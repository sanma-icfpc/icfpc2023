#include "stdafx.h"
#include "cached_compute_score.h"

#include "spec.h"

CachedComputeScore::CachedComputeScore(const Problem& problem)
    : m_problem(problem),
      m_num_attendees(problem.attendees.size()),
      m_num_musicians(problem.musicians.size()),
      m_num_pillars(problem.extension.consider_pillars ? problem.pillars.size()
                                                       : 0),
      m_harmony_cache(m_num_musicians, 0.0),
      m_score_cache(problem.attendees.size() * problem.musicians.size(), 0),
      m_audible_cache(problem.attendees.size() * problem.musicians.size() *
                          problem.musicians.size(),
                      1),
      m_blocker_count_cache(problem.attendees.size() * problem.musicians.size(),
                            0) {}

int64_t CachedComputeScore::change_musician(int k_changed,
                                            const Placement& curr_placement, bool dry_run) {
  Timer timer;

  const Placement prev_placement = m_solution.placements[k_changed];
  m_solution.placements[k_changed] = curr_placement;
  const auto& musicians = m_problem.musicians;
  const auto& attendees = m_problem.attendees;
  const auto& placements = m_solution.placements;
  const auto& pillars = m_problem.pillars;
  const int M = musicians.size();
  const int A = attendees.size();
  const int P = pillars.size();
  LOG_ASSERT(M == placements.size());

  std::vector<uint8_t> musicians_affected;
  if (m_compute_affected) {
    musicians_affected.assign(M, 0);
    m_affected_musicians_indices.clear();
  }

  auto harmony_cache_backup = m_harmony_cache;

  // スコアの更新前に、ブロック状況の更新が必要(ブロックは新旧両方を同時に利用するため)
  // pillarはblocker_countに加味されているので特別扱いする必要は無い
  int64_t influence_diff = 0;
#pragma omp parallel for reduction(+ : influence_diff)
  for (auto k_src = 0; k_src < M; ++k_src) {
    if (k_changed == k_src)
      continue;

    double old_harmony = 0.0, new_harmony = 0.0;
    if (m_problem.extension.consider_harmony) {
      if (musicians[k_changed] == musicians[k_src]) {
        old_harmony = m_harmony_cache[k_src];
        m_harmony_cache[k_src] +=
            inverse_distance(placements[k_src], curr_placement) -
            inverse_distance(placements[k_src], prev_placement);
        new_harmony = m_harmony_cache[k_src];
      }
    }

    for (int i = 0; i < A; ++i) {
      const bool old_audible = partial_audible(k_src, k_changed, i);
      const bool new_audible =
          !is_intersect(placements[k_changed], k_musician_radius,
                        placements[k_src], attendees[i]);
      if (!dry_run)
        partial_audible(k_src, k_changed, i) = new_audible;

      const auto old_blocker_count = blocker_count(k_src, i);
      auto new_blocker_count = old_blocker_count;
      if (old_audible && !new_audible)
        new_blocker_count += 1;
      if (!old_audible && new_audible)
        new_blocker_count -= 1;
      if (!dry_run)
        blocker_count(k_src, i) = new_blocker_count;

      if (m_compute_affected &&
          (old_blocker_count == 0) != (new_blocker_count == 0)) {
        musicians_affected[k_src] = 1;
      }

      int64_t partial_influence_diff = 0;
      partial_influence_diff -= old_blocker_count == 0
                           ? (int64_t)std::ceil(m_solution.volumes[k_src] *
                                                (1.0 + old_harmony) *
                                                partial_score(k_src, i))
                           : 0;
      partial_influence_diff += new_blocker_count == 0
                           ? (int64_t)std::ceil(m_solution.volumes[k_src] *
                                                (1.0 + new_harmony) *
                                                partial_score(k_src, i))
                           : 0;
      influence_diff += partial_influence_diff;
    }
  }

  if (m_compute_affected) {
    for (int k = 0; k < M; ++k) {
      if (musicians_affected[k]) {
        m_affected_musicians_indices.push_back(k);
      }
    }
  }

  const double old_harmony = m_harmony_cache[k_changed];
  if (m_problem.extension.consider_harmony) {
    for (int k_other = 0; k_other < m_num_musicians; ++k_other) {
      if (k_other != k_changed && musicians[k_changed] == musicians[k_other]) {
        m_harmony_cache[k_changed] +=
            inverse_distance(placements[k_other], curr_placement) -
            inverse_distance(placements[k_other], prev_placement);
      }
    }
  }
  const double new_harmony = m_harmony_cache[k_changed];

  // 移動したmusicianが得る効果の増減
#pragma omp parallel for reduction(+ : influence_diff)
  for (auto i = 0; i < A; ++i) {
    const int64_t old_blocker_count = blocker_count(k_changed, i);
    int64_t new_blocker_count = 0;
    for (int k_other = 0; k_other < M; ++k_other) {
      if (k_other == k_changed)
        continue;
      auto new_partial_audible = 
          !is_intersect(placements[k_other], k_musician_radius,
                        placements[k_changed], attendees[i]);
      if (!new_partial_audible)
        new_blocker_count++;
      if (!dry_run)
        partial_audible(k_changed, k_other, i) = new_partial_audible;
    }
    for (int p = 0; p < P; ++p) {
      if (is_intersect(pillars[p], pillars[p].r, placements[k_changed],
                       attendees[i])) {
        new_blocker_count++;
      }
    }

    int64_t partial_influence_diff = 0;
    partial_influence_diff -= old_blocker_count == 0
                         ? (int64_t)std::ceil(m_solution.volumes[k_changed] *
                                              (1.0 + old_harmony) *
                                              partial_score(k_changed, i))
                         : 0;
    
    const int64_t new_partial_score = (int64_t)std::ceil(
          1e6 * attendees[i].tastes[musicians[k_changed]] /
          distance_squared(placements[k_changed], attendees[i]));

    partial_influence_diff += new_blocker_count == 0
                         ? (int64_t)std::ceil(m_solution.volumes[k_changed] *
                                              (1.0 + new_harmony) *
                                              new_partial_score)
                         : 0;
    influence_diff += partial_influence_diff;

    if (!dry_run) {
      blocker_count(k_changed, i) = new_blocker_count;
      partial_score(k_changed, i) = new_partial_score;
    }
  }

  m_accum_elapsed_ms_partial += timer.elapsed_ms();
  m_call_count_partial += 1;

  if (dry_run) {
    m_solution.placements[k_changed] = prev_placement;
    m_harmony_cache = harmony_cache_backup;
  } else {
    m_score += influence_diff;
  }

  return influence_diff;
}

int64_t CachedComputeScore::change_musician_volume(int k_changed,
                                                   double curr_volume) {
  const auto& attendees = m_problem.attendees;
  const int A = attendees.size();

  const double prev_volume = m_solution.volumes[k_changed];
  m_solution.volumes[k_changed] = curr_volume;

  int64_t old_influence = 0;
  int64_t new_influence = 0;
  for (int i = 0; i < A; ++i) {
    if (blocker_count(k_changed, i) == 0) {
      const double qI =
          (1.0 + m_harmony_cache[k_changed]) * partial_score(k_changed, i);
      old_influence += (int64_t)std::ceil(prev_volume * qI);
      new_influence += (int64_t)std::ceil(curr_volume * qI);
    }
  }
  m_score += new_influence - old_influence;

  return new_influence - old_influence;
}

int64_t CachedComputeScore::full_compute(const Solution& solution) {
  Timer timer;

  m_solution = solution;
  const auto& musicians = m_problem.musicians;
  const auto& attendees = m_problem.attendees;
  const auto& placements = solution.placements;
  const auto& pillars = m_problem.pillars;
  const int M = musicians.size();
  const int A = attendees.size();
  const int P = pillars.size();

  m_score = 0;
  m_score_cache.assign(m_score_cache.size(), 0);
  m_audible_cache.assign(m_audible_cache.size(), 1);
  m_blocker_count_cache.assign(m_blocker_count_cache.size(), 0);
  m_harmony_cache.assign(m_harmony_cache.size(), 0.0);

  if (m_problem.extension.consider_harmony) {
#pragma omp parallel for
    for (int k = 0; k < M; k++) {
      double harmony = 0.0;
      for (int k_other = 0; k_other < M; ++k_other) {
        if (k != k_other && musicians[k] == musicians[k_other]) {
          harmony += inverse_distance(placements[k], placements[k_other]);
        }
      }
      m_harmony_cache[k] = harmony;
    }
  }

  int64_t score_diff = 0;
#pragma omp parallel for reduction(+ : score_diff)
  for (auto k_src = 0; k_src < M; ++k_src) {
    for (int i = 0; i < A; ++i) {
      for (int k_other = 0; k_other < M; ++k_other) {
        if (k_src != k_other) {
          partial_audible(k_src, k_other, i) =
              !is_intersect(placements[k_other], k_musician_radius,
                            placements[k_src], attendees[i]);
          if (!partial_audible(k_src, k_other, i)) {
            blocker_count(k_src, i) += 1;
          }
        }
      }
      for (int p = 0; p < P; ++p) {
        if (is_intersect(pillars[p], pillars[p].r, placements[k_src],
                         attendees[i])) {
          blocker_count(k_src, i) += 1;
        }
      }
      const bool audible = blocker_count(k_src, i) == 0;
      const int64_t influence =
          (int64_t)ceil(1e6 * attendees[i].tastes[musicians[k_src]] /
                        distance_squared(placements[k_src], attendees[i]));
      partial_score(k_src, i) = influence;
      score_diff +=
          audible ? (int64_t)ceil(m_solution.volumes[k_src] *
                                  (1.0 + m_harmony_cache[k_src]) * influence)
                  : 0;
    }
  }
  m_score += score_diff;

  m_accum_elapsed_ms_full += timer.elapsed_ms();
  m_call_count_full += 1;
  return m_score;
}
