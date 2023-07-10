#include "stdafx.h"
#include "util.h"

#include <omp.h>
#include <iostream>
#include <optional>
#include <ranges>
#include <regex>
#include <type_traits>

void wait_for_debugger() {
#ifdef _MSC_VER
  while (!IsDebuggerPresent()) {
    LOG(WARNING) << "waiting...";
    Sleep(1000);
  }
  ::DebugBreak();
#else
  LOG(WARNING) << "windows only";
#endif
}

std::optional<Solution> create_random_solution(const Problem& problem,
                                               Xorshift& rnd,
                                               double timelimit) {
  Timer timer;

  // constexpr double eps = 1e-8;
  constexpr double eps = 0;
  constexpr double margin = 10.0 + eps;
  double bb_left = problem.stage_x + margin;
  double bb_right = problem.stage_x + problem.stage_w - margin;
  double bb_bottom = problem.stage_y + margin;
  double bb_top = problem.stage_y + problem.stage_h - margin;

  std::vector<Placement> placements;

  while (timer.elapsed_ms() < timelimit &&
         placements.size() < problem.musicians.size()) {
    double x = bb_left + (bb_right - bb_left) * rnd.next_double();
    double y = bb_bottom + (bb_top - bb_bottom) * rnd.next_double();
    bool overlap = false;
    for (const auto& [px, py] : placements) {
      if ((x - px) * (x - px) + (y - py) * (y - py) < margin * margin) {
        overlap = true;
        break;
      }
    }
    if (overlap)
      continue;
    placements.emplace_back(x, y);
  }

  if (placements.size() < problem.musicians.size())
    return std::nullopt;

  Solution solution;
  solution.placements = placements;
  solution.set_default_volumes();
  return solution;
}

bool is_valid_solution(const Problem& problem,
                       const Solution& solution,
                       bool verbose) {
#define CHECK_VALID(expr, msg) \
  if (!(expr)) {               \
    if (verbose) {             \
      valid = false;           \
      LOG(INFO) << (msg);      \
    } else {                   \
      return false;            \
    }                          \
  }

  bool valid = true;
  CHECK_VALID(
      problem.musicians.size() == solution.placements.size(),
      format("size mismatch: problem.mucisians(%d) != solution.placements(%d)",
             problem.musicians.size(), solution.placements.size()));

  for (int k = 0; k < problem.musicians.size(); k++) {
    CHECK_VALID(is_musician_on_stage(problem, solution.placements[k]),
                format("bad musician position: %d (%f, %f)", k,
                       solution.placements[k].x, solution.placements[k].y));
    for (int kk = k + 1; kk < problem.musicians.size(); kk++) {
      CHECK_VALID(
          !are_musicians_too_close(solution.placements[k],
                                   solution.placements[kk]),
          format("musicians are too close: %d (%f, %f) .. %d (%f, %f) "
                 "(distance = %f)",
                 k, solution.placements[k].x, solution.placements[k].y, kk,
                 solution.placements[kk].x, solution.placements[kk].y,
                 std::sqrt(distance_squared(solution.placements[k],
                                            solution.placements[kk]))));
    }
  }

  if (problem.extension.consider_pillars) {
    for (int i = 0; i < problem.attendees.size(); i++) {
      for (int p = 0; p < problem.pillars.size(); p++) {
        CHECK_VALID(distance_squared(problem.attendees[i], problem.pillars[p]) >
                        problem.pillars[p].r * problem.pillars[p].r,
                    format("attendee is inside a pillar: %d (%f, %f) .. %d "
                           "(%f, %f, %f) (distance = %f)",
                           i, problem.attendees[i].x, problem.attendees[i].y, p,
                           problem.pillars[p].x, problem.pillars[p].y,
                           problem.pillars[p].r,
                           std::sqrt(distance_squared(problem.attendees[i],
                                                      problem.pillars[p]))));
      }
    }
  }

  CHECK_VALID(solution.volumes.empty() ||
                  solution.volumes.size() == problem.musicians.size(),
              format("invalid volumes of size %d where musicians are %d",
                     solution.volumes.size(), problem.musicians.size()));
  if (!solution.volumes.empty()) {
    for (int i = 0; i < solution.volumes.size(); ++i) {
      CHECK_VALID(
          0.0 <= solution.volumes[i] && solution.volumes[i] <= 10.0,
          format("volume[%d] = %f out of range", i, solution.volumes[i]));
    }
  }

#undef CHECK_VALID
  if (verbose) {
    LOG(INFO) << "solution is valid";
  }
  return valid;
}

int64_t compute_score(const Problem& problem, const Solution& solution) {
  const auto& musicians = problem.musicians;
  const auto& attendees = problem.attendees;
  const auto& pillars = problem.pillars;
  const auto& placements = solution.placements;
  const auto& extension = problem.extension;
  const bool is_volume_available = solution.volumes.size() == musicians.size();

  std::vector<double> harmony(musicians.size(), 1.0);
  if (extension.consider_harmony) {
#pragma omp parallel for
    for (int k = 0; k < musicians.size(); k++) {
      int instrument = musicians[k];
      auto&& p1 = placements[k];
      double harmonyk = 1;
      for (int kk = 0; kk < musicians.size(); kk++) {
        if (k != kk && instrument == musicians[kk]) {
          auto&& p2 = placements[kk];
          double distance = std::hypot(p1.x - p2.x, p1.y - p2.y);
          double q = 1.0 / distance;
          harmonyk += q;
        }
      }
      harmony[k] = harmonyk;
    }
  }

  int64_t score = 0;
#pragma omp parallel for
  for (int k = 0; k < musicians.size(); k++) {
    int t = problem.musicians[k];
    double mx = placements[k].x, my = placements[k].y;
    for (int i = 0; i < attendees.size(); i++) {
      double ax = attendees[i].x, ay = attendees[i].y;
      bool blocked = false;
      for (int kk = 0; kk < musicians.size(); kk++) {
        double cx = placements[kk].x, cy = placements[kk].y;
        if (k != kk && is_intersect(cx, cy, 5.0, mx, my, ax, ay)) {
          blocked = true;
          break;
        }
      }

      for (int p = 0; p < pillars.size(); p++) {
        if (is_intersect(pillars[p], pillars[p].r, placements[k],
                         attendees[i])) {
          blocked = true;
          break;
        }
      }
      if (blocked)
        continue;

      double d2 = (ax - mx) * (ax - mx) + (ay - my) * (ay - my);
      double taste = attendees[i].tastes[t];
      double qk = harmony[k];
#pragma omp critical(crit_sct)
      {
        double impact = ceil(1e6 * taste / d2);
        score += (int64_t)ceil(
            (is_volume_available ? solution.volumes[k] : 1.0) * qk * impact);
      }
    }
  }

  return score;
}

std::vector<int64_t> compute_score_each_musician(const Problem& problem,
                                                 const Solution& solution) {
  const auto& musicians = problem.musicians;
  const auto& attendees = problem.attendees;
  const auto& pillars = problem.pillars;
  const auto& placements = solution.placements;
  const auto& extension = problem.extension;
  const bool is_volume_available = solution.volumes.size() == musicians.size();

  std::vector<double> harmony(musicians.size(), 1.0);
  if (extension.consider_harmony) {
#pragma omp parallel for
    for (int k = 0; k < musicians.size(); k++) {
      int instrument = musicians[k];
      auto&& p1 = placements[k];
      double harmonyk = 1;
      for (int kk = 0; kk < musicians.size(); kk++) {
        if (k != kk && instrument == musicians[kk]) {
          auto&& p2 = placements[kk];
          double distance = std::hypot(p1.x - p2.x, p1.y - p2.y);
          double q = 1.0 / distance;
          harmonyk += q;
        }
      }
      harmony[k] = harmonyk;
    }
  }

  std::vector<int64_t> scores(musicians.size());
#pragma omp parallel for
  for (int k = 0; k < musicians.size(); k++) {
    int t = problem.musicians[k];
    double mx = placements[k].x, my = placements[k].y;
    for (int i = 0; i < attendees.size(); i++) {
      double ax = attendees[i].x, ay = attendees[i].y;
      bool blocked = false;
      for (int kk = 0; kk < musicians.size(); kk++) {
        double cx = placements[kk].x, cy = placements[kk].y;
        if (k != kk && is_intersect(cx, cy, 5.0, mx, my, ax, ay)) {
          blocked = true;
          break;
        }
      }

      for (int p = 0; p < pillars.size(); p++) {
        if (is_intersect(pillars[p], pillars[p].r, placements[k],
                         attendees[i])) {
          blocked = true;
          break;
        }
      }
      if (blocked)
        continue;

      double d2 = (ax - mx) * (ax - mx) + (ay - my) * (ay - my);
      double taste = attendees[i].tastes[t];
      double qk = harmony[k];
#pragma omp critical(crit_sct)
      {
        double impact = ceil(1e6 * taste / d2);
        scores[k] += (int64_t)ceil(
            (is_volume_available ? solution.volumes[k] : 1.0) * qk * impact);
      }
    }
  }

  return scores;
}

void set_optimal_volumes(const Problem& problem,
                         Solution& solution,
                         double amplitude) {
  auto scores = compute_score_each_musician(problem, solution);
  for (int musician_id = 0; musician_id < problem.musicians.size();
       musician_id++) {
    solution.volumes[musician_id] = scores[musician_id] <= 0 ? 0.0 : amplitude;
  }
}

void set_constant_volumes(const Problem &problem,
                          Solution &solution,
                          double amplitude) {
  for (int musician_id = 0; musician_id < problem.musicians.size(); musician_id++) {
    solution.volumes[musician_id] = amplitude;
  }
}

struct SegmentMap : public std::map<double, double> {
  bool included(double x) {
    auto it = upper_bound(x);
    if (it == begin() || (--it)->second < x)
      return false;
    return true;
  }

  void insert(double l, double r) {
    auto itl = upper_bound(l), itr = upper_bound(r);
    if (itl != begin()) {
      if ((--itl)->second < l)
        ++itl;
    }
    if (itl != itr) {
      l = std::min(l, itl->first);
      r = std::max(r, std::prev(itr)->second);
      erase(itl, itr);
    }
    (*this)[l] = r;
  }
};

int64_t compute_score_fast(const Problem& problem, const Solution& solution) {
  static bool warned = false;
  if (!warned) {
    LOG(WARNING) << __FUNCTION__
                 << " does not yet support full division extension";
    warned = true;
  }

  const auto& musicians = problem.musicians;
  const auto& attendees = problem.attendees;
  const auto& placements = solution.placements;

  const double PI = atan2(0.0, -1.0);

  std::vector<int64_t> attendee_score(attendees.size());
  std::vector<SegmentMap> smaps(attendees.size());

#pragma omp parallel for
  for (int i = 0; i < (int)attendees.size(); i++) {
    double ax = attendees[i].x, ay = attendees[i].y;

    std::vector<double> dist2;
    dist2.reserve(musicians.size());
    for (const auto& [mx, my] : placements) {
      dist2.push_back((ax - mx) * (ax - mx) + (ay - my) * (ay - my));
    }

    std::vector<int> ord(musicians.size());
    std::iota(ord.begin(), ord.end(), 0);
    std::sort(ord.begin(), ord.end(),
              [&](int i, int j) { return dist2[i] < dist2[j]; });

    auto& score = attendee_score[i];
    auto& smap = smaps[i];

    for (int mid : ord) {
      auto [mx, my] = placements[mid];
      mx -= ax;
      my -= ay;
      double d2 = dist2[mid];
      double rad = atan2(my, mx);
      double drad = atan2(5.0, sqrt(d2 - 25.0));
      // DUMP(rad, drad);
      //  if rad is included in smap, it is not blocked.
      if (!smap.included(rad)) {
        double taste = attendees[i].tastes[musicians[mid]];
        score += (int64_t)ceil(1e6 * taste / d2);
      }
      // merge [rad-drad, rad+drad] into smap
      double left = rad - drad, right = rad + drad;
      if (left < -PI) {
        smap.insert(left + PI * 2, PI);
        left = -PI;
      }
      if (right > PI) {
        smap.insert(-PI, right - 2 * PI);
        right = PI;
      }
      smap.insert(left, right);
    }
  }

  return std::accumulate(attendee_score.begin(), attendee_score.end(), 0LL);
}

nlohmann::json create_problem_stats(Problem& problem) {
  auto [room_width, room_height, stage_w, stage_h, stage_x, stage_y, musicians,
        attendees, pillars, extension, problem_id] = problem;
  nlohmann::json stats;
  stats["topology"] = {
      {"room_width", room_width},
      {"room_height", room_height},
      {"stage_w", stage_w},
      {"stage_h", stage_h},
      {"stage_x", stage_x},
      {"stage_y", stage_y},
      {"stage_ratio", stage_w * stage_h / room_width / room_height}};
  stats["complexity"] = {
      {"num_musicians", musicians.size()},
      {"num_attendees", attendees.size()},
      {"num_pillars", pillars.size()},
      {"complexity_index",
       musicians.size() * musicians.size() * attendees.size()}};
  {
    std::map<int, int> type_hist;
    for (int type : musicians) {
      type_hist[type]++;
    }
    int num_musician_types = type_hist.size();
    std::pair<int, int> mode;
    double sum = 0.0, sqsum = 0.0;
    for (const auto& [type, ctr] : type_hist) {
      if (mode.second < ctr) {
        mode = {type, ctr};
      }
      sum += ctr;
      sqsum += ctr * ctr;
    }
    double mean = sum / num_musician_types;
    double var = sqsum / num_musician_types - mean * mean;
    stats["musicians"] = {{"num_types", num_musician_types},
                          {"mode_type", mode.first},
                          {"mode_amount", mode.second},
                          {"mean_amount", mean},
                          {"stdev_amount", sqrt(var)}};
  }
  {
    std::vector<double> tastes;
    for (const auto& attendee : attendees) {
      for (auto taste : attendee.tastes) {
        tastes.push_back(taste);
      }
    }
    std::sort(tastes.begin(), tastes.end());
    double min = 1e20, max = -1e20, sum = 0.0, sqsum = 0.0;
    for (auto taste : tastes) {
      sum += taste;
      sqsum += taste * taste;
      chmin(min, taste);
      chmax(max, taste);
    }
    double mean = sum / tastes.size();
    double var = sqsum / tastes.size() - mean * mean;
    double median = tastes[tastes.size() / 2];
    stats["attendees"] = {{"mean_taste", mean},
                          {"stdev_taste", sqrt(var)},
                          {"median_taste", median},
                          {"min_taste", min},
                          {"max_taste", max}};
  }
  if (!pillars.empty()) {
    const double PI = abs(atan2(0, -1));
    double min = 1e20, max = -1e20, sum = 0.0, sqsum = 0.0;
    for (const auto& [x, y, r] : pillars) {
      chmin(min, r);
      chmax(max, r);
      sum += r;
      sqsum += r * r;
    }
    double mean = sum / pillars.size();
    double var = sqsum / pillars.size() - mean * mean;
    stats["pillars"] = {{"mean_radius", mean},
                        {"stdev_radius", sqrt(var)},
                        {"min_radius", min},
                        {"max_radius", max},
                        {"area_ratio", sqsum * PI / room_width / room_height}};
  }
  return stats;
}
