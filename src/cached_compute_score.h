#pragma once

#include <stdafx.h>
#ifdef _MSC_VER
#include <ppl.h>
#endif
#include <glog/logging.h>

#include "spec.h"
#include "util.h"

struct CachedComputeScore {
 public:
  const Problem& m_problem;
  const int m_num_attendees;
  const int m_num_musicians;
  const int m_num_pillars;
  Solution m_solution;

  // private:
  std::vector<double>
      m_harmony_cache;  // 1 + m_harmony_cache[k] is the harmony.
  std::vector<int64_t> m_score_cache;
  std::vector<uint8_t> m_audible_cache;
  std::vector<int16_t> m_blocker_count_cache;
  int64_t m_score;

  double m_accum_elapsed_ms_partial = 0.0;
  int m_call_count_partial = 0;
  double m_accum_elapsed_ms_full = 0.0;
  int m_call_count_full = 0;

  bool m_compute_affected = true;
  std::vector<int> m_affected_musicians_indices;

  int sample_random_affected_musician(Xorshift& rnd) const {
    if (m_affected_musicians_indices.empty()) {
      return rnd.next_int() % m_num_musicians;
    } else {
      return m_affected_musicians_indices[rnd.next_int() %
                                          m_affected_musicians_indices.size()];
    }
  }

  double get_mean_elapsed_ms_partial() const {
    return m_accum_elapsed_ms_partial / (m_call_count_partial + 1e-9);
  }
  double get_mean_elapsed_ms_full() const {
    return m_accum_elapsed_ms_full / (m_call_count_full + 1e-9);
  }
  void report() const {
    LOG(INFO) << format(
        "PID=%d: full calculation %.4f ms/eval (%d evals), partial update %.4f "
        "ms/eval (%d evals)",
        m_problem.problem_id, get_mean_elapsed_ms_full(), m_call_count_full,
        get_mean_elapsed_ms_partial(), m_call_count_partial);
  }

 public:
  CachedComputeScore(const Problem& problem);

  int64_t& partial_score(int k, int i) {
#if 0
        LOG_ASSERT(0 <= k && k < m_num_musicians);
        LOG_ASSERT(0 <= i && i < m_num_attendees);
#endif
    return m_score_cache[k * m_num_attendees + i];
  }

  int16_t& blocker_count(int k, int i) {
#if 0
        LOG_ASSERT(0 <= k && k < m_num_musicians);
        LOG_ASSERT(0 <= i && i < m_num_attendees);
#endif
    return m_blocker_count_cache[k * m_num_attendees + i];
  }

  // k_src's ray to i is audible by k_other
  uint8_t& partial_audible(int k_src, int k_other, int i) {
#if 0
        LOG_ASSERT(0 <= k_src && k_src < m_num_musicians);
        LOG_ASSERT(0 <= k_other && k_other < m_num_musicians);
        LOG_ASSERT(k_src != k_other); // not illegal but strange.
        LOG_ASSERT(0 <= i && i < m_num_attendees);
#endif
    return m_audible_cache[(k_src * m_num_musicians + k_other) *
                               m_num_attendees +
                           i];
  }
  int64_t score() const { return m_score; }

  int64_t change_musician_volume(int k_changed, double curr_volume);

  int64_t change_musician(int k_changed, const Placement& curr_placement, bool dry_run = false);

  int64_t full_compute(const Solution& solution);
};
