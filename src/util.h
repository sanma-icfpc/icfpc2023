#pragma once
#include <stdafx.h>
#ifdef _MSC_VER
#include <ppl.h>
#endif
#include <glog/logging.h>
#include <omp.h>
#include <climits>
#include <iostream>
#include <optional>
#include <ranges>
#include <regex>
#include <type_traits>
#include "spec.h"

void wait_for_debugger();

template <typename T>
bool chmax(T& a, const T& b) {
  if (a < b) {
    a = b;
    return true;
  }
  return false;
}

template <typename T>
bool chmin(T& a, const T& b) {
  if (a > b) {
    a = b;
    return true;
  }
  return false;
}

template <typename T>
T SQ(const T& x) {
  return x * x;
}

/** string formatter **/
template <typename... Ts>
std::string format(const std::string& f, Ts... t) {
  size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
  std::vector<char> b(l + 1);
  std::snprintf(&b[0], l + 1, f.c_str(), t...);
  return std::string(&b[0], &b[0] + l);
}

/** rand **/
struct Xorshift {
  static constexpr uint64_t M = INT_MAX;
  static constexpr double e = 1.0 / M;
  uint64_t x = 88172645463325252LL;
  Xorshift() {}
  Xorshift(uint64_t seed) { reseed(seed); }
  inline void reseed(uint64_t seed) {
    x = 0x498b3bc5 ^ seed;
    for (int i = 0; i < 20; i++)
      next();
  }
  inline uint64_t next() {
    x = x ^ (x << 7);
    return x = x ^ (x >> 9);
  }
  inline int next_int() { return next() & M; }
  inline int next_int(int mod) { return next() % mod; }
  inline int next_int(int l, int r) { return l + next_int(r - l + 1); }
  inline double next_double() { return next_int() * e; }
};

/** timer **/
class Timer {
  double t = 0, paused = 0, tmp;

 public:
  Timer() { reset(); }
  static double time() {
#ifdef _MSC_VER
    return __rdtsc() / 3.0e9;
#else
    unsigned long long a, d;
    __asm__ volatile("rdtsc" : "=a"(a), "=d"(d));
    return (d << 32 | a) / 3.0e9;
#endif
  }
  void reset() { t = time(); }
  void pause() { tmp = time(); }
  void restart() { paused += time() - tmp; }
  double elapsed_ms() const { return (time() - t - paused) * 1000.0; }
};

/** dump args **/
#define DUMPOUT std::cerr
static std::ostringstream DUMPBUF;
#define DUMP(...)                                                              \
  do {                                                                         \
    DUMPBUF << "  ";                                                           \
    DUMPBUF << #__VA_ARGS__ << " :[DUMP - " << __LINE__ << ":" << __FUNCTION__ \
            << ']' << std::endl;                                               \
    DUMPBUF << "    ";                                                         \
    dump_func(__VA_ARGS__);                                                    \
    DUMPOUT << DUMPBUF.str();                                                  \
    DUMPBUF.str("");                                                           \
    DUMPBUF.clear();                                                           \
  } while (0);
inline void dump_func() {
  DUMPBUF << std::endl;
}
template <class Head, class... Tail>
void dump_func(Head&& head, Tail&&... tail) {
  DUMPBUF << head;
  if (sizeof...(Tail) == 0) {
    DUMPBUF << " ";
  } else {
    DUMPBUF << ", ";
  }
  dump_func(std::move(tail)...);
}

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, std::string>
int_to_delimited_string(T number) {
  std::string number_str = std::to_string(number);
  std::string result;

  int count = 0;
  for (auto it = number_str.rbegin(); it != number_str.rend(); ++it) {
    if (count != 0 && count % 3 == 0) {
      result.push_back('_');
    }
    result.push_back(*it);
    count++;
  }

  std::ranges::reverse(result);
  return result;
}

inline double is_intersect(double cx,
                           double cy,
                           double r,
                           double x1,
                           double y1,
                           double x2,
                           double y2) {
  if (false) {
    // 1.19->1.23ms/partialeval
    // full BBOX test
    if (std::min(cx + r, std::max(x1, x2)) -
                std::max(cx - r, std::min(x1, x2)) >=
            0 &&
        std::min(cy + r, std::max(y1, y2)) -
                std::max(cy - r, std::min(y1, y2)) >=
            0) {
      return false;
    }
    // 1.19->1.08ms/partialeval
    // separable test
    if (std::max(x1, x2) > cx - r || std::max(y1, y2) > cy - r ||
        cx + r < std::min(x1, x2) || cy + r < std::min(y1, y2)) {
      return false;
    }
  }
  x1 -= cx;
  x2 -= cx;
  y1 -= cy;
  y2 -= cy;
  double dx = x2 - x1, dy = y2 - y1;
  double a = dx * dx + dy * dy, b = x1 * dx + y1 * dy,
         c = x1 * x1 + y1 * y1 - r * r;
  if (b * b - a * c < 0)
    return false;
  double f0 = c, f1 = a + 2 * b + c;
  if (f0 * f1 <= 0)
    return true;
  return 0 <= f0 && 0 <= f1 && -a <= b && b <= 0;
}

template <typename T1, typename T2, typename T3>
inline double is_intersect(const T1& c, double r, const T2& p1, const T3& p2) {
  return is_intersect(c.x, c.y, r, p1.x, p1.y, p2.x, p2.y);
}

template <typename T1, typename T2>
inline double distance_squared(const T1& p1, const T2& p2) {
  return SQ(p1.x - p2.x) + SQ(p1.y - p2.y);
}

template <typename T1, typename T2>
inline double inverse_distance(const T1& p1, const T2& p2) {
  return 1.0 / std::hypot(p1.x - p2.x, p1.y - p2.y);
}

inline bool is_musician_on_stage(const Problem& problem,
                                 const Placement& p,
                                 double eps_margin_for_safty = 0.0) {
  return problem.stage_x + k_musician_spacing_radius + eps_margin_for_safty <=
             p.x &&
         p.x + k_musician_spacing_radius + eps_margin_for_safty <=
             problem.stage_x + problem.stage_w &&
         problem.stage_y + k_musician_spacing_radius + eps_margin_for_safty <=
             p.y &&
         p.y + k_musician_spacing_radius + eps_margin_for_safty <=
             problem.stage_y + problem.stage_h;
}

inline bool are_musicians_too_close(const Placement& p1,
                                    const Placement& p2,
                                    double eps_margin_for_safty = 1e-3) {
  // To ensure they have enough room for playing, they must not
  // have any other musician or an edge of the stage in a circle of radius 10
  // centered on them. A placement where a musician is at distance exactly 10
  // from an edge or another musician will be considered as valid.
  return distance_squared(p1, p2) + eps_margin_for_safty <=
         k_musician_spacing_radius * k_musician_spacing_radius;
}

inline std::optional<Placement> suggest_random_position(
    const Problem& problem,
    const Solution& solution,
    Xorshift& rnd,
    int i,
    int max_retry = 100) {
  for (int retry = 0; retry < max_retry; ++retry) {
    Placement placement{
        problem.stage_x + k_musician_spacing_radius +
            rnd.next_double() *
                (problem.stage_w - k_musician_spacing_radius * 2),
        problem.stage_y + k_musician_spacing_radius +
            rnd.next_double() *
                (problem.stage_h - k_musician_spacing_radius * 2)};
    if (!is_musician_on_stage(problem, placement))
      continue;
    bool conflict = false;
    for (int kk = 0; kk < solution.placements.size(); ++kk) {
      if (i != kk) {
        if (are_musicians_too_close(solution.placements[kk], placement)) {
          conflict = true;
          break;
        }
      }
    }
    if (!conflict)
      return placement;
  }

  return std::nullopt;
}

std::optional<Solution> create_random_solution(const Problem& problem,
                                               Xorshift& rnd,
                                               double timelimit = 1000);

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

  int64_t change_musician_volume(int k_changed, double curr_volume) {
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

  int64_t change_musician(int k_changed, const Placement& curr_placement) {
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

    // スコアの更新前に、ブロック状況の更新が必要(ブロックは新旧両方を同時に利用するため)
    // pillarはblocker_countに加味されているので特別扱いする必要は無い
    int64_t old_influence = 0;
    int64_t new_influence = 0;
#pragma omp parallel for reduction(+ : old_influence) \
    reduction(+ : new_influence)
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
        partial_audible(k_src, k_changed, i) =
            new_audible;  // この二重ループでは全て独立

        const auto old_blocker_count = blocker_count(k_src, i);
        if (old_audible && !new_audible)
          blocker_count(k_src, i) += 1;
        if (!old_audible && new_audible)
          blocker_count(k_src, i) -= 1;
        const auto new_blocker_count = blocker_count(k_src, i);

        if (m_compute_affected &&
            (old_blocker_count == 0) != (new_blocker_count == 0)) {
          musicians_affected[k_src] = 1;
        }
        old_influence += old_blocker_count == 0
                             ? (int64_t)std::ceil(m_solution.volumes[k_src] *
                                                  (1.0 + old_harmony) *
                                                  partial_score(k_src, i))
                             : 0;
        new_influence += new_blocker_count == 0
                             ? (int64_t)std::ceil(m_solution.volumes[k_src] *
                                                  (1.0 + new_harmony) *
                                                  partial_score(k_src, i))
                             : 0;
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
        if (k_other != k_changed &&
            musicians[k_changed] == musicians[k_other]) {
          m_harmony_cache[k_changed] +=
              inverse_distance(placements[k_other], curr_placement) -
              inverse_distance(placements[k_other], prev_placement);
        }
      }
    }
    const double new_harmony = m_harmony_cache[k_changed];

    // 移動したmusicianが得る効果の増減
#pragma omp parallel for reduction(+ : old_influence) \
    reduction(+ : new_influence)
    for (auto i = 0; i < A; ++i) {
      const int64_t old_blocker_count = blocker_count(k_changed, i);
      int64_t new_blocker_count = 0;
      for (int k_other = 0; k_other < M; ++k_other) {
        if (k_other == k_changed)
          continue;
        partial_audible(k_changed, k_other, i) =
            !is_intersect(placements[k_other], k_musician_radius,
                          placements[k_changed], attendees[i]);
        if (!partial_audible(k_changed, k_other, i))
          new_blocker_count++;
      }
      for (int p = 0; p < P; ++p) {
        if (is_intersect(pillars[p], pillars[p].r, placements[k_changed],
                         attendees[i])) {
          new_blocker_count++;
        }
      }
      blocker_count(k_changed, i) = new_blocker_count;
      old_influence += old_blocker_count == 0
                           ? (int64_t)std::ceil(m_solution.volumes[k_changed] *
                                                (1.0 + old_harmony) *
                                                partial_score(k_changed, i))
                           : 0;
      partial_score(k_changed, i) = (int64_t)std::ceil(
          1e6 * attendees[i].tastes[musicians[k_changed]] /
          distance_squared(placements[k_changed], attendees[i]));
      new_influence += new_blocker_count == 0
                           ? (int64_t)std::ceil(m_solution.volumes[k_changed] *
                                                (1.0 + new_harmony) *
                                                partial_score(k_changed, i))
                           : 0;
    }

    m_score += new_influence - old_influence;

    m_accum_elapsed_ms_partial += timer.elapsed_ms();
    m_call_count_partial += 1;

    return new_influence - old_influence;
  }

  int64_t full_compute(const Solution& solution) {
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
};

bool is_valid_solution(const Problem& problem,
                       const Solution& solution,
                       bool verbose = false);

int64_t compute_score(const Problem& problem, const Solution& solution);

std::vector<int64_t> compute_score_each_musician(const Problem& problem,
                                                 const Solution& solution);

void set_optimal_volumes(const Problem& problem,
                         Solution& solution,
                         double amplitude = 10.0);

int64_t compute_score_fast(const Problem& problem, const Solution& solution);

nlohmann::json create_problem_stats(Problem& problem);