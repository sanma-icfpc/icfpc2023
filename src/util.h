#pragma once

#include <stdafx.h>
#ifdef _MSC_VER
#include <ppl.h>
#endif
#include <glog/logging.h>

#include <climits>

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
