#pragma once
#ifdef _MSC_VER
#include <ppl.h>
#elif _OPENMP
#include <omp.h>
#endif
#include <type_traits>
#include <iostream>
#include <climits>
#include <regex>
#include <optional>
#include <glog/logging.h>
#include "spec.h"

template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }
template<typename T> T SQ(const T& x) { return x * x; }

/** string formatter **/
template<typename... Ts>
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
    inline void reseed(uint64_t seed) { x = 0x498b3bc5 ^ seed; for (int i = 0; i < 20; i++) next(); }
    inline uint64_t next() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
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
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
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
#define DUMP(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<']'<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
inline void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }

template <typename T>
inline std::enable_if_t<std::is_integral_v<T>, std::string> int_to_delimited_string(T number) {
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

inline double is_intersect(double cx, double cy, double r, double x1, double y1, double x2, double y2) {
    x1 -= cx; x2 -= cx; y1 -= cy; y2 -= cy;
    double dx = x2 - x1, dy = y2 - y1;
    double a = dx * dx + dy * dy, b = x1 * dx + y1 * dy, c = x1 * x1 + y1 * y1 - r * r;
    if (b * b - a * c < 0) return false;
    double f0 = c, f1 = a * 2 * b + c;
    if (f0 * f1 <= 0) return true;
    return 0 <= f0 && 0 <= f1 && -a <= b && b <= 0 && a * c <= b * b;
}

template <typename T1, typename T2, typename T3>
inline double is_intersect(const T1& c, double r, const T2& p1, const T3& p2) {
    return is_intersect(c.x, c.y, r, p1.x, p1.y, p2.x, p2.y);
}

template <typename T1, typename T2>
inline double distance_squared(const T1& p1, const T2& p2) {
    return SQ(p1.x - p2.x) + SQ(p1.y - p2.y);
}

inline bool are_musicians_too_close(const Placement& p1, const Placement& p2, double eps_margin_for_safty = 1e-3) {
    // To ensure they have enough room for playing, they must not
    // have any other musician or an edge of the stage in a circle of radius 10 centered
    // on them. A placement where a musician is at distance exactly 10 from an edge or
    // another musician will be considered as valid. 
    return distance_squared(p1, p2) + eps_margin_for_safty <= k_musician_spacing_radius * k_musician_spacing_radius;
}

inline std::optional<Solution> create_random_solution(const Problem& problem, Xorshift& rnd, double timelimit = 1000) {

    Timer timer;

    //constexpr double eps = 1e-8;
    constexpr double eps = 0;
    constexpr double margin = 10.0 + eps;
    double bb_left = problem.stage_x + margin;
    double bb_right = problem.stage_x + problem.stage_w - margin;
    double bb_bottom = problem.stage_y + margin;
    double bb_top = problem.stage_y + problem.stage_h - margin;

    std::vector<Placement> placements;

    while (timer.elapsed_ms() < timelimit && placements.size() < problem.musicians.size()) {
        double x = bb_left + (bb_right - bb_left) * rnd.next_double();
        double y = bb_bottom + (bb_top - bb_bottom) * rnd.next_double();
        bool overlap = false;
        for (const auto& [px, py] : placements) {
            if ((x - px) * (x - px) + (y - py) * (y - py) < margin * margin) {
                overlap = true;
                break;
            }
        }
        if (overlap) continue;
        placements.emplace_back(x, y);
    }

    if (placements.size() < problem.musicians.size()) return std::nullopt;

    Solution solution;
    solution.placements = placements;
    return solution;

}

inline std::optional<int> guess_problem_id(std::string some_file_path) {
    std::regex pattern(".*-(\\d+)");
    std::smatch matches;

    if (std::regex_search(some_file_path, matches, pattern)) {
        std::string numberString = matches[1].str();
        try {
            return stoi(numberString);
        } catch (const std::invalid_argument& e) {
            LOG(ERROR) << "Invalid argument: " << e.what();
        } catch (const std::out_of_range& e) {
            LOG(ERROR) << "Out of range: " << e.what();
        }
    }
    return std::nullopt;
}

// still has some bugs..
struct CachedComputeScore {
public:
    const Problem& m_problem;
    const int m_num_attendees;
    const int m_num_musicians;
    Solution m_solution;
private:
    std::vector<int64_t> m_score_k_i;
    std::vector<uint8_t> m_block;
    int64_t m_score;
public:
    CachedComputeScore(const Problem& problem)
     : m_problem(problem)
     , m_num_attendees(problem.attendees.size())
     , m_num_musicians(problem.musicians.size())
     , m_score_k_i(problem.attendees.size() * problem.musicians.size(), 0)
     , m_block(problem.attendees.size() * problem.musicians.size() * problem.musicians.size(), 0) {
    }

    int64_t& partial_score(int k, int i) { return m_score_k_i[k * m_num_attendees + i]; }
    int64_t partial_score(int k, int i) const { return m_score_k_i[k * m_num_attendees + i]; }
    // k's ray to i is blocked by kk
    uint8_t& partial_block(int k, int kk, int i) { return m_block[(k * m_num_musicians + kk) * m_num_attendees + i]; }
    uint8_t partial_block(int k, int kk, int i) const { return m_block[(k * m_num_musicians + kk) * m_num_attendees + i]; } 
    int64_t score() const { return m_score; }

    int64_t change_musician(int k_changed, const Placement& curr_placement) {
        const auto& musicians = m_problem.musicians;
        const auto& attendees = m_problem.attendees;
        const auto& placements = m_solution.placements;
        LOG_ASSERT(musicians.size() == placements.size());

        int64_t score_gain = 0;
        { // score that musician k_changed offers.
            int t = m_problem.musicians[k_changed];
#pragma omp parallel for reduction(+:score_gain)
            for (int i = 0; i < attendees.size(); i++) {
                double ax = attendees[i].x, ay = attendees[i].y;
                bool blocked = false;
                bool prev_blocked = false;
                for (int kk = 0; kk < musicians.size(); kk++) {
                    if (k_changed == kk) continue;
                    if (partial_block(k_changed, kk, i)) {
                        prev_blocked = true;
                    }
                    if (is_intersect(placements[kk], 5.0, curr_placement, attendees[i])) {
                        blocked = true;
                        partial_block(k_changed, kk, i) = true;
                    } else {
                        partial_block(k_changed, kk, i) = false;
                    }
                }
                double d2 = distance_squared(attendees[i], curr_placement);
                double taste = attendees[i].tastes[t];
                int64_t prev = partial_score(k_changed, i);
                int64_t curr = (int64_t)ceil(1e6 * taste / d2);
                score_gain += (blocked ? 0 : curr) - (prev_blocked ? 0 : prev);
                partial_score(k_changed, i) = curr;
            }
        }

        // score that musician k_changed blocks/unblocks k -> i.
#pragma omp parallel for reduction(+:score_gain)
        for (int k = 0; k < musicians.size(); k++) {
            if (k == k_changed) continue;
            int t = m_problem.musicians[k];
            for (int i = 0; i < attendees.size(); i++) {
                bool blocked = is_intersect(curr_placement, 5.0, placements[k], attendees[i]);
                double taste = attendees[i].tastes[t];
                int64_t prev = partial_score(k, i);
                score_gain += (blocked ? 0 : prev) - (partial_block(k, k_changed, i) ? 0 : prev);
                partial_block(k, k_changed, i) = blocked;
            }
        }

        // update musician.
        m_solution.placements[k_changed] = curr_placement;

        m_score += score_gain;
        return score_gain;
    }

    int64_t full_compute(const Solution& solution) {
        m_solution = solution;
        const auto& musicians = m_problem.musicians;
        const auto& attendees = m_problem.attendees;
        const auto& placements = solution.placements;

        m_score = 0;
        m_score_k_i.assign(m_score_k_i.size(), 0);
        m_block.assign(m_block.size(), 0);

        for (int k = 0; k < musicians.size(); k++) {
            int t = m_problem.musicians[k];
            for (int i = 0; i < attendees.size(); i++) {
                double ax = attendees[i].x, ay = attendees[i].y;
                bool blocked = false;
                for (int kk = 0; kk < musicians.size(); kk++) {
                    if (k != kk && is_intersect(placements[kk], 5.0, placements[k], attendees[i])) {
                        blocked = true;
                        partial_block(k, kk, i) = true;
                    }
                }
                double d2 = distance_squared(placements[k], attendees[i]);
                double taste = attendees[i].tastes[t];
                int64_t curr = (int64_t)ceil(1e6 * taste / d2);
                m_score += blocked ? 0 : curr;
                partial_score(k, i) = curr;
            }
        }

        return m_score;
    }
};

#ifdef _PPL_H
inline int64_t compute_score(const Problem& problem, const Solution& solution) {

    const auto& musicians = problem.musicians;
    const auto& attendees = problem.attendees;
    const auto& placements = solution.placements;

    concurrency::combinable<int64_t> score;
    concurrency::parallel_for(0, (int)musicians.size(), [&](int k) {
        int t = problem.musicians[k];
        double mx = placements[k].x, my = placements[k].y;
        for (int i = 0; i < attendees.size(); i++) {
            double ax = attendees[i].x, ay = attendees[i].y;
            bool blocked = false;
            for (int kk = 0; kk < musicians.size(); kk++) {
                double cx = placements[kk].x, cy = placements[kk].y;
                if (is_intersect(cx, cy, 5.0, mx, my, ax, ay)) {
                    blocked = true;
                    break;
                }
            }
            if (blocked) continue;
            double d2 = (ax - mx) * (ax - mx) + (ay - my) * (ay - my);
            double taste = attendees[i].tastes[t];
            score.local() += (int64_t)ceil(1e6 * taste / d2);
        }
        });

    return score.combine(std::plus<int64_t>());

}
#else
inline int64_t compute_score(const Problem& problem, const Solution& solution) {

    const auto& musicians = problem.musicians;
    const auto& attendees = problem.attendees;
    const auto& placements = solution.placements;

    int64_t score = 0;
#ifdef _OPENMP
#pragma omp parallel for
#endif
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
            if (blocked) continue;
            double d2 = (ax - mx) * (ax - mx) + (ay - my) * (ay - my);
            double taste = attendees[i].tastes[t];
#ifdef _OPENMP
#pragma omp critical(crit_sct)
#endif
            {
                score += (int64_t)ceil(1e6 * taste / d2);
            }
        }
    }

    return score;

}
#endif
