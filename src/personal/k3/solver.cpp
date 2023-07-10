#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING

#include "../../cached_compute_score.h"
#include "../../spec.h"
#include "../../util.h"

#ifndef _MSC_VER
#include <bits/stdc++.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <random>
#include <unordered_set>
#include <array>
#include <bitset>
#include <algorithm>
#include <optional>
#include <regex>
#include <nlohmann/json.hpp>
#ifdef USE_OPENCV
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#pragma GCC diagnostic pop
#endif
#ifdef _MSC_VER
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
/* g++ functions */
inline int __builtin_clz(unsigned int n) { unsigned long index; _BitScanReverse(&index, n); return 31 - index; }
inline int __builtin_ctz(unsigned int n) { unsigned long index; _BitScanForward(&index, n); return index; }
namespace std { inline int __lg(int __n) { return sizeof(int) * 8 - 1 - __builtin_clz(__n); } }
/* enable __uint128_t in MSVC */
//#include <boost/multiprecision/cpp_int.hpp>
//using __uint128_t = boost::multiprecision::uint128_t;
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro io **/
namespace aux {
    template<typename T, unsigned N, unsigned L> struct tp { static void output(std::ostream& os, const T& v) { os << std::get<N>(v) << ", "; tp<T, N + 1, L>::output(os, v); } };
    template<typename T, unsigned N> struct tp<T, N, N> { static void output(std::ostream& os, const T& v) { os << std::get<N>(v); } };
}
template<typename... Ts> std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) { os << '{'; aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t); return os << '}'; } // tuple out
template<class Ch, class Tr, class Container, class = decltype(std::begin(std::declval<Container&>()))> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x); // container out (fwd decl)
template<class S, class T> std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) { return os << '{' << p.first << ", " << p.second << '}'; } // pair out
template<class S, class T> std::istream& operator>>(std::istream& is, std::pair<S, T>& p) { return is >> p.first >> p.second; } // pair in
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>::reference& v) { os << (v ? '1' : '0'); return os; } // bool (vector) out
inline std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) { bool f = true; os << '{'; for (const auto& x : v) { os << (f ? "" : ", ") << x; f = false; } os << '}'; return os; } // vector<bool> out
template<class Ch, class Tr, class Container, class> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) { bool f = true; os << '{'; for (auto& y : x) { os << (f ? "" : ", ") << y; f = false; } return os << '}'; } // container out
template<class T, class = decltype(std::begin(std::declval<T&>())), class = typename std::enable_if<!std::is_same<T, std::string>::value>::type> std::istream& operator>>(std::istream& is, T& a) { for (auto& x : a) is >> x; return is; } // container in
template<typename T> auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) { out << t.stringify(); return out; } // struct (has stringify() func) out
/** io setup **/
struct IOSetup { IOSetup(bool f) { if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); } std::cout << std::fixed << std::setprecision(15); } }
iosetup(true); // set false when solving interective problems
/** shuffle **/
template<typename T> void shuffle_vector(std::vector<T>& v, Xorshift& rnd) { int n = v.size(); for (int i = n - 1; i >= 1; i--) { int r = rnd.next_int(i); std::swap(v[i], v[r]); } }
/** split **/
inline std::vector<std::string> split(const std::string& str, const std::string& delim) {
    std::vector<std::string> res;
    std::string buf;
    for (const auto& c : str) {
        if (delim.find(c) != std::string::npos) {
            if (!buf.empty()) res.push_back(buf);
            buf.clear();
        }
        else buf += c;
    }
    if (!buf.empty()) res.push_back(buf);
    return res;
}
/** misc **/
template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) { std::fill((T*)array, (T*)(array + N), val); } // fill array
template<typename T, typename ...Args> auto make_vector(T x, int arg, Args ...args) { if constexpr (sizeof...(args) == 0)return std::vector<T>(arg, x); else return std::vector(arg, make_vector<T>(x, args...)); }



#if 1
inline double get_temp(double stemp, double etemp, double t, double T) {
    return etemp + (stemp - etemp) * (T - t) / T;
};
#else
inline double get_temp(double stemp, double etemp, double t, double T) {
    return stemp * pow(etemp / stemp, t / T);
};
#endif


struct CachedSwapComputeScore {
  public:
    const Problem &m_problem;
    const int m_num_attendees;
    const int m_num_musicians;
    const int m_num_pillars;
    Solution m_solution;

    std::vector<double> m_harmony_cache; // 1 + m_harmony_cache[k] is the harmony.
    std::vector<int64_t> m_score_cache;
    std::vector<uint8_t> m_audible_cache;
    int64_t m_score;

    double m_accum_elapsed_ms_partial = 0.0;
    int m_call_count_partial = 0;
    double m_accum_elapsed_ms_full = 0.0;
    int m_call_count_full = 0;

    double get_mean_elapsed_ms_partial() const {
        return m_accum_elapsed_ms_partial / (m_call_count_partial + 1e-9);
    }
    double get_mean_elapsed_ms_full() const {
        return m_accum_elapsed_ms_full / (m_call_count_full + 1e-9);
    }
    void report() const {
        LOG(INFO) << format("PID=%d: full calculation %.4f ms/eval (%d evals), "
                            "partial update %.4f "
                            "ms/eval (%d evals)",
                            m_problem.problem_id, get_mean_elapsed_ms_full(),
                            m_call_count_full, get_mean_elapsed_ms_partial(),
                            m_call_count_partial);
    }

  public:
    CachedSwapComputeScore(const Problem &problem)
    : m_problem(problem),
      m_num_attendees(problem.attendees.size()),
      m_num_musicians(problem.musicians.size()),
      m_num_pillars(problem.extension.consider_pillars ? problem.pillars.size() : 0),
      m_harmony_cache(m_num_musicians, 0.0),
      m_score_cache(problem.attendees.size() * problem.musicians.size(), 0),
      m_audible_cache(problem.attendees.size() * problem.musicians.size(), 1) {}

    int64_t &partial_score(int k, int i) {
#if 0
        LOG_ASSERT(0 <= k && k < m_num_musicians);
        LOG_ASSERT(0 <= i && i < m_num_attendees);
#endif
        return m_score_cache[k * m_num_attendees + i];
    }

    // k_src's ray to i is audible by k_other
    uint8_t &is_audible(int k, int i) {
        return m_audible_cache[k * m_num_attendees + i];
    }
    int64_t score() const { return m_score; }

    int64_t swap_two_musicians(int k1, int k2, bool dry_run) {
        Timer timer;

        // swap two musicians position
        std::swap(m_solution.placements[k1], m_solution.placements[k2]);

        const auto &musicians = m_problem.musicians;
        const auto &attendees = m_problem.attendees;
        const auto &placements = m_solution.placements;
        const auto &pillars = m_problem.pillars;
        const int M = musicians.size();
        const int A = attendees.size();
        const int P = pillars.size();
        LOG_ASSERT(M == placements.size());

        auto harmony_cache_backup = m_harmony_cache;

        int64_t influence_diff = 0;
        
        // k1 について
        // position と audible list は k2 のものと swap
        // partial score は 再計算 (距離をキャッシュ？)
        // 

        m_accum_elapsed_ms_partial += timer.elapsed_ms();
        m_call_count_partial += 1;

        if(dry_run) {
            std::swap(m_solution.placements[k1], m_solution.placements[k2]);
            m_harmony_cache = harmony_cache_backup;
        } else {
            m_score += influence_diff;
        }

        return influence_diff;
    }

    int64_t full_compute(const Solution &solution) {
        Timer timer;

        m_solution = solution;
        const auto &musicians = m_problem.musicians;
        const auto &attendees = m_problem.attendees;
        const auto &placements = solution.placements;
        const auto &pillars = m_problem.pillars;
        const int M = musicians.size();
        const int A = attendees.size();
        const int P = pillars.size();

        m_score = 0;
        m_score_cache.assign(m_score_cache.size(), 0);
        m_audible_cache.assign(m_audible_cache.size(), 1);
        m_harmony_cache.assign(m_harmony_cache.size(), 0.0);

        if(m_problem.extension.consider_harmony) {
    #pragma omp parallel for
            for(int k = 0; k < M; k++) {
                double harmony = 0.0;
                for(int k_other = 0; k_other < M; ++k_other) {
                    if(k != k_other && musicians[k] == musicians[k_other]) {
                        harmony +=
                            inverse_distance(placements[k], placements[k_other]);
                    }
                }
                m_harmony_cache[k] = harmony;
            }
        }

        int64_t score_diff = 0;
    #pragma omp parallel for reduction(+ : score_diff)
        for(auto k_src = 0; k_src < M; ++k_src) {
            for(int i = 0; i < A; ++i) {
                bool audible = true;
                for(int k_other = 0; k_other < M; ++k_other) {
                    if(
                        k_src != k_other &&
                        is_intersect(placements[k_other], k_musician_radius, placements[k_src], attendees[i])
                        ) {
                        audible = false;        
                        break;                
                    }
                }
                if (audible) {
                    for(int p = 0; p < P; ++p) {
                        if(is_intersect(pillars[p], pillars[p].r, placements[k_src], attendees[i])) {
                            audible = false;
                            break;
                        }
                    }
                }
                const int64_t influence = (int64_t)ceil(
                    1e6 * attendees[i].tastes[musicians[k_src]] /
                    distance_squared(placements[k_src], attendees[i]));
                partial_score(k_src, i) = influence;
                is_audible(k_src, i) = audible;
                score_diff +=
                    audible
                        ? (int64_t)ceil(m_solution.volumes[k_src] *
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

struct CachedSwapComputeScoreWithoutHarmony {
  public:
    const Problem &m_problem;
    const int m_num_attendees;
    const int m_num_musicians;
    const int m_num_pillars;
    Solution m_solution;

    std::vector<int64_t> m_score_cache;
    std::vector<uint8_t> m_audible_cache;
    int64_t m_score;

  public:
    CachedSwapComputeScoreWithoutHarmony(const Problem &problem)
    : m_problem(problem),
      m_num_attendees(problem.attendees.size()),
      m_num_musicians(problem.musicians.size()),
      m_num_pillars(problem.extension.consider_pillars ? problem.pillars.size() : 0),
      m_score_cache(problem.attendees.size() * problem.musicians.size(), 0),
      m_audible_cache(problem.attendees.size() * problem.musicians.size(), 1) {}

    int64_t &partial_score(int k, int i) {
#if 0
        LOG_ASSERT(0 <= k && k < m_num_musicians);
        LOG_ASSERT(0 <= i && i < m_num_attendees);
#endif
        return m_score_cache[k * m_num_attendees + i];
    }

    // k_src's ray to i is audible by k_other
    uint8_t &is_audible(int k, int i) {
        return m_audible_cache[k * m_num_attendees + i];
    }
    int64_t score() const { return m_score; }

    int64_t swap_two_musicians(int k1, int k2, bool dry_run) {
        Timer timer;

        const auto &musicians = m_problem.musicians;
        const auto &attendees = m_problem.attendees;
        const auto &placements = m_solution.placements;
        const auto &pillars = m_problem.pillars;
        const int M = musicians.size();
        const int A = attendees.size();
        const int P = pillars.size();
        LOG_ASSERT(M == placements.size());

        // swap position
        std::swap(m_solution.placements[k1], m_solution.placements[k2]);

        int64_t influence_diff = 0;
        
        // k1 について
        // position と audible list は k2 のものと swap
        // partial score は 再計算 (距離をキャッシュ？)

        for (int i = 0; i < A; i++) {

            auto audible_k1 = is_audible(k1, i);
            int64_t influence_k1 = (int64_t)ceil(
                1e6 * attendees[i].tastes[musicians[k1]] /
                 distance_squared(placements[k1], attendees[i]));

            auto audible_k2 = is_audible(k2, i);
            int64_t influence_k2 = (int64_t)ceil(
                1e6 * attendees[i].tastes[musicians[k2]] /
                 distance_squared(placements[k2], attendees[i]));

            influence_diff += (audible_k2 ? influence_k1 : 0) - (audible_k1 ? partial_score(k1, i) : 0)
                            + (audible_k1 ? influence_k2 : 0) - (audible_k2 ? partial_score(k2, i) : 0);

            if (!dry_run) {
                std::swap(is_audible(k1, i), is_audible(k2, i));
                partial_score(k1, i) = influence_k1;
                partial_score(k2, i) = influence_k2;
            }

        }

        if(dry_run) {
            std::swap(m_solution.placements[k1], m_solution.placements[k2]);
        } else {
            m_score += influence_diff;
        }

        return influence_diff;
    }

    int64_t full_compute(const Solution &solution) {
        Timer timer;

        m_solution = solution;
        const auto &musicians = m_problem.musicians;
        const auto &attendees = m_problem.attendees;
        const auto &placements = solution.placements;
        const auto &pillars = m_problem.pillars;
        const int M = musicians.size();
        const int A = attendees.size();
        const int P = pillars.size();

        m_score = 0;
        m_score_cache.assign(m_score_cache.size(), 0);
        m_audible_cache.assign(m_audible_cache.size(), 1);

        if(m_problem.extension.consider_harmony) {
    #pragma omp parallel for
            for(int k = 0; k < M; k++) {
                double harmony = 0.0;
                for(int k_other = 0; k_other < M; ++k_other) {
                    if(k != k_other && musicians[k] == musicians[k_other]) {
                        harmony +=
                            inverse_distance(placements[k], placements[k_other]);
                    }
                }
            }
        }

        int64_t score_diff = 0;
    #pragma omp parallel for reduction(+ : score_diff)
        for(auto k_src = 0; k_src < M; ++k_src) {
            for(int i = 0; i < A; ++i) {
                bool audible = true;
                for(int k_other = 0; k_other < M; ++k_other) {
                    if(
                        k_src != k_other &&
                        is_intersect(placements[k_other], k_musician_radius, placements[k_src], attendees[i])
                        ) {
                        audible = false;        
                        break;                
                    }
                }
                if (audible) {
                    for(int p = 0; p < P; ++p) {
                        if(is_intersect(pillars[p], pillars[p].r, placements[k_src], attendees[i])) {
                            audible = false;
                            break;
                        }
                    }
                }
                const int64_t influence = (int64_t)ceil(
                    1e6 * attendees[i].tastes[musicians[k_src]] /
                    distance_squared(placements[k_src], attendees[i]));
                partial_score(k_src, i) = influence;
                is_audible(k_src, i) = audible;
                score_diff += audible ? (int64_t)ceil(m_solution.volumes[k_src] * influence) : 0;
            }
        }
        m_score += score_diff;

        return m_score;
    }
};

void save_solution(const Solution& sol, const std::string& out_file_format, int problem_id, int64_t score) {
    std::string out_file = format(out_file_format, problem_id, score);
    std::ofstream ofs(out_file);
    ofs << sol.to_json().dump(4);
};

Solution swap_anneal(int problem_id, const Problem& problem, Solution sol, const std::string out_file_format, double duration) {
    DUMP(problem_id);
    Timer timer;

    Xorshift rnd;

    DUMP(compute_score(problem, sol));
    set_optimal_volumes(problem, sol, 1.0);
    double capped_score = compute_score(problem, sol);
    DUMP(compute_score(problem, sol));
    set_constant_volumes(problem, sol, 1.0);
    DUMP(compute_score(problem, sol));

    CachedSwapComputeScoreWithoutHarmony cache(problem);
    cache.full_compute(sol);

    int loop = 0;
    double dump_interval = 100.0;
    double start_time = timer.elapsed_ms(), now_time = start_time, end_time = start_time + duration;
    double start_temp = std::max(100.0, capped_score * 1e-5);
    DUMP(start_temp);
    double next_dump_time = start_time + dump_interval;
    while ((now_time = timer.elapsed_ms()) < end_time) {
        loop++;
        int k1, k2;
        do {
            k1 = rnd.next_int(problem.musicians.size());
            k2 = rnd.next_int(problem.musicians.size());
        } while(k1 == k2);
        auto gain = cache.swap_two_musicians(k1, k2, true);
        double temp = get_temp(start_temp, 0, now_time - start_time, end_time - start_time);
        double prob = exp(gain / temp);
        if (rnd.next_double() < prob) {
            cache.swap_two_musicians(k1, k2, false);
        }
        if (next_dump_time < now_time) {
            DUMP(now_time, loop, cache.score());
            next_dump_time += dump_interval;
        }
    }

    sol = cache.m_solution;
    auto score = compute_score(problem, sol);
    DUMP(cache.score(), score);
    set_optimal_volumes(problem, sol);
    score = compute_score(problem, sol);
    DUMP(score);

    save_solution(sol, out_file_format, problem_id, score);

    return sol;
}

Solution move_anneal(int problem_id, const Problem& problem, Solution sol, const std::string out_file_format, double duration) {
    DUMP(problem_id);
    Timer timer;

    Xorshift rnd;

    DUMP(compute_score(problem, sol));
    set_optimal_volumes(problem, sol, 1.0);
    double capped_score = compute_score(problem, sol);
    DUMP(compute_score(problem, sol));
    set_constant_volumes(problem, sol, 1.0);
    DUMP(compute_score(problem, sol));

    CachedComputeScore cache(problem);
    cache.full_compute(sol);

    int loop = 0;
    double dump_interval = 1000.0;
    double save_interval = 10000.0;
    double start_time = timer.elapsed_ms(), now_time = start_time, end_time = start_time + duration;
    double start_temp = std::max(100.0, capped_score * 1e-5);
    DUMP(start_temp);
    double next_dump_time = start_time + dump_interval;
    double next_save_time = start_time + save_interval;
    while ((now_time = timer.elapsed_ms()) < end_time) {
        loop++;
        const int i = rnd.next_int(problem.musicians.size());
        auto old_placement = cache.m_solution.placements[i];
        auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
        if (!new_placement_opt) continue;
        auto gain = cache.change_musician(i, *new_placement_opt);
        double temp = get_temp(start_temp, 0, now_time - start_time, end_time - start_time);
        double prob = exp(gain / temp);
        if (rnd.next_double() > prob) {
            cache.change_musician(i, old_placement);
        }
        if (next_dump_time < now_time) {
            DUMP(now_time, loop, cache.score());
            next_dump_time += dump_interval;
        }
    }

    sol = cache.m_solution;
    auto score = compute_score(problem, sol);
    DUMP(cache.score(), score);
    set_optimal_volumes(problem, sol);
    score = compute_score(problem, sol);
    DUMP(score);
    save_solution(sol, out_file_format, problem_id, score);

    return sol;
}

void solve(int problem_id) {

    std::string in_file = format("../data/problems/problem-%d.json", problem_id);
    auto problem = Problem::from_file(in_file);

    std::string sol_file = format("../data/solutions/bests/solution-%d.json", problem_id);
    Solution solution;
    if (!sol_file.empty() && std::filesystem::exists(sol_file)) {
        solution = Solution::from_file(sol_file);
    }
    else {
        Xorshift rnd;
        solution = *create_random_solution(problem, rnd);
    }

    std::string out_file_format = "../data/solutions/k3_v07_alternate_optimization/solution-%d_sub=%lld.json";
    
    for (int trial = 0; trial < 10; trial++) {
        solution = swap_anneal(problem_id, problem, solution, out_file_format, 20000);
        solution = move_anneal(problem_id, problem, solution, out_file_format, 60000);
    }

}


int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    int problem_id = 1;
    if (argc > 1) {
        problem_id = std::stoi(argv[1]);
    }
    solve(problem_id);

    return 0;
}
