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

void anneal_after_volume_change(int problem_id) {
    DUMP(problem_id);
    Timer timer;

    std::string in_file = format("../data/problems/problem-%d.json", problem_id);
    std::string sol_file = format("../data/solutions/bests/solution-%d.json", problem_id);
    std::string out_file_format = "../data/solutions/k3_v05_anneal_after_volume_change/solution-%d_sub=%lld.json";
    nlohmann::json data;
    {
        std::ifstream ifs(in_file);
        ifs >> data;
    }

    Problem problem(data);

    Xorshift rnd;

    auto save = [&](const Solution& sol, int problem_id, int64_t score) {
        std::string out_file = format(out_file_format, problem_id, score);
        std::ofstream ofs(out_file);
        ofs << sol.to_json().dump(4);
    };

    auto sol = Solution::from_file(sol_file);
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
    double start_time = timer.elapsed_ms(), now_time = start_time, end_time = start_time + 60000;
    double start_temp = std::max(100.0, capped_score * 1e-5);
    DUMP(start_temp);
    double next_dump_time = start_time + dump_interval;
    double next_save_time = start_time + save_interval;
    bool save_mode = false;
    while ((now_time = timer.elapsed_ms()) < end_time) {
        loop++;
        const int i = rnd.next_int(problem.musicians.size());
        auto old_placement = cache.m_solution.placements[i];
        auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
        if (!new_placement_opt) continue;
        auto gain = cache.change_musician(i, *new_placement_opt);
        double temp = get_temp(1e4, 0, now_time - start_time, end_time - start_time);
        double prob = exp(gain / temp);
        if (rnd.next_double() > prob) {
            cache.change_musician(i, old_placement);
        }
        if (next_dump_time < now_time) {
            DUMP(now_time, loop, cache.score());
            next_dump_time += dump_interval;
        }
        if (save_mode && next_save_time < now_time) {
            save(cache.m_solution, problem_id, cache.score());
            next_save_time += save_interval;
        }
    }

    sol = cache.m_solution;
    auto score = compute_score(problem, sol);
    DUMP(cache.score(), score);
    set_optimal_volumes(problem, sol);
    score = compute_score(problem, sol);
    DUMP(score);
    save(sol, problem_id, score);
}



int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    for (int problem_id = 3; problem_id <= 90; problem_id++) {
        if (problem_id == 18) continue;
        anneal_after_volume_change(problem_id);
    }

    return 0;
}
