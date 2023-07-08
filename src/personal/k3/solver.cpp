#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING

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



Solution create_trivial_solution(const Problem& problem) {

    double width = problem.stage_w;
    double height = problem.stage_h;
    int ncols = (int)floor((width - 10.0) / 10.0);
    int nrows = (int)floor((height - 10.0) / 10.0);

    Solution solution;

    for (int row = 0; row < nrows; row++) {
        double y = problem.stage_y + row * 10.0 + 10.0;
        for (int col = 0; col < ncols; col++) {
            int id = row * ncols + col;
            if (id >= problem.musicians.size()) continue;
            double x = problem.stage_x + col * 10.0 + 10.0;
            solution.placements.emplace_back(x, y);
        }
    }

    return solution;

}



void solve(int problem_id) {

    Timer timer;

    std::string in_file = format("../data/problems/problem-%d.json", problem_id);
    //std::string sol_file = format("/home/komori3/dev/icfpc2023/data/solutions/k3_v03_climbing/solution-8_1447259902.json");
    std::string out_file_format = "../data/solutions/k3_v03_climbing/solution-%d_%lld.json";
    nlohmann::json data;
    {
        std::ifstream ifs(in_file);
        ifs >> data;
    }

    Problem problem(data);
    DUMP(problem.musicians.size(), problem.attendees.size());

#if defined(_PPL_H) || defined(_OPENMP)
    constexpr int concurrency_coeff = 1;
#else
    constexpr int concurrency_coeff = 10;
#endif
    constexpr int timelimit_phase1 = 10000 * concurrency_coeff;
    constexpr int timelimit_phase2 = 10000000 * concurrency_coeff;

    DUMP(problem_id, timelimit_phase1, timelimit_phase2, concurrency_coeff);

    Xorshift rnd;
    auto solution = create_trivial_solution(problem);
    Solution best_solution;
    double best_score = -1e20;
    //auto best_solution = Solution::from_file(sol_file);
    //double best_score = compute_score_fast(problem, best_solution);
    int loop = 0;

    auto save = [&](const Solution& sol, int problem_id, int64_t score) {
        std::string out_file = format(out_file_format, problem_id, score);
        std::ofstream ofs(out_file);
        ofs << sol.to_json().dump(4);
        DUMP(loop, best_score, timer.elapsed_ms());
    };

    double save_interval = 5000.0;
    double next_save_time = save_interval;

    while (timer.elapsed_ms() < timelimit_phase1 * concurrency_coeff) {
        loop++;
        auto solution_opt = create_random_solution(problem, rnd);
        if (solution_opt) {
            auto solution = solution_opt.value();
            double score = compute_score_fast(problem, solution);
            if (chmax(best_score, score)) {
                DUMP(loop, best_score, timer.elapsed_ms());
                best_solution = solution;
            }
        }
        if (next_save_time < timer.elapsed_ms()) {
            save(best_solution, problem_id, best_score);
            next_save_time += save_interval;
        }
    }
    DUMP(loop);

    while (timer.elapsed_ms() < timelimit_phase2 * concurrency_coeff) {

        auto solution = best_solution;
        int num_musicians = solution.placements.size();

        if (!rnd.next_int(2)) {
            int i, j;
            do {
                i = rnd.next_int(num_musicians);
                j = rnd.next_int(num_musicians);
            } while (i == j);
            if (problem.musicians[i] == problem.musicians[j]) continue;
            loop++;
            std::swap(solution.placements[i].x, solution.placements[j].x);
            std::swap(solution.placements[i].y, solution.placements[j].y);
            double score = compute_score_fast(problem, solution);
            if (chmax(best_score, score)) {
                best_solution = solution;
            }
        }
        else {
            int i = rnd.next_int(num_musicians);
            Placement prev_placement = solution.placements[i];
            Placement placement {
                problem.stage_x + k_musician_spacing_radius + rnd.next_double() * (problem.stage_w - k_musician_spacing_radius * 2),
                problem.stage_y + k_musician_spacing_radius + rnd.next_double() * (problem.stage_h - k_musician_spacing_radius * 2)
                };
            if (!is_musician_on_stage(problem, placement)) continue;
            bool conflict = false;
            for (int kk = 0; kk < solution.placements.size(); ++kk) {
                if (i != kk) {
                    if (are_musicians_too_close(solution.placements[kk], placement)) {
                        conflict = true;
                        break;
                    }
                }
            }
            if (conflict) continue;
            loop++;
            solution.placements[i] = placement;
            double score = compute_score_fast(problem, solution);
            if (chmax(best_score, score)) {
                best_solution = solution;
            }
            else {
                solution.placements[i] = prev_placement;
            }
        }

        if (next_save_time < timer.elapsed_ms()) {
            save(best_solution, problem_id, best_score);
            next_save_time += save_interval;
        }
    }
    DUMP(loop, best_score, compute_score(problem, best_solution));

}


int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    int problem_id = 1;
    if (argc == 2) {
        problem_id = std::stoi(argv[1]);
    }
    DUMP(problem_id);

    solve(problem_id);

    return 0;
}