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

std::optional<Solution> create_random_solution(const Problem& problem, Xorshift& rnd, double timelimit = 1000) {

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

void solve(int problem_id) {

    Timer timer;

    std::ifstream ifs(format("../data/problems/problem-%d.json", problem_id));
    nlohmann::json data;
    ifs >> data;

    Problem problem(data);

#ifdef _PPL_H
    constexpr int concurrency_coeff = 1;
#else
    constexpr int concurrency_coeff = 10;
#endif
    constexpr int timelimit_phase1 = 10000 * concurrency_coeff;
    constexpr int timelimit_phase2 = 60000 * concurrency_coeff;

    DUMP(problem_id, timelimit_phase1, timelimit_phase2, concurrency_coeff);

    Xorshift rnd;
    //auto solution = create_trivial_solution(problem);
    Solution best_solution;
    double best_score = -1e20;
    int loop = 0;

    while (timer.elapsed_ms() < timelimit_phase1 * concurrency_coeff) {
        loop++;
        auto solution_opt = create_random_solution(problem, rnd);
        if (solution_opt) {
            auto solution = solution_opt.value();
            double score = compute_score(problem, solution);
            if (chmax(best_score, score)) {
                DUMP(loop, best_score, timer.elapsed_ms());
                best_solution = solution;
            }
        }
    }
    DUMP(loop);

    while (timer.elapsed_ms() < timelimit_phase2 * concurrency_coeff) {
        auto solution = best_solution;
        int num_musicians = solution.placements.size();
        int i, j;
        do {
            i = rnd.next_int(num_musicians);
            j = rnd.next_int(num_musicians);
        } while (i == j);
        if (problem.musicians[i] == problem.musicians[j]) continue;
        loop++;
        std::swap(solution.placements[i].x, solution.placements[j].x);
        std::swap(solution.placements[i].y, solution.placements[j].y);
        double score = compute_score(problem, solution);
        if (chmax(best_score, score)) {
            DUMP(loop, best_score, timer.elapsed_ms());
            best_solution = solution;
        }
    }
    DUMP(loop);

    if (best_score > 0) {
        std::ofstream ofs(format("../data/solutions/k3_v02_k3_v02_random_swap_after_creation/solution-%d.json", problem_id));
        ofs << best_solution.to_json().dump(4);
    }
}


int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    solve(1);

    return 0;
}