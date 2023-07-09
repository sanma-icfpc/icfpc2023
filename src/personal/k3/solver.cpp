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

struct State {

    using Pos = Placement;

    static constexpr int MM = 1536; // max musician size: 1484 (problem 33)

    Problem m_problem;
    const int m_num_musicians;
    const int m_num_attendees;

    std::vector<std::vector<std::bitset<MM>>> m_block;
    std::vector<std::vector<int>> m_block_ctr;

    int64_t m_score;

    std::bitset<MM> m_placed;
    std::array<Pos, MM> m_pos;

    // block[sk][i][tk]: sk と i を結ぶ線分と tk 中心の円が交差するか？
        // sk から i が見える条件: block[sk][i][:] が false

        // k を削除
        // 1. for sk in K\{k}: for i in I: block[sk][i][k] = false
        // 2. block[k][:][:] = false

        // k を追加
        // 1. for sk in K\{k}: for i in I: block[sk][i][k] = (交差判定)
        // 2. for i in I: for tk in K\{k}: block[k][i][tk] = (交差判定)

        // block_ctr[sk][i]: block[sk][i][:] のうち true であるものの個数
        // block_ctr[sk][i] が 0 とそれ以外を行き来する際に gain の変化が生じるはず

    State(const Problem& problem)
        : m_problem(problem)
        , m_num_musicians(problem.musicians.size())
        , m_num_attendees(problem.attendees.size())
        , m_block(m_num_musicians, std::vector<std::bitset<MM>>(m_num_attendees))
        , m_block_ctr(m_num_musicians, std::vector<int>(m_num_attendees))
        , m_score(0) {

        m_placed.reset();
        std::fill(m_pos.begin(), m_pos.begin() + m_num_musicians, Pos(-1e9, -1e9));
    }

    int64_t full_compute(const Solution& solution) {
        const auto& as = m_problem.attendees;
        auto pos = solution.placements;
        std::copy(pos.begin(), pos.end(), m_pos.begin());
        std::vector<int64_t> score_k(m_num_musicians);
#pragma omp parallel for
        for (int sk = 0; sk < m_num_musicians; sk++) {
            m_placed[sk] = true;
            int st = m_problem.musicians[sk];
            for (int i = 0; i < m_num_attendees; i++) {
                for (int tk = 0; tk < m_num_musicians; tk++) {
                    if (sk == tk || !is_intersect(m_pos[tk], k_musician_radius, m_pos[sk], as[i])) continue;
                    m_block[sk][i][tk] = true;
                    m_block_ctr[sk][i]++;
                }
                if (!m_block_ctr[sk][i]) {
                    double d2 = distance_squared(m_pos[sk], as[i]);
                    double taste = as[i].tastes[st];
                    score_k[sk] += (int64_t)ceil(1e6 * taste / d2);
                }
            }
        }
        return m_score = std::accumulate(score_k.begin(), score_k.end(), 0LL);
    }

    int64_t push(int k, const Pos& kpos) {

        const auto& ms = m_problem.musicians;
        const auto& as = m_problem.attendees;

        std::vector<int64_t> gain_k(m_num_musicians);
        std::vector<int64_t> gain_i(m_num_attendees);

#pragma omp parallel for
        for (int sk = 0; sk < m_num_musicians; sk++) {
            if (!m_placed[sk]) continue;
            for (int i = 0; i < m_num_attendees; i++) {
                bool blocked = is_intersect(kpos, k_musician_radius, m_pos[sk], as[i]);
                if (blocked) {
                    m_block[sk][i][k] = true;
                    if (!m_block_ctr[sk][i]) {
                        auto taste = as[i].tastes[ms[sk]];
                        auto d2 = distance_squared(m_pos[sk], as[i]);
                        gain_k[sk] -= (int64_t)ceil(1e6 * taste / d2);
                    }
                    m_block_ctr[sk][i]++;
                }
            }
        }

#pragma omp parallel for
        for (int i = 0; i < m_num_attendees; i++) {
            for (int tk = 0; tk < m_num_musicians; tk++) {
                if (!m_placed[tk]) continue;
                bool blocked = is_intersect(m_pos[tk], k_musician_radius, kpos, as[i]);
                if (blocked) {
                    m_block[k][i][tk] = true;
                    m_block_ctr[k][i]++;
                }
            }
            if (!m_block_ctr[k][i]) {
                auto taste = as[i].tastes[ms[k]];
                auto d2 = distance_squared(kpos, as[i]);
                gain_i[i] += (int64_t)ceil(1e6 * taste / d2);
            }
        }

        auto gain
            = std::accumulate(gain_k.begin(), gain_k.end(), 0LL)
            + std::accumulate(gain_i.begin(), gain_i.end(), 0LL);

        m_placed[k] = true;
        m_pos[k] = kpos;
        m_score += gain;

        return gain;
    }

    int64_t pop(int k) {
        
        const auto& ms = m_problem.musicians;
        const auto& as = m_problem.attendees;

        std::vector<int64_t> gain_k(m_num_musicians);
        std::vector<int64_t> gain_i(m_num_attendees);

#pragma omp parallel for
        for (int sk = 0; sk < m_num_musicians; sk++) {
            if (!m_placed[sk]) continue;
            for (int i = 0; i < m_num_attendees; i++) {
                if (m_block[sk][i][k]) {
                    m_block[sk][i][k] = false;
                    m_block_ctr[sk][i]--;
                    if (!m_block_ctr[sk][i]) {
                        auto taste = as[i].tastes[ms[sk]];
                        auto d2 = distance_squared(m_pos[sk], as[i]);
                        gain_k[sk] += (int64_t)ceil(1e6 * taste / d2);
                    }
                }
            }
        }

#pragma omp parallel for
        for (int i = 0; i < m_num_attendees; i++) {
            m_block[k][i].reset();
            if (!m_block_ctr[k][i]) {
                auto taste = as[i].tastes[ms[k]];
                auto d2 = distance_squared(m_pos[k], as[i]);
                gain_i[i] -= (int64_t)ceil(1e6 * taste / d2);
            }
            m_block_ctr[k][i] = 0;
        }

        auto gain
            = std::accumulate(gain_k.begin(), gain_k.end(), 0LL)
            + std::accumulate(gain_i.begin(), gain_i.end(), 0LL);
        
        m_placed[k] = false;
        m_pos[k] = { -1e9, -1e9 };
        m_score += gain;

        return gain;
    }

    bool can_move(int k, const Pos& kpos) const {
        if (!is_musician_on_stage(m_problem, kpos)) return false;
        for (int kk = 0; kk < m_num_musicians; kk++) {
            if (!m_placed[kk] || k == kk) continue;
            if (are_musicians_too_close(kpos, m_pos[kk], 0)) return false;
        }
        return true;
    }

    int64_t move(int k, const Pos& kpos) {
        auto gain = pop(k);
        gain += push(k, kpos);
        return gain;
    }

    int64_t swap(int k1, int k2) {
        auto p1 = m_pos[k1], p2 = m_pos[k2];
        auto gain = pop(k1);
        gain += pop(k2);
        gain += push(k1, p2);
        gain += push(k2, p1);
        return gain;
    }

    Pos sample_random_pos(Xorshift& rnd) const {
        return {
            m_problem.stage_x + k_musician_spacing_radius + rnd.next_double() * (m_problem.stage_w - k_musician_spacing_radius * 2),
            m_problem.stage_y + k_musician_spacing_radius + rnd.next_double() * (m_problem.stage_h - k_musician_spacing_radius * 2)
        };
    }

    Solution to_solution() const {
        Solution sol;
        sol.placements = std::vector<Pos>(m_pos.begin(), m_pos.begin() + m_num_musicians);
        return sol;
    }

};

#if 1
inline double get_temp(double stemp, double etemp, double t, double T) {
    return etemp + (stemp - etemp) * (T - t) / T;
};
#else
inline double get_temp(double stemp, double etemp, double t, double T) {
    return stemp * pow(etemp / stemp, t / T);
};
#endif

void solve_anneal(int problem_id) {
    Timer timer;

    std::string in_file = format("../data/problems/problem-%d.json", problem_id);
    //std::string sol_file = format("../data/solutions/bests/solution-%d.json", problem_id);
    std::string sol_file;
    std::string out_file_format = "../data/solutions/bests/solution-%d.json";
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

    auto sol = *create_random_solution(problem, rnd);
    if (!sol_file.empty()) {
        sol = Solution::from_file(sol_file);
    }
    else {
        int loop_random = 0;
        auto best_score = compute_score(problem, sol);
        while (timer.elapsed_ms() < 3000) {
            loop_random++;
            auto nsol = create_random_solution(problem, rnd);
            if (!nsol) continue;
            auto score = compute_score(problem, *nsol);
            if (chmax(best_score, score)) {
                sol = *nsol;
                DUMP(loop_random, best_score);
            }
        }
    }

    CachedComputeScore cache(problem);
    cache.full_compute(sol);

    int loop = 0;
    double dump_interval = 1000.0;
    double save_interval = 10000.0;
    double start_time = timer.elapsed_ms(), now_time = start_time, end_time = 60000;
    double next_dump_time = now_time;
    double next_save_time = now_time;
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

void solve(int problem_id) {
    Timer timer;

    std::string in_file = format("../data/problems/problem-%d.json", problem_id);
    std::string sol_file = format("../data/solutions/bests/solution-%d.json", problem_id);

    if (!std::filesystem::exists(sol_file)) return;

    std::string out_file_format = "../data/solutions/k3_v05_opt_volumes/solution-%d_sub=%lld.json";
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

    auto sol = *create_random_solution(problem, rnd);
    if (!sol_file.empty()) {
        sol = Solution::from_file(sol_file);
    }
    else {
        int loop_random = 0;
        auto best_score = compute_score(problem, sol);
        while (timer.elapsed_ms() < 3000) {
            loop_random++;
            auto nsol = create_random_solution(problem, rnd);
            if (!nsol) continue;
            auto score = compute_score(problem, *nsol);
            if (chmax(best_score, score)) {
                sol = *nsol;
                DUMP(loop_random, best_score);
            }
        }
    }

    auto prev_score = compute_score(problem, sol);
    set_optimal_volumes(problem, sol);
    auto curr_score = compute_score(problem, sol);
    auto diff = curr_score - prev_score;
    LOG(INFO) << format("problem_id=%2d, prev_score=%12lld, curr_score=%12lld, diff=%12lld", problem_id, prev_score, curr_score, diff);
}


int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    solve_anneal(16);

    return 0;
}