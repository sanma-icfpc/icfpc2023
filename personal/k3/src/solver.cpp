#define _CRT_NONSTDC_NO_WARNINGS
#define _SILENCE_CXX17_ITERATOR_BASE_CLASS_DEPRECATION_WARNING
#include <bits/stdc++.h>
#include <random>
#include <unordered_set>
#include <array>
#include <optional>
#include <regex>
#include <json.hpp>
#ifdef _MSC_VER
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#include <conio.h>
#include <ppl.h>
#include <filesystem>
#include <intrin.h>
/* g++ functions */
int __builtin_clz(unsigned int n) { unsigned long index; _BitScanReverse(&index, n); return 31 - index; }
int __builtin_ctz(unsigned int n) { unsigned long index; _BitScanForward(&index, n); return index; }
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
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x); // container out (fwd decl)
template<class S, class T> std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) { return os << '{' << p.first << ", " << p.second << '}'; } // pair out
template<class S, class T> std::istream& operator>>(std::istream& is, std::pair<S, T>& p) { return is >> p.first >> p.second; } // pair in
std::ostream& operator<<(std::ostream& os, const std::vector<bool>::reference& v) { os << (v ? '1' : '0'); return os; } // bool (vector) out
std::ostream& operator<<(std::ostream& os, const std::vector<bool>& v) { bool f = true; os << '{'; for (const auto& x : v) { os << (f ? "" : ", ") << x; f = false; } os << '}'; return os; } // vector<bool> out
template<class Ch, class Tr, class Container> std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) { bool f = true; os << '{'; for (auto& y : x) { os << (f ? "" : ", ") << y; f = false; } return os << '}'; } // container out
template<class T, class = decltype(std::begin(std::declval<T&>())), class = typename std::enable_if<!std::is_same<T, std::string>::value>::type> std::istream& operator>>(std::istream& is, T& a) { for (auto& x : a) is >> x; return is; } // container in
template<typename T> auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) { out << t.stringify(); return out; } // struct (has stringify() func) out
/** io setup **/
struct IOSetup { IOSetup(bool f) { if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); } std::cout << std::fixed << std::setprecision(15); } }
iosetup(true); // set false when solving interective problems
/** string formatter **/
template<typename... Ts> std::string format(const std::string& f, Ts... t) { size_t l = std::snprintf(nullptr, 0, f.c_str(), t...); std::vector<char> b(l + 1); std::snprintf(&b[0], l + 1, f.c_str(), t...); return std::string(&b[0], &b[0] + l); }
/** dump **/
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define DUMP(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<']'<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
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
/** shuffle **/
template<typename T> void shuffle_vector(std::vector<T>& v, Xorshift& rnd) { int n = v.size(); for (int i = n - 1; i >= 1; i--) { int r = rnd.next_int(i); std::swap(v[i], v[r]); } }
/** split **/
std::vector<std::string> split(const std::string& str, const std::string& delim) {
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
template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }
/** using **/
//using ll = long long;
//using ld = long double;
////using ld = boost::multiprecision::cpp_bin_float_quad;
//using pii = std::pair<int, int>;
//using pll = std::pair<ll, ll>;
//using std::cin, std::cout, std::cerr, std::endl, std::string, std::vector, std::array;



struct Attendee {

    double x;
    double y;
    std::vector<double> tastes;

    std::string stringify() const {
        std::string res = format("Attendee {x=%f, y=%f, tastes=[", x, y);
        res += std::to_string(tastes[0]);
        for (int i = 1; i < tastes.size(); i++) {
            res += ", " + std::to_string(tastes[i]);
        }
        res += "]}";
        return res;
    }

};

struct Problem {

    double room_width;
    double room_height;
    double stage_w;
    double stage_h;
    double stage_x;
    double stage_y;
    std::vector<int> musicians;
    std::vector<Attendee> attendees;

    Problem(const nlohmann::json& data) {
        room_width = data["room_width"];
        room_height = data["room_height"];
        stage_w = data["stage_width"];
        stage_h = data["stage_height"];
        stage_x = data["stage_bottom_left"][0];
        stage_y = data["stage_bottom_left"][1];
        for (const auto& type : data["musicians"]) {
            musicians.push_back(type);
        }
        for (const auto& attendee_json : data["attendees"]) {
            Attendee attendee;
            attendee.x = attendee_json["x"];
            attendee.y = attendee_json["y"];
            for (const auto& taste : attendee_json["tastes"]) {
                attendee.tastes.push_back(taste);
            }
            attendees.push_back(attendee);
        }
    }

#ifdef HAVE_OPENCV_HIGHGUI
    cv::Mat_<cv::Vec3b> show(int delay = 0) const {
        int H = (int)room_height;
        int W = (int)room_width;
        cv::Mat_<cv::Vec3b> img(H, W, cv::Vec3b(155, 155, 155));
        cv::Rect stage(
            (int)stage_x,  // x
            (int)stage_y,    // y
            (int)stage_w,   // w
            (int)stage_h   // h
        );
        cv::rectangle(img, stage, cv::Scalar(255, 255, 255), cv::FILLED);
        for (const auto& a : attendees) {
            int y = (int)a.y;
            int x = (int)a.x;
            cv::circle(img, cv::Point(x, y), 2, cv::Scalar(0, 0, 255), cv::FILLED);
        }
        cv::imshow("img", img);
        cv::waitKey(delay);
    }
#endif

};

struct Placement {

    double x;
    double y;

    Placement(double x = 0.0, double y = 0.0) : x(x), y(y) {}

    std::string stringify() const { return format("(%f, %f)", x, y); }

};

struct Solution {

    std::vector<Placement> placements;

    std::string stringify() const {
        std::ostringstream oss;
        oss << placements;
        return oss.str();
    }

    nlohmann::json to_json() const {
        nlohmann::json data;
        data["placements"] = {};
        for (const auto& [x, y] : placements) {
            data["placements"].push_back({ {"x", x}, {"y", y} });
        }
        return data;
    }

};

Solution create_trivial_solution(const Problem& problem) {

    double width = problem.stage_w;
    double height = problem.stage_h;
    int ncols = (int)floor((width - 10.0) / 10.0);
    int nrows = (int)floor((height - 10.0) / 10.0);
    DUMP(ncols, nrows);

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

    DUMP(solution.placements);

    return solution;

}

double is_intersect(double cx, double cy, double r, double x1, double y1, double x2, double y2) {
    x1 -= cx; x2 -= cx; y1 -= cy; y2 -= cy;
    double dx = x2 - x1, dy = y2 - y1;
    double a = dx * dx + dy * dy, b = x1 * dx + y1 * dy, c = x1 * x1 + y1 * y1 - r * r;
    if (b * b - a * c < 0) return false;
    double f0 = c, f1 = a * 2 * b + c;
    if (f0 * f1 <= 0) return true;
    return 0 <= f0 && 0 <= f1 && -a <= b && b <= 0 && a * c <= b * b;
}

int64_t compute_score(const Problem& problem, const Solution& solution) {

    const auto& musicians = problem.musicians;
    const auto& attendees = problem.attendees;
    const auto& placements = solution.placements;

    int64_t score = 0;

    for (int k = 0; k < musicians.size(); k++) {
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
            score += (int64_t)ceil(1e6 * taste / d2);
        }
    }

    return score;

}

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    std::ifstream ifs("../../in/problem-10.json");
    nlohmann::json data;
    ifs >> data;

    Problem problem(data);

#ifdef HAVE_OPENCV_HIGHGUI
    problem.show();
#endif

    auto solution = create_trivial_solution(problem);

    DUMP(compute_score(problem, solution));
    // local:  -154693724
    // remote: -154693721

    std::ofstream ofs("../../out/problem-10.json");
    ofs << solution.to_json().dump(4);

    return 0;
}