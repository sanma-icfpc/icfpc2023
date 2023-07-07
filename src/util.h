#pragma once
#include <regex>
#include <glog/logging.h>
#include "k3common.h"
#include "spec.h"

inline double is_intersect(double cx, double cy, double r, double x1, double y1, double x2, double y2) {
    x1 -= cx; x2 -= cx; y1 -= cy; y2 -= cy;
    double dx = x2 - x1, dy = y2 - y1;
    double a = dx * dx + dy * dy, b = x1 * dx + y1 * dy, c = x1 * x1 + y1 * y1 - r * r;
    if (b * b - a * c < 0) return false;
    double f0 = c, f1 = a * 2 * b + c;
    if (f0 * f1 <= 0) return true;
    return 0 <= f0 && 0 <= f1 && -a <= b && b <= 0 && a * c <= b * b;
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

inline int64_t compute_score(const Problem& problem, const Solution& solution) {

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
