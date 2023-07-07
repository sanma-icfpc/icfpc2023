#include "../../../src/k3common.h"
#include "../../../src/spec.h"
#include "../../../src/util.h"

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

    for (int problem_id = 4; problem_id <= 45; problem_id++) {
        solve(problem_id);
    }

    return 0;
}