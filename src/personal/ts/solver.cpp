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

struct Changeset {
    std::array<int, 2> swapped_placements;
    void undo(Solution& solution) const {
        const int i = swapped_placements[0];
        const int j = swapped_placements[1];
        std::swap(solution.placements[i], solution.placements[j]);
    }
};
Changeset do_random_mutation(const Problem& problem, Xorshift& rnd, Solution& solution) {
    const size_t N = solution.placements.size();
    const int i = rnd.next_int() % N;
    const int j = (i + 1 + rnd.next_int() % (N - 1)) % N;

    std::swap(solution.placements[i], solution.placements[j]);

    Changeset changeset;
    changeset.swapped_placements[0] = i;
    changeset.swapped_placements[1] = j;
    return changeset;
}

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

    std::string solution_path;
    if (argc == 2) {
        solution_path = argv[1];
    }
    int problem_id = *guess_problem_id(solution_path);
    DUMP(problem_id);

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    Problem problem = Problem::from_file(problem_id);
    Solution input_solution = Solution::from_file(solution_path);
    double input_score = compute_score(problem, input_solution);
    DUMP(input_score);

    Xorshift rnd;
    Solution best_solution = input_solution;
    Solution current_solution = input_solution;
    double best_score = input_score;
    int loop = 0;
    while (timer.elapsed_ms() < 10000) {
        loop++;
        auto changeset = do_random_mutation(problem, rnd, best_solution);
        double score = compute_score(problem, best_solution);
        if (chmax(best_score, score)) {
            DUMP(loop, best_score);
        } else {
            changeset.undo(best_solution);
        }
    }
    DUMP(loop);
    DUMP(best_score);

    if (best_score > 0) {
        std::ofstream ofs(format("../personal/ts/solutions/solution-%d.json", problem_id));
        ofs << best_solution.to_json().dump(4);
    }

    return 0;
}