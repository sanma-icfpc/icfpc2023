#include "../../../src/k3common.h"
#include "../../../src/spec.h"
#include "../../../src/util.h"

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

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    std::ifstream ifs("../data/problems/problem-10.json");
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

    std::ofstream ofs("problem-10.json");
    ofs << solution.to_json().dump(4);

    return 0;
}