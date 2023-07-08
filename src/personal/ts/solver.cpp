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
    int change_type = 0;
    int i = -1;
    int j = -1;
    Placement i_before, i_after;
    Placement j_before, j_after;
    void apply(Solution& solution) const {
        switch (change_type) {
            case 0:
            std::swap(solution.placements[i], solution.placements[j]);
            break;
            case 1:
            solution.placements[i] = i_after;
            break;
        }
    }
    //void unapply(Solution& solution) const {
    //    switch (change_type) {
    //        case 0:
    //        std::swap(solution.placements[i], solution.placements[j]);
    //        break;
    //        case 1:
    //        break;
    //    }
    //}
    static Changeset sample_random_mutation(const Problem& problem, Xorshift& rnd, Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;
        const int j = (i + 1 + rnd.next_int() % (N - 1)) % N;
        LOG_ASSERT(0 <= i && i < N && 0 <= j && j < N && i != j);

        return Changeset { 0, i, j, 
            solution.placements[i], solution.placements[j],
            solution.placements[j], solution.placements[i],
        };
    }
    static Changeset sample_random_motion(const Problem& problem, Xorshift& rnd, Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;

        Changeset chg { 1, i, -1, 
            solution.placements[i], solution.placements[i], 
            {0.0, 0.0}, {0.0, 0.0},
        };

        for (int retry = 0; retry < 100; ++retry) {
            Placement placement {
                problem.stage_x + rnd.next_double() * problem.stage_w,
                problem.stage_y + rnd.next_double() * problem.stage_h};
            bool conflict = false;
            for (int kk = 0; kk < solution.placements.size(); ++kk) {
                if (i != kk) {
                    if (are_musicians_too_close(solution.placements[kk], placement)) {
                        conflict = true;
                        break;
                    }
                }
            }
            if (!conflict) {
                chg.i_after = placement;
                break;
            }
        }

        return chg;
    }
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

    std::string solution_path;
    std::string method = "SA";
    if (argc >= 2) {
        solution_path = argv[1];
    }
    if (argc >= 3) {
        method = argv[2];
    }
    int problem_id = *guess_problem_id(solution_path);
    DUMP(problem_id);

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    Problem problem = Problem::from_file(problem_id);
    Solution input_solution = Solution::from_file(solution_path);
    int64_t input_score = compute_score(problem, input_solution);
    DUMP(input_score);

    Xorshift rnd(42);
    Solution best_solution = input_solution;

    CachedComputeScore cache(problem);
    cache.full_compute(best_solution);
    int64_t best_score = input_score;

    if (method == "HILLCLIMB") {
        int loop = 0;
        while (timer.elapsed_ms() < 10000) {
            loop++;
            auto changeset = Changeset::sample_random_mutation(problem, rnd, best_solution);
            cache.change_musician(changeset.i, changeset.i_after);
            cache.change_musician(changeset.j, changeset.j_after);
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                changeset.apply(best_solution);
                DUMP(loop, best_score);
            } else {
                cache.change_musician(changeset.i, changeset.i_before);
                cache.change_musician(changeset.j, changeset.j_before);
            }
        }
        DUMP(loop);
        DUMP(best_score);
    }
    if (method == "SA") {
        int loop = 0;
        double T_start = 1e5;
        double T_stop = 1e3;
        double t = 0.0;
        double t_max = 10000.0;
        Solution current_solution = best_solution;
        int accept = 0, reject = 0;
        while ((t = timer.elapsed_ms()) < t_max) {
            const double T = T_stop + (T_start - T_stop) * (1.0 - t / t_max);
            loop++;
            Changeset changeset;
            int64_t gain = 0;
            if (rnd.next_double() < 0) {
                changeset = Changeset::sample_random_mutation(problem, rnd, current_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after) + cache.change_musician(changeset.j, changeset.j_after);
            } else {
                changeset = Changeset::sample_random_motion(problem, rnd, current_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after);
            }
            const double p = std::exp(double(gain) / T);
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                changeset.apply(current_solution);
                best_solution = current_solution;
                ++accept;
                DUMP(loop, best_score, accept, reject);
                DUMP(t / t_max, T, gain, p);
            } else {
                if (gain > 0 || rnd.next_double() < p) {
                    // accept
                    ++accept;
                    changeset.apply(current_solution);
                } else {
                    // reject
                    ++reject;
                    if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                    if (changeset.j >= 0) cache.change_musician(changeset.j, changeset.j_before);
                }
            }
        }
        DUMP(loop, best_score, accept, reject);
    }

    // verify
    int64_t final_score = compute_score(problem, best_solution);
    DUMP(final_score == best_score); // tends to be false due to some bug..
    std::cout << format("init_score  = %I64d (%s)", input_score, int_to_delimited_string(input_score).c_str()) << std::endl;
    std::cout << format("best_score  = %I64d (%s)", best_score, int_to_delimited_string(best_score).c_str()) << std::endl;
    std::cout << format("final_score = %I64d (%s)", final_score, int_to_delimited_string(final_score).c_str()) << std::endl;

    if (best_score > 0) {
        std::ofstream ofs(format("../data/solutions/ts/solution-%d.json", problem_id));
        ofs << best_solution.to_json().dump(4);
    }

    return 0;
}