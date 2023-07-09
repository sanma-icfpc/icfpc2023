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
    static Changeset sample_random_mutation(const Problem& problem, Xorshift& rnd, const Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;
        const int j = (i + 1 + rnd.next_int() % (N - 1)) % N;
        LOG_ASSERT(0 <= i && i < N && 0 <= j && j < N && i != j);

        return Changeset { 0, i, j, 
            solution.placements[i], solution.placements[j],
            solution.placements[j], solution.placements[i],
        };
    }
    static Changeset sample_random_delta(const Problem& problem, Xorshift& rnd, const Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;

        Changeset chg { 1, i, -1, 
            solution.placements[i], solution.placements[i], 
            {0.0, 0.0}, {0.0, 0.0},
        };

        for (int retry = 0; retry < 100; ++retry) {
            Placement placement = solution.placements[i];
            placement.x += rnd.next_double() * 3.0 - 1.5;
            placement.y += rnd.next_double() * 3.0 - 1.5;
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
            chg.i_after = placement;
        }

        return chg;
    }
    static Changeset sample_random_motion(const Problem& problem, Xorshift& rnd, const Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;

        Changeset chg { 1, i, -1, 
            solution.placements[i], solution.placements[i], 
            {0.0, 0.0}, {0.0, 0.0},
        };

        for (int retry = 0; retry < 100; ++retry) {
            Placement placement {
                problem.stage_x + k_musician_spacing_radius + rnd.next_double() * (problem.stage_w - k_musician_spacing_radius * 2),
                problem.stage_y + k_musician_spacing_radius + rnd.next_double() * (problem.stage_h - k_musician_spacing_radius * 2)};
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
            if (!conflict) {
                chg.i_after = placement;
                break;
            }
        }

        return chg;
    }
};

int main([[maybe_unused]] int argc, [[maybe_unused]] char** argv) {

  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  google::SetStderrLogging(google::INFO);
  google::SetLogDestination(google::INFO, "tssolver.log.");

  std::ios::sync_with_stdio(false);
  std::cin.tie(NULL);

    std::string solution_or_problem_path;
    std::string method = "SA";
    double t_max = 10000.0;
    if (argc >= 2) {
        solution_or_problem_path = argv[1];
    }
    if (argc >= 3) {
        method = argv[2];
    }
    if (argc >= 4) {
        t_max = std::atof(argv[3]);
    }
    int problem_id = *guess_problem_id(solution_or_problem_path);
    DUMP(problem_id);

    Timer timer;

#ifdef HAVE_OPENCV_HIGHGUI
    cv::utils::logging::setLogLevel(cv::utils::logging::LogLevel::LOG_LEVEL_SILENT);
#endif

    Xorshift rnd(42);
    Problem problem = Problem::from_file(problem_id);
    Solution input_solution = Solution::from_file(solution_or_problem_path);
    if (input_solution.placements.empty()) {
        LOG(WARNING) << "provided a problem file. starting from a random state.";
        input_solution = *create_random_solution(problem, rnd);
    }
    int64_t input_score = compute_score(problem, input_solution);
    DUMP(input_score);

    Solution best_solution = input_solution;

    const int num_force_reset_iter = 3000;
    CachedComputeScore cache(problem);
    cache.full_compute(best_solution);
    int64_t best_score = input_score;
    LOG_ASSERT(cache.score() == best_score);


    if (method == "HILLCLIMB") {
        int loop = 0;
        int accept = 0, reject = 0;
        while (timer.elapsed_ms() < t_max) {
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                best_score = cache.full_compute(best_solution);
            }
            Changeset changeset;
            if (rnd.next_int() % 100 <= 1) {
                changeset = Changeset::sample_random_mutation(problem, rnd, best_solution);
                cache.change_musician(changeset.i, changeset.i_after);
                cache.change_musician(changeset.j, changeset.j_after);
            } else {
                changeset = Changeset::sample_random_motion(problem, rnd, best_solution);
                cache.change_musician(changeset.i, changeset.i_after);
            }
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                changeset.apply(best_solution);
                ++accept;
            } else {
                if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                if (changeset.j >= 0) cache.change_musician(changeset.j, changeset.j_before);
                ++reject;
            }
            LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d})", loop, best_score, accept, reject);
        }
        DUMP(loop);
        DUMP(best_score);
    }
    if (method == "SA") {
        int loop = 0;
        double T_start = 1e5;
        double T_stop = 1e3;
        double t = 0.0;
        Solution current_solution = best_solution;
        int accept = 0, reject = 0;
        while ((t = timer.elapsed_ms()) < t_max) {
            const double T = T_stop + (T_start - T_stop) * (1.0 - t / t_max);
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) cache.full_compute(current_solution);
            Changeset changeset;
            int64_t gain = 0;
            double action = rnd.next_double();
            if (action < 0.0) {
                changeset = Changeset::sample_random_mutation(problem, rnd, current_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after) + cache.change_musician(changeset.j, changeset.j_after);
            } else if (action < 0.0) {
                changeset = Changeset::sample_random_delta(problem, rnd, current_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after);
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
            } else {
                // DUMP(t / t_max, T, gain, p);
                if (gain > 0 || rnd.next_double() < p) {
                    ++accept;
                    changeset.apply(current_solution);
                } else {
                    ++reject;
                    if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                    if (changeset.j >= 0) cache.change_musician(changeset.j, changeset.j_before);
                }
            }
        }
        DUMP(loop, best_score, accept, reject);
    }
    if (method == "ILS") {
        int loop = 0;
        double t = 0.0;
        Solution current_solution = best_solution;
        int accept = 0, reject = 0;
        while ((t = timer.elapsed_ms()) < t_max) {
            // LS
            bool reset = false;
            for (int ls = 0; ls < 100 && timer.elapsed_ms() < t_max; ++ls) {
                loop++;
                if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) reset = true;
                auto changeset = Changeset::sample_random_motion(problem, rnd, current_solution);
                auto gain = cache.change_musician(changeset.i, changeset.i_after);
                int64_t score = cache.score();
                if (chmax(best_score, score)) {
                    changeset.apply(current_solution);
                    best_solution = current_solution;
                    ++accept;
                    //DUMP("LS", loop, best_score, accept, reject);
                } else {
                    ++reject;
                    if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                }
            }
            { // kick
                if (reset) cache.full_compute(current_solution);
                auto changeset = Changeset::sample_random_mutation(problem, rnd, current_solution);
                cache.change_musician(changeset.i, changeset.i_after);
                cache.change_musician(changeset.j, changeset.j_after);
                changeset.apply(current_solution);
                ++accept;
                DUMP("kick", loop, best_score, accept, reject);
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

    cache.report();

    if (best_score > 0) {
        std::ofstream ofs(format("../data/solutions/ts/solution-%d.json", problem_id));
        ofs << best_solution.to_json().dump(4);
    }

    return 0;
}