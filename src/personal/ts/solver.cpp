#ifdef USE_OPENCV
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-enum-enum-conversion"
#include <opencv2/core.hpp>
#include <opencv2/core/utils/logger.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#pragma GCC diagnostic pop
#endif

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
    double volume_before, volume_after;
    void apply(Solution& solution) const {
        switch (change_type) {
            case 0:
            std::swap(solution.placements[i], solution.placements[j]);
            break;
            case 1:
            solution.placements[i] = i_after;
            break;
            case 2:
            solution.volumes[i] = volume_after;
            break;
        }
    }
    static Changeset sample_random_volume(const Problem& problem, Xorshift& rnd, const Solution& solution) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;
        return Changeset { 2, i, -1, 
            {}, {}, {}, {},
            solution.volumes[i], solution.volumes[i] > 0.5 ? 0.2 : 1.0}; // flip
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
    static Changeset sample_random_delta(const Problem& problem, Xorshift& rnd, const Solution& solution, double max_delta) {
        const size_t N = solution.placements.size();
        const int i = rnd.next_int() % N;

        Changeset chg { 1, i, -1, 
            solution.placements[i], solution.placements[i], 
            {0.0, 0.0}, {0.0, 0.0},
        };

        for (int retry = 0; retry < 100; ++retry) {
            Placement placement = solution.placements[i];
            placement.x += (rnd.next_double() * 2.0 - 1.0) * max_delta;
            placement.y += (rnd.next_double() * 2.0 - 1.0) * max_delta;
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

        return sample_random_motion(problem, rnd, solution, i);
    }

    static Changeset sample_random_motion(const Problem& problem, Xorshift& rnd, const Solution& solution, int i) {
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
    if (input_solution.volumes.empty()) {
        input_solution.set_default_volumes();
    }
    int64_t input_score = compute_score(problem, input_solution);
    DUMP(input_score);

    Solution best_solution = input_solution;

    const int num_force_reset_iter = 10000;
    CachedComputeScore cache(problem);
    cache.full_compute(best_solution);
    int64_t best_score = input_score;
    LOG_ASSERT(cache.score() == best_score);


    int loop = 0;
    int accept = 0, reject = 0;
    if (method == "HC") {
        while (timer.elapsed_ms() < t_max) {
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(best_solution);
                adjust_score += best_score;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d, "adjust":%d})", loop, best_score, accept, reject, adjust_score);
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
    if (method == "HCAFFECTED") {
        cache.m_compute_affected = true;
        while (timer.elapsed_ms() < t_max) {
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(best_solution);
                adjust_score += best_score;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d, "adjust":%d})", loop, best_score, accept, reject, adjust_score);
            }
            Changeset changeset;
            if (rnd.next_int() % 100 <= 1) {
                changeset = Changeset::sample_random_mutation(problem, rnd, best_solution);
                cache.change_musician(changeset.i, changeset.i_after);
                cache.change_musician(changeset.j, changeset.j_after);
            } else {
                int i = 0;
                if (rnd.next_double() < 0.1) {
                    i = rnd.next_int() % problem.musicians.size();
                } else {
                    i = cache.sample_random_affected_musician(rnd);
                }
                changeset = Changeset::sample_random_motion(problem, rnd, best_solution, i);
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
    if (method == "NEARHC") {
        while (timer.elapsed_ms() < t_max) {
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(cache.m_solution);
                adjust_score += best_score;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d, "adjust":%d})", loop, best_score, accept, reject, adjust_score);
            }
            Changeset changeset;
            changeset = Changeset::sample_random_delta(problem, rnd, cache.m_solution, 100.0);
            cache.change_musician(changeset.i, changeset.i_after);
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                best_solution = cache.m_solution;
                ++accept;
            } else {
                if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                ++reject;
            }
            LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d})", loop, best_score, accept, reject);
        }
        DUMP(loop);
        DUMP(best_score);
    }
    if (method == "LINEHC") {
        while (timer.elapsed_ms() < t_max) {
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(cache.m_solution);
                adjust_score += best_score;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d, "adjust":%d})", loop, best_score, accept, reject, adjust_score);
            }
            // 適当な点を定め、そこまでの間を分割して試して最も良いものを選ぶ
            auto changeset_goal = Changeset::sample_random_motion(problem, rnd, cache.m_solution);
            const int num_split = 10;
            bool found = false;
            for (int i = 0; i < num_split; ++i) {
                double t = double(i) / double(num_split - 1);
                Changeset changeset = changeset_goal;
                changeset.i_after = {
                    changeset.i_before.x * t + changeset_goal.i_after.x * (1.0 - t),
                    changeset.i_before.y * t + changeset_goal.i_after.y * (1.0 - t),
                };
                bool conflict = false;
                for (int kk = 0; kk < cache.m_solution.placements.size(); ++kk) {
                    if (changeset.i != kk) {
                        if (are_musicians_too_close(cache.m_solution.placements[kk], changeset.i_after)) {
                            conflict = true;
                            break;
                        }
                    }
                }
                if (conflict) continue;
                cache.change_musician(changeset.i, changeset.i_after);
                int64_t score = cache.score();
                if (chmax(best_score, score)) {
                    best_solution = cache.m_solution;
                    found = true;
                } else {
                    if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                }
                break;
            }
            if (found) { ++accept; } else { ++reject; }

            if (loop % 100 == 0) {
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "accept":%d, "reject":%d})", loop, best_score, accept, reject);
                cv::Mat img = cache.m_solution.to_mat(problem);
                cv::Mat resized;
                cv::resize(img, resized, cv::Size(img.cols / 4, img.rows / 4)); // Resize the image
                cv::imshow("img", resized);
                cv::waitKey(1);
            }
        }
        DUMP(loop);
        DUMP(best_score);
    }
    if (method == "SA") {
        cache.m_compute_affected = false;
        double T_start = std::abs(best_score) * 0.01;
        double T_stop = std::abs(best_score) * 0.0001;
        double t = 0.0;
        Solution current_solution = best_solution;
        while ((t = timer.elapsed_ms()) < t_max) {
            //const double T = T_stop + (T_start - T_stop) * (1.0 - t / t_max);
            const double T = T_start * 1000 / (loop + 1);
            loop++;
            if (num_force_reset_iter > 0 &&  loop % num_force_reset_iter == 0) {
                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(best_solution);
                adjust_score += best_score;
                int64_t score = cache.full_compute(cache.m_solution);
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f, "adjust":%d})", loop, best_score, score, accept, reject, T, adjust_score);
            }
            Changeset changeset;
            int64_t gain = 0;
            double action = rnd.next_double();
            if (action < 0.01) {
                changeset = Changeset::sample_random_mutation(problem, rnd, cache.m_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after) + cache.change_musician(changeset.j, changeset.j_after);
            } else if (action < 0.0) {
                changeset = Changeset::sample_random_delta(problem, rnd, cache.m_solution, 5.0);
                gain = cache.change_musician(changeset.i, changeset.i_after);
            } else if (action < 0.00) {
                changeset = Changeset::sample_random_volume(problem, rnd, cache.m_solution);
                gain = cache.change_musician_volume(changeset.i, changeset.volume_after);
            } else if (action < 0.02) {
                changeset = Changeset::sample_random_motion(problem, rnd, cache.m_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after);
            } else {
                changeset = Changeset::sample_random_motion(problem, rnd, cache.m_solution, cache.sample_random_affected_musician(rnd));
                gain = cache.change_musician(changeset.i, changeset.i_after);
            }
            const double p = std::exp(double(gain) / T);
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                best_solution = cache.m_solution;
                ++accept;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f})", loop, best_score, score, accept, reject, T);
            } else {
                if (gain > 0 || rnd.next_double() < p) {
                    ++accept;
                } else {
                    ++reject;
                    if (changeset.change_type <= 1) {
                        if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                        if (changeset.j >= 0) cache.change_musician(changeset.j, changeset.j_before);
                    } 
                    if (changeset.change_type == 2) {
                        if (changeset.i >= 0) cache.change_musician_volume(changeset.i, changeset.volume_before);
                    } 
                }
                if (loop % 1000 == 0) {
                    LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f})", loop, best_score, score, accept, reject, T);
                }
            }
        }
        std::cout << "Volumes:";
        for (int k = 0; k < best_solution.volumes.size(); ++k) {
            std::cout << best_solution.volumes[k] << " ";
            if (best_solution.volumes[k] < 0.5) {
                best_solution.volumes[k] = 0.0;
            } else {
                best_solution.volumes[k] = 1.0; // to be compatible.
            }
        }
        std::cout << std::endl;
        DUMP(loop, best_score, accept, reject);
    }
    if (method == "SAVOL") { // SA <-> VOL
        double T_start = 1e5;
        double T_stop = 1e3;
        double t = 0.0;
        Solution current_solution = best_solution;
        while ((t = timer.elapsed_ms()) < t_max) {
            const double T = T_start * 1000 / (loop + 1);
            loop++;
            Changeset changeset;
            int64_t gain = 0;
            double action = rnd.next_double();
            if (action < 0.0) {
                changeset = Changeset::sample_random_mutation(problem, rnd, cache.m_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after) + cache.change_musician(changeset.j, changeset.j_after);
            } else if (action < 0.0) {
                changeset = Changeset::sample_random_delta(problem, rnd, cache.m_solution, 5.0);
                gain = cache.change_musician(changeset.i, changeset.i_after);
            } else {
                changeset = Changeset::sample_random_motion(problem, rnd, cache.m_solution);
                gain = cache.change_musician(changeset.i, changeset.i_after);
            }
            const double p = std::exp(double(gain) / T);
            int64_t score = cache.score();
            if (chmax(best_score, score)) {
                best_solution = cache.m_solution;
                ++accept;
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f})", loop, best_score, score, accept, reject, T);
            } else {
                if (gain > 0 || rnd.next_double() < p) {
                    ++accept;
                } else {
                    ++reject;
                    if (changeset.change_type <= 1) {
                        if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                        if (changeset.j >= 0) cache.change_musician(changeset.j, changeset.j_before);
                    } 
                    if (changeset.change_type == 2) {
                        if (changeset.i >= 0) cache.change_musician_volume(changeset.i, changeset.volume_before);
                    } 
                }
                if (loop % 1000 == 0) {
                    LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f})", loop, best_score, score, accept, reject, T);
                }
            }

            if (loop % 10000 == 0) {
                Solution dummy_solution = best_solution;
                set_optimal_volumes(problem, dummy_solution);
                for (int i = 0; i < cache.m_solution.volumes.size(); ++i) {
                    if (dummy_solution.volumes[i] < 5.0) {
                        cache.m_solution.volumes[i] = 0.001;
                    } else {
                        cache.m_solution.volumes[i] = 1.0;
                    }
                }

                int64_t adjust_score = -best_score;
                best_score = cache.full_compute(best_solution);
                adjust_score += best_score;
                int64_t score = cache.full_compute(cache.m_solution);
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d, "T":%f, "adjust":%d})", loop, best_score, score, accept, reject, T, adjust_score);
            }
        }
        std::cout << "Volumes:";
        for (int k = 0; k < best_solution.volumes.size(); ++k) {
            std::cout << best_solution.volumes[k] << " ";
            if (best_solution.volumes[k] < 0.5) {
                best_solution.volumes[k] = 0.0;
            } else {
                best_solution.volumes[k] = 1.0; // to be compatible.
            }
        }
        std::cout << std::endl;
        DUMP(loop, best_score, accept, reject);
    }
    if (method == "ILS") {
        double t = 0.0;
        Solution current_solution = best_solution;
        while ((t = timer.elapsed_ms()) < t_max) {
            // LS
            bool reset = false;
            int64_t score = 0;
            for (int ls = 0; ls < 1000 && timer.elapsed_ms() < t_max; ++ls) {
                loop++;
                auto changeset = Changeset::sample_random_motion(problem, rnd, cache.m_solution);
                auto gain = cache.change_musician(changeset.i, changeset.i_after);
                score = cache.score();
                if (chmax(best_score, score)) {
                    best_solution = cache.m_solution;
                    ++accept;
                } else {
                    ++reject;
                    if (changeset.i >= 0) cache.change_musician(changeset.i, changeset.i_before);
                }
                if (ls % 100 == 0) {
                    LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d})", loop, best_score, score, accept, reject);
                }
            }
            { // kick from the best (and it also resets the cumulative error)
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d})", loop, best_score, score, accept, reject);
                best_score = cache.full_compute(best_solution);
                //auto changeset = Changeset::sample_random_mutation(problem, rnd, cache.m_solution);
                //cache.change_musician(changeset.i, changeset.i_after);
                //cache.change_musician(changeset.j, changeset.j_after);
                for (int i = 0; i < 3; ++i) { // as a kick, repeate for a few times.
                    auto changeset = Changeset::sample_random_motion(problem, rnd, cache.m_solution);
                    cache.change_musician(changeset.i, changeset.i_after);
                }
                ++accept;
                DUMP("kick", loop, best_score, accept, reject);
                LOG(INFO) << format(R"(RECORD {"loop": %d, "best":%lld, "current":%lld, "accept":%d, "reject":%d})", loop, best_score, score, accept, reject);
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
        std::ofstream ofs(format("../data/solutions/ts_v01_sa/solution-%d.json", problem_id));
        ofs << best_solution.to_json().dump(4);
    }

    return 0;
}