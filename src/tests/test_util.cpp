#include "stdafx.h"

#include <glog/logging.h>
#include <gtest/gtest.h>

#include "../util.h"

TEST(TestUtil, guess_problem_id) {
  EXPECT_EQ(guess_problem_id("path\\to\\problem-15.json"), 15);
}

TEST(TestUtil, CachedComputeScore_inverse_operation_cancels_out_score_gain) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  int64_t gain_forward = cache.change_musician(1, solution.placements[2]);
  LOG(INFO) << gain_forward;
  EXPECT_NE(gain_forward, 0);
  int64_t gain_backward = cache.change_musician(1, solution.placements[1]);
  LOG(INFO) << gain_backward;

  EXPECT_EQ(gain_forward + gain_backward, 0);
}

TEST(TestUtil, CachedComputeScore_change_musician) {

  std::string data_txt = R"({
    "room_width": 2000.0,
    "room_height": 5000.0,
    "stage_width": 1000.0,
    "stage_height": 200.0,
    "stage_bottom_left": [500.0, 0.0],
    "musicians": [0, 1, 0],
    "attendees": [
        {"x": 100.0, "y": 500.0, "tastes": [1000.0, -1000.0]},
        {"x": 200.0, "y": 1000.0, "tastes": [200.0, 200.0]},
        {"x": 1100.0, "y": 800.0, "tastes": [800.0, 1500.0]}
    ]
})";
  auto data = nlohmann::json::parse(data_txt);
  LOG(INFO) << data["room_width"];

  Problem problem(data);

  Xorshift rnd;
  auto solution = *create_random_solution(problem, rnd);
  LOG(INFO) << solution.placements.size();
  LOG(INFO) << problem.musicians.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  EXPECT_EQ(cache.score(), compute_score(problem, solution));

  for (int trial = 0; trial < 1121; trial++) {
    auto changeset = Changeset::sample_random_motion(problem, rnd, solution);
    auto gain = cache.change_musician(changeset.i, changeset.i_after);
    changeset.apply(solution);
    if (cache.score() != compute_score(problem, solution)) {
      DUMP(trial, gain, cache.score(), compute_score(problem, solution));
      return;
    }
    EXPECT_EQ(cache.score(), compute_score(problem, solution));
  }

  {
    auto changeset = Changeset::sample_random_motion(problem, rnd, solution);
    DUMP(changeset.change_type, changeset.i, changeset.j);
    DUMP(changeset.i_before.stringify(), changeset.i_after.stringify());
    DUMP(solution.stringify());
    auto gain = cache.change_musician(changeset.i, changeset.i_after);
    changeset.apply(solution);
    DUMP(solution.stringify());
    if (cache.score() != compute_score(problem, solution)) {
      DUMP(gain, cache.score(), compute_score(problem, solution));
      auto cache_clone(cache);
      cache_clone.full_compute(solution);
      DUMP(cache_clone.score());
      cache.check_equality(cache_clone);
    }
    EXPECT_EQ(cache.score(), compute_score(problem, solution));
  }

}

TEST(TestUtil, test_compute_score_fast) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);

  for (int trial = 0; trial < 10; trial++) {
    Solution solution = *create_random_solution(problem, rnd);
    auto score_naive = compute_score(problem, solution);
    auto score_fast = compute_score_fast(problem, solution);
    LOG(INFO) << "score_naive: " << score_naive;
    LOG(INFO) << "score_fast : " << score_fast;
    EXPECT_NEAR(score_naive, score_fast, 100);
  }
}

// vim:ts=2 sw=2 sts=2 et ci
