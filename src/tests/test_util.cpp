#include "stdafx.h"

#include <glog/logging.h>
#include <gtest/gtest.h>

#include "../util.h"

TEST(TestUtil, guess_problem_id) {
  EXPECT_EQ(guess_problem_id("path\\to\\problem-15.json"), 15);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-16.json"), 16);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-17_12345.json"), 17);
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
