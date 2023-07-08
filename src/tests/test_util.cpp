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

TEST(TestUtil, regression_test_for_compute_score_functions_without_pillars) {
  Problem problem = Problem::from_file(42);
  Solution solution = Solution::from_file("../data/solutions/k3_v01_random_creation/solution-42.json");

  Extension extension { .consider_pillars = false, .consider_harmony = false };
  auto score_naive = compute_score(problem, solution, extension);
  auto score_fast = compute_score_fast(problem, solution);
  auto score_cached = CachedComputeScore(problem).full_compute(solution);
  LOG(INFO) << "score_naive: " << score_naive;
  LOG(INFO) << "score_fast : " << score_fast;
  LOG(INFO) << "score_cached : " << score_cached;
  EXPECT_EQ(score_naive, 4382334); // revision:87ac52b
  EXPECT_EQ(score_fast, 4382334); // revision:87ac52b
  EXPECT_EQ(score_cached, 4382334); // revision:87ac52b
}

TEST(TestUtil, minimal_test_to_check_effect_of_extensions_in_compute_score_functions) {
  Xorshift rnd(42);
  Problem problem = Problem::from_file(56);
  Solution solution = *create_random_solution(problem, rnd);

  Extension extension_lightning { .consider_pillars = false, .consider_harmony = false };
  Extension extension_pillars { .consider_pillars = true, .consider_harmony = false };
  Extension extension_harmony { .consider_pillars = false, .consider_harmony = true };
  auto score_naive_lightning = compute_score(problem, solution, extension_lightning);
  auto score_naive_pillars = compute_score(problem, solution, extension_pillars);
  auto score_naive_harmony = compute_score(problem, solution, extension_harmony);
  LOG(INFO) << "score_naive(lightning): " << score_naive_lightning;
  LOG(INFO) << "score_naive(pillars): " << score_naive_pillars;
  LOG(INFO) << "score_naive(harmony): " << score_naive_harmony;
  EXPECT_LT(score_naive_pillars, score_naive_lightning);
  EXPECT_GT(score_naive_harmony, score_naive_lightning);
}

TEST(TestUtil, test_compute_score_fast) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);

  Extension extension { .consider_pillars = false, .consider_harmony = false };
  for (int trial = 0; trial < 10; trial++) {
    Solution solution = *create_random_solution(problem, rnd);
    auto score_naive = compute_score(problem, solution, extension);
    auto score_fast = compute_score_fast(problem, solution);
    LOG(INFO) << "score_naive: " << score_naive;
    LOG(INFO) << "score_fast : " << score_fast;
    EXPECT_NEAR(score_naive, score_fast, 100);
  }
}

// vim:ts=2 sw=2 sts=2 et ci
