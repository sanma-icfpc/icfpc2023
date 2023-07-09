#include "stdafx.h"

#include <glog/logging.h>
#include <gtest/gtest.h>

#include "../util.h"

TEST(TestUtil, GuessProblemId) {
  EXPECT_EQ(guess_problem_id("path\\to\\problem-15.json"), 15);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-16.json"), 16);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-17_12345.json"), 17);
}

TEST(TestUtil, IsIntersect) {
  // 1 point inside, 1 point outside.
  EXPECT_TRUE(is_intersect(0, 0, 10, 0, 0, 10, 10));
  EXPECT_TRUE(is_intersect(0, 0, 10, -1, -1, 10, 10));
  EXPECT_TRUE(is_intersect(0, 0, 10, 1, 1, 10, 10));
  EXPECT_TRUE(is_intersect(0, 0, 10, 2, 2, 10, 10));
  EXPECT_TRUE(is_intersect(0, 0, 10, 5, 5, 10, 10));

  // both outside, crossing
  EXPECT_TRUE(is_intersect(0, 0, 10, -10, -10, 10, 10));
  EXPECT_TRUE(is_intersect(0, 0, 10, -10, -9, 10, -9));
  EXPECT_TRUE(is_intersect(0, 0, 10, -10, -10, 10, -10));  // Touching
  EXPECT_FALSE(is_intersect(0, 0, 10, -10, -11, 10, -11));

  // both inside
  EXPECT_FALSE(is_intersect(0, 0, 10, 0, 0, 5, 5));
}

TEST(TestUtil, CachedComputeScore_InverseOperationCancelsOutScoreGain) {
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

TEST(TestUtil, GuessExtension) {
  EXPECT_EQ(Problem::from_file(42).extension, Extension::lightning());
  EXPECT_EQ(Problem::from_file(56).extension, Extension::full());
}

TEST(TestUtil, RegressionTest_ComputeScoreFunctionsWithoutPillars) {
  Problem problem = Problem::from_file(42);
  Solution solution = Solution::from_file(
      "../data/solutions/k3_v01_random_creation/solution-42.json");

  auto score_naive = compute_score(problem, solution);
  auto score_fast = compute_score_fast(problem, solution);
  auto score_cached = CachedComputeScore(problem).full_compute(solution);
  LOG(INFO) << "score_naive: " << score_naive;
  LOG(INFO) << "score_fast : " << score_fast;
  LOG(INFO) << "score_cached : " << score_cached;
  EXPECT_EQ(score_naive, 4382334);   // revision:87ac52b
  EXPECT_EQ(score_fast, 4382334);    // revision:87ac52b
  EXPECT_EQ(score_cached, 4382334);  // revision:87ac52b
}

// Expected values are discussed in the official Discord channel #question.
TEST(TestUtil, DISABLED_ExtensionsInComputeScoreFunction) {
  Problem problem =
      Problem::from_file("../data/test_data/problem-example.json");
  Solution solution =
      Solution::from_file("../data/test_data/solution-example.json");

  const Extension extension_pillars{.consider_pillars = true,
                                    .consider_harmony = false};
  const Extension extension_harmony{.consider_pillars = false,
                                    .consider_harmony = true};

  auto score_naive_full = compute_score(problem, solution);
  problem.extension = extension_pillars;
  auto score_naive_pillars = compute_score(problem, solution);

  // `.consider_pillars` is actually ignored for performance reasons.
  // Here we clear out problem.pillars instead for linghtning and
  // harmony tests.
  problem.pillars = std::vector<Pillar>();

  problem.extension = Extension::lightning();
  auto score_naive_lightning = compute_score(problem, solution);
  problem.extension = extension_harmony;
  auto score_naive_harmony = compute_score(problem, solution);

  EXPECT_DOUBLE_EQ(score_naive_lightning, 5343.0);
  EXPECT_DOUBLE_EQ(score_naive_full, 3270.0);
  // EXPECT_DOUBLE_EQ(score_naive_pillars, ?);
  EXPECT_DOUBLE_EQ(score_naive_harmony, 5350.0);
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
