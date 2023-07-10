#include "stdafx.h"

#include <glog/logging.h>
#include <gtest/gtest.h>
#include <ranges>

#include "../cached_compute_score.h"
#include "../util.h"

TEST(CachedComputeScoreTest, IdenticalChangeWillNotChangeScore) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  int64_t score_before = cache.score();
  int64_t gain_forward = cache.change_musician(1, solution.placements[1]);
  int64_t score_after = cache.score();
  EXPECT_EQ(gain_forward, 0);
  EXPECT_EQ(score_before, score_after);
}

TEST(CachedComputeScoreTest, InverseOperationCancelsOutScoreGain) {
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

TEST(CachedComputeScoreTest, DryRunComputesExactScoreDiff) {
  Problem problem = Problem::from_file(85);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);

  CachedComputeScore cache(problem);
  cache.full_compute(solution);

  int64_t initial_score = cache.score();
  int64_t dry_run_diff = cache.change_musician(1, solution.placements[2], true);
  int64_t simple_run_diff = cache.change_musician(1, solution.placements[2], false);
  EXPECT_EQ(dry_run_diff, simple_run_diff);
}

TEST(CachedComputeScoreTest, DryRunDoesNotChangeInternalState) {
  Problem problem = Problem::from_file(85);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);

  CachedComputeScore cache(problem);
  cache.full_compute(solution);

  int64_t initial_score = cache.score();
  EXPECT_NE(cache.change_musician(1, solution.placements[2], true), 0); // 99% sure
  EXPECT_NE(cache.change_musician(2, solution.placements[3], true), 0); // 99% sure
  EXPECT_EQ(cache.score(), initial_score);
  EXPECT_EQ(cache.m_solution.placements, solution.placements);
}

TEST(CachedComputeScoreTest, RunAfterDryRunEqualsSimpleRun) {
  Problem problem = Problem::from_file(85);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);

  CachedComputeScore cache(problem);

  cache.full_compute(solution);
  EXPECT_NE(cache.change_musician(1, solution.placements[2], true), 0); // trial->revert
  EXPECT_NE(cache.change_musician(2, solution.placements[3], true), 0); // trial->revert
  EXPECT_NE(cache.change_musician(3, solution.placements[4], false), 0); // trial->accept
  int64_t score_run_after_dry_run = cache.score();

  cache.full_compute(solution);
  EXPECT_NE(cache.change_musician(3, solution.placements[4], false), 0); // trial->accept
  int64_t score_simple_run = cache.score();

  EXPECT_EQ(score_run_after_dry_run, score_simple_run);
}

TEST(CachedComputeScoreTest, ConsistentWithReferenceScoreDuringRandomChange) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);

  for (int loop = 0; loop < 100; ++loop) {
    const int i = rnd.next_int() % M;
    auto new_placement_opt =
        suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt));                        // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0);  // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_EQ(cache.score(), score_reference);
  }

  LOG(INFO) << format(
      "full calculation %.4f ms/eval, partial update %.4f ms/eval",
      cache.get_mean_elapsed_ms_full(), cache.get_mean_elapsed_ms_partial());
}

TEST(CachedComputeScoreTest,
     ConsistentWithReferenceScoreDuringRandomChangeWithPillars) {
  Problem problem =
      Problem::from_file("../data/test_data/problem-85-mod-no-harmony.json");
  problem.extension.consider_pillars = true;
  problem.extension.consider_harmony = false;
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();
  LOG(INFO) << "Pillars: " << problem.pillars.size();
  LOG(INFO) << problem.extension.stringify();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  EXPECT_EQ(cache.score(), compute_score(problem, solution));

  for (int loop = 0; loop < 100; ++loop) {
    const int i = rnd.next_int() % M;
    auto new_placement_opt =
        suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt));                        // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0);  // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference,
                std::abs(score_reference) * 1e-4);
  }
}

TEST(CachedComputeScoreTest,
     ConsistentWithReferenceScoreDuringRandomChangeWithHarmony) {
  Problem problem = Problem::from_file(85);
  problem.extension.consider_pillars = false;
  problem.pillars.clear();
  LOG(INFO) << problem.extension.stringify();

  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  EXPECT_EQ(cache.score(), compute_score(problem, solution));

  for (int loop = 0; loop < 100; ++loop) {
    const int i = rnd.next_int() % M;
    auto new_placement_opt =
        suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt));                        // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0);  // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference,
                std::abs(score_reference) * 1e-2);
    {
      CachedComputeScore cache_ref(problem);
      cache_ref.full_compute(cache.m_solution);
      for (int i = 0; i < cache.m_harmony_cache.size(); ++i) {
        EXPECT_NEAR(cache.m_harmony_cache[i], cache_ref.m_harmony_cache[i],
                    1e-3);
      }
    }
  }
}

TEST(CachedComputeScoreTest,
     ConsistentWithReferenceScoreDuringRandomChangeFullDivision) {
  Problem problem = Problem::from_file("../data/problems/problem-85.json");
  LOG(INFO) << problem.extension.stringify();
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  EXPECT_EQ(cache.score(), compute_score(problem, solution));

  for (int loop = 0; loop < 100; ++loop) {
    const int i = rnd.next_int() % M;
    auto new_placement_opt =
        suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt));                        // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0);  // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference,
                std::abs(score_reference) * 1e-2);
  }

  LOG(INFO) << format(
      "full calculation %.4f ms/eval, partial update %.4f ms/eval",
      cache.get_mean_elapsed_ms_full(), cache.get_mean_elapsed_ms_partial());
}

TEST(CachedComputeScoreTest, ConsistentWithReferenceScoreOnVolumeChange) {
  Problem problem = Problem::from_file(42);
  Solution solution = Solution::from_file(
      "../data/solutions/k3_v01_random_creation/solution-42.json");

  solution.set_default_volumes();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  EXPECT_EQ(cache.score(), compute_score(problem, solution));

  cache.change_musician_volume(0, 5.0);
  cache.change_musician_volume(1, 8.0);
  EXPECT_EQ(cache.m_solution.volumes[0], 5.0);
  EXPECT_EQ(cache.m_solution.volumes[1], 8.0);
  EXPECT_EQ(cache.score(), compute_score(problem, cache.m_solution));
}

TEST(CachedComputeScoreTest, HillClimbingForLongTimeForProfiling) {
  Problem problem = Problem::from_file(79);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);
  const auto init_score = cache.score();

  Timer timer;
  while (timer.elapsed_ms() < 20.0 * 1000.0) {
    const int i = rnd.next_int() % M;
    auto old_placement = cache.m_solution.placements[i];
    auto new_placement_opt =
        suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt));  // 99% sure..
    if (cache.change_musician(i, *new_placement_opt) < 0) {
      cache.change_musician(i, old_placement);
    }
  }

  LOG(INFO) << format("init  score = %s",
                      int_to_delimited_string(init_score).c_str());
  LOG(INFO) << format("last  score = %s",
                      int_to_delimited_string(cache.score()).c_str());
  LOG(INFO) << format(
      "final score = %s",
      int_to_delimited_string(
          CachedComputeScore(problem).full_compute(cache.m_solution))
          .c_str());

  cache.report();
}

TEST(CachedComputeScoreTest,
     RegressionTest_ComputeScoreFunctionsWithoutPillars) {
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

TEST(CachedComputeScoreTest, ComputeScoreFunctionWith10xVolumeScores10x) {
  Problem problem = Problem::from_file(42);
  Solution solution = Solution::from_file(
      "../data/solutions/k3_v01_random_creation/solution-42.json");

  auto score_naive = compute_score(problem, solution);
  auto score_cached = CachedComputeScore(problem).full_compute(solution);
  EXPECT_EQ(score_naive, score_cached);

  solution.set_default_volumes();
  solution.volumes.assign(solution.volumes.size(), 10.0);
  auto score_naive_loud = compute_score(problem, solution);
  auto score_cached_loud = CachedComputeScore(problem).full_compute(solution);
  EXPECT_EQ(score_naive * 10.0, score_naive_loud);
  EXPECT_EQ(score_cached * 10.0, score_cached_loud);
}
