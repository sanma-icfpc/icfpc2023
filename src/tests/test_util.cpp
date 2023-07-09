#include "stdafx.h"

#include <ranges>
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

TEST(TestUtil, CachedComputeScore_IdenticalChangeWillNotChangeScore) {
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

TEST(TestUtil, CachedComputeScore_ConsistentWithReferenceScoreDuringRandomChange) {
  Problem problem = Problem::from_file(42);
  Xorshift rnd(42);
  Solution solution = *create_random_solution(problem, rnd);
  const int M = solution.placements.size();

  CachedComputeScore cache(problem);
  cache.full_compute(solution);

  for (int loop = 0; loop < 100; ++loop) {
    const int i = rnd.next_int() % M;
    auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt)); // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0); // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_EQ(cache.score(), score_reference);
  }

  LOG(INFO) << format("full calculation %.4f ms/eval, partial update %.4f ms/eval",
    cache.get_mean_elapsed_ms_full(),
    cache.get_mean_elapsed_ms_partial());
}

TEST(TestUtil, CachedComputeScore_ConsistentWithReferenceScoreDuringRandomChangeWithPillars) {
  Problem problem = Problem::from_file("../data/test_data/problem-85-mod-no-harmony.json");
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
    auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt)); // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0); // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference, std::abs(score_reference) * 1e-4);
  }
}

TEST(TestUtil, CachedComputeScore_ConsistentWithReferenceScoreDuringRandomChangeWithHarmony) {
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
    auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt)); // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0); // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference, std::abs(score_reference) * 1e-2);
    {
      CachedComputeScore cache_ref(problem);
      cache_ref.full_compute(cache.m_solution);
      for (int i = 0; i < cache.m_harmony_cache.size(); ++i) {
        EXPECT_NEAR(cache.m_harmony_cache[i], cache_ref.m_harmony_cache[i], 1e-3);
      }
    }
  }
}

TEST(TestUtil, CachedComputeScore_ConsistentWithReferenceScoreDuringRandomChangeFullDivision) {
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
    auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt)); // 99% sure..
    EXPECT_NE(cache.change_musician(i, *new_placement_opt), 0); // 99% sure..
    const int64_t score_reference = compute_score(problem, cache.m_solution);
    EXPECT_NEAR(cache.score(), score_reference, std::abs(score_reference) * 1e-2);
  }

  LOG(INFO) << format("full calculation %.4f ms/eval, partial update %.4f ms/eval",
    cache.get_mean_elapsed_ms_full(),
    cache.get_mean_elapsed_ms_partial());
}

TEST(TestUtil, CachedComputeScore_HillClimbingForLongTimeForProfiling) {
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
    auto new_placement_opt = suggest_random_position(problem, cache.m_solution, rnd, i);
    EXPECT_TRUE(bool(new_placement_opt)); // 99% sure..
    if (cache.change_musician(i, *new_placement_opt) < 0) {
      cache.change_musician(i, old_placement);
    }
  }

  LOG(INFO) << format("init  score = %s", int_to_delimited_string(init_score).c_str());
  LOG(INFO) << format("last  score = %s", int_to_delimited_string(cache.score()).c_str());
  LOG(INFO) << format("final score = %s", int_to_delimited_string(CachedComputeScore(problem).full_compute(cache.m_solution)).c_str());

  LOG(INFO) << format("full calculation %.4f ms/eval, partial update %.4f ms/eval",
    cache.get_mean_elapsed_ms_full(),
    cache.get_mean_elapsed_ms_partial());
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
TEST(TestUtil, ExtensionsInComputeScoreFunction) {
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
  EXPECT_DOUBLE_EQ(score_naive_full, 5357.0);
  // EXPECT_DOUBLE_EQ(score_naive_pillars, ?);
  EXPECT_DOUBLE_EQ(score_naive_harmony, 5357.0);
}

TEST(TestUtil, ExtensionsInComputeScoreFunctionWithSimpleProblem) {
  auto data = nlohmann::json::parse(R"({
  "room_width": 1000.0,
  "room_height": 1000.0,
  "stage_width": 1000.0,
  "stage_height": 200.0,
  "stage_bottom_left": [ 0.0, 700.0 ],
  "musicians": [0, 0],
  "attendees": [
    {
      "x": 450.0,
      "y": 450.0,
      "tastes": [1.0]
    }
  ],
  "pillars": [
    { "center": [ 500.0, 500.0 ], "radius": 50.0 }
  ]
})");
  Problem problem(data);
  problem.extension = Extension::full();

  Solution solution_blocked_by_pillar = {
      .placements = {
      {451.0, 700.0},
      {700.0, 700.0},
    },
  };
  EXPECT_NEAR(compute_score(problem, solution_blocked_by_pillar), 0.0, 1e-4);

  Solution solution_audible = {
      .placements = {
      {449.0, 700.0},
      {700.0, 700.0},
    },
  };
  EXPECT_NEAR(compute_score(problem, solution_audible),
   std::ceil((1.0 + 1.0 / (700 - 449.0)) * std::ceil(1e6 * 1.0 / SQ(700.0 - 450.0))),
   1e-4);
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
