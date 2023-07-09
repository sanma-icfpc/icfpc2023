#include "stdafx.h"

#include <glog/logging.h>
#include <gtest/gtest.h>
#include <ranges>

#include "../util.h"

TEST(TestUtil, GuessProblemId) {
  EXPECT_EQ(guess_problem_id("path\\to\\problem-15.json"), 15);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-16.json"), 16);
  EXPECT_EQ(guess_problem_id("path\\to\\solution-17_12345.json"), 17);
}

TEST(TestUtil, SolutionVolumeIO) {
  Solution solution_without_volume =
      Solution::from_file("../data/test_data/solution-example.json");
  std::vector<double> default_volumes = {1.0, 1.0, 1.0};
  EXPECT_EQ(solution_without_volume.volumes, default_volumes);

  Solution solution = Solution::from_file(
      "../data/test_data/solution-example-with-volume.json");
  std::vector<double> expected_volumes = {1.5, 5.5, 10.5};
  EXPECT_EQ(solution.volumes, expected_volumes);

  solution.volumes = {3.5, 4.5, 5.5};
  auto j = solution.to_json();
  EXPECT_TRUE(j.find("volumes") != j.end());
  if (j.contains("volumes")) {
    EXPECT_EQ(j["volumes"].size(), 3);
    if (j["volumes"].size() == 3) {
      EXPECT_EQ(j["volumes"][0], 3.5);
      EXPECT_EQ(j["volumes"][1], 4.5);
      EXPECT_EQ(j["volumes"][2], 5.5);
    }
  }
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

TEST(TestUtil, IsValidSolutionChecksVolume) {
  Problem problem = Problem::from_file(42);
  Solution solution = Solution::from_file(
      "../data/solutions/k3_v01_random_creation/solution-42.json");
  EXPECT_TRUE(is_valid_solution(problem, solution));
  solution.volumes[0] = -1.0;
  EXPECT_FALSE(is_valid_solution(problem, solution));
  solution.volumes[0] = 1.0;
  solution.volumes.pop_back();
  EXPECT_FALSE(is_valid_solution(problem, solution));
}

TEST(TestUtil, GuessExtension) {
  EXPECT_EQ(Problem::from_file(42).extension, Extension::lightning());
  EXPECT_EQ(Problem::from_file(56).extension, Extension::full());
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
      .placements =
          {
              {451.0, 700.0},
              {700.0, 700.0},
          },
  };
  EXPECT_NEAR(compute_score(problem, solution_blocked_by_pillar), 0.0, 1e-4);

  Solution solution_audible = {
      .placements =
          {
              {449.0, 700.0},
              {700.0, 700.0},
          },
  };
  EXPECT_NEAR(compute_score(problem, solution_audible),
              std::ceil((1.0 + 1.0 / (700 - 449.0)) *
                        std::ceil(1e6 * 1.0 / SQ(700.0 - 450.0))),
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
