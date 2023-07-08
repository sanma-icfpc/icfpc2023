#include "stdafx.h"
#include <nlohmann/json.hpp>
#include <fmt/format.h>
#include <CLI/CLI.hpp>
#include <omp.h>

#include "spec.h"
#include "util.h"

int main(int argc, char* argv[]) {
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();
  google::SetStderrLogging(google::INFO);
  google::SetLogDestination(google::INFO, "main.log.");

  std::ios::sync_with_stdio(false);
  std::cin.tie(NULL);

  CLI::App app{ "main module" };
  app.require_subcommand(0, 1);

  auto sub_problem_to_png = app.add_subcommand("problem-to-png");
  std::string problem_file;
  std::string output_file;
  sub_problem_to_png->add_option("problem_file", problem_file, "problem file path (JSON)");
  sub_problem_to_png->add_option("output_file", output_file, "output file path (PNG)");

  auto sub_solution_to_png = app.add_subcommand("solution-to-png");
  std::string solution_file;
  sub_solution_to_png->add_option("solution_file", solution_file, "solution file path (JSON)");
  sub_solution_to_png->add_option("output_file", output_file, "output file path (PNG)");

  CLI11_PARSE(app, argc, argv);

#ifndef HAVE_OPENCV_HIGHGUI
    LOG(ERROR) << "no OpenCV highgui!";
#else
  if (sub_problem_to_png->parsed()) {
    Problem problem = Problem::from_file(problem_file);
    cv::Mat img = problem.to_mat();
    cv::imwrite(output_file, img);
  }
  if (sub_solution_to_png->parsed()) {
    if (auto problem_id = guess_problem_id(solution_file)) {
      Problem problem = Problem::from_file(*problem_id);
      Solution solution = Solution::from_file(solution_file);
      cv::Mat img = solution.to_mat(problem);
      cv::imwrite(output_file, img);
    } else {
      LOG(ERROR) << "Could not guess the problem id!";
      return 1;
    }
  }
#endif


  return 0;
}
