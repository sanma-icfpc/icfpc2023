#pragma once
#include <stdafx.h>
#include <filesystem>
#include <regex>
#include <fstream>
#include <nlohmann/json.hpp>
#include <sstream>
#include <string>
#include <vector>
#include <optional>

static constexpr int k_musician_radius = 5;
static constexpr int k_musician_spacing_radius = 10;

struct Attendee {
  double x;
  double y;
  std::vector<double> tastes;

  std::string stringify() const {
    std::string res = "Attendee {x=" + std::to_string(x) +
                      ", y=" + std::to_string(y) + ", tastes=[";
    res += std::to_string(tastes[0]);
    for (int i = 1; i < tastes.size(); i++) {
      res += ", " + std::to_string(tastes[i]);
    }
    res += "]}";
    return res;
  }
};

struct Pillar {
  double x;
  double y;
  double r;
};

struct Extension {
  bool consider_pillars = false; // Full Division, Extension 1: Obstacles in the Room
  bool consider_harmony = false; // Full Division, Extension 2: Playing Together
  auto operator<=>(const Extension&) const = default;
  std::string stringify() const {
    return "{.consider_pillars=" + std::to_string(consider_pillars) + ", .consider_harmony=" + std::to_string(consider_harmony) + "}";
  }
  static Extension from_problem_id(int problem_id) {
    if (problem_id <= 55) {
      return lightning();
    }
    return full();
  }
  static constexpr Extension lightning() { return Extension {false, false}; }
  static constexpr Extension full() { return Extension {true, true}; }
};

inline std::optional<int> guess_problem_id(std::string some_file_path) {
  std::regex pattern(".*-(\\d+)");
  std::smatch matches;

  if (std::regex_search(some_file_path, matches, pattern)) {
    std::string numberString = matches[1].str();
    try {
      return stoi(numberString);
    } catch (const std::invalid_argument& e) {
      LOG(ERROR) << "Invalid argument: " << e.what();
    } catch (const std::out_of_range& e) {
      LOG(ERROR) << "Out of range: " << e.what();
    }
  }
  return std::nullopt;
}

struct Problem {
  double room_width;
  double room_height;
  double stage_w;
  double stage_h;
  double stage_x;
  double stage_y;
  std::vector<int> musicians;
  std::vector<Attendee> attendees;
  std::vector<Pillar> pillars;
  Extension extension;
  int problem_id = -1;

  Problem(const nlohmann::json& data) {
    room_width = data["room_width"];
    room_height = data["room_height"];
    stage_w = data["stage_width"];
    stage_h = data["stage_height"];
    stage_x = data["stage_bottom_left"][0];
    stage_y = data["stage_bottom_left"][1];
    for (const auto& type : data["musicians"]) {
      musicians.push_back(type);
    }
    for (const auto& attendee_json : data["attendees"]) {
      Attendee attendee;
      attendee.x = attendee_json["x"];
      attendee.y = attendee_json["y"];
      for (const auto& taste : attendee_json["tastes"]) {
        attendee.tastes.push_back(taste);
      }
      attendees.push_back(attendee);
    }
    for (auto&& pillar_json : data["pillars"]) {
      Pillar pillar;
      auto&& center = pillar_json["center"];
      pillar.x = center[0];
      pillar.y = center[1];
      pillar.r = pillar_json["radius"];
      pillars.push_back(pillar);
    }
  }

  static Problem from_file(int problem_id) {
    // Linux: cwd = /src
    // Windows: cwd = /vs
    std::filesystem::path json_path("../data/problems/problem-" + std::to_string(problem_id) + ".json");
    LOG_ASSERT(std::filesystem::is_regular_file(json_path));
    return Problem::from_file(json_path.string());
  }

  static Problem from_file(std::string path) {
    std::ifstream ifs(path.c_str());
    nlohmann::json data;
    ifs >> data;
    Problem problem(data);
    if (auto problem_id_opt = guess_problem_id(path)) {
      problem.extension = Extension::from_problem_id(*problem_id_opt);
      problem.problem_id = *problem_id_opt;
    } else {
      LOG(WARNING) << "Failed to guess problem ID! extension may be wrong. assuming full division..";
      problem.extension = Extension::full();
      problem.problem_id = -1;
    }
    return problem;
  }

#ifdef HAVE_OPENCV_HIGHGUI
  cv::Mat_<cv::Vec3b> to_mat() const {
    int H = (int)room_height;
    int W = (int)room_width;
    cv::Mat_<cv::Vec3b> img(H, W, cv::Vec3b(155, 155, 155));
    cv::Rect stage((int)stage_x,  // x
                   (int)stage_y,  // y
                   (int)stage_w,  // w
                   (int)stage_h   // h
    );
    cv::rectangle(img, stage, cv::Scalar(255, 255, 255), cv::FILLED);
    for (const auto& a : attendees) {
      int y = (int)a.y;
      int x = (int)a.x;
      cv::circle(img, cv::Point(x, y), 2, cv::Scalar(0, 0, 255), cv::FILLED);
    }
    for (const auto& p : pillars) {
      int x = static_cast<int>(p.x);
      int y = static_cast<int>(p.y);
      int r = static_cast<int>(p.r);
      cv::circle(img, cv::Point(x, y), r, cv::Scalar(128, 0, 0), cv::FILLED);
    }
    return img;
  }
  cv::Mat_<cv::Vec3b> show(int delay = 0) const {
    cv::Mat_<cv::Vec3b> img = to_mat();
    cv::imshow("img", img);
    cv::waitKey(delay);
    return img;
  }
#endif
};

struct Placement {
  double x;
  double y;

  Placement(double x = 0.0, double y = 0.0) : x(x), y(y) {}
  auto operator<=>(const Placement&) const = default;

  std::string stringify() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
  }
};

struct Solution {
  std::vector<Placement> placements;
  std::vector<double> volumes;

  void set_default_volumes() {
      volumes.assign(placements.size(), 1.0);
  }

  std::string stringify() const {
    if (placements.empty())
      return "{}";
    std::string res = "{" + placements[0].stringify();
    for (int i = 1; i < (int)placements.size(); i++) {
      res += ", " + placements[i].stringify();
    }
    res += "}, {" + std::to_string(volumes[0]);
    for (int i = 1; i < (int)volumes.size(); i++) {
      res += ", " + std::to_string(volumes[i]);
    }
    res += "}";
    return res;
  }

  static Solution from_file(std::string path) {
    std::ifstream ifs(path.c_str());
    if (!ifs.is_open()) {
      return {};
    }

    nlohmann::json data;
    ifs >> data;

    Solution solution;
    for (const auto& placement : data["placements"]) {
      solution.placements.push_back({placement["x"], placement["y"]});
    }
    if (data.contains("volumes")) {
      for (const auto& volume : data["volumes"]) {
        solution.volumes.push_back(volume);
      }
    } else {
      solution.set_default_volumes();
    }
    return solution;
  }

  nlohmann::json to_json() const {
    nlohmann::json data;
    data["placements"] = {};
    for (const auto& [x, y] : placements) {
      data["placements"].push_back({{"x", x}, {"y", y}});
    }
    if (!volumes.empty()) {
      if (volumes.size() != placements.size()) {
        LOG(WARNING) << "volumes " << volumes.size() << " and placements " << placements.size() << " size mismatch!";
      }
      data["volumes"] = nlohmann::json::array();
      for (const auto volume : volumes) {
        data["volumes"].push_back(volume);
      }
    }
    return data;
  }

#ifdef HAVE_OPENCV_HIGHGUI
  cv::Mat_<cv::Vec3b> to_mat(const Problem& problem) const {
    cv::Mat img = problem.to_mat();
    for (const auto& [x, y] : placements) {
      cv::circle(img, cv::Point(x, y), 2, cv::Scalar(0, 255, 0), cv::FILLED);
      cv::circle(img, cv::Point(x, y), k_musician_radius, cv::Scalar(0, 255, 0),
                 1);
    }
    return img;
  }
#endif
};