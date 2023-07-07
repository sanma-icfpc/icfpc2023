#pragma once

static constexpr int k_musician_radius = 5;

struct Attendee {

    double x;
    double y;
    std::vector<double> tastes;

    std::string stringify() const {
        std::string res = format("Attendee {x=%f, y=%f, tastes=[", x, y);
        res += std::to_string(tastes[0]);
        for (int i = 1; i < tastes.size(); i++) {
            res += ", " + std::to_string(tastes[i]);
        }
        res += "]}";
        return res;
    }

};

struct Problem {

    double room_width;
    double room_height;
    double stage_w;
    double stage_h;
    double stage_x;
    double stage_y;
    std::vector<int> musicians;
    std::vector<Attendee> attendees;

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
    }

    static Problem from_file(int problem_id) {
        return from_file(format("../data/problems/problem-%d.json", problem_id));
    }

    static Problem from_file(std::string path) {
        std::ifstream ifs(path.c_str());
        nlohmann::json data;
        ifs >> data;
        return Problem(data);
    }

#ifdef HAVE_OPENCV_HIGHGUI
    cv::Mat_<cv::Vec3b> to_mat() const {
        int H = (int)room_height;
        int W = (int)room_width;
        cv::Mat_<cv::Vec3b> img(H, W, cv::Vec3b(155, 155, 155));
        cv::Rect stage(
            (int)stage_x,  // x
            (int)stage_y,    // y
            (int)stage_w,   // w
            (int)stage_h   // h
        );
        cv::rectangle(img, stage, cv::Scalar(255, 255, 255), cv::FILLED);
        for (const auto& a : attendees) {
            int y = (int)a.y;
            int x = (int)a.x;
            cv::circle(img, cv::Point(x, y), 2, cv::Scalar(0, 0, 255), cv::FILLED);
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

    std::string stringify() const { return format("(%f, %f)", x, y); }

};

struct Solution {

    std::vector<Placement> placements;

    std::string stringify() const {
        std::ostringstream oss;
        oss << placements;
        return oss.str();
    }

    static Solution from_file(std::string path) {
        std::ifstream ifs(path.c_str());
        nlohmann::json data;
        ifs >> data;

        Solution solution;
        for (const auto& placement : data["placements"]) {
            solution.placements.push_back({placement["x"], placement["y"]});
        }
        return solution;
    }

    nlohmann::json to_json() const {
        nlohmann::json data;
        data["placements"] = {};
        for (const auto& [x, y] : placements) {
            data["placements"].push_back({ {"x", x}, {"y", y} });
        }
        return data;
    }

#ifdef HAVE_OPENCV_HIGHGUI
    cv::Mat_<cv::Vec3b> to_mat(const Problem& problem) const {
        cv::Mat img = problem.to_mat();
        for (const auto& [x, y] : placements) {
            cv::circle(img, cv::Point(x, y), 2, cv::Scalar(0, 255, 0), cv::FILLED);
            cv::circle(img, cv::Point(x, y), k_musician_radius, cv::Scalar(0, 255, 0), 1);
        }
        return img;
    }
#endif

};