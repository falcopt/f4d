#ifndef _F4D_RENDERER_HPP_
#define _F4D_RENDERER_HPP_

#include <GLFW/glfw3.h>

#include "cobra/Instance.hpp"
#include "cobra/MoveGenerators.hpp"
#include "cobra/Solution.hpp"

#define FRAME_WIDTH (10.0)
#define TRAJECTORY_SECTION_HEIGHT (300.0)
#define SOLUTION_SECTION_HEIGHT (700.0)


#define XMIN (0.0)
#define XMAX (1400.0)

#define TRAJECTORY_SECTION_Y_BEGIN (FRAME_WIDTH)
#define TRAJECTORY_SECTION_Y_END (TRAJECTORY_SECTION_Y_BEGIN + TRAJECTORY_SECTION_HEIGHT)
#define TRAJECTORY_SECTION_X_BEGIN (XMIN + FRAME_WIDTH)
#define TRAJECTORY_SECTION_X_END (XMAX * 0.66 - FRAME_WIDTH / 2.0)

#define ZOOMED_TRAJECTORY_SECTION_Y_BEGIN (TRAJECTORY_SECTION_Y_BEGIN)
#define ZOOMED_TRAJECTORY_SECTION_Y_END (TRAJECTORY_SECTION_Y_END)
#define ZOOMED_TRAJECTORY_SECTION_X_BEGIN (TRAJECTORY_SECTION_X_END + FRAME_WIDTH)
#define ZOOMED_TRAJECTORY_SECTION_X_END (XMAX - FRAME_WIDTH)

#define SOLUTION_SECTION_Y_BEGIN (ZOOMED_TRAJECTORY_SECTION_Y_END + FRAME_WIDTH)
#define SOLUTION_SECTION_Y_END (SOLUTION_SECTION_Y_BEGIN + SOLUTION_SECTION_HEIGHT)
#define SOLUTION_SECTION_X_BEGIN (XMIN + FRAME_WIDTH)
#define SOLUTION_SECTION_X_END (XMAX / 2.0 - FRAME_WIDTH / 2.0)

#define OMEGA_SECTION_Y_BEGIN (SOLUTION_SECTION_Y_BEGIN)
#define OMEGA_SECTION_Y_END (SOLUTION_SECTION_Y_END)
#define OMEGA_SECTION_X_BEGIN (SOLUTION_SECTION_X_END + FRAME_WIDTH)
#define OMEGA_SECTION_X_END (XMAX - FRAME_WIDTH)


#define YMIN (0)
#define YMAX (SOLUTION_SECTION_Y_END + FRAME_WIDTH)


class Renderer {

public:
    Renderer(const cobra::Instance& instance_, float initial_cost_, const std::vector<int>& omega_)
        : instance(instance_), initial_cost(initial_cost_), omega(omega_) {

        if (!glfwInit()) {
            std::cerr << "Cannot successfully execute 'glfwInit'\n";
            abort();
        }

        window = glfwCreateWindow(XMAX, YMAX, "FILO", nullptr, nullptr);

        if (!window) {
            std::cerr << "Cannot successfully execute 'glfwCreateWindow'\n";
            abort();
        }

        glfwSetWindowCloseCallback(window, [](GLFWwindow*) -> void {
            exit(1);
        });

        glfwSetFramebufferSizeCallback(window, [](GLFWwindow* wind, int width, int height) {
            glfwMakeContextCurrent(wind);
            glViewport(0, 0, width, height);
        });

        glfwMakeContextCurrent(window);

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);

        glfwPollEvents();

        x_min = instance.get_x_coordinate(instance.get_vertices_begin());
        x_max = x_min;
        y_min = instance.get_y_coordinate(instance.get_vertices_begin());
        y_max = y_min;
        for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {
            const auto x_i = instance.get_x_coordinate(i);
            const auto y_i = instance.get_y_coordinate(i);
            x_min = std::min(x_min, x_i);
            x_max = std::max(x_max, x_i);
            y_min = std::min(y_min, y_i);
            y_max = std::max(y_max, y_i);
        }

        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable(GL_BLEND);
        glEnable(GL_LINE_SMOOTH);

        add_trajectory_point(initial_cost, initial_cost, initial_cost, true, initial_cost);
    }

    virtual ~Renderer() {
        glfwDestroyWindow(window);
    }

    void add_trajectory_point(float shaken_solution_cost, float local_optimum_cost, float current_solution_cost, bool is_current_solution_feasible,
                              float best_solution_cost) {

        auto trajectory_point = TrajectoryPoint();

        trajectory_point.current_solution_gap = 100.0f * (current_solution_cost - initial_cost) / initial_cost;
        trajectory_point.is_current_solution_feasible = is_current_solution_feasible;
        trajectory_point.best_solution_gap = 100.0f * (best_solution_cost - initial_cost) / initial_cost;
        trajectory_point.shaken_solution_gap = 100.0f * (shaken_solution_cost - initial_cost) / initial_cost;
        trajectory_point.local_optima_gap = 100.0f * (local_optimum_cost - initial_cost) / initial_cost;
        trajectory.emplace_back(trajectory_point);

        min_gap = std::min({trajectory_point.best_solution_gap, min_gap, trajectory_point.current_solution_gap, trajectory_point.shaken_solution_gap,
                            trajectory_point.local_optima_gap});

        max_gap = std::max({trajectory_point.local_optima_gap, max_gap, trajectory_point.current_solution_gap, trajectory_point.shaken_solution_gap,
                            trajectory_point.best_solution_gap});
    }

    void draw(cobra::Solution& solution, const cobra::LRUCache& cached_vertices, cobra::MoveGenerators& move_generators) {

        glfwMakeContextCurrent(window);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();

        glOrtho(XMIN, XMAX, YMIN, YMAX, -1, 1);

        glMatrixMode(GL_MODELVIEW);

        glClearColor(0.14, 0.15, 0.16, 1.0);
        glClear(GL_COLOR_BUFFER_BIT);

        draw_hseparator(YMIN, FRAME_WIDTH);

        draw_trajectory();

        draw_zoomed_trajectory();

        draw_hseparator(ZOOMED_TRAJECTORY_SECTION_Y_END, SOLUTION_SECTION_Y_BEGIN);

        draw_solution(solution, cached_vertices, move_generators);
        draw_omega(cached_vertices);

        draw_hseparator(SOLUTION_SECTION_Y_END, YMAX);

        draw_vseparator(XMIN, XMIN + FRAME_WIDTH);
        draw_vseparator(XMAX - FRAME_WIDTH, XMAX);
        draw_vseparator(SOLUTION_SECTION_X_END, OMEGA_SECTION_X_BEGIN, SOLUTION_SECTION_Y_BEGIN, SOLUTION_SECTION_Y_END);
        draw_vseparator(TRAJECTORY_SECTION_X_END, ZOOMED_TRAJECTORY_SECTION_X_BEGIN, TRAJECTORY_SECTION_Y_BEGIN, TRAJECTORY_SECTION_Y_END);

        glfwSwapBuffers(window);

        glfwPollEvents();
    }

    void draw_hseparator(double begin, double end) {

        glBegin(GL_QUADS);
        glColor4f(0.18, 0.19, 0.20, 1);
        glVertex2d(XMIN, begin);
        glVertex2d(XMAX, begin);
        glVertex2d(XMAX, end);
        glVertex2d(XMIN, end);
        glEnd();

        glColor4f(0.25, 0.26, 0.27, 1);
        glLineWidth(1);

        glBegin(GL_LINES);
        glVertex2d(XMIN, begin);
        glVertex2d(XMAX, begin);
        glEnd();

        glBegin(GL_LINES);
        glVertex2d(XMIN, end);
        glVertex2d(XMAX, end);
        glEnd();
    }

    void draw_vseparator(double xmin, double xmax, double ymin = YMIN, double ymax = YMAX) {

        glBegin(GL_QUADS);
        glColor4f(0.18, 0.19, 0.20, 1);
        glVertex2d(xmin, ymin);
        glVertex2d(xmax, ymin);
        glVertex2d(xmax, ymax);
        glVertex2d(xmin, ymax);
        glEnd();

        glColor4f(0.25, 0.26, 0.27, 1);
        glLineWidth(1);

        glBegin(GL_LINES);
        glVertex2d(xmin, ymin);
        glVertex2d(xmin, ymax);
        glEnd();

        glBegin(GL_LINES);
        glVertex2d(xmax, ymin);
        glVertex2d(xmax, ymax);
        glEnd();
    }

    void draw_trajectory() {

        auto xs = Scaler(0, trajectory.size(), TRAJECTORY_SECTION_X_BEGIN, TRAJECTORY_SECTION_X_END);
        auto ys = Scaler(min_gap - 0.5f, max_gap + 1, TRAJECTORY_SECTION_Y_BEGIN, TRAJECTORY_SECTION_Y_END);

        glColor4f(1.0f, 1.0f, 1.0f, 0.25f);
        glPushAttrib(GL_ENABLE_BIT);
        glLineStipple(5, 0xAAAA);
        glEnable(GL_LINE_STIPPLE);
        glBegin(GL_LINES);
        for (auto i = std::min(static_cast<int>(std::round(min_gap - 0.5f)), 1); i <= max_gap + 1; i++) {
            glVertex2f(xs.scale(0), ys.scale(i));
            glVertex2f(xs.scale(trajectory.size()), ys.scale(i));
        }
        glEnd();
        glPopAttrib();


        glColor4f(1.0f, 0.0f, 0.0f, 0.50f);
        glBegin(GL_POINTS);
        for (auto x = 0u; x < trajectory.size(); x++) {
            const auto y = trajectory[x].shaken_solution_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
        }
        glEnd();


        glColor4f(0.0f, 1.0f, 0.0f, 0.50f);
        glBegin(GL_POINTS);
        for (auto x = 0u; x < trajectory.size(); x++) {
            const auto y = trajectory[x].local_optima_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
        }
        glEnd();


        glBegin(GL_LINES);
        for (auto x = 0u; x < trajectory.size() - 1; x++) {
            if (trajectory[x].is_current_solution_feasible) {
                glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
            } else {
                glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
            }
            const auto y = trajectory[x].current_solution_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
            glVertex2f(xs.scale(x + 1), ys.scale(trajectory[x + 1].current_solution_gap));
        }
        glEnd();

        glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
        glBegin(GL_LINES);
        for (auto x = 0u; x < trajectory.size() - 1; x++) {
            const auto y = trajectory[x].best_solution_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
            glVertex2f(xs.scale(x + 1), ys.scale(trajectory[x + 1].best_solution_gap));
        }
        glEnd();
    }

    void draw_zoomed_trajectory() {

        const auto zoom_iter = std::min<int>(1000, trajectory.size());

        const auto beg = std::max<size_t>(0, static_cast<int>(trajectory.size()) - zoom_iter - 1);

        auto zoom_min_gap = std::numeric_limits<float>::max();
        auto zoom_max_gap = std::numeric_limits<float>::lowest();

        const auto min_gap_value = 0.01f;

        for (auto n = beg; n < trajectory.size(); n++) {
            if (trajectory[n].best_solution_gap < zoom_min_gap) {
                zoom_min_gap = trajectory[n].best_solution_gap;
            }
            if (trajectory[n].current_solution_gap > zoom_max_gap) {
                zoom_max_gap = trajectory[n].current_solution_gap;
            }
        }

        auto xs = Scaler(beg, trajectory.size(), ZOOMED_TRAJECTORY_SECTION_X_BEGIN, ZOOMED_TRAJECTORY_SECTION_X_END);
        auto ys = Scaler(zoom_min_gap - min_gap_value, zoom_max_gap + min_gap_value, ZOOMED_TRAJECTORY_SECTION_Y_BEGIN, ZOOMED_TRAJECTORY_SECTION_Y_END);

        glColor4f(1.0f, 1.0f, 1.0f, 0.25f);
        glPushAttrib(GL_ENABLE_BIT);
        glLineStipple(5, 0xAAAA);
        glEnable(GL_LINE_STIPPLE);
        glBegin(GL_LINES);
        auto g = zoom_min_gap - min_gap_value;
        while (g + min_gap_value < zoom_max_gap + min_gap_value) {
            g += min_gap_value;  // gap increment
            glVertex2f(xs.scale(beg), ys.scale(g));
            glVertex2f(xs.scale(trajectory.size()), ys.scale(g));
        }
        glEnd();
        glPopAttrib();

        glBegin(GL_LINES);
        for (auto x = beg; x < trajectory.size() - 1; x++) {
            if (trajectory[x].is_current_solution_feasible) {
                glColor4f(1.0f, 1.0f, 0.0f, 1.0f);
            } else {
                glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
            }
            const auto y = trajectory[x].current_solution_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
            glVertex2f(xs.scale(x + 1), ys.scale(trajectory[x + 1].current_solution_gap));
        }
        glEnd();

        glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
        glBegin(GL_LINES);
        for (auto x = beg; x < trajectory.size() - 1; x++) {
            const auto y = trajectory[x].best_solution_gap;
            glVertex2f(xs.scale(x), ys.scale(y));
            glVertex2f(xs.scale(x + 1), ys.scale(trajectory[x + 1].best_solution_gap));
        }
        glEnd();
    }

    void draw_solution(cobra::Solution& solution, const cobra::LRUCache& cached_vertices, cobra::MoveGenerators& move_generators) {

        auto xs = Scaler(x_min, x_max, SOLUTION_SECTION_X_BEGIN, SOLUTION_SECTION_X_END);
        auto ys = Scaler(y_min, y_max, SOLUTION_SECTION_Y_BEGIN, SOLUTION_SECTION_Y_END);

        glLineWidth(1);
        for (auto i = instance.get_vertices_begin(); i < instance.get_vertices_end(); i++) {

            for (auto move_id : move_generators.get_move_generator_indices_involving(i)) {

                const auto& move = move_generators.get(move_id);

                glBegin(GL_LINES);
                glColor4f(0.0, 0.43, 0.75, 0.10);
                glVertex2f(xs.scale(instance.get_x_coordinate(move.get_first_vertex())), ys.scale(instance.get_y_coordinate(move.get_first_vertex())));
                glVertex2f(xs.scale(instance.get_x_coordinate(move.get_second_vertex())), ys.scale(instance.get_y_coordinate(move.get_second_vertex())));
                glEnd();
            }
        }


        for (auto i = cached_vertices.begin(); i != cached_vertices.end(); i = cached_vertices.get_next(i)) {

            for (auto move_id : move_generators.get_move_generator_indices_involving(i)) {
                const auto& move = move_generators.get(move_id);
                glBegin(GL_LINES);
                // glColor4f(1.0, 1.0, 0.0, 0.50);
                glColor4f(0.07, 0.60, 0.98, 0.20);
                glVertex2d(xs.scale(instance.get_x_coordinate(move.get_first_vertex())), ys.scale(instance.get_y_coordinate(move.get_first_vertex())));
                glVertex2d(xs.scale(instance.get_x_coordinate(move.get_second_vertex())), ys.scale(instance.get_y_coordinate(move.get_second_vertex())));
                glEnd();
            }
        }


        glLineWidth(2);
        for (auto route = solution.get_first_route(); route != cobra::Solution::dummy_route; route = solution.get_next_route(route)) {

            const auto load_ratio = static_cast<float>(solution.get_route_load(route)) / static_cast<float>(instance.get_vehicle_capacity());

            glBegin(GL_LINES);

            glColor4f(load_ratio, 1 - load_ratio, 0.00, 1.0f);

            auto curr = solution.get_first_customer(route);
            auto next = curr;

            do {

                next = solution.get_next_vertex(route, curr);

                if (next == instance.get_depot()) {
                    break;
                }

                glVertex2f(xs.scale(instance.get_x_coordinate(curr)), ys.scale(instance.get_y_coordinate(curr)));
                glVertex2f(xs.scale(instance.get_x_coordinate(next)), ys.scale(instance.get_y_coordinate(next)));

                curr = next;


            } while (true);


            glEnd();

            /*glPushAttrib(GL_ENABLE_BIT);
            glLineStipple(1, 0xAAAA);
            glEnable(GL_LINE_STIPPLE);
            glColor4f(load_ratio,1-load_ratio,0.0, 0.1f);
            glBegin(GL_LINES);
            glVertex2d(instance.get_x_coordinate(solution.get_last_customer(route)), instance.get_y_coordinate(solution.get_last_customer(route)));
            glVertex2d(instance.get_x_coordinate(instance.get_depot()), instance.get_y_coordinate(instance.get_depot()));
            glEnd();
            glPopAttrib();*/
        }

        glColor4f(1.0, 0.75, 0.00, 1.0f);
        glPointSize(20);
        glBegin(GL_POINTS);
        glVertex2d(xs.scale(instance.get_x_coordinate(instance.get_depot())), ys.scale(instance.get_y_coordinate(instance.get_depot())));
        glEnd();
        glPointSize(1);
    }


    void draw_omega(const cobra::LRUCache& cached_vertices) {

        auto omega_max = std::numeric_limits<float>::lowest();

        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {
            if (omega[i] > omega_max) {
                omega_max = omega[i];
            }
        }

        auto xs = Scaler(x_min, x_max, OMEGA_SECTION_X_BEGIN, OMEGA_SECTION_X_END);
        auto ys = Scaler(y_min, y_max, OMEGA_SECTION_Y_BEGIN, OMEGA_SECTION_Y_END);

        for (auto i = instance.get_customers_begin(); i < instance.get_customers_end(); i++) {

            const auto value = static_cast<float>(omega[i]) / static_cast<float>(omega_max);

            glColor4f(value, 0, 1 - value, 1);
            glPointSize(3);
            glBegin(GL_POINTS);
            glVertex2d(xs.scale(instance.get_x_coordinate(i)), ys.scale(instance.get_y_coordinate(i)));
            glEnd();
        }

        glColor4f(1.0, 0.75, 0.00, 1.0f);
        glPointSize(20);
        glBegin(GL_POINTS);
        glVertex2d(xs.scale(instance.get_x_coordinate(instance.get_depot())), ys.scale(instance.get_y_coordinate(instance.get_depot())));
        glEnd();
        glPointSize(1);


        for (auto i = cached_vertices.begin(); i != cached_vertices.end(); i = cached_vertices.get_next(i)) {

            glColor4f(1, 1, 1, 1);
            glPointSize(3);
            glBegin(GL_POINTS);
            glVertex2d(xs.scale(instance.get_x_coordinate(i)), ys.scale(instance.get_y_coordinate(i)));
            glEnd();
        }

        glPointSize(1);
    }

private:
    GLFWwindow* window = nullptr;

    const cobra::Instance& instance;
    float initial_cost;
    const std::vector<int>& omega;

    struct TrajectoryPoint {
        float shaken_solution_gap;
        float local_optima_gap;
        float current_solution_gap;
        bool is_current_solution_feasible;
        float best_solution_gap;
    };

    float min_gap = std::numeric_limits<float>::max();
    float max_gap = std::numeric_limits<float>::min();

    std::vector<TrajectoryPoint> trajectory;

    float x_min{}, x_max{}, y_min{}, y_max{};

    class Scaler {
    public:
        Scaler() = default;

        void set(double in_min_, double in_max_, double out_min_, double out_max_) {
            in_min = in_min_;
            in_max = in_max_;
            out_min = out_min_;
            out_max = out_max_;
            scaled_min = scale(in_min);
            scaled_max = scale(in_max);
        }

        Scaler(double in_min_, double in_max_, double out_min_, double out_max_) {
            set(in_min_, in_max_, out_min_, out_max_);
        }

        double scale(double v) {
            return (v - in_min) * (out_max - out_min) / (in_max - in_min) + out_min;
        }

        double min() {
            return scaled_min;
        }

        double max() {
            return scaled_max;
        }

    private:
        double in_min = 0, in_max = 10, out_min = 0, out_max = 10, scaled_min = 0, scaled_max = 10;
    };
};

#endif