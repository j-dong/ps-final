#include "interface.h"

constexpr int WINDOW_W = WIDTH * 10, WINDOW_H = HEIGHT * 10;

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <cstdio>

#include <iostream>
#include <thread>
#include <vector>

#define ERROR(...) ::std::fprintf(stderr, __VA_ARGS__)

std::atomic_bool running;

extern void sim_main();
extern void sim_init_grid(Grid *grid);
extern void sim_init();

static bool want_reset = false;

static const char *VERTEX_SHADER = R"zzz(#version 330 core
in vec2 c;
out vec2 texcoord;
void main() {
    texcoord = 0.5 * (1.0 + c);
    gl_Position = vec4(c, 0.0, 1.0);
}
)zzz";

static const char *FRAGMENT_SHADER = R"zzz(#version 330 core
in vec2 texcoord;
out vec4 color;

uniform sampler2D tex;

uniform vec3 pressure_color, density_color, temperature_color;
uniform vec3 fire_coeff;
uniform float fire_scale;
uniform bool debug_mode;

void main() {
    vec4 data = texture(tex, texcoord);
    // vec3 col3 = mix(vec3(0.2, 0.2, 0.2), vec3(1.0, 0.0, 0.0), data.z / 40.0);
    // color = vec4(col3 * data.y, data.y);
    if (debug_mode) {
        color = vec4(vec3(0.0) + pressure_color * data.x + density_color * data.y + temperature_color * data.z, 1.0);
    } else {
        // smoke absorption + diffuse
        color = vec4(0.6 * vec3(1.0 - clamp(data.y - 1.0, 0.0, 0.2)), 1.0) * clamp(data.y, 0.0, 1.0);
        // emission
        color += vec4(exp(fire_scale * fire_coeff * max(0.0, data.z)) - 1.0, 0.0);
    }
}
)zzz";

static const float VERTICES[][2] = {
    { -1.0, -1.0 },
    {  3.0, -1.0 },
    { -1.0,  3.0 },
};

static constexpr float CLEAR_COLOR[] = {0.0, 0.0, 0.0, 1.0};

float upload_data[HEIGHT][WIDTH][4];

float pressure_color[3]    = {0.0f, 0.0f, 1.0f},
      density_color[3]     = {1.0f, 1.0f, 0.0f},
      temperature_color[3] = {0.0f, 0.0f, 0.0f};

float fire_coeff[3] = {1.0, 0.6, 0.3};
float fire_scale    = 0.02;

bool debug_graphics = false;

static void gen_upload_data(const Grid *grid) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // TODO: colormapping
            upload_data[y][x][0] = (float) grid->pressure[y][x];
            upload_data[y][x][1] = (float) grid->density[y][x];
            upload_data[y][x][2] = (float) grid->temperature[y][x];
            upload_data[y][x][3] = 0.0f;
        }
    }
}

void render_gui() {
    SimParams *params = param_buf.swap(WRITER);
    static bool show_metrics = false, show_demo = false;
    ImGui::Begin("Parameters");
    if (ImGui::Button("Reset Simulation")) {
        want_reset = true;
    }
    if (ImGui::Button("Export Velocities")) {
        params->want_to_export = true;
    }
    ImGui::InputDouble("Timestep", &params->timestep);
    if (ImGui::CollapsingHeader("Buoyancy", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputDouble("Alpha", &params->alpha);
        ImGui::InputDouble("Beta", &params->beta);
    }
    if (ImGui::CollapsingHeader("Vorticity Confinement", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputDouble("Epsilon", &params->epsilon);
    }
    if (ImGui::CollapsingHeader("Visualization", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Debug Mode", &debug_graphics);
        ImGui::ColorEdit3("Fire Coefficients", fire_coeff);
        ImGui::InputFloat("Fire Scale", &fire_scale);
        ImGui::ColorEdit3("Pressure Color", pressure_color);
        ImGui::ColorEdit3("Density Color", density_color);
        ImGui::ColorEdit3("Temperature Color", temperature_color);
    }
    if (ImGui::CollapsingHeader("Emitter", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::InputDouble("Density", &params->emitter_density);
        ImGui::InputDouble("Temperature", &params->emitter_temp);
    }
    if (ImGui::CollapsingHeader("Debug", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Checkbox("Show Metrics", &show_metrics);
        ImGui::Checkbox("Show Demo", &show_demo);
    }
    ImGui::End();
    if (show_metrics) ImGui::ShowMetricsWindow();
    if (show_demo) ImGui::ShowDemoWindow();
    params->updated = true;
}

int main() {
    if (!glfwInit()) {
        ERROR("error initializing GLFW\n");
        return 1;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SRGB_CAPABLE, GLFW_TRUE);
    GLFWwindow *window = glfwCreateWindow(
        WINDOW_W, WINDOW_H,
        "Gas Simulation",
        nullptr, nullptr
    );
    if (!window) {
        ERROR("error creating window\n");
        return 1;
    }
    glfwMakeContextCurrent(window);
    gladLoadGLLoader((GLADloadproc) glfwGetProcAddress);
    glfwSwapInterval(1);

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLint temp;

    GLuint vs, fs;
    vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &VERTEX_SHADER, nullptr);
    glCompileShader(vs);
    glGetShaderiv(vs, GL_COMPILE_STATUS, &temp);
    if (temp != GL_TRUE) {
        glGetShaderiv(vs, GL_INFO_LOG_LENGTH, &temp);
        std::vector<GLchar> log;
        log.resize((size_t) temp);
        glGetShaderInfoLog(vs, temp, nullptr, log.data());
        ERROR("error compiling vertex shader:\n%s\n", (const char *) log.data());
        return 1;
    }
    fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &FRAGMENT_SHADER, nullptr);
    glCompileShader(fs);
    glGetShaderiv(fs, GL_COMPILE_STATUS, &temp);
    if (temp != GL_TRUE) {
        glGetShaderiv(fs, GL_INFO_LOG_LENGTH, &temp);
        std::vector<GLchar> log;
        log.resize((size_t) temp);
        glGetShaderInfoLog(fs, temp, nullptr, log.data());
        ERROR("error compiling fragment shader:\n%s\n", (const char *) log.data());
        return 1;
    }

    GLuint prog;
    prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glDeleteShader(vs);
    glDeleteShader(fs);
    glLinkProgram(prog);
    glDetachShader(prog, vs);
    glDetachShader(prog, fs);
    glGetProgramiv(prog, GL_LINK_STATUS, &temp);
    if (temp != GL_TRUE) {
        glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &temp);
        std::vector<GLchar> log;
        log.resize((size_t) temp);
        glGetProgramInfoLog(prog, temp, nullptr, log.data());
        ERROR("error linking program:\n%s\n", (const char *) log.data());
        return 1;
    }

    glUseProgram(prog);

    GLint loc_tex = glGetUniformLocation(prog, "tex");
    glUniform1i(loc_tex, 0);

    GLint loc_col_p = glGetUniformLocation(prog, "pressure_color"),
          loc_col_d = glGetUniformLocation(prog, "density_color"),
          loc_col_t = glGetUniformLocation(prog, "temperature_color");

    GLint loc_fire_coeff = glGetUniformLocation(prog, "fire_coeff"),
          loc_fire_scale = glGetUniformLocation(prog, "fire_scale"),
          loc_debug = glGetUniformLocation(prog, "debug_mode");

    GLuint tex;
    glGenTextures(1, &tex);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);

    GLuint buf;
    glGenBuffers(1, &buf);
    glBindBuffer(GL_ARRAY_BUFFER, buf);
    glBufferData(GL_ARRAY_BUFFER, sizeof VERTICES, (const void *) VERTICES, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(0);

    glClearColor(CLEAR_COLOR[0],
                 CLEAR_COLOR[1],
                 CLEAR_COLOR[2],
                 CLEAR_COLOR[3]);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_FRAMEBUFFER_SRGB);

    running.store(true, std::memory_order_relaxed);

    sim_init_grid(grids.get_init());
    grids.init();
    sim_init();
    std::thread sim_thread(sim_main);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, WIDTH, HEIGHT,
                 0, GL_RGBA, GL_FLOAT,
                 nullptr);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO &io = ImGui::GetIO(); (void) io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 150");

    bool valid = false;

    while (!glfwWindowShouldClose(window)) {
        Grid *grid = grids.swap(READER);
        GLenum err;
        if (!grid->updated) goto next;
        grid->updated = false;
        valid = true;

        gen_upload_data(grid);

        glGetError();
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, WIDTH, HEIGHT,
                     GL_RGBA, GL_FLOAT, (const GLvoid *) upload_data);
        err = glGetError();
        switch (err) {
        case GL_NO_ERROR:
            break;
        case GL_INVALID_ENUM:
            std::cerr << "error: GL_INVALID_ENUM" << std::endl;
            break;
        case GL_INVALID_VALUE:
            std::cerr << "error: GL_INVALID_VALUE" << std::endl;
            break;
        case GL_INVALID_OPERATION:
            std::cerr << "error: GL_INVALID_OPERATION" << std::endl;
            break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            std::cerr << "error: GL_INVALID_FRAMEBUFFER_OPERATION" << std::endl;
            break;
        case GL_OUT_OF_MEMORY:
            std::cerr << "error: GL_OUT_OF_MEMORY" << std::endl;
            break;
        default:
            std::cerr << "error: unknown error" << std::endl;
            break;
        }

next:
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        render_gui();
        ImGui::Render();

        glClear(GL_COLOR_BUFFER_BIT);
        glUniform3fv(loc_col_p, 1, pressure_color);
        glUniform3fv(loc_col_d, 1, density_color);
        glUniform3fv(loc_col_t, 1, temperature_color);
        glUniform3fv(loc_fire_coeff, 1, fire_coeff);
        glUniform1f(loc_fire_scale, fire_scale);
        glUniform1i(loc_debug, debug_graphics);
        if (valid) { glDrawArrays(GL_TRIANGLES, 0, 3); }

        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);

        glfwPollEvents();

        if (want_reset) {
            want_reset = false;
            std::cout << "Resetting simulation!" << std::endl;
            running.store(false, std::memory_order_relaxed);
            sim_thread.join();
            sim_init_grid(grids.get_init());
            grids.init();
            running.store(true, std::memory_order_relaxed);
            sim_thread = std::thread(sim_main);
        }
    }

    running.store(false, std::memory_order_relaxed);

    glDeleteBuffers(1, &buf);
    glDeleteTextures(1, &tex);
    glDeleteProgram(prog);

    glfwDestroyWindow(window);
    glfwTerminate();

    sim_thread.join();

    return 0;
}
