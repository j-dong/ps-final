#include "interface.h"

constexpr int WINDOW_W = WIDTH * 20, WINDOW_H = HEIGHT * 20;

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <cstdio>

#include <iostream>
#include <thread>
#include <vector>

#define ERROR(...) ::std::fprintf(stderr, __VA_ARGS__)

std::atomic_bool running;

extern void sim_main();
extern void sim_init_grid(Grid *grid);

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

void main() {
    vec4 data = texture(tex, texcoord);
    vec3 col3 = mix(vec3(0.2, 0.2, 0.2), vec3(1.0, 0.0, 0.0), data.z / 40.0);
    color = vec4(col3 * data.y, data.y);
    color.r = data.x;
    color.g = data.y;
    color.b = data.z;
    color.a = 1.0;
}
)zzz";

static const float VERTICES[][2] = {
    { -1.0, -1.0 },
    {  3.0, -1.0 },
    { -1.0,  3.0 },
};

static constexpr float CLEAR_COLOR[] = {0.7, 0.7, 0.7, 1.0};

float upload_data[HEIGHT][WIDTH][4];

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

    sim_init_grid(get_init_grid());
    init_grids();
    std::thread sim_thread(sim_main);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, WIDTH, HEIGHT,
                 0, GL_RGBA, GL_FLOAT,
                 nullptr);

    while (!glfwWindowShouldClose(window)) {
        Grid *grid = grid_swap(READER);
        GLenum err;
        if (!grid->updated) goto next;
        grid->updated = false;

        // std::cout << "Got new data!" << std::endl;
        // printf("using %p\n", grid);

        gen_upload_data(grid);

        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex);
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

        glClear(GL_COLOR_BUFFER_BIT);
        glDrawArrays(GL_TRIANGLES, 0, 3);
        glfwSwapBuffers(window);

next:
        glfwPollEvents();
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
