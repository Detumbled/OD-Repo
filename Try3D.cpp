

#define GL_SILENCE_DEPRECATION
#include "utils/CSPICE/SpiceError.hpp"

#define GLFW_INCLUDE_GLCOREARB

// GLFW e ImGui
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl3.h"

// GLM per la matematica 3D
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <iostream>

// ---------------------------------------------------------
// 1. SHADER GLSL
// ---------------------------------------------------------
const char* vertexShaderSource = 
    "#version 330 core\n"
    "layout (location = 0) in vec3 aPos;\n"
    "uniform mat4 model;\n"
    "uniform mat4 view;\n"
    "uniform mat4 projection;\n"
    "void main() {\n"
    "   gl_Position = projection * view * model * vec4(aPos, 1.0);\n"
    "}\n";

const char* fragmentShaderSource = 
    "#version 330 core\n"
    "out vec4 FragColor;\n"
    "uniform vec3 lineColor;\n"
    "void main() {\n"
    "   FragColor = vec4(lineColor, 1.0f);\n"
    "}\n";

// Helper per compilare gli shader
GLuint CompileShader(GLuint type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    int success;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char infoLog[512];
        glGetShaderInfoLog(shader, 512, nullptr, infoLog);
        std::cerr << "ERRORE COMPILAZIONE SHADER\n" << infoLog << std::endl;
    }
    return shader;
}

int main() {
    // ---------------------------------------------------------
    // 2. CSPICE TRAJECTORY GENERATION
    // ---------------------------------------------------------
    const od::SpiceErrorModeGuard spice_error_mode("RETURN", "NONE");
    furnsh_c("../kernels.tm");
    od::ensureNoSpiceError("Loading kernels");

    SpiceDouble et_start, et_end;
    str2et_c("1978-02-01T00:00:00", &et_start); 
        str2et_c("1990-04-01T00:00:00", &et_end);
    od::ensureNoSpiceError("Time conversion");

    // OpenGL richiede array 1D continui di float [x1, y1, z1, x2, y2, z2...]
    std::vector<float> earth_verts, jupiter_verts, saturn_verts, vgr_verts;
    int steps = 20000;
    double step_size = (et_end - et_start) / steps;
    const double SCALE = 1000000.0; // Scaliamo di 1 milione di km per unità OpenGL

    for (int i = 0; i <= steps; ++i) {
        SpiceDouble current_et = et_start + (i * step_size);
        SpiceDouble s_earth[6], s_saturn[6], s_vgr[6], s_jupiter[6], lt;
        
        spkezr_c("EARTH", current_et, "J2000", "NONE", "SUN", s_earth, &lt);
        spkezr_c("5", current_et, "J2000", "NONE", "SUN", s_jupiter, &lt);
        spkezr_c("6", current_et, "J2000", "NONE", "SUN", s_saturn, &lt);
        spkezr_c("-31", current_et, "J2000", "NONE", "SUN", s_vgr, &lt);

        // Salviamo X, Y e Z (!)
        earth_verts.push_back(s_earth[0]/SCALE); earth_verts.push_back(s_earth[1]/SCALE); earth_verts.push_back(s_earth[2]/SCALE);
        jupiter_verts.push_back(s_jupiter[0]/SCALE); jupiter_verts.push_back(s_jupiter[1]/SCALE); jupiter_verts.push_back(s_jupiter[2]/SCALE);
        saturn_verts.push_back(s_saturn[0]/SCALE); saturn_verts.push_back(s_saturn[1]/SCALE); saturn_verts.push_back(s_saturn[2]/SCALE);
        vgr_verts.push_back(s_vgr[0]/SCALE); vgr_verts.push_back(s_vgr[1]/SCALE); vgr_verts.push_back(s_vgr[2]/SCALE);
    }

    std::cout << "Earth vertices: " << earth_verts.size() << '\n';
    std::cout << "Jupiter vertices: " << jupiter_verts.size() << '\n';
    std::cout << "Saturn vertices: " << saturn_verts.size() << '\n';
    std::cout << "Voyager vertices: " << vgr_verts.size() << '\n';

    std::cout
    << "Earth first point: "
    << earth_verts[0] << ", "
    << earth_verts[1] << ", "
    << earth_verts[2]
    << std::endl;
    kclear_c();

    // --- GENERAZIONE GRIGLIA ---
    std::vector<float> grid_verts;
    float gridSize = 2000.0f; // 2 Billion km (scaled)
    float gridStep = 200.0f;  // Lines every 200 Million km
    
    for (float i = -gridSize; i <= gridSize; i += gridStep) {
        // Linee parallele all'asse Y
        grid_verts.push_back(i); grid_verts.push_back(-gridSize); grid_verts.push_back(0.0f);
        grid_verts.push_back(i); grid_verts.push_back(gridSize); grid_verts.push_back(0.0f);
        // Linee parallele all'asse X
        grid_verts.push_back(-gridSize); grid_verts.push_back(i); grid_verts.push_back(0.0f);
        grid_verts.push_back(gridSize); grid_verts.push_back(i); grid_verts.push_back(0.0f);
    }
    // ---------------------------------------------------------
    // 3. OPENGL & GLFW INITIALIZATION
    // ---------------------------------------------------------
    if (!glfwInit()) return -1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Necessario su macOS

    GLFWwindow* window = glfwCreateWindow(1200, 900, "OD Engine 3D - Grand Tour", nullptr, nullptr);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    // ---------------------------------------------------------
    // 4. PREPARAZIONE SHADER E BUFFER GPU (VAO/VBO)
    // ---------------------------------------------------------
    GLuint vertexShader = CompileShader(GL_VERTEX_SHADER, vertexShaderSource);
    GLuint fragmentShader = CompileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    glDeleteShader(vertexShader); glDeleteShader(fragmentShader);

    // Funzione lambda per creare rapidamente i buffer in memoria video
    auto createBuffer = [](const std::vector<float>& data, GLuint& VAO, GLuint& VBO) {
        glGenVertexArrays(1, &VAO);
        glGenBuffers(1, &VBO);
        glBindVertexArray(VAO);
        glBindBuffer(GL_ARRAY_BUFFER, VBO);
        glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float), data.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);

        glBindVertexArray(0);
    };

    GLuint vaoEarth, vboEarth, vaoJup, vboJup, vaoSat, vboSat, vaoVgr, vboVgr;
    createBuffer(earth_verts, vaoEarth, vboEarth);
    createBuffer(jupiter_verts, vaoJup, vboJup);
    createBuffer(saturn_verts, vaoSat, vboSat);
    createBuffer(vgr_verts, vaoVgr, vboVgr);

    GLuint vaoGrid, vboGrid;
    createBuffer(grid_verts, vaoGrid, vboGrid);

    // Sole 

    // --- GENERAZIONE SOLE (Sfera Wireframe 3D) ---
    std::vector<float> sun_verts;
    float sunRadius = 40.0f; // Esageriamo il raggio visivamente (40 milioni di km virtuali)
    int stacks = 15; // Linee di latitudine
    int slices = 15; // Linee di longitudine
    const float PI = 3.14159265359f;

    for (int i = 0; i <= stacks; ++i) {
        float phi = PI * float(i) / float(stacks);
        for (int j = 0; j <= slices; ++j) {
            float theta = 2.0f * PI * float(j) / float(slices);
            
            // Punto corrente (x, y, z)
            float x = sunRadius * sin(phi) * cos(theta);
            float y = sunRadius * sin(phi) * sin(theta);
            float z = sunRadius * cos(phi);
            
            // Disegna linea orizzontale (verso il punto successivo lungo l'equatore)
            float thetaNext = 2.0f * PI * float(j + 1) / float(slices);
            float xNext = sunRadius * sin(phi) * cos(thetaNext);
            float yNext = sunRadius * sin(phi) * sin(thetaNext);
            sun_verts.push_back(x); sun_verts.push_back(y); sun_verts.push_back(z);
            sun_verts.push_back(xNext); sun_verts.push_back(yNext); sun_verts.push_back(z); // z resta uguale
            
            // Disegna linea verticale (verso il punto sottostante lungo il meridiano)
            if (i < stacks) {
                float phiNext = PI * float(i + 1) / float(stacks);
                float xDown = sunRadius * sin(phiNext) * cos(theta);
                float yDown = sunRadius * sin(phiNext) * sin(theta);
                float zDown = sunRadius * cos(phiNext);
                sun_verts.push_back(x); sun_verts.push_back(y); sun_verts.push_back(z);
                sun_verts.push_back(xDown); sun_verts.push_back(yDown); sun_verts.push_back(zDown);
            }
        }
    }

    // Carica il Sole nella GPU usando la lambda che hai già creato
    GLuint vaoSun, vboSun;
    createBuffer(sun_verts, vaoSun, vboSun);

    // Ottieni le locazioni delle variabili uniform nello shader
    GLint modelLoc = glGetUniformLocation(shaderProgram, "model");
    GLint viewLoc = glGetUniformLocation(shaderProgram, "view");
    GLint projLoc = glGetUniformLocation(shaderProgram, "projection");
    GLint colorLoc = glGetUniformLocation(shaderProgram, "lineColor");

    // ---------------------------------------------------------
    // 5. IMGUI INIT
    // ---------------------------------------------------------
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 150");

    // Variabili Camera
    float camRadius = 2500.0f; // Distanza
    float camYaw = 0.0f;       // Rotazione Z
    float camPitch = 0.5f;     // Elevazione

    glm::vec3 camTarget = glm::vec3(0.0f, 0.0f, 0.0f);

    glEnable(GL_DEPTH_TEST); // Abilita la profondità 3D (Z-Buffer)

    // ---------------------------------------------------------
    // 6. RENDER LOOP PRINCIPALE
    // ---------------------------------------------------------
    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // 6a. Interfaccia ImGui (Sovrapposta al 3D)
        ImGui_ImplOpenGL3_NewFrame(); ImGui_ImplGlfw_NewFrame(); ImGui::NewFrame();

        // 6a-bis. Legenda ImGui
        ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_FirstUseEver);
        ImGui::Begin("Legenda Sistema Solare", nullptr, ImGuiWindowFlags_AlwaysAutoResize);
        
        ImGui::TextColored(ImVec4(1.0f, 1.0f, 1.0f, 1.0f), "Corpo Celeste:");
        ImGui::Separator();
        ImGui::TextColored(ImVec4(0.2f, 0.6f, 1.0f, 1.0f), "Terra");
        ImGui::TextColored(ImVec4(0.8f, 0.6f, 0.3f, 1.0f), "Giove");
        ImGui::TextColored(ImVec4(0.9f, 0.8f, 0.5f, 1.0f), "Saturno");
        ImGui::TextColored(ImVec4(1.0f, 0.2f, 0.2f, 1.0f), "Voyager 1");
        ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.1f, 1.0f), "Sole");
        
        ImGui::End();

        // --- CONTROLLI MOUSE (Drag & Scroll) ---
        ImGuiIO& io = ImGui::GetIO();

        // Verifica che il mouse NON sia sopra un menu/legenda di ImGui
        if (!io.WantCaptureMouse) {
            // 1. ZOOM (Scroll)
            float scrollSensitivity = 150.0f;
            camRadius -= io.MouseWheel * scrollSensitivity;
            if (camRadius < 50.0f) camRadius = 50.0f;       // Limite minimo per non finire dentro il Sole
            if (camRadius > 10000.0f) camRadius = 10000.0f; // Limite massimo
            
            // 2. PAN/ORBIT (Drag con Tasto Sinistro)
            if (ImGui::IsMouseDown(ImGuiMouseButton_Left)) {
                float mouseSensitivity = 0.005f;
                // Usiamo il "Delta" (di quanti pixel si è mosso il mouse in questo frame)
                camYaw -= io.MouseDelta.x * mouseSensitivity;
                camPitch += io.MouseDelta.y * mouseSensitivity;

                // Blocchiamo l'elevazione a +/- 89 gradi per evitare il ribaltamento (Gimbal Lock)
                if (camPitch > 1.55f) camPitch = 1.55f;
                if (camPitch < -1.55f) camPitch = -1.55f;
            }
        }

        // 3. PAN / DRAG (Spostamento con il Tasto Destro)
            if (ImGui::IsMouseDown(ImGuiMouseButton_Right)) {
                // La sensibilità scala con lo zoom: se sei lontano, ti muovi più velocemente
                float panSensitivity = camRadius * 0.001f; 
                
                // Calcoliamo la posizione locale della camera rispetto al target
                glm::vec3 camOffset = glm::vec3(
                    camRadius * cos(camPitch) * sin(camYaw),
                    camRadius * cos(camPitch) * cos(camYaw),
                    camRadius * sin(camPitch)
                );
                
                // Calcoliamo i vettori direzionali "Destra" e "Su" relativi allo schermo
                glm::vec3 forward = glm::normalize(-camOffset);
                glm::vec3 worldUp = glm::vec3(0.0f, 0.0f, 1.0f);
                glm::vec3 right = glm::normalize(glm::cross(forward, worldUp));
                glm::vec3 up = glm::normalize(glm::cross(right, forward));
                
                // Spostiamo il target in base al movimento del mouse
                camTarget -= right * io.MouseDelta.x * panSensitivity;
                camTarget += up * io.MouseDelta.y * panSensitivity;
            }
   

        // 6b. Rendering OpenGL
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        glViewport(0, 0, width, height);
        glClearColor(0.02f, 0.02f, 0.02f, 1.0f); // Spazio scuro
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Pulisci colore e profondità

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glDisable(GL_CULL_FACE);
        glLineWidth(2.0f); // Rende le linee visibili sui display Mac Retina

        glUseProgram(shaderProgram);

        

        // Matematica della Camera GLM
        glm::mat4 projection = glm::perspective(glm::radians(45.0f), (float)width / (float)height, 0.1f, 100000.0f);
        
        glm::vec3 camOffset = glm::vec3(
                    camRadius * cos(camPitch) * sin(camYaw),
                    camRadius * cos(camPitch) * cos(camYaw),
                    camRadius * sin(camPitch)
                );
                
                // Posizione finale = Centro focale + Offset sferico
                glm::vec3 camPos = camTarget + camOffset;
                
                // Genera la matrice di vista puntando verso camTarget
                glm::mat4 view = glm::lookAt(camPos, camTarget, glm::vec3(0.0f, 0.0f, 1.0f));
        glm::mat4 model = glm::mat4(1.0f); // Identità

        glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(projection));
        glUniformMatrix4fv(viewLoc, 1, GL_FALSE, glm::value_ptr(view));
        glUniformMatrix4fv(modelLoc, 1, GL_FALSE, glm::value_ptr(model));

        glUniform3f(colorLoc, 0.15f, 0.15f, 0.15f); // Faint dark grey
        glBindVertexArray(vaoGrid);
        glDrawArrays(GL_LINES, 0, grid_verts.size() / 3);

        // --- DRAW SUN ---
        glUniform3f(colorLoc, 1.0f, 0.8f, 0.1f); // Colore: Giallo/Arancione brillante
        glBindVertexArray(vaoSun);
        glDrawArrays(GL_LINES, 0, sun_verts.size() / 3);

        // Funzione per disegnare un'orbita
        auto drawOrbit = [&](GLuint vao, float r, float g, float b, int numSteps) {
            glUniform3f(colorLoc, r, g, b);
            glBindVertexArray(vao);
            glDrawArrays(GL_LINE_STRIP, 0, numSteps + 1);
        };

        // Disegna tutto
        drawOrbit(vaoEarth, 0.2f, 0.6f, 1.0f, steps);  // Terra: Blu
        drawOrbit(vaoJup, 0.8f, 0.6f, 0.3f, steps);    // Giove: Arancione
        drawOrbit(vaoSat, 0.9f, 0.8f, 0.5f, steps);    // Saturno: Giallo pallido
        drawOrbit(vaoVgr, 1.0f, 0.2f, 0.2f, steps);    // Voyager 1: Rosso vivo

        // 6c. Renderizza ImGui sopra OpenGL
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown(); ImGui_ImplGlfw_Shutdown(); ImGui::DestroyContext();
    glfwDestroyWindow(window); glfwTerminate();
    return 0;
}
