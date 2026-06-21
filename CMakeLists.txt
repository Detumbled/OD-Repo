cmake_minimum_required(VERSION 3.20)
project(coding_v1 LANGUAGES C CXX) 

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
enable_testing()

# ------------------------------
# 1. Configuration & Standards
# ------------------------------
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(APPLE)
    add_compile_options("-gdwarf-4")
endif()

# ------------------------------
# 2. Dependencies (Eigen, GLFW, GLEW, OpenGL)
# ------------------------------
find_package(Eigen3 CONFIG REQUIRED)
find_package(glfw3 CONFIG REQUIRED)
find_package(GLEW REQUIRED)
find_package(OpenGL REQUIRED)
find_package(Ceres CONFIG REQUIRED)

message(STATUS "Graphics Libraries Found:")
message(STATUS "  GLFW: ${glfw3_DIR}")
message(STATUS "  GLEW: ${GLEW_INCLUDE_DIRS}")
message(STATUS "  OpenGL: ${OPENGL_LIBRARIES}")

# ------------------------------
# 3. Dependencies (CSPICE)
# ------------------------------
if(NOT DEFINED ENV{CSPICE_HOME})
    message(FATAL_ERROR "Errore: CSPICE_HOME non definita!")
endif()

set(CSPICE_DIR "$ENV{CSPICE_HOME}")

find_path(CSPICE_INCLUDE_DIR
    NAMES SpiceUsr.h
    PATHS "${CSPICE_DIR}/include" "${CSPICE_DIR}/cspice/include"
    NO_DEFAULT_PATH
)

find_library(CSPICE_LIB 
    NAMES cspice cspice.a
    PATHS "${CSPICE_DIR}/lib" "${CSPICE_DIR}/cspice/lib"
    NO_DEFAULT_PATH
)

# ------------------------------
# 4. ★ NEW: ImGui + ImPlot
# ------------------------------
set(IMGUI_DIR  ${CMAKE_SOURCE_DIR}/third_party/imgui)
set(IMPLOT_DIR ${CMAKE_SOURCE_DIR}/third_party/implot)

add_library(imgui_lib STATIC
    ${IMGUI_DIR}/imgui.cpp
    ${IMGUI_DIR}/imgui_draw.cpp
    ${IMGUI_DIR}/imgui_tables.cpp
    ${IMGUI_DIR}/imgui_widgets.cpp
    ${IMGUI_DIR}/backends/imgui_impl_glfw.cpp
    ${IMGUI_DIR}/backends/imgui_impl_opengl3.cpp
    ${IMPLOT_DIR}/implot.cpp
    ${IMPLOT_DIR}/implot_items.cpp
)

target_include_directories(imgui_lib PUBLIC
    ${IMGUI_DIR}
    ${IMGUI_DIR}/backends
    ${IMPLOT_DIR}
)

# ImGui needs GLFW + OpenGL — reuse what you already found
target_link_libraries(imgui_lib PUBLIC
    glfw
    ${OPENGL_LIBRARIES}
)

if(APPLE)
    target_link_libraries(imgui_lib PUBLIC
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
endif()

# ------------------------------
# 5. Target: project_dependencies (Interface)
# ------------------------------
add_library(project_dependencies INTERFACE)

target_include_directories(project_dependencies INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    "${CSPICE_INCLUDE_DIR}"
    ${GLEW_INCLUDE_DIRS}
    "${CMAKE_SOURCE_DIR}/third_party/glm"
)

target_link_libraries(project_dependencies INTERFACE 
    Eigen3::Eigen 
    "${CSPICE_LIB}"
    glfw
    GLEW::GLEW
    ${OPENGL_LIBRARIES}
    imgui_lib          
    Ceres::ceres
    m
)

if(APPLE)
    target_link_libraries(project_dependencies INTERFACE 
        "-framework Cocoa"
        "-framework IOKit"
        "-framework CoreVideo"
    )
endif()

# ------------------------------
# 6. Static Library (Astro functions)
# ------------------------------
set(BIG_FUNCTIONS_SOURCES
    src/dummy.cpp
    src/hi.cpp
)

add_library(big_functions STATIC ${BIG_FUNCTIONS_SOURCES})
target_link_libraries(big_functions PUBLIC project_dependencies)

# ------------------------------
# 7. Executables
# ------------------------------
set(APP_SOURCES
    Try3D.cpp
    example_4_2_1.cpp
    main.cpp
    test_stations.cpp
)

foreach(app_src IN LISTS APP_SOURCES)
    get_filename_component(app_name "${app_src}" NAME_WE)
    add_executable(${app_name} "${app_src}")
    target_link_libraries(${app_name} PRIVATE big_functions)

    if(app_name MATCHES "^test_")
        add_test(NAME ${app_name} COMMAND ${app_name} WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    endif()
endforeach()
