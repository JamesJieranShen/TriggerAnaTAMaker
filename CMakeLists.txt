cmake_minimum_required(VERSION 3.16)
include(FetchContent)

project(TriggerActivityMaker LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
execute_process(
    COMMAND ${CMAKE_COMMAND} -E create_symlink
    "${CMAKE_BINARY_DIR}/compile_commands.json"
    "${CMAKE_SOURCE_DIR}/compile_commands.json"
)

# Build Options
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PROJECT_ROOT_DIR ${CMAKE_SOURCE_DIR})
add_compile_options(-Wall -Wextra -Wno-error=unused-parameter)
add_compile_options($<$<CONFIG:Release>:-march=native>)
set(default_build_type "RelWithDebInfo")

# Dependencies
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist TreePlayer Physics)

# If somebody (system package or another dependency) already provided it, reuse it.
if(NOT TARGET nlohmann_json::nlohmann_json)
  # Prefer an installed package if available (avoids duplicate target definitions).
  find_package(nlohmann_json 3.12.0 CONFIG QUIET)
endif()

if(NOT TARGET nlohmann_json::nlohmann_json)
  # Optional: reduce side effects from the subproject
  set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
  set(JSON_Install    OFF CACHE BOOL "" FORCE)  # suggested by upstream docs :contentReference[oaicite:1]{index=1}

  FetchContent_Declare(
    nlohmann_json
    GIT_REPOSITORY https://github.com/nlohmann/json.git
    GIT_TAG        v3.12.0
  )
  FetchContent_MakeAvailable(nlohmann_json)
endif()

# CLI11
fetchcontent_declare(
  CLI11
  GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
  GIT_TAG        v2.5.0
)
fetchcontent_makeavailable(CLI11)

add_subdirectory(src)
