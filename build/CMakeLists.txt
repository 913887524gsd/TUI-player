cmake_minimum_required(VERSION 3.10)
project(TUI_player)

option(CMAKE_BUILD_TYPE "" release)

# C++
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -O0")

add_compile_options(-Wall -Wshadow -Werror -Wno-format)

set(PROJECT_DIR "${CMAKE_SOURCE_DIR}/..")

file(GLOB_RECURSE SRC_FILES
    "${PROJECT_DIR}/src/*.cpp"
)

include_directories("${PROJECT_DIR}/include")

set(DYN_LIB pthread asound)
set(PLAYER_NAME player)

add_executable(${PLAYER_NAME} ${SRC_FILES})
target_link_libraries(${PLAYER_NAME} ${DYN_LIB})