cmake_minimum_required(VERSION 3.24)
project(LeedsGL)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

file(GLOB cpp_files
        "LeedsGLRenderWindow/*.cpp"
        )
list(FILTER cpp_files EXCLUDE REGEX ".*(/RenderWidget\\.cpp)")

function (ListToString result delim)
    list(GET ARGV 2 temp)
    math(EXPR N "${ARGC}-1")
    foreach(IDX RANGE 3 ${N})
        list(GET ARGV ${IDX} STR)
        set(temp "${temp}${delim}${STR}")
    endforeach()
    set(${result} "${temp}" PARENT_SCOPE)
endfunction(ListToString)

ListToString(str "\n" ${cpp_files})
message(STATUS "${str}")


add_executable(LeedsGL ${cpp_files})

if (CMAKE_HOST_APPLE)

    set(CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt@5/lib/cmake" "/usr/local/lib/cmake")
    target_link_libraries(LeedsGL "-framework OpenGL")
    set_target_properties(LeedsGL PROPERTIES LINK_FLAGS "-Wl,-F/Library/Frameworks")

endif()

find_package(Qt5 COMPONENTS
        Core
        Gui
        Widgets
        OpenGL
        REQUIRED)

target_link_libraries(LeedsGL
        Qt5::Core
        Qt5::Gui
        Qt5::Widgets
        Qt5::OpenGL
        )
