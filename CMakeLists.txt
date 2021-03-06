
# Set this to the minimal version you want to support
find_package(raylib 3.0 QUIET) # Let CMake search for a raylib-config.cmake


# You could change the QUIET above to REQUIRED and remove this if() clause
# This part downloads raylib and builds it if it's not installed on your system
if (NOT raylib_FOUND) # If there's none, fetch and build raylib
  include(FetchContent)

  FetchContent_Declare(
    raylib
    URL https://github.com/raysan5/raylib/archive/master.tar.gz
  )

  FetchContent_GetProperties(raylib)
  if (NOT raylib_POPULATED) # Have we downloaded raylib yet?
    set(FETCHCONTENT_QUIET NO)
    FetchContent_Populate(raylib)

    set(BUILD_EXAMPLES OFF CACHE BOOL "" FORCE) # don't build the supplied examples

    if (EMSCRIPTEN)
        add_definitions("-pthread -s USE_PTHREADS=1")
    endif()
    # build raylib
    add_subdirectory(${raylib_SOURCE_DIR} ${raylib_BINARY_DIR})

  endif()

endif()


set(TINYGLTF_BUILD_LOADER_EXAMPLE OFF)
add_definitions(-DTINYGLTF_NO_STB_IMAGE -DTINYGLTF_NO_STB_IMAGE_WRITE)

add_subdirectory("3rdparty/glm")
add_subdirectory("3rdparty/tinygltf")
add_subdirectory("3rdparty/longeronpp")

project(ornigame)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)


file (GLOB_RECURSE C_FILES src/*.cpp)
file (GLOB_RECURSE H_FILES src/*.h)
set (SOURCE_FILES ${C_FILES} ${H_FILES})
add_executable(ornigame ${C_FILES})


target_link_libraries(ornigame raylib longeron tinygltf glm Threads::Threads)

if (EMSCRIPTEN)
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}")
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread -s USE_PTHREADS=1 -s PTHREAD_POOL_SIZE=2 -s USE_GLFW=3 -s ASSERTIONS=1 -s WASM=1 -s ASYNCIFY -s ALLOW_MEMORY_GROWTH=1 --preload-file ${PROJECT_SOURCE_DIR}/bin@./")
    set(CMAKE_EXECUTABLE_SUFFIX ".html")
endif ()
