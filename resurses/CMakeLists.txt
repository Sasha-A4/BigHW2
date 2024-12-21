set(HEADERS
        fixed.hpp
        fast_fixed.hpp
        fixed_oper.hpp
)

set(SOURCES
        ../main.cpp
)

add_executable(main ${SOURCES} ${HEADERS})

target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})