cmake_minimum_required(VERSION 3.8)
project(miruSAXS)

set(CMAKE_CXX_STANDARD 11)

find_library(SASTools_LIB libSASTools.a PATHS ~/usr/local/lib ~/lib)

include_directories(~/usr/local/include)

#
# BOOST SECTION
#
# set Boost directory
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_STATIC_RUNTIME ON)
# Boost root variable set as environment variable in CLION Preferences
MESSAGE("** BOOST ROOT => ${BOOSTROOT}")
SET(BOOST_ROOT ${BOOSTROOT})
SET(BOOST_LIBRARYDIR ${BOOSTROOT}/stage/lib)
SET(BOOST_INCLUDEDIR ${BOOSTROOT})
include_directories(${BOOST_INCLUDEDIR})

set(BOOST_FILESYSTEM ${BOOSTROOT}/stage/lib/libboost_filesystem.a )
set(BOOST_PO ${BOOSTROOT}/stage/lib/libboost_program_options.a )
set(BOOST_RANDOM ${BOOSTROOT}/stage/lib/libboost_random.a )
set(BOOST_REGEX ${BOOSTROOT}/stage/lib/libboost_regex.a )
set(BOOST_SYSTEM ${BOOSTROOT}/stage/lib/libboost_system.a )

MESSAGE("** BOOST INCLUDE => ${BOOST_INCLUDEDIR}")
MESSAGE("** BOOST SPECIAL => ${BOOST_FILESYSTEM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_PO}")
MESSAGE("** BOOST SPECIAL => ${BOOST_RANDOM}")
MESSAGE("** BOOST SPECIAL => ${BOOST_REGEX}")
MESSAGE("** BOOST SPECIAL => ${BOOST_SYSTEM}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${Boost_INCLUDE_DIR}")
MESSAGE("** BOOST INCLUDE DIRECTORY: ${BOOST_INCLUDEDIR}")


set(SOURCE_FILES main.cpp)
include_directories(src)
add_subdirectory(src)

SET (QHULL_INCLUDE_DIRS ~/src_libs/qhull-master/src ~/src_libs/qhull-master/src/libqhull ~/src_libs/qhull-master/src/libqhullcpp)
include_directories(${QHULL_INCLUDE_DIRS})
find_library(QHULL_LIB libqhullstatic.a PATHS ~/src_libs/qhull-master/lib)
find_library(QHULLCPP_LIB libqhullcpp.a PATHS ~/src_libs/qhull-master/lib)

add_executable(miruSAXS ${SOURCE_FILES})

target_link_libraries(miruSAXS
        PUBLIC
        ${BOOST_FILESYSTEM} ${BOOST_PO} ${BOOST_RANDOM} ${BOOST_REGEX} ${BOOST_SYSTEM}
        ${SASTools_LIB}  ${QHULL_LIB} ${QHULLCPP_LIB} Files
        )