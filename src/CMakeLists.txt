project(Files)

macro(addToIt name)

    list(APPEND HEADER_FILES "${name}.h")

    if(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${name}.cpp")
        list(APPEND CPP_FILES "${name}.cpp")
    endif()

endmacro(addToIt)

MESSAGE("** CURRENT SOURCE INCLUDE DIRECTORY: ${CMAKE_CURRENT_SOURCE_DIR}")

set(HEADER_FILES "")
set(CPP_FILES "")

addToIt(Engine)

MESSAGE("** HEADER FILES : ${HEADER_FILES}")
MESSAGE("** SOURCE FILES : ${CPP_FILES}")

add_library(Files STATIC ${CPP_FILES} ${HEADER_FILES})