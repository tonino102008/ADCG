set(INCLUDE_DIR "../include")
set(LIBS_DIR "../lib")

add_definitions(-DLOGGING)

add_executable(SingleTrack SingleTrack.cpp)
target_link_directories(SingleTrack PUBLIC ${LIBS_DIR})
target_include_directories(SingleTrack PUBLIC ${INCLUDE_DIR})
target_link_libraries(SingleTrack PUBLIC -lm)
target_compile_options(SingleTrack PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(SingleTrack PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(SingleTrack PUBLIC $<LINK_ONLY:MKL::MKL>)
