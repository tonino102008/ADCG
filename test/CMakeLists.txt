add_definitions(-DLOGGING)

add_executable(test1 test1.cpp)
target_link_directories(test1 PUBLIC ${LIBS_DIR})
target_include_directories(test1 PUBLIC ${INCLUDE_DIR})
target_link_libraries(test1 PUBLIC -lm)
target_compile_options(test1 PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
target_include_directories(test1 PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
target_link_libraries(test1 PUBLIC $<LINK_ONLY:MKL::MKL>)
