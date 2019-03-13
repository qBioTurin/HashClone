cmake_minimum_required(VERSION 2.8)

project(hashcheck)

add_executable(hashcheck main.cpp classIO.cpp general.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -Wall  -fmessage-length=0")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -Wall  -fmessage-length=0")
endif()
install(TARGETS hashcheck RUNTIME DESTINATION bin)

project(freqchecker)


add_executable(freqchekcer mainfreq.cpp general.cpp)

 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall  -fmessage-length=0")
 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -g -Wall  -fmessage-length=0")

install(TARGETS freqchekcer RUNTIME DESTINATION bin)



add_executable(parser mainparser.cpp )

 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall  -fmessage-length=0")
 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -g -Wall  -fmessage-length=0")

install(TARGETS parser RUNTIME DESTINATION bin)

add_executable(testfreq test_frequency.cpp )

 set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -Wall  -fmessage-length=0")
 set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -g -Wall  -fmessage-length=0")

install(TARGETS testfreq RUNTIME DESTINATION bin)



project(createKmer)

add_executable(createKmer mainKmer.cpp classIO.cpp general.cpp)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -O3 -Wall  -fmessage-length=0")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -Wall  -fmessage-length=0")
endif()
install(TARGETS createKmer RUNTIME DESTINATION bin)
