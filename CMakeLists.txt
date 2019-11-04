cmake_minimum_required(VERSION 2.8)
project(MPCHonestMajorityNoTriplesYehuda)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11 -maes -mpclmul -mbmi2")

INCLUDE_DIRECTORIES($ENV{HOME} $ENV{HOME}/libscapi/install/include
        $ENV{HOME}/libscapi/install/include/OTExtensionBristol)

set(SOURCE_FILES main.cpp
        ProtocolParty.cpp ProtocolTimer.cpp GF2_8LookupTable.cpp TemplateField.cpp HashEncrypt.cpp)


link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib )


add_executable(MPCHonestMajorityNoTriplesYehuda ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(MPCHonestMajorityNoTriplesYehuda ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp gmpxx
        boost_system boost_thread pthread crypto dl ssl z)
