add_library( mystim SHARED
     imuInitialization.cpp 
     config.cpp
     Rotation.cpp
)
target_link_libraries( mystim
    ${THIRD_PARTY_LIBS}
)
