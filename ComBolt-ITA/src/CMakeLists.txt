ROOT_GENERATE_DICTIONARY(root_admin ROOT_admin.h LINKDEF LinkDef.h)
ROOT_GENERATE_DICTIONARY(particle_Dic particle_frzout.h particle_urqmd.h LINKDEF LinkDef_particle.h)

add_library(src SHARED 
    ROOT_admin.cpp
    root_admin.cxx
    particle_frzout.cpp
    particle_urqmd.cpp
    particle_Dic.cxx
    collision_kernel.h
    initial_state.h
    Administration.cpp
    # Butcher_tableaus.h
    grid.h
    BoltzEQ_solver.h
    # nr3.h
    # interpolation.h
    # interpolation_new.h
    Tmunu_and_hyd.h
    cornelius.cpp
    # rk4.h
    error_handler.h
    # surface_finder.h
    observables.h
    # stepperdopr5.h
        
    # ALGLIB files (only include what you need)
    # alglib-cpp/interpolation.cpp
    # alglib-cpp/ap.cpp
    # alglib-cpp/alglibinternal.cpp 
    # alglib-cpp/alglibmisc.cpp 
    # alglib-cpp/linalg.cpp 
    # alglib-cpp/solvers.cpp 
    # alglib-cpp/optimization.cpp 
    # alglib-cpp/specialfunctions.cpp 
    # alglib-cpp/integration.cpp  
    # alglib-cpp/stdafx.h
    )

target_link_libraries(src ${ROOT_LIBRARIES})
