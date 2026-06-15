find_path(PETSC_DIR NAMES include/petsc.h
    PATHS ENV PETSC_DIR
    DOC "PETSc installation directory"
)

if(PETSC_DIR)
    set(PETSC_INCLUDE_DIRS
        ${PETSC_DIR}/include
        ${PETSC_DIR}/${PETSC_ARCH}/include
    )

    find_path(PETSC_CONF_DIR NAMES petscconf.h
        PATHS ${PETSC_DIR}/${PETSC_ARCH}/include
              ${PETSC_DIR}/include
        NO_DEFAULT_PATH
    )

    find_library(PETSC_LIBRARY
        NAMES petsc
        PATHS ${PETSC_DIR}/${PETSC_ARCH}/lib
              ${PETSC_DIR}/lib
        NO_DEFAULT_PATH
    )

    if(PETSC_LIBRARY)
        set(PETSC_LIBRARIES ${PETSC_LIBRARY})
        
        file(STRINGS "${PETSC_INCLUDE_DIRS}/petscversion.h" PETSC_VERSION_HEADER
             REGEX "#define PETSC_VERSION_(MAJOR|MINOR|SUBMINOR|RELEASE) "
             LIMIT_COUNT 4)
        
        foreach(comp MAJOR MINOR SUBMINOR RELEASE)
            string(REGEX REPLACE ".*#define PETSC_VERSION_${comp} +([0-9]+).*" "\\1"
                   PETSC_VERSION_${comp} "${PETSC_VERSION_HEADER}")
        endforeach()
        
        set(PETSC_VERSION "${PETSC_VERSION_MAJOR}.${PETSC_VERSION_MINOR}.${PETSC_VERSION_SUBMINOR}")
    endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PETSc
    REQUIRED_VARS PETSC_DIR PETSC_INCLUDE_DIRS PETSC_LIBRARIES
    VERSION_VAR PETSC_VERSION
)

if(PETSC_FOUND AND NOT TARGET PETSc::PETSc)
    add_library(PETSc::PETSc UNKNOWN IMPORTED)
    set_target_properties(PETSc::PETSc PROPERTIES
        IMPORTED_LOCATION "${PETSC_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${PETSC_INCLUDE_DIRS}"
    )
endif()

mark_as_advanced(PETSC_DIR PETSC_INCLUDE_DIRS PETSC_LIBRARIES PETSC_VERSION)
