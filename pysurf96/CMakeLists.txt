# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

set(module_name "surfdisp96_ext")
set(module_src "${CMAKE_CURRENT_SOURCE_DIR}/surfdisp96.f")
set(module_c "${CMAKE_CURRENT_BINARY_DIR}/${module_name}module.c")

# add extra f77 compile flags
set(F77_FLAGS "-O3 -ffixed-line-length-none -fbounds-check -m64")
set_source_files_properties(${module_src} PROPERTIES COMPILE_FLAGS "${F77_FLAGS}")

add_custom_command(
  OUTPUT ${module_c}
  DEPENDS ${module_src}
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          ${module_src} -m ${module_name}
          --lower only: surfdisp96 :)

python_add_library(${module_name} MODULE ${module_c} ${module_src} WITH_SOABI)
target_link_libraries(${module_name} PRIVATE fortranobject)

install(TARGETS ${module_name} DESTINATION pysurf96)
