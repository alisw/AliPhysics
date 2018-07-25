<<<<<<< HEAD
prefix=@CMAKE_INSTALL_PREFIX@
exec_prefix=${prefix}
includedir=${prefix}/@INCLUDE_INSTALL_ROOT_DIR@
libdir=${exec_prefix}/@LIB_INSTALL_DIR@
=======
libdir=@LIB_INSTALL_DIR@
includedir=@INCLUDE_INSTALL_ROOT_DIR@
>>>>>>> origin/master

Name: Yaml-cpp
Description: A YAML parser and emitter for C++
Version: @YAML_CPP_VERSION@
Requires:
Libs: -L${libdir} -lyaml-cpp
Cflags: -I${includedir}
