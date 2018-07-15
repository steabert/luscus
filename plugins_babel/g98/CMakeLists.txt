INSTALL(PROGRAMS babel_g98_in.sh  babel_g98_out.sh DESTINATION bin)

#file(COPY babel_g98_in.sh  babel_g98_out.sh DESTINATION bin)
file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin  extension=log  description=Gaussian98  forward=babel_g98_out.sh  backward=babel_g98_in.sh  \n")
