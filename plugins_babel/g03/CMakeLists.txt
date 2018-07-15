INSTALL(PROGRAMS babel_g03_in.sh  babel_g03_out.sh DESTINATION bin)
file(APPEND ${TMP_CONFIG_DIR}/plugin.rc "  libpath=${CMAKE_INSTALL_PREFIX}/bin  extension=log  description=Gaussian03  forward=babel_g03_out.sh  backward=babel_g03_in.sh  \n")

#file(COPY babel_g03_in.sh  babel_g03_out.sh DESTINATION bin)

