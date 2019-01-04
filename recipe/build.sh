#export CFLAGS="-Wall -g -m64 -pipe -O2  -fPIC $CFLAGS"
#export CPPFLAGS="-I${PREFIX}/include"
#export LDFLAGS="-L${PREFIX}/lib"

MAKEFILE=ezget_Makefile.gfortran

sed "s#@cdat_EXTERNALS@#${PREFIX}#g;" ${MAKEFILE}.in > ${MAKEFILE}
make  -f ${MAKEFILE}
make -f ${MAKEFILE} install
