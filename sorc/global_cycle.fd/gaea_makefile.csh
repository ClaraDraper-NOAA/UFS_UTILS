#-----------------------------------------------------
#-use standard module.
#-----------------------------------------------------
module load cray-netcdf
module load cmake/3.17.0
module load cray-hdf5
setenv NCEPLIBS /lustre/f2/pdata/ncep_shared/NCEPLIBS/lib/
module use -a $NCEPLIBS/modulefiles
module load nemsio-intel-sandybridge/2.2.3
module load w3emc-intel-sandybridge/2.2.0
module load w3nco-intel-sandybridge/2.0.6
module load bacio-intel-sandybridge/2.0.1
module load sp-intel-sandybridge/2.0.2
module load sigio-intel-sandybridge/2.0.1
module load sfcio-intel-sandybridge/1.0.0
module load crtm-intel-sandybridge/2.2.5
module load bufr-intel-sandybridge/11.0.1
module load ip-intel-sandybridge/3.0.0
module list

setenv FCMP ftn

##setenv DEBUG='-ftrapuv -check all -check nooutput_conversion -fp-stack-check -fstack-protector -traceback -g'
setenv INCS "-I$IP_INCd -I${CRAY_NETCDF_PREFIX}/include"
setenv FFLAGS "$INCS -O3 -fp-model precise -r8 -convert big_endian -traceback -g"
setenv OMPFLAG -qopenmp
setenv LDFLG -qopenmp

setenv LIBSM "${W3NCO_LIBd} ${BACIO_LIB4} ${IP_LIBd} ${SP_LIBd} -L/${CRAY_NETCDF_PREFIX}/lib -lnetcdf -lnetcdff -lhdf5"

make -f Makefile_gaea clean
make -f Makefile_gaea
