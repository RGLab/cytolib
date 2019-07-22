#untar the lib
BASEPBNAME="hdf5"
PBTGZNAME=hdf5small_cxx_1.10.3.tar.gz
if test -d ./${BASEPBNAME}; then 
	echo 'found ' $BASEPBNAME ' header sources and tar archive;using what is there.'
else
	echo "untarring $PBTGZNAME ...";
	gunzip -dc ${PBTGZNAME} | tar xf -;
fi;

echo "building the szip library...";
cd ${BASEPBNAME}/szip;
./configure --with-pic --enable-shared=no
make


cd ../..
echo "building the hdf5 library...";
H5BUILD=$(pwd)/h5_build
if test -d ${H5BUILD}; then
	echo 'found ' $H5BUILD ' ;using what is there.'
else
	mkdir ${H5BUILD}
fi
cd ${BASEPBNAME};
./configure --with-pic --enable-shared=no --enable-cxx \
    --with-szlib \
   LDFLAGS="-L`pwd -P`/szip/src/.libs" \
    CPPFLAGS="-I`pwd -P`/szip/src/" \
    --prefix="${H5BUILD}" --libdir="${H5BUILD}/lib"
    
make lib
#somehow the faulty makefile from c++/ failed to install c++ lib and header 
#make install

#cp by hand
cd ..
USER_INCLUDE=${H5BUILD}/include
USER_LIB_DIR=${H5BUILD}/lib
HDF5_INCLUDE="${BASEPBNAME}/src"
HDF5_CXX_INCLUDE="${BASEPBNAME}/c++/src"
HDF5_LIB="${BASEPBNAME}/src/.libs/libhdf5.a"
HDF5_CXX_LIB="${BASEPBNAME}/c++/src/.libs/libhdf5_cpp.a"
SZIP_LIB="${BASEPBNAME}/szip/src/.libs/libsz.a"

mkdir -p "${USER_INCLUDE}"
cp ${HDF5_INCLUDE}/*.h ${USER_INCLUDE}
cp ${HDF5_CXX_INCLUDE}/*.h ${USER_INCLUDE}
mkdir -p "${USER_LIB_DIR}"
cp ${HDF5_LIB} ${USER_LIB_DIR}
cp ${HDF5_CXX_LIB} ${USER_LIB_DIR}
cp ${SZIP_LIB} ${USER_LIB_DIR}
