#untar the lib
BASEPBNAME="protobuf-3.7.0"
PBTGZNAME=${BASEPBNAME}.tar.gz
if test -d ./${BASEPBNAME}; then 
	echo 'found ' $BASEPBNAME ' header sources and tar archive;using what is there.'
else
	echo "untarring protobuf ...";
	tar xf ${PBTGZNAME}
fi;

echo "building protobuf...";
PBBUILD=$(pwd)/pb_build
if test -d ${PBBUILD}; then
	echo 'found ' $PBBUILD ' ;using what is there.'
else
	mkdir ${PBBUILD}
fi;

cd ${BASEPBNAME};
./configure --enable-static=yes --with-pic=yes --enable-shared=no --prefix="${PBBUILD}" --libdir="${PBBUILD}/lib" 
make install 
