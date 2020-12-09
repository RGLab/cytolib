#generate pb code
protoc --cpp_out=cytolib/ GatingSet.proto
#move the source to src
mv cytolib/GatingSet.pb.cc ../../src/
#update the include path to header
sed -i 's+"GatingSet.pb.h"+<cytolib/GatingSet.pb.h>+g' ../../src/GatingSet.pb.cc