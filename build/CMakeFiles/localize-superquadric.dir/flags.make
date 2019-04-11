# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.5

# compile CXX with /usr/lib/ccache/c++
CXX_FLAGS =   -O3 -DNDEBUG   -std=c++11

CXX_DEFINES = -DHAVE_CSTDDEF -D_USE_MATH_DEFINES -DvtkDomainsChemistry_AUTOINIT="1(vtkDomainsChemistryOpenGL2)" -DvtkIOExport_AUTOINIT="1(vtkIOExportOpenGL2)" -DvtkRenderingContext2D_AUTOINIT="1(vtkRenderingContextOpenGL2)" -DvtkRenderingCore_AUTOINIT="3(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingOpenGL2)" -DvtkRenderingOpenGL2_AUTOINIT="1(vtkRenderingGL2PSOpenGL2)" -DvtkRenderingVolume_AUTOINIT="1(vtkRenderingVolumeOpenGL2)"

CXX_INCLUDES = -I/home/gvezzani/dependencies/VTK/build/Utilities/KWIML -I/home/gvezzani/dependencies/VTK/Utilities/KWIML -I/home/gvezzani/dependencies/VTK/build/Utilities/KWSys -I/home/gvezzani/dependencies/VTK/Utilities/KWSys -I/home/gvezzani/dependencies/VTK/build/Common/Core -I/home/gvezzani/dependencies/VTK/Common/Core -I/home/gvezzani/dependencies/VTK/build/Common/Math -I/home/gvezzani/dependencies/VTK/Common/Math -I/home/gvezzani/dependencies/VTK/build/Common/Misc -I/home/gvezzani/dependencies/VTK/Common/Misc -I/home/gvezzani/dependencies/VTK/build/Common/System -I/home/gvezzani/dependencies/VTK/Common/System -I/home/gvezzani/dependencies/VTK/build/Common/Transforms -I/home/gvezzani/dependencies/VTK/Common/Transforms -I/home/gvezzani/dependencies/VTK/build/Common/DataModel -I/home/gvezzani/dependencies/VTK/Common/DataModel -I/home/gvezzani/dependencies/VTK/build/Common/Color -I/home/gvezzani/dependencies/VTK/Common/Color -I/home/gvezzani/dependencies/VTK/build/Common/ExecutionModel -I/home/gvezzani/dependencies/VTK/Common/ExecutionModel -I/home/gvezzani/dependencies/VTK/build/Common/ComputationalGeometry -I/home/gvezzani/dependencies/VTK/Common/ComputationalGeometry -I/home/gvezzani/dependencies/VTK/build/Filters/Core -I/home/gvezzani/dependencies/VTK/Filters/Core -I/home/gvezzani/dependencies/VTK/build/Filters/General -I/home/gvezzani/dependencies/VTK/Filters/General -I/home/gvezzani/dependencies/VTK/build/Imaging/Core -I/home/gvezzani/dependencies/VTK/Imaging/Core -I/home/gvezzani/dependencies/VTK/build/Imaging/Fourier -I/home/gvezzani/dependencies/VTK/Imaging/Fourier -I/home/gvezzani/dependencies/VTK/build/ThirdParty/alglib -I/home/gvezzani/dependencies/VTK/ThirdParty/alglib -I/home/gvezzani/dependencies/VTK/build/Filters/Statistics -I/home/gvezzani/dependencies/VTK/Filters/Statistics -I/home/gvezzani/dependencies/VTK/build/Filters/Extraction -I/home/gvezzani/dependencies/VTK/Filters/Extraction -I/home/gvezzani/dependencies/VTK/build/Infovis/Core -I/home/gvezzani/dependencies/VTK/Infovis/Core -I/home/gvezzani/dependencies/VTK/build/Filters/Geometry -I/home/gvezzani/dependencies/VTK/Filters/Geometry -I/home/gvezzani/dependencies/VTK/build/Filters/Sources -I/home/gvezzani/dependencies/VTK/Filters/Sources -I/home/gvezzani/dependencies/VTK/build/Rendering/Core -I/home/gvezzani/dependencies/VTK/Rendering/Core -I/home/gvezzani/dependencies/VTK/build/ThirdParty/zlib -I/home/gvezzani/dependencies/VTK/ThirdParty/zlib -I/home/gvezzani/dependencies/VTK/build/ThirdParty/freetype -I/home/gvezzani/dependencies/VTK/ThirdParty/freetype -I/home/gvezzani/dependencies/VTK/build/Rendering/FreeType -I/home/gvezzani/dependencies/VTK/Rendering/FreeType -I/home/gvezzani/dependencies/VTK/build/Rendering/Context2D -I/home/gvezzani/dependencies/VTK/Rendering/Context2D -I/home/gvezzani/dependencies/VTK/build/Charts/Core -I/home/gvezzani/dependencies/VTK/Charts/Core -I/home/gvezzani/dependencies/VTK/ThirdParty/lz4/vtklz4/lib -I/home/gvezzani/dependencies/VTK/build/ThirdParty/lz4/vtklz4 -I/home/gvezzani/dependencies/VTK/build/ThirdParty/lz4 -I/home/gvezzani/dependencies/VTK/ThirdParty/lz4 -I/home/gvezzani/dependencies/VTK/build/IO/Core -I/home/gvezzani/dependencies/VTK/IO/Core -I/home/gvezzani/dependencies/VTK/build/IO/Legacy -I/home/gvezzani/dependencies/VTK/IO/Legacy -I/home/gvezzani/dependencies/VTK/build/ThirdParty/expat -I/home/gvezzani/dependencies/VTK/ThirdParty/expat -I/home/gvezzani/dependencies/VTK/build/IO/XMLParser -I/home/gvezzani/dependencies/VTK/IO/XMLParser -I/home/gvezzani/dependencies/VTK/build/IO/XML -I/home/gvezzani/dependencies/VTK/IO/XML -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libxml2/vtklibxml2 -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libxml2 -I/home/gvezzani/dependencies/VTK/ThirdParty/libxml2 -I/home/gvezzani/dependencies/VTK/build/IO/Infovis -I/home/gvezzani/dependencies/VTK/IO/Infovis -I/home/gvezzani/dependencies/VTK/build/Utilities/EncodeString -I/home/gvezzani/dependencies/VTK/Utilities/EncodeString -I/home/gvezzani/dependencies/VTK/build/ThirdParty/glew -I/home/gvezzani/dependencies/VTK/ThirdParty/glew -I/home/gvezzani/dependencies/VTK/build/Rendering/OpenGL2 -I/home/gvezzani/dependencies/VTK/Rendering/OpenGL2 -I/home/gvezzani/dependencies/VTK/build/Rendering/ContextOpenGL2 -I/home/gvezzani/dependencies/VTK/Rendering/ContextOpenGL2 -I/home/gvezzani/dependencies/VTK/build/Testing/Core -I/home/gvezzani/dependencies/VTK/Testing/Core -I/home/gvezzani/dependencies/VTK/build/Utilities/DICOMParser -I/home/gvezzani/dependencies/VTK/Utilities/DICOMParser -I/home/gvezzani/dependencies/VTK/build/Utilities/MetaIO/vtkmetaio -I/home/gvezzani/dependencies/VTK/build/Utilities/MetaIO -I/home/gvezzani/dependencies/VTK/Utilities/MetaIO -I/home/gvezzani/dependencies/VTK/build/ThirdParty/jpeg -I/home/gvezzani/dependencies/VTK/ThirdParty/jpeg -I/home/gvezzani/dependencies/VTK/build/ThirdParty/png -I/home/gvezzani/dependencies/VTK/ThirdParty/png -I/home/gvezzani/dependencies/VTK/build/ThirdParty/tiff/vtktiff/libtiff -I/home/gvezzani/dependencies/VTK/build/ThirdParty/tiff -I/home/gvezzani/dependencies/VTK/ThirdParty/tiff -I/home/gvezzani/dependencies/VTK/build/IO/Image -I/home/gvezzani/dependencies/VTK/IO/Image -I/home/gvezzani/dependencies/VTK/build/Testing/Rendering -I/home/gvezzani/dependencies/VTK/Testing/Rendering -I/home/gvezzani/dependencies/VTK/build/Imaging/Sources -I/home/gvezzani/dependencies/VTK/Imaging/Sources -I/home/gvezzani/dependencies/VTK/build/Filters/Hybrid -I/home/gvezzani/dependencies/VTK/Filters/Hybrid -I/home/gvezzani/dependencies/VTK/build/Filters/Modeling -I/home/gvezzani/dependencies/VTK/Filters/Modeling -I/home/gvezzani/dependencies/VTK/build/Imaging/Color -I/home/gvezzani/dependencies/VTK/Imaging/Color -I/home/gvezzani/dependencies/VTK/build/Imaging/General -I/home/gvezzani/dependencies/VTK/Imaging/General -I/home/gvezzani/dependencies/VTK/build/Imaging/Hybrid -I/home/gvezzani/dependencies/VTK/Imaging/Hybrid -I/home/gvezzani/dependencies/VTK/build/Interaction/Style -I/home/gvezzani/dependencies/VTK/Interaction/Style -I/home/gvezzani/dependencies/VTK/build/Rendering/Annotation -I/home/gvezzani/dependencies/VTK/Rendering/Annotation -I/home/gvezzani/dependencies/VTK/build/Rendering/Volume -I/home/gvezzani/dependencies/VTK/Rendering/Volume -I/home/gvezzani/dependencies/VTK/build/Interaction/Widgets -I/home/gvezzani/dependencies/VTK/Interaction/Widgets -I/home/gvezzani/dependencies/VTK/build/Views/Core -I/home/gvezzani/dependencies/VTK/Views/Core -I/home/gvezzani/dependencies/VTK/build/Views/Context2D -I/home/gvezzani/dependencies/VTK/Views/Context2D -I/home/gvezzani/dependencies/VTK/build/Filters/Generic -I/home/gvezzani/dependencies/VTK/Filters/Generic -I/home/gvezzani/dependencies/VTK/build/IO/Geometry -I/home/gvezzani/dependencies/VTK/IO/Geometry -I/home/gvezzani/dependencies/VTK/build/Testing/GenericBridge -I/home/gvezzani/dependencies/VTK/Testing/GenericBridge -I/home/gvezzani/dependencies/VTK/build/Domains/Chemistry -I/home/gvezzani/dependencies/VTK/Domains/Chemistry -I/home/gvezzani/dependencies/VTK/build/Domains/ChemistryOpenGL2 -I/home/gvezzani/dependencies/VTK/Domains/ChemistryOpenGL2 -I/home/gvezzani/dependencies/VTK/build/Utilities/HashSource -I/home/gvezzani/dependencies/VTK/Utilities/HashSource -I/home/gvezzani/dependencies/VTK/build/Parallel/Core -I/home/gvezzani/dependencies/VTK/Parallel/Core -I/home/gvezzani/dependencies/VTK/build/Filters/AMR -I/home/gvezzani/dependencies/VTK/Filters/AMR -I/home/gvezzani/dependencies/VTK/build/ThirdParty/hdf5/vtkhdf5 -isystem /home/gvezzani/dependencies/VTK/ThirdParty/hdf5/vtkhdf5/hl/src -isystem /home/gvezzani/dependencies/VTK/ThirdParty/hdf5/vtkhdf5/src -I/home/gvezzani/dependencies/VTK/build/ThirdParty/hdf5 -I/home/gvezzani/dependencies/VTK/ThirdParty/hdf5 -I/home/gvezzani/dependencies/VTK/build/IO/AMR -I/home/gvezzani/dependencies/VTK/IO/AMR -I/home/gvezzani/dependencies/VTK/ThirdParty/netcdf/vtknetcdf/include -I/home/gvezzani/dependencies/VTK/build/ThirdParty/netcdf/vtknetcdf -I/home/gvezzani/dependencies/VTK/build/ThirdParty/netcdf -I/home/gvezzani/dependencies/VTK/ThirdParty/netcdf -I/home/gvezzani/dependencies/VTK/build/ThirdParty/exodusII -I/home/gvezzani/dependencies/VTK/ThirdParty/exodusII -I/home/gvezzani/dependencies/VTK/build/IO/Exodus -I/home/gvezzani/dependencies/VTK/IO/Exodus -I/home/gvezzani/dependencies/VTK/build/Imaging/Math -I/home/gvezzani/dependencies/VTK/Imaging/Math -I/home/gvezzani/dependencies/VTK/build/Rendering/VolumeOpenGL2 -I/home/gvezzani/dependencies/VTK/Rendering/VolumeOpenGL2 -I/home/gvezzani/dependencies/VTK/build/Filters/FlowPaths -I/home/gvezzani/dependencies/VTK/Filters/FlowPaths -I/home/gvezzani/dependencies/VTK/build/Filters/Imaging -I/home/gvezzani/dependencies/VTK/Filters/Imaging -I/home/gvezzani/dependencies/VTK/build/Rendering/Label -I/home/gvezzani/dependencies/VTK/Rendering/Label -I/home/gvezzani/dependencies/VTK/build/Filters/HyperTree -I/home/gvezzani/dependencies/VTK/Filters/HyperTree -I/home/gvezzani/dependencies/VTK/build/Imaging/Stencil -I/home/gvezzani/dependencies/VTK/Imaging/Stencil -I/home/gvezzani/dependencies/VTK/build/Filters/Parallel -I/home/gvezzani/dependencies/VTK/Filters/Parallel -I/home/gvezzani/dependencies/VTK/build/Filters/ParallelImaging -I/home/gvezzani/dependencies/VTK/Filters/ParallelImaging -I/home/gvezzani/dependencies/VTK/build/Filters/Points -I/home/gvezzani/dependencies/VTK/Filters/Points -I/home/gvezzani/dependencies/VTK/build/Filters/Programmable -I/home/gvezzani/dependencies/VTK/Filters/Programmable -I/home/gvezzani/dependencies/VTK/build/Filters/SMP -I/home/gvezzani/dependencies/VTK/Filters/SMP -I/home/gvezzani/dependencies/VTK/build/Filters/Selection -I/home/gvezzani/dependencies/VTK/Filters/Selection -I/home/gvezzani/dependencies/VTK/build/ThirdParty/verdict -I/home/gvezzani/dependencies/VTK/ThirdParty/verdict -I/home/gvezzani/dependencies/VTK/build/Filters/Verdict -I/home/gvezzani/dependencies/VTK/Filters/Verdict -I/home/gvezzani/dependencies/VTK/build/ThirdParty/netcdfcpp -I/home/gvezzani/dependencies/VTK/ThirdParty/netcdfcpp -I/home/gvezzani/dependencies/VTK/build/IO/NetCDF -I/home/gvezzani/dependencies/VTK/IO/NetCDF -I/home/gvezzani/dependencies/VTK/build/ThirdParty/jsoncpp -I/home/gvezzani/dependencies/VTK/ThirdParty/jsoncpp -I/home/gvezzani/dependencies/VTK/build/IO/Parallel -I/home/gvezzani/dependencies/VTK/IO/Parallel -I/home/gvezzani/dependencies/VTK/build/Filters/Texture -I/home/gvezzani/dependencies/VTK/Filters/Texture -I/home/gvezzani/dependencies/VTK/build/Filters/Topology -I/home/gvezzani/dependencies/VTK/Filters/Topology -I/home/gvezzani/dependencies/VTK/build/Infovis/Layout -I/home/gvezzani/dependencies/VTK/Infovis/Layout -I/home/gvezzani/dependencies/VTK/ThirdParty/libproj4/vtklibproj4 -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libproj4/vtklibproj4 -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libproj4 -I/home/gvezzani/dependencies/VTK/ThirdParty/libproj4 -I/home/gvezzani/dependencies/VTK/build/Geovis/Core -I/home/gvezzani/dependencies/VTK/Geovis/Core -I/home/gvezzani/dependencies/VTK/build/IO/EnSight -I/home/gvezzani/dependencies/VTK/IO/EnSight -I/home/gvezzani/dependencies/VTK/build/ThirdParty/gl2ps -I/home/gvezzani/dependencies/VTK/ThirdParty/gl2ps -I/home/gvezzani/dependencies/VTK/build/Rendering/GL2PSOpenGL2 -I/home/gvezzani/dependencies/VTK/Rendering/GL2PSOpenGL2 -I/home/gvezzani/dependencies/VTK/ThirdParty/libharu/vtklibharu/include -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libharu/vtklibharu/include -I/home/gvezzani/dependencies/VTK/build/ThirdParty/libharu -I/home/gvezzani/dependencies/VTK/ThirdParty/libharu -I/home/gvezzani/dependencies/VTK/build/IO/Export -I/home/gvezzani/dependencies/VTK/IO/Export -I/home/gvezzani/dependencies/VTK/build/IO/ExportOpenGL2 -I/home/gvezzani/dependencies/VTK/IO/ExportOpenGL2 -I/home/gvezzani/dependencies/VTK/build/Interaction/Image -I/home/gvezzani/dependencies/VTK/Interaction/Image -I/home/gvezzani/dependencies/VTK/build/IO/Import -I/home/gvezzani/dependencies/VTK/IO/Import -I/home/gvezzani/dependencies/VTK/build/IO/LSDyna -I/home/gvezzani/dependencies/VTK/IO/LSDyna -I/home/gvezzani/dependencies/VTK/build/IO/MINC -I/home/gvezzani/dependencies/VTK/IO/MINC -I/home/gvezzani/dependencies/VTK/build/ThirdParty/oggtheora -I/home/gvezzani/dependencies/VTK/ThirdParty/oggtheora -I/home/gvezzani/dependencies/VTK/build/IO/Movie -I/home/gvezzani/dependencies/VTK/IO/Movie -I/home/gvezzani/dependencies/VTK/build/IO/PLY -I/home/gvezzani/dependencies/VTK/IO/PLY -I/home/gvezzani/dependencies/VTK/build/IO/ParallelXML -I/home/gvezzani/dependencies/VTK/IO/ParallelXML -I/home/gvezzani/dependencies/VTK/build/ThirdParty/sqlite -I/home/gvezzani/dependencies/VTK/ThirdParty/sqlite -I/home/gvezzani/dependencies/VTK/build/IO/SQL -I/home/gvezzani/dependencies/VTK/IO/SQL -I/home/gvezzani/dependencies/VTK/build/Testing/IOSQL -I/home/gvezzani/dependencies/VTK/Testing/IOSQL -I/home/gvezzani/dependencies/VTK/build/IO/TecplotTable -I/home/gvezzani/dependencies/VTK/IO/TecplotTable -I/home/gvezzani/dependencies/VTK/build/IO/Video -I/home/gvezzani/dependencies/VTK/IO/Video -I/home/gvezzani/dependencies/VTK/build/Imaging/Statistics -I/home/gvezzani/dependencies/VTK/Imaging/Statistics -I/home/gvezzani/dependencies/VTK/build/Rendering/Image -I/home/gvezzani/dependencies/VTK/Rendering/Image -I/home/gvezzani/dependencies/VTK/build/Imaging/Morphological -I/home/gvezzani/dependencies/VTK/Imaging/Morphological -I/home/gvezzani/dependencies/VTK/build/Rendering/LOD -I/home/gvezzani/dependencies/VTK/Rendering/LOD -I/home/gvezzani/dependencies/VTK/build/Views/Infovis -I/home/gvezzani/dependencies/VTK/Views/Infovis -I/home/gvezzani/robot-code/localize-superquadric/src -I/usr/include/coin -isystem /home/gvezzani/robot-install/include -isystem /home/gvezzani/dependencies/VTK/build/Utilities/KWSys/vtksys -isystem /home/gvezzani/dependencies/VTK/build/ThirdParty/hdf5/vtkhdf5/hl/src -isystem /home/gvezzani/dependencies/VTK/build/ThirdParty/hdf5/vtkhdf5/src -isystem /home/gvezzani/dependencies/VTK/build/ThirdParty/netcdfcpp/vtknetcdfcpp 

