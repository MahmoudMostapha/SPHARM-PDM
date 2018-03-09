/*
 * compute the spherical parametrization of a binary mask
 *
 * author:  Martin Styner
 *
 */

#include <fstream>
#include <string>
#include <string.h>
//stringstream
#include <sstream>
// std::cout
#include <iostream>
// std::vector
#include <vector>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkMeshSpatialObject.h>
#include <itkMesh.h>
#include <itkSpatialObjectWriter.h>
#include <itkSpatialObjectReader.h>
#include <itkImageFileReader.h>
#include <itkDefaultDynamicMeshTraits.h>

#include "BinaryMask3DEqualAreaParametricMeshSource.h"

#include "GenParaMeshCLPCLP.h"
#include "itkMeshTovtkPolyData.h"
#include "vtkPolyDataToitkMesh.h"

#include <vtkPolyDataWriter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>

// Conformal Mapping
#include "itkVTKPolyDataReader.h"
#include "itkVTKPolyDataWriter.h"
#include "myConformalFlatteningMeshFilter.h"

typedef itk::Image<unsigned short, 3>                             ImageType;
typedef itk::ImageFileReader<ImageType>                           VolumeReaderType;
typedef itk::BinaryMask3DEqualAreaParametricMeshSource<ImageType> MeshSourceType;
typedef MeshSourceType::OutputMeshType                            OutputMeshType;
typedef itk::MeshSpatialObject<OutputMeshType>                    MeshSpatialObjectType;
typedef itk::ImageFileReader<ImageType> VolumeReaderType;

typedef itk::DefaultDynamicMeshTraits<double, 3, 3, double, double> MeshTraitsType;
typedef itk::Mesh<double, 3, MeshTraitsType>                        itkMeshType;
typedef itk::MeshSpatialObject<itkMeshType>                         itkMeshSOType;
typedef itk::MetaMeshConverter<3, double, MeshTraitsType>           MeshConverterType;

// ConformalMapping
// For future work, make "typedef itk::Mesh< double, 3, MeshTraitsType > MeshType;" work, need some revision
// on myConformalFlatteningFilter.txx

typedef itk::Mesh< double, 3 > MeshType;
typedef itk::ConformalFlatteningMeshFilter<MeshType, MeshType>  FilterType;
typedef MeshType::CellIdentifier  CellIdentifier;
typedef itk::VTKPolyDataReader<MeshType>  ITKReaderType;
typedef itk::VTKPolyDataWriter<MeshType>  ITKWriterType;

using namespace std;
void WriteEulerFile( std::string outEulerName, int Eulernum);
vtkSmartPointer<vtkPolyData> ReadPolyData(std::string filePath);
void WritePolyData(MeshType::Pointer meshInput, std::string& fileName);
std::string generatePrefix(std::string ConfOrInit, std::string infile, std::string initParaFileName, int numIterations);

// static int debug = 0;

int main( int argc, char * argv[] )
{

  PARSE_ARGS;
  // Read the input image
  ImageType::Pointer image;
  bool               initParaFile;
  if( initParaFileName == "NULL" )
    {
        initParaFile = false;
    }
  else
    {
        initParaFile = true;
    }
  if( debug )
    {
        std::cout << "Reading Image: " << std::endl;
    }
  VolumeReaderType::Pointer imageReader = VolumeReaderType::New();
  imageReader->SetFileName(infile);
  imageReader->Update();
  image = imageReader->GetOutput();
  // Create log file if log flag is on
  ofstream log;
  if( logFile )
    {
        log.open(outLogName.c_str(), ios::out | ios::app);
    }

  // if necessary read parmesh. When --conf is set, no initParaFile allowed(--conf set confParaMesh as default initPara); otherwise, initParaFile is allowed
  MeshSourceType::Pointer meshsrc = MeshSourceType::New();
  OutputMeshType::Pointer mesh;
  OutputMeshType::Pointer mesh0;
  OutputMeshType::Pointer parmesh;
  try
  {
        if( debug )
          {
          std::cout << "Creating Para Surface Mesh: " << std::endl;
          }
        if (useConformalMapping) {

            // Run for 0 iterations
            meshsrc->SetInput(image);
            meshsrc->SetNumberOfIterations(0);
            meshsrc->SetObjectValue(label);
            meshsrc->Update();
            // Output Mesh
            mesh0 = meshsrc->GetSurfaceMesh();

            // convert surfaces to VTK
            itkMeshTovtkPolyData ITKVTKConverter;
            ITKVTKConverter.SetInput(mesh0);
            vtkSmartPointer <vtkPolyData> SurfMesh0 = vtkSmartPointer<vtkPolyData>::New();
            SurfMesh0 = ITKVTKConverter.GetOutput();
            if( debug )
            {
                std::cout << "SurfMesh0 has " << SurfMesh0->GetNumberOfPoints() << " points." << std::endl;
            }

            vtkSmartPointer <vtkPolyDataWriter> SurfMesh0writer = vtkSmartPointer<vtkPolyDataWriter>::New();

            // Writing the SurfMesh0 file
            #if VTK_MAJOR_VERSION > 5
            SurfMesh0writer->SetInputData(SurfMesh0);
            #else
            SurfMesh0writer->SetInput(SurfMesh0);
            #endif
            std::string Surf0Common = generatePrefix("ConfInitPara", infile, "", numIterations);

            //write the SurfMesh0Name filename as Surf0Common-surf0.vtk
            std::string SurfMesh0Name("dummy");
            SurfMesh0Name.erase();
            SurfMesh0Name.append(Surf0Common);
            SurfMesh0Name.append("-surf0.vtk");
            SurfMesh0writer->SetFileName(SurfMesh0Name.c_str());
            SurfMesh0writer->Write();

            // Read the SurfMesh0 file
            ITKReaderType::Pointer ITKreader = ITKReaderType::New();
            ITKreader->SetFileName(SurfMesh0Name.c_str());
            try {
                ITKreader->Update();
            }
            catch (itk::ExceptionObject &excp) {
                std::cerr << excp << std::endl;
                return EXIT_FAILURE;
            }
            MeshType::Pointer SurfMesh0_R = ITKreader->GetOutput();


            //try to make this work in the future. MeshType can't be initialized to have the same mesh
            //traits as the BinaryMask3DEqualAreaParametricMeshSource because conformal filter
            //does not work with that mesh traits
            //MeshType::Pointer SurfMesh0_R = mesh0;

            // Use the output mesh to run ConformalMapping
            FilterType::Pointer filter = FilterType::New();
            filter->SetInput(SurfMesh0_R);
            CellIdentifier polarCellId = 0; // default set to the first cell
            filter->SetPolarCellIdentifier(polarCellId);
            filter->MapToSphere();
            // Execute the filter
            if( debug )
            {
                std::cout << "Execute the filter" << std::endl;
            }

            try {
                filter->Update();
            }
            catch (itk::ExceptionObject &excp) {
                std::cerr << excp << std::endl;
                return EXIT_FAILURE;
            }
            // Get the Smart Pointer to the Filter Output
            MeshType::Pointer Confpara = filter->GetOutput();
            // Write as VTK
            std::string ConfPara0Name("dummy");
            ConfPara0Name.erase();
            ConfPara0Name.append(Surf0Common);
            ConfPara0Name.append("-confPara0.vtk");
            //write VTK
            WritePolyData(Confpara, ConfPara0Name);
            // Read VTK
            vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(ConfPara0Name);
            //make sure it reads the VTK correctly
            if( debug )
            {
                std::cout << ConfPara0Name << " has " << polydata->GetNumberOfPoints() << " points." << std::endl;
            }
            // convert to itk mesh data structure
            vtkPolyDataToitkMesh vtkItkConverter;
            vtkItkConverter.SetInput( polydata );
            if (debug) {
                std::cout << "converting mesh " << ConfPara0Name <<std::endl;
            }
            // write out the itk meta mesh file
            itkMeshSOType::Pointer meshSO = itkMeshSOType::New();
            meshSO->SetMesh( vtkItkConverter.GetOutput() );
            MeshConverterType::Pointer itkConverter = MeshConverterType::New();
            std::string ConfPara0NameMeta("dummy");
            ConfPara0NameMeta.erase();
            ConfPara0NameMeta.append(Surf0Common);
            ConfPara0NameMeta.append("-confPara0.meta");
            itkConverter->WriteMeta( meshSO, ConfPara0NameMeta.c_str() );

            typedef itk::SpatialObjectReader<3, double, OutputMeshType::MeshTraits> ReaderType;
            ReaderType::Pointer readerSH = ReaderType::New();
            try
            {
                if( debug )
                {
                    std::cout << "Reading  " << ConfPara0NameMeta << std::endl;
                }
                readerSH->SetFileName(ConfPara0NameMeta);
                readerSH->Update();
                ReaderType::SceneType::Pointer          scene1 = readerSH->GetScene();
                ReaderType::SceneType::ObjectListType * objList =  scene1->GetObjects(1, NULL);
                // TODO: plugin name if multiple object are present
                ReaderType::SceneType::ObjectListType::iterator it = objList->begin();
                itk::SpatialObject<3> *                         curObj = *it;
                MeshSpatialObjectType::Pointer                  paraSOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
                parmesh = paraSOMesh->GetMesh();
            }
            catch( itk::ExceptionObject ex )
            {
                std::cout << ex.GetDescription() << std::endl;
                return EXIT_FAILURE ;
            }
            // Use ConformalMapping output as input
            meshsrc->SetInput(image);
            meshsrc->SetNumberOfIterations(numIterations);
            meshsrc->SetObjectValue(label);
            meshsrc->SetInitParametricMesh(parmesh);
            meshsrc->Update();
            // Output Mesh
            mesh = meshsrc->GetSurfaceMesh();
            // Create the mesh Spatial Object
        }
            //not using conformal mapping
        else
        {
            // read initParaFile if set
            vtkPolyDataToitkMesh vtkItkConverter;
            typedef itk::SpatialObjectReader<3, double, OutputMeshType::MeshTraits> ReaderType;
            ReaderType::Pointer readerSH = ReaderType::New();

            if( initParaFile )
            {
                if (strstr(initParaFileName.c_str(),".meta"))
                {
                    try
                    {
                        if( debug )
                        {
                            std::cout << "Reading " << initParaFileName << std::endl;
                        }
                        readerSH->SetFileName(initParaFileName);
                        readerSH->Update();
                        ReaderType::SceneType::Pointer          scene1 = readerSH->GetScene();
                        ReaderType::SceneType::ObjectListType * objList =  scene1->GetObjects(1, NULL);
                        // TODO: plugin name if multiple object are present
                        ReaderType::SceneType::ObjectListType::iterator it = objList->begin();
                        itk::SpatialObject<3> *                         curObj = *it;
                        MeshSpatialObjectType::Pointer                  paraSOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
                        parmesh = paraSOMesh->GetMesh();
                    }
                    catch( itk::ExceptionObject ex )
                    {
                        std::cout << ex.GetDescription() << std::endl;
                        return EXIT_FAILURE ;
                    }
                }

                if (strstr(initParaFileName.c_str(),".vtk"))
                {
                    // convert surfaces to meta
                    try
                    {
                        if( debug )
                        {
                            std::cout << "Reading vtk " << initParaFileName << std::endl;
                        }

                        // Read VTK
                        vtkSmartPointer<vtkPolyData> polydata = ReadPolyData(initParaFileName);
                        //make sure it reads the VTK correctly
                        if( debug )
                        {
                            std::cout << initParaFileName << " has " << polydata->GetNumberOfPoints() << " points." << std::endl;
                        }
                        // convert to itk mesh data structure
                        vtkItkConverter.SetInput( polydata );

                        if (debug) {
                            std::cout << "converting mesh " << initParaFileName <<std::endl;
                        }

                        // write out the itk meta mesh file
                        itkMeshSOType::Pointer meshSO = itkMeshSOType::New();
                        meshSO->SetMesh( vtkItkConverter.GetOutput() );
                        MeshConverterType::Pointer itkConverter = MeshConverterType::New();

                        std::string ConfPara0NameMetaPrefix = generatePrefix("Init", infile, initParaFileName, numIterations);

                        std::string ConfPara0NameMeta("dummy");
                        ConfPara0NameMeta.erase();
                        ConfPara0NameMeta.append(ConfPara0NameMetaPrefix);
                        ConfPara0NameMeta.append("-confPara0.meta");
                        itkConverter->WriteMeta( meshSO, ConfPara0NameMeta.c_str() );


                        typedef itk::SpatialObjectReader<3, double, OutputMeshType::MeshTraits> ReaderType;
                        try
                        {
                            if( debug )
                            {
                                std::cout << "Reading InitPara.meta" << std::endl;
                            }
                            readerSH->SetFileName("InitPara.meta");
                            readerSH->Update();
                            ReaderType::SceneType::Pointer          scene1 = readerSH->GetScene();
                            ReaderType::SceneType::ObjectListType * objList =  scene1->GetObjects(1, NULL);
                            // TODO: plugin name if multiple object are present
                            ReaderType::SceneType::ObjectListType::iterator it = objList->begin();
                            itk::SpatialObject<3> *                         curObj = *it;
                            MeshSpatialObjectType::Pointer                  paraSOMesh = dynamic_cast<MeshSpatialObjectType *>(curObj);
                            parmesh = paraSOMesh->GetMesh();
                        }
                        catch( itk::ExceptionObject ex )
                        {
                            std::cout << ex.GetDescription() << std::endl;
                            return EXIT_FAILURE ;
                        }
                        if( debug )
                        {
                            std::cout << "Converting vtk done"  << std::endl;
                        }
                    }
                    catch( itk::ExceptionObject ex )
                    {
                        std::cout << ex.GetDescription() << std::endl;
                        return EXIT_FAILURE ;
                    }
                }
            }

            meshsrc->SetInput(image);
            meshsrc->SetNumberOfIterations(numIterations);
            meshsrc->SetObjectValue(label);
            if( initParaFile )
            {
              meshsrc->SetInitParametricMesh(parmesh);
            }
            meshsrc->Update();
            // Output Mesh
            mesh = meshsrc->GetSurfaceMesh();
            // Create the mesh Spatial Object
        }

        if (meshsrc->GetEulerNum() == 2 )
        {
            //duplicate writer, could make this global and use it in if/else
            vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();

            {
                // convert surfaces to VTK
                itkMeshTovtkPolyData ITKVTKConverter;
                ITKVTKConverter.SetInput(mesh);
                vtkSmartPointer<vtkPolyData> SurfMesh = vtkSmartPointer<vtkPolyData>::New();
                SurfMesh = ITKVTKConverter.GetOutput();

                // Writing the file
                #if VTK_MAJOR_VERSION > 5
                writer->SetInputData(SurfMesh);
                #else
                writer->SetInput(SurfMesh);
                #endif
                writer->SetFileName(outSurfName.c_str() );
                writer->Write();
            }

            if( EulerFile )
            {
                WriteEulerFile(outEulerName, meshsrc->GetEulerNum() );
            }

            // Output Mesh
            parmesh = meshsrc->GetParametrizationMesh();
            // Create the mesh Spatial Object
            {
                // convert surfaces to VTK
                itkMeshTovtkPolyData ITKVTKConverter2;
                ITKVTKConverter2.SetInput( parmesh );
                vtkSmartPointer<vtkPolyData> ParaMesh = vtkSmartPointer<vtkPolyData>::New();
                ParaMesh = ITKVTKConverter2.GetOutput();
                // delete (ITKVTKConverter2);
                // Writing the file
                #if VTK_MAJOR_VERSION > 5
                writer->SetInputData(ParaMesh);
                #else
                writer->SetInput(ParaMesh);
                #endif
                writer->SetFileName(outParaName.c_str() );
                writer->Write();
                if( logFile )
                {
                    log << "Computed " << infile << std::endl;
                }
            }
       }
  }

  catch( itk::ExceptionObject e )
  {
    if( EulerFile )
    {
        WriteEulerFile(outEulerName, meshsrc->GetEulerNum() );
    }
    e.Print(std::cout);
    if( logFile )
    {
        log << "Failed " << infile << " " << e.what();
    }
    return EXIT_FAILURE ;
  }

  return EXIT_SUCCESS ;
}

void WriteEulerFile( std::string outEulerName, int Eulernum)
{
  std::ofstream file( outEulerName.c_str() );

  if( file )
  {
      file << Eulernum << std::endl;
      file.close();
  }

}

vtkSmartPointer<vtkPolyData> ReadPolyData(std::string filePath)
{
    size_t found = filePath.rfind(".vtp");
    if (found != std::string::npos)
    {
        vtkSmartPointer<vtkXMLPolyDataReader> VTPreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        VTPreader->SetFileName(filePath.c_str());
        VTPreader->Update();
        return VTPreader->GetOutput();
    }
    found = filePath.rfind(".vtk");
    if ( found != std::string::npos )
    {
        vtkSmartPointer<vtkPolyDataReader> VTKreader = vtkSmartPointer<vtkPolyDataReader>::New();
        VTKreader->SetFileName(filePath.c_str());
        VTKreader->Update();
        return VTKreader->GetOutput();
    }
    return NULL;
}

void WritePolyData(MeshType::Pointer meshInput, std::string &fileName)
{
    ITKWriterType::Pointer ITKwriter = ITKWriterType::New();
    ITKwriter->SetInput(meshInput);
    ITKwriter->SetFileName(fileName);
    ITKwriter->Update();
}
//generate prefix of filename as "label(string before first dot of the original filename)-iter[numIterations]"
std::string generatePrefix(std::string ConfOrInit, std::string infile, std::string initParaFileName, int numIterations)
{
    std::string prefix("dummy");
    prefix.erase();
    std::istringstream iss;
    if (ConfOrInit.compare("ConfInit"))
    {
        iss.str(infile);
    }
    else if (ConfOrInit.compare("Init"))
    {
        iss.str(initParaFileName);
    }
    else
    {
        return NULL;
    }
    std::vector<std::string> tokens;
    std::string token;
    while (std::getline(iss, token, '.'))
    {
        if (!token.empty())
            tokens.push_back(token);
    }
    prefix.append(tokens[0]);
    prefix.append("-iter");
    std::stringstream ni;
    ni << numIterations;
    string iter = ni.str();
    prefix.append(iter);
//    std::cout<<"prefix generated is: " << prefix <<endl;
    return prefix;
}