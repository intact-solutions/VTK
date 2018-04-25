#include <vtkPolyData.h>
#include <vtkPLYReader.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkDataSetMapper.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkLookupTable.h>
#include "json.h"

int main(int argc, char *argv[])
{
  ifstream inFile;
  inFile.open("bridge_samples.json");

  if (!inFile) {
    cout << "Unable to open file";
    exit(1); // terminate with error
  }

  Json::Value root;   // 'root' will contain the root value after parsing.
  inFile >> root;
  const Json::Value stress = root["123"]["von_mises_stress"];

  inFile.close();

  // Iterate over the sequence elements.
  std::vector<double> stress_vector;
  for (int index = 0; index < stress.size(); ++index) {
    stress_vector.push_back(stress[index].asDouble());
  }

  if (argc != 2)
  {
    std::cout << "Usage: " << argv[0] << "  Filename(.ply)" << std::endl;
    return EXIT_FAILURE;
  }

  std::string inputFilename = argv[1];

  vtkSmartPointer<vtkPLYReader> reader =
    vtkSmartPointer<vtkPLYReader>::New();
  reader->SetFileName(inputFilename.c_str());
  reader->Update();

  //add color
  vtkPolyData* outputPolyData = reader->GetOutput();
  //vtkPointData* point = outputPolyData->GetPointData()->SetAttribute;
  double range[2];
  // Generate the colors for each point based on the color map
  vtkSmartPointer<vtkDoubleArray> vtkstress =
    vtkSmartPointer<vtkDoubleArray>::New();

  vtkstress->SetNumberOfComponents(1);
  vtkstress->SetName("vtkstress");

  //for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
  for (int i = 0; i < stress_vector.size(); i++)
  {
    vtkstress->InsertNextValue(stress_vector[i]);
  }
  vtkstress->GetRange(range);
  outputPolyData->GetPointData()->SetAttribute(vtkstress, vtkDataSetAttributes::SCALARS);

  //create the colormap
  // Create the color map
  vtkSmartPointer<vtkLookupTable> colorLookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  colorLookupTable->SetTableRange(range);
  colorLookupTable->SetHueRange(.667, 0);
  colorLookupTable->Build();

 // create the height representation
  vtkNew<vtkDataSetMapper> mapper;
  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SelectColorArray("vtkstress");
  //mapper->SetScalarRange(range);
  mapper->SetLookupTable(colorLookupTable);

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  //renderer->SetBackground(0.1804, 0.5451, 0.3412); // Sea green

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}