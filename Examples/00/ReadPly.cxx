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
#include <vtkCellDataToPointData.h>
#include <vtkPointDataToCellData.h>
#include <vtkExtractEdges.h>
#include <vtkFeatureEdges.h>
#include <vtkPolyDataWriter.h>
#include "json.h"

#define MAXDIST 10
vtkSmartPointer<vtkPolyData> smooth(vtkSmartPointer<vtkPolyData> polydata);
void CreateColor(vtkSmartPointer<vtkPolyData> outputPolyData);

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
  vtkSmartPointer<vtkPolyData> outputPolyData = reader->GetOutput();
  outputPolyData->BuildCells();
  outputPolyData->BuildLinks();
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

  //point to cell and then cell to point
  // cell 2 point and contour
  //vtkSmartPointer<vtkPointDataToCellData>p2c = vtkPointDataToCellData::New();
  //vtkSmartPointer<vtkCellDataToPointData> c2p = vtkCellDataToPointData::New();
  //p2c->SetInputConnection(reader->GetOutputPort());
  //p2c->PassPointDataOn();
  //c2p->SetInputConnection(p2c->GetOutputPort());
  //c2p->SetContributingCellOption(c2p->All);
  //c2p->PassCellDataOn();

  //create the colormap
  // Create the color map
  vtkSmartPointer<vtkLookupTable> colorLookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  colorLookupTable->SetTableRange(range);
  colorLookupTable->SetHueRange(.667, 0);
  colorLookupTable->SetRampToLinear();
  colorLookupTable->Build();

 // create the scalar map
  auto smoothdata = smooth(outputPolyData);


  vtkNew<vtkDataSetMapper> mapper;
  //mapper->SetInputConnection(reader->getGetOutputPort());
  mapper->SetInputData(smoothdata);
  CreateColor(outputPolyData);
  mapper->SelectColorArray("colors");
  //mapper->InterpolateScalarsBeforeMappingOn();
  //mapper->SetScalarRange(range);
  //mapper->SetLookupTable(colorLookupTable);

  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);

  //to display the mesh
  vtkSmartPointer<vtkExtractEdges> extractEdges =
    vtkSmartPointer<vtkExtractEdges>::New();
  extractEdges->SetInputConnection(reader->GetOutputPort());
  extractEdges->Update();

  // Visualize
  vtkSmartPointer<vtkPolyDataMapper> edgemapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
  edgemapper->SetInputConnection(extractEdges->GetOutputPort());
  vtkSmartPointer<vtkActor> edgeactor =
    vtkSmartPointer<vtkActor>::New();
  edgeactor->SetMapper(edgemapper);

  ////feature edge
  //vtkSmartPointer<vtkFeatureEdges> featureEdges =
  //  vtkSmartPointer<vtkFeatureEdges>::New();
  //featureEdges->SetInputConnection(reader->GetOutputPort());
  //featureEdges->BoundaryEdgesOn();
  //featureEdges->FeatureEdgesOff();
  //featureEdges->ManifoldEdgesOff();
  //featureEdges->NonManifoldEdgesOff();
  //featureEdges->Update();

  //// Visualize
  //vtkSmartPointer<vtkPolyDataMapper> fedgeMapper =
  //  vtkSmartPointer<vtkPolyDataMapper>::New();
  //fedgeMapper->SetInputConnection(featureEdges->GetOutputPort());
  //vtkSmartPointer<vtkActor> fedgeActor =
  //  vtkSmartPointer<vtkActor>::New();
  //fedgeActor->SetMapper(fedgeMapper);


  vtkSmartPointer<vtkRenderer> renderer =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  renderWindowInteractor->SetRenderWindow(renderWindow);

  renderer->AddActor(actor);
  //renderer->AddActor(fedgeActor);
  //renderer->AddActor(edgeactor);

  //renderer->SetBackground(0.1804, 0.5451, 0.3412); // Sea green

  renderWindow->Render();
  renderWindowInteractor->Start();

  return EXIT_SUCCESS;
}

vtkSmartPointer<vtkPolyData>
smooth(vtkSmartPointer<vtkPolyData> polydata) {

  int num_points = polydata->GetNumberOfPoints();

  vtkSmartPointer<vtkPolyData> polydata_out = vtkPolyData::New();
  polydata_out->DeepCopy(polydata);

  vtkSmartPointer<vtkIdList> point_cells = vtkIdList::New();
  vtkSmartPointer<vtkIdList> cell_points = vtkIdList::New();
  // This will hold a temporary copy of the points.
  vtkSmartPointer<vtkDoubleArray> tmp_DA = vtkDoubleArray::New();
  tmp_DA->SetNumberOfTuples(num_points);
  tmp_DA->SetNumberOfComponents(1);

  polydata->BuildLinks();

  int n_passes = 5;
  for (int i = 0; i < n_passes; ++i) {
    // Try to smooth the point scalars.
    for (int point_id = 0; point_id < num_points; ++point_id) {
      double point_scalar = polydata_out->GetPointData()->GetScalars()->GetTuple1(point_id);
      double *point_position = new double[3];
      point_position = polydata->GetPoint(point_id);

      double average_scalar = point_scalar;
      double den_fact = 1;

      polydata->GetPointCells(point_id, point_cells);
      vtkIdType num_cells = point_cells->GetNumberOfIds();

      vtkIdType cell;
      for (cell = 0; (cell < num_cells); ++cell) {
        vtkIdType neighbor_cell_id;
        neighbor_cell_id = point_cells->GetId(cell);
        polydata->GetCellPoints(neighbor_cell_id, cell_points);
        vtkIdType num_cell_points = cell_points->GetNumberOfIds();
        // Loop over the cell points.
        vtkIdType neighbor_point;
        for (neighbor_point = 0; neighbor_point < num_cell_points; ++neighbor_point) {
          // Get the neighbor point id.
          vtkIdType neighbor_point_id = cell_points->GetId(neighbor_point);
          double *neighbor_position = new double[3];
          neighbor_position = polydata->GetPoint(neighbor_point_id);

          double a = neighbor_position[0] - point_position[0];
          double b = neighbor_position[1] - point_position[1];
          double c = neighbor_position[2] - point_position[2];
          double distance = sqrt(a*a + b * b + c * c);
          double weight = (1 - distance / MAXDIST) * int(distance<MAXDIST);

          // Get the neighbor point position.
          double neighbor_scalar;
          neighbor_scalar = polydata_out->GetPointData()->GetScalars()->GetTuple1(neighbor_point_id);
          // Add it to the average.
          average_scalar += neighbor_scalar * weight;
          den_fact += weight;
        } // End neighbor points loop.
      } // End of cells loop
      average_scalar /= den_fact;
      tmp_DA->SetTuple1(point_id, average_scalar);
    } // End of smooth loop.
      //    tree_points->DeepCopy( tmp_points );
    polydata_out->GetPointData()->SetScalars(tmp_DA);
  } // End of itr loop
  return polydata_out;
}

void
CreateColor(vtkSmartPointer<vtkPolyData> outputPolyData)
{
  // Get the distances from the polydata
  vtkSmartPointer<vtkDoubleArray> stress = vtkDoubleArray::SafeDownCast(outputPolyData->GetPointData()->GetArray("vtkstress"));
  double minmax[2];
  stress->GetRange(minmax);

  // Create the color map
  vtkSmartPointer<vtkLookupTable> colorLookupTable =
    vtkSmartPointer<vtkLookupTable>::New();
  colorLookupTable->SetTableRange(minmax[0], minmax[1]);
  colorLookupTable->Build();

  // Generate the colors for each point based on the color map
  vtkSmartPointer<vtkUnsignedCharArray> colors =
    vtkSmartPointer<vtkUnsignedCharArray>::New();
  colors->SetNumberOfComponents(3);
  colors->SetName("Colors");

  std::cout << "There are " << outputPolyData->GetNumberOfPoints()
    << " points." << std::endl;

  for (int i = 0; i < outputPolyData->GetNumberOfPoints(); i++)
  {
    double p[3];
    outputPolyData->GetPoint(i, p);

    double dcolor[3];
    colorLookupTable->GetColor(stress->GetValue(i), dcolor);
    /*std::cout << "dcolor: "
      << dcolor[0] << " "
      << dcolor[1] << " "
      << dcolor[2] << std::endl;*/
    unsigned char color[3];
    for (unsigned int j = 0; j < 3; j++)
    {
      color[j] = static_cast<unsigned char>(255.0 * dcolor[j]);
    }
    //std::cout << "color: "
    //  << (int)color[0] << " "
    //  << (int)color[1] << " "
    //  << (int)color[2] << std::endl;

    colors->InsertNextTypedTuple(color);
  }
  outputPolyData->GetPointData()->SetAttribute(colors, vtkDataSetAttributes::SCALARS);
}