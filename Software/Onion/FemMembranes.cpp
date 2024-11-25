#include "FemMembranes.hpp"
#include <cmath>
#include <Geometry.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <unordered_map>

#include <MeshProcessPosition.hpp>
#include <MeshProcessSelection.hpp>
#include <MeshProcessHeatMap.hpp>
#include <Process.hpp>

#include <CCUtils.hpp>
#include <MeshUtils.hpp>
#include <DistMatrix.hpp>
#include <KrylovMethods.hpp>

#include <limits>
#include <QMouseEvent>
#include <QApplication>

namespace mdx
{

  // Bubble sort for finding pore length and width
  void swap(double *xp, double *yp)
  {
    double temp = *xp;
    *xp = *yp;
    *yp = temp;
  }

  void bubblesort(double arr[], int n)
  {
    int i, j;
    for (i=0; i < n-1; i++)

      for (j=0; j < n-i-1; j++)
        if (arr[j] > arr[j+1])
          swap(&arr[j], &arr[j+1]);
  }

  

  bool AreaMin(const CCStructure cs, CCIndexDataAttr indexAttr, int num_points, double& poreArea, double& poreWidth, double& poreLength)
  {
    
    // declare and initialize variables that will be used 
    const double d_theta = 2*M_PI/num_points;       // d_theta
    double dist_squared[num_points];              // variable to store the distance to the nearest triangle 
    double theta[num_points]; 
    Point3d intp_pts[num_points];
    const Point3d origin;
    Point3d ray[num_points];
    Point3d intp[num_points];

    // rotate the ray that emanates from the origin
    // can this be parallelized??
#pragma omp parallel for
    for (int n = 0; n < num_points; n++) {
      dist_squared[n] = pow(100,2);    // reset the closest distance for each theta
      theta[n] = n*d_theta;
      ray[n] = Point3d(100*cos(theta[n]),100*sin(theta[n]),0);
      intp_pts[n] = intp[n];

      // loop through all the faces and keep the closest distance
      // needs to be nested in the loop that rotates the ray.
      for (int i = 0; i < int(cs.faces().size()); i++) {
        CCIndex f = cs.faces()[i];
        CCIndexVec nbs = cs.faceVertices(f); 
        CCIndexData &vIdx = indexAttr[nbs[0]];
        CCIndexData &mIdx = indexAttr[nbs[1]];
        CCIndexData &nIdx = indexAttr[nbs[2]];

        if ((rayTriangleIntersect(origin, ray[n], vIdx.pos, mIdx.pos, nIdx.pos, intp[n]) == 1) and 
            (pow(intp[n][0] - origin[0],2) + pow(intp[n][1] - origin[1],2) + pow(intp[n][2] - origin[2],2) < dist_squared[n])) {

          // save the point in a list replacing the previous closest point
          dist_squared[n] = pow(intp[n][0] - origin[0],2) + pow(intp[n][1] - origin[1],2) + pow(intp[n][2] - origin[2],2);
          intp_pts[n] = intp[n];
        }
      }
    } 

    // Calculate the area using polygonCentroidData
    //auto poreArea = polygonCentroidData(intp_pts).measure();
    // Calculate the triangle fan using the list of intersection points.
    // variables for area calculation simplicity
    double x0 = 0;
    double y0 = 0;
    double x1 = 0;
    double y1 = 0;
    double x2 = 0;
    double y2 = 0;
    double area = 0;

    for (int n = 0; n < num_points; n++) {
      x0 = intp_pts[n][0];
      y0 = intp_pts[n][1];

      x1 = intp_pts[(n+1)%num_points][0];
      y1 = intp_pts[(n+1)%num_points][1]; 

      area += std::abs(x0 * (y1 - y2) + x1 * (y2 - y0) + x2 * (y0 - y1)) / 2;
    }

    if (area < poreArea) {
        int n = sizeof(dist_squared)/sizeof(dist_squared[0]);
        bubblesort(dist_squared, n);

        poreArea = area;
        poreWidth = pow(dist_squared[0],0.5) + pow(dist_squared[1],0.5);
        poreLength = pow(dist_squared[num_points-2],0.5) + pow(dist_squared[num_points-1],0.5);
    }

  return true;
  }

  /// Set the anisotropy direction for selected elements

  //
  // Clinton: This is a bit more complicated than normal, because it is designed to work against different element types, so
  // here we just change the parameters based on the type and call the templated method in the .hpp file.
  //
  // You wouldn't need to do this in python.
  //
  using namespace fem; // FIXME when moving to FemLib
  bool FemSetAnisoDirLines::run(Mesh &mesh, CCStructure &cs, CCStructure &csLine, const QString &elementType, const QString &elementName, 
      int dirType, int projType, double tolerance)
  {
    CCIndexDataAttr &indexAttr = mesh.indexAttr();

    ElasticElementAttr eleAttr;
    if(!eleAttr.getAttr(mesh, elementType, elementName))
      throw QString("%1::run Invalid element type %1 or name %2").arg(elementType).arg(elementName);

    if(eleAttr.linearTriangle)
      run(cs, csLine, indexAttr, *eleAttr.linearTriangle, dirType, projType, tolerance);
    else if(eleAttr.linearWedge)
      run(cs, csLine, indexAttr, *eleAttr.linearWedge, dirType, projType, tolerance);
    else if(eleAttr.linearTetra)
      run(cs, csLine, indexAttr, *eleAttr.linearTetra, dirType, projType, tolerance);
    else
      throw QString("%1::run Invalid element type").arg(name());  

    return true;
  }
  REGISTER_PROCESS(FemSetAnisoDirLines);
  
  bool DrawLines::initialize(QWidget *parent)
  {
    mesh = currentMesh();
    stack = currentStack();
    if(!mesh)
      throw QString("%1::initialize No current mesh").arg(name());

    ccName = mesh->ccName();
    newCC = false;
    if(ccName.isEmpty()) {
      newCC = true;
      ccName = "Make Line";
    }

    cs = &mesh->ccStructure(ccName);
    indexAttr = &mesh->indexAttr();

    distance = parm("Distance").toDouble(); 
    drawLines = stringToBool(parm("Draw Lines")); 

    // Set up mouse grabber
    grabsMouse = false;
    addInMouseGrabberPool();

    w = parent;

    if(newCC)
      mesh->setCCName(ccName);

    return true;
  }

  bool DrawLines::step()
  {
    if(QApplication::queryKeyboardModifiers().testFlag(Qt::AltModifier))
      grabsMouse = true;
    else
      grabsMouse = false;

    return true;
  }

  bool DrawLines::finalize(QWidget *parent)
  {
    grabsMouse = false;
    removeFromMouseGrabberPool();

    return true;
  }

  bool DrawLines::setGroupsVisible(bool nCC = true)
  {
    mesh->drawParms(ccName).setGroupVisible("Vertices", true);
    mesh->drawParms(ccName).setGroupVisible("Edges", true);
    newCC = nCC;

    return true;
  }

  bool DrawLines::addPoint(const Point3d &p)
  {
    CCIndex v = CCIndexFactory.getIndex();
    prevCell = v;
    auto &vIdx = (*indexAttr)[v];
    vIdx.pos = p;
    vIdx.label = label;
    cs->addCell(v);

    if(newCC)
      setGroupsVisible();

    prevCell = v;
    prevPos = p;

    mesh->updateAll(ccName);

    return true;
  }

  bool DrawLines::addLine(const Point3d &p)
  {
    CCIndex v = CCIndexFactory.getIndex();
    auto &vIdx = (*indexAttr)[v];
    vIdx.pos = p;
    vIdx.label = label;
    cs->addCell(v);
    if(prevCell != CCIndex::UNDEF) {
      CCIndex e = CCIndexFactory.getIndex();
      cs->addCell(e, +prevCell -v);
      if(newCC)
        setGroupsVisible(false);
    }

    prevCell = v;
    prevPos = p;

    mesh->updateAll(ccName);

    return true;
  }

    void DrawLines::mousePressEvent(QMouseEvent *const event, Camera *const cam) 
    { 
  
      drawing = true; 
  
      // Get previous point to find depth
      Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos), &stack->getFrame()));
  
      // Get current position
      Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &stack->getFrame()));
  
      if(event->modifiers().testFlag(Qt::ShiftModifier))
        addLine(currPos);
      else {
        label = mesh->nextLabel();
        addPoint(currPos);
      }
    }
  
    void DrawLines::mouseMoveEvent(QMouseEvent *const event, Camera *const cam) 
    {
      if(!drawing)
        return;
  
      if(!drawLines)
        return;
  
      Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos), &stack->getFrame()));
      Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &stack->getFrame()));
  
      if(norm(currPos - prevPos) < distance)
        return;
  
      addLine(currPos);
    }
  
    void DrawLines::mouseReleaseEvent(QMouseEvent *const event, Camera *const cam) 
    { 
      drawing = false; 
  
      Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos),&stack->getFrame()));
      Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &stack->getFrame()));
  
      addLine(currPos);
    }
  
    void DrawLines::mouseDoubleClickEvent(QMouseEvent *const event, Camera *const cam) {}
    void DrawLines::wheelEvent(QWheelEvent *const event, Camera *const cam) {}
  
    REGISTER_PROCESS(DrawLines);
  
  bool MakeLineBezier::run(CCStructure &cs, CCIndexDataAttr &indexAttr, int points)
  {
    // Get the cutting surface
    CuttingSurface* cutSurf = cuttingSurface();
    if(!cutSurf)
      throw QString("%1::run Invalid Bezier").arg(name());

    // Create first vertex 
    CCIndex pv = CCIndexFactory.getIndex();
    // Add to cell complex
    cs.addCell(pv);
    // Set pos from bezier
    indexAttr[pv].pos = cutSurf->bezier().evalCoord(0.0, 0.0);

    double du = 1.0/points;
    for(int i = 1; i < points; i++) {
      // Create next vertex
      CCIndex p = CCIndexFactory.getIndex();
      cs.addCell(p);
      // Create edge
      indexAttr[p].pos = cutSurf->bezier().evalCoord(i * du, 0.0);
      CCIndex e = CCIndexFactory.getIndex();
      cs.addCell(e, +pv -p);
      // save previous vertex
      pv = p;
    }
    return true;
  }

  REGISTER_PROCESS(MakeLineBezier);

  bool StomaDims::run()
  {
    double xmin = 0;
    double xmax = 0;
    double ymin = 0;
    double ymax = 0;
    double length = 0;
    double width = 0;

    // get the current mesh 
    mesh = currentMesh();
    if(!mesh)
      throw QString("%1::run Invalid mesh").arg(name());
    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw QString("%1::initialize Invalid cell complex").arg(name());
    CCStructure cs = mesh->ccStructure(ccName);
    CCIndexDataAttr indexAttr = mesh->indexAttr();

    // loop through all the vertices and save max/min x/y of all faces 

    for (int i = 0; i < int(cs.vertices().size()); i++) {
      CCIndex v = cs.vertices()[i];
      auto &vIdx = indexAttr[v];

      if (vIdx.pos[0] > xmax)
        xmax = vIdx.pos[0];
      if (vIdx.pos[0] < xmin)
        xmin = vIdx.pos[0];
      if (vIdx.pos[1] > ymax)
        ymax = vIdx.pos[1];
      if (vIdx.pos[1] < ymin)
        ymin = vIdx.pos[1];
    }

    if (std::abs(xmax - xmin) >  std::abs(ymax - ymin)) {
      length = std::abs(xmax-xmin);
      width = std::abs(ymax-ymin);
    } else{
      width = std::abs(xmax-xmin);
      length = std::abs(ymax-ymin);
    }

    mdxInfo << "Stoma Length: " << length << endl;
    mdxInfo << "Stoma Width: " << width << endl;

    return true;
  }

  REGISTER_PROCESS(StomaDims);

  bool HeatMapPrint::run(Mesh &mesh, const QString &heatName, const QString &labeling, bool useSelection, bool exportLabelings) {

    CCIndexDataAttr *indexAttr = 0;
    CCStructure *cs = 0;
    if(useSelection) {
      QString ccName = mesh.ccName();
      if(ccName.isEmpty())
        throw QString("%1: No cell complex").arg(name());
      cs = &mesh.ccStructure(ccName);
      indexAttr = &mesh.indexAttr();
    }
    QStringList heatAttrList = heatName == "All" ? mesh.heatAttrList(labeling) : QStringList() << heatName;

    // Get all the labels and check the type
    IntSet labels;
    std::vector<IntDoubleAttr *> heatAttrs;
    std::vector<QString> heatNames;
    for(auto &heat : heatAttrList)
      if(mesh.heatType(heat, labeling) == "Double") {
        heatNames.push_back(heat);
        auto &heatAttr = mesh.heatAttr<double>(heat, labeling);
        if(!useSelection) // If not using the selection, get labels from attributes
          for(auto &pr : heatAttr)
            labels.insert(pr.first);

        heatAttrs.push_back(&heatAttr);
      }
    if(heatNames.size() == 0)
      throw QString("%1: No double heat maps available").arg(name());

    // get labels from selection
    if(useSelection) {
      if(!cs or !indexAttr)
        throw QString("%1: Cell complex not set").arg(name());

      for(CCIndex c : selectedCells(*cs, *indexAttr))
        labels.insert((*indexAttr)[c].label);
    }

    std::vector<QString> labelings;
    std::vector<IntIntAttr *> labelingAttrs;
    std::vector<IntQStringAttr *> labelingNames;
    if(exportLabelings) {
      for(auto &labeling :  mesh.labelingAttrList()) {
        if(labeling != "Labels") {
          auto *attr = mesh.labelMap(labeling);
          if(attr->size() > 0) {
            labelings.push_back(labeling);
            labelingAttrs.push_back(attr);
            labelingNames.push_back(&mesh.labelName(labeling));
          }
        }
      }
      if(labelings.size() == 0)
        exportLabelings = false;
    }

    //    // Write header
    //    mdxInfo << QString("%1,Description").arg(labeling);
    //    if(exportLabelings)
    //      for(uint i = 0; i < labelings.size(); i++)
    //        mdxInfo << QString(",%1,%1-Label").arg(labelings[i]);
    //
//    for(auto &heat : heatNames)
//      mdxInfo << "," << heat;
//    mdxInfo << endl;
    //
    auto &labelName = mesh.labelName(labeling);
    for(int label : labels) {
      //      mdxInfo << label << "," << labelName[label];
      //      mdxInfo << label << ",";
      //
      //      if(exportLabelings) {
      //        for(uint i = 0; i < labelings.size(); i++) {
      //          auto lp = labelingAttrs[i]->find(label);
      //          if(lp != labelingAttrs[i]->end())
      //            mdxInfo << "," << (*labelingNames[i])[lp->second] << "," << lp->second;
      //          else
      //            mdxInfo << ",,";
      //        }
      //      }
      for(uint i = 0; i < heatNames.size(); i++) {
        //        mdxInfo << ",";
        auto itr = heatAttrs[i]->find(label);
        if(itr != heatAttrs[i]->end()){
          mdxInfo << itr->second;
        } 
      }
      mdxInfo << endl;
    }

    setStatus(QString("Heat map written to terminal with %1 rows").arg(labels.size()));

    return true;
  }

// Function to read parameters from a file
void readParametersFromFile(const std::string& filePath, std::unordered_map<std::string, QString>& parameters) {
    std::ifstream file(filePath);
    
    if (!file.is_open()) {
        std::cerr << "Error opening file: " << filePath << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string name;
        std::string value;
        if (iss >> name >> value) {
            parameters[name] = QString::fromStdString(value);
        }
    }

    file.close();
}

QString findParameter(const QString& parameterName, const std::unordered_map<std::string, QString>& parameters) {
    auto it = parameters.find(parameterName.toStdString());
    if (it == parameters.end()) {
        throw QString("Parameter not found: ") + parameterName;
    }
    return it->second;
}

 // Initialize the simulation with the given parameters file
bool InflateStomata::initialiseSim(QString dataFolder, QString meshFile)
{
    std::unordered_map<std::string, QString> parameters;
    readParametersFromFile(parm("Parameters File").toStdString(), parameters);
    QStringList parms;
    QStringList names;

    // Set the correct parameters from the parameter file
    QString youngE1E3 = findParameter("young_E1E3", parameters);
    QString youngE2 = findParameter("young_E2", parameters);
    QString poisson = findParameter("poisson", parameters);
    QString thickness = findParameter("thickness", parameters);
    QString scalingFactor = findParameter("z_scale", parameters);
    QString dirichlet = dataFolder + findParameter("dirichlet", parameters);

    // Reset the mesh system
    getProcessParms("Mesh/System/Reset", parms);
    if (!runProcess("Mesh/System/Reset", parms))
    {
        throw QString("Error with resetting mesh: ") + errorMessage();
    }

    // Open the mesh system
    getProcessParms("Mesh/System/Open", parms);
    MeshOpen mso(*this);
    names = mso.parmNames();
    int fileName = names.indexOf("File Name");

    // Load the mesh with anisotropy direction already set
    parms[fileName] = dataFolder + meshFile;
    if (!runProcess("Mesh/System/Open", parms))
    {
        throw QString("Error with opening mesh: ") + errorMessage();
    }

    // Rescale the mesh based on the scaling factor
    ScaleMesh scaleMesh(*this);
    scaleMesh.setParm("Z Scale", scalingFactor);

    // Select all faces in the mesh
    MeshSelectAll selectFaces(*this);
    selectFaces.setParm("Dimension", "Faces");
    selectFaces.run();

    // Set the reference configuration for the model
    getProcessParms("Model/CCF/03 Reference Configuration", parms);
    if (!runProcess("Model/CCF/03 Reference Configuration", parms))
    {
        throw QString("Error with ref config: ") + errorMessage();
    }

    // Set the stress-strain properties for the model
    getProcessParms("Model/CCF/04 StressStrain", parms);
    if (!runProcess("Model/CCF/04 StressStrain", parms))
    {
        throw QString("Error with StressStrain: ") + errorMessage();
    }

    // Set the material properties for the model
    getProcessParms("Model/CCF/05 Set Material Properties", parms);
    FemMembraneSetMaterial fsm(*this);
    names = fsm.parmNames();
    int younge1e3 = names.indexOf("Young E1E3");
    parms[younge1e3] = youngE1E3;
    int younge2 = names.indexOf("Young E2");
    parms[younge2] = youngE2;
    if (!runProcess("Model/CCF/05 Set Material Properties", parms))
    {
        throw QString("Error with Material properties: ") + errorMessage();
    }

    // Set the anisotropy direction from the polygon mesh
    getProcessParms("Model/CCF/32 Set Ansio Dir From Lines", parms);
    FemSetAnisoDirLines fsadl(*this);
    names = fsadl.parmNames();
    int linecc = names.indexOf("Line CC");
    parms[linecc] = "Polygon";
    if (!runProcess("Model/CCF/32 Set Ansio Dir From Lines", parms))
    {
        throw QString("Model/CCF/32 Set Ansio Dir From Lines") + errorMessage();
    }

    // Clear the mesh selection
    MeshClearSelection clearSelection(*this);
    clearSelection.run();

    // Set the boundary conditions
    MeshLoadSelection loadSelect(*this);
    loadSelect.setParm("File Name", dirichlet);
    loadSelect.run();
    getProcessParms("Model/CCF/06 Set Dirichlet Boundary Conditions", parms);
    FemMembraneSetDirichlet fsd(*this);
    fsd.run();
    clearSelection.run();

    return true;
}
bool InflateStomata::loadAndSet3DCellPressure(QString cellFileName, double pressure)
{
    QStringList parms;
    QStringList names;

    // Ensure we have nothing else selected
    MeshClearSelection clearSelection(*this);
    clearSelection.run();

    // Load in the cell file
    MeshLoadSelection loadSelect(*this);
    loadSelect.setParm("File Name", cellFileName);
    loadSelect.run();

    // Set 3D cell pressure
    getProcessParms("Model/CCF/90 Set 3D Cell Pressure", parms);
    FemMembraneSet3DCellPressure fm3dp(*this);
    names = fm3dp.parmNames();
    int pressureIndex = names.indexOf("Pressure");
    parms[pressureIndex] = QString::number(pressure);
    int signalIndex = names.indexOf("Pressure Signal");
    parms[signalIndex] = "Signal Name";

    if (!runProcess("Model/CCF/90 Set 3D Cell Pressure", parms))
    {
        throw QString("Error with 3D Cell Pressure: ") + errorMessage();
    }

    // Set Face Pressure from volumes
    getProcessParms("Model/CCF/91 Set Face Pressure From Volumes", parms);
    if (!runProcess("Model/CCF/91 Set Face Pressure From Volumes", parms))
    {
        throw QString("Error with set face pressure from volumes: ") + errorMessage();
    }

    return true;
}

 bool InflateStomata::run(QString dataFolder, QString meshFile, bool differential)
 {
  QStringList parms;
  QStringList names;

  // Initialise Simulation
  initialiseSim(dataFolder, meshFile);

  // If the left and right cell need to be different pressures
  if (differential){

    // The left and right cell labels should be in the dataFolder
    QString leftFile = dataFolder + "left_cell.txt";
    QString rightFile = dataFolder + "right_cell.txt";

    // Set pressure of left cell
    loadAndSet3DCellPressure(leftFile, parm("Left cell pressure").toDouble());
    //Set pressure of right cell
    loadAndSet3DCellPressure(rightFile, parm("Right cell pressure").toDouble());
  }
  else {

    MeshClearSelection clearSelection(*this);
    clearSelection.run();
    // Select all volumes
    MeshSelectAll selectVolumes(*this);
    selectVolumes.setParm("Dimension", "Volumes");
    selectVolumes.run();

  
    // Set 3D cell pressure
    getProcessParms("Model/CCF/90 Set 3D Cell Pressure", parms);
    FemMembraneSet3DCellPressure fm3dp(*this);
    names = fm3dp.parmNames();
    int pressureIndex = names.indexOf("Pressure");
    parms[pressureIndex] = parm("Left cell pressure");
    int signalIndex = names.indexOf("Pressure Signal");
    parms[signalIndex] = "Signal Name";

    if (!runProcess("Model/CCF/90 Set 3D Cell Pressure", parms))
    {
        throw QString("Error with 3D Cell Pressure: ") + errorMessage();
    }

    // Set Face Pressure from volumes
    getProcessParms("Model/CCF/91 Set Face Pressure From Volumes", parms);
    if (!runProcess("Model/CCF/91 Set Face Pressure From Volumes", parms))
    {
        throw QString("Error with set face pressure from volumes: ") + errorMessage();
    }

    
  }

  // Run the FEM model
  getProcessParms("Model/CCF/01 FEM Membranes", parms);
  if (!runProcess("Model/CCF/01 FEM Membranes", parms))
  {
    throw QString("Error with FEM Membranes: ") + errorMessage();
  }

  // Print the results
  PoreArea pa(*this);
  pa.run();
  mesh = currentMesh();
  if(!mesh)
    throw QString("%1::run Invalid mesh").arg(name());

  AttrMap<QString,double> &values = mesh->attributes().attrMap<QString,double>("savedConstants");

  mdxInfo << "PLeft: " << '\t' << parm("Left cell pressure") << endl;
  mdxInfo << "PRight: " << '\t' << parm("Left cell pressure") << endl;
  mdxInfo << "Pore area: " << '\t' << values["Area"] << endl;

  return true;
    
  }

  // Function to run the individual simulations for RunSimulationsFigure2b
 bool RunSymplasticConnectionsFigure::DifferentialPressureSims(double pressureStart, double pressureEnd, double pressureInc, QString outputFolder) {
    std::vector<double> pressure;
    QStringList parms;
    QString meshFile = "onion_aniso.mdxm";

    // Initiate the inflate stomata process
    InflateStomata is(*this);

    for (double i = pressureStart; i <= pressureEnd; i += pressureInc) {
        pressure.push_back(i);
    }

    std::vector<double> areas(pressure.size());

    for (std::vector<double>::size_type i = 0; i < pressure.size(); i++) {
        double pLeft = pressure[i];
        double pRight = pressureEnd - pLeft;

        // Update the InflateStomata parameters
        is.setParm("Left cell pressure", QString::number(pLeft));
        is.setParm("Right cell pressure", QString::number(pRight));
        is.run(parm("Data Folder"), meshFile, true);

        // Record pore area
        PoreArea pa(*this);
        pa.run();
        mesh = currentMesh();
        if (!mesh)
            throw QString("%1::run Invalid mesh").arg(name());

        AttrMap<QString, double> &values = mesh->attributes().attrMap<QString, double>("savedConstants");
        areas[i] = values["Area"];
    }

    for (std::vector<double>::size_type j = 0; j < pressure.size(); j++) {
        mdxInfo << "PLeft: " << '\t' << pressure[j] << endl;
        mdxInfo << "PRight: " << '\t' << pressureEnd - pressure[j] << endl;
        mdxInfo << "Area: " << '\t' << areas[j] << endl;
    }

    // Write the outputs into a csv file
    std::ofstream myfile;
    QString filename = outputFolder + "output_" + QString::number(pressureEnd) + ".csv";
    myfile.open(filename.toStdString());
    myfile << "Left GC Pressure,Right GC Pressure,Area\n";
    for (std::vector<double>::size_type j = 0; j < pressure.size(); j++) {
        myfile << pressure[j] << "," << pressureEnd - pressure[j] << "," << areas[j] << "\n";
    }
    myfile.close();

    return true;
}

 // Create the process to output the images for the symplastic connections paper
 bool RunSymplasticConnectionsFigure::run()
 {
  QString figure3b = QString("3b");
  QString figure3c = QString("3c");
  QString figureS2a = QString("S2a"); // Pore area vs Pressure
  QString figureS2b = QString("S2b"); // Pore area for all meshes
  // The left and right cell labels should be in the dataFolder
  QString leftFile = dataFolder + "left_cell.txt";
  QString rightFile = dataFolder + "right_cell.txt";

  // Load in the parameters from the file
  std::unordered_map<std::string, QString> parameters;
  readParametersFromFile(parm("Parameters File").toStdString(), parameters);

  // Set up loop to run each of the simulations
  int startPressure = findParameter("start_pressure", parameters).toInt();
  int endPressure = findParameter("end_pressure", parameters).toInt();
  //QString lCell = parm("Data Folder") + findParameter("left_cell", parameters);
  //QString rCell = parm("Data Folder") + findParameter("right_cell", parameters);

  // Check if folder called "output" exists, if not, create it
  QDir dir(parm("Output Folder"));
  if (!dir.exists())
  {
    dir.mkpath(".");
  }

    if (parm("Figure")==figure3b)
  {
    for (int p = startPressure; p <= endPressure; p++)
    {
      // Calculate the increments
      double inc = static_cast<double>(p) / 10;
      // Run the differential pressure simulations
    DifferentialPressureSims(inc, p, inc, parm("Output Folder"));
    }
  }

  else if(parm("Figure")==figure3c)
  {
    // Total pressure = 5MPa
  
     // Set up the left and right cell pressures
    double lCell = 5.0 * parm("Left cell pressure proportion (figure 3c)").toDouble();
    double rCell = 5.0 - lCell;

    InflateStomata is(*this);
    QString meshFile = "onion_aniso.mdxm";
    is.initialiseSim(dataFolder, meshFile);
    MeshClearSelection clearSelection(*this);
    clearSelection.run();
    // Select all volumes
    MeshSelectAll selectVolumes(*this);
    selectVolumes.setParm("Dimension", "Volumes");
    selectVolumes.run();

  
    // Set pressure of left cell
    is.loadAndSet3DCellPressure(leftFile, lCell);
    //Set pressure of right cell
    is.loadAndSet3DCellPressure(rightFile, rCell);

    //Run the model
    getProcessParms("Model/CCF/01 FEM Membranes", parms);
    if (!runProcess("Model/CCF/01 FEM Membranes", parms))
    {
        throw QString("Error with FEM Membranes: ") + errorMessage();
    }

    
  }
  // Run the FEM model
  getProcessParms("Model/CCF/01 FEM Membranes", parms);
  if (!runProcess("Model/CCF/01 FEM Membranes", parms))
    {
      throw QString("Error with FEM Membranes: ") + errorMessage();
    }

  }

  else if (parm("Figure")==figureS2a)
  {
    // Pore area vs pressure

    // Initiate the inflate stomata process
    InflateStomata is(*this);

    QString meshFile = "onion_aniso.mdxm";
    double area;
    int maxPressure = 10;

    // Write the outputs into a csv file
    std::ofstream myfile;
    QString filename = parm("Output Folder") + "S2a.csv";
    myfile.open(filename.toStdString());
    myfile << "Left GC Pressure,Right GC Pressure,Area\n";


    for (int p =0; p < maxPressure; p++)
      {

        // Update the InflateStomata parameters
        is.setParm("Left cell pressure", QString::number(p));
        is.setParm("Right cell pressure", QString::number(p));
        is.run(parm("Data Folder"), meshFile, false);

        // Record pore area
        PoreArea pa(*this);
        pa.run();
        mesh = currentMesh();
        if(!mesh)
          throw QString("%1::run Invalid mesh").arg(name());

        AttrMap<QString,double> &values = mesh->attributes().attrMap<QString,double>("savedConstants");
        area = values["Area"];
    
        myfile << p << "," << p << "," << area << "\n";
      }
      myfile.close();
        
    }

  else if (parm("Figure")==figureS2b)
  {
    // Run the equal pressure simulations for all of the meshes
    // Hard code the list in for now - ensure that the same list of meshes is used as in the paper
    QStringList meshes = {"1_2","1_3","1_4","1_5","1_6","1_8","2_1","2_3","2_4","2_6a","2_6b","2_7a","3_1","3_2","3_3","3_4","3_6","3_7"};

    // Set the pressure
    int p = 5;

    QString meshFile = "onion_aniso.mdxm";
    double area;

    // Write the outputs into a csv file
    std::ofstream myfile;
    QString filename = parm("Output Folder") + "S2b.csv";
    myfile.open(filename.toStdString());
    myfile << "Mesh,Pressure,Area\n";

    // Need to amend the data folder to include the correct mesh folder
    QString dataFolder = parm("Data Folder");
    dataFolder.replace("/1_2/", "/");

    // Initiate the inflate stomata process
    InflateStomata is(*this);

    // Loop through all of the meshes

    for (int m =0; m < meshes.count(); m++)
    {
      
      is.setParm("Left cell pressure", QString::number(p));
      is.setParm("Right cell pressure", QString::number(p));
      // Update the InflateStomata parameters
      is.run(dataFolder + meshes[m] + "/", meshFile, false);

      // Record pore area
      PoreArea pa(*this);
      pa.run();
      mesh = currentMesh();
      if(!mesh)
        throw QString("%1::run Invalid mesh").arg(name());

      AttrMap<QString,double> &values = mesh->attributes().attrMap<QString,double>("savedConstants");
      area = values["Area"];
    
      mdxInfo << "Mesh: " << '\t' << meshes[m] << endl;
      mdxInfo << "Pressure: " << '\t' << p << endl;
      mdxInfo << "Area: " << '\t' << area << endl;
  
      myfile << meshes[m] << "," << p << "," << area << "\n";
    }

    myfile.close();

  }  

  return true;

 }

  REGISTER_PROCESS(HeatMapPrint);
  REGISTER_PROCESS(FemMembranes);
  REGISTER_PROCESS(FemMembraneRefCfg);
  REGISTER_PROCESS(FemMembraneStressStrain);
  REGISTER_PROCESS(FemMembraneDerivs);
  REGISTER_PROCESS(FemMembraneSetMaterial);
  REGISTER_PROCESS(FemMembraneAnisoDir);
  REGISTER_PROCESS(FemMembraneSetPressure);
  REGISTER_PROCESS(FemMembraneSetDirichlet);
  REGISTER_PROCESS(FemMembranePressureDerivs);
  REGISTER_PROCESS(FemMembraneDirichletDerivs);
  REGISTER_PROCESS(FemMembraneVisMaterial);
  REGISTER_PROCESS(FemMembraneVisDirections);
  REGISTER_PROCESS(FemMembraneVisDirichlet);
  REGISTER_PROCESS(FemMembraneSet3DCellPressure);
  REGISTER_PROCESS(FemMembraneSetFacePressureFromVolumes);
  REGISTER_PROCESS(InflateStomata);
  REGISTER_PROCESS(RunSymplasticConnectionsFigure);
}
