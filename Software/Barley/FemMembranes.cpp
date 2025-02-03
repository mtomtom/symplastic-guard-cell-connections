#include "FemMembranes.hpp"
#include <cmath>
#include <Geometry.hpp>
#include <iostream>
#include <fstream>
#include <MeshProcessStructure.hpp>
#include <MeshProcessPosition.hpp>

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

  //  bool DrawLines::initialize(QWidget *parent)
  //  {
  //    mesh = currentMesh();
  //    if(!mesh)
  //      throw QString("%1::initialize No current mesh").arg(name());
  //
  //    ccName = mesh->ccName();
  //    newCC = false;
  //    if(ccName.isEmpty()) {
  //      newCC = true;
  //      ccName = "Make Line";
  //    }
  //
  //    cs = &mesh->ccStructure(ccName);
  //    indexAttr = &mesh->indexAttr();
  //
  //    distance = parm("Distance").toDouble(); 
  //    drawLines = stringToBool(parm("Draw Lines")); 
  //
  //    // Set up mouse grabber
  //    grabsMouse = false;
  //    addInMouseGrabberPool();
  //
  //    w = parent;
  //
  //    if(newCC)
  //      mesh->setCCName(ccName);
  //
  //    return true;
  //  }
  //
  //  bool DrawLines::step()
  //  {
  //    if(QApplication::queryKeyboardModifiers().testFlag(Qt::AltModifier))
  //      grabsMouse = true;
  //    else
  //      grabsMouse = false;
  //
  //    return true;
  //  }
  //
  //  bool DrawLines::finalize(QWidget *parent)
  //  {
  //    grabsMouse = false;
  //    removeFromMouseGrabberPool();
  //
  //    return true;
  //  }
  //
  //  bool DrawLines::setGroupsVisible(bool nCC = true)
  //  {
  //    mesh->drawParms(ccName).setGroupVisible("Vertices", true);
  //    mesh->drawParms(ccName).setGroupVisible("Edges", true);
  //    newCC = nCC;
  //
  //    return true;
  //  }
  //
  //  bool DrawLines::addPoint(const Point3d &p)
  //  {
  //    CCIndex v = CCIndexFactory.getIndex();
  //    prevCell = v;
  //    auto &vIdx = (*indexAttr)[v];
  //    vIdx.pos = p;
  //    vIdx.label = label;
  //    cs->addCell(v);
  //
  //    if(newCC)
  //      setGroupsVisible();
  //
  //    prevCell = v;
  //    prevPos = p;
  //
  //    mesh->updateAll(ccName);
  //
  //    return true;
  //  }
  //
  //  bool DrawLines::addLine(const Point3d &p)
  //  {
  //    CCIndex v = CCIndexFactory.getIndex();
  //    auto &vIdx = (*indexAttr)[v];
  //    vIdx.pos = p;
  //    vIdx.label = label;
  //    cs->addCell(v);
  //    if(prevCell != CCIndex::UNDEF) {
  //      CCIndex e = CCIndexFactory.getIndex();
  //      cs->addCell(e, +prevCell -v);
  //      if(newCC)
  //        setGroupsVisible(false);
  //    }
  //
  //    prevCell = v;
  //    prevPos = p;
  //
  //    mesh->updateAll(ccName);
  //
  //    return true;
  //  }
  //
  //  void DrawLines::mousePressEvent(QMouseEvent *const event, Camera *const cam) 
  //  { 
  //
  //    drawing = true; 
  //
  //    // Get previous point to find depth
  //    Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos), &mesh->stack()->getFrame()));
  //
  //    // Get current position
  //    Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &mesh->stack()->getFrame()));
  //
  //    if(event->modifiers().testFlag(Qt::ShiftModifier))
  //      addLine(currPos);
  //    else {
  //      label = mesh->nextLabel();
  //      addPoint(currPos);
  //    }
  //  }
  //
  //  void DrawLines::mouseMoveEvent(QMouseEvent *const event, Camera *const cam) 
  //  {
  //    if(!drawing)
  //      return;
  //
  //    if(drawLines < 0)
  //      return;
  //
  //    Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos), &mesh->stack()->getFrame()));
  //    Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &mesh->stack()->getFrame()));
  //
  //    if(drawLines != 0 and norm(currPos - prevPos) < distance)
  //      return;
  //
  //    addLine(currPos);
  //  }
  //
  //  void DrawLines::mouseReleaseEvent(QMouseEvent *const event, Camera *const cam) 
  //  { 
  //    drawing = false; 
  //
  //    Point3d prevPoint(camera()->projectedCoordinatesOf(Vec(prevPos), &mesh->stack()->getFrame()));
  //    Point3d currPos(camera()->unprojectedCoordinatesOf(Vec(event->x(), event->y(), prevPoint.z()), &mesh->stack()->getFrame()));
  //
  //    addLine(currPos);
  //  }
  //
  //  void DrawLines::mouseDoubleClickEvent(QMouseEvent *const event, Camera *const cam) {}
  //  void DrawLines::wheelEvent(QWheelEvent *const event, Camera *const cam) {}
  //
  //  REGISTER_PROCESS(DrawLines);
  //
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

  bool CalculateAreaRay::run()
  {

    const int num_points = parm("Num Points").toInt(); // number of discretization points for the circle
    double poreArea = 1000;
    double poreWidth = 1000;
    double poreLength = 1000;

    // get the current mesh 
    mesh = currentMesh();
    if(!mesh)
      throw QString("%1::run Invalid mesh").arg(name());

    // get cc name
    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw QString("%1::initialize Invalid cell complex").arg(name());
    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();

    TransformMesh tm(*this);
    tm.run(cs, indexAttr, Point3d(0,0,-2), Point3d(0,0,0), 0.0, 1.0);

    for (int i = 0; i < 41; i++) {
        mdxInfo << -2+i*0.1 << endl;
        tm.run(cs, indexAttr, Point3d(0,0,.1), Point3d(0,0,0), 0.0, 1.0);
        AreaMin(cs, indexAttr, num_points, poreArea, poreWidth, poreLength);
    }

    mdxInfo << "Width: " << '\t' << poreWidth << endl;
    mdxInfo << "Length: " << '\t' << poreLength << endl;
    mdxInfo << "Area: " << '\t' << poreArea << endl;

    return true;
  }

  REGISTER_PROCESS(CalculateAreaRay);

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

  REGISTER_PROCESS(HeatMapPrint);

  bool GeomCSV::run()
  {
    mdxInfo << "Output geometrical data to a csv and save" << endl;

    // declare and initialize variables that will be used 
    const int num_points = parm("Num Points").toInt();                   // number of discretization points for the circle
    
    double poreArea = 1000;
    double poreWidth = 1000;
    double poreLength = 1000;
    
    double cost = 0;    
    double gc_vol_total = 0;    
    double gc_sa_total = 0;    
    double sc_vol_total = 0;    
    double sc_sa_total = 0;    

    // constants for cost function (From UoS Mean Geometry Values)
    // Open stoma parameters 
    const double gc_vol = 2369.2425;
    const double gc_sa = 1500.303542;
    const double sc_vol = 4780.564028;
    const double sc_sa = 2072.665972; 
    const double pore_area = 159.2802609;
    // Closed stoma parameters      
//    const double gc_vol = 2128.453125;
//    const double gc_sa = 1439.414583;
//    const double sc_vol = 4409.701667;
//    const double sc_sa = 1934.630417; 
//    const double pore_area = 93.18431579;
    
    
    std::string fileName = parm("File Name").toStdString();
    std::ofstream myfile;

    if (parm("Write Header") == "True")
    {
      myfile.open(fileName);
      myfile << "P_GC1,P_GC2,pore_width,pore_length,pore_area,stoma_length,stoma_width,Vol_GC2,Vol_SC1,Vol_SC2,Vol_GC1,SA_GC2,SA_SC1,SA_SC2,SA_GC1, cost\n";
    } else
    {
      myfile.open(fileName, std::ios_base::app);
    }

    // get the current mesh 
    mesh = currentMesh();
    if(!mesh)
      throw QString("%1::run Invalid mesh").arg(name());
   
    // get cc name 
    ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw QString("%1::initialize Invalid cell complex").arg(name());
    CCStructure &cs = mesh->ccStructure(ccName);
    CCIndexDataAttr &indexAttr = mesh->indexAttr();

    TransformMesh tm(*this);
    tm.run(cs, indexAttr, Point3d(0,0,-2), Point3d(0,0,0), 0.0, 1.0);

    for (int i = 0; i < 41; i++) {
        mdxInfo << -2+i*0.1 << endl;
        tm.run(cs, indexAttr, Point3d(0,0,.1), Point3d(0,0,0), 0.0, 1.0);
        AreaMin(cs, indexAttr, num_points, poreArea, poreWidth, poreLength);
    }

    tm.run(cs, indexAttr, Point3d(0,0,-2), Point3d(0,0,0), 0.0, 1.0);
    
    mdxInfo << "Wid: " << '\t' << poreWidth << endl;
    mdxInfo << "Len: " << '\t' << poreLength << endl;
    mdxInfo << "Area: " << '\t' << poreArea << endl;

    double xmin = 0;
    double xmax = 0;
    double ymin = 0;
    double ymax = 0;
    double length = 0;
    double width = 0;

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

    myfile << parm("GC Pressure").toDouble() << ", " << parm("SC Pressure").toDouble() << ", " << poreWidth << ", " << poreLength << ", " << poreArea << ", " << length << ", " << width; 

    double count = 0;
    auto &heatAttr = mesh->heatAttr<double>("Volume");
    for (auto &pr : heatAttr) {
      myfile << ", " << pr.second;
      if (count == 0 or count == 3) {
        gc_vol_total += pr.second;
      } else {
        sc_vol_total += pr.second;
      }
      count += 1;  
    }

    count = 0; 
    auto &heatAttr2 = mesh->heatAttr<double>("Cell Wall Area");
    for (auto &pr : heatAttr2) {
      myfile << ", " << pr.second;
      if (count == 0 or count == 3) {
        gc_sa_total += pr.second;
      } else {
        sc_sa_total += pr.second;
      }
      count += 1;  
    }
   
   // calculate cost function
    gc_vol_total = gc_vol_total / 2; 
    sc_vol_total = sc_vol_total / 2; 
    gc_sa_total = gc_sa_total / 2; 
    sc_sa_total = sc_sa_total / 2; 
    
    cost = pow(log(gc_vol_total/gc_vol),2) + pow(log(sc_vol_total/sc_vol),2) + pow(log(gc_sa_total/gc_sa),2) + pow(log(sc_sa_total/sc_sa),2) + pow(log(poreArea/pore_area),2);
   
    myfile << ", " << cost; 
    
    myfile << "\n";
    myfile.close();

    return true;
  }

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
  REGISTER_PROCESS(GeomCSV);

}
