#ifndef FEM_MEMBRANES_HPP
#define FEM_MEMBRANES_HPP

#include <MDXProcessFem.hpp>
#include <MeshProcessSystem.hpp>

#include <Process.hpp>

namespace mdx
{
  class FemMembranes : public fem::FemSolver
  {
    public:
      FemMembranes(const Process &proc) : FemSolver(proc) 
    {
      setName("Model/CCF/01 FEM Membranes");
      setDesc("FEM Simulation using triangular membrane elements");

      // Update parameters with our own defaults
      setParmDefault("Stress-Strain", "Model/CCF/04 StressStrain");

      // Add derivatives processes
      addParm("Element Derivs", "Process for element derivatives", "Model/CCF/02 Triangle Derivs");
      addParm("Pressure Derivs", "Process for pressure derivatives", "Model/CCF/08 Pressure Derivs");
      addParm("Dirichlet Derivs", "Process for Dirichlet derivatives", "Model/CCF/10 Dirichlet Derivs");
    }

      bool rewind(QWidget *parent)
      {
        // To rewind, we'll reload the mesh
        Mesh *mesh = currentMesh();
        if(!mesh or mesh->file().isEmpty())
          throw(QString("No current mesh, cannot rewind"));
        MeshLoad meshLoad(*this);
        meshLoad.setParm("File Name", mesh->file());
        return meshLoad.run();
      }
  };

  class FemMembraneDerivs : public fem::ElementDerivs
  {
    public:
      FemMembraneDerivs(const Process &proc) : ElementDerivs(proc) 
    {
      setName("Model/CCF/02 Triangle Derivs");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneRefCfg : public fem::SetRefCfg
  {
    public:
      FemMembraneRefCfg(const Process &proc) : SetRefCfg(proc) 
    {
      setName("Model/CCF/03 Reference Configuration");

      addParm("Thickness", "Thickness of the membrane elements", "1.0");
      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneStressStrain : public fem::StressStrain
  {
    public:
      FemMembraneStressStrain(const Process &proc) : StressStrain(proc) 
    {
      setName("Model/CCF/04 StressStrain");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneSetMaterial : public fem::SetTransIsoMaterial
  {
    public:
      FemMembraneSetMaterial(const Process &proc) : SetTransIsoMaterial(proc) 
    {
      setName("Model/CCF/05 Set Material Properties");
    }
  };

  class FemMembraneAnisoDir : public fem::SetAnisoDir
  {
    public:
      FemMembraneAnisoDir(const Process &proc) : SetAnisoDir(proc) 
    {
      setName("Model/CCF/06 Set Ansio Dir");

      setParmDefault("Element Type", "Linear Triangle");
      setParmDefault("Element Attribute", "Triangle Element");
    }
  };

  class FemMembraneSetPressure : public fem::SetPressure
  {
    public:
      FemMembraneSetPressure(const Process &proc) : SetPressure(proc) 
    {
      setName("Model/CCF/07 Set Pressure");
    }
  };

  class FemMembraneSet3DCellPressure : public fem::Set3DCellPressure
  {
    public:
      FemMembraneSet3DCellPressure(const Process &proc) : Set3DCellPressure(proc) 
    {
      setName("Model/CCF/90 Set 3D Cell Pressure");
    }
  };

  class FemMembraneSetFacePressureFromVolumes : public fem::SetFacePressureFromVolumes
  {
    public:
      FemMembraneSetFacePressureFromVolumes(const Process &proc) : SetFacePressureFromVolumes(proc) 
    {
      setName("Model/CCF/91 Set Face Pressure From Volumes");
    }
  };

  class FemMembranePressureDerivs : public fem::PressureDerivs
  {
    public:
      FemMembranePressureDerivs(const Process &proc) : PressureDerivs(proc) 
    {
      setName("Model/CCF/08 Pressure Derivs");
    }
  };

  class FemMembraneSetDirichlet : public fem::SetDirichlet
  {
    public:
      FemMembraneSetDirichlet(const Process &proc) : SetDirichlet(proc)
    {
      setName("Model/CCF/09 Set Dirichlet");
    }
  };

  class FemMembraneDirichletDerivs : public fem::DirichletDerivs
  {
    public:
      FemMembraneDirichletDerivs(const Process &proc) : DirichletDerivs(proc)
    {
      setName("Model/CCF/10 Dirichlet Derivs");
    }
  };

  class FemMembraneVisMaterial : public fem::VisTransIsoMaterial
  {
    public:
      FemMembraneVisMaterial(const Process &proc) : VisTransIsoMaterial(proc) 
    {
      setName("Model/CCF/20 Visualize Material");
    }
  };

  class FemMembraneVisDirections : public fem::VisDirections
  {
    public:
      FemMembraneVisDirections(const Process &proc) : VisDirections(proc) 
    {
      setName("Model/CCF/21 Visualize Directions");
    }
  };

  class FemMembraneVisDirichlet : public fem::VisDirichlet
  {
    public:
      FemMembraneVisDirichlet(const Process &proc) : VisDirichlet(proc)
    {
      setName("Model/CCF/22 Visualize Dirichlet");
    }
  };

  //  // Process to draw lines, uses QGLViewer mousegrabber
  //  using qglviewer::Camera;
  //  using qglviewer::Vec;
  //
  //  class DrawLines : public Process, public qglviewer::MouseGrabber
  //  {
  //  public:
  //    DrawLines(const Process &process) : Process(process) 
  //    {
  //      setName("Model/CCF/30 Draw Lines");
  //      setDesc("Draw lines with the mouse");
  //      setIcon(QIcon(":/images/DrawLines.png"));
  //
  //      addParm("Distance", "Distance (um) between points when moving mouse, 0 grabs all, -1 none", "10.0");
  //      addParm("Draw Lines", "Draw lines and points", "Yes", booleanChoice());
  //    }
  //    bool initialize(QWidget *parent);
  //    bool step();
  //    bool finalize(QWidget *parent);
  //
  //    bool setGroupsVisible(bool newCC);
  //    bool addPoint(const Point3d &p);
  //    bool addLine(const Point3d &p);
  //
  //    void checkIfGrabsMouse(int x, int y, const Camera* const cam) { setGrabsMouse(grabsMouse); }
  //    void mousePressEvent(QMouseEvent *const event, Camera *const cam);
  //    void mouseDoubleClickEvent(QMouseEvent *const event, Camera *const cam);
  //    void mouseReleaseEvent(QMouseEvent *const event, Camera *const cam);
  //    void mouseMoveEvent(QMouseEvent *const event, Camera *const cam);
  //    void wheelEvent(QWheelEvent *const event, Camera *const cam);
  //
  //  private:
  //    Mesh *mesh = 0;
  //    QString ccName;
  //    CCStructure *cs = 0;
  //    CCIndexDataAttr *indexAttr = 0;
  //
  //    double distance = 0;
  //    bool drawLines = true; 
  //
  //    Point3d prevPos;
  //    CCIndex prevCell = CCIndex::UNDEF;
  //    bool drawing = false;
  //    bool newCC = false;
  //    int label = 0;
  //
  //    bool grabsMouse = false;
  //    QWidget *w = 0;
  //  };
  //
  class MakeLineBezier : public Process
  {
    public:
      MakeLineBezier(const Process &process) : Process(process) 
    {
      setName("Model/CCF/31 Make Line Bezier");
      setDesc("Make a line from the Bezier");

      addParm("CC Name", "Name of cell complex to create lines", "Bezier Line");
      addParm("Segments", "Number of segments in the line", "100");
    }
      bool run()
      {
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw QString("%1::run No mesh").arg(name());

        QString ccName = parm("CC Name");
        if(ccName.isEmpty())
          throw QString("%1::run Cell complex parameter empty").arg(name());

        auto &cs = mesh->ccStructure(ccName);
        auto &indexAttr = mesh->indexAttr();
        int segments = parm("Segments").toInt();
        if(segments < 2)
          throw QString("%1::run Points must be at least 2").arg(name());

        mesh->updateAll(ccName);
        return run(cs, indexAttr, segments);
      }
      bool run(CCStructure &cs, CCIndexDataAttr &indexAttr, int points);
  };

  using namespace fem; // FIXME when moved into FemLib
  class FemSetAnisoDirLines : public Process
  {
    public:
      FemSetAnisoDirLines(const Process &process) : Process(process) 
    {
      setName("Model/CCF/32 Set Ansio Dir From Lines");
      addParm("Element Type", "Type of the element", elementTypeToString(TRIANGLE_3), elementTypeList());
      addParm("Element Attribute", "Attribute to store nodal values", "Triangle Element", QStringList() << "Triangle Element" << "Wedge Element");
      addParm("Direction Type", "Direction in global coordinates for the special direction", dirTypeToString(DIRECTION_E2), dirTypeList());
      addParm("Projection Type", "Project the anisotropy direction as aligned parallel or orthogonal do the specified vector", 
          projTypeToString(PARALLEL), projTypeList());

      addParm("Tolerance", "Tolerance for the projection", "1e-6");
      addParm("Line CC", "Name of CC holding lines", "Make Line");
    }

      /// Set the anisotropy direction for selected elements
      bool run() 
      {
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw QString("%1::run No mesh").arg(name());

        QString ccName = mesh->ccName();
        if(ccName.isEmpty())
          throw QString("%1::run No cell complex selected").arg(name());

        QString lineName = parm("Line CC");
        if(lineName.isEmpty())
          throw QString("%1::run Line cell complex parameter empty").arg(name());

        double tolerance = parm("Tolerance").toDouble();
        if(tolerance <= 0)
          throw QString("%1::run Tolerance must be greater than 0").arg(name());

        int dirType = stringToDirType(parm("Direction Type"));
        if(dirType >= DirectionTypeCount)
          throw QString("%1::run Invalid direction %2").arg(name()).arg(parm("Direction Type"));

        int projType = stringToProjectionType(parm("Projection Type"));
        if(projType >= ProjectionTypeCount)
          throw QString("%1::run Invalid direction %2").arg(name()).arg(parm("Projection Type"));

        // Get the cell complex
        CCStructure &cs = mesh->ccStructure(ccName);
        CCStructure &csLine = mesh->ccStructure(lineName);
        return run(*mesh, cs, csLine, parm("Element Type"), parm("Element Attribute"), dirType, projType, tolerance);
      }
      bool run(Mesh &mesh, CCStructure &cs, CCStructure &csLine, const QString &elementType, const QString &elementName, 
          int dirType, int projType, double tolerance);

      /// Set the anisotopy direction for selected elements
      //
      // Clinton: Note this class is templated, so can work for different element types.
      //
      template <typename ElementAttr_T>
        void run(const CCStructure &cs, CCStructure &csLine, CCIndexDataAttr &indexAttr, ElementAttr_T &elementAttr, 
            int dirType, int projType, double tolerance)
        {
          int dim = ElementAttr_T::value_type::dimension();
          if(dim != 2 and dim != 3)
            throw QString("%1::run Invalid dimension of element %2").arg(name()).arg(dim);

          if(selectedFaces(cs, indexAttr).size() != 0 and  selectedVolumes(cs, indexAttr).size() != 0)
            throw QString("%1::run Both volumes and faces selected, please make up your mind").arg(name());

          if(dim == 2 and  selectedVolumes(cs, indexAttr).size() != 0)
            throw QString("%1::run mismatch between element type and selection, please make up your mind").arg(name());
          else if(dim == 3 and  selectedFaces(cs, indexAttr).size() != 0)
            throw QString("%1::run mismatch between element type and selection, please make up your mind").arg(name());

          CCIndexVec cells = (dim == 2 ? selectedFaces(cs, indexAttr) : selectedVolumes(cs, indexAttr));
          if (dim == 2 and cells.size() == 0) {
            for (CCIndex face : cs.cellsOfDimension(2))
              cells.push_back(face);
          } else if(dim == 3 and cells.size() == 0) {
            for (CCIndex vol : cs.cellsOfDimension(3))
              cells.push_back(vol);
          }

          // List of line segments
          typedef std::pair<Point3d, Point3d> Point3dPair;
          CCIndexVec edges = csLine.edges();
          std::vector<Point3dPair> lines(edges.size());
#pragma omp parallel for
          for(size_t i = 0; i < edges.size(); i++) {
            auto eb = csLine.edgeBounds(edges[i]);
            lines[i] = std::make_pair(indexAttr[eb.first].pos, indexAttr[eb.second].pos);
            mdxInfo << "Edge:" << indexAttr[eb.first].pos << "-" << indexAttr[eb.second].pos << endl;
          }

          // Loop through cells and assign direction
#pragma omp parallel for
          for(size_t i = 0; i < cells.size(); i++) {
            double minDist = std::numeric_limits<double>::max();
            size_t minIdx = 0;
            Point3d pos = indexAttr[cells[i]].pos;
            for(size_t j = 0; j < lines.size(); j++) {
              double d = distLinePoint(lines[j].first, lines[j].second, pos, true);
              if(minDist > d) {
                minDist = d;
                minIdx = j;
              }
            }
            Point3d dir = normalized(lines[minIdx].first - lines[minIdx].second);

            auto &cE = elementAttr[cells[i]];
            cE.setVertices(cs, indexAttr, cells[i]);
            if(dirType == DIRECTION_E2){
              cE.setE2Dir(indexAttr, dir, projType, tolerance);
              //if (cE.e2CosSin.x() == 0 and cE.e2CosSin.y() == 0)
              //  indexAttr[cE.v[0]].selected = indexAttr[cE.v[1]].selected = indexAttr[cE.v[0]].selected =true;
            }
            else if(dirType == DIRECTION_KPAR){
              cE.setKParDir(indexAttr, dir, projType, tolerance);
              //if (cE.kParCosSin.x() == 0 and cE.kParCosSin.y() == 0)
              //  indexAttr[cE.v[0]].selected = indexAttr[cE.v[1]].selected = indexAttr[cE.v[0]].selected =true;
            }
            else if(dirType == DIRECTION_KPER){
              cE.setKPerDir(indexAttr, dir, projType, tolerance);

            }
          }
          mdxInfo << QString("%1::run %2 elements updated").arg(name()).arg(cells.size()) << endl;
        }
  };

  class CalculateAreaRay : public Process
  {
    public:
      CalculateAreaRay(const Process &proc) : Process(proc)
    {
      setName("Model/CCF/101 Calculate Area Ray");
      setDesc("Area enclosed by a mesh.");
      addParm("Num Points", "Number of points to use in calculation","100");
    }
      bool run();

    private:
      Mesh *mesh = 0;
      QString ccName;
  };

  class StomaDims : public Process
  {
    public:
      StomaDims(const Process &proc) : Process(proc)
    {
      setName("Model/CCF/102 Stoma Dimensions");
      setDesc("Stoma Length and Width");
    }
      bool run();

    private:
      Mesh *mesh = 0;
      QString ccName;
  };

  class GeomCSV : public Process
  {
    public: 
      GeomCSV(const Process &proc) : Process(proc)
    {
      setName("Model/CCF/103 Geometry to CSV");
      setDesc("Print Geometry Dimensions to CSV File");
      addParm("Num Points", "Number of points to use in Area calculation","100");
      addParm("GC Pressure", "Guard Cell Pressure","");
      addParm("SC Pressure", "Subsidiary Cell Pressure","");
      addParm("File Name", "Name of CSV File","default.csv");
      addParm("Write Header", "Write header to the CSV file?","False");
    }
      bool run();

    private:
      Mesh *mesh = 0;
      QString ccName;
  };

  class HeatMapPrint : public Process
  {
    public:
      HeatMapPrint(const Process &proc) : Process(proc)
    {
      setName("Model/CCF/103 Heat Map Display");
      setDesc("Print heat map to terminal");
      addParm("Heat Name", "Name of heat map to save, empty for current, or All", "");
      addParm("Labeling", "Labeling to use, empty for current", "");
      addParm("Use Selection","Export selected cells only?\n"
          "'No' exports the whole attribute, whether the cells exist or not","No", booleanChoice());
      addParm("Export Labelings","If Labeling is 'Labels', export the other labelins as well?\n","No", booleanChoice());
    }
      bool processParms()
      {
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw QString("%1::processParms No current mesh").arg(name());

        QString heatName = parm("Heat Name");
        if(heatName.isEmpty())
          heatName = mesh->heat();
        if(heatName.isEmpty())
          throw QString("%1::processParms No heat map specified").arg(name());
        if(heatName != "All" and !mesh->heatExists(heatName))
          throw QString("%1::run No heat map named (%2) exists").arg(name()).arg(heatName);
        setParm("Heat Name", heatName);

        QString labeling = parm("Labeling");
        if(labeling.isEmpty())
          labeling = mesh->labeling();
        if(labeling.isEmpty())
          throw QString("%1::processParms No labeling specified").arg(name());
        if(!mesh->labelingExists(labeling))
          throw QString("%1::processParms No labeling named (%2) exists").arg(name()).arg(labeling);
        setParm("Labeling", labeling);

        if(labeling != "Labels" and stringToBool("Export Labelings"))
          throw QString("%1::processParms Labeling must be 'Labels' to export labelings").arg(name());
        return true;
      } 

      bool run()
      {
        Mesh *mesh = currentMesh();
        if(!mesh)
          throw QString("%1::run No current mesh").arg(name());

        processParms();
        return run(*mesh, parm("Heat Name"), parm("Labeling"), stringToBool(parm("Use Selection")), stringToBool(parm("Export Labelings")));
      }

      bool run(Mesh &mesh, const QString &heatName, const QString &labeling, bool useSelection, bool exportLabelings);
  };
}
#endif
