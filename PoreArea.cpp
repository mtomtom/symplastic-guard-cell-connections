//
// This file is part of MorphoDynamX - http://www.MorphoDynamX.org
// Copyright (C) 2012-2016 Richard S. Smith and collaborators.
//
// If you use MorphoDynamX in your work, please cite:
//   http://dx.doi.org/10.7554/eLife.05864
//
// MorphoDynamX is free software, and is licensed under under the terms of the 
// GNU General (GPL) Public License version 2.0, http://www.gnu.org/licenses.
//
#include <FemMembranes.hpp>

namespace mdx
{
  double sign(const Point2d &p1, const Point2d &p2, const Point2d &p3)
  {
    return (p1.x() - p3.x()) * (p2.y() - p3.y()) - (p2.x() - p3.x()) * (p1.y() - p3.y());
  }

  bool pointInTriangle(const Point2d &pt, const Point2d v1, const Point2d v2, const Point2d &v3)
  {
    double d1 = sign(pt, v1, v2);
    double d2 = sign(pt, v2, v3);
    double d3 = sign(pt, v3, v1);

    bool hasNeg = (d1 < 0) or (d2 < 0) or (d3 < 0);
    bool hasPos = (d1 > 0) or (d2 > 0) or (d3 > 0);

    return !hasNeg or !hasPos;
  }

  // Raster to world coordinates
  Point2d toWorld(const Point2i &p, Point2d &origin, double step)
  {
    return origin + Point2d(p) * step;
  }

  // World to raster coordinates
  Point2d toRaster(const Point2d &p, Point2d &origin, double step)
  {
    return (p - origin)/step;
  }

  // Calculate the offset into the raster array
  size_t offset(const Point2i &p, const Point2i &size)
  {
    return p.y() * size.x() + p.x();
  }

  // Make a qglviewer quaternion from a direction and an angle
  qglviewer::Quaternion qtrn(const Point3d &p, double angle)
  {
    Quaternion rotate(p, angle);
    return qglviewer::Quaternion(rotate[0], rotate[1], rotate[2], rotate[3]);
  }

  double PoreArea::run(CCStructure &cs, CCIndexDataAttr &indexAttr, qglviewer::Frame *frame, double step)
  {
    BoundingBox2d bb;
    CCIndexPoint3dAttr pos;
    for(CCIndex v : cs.vertices()) {
      auto p = indexAttr[v].pos;
      if(frame)
        p = Point3d(frame->inverseCoordinatesOf(qglviewer::Vec(p)));
      bb |= Point2d(p);
      pos[v] = p;
    }

    Point2i size = (bb.pmax() - bb.pmin())/step + Point2d(2,2);
    std::vector<char> raster(size.x() * size.y());
    for(uint i = 0; i < raster.size(); i++)
      raster[i] = 0;
    Point2d origin = bb.pmin() - Point2d(step, step); // Pad to give border
    int inCount = 0;
    for(CCIndex f : cs.faces()) {
      BoundingBox2d tb;
      auto fV = cs.faceVertices(f);
      if(fV.size() != 3) {
        mdxInfo << "Found non triangle face:" << f << endl;
        continue;
      }
      Point2dVec tri(3);
      for(int i = 0; i < 3; i++) {
        tri[i] = Point2d(pos[fV[i]]);
        tb |= tri[i];
      }
      Point2i imin = toRaster(tb.pmin(), origin, step);
      Point2i imax = toRaster(tb.pmax(), origin, step);
      for(int y = imin.y(); y <= imax.y(); y++) 
        for(int x = imin.x(); x <= imax.x(); x++)
          if(pointInTriangle(toWorld(Point2i(x, y), origin, step), tri[0], tri[1], tri[2])) {
            inCount++;
            raster[offset(Point2d(x, y), size)] = 1;
          }
    }
    int count = 0;
    for(uint i = 0; i < raster.size(); i++)
      if(raster[i] == 0)
        count++;

    // Find connected region from border
    std::set<Point2i> tryP;
    for(int i = 0; i < size.x(); i++) {
      tryP.insert(Point2i(i, 0));
      tryP.insert(Point2i(i, size.y() - 1));
    }
    for(int i = 0; i < size.y(); i++) {
      tryP.insert(Point2i(0, i));
      tryP.insert(Point2i(size.x() - 1, i));
    }
    while(tryP.size() > 0) {
      std::set<Point2i> newP;
      for(auto &p : tryP)
        if(raster[offset(p, size)] == 0) {
          raster[offset(p, size)] = 2;

          // Check neighbors
          if(p.x() > 0)
            newP.insert(Point2i(p.x() - 1, p.y()));
          if(p.x() < size.x() - 1)
            newP.insert(Point2i(p.x() + 1, p.y()));
          if(p.y() > 0)
            newP.insert(Point2i(p.x(), p.y() - 1));
          if(p.y() < size.y() - 1)
            newP.insert(Point2i(p.x(), p.y() + 1));
        }
      tryP = newP;
    }
    int porePix = 0;
    for(uint i = 0; i < raster.size(); i++)
      if(raster[i] == 0)
        porePix++;

    return porePix * step * step;
  }

  bool PoreArea::run()
  {
    Stack *stack = currentStack();
    if(!stack)
      throw QString("%1::run No current stack").arg(name());

    Mesh *mesh = currentMesh();
    if(!mesh)
      throw QString("%1::run No current mesh");

    QString ccName = mesh->ccName();
    if(ccName.isEmpty())
      throw QString("%1::run No current cell complex");

    double step = parm("Step").toDouble();
    if(step <= 0)
      throw QString("%1::run Step must be > 0");

    double qStep = parm("Q Step").toDouble();
    if(qStep <= 0)
      throw QString("%1::run Q Step must be > 0");
    double initQStep = qStep;

    bool transform = stringToBool(parm("Transform"));
    bool rotateZ = stringToBool(parm("Rotate Z"));

    qglviewer::Frame *frame = 0;
    if(transform)
      frame = &stack->frame();

    auto &cs = mesh->ccStructure(ccName);
    auto &indexAttr = mesh->indexAttr();

    double area = run(cs, indexAttr, frame, step);
    if(transform) {
      bool improved = true;
      while(improved) {
        if(!progressAdvance(1))
          userCancel();
        improved = false;

        frame->rotate(qtrn(Point3d(1.0, 0, 0), qStep));
        double newArea = run(cs, indexAttr, frame, step);
        if(newArea > area) {
          area = newArea;
          improved = true;
        } else {
          frame->rotate(qtrn(Point3d(1.0, 0, 0), -2.0 * qStep));
          newArea = run(cs, indexAttr, frame, step);
          if(newArea > area) {
            area = newArea;
            improved = true;
          } else
            frame->rotate(qtrn(Point3d(1.0, 0, 0), qStep));
        }
        frame->rotate(qtrn(Point3d(0, 1.0, 0), qStep));
        newArea = run(cs, indexAttr, frame, step);
        if(newArea > area) {
          area = newArea;
          improved = true;
        } else {
          frame->rotate(qtrn(Point3d(0, 1.0, 0), -2.0 * qStep));
          newArea = run(cs, indexAttr, frame, step);
          if(newArea > area) {
            area = newArea;
            improved = true;
          } else
            frame->rotate(qtrn(Point3d(0, 1.0, 0), qStep));
        }
        if(rotateZ) {
          frame->rotate(qtrn(Point3d(0, 0, 1.0), qStep));
          newArea = run(cs, indexAttr, frame, step);
          if(newArea > area) {
            area = newArea;
            improved = true;
          } else {
            frame->rotate(qtrn(Point3d(0, 0, 1.0), -2.0 * qStep));
            newArea = run(cs, indexAttr, frame, step);
            if(newArea > area) {
              area = newArea;
              improved = true;
            } else
              frame->rotate(qtrn(Point3d(0, 0, 1.0), qStep));
          }
        }

        // Reduce step size and repeat
        if(!improved and qStep > initQStep/100) {
          qStep /= 10.0;
          improved = true;
        }
        updateState();
        updateViewer();
        mdxInfo << "Area:" << area << endl;

        AttrMap<QString,double> &values = mesh->attributes().attrMap<QString,double>("savedConstants");
        values["Area"] = area;
      }
    }

    return true;
  }

  REGISTER_PROCESS(PoreArea);
}