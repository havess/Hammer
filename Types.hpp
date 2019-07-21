#pragma once
#include <SFML/Graphics.hpp>
#include <glm/glm.hpp>
#include <iostream>
#include <memory>
#include <queue>
#include <set>

#include "DCEL.hpp"

template <typename T>
using sptr = std::shared_ptr<T>;
using vec3 = glm::vec3;
using vec2 = glm::vec2;
using mat3 = glm::mat3;
using mat2 = glm::mat2;

/****** CONSTANTS *******/
constexpr int width = 2000, height = 1400;
constexpr int surfaceWidth = 1600, surfaceHeight = 1000;
constexpr int delta = 200, cellSize = 80;

struct Triangle;
struct Point {
  Point(float x, float y, DCEL::VertexPtr vptr = nullptr);

  Point() : Point(0, 0, nullptr) {}

  void addCompressiveForce(const vec2 &vec);

  void addTensileForce(const vec2 &vec);

  void computeSeperationTensor();

  void clearForces() {
    tensileForces.clear();
    compressiveForces.clear();
  }
  void drawForces(sf::RenderWindow &window);

  vec2 getLargestEigenRep();

  Point operator+(const Point &p);

  /****** DATA ******/
  float x, y;
  int index;
  bool isDead = false;
  vec2 mSepEVec;
  float mSepEVal = 0;
  std::vector<vec2> tensileForces;
  std::vector<vec2> compressiveForces;
  mat2 seperationTensor;
  mat2 residualTensor = mat2(0);
  // Purely for HalfEdge data structure compatibility
  DCEL::VertexPtr vertex;
};

struct HEdge {
  sptr<HEdge> pair, next;
  sptr<Point> to, from;
  sptr<Triangle> face;
};

struct Triangle {
  Triangle(sptr<Point> p0, sptr<Point> p1, sptr<Point> p2);

  void addDirectionalStress(const vec2 &direction,
                            std::function<float(const vec2 &)> stressFn);

  float area() const;

  void computeForceOnNodes(sptr<Point> pt = nullptr);

  void computeGreensStrainDerivative();

  void relax();

  void createShape();

  void draw(sf::RenderWindow &window);

  vec2 getCentroid();

  vec2 getLargestEigenRep();

  vec2 worldToLocal(const vec2 &vec);

  vec3 worldToBarycentric(const vec2 &vec);

  void replacePoint(sptr<Point> old, sptr<Point> newP);

  /****** DATA ******/
  sptr<HEdge> e01, e12, e20;
  sptr<Point> p0, p1, p2;
  sf::ConvexShape shape;
  vec2 u, v;
  mat2 worldToLocalMat;
  mat2 stressTensor;
  mat3 barycentricBasisMat;
};

class Surface {
 public:
  Surface(const std::set<sptr<Triangle>> &triangles,
          const std::vector<sptr<Point>> &points);

  void computeSeperationTensors();
  void initFractureState();
  void FractureStep(sf::RenderWindow &window);
  void Relax();
  void draw(sf::RenderWindow &window);
  void PrintDataToLog();

  /****** DATA ******/
  std::set<sptr<Triangle>> mTriangles;
  std::vector<sptr<Point>> mPoints;
  std::vector<std::pair<sptr<Point>, sptr<Point>>> mCracks;

  struct SeperationComp {
    bool operator()(const sptr<Point> &a, const sptr<Point> &b) {
      return a->mSepEVal <= b->mSepEVal;
    }
  };
  std::priority_queue<sptr<Point>, std::vector<sptr<Point>>, SeperationComp>
      mSeperationPQ;

  static constexpr float snapThreshold = 13.0f;
};

