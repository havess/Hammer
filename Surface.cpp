#define GLM_ENABLE_EXPERIMENTAL
#include <fstream>
#include <glm/gtx/string_cast.hpp>

#include "Types.hpp"
#include "helpers.hpp"

Point::Point(float x, float y, DCEL::VertexPtr vert)
    : x(x), y(y), vertex(vert) {
  seperationTensor = mat2(0);
}

Point Point::operator+(const Point &p) { return Point(x + p.x, y + p.y); }

void Point::addTensileForce(const vec2 &vec) { tensileForces.push_back(vec); }

void Point::addCompressiveForce(const vec2 &vec) {
  compressiveForces.push_back(vec);
}

void Point::computeSeperationTensor() {
  seperationTensor = mat2(0);
  mSepEVec = vec2(0);
  mSepEVal = 0;
  if (isDead) return;
  if (tensileForces.empty() && compressiveForces.empty()) return;
  auto edge = vertex->edge;
  if (!edge->twin && !edge->prev->twin && !edge->next->twin) {
    return;
  }
  // std::cout << std::endl;
  vec2 fTensile(0), fComp(0);
  for (const auto &v : tensileForces) {
    // std::cout << "TENSILE " << glm::to_string(v) << std::endl;
    fTensile += v;
  }
  for (const auto &v : compressiveForces) fComp += v;
  mat2 mTensile(0), mComp(0);
  for (const auto &v : tensileForces) mTensile += m(v);
  for (const auto &v : compressiveForces) mComp += m(v);
  /*std::cout << "COMPUTING SEPERATION TENSOR " << glm::to_string(fTensile)
            << std::endl
            << glm::to_string(mTensile) << std::endl
            << glm::to_string(fComp) << std::endl
            << glm::to_string(mComp) << std::endl;
  std::cout << "OLD SEPERATION TENSOR for " << index << std::endl
            << glm::to_string(seperationTensor) << std::endl;*/
  seperationTensor = (-m(fTensile) + mTensile + m(fComp) - mComp) / 2.0f;
  seperationTensor += residualTensor;

  // std::cout << "NEW SEPERATION TENSOR for " << index << std::endl
  //          << glm::to_string(seperationTensor) << std::endl;
  getLargestEigenRep();
  // std::cout << index << " " << mSepEVal << std::endl;
}

vec2 Point::getLargestEigenRep() {
  if (tensileForces.empty() && compressiveForces.empty()) return vec2(0, 0);
  vec2 v1, v2;
  float l1, l2;
  eigenDecomp(seperationTensor, v1, v2, l1, l2);
  // auto ret = vec2(0);
  mSepEVal = 0;
  mSepEVec = vec2(0);
  if (l1 > 0) {
    mSepEVal = l1;
    mSepEVec = v1;
    // ret = l1 * v1;
  }
  if (l2 > 0 && l2 > l1) {
    mSepEVal = l2;
    mSepEVec = v2;
    // ret = l2 * v2;
  }
  // mSepEVec = ret;
  return mSepEVec;
}

void Point::drawForces(sf::RenderWindow &window) {
  for (auto &v : tensileForces) {
    sf::RectangleShape line(sf::Vector2f(0.1f * glm::length(v), 4));
    line.setPosition(x, y);
    line.setFillColor(sf::Color(0, 240, 0));
    auto angle = vangle(v);
    line.rotate(angle);
    // window.draw(line);
  }
  for (auto &v : compressiveForces) {
    sf::RectangleShape line(sf::Vector2f(10, 4));
    line.setPosition(x, y);
    line.setFillColor(sf::Color(240, 0, 0));
    auto angle = vangle(v);
    line.rotate(angle);
    // window.draw(line);
  }
}

Triangle::Triangle(sptr<Point> p0, sptr<Point> p1, sptr<Point> p2)
    : p0(p0), p1(p1), p2(p2), stressTensor(0) {
  vec3 a(p0->x, p0->y, 0);
  vec3 b(p1->x, p1->y, 0);
  vec3 c(p2->x, p2->y, 0);
  vec3 tempU = glm::normalize(b - a);
  auto n = glm::normalize(glm::cross(tempU, c - a));
  v = vec2(glm::cross(n, tempU));
  u = vec2(tempU);
  worldToLocalMat = glm::inverse(mat2(u, v));
  vec2 m1(0);
  vec2 m2 = worldToLocal(vec2(b - a));
  vec2 m3 = worldToLocal(vec2(c - a));
  barycentricBasisMat =
      glm::inverse(mat3(vec3(m1, 1), vec3(m2, 1), vec3(m3, 1)));
}

void Triangle::createShape() {
  shape.setPointCount(3);
  shape.setPoint(0, sf::Vector2f(p0->x, p0->y));
  shape.setPoint(1, sf::Vector2f(p1->x, p1->y));
  shape.setPoint(2, sf::Vector2f(p2->x, p2->y));
  shape.setOutlineThickness(1);
  shape.setFillColor(sf::Color::Transparent);
  shape.setOutlineColor(sf::Color(20, 20, 20));
}

void Triangle::draw(sf::RenderWindow &window) { window.draw(shape); }

vec2 Triangle::getCentroid() {
  vec2 v0(p0->x, p0->y);
  vec2 v1(p1->x, p1->y);
  vec2 v2(p2->x, p2->y);
  return (v0 + v1 + v2) / 3.0f;
}

void Triangle::addDirectionalStress(
    const vec2 &direction, std::function<float(const vec2 &)> stressFn) {
  auto d = worldToLocal(direction);
  stressTensor += stressFn(getCentroid()) * glm::outerProduct(d, d);
}

vec2 Triangle::getLargestEigenRep() {
  vec2 v1, v2;
  float l1, l2;
  eigenDecomp(stressTensor, v1, v2, l1, l2);
  auto ret = vec2(0);
  if (l1 > 0) ret = l1 * v1;
  if (l2 > 0 && l2 > l1) ret = l2 * v2;
  return ret.x * u + ret.y * v;
}

float Triangle::area() const {
  vec3 a(p0->x, p0->y, 0);
  vec3 b(p1->x, p1->y, 0);
  vec3 c(p2->x, p2->y, 0);
  return glm::length(glm::cross(b - a, c - a)) / 2.0f;
}

void Triangle::replacePoint(sptr<Point> old, sptr<Point> newP) {
  if (p0 == old)
    p0 = newP;
  else if (p1 == old)
    p1 = newP;
  else if (p2 == old)
    p2 = newP;
}

void Triangle::computeForceOnNodes(sptr<Point> point) {
  if (stressTensor == mat2(0)) return;
  // std::cout << glm::to_string(stressTensor) << std::endl;
  if (point && point != p0 && point != p1 && point != p2) return;

  vec2 a(p0->x, p0->y);
  vec2 b(p1->x, p1->y);
  vec2 c(p2->x, p2->y);

  vec2 v1, v2;
  float l1, l2;
  eigenDecomp(stressTensor, v1, v2, l1, l2);
  v1 = glm::normalize(v1);
  v2 = glm::normalize(v2);
  // Get tensile tensor
  auto tensileTensor = std::max(0.0f, l1) * m(v1) + std::max(0.0f, l2) * m(v2);

  // Get compressive tensor
  auto compressiveTensor =
      std::min(0.0f, l1) * m(v1) + std::min(0.0f, l2) * m(v2);

  auto A = area();

  for (int i = 0; i < 3; i++) {
    auto p = point ? point : (i == 0 ? p0 : (i == 1 ? p1 : p2));
    // Get tensile force
    vec2 tensileForce(0);
    for (int j = 0; j < 3; j++) {
      float BjlBikSigkl = 0.0f;
      for (int k = 0; k < 2; k++) {
        // !!!!!!CAREFUL!!!!!!! GLM uses column major rep, so
        // we need to flip the indices from the force equation
        // otherwise youll get *super* wrong stuff.
        for (int l = 0; l < 2; l++) {
          BjlBikSigkl += barycentricBasisMat[l][j] * barycentricBasisMat[k][i] *
                         tensileTensor[l][k];
        }
      }
      if (j == 0)
        tensileForce += a * BjlBikSigkl;
      else if (j == 1)
        tensileForce += b * BjlBikSigkl;
      else
        tensileForce += c * BjlBikSigkl;
    }
    tensileForce *= -A;
    // std::endl;

    // Get compressive force
    vec2 compressiveForce(0);
    for (int j = 0; j < 3; j++) {
      float BjlBikSigkl = 0.0f;
      for (int k = 0; k < 2; k++) {
        for (int l = 0; l < 2; l++) {
          BjlBikSigkl += barycentricBasisMat[l][j] * barycentricBasisMat[k][i] *
                         compressiveTensor[l][k];
        }
      }
      if (j == 0)
        compressiveForce += a * BjlBikSigkl;
      else if (j == 1)
        compressiveForce += b * BjlBikSigkl;
      else
        compressiveForce += c * BjlBikSigkl;
    }
    compressiveForce *= -A;
    if (tensileForce != vec2(0)) p->addTensileForce(tensileForce);
    if (compressiveForce != vec2(0)) p->addCompressiveForce(compressiveForce);
  }
}

void Triangle::computeGreensStrainDerivative() {}

void Triangle::relax() {
  static constexpr float timestep = 0.03f;
  for (int n = 0; n < 3; n++) {
    auto pn = n == 0 ? p0 : (n == 1 ? p1 : p2);
    vec3 totalForce(0);
    for (auto &f : pn->tensileForces) {
      totalForce += vec3(f, 0);
    }
    for (auto &f : pn->compressiveForces) {
      totalForce += vec3(f, 0);
    }
    vec3 deltaPn = timestep * totalForce;

    // Derivative of Green's strain tensor with respect to node positions.
    // Note indices for this are backwards [depth][column][row]
    mat3 dedpn[2];
    for (int r = 0; r < 2; r++) {
      dedpn[r] = mat3(0);
      for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
          float firstTerm = 0, secondTerm = 0;
          for (int m = 0; m < 3; m++) {
            auto pointm = m == 0 ? p0 : (m == 1 ? p1 : p2);
            auto pm = vec2(pointm->x, pointm->y);
            firstTerm +=
                pm[r] * barycentricBasisMat[i][n] * barycentricBasisMat[j][m];
            secondTerm +=
                pm[r] * barycentricBasisMat[i][m] * barycentricBasisMat[j][n];
          }
          dedpn[r][j][i] = 0.5f * (firstTerm + secondTerm);
        }
      }
    }

    mat3 deltaStressTensor(0);
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        for (int r = 0; r < 2; r++) {
          deltaStressTensor[i][j] += dedpn[r][i][j] * deltaPn[r];
        }
      }
    }
    stressTensor += mat2(deltaStressTensor);
  }
}

vec2 Triangle::worldToLocal(const vec2 &vec) { return worldToLocalMat * vec; }
vec3 Triangle::worldToBarycentric(const vec2 &vec) {
  vec3 v3(worldToLocal(vec - vec2(p0->x, p0->y)), 1);
  auto ret = barycentricBasisMat * v3;
  return ret;
}

Surface::Surface(const std::set<sptr<Triangle>> &triangles,
                 const std::vector<sptr<Point>> &points)
    : mTriangles(triangles), mPoints(points) {
  auto linearStressFn = [](const vec2 &v) {
    auto diff = std::abs(v.x - (delta + surfaceWidth / 2.0f));
    return 2.0f * std::clamp(1.0f - diff / 300.0f, 0.0f, 1.0f);
  };

  auto radialStressFn = [](const vec2 &v) {
    auto diff = v - vec2(delta + 100, delta + 100);
    auto dist = std::min(float(M_PI * 150.0f), glm::length(diff));
    auto wi = 2.0f * std::cos(dist / 300.0f);
    return wi;
  };

  /************** CREATE SHAPE, ADD STRESS, ADD FORCE TO NODES **********/
  for (auto &t : mTriangles) {
    t->createShape();
    auto centroid = t->getCentroid();
    auto dir = glm::normalize(centroid - vec2(delta + 100, delta + 100));
    t->addDirectionalStress(vec2(1, 0), linearStressFn);
    vec2 perp(dir.y, -dir.x);
    // t->addDirectionalStress(perp, [radialStressFn](const vec2 &v) {
    //   return 1.25f * radialStressFn(v);
    // });
    t->computeForceOnNodes();
  }
}

void Surface::initFractureState() {
  mSeperationPQ = std::priority_queue<sptr<Point>, std::vector<sptr<Point>>,
                                      SeperationComp>();
  for (auto &p : mPoints) mSeperationPQ.push(p);
}

bool getIntersection(const vec2 &normal, const vec2 &p, DCEL::HalfEdgePtr edge,
                     vec2 &inter) {
  auto p0 = edge->prev->vertex->point;
  vec2 a(p0.x, p0.y);
  auto p1 = edge->vertex->point;
  vec2 b(p1.x, p1.y);
  bool aboveA = glm::dot(a - p, normal) > 0;
  bool aboveB = glm::dot(b - p, normal) > 0;
  if ((aboveA && aboveB) || (!aboveA && !aboveB)) return false;
  vec3 p2(a, 0);
  vec3 p4(b, 0);
  vec3 p3(p, 0);
  vec3 n(normal, 0);
  float prod1 = glm::dot(p2 - p3, n);
  float prod2 = glm::dot(p4 - p2, n);
  float prod3 = prod1 / prod2;
  auto ret = p2 - (p4 - p2) * prod3;
  inter.x = ret.x;
  inter.y = ret.y;
  return true;
}

void Surface::FractureStep(sf::RenderWindow &window) {
  auto isBoundaryPoint = [](const DCEL::VertexPtr &v) {
    auto boundRunner = v->edge;
    bool first = true;
    while (boundRunner && (boundRunner != v->edge || first)) {
      first = false;
      boundRunner = boundRunner->vertexNextCW();
    }
    return boundRunner == nullptr;
  };

  // Get point with highest seperation tensor
  auto top = mSeperationPQ.top();
  std::cout << "FORCES " << top->tensileForces.size() << std::endl;
  std::cout << "STEP " << top->index << " " << top->mSepEVal << std::endl;
  top->residualTensor = mat2(0);
  auto topEigenVector = top->mSepEVec;
  auto topEigenValue = top->mSepEVal;
  mSeperationPQ.pop();
  auto p = top->vertex->point;
  /****** WORKAROUND FOR INSTABILITY ******/
  static int lastIndex = top->index;
  static int count = 0;
  static float lastVal = topEigenValue;
  if (lastIndex != top->index)
    count = 0;
  else if (std::abs(lastVal - topEigenValue) < 0.1)
    count++;

  if (topEigenValue > 1e6) {
    top->isDead = true;
    return;
  }
  lastIndex = top->index;
  lastVal = topEigenValue;

  if (count == 15) {
    top->isDead = true;
  }

  auto normal = glm::normalize(vec3(top->mSepEVec, 0));

  // Find number of free boundaries.
  assert(top->vertex->edge && "EDGE CANT BE NULL");
  auto runner = top->vertex->edge;
  bool first = true;
  bool firstIntersection = false;
  while (runner && (runner != top->vertex->edge || first)) {
    first = false;
    auto next = runner->vertexNextCW();
    if (!next) {
      firstIntersection = true;
      break;
    }
    runner = next;
  }

  std::vector<sptr<Point>> addResidualTo;
  while (runner) {
    auto next = runner->twin ? runner->vertexNextCCW() : nullptr;
    auto backwardsRunner = runner->vertexNextCW();
    // Check if the current triangle 'runner' belongs to
    // intersects with the plane.
    vec2 inter;
    if (getIntersection(normal, vec2(p.x, p.y), runner->prev, inter)) {
      // Check if we should snap.
      vec2 crack(inter.x - p.x, inter.y - p.y);
      auto p1 = runner->prev->vertex->point;
      auto p2 = runner->next->vertex->point;
      vec2 e1 = vec2(p1.x - p.x, p1.y - p.y);
      vec2 e2 = vec2(p2.x - p.x, p2.y - p.y);

      if (std::abs(vangle(crack) - vangle(e1)) < snapThreshold) {
        inter = vec2(p1.x, p1.y);
        // Check if it is a free boundary edge.
        if (runner->twin) {
          mCracks.push_back(
              std::make_pair(mPoints[runner->prev->vertex->pointIndex],
                             mPoints[runner->vertex->pointIndex]));
          if (!isBoundaryPoint(runner->prev->vertex))
            addResidualTo.push_back(mPoints[runner->prev->vertex->pointIndex]);
          else {
            std::cout << "CASE 1" << std::endl;
            // If we hit a boundary point we need to split this point into 2
            // now.
            auto newVertex = std::make_shared<DCEL::Vertex>(p1, mPoints.size());
            mPoints.push_back(std::make_shared<Point>(p1.x, p1.y, newVertex));
            mPoints.back()->index = mPoints.size() - 1;
            if (runner->twin) {
              auto fixRunner = runner->twin;
              while (fixRunner) {
                fixRunner->face->replacePoint(
                    mPoints[fixRunner->vertex->pointIndex], mPoints.back());
                fixRunner->vertex = newVertex;
                newVertex->edge = fixRunner;
                fixRunner = fixRunner->vertexNextCW();
              }
            }
          }
          backwardsRunner = runner;
          runner->twin->twin = nullptr;
          runner->twin = nullptr;

          //
        } else {
          runner = next;
          std::cout << "CONTINUE" << std::endl;
          continue;
        }
      } else if (std::abs(vangle(crack) - vangle(e2)) < snapThreshold) {
        inter = vec2(p2.x, p2.y);
        // Check if it is a free boundary edge.
        if (runner->next->twin) {
          mCracks.push_back(
              std::make_pair(mPoints[runner->next->vertex->pointIndex],
                             mPoints[runner->vertex->pointIndex]));
          if (!isBoundaryPoint(runner->next->vertex))
            addResidualTo.push_back(mPoints[runner->next->vertex->pointIndex]);
          else {
            std::cout << "CASE 2" << std::endl;
            auto newVertex = std::make_shared<DCEL::Vertex>(p2, mPoints.size());
            mPoints.push_back(std::make_shared<Point>(p2.x, p2.y, newVertex));
            mPoints.back()->index = mPoints.size() - 1;
            if (runner->next->twin) {
              auto fixRunner = runner->next->twin->prev;
              while (fixRunner) {
                fixRunner->face->replacePoint(
                    mPoints[fixRunner->vertex->pointIndex], mPoints.back());
                fixRunner->vertex = newVertex;
                newVertex->edge = fixRunner;
                if (!fixRunner->twin) break;
                fixRunner = fixRunner->vertexNextCCW();
              }
            }
          }
          runner->next->twin->twin = nullptr;
          runner->next->twin = nullptr;
        } else {
          std::cout << "CONTINUE" << std::endl;
          runner = next;
          continue;
        }
      } else {
        // Create new point at the tip of the crack edge.
        auto newVertex = std::make_shared<DCEL::Vertex>(
            Point2D(inter.x, inter.y), mPoints.size());
        mPoints.push_back(std::make_shared<Point>(inter.x, inter.y, newVertex));
        mPoints.back()->index = mPoints.size() - 1;

        // Get 4 points belonging to new split triangle.
        auto t0 = mPoints.back();
        auto t1 = mPoints[runner->prev->vertex->pointIndex];
        auto t2 = mPoints[runner->vertex->pointIndex];
        auto t3 = mPoints[runner->next->vertex->pointIndex];

        // Add first triangle.
        auto tri0 = std::make_shared<Triangle>(t0, t1, t2);
        tri0->createShape();
        mTriangles.insert(tri0);
        auto angle = vangle(tri0->u) - vangle(runner->face->u);
        auto rot = rotate2(angle);
        tri0->stressTensor =
            rot * runner->face->stressTensor * glm::transpose(rot);

        // Add second triangle.
        auto tri1 = std::make_shared<Triangle>(t0, t2, t3);
        tri1->createShape();
        mTriangles.insert(tri1);
        angle = vangle(tri1->u) - vangle(runner->face->u);
        rot = rotate2(angle);
        tri1->stressTensor =
            rot * runner->face->stressTensor * glm::transpose(rot);
        mTriangles.erase(runner->face);

        // Create new halfedges, we dont care what the new face indexes
        // are since we arent going to relax the surface.
        // The variable names here correspond to a drawing I had, this is
        // going to look like nonsense without the drawing and my ASCII
        // art game is weak.
        auto a = std::make_shared<DCEL::HalfEdge>(-1, -1);
        auto b = std::make_shared<DCEL::HalfEdge>(-1, -1);
        auto d = std::make_shared<DCEL::HalfEdge>(-1, -1);
        auto e = std::make_shared<DCEL::HalfEdge>(-1, -1);
        auto f = runner->prev->twin;
        auto g = f ? std::make_shared<DCEL::HalfEdge>(-1, -1) : nullptr;
        auto h = f ? std::make_shared<DCEL::HalfEdge>(-1, -1) : nullptr;
        auto k = f ? std::make_shared<DCEL::HalfEdge>(-1, -1) : nullptr;
        auto l = f ? std::make_shared<DCEL::HalfEdge>(-1, -1) : nullptr;
        if (f) addResidualTo.push_back(mPoints.back());

        // a is a free boundary
        a->twin = nullptr;
        a->prev = e;
        a->next = runner->next;
        a->vertex = runner->vertex;
        a->face = tri1;

        // b is a free boundary
        b->twin = nullptr;
        b->prev = runner;
        b->next = d;
        b->vertex = newVertex;
        b->face = tri0;
        newVertex->edge = b;

        d->twin = g;
        d->prev = b;
        d->next = runner;
        d->vertex = runner->prev->vertex;
        d->face = tri0;

        e->twin = h;
        e->prev = runner->next;
        e->next = a;
        e->vertex = newVertex;
        e->face = tri0;

        runner->next->prev = a;
        runner->next->next = e;
        runner->next->face = tri1;

        runner->prev->vertex->edge = d;
        runner->prev = d;
        runner->next = b;
        mTriangles.erase(runner->face);
        runner->face = tri0;

        mCracks.push_back(std::make_pair(mPoints[b->vertex->pointIndex],
                                         mPoints[a->vertex->pointIndex]));

        // If the split edge is not a free boundary fix the other adjacent
        // triangle.
        if (f) {
          // Get 4th point for triangle we are fixing.
          auto t4 = mPoints[f->next->vertex->pointIndex];

          // Create split triangles
          auto tri2 = std::make_shared<Triangle>(t0, t4, t1);
          tri2->createShape();
          mTriangles.insert(tri2);
          angle = vangle(tri2->u) - vangle(f->face->u);
          rot = rotate2(angle);
          tri2->stressTensor =
              rot * f->face->stressTensor * glm::transpose(rot);

          auto tri3 = std::make_shared<Triangle>(t0, t3, t4);
          tri3->createShape();
          mTriangles.insert(tri3);
          angle = vangle(tri3->u) - vangle(f->face->u);
          rot = rotate2(angle);
          tri3->stressTensor =
              rot * f->face->stressTensor * glm::transpose(rot);

          g->twin = d;
          g->prev = f->prev;
          g->next = k;
          g->vertex = newVertex;
          g->face = tri2;

          h->twin = e;
          h->prev = l;
          h->next = f->next;
          h->vertex = f->vertex;
          h->face = tri3;
          f->vertex->edge = h;

          k->twin = l;
          k->prev = g;
          k->next = f->prev;
          k->vertex = f->next->vertex;
          k->face = tri2;

          l->twin = k;
          l->prev = f->next;
          l->next = h;
          l->vertex = newVertex;
          l->face = tri3;

          f->next->prev = h;
          f->next->next = l;

          f->prev->prev = k;
          f->prev->next = g;

          mTriangles.erase(f->face);
        } else {
          auto newBoundaryVertex = std::make_shared<DCEL::Vertex>(
              Point2D(inter.x, inter.y), mPoints.size());
          mPoints.push_back(
              std::make_shared<Point>(inter.x, inter.y, newVertex));
          mPoints.back()->index = mPoints.size() - 1;
          e->vertex = newBoundaryVertex;
          e->face->replacePoint(mPoints[mPoints.size() - 2], mPoints.back());
          newBoundaryVertex->edge = e;
        }
        backwardsRunner = a;
      }
      if (firstIntersection) {
        // Create new point at the point we're cracking.
        auto splitVertex = std::make_shared<DCEL::Vertex>(p, mPoints.size());
        mPoints.push_back(std::make_shared<Point>(p.x, p.y, splitVertex));
        mPoints.back()->index = mPoints.size() - 1;

        // Move backwards and fix
        while (backwardsRunner) {
          backwardsRunner->face->replacePoint(
              mPoints[backwardsRunner->vertex->pointIndex], mPoints.back());
          backwardsRunner->vertex = splitVertex;
          splitVertex->edge = backwardsRunner;
          backwardsRunner = backwardsRunner->vertexNextCW();
        }
      }
      firstIntersection = true;
    }

    runner = next;
  }

  for (auto &p : mPoints) p->clearForces();
  for (auto &t : mTriangles) t->computeForceOnNodes();
  auto residualVal = std::max(topEigenValue - 1.5f, 0.0f);
  for (auto &p : addResidualTo) {
    std::cout << "ADDING RESIDUAL " << residualVal << " TO " << p->index
              << std::endl;
    p->residualTensor = 1.0f * m(glm::normalize(topEigenVector)) * residualVal;
  }
}

void Surface::Relax() {
  for (auto &t : mTriangles) {
    t->relax();
  }
  for (auto &p : mPoints) {
    p->clearForces();
  }
  for (auto &t : mTriangles) {
    t->computeForceOnNodes();
  }
}

void Surface::computeSeperationTensors() {
  std::cout << "COMPUTING SEP TENSORS" << std::endl;
  for (auto &p : mPoints) p->computeSeperationTensor();
}

void Surface::draw(sf::RenderWindow &window) {
  for (auto &t : mTriangles) {
    //  sf::CircleShape middle(10);
    //  middle.setPosition(width / 2 - 10, height / 2 - 10);
    // window.draw(middle);
    t->draw(window);
    /*  sf::CircleShape shape(5);
      auto testW = t->getLargestEigenRep();
      auto centroid = t->getCentroid() - vec2(5, 5);
      shape.setPosition(centroid.x, centroid.y);
      centroid = t->getCentroid();
      shape.setFillColor(sf::Color(250, 0, 0));
      // window.draw(shape);
      auto angle = vangle(t->u);
      sf::RectangleShape line(sf::Vector2f(20, 3));
      line.setPosition(centroid.x, centroid.y);
      line.setFillColor(sf::Color(0, 0, 255));
      line.rotate(angle);
      // window.draw(line);
      line.rotate(-angle);
      angle = vangle(t->v);
      line.setFillColor(sf::Color(0, 255, 0));
      line.rotate(angle);
      // window.draw(line);
      auto stress = glm::length(testW);
      line.setSize(sf::Vector2f(10 * stress, 1));
      line.setFillColor(sf::Color(255, 255, 0));
      line.rotate(-angle);
      angle = vangle(testW);
      line.rotate(angle);
      // window.draw(line);
      line.rotate(180);
      // window.draw(line);*/
  }
  for (auto &p : mPoints) {
    auto v = p->getLargestEigenRep();
    vec2 v1, v2;
    float l1, l2;
    eigenDecomp(p->seperationTensor, v1, v2, l1, l2);
    if (l1 < 20 && l2 < 20) continue;
    auto l0 = std::max(l1, l2);
    vec2 v0 = l0 == l1 ? v1 : v2;
    if (v == vec2(0)) continue;
    sf::RectangleShape line(sf::Vector2f(l0, 4));
    line.setFillColor(sf::Color(199, 21, 133));
    line.setPosition(p->x, p->y);
    auto angle = vangle(v0);
    line.rotate(angle);
    window.draw(line);
    line.rotate(180);
    window.draw(line);
    // p->drawForces(window);
  }

  for (auto &pair : mCracks) {
    auto p0 = pair.first, p1 = pair.second;
    sf::VertexArray line(sf::LineStrip, 2);
    line[0].position = sf::Vector2f(p0->x, p0->y);
    line[0].color = sf::Color(255, 255, 0);
    line[1].position = sf::Vector2f(p1->x, p1->y);
    line[1].color = sf::Color(255, 255, 0);
    window.draw(line);
  }
  /* sf::Font font;
   if (!font.loadFromFile("arial.ttf")) {
     std::cout << " error... " << std::endl;
   }
   std::map<int64_t, int> existing;
   for (auto &p : mPoints) {
     sf::Text text;

     // select the font
     text.setFont(font);  // font is a sf::Font

     // set the string to display
     text.setString(std::to_string(p->index));
     int64_t key = (int64_t(p->x) << 32) + int64_t(p->y);
     auto y = p->y + existing[key] * 26;
     text.setPosition(sf::Vector2f(p->x, y));
     existing[key]++;

     // set the character size
     text.setCharacterSize(24);  // in pixels, not points!
     // window.draw(text);
   }*/
}

void Surface::PrintDataToLog() {
  std::ofstream file;
  file.open("Log.txt");
  file << "------------------ POINTS -----------------" << std::endl
       << std::endl;
  for (int i = 0; i < mPoints.size(); i++) {
    auto p = mPoints[i];
    file << i << "\tX: " << int(p->x) << "\tY: " << int(p->y) << std::endl;
  }

  file << std::endl
       << "------------------ CRACKS -----------------" << std::endl
       << std::endl;
  for (int i = 0; i < mCracks.size(); i++) {
    auto c = mCracks[i];
    auto p0 = c.first, p1 = c.second;
    file << p0->index << " " << p1->index << std::endl;
  }

  file << std::endl
       << "------------------ TRIANGLES -----------------" << std::endl
       << std::endl;
  int tr = 0;
  for (auto &t : mTriangles) {
    file << tr++ << " " << t->p0->index << " " << t->p1->index << " "
         << t->p2->index << std::endl;
  }
}
