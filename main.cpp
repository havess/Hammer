#include <SFML/Audio.hpp>
#include <SFML/Graphics.hpp>
#include <algorithm>
#include <glm/glm.hpp>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <iostream>
#include <memory>
#include <random>
#include <vector>

#include "Beachline.hpp"
#include "Point2D.h"
#include "Types.hpp"
#include "VoronoiDiagram.hpp"

/**
 Converts vector of Point2D into vector of double.
 If coord_id == 0, function returns x-coordinate
 If coord_id == 1, function returns y-coordinate
 */
std::vector<double> get_coordinate(const std::vector<Point2D> &points,
                                   int coord_id) {
  std::vector<double> c(points.size(), 0.0);
  if (coord_id == 0) {
    for (size_t i = 0; i < points.size(); ++i) {
      c[i] = points[i].x;
    }
  } else {
    for (size_t i = 0; i < points.size(); ++i) {
      c[i] = points[i].y;
    }
  }
  return c;
}

/**
 Converts vector of Vertex pointers into vector of double.
 If coord_id == 0, function returns x-coordinate
 If coord_id == 1, function returns y-coordinate
 */
std::vector<double> get_coordinate(const std::vector<VertexPtr> &points,
                                   int coord_id) {
  std::vector<double> c(points.size(), 0.0);
  if (coord_id == 0) {
    for (size_t i = 0; i < points.size(); ++i) {
      c[i] = points[i]->point.x;
    }
  } else {
    for (size_t i = 0; i < points.size(); ++i) {
      c[i] = points[i]->point.y;
    }
  }
  return c;
}

double random_number(bool width) {
  return double(rand()) / double(RAND_MAX) *
             (width ? surfaceWidth : surfaceHeight) +
         delta;
  // return double(rand() % ((width ? surfaceWidth : surfaceHeight) + 1) +
  // delta);
}

std::vector<Point2D> randomPoint(int number) {
  srand(static_cast<unsigned int>(time(0)));
  std::vector<Point2D> points;
  for (int i = 0; i < number; ++i) {
    double x = random_number(true), y = random_number(false);
    points.push_back(Point2D(x, y));
  }
  std::sort(points.begin(), points.end(),
            [](const Point2D &p1, const Point2D &p2) {
              return (fabs(p1.y - p2.y) < POINT_EPSILON && p1.x < p2.x) ||
                     (fabs(p1.y - p2.y) >= POINT_EPSILON && p1.y < p1.y);
            });
  for (int i = 1; i < number; ++i) {
    if ((points[i - 1] - points[i]).norm() < POINT_EPSILON) {
      points[i - 1].x = random_number(true);
    }
  }
  for (int i = number - 1; i >= 0; --i) {
    int j = rand() % (i + 1);
    std::swap(points[i], points[j]);
  }
  return points;
}

void initEdgePointsVis(HalfEdgePtr h, std::vector<double> &x,
                       std::vector<double> &y,
                       const std::vector<Point2D> &points) {
  if (h->vertex != nullptr && h->twin->vertex != nullptr) {
    x[0] = h->vertex->point.x;
    x[1] = h->twin->vertex->point.x;

    y[0] = h->vertex->point.y;
    y[1] = h->twin->vertex->point.y;

  } else if (h->vertex != nullptr) {
    x[0] = h->vertex->point.x;
    y[0] = h->vertex->point.y;

    Point2D norm = (points[h->l_index] - points[h->r_index])
                       .normalized()
                       .getRotated90CCW();
    x[1] = x[0] + norm.x * 1000;
    y[1] = y[0] + norm.y * 1000;

  } else if (h->twin->vertex != nullptr) {
    x[0] = h->twin->vertex->point.x;
    y[0] = h->twin->vertex->point.y;

    Point2D norm = (points[h->twin->l_index] - points[h->twin->r_index])
                       .normalized()
                       .getRotated90CCW();
    x[1] = x[0] + norm.x * 1000;
    y[1] = y[0] + norm.y * 1000;

  } else {
    Point2D p1 = points[h->l_index], p2 = points[h->r_index];

    Point2D norm = (p1 - p2).normalized().getRotated90CCW();
    Point2D c = 0.5 * (p1 + p2);

    x[0] = c.x + norm.x * 1000;
    x[1] = c.x - norm.x * 1000;

    y[0] = c.y + norm.y * 1000;
    y[1] = c.y - norm.y * 1000;
  }
}

Surface GenerateRegularRectangle() {
  /****** INITIAL GRID ******/
  std::set<sptr<Triangle>> triangles;
  std::vector<sptr<Point>> vpoints;
  sptr<Point> **points = new sptr<Point> *[(surfaceWidth / cellSize) + 1];
  for (int i = 0; i <= surfaceWidth / cellSize; i++)
    points[i] = new sptr<Point>[(surfaceHeight / cellSize) + 1];

  std::cout << "SIZE: " << surfaceWidth / cellSize + 1 << " "
            << surfaceHeight / cellSize + 1 << std::endl;

  /********** CREATE GRID POINTS **********/
  for (int x = delta; x <= surfaceWidth + delta; x += cellSize) {
    for (int y = delta; y <= surfaceHeight + delta; y += cellSize) {
      int indX = (x - delta) / cellSize;
      int indY = (y - delta) / cellSize;
      points[indX][indY] = std::make_shared<Point>(x, y);
      vpoints.push_back(points[indX][indY]);
    }
  }

  /************ CREATE TRIANGLES ***********/
  for (int i = 0; i <= surfaceWidth / cellSize - 1; i++) {
    for (int j = 0; j <= surfaceHeight / cellSize - 1; j++) {
      auto t1 = std::make_shared<Triangle>(points[i][j], points[i + 1][j],
                                           points[i + 1][j + 1]);
      auto t2 = std::make_shared<Triangle>(points[i][j], points[i + 1][j + 1],
                                           points[i][j + 1]);
      triangles.insert(t1);
      triangles.insert(t2);
    }
  }
  return Surface(triangles, vpoints);
}

Surface GenerateDiscreteRectangle(sf::RenderWindow &window, int nPts) {
  /****** INITIAL GRID ******/
  std::set<sptr<Triangle>> triangles;
  std::vector<sptr<Point>> points;
  std::vector<Point2D> pts = randomPoint(nPts);
  std::vector<HalfEdgePtr> halfedges, faces;
  std::vector<VertexPtr> vertices;

  // Add edge points so that the edges arent poorly formed triangles.
  for (int x = delta - 20; x <= surfaceWidth + delta + 20; x += 80) {
    pts.push_back(Point2D(x, delta - 20));
    pts.push_back(Point2D(x, surfaceHeight + delta + 20));
  }
  for (int y = delta; y <= surfaceHeight + delta; y += 80) {
    pts.push_back(Point2D(delta - 20, y));
    pts.push_back(Point2D(surfaceWidth + delta + 20, y));
  }

  // Construct initial Voronoi diagram
  build_voronoi(pts, halfedges, vertices, faces);
  /********** DRAW VORONOI DIAGRAM *********/
  /*for (size_t i = 0; i < halfedges.size(); ++i) {
    HalfEdgePtr h = halfedges[i];

    std::vector<double> x(2, 0.0), y(2, 0.0);
    initEdgePointsVis(h, x, y, pts);
    sf::VertexArray line(sf::LineStrip, 2);
    line[0].position = sf::Vector2f(x[0], y[0]);
    line[1].position = sf::Vector2f(x[1], y[1]);
    window.draw(line);
  }
  for (auto &p : pts) {
    sf::CircleShape pt(3);
    pt.setFillColor(sf::Color::Red);
    pt.setPosition(sf::Vector2f(p.x - 3, p.y - 3));
    window.draw(pt);
  }*/
  /******* RELAX MESH *********/
  for (int i = 0; i < 3; i++) {
    for (auto &e : faces) {
      if (!e->is_finite()) continue;
      auto runner = e->next;
      int numSides = 0;
      Point2D centroid;
      auto a = e->twin->vertex->point;
      bool open = false;
      while (runner && runner != e) {
        if (!runner->is_finite()) {
          open = true;
          break;
        }
        auto b = runner->vertex->point;
        auto c = runner->twin->vertex->point;
        centroid += (a + b + c) / 3.0f;
        numSides++;
        runner = runner->next;
      }
      if (open) continue;
      centroid = centroid / numSides;
      if (centroid.x > delta + surfaceWidth || centroid.x < delta) continue;
      if (centroid.y > delta + surfaceHeight || centroid.y < delta) continue;
      pts[e->l_index].setX(centroid.x);
      pts[e->l_index].setY(centroid.y);
    }

    halfedges.clear();
    vertices.clear();
    faces.clear();
    build_voronoi(pts, halfedges, vertices, faces);
  }

  std::map<int, int> added;
  std::map<int64_t, HalfEdgePtr> twinLookupTable;
  auto createKey = [](int a, int b) {
    int64_t key = a < b ? a : b;
    key = key << 32;
    key += a < b ? b : a;
    return key;
  };

  /*********** CREATE DELAUNAY TRIANGULATION *********/
  for (auto &vert : vertices) {
    auto edge = vert->edge;
    auto a = pts[edge->l_index];
    if (!added.count(edge->l_index)) {
      points.push_back(std::make_shared<Point>(
          a.x, a.y,
          std::make_shared<Vertex>(Point2D(a.x, a.y), points.size())));
      points.back()->index = points.size() - 1;
      added[edge->l_index] = points.size() - 1;
    }
    auto b = pts[edge->r_index];
    if (!added.count(edge->r_index)) {
      points.push_back(std::make_shared<Point>(
          b.x, b.y,
          std::make_shared<Vertex>(Point2D(b.x, b.y), points.size())));
      points.back()->index = points.size() - 1;
      added[edge->r_index] = points.size() - 1;
    }
    auto ind = edge->next ? edge->next->r_index : edge->prev->r_index;
    auto c = pts[ind];
    if (!added.count(ind)) {
      points.push_back(std::make_shared<Point>(
          c.x, c.y,
          std::make_shared<Vertex>(Point2D(c.x, c.y), points.size())));
      points.back()->index = points.size() - 1;
      added[ind] = points.size() - 1;
    }

    auto tri = std::make_shared<Triangle>(points[added[edge->l_index]],
                                          points[added[edge->r_index]],
                                          points[added[ind]]);
    triangles.insert(tri);
    // Create half edges for this triangle
    auto halfEdgeA =
        std::make_shared<DCEL::HalfEdge>(edge->l_index, edge->r_index, tri);
    auto halfEdgeB = std::make_shared<DCEL::HalfEdge>(edge->r_index, ind, tri);
    auto halfEdgeC = std::make_shared<DCEL::HalfEdge>(edge->l_index, ind, tri);

    // Get the key for each new edge.
    auto keyA = createKey(edge->l_index, edge->r_index);
    auto keyB = createKey(edge->r_index, ind);
    auto keyC = createKey(edge->l_index, ind);

    halfEdgeA->twin = twinLookupTable[keyA];
    if (!halfEdgeA->twin) {
      twinLookupTable[keyA] = halfEdgeA;
    } else {
      twinLookupTable[keyA]->twin = halfEdgeA;
    }

    halfEdgeB->twin = twinLookupTable[keyB];
    if (!halfEdgeB->twin) {
      twinLookupTable[keyB] = halfEdgeB;
    } else {
      twinLookupTable[keyB]->twin = halfEdgeB;
    }

    halfEdgeC->twin = twinLookupTable[keyC];
    if (!halfEdgeC->twin) {
      twinLookupTable[keyC] = halfEdgeC;
    } else {
      twinLookupTable[keyC]->twin = halfEdgeC;
    }

    halfEdgeA->next = halfEdgeB;
    halfEdgeA->prev = halfEdgeC;

    halfEdgeB->next = halfEdgeC;
    halfEdgeB->prev = halfEdgeA;

    halfEdgeC->next = halfEdgeA;
    halfEdgeC->prev = halfEdgeB;

    halfEdgeA->vertex = points[added[edge->r_index]]->vertex;
    points[added[edge->r_index]]->vertex->edge = halfEdgeA;

    halfEdgeB->vertex = points[added[ind]]->vertex;
    points[added[ind]]->vertex->edge = halfEdgeB;

    halfEdgeC->vertex = points[added[edge->l_index]]->vertex;
    points[added[edge->l_index]]->vertex->edge = halfEdgeC;
  }

  /*********** CREATE SURFACE ************/
  return Surface(triangles, points);
}

int main() {
  /************ SFML SETUP **************/
  sf::ContextSettings settings;
  settings.antialiasingLevel = 8;
  // Create the main window
  sf::RenderWindow window(sf::VideoMode(width, height), "SFML window",
                          sf::Style::Default, settings);
  sf::Font font;
  if (!font.loadFromFile("/Library/Fonts/Arial.ttf")) return EXIT_FAILURE;
  window.clear(sf::Color(40, 40, 40));

  auto s = GenerateDiscreteRectangle(window, 3000);
  auto sR = GenerateRegularRectangle();
  /*********** COMPUTE POINT SEPERATION TENSORS ***********/
  std::cout << "COMPUTING SEPERATION TENSORS" << std::endl;
  s.computeSeperationTensors();
  s.initFractureState();

  bool keepGoing = false;
  while (window.isOpen()) {
    // Process events
    sf::Event event;
    while (window.pollEvent(event)) {
      // Close window: exit
      if (event.type == sf::Event::Closed) window.close();
      if (sf::Keyboard::isKeyPressed(sf::Keyboard::Left)) {
        keepGoing = true;
      }
      if (sf::Keyboard::isKeyPressed(sf::Keyboard::Right)) {
        s.PrintDataToLog();
      }

      //}
    }
    if (keepGoing) {
      s.initFractureState();
      s.FractureStep(window);
      for (int i = 0; i < 3; i++) s.Relax();
      s.computeSeperationTensors();
      // keepGoing = false;
    }
    // Clear screen
    window.clear(sf::Color(40, 40, 40));
    s.draw(window);
    // Draw the string
    window.display();
  }
  return EXIT_SUCCESS;
}
