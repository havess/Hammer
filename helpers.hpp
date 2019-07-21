#pragma once
#include "Eigen/Dense"
#include "Types.hpp"

float vangle(const vec2 &vec) { return glm::degrees(std::atan2(vec.y, vec.x)); }

mat2 rotate2(float degrees) {
  auto rad = glm::radians(degrees);
  auto s = std::sin(rad), c = std::cos(rad);
  return mat2(vec2(c, -s), vec2(s, c));
}

void eigenDecomp(const mat2 &A, vec2 &v1, vec2 &v2, float &l1, float &l2) {
  auto a = A[0][0], b = A[1][0], c = A[0][1], d = A[1][1];
  Eigen::Matrix2f m;
  m(0, 0) = a;
  m(0, 1) = b;
  m(1, 0) = c;
  m(1, 1) = d;
  Eigen::EigenSolver<Eigen::Matrix2f> s(m);
  auto ev = s.eigenvalues();
  auto evec = s.eigenvectors();
  l1 = ev(0, 0).real();
  l2 = ev(1, 0).real();
  if (std::abs(l1) < 0.00001) l1 = 0;
  if (std::abs(l2) < 0.00001) l2 = 0;
  v1 = vec2(evec(0, 0).real(), evec(1, 0).real());
  v2 = vec2(evec(0, 1).real(), evec(1, 1).real());
  // v1 = glm::normalize(v1);
  // v2 = glm::normalize(v2);
}

mat2 m(const vec2 &v) {
  if (glm::length(v) < 0.0001) return mat2(0);
  return glm::outerProduct(v, v) / glm::length(v);
}
