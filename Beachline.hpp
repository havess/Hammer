//
//  Beachline.h
//  FortuneAlgo
//
//  Created by Dmytro Kotsur on 20/04/2018.
//  Copyright © 2018 Dmytro Kotsur. All rights reserved.
//
//  Modified by Me (Sam Haves), I did not want to write a custom balancing BST.
//  Changes are (mostly) just stylistic.

#pragma once

#include <time.h>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

#include "DCEL.hpp"
#include "Parabola.hpp"
#include "Types.hpp"

struct Event;

using namespace DCEL;

class Node;
using NodePtr = sptr<Node>;

class Node {
 public:
  // Height of the tree
  int height;

  // Pointer to a position of a sweepline
  double *sweepline;

  // Pointer to a vector of input points
  const std::vector<Point2D> *points;

  // Indices of the points
  std::pair<int, int> indices;

  // Pointers to left, right children and parent node
  NodePtr left, right, parent;

  // Pointer to a circle event for a leaf node or halfedge for an internal node
  sptr<Event> circle_event;
  sptr<HalfEdge> edge;

  // Constructor
  Node(const std::pair<int, int> &_indices, double *_sweepline = nullptr,
       const std::vector<Point2D> *_points = nullptr, NodePtr _left = nullptr,
       NodePtr _right = nullptr, NodePtr _parent = nullptr, int _height = 1);

  // Pointers to a next and previous arc-nodes
  NodePtr next, prev;

  // Leaf is defined as <p_i,p_i>
  inline bool is_leaf() { return indices.first == indices.second; }

  // Returns id of the node (only for leafs)
  inline int get_id() { return indices.first; }

  // Check if indices are equal
  inline bool has_indices(int a, int b) {
    return indices.first == a && indices.second == b;
  }

  // Check if indices are equal
  inline bool has_indices(const std::pair<int, int> &p) {
    return indices.first == p.first && indices.second == p.second;
  }

  // Return x-coordinate of:
  //  - in case of leaf node - corresponding focus of parabola;
  //  - in case of internal node - breakpoint;
  double value();
};

/**
 Connect as a list
 */
void connect(NodePtr prev, NodePtr next);

/**
 Check if the node is a root node
 */
bool is_root(NodePtr node);

/**
 Get height of the node
 */
int get_height(NodePtr node);

/**
 Update height of the node
 */
void update_height(NodePtr node);

/**
 Get balance of the node (difference between the height of left and right
 subtrees)
 */
int get_balance(NodePtr node);

/**
 Performs rotation of a tree around `node` such that it goes to the left subtree
 */
NodePtr rotate_left(NodePtr node);

/**
 Performs rotation of a tree around `node` such that it goes to the right
 subtree
 */
NodePtr rotate_right(NodePtr node);

/**
 Find a leaf in a tree such that x is under the parabolic arc,
 which corresponds to this leaf.
 */
NodePtr find(NodePtr root, double x);

/**
 Replace a leaf `node` with a new subtree, which has root `new_node`.
 The function rebalances the tree and returns the pointer to a new root node.
 */
NodePtr replace(NodePtr node, NodePtr new_node);

/**
 Remove a disappearing arc related to a circle event.
 The function rebalances the tree and returns the pointer to a new root node.
 */
NodePtr remove(NodePtr leaf);

/**
 Returns breakpoints for a given arc
 */
std::pair<NodePtr, NodePtr> breakpoints(NodePtr leaf);

NodePtr make_subtree(int index, int index_behind, double *sweepline,
                     const std::vector<Point2D> *points,
                     std::vector<HalfEdgePtr> &edges);

NodePtr make_simple_subtree(int index, int index_behind, double *sweepline,
                            const std::vector<Point2D> *points,
                            std::vector<HalfEdgePtr> &edges);

bool _validate(NodePtr node);

bool _check_balance(NodePtr node);

/**
 Print tree
 */
void print_tree(NodePtr root, int width = 7);

