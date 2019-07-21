#pragma once
//
//  VoronoiDiagram.hpp
//  FortuneAlgo
//
//  Created by Dmytro Kotsur on 06/05/2018.
//  Copyright © 2018 Dmytro Kotsur. All rights reserved.
//

#include "Beachline.hpp"
#include "Point2D.h"

void build_voronoi(const std::vector<Point2D> &points,
                   std::vector<HalfEdgePtr> &halfedges,
                   std::vector<VertexPtr> &vertices,
                   std::vector<HalfEdgePtr> &faces);

// std::vector<HalfEdgePtr> init
//
