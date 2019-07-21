#include "Voronoi.hpp"

/****** This file contains the code that generates a voronoi diagram ******/
/// I use Fortunes Algorith, also known as the beach line algorithm to compute
/// the voronoi diagram of a set of 2D points. This is going to be hard,
/// lord help me.

/***** Pseudocode *****/

/// INPUT: A set P of points on the 2D plane
/// OUTPUT: The Voronoi Diagram inside a bounding box given using a double
/// connected edge list.

/// 1. Initialize an empty event queue Q, with all site events. Initialize a
/// empty status structure T, and an empty DC Edge List D.

/// 2. while Q is not empty:
/// 3.  do Remove the event with the largest y coordinate from Q.
/// 4.    if the event is a site event, occuring at site pi.
/// 5.      HANDLE_SITE_EVENT(pi);
/// 6.    else
/// 7.      HANDLE_CIRCLE_EVENT(L) where L is the leaf of T representing the
///         arc that will disappear.
/// 8. The internal nodes still present in T correspond to the half-infinite
/// edges of the Voronoi diagram. Computing the bounding box that contains all
/// vertices fo the Voronoi Diagram and attach the half-infinite edges to the
/// bounding box by updating D appropriately.
/// 9. Traverse the half edges of D to form cells.

void Voronoi::Compute() {
  // Initialize Q with all site events.
  for (const auto &v : mPoints) {
    mQ.push(std::make_shared<SiteEvent>(v));
  }
  double sweepPos = 0;

  while (!mQ.empty()) {
    auto event = mQ.top();
    mQ.pop();

    sweepPos = event->p.y;

    if (event->type == EventType::Site) {
      HandleSiteEvent(event);
    } else if (event->type == EventType::Circle) {
      HandleCircleEvent(event);
    }
  }

  /// 8. The internal nodes still present in T correspond to the half-infinite
  /// edges of the Voronoi diagram. Computing the bounding box that contains all
  /// vertices fo the Voronoi Diagram and attach the half-infinite edges to the
  /// bounding box by updating D appropriately.
  /// 9. Traverse the half edges of D to form cells.
}

/// 1. If T is empty, insert pi into it and return.
/// 2. Search in T for the arc alpha above pi. If the leaf representing
/// alpha has a pointer to a circle event in Q, then that event is a false
/// alarm and must be deleted from Q.
/// 3. Replace the leaf in T that represents alpha with a subtree having three
/// leaves. The middle leaf stores the new site pi, and the other two store the
/// site pj that was previously stored in alpha. Store the tuples <pi, pj> and
/// <pj, pi> representing the two new breakpoints at the two new internal nodes.
/// Rebalance T if necessary.
/// 4. Create new half-edge records in the Voronoi Diagram structure for
/// the new edge seperating V(pi) and V(pj) which will be traced out by the two
/// new breakpoints.
/// 5. Check the triple of consecutive arcs where the new arc for pi is the left
/// arc to see if the breakpoints converge. If so, insert the circle event into
/// Q and add pointers between node in T and node in Q. Do the same for the
/// triple where the new arc is the right arc.
void Voronoi::HandleSiteEvent(sptr<Event> siteEvent) {
  auto event = std::static_cast<sptr<SiteEvent>>(siteEvent);
  if (mT == nullptr) {
    mT = std::make_shared<Node>();
  } else {
    auto arc = BeachLine::Find(mT, event->p.x);
    sptr<Node> subtree, leftLeaf, rightLeaf;
    if (arc->circleEvent) {
      arc->circleEvent->skip = true;
    }

    // TODO Implement ComputeIntersections
    auto numIntersections = ComputeIntersections(event);
  }
}

/// 1. Delete leaf that represent disappearing arc alpha from T. Update the
/// tuples representing the breakpoints at the internal nodes. Rebalance T if
/// required. Delete all circle events involving alpha from Q. These can be
/// found by using the successor and predecessor pointers of leaf in T.
/// 2. Add the center of the circle causing the event as a vertex record to D.
/// Create 2 HEdge records corresponding to the new breakpoint of the beachline.
/// Set ptrs between them appropriately. Attach the three new records to the
/// HEdge records that end at the vertex.
/// 3. Check the new triple of consecutive arcs that has the former left
/// neighbor of alpha as its middle arc to see if the two breakpoints of the
/// triple converge. If so, insert the corresponding circle event into Q, and
/// set pointers between the new circle event in Q and the corresponding leaf of
/// T. Do the same for the triple where the former right neighbor is the middle
/// arc.
///
void Voronoi::HandleCircleEvent(const Node *leaf) {}
