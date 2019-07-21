#pragma once
#include <algorithm>
#include <glm/glm.hpp>
#include <iostream>
#include <memory>
#include <queue>
#include <vector>
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>

enum class EventType {
  Site,
  Circle,
};

struct Event {
  Event(const vec2 &p, EventType type) : p(p), type(type) {}
  vec2 p;
  EventType type;
};

struct EventPtrComparator {
  bool operator()(const sptr<Event> &e1, const sptr<Event> &e2) {
    return e1->p.y > e2->p.y || (e1->p.y == e2->p.y && e1->p.x > e2->p.x);
  }
};

struct SiteEvent : public Event {
  SiteEvent(const vec2 &p) : Event(p, EventType::Site) {}
}

struct CircleEvent : public Event {
  CircleEvent(const vec2 &p) : Event(p, EventType::Circle) {}
  bool skip = false;
}

class Voronoi {
 public:
  Voronoi(const std::vector<vec2> &pts) : mPoints(pts) {}
  void Compute();

 private:
  void HandleCircleEvent();
  void HandleSiteEvent();

  std::vector<vec2> mPoints;
  std::priority_queue<sptr<Event>, std::vector<sptr<Event>>, EventPtrComparator>
      mQ;
};
