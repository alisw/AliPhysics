#ifndef __FASTJET__VORONOI_H__
#define __FASTJET__VORONOI_H__

//STARTHEADER
// $Id: Voronoi.hh 1161 2008-04-01 20:29:38Z soyez $
//
// Copyright (c) 1994 by AT&T Bell Laboratories (see below)
//
//
//----------------------------------------------------------------------
// This file is included as part of FastJet but was mostly written by
// S. Fortune in C, put into C++ with memory management by S
// O'Sullivan, and with further interface and memeory management
// modifications by Gregory Soyez.
//
// Permission to use, copy, modify, and distribute this software for
// any purpose without fee is hereby granted, provided that this
// entire notice is included in all copies of any software which is or
// includes a copy or modification of this software and in all copies
// of the supporting documentation for such software. THIS SOFTWARE IS
// BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED WARRANTY.
// IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY REPRESENTATION
// OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY OF THIS
// SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
//
//----------------------------------------------------------------------
//ENDHEADER


/*
* The author of this software is Steven Fortune.  
* Copyright (c) 1994 by AT&T Bell Laboratories.
* Permission to use, copy, modify, and distribute this software for any
* purpose without fee is hereby granted, provided that this entire notice
* is included in all copies of any software which is or includes a copy
* or modification of this software and in all copies of the supporting
* documentation for such software.
* THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
* OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
*/

/* 
* This code was originally written by Stephan Fortune in C code.  I,
* Shane O'Sullivan, have since modified it, encapsulating it in a C++
* class and, fixing memory leaks and adding accessors to the Voronoi
* Edges.  Permission to use, copy, modify, and distribute this
* software for any purpose without fee is hereby granted, provided
* that this entire notice is included in all copies of any software
* which is or includes a copy or modification of this software and in
* all copies of the supporting documentation for such software.  THIS
* SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
* WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
* REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE
* MERCHANTABILITY OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR
* PURPOSE.
*/

#include "fastjet/ClusterSequenceWithArea.hh"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DELETED -2
#define le 0
#define re 1

//using namespace std;

namespace fastjet {      // defined in fastjet/internal/base.hh

/**
 * \class Point
 * class to handle a 2d point
 */
class Point{
public:
  /// defailt ctor
  Point() : x(0.0), y(0.0) {}

  /// ctor with initialisation
  Point(double _x, double _y) : x(_x), y(_y) {}

  /// addition
  inline Point operator + (const Point &p) const{
    return Point(x+p.x, y+p.y);
  }

  /// subtraction
  inline Point operator - (const Point &p) const{
    return Point(x-p.x, y-p.y);
  }

  /// scalar multiplication
  inline Point operator * (const double t) const{
    return Point(x*t, y*t);
  }

  /// vector coordinates
  double x,y;
};


/// norm of a vector
inline double norm(const Point p){
  return p.x*p.x+p.y*p.y;
}


/// 2D vector product
inline double vector_product(const Point &p1, const Point &p2){
  return p1.x*p2.y-p1.y*p2.x;
}


/// scalar product
inline double scalar_product(const Point &p1, const Point &p2){
  return p1.x*p2.x+p1.y*p2.y;
}


/**
 * \class GraphEdge
 * handle an edge of the Voronoi Diagram.
 */
class GraphEdge{
public:
  /// coordinates of the extreme points
  double x1,y1,x2,y2;

  /// indices of the parent sites that define the edge
  int point1, point2;

  /// pointer to the next edge
  GraphEdge* next;
};


/**
 * \class Site
 * structure used both for particle sites and for vertices.
 */
class Site{
 public:
  Point	coord;
  int sitenbr;
  int refcnt;
};



class Freenode{
public:
  Freenode *nextfree;
};


class FreeNodeArrayList{
public:
  Freenode* memory;
  FreeNodeArrayList* next;
};


class Freelist{
public:
  Freenode *head;
  int nodesize;
};

class Edge{
public:
  double a,b,c;
  Site *ep[2];
  Site *reg[2];
  int edgenbr;
};


class Halfedge{
public:
  Halfedge *ELleft, *ELright;
  Edge *ELedge;
  int ELrefcnt;
  char ELpm;
  Site *vertex;
  volatile double ystar;
  Halfedge *PQnext;
};


class VoronoiDiagramGenerator{
public:
  VoronoiDiagramGenerator();
  ~VoronoiDiagramGenerator();

  bool generateVoronoi(std::vector<Point> *_parent_sites,
		       double minX, double maxX, double minY, double maxY, 
		       double minDist=0);

  inline void resetIterator(){
    iteratorEdges = allEdges;
  }

  bool getNext(GraphEdge **e){
    if(iteratorEdges == 0)
      return false;
    
    *e = iteratorEdges;
    iteratorEdges = iteratorEdges->next;
    return true;
  }
  
  std::vector<Point> *parent_sites;
  int n_parent_sites;

private:
  void cleanup();
  void cleanupEdges();
  char *getfree(Freelist *fl);	
  Halfedge *PQfind();
  int PQempty();
	
  Halfedge **ELhash;
  Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
  Halfedge *HEcreate(Edge *e,int pm);
  
  Point PQ_min();
  Halfedge *PQextractmin();	
  void freeinit(Freelist *fl,int size);
  void makefree(Freenode *curr,Freelist *fl);
  void geominit();
  void plotinit();
  bool voronoi(int triangulate);
  void ref(Site *v);
  void deref(Site *v);
  void endpoint(Edge *e,int lr,Site * s);

  void ELdelete(Halfedge *he);
  Halfedge *ELleftbnd(Point *p);
  Halfedge *ELright(Halfedge *he);
  void makevertex(Site *v);
  void out_triple(Site *s1, Site *s2,Site * s3);
  
  void PQinsert(Halfedge *he,Site * v, double offset);
  void PQdelete(Halfedge *he);
  bool ELinitialize();
  void ELinsert(Halfedge *lb, Halfedge *newHe);
  Halfedge * ELgethash(int b);
  Halfedge *ELleft(Halfedge *he);
  Site *leftreg(Halfedge *he);
  void out_site(Site *s);
  bool PQinitialize();
  int PQbucket(Halfedge *he);
  void clip_line(Edge *e);
  char *myalloc(unsigned n);
  int right_of(Halfedge *el,Point *p);

  Site *rightreg(Halfedge *he);
  Edge *bisect(Site *s1, Site *s2);
  double dist(Site *s,Site *t);
  Site *intersect(Halfedge *el1, Halfedge *el2, Point *p=0);

  void out_bisector(Edge *e);
  void out_ep(Edge *e);
  void out_vertex(Site *v);
  Site *nextone();

  void pushGraphEdge(double x1, double y1, double x2, double y2, 
		     Site *s1, Site *s2);

  void openpl();
  void circle(double x, double y, double radius);
  void range(double minX, double minY, double maxX, double maxY);

  Freelist hfl;
  Halfedge *ELleftend, *ELrightend;
  int ELhashsize;
  
  int triangulate, sorted, plot, debug;
  double xmin, xmax, ymin, ymax, deltax, deltay;
  
  Site *sites;
  int nsites;
  int siteidx;
  int sqrt_nsites;
  int nvertices;
  Freelist sfl;
  Site *bottomsite;
  
  int nedges;
  Freelist efl;
  int PQhashsize;
  Halfedge *PQhash;
  int PQcount;
  int PQmin;
  
  int ntry, totalsearch;
  double pxmin, pxmax, pymin, pymax, cradius;
  int total_alloc;
  
  double borderMinX, borderMaxX, borderMinY, borderMaxY;
  
  FreeNodeArrayList* allMemoryList;
  FreeNodeArrayList* currentMemoryBlock;
  
  GraphEdge* allEdges;
  GraphEdge* iteratorEdges;
  
  double minDistanceBetweenSites;
};

int scomp(const void *p1,const void *p2);


} // fastjet namespace 

#endif
