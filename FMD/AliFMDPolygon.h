// -*- mode: C++ -*-
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * Latest changes by Christian Holm Christensen <cholm@nbi.dk>
 *
 * See cxx source for full Copyright notice                               
 */
#ifndef ALIFMDPOLYGON_H
#define ALIFMDPOLYGON_H
#ifndef ROOT_TVector2
# include <TVector2.h>
#endif
#ifndef ROOT_TObjArray
# include <TObjArray.h>
#endif

class AliFMDPolygon : public TObject 
{
private:
  enum {
    kUnknown, 
    kConvex, 
    kConcave
  };
  mutable Int_t fState;
  // List of coordinates 
  TObjArray fVerticies;
  // Force convexity check 
  bool ConvexCheck() const;
  // Check if a point is at the right-hand side of a segment 
  bool IsOnLeftHand(const TVector2* c, size_t i1, size_t i2) const;
public:
  // Construct a alipolygon with N sides
  AliFMDPolygon();
  virtual ~AliFMDPolygon();

  // Clear the polygon 
  void Clear(Option_t* option="");
  
  // Add a vertex 
  bool AddVertex(TVector2* c);
  bool AddVertex(double x, double y);

  // Get a vertex point 
  const TVector2& GetVertex(size_t i) const;
  
  // Check if a point is inside the polygon 
  bool Contains(const TVector2* c) const;
  bool Contains(double x, double y) const;
  
  // Get the number of verticies 
  size_t GetNVerticies() const { return fVerticies.GetEntries(); }  
  // Get the coordinates 
  const TObjArray& GetVerticies() const { return fVerticies; }
  
  void Draw(const char* option="PL", const char* name=0) const;

  ClassDef(AliFMDPolygon,1) // Polygon parameters
};

#endif
//
// EOF
//

  
      
