#ifndef REVE_NLTProjector
#define REVE_NLTProjector

#include "TNamed.h"
#include "TString.h"

namespace std {
template<typename _Tp> class allocator;
template<typename _Tp, typename _Alloc > class list;
}

class TGeoManager;
class TGeoMatrix;
class TBuffer3D;

namespace Reve {
class Vector;
class PointSet;
class NLTPolygonSet;
class NLTPolygon;
class TrackList;
/**************************************************************************/
//  NLTProjections
/**************************************************************************/
class NLTProjection {
public: 
  enum Type_e       { CFishEye, RhoZ, RhoPhi}; 
 
  Float_t             fDistortion; // sensible values from 0 to 0.01
 
  virtual   Bool_t    AcceptSegment(Vector& /*v1*/, Vector& /*v2*/) {return kTRUE;}
  virtual   Vector*   Project(Vector* pnts, Int_t npnts, Bool_t create_new = kTRUE);
  NLTProjection():fDistortion(0.) {}
  ClassDef(NLTProjection, 0)
};

class PhiZ: public   NLTProjection {
};

class RhoZ: public NLTProjection {
public:
  virtual   Bool_t    AcceptSegment(Vector& v1, Vector& v2); 
  virtual   Vector*   Project(Vector* pnts, Int_t npnts, Bool_t copy = kTRUE);  
  ClassDef(RhoZ, 0)
};

class CircularFishEye: public NLTProjection {
public:
  virtual   Vector*   Project(Vector* pnts, Int_t npnts, Bool_t copy = kTRUE); 
  ClassDef(CircularFishEye, 0)
};

/**************************************************************************/
//  NLTProjector
/**************************************************************************/
class NLTProjector : public TNamed
{ 
protected:
  void   CheckPoint(Int_t idx, Float_t* bbox);
  Bool_t IsFirstIdxHead(Int_t s0, Int_t s1, TBuffer3D* buff);

private:
  NLTProjection*  fProjection;

  Float_t         fEps;    // distance accounted in reducing the ponts

  // temporary variables cashed 
  Int_t*          fIdxMap; // map from original to projected and reduced point needed oly for geometry
  Int_t           fNRPnts; // number of reduced and projected points
  Vector*         fRPnts;  // reduced and projected points

  void            ReducePoints(Vector* p, Int_t N);
  void            MakePolygonsFromBP(TBuffer3D* buff, std::list<NLTPolygon, std::allocator<NLTPolygon> >& pols);
  void            MakePolygonsFromBS(TBuffer3D* buff, std::list<NLTPolygon, std::allocator<NLTPolygon> >& pols);
  void            CleanUp();

public:
  NLTProjector();
  virtual ~NLTProjector();

  NLTPolygonSet*  ProjectGeoShape(TBuffer3D* buff, Int_t useBuffPols=-1);
  void            ProjectPointSet(PointSet* ps);
  void            ProjectTrackList(TrackList* tl);
  
  void            SetProjection(NLTProjection* p){ fProjection = p;}
  void            SetProjection(NLTProjection::Type_e type, Float_t distort = 0.);
  NLTProjection*  GetProjection(){ return fProjection; }

  void            DumpBuffer(TBuffer3D* b);
  ClassDef(NLTProjector, 0) //GUI for editing TGLViewer attributes
  };
}
#endif
