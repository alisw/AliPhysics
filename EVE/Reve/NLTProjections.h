#ifndef REVE_NLTProjections
#define REVE_NLTProjections

#include <Reve/PODs.h>

namespace Reve {

////////////////////////////////////////////////////////////////
//                                                            //
// NLTProjection                                              //
//                                                            //
////////////////////////////////////////////////////////////////

class NLTProjection
{
public:
  enum PType_e   { PT_Unknown, PT_CFishEye, PT_RhoZ };     // type
  enum PProc_e   { PP_Plane, PP_Distort, PP_Full };        // procedure
  enum GeoMode_e { GM_Unknown, GM_Polygons, GM_Segments }; // reconstruction of geometry

protected:
  PType_e             fType;          // type
  GeoMode_e           fGeoMode;       // way of polygon reconstruction
  const char*         fName;          // name

  Vector              fCenter;        // center of distortion
  Vector              fZeroPosVal;    // projected origin (0, 0, 0)

  Float_t             fDistortion;    // distortion
  Float_t             fFixedRadius;   // projected radius independent of distortion
  Float_t             fScale;         // scale factor to keep projected radius fixed
  Vector              fUpLimit;       // convergence of point +infinity
  Vector              fLowLimit;      // convergence of point -infinity

public:
  NLTProjection(Vector& center);
  virtual ~NLTProjection(){}

  virtual   void      ProjectPoint(Float_t&, Float_t&, Float_t&, PProc_e p = PP_Full ) = 0;
  virtual   void      ProjectPointFv(Float_t* v){ ProjectPoint(v[0], v[1], v[2]); }
  virtual   void      ProjectVector(Vector& v);

  const     char*     GetName(){return fName;}
  void                SetName(const char* txt){ fName = txt; }

  virtual void        SetCenter(Vector& v){ fCenter = v; UpdateLimit();}
  virtual Float_t*    GetProjectedCenter() { return fCenter.c_vec(); }

  void                SetType(PType_e t){fType = t;}
  PType_e             GetType(){return fType;}

  void                SetGeoMode(GeoMode_e m){fGeoMode = m;}
  GeoMode_e           GetGeoMode(){return fGeoMode;}

  void                UpdateLimit();
  void                SetDistortion(Float_t d);
  Float_t             GetDistortion(){return fDistortion;}
  void                SetFixedRadius(Float_t x);
  Float_t             GetFixedRadius(){return fFixedRadius;}

  virtual   Bool_t    AcceptSegment(Vector&, Vector&, Float_t /*tolerance*/) { return kTRUE; }
  virtual   void      SetDirectionalVector(Int_t screenAxis, Vector& vec);

  // utils to draw axis
  virtual Float_t     GetValForScreenPos(Int_t ax, Float_t value);
  virtual Float_t     GetScreenVal(Int_t ax, Float_t value);
  Float_t             GetLimit(Int_t i, Bool_t pos) { return pos ? fUpLimit[i] : fLowLimit[i]; }

  static   Float_t    fgEps;  // resolution of projected points

  ClassDef(NLTProjection, 0); // Base-class for non-linear projection.
}; // endclass NLTProjection

////////////////////////////////////////////////////////////////
//                                                            //
// NLTRhoZ                                                    //
//                                                            //
////////////////////////////////////////////////////////////////

class NLTRhoZ: public NLTProjection
{
private:
  Vector   fProjectedCenter; // projected center of distortion.
public:
  NLTRhoZ(Vector& center) : NLTProjection(center) { fType = PT_RhoZ; fName="RhoZ"; }
  virtual ~NLTRhoZ() {}

  virtual   Bool_t    AcceptSegment(Vector& v1, Vector& v2, Float_t tolerance);
  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc = PP_Full);
  virtual   void      SetDirectionalVector(Int_t screenAxis, Vector& vec);

  virtual   void      SetCenter(Vector& center);
  virtual Float_t*    GetProjectedCenter() { return fProjectedCenter.c_vec(); }
  ClassDef(NLTRhoZ, 0);  // Rho/Z non-linear projection.
}; // endclass NLTRhoZ

////////////////////////////////////////////////////////////////
//                                                            //
// NLTCircularFishEye                                         //
//                                                            //
////////////////////////////////////////////////////////////////

class NLTCircularFishEye : public NLTProjection
{
public:
  NLTCircularFishEye(Vector& center):NLTProjection(center) { fType = PT_CFishEye; fGeoMode = GM_Polygons; fName="CircularFishEye"; }
  virtual ~NLTCircularFishEye() {}

  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc = PP_Full);

  ClassDef(NLTCircularFishEye, 0); // XY non-linear projection.
}; // endclass NLTCircularFishEye

}
#endif
