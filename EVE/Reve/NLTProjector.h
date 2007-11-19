#ifndef REVE_NLTProjector
#define REVE_NLTProjector

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>
#include <Reve/PODs.h>

namespace Reve {

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
};

class NLTCircularFishEye : public NLTProjection
{
public:
  NLTCircularFishEye(Vector& center):NLTProjection(center) { fType = PT_CFishEye; fGeoMode = GM_Polygons; fName="CircularFishEye"; }
  virtual ~NLTCircularFishEye() {}

  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc = PP_Full);

  ClassDef(NLTCircularFishEye, 0); // XY non-linear projection.
};

/**************************************************************************/
//  NLTProjector
/**************************************************************************/
class NLTProjector : public RenderElementList,
		     public TAttBBox,
                     public TAtt3D
{
private:
  NLTProjector(const NLTProjector&);            // Not implemented
  NLTProjector& operator=(const NLTProjector&); // Not implemented

  NLTProjection*  fProjection;  // projection

  Bool_t          fDrawCenter;  // draw center of distortion
  Bool_t          fDrawOrigin;  // draw origin
  Vector          fCenter;      // center of distortion

  Int_t           fSplitInfoMode;  // tick-mark position
  Int_t           fSplitInfoLevel; // tick-mark density
  Color_t         fAxisColor;      // color of axis

  Float_t         fCurrentDepth;   // z depth of object being projected

  virtual Bool_t  ShouldImport(RenderElement* rnr_el);

public:
  NLTProjector();
  virtual ~NLTProjector();

  void            SetProjection(NLTProjection::PType_e type, Float_t distort=0);
  NLTProjection*  GetProjection() { return fProjection; }

  virtual void    UpdateName();

  void            SetAxisColor(Color_t col)  { fAxisColor = col;       }
  Color_t         GetAxisColor()       const { return fAxisColor;      }
  void            SetSplitInfoMode(Int_t x)  { fSplitInfoMode = x;     }
  Int_t           GetSplitInfoMode()   const { return fSplitInfoMode;  }
  void            SetSplitInfoLevel(Int_t x) { fSplitInfoLevel = x;    }
  Int_t           GetSplitInfoLevel()  const { return fSplitInfoLevel; }

  void            SetDrawCenter(Bool_t x){ fDrawCenter = x; }
  Bool_t          GetDrawCenter(){ return fDrawCenter; }
  void            SetDrawOrigin(Bool_t x){ fDrawOrigin = x; }
  Bool_t          GetDrawOrigin(){ return fDrawOrigin; }

  void            SetCenter(Float_t x, Float_t y, Float_t z);
  Vector&         GetCenter(){return fCenter;}

  void            SetCurrentDepth(Float_t d) { fCurrentDepth = d;      }
  Float_t         GetCurrentDepth()    const { return fCurrentDepth;   }

  virtual Bool_t  HandleElementPaste(RenderElement* el);
  virtual void    ImportElementsRecurse(RenderElement* rnr_el, RenderElement* parent);
  virtual void    ImportElements(RenderElement* rnr_el);
  virtual void    ProjectChildren();
  virtual void    ProjectChildrenRecurse(RenderElement* rnr_el);

  virtual void    ComputeBBox();
  virtual void    Paint(Option_t* option = "");

  ClassDef(NLTProjector, 0); // Project NLTProjectable object.
};

}
#endif
