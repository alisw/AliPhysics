#ifndef REVE_NLTProjector
#define REVE_NLTProjector

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>
#include <Reve/PODs.h>

namespace Reve {


/**************************************************************************/
//  NLTProjections
/**************************************************************************/

class NLTProjection
{
public:
  enum PType_e   { PT_Unknown, PT_CFishEye, PT_RhoZ };  
  enum PProc_e  { PP_Plane, PP_Distort, PP_Full };
  enum GeoMode_e { GM_Unknown, GM_Polygons, GM_Segments };  

protected:
  PType_e             fType;
  GeoMode_e           fGeoMode;
  const char*         fName;

  Vector              fCenter;
  Vector              fProjectedZero;
  Vector              fZeroPosVal;

  Float_t             fDistortion; // sensible values from 0 to 0.01
  Float_t             fFixedRadius;
  Float_t             fScale;
  Vector              fUpLimit;  
  Vector              fLowLimit;

  
public:
  NLTProjection(Vector& center);
  virtual ~NLTProjection(){}

  // virtual   void      ProjectPoint(Float_t&, Float_t&, Float_t&, PProc_e p = PP_Full ){p =p;}
  virtual   void      ProjectPoint(Float_t&, Float_t&, Float_t&, PProc_e p = PP_Full ) = 0;
  virtual   void      ProjectPointFv(Float_t* v){ ProjectPoint(v[0], v[1], v[2]); }
  virtual   void      ProjectVector(Vector& v);
  virtual   Vector*   Project(Vector* pnts, Int_t npnts, Bool_t create_new = kTRUE);

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

  static   Float_t    fgEps;

  ClassDef(NLTProjection, 0);
};

class RhoZ: public NLTProjection
{
private:
  Float_t  fCenterR;
  Vector   fProjectedCenter;
public:
  RhoZ(Vector& center) : NLTProjection(center), fCenterR(0) { fType = PT_RhoZ; fName="RhoZ"; }
  virtual ~RhoZ() {}

  virtual   Bool_t    AcceptSegment(Vector& v1, Vector& v2, Float_t tolerance); 
  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc = PP_Full);
  virtual   void      SetDirectionalVector(Int_t screenAxis, Vector& vec);

  virtual   void      SetCenter(Vector& center); 
  virtual Float_t*    GetProjectedCenter() { return fProjectedCenter.c_vec(); }
  ClassDef(RhoZ, 0);
};

class CircularFishEye : public NLTProjection
{
public:
  CircularFishEye(Vector& center):NLTProjection(center) { fType = PT_CFishEye; fGeoMode = GM_Polygons; fName="CircularFishEye"; }
  virtual ~CircularFishEye() {}

  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z, PProc_e proc = PP_Full); 

  ClassDef(CircularFishEye, 0);
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
  
  NLTProjection*  fProjection;

  Bool_t          fDrawCenter;
  Bool_t          fDrawOrigin;
  Vector          fCenter;

  Int_t           fSplitInfoMode;
  Int_t           fSplitInfoLevel;
  Color_t         fAxisColor;  

  Float_t         fCurrentDepth;

public:
  NLTProjector();
  virtual ~NLTProjector();

  void            SetProjection(NLTProjection::PType_e type, Float_t distort=0);
  void            SetProjection(NLTProjection* p);
  NLTProjection*  GetProjection() { return fProjection; }

  virtual void    UpdateName();
  // scale info 
  void            SetSplitInfoMode(Int_t x)  { fSplitInfoMode = x;     }
  Int_t           GetSplitInfoMode()   const { return fSplitInfoMode;  }

  void            SetSplitInfoLevel(Int_t x) { fSplitInfoLevel = x;    }
  Int_t           GetSplitInfoLevel()  const { return fSplitInfoLevel; }

  void            SetAxisColor(Color_t col)  { fAxisColor = col;       }
  Color_t         GetAxisColor()       const { return fAxisColor;      }

  void            SetCurrentDepth(Float_t d) { fCurrentDepth = d;      }
  Float_t         GetCurrentDepth()    const { return fCurrentDepth;   }

  void            SetCenter(Float_t x, Float_t y, Float_t z);
  Vector&         GetCenter(){return fCenter;}

  void            SetDrawCenter(Bool_t x){ fDrawCenter = x; }
  Bool_t          GetDrawCenter(){ return fDrawCenter; }

  void            SetDrawOrigin(Bool_t x){ fDrawOrigin = x; }
  Bool_t          GetDrawOrigin(){ return fDrawOrigin; }

  virtual Bool_t  HandleElementPaste(RenderElement* el);

  virtual Bool_t  ShouldImport(RenderElement* rnr_el);
  virtual void    ImportElementsRecurse(RenderElement* rnr_el, RenderElement* parent);
  virtual void    ImportElements(RenderElement* rnr_el);

  virtual void    ProjectChildrenRecurse(RenderElement* rnr_el);
  virtual void    ProjectChildren();

  virtual void    ComputeBBox();
  virtual void    Paint(Option_t* option = "");

  ClassDef(NLTProjector, 0); //GUI for editing TGLViewer attributes
};

}
#endif
