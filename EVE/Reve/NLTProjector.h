#ifndef REVE_NLTProjector
#define REVE_NLTProjector

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>

namespace Reve {

class Vector;

/**************************************************************************/
//  NLTProjections
/**************************************************************************/

class NLTProjection
{
public:
  enum PType_e { PT_Unknown, PT_CFishEye, PT_RhoZ }; // , PT_RhoPhi}; 

protected:
  PType_e             fType;
  const char*         fName;
  Float_t             fDistortion; // sensible values from 0 to 0.01
  Float_t             fFixedRadius;
  Float_t             fScale;
 
public:
  NLTProjection() : fType(PT_Unknown), fName(0), fDistortion(0), fFixedRadius(300), fScale(1.0f) {}
  virtual   ~NLTProjection() {}

  virtual   void      ProjectPoint(Float_t&, Float_t&, Float_t&){}
  virtual   void      ProjectPointFv(Float_t* v){ ProjectPoint(v[0], v[1], v[2]); }
  virtual   void      ProjectVector(Vector& v);
  virtual   Vector*   Project(Vector* pnts, Int_t npnts, Bool_t create_new = kTRUE);

  const     char*     GetName(){return fName;}
  void                SetName(const char* txt){ fName = txt; }

  void                SetType(PType_e t){fType = t;}
  PType_e             GetType(){return fType;}
  void                SetDistortion(Float_t d);
  Float_t             GetDistortion(){return fDistortion;}
  void                SetFixedRadius(Float_t x){fFixedRadius = x;}
  Float_t             GetFixedRadius(){return fFixedRadius;}

  virtual   Bool_t    AcceptSegment(Vector&, Vector&, Float_t /*tolerance*/) { return kTRUE; } 
  virtual   void      SetDirectionalVector(Int_t screenAxis, Vector& vec);

  // utils to draw axis
  virtual Float_t     GetValForScreenPos(Int_t ax, Float_t value);
  virtual Float_t     GetScreenVal(Int_t ax, Float_t value);

  static   Float_t    fgEps;

  ClassDef(NLTProjection, 0);
};

class RhoZ: public NLTProjection
{
public:
  RhoZ() : NLTProjection() { fType = PT_RhoZ; fName="RhoZ";}
  virtual ~RhoZ() {}

  virtual   Bool_t    AcceptSegment(Vector& v1, Vector& v2, Float_t tolerance); 
  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z);
  virtual   void      SetDirectionalVector(Int_t screenAxis, Vector& vec);

  ClassDef(RhoZ, 0);
};

class CircularFishEye : public NLTProjection
{
public:
  CircularFishEye():NLTProjection() { fType = PT_CFishEye; fName="CircularFishEye"; }
  virtual ~CircularFishEye() {}

  virtual   void      ProjectPoint(Float_t& x, Float_t& y, Float_t& z); 

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
