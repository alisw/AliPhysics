#ifndef REVE_NLTProjector
#define REVE_NLTProjector

#include <TAtt3D.h>
#include <TAttBBox.h>

#include <Reve/RenderElement.h>
#include <Reve/NLTProjections.h>
#include <Reve/PODs.h>

namespace Reve {

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

  ClassDef(NLTProjector, 0); // Manages and steers NLT projections.
};

}
#endif
