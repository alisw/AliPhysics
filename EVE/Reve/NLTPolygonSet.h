#ifndef REVE_NLTPolygonSet_H
#define REVE_NLTPolygonSet_H

#include <Reve/RenderElement.h>
#include <Reve/NLTBases.h>

#include "TNamed.h"
#include "TAtt3D.h"
#include "TAttBBox.h"
#include "TColor.h"
#include "PODs.h"

class TBuffer3D;

namespace std {
  template<typename _Tp> class allocator;
  template<typename _Tp, typename _Alloc > class list;
}

namespace Reve {

class Vector;

class NLTPolygon
{
public:
  Int_t     fNPnts;
  Int_t*    fPnts; 

  NLTPolygon() : fNPnts(0), fPnts(0) {}
  NLTPolygon(Int_t n, Int_t* p) : fNPnts(n), fPnts(p) {}
  NLTPolygon(const NLTPolygon& x) : fNPnts(x.fNPnts), fPnts(x.fPnts) {}
  virtual ~NLTPolygon() {}

  NLTPolygon& operator=(const NLTPolygon& x)
  { fNPnts = x.fNPnts; fPnts = x.fPnts; return *this; }

  Int_t FindPoint(Int_t pi)
  { for (Int_t i=0; i<fNPnts; ++i) if (fPnts[i] == pi) return i; return -1; }

  ClassDef(NLTPolygon, 0)
};


class NLTPolygonSet :  public RenderElementList,
                       public NLTProjected,
		       public TAtt3D, 
                       public TAttBBox
{
  friend class NLTPolygonSetGL;
  friend class NLTPolygonSetEditor;
private:
  NLTPolygonSet(const NLTPolygonSet&);            // Not implemented
  NLTPolygonSet& operator=(const NLTPolygonSet&); // Not implemented

public:
  typedef std::list<Reve::NLTPolygon>     vpPolygon_t;
  typedef vpPolygon_t::iterator           vpPolygon_i;
  typedef vpPolygon_t::const_iterator     vpPolygon_ci;

private:
  TBuffer3D*   fBuff;
  Int_t*       fIdxMap; // map from original to projected and reduced point needed oly for geometry

  Bool_t       IsFirstIdxHead(Int_t s0, Int_t s1);
  void         AddPolygon(std::list<Int_t, std::allocator<Int_t> >& pp);

  void         ProjectAndReducePoints();
  void         MakePolygonsFromBP();
  void         MakePolygonsFromBS();
  void         ClearPolygonSet(); 

protected: 
  vpPolygon_t  fPols;   // NLT polygons

  Float_t      fEps;    // distance accounted in reducing the ponts
  Int_t        fNPnts;  // number of reduced and projected points
  Vector*      fPnts;   // reduced and projected points

  Color_t      fFillColor;  
  Color_t      fLineColor;
  Float_t      fLineWidth;

  UChar_t      fTransparency;

public:
  NLTPolygonSet(const Text_t* n="NLTPolygonSet", const Text_t* t="");
  virtual ~NLTPolygonSet();

  virtual void    SetProjection(NLTProjector* proj, NLTProjectable* model);
  virtual void    UpdateProjection();

  void            ProjectBuffer3D(); 

  virtual void    ComputeBBox();
  virtual void    Paint(Option_t* option = "");

  virtual void    DumpPolys() const; 
  void            DumpBuffer3D();

  //rendering
  virtual Bool_t  CanEditMainColor()        { return kTRUE; }
  virtual Color_t GetLineColor() const { return fLineColor; }
  
  virtual Bool_t  CanEditMainTransparency()      { return kTRUE; }
  virtual UChar_t GetMainTransparency() const    { return fTransparency; }
  virtual void    SetMainTransparency(UChar_t t) { fTransparency = t; }  

  virtual void    SetFillColor(Pixel_t pixel) { fFillColor = Color_t(TColor::GetColor(pixel));}
  virtual void    SetLineColor(Pixel_t pixel) { fLineColor = Color_t(TColor::GetColor(pixel));}

  virtual void    SetFillColor(Color_t c) { fFillColor = c; }
  virtual void    SetLineColor(Color_t c) { fLineColor = c; }
  virtual void    SetLineWidth(Double_t lw){fLineWidth = lw;}

  ClassDef(NLTPolygonSet,0) 

  }; // endclass NLTPolygonSet
} // namespace Reve

#endif
