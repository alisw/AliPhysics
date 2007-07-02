#ifndef REVE_NLTPolygonSet_H
#define REVE_NLTPolygonSet_H

#include <Reve/RenderElement.h>

#include "TNamed.h"
#include "TAtt3D.h"
#include "TAttBBox.h"
#include "TColor.h"

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
  NLTPolygon& operator=(const NLTPolygon& x)
  { fNPnts = x.fNPnts; fPnts = x.fPnts; return *this; }

  virtual ~NLTPolygon() {}

  ClassDef(NLTPolygon, 0)
};


class NLTPolygonSet :  public Reve::RenderElement,
                       public TNamed, 
		       public TAtt3D, 
                       public TAttBBox
{
  friend class NLTPolygonSetGL;
  friend class NLTPolygonSetEditor;

  NLTPolygonSet(const NLTPolygonSet&);            // Not implemented
  NLTPolygonSet& operator=(const NLTPolygonSet&); // Not implemented

protected:
  Int_t        fNPnts;
  Vector*      fPnts;

  Int_t        fNPols; // number of polygons
  NLTPolygon*  fPols;  // vector of polygon structs

  Color_t      fFillColor;  
  Color_t      fLineColor;
  Float_t      fLineWidth;
  Float_t      fZDepth;

public:
  NLTPolygonSet(const Text_t* n="NLTPolygonSet", const Text_t* t="");
  virtual ~NLTPolygonSet();

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option = "");

  void SetPoints(Vector* p, Int_t n) {fPnts = p; fNPnts = n;}
  void SetPolygons(NLTPolygon* p, Int_t n){fPols=p; fNPols=n;}

  virtual Color_t GetFillColor() const { return fFillColor; }
  virtual Color_t GetLineColor() const { return fLineColor; }
  
  Int_t GetSize(){return fNPols;}
  virtual void Dump() const;

  virtual void SetFillColor(Pixel_t pixel) { fFillColor = Color_t(TColor::GetColor(pixel));}
  virtual void SetLineColor(Pixel_t pixel) { fLineColor = Color_t(TColor::GetColor(pixel));}

  virtual void SetFillColor(Color_t c) { fFillColor = c; }
  virtual void SetLineColor(Color_t c) { fLineColor = c; }

  virtual void SetZDepth(Float_t z){fZDepth = z;}
  virtual void SetLineWidth(Double_t lw){fLineWidth = lw;}


  ClassDef(NLTPolygonSet,0) 
    }; // endclass NLTPolygonSet
} // namespace Reve

#endif
