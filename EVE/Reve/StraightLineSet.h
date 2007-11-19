// $Header$

#ifndef REVE_StraightLineSet_H
#define REVE_StraightLineSet_H

#include <Reve/Reve.h>

#include <Gtypes.h>
#include <TNamed.h>
#include <TQObject.h>
#include <TAtt3D.h>
#include <TAttMarker.h>
#include <TAttLine.h>
#include <TAttBBox.h>

#include <Reve/Reve.h>
#include <Reve/RenderElement.h>
#include <Reve/NLTBases.h>
#include <Reve/Plex.h>
#include <Reve/ZTrans.h>
class TRandom;

namespace Reve {

class StraightLineSet : public RenderElement,
                        public NLTProjectable,
		        public TNamed, 
                        public TQObject,
		        public TAtt3D,
                        public TAttLine,
                        public TAttMarker,
		        public TAttBBox
{
private:
  StraightLineSet(const StraightLineSet&);            // Not implemented
  StraightLineSet& operator=(const StraightLineSet&); // Not implemented

public:
  struct Line
  {
    Float_t        fV1[3];
    Float_t        fV2[3];
    TRef           fRef;

    Line(Float_t x1, Float_t y1, Float_t z1,Float_t x2, Float_t y2, Float_t z2)
    {
      fV1[0] = x1, fV1[1] = y1, fV1[2] = z1;
      fV2[0] = x2, fV2[1] = y2, fV2[2] = z2;
    }
  };
 
  struct Marker
  {
    Int_t   	   fLineID;
    Float_t 	   fPos;
    TRef           fRef;

    Marker(Int_t lineID, Float_t pos) : fLineID(lineID), fPos(pos) {};
  };

protected:
  VoidCPlex         fLinePlex;
  VoidCPlex         fMarkerPlex;

  Bool_t            fOwnLinesIds;       //Flag specifying if id-objects are owned by the QuadSet
  Bool_t            fOwnMarkersIds;     //Flag specifying if id-objects are owned by the QuadSet

  Bool_t            fRnrMarkers;
  Bool_t            fRnrLines;

  Line*             fLastLine; //!          

  Bool_t            fTrans;
  ZTrans            fHMTrans;
public:
  StraightLineSet(const Text_t* n="StraightLine", const Text_t* t="");
  virtual ~StraightLineSet() {}

  virtual Bool_t  CanEditMainHMTrans() { return  kTRUE; }
  virtual ZTrans* PtrMainHMTrans()     { return &fHMTrans; }

  ZTrans& RefHMTrans() { return fHMTrans; }
  void SetTransMatrix(Double_t* carr)        { fHMTrans.SetFrom(carr); }
  void SetTransMatrix(const TGeoMatrix& mat) { fHMTrans.SetFrom(mat);  }


  virtual void SetLineColor(Color_t col) { SetMainColor(col); }


  void    AddLine(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
  void    AddMarker(Int_t lineID, Float_t pos);

  VoidCPlex& GetLinePlex()   { return fLinePlex;   }
  VoidCPlex& GetMarkerPlex() { return fMarkerPlex; }

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option="");

  virtual void SetRnrMarkers(Bool_t x) {fRnrMarkers = x;}
  virtual Bool_t GetRnrMarkers(){return fRnrMarkers;}

  virtual void SetRnrLines(Bool_t x) {fRnrLines = x;}
  virtual Bool_t GetRnrLines(){return fRnrLines;}

  virtual TClass* ProjectedClass() const;

  ClassDef(StraightLineSet, 1); // Set of lines and optional markers.
}; // endclass StraightLineSet



/**************************************************************************/

class NLTSLineSet : public StraightLineSet,
                    public NLTProjected
{
private:
  NLTSLineSet(const NLTSLineSet&);            // Not implemented
  NLTSLineSet& operator=(const NLTSLineSet&); // Not implemented

protected:

public:
  NLTSLineSet();
  virtual ~NLTSLineSet() {}

  virtual Bool_t  CanEditMainHMTrans() { return  kFALSE; }

  virtual void SetProjection(NLTProjector* proj, NLTProjectable* model);

  virtual void UpdateProjection();

  ClassDef(NLTSLineSet, 1); // NLT projected StraightLineSet.
}; // endclass NLTSLineSet

} // Reve

#endif
