// $Header$

#ifndef REVE_StraightLineSet_H
#define REVE_StraightLineSet_H

#include <Reve/Reve.h>

#include <Gtypes.h>
#include <TNamed.h>
#include <TQObject.h>
#include <TAtt3D.h>
#include <TAttMarker.h>
#include <TAttBBox.h>

#include <Reve/Reve.h>
#include <Reve/RenderElement.h>
#include <Reve/Plex.h>
class TRandom;

namespace Reve {

class StraightLineSet : public RenderElement,
		        public TNamed, public TQObject,
		        public TAtt3D,
                        public TAttMarker,
		        public TAttBBox
{
  friend class StraightLineSetGL; 
  friend class StraightLineSetEditor;

private:
  StraightLineSet(const StraightLineSet&);            // Not implemented
  StraightLineSet& operator=(const StraightLineSet&); // Not implemented

protected:
  struct Line {
    Float_t        fV1[3];
    Float_t        fV2[3];
    TRef           fRef;

    Line(Float_t x1, Float_t y1, Float_t z1,Float_t x2, Float_t y2, Float_t z2)
    {
      fV1[0] = x1, fV1[1] = y1, fV1[2] = z1;
      fV2[0] = x2, fV2[1] = y2, fV2[2] = z2;
    }
  };
 
  struct Marker {
    Int_t   	   fLineID;
    Float_t 	   fPos;
    TRef           fRef;

    Marker(Int_t lineID, Float_t pos) : fLineID(lineID), fPos(pos) {};
  };

protected:
  Color_t           fColor;

  Bool_t            fOwnLinesIds;       //Flag specifying if id-objects are owned by the QuadSet
  Bool_t            fOwnMarkersIds;       //Flag specifying if id-objects are owned by the QuadSet
  VoidCPlex         fLinePlex;
  Line*             fLastLine;     //!
  VoidCPlex         fMarkerPlex;

  Bool_t            fRnrMarkers;
  Bool_t            fRnrLines;

public:
  StraightLineSet(const Text_t* n="StraightLine", const Text_t* t="");
  virtual ~StraightLineSet() {}

  virtual Bool_t CanEditMainColor() { return kTRUE; }

  void AddLine(Float_t x1, Float_t y1, Float_t z1, Float_t x2, Float_t y2, Float_t z2);
  void AddMarker(Int_t lineID, Float_t pos);

  virtual void ComputeBBox();
  virtual void Paint(Option_t* option="");

  Color_t GetColor(){ return fColor; };
  void    SetColor(Color_t c){ fColor=c; }

  ClassDef(StraightLineSet, 1);
}; // endclass StraightLineSet

}

#endif
