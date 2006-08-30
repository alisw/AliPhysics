// $Header$

#ifndef REVE_TrackEditors_H
#define REVE_TrackEditors_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Reve {

class RGValuator;
class RGDoubleValuator;

class TrackList;

class TrackListEditor : public TGedFrame
{
  TrackListEditor(const TrackListEditor&);            // Not implemented
  TrackListEditor& operator=(const TrackListEditor&); // Not implemented

protected:
  TrackList* fTC; // fModel dynamic-casted to TrackListEditor

  TGNumberEntry*     fMaxR;
  TGNumberEntry*     fMaxZ;
  TGNumberEntry*     fMaxOrbits;
  TGNumberEntry*     fMinAng;
  TGNumberEntry*     fDelta;

  TGCheckButton*     fRnrTracks;
  TGCheckButton*     fRnrMarkers;

  TGCheckButton*     fFitDaughters;
  TGCheckButton*     fFitDecay;

  RGDoubleValuator*  fPtRange;

public:
  TrackListEditor(const TGWindow* p, Int_t id, Int_t width = 170, Int_t height = 30, UInt_t options = kChildFrame, Pixel_t back = GetDefaultFrameBackground());
  ~TrackListEditor();

  virtual void SetModel(TVirtualPad* pad, TObject* obj, Int_t event);

  void DoMaxR();
  void DoMaxZ();
  void DoMaxOrbits();
  void DoMinAng();
  void DoDelta();

  void DoRnrTracks();
  void DoRnrMarkers();

  void DoFitDaughters();
  void DoFitDecay();

  void DoPtRange();

  ClassDef(TrackListEditor, 1); // Editor for TrackList
}; // endclass TrackListEditor

}

#endif
