// $Header$

#ifndef REVE_TrackEditors_H
#define REVE_TrackEditors_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TGComboBox;
class TGLineWidthComboBox;

namespace Reve {

class RGValuator;
class RGDoubleValuator;

class TrackList;

/**************************************************************************/
// TrackListEditor
/**************************************************************************/

class TrackListEditor : public TGedFrame
{
  TrackListEditor(const TrackListEditor&);            // Not implemented
  TrackListEditor& operator=(const TrackListEditor&); // Not implemented

protected:
  TrackList* fTC; // fModel dynamic-casted to TrackListEditor

  Reve::RGValuator*  fMaxR;
  Reve::RGValuator*  fMaxZ;
  TGNumberEntry*     fMaxOrbits;
  TGNumberEntry*     fMinAng;
  TGNumberEntry*     fDelta;

  TGLineWidthComboBox* fWidthCombo;

  TGCheckButton*     fRnrTracks;
  TGCheckButton*     fRnrMarkers;

  TGCheckButton*     fFitDaughters;
  TGCheckButton*     fFitDecay;

  RGDoubleValuator*  fPtRange;

public:
  TrackListEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		  UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  ~TrackListEditor();

  virtual void SetModel(TObject* obj);

  void DoMaxR();
  void DoMaxZ();
  void DoMaxOrbits();
  void DoMinAng();
  void DoDelta();

  void DoLineWidth(Int_t width);

  void DoRnrTracks();
  void DoRnrMarkers();

  void DoFitDaughters();
  void DoFitDecay();

  void DoPtRange();

  ClassDef(TrackListEditor, 1); // Editor for TrackList
}; // endclass TrackListEditor



/**************************************************************************/
// 
/**************************************************************************/

class TrackCounter;

class TrackCounterEditor : public TGedFrame
{
private:
  TrackCounterEditor(const TrackCounterEditor&);            // Not implemented
  TrackCounterEditor& operator=(const TrackCounterEditor&); // Not implemented

protected:
  TrackCounter* fM; // fModel dynamic-casted to TrackCounter

  // Declare widgets
  TGComboBox*    fClickAction;
  TGLabel*       fInfoLabel;
  TGNumberEntry* fEventId;

public:
  TrackCounterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		     UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TrackCounterEditor();

  virtual void SetModel(TObject* obj);

  void DoOrtoXY();
  void DoOrtoZY();
  void DoPersp();

  void DoPrev();
  void DoNext();
  void DoSetEvent();

  void DoPrintReport();
  void DoFileReport();
  void DoShowHistos();

  void DoClickAction(Int_t);

  ClassDef(TrackCounterEditor, 1); // Editor for TrackCounter
}; // endclass TrackCounterEditor

}

#endif
