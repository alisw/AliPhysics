// $Header$

#ifndef ALIEVE_TrackFitterEditor_H
#define ALIEVE_TrackFitterEditor_H

#include <TGedFrame.h>

class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

namespace Alieve {

class TrackFitter;

class TrackFitterEditor : public TGedFrame
{
private:
  TrackFitterEditor(const TrackFitterEditor&);            // Not implemented
  TrackFitterEditor& operator=(const TrackFitterEditor&); // Not implemented

protected:
  TrackFitter* fM; // fModel dynamic-casted to TrackFitterEditor

  TGTextButton* fFit;   // button to fit selection
  TGTextButton* fReset; // button to reset selection
  TGTextButton* fStart; // button to connect to signal
  TGTextButton* fStop;  // button to disconnect from signal
  TGTextButton* fGraph; // button to draw graph

public:
  TrackFitterEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30, UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~TrackFitterEditor() {}

  virtual void SetModel(TObject* obj);

  void DoStart();
  void DoFit();
  void DoReset(); 
  void DoStop();
  void DoGraph();

  ClassDef(TrackFitterEditor, 0); // Editor for TrackFitter class.
}; // endclass TrackFitterEditor

}

#endif
