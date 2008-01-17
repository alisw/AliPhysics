// $Header$

#ifndef ALIEVE_JetPlaneEditor_H
#define ALIEVE_JetPlaneEditor_H

#include <TGedFrame.h>
#include <RQ_OBJECT.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TEveGValuator;

namespace Alieve {

// class TEveGValuator;

class JetPlane;

class JetPlaneEditor : public TGedFrame
{
private:
  JetPlaneEditor(const JetPlaneEditor&);            // Not implemented
  JetPlaneEditor& operator=(const JetPlaneEditor&); // Not implemented

protected:
  JetPlane          *fM; // fModel dynamic-casted to JetPlaneEditor

  // Declare widgets
  // TGSomeWidget*   fXYZZ;
  TGCheckButton*     fRnrJets;
  TGCheckButton*     fRnrTracks;
  TEveGValuator*  fEnergyScale;
  TEveGValuator*  fEnergyColorScale;
  TGButton          *fOneSelection, *fTwoSelection;
  TGButton          *fInformationSetup;

public:
  JetPlaneEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
		 UInt_t options = kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~JetPlaneEditor();

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DoRnrJets();
  void DoRnrTracks();
  void DoEnergyColorScale();
  void DoEnergyScale();
  void DoOneSelection();
  void DoTwoSelection();
  void DoStaticDataWindow();


  // --- Internal class for common settings
public:
  class StaticDataWindow : public TGTransientFrame
  {
  private:
    TGCompositeFrame    *fFrame1, *fF2;
    TGButton            *fOkButton, *fCancelButton;
    TGLayoutHints       *fL1, *fL2, *fL3, *fL5;
    TGTab               *fTab;
    TGButton            *fChk1, *fChk2,*fChk3, *fChk4,*fChk5;

  public:
    StaticDataWindow(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
		     UInt_t options = kVerticalFrame);
    virtual ~StaticDataWindow();

    // slots
    void DoClose();
    void DoOK();
    void DoCancel();
    void DoTab(Int_t id);

    ClassDef(StaticDataWindow, 0);
  };

protected:
  static StaticDataWindow* fgStaticWindow;

  ClassDef(JetPlaneEditor, 1); // Editor for JetPlane
}; // endclass JetPlaneEditor

}

#endif
