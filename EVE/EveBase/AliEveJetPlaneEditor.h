// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveJetPlaneEditor_H
#define AliEveJetPlaneEditor_H

#include <TGedFrame.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;
class TEveGValuator;

class AliEveJetPlane;

class AliEveJetPlaneEditor : public TGedFrame
{
public:
  AliEveJetPlaneEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                       UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveJetPlaneEditor() {}

  virtual void SetModel(TObject* obj);

  // Declare callback/slot methods
  // void DoXYZZ();
  void DoRnrJets();
  void DoRnrTracks();
  void DoArrowJetScale();
  void DoArrowTrackScale();
	void DoEnergyScale();
  void DoOneSelection();
  void DoTwoSelection();
  void DoStaticDataWindow();


  // --- Internal class for common settings

  class StaticDataWindow : public TGTransientFrame
  {
  public:
    StaticDataWindow(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
		     UInt_t options = kVerticalFrame);
    virtual ~StaticDataWindow();

    // slots
    void DoClose();
    void DoOK();
    void DoCancel();
    void DoTab(Int_t id);

  private:
    StaticDataWindow(const StaticDataWindow&);            // Not implemented
    StaticDataWindow& operator=(const StaticDataWindow&); // Not implemented

    TGCompositeFrame    *fFrame1, *fF2;             // Frames.
    TGButton            *fOkButton, *fCancelButton; // Ok, cancel buttons.
    TGLayoutHints       *fL1, *fL2, *fL3, *fL5;     // Layout hints.
    TGTab               *fTab;                      // Tab container.
    TGButton            *fChk1, *fChk2,*fChk3, *fChk4,*fChk5; // Check-buttons.

    ClassDef(StaticDataWindow, 0); // Common settings for all AliEveJetPlane objects.
  };

protected:
  AliEveJetPlane   *fM; // Model object.

  TGCheckButton    *fRnrJets;          // Widget for flag RnrJets.
  TGCheckButton    *fRnrTracks;        // Widget for flag RnrTracks.
  TEveGValuator    *fEnergyScale;      // Widget for EnergyScale.
  TEveGValuator    *fArrowJetScale;    // Widget for ArrowJetScale.
  TEveGValuator    *fArrowTrackScale;  // Widget for ArrowTrackScale.
  TGButton         *fOneSelection, *fTwoSelection;  // Widgets for one/two selection flags.
  TGButton         *fInformationSetup; // Widget for InformationSetup.

  static StaticDataWindow* fgStaticWindow; // Common window for global settings.

private:
  AliEveJetPlaneEditor(const AliEveJetPlaneEditor&);            // Not implemented
  AliEveJetPlaneEditor& operator=(const AliEveJetPlaneEditor&); // Not implemented

  ClassDef(AliEveJetPlaneEditor, 0); // Editor for AliEveJetPlane.
};

#endif
