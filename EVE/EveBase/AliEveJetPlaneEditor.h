// $Id$
// Main authors: Matevz Tadel & Alja Mrak-Tadel: 2006, 2007

/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef ALIEVE_JetPlaneEditor_H
#define ALIEVE_JetPlaneEditor_H

#include <TGedFrame.h>
#include <RQ_OBJECT.h>

class TGButton;
class TGCheckButton;
class TGNumberEntry;
class TGColorSelect;

class TEveGValuator;


// class TEveGValuator;

class AliEveJetPlane;

class AliEveJetPlaneEditor : public TGedFrame
{
private:
  AliEveJetPlaneEditor(const AliEveJetPlaneEditor&);            // Not implemented
  AliEveJetPlaneEditor& operator=(const AliEveJetPlaneEditor&); // Not implemented

protected:
  AliEveJetPlane   *fM; // Model object.

  TGCheckButton    *fRnrJets;          // Widget for flag RnrJets.
  TGCheckButton    *fRnrTracks;        // Widget for flag RnrTracks.
  TEveGValuator    *fEnergyScale;      // Widget for EnergyScale.
  TEveGValuator    *fEnergyColorScale; // Widget for EnergyColorScale.
  TGButton         *fOneSelection, *fTwoSelection;  // Widgets for one/two selection flags.
  TGButton         *fInformationSetup; // Widget for InformationSetup.

public:
  AliEveJetPlaneEditor(const TGWindow* p=0, Int_t width=170, Int_t height=30,
                       UInt_t options=kChildFrame, Pixel_t back=GetDefaultFrameBackground());
  virtual ~AliEveJetPlaneEditor() {}

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
    StaticDataWindow(const StaticDataWindow&);            // Not implemented
    StaticDataWindow& operator=(const StaticDataWindow&); // Not implemented

    TGCompositeFrame    *fFrame1, *fF2;             // Frames.
    TGButton            *fOkButton, *fCancelButton; // Ok, cancel buttons.
    TGLayoutHints       *fL1, *fL2, *fL3, *fL5;     // Layout hints.
    TGTab               *fTab;                      // Tab container.
    TGButton            *fChk1, *fChk2,*fChk3, *fChk4,*fChk5; // Check-buttons.

  public:
    StaticDataWindow(const TGWindow *p, const TGWindow *main, UInt_t w, UInt_t h,
		     UInt_t options = kVerticalFrame);
    virtual ~StaticDataWindow();

    // slots
    void DoClose();
    void DoOK();
    void DoCancel();
    void DoTab(Int_t id);

    ClassDef(StaticDataWindow, 0); // Common settings for all AliEveJetPlane objects.
  };

protected:
  static StaticDataWindow* fgStaticWindow; // Common window for global settings.

  ClassDef(AliEveJetPlaneEditor, 1); // Editor for AliEveJetPlane.
}; // endclass AliEveJetPlaneEditor

#endif
