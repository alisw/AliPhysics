#ifndef ALITRDQAGUIBLACKCHAMBER_H
#define ALITRDQAGUIBLACKCHAMBER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiBlackChamber.h 23387 2008-01-17 17:25:16Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                    
// This class is a Graphical User Interface for the Quality Monitorig 
// of black (non zero zuppresed) events from TRD. 
// It lets display and browse throu histograms created by the class 
// AliTRDqaBlackEvents.
// The class works in cooperation with AliTRDqaGuiMainBlack.
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
//                                                                    
////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"  

class TH1;
class TString;
class TGLabel;
class TGComboBox;
class TGTextButton;
class TRootEmbeddedCanvas;

class AliTRDqaGuiBlackChamber : public TGCompositeFrame {

 public:

  AliTRDqaGuiBlackChamber();
  AliTRDqaGuiBlackChamber(TGWindow *parent);
  ~AliTRDqaGuiBlackChamber() {};
  
  void SetQAFile(const char *filename);
  void SetSM(Int_t idxSM);
  void SetChamber(Int_t idxChamber);
  void SetView(Int_t idxView);

  void SetRangePed(Double_t min, Double_t max) {
    fSetRangePed = 1;
    fRangePed[0] = min;
    fRangePed[1] = max;
  }

  void SetRangeNoise(Double_t min, Double_t max) {
    fSetRangeNoise = 1;
    fRangeNoise[0] = min;
    fRangeNoise[1] = max;
  }


  // void Play();  // *SLOT*
  void PreviusChamber() {if (fIdxChamber > 0) SetChamber(fIdxChamber-1);}   // *SLOT*
  void NextChamber()    {if (fIdxChamber < fgknChamber-1) SetChamber(fIdxChamber+1);}   // *SLOT*

  void PreviusSM() {if (fIdxSM > 0) SetSM(fIdxSM-1);}      // *SLOT*
  void NextSM()    {if (fIdxSM < fgknSM-1) SetSM(fIdxSM+1);}   // *SLOT*
  
  void SelectChamber(Int_t idx) {SetChamber(idx);} // *SLOT*
  void SelectSM(Int_t idx) {SetSM(idx);}       // *SLOT*
  void SelectView(Int_t idx) {SetView(idx);}   // *SLOT*


 protected:
  
  Int_t fView;
  
  static const Int_t fgknSM;       // number of supermodules
  static const Int_t fgknChamber;  // number of chamberd (30)

  Int_t    fSetRangePed;           // flag for range in pedestals
  Double_t fRangePed[2];           // range for pedelstals 
  
  Int_t    fSetRangeNoise;         // flag for range in noise
  Double_t fRangeNoise[2];         // range in noise
  
  Int_t fIdxSM;                    // active super module
  Int_t fIdxChamber;               // active chamber 
  //Int_t fView;

  //char fFileName[256];                  // file with histograms
  const Char_t *fFileName;              // file with histograms
  
  TRootEmbeddedCanvas *fCanvasList[5];  // canvases
  TH1 *fHistList[5];                    // and histos
  
  // 
  TGCompositeFrame *fGPanel;            // panel with buttons
  TGCompositeFrame *fGCanvas;           // canvas
  TGCompositeFrame *fGCanvasUp;         // canvas
  TGCompositeFrame *fGCanvasDown;       // canvas
  
  // steering panel
  // TGLabel      *fGLabel;
  TGComboBox   *fGSelectSM;             // selector for Super Module
  TGComboBox   *fGSelectChamber;        // selector for Chamber 
  TGComboBox   *fGSelectView;           // select view

  TGTextButton *fGPrevSM;        // button
  TGTextButton *fGPrevChamber;   // button
  TGTextButton *fGNextSM;        // button
  TGTextButton *fGNextChamber;   // button
  //TGTextButton *fGPlay;
 
  //void CreateHistAmplitude();
  // void CreateHistTimeCharge();
  //void CreateHistTimeMPV();

 private:
  AliTRDqaGuiBlackChamber& operator = (const AliTRDqaGuiBlackChamber& /*g*/) { return *this; };
  AliTRDqaGuiBlackChamber(const AliTRDqaGuiBlackChamber &);

  ClassDef(AliTRDqaGuiBlackChamber,1) // 

};

#endif
