#ifndef ALITRDQAGUIBLACKSM_H 
#define ALITRDQAGUIBLACKSM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiBlackSM.h 23387 2008-01-17 17:25:16Z cblume $ */

//////////////////////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"  

class TH1;
class TString;
class TGLabel;
class TGComboBox;
class TGTextButton;
class TRootEmbeddedCanvas;

class AliTRDqaGuiBlackSM : public TGCompositeFrame {
  
 public:

  AliTRDqaGuiBlackSM();
  AliTRDqaGuiBlackSM(TGWindow *parent);
  ~AliTRDqaGuiBlackSM() {}
  
  void SetQAFile(const char *filename);
  void SetSM(Int_t idx);

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

  //void Play();  // *SLOT*
  void PreviusSM() {if (fIdx > 0) SetSM(fIdx-1);}   // *SLOT*
  void NextSM()    {if (fIdx < 17) SetSM(fIdx+1);}  // *SLOT*
  void SelectSM(Int_t idx) {SetSM(idx);} // *SLOT*
  void SelectType(Int_t idx); // *SLOT*

 protected:

  Int_t fIdx;                   // SuperModule Index
  Int_t fIdxType;               // data type index
  const char *fNameList[5];     // list of possible data types

  Int_t    fSetRangePed;        // flag if use range for pedestals
  Double_t fRangePed[2];        // range for pedestals
  
  Int_t    fSetRangeNoise;      // flag if use range for noise 
  Double_t fRangeNoise[2];      // range for noise


  char fFileName[265];            // file with histograms
  
  TRootEmbeddedCanvas *fCanvasList[30];  // list of canvases
  TH1 *fHistList[30];                    // and histograms
  
  // 
  TGCompositeFrame *fGPanel;     // panel with buttons
  TGCompositeFrame *fGCanvas;    // and with canvases
  
  // steering panel
  // TGLabel      *fGLabel;
  TGComboBox   *fGSelect;       // sm selection
  TGTextButton *fGPrev;         // previus button
  TGTextButton *fGNext;         // next button
  //TGTextButton *fGPlay;

  TGComboBox *fGSelectType;     // data type selection

 private:
  AliTRDqaGuiBlackSM& operator = (const AliTRDqaGuiBlackSM& /*g*/) { return *this; };
  AliTRDqaGuiBlackSM(const AliTRDqaGuiBlackSM&);

  ClassDef(AliTRDqaGuiBlackSM,1) // 
};

#endif
