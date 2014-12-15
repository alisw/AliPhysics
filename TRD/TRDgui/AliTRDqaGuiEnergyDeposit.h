#ifndef ALITRDQAGUIENERGYDEPOSIT_H
#define ALITRDQAGUIENERGYDEPOSIT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
//
// S. Radomski 
// Uni-Heidelberg
// April 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"  

class TH1D;
class TString;
class TGLabel;
class TGComboBox;
class TGTextButton;
class TRootEmbeddedCanvas;

class AliTRDqaGuiEnergyDeposit : public TGCompositeFrame {

 public:

  AliTRDqaGuiEnergyDeposit();
  AliTRDqaGuiEnergyDeposit(TGWindow *parent);
  ~AliTRDqaGuiEnergyDeposit() {}
  
  void SetQAFile(const char *filename);
  void SetType(Int_t idx);

  void PreviusType()         {if (fIdx > 0) SetType(fIdx-1);}   // *SLOT*
  void NextType()            {if (fIdx < 17) SetType(fIdx+1);}  // *SLOT*
  void SelectType(Int_t idx) {SetType(idx);} // *SLOT*
  
 protected:

  Int_t fIdx;                          // super module index
  const char *fNameList[6];            // names of histograms

  //char fFileName[256];                // file with histograms
  const Char_t *fFileName;              // file with histograms
 
  TRootEmbeddedCanvas *fCanvasList[6];  // canvases
  TH1D *fHistList[6];                   // and histograms
  
  // 
  TGCompositeFrame *fGPanel;            // panel with buttons
  TGCompositeFrame *fGCanvas;           // and with histograms
  
  // steering panel
  TGComboBox   *fGSelect;    // step selector button 
  TGTextButton *fGPrev;      // previus step
  TGTextButton *fGNext;      // next step

 private:
  AliTRDqaGuiEnergyDeposit& operator = (const AliTRDqaGuiEnergyDeposit& /*g*/) { return *this; };
  AliTRDqaGuiEnergyDeposit(const AliTRDqaGuiEnergyDeposit&);

  ClassDef(AliTRDqaGuiEnergyDeposit,2) // 
};

#endif
