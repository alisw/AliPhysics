#ifndef ALITRDQAGUIJPSISTEPS_H
#define ALITRDQAGUIJPSISTEPS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

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

class AliTRDqaGuiJPsi : public TGCompositeFrame {

 public:

  AliTRDqaGuiJPsi();
  AliTRDqaGuiJPsi(TGWindow *parent);
  ~AliTRDqaGuiJPsi() {}
  
  void SetQAFile(const char *filename);
  void SetStep(Int_t idx);

  void PreviusStep()         {if (fIdx > 0) SetStep(fIdx-1);}   // *SLOT*
  void NextStep()            {if (fIdx < 17) SetStep(fIdx+1);}  // *SLOT*
  void SelectStep(Int_t idx) {SetStep(idx);} // *SLOT*
  
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
  AliTRDqaGuiJPsi& operator = (const AliTRDqaGuiJPsi& /*g*/) { return *this; };
  AliTRDqaGuiJPsi(const AliTRDqaGuiJPsi&);

  ClassDef(AliTRDqaGuiJPsi,2) // 
};

#endif
