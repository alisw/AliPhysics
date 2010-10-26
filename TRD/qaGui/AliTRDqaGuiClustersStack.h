#ifndef ALITRDQAGUICLUSTERSSTACK_H
#define ALITRDQAGUICLUSTERSSTACK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiClustersStack.h 23387 2008-01-17 17:25:16Z cblume $ */

#include "TGFrame.h"  

class TH1D;
class TString;
class TGLabel;
class TGComboBox;
class TGTextButton;
class TRootEmbeddedCanvas;

class AliTRDqaGuiClustersStack : public TGCompositeFrame {

 public:

  AliTRDqaGuiClustersStack();
  AliTRDqaGuiClustersStack(TGWindow *parent);
  ~AliTRDqaGuiClustersStack() {}

  
  void SetQAFile(const char *filename);
  void SetSM(Int_t idxSM);
  void SetStack(Int_t idxStack);
  void SetView(Int_t idxView);

  // void Play();  // *SLOT*
  void PreviusStack() {if (fIdxStack > 0) SetStack(fIdxStack-1);}   // *SLOT*
  void NextStack()    {if (fIdxStack < fgknStack-1) SetStack(fIdxStack+1);}   // *SLOT*

  void PreviusSM() {if (fIdxSM > 0) SetSM(fIdxSM-1);}      // *SLOT*
  void NextSM()    {if (fIdxSM < fgknSM-1) SetSM(fIdxSM+1);}   // *SLOT*
  
  void SelectStack(Int_t idx) {SetStack(idx);} // *SLOT*
  void SelectSM(Int_t idx) {SetSM(idx);}       // *SLOT*
  void SelectView(Int_t idx) {SetView(idx);}   // *SLOT*


 protected:

  static const Int_t fgknSM;     // number of Super Modules
  static const Int_t fgknStack;  // number of stacks
  static const Int_t fgknCh;     // number of chambers per stack

  Int_t fIdxSM;      // active super module
  Int_t fIdxStack;   // active stack
  Int_t fView;       // data type 

  //char fFileName[256];      // file with histograms
  const Char_t *fFileName;    // file with histograms
  
  TRootEmbeddedCanvas *fCanvasList[6];   // canvases
  TH1D *fHistList[6];                    // and histos
  
  // 
  TGCompositeFrame *fGPanel;   // panel with buttons
  TGCompositeFrame *fGCanvas;  // and with canvases
  
  // steering panel
  // TGLabel      *fGLabel;
  TGComboBox   *fGSelectSM;        // select super module
  TGComboBox   *fGSelectStack;     // select stack
  TGComboBox   *fGSelectView;      // select data type

  TGTextButton *fGPrevSM;       // button
  TGTextButton *fGPrevStack;    // button
  TGTextButton *fGNextSM;       // buton
  TGTextButton *fGNextStack;    // button
  TGTextButton *fGPlay;         // inactive button
 
  void CreateHistAmplitude();
  void CreateHistTimeCharge();
  void CreateHistTimeMPV();

 private:  
  AliTRDqaGuiClustersStack& operator = (const AliTRDqaGuiClustersStack& /*g*/) { return *this; };
  AliTRDqaGuiClustersStack(const AliTRDqaGuiClustersStack &);
  ClassDef(AliTRDqaGuiClustersStack,2) // 
};

#endif
