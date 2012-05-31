#ifndef ALITRDQAGUICLUSTERSSM_H
#define ALITRDQAGUICLUSTERSSM_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiClustersSM.h 23387 2008-01-17 17:25:16Z cblume $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of clusters on the full detector level.
// It displays histograms created by the AliTRDQADataMakerRec 
// run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"  

class TH1D;
class TString;
class TGLabel;
class TGComboBox;
class TGTextButton;
class TRootEmbeddedCanvas;

class AliTRDqaGuiClustersSM : public TGCompositeFrame {

 public:

  AliTRDqaGuiClustersSM();
  AliTRDqaGuiClustersSM(TGWindow *parent);
  ~AliTRDqaGuiClustersSM() {}
  
  void SetQAFile(const char *filename);
  void SetSM(Int_t idx);

  void Play();  // *SLOT*
  void PreviusSM() {if (fIdx > 0) SetSM(fIdx-1);}   // *SLOT*
  void NextSM()    {if (fIdx < 17) SetSM(fIdx+1);}  // *SLOT*
  void SelectSM(Int_t idx) {SetSM(idx);} // *SLOT*
  
 protected:

  Int_t fIdx;                          // super module index
  const char *fNameList[4];            // names of histograms
  static const Int_t fgkLogList[4];    // flag for logaritmic scale
   
  //char fFileName[256];                  // file with histograms
  const Char_t *fFileName;              // file with histograms
 
  TRootEmbeddedCanvas *fCanvasList[4];  // canvases
  TH1D *fHistList[4];                   // and histograms
  
  // 
  TGCompositeFrame *fGPanel;            // panel with buttons
  TGCompositeFrame *fGCanvas;           // and with histograms
  
  // steering panel
  // TGLabel      *fGLabel;
  TGComboBox   *fGSelect;    // sm selector button 
  TGTextButton *fGPrev;      // previus sm
  TGTextButton *fGNext;      // next sm
  TGTextButton *fGPlay;      // loop throu sm

 private:
  AliTRDqaGuiClustersSM& operator = (const AliTRDqaGuiClustersSM& /*g*/) { return *this; };
  AliTRDqaGuiClustersSM(const AliTRDqaGuiClustersSM&);

  ClassDef(AliTRDqaGuiClustersSM,2) // 
};

#endif
