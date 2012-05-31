#ifndef ALITRDQAGUIESDS_H
#define ALITRDQAGUIESDS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiESDs.h 23387 2008-01-17 17:25:16Z cblume $ */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
// of ESD (Event Summary Data)
// It displays histograms created by 
// the AliTRDQADataMakerRec run during the reconstruction 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"  

class TH1D;
class TRootEmbeddedCanvas;

class AliTRDqaGuiESDs : public TGCompositeFrame {

 public:

  AliTRDqaGuiESDs();
  AliTRDqaGuiESDs(TGWindow *parent, Int_t page);
  ~AliTRDqaGuiESDs() {};

  void SetPage(Int_t page) {fPage = page;}
  void SetQAFile(const char *filename);

 protected:

  Int_t fPage;                          // histogram set
  const char *fNameList[4*6];            // list of histograms
  static const Int_t fgkLogList[4*6];    // flag for log scale
 
  TRootEmbeddedCanvas *fCanvasList[6];  // canvas list
  TH1D *fHistList[6];                   // and histograms

 private:
  AliTRDqaGuiESDs& operator = (const AliTRDqaGuiESDs& /*g*/) { return *this; };
  AliTRDqaGuiESDs(const AliTRDqaGuiESDs &);

  ClassDef(AliTRDqaGuiESDs,1) // 
};

#endif
