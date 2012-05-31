#ifndef ALITRDQAGUICLUSTERS_H
#define ALITRDQAGUICLUSTERS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiClusters.h 23387 2008-01-17 17:25:16Z cblume $ */

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
class TH2D;
class TRootEmbeddedCanvas;

class AliTRDqaGuiClusters : public TGCompositeFrame {
  
 public:

  AliTRDqaGuiClusters();
  AliTRDqaGuiClusters(TGWindow *parent);
  ~AliTRDqaGuiClusters() {};
  
  void SetQAFile(const char *filename);

 protected:

  const char *fNameList[4];           // list of histogram names
  static const Int_t fgkLogList[4];   // flag for log scale
 
  TRootEmbeddedCanvas *fCanvasList[4]; // list of canvases
  TH1D *fHistList[4];                  // and histograms
  TH1D *fHistRefs[4][3];               // histograms with colors

  // TH2D *BuildHisto(TH1D *ref, TH1D *data);
  void BuildColor(Int_t i, TH1D *ref);

 private:
  AliTRDqaGuiClusters& operator = (const AliTRDqaGuiClusters& /*g*/) { return *this; };
  AliTRDqaGuiClusters(const AliTRDqaGuiClusters &);

  ClassDef(AliTRDqaGuiClusters,1) // 
};

#endif
