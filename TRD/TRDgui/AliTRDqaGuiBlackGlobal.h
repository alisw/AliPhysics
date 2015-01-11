#ifndef ALITRDQAGUIBLACKGLOBAL_H 
#define ALITRDQAGUIBLACKGLOBAL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTRDqaGuiBlackGlobal.h 23387 2008-01-17 17:25:16Z cblume $ */

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
class TGraph;
class TString;
class TRootEmbeddedCanvas;

class AliTRDqaGuiBlackGlobal : public TGCompositeFrame {
  
 public:

  AliTRDqaGuiBlackGlobal();
  AliTRDqaGuiBlackGlobal(TGWindow *parent);
  ~AliTRDqaGuiBlackGlobal() {}
  
  void SetQAFile(const char *filename);
 
 protected:

  char fFileName[265];                    // file with histograms
  
  TRootEmbeddedCanvas *fCanvasList[6];    // list of canvases
  TH1    *fHistList[3];                   // and histograms
  TGraph *fGraphList[3];                  // trend graphs

 private:

  AliTRDqaGuiBlackGlobal& operator = (const AliTRDqaGuiBlackGlobal& /*g*/) { return *this; };
  AliTRDqaGuiBlackGlobal(const AliTRDqaGuiBlackGlobal&);

  ClassDef(AliTRDqaGuiBlackGlobal,2)      // Gui class for black events 

};

#endif
