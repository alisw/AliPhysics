#ifndef ALITRDQAGUIMAINANALYSIS_H
#define ALITRDQAGUIMAINANALYSIS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////////////////////////////
//
// This class is a Graphical User Interface for the Quality Monitorig 
//
// S. Radomski 
// Uni-Heidelberg
// Feb. 2008
// 
//////////////////////////////////////////////////////////////////////////////////

#include "TGFrame.h"

class TGWindow;
class TGTab;

class AliTRDqaGuiJPsi;
class AliTRDqaGuiEnergyDeposit;

class AliTRDqaGuiMainAnalysis : public TGMainFrame {
 
 public:

  AliTRDqaGuiMainAnalysis();
  AliTRDqaGuiMainAnalysis(TGWindow *parent);
  ~AliTRDqaGuiMainAnalysis() {}

  void SetQAFile();
  
 protected:
  
  TGTab *fGTabPanel;                  // main tab panel
  AliTRDqaGuiJPsi          *fGJPsi;       // 
  AliTRDqaGuiEnergyDeposit *fGED;

 private:
  AliTRDqaGuiMainAnalysis& operator = (const AliTRDqaGuiMainAnalysis& /*g*/) { return *this; };
  AliTRDqaGuiMainAnalysis (const AliTRDqaGuiMainAnalysis&);

  ClassDef(AliTRDqaGuiMainAnalysis,1) //
};

#endif
