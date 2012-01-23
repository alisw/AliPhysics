#ifndef ALIANALYSISTASKOMEGAPI0PIPI_H
#define ALIANALYSISTASKOMEGAPI0PIPI_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//---------------------------------------------------------------------------// 
// omega(782) -> pi0 pi+ pi- analysis.                                       //
//---------------------------------------------------------------------------//

#include "AliAnalysisTaskSE.h"

class TH1F;
class TH2F;
class TString;

class AliAnalysisTaskOmegaPi0PiPi : public AliAnalysisTaskSE
{
public:

  AliAnalysisTaskOmegaPi0PiPi();
  AliAnalysisTaskOmegaPi0PiPi(const char* name);
  virtual ~AliAnalysisTaskOmegaPi0PiPi();

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t * opt);
  
  void AnalyzeModules(const char* modules="123") { fModules = modules; }
  TString GetAnalyzedModules() const { return fModules; }
  
private:

  AliAnalysisTaskOmegaPi0PiPi(const AliAnalysisTaskOmegaPi0PiPi&); 
  AliAnalysisTaskOmegaPi0PiPi& operator=(const AliAnalysisTaskOmegaPi0PiPi&); 

private:
  
  TList* fOutputContainer; // List of output histograms
  TH1F*  fhM2piSel;        // V0 inv. mass, Dxy cut applied
  TH1F*  fhDxy;            // dxy of V0s
  TH1F*  fhMggSel;         // two-cluster inv. mass
  TH1F*  fhImpXY;          // track XY-impact parameters

  TH2F*  fhM3pi0to2;   // M(3pi) vs pT(3pi), 0 < ptrack < 2 GeV
  TH2F*  fhM3pi2to4;   // M(3pi) vs pT(3pi), 2 < ptrack < 4 GeV
  TH2F*  fhM3pi4to6;   // M(3pi) vs pT(3pi), 4 < ptrack < 6 GeV
  TH2F*  fhM3pi6to8;   // M(3pi) vs pT(3pi), 6 < ptrack < 8 GeV
  TH2F*  fhM3pi8to10;  // M(3pi) vs pT(3pi), 8 < ptrack < 10 GeV
  TH2F*  fhM3pi10to12; // M(3pi) vs pT(3pi), 10< ptrack < 12 GeV
  TH2F*  fhM3pi12;     // M(3pi) vs pT(3pi), ptrack > 12 GeV
  TH2F*  fhM3pi0to8;   // M(3pi) vs pT(3pi), 0 < ptrack < 8 GeV

  TString fModules;

  ClassDef(AliAnalysisTaskOmegaPi0PiPi,2);
  
};

#endif
