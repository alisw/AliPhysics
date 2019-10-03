/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/
//
//
// VZERO event plane task for 2010
// Gain equalization + Recentering 
// Need a root file
//
//
#ifndef ALIHFEVZEROEVENTPLANE_H
#define ALIHFEVZEROEVENTPLANE_H

#include "TProfile.h"
#include <AliVEvent.h>
#include "TString.h"
//#include "TH1F.h"
//#include "TList.h"


class AliHFEVZEROEventPlane : public TNamed {
 public:
  AliHFEVZEROEventPlane();
  AliHFEVZEROEventPlane(const char *name, const Char_t *title);
  AliHFEVZEROEventPlane(const AliHFEVZEROEventPlane &ref);
  AliHFEVZEROEventPlane& operator=(const AliHFEVZEROEventPlane &ref);
  virtual void Copy(TObject &o) const;
  ~AliHFEVZEROEventPlane();

  void ProcessEvent(AliVEvent *event);
  
  void  SetNameFile(TString namefile) {fnamefile = namefile;};
  Bool_t OpenInfoCalbration(Int_t run);

  Double_t GetEventPlaneV0A() const {return fEventPlaneV0A;};
  Double_t GetEventPlaneV0C() const {return fEventPlaneV0C;};
  Double_t GetEventPlaneV0()  const {return fEventPlaneV0;};

  TList *GetOutputList() const {return fOutputList;};

 private:
  virtual void Analyze(AliVEvent* esdEvent); 

  Double_t fEventPlaneV0A;          // Corrected event plane V0A
  Double_t fEventPlaneV0C;          // Corrected event plane V0C
  Double_t fEventPlaneV0;           // Corrected event plane V0

  AliVEvent* fESD;                  //! ESD object
  Int_t        fRun;                // Run number
  TProfile *fMultV0;                //! fMultiplicityV0
  Float_t fV0Cpol,fV0Apol;          // fV0Cpol, fV0Apol
  static const Int_t fgknCentrBin = 9; // Centrality bins
  Float_t fMeanQ[fgknCentrBin][2][2];  // mean for centering
  Float_t fWidthQ[fgknCentrBin][2][2]; // rms for centering
  TString fnamefile;                // name of the file with the coefficient
  TList       *fOutputList;         //! Output list
  TProfile    *fMultV0Before;       //! fMultiplicityV0 Before
  TProfile    *fMultV0After;        //! fMultiplicityV0 After
  TH1F *fQBefore[fgknCentrBin][2][2]; //! Q centering before
  TH1F *fQAfter[fgknCentrBin][2][2];  //! Q centering after
 

  ClassDef(AliHFEVZEROEventPlane, 1);    //Analysis task for high pt analysis 
};

#endif
