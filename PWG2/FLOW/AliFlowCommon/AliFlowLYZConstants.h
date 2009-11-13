/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliFlowLYZConstants.h 25556 2008-05-02 08:22:51Z snelling $ */

#ifndef ALIFLOWLYZCONSTANTS_H
#define ALIFLOWLYZCONSTANTS_H

#include <TROOT.h>

// Description: constants for the LYZ flow makers

#include <TNamed.h>

class AliFlowLYZConstants : public TNamed {

public:
  AliFlowLYZConstants();
  virtual ~AliFlowLYZConstants();

  static AliFlowLYZConstants* GetMaster();

  const Int_t GetNtheta() const {return fNtheta;}
  const Int_t GetNbins() const {return fNbins;}
  const Double_t GetMaxSUM() const {return fMaxSUM;}
  const Double_t GetMaxPROD() const {return fMaxPROD;}

private:
  AliFlowLYZConstants& operator= (const AliFlowLYZConstants& c);
  AliFlowLYZConstants(const AliFlowLYZConstants& a);

  //data members
  Int_t fNtheta;     // number of reference angles theta
  Int_t fNbins;      // number of bins in fHistGtheta (AliFlowLYZHist1)
  Double_t  fMaxSUM;   // upper limit for fHistGtheta (AliFlowLYZHist1)
  Double_t  fMaxPROD;   // upper limit for fHistGtheta (AliFlowLYZHist1)

  static AliFlowLYZConstants* fgPMasterConfig;
  
  ClassDef(AliFlowLYZConstants,1)  // macro for rootcint 
};

#endif

