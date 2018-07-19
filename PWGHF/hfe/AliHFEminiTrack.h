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
// Originator:  M.Fasel <M.Fasel@gsi.de>
// Base class for AliHFEminiEventCreator.cxx/h
// Contributors: Nirbhay K. Behera, Jiyeon Kwon, Jonghan Park

#ifndef ALIHFEMINITRACK_H
#define ALIHFEMINITRACK_H

#include <TObject.h>
#include <TMath.h>
#include <TBits.h>

class AliHFEminiTrack : public TObject{
 public:
  AliHFEminiTrack();
  AliHFEminiTrack(const AliHFEminiTrack &ref);
  AliHFEminiTrack &operator=(const AliHFEminiTrack &ref);
  ~AliHFEminiTrack() {}
  
  // -------------- Getters ------------------------
  Double_t Pt() const { return TMath::Abs(fSignedPt); }
  Int_t Charge() const {
    if(fSignedPt >= 0.) return 1.;
    return -1;
  }
  
  Double_t MCPt() const { return TMath::Abs(fMCSignedPt); }
  Int_t    MCPDG() const { return fMCPDG; }
  Int_t    MCMotherPdg() const { return fMCMotherPdg; }
  Int_t    MCSource() const { return fMCSource; }
  Int_t    MCElectronSource() const { return static_cast<Int_t>(fMCEleSource); }
  Double_t MCElectronSourcePt() const { return fMCEleSourcePt; }
  Double_t HFEImpactParameter() const { return fHFEImpactParam[0]; }
  Double_t HFEImpactParameterResolution() const { return fHFEImpactParam[1]; }
  

  // -------------- Setters ------------------------
  void SetSignedPt(Double_t abspt, Bool_t positivecharge) { 
    Double_t charge = positivecharge ? 1. : -1;
    fSignedPt = abspt * charge;
  }
  
  void SetMCSignedPt(Double_t abspt, Bool_t positivecharge){
    Double_t charge = positivecharge ? 1. : -1;
    fMCSignedPt = abspt * charge;
  }
  void SetMCPDG(Int_t mcpdg) {fMCPDG = mcpdg; }
  void SetMCMotherPdg(Int_t pdg) { fMCMotherPdg = pdg; }
  void SetMCSource(Int_t mcsource) { fMCSource = mcsource; }
  void SetMCElectronSource(Int_t source) { fMCEleSource = static_cast<UChar_t>(source); }
  void SetMCElectronSourcePt(Double_t sourcept) { fMCEleSourcePt = sourcept; }
  void SetHFEImpactParam(Double_t impactParam, Double_t impactParamResolution){
    fHFEImpactParam[0] = impactParam;
    fHFEImpactParam[1] = impactParamResolution;
  }
  
 private:
  Double_t fSignedPt;                     // signed pt
  Double_t fMCSignedPt;                   // MCSignedPt
  Int_t    fMCPDG;                        // MCPDG
  Int_t    fMCMotherPdg;                  // MCMP
  Int_t    fMCGrandMotherPdg;             // MCGMP
  Int_t    fMCSource;                     // MCSource
  UChar_t  fMCEleSource;                  // MC Electron Source (AliHFEmcQA)
  Double_t fMCEleSourcePt;                // MC Electron Source's pt (AliHFEmcQA)
  Double_t fHFEImpactParam[2];            // HFE impact paramter (value, resolution) for beauty analysis

  ClassDef(AliHFEminiTrack, 4)
};
#endif
