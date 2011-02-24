#ifndef ALIHFSYSTERR_H
#define ALIHFSYSTERR_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliRDHFSystErr
// to handle systematic errors for charm hadrons
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include <TNamed.h>
#include <TH1F.h>
#include "TGraphAsymmErrors.h"


class AliHFSystErr : public TNamed 
{
 public:

  AliHFSystErr(const Char_t* name="HFSystErr", const Char_t* title="");
  AliHFSystErr(Int_t decay,const Char_t* name="HFSystErr", const Char_t* title="");
  
  virtual ~AliHFSystErr();
  
  void DrawErrors(TGraphAsymmErrors *grErrFeeddown=0) const; 

  Double_t GetNormErr() const {return (fNorm ? fNorm->GetBinContent(0) : 0.);}
  Double_t GetBRErr() const {return (fBR ? fBR->GetBinContent(0) : 0.);}
  Double_t GetCutsEffErr(Double_t pt) const;
  Double_t GetMCPtShapeErr(Double_t pt) const;
  Double_t GetSeleEffErr(Double_t pt) const;
  Double_t GetPartAntipartErr(Double_t pt) const;
  Double_t GetPIDEffErr(Double_t pt) const;
  Double_t GetRawYieldErr(Double_t pt) const;
  Double_t GetTrackingEffErr(Double_t pt) const;
  Double_t GetTotalSystErr(Double_t pt,Double_t feeddownErr=0) const;


 private:

  AliHFSystErr(const AliHFSystErr& source);
  AliHFSystErr& operator=(const AliHFSystErr& source); 
 

  void InitD0toKpi();
  void InitDplustoKpipi();
  void InitDstartoD0pi();

  TH1F* ReflectHisto(TH1F *hin) const;

  TH1F *fNorm;            // normalization
  TH1F *fRawYield;        // raw yield 
  TH1F *fTrackingEff;     // tracking efficiency
  TH1F *fBR;              // branching ratio
  TH1F *fCutsEff;         // cuts efficiency
  TH1F *fPIDEff;          // PID efficiency
  TH1F *fMCPtShape;       // MC dNdpt
  TH1F *fPartAntipart;    // particle=antiparticle

  ClassDef(AliHFSystErr,1);  // class for systematic errors of charm hadrons
};

#endif

