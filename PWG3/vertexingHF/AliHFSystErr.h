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
#include <TString.h>
#include <TH1F.h>
#include "AliLog.h"
#include "TGraphAsymmErrors.h"


class AliHFSystErr : public TNamed 
{
 public:

  AliHFSystErr(const Char_t* name="HFSystErr", const Char_t* title="");
  
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

  // Setting  the run number
  //  set the two last numbers of the year (is 10 for 2010)
  void SetRunNumber(Int_t number) { 
    fRunNumber = number; 
    AliInfo(Form(" Settings for run year 20%2d",fRunNumber));
  }
  // Setting the collision type
  //  0 is pp, 1 is PbPb
  void SetCollisionType(Int_t type) { 
    fCollisionType = type; 
    if (fCollisionType==0) { AliInfo(" Settings for p-p collisions"); }
    else if(fCollisionType==1) { AliInfo(" Settings for Pb-Pb collisions"); }
  }
  // Setting for the centrality class
  //  0100 for MB, 020 (4080) for 0-20 (40-80) CC and so on
  void SetCentrality(TString centrality) { 
    fCentralityClass = centrality; 
    AliInfo(Form(" Settings for centrality class %s",fCentralityClass.Data()));
  }
  void SetIsLowEnergy(Bool_t flag) { 
    fIsLowEnergy = flag; 
    AliInfo(" Settings for the low energy run");
  }

  // Function to initialize the variables/histograms
  void Init(Int_t decay);

 private:

  AliHFSystErr(const AliHFSystErr& source);
  AliHFSystErr& operator=(const AliHFSystErr& source); 
 
  void InitD0toKpi2010pp();
  void InitDplustoKpipi2010pp();
  void InitDstartoD0pi2010pp();
  void InitDstoKKpi2010pp();

  void InitD0toKpi2010PbPb020();
  void InitDplustoKpipi2010PbPb020();
  void InitDstartoD0pi2010PbPb020();

  void InitD0toKpi2010PbPb4080();
  void InitDplustoKpipi2010PbPb4080();
  void InitDstartoD0pi2010PbPb4080();

  void InitD0toKpi2010ppLowEn();
  void InitDplustoKpipi2010ppLowEn();
 void InitDstartoD0pi2010ppLowEn();

 void InitDstartoD0pi2010PbPb2040();

  TH1F* ReflectHisto(TH1F *hin) const;

  TH1F *fNorm;            // normalization
  TH1F *fRawYield;        // raw yield 
  TH1F *fTrackingEff;     // tracking efficiency
  TH1F *fBR;              // branching ratio
  TH1F *fCutsEff;         // cuts efficiency
  TH1F *fPIDEff;          // PID efficiency
  TH1F *fMCPtShape;       // MC dNdpt
  TH1F *fPartAntipart;    // particle=antiparticle

  Int_t fRunNumber;        // Run Number (year)
  Int_t fCollisionType;    // Collision type: pp=0, PbPb=1
  TString fCentralityClass;  // Centrality class
                           // MB:0100, 0-10:010, 0-20:020 ...40-80:4080...
  Bool_t fIsLowEnergy;     // flag for the low energy (2.76TeV) run

  ClassDef(AliHFSystErr,3);  // class for systematic errors of charm hadrons
};

#endif

