#ifndef ALIHFSYSTERR_H
#define ALIHFSYSTERR_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFSystErr
/// \brief to handle systematic errors for charm hadrons
/// \author Author: A.Dainese, andrea.dainese@pd.infn.it
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
  Double_t GetDataDrivenFDErr(Double_t pt) const;
  Double_t GetRawYieldFDCorr(Double_t pt) const;
  Double_t GetTotalSystErr(Double_t pt,Double_t feeddownErr=0) const;

  void  ResetRawYieldErr(Double_t pt, Double_t val){
    fRawYield->SetBinContent(fRawYield->FindBin(pt),val);
  }
  void  ResetCutEfficErr(Double_t pt, Double_t val){
    fCutsEff->SetBinContent(fCutsEff->FindBin(pt),val);
  }
  void  ResetPIDEfficErr(Double_t pt, Double_t val){
    fPIDEff->SetBinContent(fPIDEff->FindBin(pt),val);
  }
  void  ResetMCPtShapeErr(Double_t pt, Double_t val){
    fMCPtShape->SetBinContent(fMCPtShape->FindBin(pt),val);
  }
  void  ResetTrackEfficErr(Double_t pt, Double_t val){
    fTrackingEff->SetBinContent(fTrackingEff->FindBin(pt),val);
  }
  /// Setting  the run number
  ///  set the two last numbers of the year (is 10 for 2010)
  void SetRunNumber(Int_t number) {
    fRunNumber = number;
    AliInfo(Form(" Settings for run year 20%2d",fRunNumber));
  }
  Int_t GetRunNumber() const { return fRunNumber; }
  /// Setting the collision type
  ///  0 is pp, 1 is PbPb, 2 is pPb
  void SetCollisionType(Int_t type) {
    fCollisionType = type;
    if (fCollisionType==0) { AliInfo(" Settings for p-p collisions"); }
    else if(fCollisionType==1) { AliInfo(" Settings for Pb-Pb collisions"); }
    else if(fCollisionType==2) { AliInfo(" Settings for p-Pb collisions"); }
  }
  Int_t GetCollisionType() const { return fCollisionType; }
  /// Setting for the centrality class
  ///  0100 for MB, 020 (4080) for 0-20 (40-80) CC and so on
  void SetCentrality(TString centrality) {
    fCentralityClass = centrality;
    AliInfo(Form(" Settings for centrality class %s",fCentralityClass.Data()));
  }
  void SetIsLowEnergy(Bool_t flag) {
    fIsLowEnergy = flag;
    if(flag) AliInfo(" Settings for the low energy run");
  }
  void SetIsLowPtAnalysis(Bool_t flag){
    fIsLowPtAnalysis = flag;
    if(flag) AliInfo("Settings for the low pt analysis");
  }
  void SetIsPass4Analysis(Bool_t flag){
    fIsPass4Analysis = flag;
    if(flag) AliInfo("Settings for the pass4 analysis");
  }
  void SetIs5TeVAnalysis(Bool_t flag){
    fIs5TeVAnalysis = flag;
    if(flag) AliInfo("Settings for the 5TeV analysis");
  }
  void SetStandardBins(Bool_t flag){
    fStandardBins= flag;
  }

  void SetIsBDTAnalysis(Bool_t flag){
    fIsBDTAnalysis = flag;
    if(flag) AliInfo("Settings for the Lc and Ds BDT analysis");
  }

  void SetIsMLAnalysis(Bool_t flag){
    fIsMLAnalysis = flag;
    if(flag) AliInfo("Settings for the Lc ML analysis");
  }

  void SetIsDataDrivenFDAnalysis(Bool_t flag){
    fIsDataDrivenFDAnalysis = flag;
    if(flag) AliInfo("Settings for the FD data-driven analyses");
  }

  void SetIsPbPb2010EnergyScan(Bool_t flag) {
    fIsCentScan = flag;
    if(flag) AliInfo(" Settings for the PbPb 2010 energy scan");
  }

  /// Settings of rapidity ranges for pPb 0-100% CC
  void SetRapidity(TString rapidity) {
    fRapidityRange = rapidity;
    AliInfo(Form(" Settings for rapidity interval %s",fRapidityRange.Data()));
  }
  void SetIspPb2011RapidityScan(Bool_t flag){
    fIsRapidityScan = flag;
    if(flag) AliInfo("Settings for the pPb vs y measurement");
  }

  /// Function to initialize the variables/histograms
  void Init(Int_t decay);

  void InitD0toKpi2010PbPb010CentScan();
  void InitD0toKpi2010PbPb1020CentScan();
  void InitD0toKpi2010PbPb2040CentScan();
  void InitD0toKpi2010PbPb4060CentScan();
  void InitD0toKpi2010PbPb6080CentScan();

  void InitD0toKpi2011PbPb3050InPlane();
  void InitD0toKpi2011PbPb3050OutOfPlane();

  void InitDplustoKpipi2010PbPb010CentScan();
  void InitDplustoKpipi2010PbPb1020CentScan();
  void InitDplustoKpipi2010PbPb2040CentScan();
  void InitDplustoKpipi2010PbPb4060CentScan();
  void InitDplustoKpipi2010PbPb6080CentScan();

  void InitDstartoD0pi2010PbPb010CentScan();
  void InitDstartoD0pi2010PbPb1020CentScan();
  void InitDstartoD0pi2010PbPb2040CentScan();
  void InitDstartoD0pi2010PbPb4060CentScan();
  void InitDstartoD0pi2010PbPb6080CentScan();

  void InitD0toKpi2011PbPb010CentScan();
  void InitD0toKpi2011PbPb1020CentScan();
  void InitD0toKpi2011PbPb2030CentScan();
  void InitD0toKpi2011PbPb3040CentScan();
  void InitD0toKpi2011PbPb4050CentScan();
  void InitD0toKpi2010PbPb5080CentScan();

  void InitDplustoKpipi2011PbPb010CentScan();
  void InitDplustoKpipi2011PbPb1020CentScan();
  void InitDplustoKpipi2011PbPb2030CentScan();
  void InitDplustoKpipi2011PbPb3040CentScan();
  void InitDplustoKpipi2011PbPb4050CentScan();
  void InitDplustoKpipi2010PbPb5080CentScan();

  void InitDstartoD0pi2011PbPb010CentScan();
  void InitDstartoD0pi2011PbPb1020CentScan();
  void InitDstartoD0pi2011PbPb2030CentScan();
  void InitDstartoD0pi2011PbPb3040CentScan();
  void InitDstartoD0pi2011PbPb4050CentScan();
  void InitDstartoD0pi2010PbPb5080CentScan();

  void InitD0toKpi2013pPb0100RapScan0804();
  void InitD0toKpi2013pPb0100RapScan0401();
  void InitD0toKpi2013pPb0100RapScan0101();
  void InitD0toKpi2013pPb0100RapScan0104();
  void InitD0toKpi2013pPb0100RapScan0408();

  void InitDplustoKpipi2013pPb0100RapScan0804();
  void InitDplustoKpipi2013pPb0100RapScan0401();
  void InitDplustoKpipi2013pPb0100RapScan0101();
  void InitDplustoKpipi2013pPb0100RapScan0104();
  void InitDplustoKpipi2013pPb0100RapScan0408();

  void InitDstartoD0pi2013pPb0100RapScan0804();
  void InitDstartoD0pi2013pPb0100RapScan0401();
  void InitDstartoD0pi2013pPb0100RapScan0101();
  void InitDstartoD0pi2013pPb0100RapScan0104();
  void InitDstartoD0pi2013pPb0100RapScan0408();

  void InitD0toKpi2013pPb020V0A();
  void InitD0toKpi2013pPb2040V0A();
  void InitD0toKpi2013pPb4060V0A();
  void InitD0toKpi2013pPb60100V0A();

  void InitD0toKpi2013pPb020ZNA();
  void InitD0toKpi2013pPb2040ZNA();
  void InitD0toKpi2013pPb4060ZNA();
  void InitD0toKpi2013pPb60100ZNA();

  void InitD0toKpi2016pPb010ZNA();
  void InitD0toKpi2016pPb1020ZNA();
  void InitD0toKpi2016pPb2040ZNA();
  void InitD0toKpi2016pPb4060ZNA();
  void InitD0toKpi2016pPb60100ZNA();

  void InitD0toKpi2013pPb020CL1();
  void InitD0toKpi2013pPb2040CL1();
  void InitD0toKpi2013pPb4060CL1();
  void InitD0toKpi2013pPb60100CL1();

  void InitDstartoD0pi2013pPb020V0A();
  void InitDstartoD0pi2013pPb2040V0A();
  void InitDstartoD0pi2013pPb4060V0A();
  void InitDstartoD0pi2013pPb60100V0A();

  void InitDstartoD0pi2013pPb020ZNA();
  void InitDstartoD0pi2013pPb2040ZNA();
  void InitDstartoD0pi2013pPb4060ZNA();
  void InitDstartoD0pi2013pPb60100ZNA();

   void InitDstartoD0pi2016pPb010ZNA();
   void InitDstartoD0pi2016pPb1020ZNA();
   void InitDstartoD0pi2016pPb2040ZNA();
   void InitDstartoD0pi2016pPb4060ZNA();
   void InitDstartoD0pi2016pPb60100ZNA();

  void InitDstartoD0pi2013pPb020CL1();
  void InitDstartoD0pi2013pPb2040CL1();
  void InitDstartoD0pi2013pPb4060CL1();
  void InitDstartoD0pi2013pPb60100CL1();

  void InitDplustoKpipi2013pPb020V0A();
  void InitDplustoKpipi2013pPb2040V0A();
  void InitDplustoKpipi2013pPb4060V0A();
  void InitDplustoKpipi2013pPb60100V0A();

  void InitDplustoKpipi2013pPb020ZNA();
  void InitDplustoKpipi2013pPb2040ZNA();
  void InitDplustoKpipi2013pPb4060ZNA();
  void InitDplustoKpipi2013pPb60100ZNA();

  void InitDplustoKpipi2013pPb020CL1();
  void InitDplustoKpipi2013pPb2040CL1();
  void InitDplustoKpipi2013pPb4060CL1();
  void InitDplustoKpipi2013pPb60100CL1();

  void InitDplustoKpipi2016pPb140trkl();
  void InitDplustoKpipi2016pPb4070trkl();
  void InitDplustoKpipi2016pPb70200trkl();
  void InitDplustoKpipi2016pPb010ZNA();
  void InitDplustoKpipi2016pPb1020ZNA();
  void InitDplustoKpipi2016pPb2040ZNA();
  void InitDplustoKpipi2016pPb4060ZNA();
  void InitDplustoKpipi2016pPb60100ZNA();

 private:

  AliHFSystErr(const AliHFSystErr& source);
  AliHFSystErr& operator=(const AliHFSystErr& source);

  void InitD0toKpi2010pp();
  void InitD0toKpi2010ppLowEn();
  void InitD0toKpi2010ppLowPtAn();
  void InitD0toKpi2010ppPass4();
  void InitD0toKpi2015pp5TeV();
  void InitD0toKpi2017pp5TeV();
  void InitD0toKpi2017pp5TeV_finebins();
  void InitDplustoKpipi2017pp5TeVML();
  void InitD0toKpi2017pp5TeVLowPtAn();
  void InitD0toKpi2017pp5TeVLowPtAn_finebins();
  void InitD0toKpi2016pp13TeV();
  void InitD0toKpi20161718pp13TeVmb();
  void InitD0toKpi20161718pp13TeVlm();
  void InitD0toKpi20161718pp13TeVhm();
  void InitD0toKpi20161718pp13TeVFineBins();
  void InitD0toKpi2011PbPb07half();
  void InitD0toKpi2010PbPb020();
  void InitD0toKpi2010PbPb4080();
  void InitD0toKpi2011PbPb3050();
  void InitD0toKpi2011PbPb010();
  void InitD0toKpi2013pPb0100();
  void InitD0toKpi2016pPb0100();
  void InitD0toKpi2016pPb5TeV_finebins();
  void InitD0toKpi2013pPb0100LowPtAn();
  void InitD0toKpi2016pPb0100LowPtAn();

  void InitDplustoKpipi2010pp();
  void InitDplustoKpipi2010ppPass4();
  void InitDplustoKpipi2010ppLowEn();
  void InitDplustoKpipi2012pp();
  void InitDplustoKpipi2015pp5TeV();
  void InitDplustoKpipi2017pp5TeV();
  void InitDplustoKpipi2017pp5TeV_finebins();
  void InitDplustoKpipi2016pp13TeV();
  void InitDplustoKpipi2011PbPb07half();
  void InitDplustoKpipi2010PbPb020();
  void InitDplustoKpipi2010PbPb4080();
  void InitDplustoKpipi2011PbPb3050();
  void InitDplustoKpipi2011PbPb010();
  void InitDplustoKpipi2013pPb0100();
  void InitDplustoKpipi2016pPb0100();
  void InitDplustoKpipi2016pPb5TeV_finebins();

  void InitDstartoD0pi2010pp();
  void InitDstartoD0pi2010ppLowEn();
  void InitDstartoD0pi2012pp();
  void InitDstartoD0pi2011PbPb07half();
  void InitDstartoD0pi2010PbPb020();
  void InitDstartoD0pi2010PbPb2040();
  void InitDstartoD0pi2010PbPb4080();
  void InitDstartoD0pi2011PbPb3050();
  void InitDstartoD0pi2011PbPb010();
  void InitDstartoD0pi2013pPb0100();
  void InitDstartoD0pi2016pPb0100();
  void InitDstartoD0pi2016pPb0100_fb();
  void InitDstartoD0pi2010ppPass4();
  void InitDstartoD0pi2017pp5TeV();
  void InitDstartoD0pi2017pp5TeV_finebins();
  void InitDstartoKpipi2016pp13TeV();

  void InitDstoKKpi2010pp();
  void InitDstoKKpi2010ppPass4();
  void InitDstoKKpi2017pp5TeV();
  void InitDstoKKpi2017pp5TeVBDT();
  void InitDstoKKpi2011PbPb07half();
  void InitDstoKKpi2011PbPb010();
  void InitDstoKKpi2011PbPb2050();
  void InitDstoKKpi2013pPb0100();
  void InitDstoKKpi2016pPb0100();
  void InitDstoKKpi2016pPb140trkl();
  void InitDstoKKpi2016pPb4070trkl();
  void InitDstoKKpi2016pPb70200trkl();
  void InitDstoKKpi2016pp13TeV();

  void InitLctopKpi2010pp();
  void InitLctopKpi2010ppBDT();
  void InitLctopKpi2013pPb();
  void InitLctopKpi2013pPbBDT();
  void InitLctopKpi2016pPb();
  void InitLctopKpi2017pp();
  void InitLctopKpi20161718pp13TeV();
  void InitLctopKpi20161718pp13TeVFineBins();

  void InitLctopK0S2010pp();
  void InitLctopK0S2013pPb();
  void InitLctopK0S2013pPbBDT();
  void InitLctopK0S2016pPb();
  void InitLctopK0S2016pPbBDT();
  void InitLctopK0S2017pp5TeV();
  void InitLctopK0S20161718pp13TeVBDT();

  void InitLctopK0S2018PbPb010BDT();
  void InitLctopK0S2018PbPb3050BDT();
  void InitLctopK0S2018PbPb010ML();
  void InitLctopK0S2018PbPb3050ML();
  void InitLctopK0S2018PbPb010();
  void InitLctopK0S2018PbPb3050();

  void InitD0toKpi2015PbPb010();
  void InitD0toKpi2015PbPb3050();
  void InitD0toKpi2015PbPb6080();

  void InitDplustoKpipi2015PbPb010();
  void InitDplustoKpipi2015PbPb3050();
  void InitDplustoKpipi2015PbPb6080();

  void InitDplustoKpipi2018PbPb010();
  void InitDplustoKpipi2018PbPb3050();

  void InitDstoKKpi2015PbPb010();
  void InitDstoKKpi2015PbPb3050();
  void InitDstoKKpi2015PbPb6080();
  void InitDstoKKpi2018PbPb010();
  void InitDstoKKpi2018PbPb010BDT();
  void InitDstoKKpi2018PbPb3050();
  void InitDstoKKpi2018PbPb3050BDT();

  void InitDstartoD0pi2015PbPb010();
  void InitDstartoD0pi2015PbPb3050();
  void InitDstartoD0pi2015PbPb6080();

  void InitDstartoD0pi2018PbPb010();
  void InitDstartoD0pi2018PbPb3050();

  void InitD0toKpi2018PbPb010();
  void InitD0toKpi2018PbPb3050();
  void InitD0toKpi2018PbPb010LowPtAn();

  void InitLctopKpiFromScpp13TeV201620172018(); // Lc(<-Sc)
  void InitScpp13TeV201620172018(); // Sc
  void InitLctopK0SFromScpp13TeV201620172018BDT(); // Lc(<-Sc), Lc->pK0S, BDT
  void InitScpp13TeV201620172018BDT(); // Sc, Lc->pK0S, BDT

  // data-driven non-prompt analyses
  void InitNonPromptDplustoKpipi2017pp5TeVML();
  void InitNonPromptDstoKKpi2017pp5TeVML();

  TH1F* ReflectHisto(TH1F *hin) const;

  TH1F *fNorm;            /// normalization
  TH1F *fRawYield;        /// raw yield
  TH1F *fTrackingEff;     /// tracking efficiency
  TH1F *fBR;              /// branching ratio
  TH1F *fCutsEff;         /// cuts efficiency
  TH1F *fPIDEff;          /// PID efficiency
  TH1F *fMCPtShape;       /// MC dNdpt
  TH1F *fPartAntipart;    /// particle=antiparticle
  TH1F *fDataDrivenFD;    /// prompt/FD fraction in case of data-driven analysis
  TH1F *fRawYieldFDCorr;  /// correlation between raw yield and prompt/FD fraction syst. unc. in case of data-driven analysis

  Int_t fRunNumber;        /// Run Number (year)
  Int_t fCollisionType;    /// Collision type: pp=0, PbPb=1
  TString fCentralityClass;  /// Centrality class
  /// MB:0100, 0-10:010, 0-20:020 ...40-80:4080...
  TString fRapidityRange;  /// Rapidity range fot y measurements

  Bool_t fIsLowEnergy;     /// flag for the low energy (2.76TeV) run
  Bool_t fIsLowPtAnalysis; /// flag for the low pt analysis (no topological cuts)
  Bool_t fIsPass4Analysis; /// flag for the pass4 analysis
  Bool_t fIs5TeVAnalysis; /// flag for the pp5TeV analysis
  Bool_t fIsBDTAnalysis;   /// flag for the Lc BDT analysis and Ds BDT analysis
  Bool_t fIsCentScan;      /// flag fot the PbPb centrality scan
  Bool_t fStandardBins;    /// flag for the standard bins in pp@5TeV and pPb@5TeV
  Bool_t fIsRapidityScan;  /// flag for the pPb vs y measurement
  Bool_t fIsMLAnalysis;   /// flag for the Lc ML analysis
  Bool_t fIsDataDrivenFDAnalysis;   /// flag for the non-prompt data-driven analyses 

  /// \cond CLASSIMP
  ClassDef(AliHFSystErr,14);  /// class for systematic errors of charm hadrons
  /// \endcond
};

#endif
