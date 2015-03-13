#ifndef AliAnalysisTaskElectronEfficiency_h
#define AliAnalysisTaskElectronEfficiency_h
/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//###########################################################
//#                                                         #
//#             Single Electron Efficiency Task             #
//#        and  Pair-Prefilter Efficiency Task              #
//#                                                         #
//#  Authors:                                               #
//#   Patrick Reichelt, Uni Ffm / Patrick.Reichelt@cern.ch  #
//#   Theo Broeker, Uni Ffm / Theo.Broeker@cern.ch          #
//#   Markus Koehler, GSI / M.Koehler@gsi.de                #
//#                                                         #
//###########################################################
/**
 Fills 3D histograms (mcPt, mcEta, mcPhi) for generated and reconstructed electrons.
 Ratios of these for given acceptance regions give 3D track efficiencies, to be then used in a pair efficiency generator.
 Cut instances are defined by adding AliAnalysisFilters via a Config file. Therefore also an LMEECutLib is supported.
 Additional functionality:
 Determination of random electron rejection efficiency due to pair-prefiltering (used for photon conversion + Dalitz rejection).
 It is estimated by pairing primary, non-injected, charged pions with the selected electrons (so the pair has no real correlation)
 and applying the prefilter pair cuts to these random pairs. All and rejected pions are stored in 3D histograms. 
 ---
 As examples to set up the task, see PWGDQ/dielectron/macrosLMEE/Config_reichelt_ElectronEfficiency.C or Config_tbroeker_ElectronEfficiency.C
**/


#include "AliAnalysisTaskSE.h"
#include "TTreeStream.h"

#include <vector>

#include "TParticlePDG.h"
#include "TDatabasePDG.h"

class AliPIDResponse;
class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDEvent;
class AliMCEvent;
class AliESDtrackCuts;
class AliDielectronPID;
class TString;

class AliAnalysisTaskElectronEfficiency : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskElectronEfficiency(const char *name = "AliAnalysisTaskElectronEfficiency");
  virtual ~AliAnalysisTaskElectronEfficiency();
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(const Option_t*);
  
  void          SetMC(Bool_t hasMC)                           {fIsMC = hasMC;}
  void          SetRequireVertex(Bool_t bReqVtx)              {fRequireVtx=bReqVtx;}
  void          SetMaxVertexZ(Double_t maxZ)                  {fMaxVtxZ=maxZ;}
  void          SetCentralityRange(Double_t min, Double_t max){fCentMin=min; fCentMax=max;}
  void          SetCheckV0daughterElectron(Bool_t bCheckV0)   {fCheckV0daughterElectron=bCheckV0;}
  void          SetCutInjectedSignal(Bool_t bCutInj)          {fCutInjectedSignal=bCutInj;}
  void          SetEtaRangeGEN(Double_t min, Double_t max)    {fEtaMinGEN=min; fEtaMaxGEN=max;}
  void          SetPtRangeGEN(Double_t min, Double_t max)     {fPtMinGEN=min; fPtMaxGEN=max;}
  void          SetSupportedCutInstance(Int_t supp)           {fSupportedCutInstance=supp;}
  void          SetWriteTree(Bool_t write)                    {fWriteTree=write;}
  void          SetPIDResponse(AliPIDResponse *fPIDRespIn)    {fPIDResponse=fPIDRespIn;}
  
  void          SetBins(Int_t Nptbins, Double_t *PtBins, Int_t Netabins, Double_t *EtaBins, Int_t Nphibins, Double_t *PhiBins) {
    /**/          fPtBins=PtBins;   fEtaBins=EtaBins;   fPhiBins=PhiBins;
    /**/          fNptBins=Nptbins; fNetaBins=Netabins; fNphiBins=Nphibins;
    /**/        }
  void          SetRunBins(TString runs)                      { fsRunBins=runs; }
  void          AttachTrackCuts(AliAnalysisFilter *cuts)      { fvTrackCuts.push_back(cuts); }
  void          AttachDoPrefilterEff(Bool_t doPref)           { fvDoPrefilterEff.push_back(doPref); }
  void          AttachRejCutMee(Double_t rejcut)              { fvRejCutMee.push_back(rejcut); }
  void          AttachRejCutTheta(Double_t rejcut)            { fvRejCutTheta.push_back(rejcut); }
  void          AttachRejCutPhiV(Double_t rejcut)             { fvRejCutPhiV.push_back(rejcut); }
  
  
  virtual void  CreateHistograms(TString names, Int_t cutInstance);
  void          CreateHistoGen();
  void          CreateSupportHistos();
  
  Int_t         GetNCutsets() const { return fvReco_Ele.size(); }
  //AliPIDResponse* GetPIDResponse() { return fPIDResponse; }
  
 private:
  Bool_t        IsInjectedSignal(AliMCEvent* mcEventLocal, Int_t tracklabel);
  void          CalcPrefilterEff(AliMCEvent* mcEventLocal, const std::vector< std::vector<Int_t> > & vvEleCand);
  Double_t      PhivPair(Double_t MagField, Int_t charge1, Int_t charge2, TVector3 fD1, TVector3 fD2);
  const char*   GetParticleName(Int_t pdg) {
    /**/          TParticlePDG* p1 = TDatabasePDG::Instance()->GetParticle(pdg);
    /**/          if(p1) return p1->GetName();
    /**/          return Form("%d", pdg);
    /**/        }
  
  AliESDEvent*      fESD;
  AliMCEvent*       mcEvent;
  AliPIDResponse*   fPIDResponse;
  Bool_t            fIsMC;
  Bool_t            fRequireVtx;
  Bool_t            fCheckV0daughterElectron;
  Bool_t            fCutInjectedSignal;
  Double_t          fMaxVtxZ;
  Double_t          fCentMin;
  Double_t          fCentMax;
  Double_t          fEtaMinGEN;
  Double_t          fEtaMaxGEN;
  Double_t          fPtMinGEN;
  Double_t          fPtMaxGEN;
  Int_t             fSupportedCutInstance;
  //Int_t             fEventcount; // for debugging
  
  Int_t             fNptBins;
  Int_t             fNetaBins;
  Int_t             fNphiBins;
  Double_t*         fPtBins;
  Double_t*         fEtaBins;
  Double_t*         fPhiBins;
  TString           fsRunBins; // for run dependency histogram
  
  //Cut Settings
  std::vector<AliAnalysisFilter*> fvTrackCuts;
  std::vector<Bool_t>             fvDoPrefilterEff;
  std::vector<Double_t>           fvRejCutMee;
  std::vector<Double_t>           fvRejCutTheta;
  std::vector<Double_t>           fvRejCutPhiV;
  //Efficiency Histograms
  TH3F*                           fNgen;
  std::vector<TH3F*>              fvReco_Ele;
  std::vector<TH3F*>              fvReco_Ele_poslabel; // store also result when using only positive label tracks.
  TH3F*                           fAllPions;
  std::vector<TH3F*>              fvPionsRejByAllSigns;
  std::vector<TH3F*>              fvPionsRejByUnlike;
  //std::vector<TH3F*>             fvReco_Pio; // be really careful if you need to implement this (see comments in UserExec).
  //std::vector<TH3F*>             fvReco_Kao; // be really careful if you need to implement this (see comments in UserExec).
  //std::vector<TH3F*>             fvReco_Pro; // be really careful if you need to implement this (see comments in UserExec).
  
  TList*                          fOutputList; // ! output data container
  TList*                          fOutputListSupportHistos; // ! output data container   
  
  //Output Tree with Tracks
  TTree*                          tracksT;
  Bool_t                          fWriteTree;
  
  //Track Variables
  // these are for data and MC:
  Float_t   pxESD;
  Float_t   pyESD;
  Float_t   pzESD;
  Float_t   pTPC;
  Int_t     chargeT;
  Float_t   signalITS;
  Float_t   signalTPC;
  Float_t   beta;
  Float_t   kchi2ITS;
  Int_t     kNclsITS;
  Float_t   kITSchi2Cl;//redundant
  Int_t     kNclsTPC;
  Float_t   kTPCchi2Cl;
  Int_t     kNclsTPCdEdx;
  Float_t   kNFclsTPCr;
  Float_t   kNFclsTPCfCross;
  Int_t     kNtrkltsTRD;
  Int_t     kNtrkltsTRDPID;
  Float_t   sigmaEleITS;
  Float_t   sigmaEleTPC;
  Float_t   sigmaEleTOF;
  //Float_t   sigmaEleTRD;
  Float_t   probEleTRD;
  Float_t   sigmaPioITS;
  Float_t   sigmaPioTPC;
  Float_t   sigmaPioTOF;
  Float_t   sigmaKaoITS;
  Float_t   sigmaKaoTPC;
  Float_t   sigmaProITS;
  Float_t   sigmaProTPC;
  Bool_t    isGlobalT;
  Bool_t    isGlobalSDD;
  // these are only for MC:
  Int_t     labelT;
  Int_t     pdgT;
  Int_t     labelmotherT;
  Int_t     pdgmotherT;
  Int_t     labelgrandmotherT;
  //Float_t   pMC = 0;
  Float_t   pxMC;
  Float_t   pyMC;
  Float_t   pzMC;
  UInt_t    selectedByCut; // bit mask
  
  AliAnalysisTaskElectronEfficiency(const AliAnalysisTaskElectronEfficiency&); // not implemented
  AliAnalysisTaskElectronEfficiency& operator=(const AliAnalysisTaskElectronEfficiency&); // not implemented
  
  ClassDef(AliAnalysisTaskElectronEfficiency, 1);
};

#endif
