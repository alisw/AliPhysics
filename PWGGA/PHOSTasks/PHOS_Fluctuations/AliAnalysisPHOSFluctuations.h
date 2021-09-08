#ifndef AliAnalysisPHOSFluctuations_cxx
#define AliAnalysisPHOSFluctuations_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for for gamma-hadron fluctuation study
// Authors: Dmitri Peresunko, Evgenia Nekrasova
// 02-Sept-2021

class AliPPVsMultUtils ;
class AliPIDResponse;
class AliAnalysisUtils;
class AliCaloPhoton;
//#include "AliPPVsMultUtils.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisPHOSFluctuations : public AliAnalysisTaskSE {

public:

  enum etaCut{kPhosEta,kTpcEta} ;
  enum phiCut{kPhosPhi,kFullPhi} ;
  enum runtype{kpp,kpPb,kPbPb} ;

  AliAnalysisPHOSFluctuations(const char *name = "AliAnalysisPHOSFluctuations");
  virtual ~AliAnalysisPHOSFluctuations() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void SetRunType(runtype a){fRunType=a;  } 
  void SetChargedCut(etaCut ec, phiCut pc){
    if(ec == kPhosEta)fChEtaCutMax=0.15;
    if(ec == kTpcEta)fChEtaCutMax=0.8; 
    if(pc == kPhosPhi){fChPhiMin=4.1887902; fChPhiMax=5.5850536;}
    if(pc == kFullPhi){fChPhiMin=0.; fChPhiMax=TMath::TwoPi();}
  }
  void SetChargedPtCut(float ptMin, float ptMax){fChPtMin = ptMin; fChPtMax = ptMax; }
  void SetPhotonCut(etaCut ec, phiCut pc){
    if(ec == kPhosEta)fPhEtaCutMax=0.15;
    if(ec == kTpcEta) fPhEtaCutMax=0.8; 
    if(pc == kPhosPhi){fPhPhiMin=4.1887902; fPhPhiMax=5.5850536;}
    if(pc == kFullPhi){fPhPhiMin=0.; fPhPhiMax=TMath::TwoPi();}
  }
  void SetPhotonPtCut(float ptMin, float ptMax){fPhPtMin = ptMin; fPhPtMax = ptMax; }
  void SetPi0Cut(etaCut ec, phiCut pc){
    if(ec == kPhosEta)fPi0EtaCutMax=0.15;
    if(ec == kTpcEta) fPi0EtaCutMax=0.8; 
    if(pc == kPhosPhi){fPi0PhiMin=4.1887902; fPi0PhiMax=5.5850536;}
    if(pc == kFullPhi){fPi0PhiMin=0.; fPi0PhiMax=TMath::TwoPi();}
  }
  void SetPi0PtCut(float ptMin, float ptMax){fPi0PtMin = ptMin; fPi0PtMax = ptMax; }

  
private:
  AliAnalysisPHOSFluctuations(const AliAnalysisPHOSFluctuations&); // not implemented
  AliAnalysisPHOSFluctuations& operator=(const AliAnalysisPHOSFluctuations&); // not implemented

protected: 
  void SelectPhotons() ;
  void SelectTracks() ;
  void ProcessMC() ; //Analyze MC info
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  int  CommonAnsestor(int prim1, int prim2);

protected:
  static constexpr int NCENT = 20; 
  static constexpr int NPID = 4 ; 

  float fPhPtMin = 0.3 ; 
  float fPhPtMax = 1.  ;
  float fPi0PtMin = 0.3 ; 
  float fPi0PtMax = 1.  ;
  float fChPtMin = 0.3 ; 
  float fChPtMax = 1.  ;

  float fChEtaCutMax = 0.8 ; 
  float fPhEtaCutMax = 0.8 ; 
  float fPi0EtaCutMax = 0.8 ; 

  float fChPhiMin = 0. ; 
  float fChPhiMax = 6.2831853  ;
  float fPhPhiMin = 0. ; 
  float fPhPhiMax = 6.2831853  ;
  float fPi0PhiMin = 0. ; 
  float fPi0PhiMax = 6.2831853  ;

  Int_t fRunNumber =0 ;
  double fCentrality =0.; //
  int fCentBin =0 ;
  runtype   fRunType = kpp;
  AliPPVsMultUtils * fPPUtils = nullptr;
  THashList * fOutputContainer = nullptr;         //final histogram container 
  TList   * fCurrentMixedList = nullptr;          //! list of previous evetns for given centrality
  TList   * fPHOSEvents[20][NCENT] ;    //!Previous events for mixing

  TClonesArray       *fStack = nullptr;           //! Pointer to MC stack (AOD)
  AliPIDResponse     *fPIDResponse  = nullptr;    //! PID response 
  AliAnalysisUtils   *fUtils = nullptr;           //!
  TClonesArray       *fArrGamma  = nullptr ;      //! list of photons in event

  std::vector<int> fpi0list ;
  std::vector<int> fChPrimaryList1 ;
  std::vector<int> fChPrimaryList2 ;
  std::vector<int> fPhPrimaryList ;
  std::vector<int> fPi0PrimaryList ;
  std::vector<int> frecpi0 ;

  int fRecPhot[NPID] ={0} ;
  int fRecPhotTrue[NPID]={0} ;
  int fRecPhotTruePi0[NPID]={0} ;
  int fRecPhotTruePi0Single[NPID]={0} ;

  int fRecPipm = 0;
  int fRecPipmTrue = 0;

  //QA histos
  TH1F * fhSelEvents = nullptr; 
  TH1F * fhTrackPt = nullptr ;
  TH1F * fhPionPt = nullptr ;
  TH1F * fhKaonPt = nullptr ;
  TH1F * fhProtonPt = nullptr ;
  TH1F * fhUndefPt = nullptr ;
  TH1F * fhPhotonPt = nullptr ;
  TH1F * fhMCPhotonPt = nullptr ;

  //MC histos
  TH1F * fhMCPrimPi0N = nullptr ;
  TH1F * fhMCPrimPi01N  = nullptr;
  TH1F * fhMCPrimPi0NoresN = nullptr ;
  TH1F * fhMCPrimPi0Nores1N = nullptr ;

  TH1F * fhMCPrimGammaN = nullptr ;
  TH1F * fhMCPrimGamma1N = nullptr;
  TH1F * fhMCPrimGammaPi0N = nullptr ;
  TH1F * fhMCPrimGammaPi01N = nullptr ;
  TH1F * fhMCPrimGammaPi0SingleN = nullptr ;
  TH1F * fhMCPrimGammaPi0Single1N = nullptr ;
  TH1F * fhMCPrimGammaAllSingleN = nullptr ;
  TH1F * fhMCPrimGammaAllSingle1N = nullptr ;
  TH1F * fhMCPrimGammaPi0SingleNoresN = nullptr ;
  TH1F * fhMCPrimGammaPi0SingleNores1N = nullptr ;

  TH1F * fhMCPrimPipmN = nullptr ;
  TH1F * fhMCPrimPipm1N  = nullptr;
  TH1F * fhMCPrimPipmNoresNa = nullptr ;
  TH1F * fhMCPrimPipmNores1Na  = nullptr;
  TH1F * fhMCPrimPipmNoresNb = nullptr ;
  TH1F * fhMCPrimPipmNores1Nb  = nullptr;
  TH1F * fhMCPrimPipmPi0 = nullptr ;
  TH1F * fhMCPrimPipmPi0Nores = nullptr ;
  TH1F * fhMCPrimPipmGamma = nullptr ;
  TH1F * fhMCPrimPipmGammaPi0  = nullptr;
  TH1F * fhMCPrimPipmGammaPi0Single  = nullptr;
  TH1F * fhMCPrimPipmGammaAllSingle  = nullptr;
  TH1F * fhMCPrimPipmGammaPi0SingleNores = nullptr ;

  TH1F * fhPipmN = nullptr ;
  TH1F * fhPipm1N  = nullptr;
  TH1F * fhPipmTrueN = nullptr ;
  TH1F * fhPipmTrue1N = nullptr ;

  TH1F * fhGammaN[NPID] = {nullptr} ;
  TH1F * fhGamma1N[NPID] = {nullptr} ;
  TH1F * fhGammaTrueN[NPID] = {nullptr} ;
  TH1F * fhGammaTrue1N[NPID] = {nullptr} ;
  TH1F * fhGammaTruePi0N[NPID] = {nullptr} ;
  TH1F * fhGammaTruePi01N[NPID] = {nullptr} ;
  TH1F * fhGammaTruePi0SingleN[NPID] = {nullptr} ;
  TH1F * fhGammaTruePi0Single1N[NPID] = {nullptr} ;

  TH1F * fhGammaPipm[NPID] = {nullptr} ;
  TH1F * fhGammaPipmTrue[NPID] = {nullptr} ;
  TH1F * fhGammaPipmTruePi0[NPID] = {nullptr} ;
  TH1F * fhGammaPipmTruePi0Single[NPID] = {nullptr} ;

  //Inv masses
  TH2F * fhReal[NPID] = {nullptr} ;
  TH2F * fhRealTrue[NPID]  = {nullptr};
  TH2F * fhRealCommon[NPID]  = {nullptr};
  TH2F * fhMixed[NPID] = {nullptr} ;
  TH2F * fEgammaEpi0 = nullptr;

  ClassDef(AliAnalysisPHOSFluctuations, 1); // PHOS analysis task
};

#endif
