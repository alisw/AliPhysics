#ifndef ALIANALYSISPHOSFLUCTUATIONS_H
#define ALIANALYSISPHOSFLUCTUATIONS_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

// Analysis task for for gamma-hadron fluctuation study
// Authors: Dmitri Peresunko, Evgenia Nekrasova
// 02-Sept-2021

class AliPPVsMultUtils ;
class AliPIDResponse;
class AliAnalysisUtils;
class AliCaloPhoton;

#include "AliAnalysisTaskSE.h"

class AliAnalysisPHOSFluctuations : public AliAnalysisTaskSE {

public:

  enum etaCut{kPhosEta,kTpcEta} ;
  enum phiCut{kPhosPhi,kFullPhi} ;
  enum runtype{kpp,kpPb,kPbPb} ;

  AliAnalysisPHOSFluctuations(){} 
  AliAnalysisPHOSFluctuations(const char *name);
  AliAnalysisPHOSFluctuations(const AliAnalysisPHOSFluctuations&){} // not implemented
  AliAnalysisPHOSFluctuations& operator=(const AliAnalysisPHOSFluctuations&ap){   // not implemented
    this->~AliAnalysisPHOSFluctuations();
    new(this) AliAnalysisPHOSFluctuations(ap);
    return *this;
  }
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

protected: 
  void SelectPhotons() ;
  void SelectTracks() ;
  void ProcessMC() ; //Analyze MC info
  void FillHistogram(const char * key,Double_t x) const ; //Fill 1D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y) const ; //Fill 2D histogram witn name key
  void FillHistogram(const char * key,Double_t x, Double_t y, Double_t z) const ; //Fill 3D histogram witn name key
  int  CommonAnsestor(int prim1, int prim2);

protected:
  static const int NCENT = 20; 
  static const int NPID = 4 ; 

  float fPhPtMin  ; 
  float fPhPtMax ;
  float fPi0PtMin  ; 
  float fPi0PtMax   ;
  float fChPtMin  ; 
  float fChPtMax   ;

  float fChEtaCutMax  ; 
  float fPhEtaCutMax  ; 
  float fPi0EtaCutMax  ; 

  float fChPhiMin ; 
  float fChPhiMax ;
  float fPhPhiMin ; 
  float fPhPhiMax ;
  float fPi0PhiMin ; 
  float fPi0PhiMax ;

  Int_t fRunNumber ;                           //!
  double fCentrality ;                         //!
  int fCentBin ;
  runtype   fRunType ;
  AliPPVsMultUtils * fPPUtils ;          //!
  THashList * fOutputContainer ;         //! final histogram container 
  TList   * fCurrentMixedList ;          //! list of previous evetns for given centrality
  TList   * fPHOSEvents[20][NCENT] ;              //!Previous events for mixing

  TClonesArray       *fStack ;           //! Pointer to MC stack (AOD)
  AliPIDResponse     *fPIDResponse ;    //! PID response 
  AliAnalysisUtils   *fUtils ;           //!
  TClonesArray       *fArrGamma ;      //! list of photons in event

  std::vector<int> fpi0list ;               //!
  std::vector<int> fChPrimaryList1 ;        //!
  std::vector<int> fChPrimaryList2 ;        //!
  std::vector<int> fPhPrimaryList ;         //! 
  std::vector<int> fPi0PrimaryList ;        //! 
  std::vector<int> frecpi0 ;                //!

  int fRecPhot[NPID] ;                 //!
  int fRecPhotTrue[NPID] ;              //!  
  int fRecPhotTruePi0[NPID];           //!
  int fRecPhotTruePi0Single[NPID];     //!

  int fRecPipm ;
  int fRecPipmTrue ;

  //QA histos
  TH1F * fhSelEvents ;     //!
  TH1F * fhTrackPt ;      //!
  TH1F * fhPionPt  ;       //!
  TH1F * fhKaonPt  ;       //!
  TH1F * fhProtonPt  ;     //!
  TH1F * fhUndefPt  ;      //!
  TH1F * fhPhotonPt  ;     //!
  TH1F * fhMCPhotonPt  ;   //!

  //MC histos
  TH1F * fhMCPrimPi0N  ;        //!
  TH1F * fhMCPrimPi01N  ;       //!
  TH1F * fhMCPrimPi0NoresN  ;   //!
  TH1F * fhMCPrimPi0Nores1N  ;  //!

  TH1F * fhMCPrimGammaN  ;      //!
  TH1F * fhMCPrimGamma1N ;      //!
  TH1F * fhMCPrimGammaPi0N  ;   //!
  TH1F * fhMCPrimGammaPi01N  ;  //!
  TH1F * fhMCPrimGammaPi0SingleN  ;        //!
  TH1F * fhMCPrimGammaPi0Single1N  ;       //!
  TH1F * fhMCPrimGammaAllSingleN  ;        //!
  TH1F * fhMCPrimGammaAllSingle1N  ;       //!
  TH1F * fhMCPrimGammaPi0SingleNoresN  ;   //!
  TH1F * fhMCPrimGammaPi0SingleNores1N  ;  //!

  TH1F * fhMCPrimPipmN  ;                    //!  
  TH1F * fhMCPrimPipm1N  ;                   //!
  TH1F * fhMCPrimPipmNoresNa  ;              //!
  TH1F * fhMCPrimPipmNores1Na  ;             //!
  TH1F * fhMCPrimPipmNoresNb  ;              //!
  TH1F * fhMCPrimPipmNores1Nb  ;             //!
  TH1F * fhMCPrimPipmPi0  ;                  //!
  TH1F * fhMCPrimPipmPi0Nores  ;             //!
  TH1F * fhMCPrimPipmGamma  ;                //!
  TH1F * fhMCPrimPipmGammaPi0  ;             //!
  TH1F * fhMCPrimPipmGammaPi0Single  ;       //!
  TH1F * fhMCPrimPipmGammaAllSingle  ;       //!
  TH1F * fhMCPrimPipmGammaPi0SingleNores  ;  //!

  TH1F * fhPipmN  ;                           //!
  TH1F * fhPipm1N  ;                          //!
  TH1F * fhPipmTrueN  ;                       //!
  TH1F * fhPipmTrue1N  ;                      //!

  TH1F * fhGammaN[NPID] ;                  //!
  TH1F * fhGamma1N[NPID] ;                 //!
  TH1F * fhGammaTrueN[NPID] ;              //!
  TH1F * fhGammaTrue1N[NPID] ;             //!
  TH1F * fhGammaTruePi0N[NPID] ;           //!
  TH1F * fhGammaTruePi01N[NPID] ;          //!
  TH1F * fhGammaTruePi0SingleN[NPID] ;     //!
  TH1F * fhGammaTruePi0Single1N[NPID] ;    //!

  TH1F * fhGammaPipm[NPID] ;               //!
  TH1F * fhGammaPipmTrue[NPID] ;           //!
  TH1F * fhGammaPipmTruePi0[NPID] ;        //!
  TH1F * fhGammaPipmTruePi0Single[NPID] ;  //!

  //Inv masses
  TH2F * fhReal[NPID] ;                    //! 
  TH2F * fhRealTrue[NPID] ;                //!   
  TH2F * fhRealCommon[NPID] ;              //!
  TH2F * fhMixed[NPID] ;                   //!  
  TH2F * fEgammaEpi0 ;                        //!

  ClassDef(AliAnalysisPHOSFluctuations, 1); // PHOS analysis task
};
#endif
