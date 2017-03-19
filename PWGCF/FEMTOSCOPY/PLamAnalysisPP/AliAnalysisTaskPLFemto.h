#ifndef ALIANALYSISTASKPLFEMTO_H
#define ALIANALYSISTASKPLFEMTO_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */ 

//#include <TH2F.h>
#include "AliAnalysisTaskSE.h"
#include <vector>
#include "AliFemtoLambdaEventCollection2.h"
#include "AliFemtoProtonParticle.h"
#include "AliFemtoLambdaParticle.h"
#include "AliFemtoXiParticle.h"
#include "AliAnalysisUtils.h"
#include "AliFemtoCutValues.h"
#include "AliAODv0.h"
#include "AliAODpidUtil.h"
//#include "TRandom3.h"

//class AliFemtoLambdaEventCollection2;
//class AliFemtoLambdaEvent;

class AliAODEvent;
class AliAODv0;
class AliAODMCParticle;

class AliAnalysisTaskPLFemto : public AliAnalysisTaskSE
{

  
 public:
  

  AliAnalysisTaskPLFemto();
  //AliAnalysisTaskPLFemto(const Char_t* name,AliRDHFCutsDStartoKpipi* cuts);
  AliAnalysisTaskPLFemto(const Char_t* name,Bool_t OnlineCase, TString whichV0, TString whichV0region,const int whichfilterbit);
  virtual ~AliAnalysisTaskPLFemto();

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  

  enum pairType
  {
    kProtonProton,
    kProtonProtonProton,
    kAntiProtonAntiProton,
    kAntiProtonProton,
    kProtonV0,
    kProtonProtonV0,
    kProtonProtonXi,
    kAntiProtonV0,
    kProtonAntiV0,
    kAntiProtonAntiV0,
    kV0V0,
    kAntiV0AntiV0,
    kAntiV0V0,
    kProtonXi,
    kAntiProtonXi
  };


  
  Bool_t SelectPID(const AliAODTrack *track,const AliAODEvent *aodevent, Int_t type);
  Bool_t SelectV0PID(const AliAODv0 *v0,TString ParOrAPar);
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *track);
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,const AliAODTrack *nTrack);
  Bool_t AcceptTrack(const AliAODTrack *track,int type);
  //Bool_t GetMC() const {return fUseMCInfo;}
  Bool_t GoodDCA(const AliAODTrack *track,const AliAODEvent *aodEvent,Int_t type);
  Bool_t SingleParticleQualityCuts(const AliAODTrack *aodtrack,const AliAODEvent *aodEvent);
  Bool_t SingleParticleQualityCuts(const AliAODv0 *v0);
  Bool_t TestPIDHypothesis(const AliAODTrack *track,Int_t type);
  Bool_t ClosePairRejecter(Float_t deltaPhiS,Float_t deltaEta);
  Bool_t CheckIfRealV0(AliAODv0 *v0,AliAODEvent *aodEvent,Int_t PDGmother,Int_t PDGdaugh1, Int_t PDGdaugh2);
  Bool_t V0TopologicalSelection(AliAODv0 *v0,AliAODEvent *aodevent,TString ParOrAPar);
  Bool_t PileUpRejection(AliAODEvent* ev);
  Bool_t CheckMCPID(AliAODEvent *aodEvent, AliAODTrack* aodtrack,int pdg);
  Bool_t CheckMCPID(AliAODv0 *v0,AliAODEvent* aodEvent,int pdg);

  void DefineHistograms(TString whichV0);

  void FillV0Selection(AliAODv0 *v0,AliAODEvent* aodEvent,TString ParOrAPar);
  void SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliAODTrack *track);
  void PairAnalysis(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton,Int_t REorME,pairType inputPair);
  void PairAnalysis(const AliFemtoLambdaParticle &v01,const AliFemtoLambdaParticle &v02,Int_t REorME,pairType inputPair);
  void PairAnalysis(const AliFemtoXiParticle &xi,const AliFemtoProtonParticle &proton,Int_t REorME,pairType inputPair);
  void PairAnalysis(const AliFemtoProtonParticle &protonA,const AliFemtoProtonParticle &protonB,Int_t REorME,pairType inputPair);
  void TripleAnalysis(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2,Int_t REorME,pairType inputPair);
  void TripleAnalysis(const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2,const AliFemtoProtonParticle &proton3,Int_t REorME,pairType inputPair);
  void ProtonSelector(AliAODEvent *aodEvent);
  void V0Selector(AliAODEvent *aodEvent);
  void XiSelector(AliAODEvent *aodEvent);
  void FIFOShifter(Int_t zBin,Int_t MultBin);
  void TrackCleaner();
  void ParticlePairer();
  void GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield,double globalPositionsAtRadii[9][3], double PrimaryVertex[3]);
  void GetV0Origin(AliAODv0 *v0,AliAODEvent *aodEvent);
  void GetProtonOrigin(AliAODTrack *proton,AliAODEvent *aodEvent,Int_t type);
  void GetMomentumMatrix(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton);
  void GetMomentumMatrix(const AliFemtoProtonParticle &protonA,const AliFemtoProtonParticle &protonB);
  void GetMCMomentum(AliAODTrack *aodTrack,double *Protontruth,double *ProtontruthMother,double *ProtontruthMotherParton,int *PDGcodes,AliAODEvent *aodEvent,Int_t type);
  void GetMCMomentum(AliAODv0 *v0,double *V0momtruth,double *V0momtruthMother,AliAODEvent *aodEvent);
  void GetAverageSeparation(const double globalPositions1st[9][3],const double globalPositions2nd[9][3],Int_t REorME,TString whichpair,double *averageSeparation);
  void CalculateDecayMatrix();
  
  Int_t *BufferFiller();
  Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev,Double_t mineta,Double_t maxeta);
  
  Float_t PhiS(AliAODTrack *track,const Float_t bfield,const Float_t radius,const Float_t decVtx[2]);
  Float_t PhiS(AliAODTrack *track,const Float_t bfield,const Float_t radius);
  Float_t EtaS(const Float_t zVertex, const Float_t radius);
  Double_t *DCAxy(const AliAODTrack *track, const AliVEvent *evt);
  Double_t Qinv(TLorentzVector v0,TLorentzVector proton);
  Double_t relKcalc(TLorentzVector track1,TLorentzVector track2);
  Double_t LpCorrfunc(Double_t relk,Double_t source);
  Double_t ppCorrfunc(Double_t relk,Double_t source);

    
 private:
  
  AliAnalysisTaskPLFemto(const AliAnalysisTaskPLFemto &source);
  AliAnalysisTaskPLFemto& operator=(const AliAnalysisTaskPLFemto& source); 
  
  Int_t  fEvents;                //  n. of events
  Bool_t fUseMCInfo;             //  Use MC info
  Bool_t fOnlineV0;
  TList *fOutput;               //!  User output
  TList *fOutputSP;             //!  User output2
  TList *fOutputTP;             //!  User output4
  TList *fOutputPID;            //!  User output3
  AliMCEvent* fAODMCEvent; //!
  TH1F *fCEvents;          //!

  Double_t *fTPCradii; //!
  Double_t fSphericityvalue;
  Double_t fPrimVertex[3];
    
  TString fwhichV0,fwhichAntiV0,fwhichV0region;
  Int_t fV0Counter;
  Int_t fAntiV0Counter;
  Int_t fProtonCounter;
  Int_t fAntiProtonCounter;
  Int_t fXiCounter;
  
  Int_t fEventNumber;
  int fWhichfilterbit;
  //PID:
  AliPIDResponse  *fPIDResponse; //! PID response object
  //AliTPCPIDResponse *fTPCResponse;
  
  //for tracking mixed events etc.
  AliAODTrack     **fGTI;                  //! Array of pointers which stores global track infos
  const UShort_t  fTrackBuffSize;          //! Size fo the above arra

  AliFemtoLambdaEventCollection2 ***fEC; //!
  AliFemtoLambdaEvent *fEvt; //!
  AliFemtoLambdaParticle *fV0cand; //!
  AliFemtoLambdaParticle *fAntiV0cand; //!
  AliFemtoProtonParticle *fProtoncand; //!
  AliFemtoProtonParticle *fAntiProtoncand; //!
  AliFemtoXiParticle *fXicand; //!
  AliFemtoCutValues *fCuts; //!
  AliAnalysisUtils *fAnaUtils; //! Object to use analysis utils like pile-up rejection


  //Functions:
  Bool_t CheckGlobalTrack(const AliAODTrack *track);
  Int_t *MultiplicityBinFinderzVertexBinFinder(AliAODEvent *aodEvent,Double_t zVertex);
  Float_t GetCorrectedTOFSignal(const AliVTrack *track,const AliAODEvent *aodevent); //!
  
  Double_t CalcSphericityEvent(AliAODEvent *AODevent);
  Double_t DEtacalc(const double zprime1, const double zprime2, const double radius);
  Double_t DPhicalc(const double dxprime, const double dyprime, const double radius);

  enum 
  {
    kZVertexBins = 10,//10
    kEventsToMix =  10,//15
    kV0TrackLimit   = 100,              //maximum number of v0s, array size
    kProtonTrackLimit = 100, //maximum number of protons in array
    kXiTrackLimit = 100, //maximum number of Xis in array
    kZVertexStart = -10,
    kZVertexStop = 10,
    kMultiplicityBins = 13,//5
  };
  
#ifdef __ROOT__
  ClassDef(AliAnalysisTaskPLFemto,1);
#endif
};

#endif

