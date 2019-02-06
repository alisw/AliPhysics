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
#include <deque>
#include "AliFemtoLambdaEventCollection2.h"
#include "AliFemtoProtonParticle.h"
#include "AliFemtoLambdaParticle.h"
#include "AliFemtoXiParticle.h"
#include "AliAnalysisUtils.h"
#include "AliFemtoCutValues.h"
#include "AliAODv0.h"
#include "AliAODpidUtil.h"
#include "AliEventCuts.h"


//class AliFemtoLambdaEventCollection2;
//class AliFemtoLambdaEvent;

class AliAODEvent;
class AliAODv0;
class AliAODMCParticle;
class AliFemtoCutValues;

class AliAnalysisTaskPLFemto : public AliAnalysisTaskSE
{

  
 public:
  
  AliAnalysisTaskPLFemto();
  AliAnalysisTaskPLFemto(const Char_t* name,Bool_t OnlineCase, TString whichV0, TString whichV0region);
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
    kProtonProtonQA,
    kProtonProtonProton,
    kAntiProtonAntiProton,
    kAntiProtonAntiProtonQA,
    kAntiProtonProton,
    kProtonV0,
    kProtonV0QA,
    kProtonProtonV0,
    kProtonProtonXi,
    kAntiProtonV0,
    kProtonAntiV0,
    kAntiProtonAntiV0,
    kAntiProtonAntiV0QA,
    kV0V0,
    kV0V0QA,
    kAntiV0AntiV0,
    kAntiV0AntiV0QA,
    kAntiV0V0,
    kProtonXi,
    kAntiProtonXi
  };

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
  void TrackCleaner();
  void ParticlePairer(Int_t zBin,Int_t multBin);
  void GetV0Origin(AliAODv0 *v0,AliAODEvent *aodEvent);
  void GetProtonOrigin(AliAODTrack *AODtrack,AliAODEvent *aodEvent,Int_t type);
  void GetMomentumMatrix(const AliFemtoLambdaParticle &v01,const AliFemtoLambdaParticle &v02);
  void GetMomentumMatrix(const AliFemtoLambdaParticle &v0,const AliFemtoProtonParticle &proton);
  void GetMomentumMatrix(const AliFemtoProtonParticle &proton1,const AliFemtoProtonParticle &proton2);
  void GetMCMomentum(AliAODTrack *aodtrack,double *Protontruth,double *ProtontruthMother,double *ProtontruthMotherParton,int *PDGcodes,AliAODEvent *aodEvent, Int_t type);
  void GetMCMomentum(AliAODv0 *v0,double *V0momtruth,double *V0momtruthMother,AliAODEvent *aodEvent);
  void GetNumberOfTPCtracksInEtaRange(AliAODEvent* aodEvent,Double_t mineta,Double_t maxeta);
  void DCAxy(const AliAODTrack *track, const AliVEvent *evt,Double_t *dcaxy);
  void BufferFiller(Int_t zBin,Int_t multBin,Int_t *NumberOfParticles);
  
  TVector3 GetGlobalPositionAtGlobalRadiiThroughTPC(const AliAODTrack *track, const Float_t bfield,double Rwanted,float PrimaryVertex[3]);
  
  void SetRun1(bool run1) { fIsRun1 = run1; }
  void SetIsLightweight(bool isLight) { fIsLightweight = isLight; }
  Bool_t SelectPID(const AliAODTrack *track,const AliAODEvent *aodevent,Int_t type);
  Bool_t SelectV0PID(const AliAODv0 *v0,TString ParOrAPar);
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *track);
  Bool_t GoodTPCFitMapSharedMap(const AliAODTrack *pTrack,const AliAODTrack *nTrack);
  Bool_t AcceptTrack(const AliAODTrack *track,int type);
  Bool_t GoodDCA(const AliAODTrack *AODtrack,const AliAODEvent *aodEvent,Int_t type);
  Bool_t SingleParticleQualityCuts(const AliAODTrack *aodtrack,const AliAODEvent *aodEvent);
  Bool_t SingleParticleQualityCuts(const AliAODv0 *v0);
  Bool_t TestPIDHypothesis(const AliAODTrack *track,Int_t type);
  Bool_t ClosePairRejecter(Float_t deltaPhiS,Float_t deltaEta);
  Bool_t CheckIfRealV0(AliAODv0 *v0,AliAODEvent *aodEvent,Int_t PDGmother,Int_t PDGdaugh1, Int_t PDGdaugh2);
  Bool_t V0TopologicalSelection(AliAODv0 *v0,AliAODEvent *aodevent,TString ParOrAPar);
  Bool_t PileUpRejection(AliAODEvent* ev);
  Bool_t CheckMCPID(AliAODEvent *aodEvent, AliAODTrack* aodtrack,int pdg);
  Bool_t CheckMCPID(AliAODv0 *v0,AliAODEvent* aodEvent,int pdg);

  Int_t GetNumberOfTrackletsInEtaRange(AliAODEvent* ev,Double_t mineta,Double_t maxeta);

  Float_t PhiS(AliAODTrack *track,const Float_t bfield,const Float_t radius,const Float_t decVtx[2]);
  Float_t PhiS(AliAODTrack *track,const Float_t bfield,const Float_t radius);
  Float_t EtaS(const Float_t zVertex, const Float_t radius);

  Double_t Qinv(TLorentzVector trackV0,TLorentzVector trackProton);
  Double_t relKcalc(TLorentzVector track1,TLorentzVector track2);
  Double_t LpCorrfunc(Double_t relk,Double_t source);
  Double_t ppCorrfunc(Double_t relk,Double_t source);
  Double_t GetAverageSeparation(const TVector3 globalPositions1st[],const TVector3 globalPositions2nd[]);

  void SetSystematics(AliFemtoCutValues::systematics cutType){fcutType = cutType;}

  AliEventCuts fAliEventCuts;
  void SetTrigger(UInt_t trigger) { fTrigger = trigger; }
  UInt_t GetTrigger() const { return fTrigger; }
  void SetV0Percentile(float v0perc) { fV0PercentileMax = v0perc; }

 private:
  
  AliAnalysisTaskPLFemto(const AliAnalysisTaskPLFemto &source);
  AliAnalysisTaskPLFemto& operator=(const AliAnalysisTaskPLFemto& source); 
  
  UInt_t fTrigger;
  bool fIsRun1;
  bool fIsLightweight;
  float fV0PercentileMax;
  Bool_t fUseMCInfo;             //  Use MC info
  Bool_t fOnlineV0;
  TList *fOutput;               //!  User output
  TList *fOutputAliEvent;       //! User output AliEventCuts
  TList *fOutputSP;             //!  User output2
  TList *fOutputTP;             //!  User output4
  TList *fOutputPID;            //!  User output3
  AliMCEvent* fAODMCEvent; //!
  TH1F *fCEvents;          //!


  TString fwhichV0,fwhichAntiV0,fwhichV0region;
  Int_t fV0Counter;
  Int_t fAntiV0Counter;
  Int_t fProtonCounter;
  Int_t fAntiProtonCounter;
  Int_t fXiCounter;

  Int_t fEventNumber;
  int fWhichfilterbit;
  //PID:
  
  Float_t *fTPCradii; //!
  Float_t fPrimVertex[3];

  AliFemtoCutValues::systematics fcutType;

  AliPIDResponse  *fPIDResponse; //! PID response object

  //for tracking mixed events etc.
  AliAODTrack     **fGTI;                  //! Array of pointers which stores global track infos
  const UShort_t  fTrackBuffSize;          //! Size fo the above arra

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

  Float_t CalcSphericityEvent(AliAODEvent *AODevent);
  Float_t DEtacalc(const double zprime1, const double zprime2, const double radius);
  Float_t DPhicalc(const double dxprime, const double dyprime, const double radius);

  enum 
  {
    kZVertexBins = 10,//10
    kEventsToMix =  10,//15
    kV0TrackLimit   = 100, //maximum number of v0s, array size
    kProtonTrackLimit = 100, //maximum number of protons in array
    kXiTrackLimit = 100, //maximum number of Xis in array
    kMultiplicityBins = 13,//5
  };


  //Check buffer with deque (mainly to avoid empty evt mixing):

  enum
  {
    kProton,
    kAntiProton,
    kV0,
    kAntiV0,
    kXi
  };

  Bool_t fSEPairAnalysisDecider[5];

  //every particle has its own vector
  std::vector<AliFemtoProtonParticle> fProtonTrackVector;
  std::vector<AliFemtoProtonParticle> fAntiProtonTrackVector;
  std::vector<AliFemtoLambdaParticle> fLambdaTrackVector;
  std::vector<AliFemtoLambdaParticle> fAntiLambdaTrackVector;
  std::vector<AliFemtoXiParticle>     fXiTrackVector;

  //every particle has also its own buffer
  std::deque< std::vector< AliFemtoProtonParticle > > fProtonEvtBuffer[kZVertexBins][kMultiplicityBins];
  std::deque< std::vector< AliFemtoProtonParticle > > fAntiProtonEvtBuffer[kZVertexBins][kMultiplicityBins];
  std::deque< std::vector< AliFemtoLambdaParticle > > fLambdaEvtBuffer[kZVertexBins][kMultiplicityBins];
  std::deque< std::vector< AliFemtoLambdaParticle > > fAntiLambdaEvtBuffer[kZVertexBins][kMultiplicityBins];
  std::deque< std::vector< AliFemtoXiParticle > >     fXiEvtBuffer[kZVertexBins][kMultiplicityBins];

  ClassDef(AliAnalysisTaskPLFemto,3)
};

#endif

