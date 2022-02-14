/*Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */
//NUOVO
// Author:Ramona Lea based on M. Nicassio m.nicassio@gsi.de maria.nicassio@cern.ch

#ifndef ALIANALYSISTASKK0SPFEMTO_H
#define ALIANALYSISTASKK0SPFEMTO_H
class AliPIDResponse;
class AliMultSelection;
class AliAODEvent;
class AliAODtrackCuts;
//class THnSparse;
//class AliAnalysisK0SPEventCollection;
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisK0SPEventCollection.h"
#include "AliEventCuts.h"

class AliAnalysisTaskK0SPFemto : public AliAnalysisTaskSE  
{
 public:
  AliAnalysisTaskK0SPFemto();
  AliAnalysisTaskK0SPFemto(const char *name);
  virtual                 ~AliAnalysisTaskK0SPFemto();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);

  // void                    ProcessMCParticles();
  void                    SetCentrality   (Float_t lowlimcent = 0., Float_t uplimcent = 90.) { fCentrLowLim = lowlimcent;  fCentrUpLim = uplimcent; }
  void                    SetFilterBit(Int_t filterbit = 4){fFilterBit=filterbit;};
  void                    SetMC(Bool_t useMC){fIsMC=useMC;};
  Bool_t                  GetMC(){return fIsMC;};
  void                    SetCollidingSystem (const char* collSystem = "pp") { fCollidingSystem = collSystem ;};

  void DoPairsh1h2    (const Float_t lcentrality, Int_t fieldsign, const Double_t fSphericityvalue, const Double_t fSpherocityvalue);
  double CalculateKstar(double momentum1[3], double momentum2[3], double mass1, double mass2); 
  double CalculateSphericityofEvent(AliAODEvent *aodEvent);
  double CalculateSpherocityEvent(AliAODEvent *aodEvent);

  AliEventCuts fEventCuts;

 private:
  
  TString     fCollidingSystem;
  Int_t Neventi;
  AliAODEvent*  fAOD;             //! input event
  Bool_t        fIsMC;

  Double_t fNominalMassK0s=0.4976;
  //  AliAODtrackCuts    *fAODtrackCuts;              //! basic cut variables for tracks added ! not sure
  AliPIDResponse* fPIDResponse;	     //!

  Float_t fCentrLowLim;                           // centrality settings: lower limit
  Float_t fCentrUpLim;                            // upper limit


  TList*        fOutputContainer;    //! output list
  TTree*        fHistSparseSignal;   //! output tree
  TTree*        fHistSparseBkg;      //! output tree

  Int_t     *farrGT; //!
  UShort_t  fTrackBufferSize;          // Size fo the above array, ~12000 for PbPb

  AliAnalysisK0SPEventCollection ***fEventColl; //
  AliAnalysisK0SPEvent *fEvt;                   // 
    
  Int_t     fMaxFirstMult;
  Int_t     fMaxSecondMult;
  Int_t fzVertexBins;
  Int_t fnCentBins;
  //  Int_t fnSphericityBins;
  Short_t fnEventsToMix;
  Int_t       fFilterBit;
  Bool_t      fHMtrigger;
  Float_t fPDGMassFirst;
  Int_t fPDGcodeFirst;
  Float_t fPDGMassSecond;
  Int_t fPDGcodeSecond;

  Float_t tCentrality;
  Float_t tSphericity;
  Float_t tSpherocity;
  Float_t tKtpair;
  Float_t tkStar;
  Bool_t  tIsCommonParton;

  Float_t tPtV0;
  Float_t tPTotV0;
  Float_t tThetaV0;
  Float_t tPhiV0;
  Float_t tDcaPosV0;
  Float_t tDcaNegV0;
  Float_t tInvMassK0s;
  Float_t tInvMassLambda;
  Float_t tInvMassAntiLambda;
  Float_t tCosPointingAngleV0;


  Float_t tPtP;
  Float_t tPTotP;
  Float_t tThetaP;
  Float_t tPhiP;
  Int_t   tSignP;
  Float_t tDCAxyP;
  Float_t tDCAzP;
  Float_t tMassTOFP;


  Int_t tMCtruepair;
  Int_t tMCSameMother;
  Int_t tMCMotherV0;
  Int_t tMCMotherP;
  Int_t tMCptcTypeV0;
  Int_t tMCptcTypeP;
  Int_t tMCSameGM;
  Int_t tMotherPDG;
  Int_t tpdgcodeV0;
  Int_t tpdgcodeP;
  Float_t tKstarGen;

  TH1F  *fHistEventMultiplicity;                //! Counts Event multiplicity
  //  TH1I  *hmult;                                 //! Counts tracklets
  TH1F  *fHistCentrality;                       //! Counts events in centrality
  TH1F  *fHistVertexDistribution;               //! Vertex distribution
  TH1F  *fHistSphericity; 
  TH1F  *fHistSpherocity; 

  TH1F  *fHistMassK0S;
  TH2F  *fHistFirstNPionTPCdEdx;
  TH2F  *fHistFirstPPionTPCdEdx;
  TH2F  *fHistSecondTPCdEdx;
  TH2F  *fHistSecondMassTOFvsPt3sTPC;

  AliAnalysisTaskK0SPFemto(const AliAnalysisTaskK0SPFemto &);
  AliAnalysisTaskK0SPFemto &operator=(const AliAnalysisTaskK0SPFemto &);
  
  ClassDef(AliAnalysisTaskK0SPFemto, 3);
};

#endif
