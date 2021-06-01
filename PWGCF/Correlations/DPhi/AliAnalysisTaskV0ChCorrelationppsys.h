/* Copyright(c) 1998-2013, ALICE Experiment at CERN, All rights reserved.  
 * See cxx source for full Copyright notice  
 *
 * AliAnalysisTaskV0ChCorrelationppsys class
 *
 * The task selects candidates for K0s, Lambdas and AntiLambdas(trigger particles)
 * and calculates correlations with charged unidentified particles (associated particles) in phi and eta. 
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 * Authored by Zhongabo Yin,  Zhongbao.Yin@cern.ch
 *The task has been updated by mustafa.naji.anaam@crn.ch//last updated for sys 2020
 */

//#ifndef ALIANALYSISTASKXICHCORRELATIONS_H
#ifndef AliAnalysisTaskV0ChCorrelationppsys_H
#define AliAnalysisTaskV0ChCorrelationppsys_H

//#define ALIANALYSISTASKXICHCORRELATIONS_H
#include "AliVParticle.h"
class TH1F;
class TH1D;
class TH2F;
class THnSparse;
class TList;
class AliPIDResponse;
class AliEventPoolManager;
class TH1I;
class AliAODv0;
class AliVParticle;
class AliEventCuts;

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"    
#endif

class AliAnalysisTaskV0ChCorrelationppsys : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskV0ChCorrelationppsys();//
   AliEventCuts *fEventCuts; /// Event cuts
   AliAnalysisTaskV0ChCorrelationppsys(const char *name,
                                        Double_t centMin, Double_t centMax,
                                        Bool_t effCorr);
   virtual ~AliAnalysisTaskV0ChCorrelationppsys();

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     Terminate(Option_t *);
   
   void SetAnalysisMC(Bool_t AnalysisMC = kFALSE) {fAnalysisMC = AnalysisMC;}
   void SetTrackPileUpCut(Bool_t RejectTrackPileUp = kTRUE) {fRejectTrackPileUp = RejectTrackPileUp;}
  void SetV0PileUpCut(Bool_t RejectV0PileUp = kTRUE) {fRejectV0PileUp = RejectV0PileUp;}
   void SetEventPileUpCut(Bool_t cut) { fRejectEventPileUp = cut; }

   void SetOStatus(Int_t OStatus = 1) {fOStatus = OStatus;}
   //void SetEfficiencyCorrection(Bool_t effCorr = kTRUE){fEffCorr = effCorr;}

  //----------------------Mixing part----------------------------------
   void SetMixingTracks(Double_t value){fMixingTracks = value;}
   void SetPoolSize(Double_t value){fPoolSize = value;}
  //-----------------------Variable----------------------------
   void SetVtxCut(Double_t value){fPrimaryVertexCut = value;}  
   void SetNumOfVzBins(Double_t value){fNumOfVzBins = value;}  

   void SetVtxXMin(Double_t value){fVtxXMin = value;}
   void SetVtxYMin(Double_t value){fVtxYMin = value;}
   void SetVtxZMin(Double_t value){fVtxZMin = value;}
   void SetCentMin(Double_t value){fCentMin = value;}
   void SetCentMax(Double_t value){fCentMax = value;}
   //----------------------Track---------------------------
   void SetTrackMCPtMin(Double_t value){fTrackMCPtMin = value;} 
   void SetTrackPtMin(Double_t value){fTrackPtMin = value;} 
   void SetTrackPtMax(Double_t value){fTrackPtMax = value;}
   void SetTrackEta(Double_t value){fTrackEta = value;}
   void SetFilterBit(Int_t value){fFilterBit = value;}
   void SetAssocNcls(Double_t value){fAssocNcls = value;}
   //---------------------V0---------------------------------
   void SetV0MCPtMin(Double_t value){fV0MCPtMin = value;}
   void SetV0PtMin(Double_t value){fV0PtMin = value;}
   void SetV0PtMax(Double_t value){fV0PtMax = value;}
   void SetV0Eta(Double_t value){fV0Eta = value;}
   void SetK0sLifeTimeMin(Double_t value){fK0sLifeTimeMin = value;}
   void SetK0sLifeTimeMax(Double_t value){fK0sLifeTimeMax = value;}
   void SetLambdaLifeTimeMin(Double_t value){fLambdaLifeTimeMin = value;}
   void SetLambdaLifeTimeMax(Double_t value){fLambdaLifeTimeMax = value;}
   
   void SetV0DaughterPtMinCut(Double_t value){fV0DaughterPtMinCut = value;}
   void SetDCANegtoPrimVertex(Double_t value){fDCANegtoPrimVertex = value;}
   void SetDCAPostoPrimVertex(Double_t value){fDCAPostoPrimVertex = value;}
   void SetDCAV0DaughtersMax(Double_t value){fDCAV0DaughtersMax = value;}

  void Setk0sCPA(Double_t value){fk0sCPA = value;}
  void SetLambdaCPA(Double_t value){fLambdaCPA = value;}
   void SetCosPointingAngleMin(Double_t value){fCosPointingAngleMin = value;}
   void Set2DFiducialMin(Double_t value){f2DFiducialMin = value;}
   
   void SetV0DaughterTrackTPCCluster(Double_t value){fV0DaughterTrackTPCCluster = value;}
   void SetNCrossedRowsTPCfindable(Double_t value){fNCrossedRowsTPCfindable = value;}
   
   void SetK0sMassWindow(Double_t value){fK0sMassWindow = value;}
   void SetLambdaMassWindow(Double_t value){fLambdaMassWindow = value;}
   void SetPtArmV0AlphaV0(Double_t value){fPtArmV0AlphaV0 = value;}
   void SetLambdaCosPointingAngleMin(Double_t value){fLambdaCosPointingAngleMin = value;}
   void SetAntiLambdaCosPointingAngleMin(Double_t value){fAntiLambdaCosPointingAngleMin = value;}
   void SetLambdaAlphaV0Min(Double_t value){fLambdaAlphaV0Min = value;}
   void SetAntiLambdaAlphaV0Max(Double_t value){fAntiLambdaAlphaV0Max = value;}
   void SetLambdaDCA2PVMax(Double_t value){fLambdaDCA2PVMax = value;}
   void SetAntiLambdaDCA2PVMax(Double_t value){fAntiLambdaDCA2PVMax = value;}
  // void SetCorrelationsAnalysis(Bool_t corr) { fCorrelations = corr; }//////////////////////////////

   //-----------------------------PID---------------------------------
   void SetV0PIDSigma(Double_t value){fV0PIDSigma = value;}
   //-----------------------------------------------------------------

   Bool_t IsGoodPrimaryTrack(const AliAODTrack* t);
	//Bool_t IsK0InvMass(const Double_t mass)const;////////////////////
   Bool_t IsGoodDaughterTrack(const AliAODTrack* t);

   Int_t GetOStatus() { return fOStatus; }

   Bool_t IsGoodV0(AliAODv0* aodV0 ,Int_t oSta);

   Bool_t IsMCV0Primary(AliAODv0 *v0, Int_t specie);//0 for K0S, 1 for Lambda, 2 for AntiLambda
   Bool_t IsMCV0FromXi(AliAODv0 *v0, Int_t specie);
   Int_t IsMcV0(AliAODv0 *v0, Int_t specie)const;
   Int_t GetV0Label(Int_t lab1, Int_t lab2, Int_t specie)const;


  
 // void RemovingInjectedSignal(TObjArray* tracks, TObject* mcObj, Int_t maxLabel);

private:

   AliAnalysisTaskV0ChCorrelationppsys(const AliAnalysisTaskV0ChCorrelationppsys&);            //not implemented
   AliAnalysisTaskV0ChCorrelationppsys& operator=(const AliAnalysisTaskV0ChCorrelationppsys&); //not implemented 
 
   void AddQATrackCandidates();
   void AddQAV0Candidates();
   void AddQAEvent();
   void AddQAAnalysisK0s();
   void AddQAAnalysisLambda();
   void AddQAAnalysisAntiLambda();


   static const int    kNVtxZ = 8;
  
   Double_t        fMixingTracks;  // size of track buffer for event mixing
   Double_t        fPoolSize;
   AliEventPoolManager     *fPoolMgr;  //! event pool manager
   TList           *fOutput;  //! Output list
   TList           *fOutput2; //! Output list
   TList           *fOutput3; //! Output list
   TList           *fOutput4; //! Output list
   TList           *fOutput5; //! Output list
   TList           *fOutput6; //! Output list
   TList           *fOutput7; //! Output list
   TList           *fOutput8; //! Output list
   


   //TF1             *fMultiplicityV0McorrCut;
   AliPIDResponse  *fPIDResponse;  // PID response
   Double_t        fPrimaryVertexCut;//fNumOfVzBins
   Double_t        fNumOfVzBins;//
   Double_t        fVtxXMin;
   Double_t        fVtxYMin;
   Double_t        fVtxZMin;
   Double_t        fCentMin;
   Double_t        fCentMax;
   //---------------------------------Track-------------------------------------
   Double_t        fTrackMCPtMin;
   Double_t        fTrackPtMin;
   Double_t        fTrackPtMax;
   Double_t        fTrackEta;
   Int_t           fFilterBit;
   Double_t        fAssocNcls;
   //---------------------------------V0----------------------------------------
   Double_t        fV0MCPtMin;
   Double_t        fV0PtMin;
   Double_t        fV0PtMax;
   Double_t        fV0Eta;
   Double_t        fK0sLifeTimeMin;
   Double_t        fK0sLifeTimeMax;
   Double_t        fLambdaLifeTimeMin;
   Double_t        fLambdaLifeTimeMax;
   
   Double_t        fV0DaughterPtMinCut;
   Double_t        fDCANegtoPrimVertex;
   Double_t        fDCAPostoPrimVertex;
   Double_t        fDCAV0DaughtersMax;
   Double_t        fCosPointingAngleMin;

   Double_t        fCPA;
   Double_t        fk0sCPA;                
   Double_t        fLambdaCPA;
   Double_t        f2DFiducialMin;
   Double_t        fV0DaughterTrackTPCCluster;
   Double_t        fNCrossedRowsTPCfindable;
   
   Double_t        fK0sMassWindow;
   Double_t        fLambdaMassWindow;
   Double_t        fPtArmV0AlphaV0;
   Double_t        fLambdaCosPointingAngleMin;
   Double_t        fAntiLambdaCosPointingAngleMin;
   Double_t        fLambdaAlphaV0Min;
   Double_t        fAntiLambdaAlphaV0Max;
   Double_t        fLambdaDCA2PVMax;
   Double_t        fAntiLambdaDCA2PVMax;
   Double_t        massK0s;

   TObjArray       *selectedK0s;//!
   TObjArray       *selectedLambda;//!
   TObjArray       *selectedAntiLambda;//!
   TObjArray       *selectedTracks;//!
   TObjArray       *trigParticles;//!

   //-----------------------------PID-------------------------------------------
   Double_t        fV0PIDSigma;
   //-----------------------------MC--------------------------------------------
   TH2F             *fHistRCK0sPt;//!
   TH2F             *fHistRCLambdaPt;//!
   TH2F             *fHistRCPrimLambdaDCAPt;//!
   TH2F             *fHistRCSecLambdaDCAPt;//!
   TH2F             *fHistRCAntiLambdaPt;//!
   TH2F             *fHistRCPrimAntiLambdaDCAPt;//!
   TH2F             *fHistRCSecAntiLambdaDCAPt;//!
   TH2F             *fHistRCTrkPtEta;//!
   TH2F             *fHistRCPriTrkPtEta;//!
   TH2F             *fHistRCSecTrkPtEta;//!
   TH2F             *fHistMCtruthTrkPtEta;//!
   TH2F             *fHistMCtruthK0sPt;//!
   TH2F             *fHistMCtruthLambdaPt;//!
   TH2F             *fHistMCtruthAntiLambdaPt;//!
                       
   Int_t           fOStatus;   // checks for online and offline V0s 

   Bool_t            fAnalysisMC;
   Bool_t            fRejectTrackPileUp;// 
   Bool_t            fRejectV0PileUp;// 
   Bool_t            fRejectEventPileUp;
   TClonesArray     *fMCArray;//! MC array for AOD
   //----------------------------Correction-------------------------------
   Bool_t            fEffCorr;
   TList            *fEffList;  //!  
   TH2F             *fHistEffEtaPtK0s;//!
   TH2F             *fHistEffEtaPtLambda;//!
   TH2F             *fHistEffEtaPtAntiLambda; //!
   TH2F             *fHistEffEtaPtTrack;//!

      //----------------------------------------------------------------------
   Double_t         fMassMean[2]; //mass mean for K0s, Lambda 
   Double_t         fMassRes[2];  //mass resolution for K0s, Lambda

   Double_t         fMassLowK0s[3];
   Double_t         fMassHighK0s[3];
   Double_t         fMassLowLambda[3];
   Double_t         fMassHighLambda[3];
   Double_t         fMassLowAntiLambda[3];
   Double_t         fMassHighAntiLambda[3];

   TH3F             *fHistRCTrk[kNVtxZ];//!
   TH3F             *fHistRCPriTrk[kNVtxZ];//!
   TH3F             *fHistRCSecTrk[kNVtxZ];//!
   TH3F             *fHistMCtruthTrk[kNVtxZ];//!
   TH3F             *fHistMCtruthK0s[kNVtxZ];//!
   TH3F             *fHistMCtruthLambda[kNVtxZ];//!
   TH3F             *fHistMCtruthAntiLambda[kNVtxZ];//!
   TH3F             *fHistRCK0s[kNVtxZ];//!
   TH3F             *fHistRCLambda[kNVtxZ];//!
   TH3F             *fHistRCAntiLambda[kNVtxZ];//!



   Double_t         fBestPrimaryVtxPos[3];
   ClassDef(AliAnalysisTaskV0ChCorrelationppsys, 1); // class for CascadeCh correlation analysis
};


class AliV0hParticles : public AliVParticle
{
    public:




    
    AliV0hParticles(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t iDh)//, Bool_t status)
    : fEta(eta), fPhi(phi), fPt(pt), fCandidate(candidate), fLabel(label), fIDh(iDh)//, fRecStatus(status)
    {
    }
    AliV0hParticles(Float_t eta, Float_t phi, Float_t pt, Short_t candidate, Int_t label,Int_t idpos, Int_t idneg, /*Bool_t status,*/Double_t mass)
    : fEta(eta), fPhi(phi), fPt(pt), fCandidate(candidate), fLabel(label), fIDpos(idpos), fIDneg(idneg),/* fRecStatus(status), */fMass(mass)
    {
    }



AliV0hParticles(Float_t eta, Float_t phi, Float_t pt,Double_t mass)//, Short_t candidate, Int_t label,Int_t iDh)
    : fEta(eta), fPhi(phi), fPt(pt),fMass(mass)//, fCandidate(candidate), fLabel(label), fIDh(iDh)
    {
    }



  AliV0hParticles(Double_t pt,Double_t phi,Double_t eta,  Double_t mass,Short_t candidate)
    :fEta(eta),fPhi(phi),fPt(pt), fMass(mass),fCandidate(candidate)
   {
    }



   virtual ~AliV0hParticles() {}

    // kinematics                                             
    virtual Double_t Px() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Py() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Pz() const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Pt() const { return fPt; }
    virtual Double_t P()  const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   PxPyPz(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Xv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Yv() const { AliFatal("Not implemented"); return 0; }
    virtual Double_t Zv() const { AliFatal("Not implemented"); return 0; }
    virtual Bool_t   XvYvZv(Double_t[3]) const { AliFatal("Not implemented"); return 0; }

    virtual Double_t OneOverPt()  const { AliFatal("Not implemented"); return 0; }

    virtual Double_t Phi()        const { return fPhi; }
    virtual Double_t Theta()      const { AliFatal("Not implemented"); return 0; }
    virtual Double_t E()          const { AliFatal("Not implemented"); return 0; }
    virtual Double_t M()          const { return fMass; }

    virtual Double_t Eta()        const { return fEta; }
    virtual Double_t Y()          const { AliFatal("Not implemented"); return 0; }

    virtual Short_t Charge()      const { AliFatal("Not implemented"); return 0; }
    virtual Int_t   GetLabel()    const { AliFatal("Not implemented"); return 0; }
    // PID                       
    virtual Int_t   PdgCode()     const { AliFatal("Not implemented"); return 0; }
    virtual const Double_t *PID() const { AliFatal("Not implemented"); return 0; }
        


    //    virtual Double_t Dca2PV() const {return fDca2PV;}
    virtual Short_t WhichCandidate() const { return fCandidate; }

   Int_t   MyLabel()             const { return fLabel; }
    Int_t   GetIDCh()             const { return fIDh; }
    Int_t   GetIDPos()            const { return fIDpos; }
    Int_t   GetIDNeg()            const { return fIDneg; }
   // Bool_t  GetRecStatus()        const { return fRecStatus; }
    Double_t GetMass()            const { return fMass; }
    

      private:
      Double_t fEta;      // eta
      Double_t fPhi;      // phi
      Double_t fPt;       // pT
      Short_t fCandidate;   // V0 candidate: 1 - K0, 2 - Lambda, 3 - Antilambda, 4 - ChTrack
      Int_t fLabel;   // Label of MC particles
      Int_t fIDh;   // Label
      Int_t fIDpos;   // Label of possitive charged daughter
      Int_t fIDneg;   // Label of negative charged daughter
      Bool_t fRecStatus;   // reconstruction status
      Double_t fMass; // mass
  
      ClassDef( AliV0hParticles, 1) // class required for correlatios calculation and event mixing
 };

#endif
