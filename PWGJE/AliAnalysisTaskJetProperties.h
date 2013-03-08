// *************************************************************************************
// * Task for Jetproperties and jet shape analysis in PWG4 Jet Task Force Train for pp *
// *************************************************************************************

#ifndef ALIANALYSISTASKJETPROPERTIES_H
#define ALIANALYSISTASKJETPROPERTIES_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class AliESDEvent;
class AliAODEvent;
class AliAODExtension;
class TList;
class TH1F;
class TH2F;
class TH3F;
class TProfile;
class THnSparse; 
class TRandom3;
class TArrayS;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskJetProperties : public AliAnalysisTaskSE {
  
 public:
  AliAnalysisTaskJetProperties(); 
  AliAnalysisTaskJetProperties(const char *name);
  virtual ~AliAnalysisTaskJetProperties();
  
  virtual Bool_t Notify();
  virtual void   UserCreateOutputObjects();
  virtual void   Init();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t* );
  
  virtual void   SetJetBranch(const char* c){fBranchJets = c;}
  virtual void   SetNonStdFile(char* c){fNonStdFile = c;}
  virtual void   SetTrackType(Int_t i){fTrackType = i;}
  virtual void   SetEventCuts(Float_t VtxZ=10.,Int_t nContributors=2)
  {fMaxVertexZ = VtxZ; fNContributors = nContributors;}
  virtual void   SetTrackCuts(Float_t trackPt = 0.15, Float_t trackEtaMin = -0.9, Float_t trackEtaMax = 0.9)
  {fTrackPtCut = trackPt; fTrackEtaMin = trackEtaMin; fTrackEtaMax = trackEtaMax;}
  virtual void   SetJetCuts(Float_t jetPt = 5., Float_t jetEtaMin = -0.5, Float_t jetEtaMax = 0.5) 
  {fJetPtCut = jetPt; fJetEtaMin = jetEtaMin; fJetEtaMax = jetEtaMax;}
  virtual void   SetJetRejectType(Int_t i){fJetRejectType = i;}
  virtual void   SetFilterMask(UInt_t i) {fFilterMask = i;}
  virtual void   UsePhysicsSelection(Bool_t b) {fUsePhysicsSelection = b;}
  virtual void   SetJetRadius(Float_t r){fJetRadius = r;}
  
  enum {kTrackUndef=0, kTrackAOD, kTrackKine,kTrackAODMC};//for track selection
  enum {kNoReject=0, kReject1Track};//for jet rejection
  
 protected:
  Int_t	   GetListOfJetTracks(TList* l, const AliAODJet* j);
  Int_t	   GetListOfJets(TList* list);
  void     FillJetProperties(TList *jetlist);
  void     FillJetShape(TList *jetlist);
  void     FillJetShapeUE(TList *jetlist);
  void     FillFFCorr(TList *jetlist);

  void     GetTracksTiltedwrpJetAxis(Float_t alpha, TList* inputlist, TList* outputlist, const AliAODJet* jet, Double_t radius,Double_t& sumPt);
  Int_t    GetListOfTracks(TList *list, Int_t type);
  
  AliESDEvent*     fESD;          // ESD event
  AliAODEvent*     fAOD;          // AOD event
  AliAODEvent*     fAODJets;      // AOD event with jet branch (case we have AOD both in input and output)
  AliAODExtension* fAODExtension; //! where we take the jets from can be input or output AOD
  //AliMCEvent*  fMCEvent;  // MC event
  
  TString fNonStdFile;          // name of delta aod file to catch the extension
  TString fBranchJets;          // branch name for reconstructed jets
  Int_t   fTrackType;           // type of generated tracks
  Int_t   fJetRejectType;       // type of jets rejected
  Bool_t  fUseAODInputJets;     // take jets from in/output - only relevant if AOD event both in input AND output and we want to use output
  UInt_t  fFilterMask;	        // filter bit for selected tracks
  Bool_t  fUsePhysicsSelection; // switch for event selection
  Float_t fMaxVertexZ;          // maximum abs(z) position of primiary vertex [cm]
  Int_t   fNContributors;       // contributors to primary vertex
  // track cuts
  Float_t fTrackPtCut;          // track transverse momentum cut
  Float_t fTrackEtaMin;         // track eta cut
  Float_t fTrackEtaMax;         // track eta cut
  // jet cuts
  Float_t fJetPtCut;            // jet transverse momentum cut
  Float_t fJetEtaMin;           // jet eta cut
  Float_t fJetEtaMax;           // jet eta cut
  Float_t fAvgTrials;           // average number of trials per event
  Float_t fJetRadius;
  
  TList	    *fJetList;                //! List of jets
  TList	    *fTrackList;              //! List of selected tracks
  TList	    *fTrackListUE;            //! List of tracks in jet cone UE
  TList	    *fTrackListJet;           //! List of tracks in a jet
  TList	    *fCommonHistList;         //! List of common histos
  TH1F      *fh1EvtSelection;         //! event cuts 
  TH1F	    *fh1VertexNContributors;  //! NContributors to prim vertex
  TH1F	    *fh1VertexZ;              //! prim vertex z distribution
  TProfile  *fh1Xsec;                 //! pythia cross section and trials
  TH1F*     fh1Trials;                //! sum of trials
  TH1F*     fh1PtHard;                //! pt hard of the event
  TH1F*     fh1PtHardTrials;          //! pt hard of the event
 
  TH2F*     fh2EtaJet;                //!jet eta distribution
  TH2F*     fh2PhiJet;                //!jet phi distribution
  TH2F*     fh2PtJet;                 //!jet pt distribution
  TH1F*     fh1PtJet;                 //!jet pt distribution 1D
  TH2F*     fh2NtracksJet;            //!number of tracks in jet
  TProfile* fProNtracksJet;           //!number of tracks in jet
  TH2F*     fh2EtaTrack;              //!track eta distribution
  TH2F*     fh2PhiTrack;              //!track phi distribution
  TH2F*     fh2PtTrack;               //!track pt distribution
  TH2F*     fh2FF;                    //!fragmentation function (z)
  TH2F*     fh2Ksi;                   //!fragmentation function (Ksi)
  TH2F*     fh2DelEta;                //!delta eta distribution
  TH2F*     fh2DelPhi;                //!delta phi distribution
  TH2F*     fh2DelR;                  //!delta R distribution
  
  TH1F*     fh1PtLeadingJet;          //!highest jet pt
  TH2F*     fh2NtracksLeadingJet;     //!number of tracks in jet
  TProfile* fProNtracksLeadingJet;    //!number of tracks in jet
  TH2F*     fh2DelR80pcNch;           //!R containing 80% of Nch vs jet pt
  TProfile* fProDelR80pcNch;          //!R containing 80% of Nch vs jet pt
  TH2F*     fh2DelR80pcPt;            //!R containing 80% of pT vs jet pt
  TProfile* fProDelR80pcPt;           //!R containing 80% of pT vs jet pt
  TH2F*     fh2AreaCh;                //!charged jet area vs jet pT
  TProfile* fProAreaCh;               //!charged jet area vs jet pT
  TH3F*     fh3PtDelRNchSum;          //!Nch sum vs R
  TH3F*     fh3PtDelRPtSum;           //!Pt sum vs R
  TProfile* fProDiffJetShape;         //!Diff jet shape pT=20-100
  TProfile* fProIntJetShape;          //!Int jet shape pT=20-100
  TProfile* fProDelRNchSum[13];        //!Nch sum vs R
  TProfile* fProDelRPtSum[13];         //!Pt sum vs R
  TProfile* fProDiffJetShapeA[13];     //!Diff jet shape pT
  TProfile* fProIntJetShapeA[13];      //!Int jet shape pT

  TH1F*     fh1PtSumInJetConeUE;        //!pt sum of UE in cone R
  TH2F*     fh2NtracksLeadingJetUE;     //!number of tracks in jet (UE)
  TProfile* fProNtracksLeadingJetUE;    //!number of tracks in jet (UE)
  TH2F*     fh2DelR80pcNchUE;           //!R containing 80% of Nch vs jet pt (UE)
  TProfile* fProDelR80pcNchUE;          //!R containing 80% of Nch vs jet pt (UE)
  TH2F*     fh2DelR80pcPtUE;            //!R containing 80% of pT vs jet pt (UE)
  TProfile* fProDelR80pcPtUE;           //!R containing 80% of pT vs jet pt (UE)
  TH2F*     fh2AreaChUE;                //!charged jet area vs jet pT (UE)
  TProfile* fProAreaChUE;               //!charged jet area vs jet pT (UE)
  TH3F*     fh3PtDelRNchSumUE;          //!Nch sum vs R (UE)
  TH3F*     fh3PtDelRPtSumUE;           //!Pt sum vs R (UE)
  TProfile* fProDiffJetShapeUE;         //!Diff jet shape pT=20-100 (UE)
  TProfile* fProIntJetShapeUE;          //!Int jet shape pT=20-100 (UE)
  TProfile* fProDelRNchSumUE[13];        //!Nch sum vs R (UE)
  TProfile* fProDelRPtSumUE[13];         //!Pt sum vs R (UE)
  TProfile* fProDiffJetShapeAUE[13];     //!Diff jet shape pT (UE)
  TProfile* fProIntJetShapeAUE[13];      //!Int jet shape pT (UE)

  TH1F* fh1CorrJetPt;                    //!jet pt in FillCorr
  TH2F* fh2CorrPtTrack1;                 //!trk pt in FillCorr
  TH2F* fh2CorrFF1;                      //!FF in FillCorr
  TH2F* fh2CorrKsi1;                     //!Ksi in FillCorr
  TH2F* fh2CorrjT1;                      //!jT in FillCorr
  TH1F* fh1JetPtvsTrkSum;                //!QA plots to check jetPt-trkPtsum
  TH2F* fh2CorrPt1Pt2[6];                //!pt1-pt2 in FillCorr
  TH2F* fh2CorrZ1Z2[6];                  //!z1-z2 in FillCorr
  TH2F* fh2CorrKsi1Ksi2[6];              //!ksi1-ksi2 in FillCorr
  TH2F* fh2CorrjT1jT2[6];                //!jt1-jt2 in FillCorr
  
  AliAnalysisTaskJetProperties(const AliAnalysisTaskJetProperties&);// not implemented
  AliAnalysisTaskJetProperties& operator=(const AliAnalysisTaskJetProperties&);// not implemented
  ClassDef(AliAnalysisTaskJetProperties, 2);
};

#endif
