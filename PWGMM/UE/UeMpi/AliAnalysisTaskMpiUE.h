/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* Add a description of your MPI analysis */

#ifndef AliAnalysisTaskMpiUE_H
#define AliAnalysisTaskMpiUE_H

class AliESDtrackCuts;
class AliESDEvent;
class TList;
class TH1D;
class TH2D;
class TH3D;
class TH1I;
class TProfile;
class THnSparse;

#include "AliAnalysisTaskSE.h"


#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliGenEventHeader.h"



class AliAnalysisTaskMpiUE : public AliAnalysisTaskSE
{
public:
  AliAnalysisTaskMpiUE();
  AliAnalysisTaskMpiUE(const char *name);
  virtual                 ~AliAnalysisTaskMpiUE();

  virtual void            UserCreateOutputObjects();
  virtual void            UserExec(Option_t* option);
  virtual void            Terminate(Option_t* option);

  void   TracksLoop(AliESDEvent *);
  static Bool_t IsMinimumBias(AliVEvent *);
  static Bool_t IsINELgtZERO(AliVEvent *);
  static Bool_t HasAcceptedVertexPosition(AliVEvent *);
  static Bool_t HasNoPileupSPDInMultBins(AliVEvent *);
  static Bool_t HasNoInconsistentSPDandTrackVertices(AliVEvent *);
  static Bool_t IsEventSelected(AliVEvent *);
  static Bool_t IsInCompleteEvent(AliVEvent *);
  static Bool_t IsSPDClusterVsTrackletBG(AliVEvent *);
  static Bool_t IsOutOfBunchPileup(AliVEvent *);
  static Bool_t IsOutOfBunchPileupFastor(AliVEvent *);
  static Bool_t HasBCMod4(AliVEvent *);

  static Bool_t V0Decision(AliVEvent *);
  static Bool_t V0Asymmetry(AliVEvent *);

  void       SetUseMC(Bool_t mc = kFALSE)              {fUseMC = mc;}   // use to analyse MC data
  void       SetTrackingEfficiency(TH1D* hist) {fTrackingEfficiency = hist;}       // Set the tracking efficiency correction
  virtual Double_t deltaPhi(Double_t phi1, Double_t phi2);                               // difference in the phi of GenLeadPart and RecLeadtrack
  virtual Double_t DeltaPhi(Double_t phi1, Double_t phi2);                               // difference in the phi of lead and track
  virtual Double_t DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);   // distance between lead and track

protected:

  void       FillMCGen();
  


private:
  AliESDtrackCuts*        fCuts;                                       //!
  AliESDtrackCuts*        fCuts2;                                      //!
  AliESDEvent*            fESD;                                        //! input ESD event
  AliStack*    fStack;                                                 //! MC stack
  AliMCEvent*  fMCEvent;                                               //! MC Event
  Bool_t       fUseMC;                // analyze MC events


  TList*                  fOutputList;                                      //! output list in the root file
  TH1I*                   fEventCounter;                                    //! Event Summary with cuts


  TH1I*                   fHistEvents;                 //!

  TH1D*                   fHistAllTracksDCAxy;         //!
  TH1D*                   fHistPrimaryDCAxy;           //!
  TH1D*                   fHistStrangeDCAxy;           //!
  TH1D*                   fHistPionDecayDCAxy;         //!
  TH1D*                   fHistMuonDecayDCAxy;         //!
  TH1D*                   fHistMaterialDCAxy;          //!
  TH1D*                   fHistOtherDCAxy;             //!

  TH1D*                   fHistAllTracksPt;            //!
  TH1D*                   fHistPrimaryPt;              //!
  TH1D*                   fHistStrangePt;              //!
  TH1D*                   fHistMaterialPt;             //!

  TH1D*                   fHistDCAxy;                  //!  for primaries in central DCAxy regions |DCAxy|<0.3cm



  TH1I*                   fHistMultiplicityStd;                                  //! Std. Reference Multiplicity distribution
  TH1D*                   fHistPt;                                               //! pT distribution of Reconstructed tracks
  TH1D*                   fHistLeadPt;                                           //! Rec Lead pT distribution of tracks
  TH1D*                   fHistMatchLeadPt;                                      //! Rec Lead pT distribution of tracks matched with Gen lead tracks
  TH1D*                   fHistNoMatchLeadPt;                                    //! Rec Lead pT distribution of tracks not matched with Gen lead tracks
  TH1D*                   fHistSumPt;                                            //! Rec Sum pT distribution of tracks
  TH1D*                   fHistAvgPt;                                            //! Avg pT distribution of tracks
  TH1D*                   fHistEta;                                              //! Rec Eta distribution of tracks
  TH1D*                   fHistPhi;                                              //! Rec Phi distribution of tracks
  TH2D*                   fHistEta_Phi;                                          //! Rec Eta-Phi distribution of tracks
    TH2D*                   fHistEta_Pt;                                          //! Rec Eta-Pt distribution of tracks
  TH1D*                   fHistZvtx;                                             //! Rec Vertex Z
  TH1D*                   fHistZvtxCut;                                          //! Rec Vertex Z With Cut (-10,10)
  TH1D*                   fHistTracks;                                           //! Number distribution of track
  TH2D*                   fHistZvtxMult;                                         //! with vertex cut  (zvtx-mult)
  TH2D*                   fHistZvtxMultTrue;                                     //! with vertex cut  (zvtx-multTrue)
  TH2D*                   fHistZvtxMultNcon1;                                    //! with vertex cut & Ncontribytors > 1  (zvtx-mult)
  TH2D*                   fHistZvtxMultNcon1True;                               //! with vertex cut & Ncontribytors > 1  (zvtx-multTrue)
  TH2D*                   fHistZvtxMultNcon2;                                    //! with vertex cut & Ncontribytors > 2  (zvtx-mult)

  TH1D*                   fHistGenPt;                                            //! pT distribution of Generated tracks
  TH1D*                   fHistGenLeadPt;                                        //! Gen Lead pT distribution of tracks
  TH1D*                   fHistGenSumPt;                                         //! Gen Sum pT distribution of tracks
  TH1D*                   fHistGenAvgPt;                                         //! Gen Avg pT distribution of tracks
  TH1D*                   fHistGenEta;                                           //! Eta distribution of Generated tracks
  TH1D*                   fHistGenPhi;                                           //! Phi distribution of Generated tracks
  TH1D*                   fHistGenZvtx;                                          //! Gen vertex Z
  TH1D*                   fHistGenZvtxCut;                                       //! Gen vertex Z with Cut (-10,10)
  TH1D*                   fHistGenTracks;                                        //! Number Generated tracks
  TH2D*                   fHistGenZvtxMult;                                      //!
  TH2D*                   fHistMultResponseMat;                                  //!

  THnSparse*                   fHistTrack;                                                 //! pT,eta,Phi distribution of tracks
  THnSparse*                   fHistLeadTrack;                                             //! Lead pT,Sum pT , Avg pT distribution of tracks
  THnSparse*                   fHistGenTrack;                                              //!Gen pT,eta,Phi distribution of tracks
  THnSparse*                   fHistGenLeadTrack;                                          //! Gen Lead pT,Sum pT , Avg pT distribution of tracks

  TH1D*                       fHistLeadDPhi;                                  //!
  TH1D*                       fHistLeadDPhiTr;                                //!
  TH1D*                       fHistLeadDPhiTo;                                //!
  TH1D*                       fHistLeadDPhiAw;                                //!
  TProfile*                   fHistMatchParticleDensityTransverse;            //!
  TProfile*                   fHistMatchEnergyDensityTransverse;              //!

  TProfile*                   fHistParticleDensityTowards;                    //!
  TProfile*                   fHistParticleDensityTransverse;                 //!
  TProfile*                   fHistParticleDensityAway;                       //!
  TProfile*                   fHistEnergyDensityTowards;                      //!
  TProfile*                   fHistEnergyDensityTransverse;                   //!
  TProfile*                   fHistEnergyDensityAway;                         //!

  TProfile*                   fHistGenParticleDensityTowards;                 //!
  TProfile*                   fHistGenParticleDensityTransverse;              //!
  TProfile*                   fHistGenParticleDensityAway;                    //!
  TProfile*                   fHistGenEnergyDensityTowards;                   //!
  TProfile*                   fHistGenEnergyDensityTransverse;                //!
  TProfile*                   fHistGenEnergyDensityAway;                      //!
  TF1*                        f1;                                             //!
  TH1D*                       fTrackingEfficiency;                           

  THnSparse               *fHnDelta;                                                       //! multi-dimentiona histogram for delta properties
  THnSparse               *fHnDeltaTo;                                                     //! multi-dimentiona histogram for Towardss delta properties
  THnSparse               *fHnDeltaTr;                                                     //! multi-dimentiona histogram for Transvere delta properties
  THnSparse               *fHnDeltaAw;                                                     //! multi-dimentiona histogram for Away delta properties

  TH2D*                   fHistDPhiEta;                                                    //! DeltaEta-DeltaPhi of the Leading track and associated track
  TH2D*                   fHistGenDPhiEta;                                                 //! DeltaEta-DeltaPhi for Genenated particles of the Leading track and associated track


    TH1D*                   fHistGenMultTrans;                                        //! Number Generated tracks in Transvere
    TH1D*                   fHistGenMultTowds;
    TH1D*                   fHistGenMultAway;
    TH1D*                   fHistGenEngyTrans;
    TH1D*                   fHistGenEngyTowds;
    TH1D*                   fHistGenEngyAway;
    TH1D*                   fHistMultTrans;                                        //! Number of tracks in Transvere
    TH1D*                   fHistMultTowds;
    TH1D*                   fHistMultAway;
    TH1D*                   fHistEngyTrans;
    TH1D*                   fHistEngyTowds;
    TH1D*                   fHistEngyAway;


  AliMCParticle*          fGenLeadPart;                                         //!
  Double_t                fEtaCut;                                             //!
  Double_t                fPtMin;                                              //!
    Float_t                fLeadPtCutMin;                                             //!
    Float_t                fLeadPtCutMax;                                             //!

    
  Double_t                 fVertexMC[3];                  //! Gen vertex X,Y,Z
  Float_t                 fGenNTracks;                                         //
  Float_t                 fRecNTracks;                                         //

  Double_t                fGenLeadPhi;
  Double_t                fRecLeadPhi;


  AliAnalysisTaskMpiUE(const AliAnalysisTaskMpiUE&);                  // not implemented
  AliAnalysisTaskMpiUE& operator=(const AliAnalysisTaskMpiUE&);       // not implemented

  ClassDef(AliAnalysisTaskMpiUE, 3);
};

#endif
