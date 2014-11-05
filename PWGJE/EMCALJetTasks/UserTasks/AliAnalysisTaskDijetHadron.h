#ifndef ALIANALYSISTASKDIJETHADRON_H
#define ALIANALYSISTASKDIJETHADRON_H

// $Id$
#include <vector>

class TH1;
class TH2;
class TH3;
class TH2F;
class TH1F;
class TF1;
class TH3F;
class THnSparse;
class TClonesArray;
class TObject;
class TString;
class AliNamedString;
class AliAODEvent;
class AliESDEvent;
class AliMCEvent;
class AliRhoParameter;
class TRandom3;
class AliEmcalJet;
class AliVTrack;
class AliNamedArrayI;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;

template<class T> class TParameter;

#include "AliAnalysisTaskSE.h"
#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDijetHadron : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskDijetHadron();
  AliAnalysisTaskDijetHadron(const char *name);
  virtual ~AliAnalysisTaskDijetHadron() {;}

  void                        UserCreateOutputObjects();

  void                        SetMCJetPtThreshold(Double_t t)                      { fMCJetPtThreshold        = t          ; }
  void                        SetLeadingHadronPtThreshold1(Double_t u1)            { fleadingHadronPtcut1     = u1         ; }
  void                        SetLeadingHadronPtThreshold2(Double_t u2)            { fleadingHadronPtcut2     = u2         ; }
  void                        SetLeadingHadronPtThreshold3(Double_t u3)            { fleadingHadronPtcut3     = u3         ; }
  void                        SetJet1PtThreshold1(Double_t v1)                     { fJet1Ptcut1     = v1                  ; }
  void                        SetJet1PtThreshold2(Double_t v2)                     { fJet1Ptcut2     = v2                  ; }
  void                        SetJet1PtThreshold3(Double_t v3)                     { fJet1Ptcut3     = v3                  ; }
  void                        SetJet2PtThreshold1(Double_t w1)                     { fJet2Ptcut1     = w1                  ; }
  void                        SetJet2PtThreshold2(Double_t w2)                     { fJet2Ptcut2     = w2                  ; }
  void                        SetJet2PtThreshold3(Double_t w3)                     { fJet2Ptcut3     = w3                  ; }
  void                        SetConeRadius(Double_t r)                            { fConeRadius     = r                   ; }
  void                        SetConeEtaPhiEMCAL() ;
  void                        SetConeEtaPhiTPC()   ;
  void                        SetConeEtaLimits(Float_t min, Float_t max)           { fConeMinEta = min, fConeMaxEta = max  ; }
  void                        SetConePhiLimits(Float_t min, Float_t max)           { fConeMinPhi = min, fConeMaxPhi = max  ; }

 protected:
  void                        AllocateHistogramArrays()                                                                     ;
  void                        ExecOnce()                                                                                    ;
  Bool_t                      FillHistograms()                                                                              ;
  void                        GetLeadingJets(Int_t &maxJetIndex, Int_t &max2JetIndex)                                       ;
  AliEmcalJet*                NextEmbeddedJet(Bool_t reset=kFALSE)                                                          ;
  void                        DoEmbTrackLoop()                                                                              ;
  void                        DoEmbClusterLoop()                                                                            ;
  Double_t                    GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz);
  Double_t                    GetDPhi(Double_t mphi,Double_t vphi);

  Double_t                    fMCJetPtThreshold;                               // threshold for MC jets
  Double_t                    fleadingHadronPtcut1;                            // threshold for leading hadron pT NO1
  Double_t                    fleadingHadronPtcut2;                            // threshold for leading hadron pT NO2
  Double_t                    fleadingHadronPtcut3;                            // threshold for leading hadron pT NO3
  Double_t                    fJet1Ptcut1;                                     // threshold for leading Jet pT NO1
  Double_t                    fJet1Ptcut2;                                     // threshold for leading Jet pT NO2
  Double_t                    fJet1Ptcut3;                                     // threshold for leading Jet pT NO3
  Double_t                    fJet2Ptcut1;                                     // threshold for subleading Jet pT NO1
  Double_t                    fJet2Ptcut2;                                     // threshold for subleading Jet pT NO2
  Double_t                    fJet2Ptcut3;                                     // threshold for subleading Jet pT NO3
  Double_t                    fConeRadius;                                     // Radius of the jet cones
  Float_t                     fConeMinEta;                                     // Minimum eta of the jet cones
  Float_t                     fConeMaxEta;                                     // Maximum eta of the jet cones
  Float_t                     fConeMinPhi;                                     // Minimum phi of the jet cones
  Float_t                     fConeMaxPhi;                                     // Maximum phi of the jet cones

  AliJetContainer            *fJetsCont;                                       //!PbPb Jets
  AliParticleContainer       *fTracksCont;                                     //!PbPb Tracks
  AliClusterContainer        *fCaloClustersCont;                               //!PbPb Clusters  
  AliJetContainer            *fMCJetsCont;                                     //!MC jets
  AliParticleContainer       *fMCTracksCont;                                   //!MC tracks
  AliClusterContainer        *fMCCaloClustersCont;                             //!MC clusters  
  AliJetContainer            *fEmbJetsCont;                                    //!EMB jets
  AliParticleContainer       *fEmbTracksCont;                                  //!EMB tracks
  AliClusterContainer        *fEmbCaloClustersCont;                            //!EMB clusters  

  //User Task
  TH1                        *fCent_V0;                                        //!Centrality
  TH1                        *fVertex_z_cut;                                   //!z_vertex_cut
  TH1                        *fEP2;                                   //!z_vertex_cut
  TH1                        *fJetBG_rho;                                      //!rhoValue
  TH2                        *fJetBG_rho_Cent;                                 //!rho vs. Centrality
  TH1                        **fTrackPt_PbPb;                                  //!PbPb, trackPt
  TH1                        **fTrackPhi_PbPb;                                 //!PbPb, trackPhi
  TH1                        **fTrackEta_PbPb;                                 //!PbPb, trackEta
  TH2                        **fTrack_Phi_Eta_PbPb;                            //!PbPb, trackPhi vs. trackEta
  TH1                        **fTrackPt_MC;                                    //!MC, trackPt
  TH1                        **fTrackPhi_MC;                                   //!MC, trackPhi
  TH1                        **fTrackEta_MC;                                   //!MC, trackEta
  TH2                        **fTrack_Phi_Eta_MC;                              //!MC, trackPhi vs. trackEta
  TH1                        **fTrackPt_EMB;                                   //!EMB, trackPt
  TH1                        **fTrackPhi_EMB;                                  //!EMB, trackPhi
  TH1                        **fTrackEta_EMB;                                  //!EMB, trackEta
  TH2                        **fTrack_Phi_Eta_EMB;                             //!EMB, trackPhi vs. trackEta
  TH1                        *fJetPt_PbPb[4][3];                               //!PbPb, jetPt
  TH1                        *fJetPhi_PbPb[4][3];                              //!PbPb, jetPhi
  TH1                        *fJetEta_PbPb[4][3];                              //!PbPb, jetEta
  TH2                        *fJet_Phi_Eta_PbPb[4][3];                         //!PbPb, jetPhi vs. jetEta
  TH1                        *fJetPt_BG_PbPb[4][3];                            //!PbPb, pT - rho * area
  TH2                        *fJetDeltaEP_PbPb[4][3];                          //!PbPb, jetDeltaEP
  TH1                        *fJet1Pt_PbPb[4][3][4][4];                        //!PbPb, leadingjetPt
  TH1                        *fJet2Pt_PbPb[4][3][4][4];                        //!PbPb, subleadingjetPt
  TH1                        *fJet1Pt_BG_PbPb[4][3][4][4];                     //!PbPb, pT - rho * area, jet1
  TH1                        *fJet2Pt_BG_PbPb[4][3][4][4];                     //!PbPb, pT - rho * area, jet2
  TH1                        *fJetDeltaPhi_PbPb[4][3][4][4];                   //!PbPb, jetDeltaPhi
  TH1                        *fJetDeltaEta_PbPb[4][3][4][4];                   //!PbPb, jetDeltaEta
  TH1                        *fJet1SelectPt_BG_PbPb[4][3][4][4];               //!PbPb, selectleadingjetPt
  TH1                        *fJet2SelectPt_BG_PbPb[4][3][4][4];               //!PbPb, selectsubleadingjetPt
  TH2                        *fJet1EP_PbPb[4][3][4][4];                        //!PbPb, jet1DeltaEP
  TH1                        *fAj_PbPb[4][3][4][4];                            //!PbPb, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fJetPt_MC[4][3];                                 //!MC, jetPt
  TH1                        *fJetPhi_MC[4][3];                                //!MC, jetPhi
  TH1                        *fJetEta_MC[4][3];                                //!MC, jetEta
  TH2                        *fJet_Phi_Eta_MC[4][3];                           //!MC, jetPhi vs. jetEta
  TH2                        *fJetDeltaEP_MC[4][3];                            //!MC, jetDeltaEP
  TH1                        *fJet1Pt_MC[4][3][4][4];                          //!MC, leadingjetPt
  TH1                        *fJet2Pt_MC[4][3][4][4];                          //!MC, subleadingjetPt
  TH1                        *fJetDeltaPhi_MC[4][3][4][4];                     //!MC, jetDeltaPhi
  TH1                        *fJetDeltaEta_MC[4][3][4][4];                     //!MC, jetDeltaEta
  TH2                        *fJet1EP_MC[4][3][4][4];                          //!MC, jet1DeltaEP
  TH1                        *fAj_MC[4][3][4][4];                              //!MC, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fJetPt_EMB[4][3];                                //!EMB, jetPt
  TH1                        *fJetPhi_EMB[4][3];                               //!EMB, jetPhi
  TH1                        *fJetEta_EMB[4][3];                               //!EMB, jetEta
  TH2                        *fJet_Phi_Eta_EMB[4][3];                          //!EMB, jetPhi vs. jetEta
  TH1                        *fJetPt_BG_EMB[4][3];                             //!EMB, pT - rho * area
  TH1                        *fJetDeltaPt[4][3];                               //!EMB, pT - rho * area - pT(embtrack)
  TH2                        *fJetDeltaEP_EMB[4][3];                           //!EMB, jetDeltaEP
  TH1                        *fJet1Pt_EMB[4][3][4][4];                         //!EMB, leadingjetPt
  TH1                        *fJet2Pt_EMB[4][3][4][4];                         //!EMB, subleadingjetPt
  TH1                        *fJet1Pt_BG_EMB[4][3][4][4];                      //!EMB, pT - rho * area, jet1
  TH1                        *fJet2Pt_BG_EMB[4][3][4][4];                      //!EMB, pT - rho * area, jet2
  TH1                        *fJet1DeltaPt[4][3][4][4];                        //!EMB, pT - rho * area - pT(embtrack), jet1
  TH1                        *fJet2DeltaPt[4][3][4][4];                        //!EMB, pT - rho * area - pT(embtrack), jet2
  TH1                        *fJetDeltaPhi_EMB[4][3][4][4];                    //!EMB, jetDeltaPhi
  TH1                        *fJetDeltaEta_EMB[4][3][4][4];                    //!EMB, jetDeltaEta
  TH1                        *fJet1SelectPt_BG_EMB[4][3][4][4];                //!EMB, selectleadingjetPt
  TH1                        *fJet2SelectPt_BG_EMB[4][3][4][4];                //!EMB, selectsubleadingjetPt
  TH1                        *fJet1SelectDeltaPt[4][3][4][4];                  //!EMB, selectleadingjetPt
  TH1                        *fJet2SelectDeltaPt[4][3][4][4];                  //!EMB, selectsubleadingjetPt
  TH2                        *fJet1EP_EMB[4][3][4][4];                         //!EMB, jet1DeltaEP
  TH1                        *fAj_EMB[4][3][4][4];                             //!EMB, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fHJetDeltaPhi_Aj0_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJet_EP_Aj0_PbPb[4][3][4][4][4];                //!PbPb, HjetEP, no Aj cut
  TH1                        *fHJet_EP_Aj1_PbPb[4][3][4][4][4];                //!PbPb, HjetEP, Aj1
  TH1                        *fHJet_EP_Aj2_PbPb[4][3][4][4][4];                //!PbPb, HjetEP, Aj2
  TH1                        *fHJet_EP_Aj3_PbPb[4][3][4][4][4];                //!PbPb, HjetEP, Aj3
  TH1                        *fHJet_EP_Aj4_PbPb[4][3][4][4][4];                //!PbPb, HjetEP, Aj4
  TH1                        *fHJetDeltaPhi_Aj0_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJet_EP_Aj0_MC[4][3][4][4][4];                  //!MC, HjetEP, no Aj cut
  TH1                        *fHJet_EP_Aj1_MC[4][3][4][4][4];                  //!MC, HjetEP, Aj1
  TH1                        *fHJet_EP_Aj2_MC[4][3][4][4][4];                  //!MC, HjetEP, Aj2
  TH1                        *fHJet_EP_Aj3_MC[4][3][4][4][4];                  //!MC, HjetEP, Aj3
  TH1                        *fHJet_EP_Aj4_MC[4][3][4][4][4];                  //!MC, HjetEP, Aj4
  TH1                        *fHJetDeltaPhi_Aj0_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJet_EP_Aj0_EMB[4][3][4][4][4];                 //!EMB, HjetEP, no Aj cut
  TH1                        *fHJet_EP_Aj1_EMB[4][3][4][4][4];                 //!EMB, HjetEP, Aj1
  TH1                        *fHJet_EP_Aj2_EMB[4][3][4][4][4];                 //!EMB, HjetEP, Aj2
  TH1                        *fHJet_EP_Aj3_EMB[4][3][4][4][4];                 //!EMB, HjetEP, Aj3
  TH1                        *fHJet_EP_Aj4_EMB[4][3][4][4][4];                 //!EMB, HjetEP, Aj4

 private:
  AliVEvent                  *fEvent;
  Double_t                    fCentrality;                                     //! V0M for current event

  AliAnalysisTaskDijetHadron(const AliAnalysisTaskDijetHadron&);                     // not implemented
  AliAnalysisTaskDijetHadron &operator=(const AliAnalysisTaskDijetHadron&);          // not implemented

  ClassDef(AliAnalysisTaskDijetHadron, 5)                                         // Jet-Hadron analysis task
};
#endif
