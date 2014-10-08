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

  void                        SetJetMinRC2LJ(Float_t d)                            { fMinRC2LJ                = d          ; }
  void                        SetRCperEvent(Int_t n)                               { fRCperEvent              = n          ; }
  void                        SetMCJetPtThreshold(Double_t t)                      { fMCJetPtThreshold        = t          ; }
  void                        SetConeRadius(Double_t r)                            { fConeRadius              = r          ; }
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
  void                        GetRandomCone(Float_t &pt, Float_t &eta, Float_t &phi, AliParticleContainer* tracks, AliClusterContainer* clusters,
					    AliEmcalJet *jet = 0, Bool_t bPartialExclusion = 0) const;
  Double_t                    GetZ(const Double_t trkPx, const Double_t trkPy, const Double_t trkPz, const Double_t jetPx, const Double_t jetPy, const Double_t jetPz);
  Double_t                    GetNColl() const;


  Double_t                    fMCJetPtThreshold;                               // threshold for MC jets
  Float_t                     fMinRC2LJ;                                       // Minimum distance random cone to leading jet
  Int_t                       fRCperEvent;                                     // No. of random cones per event
  Double_t                    fConeRadius;                                     // Radius of the random cones
  Float_t                     fConeMinEta;                                     // Minimum eta of the random cones
  Float_t                     fConeMaxEta;                                     // Maximum eta of the random cones
  Float_t                     fConeMinPhi;                                     // Minimum phi of the random cones
  Float_t                     fConeMaxPhi;                                     // Maximum phi of the random cones

  AliJetContainer            *fJetsCont;                                       //!PbPb Jets
  AliParticleContainer       *fTracksCont;                                     //!PbPb Tracks
  AliClusterContainer        *fCaloClustersCont;                               //!PbPb Clusters  
  AliJetContainer            *fMCJetsCont;                                     //!MC jets
  AliParticleContainer       *fMCTracksCont;                                   //!MC tracks
  AliClusterContainer        *fMCCaloClustersCont;                             //!MC clusters  
  AliJetContainer            *fEmbJetsCont;                                    //!EMB jets
  AliParticleContainer       *fEmbTracksCont;                                  //!EMB tracks
  AliClusterContainer        *fEmbCaloClustersCont;                            //!EMB clusters  
  //AliParticleContainer       *fRandTracksCont;                               //!Randomized tracks
  //AliClusterContainer        *fRandCaloClustersCont;                         //!Randomized clusters

  // Random cones
  TH2                        *fHistRCPhiEta;                                   //!Phi-Eta distribution of random cones
  TH1                       **fHistRCPt;                                       //!Random cone pt
  TH1                       **fHistRCPtExLJ;                                   //!Random cone pt, imposing min distance from leading jet
  TH1                       **fHistRCPtExPartialLJ;                            //!Random cone pt, imposing min distance from leading jet with 1/ncoll probability
  //TH1                       **fHistRCPtRand;                                 //!Random cone pt, randomized particles
  TH2                       **fHistRhoVSRCPt;                                  //!Area(RC) * rho vs. Pt(RC)
  TH2                       **fHistDeltaPtRCvsEP;                              //!deltaPt = Pt(RC) - A * rho vs. event plane
  TH1                       **fHistDeltaPtRCExLJ;                              //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet
  TH1                       **fHistDeltaPtRCExPartialLJ;                       //!deltaPt = Pt(RC) - A * rho, imposing min distance from leading jet with 1/ncoll probability
  //TH1                       **fHistDeltaPtRCRand;                            //!deltaPt = Pt(RC) - A * rho, randomzied particles

  // Jet embedding
  TH3                       **fHistEmbJetsPtArea;                              //!Pt vs. area of EMB jets
  TH3                       **fHistEmbJetsCorrPtArea;                          //!Pt-rho*A vs. area of EMB jets
  TH2                       **fHistEmbPartPtvsJetPt;                           //!MC jet pt total jet pt
  TH2                       **fHistEmbPartPtvsJetCorrPt;                       //!MC jet pt total jet pt - rho*A
  TH2                       **fHistJetPtvsJetCorrPt;                           //!Pt vs jet pt - rho*A
  TH1                       **fHistDistLeadPart2JetAxis;                       //!Distance between leading particle and jet axis
  TH2                       **fHistEmbBkgArea;                                 //!Pt(embjet) - Pt(embtrack) vs. area of EMB jets
  TH2                       **fHistRhoVSEmbBkg;                                //!Area(embjet) * rho vs. Pt(embjet) - Pt(embtrack)
  TH2                       **fHistDeltaPtEmbArea;                             //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. Area(embjet)
  TH2                       **fHistDeltaPtEmbvsEP;                             //!deltaPt = Pt(embjet) - Area(embjet) * rho - Pt(embtrack) vs. event plane
  TH2                        *fHistRCPtExLJVSDPhiLJ;                           //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet
  TH2                        *fHistRCPtExPartialLJVSDPhiLJ;                    //!Random cone pt, imposing min distance from leading jet, vs. deltaPhi leading jet with 1/ncoll probability
  TH2                        *fHistEmbJetsPhiEta;                              //!Phi-Eta distribution of EMB jets
  TH2                        *fHistLeadPartPhiEta;                             //!Phi-Eta distribution of the leading particle of EMB jets

  //User Task
  TH1                        *fCent_V0;                                        //!Centrality
  TH1                        *fVertex_z_cut;                                   //!z_vertex_cut
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
  TH1                        *fJet1Pt_PbPb[4][3][4][4];                        //!PbPb, leadingjetPt
  TH1                        *fJet2Pt_PbPb[4][3][4][4];                        //!PbPb, subleadingjetPt
  TH1                        *fJet1Pt_BG_PbPb[4][3][4][4];                     //!PbPb, pT - rho * area, jet1
  TH1                        *fJet2Pt_BG_PbPb[4][3][4][4];                     //!PbPb, pT - rho * area, jet2
  TH1                        *fJetDeltaPhi_PbPb[4][3][4][4];                   //!PbPb, jetDeltaPhi
  TH1                        *fJetDeltaEta_PbPb[4][3][4][4];                   //!PbPb, jetDeltaEta
  TH1                        *fJetDeltaEP_PbPb[4][3][4][4];                    //!PbPb, jetDeltaEP
  TH1                        *fJet1SelectPt_BG_PbPb[4][3][4][4];               //!PbPb, selectleadingjetPt
  TH1                        *fJet2SelectPt_BG_PbPb[4][3][4][4];               //!PbPb, selectsubleadingjetPt
  TH1                        *fAj_PbPb[4][3][4][4];                            //!PbPb, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fJetPt_MC[4][3];                                 //!MC, jetPt
  TH1                        *fJetPhi_MC[4][3];                                //!MC, jetPhi
  TH1                        *fJetEta_MC[4][3];                                //!MC, jetEta
  TH2                        *fJet_Phi_Eta_MC[4][3];                           //!MC, jetPhi vs. jetEta
  TH1                        *fJet1Pt_MC[4][3][4][4];                          //!MC, leadingjetPt
  TH1                        *fJet2Pt_MC[4][3][4][4];                          //!MC, subleadingjetPt
  TH1                        *fJetDeltaPhi_MC[4][3][4][4];                     //!MC, jetDeltaPhi
  TH1                        *fJetDeltaEta_MC[4][3][4][4];                     //!MC, jetDeltaEta
  TH1                        *fJetDeltaEP_MC[4][3][4][4];                      //!MC, jetDeltaEP
  TH1                        *fAj_MC[4][3][4][4];                              //!MC, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fJetPt_EMB[4][3];                                //!EMB, jetPt
  TH1                        *fJetPhi_EMB[4][3];                               //!EMB, jetPhi
  TH1                        *fJetEta_EMB[4][3];                               //!EMB, jetEta
  TH2                        *fJet_Phi_Eta_EMB[4][3];                          //!EMB, jetPhi vs. jetEta
  TH1                        *fJetPt_BG_EMB[4][3];                             //!EMB, pT - rho * area
  TH1                        *fJetDeltaPt[4][3];                               //!EMB, pT - rho * area - pT(embtrack)
  TH1                        *fJet1Pt_EMB[4][3][4][4];                         //!EMB, leadingjetPt
  TH1                        *fJet2Pt_EMB[4][3][4][4];                         //!EMB, subleadingjetPt
  TH1                        *fJet1Pt_BG_EMB[4][3][4][4];                      //!EMB, pT - rho * area, jet1
  TH1                        *fJet2Pt_BG_EMB[4][3][4][4];                      //!EMB, pT - rho * area, jet2
  TH1                        *fJet1DeltaPt[4][3][4][4];                        //!EMB, pT - rho * area - pT(embtrack), jet1
  TH1                        *fJet2DeltaPt[4][3][4][4];                        //!EMB, pT - rho * area - pT(embtrack), jet2
  TH1                        *fJetDeltaPhi_EMB[4][3][4][4];                    //!EMB, jetDeltaPhi
  TH1                        *fJetDeltaEta_EMB[4][3][4][4];                    //!EMB, jetDeltaEta
  TH1                        *fJetDeltaEP_EMB[4][3][4][4];                     //!EMB, jetDeltaEP
  TH1                        *fJet1SelectPt_BG_EMB[4][3][4][4];                //!EMB, selectleadingjetPt
  TH1                        *fJet2SelectPt_BG_EMB[4][3][4][4];                //!EMB, selectsubleadingjetPt
  TH1                        *fJet1SelectDeltaPt[4][3][4][4];                  //!EMB, selectleadingjetPt
  TH1                        *fJet2SelectDeltaPt[4][3][4][4];                  //!EMB, selectsubleadingjetPt
  TH1                        *fAj_EMB[4][3][4][4];                             //!EMB, Aj(energy balance) -> Aj = (jet1-jet2)/(jet1+jet2)

  TH1                        *fHJetDeltaPhi_Aj0_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_PbPb[4][3][4][4][4];           //!PbPb, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJetPt_Aj0_PbPb[4][3][4][4][4];                 //!PbPb, HjetPt, no Aj cut
  TH1                        *fHJetPt_Aj1_PbPb[4][3][4][4][4];                 //!PbPb, HjetPt, Aj1
  TH1                        *fHJetPt_Aj2_PbPb[4][3][4][4][4];                 //!PbPb, HjetPt, Aj2
  TH1                        *fHJetPt_Aj3_PbPb[4][3][4][4][4];                 //!PbPb, HjetPt, Aj3
  TH1                        *fHJetPt_Aj4_PbPb[4][3][4][4][4];                 //!PbPb, HjetPt, Aj4
  TH1                        *fHJetDeltaPhi_Aj0_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_MC[4][3][4][4][4];             //!MC, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJetPt_Aj0_MC[4][3][4][4][4];                   //!MC, HjetPt, no Aj cut
  TH1                        *fHJetPt_Aj1_MC[4][3][4][4][4];                   //!MC, HjetPt, Aj1
  TH1                        *fHJetPt_Aj2_MC[4][3][4][4][4];                   //!MC, HjetPt, Aj2
  TH1                        *fHJetPt_Aj3_MC[4][3][4][4][4];                   //!MC, HjetPt, Aj3
  TH1                        *fHJetPt_Aj4_MC[4][3][4][4][4];                   //!MC, HjetPt, Aj4
  TH1                        *fHJetDeltaPhi_Aj0_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, no Aj cut
  TH1                        *fHJetDeltaPhi_Aj1_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj1(0.0 to 0.2)
  TH1                        *fHJetDeltaPhi_Aj2_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj2(0.2 to 0.4)
  TH1                        *fHJetDeltaPhi_Aj3_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj3(0.4 to 0.6)
  TH1                        *fHJetDeltaPhi_Aj4_EMB[4][3][4][4][4];            //!EMB, HjetDeltaPhi, Aj4(0.6 to 0.8)
  TH1                        *fHJetPt_Aj0_EMB[4][3][4][4][4];                  //!EMB, HjetPt, no Aj cut
  TH1                        *fHJetPt_Aj1_EMB[4][3][4][4][4];                  //!EMB, HjetPt, Aj1
  TH1                        *fHJetPt_Aj2_EMB[4][3][4][4][4];                  //!EMB, HjetPt, Aj2
  TH1                        *fHJetPt_Aj3_EMB[4][3][4][4][4];                  //!EMB, HjetPt, Aj3
  TH1                        *fHJetPt_Aj4_EMB[4][3][4][4][4];                  //!EMB, HjetPt, Aj4

  TH1                        *fHJetDeltaPhiasEP_Aj0_PbPb[4][4][4][4][4];       //!PbPb, HjetDeltaPhi, asEP, no Aj cut
  TH1                        *fHJetDeltaPhiasEP_Aj1_PbPb[4][4][4][4][4];       //!PbPb, HjetDeltaPhi, asEP, Aj1
  TH1                        *fHJetDeltaPhiasEP_Aj2_PbPb[4][4][4][4][4];       //!PbPb, HjetDeltaPhi, asEP, Aj2
  TH1                        *fHJetDeltaPhiasEP_Aj3_PbPb[4][4][4][4][4];       //!PbPb, HjetDeltaPhi, asEP, Aj3
  TH1                        *fHJetDeltaPhiasEP_Aj4_PbPb[4][4][4][4][4];       //!PbPb, HjetDeltaPhi, asEP, Aj4
  TH1                        *fHJetPtasEP_Aj0_PbPb[4][4][4][4][4];             //!PbPb, HjetPt, asEP, no Aj cut
  TH1                        *fHJetPtasEP_Aj1_PbPb[4][4][4][4][4];             //!PbPb, HjetPt, asEP, Aj1
  TH1                        *fHJetPtasEP_Aj2_PbPb[4][4][4][4][4];             //!PbPb, HjetPt, asEP, Aj2
  TH1                        *fHJetPtasEP_Aj3_PbPb[4][4][4][4][4];             //!PbPb, HjetPt, asEP, Aj3
  TH1                        *fHJetPtasEP_Aj4_PbPb[4][4][4][4][4];             //!PbPb, HjetPt, asEP, Aj4
  TH1                        *fHJetDeltaPhiasEP_Aj0_MC[4][4][4][4][4];         //!MC, HjetDeltaPhi, asEP, no Aj cut
  TH1                        *fHJetDeltaPhiasEP_Aj1_MC[4][4][4][4][4];         //!MC, HjetDeltaPhi, asEP, Aj1
  TH1                        *fHJetDeltaPhiasEP_Aj2_MC[4][4][4][4][4];         //!MC, HjetDeltaPhi, asEP, Aj2
  TH1                        *fHJetDeltaPhiasEP_Aj3_MC[4][4][4][4][4];         //!MC, HjetDeltaPhi, asEP, Aj3
  TH1                        *fHJetDeltaPhiasEP_Aj4_MC[4][4][4][4][4];         //!MC, HjetDeltaPhi, asEP, Aj4
  TH1                        *fHJetPtasEP_Aj0_MC[4][4][4][4][4];               //!MC, HjetPt, asEP, no Aj cut
  TH1                        *fHJetPtasEP_Aj1_MC[4][4][4][4][4];               //!MC, HjetPt, asEP, Aj1
  TH1                        *fHJetPtasEP_Aj2_MC[4][4][4][4][4];               //!MC, HjetPt, asEP, Aj2
  TH1                        *fHJetPtasEP_Aj3_MC[4][4][4][4][4];               //!MC, HjetPt, asEP, Aj3
  TH1                        *fHJetPtasEP_Aj4_MC[4][4][4][4][4];               //!MC, HjetPt, asEP, Aj4
  TH1                        *fHJetDeltaPhiasEP_Aj0_EMB[4][4][4][4][4];        //!EMB, HjetDeltaPhi, asEP, no Aj cut
  TH1                        *fHJetDeltaPhiasEP_Aj1_EMB[4][4][4][4][4];        //!EMB, HjetDeltaPhi, asEP, Aj1
  TH1                        *fHJetDeltaPhiasEP_Aj2_EMB[4][4][4][4][4];        //!EMB, HjetDeltaPhi, asEP, Aj2
  TH1                        *fHJetDeltaPhiasEP_Aj3_EMB[4][4][4][4][4];        //!EMB, HjetDeltaPhi, asEP, Aj3
  TH1                        *fHJetDeltaPhiasEP_Aj4_EMB[4][4][4][4][4];        //!EMB, HjetDeltaPhi, asEP, Aj4
  TH1                        *fHJetPtasEP_Aj0_EMB[4][4][4][4][4];              //!EMB, HjetPt, asEP, no Aj cut
  TH1                        *fHJetPtasEP_Aj1_EMB[4][4][4][4][4];              //!EMB, HjetPt, asEP, Aj1
  TH1                        *fHJetPtasEP_Aj2_EMB[4][4][4][4][4];              //!EMB, HjetPt, asEP, Aj2
  TH1                        *fHJetPtasEP_Aj3_EMB[4][4][4][4][4];              //!EMB, HjetPt, asEP, Aj3
  TH1                        *fHJetPtasEP_Aj4_EMB[4][4][4][4][4];              //!EMB, HjetPt, asEP, Aj4


 private:
  AliVEvent                  *fEvent;
  Double_t                    fCentrality;                                     //! V0M for current event
  //AliNamedString             *fPtHardBinName;                                //!Pt hard bin param
  //Int_t                       fPtHardBin;                                    //!        
  //TH1F                        *fhPtHardBins;                                 //!

  AliAnalysisTaskDijetHadron(const AliAnalysisTaskDijetHadron&);                     // not implemented
  AliAnalysisTaskDijetHadron &operator=(const AliAnalysisTaskDijetHadron&);          // not implemented

  ClassDef(AliAnalysisTaskDijetHadron, 5)                                         // Jet-Hadron analysis task
};
#endif
