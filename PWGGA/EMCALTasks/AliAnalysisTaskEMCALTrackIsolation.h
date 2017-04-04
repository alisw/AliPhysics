#ifndef ALIANALYSISTASKEMCALTRACKISOLATION_H
#define ALIANALYSISTASKEMCALTRACKISOLATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

  ///////////////////////////////////////////////////////////////////////////
  ///\class AliAnalysisTaskEMCALTrackIsolation
  ///\brief Task for Isolated Gamma in p-p,p-Pb and eventually g-h Correlation
  ///
  /// Task designed to study isolated photon using different methods to study underlying events and perform isolation
  /// It is based on the jet framework
  ///
  /// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
  ///////////////////////////////////////////////////////////////////////////

  //ROOT System
class TH1D;
class TH2D;
class TH3D;
class TF1;
class THnSparse;
class TList;
class TObjArray;
class AliEMCALGeometry;
class AliESDEvent;
class AliESDtrack;
class TClonesArray;
class TList;
class TString;
class AliVParticle;
class AliESDtrackCuts;
class AliAODEvent;
class AliMCEvent;
class AliStack;
class TParticle;
class AliTrackerContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliTrackContainer;
class AliEmcalParticle;
  //AliRoot Syste
class AliEMCALTrack;
class AliVCluster;
class AliAODCaloCluster;
  //class AliMagF;
class AliEMCALRecoUtils;
  //class AliAnalysisFilter;
class AliAODTrack;
class AliAODMCParticle;
class AliGenPythiaEventHeader;
class AliAODMCHeader;
  //class AliEventPoolManager;

#include "AliAnalysisTaskEmcal.h"
#include <vector>
using std::vector;

class AliAnalysisTaskEMCALTrackIsolation: public AliAnalysisTaskEmcal {
public:
  AliAnalysisTaskEMCALTrackIsolation();
  AliAnalysisTaskEMCALTrackIsolation(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskEMCALTrackIsolation();
  
  void         UserCreateOutputObjects();
  
  void         SetIsoConeRadius(Float_t r)                    { fIsoConeRadius = r ;}
  void         SetPtIsoThreshold(Float_t r)                   { fptIsoThreshold = r ;}
  void         SetIsoMethod(Int_t r)                          { fIsoMethod = r ;}
  void         SetPtIsoMethod (Int_t r )                      { fptIsoMethod = r ;}
  void         SetUEMethod (Int_t r )                         { fUEMethod = r ;}
  void         SetQA (Bool_t QA)                              { fQA = QA;}
  void         SetMC (Bool_t MC)                              { fIsMC = MC;}
  void         SetUSEofTPC (Bool_t TPC)                       { fTPC4Iso = TPC;}
  void         SetLTAnalysis (Bool_t LT)                      { fisLTAnalysis = LT;}
  void         SetRejectEventWithoutTracks(Bool_t revwotr)    { fRejectionEventWithoutTracks = revwotr;}
  void         SetAnalysispPb(Bool_t ana)                     { fAnalysispPb = ana;}
  void         SetTriggerLevel1(Int_t L)                      { fTriggerLevel1 = L;}
  void			   SetPtBinning(vector<Double_t> binedges)        { fBinsPt = binedges; }
  void			   SetPtisoBinning(vector<Double_t> binedges)			{ fBinsPtiso = binedges; }
  void			   SetPtueBinning(vector<Double_t> binedges)			{ fBinsPtue = binedges; }
  void			   SetEtaBinning(vector<Double_t> binedges)			  { fBinsEta = binedges; }
  void			   SetPhiBinning(vector<Double_t> binedges)			  { fBinsPhi = binedges; }
  void			   SetLabelBinning(vector<Double_t> binedges)			{ fBinsLabel = binedges; }
  void			   SetPDGBinning(vector<Double_t> binedges)			  { fBinsPDG = binedges; }
  void			   SetMomPDGBinning(vector<Double_t> binedges)		{ fBinsMomPDG = binedges; }
  void			   SetTrackPDGBinning(vector<Double_t> binedges)	{ fBinsTrackPDG = binedges; }
  void			   SetDxBinning(vector<Double_t> binedges)			   { fBinsDx = binedges; }
  void			   SetDzBinning(vector<Double_t> binedges)			   { fBinsDz = binedges; }
  void			   SetDecayBinning(vector<Double_t> binedges)			{ fBinsDecay = binedges; }
  
  void         SetMCtruth(Bool_t mctruth)                     {fMCtruth=mctruth;}
  
protected:
  
  void         FillQAHistograms(AliAODTrack *coi); // Fill some QA histograms
  void         EtIsoClusOnlyEtaBand(AliAODTrack *t, Double_t &ptIso, Double_t &etaBand, Int_t index);   //EIsoCone via ClusterOnly 
  void         PtIsoTrackPhiBand(AliAODTrack *t, Double_t &ptIso, Double_t &phiBand, Int_t index);   //PIsoCone via Track UE via PhiBand TPC
  void         PtIsoTrackEtaBand(AliAODTrack *t, Double_t &ptIso, Double_t &etaBand, Int_t index);   //PIsoCone via Track UE via EtaBand TPC
  void         PtIsoTrackOrthCones(AliAODTrack *t, Double_t &ptIso, Double_t &cones, Int_t index);   //PIsoCone via Tracks UE via Orthogonal Cones in Phi
  void         PtIsoTrackFullTPC(AliAODTrack *t, Double_t &ptIso, Double_t &full, Int_t index);      //PIsoCone via Tracks UE via FullTPC - IsoCone - B2BEtaBand
  
  Bool_t       CheckBoundaries(AliAODTrack *coi);
  Double_t*    GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const;
  void         ExecOnce();
  Bool_t       Run();
  Bool_t       SelectCandidate(AliAODTrack* );
  void         AnalyzeMC();
  void         LookforParticle(Int_t, Double_t, Double_t, Double_t, Double_t);
  void         IsolationAndUEinEMCAL(AliAODTrack *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index);
  void         IsolationAndUEinTPC(AliAODTrack *coi, Double_t& isolation,Double_t& ue,Double_t eTThreshold, Int_t index);
  void         AddParticleToUEMC(Double_t& sumUE,AliAODMCParticle* mcpp,Double_t eta,Double_t phi);
  void         CalculateUEDensityMC(Double_t& sumUE);
  using        AliAnalysisTaskEmcal::FillGeneralHistograms;
  Bool_t       FillGeneralHistograms(AliAODTrack *TOI, Int_t index);
  
  AliAODEvent         *fAOD;                           //!<!
  AliVEvent           *fVevent;                        //!<! AliVEvent
  
  TClonesArray        *fNTracks;                       //
  TClonesArray        *fAODMCParticles;                //!<!
  AliAODMCHeader      *fmcHeader;                      //!<!
  TClonesArray        *fTracksAna;                     //!<! hybrid track array in
  TClonesArray        *fClusters;                      //!<! clusters container
  AliStack            *fStack;                         //!<!
  AliEMCALRecoUtils   *fEMCALRecoUtils;                //!<!  EMCAL utils for cluster rereconstruction.
  
  Float_t             fIsoConeRadius;                  // Radius for the Isolation Cont
  Int_t               fptIsoMethod;                    // Isolation definition 0=SumEt<EtThr, 1=SumEt<%Ephoton, 2=Etmax<EtThr
  Double_t            fptIsoThreshold;                 // Et isolation threshold, supposed to be % if method one is choosed
  Int_t               fIsoMethod;                      // IsoMethod =0 ->Tracks; IsoMethod =1 Clusters Only
  Double_t            fdetacut;                        // cut on deta between track and MCparticle
  Double_t            fdphicut;                        // cut on dphi between track and MCparticle

  Bool_t              fQA;                             // Flag for few further QA plots
  Int_t               fUEMethod;                       //
  Bool_t              fIsMC;                           // Flag for MC Truth Analysis
  Bool_t              fTPC4Iso;                        // 0=EMCAL_ONLY; 1=Candidate in EMCAL+ TPC for Isolation and UE
  Int_t               fNDimensions;                    //!<!number of Dimensions for the THnSPARSE Reconstruction
  Int_t               fMCDimensions;                   //!<!number of Dimensions for the THnSPARSE Truth
  Int_t               fMCQAdim;                        //!<!number of Dimensions for the THnSPARSE Mix
  Bool_t              fisLTAnalysis;                   // Flag to pass from Leading Clusters Analysis to a NC One
  Bool_t              fRejectionEventWithoutTracks;    // able/disable rejction of events without tracks
  Bool_t              fAnalysispPb;                    // able/disable the pPb analysis facilities
  Int_t               fTriggerLevel1;                  // enable to "simulate" the L1 trigger in MC: 1 = EG1 and 2 = EG2
  Int_t               fTest1;
  Int_t               fTest2;
  Bool_t              fMCtruth;
  
    // Initialization of variables for THnSparse
  
  std::vector<Double_t>    fBinsPt;
  std::vector<Double_t>    fBinsPtiso;
  std::vector<Double_t>    fBinsPtue;
  std::vector<Double_t>    fBinsEta;
  std::vector<Double_t>    fBinsPhi;
  std::vector<Double_t>    fBinsLabel;
  std::vector<Double_t>    fBinsPDG;
  std::vector<Double_t>    fBinsMomPDG;
  std::vector<Double_t>    fBinsTrackPDG;
  std::vector<Double_t>    fBinsDx;
  std::vector<Double_t>    fBinsDz;
  std::vector<Double_t>    fBinsDecay;

    //IMPLEMENT ALL THE HISTOGRAMS AND ALL THE OUTPUT OBJECTS WE WANT!!!
  TH1D        *fTrackMult;                      ///< Track Multiplicity ---QA
  TH2D        *fEtaPhiTrack;                     ///< EMCAL Track Distribution EtaPhi ---QA
  TH1D        *fPT;                             //!<! Pt distribution
  TH1D        *fP;                              //!<! P distribution
  TH1D        *fVz;                             //!<! Veretex Z distribution
  TH1D        *fEvents;                         //!<! Number of Events
  TH1D        *fPtaftTime;                      //!<! E distribution for clusters after cluster time cut
  TH1D        *fPtaftFC;                        //!<! E distribution for clusters after fiducial cut
  TH1D        *fTrackTime;                       //!<! Time distribution for clusters
  TH2D        *fPtIsoTrack;                     //!<! Isolation Pt with Tracks
  
  TH2D        *fPhiBandUETracks;                //!<! UE with Phi Band (Tracks)
  TH2D        *fEtaBandUETracks;                //!<! UE with Eta Band (Tracks)
  TH2D        *fPerpConesUETracks;              //!<! UE with Cones (Tracks ONLY)
  TH2D        *fTPCWithoutIsoConeB2BbandUE;     //!<! UE with Full TPC except IsoCone and EtaBand in Back2Back
  TH1D        *fPtIsolatedNTracks;               //!<! Isolated neutral clusters

  TH2D        *fTestIndex;                      //!<! Index and local index test
  TH2D        *fTestIndexPt;
  TH2D        *fTestLocalIndexPt;
  TH2D        *fTrackMultvsPt;
  TH1F        *fHistXsection;
  TH1F        *fHistTrials;
  TH2F        *fPtTracksVSpTTR;                 //!<!
  TH2F        *fPtTracksVSpTTR_MC;              //!<!
  TH2F        *fPhiTracksVStrackPt;             //!<!
  TH2F        *fEtaTracksVStrackPt;             //!<!
 
  THnSparse   *fOutputTHnS;                     //!<! 1st Method 4 Output
  THnSparse   *fOutMCTruth;                     //!<! 1st Method 4 MC truth Output //Isolation on pTMax
  THnSparse   *fOutTrackMC;                     //!<! 1st Method 4 MC+Truth Output via Track label

private:
  AliAnalysisTaskEMCALTrackIsolation(const AliAnalysisTaskEMCALTrackIsolation&);            // not implemented
  AliAnalysisTaskEMCALTrackIsolation&operator=(const AliAnalysisTaskEMCALTrackIsolation&); // not implemented
  
    /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALTrackIsolation, 3);    //EMCAL Neutrals base analysis task
                                                        /// \endcond
};
#endif



