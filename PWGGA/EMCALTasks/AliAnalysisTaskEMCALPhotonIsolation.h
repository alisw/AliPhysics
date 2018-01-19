#ifndef ALIANALYSISTASKEMCALPHOTONISOLATION_H
#define ALIANALYSISTASKEMCALPHOTONISOLATION_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/// \class AliAnalysisTaskEMCALPhotonIsolation
/// \ingroup EMCALTasks
/// \brief Task for isolated photons in pp, p-Pb collisions and eventually g-h correlations.
///
/// Task designed to study isolated photons using different methods to study underlying events and perform isolation.
/// It is based on the jet framework.
///
/// \author Lucile Ronflette <lucile.ronflette@cern.ch>, Subatech, Nantes
/// \author Davide Francesco Lodato <davide.francesco.lodato@cern.ch>, Utrecht University
/// \author Marco Marquard <marco.marquard@cern.ch>, University Frankfurt am Main
/// \date Jun 26, 2014

// ROOT System
class TH1D;
class TH2D;
class TH3D;
class TF1;
class THnSparse;
class TList;
class TObjArray;
class AliEMCALGeometry;
class AliESDCaloCells;
class AliESDEvent;
class AliESDtrack;
class TClonesArray;
class TList;
class TString;
class AliVCluster;
class AliVParticle;
class AliESDtrackCuts;
class AliAODEvent;
class AliAODCaloCells;
class AliVCluster;
class AliMCEvent;
class AliStack;
class TParticle;
class AliClusterContainer;
class AliParticleContainer;
class AliTrackContainer;
class AliEmcalParticle;
// AliRoot Syste
class AliEMCALTrack;
// class AliMagF;
class AliEMCALRecoUtils;
// class AliAnalysisFilter;
class AliAODTrack;
class AliAODCaloCluster;
class AliESDCaloCluster;
class AliVCaloCells;
class AliAODMCParticle;
class AliGenPythiaEventHeader;
class AliAODMCHeader;
// class AliEventPoolManager;

#include "AliAnalysisTaskEmcal.h"
#include "AliHistogramRanges.h"
#include <vector>
using std::vector;

class AliAnalysisTaskEMCALPhotonIsolation: public AliAnalysisTaskEmcal {

 public:

  AliAnalysisTaskEMCALPhotonIsolation();
  AliAnalysisTaskEMCALPhotonIsolation ( const char * name, Bool_t histo = kFALSE );
  virtual ~AliAnalysisTaskEMCALPhotonIsolation();
  
  void                         UserCreateOutputObjects();
  virtual AliHistogramRanges * GetHistogramRangesAndBinning()                             { if(!fHistoRangeContainer) fHistoRangeContainer = new AliHistogramRanges(); return fHistoRangeContainer; }
  
  void                         SetIsoConeRadius ( Float_t r )                             { fIsoConeRadius = r ; }
  void                         SetEtIsoThreshold ( Float_t r )                            { fEtIsoThreshold = r ; }
  void                         SetCTMdeltaEta ( Float_t r )                               { fdetacut = r ; }
  void                         SetCTMdeltaPhi ( Float_t r )                               { fdphicut = r ; }
  void                         SetCTMdeltaEtaIso ( Float_t r )                            { fdetacutIso = r ; }
  void                         SetCTMdeltaPhiIso ( Float_t r )                            { fdphicutIso = r ; }
  void                         SetIsoMethod ( Int_t r )                                   { fIsoMethod = r ; }
  void                         SetEtIsoMethod ( Int_t r )                                 { fEtIsoMethod = r ; }
  void                         SetUEMethod ( Int_t UE )                                   { fUEMethod = UE; }
  void                         SetOutputFormat ( Int_t iOut )                             { fWho = iOut; }
  void                         SetQA ( Bool_t QA )                                        { fQA = QA; }
  void                         SetMC ( Bool_t MC )                                        { fIsMC = MC; }
  void                         SetUSEofTPC ( Bool_t TPC )                                 { fTPC4Iso = TPC; }
  void                         SetLCAnalysis ( Bool_t LC )                                { fisLCAnalysis = LC; }
  void                         SetNLMCut ( Bool_t isNLMCut, Int_t NLMCut, Int_t NLMmin )  { fIsNLMCut = isNLMCut; fNLMCut = NLMCut; fNLMmin = NLMmin; }
  void                         SetSmearForClusters ( Int_t whichNLM )                     { fWhich = whichNLM; }
  void                         SetTMClusterRejection ( Bool_t tm )                        { fTMClusterRejected = tm; }
  void                         SetTMClusterRejectioninCone ( Bool_t tm )                  { fTMClusterInConeRejected = tm; }
  void                         SetRejectEventWithoutTracks ( Bool_t revwotr )             { fRejectionEventWithoutTracks = revwotr; }
  void                         SetAnalysispPb ( Bool_t ana )                              { fAnalysispPb = ana; }
  void                         SetTriggerLevel1 ( Int_t L )                               { fTriggerLevel1 = L; }
  void                         SetM02Smearing ( Bool_t smear )                            { fSSsmearing = smear; }
  void                         SetWidth4Smear ( Float_t width )                           { fSSsmearwidth = width; }
  void                         SetMean4Smear ( Float_t mean )                             { fSSsmear_mean = mean; }
  void                         SetExtraIsoCuts ( Bool_t bExtraIsoCuts )                   { fExtraIsoCuts = bExtraIsoCuts; }
  void			       SetPtBinning ( vector<Double_t> binedges )                 { fBinsPt = binedges; }
  void			       SetM02Binning ( vector<Double_t> binedges )                { fBinsM02 = binedges; }
  void			       SetEtisoBinning ( vector<Double_t> binedges )              { fBinsEtiso = binedges; }
  void			       SetEtueBinning ( vector<Double_t> binedges )               { fBinsEtue = binedges; }
  void			       SetEtaBinning ( vector<Double_t> binedges )                { fBinsEta = binedges; }
  void			       SetPhiBinning ( vector<Double_t> binedges )                { fBinsPhi = binedges; }
  void			       SetLabelBinning ( vector<Double_t> binedges )              { fBinsLabel = binedges; }
  void			       SetPDGBinning ( vector<Double_t> binedges )                { fBinsPDG = binedges; }
  void			       SetMomPDGBinning ( vector<Double_t> binedges )             { fBinsMomPDG = binedges; }
  void			       SetClustPDGBinning ( vector<Double_t> binedges )           { fBinsClustPDG = binedges; }
  void			       SetDxBinning ( vector<Double_t> binedges )                 { fBinsDx = binedges; }
  void			       SetDzBinning ( vector<Double_t> binedges )                 { fBinsDz = binedges; }
  void			       SetDecayBinning ( vector<Double_t> binedges )              { fBinsDecay = binedges; }
  void                         SetMCtruth ( Bool_t mctruth )                              { fMCtruth = mctruth; }
  void                         SetPeriod ( const char * period )                          { fPeriod = period; }
  void                         SetRejectPileUpEvent ( Bool_t rpue )                       { fRejectPileUpEvent = rpue; }
  void                         SetNcontributorsToPileUp ( Int_t nCtoPU )                  { fNContrToPileUp = nCtoPU; }
  void                         SetFiducialCut ( Float_t fiducial )                        { fFiducialCut = fiducial; }
  void                         SetComputeAreasPerEvent ( Bool_t eventAreas )              { fAreasPerEvent = eventAreas; }
  void                         Set2012L1Analysis ( Bool_t is2012L1 )                      { f2012EGA = is2012L1; }
  void                         SetANWithNoSameTcard ( Bool_t iNoSameTCard )               { fANnoSameTcard=iNoSameTCard; }
  
 protected:
  
  void                         FillQAHistograms      ( AliVCluster * coi, TLorentzVector vecCOI );                          // Fill some QA histograms
  void                         EtIsoCellPhiBand      ( TLorentzVector c, Double_t &etIso, Double_t &phiBand );              // EIsoCone via Cells UE via PhiBand EMCal
  void                         EtIsoCellEtaBand      ( TLorentzVector c, Double_t &etIso, Double_t &etaBand );              // EIsoCone via Cells UE via EtaBand EMCal
  void                         EtIsoClusPhiBand      ( TLorentzVector c, Double_t &etIso, Double_t &etaBand, Int_t index ); // EIsoCone via Clusters + Track UE via EtaBand EMCal
  void                         EtIsoClusEtaBand      ( TLorentzVector c, Double_t &etIso, Double_t &etaBand, Int_t index ); // EIsoCone via Clusters + Track UE via EtaBand EMCal
  void                         PtIsoTrackPhiBand     ( TLorentzVector c, Double_t &ptIso, Double_t &phiBand );              // PIsoCone via Track UE via PhiBand TPC
  void                         PtIsoTrackEtaBand     ( TLorentzVector c, Double_t &ptIso, Double_t &etaBand );              // PIsoCone via Track UE via EtaBand TPC
  void                         PtIsoTrackOrthCones   ( TLorentzVector c, Double_t &ptIso, Double_t &cones );                // PIsoCone via Tracks UE via Orthogonal Cones in Phi
  void                         PtIsoTrackFullTPC     ( TLorentzVector c, Double_t &ptIso, Double_t &full );                 // PIsoCone via Tracks UE via FullTPC - IsoCone - B2BEtaBand
  void                         ComputeConeArea       ( TLorentzVector c, Double_t &coneArea );                              // Isolation cone area depending on the cluster position
  void                         ComputeEtaBandArea    ( TLorentzVector c, Double_t &etaBandArea_InclCone );                  // Eta-band area depending on the cluster position
  void                         ApplySmearing         ( AliVCluster * coi, Double_t &m02COI );                               // Applying smearing on MC

  Bool_t                       ClustTrackMatching    ( AliVCluster * emccluster, Bool_t candidate );
  Int_t                        GetNLM                ( AliVCluster * coi, AliVCaloCells * cells );
  Int_t                        GetNLM                ( AliVCluster * coi, AliVCaloCells * cells, Int_t * absIdList, Float_t * maxEList );
  Bool_t                       AreNeighbours         ( Int_t abscell1, Int_t abscell2 ) const;
  Float_t                      RecalEnClust          ( AliVCluster * cluster, AliVCaloCells * cells );
  void                         RecalAmpCell          ( Float_t &amp, Int_t absId ) const ;
  Bool_t                       IsAbsIDsFromTCard     ( Int_t absId1, Int_t absId2, Int_t & rowDiff, Int_t & colDiff ) const ;

  Bool_t                       CheckBoundaries       ( TLorentzVector vecCOI );
  void                         FillInvMassHistograms ( Bool_t iso, Double_t m02COI, TLorentzVector c, Int_t index, Double_t isolation );
  
  Double_t*                    GenerateFixedBinArray ( Int_t n, Double_t min, Double_t max ) const;
  void                         ExecOnce();
  Bool_t                       Run();
  Bool_t                       SelectCandidate       ( AliVCluster * coi );
  void                         AnalyzeMC();
  void                         LookforParticle       ( Int_t, Double_t, Double_t, Double_t, Double_t, Double_t, Double_t );
  Bool_t                       MCSimTrigger          ( AliVEvent * eventIn, Int_t triggerLevel );
  void                         IsolationAndUEinEMCAL ( AliVCluster * coi, Double_t &isolation, Double_t &ue, Double_t eTThreshold, Int_t index );
  void                         IsolationAndUEinTPC   ( AliVCluster * coi, Double_t &isolation, Double_t &ue, Double_t eTThreshold, Int_t index );
  void                         AddParticleToUEMC     ( Double_t &sumUE, AliAODMCParticle * mcpp, Double_t eta, Double_t phi );
  void                         CalculateUEDensityMC  ( Double_t &sumUE );


  using AliAnalysisTaskEmcal::FillGeneralHistograms;
  Bool_t                       FillGeneralHistograms ( AliVCluster * COI, TLorentzVector VecCOI, Int_t index );
  
  AliAODEvent                * fAOD;                            ///<
  AliVEvent                  * fVevent;                         ///< AliVEvent
  
  TClonesArray               * fNCluster;                       ///< Neutral clusters
  TClonesArray               * fAODMCParticles;                 ///<
  AliAODMCHeader             * fmcHeader;                       ///<
  TClonesArray               * fTracksAna;                      ///< Hybrid track array in
  AliStack                   * fStack;                          ///<
  AliEMCALRecoUtils          * fEMCALRecoUtils;                 ///< EMCal utils for cluster rereconstruction.
  
  Int_t                        fWho;                            ///< Mode for the output objects ( 0 = TTree, 1 = THnSparse, 2 = TH * D/TH * F )
  Bool_t                       fSSsmearing;                     ///<
  Float_t                      fSSsmearwidth;                   ///<
  Float_t                      fSSsmear_mean;                   ///<
  Int_t                        fWhich;                          ///<
  Bool_t                       fRejectPileUpEvent;              ///<
  Int_t                        fNContrToPileUp;                 ///<
  
  /* TList       * fOutputList;                                     ///< Output list */
  /* TGeoHMatrix * fGeomMatrix[12];                                 ///< Geometry misalignment matrices for EMCal */
  
  Float_t                      fIsoConeRadius;                  ///< Radius for the isolation cone
  Int_t                        fEtIsoMethod;                    ///< Isolation definition ( 0 = SumEt<EtThr, 1 = SumEt<%Ephoton, 2 = Etmax<EtThr )
  Double_t                     fEtIsoThreshold;                 ///< Et isolation threshold, supposed to be a percentage ( < 1 ) if method one is chosen ( fEtIsoMethod = 1 )
  Double_t                     fdetacut;                        ///< Cut on deta between track and cluster
  Double_t                     fdphicut;                        ///< Cut on dphi between track and cluster
  Double_t                     fdetacutIso;                     ///< Cut on deta between track and cluster for isolation
  Double_t                     fdphicutIso;                     ///< Cut on dphi between track and cluster for isolation
  Double_t                     fM02mincut;                      ///< lambda0^2 ( sigma_long^2 ) minimum cut
  Double_t                     fM02maxcut;                      ///< lambda0^2 ( sigma_long^2 ) maximum cut
  Bool_t                       fExtraIsoCuts;                   ///< Enable/disable cuts on Ncell and DTBC for clusters in Eiso calculation
  Bool_t                       fQA;                             ///< Enable/disable a few further QA plots wrt the ones already done in the EMCalTask
  Bool_t                       fIsMC;                           ///< Enable/disable MC analysis
  Bool_t                       fTPC4Iso;                        ///< Acceptance for isolation and UE studies ( 0 = candidate in EMCal acceptance, 1 = candidate in TPC acceptance )
  Int_t                        fIsoMethod;                      ///< Isolation method ( 0 = cells, 1 = clusters + tracks, 2 = tracks only ( within EMCal or TPC acceptance ), 3 = clusters only )
  Int_t                        fUEMethod;                       ///< UE method ( within EMCal or TPC acceptance: 0 = phi-band, 1 = eta-band; only with TPC: 2 = orthogonal cones, 3 = full TPC )
  Int_t                        fNDimensions;                    ///< Number of dimensions for the THnSPARSE reconstruction
  Int_t                        fMCDimensions;                   ///< Number of dimensions for the THnSPARSE truth
  Int_t                        fMCQAdim;                        ///< Number of dimensions for the THnSPARSE mix
  Bool_t                       fisLCAnalysis;                   ///< Flag to pass from Leading Cluster ( LC ) analysis to a NC One
  Bool_t                       fIsNLMCut;                       ///< Enable/disable cut on NLM
  Int_t                        fNLMCut;                         ///< Max value for NLM cut
  Int_t                        fNLMmin;                         ///< Min value of NLM
  Bool_t                       fTMClusterRejected;              ///< Enable/disable TM cluster rejection
  Bool_t                       fTMClusterInConeRejected;        ///< Enable/disable TM cluster rejection in isolation cone
  Bool_t                       fRejectionEventWithoutTracks;    ///< Enable/disable rejction of events without tracks
  Bool_t                       fAnalysispPb;                    ///< Enable/disable the p-Pb analysis facilities
  Int_t                        fTriggerLevel1;                  ///< Choice of the L1 gamma trigger to "simulate" for the MC ( 0 = no simulation ( "MB" case ), 1 = EMCEGA1, 2 = EMCEGA2 )
  Bool_t                       fMCtruth;                        ///< Enable/disable MC truth analysis
  TString                      fPeriod;                         ///< String containing the LHC period
  Float_t                      fFiducialCut;                    ///< Variable fiducial cut from the border of the EMCal/TPC acceptance
  Bool_t                       fAreasPerEvent;                  ///< Enable/disable the event-by-event cone area computation
  
  // Initialization for TTree variables
  Double_t                     fEClustersT;                     ///< E for all clusters
  Double_t                     fPtClustersT;                    ///< Pt for all clusters
  Double_t                     fEtClustersT;                    ///< Et for all clusters
  Double_t                     fEtaClustersT;                   ///< Eta for all clusters
  Double_t                     fPhiClustersT;                   ///< Phi for all clusters
  Double_t                     fM02ClustersT;                   ///< lambda0^2 (sigma_long^2) for all clusters
  Int_t                        fevents;                         ///< Number of events
  Int_t                        fNClustersT;                     ///< Clusters multiplicity
  Double_t                     flambda0T;                       ///< M02 for considered clusters (leading one or all depending on flag)
  Double_t                     fM02isoT;                        ///< M02 for isolated clusters
  Double_t                     fM02noisoT;                      ///< M02 for non isolated clusters
  Double_t                     fPtnoisoT;                       ///< Pt for non isolated clusters
  Double_t                     fEtT;                            ///< Et for considered clusters (leading one or all depending on flag)
  Double_t                     fPtT;                            ///< Pt for considered clusters (leading one or all depending on flag)
  Double_t                     fPtisoT;                         ///< Pt for all isolated neutral clusters
  Double_t                     fEtisolatedT;                    ///< Et for isolated clusters
  Double_t                     fPtisolatedT;                    ///< Pt for isolated clusters
  Double_t                     fetaT;                           ///< Eta for considered clusters
  Double_t                     fphiT;                           ///< Phi for considered clusters
  Double_t                     fsumEtisoconeT;                  ///< Sum Et in cone
  Double_t                     fsumEtUE;                        ///< Sum UE
  Bool_t                       fANnoSameTcard;                  ///<
  
  // Initialization of variables for THnSparse
  std::vector<Double_t>        fBinsPt;                         ///<
  std::vector<Double_t>        fBinsM02;                        ///<
  std::vector<Double_t>        fBinsEtiso;                      ///<
  std::vector<Double_t>        fBinsEtue;                       ///<
  std::vector<Double_t>        fBinsEta;                        ///<
  std::vector<Double_t>        fBinsPhi;                        ///<
  std::vector<Double_t>        fBinsLabel;                      ///<
  std::vector<Double_t>        fBinsPDG;                        ///<
  std::vector<Double_t>        fBinsMomPDG;                     ///<
  std::vector<Double_t>        fBinsClustPDG;                   ///<
  std::vector<Double_t>        fBinsDx;                         ///<
  std::vector<Double_t>        fBinsDz;                         ///<
  std::vector<Double_t>        fBinsDecay;                      ///<
  
  //IMPLEMENT ALL THE HISTOGRAMS AND ALL THE OUTPUT OBJECTS WE WANT!!!

  TH1D                       * fTrackMult;                      ///<  Track Multiplicity ---QA
  TH2D                       * fPtvsSumUE_MC;                   //!<!
  TH2D                       * fEtaPhiClus;                     ///<  EMCal Cluster Distribution EtaPhi ---QA
  TH2D                       * fClusEvsClusT;                   //!<! Cluster Energy vs Cluster Time ---QA
  TH1D                       * fPT;                             //!<! Pt distribution
  TH1D                       * fE;                              //!<! E distribution
  TH2D                       * fNLM;                            //!<! NLM distribution
  TH2D                       * fNLM2_NC_Acc;                    //!<! NLM (1,2) distribution for Neutral Clusters in Acceptance
  TH2D                       * fNLM2_NC_Acc_noTcard;            //!<! NLM (1,2) distribution for Neutral Clusters in Acceptance w/o NLM=2 in SameT
  TH1D                       * fVz;                             //!<! Vertex Z distribution
  TH1D                       * fEvents;                         //!<! Number of Events
  TH1D                       * fPtaftTime;                      //!<! E distribution for clusters after Cluster Time cut
  TH1D                       * fPtaftCell;                      //!<! Pt distribution for clusters after NCells cut
  TH1D                       * fPtaftNLM;                       //!<! Pt distribution for clusters after NLM cut
  TH3F                       * fClusEtVsEtaPhiMatched;          //!<! Track-matched cluster eta vs. phi vs. E_T
  TH3F                       * fClusEtVsEtaPhiUnmatched;        //!<! Not track-matched cluster eta vs. phi vs. E_T
  TH1D                       * fPtaftTM;                        //!<! E distribution for neutral clusters
  TH1D                       * fPtaftDTBC;                      //!<! E distribution for NC after DistanceToBadChannel cut
  TH1D                       * fPtaftFC;                        //!<! E distribution for clusters after Fiducial cut
  TH1D                       * fPtaftM02C;                      //!<! E distribution for clusters after Shower Shape cut
  TH1D                       * fClusTime;                       //!<! Time distribution for clusters
  TH2D                       * fM02;                            //!<! Squared_Lambda0 (Squared_sigma_long) distribution
  TH3F                       * fEtaPhiClusVsM02;                //!<! Cluster eta vs. phi vs. sigma_long squared (cluster energy from 14 to 16 GeV)
  TH3F                       * fEtaPhiClusVsEtIsoClus;          //!<! Cluster eta vs. phi vs. neutral contribution to the energy in isolation cone (cluster energy from 14 to 16 GeV)
  TH3F                       * fEtaPhiClusVsPtIsoTrack;         //!<! Cluster eta vs. phi vs. charged contribution to the energy in isolation cone (cluster energy from 14 to 16 GeV)
  TH3F                       * fEtaPhiClusVsPtUETrackCside;     //!<! Cluster eta vs. phi vs. charged contribution to the energy in C-side UE (cluster energy from 14 to 16 GeV)
  TH3F                       * fEtaPhiClusVsPtUETrackAside;     //!<! Cluster eta vs. phi vs. charged contribution to the energy in A-side UE (cluster energy from 14 to 16 GeV)
  TH1D                       * fDeltaETAClusTrack;              //!<! dEta Cluster-Track
  TH1D                       * fDeltaPHIClusTrack;              //!<! dPhi Cluster-Track
  TH1D                       * fDeltaETAClusTrackMatch;         //!<! dEta Cluster-Track matched
  TH1D                       * fDeltaPHIClusTrackMatch;         //!<! dPhi Cluster-Track matched
  TH1D                       * fEtIsoCells;                     //!<! Isolation Energy with EMCal Cells
  TH2D                       * fEtIsoClust;                     //!<! Isolation Energy with EMCal Clusters
  TH2D                       * fPtIsoTrack;                     //!<! Isolation Pt with Tracks
  TH1D                       * fPtEtIsoTC;                      //!<! Isolation with Pt from Tracks and Et from NON-Matched Clusters
  TH2D                       * fPhiBandUEClust;                 //!<! UE with Phi Band (Clusters)
  TH2D                       * fEtaBandUEClust;                 //!<! UE with Eta Band (Clusters)
  TH2D                       * fPhiBandUECells;                 //!<! UE with Phi Band (Cells)
  TH2D                       * fEtaBandUECells;                 //!<! UE with Eta Band (Cells)
  TH2D                       * fPhiBandUETracks;                //!<! UE with Phi Band (Tracks)
  TH2D                       * fEtaBandUETracks;                //!<! UE with Eta Band (Tracks)
  TH2D                       * fPerpConesUETracks;              //!<! UE with Cones (Tracks ONLY)
  TH2D                       * fTPCWithoutIsoConeB2BbandUE;     //!<! UE with Full TPC except IsoCone and EtaBand in Back2Back
  TH1D                       * fNTotClus10GeV;                  //!<! Number of TOTAL clusters with Energy bigger than 10 GeV
  TH1D                       * fEtIsolatedCells;                //!<! Isolated photons, isolation with cells
  TH1D                       * fEtIsolatedClust;                //!<! Isolated photons, isolation with clusters
  TH1D                       * fPtIsolatedNClust;               //!<! Isolated neutral clusters
  TH1D                       * fPtIsolatedNTracks;              //!<! Isolated neutral clusters with tracks
  TH1D                       * fEtIsolatedTracks;               //!<! Isolated photons, isolation with tracks
  TH2D                       * fPtvsM02iso;                     //!<! Isolated clusters, Pt distribution vs M02
  TH2D                       * fPtvsM02noiso;                   //!<! Non isolated clusters, pt distribution vs M02
  TH2D                       * fTestIndex;                      //!<! Index and local index test
  TH2D                       * fTestIndexE;                     //!<! Index vs cluster energy test
  TH2D                       * fTestLocalIndexE;                //!<! Local index vs cluster energy test
  TH3F                       * fTestEnergyCone;                 //!<! Energy cone clusters vs tracks test
  TH2F                       * fEtaBandVsConeArea;              //!<! Eta-band vs. cone area distribution (depending on the cluster position)
  TH3F                       * fPtVsConeVsEtaBand;              //!<! Energy cone clusters vs tracks test (not normalised)
  TH3F                       * fPtVsNormConeVsNormEtaBand;      //!<! Energy cone clusters vs tracks test (area normalised)
  TH3F                       * fPtvsM02vsSumUE_Norm;
  TH2D                       * fTestEtaPhiCone;                 //!<! Eta vs phi test for clusters in cone
  TH3D                       * fInvMassM02iso;                  //!<!
  TH3D                       * fInvMassM02noiso;                //!<!
  TH3F                       * fPtvsM02vsSum;                   //!<!
  TH3F                       * fPtvsM02vsSumUE;                 //!<!
  TH3D                       * fTrackMultvsSumChargedvsUE;      //!<!
  TH2D                       * fTrackMultvsPt;                  //!<!
  TH3D                       * fTracksConeEtaPt;                //!<!
  TH3D                       * fTracksConeEtaM02;               //!<!
  TH1F                       * fHistXsection;                   //!<!
  TH1F                       * fHistTrials;                     //!<!
  TH2F                       * fPtTracksVSpTNC;                 //!<!
  TH2F                       * fCTdistVSpTNC;                   //!<!
  TH2F                       * fPtTracksVSpTNC_MC;              //!<!
  TH3F                       * fpi0VSclusterVSIsolation;        //!<!
  TH3F                       * fpi0VSclusterVSM02;              //!<!
  TH3F                       * fpi0VSM02VSIsolation;            //!<!
  TH3F                       * fEtVSM02VSEisoclust;             //!<!
  TH3F                       * fEtVSM02VSPisotrack;             //!<!
  TH2F                       * fPhiTracksVSclustPt;             //!<!
  TH2F                       * fEtaTracksVSclustPt;             //!<!
  TH2F                       * fTracksPhiVsPt;                  //!<!
  TH2F                       * fTracksEtaVsPt;                  //!<!
  TH2F                       * fTrackResolutionPtMC;            //!<!
  TH1D                       * fVzBeforecut;                    //!<!
  
  THnSparse                  * fOutputTHnS;                     ///< 1st Method 4 Output
  THnSparse                  * fOutMCTruth;                     ///< 1st Method 4 MC truth Output // Isolation on pTMax
  THnSparse                  * fOutClustMC;                     ///< 1st Method 4 MC+Truth Output via Clusterlabel
  
  TTree                      * fOutputQATree;                   ///< 2nd method 4 QA Output
  TTree                      * fOutputTree;                     ///< 2nd Method 4 Output
  
  TH3D                       * fphietaPhotons;                  //!<!
  TH3D                       * fphietaOthers;                   //!<!
  TH3D                       * fphietaOthersBis;                //!<!
  TH2F                       * fSPDclustVsSPDtracklets;         //!<!
  TH1F                       * fnPUevents;                      //!<!
  Bool_t                       f2012EGA;                        // Analyze only Events with EGA recalc patches above threshold
  Int_t                        fAbsIDNLM[2];                    //!<!
  // TH1                        * fPDGM02;                         //!<! check for zeroM02 clusters
  // TH2                        * fEtrueEclustM02;                 //!<! check for zeroM02 clusters
  // TH2                        * fDphiDetaM02;                    //!<! check for zeroM02 clusters
  // TH1D                       * fMomPDGM02;                      //!<!
  // TH2D                       * fTvsE_MismatchEM02;              //!<!
  
  /* AliParticleContainer       * fTracksCont;                     ///< Tracks */
  /* AliParticleContainer       * fclusters;                       ///< Container for Particle container 4 clusters */
  
 private:

  AliHistogramRanges         * fHistoRangeContainer;            ///< Histogram ranges container

  AliAnalysisTaskEMCALPhotonIsolation ( const AliAnalysisTaskEMCALPhotonIsolation & );             // Not implemented
  AliAnalysisTaskEMCALPhotonIsolation&operator = ( const AliAnalysisTaskEMCALPhotonIsolation & ); // Not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEMCALPhotonIsolation, 22);            // EMCal neutrals base analysis task
  /// \endcond
};
#endif
