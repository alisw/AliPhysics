#ifndef ALIANALYSISTASKREVEA_H
#define ALIANALYSISTASKREVEA_H


class TH1I;
class TH1F;
class TF1;
class TH2F;
class TH2D;
class TH1D;
class TH3D;
class TRandom3;
class TLorentzVector;
class TArrayD;
class TArrayF;
class TArrayL;
class THnSparse;
class TProfile;
class TList;
class TClonesArray;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVParticle;
class AliLog;
class AliAnalysisUtils;
class AliJetContainer;
class AliParticleContainer;
class AliMultSelection;


#include <vector>
using std::vector;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"


// ANALYSIS OF H+JET IN PP 
// Author Filip Krizek   (25 JAN 2022)

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisTaskRevEA : public AliAnalysisTaskEmcalJet {
   public:

   enum Analysis {
      kDetLevel=0,
      kPartLevel=1,
      fkTTbins = 20 //maximum number of TT bins
   };

   enum {fkV0Mnorm, fkCE}; //detector level   norm : divided by mean V0M    norm2: average  A/mean A  +C/ mean C


   enum {kNormal=0, kMC=1}; //type of analysis
   enum {kMB=0, kHM=1, kTG};  //triggers   MB, HM, GA


   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskRevEA();
   AliAnalysisTaskRevEA(const char *name);
   virtual  ~AliAnalysisTaskRevEA();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   static AliAnalysisTaskRevEA* AddTaskRevEA(
       Int_t       mode               = AliAnalysisTaskRevEA::kNormal, // analysis mode   normal=0, mc=1 or embedded=2
       const char* jetarrayname       = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets (or combined event jets when emb)
       const char* jetarraynamePartMC = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
       const char* trackarrayname     = "tracks",            //name of track TClonesArray for detector level tracks  (or tracks in the combined event when embedding)
       const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
       const char* ktjetarrayname     = "", //name of rho task real data
       const char* ktjetarraynamePartMC = "", //name of rho task mc data
       Double_t    jetRadius          = 0.4,  //radius of analyzed jets
       UInt_t      trigger            = AliVEvent::kAny,  //trigger
       Double_t    trackEtaWindow     = 0.9,   //pseudorapidity range for tracks
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       Double_t    acut               = 0.,   //cut on relative jet area
       const char* suffix = ""                              //SUBWAGON has to be the last parameter
  );


  // ######### SETTERS/GETTERS

  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;}

  void        SetAcceptanceWindows(Double_t trackEta){ fTrackEtaWindow  = trackEta; }


  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void        SetTrackContainerName(const char* name){ fMyTrackContainerName = name;}
  void        SetMCParticleContainerName(const char* name){ fMyParticleContainerName = name;}
  void        SetMCDetLevelContainerName(const char* name){ fMyDetLevelContainerName = name;}

  void        SetJetContainerName(const char* name){ fMyJetContainerName = name;}
  void        SetMCPartJetContainerName(const char* name){ fMyJetParticleContainerName = name;}

  void        SetKTJetContainerName(const char* name){ fMyKTJetContainerName = name;}
  void        SetKTMCPartJetContainerName(const char* name){ fMyKTJetParticleContainerName = name;}

  void        SetHadronTT(Int_t tl,Int_t th){fHadronTTLowPt[fnHadronTTBins] = tl; fHadronTTHighPt[fnHadronTTBins] = th; fnHadronTTBins++; }

  void        SetMode(Int_t mode){ fMode = mode;}       //Analysi mode normal=0 or embedded event =1

  void        SetPhiCut(Double_t pcut){ fPhiCut = TMath::Pi()-pcut; }


  void        SetJetRadius(Double_t jr) { fJetR = jr;}

  void        SetJetAcut(Double_t ac){ fJetAcut = ac;}



  Bool_t      PassedMinBiasTrigger();
  Bool_t      PassedHighMultTrigger();



  Bool_t      RetrieveEventObjects();
  Bool_t      Run();
  Bool_t      FillHistograms();

  void AnalyzeParticleLevel();
  void InitEventProperties();
  void FindDetectorLevelTT();
  void FindParticleLevelTT();
  void FillResponseMatrix();
  void FillResponseMatrix2D();
  void GeneralTrackProperties();
  void AnalyzeRawData();

  void SetMaxFacPtHard(Float_t maxfacpthard){ fMaxFacPtHard = maxfacpthard;} 

 private:

  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isGen=0);
  Bool_t      IsEventInAcceptance(AliVEvent* event);
  Bool_t      IsOutlier(); 

  // ######### STANDARD FUNCTIONS
  void        ExecOnceLocal();

  Double_t    GetMyRho(AliJetContainer* ktjets);


  Double_t    GetDeltaPt(Double_t phiTT, Double_t etaTT, Double_t phiLJ, Double_t etaLJ, Double_t phiSJ, Double_t etaSJ, Double_t rho, Int_t level);

  // ########## USAGE TRIGGERS
  Bool_t      fUseDefaultVertexCut;                   // trigger if automatic vertex cut from helper class should be done
  Bool_t      fUsePileUpCut;                          // trigger if pileup cut should be done


  // ########## SOURCE INFORMATION
  TString     fMyTrackContainerName;                  // name of detector level track container or  tracks in combined embedded event
  TString     fMyParticleContainerName;               // name of particle level MC particle container
  TString     fMyDetLevelContainerName;               // name of detector level MC track container name used while embeding
  TString     fMyJetContainerName;                    // name of detector level jet container
  TString     fMyJetParticleContainerName;            // name of particle level MC jet container
  TString     fMyKTJetContainerName;                  // name of KT detector level jet container
  TString     fMyKTJetParticleContainerName;          // name of KT particle level MC jet container


  AliTrackContainer    *fTrkContainerDetLevel;        //! detector level track container   (or tracks in combined embedded events)
  AliParticleContainer *fParticleContainerPartLevel;  //! particle level container with pythia particles
  AliJetContainer      *fJetContainerDetLevel;        //! detector level jet container   (or jets in combined events after embedding)
  AliJetContainer      *fJetContainerPartLevel;       //! particle level jet container

  AliJetContainer      *fKTJetContainerDetLevel;        //! KT detector level jet container   (or jets in combined events after embedding)
  AliJetContainer      *fKTJetContainerPartLevel;       //! KT particle level jet container

  // ########## CENTRALITY
  AliMultSelection*  fMultSelection;                  //! object which handels centrality

  Bool_t  fIsMinBiasTrig;                             //! triggered by Min Bias Trig
  Bool_t  fIsHighMultTrig;                            //! triggered by high multiplicity trigger

  Float_t fCentralityV0M;                             //! Centrality from V0M

  Double_t  fMultV0Mnorm;                              //!  mult. V0M/mean

  // ########## CUTS
  Double_t    fTrackEtaWindow;                        // +- window in eta for tracks
  Double_t    fMinTrackPt;                            // Min track pt to be accepted

  // ########## GENERAL ////VARS
  AliAnalysisUtils*   fHelperClass;                   //! Vertex selection helper
  Bool_t              fInitializedLocal;              //! trigger if tracks/jets are loaded  initiates calling   ExecOnce


   TH1D     *fHistEvtSelection;                       //! gc event statistics

   TH1D     *fhVertexZall;                            //! gc vertexZ inclusive before cut
   TH1D     *fhVertexZ;                               //! gc vertexZ inclusive after cut

   TH2D     *fhTrackPhiIncl[kTG];                        //!  phi inclusive ch hadron
   TH2D     *fhTrackEtaIncl[kTG];                        //!  eta inclusive minimum bias
   TH2D     *fhJetPhiIncl[kTG];                          //!  phi inclusive ch jet
   TH2D     *fhJetEtaIncl[kTG];                          //!  eta inclusive

   TH1D     *fhRho[kTG];                          //! minimum bias rho inclusive
   TH1D     *fhRhoTTH[kTG][fkTTbins];             //! in events MB with hadron TT

   TH1D     *fhRhoMBpart;                         //! minimum bias rho inclusive particle level
   TH1D     *fhRhoTTHinMBpart[fkTTbins];          //! in events MB with hadron TT particle level

   TH2D* fhCentrality[kTG];                      //! estimated centrality based on   V0M

   TH1D* fhSignal[kTG][fkCE];                          //! centrality estimators:   V0Mnorm  in MB
   TH1D* fhSignalTTH[kTG][fkCE][fkTTbins];             //! distributions of centrality estimators biased with hadron TT in min bias

   TH1D* fhSignal_PartLevel[fkCE];                 //! particle level centrality estimators: V0Mnorm  in MB
   TH1D* fhSignalTTH_PartLevel[fkCE][fkTTbins];  //! particle level distributions of centrality estimators biased with hadron TT in min bias


   TH1D* fhMultTTH[kTG][fkTTbins];                       //! multiplicity of hadron TT


   //hadron TT
   TH2D* fhTTH_V0Mnorm[kTG][fkTTbins];                   //! counter of semi-inclusive hadron TT versus V0M/mean

   TH2D* fhTTH_V0Mnorm_PartLevel[fkTTbins];           //! counter of semi-inclusive hadron TT   V0M/mean particle level with V0 coincidence
   TH2D* fhTT_Corresp[fkTTbins];                       //! counter of "corresponding TT" // KA


   //recoil jet yields with hadron TT
   TH2D* fhRecoilJetPtTTH_V0Mnorm[kTG][fkTTbins];           //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M/mean
   TH2D* fhRecoilJetPtTTH_V0Mnorm_PartLevel[fkTTbins];      //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M/mean
   TH2D* fhRecoilJetPtZero_TTH_V0Mnorm_PartLevel[fkTTbins]; //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M/mean. Jet pT is NOT CORRECTED on RhokT (added by KA)

   TH3D* fhRecoilJetPhiTTH_V0Mnorm[kTG][fkTTbins];                   //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)
   TH3D* fhRecoilJetPhiTTH_V0Mnorm_PartLevel[fkTTbins];              //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  minimum bias particle level
   TH3D* fhRecoilJetPtZero_DeltaPhi_TTH_V0Mnorm_PartLevel[fkTTbins]; //! recoil jet  (V0M/mean , recoil jet pT,  delta phi). Jet pT is NOT CORRECTED on RhokT (added by KA)


   TH2D* fhDeltaPtTTH_RC_V0Mnorm[kTG][fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT  versus V0M/mean V0M

   //_______________________________________________________________________________________________________________________________________________
   // 1D Unfolding

   // Inclusive jets
   TH1D* fhJetPtPartLevelCorr;                      //! response matrix normalization spectrum, jet pT corrected on rho
   TH1D* fhJetPtPartLevelZero;                      //! response matrix normalization spectrum, jet pT is not corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorr;       //! response matrix, jet pT corrected on rho
   TH2D* fhJetPtZeroPartLevel_Vs_JetPtDetLevelCorr; //! response matrix, part level jet pT is not corrected on rho, detectot level jet pT is corrected on rho (added by KA)
   TH2D* fhJetPtZeroPartLevelVsJetPtZeroDetLevel;   //! response matrix, jet pT not corrected on rho

   // Recoil jets with TT on detector level
   TH1D* fhRecoilJetPtPartLevelCorr[fkTTbins];                        //! response matrix normalization spectrum, jet pT corrected on rho, built from recoil jets //FF
   TH1D* fhRecoilJetPtZeroPartLevel[fkTTbins];                        //! response matrix normalization spectrum, jet pT not corrected on rho, built from recoil jets
   TH2D* fhRecoilJetPtPartLevelVsJetPtDetLevelCorr[fkTTbins];         //! response matrix, jet pT corrected on rho, built from recoil jets //FF
   TH2D* fhRecoilJetPtZeroPartLevelVsJetPtDetLevelCorr[fkTTbins];     //! response matrix, part level jet pT not corrected on rho, built from recoil jets //FF
   TH2D* fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevelCorr[fkTTbins]; //! response matrix, jet pT not corrected on rho, built from recoil jets //FF

   // Recoil jets with corresponding TT
   TH1D* fhRecoilJetPtPartLevel_CorrespTT[fkTTbins];                        //! response matrix normalization spectrum, jet pT corrected on rho, built from recoil jets on particle level //KA
   TH1D* fhRecoilJetPtZeroPartLevel_CorrespTT[fkTTbins];                    //! response matrix normalization spectrum, jet pT not corrected on rho, built from recoil jets on particle level
   TH2D* fhRecoilJetPtPartLevelVsJetPtDetLevel_CorrespTT[fkTTbins];         //! response matrix jet pT corrected on rho built from jets recoil from corresponding TT // KA
   TH2D* fhRecoilJetPtZeroPartLevelVsJetPtDetLevel_CorrespTT[fkTTbins];     //! response matrix jet pT corrected on rho built from jets recoil from corresponding TT // KA
   TH2D* fhRecoilJetPtZeroPartLevelVsJetPtZeroDetLevel_CorrespTT[fkTTbins]; //! response matrix jet pT corrected on rho built from jets recoil from corresponding TT // KA

   //_________________________________________________________________________________________________________________________________________________________________________________________
   //2D unfolding

   // Inclusive jets
   TH2D *fhPhi_JetPtPartLevel_InclusiveJets;                                   //! histogram with missed events -> no matching with DetLevel jets. Inclusive Jets
   TH2D *fhPhi_JetPtZeroPartLevel_InclusiveJets;                               //! histogram with missed events -> no matching with DetLevel jets. Inclusive Jets
   THnSparse *fhPhi_JetPtDetLevel_Vs_Phi_JetPtPartLevel_InclusiveJets;         //! 4D-response matrix delta phi vs jet pT.  Inclusive Jets
   THnSparse *fhPhi_JetPtDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets;     //! 4D-response matrix delta phi vs jet pT.  Inclusive Jets, part level jet pT not corrected on rho
   THnSparse *fhPhi_JetPtZeroDetLevel_Vs_Phi_JetPtZeroPartLevel_InclusiveJets; //! 4D-response matrix delta phi vs jet pT.  Inclusive Jets, jet pT not corrected on rho

   // Recoil jets with TT on detector level
   TH2D *fhDeltaPhi_JetPtPartLevel[fkTTbins];                                        //! histogram with missed events -> no matching with DetLevel jets. TT events
   TH2D *fhDeltaPhi_JetPtZero_PartLevel[fkTTbins];                                   //! histogram with missed events -> no matching with DetLevel jets. TT events
   THnSparse *fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel[fkTTbins];         //! 4D-response matrix delta phi vs jet pT. TT events
   THnSparse *fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[fkTTbins];     //! 4D-response matrix delta phi vs jet pT. TT events, part level jet pT not corrected on rho
   THnSparse *fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel[fkTTbins]; //! 4D-response matrix delta phi vs jet pT. TT events, jet pT not corrected on rho

   // Recoil jets with corresponding TT
   TH2D *fhDeltaPhi_JetPtPartLevel_CorrespTT[fkTTbins];                                        //! histogram with missed events -> no matching with DetLevel jets. TT events
   TH2D *fhDeltaPhi_JetPtZeroPartLevel_CorrespTT[fkTTbins];                                    //! histogram with missed events -> no matching with DetLevel jets. TT events
   THnSparse *fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtPartLevel_CorrespTT[fkTTbins];         //! 4D-response matrix delta phi vs jet pT. TT events
   THnSparse *fhDeltaPhi_JetPtDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[fkTTbins];     //! 4D-response matrix delta phi vs jet pT. TT events, part level jet pT not corrected on rho
   THnSparse *fhDeltaPhi_JetPtZeroDetLevel_Vs_DeltaPhi_JetPtZeroPartLevel_CorrespTT[fkTTbins]; //! 4D-response matrix delta phi vs jet pT. TT events, jet pT not corrected on rho

   // TH1D *fhNumberOf_ChosenTT_PartLevel[fkTTbins];                             //! number of particle level TT (I suggest to delete it, since it is not needed) added by KA

   Double_t fArray_for_filling[4]; // Array for filling 4D Sparse object
   ///////////////

   TH2D* fhJetPtResolutionVsPtPartLevel;                //! resolution of jet pT

   THnSparse* fhNumberOfHighPtJetsRecoil[kTG][fkTTbins]; //! number of jets with pT larger than X in recoil region of physical TT
   TH1D* fhRecoilJetPtEvtByEvent[fkTTbins];              //! event by event pt spectrum of jets with physical TT

   THnSparse* fhNumberOfHighPtJetsRecoilPartLevel[fkTTbins]; //! number of jets with pT larger than X in recoil region of physical TT
   TH1D* fhRecoilJetPtEvtByEventPartLevel[fkTTbins];         //! event by event pt spectrum of jets with physical TT

   Double_t fZVertexCut;                                 // vertex cut in z
                                                         
   Int_t    fnHadronTTBins;                              // number of TT bins charged hadron
   Int_t    fHadronTTLowPt[fkTTbins];                    // low pt TT range charged hadron
   Int_t    fHadronTTHighPt[fkTTbins];                   // high pt TT range charged hadron


   vector<Int_t> fHadronTT_Labels[fkTTbins];            //! array which stores the label of hadron TT in given event 
   vector<Int_t> fHadronTT_Labels_PartLevel[fkTTbins];  //! array which stores the label of hadron TT in given event 

   Int_t    fMode;                                      // Analysis mode   0=real data, 1=mc pythia, 2=embedded
   
   Int_t fIndexTTH[fkTTbins];                      //! index of the chosen hadron trigger
   Int_t fIndexTTH_PartLevel[fkTTbins];             //! index of the chosen hadron trigger particle level

   vector<TLorentzVector> fTTH[fkTTbins];          //!  hadron trigger candidate

   vector<TLorentzVector> fTTH_PartLevel[fkTTbins];  //!  hadron trigger candidate particle level

   Bool_t fFillSigTT;                            //! flag which labels whether the event should be filled for ref or sig class

   Double_t fPhiCut;                             // phi angle cut on the recoil jet  pi-fPhiCut
   TRandom3* fRandom;                            //! Radom

   Double_t frhovec[999];                        //! auxiliary array to store pT/A of kT jets for kt

   Double_t fJetR;                                 // jet radius
   Double_t fJetAcut;                              // jet area cut


   Double_t fRho;                                  //! underlying event density real events
   Double_t fRhoMC;                                //! underlying event density detector level events events

   Bool_t fTrigflag[2];                            //! trigger flags
   Int_t fdeltapT[fkTTbins];                       //! delta pT detector level

   Float_t  fMaxFacPtHard;   // Cut on  pthat events. How many times can be jet pT larger than pthat //FK

   TH3D* fhNjetReMx_V0MnormDetLev_15GeV;  //! FILIP response matrix for distribution counting jets pT gt 15 GeV jets in events  
   TH2D* fhNjetNorm_V0MnormDetLev_15GeV;  //! FILIP normalization for ReMx for distribution counting jets pT gt 15 GeV jets in events  

   TH3D* fhNjetReMx_V0MnormDetLev_20GeV;  //! FILIP response matrix for distribution counting jets pT gt 20 GeV jets in events  
   TH2D* fhNjetNorm_V0MnormDetLev_20GeV;  //! FILIP normalization for ReMx for distribution counting jets pT gt 20 GeV jets in events 

   AliAnalysisTaskRevEA(const AliAnalysisTaskRevEA&);
   AliAnalysisTaskRevEA& operator=(const AliAnalysisTaskRevEA&);

   ClassDef(AliAnalysisTaskRevEA, 2); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
}
}
#endif
