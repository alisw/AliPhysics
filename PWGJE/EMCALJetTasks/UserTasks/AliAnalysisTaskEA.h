#ifndef ALIANALYSISTASKEA_H
#define ALIANALYSISTASKEA_H


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
class TRandom3;
class AliJetContainer;
class AliParticleContainer;
class AliClusterContainer;
class AliMultSelection;

namespace PWGJE {
  namespace EMCALJetTasks {
    class AliAnalysisEmcalJetHelperEA;
  }
}

#include <vector>
using std::vector;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"


// ANALYSIS OF EVENT ACTIVITY WITH HIGH PT HADRON BIAS
// Author Filip Krizek   (8 AUG 2019)

namespace PWGJE {
 
namespace EMCALJetTasks {

class AliAnalysisTaskEA : public AliAnalysisTaskEmcalJet {
   public:

   enum Analysis {
      kDetLevel=0, 
      kPartLevel=1,
      kEmbLevel=2,
      fkVtx=3,
      fkTTbins = 20 //maximum number of TT bins
   };

   enum {fkV0A, fkV0C, fkV0M, fkV0Mnorm1, fkCE}; //detector level   norm1 : divided by mean V0M    norm2: average  A/mean A  +C/ mean C


   enum {kNormal=0, kMC=1, kEmbedding=2}; //type of analysis    
   enum {kMB=0, kHM=1, kGA=2, kTG};  //triggers   MB, HM, GA
   enum {krhokt =0, krhocms=1, kRho};  // rho kt,  rho cms


   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskEA();
   AliAnalysisTaskEA(const char *name);
   virtual  ~AliAnalysisTaskEA();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   static AliAnalysisTaskEA* AddTaskEA(
       Int_t       mode               = AliAnalysisTaskEA::kNormal, // analysis mode   normal=0, mc=1 or embedded=2 
       const char* jetarrayname       = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets (or combined event jets when emb)
       const char* jetarraynamePartMC = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
       const char* jetarraynameDetMC  = "",                  //name of jet TClones array for detector level MC particle level jets for embedding
       const char* trackarrayname     = "tracks",            //name of track TClonesArray for detector level tracks  (or tracks in the combined event when embedding)  
       const char* mcpariclearraynamePartMC = "mcparticles", //name of mcparticle TClonesArray array for MC particle level jets
       const char* tracknameDetMC  = "",            //name of track TClonesArray array for detector level MC particle level jets  when doing embedding
       const char* clusterarrayname   = "caloClusters", //name of EMCAL cluster TClonesArray array  for detector level
       const char* ktjetarrayname     = "", //name of rho task real data 
       const char* ktjetarraynamePartMC = "", //name of rho task mc data
       const char* ktjetarraynameDetMC  = "", //name of rho task pythia detector level data
       Double_t    jetRadius          = 0.4,  //radius of analyzed jets
       UInt_t      trigger            = AliVEvent::kAny,  //trigger
       Double_t    trackEtaWindow     = 0.9,   //pseudorapidity range for tracks
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       Double_t    acut               = 0.,   //cut on relative jet area
       Double_t    emcaltofcut        = 30e-9,   //cut on relative jet area
       const char* suffix = ""                              //SUBWAGON has to be the last parameter
  );

   
  // ######### SETTERS/GETTERS

  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}  
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;} 

  void        SetAcceptanceWindows(Double_t trackEta){ 
                 fTrackEtaWindow  = trackEta; 
              } 


  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void        SetTrackContainerName(const char* name){ fMyTrackContainerName = name;}
  void        SetMCParticleContainerName(const char* name){ fMyParticleContainerName = name;}
  void        SetMCDetLevelContainerName(const char* name){ fMyDetLevelContainerName = name;}
  void        SetClusterContainerName(const char* name){ fMyClusterContainerName = name;}

  void        SetJetContainerName(const char* name){ fMyJetContainerName = name;} 
  void        SetMCPartJetContainerName(const char* name){ fMyJetParticleContainerName = name;}
  void        SetMCDetJetContainerName(const char* name){ fMyJetDetLevelContainerName = name;}

  void        SetKTJetContainerName(const char* name){ fMyKTJetContainerName = name;} 
  void        SetKTMCPartJetContainerName(const char* name){ fMyKTJetParticleContainerName = name;}
  void        SetKTMCDetJetContainerName(const char* name){ fMyKTJetDetLevelContainerName = name;}

  void        SetHadronTT(Int_t tl,Int_t th){ fHadronTTLowPt[fnHadronTTBins] = tl; fHadronTTHighPt[fnHadronTTBins] = th; fnHadronTTBins++; }
  void        SetJetChTT(Int_t tl,Int_t th){  fJetChTTLowPt[fnJetChTTBins] = tl;   fJetChTTHighPt[fnJetChTTBins] = th;   fnJetChTTBins++;  }
  void        SetClusterTT(Int_t tl,Int_t th){  fClusterTTLowPt[fnClusterTTBins] = tl;   fClusterTTHighPt[fnClusterTTBins] = th;   fnClusterTTBins++;  }

  void        SetMode(Int_t mode){ fMode = mode;}       //Analysi mode normal=0 or embedded event =1 

  void        SetPhiCut(Double_t pcut){ fPhiCut = TMath::Pi()-pcut; }

  void        SetMinFractionShared(Double_t fr){ fMinFractionShared = fr;} 

  void        SetJetRadius(Double_t jr) { fJetR = jr;}                               

  void        SetJetAcut(Double_t ac){ fJetAcut = ac;} 



  Bool_t      PassedGATrigger();
  Bool_t      PassedMinBiasTrigger();
  Bool_t      PassedHighMultTrigger();
  Int_t       GetMaxDistanceFromBorder(AliVCluster* cluster);
  Bool_t      FinalClusterCuts(AliVCluster* cluster);
//  std::string MatchTrigger(const std::string &triggerstring);

  Bool_t      RetrieveEventObjects();
  Bool_t      Run();
  Bool_t      FillHistograms();


 private:


  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isGen=0);  
  Bool_t      IsEventInAcceptance(AliVEvent* event);     


  // ######### STANDARD FUNCTIONS
  void        ExecOnceLocal();                    

  void        GetMyRho(AliJetContainer* ktjets, Double_t &rhokt, Double_t &rhocms);
  

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
  TString     fMyJetDetLevelContainerName;            // name of detector level level MC jet container which is used while embeding
  TString     fMyClusterContainerName;                // name of detector level jet container 
  TString     fMyKTJetContainerName;                  // name of KT detector level jet container 
  TString     fMyKTJetParticleContainerName;          // name of KT particle level MC jet container 
  TString     fMyKTJetDetLevelContainerName;          // name of KT detector level level MC jet container which is used while embeding


  AliTrackContainer    *fTrkContainerDetLevel;        //! detector level track container   (or tracks in combined embedded events)
  AliParticleContainer *fParticleContainerPartLevel;  //! particle level container with pythia particles
  AliTrackContainer    *fTrkContainerDetLevelEMB;     //! detector level container with pythia tracks  used for embedding
  AliJetContainer      *fJetContainerDetLevel;        //! detector level jet container   (or jets in combined events after embedding) 
  AliJetContainer      *fJetContainerPartLevel;       //! particle level jet container
  AliJetContainer      *fJetContainerDetLevelEMB;     //! detector level jet container with pythia detector level jets  

  AliClusterContainer  *fClusterContainerDetLevel;    //! detector level EMCAL cluster container   


  AliJetContainer      *fKTJetContainerDetLevel;        //! KT detector level jet container   (or jets in combined events after embedding) 
  AliJetContainer      *fKTJetContainerPartLevel;       //! KT particle level jet container
  AliJetContainer      *fKTJetContainerDetLevelEMB;     //! KT detector level jet container with pythia detector level jets  

  // ########## CENTRALITY
  AliMultSelection*  fMultSelection;                  //! object which handels centrality 

  Bool_t  fIsMinBiasTrig;                             //! triggered by Min Bias Trig
  Bool_t  fIsEmcalTrig;                               //! triggered by EMCAL
  Bool_t  fIsHighMultTrig;                            //! triggered by high multiplicity trigger 
                                                      
  Float_t fCentralityV0A;                             //! Centrality from V0A
  Float_t fCentralityV0C;                             //! Centrality from V0C
  Float_t fCentralityV0M;                             //! Centrality from V0M
                                                      
  Double_t fxVertex;                                  //!  X vertex from ITS
  Double_t fyVertex;                                  //!  Y vertex from ITS
  Double_t fzVertex;                                  //!  Z vertex from ITS

  Int_t    fNTracklets;                               //!  no. tracklets
 
  Double_t  fMultV0A;                                  //!  mult. V0A
  Double_t  fMultV0C;                                  //!  mult. V0C
  Double_t  fMultV0M;                                  //!  mult. V0A+V0C

  Double_t  fMultV0Mnorm;                              //!  mult. V0M/mean
  Double_t  fAsymV0M;                                  //!  V0A-V0C / V0A + V0C


  Double_t  fMultV0A_PartLevel;                        //!  mult. V0A       particle level
  Double_t  fMultV0C_PartLevel;                        //!  mult. V0C       particle level
  Double_t  fMultV0M_PartLevel;                        //!  mult. V0A+V0C   particle level

  Double_t  fMultV0Mnorm_PartLevel;                    //!  mult. V0M normalized by mean V0M  particle level
  Double_t  fAsymV0M_PartLevel;                        //!  V0A-V0C / V0A + V0C


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
   TH2D     *fhClusterPhiIncl[kTG];                      //!  phi inclusive cluster
   TH2D     *fhClusterEtaIncl[kTG];                      //!  eta inclusive

   THnSparse *fhTrackPtEtaPhiV0norm[kTG];               //!  pt, eta, phi, V0M for inclusive tracks 
   THnSparse *fhJetPtEtaPhiV0norm[kTG];                  //!  pt, eta, phi, V0M for inclusive tracks 
   THnSparse *fhJetPtEtaPhiV0normTTH[kTG][fkTTbins];    //!  pt, eta, phi, V0M for inclusive tracks 
   THnSparse *fhJetPtAreaV0norm[kTG][kRho];             //!  pt, area V0M for inclusive jets 
   THnSparse *fhJetPtAreaV0norm_PartLevel[kRho];        //!  pt, area V0M for inclusive jets 
 
   TH1D     *fhRho[kTG][kRho];                          //! minimum bias rho inclusive
   TH1D     *fhRhoTTH[kTG][fkTTbins][kRho];             //! in events MB with hadron TT
   TH1D     *fhRhoTTJ[kTG][fkTTbins];                   //! in events MB with jet TT
   TH1D     *fhRhoTTC[kTG][fkTTbins][kRho];             //! in events MB with cluster TT

   TH1D     *fhRhoMBpart[kRho];                         //! minimum bias rho inclusive particle level
   TH1D     *fhRhoTTHinMBpart[fkTTbins][kRho];          //! in events MB with hadron TT particle level
   TH1D     *fhRhoTTCinMBpart[fkTTbins][kRho];          //! in events MB with cluster TT particle level

   TH1D* fhVertex[fkVtx];                             //! vertex distribution 
                                                      
   TH2D* fhCentrality[kTG][fkCE];                      //! estimated centrality based on  mult V0, mult VC, V0M 
   TH2D* fhCentralityTTH[kTG][fkCE][fkTTbins];         //! estimated centrality  biased with hadron TT 
   TH2D* fhCentralityTTJ[kTG][fkCE][fkTTbins];         //! estimated centrality  biased with ch jet TT 
   TH2D* fhCentralityTTC[kTG][fkCE][fkTTbins];         //! estimated centrality  biased with cluster TT
                                                      
   TH1D* fhSignal[kTG][fkCE];                          //! centrality estimators:  mult V0, mult VC, tracklets, znatower0, znctower0, V0M, V0Mnorm  in MB
   TH1D* fhSignalTTH[kTG][fkCE][fkTTbins];             //! distributions of centrality estimators biased with hadron TT in min bias
   TH1D* fhSignalTTJ[kTG][fkCE][fkTTbins];             //! distributions of centrality estimators biased with ch jet TT in min bias
   TH1D* fhSignalTTC[kTG][fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in min bias 


   TH1D* fhSignal_PartLevel[fkCE];                 //! particle level centrality estimators:  mult V0, mult VC, V0M, V0Mnorm  in MB
   TH1D* fhSignalTTH_PartLevel[fkCE][fkTTbins];  //! particle level distributions of centrality estimators biased with hadron TT in min bias
   TH1D* fhSignalTTC_PartLevel[fkCE][fkTTbins];  //! particle level distributions of centrality estimators biased with cluster TT in min bias

   TH2D* fhV0MAssymVsV0Mnorm[kTG];                    //! V0AC asymmetry versus V0Mnorm  in inclusive events
   TH2D* fhV0MAssymVsV0Mnorm_PartLevel;               //! V0AC asymmetry versus V0Mnorm  in inclusive events

   TH2D* fhV0MAssymVsV0MnormTTH[kTG][fkTTbins];       //! V0AC asymmetry versus V0Mnorm  in TTH events 
   TH2D* fhV0MAssymVsV0MnormTTH_PartLevel[fkTTbins];  //! V0AC asymmetry versus V0Mnorm  in TTH events


   TH2D* fhV0AvsV0C;                                   //! V0A vs V0C in MB 
//   TH2D* fhV0MvsV0Mnorm;                               //! V0M vs V0Mnorm in MB 
   TH2D* fhV0AvsSPD;                                   //! V0A vs SPD in MB 
   TH2D* fhV0CvsSPD;                                   //! V0C vs SPD in MB 
   TH2D* fhV0AvsV0CTTH[fkTTbins];                      //! V0A vs V0C biased with hadron TT
   TH2D* fhV0AvsV0CTTJ[fkTTbins];                      //! V0A vs V0C biased with ch jet TT
   TH2D* fhV0AvsV0CTTCinMB[fkTTbins];                  //! V0A vs V0C biased with cluster TT in min bias 
   TH2D* fhV0AvsV0CTTCinGA[fkTTbins];                  //! V0A vs V0C biased with cluster TT in Gamma trigger
                                                      
   TH1D* fhMultTTH[kTG][fkTTbins];                       //! multiplicity of hadron TT 
   TH1D* fhMultTTJ[kTG][fkTTbins];                       //! multiplicity of charged jet TT 
   TH1D* fhMultTTC[kTG][fkTTbins];                       //! multiplicity of cluster TT 

   TH2D* fhTrackMult[kTG];                                 //! multiplicity of midrapidity charged tracks 
   TH2D* fhMeanTrackPt[kTG];                               //! mean track pT 

   //hadron TT
//   TH2D* fhTTH_CentV0M[kTG][fkTTbins];                    //! counter of semi-inclusive hadron TT versus V0M    centrality
   TH2D* fhTTH_V0Mnorm1[kTG][fkTTbins];                   //! counter of semi-inclusive hadron TT versus V0M/mean 

   TH2D* fhTTH_V0Mnorm1_PartLevel[fkTTbins];         //! counter of semi-inclusive hadron TT   V0M/mean particle level

   TH3D* fhTTH_3D_V0Mnorm1[kTG][fkTTbins];           //! counter of semi-inclusive hadron l TT in MB versus V0M/mean

   TH3D* fhTTH_3D_V0Mnorm1_PartLevel[fkTTbins];      //! counter of semi-inclusive hadron TT in MB versus V0M/mean particle level



   //EMCAL cluster TT
//   TH2D* fhTTC_CentV0M[kTG][fkTTbins];                    //! counter of semi-inclusive emcal TT in MB versus V0M  centrality 
   TH2D* fhTTC_V0Mnorm1[kTG][fkTTbins];                   //! counter of semi-inclusive emcal TT in MB versus V0M/mean

   TH2D* fhTTC_V0Mnorm1_PartLevel[fkTTbins];            //! counter of semi-inclusive emcal TT in MB versus V0M/mean particle level


   //recoil jet yields with hadron TT
//   TH2D* fhRecoilJetPtTTH_CentV0M[kTG][fkTTbins];         //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M centrality  
   TH2D* fhRecoilJetPtTTH_V0Mnorm1[kTG][fkTTbins][kRho];       //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M/mean
   TH2D* fhRecoilJetPtTTH_V0Mnorm1_PartLevel[fkTTbins][kRho];  //! pT spectrum of recoil jets associated to semi-inclusive hadron TT versus V0M/mean

   TH3D* fhRecoilJetPhiTTH_V0Mnorm1[kTG][fkTTbins][kRho];        //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  
   TH3D* fhRecoilJetPhiTTH_V0Mnorm1_PartLevel[fkTTbins][kRho];   //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  minimum bias particle level 

   THnSparse *fhRecoilJetTTH_V0Mnorm1[kTG][fkTTbins][kRho];      //! recoil jet  (V0M/mean, V0 assymetry , recoil jet pT,  delta phi)
   THnSparse *fhRecoilJetTTH_V0Mnorm1_PartLevel[fkTTbins][kRho]; //! recoil jet  (V0M/mean, V0 assymetry , recoil jet pT,  delta phi)


   TH2D* fhRecoilJetPtTTHref_V0Mnorm1_rhoShift1[kTG][kRho];   //! reference pT spectrum of recoil jets TT67 vs V0M/mean with rho shifted for TT12,20
   TH2D* fhRecoilJetPtTTHref_V0Mnorm1_rhoShift2[kTG][kRho];   //! reference pT spectrum of recoil jets TT67 vs V0M/mean with rho shifted for TT20,30
   
   //recoil jet yields with EMCAL cluster TT
//   TH2D* fhRecoilJetPtTTC_CentV0M[kTG][fkTTbins];         //! pT spectrum of recoil jets associated to semi-inclusive emcal TT versus V0M    centrality 
   TH2D* fhRecoilJetPtTTC_V0Mnorm1[kTG][fkTTbins][kRho];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT versus V0M/mean

   TH2D* fhRecoilJetPtTTC_V0Mnorm1_PartLevel[fkTTbins][kRho];   //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M/mean


//   TH2D* fhDeltaPtTTH_RC_CentV0M[kTG][fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT versus V0M centrality  
//   TH2D* fhDeltaPtTTC_RC_CentV0M[kTG][fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT versus V0M centrality   

   TH2D* fhDeltaPtTTH_RC_V0Mnorm1[kTG][fkTTbins][kRho];         //! delta pT spectrum from random cones  in events with hadron TT  versus V0M/mean V0M  
   TH2D* fhDeltaPtTTC_RC_V0Mnorm1[kTG][fkTTbins][kRho];         //! delta pT spectrum from random cones  in events with emcal  TT  versus V0M/mean V0M  



   TH2D* fhDeltaPtTTH_RC_V0Mnorm1_PartLevel[fkTTbins][kRho];  //! delta pT spectrum from random cones  in events with  hadron TT in MB versus V0M/mean V0M   PARTICLE LEVEL
   TH2D* fhDeltaPtTTC_RC_V0Mnorm1_PartLevel[fkTTbins][kRho];  //! delta pT spectrum from random cones  in events with  emcal  TT in MB versus V0M/mean V0M   PARTICLE LEVEL


   TH3D* fhPtTrkTruePrimGen;                            //! physical primary mc particle eta vs pT  vs V0Mnorm
   TH3D* fhPtTrkTruePrimRec;                            //! physical primary detector level track eta vs pT vs V0Mnorm
   TH3D* fhPtTrkSecOrFakeRec;                           //! secondary tracks eta vs pT vs V0Mnorm

   TH1D* fhJetPtPartLevelCorr[kRho];                          //! response matrix normalization spectrum, jet pT corrected on rho
   TH1D* fhJetPtPartLevelZero;                          //! response matrix normalization spectrum, jet pT is not corrected on rho
   TH1D* fhJetPtPartLevelCorrTTHdl[fkTTbins][kRho];           //! response matrix normalization spectrum, events with det. level TTH 

   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorr[kRho];           //! response matrix jet pT corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelZero;           //! response matrix jet pT not corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[fkTTbins][kRho];  //! response matrix events with detector level TTH

   TH2D* fhJetPtResolutionVsPtPartLevel[kRho];                //! resolution of jet pT

   TH2D* fhOneOverPtVsPhiNeg;                           //! 1/p_T,track  versus phi for negative tracks //AID//             
   TH2D* fhOneOverPtVsPhiPos;                           //! 1/p_T,track  versus phi for positive tracks //AID//             
   TH2D* fhSigmaPtOverPtVsPt;                           //! resolution of 1/p_T,track  versus p_T,track //AID//             
   TH2D* fhDCAinXVsPt;                                  //! X DCA versus pT  //AID// 
   TH2D* fhDCAinYVsPt;                                  //! Y DCA versus pT  //AID//
   TH2D* fhDCAinXVsPtPhysPrimary;                       //! X DCA versus pT for physical primaries  //AID// 
   TH2D* fhDCAinYVsPtPhysPrimary;                       //! Y DCA versus pT for physical primaries  //AID//
   TH2D* fhDCAinXVsPtSecondary;                         //! X DCA versus pT for secondaries //AID// 
   TH2D* fhDCAinYVsPtSecondary;                         //! Y DCA versus pT for secondaries //AID//
   TH2D* fhFractionOfSecInJet[kRho];                    //! Fraction of jet pT carried by secondaries //AID//

   TH2D* fhV0ARunByRunMB;                               //! run by run V0M
   TH2D* fhV0CRunByRunMB;                               //! run by run V0M
   TH2D* fhV0MRunByRunMB;                               //! run by run V0M
   TH2D* fhV0MnormRunByRunMB;                           //! run by run V0M norm

   //EMBEDDING
   TH2D* fhTrackEtaInclEMB;                              //!  Eta dist inclusive embedded tracks vs pT 
   TH3D* fhRecoilJetPhiTTH_EMB_V0Mnorm1[kTG][fkTTbins][kRho];  //!  filled with any detector level pythia recoil jet 
   TH3D* fhRecoilJetPhiTTH_TAG_V0Mnorm1[kTG][fkTTbins][kRho];  //!  filled  tagged closest detector level pythia recoil jet 

   TH1D* fhJetPtPartLevelCorr_EMB[kTG][kRho];                   //! response matrix normalization spectrum, jet pT corrected on rho
   TH1D* fhJetPtPartLevelZero_EMB[kTG];                   //! response matrix normalization spectrum, jet pT is not corrected on rho
   TH1D* fhJetPtPartLevelCorrTTHdl_EMB[kTG][fkTTbins][kRho];    //! response matrix normalization spectrum, events with det. level TTH 
   TH1D* fhJetPtPartLevelZeroTTHdl_EMB[kTG][fkTTbins];    //! response matrix normalization spectrum, events with det. level TTH 


   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorr_EMB[kTG][kRho];           //! response matrix jet pT corrected on rho    embedded to minimum bias events
   TH2D* fhJetPtPartLevelVsJetPtDetLevelZero_EMB[kTG][kRho];           //! response matrix jet pT not corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl_EMB[kTG][fkTTbins][kRho];  //! response matrix events with detector level TTH
   TH2D* fhJetPtPartLevelVsJetPtDetLevelZeroTTHdl_EMB[kTG][fkTTbins][kRho];  //! response matrix events with detector level TTH

   Double_t fMinFractionShared;     // cut on shared fraction

   TH2D *fhSharedJetFraction[kTG];  //!  shared fraction

   TH1F *fhTrialsEMBtot[kTG];        //! number of trials    after  event selection in      embedding
   TProfile *fhXsectionEMBtot[kTG];  //! Xsection     after  event selection in      embedding
   TH1F *fhTrialsEMB[kTG];        //! number of trials    after  event selection in      embedding
   TProfile *fhXsectionEMB[kTG];  //! Xsection     after  event selection in      embedding
   TH1F *fhPtHardEMB[kTG];        //! pthard

   Double_t fZVertexCut;                              // vertex cut in z 
                                                     
   Int_t    fnHadronTTBins;                           // number of TT bins charged hadron 
   Int_t    fnJetChTTBins;                            // number of TT bins charged jet
   Int_t    fnClusterTTBins;                          // number of TT bins gamma cluster
   Int_t    fHadronTTLowPt[fkTTbins];                 // low pt TT range charged hadron 
   Int_t    fHadronTTHighPt[fkTTbins];                // high pt TT range charged hadron 
   Int_t    fJetChTTLowPt[fkTTbins];                  // low pt TT range charged jet
   Int_t    fJetChTTHighPt[fkTTbins];                 // high pt TT range charged jet
   Int_t    fClusterTTLowPt[fkTTbins];                // low pt TT range charged jet
   Int_t    fClusterTTHighPt[fkTTbins];               // high pt TT range charged jet
                                                     
                                                     
   Int_t    fHadronTT[fkTTbins];                      //! array which stores the number of hadron TT in given event 
   Int_t    fJetChTT[fkTTbins];                       //! array which stores the number of jets TT in given event 
   Int_t    fClusterTT[fkTTbins];                     //! array which stores the number of cluster TT in given event 

   Int_t    fHadronTT_PartLevel[fkTTbins];             //! array which stores the number of hadron TT in given event particle level 
   Int_t    fClusterTT_PartLevel[fkTTbins];            //! array which stores the number of gamma TT in given event particle level 
                                                      
   Int_t    fMode;                                    // Analysis mode   0=real data, 1=mc pythia, 2=embedded 
   AliEMCALRecoUtils          *fFiducialCellCut;      //!

   PWGJE::EMCALJetTasks::AliAnalysisEmcalJetHelperEA  *fHelperEA;                   // wrapper for  mean V0 multiplicities
   Double_t fMeanV0M;                                 //! mean V0C signal in incl. MB  run by run
   Double_t fMeanV0M_PartLevel;                       // mean V0M signal in incl. MB particle level 

   Int_t fIndexTTC[fkTTbins];                      //! index of the chosen EMCAL cluster trigger 
   Int_t fIndexTTH[fkTTbins];                      //! index of the chosen hadron trigger 
   Int_t fIndexTTJ[fkTTbins];                      //! index of the chosen jet trigger 

   Int_t fdeltapT[fkTTbins][kRho];                       //! delta pT detector level
   Int_t fdeltapT_PartLevel[fkTTbins][kRho];             //! delta pT particle level

   Int_t fIndexTTH_PartLevel[fkTTbins];             //! index of the chosen hadron trigger particle level 
   Int_t fIndexTTC_PartLevel[fkTTbins];             //! index of the chosen gamma trigger particle level 

   vector<TLorentzVector> fTTC[fkTTbins];        //!  EMCAL cluster trigger candidate 
   vector<TLorentzVector> fTTH[fkTTbins];        //!  hadron trigger candidate 
   vector<TLorentzVector> fTTJ[fkTTbins];        //!  jet trigger candidate 

   vector<TLorentzVector> fTTH_PartLevel[fkTTbins];  //!  hadron trigger candidate particle level 
   vector<TLorentzVector> fTTC_PartLevel[fkTTbins];  //!  gamma trigger candidate particle level 
  
   Bool_t fFillSigTT;                            //! flag which labels whether the event should be filled for ref or sig class 

   Double_t fPhiCut;                             // phi angle cut on the recoil jet  pi-fPhiCut
   TRandom3* fRandom;                            //! Radom 

   Double_t frhoveckt[999];                        //! auxiliary array to store pT/A of kT jets for kt
   Double_t frhovecms[999];                        //! auxiliary array to store pT/A of kT jets for cms

   Double_t fJetR;                                 // jet radius
   Double_t fJetAcut;                              // jet area cut

   AliAnalysisTaskEA(const AliAnalysisTaskEA&);
   AliAnalysisTaskEA& operator=(const AliAnalysisTaskEA&);

   ClassDef(AliAnalysisTaskEA, 21); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
}
}
#endif
