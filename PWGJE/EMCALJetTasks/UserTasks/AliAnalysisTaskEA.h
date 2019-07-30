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

#include <vector>
using std::vector;

#include "AliAnalysisTaskEmcalJet.h"
#include "AliEventCuts.h"
#include "AliFiducialCut.h"
#include "AliEMCALRecoUtils.h"


// ANALYSIS OF EVENT ACTIVITY WITH HIGH PT HADRON BIAS
// Author Filip Krizek   (25 OCT 2018)

class AliAnalysisTaskEA : public AliAnalysisTaskEmcalJet {
   public:

   enum Analysis {
      kDetLevel=0, 
      kPartLevel=1,
      fkVtx=3,
      fkTTbins = 20 //maximum number of TT bins
   };

   enum {fkV0A, fkV0C, fkV0M, fkV0Mnorm1, fkV0Mnorm2,fkCE}; //detector level   norm1 : divided by mean V0M    norm2: average  A/mean A  +C/ mean C


   enum {kpp=0, kpPb=1};


   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskEA();
   AliAnalysisTaskEA(const char *name);
   virtual  ~AliAnalysisTaskEA();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   static AliAnalysisTaskEA* AddTaskEA(
       Int_t       system             = AliAnalysisTaskEA::kpPb, // collision system 
       const char* jetarrayname       = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets
       const char* jetarraynameMC     = "Jet_AKTChargedR040_mcparticles_pT0150_pt_scheme", //name of jet TClones array for MC particle level jets
       const char* trackarrayname     = "tracks", //name of track TClonesArray for detector level jets
       const char* mcpariclearrayname = "mcparticles", //name of track TClonesArray array for MC particle level jets
       const char* clusterarrayname   = "caloClusters", //name of EMCAL cluster TClonesArray array  for detector level
       const char* rhoname            = "", //name of track TClonesArray for detector level jets
       const char* mcrhoname          = "", //name of track TClonesArray array for MC particle level jets
       Double_t    jetRadius          = 0.4,  //radius of analyzed jets
       UInt_t      trigger            = AliVEvent::kAny,  //trigger
       Int_t       isMC               = 0,     // 0=real data    , 1= particle+detector level simulation 
       Double_t    trackEtaWindow     = 0.9,   //pseudorapidity range for tracks
       Bool_t      useVertexCut       = kTRUE,  // vertex cut
       Bool_t      usePileUpCut       = kTRUE, // discard pile up event
       Double_t    acut               = 0.6,   //cut on relative jet area
       Double_t    emcaltofcut        = 30e-9,   //cut on relative jet area
       const char* suffix = ""                              //SUBWAGON has to be the last parameter
  );

   
  // ######### SETTERS/GETTERS
  void        SetMC(Int_t bMC){
                 fMC = bMC;
              }

  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}  
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;} 

  void        SetAcceptanceWindows(Double_t trackEta){ 
                 fTrackEtaWindow  = trackEta; 
              } 


  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void        SetExternalRhoTaskName(const char* name) {fRhoTaskName = name;}
  void        SetExternalRhoTaskNameMC(const char* name) {fRhoTaskNameMC = name;}

  void        SetTrackContainerName(const char* name){ fMyTrackContainerName = name;}
  void        SetMCParticleContainerName(const char* name){ fMyParticleContainerName = name;}
  void        SetClusterContainerName(const char* name){ fMyClusterContainerName = name;}

  void        SetJetContainerName(const char* name){ fMyJetContainerName = name;}
  void        SetMCJetContainerName(const char* name){ fMyJetParticleContainerName = name;}
  void        SetHadronTT(Int_t tl,Int_t th){ fHadronTTLowPt[fnHadronTTBins] = tl; fHadronTTHighPt[fnHadronTTBins] = th; fnHadronTTBins++; }
  void        SetJetChTT(Int_t tl,Int_t th){  fJetChTTLowPt[fnJetChTTBins] = tl;   fJetChTTHighPt[fnJetChTTBins] = th;   fnJetChTTBins++;  }
  void        SetClusterTT(Int_t tl,Int_t th){  fClusterTTLowPt[fnClusterTTBins] = tl;   fClusterTTHighPt[fnClusterTTBins] = th;   fnClusterTTBins++;  }
  void        SetFillTTree(Bool_t b){ fFillTTree = b; } //fill output TTree
  void        SetSystem(Int_t sys){ fSystem = sys;}     // Collision system pp or pP pp or pPb 

  void        SetMeanV0APartLevel(Double_t mva){ fMeanV0A_PartLevel = mva; }
  void        SetMeanV0CPartLevel(Double_t mvc){ fMeanV0C_PartLevel = mvc; }
  void        SetMeanV0MPartLevel(Double_t mvm){ fMeanV0M_PartLevel = mvm; }
 
  void        SetPhiCut(Double_t pcut){ fPhiCut = TMath::Pi()-pcut; }

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

  Double_t    GetExternalRho(Bool_t isMC); 

  Double_t    GetDeltaPt(Double_t phiTT, Double_t etaTT, Double_t phiLJ, Double_t etaLJ, Double_t phiSJ, Double_t etaSJ, Double_t rho, Int_t level);
  
  // ########## USAGE TRIGGERS 
  Bool_t      fUseDefaultVertexCut;                   // trigger if automatic vertex cut from helper class should be done 
  Bool_t      fUsePileUpCut;                          // trigger if pileup cut should be done
  

  // ########## SOURCE INFORMATION
  TString     fMyTrackContainerName;                  // name of detector level track container 
  TString     fMyParticleContainerName;               // name of particle level MC particle container
  TString     fMyJetContainerName;                    // name of detector level jet container 
  TString     fMyJetParticleContainerName;            // name of particle level MC jet container 
  TString     fMyClusterContainerName;                // name of detector level jet container 

  AliTrackContainer    *fTrkContainerDetLevel;        //! detector level track container
  AliParticleContainer *fParticleContainerPartLevel;  //! particle level container with particles
  AliJetContainer      *fJetContainerDetLevel;        //! detector level jet container  
  AliJetContainer      *fJetContainerPartLevel;       //! particle level jet container
  AliClusterContainer  *fClusterContainerDetLevel;    //! detector level EMCAL cluster container   

  TString     fRhoTaskName;                           // name of rho CMS bg task for this analysis
  TString     fRhoTaskNameMC;                         // MC name of rho CMS bg task for this analysis


  // ########## CENTRALITY
  TTree   *fCentralityTree;                           //! output tree
  AliMultSelection*  fMultSelection;                  //! object which handels centrality 

  char    fTrigClass[1000];                           //! fired trigger classes
  Bool_t  fIsMinBiasTrig;                             //! triggered by Min Bias Trig
  Bool_t  fIsEmcalTrig;                               //! triggered by EMCAL
  Bool_t  fIsHighMultTrig;                            //! triggered by high multiplicity trigger 
                                                      
  Float_t fCentralityV0A;                             //! Centrality from V0A
  Float_t fCentralityV0C;                             //! Centrality from V0C
  Float_t fCentralityV0M;                             //! Centrality from V0M
  Float_t fCentralityCL1;                             //! Centrality from Clusters in layer 1
  Float_t fCentralityZNA;                             //! Centrality from ZNA
  Float_t fCentralityZNC;                             //! Centrality from ZNC
                                                      
  Double_t fxVertex;                                  //!  X vertex from ITS
  Double_t fyVertex;                                  //!  Y vertex from ITS
  Double_t fzVertex;                                  //!  Z vertex from ITS
  Bool_t   fVertexer3d;                               //!  Is vertex from 3d vertexer?
  //                                                  
  Int_t    fNTracklets;                               //!  no. tracklets
  Int_t    fNClusters[2];                             //!  no. clusters on SPD layers
  //                                                  
  Int_t    fIsV0ATriggered;                           //!  VOA decision
  Int_t    fIsV0CTriggered;                           //!  V0C decision
  Double_t  fMultV0A;                                  //!  mult. V0A
  Double_t  fMultV0C;                                  //!  mult. V0C
  Double_t  fMultV0M;                                  //!  mult. V0A+V0C
  Double_t  fMultV0Anorm;                              //!  mult. V0A normalized by mean V0A
  Double_t  fMultV0Cnorm;                              //!  mult. V0C normalized by mean V0C
  Double_t  fMultV0Mnorm;                              //!  mult. V0M/mean
  Double_t  fMultV0AV0Cnorm;                           //!  mult. (V0A/mean+V0C/mean)/2 

  Double_t  fMultV0A_PartLevel;                        //!  mult. V0A       particle level
  Double_t  fMultV0C_PartLevel;                        //!  mult. V0C       particle level
  Double_t  fMultV0M_PartLevel;                        //!  mult. V0A+V0C   particle level
  Double_t  fMultV0Anorm_PartLevel;                    //!  mult. V0A normalized by mean V0A   particle level
  Double_t  fMultV0Cnorm_PartLevel;                    //!  mult. V0C normalized by mean V0C   particle level
  Double_t  fMultV0Mnorm_PartLevel;                    //!  mult. V0M normalized by mean V0M  particle level
  Double_t  fMultV0AV0Cnorm_PartLevel;                 //!  mult. (V0A/mean+V0C/mean)/2   particle level

  Float_t  fRingMultV0[8];                            //!  V0 ring mult.
  //                                                  
  Float_t  fZEM1Energy;                               //!  ZEM1 Energy
  Float_t  fZEM2Energy;                               //!  ZEM2 Energy
                                                      
  Float_t  fZNCtower[5];                              //!  ZNC 5 tower signals
  Float_t  fZPCtower[5];                              //!  ZPC 5 tower signals
  Float_t  fZNAtower[5];                              //!  ZNA 5 tower signals
  Float_t  fZPAtower[5];                              //!  ZPA 5 tower signals
  Float_t  fZNCtowerLG[5];                            //!  ZNC 5 tower signals
  Float_t  fZPCtowerLG[5];                            //!  ZPC 5 tower signals
  Float_t  fZNAtowerLG[5];                            //!  ZNA 5 tower signals
  Float_t  fZPAtowerLG[5];                            //!  ZPA 5 tower signals


  // ########## CUTS 
  Double_t    fTrackEtaWindow;                        // +- window in eta for tracks  
  Double_t    fMinTrackPt;                            // Min track pt to be accepted  

  // ########## GENERAL ////VARS
  Bool_t              fMC;                            //  real data  or MC  flag
  AliAnalysisUtils*   fHelperClass;                   //! Vertex selection helper
  Bool_t              fInitializedLocal;              //! trigger if tracks/jets are loaded  initiates calling   ExecOnce 


   TH1D               *fHistEvtSelection;             //! gc event statistics

   TH1D     *fhVertexZall;                            //! gc vertexZ inclusive before cut
   TH1D     *fhVertexZ;                               //! gc vertexZ inclusive after cut
                                                      
   TH2D     *fhTrackPhiInclMB;                        //! minimum bias phi inclusive ch hadron
   TH2D     *fhTrackEtaInclMB;                        //! minimum bias eta inclusive minimum bias
   TH2D     *fhTrackEtaInclHM;                        //! minimum bias eta inclusive high multiplicity
   TH2D     *fhJetPhiIncl;                            //! minimum bias phi inclusive ch jet
   TH2D     *fhJetEtaIncl;                            //! minimum bias eta inclusive
   TH2D     *fhClusterPhiInclMB;                      //! minimum bias phi inclusive cluster
   TH2D     *fhClusterEtaInclMB;                      //! minimum bias eta inclusive
   TH2D     *fhClusterPhiInclGA;                      //! minimum bias phi inclusive cluster
   TH2D     *fhClusterEtaInclGA;                      //! minimum bias eta inclusive
 
   TH1D     *fhRhoMB;                                 //! minimum bias rho inclusive
   TH1D     *fhRhoHM;                                 //! high mult rho inclusive
   TH1D     *fhRhoTTHinMB[fkTTbins];                  //! in events MB with hadron TT
   TH1D     *fhRhoTTHinHM[fkTTbins];                  //! in events HM with hadron TT
   TH1D     *fhRhoTTJinMB[fkTTbins];                  //! in events MB with jet TT
   TH1D     *fhRhoTTJinHM[fkTTbins];                  //! in events HM with jet TT
   TH1D     *fhRhoTTCinMB[fkTTbins];                  //! in events MB with cluster TT
   TH1D     *fhRhoTTCinHM[fkTTbins];                  //! in events HM with cluster TT
   TH1D     *fhRhoTTCinGA[fkTTbins];                  //! in events GA with cluster TT

   TH1D     *fhRhoMBpart;                             //! minimum bias rho inclusive particle level
   TH1D     *fhRhoTTHinMBpart[fkTTbins];              //! in events MB with hadron TT particle level
   TH1D     *fhRhoTTCinMBpart[fkTTbins];              //! in events MB with cluster TT particle level

   TH1D* fhVertex[fkVtx];                             //! vertex distribution 
   TH1D* fhVertexTTH[fkVtx][fkTTbins];                //! vertex distribution in events biased with hadron TT
                                                      
   TH2D* fhCentralityMB[fkCE];                        //! estimated centrality based on  mult V0, mult VC, V0M 
   TH2D* fhCentralityHM[fkCE];                        //! estimated centrality based on  mult V0, mult VC, V0M 
   TH2D* fhCentralityTTHinMB[fkCE][fkTTbins];         //! estimated centrality in MB events biased with hadron TT 
   TH2D* fhCentralityTTHinHM[fkCE][fkTTbins];         //! estimated centrality in HM events biased with hadron TT 
   TH2D* fhCentralityTTJinMB[fkCE][fkTTbins];         //! estimated centrality in MB events biased with ch jet TT 
   TH2D* fhCentralityTTJinHM[fkCE][fkTTbins];         //! estimated centrality in HM events biased with ch jet TT 
   TH2D* fhCentralityTTCinMB[fkCE][fkTTbins];         //! estimated centrality in MB events biased with cluster TT
   TH2D* fhCentralityTTCinGA[fkCE][fkTTbins];         //! estimated centrality in GA events biased with cluster TT
                                                      
   TH1D* fhSignalMB[fkCE];                            //! centrality estimators:  mult V0, mult VC, tracklets, znatower0, znctower0, V0M, V0Mnorm  in MB
   TH1D* fhSignalHM[fkCE];                            //! distributions of centrality estimators in HM events
   TH1D* fhSignalTTHinMB[fkCE][fkTTbins];             //! distributions of centrality estimators biased with hadron TT in min bias
   TH1D* fhSignalTTHinHM[fkCE][fkTTbins];             //! distributions of centrality estimators biased with hadron TT in high mult
   TH1D* fhSignalTTJinMB[fkCE][fkTTbins];             //! distributions of centrality estimators biased with ch jet TT in min bias
   TH1D* fhSignalTTJinHM[fkCE][fkTTbins];             //! distributions of centrality estimators biased with ch jet TT in high mult
   TH1D* fhSignalTTCinMB[fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in min bias 
   TH1D* fhSignalTTCinHM[fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in high mult 
   TH1D* fhSignalTTCinGA[fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in Gamma trigger


   TH1D* fhSignalMB_PartLevel[fkCE];                 //! particle level centrality estimators:  mult V0, mult VC, V0M, V0Mnorm  in MB
   TH1D* fhSignalTTHinMB_PartLevel[fkCE][fkTTbins];  //! particle level distributions of centrality estimators biased with hadron TT in min bias
   TH1D* fhSignalTTCinMB_PartLevel[fkCE][fkTTbins];  //! particle level distributions of centrality estimators biased with cluster TT in min bias

   TH2D* fhV0AvsV0C;                                   //! V0A vs V0C in MB 
   TH2D* fhV0MvsV0Mnorm;                               //! V0M vs V0Mnorm in MB 
   TH2D* fhV0AvsSPD;                                   //! V0A vs SPD in MB 
   TH2D* fhV0CvsSPD;                                   //! V0C vs SPD in MB 
   TH2D* fhV0AvsV0CTTH[fkTTbins];                      //! V0A vs V0C biased with hadron TT
   TH2D* fhV0AvsV0CTTJ[fkTTbins];                      //! V0A vs V0C biased with ch jet TT
   TH2D* fhV0AvsV0CTTCinMB[fkTTbins];                  //! V0A vs V0C biased with cluster TT in min bias 
   TH2D* fhV0AvsV0CTTCinGA[fkTTbins];                  //! V0A vs V0C biased with cluster TT in Gamma trigger
                                                      
   TH1D* fhMultTTHinMB[fkTTbins];                       //! multiplicity of hadron TT in MB event
   TH1D* fhMultTTHinHM[fkTTbins];                       //! multiplicity of hadron TT in HM event
   TH1D* fhMultTTJinMB[fkTTbins];                       //! multiplicity of charged jet TT in MB event
   TH1D* fhMultTTJinHM[fkTTbins];                       //! multiplicity of charged jet TT in HM event
   TH1D* fhMultTTCinMB[fkTTbins];                       //! multiplicity of cluster TT in MB event
   TH1D* fhMultTTCinHM[fkTTbins];                       //! multiplicity of cluster TT in HM event
   TH1D* fhMultTTCinGA[fkTTbins];                       //! multiplicity of cluster TT in Gamma trigger event

   TH2D* fhTrackMultMB;                                 //! multiplicity of midrapidity charged tracks for MB events 
   TH2D* fhTrackMultHM;                                 //! multiplicity of  midrapidity charged tracks for HM events 
   TH2D* fhMeanTrackPtMB;                               //! mean track pT for MB events 
   TH2D* fhMeanTrackPtHM;                               //! mean track pT for HM events 

   //hadron TT
   TH2D* fhTTHinMB_V0M[fkTTbins];                        //! counter of semi-inclusive hadron TT  in MB versus V0M   
   TH2D* fhTTHinMB_CentV0M[fkTTbins];                    //! counter of semi-inclusive hadron TT  in MB versus V0M    centrality
   TH2D* fhTTHinMB_V0Mnorm1[fkTTbins];                   //! counter of semi-inclusive hadron TT  in MB versus V0M/mean 
   TH2D* fhTTHinMB_V0Mnorm2[fkTTbins];                   //! counter of semi-inclusive hadron TT  in MB versus (V0A/mean + V0C/mean)/2

   TH2D* fhTTHinMB_V0M_PartLevel[fkTTbins];              //! counter of semi-inclusive hadron TT  in MB versus V0M  particle level 
   TH2D* fhTTHinMB_V0Mnorm1_PartLevel[fkTTbins];         //! counter of semi-inclusive hadron TT  in MB versus V0M/mean particle level
   TH2D* fhTTHinMB_V0Mnorm2_PartLevel[fkTTbins];         //! counter of semi-inclusive hadron TT  in MB versus V0A/<V0A>+V0C/<V0C> particle level

   TH2D* fhTTHinHM_V0M[fkTTbins];                        //! counter of semi-inclusive hadron TT in HM versus V0M   
   TH2D* fhTTHinHM_CentV0M[fkTTbins];                    //! counter of semi-inclusive hadron TT in HM versus V0M centrality   
   TH2D* fhTTHinHM_V0Mnorm1[fkTTbins];                   //! counter of semi-inclusive hadron TT in HM versus V0M/mean V0M
   TH2D* fhTTHinHM_V0Mnorm2[fkTTbins];                   //! counter of semi-inclusive hadron TT in HM versus (V0A/mean +V0C/mean)/2

   //EMCAL cluster TT
   TH2D* fhTTCinMB_V0M[fkTTbins];                        //! counter of semi-inclusive emcal TT in MB versus V0M   
   TH2D* fhTTCinMB_CentV0M[fkTTbins];                    //! counter of semi-inclusive emcal TT in MB versus V0M  centrality 
   TH2D* fhTTCinMB_V0Mnorm1[fkTTbins];                   //! counter of semi-inclusive emcal TT in MB versus V0M/mean
   TH2D* fhTTCinMB_V0Mnorm2[fkTTbins];                   //! counter of semi-inclusive emcal TT in MB versus (V0A/mean+V0C/mean)/2

   TH2D* fhTTCinMB_V0M_PartLevel[fkTTbins];              //! counter of semi-inclusive emcal TT in MB versus V0M    particle level
   TH2D* fhTTCinMB_V0Mnorm1_PartLevel[fkTTbins];         //! counter of semi-inclusive emcal TT in MB versus V0M/mean particle level
   TH2D* fhTTCinMB_V0Mnorm2_PartLevel[fkTTbins];         //! counter of semi-inclusive emcal TT in MB versus (V0A/mean +V0C/mean)/2 particle level

   TH2D* fhTTCinHM_V0M[fkTTbins];                        //! counter of semi-inclusive emcal TT in HM versus V0M   
   TH2D* fhTTCinHM_CentV0M[fkTTbins];                    //! counter of semi-inclusive emcal TT in HM versus V0M  centrality  
   TH2D* fhTTCinHM_V0Mnorm1[fkTTbins];                   //! counter of semi-inclusive emcal TT in HM versus V0M/mean
   TH2D* fhTTCinHM_V0Mnorm2[fkTTbins];                   //! counter of semi-inclusive emcal TT in HM versus (V0A/mean +V0C/mean)/2

   TH2D* fhTTCinGA_V0M[fkTTbins];                        //! counter of semi-inclusive emcal TT in GA versus V0M   
   TH2D* fhTTCinGA_CentV0M[fkTTbins];                    //! counter of semi-inclusive emcal TT in GA versus V0M centrality   
   TH2D* fhTTCinGA_V0Mnorm1[fkTTbins];                   //! counter of semi-inclusive emcal TT in GA versus V0M/mean
   TH2D* fhTTCinGA_V0Mnorm2[fkTTbins];                   //! counter of semi-inclusive emcal TT in GA versus (V0A/mean+V0C/mean)/2

   //recoil jet yields with hadron TT
   TH2D* fhRecoilJetPtTTHinMB_V0M[fkTTbins];             //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus V0M   
   TH2D* fhRecoilJetPtTTHinMB_CentV0M[fkTTbins];         //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus V0M centrality  
   TH2D* fhRecoilJetPtTTHinMB_V0Mnorm1[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus V0M/mean
   TH2D* fhRecoilJetPtTTHinMB_V0Mnorm2[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus (V0A/mean+V0C/mean)/2

   TH3D* fhRecoilJetPhiTTHinMB_V0Mnorm1[fkTTbins];           //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  minimum bias 
   TH3D* fhRecoilJetPhiTTHinMB_V0Mnorm1_PartLevel[fkTTbins]; //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  minimum bias particle level 

   TH2D* fhRecoilJetPtTTHinMB_V0M_PartLevel[fkTTbins];     //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus V0M   
   TH2D* fhRecoilJetPtTTHinMB_V0Mnorm1_PartLevel[fkTTbins]; //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus V0M/mean
   TH2D* fhRecoilJetPtTTHinMB_V0Mnorm2_PartLevel[fkTTbins]; //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in MB versus (V0A/mean+V0C/mean)/2

   TH2D* fhRecoilJetPtTTHinHM_V0M[fkTTbins];            //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in HM versus V0M   
   TH2D* fhRecoilJetPtTTHinHM_CentV0M[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in HM versus V0M centrality  
   TH2D* fhRecoilJetPtTTHinHM_V0Mnorm1[fkTTbins];       //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in HM versus V0M/mean
   TH2D* fhRecoilJetPtTTHinHM_V0Mnorm2[fkTTbins];       //! pT spectrum of recoil jets associated to semi-inclusive hadron TT in HM versus (V0A/mean+V0C/mean)/2

   TH3D* fhRecoilJetPhiTTHinHM_V0Mnorm1[fkTTbins];       //! recoil jet  (V0M/mean , recoil jet pT,  delta phi)  high multiplicity

   //recoil jet yields with EMCAL cluster TT
   TH2D* fhRecoilJetPtTTCinMB_V0M[fkTTbins];             //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M   
   TH2D* fhRecoilJetPtTTCinMB_CentV0M[fkTTbins];         //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M    centrality 
   TH2D* fhRecoilJetPtTTCinMB_V0Mnorm1[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M/mean
   TH2D* fhRecoilJetPtTTCinMB_V0Mnorm2[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus (V0A/mean+V0C/mean)/2

   TH2D* fhRecoilJetPtTTCinMB_V0M_PartLevel[fkTTbins];      //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M   
   TH2D* fhRecoilJetPtTTCinMB_V0Mnorm1_PartLevel[fkTTbins]; //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0M/mean
   TH2D* fhRecoilJetPtTTCinMB_V0Mnorm2_PartLevel[fkTTbins]; //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in MB versus V0A/<V0A>+V0C/<V0C>

   TH2D* fhRecoilJetPtTTCinHM_V0M[fkTTbins];             //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in HM versus V0M   
   TH2D* fhRecoilJetPtTTCinHM_CentV0M[fkTTbins];         //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in HM versus V0M centrality  
   TH2D* fhRecoilJetPtTTCinHM_V0Mnorm1[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in HM versus V0M/mean
   TH2D* fhRecoilJetPtTTCinHM_V0Mnorm2[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in HM versus (V0A/mean+V0C/mean)/2

   TH2D* fhRecoilJetPtTTCinGA_V0M[fkTTbins];            //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in GA versus V0M   
   TH2D* fhRecoilJetPtTTCinGA_CentV0M[fkTTbins];        //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in GA versus V0M centrality  
   TH2D* fhRecoilJetPtTTCinGA_V0Mnorm1[fkTTbins];       //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in GA versus V0M/mean
   TH2D* fhRecoilJetPtTTCinGA_V0Mnorm2[fkTTbins];       //! pT spectrum of recoil jets associated to semi-inclusive emcal TT in GA versus (V0A/mean+V0C/mean)/2



   TH2D* fhDeltaPtTTHinMB_RC_CentV0M[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in MB versus V0M centrality  
   TH2D* fhDeltaPtTTHinHM_RC_CentV0M[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in HM versus V0M centrality  
   TH2D* fhDeltaPtTTCinMB_RC_CentV0M[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in MB versus V0M centrality   
   TH2D* fhDeltaPtTTCinHM_RC_CentV0M[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in HM versus V0M centrality  
   TH2D* fhDeltaPtTTCinGA_RC_CentV0M[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in GA versus V0M centrality  

   TH2D* fhDeltaPtTTHinMB_RC_V0Mnorm1[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in MB versus V0M/mean V0M  
   TH2D* fhDeltaPtTTHinHM_RC_V0Mnorm1[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in HM versus V0M/mean V0M
   TH2D* fhDeltaPtTTCinMB_RC_V0Mnorm1[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in MB versus V0M/mean V0M  
   TH2D* fhDeltaPtTTCinHM_RC_V0Mnorm1[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in HM versus V0M/mean V0M 
   TH2D* fhDeltaPtTTCinGA_RC_V0Mnorm1[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in GA versus V0M/mean V0M 

   TH2D* fhDeltaPtTTHinMB_RC_V0Mnorm2[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in MB versus (V0A/mean + V0C/mean)/2
   TH2D* fhDeltaPtTTHinHM_RC_V0Mnorm2[fkTTbins];         //! delta pT spectrum from random cones  in events with hadron TT in HM versus (V0A/mean + V0C/mean)/2
   TH2D* fhDeltaPtTTCinMB_RC_V0Mnorm2[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in MB versus (V0A/mean + V0C/mean)/2
   TH2D* fhDeltaPtTTCinHM_RC_V0Mnorm2[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in HM versus (V0A/mean + V0C/mean)/2
   TH2D* fhDeltaPtTTCinGA_RC_V0Mnorm2[fkTTbins];         //! delta pT spectrum from random cones  in events with emcal  TT in GA versus (V0A/mean + V0C/mean)/2


   TH2D* fhDeltaPtTTHinMB_RC_V0Mnorm1_PartLevel[fkTTbins];  //! delta pT spectrum from random cones  in events with  hadron TT in MB versus V0M/mean V0M   PARTICLE LEVEL
   TH2D* fhDeltaPtTTHinMB_RC_V0Mnorm2_PartLevel[fkTTbins];  //! delta pT spectrum from random cones  in events with  hadron TT in MB versus (V0A/mean + V0C/mean)/2 P.L.
   TH2D* fhDeltaPtTTCinMB_RC_V0Mnorm1_PartLevel[fkTTbins];  //! delta pT spectrum from random cones  in events with  emcal  TT in MB versus V0M/mean V0M   PARTICLE LEVEL
   TH2D* fhDeltaPtTTCinMB_RC_V0Mnorm2_PartLevel[fkTTbins];  //! delta pT spectrum from random cones  in events with  emcal  TT in MB versus (V0A/mean + V0C/mean)/2 P.L.


   TH2D* fhPtTrkTruePrimGen;                            //! physical primary mc particle eta vs pT
   TH2D* fhPtTrkTruePrimRec;                            //! physical primary detector level track eta vs pT
   TH2D* fhPtTrkSecOrFakeRec;                           //! secondary tracks eta vs pT

   TH1D* fhJetPtPartLevelCorr;                          //! response matrix normalization spectrum, jet pT corrected on rho
   TH1D* fhJetPtPartLevelZero;                          //! response matrix normalization spectrum, jet pT is not corrected on rho
   TH1D* fhJetPtPartLevelCorrTTHpl[fkTTbins];           //! response matrix normalization spectrum, events with part. level TTH 
   TH1D* fhJetPtPartLevelCorrTTHdl[fkTTbins];           //! response matrix normalization spectrum, events with det. level TTH 

   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorr;           //! response matrix jet pT corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelZero;           //! response matrix jet pT not corrected on rho
   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorrTTHpl[fkTTbins];  //! response matrix events with particle level TTH
   TH2D* fhJetPtPartLevelVsJetPtDetLevelCorrTTHdl[fkTTbins];  //! response matrix events with detector level TTH

   TH2D* fhJetPtResolutionVsPtPartLevel;                //! resolution of jet pT

   TH2D* fhOneOverPtVsPhiNeg;                           //! 1/p_T,track  versus phi for negative tracks //AID//             
   TH2D* fhOneOverPtVsPhiPos;                           //! 1/p_T,track  versus phi for positive tracks //AID//             
   TH2D* fhSigmaPtOverPtVsPt;                           //! resolution of 1/p_T,track  versus p_T,track //AID//             
   TH2D* fhDCAinXVsPt;                                  //! X DCA versus pT  //AID// 
   TH2D* fhDCAinYVsPt;                                  //! Y DCA versus pT  //AID//
   TH2D* fhDCAinXVsPtPhysPrimary;                       //! X DCA versus pT for physical primaries  //AID// 
   TH2D* fhDCAinYVsPtPhysPrimary;                       //! Y DCA versus pT for physical primaries  //AID//
   TH2D* fhDCAinXVsPtSecondary;                         //! X DCA versus pT for secondaries //AID// 
   TH2D* fhDCAinYVsPtSecondary;                         //! Y DCA versus pT for secondaries //AID//
   TH2D* fhFractionOfSecInJet;                          //! Fraction of jet pT carried by secondaries //AID//

   TH2D* fhV0ARunByRunMB;                               //! run by run V0M
   TH2D* fhV0CRunByRunMB;                               //! run by run V0M
   TH2D* fhV0MRunByRunMB;                               //! run by run V0M
   TH2D* fhV0MnormRunByRunMB;                           //! run by run V0M norm

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
                                                      
   Bool_t   fFillTTree;                               // Fill output TTree
   Int_t    fSystem;                                  // Collision system 
   AliEMCALRecoUtils          *fFiducialCellCut;      //!

   Double_t fMeanV0A[2000];                           //! mean V0A signal in incl. MB  run by run 
   Double_t fMeanV0C[2000];                           //! mean V0C signal in incl. MB  run by run
   Double_t fMeanV0M[2000];                           //! mean V0C signal in incl. MB  run by run
   Int_t    fRuns[2000];                              //! run numbers
   Int_t    fnRun;                                    //!  the number of runs 
   Double_t fMeanV0A_PartLevel;                       // mean V0A signal in incl. MB particle level 
   Double_t fMeanV0C_PartLevel;                       // mean V0C signal in incl. MB particle level 
   Double_t fMeanV0M_PartLevel;                       // mean V0M signal in incl. MB particle level 
 


   Int_t fIndexTTC[fkTTbins];                      //! index of the chosen EMCAL cluster trigger 
   Int_t fIndexTTH[fkTTbins];                      //! index of the chosen hadron trigger 
   Int_t fIndexTTJ[fkTTbins];                      //! index of the chosen jet trigger 

   Int_t fdeltapT[fkTTbins];                       //! delta pT detector level
   Int_t fdeltapT_PartLevel[fkTTbins];             //! delta pT particle level

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

   AliAnalysisTaskEA(const AliAnalysisTaskEA&);
   AliAnalysisTaskEA& operator=(const AliAnalysisTaskEA&);

   ClassDef(AliAnalysisTaskEA, 14); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
#endif
