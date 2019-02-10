#ifndef ALIANALYSISTASKEA_H
#define ALIANALYSISTASKEA_H


class TH1I;
class TH1F;
class TF1;
class TH2F;
class TH2D;
class TH1D;
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

   enum {fkV0A, fkV0C, fkSPD, fkZNA, fkZNC, fkCE}; 

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

  void        SetMeanV0A(Double_t mva){ fMeanV0A = mva; }
  void        SetMeanV0C(Double_t mvc){ fMeanV0C = mvc; }

  Bool_t      PassedGATrigger();
  Bool_t      PassedMinBiasTrigger();
  Int_t       GetMaxDistanceFromBorder(AliVCluster* cluster);
  Bool_t      FinalClusterCuts(AliVCluster* cluster);

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
                                                      
  Float_t fCentralityV0A;                             //! Centrality from V0A
  Float_t fCentralityV0C;                             //! Centrality from V0C
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
  Float_t  fMultV0A;                                  //!  mult. V0A
  Float_t  fMultV0C;                                  //!  mult. V0C
  Float_t  fMultV0Anorm;                              //!  mult. V0A normalized by mean V0A
  Float_t  fMultV0Cnorm;                              //!  mult. V0C normalized by mean V0C
  Float_t  fMultV0AV0Cnorm;                           //!  mult. V0A+V0C normalized by means 
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


   TH1I               *fHistEvtSelection;             //! gc event statistics

   TH1F     *fhVertexZ;                               //! gc vertexZ inclusive
                                                      
   TH2F     *fhTrackPhiIncl;                          //! minimum bias phi inclusive ch hadron
   TH2F     *fhTrackEtaIncl;                          //! minimum bias eta inclusive
   TH2F     *fhJetPhiIncl;                            //! minimum bias phi inclusive ch jet
   TH2F     *fhJetEtaIncl;                            //! minimum bias eta inclusive
   TH2F     *fhClusterPhiInclMB;                        //! minimum bias phi inclusive cluster
   TH2F     *fhClusterEtaInclMB;                        //! minimum bias eta inclusive
   TH2F     *fhClusterPhiInclGA;                      //! minimum bias phi inclusive cluster
   TH2F     *fhClusterEtaInclGA;                      //! minimum bias eta inclusive
 
   TH1F     *fhRhoIncl;                               //! minimum bias rho inclusive
   TH1F     *fhRhoTTH[fkTTbins];                      //! in events MB with hadron TT
   TH1F     *fhRhoTTJ[fkTTbins];                      //! in events MB with jet TT
   TH1F     *fhRhoTTCinMB[fkTTbins];                  //! in events MB with cluster TT
   TH1F     *fhRhoTTCinGA[fkTTbins];                  //! in events GA with cluster TT

   TH1D* fhVertex[fkVtx];                             //! vertex distribution 
   TH1D* fhVertexTTH[fkVtx][fkTTbins];                //! vertex distribution in events biased with hadron TT
                                                      
   TH1D* fhCentralityMB[fkCE];                        //! estimated centrality based on  mult V0, mult VC, tracklets, znatower0, znctower0
   TH1D* fhCentralityTTH[fkCE][fkTTbins];             //! estimated centrality in events biased with hadron TT
   TH1D* fhCentralityTTJ[fkCE][fkTTbins];             //! estimated centrality in events biased with ch jet TT
   TH1D* fhCentralityTTCinMB[fkCE][fkTTbins];         //! estimated centrality in MB events biased with cluster TT
   TH1D* fhCentralityTTCinGA[fkCE][fkTTbins];         //! estimated centrality in GA events biased with cluster TT
                                                      
   TH1D* fhSignalMB[fkCE];                            //! distributions of centrality estimators:  mult V0, mult VC, tracklets, znatower0, znctower0
   TH1D* fhNormSumV0AV0CMB;                           //! distributions of (mult V0/mean V0A) + (mult VC/mean V0C)  in MB
   TH1D* fhSignalTTH[fkCE][fkTTbins];                 //! distributions of centrality estimators biased with hadron TT
   TH1D* fhSignalTTJ[fkCE][fkTTbins];                 //! distributions of centrality estimators biased with ch jet TT
   TH1D* fhSignalTTCinMB[fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in min bias 
   TH1D* fhSignalTTCinGA[fkCE][fkTTbins];             //! distributions of centrality estimators biased with cluster TT in Gamma trigger
   TH1D* fhNormSumV0AV0CTTH[fkTTbins];          //! distributions of (mult V0/mean V0A) + (mult VC/mean V0C) with hadron TT
   TH1D* fhNormSumV0AV0CTTJ[fkTTbins];          //! distributions of (mult V0/mean V0A) + (mult VC/mean V0C) with ch jet TT
   TH1D* fhNormSumV0AV0CTTCinMB[fkTTbins];      //! distributions of (mult V0/mean V0A) + (mult VC/mean V0C) with cluster TT in min bias 
   TH1D* fhNormSumV0AV0CTTCinGA[fkTTbins];      //! distributions of (mult V0/mean V0A) + (mult VC/mean V0C) with cluster TT in Gamma trigger

   TH2F* fhV0AvsV0C;                                   //! V0A vs V0C in MB 
   TH2F* fhV0AvsSPD;                                   //! V0A vs SPD in MB 
   TH2F* fhV0CvsSPD;                                   //! V0C vs SPD in MB 
   TH2F* fhV0AvsV0CTTH[fkTTbins];                      //! V0A vs V0C biased with hadron TT
   TH2F* fhV0AvsV0CTTJ[fkTTbins];                      //! V0A vs V0C biased with ch jet TT
   TH2F* fhV0AvsV0CTTCinMB[fkTTbins];                  //! V0A vs V0C biased with cluster TT in min bias 
   TH2F* fhV0AvsV0CTTCinGA[fkTTbins];                  //! V0A vs V0C biased with cluster TT in Gamma trigger
                                                      
   TH1D* fhMultTTHinMB[fkTTbins];                       //! multiplicity of hadron TT in MB event
   TH1D* fhMultTTJinMB[fkTTbins];                       //! multiplicity of charged jet TT in MB event
   TH1D* fhMultTTCinMB[fkTTbins];                       //! multiplicity of cluster TT in MB event
   TH1D* fhMultTTCinGA[fkTTbins];                       //! multiplicity of cluster TT in Gamma trigger event


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
                                                     
                                                     
   Int_t    fHadronTT[fkTTbins];                      //! array which stores the number of triggers in given event 
   Int_t    fJetChTT[fkTTbins];                       //! array which stores the number of jets in given event 
   Int_t    fClusterTT[fkTTbins];                     //! array which stores the number of jets in given event 
                                                      
   Bool_t   fFillTTree;                               // Fill output TTree
   Int_t    fSystem;                                  // Collision system 
   AliEMCALRecoUtils          *fFiducialCellCut;      //!<!

   Double_t fMeanV0A;                                 // mean V0A signal in incl. MB 
   Double_t fMeanV0C;                                 // mean V0C signal in incl. MB 


   AliAnalysisTaskEA(const AliAnalysisTaskEA&);
   AliAnalysisTaskEA& operator=(const AliAnalysisTaskEA&);

   ClassDef(AliAnalysisTaskEA, 6); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
#endif
