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



   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskEA();
   AliAnalysisTaskEA(const char *name);
   virtual  ~AliAnalysisTaskEA();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   
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

  void        SetJetContainerName(const char* name){ fMyJetContainerName = name;}
  void        SetMCJetContainerName(const char* name){ fMyJetParticleContainerName = name;}
  void        SetHadronTT(Int_t tl,Int_t th){ fHadronTTLowPt[fnHadronTTBins] = tl; fHadronTTHighPt[fnHadronTTBins] = th; fnHadronTTBins++; }
  void        SetJetChTT(Int_t tl,Int_t th){  fJetChTTLowPt[fnJetChTTBins] = tl;   fJetChTTHighPt[fnJetChTTBins] = th;   fnJetChTTBins++;  }
  void        SetFillTTree(Bool_t b){ fFillTTree = b; } //fill output TTree

  Bool_t      PassedGATrigger();
  Bool_t      PassedMinBiasTrigger();
  Bool_t      PreSelection(AliVCluster* cluster);
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

  AliTrackContainer    *fTrkContainerDetLevel;        //! detector level track container
  AliParticleContainer *fParticleContainerPartLevel;  //! particle level container with particles
  AliJetContainer      *fJetContainerDetLevel;        //! detector level jet container  
  AliJetContainer      *fJetContainerPartLevel;       //! particle level jet container

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
                                                      
   TH2F     *fhTrackPhiIncl;                          //! minimum bias phi inclusive
   TH2F     *fhTrackEtaIncl;                          //! minimum bias eta inclusive
   TH2F     *fhJetPhiIncl;                            //! minimum bias phi inclusive
   TH2F     *fhJetEtaIncl;                            //! minimum bias eta inclusive
   TH1F     *fhRhoIncl;                               //! minimum bias rho inclusive

   TH1D* fhVertex[fkVtx];                             //! vertex distribution 
   TH1D* fhVertexTTH[fkVtx][fkTTbins];                //! vertex distribution in events biased with hadron TT
                                                      
   TH1D* fhCentralityMB[fkCE];                        //! estimated centrality based on  mult V0, mult VC, tracklets, znatower0, znctower0
   TH1D* fhCentralityTTH[fkCE][fkTTbins];             //! estimated centrality in events biased with hadron TT
   TH1D* fhCentralityTTJ[fkCE][fkTTbins];             //! estimated centrality in events biased with hadron TT
                                                      
   TH1D* fhSignalMB[fkCE];                            //! distributions of centrality estimators:  mult V0, mult VC, tracklets, znatower0, znctower0
   TH1D* fhSignalTTH[fkCE][fkTTbins];                 //! distributions of centrality estimators biased with hadron TT
   TH1D* fhSignalTTJ[fkCE][fkTTbins];                 //! distributions of centrality estimators biased with ch jet TT
                                                      
   TH1D* fhMultTTHinMB[fkTTbins];                       //! multiplicity of hadron TT in MB event
   TH1D* fhMultTTJinMB[fkTTbins];                       //! multiplicity of charged jet TT in MB event


   Double_t fZVertexCut;                              // vertex cut in z 
                                                     
   Int_t    fnHadronTTBins;                           // number of TT bins
   Int_t    fnJetChTTBins;                            // number of TT bins
   Int_t    fHadronTTLowPt[fkTTbins];                 // low pt TT range
   Int_t    fHadronTTHighPt[fkTTbins];                // low pt TT range
   Int_t    fJetChTTLowPt[fkTTbins];                  // low pt TT range
   Int_t    fJetChTTHighPt[fkTTbins];                 // low pt TT range
                                                     
                                                     
   Int_t    fHadronTT[fkTTbins];                      //! array which stores the number of triggers in given event 
   Int_t    fJetChTT[fkTTbins];                       //! array which stores the number of jets in given event 
                                                      
   Bool_t   fFillTTree;                               // Fill output TTree

   AliAnalysisTaskEA(const AliAnalysisTaskEA&);
   AliAnalysisTaskEA& operator=(const AliAnalysisTaskEA&);

   ClassDef(AliAnalysisTaskEA, 2); // Charged jet analysis for pAliAnalysisTaskHJetSpectra/home/fkrizek/z501.ALIC

};
#endif
