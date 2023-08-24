#ifndef ALIANALYSISTASKHJETSPECTRA_H
#define ALIANALYSISTASKHJETSPECTRA_H


class TH1I;
class TH1F;
class TF1;
class TH2F;
class TH2D;
class TH1D;
class TArrayD;
class TArrayF;
class TLorentzVector;
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

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (2 Jan 2021)

class AliAnalysisTaskHJetSpectra : public AliAnalysisTaskEmcalJet {
   public:

  enum MyAnalType {
    kRec    = 0,  // reconstructed real data 
    kEff    = 1   // MC true+recontructed 
  };

  enum MySystem {  //collision system
    kpp    = 0,  
    kpPb   = 1,  
    kPbPb  = 2   
  };


  enum {kRef=0, kSig=1, kTT=2}; //trigger track bins

   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskHJetSpectra();
   AliAnalysisTaskHJetSpectra(const char *name);
   virtual ~AliAnalysisTaskHJetSpectra();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

 
   static AliAnalysisTaskHJetSpectra* AddTaskHJetSpectra(
      Int_t collisionSystem        = 0,  // 0=pp,  1=pPb
      Int_t typeOfAnal             = 0,  // 0= Realdata, 1 = Eff with PYTHIA MC
      const char* jetarrayname   = "Jet_AKTChargedR040_tracks_pT0150_pt_scheme", //name of jet TClones array for detector level jets
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
      Double_t    acut               = 0.   //cut on relative jet area
   );


   
  // ######### SETTERS/GETTERS
  void        SetAnalysisType(Int_t sys, Int_t typeOfAnal){
                 fTypeOfAnal = typeOfAnal;
                 fCollisionSystem = sys;
              }

  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;}  
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;} 
  void        SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;} 


  void        SetSignalJetRadius(Double_t radius) {
                 fSignalJetRadius        = radius;
                 fSignalJetRadiusSquared = fSignalJetRadius*fSignalJetRadius;
                 fSignalJetEtaWindow     = fTrackEtaWindow-fSignalJetRadius; 
              }

  void        SetAcceptanceWindows(Double_t trackEta, Double_t signalJetRadius){ 
                 fTrackEtaWindow  = trackEta; 
                 fSignalJetRadius = signalJetRadius; 
                 fSignalJetEtaWindow = fTrackEtaWindow-fSignalJetRadius; 
                 fSignalJetRadiusSquared = fSignalJetRadius*fSignalJetRadius;
              } 


  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void SetTT(Double_t tlr, Double_t thr,Double_t tls, Double_t ths);
  void SetDphi(Double_t dphi){ fDphiCut = TMath::Pi() - dphi;} 

  void        SetTrackContainerName(const char* name){ fMyTrackContainerName = name;}
  void        SetMCParticleContainerName(const char* name){ fMyParticleContainerName = name;}
  void        SetJetContainerName(const char* name){ fMyJetContainerName = name;}
  void        SetMCPartJetContainerName(const char* name){ fMyJetParticleContainerName = name;}
  void        SetKTJetContainerName(const char* name){ fMyKTJetContainerName = name;}
  void        SetKTMCPartJetContainerName(const char* name){ fMyKTJetParticleContainerName = name;}



  Bool_t   RetrieveEventObjects();
  Bool_t   Run();
  Bool_t   FillHistograms();

 private:

  // ######### MAIN CALCULATION FUNCTIONS
  Double_t GetDeltaPt(Double_t rho, Double_t ttPhi, Double_t ttEta, AliParticleContainer *recTrkCont);
                     
  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius, AliParticleContainer *trkArray); 

  Double_t    GetSimPrimaryVertex(); 


  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isGen=0);  
  Bool_t      IsEventInAcceptance(AliVEvent* event);     
  Bool_t      IsMCEventInAcceptance(AliVEvent* event);   
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet, Bool_t suppressGhost=1); 

  Double_t EstimateBgKT(AliJetContainer *jetCont, AliParticleContainer  *trkArray, TLorentzVector myTT);  // median p/A of kt jets
   

  Double_t GetDeltaR(Double_t phi1, Double_t phi2, Double_t eta1, Double_t eta2); //angular distance between phi1,eta1 and phi2,eta2

  // ######### STANDARD FUNCTIONS
  void      ExecOnceLocal();                    

  // ########## USAGE TRIGGERS 
  Int_t               fCollisionSystem;      // collision system MySystem
  Int_t               fTypeOfAnal;           //kind of analysis MyAnalType

  Bool_t              fUseDefaultVertexCut;   // trigger if automatic vertex cut from helper class should be done 
  Bool_t              fUsePileUpCut;          // trigger if pileup cut should be done
  

  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Double_t            fSignalJetRadiusSquared;       // Radius for the signal jets
  // ########## CUTS 
  Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets 
  Double_t            fTrackEtaWindow;        //gc +- window in eta for tracks  
  Double_t            fMinTrackPt;            //gc Min track pt to be accepted  
  Double_t            fMinJetArea;            // Min jet area to be accepted

  // ########## GENERAL ////VARS
  TRandom3*           fRandom;                //! A random number
  AliAnalysisUtils*   fHelperClass;           //! gc Vertex selection helper
  Bool_t              fInitializedLocal;           //! gc trigger if tracks/jets are loaded  initiates calling   ExecOnce 


  Double_t            fTTlow[kTT];  //gc trigger particles TT bin lower boundary
  Double_t            fTThigh[kTT]; //gc trigger particles TT bin upper boundary
  Double_t            fDphiCut; //minimal azimuthal angle between trigger and assoc jet 


   TH1I               *fHistEvtSelection;        //! gc event statistics
   TH2D               *fhTTDet_V0A[kTT];             //! trigger counter at detector level vs centrality
   TH2D               *fhTTDet_ZNA[kTT];             //! trigger counter at detector level vs centrality
   TH2F               *fhTTMultDet_V0A[kTT];      //! trigger multiplicity in event with V0A 
   TH2F               *fhTTMultDet_ZNA[kTT];      //! trigger multiplicity in event with ZNA centrality
   TH1F               *fhTTPart[kTT];            //! trigger counter at particle level
   TH1F               *fhTTMultPart[kTT];        //! trigger multiplicity in event
   TH2D               *fHJetSpecDet_V0A[kTT];    //! TT associated spectrum of detector level jets with V0A centrality
   TH2D               *fHJetSpecDet_ZNA[kTT];    //! TT associated spectrum of detector level jets with ZNA centrality
   TH1D               *fHJetSpecPart[kTT];       //! TT associated spectrum of particle level jets

   Double_t  fRhoDet[kTT];   // pT density at detector level 
   TH2F    *fhRhoTT[kTT]; //! gc X=centrality, Y=rho from perp cone

   TH2D    *fhDeltaPt_V0A[kTT]; //! X=centrality, Y = delta pT 
   TH2D    *fhDeltaPt_ZNA[kTT]; //! X=centrality, Y = delta pT 

   TH2F     *fhTrackPhi; //! gc track phi vs track pT
   TH2F     *fhTrackEta; //! track eta vs track pT
   TH2F     *fhJetPhiDet;//!minimum bias phi inclusive
   TH2F     *fhJetEtaDet;//!minimum bias eta inclusive
   TH2F     *fhJetEtaPart; //! eta vs pt distribution of generator level jets
   TH1D     *fhJetPtPart; //! pt distribution of generator level jets  normalization of response matrix

   TH1F     *fhVertexZ;  //! gc vertexZ inclusive
   TH1F     *fhVertexXAccept; //! gc vertexZ accepted after vtx cut
   TH1F     *fhVertexYAccept; //! gc vertexZ accepted after vtx cut
   TH1F     *fhVertexZAccept; //! gc vertexZ accepted after vtx cut
   TH1F     *fhVertexZMC;  //! gc vertexZ inclusive in MC
   TH1F     *fhVertexZAcceptMC; //! gc vertexZ accepted after vtx cut in MC

   TH1F     *fhCentralityV0A;  //! centrality from V0A
   TH1F     *fhCentralityZNA;  //! centrality from ZNA 
   TH1F     *fhCentralityTT_V0A[kTT];  //! centrality V0 multiplicity when TT is present
   TH1F     *fhCentralityTT_ZNA[kTT];  //! centrality ZNA energy when TT is present
 
   TH2F     *fhVzeroATotMult; //! V0A multiplicity for given V0A centrality selection
   TH2F     *fhZNAEnergy; //! ZDC A neutral energy for given V0A centrality selection


   TH2D*  fhJetPtPartVsJetPtDet; //! pt jet gen level vs pT jet rec level  response matrix
   TH1D*  fhPhysPrimaryPtDet; //! pt spectrum of true reconstructed primary tracks    
   TH1D*  fhPhysPrimaryPtPart; //! pt spectrum of true generated primary track    
   TH1D*  fhPtTrkSecOrFakeDet; //! pt spectrum of reconstructed fake or secondary tracks    

   TH2D* fhJetAreaVsPt;            //dch: 2D jet area vs pT
 
   Double_t fZVertexCut; // vertex cut in z 
   Double_t fCutPhi;     // azimuthal cat around TT  to exclude TTjet + recoil jet in perp rho estimate

   vector<TLorentzVector> fTTH_Det[kTT];   //! trigger track candidates at detector level
   vector<TLorentzVector> fTTH_Part[kTT];  //! trigger track candidates at particle level


   Double_t frhovec[999]; //auxiliary array to store pT/A of kT jets

   AliMultSelection* fMultSelection;  //! object which handels centrality 

   TString fMyTrackContainerName; // name of detector level track container 
   TString fMyParticleContainerName;// name of particle level MC particle container
   TString fMyJetContainerName; // name of detector level jet container
   TString fMyJetParticleContainerName; // name of particle level MC jet container
   TString fMyKTJetContainerName; // name of KT detector level jet container
   TString fMyKTJetParticleContainerName; // name of KT particle level MC jet container

   AliAnalysisTaskHJetSpectra(const AliAnalysisTaskHJetSpectra&);
   AliAnalysisTaskHJetSpectra& operator=(const AliAnalysisTaskHJetSpectra&);

   ClassDef(AliAnalysisTaskHJetSpectra, 28); // Charged jet analysis for pA

};
#endif
