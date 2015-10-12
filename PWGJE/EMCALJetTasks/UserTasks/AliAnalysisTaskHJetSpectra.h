#ifndef ALIANALYSISTASKHJETSPECTRA_H
#define ALIANALYSISTASKHJETSPECTRA_H


class TH1I;
class TH1F;
class TH2F;
class TH2D;
class TH1D;
class TArrayD;
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

#include <vector>
#include "AliAnalysisTaskEmcalJet.h"

using std::vector;
// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (7.Oct. 2015)

class AliAnalysisTaskHJetSpectra : public AliAnalysisTaskEmcalJet {
   public:

  enum MyContainer {
     kContainerOne = 0, //analyze real data 
     kContainerTwo = 1, //analyze monte carlo
     kContainerThree,   //KT real 
     kContainerFour    //KT MC
  };

  enum MyDataType {
    kReal   = 0,  // reconstructed real data 
    kPythia = 1,  // pythia simulation 
    kHijing = 2   // hijing simulation
  };

  enum MyAnalType {
    kRec    = 0,  // reconstructed real data 
    kEff    = 1,  // MC true+recontructed 
    kEmb    = 2,  // embedding pythia jet
    kEmbSingl = 3,// embedding single track
    kKine   = 4   // kine 
  };

  enum MyRho {
    kConeRho=0, 
    kCMSRho, 
    kKtRho,
    kZeroRho,   //WITHOUT UE SUBTRACTION 
    kRho
  };

  enum MySystem {  //collision system
    kpp    = 0,  
    kpPb   = 1,  
    kPbPb  = 2   
  };

   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskHJetSpectra();
   AliAnalysisTaskHJetSpectra(const char *name);
   virtual ~AliAnalysisTaskHJetSpectra();
   void     UserCreateOutputObjects();
   void     Terminate(Option_t *);

   
  // ######### SETTERS/GETTERS
  void        SetAnalysisType(Int_t sys, Int_t typeOfData, Int_t typeOfAnal){
                 fTypeOfData = typeOfData;
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

  void        SetCentralityType(const char* type, Float_t cl=0., Float_t cm=100.){
                  fCentralityType = type;    
                  fCentPercMin    = (Double_t) cl;    
                  fCentPercMax    = (Double_t) cm;    
              } 

  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void        SetExternalRhoTaskName(const char* name) {fRhoTaskName = name;}
  void        SetExternalRhoTaskNameMC(const char* name) {fRhoTaskNameMC = name;}

  void SetTT(Double_t ttlow, Double_t tthigh){ fTTlow = ttlow; fTThigh = tthigh; } 
  void SetTTType(Int_t tttype){ fTTtype = tttype;} 
  void SetDphi(Double_t dphi){ fDphiCut = TMath::Pi() - dphi;} 
  void SetDoubleBinPrecision(Bool_t db){ fUseDoubleBinPrecision = db;} 


  void SetNofRandomCones(Int_t nrc){ fNofRandomCones = nrc;}
  void SetMinFractionShared(Double_t f)  { fMinFractionShared = f; }

  Bool_t   RetrieveEventObjects();
  Bool_t   Run();
  Bool_t   FillHistograms();

 private:

  // ######### MAIN CALCULATION FUNCTIONS
  void    GetDeltaPt(Int_t nrho,  TArrayD &rho, Double_t *dpt, 
                     Double_t ttPhi, Double_t ttEta, TClonesArray *trkArray, Bool_t isGen,
                     Double_t leadingJetExclusionProbability = 0); 
                     


  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius, TClonesArray *trkArray, Bool_t isGen); 
  //Double_t    GetPtHard();             
  Double_t    GetImpactParameter();   
  Double_t    GetSimPrimaryVertex(); 


  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isGen=0);  
  Bool_t      IsEventInAcceptance(AliVEvent* event);     
  Bool_t      IsMCEventInAcceptance(AliVEvent* event);   
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet, Bool_t suppressGhost=1); 
  

   Double_t RelativePhi(Double_t mphi,Double_t vphi); 
   Double_t EstimateBgCone(AliJetContainer *jetCont, TClonesArray *trkArray, AliVParticle* triggerHadron, Bool_t isGen=kFALSE);   
   Double_t EstimateBgKT(AliJetContainer *jetCont, TClonesArray *trkArray, AliVParticle* trackTT);  // median p/A of kt jets
   Double_t EstimateBgKTcms(AliJetContainer *jetCont, TClonesArray *trkArray, AliVParticle* triggerHadron); //CMS background
   


  Double_t GetNcoll(Double_t centr);  //gen Ncoll for given centrality
  Double_t GetDeltaR(Double_t phi1, Double_t phi2, Double_t eta1, Double_t eta2);
  Double_t GetFractionSharedPt(AliEmcalJet *jRec, AliJetContainer *jconRec, AliEmcalJet *jGen, AliJetContainer *jconGen);

  // ######### STANDARD FUNCTIONS
  void      ExecOnceLocal();                    

  // ########## USAGE TRIGGERS 
  Int_t               fCollisionSystem;      // collision system MySystem
  Int_t               fTypeOfData;           //kind of input data   MyDataType
  Int_t               fTypeOfAnal;           //kind of analysis MyAnalType

  Bool_t              fUseDefaultVertexCut;   // trigger if automatic vertex cut from helper class should be done 
  Bool_t              fUsePileUpCut;          // trigger if pileup cut should be done
  

  // ########## SOURCE INFORMATION
  TString             fRhoTaskName;           // name of rho CMS bg task for this analysis
  TString             fRhoTaskNameMC;         // MC name of rho CMS bg task for this analysis
  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Double_t            fSignalJetRadiusSquared;       // Radius for the signal jets
  // ########## CUTS 
  Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets 
  Double_t            fTrackEtaWindow;        //gc +- window in eta for tracks  
  Double_t            fMinTrackPt;            //gc Min track pt to be accepted  
  Double_t            fMinJetArea;            // Min jet area to be accepted
  TString             fCentralityType;        //gc Used centrality estimate (V0A, V0C, V0M, ...) 
  Double_t            fCentPercMin;           //centrality range lower cut    
  Double_t            fCentPercMax;           //centrality range upper cut
  Double_t            fMinFractionShared;     //Minimal fraction shared by embedded and rec jet

  // ########## EVENT PROPERTIES
  Double_t            fCrossSection;          //! gc value is filled, if pythia header is accessible 
  Double_t            fTrials;                //! gc value is filled, if pythia header is accessible 
  Double_t            fImpParam;              //! impact parameter from hijing

  // ########## GENERAL ////VARS
  TRandom3*           fRandom;                //! A random number
  AliAnalysisUtils*   fHelperClass;           //! gc Vertex selection helper
  Bool_t              fInitializedLocal;           //! gc trigger if tracks/jets are loaded  initiates calling   ExecOnce 


  Double_t            fTTlow;  //gc trigger particles TT bin lower boundary
  Double_t            fTThigh; //gc trigger particles TT bin upper boundary
  Int_t               fTTtype; //trigger particle type 0=single inclusive, 2 = inclusive  
  Double_t            fDphiCut; //minimal azimuthal angle between trigger and assoc jet 
  Bool_t              fUseDoubleBinPrecision; //use double bin precision

   TH1I               *fHistEvtSelection;     //! gc event statistics
   TH1D               *fh1Ntriggers;  //! trigger counter
   TH1D               *fh1TriggerMult; //! tirgger multiplicity in event
   TH1D               *fh1NtriggersGen;  //! trigger counter
   TH1D               *fh1TriggerMultGen;  //! trigger multiplicity in event
   THnSparse          *fHJetSpec[kRho];//!  TT associated spectrum of jets
   THnSparse          *fHJetSpecGen[kRho];//!TT associated spectrum of jets

   TH1F    *fhRhoTT[kRho-1]; //! gc X=rho from perp cone, Y=centrality
   TH1F    *fhRhoIncl[kRho-1]; //! gc X=rho from perp cone, Y=centrality
 
   TH1F    *fARhoTT[kRho-1]; //! jet area times rho from perp cone
   TH1F    *fARhoTTGen[kRho-1]; //! #### jet area times rho from perp cone

   TH1D    *fhDeltaPt[kRho-1]; //!  delta pT 
   TH1D    *fhDeltaPtEmb[kRho-1]; //! embedded delta pT 
   TH2D    *fhDeltaPtEmb2D[kRho-1]; //! embedded delta pT versus pT of the embedded jet 
   TH1D    *fhDeltaPtIncl[kRho-1]; //!  delta pT from RndCone using rho from perp cone inclusive event

   TH2F    *fhKTAreaPt;//!KT jets area versus PT

   TH2F     *fhJetPhi;   //! gc jet phi vs jet pT
   TH2F     *fhJetPhiGen;   //! gc jet phi vs jet pT
   TH2F     *fhTrackPhi; //! gc track phi vs track pT
   TH2F     *fhJetEta;   //! jet eta vs jet pT 
   TH2F     *fhJetEtaGen;   //! jet eta vs jet pT 
   TH2F     *fhTrackEta; //! track eta vs track pT
   TH1F     *fhTrackPt; //! gc X=centrality; Y= track pT
   TH1F     *fhTrackPtGen; //!   gc X=centrality; Y= track pT
   TH1F     *fhVertexZ;  //! gc vertexZ inclusive
   TH1F     *fhVertexZAccept; //! gc vertexZ accepted after vtx cut
   TH1F     *fhVertexZMC;  //! gc vertexZ inclusive in MC
   TH1F     *fhVertexZAcceptMC; //! gc vertexZ accepted after vtx cut in MC
   TH2F     *fhDphiTriggerJet[kRho]; //! gc Delta phi versus jet pT
   TH2F     *fhDphiTriggerJetGen[kRho]; //! gc Delta phi versus jet pT
   TH1F     *fhDphiTriggerJetAccept; //!Dphi of accepted jets after dphi cut

   TH1F     *fhCentrality;     //! centrality 
   TH1F     *fhCentralityV0M;  //! centrality V0 multiplicity A+C
   TH1F     *fhCentralityV0A;  //! centrality from V0A
   TH1F     *fhCentralityV0C;  //! centrality from V0C
   TH1F     *fhCentralityZNA;  //! centrality from ZNA


   //TProfile*     fh1Xsec;   //! gc pythia cross section and trials
   //TH1F*         fh1Trials; //! gc trials are added
   //TH1F*         fh1PtHard;  //! Pt har of the event...      
   TH1D*         fhImpactParameter; //! impact parameter distribution hijing
   TH1D*         fhImpactParameterTT; //! impact parameter distribution hijing versus TT

   TH1D*  fhJetPtGen[kRho];
   TH2D*  fhJetPtGenVsJetPtRec[kRho];
   TH2D*  fhJetPtResolutionVsPtGen[kRho];
   TH2D*  fhPtTrkTruePrimRec; // pt spectrum of true reconstructed primary tracks    
   TH2D*  fhPtTrkTruePrimGen; // pt spectrum of true generated primary track    
   TH2D*  fhPtTrkSecOrFakeRec; // pt spectrum of reconstructed fake or secondary tracks    

 
   TArrayD  fRhoRec;   // labels of particles on reconstructed track level
   TArrayD  fRhoMC;   // labels of particles on reconstructed track level


   Int_t  fNofRandomCones; // the number of random cones per event
 
   Double_t fZVertexCut; // vertex cut in z 
   Double_t fCutPhi;     // azimuthal cat around TT  to exclude TTjet + recoil jet in perp rho estimate

   std::vector<int> fTrigTracksGen; //list of trigger particle indices true MC
   std::vector<int> fTrigTracks; //list pf trigger particle indices

  AliAnalysisTaskHJetSpectra(const AliAnalysisTaskHJetSpectra&);
  AliAnalysisTaskHJetSpectra& operator=(const AliAnalysisTaskHJetSpectra&);

  ClassDef(AliAnalysisTaskHJetSpectra, 5); // Charged jet analysis for pA

};
#endif
