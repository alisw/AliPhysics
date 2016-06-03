#ifndef ALIANALYSISTASKHJETSPECTRA_H
#define ALIANALYSISTASKHJETSPECTRA_H


class TH1I;
class TH1F;
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

#include "AliAnalysisTaskEmcalJet.h"

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (5.Nov. 2015)

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

  enum MyCentBins {  //collision 
    kC0100 = 0,  
    kC020   = 1,  
    kC2050,  
    kC50100,
    kCOverflow, 
    kCAll
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

  void        SetCentralityType(const char* type){
                  fCentralityType = type;    
              } 

  void        SetVertexCut(Double_t vz){ fZVertexCut = vz; }   
  void        SetMinTrackPt(Double_t mpt){ fMinTrackPt = mpt;}

  void        SetExternalRhoTaskName(const char* name) {fRhoTaskName = name;}
  void        SetExternalRhoTaskNameMC(const char* name) {fRhoTaskNameMC = name;}

  void SetTT(Double_t ttlow, Double_t tthigh){ fTTlow = ttlow; fTThigh = tthigh; } 
  void SetTTType(Int_t tttype){ fTTtype = tttype;} 
  void SetDphi(Double_t dphi){ fDphiCut = TMath::Pi() - dphi;} 
  void SetDoubleBinPrecision(Bool_t db){ fUseDoubleBinPrecision = db;} 

  void SkipEventsWithSignalTT(Double_t ttlow, Double_t tthigh){ fTTLowToSkip = ttlow; fTTHighToSkip = tthigh; } //skip TT range

  void SetNofRandomCones(Int_t nrc){ fNofRandomCones = nrc;}
  void SetMinFractionShared(Double_t f)  { fMinFractionShared = f; }

  Bool_t   RetrieveEventObjects();
  Bool_t   Run();
  Bool_t   FillHistograms();

 private:

  // ######### MAIN CALCULATION FUNCTIONS
  void    GetDeltaPt(Int_t nrho,  TArrayD &rho, Double_t *dpt, 
                     Double_t ttPhi, Double_t ttEta, AliParticleContainer *trkArray, Bool_t isGen);
                     //Double_t leadingJetExclusionProbability = 0); 
                     


  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius, AliParticleContainer *trkArray, Bool_t isGen); 
  //Double_t    GetPtHard();             
  Double_t    GetImpactParameter();   
  Double_t    GetSimPrimaryVertex(); 


  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track, Bool_t isGen=0);  
  Bool_t      IsEventInAcceptance(AliVEvent* event);     
  Bool_t      IsMCEventInAcceptance(AliVEvent* event);   
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet, Bool_t suppressGhost=1); 
  

   Double_t RelativePhi(Double_t mphi,Double_t vphi); 
   Double_t EstimateBgCone(AliJetContainer *jetCont, AliParticleContainer *trkArray, AliVParticle* triggerHadron, Bool_t isGen=kFALSE);   
   Double_t EstimateBgKT(AliJetContainer *jetCont, AliParticleContainer  *trkArray, AliVParticle* trackTT);  // median p/A of kt jets
   Double_t EstimateBgKTcms(AliJetContainer *jetCont, AliParticleContainer *trkArray, AliVParticle* triggerHadron); //CMS background
   


  //Double_t GetNcoll(Double_t centr);  //gen Ncoll for given centrality
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

   Double_t           fTTLowToSkip;//  skip event that contains particle in the range fTTLowToSkip to fTTHighToSkip
   Double_t           fTTHighToSkip;// skip event that contains particle in the range fTTLowToSkip to fTTHighToSkip

   TH1I               *fHistEvtSelection;     //! gc event statistics
   TH1D               *fh1Ntriggers[kCAll];  //! trigger counter
   TH1D               *fh1TriggerMult[kCAll]; //! tirgger multiplicity in event
   TH1D               *fh1NtriggersGen[kCAll];  //! trigger counter
   TH1D               *fh1TriggerMultGen[kCAll];  //! trigger multiplicity in event
   TH2D               *fHJetSpec[kCAll][kRho];//!  TT associated spectrum of jets
   TH2D               *fHJetSpecGen[kCAll][kRho];//!TT associated spectrum of jets

   TH1F    *fhRhoTT[kCAll][kRho-1]; //! gc X=rho from perp cone, Y=centrality
   TH1F    *fhRhoIncl[kCAll][kRho-1]; //! gc X=rho from perp cone, Y=centrality
 
   TH1F    *fARhoTT[kCAll][kRho-1]; //! jet area times rho from perp cone
   //TH1F    *fARhoTTGen[kRho-1]; //! #### jet area times rho from perp cone

   TH1D    *fhDeltaPt[kCAll][kRho-1]; //!  delta pT 
   TH1D    *fhDeltaPtEmb[kCAll][kRho-1]; //! embedded delta pT 
   TH2D    *fhDeltaPtEmb2D[kCAll][kRho-1]; //! embedded delta pT versus pT of the embedded jet 
   TH1D    *fhDeltaPtEmbPerp[kCAll][kRho-1]; //! embedded delta pT (emb track is perp to TT)
   TH2D    *fhDeltaPtEmbPerp2D[kCAll][kRho-1]; //! embedded delta pT versus pT of the embedded jet (emb track is perp to TT)
   TH1D    *fhDeltaPtEmbBc2Bc[kCAll][kRho-1]; //! embedded delta pT (emb track is back-to-back in azimuth to TT)
   TH2D    *fhDeltaPtEmbBc2Bc2D[kCAll][kRho-1]; //! embedded delta pT versus pT of the embedded jet (emb track is backtoback in azimtuh w.r.t to TT)

   TH1D    *fhDeltaPtIncl[kCAll][kRho-1]; //!  delta pT from RndCone using rho from perp cone inclusive event

   TH2F    *fhKTAreaPt;//!KT jets area versus PT

   TH2F     *fhJetPhi[kCAll];   //! gc jet phi vs jet pT
   TH2F     *fhJetPhiGen[kCAll];   //! gc jet phi vs jet pT
   TH2F     *fhTrackPhi[kCAll]; //! gc track phi vs track pT
   TH2F     *fhJetEta[kCAll];   //! jet eta vs jet pT 
   TH2F     *fhJetEtaGen[kCAll];   //! jet eta vs jet pT 
   TH2F     *fhTrackEta[kCAll]; //! track eta vs track pT
   TH1F     *fhTrackPt[kCAll]; //! gc X=centrality; Y= track pT
   TH1F     *fhTrackPtGen[kCAll]; //!   gc X=centrality; Y= track pT
   TH1F     *fhVertexZ;  //! gc vertexZ inclusive
   TH1F     *fhVertexZAccept; //! gc vertexZ accepted after vtx cut
   TH1F     *fhVertexZMC;  //! gc vertexZ inclusive in MC
   TH1F     *fhVertexZAcceptMC; //! gc vertexZ accepted after vtx cut in MC
   TH2F     *fhDphiTriggerJet[kCAll][kRho]; //! gc Delta phi versus jet pT
   TH2F     *fhDphiTriggerJetGen[kCAll][kRho]; //! gc Delta phi versus jet pT
   TH1F     *fhDphiTriggerJetAccept; //!Dphi of accepted jets after dphi cut

   TH1F     *fhCentrality[kCAll];     //! centrality 
   TH1F     *fhCentralityV0M;  //! centrality V0 multiplicity A+C
   TH1F     *fhCentralityV0A;  //! centrality from V0A
   TH1F     *fhCentralityV0C;  //! centrality from V0C
   TH1F     *fhCentralityZNA;  //! centrality from ZNA
   TH1F     *fhCentralityV0MTT;  //! centrality V0 multiplicity A+C when TT is present
   TH1F     *fhCentralityV0ATT;  //! centrality from V0A when TT is present
   TH1F     *fhCentralityV0CTT;  //! centrality from V0C when TT is present
   TH1F     *fhCentralityZNATT;  //! centrality from ZNA when TT is present
  
   TH1F     *fhVzeroATotMult[kCAll]; //! V0A multiplicity for given V0A centrality selection
   TH1F     *fhVzeroATotMultTT[kCAll];   //! V0A multiplicity 

   TH1F     *fhZNAEnergy[kCAll]; //! ZDC A neutral energy for given V0A centrality selection
   TH1F     *fhZNAEnergyTT[kCAll];   //! ZDC A neutral energy 

   TH1D     *fhTrackMultiplicity[kCAll]; //! multiplicity of tracks
   TH1D     *fhTrackMultiplicityTT[kCAll]; //! multiplicity of tracks in event with TT track

   THnSparse  *fhZNAVzeroATrack[kCAll]; //! ZNA energy versus Vzero A mult versus track mult.
   THnSparse  *fhZNAVzeroATrackTT[kCAll]; //! ZNA energy versus Vzero mult. versus track mult. in events with TT
   

   //TProfile*     fh1Xsec;   //! gc pythia cross section and trials
   //TH1F*         fh1Trials; //! gc trials are added
   //TH1F*         fh1PtHard;  //! Pt har of the event...      
   TH1D*         fhImpactParameter[kCAll]; //! impact parameter distribution hijing
   TH1D*         fhImpactParameterTT[kCAll]; //! impact parameter distribution hijing versus TT

   TH1D*  fhJetPtGen[kCAll][kRho];
   TH2D*  fhJetPtGenVsJetPtRec[kCAll][kRho];
   TH2D*  fhJetPtResolutionVsPtGen[kCAll][kRho];
   TH2D*  fhPtTrkTruePrimRec[kCAll]; // pt spectrum of true reconstructed primary tracks    
   TH2D*  fhPtTrkTruePrimGen[kCAll]; // pt spectrum of true generated primary track    
   TH2D*  fhPtTrkSecOrFakeRec[kCAll]; // pt spectrum of reconstructed fake or secondary tracks    

 
   TArrayD  fRhoRec;   // labels of particles on reconstructed track level
   TArrayD  fRhoMC;   // labels of particles on reconstructed track level
   TArrayD  fCentralityBins; //bin boaders


   Int_t  fNofRandomCones; // the number of random cones per event
 
   Double_t fZVertexCut; // vertex cut in z 
   Double_t fCutPhi;     // azimuthal cat around TT  to exclude TTjet + recoil jet in perp rho estimate

   AliVParticle* fTrigTracksGen[999]; //list of trigger particle indices true MC
   AliVParticle* fTrigTracks[999]; //list pf trigger particle indices

   Int_t ficb[2];  //centrality bin 0=MB 1=CENT bin
   Double_t ftmpArray[2]; //tmp array
   Double_t ftmpArrayX[3]; //tmp array
   TArrayF fpyVtx;   //primaru vertex
   Double_t frhovec[999]; //auxiliary array to store pT/A of kT jets

   AliAnalysisTaskHJetSpectra(const AliAnalysisTaskHJetSpectra&);
   AliAnalysisTaskHJetSpectra& operator=(const AliAnalysisTaskHJetSpectra&);

   ClassDef(AliAnalysisTaskHJetSpectra, 16); // Charged jet analysis for pA

};
#endif
