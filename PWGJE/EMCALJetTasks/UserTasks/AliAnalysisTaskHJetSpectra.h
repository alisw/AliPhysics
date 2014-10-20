#ifndef ALIANALYSISTASKHJETSPECTRA_H
#define ALIANALYSISTASKHJETSPECTRA_H


class TH1I;
class TH1F;
class TH2F;
class TH2D;
class TH1D;
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

// ANALYSIS OF HIGH PT HADRON TRIGGER ASSOCIATED SPECTRUM OF RECOIL JETS IN P+PB
// Author Filip Krizek   (17.May. 2014)

class AliAnalysisTaskHJetSpectra : public AliAnalysisTaskSE {
   public:
   // ######### CONTRUCTORS/DESTRUCTORS AND STD FUNCTIONS
   AliAnalysisTaskHJetSpectra();
   AliAnalysisTaskHJetSpectra(const char *name, const char* trackArrayName, const char* jetArrayName, const char* backgroundJetArrayName);
   virtual ~AliAnalysisTaskHJetSpectra();
   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual Bool_t   UserNotify();
   virtual void     Terminate(Option_t *);

  // ######### SETTERS/GETTERS
  void        SetAnalyzePythia(Bool_t val) { if(val) fAnalyzePythia = kTRUE;}
  void        SetAnalyzeMC(Int_t val); 
  void        SetUseDefaultVertexCut (Bool_t val) {fUseDefaultVertexCut = val;} 
  void        SetUsePileUpCut (Bool_t val) {fUsePileUpCut = val;} 
  void        SetNumberOfCentralityBins(Int_t val) {fNumberOfCentralityBins = val;} 
  void        SetSignalJetMinArea(Double_t minArea) {fMinJetArea = minArea;}
  void        SetRandConeRadius(Double_t radius) {fRandConeRadius = radius;}
  void        SetSignalJetRadius(Double_t radius) {fSignalJetRadius = radius;}
  void        SetBackgroundJetRadius(Double_t radius) {fBackgroundJetRadius = radius;}
  void        SetMinPtOfJetsToBeRemovedInBg(Double_t minPt) {fBackgroundJetPtMin = minPt;}
  void        SetCentralityType(const char* type) {fCentralityType = type;} 
  void        SetExternalRhoTaskName(const char* name) {fRhoTaskName = name;}
  void        SetAcceptanceWindows(Double_t trackEta, Double_t signalJetRadius, Double_t bgrdJetRadius){
                 fTrackEtaWindow = trackEta; 
                 fSignalJetRadius = signalJetRadius; 
                 fBackgroundJetRadius = bgrdJetRadius; 
                 fSignalJetEtaWindow = fTrackEtaWindow-fSignalJetRadius; 
                 fBackgroundJetEtaWindow = fTrackEtaWindow-fBackgroundJetRadius;} 

  void SetTT(Double_t ttlow, Double_t tthigh){ fTTlow = ttlow; fTThigh = tthigh; }
  void SetTTType(Int_t tttype){ fTTtype = tttype;} 
  void SetDphi(Double_t dphi){ fDphiCut = TMath::Pi() - dphi;} 
  void SetDoubleBinPrecision(Bool_t db){ fUseDoubleBinPrecision = db;} 
  //void SetMC(Bool_t mc){ fIsMC = mc; };
  void SetNofRandomCones(Int_t nrc){ fNofRandomCones = nrc;}

 private:

  // ######### MAIN CALCULATION FUNCTIONS
  void    GetDeltaPt(Double_t rho1, Double_t &dpt1, 
                     Double_t rho2, Double_t &dpt2, 
                     Double_t rho3, Double_t &dpt3, 
                     Double_t &rcPhi, Double_t &rcEta,
                     Double_t leadingJetExclusionProbability = 0); 
                     


  Double_t    GetConePt(Double_t eta, Double_t phi, Double_t radius); 
  Double_t    GetPtHard();
  Double_t    GetImpactParameter();
  Double_t    GetSimPrimaryVertex();
//FK//  Double_t    GetPythiaTrials();

  void        GetPerpendicularCone(Double_t vecPhi, Double_t vecTheta, Double_t& conePt);

  // ######### CHECK FUNCTIONS
  Bool_t      IsTrackInAcceptance(AliVParticle* track); 
  Bool_t      IsEventInAcceptance(AliVEvent* event); 
  Bool_t      IsBackgroundJetInAcceptance(AliEmcalJet* jet); 
  Bool_t      IsSignalJetInAcceptance(AliEmcalJet* jet); 
  

   Double_t RelativePhi(Double_t mphi,Double_t vphi); 
   Double_t EstimateBgRhoMedian(); 
   Double_t EstimateBgCone();  
   Double_t GetExternalRho(); 

   Bool_t DistantCones(Double_t phi1, Double_t eta1, Double_t r1, Double_t phi2, Double_t eta2, Double_t r2);

  // ######### STANDARD FUNCTIONS
  void      Calculate(AliVEvent* event);   
  void      ExecOnce();                    

  TList*              fOutputList;            //! Output list
  // ########## USAGE TRIGGERS 
  Bool_t              fAnalyzePythia;         // trigger if pythia properties should be processed
  Bool_t              fAnalyzeHijing;         // trigger if pythia properties should be processed
  Bool_t              fIsKinematics;          // trigger if data is kinematics only (for naming reasons)
  Bool_t              fUseDefaultVertexCut;   // trigger if automatic vertex cut from helper class should be done 
  Bool_t              fUsePileUpCut;          // trigger if pileup cut should be done
  

  // ########## SOURCE INFORMATION
  TClonesArray*       fJetArray;              //! object containing the jets   
  TClonesArray*       fTrackArray;            //! object containing the tracks 
  TClonesArray*       fBackgroundJetArray;    //! object containing background jets
  TString*            fJetArrayName;          // name of object containing the jets
  TString*            fTrackArrayName;        // name of object containing the tracks 
  TString*            fBackgroundJetArrayName;// name of object containing event wise bckgrds
  TString             fRhoTaskName;           // name of rho CMS bg task for this analysis
  // ########## JET/DIJET/RC PROPERTIES
  Double_t            fRandConeRadius;        // Radius for the random cones
  Double_t            fRandConeRadiusSquared; // Radius for the random cones squared
  Double_t            fSignalJetRadius;       // Radius for the signal jets
  Double_t            fBackgroundJetRadius;   // Radius for the jets to be removed from bg 
  Double_t            fBackgroundJetPtMin;    // Minimum pt of jets which are ignored during bg calculation
  // ########## CUTS 
  Double_t            fSignalJetEtaWindow;    // +- window in eta for signal jets 
  Double_t            fBackgroundJetEtaWindow;// +- window in eta for background jets 
  Double_t            fTrackEtaWindow;        // +- window in eta for tracks  
  Double_t            fMinTrackPt;            // Min track pt to be accepted  
  Double_t            fMinJetArea;            // Min jet area to be accepted
  Int_t               fNumberOfCentralityBins;// Number of centrality bins used for histograms
  TString             fCentralityType;        // Used centrality estimate (V0A, V0C, V0M, ...) 

  // ########## EVENT PROPERTIES
  Double_t            fCrossSection;          //! value is filled, if pythia header is accessible 
  Double_t            fTrials;                //! value is filled, if pythia header is accessible 
  Double_t            fImpParam;              //! impact parameter from hijing

  // ########## GENERAL ////VARS
  TRandom3*           fRandom;                //! A random number
  AliAnalysisUtils*   fHelperClass;           //! Vertex selection helper
  Bool_t              fInitialized;           //! trigger if tracks/jets are loaded  initiates calling   ExecOnce 


  Double_t            fTTlow;  //trigger particles TT bin lower boundary
  Double_t            fTThigh; //trigger particles TT bin upper boundary
  Int_t               fTTtype; //trigger particle type 0=single inclusive, 2 = inclusive  
  Double_t            fDphiCut; //minimal azimuthal angle between trigger and assoc jet 
  Bool_t              fUseDoubleBinPrecision; //use double bin precision

   TH1I               *fHistEvtSelection;     //!  event statistics
   TH2F               *fh2Ntriggers;  //! trigger counter
   THnSparse         *fHJetSpec;//!  TT associated spectrum of jets
   THnSparse         *fHJetSpecSubUeMedian;//! TT associated spectrum of jets, jetPT corredted for UE cell median
   THnSparse         *fHJetSpecSubUeCone;//! TT associated spectrum of jets, jetPT corredted for UE perp cone
   THnSparse         *fHJetSpecSubUeCMS; //! TT associated spectrum of jets, jetPT corredted for UE CMS

   TH2F    *fhRhoCellMedian; //! X=rho from cell median Y=centrality
   TH2F    *fhRhoCone; //! X=rho from perp cone, Y=centrality
   TH2F    *fhRhoCMS;  //! X=rho from CMS, Y=centrality
   TH2F    *fhRhoCellMedianIncl; //! X=rho from cell median Y=centrality
   TH2F    *fhRhoConeIncl; //! X=rho from perp cone, Y=centrality
   TH2F    *fhRhoCMSIncl;  //! X=rho from CMS, Y=centrality
 
   TH1F    *fARhoCellMedian;//! jet area times rho from cell median
   TH1F    *fARhoCone; //! jet area times rho from perp cone
   TH1F    *fARhoCMS;//! jet area times rho from CMS

   TH2D    *fhDeltaPtMedian; //! delta pT from RndCone using rho from cell median high pT particle in event 
   TH2D    *fhDeltaPtCone; //! delta pT from RndCone using rho from perp cone high pT particle in event
   TH2D    *fhDeltaPtCMS; //! delta pT from RndCone using rho CMS high pT particle in event
   TH2D    *fhDeltaPtMedianIncl; //! delta pT from RndCone using rho from cell median inclusive event
   TH2D    *fhDeltaPtConeIncl; //! delta pT from RndCone using rho from perp cone inclusive event
   TH2D    *fhDeltaPtCMSIncl; //! delta pT from RndCone using rho CMS inclusive event

   TH2D    *fhDeltaPtMedianNearSide; //!  delta pt fluctuations from near side w.r.t. trigger
   TH2D    *fhDeltaPtMedianAwaySide;//! delta pt from away side
   TH2D    *fhDeltaPtCMSNearSide;//! delta pt fluctuations from near side w.r.t. trigger
   TH2D    *fhDeltaPtCMSAwaySide;//! delta pt from away side

   TH2D    *fhDeltaPtMedianExclTrigCone;//! delta pt exclude a cone around trigger
   TH2D    *fhDeltaPtCMSExclTrigCone;//!  delta pt exclude a cone around trigger

   TH2D    *fhDeltaPtMedianExclAwayJet;//! delta pt exclude a cone around leading jet on away side 
   TH2D    *fhDeltaPtCMSExclAwayJet;//!  delta pt exclude a cone around leading jet on away side



   TH2F     *fhJetPhi;   //! jet phi vs jet pT
   TH2F     *fhTrackPhi; //! track phi vs track pT
   TH2F     *fhJetEta;   //! jet eta vs jet pT 
   TH2F     *fhTrackEta; //! track eta vs track pT
   TH2F     *fhTrackCentVsPt; //!  X=centrality; Y= track pT
   TH1F     *fhVertexZ;  //! vertexZ inclusive
   TH1F     *fhVertexZAccept; //! vertexZ accepted after vtx cut
   TH2F     *fhDphiTriggerJetMinBias; //! Delta phi versus jet pT
   TH2F     *fhDphiTriggerJetCent20; //! Delta phi versus jet pT
   TH1F     *fhDphiTriggerJetAccept; //!Dphi of accepted jets after dphi cut

   TH1F     *fhCentrality;     //! centrality 
   TH1F     *fhCentralityV0M;  //! centrality V0 multiplicity A+C
   TH1F     *fhCentralityV0A;  //! centrality from V0A
   TH1F     *fhCentralityV0C;  //! centrality from V0C
   TH1F     *fhCentralityZNA;  //! centrality from ZNA

   Int_t    fNofRndTrials;     //! number of random trials for cell area estimate
   Double_t fJetFreeAreaFrac;  //! minimal fraction of cell area to be accepted to cell median
   Int_t    fnEta;    //! the number of cells in eta direction
   Int_t    fnPhi;    //! the number of cell in phi direction 
   Double_t fEtaSize;  //! size of cell in eta
   Double_t fPhiSize;  //! size of cell in phi
   Double_t fCellArea; //! cell area  

   TProfile*     fh1Xsec;   //! pythia cross section and trials
   TH1F*         fh1Trials; //! trials are added
   TH1F*         fh1PtHard;  //! Pt har of the event...      
   TH1D*         fhImpactParameter; //! impact parameter distribution hijing
   TH1D*         fhImpactParameterTT; //! impact parameter distribution hijing versus TT

   Int_t  fNofRandomCones; // the number of random cones per event
   
   Double_t fRConesR;      // small random cone of radius R=0.1 
   Double_t fRConesRSquared;      // small random cone of radius R=0.1 
   Int_t    fnRCones;      // the number of small random cones R=0.1 
   Double_t fRConePhi[50]; //! phi of small R=0.1 random cone 
   Double_t fRConeEta[50]; //! eta of small R=0.1 random cone 
 

   //Bool_t fIsMC;   

  AliAnalysisTaskHJetSpectra(const AliAnalysisTaskHJetSpectra&);
  AliAnalysisTaskHJetSpectra& operator=(const AliAnalysisTaskHJetSpectra&);

  ClassDef(AliAnalysisTaskHJetSpectra, 3); // Charged jet analysis for pA

};
#endif
