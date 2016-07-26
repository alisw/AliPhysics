#ifndef ALIANALYSISTASKTWOPLUSONE_H
#define ALIANALYSISTASKTWOPLUSONE_H

class TList;
class TH1F;
class TH2F;
class TProfile;
class THnSparse;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"
#include "AliTwoPlusOneContainer.h"
#include "TString.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "THn.h"

class AliAODEvent;
class AliAnalyseLeadingTrackUE;
class AliInputEventHandler;
class AliMCEvent;
class AliMCEventHandler;
class AliTwoPlusOneContainer;
class AliVParticle;
class TH1;
class TObjArray;
class AliEventPoolManager;
class AliESDEvent;
class AliHelperPID;
class AliAnalysisUtils;
class TFormula;
class TMap;
class AliGenEventHeader;

class AliAnalysisTaskTwoPlusOne : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskTwoPlusOne(const char *name="AliAnalysisTaskTwoPlusOne");
  virtual ~AliAnalysisTaskTwoPlusOne();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t *);

  // Setters/Getters
  // general configuration
  void   SetMixingTracks(Int_t tracks) { fMixingTracks = tracks; }
  void   SetTracksInVertex( Int_t val ){ fnTracksVertex = val; }
  void   SetfMode( Int_t val ){ fMode = val; }
  void   SetfIsNano( Int_t val ){ fIsNano = val; }
  void   SetZVertex( Double_t val )    { fZVertex = val; }
  
  // track cuts
  void   SetTrackEtaCut( Double_t val )    { fTrackEtaCut = val; }
  void   SetTrackEtaCutMin( Double_t val )    { fTrackEtaCutMin = val; }
  void   SetPtMin(Double_t val)            { fPtMin = val; }
  void   SetFilterBit( UInt_t val )        { fFilterBit = val;  }
  void   SetDCAXYCut(TFormula* value)      { fDCAXYCut = value; }
  void   SetSharedClusterCut(Float_t value) { fSharedClusterCut = value; }
  void   SetCrossedRowsCut(Int_t value)    { fCrossedRowsCut = value; }
  void   SetFoundFractionCut(Double_t value) { fFoundFractionCut = value; }
  void   SetTrackStatus(UInt_t status)     { fTrackStatus = status; }

  void   SetThreeParticleMixed(Bool_t flag) { fThreeParticleMixed = flag; }
  void   SetUseEventCombination(Bool_t flag) { fUseEventCombination = flag; }
  void   SetUsePP(Bool_t flag)             { fUsePP = flag; }
  
  void   SetCentralityMethod(const char* method) { fCentralityMethod = method; }

  void   SetCustomBinning(const char* binningStr) { fCustomBinning = binningStr; }
  void   SetUEHist_name(const char* ueHist_name) { fUEHist_name = ueHist_name; }
  
  void   SetAlpha(Double_t val){fAlpha = val; }
  void   SetUseLeadingPt(Bool_t flag) { fUseLeadingPt = flag; }
  void   SetUseAllT1(Bool_t flag) { fUseAllT1 = flag; }
  void   SetUseBackgroundSameOneSide(Bool_t flag) { fUseBackgroundSameOneSide = flag; }
  void   SetUseBackgroundSameFromMixedComb(Bool_t flag) { fUseBackgroundSameFromMixedComb  = flag; }
  void   SetUseSmallerPtAssoc(Bool_t flag) { fUseSmallerPtAssoc = flag; }
  void   SetRunCorrelations(Bool_t flag) { fRunCorrelations = flag; }
  void   SetRunIfPoolReady(Bool_t flag) { fRunIfPoolReady = flag; }
  void   SetRandomPosition(Bool_t flag) { fRandomPosition = flag; }
  void   SetSelectCentrality(Bool_t flag) { fSelectCentrality = flag; }

  void SetEfficiencyCorrection(THnF* hist) { fEfficiencyCorrection = hist; }

 private:
  void            AddSettingsTree();                                  // add list of settings to output list

  TObjArray* CloneAndReduceTrackList(TObjArray* tracks);
  void AddEventCombination(TObjArray* tracks);
  AliGenEventHeader* GetFirstHeader();

  //general configuration
  Int_t  		fMixingTracks;		// size of track buffer for event mixing
  Int_t                 fMode; // fMode = 0; data mode
                               // fMode = 1; mc mode

  Int_t                 fIsNano; //fIsNano = 0; normal mode
                                 //fIsNano = 1; use nano AODs as input

  // Pointers to external UE classes
  AliAnalyseLeadingTrackUE*     fAnalyseUE;      //! points to class containing common analysis algorithms
  AliTwoPlusOneContainer* fHistos;// Histogram class based on UEHist
  
  // Handlers and events
  AliAODEvent*             fAOD;             //! AOD Event 
  AliEventPoolManager*     fPoolMgr;         //! event pool manager
  AliMCEvent*              fMcEvent;         //! MC event
  AliInputEventHandler*    fMcHandler;       //! MCEventHandler 

  TObjArray*     fEventCombination;         //reduced tracklist which contains 4 semi central events which have the same multiplicity as 1 central event
  Int_t          fUsedEvents;               //used events in fEventCombination

  // Histogram settings
  TList*              fListOfHistos;    //  Output list of containers 

  // Event QA cuts
  Int_t               fnTracksVertex;        // QA tracks pointing to principal vertex
  Double_t            fZVertex;              // Position of Vertex in Z direction
  TString             fCentralityMethod;     // Method to determine centrality

    // Track cuts
    Double_t      	fTrackEtaCut;          // Maximum Eta cut on particles
    Double_t      	fTrackEtaCutMin;       // Minimum Eta cut on particles
    Double_t            fPtMin;                // Min pT to start correlations
    TFormula*           fDCAXYCut;             // additional pt dependent cut on DCA XY (only for AOD)
    Double_t            fSharedClusterCut;  // cut on shared clusters (only for AOD)
    Int_t		fCrossedRowsCut;   // cut on crossed rows (only for AOD)
    Double_t	 	fFoundFractionCut;     // cut on crossed rows/findable clusters (only for AOD)
    UInt_t           	fFilterBit;            // Select tracks from an specific track cut 
    UInt_t         	fTrackStatus;          // if non-0, the bits set in this variable are required for each track
    
    Bool_t              fThreeParticleMixed;   //0 use trigger from one event and mixed particles from another; 1 use trigger particles from two different events and mixed event from a third event
    Bool_t              fUseEventCombination;   //0 normal analysis run; 1 add 4 30-50% events up to the multiplicity of an 0-5% event
    Bool_t              fUsePP;                 //0 PbPb collisions; 1 pp collisions

    TString fCustomBinning;	   // supersedes default binning if set, see AliUEHist::GetBinning or AliUEHistograms::AliUEHistograms for syntax and examples
    TString fUEHist_name;	   // name of the AliUEHist in the AliTwoPlusOneContainer
    Double_t fAlpha;            //sets the alpha parameter in the container
    Bool_t fUseLeadingPt;        //decides if all particles of a cone are used as trigger particles or only the leading particles within alpha (apply this on near and away side)
    Bool_t fUseAllT1;            //decides if the near side yield is filled for all away side yields or only for the highest one
    Bool_t fUseBackgroundSameOneSide;            //decides if background Same is searched on one side or on both with half the alpha
    Bool_t fUseBackgroundSameFromMixedComb; //decides if the background same method is filled with the same event or with the mixed combinatorics events
    Bool_t fUseSmallerPtAssoc;                  //uses only pT which is smaller than the trigger pT
    Bool_t fRunCorrelations;                    //run correlation analysis, otherwise only the particle distribution is analyzed
    Bool_t fRunIfPoolReady;                    //run whole analysis only if the pool exists and it is ready
    Bool_t fRandomPosition;                    //use randomized position of the particles in the eveng
    Bool_t fSelectCentrality;                  //sort out centralities 7.5% - 30% and >50% because they are not used for the analysis anyway
    
    THnF* fEfficiencyCorrection;     // if non-0 this efficiency correction is applied on the fly to the filling for all particles. The factor is multiplicative, i.e. should contain 1/efficiency. Axes: eta, pT, centrality, z-vtx

    AliAnalysisTaskTwoPlusOne(const AliAnalysisTaskTwoPlusOne&); // not implemented
    AliAnalysisTaskTwoPlusOne& operator=(const AliAnalysisTaskTwoPlusOne&); // not implemented


    ClassDef(AliAnalysisTaskTwoPlusOne, 9); // two plus one analysis with two trigger particles and particle correlations to these triggers
};

#endif
