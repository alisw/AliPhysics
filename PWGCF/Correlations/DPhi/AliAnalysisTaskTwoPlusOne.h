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
  
  void   SetEventSelectionBit( UInt_t val )        { fSelectBit = val;  }
  
  void   SetCentralityMethod(const char* method) { fCentralityMethod = method; }

  void   SetCustomBinning(const char* binningStr) { fCustomBinning = binningStr; }
  
  void   SetAlpha(Double_t val){fAlpha = val; }

 private:
  void            AddSettingsTree();                                  // add list of settings to output list

  TObjArray* CloneAndReduceTrackList(TObjArray* tracks);

  //general configuration
  Int_t  		fMixingTracks;		// size of track buffer for event mixing

  // Pointers to external UE classes
  AliAnalyseLeadingTrackUE*     fAnalyseUE;      //! points to class containing common analysis algorithms
  AliTwoPlusOneContainer* fHistos;// Histogram class based on UEHist
  
  // Handlers and events
  AliAODEvent*             fAOD;             //! AOD Event 
  AliEventPoolManager*     fPoolMgr;         //! event pool manager

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
    UInt_t         	fSelectBit;            // Select events according to AliAnalysisTaskJetServices bit maps 

    TString fCustomBinning;	   // supersedes default binning if set, see AliUEHist::GetBinning or AliUEHistograms::AliUEHistograms for syntax and examples
    Double_t fAlpha;            //sets the alpha parameter in the container
    
    AliAnalysisTaskTwoPlusOne(const AliAnalysisTaskTwoPlusOne&); // not implemented
    AliAnalysisTaskTwoPlusOne& operator=(const AliAnalysisTaskTwoPlusOne&); // not implemented


  ClassDef(AliAnalysisTaskTwoPlusOne, 1); // two plus one analysis with two trigger particles and particle correlations to these triggers
};

#endif
