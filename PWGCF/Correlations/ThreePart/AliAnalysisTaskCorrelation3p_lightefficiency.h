//-*- Mode: C++ -*-
// $Id$

#ifndef ALIANALYSISTASKCORRELATION3P_LIGHTEFFICIENCY_H
#define ALIANALYSISTASKCORRELATION3P_LIGHTEFFICIENCY_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TArrayD.h"
class AliCorrelation3p;
class AliESDtrackCuts;
class AliAODTrack;
class AliCentrality;
class TH1;
class TH2;
class TText;
class TList;
class TTree;
class TRandom3;
class TF1;

class AliAnalysisTaskCorrelation3p_lightefficiency : public AliAnalysisTaskSE {
  public:
  /// default constructor
  AliAnalysisTaskCorrelation3p_lightefficiency();
  /// constructor with options
  AliAnalysisTaskCorrelation3p_lightefficiency(const char *name,const char* opt);
  /// destructor
  virtual ~AliAnalysisTaskCorrelation3p_lightefficiency();
  /// inherited from AliAnalysisTask: connect tree branches at input slots
  virtual void ConnectInputData(Option_t *option="") {
    return AliAnalysisTaskSE::ConnectInputData(option);
  }
  /// inherited from AliAnalysisTaskSE: create output objects
  virtual void UserCreateOutputObjects();
  /// inherited from AliAnalysisTaskSE: event processing
  virtual void UserExec(Option_t*);
  /// inherited from AliAnalysisTask: called in SlaveTerminate phase for each task
  virtual void FinishTaskOutput();
  /// inherited from AliAnalysisTask: final step
  virtual void Terminate(Option_t*);
  //Setters for controlling the Task through an AddTaskMacro:
  void SetOption(const char* opt) { fOption = opt; }
  //General setters:
  void SetCentralityEstimator(const char * centr) {fCentralityEstimator = centr;}
  void SetMaxNumberOfTracks(Double_t MaxN){fMaxNumberOfTracksInPPConsidered = MaxN;}
  void SetMixingScheme(Int_t MaxNEventMix,Int_t MinNofTracksMix,TArrayD MBinEdges,TArrayD ZBinEdges);
  //General cut variables:
  void SetMaxVertex(Double_t Vertex){fMaxVz = Vertex;}
  void SetAcceptanceCut(float Acceptance){fAcceptancecut=Acceptance;}
  void SetMinPt(Double_t MinTriggerPt){fMinPt = MinTriggerPt;}
  void SetMaxPt(Double_t MaxTriggerPt){fMaxPt = MaxTriggerPt;}
  void SetNClustersTPC(Int_t NClusters){fMinNClustersTPC=NClusters;}
  enum CollisionType{pp,PbPb,pPb};
  void SetCollisionType(AliAnalysisTaskCorrelation3p_lightefficiency::CollisionType type){fCollisionType=type;}

  //Cuts for clusters and pi0s.
  
 protected:
  TList*  fOutput;                  //! list send on output slot 1
  TText*  fTextBox;                 //! text information
  TString fOption;                 // option string  

 private:
  AliAnalysisTaskCorrelation3p_lightefficiency(const AliAnalysisTaskCorrelation3p_lightefficiency&);
  AliAnalysisTaskCorrelation3p_lightefficiency& operator=(const AliAnalysisTaskCorrelation3p_lightefficiency&);
  //initializers used:
  int 			DefineSlots();
  void 		    	InitializeQAhistograms();
  void 			InitializeEffHistograms();
  void 			FillHistogram(const char * key,Double_t x);
  void 			FillHistogram(const char * key,Double_t x,Double_t y);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z,Double_t a,Double_t b);
  Bool_t 		SelectEvent();
  //Functions to make the array of associated and triggers:
  Int_t 	    	GetTracks(TObjArray *allrelevantParticles, AliVEvent* pEvent);
  void 			GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);
  void 			GetMCArray();
  void 			GetCentralityAndVertex();
  Bool_t	    	IsSelected(AliVParticle * p);
  Bool_t 	    	IsSelectedTrack(AliVParticle * p);
  Bool_t 	    	IsSelectedTrackAOD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackESD(AliVParticle* p);


  //Event characteristics:
  AliCentrality*    fCentrality;   //! centrality object
  const AliVVertex* fVertexobj; //! Vertex object
  Double_t 	    fVertex[3];//vertex
  CollisionType     fCollisionType;
  Bool_t 	    fisESD;
  Bool_t 	    fisAOD;
  TClonesArray*     fMcArray;
  //Objects that contain needed/used objects for the task:
  TArrayD 	    fMBinEdges; //Contains bin edges in centrality.
  TArrayD 	    fZBinEdges; //Edges for vZ binning.
  Int_t 	    fMaxNEventMix; //Maximum number of events per bin for event mixing.
  Int_t		    fMinNofTracksMix;//Minimum number of mixing tracks.
  TString 	    fCentralityEstimator; //! Centrality estimator ("V0M", "ZNA")
  Double_t 	    fCentralityPercentile;	
  Double_t 	    fMultiplicity;
  //cut variables:
  ////event:
  Double_t	    fMaxVz; //Vertex cut variable.
  Double_t 	    fMaxMult; //Is set automatically when setting the binning.
  Double_t 	    fMaxNumberOfTracksInPPConsidered; //Maximum number of tracks in an pp event used for correlations.
  Int_t 	    fNTriggers;//Contains the number of triggers in a given event.
  Int_t 	    fNAssociated;//Number of associated in the event.
  ////tracks:
  float 	    fAcceptancecut; //maximum eta
  Double_t 	    fMinPt;
  Double_t 	    fMaxPt;
  Int_t 	    fMinNClustersTPC;

  //Class definition.
  ClassDef(AliAnalysisTaskCorrelation3p_lightefficiency, 1);
};

#endif