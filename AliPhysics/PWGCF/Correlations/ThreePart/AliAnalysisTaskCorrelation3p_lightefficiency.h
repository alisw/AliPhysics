//-*- Mode: C++ -*-
// $Id$

#ifndef ALIANALYSISTASKCORRELATION3P_LIGHTEFFICIENCY_H
#define ALIANALYSISTASKCORRELATION3P_LIGHTEFFICIENCY_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TArrayD.h"
#include "AliLog.h"
class AliCorrelation3p;
class AliESDtrackCuts;
class AliAODTrack;
class AliFilteredTrack;
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
  enum Period {P10b,P10c,P10d,P10e,P10h,P11a,P11h, Nperiods = P11h};
  void SetPeriod(Period period){fperiod = period;}
  void SetCollisionType(AliAnalysisTaskCorrelation3p_lightefficiency::CollisionType type){fCollisionType=type;}
  void SetRemoveSignals(Bool_t rem = kTRUE){fRemoveSignals = rem;}
  //Cuts for clusters and pi0s.
  void SetTrackCut(const char* cutmask){
    fCutMask = 0;//Defaults to GlobalHybrid tracks.
    //GlobalHybrid - 0
    if(TString(cutmask).CompareTo("GlobalHybrid")==0){
      fCutMask = 0;
      
      AliWarning("Track cut set to GlobalHybrid");
    }
    //Filter Bit 4 - 1
    if(TString(cutmask).CompareTo("BIT4")==0){
      fCutMask = 1;
      AliWarning("Track cut set to Bit4");
    }
    //Filter Bit 5 - 2
    if(TString(cutmask).CompareTo("BIT5")==0){
      fCutMask = 2;
      AliWarning("Track cut set to Bit5");
    }      
    //Filter Bit 6 - 3
    if(TString(cutmask).CompareTo("BIT6")==0){
      fCutMask = 3;
      AliWarning("Track cut set to Bit6");
    }
    //Filter Bit 5 | 6 - 4
    if(TString(cutmask).CompareTo("BIT56")==0){
      fCutMask = 4;
      AliWarning("Track cut set to Bit5|Bit6");
    }    
  }
  void SetTree(bool istree){fisTree=istree;}
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
  void 			IncrementHist(const char * key,int bin);
  void 			FillHistogram(const char * key,Double_t x,Double_t y);
  void 			FillHistogramGenPar(const char * key,Int_t x,Int_t y,Double_t z);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z,Double_t a,Double_t b);
  Bool_t 		SelectEvent();
  //Functions to make the array of associated and triggers:
  Int_t 	    	GetTracks(TObjArray *allrelevantParticles, AliVEvent* pEvent);
  void 			GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);
  void 			GetMCArray();
  void 			GetCentralityAndVertex();
  void 			GetCentralityAndVertex(AliVEvent * pevent);
  Bool_t	    	IsSelected(AliVParticle * p);
  Bool_t 	    	IsSelectedTrack(AliVParticle * p);
  Bool_t 	    	IsSelectedTrackAOD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackESD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackFiltered(AliVParticle* p);
  void 			FillFilterBit(AliFilteredTrack* p);
  Bool_t		IsMCFilteredTrack(AliVParticle* p);
  Bool_t 		IsAddedSignal(AliVParticle* p);
  void 			GeneratorStat(AliVParticle* p);
  Int_t 		FillGen(Int_t GetGenereratorIndex);
  const char* 		NameGen(Int_t GetGenereratorIndex);
  Int_t 		FillPDG(Int_t GetPDG);
  const char* 		NamePDG(Int_t GetPDG);

  //Event characteristics:
  AliCentrality*    fCentrality;   //! centrality object
  const AliVVertex* fVertexobj; //! Vertex object
  Double_t 	    fVertex[3];//vertex
  CollisionType     fCollisionType;
  Period 	    fperiod;
  Bool_t 	    fisESD;
  Bool_t 	    fisAOD;
  Bool_t 	    fRemoveSignals;
  Bool_t 	    fisTree;
  TClonesArray*     fMcArray;
  //Objects that contain needed/used objects for the task:
  TArrayD 	    fMBinEdges; //Contains bin edges in centrality.
  TArrayD 	    fZBinEdges; //Edges for vZ binning.
  Int_t 	    fMaxNEventMix; //Maximum number of events per bin for event mixing.
  Int_t		    fMinNofTracksMix;//Minimum number of mixing tracks.
  TString 	    fCentralityEstimator; //! Centrality estimator ("V0M", "ZNA")
  Double_t 	    fCentralityPercentile;	
  Double_t 	    fMultiplicity;
  Int_t*	    fRunNumberList;  //! List containing the lists over runnumbers.
  Int_t 	    fNruns; // Number of runs in the choosen period
  Int_t 	    fRun;//the run we are in
  Double_t	    fRunFillValue;//Fill value for the relevant histograms
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
  Int_t 	    fCutMask; // To correspond to different cut masks.
  //Objects to adress the correct run number lists.
  static const Int_t fNRunsP10b = 57;
  static const Int_t fNRunsP10c = 46;
  static const Int_t fNRunsP10d = 62;
  static const Int_t fNRunsP10e = 126;
  static const Int_t fNRunsP10h = 93;
  static const Int_t fNRunsP11a = 58;
  static const Int_t fNRunsP11h = 108;
  //Class definition.
  ClassDef(AliAnalysisTaskCorrelation3p_lightefficiency, 1);
};

#endif