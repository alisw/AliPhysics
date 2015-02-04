//-*- Mode: C++ -*-
// $Id$

#ifndef ALIANALYSISTASKCORRELATION3P_H
#define ALIANALYSISTASKCORRELATION3P_H

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

class AliAnalysisTaskCorrelation3p : public AliAnalysisTaskSE {
  public:
  /// default constructor
  AliAnalysisTaskCorrelation3p();
  /// constructor with options
  AliAnalysisTaskCorrelation3p(const char *name,const char* opt);
  /// destructor
  virtual ~AliAnalysisTaskCorrelation3p();
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
  enum Period {P10b,P10c,P10d,P10e,P10h,P11a,P11h, Nperiods = P10h};
  void SetPeriod(Period period){fperiod = period;}
  enum Trigger { pi0, tracks , pi0MC, tracksMC};
  void SetTrigger(Trigger trigger){ftrigger = trigger;}
  void SetCentralityEstimator(const char * centr) {fCentralityEstimator = centr;}
  void SetMaxNumberOfTracks(Double_t MaxN){fMaxNumberOfTracksInPPConsidered = MaxN;}
  void SetMixingScheme(Int_t MaxNEventMix,Int_t MinNofTracksMix,TArrayD MBinEdges,TArrayD ZBinEdges);
  //General cut variables:
  void SetMaxVertex(Double_t Vertex){fMaxVz = Vertex;}
  void SetAcceptanceCut(float Acceptance){fAcceptancecut=Acceptance;}
  void SetMinTriggerPt(Double_t MinTriggerPt){fMinTriggerPt = MinTriggerPt;}
  void SetMaxTriggerPt(Double_t MaxTriggerPt){fMaxTriggerPt = MaxTriggerPt;}
  void SetMinAssociatedPt(Double_t MinAssociatedPt){fMinAssociatedPt = MinAssociatedPt;}
  void SetMaxAssociatedPt(Double_t MaxAssociatedPt){fMaxAssociatedPt = MaxAssociatedPt;}
  void SetNClustersTPC(Int_t NClusters){fMinNClustersTPC=NClusters;}
  //Cuts for clusters and pi0s.
  void EnableTOFCut(Bool_t enable = kTRUE, Double_t TOFCut = 100.e-9){fTOFCutEnabled=enable; fTOFCut=TOFCut;}
  void SetfMinClusterEnergy(Double_t MinClusterEnergy = 0.3){fMinClusterEnergy = MinClusterEnergy;}
  void SetfMinBCDistance(Double_t MinBCDistance = 0.0){fMinBCDistance = MinBCDistance;}
  void SetMinNCells(Int_t MinNCells = 3){fMinNCells = MinNCells;}
  void SetfMinM02(Double_t MinM02 = 0.2){fMinM02 = MinM02;}
  void SetMassWindow(Double_t massMean = 0.135, Double_t massSigma = 0.01) { fMassInvMean = massMean; fMassInvSigma = massSigma; }
  void Setphospi0s(bool phos) {fphospions = phos;}
  void Setemcalpi0s(bool emcal) {femcalpions = emcal;}
  void SetGenerate(){fgenerate=kTRUE;}
  void Askforgensettings();
  TF1* pTdistribution(TH1D* hist, const char* name);
  
 protected:
  TList*  fOutput;                  //! list send on output slot 1
  TText*  fTextBox;                 //! text information
  TString fOption;                 // option string  

 private:
  AliAnalysisTaskCorrelation3p(const AliAnalysisTaskCorrelation3p&);
  AliAnalysisTaskCorrelation3p& operator=(const AliAnalysisTaskCorrelation3p&);
  //initializers used:
  int 			DefineSlots();
  void 		    	UsePeriod();
  void 		    	InitializeQAhistograms();
  void 		    	MakeRunNumbers();
  void 			FillHistogram(const char * key,Double_t x);
  void 			FillHistogram(const char * key,Double_t x,Double_t y);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z);
  Bool_t 		SelectEvent();
  void 			execgenerate();
  //Functions to make the array of associated and triggers:
  Int_t 	    	GetTracks(TObjArray *allrelevantParticles, AliVEvent* pEvent);
  void 		    	GetPi0s(TObjArray* allrelevantParticles, AliVEvent* pEvent);
  void 			GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);
//   void 			GetMCArray();
  void 			GetCentralityAndVertex();
  Bool_t	    	GoodCluster(AliVCluster *clu);
  Bool_t	    	IsSelected(AliVParticle * p);
  Bool_t 	    	IsSelectedTrigger(AliVParticle * p);
  Bool_t 	    	IsSelectedAssociated(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackAOD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackESD(AliVParticle* p);


  enum CollisionType{pp,PbPb,pPb};
  //Event characteristics:
  AliCentrality*    fCentrality;   //! centrality object
  const AliVVertex* fVertexobj; //! Vertex object
  Int_t 	    fRun; //The run we are in.
  Double_t 	    fVertex[3];//vertex
  Period 	    fperiod;
  CollisionType     fCollisionType;
  Bool_t 	    fisESD;
  Bool_t 	    fisAOD;
  Bool_t 	    fgenerate;//if true, no event is opened and the particles are created on the fly.
  TRandom3 *	    fRandom;
  TClonesArray*     fMcArray;
  //Objects that contain needed/used objects for the task:
  AliESDtrackCuts*  fTrackCuts;    //! track cuts object
  TObject*          fCorrelator;   //! correlation steering class
  Int_t*	    fRunNumberList;  //! List containing the lists over runnumbers.
  Int_t 	    fNruns; // Number of runs in the choosen period
  Double_t 	    fRunFillValue; //The bin of the run we are in.
  TArrayD 	    fMBinEdges; //Contains bin edges in centrality.
  TArrayD 	    fZBinEdges; //Edges for vZ binning.
  Int_t 	    fMaxNEventMix; //Maximum number of events per bin for event mixing.
  Int_t		    fMinNofTracksMix;//Minimum number of mixing tracks.
  TString 	    fCentralityEstimator; //! Centrality estimator ("V0M", "ZNA")
  Trigger 	    ftrigger;
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
  Double_t 	    fMinTriggerPt;
  Double_t 	    fMaxTriggerPt;
  Double_t 	    fMinAssociatedPt;
  Double_t 	    fMaxAssociatedPt;
  Int_t 	    fMinNClustersTPC;
  ////clusters:
  Double_t 	    fMinClusterEnergy;
  Double_t 	    fMinBCDistance;  //distance to nearest bad channel
  Int_t    	    fMinNCells;
  Double_t 	    fMinM02;
  Bool_t   	    fTOFCutEnabled;
  Double_t 	    fTOFCut;
  Double_t	    fMassInvMean ;
  Double_t	    fMassInvSigma ;
  bool    	    fphospions;
  bool    	    femcalpions;
  //Objects to adress the correct run number lists.
  static const Int_t fNRunsP10b = 57;
  static const Int_t fNRunsP10c = 46;
  static const Int_t fNRunsP10d = 62;
  static const Int_t fNRunsP10e = 126;
  static const Int_t fNRunsP10h = 93;
  static const Int_t fNRunsP11a = 58;
  static const Int_t fNRunsP11h = 108;
  //Class definition.
  ClassDef(AliAnalysisTaskCorrelation3p, 1);
};

#endif