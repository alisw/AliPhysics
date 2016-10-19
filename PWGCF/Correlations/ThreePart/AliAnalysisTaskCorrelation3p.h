//-*- Mode: C++ -*-
// $Id$

#ifndef ALIANALYSISTASKCORRELATION3P_H
#define ALIANALYSISTASKCORRELATION3P_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TArrayD.h"
#include "TFile.h"
#include "THn.h"
#include "TH3D.h"
#include "TF1.h"
#include <iosfwd>

class AliCorrelation3p;
class AliESDtrackCuts;
class AliAODTrack;
class AliCentrality;
class TH1;
class TH2;
class TH3;
class TText;
class TList;
class TTree;
class TRandom3;
using namespace std;
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
  enum Period {P10b,P10c,P10d,P10e,P10h,P11a,P11h, Nperiods = P11h};
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
  void SetQA(){fQA=kTRUE;}
  void SetQAtask(bool qatask){fqatask=qatask;}
  void SetBinVer(int binver){fBinVer=binver;}
  void SetExtraMixed(bool extramixed){fExtraMixed=extramixed;}
  void SetLeading(bool leading){fLeading=leading;}
  void SetWeights(const char* file){
    TFile* wfile = TFile::Open(file,"OLD");
    if(wfile){
      fWeights = dynamic_cast<TH3D*>(wfile->Get("hnWeight")->Clone("hnWeight_task"));
      fWeights->SetTitle("Weight for low pT");
      fWeights->GetXaxis()->SetTitle("eta []");
      fWeights->GetYaxis()->SetTitle("Vertex [cm]");
      fWeights->GetZaxis()->SetTitle("pT [GeV/c]");
      fWeights->SetDirectory(0x0);
      fWeightshpt = dynamic_cast<TH2D*>(wfile->Get("hnWeight_highpt")->Clone("hnWeight_highpt_task"));
      fWeightshpt->SetTitle("Weight for high pT");
      fWeightshpt->GetXaxis()->SetTitle("eta [%]");
      fWeightshpt->GetYaxis()->SetTitle("Vertex [cm]");
      fWeightshpt->SetDirectory(0x0);
      fpTfunction= dynamic_cast<TF1*>(wfile->Get("pT_function")->Clone("pT_function_task"));
      fpTfunction->SetTitle("Function dependence of the Weight in pT.");
      fpTfunction->GetXaxis()->SetTitle("pT [GeV/c]");
    }
    else{fWeights = NULL;}
    
    wfile->Close();   
  }
  void SetTrackCut(const char* cutmask){
    fCutMask = 0;//Defaults to GlobalHybrid tracks.
    //GlobalHybrid - 0
    if(TString(cutmask).CompareTo("GlobalHybrid")==0) fCutMask = 0;
    //Filter Bit 4 - 1
    if(TString(cutmask).CompareTo("BIT4")==0) fCutMask = 1;
    //Filter Bit 5 - 2
    if(TString(cutmask).CompareTo("BIT5")==0) fCutMask = 2;
    //Filter Bit 6 - 3
    if(TString(cutmask).CompareTo("BIT6")==0) fCutMask = 3;
    //Filter Bit5 | Bit 6 - 4
    if(TString(cutmask).CompareTo("BIT5|BIT6")==0) fCutMask = 4;    
    //Filter Global & !(Bit 5 |Bit 6) - 5
    if(TString(cutmask).CompareTo("ExclusiveGlobal")==0) fCutMask = 5;    
  }
  void SetDstTree(bool tree = false){fisDstTree = tree;}
  void SetMaxTracksPerEvent(int i){fMaxTracksperEvent = i;}
//   void SetEfficiencies(){fefficiencies=kTRUE;}
  void Askforgensettings();
  void SetMoreOutputs(bool out = true){fMoreOutputs = out;}
  void SetStartEvent(int startevent){fStartAtEvent = startevent;}
  void SetNEvents(int nevents){fNEventsToProcess = nevents;}
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
  void 			InitializeEffHistograms();
  void 		    	MakeRunNumbers();
  void 			FillHistogram(const char * key,Double_t x);
  void 			FillHistogram(const char * key,Double_t x,Double_t y);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z);
  void 			FillHistogram(const char * key,Double_t x,Double_t y,Double_t z,Double_t a,Double_t b);
  Bool_t 		SelectEvent();
  void 			execgenerate();
  //Functions to make the array of associated and triggers:
  Int_t 	    	GetTracks(TObjArray *allrelevantParticles, AliVEvent* pEvent);
  void 		    	GetPi0s(TObjArray* allrelevantParticles, AliVEvent* pEvent);
  void 			GetDCA(Double_t& DCAtang, Double_t& DCAlong, AliAODTrack* AODt);
  void 			GetMCArray();
  void 			GetCentralityAndVertex();
  void 			GetCentralityAndVertex(AliVEvent * pevent);
  Bool_t	    	GoodCluster(AliVCluster *clu);
  Bool_t	    	IsSelected(AliVParticle * p);
  Bool_t 	    	IsSelectedTrigger(AliVParticle * p);
  Bool_t 	    	IsSelectedAssociated(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackAOD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackESD(AliVParticle* p);
  Bool_t 	    	IsSelectedTrackFiltered(AliVParticle* p);


  enum CollisionType{pp,PbPb,pPb};
  //Task control:
  Int_t fstarttime;
  //Event characteristics:
  AliCentrality*    fCentrality;//! centrality object
  const AliVVertex* fVertexobj;//! Vertex object
  Int_t 	    fRun;//The run we are in.
  Int_t		    fNEventsProcessed;// number of events Processed
  Int_t 	    fNEventsParsed;// number of events the job has Parsed
  Int_t 	    fNEventsToProcess;// number of events before the task stops.
  Int_t 	    fStartAtEvent;// start to do correlations after event x.
  Double_t 	    fVertex[3];//vertex
  Period 	    fperiod;
  CollisionType     fCollisionType;
  Bool_t 	    fisESD;
  Bool_t 	    fisAOD;
  Bool_t	    fisDstTree;
  Bool_t 	    fgenerate;//if true, no event is opened and the particles are created on the fly.
  Bool_t 	    fQA;//if true, correlations are not build.
  Bool_t 	    fqatask;//if true AliCorrelation3p_noQA is used.
  Bool_t	    fMoreOutputs;//If true, more outputs are given.
  Bool_t	    fExtraMixed;//If true, META META2 and METrigger are explicitly correlated.
  Bool_t	    fLeading;//If true only the leading pT track is considered for trigger.
  TH3D *            fWeights;//TH3D to hold the correction weights Axis: 0 = centrality, 1 = vertex,2 = pT. for pT<4GeV/c
  TH2D * 	    fWeightshpt;//TH2D to hold the correction weights for high pT>4GeV/c: 0 = centrality, 1 = vertex
  TF1  * 	    fpTfunction;//TF1 to hold the pT dependence over pT = 4GeV/c.
  TRandom3 *	    fRandom;//!
  TClonesArray*     fMcArray;//!
  //Objects that contain needed/used objects for the task:
  AliESDtrackCuts*  fTrackCuts;//! track cuts object
  TObject*          fCorrelator;//! correlation steering class
  Int_t*	    fRunNumberList;//! List containing the lists over runnumbers.
  Int_t 	    fNruns;// Number of runs in the choosen period
  Double_t 	    fRunFillValue;//The bin of the run we are in.
  TArrayD 	    fMBinEdges;//Contains bin edges in centrality.
  TArrayD 	    fZBinEdges;//Edges for vZ binning.
  Int_t 	    fMaxNEventMix;//Maximum number of events per bin for event mixing.
  Int_t		    fMinNofTracksMix;//Minimum number of mixing tracks.
  Int_t		    fMaxTracksperEvent;//
  TString 	    fCentralityEstimator;//! Centrality estimator ("V0M", "ZNA")
  Trigger 	    ftrigger;
  Double_t 	    fCentralityPercentile;	
  Double_t 	    fMultiplicity;
  Int_t 	    fBinVer;
  //cut variables:
  ////event:
  Double_t	    fMaxVz;//Vertex cut variable.
  Double_t 	    fMaxMult;//Is set automatically when setting the binning.
  Double_t 	    fMaxNumberOfTracksInPPConsidered;//Maximum number of tracks in an pp event used for correlations.
  Double_t 	    fNTriggers;//Contains the number of triggers in a given event.
  Double_t 	    fNAssociated;//Number of associated in the event.
  ////tracks:
  float 	    fAcceptancecut;//maximum eta
  Double_t 	    fMinTriggerPt;
  Double_t 	    fMaxTriggerPt;
  Double_t 	    fMinAssociatedPt;
  Double_t 	    fMaxAssociatedPt;
  Int_t 	    fMinNClustersTPC;
  Int_t 	    fCutMask;//To correspond to different cut masks.
  ////clusters:
  Double_t 	    fMinClusterEnergy;
  Double_t 	    fMinBCDistance;//distance to nearest bad channel
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
  ClassDef(AliAnalysisTaskCorrelation3p, 5);
};

#endif