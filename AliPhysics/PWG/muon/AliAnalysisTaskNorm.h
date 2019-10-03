#ifndef ALIANALYSISTASKNORM_H
#define ALIANALYSISTASKNORM_H
 
/* $Id$ */ 
//
// AliAnalysisTaskNorm
// Analysis task for muon trigger normalization 
// 
// author: C. Hadjidakis 
//

#ifndef ALIANALYSISTASKSE_H
#  include "AliAnalysisTaskSE.h"
#endif

class AliCounterCollection;
class AliMuonEventCuts;
class AliVEvent;

class AliAnalysisTaskNorm : public AliAnalysisTaskSE
{
 public:

  AliAnalysisTaskNorm(const char *name = "<default name>");
  virtual ~AliAnalysisTaskNorm();

  virtual void Print(Option_t *opt="") const;
  void Print(TObjArray *obj) const;
  virtual void Terminate(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* opt);

  //
  static Bool_t IsAODEvent( const AliVEvent *event);
  static TString GetFiredTriggerClasses( const AliVEvent *event);

  //build list of counters
  TList *BuildListOfNormFactor(const TObjArray*);
  TList *BuildListOfTrigger(const TObjArray*);
  TList *BuildListOfCentrality(AliCentrality *centrality);  
  TList *BuildListOfPileUp( const AliVEvent *event);
  TList *BuildListOfTracklets( const AliVEvent *event);
  TList *BuildListOfV0AMult( const AliVEvent *event);

  //fill histograms
  void FillHistoPileUpVertices( const AliVEvent *event, const TObjArray*);
  void FillHistoMult( const AliVEvent *event, const TObjArray*);

  Bool_t CheckPattern(TString, TObjArray*, TObjArray*);
  TObjArray *BuildArrayOfTrigger(const TObjArray*, TString keepPattern="", TString rejectPattern="OTHER,TRUE,PHI,ANY,EMC,-ACE-,-ABCE-,WU,MUP,SPI,SHM");
  //
  Bool_t IsPileupFromSPDInMultBins(const AliVEvent *event) const; 

  TString GetOCDBPath() {return "raw://";}; // OCDB to be used (raw:// by default)

  //Fill object methods
  void FillEventCounters(Int_t, TList*,TList*,TList*,TList*,TList*,TList*,Bool_t,Bool_t,TString,Bool_t);

  //
  AliMuonEventCuts* GetMuonEventCuts() { return fMuonEventCuts;}
  void SetIsMC(Bool_t isMC) { fIsMC = isMC;};
  void SetIsESD(Bool_t isESD) { fIsESD = isESD;};
  void SetBeamConf(TString sBeamConf) { fBeamConf = sBeamConf;};

 private:

  enum eListV0A {
    kV0AMB  = 0, ///<V0A mult for CINT7 trigger
    kV0AMUL = 1, ///<V0A mult for CMUL7 trigger
    kV0AMultvsCentMB = 2, ///<V0A mult vs centrality for CMUL7 trigger
    kV0ACentvsV0CCentMUL = 3, ///<V0A vs V0C centrality for CMUL7 trigger
    kV0ACentvsCL1CentMUL = 4, ///<V0A vs CL1 centrality for CMUL7 trigger
    kV0ACentvsV0CCentMB = 5, ///<V0A vs V0C centrality for CINT7 trigger
    kV0ACentvsCL1CentMB = 6, ///<V0A vs CL1 centrality for CINT7 trigger
    kV0CCentvsCL1CentMB = 7, ///<V0A vs CL1 centrality for CINT7 trigger
    kV0ACentMB  = 8, ///<V0A centrality for CINT7 trigger
    kV0ACentMUL = 9, ///<V0A centrality for CMUL7 trigger
    kV0AvsSPDTrackletsMB = 10, ///<V0A mult vs SPD tracklets nr for CINT7 trigger
    kV0AvsSPDTrackletsMUL = 11, ///<V0A mult vs SPD tracklets nr for CMUL7 trigger
    kV0AvsSPDTrackletsMBPuCut = 12, ///<V0A mult vs SPD tracklets nr for CINT7 trigger
    kV0AvsSPDTrackletsMULPuCut = 13, ///<V0A mult vs SPD tracklets nr for CMUL7 trigger
    kV0AvsSPDTrackletsMBPuCut2 = 14, ///<V0A mult vs SPD tracklets nr for CINT7 trigger
    kV0AvsSPDTrackletsMULPuCut2 = 15, ///<V0A mult vs SPD tracklets nr for CMUL7 trigger
    kV0ACentMBPuCut  = 16, ///<V0A centrality for CINT7 trigger
    kV0ACentMULPuCut = 17, ///<V0A centrality for CMUL7 trigger
    kV0ACentMBPuCut2  = 18, ///<V0A centrality for CINT7 trigger
    kV0ACentMULPuCut2 = 19 ///<V0A centrality for CMUL7 trigger
  };
  
  enum eListZN {
    kZNMB  = 0, ///<ZN energy for CINT7 trigger
    kZNMUL = 1, ///<ZN energy for CMUL7 trigger
    kZNCentMB  = 2, ///<ZN cent for CINT7 trigger
    kZNCentMUL = 3, ///<ZN cent for CMUL7 trigger
    kZNMultvsCentMB = 4, ///<ZN mult vs centrality for CINT7 trigger
    kZNvsSPDTrackletsMB = 5, ///<ZN energy vs SPD tracklets nr for CINT7 trigger
    kZNvsSPDTrackletsMUL = 6, ///<ZN energy vs SPD tracklets nr for CMUL7 trigger
    kZNvsSPDTrackletsMBPuCut = 7, ///<ZN energy vs SPD tracklets nr for CINT7 trigger w SPD pu
    kZNvsSPDTrackletsMULPuCut = 8, ///<ZN energy vs SPD tracklets nr for CMUL7 trigger w SPD pu
    kZNvsSPDTrackletsMBPuCut2 = 9, ///<ZN energy vs SPD tracklets nr for CINT7 trigger w MV pu
    kZNvsSPDTrackletsMULPuCut2 = 10, ///<ZN energy vs SPD tracklets nr for CMUL7 trigger w MV pu
    kZNCentMBPuCut  = 11, ///<ZN cent for CINT7 trigger with SPD pu
    kZNCentMULPuCut = 12, ///<ZN cent for CMUL7 trigger with SPD pu
    kZNCentMBPuCut2  = 13, ///<ZN cent for CINT7 trigger with MV pu
    kZNCentMULPuCut2 = 14 ///<ZN cent for CMUL7 trigger with MV pu
   };

  enum eListVertex {
    kVZMB = 0, ///< primary vertex Z
    kVnCMB = 1, ///<primary vertex nContributors
    kPileupVZMB = 2, ///< pileup vertices Z
    kPileupnCMB = 3, ///< pileup nContributors
    kVZMUL = 4, ///< primary vertex Z
    kVnCMUL = 5, ///<primary vertex nContributors
    kPileupVZMUL = 6, ///< pileup vertices Z
    kPileupnCMUL = 7, ///< pileup nContributors
    kDeltaZMB = 8,///<distrance between primary and pileupvertices
    kDeltaZMUL = 9,///<distrance between primary and pileupvertices
    kNPileupMB = 10,///<number of pileupvertices
    kNPileupMUL = 11///<number of pileupvertices
  };

  AliCounterCollection *fRunCounters;//! run counters
  AliCounterCollection *fEventCounters;//! event counters
  TObjArray *fListVertex; //!list of objects vertex related
  TObjArray *fListV0A; //!list of objects V0A related
  TObjArray *fListZN; //!list of objects ZN related
  Bool_t fIsESD; //ESD or AOD
  Bool_t fIsMC; //MC or data
  TString fBeamConf; //Beam configuration (p-Pb,Pb-p)
  AliMuonEventCuts *fMuonEventCuts; // cuts for muon events 

  TObjArray *fSCentEst; // array of centrality estimators
  TArrayI fCentBin; //centrality range for event counters
  TObjArray *fSCentBin; // array of centrality bins
  TArrayI fTrackletsBin; //ntracklets range for event counters
  TObjArray *fSTrackletsBin; //array of ntracklets bins
  TArrayI fV0AMultBin; //V0Amult range for event counters
  TObjArray *fSV0AMultBin; //array of V0Amult bins
  
  AliAnalysisTaskNorm(const AliAnalysisTaskNorm&);// not implemented
  AliAnalysisTaskNorm& operator=(const AliAnalysisTaskNorm&); //not implemented;

  ClassDef(AliAnalysisTaskNorm,4); 

};

#endif
