#ifndef  ALIANALYSISHELPERJETTASKS_H
#define  ALIANALYSISHELPERJETTASKS_H
#include "TObject.h"

class AliMCEvent;
class AliAODJet;
class AliVEvent;
class TString;
class AliGenPythiaEventHeader;
class TVector3;

// Helper Class that contains a lot of usefull static functions (i.e. for Flavor selection.

class AliAnalysisHelperJetTasks : public TObject {
 public:
  AliAnalysisHelperJetTasks() : TObject() {;}
  virtual ~AliAnalysisHelperJetTasks(){;}
  
  static AliGenPythiaEventHeader*  GetPythiaEventHeader(AliMCEvent *mcEvent);
  static void PrintStack(AliMCEvent *mcEvent,Int_t iFirst = 0,Int_t iLast = 0,Int_t iMaxPrint = 10);
  static void GetClosestJets(AliAODJet *genJets,
			     const Int_t &kGenJets,
			     AliAODJet *recJets,
			     const Int_t &kRecJets,
			     Int_t *iGenIndex,
			     Int_t *iRecIndex,
			     Int_t iDebug, Float_t maxDist = 0.5);

  static void MergeOutput(char* cFiles, char* cDir = "",char *cList = "",char* cOutFile ="allpt.root",Bool_t bUpdate = false); // Merges the files in the input text file  needs the two histograms fh1PtHard_Trials, fh1Xsec and the name of the input list
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  static Bool_t PrintDirectorySize(const char* currFile); // print the size of the output on a given file
  static Bool_t GetEventShapes(TVector3 &n01, TVector3 * pTrack, Int_t nTracks, Double_t * eventShapes);


  enum {kMaxJets = 6}; //  needed for array size not to fragemnt memory on the heap by many new/delete 
  enum { kNone = 1<<0,kBunchBunch = 1<<1,kBunchEmpty = 1<<2,kEmptyEmpty= 1<<3,
	 kPhysicsSelection = 1<<4, kVertexIn = 1<<5, kIsCosmic = 1<<6, kIsPileUp = 1<<7,kTotalSelections = (1<<8) - 1};

  enum Trigger {kAcceptAll = 0,kMB1,kMB2,kMB3,kSPDGFO,kTrigger}; // 

  static Bool_t Selected(Bool_t bSet = kFALSE,Bool_t bNew = kTRUE); // static function to store the state of selection from service task

  static Bool_t IsPileUp(); // Wrapper for SelectInfo with PileUp
  static Bool_t IsCosmic(); // Wrapper for SelectInfo with cosmic

  static UInt_t SelectInfo(Bool_t bSet = kFALSE,UInt_t iNew = 0); // static function to store the state bitmask of the selection from service task
  
  // these methods have been essentially copied from PWG0/AliTriggerAnalysis and expanded to use with AOD
  static Bool_t IsTriggerFired(const AliVEvent* aEsd, Trigger trigger);

  private:
  
  ClassDef(AliAnalysisHelperJetTasks, 2) 
};

#endif // ALIANALYSISHELPERJETTASKS_H
