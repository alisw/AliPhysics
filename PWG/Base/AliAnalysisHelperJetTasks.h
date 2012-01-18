#ifndef  ALIANALYSISHELPERJETTASKS_H
#define  ALIANALYSISHELPERJETTASKS_H
#include "TObject.h"

class AliMCEvent;
class AliAODJet;
class AliVEvent;
class TString;
class TArrayI;
class TArrayF;
class AliGenPythiaEventHeader;
class TVector3;
class AliGenEventHeader;


//
// Helper Class that contains a lot of 
// usefull static functions jet matchin pythia access etc.
//


class AliAnalysisHelperJetTasks : public TObject {
 public:
  AliAnalysisHelperJetTasks() : TObject() {;}
  virtual ~AliAnalysisHelperJetTasks(){;}


  enum {kMaxJets = 6}; //  needed for array size not to fragemnt memory on the heap by many new/delete 


  enum { kNone = 1<<0,
	 kBunchBunch = 1<<1,
	 kBunchEmpty = 1<<2,
	 kEmptyEmpty= 1<<3,
	 kV0A=1<<4,
	 kV0C=1<<5,
	 kNoV0BG=1<<6,
	 kSPDFO=1<<7,
	 kPhysicsSelection = 1<<8, 
	 kVertexIn = 1<<9, 
	 kIsCosmic = 1<<10, 
	 kIsPileUp = 1<<11,
	 kIsMCND=1<<12,
	 kIsMCDD=1<<13,
	 kIsMCSD=1<<14,
	 kTotalSelections = (1<<15) - 1};

  enum Trigger {kAcceptAll = 0,kMB1,kMB2,kMB3,kSPDGFO,kTrigger}; // 
  // same as in PWG0Helper
  enum MCProcessType { kInvalidProcess = -1, kND = 0x1, kDD = 0x2, kSD = 0x4, kOnePart = 0x8 };

  static AliGenPythiaEventHeader*  GetPythiaEventHeader(const AliMCEvent *mcEvent);
  static void PrintStack(AliMCEvent *mcEvent,Int_t iFirst = 0,Int_t iLast = 0,Int_t iMaxPrint = 10);
  static void GetClosestJets(AliAODJet *genJets,
			     const Int_t &kGenJets,
			     const AliAODJet *recJets,
			     const Int_t &kRecJets,
			     Int_t *iGenIndex,
			     Int_t *iRecIndex,
			     Int_t iDebug = 0, Float_t maxDist = 0.3);

  static void GetClosestJets(const TList *genJetsList,const Int_t &kGenJets,
			     const TList *recJetsList,const Int_t &kRecJets,
			     TArrayI &iGenIndex,TArrayI &iRecIndex,
			     Int_t iDebug = 0,Float_t maxDist = 0.3);
				 
  static void GetJetMatching(const TList *genJetsList, const Int_t &kGenJets,
                             const TList *recJetsList, const Int_t &kRecJets,
                             TArrayI &iMatchIndex, TArrayF &fPtFraction,
                             Int_t iDebug = 0, Float_t maxDist = 0.3, Int_t mode=1);
							 
  static Double_t GetFractionOfJet(const AliAODJet *recJet,const AliAODJet *genJet,Int_t mode=1);


  static void MergeOutputDirs(const char* cFiles,const char* cPattern,const char *cOutFile,Bool_t bUpdate = false); // merges all directories containing the pattern

  static void MergeOutput(const char* cFiles,const char* cDir = "",const char *cList = "",const char* cOutFile ="allpt.root",Bool_t bUpdate = false); // Merges the files in the input text file  needs the two histograms fh1PtHard_Trials, fh1Xsec and the name of the input list
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  static Bool_t PrintDirectorySize(const char* currFile,Int_t iDetail = -1); // print the size of the output on a given file
  static Bool_t GetEventShapes(TVector3 &n01,const TVector3 * pTrack, Int_t nTracks, Double_t * eventShapes);

  static MCProcessType GetPythiaEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
  static MCProcessType GetDPMjetEventProcessType(AliGenEventHeader* aHeader, Bool_t adebug = kFALSE);
  static Int_t GetLastProcessType() { return fgLastProcessType; }

  static Bool_t Selected(Bool_t bSet = kFALSE,Bool_t bNew = kTRUE); // static function to store the state of selection from service task
  static Double_t ReactionPlane(Bool_t bSet = kFALSE,Double_t fNew = 0); 
  static Int_t GetPhiBin(Double_t phi,Int_t fNRPbins);

  static Bool_t IsPileUp(); // Wrapper for SelectInfo with PileUp
  static Bool_t IsCosmic(); // Wrapper for SelectInfo with cosmic
  static Bool_t TestSelectInfo(UInt_t iMask); // Wrapper for testing the SelectInfo bitmask
  static Bool_t TestEventClass(Int_t iClass); // Wrapper for testing the SelectInfo bitmask

  static UInt_t SelectInfo(Bool_t bSet = kFALSE,UInt_t iNew = 0); // static function to store the state bitmask of the selection from service task
  static Int_t  EventClass(Bool_t bSet = kFALSE,Int_t iNew = 0); // static function to store the event class of the selection from service task
  
  // these methods have been essentially copied from PWG0/AliTriggerAnalysis and expanded to use with AOD
  static Bool_t IsTriggerFired(const AliVEvent* aEsd, Trigger trigger);

  private:
  
  static Int_t fgLastProcessType;    // stores the raw value of the last process type extracted
 
  ClassDef(AliAnalysisHelperJetTasks, 7) 
};

#endif // ALIANALYSISHELPERJETTASKS_H
