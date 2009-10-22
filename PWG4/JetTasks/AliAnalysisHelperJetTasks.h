#ifndef  ALIANALYSISHELPERJETTASKS_H
#define  ALIANALYSISHELPERJETTASKS_H
#include "TObject.h"

class AliMCEvent;
class AliAODJet;
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

  static void MergeOutput(char* cFiles, char* cList = "pwg4spec"); // Merges the files in the input text file  needs the two histograms fh1PtHard_Trials, fh1Xsec and the name of the input list
  static Bool_t PythiaInfoFromFile(const char* currFile,Float_t &fXsec,Float_t &fTrials);// get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  static Bool_t GetThrustAxis(TVector3 &n01, TVector3 * pTrack, const Int_t &nTracks);
  static Bool_t GetEventShapes(TVector3 &n01, TVector3 * pTrack, Int_t nTracks, Double_t * eventShapes);
  enum {kMaxJets = 6}; //  needed for array size not to fragemnt memory on the heap by many new/delete 
  private:
  
  ClassDef(AliAnalysisHelperJetTasks, 1) 
};

#endif // ALIANALYSISHELPERJETTASKS_H
