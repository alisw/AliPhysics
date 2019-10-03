#ifndef ALIMUONQAMERGESUBMITTER_H
#define ALIMUONQAMERGESUBMITTER_H

#include "AliMuonGridSubmitter.h"
#include "TString.h"

#include <vector>

/**
  \ingroup pwg_muondep_submitter
  \class AliMuonQAMergeSubmitter
  \brief QA merging Grid jobs submitter
  \author Laurent Aphecetche, Subatech

  @todo This class has not be completed, thus there is no guarantee it's working in any way...
  */

class AliMuonQAMergeSubmitter : public AliMuonGridSubmitter
{
public:
  AliMuonQAMergeSubmitter(const char* period, const char* pass);
  virtual ~AliMuonQAMergeSubmitter();
  
  Bool_t Run(const char* mode);

  Bool_t Submit(Int_t runNumber, Bool_t dryRun);
  
  Int_t Submit(Bool_t dryRun);

  TString MergeJDLName(Bool_t final) const { return (final ? "QAMerge_final.jdl" : "QAMerge.jdl"); }

  Bool_t Generate(const char* name) const;

  Bool_t SetRemoteDir(const char* dir);

  UInt_t MakeXMLCollectionForRun(Int_t runNumber, Int_t stage);

  UInt_t GetSplitMaxInputFileNumber() const { return fSplitMaxInputFileNumber; }
  
  void SetSplitMaxInputFileNumber(UInt_t n) { fSplitMaxInputFileNumber=n; }
  
  virtual void Print(Option_t* opt="") const;
  
  void ShowStages();
  
  void ShowStage(Int_t runNumber);

private:

  TString fPeriod;
  TString fPass;
  TString fWhatToMerge; // file to be merged
  UInt_t fSplitMaxInputFileNumber;
  
/// \cond CLASSIMP
  ClassDef(AliMuonQAMergeSubmitter,1);
/// \endcond
};

#endif
