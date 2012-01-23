#ifndef AliAnalysisTaskAODvsESD_cxx
#define AliAnalysisTaskAODvsESD_cxx

/* $Id$ */ 

class TList;
class TH2F;
class AliESDEvent;
class AliAODEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskAODvsESD : public AliAnalysisTask {
 public:
  AliAnalysisTaskAODvsESD(const char *name = "AliAnalysisTaskAODvsESD");
  virtual ~AliAnalysisTaskAODvsESD() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(const Option_t*);
  
 private:
  AliESDEvent *fESD;         // ESD object
  AliAODEvent *fAOD;         // AOD object

  TList   *fList;             // list of ntuples
  TNtuple *fMuonNtuple;       // NTuple for single muons ESD
  TNtuple *fMuonNtupleAOD;    // NTuple for single muons AOD
  TH1F    *fInvMass;          // ESD inv. mass histo
  TH1F    *fInvMassAOD;       // AOD inv. mass histo
   
  ClassDef(AliAnalysisTaskAODvsESD, 1); // example of analysis
};

#endif

