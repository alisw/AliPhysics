#include "TH2.h"
#include "TObjArray.h"

#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"

class AliAnalysisTaskSingleMu : public AliAnalysisTask {
 public:
  AliAnalysisTaskSingleMu(const char *name = "AliAnalysisTaskSingleMu");
  virtual ~AliAnalysisTaskSingleMu() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

protected:
  Bool_t MuonPassesCuts(AliAODTrack &muonTrack,
			TLorentzVector &lorVec,
			Int_t &trigMatch);

  const AliAODVertex* GetVertex();
  void ResetHistos();

private:
  AliAODEvent *fAOD; //!< ESDevent object

  static const Int_t fgkNhistos = 1;
  static const Int_t fgkNTrigCuts = 4;

  enum
  {
    kNoMatchTrig,
    kAllPtTrig,
    kLowPtTrig,
    kHighPtTrig
  };

  TString trigName[fgkNTrigCuts]; //!< trigger cut names 

  TObjArray * fOutputContainer; //!< output data container

  TH2F *fVzVsPt[fgkNTrigCuts]; //!< Single muon spectrum
     
  ClassDef(AliAnalysisTaskSingleMu, 0); // Single muon analysis
};

