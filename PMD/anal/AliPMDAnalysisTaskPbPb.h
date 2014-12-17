#ifndef AliPMDAnalysisTaskPbPb_cxx
#define AliPMDAnalysisTaskPbPb_cxx

 /**************************************************************************

                   A template class to read tracks (PMD Cluster)
		         Runs in Local and Grid Modes
			Can be used for PbPb PMD analysis
		      Origin: Satyajit Jena <sjena@cern.ch>

 **************************************************************************/


class TH1F;
class AliESDEvent;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliPMDAnalysisTaskPbPb : public AliAnalysisTaskSE {
 public:
 AliPMDAnalysisTaskPbPb(): AliAnalysisTaskSE(), fOutputList(0), fTrackCuts(0),fESD(0), fHistPt(0), fHistEta(0) {}
  AliPMDAnalysisTaskPbPb(const char *name);
  virtual ~AliPMDAnalysisTaskPbPb() {}
  
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetIsMC(Bool_t isMC) { fIsMC = isMC; }  
 
 private:
  TList *fOutputList;
  
  AliESDtrackCuts *fTrackCuts;
 
  AliESDEvent *fESD;     // ESD object
 
  TH1F    *fHistPt;  // Pt spectrum
  TH1F   *fHistEta; // Pt spectrum
  TH2F   *fhEsdXYP; //
  TH2F   *fhEsdXYC; //

  Bool_t fIsMC; // MC truth 

  AliPMDAnalysisTaskPbPb(const AliPMDAnalysisTaskPbPb&); // not implemented
  AliPMDAnalysisTaskPbPb& operator=(const AliPMDAnalysisTaskPbPb&); // not implemented
  
  ClassDef(AliPMDAnalysisTaskPbPb, 1); // example of analysis
};

#endif
