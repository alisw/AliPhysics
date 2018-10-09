#ifndef AliAnalysisTaskHFEIPDistribution_cxx
#define AliAnalysisTaskHFEIPDistribution_cxx

#include "AliAnalysisTaskSE.h"
#include "AliHFEcuts.h"

class TH1F;
class TH2F;
class TH2D;
class TH3D;
class TTree;
class THnSparse;
class AliESDEvent;
class AliESDpid;
class AliAODEvent;
class AliAODpid;
class AliAODTrack;
class AliHFEcollection;
class TArrayD;




class AliAnalysisTaskHFEIPDistribution : public AliAnalysisTaskSE {
 public:
     
  AliAnalysisTaskHFEIPDistribution();
  AliAnalysisTaskHFEIPDistribution(const char *name);
    virtual ~AliAnalysisTaskHFEIPDistribution();// {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *);
  virtual void   Process(AliAODEvent *const aodEvent=0);
  virtual void   Terminate(const Option_t*);
//  virtual void   ConnectInputData(Option_t *);
  

 private:
  //AliESDEvent *fESD;    //ESD object
  AliAODEvent *fAOD;    //ESD object
  //AliESDpid   *fESDpid;           // basic TPC object for n-sigma cuts
  TObjArray * fOutputContainer; // ! output data container
  
  Bool_t PassesTrackCuts(AliAODTrack *track);
  Bool_t PassesElectronPID(AliAODTrack *track, AliPIDResponse *pid);
  
  AliAnalysisTaskHFEIPDistribution(const AliAnalysisTaskHFEIPDistribution&); // not implemented
  AliAnalysisTaskHFEIPDistribution& operator=(const AliAnalysisTaskHFEIPDistribution&); // not implemented
  TH2D * EPCent;
  TH2D * EPCentCorrected;
  TH2D * EPCentV0A;
  TH2D * EPCentV0C;
  
  
  TH2D * fIPData;
  
  TH1D * EP2040;
  TH1D * EP2040Corrected;
  TH1D * EP2040V0A;
  TH1D * EP2040V0C;
  
  TH2D * TPCnSigma;
  
  TH1D * DeltaPhi;

  TH2D * fpTIP2040IP;
  TH2D * fpTIP2040OOP;
  TH2D * fpTIP3050IP;
  TH2D * fpTIP3050OOP;
  
  TH1D * EventSelectionSteps;
  
  //AliHFEcuts * hfetrackCuts;           // Track cuts
  AliHFEextraCuts * fExtraCuts;
  
  
  ClassDef(AliAnalysisTaskHFEIPDistribution, 1); // example of analysis
};

#endif
