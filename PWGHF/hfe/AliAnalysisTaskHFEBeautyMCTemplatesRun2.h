#ifndef AliAnalysisTaskHFEBeautyMCTemplatesRun2_cxx
#define AliAnalysisTaskHFEBeautyMCTemplatesRun2_cxx

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
class AliHFEsignalCuts;
class AliAODv0KineCuts;
class TRandom3;




class AliAnalysisTaskHFEBeautyMCTemplatesRun2 : public AliAnalysisTaskSE {
 public:
     
  AliAnalysisTaskHFEBeautyMCTemplatesRun2();
  AliAnalysisTaskHFEBeautyMCTemplatesRun2(const char *name);
  virtual ~AliAnalysisTaskHFEBeautyMCTemplatesRun2();// {}
  
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
  
  
  AliAnalysisTaskHFEBeautyMCTemplatesRun2(const AliAnalysisTaskHFEBeautyMCTemplatesRun2&); // not implemented
  AliAnalysisTaskHFEBeautyMCTemplatesRun2& operator=(const AliAnalysisTaskHFEBeautyMCTemplatesRun2&); // not implemented
  Int_t FindSource(AliMCParticle * mcple, AliMCEvent* fMCEvent, Double_t &motherPt, Double_t &GroundStateMotherPt);
  Bool_t PassesTrackCuts(AliAODTrack *track);
  Bool_t PassesMinimalTrackCuts(AliAODTrack *track);
  Bool_t PassesElectronPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesPionPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesKaonPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesITSTrackCuts(AliAODTrack *track);
  Bool_t IsAddedSignal(AliMCParticle * mcple);
  void PrintHierarchy(AliMCParticle * mcple, AliMCEvent* fMCEvent);
  Int_t MotherPDG(AliMCParticle * mcple, AliMCEvent* fMCEvent);
  Int_t CharmSource(AliMCParticle * mcple, AliMCEvent* fMCEvent);
  Int_t BeautySource(AliMCParticle * mcple, AliMCEvent* fMCEvent);

  TString GetPeriodNameByLPM(TString lTag);
  
  TH1D * fCentrality;
  //TH2D * fSourceGenerator;
  // A lot of different histograms for cross checks and systematics so I only need to run the train once
  TH2D * fDCACharm;
  TH2D * fDCACharm3050;
  TH2D * fDCACharm3050IP;
  TH2D * fDCACharm3050OOP;
  TH2D * fDCABeauty;
  TH2D * fDCABeautyHalfRAA;
  TH2D * fDCABeautyHalfRAAIP;
  TH2D * fDCABeautyHalfRAAOOP;
  TH2D * fDCABeautyRAA;
  TH2D * fDCABeautyRAAIP;
  TH2D * fDCABeautyRAAOOP;
  TH2D * fDCAConversion;
  TH2D * fDCADalitz;
  TH2D * fDCACharmNew;
  TH3D * fDCACharmNewDPlus;
  TH3D * fDCACharmNewDZero;
  TH3D * fDCACharmNewDS;
  TH3D * fDCACharmNewLambdaC;
  TH3D * fDCACharmNewOtherC;
  TH2D * fDCACharmNew3050;
  TH2D * fDCACharmNew3050IP;
  TH2D * fDCACharmNew3050OOP;
  TH3D * fDCACharmWeightedNew3050;
  TH3D * fDCACharmWeightedNew3050IP;
  TH3D * fDCACharmWeightedNew3050OOP;
  TH2D * fDCABeautyNew;
  TH3D * fDCABeautyNewBZero;
  TH3D * fDCABeautyNewBPlus;
  TH3D * fDCABeautyNewBS;
  TH3D * fDCABeautyNewLambdaB;
  TH3D * fDCABeautyNewOtherB;
  TH2D * fDCABeautyNewHalfRAA;
  TH2D * fDCABeautyNewHalfRAAIP;
  TH2D * fDCABeautyNewHalfRAAOOP;
  TH3D * fDCABeautyWeightedNewHalfRAA;
  TH3D * fDCABeautyWeightedNewHalfRAAIP;
  TH3D * fDCABeautyWeightedNewHalfRAAOOP;
  TH2D * fDCABeautyNewRAA;
  TH2D * fDCABeautyNewRAAIP;
  TH2D * fDCABeautyNewRAAOOP;
  TH2D * fDCAConversionNew;
  TH3D * fDCAConversionNewCent;
  TH2D * fDCADalitzNew;
  TH2D * fDCADalitzCharm; // c -> x -> e
  TH2D * fDCADalitzBeauty;
  TH2D * fDCAConversionCharm; // c -> gamma  -> e
  TH2D * fDCAConversionBeauty;
  TH1D * fBeautyMotherpT;
  TH1D * fCharmMotherpT;
  TH1D * fGroundStateBeautyMotherpT;
  TH1D * fGroundStateCharmMotherpT;
  TH2D * fDCAHadrons;
  TH3D * fDCAWErrHadrons; // Pions, mostly
  TH3D * fDCAHadronsFineBins;
  TH2D * fDCAKaons; // Should have less contamination, but have higher mass
  TH3D * fDCAWErrKaons;
  TH3D * fDCAKaonsFineBins;
  
  TH2D * fDCAHadronsCorrected;
  TH3D * fPionV0pTRNoCuts;  // pt, R, cent
  TH3D * fPionV0pTRWithCuts;
  TH3D * fPionTPCSignal; // To check TPC signal shape in MC
  
  //AliHFEcuts * hfetrackCuts;           // Track cuts
  AliHFEsignalCuts * fSignalCuts;
  AliHFEextraCuts * fExtraCuts;
  TClonesArray *fAODArrayMCInfo;    // ! MC info particle AOD
  AliAODv0KineCuts * fAODV0Cuts;
  TRandom3 * fRd;
  
  ClassDef(AliAnalysisTaskHFEBeautyMCTemplatesRun2, 2); // example of analysis
};

#endif
