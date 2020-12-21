#ifndef AliAnalysisTaskHFEBeautyMCTemplates_cxx
#define AliAnalysisTaskHFEBeautyMCTemplates_cxx

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




class AliAnalysisTaskHFEBeautyMCTemplates : public AliAnalysisTaskSE {
 public:
     
  AliAnalysisTaskHFEBeautyMCTemplates();
  AliAnalysisTaskHFEBeautyMCTemplates(const char *name);
  virtual ~AliAnalysisTaskHFEBeautyMCTemplates();// {}
  
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
  
  
  AliAnalysisTaskHFEBeautyMCTemplates(const AliAnalysisTaskHFEBeautyMCTemplates&); // not implemented
  AliAnalysisTaskHFEBeautyMCTemplates& operator=(const AliAnalysisTaskHFEBeautyMCTemplates&); // not implemented
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
  TH2D * fDCACharmNew3050;
  TH2D * fDCACharmNew3050IP;
  TH2D * fDCACharmNew3050OOP;
  TH2D * fDCABeautyNew;
  TH2D * fDCABeautyNewHalfRAA;
  TH2D * fDCABeautyNewHalfRAAIP;
  TH2D * fDCABeautyNewHalfRAAOOP;
  TH2D * fDCABeautyNewRAA;
  TH2D * fDCABeautyNewRAAIP;
  TH2D * fDCABeautyNewRAAOOP;
  TH2D * fDCAConversionNew;
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
  
  //AliHFEcuts * hfetrackCuts;           // Track cuts
  AliHFEsignalCuts * fSignalCuts;
  AliHFEextraCuts * fExtraCuts;
  TClonesArray *fAODArrayMCInfo;    // ! MC info particle AOD
  AliAODv0KineCuts * fAODV0Cuts;
  TRandom3 * fRd;
  
  ClassDef(AliAnalysisTaskHFEBeautyMCTemplates, 1); // example of analysis
};

#endif
