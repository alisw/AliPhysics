#ifndef AliAnalysisTaskHFEIPCorrection_cxx
#define AliAnalysisTaskHFEIPCorrection_cxx

#include "AliAnalysisTaskSE.h"
#include "AliHFEcuts.h"

class TH1F;
class TH2F;
class TH2D;
class TH3D;
class TString;
class TTree;
class THnSparse;
class AliESDEvent;
class AliESDpid;
class AliAODEvent;
class AliAODpid;
class AliAODTrack;
class AliHFEcollection;
class TArrayD;
class AliAODVertex;
class TRandom3;
class TSpline3;




class AliAnalysisTaskHFEIPCorrection : public AliAnalysisTaskSE {
 public:
     
  AliAnalysisTaskHFEIPCorrection();
  AliAnalysisTaskHFEIPCorrection(const char *name);
    virtual ~AliAnalysisTaskHFEIPCorrection();// {}
  
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
  Bool_t PassesTrackCutsNoFirst(AliAODTrack *track);
  Bool_t PassesElectronPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesWeakerElectronPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesPionPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesKaonPID(AliAODTrack *track, AliPIDResponse *pid);
  Bool_t PassesMinimalTrackCuts(AliAODTrack *track);
  Bool_t PassesITSTrackCuts(AliAODTrack *track);
  Int_t ReturnRunBin(Int_t RunNr);

  Int_t IsInMisalignedRegion(AliAODTrack *track, double vtxz); // 15o pass1 misalignement
  AliAODVertex * CorrectVertex(AliAODEvent *aodEvent, double vtxz); // Vertex without using excluded regions
  void GetTrackImpactParameter(AliAODEvent *aodEvent, AliAODTrack *track, AliAODVertex * pvtx, Double_t &dcaxy); // Calculate IP from other vertex

  void GetCorrectedImpactParameter(AliAODEvent *aodEvent, AliAODTrack *track, Double_t primVertexZ, Double_t &dcaxy); // correct for effects in phi, z, and pt
  TString GetPeriodNameByLPM(TString lTag);
  
  AliAnalysisTaskHFEIPCorrection(const AliAnalysisTaskHFEIPCorrection&); // not implemented
  AliAnalysisTaskHFEIPCorrection& operator=(const AliAnalysisTaskHFEIPCorrection&); // not implemented
  
  
  TH2D * fIPData;
  
  TH1D * fCentrality;
  TH1D * EP2040;
  TH1D * EP2040Corrected;
  TH1D * EP2040V0A;
  TH1D * EP2040V0C;
  
  TH2D * TPCnSigma;
  TH3D * fTPCnSigmaCentIP;
  TH3D * fTPCnSigmaCentOOP;
  TH3D * fEleV0TPCnSigmaCentIP;
  TH3D * fEleV0TPCnSigmaCentOOP;

  TH2D * EPCent;
  TH2D * EPCentUncorrected;
  TH2D * EPCentV0A;
  TH2D * EPCentV0C;

  TH1D * DeltaPhi;
  AliHFEextraCuts * fExtraCuts;

  TH2D * fpTIP2040IP;
  TH2D * fpTIP2040OOP;
  TH2D * fpTIP3050IP;
  TH2D * fpTIP3050OOP;
  TH2D * fpTIP3050IPAlternativeCut; // TPC cut of -1 instead of -0.5 to combat
  TH2D * fpTIP3050OOPAlternativeCut; // loss of efficiency due to shifted center
  
  TH1D * EventSelectionSteps;

  TH3D * fPionV0pTRNoCuts;
  TH3D * fPionV0pTRWithCuts;
  TH3D * fPionV0pTRNoCutsIP;
  TH3D * fPionV0pTRWithCutsIP;
  TH3D * fPionV0pTRNoCutsOOP;
  TH3D * fPionV0pTRWithCutsOOP;
  TH2D * fPionV0pTTPC;
  TH2D * fPionV0pTTPCWithCuts;
  TH2D * fPionV0pTTPCIP;
  TH2D * fPionV0pTTPCOOP;
  TH2D * fPionV0pTTPCIPWTOF; // These are actually electron V0s
  TH2D * fPionV0pTTPCOOPWTOF;
  TH2D * fPionV0pTTPCIPnoFirst;
  TH2D * fPionV0pTTPCOOPnoFirst;
  TH2D * fPionV0pTTPCIPWTOFnoFirst;
  TH2D * fPionV0pTTPCOOPWTOFnoFirst;

  TH3D * fEPLowHighCent;
  TH3D * fEPLowVZEROCent;
  TH3D * fEPHighVZEROCent;
  TH3D * fEPLowHighCent2;
  TH3D * fEPLowVZEROCent2;
  TH3D * fEPHighVZEROCent2;
  
  TH3D * fDCARegionRun;
  TH3D * fDCAPhiZHadrons;
  TH3D * fDCAPhiZHadronsEarlyRuns;
  TH3D * fDCAPhiZHadronsLateRuns;
  TH3D * fDCAPhiZHadronsC;
  TH3D * fDCAPhipTHadrons;
  TH3D * fDCAPhipTHadronsEarlyRuns;
  TH3D * fDCAPhipTHadronsLateRuns;
  TH3D * fDCAPhipTHadronsC;
  TH3D * fDCAPhiZKaons;
  TH3D * fDCAPhiZKaonsC;
  TH3D * fDCAPhipTKaons;
  TH3D * fDCAPhipTKaonsC;
  TH3D * fpTPhiZHadrons;
  TH3D * fDCAWErrHadrons;
  TH2D * fDCAHadrons;
  TH3D * fDCAHadronsFineBins;
  TH2D * fDCAKaons; // Should have less contamination, but have higher mass
  TH3D * fDCAWErrKaons;
  TH3D * fDCAKaonsFineBins;


  
  //AliHFEcuts * hfetrackCuts;           // Track cuts
  AliAODv0KineCuts * fAODV0Cuts;
  TRandom3 * fRd;
  TSpline3 * fSplineCorr;
  
  ClassDef(AliAnalysisTaskHFEIPCorrection, 5);
};

#endif
