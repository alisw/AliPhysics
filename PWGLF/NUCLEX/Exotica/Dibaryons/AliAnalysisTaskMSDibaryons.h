#ifndef ALIANALYSISTASKMSDIBARYONS_H
#define ALIANALYSISTASKMSDIBARYONS_H

// ROOT includes
#include <THashList.h>
#include <TH1.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisUtils.h>
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliEventPoolManager.h>

#include <AliMultSelection.h>
#include <AliEventplane.h>
#include <AliQnCorrectionsManager.h>
#include <AliAnalysisTaskFlowVectorCorrections.h>
#include <AliQnCorrectionsQnVector.h>

class AliAnalysisTaskMSDibaryons : public AliAnalysisTaskSE {
 public:
  enum FlowMethod{
    kOFF = -1,
    kEP  = 0,
      kSP  = 1
  };

  enum QnDetector{
    kNone      = -1,
    kFullTPC   = 0,
    kTPCNegEta = 1,
    kTPCPosEta = 2,
    kFullV0    = 3,
    kV0A       = 4,
    kV0C       = 5
  };
  AliAnalysisTaskMSDibaryons();
  AliAnalysisTaskMSDibaryons(const char *name);
  virtual ~AliAnalysisTaskMSDibaryons();
  
  virtual void  SetTrigger(UInt_t ktriggerInt=AliVEvent::kINT7){ ftrigBit=ktriggerInt; }
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);

  //setter for event selection
  void SetCentralityEstimator(TString estimator) {fEstimator = estimator;}
  void SetCentralityMin(Float_t min) {fCentralityMin = min;}
  void SetCentralityMax(Float_t max) {fCentralityMax = max;}

 protected:
  Bool_t   EventSelection(AliAODEvent *data);
  Bool_t   ProtonSelection(AliAODTrack *trk);
  Double_t InvMassLambda(AliAODcascade *casc);
  Double_t InvMassAntiLambda(AliAODcascade *casc);
  Double_t InvMassXi(AliAODcascade *casc);
  Double_t InvMassOmega(AliAODcascade *casc);
  Double_t LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const;
  Double_t xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const;
  Double_t InvMasslambda(AliAODv0 *v0);
  Double_t InvMassK0(AliAODv0 *v0);
  Double_t InvMassAntilambda(AliAODv0 *v0);
  Double_t CalculateInvMassAntilambda(AliAODv0 *v0,AliAODTrack *antiprotontrk,AliAODTrack *piontrk);
  Double_t CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
		     Double_t px2,Double_t py2,Double_t pz2);
  Double_t InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
			 Double_t px2,Double_t py2,Double_t pz2,Double_t energysum);
  //Bool_t ExtractQnVector();
  Double_t ExtractQnVector();
  const AliQnCorrectionsQnVector *GetQnVectorFromList(const TList *qnlist, const char* subdetector, const char *expcorr, const char *altcorr);

 private:
  AliVEvent       *fEvent;    //! ESD object
  AliESDEvent     *fESD;    //! ESD object
  AliAODEvent     *fAOD;    //! AOD object
  AliAODHeader    *fHeader; //! AOD header
  AliPIDResponse  *fPIDResponse;
  UInt_t          ftrigBit;
  TString fEstimator;//V0[M|A|C], ZN[A|C], CL[0|1]
  AliMultSelection *fMultSelection;
  Float_t fCentralityMain;
  Float_t fCentralityMin;
  Float_t fCentralityMax;

  AliEventPoolManager *fPoolManagerlambda;
  AliEventPoolManager *fPoolManagerantilambda;
  AliEventPoolManager *fPoolManagerxi;
  AliEventPoolManager *fPoolManagerxip;
  AliEventPoolManager *fPoolManagerproton;
  AliEventPoolManager *fPoolManagerantiproton;
  AliEventPoolManager *fPoolManagerlam;
  AliEventPoolManager *fPoolManagerantilam;
  AliEventPoolManager *fPoolManagerx;
  AliEventPoolManager *fPoolManagerxp;
  AliEventPoolManager *fPoolManageromega;
  AliEventPoolManager *fPoolManageromegap;
  AliEventPoolManager *fPoolManagerp;
  AliEventPoolManager *fPoolManagerantip;

  Bool_t fIsFlowTask;
  AliQnCorrectionsManager *fFlowQnVectorMgr;
  THashList             *fOutputList; //! Output list

  AliAnalysisTaskMSDibaryons(const AliAnalysisTaskMSDibaryons&); // not implemented
  AliAnalysisTaskMSDibaryons& operator=(const AliAnalysisTaskMSDibaryons&); // not implemented

  ClassDef(AliAnalysisTaskMSDibaryons,2);

};

#endif
