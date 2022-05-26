#ifndef ALIANALYSISTASKMSDIBARYONS_H
#define ALIANALYSISTASKMSDIBARYONS_H

// ROOT includes
#include <TList.h>
#include <TH1.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisUtils.h>
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliEventPoolManager.h>

class AliAnalysisTaskMSDibaryons : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMSDibaryons();
  AliAnalysisTaskMSDibaryons(const char *name);
  virtual ~AliAnalysisTaskMSDibaryons();
  
  virtual void  SetTrigger(UInt_t ktriggerInt=AliVEvent::kINT7){ ftrigBit=ktriggerInt; }
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);

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

 private:
  AliAODEvent       *fAOD;    //! AOD object
  AliAODHeader      *fHeader; //! AOD header
  AliAnalysisUtils  *fUtils;
  AliPIDResponse    *fPIDResponse;

  UInt_t            ftrigBit;

  Int_t             nevt;
  Int_t             lambdacounter;
  Int_t             eventdepth[5][10];
  Int_t             eventdepth2[5][10];
  Int_t             mnLambda[5][10][5];
  Double_t          mEnergy[5][10][5][10];
  Double_t          mPx[5][10][5][10];
  Double_t          mPy[5][10][5][10];
  Double_t          mPz[5][10][5][10];  

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

  TList             *fOutputList; //! Output list

  AliAnalysisTaskMSDibaryons(const AliAnalysisTaskMSDibaryons&); // not implemented
  AliAnalysisTaskMSDibaryons& operator=(const AliAnalysisTaskMSDibaryons&); // not implemented

  ClassDef(AliAnalysisTaskMSDibaryons,1);

};

#endif
