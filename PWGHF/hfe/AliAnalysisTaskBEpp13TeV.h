/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskBEpp13TeV_H
#define AliAnalysisTaskBEpp13TeV_H

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliHFEextraCuts;

class AliAnalysisTaskBEpp13TeV : public AliAnalysisTaskSE  
{
  public:
                            AliAnalysisTaskBEpp13TeV();
                            AliAnalysisTaskBEpp13TeV(const char *name);
	virtual                 ~AliAnalysisTaskBEpp13TeV();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);

	//Setter
	void SetMCanalysis() { fIsMC = true; }
	void SetMinTPCnCrossedRow(int minTPCncls) { fMinTPCnCrossedRow = minTPCncls; }
	void SetMinTPCNclsPID(int minTPCnclsPID) { fMinTPCNclsPID = minTPCnclsPID; }
	void SetMaxTPCchi2(double maxTPCchi2) { fMaxTPCchi2 = maxTPCchi2; }
	void SetMinTPCclsRatio(double minTPCclsratio) { fMinTPCclsRatio = minTPCclsratio; }
	void SetMinITSNcls(int minITSncls) { fMinITSNcls = minITSncls; }
	void SetITSlayer(int itsLayer) { fITSlayer = itsLayer; }
	void SetPIDCuts(double tpcPIDlow, double tpcPIDhigh, double tofPID) { fTPCnsigmaLow = tpcPIDlow;fTPCnsigmaHigh = tpcPIDhigh; fTOFnsigma = tofPID; }

  private:

	bool PassEventCuts(AliAODEvent *event);
	bool PassPileUpEvent(AliAODEvent *event);
	bool PassTrackCuts(AliAODTrack *track);

	bool		fIsMC;
	int			fMinTPCnCrossedRow;
	int			fMinTPCNclsPID;
	double		fMaxTPCchi2;
	double		fMinTPCclsRatio;
	int			fMinITSNcls;
	int			fITSlayer;
	double		fTPCnsigmaLow;
	double		fTPCnsigmaHigh;
	double		fTOFnsigma;

    AliAODEvent				*fAOD;           //! input event
    TList					*fOutputList;    //! output list
    AliPIDResponse			*fPIDResponse;
	AliHFEextraCuts			*fExtraCuts;
	
	TH1F	*fHistPt;        //! dummy histogram
	TH1F	*hNrEvents;
	
	// track cut QA
	TH1F	*hFilterMask;
	TH1F	*hTPCnCrossedRow;
	TH1F	*hTPCclsPID;
	TH1F	*hTPCchi2;
	TH1F	*hTPCclsRatio;
	TH1F	*hITSNcls;
	TH1F	*hITSlayer;
	TH1F	*hDCAxy;
	TH1F	*hDCAz;
	TH1F	*hPt;
	TH1F	*hEta;
	TH1F	*hPhi;

	// pid cut
	TH2F	*hTPCnsigma;
	TH2F	*hTPCnsigmaTOFcut;
	TH2F	*hTPCnsigmaTOFcutPt;
	TH2F	*hTPCnsigmaQA;
	TH2F	*hTPCnsigmaPiQA;
	TH2F	*hTOFnsigma;
	TH2F	*hTOFnsigmaQA;

	// dca
	TH2F	*dcaTrack;
	TH2F	*dcaPion;
	
    AliAnalysisTaskBEpp13TeV(const AliAnalysisTaskBEpp13TeV&); // not implemented
    AliAnalysisTaskBEpp13TeV& operator=(const AliAnalysisTaskBEpp13TeV&); // not implemented

    ClassDef(AliAnalysisTaskBEpp13TeV, 1);
};

#endif
