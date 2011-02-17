/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial pures is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any pure. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// AliAnalysisTaskPIDflowQA:
// QA for pid
//
//origin: Marek Chojnacki, Marek.Chojnacki@cern.ch
//modified: Mikolaj Krzewicki, Mikolaj.Krzewicki@cern.ch

#include "AliAnalysisTaskPIDflowQA.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "AliVEvent.h"
#include "AliESDtrackCuts.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliFlowEventCuts.h"

ClassImp( AliAnalysisTaskPIDflowQA)

//________________________________________________________________________
AliAnalysisTaskPIDflowQA:: AliAnalysisTaskPIDflowQA():
  AliAnalysisTaskSE("AliAnalysisTaskPIDflowQA"),
  fESD(NULL),
  fCuts(NULL),
  fEventCuts(NULL),
  fESDpid(NULL),
  fMC(kFALSE),
  fITSsignal(NULL),
  fTPCsignal(NULL),
  fTOFsignal(NULL),
  fITSsignalpi(NULL),
  fTPCsignalpi(NULL),
  fTOFsignalpi(NULL),
  fITSsignalK(NULL),
  fTPCsignalK(NULL),
  fTOFsignalK(NULL),
  fITSsignalp(NULL),
  fTPCsignalp(NULL),
  fTOFsignalp(NULL),
  fITSsignalpiMC(NULL),
  fTPCsignalpiMC(NULL),
  fTOFsignalpiMC(NULL),
  fITSsignalKMC(NULL),
  fTPCsignalKMC(NULL),
  fTOFsignalKMC(NULL),
  fITSsignalpMC(NULL),
  fTPCsignalpMC(NULL),
  fTOFsignalpMC(NULL),
  fITSsignalpip(NULL),
  fTPCsignalpip(NULL),
  fTOFsignalpip(NULL),
  fITSsignalKp(NULL),
  fTPCsignalKp(NULL),
  fTOFsignalKp(NULL),
  fITSsignalpp(NULL),
  fTPCsignalpp(NULL),
  fTOFsignalpp(NULL),
  fITSsignalpiMCp(NULL),
  fTPCsignalpiMCp(NULL),
  fTOFsignalpiMCp(NULL),
  fITSsignalKMCp(NULL),
  fTPCsignalKMCp(NULL),
  fTOFsignalKMCp(NULL),
  fITSsignalpMCp(NULL),
  fTPCsignalpMCp(NULL),
  fTOFsignalpMCp(NULL),
  fTOFsignalPiExpKvsPt(NULL),
  fTOFsignalPiExpPvsPt(NULL),
  fTOFsignalKExpPivsPt(NULL),
  fTOFsignalKExpPvsPt(NULL),
  fTOFsignalPExpPivsPt(NULL),
  fTOFsignalPExpKvsPt(NULL),
  fTOFsignalPiExpKvsP(NULL),
  fTOFsignalPiExpPvsP(NULL),
  fTOFsignalKExpPivsP(NULL),
  fTOFsignalKExpPvsP(NULL),
  fTOFsignalPExpPivsP(NULL),
  fTOFsignalPExpKvsP(NULL),
  fTOFsignalBeta(NULL),
  fTOFsignalPiBeta(NULL),
  fTOFsignalKBeta(NULL),
  fTOFsignalPBeta(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fMeanPtvsPt(NULL),
  fOutputList(NULL)
{
  //def ctor
}

//________________________________________________________________________
AliAnalysisTaskPIDflowQA:: AliAnalysisTaskPIDflowQA(const char *name):
  AliAnalysisTaskSE(name),
  fESD(NULL),
  fCuts(NULL),
  fEventCuts(NULL),
  fESDpid(NULL),
  fMC(kFALSE),
  fITSsignal(NULL),
  fTPCsignal(NULL),
  fTOFsignal(NULL),
  fITSsignalpi(NULL),
  fTPCsignalpi(NULL),
  fTOFsignalpi(NULL),
  fITSsignalK(NULL),
  fTPCsignalK(NULL),
  fTOFsignalK(NULL),
  fITSsignalp(NULL),
  fTPCsignalp(NULL),
  fTOFsignalp(NULL),
  fITSsignalpiMC(NULL),
  fTPCsignalpiMC(NULL),
  fTOFsignalpiMC(NULL),
  fITSsignalKMC(NULL),
  fTPCsignalKMC(NULL),
  fTOFsignalKMC(NULL),
  fITSsignalpMC(NULL),
  fTPCsignalpMC(NULL),
  fTOFsignalpMC(NULL),
  fITSsignalpip(NULL),
  fTPCsignalpip(NULL),
  fTOFsignalpip(NULL),
  fITSsignalKp(NULL),
  fTPCsignalKp(NULL),
  fTOFsignalKp(NULL),
  fITSsignalpp(NULL),
  fTPCsignalpp(NULL),
  fTOFsignalpp(NULL),
  fITSsignalpiMCp(NULL),
  fTPCsignalpiMCp(NULL),
  fTOFsignalpiMCp(NULL),
  fITSsignalKMCp(NULL),
  fTPCsignalKMCp(NULL),
  fTOFsignalKMCp(NULL),
  fITSsignalpMCp(NULL),
  fTPCsignalpMCp(NULL),
  fTOFsignalpMCp(NULL),
  fTOFsignalPiExpKvsPt(NULL),
  fTOFsignalPiExpPvsPt(NULL),
  fTOFsignalKExpPivsPt(NULL),
  fTOFsignalKExpPvsPt(NULL),
  fTOFsignalPExpPivsPt(NULL),
  fTOFsignalPExpKvsPt(NULL),
  fTOFsignalPiExpKvsP(NULL),
  fTOFsignalPiExpPvsP(NULL),
  fTOFsignalKExpPivsP(NULL),
  fTOFsignalKExpPvsP(NULL),
  fTOFsignalPExpPivsP(NULL),
  fTOFsignalPExpKvsP(NULL),
  fTOFsignalBeta(NULL),
  fTOFsignalPiBeta(NULL),
  fTOFsignalKBeta(NULL),
  fTOFsignalPBeta(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fMeanPtvsPt(NULL),
  fOutputList(NULL)
{
  //Constructor
  fESDpid=new AliESDpid();

  //old
  fESDpid->GetTPCResponse().SetBetheBlochParameters(0.0283086,
      2.63394e+01,
      5.04114e-11,
      2.12543e+00,
      4.88663e+00 );
  //new
  //fESDpid->GetTPCResponse().SetBetheBlochParameters(1.28949/50.,
  //                                                  2.74095e+01,
  //                                                  TMath::Exp(-3.21763e+01),
  //                                                  2.44026,
  //                                                  6.58800);

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::UserCreateOutputObjects()
{
  //UserCreateOutputObject
  if (fOutputList) fOutputList->Delete();
  delete fOutputList;
  fOutputList=new TList();

  const  Int_t ndec=2;
  Int_t startvalue=-1;
  const  Int_t npredec=50;
  Double_t tabx[ndec*npredec+1];
  for (Int_t i=0; i<ndec; i++)
  {
    for (Int_t j=0; j<npredec; j++)
    {
      tabx[npredec*i+j]=TMath::Power(10,((Double_t)i)+((Double_t)startvalue)+((Double_t)j)/((Double_t)npredec));
    }
  }
  tabx[ndec*npredec]=TMath::Power(10,ndec+startvalue);


  fITSsignal=new TH2F("fITSsignal","fITSsignal;dEdx;p[GeV/c]",ndec*npredec,tabx,900,0,900);
  fOutputList->Add(fITSsignal);
  fTPCsignal=new TH2F("fTPCsignal","fTPCsignal;dEdx;p[GeV/c]",ndec*npredec,tabx,900,0,900);
  fOutputList->Add(fTPCsignal);
  fTOFsignal=new TH2F("fTOFsignal","fTOFsignal;t-t_{#pi};p[GeV/c]",ndec*npredec,tabx,1200,-2000,10000);
  fOutputList->Add(fTOFsignal);

  Int_t kPtBins=60;
  Double_t binsPtDummy[kPtBins+1];
  binsPtDummy[0]=0.0;
  for(int i=1; i<=kPtBins+1; i++)
  {
    if(binsPtDummy[i-1]+0.05<1.01)
      binsPtDummy[i]=binsPtDummy[i-1]+0.05;
    else
      binsPtDummy[i]=binsPtDummy[i-1]+0.1;
  }

  Int_t kPBins=60;
  Double_t binsPDummy[kPBins+1];
  binsPDummy[0]=0.0;
  for(int i=1; i<=kPBins+1; i++)
  {
    if(binsPDummy[i-1]+0.05<1.01)
      binsPDummy[i]=binsPDummy[i-1]+0.05;
    else
      binsPDummy[i]=binsPDummy[i-1]+0.1;
  }

  fITSsignalpi=new TH2F("fITSsignalpi",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for pi+
  fTPCsignalpi=new TH2F("fTPCsignalpi",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for pi+
  fTOFsignalpi=new TH2F("fTOFsignalpi",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for pi+
  fOutputList->Add(fITSsignalpi);
  fOutputList->Add(fTPCsignalpi);
  fOutputList->Add(fTOFsignalpi);

  fITSsignalK=new TH2F("fITSsignalK",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for K+
  fTPCsignalK=new TH2F("fTPCsignalK",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for K+
  fTOFsignalK=new TH2F("fTOFsignalK",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for K+
  fOutputList->Add(fITSsignalK);
  fOutputList->Add(fTPCsignalK);
  fOutputList->Add(fTOFsignalK);

  fITSsignalp=new TH2F("fITSsignalp",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for p
  fTPCsignalp=new TH2F("fTPCsignalp",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for p
  fTOFsignalp=new TH2F("fTOFsignalp",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for p
  fOutputList->Add(fITSsignalp);
  fOutputList->Add(fTPCsignalp);
  fOutputList->Add(fTOFsignalp);

  fTOFsignalPiExpKvsPt=new TH2F("fTOFsignalPiExpKvsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPiExpPvsPt=new TH2F("fTOFsignalPiExpPvsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalKExpPivsPt=new TH2F("fTOFsignalKExpPivsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalKExpPvsPt=new TH2F("fTOFsignalKExpPvsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPExpPivsPt=new TH2F("fTOFsignalPExpPivsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPExpKvsPt=new TH2F("fTOFsignalPExpKvsPt",";pt[GeV/c];expected signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fOutputList->Add(fTOFsignalPiExpKvsPt);
  fOutputList->Add(fTOFsignalPiExpPvsPt);
  fOutputList->Add(fTOFsignalKExpPivsPt);
  fOutputList->Add(fTOFsignalKExpPvsPt);
  fOutputList->Add(fTOFsignalPExpPivsPt);
  fOutputList->Add(fTOFsignalPExpKvsPt);

  //p
  fITSsignalpip=new TH2F("fITSsignalpip",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for pi+
  fTPCsignalpip=new TH2F("fTPCsignalpip",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//TPC PID signal as function of p for pi+
  fTOFsignalpip=new TH2F("fTOFsignalpip",";p[GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for pi+
  fOutputList->Add(fITSsignalpip);
  fOutputList->Add(fTPCsignalpip);
  fOutputList->Add(fTOFsignalpip);

  fITSsignalKp=new TH2F("fITSsignalKp",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for K+
  fTPCsignalKp=new TH2F("fTPCsignalKp",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//TPC PID signal as function of p for K+
  fTOFsignalKp=new TH2F("fTOFsignalKp",";p[GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for K+
  fOutputList->Add(fITSsignalKp);
  fOutputList->Add(fTPCsignalKp);
  fOutputList->Add(fTOFsignalKp);

  fITSsignalpp=new TH2F("fITSsignalpp",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//ITS PID signal as function of p for p
  fTPCsignalpp=new TH2F("fTPCsignalpp",";p[GeV/c];signal",kPBins,binsPDummy,600,-4,4);//TPC PID signal as function of p for p
  fTOFsignalpp=new TH2F("fTOFsignalpp",";p[GeV/c];signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of p for p
  fOutputList->Add(fITSsignalpp);
  fOutputList->Add(fTPCsignalpp);
  fOutputList->Add(fTOFsignalpp);

  fTOFsignalPiExpKvsP=new TH2F("fTOFsignalPiExpKvsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPiExpPvsP=new TH2F("fTOFsignalPiExpPvsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalKExpPivsP=new TH2F("fTOFsignalKExpPivsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalKExpPvsP=new TH2F("fTOFsignalKExpPvsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPExpPivsP=new TH2F("fTOFsignalPExpPivsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fTOFsignalPExpKvsP=new TH2F("fTOFsignalPExpKvsP",";p[GeV/c];expected signal",kPBins,binsPDummy,1000,-8000,8000);//TOF PID signal as function of pt
  fOutputList->Add(fTOFsignalPiExpKvsP);
  fOutputList->Add(fTOFsignalPiExpPvsP);
  fOutputList->Add(fTOFsignalKExpPivsP);
  fOutputList->Add(fTOFsignalKExpPvsP);
  fOutputList->Add(fTOFsignalPExpPivsP);
  fOutputList->Add(fTOFsignalPExpKvsP);

  fTOFsignalBeta=new TH2F("fTOFsignalBeta",";p[GeV/c];#beta",kPBins,binsPDummy,1000, 0.2, 1.1);//
  fTOFsignalPiBeta=new TH2F("fTOFsignalPiBeta",";p[GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,1000, -1.0, 1.0);//
  fTOFsignalKBeta=new TH2F("fTOFsignalKBeta",";p[GeV/c];#beta-#beta_{K}",kPBins,binsPDummy,1000, -1.0, 1.0);//
  fTOFsignalPBeta=new TH2F("fTOFsignalPBeta",";p[GeV/c];#beta-#beta_{p}",kPBins,binsPDummy,1000, -1.0, 1.0);//
  fOutputList->Add(fTOFsignalBeta);
  fOutputList->Add(fTOFsignalPiBeta);
  fOutputList->Add(fTOFsignalKBeta);
  fOutputList->Add(fTOFsignalPBeta);

  if(fMC)
  {
    fITSsignalpiMC=new TH2F("fITSsignalpiMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for pi+
    fTPCsignalpiMC=new TH2F("fTPCsignalpiMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for pi+
    fTOFsignalpiMC=new TH2F("fTOFsignalpiMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for pi+
    fOutputList->Add(fITSsignalpiMC);
    fOutputList->Add(fTPCsignalpiMC);
    fOutputList->Add(fTOFsignalpiMC);

    fITSsignalKMC=new TH2F("fITSsignalKMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for K+
    fTPCsignalKMC=new TH2F("fTPCsignalKMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for K+
    fTOFsignalKMC=new TH2F("fTOFsignalKMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for K+
    fOutputList->Add(fITSsignalKMC);
    fOutputList->Add(fTPCsignalKMC);
    fOutputList->Add(fTOFsignalKMC);

    fITSsignalpMC=new TH2F("fITSsignalpMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for p
    fTPCsignalpMC=new TH2F("fTPCsignalpMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for p
    fTOFsignalpMC=new TH2F("fTOFsignalpMC",";pt[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for p
    fOutputList->Add(fITSsignalpMC);
    fOutputList->Add(fTPCsignalpMC);
    fOutputList->Add(fTOFsignalpMC);

    fITSsignalpiMCp=new TH2F("fITSsignalpiMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for pi+
    fTPCsignalpiMCp=new TH2F("fTPCsignalpiMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for pi+
    fTOFsignalpiMCp=new TH2F("fTOFsignalpiMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for pi+
    fOutputList->Add(fITSsignalpiMCp);
    fOutputList->Add(fTPCsignalpiMCp);
    fOutputList->Add(fTOFsignalpiMCp);

    fITSsignalKMCp=new TH2F("fITSsignalKMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for K+
    fTPCsignalKMCp=new TH2F("fTPCsignalKMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for K+
    fTOFsignalKMCp=new TH2F("fTOFsignalKMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for K+
    fOutputList->Add(fITSsignalKMCp);
    fOutputList->Add(fTPCsignalKMCp);
    fOutputList->Add(fTOFsignalKMCp);

    fITSsignalpMCp=new TH2F("fITSsignalpMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//ITS PID signal as function of pt for p
    fTPCsignalpMCp=new TH2F("fTPCsignalpMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,600,-4,4);//TPC PID signal as function of pt for p
    fTOFsignalpMCp=new TH2F("fTOFsignalpMCp",";p[GeV/c];signal",kPtBins,binsPtDummy,1000,-8000,8000);//TOF PID signal as function of pt for p
    fOutputList->Add(fITSsignalpMCp);
    fOutputList->Add(fTPCsignalpMCp);
    fOutputList->Add(fTOFsignalpMCp);
  }

  fPvsPt=new TH2F("fPvsPt","p vs p_{t}",kPBins,binsPDummy,kPtBins,binsPtDummy);
  fOutputList->Add(fPvsPt);

  fMeanPvsP = new TProfile("fMeanPvsP","Mean P vs P",kPBins,binsPDummy);
  fMeanPtvsPt = new TProfile("fMeanPtvsPt","Mean Pt vs Pt",kPtBins,binsPtDummy);
  fOutputList->Add(fMeanPvsP);
  fOutputList->Add(fMeanPtvsPt);

  PostData(1,  fOutputList);
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::UserExec(Option_t *)
{
  fESD = dynamic_cast<AliESDEvent*> (InputEvent());
  if (!fESD) return;

  //do the calibration bit
  fESDpid->SetTOFResponse(fESD,AliESDpid::kTOF_T0); // to use T0-TOF 
  fESDpid->MakePID(fESD,kFALSE);

  if(!fCuts || !fEventCuts)
  {
    Printf("No CUTS Defined.........\n");
    PostData(1,  fOutputList);
    return;
  }

  if (!(fEventCuts->IsSelected(fESD)))
  {
    return;
  }

  AliStack* stack=0x0;
  if(fMC)
  {
    AliMCEvent* mcEvent  = (AliMCEvent*) MCEvent();
    Printf("MC particles: %d", mcEvent->GetNumberOfTracks());
    stack = mcEvent->Stack();
  }

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Int_t nTracks=fESD->GetNumberOfTracks();

  AliESDtrack *trackESD=0;

  for(int tr1=0; tr1<nTracks; tr1++)
  {
    trackESD=fESD->GetTrack(tr1);
    if (!trackESD) continue;

    if(!(fCuts->AcceptTrack(trackESD))) continue;

    Int_t label=-1;
    if(fMC) label=trackESD->GetLabel();

    Int_t pdgcode=0;
    if(stack&&fMC)
    {
      TParticle* particle2 = stack->Particle(TMath::Abs(label));
      pdgcode=particle2->GetPdgCode();
    }

    Double_t p=trackESD->GetP();
    Double_t pt=trackESD->Pt();
    fPvsPt->Fill(p,pt);
    fMeanPvsP->Fill(p,p);
    fMeanPtvsPt->Fill(pt,pt);

    pidITS(trackESD,pdgcode);
    pidTPC(trackESD,pdgcode);
    pidTOF(trackESD,pdgcode);
  }

  // Post output data.
  PostData(1,  fOutputList);
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::Terminate(Option_t *)
{
  //Terminate
  if(fCuts)
    fCuts->Dump();
  if(fMC)
    Printf("MC On\n");

  Printf("AliAnalysisTaskPIDflowQA: end of Terminate");
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidITS(AliESDtrack* t, Int_t pdgcode)
{
  Int_t ngoodSDDSSD=0;
  Double_t sample[4]= {0.0,0.0,0.0,0.0};
  t->GetITSdEdxSamples(sample);
  for(int i=0; i<4; i++)
  {
    if(sample[i]>50.0)
      ngoodSDDSSD++;
  }
  if(ngoodSDDSSD<3) return;

  Float_t dedx=(Float_t)t->GetITSsignal();
  if(dedx<1.0) return;

  Bool_t ifSA=!(t->IsOn(AliESDtrack::kTPCin));
  Float_t p=t->GetP();
  Float_t pt=t->Pt();

  Float_t signalpi=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(2),ifSA);
  Float_t signalK=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(3),ifSA);
  Float_t signalp=fESDpid->GetITSResponse().Bethe(p,AliPID::ParticleMass(4),ifSA);
  if(signalpi<1.0||signalK<1.0||signalp<1.0) return;

  fITSsignal->Fill(p,dedx);

  fITSsignalpip->Fill(p,TMath::Log(dedx)-TMath::Log(signalpi));
  fITSsignalKp->Fill(p,TMath::Log(dedx)-TMath::Log(signalK));
  fITSsignalpp->Fill(p,TMath::Log(dedx)-TMath::Log(signalp));

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fITSsignalpiMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalpi));
    else if(TMath::Abs(pdgcode)==321)
      fITSsignalKMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalK));
    else if (TMath::Abs(pdgcode)==2212)
      fITSsignalpMCp->Fill(p,TMath::Log(dedx)-TMath::Log(signalp));
  }

  fITSsignalpi->Fill(pt,TMath::Log(dedx)-TMath::Log(signalpi));
  fITSsignalK->Fill(pt,TMath::Log(dedx)-TMath::Log(signalK));
  fITSsignalp->Fill(pt,TMath::Log(dedx)-TMath::Log(signalp));

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fITSsignalpiMC->Fill(pt,TMath::Log(dedx)-TMath::Log(signalpi));
    else if(TMath::Abs(pdgcode)==321)
      fITSsignalKMC->Fill(pt,TMath::Log(dedx)-TMath::Log(signalK));
    else if (TMath::Abs(pdgcode)==2212)
      fITSsignalpMC->Fill(pt,TMath::Log(dedx)-TMath::Log(signalp));
  }
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTPC(AliESDtrack* t, Int_t pdgcode)
{
  Double_t pinTPCglobal=t->GetInnerParam()->GetP();
  Float_t sigPion     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kPion);
  Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
  Float_t sigProton   = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
  Float_t pt=t->Pt();
  Float_t p=t->GetP();
  Double_t tpcSignal =t ->GetTPCsignal();
  if(!(sigPion>0.0&&sigKaon>0.0&&sigProton>0.0))
    return;

  fTPCsignal->Fill(pinTPCglobal,tpcSignal);

  fTPCsignalpip->Fill(p,(tpcSignal-sigPion)/sigPion);
  fTPCsignalKp->Fill(p,(tpcSignal-sigKaon)/sigKaon);
  fTPCsignalpp->Fill(p,(tpcSignal-sigProton)/sigProton);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTPCsignalpiMCp->Fill(p,(tpcSignal-sigPion)/sigPion);
    else if(TMath::Abs(pdgcode)==321)
      fTPCsignalKMCp->Fill(p,(tpcSignal-sigKaon)/sigKaon);
    else if (TMath::Abs(pdgcode)==2212)
      fTPCsignalpMCp->Fill(p,(tpcSignal-sigProton)/sigProton);
  }

  fTPCsignalpi->Fill(pt,(tpcSignal-sigPion)/sigPion);
  fTPCsignalK->Fill(pt,(tpcSignal-sigKaon)/sigKaon);
  fTPCsignalp->Fill(pt,(tpcSignal-sigProton)/sigProton);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTPCsignalpiMC->Fill(pt,(tpcSignal-sigPion)/sigPion);
    else if(TMath::Abs(pdgcode)==321)
      fTPCsignalKMC->Fill(pt,(tpcSignal-sigKaon)/sigKaon);
    else if (TMath::Abs(pdgcode)==2212)
      fTPCsignalpMC->Fill(pt,(tpcSignal-sigProton)/sigProton);
  }
}

//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTOF(AliESDtrack* t, Int_t pdgcode)
{
  Bool_t goodtrack = (t) &&
                     (t->GetStatus() & AliESDtrack::kTOFpid) &&
                     (t->GetTOFsignal() > 12000) &&
                     (t->GetTOFsignal() < 100000) &&
                     (t->GetIntegratedLength() > 365) &&
                     !(t->GetStatus() & AliESDtrack::kTOFmismatch);

  if (!goodtrack) return;

  const Float_t c = 2.99792457999999984e-02;
  Float_t pt=t->Pt();
  Float_t p=t->GetP();
  Float_t L = t->GetIntegratedLength();
  Float_t fT0track=fESDpid->GetTOFResponse().GetStartTime(p);
  Float_t timeTOF=t->GetTOFsignal()- fT0track;

  //calculate beta for the track
  Float_t beta = L/timeTOF/c;
  
  //2=pion 3=kaon 4=protons
  Double_t inttimes[5]= {-1.0,-1.0,-1.0,-1.0,-1.0};
  t->GetIntegratedTimes(inttimes);
  Float_t betaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    betaHypothesis[i] = L/inttimes[i]/c;
  }

  fTOFsignal->Fill(p,timeTOF-inttimes[2]);

  //beta part
  fTOFsignalBeta->Fill(p,beta);
  fTOFsignalPiBeta->Fill(p,beta-betaHypothesis[2]);
  fTOFsignalKBeta->Fill(p,beta-betaHypothesis[3]);
  fTOFsignalPBeta->Fill(p,beta-betaHypothesis[4]);

  //P part
  fTOFsignalpip->Fill(p,timeTOF-inttimes[2]);
  fTOFsignalKp->Fill(p,timeTOF-inttimes[3]);
  fTOFsignalpp->Fill(p,timeTOF-inttimes[4]);

  fTOFsignalPiExpKvsPt->Fill(pt,-inttimes[2]+inttimes[3]);
  fTOFsignalPiExpPvsPt->Fill(pt,-inttimes[2]+inttimes[4]);
  fTOFsignalKExpPivsPt->Fill(pt,-inttimes[3]+inttimes[2]);
  fTOFsignalKExpPvsPt->Fill(pt,-inttimes[3]+inttimes[4]);
  fTOFsignalPExpPivsPt->Fill(pt,-inttimes[4]+inttimes[2]);
  fTOFsignalPExpKvsPt->Fill(pt,-inttimes[4]+inttimes[3]);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTOFsignalpiMCp->Fill(p,timeTOF-inttimes[2]);
    else if(TMath::Abs(pdgcode)==321)
      fTOFsignalKMCp->Fill(p,timeTOF-inttimes[3]);
    else if (TMath::Abs(pdgcode)==2212)
      fTOFsignalpMCp->Fill(p,timeTOF-inttimes[4]);
  }

  //Pt part
  fTOFsignalpi->Fill(pt,timeTOF-inttimes[2]);
  fTOFsignalK->Fill(pt,timeTOF-inttimes[3]);
  fTOFsignalp->Fill(pt,timeTOF-inttimes[4]);

  fTOFsignalPiExpKvsP->Fill(p,-inttimes[2]+inttimes[3]);
  fTOFsignalPiExpPvsP->Fill(p,-inttimes[2]+inttimes[4]);
  fTOFsignalKExpPivsP->Fill(p,-inttimes[3]+inttimes[2]);
  fTOFsignalKExpPvsP->Fill(p,-inttimes[3]+inttimes[4]);
  fTOFsignalPExpPivsP->Fill(p,-inttimes[4]+inttimes[2]);
  fTOFsignalPExpKvsP->Fill(p,-inttimes[4]+inttimes[3]);

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTOFsignalpiMC->Fill(pt,timeTOF-inttimes[2]);
    else if(TMath::Abs(pdgcode)==321)
      fTOFsignalKMC->Fill(pt,timeTOF-inttimes[3]);
    else if (TMath::Abs(pdgcode)==2212)
      fTOFsignalpMC->Fill(pt,timeTOF-inttimes[4]);
  }
}

Float_t AliAnalysisTaskPIDflowQA::Beta(Float_t m, Float_t p) 
{
  //get theoretical beta
  return TMath::Sqrt(1. / (1. + m * m / (p * p)));
}
 
