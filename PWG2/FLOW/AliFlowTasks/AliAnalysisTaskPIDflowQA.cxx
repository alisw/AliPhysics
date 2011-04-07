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

#include <stdio.h>

#include "AliAnalysisTaskPIDflowQA.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TMath.h"
#include "AliVEvent.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliFlowEventCuts.h"
#include "AliFlowTrackCuts.h"
#include "AliVEventHandler.h"
#include "AliInputEventHandler.h"
#include "TTree.h"
#include "TFile.h"

ClassImp( AliAnalysisTaskPIDflowQA)

//________________________________________________________________________
AliAnalysisTaskPIDflowQA:: AliAnalysisTaskPIDflowQA():
  AliAnalysisTaskSE("AliAnalysisTaskPIDflowQA"),
  fESD(NULL),
  fCuts(NULL),
  fEventCuts(NULL),
  fESDpid(NULL),
  fUseDebugFile(kFALSE),
  fFile(NULL),
  fTPCsignal(NULL),
  fTPCsignalPi(NULL),
  fTPCsignalK(NULL),
  fTPCsignalP(NULL),
  fTPCsignalPimc(NULL),
  fTPCsignalKmc(NULL),
  fTPCsignalPmc(NULL),
  fTOFtime(NULL),
  fTOFtimeE(NULL),
  fTOFtimePi(NULL),
  fTOFtimeK(NULL),
  fTOFtimeP(NULL),
  fTOFbeta(NULL),
  fTOFbetaE(NULL),
  fTOFbetaPi(NULL),
  fTOFbetaK(NULL),
  fTOFbetaP(NULL),
  fTOFinvbeta(NULL),
  fTOFinvbetaE(NULL),
  fTOFinvbetaPi(NULL),
  fTOFinvbetaK(NULL),
  fTOFinvbetaP(NULL),
  fTOFrawtime(NULL),
  fTOFrawtimeE(NULL),
  fTOFrawtimePi(NULL),
  fTOFrawtimeK(NULL),
  fTOFrawtimeP(NULL),
  fTOFrawbeta(NULL),
  fTOFrawbetaE(NULL),
  fTOFrawbetaPi(NULL),
  fTOFrawbetaK(NULL),
  fTOFrawbetaP(NULL),
  fTOFrawinvbeta(NULL),
  fTOFrawinvbetaE(NULL),
  fTOFrawinvbetaPi(NULL),
  fTOFrawinvbetaK(NULL),
  fTOFrawinvbetaP(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fCutsTOFbetaElectrons(NULL),
  fCutsTOFbetaPions(NULL),
  fCutsTOFbetaKaons(NULL),
  fCutsTOFbetaProtons(NULL),
  fCutsTOFbetaSimpleElectrons(NULL),
  fCutsTOFbetaSimplePions(NULL),
  fCutsTOFbetaSimpleKaons(NULL),
  fCutsTOFbetaSimpleProtons(NULL),
  fCutsTPCdedxElectrons(NULL),
  fCutsTPCdedxPions(NULL),
  fCutsTPCdedxKaons(NULL),
  fCutsTPCdedxProtons(NULL),
  fCutsTPCpidElectrons(NULL),
  fCutsTPCpidPions(NULL),
  fCutsTPCpidKaons(NULL),
  fCutsTPCpidProtons(NULL),
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
  fUseDebugFile(kFALSE),
  fFile(NULL),
  fTPCsignal(NULL),
  fTPCsignalPi(NULL),
  fTPCsignalK(NULL),
  fTPCsignalP(NULL),
  fTPCsignalPimc(NULL),
  fTPCsignalKmc(NULL),
  fTPCsignalPmc(NULL),
  fTOFtime(NULL),
  fTOFtimeE(NULL),
  fTOFtimePi(NULL),
  fTOFtimeK(NULL),
  fTOFtimeP(NULL),
  fTOFbeta(NULL),
  fTOFbetaE(NULL),
  fTOFbetaPi(NULL),
  fTOFbetaK(NULL),
  fTOFbetaP(NULL),
  fTOFinvbeta(NULL),
  fTOFinvbetaE(NULL),
  fTOFinvbetaPi(NULL),
  fTOFinvbetaK(NULL),
  fTOFinvbetaP(NULL),
  fTOFrawtime(NULL),
  fTOFrawtimeE(NULL),
  fTOFrawtimePi(NULL),
  fTOFrawtimeK(NULL),
  fTOFrawtimeP(NULL),
  fTOFrawbeta(NULL),
  fTOFrawbetaE(NULL),
  fTOFrawbetaPi(NULL),
  fTOFrawbetaK(NULL),
  fTOFrawbetaP(NULL),
  fTOFrawinvbeta(NULL),
  fTOFrawinvbetaE(NULL),
  fTOFrawinvbetaPi(NULL),
  fTOFrawinvbetaK(NULL),
  fTOFrawinvbetaP(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fCutsTOFbetaElectrons(NULL),
  fCutsTOFbetaPions(NULL),
  fCutsTOFbetaKaons(NULL),
  fCutsTOFbetaProtons(NULL),
  fCutsTOFbetaSimpleElectrons(NULL),
  fCutsTOFbetaSimplePions(NULL),
  fCutsTOFbetaSimpleKaons(NULL),
  fCutsTOFbetaSimpleProtons(NULL),
  fCutsTPCdedxElectrons(NULL),
  fCutsTPCdedxPions(NULL),
  fCutsTPCdedxKaons(NULL),
  fCutsTPCdedxProtons(NULL),
  fCutsTPCpidElectrons(NULL),
  fCutsTPCpidPions(NULL),
  fCutsTPCpidKaons(NULL),
  fCutsTPCpidProtons(NULL),
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
  fOutputList->SetOwner(kTRUE);

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

  fTPCsignal=new TH2F("fTPCsignal",";p [GeV/c];dEdx",kPBins,binsPDummy,500,0,500);
  fOutputList->Add(fTPCsignal);
  fTPCsignalPi=new TH2F("fTPCsignalPi",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for pi+
  fTPCsignalK=new TH2F("fTPCsignalK",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for K+
  fTPCsignalP=new TH2F("fTPCsignalP",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for p
  fOutputList->Add(fTPCsignalPi);
  fOutputList->Add(fTPCsignalK);
  fOutputList->Add(fTPCsignalP);

  fTOFtime=new TH2F("fTOFtime",";p[GeV/c];#time",kPBins,binsPDummy,1000, 12000, 80000);//
  fOutputList->Add(fTOFtime);
  fTOFtimeE=new TH2F("fTOFtimeE",";p [GeV/c];#time-#time_{#pi}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFtimePi=new TH2F("fTOFtimePi",";p [GeV/c];#time-#time_{#pi}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFtimeK=new TH2F("fTOFtimeK",";p [GeV/c];#time-#time_{K}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFtimeP=new TH2F("fTOFtimeP",";p [GeV/c];#time-#time_{p}",kPBins,binsPDummy,500, -8000, 8000);//
  fOutputList->Add(fTOFtimeE);
  fOutputList->Add(fTOFtimePi);
  fOutputList->Add(fTOFtimeK);
  fOutputList->Add(fTOFtimeP);

  fTOFbeta=new TH2F("fTOFbeta",";p[GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fOutputList->Add(fTOFbeta);
  fTOFbetaE=new TH2F("fTOFbetaE",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFbetaPi=new TH2F("fTOFbetaPi",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFbetaK=new TH2F("fTOFbetaK",";p [GeV/c];#beta-#beta_{K}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFbetaP=new TH2F("fTOFbetaP",";p [GeV/c];#beta-#beta_{p}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fOutputList->Add(fTOFbetaE);
  fOutputList->Add(fTOFbetaPi);
  fOutputList->Add(fTOFbetaK);
  fOutputList->Add(fTOFbetaP);

  fTOFinvbeta=new TH2F("fTOFinvbeta",";p[GeV/c];1/#beta",kPBins,binsPDummy,1000, 0.90, 2.5);//
  fOutputList->Add(fTOFinvbeta);
  fTOFinvbetaE=new TH2F("fTOFinvbetaE",";p [GeV/c];1/#beta-1/#beta_{#pi}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFinvbetaPi=new TH2F("fTOFinvbetaPi",";p [GeV/c];1/#beta-1/#beta_{#pi}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFinvbetaK=new TH2F("fTOFinvbetaK",";p [GeV/c];1/#beta-1/#beta_{K}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFinvbetaP=new TH2F("fTOFinvbetaP",";p [GeV/c];1/#beta-1/#beta_{p}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fOutputList->Add(fTOFinvbetaE);
  fOutputList->Add(fTOFinvbetaPi);
  fOutputList->Add(fTOFinvbetaK);
  fOutputList->Add(fTOFinvbetaP);

  fTOFrawtime=new TH2F("fTOFrawtime",";p[GeV/c];#time",kPBins,binsPDummy,1000, 12000, 80000);//
  fOutputList->Add(fTOFrawtime);
  fTOFrawtimeE=new TH2F("fTOFrawtimeE",";p [GeV/c];#time-#time_{#pi}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFrawtimePi=new TH2F("fTOFrawtimePi",";p [GeV/c];#time-#time_{#pi}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFrawtimeK=new TH2F("fTOFrawtimeK",";p [GeV/c];#time-#time_{K}",kPBins,binsPDummy,500, -8000, 8000);//
  fTOFrawtimeP=new TH2F("fTOFrawtimeP",";p [GeV/c];#time-#time_{p}",kPBins,binsPDummy,500, -8000, 8000);//
  fOutputList->Add(fTOFrawtimeE);
  fOutputList->Add(fTOFrawtimePi);
  fOutputList->Add(fTOFrawtimeK);
  fOutputList->Add(fTOFrawtimeP);

  fTOFrawbeta=new TH2F("fTOFrawbeta",";p[GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fOutputList->Add(fTOFrawbeta);
  fTOFrawbetaE=new TH2F("fTOFrawbetaE",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFrawbetaPi=new TH2F("fTOFrawbetaPi",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFrawbetaK=new TH2F("fTOFrawbetaK",";p [GeV/c];#beta-#beta_{K}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fTOFrawbetaP=new TH2F("fTOFrawbetaP",";p [GeV/c];#beta-#beta_{p}",kPBins,binsPDummy,500, -0.25, 0.25);//
  fOutputList->Add(fTOFrawbetaE);
  fOutputList->Add(fTOFrawbetaPi);
  fOutputList->Add(fTOFrawbetaK);
  fOutputList->Add(fTOFrawbetaP);

  fTOFrawinvbeta=new TH2F("fTOFrawinvbeta",";p[GeV/c];1/#beta",kPBins,binsPDummy,1000, 0.90, 2.5);//
  fOutputList->Add(fTOFrawinvbeta);
  fTOFrawinvbetaE=new TH2F("fTOFrawinvbetaE",";p [GeV/c];1/#beta-1/#beta_{#pi}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFrawinvbetaPi=new TH2F("fTOFrawinvbetaPi",";p [GeV/c];1/#beta-1/#beta_{#pi}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFrawinvbetaK=new TH2F("fTOFrawinvbetaK",";p [GeV/c];1/#beta-1/#beta_{K}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fTOFrawinvbetaP=new TH2F("fTOFrawinvbetaP",";p [GeV/c];1/#beta-1/#beta_{p}",kPBins,binsPDummy,600, -0.3, 0.3);//
  fOutputList->Add(fTOFrawinvbetaE);
  fOutputList->Add(fTOFrawinvbetaPi);
  fOutputList->Add(fTOFrawinvbetaK);
  fOutputList->Add(fTOFrawinvbetaP);

  fPvsPt=new TH2F("fPvsPt","p vs p_{t};p [GeV/c];p_{t} [GeV/c]",kPBins,binsPDummy,kPtBins,binsPtDummy);
  fOutputList->Add(fPvsPt);

  fMeanPvsP = new TProfile("fMeanPvsP","Mean P vs P;p [Gev/c];<p> [GeV/c]",kPBins,binsPDummy);
  fOutputList->Add(fMeanPvsP);

  fTPCvsGlobalMult = new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",500,0,2500,500,0,3500);
  fOutputList->Add(fTPCvsGlobalMult);

  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  fStandardTPCCuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();

  fCutsTOFbetaElectrons = new AliFlowTrackCuts("TOFbeta e");
  fCutsTOFbetaElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFbetaElectrons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaElectrons->SetQA();
  fCutsTOFbetaPions = new AliFlowTrackCuts("TOFbeta pi");
  fCutsTOFbetaPions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFbetaPions->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaPions->SetRejectElectronsWithTPCpid();
  fCutsTOFbetaPions->SetQA();
  fCutsTOFbetaKaons = new AliFlowTrackCuts("TOFbeta K");
  fCutsTOFbetaKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFbetaKaons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaKaons->SetQA();
  fCutsTOFbetaProtons = new AliFlowTrackCuts("TOFbeta p");
  fCutsTOFbetaProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFbetaProtons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaProtons->SetQA();

  fCutsTOFbetaSimpleElectrons = new AliFlowTrackCuts("TOFbetaSimple e");
  fCutsTOFbetaSimpleElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTOFbetaSimple);
  fCutsTOFbetaSimpleElectrons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaSimpleElectrons->SetQA();
  fCutsTOFbetaSimplePions = new AliFlowTrackCuts("TOFbetaSimple pi");
  fCutsTOFbetaSimplePions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTOFbetaSimple);
  fCutsTOFbetaSimplePions->SetRejectElectronsWithTPCpid();
  fCutsTOFbetaSimplePions->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaSimplePions->SetQA();
  fCutsTOFbetaSimpleKaons = new AliFlowTrackCuts("TOFbetaSimple K");
  fCutsTOFbetaSimpleKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTOFbetaSimple);
  fCutsTOFbetaSimpleKaons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaSimpleKaons->SetQA();
  fCutsTOFbetaSimpleProtons = new AliFlowTrackCuts("TOFbetaSimple p");
  fCutsTOFbetaSimpleProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTOFbetaSimple);
  fCutsTOFbetaSimpleProtons->SetRequireStrictTOFTPCagreement();
  fCutsTOFbetaSimpleProtons->SetQA();

  fCutsTPCdedxElectrons = new AliFlowTrackCuts("TPCdedx e");
  fCutsTPCdedxElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCdedxElectrons->SetQA();
  fCutsTPCdedxPions = new AliFlowTrackCuts("TPCdedx Pi");
  fCutsTPCdedxPions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCdedxPions->SetQA();
  fCutsTPCdedxKaons = new AliFlowTrackCuts("TPCdedx K");
  fCutsTPCdedxKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCdedxKaons->SetQA();
  fCutsTPCdedxProtons = new AliFlowTrackCuts("TPCdedx p");
  fCutsTPCdedxProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCdedxProtons->SetQA();

  fCutsTPCpidElectrons = new AliFlowTrackCuts("TPCpid e");
  fCutsTPCpidElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTPCpid);
  fCutsTPCpidElectrons->SetQA();
  fCutsTPCpidPions = new AliFlowTrackCuts("TPCpid Pi");
  fCutsTPCpidPions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTPCpid);
  fCutsTPCpidPions->SetQA();
  fCutsTPCpidKaons = new AliFlowTrackCuts("TPCpid K");
  fCutsTPCpidKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTPCpid);
  fCutsTPCpidKaons->SetQA();
  fCutsTPCpidProtons = new AliFlowTrackCuts("TPCpid p");
  fCutsTPCpidProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTPCpid);
  fCutsTPCpidProtons->SetQA();

  //fOutputList->Add(fESDpid);

  fOutputList->Add(fCutsTPCdedxElectrons->GetQA());
  fOutputList->Add(fCutsTPCdedxPions->GetQA());
  fOutputList->Add(fCutsTPCdedxKaons->GetQA());
  fOutputList->Add(fCutsTPCdedxProtons->GetQA());
  fOutputList->Add(fCutsTPCpidElectrons->GetQA());
  fOutputList->Add(fCutsTPCpidPions->GetQA());
  fOutputList->Add(fCutsTPCpidKaons->GetQA());
  fOutputList->Add(fCutsTPCpidProtons->GetQA());
  fOutputList->Add(fCutsTOFbetaElectrons->GetQA());
  fOutputList->Add(fCutsTOFbetaPions->GetQA());
  fOutputList->Add(fCutsTOFbetaKaons->GetQA());
  fOutputList->Add(fCutsTOFbetaProtons->GetQA());
  fOutputList->Add(fCutsTOFbetaSimpleElectrons->GetQA());
  fOutputList->Add(fCutsTOFbetaSimplePions->GetQA());
  fOutputList->Add(fCutsTOFbetaSimpleKaons->GetQA());
  fOutputList->Add(fCutsTOFbetaSimpleProtons->GetQA());

  if (fUseDebugFile) fFile = fopen("debug.txt","w");

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

  AliStack* stack=NULL;
  AliMCEvent* mcEvent = MCEvent();
  if (mcEvent) stack = mcEvent->Stack();
  if (mcEvent) Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

  Printf("There are %d tracks in this event", fESD->GetNumberOfTracks());
  Int_t nTracks=fESD->GetNumberOfTracks();

  AliESDtrack *trackESD=0;

  fCuts->SetEvent(fESD,mcEvent);
  fCutsTPCdedxElectrons->SetEvent(fESD,mcEvent);
  fCutsTPCdedxPions->SetEvent(fESD,mcEvent);
  fCutsTPCdedxKaons->SetEvent(fESD,mcEvent);
  fCutsTPCdedxProtons->SetEvent(fESD,mcEvent);
  fCutsTPCpidElectrons->SetEvent(fESD,mcEvent);
  fCutsTPCpidPions->SetEvent(fESD,mcEvent);
  fCutsTPCpidKaons->SetEvent(fESD,mcEvent);
  fCutsTPCpidProtons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaElectrons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaPions->SetEvent(fESD,mcEvent);
  fCutsTOFbetaKaons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaProtons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaSimpleElectrons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaSimplePions->SetEvent(fESD,mcEvent);
  fCutsTOFbetaSimpleKaons->SetEvent(fESD,mcEvent);
  fCutsTOFbetaSimpleProtons->SetEvent(fESD,mcEvent);

  for(int tr1=0; tr1<nTracks; tr1++)
  {
    trackESD=fESD->GetTrack(tr1);
    if (!trackESD) continue;

    Double_t p=trackESD->GetP();
    Double_t pt=trackESD->Pt();

    if(!(fCuts->IsSelected(trackESD))) continue;

    Int_t label=-1;
    if(mcEvent) label=trackESD->GetLabel();

    Int_t pdgcode=0;
    if(stack)
    {
      TParticle* particle2 = stack->Particle(TMath::Abs(label));
      pdgcode=particle2->GetPdgCode();
    }

    fPvsPt->Fill(p,pt);
    fMeanPvsP->Fill(p,p);

    pidTPC(trackESD,pdgcode);
    pidTOF(trackESD,pdgcode);
  }

  //check the correlation between the global and TPConly number of tracks
  fStandardGlobalCuts->SetEvent(fESD);
  fStandardTPCCuts->SetEvent(fESD);
  Int_t multGlobal = fStandardGlobalCuts->Count();
  Int_t multTPC = fStandardTPCCuts->Count();
  fTPCvsGlobalMult->Fill(multGlobal,multTPC);

  if (fFile)
  {
    const AliESDVertex* pvtx = fESD->GetPrimaryVertex();
    const AliESDVertex* tpcvtx = fESD->GetPrimaryVertexTPC();
    const AliESDVertex* spdvtx = fESD->GetPrimaryVertexSPD();
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliVEventHandler* handler = mgr->GetInputEventHandler();
    TTree* tree = handler->GetTree();
    TFile* file = tree->GetCurrentFile();
    if (multTPC>(23+1.216*multGlobal) || multTPC<(-20+1.087*multGlobal))
    {
      fprintf(fFile, "%i %i %s %i\n",multTPC,multGlobal,file->GetName(),fESD->GetEventNumberInFile());
      fprintf(fFile, "  primary vertex: x: %.2f, y: %.2f, z: %.2f, n: %i\n", pvtx->GetX(), pvtx->GetY(), pvtx->GetZ(), pvtx->GetNContributors());
      fprintf(fFile, "      SPD vertex: x: %.2f, y: %.2f, z: %.2f, n: %i\n", spdvtx->GetX(), spdvtx->GetY(), spdvtx->GetZ(), spdvtx->GetNContributors());
      fprintf(fFile, "      TPC vertex: x: %.2f, y: %.2f, z: %.2f, n: %i\n", tpcvtx->GetX(), tpcvtx->GetY(), tpcvtx->GetZ(), tpcvtx->GetNContributors());
    }
  }
}

//________________________________________________________________________
void  AliAnalysisTaskPIDflowQA::Terminate(Option_t *)
{
  //Terminate
  if(fCuts)
    fCuts->Dump();

  Printf("AliAnalysisTaskPIDflowQA: end of Terminate");
}


//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTPC(AliESDtrack* t, Int_t)
{
  //do TPC pid
  const AliExternalTrackParam* innerParam = t->GetInnerParam();
  if (!innerParam) return;
  Double_t pinTPCglobal=innerParam->GetP();
  Double_t tpcSignal =t->GetTPCsignal();
  Float_t p=innerParam->P();
  Float_t sigPion     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kPion);
  Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
  Float_t sigProton   = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
  if(!(sigPion>0.0&&sigKaon>0.0&&sigProton>0.0))
    return;

  fTPCsignal->Fill(pinTPCglobal,tpcSignal);

  fTPCsignalPi->Fill(p,(tpcSignal-sigPion)/sigPion);
  fTPCsignalK->Fill(p,(tpcSignal-sigKaon)/sigKaon);
  fTPCsignalP->Fill(p,(tpcSignal-sigProton)/sigProton);

  fCutsTPCdedxElectrons->IsSelected(t);
  fCutsTPCdedxPions->IsSelected(t);
  fCutsTPCdedxKaons->IsSelected(t);
  fCutsTPCdedxProtons->IsSelected(t);
  fCutsTPCpidElectrons->IsSelected(t);
  fCutsTPCpidPions->IsSelected(t);
  fCutsTPCpidKaons->IsSelected(t);
  fCutsTPCpidProtons->IsSelected(t);
}

//______________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTOF(AliESDtrack* track, Int_t)
{
  //do TOF pid
  Bool_t goodtrack = (track) &&
                     (track->GetStatus() & AliESDtrack::kTOFpid) &&
                     (track->GetTOFsignal() > 12000) &&
                     (track->GetTOFsignal() < 100000) &&
                     (track->GetIntegratedLength() > 365);
  
  if (!goodtrack) return;

  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t l = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid->GetTOFResponse().GetStartTime(p);
  //time
  Float_t timeTOF = track->GetTOFsignal()- trackT0; 
  Double_t integratedTimes[5] = {-1.0,-1.0,-1.0,-1.0,-1.0};
  track->GetIntegratedTimes(integratedTimes);
  //beta
  Float_t beta = l/timeTOF/c;
  Float_t betaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  Float_t betadiff[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    betaHypothesis[i] = l/integratedTimes[i]/c;
    betadiff[i] = beta-betaHypothesis[i];
  }
  //inverse beta
  Float_t invbeta = 1/beta;
  Float_t invbetaHypothesis[5] = {0.0,0.0,0.0,0.0,0.0};
  Float_t invbetadiff[5] = {0.0,0.0,0.0,0.0,0.0};
  for (Int_t i=0;i<5;i++)
  {
    invbetaHypothesis[i] = 1/betaHypothesis[i];
    invbetadiff[i] = invbeta-invbetaHypothesis[i];
  }

  Double_t tpcpid[AliPID::kSPECIES];
  track->GetTPCpid(tpcpid);

  //base hists
  fTOFrawtime->Fill(p,timeTOF);
  fTOFrawtimeE->Fill(p,timeTOF-integratedTimes[0]);
  fTOFrawtimePi->Fill(p,timeTOF-integratedTimes[2]);
  fTOFrawtimeK->Fill(p,timeTOF-integratedTimes[3]);
  fTOFrawtimeP->Fill(p,timeTOF-integratedTimes[4]);

  fTOFrawbeta->Fill(p,beta);
  fTOFrawbetaE->Fill(p,betadiff[0]);
  fTOFrawbetaPi->Fill(p,betadiff[2]);
  fTOFrawbetaK->Fill(p,betadiff[3]);
  fTOFrawbetaP->Fill(p,betadiff[4]);

  fTOFrawinvbeta->Fill(p,invbeta);
  fTOFrawinvbetaE->Fill(p,invbetadiff[0]);
  fTOFrawinvbetaPi->Fill(p,invbetadiff[2]);
  fTOFrawinvbetaK->Fill(p,invbetadiff[3]);
  fTOFrawinvbetaP->Fill(p,invbetadiff[4]);

  //cleanup with TPC
  if (track->GetStatus() & AliESDtrack::kTOFmismatch) return;
  if (!TPCTOFagree(track)) return;

  //responses
  fTOFtime->Fill(p,timeTOF);
  fTOFtimeE->Fill(p,timeTOF-integratedTimes[0]);
  fTOFtimePi->Fill(p,timeTOF-integratedTimes[2]);
  fTOFtimeK->Fill(p,timeTOF-integratedTimes[3]);
  fTOFtimeP->Fill(p,timeTOF-integratedTimes[4]);

  fTOFbeta->Fill(p,beta);
  fTOFbetaE->Fill(p,betadiff[0]);
  fTOFbetaPi->Fill(p,betadiff[2]);
  fTOFbetaK->Fill(p,betadiff[3]);
  fTOFbetaP->Fill(p,betadiff[4]);

  fTOFinvbeta->Fill(p,invbeta);
  fTOFinvbetaE->Fill(p,invbetadiff[0]);
  fTOFinvbetaPi->Fill(p,invbetadiff[2]);
  fTOFinvbetaK->Fill(p,invbetadiff[3]);
  fTOFinvbetaP->Fill(p,invbetadiff[4]);

  fCutsTOFbetaElectrons->IsSelected(track);
  fCutsTOFbetaPions->IsSelected(track);
  fCutsTOFbetaKaons->IsSelected(track);
  fCutsTOFbetaProtons->IsSelected(track);
  fCutsTOFbetaSimpleElectrons->IsSelected(track);
  fCutsTOFbetaSimplePions->IsSelected(track);
  fCutsTOFbetaSimpleKaons->IsSelected(track);
  fCutsTOFbetaSimpleProtons->IsSelected(track);
}

//______________________________________________________________________________
Float_t AliAnalysisTaskPIDflowQA::Beta(Float_t m, Float_t p) 
{
  //get theoretical beta
  return TMath::Sqrt(1. / (1. + m * m / (p * p)));
}

//---------------------------------------------------------------//
Bool_t AliAnalysisTaskPIDflowQA::TPCTOFagree(const AliESDtrack *track)
{
  Bool_t status = kFALSE;
  
  Float_t mass[5] = {5.10998909999999971e-04,1.05658000000000002e-01,1.39570000000000000e-01,4.93676999999999977e-01,9.38271999999999995e-01};
  

  Double_t exptimes[5];
  track->GetIntegratedTimes(exptimes);
  
  Float_t dedx = track->GetTPCsignal();

  Float_t p = track->P();
  Float_t time = track->GetTOFsignal()- fESDpid->GetTOFResponse().GetStartTime(p);
  Float_t tl = track->GetIntegratedLength();

  Float_t betagammares =  fESDpid->GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);

  Float_t betagamma1 = tl/(time-5 *betagammares) * 33.3564095198152043;

//  printf("betagamma1 = %f\n",betagamma1);

  if(betagamma1 < 0.1) betagamma1 = 0.1;

  if(betagamma1 < 0.99999) betagamma1 /= TMath::Sqrt(1-betagamma1*betagamma1);
  else betagamma1 = 100;

  Float_t betagamma2 = tl/(time+5 *betagammares) * 33.3564095198152043;
//  printf("betagamma2 = %f\n",betagamma2);

  if(betagamma2 < 0.1) betagamma2 = 0.1;

  if(betagamma2 < 0.99999) betagamma2 /= TMath::Sqrt(1-betagamma2*betagamma2);
  else betagamma2 = 100;


  Double_t ptpc[3];
  track->GetInnerPxPyPz(ptpc);
  Float_t momtpc=TMath::Sqrt(ptpc[0]*ptpc[0] + ptpc[1]*ptpc[1] + ptpc[2]*ptpc[2]);
 
  for(Int_t i=0;i < 5;i++){
    Float_t resolutionTOF =  fESDpid->GetTOFResponse().GetExpectedSigma(p, exptimes[i], mass[i]);
    if(TMath::Abs(exptimes[i] - time) < 5 * resolutionTOF){
      Float_t dedxExp = 0;
      if(i==0) dedxExp =  fESDpid->GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kElectron);
      else if(i==1) dedxExp =  fESDpid->GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kMuon);
      else if(i==2) dedxExp =  fESDpid->GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kPion);
      else if(i==3) dedxExp =  fESDpid->GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kKaon);
      else if(i==4) dedxExp =  fESDpid->GetTPCResponse().GetExpectedSignal(momtpc,AliPID::kProton);

      Float_t resolutionTPC = 2;
      if(i==0) resolutionTPC =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kElectron); 
      else if(i==1) resolutionTPC =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kMuon);
      else if(i==2) resolutionTPC =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kPion);
      else if(i==3) resolutionTPC =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kKaon);
      else if(i==4) resolutionTPC =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);

      if(TMath::Abs(dedx - dedxExp) < 3 * resolutionTPC){
	status = kTRUE;
      }
    }
  }

  Float_t bb1 =  fESDpid->GetTPCResponse().Bethe(betagamma1);
  Float_t bb2 =  fESDpid->GetTPCResponse().Bethe(betagamma2);
  Float_t bbM =  fESDpid->GetTPCResponse().Bethe((betagamma1+betagamma2)*0.5);


  //  status = kFALSE;
  // for nuclei
  Float_t resolutionTOFpr =   fESDpid->GetTOFResponse().GetExpectedSigma(p, exptimes[4], mass[4]);
  Float_t resolutionTPCpr =   fESDpid->GetTPCResponse().GetExpectedSigma(momtpc,track->GetTPCsignalN(),AliPID::kProton);
  if(TMath::Abs(dedx-bb1) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bb2) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  else if(TMath::Abs(dedx-bbM) < resolutionTPCpr*3 && exptimes[4] < time-7*resolutionTOFpr){
     status = kTRUE;
  }
  
  return status;
}
 
