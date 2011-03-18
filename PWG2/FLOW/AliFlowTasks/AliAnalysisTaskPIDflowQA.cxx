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
  fMC(kFALSE),
  fUseDebugFile(kFALSE),
  fFile(NULL),
  fTPCsignal(NULL),
  fTPCsignalPi(NULL),
  fTPCsignalK(NULL),
  fTPCsignalP(NULL),
  fTPCsignalPimc(NULL),
  fTPCsignalKmc(NULL),
  fTPCsignalPmc(NULL),
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
  fTOFbetaAfterElectronsCuts(NULL),
  fTOFbetaAfterPionCuts(NULL),
  fTOFbetaAfterKaonCuts(NULL),
  fTOFbetaAfterProtonCuts(NULL),
  fTPCsignalAfterPionCuts(NULL),
  fTPCsignalAfterKaonCuts(NULL),
  fTPCsignalAfterProtonCuts(NULL),
  fTOFbetaAfterElectronsCuts1(NULL),
  fTOFbetaAfterPionCuts1(NULL),
  fTOFbetaAfterKaonCuts1(NULL),
  fTOFbetaAfterProtonCuts1(NULL),
  fTPCsignalAfterPionCuts1(NULL),
  fTPCsignalAfterKaonCuts1(NULL),
  fTPCsignalAfterProtonCuts1(NULL),
  fTOFbetaEafter(NULL),
  fTOFbetaPiafter(NULL),
  fTOFbetaKafter(NULL),
  fTOFbetaPafter(NULL),
  fTPCsignalPiafter(NULL),
  fTPCsignalKafter(NULL),
  fTPCsignalPafter(NULL),
  fTOFyieldSelEmcE(NULL),
  fTOFyieldSelPimcE(NULL),
  fTOFyieldSelKmcE(NULL),
  fTOFyieldSelPmcE(NULL),
  fTOFyieldSelEmcM(NULL),
  fTOFyieldSelPimcM(NULL),
  fTOFyieldSelKmcM(NULL),
  fTOFyieldSelPmcM(NULL),
  fTOFyieldSelEmcPi(NULL),
  fTOFyieldSelPimcPi(NULL),
  fTOFyieldSelKmcPi(NULL),
  fTOFyieldSelPmcPi(NULL),
  fTOFyieldSelEmcK(NULL),
  fTOFyieldSelPimcK(NULL),
  fTOFyieldSelKmcK(NULL),
  fTOFyieldSelPmcK(NULL),
  fTOFyieldSelEmcP(NULL),
  fTOFyieldSelPimcP(NULL),
  fTOFyieldSelKmcP(NULL),
  fTOFyieldSelPmcP(NULL),
  fTPCyieldSelEmcE(NULL),
  fTPCyieldSelPimcE(NULL),
  fTPCyieldSelKmcE(NULL),
  fTPCyieldSelPmcE(NULL),
  fTPCyieldSelEmcM(NULL),
  fTPCyieldSelPimcM(NULL),
  fTPCyieldSelKmcM(NULL),
  fTPCyieldSelPmcM(NULL),
  fTPCyieldSelEmcPi(NULL),
  fTPCyieldSelPimcPi(NULL),
  fTPCyieldSelKmcPi(NULL),
  fTPCyieldSelPmcPi(NULL),
  fTPCyieldSelEmcK(NULL),
  fTPCyieldSelPimcK(NULL),
  fTPCyieldSelKmcK(NULL),
  fTPCyieldSelPmcK(NULL),
  fTPCyieldSelEmcP(NULL),
  fTPCyieldSelPimcP(NULL),
  fTPCyieldSelKmcP(NULL),
  fTPCyieldSelPmcP(NULL),
  fTPCdedxAfterTOFpidPions(NULL),
  fTPCdedxAfterTOFpidKaons(NULL),
  fTPCdedxAfterTOFpidProtons(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fCutsTOFElectrons(NULL),
  fCutsTOFPions(NULL),
  fCutsTOFKaons(NULL),
  fCutsTOFProtons(NULL),
  fCutsTPCElectrons(NULL),
  fCutsTPCPions(NULL),
  fCutsTPCKaons(NULL),
  fCutsTPCProtons(NULL),
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
  fUseDebugFile(kFALSE),
  fFile(NULL),
  fTPCsignal(NULL),
  fTPCsignalPi(NULL),
  fTPCsignalK(NULL),
  fTPCsignalP(NULL),
  fTPCsignalPimc(NULL),
  fTPCsignalKmc(NULL),
  fTPCsignalPmc(NULL),
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
  fTOFbetaAfterElectronsCuts(NULL),
  fTOFbetaAfterPionCuts(NULL),
  fTOFbetaAfterKaonCuts(NULL),
  fTOFbetaAfterProtonCuts(NULL),
  fTPCsignalAfterPionCuts(NULL),
  fTPCsignalAfterKaonCuts(NULL),
  fTPCsignalAfterProtonCuts(NULL),
  fTOFbetaAfterElectronsCuts1(NULL),
  fTOFbetaAfterPionCuts1(NULL),
  fTOFbetaAfterKaonCuts1(NULL),
  fTOFbetaAfterProtonCuts1(NULL),
  fTPCsignalAfterPionCuts1(NULL),
  fTPCsignalAfterKaonCuts1(NULL),
  fTPCsignalAfterProtonCuts1(NULL),
  fTOFbetaEafter(NULL),
  fTOFbetaPiafter(NULL),
  fTOFbetaKafter(NULL),
  fTOFbetaPafter(NULL),
  fTPCsignalPiafter(NULL),
  fTPCsignalKafter(NULL),
  fTPCsignalPafter(NULL),
  fTOFyieldSelEmcE(NULL),
  fTOFyieldSelPimcE(NULL),
  fTOFyieldSelKmcE(NULL),
  fTOFyieldSelPmcE(NULL),
  fTOFyieldSelEmcM(NULL),
  fTOFyieldSelPimcM(NULL),
  fTOFyieldSelKmcM(NULL),
  fTOFyieldSelPmcM(NULL),
  fTOFyieldSelEmcPi(NULL),
  fTOFyieldSelPimcPi(NULL),
  fTOFyieldSelKmcPi(NULL),
  fTOFyieldSelPmcPi(NULL),
  fTOFyieldSelEmcK(NULL),
  fTOFyieldSelPimcK(NULL),
  fTOFyieldSelKmcK(NULL),
  fTOFyieldSelPmcK(NULL),
  fTOFyieldSelEmcP(NULL),
  fTOFyieldSelPimcP(NULL),
  fTOFyieldSelKmcP(NULL),
  fTOFyieldSelPmcP(NULL),
  fTPCyieldSelEmcE(NULL),
  fTPCyieldSelPimcE(NULL),
  fTPCyieldSelKmcE(NULL),
  fTPCyieldSelPmcE(NULL),
  fTPCyieldSelEmcM(NULL),
  fTPCyieldSelPimcM(NULL),
  fTPCyieldSelKmcM(NULL),
  fTPCyieldSelPmcM(NULL),
  fTPCyieldSelEmcPi(NULL),
  fTPCyieldSelPimcPi(NULL),
  fTPCyieldSelKmcPi(NULL),
  fTPCyieldSelPmcPi(NULL),
  fTPCyieldSelEmcK(NULL),
  fTPCyieldSelPimcK(NULL),
  fTPCyieldSelKmcK(NULL),
  fTPCyieldSelPmcK(NULL),
  fTPCyieldSelEmcP(NULL),
  fTPCyieldSelPimcP(NULL),
  fTPCyieldSelKmcP(NULL),
  fTPCyieldSelPmcP(NULL),
  fTPCdedxAfterTOFpidPions(NULL),
  fTPCdedxAfterTOFpidKaons(NULL),
  fTPCdedxAfterTOFpidProtons(NULL),
  fPvsPt(NULL),
  fMeanPvsP(NULL),
  fTPCvsGlobalMult(NULL),
  fStandardGlobalCuts(NULL),
  fStandardTPCCuts(NULL),
  fCutsTOFElectrons(NULL),
  fCutsTOFPions(NULL),
  fCutsTOFKaons(NULL),
  fCutsTOFProtons(NULL),
  fCutsTPCElectrons(NULL),
  fCutsTPCPions(NULL),
  fCutsTPCKaons(NULL),
  fCutsTPCProtons(NULL),
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
  fTPCdedxAfterTOFpidPions=new TH2F("fTPCsignalAfterTOFpions",";p [GeV/c];dEdx",kPBins,binsPDummy,500,0,500);
  fOutputList->Add(fTPCdedxAfterTOFpidPions);
  fTPCdedxAfterTOFpidKaons=new TH2F("fTPCsignalAfterTOFkaons",";p [GeV/c];dEdx",kPBins,binsPDummy,500,0,500);
  fOutputList->Add(fTPCdedxAfterTOFpidKaons);
  fTPCdedxAfterTOFpidProtons=new TH2F("fTPCsignalAfterTOFprotons",";p [GeV/c];dEdx",kPBins,binsPDummy,500,0,500);
  fOutputList->Add(fTPCdedxAfterTOFpidProtons);
  fTPCsignalAfterPionCuts=new TH2F("fTPCsignalAfterPionCuts",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fTPCsignalAfterKaonCuts=new TH2F("fTPCsignalAfterKaonCuts",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fTPCsignalAfterProtonCuts=new TH2F("fTPCsignalAfterProtonCuts",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fOutputList->Add(fTPCsignalAfterPionCuts);
  fOutputList->Add(fTPCsignalAfterKaonCuts);
  fOutputList->Add(fTPCsignalAfterProtonCuts);
  fTPCsignalAfterPionCuts1=new TH2F("fTPCsignalAfterPionCuts1",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fTPCsignalAfterKaonCuts1=new TH2F("fTPCsignalAfterKaonCuts1",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fTPCsignalAfterProtonCuts1=new TH2F("fTPCsignalAfterProtonCuts1",";p [GeV/c];dE/dx",kPBins,binsPDummy,500, 0., 500.);//
  fOutputList->Add(fTPCsignalAfterPionCuts1);
  fOutputList->Add(fTPCsignalAfterKaonCuts1);
  fOutputList->Add(fTPCsignalAfterProtonCuts1);

  fTPCsignalPi=new TH2F("fTPCsignalPi",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for pi+
  fTPCsignalK=new TH2F("fTPCsignalK",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for K+
  fTPCsignalP=new TH2F("fTPCsignalP",";p [GeV/c];signal",kPBins,binsPDummy,300,-2,2);//TPC PID signal as function of p for p
  fOutputList->Add(fTPCsignalPi);
  fOutputList->Add(fTPCsignalK);
  fOutputList->Add(fTPCsignalP);
  fTPCsignalPiafter=new TH2F("fTPCsignalPiafter",";p [GeV/c];(dE/dx-dE/dx_{#pi})/dE/dx_{#pi}",kPBins,binsPDummy,300, -2., 2.);//
  fTPCsignalKafter=new TH2F("fTPCsignalKafter",";p [GeV/c];(dE/dx-dE/dx_{K})/dE/dx_{K}",kPBins,binsPDummy,300, -2., 2.);//
  fTPCsignalPafter=new TH2F("fTPCsignalPafter",";p [GeV/c];(dE/dx-dE/dx_{p})/dE/dx_{p}",kPBins,binsPDummy,300, -2., 2.);//
  fOutputList->Add(fTPCsignalPiafter);
  fOutputList->Add(fTPCsignalKafter);
  fOutputList->Add(fTPCsignalPafter);

  if(fMC)
  {
    fTPCsignalPimc=new TH2F("fTPCsignalPionsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-2,2);//TPC PID signal as function of pt for pi+
    fTPCsignalKmc=new TH2F("fTPCsignalKaonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-2,2);//TPC PID signal as function of pt for K+
    fTPCsignalPmc=new TH2F("fTPCsignalProtonsMC",";p [GeV/c];signal",kPBins,binsPDummy,600,-2,2);//TPC PID signal as function of pt for p
    fOutputList->Add(fTPCsignalPimc);
    fOutputList->Add(fTPCsignalKmc);
    fOutputList->Add(fTPCsignalPmc);
  }

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

  fTOFbetaAfterElectronsCuts=new TH2F("fTOFbetaAfterElectronsCuts",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterPionCuts=new TH2F("fTOFbetaAfterPionCuts",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterKaonCuts=new TH2F("fTOFbetaAfterKaonCuts",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterProtonCuts=new TH2F("fTOFbetaAfterProtonCuts",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fOutputList->Add(fTOFbetaAfterElectronsCuts);
  fOutputList->Add(fTOFbetaAfterPionCuts);
  fOutputList->Add(fTOFbetaAfterKaonCuts);
  fOutputList->Add(fTOFbetaAfterProtonCuts);
  fTOFbetaAfterElectronsCuts1=new TH2F("fTOFbetaAfterElectronsCuts1",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterPionCuts1=new TH2F("fTOFbetaAfterPionCuts1",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterKaonCuts1=new TH2F("fTOFbetaAfterKaonCuts1",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fTOFbetaAfterProtonCuts1=new TH2F("fTOFbetaAfterProtonCuts1",";p [GeV/c];#beta",kPBins,binsPDummy,1000, 0.4, 1.1);//
  fOutputList->Add(fTOFbetaAfterElectronsCuts1);
  fOutputList->Add(fTOFbetaAfterPionCuts1);
  fOutputList->Add(fTOFbetaAfterKaonCuts1);
  fOutputList->Add(fTOFbetaAfterProtonCuts1);

  fTOFbetaEafter=new TH2F("fTOFbetaEafter",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25 );//
  fTOFbetaPiafter=new TH2F("fTOFbetaPiafter",";p [GeV/c];#beta-#beta_{#pi}",kPBins,binsPDummy,500, -0.25, 0.25 );//
  fTOFbetaKafter=new TH2F("fTOFbetaKafter",";p [GeV/c];#beta-#beta_{K}",kPBins,binsPDummy,500, -0.25, 0.25 );//
  fTOFbetaPafter=new TH2F("fTOFbetaPafter",";p [GeV/c];#beta-#beta_{p}",kPBins,binsPDummy,500, -0.25, 0.25 );//
  fOutputList->Add(fTOFbetaEafter);
  fOutputList->Add(fTOFbetaPiafter);
  fOutputList->Add(fTOFbetaKafter);
  fOutputList->Add(fTOFbetaPafter);

  if (fMC)
  {
    fTOFyieldSelEmcE = new TH1F("fTOFyieldSelEmcE",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPimcE = new TH1F("fTOFyieldSelPimcE",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelKmcE = new TH1F("fTOFyieldSelKmcE",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPmcE = new TH1F("fTOFyieldSelPmcE",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTOFyieldSelEmcE);
    fOutputList->Add(fTOFyieldSelPimcE);
    fOutputList->Add(fTOFyieldSelKmcE);
    fOutputList->Add(fTOFyieldSelPmcE);
    fTOFyieldSelEmcM = new TH1F("fTOFyieldSelEmcM",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPimcM = new TH1F("fTOFyieldSelPimcM",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelKmcM = new TH1F("fTOFyieldSelKmcM",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPmcM = new TH1F("fTOFyieldSelPmcM",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTOFyieldSelEmcM);
    fOutputList->Add(fTOFyieldSelPimcM);
    fOutputList->Add(fTOFyieldSelKmcM);
    fOutputList->Add(fTOFyieldSelPmcM);
    fTOFyieldSelEmcPi = new TH1F("fTOFyieldSelEmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPimcPi = new TH1F("fTOFyieldSelPimcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelKmcPi = new TH1F("fTOFyieldSelKmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPmcPi = new TH1F("fTOFyieldSelPmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTOFyieldSelEmcPi);
    fOutputList->Add(fTOFyieldSelPimcPi);
    fOutputList->Add(fTOFyieldSelKmcPi);
    fOutputList->Add(fTOFyieldSelPmcPi);
    fTOFyieldSelEmcK = new TH1F("fTOFyieldSelEmcK",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPimcK = new TH1F("fTOFyieldSelPimcK",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelKmcK = new TH1F("fTOFyieldSelKmcK",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPmcK = new TH1F("fTOFyieldSelPmcK",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTOFyieldSelEmcK);
    fOutputList->Add(fTOFyieldSelPimcK);
    fOutputList->Add(fTOFyieldSelKmcK);
    fOutputList->Add(fTOFyieldSelPmcK);
    fTOFyieldSelEmcP = new TH1F("fTOFyieldSelEmcP",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPimcP = new TH1F("fTOFyieldSelPimcP",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelKmcP = new TH1F("fTOFyieldSelKmcP",";p [Gev/c];",kPBins,binsPDummy);
    fTOFyieldSelPmcP = new TH1F("fTOFyieldSelPmcP",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTOFyieldSelEmcP);
    fOutputList->Add(fTOFyieldSelPimcP);
    fOutputList->Add(fTOFyieldSelKmcP);
    fOutputList->Add(fTOFyieldSelPmcP);
    fTPCyieldSelEmcE = new TH1F("fTPCyieldSelEmcE",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPimcE = new TH1F("fTPCyieldSelPimcE",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelKmcE = new TH1F("fTPCyieldSelKmcE",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPmcE = new TH1F("fTPCyieldSelPmcE",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTPCyieldSelEmcE);
    fOutputList->Add(fTPCyieldSelPimcE);
    fOutputList->Add(fTPCyieldSelKmcE);
    fOutputList->Add(fTPCyieldSelPmcE);
    fTPCyieldSelEmcM = new TH1F("fTPCyieldSelEmcM",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPimcM = new TH1F("fTPCyieldSelPimcM",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelKmcM = new TH1F("fTPCyieldSelKmcM",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPmcM = new TH1F("fTPCyieldSelPmcM",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTPCyieldSelEmcM);
    fOutputList->Add(fTPCyieldSelPimcM);
    fOutputList->Add(fTPCyieldSelKmcM);
    fOutputList->Add(fTPCyieldSelPmcM);
    fTPCyieldSelEmcPi = new TH1F("fTPCyieldSelEmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPimcPi = new TH1F("fTPCyieldSelPimcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelKmcPi = new TH1F("fTPCyieldSelKmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPmcPi = new TH1F("fTPCyieldSelPmcPi",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTPCyieldSelEmcPi);
    fOutputList->Add(fTPCyieldSelPimcPi);
    fOutputList->Add(fTPCyieldSelKmcPi);
    fOutputList->Add(fTPCyieldSelPmcPi);
    fTPCyieldSelEmcK = new TH1F("fTPCyieldSelEmcK",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPimcK = new TH1F("fTPCyieldSelPimcK",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelKmcK = new TH1F("fTPCyieldSelKmcK",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPmcK = new TH1F("fTPCyieldSelPmcK",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTPCyieldSelEmcK);
    fOutputList->Add(fTPCyieldSelPimcK);
    fOutputList->Add(fTPCyieldSelKmcK);
    fOutputList->Add(fTPCyieldSelPmcK);
    fTPCyieldSelEmcP = new TH1F("fTPCyieldSelEmcP",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPimcP = new TH1F("fTPCyieldSelPimcP",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelKmcP = new TH1F("fTPCyieldSelKmcP",";p [Gev/c];",kPBins,binsPDummy);
    fTPCyieldSelPmcP = new TH1F("fTPCyieldSelPmcP",";p [Gev/c];",kPBins,binsPDummy);
    fOutputList->Add(fTPCyieldSelEmcP);
    fOutputList->Add(fTPCyieldSelPimcP);
    fOutputList->Add(fTPCyieldSelKmcP);
    fOutputList->Add(fTPCyieldSelPmcP);
  }

  fPvsPt=new TH2F("fPvsPt","p vs p_{t};p [GeV/c];p_{t} [GeV/c]",kPBins,binsPDummy,kPtBins,binsPtDummy);
  fOutputList->Add(fPvsPt);

  fMeanPvsP = new TProfile("fMeanPvsP","Mean P vs P;p [Gev/c];<p> [GeV/c]",kPBins,binsPDummy);
  fOutputList->Add(fMeanPvsP);

  fTPCvsGlobalMult = new TH2F("fTPCvsGlobalMult","TPC only vs Global track multiplicity;global;TPC only",500,0,2500,500,0,3500);
  fOutputList->Add(fTPCvsGlobalMult);

  fStandardGlobalCuts = AliFlowTrackCuts::GetStandardGlobalTrackCuts2010();
  fStandardTPCCuts = AliFlowTrackCuts::GetStandardTPCStandaloneTrackCuts2010();

  fCutsTOFElectrons = new AliFlowTrackCuts("tof electron cuts");
  fCutsTOFElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFElectrons->SetQA();
  fCutsTOFPions = new AliFlowTrackCuts("tof pion cuts");
  fCutsTOFPions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFPions->SetQA();
  fCutsTOFKaons = new AliFlowTrackCuts("tof kaon cuts");
  fCutsTOFKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFKaons->SetQA();
  fCutsTOFProtons = new AliFlowTrackCuts("tof proton cuts");
  fCutsTOFProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTOFbeta);
  fCutsTOFProtons->SetQA();
  fCutsTPCElectrons = new AliFlowTrackCuts("tpc electron cuts");
  fCutsTPCElectrons->SetPID(AliPID::kElectron, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCElectrons->SetQA();
  fCutsTPCPions = new AliFlowTrackCuts("tpc pion cuts");
  fCutsTPCPions->SetPID(AliPID::kPion, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCPions->SetQA();
  fCutsTPCKaons = new AliFlowTrackCuts("tpc kaon cuts");
  fCutsTPCKaons->SetPID(AliPID::kKaon, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCKaons->SetQA();
  fCutsTPCProtons = new AliFlowTrackCuts("tpc proton cuts");
  fCutsTPCProtons->SetPID(AliPID::kProton, AliFlowTrackCuts::kTPCdedx);
  fCutsTPCProtons->SetQA();

  //fOutputList->Add(fESDpid);

  TH1* h=NULL; 
  h = static_cast<TH1*>(fCutsTOFPions->GetQA()->At(0));
  h->SetName(Form("pion direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFKaons->GetQA()->At(0));
  h->SetName(Form("kaon direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFProtons->GetQA()->At(0));
  h->SetName(Form("proton direct %s",h->GetName()));
  fOutputList->Add(fCutsTOFPions->GetQA()->At(0));
  fOutputList->Add(fCutsTOFKaons->GetQA()->At(0));
  fOutputList->Add(fCutsTOFProtons->GetQA()->At(0));

  h = static_cast<TH1*>(fCutsTOFPions->GetQA()->At(1));
  h->SetName(Form("pion direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFKaons->GetQA()->At(1));
  h->SetName(Form("kaon direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFProtons->GetQA()->At(1));
  h->SetName(Form("proton direct %s",h->GetName()));
  fOutputList->Add(fCutsTOFPions->GetQA()->At(1));
  fOutputList->Add(fCutsTOFKaons->GetQA()->At(1));
  fOutputList->Add(fCutsTOFProtons->GetQA()->At(1));

  h = static_cast<TH1*>(fCutsTOFPions->GetQA()->At(2));
  h->SetName(Form("pion direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFKaons->GetQA()->At(2));
  h->SetName(Form("kaon direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFProtons->GetQA()->At(2));
  h->SetName(Form("proton direct %s",h->GetName()));
  fOutputList->Add(fCutsTOFPions->GetQA()->At(2));
  fOutputList->Add(fCutsTOFKaons->GetQA()->At(2));
  fOutputList->Add(fCutsTOFProtons->GetQA()->At(2));

  h = static_cast<TH1*>(fCutsTOFPions->GetQA()->At(3));
  h->SetName(Form("pion direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFKaons->GetQA()->At(3));
  h->SetName(Form("kaon direct %s",h->GetName()));
  h = static_cast<TH1*>(fCutsTOFProtons->GetQA()->At(3));
  h->SetName(Form("proton direct %s",h->GetName()));
  fOutputList->Add(fCutsTOFPions->GetQA()->At(3));
  fOutputList->Add(fCutsTOFKaons->GetQA()->At(3));
  fOutputList->Add(fCutsTOFProtons->GetQA()->At(3));

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

  AliStack* stack=0x0;
  AliMCEvent* mcEvent=NULL;
  if(fMC)
  {
    mcEvent = (AliMCEvent*) MCEvent();
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

    Double_t p=trackESD->GetP();
    Double_t pt=trackESD->Pt();

    if(!(fCuts->IsSelected(trackESD))) continue;

    Int_t label=-1;
    if(fMC) label=trackESD->GetLabel();

    Int_t pdgcode=0;
    if(stack&&fMC)
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
  if(fMC)
    Printf("MC On\n");

  Printf("AliAnalysisTaskPIDflowQA: end of Terminate");
}


//________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTPC(AliESDtrack* t, Int_t pdgcode)
{
  //do TPC pid
  const AliExternalTrackParam* innerParam = t->GetInnerParam();
  if (!innerParam) return;
  Double_t pinTPCglobal=innerParam->GetP();
  Double_t tpcSignal =t ->GetTPCsignal();
  Float_t p=innerParam->P();
  Float_t pt=innerParam->Pt();
  Float_t sigPion     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kPion);
  Float_t sigKaon     = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kKaon);
  Float_t sigProton   = fESDpid->GetTPCResponse().GetExpectedSignal(pinTPCglobal, AliPID::kProton);
  if(!(sigPion>0.0&&sigKaon>0.0&&sigProton>0.0))
    return;

  fTPCsignal->Fill(pinTPCglobal,tpcSignal);

  fTPCsignalPi->Fill(p,(tpcSignal-sigPion)/sigPion);
  fTPCsignalK->Fill(p,(tpcSignal-sigKaon)/sigKaon);
  fTPCsignalP->Fill(p,(tpcSignal-sigProton)/sigProton);

  if (fCutsTPCPions->IsSelected(t)) 
  {
    fTPCsignalAfterPionCuts->Fill(p,tpcSignal);
    fTPCsignalPiafter->Fill(p,(tpcSignal-sigPion)/sigPion);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTPCyieldSelPimcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTPCyieldSelPimcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTPCyieldSelPimcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTPCyieldSelPimcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTPCyieldSelPimcP->Fill(pt);
      }
    }
  }
  if (fCutsTPCKaons->IsSelected(t)) 
  {
    fTPCsignalAfterKaonCuts->Fill(p,tpcSignal);
    fTPCsignalKafter->Fill(p,(tpcSignal-sigKaon)/sigKaon);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTPCyieldSelKmcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTPCyieldSelKmcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTPCyieldSelKmcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTPCyieldSelKmcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTPCyieldSelKmcP->Fill(pt);
      }
    }
  }
  if (fCutsTPCProtons->IsSelected(t)) 
  {
    fTPCsignalAfterProtonCuts->Fill(p,tpcSignal);
    fTPCsignalPafter->Fill(p,(tpcSignal-sigProton)/sigProton);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTPCyieldSelPmcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTPCyieldSelPmcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTPCyieldSelPmcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTPCyieldSelPmcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTPCyieldSelPmcP->Fill(pt);
      }
    }
  }

  if(fMC)
  {
    if(TMath::Abs(pdgcode)==211)
      fTPCsignalPimc->Fill(p,(tpcSignal-sigPion)/sigPion);
    else if(TMath::Abs(pdgcode)==321)
      fTPCsignalKmc->Fill(p,(tpcSignal-sigKaon)/sigKaon);
    else if (TMath::Abs(pdgcode)==2212)
      fTPCsignalPmc->Fill(p,(tpcSignal-sigProton)/sigProton);
  }
}

//______________________________________________________________________________
void AliAnalysisTaskPIDflowQA::pidTOF(AliESDtrack* track, Int_t pdgcode)
{
  //do TOF pid
  Bool_t goodtrack = (track) &&
                     (track->GetStatus() & AliESDtrack::kTOFpid) &&
                     (track->GetTOFsignal() > 12000) &&
                     (track->GetTOFsignal() < 100000) &&
                     (track->GetIntegratedLength() > 365) &&
                     !(track->GetStatus() & AliESDtrack::kTOFmismatch);

  if (!goodtrack) return;

  const AliExternalTrackParam* innerParam = track->GetInnerParam();
  if (!innerParam) return;
  Double_t pinTPCglobal=innerParam->GetP();
  Double_t tpcSignal =track->GetTPCsignal();
  const Float_t c = 2.99792457999999984e-02;  
  Float_t p = track->GetP();
  Float_t pt = track->Pt();
  Float_t l = track->GetIntegratedLength();  
  Float_t trackT0 = fESDpid->GetTOFResponse().GetStartTime(p);
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

  /////////////simple cuts
  Bool_t isPion   = ( (betadiff[2]<0.015) && (betadiff[2]>-0.015) &&
                      (betadiff[3]>0.025) &&
                      (betadiff[4]>0.03) );

  Bool_t isKaon   = ( (betadiff[3]<0.015) && (betadiff[3]>-0.015) &&
                      (betadiff[2]<-0.03) &&
                      (betadiff[4]>0.03) );

  Bool_t isProton = ( (betadiff[4]<0.015) && (betadiff[4]>-0.015) &&
                      (betadiff[3]<-0.025) &&
                      (betadiff[2]<-0.025) );

  if (isPion)     fTOFbetaAfterPionCuts1->Fill(p,beta);
  if (isKaon)     fTOFbetaAfterKaonCuts1->Fill(p,beta);
  if (isProton)   fTOFbetaAfterProtonCuts1->Fill(p,beta);

  //responses
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

  if (fCutsTOFElectrons->PassesTOFbetaCut(track)) 
  {
    fTOFbetaAfterElectronsCuts->Fill(p,beta);
    fTOFbetaEafter->Fill(p,beta-betaHypothesis[0]);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTOFyieldSelEmcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTOFyieldSelEmcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTOFyieldSelEmcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTOFyieldSelEmcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTOFyieldSelEmcP->Fill(pt);
      }
    }
  }
  if (fCutsTOFPions->PassesTOFbetaCut(track)) 
  {
    fTPCdedxAfterTOFpidPions->Fill(pinTPCglobal,tpcSignal);
    fTOFbetaAfterPionCuts->Fill(p,beta);
    fTOFbetaPiafter->Fill(p,beta-betaHypothesis[2]);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTOFyieldSelPimcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTOFyieldSelPimcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTOFyieldSelPimcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTOFyieldSelPimcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTOFyieldSelPimcP->Fill(pt);
      }
    }
  }
  if (fCutsTOFKaons->PassesTOFbetaCut(track)) 
  {
    fTPCdedxAfterTOFpidKaons->Fill(pinTPCglobal,tpcSignal);
    fTOFbetaAfterKaonCuts->Fill(p,beta);
    fTOFbetaKafter->Fill(p,beta-betaHypothesis[3]);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTOFyieldSelKmcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTOFyieldSelKmcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTOFyieldSelKmcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTOFyieldSelKmcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTOFyieldSelKmcP->Fill(pt);
      }
    }
  }
  if (fCutsTOFProtons->PassesTOFbetaCut(track)) 
  {
    fTPCdedxAfterTOFpidProtons->Fill(pinTPCglobal,tpcSignal);
    fTOFbetaAfterProtonCuts->Fill(p,beta);
    fTOFbetaPafter->Fill(p,beta-betaHypothesis[4]);
    if(fMC)
    {
      if (TMath::Abs(pdgcode)==11)
      {
        fTOFyieldSelPmcE->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==13)
      {
        fTOFyieldSelPmcM->Fill(pt);
      }
      if (TMath::Abs(pdgcode)==211)
      {
        fTOFyieldSelPmcPi->Fill(pt);
      }
      else if(TMath::Abs(pdgcode)==321)
      {
        fTOFyieldSelPmcK->Fill(pt);
      }
      else if (TMath::Abs(pdgcode)==2212)
      {
        fTOFyieldSelPmcP->Fill(pt);
      }
    }
  }

}

//______________________________________________________________________________
Float_t AliAnalysisTaskPIDflowQA::Beta(Float_t m, Float_t p) 
{
  //get theoretical beta
  return TMath::Sqrt(1. / (1. + m * m / (p * p)));
}
 
