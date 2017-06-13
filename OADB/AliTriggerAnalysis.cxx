/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//-------------------------------------------------------------------------
//                      Implementation of   Class AliTriggerAnalysis
// This class provides function to check if events have been triggered based 
// on the data in ESD and AODs. The trigger bits, trigger class inputs and 
// only the data (offline trigger) can be used
// Origin: Jan Fiete Grosse-Oetringhaus, CERN
// Current support and development: Evgeny Kryshen, PNPI
//-------------------------------------------------------------------------

#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TIterator.h"
#include "TParameter.h"
#include "TMap.h"
#include "TRandom.h"
#include "TEllipse.h"
#include "AliTriggerAnalysis.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliMultiplicity.h"
#include "AliESDAD.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliESDFMD.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliDAQ.h"
#include "AliESDTrdTrack.h"
#include "AliVCaloTrigger.h"
#include "AliAODTZERO.h"
#include "AliAODEvent.h"
ClassImp(AliTriggerAnalysis)

AliTriggerAnalysis::AliTriggerAnalysis(TString name) :
AliOADBTriggerAnalysis(name.Data()),
fSPDGFOEfficiency(0),
fMC(kFALSE),
fPileupCutsEnabled(kFALSE),
fDoFMD(kFALSE),
fHistList(new TList()),
fHistStat(0),
fHistFiredBitsSPD(0),
fHistSPDClsVsTklAll(0),
fHistSPDClsVsTklCln(0),
fHistV0C012vsTklAll(0),
fHistV0C012vsTklCln(0),
fHistV0MOnVsOfAll(0),
fHistV0MOnVsOfCln(0),
fHistSPDOnVsOfAll(0),
fHistSPDOnVsOfCln(0),
fHistV0C3vs012All(0),
fHistV0C3vs012Cln(0),
fHistSPDVtxPileupAll(0),
fHistSPDVtxPileupCln(0),
fHistV0MOnAll(0),
fHistV0MOnAcc(0),
fHistV0MOnVHM(0),
fHistV0MOfAll(0),
fHistV0MOfAcc(0),
fHistOFOAll(0),
fHistOFOAcc(0),
fHistTKLAll(0),
fHistTKLAcc(0),
fHistVIRvsBCmod4pup(0),
fHistVIRvsBCmod4acc(0),
fHistVIRCln(0),
fHistBBAflagsAll(0),
fHistBBAflagsAcc(0),
fHistBBCflagsAll(0),
fHistBBCflagsAcc(0),
fHistBGAflagsAll(0),
fHistBGAflagsAcc(0),
fHistBGCflagsAll(0),
fHistBGCflagsAcc(0),
fHistAD(0),
fHistADAAll(0),
fHistADAAcc(0),
fHistADCAll(0),
fHistADCAcc(0),
fHistV0AAll(0),
fHistV0AAcc(0),
fHistV0CAll(0),
fHistV0CAcc(0),
fHistTimeZNA(0),
fHistTimeZNC(0),
fHistZDC(0),
fHistTDCZDC(0),
fHistTimeZNSumVsDif(0),
fHistTimeCorrZDC(0),
fHistFMDA(0),
fHistFMDC(0),
fHistFMDSingle(0),
fHistFMDSum(0),
fHistT0(0),
fHistOFOvsTKLAcc(0),
fHistV0MOnVsOfAcc(0),
fTriggerClasses(new TMap)
{
  // constructor
  fHistList->SetName("histos");
  fHistList->SetOwner();
  fTriggerClasses->SetOwner();
}

//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::SetParameters(AliOADBTriggerAnalysis* oadb){
  fZDCCutRefSumCorr     = oadb->GetZDCCutRefSumCorr();
  fZDCCutRefDeltaCorr   = oadb->GetZDCCutRefDeltaCorr();
  fZDCCutSigmaSumCorr   = oadb->GetZDCCutSigmaSumCorr();
  fZDCCutSigmaDeltaCorr = oadb->GetZDCCutSigmaDeltaCorr();
  fZDCCutZNATimeCorrMax = oadb->GetZDCCutZNATimeCorrMax();
  fZDCCutZNATimeCorrMin = oadb->GetZDCCutZNATimeCorrMin();
  fZDCCutZNCTimeCorrMax = oadb->GetZDCCutZNCTimeCorrMax();
  fZDCCutZNCTimeCorrMin = oadb->GetZDCCutZNCTimeCorrMin();
  fSPDClsVsTklA         = oadb->GetSPDClsVsTklA();
  fSPDClsVsTklB         = oadb->GetSPDClsVsTklB();
  fV0C012vsTklA         = oadb->GetV0C012vsTklA();
  fV0C012vsTklB         = oadb->GetV0C012vsTklB();
  fV0MOnVsOfA           = oadb->GetV0MOnVsOfA();
  fV0MOnVsOfB           = oadb->GetV0MOnVsOfB();
  fSPDOnVsOfA           = oadb->GetSPDOnVsOfA();
  fSPDOnVsOfB           = oadb->GetSPDOnVsOfB();
  fVtxMinContributors   = oadb->GetVtxMinContributors();
  fVtxMinZdist          = oadb->GetVtxMinZdist();
  fVtxNSigmaZdist       = oadb->GetVtxNSigmaZdist();
  fVtxNSigmaDiamXY      = oadb->GetVtxNSigmaDiamXY();
  fVtxNSigmaDiamZ       = oadb->GetVtxNSigmaDiamZ();
  fV0CasymA             = oadb->GetV0CasymA();
  fV0CasymB             = oadb->GetV0CasymB();
  fNBCsPast             = oadb->GetNBCsPast();
  fNBCsFuture           = oadb->GetNBCsFuture();
  fVIRBBAflags          = oadb->GetVIRBBAflags();
  fVIRBBCflags          = oadb->GetVIRBBCflags();
  fVIRBGAflags          = oadb->GetVIRBGAflags();
  fVIRBGCflags          = oadb->GetVIRBGCflags();
  fVHMBBAflags          = oadb->GetVHMBBAflags();
  fVHMBBCflags          = oadb->GetVHMBBCflags();
  fVHMBGAflags          = oadb->GetVHMBGAflags();
  fVHMBGCflags          = oadb->GetVHMBGCflags();
  fV0MOnThreshold       = oadb->GetV0MOnThreshold();
  fV0MOfThreshold       = oadb->GetV0MOfThreshold();
  fSPDGFOThreshold      = oadb->GetSPDGFOThreshhold();
  fSH1OuterThreshold    = oadb->GetSH1OuterThreshold();
  fSH2OuterThreshold    = oadb->GetSH2OuterThreshold();
  fTklThreshold         = oadb->GetTklThreshold();
  fFMDLowCut            = oadb->GetFMDLowThreshold();
  fFMDHitCut            = oadb->GetFMDHitThreshold();
  fTRDptHSE             = oadb->GetTRDptHSE();
  fTRDpidHSE            = oadb->GetTRDpidHSE();
  fTRDptHQU             = oadb->GetTRDptHQU();
  fTRDpidHQU            = oadb->GetTRDpidHQU();
  fTRDptHEE             = oadb->GetTRDptHEE();
  fTRDpidHEE            = oadb->GetTRDpidHEE();
  fTRDminSectorHEE      = oadb->GetTRDminSectorHEE();
  fTRDmaxSectorHEE      = oadb->GetTRDmaxSectorHEE();
  fTRDptHJT             = oadb->GetTRDptHJT();
  fTRDnHJT              = oadb->GetTRDnHJT();
}

//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::~AliTriggerAnalysis(){
  delete fHistList;
  delete fTriggerClasses;
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::EnableHistograms(Bool_t isLowFlux){
  // creates the monitoring histograms 
  // dynamical range of histograms can be adapted for pp and pPb via isLowFlux flag)
  
  // do not add these hists to the directory
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);
  fHistStat            = new TH1F("fHistStat","Accepted events;;",262144,-0.5,262143.5);
  fHistFiredBitsSPD    = new TH1F("fHistFiredBitsSPD","SPD GFO Hardware;chip number;events", 1200, -0.5, 1199.5);
  fHistSPDClsVsTklAll  = new TH2F("fHistSPDClsVsTklAll",                  "All events;n tracklets;n clusters",200,0,isLowFlux?200:6000,500,0,isLowFlux?1000:20000);
  fHistSPDClsVsTklCln  = new TH2F("fHistSPDClsVsTklCln","Events cleaned by other cuts;n tracklets;n clusters",200,0,isLowFlux?200:6000,500,0,isLowFlux?1000:20000);
  fHistV0C012vsTklAll  = new TH2F("fHistV0C012vsTklAll",                  "All events;n tracklets;V0C012 multiplicity",150,0,150,150,0,600);
  fHistV0C012vsTklCln  = new TH2F("fHistV0C012vsTklCln","Events cleaned by other cuts;n tracklets;V0C012 multiplicity",150,0,150,150,0,600);
  fHistV0MOnVsOfAll    = new TH2F("fHistV0MOnVsOfAll",                  "All events;Offline V0M;Online V0M",200,0,isLowFlux?1000:50000,400,0,isLowFlux?8000:40000);
  fHistV0MOnVsOfCln    = new TH2F("fHistV0MOnVsOfCln","Events cleaned by other cuts;Offline V0M;Online V0M",200,0,isLowFlux?1000:50000,400,0,isLowFlux?8000:40000);
  fHistSPDOnVsOfAll    = new TH2F("fHistSPDOnVsOfAll",                  "All events;Offline FOR;Online FOR",300,0,isLowFlux?300:1200 ,300,0,isLowFlux?300:1200);
  fHistSPDOnVsOfCln    = new TH2F("fHistSPDOnVsOfCln","Events cleaned by other cuts;Offline FOR;Online FOR",300,0,isLowFlux?300:1200 ,300,0,isLowFlux?300:1200);
  fHistV0C3vs012All    = new TH2F("fHistV0C3vs012All",                  "All events;V0C012 multiplicity;V0C3 multiplicity",200,0,800,300,0,300);
  fHistV0C3vs012Cln    = new TH2F("fHistV0C3vs012Cln","Events cleaned by other cuts;V0C012 multiplicity;V0C3 multiplicity",200,0,800,300,0,300);
  fHistSPDVtxPileupAll = new TH1F("fHistSPDVtxPileupAll",";SPD Vtx pileup",2,0,2);
  fHistSPDVtxPileupCln = new TH1F("fHistSPDVtxPileupCln",";SPD Vtx pileup",2,0,2);
  fHistVIRvsBCmod4pup  = new TH2F("fHistVIRvsBCmod4pup","VIR vs BC%4 for events identified as SPD or V0 pileup;VIR;BC%4",21,-10.5,10.5,4,-0.5,3.5);
  fHistVIRvsBCmod4acc  = new TH2F("fHistVIRvsBCmod4acc","VIR vs BC%4 for accepted events;VIR;BC%4",21,-10.5,10.5,4,-0.5,3.5);
  fHistVIRCln          = new TH1F("fHistVIRCln","Events cleaned by other cuts",21,-10.5,10.5);
  fHistBBAflagsAll     = new TH1F("fHistBBAflagsAll",";BBA flags;",33,-0.5,32.5);
  fHistBBAflagsAcc     = new TH1F("fHistBBAflagsAcc",";BBA flags;",33,-0.5,32.5);
  fHistBBCflagsAll     = new TH1F("fHistBBCflagsAll",";BBC flags;",33,-0.5,32.5);
  fHistBBCflagsAcc     = new TH1F("fHistBBCflagsAcc",";BBC flags;",33,-0.5,32.5);
  fHistBGAflagsAll     = new TH1F("fHistBGAflagsAll",";BGA flags;",33,-0.5,32.5);
  fHistBGAflagsAcc     = new TH1F("fHistBGAflagsAcc",";BGA flags;",33,-0.5,32.5);
  fHistBGCflagsAll     = new TH1F("fHistBGCflagsAll",";BGC flags;",33,-0.5,32.5);
  fHistBGCflagsAcc     = new TH1F("fHistBGCflagsAcc",";BGC flags;",33,-0.5,32.5);
  fHistV0MOnAll        = new TH1F("fHistV0MOnAll",              "All events;Online V0M;",isLowFlux?8000:40000,0,isLowFlux?8000:40000);
  fHistV0MOnAcc        = new TH1F("fHistV0MOnAcc",         "Accepted events;Online V0M;",isLowFlux?8000:40000,0,isLowFlux?8000:40000);
  fHistV0MOnVHM        = new TH1F("fHistV0MOnVHM","Events with VHM clean up;Online V0M;",isLowFlux?8000:40000,0,isLowFlux?8000:40000);
  fHistV0MOfAll        = new TH1F("fHistV0MOfAll",     "All events;Offline V0M;",isLowFlux?1000:50000,0,isLowFlux?1000:50000);
  fHistV0MOfAcc        = new TH1F("fHistV0MOfAcc","Accepted events;Offline V0M;",isLowFlux?1000:50000,0,isLowFlux?1000:50000);
  fHistOFOAll          = new TH1F("fHistOFOAll"  ,     "All events;Online outer FO chips",800,0,800);
  fHistOFOAcc          = new TH1F("fHistOFOAcc"  ,"Accepted events;Online outer FO chips",800,0,800);
  fHistTKLAll          = new TH1F("fHistTKLAll"  ,     "All events;n tracklets;",isLowFlux?200:6000,0,isLowFlux?200:6000);
  fHistTKLAcc          = new TH1F("fHistTKLAcc"  ,"Accepted events;n tracklets;",isLowFlux?200:6000,0,isLowFlux?200:6000);
  fHistAD              = new TH2F("fHistAD", "ADC+ADA vs ADC+ADA;ADC-ADA time (ns);ADC+ADA time (ns)", 300, -150, 150, 300, -50, 250);
  fHistADAAll          = new TH1F("fHistADAAll",     "All events;ADA mean time (ns);", 2000, -100, 100);
  fHistADAAcc          = new TH1F("fHistADAAcc","Accepted events;ADA mean time (ns);", 2000, -100, 100);
  fHistADCAll          = new TH1F("fHistADCAll",     "All events;ADC mean time (ns);", 2000, -100, 100);
  fHistADCAcc          = new TH1F("fHistADCAcc","Accepted events;ADC mean time (ns);", 2000, -100, 100);
  fHistV0AAll          = new TH1F("fHistV0AAll",     "All events;V0A mean time (ns);", 400, -100, 100);
  fHistV0AAcc          = new TH1F("fHistV0AAcc","Accepted events;V0A mean time (ns);", 400, -100, 100);
  fHistV0CAll          = new TH1F("fHistV0CAll",     "All events;V0C mean time (ns);", 400, -100, 100);
  fHistV0CAcc          = new TH1F("fHistV0CAcc","Accepted events;V0C mean time (ns);", 400, -100, 100);
  fHistZDC             = new TH1F("fHistZDC", "ZDC;trigger bits;events", 8, -1.5, 6.5);
  fHistTDCZDC          = new TH1F("fHistTDCZDC", "ZDC;TDC bits;events", 32, -0.5, 32-0.5);
  fHistTimeZNA         = new TH1F("fHistTimeZNA",";Corrected ZNA time (ns);",200,-10,10);
  fHistTimeZNC         = new TH1F("fHistTimeZNC",";Corrected ZNC time (ns);",200,-10,10);
  fHistTimeZNSumVsDif  = new TH2F("fHistTimeZNSumVsDif",";Corrected ZNC-ZNA time (ns);Corrected ZNC+ZNA time (ns)", 100, -5, 5,100,  -5, 5);
  fHistTimeCorrZDC     = new TH2F("fHistTimeCorrZDC"   ,";Corrected ZNC-ZNA time (ns);Corrected ZNC+ZNA time (ns)", 160,-320,320,160,-320,320);
  fHistFMDA            = new TH1F("fHistFMDA", "FMDA;combinations above threshold;", 102, -1.5, 100.5);
  fHistFMDC            = new TH1F("fHistFMDC", "FMDC;combinations above threshold;", 102, -1.5, 100.5);
  fHistFMDSingle       = new TH1F("fHistFMDSingle", "FMD single;multiplicity value;", 1000, 0, 10);
  fHistFMDSum          = new TH1F("fHistFMDSum", "FMD sum;multiplicity value;counts", 1000, 0, 10);
  fHistT0              = new TH1F("fHistT0", ";T0 time (ns);", 100, -25, 25);
  fHistOFOvsTKLAcc     = new TH2F("fHistOFOvsTKLAcc","Accepted events; n tracklets; Online outer FO chips",200,0,isLowFlux?200:6000,200,0,isLowFlux?200:800);
  fHistV0MOnVsOfAcc    = new TH2F("fHistV0MOnVsOfAcc","Accepted events; Offline V0M; Online V0M",200,0,isLowFlux?1000:50000,400,0,isLowFlux?8000:40000);

  TF1* fFuncSPDClsVsTkl = new TF1("fFuncSPDClsVsTkl","[0]+[1]*x",0,fHistSPDClsVsTklCln->GetXaxis()->GetXmax());
  fFuncSPDClsVsTkl->SetParameters(fSPDClsVsTklA,fSPDClsVsTklB);
  fHistSPDClsVsTklCln->GetListOfFunctions()->Add(fFuncSPDClsVsTkl);

  TF1* fFuncV0C012vsTkl = new TF1("fFuncV0C012vsTkl","[0]+[1]*x",0,6);
  fFuncV0C012vsTkl->SetParameters(fV0C012vsTklA,fV0C012vsTklB);
  fHistV0C012vsTklCln->GetListOfFunctions()->Add(fFuncV0C012vsTkl);
  
  TF1* fFuncV0MOnVsOf = new TF1("fFuncV0MOnVsOf","[0]+[1]*x",0,fHistV0MOnVsOfCln->GetXaxis()->GetXmax());
  fFuncV0MOnVsOf->SetParameters(fV0MOnVsOfA,fV0MOnVsOfB);
  fHistV0MOnVsOfCln->GetListOfFunctions()->Add(fFuncV0MOnVsOf);

  TF1* fFuncSPDOnVsOf = new TF1("fFuncSPDOnVsOf","[0]+[1]*x",0,fHistSPDOnVsOfCln->GetXaxis()->GetXmax());
  fFuncSPDOnVsOf->SetParameters(fSPDOnVsOfA,fSPDOnVsOfB);
  fHistSPDOnVsOfCln->GetListOfFunctions()->Add(fFuncSPDOnVsOf);

  TF1* fFuncV0C3vs012 = new TF1("fFuncV0C3vs012","[0]+[1]*x",0,fHistV0C3vs012Cln->GetXaxis()->GetXmax());
  fFuncV0C3vs012->SetParameters(fV0CasymA,fV0CasymB);
  fHistV0C3vs012Cln->GetListOfFunctions()->Add(fFuncV0C3vs012);
  
  TEllipse* ellipse = new TEllipse(fZDCCutRefDeltaCorr,fZDCCutRefSumCorr,fZDCCutSigmaDeltaCorr,fZDCCutSigmaSumCorr);
  ellipse->SetFillStyle(0);
  ellipse->SetLineColor(kMagenta);
  fHistTimeZNSumVsDif->GetListOfFunctions()->Add(ellipse);
  
  fHistList->Add(fHistStat);
  fHistList->Add(fHistFiredBitsSPD);
  fHistList->Add(fHistSPDClsVsTklAll);
  fHistList->Add(fHistSPDClsVsTklCln);
  fHistList->Add(fHistV0C012vsTklAll);
  fHistList->Add(fHistV0C012vsTklCln);
  fHistList->Add(fHistV0MOnVsOfAll);
  fHistList->Add(fHistV0MOnVsOfCln);
  fHistList->Add(fHistSPDOnVsOfAll);
  fHistList->Add(fHistSPDOnVsOfCln);
  fHistList->Add(fHistV0C3vs012All);
  fHistList->Add(fHistV0C3vs012Cln);
  fHistList->Add(fHistSPDVtxPileupAll);
  fHistList->Add(fHistSPDVtxPileupCln);
  fHistList->Add(fHistVIRvsBCmod4pup);
  fHistList->Add(fHistVIRvsBCmod4acc);
  fHistList->Add(fHistVIRCln);
  fHistList->Add(fHistBBAflagsAll);
  fHistList->Add(fHistBBAflagsAcc);
  fHistList->Add(fHistBBCflagsAll);
  fHistList->Add(fHistBBCflagsAcc);
  fHistList->Add(fHistBGAflagsAll);
  fHistList->Add(fHistBGAflagsAcc);
  fHistList->Add(fHistBGCflagsAll);
  fHistList->Add(fHistBGCflagsAcc);
  fHistList->Add(fHistV0MOnAll);
  fHistList->Add(fHistV0MOnAcc);
  fHistList->Add(fHistV0MOnVHM);
  fHistList->Add(fHistV0MOfAll);
  fHistList->Add(fHistV0MOfAcc);
  fHistList->Add(fHistOFOAll);
  fHistList->Add(fHistOFOAcc);
  fHistList->Add(fHistTKLAll);
  fHistList->Add(fHistTKLAcc);
  fHistList->Add(fHistAD);
  fHistList->Add(fHistADAAll);
  fHistList->Add(fHistADAAcc);
  fHistList->Add(fHistADCAll);
  fHistList->Add(fHistADCAcc);
  fHistList->Add(fHistV0AAll);
  fHistList->Add(fHistV0AAcc);
  fHistList->Add(fHistV0CAll);
  fHistList->Add(fHistV0CAcc);
  fHistList->Add(fHistZDC);
  fHistList->Add(fHistTDCZDC);
  fHistList->Add(fHistTimeZNA);
  fHistList->Add(fHistTimeZNC);
  fHistList->Add(fHistTimeZNSumVsDif);
  fHistList->Add(fHistTimeCorrZDC);
  fHistList->Add(fHistFMDA);
  fHistList->Add(fHistFMDC);
  fHistList->Add(fHistFMDSingle);
  fHistList->Add(fHistFMDSum);
  fHistList->Add(fHistT0);
  fHistList->Add(fHistOFOvsTKLAcc);
  fHistList->Add(fHistV0MOnVsOfAcc);
  
  TH1::AddDirectory(oldStatus);
}

TObject* AliTriggerAnalysis::GetHistogram(const char* histName) { 
  return fHistList->FindObject(histName); 
}


//-------------------------------------------------------------------------------------------------
const char* AliTriggerAnalysis::GetTriggerName(Trigger trigger){
  // returns the name of the requested trigger
  // the returned string will only be valid until the next call to this function [not thread-safe]
  
  static TString str;
  
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  
  switch (triggerNoFlags)  {
    case kAcceptAll :      str = "ACCEPT ALL (bypass!)";      break;
    case kMB1 :            str = "MB1";                       break;
    case kMB2 :            str = "MB2";                       break;
    case kMB3 :            str = "MB3";                       break;
    case kSPDGFO :         str = "SPD GFO";                   break;
    case kSPDGFOBits :     str = "SPD GFO Bits";              break;
    case kSPDGFOL0 :       str = "SPD GFO L0 (first layer)";  break;
    case kSPDGFOL1 :       str = "SPD GFO L1 (second layer)"; break;
    case kSPDClsVsTrkBG :  str = "Cluster vs Tracklets";      break;
    case kV0MOnVsOfPileup: str = "V0M on-cs-of pileup";       break;
    case kSPDOnVsOfPileup: str = "SPD on-cs-of pileup";       break;
    case kV0PFPileup:      str = "V0 PF pileup";              break;
    case kSPDVtxPileup:    str = "SPD vertex pileup";         break;
    case kV0Casym:         str = "V0C012 vs V0C3 asymmetry";  break;
    case kADA :            str = "AD A BB";                   break;
    case kADC :            str = "AD C BB";                   break;
    case kADABG :          str = "AD A BG";                   break;
    case kADCBG :          str = "AD C BG";                   break;
    case kV0A :            str = "V0 A BB";                   break;
    case kV0C :            str = "V0 C BB";                   break;
    case kV0OR :           str = "V0 OR BB";                  break;
    case kV0AND :          str = "V0 AND BB";                 break;
    case kV0ABG :          str = "V0 A BG";                   break;
    case kV0CBG :          str = "V0 C BG";                   break;
    case kVHM :            str = "VHM";                       break;
    case kV0M :            str = "V0M";                       break;
    case kTKL :            str = "TKL";                       break;
    case kSH1 :            str = "SH1";                       break;
    case kSH2 :            str = "SH2";                       break;
    case kZDC :            str = "ZDC";                       break;
    case kZDCA :           str = "ZDC A";                     break;
    case kZDCC :           str = "ZDC C";                     break;
    case kZNA :            str = "ZN A";                      break;
    case kZNC :            str = "ZN C";                      break;
    case kZNABG :          str = "ZN A BG";                   break;
    case kZNCBG :          str = "ZN C BG";                   break;
    case kFMDA :           str = "FMD A";                     break;
    case kFMDC :           str = "FMD C";                     break;
    case kFPANY :          str = "SPD GFO | V0 | ZDC | FMD";  break;
    case kNSD1 :           str = "NSD1";                      break;
    case kMB1Prime:        str = "MB1prime";                  break;
    case kZDCTDCA :        str = "ZDC TDC A";                 break;
    case kZDCTDCC :        str = "ZDC TDC C";                 break;
    case kZDCTime :        str = "ZDC Time Cut";              break;
    case kCentral :        str = "V0 Central";                break;
    case kSemiCentral :    str = "V0 Semi-central";           break;
    case kEmcalL0 :        str = "EMCAL";                     break;
    case kTRDHCO :         str = "TRD cosmics";               break;
    case kTRDHJT :         str = "TRD Jet";                   break;
    case kTRDHSE :         str = "TRD Single Electron";       break;
    case kTRDHQU :         str = "TRD Quarkonia";             break;
    case kTRDHEE :         str = "TRD Dielectron";            break;
    default:               str = "";                          break;
  }
  if (trigger & kOfflineFlag) str += " OFFLINE";
  return str;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerFired(const AliVEvent* event, Trigger trigger){
  // checks if an event has been triggered
  if (trigger & kOfflineFlag) return IsOfflineTriggerFired(event, trigger);
  return IsTriggerBitFired(event, trigger);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsTriggerBitFired(const AliVEvent* event, ULong64_t tclass) const {
  // Checks if corresponding bit in mask is on
  // TODO Why we need this function
  ULong64_t trigmask = event->GetTriggerMask();
  return (trigmask & (1ull << (tclass-1)));
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::EvaluateTrigger(const AliVEvent* event, Trigger trigger){
  // evaluates a given trigger
  // trigger combinations are not supported, for that see IsOfflineTriggerFired

  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  Bool_t offline = trigger & kOfflineFlag;
  
  if (!offline) {
    if ( triggerNoFlags==kT0BG
      || triggerNoFlags==kT0Pileup
      || triggerNoFlags==kSPDClsVsTrkBG
      || triggerNoFlags==kV0MOnVsOfPileup
      || triggerNoFlags==kSPDOnVsOfPileup
      || triggerNoFlags==kV0PFPileup
      || triggerNoFlags==kSPDVtxPileup
      || triggerNoFlags==kV0Casym
      || triggerNoFlags==kV0C012vsTklBG
      || triggerNoFlags==kTKL
      || triggerNoFlags==kZDCA
      || triggerNoFlags==kZDCC
      || triggerNoFlags==kZDCTDCA
      || triggerNoFlags==kZDCTDCC
      || triggerNoFlags==kZDCTime
      || triggerNoFlags==kZNA
      || triggerNoFlags==kZNC
      || triggerNoFlags==kZNABG
      || triggerNoFlags==kZNCBG
      || triggerNoFlags==kFMDA
      || triggerNoFlags==kFMDC
      || triggerNoFlags==kTPCLaserWarmUp
      || triggerNoFlags==kTPCHVdip
      || triggerNoFlags==kIncompleteEvent
      || triggerNoFlags==kEMCAL
      || triggerNoFlags==kEmcalL0
      || triggerNoFlags==kEmcalL1GammaHigh
      || triggerNoFlags==kEmcalL1GammaLow
      || triggerNoFlags==kEmcalL1JetHigh
      || triggerNoFlags==kEmcalL1JetLow
      || triggerNoFlags==kTRDHCO
      || triggerNoFlags==kTRDHJT
      || triggerNoFlags==kTRDHSE
      || triggerNoFlags==kTRDHQU
      || triggerNoFlags==kTRDHEE
      ) AliFatal(Form("Online trigger not available for trigger %d", triggerNoFlags));
  } else {
    if (  triggerNoFlags==kCTPV0A 
        ||triggerNoFlags==kCTPV0C
        ||triggerNoFlags==kVHM
        ||triggerNoFlags==kSH1
        ||triggerNoFlags==kSH2
        ||triggerNoFlags==kCentral
        ||triggerNoFlags==kSemiCentral
      ) AliFatal(Form("Offline trigger not available for trigger %d", triggerNoFlags));
  }
  
  switch (triggerNoFlags) {
    case kCTPV0A:          return event->GetHeader()->IsTriggerInputFired("V0A");
    case kCTPV0C:          return event->GetHeader()->IsTriggerInputFired("V0C");
    case kSPDGFO:          return SPDFiredChips(event, !offline, kFALSE, 0); 
    case kSPDGFOL0:        return SPDFiredChips(event, !offline, kFALSE, 1);
    case kSPDGFOL1:        return SPDFiredChips(event, !offline, kFALSE, 2);
    case kSPDClsVsTrkBG:   return IsSPDClusterVsTrackletBG(event);
    case kV0MOnVsOfPileup: return IsV0MOnVsOfPileup(event);
    case kSPDOnVsOfPileup: return IsSPDOnVsOfPileup(event);
    case kV0PFPileup:      return IsV0PFPileup(event);
    case kSPDVtxPileup:    return IsSPDVtxPileup(event);
    case kV0Casym:         return IsV0Casym(event);
    case kV0C012vsTklBG:   return IsV0C012vsTklBG(event);
    case kADA:             return ADTrigger(event, kASide, !offline) == kADBB; 
    case kADC:             return ADTrigger(event, kCSide, !offline) == kADBB;
    case kADABG:           return ADTrigger(event, kASide, !offline) == kADBG;
    case kADCBG:           return ADTrigger(event, kCSide, !offline) == kADBG;
    case kV0A:             return V0Trigger(event, kASide, !offline) == kV0BB; 
    case kV0C:             return V0Trigger(event, kCSide, !offline) == kV0BB;
    case kV0ABG:           return V0Trigger(event, kASide, !offline) == kV0BG;
    case kV0CBG:           return V0Trigger(event, kCSide, !offline) == kV0BG;
    case kVHM:             return VHMTrigger(event,!offline);
    case kV0M:             return V0MTrigger(event,!offline);
    case kTKL:             return TKLTrigger(event);
    case kSH1:             return SH1Trigger(event);
    case kSH2:             return SH2Trigger(event);
    case kT0:              return T0Trigger(event, !offline) == kT0BB;
    case kT0BG:            return T0Trigger(event, !offline) == kT0DecBG;
    case kT0Pileup:        return T0Trigger(event, !offline) == kT0DecPileup;
    case kZDCA:            return ZDCTrigger(event, kASide);
    case kZDCC:            return ZDCTrigger(event, kCSide);
    case kZDCTDCA:         return ZDCTDCTrigger(event, kASide);
    case kZDCTDCC:         return ZDCTDCTrigger(event, kCSide);
    case kZDCTime:         return ZDCTimeTrigger(event);
    case kZNA:             return ZDCTDCTrigger(event,kASide,kTRUE,kFALSE,kFALSE);
    case kZNC:             return ZDCTDCTrigger(event,kCSide,kTRUE,kFALSE,kFALSE);
    case kZNABG:           return ZDCTimeBGTrigger(event,kASide);
    case kZNCBG:           return ZDCTimeBGTrigger(event,kCSide);
    case kFMDA:            return FMDTrigger(event, kASide);
    case kFMDC:            return FMDTrigger(event, kCSide);
    case kTPCLaserWarmUp:  return IsLaserWarmUpTPCEvent(event);
    case kTPCHVdip:        return IsHVdipTPCEvent(event);
    case kIncompleteEvent: return IsIncompleteEvent(event);
    case kEMCAL:           return EMCALCellsTrigger(event);
    case kEmcalL0:         return EMCALTrigger(event,kEmcalL0);
    case kEmcalL1GammaHigh:return EMCALTrigger(event,kEmcalL1GammaHigh);
    case kEmcalL1GammaLow: return EMCALTrigger(event,kEmcalL1GammaLow);
    case kEmcalL1JetHigh:  return EMCALTrigger(event,kEmcalL1JetHigh);
    case kEmcalL1JetLow:   return EMCALTrigger(event,kEmcalL1JetLow);
    case kTRDHCO:          return TRDTrigger(event,kTRDHCO);
    case kTRDHJT:          return TRDTrigger(event,kTRDHJT);
    case kTRDHSE:          return TRDTrigger(event,kTRDHSE);
    case kTRDHQU:          return TRDTrigger(event,kTRDHQU);
    case kTRDHEE:          return TRDTrigger(event,kTRDHEE);
    case kCentral: {
      if (!event->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!event->GetVZEROData()->TestBit(AliVVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return event->GetVZEROData()->GetTriggerBits() & (1<<AliVVZERO::kCTA2andCTC2);
    }
    case kSemiCentral: {
      if (!event->GetVZEROData()) { AliWarning("V0 centrality trigger bits were not filled!"); return kFALSE; }
      if (!event->GetVZEROData()->TestBit(AliVVZERO::kTriggerChargeBitsFilled)) return kFALSE;
      return event->GetVZEROData()->GetTriggerBits() & (1<<AliVVZERO::kCTA1andCTC1);
    }
    default: AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsOfflineTriggerFired(const AliVEvent* event, Trigger trigger){
  // checks if an event has been triggered "offline"
  UInt_t triggerNoFlags = (UInt_t) trigger % (UInt_t) kStartOfFlags;
  if (trigger & kOneParticle) AliError("AliTriggerAnalysis::kOneParticle functionality is obsolete");
  if (trigger & kOneTrack)    AliError("AliTriggerAnalysis::kOneTrack functionality is obsolete");

  Bool_t decision = kFALSE;
  switch (triggerNoFlags) {
    case kAcceptAll:        return kTRUE; 
    case kMB1:              return SPDGFOTrigger(event, 0) ||  V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kMB2:              return SPDGFOTrigger(event, 0) && (V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB);
    case kMB3:              return SPDGFOTrigger(event, 0) &&  V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kSPDGFO:           return SPDGFOTrigger(event, 0);
    case kSPDGFOBits:       return SPDGFOTrigger(event, 1);
    case kV0MOnVsOfPileup:  return IsV0MOnVsOfPileup(event);
    case kSPDOnVsOfPileup:  return IsSPDOnVsOfPileup(event);
    case kV0PFPileup:       return IsV0PFPileup(event);
    case kSPDVtxPileup:     return IsSPDVtxPileup(event);
    case kV0Casym:          return IsV0Casym(event);
    case kV0C012vsTklBG:    return IsV0C012vsTklBG(event);
    case kADA:              return ADTrigger(event, kASide, kFALSE) == kADBB;
    case kADC:              return ADTrigger(event, kCSide, kFALSE) == kADBB;
    case kADABG:            return ADTrigger(event, kASide, kFALSE) == kADBG;
    case kADCBG:            return ADTrigger(event, kCSide, kFALSE) == kADBG;
    case kV0A:              return V0Trigger(event, kASide, kFALSE) == kV0BB;
    case kV0C:              return V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0OR:             return V0Trigger(event, kASide, kFALSE) == kV0BB || V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0AND:            return V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
    case kV0ABG:            return V0Trigger(event, kASide, kFALSE) == kV0BG;
    case kV0CBG:            return V0Trigger(event, kCSide, kFALSE) == kV0BG;
    case kV0M:              return V0MTrigger(event,kFALSE);
    case kTKL:              return TKLTrigger(event);
    case kZDC:              return ZDCTrigger(event, kASide) || ZDCTrigger(event, kCentralBarrel) || ZDCTrigger(event, kCSide);
    case kZDCA:             return ZDCTrigger(event, kASide);
    case kZDCC:             return ZDCTrigger(event, kCSide);
    case kZNA:              return ZDCTDCTrigger(event,kASide,kTRUE,kFALSE,kFALSE);
    case kZNC:              return ZDCTDCTrigger(event,kCSide,kTRUE,kFALSE,kFALSE);
    case kZNABG:            return ZDCTimeBGTrigger(event,kASide);
    case kZNCBG:            return ZDCTimeBGTrigger(event,kCSide);
    case kFMDA:             return FMDTrigger(event, kASide);
    case kFMDC:             return FMDTrigger(event, kCSide);
    case kEMCAL:            return EMCALCellsTrigger(event);
    case kEmcalL0:          return EMCALTrigger(event,kEmcalL0);
    case kEmcalL1GammaHigh: return EMCALTrigger(event,kEmcalL1GammaHigh);
    case kEmcalL1GammaLow:  return EMCALTrigger(event,kEmcalL1GammaLow);
    case kEmcalL1JetHigh:   return EMCALTrigger(event,kEmcalL1JetHigh);
    case kEmcalL1JetLow:    return EMCALTrigger(event,kEmcalL1JetLow);
    case kTRDHCO:           return TRDTrigger(event,kTRDHCO);
    case kTRDHJT:           return TRDTrigger(event,kTRDHJT);
    case kTRDHSE:           return TRDTrigger(event,kTRDHSE);
    case kTRDHQU:           return TRDTrigger(event,kTRDHQU);
    case kTRDHEE:           return TRDTrigger(event,kTRDHEE);
    case kNSD1:             return SPDFiredChips(event, 0) >= 5 || (V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB);
    case kFPANY:            decision |= SPDGFOTrigger(event, 0); 
                            decision |= V0Trigger(event, kASide, kFALSE) == kV0BB;
                            decision |= V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            decision |= ZDCTrigger(event, kASide);
                            decision |= ZDCTrigger(event, kCentralBarrel);
                            decision |= ZDCTrigger(event, kCSide);
                            decision |= FMDTrigger(event, kASide);
                            decision |= FMDTrigger(event, kCSide);
                            return decision; 
    case kMB1Prime:         decision |= SPDGFOTrigger(event, 0) && V0Trigger(event, kASide, kFALSE) == kV0BB;
                            decision |= SPDGFOTrigger(event, 0) && V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            decision |= V0Trigger(event, kASide, kFALSE) == kV0BB && V0Trigger(event, kCSide, kFALSE) == kV0BB;
                            return decision;
    default:                AliFatal(Form("Trigger type %d not implemented", triggerNoFlags));
  }
  
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::SPDFiredChips(const AliVEvent* event, Int_t origin, Int_t fillHists, Int_t layer){
  // returns the number of fired chips in the SPD
  //
  // origin = 0 --> event->GetMultiplicity()->GetNumberOfFiredChips() (filled from clusters)
  // origin = 1 --> event->GetMultiplicity()->TestFastOrFiredChips() (from hardware bits)
  // layer  = 0 --> both layers
  // layer  = 1 --> inner
  // layer  = 2 --> outer
  
  const AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) { 
    AliError("AliVMultiplicity not available"); 
    return -1; 
  }
  
  if (origin == 0) {
    if (layer == 0) return mult->GetNumberOfFiredChips(0) + mult->GetNumberOfFiredChips(1);
    return mult->GetNumberOfFiredChips(layer-1); 
  }
  
  if (origin == 1) {
    Int_t nChips = 0;
    Int_t firstChip = 0;
    Int_t lastChip  = 1200;
    if(layer == 1) lastChip  = 400;
    if(layer == 2) firstChip = 400;

    for (Int_t i=firstChip; i<lastChip; i++) {
      if (mult->TestFastOrFiredChips(i)) {
        // efficiency simulation (if enabled)
        if (fSPDGFOEfficiency) if (gRandom->Uniform() > fSPDGFOEfficiency->GetBinContent(i+1)) continue;
        nChips++;
        if (fillHists) fHistFiredBitsSPD->Fill(i);
      }
    }
    return nChips;
  }
  
  return -1;
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::ADDecision AliTriggerAnalysis::ADTrigger(const AliVEvent* event, AliceSide side, Bool_t online, Int_t fillHists){
  // Returns the AD trigger decision 
  // argument 'online' is used as a switch between online and offline trigger algorithms
  
  const AliVAD* ad = event->GetADData();
  if (!ad) { 
    // print error only for runs from >=2015
    if (event->GetRunNumber()>=208505) AliError("AliVAD not available");
    return kADInvalid; 
  }
  if (side != kASide && side != kCSide) {
    AliError("Invalid AD side argument");
    return kADInvalid;
  }
  
  AliDebug(2,Form("In ADTrigger: %f %f",ad->GetADATime(),ad->GetADCTime()));

  if (online) {
    UShort_t bits = ad->GetTriggerBits();
    if (side==kASide && (bits & 1<<12)) return kADBB;
    if (side==kCSide && (bits & 1<<13)) return kADBB;
    if (side==kASide && (bits & 1<< 3)) return kADBG;
    if (side==kCSide && (bits & 1<< 5)) return kADBG;
  } else {
    if (fillHists==1) {
      fHistAD->Fill(ad->GetADCTime()-ad->GetADATime(),ad->GetADATime()+ad->GetADCTime());
      if (side == kASide) fHistADAAll->Fill(ad->GetADATime());
      if (side == kCSide) fHistADCAll->Fill(ad->GetADCTime());
    } else if (fillHists==2) {
      if (side == kASide) fHistADAAcc->Fill(ad->GetADATime());
      if (side == kCSide) fHistADCAcc->Fill(ad->GetADCTime());
    }
    if      (side == kASide) return (ADDecision) ad->GetADADecision();
    else if (side == kCSide) return (ADDecision) ad->GetADCDecision();
  }
  
  return kADEmpty;
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::V0Decision AliTriggerAnalysis::V0Trigger(const AliVEvent* event, AliceSide side, Bool_t online, Int_t fillHists){
  // Returns the V0 trigger decision 
  // argument 'online' is used as a switch between online and offline trigger algorithms
  
  const AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) { 
    AliError("AliVVZERO not available");  
    return kV0Invalid; 
  }
  if (!vzero->TestBit(AliVVZERO::kDecisionFilled)) {
    AliError("V0 decisions not filled");
    return kV0Invalid;
  }
  if (!vzero->TestBit(AliVVZERO::kOnlineBitsFilled)) {
    AliError("V0 online trigger bits not filled");
  }
  if (side != kASide && side != kCSide) {
    AliError("Invalid V0 side argument");
    return kV0Invalid;
  }
  
  AliDebug(2,Form("In V0Trigger: %f %f",vzero->GetV0ATime(),vzero->GetV0CTime()));
  
  if (online) {
    // Workaround for high multiplicity in V0C trigger (Pb-Pb 2015):
    // high-mult events drop out from online beam-beam trigger window in MC
    if (fMC && side==kCSide && vzero->GetMTotV0C()>1000) return kV0BB;
    
    Int_t begin = (side == kASide) ? 32 :  0;
    Int_t end   = (side == kASide) ? 64 : 32;
    for (Int_t i=begin; i<end; i++) if (vzero->GetBBFlag(i)) return kV0BB;
    for (Int_t i=begin; i<end; i++) if (vzero->GetBGFlag(i)) return kV0BG;
  } else {
    if (fillHists==1) {
      if (side == kASide) fHistV0AAll->Fill(vzero->GetV0ATime());
      if (side == kCSide) fHistV0CAll->Fill(vzero->GetV0CTime());
    } else if (fillHists==2){
      if (side == kASide) fHistV0AAcc->Fill(vzero->GetV0ATime());
      if (side == kCSide) fHistV0CAcc->Fill(vzero->GetV0CTime());
    }
    if      (side == kASide) return (V0Decision) vzero->GetV0ADecision();
    else if (side == kCSide) return (V0Decision) vzero->GetV0CDecision();
  }
  
  return kV0Empty;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTDCTrigger(const AliVEvent* event, AliceSide side, Bool_t useZN, Bool_t useZP, Int_t fillHists) const{
  // Returns if ZDC triggered, based on TDC information 
  
  Bool_t zdcNA = kFALSE;
  Bool_t zdcNC = kFALSE;
  Bool_t zdcPA = kFALSE;
  Bool_t zdcPC = kFALSE;
  
  if (fMC) { // If it's MC, we use the energy
    Double_t minEnergy = 0;
    zdcNA = event->GetZDCN2Energy()>minEnergy;
    zdcNC = event->GetZDCN1Energy()>minEnergy;
    zdcPA = event->GetZDCP2Energy()>minEnergy;
    zdcPC = event->GetZDCP1Energy()>minEnergy;
  } else if (event->GetDataLayoutType()==AliVEvent::kESD){
    const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(event);
    AliESDZDC* esdZDC = esd->GetESDZDC();
    for (Int_t i=0;i<4;i++){
      zdcNA|= esdZDC->GetZDCTDCData(esdZDC->GetZNATDCChannel(),i)!=0;
      zdcNC|= esdZDC->GetZDCTDCData(esdZDC->GetZNCTDCChannel(),i)!=0;
      zdcPA|= esdZDC->GetZDCTDCData(esdZDC->GetZPATDCChannel(),i)!=0;
      zdcPC|= esdZDC->GetZDCTDCData(esdZDC->GetZPCTDCChannel(),i)!=0;
    }
  } else if (event->GetDataLayoutType()==AliVEvent::kAOD){
    const AliAODEvent* aod = dynamic_cast<const AliAODEvent*>(event);
    AliAODZDC* aodZDC = aod->GetZDCData();
    for (Int_t i=0;i<4;i++){
      // 999 is set if corresponding esdZDC->GetZDCTDCData(ch,i) is 0
      zdcNA|= aodZDC->GetZNATDCm(i)<998;
      zdcNC|= aodZDC->GetZNCTDCm(i)<998;
      zdcPA|= aodZDC->GetZPATDCm(i)<998;
      zdcPC|= aodZDC->GetZPCTDCm(i)<998;
    }
  } else {
    return kFALSE;
  }
  
  if (side == kASide) return ((useZP && zdcPA) || (useZN && zdcNA));
  if (side == kCSide) return ((useZP && zdcPC) || (useZN && zdcNC));
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTimeTrigger(const AliVEvent* event, Int_t fillHists) const {
  // This method implements a selection based on the timing in both sides of zdcN
  // It can be used in order to eliminate parasitic collisions
  // usage of uncorrected timings is deprecated
  // TODO: implement selection on AOD in MC
  if(fMC) {
    if (event->GetDataLayoutType()==AliVEvent::kESD) {
      const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(event);
      AliESDZDC *esdZDC = esd->GetESDZDC();
      UInt_t esdFlag = esdZDC->GetESDQuality();
      Bool_t znaFired = (esdFlag & 0x01) == 0x01;
      Bool_t zncFired = (esdFlag & 0x10) == 0x10;
      return znaFired | zncFired;
    } else {
      return kTRUE;
    }
  }
  
  Float_t zna[4]={0};
  Float_t znc[4]={0};

  if (event->GetDataLayoutType()==AliVEvent::kESD) {
    const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(event);
    AliESDZDC* esdZDC = esd->GetESDZDC();
    Int_t detChZNA  = esdZDC->GetZNATDCChannel();
    Int_t detChZNC  = esdZDC->GetZNCTDCChannel();
    if (esd->GetRunNumber()>=245726 && esd->GetRunNumber()<=245793) detChZNA = 10; // use  timing from the common ZNA PMT
    for (Int_t i=0;i<4;i++) zna[i] = esdZDC->GetZDCTDCCorrected(detChZNA,i);
    for (Int_t i=0;i<4;i++) znc[i] = esdZDC->GetZDCTDCCorrected(detChZNC,i);
  } else if (event->GetDataLayoutType()==AliVEvent::kAOD){
    const AliAODEvent* aod = dynamic_cast<const AliAODEvent*>(event);
    AliAODZDC* aodZDC = aod->GetZDCData();
    for (Int_t i=0;i<4;i++) zna[i]=aodZDC->GetZNATDCm(i);
    for (Int_t i=0;i<4;i++) znc[i]=aodZDC->GetZNCTDCm(i);
  } else {
    return kFALSE;
  }
  
  if(fillHists) {
    for (Int_t i=0;i<4;i++) {
      fHistTimeZNA->Fill(zna[i]);
      fHistTimeZNC->Fill(znc[i]);
      for (Int_t j=0;j<4;j++) {
        fHistTimeCorrZDC->Fill(znc[i]-zna[j],znc[i]+zna[j]);
        fHistTimeZNSumVsDif->Fill(znc[i]-zna[j],znc[i]+zna[j]);
      }
    }
  }
  
  for (Int_t i=0;i<4;i++) {
    for (Int_t j=0;j<4;j++) {
      if (TMath::Power((znc[i]-zna[j]-fZDCCutRefDeltaCorr)/fZDCCutSigmaDeltaCorr,2)+
          TMath::Power((znc[i]+zna[j]-fZDCCutRefSumCorr  )/fZDCCutSigmaSumCorr  ,2)<1.0) return kTRUE;
    }
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTimeBGTrigger(const AliVEvent* event, AliceSide side) const{
  // This method implements a selection based on the timing in zdcN
  // It can be used in order to flag background
  if(fMC) return kFALSE;

  Float_t zna[4]={0};
  Float_t znc[4]={0};
  
  if (event->GetDataLayoutType()==AliVEvent::kESD) {
    const AliESDEvent* esd = dynamic_cast<const AliESDEvent*>(event);
    AliESDZDC* esdZDC = esd->GetESDZDC();
    Int_t detChZNA  = esdZDC->GetZNATDCChannel();
    Int_t detChZNC  = esdZDC->GetZNCTDCChannel();
    for (Int_t i=0;i<4;i++) zna[i] = esdZDC->GetZDCTDCCorrected(detChZNA,i);
    for (Int_t i=0;i<4;i++) znc[i] = esdZDC->GetZDCTDCCorrected(detChZNC,i);
  } else if (event->GetDataLayoutType()==AliVEvent::kAOD){
    const AliAODEvent* aod = dynamic_cast<const AliAODEvent*>(event);
    AliAODZDC* aodZDC = aod->GetZDCData();
    for (Int_t i=0;i<4;i++) zna[i] = aodZDC->GetZNATDCm(i);
    for (Int_t i=0;i<4;i++) znc[i] = aodZDC->GetZNCTDCm(i);
  } else {
    return kFALSE;
  }
  
  Bool_t znabadhit = kFALSE;
  Bool_t zncbadhit = kFALSE;
  
  for(Int_t i = 0; i < 4; ++i) {
    Float_t absZNA = TMath::Abs(zna[i]);
    Float_t absZNC = TMath::Abs(znc[i]);
    if(absZNA<fZDCCutZNATimeCorrMax && absZNA>fZDCCutZNATimeCorrMin) znabadhit = kTRUE;
    if(absZNC<fZDCCutZNCTimeCorrMax && absZNC>fZDCCutZNCTimeCorrMin) zncbadhit = kTRUE;
  }
  
  if (side == kASide) return znabadhit;
  if (side == kCSide) return zncbadhit;

  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::ZDCTrigger(const AliVEvent* event, AliceSide side) const {
  // Returns if ZDC triggered

  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
    AliError("ZDCTrigger method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);
  
  AliESDZDC* zdcData = aEsd->GetESDZDC();
  UInt_t quality = zdcData->GetESDQuality();

  // from Nora's presentation, general first physics meeting 16.10.09
  static UInt_t zpc  = 0x20;
  static UInt_t znc  = 0x10;
  static UInt_t zem1 = 0x08;
  static UInt_t zem2 = 0x04;
  static UInt_t zpa  = 0x02;
  static UInt_t zna  = 0x01;
  
  if (side == kASide         && ((quality & zpa)  || (quality & zna ))) return kTRUE;
  if (side == kCentralBarrel && ((quality & zem1) || (quality & zem2))) return kTRUE;
  if (side == kCSide         && ((quality & zpc)  || (quality & znc ))) return kTRUE;
  
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Int_t AliTriggerAnalysis::FMDHitCombinations(const AliESDEvent* aEsd, AliceSide side, Int_t fillHists){
  // returns number of hit combinations above threshold
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)
  if (!fDoFMD) return -1;
  
  // Workaround for AliESDEvent::GetFMDData is not const!
  const AliESDFMD* fmdData = (const_cast<AliESDEvent*>(aEsd))->GetFMDData();
  if (!fmdData) {
    AliError("AliESDFMD not available");
    return -1;
  }
  
  Int_t detFrom = (side == kASide) ? 1 : 3;
  Int_t detTo   = (side == kASide) ? 2 : 3;
  
  Int_t triggers = 0;
  Float_t totalMult = 0;
  for (UShort_t det=detFrom;det<=detTo;det++) {
    Int_t nRings = (det == 1 ? 1 : 2);
    for (UShort_t ir = 0; ir < nRings; ir++) {
      Char_t   ring = (ir == 0 ? 'I' : 'O');
      UShort_t nsec = (ir == 0 ? 20  : 40);
      UShort_t nstr = (ir == 0 ? 512 : 256);
      for (UShort_t sec =0; sec < nsec;  sec++) {
        for (UShort_t strip = 0; strip < nstr; strip++) {
          Float_t mult = fmdData->Multiplicity(det,ring,sec,strip);
          if (mult == AliESDFMD::kInvalidMult) continue;
          if (fillHists) fHistFMDSingle->Fill(mult);
          if (mult > fFMDLowCut)
            totalMult = totalMult + mult;
          else {
            if (totalMult > fFMDHitCut) triggers++;
            if (fillHists) fHistFMDSum->Fill(totalMult);
            totalMult = 0;
          }
        }
      }
    }
  }
  return triggers;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::FMDTrigger(const AliVEvent* event, AliceSide side){
  // Returns if the FMD triggered
  // Authors: FMD team, Hans Dalsgaard (code merged from FMD/AliFMDOfflineTrigger)
  if (event->GetDataLayoutType()!=AliVEvent::kESD) return 0;
  return FMDHitCombinations((AliESDEvent*) event, side, kFALSE);
}


//-------------------------------------------------------------------------------------------------
AliTriggerAnalysis::T0Decision AliTriggerAnalysis::T0Trigger(const AliVEvent* event, Bool_t online, Int_t fillHists){
  // Returns the T0 TVDC trigger decision
  //  
  // argument 'online' is used as a switch between online and offline trigger algorithms
  // in online mode return 0TVX 
  // in offline mode in addition check pile-up and background :
  // pile-up read from ESD: check if TVDC (0TVX module name) has more 1 hit;
  // background flag read from ESD : check in given time interval OrA and OrC were correct but TVDC not
  // 
  // Based on an algorithm by Alla Maevskaya
  // TODO: read vtx thresholds from OCDB
  
  if (event->GetDataLayoutType()==AliVEvent::kAOD) {
    // AOD analysis
    const AliAODTZERO* tzero = dynamic_cast<const AliAODEvent*>(event)->GetTZEROData();
    if (!tzero) {
      AliError("AliAODTZERO not available");
      return kT0Invalid;
    }
    if (online) {
      UInt_t input0TVX = 2; 
      if (event->GetRunNumber()<229355) input0TVX = 6;
      if (event->GetRunNumber()<224944) input0TVX = 3;
      if (event->GetHeader()->GetL0TriggerInputs() & 1<<input0TVX) return kT0BB;
    } else {
      if (tzero->GetPileupFlag()) return kT0DecPileup;
      if (tzero->GetBackgroundFlag()) return kT0DecBG;
      if (tzero->GetT0zVertex()>-12.3 && tzero->GetT0zVertex() < 10.3) return kT0BB;
    }
  } 
  else if (event->GetDataLayoutType()==AliVEvent::kESD) {
    // ESD analysis
    const AliESDTZERO* tzero = dynamic_cast<const AliESDEvent*>(event)->GetESDTZERO();
    if (!tzero) {
      AliError("AliESDTZERO not available");
      return kT0Invalid;
    }
    if (online) {
      if (event->GetHeader()->IsTriggerInputFired("0TVX")) return kT0BB;
    } else {
      Float_t tvdc0 = tzero->GetTVDC(0);
      if(fillHists) fHistT0->Fill(tvdc0);
      if (tzero->GetPileupFlag()) return kT0DecPileup;
      if (tzero->GetBackgroundFlag()) return kT0DecBG;
      if (tvdc0>-5 && tvdc0<5 && tvdc0!=0) return kT0BB;
    }
    if (fMC) if (tzero->GetT0zVertex()>-12.3 && tzero->GetT0zVertex() < 10.3) return kT0BB;
  }
  
  return kT0Empty;
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::EMCALCellsTrigger(const AliVEvent* event){
  //
  // Returns the EMCAL trigger decision
  // so far only implemented for LHC11a data
  // see http://alisoft.cern.ch/viewvc/trunk/PWGGA/EMCALTasks/AliEmcalPhysicsSelection.cxx?view=markup&root=AliRoot Revision 56136
  //
  
  Bool_t isFired = kTRUE;
  const Int_t runNumber = event->GetRunNumber();
  
  // Get EMCAL cells
  AliVCaloCells *cells = event->GetEMCALCells();
  const Short_t nCells = cells->GetNumberOfCells();
  
  // count cells above threshold per sm
  Int_t nCellCount[10] = {0,0,0,0,0,0,0,0,0,0};
  for(Int_t iCell=0; iCell<nCells; ++iCell) {
    Short_t cellId = cells->GetCellNumber(iCell);
    Double_t cellE = cells->GetCellAmplitude(cellId);
    Int_t sm       = cellId / (24*48);
    if (cellE>0.1)
      ++nCellCount[sm];
  }
  
  // Trigger decision for LHC11a
  Bool_t isLedEvent = kFALSE;
  if ((runNumber>=144871) && (runNumber<=146860)) {
    if (nCellCount[4] > 100)
      isLedEvent = kTRUE;
    else {
      if ((runNumber>=146858) && (runNumber<=146860)) {
        if (nCellCount[3]>=35)
          isLedEvent = kTRUE;
      }
    }
  }
  
  if (isLedEvent) {
    isFired = kFALSE;
  }
  
  return isFired;
}


//----------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::TRDTrigger(const AliVEvent* event, Trigger trigger){
  // evaluate the TRD trigger conditions,
  // so far HCO, HSE, HQU, HJT, HEE
  if(trigger!=kTRDHCO && trigger!=kTRDHJT && trigger!=kTRDHSE && trigger!=kTRDHQU && trigger!=kTRDHEE) {
    AliWarning("Beware you are erroneously trying to use this function (wrong trigger)");
    return kFALSE;
  }
  
  Int_t nTrdTracks = event->GetNumberOfTrdTracks();
  if (nTrdTracks<=0) return kFALSE;
  if      (trigger==kTRDHCO) return kTRUE;
  else if (trigger!=kTRDHJT) {
    for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
      AliVTrdTrack* trdTrack = event->GetTrdTrack(iTrack);
      if (!trdTrack) continue;
      // for the electron triggers we only consider matched tracks
      if(trigger==kTRDHQU) if (TMath::Abs(trdTrack->Pt())>fTRDptHQU && trdTrack->GetPID()>fTRDpidHQU) return kTRUE;
      if(trigger==kTRDHSE) if (TMath::Abs(trdTrack->Pt())>fTRDptHSE && trdTrack->GetPID()>fTRDpidHSE) return kTRUE; 
      if(trigger==kTRDHEE) if (TMath::Abs(trdTrack->Pt())>fTRDptHSE && trdTrack->GetPID()>fTRDpidHSE && trdTrack->GetSector()>=fTRDminSectorHEE && trdTrack->GetSector()<=fTRDmaxSectorHEE) return kTRUE;
    }
  } 
  else if (trigger==kTRDHJT) {
    Int_t nTracks[90] = { 0 }; // stack-wise counted number of tracks above pt threshold
    for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
      AliVTrdTrack *trdTrack = event->GetTrdTrack(iTrack);    
      if (!trdTrack) continue;
      Int_t globalStack = 5*trdTrack->GetSector() + trdTrack->GetStack();
      // stack-wise counting of tracks above pt threshold for jet trigger
      if (TMath::Abs(trdTrack->GetPt()) >= fTRDptHJT) ++nTracks[globalStack];
    }
    // check if HJT condition is fulfilled in any stack
    for (Int_t iStack = 0; iStack < 90; iStack++) if (nTracks[iStack] >= fTRDnHJT) return kTRUE;
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::EMCALTrigger(const AliVEvent* event, Trigger trigger){
  AliVCaloTrigger* emcalTrigger = event->GetCaloTrigger("EMCAL");
  if (!emcalTrigger) return kFALSE;
  Int_t emcalTriggerBits = 0;
  emcalTrigger->GetTriggerBits(emcalTriggerBits);
  if      (trigger==kEmcalL0         ) { return emcalTriggerBits & 1<<0; }
  else if (trigger==kEmcalL1GammaHigh) { return emcalTriggerBits & 1<<1; }
  else if (trigger==kEmcalL1GammaLow ) { return emcalTriggerBits & 1<<2; }
  else if (trigger==kEmcalL1JetHigh  ) { return emcalTriggerBits & 1<<3; }
  else if (trigger==kEmcalL1JetLow   ) { return emcalTriggerBits & 1<<4; }
  else {
    AliWarning("Beware you are erroneously trying to use this function (wrong trigger)");
    return kFALSE;
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsLaserWarmUpTPCEvent(const AliVEvent* event){
  // This function flags noisy TPC events which can happen during laser warm-up.
  Int_t trackCounter = 0;
  for (Int_t i=0; i<event->GetNumberOfTracks(); i++) {
    AliVTrack *track = dynamic_cast<AliVTrack*>(event->GetTrack(i));
    if (!track) continue;
    if (track->GetTPCNcls() < 30) continue;
    if (TMath::Abs(track->Eta()) > 0.005) continue;
    if (track->Pt() < 4) continue;
    if (track->GetKinkIndex(0) > 0) continue;
    UInt_t status = track->GetStatus();
    if ((status&AliESDtrack::kITSrefit)==AliESDtrack::kITSrefit) continue; // explicitly ask for tracks without ITS refit
    if ((status&AliESDtrack::kTPCrefit)!=AliESDtrack::kTPCrefit) continue;
    if (track->GetTPCsignal() > 10) continue;          // explicitly ask for tracks without dE/dx
    trackCounter++;
  }
  if (trackCounter > 15) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsHVdipTPCEvent(const AliVEvent* event) {
  // This function flags events in which the TPC chamber HV is not at its nominal value
  if (fMC) return kFALSE; // there are no dip events in MC
  if (event->GetRunNumber()>197692) return kFALSE; // no dip events in run2
  if (event->GetDataLayoutType()!=AliVEvent::kESD) {
    AliWarning("IsHVdipTPCEvent method implemented for ESDs only");
    return kFALSE;
  }
  const AliESDEvent* aEsd = dynamic_cast<const AliESDEvent*>(event);

  if (!aEsd->IsDetectorOn(AliDAQ::kTPC)) return kTRUE;
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsIncompleteEvent(const AliVEvent* event){
  // Check whether the event is incomplete 
  // (due to DAQ-HLT issues, it could be only part of the event was saved)
  if (fMC) return kFALSE; // there are no incomplete events on MC
  return const_cast<AliVEvent*>(event)->IsIncompleteDAQ(); 
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsSPDClusterVsTrackletBG(const AliVEvent* event, Int_t fillHists){
  // rejects BG based on the cluster vs tracklet correlation
  // returns true if the event is BG
  if (!fPileupCutsEnabled) return kFALSE;
  const AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) { 
    AliError("No multiplicity object"); 
    return kFALSE; 
  }
  Int_t nTkl = mult->GetNumberOfTracklets();
  Int_t nCls = event->GetNumberOfITSClusters(0) + event->GetNumberOfITSClusters(1);
  if      (fillHists==1) fHistSPDClsVsTklAll->Fill(nTkl,nCls);
  else if (fillHists==2) fHistSPDClsVsTklCln->Fill(nTkl,nCls);
  return nCls > fSPDClsVsTklA + nTkl*fSPDClsVsTklB;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0C012vsTklBG(const AliVEvent* event, Int_t fillHists){
  // rejects BG based on V0C012 vs tracklet correlation
  // returns true if the event is BG
  const AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) { 
    AliError("AliVMultiplicity not available"); 
    return kFALSE; 
  }
  
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }

  Float_t nTkl       = mult->GetNumberOfTracklets();
  Float_t multV0C012 = vzero->GetMTotV0C()-vzero->GetMRingV0C(3);
  if      (fillHists==1) fHistV0C012vsTklAll->Fill(nTkl,multV0C012);
  else if (fillHists==2) fHistV0C012vsTklCln->Fill(nTkl,multV0C012);
  return nTkl < 6 && multV0C012 > fV0C012vsTklA + nTkl*fV0C012vsTklB;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0MOnVsOfPileup(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  if (!fPileupCutsEnabled) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  if (!vzero->TestBit(AliVVZERO::kTriggerChargeBitsFilled)){
    if (!fillHists) AliWarning("V0 trigger charge info not found");
    return kFALSE;
  }
  // V0A0 excluded from online V0A charge sum => excluding also from offline sum for consistency
  Float_t on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
  Float_t of = vzero->GetMTotV0A()-vzero->GetMRingV0A(0)+vzero->GetMTotV0C();  
  if      (fillHists==1) fHistV0MOnVsOfAll->Fill(of,on);
  else if (fillHists==2) fHistV0MOnVsOfCln->Fill(of,on);
  return (on < fV0MOnVsOfA + fV0MOnVsOfB*of);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsSPDOnVsOfPileup(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  if (!fPileupCutsEnabled) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  TBits ofMap = mult->GetFiredChipMap();
  Int_t on = onMap.CountBits(0);
  Int_t of = ofMap.CountBits(0);
  if      (fillHists==1) fHistSPDOnVsOfAll->Fill(of,on);
  else if (fillHists==2) fHistSPDOnVsOfCln->Fill(of,on);
  return (on < fSPDOnVsOfA + fSPDOnVsOfB*of);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0PFPileup(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  if (!fPileupCutsEnabled) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  if (!vzero->TestBit(AliVVZERO::kPastFutureFlagsFilled)){
    if (!fillHists) AliWarning("V0 past future info not found");
    return kFALSE;
  }

  Bool_t vir[21] = {0};
  UChar_t bcMod4 = event->GetBunchCrossNumber()%4;

  for (Int_t bc=0;bc<=20;bc++) {
    UChar_t nBBA=0;
    UChar_t nBBC=0;
    UChar_t nBGA=0;
    UChar_t nBGC=0;
    if (fVIRBBAflags<33) for (Int_t i=0;i<32;i++) nBBA+=vzero->GetPFBBFlag(i+32,bc);
    if (fVIRBBCflags<33) for (Int_t i=0;i<32;i++) nBBC+=vzero->GetPFBBFlag(i   ,bc);
    if (fVIRBGAflags<33) for (Int_t i=0;i<32;i++) nBGA+=vzero->GetPFBGFlag(i+32,bc);
    if (fVIRBGCflags<33) for (Int_t i=0;i<32;i++) nBGC+=vzero->GetPFBGFlag(i   ,bc);
    vir[bc] |= nBBA>=fVIRBBAflags;
    vir[bc] |= nBBC>=fVIRBBCflags;
    vir[bc] |= nBGA>=fVIRBGAflags;
    vir[bc] |= nBGC>=fVIRBGCflags;
    if (fillHists==1) {
      if (bc==10) continue;
      if (!vir[bc]) continue;
      if (!IsSPDOnVsOfPileup(event) && !IsV0MOnVsOfPileup(event)) continue;
      fHistVIRvsBCmod4pup->Fill(10-bc,bcMod4);
    }
  }
  
  // clock index is counting from future to past
  Int_t bcMin = 10 - fNBCsFuture + bcMod4;
  Int_t bcMax = 10 + fNBCsPast   + bcMod4;
  for (Int_t bc=bcMin;bc<=bcMax;bc++) {
    if (bc==10) continue; // skip current bc
    if (bc < 0) continue;
    if (bc >20) continue;
    if (vir[bc]) return kTRUE;
  }

  if (fillHists) {
    for (Int_t bc=0;bc<=20;bc++) {
      if (!vir[bc]) continue;
      if (fillHists==2) fHistVIRCln->Fill(10-bc);
      if (bc==10) continue;
      if (fillHists==1) fHistVIRvsBCmod4acc->Fill(10-bc,bcMod4);
    }
  }
  return kFALSE;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsSPDVtxPileup(const AliVEvent* event, Int_t fillHists) {
  if (!fPileupCutsEnabled) return kFALSE;
  Bool_t pileup = event->IsPileupFromSPD(fVtxMinContributors,fVtxMinZdist,fVtxNSigmaZdist,fVtxNSigmaDiamXY,fVtxNSigmaDiamZ);
  if      (fillHists==1) fHistSPDVtxPileupAll->Fill(pileup);
  else if (fillHists==2) fHistSPDVtxPileupCln->Fill(pileup);
  return pileup;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::IsV0Casym(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Float_t multV0C012 = vzero->GetMRingV0C(0)+vzero->GetMRingV0C(1)+vzero->GetMRingV0C(2);
  Float_t multV0C3   = vzero->GetMRingV0C(3);
  
  if      (fillHists==1) fHistV0C3vs012All->Fill(multV0C012,multV0C3);
  else if (fillHists==2) fHistV0C3vs012Cln->Fill(multV0C012,multV0C3);
  return (multV0C3 < fV0CasymA + fV0CasymB*multV0C012);
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::VHMTrigger(const AliVEvent* event, Int_t fillHists){
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  Int_t nBBA = 0;
  Int_t nBBC = 0;
  Int_t nBGA = 0;
  Int_t nBGC = 0;
  
  for (Int_t i=0;i<32;i++) {
    nBBA += vzero->GetBBFlag(i+32);
    nBBC += vzero->GetBBFlag(i   );
    nBGA += vzero->GetBGFlag(i+32);
    nBGC += vzero->GetBGFlag(i   );
  }
  if (fillHists==1){
    fHistBBAflagsAll->Fill(nBBA);
    fHistBBCflagsAll->Fill(nBBC);
    fHistBGAflagsAll->Fill(nBGA);
    fHistBGCflagsAll->Fill(nBGC);
  } else if (fillHists==2) {
    fHistBBAflagsAcc->Fill(nBBA);
    fHistBBCflagsAcc->Fill(nBBC);
    fHistBGAflagsAcc->Fill(nBGA);
    fHistBGCflagsAcc->Fill(nBGC);
  }
  
  Bool_t vhm = 1;
  vhm *= nBBA>=fVHMBBAflags;
  vhm *= nBBC>=fVHMBBCflags;
  vhm *= nBGA<=fVHMBGAflags;
  vhm *= nBGC<=fVHMBGCflags;
  if (fillHists==1 && vhm) {
    Float_t on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
    fHistV0MOnVHM->Fill(on);
  }
  
  return vhm;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::V0MTrigger(const AliVEvent* event, Bool_t online, Int_t fillHists){
  AliVVZERO* vzero = event->GetVZEROData();
  if (!vzero) {
    AliError("AliVVZERO not available");
    return kFALSE;
  }
  if (online && !vzero->TestBit(AliVVZERO::kTriggerChargeBitsFilled)){
    if (!fillHists) AliWarning("V0 trigger charge info not found");
    return kFALSE;
  }

  Int_t   on = vzero->GetTriggerChargeA()+vzero->GetTriggerChargeC();
  Float_t of = vzero->GetMTotV0A()+vzero->GetMTotV0C();

  if (fillHists==1) {
    if (online) fHistV0MOnAll->Fill(on); 
    else        fHistV0MOfAll->Fill(of); 
  } else if (fillHists==2) {
    if (online) fHistV0MOnAcc->Fill(on);
    else        fHistV0MOfAcc->Fill(of);
    if (!online) fHistV0MOnVsOfAcc->Fill(of,on);
  }
  
  return online ? on>=fV0MOnThreshold: of>=fV0MOfThreshold;
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::SH1Trigger(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  Int_t on = onMap.CountBits(400);
  if      (fillHists==1) fHistOFOAll->Fill(on);
  else if (fillHists==2) fHistOFOAcc->Fill(on);
  return (on>=fSH1OuterThreshold);
}


//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::SH2Trigger(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  TBits onMap = mult->GetFastOrFiredChips();
  Int_t on = onMap.CountBits(400);
  if (fillHists) AliWarning("Histograms not filled. Use SH1Trigger for OFO histograms");
  return (on>=fSH2OuterThreshold);
}

//-------------------------------------------------------------------------------------------------
Bool_t AliTriggerAnalysis::TKLTrigger(const AliVEvent* event, Int_t fillHists){
  if (fMC) return kFALSE;
  AliVMultiplicity* mult = event->GetMultiplicity();
  if (!mult) {
    AliError("AliVMultiplicity not available");
    return kFALSE;
  }
  Int_t nTkl = mult->GetNumberOfTracklets();
  if      (fillHists==1) fHistTKLAll->Fill(nTkl);
  else if (fillHists==2) {
    fHistTKLAcc->Fill(nTkl);
    fHistOFOvsTKLAcc->Fill(mult->GetFastOrFiredChips().CountBits(400),nTkl);
  }
  return (nTkl>=fTklThreshold);
}


//-------------------------------------------------------------------------------------------------
Long64_t AliTriggerAnalysis::Merge(TCollection* list){
  // Merge a list of objects with this (needed for PROOF).
  // Returns the number of merged objects (including this).
  if (!list) return 0;
  if (list->IsEmpty()) return 0;
  TIterator* iter = list->MakeIterator();
  TObject* obj;
  TList histListCollection; 
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliTriggerAnalysis* entry = dynamic_cast<AliTriggerAnalysis*> (obj);
    if (entry == 0) continue;
    histListCollection.Add(entry->fHistList);
    count++;
  }
  fHistList->Merge(&histListCollection);
  delete iter;
  return count+1;
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillHistograms(const AliVEvent* event,Bool_t onlineDecision, Bool_t offlineDecision){
  Bool_t pileupCutsStatus = fPileupCutsEnabled;
  fPileupCutsEnabled = kTRUE;

  SPDFiredChips(event,1,kTRUE,0);
  Int_t decisionADA        = ADTrigger(event, kASide, kFALSE, 1);
  Int_t decisionADC        = ADTrigger(event, kCSide, kFALSE, 1);
  Int_t decisionV0A        = V0Trigger(event, kASide, kFALSE, 1);
  Int_t decisionV0C        = V0Trigger(event, kCSide, kFALSE, 1);
  Bool_t isSPDClsVsTklBG   = IsSPDClusterVsTrackletBG(event,1);
  Bool_t isV0C012vsTklBG   = IsV0C012vsTklBG(event,1);
  Bool_t isV0MOnVsOfPileup = IsV0MOnVsOfPileup(event,1);
  Bool_t isSPDOnVsOfPileup = IsSPDOnVsOfPileup(event,1);
  Bool_t isV0PFPileup      = IsV0PFPileup(event,1);
  Bool_t isSPDVtxPileup    = IsSPDVtxPileup(event,1);
  Bool_t isV0Casym         = IsV0Casym(event,1);
  Bool_t isVHMTrigger      = VHMTrigger(event,1);
  Bool_t isV0MOnTrigger    = V0MTrigger(event,kTRUE,1);
  Bool_t isV0MOfTrigger    = V0MTrigger(event,kFALSE,1);
  Bool_t isSH1Trigger      = SH1Trigger(event,1);
  Bool_t isTKLTrigger      = TKLTrigger(event,1);
  Bool_t isZDCTimeTrigger  = ZDCTimeTrigger(event,1);
  Bool_t isZNATimeBG       = ZDCTimeBGTrigger(event,AliTriggerAnalysis::kASide);
  Bool_t isZNCTimeBG       = ZDCTimeBGTrigger(event,AliTriggerAnalysis::kCSide);
  Bool_t isV0A             = decisionV0A==kV0BB;
  Bool_t isV0C             = decisionV0C==kV0BB;
  
  fHistStat->AddBinContent(1);
  if (onlineDecision)  fHistStat->AddBinContent(2);
  if (offlineDecision) fHistStat->AddBinContent(3);
  if (onlineDecision & offlineDecision) fHistStat->AddBinContent(4);
  Int_t accept = 0;
  if (isV0A)              accept |= 1 << 3;
  if (isV0C)              accept |= 1 << 4;
  if (!isSPDClsVsTklBG)   accept |= 1 << 5;
  if (!isV0C012vsTklBG)   accept |= 1 << 6;
  if (!isV0MOnVsOfPileup) accept |= 1 << 7;
  if (!isSPDOnVsOfPileup) accept |= 1 << 8;
  if (!isSPDVtxPileup)    accept |= 1 << 9;
  if (!isV0PFPileup)      accept |= 1 <<10;
  if (!isV0Casym)         accept |= 1 <<11;
  if (isVHMTrigger)       accept |= 1 <<12;
  if (isV0MOfTrigger)     accept |= 1 <<13;
  if (isSH1Trigger)       accept |= 1 <<14;
  if (isZDCTimeTrigger)   accept |= 1 <<15;
  if (!isZNATimeBG)       accept |= 1 <<16;
  if (!isZNCTimeBG)       accept |= 1 <<17;
  if (accept) fHistStat->Fill(accept);
  
  Bool_t acceptDefault = isV0A & isV0C;
  acceptDefault &= !isSPDClsVsTklBG;
  acceptDefault &= !isV0C012vsTklBG;
  acceptDefault &= !isV0PFPileup;
  acceptDefault &= !isSPDVtxPileup;
  acceptDefault &= !isV0MOnVsOfPileup;
  acceptDefault &= !isSPDOnVsOfPileup;
  acceptDefault &= !isV0Casym;

  // Fill histograms for cleaned events
  if (isV0A & isV0C                    & !isV0C012vsTklBG & !isV0PFPileup & !isSPDVtxPileup & !isV0MOnVsOfPileup & !isSPDOnVsOfPileup & !isV0Casym) IsSPDClusterVsTrackletBG(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG                    & !isV0PFPileup & !isSPDVtxPileup & !isV0MOnVsOfPileup & !isSPDOnVsOfPileup & !isV0Casym) IsV0C012vsTklBG(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG & !isV0C012vsTklBG                 & !isSPDVtxPileup & !isV0MOnVsOfPileup & !isSPDOnVsOfPileup & !isV0Casym) IsV0PFPileup(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG & !isV0C012vsTklBG & !isV0PFPileup                   & !isV0MOnVsOfPileup & !isSPDOnVsOfPileup & !isV0Casym) IsSPDVtxPileup(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG & !isV0C012vsTklBG & !isV0PFPileup & !isSPDVtxPileup                      & !isSPDOnVsOfPileup & !isV0Casym) IsV0MOnVsOfPileup(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG & !isV0C012vsTklBG & !isV0PFPileup & !isSPDVtxPileup & !isV0MOnVsOfPileup                      & !isV0Casym) IsSPDOnVsOfPileup(event,2);
  if (isV0A & isV0C & !isSPDClsVsTklBG & !isV0C012vsTklBG & !isV0PFPileup & !isSPDVtxPileup & !isV0MOnVsOfPileup & !isSPDOnVsOfPileup             ) IsV0Casym(event,2);
  
  // Fill distributions for events accepted by general cuts
  if (acceptDefault){
    ADTrigger(event, kASide, kFALSE, 2);
    ADTrigger(event, kCSide, kFALSE, 2);
    V0Trigger(event, kASide, kFALSE, 2);
    V0Trigger(event, kCSide, kFALSE, 2);
    VHMTrigger(event, 2);
    V0MTrigger(event,kTRUE , 2);
    V0MTrigger(event,kFALSE, 2);
    SH1Trigger(event,2);
    TKLTrigger(event,2);
  }
  fPileupCutsEnabled = pileupCutsStatus;

//  TODO: Adjust for AOD
//  AliESDZDC* zdcData = event->GetESDZDC();
//  if (zdcData)  {
//    UInt_t quality = zdcData->GetESDQuality();
//    
//    // from Nora's presentation, general first physics meeting 16.10.09
//    static UInt_t zpc  = 0x20;
//    static UInt_t znc  = 0x10;
//    static UInt_t zem1 = 0x08;
//    static UInt_t zem2 = 0x04;
//    static UInt_t zpa  = 0x02;
//    static UInt_t zna  = 0x01;
//    
//    fHistZDC->Fill(1, (quality & zna)  ? 1 : 0);
//    fHistZDC->Fill(2, (quality & zpa)  ? 1 : 0);
//    fHistZDC->Fill(3, (quality & zem2) ? 1 : 0);
//    fHistZDC->Fill(4, (quality & zem1) ? 1 : 0);
//    fHistZDC->Fill(5, (quality & znc)  ? 1 : 0);
//    fHistZDC->Fill(6, (quality & zpc)  ? 1 : 0);
//  }
//  else {
//    fHistZDC->Fill(-1);
//    AliError("AliESDZDC not available");
//  }
  
//  if (fDoFMD) {
//    fHistFMDA->Fill(FMDHitCombinations(event, kASide, kTRUE));
//    fHistFMDC->Fill(FMDHitCombinations(event, kCSide, kTRUE));
//  }
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::SaveHistograms() const {
  // write histograms to current directory
  if (fSPDGFOEfficiency)   fSPDGFOEfficiency->Write();
  fTriggerClasses->Write("fTriggerClasses", TObject::kSingleKey);
  fHistList->Write("histos",TObject::kSingleKey);
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::FillTriggerClasses(const AliVEvent* event){
  // fills trigger classes map
  TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(event->GetFiredTriggerClasses().Data()));
  if (!count) {
    count = new TParameter<Long64_t>(event->GetFiredTriggerClasses(), 0);
    fTriggerClasses->Add(new TObjString(event->GetFiredTriggerClasses().Data()), count);
  }
  count->SetVal(count->GetVal() + 1);
}


//-------------------------------------------------------------------------------------------------
void AliTriggerAnalysis::PrintTriggerClasses() const {
  // print trigger classes
  
  Printf("Trigger Classes:");
  Printf("Event count for trigger combinations:");
  TMap singleTrigger;
  singleTrigger.SetOwner();
  TIterator* iter = fTriggerClasses->MakeIterator();
  TObjString* obj = 0;
  while ((obj = dynamic_cast<TObjString*> (iter->Next()))) {
    TParameter<Long64_t>* param = static_cast<TParameter<Long64_t>*> (fTriggerClasses->GetValue(obj));
    Printf(" %s: %ld triggers", obj->String().Data(), (Long_t)param->GetVal());
    TObjArray* tokens = obj->String().Tokenize(" ");
    for (Int_t i=0; i<tokens->GetEntries(); i++) {
      TParameter<Long64_t>* count = dynamic_cast<TParameter<Long64_t>*> (singleTrigger.GetValue(((TObjString*) tokens->At(i))->String().Data()));
      if (!count) {
        count = new TParameter<Long64_t>(((TObjString*) tokens->At(i))->String().Data(), 0);
        singleTrigger.Add(new TObjString(((TObjString*) tokens->At(i))->String().Data()), count);
      }
      count->SetVal(count->GetVal() + param->GetVal());
    }
    delete tokens;
  }
  delete iter;
  
  Printf("Event count for single trigger:");
  iter = singleTrigger.MakeIterator();
  while ((obj = dynamic_cast<TObjString*> (iter->Next()))) {
    TParameter<Long64_t>* param = static_cast<TParameter<Long64_t>*> (singleTrigger.GetValue(obj));
    Printf("  %s: %ld triggers", obj->String().Data(), (Long_t)param->GetVal());
  }
  delete iter;
  singleTrigger.DeleteAll();
}

void AliTriggerAnalysis::Browse(TBrowser *b){
   // Browse this object.
   // If b=0, there is no Browse call TObject::Browse(0) instead.
   //         This means TObject::Inspect() will be invoked indirectly
  AliOADBTriggerAnalysis::Browse(b);
  if (b) b->Add(fHistList);
  else   TObject::Browse(b);
}
