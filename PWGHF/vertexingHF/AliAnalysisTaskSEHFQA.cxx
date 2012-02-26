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

/* $Id$ */

/////////////////////////////////////////////////////////////
//
// AliAnalysisTaskSE for HF quality assurance
//
// Author: Chiara Bianchin, chiara.bianchin@pd.infn.it
/////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliVertexerTracks.h"
#include "AliPID.h"
#include "AliTPCPIDResponse.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliCounterCollection.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"
#include "AliInputEventHandler.h"

#include "AliAnalysisTaskSEHFQA.h"

ClassImp(AliAnalysisTaskSEHFQA)

//____________________________________________________________________________

AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA():AliAnalysisTaskSE(),
  fNEntries(0x0),
  fOutputPID(0x0),
  fOutputTrack(0x0),
  fOutputCounters(0x0),
  fOutputCheckCentrality(0x0),
  fOutputEvSelection(0x0),
  fDecayChannel(AliAnalysisTaskSEHFQA::kD0toKpi),
  fCuts(0x0),
  fEstimator(AliRDHFCuts::kCentTRK),
  fReadMC(kFALSE),
  fSimpleMode(kFALSE),
  fOnOff()
{
  //default constructor
  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;
}

//____________________________________________________________________________
AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA(const char *name, AliAnalysisTaskSEHFQA::DecChannel ch,AliRDHFCuts* cuts):
  AliAnalysisTaskSE(name),
  fNEntries(0x0),
  fOutputPID(0x0),
  fOutputTrack(0x0),
  fOutputCounters(0x0),
  fOutputCheckCentrality(0x0),
  fOutputEvSelection(0x0),
  fDecayChannel(ch),
  fCuts(0x0),
  fEstimator(AliRDHFCuts::kCentTRK),
  fReadMC(kFALSE),
  fSimpleMode(kFALSE),
  fOnOff()
{
  //constructor

  //SetCutObject(cuts);
  fCuts=cuts;

  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;

  // Output slot #1 writes into a TH1F container (number of events)
  DefineOutput(1,TH1F::Class());  //My private output
  // Output slot #2 writes into a TList container (PID)
  if (fOnOff[1]) DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TList container (Tracks)
  if (fOnOff[0]) DefineOutput(3,TList::Class());  //My private output
  // Output slot #4 writes into a AliRDHFCuts container (cuts)
  switch(fDecayChannel){
  case 0:
    DefineOutput(4,AliRDHFCutsDplustoKpipi::Class());  //My private output
    break;
  case 1:
    DefineOutput(4,AliRDHFCutsD0toKpi::Class());  //My private output
    break;
  case 2:
    DefineOutput(4,AliRDHFCutsDStartoKpipi::Class());  //My private output
    break;
  case 3:
    DefineOutput(4,AliRDHFCutsDstoKKpi::Class());  //My private output
    break;
  case 4:
    DefineOutput(4,AliRDHFCutsD0toKpipipi::Class());  //My private output
    break;
  case 5:
    DefineOutput(4,AliRDHFCutsLctopKpi::Class());  //My private output
    break;
  }
  if (fOnOff[2]) {
    // Output slot #5 writes into a TList container (AliCounterCollection)
    DefineOutput(5,TList::Class());  //My private output
    // Output slot #6 writes into a TList container (TH1F)
    DefineOutput(6,TList::Class());  //My private output
  }

  if(fOnOff[3]) DefineOutput(7,TList::Class());  //My private output

}

//___________________________________________________________________________
AliAnalysisTaskSEHFQA::~AliAnalysisTaskSEHFQA()
{
  //destructor
  delete fNEntries;

  delete fOutputPID;

  delete fOutputTrack;

  delete fOutputCounters;

  delete fOutputCheckCentrality;

  delete fOutputEvSelection;

}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::Init(){

  //initialization
  if(fDebug > 1) printf("AnalysisTaskSEHFQA::Init() \n");

  switch(fDecayChannel){
  case 0:
    {
      AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;
  case 1:
    {
      AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;
  case 2:
    {
      AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;
  case 3:
    {
      AliRDHFCutsDstoKKpi* copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;
  case 4:
    {
      AliRDHFCutsD0toKpipipi* copycut=new AliRDHFCutsD0toKpipipi(*(static_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;
  case 5:
    {
      AliRDHFCutsLctopKpi* copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fCuts)));
      // Post the data
      PostData(4,copycut);
    }
    break;

  default:
    return;
  }



}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::UserCreateOutputObjects()
{

  //create the output container
  if(fDebug > 1) printf("AnalysisTaskSEHFQA::UserCreateOutputObjects() \n");

  //count events

  fNEntries=new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Counts the number of events", 10,-0.5,9.5);
  fNEntries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNEntries->GetXaxis()->SetBinLabel(2,"Pile-up Rej");
  fNEntries->GetXaxis()->SetBinLabel(3,"No VertexingHF");
  fNEntries->GetXaxis()->SetBinLabel(4,"nCandidates(AnCuts)");
  fNEntries->GetXaxis()->SetBinLabel(5,"EventsWithGoodVtx");
  //fNEntries->GetXaxis()->SetBinLabel(6,"N. of 0SMH");
  fNEntries->GetXaxis()->SetBinLabel(6,"N. of CSH1");
  if(fReadMC){
    fNEntries->GetXaxis()->SetBinLabel(7,"MC Cand from c");
    fNEntries->GetXaxis()->SetBinLabel(8,"MC Cand from b");
    fNEntries->GetXaxis()->SetBinLabel(9,"N fake Trks");
    fNEntries->GetXaxis()->SetBinLabel(10,"N true Trks");
  } else{
    fNEntries->GetXaxis()->SetBinLabel(7,"N candidates");
  }

  fNEntries->GetXaxis()->SetNdivisions(1,kFALSE);

  //PID
  if(fOnOff[1]){
    fOutputPID=new TList();
    fOutputPID->SetOwner();
    fOutputPID->SetName(GetOutputSlot(2)->GetContainer()->GetName());

    //TOF pid
    TH1F* hTOFflags=new TH1F("hTOFflags","TOF flags",6,-0.5,5.5);
    hTOFflags->SetMinimum(0.);
    hTOFflags->GetXaxis()->SetBinLabel(1,"All Tracks");
    hTOFflags->GetXaxis()->SetBinLabel(2,"kTPCout");
    hTOFflags->GetXaxis()->SetBinLabel(3,"kTOFout");
    hTOFflags->GetXaxis()->SetBinLabel(4,"kTIME");
    hTOFflags->GetXaxis()->SetBinLabel(5,"kTOFpid");
    hTOFflags->GetXaxis()->SetBinLabel(6,"kTOFmismatch");

    TString hname="hTOFsig";
    TH1F* hTOFsig=new TH1F(hname.Data(),"Distribution of TOF signal;TOF time [ps];Entries", 100, -2.e3,40.e3);

    hname="hTOFtime";
    TH1F* hTOFtime=new TH1F(hname.Data(),"Distribution of TOF time Kaon;TOF time(Kaon) [ps];Entries", 1000, 0.,50000.);

    hname="hTOFtimeKaonHyptime";
    TH2F* hTOFtimeKaonHyptime=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",500,0.,10.,1000,-20000.,20000.);

    hname="hTOFtimeKaonHyptimeAC";
    TH2F* hTOFtimeKaonHyptimeAC=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",500,0.,10.,1000,-20000.,20000.);

    hname="hTOFsigmaKSigPid";
    TH2F* hTOFsigmaKSigPid=new TH2F(hname.Data(),"(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigmaPionSigPid";
    TH2F* hTOFsigmaPionSigPid=new TH2F(hname.Data(),"(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigmaProtonSigPid";
    TH2F* hTOFsigmaProtonSigPid=new TH2F(hname.Data(),"(TOFsignal-timep)/tofSigPid;p[GeV/c];(TOFsignal-time p)/tofSigPid",500,0.,10.,400,-20,20);

    hname="hTOFsigPid3sigPion";
    TH1F* hTOFsigPid3sigPion=new TH1F(hname.Data(),"TOF PID resolution (#pi) [ps]",500,0.,1000.);

    hname="hTOFsigPid3sigKaon";
    TH1F* hTOFsigPid3sigKaon=new TH1F(hname.Data(),"TOF PID resolution (K) [ps]",500,0.,1000.);

    hname="hTOFsigPid3sigProton";
    TH1F* hTOFsigPid3sigProton=new TH1F(hname.Data(),"TOF PID resolution (p) [ps]",500,0.,1000.);


    //TPC pid
    hname="hTPCsig";
    TH1F* hTPCsig=new TH1F(hname.Data(),"Distribution of TPC signal;TPC sig;Entries", 100, 35.,100.);

    hname="hTPCsigvsp";
    TH2F* hTPCsigvsp=new TH2F(hname.Data(),"TPCsig vs p;TPC p[GeV/c];TPCsig",500,0.,10.,1000,35.,100.);
 
    hname="hTPCsigvspAC";
    TH2F* hTPCsigvspAC=new TH2F(hname.Data(),"TPCsig vs p;TPCp[GeV/c];TPCsig",500,0.,10.,1000,35.,100.);

    hname="hTPCsigmaK";
    TH2F* hTPCsigmaK=new TH2F(hname.Data(),"TPC Sigma for K as a function of momentum;p[GeV/c];Sigma Kaon",500,0.,10.,400,-20,20);

    hname="hTPCsigmaPion";
    TH2F* hTPCsigmaPion=new TH2F(hname.Data(),"TPC Sigma for #pi as a function of momentum;p[GeV/c];Sigma #pi",500,0.,10.,400,-20,20);

    hname="hTPCsigmaProton";
    TH2F* hTPCsigmaProton=new TH2F(hname.Data(),"TPC Sigma for proton as a function of momentum;p[GeV/c];Sigma Proton",500,0.,10.,400,-20,20);

    fOutputPID->Add(hTOFflags);
    fOutputPID->Add(hTOFsig);
    fOutputPID->Add(hTPCsig);
    fOutputPID->Add(hTOFtime);
    fOutputPID->Add(hTOFtimeKaonHyptime);
    fOutputPID->Add(hTOFtimeKaonHyptimeAC);
    fOutputPID->Add(hTOFsigmaKSigPid);
    fOutputPID->Add(hTOFsigmaPionSigPid);
    fOutputPID->Add(hTOFsigmaProtonSigPid);
    fOutputPID->Add(hTOFsigPid3sigPion);
    fOutputPID->Add(hTOFsigPid3sigKaon);
    fOutputPID->Add(hTOFsigPid3sigProton);
    fOutputPID->Add(hTPCsigvsp);
    fOutputPID->Add(hTPCsigvspAC);
    fOutputPID->Add(hTPCsigmaK);
    fOutputPID->Add(hTPCsigmaPion);
    fOutputPID->Add(hTPCsigmaProton);
  }

  //quality of the tracks
  if(fOnOff[0]){
    fOutputTrack=new TList();
    fOutputTrack->SetOwner();
    fOutputTrack->SetName(GetOutputSlot(3)->GetContainer()->GetName());

    TString hname="hnClsITS";
    TH1F* hnClsITS=new TH1F(hname.Data(),"Distribution of number of ITS clusters;nITScls;Entries",7,-0.5,6.5);

    hname="hnClsITSselTr";
    TH1F* hnClsITSselTr=new TH1F(hname.Data(),"Distribution of number of ITS clusters selected tracks;nITScls;Entries",7,-0.5,6.5);

    hname="hnClsITS-SA";
    TH1F* hnClsITSSA=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA);nITScls;Entries",7,-0.5,6.5);


    hname="hnLayerITS";
    TH1F* hnLayerITS=new TH1F(hname.Data(),"Number of tracks with point in layer;ITS layer;",7,-1.5,5.5);
    hnLayerITS->GetXaxis()->SetBinLabel(1,"n tracks");

    hname="hnLayerITSsa";
    TH1F* hnLayerITSsa=new TH1F(hname.Data(),"Number of tracks with point in layer;ITS layer;",7,-1.5,5.5);
    hnLayerITSsa->GetXaxis()->SetBinLabel(1,"n tracks");
   
    hname="hnClsSPD";
    TH1F* hnClsSPD=new TH1F(hname.Data(),"Distribution of number of SPD clusters;nSPDcls;Entries",3,-0.5,2.5);

    hname="hptGoodTr";
    TH1F* hptGoodTr=new TH1F(hname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
    hptGoodTr->SetTitleOffset(1.3,"Y");

    hname="hdistrGoodTr";
    TH1F* hdistrGoodTr=new TH1F(hname.Data(),"Distribution of number of 'good' tracks per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
    hdistrGoodTr->SetTitleOffset(1.3,"Y");

    hname="hd0";
    TH1F* hd0=new TH1F(hname.Data(),"Impact parameter (rphi) distribution of 'good' tracks;d_{0rphi}[cm];Entries/10^{3} cm",200,-0.1,0.1);

    hname="hd0z";
    TH1F* hd0z=new TH1F(hname.Data(),"Impact parameter (z) distribution of 'good' tracks;d_{0z}[cm];Entries/10^{3} cm",200,-0.1,0.1);

    fOutputTrack->Add(hnClsITS);
    fOutputTrack->Add(hnClsITSselTr);
    fOutputTrack->Add(hnClsITSSA);
    fOutputTrack->Add(hnLayerITS);
    fOutputTrack->Add(hnLayerITSsa);
    fOutputTrack->Add(hnClsSPD);
    fOutputTrack->Add(hptGoodTr);
    fOutputTrack->Add(hdistrGoodTr);
    fOutputTrack->Add(hd0);
    fOutputTrack->Add(hd0z);

    if(fReadMC){
      hname="hdistrFakeTr";
      TH1F* hdistrFakeTr=new TH1F(hname.Data(),"Distribution of number of fake tracks per event;no.fake-tracks/ev;Entries",4000,-0.5,3999.5);
      hdistrGoodTr->SetTitleOffset(1.3,"Y");

      hname="hd0f";
      TH1F* hd0f=new TH1F(hname.Data(),"Impact parameter distribution of fake tracks;d_{0}[cm];Entries/10^{3} cm",200,-0.1,0.1);

      hname="hptFakeTr";
      TH1F* hptFakeTr=new TH1F(hname.Data(),"Pt distribution of fake tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      hptFakeTr->SetTitleOffset(1.3,"Y");

      fOutputTrack->Add(hptFakeTr);
      fOutputTrack->Add(hdistrFakeTr);
      fOutputTrack->Add(hd0f);
   
    }
  }

  
  if(fOnOff[2] && fCuts->GetUseCentrality()){

    //Centrality (Counters)
    fOutputCounters=new TList();
    fOutputCounters->SetOwner();
    fOutputCounters->SetName(GetOutputSlot(5)->GetContainer()->GetName());

    AliCounterCollection *stdEstimator=new AliCounterCollection("stdEstimator");
    stdEstimator->AddRubric("run",500000);
    stdEstimator->AddRubric("centralityclass","-10_0/0_10/10_20/20_30/30_40/40_50/50_60/60_70/70_80/80_90/90_100/-990_-980");
    stdEstimator->Init();
    AliCounterCollection *secondEstimator=new AliCounterCollection("secondEstimator");
    secondEstimator->AddRubric("run",500000);
    secondEstimator->AddRubric("centralityclass","-10_0/0_10/10_20/20_30/30_40/40_50/50_60/60_70/70_80/80_90/90_100/-990_-980");
    secondEstimator->Init();
 
    fOutputCounters->Add(stdEstimator);
    fOutputCounters->Add(secondEstimator);
 
    //Centrality (Checks)
    fOutputCheckCentrality=new TList();
    fOutputCheckCentrality->SetOwner();
    fOutputCheckCentrality->SetName(GetOutputSlot(6)->GetContainer()->GetName());

    TString hname="hNtrackletsIn";
    TH1F* hNtrackletsIn=new TH1F(hname.Data(),"Number of tracklets in Centrality range;ntracklets;Entries",5000,-0.5,4999.5);

    hname="hMultIn";
    TH1F* hMultIn=new TH1F(hname.Data(),"Multiplicity;multiplicity in Centrality range;Entries",10000,-0.5,9999.5);

    hname="hNtrackletsOut";
    TH1F* hNtrackletsOut=new TH1F(hname.Data(),"Number of tracklets out of Centrality range;ntracklets;Entries",5000,-0.5,4999.5);

    hname="hMultOut";
    TH1F* hMultOut=new TH1F(hname.Data(),"Multiplicity out of Centrality range;multiplicity;Entries",10000,-0.5,9999.5);

    hname="hMultvsPercentile";
    TH2F* hMultvsPercentile=new TH2F(hname.Data(),"Multiplicity vs Percentile;multiplicity;percentile",10000,-0.5,9999.5,12,-10.,110);


    fOutputCheckCentrality->Add(hNtrackletsIn);
    fOutputCheckCentrality->Add(hNtrackletsOut);
    fOutputCheckCentrality->Add(hMultIn);
    fOutputCheckCentrality->Add(hMultOut);
    fOutputCheckCentrality->Add(hMultvsPercentile);

    PostData(6,fOutputCheckCentrality);
  
  } else{
    if(fOnOff[0]){
      TString hname="hNtracklets";
      TH1F* hNtracklets=new TH1F(hname.Data(),"Number of tracklets;ntracklets;Entries",5000,-0.5,4999.5);

      hname="hMult";
      TH1F* hMult=new TH1F(hname.Data(),"Multiplicity;multiplicity;Entries",10000,-0.5,9999.5);
      fOutputTrack->Add(hNtracklets);
      fOutputTrack->Add(hMult);
    }
  }

  //event selection (z vertex for the moment)
  if(fOnOff[3]){
    fOutputEvSelection=new TList();
    fOutputEvSelection->SetOwner();
    fOutputEvSelection->SetName(GetOutputSlot(7)->GetContainer()->GetName());
    AliCounterCollection *evselection=new AliCounterCollection("evselection");
    evselection->AddRubric("run",500000);
    evselection->AddRubric("evnonsel","zvtx");
    evselection->Init();

    TH1F* hzvtx=new TH1F("hzvtx", "Distribution of z_{VTX};z_{VTX} [cm];Entries",100,-20,20);


    TH2F* hTrigCent=new TH2F("hTrigCent","Centrality vs. Trigger types",12,-1.5,10.5,12,-10,110);
    hTrigCent->GetXaxis()->SetBinLabel(1,"All");
    hTrigCent->GetXaxis()->SetBinLabel(2,"kAny");
    hTrigCent->GetXaxis()->SetBinLabel(3,"kMB");
    hTrigCent->GetXaxis()->SetBinLabel(4,"kINT7");
    hTrigCent->GetXaxis()->SetBinLabel(5,"kCINT5");
    hTrigCent->GetXaxis()->SetBinLabel(6,"kCent");
    hTrigCent->GetXaxis()->SetBinLabel(7,"kSemiCent");
    hTrigCent->GetXaxis()->SetBinLabel(8,"kEMC1+7");
    hTrigCent->GetXaxis()->SetBinLabel(9,"kEMCJET+GAMMA");
    hTrigCent->GetXaxis()->SetBinLabel(10,"Muons");
    hTrigCent->GetXaxis()->SetBinLabel(11,"PHOS");
    hTrigCent->GetXaxis()->SetBinLabel(12,"Others");

    TH2F* hTrigMul=new TH2F("hTrigMul","Multiplicity vs. Trigger types",12,-1.5,10.5,100,0.,10000.);
    hTrigMul->GetXaxis()->SetBinLabel(1,"All");
    hTrigMul->GetXaxis()->SetBinLabel(2,"kAny");
    hTrigMul->GetXaxis()->SetBinLabel(3,"kMB");
    hTrigMul->GetXaxis()->SetBinLabel(4,"kINT7");
    hTrigMul->GetXaxis()->SetBinLabel(5,"kCINT5");
    hTrigMul->GetXaxis()->SetBinLabel(6,"kCent");
    hTrigMul->GetXaxis()->SetBinLabel(7,"kSemiCent");
    hTrigMul->GetXaxis()->SetBinLabel(8,"kEMC1+7");
    hTrigMul->GetXaxis()->SetBinLabel(9,"kEMCJET+GAMMA");
    hTrigMul->GetXaxis()->SetBinLabel(10,"Muons");
    hTrigMul->GetXaxis()->SetBinLabel(11,"PHOS");
    hTrigMul->GetXaxis()->SetBinLabel(12,"Others");

    TH2F* hTrigCentSel=new TH2F("hTrigCentSel","Trigger types",12,-1.5,10.5,12,-10,110);
    hTrigCentSel->GetXaxis()->SetBinLabel(1,"All");
    hTrigCentSel->GetXaxis()->SetBinLabel(2,"kAny");
    hTrigCentSel->GetXaxis()->SetBinLabel(3,"kMB");
    hTrigCentSel->GetXaxis()->SetBinLabel(4,"kINT7");
    hTrigCentSel->GetXaxis()->SetBinLabel(5,"kCINT5");
    hTrigCentSel->GetXaxis()->SetBinLabel(6,"kCent");
    hTrigCentSel->GetXaxis()->SetBinLabel(7,"kSemiCent");
    hTrigCentSel->GetXaxis()->SetBinLabel(8,"kEMC1+7");
    hTrigCentSel->GetXaxis()->SetBinLabel(9,"kEMCJET+GAMMA");
    hTrigCentSel->GetXaxis()->SetBinLabel(10,"Muons");
    hTrigCentSel->GetXaxis()->SetBinLabel(11,"PHOS");
    hTrigCentSel->GetXaxis()->SetBinLabel(12,"Others");

    AliCounterCollection *trigCounter=new AliCounterCollection("trigCounter");
    trigCounter->AddRubric("run",500000);
    trigCounter->AddRubric("triggerType","Any/MB/Cent/SemiCent/EMCAL");
    trigCounter->Init();

    fOutputEvSelection->Add(evselection);
    fOutputEvSelection->Add(hzvtx);
    fOutputEvSelection->Add(hTrigCent);
    fOutputEvSelection->Add(hTrigMul);
    fOutputEvSelection->Add(hTrigCentSel);
    fOutputEvSelection->Add(trigCounter);
  }
//  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
//  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
//  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
//  fCuts->GetPidHF()->SetPidResponse(pidResp);
  // Post the data
  PostData(1,fNEntries);
  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) PostData(5,fOutputCounters);
  if(fOnOff[3]) PostData(7,fOutputEvSelection);

  if(!fOnOff[0] && !fOnOff[1] && !fOnOff[2]) AliError("Nothing will be filled!");
}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(fDebug>2) printf("Analysing decay %d\n",fDecayChannel);
  // Post the data already here
  PostData(1,fNEntries);
  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) {
    PostData(5,fOutputCounters);
    if(fCuts->GetUseCentrality()) PostData(6,fOutputCheckCentrality);
  }

  TClonesArray *arrayProng =0;
  Int_t pdg=0;
  Int_t *pdgdaughters=0x0;

  if(!aod && AODEvent() && IsStandardAOD()) { 
    // In case there is an AOD handler writing a standard AOD, use the AOD 
    // event in memory rather than the input (ESD) event.    
    aod = dynamic_cast<AliAODEvent*> (AODEvent());
    // in this case the braches in the deltaAOD (AliAOD.VertexingHF.root)
    // have to taken from the AOD event hold by the AliAODExtension
    AliAODHandler* aodHandler = (AliAODHandler*) 
      ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      
      AliAODExtension *ext = (AliAODExtension*)aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
   
   
      
      switch(fDecayChannel){
      case 0:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	pdg=411;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=211;//pi
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break; 
      case 1:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
	pdg=421;
	if(fReadMC){
	  pdgdaughters =new Int_t[2];
	  pdgdaughters[0]=211;//pi 
	  pdgdaughters[1]=321;//K
	}
	break; 
      case 2:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
	pdg=413;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[1]=211;//pi
	  pdgdaughters[0]=321;//K
	  pdgdaughters[2]=211;//pi (soft?)
	}
	break; 
      case 3:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	pdg=431;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=321;//K
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break; 
      case 4:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
	pdg=421;
	if(fReadMC){
	  pdgdaughters =new Int_t[4];
	  pdgdaughters[0]=321;
	  pdgdaughters[1]=211;
	  pdgdaughters[2]=211;
	  pdgdaughters[3]=211;
	}
	break; 
      case 5:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	pdg=4122;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=2212;//p
	  pdgdaughters[1]=321;//K
	  pdgdaughters[2]=211;//pi
	}
	break; 
      }
    }
  } else if(aod) {
    switch(fDecayChannel){
    case 0:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=411;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=211;//pi
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break; 
    case 1:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
      pdg=421;
      if(fReadMC){
	pdgdaughters =new Int_t[2];
	pdgdaughters[0]=211;//pi 
	pdgdaughters[1]=321;//K
      }
      break; 
    case 2:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
      pdg=413;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[1]=211;//pi
	pdgdaughters[0]=321;//K
	pdgdaughters[2]=211;//pi (soft?)
      }
      break; 
    case 3:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=431;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=321;//K
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break; 
    case 4:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm4Prong");
      pdg=421;
      if(fReadMC){
	pdgdaughters =new Int_t[4];
	pdgdaughters[0]=321;
	pdgdaughters[1]=211;
	pdgdaughters[2]=211;
	pdgdaughters[3]=211;
      }
      break; 
    case 5:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=4122;
      if(fReadMC){
	pdgdaughters =new Int_t[3];
	pdgdaughters[0]=2212;//p
	pdgdaughters[1]=321;//K
	pdgdaughters[2]=211;//pi
      }
      break; 
    }
  }
  Bool_t isSimpleMode=fSimpleMode;
  if(!arrayProng) {
    AliInfo("Branch not found! The output will contain only trak related histograms\n");
    isSimpleMode=kTRUE;
    fNEntries->Fill(2);
  }
  
  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;

  if(!aod) {
    delete [] pdgdaughters;
    return;
  }

  //check if MC
  if(fReadMC) {
    // load MC particles
    mcArray = (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!mcArray) {
      printf("AliAnalysisTaskSEHFQA::UserExec: MC particles branch not found!\n");
      delete [] pdgdaughters;
      return;
    }
    
    // load MC header
    mcHeader = (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      printf("AliAnalysisTaskSEHFQA::UserExec: MC header branch not found!\n");
      delete [] pdgdaughters;
      return;
    }
  }

  
  UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  Double_t centrality=fCuts->GetCentrality(aod);
  Double_t multiplicity=aod->GetHeader()->GetRefMultiplicity();
  Int_t runNumber = aod->GetRunNumber();
  if(fOnOff[3]){
    TH2F* hTrigC=(TH2F*)fOutputEvSelection->FindObject("hTrigCent");
    TH2F* hTrigM=(TH2F*)fOutputEvSelection->FindObject("hTrigMul");
    AliCounterCollection* trigCount=(AliCounterCollection*)fOutputEvSelection->FindObject("trigCounter");

    hTrigC->Fill(-1.,centrality);
    hTrigM->Fill(-1.,multiplicity);
    
    if(evSelMask & AliVEvent::kAny){
      hTrigC->Fill(0.,centrality);
      hTrigM->Fill(0.,multiplicity);
      trigCount->Count(Form("triggerType:Any/Run:%d",runNumber));
    }  
    if(evSelMask & AliVEvent::kMB){
      hTrigC->Fill(1.,centrality);
      hTrigM->Fill(1.,multiplicity);
      trigCount->Count(Form("triggerType:MB/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kINT7){ 
      hTrigC->Fill(2.,centrality);
      hTrigM->Fill(2.,multiplicity);
    }
    if(evSelMask & AliVEvent::kCINT5){ 
      hTrigC->Fill(3.,centrality);
      hTrigM->Fill(3.,multiplicity);
    }
    if(evSelMask & AliVEvent::kCentral){
      hTrigC->Fill(4.,centrality);
      hTrigM->Fill(4.,multiplicity);
      trigCount->Count(Form("triggerType:Cent/Run:%d",runNumber));
    }
    if(evSelMask & AliVEvent::kSemiCentral){ 
      hTrigC->Fill(5.,centrality);
      hTrigM->Fill(5.,multiplicity);
      trigCount->Count(Form("triggerType:SemiCent/Run:%d",runNumber));
    }
    if(evSelMask & (AliVEvent::kEMC1 | AliVEvent::kEMC7)){
      hTrigC->Fill(6.,centrality);
      hTrigM->Fill(6.,multiplicity);
    }
    if(evSelMask & (AliVEvent::kEMCEJE | AliVEvent::kEMCEGA)){
      hTrigC->Fill(7.,centrality);
      hTrigM->Fill(7.,multiplicity);
      trigCount->Count(Form("triggerType:EMCAL/Run:%d",runNumber));
    }
    if(evSelMask & (((AliVEvent::kCMUS5 | AliVEvent::kMUSH7) | (AliVEvent::kMUL7 | AliVEvent::kMUU7)) |  (AliVEvent::kMUS7 | AliVEvent::kMUON))){
      hTrigC->Fill(8.,centrality);
      hTrigM->Fill(8.,multiplicity);
    }
    if(evSelMask & (AliVEvent::kPHI1 | AliVEvent::kPHI7)){ 
      hTrigC->Fill(9.,centrality);
      hTrigM->Fill(9.,multiplicity);
    }
    if(evSelMask & (AliVEvent::kDG5 | AliVEvent::kZED)){
      hTrigC->Fill(10.,centrality);
      hTrigM->Fill(10.,multiplicity);
    }
  }
  

  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) {
    delete [] pdgdaughters;
    return;
  }

  // count event
  fNEntries->Fill(0); 

  //count events with good vertex
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) fNEntries->Fill(4);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  //TString trigclass=aod->GetFiredTriggerClasses();
  //if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNEntries->Fill(5); //tmp




  Bool_t evSelbyCentrality=kTRUE,evSelected=kTRUE,evSelByVertex=kTRUE,evselByPileup=kTRUE,evSelByPS=kTRUE;
  //select event
  if(!fCuts->IsEventSelected(aod)) {
    evSelected=kFALSE;
    if(fCuts->GetWhyRejection()==1) {fNEntries->Fill(1); evselByPileup=kFALSE;}// rejected for pileup
    if(fCuts->GetWhyRejection()==2 || fCuts->GetWhyRejection()==3) evSelbyCentrality=kFALSE; //rejected by centrality
    if(fCuts->GetWhyRejection()==4) evSelByVertex=kFALSE; //rejected by vertex
    if(fCuts->GetWhyRejection()==5) fNEntries->Fill(5);//tmp
    if(fCuts->GetWhyRejection()==6 && fOnOff[3]) ((AliCounterCollection*)fOutputEvSelection->FindObject("evselection"))->Count(Form("evnonsel:zvtx/Run:%d",runNumber));
    if(fCuts->GetWhyRejection()==7) { evSelByPS=kFALSE; }
  }
  if(evSelected){
    TH2F* hTrigS=(TH2F*)fOutputEvSelection->FindObject("hTrigCentSel");
    hTrigS->Fill(-1.,centrality);

    if(evSelMask & AliVEvent::kAny) hTrigS->Fill(0.,centrality);
    if(evSelMask & AliVEvent::kMB) hTrigS->Fill(1.,centrality);
    if(evSelMask & AliVEvent::kINT7) hTrigS->Fill(2.,centrality);
    if(evSelMask & AliVEvent::kCINT5) hTrigS->Fill(3.,centrality);
    if(evSelMask & AliVEvent::kCentral) hTrigS->Fill(4.,centrality);
    if(evSelMask & AliVEvent::kSemiCentral) hTrigS->Fill(5.,centrality);
    if(evSelMask & (AliVEvent::kEMC1 | AliVEvent::kEMC7)) hTrigS->Fill(6.,centrality);
    if(evSelMask & (AliVEvent::kEMCEJE | AliVEvent::kEMCEGA)) hTrigS->Fill(7.,centrality);
    if(evSelMask & (((AliVEvent::kCMUS5 | AliVEvent::kMUSH7) | (AliVEvent::kMUL7 | AliVEvent::kMUU7)) |  (AliVEvent::kMUS7 | AliVEvent::kMUON))) hTrigS->Fill(8.,centrality);
    if(evSelMask & (AliVEvent::kPHI1 | AliVEvent::kPHI7)) hTrigS->Fill(9.,centrality);
    if(evSelMask & (AliVEvent::kDG5 | AliVEvent::kZED)) hTrigS->Fill(10.,centrality);
  }
  
  if(evSelected || (!evSelbyCentrality && evSelByVertex && evselByPileup && evSelByPS)){ //events selected or not selected because of centrality
    if(fOnOff[2] && fCuts->GetUseCentrality()){

      Float_t stdCentf=fCuts->GetCentrality(aod);
      Int_t stdCent = (Int_t)(stdCentf+0.5);
      Float_t secondCentf =fCuts->GetCentrality(aod,fEstimator);
      Int_t secondCent = (Int_t)(secondCentf+0.5);
      Int_t mincent=stdCent-stdCent%10;
      if(stdCentf==-1) {
	mincent=-10; 
	stdCent=-1;
      }
      if(mincent==100)mincent--;
      ((AliCounterCollection*)fOutputCounters->FindObject("stdEstimator"))->Count(Form("centralityclass:%d_%d/Run:%d",mincent,mincent+10,runNumber));

      mincent=secondCent-secondCent%10;
      if(secondCentf==-1) {
	mincent=-10;
	secondCent=-1;
      }
      if(mincent==100)mincent--;
      ((AliCounterCollection*)fOutputCounters->FindObject("secondEstimator"))->Count(Form("centralityclass:%d_%d/Run:%d",mincent,mincent+10,runNumber));

      if(stdCent<fCuts->GetMinCentrality() || stdCent>fCuts->GetMaxCentrality()){
	((TH1F*)fOutputCheckCentrality->FindObject("hNtrackletsOut"))->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	((TH1F*)fOutputCheckCentrality->FindObject("hMultOut"))->Fill(aod->GetHeader()->GetRefMultiplicity());
      }else{
	((TH1F*)fOutputCheckCentrality->FindObject("hNtrackletsIn"))->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	((TH1F*)fOutputCheckCentrality->FindObject("hMultIn"))->Fill(aod->GetHeader()->GetRefMultiplicity());
      }
      ((TH2F*)fOutputCheckCentrality->FindObject("hMultvsPercentile"))->Fill(aod->GetHeader()->GetRefMultiplicity(),stdCentf);

      PostData(6,fOutputCheckCentrality);

    } else{
      if(fOnOff[0]){
	((TH1F*)fOutputTrack->FindObject("hNtracklets"))->Fill(aod->GetTracklets()->GetNumberOfTracklets());
	((TH1F*)fOutputTrack->FindObject("hMult"))->Fill(aod->GetHeader()->GetRefMultiplicity());
      }
    }
  }

  if(fOnOff[3]){
    const AliVVertex *vertex = aod->GetPrimaryVertex();
    Double_t zvtx=vertex->GetZ();
    if(evSelected || (!evSelected && TMath::Abs(zvtx) > 10.))
    ((TH1F*)fOutputEvSelection->FindObject("hzvtx"))->Fill(zvtx);
  }

  if(!evSelected) {
    delete [] pdgdaughters;
    return; //discard all events not selected (vtx and/or centrality)
  }


  AliAODPidHF* pidHF=fCuts->GetPidHF();
  if(!pidHF) {
    delete [] pdgdaughters;
    return;
  }
  AliPIDResponse* respF=pidHF->GetPidResponse();
  AliTPCPIDResponse* tpcres=new AliTPCPIDResponse();
  if(pidHF->GetOldPid()) pidHF->SetBetheBloch(*tpcres);
  Bool_t oldPID=pidHF->GetOldPid();


  Int_t ntracks=0;
  Int_t isGoodTrack=0, isFakeTrack=0;

  if(aod) ntracks=aod->GetNTracks();

  if(fOnOff[0] || fOnOff[1]){
    //loop on tracks in the event
    for (Int_t k=0;k<ntracks;k++){
      AliAODTrack* track=aod->GetTrack(k);
      AliAODPid *pid = track->GetDetPid();


      if(fOnOff[1]){
	if(!pid)  {if (fDebug>1)cout<<"No AliAODPid found"<<endl; continue;}
	Double_t times[AliPID::kSPECIES];
	pid->GetIntegratedTimes(times);
    
	Double_t tofRes[AliPID::kSPECIES];
	pid->GetTOFpidResolution(tofRes);

	//check TOF
	TH1F* htmpfl=((TH1F*)fOutputPID->FindObject("hTOFflags"));
	htmpfl->Fill(0.);
	if (track->GetStatus()&AliESDtrack::kTPCout) htmpfl->Fill(1.);
	if (track->GetStatus()&AliESDtrack::kTOFout) htmpfl->Fill(2.);
	if (track->GetStatus()&AliESDtrack::kTIME) htmpfl->Fill(3.);
	if (track->GetStatus()&AliESDtrack::kTOFpid) htmpfl->Fill(4.);
	if (track->GetStatus()&AliESDtrack::kTOFmismatch) htmpfl->Fill(5.);

	if(pidHF && pidHF->CheckStatus(track,"TOF")){
	  ((TH1F*)fOutputPID->FindObject("hTOFtime"))->Fill(times[AliPID::kProton]);
	  ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptime"))->Fill(track->P(),pid->GetTOFsignal()-times[3]); //3 is kaon
	  ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(pid->GetTOFsignal());
	  if (pid->GetTOFsignal()< 0) ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(-1);

	  
	  // test TOF sigma PID
	  if (tofRes[2] != 0.) {   // protection against 'old' AODs...
	    ((TH2F*)fOutputPID->FindObject("hTOFsigmaKSigPid"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kKaon])/tofRes[3]);
	    ((TH2F*)fOutputPID->FindObject("hTOFsigmaPionSigPid"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kPion])/tofRes[2]);
	    ((TH2F*)fOutputPID->FindObject("hTOFsigmaProtonSigPid"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kProton])/tofRes[4]);
	    for (Int_t iS=2; iS<5; iS++){ //we plot TOF Pid resolution for 3-sigma identified particles
	      if ( (TMath::Abs(times[iS]-pid->GetTOFsignal())/tofRes[iS])<3.){
		switch (iS) {
		case AliPID::kPion:
		  ((TH1F*)fOutputPID->FindObject("hTOFsigPid3sigPion"))->Fill(tofRes[iS]);
		  break;
		case AliPID::kKaon:
		  ((TH1F*)fOutputPID->FindObject("hTOFsigPid3sigKaon"))->Fill(tofRes[iS]);
		  break;
		case AliPID::kProton:
		  ((TH1F*)fOutputPID->FindObject("hTOFsigPid3sigProton"))->Fill(tofRes[iS]);
		  break;
		default:
		  break;
		}
	      }
	    }
	  }
	}//if TOF status

	if(pidHF && pidHF->CheckStatus(track,"TPC")){ 

	  Double_t TPCp=pid->GetTPCmomentum();
	  Double_t TPCsignal=pid->GetTPCsignal();
	  ((TH1F*)fOutputPID->FindObject("hTPCsig"))->Fill(TPCsignal);
	  ((TH1F*)fOutputPID->FindObject("hTPCsigvsp"))->Fill(TPCp,TPCsignal);
	  //if (pidHF->IsKaonRaw(track, "TOF"))
	  if(!oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaK"))->Fill(TPCp,respF->NumberOfSigmasTPC(track,AliPID::kKaon));
	  if (oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaK"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kKaon));
	  //if (pidHF->IsPionRaw(track, "TOF"))
	  if(oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaPion"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kPion));
	  if(!oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaPion"))->Fill(TPCp,respF->NumberOfSigmasTPC(track,AliPID::kPion));
	  //if (pidHF->IsProtonRaw(track,"TOF"))
	  if(oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaProton"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kProton));
	  if(!oldPID) ((TH2F*)fOutputPID->FindObject("hTPCsigmaProton"))->Fill(TPCp,respF->NumberOfSigmasTPC(track,AliPID::kProton));
	}//if TPC status
      } //end PID histograms

      Int_t nclsTot=0,nclsSPD=0;

      //check clusters of the tracks
      if(fOnOff[0]){

	((TH1F*)fOutputTrack->FindObject("hnLayerITS"))->Fill(-1);
	for(Int_t l=0;l<6;l++) {
	  if(TESTBIT(track->GetITSClusterMap(),l)) {
	    ((TH1F*)fOutputTrack->FindObject("hnLayerITS"))->Fill(l);
	    nclsTot++; if(l<2) nclsSPD++;
	  }
	}
	((TH1F*)fOutputTrack->FindObject("hnClsITS"))->Fill(nclsTot);
	((TH1F*)fOutputTrack->FindObject("hnClsSPD"))->Fill(nclsSPD);
	if(track->Pt()>0.3 &&
	   TMath::Abs(track->Eta())<0.8 &&
	   track->GetStatus()&AliESDtrack::kITSrefit &&
	   track->GetStatus()&AliESDtrack::kTPCrefit &&
	   nclsSPD>0){
	  ((TH1F*)fOutputTrack->FindObject("hnClsITSselTr"))->Fill(nclsTot);
	}
	if(!(track->GetStatus()&AliESDtrack::kTPCin) && track->GetStatus()&AliESDtrack::kITSrefit && !(track->GetStatus()&AliESDtrack::kITSpureSA)){//tracks retrieved in the ITS and not reconstructed in the TPC
	  ((TH1F*)fOutputTrack->FindObject("hnClsITS-SA"))->Fill(nclsTot);
	  ((TH1F*)fOutputTrack->FindObject("hnLayerITS"))->Fill(-1);
	  for(Int_t l=0;l<6;l++) {
	    if(TESTBIT(track->GetITSClusterMap(),l)) {
	      ((TH1F*)fOutputTrack->FindObject("hnLayerITSsa"))->Fill(l);
	    }
	  }
	}
	Int_t label=0;
	if(fReadMC){
	  label=track->GetLabel();
	  if (label<0)fNEntries->Fill(8);
	  else fNEntries->Fill(9); 
	}

	if(isSimpleMode){

	  if (track->Pt()>0.3 &&
	      track->GetStatus()&AliESDtrack::kTPCrefit &&
	      track->GetStatus()&AliESDtrack::kITSrefit &&
	      /*nclsTot>3 &&*/
	      nclsSPD>0) {//count good tracks

	    
	    if(fReadMC && label<0) {
	      ((TH1F*)fOutputTrack->FindObject("hptFakeTr"))->Fill(track->Pt());
	      isFakeTrack++;	
	    } else {
	      ((TH1F*)fOutputTrack->FindObject("hptGoodTr"))->Fill(track->Pt());
	      isGoodTrack++;
	    }
	  }
	} //simple mode: no IsSelected on tracks: use "manual" cuts   
      } //fill track histos
    } //end loop on tracks

      //fill once per event
    if(fOnOff[0]){
      if (fReadMC) ((TH1F*)fOutputTrack->FindObject("hdistrFakeTr"))->Fill(isFakeTrack);
      ((TH1F*)fOutputTrack->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
    }

    if(!isSimpleMode){
      // loop over candidates
      Int_t nCand = arrayProng->GetEntriesFast();
      Int_t ndaugh=3;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi) ndaugh=2;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpipipi) ndaugh=4;

      for (Int_t iCand = 0; iCand < nCand; iCand++) {
	AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
	if(d->GetSelectionMap()) {
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi && !d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) continue; //skip the D0 from Dstar
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kDplustoKpipi && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) continue; //skip the 3 prong !D+
	}

	if(fReadMC){ 
	  Int_t labD = d->MatchToMC(pdg,mcArray,ndaugh,pdgdaughters);
	  if(labD>=0){
	    AliAODMCParticle *partD = (AliAODMCParticle*)mcArray->At(labD);
	    Int_t label=partD->GetMother();
	    AliAODMCParticle *mot = (AliAODMCParticle*)mcArray->At(label);
	    while(label>=0){//get first mother
	      mot = (AliAODMCParticle*)mcArray->At(label);
	      label=mot->GetMother();
	    }
	    Int_t pdgMotCode = mot->GetPdgCode();
	
	    if(TMath::Abs(pdgMotCode)==4) fNEntries->Fill(6); //from primary charm
	    if(TMath::Abs(pdgMotCode)==5) fNEntries->Fill(7); //from beauty

	  }
	}//end MC
	else fNEntries->Fill(6); //count the candidates (data)

	for(Int_t id=0;id<ndaugh;id++){

	  //other histograms to be filled when the cut object is given
	  AliAODTrack* track=(AliAODTrack*)d->GetDaughter(id);

	  //track quality

	  if (fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(pdg)) && fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod)) {
	    
	    Int_t label=0;
	    if(fReadMC)label=track->GetLabel();
	    if(fOnOff[0]){
	      
	      if(fReadMC && label<0) {
		isFakeTrack++;
		((TH1F*)fOutputTrack->FindObject("hptFakeTr"))->Fill(track->Pt());
	   
		((TH1F*)fOutputTrack->FindObject("hd0f"))->Fill(d->Getd0Prong(id));
	      } else {
		isGoodTrack++;
		((TH1F*)fOutputTrack->FindObject("hptGoodTr"))->Fill(track->Pt());
		((TH1F*)fOutputTrack->FindObject("hd0"))->Fill(d->Getd0Prong(id));
		Double_t d0rphiz[2],covd0[3];
		Bool_t isDCA=track->PropagateToDCA(aod->GetPrimaryVertex(),aod->GetMagneticField(),9999.,d0rphiz,covd0);
		if(isDCA) ((TH1F*)fOutputTrack->FindObject("hd0z"))->Fill(d0rphiz[1]);
	      }
	    }
	    if (fCuts->IsSelected(d,AliRDHFCuts::kAll,aod) && fOnOff[1]){
	  
	      AliAODPid *pid = track->GetDetPid();
	      Double_t times[5];
	      pid->GetIntegratedTimes(times);
	      if(pidHF && pidHF->CheckStatus(track,"TOF")) ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptimeAC"))->Fill(track->P(),pid->GetTOFsignal()-times[AliPID::kKaon]);
	      if(pidHF && pidHF->CheckStatus(track,"TPC")) ((TH2F*)fOutputPID->FindObject("hTPCsigvspAC"))->Fill(pid->GetTPCmomentum(),pid->GetTPCsignal());

	      fNEntries->Fill(3);
	    } //end analysis cuts
	  } //end acceptance and track cuts
	} //end loop on tracks in the candidate
      } //end loop on candidates
      if(fOnOff[0]){
	if(fReadMC) ((TH1F*)fOutputTrack->FindObject("hdistrFakeTr"))->Fill(isFakeTrack);
	((TH1F*)fOutputTrack->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
      }
    }
  } //end if on pid or track histograms

  delete tpcres;
  delete [] pdgdaughters;
  PostData(1,fNEntries);
  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) PostData(5,fOutputCounters);
  //Post data 6 done in case of centrality on   

}

//____________________________________________________________________________
void AliAnalysisTaskSEHFQA::Terminate(Option_t */*option*/){
  //terminate analysis

  fNEntries = dynamic_cast<TH1F*>(GetOutputData(1));
  if(!fNEntries){
    printf("ERROR: %s not available\n",GetOutputSlot(1)->GetContainer()->GetName());
    return;
  }

  fOutputPID = dynamic_cast<TList*> (GetOutputData(2));
  if (!fOutputPID && fOnOff[1]) {     
    printf("ERROR: %s not available\n",GetOutputSlot(2)->GetContainer()->GetName());
    return;
  }

  fOutputTrack = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputTrack && fOnOff[0]) {     
    printf("ERROR: %s not available\n",GetOutputSlot(3)->GetContainer()->GetName());
    return;
  }

}

