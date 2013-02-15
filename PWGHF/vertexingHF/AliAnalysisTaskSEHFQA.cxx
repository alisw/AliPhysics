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
#include <TH3F.h>
#include <TProfile2D.h>
#include <TDatabasePDG.h>

#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliESDtrack.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
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
#include "AliRDHFCutsLctoV0.h"
#include "AliInputEventHandler.h"

#include "AliFlowEvent.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowVector.h"

#include "AliAnalysisTaskSEHFQA.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskSEHFQA)

//____________________________________________________________________________

AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA():AliAnalysisTaskSE(),
  fNEntries(0x0),
  fOutputPID(0x0),
  fOutputTrack(0x0),
  fOutputCounters(0x0),
  fOutputCheckCentrality(0x0),
  fOutputEvSelection(0x0),
  fOutputFlowObs(0x0),
  fDecayChannel(AliAnalysisTaskSEHFQA::kD0toKpi),
  fCuts(0x0),
  fFlowEvent(0x0),
  fRFPcuts(0x0),
  fEstimator(AliRDHFCuts::kCentTRK),
  fReadMC(kFALSE),
  fSimpleMode(kFALSE),
  fUseSelectionBit(kTRUE),
  fOnOff()
{
  //default constructor
  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;
  fOnOff[4]=kTRUE;
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
  fOutputFlowObs(0x0),
  fDecayChannel(ch),
  fCuts(0x0),
  fFlowEvent(0x0),
  fRFPcuts(0x0),
  fEstimator(AliRDHFCuts::kCentTRK),
  fReadMC(kFALSE),
  fSimpleMode(kFALSE),
  fUseSelectionBit(kTRUE),
  fOnOff()
{
  //constructor

  //SetCutObject(cuts);
  fCuts=cuts;

  fOnOff[0]=kTRUE;
  fOnOff[1]=kTRUE;
  fOnOff[2]=kTRUE;
  fOnOff[3]=kTRUE;
  fOnOff[4]=kTRUE;

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
  case kLambdactoV0:
    DefineOutput(4,AliRDHFCutsLctoV0::Class());  //My private output
    break;
  }
  if (fOnOff[2]) {
    // Output slot #5 writes into a TList container (AliCounterCollection)
    DefineOutput(5,TList::Class());  //My private output
    // Output slot #6 writes into a TList container (TH1F)
    DefineOutput(6,TList::Class());  //My private output
  }

  if(fOnOff[3]) DefineOutput(7,TList::Class());  //My private output
  if(fOnOff[4]) DefineOutput(8,TList::Class());  //My private output

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

  if(fOnOff[4]) {
    delete fOutputFlowObs;
    delete fFlowEvent;
  }
}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::Init(){

  //initialization
  if(fDebug > 1) printf("AnalysisTaskSEHFQA::Init() \n");
  AliRDHFCuts *copycut = 0x0;

  switch(fDecayChannel){
  case 0:
    {
      copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fCuts)));
    }
    break;
  case 1:
    {
      copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
    }
    break;
  case 2:
    {
      copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
    }
    break;
  case 3:
    {
      copycut=new AliRDHFCutsDstoKKpi(*(static_cast<AliRDHFCutsDstoKKpi*>(fCuts)));
    }
    break;
  case 4:
    {
      copycut=new AliRDHFCutsD0toKpipipi(*(static_cast<AliRDHFCutsD0toKpipipi*>(fCuts)));
    }
    break;
  case 5:
    {
      copycut=new AliRDHFCutsLctopKpi(*(static_cast<AliRDHFCutsLctopKpi*>(fCuts)));
    }
    break;
  case kLambdactoV0:
    {
      copycut=new AliRDHFCutsLctoV0(*(static_cast<AliRDHFCutsLctoV0*>(fCuts)));
    }
    break;
  default:
    AliFatal("Bad initialization for the decay channe - Exiting...");
    break;
  }  

  const char* nameoutput=GetOutputSlot(4)->GetContainer()->GetName();
  if (copycut){
    copycut->SetName(nameoutput);
    
    // Post the data
    PostData(4,copycut);
  }else{
    AliFatal("Failing initializing AliRDHFCuts object - Exiting...");
  }	

  return;

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
  fNEntries->GetXaxis()->SetBinLabel(6,"N candidates");
  if(fReadMC){
    fNEntries->GetXaxis()->SetBinLabel(7,"MC Cand from c");
    fNEntries->GetXaxis()->SetBinLabel(8,"MC Cand from b");
    fNEntries->GetXaxis()->SetBinLabel(9,"N fake Trks");
    fNEntries->GetXaxis()->SetBinLabel(10,"N true Trks");
  }

  fNEntries->GetXaxis()->SetNdivisions(1,kFALSE);

  //PID
  if(fOnOff[1]){
    fOutputPID=new TList();
    fOutputPID->SetOwner();
    fOutputPID->SetName(GetOutputSlot(2)->GetContainer()->GetName());

    //TOF pid
    TH1F* hTOFflags=new TH1F("hTOFflags","TOF flags",7,-0.5,6.5);
    hTOFflags->SetMinimum(0.);
    hTOFflags->GetXaxis()->SetBinLabel(1,"All Tracks");
    hTOFflags->GetXaxis()->SetBinLabel(2,"kTPCout");
    hTOFflags->GetXaxis()->SetBinLabel(3,"kTOFout");
    hTOFflags->GetXaxis()->SetBinLabel(4,"kTIME");
    hTOFflags->GetXaxis()->SetBinLabel(5,"kTOFpid");
    hTOFflags->GetXaxis()->SetBinLabel(6,"kTOFmismatch");
    hTOFflags->GetXaxis()->SetBinLabel(7,"kDetPidOK");

    TString hname="hTOFsig";
    TH1F* hTOFsig=new TH1F(hname.Data(),"Distribution of TOF signal;TOF time [ps];Entries", 100, -2.e3,40.e3);

    hname="hTOFstartTimeMask";
    TH1F* hTOFstartTimeMask=new TH1F(hname.Data(),"TOF start time mask; Mask ;Entries", 8, -0.5,7.5);
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(1,"FILL");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(2,"TOF");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(3,"T0A");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(4,"TOF.and.T0A");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(5,"T0C");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(6,"TOF.and.T0C");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(7,"T0AC");
    hTOFstartTimeMask->GetXaxis()->SetBinLabel(8,"TOF.and.T0AC");

    hname="hTOFstartTimeRes";
    TH1F* hTOFstartTimeRes=new TH1F(hname.Data(),"TOF start time resolution; Resolution (ps) ;Entries", 100, 0.,300.);

    hname="hTOFstartTimeDistrib";
    TH1F* hTOFstartTimeDistrib=new TH1F(hname.Data(),"TOF start time distribution; Start time ;Entries", 400, -1000.,1000.);

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
    fOutputPID->Add(hTOFstartTimeMask);
    fOutputPID->Add(hTOFstartTimeRes);
    fOutputPID->Add(hTOFstartTimeDistrib);
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

    if(fReadMC){
      //TOF
      hname="hTOFsigmaMCKSigPid";
      TH2F* hTOFsigmaMCKSigPid=new TH2F(hname.Data(),"(TOFsignal-timeK)/tofSigPid;p[GeV/c];(TOFsignal-timeK)/tofSigPid",500,0.,10.,400,-20,20);

      hname="hTOFsigmaMCPionSigPid";
      TH2F* hTOFsigmaMCPionSigPid=new TH2F(hname.Data(),"(TOFsignal-time#pi)/tofSigPid;p[GeV/c];(TOFsignal-time#pi)/tofSigPid",500,0.,10.,400,-20,20);

      hname="hTOFsigmaMCProtonSigPid";
      TH2F* hTOFsigmaMCProtonSigPid=new TH2F(hname.Data(),"(TOFsignal-timep)/tofSigPid;p[GeV/c];(TOFsignal-time p)/tofSigPid",500,0.,10.,400,-20,20);

      //TPC
      hname="hTPCsigmaMCK";
      TH2F* hTPCsigmaMCK=new TH2F(hname.Data(),"TPC Sigma for K as a function of momentum;p[GeV/c];Sigma Kaon",500,0.,10.,400,-20,20);

      hname="hTPCsigmaMCPion";
      TH2F* hTPCsigmaMCPion=new TH2F(hname.Data(),"TPC Sigma for #pi as a function of momentum;p[GeV/c];Sigma #pi",500,0.,10.,400,-20,20);

      hname="hTPCsigmaMCProton";
      TH2F* hTPCsigmaMCProton=new TH2F(hname.Data(),"TPC Sigma for proton as a function of momentum;p[GeV/c];Sigma Proton",500,0.,10.,400,-20,20);

      fOutputPID->Add(hTOFsigmaMCKSigPid);
      fOutputPID->Add(hTOFsigmaMCPionSigPid);
      fOutputPID->Add(hTOFsigmaMCProtonSigPid);
      fOutputPID->Add(hTPCsigmaMCK);
      fOutputPID->Add(hTPCsigmaMCPion);
      fOutputPID->Add(hTPCsigmaMCProton);

    }
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
    hnLayerITS->GetXaxis()->SetBinLabel(2,"SPDin");
    hnLayerITS->GetXaxis()->SetBinLabel(3,"SPDout");
    hnLayerITS->GetXaxis()->SetBinLabel(4,"SDDin");
    hnLayerITS->GetXaxis()->SetBinLabel(5,"SDDout");
    hnLayerITS->GetXaxis()->SetBinLabel(6,"SSDin");
    hnLayerITS->GetXaxis()->SetBinLabel(7,"SSDout");

    hname="hnLayerITSsa";
    TH1F* hnLayerITSsa=new TH1F(hname.Data(),"Number of tracks with point in layer;ITS layer;",7,-1.5,5.5);
    hnLayerITSsa->GetXaxis()->SetBinLabel(1,"n tracks");
    hnLayerITSsa->GetXaxis()->SetBinLabel(2,"SPDin");
    hnLayerITSsa->GetXaxis()->SetBinLabel(3,"SPDout");
    hnLayerITSsa->GetXaxis()->SetBinLabel(4,"SDDin");
    hnLayerITSsa->GetXaxis()->SetBinLabel(5,"SDDout");
    hnLayerITSsa->GetXaxis()->SetBinLabel(6,"SSDin");
    hnLayerITSsa->GetXaxis()->SetBinLabel(7,"SSDout");
   
    hname="hnClsSPD";
    TH1F* hnClsSPD=new TH1F(hname.Data(),"Distribution of number of SPD clusters;nSPDcls;Entries",3,-0.5,2.5);

    hname="hptGoodTr";
    TH1F* hptGoodTr=new TH1F(hname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
    hptGoodTr->SetTitleOffset(1.3,"Y");

    if(!fSimpleMode){
      hname="hptGoodTrFromDaugh";
      TH1F* hptGoodTrFromDaugh=new TH1F(hname.Data(),"Pt distribution of 'good' candidate's daughters;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      hptGoodTrFromDaugh->SetTitleOffset(1.3,"Y");
      fOutputTrack->Add(hptGoodTrFromDaugh);
    }

    hname="hdistrGoodTr";
    TH1F* hdistrGoodTr=new TH1F(hname.Data(),"Distribution of number of 'good' candidate's daughters per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
    hdistrGoodTr->SetTitleOffset(1.3,"Y");

    hname="hdistrSelTr";
    TH1F* hdistrSelTr=new TH1F(hname.Data(),"Distribution of number of Selected tracks per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
    hdistrSelTr->SetTitleOffset(1.3,"Y");

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
    fOutputTrack->Add(hdistrSelTr);
    fOutputTrack->Add(hd0);
    fOutputTrack->Add(hd0z);

    if(fReadMC){
      hname="hdistrFakeTr";
      TH1F* hdistrFakeTr=new TH1F(hname.Data(),"Distribution of number of fake tracks per event;no.fake-tracks/ev;Entries",4000,-0.5,3999.5);
      hdistrFakeTr->SetTitleOffset(1.3,"Y");

      hname="hd0f";
      TH1F* hd0f=new TH1F(hname.Data(),"Impact parameter distribution of fake tracks;d_{0}[cm];Entries/10^{3} cm",200,-0.1,0.1);

      hname="hptFakeTr";
      TH1F* hptFakeTr=new TH1F(hname.Data(),"Pt distribution of fake tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
      hptFakeTr->SetTitleOffset(1.3,"Y");
      if(!fSimpleMode){
	hname="hptFakeTrFromDaugh";
	TH1F* hptFakeTrFromDaugh=new TH1F(hname.Data(),"Pt distribution of fake tracks from daughters;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
	hptFakeTrFromDaugh->SetTitleOffset(1.3,"Y");
	fOutputTrack->Add(hptFakeTrFromDaugh);
      }

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
    TH2F* hMultvsPercentile=new TH2F(hname.Data(),"Multiplicity vs Percentile;multiplicity;percentile",10000,-0.5,9999.5,240,-10.,110);

    hname="hntrklvsPercentile";
    TH2F* hntrklvsPercentile=new TH2F(hname.Data(),"N tracklets vs Percentile;ntracklets;percentile",5000,-0.5,4999.5,240,-10.,110);

    hname="hnTPCTracksvsPercentile";
    TH2F* hnTPCTracksvsPercentile=new TH2F(hname.Data(),"N TPC tracks vs Percentile;nTPCTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hnTPCITSTracksvsPercentile";
    TH2F* hnTPCITSTracksvsPercentile=new TH2F(hname.Data(),"N TPC+ITS tracks vs Percentile;nTPCITSTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hnTPCITS1SPDTracksvsPercentile";
    TH2F* hnTPCITS1SPDTracksvsPercentile=new TH2F(hname.Data(),"N TPC+ITS+1SPD tracks vs Percentile;nTPCITS1SPDTracks;percentile",5000,-0.5,9999.5,240,-10.,110);

    hname="hV0MultiplicityPercentile";
    TH2F*hV0MultiplicityPercentile = new TH2F(hname.Data(),"V0 Multiplicity vs Percentile;V0 multiplicity;percentile",1000,-0.5,9999.5,120,-10.,110);

    hname="hV0MultiplicityNtrackletsIn";
    TH2F*hV0MultiplicityNtrackletsIn = new TH2F(hname.Data(),"V0 Multiplicity vs Number of tracklets in the CC;V0 multiplicity;percentile",1000,-0.5,9999.5,5000,-0.5,4999.5);

    hname="hStdPercentileSPDPercentile";
    TH2F* hStdPercentileSPDPercentile = new TH2F(hname.Data(),"Std estimator Percentile Vs SPD Percentile;Std estimator percentile;SPD percentile",120,-10.,110,120,-10.,110);

    fOutputCheckCentrality->Add(hNtrackletsIn);
    fOutputCheckCentrality->Add(hNtrackletsOut);
    fOutputCheckCentrality->Add(hMultIn);
    fOutputCheckCentrality->Add(hMultOut);
    fOutputCheckCentrality->Add(hMultvsPercentile);
    fOutputCheckCentrality->Add(hntrklvsPercentile);
    fOutputCheckCentrality->Add(hnTPCTracksvsPercentile);
    fOutputCheckCentrality->Add(hnTPCITSTracksvsPercentile);
    fOutputCheckCentrality->Add(hnTPCITS1SPDTracksvsPercentile);
    fOutputCheckCentrality->Add(hV0MultiplicityPercentile);
    fOutputCheckCentrality->Add(hV0MultiplicityNtrackletsIn);
    fOutputCheckCentrality->Add(hStdPercentileSPDPercentile);

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

    TH1F* hxvtx=new TH1F("hxvtx", "Distribution of x_{VTX};x_{VTX} [cm];Entries",100,-1,1);
    TH1F* hyvtx=new TH1F("hyvtx", "Distribution of y_{VTX};y_{VTX} [cm];Entries",100,-1,1);
    TH1F* hzvtx=new TH1F("hzvtx", "Distribution of z_{VTX};z_{VTX} [cm];Entries",100,-30,30);
    TH1F* hxvtxSelEv=new TH1F("hxvtxSelEv", "Distribution of x_{VTX} Selected Ev;x_{VTX} [cm];Entries",100,-1,1);
    TH1F* hyvtxSelEv=new TH1F("hyvtxSelEv", "Distribution of y_{VTX} Selected Ev;y_{VTX} [cm];Entries",100,-1,1);
    TH1F* hzvtxSelEv=new TH1F("hzvtxSelEv", "Distribution of z_{VTX} Selected Ev;z_{VTX} [cm];Entries",100,-30,30);
    TH1F* hWhichVert=new TH1F("hWhichVert","Vertex Type",4,-1.5,2.5);
    hWhichVert->GetXaxis()->SetBinLabel(1,"Not found");
    hWhichVert->GetXaxis()->SetBinLabel(2,"Track");
    hWhichVert->GetXaxis()->SetBinLabel(3,"SPD-3D");
    hWhichVert->GetXaxis()->SetBinLabel(4,"SPD-z");
    TH1F* hWhichVertSelEv=new TH1F("hWhichVertSelEv","Vertex Type",4,-1.5,2.5);
    hWhichVertSelEv->GetXaxis()->SetBinLabel(1,"Not found");
    hWhichVertSelEv->GetXaxis()->SetBinLabel(2,"Track");
    hWhichVertSelEv->GetXaxis()->SetBinLabel(3,"SPD-3D");
    hWhichVertSelEv->GetXaxis()->SetBinLabel(4,"SPD-z");

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
    trigCounter->AddRubric("triggerType","All/Any/MB/Cent/SemiCent/EMCAL/MUON/NoPhysSelMUON/NoPhysSelEvNot7/NoPhysSelCMUP1/NoPhysSelMB/NoPhysSelCent/NoPhysSelSemiCent/CINT7");
    trigCounter->Init();

    TH1F* hWhyEvRejected=new TH1F("hWhyEvRejected", "Why Event rejected",6,-0.5,5.5);

    hWhyEvRejected->GetXaxis()->SetBinLabel(1,"pileup");
    hWhyEvRejected->GetXaxis()->SetBinLabel(2,"centrality");
    hWhyEvRejected->GetXaxis()->SetBinLabel(3,"Vertex (more reasons)");
    hWhyEvRejected->GetXaxis()->SetBinLabel(4,"trigger");
    hWhyEvRejected->GetXaxis()->SetBinLabel(5,"z vertex out");
    hWhyEvRejected->GetXaxis()->SetBinLabel(6,"physics sel");


    fOutputEvSelection->Add(evselection);
    fOutputEvSelection->Add(hxvtx);
    fOutputEvSelection->Add(hyvtx);
    fOutputEvSelection->Add(hzvtx);
    fOutputEvSelection->Add(hxvtxSelEv);
    fOutputEvSelection->Add(hyvtxSelEv);
    fOutputEvSelection->Add(hzvtxSelEv);
    fOutputEvSelection->Add(hWhichVert);
    fOutputEvSelection->Add(hWhichVertSelEv);
    fOutputEvSelection->Add(hTrigCent);
    fOutputEvSelection->Add(hTrigMul);
    fOutputEvSelection->Add(hTrigCentSel);
    fOutputEvSelection->Add(trigCounter);
    fOutputEvSelection->Add(hWhyEvRejected);

  }
  if(fOnOff[4]){ // FLOW OBSERVABLES
    fOutputFlowObs=new TList();
    fOutputFlowObs->SetOwner();
    fOutputFlowObs->SetName(GetOutputSlot(8)->GetContainer()->GetName());

    fFlowEvent = new AliFlowEvent(3000);
    fRFPcuts = new AliFlowTrackCuts("rfpCuts");

    TH2F *hFEvents = new TH2F("hFlowEvents","FlowEvent Selection",7,0,7,7,-10,60);
    hFEvents->GetXaxis()->SetBinLabel(1,"REACHED");
    hFEvents->GetXaxis()->SetBinLabel(2,"TRIGGERED");
    hFEvents->GetXaxis()->SetBinLabel(3,"kMB");
    hFEvents->GetXaxis()->SetBinLabel(4,"kCent");
    hFEvents->GetXaxis()->SetBinLabel(5,"kSemiC");
    hFEvents->GetXaxis()->SetBinLabel(6,"Triggered + vtx cut");
    hFEvents->GetXaxis()->SetBinLabel(7,"UnexpectedBehaviour");
    fOutputFlowObs->Add(hFEvents);

    TProfile2D *hQ[3];
    TH2F *hAngleQ[3];
    TH3F *hPhiEta[3];
    TString ref[3] = {"FB1","FB128","VZE"};
    Int_t etabin[3] = {40,40,20};
    Int_t etamax[3] = { 1, 1, 5};
    for(Int_t i=0; i<3; ++i) {
      hQ[i]= new TProfile2D( Form("h%s_Q",ref[i].Data()),
			     Form("Q_{2} components for %s",ref[i].Data()),
			     4,0,4,12,0,60,"s");
      hQ[i]->GetXaxis()->SetBinLabel(1,"Qx^{-}");
      hQ[i]->GetXaxis()->SetBinLabel(2,"Qy^{-}");
      hQ[i]->GetXaxis()->SetBinLabel(3,"Qx^{+}");
      hQ[i]->GetXaxis()->SetBinLabel(4,"Qy^{+}");
      hQ[i]->GetYaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(hQ[i]);

      hAngleQ[i] = new TH2F( Form("h%s_AngleQ",ref[i].Data()),
			     Form("#Psi_{2} for %s",ref[i].Data()),
			     72,0,TMath::Pi(),12,0,60);
      hAngleQ[i]->GetXaxis()->SetTitle( Form("#Psi_{2}^{%s}",ref[i].Data()) );
      hAngleQ[i]->GetYaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(hAngleQ[i]);

      hPhiEta[i] = new TH3F( Form("h%s_PhiEta",ref[i].Data()),
			     Form("Eta vs Phi for %s",ref[i].Data()),
			     144,0,TMath::TwoPi(),etabin[i],-1.0*etamax[i],+1.0*etamax[i],12,0,60);
      hPhiEta[i]->GetXaxis()->SetTitle("Phi");
      hPhiEta[i]->GetYaxis()->SetTitle("Eta");
      hPhiEta[i]->GetZaxis()->SetTitle("Centrality");
      fOutputFlowObs->Add(hPhiEta[i]);

    }
    TH3F *hTPCVZE_AngleQ = new TH3F("hTPCVZE_AngleQ","#Psi_{2}^{VZERO} vs #Psi_{2}^{TPC}",   72,0,TMath::Pi(),72,0,TMath::Pi(),12,0,60);
    hTPCVZE_AngleQ->GetXaxis()->SetTitle("#Psi_{2}^{TPC}");
    hTPCVZE_AngleQ->GetYaxis()->SetTitle("#Psi_{2}^{VZE}");
    hTPCVZE_AngleQ->GetZaxis()->SetTitle("Centrality");
    fOutputFlowObs->Add(hTPCVZE_AngleQ);

    TH2F *hCentVsMultRPS = new TH2F("hCentVsMultRPS", " Centrality Vs. Multiplicity RPs",5000, 0, 5000.,12,0,60 );
    hCentVsMultRPS->GetXaxis()->SetTitle("Multiplicity RPs");
    hCentVsMultRPS->GetYaxis()->SetTitle("Centrality");
    fOutputFlowObs->Add(hCentVsMultRPS);
  }

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();
  if (fCuts->GetIsUsePID() && fDecayChannel==kLambdactoV0) {
    fCuts->GetPidHF()->SetPidResponse(pidResp);
    AliRDHFCutsLctoV0* lccuts=dynamic_cast<AliRDHFCutsLctoV0*>(fCuts);
    if(lccuts){
      lccuts->GetPidV0pos()->SetPidResponse(pidResp);
      lccuts->GetPidV0neg()->SetPidResponse(pidResp);
      fCuts->GetPidHF()->SetOldPid(kFALSE);
      lccuts->GetPidV0pos()->SetOldPid(kFALSE);
      lccuts->GetPidV0neg()->SetOldPid(kFALSE);
    }
  }

  // Post the data
  PostData(1,fNEntries);

  if(fOnOff[1]) PostData(2,fOutputPID);
  if(fOnOff[0]) PostData(3,fOutputTrack);
  PostData(4,fCuts);
  if(fOnOff[2]) PostData(5,fOutputCounters);
  if(fOnOff[3]) PostData(7,fOutputEvSelection);
  if(fOnOff[4]) PostData(8,fOutputFlowObs);

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
      case kLambdactoV0:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("CascadesHF");
	pdg=4122;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=2212;//p
	  pdgdaughters[1]=211;//pi
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
      case kLambdactoV0:
	arrayProng=(TClonesArray*)aod->GetList()->FindObject("CascadesHF");
	pdg=4122;
	if(fReadMC){
	  pdgdaughters =new Int_t[3];
	  pdgdaughters[0]=2212;//p
	  pdgdaughters[1]=211;//pi
	  pdgdaughters[2]=211;//pi
	}
	break; 
    }
  }
  Bool_t isSimpleMode=fSimpleMode;
  if(!arrayProng) {
    AliInfo("Branch not found! The output will contain only track related histograms\n");
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
  TString trigClass=aod->GetFiredTriggerClasses();
  Int_t nAODtracks=aod->GetNTracks();
  Int_t nSelTracksTPCOnly=0;
  Int_t nSelTracksTPCITS=0;
  Int_t nSelTracksTPCITS1SPD=0;

  for (Int_t k=0;k<nAODtracks;k++){
    AliAODTrack* track=aod->GetTrack(k);
    if(track->GetID()<0) continue;
    Int_t nclsTot=0,nclsSPD=0;
    for(Int_t l=0;l<6;l++) {
      if(TESTBIT(track->GetITSClusterMap(),l)) {
	nclsTot++; if(l<2) nclsSPD++;
      }
    }
    UShort_t nTPCClus=track->GetTPCClusterMap().CountBits();
    if(TMath::Abs(track->Eta())<0.8 && nTPCClus>=70 && track->GetStatus()&AliESDtrack::kTPCrefit){
      if(track->TestFilterBit(1))  nSelTracksTPCOnly++;
      if(track->GetStatus()&AliESDtrack::kITSrefit){
	nSelTracksTPCITS++;
	if(nclsSPD>0) nSelTracksTPCITS1SPD++;      
      }
    }
  }

  if(fOnOff[4]) {
    FillFlowObs(aod);
    PostData(8,fOutputFlowObs);
  }
  if(fOnOff[3]){
    TH2F* hTrigC=(TH2F*)fOutputEvSelection->FindObject("hTrigCent");
    TH2F* hTrigM=(TH2F*)fOutputEvSelection->FindObject("hTrigMul");
    AliCounterCollection* trigCount=(AliCounterCollection*)fOutputEvSelection->FindObject("trigCounter");

    hTrigC->Fill(-1.,centrality);
    hTrigM->Fill(-1.,multiplicity);
    trigCount->Count(Form("triggerType:All/Run:%d",runNumber));
    if(evSelMask==0){
      if(aod->GetEventType()!=7){
	trigCount->Count(Form("triggerType:NoPhysSelEvNot7/Run:%d",runNumber));
      }else if(trigClass.Contains("CMUP1")){
	trigCount->Count(Form("triggerType:NoPhysSelCMUP1/Run:%d",runNumber));
      }else if(trigClass.Contains("MUON")){
	trigCount->Count(Form("triggerType:NoPhysSelMUON/Run:%d",runNumber));
      }else if(trigClass.Contains("CPBI2_B1-B") || trigClass.Contains(" CPBI2WU_B1-B")){
	trigCount->Count(Form("triggerType:NoPhysSelMB/Run:%d",runNumber));
      }else if(trigClass.Contains("CCENT") || trigClass.Contains("CVHN")){
	trigCount->Count(Form("triggerType:NoPhysSelCent/Run:%d",runNumber));
      }else if(trigClass.Contains("CSEMI") || trigClass.Contains("CVLN")){
	trigCount->Count(Form("triggerType:NoPhysSelSemiCent/Run:%d",runNumber));
      }
    }
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
      trigCount->Count(Form("triggerType:CINT7/Run:%d",runNumber));
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
      trigCount->Count(Form("triggerType:MUON/Run:%d",runNumber));
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

  Double_t pos[3],cov[6];
  vtx1->GetXYZ(pos);
  vtx1->GetCovarianceMatrix(cov);
  const AliESDVertex vESD(pos,cov,100.,100);

  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) fNEntries->Fill(4);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  //TString trigclass=aod->GetFiredTriggerClasses();
  //if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNEntries->Fill(5); //tmp




  Bool_t evSelbyCentrality=kTRUE,evSelected=kTRUE,evSelByVertex=kTRUE,evselByPileup=kTRUE,evSelByPS=kTRUE;

  TH1F* hWhyEvRejected=(TH1F*)fOutputEvSelection->FindObject("hWhyEvRejected");

  //select event
  if(!fCuts->IsEventSelected(aod)) {
    evSelected=kFALSE;
    if(fCuts->IsEventRejectedDueToPileupSPD()) {hWhyEvRejected->Fill(1); evselByPileup=kFALSE;}// rejected for pileup
    if(fCuts->IsEventRejectedDueToCentrality()) {hWhyEvRejected->Fill(2); evSelbyCentrality=kFALSE; //rejected by centrality
    }
    if(fCuts->IsEventRejectedDueToNotRecoVertex() ||
       fCuts->IsEventRejectedDueToVertexContributors() ||
       fCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion()){ 
      evSelByVertex=kFALSE; 
      hWhyEvRejected->Fill(3);
    }
    if(fCuts->IsEventRejectedDueToTrigger()) hWhyEvRejected->Fill(4);//tmp
    if(fCuts->IsEventRejectedDueToZVertexOutsideFiducialRegion() && fOnOff[3]) {((AliCounterCollection*)fOutputEvSelection->FindObject("evselection"))->Count(Form("evnonsel:zvtx/Run:%d",runNumber)); hWhyEvRejected->Fill(5);
    }
    if(fCuts->IsEventRejectedDuePhysicsSelection()) { evSelByPS=kFALSE;hWhyEvRejected->Fill(6); }
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
      AliAODVZERO *vzeroAOD = (AliAODVZERO*)aod->GetVZEROData();
      Float_t vzeroMult = vzeroAOD->GetMTotV0A() +  vzeroAOD->GetMTotV0C();
      AliCentrality *aodcent = aod->GetCentrality();
      Float_t spdCentf = aodcent->GetCentralityPercentile("CL1");
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
      ((TH2F*)fOutputCheckCentrality->FindObject("hntrklvsPercentile"))->Fill(aod->GetTracklets()->GetNumberOfTracklets(),stdCentf);
      ((TH2F*)fOutputCheckCentrality->FindObject("hnTPCTracksvsPercentile"))->Fill(nSelTracksTPCOnly,stdCentf);
      ((TH2F*)fOutputCheckCentrality->FindObject("hnTPCITSTracksvsPercentile"))->Fill(nSelTracksTPCITS,stdCentf);
      ((TH2F*)fOutputCheckCentrality->FindObject("hnTPCITS1SPDTracksvsPercentile"))->Fill(nSelTracksTPCITS1SPD,stdCentf);
      ((TH2F*)fOutputCheckCentrality->FindObject("hV0MultiplicityPercentile"))->Fill(vzeroMult,stdCentf);
      ((TH2F*)fOutputCheckCentrality->FindObject("hV0MultiplicityNtrackletsIn"))->Fill(vzeroMult,aod->GetTracklets()->GetNumberOfTracklets());
      ((TH2F*)fOutputCheckCentrality->FindObject("hStdPercentileSPDPercentile"))->Fill(stdCentf,spdCentf);

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
    Double_t xvtx=vertex->GetX();
    Double_t yvtx=vertex->GetY();
    Double_t zvtx=vertex->GetZ();
    Int_t vtxTyp=0;
    if(vertex->GetNContributors()<=0) vtxTyp=-1;
    TString title=vertex->GetTitle();
    if(title.Contains("Z")) vtxTyp=3;
    if(title.Contains("3D")) vtxTyp=2;    
    ((TH1F*)fOutputEvSelection->FindObject("hxvtx"))->Fill(xvtx);
    ((TH1F*)fOutputEvSelection->FindObject("hyvtx"))->Fill(yvtx);
    ((TH1F*)fOutputEvSelection->FindObject("hzvtx"))->Fill(zvtx);
    ((TH1F*)fOutputEvSelection->FindObject("hWhichVert"))->Fill(vtxTyp);
    if(evSelected){
      ((TH1F*)fOutputEvSelection->FindObject("hxvtxSelEv"))->Fill(xvtx);
      ((TH1F*)fOutputEvSelection->FindObject("hyvtxSelEv"))->Fill(yvtx);
      ((TH1F*)fOutputEvSelection->FindObject("hzvtxSelEv"))->Fill(zvtx);
      ((TH1F*)fOutputEvSelection->FindObject("hWhichVertSelEv"))->Fill(vtxTyp);
    }
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
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse *pidResp=inputHandler->GetPIDResponse();

  //AliPIDResponse* respF=pidHF->GetPidResponse();
  AliTPCPIDResponse* tpcres=new AliTPCPIDResponse();
  Bool_t oldPID=pidHF->GetOldPid();
  if(oldPID){ 
    Double_t alephParameters[5];
    pidHF->GetTPCBetheBlochParams(alephParameters);
    tpcres->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);
  }


  Int_t ntracks=0;
  Int_t isGoodTrack=0, isFakeTrack=0, isSelTrack=0;

  if(aod) ntracks=aod->GetNTracks();

  if(fOnOff[0] || fOnOff[1]){
    //loop on tracks in the event
    for (Int_t k=0;k<ntracks;k++){
      AliAODTrack* track=aod->GetTrack(k);

      // Track selection cuts
      if(track->GetID()<0) continue;
      Bool_t selTrack=kTRUE;
      ULong_t trStatus=track->GetStatus();
      if (!((trStatus & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
	  !((trStatus & AliVTrack::kITSrefit) == AliVTrack::kITSrefit)){
	selTrack=kFALSE;
      }
      Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
      Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
      if (track->GetTPCNclsF()>0) {
	ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
      }
      if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ){
	selTrack=kFALSE;	
      }

      AliAODPid *pid = track->GetDetPid();
      if(!pid && fDebug>1) cout<<"No AliAODPid found"<<endl;

      if(pid && fOnOff[1]){
	Double_t times[AliPID::kSPECIES];
	pid->GetIntegratedTimes(times);
    
	Double_t tofRes[AliPID::kSPECIES];
	pid->GetTOFpidResolution(tofRes);

	//check TOF
	TH1F* htmpfl=((TH1F*)fOutputPID->FindObject("hTOFflags"));
	htmpfl->Fill(0.);
	if (trStatus&AliESDtrack::kTPCout) htmpfl->Fill(1.);
	if (trStatus&AliESDtrack::kTOFout) htmpfl->Fill(2.);
	if (trStatus&AliESDtrack::kTIME) htmpfl->Fill(3.);
	if (trStatus&AliESDtrack::kTOFpid) htmpfl->Fill(4.);
	if (trStatus&AliESDtrack::kTOFmismatch) htmpfl->Fill(5.);

	Bool_t isTOFok=kFALSE;
	if(pidResp){
	  Double_t prob[AliPID::kSPECIES];
	  if(pidResp->ComputeTOFProbability(track,AliPID::kSPECIES,prob)==AliPIDResponse::kDetPidOk){
	    isTOFok=kTRUE;
	    htmpfl->Fill(6.);
	  }
	}

	if(selTrack && isTOFok){
	  Double_t tofTime=pid->GetTOFsignal();
	  AliTOFHeader* tofH=(AliTOFHeader*)aod->GetTOFHeader();
	  if (tofH && (TMath::Abs(tofRes[0]) <= 1.E-16) ) { // new AOD
            // with new AOD we need to retrieve startTime, subtract it and retrieve correctly TOF PID resolutions  *PA*
	    AliTOFPIDResponse tofResp=pidResp->GetTOFResponse();
	    Double_t startTime = tofResp.GetStartTime(track->P());
	    Float_t startTimeRes = tofResp.GetStartTimeRes(track->P());  
	    Int_t startTimeMask = tofResp.GetStartTimeMask(track->P());  
	    ((TH1F*)fOutputPID->FindObject("hTOFstartTimeDistrib"))->Fill(startTime);
	    ((TH1F*)fOutputPID->FindObject("hTOFstartTimeMask"))->Fill(startTimeMask);
	    ((TH1F*)fOutputPID->FindObject("hTOFstartTimeRes"))->Fill(startTimeRes);
	    tofTime-=startTime;
	    for (Int_t type=0;type<AliPID::kSPECIES;type++) tofRes[type]=tofResp.GetExpectedSigma(track->P(),times[type],AliPID::ParticleMassZ(type)); 
	  }
	  ((TH1F*)fOutputPID->FindObject("hTOFtime"))->Fill(times[AliPID::kProton]);
	  ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptime"))->Fill(track->P(),tofTime-times[3]); //3 is kaon
	  ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(tofTime);
	  if (pid->GetTOFsignal()< 0) ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(-1);

	  Double_t nsigma[3]={-10,-10,-10};
	  nsigma[0]=pidResp->NumberOfSigmasTOF(track,AliPID::kPion);
	  nsigma[1]=pidResp->NumberOfSigmasTOF(track,AliPID::kKaon);
	  nsigma[2]=pidResp->NumberOfSigmasTOF(track,AliPID::kProton);

	  ((TH2F*)fOutputPID->FindObject("hTOFsigmaKSigPid"))->Fill(track->P(),nsigma[1]);
	  ((TH2F*)fOutputPID->FindObject("hTOFsigmaPionSigPid"))->Fill(track->P(),nsigma[0]);
	  ((TH2F*)fOutputPID->FindObject("hTOFsigmaProtonSigPid"))->Fill(track->P(),nsigma[2]);
	  if(fReadMC){
	    Int_t label=track->GetLabel();
	    if(label<=0) continue;
	    AliMCParticle* mcpart=(AliMCParticle*)mcArray->At(label);
	    if(mcpart){
	      Int_t abspdgcode=TMath::Abs(mcpart->PdgCode());
	      if(abspdgcode==211) ((TH2F*)fOutputPID->FindObject("hTOFsigmaMCPionSigPid"))->Fill(track->P(),nsigma[0]);
	      if(abspdgcode==321) ((TH2F*)fOutputPID->FindObject("hTOFsigmaMCKSigPid"))->Fill(track->P(),nsigma[1]);
	      if(abspdgcode==2212) ((TH2F*)fOutputPID->FindObject("hTOFsigmaMCProtonSigPid"))->Fill(track->P(),nsigma[2]);

	    }
	  }

	  for (Int_t iS=2; iS<5; iS++){ //we plot TOF Pid resolution for 3-sigma identified particles
	    if ( TMath::Abs(nsigma[iS-2])<3.){
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
	}//if TOF status
	//}

	if(pidHF && pidHF->CheckStatus(track,"TPC")){ 

	  Double_t TPCp=pid->GetTPCmomentum();
	  Double_t TPCsignal=pid->GetTPCsignal();
	  ((TH1F*)fOutputPID->FindObject("hTPCsig"))->Fill(TPCsignal);
	  ((TH1F*)fOutputPID->FindObject("hTPCsigvsp"))->Fill(TPCp,TPCsignal);
	  //if (pidHF->IsKaonRaw(track, "TOF"))
	  Double_t nsigma[3]={-10,-10,-10};
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kPion,nsigma[0]);	 
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kKaon,nsigma[1]);	 
	  pidHF->GetnSigmaTPC(track,(Int_t)AliPID::kProton,nsigma[2]);	 

	  ((TH2F*)fOutputPID->FindObject("hTPCsigmaK"))->Fill(TPCp,nsigma[1]);
	  
	  ((TH2F*)fOutputPID->FindObject("hTPCsigmaPion"))->Fill(TPCp,nsigma[0]);	  
	  ((TH2F*)fOutputPID->FindObject("hTPCsigmaProton"))->Fill(TPCp,nsigma[2]);
	  
	  if(fReadMC){
	    Int_t label=track->GetLabel();
	    if(label<=0) continue;
	    AliMCParticle* mcpart=(AliMCParticle*)mcArray->At(label);
	    if(mcpart){
	      Int_t abspdgcode=TMath::Abs(mcpart->PdgCode());
	      if(abspdgcode==211) ((TH2F*)fOutputPID->FindObject("hTPCsigmaMCPion"))->Fill(track->P(),nsigma[0]);
	      if(abspdgcode==321) ((TH2F*)fOutputPID->FindObject("hTPCsigmaMCK"))->Fill(track->P(),nsigma[1]);
	      if(abspdgcode==2212) ((TH2F*)fOutputPID->FindObject("hTPCsigmaMCProton"))->Fill(track->P(),nsigma[2]);

	    }

	  }


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
	if(fCuts->IsDaughterSelected(track,&vESD,fCuts->GetTrackCuts())){
	  isSelTrack++;
	}//select tracks for our analyses
      } //fill track histos
    } //end loop on tracks

    //fill once per event
    if(fOnOff[0]){
      if (fReadMC) ((TH1F*)fOutputTrack->FindObject("hdistrFakeTr"))->Fill(isFakeTrack);
      ((TH1F*)fOutputTrack->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
      ((TH1F*)fOutputTrack->FindObject("hdistrSelTr"))->Fill(isSelTrack);
    }

    if(!isSimpleMode){
      // loop over candidates
      Int_t nCand = arrayProng->GetEntriesFast();
      Int_t ndaugh=3;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi) ndaugh=2;
      if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpipipi) ndaugh=4;

      for (Int_t iCand = 0; iCand < nCand; iCand++) {
	AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
	if(fUseSelectionBit && d->GetSelectionMap()) {
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi && !d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts)) continue; //skip the D0 from Dstar
	  if(fDecayChannel==AliAnalysisTaskSEHFQA::kDplustoKpipi && !d->HasSelectionBit(AliRDHFCuts::kDplusCuts)) continue; //skip the 3 prong !D+
	}

	if(fReadMC){ 

	  Int_t labD = -1;
	  if (fDecayChannel==AliAnalysisTaskSEHFQA::kLambdactoV0 && (dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0()) {

	    Int_t pdgDgLctoV0bachelor[2]={310,2212};
	    Int_t pdgDgV0toDaughters[2]={211,211};
	    Int_t mcLabelK0S = (dynamic_cast<AliAODRecoCascadeHF*>(d))->MatchToMC(pdg,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE); // Lc->K0S+p and cc
	    pdgDgLctoV0bachelor[0]=3122, pdgDgLctoV0bachelor[1]=211;
	    pdgDgV0toDaughters[0]=2212,  pdgDgV0toDaughters[1]=211;
	    Int_t mcLabelLambda = (dynamic_cast<AliAODRecoCascadeHF*>(d))->MatchToMC(pdg,pdgDgLctoV0bachelor[0],pdgDgLctoV0bachelor,pdgDgV0toDaughters,mcArray,kTRUE); // Lc->Lambda+pi and cc
	    if (mcLabelK0S!=-1 || mcLabelLambda!=-1) AliInfo(Form("mcLabelK0S=%d - mcLabelLambda=%d",mcLabelK0S,mcLabelLambda));

	    if (mcLabelK0S!=-1 && mcLabelLambda!=-1)
	      AliInfo("Strange: current Lc->V0+bachelor candidate has two MC different labels!");
	    else if (mcLabelK0S>-1 && mcLabelLambda==-1)
	      labD = mcLabelK0S;
	    else if (mcLabelLambda>-1 && mcLabelK0S==-1)
	      labD = mcLabelLambda;
	  }
	  else
	    labD = d->MatchToMC(pdg,mcArray,ndaugh,pdgdaughters);

	  if(labD>=0){
	    AliAODMCParticle *partD = (AliAODMCParticle*)mcArray->At(labD);
	    Int_t label=partD->GetMother();
	    AliAODMCParticle *mot = (AliAODMCParticle*)mcArray->At(label);
	    while(label>=0){//get first mother
	      mot = (AliAODMCParticle*)mcArray->At(label);
	      label=mot->GetMother();
	    }
	    if(mot){
	      Int_t pdgMotCode = mot->GetPdgCode();
	
	      if(TMath::Abs(pdgMotCode)==4) fNEntries->Fill(6); //from primary charm
	      if(TMath::Abs(pdgMotCode)==5) fNEntries->Fill(7); //from beauty
	    }
	  }
	}//end MC
	fNEntries->Fill(5); //count the candidates (data and MC)

	for(Int_t id=0;id<ndaugh;id++){
	  //other histograms to be filled when the cut object is given
	  AliAODTrack* track=0;

	  if (fDecayChannel==AliAnalysisTaskSEHFQA::kLambdactoV0 && (dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0()) {
	    if (id==0)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->GetBachelor();
	    else if (id==1)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0PositiveTrack();
	    else if (id==2)
	      track=(AliAODTrack*)(dynamic_cast<AliAODRecoCascadeHF*>(d))->Getv0NegativeTrack();
	  }
	  else 
	    track=(AliAODTrack*)d->GetDaughter(id);

	  //track quality

	  if (fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(pdg)) && fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod)) {
	    
	    Int_t label=0;
	    if(fReadMC)label=track->GetLabel();
	    if(fOnOff[0]){
	      
	      if(fReadMC && label<0) {
		isFakeTrack++;
		((TH1F*)fOutputTrack->FindObject("hptFakeTrFromDaugh"))->Fill(track->Pt());
	   
		((TH1F*)fOutputTrack->FindObject("hd0f"))->Fill(d->Getd0Prong(id));
	      } else {
		((TH1F*)fOutputTrack->FindObject("hptGoodTrFromDaugh"))->Fill(track->Pt());
		((TH1F*)fOutputTrack->FindObject("hd0"))->Fill(d->Getd0Prong(id));
		Double_t d0rphiz[2],covd0[3];
		Bool_t isDCA=track->PropagateToDCA(aod->GetPrimaryVertex(),aod->GetMagneticField(),9999.,d0rphiz,covd0);
		if(isDCA) ((TH1F*)fOutputTrack->FindObject("hd0z"))->Fill(d0rphiz[1]);
	      }
	    }
	    if (fCuts->IsSelected(d,AliRDHFCuts::kAll,aod) && fOnOff[1]){
	      fNEntries->Fill(3); //candidates passing analysis cuts
	      AliAODPid *pid = track->GetDetPid();
	      if(pid){
		Double_t times[5];
		pid->GetIntegratedTimes(times);
		if(pidHF && pidHF->CheckStatus(track,"TOF")){
		  Double_t tofTime=pid->GetTOFsignal();
		  AliTOFHeader* tofH=(AliTOFHeader*)aod->GetTOFHeader();
		  Double_t tofRes[AliPID::kSPECIES];
		  pid->GetTOFpidResolution(tofRes);
		  if (tofH && (TMath::Abs(tofRes[0]) <= 1.E-16) ) { // new AOD
		    AliTOFPIDResponse tofResp=pidHF->GetPidResponse()->GetTOFResponse();
		    Double_t startTime=tofResp.GetStartTime(track->P());
		    tofTime-=startTime;
		  }
		  ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptimeAC"))->Fill(track->P(),tofTime-times[AliPID::kKaon]);
		}
		if(pidHF && pidHF->CheckStatus(track,"TPC")) ((TH2F*)fOutputPID->FindObject("hTPCsigvspAC"))->Fill(pid->GetTPCmomentum(),pid->GetTPCsignal());
	      }
	      
	    } //end analysis cuts
	  } //end acceptance and track cuts
	} //end loop on tracks in the candidate
      } //end loop on candidates
      
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
void AliAnalysisTaskSEHFQA::FillFlowObs(AliAODEvent *aod){
  //fills the flow observables
  Double_t cc;
  cc = fCuts->GetCentrality(aod);
  ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(0., cc);

  UInt_t mask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  UInt_t trigger=AliVEvent::kMB | AliVEvent::kCentral | AliVEvent::kSemiCentral;
  if(mask & trigger) {
    ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(1.,cc); // fired
    if (mask & AliVEvent::kMB) ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(2.,cc);
    if (mask & AliVEvent::kCentral) ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(3.,cc);
    if (mask & AliVEvent::kSemiCentral) ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(4.,cc);
    Bool_t rejected=false;
    if(cc<0 || cc>60) rejected=true;
    const AliVVertex *vertex = aod->GetPrimaryVertex();
    Double_t zvtx=vertex->GetZ();
    if(TMath::Abs(zvtx)>fCuts->GetMaxVtxZ()) rejected=true;
    if(rejected) return; //not interesting for flow QA
  } else {
    return;
  }

  // event accepted
  ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(5.,cc);
  fRFPcuts->SetParamType(AliFlowTrackCuts::kGlobal);
  fRFPcuts->SetPtRange(0.2,5.);
  fRFPcuts->SetEtaRange(-0.8,0.8);
  fRFPcuts->SetMinNClustersTPC(70);
  fRFPcuts->SetMinChi2PerClusterTPC(0.2);
  fRFPcuts->SetMaxChi2PerClusterTPC(4.0);
  fRFPcuts->SetAcceptKinkDaughters(kFALSE);
  fRFPcuts->SetEvent(aod);

  TString ref[3] = {"FB1","FB128","VZE"};
  Double_t psi[3];
  for(Int_t i=0; i!=3; ++i) {
    if(i==0) { // switching to bit 1
      fRFPcuts->SetMinimalTPCdedx(10.);
      fRFPcuts->SetAODfilterBit(1);
    } else { // switching to bit 128
      fRFPcuts->SetMinimalTPCdedx(-1);
      fRFPcuts->SetAODfilterBit(128);
    }
    if(i>1) {
      fRFPcuts->SetParamType(AliFlowTrackCuts::kV0);
      fRFPcuts->SetEtaRange(-5,+5);
      fRFPcuts->SetPhiMin(0);
      fRFPcuts->SetPhiMax(TMath::TwoPi());
    }
    fFlowEvent->Fill(fRFPcuts,fRFPcuts);
    fFlowEvent->TagSubeventsInEta(-5,0,0,+5);
    // getting informationt
    AliFlowVector vQ, vQaQb[2];
    fFlowEvent->Get2Qsub(vQaQb,2);
    vQ = vQaQb[0]+vQaQb[1];
    Double_t dMa=vQaQb[0].GetMult();
    Double_t dMb=vQaQb[1].GetMult();
    if( dMa<2 || dMb<2 ) {
      ((TH2F*) fOutputFlowObs->FindObject("hFlowEvents"))->Fill(6.,cc); //???
      continue;
    }
    psi[i] = vQ.Phi()/2;
    // publishing
    ((TProfile2D*) fOutputFlowObs->FindObject( Form("h%s_Q",ref[i].Data())))->Fill(0,cc,vQaQb[0].X()/dMa,dMa); // Qx-
    ((TProfile2D*) fOutputFlowObs->FindObject( Form("h%s_Q",ref[i].Data())))->Fill(1,cc,vQaQb[0].Y()/dMa,dMa); // Qy-
    ((TProfile2D*) fOutputFlowObs->FindObject( Form("h%s_Q",ref[i].Data())))->Fill(2,cc,vQaQb[1].X()/dMb,dMb); // Qx+
    ((TProfile2D*) fOutputFlowObs->FindObject( Form("h%s_Q",ref[i].Data())))->Fill(3,cc,vQaQb[1].Y()/dMb,dMb); // Qy+
    ((TH2F*) fOutputFlowObs->FindObject( Form("h%s_AngleQ",ref[i].Data()) ))->Fill(psi[i],cc); // Psi
    AliFlowTrackSimple *track;
    for(Int_t t=0; t!=fFlowEvent->NumberOfTracks(); ++t) {
      track = (AliFlowTrackSimple*) fFlowEvent->GetTrack(t);
      if(!track) continue;
      if(!track->InRPSelection()) continue;
      ((TH3F*) fOutputFlowObs->FindObject( Form("h%s_PhiEta",ref[i].Data()) ))->Fill(track->Phi(),track->Eta(),cc,track->Weight()); //PhiEta
    }
  
  //histo filled only for TPCFB1
  if (i==0) {
    ((TH2F*) fOutputFlowObs->FindObject("hCentVsMultRPS"))->Fill(fFlowEvent->GetNumberOfRPs(),cc);
  }
  }
  // TPC vs VZERO
  ((TH3F*) fOutputFlowObs->FindObject( "hTPCVZE_AngleQ" ))->Fill(psi[0],psi[2],cc);
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

