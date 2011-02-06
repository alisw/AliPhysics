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
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAnalysisVertexingHF.h"
#include "AliAnalysisTaskSE.h"
#include "AliRDHFCuts.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"

#include "AliAnalysisTaskSEHFQA.h"


ClassImp(AliAnalysisTaskSEHFQA)

//____________________________________________________________________________

AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA():AliAnalysisTaskSE(),
  fNEntries(0x0),
  fOutputPID(0x0),
  fOutputTrack(0x0),
  fDecayChannel(AliAnalysisTaskSEHFQA::kD0toKpi),
  fCuts(0x0)
{
  //default constructor
}

//____________________________________________________________________________
AliAnalysisTaskSEHFQA::AliAnalysisTaskSEHFQA(const char *name, AliAnalysisTaskSEHFQA::DecChannel ch,AliRDHFCuts* cuts):
AliAnalysisTaskSE(name),
fNEntries(0x0),
fOutputPID(0x0),
fOutputTrack(0x0),
fDecayChannel(ch),
fCuts(0x0)
{
  //constructor

  //SetCutObject(cuts);
  fCuts=cuts;

  // Output slot #1 writes into a TH1F container (number of events)
  DefineOutput(1,TH1F::Class());  //My private output
  // Output slot #2 writes into a TList container (PID)
  DefineOutput(2,TList::Class());  //My private output
  // Output slot #3 writes into a TList container (Tracks)
  DefineOutput(3,TList::Class());  //My private output
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
}

//___________________________________________________________________________
AliAnalysisTaskSEHFQA::~AliAnalysisTaskSEHFQA()
{
  //destructor
  if(fNEntries){
    delete fNEntries;
    fNEntries=0;
  }
  if(fOutputPID){
    delete fOutputPID;
    fOutputPID=0;
  }
  if(fOutputTrack){
    delete fOutputTrack;
    fOutputTrack=0;
  }
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

  fNEntries=new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Counts the number of events", 6,-0.5,5.5);
  fNEntries->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fNEntries->GetXaxis()->SetBinLabel(2,"Pile-up Rej");
  fNEntries->GetXaxis()->SetBinLabel(3,"No VertexingHF");
  fNEntries->GetXaxis()->SetBinLabel(4,"nCandidates(AnCuts)");
  fNEntries->GetXaxis()->SetBinLabel(5,"EventsWithGoodVtx");
  fNEntries->GetXaxis()->SetBinLabel(6,"N. of 0SMH");
  fNEntries->GetXaxis()->SetNdivisions(1,kFALSE);

  //PID

  fOutputPID=new TList();
  fOutputPID->SetOwner();
  fOutputPID->SetName(GetOutputSlot(2)->GetContainer()->GetName());

  //TOF pid
  TString hname="hTOFsig";
  TH1F* hTOFsig=new TH1F(hname.Data(),"Distribution of TOF signal;TOF time [ps];Entries", 100, -2.e3,40.e3);

  hname="hTOFtime";
  TH1F* hTOFtime=new TH1F(hname.Data(),"Distribution of TOF time Kaon;TOF time(Kaon) [ps];Entries", 1000, 0.,50000.);

  hname="hTOFtimeKaonHyptime";
  TH2F* hTOFtimeKaonHyptime=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",200,0.,4.,1000,-20000.,20000.);

  hname="hTOFtimeKaonHyptimeAC";
  TH2F* hTOFtimeKaonHyptimeAC=new TH2F(hname.Data(),"TOFtime - timeHypothesisForKaon;p[GeV/c];TOFtime - timeHypothesisForKaon [ps]",200,0.,4.,1000,-20000.,20000.);

  hname="hTOFsigmaK160";
  TH2F* hTOFsigmaK160=new TH2F(hname.Data(),"(TOFsignal-timeK)/160 ps;p[GeV/c];(TOFsignal-timeK)/160 ps",200,0.,4.,100,-5,5);

  hname="hTOFsigmaPion160";
  TH2F* hTOFsigmaPion160=new TH2F(hname.Data(),"(TOFsignal-time#pi)/160 ps;p[GeV/c];(TOFsignal-time#pi)/160 ps",200,0.,4.,100,-5,5);

  hname="hTOFsigmaProton160";
  TH2F* hTOFsigmaProton160=new TH2F(hname.Data(),"(TOFsignal-timep)/160 ps;p[GeV/c];(TOFsignal-time p)/160 ps",200,0.,4.,100,-5,5);

  hname="hTOFsigmaK120";
  TH2F* hTOFsigmaK120=new TH2F(hname.Data(),"(TOFsignal-timeK)/120 ps;p[GeV/c];(TOFsignal-timeK)/120 ps",200,0.,4.,100,-5,5);

  hname="hTOFsigmaPion120";
  TH2F* hTOFsigmaPion120=new TH2F(hname.Data(),"(TOFsignal-time#pi)/120 ps;p[GeV/c];(TOFsignal-time#pi)/120 ps",200,0.,4.,100,-5,5);

  hname="hTOFsigmaProton120";
  TH2F* hTOFsigmaProton120=new TH2F(hname.Data(),"(TOFsignal-timep)/120 ps;p[GeV/c];(TOFsignal-time p)/120 ps",200,0.,4.,100,-5,5);

  //TPC pid
  hname="hTPCsig";
  TH1F* hTPCsig=new TH1F(hname.Data(),"Distribution of TPC signal;TPC sig;Entries", 100, 35.,100.);

  hname="hTPCsigvsp";
  TH2F* hTPCsigvsp=new TH2F(hname.Data(),"TPCsig vs p;TPC p[GeV/c];TPCsig",200,0.,4.,1000,35.,100.);
 
  hname="hTPCsigvspAC";
  TH2F* hTPCsigvspAC=new TH2F(hname.Data(),"TPCsig vs p;TPCp[GeV/c];TPCsig",200,0.,4.,1000,35.,100.);

  hname="hTPCsigmaK";
  TH2F* hTPCsigmaK=new TH2F(hname.Data(),"TPC Sigma for K as a function of momentum;p[GeV/c];Sigma Kaon",200,0.,4.,200,-5,5);

  hname="hTPCsigmaPion";
  TH2F* hTPCsigmaPion=new TH2F(hname.Data(),"TPC Sigma for #pi as a function of momentum;p[GeV/c];Sigma #pi",200,0.,4.,200,-5,5);

  hname="hTPCsigmaProton";
  TH2F* hTPCsigmaProton=new TH2F(hname.Data(),"TPC Sigma for proton as a function of momentum;p[GeV/c];Sigma Proton",200,0.,4.,200,-5,5);

  fOutputPID->Add(hTOFsig);
  fOutputPID->Add(hTPCsig);
  fOutputPID->Add(hTOFtime);
  fOutputPID->Add(hTOFtimeKaonHyptime);
  fOutputPID->Add(hTOFtimeKaonHyptimeAC);
  fOutputPID->Add(hTOFsigmaK160);
  fOutputPID->Add(hTOFsigmaPion160);
  fOutputPID->Add(hTOFsigmaProton160);
  fOutputPID->Add(hTOFsigmaK120);
  fOutputPID->Add(hTOFsigmaPion120);
  fOutputPID->Add(hTOFsigmaProton120);
  fOutputPID->Add(hTPCsigvsp);
  fOutputPID->Add(hTPCsigvspAC);
  fOutputPID->Add(hTPCsigmaK);
  fOutputPID->Add(hTPCsigmaPion);
  fOutputPID->Add(hTPCsigmaProton);

  //quality of the tracks

  fOutputTrack=new TList();
  fOutputTrack->SetOwner();
  fOutputTrack->SetName(GetOutputSlot(3)->GetContainer()->GetName());

  hname="hnClsITS";
  TH1F* hnClsITS=new TH1F(hname.Data(),"Distribution of number of ITS clusters;nITScls;Entries",7,-0.5,6.5);

  hname="hnClsITS-SA";
  TH1F* hnClsITSSA=new TH1F(hname.Data(),"Distribution of number of ITS clusters(ITS-SA);nITScls;Entries",7,-0.5,6.5);

  hname="hnClsSPD";
  TH1F* hnClsSPD=new TH1F(hname.Data(),"Distribution of number of SPD clusters;nSPDcls;Entries",3,-0.5,2.5);

  hname="hptGoodTr";
  TH1F* hptGoodTr=new TH1F(hname.Data(),"Pt distribution of 'good' tracks;p_{t}[GeV];Entries/0.05 GeV/c",400,0.,20.);
  hptGoodTr->SetTitleOffset(1.3,"Y");

  hname="hdistrGoodTr";
  TH1F* hdistrGoodTr=new TH1F(hname.Data(),"Distribution of number of 'good' tracks per event;no.good-tracks/ev;Entries",4000,-0.5,3999.5);
  hdistrGoodTr->SetTitleOffset(1.3,"Y");

  hname="hNtracklets";
  TH1F* hNtracklets=new TH1F(hname.Data(),"Number of tracklets;ntracklets;Entries",5000,-0.5,4999.5);

  hname="hMult";
  TH1F* hMult=new TH1F(hname.Data(),"Multiplicity;multiplicity;Entries",10000,-0.5,9999.5);

  hname="hd0";
  TH1F* hd0=new TH1F(hname.Data(),"Impact parameter distribution of 'good' tracks;d_{0}[cm];Entries/10^{3} cm",200,-0.1,0.1);

  fOutputTrack->Add(hnClsITS);
  fOutputTrack->Add(hnClsITSSA);
  fOutputTrack->Add(hnClsSPD);
  fOutputTrack->Add(hptGoodTr);
  fOutputTrack->Add(hdistrGoodTr);
  fOutputTrack->Add(hNtracklets);
  fOutputTrack->Add(hMult);
  fOutputTrack->Add(hd0);
  
  // Post the data
  PostData(1,fNEntries);
  PostData(2,fOutputPID);
  PostData(3,fOutputTrack);
  PostData(4,fCuts);
}

//___________________________________________________________________________
void AliAnalysisTaskSEHFQA::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event

  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(fDebug>2) printf("Analysing decay %d\n",fDecayChannel);
  // Post the data already here
  PostData(1,fNEntries);
  PostData(2,fOutputPID);
  PostData(3,fOutputTrack);
  PostData(4,fCuts);

  TClonesArray *arrayProng =0;
  Int_t pdg=0;
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
	break; 
      case 1:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
	pdg=421;
	break; 
      case 2:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
	pdg=413;
	break; 
      case 3:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	pdg=431;
	break; 
      case 4:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm4Prong");
	pdg=421;
	break; 
      case 5:
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
	pdg=4122;
	break; 
      }
    }
  } else if(aod) {
    switch(fDecayChannel){
    case 0:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=411;
      break; 
    case 1:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
      pdg=421;
      break; 
    case 2:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
      pdg=413;
      break; 
    case 3:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=431;
      break; 
    case 4:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm4Prong");
      pdg=421;
      break; 
    case 5:
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
      pdg=4122;
      break; 
    }
  }
  Bool_t isSimpleMode=kFALSE;
  if(!arrayProng) {
    AliInfo("Branch not found! The output will contain only trak related histograms\n");
    isSimpleMode=kTRUE;
    fNEntries->Fill(2);
  }
  
  if(!aod) return;
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  // count event
  fNEntries->Fill(0); 

  //count events with good vertex
  // AOD primary vertex
  AliAODVertex *vtx1 = (AliAODVertex*)aod->GetPrimaryVertex();
  TString primTitle = vtx1->GetTitle();
  if(primTitle.Contains("VertexerTracks") && vtx1->GetNContributors()>0) fNEntries->Fill(4);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fNEntries->Fill(5);


  //select event
  if(!fCuts->IsEventSelected(aod)) {
    // rejected for pileup
    if(fCuts->GetWhyRejection()==1) fNEntries->Fill(1);
    return;
  }

  Int_t ntracks=0;
  Int_t isGoodTrack=0;

  if(aod) ntracks=aod->GetNTracks();


  ((TH1F*)fOutputTrack->FindObject("hNtracklets"))->Fill(aod->GetTracklets()->GetNumberOfTracklets());
  ((TH1F*)fOutputTrack->FindObject("hMult"))->Fill(aod->GetHeader()->GetRefMultiplicity());


  //loop on tracks in the event
  for (Int_t k=0;k<ntracks;k++){
    AliAODTrack* track=aod->GetTrack(k);
    AliAODPidHF* pidHF=fCuts->GetPidHF();
    AliAODPid *pid = track->GetDetPid();
    if(!pid)  {if (fDebug>1)cout<<"No AliAODPid found"<<endl; continue;}

    Double_t times[AliPID::kSPECIES];
    pid->GetIntegratedTimes(times);
    
    //check TOF
    if(pidHF && pidHF->CheckStatus(track,"TOF")){
      ((TH1F*)fOutputPID->FindObject("hTOFtime"))->Fill(times[AliPID::kProton]);
      ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptime"))->Fill(track->P(),pid->GetTOFsignal()-times[3]); //3 is kaon
      ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(pid->GetTOFsignal());
      if (pid->GetTOFsignal()< 0) ((TH1F*)fOutputPID->FindObject("hTOFsig"))->Fill(-1);

      ((TH2F*)fOutputPID->FindObject("hTOFsigmaK160"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kKaon])/160);
      ((TH2F*)fOutputPID->FindObject("hTOFsigmaPion160"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kPion])/160);
      ((TH2F*)fOutputPID->FindObject("hTOFsigmaProton160"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kProton])/160);
      ((TH2F*)fOutputPID->FindObject("hTOFsigmaK120"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kKaon])/120);
      ((TH2F*)fOutputPID->FindObject("hTOFsigmaPion120"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kPion])/120);
      ((TH2F*)fOutputPID->FindObject("hTOFsigmaProton120"))->Fill(track->P(),(pid->GetTOFsignal()-times[AliPID::kProton])/120);

    }//if TOF status

    if(pidHF && pidHF->CheckStatus(track,"TPC")){ 
      Double_t alephParameters[5];
      //this is recommended for LHC10d
      alephParameters[0] = 1.34490e+00/50.;
      alephParameters[1] =  2.69455e+01;
      alephParameters[2] =  TMath::Exp(-2.97552e+01);
      alephParameters[3] = 2.35339e+00;
      alephParameters[4] = 5.98079e+00;
      /*
      //this is recommended for LHC10bc
      alephParameters[0] = 0.0283086/0.97;
      alephParameters[1] = 2.63394e+01;
      alephParameters[2] = 5.04114e-11;
      alephParameters[3] = 2.12543e+00;
      alephParameters[4] = 4.88663e+00;
      */
      AliTPCPIDResponse* tpcres=new AliTPCPIDResponse();
      tpcres->SetBetheBlochParameters(alephParameters[0],alephParameters[1],alephParameters[2],alephParameters[3],alephParameters[4]);
      Double_t TPCp=pid->GetTPCmomentum();
      Double_t TPCsignal=pid->GetTPCsignal();
      ((TH1F*)fOutputPID->FindObject("hTPCsig"))->Fill(TPCsignal);
      ((TH1F*)fOutputPID->FindObject("hTPCsigvsp"))->Fill(TPCp,TPCsignal);
      ((TH2F*)fOutputPID->FindObject("hTPCsigmaK"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kKaon));
      ((TH2F*)fOutputPID->FindObject("hTPCsigmaPion"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kPion));
      ((TH2F*)fOutputPID->FindObject("hTPCsigmaProton"))->Fill(TPCp,tpcres->GetNumberOfSigmas(TPCp,TPCsignal,track->GetTPCNcls(),AliPID::kProton));
      delete tpcres;
    }//if TPC status

    //check clusters of the tracks
    Int_t nclsTot=0,nclsSPD=0;
    
    for(Int_t l=0;l<6;l++) {
      if(TESTBIT(track->GetITSClusterMap(),l)) {
	nclsTot++; if(l<2) nclsSPD++;
      }
    }
    ((TH1F*)fOutputTrack->FindObject("hnClsITS"))->Fill(nclsTot);
    ((TH1F*)fOutputTrack->FindObject("hnClsSPD"))->Fill(nclsSPD);

    if(!(track->GetStatus()&AliESDtrack::kTPCin) && track->GetStatus()&AliESDtrack::kITSrefit && !(track->GetStatus()&AliESDtrack::kITSpureSA)){//tracks retrieved in the ITS and not reconstructed in the TPC
      ((TH1F*)fOutputTrack->FindObject("hnClsITS-SA"))->Fill(nclsTot);
    }

    if(isSimpleMode){

      if (track->Pt()>0.3 &&
	  track->GetStatus()&AliESDtrack::kTPCrefit &&
	  track->GetStatus()&AliESDtrack::kITSrefit &&
	  /*nclsTot>3 &&*/
	  nclsSPD>0) {//fill hist good tracks

	((TH1F*)fOutputTrack->FindObject("hptGoodTr"))->Fill(track->Pt());
	
	isGoodTrack++;
      
	((TH1F*)fOutputTrack->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
      }
    }//simple mode: no IsSelected on tracks: use "manual" cuts
      
  } //end loop on tracks

  if(!isSimpleMode){
  // loop over candidates
  Int_t nCand = arrayProng->GetEntriesFast();
  Int_t ndaugh=3;
  if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpi) ndaugh=2;
  if(fDecayChannel==AliAnalysisTaskSEHFQA::kD0toKpipipi) ndaugh=4;

  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    AliAODRecoDecayHF *d = (AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);

    for(Int_t id=0;id<ndaugh;id++){

      //other histograms to be filled when the cut object is given
      AliAODTrack* track=(AliAODTrack*)d->GetDaughter(id);

      //track quality

      if (fCuts->IsInFiducialAcceptance(d->Pt(),d->Y(pdg)) && fCuts->IsSelected(d,AliRDHFCuts::kTracks,aod)) {
	
	((TH1F*)fOutputTrack->FindObject("hptGoodTr"))->Fill(track->Pt());
	isGoodTrack++;
      
	((TH1F*)fOutputTrack->FindObject("hdistrGoodTr"))->Fill(isGoodTrack);
      
	((TH1F*)fOutputTrack->FindObject("hd0"))->Fill(d->Getd0Prong(id));
  
	if (fCuts->IsSelected(d,AliRDHFCuts::kAll,aod)){
	  
	  AliAODPid *pid = track->GetDetPid();
	  AliAODPidHF* pidHF=fCuts->GetPidHF();
	  Double_t times[5];
	  pid->GetIntegratedTimes(times);
	  if(pidHF && pidHF->CheckStatus(track,"TOF")) ((TH2F*)fOutputPID->FindObject("hTOFtimeKaonHyptimeAC"))->Fill(track->P(),pid->GetTOFsignal()-times[AliPID::kKaon]);
	  if(pidHF && pidHF->CheckStatus(track,"TPC")) ((TH2F*)fOutputPID->FindObject("hTPCsigvspAC"))->Fill(pid->GetTPCmomentum(),pid->GetTPCsignal());

	  fNEntries->Fill(3);
	} //end analysis cuts
      } //end acceptance and track cuts
    } //end loop on tracks in the candidate
  } //end loop on candidates  
  }
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
  if (!fOutputPID) {     
    printf("ERROR: %s not available\n",GetOutputSlot(2)->GetContainer()->GetName());
    return;
  }

  fOutputTrack = dynamic_cast<TList*> (GetOutputData(3));
  if (!fOutputTrack) {     
    printf("ERROR: %s not available\n",GetOutputSlot(3)->GetContainer()->GetName());
    return;
  }

}

