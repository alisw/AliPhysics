/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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
// AliAnalysisTaskSEHFv2 gives the needed tools for the D 
// mesons v2 analysis with event plane method
// Authors: Chiara Bianchin, cbianchi@pd.infn.it, 
// Robert Grajcarek, grajcarek@physi.uni-heidelberg.de
// Giacomo Ortona, ortona@to.infn.it
// Carlos Eugenio Perez Lara, carlos.eugenio.perez.lara@cern.ch
// 
/////////////////////////////////////////////////////////////

/* $Id$ */

#include <Riostream.h>
#include <TClonesArray.h>
#include <TCanvas.h>
#include <TList.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TDatabasePDG.h>
#include <TRandom3.h>
#include <TVector2.h>
#include <TArrayF.h>

#include <AliLog.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliAnalysisManager.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF4Prong.h"
#include "AliAODRecoCascadeHF.h"

#include "AliAnalysisTaskSE.h"
#include "AliRDHFCutsDplustoKpipi.h"
#include "AliRDHFCutsD0toKpipipi.h"
#include "AliRDHFCutsDstoKKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsLctopKpi.h"

#include "AliHFMassFitter.h"
#include "AliEventplane.h"
#include "AliFlowTrack.h"
#include "AliFlowVector.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"  

#include "AliAnalysisTaskSEHFv2.h"

ClassImp(AliAnalysisTaskSEHFv2)


//________________________________________________________________________
AliAnalysisTaskSEHFv2::AliAnalysisTaskSEHFv2():
AliAnalysisTaskSE(),
  fhEventsInfo(0),
  fOutput(0),
  fRDCuts(0),
  fParHist(0),
  fLowmasslimit(1.765),
  fUpmasslimit(1.965),
  fNPtBins(1),
  fNPhiBinLims(2),
  fPhiBins(0),
  fCentLowLimit(0),
  fCentUpLimit(100),
  fNMassBins(200),
  fReadMC(kFALSE),    
  fUseAfterBurner(kFALSE),
  fDecChannel(0),
  fAfterBurner(0),
  fUseV0EP(kFALSE),
  fV0EPorder(2),
  fMinCentr(10),
  fMaxCentr(80),
  fEtaGap(kFALSE)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEHFv2::AliAnalysisTaskSEHFv2(const char *name,AliRDHFCuts *rdCuts,Int_t decaychannel, Int_t nlimsphibin, Float_t *phibinlimits):
  AliAnalysisTaskSE(name),
  fhEventsInfo(0),
  fOutput(0),
  fRDCuts(rdCuts),
  fParHist(0),
  fLowmasslimit(0),
  fUpmasslimit(0),
  fNPtBins(1),
  fNPhiBinLims(2),
  fPhiBins(0),
  fCentLowLimit(0),
  fCentUpLimit(100),
  fNMassBins(200),
  fReadMC(kFALSE),
  fUseAfterBurner(kFALSE),
  fDecChannel(decaychannel),
  fAfterBurner(0),
  fUseV0EP(kFALSE),
  fV0EPorder(2),
  fMinCentr(10),
  fMaxCentr(80),
  fEtaGap(kFALSE)
{

  Int_t pdg=421;
  switch(fDecChannel){
  case 0:
    pdg=411;
    break;
  case 1:
    pdg=421;
    break;
  case 2:
    pdg=413;
    break;
  case 3:
    pdg=431;
    break;
  case 4:
    pdg=421;
    break;
  case 5:
    pdg=4122;
    break;
  }
  fAfterBurner = new AliHFAfterBurner(fDecChannel);
  if(pdg==413) SetMassLimits((Float_t)0.135,(Float_t)0.165);
  else SetMassLimits((Float_t)0.2,pdg); //check range
  fNPtBins=fRDCuts->GetNPtBins();
  if(nlimsphibin>2) fNPhiBinLims=nlimsphibin;
  else AliInfo("At least 2 limits in Delta phi needed");

  fPhiBins=new Float_t[fNPhiBinLims];
  for(Int_t i=0;i<fNPhiBinLims;i++) fPhiBins[i]=phibinlimits[i];

  if(fDebug>1)fRDCuts->PrintAll();
  // Output slot #1 writes into a TH1F container
  DefineOutput(1,TH1F::Class());   //Info on the number of events etc.
  // Output slot #2 writes into a TList container
  DefineOutput(2,TList::Class());  //Main output
  // Output slot #3 writes into a AliRDHFCuts container (cuts)
  switch(fDecChannel){
  case 0:
    DefineOutput(3,AliRDHFCutsDplustoKpipi::Class());  //Cut object for Dplus
    break;
  case 1:
    DefineOutput(3,AliRDHFCutsD0toKpi::Class());  //Cut object for D0
    break;
  case 2:
    DefineOutput(3,AliRDHFCutsDStartoKpipi::Class());  //Cut object for D*
    break;
  }
  //DefineOutput(4,AliFlowEventSimple::Class());
  //DefineOutput(4,TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskSEHFv2::~AliAnalysisTaskSEHFv2()
{
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if(fhEventsInfo){
    delete fhEventsInfo;
    fhEventsInfo=0;
  }

  if(fRDCuts) {
    delete fRDCuts;
    fRDCuts= 0;
  } 

  if(fParHist) {
    delete fParHist;
    fParHist= 0;
  } 
  for(Int_t i=0;i<6;i++){
    if(fHistvzero[i]) {
      delete fHistvzero[i];
      fHistvzero[i]=0x0;
    }
  } 
  if(fAfterBurner){
    delete fAfterBurner;
    fAfterBurner=0;  
  }
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFv2::SetVZEROParHist(TH2D** histPar){
  for(Int_t i=0;i<6;i++)fHistvzero[i]=(TH2D*)histPar[i]->Clone();
  for(Int_t i=0;i<6;i++){
    if(!fHistvzero[i]){
      printf("No VZERO histograms!\n");
      fUseV0EP=kFALSE;
      return;
    }
  }
  DefineOutput(4,TList::Class());
  fUseV0EP=kTRUE;
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFv2::SetMassLimits(Float_t range, Int_t pdg){
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  fUpmasslimit = mass+range;
  fLowmasslimit = mass-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFv2::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  if(uplimit>lowlimit)
    {
      fUpmasslimit = uplimit;
      fLowmasslimit = lowlimit;
    }
}


//________________________________________________________________________
void AliAnalysisTaskSEHFv2::LocalInit()
{
  // Initialization

  if(fDebug > 1) printf("AnalysisTaskSEHFv2::Init() \n");

  fRDCuts->SetMinCentrality(fMinCentr);
  fRDCuts->SetMaxCentrality(fMaxCentr);

  switch(fDecChannel){
  case 0:
    {
      AliRDHFCutsDplustoKpipi* copycut=new AliRDHFCutsDplustoKpipi(*(static_cast<AliRDHFCutsDplustoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 1:
    {
      AliRDHFCutsD0toKpi* copycut=new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  case 2:
    {
      AliRDHFCutsDStartoKpipi* copycut=new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fRDCuts)));
      // Post the data
      PostData(3,copycut);
    }
    break;
  default:
    return;
  }
  return;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFv2::UserCreateOutputObjects()
{
  // Create the output container
 
  if(fDebug > 1) printf("AnalysisTaskSEHFv2::UserCreateOutputObjects() \n");

  fhEventsInfo=new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Number of AODs scanned",8,-0.5,7.5);
  fhEventsInfo->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fhEventsInfo->GetXaxis()->SetBinLabel(2,"nEvSelected");
  fhEventsInfo->GetXaxis()->SetBinLabel(3,"nCandidatesSelected");
  fhEventsInfo->GetXaxis()->SetBinLabel(4,"out of pt bounds");
  fhEventsInfo->GetXaxis()->SetBinLabel(5,"Pile-up Rej");
  fhEventsInfo->GetXaxis()->SetBinLabel(6,"N. of 0SMH");
  fhEventsInfo->GetXaxis()->SetBinLabel(7,Form("Ev Sel in Centr %.0f-%.0f%s",fRDCuts->GetMinCentrality(),fRDCuts->GetMaxCentrality(),"%"));
  fhEventsInfo->GetXaxis()->SetBinLabel(8,"mismatch lab");
  fhEventsInfo->GetXaxis()->SetNdivisions(1,kFALSE);


  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("MainOutput");
  
  for(Int_t icentr=fMinCentr+5;icentr<=fMaxCentr;icentr=icentr+5){
    TString centrname;centrname.Form("centr%d_%d",icentr-5,icentr);
    for(Int_t i=0;i<fNPtBins;i++){
      
      TString hname;
      TString title;
      hname.Form("hPhi_pt%d",i);hname.Append(centrname.Data());
      title.Form("Phi distribution (Pt bin %d %s);#phi;Entries",i,centrname.Data());
     
      TH1F* hPhi=new TH1F(hname.Data(),title.Data(),96,0.,2*TMath::Pi());
      hPhi->Sumw2();
      fOutput->Add(hPhi);
     
      for(Int_t j=0;j<fNPhiBinLims-1;j++){
       
	hname.Form("hMass_pt%dphi%d",i,j);hname.Append(centrname.Data());
	title.Form("Invariant mass (Pt bin %d, Phi bin %d, %s);M(GeV/c^{2});Entries",i,j,centrname.Data());
       
	TH1F* hMass=new TH1F(hname.Data(),title.Data(),fNMassBins,fLowmasslimit,fUpmasslimit);
	hMass->Sumw2();
	fOutput->Add(hMass);
       
       
	if(fReadMC){
	  hname.Form("hSgn_pt%dphi%d",i,j);hname.Append(centrname.Data());
	  title.Form("Invariant mass S (Pt bin %d, Phi bin %d, %s);M(GeV/c^{2});Entries",i,j,centrname.Data());
	  TH1F* hSgn=new TH1F(hname.Data(),title.Data(),fNMassBins,fLowmasslimit,fUpmasslimit);
	  hSgn->Sumw2();
	  fOutput->Add(hSgn);
	 
	  hname.Form("hBkg_pt%dphi%d",i,j);hname.Append(centrname.Data());
	  title.Form("Invariant mass B (Pt bin %d, Phi bin %d, %s);M(GeV/c^{2});Entries",i,j,centrname.Data());
	  TH1F* hBkg=new TH1F(hname.Data(),title.Data(),fNMassBins,fLowmasslimit,fUpmasslimit);
	  hBkg->Sumw2();
	  fOutput->Add(hBkg);
	 
	  if((fDecChannel != AliAnalysisTaskSEHFv2::kDplustoKpipi) && (fDecChannel != AliAnalysisTaskSEHFv2::kDstartoKpipi)){
	    hname.Form("hRfl_pt%dphi%d",i,j);hname.Append(centrname.Data());
	    title.Form("Invariant mass Reflections (Pt bin %d, Phi bin %d %s);M(GeV/c^{2});Entries",i,j,centrname.Data());
	    TH1F* hRfl=new TH1F(hname.Data(),title.Data(),fNMassBins,fLowmasslimit,fUpmasslimit);
	    hRfl->Sumw2();
	    fOutput->Add(hRfl);
	  }
	}
      }

      TH2F* hMc2phi=new TH2F(Form("hMc2phi_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
      fOutput->Add(hMc2phi);

      if (fReadMC){
	TH2F* hMc2phiS=new TH2F(Form("hMc2phiS_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
	fOutput->Add(hMc2phiS);
	TH2F * hMphiS=new TH2F(Form("hMphiS_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
	fOutput->Add(hMphiS);
	TH2F* hMc2phiB=new TH2F(Form("hMc2phiB_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
	fOutput->Add(hMc2phiB);
	TH2F * hMphiB=new TH2F(Form("hMphiB_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
	fOutput->Add(hMphiB);
	if((fDecChannel != AliAnalysisTaskSEHFv2::kDplustoKpipi) &&(fDecChannel != AliAnalysisTaskSEHFv2::kDstartoKpipi)){
	  TH2F* hMc2phiR=new TH2F(Form("hMc2phiR_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
	  fOutput->Add(hMc2phiR);
	  TH2F* hMphiR=new TH2F(Form("hMphiR_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi (p_{t} bin %d %s);#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),96,0,2*TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
	  fOutput->Add(hMphiR);
	}
      }
    }
  
    TH2F* hMphi=new TH2F(Form("hMphi%s",centrname.Data()),Form("Mass vs #Delta#phi %s;#Delta#phi;M (GeV/c^{2})",centrname.Data()),96,0,TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
    fOutput->Add(hMphi);

    TH1F* hEvPlaneneg=new TH1F(Form("hEvPlaneneg%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
    fOutput->Add(hEvPlaneneg);
	
	TH1F* hEvPlanepos=new TH1F(Form("hEvPlanepos%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
    fOutput->Add(hEvPlanepos);

    //TH1F* hEvPlaneCheck=new TH1F(Form("hEvPlaneCheck%s",centrname.Data()),Form("Event plane angle - Event plane angle per candidate %s;(#phi(Ev Plane) - #phi(Ev Plane Candidate))/#phi(EvPlane);Entries",centrname.Data()),200,-0.2,0.2);
    //fOutput->Add(hEvPlaneCheck);

    TH1F* hEvPlaneCand=new TH1F(Form("hEvPlaneCand%s",centrname.Data()),Form("Event plane angle - Event plane angle per candidate %s;#phi(Ev Plane Candidate);Entries",centrname.Data()),200,-TMath::Pi(),TMath::Pi());
    fOutput->Add(hEvPlaneCand);

    TH1F* hEvPlaneReso=new TH1F(Form("hEvPlaneReso%s",centrname.Data()),Form("Event plane angle Resolution %s;cos2(#psi_{A}-#psi_{B});Entries",centrname.Data()),220,-1.1,1.1);
    fOutput->Add(hEvPlaneReso);
  }
  
  TH1F* hPhiBins=new TH1F("hPhiBins","Bins in #Delta#phi used in this analysis;#phi bin;n jobs",fNPhiBinLims-1,fPhiBins);
  fOutput->Add(hPhiBins);
  for(Int_t k=0;k<fNPhiBinLims-1;k++)hPhiBins->SetBinContent(k+1,1);
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);
  if(fUseV0EP){
    fParHist = new TList();
    fParHist->SetOwner();
    fParHist->SetName("VZEROcorr");
    for(Int_t i=0;i<6;i++){
      fParHist->Add((TH2D*)fHistvzero[i]);
    }
    PostData(4,fParHist);
  }
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEHFv2::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:
  // heavy flavor candidates association to MC truth
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());
  if(fDebug>2) printf("Analysing decay %d\n",fDecChannel);
  // Post the data already here
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

  TClonesArray *arrayProng =0;
  Int_t absPdgMom=0;
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
   
      
      if(fDecChannel==0){
	absPdgMom=411;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Charm3Prong");
      }
      if(fDecChannel==1){
	absPdgMom=421;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("D0toKpi");
      }
      if(fDecChannel==2){
	absPdgMom=413;
	arrayProng=(TClonesArray*)aodFromExt->GetList()->FindObject("Dstar");
      }
    }
  } else {
    if(fDecChannel==0){
      absPdgMom=411;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Charm3Prong");
    }
    if(fDecChannel==1){
      absPdgMom=421;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("D0toKpi");
    }
    if(fDecChannel==2){
      absPdgMom=413;
      arrayProng=(TClonesArray*)aod->GetList()->FindObject("Dstar");
    }
  }

  if(!arrayProng) {
    AliError("AliAnalysisTaskSEHFv2::UserExec:Branch not found!\n");
    return;
  }
  
  // fix for temporary bug in ESDfilter 
  // the AODs with null vertex pointer didn't pass the PhysSel
  if(!aod->GetPrimaryVertex() || TMath::Abs(aod->GetMagneticField())<0.001) return;

  TClonesArray *arrayMC=0;
  AliAODMCHeader *mcHeader=0;
  
  // load MC particles
  if(fReadMC){
    
    arrayMC =  (TClonesArray*)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!arrayMC) {
      AliWarning("AliAnalysisTaskSEHFv2::UserExec:MC particles branch not found!\n");
      //    return;
    }
    
    // load MC header
    mcHeader =  (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if(!mcHeader) {
      AliError("AliAnalysisTaskSEHFv2::UserExec:MC header branch not found!\n");
      return;
    }
  }

  fhEventsInfo->Fill(0); // count event

  AliAODRecoDecayHF *d=0;

  Int_t nCand = arrayProng->GetEntriesFast();
  if(fDebug>2) printf("Number of D2H: %d\n",nCand);

  // trigger class for PbPb C0SMH-B-NOPF-ALLNOTRD, C0SMH-B-NOPF-ALL
  TString trigclass=aod->GetFiredTriggerClasses();
  if(trigclass.Contains("C0SMH-B-NOPF-ALLNOTRD") || trigclass.Contains("C0SMH-B-NOPF-ALL")) fhEventsInfo->Fill(5);

  if(fRDCuts->IsEventSelectedInCentrality(aod)>0) return;
  else fhEventsInfo->Fill(6);
  if(fRDCuts->IsEventSelected(aod)) fhEventsInfo->Fill(1);
  else{
    if(fRDCuts->GetWhyRejection()==1) // rejected for pileup
      fhEventsInfo->Fill(4);
    return;
  }
 
  AliEventplane *pl=0x0;
  TVector2* q=0x0;
  Double_t rpangleevent=0;
  Double_t rpangleeventneg=0;
  Double_t rpangleeventpos=0;
  Double_t eventplane=0;
  TVector2 *qsub1=0x0;
  TVector2 *qsub2=0x0;
  
  //determine centrality bin
  Float_t centr=fRDCuts->GetCentrality(aod);
  Int_t icentr=0;
  for(Int_t ic=fMinCentr+5;ic<=fMaxCentr;ic=ic+5){
    if(ic>centr){
      icentr=ic;
      break;
    }
  }
  TString centrbinname=Form("centr%d_%d",icentr-5,icentr);
  if(fReadMC){
    fUseV0EP=kTRUE;
    TRandom3 *g = new TRandom3(0);
    rpangleevent=g->Rndm()*TMath::Pi();
    delete g;g=0x0;
    eventplane=rpangleevent;
	((TH1F*)fOutput->FindObject(Form("hEvPlanepos%s",centrbinname.Data())))->Fill(rpangleevent);
    if(fUseAfterBurner)fAfterBurner->SetEventPlane((Double_t)rpangleevent);
  }else{
    if(fUseV0EP){
      rpangleevent=GetEventPlaneFromV0(aod);
      eventplane=rpangleevent;
	  ((TH1F*)fOutput->FindObject(Form("hEvPlanepos%s",centrbinname.Data())))->Fill(rpangleevent);
    }else{
      // event plane and resolution 
      //--------------------------------------------------------------------------
      // extracting Q vectors and calculating v2 + resolution
      pl = aod->GetHeader()->GetEventplaneP();
      if(!pl){
	AliError("AliAnalysisTaskSEHFv2::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
	return;
      }
      if(fEtaGap){
        qsub1 = pl->GetQsub1();
   	qsub2 = pl->GetQsub2();
	rpangleeventpos = qsub1->Phi()/2.;
	rpangleeventneg = qsub2->Phi()/2.;
	((TH1F*)fOutput->FindObject(Form("hEvPlanepos%s",centrbinname.Data())))->Fill(rpangleeventpos);
	((TH1F*)fOutput->FindObject(Form("hEvPlaneneg%s",centrbinname.Data())))->Fill(rpangleeventneg);
	}
      else if(!fEtaGap){
	q = pl->GetQVector();
        rpangleevent = pl->GetEventplane("Q");
	((TH1F*)fOutput->FindObject(Form("hEvPlanepos%s",centrbinname.Data())))->Fill(rpangleevent); // reaction plane angle without autocorrelations removal
        } 
       
      Double_t deltaPsi = pl->GetQsubRes();
      if(TMath::Abs(deltaPsi)>TMath::Pi()/2.){
	if(deltaPsi>0.) deltaPsi-=TMath::Pi();
	else deltaPsi +=TMath::Pi();
      } // difference of subevents reaction plane angle cannot be bigger than phi/2
      Double_t planereso = TMath::Cos(2.*deltaPsi); // reaction plane resolution
      //--------------------------------------------------------------------------
      ((TH1F*)fOutput->FindObject(Form("hEvPlaneReso%s",centrbinname.Data())))->Fill(planereso);    
    }
  }
  

  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
    
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(absPdgMom));
    Int_t isSelected= fRDCuts->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    Bool_t isSelBit=kTRUE;
    if(fDecChannel==0) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
    if(fDecChannel==1) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
    if(fDecChannel==2) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0fromDstarCuts);
    if(isSelected && isFidAcc && isSelBit) {
      fhEventsInfo->Fill(2); // candidate selected
      if(fDebug>3) printf("+++++++Is Selected\n");
      
      Float_t* invMass=0x0;
      Int_t nmasses;
      CalculateInvMasses(d,invMass,nmasses);

      Int_t ptbin=fRDCuts->PtBin(d->Pt());
      if(ptbin==-1) {
	fhEventsInfo->Fill(3);
	continue;
      }

      if(!fUseV0EP) {
	eventplane = GetEventPlaneForCandidate(d,q,pl,qsub1,qsub2); // remove autocorrelations
	((TH1F*)fOutput->FindObject(Form("hEvPlaneCand%s",centrbinname.Data())))->Fill(rpangleevent-eventplane);
	//((TH1F*)fOutput->FindObject(Form("hEvPlaneCheck%s",centrbinname.Data())))->Fill((rpangleevent-eventplane)/100.*rpangleevent);
      }

      Float_t phi=d->Phi();
      ((TH1F*)fOutput->FindObject(Form("hPhi_pt%d%s",ptbin,centrbinname.Data())))->Fill(phi);
      
      if(fReadMC&&fUseAfterBurner)phi=fAfterBurner->GetNewAngle(d,arrayMC);
      Float_t deltaphi=GetPhi0Pi(phi-eventplane);

      //fill the histograms with the appropriate method
      if(fDecChannel==0)FillDplus(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
      if(fDecChannel==1)FillD02p(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
      if(fDecChannel==2)FillDstar(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
    }// end if selected
  }
  
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);
  //PostData(4,flowEvent);
  return;
}

//***************************************************************************

// Methods used in the UserExec

void AliAnalysisTaskSEHFv2::CalculateInvMasses(AliAODRecoDecayHF* d,Float_t*& masses,Int_t& nmasses){
  //Calculates all the possible invariant masses for each candidate
  //NB: the implementation for each candidate is responsibility of the corresponding developer

  if(fDecChannel==0){
    //D+ -- Giacomo
    nmasses=1;
    masses=new Float_t[nmasses];
    Int_t pdgdaughters[3] = {211,321,211};
    masses[0]=d->InvMass(3,(UInt_t*)pdgdaughters);
  }
  if(fDecChannel==1){
    //D0 (Kpi)  -- Chiara
    const Int_t ndght=2;
    nmasses=2;
    masses=new Float_t[nmasses];
    Int_t pdgdaughtersD0[ndght]={211,321};//pi,K 
    masses[0]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0); //D0
    Int_t pdgdaughtersD0bar[ndght]={321,211};//K,pi 
    masses[1]=d->InvMass(ndght,(UInt_t*)pdgdaughtersD0bar); //D0bar
  }
  if(fDecChannel==2){
    //D* -- Robert,Yifei, Alessandro
    nmasses=1;
    masses=new Float_t[nmasses];
    masses[0]=((AliAODRecoCascadeHF*)d)->DeltaInvMass();
  } 
}

//******************************************************************************

//Methods to fill the istograms with MC information, one for each candidate
//NB: the implementation for each candidate is responsibility of the corresponding developer

void AliAnalysisTaskSEHFv2::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi,Float_t* masses,Int_t isSel,Int_t icentr){
  //D+ channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
 
  Int_t phibin=GetPhiBin(deltaphi);

  ((TH1F*)fOutput->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[0]);
  Int_t pdgdaughters[3] = {211,321,211};

  if(fReadMC){
    Int_t lab=-1;
    if(fUseAfterBurner){
      Bool_t isSignal=fAfterBurner->GetIsSignal();
      if(isSignal)lab=10;
    }else {
      lab = d->MatchToMC(411,arrayMC,3,pdgdaughters);
    }
    if(lab>=0){ //signal
      ((TH1F*)fOutput->FindObject(Form("hSgn_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH1F*)fOutput->FindObject(Form("hBkg_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }   
}


void AliAnalysisTaskSEHFv2::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi,Float_t* masses,Int_t isSel,Int_t icentr){

  //D0->Kpi channel

  //mass histograms
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t phibin=GetPhiBin(deltaphi);
  if(isSel==1 || isSel==3) {
    ((TH1F*)fOutput->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[0]);
  }
  if(isSel>=2) {
    ((TH1F*)fOutput->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[1]);
    ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
    ((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[1]);
  }


  //MC histograms
  if(fReadMC){

    Int_t matchtoMC=-1;

    //D0
    Int_t pdgdaughters[2];
    pdgdaughters[0]=211;//pi 
    pdgdaughters[1]=321;//K
    Int_t nprongs=2;
    Int_t absPdgMom=421;

    matchtoMC = d->MatchToMC(absPdgMom,arrayMC,nprongs,pdgdaughters);

    Int_t prongPdgPlus=421,prongPdgMinus=(-1)*421;
    if((isSel==1 || isSel==3)){ //D0
      if(matchtoMC>=0){
	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	
	if(pdgMC==prongPdgPlus) {
	  ((TH1F*)fOutput->FindObject(Form("hSgn_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
	}
	else {
	  ((TH1F*)fOutput->FindObject(Form("hRfl_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiR_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiRcentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[0]);
	}
      } else {
	((TH1F*)fOutput->FindObject(Form("hBkg_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
	((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[0]);
      }
    }
    if(isSel>=2){ //D0bar
      if(matchtoMC>=0){
	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	
	if(pdgMC==prongPdgMinus) {
	  ((TH1F*)fOutput->FindObject(Form("hSgn_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
	}
	else {
	  ((TH1F*)fOutput->FindObject(Form("hRfl_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiR_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[1]);
	}
      } else {
	((TH1F*)fOutput->FindObject(Form("hBkg_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[1]);
	((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
      }
    }
  }
}

void AliAnalysisTaskSEHFv2::FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, Float_t* masses,Int_t isSel,Int_t icentr){
  //D* channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  Int_t phibin=GetPhiBin(deltaphi);
  
  ((TH1F*)fOutput->FindObject(Form("hMass_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMphicentr%d_%d",icentr-5,icentr)))->Fill(deltaphi,masses[0]);
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  if(fReadMC){
    Int_t lab=-1;
    lab = ((AliAODRecoCascadeHF*)d)->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,arrayMC);
    if(lab>=0){ //signal
      ((TH1F*)fOutput->FindObject(Form("hSgn_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH1F*)fOutput->FindObject(Form("hBkg_pt%dphi%dcentr%d_%d",ptbin,phibin,icentr-5,icentr)))->Fill(masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }
  

}


//________________________________________________________________________
Int_t AliAnalysisTaskSEHFv2::GetPhiBin(Float_t deltaphi){

  //give the bin corresponding to the value of deltaphi according to the binning requested in the constructor
  Int_t phibin=0;
  for(Int_t i=0;i<fNPhiBinLims-1;i++) {
    if(deltaphi>=fPhiBins[i] && deltaphi<fPhiBins[i+1]) {
      phibin=i;
      break;
    }
  }
  return phibin;
}
//________________________________________________________________________
// Float_t AliAnalysisTaskSEHFv2::GetPhi02Pi(Float_t phi){
//   Float_t result=phi;
//   while(result<0){
//     result=result+2*TMath::Pi();
//   }
//   while(result>TMath::Pi()*2){
//     result=result-2*TMath::Pi();
//   }
//   return result;
// }
//________________________________________________________________________
Float_t AliAnalysisTaskSEHFv2::GetPhi0Pi(Float_t phi){
  Float_t result=phi;
  while(result<0){
    result=result+TMath::Pi();
  }
  while(result>TMath::Pi()){
    result=result-TMath::Pi();
  }
  return result;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSEHFv2::GetEventPlaneForCandidate(AliAODRecoDecayHF* d, TVector2* q,AliEventplane *pl,TVector2* qsub1,TVector2* qsub2){
  // remove autocorrelations 
 
  TArrayF* qx = 0x0;
  TArrayF* qy = 0x0;
  TVector2 qcopy; 
  if(!fEtaGap){
    qx = pl->GetQContributionXArray();
    qy = pl->GetQContributionYArray();
    qcopy = *q;
    }
  else {
    if(d->Eta()>0.){
      qx = pl->GetQContributionXArraysub1();
      qy = pl->GetQContributionYArraysub1();
      qcopy = *qsub1;
    }
    else{
      qx = pl->GetQContributionXArraysub2();
      qy = pl->GetQContributionYArraysub2();
      qcopy = *qsub2;
    }
  }
  
  
 
  if(fDecChannel==2){
    //D* -- Yifei, Alessandro,Robert
    AliAODRecoDecayHF2Prong* theD0particle = ((AliAODRecoCascadeHF*)d)->Get2Prong();
    AliAODTrack *track0 = (AliAODTrack*)theD0particle->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)theD0particle->GetDaughter(1);  
    AliAODTrack *track2 = ((AliAODRecoCascadeHF*)d)->GetBachelor(); 
    // reduce global q vector

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    TVector2 q2;
    if((track2->GetID()) < qx->fN){
      q2.Set(qx->At(track2->GetID()),qy->At(track2->GetID()));}
      
    qcopy = qcopy -(q0+q1+q2);
  
  }
  
  // reduce Q vector for D+ and D0
  
  if(fDecChannel==1){    
    //D0 -- Chiara
    AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1);  

    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    qcopy = qcopy -(q0+q1);
  }

  if(fDecChannel==0){
    //D+ -- Giacomo
    AliAODTrack *track0 = (AliAODTrack*)d->GetDaughter(0);
    AliAODTrack *track1 = (AliAODTrack*)d->GetDaughter(1);  
    AliAODTrack *track2 = (AliAODTrack*)d->GetDaughter(2);  
    
    TVector2 q0;
    if((track0->GetID()) < qx->fN){
      q0.Set(qx->At(track0->GetID()),qy->At(track0->GetID()));}
	
    TVector2 q1;
    if((track1->GetID()) < qx->fN){
      q1.Set(qx->At(track1->GetID()),qy->At(track1->GetID()));}
	
    TVector2 q2;
    if((track2->GetID()) < qx->fN){
      q2.Set(qx->At(track2->GetID()),qy->At(track2->GetID()));}
      
    qcopy = qcopy -(q0+q1+q2);
	
  }

  return qcopy.Phi()/2.;

}
//________________________________________________________________________
Float_t AliAnalysisTaskSEHFv2::GetEventPlaneFromV0(AliAODEvent *aodEvent){

  Int_t centr=fRDCuts->GetCentrality(aodEvent);
  centr=centr-centr%10;
  //temporary fix
  if(centr<20)centr=20;
  if(centr>70)centr=70;
  //end temporary fix
  Int_t binx=0;
  Int_t iParHist=(centr-20)/10;

  TString name;name.Form("parhist%d_%d",centr,centr+10);

  if(fDebug>15)printf("EPfromV0 centr %d, iparhist %d (%p-%p)\n",centr,iParHist,fParHist->FindObject(name.Data()),fParHist->At(iParHist));

  Int_t runnumber=aodEvent->GetRunNumber();
  if(fParHist->At(iParHist)){
    for(Int_t i=1;i<=((TH2D*)fParHist->At(iParHist))->GetNbinsX()&&binx<=0;i++){
      Int_t run=atoi(((TH2D*)fParHist->At(iParHist))->GetXaxis()->GetBinLabel(i));
      if(run>=runnumber)binx=i;
    }
  }else{
    fhEventsInfo->Fill(7);
  }

  AliFlowTrackCuts* cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
  cutsRP->SetEvent(aodEvent, MCEvent());//, 0x0);
  cutsRP->SetName( Form("rp_cuts") );
  AliFlowTrackCuts* dummy = new AliFlowTrackCuts("null_cuts");
  dummy->SetParamType(AliFlowTrackCuts::kGlobal);
  dummy->SetPtRange(+1,-1); // select nothing QUICK
  dummy->SetEtaRange(+1,-1); // select nothing VZERO
  dummy->SetEvent(aodEvent,MCEvent());

  //////////////// construct the flow event container ////////////
  AliFlowEvent flowEvent(cutsRP,dummy);
  flowEvent.SetReferenceMultiplicity( 64 );
  for(Int_t i=0;i<64&&binx>0;i++){
    AliFlowTrack *flowTrack=flowEvent.GetTrack(i);
    Double_t inte=((TH2D*)fParHist->At(iParHist))->Integral(binx,binx,i+1,i+1);
    if(inte>0)flowTrack->SetWeight(flowTrack->Weight()/inte);
  }
  if(fDebug>15)printf("EPfromV0 flow tracks weights done\n");

  AliFlowVector qvec=flowEvent.GetQ(fV0EPorder);
  Double_t angleEP=(1./(Double_t)fV0EPorder)*qvec.Phi();
  if(fDebug>15)printf("EPfromV0 phi %f\n",angleEP);
  return angleEP;
}
//________________________________________________________________________
void AliAnalysisTaskSEHFv2::Terminate(Option_t */*option*/)
{
  
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEHFv2: Terminate() \n");
  /*
    fhEventsInfo = dynamic_cast<TH1F*> (GetOutputData(1));
    if(!fhEventsInfo){
    printf("Error: hEventsInfo not available\n");
    return;
    }

    fOutput = dynamic_cast<TList*> (GetOutputData(2));
    if (!fOutput) {     
    printf("ERROR: fOutput not available\n");
    return;
    }
    switch(fDecChannel){
    case 0:
    fRDCuts = dynamic_cast<AliRDHFCutsDplustoKpipi*> (GetOutputData(3));
    break;
    case 1:
    fRDCuts = dynamic_cast<AliRDHFCutsD0toKpi*> (GetOutputData(3));
    break;
    case 2:
    fRDCuts = dynamic_cast<AliRDHFCutsDStartoKpipi*> (GetOutputData(3));
    break;
    default:
    break;
    }
    if (!fRDCuts) {     
    printf("ERROR: fRDCuts not available\n");
    return;
    }
   
    // TCanvas* cvex=new TCanvas("example","Example Mass");
    // cvex->cd();
    // ((TH1F*)fOutput->FindObject("hMass_pt0phi0"))->Draw();
    // TCanvas* cv2d=new TCanvas("2d","mass-cos2phi");
    // cv2d->cd();
    // ((TH2F*)fOutput->FindObject("hMc2phi"))->Draw("colz");
    // TCanvas *cstat=new TCanvas("cstat","Stat");
    // cstat->SetGridy();
    // cstat->cd();
    // fhEventsInfo->Draw("htext0");
   
    //  TMultiGraph *multig = new TMultiGraph();
    if(fDecChannel==2)return;
    TGraphErrors *g[fNPtBins];
    TH1F *h = new TH1F("h","h",100,0.,1.);
    TString hname;
    TString gname;
    for(Int_t ipt = 0;ipt<fNPtBins;ipt++){
    g[ipt] = new TGraphErrors(fNPhiBinLims);
    gname.Form("hMass_pt%d",ipt);
    g[ipt]->SetTitle(gname.Data());
    for(Int_t iphi = 0;iphi<fNPhiBinLims;iphi++){  
    hname.Form("hMass_pt%dphi%d",ipt,iphi);
    h = (TH1F*)fOutput->FindObject("hMass_pt0phi0");
    AliHFMassFitter fitter(h,fLowmasslimit,fUpmasslimit,2,0);
    Int_t pdg=0;
    switch(fDecChannel){
    case 0:
    pdg=411;
    break;
    case 1:
    pdg=421;
    break;
    case 2:
    pdg=413;
    break;
    default:
    break;
    }
    fitter.SetInitialGaussianMean(TDatabasePDG::Instance()->GetParticle(pdg)->Mass());
    fitter.SetInitialGaussianSigma(0.012);
    fitter.InitNtuParam(Form("ntuPtbin%d",ipt));
    Double_t signal=0, errSignal=0;
    if(fitter.MassFitter(kFALSE)){
    fitter.Signal(3,signal,errSignal);
    }
    g[ipt]->SetPoint(iphi,fPhiBins[iphi],signal);
    g[ipt]->SetPointError(iphi,fPhiBins[iphi],errSignal);
    }//end loop over phi
     //     multig->Add(g[ipt],"ap");
     }//end loop on pt
     TCanvas *cdndphi = new TCanvas("dN/d#phi","dN/d#phi");
     cdndphi->Divide(1,fNPtBins); 
     for(Int_t ipt = 0;ipt<fNPtBins;ipt++){
     cdndphi->cd(ipt+1);
     g[ipt]->Draw("AP");
     }
  */
  return;
}

//-------------------------------------------
/*
  Float_t GetEventPlaneFromVZERO(){
  AliAODHFUtil *info = (AliAODHFUtil*) aodEvent->GetList()->FindObject("fHFUtilInfo");
  if (!info) return -999.;
  //============= FIX ONLY FOR AOD033
  Double_t *par0;
  Double_t par0_137161[64] = { 6.71e-02 , 6.86e-02 , 7.06e-02 , 6.32e-02 , 
  5.91e-02 , 6.07e-02 , 5.78e-02 , 5.73e-02 , 5.91e-02 , 6.22e-02 , 
  5.90e-02 , 6.11e-02 , 5.55e-02 , 5.29e-02 , 5.19e-02 , 5.56e-02 , 
  6.25e-02 , 7.03e-02 , 5.64e-02 , 5.81e-02 , 4.57e-02 , 5.30e-02 , 
  5.13e-02 , 6.43e-02 , 6.27e-02 , 6.48e-02 , 6.07e-02 , 1.01e-01 , 
  6.68e-02 , 7.16e-02 , 6.36e-02 , 5.95e-02 , 2.52e-02 , 2.82e-02 , 
  2.56e-02 , 2.86e-02 , 2.82e-02 , 2.10e-02 , 2.13e-02 , 2.32e-02 , 
  2.75e-02 , 4.34e-02 , 3.78e-02 , 4.52e-02 , 4.11e-02 , 3.89e-02 , 
  4.10e-02 , 3.73e-02 , 4.51e-02 , 5.07e-02 , 5.42e-02 , 4.74e-02 , 
  4.33e-02 , 4.44e-02 , 4.64e-02 , 3.01e-02 , 6.38e-02 , 5.26e-02 , 
  4.99e-02 , 5.26e-02 , 5.47e-02 , 3.84e-02 , 5.00e-02 , 5.20e-02 };
  Double_t par0_137366[64] = { 7.12e-02 , 7.34e-02 , 7.39e-02 , 6.54e-02 , 6.11e-02 , 6.31e-02 , 6.15e-02 , 
  6.00e-02 , 6.10e-02 , 6.49e-02 , 6.17e-02 , 6.33e-02 , 6.00e-02 , 5.48e-02 , 
  5.44e-02 , 5.81e-02 , 6.49e-02 , 7.07e-02 , 5.91e-02 , 6.18e-02 , 4.82e-02 , 
  5.67e-02 , 5.36e-02 , 6.60e-02 , 6.37e-02 , 6.78e-02 , 6.31e-02 , 1.04e-01 , 
  6.91e-02 , 7.32e-02 , 6.61e-02 , 6.16e-02 , 2.64e-02 , 2.81e-02 , 2.64e-02 , 
  2.85e-02 , 2.87e-02 , 2.18e-02 , 2.19e-02 , 2.43e-02 , 2.81e-02 , 4.37e-02 , 
  3.90e-02 , 4.66e-02 , 4.24e-02 , 4.09e-02 , 4.21e-02 , 3.88e-02 , 4.83e-02 , 
  5.23e-02 , 5.44e-02 , 4.85e-02 , 4.42e-02 , 4.58e-02 , 4.74e-02 , 3.14e-02 , 
  6.31e-02 , 5.30e-02 , 5.01e-02 , 5.33e-02 , 5.70e-02 , 3.95e-02 , 4.98e-02 , 5.31e-02 };
  if (aodEvent->GetRunNumber() <= 137165)
  par0=par0_137161;
  else
  par0=par0_137366;
  Float_t vChCorr[64];
  for(int i=0; i!=64; ++i)
  vChCorr[i] = (info->GetVZEROChannel(i))/par0[i]/64.;
  //============= END OF FIX AOD033
  Float_t multR[8];
  for(int i=0; i!=8; ++i) multR[i] = 0;
  for(int i=0; i!=64; ++i)
  multR[i/8] += vChCorr[i];
  for(int i=0; i!=8; ++i) 
  if(multR[i]) {
  double Qx=0, Qy=0;
  for(int j=0; j!=8; ++j) {
  double phi = TMath::PiOver4()*(0.5+j);
  Qx += TMath::Cos(2*phi)*vChCorr[i*8+j]/multR[i];
  Qy += TMath::Sin(2*phi)*vChCorr[i*8+j]/multR[i];
  }
  }
  return 0.5*TMath::ATan2(Qy,Qx)+TMath::PiOver2();
  }

*/
