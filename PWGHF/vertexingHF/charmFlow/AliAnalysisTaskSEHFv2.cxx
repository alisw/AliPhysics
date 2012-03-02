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
  fLowmasslimit(1.669),
  fUpmasslimit(2.069),
  fNPtBins(1),
  fNMassBins(200),
  fReadMC(kFALSE),    
  fUseAfterBurner(kFALSE),
  fDecChannel(0),
  fAfterBurner(0),
  fEventPlaneMeth(kTPCVZERO),
  fEventPlanesComp(10),
  fV0EPorder(2),
  fMinCentr(20),
  fMaxCentr(80),
  fEtaGap(kFALSE)
{
  // Default constructor
}

//________________________________________________________________________
AliAnalysisTaskSEHFv2::AliAnalysisTaskSEHFv2(const char *name,AliRDHFCuts *rdCuts,Int_t decaychannel):
  AliAnalysisTaskSE(name),
  fhEventsInfo(0),
  fOutput(0),
  fRDCuts(rdCuts),
  fLowmasslimit(0),
  fUpmasslimit(0),
  fNPtBins(1),
  fNMassBins(200),
  fReadMC(kFALSE),
  fUseAfterBurner(kFALSE),
  fDecChannel(decaychannel),
  fAfterBurner(0),
  fEventPlaneMeth(kTPCVZERO),
  fEventPlanesComp(10),
  fV0EPorder(2),
  fMinCentr(20),
  fMaxCentr(80),
  fEtaGap(kFALSE)
{
  // standard constructor
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
  delete fOutput;
  delete fhEventsInfo;
  delete fRDCuts;
  delete fAfterBurner;
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFv2::SetMassLimits(Float_t range, Int_t pdg){
  // Set limits for mass spectra plots
  Float_t mass=0;
  Int_t abspdg=TMath::Abs(pdg);
  mass=TDatabasePDG::Instance()->GetParticle(abspdg)->Mass();
  fUpmasslimit = mass+range;
  fLowmasslimit = mass-range;
}
//_________________________________________________________________
void  AliAnalysisTaskSEHFv2::SetMassLimits(Float_t lowlimit, Float_t uplimit){
  // Set limits for mass spectra plots
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

  fhEventsInfo=new TH1F(GetOutputSlot(1)->GetContainer()->GetName(), "Number of AODs scanned",6,-0.5,5.5);
  fhEventsInfo->GetXaxis()->SetBinLabel(1,"nEventsAnal");
  fhEventsInfo->GetXaxis()->SetBinLabel(2,"nEvSelected");
  fhEventsInfo->GetXaxis()->SetBinLabel(3,"nCandidatesSelected");
  fhEventsInfo->GetXaxis()->SetBinLabel(4,"out of pt bounds");
  fhEventsInfo->GetXaxis()->SetBinLabel(5,Form("Ev Sel in Centr %.0f-%.0f%s",fRDCuts->GetMinCentrality(),fRDCuts->GetMaxCentrality(),"%"));
  fhEventsInfo->GetXaxis()->SetBinLabel(6,"mismatch lab");
  fhEventsInfo->GetXaxis()->SetNdivisions(1,kFALSE);


  // Several histograms are more conveniently managed in a TList
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("MainOutput");
  
  for(Int_t icentr=fMinCentr+5;icentr<=fMaxCentr;icentr=icentr+5){
    TString centrname;centrname.Form("centr%d_%d",icentr-5,icentr);

      TH2F* hMPtCand=new TH2F(Form("hMPtCand%s",centrname.Data()),Form("Mass vs pt %s;pt (GeV);M (GeV/c^{2})",centrname.Data()),120,0,24.,fNMassBins,fLowmasslimit,fUpmasslimit);
      fOutput->Add(hMPtCand);//For <pt> calculation


    //Candidate distributions
    for(Int_t i=0;i<fNPtBins;i++){ 
      TH2F* hMc2phi=new TH2F(Form("hMc2phi_pt%d%s",i,centrname.Data()),Form("Mass vs cos2#Delta#phi (p_{t} bin %d %s);cos2#Delta#phi;M (GeV/c^{2})",i,centrname.Data()),100,-1.,1.,fNMassBins,fLowmasslimit,fUpmasslimit);
      fOutput->Add(hMc2phi);//for 2D analysis
      
      TH2F* hMphi=new TH2F(Form("hMphi_pt%d%s",i,centrname.Data()),Form("Mass vs #Delta#phi %s;#Delta#phi;M (GeV/c^{2})",centrname.Data()),96,0,TMath::Pi(),fNMassBins,fLowmasslimit,fUpmasslimit);
      fOutput->Add(hMphi);//for phi bins analysis
      
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


    //Event Plane
    TH2F* hEvPlane=new TH2F(Form("hEvPlane%s",centrname.Data()),Form("VZERO/TPC Event plane angle %s;#phi Ev Plane (TPC);#phi Ev Plane (VZERO);Entries",centrname.Data()),200,0.,TMath::Pi(),200,0.,TMath::Pi());
    fOutput->Add(hEvPlane);

    TH1F* hEvPlaneneg=new TH1F(Form("hEvPlaneneg%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
    fOutput->Add(hEvPlaneneg);
	
    TH1F* hEvPlanepos=new TH1F(Form("hEvPlanepos%s",centrname.Data()),Form("Event plane angle %s;#phi Ev Plane;Entries",centrname.Data()),200,0.,TMath::Pi());
    fOutput->Add(hEvPlanepos);

    TH1F* hEvPlaneCand=new TH1F(Form("hEvPlaneCand%s",centrname.Data()),Form("Event plane angle - Event plane angle per candidate %s;#phi(Ev Plane Candidate);Entries",centrname.Data()),200,-TMath::Pi(),TMath::Pi());
    fOutput->Add(hEvPlaneCand);

    TH1F* hEvPlaneReso=new TH1F(Form("hEvPlaneReso%s",centrname.Data()),Form("Event plane angle Resolution %s;cos2(#psi_{A}-#psi_{B});Entries",centrname.Data()),220,-1.1,1.1);
    fOutput->Add(hEvPlaneReso);
  }
  
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

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

  Bool_t isEvSel=fRDCuts->IsEventSelected(aod);
  if(!isEvSel){
    if(!fRDCuts->IsEventRejectedDueToCentrality())fhEventsInfo->Fill(4);
    return;
  }

  fhEventsInfo->Fill(1);
  AliEventplane *pl=aod->GetEventplane();
  if(!pl){
    AliError("AliAnalysisTaskSEHFv2::UserExec:no eventplane! v2 analysis without eventplane not possible!\n");
    return;
  }
  
  //Event plane
  Double_t eventplane=0;
  Double_t rpangleTPC=0;
  Double_t rpangleVZERO=0;
  Double_t planereso=0;
  Double_t deltaPsi=0;
  Double_t rpangleeventneg=0;
  Double_t rpangleeventpos=0;

  //For candidate removal from TPC EP
  TVector2* q=0x0;
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
    fEventPlaneMeth=kVZERO;
    TRandom3 *g = new TRandom3(0);
    rpangleVZERO=g->Rndm()*TMath::Pi();
    delete g;g=0x0;
    eventplane=rpangleVZERO;
    if(fUseAfterBurner)fAfterBurner->SetEventPlane((Double_t)eventplane);
  }else{
    if(fEventPlaneMeth!=kTPC){
      //VZERO EP and resolution
      rpangleVZERO=pl->GetEventplane("V0",aod,fV0EPorder);
      if(fEventPlaneMeth!=kTPCVZERO){
	//Using V0A/C to determine resolution
	rpangleeventneg= pl->GetEventplane("V0A",aod,fV0EPorder);
	rpangleeventpos= pl->GetEventplane("V0C",aod,fV0EPorder);
	deltaPsi =rpangleeventneg-rpangleeventpos;
	if(TMath::Abs(deltaPsi)>TMath::Pi()/2.){// difference of subevents reaction plane angle cannot be bigger than phi/2
	  if(deltaPsi>0.) deltaPsi-=TMath::Pi();
	  else deltaPsi +=TMath::Pi();
	} 
	eventplane=rpangleVZERO;
      }
    }
    if(fEventPlaneMeth!=kVZERO){
      // TPC event plane and resolution 

      // extracting Q vectors and calculating v2 + resolution
      if(fEtaGap){
        qsub1 = pl->GetQsub1();
   	qsub2 = pl->GetQsub2();
	rpangleeventpos = qsub1->Phi()/2.;
	rpangleeventneg = qsub2->Phi()/2.;
      }
      else if(!fEtaGap){
	q = pl->GetQVector();
        rpangleTPC = pl->GetEventplane("Q");
      } 
      if(fEventPlaneMeth!=kVZEROTPC){
	deltaPsi = pl->GetQsubRes();
	if(TMath::Abs(deltaPsi)>TMath::Pi()/2.){
	  if(deltaPsi>0.) deltaPsi-=TMath::Pi();
	  else deltaPsi +=TMath::Pi();
	} // difference of subevents reaction plane angle cannot be bigger than phi/2
	planereso = TMath::Cos(2.*deltaPsi); // reaction plane resolution
      }
    }
  }
  if(TMath::Abs(rpangleTPC-rpangleVZERO)>fEventPlanesComp)return;
  //Filling EP-related histograms
  planereso = TMath::Cos(2.*deltaPsi); // reaction plane resolution
  ((TH2F*)fOutput->FindObject(Form("hEvPlane%s",centrbinname.Data())))->Fill(rpangleTPC,rpangleVZERO); // reaction plane angle without autocorrelations removal
  ((TH1F*)fOutput->FindObject(Form("hEvPlaneReso%s",centrbinname.Data())))->Fill(planereso); //RP resolution   
  ((TH1F*)fOutput->FindObject(Form("hEvPlanepos%s",centrbinname.Data())))->Fill(rpangleeventpos); //Angle of first subevent
  ((TH1F*)fOutput->FindObject(Form("hEvPlaneneg%s",centrbinname.Data())))->Fill(rpangleeventneg); //Angle of second subevent

  //Loop on D candidates
  for (Int_t iCand = 0; iCand < nCand; iCand++) {
    d=(AliAODRecoDecayHF*)arrayProng->UncheckedAt(iCand);
    Bool_t isSelBit=kTRUE;
    if(fDecChannel==0) isSelBit=d->HasSelectionBit(AliRDHFCuts::kDplusCuts);
    if(fDecChannel==1) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0toKpiCuts);
    if(fDecChannel==2) isSelBit=d->HasSelectionBit(AliRDHFCuts::kD0fromDstarCuts);
    if(!isSelBit)continue;
    Int_t ptbin=fRDCuts->PtBin(d->Pt());
    if(ptbin<0) {
      fhEventsInfo->Fill(3);
      continue;
    }
    Bool_t isFidAcc = fRDCuts->IsInFiducialAcceptance(d->Pt(),d->Y(absPdgMom));
    if(!isFidAcc)continue;    
    Int_t isSelected= fRDCuts->IsSelected(d,AliRDHFCuts::kCandidate,aod);
    if(!isSelected)continue;

    fhEventsInfo->Fill(2); // candidate selected
    if(fDebug>3) printf("+++++++Is Selected\n");
      
    Float_t* invMass=0x0;
    Int_t nmasses;
    CalculateInvMasses(d,invMass,nmasses);

    if(fEventPlaneMeth%2==0){
      eventplane = GetEventPlaneForCandidate(d,q,pl,qsub1,qsub2); // remove autocorrelations
      ((TH1F*)fOutput->FindObject(Form("hEvPlaneCand%s",centrbinname.Data())))->Fill(rpangleTPC-eventplane);
    }

    Double_t phi=d->Phi(); 
    if(fReadMC&&fUseAfterBurner)phi=fAfterBurner->GetNewAngle(d,arrayMC);
    Float_t deltaphi=GetPhi0Pi(phi-eventplane);
    
    //fill the histograms with the appropriate method
    if(fDecChannel==0)FillDplus(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
    else if(fDecChannel==1)FillD02p(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
    else if(fDecChannel==2)FillDstar(d,arrayMC,ptbin,deltaphi,invMass,isSelected,icentr);
    
    delete invMass;
  }
  
  PostData(1,fhEventsInfo);
  PostData(2,fOutput);

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

//Methods to fill the histograms, one for each channel
//NB: the implementation for each candidate is responsibility of the corresponding developer

//******************************************************************************
void AliAnalysisTaskSEHFv2::FillDplus(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr){
  //D+ channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
 
  ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentr-5,icentr)))->Fill(d->Pt(),masses[0]);
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
      ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }   
}

//******************************************************************************
void AliAnalysisTaskSEHFv2::FillD02p(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr){

  //D0->Kpi channel

  //mass histograms
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }
  if(isSel==1 || isSel==3) {
    ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentr-5,icentr)))->Fill(d->Pt(),masses[0]);
  }
  if(isSel>=2) {
    ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
    ((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
    ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentr-5,icentr)))->Fill(d->Pt(),masses[1]);
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
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
	}
	else {
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiR_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiR_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
	}
      } else {
	((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
	((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
      }
    }
    if(isSel>=2){ //D0bar
      if(matchtoMC>=0){
	AliAODMCParticle *dMC = (AliAODMCParticle*)arrayMC->At(matchtoMC);
	Int_t pdgMC = dMC->GetPdgCode();
	
	if(pdgMC==prongPdgMinus) {
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
	}
	else {
	  ((TH2F*)fOutput->FindObject(Form("hMc2phiR_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	  ((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
	}
      } else {
	((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[1]);
	((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[1]);
      }
    }
  }
}
//******************************************************************************
void AliAnalysisTaskSEHFv2::FillDstar(AliAODRecoDecayHF* d,TClonesArray *arrayMC,Int_t ptbin,Float_t deltaphi, const Float_t* masses,Int_t isSel,Int_t icentr){
  //D* channel
  if(!isSel){
    if(fDebug>3)AliWarning("Candidate not selected\n");
    return;
  }
  if(!masses){
    if(fDebug>3)AliWarning("Masses not calculated\n");
    return;
  }

  ((TH2F*)fOutput->FindObject(Form("hMc2phi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMphi_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
  ((TH2F*)fOutput->FindObject(Form("hMPtCandcentr%d_%d",icentr-5,icentr)))->Fill(d->Pt(),masses[0]);
  Int_t pdgDgDStartoD0pi[2]={421,211};
  Int_t pdgDgD0toKpi[2]={321,211};
  
  if(fReadMC){
    Int_t lab=-1;
    lab = ((AliAODRecoCascadeHF*)d)->MatchToMC(413,421,pdgDgDStartoD0pi,pdgDgD0toKpi,arrayMC);
    if(lab>=0){ //signal
      ((TH2F*)fOutput->FindObject(Form("hMc2phiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiS_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } else{ //background
      ((TH2F*)fOutput->FindObject(Form("hMc2phiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(TMath::Cos(2*deltaphi),masses[0]);
      ((TH2F*)fOutput->FindObject(Form("hMphiB_pt%dcentr%d_%d",ptbin,icentr-5,icentr)))->Fill(deltaphi,masses[0]);
    } 
  }
}

//________________________________________________________________________
void AliAnalysisTaskSEHFv2::SetEventPlaneMethod(Int_t method){
  if(method>3||method<0){
    AliWarning("No EP method associated to the selection, setting to TPC EP\n");
    method=kTPCVZERO;
  }
  fEventPlaneMeth=method;
}

//________________________________________________________________________
Float_t AliAnalysisTaskSEHFv2::GetPhi0Pi(Float_t phi){
  // Sets the phi angle in the range 0-pi
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
Float_t AliAnalysisTaskSEHFv2::GetEventPlaneForCandidate(AliAODRecoDecayHF* d, const TVector2* q,AliEventplane *pl, const TVector2* qsub1, const TVector2* qsub2){
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
// //________________________________________________________________________
// Float_t AliAnalysisTaskSEHFv2::GetEventPlaneFromV0(AliAODEvent *aodEvent){
//   // Compute event plane for VZERO - Obsolete, used for 2010 data

//   Int_t centr=fRDCuts->GetCentrality(aodEvent);
//   centr=centr-centr%10;
//   //temporary fix
//   if(centr<20)centr=20;
//   if(centr>70)centr=70;
//   //end temporary fix
//   Int_t binx=0;
//   Int_t iParHist=(centr-20)/10;

//   TString name;name.Form("parhist%d_%d",centr,centr+10);

//   if(fDebug>15)printf("EPfromV0 centr %d, iparhist %d (%p-%p)\n",centr,iParHist,fParHist->FindObject(name.Data()),fParHist->At(iParHist));

//   Int_t runnumber=aodEvent->GetRunNumber();
//   if(fParHist->At(iParHist)){
//     for(Int_t i=1;i<=((TH2D*)fParHist->At(iParHist))->GetNbinsX()&&binx<=0;i++){
//       Int_t run=atoi(((TH2D*)fParHist->At(iParHist))->GetXaxis()->GetBinLabel(i));
//       if(run>=runnumber)binx=i;
//     }
//   }else{
//     fhEventsInfo->Fill(7);
//   }

//   AliFlowTrackCuts* cutsRP = AliFlowTrackCuts::GetStandardVZEROOnlyTrackCuts();
//   cutsRP->SetEvent(aodEvent, MCEvent());//, 0x0);
//   cutsRP->SetName( Form("rp_cuts") );
//   AliFlowTrackCuts* dummy = new AliFlowTrackCuts("null_cuts");
//   dummy->SetParamType(AliFlowTrackCuts::kGlobal);
//   dummy->SetPtRange(+1,-1); // select nothing QUICK
//   dummy->SetEtaRange(+1,-1); // select nothing VZERO
//   dummy->SetEvent(aodEvent,MCEvent());

//   //////////////// construct the flow event container ////////////
//   AliFlowEvent flowEvent(cutsRP,dummy);
//   flowEvent.SetReferenceMultiplicity( 64 );
//   for(Int_t i=0;i<64&&binx>0;i++){
//     AliFlowTrack *flowTrack=flowEvent.GetTrack(i);
//     Double_t inte=((TH2D*)fParHist->At(iParHist))->Integral(binx,binx,i+1,i+1);
//     if(inte>0)flowTrack->SetWeight(flowTrack->Weight()/inte);
//   }
//   if(fDebug>15)printf("EPfromV0 flow tracks weights done\n");

//   AliFlowVector qvec=flowEvent.GetQ(fV0EPorder);
//   Double_t angleEP=(1./(Double_t)fV0EPorder)*qvec.Phi();
//   if(fDebug>15)printf("EPfromV0 phi %f\n",angleEP);
//   return angleEP;
// }
//________________________________________________________________________
void AliAnalysisTaskSEHFv2::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
  if(fDebug > 1) printf("AnalysisTaskSEHFv2: Terminate() \n");
 
  return;
}

