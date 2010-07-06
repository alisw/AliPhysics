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

//*************************************************************************
// Class AliAnalysisTaskVertexESD
// AliAnalysisTask to extract from ESD the information for the analysis
// of the primary vertex reconstruction efficiency and resolution
// (for MC events) and distributions (for real data). Three vertices:
// - SPD tracklets
// - ITS+TPC tracks
// - TPC-only tracks
//
// Author: A.Dainese, andrea.dainese@pd.infn.it
//*************************************************************************

#include <TChain.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>  
#include <TCanvas.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliMultiplicity.h"
#include "AliVertexerTracks.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"
//#include "AliTriggerAnalysis.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"

#include "AliGenEventHeader.h" 
#include "AliAnalysisTaskVertexESD.h"


ClassImp(AliAnalysisTaskVertexESD)

//________________________________________________________________________
AliAnalysisTaskVertexESD::AliAnalysisTaskVertexESD(const char *name) : 
AliAnalysisTaskSE(name), 
fCheckEventType(kTRUE),
fReadMC(kFALSE),
fRecoVtxTPC(kFALSE),
fRecoVtxITSTPC(kFALSE),
fRecoVtxITSTPCHalfEvent(kFALSE),
fOnlyITSTPCTracks(kFALSE),
fOnlyITSSATracks(kFALSE),
fFillNtuple(kFALSE),
fFillTreeBeamSpot(kFALSE),
fESD(0), 
fOutput(0), 
fNtupleVertexESD(0),
fhSPDVertexX(0),
fhSPDVertexY(0),
fhSPDVertexZ(0),
fhTRKVertexX(0),
fhTRKVertexY(0),
fhTRKVertexZ(0),
fhTPCVertexX(0),
fhTPCVertexY(0),
fhTPCVertexZ(0),
fhTrackRefs(0),
fTreeBeamSpot(0)
{
  // Constructor

  // Define input and output slots here
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());  //My private output
}
//________________________________________________________________________
AliAnalysisTaskVertexESD::~AliAnalysisTaskVertexESD()
{
  // Destructor

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
}


//________________________________________________________________________
void AliAnalysisTaskVertexESD::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();

  fNtupleVertexESD = new TNtuple("fNtupleVertexESD","vertices","run:tstamp:bunchcross:triggered:dndygen:xtrue:ytrue:ztrue:xSPD:xerrSPD:ySPD:yerrSPD:zSPD:zerrSPD:ntrksSPD:SPD3D:dphiSPD:xTPC:xerrTPC:yTPC:yerrTPC:zTPC:zerrTPC:ntrksTPC:constrTPC:xTRK:xerrTRK:yTRK:yerrTRK:zTRK:zerrTRK:ntrksTRK:constrTRK:ntrklets:nESDtracks:nITSrefit5or6:nTPCin:nTPCinEta09:SPD0cls:xSPDp:xerrSPDp:ySPDp:yerrSPDp:zSPDp:zerrSPDp:ntrksSPDp:xTPCnc:xerrTPCnc:yTPCnc:yerrTPCnc:zTPCnc:zerrTPCnc:ntrksTPCnc:xTPCc:xerrTPCc:yTPCc:yerrTPCc:zTPCc:zerrTPCc:ntrksTPCc:xTRKnc:xerrTRKnc:yTRKnc:yerrTRKnc:zTRKnc:zerrTRKnc:ntrksTRKnc:xTRKc:xerrTRKc:yTRKc:yerrTRKc:zTRKc:zerrTRKc:ntrksTRKc:deltaxTRKnc:deltayTRKnc:deltazTRKnc:deltaxerrTRKnc:deltayerrTRKnc:deltazerrTRKnc:ntrksEvenTRKnc:ntrksOddTRKnc");

  fOutput->Add(fNtupleVertexESD);

  fhSPDVertexX = new TH1F("fhSPDVertexX","SPDVertex x; x vertex [cm]; events",200,-1,1);
  fOutput->Add(fhSPDVertexX);
  fhSPDVertexY = new TH1F("fhSPDVertexY","SPDVertex y; y vertex [cm]; events",200,-1,1);
  fOutput->Add(fhSPDVertexY);
  fhSPDVertexZ = new TH1F("fhSPDVertexZ","SPDVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(fhSPDVertexZ);
  fhTRKVertexX = new TH1F("fhTRKVertexX","TRKVertex x; x vertex [cm]; events",200,-1,1);
  fOutput->Add(fhTRKVertexX);
  fhTRKVertexY = new TH1F("fhTRKVertexY","TRKVertex y; y vertex [cm]; events",200,-1,1);
  fOutput->Add(fhTRKVertexY);
  fhTRKVertexZ = new TH1F("fhTRKVertexZ","TRKVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(fhTRKVertexZ);
  fhTPCVertexX = new TH1F("fhTPCVertexX","TPCVertex x; x vertex [cm]; events",200,-3,3);
  fOutput->Add(fhTPCVertexX);
  fhTPCVertexY = new TH1F("fhTPCVertexY","TPCVertex y; y vertex [cm]; events",200,-3,3);
  fOutput->Add(fhTPCVertexY);
  fhTPCVertexZ = new TH1F("fhTPCVertexZ","TPCVertex z; z vertex [cm]; events",200,-20,20);
  fOutput->Add(fhTPCVertexZ);

  fhTrackRefs = new TH2F("fhTrackRefs","Track references; x; y",1000,-4,4,1000,-4,4);
  fOutput->Add(fhTrackRefs);

  fTreeBeamSpot = new TTree("fTreeBeamSpot", "BeamSpotTree");
  UShort_t triggered, ntrkletsS, ntrksTRKnc;
  UInt_t run, bx;
  Float_t cetTimeLHC,xTRKnc, yTRKnc, zTRKnc;
  fTreeBeamSpot->Branch("run", &run, "run/i");
  fTreeBeamSpot->Branch("cetTimeLHC", &cetTimeLHC, "cetTimeLHC/F");
  fTreeBeamSpot->Branch("bx", &bx, "bx/i");
  fTreeBeamSpot->Branch("triggered", &triggered, "triggered/s");
  fTreeBeamSpot->Branch("ntrkletsS", &ntrkletsS, "ntrkletsS/s");
  fTreeBeamSpot->Branch("xTRKnc", &xTRKnc, "xTRKnc/F");
  fTreeBeamSpot->Branch("yTRKnc", &yTRKnc, "yTRKnc/F");
  fTreeBeamSpot->Branch("zTRKnc", &zTRKnc, "zTRKnc/F");
  fTreeBeamSpot->Branch("ntrksTRKnc", &ntrksTRKnc, "ntrksTRKnc/s");
  fOutput->Add(fTreeBeamSpot);

  PostData(1, fOutput);

  return;
}

//________________________________________________________________________
void AliAnalysisTaskVertexESD::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!InputEvent()) {
    Printf("ERROR: fESD not available");
    return;
  }
  AliESDEvent* esdE = (AliESDEvent*) InputEvent();
  
  // Select PHYSICS events (type=7, for data) and MC events (type=0)
  // fCheckEventType is kFALSE if fReadMC is kTRUE, hence check is skipped
  if(fCheckEventType) {
    if(esdE->GetEventType()!=7 && esdE->GetEventType()!=0) return; 
  }


  TArrayF mcVertex(3);
  mcVertex[0]=9999.; mcVertex[1]=9999.; mcVertex[2]=9999.;
  Float_t dNchdy=-999.;
  // ***********  MC info ***************
  if (fReadMC) {
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      return;
    }
    
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      return;
    }
    
    AliStack* stack = mcEvent->Stack();
    if (!stack) {
      AliDebug(AliLog::kError, "Stack not available");
      return;
    }
    
    AliHeader* header = mcEvent->Header();
    if (!header) {
      AliDebug(AliLog::kError, "Header not available");
      return;
    }
    AliGenEventHeader* genHeader = header->GenEventHeader();
    genHeader->PrimaryVertex(mcVertex);

    Int_t ngenpart = (Int_t)stack->GetNtrack();
    //printf("# generated particles = %d\n",ngenpart);
    dNchdy=0;
    for(Int_t ip=0; ip<ngenpart; ip++) {
      TParticle* part = (TParticle*)stack->Particle(ip);
      // keep only electorns, muons, pions, kaons and protons
      Int_t apdg = TMath::Abs(part->GetPdgCode());
      if(apdg!=11 && apdg!=13 && apdg!=211 && apdg!=321 && apdg!=2212) continue;      
      // reject secondaries
      if(TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]))>0.0010) continue;
      // reject incoming protons
      Double_t energy  = part->Energy();
      if(energy>900.) continue;
      Double_t pz = part->Pz();
      Double_t y = 0.5*TMath::Log((energy+pz+1.e-13)/(energy-pz+1.e-13));
      if(TMath::Abs(y)<1.0) dNchdy += 0.5; // count 1/2 of particles in |y|<1
      // tracks refs
      TClonesArray *trefs=0;
      Int_t ntrefs = mcEvent->GetParticleAndTR(ip,part,trefs);
      if(ntrefs<0) continue;
      for(Int_t iref=0; iref<ntrefs; iref++) {
	AliTrackReference *tref = (AliTrackReference*)trefs->At(iref);
	if(tref->R()>10.) continue;
	fhTrackRefs->Fill(tref->X(),tref->Y());
      }
    }
    //printf("# primary particles = %7.1f\n",dNchdy);
  } 
  // ***********  MC info ***************

    
  // Trigger
  //ULong64_t triggerMask;
  //ULong64_t spdFO = (1 << 14);
  //ULong64_t v0left = (1 << 10);
  //ULong64_t v0right = (1 << 11);
  
  //triggerMask=esdE->GetTriggerMask();
  // MB1: SPDFO || V0L || V0R
  //Bool_t eventTriggered = (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right))); 
  //MB2: GFO && V0R
  //triggerMask & spdFO && ((triggerMask&v0left) || (triggerMask&v0right))
  //Bool_t eventTriggered = (triggerMask & spdFO); 
 
  //static AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis();
  Bool_t eventTriggered = 0;//triggerAnalysis->IsTriggerFired(esdE, AliTriggerAnalysis::kSPDGFO /*| AliTriggerAnalysis::kOfflineFlag*/); 

  // use response of AliPhysicsSelection
  eventTriggered = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();


  Int_t nInFile = esdE->GetEventNumberInFile();

 const AliMultiplicity *alimult = esdE->GetMultiplicity();
  Int_t ntrklets=0,spd0cls=0;
  if(alimult) {
    ntrklets = alimult->GetNumberOfTracklets();

    for(Int_t l=0;l<alimult->GetNumberOfTracklets();l++){
      if(alimult->GetDeltaPhi(l)<-9998.) ntrklets--;
    }
    spd0cls = alimult->GetNumberOfSingleClusters()+ntrklets;
  }


  UShort_t triggered, ntrkletsS, ntrksTRKnc;
  UInt_t run, bx;
  Float_t cetTimeLHC,xTRKnc, yTRKnc, zTRKnc;
  fTreeBeamSpot->SetBranchAddress("run", &run);
  fTreeBeamSpot->SetBranchAddress("cetTimeLHC", &cetTimeLHC);
  fTreeBeamSpot->SetBranchAddress("bx", &bx);
  fTreeBeamSpot->SetBranchAddress("triggered", &triggered);
  fTreeBeamSpot->SetBranchAddress("ntrkletsS", &ntrkletsS);
  fTreeBeamSpot->SetBranchAddress("xTRKnc", &xTRKnc);
  fTreeBeamSpot->SetBranchAddress("yTRKnc", &yTRKnc);
  fTreeBeamSpot->SetBranchAddress("zTRKnc", &zTRKnc);
  fTreeBeamSpot->SetBranchAddress("ntrksTRKnc", &ntrksTRKnc);


  Double_t tstamp = esdE->GetTimeStamp();
  Float_t cetTime =(tstamp-1262304000.)+7200.;

  Float_t cetTime1h =(tstamp-1262304000.)+3600.;

  cetTimeLHC = (Float_t)cetTime1h;

  Int_t ntracks = esdE->GetNumberOfTracks();
  Int_t nITS5or6=0,nTPCin=0,nTPCinEta09=0;
  //printf("Tracks # = %d\n",esdE->GetNumberOfTracks());
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *t = esdE->GetTrack(itr);
    if(t->GetNcls(0)>=5) nITS5or6++;
    Double_t z0; t->GetZAt(0,esdE->GetMagneticField(),z0);
    if(t->GetNcls(1)>0 && TMath::Abs(t->GetD(0,0,esdE->GetMagneticField()))<2.8 && TMath::Abs(z0)<20) {
      nTPCin++;
      if(TMath::Abs(t->GetTgl())<1.5) nTPCinEta09++;
    }
  }


  const AliESDVertex *spdv=esdE->GetPrimaryVertexSPD();
  const AliESDVertex *spdvp=esdE->GetPileupVertexSPD(0);
  const AliESDVertex *tpcv=esdE->GetPrimaryVertexTPC();
  const AliESDVertex *trkv=esdE->GetPrimaryVertexTracks();
  
  // fill histos
  
  if(spdv) {
    if(spdv->GetNContributors()>0) {
      TString title=spdv->GetTitle();
      if(title.Contains("3D")) {
	fhSPDVertexX->Fill(spdv->GetXv());
	fhSPDVertexY->Fill(spdv->GetYv());
      }
      fhSPDVertexZ->Fill(spdv->GetZv());
    }
  }
  
  if(trkv) {
    if(trkv->GetNContributors()>0) {
      fhTRKVertexX->Fill(trkv->GetXv());
      fhTRKVertexY->Fill(trkv->GetYv());
      fhTRKVertexZ->Fill(trkv->GetZv());
    }
  }
  
  if(tpcv) {
    if(tpcv->GetNContributors()>0) {
      fhTPCVertexX->Fill(tpcv->GetXv());
      fhTPCVertexY->Fill(tpcv->GetYv());
      fhTPCVertexZ->Fill(tpcv->GetZv());
    }
  } 

  Float_t xpile=-999.;
  Float_t ypile=-999.;
  Float_t zpile=-999.;
  Float_t expile=-999.;
  Float_t eypile=-999.;
  Float_t ezpile=-999.;
  Int_t ntrkspile=-1;

  if(esdE->GetNumberOfPileupVerticesSPD()>0 && spdvp){
    xpile=spdvp->GetXv();
    expile=spdvp->GetXRes();
    ypile=spdvp->GetYv();
    eypile=spdvp->GetYRes();
    zpile=spdvp->GetZv();
    ezpile=spdvp->GetZRes();
    ntrkspile=spdvp->GetNContributors();  
  }

  // fill ntuple
  Int_t isize=82;
  Float_t xnt[82]; for(Int_t iii=0;iii<isize;iii++) xnt[iii]=0.;

  Int_t isizeSecNt=9;
  Float_t secnt[9]; for(Int_t iii=0;iii<isizeSecNt;iii++) secnt[iii]=0.;

  Int_t index=0;
  Int_t indexSecNt=0;

  xnt[index++]=(Float_t)esdE->GetRunNumber();
  secnt[indexSecNt++]=(Float_t)esdE->GetRunNumber();
  run = (Int_t)esdE->GetRunNumber();

  xnt[index++]=cetTime; //(Float_t)esdE->GetTimeStamp();
  //secnt[indexSecNt++]=cetTime;
  secnt[indexSecNt++]=cetTime1h;

  xnt[index++]=(Float_t)esdE->GetBunchCrossNumber();
  secnt[indexSecNt++]=(Float_t)esdE->GetBunchCrossNumber();
  bx = (Int_t)esdE->GetBunchCrossNumber();

  xnt[index++]=(eventTriggered ? 1. : 0.);
  secnt[indexSecNt++]=(eventTriggered ? 1. : 0.);
  triggered = (UShort_t)(eventTriggered ? 1 : 0);

  xnt[index++]=(Float_t)dNchdy;
  
  xnt[index++]=mcVertex[0];
  xnt[index++]=mcVertex[1];
  xnt[index++]=mcVertex[2];
  
  xnt[index++]=spdv->GetXv();
  xnt[index++]=spdv->GetXRes();
  xnt[index++]=spdv->GetYv();
  xnt[index++]=spdv->GetYRes();
  xnt[index++]=spdv->GetZv();
  xnt[index++]=spdv->GetZRes();
  xnt[index++]=spdv->GetNContributors();
  TString spdtitle = spdv->GetTitle();
  xnt[index++]=(spdtitle.Contains("vertexer: 3D") ? 1. : 0.);
  xnt[index++]=spdv->GetDispersion();
  
  xnt[index++]=tpcv->GetXv();
  xnt[index++]=tpcv->GetXRes();
  xnt[index++]=tpcv->GetYv();
  xnt[index++]=tpcv->GetYRes();
  xnt[index++]=tpcv->GetZv();
  xnt[index++]=tpcv->GetZRes();
  xnt[index++]=tpcv->GetNContributors();
  TString tpctitle = tpcv->GetTitle();
  xnt[index++]=(tpctitle.Contains("WithConstraint") ? 1. : 0.);
  
  xnt[index++]=trkv->GetXv();
  xnt[index++]=trkv->GetXRes();
  xnt[index++]=trkv->GetYv();
  xnt[index++]=trkv->GetYRes();
  xnt[index++]=trkv->GetZv();
  xnt[index++]=trkv->GetZRes();
  xnt[index++]=trkv->GetNContributors();// tpccontrorig;
  TString trktitle = trkv->GetTitle();
  xnt[index++]=(trktitle.Contains("WithConstraint") ? 1. : 0.);  

  xnt[index++]=float(ntrklets);
  secnt[indexSecNt++]=float(ntrklets);
  ntrkletsS = (UShort_t)ntrklets;

  xnt[index++]=float(ntracks);
  xnt[index++]=float(nITS5or6);

  xnt[index++]=float(nTPCin);
  xnt[index++]=float(nTPCinEta09);

  xnt[index++]=spd0cls;
    
  xnt[index++]=xpile;
  xnt[index++]=expile;
  xnt[index++]=ypile;
  xnt[index++]=eypile;
  xnt[index++]=zpile;
  xnt[index++]=ezpile;
  xnt[index++]=(Float_t)ntrkspile;  
    
    
 
  // add recalculated vertices TRK and TPC

  if(fRecoVtxTPC) {
    AliESDVertex *tpcvnc = ReconstructPrimaryVertexTPC(kFALSE);
    xnt[index++]=tpcvnc->GetXv();
    xnt[index++]=tpcvnc->GetXRes();
    xnt[index++]=tpcvnc->GetYv();
    xnt[index++]=tpcvnc->GetYRes();
    xnt[index++]=tpcvnc->GetZv();
    xnt[index++]=tpcvnc->GetZRes();
    xnt[index++]=tpcvnc->GetNContributors();
    delete tpcvnc; tpcvnc=0;
    
    AliESDVertex *tpcvc = ReconstructPrimaryVertexTPC(kTRUE);
    xnt[index++]=tpcvc->GetXv();
    xnt[index++]=tpcvc->GetXRes();
    xnt[index++]=tpcvc->GetYv();
    xnt[index++]=tpcvc->GetYRes();
    xnt[index++]=tpcvc->GetZv();
    xnt[index++]=tpcvc->GetZRes();
    xnt[index++]=tpcvc->GetNContributors();
    delete tpcvc; tpcvc=0;
  } else index+=14;

  if(fRecoVtxITSTPC) {
    AliESDVertex *trkvnc = ReconstructPrimaryVertexITSTPC(kFALSE);
    xnt[index++]=trkvnc->GetXv();
    xnt[index++]=trkvnc->GetXRes();
    xnt[index++]=trkvnc->GetYv();
    xnt[index++]=trkvnc->GetYRes();
    xnt[index++]=trkvnc->GetZv();
    xnt[index++]=trkvnc->GetZRes();
    xnt[index++]=trkvnc->GetNContributors();
  
    secnt[indexSecNt++]=trkvnc->GetXv();
    secnt[indexSecNt++]=trkvnc->GetYv();
    secnt[indexSecNt++]=trkvnc->GetZv();
    secnt[indexSecNt++]=trkvnc->GetNContributors();

    xTRKnc = (Float_t)trkvnc->GetXv();
    yTRKnc = (Float_t)trkvnc->GetYv();
    zTRKnc = (Float_t)trkvnc->GetZv();
    ntrksTRKnc = (trkvnc->GetNContributors()<0 ? 0 : (UShort_t)trkvnc->GetNContributors());

    delete trkvnc; trkvnc=0;


    AliESDVertex *trkvc = ReconstructPrimaryVertexITSTPC(kTRUE);
    xnt[index++]=trkvc->GetXv();
    xnt[index++]=trkvc->GetXRes();
    xnt[index++]=trkvc->GetYv();
    xnt[index++]=trkvc->GetYRes();
    xnt[index++]=trkvc->GetZv();
    xnt[index++]=trkvc->GetZRes();
    xnt[index++]=trkvc->GetNContributors();
    delete trkvc; trkvc=0;
  } else index+=14;
  
  if(fRecoVtxITSTPCHalfEvent) {
    AliESDVertex *trkvncodd  = ReconstructPrimaryVertexITSTPC(kFALSE,1);
    AliESDVertex *trkvnceven = ReconstructPrimaryVertexITSTPC(kFALSE,2);
    if(trkvncodd->GetNContributors()>0 && 
       trkvnceven->GetNContributors()>0) {
      xnt[index++]=trkvncodd->GetXv()-trkvnceven->GetXv();
      xnt[index++]=trkvncodd->GetYv()-trkvnceven->GetYv();
      xnt[index++]=trkvncodd->GetZv()-trkvnceven->GetZv();
      xnt[index++]=TMath::Sqrt(trkvncodd->GetXRes()*trkvncodd->GetXRes()+trkvnceven->GetXRes()*trkvnceven->GetXRes());
      xnt[index++]=TMath::Sqrt(trkvncodd->GetYRes()*trkvncodd->GetYRes()+trkvnceven->GetYRes()*trkvnceven->GetYRes());
      xnt[index++]=TMath::Sqrt(trkvncodd->GetZRes()*trkvncodd->GetZRes()+trkvnceven->GetZRes()*trkvnceven->GetZRes());
      xnt[index++]=trkvnceven->GetNContributors();
      xnt[index++]=trkvncodd->GetNContributors();
    } else {
      xnt[index++]=0.;
      xnt[index++]=0.;
      xnt[index++]=0.;
      xnt[index++]=0.;
      xnt[index++]=0.;
      xnt[index++]=0.;
      xnt[index++]=-1;
      xnt[index++]=-1;
    }
    delete trkvncodd; trkvncodd=0;
    delete trkvnceven; trkvnceven=0;    
  } else index+=8;
  

  if(index>isize) printf("AliAnalysisTaskVertexESD: ERROR, index!=isize\n");
  if(fFillNtuple) fNtupleVertexESD->Fill(xnt);

  if(indexSecNt>isizeSecNt) printf("AliAnalysisTaskVertexESD: ERROR, indexSecNt!=isizeSecNt\n");

  // if(indexTree>isizeTree) printf("AliAnalysisTaskVertexESD: ERROR, indexTree!=isizeTree\n");
  // only every second event (to reduce output size)
  if(fFillTreeBeamSpot && (nInFile%2 == 0)) fTreeBeamSpot->Fill();

  
  // Post the data already here
  PostData(1, fOutput);

  return;
}      

//________________________________________________________________________
void AliAnalysisTaskVertexESD::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }
  
  if (!fNtupleVertexESD){
    Printf("ERROR: fNtuple not available");
    return;
  }
  
  fNtupleVertexESD = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleVertexESD"));


  return;
}

//_________________________________________________________________________
AliESDVertex* AliAnalysisTaskVertexESD::ReconstructPrimaryVertexTPC(Bool_t constr) const {
  // On the fly reco of TPC vertex from ESD
  AliESDEvent* evt = (AliESDEvent*) fInputEvent;
  AliVertexerTracks vertexer(evt->GetMagneticField());
  vertexer.SetTPCMode(); // defaults
  Float_t diamondcovxy[3]; evt->GetDiamondCovXY(diamondcovxy);
  Double_t pos[3]={evt->GetDiamondX(),evt->GetDiamondY(),0}; 
  Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
  AliESDVertex *initVertex = new AliESDVertex(pos,cov,1.,1);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  if(!constr) vertexer.SetConstraintOff();

  return vertexer.FindPrimaryVertex(evt);
}

//_________________________________________________________________________
AliESDVertex* AliAnalysisTaskVertexESD::ReconstructPrimaryVertexITSTPC(Bool_t constr,Int_t mode) const {
  // On the fly reco of ITS+TPC vertex from ESD
  // mode=0 use all tracks
  // mode=1 use odd-number tracks
  // mode=2 use even-number tracks

  AliESDEvent* evt = (AliESDEvent*) fInputEvent;
  AliVertexerTracks vertexer(evt->GetMagneticField());
  vertexer.SetITSMode(); // defaults
  vertexer.SetMinClusters(4); // default is 5
  //vertexer.SetITSpureSA();
  Float_t diamondcovxy[3]; evt->GetDiamondCovXY(diamondcovxy);
  Double_t pos[3]={evt->GetDiamondX(),evt->GetDiamondY(),0}; 
  Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
  AliESDVertex *initVertex = new AliESDVertex(pos,cov,1.,1);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  if(!constr) vertexer.SetConstraintOff();

  // use only ITS-TPC or only ITS-SA tracks
  if(fOnlyITSTPCTracks || fOnlyITSSATracks || mode>0) {
    Int_t iskip=0;
    Int_t *skip = new Int_t[evt->GetNumberOfTracks()];
    for(Int_t itr=0;itr<evt->GetNumberOfTracks(); itr++) {
      AliESDtrack* track = evt->GetTrack(itr);
      if(fOnlyITSTPCTracks && track->GetNcls(1)==0) { // skip ITSSA
	skip[iskip++]=itr;
	continue;
      }
      if(fOnlyITSSATracks && track->GetNcls(1)>0) { // skip ITSTPC
	skip[iskip++]=itr;
	continue;
      }
      if(mode==1 && itr%2==0) skip[iskip++]=itr;
      if(mode==2 && itr%2==1) skip[iskip++]=itr;
    }
    vertexer.SetSkipTracks(iskip,skip);
    delete [] skip; skip=NULL;
  }

  return vertexer.FindPrimaryVertex(evt);
}
