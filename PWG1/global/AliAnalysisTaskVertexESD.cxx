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
fOnlyITSTPCTracks(kFALSE),
fOnlyITSSATracks(kFALSE),
fFillNtuple(kFALSE),
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
fhTrackRefs(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
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

  fNtupleVertexESD = new TNtuple("fNtupleVertexESD","vertices","run:tstamp:xtrue:ytrue:ztrue:xSPD:xerrSPD:ySPD:yerrSPD:zSPD:zerrSPD:ntrksSPD:xTPC:xerrTPC:yTPC:yerrTPC:zTPC:zerrTPC:ntrksTPC:xTRK:xerrTRK:yTRK:yerrTRK:zTRK:zerrTRK:ntrksTRK:ntrklets:nESDtracks:nITSrefit5or6:nTPCin:nTPCinEta09:dndygen:triggered:SPD3D:SPD0cls:constrTRK:constrTPC");

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
  
  if(fCheckEventType && (esdE->GetEventType())!=7) return; 


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
  const AliESDVertex *tpcv=esdE->GetPrimaryVertexTPC();
  const AliESDVertex *trkv=esdE->GetPrimaryVertexTracks();

  //Float_t tpccontrorig=tpcv->GetNContributors();

  if(fRecoVtxTPC) {
    tpcv = 0;
    tpcv = ReconstructPrimaryVertexTPC();
  }
  if(fRecoVtxITSTPC) {
    trkv = 0;
    trkv = ReconstructPrimaryVertexITSTPC();
  }

  const AliMultiplicity *alimult = esdE->GetMultiplicity();
  Int_t ntrklets=0,spd0cls=0;
  if(alimult) {
    ntrklets = alimult->GetNumberOfTracklets();
    for(Int_t l=0;l<alimult->GetNumberOfTracklets();l++){
      if(alimult->GetDeltaPhi(l)<-9998.) ntrklets--;
    }
    spd0cls = alimult->GetNumberOfSingleClusters()+ntrklets;
  }
  
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
  

  // fill ntuple
  Int_t isize=37;
  Float_t xnt[37];
  
  Int_t index=0;

  xnt[index++]=(Float_t)esdE->GetRunNumber();
  xnt[index++]=(Float_t)esdE->GetTimeStamp();

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
  
  xnt[index++]=tpcv->GetXv();
  xnt[index++]=tpcv->GetXRes();
  xnt[index++]=tpcv->GetYv();
  xnt[index++]=tpcv->GetYRes();
  xnt[index++]=tpcv->GetZv();
  xnt[index++]=tpcv->GetZRes();
  xnt[index++]=tpcv->GetNContributors();
  
  xnt[index++]=trkv->GetXv();
  xnt[index++]=trkv->GetXRes();
  xnt[index++]=trkv->GetYv();
  xnt[index++]=trkv->GetYRes();
  xnt[index++]=trkv->GetZv();
  xnt[index++]=trkv->GetZRes();
  xnt[index++]=trkv->GetNContributors();// tpccontrorig;
  

  xnt[index++]=float(ntrklets);
  xnt[index++]=float(ntracks);
  xnt[index++]=float(nITS5or6);
  xnt[index++]=float(nTPCin);
  xnt[index++]=float(nTPCinEta09);

  xnt[index++]=float(dNchdy);

  xnt[index++]=(eventTriggered ? 1. : 0.);

  TString spdtitle = spdv->GetTitle();
  xnt[index++]=(spdtitle.Contains("vertexer: 3D") ? 1. : 0.);

  xnt[index++]=spd0cls;

  TString trktitle = trkv->GetTitle();
  xnt[index++]=(trktitle.Contains("WithConstraint") ? 1. : 0.);

  TString tpctitle = tpcv->GetTitle();
  xnt[index++]=(tpctitle.Contains("WithConstraint") ? 1. : 0.);


  if(index!=isize) printf("AliAnalysisTaskVertexESD: ERROR, index!=isize\n");

  if(fFillNtuple) fNtupleVertexESD->Fill(xnt);
  
  // Post the data already here
  PostData(1, fOutput);

  return;
}      

//________________________________________________________________________
void AliAnalysisTaskVertexESD::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(0));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }

  fNtupleVertexESD = dynamic_cast<TNtuple*>(fOutput->FindObject("fNtupleVertexESD"));

  return;
}

//_________________________________________________________________________
AliESDVertex* AliAnalysisTaskVertexESD::ReconstructPrimaryVertexTPC() const {
  // On the fly reco of TPC vertex from ESD
  AliESDEvent* evt = (AliESDEvent*) fInputEvent;
  AliVertexerTracks vertexer(evt->GetMagneticField());
  vertexer.SetTPCMode(); // defaults
  //vertexer.SetTPCMode(0.1,1.0,5.,0,1,3.,0.1,1.5);
  Double_t pos[3]={+0.0220,-0.0340,+0.270}; 
  Double_t err[3]={0.0200,0.0200,7.5};
  AliESDVertex *initVertex = new AliESDVertex(pos,err);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  vertexer.SetConstraintOff();

  return vertexer.FindPrimaryVertex(evt);
}

//_________________________________________________________________________
AliESDVertex* AliAnalysisTaskVertexESD::ReconstructPrimaryVertexITSTPC() const {
  // On the fly reco of ITS+TPC vertex from ESD
  AliESDEvent* evt = (AliESDEvent*) fInputEvent;
  AliVertexerTracks vertexer(evt->GetMagneticField());
  vertexer.SetITSMode(); // defaults
  //vertexer.SetTPCMode(0.1,1.0,5.,0,1,3.,0.1,1.5);
  Double_t pos[3]={+0.0220,-0.0340,+0.270}; 
  Double_t err[3]={0.0200,0.0200,7.5};
  AliESDVertex *initVertex = new AliESDVertex(pos,err);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  vertexer.SetConstraintOff();

  // use only ITS-TPC or only ITS-SA tracks
  if(fOnlyITSTPCTracks || fOnlyITSSATracks) {
    Int_t iskip=0;
    Int_t *skip = new Int_t[evt->GetNumberOfTracks()];
    for(Int_t itr=0;itr<evt->GetNumberOfTracks(); itr++) {
      AliESDtrack* track = evt->GetTrack(itr);
      if(fOnlyITSTPCTracks && track->GetNcls(1)==0) { // skip ITSSA
	skip[iskip++]=itr;
      }
      if(fOnlyITSSATracks && track->GetNcls(1)>0) { // skip ITSTPC
	skip[iskip++]=itr;
      }
    }
    vertexer.SetSkipTracks(iskip,skip);
    delete [] skip; skip=NULL;
  }

  return vertexer.FindPrimaryVertex(evt);
}
