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

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliLog.h"

#include "AliGenEventHeader.h" 
#include "AliAnalysisTaskVertexESD.h"


ClassImp(AliAnalysisTaskVertexESD)

//________________________________________________________________________
AliAnalysisTaskVertexESD::AliAnalysisTaskVertexESD(const char *name) : 
AliAnalysisTask(name, "VertexESDTask"), 
fReadMC(kFALSE),
fESD(0), 
fOutput(0), 
fNtupleVertexESD(0),
fhTrackRefs(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(0, TList::Class());  //My private output
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
void AliAnalysisTaskVertexESD::ConnectInputData(Option_t *) 
{
  // Connect ESD or AOD here
  // Called once

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    
    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else {
      fESD = esdH->GetEvent();
    }

  }
  
  return;
}

//________________________________________________________________________
void AliAnalysisTaskVertexESD::CreateOutputObjects()
{
  // Create histograms
  // Called once

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();

  fNtupleVertexESD = new TNtuple("fNtupleVertexESD","vertices","xtrue:ytrue:ztrue:xSPD:xdiffSPD:xerrSPD:ySPD:ydiffSPD:yerrSPD:zSPD:zdiffSPD:zerrSPD:ntrksSPD:xTPC:xdiffTPC:xerrTPC:yTPC:ydiffTPC:yerrTPC:zTPC:zdiffTPC:zerrTPC:ntrksTPC:xTRK:xdiffTRK:xerrTRK:yTRK:ydiffTRK:yerrTRK:zTRK:zdiffTRK:zerrTRK:ntrksTRK:ntrklets:nESDtracks:nITSrefit5or6:nTPCin:nTPCinEta09:dndygen:triggered:SPD3D:SPD0cls:constrTRK");

  fOutput->Add(fNtupleVertexESD);

  fhTrackRefs = new TH2F("fhTrackRefs","Track references; x; y",1000,-4,4,1000,-4,4);
  fOutput->Add(fhTrackRefs);

 return;
}

//________________________________________________________________________
void AliAnalysisTaskVertexESD::Exec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
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
  ULong64_t triggerMask;
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);
  
  triggerMask=fESD->GetTriggerMask();
  // MB1: SPDFO || V0L || V0R
  Bool_t eventTriggered = (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right))); 
  //MB2: GFO && V0R
  //triggerMask & spdFO && ((triggerMask&v0left) || (triggerMask&v0right))
 

  Int_t ntracks = fESD->GetNumberOfTracks();
  Int_t nITS5or6=0,nTPCin=0,nTPCinEta09=0;
  //printf("Tracks # = %d\n",fESD->GetNumberOfTracks());
  for(Int_t itr=0; itr<ntracks; itr++) {
    AliESDtrack *t = fESD->GetTrack(itr);
    if(t->GetNcls(0)>=5) nITS5or6++;
    Double_t z0; t->GetZAt(0,fESD->GetMagneticField(),z0);
    if(t->GetNcls(1)>0 && TMath::Abs(t->GetD(0,0,fESD->GetMagneticField()))<2.8 && TMath::Abs(z0)<20) {
      nTPCin++;
      if(TMath::Abs(t->GetTgl())<1.5) nTPCinEta09++;
    }
  }

    
  const AliESDVertex *spdv=fESD->GetPrimaryVertexSPD();
  const AliESDVertex *tpcv=fESD->GetPrimaryVertexTPC();
  const AliESDVertex *trkv=fESD->GetPrimaryVertex();

  //Float_t tpccontrorig=tpcv->GetNContributors();
  //tpcv = 0;
  //tpcv = ReconstructPrimaryVertexTPC();

  const AliMultiplicity *alimult = fESD->GetMultiplicity();
  Int_t ntrklets=0,spd0cls=0;
  if(alimult) {
    ntrklets = alimult->GetNumberOfTracklets();
    for(Int_t l=0;l<alimult->GetNumberOfTracklets();l++){
      if(alimult->GetDeltaPhi(l)<-9998.) ntrklets--;
    }
    spd0cls = alimult->GetNumberOfSingleClusters()+ntrklets;
  }


  Float_t xnt[43];
  
  xnt[0]=mcVertex[0];
  xnt[1]=mcVertex[1];
  xnt[2]=mcVertex[2];
  
  xnt[3]=spdv->GetXv();
  xnt[4]=(spdv->GetXv()-mcVertex[0])*10000.;
  xnt[5]=spdv->GetXRes();
  xnt[6]=spdv->GetYv();
  xnt[7]=(spdv->GetYv()-mcVertex[1])*10000.;
  xnt[8]=spdv->GetYRes();
  xnt[9]=spdv->GetZv();
  xnt[10]=(spdv->GetZv()-mcVertex[2])*10000.;
  xnt[11]=spdv->GetZRes();
  xnt[12]=spdv->GetNContributors();
  
  xnt[13]=tpcv->GetXv();
  xnt[14]=(tpcv->GetXv()-mcVertex[0])*10000.;
  xnt[15]=tpcv->GetXRes();
  xnt[16]=tpcv->GetYv();
  xnt[17]=(tpcv->GetYv()-mcVertex[1])*10000.;
  xnt[18]=tpcv->GetYRes();
  xnt[19]=tpcv->GetZv();
  xnt[20]=(tpcv->GetZv()-mcVertex[2])*10000.;
  xnt[21]=tpcv->GetZRes();
  xnt[22]=tpcv->GetNContributors();
  
  xnt[23]=trkv->GetXv();
  xnt[24]=(trkv->GetXv()-mcVertex[0])*10000.;
  xnt[25]=trkv->GetXRes();
  xnt[26]=trkv->GetYv();
  xnt[27]=(trkv->GetYv()-mcVertex[1])*10000.;
  xnt[28]=trkv->GetYRes();
  xnt[29]=trkv->GetZv();
  xnt[30]=(trkv->GetZv()-mcVertex[2])*10000.;
  xnt[31]=trkv->GetZRes();
  xnt[32]=trkv->GetNContributors();// tpccontrorig;
  

  xnt[33]=float(ntrklets);
  xnt[34]=float(ntracks);
  xnt[35]=float(nITS5or6);
  xnt[36]=float(nTPCin);
  xnt[37]=float(nTPCinEta09);

  xnt[38]=float(dNchdy);

  xnt[39]=0.; 
  if(eventTriggered) xnt[39]=1.;

  xnt[40]=0.; 
  TString spdtitle = spdv->GetTitle();
  if(spdtitle.Contains("vertexer: 3D")) xnt[40]=1.;

  xnt[42]=0.; 
  TString trktitle = trkv->GetTitle();
  if(trktitle.Contains("WithConstraint")) xnt[42]=1.;

  xnt[41]=spd0cls;

  fNtupleVertexESD->Fill(xnt);
    
  // Post the data already here
  PostData(0, fOutput);

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
AliESDVertex* AliAnalysisTaskVertexESD::ReconstructPrimaryVertexTPC() {
  // On the fly reco of TPC vertex from ESD

  AliVertexerTracks vertexer(fESD->GetMagneticField());
  vertexer.SetTPCMode(0.1,1.0,5.,0,1,3.,0.1,1.5); // defaults
  //vertexer.SetTPCMode(0.1,1.0,5.,0,1,3.,0.1,1.5); // defaults
  Double_t pos[3]={+0.0220,-0.0340,+0.270}; 
  Double_t err[3]={0.0200,0.0200,7.5};
  AliESDVertex *initVertex = new AliESDVertex(pos,err);
  vertexer.SetVtxStart(initVertex);
  delete initVertex;
  //vertexer.SetConstraintOff();

  return vertexer.FindPrimaryVertex(fESD);
}
