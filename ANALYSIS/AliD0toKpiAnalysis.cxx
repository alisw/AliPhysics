/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

//----------------------------------------------------------------------------
//    Implementation of the D0toKpi reconstruction and analysis class
// Note: the two decay tracks are labelled: 0 (positive track)
//                                          1 (negative track)
// An example of usage can be found in the macro AliD0toKpiTest.C
//            Origin: A. Dainese    andrea.dainese@pd.infn.it            
//----------------------------------------------------------------------------
#include <TKey.h>
#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TSystem.h>
#include <TParticle.h>
#include "AliESD.h"
#include "AliStack.h"
#include "AliRunLoader.h"
#include "AliITSVertexerTracks.h"
#include "AliESDVertex.h"
#include "AliESDv0.h"
#include "AliD0toKpi.h"
#include "AliD0toKpiAnalysis.h"
#include "AliLog.h"

typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumlab;
  Int_t mumpdg;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} REFTRACK;

ClassImp(AliD0toKpiAnalysis)

//----------------------------------------------------------------------------
AliD0toKpiAnalysis::AliD0toKpiAnalysis() {
  // Default constructor
  SetPtCut();
  Setd0Cut();
  SetMassCut();
  SetD0Cuts();
  SetVertex1();
  SetPID();
  fVertexOnTheFly = kFALSE;
  fSim = kFALSE;
  fOnlySignal = kFALSE;
}
//----------------------------------------------------------------------------
AliD0toKpiAnalysis::~AliD0toKpiAnalysis() {}
//----------------------------------------------------------------------------
void AliD0toKpiAnalysis::ApplySelection(const Char_t *inName,const Char_t *outName) const {
  // select candidates that pass fD0Cuts and write them to a new file

  TFile *inFile = TFile::Open(inName);

  TTree *treeD0in=(TTree*)inFile->Get("TreeD0");
  AliD0toKpiAnalysis *inAnalysis = (AliD0toKpiAnalysis*)inFile->Get("D0toKpiAnalysis");
  printf("+++\n+++  I N P U T   S T A T U S:\n+++\n");
  inAnalysis->PrintStatus();


  AliD0toKpi *d = 0; 
  treeD0in->SetBranchAddress("D0toKpi",&d);
  Int_t entries = (Int_t)treeD0in->GetEntries();

  printf("+++\n+++ Number of D0 in input tree:  %d\n+++\n",entries);

  TTree *treeD0out = new TTree("TreeD0","Tree with selected D0 candidates");
  treeD0out->Branch("D0toKpi","AliD0toKpi",&d,200000,0);


  Int_t okD0=0,okD0bar=0;
  Int_t nSel = 0;

  for(Int_t i=0; i<entries; i++) {
    // get event from tree
    treeD0in->GetEvent(i);

    if(fSim && fOnlySignal && !d->IsSignal()) continue; 

    // check if candidate passes selection (as D0 or D0bar)
    if(d->Select(fD0Cuts,okD0,okD0bar)) {
      nSel++;
      treeD0out->Fill();
    }

  }

  AliD0toKpiAnalysis *outAnalysis = (AliD0toKpiAnalysis*)inAnalysis->Clone("D0toKpiAnalysis");
  outAnalysis->SetD0Cuts(fD0Cuts);
  printf("------------------------------------------\n");
  printf("+++\n+++  O U T P U T   S T A T U S:\n+++\n");
  outAnalysis->PrintStatus();

  printf("+++\n+++ Number of D0 in output tree:  %d\n+++\n",nSel);

  TFile* outFile = new TFile(outName,"recreate");
  treeD0out->Write();
  outAnalysis->Write();
  outFile->Close();

  return;
}
//----------------------------------------------------------------------------
Double_t AliD0toKpiAnalysis::CalculateTOFmass(Double_t mom,Double_t length,
					      Double_t time) const {
  // calculated the mass from momentum, track length from vertex to TOF
  // and time measured by the TOF
  if(length==0.) return -1000.;
  Double_t a = time*time/length/length;
  if(a > 1.) {
    a = TMath::Sqrt(a-1.);
  } else {
    a = -TMath::Sqrt(1.-a);
  }

  return mom*a;
}
/*
//----------------------------------------------------------------------------
void AliD0toKpiAnalysis::FindCandidates(Int_t evFirst,Int_t evLast,
					const Char_t *outName) {
  // Find D0 candidates and calculate parameters

  TString trkName("AliITStracksV2.root");
  if(gSystem->AccessPathName(trkName.Data(),kFileExists)) {
    printf("AliD0toKpiAnalysis::FindCandidates(): No tracks file found!\n"); 
    return;
  }

  TString outName1=outName;
  TString outName2("nTotEvents.dat");

  Int_t    ev;
  Int_t    nTotEv=0,nD0rec=0,nD0rec1ev=0;
  Double_t dca;
  Double_t v2[3],mom[6],d0[2];
  Double_t alphaP,alphaN,ptP,ptN,phiP,phiN;
  Int_t    iTrkP,iTrkN,trkEntries;
  Int_t    nTrksP=0,nTrksN=0;
  Int_t    trkNum[2];
  Int_t    okD0=0,okD0bar=0;
  Char_t   trkTreeName[100],vtxName[100];
  AliITStrackV2 *postrack = 0;
  AliITStrackV2 *negtrack = 0;

  // create the AliITSVertexerTracks object
  // (it will be used only if fVertexOnTheFly=kTrue)
  AliITSVertexerTracks *vertexer1 = new AliITSVertexerTracks;
  vertexer1->SetMinTracks(2);
  Int_t  skipped[2];
  Bool_t goodVtx1;
  

  // define the cuts for vertexing
  Double_t vtxcuts[]={50., // max. allowed chi2
		      0.0, // min. allowed negative daughter's impact param 
		      0.0, // min. allowed positive daughter's impact param 
		      1.0, // max. allowed DCA between the daughter tracks
        	     -1.0, // min. allowed cosine of pointing angle
		      0.0, // min. radius of the fiducial volume
		      2.9};// max. radius of the fiducial volume
  
  // create the AliV0vertexer object
  AliV0vertexer *vertexer2 = new AliV0vertexer(vtxcuts);

  // create tree for reconstructed D0s
  AliD0toKpi *ioD0toKpi=0;
  TTree *treeD0 = new TTree("TreeD0","Tree with D0 candidates");
  treeD0->Branch("D0toKpi","AliD0toKpi",&ioD0toKpi,200000,0);

  // open file with tracks
  TFile *trkFile = TFile::Open(trkName.Data());

  // loop on events in file
  for(ev=evFirst; ev<=evLast; ev++) {
    printf(" --- Looking for D0->Kpi in event  %d ---\n",ev);

    // retrieve primary vertex from file
    sprintf(vtxName,"VertexTracks_%d",ev);
    AliESDVertex *vertex1stored = (AliESDVertex*)trkFile->Get(vtxName);
    vertex1stored->GetXYZ(fV1);
    delete vertex1stored;

    // retrieve tracks from file
    sprintf(trkTreeName,"TreeT_ITS_%d",ev);
    TTree *trkTree=(TTree*)trkFile->Get(trkTreeName);
    if(!trkTree) { 
      printf("AliD0toKpiAnalysis::FindCandidates():\n  tracks tree not found for evet %d\n",ev); 
      continue; 
    }
    trkEntries = (Int_t)trkTree->GetEntries();
    printf(" Number of tracks: %d\n",trkEntries);

    // count the total number of events
    nTotEv++;

    // call function which applies sigle-track selection and
    // separetes positives and negatives
    TObjArray trksP(trkEntries/2);
    Int_t *trkEntryP = new Int_t[trkEntries];
    TObjArray trksN(trkEntries/2);
    Int_t *trkEntryN = new Int_t[trkEntries];
    SelectTracks(*trkTree,trksP,trkEntryP,nTrksP,trksN,trkEntryN,nTrksN);
      
    nD0rec1ev = 0;

    // loop on positive tracks
    for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
      if(iTrkP%10==0) printf("  Processing positive track number %d of %d\n",iTrkP,nTrksP);
	  
      // get track from track array
      postrack = (AliITStrackV2*)trksP.At(iTrkP);
      trkNum[0] = trkEntryP[iTrkP];      

      // loop on negative tracks 
      for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {
	// get track from tracks array
	negtrack = (AliITStrackV2*)trksN.At(iTrkN);
	trkNum[1] = trkEntryN[iTrkN];      

	AliITStrackV2 nt(*negtrack), pt(*postrack), *pnt=&nt, *ppt=&pt;

	//
	// ----------- DCA MINIMIZATION ------------------
	//
	// find the DCA and propagate the tracks to the DCA 
	dca = vertexer2->PropagateToDCA(pnt,ppt);

	// define the AliV0vertex object
	AliV0vertex *vertex2 = new AliV0vertex(*pnt,*ppt);
	  
	// get position of the secondary vertex
	vertex2->GetXYZ(v2[0],v2[1],v2[2]);
	
	delete vertex2;
  
	// momenta of the tracks at the vertex
	ptP = 1./TMath::Abs(ppt->Get1Pt());
	alphaP = ppt->GetAlpha();
	phiP = alphaP+TMath::ASin(ppt->GetSnp());
	mom[0] = ptP*TMath::Cos(phiP); 
	mom[1] = ptP*TMath::Sin(phiP);
	mom[2] = ptP*ppt->GetTgl();
	  
	ptN = 1./TMath::Abs(pnt->Get1Pt());
	alphaN = pnt->GetAlpha();
	phiN = alphaN+TMath::ASin(pnt->GetSnp());
	mom[3] = ptN*TMath::Cos(phiN); 
	mom[4] = ptN*TMath::Sin(phiN);
	mom[5] = ptN*pnt->GetTgl();
	  
	goodVtx1 = kTRUE;
	// no vertexing if DeltaMass > fMassCut 
	if(fVertexOnTheFly) {
	  goodVtx1 = kFALSE;
	  if(SelectInvMass(mom)) {
	    // primary vertex from *other* tracks in event
	    vertexer1->SetVtxStart(fV1[0],fV1[1]);
	    skipped[0] = trkEntryP[iTrkP];
	    skipped[1] = trkEntryN[iTrkN];
	    vertexer1->SetSkipTracks(2,skipped);
	    AliESDVertex *vertex1onfly = 
	      (AliESDVertex*)vertexer1->VertexOnTheFly(*trkTree); 
	    if(vertex1onfly->GetNContributors()>0) goodVtx1 = kTRUE;
	    vertex1onfly->GetXYZ(fV1);
	    //vertex1onfly->PrintStatus();
	    delete vertex1onfly;	
	  }
	}	  

	// impact parameters of the tracks w.r.t. the primary vertex
	d0[0] =  10000.*ppt->GetD(fV1[0],fV1[1]);
	d0[1] = -10000.*pnt->GetD(fV1[0],fV1[1]);

	// create the object AliD0toKpi
	AliD0toKpi theD0(ev,trkNum,fV1,v2,dca,mom,d0);

	// select D0s
	if(goodVtx1 && theD0.Select(fD0Cuts,okD0,okD0bar)) {
	      
	  // fill the tree
	  ioD0toKpi=&theD0;
	  treeD0->Fill();

	  nD0rec++; nD0rec1ev++;
	  ioD0toKpi=0;  
	}
	
	negtrack = 0;
      } // loop on negative tracks
      postrack = 0;
    }   // loop on positive tracks

    trksP.Delete();
    trksN.Delete();
    delete [] trkEntryP;
    delete [] trkEntryN;
    delete trkTree;

    printf(" Number of D0 candidates: %d\n",nD0rec1ev);
  }    // loop on events in file


  printf("\n+++\n+++ Total number of events: %d\n+++\n",nTotEv);
  printf("\n+++\n+++ Total number of D0 candidates: %d\n+++\n",nD0rec);

  delete vertexer1;
  delete vertexer2;

  trkFile->Close();

  // create a copy of this class to be written to output file
  AliD0toKpiAnalysis *copy = (AliD0toKpiAnalysis*)this->Clone("D0toKpiAnalysis");
  copy->PrintStatus();

  // add PDG codes to decay tracks in found candidates (in simulation mode)
  // and store tree in the output file
  if(!fSim) {
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeD0->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  } else {
    printf(" Now adding information from simulation (PDG codes) ...\n");
    TTree *treeD0sim = new TTree("TreeD0","Tree with D0 candidates");
    MakeTracksRefFile(evFirst,evLast);
    SimulationInfo(treeD0,treeD0sim);
    delete treeD0;
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeD0sim->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  }

  // write to a file the total number of events
  FILE *outdat = fopen(outName2.Data(),"w");
  fprintf(outdat,"%d\n",nTotEv);
  fclose(outdat);

  return;
}
*/
//----------------------------------------------------------------------------
void AliD0toKpiAnalysis::FindCandidatesESD(Int_t evFirst,Int_t evLast,
					   const Char_t *outName) {
  // Find D0 candidates and calculate parameters

  TString esdName("AliESDs.root");
  if(gSystem->AccessPathName(esdName.Data(),kFileExists)) {
    printf("AliD0toKpiAnalysis::FindCandidatesESD(): No ESDs file found!\n"); 
    return;
  }

  TString outName1=outName;
  TString outName2("nTotEvents.dat");

  Int_t    nTotEv=0,nD0rec=0,nD0rec1ev=0;
  Double_t dca;
  Double_t v2[3],mom[6],d0[2];
  //Double_t alphaP,alphaN,ptP,ptN,phiP,phiN;
  Int_t    iTrkP,iTrkN,trkEntries;
  Int_t    nTrksP=0,nTrksN=0;
  Int_t    trkNum[2];
  Double_t tofmass[2];
  Int_t    okD0=0,okD0bar=0;
  AliESDtrack *postrack = 0;
  AliESDtrack *negtrack = 0;

  // create the AliITSVertexerTracks object
  // (it will be used only if fVertexOnTheFly=kTrue)
  AliITSVertexerTracks *vertexer1 = new AliITSVertexerTracks;
  vertexer1->SetMinTracks(2);
  Int_t  skipped[2];
  Bool_t goodVtx1;
  
  /*
  // define the cuts for vertexing
  Double_t vtxcuts[]={50., // max. allowed chi2
		      0.0, // min. allowed negative daughter's impact param 
		      0.0, // min. allowed positive daughter's impact param 
		      1.0, // max. allowed DCA between the daughter tracks
        	     -1.0, // min. allowed cosine of pointing angle
		      0.0, // min. radius of the fiducial volume
		      2.9};// max. radius of the fiducial volume
  
  // create the AliV0vertexer object
  AliV0vertexer *vertexer2 = new AliV0vertexer(vtxcuts);
  */

  // create tree for reconstructed D0s
  AliD0toKpi *ioD0toKpi=0;
  TTree *treeD0 = new TTree("TreeD0","Tree with D0 candidates");
  treeD0->Branch("D0toKpi","AliD0toKpi",&ioD0toKpi,200000,0);

  // open file with tracks
  TFile *esdFile = TFile::Open(esdName.Data());
  AliESD* event = new AliESD;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    Error("FindCandidatesESD", "no ESD tree found");
    return;
  }
  tree->SetBranchAddress("ESD", &event);

  // loop on events in file
  for (Int_t iEvent = evFirst; iEvent < tree->GetEntries(); iEvent++) {
    if (iEvent > evLast) break;
    tree->GetEvent(iEvent);
    Int_t ev = (Int_t)event->GetEventNumber();
    printf("--- Finding D0 -> Kpi in event %d\n",ev);

    // retrieve primary vertex from file
    //sprintf(vtxName,"Vertex_%d",ev);
    //AliESDVertex *vertex1stored = (AliESDVertex*)trkFile->Get(vtxName);
    //vertex1stored->GetXYZ(fV1);
    //delete vertex1stored;
    event->GetVertex()->GetXYZ(fV1);

    trkEntries = event->GetNumberOfTracks();
    printf(" Number of tracks: %d\n",trkEntries);

    // count the total number of events
    nTotEv++;

    // call function which applies sigle-track selection and
    // separetes positives and negatives
    TObjArray trksP(trkEntries/2);
    Int_t    *trkEntryP   = new Int_t[trkEntries];
    TObjArray trksN(trkEntries/2);
    Int_t    *trkEntryN   = new Int_t[trkEntries];
    TTree *trkTree = new TTree();
    if(fVertexOnTheFly) {
      SelectTracksESDvtx(*event,trkTree,trksP,trkEntryP,nTrksP,
	        	                trksN,trkEntryN,nTrksN);
    } else {
      SelectTracksESD(*event,trksP,trkEntryP,nTrksP,
	        	     trksN,trkEntryN,nTrksN);
    }      

    AliDebugClass(1,Form(" pos. tracks: %d    neg .tracks: %d",nTrksP,nTrksN));

    nD0rec1ev = 0;

    // loop on positive tracks
    for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
      if(iTrkP%10==0) AliDebugClass(1,Form("  Processing positive track number %d of %d",iTrkP,nTrksP));
	  
      // get track from track array
      postrack = (AliESDtrack*)trksP.UncheckedAt(iTrkP);
      trkNum[0] = trkEntryP[iTrkP];      

      // loop on negative tracks 
      for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {

	// get track from tracks array
	negtrack = (AliESDtrack*)trksN.UncheckedAt(iTrkN);
	trkNum[1] = trkEntryN[iTrkN];      

        {
	//
	// ----------- DCA MINIMIZATION ------------------
	//
	// find the DCA and propagate the tracks to the DCA
	Double_t b=event->GetMagneticField(); 
	AliESDtrack nt(*negtrack), pt(*postrack);
	dca = nt.PropagateToDCA(&pt,b);

	// define the AliESDv0 object
	AliESDv0 vertex2(nt,trkNum[0],pt,trkNum[1]);
	  
	// get position of the secondary vertex
	vertex2.GetXYZ(v2[0],v2[1],v2[2]);
        vertex2.GetPPxPyPz(mom[0],mom[1],mom[2]);
        vertex2.GetNPxPyPz(mom[3],mom[4],mom[5]);
	// impact parameters of the tracks w.r.t. the primary vertex

	d0[0] =  10000.*pt.GetD(fV1[0],fV1[1],b);
	d0[1] = -10000.*nt.GetD(fV1[0],fV1[1],b);
	}
        /*
	// momenta of the tracks at the vertex
        //Double_t x,par[5]; postrack->GetExternalParameters(x,par); 
	//ptP = 1./TMath::Abs(par[4]);
	//alphaP = postrack->GetAlpha();
	//phiP = alphaP+TMath::ASin(par[2]);
	  postrack->GetPxPyPz();
	mom[0] = ptP*TMath::Cos(phiP); 
	mom[1] = ptP*TMath::Sin(phiP);
	mom[2] = ptP*par[3];
	  
	ptN = 1./TMath::Abs(pnt->Get1Pt());
	alphaN = pnt->GetAlpha();
	phiN = alphaN+TMath::ASin(pnt->GetSnp());
	mom[3] = ptN*TMath::Cos(phiN); 
	mom[4] = ptN*TMath::Sin(phiN);
	mom[5] = ptN*pnt->GetTgl();
	*/
	goodVtx1 = kTRUE;
	// no vertexing if DeltaMass > fMassCut 
	if(fVertexOnTheFly) {
	  goodVtx1 = kFALSE;
	  if(SelectInvMass(mom)) {
	    // primary vertex from *other* tracks in event
	    vertexer1->SetVtxStart(fV1[0],fV1[1]);
	    skipped[0] = trkEntryP[iTrkP];
	    skipped[1] = trkEntryN[iTrkN];
	    vertexer1->SetSkipTracks(2,skipped);
	    AliESDVertex *vertex1onfly = 
	      (AliESDVertex*)vertexer1->VertexOnTheFly(*trkTree); 
	    if(vertex1onfly->GetNContributors()>0) goodVtx1 = kTRUE;
	    vertex1onfly->GetXYZ(fV1);
	    //vertex1onfly->PrintStatus();
	    delete vertex1onfly;	
	  }
	}	  

	// create the object AliD0toKpi
	AliD0toKpi theD0(ev,trkNum,fV1,v2,dca,mom,d0);

	// select D0s
	if(goodVtx1 && theD0.Select(fD0Cuts,okD0,okD0bar)) {
	      
	  // get PID info from ESD
	  AliESDtrack *t0 = (AliESDtrack*)event->GetTrack(trkNum[0]);
	  Double_t esdpid0[5];
	  t0->GetESDpid(esdpid0);
	  if(t0->GetStatus()&AliESDtrack::kTOFpid) {
	    tofmass[0] = CalculateTOFmass(t0->GetP(),
					  t0->GetIntegratedLength(),
					  t0->GetTOFsignal());
	  } else {
	    tofmass[0] = -1000.;
	  }
	  AliESDtrack *t1 = (AliESDtrack*)event->GetTrack(trkNum[1]);
	  Double_t esdpid1[5];
	  t1->GetESDpid(esdpid1);
	  if(t1->GetStatus()&AliESDtrack::kTOFpid) {
	    tofmass[1] = CalculateTOFmass(t1->GetP(),
					  t1->GetIntegratedLength(),
					  t1->GetTOFsignal());
	  } else {
	    tofmass[1] = -1000.;
	  }

	  theD0.SetPIDresponse(esdpid0,esdpid1);
	  theD0.SetTOFmasses(tofmass);

	  // fill the tree
	  ioD0toKpi=&theD0;
	  treeD0->Fill();

	  nD0rec++; nD0rec1ev++;
	  ioD0toKpi=0;  
	}
	
	negtrack = 0;
      } // loop on negative tracks
      postrack = 0;
    }   // loop on positive tracks

    trksP.Delete();
    trksN.Delete();
    delete [] trkEntryP;
    delete [] trkEntryN;
    delete trkTree;


    printf(" Number of D0 candidates: %d\n",nD0rec1ev);
  }    // loop on events in file


  printf("\n+++\n+++ Total number of events: %d\n+++\n",nTotEv);
  printf("\n+++\n+++ Total number of D0 candidates: %d\n+++\n",nD0rec);

  delete vertexer1;
  //delete vertexer2;

  esdFile->Close();

  // create a copy of this class to be written to output file
  AliD0toKpiAnalysis *copy = (AliD0toKpiAnalysis*)this->Clone("D0toKpiAnalysis");

  // add PDG codes to decay tracks in found candidates (in simulation mode)
  // and store tree in the output file
  if(!fSim) {
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeD0->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  } else {
    printf(" Now adding information from simulation (PDG codes) ...\n");
    TTree *treeD0sim = new TTree("TreeD0","Tree with D0 candidates");
    MakeTracksRefFileESD();
    SimulationInfo(treeD0,treeD0sim);
    delete treeD0;
    TFile *outroot = new TFile(outName1.Data(),"recreate");
    treeD0sim->Write();
    copy->Write();
    outroot->Close();
    delete outroot;
  }

  // write to a file the total number of events
  FILE *outdat = fopen(outName2.Data(),"w");
  fprintf(outdat,"%d\n",nTotEv);
  fclose(outdat);

  return;
}
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::PrintStatus() const {
  // Print parameters being used

  printf("Preselections:\n");
  printf("    fPtCut   = %f GeV\n",fPtCut);
  printf("    fd0Cut   = %f micron\n",fd0Cut);
  printf("    fMassCut = %f GeV\n",fMassCut);
  if(fVertexOnTheFly) printf("Primary vertex on the fly\n");
  if(fSim) { 
    printf("Simulation mode\n");
    if(fOnlySignal) printf("  Only signal goes to file\n");
  }
  printf("Cuts on candidates:\n");
  printf("    |M-MD0| [GeV]    < %f\n",fD0Cuts[0]);
  printf("    dca    [micron]  < %f\n",fD0Cuts[1]);
  printf("    cosThetaStar     < %f\n",fD0Cuts[2]);
  printf("    pTK     [GeV]    > %f\n",fD0Cuts[3]);
  printf("    pTpi    [GeV]    > %f\n",fD0Cuts[4]);
  printf("    |d0K|  [micron]  < %f\n",fD0Cuts[5]);
  printf("    |d0pi| [micron]  < %f\n",fD0Cuts[6]);
  printf("    d0d0  [micron^2] < %f\n",fD0Cuts[7]);
  printf("    cosThetaPoint    > %f\n",fD0Cuts[8]);

  return;
}
//-----------------------------------------------------------------------------
Bool_t AliD0toKpiAnalysis::SelectInvMass(const Double_t p[6]) const {
  // Apply preselection in the invariant mass of the pair

  Double_t mD0 = 1.8645;
  Double_t mPi = 0.13957;
  Double_t mKa = 0.49368;

  Double_t energy[2];
  Double_t mom2[2],momTot2;

  mom2[0] = p[0]*p[0] + p[1]*p[1] + p[2]*p[2];
  mom2[1] = p[3]*p[3] + p[4]*p[4] + p[5]*p[5];

  momTot2 = (p[0]+p[3])*(p[0]+p[3])+
            (p[1]+p[4])*(p[1]+p[4])+
            (p[2]+p[5])*(p[2]+p[5]);

  // D0 -> K- pi+
  energy[1] = TMath::Sqrt(mKa*mKa+mom2[1]);
  energy[0] = TMath::Sqrt(mPi*mPi+mom2[0]);

  Double_t minvD0 = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);
    
  // D0bar -> K+ pi-
  energy[0] = TMath::Sqrt(mKa*mKa+mom2[0]);
  energy[1] = TMath::Sqrt(mPi*mPi+mom2[1]);

  Double_t minvD0bar = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);

  if(TMath::Abs(minvD0-mD0)    < fMassCut) return kTRUE;
  if(TMath::Abs(minvD0bar-mD0) < fMassCut) return kTRUE;
  return kFALSE;
}
/*
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SelectTracks(TTree &trkTree,
		  TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
		  TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const {
  // Create two TObjArrays with positive and negative tracks and 
  // apply single-track preselection

  nTrksP=0,nTrksN=0;
 
  Int_t entr = (Int_t)trkTree.GetEntries();

  // trasfer tracks from tree to arrays
  for(Int_t i=0; i<entr; i++) {

    AliITStrackV2 *track = new AliITStrackV2; 
    trkTree.SetBranchAddress("tracks",&track);

    trkTree.GetEvent(i);

    // single track selection
    if(!SingleTrkCuts(*track)) { delete track; continue; }

    if(track->Get1Pt()>0.) { // negative track
      trksN.AddLast(track);
      trkEntryN[nTrksN] = i;
      nTrksN++;
    } else {                 // positive track
      trksP.AddLast(track);
      trkEntryP[nTrksP] = i;
      nTrksP++;
    }

  }

  return;
}
*/
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SelectTracksESD(AliESD &event,
        TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
        TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const {
  // Create two TObjArrays with positive and negative tracks and 
  // apply single-track preselection

  nTrksP=0,nTrksN=0;

  Int_t entr = event.GetNumberOfTracks();
 
  // transfer ITS tracks from ESD to arrays and to a tree
  for(Int_t i=0; i<entr; i++) {

    AliESDtrack *esdtrack = event.GetTrack(i);
    UInt_t status = esdtrack->GetStatus();

    if(!(status&AliESDtrack::kITSrefit)) continue;

    // single track selection
    if(!SingleTrkCuts(*esdtrack,event.GetMagneticField())) continue;

    if(esdtrack->GetSign()<0) { // negative track
      trksN.AddLast(esdtrack);
      trkEntryN[nTrksN] = i;
      nTrksN++;
    } else {                 // positive track
      trksP.AddLast(esdtrack);
      trkEntryP[nTrksP] = i;
      nTrksP++;
    }

  } // loop on ESD tracks

  return;
}
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SelectTracksESDvtx(AliESD &event,TTree *trkTree,
        TObjArray &trksP,Int_t *trkEntryP,Int_t &nTrksP,
        TObjArray &trksN,Int_t *trkEntryN,Int_t &nTrksN) const {
  // Create two TObjArrays with positive and negative tracks and 
  // apply single-track preselection

  nTrksP=0,nTrksN=0;

  Int_t entr = event.GetNumberOfTracks();
 
  AliESDtrack *esdtrackfortree = 0;
  trkTree->Branch("tracks","AliESDtrack",&esdtrackfortree,entr,0);


  // transfer the tracks from ESD to arrays and to a tree
  for(Int_t i=0; i<entr; i++) {

    AliESDtrack *esdtrack = event.GetTrack(i);
    UInt_t status = esdtrack->GetStatus();

    if(!(status&AliESDtrack::kITSrefit)) continue;

    // store track in the tree to be used for primary vertex finding
    esdtrackfortree = new AliESDtrack(*esdtrack);
    trkTree->Fill();
    //itstrackfortree = 0;
    delete esdtrackfortree;

    // single track selection
    if(!SingleTrkCuts(*esdtrack,event.GetMagneticField())) continue;

    if(esdtrack->GetSign()<0) { // negative track
      trksN.AddLast(esdtrack);
      trkEntryN[nTrksN] = i;
      nTrksN++;
    } else {                 // positive track
      trksP.AddLast(esdtrack);
      trkEntryP[nTrksP] = i;
      nTrksP++;
    }

  } // loop on esd tracks

  //delete itstrackfortree;

  return;
}
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SetD0Cuts(Double_t cut0,Double_t cut1,
				   Double_t cut2,Double_t cut3,Double_t cut4,
				   Double_t cut5,Double_t cut6,
				   Double_t cut7,Double_t cut8) {
  // Set the cuts for D0 selection
  fD0Cuts[0] = cut0;
  fD0Cuts[1] = cut1;
  fD0Cuts[2] = cut2;
  fD0Cuts[3] = cut3;
  fD0Cuts[4] = cut4;
  fD0Cuts[5] = cut5;
  fD0Cuts[6] = cut6;
  fD0Cuts[7] = cut7;
  fD0Cuts[8] = cut8;

  return;
}
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SetD0Cuts(const Double_t cuts[9]) {
  // Set the cuts for D0 selection

  for(Int_t i=0; i<9; i++) fD0Cuts[i] = cuts[i];

  return;
}
//-----------------------------------------------------------------------------
Bool_t 
AliD0toKpiAnalysis::SingleTrkCuts(const AliESDtrack& trk, Double_t b) const {
  // Check if track passes some kinematical cuts  
  // Magnetic field "b" (kG)

  if(TMath::Abs(1./trk.GetParameter()[4]) < fPtCut) 
    return kFALSE;
  if(TMath::Abs(10000.*trk.GetD(fV1[0],fV1[1],b)) < fd0Cut) 
    return kFALSE;

  return kTRUE;
}
/*
//----------------------------------------------------------------------------
void AliD0toKpiAnalysis::MakeTracksRefFile(Int_t evFirst,Int_t evLast) 
  const {
  // Create a file with simulation info for the reconstructed tracks
  
  TFile *out = TFile::Open("ITStracksRefFile.root","recreate");
  TFile *trk = TFile::Open("AliITStracksV2.root");
  AliRunLoader *rl = AliRunLoader::Open("galice.root");

  // load kinematics
  rl->LoadKinematics();
  
  Int_t      label;
  TParticle *part;  
  TParticle *mumpart;
  REFTRACK   reftrk;
  
  for(Int_t ev=evFirst; ev<=evLast; ev++){
    rl->GetEvent(ev);  
    AliStack *stack = rl->Stack();

    trk->cd();

    // Tree with tracks
    char tname[100];
    sprintf(tname,"TreeT_ITS_%d",ev);

    TTree *tracktree=(TTree*)trk->Get(tname);
    if(!tracktree) continue;
    AliITStrackV2 *track = new AliITStrackV2; 
    tracktree->SetBranchAddress("tracks",&track);
    Int_t nentr=(Int_t)tracktree->GetEntries();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",ev);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&reftrk,"lab/I:pdg:mumlab:mumpdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for(Int_t i=0; i<nentr; i++) {
      tracktree->GetEvent(i);
      label = TMath::Abs(track->GetLabel());

      part = (TParticle*)stack->Particle(label);
      reftrk.lab = label;
      reftrk.pdg = part->GetPdgCode();
      reftrk.mumlab = part->GetFirstMother();
      if(part->GetFirstMother()>=0) {
	mumpart = (TParticle*)stack->Particle(part->GetFirstMother());
	reftrk.mumpdg = mumpart->GetPdgCode();
      } else {
	reftrk.mumpdg=-1;
      }
      reftrk.Vx = part->Vx();
      reftrk.Vy = part->Vy();
      reftrk.Vz = part->Vz();
      reftrk.Px = part->Px();
      reftrk.Py = part->Py();
      reftrk.Pz = part->Pz();
      
      reftree->Fill();
    } // loop on tracks   

    out->cd();
    reftree->Write();

    delete track;
    delete reftree;
    delete tracktree;
    delete stack;
  } // loop on events

  trk->Close();
  out->Close();

  return;
}
*/
//----------------------------------------------------------------------------
void AliD0toKpiAnalysis::MakeTracksRefFileESD() const {
  // Create a file with simulation info for the reconstructed tracks
  
  TFile *outFile = TFile::Open("ITStracksRefFile.root","recreate");
  TFile *esdFile = TFile::Open("AliESDs.root");
  AliRunLoader *rl = AliRunLoader::Open("galice.root");

  // load kinematics
  rl->LoadKinematics();
  
  Int_t      label;
  TParticle *part;  
  TParticle *mumpart;
  REFTRACK   reftrk;
  
  TKey *key=0;
  TIter next(esdFile->GetListOfKeys());
  // loop on events in file
  while ((key=(TKey*)next())!=0) {
    AliESD *event=(AliESD*)key->ReadObj();
    Int_t ev = (Int_t)event->GetEventNumber();

    rl->GetEvent(ev);  
    AliStack *stack = rl->Stack();

    Int_t nentr=(Int_t)event->GetNumberOfTracks();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",ev);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&reftrk,"lab/I:pdg:mumlab:mumpdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for(Int_t i=0; i<nentr; i++) {
      AliESDtrack *esdtrack = (AliESDtrack*)event->GetTrack(i);
      label = TMath::Abs(esdtrack->GetLabel());

      part = (TParticle*)stack->Particle(label);
      reftrk.lab = label;
      reftrk.pdg = part->GetPdgCode();
      reftrk.mumlab = part->GetFirstMother();
      if(part->GetFirstMother()>=0) {
	mumpart = (TParticle*)stack->Particle(part->GetFirstMother());
	reftrk.mumpdg = mumpart->GetPdgCode();
      } else {
	reftrk.mumpdg=-1;
      }
      reftrk.Vx = part->Vx();
      reftrk.Vy = part->Vy();
      reftrk.Vz = part->Vz();
      reftrk.Px = part->Px();
      reftrk.Py = part->Py();
      reftrk.Pz = part->Pz();
      
      reftree->Fill();

    } // loop on tracks   

    outFile->cd();
    reftree->Write();

    delete reftree;
    delete event;
    delete stack;
  } // loop on events

  esdFile->Close();
  outFile->Close();

  return;
}
//-----------------------------------------------------------------------------
void AliD0toKpiAnalysis::SimulationInfo(TTree *treeD0in,TTree *treeD0out) const {
  // add pdg codes to candidate decay tracks (for sim)

  TString refFileName("ITStracksRefFile.root");
  if(fSim && gSystem->AccessPathName(refFileName.Data(),kFileExists)) { 
    printf("AliD0toKpiAnalysis::SimulationInfo: no reference file found!\n"); 
    return;
  }
  TFile *refFile = TFile::Open(refFileName.Data());

  Char_t refTreeName[100];
  Int_t  event;
  Int_t  pdg[2],mumpdg[2],mumlab[2];
  REFTRACK reftrk;

  // read-in reference tree for event 0 (the only event for Pb-Pb)
  sprintf(refTreeName,"Tree_Ref_%d",0);
  TTree *refTree0 = (TTree*)refFile->Get(refTreeName);
  refTree0->SetBranchAddress("rectracks",&reftrk);

  AliD0toKpi *theD0 = 0; 
  treeD0in->SetBranchAddress("D0toKpi",&theD0);
  treeD0out->Branch("D0toKpi","AliD0toKpi",&theD0,200000,0);

  Int_t entries = (Int_t)treeD0in->GetEntries();

  for(Int_t i=0; i<entries; i++) {
    if(i%100==0) printf("  done %d candidates of %d\n",i,entries);    

    treeD0in->GetEvent(i);
    event = theD0->EventNo();

    if(event==0) { // always true for Pb-Pb (avoid to read-in tree every time)
      refTree0->GetEvent(theD0->GetTrkNum(0));
      pdg[0]    = reftrk.pdg;
      mumpdg[0] = reftrk.mumpdg;
      mumlab[0] = reftrk.mumlab;
      refTree0->GetEvent(theD0->GetTrkNum(1));
      pdg[1]    = reftrk.pdg;
      mumpdg[1] = reftrk.mumpdg;
      mumlab[1] = reftrk.mumlab;
    } else {
      sprintf(refTreeName,"Tree_Ref_%d",event);
      TTree *refTree = (TTree*)refFile->Get(refTreeName);
      refTree->SetBranchAddress("rectracks",&reftrk);
      refTree->GetEvent(theD0->GetTrkNum(0));
      pdg[0]    = reftrk.pdg;
      mumpdg[0] = reftrk.mumpdg;
      mumlab[0] = reftrk.mumlab;
      refTree->GetEvent(theD0->GetTrkNum(1));
      pdg[1]    = reftrk.pdg;
      mumpdg[1] = reftrk.mumpdg;
      mumlab[1] = reftrk.mumlab;
      delete refTree;
    }
    
    theD0->SetPdgCodes(pdg);
    theD0->SetMumPdgCodes(mumpdg);
    
    if(TMath::Abs(mumpdg[0])==421 && TMath::Abs(mumpdg[1])==421 
       && mumlab[0]==mumlab[1]) theD0->SetSignal();
    
    if(!fOnlySignal || theD0->IsSignal()) treeD0out->Fill();

  }

  delete refTree0;

  refFile->Close();

  return;
}






