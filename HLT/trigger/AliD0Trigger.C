/****************************************************************************
 *                                                                          *
 * This macro computes the position of the D0 decay vertex using            *
 * helix DCA minimization by J.Belikov                                      *
 *                                                                          *
 * Reconstructed D0 are stored in a tree as AliD0toKpi objects              *
 *                                                                          *
 ****************************************************************************/ 
//#define __COMPILE__
//#ifdef __COMPILE__
#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <iostream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TVector3.h>
#include <TTree.h>
#include <TParticle.h>
#include <TArray.h>
//----- AliRoot headers ---------------
#include "alles.h"
#include "AliRun.h"
#include "AliKalmanTrack.h"
#include "AliITStrackV2.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliV0vertex.h"
#include "AliV0vertexer.h"
#include "AliITSVertex.h"
#include "AliITSVertexer.h"
#include "AliITSVertexerTracks.h"
#include "AliD0toKpi.h"
#endif
//-------------------------------------

typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumlab;
  Int_t mumpdg;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} RECTRACK;

// field (T)
const Double_t kBz = 0.4;

// primary vertex
Double_t gv1[3] = {0.,0.,0,};

// sigle track cuts
const Double_t kPtCut = 0.5;  // GeV/c
const Double_t kd0Cut = 50.; // micron
//const Double_t kPtCut = 0.;  // GeV/c
//const Double_t kd0Cut = 0.; // micron

// cuts on D0 candidate (to be passed to function AliD0toKpi::Select())
//
// cuts[0] = inv. mass half width [GeV]   
// cuts[1] = dca [micron]
// cuts[2] = cosThetaStar 
// cuts[3] = pTK [GeV/c]
// cuts[4] = pTPi [GeV/c]
// cuts[5] = d0K [micron]   upper limit!
// cuts[6] = d0Pi [micron]  upper limit!
// cuts[7] = d0d0 [micron^2]
// cuts[8] = cosThetaPoint
const Double_t kD0cuts[9] = {0.012,
                             500.,
			     0.8,
			     0.5,
			     0.5,
			     500.,
			     500.,
			    -5000.,
			     0.8};
const Double_t kD0cutsAll[9] = {.1,
				1.e5,
				1.1,
				0.,
				0.,
				1.e10,
				1.e10,
				1.e10,
				-1.1};
// mass cut for 1st level rejection
const Double_t kMassCut = 0.1; // GeV

// this function ckecks files existence
Bool_t GetInputFiles(char* path[1024]);

// 1st level rejection based on inv. mass
Bool_t SelectInvMass(const Double_t p[6]);

// this function applies single track cuts
Bool_t TrkCuts(const AliITStrackV2& trk);

// this function creates TObjArrays with positive and negative tracks
void   SelectTracks(TTree& itsTree,
		      TObjArray& trksP,Int_t* itsEntryP,Int_t& nTrksP,
		      TObjArray& trksN,Int_t* itsEntryN,Int_t& nTrksN);


void   AliD0Trigger(Int_t evFirst=0,Int_t evLast=1,Char_t* path="./") {

  const Char_t *name="AliD0Trigger";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);  
  TStopwatch timer;
  AliKalmanTrack::SetConvConst(100/0.299792458/kBz);

  // load library for D0 candidates class
  //gSystem->Load("libD0toKpi.so");

  // check existence of input files 
  if(!GetInputFiles(path)) { cerr<<"No tracks file found"<<endl; return; }

  Char_t fname[1024];
  sprintf(fname,"%s/AliD0toKpiBkg.root",path);
  TString outName(fname);

  Int_t ev;
  Int_t nD0rec=0,nD0rec1ev=0;
  Double_t dca;
  Double_t v2[3],mom[6],d0[2];
  Int_t pdg[2],mum[2];
  Double_t alphaP,alphaN,ptP,ptN,phiP,phiN;
  Int_t iTrkP,iTrkN,itsEntries;
  Int_t nTrksP=0,nTrksN=0;
  Int_t okD0=0,okD0bar=0;
  Char_t trksName[100],refsName[100],vtxName[100];
  RECTRACK rectrk;  
  AliITStrackV2 *postrack = 0;
  AliITStrackV2 *negtrack = 0;

  // create the AliITSVertexerTracks object
  AliITSVertexerTracks *vertexer1 = new AliITSVertexerTracks;
  //vertexer1->SetUseThrustFrame(1);
  vertexer1->SetMinTracks(2);
  vertexer1->SetDebug(0);
  Int_t skipped[2];

  // define the cuts for vertexing
  Double_t vtxcuts[]={33., // max. allowed chi2
		      0.0, // min. allowed negative daughter's impact param 
		      0.0, // min. allowed positive daughter's impact param 
		      1.0, // max. allowed DCA between the daughter tracks
        	     -1.0, // min. allowed cosine of V0's pointing angle
		      0.0, // min. radius of the fiducial volume
		      2.9};// max. radius of the fiducial volume
  
  // create the AliV0vertexer object
  AliV0vertexer *vertexer2 = new AliV0vertexer(vtxcuts);


  // Create tree for reconstructed D0s
  TTree D0Tree("TreeD0","Tree with D0 candidates");
  AliD0toKpi *ioD0toKpi=0;
  D0Tree.Branch("D0toKpi","AliD0toKpi",&ioD0toKpi,200000,0);

  // Open file with ITS tracks
  Char_t fnameTrack[1024];
  sprintf(fnameTrack,"%s/AliITStracksV2.root",path);
  TFile* itstrks = TFile::Open(fnameTrack);

  // Open file with ITS track references
  Char_t fnameRef[1024];
  sprintf(fnameRef,"%s/ITStracksRefFile.root",path);
  TFile* itsrefs = TFile::Open(fnameRef);

  // Open Primary Vertex File
  Char_t fnameVertex[1024];
  sprintf(fnameVertex,"%s/Vtx4Tracking.root",path);
  TFile* vertexFile = TFile::Open(fnameVertex);

  // loop on events in file
  for(ev=evFirst; ev<=evLast; ev++) {
    printf(" --- Processing event  %d ---\n",ev);
    sprintf(vtxName,"VertexTracks_%d",ev);
    sprintf(trksName,"TreeT_ITS_%d",ev);
    sprintf(refsName,"Tree_Ref_%d",ev);


    // tracks from ITS
    TTree *itsTree=(TTree*)itstrks->Get(trksName);
    if(!itsTree) continue;
    itsEntries = (Int_t)itsTree->GetEntries();
    printf("+++\n+++ Number of tracks in ITS: %d\n+++\n\n",itsEntries);
    
    // tree from reference file
    TTree *refTree=(TTree*)itsrefs->Get(refsName);
    refTree->SetBranchAddress("rectracks",&rectrk);
    
    //get the primary vertex from headerfile

    AliITSVertex *vertex = vertexFile->Get("Vertex_0");
    vertex->GetXYZ(gv1);

    cout<<gv1[2]<<endl;
    
    // call function which applies sigle track selection and
    // separetes positives and negatives
    
    TObjArray trksP(itsEntries/2);
    Int_t *itsEntryP = new Int_t[itsEntries];
    TObjArray trksN(itsEntries/2);
    Int_t *itsEntryN = new Int_t[itsEntries];
    SelectTracks(*itsTree,trksP,itsEntryP,nTrksP,trksN,itsEntryN,nTrksN);

    nD0rec1ev = 0;

    cout<<"NnegTracks: "<<nTrksN<<endl;

    // loop on positive tracks
    for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
      //if(iTrkP % 10==0) cerr<<"--- Processing positive track number: "<<iTrkP<<" of "<<nTrksP<<"\r";
	  
      // get track from track array
      postrack = (AliITStrackV2*)trksP.At(iTrkP);

      // get info on tracks PDG and mothers PDG from reference file
      refTree->GetEvent(itsEntryP[iTrkP]);
      pdg[0] = rectrk.pdg;
      mum[0] = rectrk.mumpdg;

      // loop on negative tracks 
      for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {
	
	// get track from track array
	negtrack = (AliITStrackV2*)trksN.At(iTrkN);
      
	// get info on tracks PDG and mothers PDG from reference file
	refTree->GetEvent(itsEntryN[iTrkN]);
	pdg[1] = rectrk.pdg;
	mum[1] = rectrk.mumpdg;
 
	//
	// ----------- DCA MINIMIZATION ------------------
	//
	// find the DCA and propagate the tracks to the DCA 
	dca = vertexer2->PropagateToDCA(negtrack,postrack);

	// define the AliV0vertex object
	AliV0vertex *vertex2 = new AliV0vertex(*negtrack,*postrack);
	  
	// get position of the vertex
	vertex2->GetXYZ(v2[0],v2[1],v2[2]);
	
	delete vertex2;
  
	// momenta of the tracks at the vertex
	ptP = 1./TMath::Abs(postrack->Get1Pt());
	alphaP = postrack->GetAlpha();
	phiP = alphaP+TMath::ASin(postrack->GetSnp());
	mom[0] = ptP*TMath::Cos(phiP); 
	mom[1] = ptP*TMath::Sin(phiP);
	mom[2] = ptP*postrack->GetTgl();
	  
	ptN = 1./TMath::Abs(negtrack->Get1Pt());
	alphaN = negtrack->GetAlpha();
	phiN = alphaN+TMath::ASin(negtrack->GetSnp());
	mom[3] = ptN*TMath::Cos(phiN); 
	mom[4] = ptN*TMath::Sin(phiN);
	mom[5] = ptN*negtrack->GetTgl();
	  

	Bool_t goodVtx = kFALSE;
	// no vertexing if DeltaMass > kMassCut 
	if(SelectInvMass(mom)) {
	  // primary vertex from other tracks in event
	  vertexer1->SetVtxStart(0.,0.);
	  skipped[0] = itsEntryP[iTrkP];
	  skipped[1] = itsEntryN[iTrkN];
	  vertexer1->SetSkipTracks(2,skipped);
	  // impact parameters of the tracks w.r.t. the smeared primary vertex
	  d0[0] =  10000.*postrack->GetD(gv1[0],gv1[1]);
	  d0[1] = -10000.*negtrack->GetD(gv1[0],gv1[1]);
	}
	  
	// create the object TD0rec and store it in the tree
	AliD0toKpi theD0(kFALSE,ev,gv1,v2,dca,mom,d0,pdg,mum);

	
	// select D0s
	if(goodVtx && theD0.Select(kD0cuts,okD0,okD0bar)) {
	      
	  // compute the weights
	  theD0.ComputeWgts();

	  // fill the tree
	  ioD0toKpi=&theD0;
	  D0Tree.Fill();

	  nD0rec++; nD0rec1ev++;
	  ioD0toKpi=0;  
	}
	
	negtrack = 0;
      } // loop on negative tracks
      postrack = 0;
    }   // loop on positive tracks

    delete [] itsEntryP;
    delete [] itsEntryN;
    delete itsTree;
    delete refTree;

    printf("\n+++\n+++ Number of D0 candidates: %d\n+++\n",nD0rec1ev);
    printf("----------------------------------------------------------\n");

  }    // loop on events in file

  printf("\n+++\n+++ Total number of D0 candidates: %d\n+++\n",nD0rec);

  delete vertexer1;
  delete vertexer2;

  itstrks->Close();
  itsrefs->Close();

  // store tree in a file
  TFile* outroot = new TFile(outName.Data(),"recreate");
  D0Tree.Write();
  outroot->Close();
  delete outroot;

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return;
}
//___________________________________________________________________________
Bool_t   GetInputFiles(Char_t* path2="./") {
  Char_t fnameTrack[1024];
  sprintf(fnameTrack,"%s/AliITStracksV2.root",path2);
  TString itsName(fnameTrack);
  Char_t fnameRef[1024];
  sprintf(fnameRef,"%s/ITStracksRefFile.root",path2);
  TString refName(fnameRef);

  if(gSystem->AccessPathName(itsName.Data(),kFileExists) ||
     gSystem->AccessPathName(refName.Data(),kFileExists)) return kFALSE;

  return kTRUE;
}
//___________________________________________________________________________
Bool_t SelectInvMass(const Double_t p[6]) {

  Double_t mD0 = 1.8645;
  Double_t mPi = 0.13957;
  Double_t mKa = 0.49368;

  Double_t energy[2];
  Double_t mom2[2],momTot2;

  mom2[0] = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
  mom2[1] = p[3]*p[3]+p[4]*p[4]+p[5]*p[5];

  momTot2 = (p[0]+p[3])*(p[0]+p[3])+
            (p[1]+p[4])*(p[1]+p[4])+
            (p[2]+p[5])*(p[2]+p[5]);

  // D0 -> K- Pi+
  energy[1] = TMath::Sqrt(mKa*mKa+mom2[1]);
  energy[0] = TMath::Sqrt(mPi*mPi+mom2[0]);

  Double_t minvD0 = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);
    

  // D0bar -> K+ Pi-
  energy[0] = TMath::Sqrt(mKa*mKa+mom2[0]);
  energy[1] = TMath::Sqrt(mPi*mPi+mom2[1]);

  Double_t minvD0bar = TMath::Sqrt((energy[0]+energy[1])*(energy[0]+energy[1])-momTot2);

  if(TMath::Abs(minvD0-mD0) < kMassCut)    return kTRUE;
  if(TMath::Abs(minvD0bar-mD0) < kMassCut) return kTRUE;
  return kFALSE;
}
//___________________________________________________________________________
void   SelectTracks(TTree& itsTree,
		    TObjArray& trksP,Int_t* itsEntryP,Int_t& nTrksP,
		    TObjArray& trksN,Int_t* itsEntryN,Int_t& nTrksN) {
//
// this function creates two TObjArrays with positive and negative tracks
//
  nTrksP=0,nTrksN=0;

 
  Int_t entr = (Int_t)itsTree.GetEntries();

  // trasfer tracks from tree to arrays
  for(Int_t i=0; i<entr; i++) {

    AliITStrackV2 *itstrack = new AliITStrackV2; 
    itsTree.SetBranchAddress("tracks",&itstrack);

    itsTree.GetEvent(i);

    // single track selection
    if(!TrkCuts(*itstrack)) { delete itstrack; continue; }

    if(itstrack->Get1Pt()>0.) { // negative track
      trksN.AddLast(itstrack);
      itsEntryN[nTrksN] = i;
      nTrksN++;
    } else {                    // positive track
      trksP.AddLast(itstrack);
      itsEntryP[nTrksP] = i;
      nTrksP++;
    }

  }

  return;
}
//____________________________________________________________________________
Bool_t TrkCuts(const AliITStrackV2& trk) {
// 
// this function tells if track passes some kinematical cuts  
//
  if(TMath::Abs(1./trk.Get1Pt()) < kPtCut)                return kFALSE;
  if(TMath::Abs(10000.*trk.GetD(gv1[0],gv1[1])) < kd0Cut) return kFALSE;

  return kTRUE;
}








