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
#include "AliD0Trigger.h"
#endif
//-------------------------------------
// field (T)
const Double_t kBz = 0.4;

// primary vertex
Double_t primaryvertex[3] = {0.,0.,0,};

//sec. vertex
double v2[3]={0,0,0};

// sigle track cuts
const Double_t kPtCut = 0.5;  // GeV/c
const Double_t kd0Cut = 50.; // micron
const Double_t kd0CutHigh = 200.; // micron


//cuts for combined tracks
const Double_t cuts[7] = {0.005,     // cuts[0] = lowest V0 cut  (cm)
			  0.015,     // cuts[1] = highest V0 cut (cm)
			  0.05,      // cuts[2] = inv. mass cut (diferense) (Gev/c)
			  0.95,      // cuts[3] = max cosine value for pointing angle
			  -5000,     // cuts[4] = d0d0
			  0.8,       // cuts[5] = costhetastar
			  0.5};      // cuts[6] = ptchild
//cut for distance of closest aprach
double cutDCA=0.01;

// this function applies single track cuts
Bool_t TrkCuts(const AliITStrackV2& trk);

// this function creates TObjArrays with positive and negative tracks
void   SelectTracks(TTree& itsTree,
                      TObjArray& trksP,Int_t* itsEntryP,Int_t& nTrksP,
                      TObjArray& trksN,Int_t* itsEntryN,Int_t& nTrksN);

//void GetPrimaryVertex(int i,Char_t* path="./");

Int_t iTrkP,iTrkN,itsEntries;
Char_t trksName[100];
Int_t nTrksP=0,nTrksN=0;
Int_t nD0=0;
int ev=0;
double mom[6];

void RunD0offline(Int_t evFirst=0,Int_t evLast=1,Char_t* path="./") {

  const Char_t *name="AliD0offline";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name); 

  AliKalmanTrack::SetConvConst(100/0.299792458/kBz);

  // Open file with ITS tracks
  Char_t fnameTrack[1024];
  sprintf(fnameTrack,"%s/AliITStracksV2.root",path);
  TFile* itstrks = TFile::Open(fnameTrack);

   // tracks from ITS
  sprintf(trksName,"TreeT_ITS_%d",ev);
  TTree *itsTree=(TTree*)itstrks->Get(trksName);
  if(!itsTree) continue;
  itsEntries = (Int_t)itsTree->GetEntries();
  printf("+++\n+++ Number of tracks in ITS: %d\n+++\n\n",itsEntries);

  //Getting primary Vertex
  GetPrimaryVertex(0,path);

  // call function which applies sigle track selection and
  // separetes positives and negatives
  TObjArray trksP(itsEntries/2);
  Int_t *itsEntryP = new Int_t[itsEntries];
  TObjArray trksN(itsEntries/2);
  Int_t *itsEntryN = new Int_t[itsEntries];
  SelectTracks(*itsTree,trksP,itsEntryP,nTrksP,trksN,itsEntryN,nTrksN); 

  cout<<"#pos: "<<nTrksP<<endl;
  cout<<"#neg: "<<nTrksN<<endl;

  //the offline stuff
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

  AliD0Trigger * D0 = new AliD0Trigger(cuts,kBz,primaryvertex);

  double ptP,alphaP,phiP,ptN,alphaN,phiN,dca;

  for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
    postrack = (AliITStrackV2*)trksP.At(iTrkP);
    for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {
      negtrack = (AliITStrackV2*)trksN.At(iTrkN);
      D0.SetTracks(postrack,negtrack);
      //
      // ----------- DCA MINIMIZATION ------------------
      //
      // find the DCA and propagate the tracks to the DCA 
      double dca = vertexer2->PropagateToDCA(negtrack,postrack);
  
      if(dca<cutDCA){
	// define the AliV0vertex object
	AliV0vertex *vertex2 = new AliV0vertex(*negtrack,*postrack);
	// get position of the vertex
	vertex2->GetXYZ(v2[0],v2[1],v2[2]);
	delete vertex2;	
	if(D0.pTchild()){
	  if(D0.d0d0()){
	    if(D0.FindV0offline(v2)){
	      
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
	      
	      D0.SetMomenta(mom);
	      
	      if(D0.FindInvMass()){
		if(D0.CosThetaStar()){
		  if(D0.PointingAngle()){
		    nD0++;
		  }
		}
	      }
	    } 
	  }
	}
      }
    }
  } 
  cout<<"#D0: "<<nD0<<endl;
  gBenchmark->Stop(name);
  gBenchmark->Show(name);
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
  if(TMath::Abs(10000.*trk.GetD(primaryvertex[0],primaryvertex[1])) < kd0Cut) return kFALSE;
  if(TMath::Abs(10000.*trk.GetD(primaryvertex[0],primaryvertex[1])) > kd0CutHigh) return kFALSE;

  return kTRUE;
}
//____________________________________________________________________________
void GetPrimaryVertex(int i,Char_t* path="./") {

  int event=i;

  Char_t falice[1024];
  sprintf(falice,"%s/galice.root",path);
  TFile * galice = new TFile(falice);
  
  TDirectory * curdir;  

  Char_t vname[20];
  galice->cd();
  
  sprintf(vname,"Vertex_%d",event);
  TArrayF o = 0;
  o.Set(3);
  AliHeader * header = 0;
  
  TTree * treeE = (TTree*)gDirectory->Get("TE");
  treeE->SetBranchAddress("Header",&header);
  treeE->GetEntry(event);
  AliGenEventHeader* genHeader = header->GenEventHeader();
  if(genHeader){
    genHeader->PrimaryVertex(o);
    primaryvertex[0] = (Double_t)o[0];
    primaryvertex[1] = (Double_t)o[1];
    primaryvertex[2] = (Double_t)o[2];
  }
  else{
    printf("Can't find Header");
  }
  delete header;
  delete galice;
}
