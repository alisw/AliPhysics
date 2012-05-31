//#define __COMPILE__
//#ifdef __COMPILE__
#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include <iostream.h>
#include <fstream.h>
//--------Root headers ---------------
#include <TSystem.h>
#include <TFile.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TVector3.h>
#include <TTree.h>
#include <TObjArray>
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
double primaryvertex[3]={0.,0.,0.};

//sec. vertex
double v2[3]={0,0,0};

// sigle track cuts
const Double_t kPtCut = 0.5;        // 0.5 GeV/c
const Double_t kd0Cut = 50.;        // 50  micron
const Double_t kd0CutHigh = 400.;   // 200 micron


//cuts for combined tracks
const Double_t cuts[7] = {0.005,     // 0.005 cuts[0] = lowest V0 cut  (cm)
			  0.800,     // 0.015 cuts[1] = highest V0 cut (cm)
			  0.012,     // 0.012 cuts[2] = inv. mass cut (diferense) (Gev/c)
			  0.8,       // 0.8   cuts[3] = min. cosine value for pointing angle
			  -5000,     // -5000 cuts[4] = d0d0
			  0,         // 0.8   cuts[5] = costhetastar
			  0.5};      // 0.5   cuts[6] = ptchild
//cut for distance of closest aprach
double cutDCA=0.05;   //0.05

// this function applies single track cuts
Bool_t TrkCuts(const AliITStrackV2& trk);

// this function creates TObjArrays with positive and negative tracks
void   SelectTracks(TTree& itsTree,
                      TObjArray& trksP,Int_t* itsEntryP,Int_t& nTrksP,
                      TObjArray& trksN,Int_t* itsEntryN,Int_t& nTrksN);

//void GetPrimaryVertex(int i,Char_t* path="./");

//void PtD0(Char_t* path="./");

Int_t iTrkP,iTrkN,itsEntries;
Char_t trksName[100],refsName[100];
Int_t nTrksP=0,nTrksN=0;
Int_t nD0=0;
int ev=0;
double mom[6];
int event[10000];
int index=0;
Bool_t isSignal;
Int_t nTotEv=0,nD0recSgn=0,nD0recBkgS=0,nD0recBkg=0,nD0rec1ev=0;
Int_t pdg[2],mum[2],mumlab[2];
RECTRACK rectrk;

void RunD0offline(Char_t* path="./",bool h=false,bool PtD0=false) {
    
  AliKalmanTrack::SetConvConst(100/0.299792458/kBz);
  
  // Open file with ITS tracks
  Char_t fnameTrack[1024];
  //sprintf(fnameTrack,"%s/AliITStracksV2.root",path);
  sprintf(fnameTrack,"%s/AliITStracksV2.root",path);
  TFile* itstrks = TFile::Open(fnameTrack);
  
  Char_t refFile[1024];
  //sprintf(fnameTrack,"%s/ITStracksRefFile.root",path);
  sprintf(refFile,"%s/ITStracksRefFile.root",path);
  TFile* itsrefs = TFile::Open(refFile);
  

  //the offline stuff
  // define the cuts for vertexing
  Double_t vtxcuts[]={33., // max. allowed chi2
		      0.0, // min. allowed negative daughter's impact param 
		      0.0, // min. allowed positive daughter's impact param 
		      1.0, // max. allowed DCA between the daughter tracks
		      -1.0, // min. allowed cosine of V0's pointing angle
		      0.0, // min. radius of the fiducial volume
		      2.9};// max. radius of the fiducial volume
  
  TH1F *h1 = new TH1F("h1","Transvers momentun of reconstructed D0",100,0,10);
  TH1F *h2 = new TH1F("h2","Transvers momentun of D0 with |Eta|<0.9",100,0,10);
  TH1F *h3 = new TH1F("h3","Eta reconstructed of D0",100,-5,5);
  TH1F *h4 = new TH1F("h4","Eta of D0",100,-5,5);
  
  Char_t falice[1024];
  sprintf(falice,"%s/galice.root",path);
  TFile *f = new TFile(falice);   
  gAlice=(AliRun*)f->Get("gAlice");
  int nEvent=gAlice->GetEventsPerRun();
  //int nEvent=20;
  cout<<"#Events: "<<nEvent<<endl;
  delete gAlice;

  const Char_t *name="AliD0offline";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name); 
  
  for(ev=0;ev<nEvent;ev++) {
    
    cout<<"\n Event: "<<ev<<endl;
    
    // tracks from ITS
    sprintf(trksName,"TreeT_ITS_%d",ev);
    TTree *itsTree=(TTree*)itstrks->Get(trksName);
    if(!itsTree) continue;
    itsEntries = (Int_t)itsTree->GetEntries();
    printf("+++\n+++ Number of tracks in ITS: %d\n+++\n\n",itsEntries);
    
    // tree from reference file
    
    sprintf(refsName,"Tree_Ref_%d",ev);
    TTree *refTree=(TTree*)itsrefs->Get(refsName);
    refTree->SetBranchAddress("rectracks",&rectrk);

    //Getting primary Vertex
    GetPrimaryVertex(ev,path);
    
    // count the total number of events
    nTotEv++;

    // call function which applies sigle track selection and
    // separetes positives and negatives
    TObjArray trksP(itsEntries/2);
    Int_t *itsEntryP = new Int_t[itsEntries];
    TObjArray trksN(itsEntries/2);
    Int_t *itsEntryN = new Int_t[itsEntries];
    SelectTracks(*itsTree,trksP,itsEntryP,nTrksP,trksN,itsEntryN,nTrksN); 
    
    cout<<"#pos: "<<nTrksP<<endl;
    cout<<"#neg: "<<nTrksN<<endl;
    
    // create the AliV0vertexer object
    AliV0vertexer *vertexer2 = new AliV0vertexer(vtxcuts);
    
    AliD0Trigger * D0 = new AliD0Trigger(cuts,kBz,primaryvertex);
    
    double ptP,alphaP,phiP,ptN,alphaN,phiN,dca;
    
    for(iTrkP=0; iTrkP<nTrksP; iTrkP++) {
      postrack = (AliITStrackV2*)trksP.At(iTrkP);

      // get info on tracks PDG and mothers PDG from reference file
      refTree->GetEvent(itsEntryP[iTrkP]);
      pdg[0] = rectrk.pdg;
      mum[0] = rectrk.mumpdg;
      mumlab[0] = rectrk.mumlab;

      for(iTrkN=0; iTrkN<nTrksN; iTrkN++) {
	negtrack = (AliITStrackV2*)trksN.At(iTrkN);

	// get info on tracks PDG and mothers PDG from reference file
	refTree->GetEvent(itsEntryN[iTrkN]);
	pdg[1] = rectrk.pdg;
	mum[1] = rectrk.mumpdg;
	mumlab[1] = rectrk.mumlab;

	D0->SetTracks(postrack,negtrack);

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
	  //if(D0->FindV0offline(v2) && D0->d0d0()){
	  if(D0->d0d0()){
	    D0->SetV0(v2);
	    D0->FindMomentaOffline();
	    if(D0->FindInvMass() && D0->CosThetaStar() && D0->PointingAngle() && D0->pTchild()){
	      nD0++;
	      event[index]=ev; index++;
	      if(h){
		h1->Fill(D0->Pt());
		h3->Fill(D0->Eta());
	      }
	      
	      if(mumlab[0]==mumlab[1] && TMath::Abs(mum[0])==421 && TMath::Abs(mum[1])==421) { 
		nD0recSgn++;
	      } 
	      else { 
		if(TMath::Abs(mum[0])==421 || TMath::Abs(mum[0])==411 ||
		   TMath::Abs(mum[1])==421 || TMath::Abs(mum[1])==411) {
		  nD0recBkgS++;
		} else {
		  nD0recBkg++;
		}  
	      }
	    }
	  } 
	}
      }
    } 
    
    //delete D0;
    //delete itstrks;
    //delete itsTree;
    //delete trksP;
    //delete itsEntryP;
    //delete trksN;
    //delete itsEntryN;
  }
  
  cout<<"\nMy #D0: "<<nD0<<"\n"<<endl;
  printf("\n+++\n+++ Total number of events: %d\n+++\n",nTotEv);
  printf("\n+++\n+++ Total number of D0 candidates:\n
                 +++    Sgn:   %d\n
                 +++    BkgS:  %d\n
                 +++    Bkg:   %d\n+++\n",nD0recSgn,nD0recBkgS,nD0recBkg);

  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  
  if(h){
    if(!PtD0){
      TCanvas *c = new TCanvas("c","c",700,1000);
      c->Divide(1,2);
      c->cd(1);
      h1->Draw();
      pt = new TPaveText(7.3,0.4,11,1.8,"br");
      pt->SetFillColor(0);
      pt->SetBorderSize(1);
      pt->AddText("Cuts:");
      Char_t st[1024];
      sprintf(st,"First Pt:         %g",kPtCut);
      pt->AddText(st);
      sprintf(st,"d0 low:           %g",kd0Cut);
      pt->AddText(st);
      sprintf(st,"d0 high:          %g",kd0CutHigh);
      pt->AddText(st);
      sprintf(st,"V0 low:           %g",cuts[0]);
      pt->AddText(st);
      sprintf(st,"V0 high:          %g",cuts[1]);
      pt->AddText(st);
      sprintf(st,"InvMass Diff:     %g",cuts[2]);
      pt->AddText(st);
      sprintf(st,"cosPointingAngle: %g",cuts[3]);
      pt->AddText(st);
      sprintf(st,"d0d0:             %g",cuts[4]);
      pt->AddText(st);
      sprintf(st,"cosThetaStar:     %g",cuts[5]);
      pt->AddText(st);
      sprintf(st,"PtChild:          %g",cuts[6]);
      pt->AddText(st);
      sprintf(st,"DCA:              %g",cutDCA);
      pt->AddText(st);
      pt->Draw();
      c_1->Modified();
      c->cd();
      c->cd(2);
      h3->Draw();
    }
    if(PtD0){
      PtD0(path,h2,h4);
      TCanvas *c = new TCanvas("c","c",1000,700);
      c->Divide(2,2);
      c->cd(1);
      h1->Draw();
      pt = new TPaveText(7.3,0.4,11,1.8,"br");
      pt->SetFillColor(0);
      pt->SetBorderSize(1);
      pt->AddText("Cuts:");
      Char_t st[1024];
      sprintf(st,"First Pt:         %g",kPtCut);
      pt->AddText(st);
      sprintf(st,"d0 low:           %g",kd0Cut);
      pt->AddText(st);
      sprintf(st,"d0 high:          %g",kd0CutHigh);
      pt->AddText(st);
      sprintf(st,"V0 low:           %g",cuts[0]);
      pt->AddText(st);
      sprintf(st,"V0 high:          %g",cuts[1]);
      pt->AddText(st);
      sprintf(st,"InvMass Diff:     %g",cuts[2]);
      pt->AddText(st);
      sprintf(st,"cosPointAng: %g",cuts[3]);
      pt->AddText(st);
      sprintf(st,"d0d0:             %g",cuts[4]);
      pt->AddText(st);
      sprintf(st,"cosTheta*:     %g",cuts[5]);
      pt->AddText(st);
      sprintf(st,"PtChild:          %g",cuts[6]);
      pt->AddText(st);
      sprintf(st,"DCA:              %g",cutDCA);
      pt->AddText(st);
      pt->Draw();
      c_1->Modified();
      c->cd();
      c->cd(2);
      h3->Draw();
      c->cd(3);
      h2->Draw();
      c->cd(4);
      h4->Draw();
    }    
  }
  if(h){
    if(!PtD0){
      Char_t outName[1024];
      sprintf(outName,"%s/ReconstructedD0.root",path);
      TFile* outroot = new TFile(outName,"recreate");
      h1->Write();
      h3->Write();
      outroot->Close();
      delete outroot;
    }
  }
  if(PtD0){
    Char_t outName[1024];
    sprintf(outName,"%s/ReconstructedD0.root",path);
    TFile* outroot = new TFile(outName,"recreate");
    h1->Write();
    h2->Write();
    h3->Write();
    h4->Write();
    outroot->Close();
    delete outroot;
  }
  Char_t foutName[1024];
  sprintf(foutName,"%s/Cuts",path);
  ofstream fout(foutName);
  Char_t st2[1024];
  sprintf(st2,"First Pt:       %g",kPtCut);
  fout<<st2<<endl;
  sprintf(st2,"d0 low:         %g",kd0Cut);
  fout<<st2<<endl;
  sprintf(st2,"d0 high:        %g",kd0CutHigh);
  fout<<st2<<endl;
  sprintf(st2,"V0 low:         %g",cuts[0]);
  fout<<st2<<endl;
  sprintf(st2,"V0 high:        %g",cuts[1]);
  fout<<st2<<endl;
  sprintf(st2,"InvMass Diff:   %g",cuts[2]);
  fout<<st2<<endl;
  sprintf(st2,"cosPointAng:    %g",cuts[3]);
  fout<<st2<<endl;
  sprintf(st2,"d0d0:           %g",cuts[4]);
  fout<<st2<<endl;
  sprintf(st2,"cosTheta*:      %g",cuts[5]);
  fout<<st2<<endl;
  sprintf(st2,"PtChild:        %g",cuts[6]);
  fout<<st2<<endl;
  sprintf(st2,"DCA:            %g",cutDCA);
  fout<<st2<<endl;
  fout.close();

  Char_t fName[1024];
  sprintf(fName,"%s/Events",path);
  ofstream fevent(fName);
  for(int i=0;i<nD0;i++){fevent<<event[i]<<endl;}
  fevent.close();
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
  //if(TMath::Abs(10000.*trk.GetD(primaryvertex[0],primaryvertex[1])) > kd0CutHigh) return kFALSE;

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
//____________________________________________________________________________
PtD0(Char_t* path="./",TH1F * h1,TH1F * h4){
  
  Char_t falice[1024];
  sprintf(falice,"%s/galice.root",path);  
  TFile *f = new TFile(falice);
  gAlice=(AliRun*)f->Get("gAlice");
  
  TParticle *p;
  int nd0=0;
  int nkminus =0;
  int npipluss = 0;
  int nEvent=gAlice->GetEventsPerRun();
  
  for (Int_t i = 0; i <nEvent; i++) {
    cout<<"Event: "<<i<<endl;
    gAlice->GetEvent(i);
    Int_t nPart = gAlice->GetNtrack();
    for (Int_t iPart = 0; iPart < nPart; iPart++) {
      //cout<<"Particlenr.: "<<iPart<<endl;
      p = (TParticle*)gAlice->Particle(iPart);
      if (p->GetPdgCode()==421){
	if(fabs(p->Eta())<0.9) h1->Fill(p.Pt());
	h4->Fill(p.Eta());
	nd0++;
      }
      if (p->GetPdgCode()==-321){
	nkminus++;
      }
      if (p->GetPdgCode()==211){
	npipluss++;
      }
    }
  }
  delete gAlice;
}
