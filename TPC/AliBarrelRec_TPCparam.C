/****************************************************************************
 * This macro is supposed to do reconstruction in the ITS via Kalman        *
 * tracker V2. The ITStracker is feeded with parametrized TPC tracks        * 
 *                                                                          *
 * It does the following steps:                                             *
 *             1) TPC tracking parameterization                             *
 *             2) Fast points in ITS                                        *
 *             3) ITS cluster finding V2                                    *
 *             4) Determine z position of primary vertex                    *
 *                - read from event header in PbPb events                   *
 *                - determined using points on SPD in pp events (M.Masera)  *
 *             5) ITS track finding V2                                      *
 *                * in pp, redetermine the z position of the primary vertex *
 *                  using the reconstructed tracks
 *             6) Create a reference file with simulation info (p,PDG...)   *
 *                                                                          *
 * (Origin: A.Dainese, Padova, andrea.dainese@pd.infn.it                    * 
 *  from AliTPCtest.C & AliITStestV2.C by I.Belikov                         *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
//-- --- standard headers------------- 
#include "Riostream.h"
//--------Root headers ---------------
#include <TFile.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TSystem.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TParticle.h>
#include <TRandom.h>
//----- AliRoot headers ---------------
#include "alles.h"
#include "AliRun.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliArrayI.h"
#include "AliDigits.h"
#include "AliITS.h"
#include "AliTPC.h"
#include "AliITSgeom.h"
#include "AliITSRecPoint.h"
#include "AliITSclusterV2.h"
#include "AliITSsimulationFastPoints.h" 
#include "AliITStrackerV2.h"
#include "AliKalmanTrack.h"
#include "AliTPCtrackerParam.h"
//-------------------------------------
#endif

// structure for track references
typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumlab;
  Int_t mumpdg;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} RECTRACK;

//===== Functions definition ================================================= 

Int_t MarkEvtsToSkip(const Char_t *evtsName,
		           Bool_t *skipEvt);

Int_t TPCParamTracks(const Char_t  *galiceName,
		     const Char_t  *outName,
		     const Int_t    coll,
		     const Double_t Bfield,
		           Int_t    n);

Int_t ITSHits2FastRecPoints(const Char_t *galiceName,
			          Bool_t *skipEvt,
			          Int_t  n);

Int_t ITSFindClustersV2(const Char_t *galiceName,
			const Char_t *outName,
			      Bool_t *skipEvt,
			      Int_t   n);

Int_t ZvtxFromHeader(const Char_t *galiceName,
		           Bool_t *skipEvt,
		           Int_t   n);

Int_t ZvtxFromSPD(const Char_t *galiceName,
		        Bool_t *skipEvt,
		        Int_t   n);

Int_t ZvtxFromTracks(const Char_t *trkame,
		           Bool_t *skipEvt,
		           Int_t  n);

Int_t ZvtxFastpp(const Char_t *galiceName,
		       Bool_t *skipEvt,
		        Int_t  n);

void EvalZ(TH1F    *hist,
	   Int_t    sepa,
	  Float_t  &av,
	  Float_t  &sig,
	  Int_t     ncoinc,
	  TArrayF  *zval,
	  ofstream *deb);

Bool_t VtxTrack(Double_t pt,
		Double_t d0rphi);

Double_t d0zRes(Double_t pt);

Int_t ITSFindTracksV2(const Char_t   *galiceName,
		      const Char_t   *inName,
		      const Char_t   *inName2,
		      const Char_t   *outName,
		            Bool_t   *skipEvt,
			    Option_t *vtxMode,
			    Int_t     n);

Int_t ITSMakeRefFile(const Char_t *galiceName,
		     const Char_t *inName, 
		     const Char_t *outName,
		           Bool_t *skipEvt,
		            Int_t  n);

void WriteZvtx(const Char_t   *name,
		     Double_t *zvtx,
		     Int_t     n);

void ReadZvtx(const Char_t   *name,
	            Double_t *zvtx,
	            Int_t     n);

//=============================================================================

Int_t AliBarrelRec_TPCparam(Int_t n=1) {

  const Char_t *name=" AliBarrelRec_TPCparam";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  const Char_t *evtsName    = "EvtsToSkip.dat";
  const Char_t *TPCtrkNameS = "AliTPCtracksParam.root";
  const Char_t *galiceName  = "galice.root";
  const Char_t *ITSclsName  = "AliITSclustersV2.root";
  const Char_t *ITStrkName  = "AliITStracksV2.root";
  const Char_t *ITStrkName2 = "AliITStracksV2_2.root";
  const Char_t *ITSrefName  = "ITStracksRefFile.root";

  // set here the code for the type of collision (needed for TPC tracking
  // parameterization). available collisions:
  //
  // coll = 0 ->   PbPb6000 (HIJING with b<2fm) 
  // coll = 1 ->   pp 
  const Int_t    collcode    = 1;  
  const Bool_t   slowVtx     = kFALSE;
  const Bool_t   retrack     = kFALSE;
  // set here the value of the magnetic field
  const Double_t BfieldValue = 0.4;
 
  AliKalmanTrack::SetConvConst(100/0.299792458/BfieldValue);

  Bool_t *skipEvt = new Bool_t[n];

  // Mark events that have to be skipped (read from file ascii evtsName)
  for(Int_t i=0;i<n;i++) skipEvt[i] = kFALSE;
  if(!gSystem->AccessPathName(evtsName,kFileExists)) { 
    MarkEvtsToSkip(evtsName,skipEvt);
  }
  
  // ********** Build TPC tracks with parameterization *********** //
  TPCParamTracks(galiceName,TPCtrkNameS,collcode,BfieldValue,n);
  
    
  // ********** ITS RecPoints ************************************ //
  ITSHits2FastRecPoints(galiceName,skipEvt,n);
  
  
  // ********** Find ITS clusters ******************************** //
  ITSFindClustersV2(galiceName,ITSclsName,skipEvt,n);
  

  // ********** Tracking in ITS ********************************** //
  Char_t *vtxMode;
  //  Pb-Pb
  if(collcode==0) {
    ZvtxFromHeader(galiceName,skipEvt,n);
    vtxMode="Header";
    ITSFindTracksV2(galiceName,TPCtrkNameS,ITSclsName,ITStrkName,skipEvt,vtxMode,n); 
  }
  //   pp
  if(collcode==1 && slowVtx) {
    ZvtxFromSPD(galiceName,skipEvt,n);  
    vtxMode="SPD";
    ITSFindTracksV2(galiceName,TPCtrkNameS,ITSclsName,ITStrkName,skipEvt,vtxMode,n); 
    ZvtxFromTracks(ITStrkName,skipEvt,n);
    if(retrack) {
      vtxMode="Tracks"; 
      ITSFindTracksV2(galiceName,TPCtrkNameS,ITSclsName,ITStrkName2,skipEvt,vtxMode,n);
    }
  }  
  if(collcode==1 && !slowVtx) {
    ZvtxFastpp(galiceName,skipEvt,n);  
    vtxMode="Fast";
    ITSFindTracksV2(galiceName,TPCtrkNameS,ITSclsName,ITStrkName,skipEvt,vtxMode,n); 
  }

  // ********** Make ITS tracks reference file ******************* //
  ITSMakeRefFile(galiceName,ITStrkName,ITSrefName,skipEvt,n);
  
  

  delete [] skipEvt;

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t MarkEvtsToSkip(const Char_t *evtsName,Bool_t *skipEvt) {

  cerr<<"\n*******************************************************************\n";
  cerr<<"\nChecking for events to skip...\n";

  Int_t evt,ncol;

  FILE *f = fopen(evtsName,"r");
  while(1) {
    ncol = fscanf(f,"%d",&evt);
    if(ncol<1) break;
    skipEvt[evt] = kTRUE;
    cerr<<" ev. "<<evt<<" will be skipped\n";
  }
  fclose(f);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t TPCParamTracks(const Char_t *galiceName,const Char_t *outName,
		     const Int_t coll,const Double_t Bfield,Int_t n) {

  cerr<<"\n*******************************************************************\n";

  const Char_t *name="TPCParamTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  TFile *outFile=TFile::Open(outName,"recreate");
  TFile *inFile =TFile::Open(galiceName);
 
  AliTPCtrackerParam tracker(coll,Bfield,n);
  tracker.BuildTPCtracks(inFile,outFile);

  delete gAlice; gAlice=0;

  inFile->Close();
  outFile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t ITSHits2FastRecPoints(const Char_t *galiceName,Bool_t *skipEvt,Int_t n) {
 
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ITSHits2FastRecPoints";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
 
  Int_t nsignal=25;
  Int_t size=-1;

  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = TFile::Open(galiceName,"UPDATE");

  gAlice = (AliRun*)file->Get("gAlice");


  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  if(!ITS) return 1;

  // Set the simulation model
  for(Int_t i=0;i<3;i++) {
    ITS->SetSimulationModel(i,new AliITSsimulationFastPoints());
  }
   

  //
  // Event Loop
  //
  for(Int_t ev=0; ev<n; ev++) {
    if(skipEvt[ev]) continue;
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;
    Int_t nparticles = gAlice->GetEvent(ev);
    cerr<<"Number of particles: "<<nparticles<<endl;
    gAlice->SetEvent(ev);
    if(!gAlice->TreeR()) gAlice->MakeTree("R");
    if(!gAlice->TreeR()->GetBranch("ITSRecPointsF")) {  
      ITS->MakeBranch("RF");
      if(nparticles <= 0) return 1;

      Int_t bgr_ev=Int_t(ev/nsignal);
      //printf("bgr_ev %d\n",bgr_ev);
      ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
    }
  } // event loop 

  delete gAlice; gAlice=0;
  file->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t ITSFindClustersV2(const Char_t *galiceName,const Char_t *outName,
			Bool_t *skipEvt,Int_t n) {
  //
  // This function converts AliITSRecPoint(s) to AliITSclusterV2
  //
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ITSFindClustersV2";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);


  TFile *inFile=TFile::Open(galiceName);
  if(!inFile->IsOpen()) { cerr<<"Can't open galice.root !\n"; return 1; }
   
  if(!(gAlice=(AliRun*)inFile->Get("gAlice"))) {
    cerr<<"Can't find gAlice !\n";
    return 1;
  }
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  if(!ITS) { cerr<<"Can't find the ITS !\n"; return 1; }
  AliITSgeom *geom=ITS->GetITSgeom();

  TFile *outFile = TFile::Open(outName,"recreate");
  geom->Write();

  // loop on events
  for(Int_t ev=0; ev<n; ev++) {
    if(skipEvt[ev]) continue;
    cerr<<"--- Processing event "<<ev<<" ---"<<endl;
    inFile->cd();
    gAlice->GetEvent(ev);

    TClonesArray *clusters = new TClonesArray("AliITSclusterV2",10000);
    Char_t cname[100];
    sprintf(cname,"TreeC_ITS_%d",ev);
    TTree *cTree = new TTree(cname,"ITS clusters");
    cTree->Branch("Clusters",&clusters);

    TTree *pTree=gAlice->TreeR();
    if(!pTree) { cerr<<"Can't get TreeR !\n"; return 1; }
    TBranch *branch = pTree->GetBranch("ITSRecPointsF");
    if(!branch) { cerr<<"Can't get ITSRecPoints branch !\n"; return 1; }
    TClonesArray *points = new TClonesArray("AliITSRecPoint",10000);
    branch->SetAddress(&points);

    TClonesArray &cl=*clusters;
    Int_t nclusters=0;
    Int_t nentr=(Int_t)branch->GetEntries();

    cerr<<"Number of entries in TreeR_RF: "<<nentr<<endl;

    Float_t *lp  = new Float_t[5]; 
    Int_t   *lab = new Int_t[6];

    for(Int_t i=0; i<nentr; i++) {
      points->Clear();
      branch->GetEvent(i);
      Int_t ncl=points->GetEntriesFast(); if(ncl==0){cTree->Fill();continue;}
      Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
      Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
      Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
      Double_t yshift = x*rot[0] + y*rot[1];
      Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
      nclusters+=ncl;
      
      Float_t kmip=1; // ADC->mip normalization factor for the SDD and SSD 
      if(lay==4 || lay==3){kmip=280.;};
      if(lay==6 || lay==5){kmip=38.;};

      for(Int_t j=0; j<ncl; j++) {
	AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
	lp[0]=-p->GetX()-yshift; if(lay==1) lp[0]=-lp[0];
	lp[1]=p->GetZ()+zshift;
	lp[2]=p->GetSigmaX2();
	lp[3]=p->GetSigmaZ2();
	lp[4]=p->GetQ(); lp[4]/=kmip;
	lab[0]=p->GetLabel(0);lab[1]=p->GetLabel(1);lab[2]=p->GetLabel(2);
	lab[3]=ndet;
	
	Int_t label=lab[0];
	if(label>=0) {
	  TParticle *part=(TParticle*)gAlice->Particle(label);
	  label=-3;
	  while (part->P() < 0.005) {
	    Int_t m=part->GetFirstMother();
	    if(m<0) {cerr<<"Primary momentum: "<<part->P()<<endl; break;}
	    label=m;
	    part=(TParticle*)gAlice->Particle(label);
	  }
	  if      (lab[1]<0) lab[1]=label;
	  else if(lab[2]<0) lab[2]=label;
	  else cerr<<"No empty labels !\n";
	}
	
	new(cl[j]) AliITSclusterV2(lab,lp);
      }
      cTree->Fill(); clusters->Delete();
      points->Delete();
    }

    outFile->cd();
    cTree->Write();

    cerr<<"Number of clusters: "<<nclusters<<endl;
    delete [] lp;
    delete [] lab;
    delete cTree; delete clusters; delete points;

  } // end loop on events
    

  delete gAlice; gAlice=0;
  
  inFile->Close();
  outFile->Close();
  
  gBenchmark->Stop(name);
  gBenchmark->Show(name);

   return 0;
}
//-----------------------------------------------------------------------------
Int_t ZvtxFromHeader(const Char_t *galiceName,
		     Bool_t *skipEvt,Int_t n) {

  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ZvtxFromHeader";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
 
  TFile *galice = TFile::Open(galiceName);  

  Double_t *zvtx    = new Double_t[n];

  for(Int_t ev=0; ev<n; ev++){
    if(skipEvt[ev]) continue;
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;
    
    TArrayF o = 0;
    o.Set(3);
    AliHeader* header = 0;
    TTree* treeE = (TTree*)gDirectory->Get("TE");
    treeE->SetBranchAddress("Header",&header);
    treeE->GetEntry(ev);
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if(genHeader) {
      // get primary vertex position
      genHeader->PrimaryVertex(o);
      // set position of primary vertex
      zvtx[ev] = (Double_t)o[2];  
    } else {
      cerr<<" ! event header not found : setting z vertex to 0 !"<<endl;
      zvtx[ev] = 0.;
    }    
    delete header;
  }

  galice->Close();

  // Write vertices to file
  WriteZvtx("zvtxHeader.dat",zvtx,n);

  delete [] zvtx;
 
  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t ZvtxFastpp(const Char_t *galiceName,
		 Bool_t *skipEvt,Int_t n) {

  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ZvtxFastpp";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
 
  TFile *galice = TFile::Open(galiceName);  
  AliRun *gAlice = (AliRun*)galice->Get("gAlice");

  Double_t *zvtx    = new Double_t[n];

  TParticle *part;
  Int_t      nPart,b;
  Double_t dNchdy,E,theta,eta,zvtxTrue,sigmaz,probVtx;

  for(Int_t ev=0; ev<n; ev++) {
    if(skipEvt[ev]) continue;
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;

    TArrayF o = 0;
    o.Set(3);
    AliHeader* header = 0;
    TTree* treeE = (TTree*)galice->Get("TE");
    treeE->SetBranchAddress("Header",&header);
    treeE->GetEntry(ev);
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if(genHeader) {
      // get primary vertex position
      genHeader->PrimaryVertex(o);
      // set position of primary vertex
      zvtxTrue = (Double_t)o[2];  
    } else {
      cerr<<" ! event header not found : setting z vertex to 0 !"<<endl;
      zvtxTrue = 0.;
    }    
    delete header;


    // calculate dNch/dy pf the event (charged primaries in |eta|<0.5)
    nPart = (Int_t)gAlice->GetEvent(ev);
    dNchdy = 0.;
    for(b=0; b<nPart; b++) {
      part = (TParticle*)gAlice->Particle(b);
      if(TMath::Abs(part->GetPdgCode())<10) continue; // reject partons
      if(TMath::Abs(part->Vx()*part->Vx()+part->Vy()*part->Vy())>0.01) continue; // reject secondaries
      E  = part->Energy();
      if(E>6900.) continue; // reject incoming protons
      theta = part->Theta();
      if(theta<.1 || theta>3.) continue;
      eta = -TMath::Log(TMath::Tan(theta/2.));
      if(TMath::Abs(eta)>0.5) continue; // count particles in |eta|<0.5
      dNchdy+=1.;
    }

    // get sigma(z) corresponding to the mult of this event
    if(dNchdy>1.)  {
      sigmaz   = 0.0417/TMath::Sqrt(dNchdy);
    } else {
      sigmaz = 0.0500;
    }
    // smear the original position of the primary vertex
    zvtx[ev] = gRandom->Gaus(zvtxTrue,sigmaz);

    // compute the probability that the vertex is found
    probVtx = 1.;
    if(dNchdy<24.) probVtx = 0.85+0.00673*dNchdy;
    if(dNchdy<14.) probVtx = 0.85;

    if(gRandom->Rndm()>probVtx) zvtx[ev] = -1000.;// no zvtx for this event

  }
  galice->Close();

  // Write vertices to file
  WriteZvtx("zvtxFastpp.dat",zvtx,n);

  delete [] zvtx;
  //delete gAlice; gAlice=0;

 
  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//-----------------------------------------------------------------------------
Int_t ZvtxFromSPD(const Char_t *galiceName,Bool_t *skipEvt,Int_t n) {


  Double_t *zvtx    = new Double_t[n];

  Bool_t debug = kFALSE;
  ofstream fildeb;
  if(debug)fildeb.open("FindZV.out");
  ofstream filou;
  //filou.open("zcoor.dat");
  const Float_t kPi=TMath::Pi();
  const Int_t kFirstL1=0;
  const Int_t kLastL1=79;
  const Int_t kFirstL2=80;
  const Int_t kLastL2=239;
  //  const Float_t kDiffPhiMax=0.0005;
  const Float_t kDiffPhiMax=0.01;
  const Float_t kphl=0.;
  const Float_t kphh=kPi/18.;
  TDirectory *currdir = gDirectory;
  TH2F *hz2z1 = 0;
  TH1F *zvdis = 0;
  TH2F *dpvsz = 0;
  TProfile *scoc = new TProfile("scoc","Zgen - Zmeas vs occ lay 1",5,10.,60.);
  TProfile *prof = new TProfile("prof","z2 vs z1",50,-15.,15.);
  TH1F *phdiff = new TH1F("phdiff","#Delta #phi",100,kphl,kphh);
  TH1F *phdifsame = new TH1F("phdifsame","#Delta #phi same track",100,kphl,kphh);
  if(gAlice){delete gAlice; gAlice = 0;}
  TFile *file = TFile::Open(galiceName);
  gAlice = (AliRun*)file->Get("gAlice");

  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  AliITSgeom *geom = ITS->GetITSgeom();
  TTree *TR=0;
  Float_t lc[3]; for(Int_t ii=0; ii<3; ii++) lc[ii]=0.;
  Float_t gc[3]; for(Int_t ii=0; ii<3; ii++) gc[ii]=0.;
  Float_t lc2[3]; for(Int_t ii=0; ii<3; ii++) lc2[ii]=0.;
  Float_t gc2[3]; for(Int_t ii=0; ii<3; ii++) gc2[ii]=0.;
  char name[30];
  Float_t avesca=0;
  Float_t avesca2=0;
  Int_t N=0;
  for(Int_t event=0; event<n; event++){
    if(skipEvt[event]) continue;
    sprintf(name,"event_%d",event);
    hz2z1=new TH2F(name,"z2 vs z1",50,-15.,15.,50,-15.,15.);
    hz2z1->SetDirectory(currdir);
    gAlice->GetEvent(event);
    TParticle * part = gAlice->Particle(0);
    Float_t oriz=part->Vz();
    cout<<"\n ==============================================================\n";
    cout<<"   Processing event "<<event<<" "<<oriz<<endl;
    cout<<" ==============================================================\n";
    if(debug){
    fildeb<<"\n ==============================================================\n";
    fildeb<<"   Processing event "<<event<<" "<<oriz<<endl;
    fildeb<<" ==============================================================\n";
    }
    file->cd();
    TR = gAlice->TreeR();
    if(!TR) cerr<<"TreeR not found\n";
    TBranch *branch=TR->GetBranch("ITSRecPointsF");
    if(!branch) cerr<<"Branch ITSRecPointsF not found\n"; 
    TClonesArray *points  = new TClonesArray("AliITSRecPoint",10000);
    branch->SetAddress(&points);
    Int_t nentr=(Int_t)branch->GetEntries();
    Float_t zave=0;
    Float_t rmszav=0;
    Float_t zave2=0;
    Int_t firipixe=0;
    cout<<endl;
    for(Int_t module= kFirstL1; module<=kLastL1;module++){
      points->Clear();
      branch->GetEvent(module);
      Int_t nrecp1 = points->GetEntriesFast();
      //cerr<<nrecp1<<endl;
      for(Int_t i=0; i<nrecp1;i++){
        AliITSRecPoint *current = (AliITSRecPoint*)points->UncheckedAt(i);
        lc[0]=current->GetX();
        lc[2]=current->GetZ();
        geom->LtoG(module,lc,gc);
        zave+=gc[2];
        zave2+=gc[2]*gc[2];
        firipixe++;
      }
      points->Delete();
    }
    if(firipixe>1){
      rmszav=TMath::Sqrt(zave2/(firipixe-1)-zave*zave/firipixe/(firipixe-1));
      zave=zave/firipixe;
      cout<<"Z aver first layer = "<<zave<<" RMS= "<<rmszav<<" pixels = "<<firipixe<<endl;
      if(debug)fildeb<<"Z aver first layer = "<<zave<<" RMS= "<<rmszav<<" pixels = "<<firipixe<<endl;
    }
    else {
      cout<<"No rec points on first layer for this event\n";
      if(debug)fildeb<<"No rec points on first layer for this event\n";
    }
    Float_t zlim1=zave-rmszav;
    Float_t zlim2=zave+rmszav;
    Int_t sepa=(Int_t)((zlim2-zlim1)*10.+1.);
    zlim2=zlim1 + sepa/10.;
    sprintf(name,"pz_ev_%d",event);
    dpvsz=new TH2F(name,"#Delta #phi vs Z",sepa,zlim1,zlim2,100,0.,0.03);
    dpvsz->SetDirectory(currdir);
    sprintf(name,"z_ev_%d",event);
    zvdis=new TH1F(name,"zv distr",sepa,zlim1,zlim2);
    zvdis->SetDirectory(currdir);
    cout<<"Z limits: "<<zlim1<<" "<<zlim2<<endl;
    if(debug)fildeb<<"Z limits: "<<zlim1<<" "<<zlim2<<endl;
    Int_t sizarr=100;
    TArrayF *zval = new TArrayF(sizarr);
    Int_t ncoinc=0;
    for(Int_t module= kFirstL1; module<=kLastL1;module++){
      //cout<<"processing module   "<<module<<"                  \r";
      branch->GetEvent(module);
      Int_t nrecp1 = points->GetEntriesFast();
      TObjArray *poiL1 = new TObjArray(nrecp1);
      for(Int_t i=0; i<nrecp1;i++)poiL1->AddAt(points->UncheckedAt(i),i);
      //ITS->ResetRecPoints();
      points->Delete();
      for(Int_t i=0; i<nrecp1;i++){
        AliITSRecPoint *current = (AliITSRecPoint*)poiL1->UncheckedAt(i);
        lc[0]=current->GetX();
        lc[2]=current->GetZ();
        geom->LtoG(module,lc,gc);
        Float_t r1=TMath::Sqrt(gc[0]*gc[0]+gc[1]*gc[1]);
        Float_t phi1 = TMath::ATan2(gc[1],gc[0]);
        if(phi1<0)phi1=2*kPi+phi1;
        Int_t lab1 = current->GetLabel(0);
        Float_t phi1d=phi1*180./kPi;
        //cerr<<"module "<<module<<" "<<gc[0]<<" "<<gc[1]<<" "<<gc[2]<<" "<<phi1d<<" "<<lab1<<"     \n";
        for(Int_t modul2=kFirstL2; modul2<=kLastL2; modul2++){
          branch->GetEvent(modul2);
          Int_t nrecp2 = points->GetEntriesFast();
          for(Int_t j=0; j<nrecp2;j++){
            AliITSRecPoint *recp = (AliITSRecPoint*)points->UncheckedAt(j);
            lc2[0]=recp->GetX();
            lc2[2]=recp->GetZ();
            geom->LtoG(modul2,lc2,gc2);
            Float_t r2=TMath::Sqrt(gc2[0]*gc2[0]+gc2[1]*gc2[1]);
            Float_t zr0=(r2*gc[2]-r1*gc2[2])/(r2-r1);
            Int_t lab2 = recp->GetLabel(0);
            Float_t phi2 = TMath::ATan2(gc2[1],gc2[0]);
            if(phi2<0)phi2=2.*kPi+phi2;
            Float_t diff = TMath::Abs(phi2-phi1);
            if(diff>kPi)diff=2.*kPi-diff;
            if(lab1==lab2)phdifsame->Fill(diff);
            if(zr0>zlim1 && zr0<zlim2){
              phdiff->Fill(diff);
              dpvsz->Fill(zr0,diff);
              if(diff<kDiffPhiMax ){
                zvdis->Fill(zr0);
                zval->AddAt(zr0,ncoinc);
                ncoinc++;
                if(ncoinc==(sizarr-1)){
                  sizarr+=100;
                  zval->Set(sizarr);
                }
                Float_t phi2d=phi2*180./kPi;
                Float_t diffd=diff*180./kPi;
                //cerr<<"module2 "<<modul2<<" "<<gc2[0]<<" "<<gc2[1]<<" "<<gc2[2]<<" "<<phi2d<<" "<<diffd<<" "<<lab2<<"     \n";
                //                cout<<"zr0= "<<zr0<<endl;
                hz2z1->Fill(gc[2],gc2[2]);
                prof->Fill(gc[2]-oriz,gc2[2]-oriz);
              }
            }
          }
	  points->Delete();
        }
      }
      delete poiL1;
    }         // loop on modules
    cout<<endl;
    Float_t zcmp=0;
    Float_t zsig=0;
    EvalZ(zvdis,sepa,zcmp,zsig,ncoinc,zval,&fildeb);
    if(zcmp!=-100 && zsig!=-100){
      N++;
      Float_t scarto = (oriz-zcmp)*10000.;
      Float_t ascarto=TMath::Abs(scarto);
      scoc->Fill(firipixe,ascarto);
      avesca+=ascarto;
      avesca2+=scarto*scarto;
      if(debug)fildeb<<"Mean value: "<<zcmp<<"  +/-"<<zsig<<"  Zgen- Zmeas (micron)= "<<scarto<<endl;
      cout<<"Mean value: "<<zcmp<<"  RMS: "<<zsig<<"  Zgen- Zmeas (micron)= "<<scarto<<endl;
      //filou<<event<<" "<<zcmp<<" "<<zsig<<" "<<oriz<<" "<<scarto<<endl;
      zvtx[event] = (Double_t)zcmp;
    }
    else {
      cout<<"Not enough points in event "<<event<<endl;
      if(debug)fildeb<<"Not enough points\n";
      //filou<<event<<" "<<-1000.<<" "<<0.<<" "<<oriz<<" "<<0.<<endl;
      zvtx[event] = -1000.;
    }
    delete zval;
  }           // loop on events
  file->Close();
  //hz2z1->Draw("box");

  // Write vertices to file
  WriteZvtx("zvtxSPD.dat",zvtx,n);
  delete [] zvtx;

  if(N>1) {
    avesca2=TMath::Sqrt((avesca2/(N-1)-avesca*avesca/N/(N-1))/N);
    avesca=avesca/N;
  } else {
    avesca2 = 0;
  }
  cout<<"\n \n ==========================================================\n";
  cout<<"Number of events with vertex estimate  "<<N<<endl;
  cout<<"Average of the (abs) Zgen-Zmeas   :   "<<avesca;
  cout<<" +/- "<<avesca2<<" micron"<<endl;
  if(debug){
    fildeb<<"\n \n ==========================================================\n";
    fildeb<<"Number of events with vertex estimate  "<<N<<endl;
    fildeb<<"Average of the (abs) Zgen-Zmeas   :   "<<avesca;
    fildeb<<" +/- "<<avesca2<<endl;
  }
  return 0;
}


//-----------------------------------------------------------------------------
void EvalZ(TH1F *hist,Int_t sepa,Float_t &av, Float_t &sig,Int_t ncoinc, 
	   TArrayF *zval,ofstream *deb) {

  av=0;
  sig=0;
  Int_t N=0;
  Int_t NbinNotZero=0;
  for(Int_t i=1;i<=sepa;i++){
    Float_t cont=hist->GetBinContent(i);
    av+=cont;
    sig+=cont*cont;
    N++;
    if(cont!=0)NbinNotZero++;
  }
  if(NbinNotZero==0){av=-100; sig=-100; return;}
  if(NbinNotZero==1){sig=-100; return;}
  sig=TMath::Sqrt(sig/(N-1)-av*av/N/(N-1));
  av=av/N;
  Float_t val1=hist->GetBinLowEdge(sepa); 
  Float_t val2=hist->GetBinLowEdge(1);
  for(Int_t i=1;i<=sepa;i++){
    Float_t cont=hist->GetBinContent(i);
    if(cont>(av+sig*2.)){
      Float_t curz=hist->GetBinLowEdge(i);
      if(curz<val1)val1=curz;
      if(curz>val2)val2=curz;
    }
  }
  val2+=hist->GetBinWidth(1);
  cout<<"Values for Z finding: "<<val1<<" "<<val2<<endl;
  if(deb)(*deb)<<"Values for Z finding: "<<val1<<" "<<val2<<endl;
  av=0;
  sig=0;
  N=0;
  for(Int_t i=0; i<ncoinc; i++){
    Float_t z=zval->At(i);
    if(z<val1)continue;
    if(z>val2)continue;
    av+=z;
    sig+=z*z;
    N++;
  }
  if(N<=1){av=-100; sig=-100; return;}
  sig=TMath::Sqrt((sig/(N-1)-av*av/N/(N-1))/N);
  av=av/N;
  return;
}
//-----------------------------------------------------------------------------
Int_t ZvtxFromTracks(const Char_t *trkName,Bool_t *skipEvt,Int_t n) {

  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ZvtxFromTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  Double_t *zvtx = new Double_t[n];
  Double_t *zvtxSPD = new Double_t[n];
  ReadZvtx("zvtxSPD.dat",zvtxSPD,n);

  TFile *itstrk = TFile::Open(trkName);


  Int_t    nTrks,i,nUsedTrks;
  Double_t pt,d0rphi,d0z; 
  Double_t zvtxTrks,ezvtxTrks,sumWeights;

  for(Int_t ev=0; ev<n; ev++){
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;

    // Tree with ITS tracks
    Char_t tname[100];
    sprintf(tname,"TreeT_ITS_%d",ev);

    TTree *tracktree=(TTree*)itstrk->Get(tname);
    if(!tracktree) continue;
    AliITStrackV2 *itstrack=new AliITStrackV2; 
    tracktree->SetBranchAddress("tracks",&itstrack);
    nTrks = (Int_t)tracktree->GetEntries();
    nUsedTrks = 0;
    zvtxTrks = 0.;
    ezvtxTrks = 0.;
    sumWeights = 0.;

    //cerr<<" ITS tracks: "<<nTrks<<endl;
    // loop on tracks
    for(i=0; i<nTrks; i++) {
      tracktree->GetEvent(i);
      // get pt and impact parameters
      pt = 1./TMath::Abs(itstrack->Get1Pt());
      d0rphi = TMath::Abs(itstrack->GetD());
      itstrack->PropagateToVertex();
      d0z = itstrack->GetZ();
      // check if the track is from (0,0) in (x,y)
      if(!VtxTrack(pt,d0rphi)) continue;
      nUsedTrks++;
      zvtxTrks   += d0z/d0zRes(pt)/d0zRes(pt);
      sumWeights += 1./d0zRes(pt)/d0zRes(pt);
 
    } // loop on tracks  

    if(nUsedTrks>1) {
      // estimated position in z of vertex
      zvtxTrks /= sumWeights;
      // estimated error
      ezvtxTrks = TMath::Sqrt(1./sumWeights);
    } else {
      zvtxTrks = zvtxSPD[ev];
    }

    zvtx[ev] = zvtxTrks;

  } // loop on events


  // Write vertices to file
  WriteZvtx("zvtxTracks.dat",zvtx,n);
  delete [] zvtx;
  delete [] zvtxSPD;

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return 0;
}
//----------------------------------------------------------------------------
Bool_t VtxTrack(Double_t pt,Double_t d0rphi) {

  Double_t d0rphiRes = TMath::Sqrt(11.59*11.59+65.76*65.76/TMath::Power(pt,1.878))/10000.;
  Double_t nSigma = 3.;
  if(d0rphi<nSigma*d0rphiRes) { return kTRUE; } else { return kFALSE; }
}
//----------------------------------------------------------------------------
Double_t d0zRes(Double_t pt) {
  Double_t res = TMath::Sqrt(34.05*34.05+170.1*170.1/TMath::Power(pt,1.226));
  return res/10000.;
}
//-----------------------------------------------------------------------------
Int_t ITSFindTracksV2(const Char_t *galiceName,const Char_t *inName, 
		      const Char_t *inName2,const Char_t *outName, 
		      Bool_t *skipEvt,Option_t *vtxMode,Int_t n) {

  
  cerr<<"\n*******************************************************************\n";

  const Char_t *name="ITSFindTracksV2";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  // Read vertices from file 
  Double_t *zvtx = new Double_t[n];
  Char_t *vtxfile="zvtxHeader.dat";

  const Char_t *header = strstr(vtxMode,"Header");
  const Char_t *fastpp = strstr(vtxMode,"Fast");
  const Char_t *spd    = strstr(vtxMode,"SPD");
  const Char_t *tracks = strstr(vtxMode,"Tracks");

  if(header) vtxfile = "zvtxHeader.dat";
  if(fastpp) vtxfile = "zvtxFastpp.dat";
  if(spd)    vtxfile = "zvtxSPD.dat";
  if(tracks) vtxfile = "zvtxTracks.dat";

  ReadZvtx(vtxfile,zvtx,n);


  TFile *outFile = TFile::Open(outName,"recreate");
  TFile *inFile  = TFile::Open(inName);
  TFile *inFile2 = TFile::Open(inName2);
  
  AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
  if(!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}
  
  Int_t flag1stPass,flag2ndPass;

  for(Int_t ev=0; ev<n; ev++){
    if(skipEvt[ev]) continue;
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;
    AliITStrackerV2 tracker(geom,ev);

    // set position of primary vertex
    Double_t vtx[3];
    vtx[0]=0.; vtx[1]=0.; vtx[2]=zvtx[ev];

    flag1stPass=1; // vtx constraint
    flag2ndPass=0; // no vtx constraint

    // no vtx constraint if vertex not found
    if(vtx[2]<-999.) {
      flag1stPass=0;
      vtx[2]=0.;
    }

    cerr<<"+++\n+++ Setting primary vertex z = "<<vtx[2]<<
      " cm for ITS tracking\n+++\n";
    tracker.SetVertex(vtx);  

    // setup vertex constraint in the two tracking passes
    Int_t flags[2];
    flags[0]=flag1stPass;
    tracker.SetupFirstPass(flags);
    flags[0]=flag2ndPass;
    tracker.SetupSecondPass(flags);
    
    tracker.Clusters2Tracks(inFile,outFile);
  }

  inFile->Close();
  inFile2->Close();
  outFile->Close();
 
  delete [] zvtx;

  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  
  return 0;
}
//-----------------------------------------------------------------------------
Int_t ITSMakeRefFile(const Char_t *galice,const Char_t *inname, 
		     const Char_t *outname,Bool_t *skipEvt,Int_t n) {

 
  cerr<<"\n*******************************************************************\n";

  Int_t rc=0;
  const Char_t *name="ITSMakeRefFile";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
  
  
  TFile *out = TFile::Open(outname,"recreate");
  TFile *trk = TFile::Open(inname);
  TFile *kin = TFile::Open(galice);

  
  // Get gAlice object from file
  if(!(gAlice=(AliRun*)kin->Get("gAlice"))) {
    cerr<<"gAlice has not been found on galice.root !\n";
    return 1;
  }
  
  Int_t label;
  TParticle *Part;  
  TParticle *Mum;
  static RECTRACK rectrk;
  

  for(Int_t ev=0; ev<n; ev++){
    if(skipEvt[ev]) continue;
    cerr<<" --- Processing event "<<ev<<" ---"<<endl;

    gAlice->GetEvent(ev);  

    trk->cd();

    // Tree with ITS tracks
    char tname[100];
    sprintf(tname,"TreeT_ITS_%d",ev);

    TTree *tracktree=(TTree*)trk->Get(tname);
    if(!tracktree) continue;
    AliITStrackV2 *itstrack=new AliITStrackV2; 
    tracktree->SetBranchAddress("tracks",&itstrack);
    Int_t nentr=(Int_t)tracktree->GetEntries();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",ev);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&rectrk,"lab/I:pdg:mumlab:mumpdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for(Int_t i=0; i<nentr; i++) {
      tracktree->GetEvent(i);
      label = TMath::Abs(itstrack->GetLabel());

      Part = (TParticle*)gAlice->Particle(label);
      rectrk.lab=label;
      rectrk.pdg=Part->GetPdgCode();
      rectrk.mumlab = Part->GetFirstMother();
      if(Part->GetFirstMother()>=0) {
	Mum = (TParticle*)gAlice->Particle(Part->GetFirstMother());
	rectrk.mumpdg=Mum->GetPdgCode();
      } else {
	rectrk.mumpdg=-1;
      }
      rectrk.Vx=Part->Vx();
      rectrk.Vy=Part->Vy();
      rectrk.Vz=Part->Vz();
      rectrk.Px=Part->Px();
      rectrk.Py=Part->Py();
      rectrk.Pz=Part->Pz();
      
      reftree->Fill();
    } // loop on tracks   

    out->cd();
    reftree->Write();

    delete itstrack;
    delete reftree;
  } // loop on events

  trk->Close();
  kin->Close();
  out->Close();
  
  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  

  return rc;
}
//-----------------------------------------------------------------------------
void WriteZvtx(const Char_t *name,Double_t *zvtx,Int_t n) {

  FILE *f = fopen(name,"w");
  for(Int_t i=0;i<n;i++) {
    fprintf(f,"%f\n",zvtx[i]);
  }
  fclose(f);

  return;
}
//-----------------------------------------------------------------------------
void ReadZvtx(const Char_t *name,Double_t *zvtx,Int_t n) {

  FILE *f = fopen(name,"r");
  for(Int_t i=0;i<n;i++) {
    fscanf(f,"%lf",&zvtx[i]);
  }
  fclose(f);

  return;
}










