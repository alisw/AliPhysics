/****************************************************************************
 * This macro is supposed to do reconstruction in the ITS via Kalman        *
 * tracker V2. The ITStracker is feeded with parametrized TPC tracks        * 
 *                                                                          *
 * It does the following steps:                                             *
 *             1) TPC tracking parameterization                             *
 *             2) ITS cluster finding V2 (via fast points !)                *
 *             3) ITS track finding V2                                      *
 *             4) Create a reference file with simulation info (p,PDG...)   *
 *                                                                          *
 * (Origin: A.Dainese, Padova, andrea.dainese@pd,infn.it                    * 
 *  from AliBarrelReconstruction.C I.Belikov, CERN, Jouri.Belikov@cern.ch)  *
 ****************************************************************************/

#ifndef __CINT__
  #include "alles.h"
  #include "AliMagF.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSRecPoint.h"
  #include "AliITSclusterV2.h"
  #include "AliITSsimulationFastPoints.h"
  #include "AliITStrackerV2.h"
  #include "AliTPCtrackerParam.h"
#endif

typedef struct {
  Int_t lab;
  Int_t pdg;
  Int_t mumpdg;
  Float_t Vx,Vy,Vz;
  Float_t Px,Py,Pz;
} RECTRACK;


Int_t TPCParamTracks(const Char_t *galice,const Char_t *outname,const Int_t coll,const Double_t Bfield,Int_t n);
Int_t ITSFindClusters(const Char_t *inname,const Char_t *outname,Int_t n);
Int_t ITSFindTracks(const Char_t *galice,const Char_t *inname,const Char_t *inname2,const Char_t *outname,Int_t n);
Int_t ITSMakeRefFile(const Char_t *galice, const Char_t *inname, const Char_t *outname, Int_t n);

Int_t AliBarrelRec_TPCparam(Int_t n=1) {
  const Char_t *TPCtrkNameS="AliTPCtracksParam.root";
  const Char_t *galiceName="galice.root";
  const Char_t *ITSclsName="AliITSclustersV2.root";
  const Char_t *ITStrkName="AliITStracksV2.root";
  const Char_t *ITSrefName="ITStracksRefFile.root";

  // set here the code for the type of collision (needed for TPC tracking
  // parameterization). available collisions:
  //
  // coll = 0 ->   PbPb6000 (HIJING with b<2fm) 
  const Int_t    collcode = 0;  
  // set here the value of the magnetic field
  const Double_t BfieldValue = 0.4;



  AliKalmanTrack::SetConvConst(100/0.299792458/BfieldValue);

  
  // ********** Build TPC tracks with parameterization *********** //
  if (TPCParamTracks(galiceName,TPCtrkNameS,collcode,BfieldValue,n)) {
    cerr<<"Failed to get TPC hits !\n";
    return 1;
  }

  
  // ********** Find ITS clusters *********** //
  if (ITSFindClusters(galiceName,ITSclsName,n)) {
    cerr<<"Failed to get ITS clusters !\n";
    return 1;
  } 
  

  // ********* Find ITS tracks *********** //
  if (ITSFindTracks(galiceName,TPCtrkNameS,ITSclsName,ITStrkName,n)) {
    cerr<<"Failed to get ITS tracks !\n";
    return 1;
  } 
  
  
  // ********* Make ITS tracks reference file *********** //
  if (ITSMakeRefFile(galiceName,ITStrkName,ITSrefName,n)) {
    cerr<<"Failed to get ITS tracks ref file!\n";
    return 1;
  } 
  

  return 0;
}


Int_t TPCParamTracks(const Char_t *galice, const Char_t *outname,
		     const Int_t coll,const Double_t Bfield,Int_t n) {

  cerr<<"\n*******************************************************************\n";
  Int_t rc;

  const Char_t *name="TPCParamTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);

  TFile *outfile=TFile::Open(outname,"recreate");
  TFile *infile =TFile::Open(galice);

  AliTPCtrackerParam tracker(coll,Bfield);
  rc = tracker.BuildTPCtracks(infile,outfile,n);

  delete gAlice; gAlice=0;

  infile->Close();
  outfile->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return rc;
}

Int_t ITSFindClusters(const Char_t *inname, const Char_t *outname, Int_t n) {

 
  cerr<<"\n*******************************************************************\n";

  Int_t rc=0;
  const Char_t *name="ITSFindClusters";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
  TFile *out=TFile::Open(outname,"recreate");
  TFile *in =TFile::Open(inname,"update");

  
  if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
    cerr<<"Can't get gAlice !\n";
    return 1;
  }
  
  
  AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
  if (!ITS) { cerr<<"Can't get the ITS !\n"; return 1;}
  AliITSgeom *geom=ITS->GetITSgeom();
  out->cd();   
  geom->Write();
  
  Int_t ev=0;
  for (ev = 0; ev<n; ev++){
    in->cd();   // !!!! MI directory must point to galice. - othervise problem with Tree -connection
    gAlice->GetEvent(ev);
    //gAlice->TreeR()->Reset();   //reset reconstructed tree
    
     
    TTree *pTree=gAlice->TreeR();
    if (!pTree){
      gAlice->MakeTree("R");
      pTree = gAlice->TreeR();
    }
    TBranch *branch=pTree->GetBranch("ITSRecPoints");
    if (!branch) {
      //if not reconstructed ITS branch do reconstruction 
      ITS->MakeBranch("R",0);
      //////////////// Taken from ITSHitsToFastPoints.C ///////////////////////
      for (Int_t i=0;i<3;i++) { 
	ITS->SetSimulationModel(i,new AliITSsimulationFastPoints()); 
      }
      Int_t nsignal=25;
      Int_t size=-1;
      Int_t bgr_ev=Int_t(ev/nsignal);
      ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
      ///////////////////////////////////////////////////////////////////////
      gAlice->GetEvent(ev);   //MI comment  - in HitsToFast... they reset treeR to 0 
      //they overwrite full reconstructed event ???? ... so lets connect TreeR one more
      //time
    }


     
    out->cd();
    TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
    char   cname[100];
    sprintf(cname,"TreeC_ITS_%d",ev);
    
    TTree *cTree=new TTree(cname,"ITS clusters");
    cTree->Branch("Clusters",&clusters);
     
    pTree=gAlice->TreeR();
    if (!pTree) { cerr<<"Can't get TreeR !\n"; return 1; }
    branch=pTree->GetBranch("ITSRecPoints");
    if (!branch) { cerr<<"Can't get ITSRecPoints branch !\n"; return 1;}
    TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
    branch->SetAddress(&points);
     
    TClonesArray &cl=*clusters;
    Int_t nclusters=0;
    Int_t nentr=(Int_t)pTree->GetEntries();
    AliITSgeom *geom=ITS->GetITSgeom();

    for (Int_t i=0; i<nentr; i++) {
      if (!pTree->GetEvent(i)) {cTree->Fill(); continue;}
      Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
      Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
      Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
      Double_t yshift = x*rot[0] + y*rot[1];
      Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
      Int_t ncl=points->GetEntriesFast();
      nclusters+=ncl;
      Float_t lp[5];
      Int_t lab[6]; 
      for (Int_t j=0; j<ncl; j++) {
	AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
	lp[0]=-p->GetX()-yshift; if (lay==1) lp[0]=-lp[0];
	lp[1]=p->GetZ()+zshift;
	lp[2]=p->GetSigmaX2();
	lp[3]=p->GetSigmaZ2();
	lp[4]=p->GetQ();
	lab[0]=p->GetLabel(0);lab[1]=p->GetLabel(1);lab[2]=p->GetLabel(2);
	lab[3]=ndet;
	
	Int_t label=lab[0];
	TParticle *part=(TParticle*)gAlice->Particle(label);
	label=-3;
	while (part->P() < 0.005) {
	  Int_t m=part->GetFirstMother();
	  if (m<0) {cerr<<"Primary momentum: "<<part->P()<<endl; break;}
	  label=m;
	  part=(TParticle*)gAlice->Particle(label);
	}
	if      (lab[1]<0) lab[1]=label;
	else if (lab[2]<0) lab[2]=label;
	else cerr<<"No empty labels !\n";
	
	new(cl[j]) AliITSclusterV2(lab,lp);
      }
      cTree->Fill(); clusters->Delete();
      points->Delete();
    }
    cTree->Write();
    cerr<<"Number of clusters: "<<nclusters<<endl;
    delete cTree; delete clusters; delete points;
    
  }

  
  delete gAlice; gAlice=0;
  in->Close();
  out->Close();
  gBenchmark->Stop(name);
  gBenchmark->Show(name);

  return rc;
}

Int_t ITSFindTracks(const Char_t *galice, const Char_t * inname, 
                    const Char_t *inname2, const Char_t *outname, 
                    Int_t n) {

 
  cerr<<"\n*******************************************************************\n";

  Int_t rc=0;
  const Char_t *name="ITSFindTracks";
  cerr<<'\n'<<name<<"...\n";
  gBenchmark->Start(name);
  
  
  TFile *out=TFile::Open(outname,"recreate");
  TFile *in =TFile::Open(inname);
  TFile *in2 =TFile::Open(inname2);
  
  AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
  if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}

  
  for (Int_t ev=0; ev<n; ev++){
    AliITStrackerV2 tracker(geom,ev);
    
    TArrayF o;
    o.Set(3);
    TFile* vtxFile = new TFile(galice);
    vtxFile->cd();
    AliHeader* header = 0;
    TTree* treeE = (TTree*)gDirectory->Get("TE");
    treeE->SetBranchAddress("Header",&header);
    treeE->GetEntry(ev);
    AliGenEventHeader* genHeader = header->GenEventHeader();
    if(genHeader) {
      // get primary vertex position
      genHeader->PrimaryVertex(o);
      vtxFile->Close();
      delete header;
    // set position of primary vertex
      Double_t vtx[3];
      vtx[0]=o[0]; vtx[1]=o[1]; vtx[2]=o[2];
      cerr<<"+++\n+++ Reading primary vertex position from galice.root\n+++\n+++ Setting primary vertex z = "<<vtx[2]<<" cm for ITS tracking\n+++\n";
      tracker.SetVertex(vtx);  
    }else {
      cerr<<"+++\n+++ Event header not found in galice.root:\n+++ Primary vertex in (0,0,0) [default]\n+++\n";
    }    

    // setup vertex constraint in the two tracking passes
    Int_t flags[2];
    flags[0]=0;
    tracker.SetupFirstPass(flags);
    flags[0]=-1;
    tracker.SetupSecondPass(flags);
    
    rc=tracker.Clusters2Tracks(in,out);
 
  }

  delete gAlice;  gAlice=0;
  in->Close();
  in2->Close();
  out->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  
  return rc;
}


Int_t ITSMakeRefFile(const Char_t *galice, const Char_t *inname, 
		     const Char_t *outname, Int_t n) {

 
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
  

  for(Int_t event=0; event<n; event++){

    AliITStrackV2 *itstrack=0;

    Int_t nparticles=gAlice->GetEvent(event);  

    trk->cd();

    // Tree with ITS tracks
    char tname[100];
    sprintf(tname,"TreeT_ITS_%d",event);

    TTree *tracktree=(TTree*)trk->Get(tname);
    TBranch *tbranch=tracktree->GetBranch("tracks");
    Int_t nentr=(Int_t)tracktree->GetEntries();

    // Tree for true track parameters
    char ttname[100];
    sprintf(ttname,"Tree_Ref_%d",event);
    TTree *reftree = new TTree(ttname,"Tree with true track params");
    reftree->Branch("rectracks",&rectrk,"lab/I:pdg:Vx/F:Vy:Vz:Px:Py:Pz");

    for (Int_t i=0; i<nentr; i++) {
      itstrack=new AliITStrackV2;
      tbranch->SetAddress(&itstrack);
      tracktree->GetEvent(i);
      label = TMath::Abs(itstrack->GetLabel());

      Part = (TParticle*)gAlice->Particle(label);
      rectrk.lab=label;
      rectrk.pdg=Part->GetPdgCode();
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
    }   

    out->cd();
    reftree->Write();

  }

  trk->Close();
  kin->Close();
  out->Close();

  gBenchmark->Stop(name);
  gBenchmark->Show(name);
  

  return rc;

}












