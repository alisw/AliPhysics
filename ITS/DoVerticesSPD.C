#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TFile.h"
#include "TGeoManager.h"
#include "TInterpreter.h"
#include "TClassTable.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TArrayF.h"
#include "TParticle.h"
#include "AliESD.h"
#include "AliGenPythiaEventHeader.h"
#include "AliTracker.h"
#include "AliHeader.h"
#include "AliITSLoader.h"
#include "AliITSsegmentationSPD.h"
#include "AliVertexerTracks.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliGRPManager.h"
#include "AliITSDetTypeRec.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#endif

/*  $Id$    */

Bool_t DoVerticesSPD(Bool_t isMC=kFALSE, Int_t pileupalgo=1, Int_t optdebug=0){

  TFile *fint = new TFile("VertexSPD.root","recreate");
  TNtuple *nt = new TNtuple("ntvert","vertices","xtrue:ytrue:ztrue:zZ:zdiffZ:zerrZ:ntrksZ:x3D:xdiff3D:xerr3D:y3D:ydiff3D:yerr3D:z3D:zdiff3D:zerr3D:ntrks3D:dndy:ntrklets:nrp1:ptyp:is3D:isTriggered");

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
    //    man->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
    man->SetRun(0);
  }
  else {
    printf("Using deafult storage \n");
  }
 
  // retrives geometry 
  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  }
  AliGeomManager::ApplyAlignObjsFromCDB("ITS");

  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    Error("DoVertices", "getting run loader from file %s failed", 
	  "galice.root");
    return kFALSE;
  }
  
  if(isMC){
    runLoader->LoadgAlice();
    gAlice = runLoader->GetAliRun();
    if (!gAlice) {
      Error("DoVertices", "no galice object found");
      return kFALSE;
    }
    runLoader->LoadKinematics();
    runLoader->LoadHeader();
  }
  AliITSLoader* ITSloader =  (AliITSLoader*) runLoader->GetLoader("ITSLoader");
  ITSloader->LoadRecPoints("read");

  Int_t totev=runLoader->GetNumberOfEvents();
  if(optdebug)  printf("Number of events= %d\n",totev);

  TFile* esdFile = TFile::Open("AliESDs.root");
  if (!esdFile || !esdFile->IsOpen()) {
    printf("Error in opening ESD file");
    return kFALSE;
  }

  AliESDEvent * esd = new AliESDEvent;
  TTree* tree = (TTree*) esdFile->Get("esdTree");
  if (!tree) {
    printf("Error: no ESD tree found");
    return kFALSE;
  }
  esd->ReadFromTree(tree);
  AliGRPManager * grpMan=new AliGRPManager();
  grpMan->ReadGRPEntry();
  grpMan->SetMagField();
  printf("Magnetic field set to %f T\n",AliTracker::GetBz()/10.);

  AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec();
  AliITSsegmentation* seg = new AliITSsegmentationSPD();  
  detTypeRec->SetSegmentationModel(0,seg);

  Double_t xnom=0.,ynom=0.;
  AliITSVertexerZ *vertz = new AliITSVertexerZ(xnom,ynom);
  vertz->Init("default");
  AliITSVertexer3D *vert3d = new AliITSVertexer3D();  
  vert3d->Init("default");
  vert3d->SetWideFiducialRegion(40.,1.);
  vert3d->SetPileupAlgo(pileupalgo);
  vert3d->PrintStatus();
  vertz->SetDetTypeRec(detTypeRec);
  vert3d->SetDetTypeRec(detTypeRec);
  vert3d->SetComputeMultiplicity(kTRUE);
  /* uncomment these lines to use diamond constrain */
//   Double_t posdiam[3]={0.03,0.1,0.};
//   Double_t sigdiam[3]={0.01,0.01,10.0};
//   AliESDVertex* diam=new AliESDVertex(posdiam,sigdiam);
//   vertz->SetVtxStart(diam);
//   vert3d->SetVtxStart(diam);
  /* end lines to be uncommented to use diamond constrain */

  Int_t goodz=0,good3d=0;
  // Trigger mask
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
    TArrayF mcVertex(3); 
    runLoader->GetEvent(iEvent);
    Double_t dNchdy = 0.;
    Int_t ptype = 0;
    if(optdebug){
      printf("==============================================================\n");
      printf("Event: %d \n",iEvent);
    }
    if(isMC){
      runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
      AliGenPythiaEventHeader *evh=(AliGenPythiaEventHeader*)runLoader->GetHeader()->GenEventHeader();
      ptype = evh->ProcessType();

      AliStack* stack = runLoader->Stack();
      TTree *treek=(TTree*)runLoader->TreeK();
      Int_t npart = (Int_t)treek->GetEntries();
      if(optdebug) printf("Process Type = %d --- Particles  %d\n",ptype,npart);
      
      // loop on particles to get generated dN/dy
      for(Int_t iPart=0; iPart<npart; iPart++) {
	if(!stack->IsPhysicalPrimary(iPart)) continue;
	TParticle* part = (TParticle*)stack->Particle(iPart);
	if(part->GetPDG()->Charge() == 0) continue;
	Double_t eta=part->Eta();
	
	if(TMath::Abs(eta)<1.5) dNchdy+=1.; 
      }
      if(optdebug) printf(" dNch/dy = %f\n",dNchdy);
    }
    tree->GetEvent(iEvent);

    TTree* cltree = ITSloader->TreeR();
    ULong64_t triggerMask=esd->GetTriggerMask();
    // MB1: SPDFO || V0L || V0R
    Bool_t eventTriggered = (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)));
    AliESDVertex* vtxz = vertz->FindVertexForCurrentEvent(cltree);
    AliESDVertex* vtx3d = vert3d->FindVertexForCurrentEvent(cltree);
    AliMultiplicity *alimult = vert3d->GetMultiplicity();
    Int_t  ntrklets=0;
    Int_t nrecp1=0;
    if(alimult){
      ntrklets=alimult->GetNumberOfTracklets() ;
      nrecp1=ntrklets+alimult->GetNumberOfSingleClusters();
    }


    TDirectory *current = gDirectory;
    fint->cd();
    Float_t xnt[23];

    Double_t zz = 0.;
    Double_t zresz = 0.;
    Double_t zdiffz = 0.;

    Int_t ntrkz = 0;
    if(vtxz){
 
      ntrkz = vtxz->GetNContributors();
      if(ntrkz>0)goodz++;
      zz=vtxz->GetZv();
      zresz =vtxz->GetZRes(); // microns
      zdiffz = 10000.*(zz-mcVertex[2]); // microns
    }

    Double_t x3d = 0.;
    Double_t xerr3d = 0.;
    Double_t xdiff3d = 0.;

    Double_t y3d = 0.;
    Double_t yerr3d = 0.;
    Double_t ydiff3d = 0.;

    Double_t z3d = 0.;
    Double_t zerr3d = 0.;
    Double_t zdiff3d = 0.;

    Int_t ntrk3d = 0;
    Bool_t is3d=kFALSE;
    if(vtx3d){

      if(vtx3d->IsFromVertexer3D()) is3d=kTRUE;
      ntrk3d = vtx3d->GetNContributors();
      if(is3d && ntrk3d>0)good3d++;
      x3d = vtx3d->GetXv();
      xerr3d=vtx3d->GetXRes();
      xdiff3d = 10000.*(x3d-mcVertex[0]);  // microns

      y3d = vtx3d->GetYv();
      yerr3d=vtx3d->GetYRes();
      ydiff3d = 10000.*(y3d-mcVertex[1]);  // microns

      z3d = vtx3d->GetZv();
      zerr3d=vtx3d->GetZRes();
      zdiff3d = 10000.*(z3d-mcVertex[2]);  // microns

    }

    xnt[0]=mcVertex[0];//x
    xnt[1]=mcVertex[1];//y
    xnt[2]=mcVertex[2];//z

    xnt[3]=zz;
    xnt[4]=zdiffz;
    xnt[5]=zresz;
    xnt[6]=ntrkz;

    xnt[7]=x3d;
    xnt[8]=xdiff3d;
    xnt[9]=xerr3d;
    xnt[10]=y3d;
    xnt[11]=ydiff3d;
    xnt[12]=yerr3d;
    xnt[13]=z3d;
    xnt[14]=zdiff3d;
    xnt[15]=zerr3d;
    xnt[16]=ntrk3d;

    xnt[17]=dNchdy;
    xnt[18]=float(ntrklets);
    xnt[19]=float(nrecp1);
    xnt[20]=float(ptype);
    xnt[21]=float(is3d);
    xnt[22]=float(eventTriggered);
    nt->Fill(xnt);
    current->cd();
    
    if(optdebug){
      printf("Vertexer3D: \tx=%.5f \ty=%.5f \tz(cm)=%.5f \tContributors=%d \n",x3d,y3d,z3d,ntrk3d);
      printf("VertexerZ:  \tx=%.5f \ty=%.5f \tz(cm)=%.5f \tContributors=%d \n",0.,0.,zz,ntrkz);
      if(isMC) printf("True Pos.: \tx=%.5f \ty=%.5f \tz(cm)=%.5f \tdN/dy=%.1f\n",mcVertex[0],mcVertex[1],mcVertex[2],dNchdy);
      printf("Multiplicity: Tracklets=%d  ClustersLay1=%d\n",ntrklets,nrecp1);
    }
    
  }
  fint->cd();
  nt->Write();
  fint->Close();
  delete fint;
  if(optdebug){
   printf("***********************************************\n");
   printf("Number of Z vertices found= %d\n",goodz);
   printf("efficiency=%f\n",float(goodz)/float(totev));
   printf("Number of 3D vertices found= %d\n",good3d);
   printf("efficiency=%f\n",float(good3d)/float(totev));
  }
  return kTRUE;
}
