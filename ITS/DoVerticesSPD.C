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
#include "AliVertexerTracks.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"
#include "AliITSVertexer3D.h"
#include "AliITSVertexerZ.h"
#include "AliESDVertex.h"
#include "AliRun.h"
#include "AliRunLoader.h"
#include "AliGenEventHeader.h"
#include "AliStack.h"
#endif

/*  $Id$    */

Bool_t DoVerticesSPD(Int_t optdebug=1){
  TFile *fint = new TFile("VertexSPD.root","recreate");
  TNtuple *nt = new TNtuple("ntvert","vertices","xtrue:ytrue:ztrue:zZ:zdiffZ:zerrZ:ntrksZ:x3D:xdiff3D:xerr3D:y3D:ydiff3D:yerr3D:z3D:zdiff3D:zerr3D:ntrks3D:dndy:ntrklets:nrp1:ptyp");

  if (gClassTable->GetID("AliRun") < 0) {
    gInterpreter->ExecuteMacro("loadlibs.C");
  }
  // Set OCDB if needed
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man->IsDefaultStorageSet()) {
    printf("Setting a local default storage and run number 0\n");
    man->SetDefaultStorage("local://$ALICE_ROOT");
    man->SetRun(0);
  }
  else {
    printf("Using deafult storage \n");
  }
 
  // retrives geometry 
  if(!gGeoManager){
    AliGeomManager::LoadGeometry("geometry.root");
  }

  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");
  if (!runLoader) {
    Error("DoVertices", "getting run loader from file %s failed", 
	  "galice.root");
    return kFALSE;
  }
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();
  if (!gAlice) {
    Error("DoVertices", "no galice object found");
    return kFALSE;
  }
  runLoader->LoadKinematics();
  runLoader->LoadHeader();
  AliITSLoader* ITSloader =  (AliITSLoader*) runLoader->GetLoader("ITSLoader");
  ITSloader->LoadRecPoints("read");

  AliMagF * magf = gAlice->Field();
  AliTracker::SetFieldMap(magf,kTRUE);
  if(optdebug) printf("MagneticField=%f\n",AliTracker::GetBz());
  
  Int_t totev=runLoader->GetNumberOfEvents();
  if(optdebug)  printf("Number of events= %d\n",totev);

  
//   TFile* esdFile = TFile::Open("AliESDs.root");
//   if (!esdFile || !esdFile->IsOpen()) {
//     Error("DoVertices", "opening ESD file %s failed", "AliESDs.root");
//     return kFALSE;
//   }
//   AliESD* esd = new AliESD;
//   TTree* tree = (TTree*) esdFile->Get("esdTree");
//   if (!tree) {
//     Error("DoVertices", "no ESD tree found");
//     return kFALSE;
//   }
//   tree->SetBranchAddress("ESD", &esd);

  Double_t xnom=0.,ynom=0.;
  AliITSVertexerZ *vertz = new AliITSVertexerZ(xnom,ynom);
  vertz->Init("default");
  AliITSVertexer3D *vert3d = new AliITSVertexer3D();
  vert3d->Init("default");
  //  vert3d->SetDebug(10);
  //  vertz->ConfigIterations(5);

  Int_t goodz=0,good3d=0;

  for (Int_t iEvent = 0; iEvent < totev; iEvent++) {
    TArrayF mcVertex(3); 
    runLoader->GetEvent(iEvent);
    runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(mcVertex);
    AliGenPythiaEventHeader *evh=(AliGenPythiaEventHeader*)runLoader->GetHeader()->GenEventHeader();
    Int_t ptype = evh->ProcessType();
    if(optdebug){
      printf("==============================================================\n");
      printf("\nEvent: %d ---- Process Type = %d \n",iEvent,ptype);
    }

    AliStack* stack = runLoader->Stack();
    TTree *treek=(TTree*)runLoader->TreeK();
    Int_t npart = (Int_t)treek->GetEntries();
    if(optdebug) printf("particles  %d\n",npart);

    Double_t dNchdy = 0.;

   // loop on particles
    for(Int_t pa=0; pa<npart; pa++) {
      TParticle* part = (TParticle*)stack->Particle(pa);
      Int_t pdg = part->GetPdgCode();
      Int_t apdg = TMath::Abs(pdg);
      Double_t energy  = part->Energy();
      if(energy>6900.) continue; // reject incoming protons
      Double_t pz = part->Pz();
      Double_t y = 0.5*TMath::Log((energy+pz+1.e-13)/(energy-pz+1.e-13));


      if(apdg!=11 && apdg!=13 && apdg!=211 && apdg!=321 && apdg!=2212) continue;      // reject secondaries
      if(TMath::Sqrt((part->Vx()-mcVertex[0])*(part->Vx()-mcVertex[0])+(part->Vy()-mcVertex[1])*(part->Vy()-mcVertex[1]))>0.0010) continue;
      if(TMath::Abs(y)<1.0) dNchdy += 0.5; // count 1/2 of particles in |y|<1
    }
    if(optdebug) printf(" dNch/dy = %f\n",dNchdy);
 
    TTree* cltree = ITSloader->TreeR();

    AliESDVertex* vtxz = vertz->FindVertexForCurrentEvent(cltree);
    AliMultiplicity *alimult = vertz->GetMultiplicity();
    Int_t ntrklets=0,nrecp1=0;
    if(alimult) {
      nrecp1=alimult->GetNumberOfTracklets() ;
      ntrklets = alimult->GetNumberOfTracklets();
      for(Int_t l=0;l<alimult->GetNumberOfTracklets();l++){
	if(alimult->GetDeltaPhi(l)<-9998.) ntrklets--;
      }
    }

    AliESDVertex* vtx3d = vert3d->FindVertexForCurrentEvent(cltree);

    TDirectory *current = gDirectory;
    fint->cd();
    Float_t xnt[21];

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
    if(vtx3d){

      ntrk3d = vtx3d->GetNContributors();
      if(ntrk3d>0)good3d++;
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
    xnt[11]=y3d;
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
    nt->Fill(xnt);
    current->cd();
    
    if(optdebug){
      printf("\nVertexerZ:  \tz(cm)=%9.5f \tTracklets=%d \tztrue(cm)=%9.5f \tzdiff(um)=%8.2f\n",zz,ntrkz,mcVertex[2],zdiffz);
      printf("Vertexer3D: \tz(cm)=%9.5f \tTracklets=%d \tztrue(cm)=%9.5f \tzdiff(um)=%8.2f\n",z3d,ntrk3d,mcVertex[2],zdiff3d);
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
