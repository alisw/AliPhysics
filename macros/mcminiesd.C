// This is an example macro which loops on all ESD events in a file set
// and stores the selected information in form of mini ESD tree.
// It also extracts the corresponding MC information in a parallel MC tree.
//
// Input: galice.root,Kinematics.root,AliESDs.root (read only)
// Output: miniesd.root (recreate)
//
// Usage: faster if compiled
// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT/TPC");  
// gROOT->LoadMacro("mcminiesd.C+");
// mcminiesd()
//
// Author: Peter Hristov
// e-mail: Peter.Hristov@cern.ch

#if !defined(__CINT__) || defined(__MAKECINT__)

// Root include files
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TStopwatch.h>
#include <TObject.h>
#include <TParticle.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TArrayF.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>

// AliRoot include files
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliStack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliTPCtrack.h"
#include "AliTracker.h"

#endif

void selectRecTracks(AliESD* esdOut, Int_t * select, Int_t &nkeep0) {
  // Select also all the reconstructed tracks in case it was not done before
  Int_t nrec = esdOut->GetNumberOfTracks();
  for(Int_t irec=0; irec<nrec; irec++) {
    AliESDtrack * track = esdOut->GetTrack(irec);
    UInt_t label = TMath::Abs(track->GetLabel());
    if (label>=10000000) continue; // Underlying event. 10000000 is the
                                   // value of fkMASKSTEP in AliRunDigitizer
    if (!select[label]) {
      select[label] = 1;
      nkeep0++;
    }
  }
}

void copyGeneralESDInfo(AliESD* esdIn, AliESD* esdOut) {
  // Event and run number
  esdOut->SetEventNumber(esdIn->GetEventNumber());
  esdOut->SetRunNumber(esdIn->GetRunNumber());
  
  // Trigger
  esdOut->SetTriggerMask(esdIn->GetTriggerMask());

  // Magnetic field
  esdOut->SetMagneticField(esdIn->GetMagneticField());

  // Copy ESD vertex
  const AliESDVertex * vtxIn = esdIn->GetVertex();
  AliESDVertex * vtxOut = 0x0;
  if (vtxIn) {
    Double_t pos[3];
    vtxIn->GetXYZ(pos);
    Double_t cov[6];
    vtxIn->GetCovMatrix(cov);
  
    vtxOut = new AliESDVertex(pos,cov,
					     vtxIn->GetChi2(),
					     vtxIn->GetNContributors());
    Double_t tp[3];
    vtxIn->GetTruePos(tp);
    vtxOut->SetTruePos(tp);
  }
  else
    vtxOut = new AliESDVertex();
  
  esdOut->SetVertex(vtxOut);
}

void selectMiniESD(AliESD* esdIn, AliESD* &esdOut) {
  // This function copies the general ESD information
  // and selects the reconstructed tracks of interest
  // The second argument is a reference to pointer since we 
  // want to operate not with the content, but with the pointer itself

  printf("--------------------\n");
  printf("Selecting data mini ESD\n");
  TStopwatch timer;
  timer.Start();

  // Create the new output ESD
  esdOut = new AliESD();

  // Copy the general information
  copyGeneralESDInfo(esdIn, esdOut);

  // Select tracks
  Int_t ntrk = esdIn->GetNumberOfTracks();
  
  // Loop on tracks
  for (Int_t itrk=0; itrk<ntrk; itrk++) {
    
    AliESDtrack * trackIn = esdIn->GetTrack(itrk);
    UInt_t status=trackIn->GetStatus();

    //select only tracks with TPC or ITS refit
    if ((status&AliESDtrack::kTPCrefit)==0
	&& (status&AliESDtrack::kITSrefit)==0
	) continue;

    //select only tracks with "combined" PID
    if ((status&AliESDtrack::kESDpid)==0) continue;
    
    AliESDtrack * trackOut = new AliESDtrack(*trackIn);
    
    // Reset most of the redundant information
    trackOut->MakeMiniESDtrack();
    
    esdOut->AddTrack(trackOut);
    
    // Do not delete trackIn, it is a second pointer to an
    // object belonging to the TClonesArray
    delete trackOut;
  }

  printf("%d tracks selected\n",esdOut->GetNumberOfTracks());
  esdOut->Print();
  timer.Stop();
  timer.Print();
  printf("--------------------\n");
}

void fixMotherDaughter(TClonesArray& particles, Int_t * map) {
  // Fix mother-daughter relationship
  // using the map of old to new indexes
  Int_t nkeep = particles.GetEntriesFast();
  for (Int_t i=0; i<nkeep; i++) {
    TParticle * part = (TParticle*)particles[i];
    
    // mother
    Int_t mum = part->GetFirstMother();      
    if (mum>-1 && i>map[mum]) 
      part->SetFirstMother(map[mum]);
    
    // old indexes
    Int_t fd = part->GetFirstDaughter();
    Int_t ld = part->GetLastDaughter();
    
    // invalidate daughters
    part->SetFirstDaughter(-1);
    part->SetLastDaughter(-1);
    
    for (Int_t id=TMath::Max(fd,0); id<=ld; id++) {
      if (map[id]>-1) {
	// this daughter was selected
	if(part->GetFirstDaughter() < 0) part->SetFirstDaughter(map[id]);
	part->SetLastDaughter(map[id]);
      }
    }
  }
}

void checkConsistency(TClonesArray& particles) {
  // Check consistency
  // each mother is before the daughters,
  // each daughter from the mother's list 
  // points back to its mother
  Int_t nkeep = particles.GetEntriesFast();
  for (Int_t i=0; i<nkeep; i++) {
    TParticle * part = (TParticle*)particles[i];
    
    Int_t mum = part->GetFirstMother();
    
    if (mum>-1 && i<mum) {
      cout << "Problem: mother " << mum << " after daughter " << i << endl;
    }
    
    Int_t fd = part->GetFirstDaughter();
    Int_t ld = part->GetLastDaughter();
    
    if (fd > ld ) cout << "Problem " << fd << " > " << ld << endl;
    if (fd > -1 && fd < i ) 
      cout << "Problem: mother " << i << " after daughter " << ld << endl;
    
    for (Int_t id=TMath::Max(fd,0); id<=ld; id++) {
      TParticle * daughter = (TParticle*)particles[id];
      if (daughter->GetFirstMother() != i) {
	cout << "Problem "<< i << " not " 
	     << daughter->GetFirstMother()  << endl;
	daughter->Print();
      }
    }
    
  }
}


void kinetree(AliRunLoader* runLoader, AliESD* &esdOut, TClonesArray* &ministack) {

  printf("--------------------\n");
  printf("Selecting MC mini ESD\n");
  TStopwatch timer;
  timer.Start();

  // Particle properties
  TDatabasePDG * pdgdb = TDatabasePDG::Instance();

  // Particle stack
  AliStack * stack = runLoader->Stack();

  Int_t npart = stack->GetNtrack();

  // Particle momentum and vertex. Will be extracted from TParticle
  TLorentzVector momentum;
  TArrayF vertex(3);
  runLoader->GetHeader()->GenEventHeader()->PrimaryVertex(vertex);
  TVector3 dvertex;

  // Counter for selected particles
  Int_t nkeep0 = 0;

  Int_t * select = new Int_t[npart];

  // Loop on particles: physics selection
  for (Int_t ipart=0; ipart<npart; ipart++) {
      
    select[ipart] = 0;

    TParticle * part = stack->Particle(ipart);

    if (!part) {
      cout << "Empty particle " << ipart << endl;
      continue;
    }

    // Particle momentum and vertex of origin
    part->Momentum(momentum);
    dvertex.SetXYZ(part->Vx()-vertex[0],
		   part->Vy()-vertex[1],
		   part->Vz()-vertex[2]);

    // Now select only the "interesting" particles

    // Select all particles from the primary vertex:
    // resonances and products of resonance decays
    if(dvertex.Mag()<0.0001) {
      select[ipart] = 1;
      nkeep0++;
      continue;
    }

    // Reject particles born outside ITS
    if (dvertex.Perp()>100) continue;

    // Reject particles with too low momentum, they don't rich TPC
    if (part->Pt()<0.1) continue;

    // Reject "final state" neutral particles
    // This has to be redone after this loop
    Int_t pdgcode = part->GetPdgCode();
    TParticlePDG * pdgpart = pdgdb->GetParticle(pdgcode);
    Int_t charge  = 0;
    if (pdgpart) charge = (Int_t) pdgpart->Charge();

    if (!charge) {
      if (part->GetFirstDaughter()<0) continue;
    }
      
    // Select the rest of particles except the case
    // particle -> particle + X
    // for example bremstrahlung and delta electrons
    Int_t mumid = part->GetFirstMother();
    TParticle * mother = stack->Particle(mumid);
    if (mother) {
      Int_t mumpdg = mother->GetPdgCode();
      Bool_t skip = kFALSE;
      for (Int_t id=mother->GetFirstDaughter();
	   id<=mother->GetLastDaughter();
	   id++) {
	TParticle * daughter=stack->Particle(id);
	if (!daughter) continue;
	if (mumpdg == daughter->GetPdgCode()) {
	  skip=kTRUE;
	  break;
	}
      }
      if (skip) continue;
    }
    
    // All the criteria are OK, this particle is selected
    select[ipart] = 1;
    nkeep0++;

  } // end loop over particles in the current event

  selectRecTracks(esdOut, select, nkeep0);

  // Container for the selected particles
  TClonesArray * pparticles = new TClonesArray("TParticle",nkeep0);
  TClonesArray &particles = *pparticles;

  // Map of the indexes in the stack and the new ones in the TClonesArray
  Int_t * map = new Int_t[npart];

  // Counter for selected particles
  Int_t nkeep = 0;

  for (Int_t ipart=0; ipart<npart; ipart++) {
      
    map[ipart] = -99; // Reset the map
      
    if (select[ipart]) {

      TParticle * part = stack->Particle(ipart);
      map[ipart] = nkeep;
      new(particles[nkeep++]) TParticle(*part);

    }
  }

  if (nkeep0 != nkeep) printf("Strange %d is not %d\n",nkeep0,nkeep);

  delete [] select;

  // Fix mother-daughter relationship
  fixMotherDaughter(particles,map);

  // Check consistency
  checkConsistency(particles);

  // Now remove the "final state" neutral particles

  TClonesArray * particles1 = new TClonesArray("TParticle",nkeep);
  Int_t * map1 = new Int_t[nkeep];
  Int_t nkeep1 = 0;
    
  // Loop on particles
  for (Int_t ipart=0; ipart<nkeep; ipart++) {
      
    map1[ipart] = -99; // Reset the map

    TParticle * part = (TParticle*)particles[ipart];
    
    // Reject "final state" neutral particles
    // This has to be redone after this loop
    Int_t pdgcode = part->GetPdgCode();
    TParticlePDG * pdgpart = pdgdb->GetParticle(pdgcode);
    Int_t charge  =  0;
    if (pdgpart) charge = (Int_t) pdgpart->Charge();
    
    if (!charge) {
      if (part->GetFirstDaughter()<0) continue;
    }

    // Select the particle
    map1[ipart] = nkeep1;
    TClonesArray &ar = *particles1;
    new(ar[nkeep1++]) TParticle(*part);
  }

  particles.Delete();
  delete pparticles;
  cout << nkeep1 << " particles finally selected" << endl;

  fixMotherDaughter(*particles1,map1);
  checkConsistency(*particles1);

  // Remap the labels of reconstructed tracks

  AliESD * esdNew = new AliESD();

  copyGeneralESDInfo(esdOut,esdNew);

  // Tracks
  Int_t nrec = esdOut->GetNumberOfTracks();
  for(Int_t irec=0; irec<nrec; irec++) {
    AliESDtrack * track = esdOut->GetTrack(irec);
    UInt_t label = TMath::Abs(track->GetLabel());
    if (label<10000000) { // Signal event. 10000000 is the
                          // value of fkMASKSTEP in AliRunDigitizer
      track->SetLabel(TMath::Sign(map1[map[label]],track->GetLabel()));
    }
    esdNew->AddTrack(track);
  }

  delete esdOut;
  esdOut = esdNew;

  ministack = particles1;


  delete [] map;
  delete [] map1;


  timer.Stop();
  timer.Print();
  printf("--------------------\n");
}



AliTPCtrack * particle2track(TParticle * part) {

  // Converts TParticle to AliTPCtrack

  UInt_t index;
  Double_t xx[5];
  Double_t cc[15];
  Double_t xref;
  Double_t alpha;

  // Calculate alpha: the rotation angle of the corresponding TPC sector
  alpha = part->Phi()*180./TMath::Pi();
  if (alpha<0) alpha+= 360.;
  if (alpha>360) alpha -= 360.;

  Int_t sector = (Int_t)(alpha/20.);
  alpha = 10. + 20.*sector;
  alpha /= 180;
  alpha *= TMath::Pi();


  // Reset the covariance matrix
  for (Int_t i=0; i<15; i++) cc[i]=0.;


  // Index
  index = part->GetUniqueID();

  // Get the vertex of origin and the momentum
  TVector3 ver(part->Vx(),part->Vy(),part->Vz());
  TVector3 mom(part->Px(),part->Py(),part->Pz());


  // Rotate to the local (sector) coordinate system
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);

  // X of the referense plane
  xref = ver.X();

  // Track parameters
  // fX     = xref         X-coordinate of this track (reference plane)
  // fAlpha = Alpha        Rotation angle the local (TPC sector) 
  // fP0    = YL           Y-coordinate of a track
  // fP1    = ZG           Z-coordinate of a track
  // fP2    = C*x0         x0 is center x in rotated frame
  // fP3    = Tgl          tangent of the track momentum dip angle
  // fP4    = C            track curvature

  // Magnetic field
  TVector3 field(0.,0.,1/AliKalmanTrack::GetConvConst());
  // radius [cm] of track projection in (x,y) 
  Double_t rho =
    mom.Pt()*AliKalmanTrack::GetConvConst();

  Double_t charge = 
    TDatabasePDG::Instance()->GetParticle(part->GetPdgCode())->Charge();
  charge /=3.;

  if (TMath::Abs(charge) < 0.9) charge=1.e-9;
 
  TVector3 vrho = mom.Cross(field);
  vrho *= charge;
  vrho.SetMag(rho);

  vrho += ver;
  Double_t x0 = vrho.X();

  xx[0] = ver.Y();
  xx[1] = ver.Z();
  xx[3] = mom.Pz()/mom.Pt();
  xx[4] = -charge/rho;
  xx[2] = xx[4]*x0;
  //  if (TMath::Abs(charge) < 0.9) xx[4] = 0;

  AliTPCtrack * tr = new AliTPCtrack(index, xx, cc, xref, alpha);
  tr->SetLabel(index);
  return tr;

}

void mcminiesd() {

  printf("====================\n");
  printf("Main program\n");
  TStopwatch timer;
  timer.Start();

  // Input "data" file, tree, and branch
  TFile * fIn = TFile::Open("AliESDs.root");
  TTree * tIn = (TTree*)fIn->Get("esdTree");
  TBranch * bIn = tIn->GetBranch("ESD");  

  AliESD * esdIn = 0; // The input ESD object is put here
  bIn->SetAddress(&esdIn);

  // Output file, trees, and branches.
  // Both "data" and MC are stored in the same file
  TFile * fOut   = TFile::Open("miniesd.root","recreate");
  TTree * tOut   = new TTree("esdTree","esdTree");
  TTree * tMC = new TTree("mcStackTree","mcStackTree");

  AliESD * esdOut = 0; // The newly created ESD object
  TBranch * bOut = tOut->Branch("ESD","AliESD",&esdOut);
  bOut->SetCompressionLevel(9);

  // TParticles are stored in TClonesArray
  // The corresponding branch is created using "dummy" TClonesArray
  TClonesArray * ministack = new TClonesArray("TParticle");
  TBranch * bMC = tMC->Branch("Stack",&ministack);
  ministack->Delete();
  delete ministack;
  bMC->SetCompressionLevel(9);

  // Event header
  AliHeader * headerIn = 0; //
  TBranch * bHeader = tMC->Branch("Header","AliHeader",&headerIn);
  bHeader->SetCompressionLevel(9);

 // Run loader
  AliRunLoader* runLoader = AliRunLoader::Open("galice.root");

  // gAlice
  runLoader->LoadgAlice();
  gAlice = runLoader->GetAliRun();

  // Magnetic field
  AliTracker::SetFieldMap(gAlice->Field(),1); // 1 means uniform magnetic field

  // Now load kinematics and event header
  runLoader->LoadKinematics();
  runLoader->LoadHeader();

  // Loop on events: check that MC and data contain the same number of events
  Int_t nevESD = bIn->GetEntries();
  Int_t nevMC = (Int_t)runLoader->GetNumberOfEvents();

  if (nevESD != nevMC)
    cout << "Inconsistent number of ESD("<<nevESD<<") and MC("<<nevMC<<") events" << endl;

  // Loop on events
  Int_t nev=TMath::Min(nevMC,nevESD);

  cout << nev << " events to process" << endl;

  for (Int_t iev=0; iev<nev; iev++) {
    cout << "Event " << iev << endl;
    // Get "data" event
    bIn->GetEntry(iev);

    // "Data" selection
    selectMiniESD(esdIn,esdOut);

    // Get MC event
    runLoader->GetEvent(iev);

    // MC selection
    kinetree(runLoader,esdOut,ministack);
    bMC->SetAddress(&ministack); // needed because ministack has changed

    // Header
    headerIn = runLoader->GetHeader();

    // Fill the trees
    tOut->Fill();
    tMC->Fill();
    delete esdOut;
    esdOut = 0x0;
    //    delete esdIn;
    //    esdIn = 0x0;
    ministack->Delete();
    delete ministack;

  }

  fOut = tOut->GetCurrentFile();
  fOut->Write();
  fOut->Close();

  fIn->Close();

  timer.Stop();
  timer.Print();
  printf("====================\n");
}
