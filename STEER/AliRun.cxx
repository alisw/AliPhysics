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

/*
$Log$
Revision 1.21  1999/11/25 10:40:08  fca
Fixing daughters information also in primary tracks

Revision 1.20  1999/10/04 18:08:49  fca
Adding protection against inconsistent Euclid files

Revision 1.19  1999/09/29 07:50:40  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Control class for Alice C++                                              //
//  Only one single instance of this class exists.                           //
//  The object is created in main program aliroot                            //
//  and is pointed by the global gAlice.                                     //
//                                                                           //
//   -Supports the list of all Alice Detectors (fModules).                 //
//   -Supports the list of particles (fParticles).                           //
//   -Supports the Trees.                                                    //
//   -Supports the geometry.                                                 //
//   -Supports the event display.                                            //
//Begin_Html
/*
<img src="picts/AliRunClass.gif">
*/
//End_Html
//Begin_Html
/*
<img src="picts/alirun.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TRandom.h>
#include <TBRIK.h> 
#include <TNode.h> 
#include <TCint.h> 
#include <TSystem.h>
#include <TObjectTable.h>

#include "TParticle.h"
#include "AliRun.h"
#include "AliDisplay.h"

#include "AliCallf77.h" 
 
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
 
AliRun *gAlice;

static AliHeader *header;

#ifndef WIN32

# define rxgtrak rxgtrak_
# define rxstrak rxstrak_
# define rxkeep  rxkeep_ 
# define rxouth  rxouth_
#else

# define rxgtrak RXGTRAK 
# define rxstrak RXSTRAK 
# define rxkeep  RXKEEP  
# define rxouth  RXOUTH
#endif

static TArrayF sEventEnergy;
static TArrayF sSummEnergy;
static TArrayF sSum2Energy;

ClassImp(AliRun)

//_____________________________________________________________________________
AliRun::AliRun()
{
  //
  // Default constructor for AliRun
  //
  header=&fHeader;
  fRun       = 0;
  fEvent     = 0;
  fCurrent   = -1;
  fModules = 0;
  fGenerator = 0;
  fTreeD     = 0;
  fTreeK     = 0;
  fTreeH     = 0;
  fTreeE     = 0;
  fTreeR     = 0;
  fParticles = 0;
  fGeometry  = 0;
  fDisplay   = 0;
  fField     = 0;
  fMC        = 0;
  fNdets     = 0;
  fImedia    = 0;
  fTrRmax    = 1.e10;
  fTrZmax    = 1.e10;
  fInitDone  = kFALSE;
  fLego      = 0;
  fPDGDB     = 0;        //Particle factory object!
}

//_____________________________________________________________________________
AliRun::AliRun(const char *name, const char *title)
  : TNamed(name,title)
{
  //
  //  Constructor for the main processor.
  //  Creates the geometry
  //  Creates the list of Detectors.
  //  Creates the list of particles.
  //
  Int_t i;
  
  gAlice     = this;
  fTreeD     = 0;
  fTreeK     = 0;
  fTreeH     = 0;
  fTreeE     = 0;
  fTreeR     = 0;
  fTrRmax    = 1.e10;
  fTrZmax    = 1.e10;
  fGenerator = 0;
  fInitDone  = kFALSE;
  fLego      = 0;
  fField     = 0;
  
  gROOT->GetListOfBrowsables()->Add(this,name);
  //
  // create the support list for the various Detectors
  fModules = new TObjArray(77);
  //
  // Create the TNode geometry for the event display
  
  BuildSimpleGeometry();
  

  fNtrack=0;
  fHgwmk=0;
  fCurrent=-1;
  header=&fHeader;
  fRun       = 0;
  fEvent     = 0;
  //
  // Create the particle stack
  fParticles = new TClonesArray("TParticle",100);
  
  fDisplay = 0;
  //
  // Create default mag field
  SetField();
  //
  fMC      = gMC;
  //
  // Prepare the tracking medium lists
  fImedia = new TArrayI(1000);
  for(i=0;i<1000;i++) (*fImedia)[i]=-99;
  //
  // Make particles
  fPDGDB     = TDatabasePDG::Instance();        //Particle factory object!
}

//_____________________________________________________________________________
AliRun::~AliRun()
{
  //
  // Defaullt AliRun destructor
  //
  delete fImedia;
  delete fField;
  delete fMC;
  delete fGeometry;
  delete fDisplay;
  delete fGenerator;
  delete fLego;
  delete fTreeD;
  delete fTreeK;
  delete fTreeH;
  delete fTreeE;
  delete fTreeR;
  if (fModules) {
    fModules->Delete();
    delete fModules;
  }
  if (fParticles) {
    fParticles->Delete();
    delete fParticles;
  }
}

//_____________________________________________________________________________
void AliRun::AddHit(Int_t id, Int_t track, Int_t *vol, Float_t *hits) const
{
  //
  //  Add a hit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) ((AliModule*) dets[id])->AddHit(track,vol,hits);
}

//_____________________________________________________________________________
void AliRun::AddDigit(Int_t id, Int_t *tracks, Int_t *digits) const
{
  //
  // Add digit to detector id
  //
  TObjArray &dets = *fModules;
  if(dets[id]) ((AliModule*) dets[id])->AddDigit(tracks,digits);
}

//_____________________________________________________________________________
void AliRun::Browse(TBrowser *b)
{
  //
  // Called when the item "Run" is clicked on the left pane
  // of the Root browser.
  // It displays the Root Trees and all detectors.
  //
  if (fTreeK) b->Add(fTreeK,fTreeK->GetName());
  if (fTreeH) b->Add(fTreeH,fTreeH->GetName());
  if (fTreeD) b->Add(fTreeD,fTreeD->GetName());
  if (fTreeE) b->Add(fTreeE,fTreeE->GetName());
  if (fTreeR) b->Add(fTreeR,fTreeR->GetName());
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    b->Add(detector,detector->GetName());
  }
}

//_____________________________________________________________________________
void AliRun::Build()
{
  //
  // Initialize Alice geometry
  // Dummy routine
  //
}
 
//_____________________________________________________________________________
void AliRun::BuildSimpleGeometry()
{
  //
  // Create a simple TNode geometry used by Root display engine
  //
  // Initialise geometry
  //
  fGeometry = new TGeometry("AliceGeom","Galice Geometry for Hits");
  new TMaterial("void","Vacuum",0,0,0);  //Everything is void
  TBRIK *brik = new TBRIK("S_alice","alice volume","void",2000,2000,3000);
  brik->SetVisibility(0);
  new TNode("alice","alice","S_alice");
}

//_____________________________________________________________________________
void AliRun::CleanDetectors()
{
  //
  // Clean Detectors at the end of event
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->FinishEvent();
  }
}

//_____________________________________________________________________________
void AliRun::CleanParents()
{
  //
  // Clean Particles stack.
  // Set parent/daughter relations
  //
  TClonesArray &particles = *(gAlice->Particles());
  TParticle *part;
  int i;
  for(i=0; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    if(!part->TestBit(Daughters_Bit)) {
      part->SetFirstDaughter(-1);
      part->SetLastDaughter(-1);
    }
  }
}

//_____________________________________________________________________________
Int_t AliRun::DistancetoPrimitive(Int_t, Int_t)
{
  //
  // Return the distance from the mouse to the AliRun object
  // Dummy routine
  //
  return 9999;
}

//_____________________________________________________________________________
void AliRun::DumpPart (Int_t i)
{
  //
  // Dumps particle i in the stack
  //
  TClonesArray &particles = *fParticles;
  ((TParticle*) particles[i])->Print();
}

//_____________________________________________________________________________
void AliRun::DumpPStack ()
{
  //
  // Dumps the particle stack
  //
  TClonesArray &particles = *fParticles;
  printf(
	 "\n\n=======================================================================\n");
  for (Int_t i=0;i<fNtrack;i++) 
    {
      printf("-> %d ",i); ((TParticle*) particles[i])->Print();
      printf("--------------------------------------------------------------\n");
    }
  printf(
	 "\n=======================================================================\n\n");
}

//_____________________________________________________________________________
void AliRun::SetField(Int_t type, Int_t version, Float_t scale,
		      Float_t maxField, char* filename)
{
  //
  //  Set magnetic field parameters
  //  type      Magnetic field transport flag 0=no field, 2=helix, 3=Runge Kutta
  //  version   Magnetic field map version (only 1 active now)
  //  scale     Scale factor for the magnetic field
  //  maxField  Maximum value for the magnetic field

  //
  // --- Sanity check on mag field flags
  if(type<0 || type > 2) {
    Warning("SetField",
	    "Invalid magnetic field flag: %5d; Helix tracking chosen instead\n"
	   ,type);
    type=2;
  }
  if(fField) delete fField;
  if(version==1) {
    fField = new AliMagFC("Map1"," ",type,version,scale,maxField);
  } else if(version<=3) {
    fField = new AliMagFCM("Map2-3",filename,type,version,scale,maxField);
    fField->ReadField();
  } else {
    Warning("SetField","Invalid map %d\n",version);
  }
}

//_____________________________________________________________________________
void AliRun::FillTree()
{
  //
  // Fills all AliRun TTrees
  //
  if (fTreeK) fTreeK->Fill();
  if (fTreeH) fTreeH->Fill();
  if (fTreeD) fTreeD->Fill();
  if (fTreeR) fTreeR->Fill();
}
 
//_____________________________________________________________________________
void AliRun::FinishPrimary()
{
  //
  // Called  at the end of each primary track
  //
  
  //  static Int_t count=0;
  //  const Int_t times=10;
  // This primary is finished, purify stack
  gAlice->PurifyKine();

  // Write out hits if any
  if (gAlice->TreeH()) {
    gAlice->TreeH()->Fill();
  }
  
  // Reset Hits info
  gAlice->ResetHits();

  //
  //  if(++count%times==1) gObjectTable->Print();
}

//_____________________________________________________________________________
void AliRun::FinishEvent()
{
  //
  // Called at the end of the event.
  //
  
  //Update the energy deposit tables
  Int_t i;
  for(i=0;i<sEventEnergy.GetSize();i++) {
    sSummEnergy[i]+=sEventEnergy[i];
    sSum2Energy[i]+=sEventEnergy[i]*sEventEnergy[i];
  }
  sEventEnergy.Reset();
  
  // Clean detector information
  CleanDetectors();
  
  // Write out the kinematics
  if (fTreeK) {
    CleanParents();
    fTreeK->Fill();
  }
  
  // Write out the digits
  if (fTreeD) {
    fTreeD->Fill();
    ResetDigits();
  }
  
  // Write out reconstructed clusters  
  if (fTreeR) {
    fTreeR->Fill();
  }

  // Write out the event Header information
  if (fTreeE) fTreeE->Fill();
  
  // Reset stack info
  ResetStack();
  
  // Write Tree headers
  //  Int_t ievent = fHeader.GetEvent();
  //  char hname[30];
  //  sprintf(hname,"TreeK%d",ievent);
  if (fTreeK) fTreeK->Write();
  //  sprintf(hname,"TreeH%d",ievent);
  if (fTreeH) fTreeH->Write();
  //  sprintf(hname,"TreeD%d",ievent);
  if (fTreeD) fTreeD->Write();
  //  sprintf(hname,"TreeR%d",ievent);
  if (fTreeR) fTreeR->Write();
}

//_____________________________________________________________________________
void AliRun::FinishRun()
{
  //
  // Called at the end of the run.
  //

  // Clean detector information
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->FinishRun();
  }
  
  //Output energy summary tables
  EnergySummary();
  
  // file is retrieved from whatever tree
  TFile *File = 0;
  if (fTreeK) File = fTreeK->GetCurrentFile();
  if ((!File) && (fTreeH)) File = fTreeH->GetCurrentFile();
  if ((!File) && (fTreeD)) File = fTreeD->GetCurrentFile();
  if ((!File) && (fTreeE)) File = fTreeE->GetCurrentFile();
  if( NULL==File ) {
    Error("FinishRun","There isn't root file!");
    exit(1);
  }
  File->cd();
  fTreeE->Write();
  
  // Clean tree information
  delete fTreeK; fTreeK = 0;
  delete fTreeH; fTreeH = 0;
  delete fTreeD; fTreeD = 0;
  delete fTreeR; fTreeR = 0;
  delete fTreeE; fTreeE = 0;
  
  // Write AliRun info and all detectors parameters
  Write();
  
  // Close output file
  File->Write();
  File->Close();
}

//_____________________________________________________________________________
void AliRun::FlagTrack(Int_t track)
{
  //
  // Flags a track and all its family tree to be kept
  //
  int curr;
  TParticle *particle;

  curr=track;
  while(1) {
    particle=(TParticle*)fParticles->UncheckedAt(curr);
    
    // If the particle is flagged the three from here upward is saved already
    if(particle->TestBit(Keep_Bit)) return;
    
    // Save this particle
    particle->SetBit(Keep_Bit);
    
    // Move to father if any
    if((curr=particle->GetFirstMother())==-1) return;
  }
}
 
//_____________________________________________________________________________
void AliRun::EnergySummary()
{
  //
  // Print summary of deposited energy
  //

  Int_t ndep=0;
  Float_t edtot=0;
  Float_t ed, ed2;
  Int_t kn, i, left, j, id;
  const Float_t zero=0;
  Int_t ievent=fHeader.GetEvent()+1;
  //
  // Energy loss information
  if(ievent) {
    printf("***************** Energy Loss Information per event (GEV) *****************\n");
    for(kn=1;kn<sEventEnergy.GetSize();kn++) {
      ed=sSummEnergy[kn];
      if(ed>0) {
	sEventEnergy[ndep]=kn;
	if(ievent>1) {
	  ed=ed/ievent;
	  ed2=sSum2Energy[kn];
	  ed2=ed2/ievent;
	  ed2=100*TMath::Sqrt(TMath::Max(ed2-ed*ed,zero))/ed;
	} else 
	  ed2=99;
	sSummEnergy[ndep]=ed;
	sSum2Energy[ndep]=TMath::Min((Float_t) 99.,TMath::Max(ed2,zero));
	edtot+=ed;
	ndep++;
      }
    }
    for(kn=0;kn<(ndep-1)/3+1;kn++) {
      left=ndep-kn*3;
      for(i=0;i<(3<left?3:left);i++) {
	j=kn*3+i;
        id=Int_t (sEventEnergy[j]+0.1);
	printf(" %s %10.3f +- %10.3f%%;",gMC->VolName(id),sSummEnergy[j],sSum2Energy[j]);
      }
      printf("\n");
    }
    //
    // Relative energy loss in different detectors
    printf("******************** Relative Energy Loss per event ********************\n");
    printf("Total energy loss per event %10.3f GeV\n",edtot);
    for(kn=0;kn<(ndep-1)/5+1;kn++) {
      left=ndep-kn*5;
      for(i=0;i<(5<left?5:left);i++) {
	j=kn*5+i;
        id=Int_t (sEventEnergy[j]+0.1);
	printf(" %s %10.3f%%;",gMC->VolName(id),100*sSummEnergy[j]/edtot);
      }
      printf("\n");
    }
    for(kn=0;kn<75;kn++) printf("*"); 
    printf("\n");
  }
  //
  // Reset the TArray's
  sEventEnergy.Set(0);
  sSummEnergy.Set(0);
  sSum2Energy.Set(0);
}

//_____________________________________________________________________________
AliModule *AliRun::GetModule(const char *name)
{
  //
  // Return pointer to detector from name
  //
  return (AliModule*)fModules->FindObject(name);
}
 
//_____________________________________________________________________________
AliDetector *AliRun::GetDetector(const char *name)
{
  //
  // Return pointer to detector from name
  //
  return (AliDetector*)fModules->FindObject(name);
}
 
//_____________________________________________________________________________
Int_t AliRun::GetModuleID(const char *name)
{
  //
  // Return galice internal detector identifier from name
  //
  Int_t i=-1;
  TObject *mod=fModules->FindObject(name);
  if(mod) i=fModules->IndexOf(mod);
  return i;
}
 
//_____________________________________________________________________________
Int_t AliRun::GetEvent(Int_t event)
{
  //
  // Connect the Trees Kinematics and Hits for event # event
  // Set branch addresses
  //

  // Reset existing structures
  ResetStack();
  ResetHits();
  ResetDigits();
  
  // Delete Trees already connected
  if (fTreeK) delete fTreeK;
  if (fTreeH) delete fTreeH;
  if (fTreeD) delete fTreeD;
  if (fTreeR) delete fTreeR;

  // Get header from file
  if(fTreeE) fTreeE->GetEntry(event);
  else Error("GetEvent","Cannot file Header Tree\n");
  
  // Get Kine Tree from file
  char treeName[20];
  sprintf(treeName,"TreeK%d",event);
  fTreeK = (TTree*)gDirectory->Get(treeName);
  if (fTreeK) fTreeK->SetBranchAddress("Particles", &fParticles);
  else    Error("GetEvent","cannot find Kine Tree for event:%d\n",event);
  
  // Get Hits Tree header from file
  sprintf(treeName,"TreeH%d",event);
  fTreeH = (TTree*)gDirectory->Get(treeName);
  if (!fTreeH) {
    Error("GetEvent","cannot find Hits Tree for event:%d\n",event);
  }
  
  // Get Digits Tree header from file
  sprintf(treeName,"TreeD%d",event);
  fTreeD = (TTree*)gDirectory->Get(treeName);
  if (!fTreeD) {
    Warning("GetEvent","cannot find Digits Tree for event:%d\n",event);
  }
  
  
  // Get Reconstruct Tree header from file
  sprintf(treeName,"TreeR%d",event);
  fTreeR = (TTree*)gDirectory->Get(treeName);
  if (!fTreeR) {
    //    printf("WARNING: cannot find Reconstructed Tree for event:%d\n",event);
  }
 
  // Set Trees branch addresses
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->SetTreeAddress();
  }
  
  if (fTreeK) fTreeK->GetEvent(0);
  fNtrack = Int_t (fParticles->GetEntries());
  return fNtrack;
}

//_____________________________________________________________________________
TGeometry *AliRun::GetGeometry()
{
  //
  // Import Alice geometry from current file
  // Return pointer to geometry object
  //
  if (!fGeometry) fGeometry = (TGeometry*)gDirectory->Get("AliceGeom");
  //
  // Unlink and relink nodes in detectors
  // This is bad and there must be a better way...
  //
  
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    detector->SetTreeAddress();
    TList *dnodes=detector->Nodes();
    Int_t j;
    TNode *node, *node1;
    for ( j=0; j<dnodes->GetSize(); j++) {
      node = (TNode*) dnodes->At(j);
      node1 = fGeometry->GetNode(node->GetName());
      dnodes->Remove(node);
      dnodes->AddAt(node1,j);
    }
  }
  return fGeometry;
}

//_____________________________________________________________________________
void AliRun::GetNextTrack(Int_t &mtrack, Int_t &ipart, Float_t *pmom,
			  Float_t &e, Float_t *vpos, Float_t *polar,
			  Float_t &tof)
{
  //
  // Return next track from stack of particles
  //
  TVector3 pol;
  fCurrent=-1;
  TParticle *track;
  for(Int_t i=fNtrack-1; i>=0; i--) {
    track=(TParticle*) fParticles->UncheckedAt(i);
    if(!track->TestBit(Done_Bit)) {
      //
      // The track has not yet been processed
      fCurrent=i;
      ipart=track->GetPdgCode();
      pmom[0]=track->Px();
      pmom[1]=track->Py(); 
      pmom[2]=track->Pz();
      e     =track->Energy();
      vpos[0]=track->Vx();
      vpos[1]=track->Vy();
      vpos[2]=track->Vz();
      track->GetPolarisation(pol);
      polar[0]=pol.X();
      polar[1]=pol.Y();
      polar[2]=pol.Z();
      tof=track->T();
      track->SetBit(Done_Bit);
      break;
    }
  }
  mtrack=fCurrent;
  //
  // stop and start timer when we start a primary track
  Int_t nprimaries = fHeader.GetNprimary();
  if (fCurrent >= nprimaries) return;
  if (fCurrent < nprimaries-1) {
    fTimer.Stop();
    track=(TParticle*) fParticles->UncheckedAt(fCurrent+1);
    //    track->SetProcessTime(fTimer.CpuTime());
  }
  fTimer.Start();
}

//_____________________________________________________________________________
Int_t AliRun::GetPrimary(Int_t track)
{
  //
  // return number of primary that has generated track
  //
  int current, parent;
  TParticle *part;
  //
  parent=track;
  while (1) {
    current=parent;
    part = (TParticle *)fParticles->UncheckedAt(current);
    parent=part->GetFirstMother();
    if(parent<0) return current;
  }
}
 
//_____________________________________________________________________________
void AliRun::Init(const char *setup)
{
  //
  // Initialize the Alice setup
  //

  gROOT->LoadMacro(setup);
  gInterpreter->ProcessLine("Config();");

  gMC->DefineParticles();  //Create standard MC particles

  TObject *objfirst, *objlast;

  fNdets = fModules->GetLast()+1;

  //
  //=================Create Materials, geometry, histograms, etc
   TIter next(fModules);
   AliModule *detector;
   while((detector = (AliModule*)next())) {
      detector->SetTreeAddress();
      objlast = gDirectory->GetList()->Last();
      
      // Initialise detector materials, geometry, histograms,etc
      detector->CreateMaterials();
      detector->CreateGeometry();
      detector->BuildGeometry();
      detector->Init();
      
      // Add Detector histograms in Detector list of histograms
      if (objlast) objfirst = gDirectory->GetList()->After(objlast);
      else         objfirst = gDirectory->GetList()->First();
      while (objfirst) {
	detector->Histograms()->Add(objfirst);
	objfirst = gDirectory->GetList()->After(objfirst);
      }
   }
   SetTransPar(); //Read the cuts for all materials
   
   MediaTable(); //Build the special IMEDIA table
   
   //Close the geometry structure
   gMC->Ggclos();
   
   //Initialise geometry deposition table
   sEventEnergy.Set(gMC->NofVolumes()+1);
   sSummEnergy.Set(gMC->NofVolumes()+1);
   sSum2Energy.Set(gMC->NofVolumes()+1);
   
   //Create the color table
   gMC->SetColors();
   
   //Compute cross-sections
   gMC->Gphysi();
   
   //Write Geometry object to current file.
   fGeometry->Write();
   
   fInitDone = kTRUE;
}

//_____________________________________________________________________________
void AliRun::MediaTable()
{
  //
  // Built media table to get from the media number to
  // the detector id
  //
  Int_t kz, nz, idt, lz, i, k, ind;
  //  Int_t ibeg;
  TObjArray &dets = *gAlice->Detectors();
  AliModule *det;
  //
  // For all detectors
  for (kz=0;kz<fNdets;kz++) {
    // If detector is defined
    if((det=(AliModule*) dets[kz])) {
        TArrayI &idtmed = *(det->GetIdtmed()); 
        for(nz=0;nz<100;nz++) {
	// Find max and min material number
	if((idt=idtmed[nz])) {
	  det->LoMedium() = det->LoMedium() < idt ? det->LoMedium() : idt;
	  det->HiMedium() = det->HiMedium() > idt ? det->HiMedium() : idt;
	}
      }
      if(det->LoMedium() > det->HiMedium()) {
	det->LoMedium() = 0;
	det->HiMedium() = 0;
      } else {
	if(det->HiMedium() > fImedia->GetSize()) {
	  Error("MediaTable","Increase fImedia from %d to %d",
		fImedia->GetSize(),det->HiMedium());
	  return;
	}
	// Tag all materials in rage as belonging to detector kz
	for(lz=det->LoMedium(); lz<= det->HiMedium(); lz++) {
	  (*fImedia)[lz]=kz;
	}
      }
    }
  }
  //
  // Print summary table
  printf(" Traking media ranges:\n");
  for(i=0;i<(fNdets-1)/6+1;i++) {
    for(k=0;k< (6<fNdets-i*6?6:fNdets-i*6);k++) {
      ind=i*6+k;
      det=(AliModule*)dets[ind];
      if(det)
	printf(" %6s: %3d -> %3d;",det->GetName(),det->LoMedium(),
	       det->HiMedium());
      else
	printf(" %6s: %3d -> %3d;","NULL",0,0);
    }
    printf("\n");
  }
}

//____________________________________________________________________________
void AliRun::SetGenerator(AliGenerator *generator)
{
  //
  // Load the event generator
  //
  if(!fGenerator) fGenerator = generator;
}

//____________________________________________________________________________
void AliRun::SetTransPar(char* filename)
{
  //
  // Read filename to set the transport parameters
  //


  const Int_t ncuts=10;
  const Int_t nflags=11;
  const Int_t npars=ncuts+nflags;
  const char pars[npars][7] = {"CUTGAM" ,"CUTELE","CUTNEU","CUTHAD","CUTMUO",
			       "BCUTE","BCUTM","DCUTE","DCUTM","PPCUTM","ANNI",
			       "BREM","COMP","DCAY","DRAY","HADR","LOSS",
			       "MULS","PAIR","PHOT","RAYL"};
  char line[256];
  char detName[7];
  char* filtmp;
  Float_t cut[ncuts];
  Int_t flag[nflags];
  Int_t i, itmed, iret, ktmed, kz;
  FILE *lun;
  //
  // See whether the file is there
  filtmp=gSystem->ExpandPathName(filename);
  lun=fopen(filtmp,"r");
  delete [] filtmp;
  if(!lun) {
    Warning("SetTransPar","File %s does not exist!\n",filename);
    return;
  }
  //
  printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
  printf(" *%59s\n","*");
  printf(" *       Please check carefully what you are doing!%10s\n","*");
  printf(" *%59s\n","*");
  //
  while(1) {
    // Initialise cuts and flags
    for(i=0;i<ncuts;i++) cut[i]=-99;
    for(i=0;i<nflags;i++) flag[i]=-99;
    itmed=0;
    for(i=0;i<256;i++) line[i]='\0';
    // Read up to the end of line excluded
    iret=fscanf(lun,"%[^\n]",line);
    if(iret<0) {
      //End of file
      fclose(lun);
      printf(" *%59s\n","*");
      printf(" "); for(i=0;i<60;i++) printf("*"); printf("\n");
      return;
    }
    // Read the end of line
    fscanf(lun,"%*c");
    if(!iret) continue;
    if(line[0]=='*') continue;
    // Read the numbers
    iret=sscanf(line,"%s %d %f %f %f %f %f %f %f %f %f %f %d %d %d %d %d %d %d %d %d %d %d",
		detName,&itmed,&cut[0],&cut[1],&cut[2],&cut[3],&cut[4],&cut[5],&cut[6],&cut[7],&cut[8],
		&cut[9],&flag[0],&flag[1],&flag[2],&flag[3],&flag[4],&flag[5],&flag[6],&flag[7],
		&flag[8],&flag[9],&flag[10]);
    if(!iret) continue;
    if(iret<0) {
      //reading error
      Warning("SetTransPar","Error reading file %s\n",filename);
      continue;
    }
    // Check that the module exist
    AliModule *mod = GetModule(detName);
    if(mod) {
      // Get the array of media numbers
      TArrayI &idtmed = *mod->GetIdtmed();
      // Check that the tracking medium code is valid
      if(0<=itmed && itmed < 100) {
	ktmed=idtmed[itmed];
	if(!ktmed) {
	  Warning("SetTransPar","Invalid tracking medium code %d for %s\n",itmed,mod->GetName());
	  continue;
	}
	// Set energy thresholds
	for(kz=0;kz<ncuts;kz++) {
	  if(cut[kz]>=0) {
	    printf(" *  %-6s set to %10.3E for tracking medium code %4d for %s\n",
		   pars[kz],cut[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,pars[kz],cut[kz]);
	  }
	}
	// Set transport mechanisms
	for(kz=0;kz<nflags;kz++) {
	  if(flag[kz]>=0) {
	    printf(" *  %-6s set to %10d for tracking medium code %4d for %s\n",
		   pars[ncuts+kz],flag[kz],itmed,mod->GetName());
	    gMC->Gstpar(ktmed,pars[ncuts+kz],Float_t(flag[kz]));
	  }
	}
      } else {
	Warning("SetTransPar","Invalid medium code %d *\n",itmed);
	continue;
      }
    } else {
      Warning("SetTransPar","Module %s not present\n",detName);
      continue;
    }
  }
}

//_____________________________________________________________________________
void AliRun::MakeTree(Option_t *option)
{
  //
  //  Create the ROOT trees
  //  Loop on all detectors to create the Root branch (if any)
  //

  //
  // Analyse options
  char *K = strstr(option,"K");
  char *H = strstr(option,"H");
  char *E = strstr(option,"E");
  char *D = strstr(option,"D");
  char *R = strstr(option,"R");
  //
  if (K && !fTreeK) fTreeK = new TTree("TreeK0","Kinematics");
  if (H && !fTreeH) fTreeH = new TTree("TreeH0","Hits");
  if (D && !fTreeD) fTreeD = new TTree("TreeD0","Digits");
  if (E && !fTreeE) fTreeE = new TTree("TE","Header");
  if (R && !fTreeR) fTreeR = new TTree("TreeR0","Reconstruction");
  if (fTreeH) fTreeH->SetAutoSave(1000000000); //no autosave
  //
  // Create a branch for hits/digits for each detector
  // Each branch is a TClonesArray. Each data member of the Hits classes
  // will be in turn a subbranch of the detector master branch
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     if (H || D || R) detector->MakeBranch(option);
  }
  //  Create a branch for particles
  if (fTreeK && K) fTreeK->Branch("Particles",&fParticles,4000);
  
  //  Create a branch for Header
  if (fTreeE && E) fTreeE->Branch("Header","AliHeader",&header,4000);
}

//_____________________________________________________________________________
Int_t AliRun::PurifyKine(Int_t lastSavedTrack, Int_t nofTracks)
{
  //
  // PurifyKine with external parameters
  //
  fHgwmk = lastSavedTrack;
  fNtrack = nofTracks;
  PurifyKine();
  return fHgwmk;
}

//_____________________________________________________________________________
void AliRun::PurifyKine()
{
  //
  // Compress kinematic tree keeping only flagged particles
  // and renaming the particle id's in all the hits
  //
  TClonesArray &particles = *fParticles;
  int nkeep=fHgwmk+1, parent, i;
  TParticle *part, *partnew, *father;
  AliHit *OneHit;
  int *map = new int[particles.GetEntries()];

  // Save in Header total number of tracks before compression
  fHeader.SetNtrack(fHeader.GetNtrack()+fNtrack-fHgwmk);

  // Preset map, to be removed later
  for(i=0; i<fNtrack; i++) {
    if(i<=fHgwmk) map[i]=i ; else map[i] = -99 ;}
  // Second pass, build map between old and new numbering
  for(i=fHgwmk+1; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    if(part->TestBit(Keep_Bit)) {
      
      // This particle has to be kept
      map[i]=nkeep;
      if(i!=nkeep) {
	
	// Old and new are different, have to copy
	partnew = (TParticle *)particles.UncheckedAt(nkeep);
	*partnew = *part;
      } else partnew = part;
      
      // as the parent is always *before*, it must be already
      // in place. This is what we are checking anyway!
      if((parent=partnew->GetFirstMother())>fHgwmk) {
	if(map[parent]==-99) printf("map[%d] = -99!\n",parent);
	partnew->SetFirstMother(map[parent]);
      }
      nkeep++;
    }
  }
  fNtrack=nkeep;
  
  // Fix daughters information
  for (i=0; i<fNtrack; i++) {
    part = (TParticle *)particles.UncheckedAt(i);
    parent = part->GetFirstMother();
    if(parent>=0) {
      father = (TParticle *)particles.UncheckedAt(parent);
      if(father->TestBit(Daughters_Bit)) {
      
	if(i<father->GetFirstDaughter()) father->SetFirstDaughter(i);
	if(i>father->GetLastDaughter())  father->SetLastDaughter(i);
      } else {
	// Iitialise daughters info for first pass
	father->SetFirstDaughter(i);
	father->SetLastDaughter(i);
	father->SetBit(Daughters_Bit);
      }
    }
  }
  
  // Now loop on all detectors and reset the hits
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
    if (!detector->Hits()) continue;
    TClonesArray &vHits=*(detector->Hits());
    if(vHits.GetEntries() != detector->GetNhits())
      printf("vHits.GetEntries()!=detector->GetNhits(): %d != %d\n",
	     vHits.GetEntries(),detector->GetNhits());
    for (i=0; i<detector->GetNhits(); i++) {
      OneHit = (AliHit *)vHits.UncheckedAt(i);
      OneHit->SetTrack(map[OneHit->GetTrack()]);
    }
  }

  fHgwmk=nkeep-1;
  particles.SetLast(fHgwmk);
  delete [] map;
}

//_____________________________________________________________________________
void AliRun::Reset(Int_t run, Int_t idevent)
{
  //
  //  Reset all Detectors & kinematics & trees
  //
  char hname[30];
  //
  ResetStack();
  ResetHits();
  ResetDigits();

  // Initialise event header
  fHeader.Reset(run,idevent);

  if(fTreeK) {
    fTreeK->Reset();
    sprintf(hname,"TreeK%d",idevent);
    fTreeK->SetName(hname);
  }
  if(fTreeH) {
    fTreeH->Reset();
    sprintf(hname,"TreeH%d",idevent);
    fTreeH->SetName(hname);
  }
  if(fTreeD) {
    fTreeD->Reset();
    sprintf(hname,"TreeD%d",idevent);
    fTreeD->SetName(hname);
  }
  if(fTreeR) {
    fTreeR->Reset();
    sprintf(hname,"TreeR%d",idevent);
    fTreeR->SetName(hname);
  }
}

//_____________________________________________________________________________
void AliRun::ResetDigits()
{
  //
  //  Reset all Detectors digits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetDigits();
  }
}

//_____________________________________________________________________________
void AliRun::ResetHits()
{
  //
  //  Reset all Detectors hits
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetHits();
  }
}

//_____________________________________________________________________________
void AliRun::ResetPoints()
{
  //
  // Reset all Detectors points
  //
  TIter next(fModules);
  AliModule *detector;
  while((detector = (AliModule*)next())) {
     detector->ResetPoints();
  }
}

//_____________________________________________________________________________
void AliRun::Run(Int_t nevent, const char *setup)
{
  //
  // Main function to be called to process a galice run
  // example
  //    Root > gAlice.Run(); 
  // a positive number of events will cause the finish routine
  // to be called
  //

  Int_t i, todo;
  // check if initialisation has been done
  if (!fInitDone) Init(setup);
  
  // Create the Root Tree with one branch per detector
  if(!fEvent) {
    gAlice->MakeTree("KHDER");
  }

  todo = TMath::Abs(nevent);
  for (i=0; i<todo; i++) {
  // Process one run (one run = one event)
     gAlice->Reset(fRun, fEvent);
     gMC->Gtrigi();
     gMC->Gtrigc();
     gMC->Gtrig();
     gAlice->FinishEvent();
     fEvent++;
  }
  
  // End of this run, close files
  if(nevent>0) gAlice->FinishRun();
}

//_____________________________________________________________________________
void AliRun::RunLego(const char *setup,Int_t ntheta,Float_t themin,
		     Float_t themax,Int_t nphi,Float_t phimin,Float_t phimax,
		     Float_t rmin,Float_t rmax,Float_t zmax)
{
  //
  // Generates lego plots of:
  //    - radiation length map phi vs theta
  //    - radiation length map phi vs eta
  //    - interaction length map
  //    - g/cm2 length map
  //
  //  ntheta    bins in theta, eta
  //  themin    minimum angle in theta (degrees)
  //  themax    maximum angle in theta (degrees)
  //  nphi      bins in phi
  //  phimin    minimum angle in phi (degrees)
  //  phimax    maximum angle in phi (degrees)
  //  rmin      minimum radius
  //  rmax      maximum radius
  //  
  //
  //  The number of events generated = ntheta*nphi
  //  run input parameters in macro setup (default="Config.C")
  //
  //  Use macro "lego.C" to visualize the 3 lego plots in spherical coordinates
  //Begin_Html
  /*
    <img src="picts/AliRunLego1.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego2.gif">
  */
  //End_Html
  //Begin_Html
  /*
    <img src="picts/AliRunLego3.gif">
  */
  //End_Html
  //

  // check if initialisation has been done
  if (!fInitDone) Init(setup);
  
  fLego = new AliLego("lego","lego");
  fLego->Init(ntheta,themin,themax,nphi,phimin,phimax,rmin,rmax,zmax);
  fLego->Run();
  
  // Create only the Root event Tree
  gAlice->MakeTree("E");
  
  // End of this run, close files
  gAlice->FinishRun();
}

//_____________________________________________________________________________
void AliRun::SetCurrentTrack(Int_t track)
{ 
  //
  // Set current track number
  //
  fCurrent = track; 
}
 
//_____________________________________________________________________________
void AliRun::SetTrack(Int_t done, Int_t parent, Int_t pdg, Float_t *pmom,
		      Float_t *vpos, Float_t *polar, Float_t tof,
		      const char *mecha, Int_t &ntr, Float_t weight)
{ 
  //
  // Load a track on the stack
  //
  // done     0 if the track has to be transported
  //          1 if not
  // parent   identifier of the parent track. -1 for a primary
  // pdg    particle code
  // pmom     momentum GeV/c
  // vpos     position 
  // polar    polarisation 
  // tof      time of flight in seconds
  // mecha    production mechanism
  // ntr      on output the number of the track stored
  //
  TClonesArray &particles = *fParticles;
  TParticle *particle;
  Float_t mass;
  const Int_t firstdaughter=-1;
  const Int_t lastdaughter=-1;
  const Int_t KS=0;
  //  const Float_t tlife=0;
  
  //
  // Here we get the static mass
  // For MC is ok, but a more sophisticated method could be necessary
  // if the calculated mass is required
  // also, this method is potentially dangerous if the mass
  // used in the MC is not the same of the PDG database
  //
  mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  Float_t e=TMath::Sqrt(mass*mass+pmom[0]*pmom[0]+
			pmom[1]*pmom[1]+pmom[2]*pmom[2]);
  
  //printf("Loading particle %s mass %f ene %f No %d ip %d pos %f %f %f mom %f %f %f KS %d m %s\n",
  //pname,mass,e,fNtrack,pdg,vpos[0],vpos[1],vpos[2],pmom[0],pmom[1],pmom[2],KS,mecha);
  
  particle=new(particles[fNtrack]) TParticle(pdg,KS,parent,-1,firstdaughter,
					     lastdaughter,pmom[0],pmom[1],pmom[2],
					     e,vpos[0],vpos[1],vpos[2],tof);
  //					     polar[0],polar[1],polar[2],tof,
  //					     mecha,weight);
  ((TParticle*)particles[fNtrack])->SetPolarisation(TVector3(polar[0],polar[1],polar[2]));
  ((TParticle*)particles[fNtrack])->SetWeight(weight);
  if(!done) particle->SetBit(Done_Bit);
  
  if(parent>=0) {
    particle=(TParticle*) fParticles->UncheckedAt(parent);
    particle->SetLastDaughter(fNtrack);
    if(particle->GetFirstDaughter()<0) particle->SetFirstDaughter(fNtrack);
  } else { 
    //
    // This is a primary track. Set high water mark for this event
    fHgwmk=fNtrack;
    //
    // Set also number if primary tracks
    fHeader.SetNprimary(fHgwmk+1);
    fHeader.SetNtrack(fHgwmk+1);
  }
  ntr = fNtrack++;
}

//_____________________________________________________________________________
void AliRun::KeepTrack(const Int_t track)
{ 
  //
  // flags a track to be kept
  //
  TClonesArray &particles = *fParticles;
  ((TParticle*)particles[track])->SetBit(Keep_Bit);
}
 
//_____________________________________________________________________________
void AliRun::StepManager(Int_t id) const
{
  //
  // Called at every step during transport
  //

  Int_t copy;
  //
  // --- If lego option, do it and leave 
  if (fLego) {
    fLego->StepManager();
    return;
  }
  //Update energy deposition tables
  sEventEnergy[gMC->CurrentVolID(copy)]+=gMC->Edep();
  
  //Call the appropriate stepping routine;
  AliModule *det = (AliModule*)fModules->At(id);
  if(det) det->StepManager();
}

//_____________________________________________________________________________
void AliRun::ReadEuclid(const char* filnam, const AliModule *det, char* topvol)
{
  //                                                                     
  //       read in the geometry of the detector in euclid file format    
  //                                                                     
  //        id_det : the detector identification (2=its,...)            
  //        topvol : return parameter describing the name of the top    
  //        volume of geometry.                                          
  //                                                                     
  //            author : m. maire                                        
  //                                                                     
  //     28.07.98
  //     several changes have been made by miroslav helbich
  //     subroutine is rewrited to follow the new established way of memory
  //     booking for tracking medias and rotation matrices.
  //     all used tracking media have to be defined first, for this you can use
  //     subroutine  greutmed.
  //     top volume is searched as only volume not positioned into another 
  //

  Int_t i, nvol, iret, itmed, irot, numed, npar, ndiv, iaxe;
  Int_t ndvmx, nr, flag;
  char key[5], card[77], natmed[21];
  char name[5], mother[5], shape[5], konly[5], volst[7000][5];
  char *filtmp;
  Float_t par[50];
  Float_t teta1, phi1, teta2, phi2, teta3, phi3, orig, step;
  Float_t xo, yo, zo;
  const Int_t maxrot=5000;
  Int_t idrot[maxrot],istop[7000];
  FILE *lun;
  //
  // *** The input filnam name will be with extension '.euc'
  filtmp=gSystem->ExpandPathName(filnam);
  lun=fopen(filtmp,"r");
  delete [] filtmp;
  if(!lun) {
    Error("ReadEuclid","Could not open file %s\n",filnam);
    return;
  }
  //* --- definition of rotation matrix 0 ---  
  TArrayI &idtmed = *(det->GetIdtmed());
  for(i=1; i<maxrot; ++i) idrot[i]=-99;
  idrot[0]=0;
  nvol=0;
 L10:
  for(i=0;i<77;i++) card[i]=0;
  iret=fscanf(lun,"%77[^\n]",card);
  if(iret<=0) goto L20;
  fscanf(lun,"%*c");
  //*
  strncpy(key,card,4);
  key[4]='\0';
  if (!strcmp(key,"TMED")) {
    sscanf(&card[5],"%d '%[^']'",&itmed,natmed);
    if( itmed<0 || itmed>=100 ) {
      Error("ReadEuclid","TMED illegal medium number %d for %s\n",itmed,natmed);
      exit(1);
    }
    //Pad the string with blanks
    i=-1;
    while(natmed[++i]);
    while(i<20) natmed[i++]=' ';
    natmed[i]='\0';
    //
    if( idtmed[itmed]<=0 ) {
      Error("ReadEuclid","TMED undefined medium number %d for %s\n",itmed,natmed);
      exit(1);
    }
    gMC->Gckmat(idtmed[itmed],natmed);
    //*
  } else if (!strcmp(key,"ROTM")) {
    sscanf(&card[4],"%d %f %f %f %f %f %f",&irot,&teta1,&phi1,&teta2,&phi2,&teta3,&phi3);
    if( irot<=0 || irot>=maxrot ) {
      Error("ReadEuclid","ROTM rotation matrix number %d illegal\n",irot);
      exit(1);
    }
    det->AliMatrix(idrot[irot],teta1,phi1,teta2,phi2,teta3,phi3);
    //*
  } else if (!strcmp(key,"VOLU")) {
    sscanf(&card[5],"'%[^']' '%[^']' %d %d", name, shape, &numed, &npar);
    if (npar>0) {
      for(i=0;i<npar;i++) fscanf(lun,"%f",&par[i]);
      fscanf(lun,"%*c");
    }
    gMC->Gsvolu( name, shape, idtmed[numed], par, npar);
    //*     save the defined volumes
    strcpy(volst[++nvol],name);
    istop[nvol]=1;
    //*
  } else if (!strcmp(key,"DIVN")) {
    sscanf(&card[5],"'%[^']' '%[^']' %d %d", name, mother, &ndiv, &iaxe);
    gMC->Gsdvn  ( name, mother, ndiv, iaxe );
    //*
  } else if (!strcmp(key,"DVN2")) {
    sscanf(&card[5],"'%[^']' '%[^']' %d %d %f %d",name, mother, &ndiv, &iaxe, &orig, &numed);
    gMC->Gsdvn2( name, mother, ndiv, iaxe, orig,idtmed[numed]);
    //*
  } else if (!strcmp(key,"DIVT")) {
    sscanf(&card[5],"'%[^']' '%[^']' %f %d %d %d", name, mother, &step, &iaxe, &numed, &ndvmx);
    gMC->Gsdvt ( name, mother, step, iaxe, idtmed[numed], ndvmx);
    //*
  } else if (!strcmp(key,"DVT2")) {
    sscanf(&card[5],"'%[^']' '%[^']' %f %d %f %d %d", name, mother, &step, &iaxe, &orig, &numed, &ndvmx);
    gMC->Gsdvt2 ( name, mother, step, iaxe, orig, idtmed[numed], ndvmx );
    //*
  } else if (!strcmp(key,"POSI")) {
    sscanf(&card[5],"'%[^']' %d '%[^']' %f %f %f %d '%[^']'", name, &nr, mother, &xo, &yo, &zo, &irot, konly);
    if( irot<0 || irot>=maxrot ) {
      Error("ReadEuclid","POSI %s#%d rotation matrix number %d illegal\n",name,nr,irot);
      exit(1);
    }
    if( idrot[irot] == -99) {
      Error("ReadEuclid","POSI %s#%d undefined matrix number %d\n",name,nr,irot);
      exit(1);
    }
    //*** volume name cannot be the top volume
    for(i=1;i<=nvol;i++) {
      if (!strcmp(volst[i],name)) istop[i]=0;
    }
    //*
    gMC->Gspos  ( name, nr, mother, xo, yo, zo, idrot[irot], konly );
    //*
  } else if (!strcmp(key,"POSP")) {
    sscanf(&card[5],"'%[^']' %d '%[^']' %f %f %f %d '%[^']' %d", name, &nr, mother, &xo, &yo, &zo, &irot, konly, &npar);
    if( irot<0 || irot>=maxrot ) {
      Error("ReadEuclid","POSP %s#%d rotation matrix number %d illegal\n",name,nr,irot);
      exit(1);
    }
    if( idrot[irot] == -99) {
      Error("ReadEuclid","POSP %s#%d undefined matrix number %d\n",name,nr,irot);
      exit(1);
    }
    if (npar > 0) {
      for(i=0;i<npar;i++) fscanf(lun,"%f",&par[i]);
      fscanf(lun,"%*c");
    }
    //*** volume name cannot be the top volume
    for(i=1;i<=nvol;i++) {
      if (!strcmp(volst[i],name)) istop[i]=0;
    }
    //*
    gMC->Gsposp ( name, nr, mother, xo,yo,zo, idrot[irot], konly, par, npar);
  }
  //*
  if (strcmp(key,"END")) goto L10;
  //* find top volume in the geometry
  flag=0;
  for(i=1;i<=nvol;i++) {
    if (istop[i] && flag) {
      Warning("ReadEuclid"," %s is another possible top volume\n",volst[i]);
    }
    if (istop[i] && !flag) {
      strcpy(topvol,volst[i]);
      printf(" *** GREUCL *** volume %s taken as a top volume\n",topvol);
      flag=1;
    }
  }
  if (!flag) {
    Warning("ReadEuclid","top volume not found\n");
  }
  fclose (lun);
  //*
  //*     commented out only for the not cernlib version
  printf(" *** GREUCL *** file: %s is now read in\n",filnam);
  //
  return;
  //*
  L20:
  Error("ReadEuclid","reading error or premature end of file\n");
}

//_____________________________________________________________________________
void AliRun::ReadEuclidMedia(const char* filnam, const AliModule *det)
{
  //                                                                     
  //       read in the materials and tracking media for the detector     
  //                   in euclid file format                             
  //                                                                     
  //       filnam: name of the input file                                
  //       id_det: id_det is the detector identification (2=its,...)     
  //                                                                     
  //            author : miroslav helbich                                
  //
  Float_t sxmgmx = gAlice->Field()->Max();
  Int_t   isxfld = gAlice->Field()->Integ();
  Int_t end, i, iret, itmed;
  char key[5], card[130], natmed[21], namate[21];
  Float_t ubuf[50];
  char* filtmp;
  FILE *lun;
  Int_t imate;
  Int_t nwbuf, isvol, ifield, nmat;
  Float_t a, z, dens, radl, absl, fieldm, tmaxfd, stemax, deemax, epsil, stmin;
  //
  end=strlen(filnam);
  for(i=0;i<end;i++) if(filnam[i]=='.') {
    end=i;
    break;
  }
  //
  // *** The input filnam name will be with extension '.euc'
  printf("The file name is %s\n",filnam); //Debug
  filtmp=gSystem->ExpandPathName(filnam);
  lun=fopen(filtmp,"r");
  delete [] filtmp;
  if(!lun) {
    Warning("ReadEuclidMedia","Could not open file %s\n",filnam);
    return;
  }
  //
  // Retrieve Mag Field parameters
  Int_t ISXFLD=gAlice->Field()->Integ();
  Float_t SXMGMX=gAlice->Field()->Max();
  //  TArrayI &idtmed = *(det->GetIdtmed());
  //
 L10:
  for(i=0;i<130;i++) card[i]=0;
  iret=fscanf(lun,"%4s %[^\n]",key,card);
  if(iret<=0) goto L20;
  fscanf(lun,"%*c");
  //*
  //* read material
  if (!strcmp(key,"MATE")) {
    sscanf(card,"%d '%[^']' %f %f %f %f %f %d",&imate,namate,&a,&z,&dens,&radl,&absl,&nwbuf);
    if (nwbuf>0) for(i=0;i<nwbuf;i++) fscanf(lun,"%f",&ubuf[i]);
    //Pad the string with blanks
    i=-1;
    while(namate[++i]);
    while(i<20) namate[i++]=' ';
    namate[i]='\0';
    //
    det->AliMaterial(imate,namate,a,z,dens,radl,absl,ubuf,nwbuf);
    //* read tracking medium
  } else if (!strcmp(key,"TMED")) {
    sscanf(card,"%d '%[^']' %d %d %d %f %f %f %f %f %f %d",
	   &itmed,natmed,&nmat,&isvol,&ifield,&fieldm,&tmaxfd,
	   &stemax,&deemax,&epsil,&stmin,&nwbuf);
    if (nwbuf>0) for(i=0;i<nwbuf;i++) fscanf(lun,"%f",&ubuf[i]);
    if (ifield<0) ifield=isxfld;
    if (fieldm<0) fieldm=sxmgmx;
    //Pad the string with blanks
    i=-1;
    while(natmed[++i]);
    while(i<20) natmed[i++]=' ';
    natmed[i]='\0';
    //
    det->AliMedium(itmed,natmed,nmat,isvol,ISXFLD,SXMGMX,tmaxfd,
		   stemax,deemax,epsil,stmin,ubuf,nwbuf);
    //    (*fImedia)[idtmed[itmed]-1]=id_det;
    //*
  }
  //*
  if (strcmp(key,"END")) goto L10;
  fclose (lun);
  //*
  //*     commented out only for the not cernlib version
  Warning("ReadEuclidMedia","file: %s is now read in\n",filnam);
  //*
  return;
  //*
 L20:
  Warning("ReadEuclidMedia","reading error or premature end of file\n");
} 
 
//_____________________________________________________________________________
void AliRun::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliRun.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    if (!gAlice) gAlice = this;
    gROOT->GetListOfBrowsables()->Add(this,"Run");
    fTreeE = (TTree*)gDirectory->Get("TE");
    if (fTreeE) fTreeE->SetBranchAddress("Header", &header);
    else    Error("Streamer","cannot find Header Tree\n");
    R__b >> fNtrack;
    R__b >> fHgwmk;
    R__b >> fDebug;
    fHeader.Streamer(R__b);
    R__b >> fModules;
    R__b >> fParticles;
    R__b >> fField; 
    //    R__b >> fMC;
    R__b >> fNdets;
    R__b >> fTrRmax;
    R__b >> fTrZmax;
    R__b >> fGenerator;
    if(R__v>1) {
      R__b >> fPDGDB;        //Particle factory object!
      fTreeE->GetEntry(0);
    } else {
      fHeader.SetEvent(0);
      fPDGDB     = TDatabasePDG::Instance();        //Particle factory object!
    }
  } else {
    R__b.WriteVersion(AliRun::IsA());
    TNamed::Streamer(R__b);
    R__b << fNtrack;
    R__b << fHgwmk;
    R__b << fDebug;
    fHeader.Streamer(R__b);
    R__b << fModules;
    R__b << fParticles;
    R__b << fField;
    //    R__b << fMC;
    R__b << fNdets;
    R__b << fTrRmax;
    R__b << fTrZmax;
    R__b << fGenerator;
    R__b << fPDGDB;        //Particle factory object!
  }
} 


//_____________________________________________________________________________
//
//                 Interfaces to Fortran
//
//_____________________________________________________________________________

extern "C" void type_of_call  rxgtrak (Int_t &mtrack, Int_t &ipart, Float_t *pmom, 
				       Float_t &e, Float_t *vpos, Float_t *polar,
				       Float_t &tof)
{
  //
  //     Fetches next track from the ROOT stack for transport. Called by the
  //     modified version of GTREVE.
  //
  //              Track number in the ROOT stack. If MTRACK=0 no
  //      mtrack  more tracks are left in the stack to be
  //              transported.
  //      ipart   Particle code in the GEANT conventions.
  //      pmom[3] Particle momentum in GeV/c
  //      e       Particle energy in GeV
  //      vpos[3] Particle position
  //      tof     Particle time of flight in seconds
  //
  Int_t pdg;
  gAlice->GetNextTrack(mtrack, pdg, pmom, e, vpos, polar, tof);
  ipart = gMC->IdFromPDG(pdg);
  mtrack++;
}

//_____________________________________________________________________________
extern "C" void type_of_call 
#ifndef WIN32
rxstrak (Int_t &keep, Int_t &parent, Int_t &ipart, Float_t *pmom, 
	       Float_t *vpos, Float_t &tof, const char* cmech, Int_t &ntr, const int cmlen)
#else
rxstrak (Int_t &keep, Int_t &parent, Int_t &ipart, Float_t *pmom,
	 Float_t *vpos, Float_t &tof, const char* cmech, const int cmlen,
	 Int_t &ntr)
#endif
{
  //
  //     Fetches next track from the ROOT stack for transport. Called by GUKINE
  //     and GUSTEP.
  //
  //              Status of the track. If keep=0 the track is put
  //      keep    on the ROOT stack but it is not fetched for
  //              transport.
  //      parent  Parent track. If parent=0 the track is a primary.
  //              In GUSTEP the routine is normally called to store
  //              secondaries generated by the current track whose
  //              ROOT stack number is MTRACK (common SCKINE.
  //      ipart   Particle code in the GEANT conventions.
  //      pmom[3] Particle momentum in GeV/c
  //      vpos[3] Particle position
  //      tof     Particle time of flight in seconds
  //
  //      cmech   (CHARACTER*10) Particle origin. This field is user
  //              defined and it is not used inside the GALICE code.
  //      ntr     Number assigned to the particle in the ROOT stack.
  //
  char mecha[11];
  Float_t polar[3]={0.,0.,0.};
  for(int i=0; i<10 && i<cmlen; i++) mecha[i]=cmech[i];
  mecha[10]=0;
  Int_t pdg=gMC->PDGFromId(ipart);
  gAlice->SetTrack(keep, parent-1, pdg, pmom, vpos, polar, tof, mecha, ntr);
  ntr++;
}

//_____________________________________________________________________________
extern "C" void type_of_call  rxkeep(const Int_t &n)
{
  if( NULL==gAlice ) exit(1);
  
  if( n<=0 || n>gAlice->Particles()->GetEntries() )
    {
      printf("  Bad index n=%d must be 0<n<=%d\n",
	     n,gAlice->Particles()->GetEntries());
      exit(1);
    }
  
  ((TParticle*)(gAlice->Particles()->UncheckedAt(n-1)))->SetBit(Keep_Bit);
}

//_____________________________________________________________________________
extern "C" void type_of_call  rxouth ()
{
  //
  // Called by Gtreve at the end of each primary track
  //
  gAlice->FinishPrimary();
}

