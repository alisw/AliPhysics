///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base class for ALICE modules. Both sensitive modules (detectors) and      //
// non-sensitive ones are described by this base class. This class           //
// supports the hit and digit trees produced by the simulation and also      //
// the objects produced by the reconstruction.                               //
//                                                                           //
// This class is also responsible for building the geometry of the           //
// detectors.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliDetectorClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliDetector.h"
#include "AliRun.h"
#include "AliHit.h"
#include "AliPoints.h"
#include "AliMC.h"
#include <TClass.h>
#include <TNode.h>
#include <TRandom.h>

// Static variables for the hit iterator routines
static Int_t sMaxIterHit=0;
static Int_t sCurIterHit=0;

ClassImp(AliDetector)
 
//_____________________________________________________________________________
AliDetector::AliDetector()
{
  //
  // Default constructor for the AliDetector class
  //
  fNhits      = 0;
  fNdigits    = 0;
  fHistograms = 0;
  fNodes      = 0;
  fPoints     = 0;
  fHits       = 0;
  fDigits     = 0;
  fTimeGate   = 200.e-9;
  fActive     = kTRUE;
  fBufferSize = 16000;
}
 
//_____________________________________________________________________________
AliDetector::AliDetector(const char* name,const char *title):TNamed(name,title)
{
  //
  // Normal constructor invoked by all Detectors.
  // Create the list for detector specific histograms
  // Add this Detector to the global list of Detectors in Run.
  //

  fTimeGate   = 200.e-9;
  fActive     = kTRUE;
  fNhits      = 0;
  fHits       = 0;
  fDigits     = 0;
  fNdigits    = 0;
  fPoints     = 0;
  fBufferSize = 16000;
  //
  // Initialises the histogram list
  fHistograms = new TList();
  //
  // Initialises the list of ROOT TNodes
  fNodes      = new TList();
  //  
  // Get the detector numeric ID
  Int_t id = gAlice->GetDetectorID(name);
  if (id < 0) {
    // Unknown detector !
     printf(" * AliRun::Ctor * ERROR Unknown detector: %s\n",name);
     return;
  }
  //
  // Add this detector to the list of detectors
  gAlice->Detectors()->AddAtAndExpand(this,id);
  //
  //
  SetMarkerColor(3);
  //
  // Allocate space for tracking media and material indexes
  fIdtmed = new TArrayI(100);
  fIdmate  = new TArrayI(100);
  for(Int_t i=0;i<100;i++) (*fIdmate)[i]=(*fIdtmed)[i]=0;
  //
  // Prepare to find the tracking media range
  fLoMedium = 65536;
  fHiMedium = 0;
}
 
//_____________________________________________________________________________
AliDetector::~AliDetector()
{
  //
  // Destructor
  //
  fNhits      = 0;
  fNdigits    = 0;
  fHistograms = 0;
  //
  // Delete ROOT geometry
  fNodes->Clear();
  delete fNodes;
  //
  // Delete space point structure
  if (fPoints) fPoints->Delete();
  delete fPoints;
  fPoints     = 0;
  //
  // Delete TArray objects
  delete fIdtmed;
  delete fIdmate;
}
 
//_____________________________________________________________________________
void AliDetector::Browse(TBrowser *b)
{
  //
  // Insert Detector objects in the list of objects to be browsed
  //
  char name[64];
  if( fHits == 0) return;
  TObject *obj;
  Int_t i, nobjects;
  //
  nobjects = fHits->GetEntries();
  for (i=0;i<nobjects;i++) {
    obj = fHits->At(i);
    sprintf(name,"%s_%d",obj->GetName(),i);
    b->Add(obj, &name[0]);
  }
}

//_____________________________________________________________________________
void AliDetector::Disable()
{
  //
  // Disable detector on viewer
  //
  fActive = kFALSE;
  TIter next(fNodes);
  TNode *node;
  //
  // Loop through geometry to disable all
  // nodes for this detector
  while((node = (TNode*)next())) {
    node->SetVisibility(0);
  }   
}

//_____________________________________________________________________________
Int_t AliDetector::DistancetoPrimitive(Int_t, Int_t)
{
  //
  // Return distance from mouse pointer to object
  // Dummy routine for the moment
  //
  return 9999;
}

//_____________________________________________________________________________
void AliDetector::Enable()
{
  //
  // Enable detector on the viewver
  //
  fActive = kTRUE;
  TIter next(fNodes);
  TNode *node;
  //
  // Loop through geometry to enable all
  // nodes for this detector
  while((node = (TNode*)next())) {
    node->SetVisibility(1);
  }   
}

//_____________________________________________________________________________
void AliDetector::FinishRun()
{
  //
  // Procedure called at the end of a run.
  //
}

//_____________________________________________________________________________
AliHit* AliDetector::FirstHit(Int_t track)
{
  //
  // Initialise the hit iterator
  // Return the address of the first hit for track
  // If track>=0 the track is read from disk
  // while if track<0 the first hit of the current
  // track is returned
  // 
  if(track>=0) {
    gAlice->ResetHits();
    gAlice->TreeH()->GetEvent(track);
  }
  //
  sMaxIterHit=fHits->GetEntriesFast();
  sCurIterHit=0;
  if(sMaxIterHit) return (AliHit*) fHits->UncheckedAt(0);
  else            return 0;
}

//_____________________________________________________________________________
AliHit* AliDetector::NextHit()
{
  //
  // Return the next hit for the current track
  //
  if(sMaxIterHit) {
    if(++sCurIterHit<sMaxIterHit) 
      return (AliHit*) fHits->UncheckedAt(sCurIterHit);
    else        
      return 0;
  } else {
    printf("* AliDetector::NextHit * Hit Iterator called without calling FistHit before\n");
    return 0;
  }
}

//_____________________________________________________________________________
void AliDetector::LoadPoints(Int_t)
{
  //
  // Store x, y, z of all hits in memory
  //
  if (fHits == 0) return;
  //
  if (fPoints == 0) fPoints = new TObjArray(gAlice->GetNtrack());
  Int_t nhits = fHits->GetEntriesFast();
  if (nhits == 0) return;
  AliHit *ahit;
  //
  AliPoints *points = 0;
  Int_t trko=-99, trk;
  //
  // Loop over all the hits and store their position
  for (Int_t hit=0;hit<nhits;hit++) {
    ahit = (AliHit*)fHits->UncheckedAt(hit);
    if(trko!=(trk=ahit->GetTrack())) {
      //
      // Initialise a new track
      trko=trk;
      points = new AliPoints(nhits);
      fPoints->AddAt(points,trk);
      points->SetMarkerColor(GetMarkerColor());
      points->SetMarkerStyle(GetMarkerStyle());
      points->SetMarkerSize(GetMarkerSize());
      points->SetDetector(this);
      points->SetParticle(trk);
    }
    points->SetPoint(hit,ahit->fX,ahit->fY,ahit->fZ);
  }
}

//_____________________________________________________________________________
void AliDetector::MakeBranch(Option_t *option)
{
  //
  // Create a new branch in the current Root Tree
  // The branch of fHits is automatically split
  //
  char branchname[10];
  sprintf(branchname,"%s",GetName());
  //
  // Get the pointer to the header
  char *H = strstr(option,"H");
  //
  if (fHits   && gAlice->TreeH() && H) {
    gAlice->TreeH()->Branch(branchname,&fHits, fBufferSize);
    printf("* AliDetector::MakeBranch * Making Branch %s for hits\n",branchname);
  }	
}

//_____________________________________________________________________________
void AliDetector::AliMaterial(Int_t imat, const char* name, Float_t a, 
			      Float_t z, Float_t dens, Float_t radl,
			      Float_t absl, Float_t *buf, Int_t nwbuf) const
{
  //
  // Store the parameters for a material
  //
  // imat        the material index will be stored in (*fIdmate)[imat]
  // name        material name
  // a           atomic mass
  // z           atomic number
  // dens        density
  // radl        radiation length
  // absl        absorbtion length
  // buf         adress of an array user words
  // nwbuf       number of user words
  //
  Int_t kmat;
  AliMC::GetMC()->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
  (*fIdmate)[imat]=kmat;
}
  

//_____________________________________________________________________________
void AliDetector::AliMixture(Int_t imat, const char *name, Float_t *a,
			     Float_t *z, Float_t dens, Int_t nlmat,
			     Float_t *wmat) const
{ 
  //
  // Defines mixture or compound imat as composed by 
  // nlmat materials defined by arrays a, z and wmat
  // 
  // If nlmat > 0 wmat contains the proportion by
  // weights of each basic material in the mixture  
  // 
  // If nlmat < 0 wmat contains the number of atoms 
  // of eack kind in the molecule of the compound
  // In this case, wmat is changed on output to the relative weigths.
  //
  // imat        the material index will be stored in (*fIdmate)[imat]
  // name        material name
  // a           array of atomic masses
  // z           array of atomic numbers
  // dens        density
  // nlmat       number of components
  // wmat        array of concentrations
  //
  Int_t kmat;
  AliMC::GetMC()->Mixture(kmat, name, a, z, dens, nlmat, wmat);
  (*fIdmate)[imat]=kmat;
} 
 
//_____________________________________________________________________________
void AliDetector::AliMedium(Int_t numed, const char *name, Int_t nmat,
			    Int_t isvol, Int_t ifield, Float_t fieldm,
			    Float_t tmaxfd, Float_t stemax, Float_t deemax,
			    Float_t epsil, Float_t stmin, Float_t *ubuf,
			    Int_t nbuf) const
{ 
  //
  // Store the parameters of a tracking medium
  //
  // numed       the medium number is stored into (*fIdtmed)[numed-1]
  // name        medium name
  // nmat        the material number is stored into (*fIdmate)[nmat]
  // isvol       sensitive volume if isvol!=0
  // ifield      magnetic field flag (see below)
  // fieldm      maximum magnetic field
  // tmaxfd      maximum deflection angle due to magnetic field
  // stemax      maximum step allowed
  // deemax      maximum fractional energy loss in one step
  // epsil       tracking precision in cm
  // stmin       minimum step due to continuous processes
  //
  // ifield =  0       no magnetic field
  //        = -1       user decision in guswim
  //        =  1       tracking performed with Runge Kutta
  //        =  2       tracking performed with helix
  //        =  3       constant magnetic field along z
  //  
  Int_t kmed;
  Int_t *idtmed = gAlice->Idtmed();
  AliMC::GetMC()->Medium(kmed,name, (*fIdmate)[nmat], isvol, ifield, fieldm,
			 tmaxfd, stemax, deemax, epsil,	stmin, ubuf, nbuf); 
  idtmed[numed-1]=kmed;
} 
 
//_____________________________________________________________________________
void AliDetector::AliMatrix(Int_t &nmat, Float_t theta1, Float_t phi1,
			    Float_t theta2, Float_t phi2, Float_t theta3,
			    Float_t phi3) const
{
  // 
  // Define a rotation matrix. Angles are in degrees.
  //
  // nmat        on output contains the number assigned to the rotation matrix
  // theta1      polar angle for axis I
  // phi1        azimuthal angle for axis I
  // theta2      polar angle for axis II
  // phi2        azimuthal angle for axis II
  // theta3      polar angle for axis III
  // phi3        azimuthal angle for axis III
  //
  AliMC::GetMC()->Matrix(nmat, theta1, phi1, theta2, phi2, theta3, phi3); 
} 

//_____________________________________________________________________________
void AliDetector::ResetDigits()
{
  //
  // Reset number of digits and the digits array
  //
  fNdigits   = 0;
  if (fDigits)   fDigits->Clear();
}

//_____________________________________________________________________________
void AliDetector::ResetHits()
{
  //
  // Reset number of hits and the hits array
  //
  fNhits   = 0;
  if (fHits)   fHits->Clear();
}

//_____________________________________________________________________________
void AliDetector::ResetPoints()
{
  //
  // Reset array of points
  //
  if (fPoints) {
    fPoints->Delete();
    delete fPoints;
    fPoints = 0;
  }
}

  
//_____________________________________________________________________________
void AliDetector::StepManager()
{
  //
  // Procedure called at every step inside the detector
  //
  printf("* AliDetector::StepManager * Generic Step Manager called for Detector: %s\n",fName.Data());
}

//_____________________________________________________________________________
void AliDetector::SetEuclidFile(char* material, char* geometry)
{
  //
  // Sets the name of the Euclid file
  //
  fEuclidMaterial=material;
  if(geometry) {
    fEuclidGeometry=geometry;
  } else {
    char* name = new char[strlen(material)];
    strcpy(name,material);
    strcpy(&name[strlen(name)-4],".euc");
    fEuclidGeometry=name;
    delete [] name;
  }
}
 
//_____________________________________________________________________________
void AliDetector::SetTreeAddress()
{
  //
  // Set branch address for the Hits and Digits Trees
  //
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  //
  // Branch address for hit tree
  TTree *treeH = gAlice->TreeH();
  if (treeH && fHits) {
    branch = treeH->GetBranch(branchname);
    if (branch) branch->SetAddress(&fHits);
  }
  //
  // Branch address for digit tree
  TTree *treeD = gAlice->TreeD();
  if (treeD && fDigits) {
    branch = treeD->GetBranch(branchname);
    if (branch) branch->SetAddress(&fDigits);
  }
}

//_____________________________________________________________________________
void AliDetector::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class Detector.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    fEuclidMaterial.Streamer(R__b);
    fEuclidGeometry.Streamer(R__b);
    R__b >> fTimeGate;
    R__b >> fActive;
    R__b >> fIshunt;
    //R__b >> fNhits;
    R__b >> fHistograms;
    //
    // Stream the pointers but not the TClonesArrays
    R__b >> fNodes; // diff

    R__b >> fHits; // diff
    R__b >> fDigits; // diff
    
  } else {
    R__b.WriteVersion(AliDetector::IsA());
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    fEuclidMaterial.Streamer(R__b);
    fEuclidGeometry.Streamer(R__b);
    R__b << fTimeGate;
    R__b << fActive;
    R__b << fIshunt;
    //R__b << fNhits;
    R__b << fHistograms;
    //
    // Stream the pointers but not the TClonesArrays
    R__b << fNodes; // diff
    
    R__b << fHits; // diff
    R__b << fDigits; // diff
  }
}
 
