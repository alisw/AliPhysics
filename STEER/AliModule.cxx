///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Base class for ALICE modules. Both sensitive modules (Modules) and      //
// non-sensitive ones are described by this base class. This class           //
// supports the hit and digit trees produced by the simulation and also      //
// the objects produced by the reconstruction.                               //
//                                                                           //
// This class is also responsible for building the geometry of the           //
// Modules.                                                                //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliModuleClass.gif">
*/
//End_Html
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
#include "AliModule.h"
#include "AliRun.h"
#include "AliHit.h"
#include "AliPoints.h"
#include <TClass.h>
#include <TNode.h>
#include <TRandom.h>

ClassImp(AliModule)
 
//_____________________________________________________________________________
AliModule::AliModule()
{
  //
  // Default constructor for the AliModule class
  //
  fHistograms = 0;
  fNodes      = 0;
}
 
//_____________________________________________________________________________
AliModule::AliModule(const char* name,const char *title):TNamed(name,title)
{
  //
  // Normal constructor invoked by all Modules.
  // Create the list for Module specific histograms
  // Add this Module to the global list of Modules in Run.
  //
  //
  // Initialises the histogram list
  fHistograms = new TList();
  //
  // Initialises the list of ROOT TNodes
  fNodes      = new TList();
  //  
  // Get the Module numeric ID
  Int_t id = gAlice->GetModuleID(name);
  if (id < 0) {
    // Unknown Module !
     Warning("AliRun::Ctor","ERROR Unknown Module: %s\n",name);
     return;
  }
  //
  // Add this Module to the list of Modules
  gAlice->Modules()->AddAtAndExpand(this,id);
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
AliModule::~AliModule()
{
  //
  // Destructor
  //
  fHistograms = 0;
  //
  // Delete ROOT geometry
  fNodes->Clear();
  delete fNodes;
  //
  // Delete TArray objects
  delete fIdtmed;
  delete fIdmate;
}
 
//_____________________________________________________________________________
void AliModule::Disable()
{
  //
  // Disable Module on viewer
  //
  fActive = kFALSE;
  TIter next(fNodes);
  TNode *node;
  //
  // Loop through geometry to disable all
  // nodes for this Module
  while((node = (TNode*)next())) {
    node->SetVisibility(0);
  }   
}

//_____________________________________________________________________________
Int_t AliModule::DistancetoPrimitive(Int_t, Int_t)
{
  //
  // Return distance from mouse pointer to object
  // Dummy routine for the moment
  //
  return 9999;
}

//_____________________________________________________________________________
void AliModule::Enable()
{
  //
  // Enable Module on the viewver
  //
  fActive = kTRUE;
  TIter next(fNodes);
  TNode *node;
  //
  // Loop through geometry to enable all
  // nodes for this Module
  while((node = (TNode*)next())) {
    node->SetVisibility(1);
  }   
}

//_____________________________________________________________________________
void AliModule::AliMaterial(Int_t imat, const char* name, Float_t a, 
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
void AliModule::AliMixture(Int_t imat, const char *name, Float_t *a,
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
void AliModule::AliMedium(Int_t numed, const char *name, Int_t nmat,
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
void AliModule::AliMatrix(Int_t &nmat, Float_t theta1, Float_t phi1,
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
void AliModule::SetEuclidFile(char* material, char* geometry)
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
void AliModule::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class Module.
  //
  if (R__b.IsReading()) {
    Version_t R__v = R__b.ReadVersion(); if (R__v) { }
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    fEuclidMaterial.Streamer(R__b);
    fEuclidGeometry.Streamer(R__b);
    R__b >> fActive;
    R__b >> fHistograms;
    //
    // Stream the pointers but not the TClonesArrays
    R__b >> fNodes; // diff
  } else {
    R__b.WriteVersion(AliModule::IsA());
    TNamed::Streamer(R__b);
    TAttLine::Streamer(R__b);
    TAttMarker::Streamer(R__b);
    fEuclidMaterial.Streamer(R__b);
    fEuclidGeometry.Streamer(R__b);
    R__b << fActive;
    R__b << fHistograms;
    //
    // Stream the pointers but not the TClonesArrays
    R__b << fNodes; // diff
  }
}
 
