// $Id$
// Category: geometry
//
// C++ interface to Geant3 basic routines 
// for building Geant4 geometry
//
// by V. Berejnoi, 25.2.1999
// materials, tracking media support 
// added by I.Hrivnacova, 27.5.1999

#include "TG4GeometryManager.h"
#include "TG4GeometryOutputManager.h"
#include "TG4PhysicsManager.h"
#include "TG4VSensitiveDetector.h"
#include "TG4Limits.h"
#include "TG4Globals.h"
#include "TG3Units.h"

#include <G3toG4.hh> 
#include <G3toG4BuildTree.hh>
#include <G3VolTable.hh>
#include <G3RotTable.hh>
#include <G3EleTable.hh>
#include <G3MatTable.hh>
#include <G3MedTable.hh>
#include <G3SensVolVector.hh>

#include <G4SDManager.hh>
#include <G4VSensitiveDetector.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Element.hh> 

#include <math.h>
#include <ctype.h>

// extern global method from g3tog4
void G3CLRead(G4String &, char *);

TG4GeometryManager* TG4GeometryManager::fgInstance = 0;

TG4GeometryManager::TG4GeometryManager() 
  : fMediumCounter(0),
    fMaterialCounter(0),
    fMatrixCounter(0),
    fUseG3TMLimits(false),
    fWriteGeometry(true)
{
//
  if (fgInstance) {
    TG4Globals::Exception(
      "TG4GeometryManager: attempt to create two instances of singleton.");
  }

  fOutputManager = new TG4GeometryOutputManager();

  fgInstance = this;
      
  // instantiate the default element table
  //TG4ElementTable::Instance();
}

TG4GeometryManager::TG4GeometryManager(const TG4GeometryManager& right) {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4GeometryManager singleton.");
}


TG4GeometryManager::~TG4GeometryManager() {
//
  delete fOutputManager;
}

//==================================================================== =========
//
// operators
//
//==================================================================== =========

TG4GeometryManager& 
TG4GeometryManager::operator=(const TG4GeometryManager& right)
{
  // check assignement to self
  if (this == &right) return *this;

  TG4Globals::Exception(
    "Attempt to assign TG4GeometryManager singleton.");
    
  return *this;  
}    
          

//=============================================================================
//
// private methods
//
//=============================================================================


G4double* TG4GeometryManager::CreateG4doubleArray(Float_t* array, 
               G4int size) const
{
// Converts Float_t* array to G4double*,
// !! The new array has to be deleted by user.
// ---

  G4double* doubleArray;
  if (size>0) {
    doubleArray = new G4double[size]; 
    for (G4int i=0; i<size; i++) doubleArray[i] = array[i];
  }
  else {
    doubleArray = 0; 
  }  
  return doubleArray;
}


G4String TG4GeometryManager::CutName(const char* name) const
{
// Removes spaces after the name if present.
// ---

  G4String cutName = name;
  G4int i = cutName.length();
  while (cutName(--i) == ' ') cutName = cutName(0,i);

  return cutName;
}  


void TG4GeometryManager::GstparCut(G4int itmed, TG3Cut par, G4double parval)
{
// Sets special tracking medium parameter. 
// It is applied to all logical volumes that use the specified 
// tracking medium.
// ---

  // get medium from table
  G3MedTableEntry* medium = G3Med.get(itmed);
  if (!medium) {
    G4String text = "TG4GeometryManager::GstparCut: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  

  // get/create user limits
  G4UserLimits* limits = medium->GetLimits();
  TG4Limits* tg4Limits;
  if (limits) {
    tg4Limits = dynamic_cast<TG4Limits*> (limits);
    if (!tg4Limits)
      G4Exception("TG4GeometryManager::GstparCut: Wrong limits type.");
  }    
  else {     
    tg4Limits = new TG4Limits();
    medium->SetLimits(tg4Limits);
    // add verbose 
    G4cout << "TG4GeometryManager::GstparCut: new TG4Limits() for medium " 
           << itmed << " has been created." << endl;  
  }	   
  // set parameter
  tg4Limits->SetG3Cut(par, parval*GeV);
}


void TG4GeometryManager::GstparFlag(G4int itmed, TG3Flag par, G4double parval)
{
// Sets special tracking medium parameter. 
// It is applied to all logical volumes that use the specified 
// tracking medium.
// ---

  // get medium from table
  G3MedTableEntry* medium = G3Med.get(itmed);
  if (!medium) {
    G4String text = "TG4GeometryManager::GstparFlag: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  

  // get/create user limits
  G4UserLimits* limits = medium->GetLimits();
  TG4Limits* tg4Limits;
  if (limits) {
    tg4Limits = dynamic_cast<TG4Limits*> (limits);
    if (!tg4Limits)
      G4Exception("TG4GeometryManager::GstparFlag: Wrong limits type.");
  }    
  else {     
    tg4Limits = new TG4Limits();
    medium->SetLimits(tg4Limits);
    // add verbose 
    G4cout << "TG4GeometryManager::GstparFlag: new TG4Limits() for medium " 
           << itmed << " has been created." << endl;  
  }
  // set parameter
  tg4Limits->SetG3Flag(par, parval);
}

 
void TG4GeometryManager::FillMediumIdVector()
{
// The second index for materials (having its origin in
// G4 tracking media concept) is stored in a vector of G4int
// parallel to G4MaterialTable.
// ---

  // initialize vector 
  G4int nofMaterials = G4Material::GetNumberOfMaterials();
  G4int i;
  for (i=0; i<nofMaterials; i++) 
    fMediumIdVector.push_back(0);
  
  // fill vector
  for (i=0; i<fMediumCounter; i++) {
    // material index (index in G4Material table)
    G3MedTableEntry* mte = G3Med.get(i+1);
    G4int materialIndex = mte->GetMaterial()->GetIndex();

    // medium index (ID of G3MedTableEntry)
    G4int mediumIndex = mte->GetID();
    
    // store medium index in the vector
    fMediumIdVector[materialIndex] = mediumIndex;
  }  

  // add verbose
  G4cout << "Total nof materials: " << nofMaterials << endl;
  G4cout << "Total nof tracking medias: " << fMediumCounter << endl;  
}    


//=============================================================================
//
// public methods - AliMC implementation
//
//=============================================================================

 
void TG4GeometryManager::Material(Int_t& kmat, const char* name, Float_t a, 
          Float_t z, Float_t dens, Float_t radl, Float_t absl, Float_t* buf, 
          Int_t nwbuf)
{
// Creates G4Material.
// !! Parameters radl, absl, buf, nwbuf are ignored in G4gsmate
// Comment: 
// absl - this parameter is ignored by GEANT3, too
// ---

    kmat = ++fMaterialCounter;
    G4double* bufin = CreateG4doubleArray(buf, nwbuf); 

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsmate(kmat, name, a, z, dens, radl, nwbuf, bufin); 

    G4gsmate(kmat, name, a, z, dens, radl, nwbuf, bufin); 

    delete [] bufin;

    if (nwbuf > 0) {  
      G4String matName = name;
      G4String text 
        = "TG4GeometryManager: user defined parameters for material ";
      text = text + matName;
      text = text + " are ignored by Geant4.";	
      TG4Globals::Warning(text);
    }
}
  
 
void TG4GeometryManager::Mixture(Int_t& kmat, const char *name, Float_t *a, 
          Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat)
{ 
// Creates G4Material composed of more elements.
// !! Parameters radl, absl, buf, nwbuf are ignored in G4gsmate
// Comment: 
// absl - this parameter is ignored by GEANT3, too
// ---

   Int_t npar = abs(nlmat);
   G4double *ain = CreateG4doubleArray(a, npar); 
   G4double *zin = CreateG4doubleArray(z, npar); 
   G4double *wmatin = CreateG4doubleArray(wmat, npar); 

   kmat = ++fMaterialCounter;

   // write token to the output file
   if (fWriteGeometry) 
     fOutputManager->WriteGsmixt(kmat, name, ain, zin, dens, nlmat, wmatin);

   G4gsmixt(kmat, name, ain, zin, dens, nlmat, wmatin);
   
   // !!! in Geant3:
   // After a call with ratios by number (negative number of elements), 
   // the ratio array is changed to the ratio by weight, so all successive 
   // calls with the same array must specify the number of elements as 
   // positive
   
   // wmatin may be modified
   for (G4int i=0; i<npar; i++) wmat[i] = wmatin[i]; 

   delete [] ain;
   delete [] zin;
   delete [] wmatin;
} 

void TG4GeometryManager::Medium(Int_t& kmed, const char *name, Int_t nmat, 
          Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
          Float_t stemax, Float_t deemax, Float_t epsil, 
          Float_t stmin, Float_t* ubuf, Int_t nbuf)
{ 
// Creates a temporary "medium" that is used for 
// assigning corresponding parameters to G4 objects:
// NTMED is stored as a second material index;
// ISVOL is used for builing TG3SensVolVector;
// STEMAX is passed in TG4Limits (if fUseG3TMLimits is set true);
// !! The other parameters (IFIELD, FIELDM, TMAXFD, DEEMAX, EPSIL, STMIN)
// are ignored by Geant4.
// ---

//  Geant3 desription:
//  ==================
//  NTMED  Tracking medium number
//  NAME   Tracking medium name
//  NMAT   Material number
//  ISVOL  Sensitive volume flag
//  IFIELD Magnetic field
//  FIELDM Max. field value (Kilogauss)
//  TMAXFD Max. angle due to field (deg/step)
//  STEMAX Max. step allowed
//  DEEMAX Max. fraction of energy lost in a step
//  EPSIL  Tracking precision (cm)
//  STMIN  Min. step due to continuos processes (cm)
//
//  IFIELD = 0 if no magnetic field; IFIELD = -1 if user decision in GUSWIM;
//  IFIELD = 1 if tracking performed with GRKUTA; IFIELD = 2 if tracking
//  performed with GHELIX; IFIELD = 3 if tracking performed with GHELX3.  
// ---

  kmed = ++fMediumCounter;

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGstmed(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, 
        stemax, deemax, epsil, stmin, 0, 0);

  G4gstmed(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
       epsil, stmin, 0, fUseG3TMLimits);
     // !! instead of the nbuf argument the bool fIsG3Default is passed

  if (nbuf > 0) {  
    G4String medName = name;
    G4String text
      = "TG4GeometryManager: user defined parameters for medium ";
    text = text + medName;
    text = text + " are ignored by Geant4.";  
    TG4Globals::Warning(text);
  }
} 


void TG4GeometryManager::Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
           Float_t thetaY, Float_t phiY, Float_t thetaZ, Float_t phiZ)
{
// Creates G4RotationMatrix.
// ---

  krot = ++fMatrixCounter;

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGsrotm(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ);

  G4gsrotm(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ);
}
  

G4Material* TG4GeometryManager::MixMaterials(G4String name, G4double density, 
        TG4StringVector* matNames, TG4doubleVector* matWeights)
{
// Creates a mixture of selected materials
// ---

  // number of materials to be mixed  
  G4int nofMaterials = matNames->entries();
  if (nofMaterials != matWeights->entries()) {
    G4String text = "TG4GeometryManager::MixMaterials: ";
    text = text +  "different number of material names and weigths.";
    TG4Globals::Exception(text);
  }    
  // add verbose
  // G4cout << "Nof of materials to be mixed: " << nofMaterials << endl;

  // fill vector of materials
  TG4MaterialVector matVector;  
  G4int im;
  for (im=0; im< nofMaterials; im++) {
    // material
    G4Material* material = G4Material::GetMaterial((*matNames)[im]);
    matVector.insert(material);
  } 

  // create the mixed material
  G4Material* mixture = new G4Material(name, density, nofMaterials);
  for (im=0; im< nofMaterials; im++) {
    G4Material* material = matVector[im];
    G4double fraction = (*matWeights)[im];
    mixture->AddMaterial(material, fraction);
  }

  return mixture;
}  

 
void TG4GeometryManager::Ggclos() 
{ 
// Sets the top VTE in temporary G3 volume table.
// Close geometry output file (if fWriteGeometry is set true).
//
//  Geant3 desription:
//  ==================
// close out the geometry
// ---

  if (fWriteGeometry) fOutputManager->WriteGgclos();

  G4ggclos();    
} 
 

void TG4GeometryManager::Gfmate(Int_t imat, char *name, Float_t &a, 
          Float_t &z, Float_t &dens, Float_t &radl, Float_t &absl,
          Float_t* ubuf, Int_t& nbuf) 
{ 
//  Geant3 desription:
//  ==================
// Return parameters for material IMAT 
// ---

  G4Material* material = G3Mat.get(imat);
  
  if (material) {
    // to do: change this correctly 
    // !! unsafe conversion
    const char* chName = material->GetName();
    name = (char*)chName;
    a = GetEffA(material);
    z = GetEffZ(material);
    
    dens = material->GetDensity();
    dens /= TG3Units::MassDensity();

    radl = material->GetRadlen();
    radl /= TG3Units::Length();

    // the following parameters are not defined in Geant4
    absl = 0.; 
    ubuf = 0;
    nbuf = 0;
  }
  else {
    TG4Globals::Exception(
     "TG4GeometryManager::Gfmate: material has not been found.");
  }
} 

 
void  TG4GeometryManager::Gstpar(Int_t itmed, const char *param, 
           Float_t parval) 
{ 
// Passes the tracking medium parameter to TG4Limits.
// The tracking medium parameter is set only in case
// its value is different from the "global" physics setup.
// (If: CheckCut/FlagWithG3Defaults is used checking
//  is performed with respect to G3 default values.)
// When any cut/flag parameter is set in limits
// the physics manager is locked and the physics setup
// cannot be changed.
// Applying this TG4Limits to particles and physical
// processes is still in development.
//
//  Geant3 desription:
//  ==================
//  To change the value of cut  or mechanism "CHPAR"
//  to a new value PARVAL  for tracking medium ITMED
//  The  data   structure  JTMED   contains  the   standard  tracking
//  parameters (CUTS and flags to control the physics processes)  which
//  are used  by default  for all  tracking media.   It is  possible to
//  redefine individually  with GSTPAR  any of  these parameters  for a
//  given tracking medium. 
//  ITMED     tracking medium number 
//  CHPAR     is a character string (variable name) 
//  PARVAL    must be given as a floating point.
// ---

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGstpar(itmed, param, parval); 

  // get physics setup
  TG4PhysicsManager* physicsManager = TG4PhysicsManager::Instance();
  //TG4CutVector* cutVector = physicsManager->GetCutVector();
  //TG4FlagVector* flagVector = physicsManager->GetFlagVector();

  G4String name = CutName(param); 
  TG3Cut g3Cut;
  if (physicsManager->CheckCutWithCutVector(name, parval, g3Cut)) {
      GstparCut(itmed, g3Cut, parval);
      physicsManager->Lock();
  }  
  else {
    TG3Flag g3Flag;
    if (physicsManager->CheckFlagWithFlagVector(name, parval, g3Flag)) {
      GstparFlag(itmed, g3Flag, parval);
      physicsManager->Lock();
    } 
    else if (g3Cut==kNoG3Cuts && g3Flag==kNoG3Flags) { 
      G4String text = "TG4GeometryManager::Gstpar:";
      text = text + name;
      text = text + " parameter is not yet implemented.";
      TG4Globals::Warning(text);
    }	
  } 
} 
 
 
void  TG4GeometryManager::Gsckov(Int_t itmed, Int_t npckov, Float_t* ppckov,
			 Float_t* absco, Float_t* effic, Float_t* rindex)
{
//
//  Geant3 desription:
//  ==================
//
//    Stores the tables for UV photon tracking in medium ITMED 
//    Please note that it is the user's responsability to 
//    provide all the coefficients:
//
//
//       ITMED       Tracking medium number
//       NPCKOV      Number of bins of each table
//       PPCKOV      Value of photon momentum (in GeV)
//       ABSCO       Absorbtion coefficients 
//                   dielectric: absorbtion length in cm
//                   metals    : absorbtion fraction (0<=x<=1)
//       EFFIC       Detection efficiency for UV photons 
//       RINDEX      Refraction index (if=0 metal)
// ---

  G4double* ppckovDbl = CreateG4doubleArray(ppckov, npckov); 
  G4double* abscoDbl  = CreateG4doubleArray(absco, npckov); 
  G4double* efficDbl  = CreateG4doubleArray(effic, npckov); 
  G4double* rindexDbl = CreateG4doubleArray(rindex, npckov); 
  
  // add units
  G4int i;
  for (i=0; i<npckov; i++) {
    ppckovDbl[i] = ppckovDbl[i]*TG3Units::Energy();
    abscoDbl[i]  = abscoDbl[i]*TG3Units::Length();
  }  

  // create material properties table
  G4MaterialPropertiesTable* table = new G4MaterialPropertiesTable(); 
  table->AddProperty("ABSLENGTH", ppckovDbl, abscoDbl, npckov);
                    // used in G4OpAbsorption process
  table->AddProperty("EFFICIENCY", ppckovDbl, efficDbl, npckov);
                    // used in G4OpBoundary process
  table->AddProperty("RINDEX", ppckovDbl, rindexDbl, npckov);
                    // used in G4Cerenkov, G4OpRayleigh, G4OpBoundary

  // get material of medium from table
  G3MedTableEntry* medium = G3Med.get(itmed);
  if (!medium) {
    G4String text = "TG4GeometryManager::Gsckov: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  
  G4Material* material = medium->GetMaterial();
  
  // set material properties table 
  material->SetMaterialPropertiesTable(table);

  G4cout << "The tables for UV photon tracking set for "
         << material->GetName() << endl;
  for (i=0; i<npckov; i++)
    G4cout << ppckovDbl[i] << " " << rindexDbl[i] << endl;
	 
  delete ppckovDbl;
  delete abscoDbl;
  delete efficDbl;
  delete rindexDbl;	 
}			 

 
void  TG4GeometryManager::Gsdvn(const char *name, const char *mother, 
           Int_t ndiv, Int_t iaxis) 
{ 
//  Geant3 desription:
//  ==================
//  NAME   Volume name
//  MOTHER Mother volume name
//  NDIV   Number of divisions
//  IAXIS  Axis value
//
//  X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
//  It divides a previously defined volume.
//  ---

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsdvn(name, mother, ndiv, iaxis);

    G4gsdvn(CutName(name), CutName(mother), ndiv, iaxis);

    // register name in name map
    fNameMap.AddName(CutName(name));
} 
 
 
void  TG4GeometryManager::Gsdvn2(const char *name, const char *mother, 
           Int_t ndiv, Int_t iaxis, Float_t c0i, Int_t numed) 
{ 
//  Geant3 desription:
//  ==================
//  DIVIDES MOTHER INTO NDIV DIVISIONS CALLED NAME
//  ALONG AXIS IAXIS STARTING AT COORDINATE VALUE C0.
//  THE NEW VOLUME CREATED WILL BE MEDIUM NUMBER NUMED.
// ---

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsdvn2(name, mother, ndiv, iaxis, c0i, numed);

    G4gsdvn2(CutName(name),CutName(mother), ndiv, iaxis, c0i, numed);

    // register name in name map
    fNameMap.AddName(CutName(name));
} 
 
 
void  TG4GeometryManager::Gsdvt(const char *name, const char *mother, 
           Float_t step, Int_t iaxis, Int_t numed, Int_t ndvmx)
{ 
//  Geant3 desription:
//  ==================
//       Divides MOTHER into divisions called NAME along
//       axis IAXIS in steps of STEP. If not exactly divisible 
//       will make as many as possible and will centre them 
//       with respect to the mother. Divisions will have medium 
//       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
//       NDVMX is the expected maximum number of divisions
//          (If 0, no protection tests are performed) 
// ---

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsdvt(name, mother, step, iaxis, numed, ndvmx);

    G4gsdvt(CutName(name), CutName(mother), step, iaxis, numed, ndvmx);

    // register name in name map
    fNameMap.AddName(CutName(name));
} 
 
 
void  TG4GeometryManager::Gsdvt2(const char *name, const char *mother, 
           Float_t step, Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx)
{ 
//  Geant3 desription:
//  ==================
// Create a new volume by dividing an existing one
//                                                                    
//           Divides MOTHER into divisions called NAME along          
//            axis IAXIS starting at coordinate value C0 with step    
//            size STEP.                                              
//           The new volume created will have medium number NUMED.    
//           If NUMED is 0, NUMED of mother is taken.                 
//           NDVMX is the expected maximum number of divisions        
//             (If 0, no protection tests are performed)              
// ---

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsdvt2(name, mother, step, iaxis, c0, numed, ndvmx);

    G4gsdvt2(CutName(name), CutName(mother), step, iaxis, c0, numed, ndvmx);

    // register name in name map
    fNameMap.AddName(CutName(name));
} 
 
 
void  TG4GeometryManager::Gsord(const char *name, Int_t iax) 
{ 
// No corresponding action in G4.
//
//  Geant3 desription:
//  ==================
//    Flags volume CHNAME whose contents will have to be ordered 
//    along axis IAX, by setting the search flag to -IAX
//           IAX = 1    X axis 
//           IAX = 2    Y axis 
//           IAX = 3    Z axis 
//           IAX = 4    Rxy (static ordering only  -> GTMEDI)
//           IAX = 14   Rxy (also dynamic ordering -> GTNEXT)
//           IAX = 5    Rxyz (static ordering only -> GTMEDI)
//           IAX = 15   Rxyz (also dynamic ordering -> GTNEXT)
//           IAX = 6    PHI   (PHI=0 => X axis)
//           IAX = 7    THETA (THETA=0 => Z axis)
// ---

  TG4Globals::Warning("TG4GeometryManager::Gsord: dummy method.");
} 
 
 
void  TG4GeometryManager::Gspos(const char *vname, Int_t num, 
          const char *vmoth, Float_t x, Float_t y, Float_t z, Int_t irot, 
          const char *vonly) 
{ 
//  Geant3 desription:
//  ==================
// Position a volume into an existing one
//
//  NAME   Volume name
//  NUMBER Copy number of the volume
//  MOTHER Mother volume name
//  X      X coord. of the volume in mother ref. sys.
//  Y      Y coord. of the volume in mother ref. sys.
//  Z      Z coord. of the volume in mother ref. sys.
//  IROT   Rotation matrix number w.r.t. mother ref. sys.
//  ONLY   ONLY/MANY flag
//
//  It positions a previously defined volume in the mother.
// ---  

   // write token to the output file
   if (fWriteGeometry) 
     fOutputManager->WriteGspos(vname, num, vmoth, x, y, z, irot, vonly);

   G4gspos(CutName(vname), num, CutName(vmoth), x, y, z, irot, vonly);

   // register name in name map
   fNameMap.AddName(CutName(vname));
} 
 
 
void  TG4GeometryManager::Gsposp(const char *name, Int_t nr, 
           const char *mother, Float_t x, Float_t y, Float_t z, Int_t irot, 
           const char *konly, Float_t *upar, Int_t np ) 
{ 
//  Geant3 desription:
//  ==================
//      Place a copy of generic volume NAME with user number
//      NR inside MOTHER, with its parameters UPAR(1..NP)
// ---

   G4double* parin = CreateG4doubleArray(upar, np); 

   // write token to the output file
   if (fWriteGeometry) 
     fOutputManager->WriteGsposp(name, nr, mother, x, y, z, irot, konly, parin, np);

   G4gsposp(CutName(name), nr, CutName(mother), x, y, z, irot, konly, 
             parin, np);

   delete [] parin;

   // register name in name map
   fNameMap.AddName(CutName(name));
} 
 
 
Int_t TG4GeometryManager::Gsvolu(const char *name, const char *shape, 
          Int_t nmed, Float_t *upar, Int_t npar) 
{ 
//  Geant3 desription:
//  ==================
//  NAME   Volume name
//  SHAPE  Volume type
//  NUMED  Tracking medium number
//  NPAR   Number of shape parameters
//  UPAR   Vector containing shape parameters
//
//  It creates a new volume in the JVOLUM data structure.
// ---  

  G4double* parin = CreateG4doubleArray(upar, npar); 

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGsvolu(name, shape, nmed, parin, npar);    

  G4gsvolu(CutName(name), CutName(shape), nmed, parin, npar);

  delete [] parin;
  
  // register name in name map
  fNameMap.AddName(CutName(name));

  return 0;
} 


void TG4GeometryManager::WriteEuclid(const char* fileName, 
          const char* topVolName, Int_t number, Int_t nlevel)
{
//  Geant3 desription:
//  ==================
//
//     ******************************************************************
//     *                                                                *
//     *  Write out the geometry of the detector in EUCLID file format  *
//     *                                                                *
//     *       filnam : will be with the extension .euc                 *
//     *       topvol : volume name of the starting node                *
//     *       number : copy number of topvol (relevant for gsposp)     *
//     *       nlevel : number of  levels in the tree structure         *
//     *                to be written out, starting from topvol         *
//     *                                                                *
//     *       Author : M. Maire                                        *
//     *                                                                *
//     ******************************************************************
//
//     File filnam.tme is written out with the definitions of tracking
//     medias and materials.
//     As to restore original numbers for materials and medias, program
//     searches in the file euc_medi.dat and comparing main parameters of
//     the mat. defined inside geant and the one in file recognizes them
//     and is able to take number from file. If for any material or medium,
//     this procedure fails, ordering starts from 1.
//     Arrays IOTMED and IOMATE are used for this procedure
// ---

  TG4Globals::Warning(
    "TG4GeometryManager::WriteEuclid(..) is not yet implemented.");
}

 
Int_t TG4GeometryManager::VolId(const Text_t* volName) const
{ 
// Returns the sensitive detector identifier.
// !! Gives exception in case logical volume is not associated with 
// a sensitive detector.
// ---

  G4String g4VolName = CutName(volName);
  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
  
    if ((sd) && (sd->GetName()==g4VolName)) {
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd)
        return tsd->GetID();
      else {
        TG4Globals::Exception(
          "TG4GeometryManager::VolId: Unknown sensitive detector type");
        return 0;
      }   	
    }   
  }

  G4String text = "TG4GeometryManager::VolId: Sensitive detector ";
  text = text + g4VolName;
  text = text + " is not defined.\n"; 
  text = text + "    Set /alDet/setAllSensitive true in PreInit.";
  TG4Globals::Exception(text);
  return 0;
}


const char* TG4GeometryManager::VolName(Int_t id) const
{
// Returns the name of the sensitive detector with the given identifier.
// !! Gives exception in case logical volume is not associated with 
// a sensitive detector.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = pLVStore->at(i);
    G4VSensitiveDetector* sd = lv->GetSensitiveDetector();
    
    if (sd) {
      G4int sdID;
      TG4VSensitiveDetector* tsd = dynamic_cast<TG4VSensitiveDetector*>(sd);
      if (tsd)
        sdID = tsd->GetID();
      else {
        TG4Globals::Exception(
          "TG4GeometryManager::VolId: Unknown sensitive detector type");
        return 0;
      }   
      if (sdID == id) return sd->GetName();
    }  
  }

  G4String text = "TG4GeometryManager::VolName:\n";
  text = text + "    Sensitive detector with given id is not defined. \n";
  text = text + "    Set /alDet/setAllSensitive true in PreInit.";
  TG4Globals::Exception(text);
  return "";	       	         
}


Int_t TG4GeometryManager::NofVolumes() const
{
// Returns the total number of sensitive detectors.
// ---

  return NofSensitiveDetectors();
}

 
void TG4GeometryManager::ReadG3Geometry(G4String filePath)
{
// Processes g3calls.dat file and fills G3 tables.
// ---

  // add verbose
  G4cout << "Reading the call list file " << filePath << "..." << endl;
  G3CLRead(filePath, NULL);
  G4cout << "Call list file read completed. Build geometry" << endl;
}

 
//=============================================================================
//
// public methods - Geant4 only
//
//=============================================================================

 
G4VPhysicalVolume* TG4GeometryManager::CreateG4Geometry()
{
// Creates G4 geometry objects according to the G3VolTable 
// and returns the top physical volume in case it was created
// (return zero otherwise).
// ---

  // set the first entry in the G3Vol table
  Ggclos();
  G3VolTableEntry* first = G3Vol.GetFirstVTE();
  
  // close g3calls.dat
  if (fWriteGeometry) fOutputManager->CloseFile();  

  // create G4 geometry
  G3toG4BuildTree(first,0);  

  // print G3 volume table statistics
  G3Vol.VTEStat();

  // print G4 geometry statistics
  G4cout << "G4 Stat: instantiated " 
         << NofG4LogicalVolumes()  << " logical volumes \n"
	 << "                      " 
	 << NofG4PhysicalVolumes() << " physical volumes" << endl;

  // position the first entry 
  // (in Geant3 the top volume cannot be positioned)
  // 
  G4VPhysicalVolume* top = 0;
  if (first->GetLV()->GetNoDaughters() == 0) {
    top = new G4PVPlacement(0, G4ThreeVector(), first->GetName(), 
                            first->GetLV(), 0, false, 0);
  }
  return top;		      
}

 
void TG4GeometryManager::UseG3TrackingMediaLimits()
{
// Sets fUseG3TMLimits option.
// !! This method has to be called only before starting
// creating geometry.
// ---

  if (fMediumCounter == 0) {
    fUseG3TMLimits = true;
  }
  else {
    G4String text = "TG4GeometryManager::UseG3TMLimits: \n";
    text = text + "    It is too late to set G3 defaults. \n";
    text = text + "    Some media has been already processed.";
    TG4Globals::Exception(text);
  }
}

 
void TG4GeometryManager::ClearG3Tables()
{ 
// Clears G3 volumes, materials, rotations(?) tables
// and sensitive volumes vector.
// The top volume is kept in the vol table.
// ---

  // clear volume table 
  // but keep the top volume in the table 
  G3VolTableEntry* top = G3Vol.GetFirstVTE();
  G4String name = top->GetName();
  G4String shape = top->GetShape(); 
  G3VolTableEntry* keep 
    = new G3VolTableEntry(name, shape, top->GetRpar(), top->GetNpar(), 
                           top->GetNmed(), top->GetSolid(), false);
  keep->SetLV(top->GetLV());
  G3Vol.Clear();  
  G3Vol.PutVTE(keep);
  
  // clear other tables
  G3Mat.Clear();
  //G3Rot.Clear();
  G3SensVol.clear(); 
}

 
void TG4GeometryManager::ClearG3TablesFinal()
{
// Clears G3 medias and volumes tables
// (the top volume is removed from the vol table)
// ---

  // fill medium id vector
  FillMediumIdVector();

  G3Med.Clear();
  G3Vol.Clear();  

  // reset medium counter
  //fMaterialCounter = 0;
  fMediumCounter = 0;  
}

 
void TG4GeometryManager::OpenOutFile(G4String filePath)
{ 
// Opens output file.
// ---

  fOutputManager->OpenFile(filePath);
}

 
void TG4GeometryManager::PrintNameMap()
{
// Prints the map of volumes names to second names.
// ---

  fNameMap.PrintAll();
}

 
void TG4GeometryManager::SetWriteGeometry(G4bool writeGeometry)
{ 
// Controls geometry output.
// ---

  fWriteGeometry = writeGeometry; 
}

 
void TG4GeometryManager::SetMapSecond(const G4String& name)
{
// Sets the second name for the map of volumes names.
// ---

  fNameMap.SetSecond(name);
}


Int_t TG4GeometryManager::NofG3Volumes() const
{
// Returns the total number of logical volumes corresponding
// to G3 volumes. (
// The logical volume that were created by Gsposp method 
// with a generic name (name_copyNo) are NOT included.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();

  G4int counter = 0;  
  for (G4int i=0; i<pLVStore->entries(); i++) {
    G4LogicalVolume* lv = (*pLVStore)[i];
    if (IsG3Volume(lv->GetName())) counter++;
  }
  
  return counter;  
}


Int_t TG4GeometryManager::NofG4LogicalVolumes() const
{
// Returns the total number of logical volumes in the geometry.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();
  return pLVStore->entries();
}


Int_t TG4GeometryManager::NofG4PhysicalVolumes() const
{
// Returns the total number of physical volumes in the geometry.
// ---

  G4LogicalVolumeStore* pLVStore = G4LogicalVolumeStore::GetInstance();

  G4int counter = 0;
  for (G4int i=0; i<pLVStore->entries(); i++) {
    counter += ((*pLVStore)[i])->GetNoDaughters();
  }
  
  return counter;  
}


Int_t TG4GeometryManager::NofSensitiveDetectors() const
{
// Returns the total number of sensitive detectors.
// ---

  return TG4VSensitiveDetector::GetTotalNofSensitiveDetectors();
}

 
G4bool TG4GeometryManager::IsG3Volume(G4String lvName) const
{
// Returns true if the logical volume of given volumeName
// was not created by Gsposp method with a generic name 
// (name_copyNo).
// ---

  if (lvName.contains(gSeparator))
    return false;  
  else
    return true;   
}

 
void TG4GeometryManager::G4ToG3VolumeName(G4String& name) const
{
// Cuts _copyNo extension added to logical volume name in case 
// the logical volume was created by Gsposp method.
// ---

  if (name.contains(gSeparator)) 
  name = name(0,name.first(gSeparator));
}

 
const G4String& TG4GeometryManager::GetMapSecond(const G4String& name)
{ 
// Returns the second string associated with the name in
// the name map.
// ---

  return fNameMap.GetSecond(name); 
}

 
G4int TG4GeometryManager::GetMediumId(G4Material* material) const
{
// Returns the second index for materials (having its origin in
// G4 tracking media concept)
// ---

  return fMediumIdVector[material->GetIndex()];
}  


G4double TG4GeometryManager::GetEffA(G4Material* material) const
{
// Returns A or effective A=sum(pi*Ai) (if compound/mixture)
// of given material.
// ---

  G4double a = 0.;
  G4int nofElements = material->GetNumberOfElements();
  if (nofElements > 1) {
    G4String text = "Effective A for material mixture (";
    text = text + material->GetName();
    text = text + ") is used.";
    //TG4Globals::Warning(text);

    for (G4int i=0; i<nofElements; i++) {
      G4double aOfElement = material->GetElement(i)->GetA();
      G4double massFraction = material->GetFractionVector()[i];      
      a += aOfElement*massFraction /(TG3Units::AtomicWeight());
    }
  }
  else { 
    a = material->GetA();
    a /= TG3Units::AtomicWeight();
  }
  return a;
}


G4double TG4GeometryManager::GetEffZ(G4Material* material) const
{
// Returns Z or effective Z=sum(pi*Zi) (if compound/mixture)
// of given material.
// ---

  G4double z = 0.;
  G4int nofElements = material->GetNumberOfElements();
  if (nofElements > 1) {
    G4String text = "Effective Z for material mixture (";
    text = text + material->GetName();
    text = text + ") is used.";
    //TG4Globals::Warning(text);

    for (G4int i=0; i<nofElements; i++) {
      G4double zOfElement = material->GetElement(i)->GetZ();
      G4double massFraction = material->GetFractionVector()[i];
      z += zOfElement*massFraction;
    }
  }
  else { 
    z = material->GetZ(); 
  }  
  return z;
}
