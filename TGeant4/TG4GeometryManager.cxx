// $Id$
// Category: geometry
//
// Author: V. Berejnoi, I. Hrivnacova
//
// Class TG4GeometryManager
// ------------------------
// See the class description in the header file.
// C++ interface to Geant3 basic routines for building Geant4 geometry
// by V. Berejnoi, 25.2.1999;
// materials, tracking media support 
// added by I.Hrivnacova, 27.5.1999.

#include "TG4GeometryManager.h"
#include "TG4GeometryOutputManager.h"
#include "TG4GeometryServices.h"
#include "TG4Limits.h"
#include "TG4G3Units.h"
#include "TG4G3CutVector.h"
#include "TG4G3ControlVector.h"
#include "TG4Globals.h"

#include <G3toG4.hh> 
#include <G3toG4MANY.hh>
#include <G3toG4BuildTree.hh>
#include <G3VolTable.hh>
#include <G3RotTable.hh>
#include <G3EleTable.hh>
#include <G3MatTable.hh>
#include <G3MedTable.hh>
#include <G3SensVolVector.hh>

#include <G4LogicalVolumeStore.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4Element.hh> 

// extern global method from g3tog4
void G3CLRead(G4String &, char *);

TG4GeometryManager* TG4GeometryManager::fgInstance = 0;

//_____________________________________________________________________________
TG4GeometryManager::TG4GeometryManager() 
  : TG4Verbose("geometryManager"),
    fMediumCounter(0),
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

  fOutputManager 
    = new TG4GeometryOutputManager();

  fGeometryServices 
    = new TG4GeometryServices(&fMediumMap, &fNameMap);

  fgInstance = this;
      
  // instantiate the default element table
  //TG4ElementTable::Instance();
}

//_____________________________________________________________________________
TG4GeometryManager::TG4GeometryManager(const TG4GeometryManager& right)
  : TG4Verbose("geometryManager") {
// 
  TG4Globals::Exception(
    "Attempt to copy TG4GeometryManager singleton.");
}


//_____________________________________________________________________________
TG4GeometryManager::~TG4GeometryManager() {
//
  delete fOutputManager;
  delete fGeometryServices;
}

//=============================================================================
//
// operators
//
//=============================================================================


//_____________________________________________________________________________
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
 
//_____________________________________________________________________________
void TG4GeometryManager::FillMediumMap()
{
// Maps G3 tracking medium IDs to volumes names.
// ---


  static G4int done = 0;
  
  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

  for (G4int i=done; i<lvStore->size(); i++) {
    G4String name  = ((*lvStore)[i])->GetName();
    
    G4String g3Name(name);
    if (name.find("_refl")) g3Name = g3Name.substr(0, g3Name.find("_refl"));

    G4int mediumID = G3Vol.GetVTE(g3Name)->GetNmed();
    fMediumMap.Add(name, mediumID);     
  }
  
  done = lvStore->size();
}    


//=============================================================================
//
// public methods - AliMC implementation
//
//=============================================================================

 
//_____________________________________________________________________________
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
    G4double* bufin = fGeometryServices->CreateG4doubleArray(buf, nwbuf); 
    G4String namein = fGeometryServices->CutMaterialName(name);

    // write token to the output file
    if (fWriteGeometry) 
      fOutputManager->WriteGsmate(kmat, namein, a, z, dens, radl, nwbuf, bufin); 

    // create new material only if it does not yet exist
    G4Material* material = fGeometryServices->FindMaterial(a, z, dens); 
    if (material) {
      // verbose
      if (VerboseLevel() > 1) {
        G4cout << "!!! Material " << namein << " already exists as "
               << material->GetName() << G4endl;
      }	       
      G3Mat.put(kmat, material);	    
    }  	       
    else
      G4gsmate(kmat, namein, a, z, dens, radl, nwbuf, bufin); 
      
    // save the original material name
    fMaterialNameVector.push_back(namein);  

    delete [] bufin;

    if (nwbuf > 0) {  
      G4String text 
        = "TG4GeometryManager: user defined parameters for material ";
      text = text + namein;
      text = text + " are ignored by Geant4.";	
      TG4Globals::Warning(text);
    }
}
  
 
//_____________________________________________________________________________
void TG4GeometryManager::Mixture(Int_t& kmat, const char *name, Float_t *a, 
          Float_t *z, Float_t dens, Int_t nlmat, Float_t *wmat)
{ 
// Creates G4Material composed of more elements.
// !! Parameters radl, absl, buf, nwbuf are ignored in G4gsmate
// Comment: 
// absl - this parameter is ignored by GEANT3, too
// ---

   Int_t npar = abs(nlmat);
   G4double *ain = fGeometryServices->CreateG4doubleArray(a, npar); 
   G4double *zin = fGeometryServices->CreateG4doubleArray(z, npar); 
   G4double *wmatin = fGeometryServices->CreateG4doubleArray(wmat, npar); 
   G4String namein = fGeometryServices->CutMaterialName(name);

   kmat = ++fMaterialCounter;

   // write token to the output file
   if (fWriteGeometry) 
     fOutputManager->WriteGsmixt(kmat, namein, ain, zin, dens, nlmat, wmatin);

   // create new material only if it does not yet exist
   G4Material* material 
     = fGeometryServices->FindMaterial(ain, zin, dens, nlmat, wmatin); 
   if (material) {
     // verbose
     if (VerboseLevel() > 1) {
       G4cout << "!!! Material " << namein << " already exists as "
              << material->GetName() << G4endl;
     }	      
     G3Mat.put(kmat, material);	    
   }  
   else
     G4gsmixt(kmat, namein, ain, zin, dens, nlmat, wmatin);
     
    // save the original material name
    fMaterialNameVector.push_back(namein);  

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

//_____________________________________________________________________________
void TG4GeometryManager::Medium(Int_t& kmed, const char *name, Int_t nmat, 
          Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, 
          Float_t stemax, Float_t deemax, Float_t epsil, 
          Float_t stmin, Float_t* ubuf, Int_t nbuf)
{ 
// Creates a temporary "medium" that is used for 
// assigning corresponding parameters to G4 objects:
// NTMED is stored as a second material index;
// ISVOL is used for builing G3SensVolVector;
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

  G4String namein = fGeometryServices->CutMaterialName(name);

  kmed = ++fMediumCounter;

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGstmed(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, 
        stemax, deemax, epsil, stmin, 0, 0);

  G4gstmed(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
       epsil, stmin, 0, fUseG3TMLimits);
     // !! instead of the nbuf argument the bool fIsG3Default is passed

  // generate new unique name  
  G4String newName 
    = fGeometryServices
        ->GenerateLimitsName(kmed, namein, fMaterialNameVector[nmat-1]);
  fMediumNameVector.push_back(newName);
  
  if (nbuf > 0) {  
    G4String medName = name;
    G4String text
      = "TG4GeometryManager: user defined parameters for medium ";
    text = text + medName;
    text = text + " are ignored by Geant4.";  
    TG4Globals::Warning(text);
  }
} 


//_____________________________________________________________________________
void TG4GeometryManager::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
           Double_t thetaY, Double_t phiY, Double_t thetaZ, Double_t phiZ)
{
// Creates G4RotationMatrix.
// ---

  krot = ++fMatrixCounter;

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGsrotm(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ);

  G4gsrotm(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ);
}
  

//_____________________________________________________________________________
void TG4GeometryManager::Matrix(Int_t& krot, Float_t thetaX, Float_t phiX, 
           Float_t thetaY, Float_t phiY, Float_t thetaZ, Float_t phiZ)
{
// Single precision interface.
// ---

  //TG4Globals::Warning("TG4GeometryManager::Matrix in single precision.");
  
  Double_t dthetaX = thetaX;
  Double_t dphiX   = phiX; 
  Double_t dthetaY = thetaY; 
  Double_t dphiY   = phiY;
  Double_t dthetaZ = thetaZ;
  Double_t dphiZ   = phiZ;

  Matrix(krot, dthetaX, dphiX, dthetaY, dphiY, dthetaZ, dphiZ);
}
  

//_____________________________________________________________________________
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
 

//_____________________________________________________________________________
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
    a = fGeometryServices->GetEffA(material);
    z = fGeometryServices->GetEffZ(material);
    
    dens = material->GetDensity();
    dens /= TG4G3Units::MassDensity();

    radl = material->GetRadlen();
    radl /= TG4G3Units::Length();

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

 
//_____________________________________________________________________________
void  TG4GeometryManager::Gstpar(Int_t itmed, const char *param, 
           Float_t parval) 
{ 
// Write token to the output file only,
// the method is performed by TG4PhysicsManager.
// ---

  if (fWriteGeometry) 
    fOutputManager->WriteGstpar(itmed, param, parval); 
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::SetCerenkov(Int_t itmed, Int_t npckov, 
                             Float_t* ppckov, Float_t* absco, Float_t* effic, 
			     Float_t* rindex)
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

  G4double* ppckovDbl = fGeometryServices->CreateG4doubleArray(ppckov, npckov); 
  G4double* abscoDbl  = fGeometryServices->CreateG4doubleArray(absco, npckov); 
  G4double* efficDbl  = fGeometryServices->CreateG4doubleArray(effic, npckov); 
  G4double* rindexDbl = fGeometryServices->CreateG4doubleArray(rindex, npckov); 
  
  // add units
  G4int i;
  for (i=0; i<npckov; i++) {
    ppckovDbl[i] = ppckovDbl[i]*TG4G3Units::Energy();
    abscoDbl[i]  = abscoDbl[i]*TG4G3Units::Length();
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
    G4String text = "TG4GeometryManager::SetCerenkov: \n";
    text = text + "    Medium not found."; 
    G4Exception(text);
  }  
  G4Material* material = medium->GetMaterial();
  
  // set material properties table 
  material->SetMaterialPropertiesTable(table);

  // verbose
  if (VerboseLevel() > 0) {
    G4cout << "The tables for UV photon tracking set for "
           << material->GetName() << G4endl;
  }	   
  for (i=0; i<npckov; i++)
    G4cout << ppckovDbl[i] << " " << rindexDbl[i] << G4endl;
	 
  delete ppckovDbl;
  delete abscoDbl;
  delete efficDbl;
  delete rindexDbl;	 
}			 

 
//_____________________________________________________________________________
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

    G4gsdvn(fGeometryServices->CutName(name), 
            fGeometryServices->CutName(mother), ndiv, iaxis);

    // register name in name map
    fNameMap.AddName(fGeometryServices->CutName(name));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvn2(const char *name, const char *mother, 
           Int_t ndiv, Int_t iaxis, Double_t c0i, Int_t numed) 
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

    G4gsdvn2(fGeometryServices->CutName(name),
             fGeometryServices->CutName(mother), ndiv, iaxis, c0i, numed);

    // register name in name map
    fNameMap.AddName(fGeometryServices->CutName(name));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvn2(const char *name, const char *mother, 
           Int_t ndiv, Int_t iaxis, Float_t c0i, Int_t numed) 
{ 
// Single precision interface.
// ---

    //TG4Globals::Warning("TG4GeometryManager::Gsdvn2 in single precision.");

    G4double dc0i = c0i;
    
    Gsdvn2(name, mother, ndiv, iaxis, dc0i, numed);
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvt(const char *name, const char *mother, 
           Double_t step, Int_t iaxis, Int_t numed, Int_t ndvmx)
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

    G4gsdvt(fGeometryServices->CutName(name), 
            fGeometryServices->CutName(mother), step, iaxis, numed, ndvmx);

    // register name in name map
    fNameMap.AddName(fGeometryServices->CutName(name));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvt(const char *name, const char *mother, 
           Float_t step, Int_t iaxis, Int_t numed, Int_t ndvmx)
{ 
// Single precision interface.
// ---

    //TG4Globals::Warning("TG4GeometryManager::Gsdvt in single precision.");

    G4double dstep = step;

    Gsdvt(name, mother, dstep, iaxis, numed, ndvmx);
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvt2(const char *name, const char *mother, 
           Double_t step, Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx)
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

    G4gsdvt2(fGeometryServices->CutName(name), 
             fGeometryServices->CutName(mother), step, iaxis, c0, numed, ndvmx);

    // register name in name map
    fNameMap.AddName(fGeometryServices->CutName(name));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsdvt2(const char *name, const char *mother, 
           Float_t step, Int_t iaxis, Float_t c0, Int_t numed, Int_t ndvmx)
{ 
// Single precision interface.
// ---

    //TG4Globals::Warning("TG4GeometryManager::Gsdvt2 in single precision.");
    
    G4double dstep = step;
    G4double dc0   = c0;

    Gsdvt2(name, mother, dstep, iaxis, dc0, numed, ndvmx);
} 
 
 
//_____________________________________________________________________________
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
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gspos(const char *vname, Int_t num, 
          const char *vmoth, Double_t x, Double_t y, Double_t z, Int_t irot, 
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

   G4gspos(fGeometryServices->CutName(vname), num,
           fGeometryServices->CutName(vmoth), x, y, z, irot, vonly);

   // register name in name map
   fNameMap.AddName(fGeometryServices->CutName(vname));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gspos(const char *vname, Int_t num, 
          const char *vmoth, Float_t x, Float_t y, Float_t z, Int_t irot, 
          const char *vonly) 
{ 
// Single precision interface.
// ---

   //TG4Globals::Warning("TG4GeometryManager::Gspos in single precision.");

   G4double dx = x; 
   G4double dy = y;
   G4double dz = z;
   
   Gspos(vname, num, vmoth, dx, dy, dz, irot, vonly);
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsposp(const char *name, Int_t nr, 
           const char *mother, Double_t x, Double_t y, Double_t z, Int_t irot, 
           const char *konly, Double_t *upar, Int_t np ) 
{ 
//  Geant3 desription:
//  ==================
//      Place a copy of generic volume NAME with user number
//      NR inside MOTHER, with its parameters UPAR(1..NP)
// ---

   // write token to the output file
   if (fWriteGeometry) 
     fOutputManager->WriteGsposp(name, nr, mother, x, y, z, irot, konly, upar, np);

   G4gsposp(fGeometryServices->CutName(name), nr, 
            fGeometryServices->CutName(mother), x, y, z, irot, konly, 
             upar, np);

   // register name in name map
   fNameMap.AddName(fGeometryServices->CutName(name));
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsposp(const char *name, Int_t nr, 
           const char *mother, Float_t x, Float_t y, Float_t z, Int_t irot, 
           const char *konly, Float_t *upar, Int_t np ) 
{ 
// Single precision interface.
// ---

   //TG4Globals::Warning("TG4GeometryManager::Gsposp in single precision.");

   G4double dx = x; 
   G4double dy = y;
   G4double dz = z;
   G4double* parin = fGeometryServices->CreateG4doubleArray(upar, np); 

   Gsposp(name, nr, mother, dx, dy, dz, irot, konly, parin, np);

   delete [] parin;
} 
 
 
//_____________________________________________________________________________
void  TG4GeometryManager::Gsbool(const char* onlyVolName, 
                                 const char* manyVolName)
{ 
// Helps for resolving MANY.
// Specifies the ONLY volume that overlaps with the 
// specified MANY and has to be substracted.
// ---  

   // write token to the output file
   //if (fWriteGeometry) 
   //  fOutputManager->WriteGsbool(onlyVolName, manyVolName);

   G4gsbool(onlyVolName, manyVolName);
} 
 
 
//_____________________________________________________________________________
Int_t TG4GeometryManager::Gsvolu(const char *name, const char *shape, 
          Int_t nmed, Double_t *upar, Int_t npar) 
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

  // write token to the output file
  if (fWriteGeometry) 
    fOutputManager->WriteGsvolu(name, shape, nmed, upar, npar);    

  G4gsvolu(fGeometryServices->CutName(name), 
           fGeometryServices->CutName(shape), nmed, upar, npar);
  
  // register name in name map
  fNameMap.AddName(fGeometryServices->CutName(name));

  return 0;
} 


//_____________________________________________________________________________
Int_t TG4GeometryManager::Gsvolu(const char *name, const char *shape, 
          Int_t nmed, Float_t *upar, Int_t npar) 
{ 
// Single precision interface.
// ---

  //TG4Globals::Warning("TG4GeometryManager::Gsvolu in single precision.");

  G4double* parin = fGeometryServices->CreateG4doubleArray(upar, npar); 

  G4int result
   = Gsvolu(name, shape, nmed, parin, npar);

  delete [] parin;

  return result;
} 


//_____________________________________________________________________________
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

 
//=============================================================================
//
// public methods - Geant4 only
//
//=============================================================================

 
//_____________________________________________________________________________
G4VPhysicalVolume* TG4GeometryManager::CreateG4Geometry()
{
// Creates G4 geometry objects according to the G3VolTable 
// and returns the world physical volume.
// ---

  // set the first entry in the G3Vol table
  Ggclos();
  G3VolTableEntry* first = G3Vol.GetFirstVTE();
  
  // transform MANY to Boolean solids
  G3toG4MANY(first);
  
  // create G4 geometry
  G3toG4BuildTree(first,0);  
  
  // fill medium map
  FillMediumMap();

  // print G3 volume table statistics
  G3Vol.VTEStat();

  // print G4 geometry statistics
  if (VerboseLevel() > 0) {
    G4cout << "G4 Stat: instantiated " 
           << fGeometryServices->NofG4LogicalVolumes()  
	   << " logical volumes \n"
	   << "                      " 
	   << fGeometryServices->NofG4PhysicalVolumes() 
	   << " physical volumes" << G4endl;
  }	   

  // position the first entry 
  // (in Geant3 the top volume cannot be positioned)
  // 
  if (!fGeometryServices->GetWorld()) {
    G4VPhysicalVolume* world
       = new G4PVPlacement(0, G4ThreeVector(), first->GetName(), 
                           first->GetLV(), 0, false, 0);
    fGeometryServices->SetWorld(world);			    
  }
  return fGeometryServices->GetWorld();		      
}

 
//_____________________________________________________________________________
void TG4GeometryManager::SetUserLimits(const TG4G3CutVector& cuts,
                               const TG4G3ControlVector& controls) const
{
// Sets user limits defined in G3MedTable for all logical volumes.
// ---

  G4LogicalVolumeStore* lvStore = G4LogicalVolumeStore::GetInstance();

  for (G4int i=0; i<lvStore->size(); i++) {
    G4LogicalVolume* lv = (*lvStore)[i];

    // get limits from G3Med
    G4int mediumIndex = fGeometryServices->GetMediumId(lv);    
    G4UserLimits* limits = G3Med.get(mediumIndex)->GetLimits();
    TG4Limits* tg4Limits = fGeometryServices->GetLimits(limits);

    // get tracking medium name
    G4String name = fMediumNameVector[mediumIndex-1];
    
    if (tg4Limits) 
      tg4Limits->SetName(name);
    else {
      tg4Limits = fGeometryServices->FindLimits(name, true);  
      if (!tg4Limits) 
        tg4Limits = new TG4Limits(name, cuts, controls); 
    }
    
    // update controls in limits according to the setup 
    // in the passed vector
    tg4Limits->Update(controls);

    // set limits to logical volume
    lv->SetUserLimits(tg4Limits);
    //NO TG4 lv->SetUserLimits(0);
  } 
}

//_____________________________________________________________________________
void TG4GeometryManager::ReadG3Geometry(G4String filePath)
{
// Processes g3calls.dat file and fills G3 tables.
// ---

  // verbose
  if (VerboseLevel() > 0) {
    G4cout << "Reading the call list file " << filePath << "..." << G4endl;
  }
    
  G3CLRead(filePath, NULL);

  if (VerboseLevel() > 0) {
    G4cout << "Call list file read completed. Build geometry" << G4endl;
  }  
}

 
//_____________________________________________________________________________
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

 
//_____________________________________________________________________________
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
  //G3Rot.Clear();
  G3SensVol.clear(); 
}

 
//_____________________________________________________________________________
void TG4GeometryManager::ClearG3TablesFinal()
{
// Clears G3 medias and volumes tables
// (the top volume is removed from the vol table)
// ---

  G3Mat.Clear();
  G3Med.Clear();
  G3Vol.Clear();  
}

 
//_____________________________________________________________________________
void TG4GeometryManager::OpenOutFile(G4String filePath)
{ 
// Opens output file.
// ---

  fOutputManager->OpenFile(filePath);
}

 
//_____________________________________________________________________________
void TG4GeometryManager::CloseOutFile()
{ 
// Closes output file.
// ---

  fOutputManager->CloseFile();
}

 
//_____________________________________________________________________________
void TG4GeometryManager::SetWriteGeometry(G4bool writeGeometry)
{ 
// Controls geometry output.
// ---

  fWriteGeometry = writeGeometry; 
}

 
//_____________________________________________________________________________
void TG4GeometryManager::SetMapSecond(const G4String& name)
{
// Sets the second name for the map of volumes names.
// ---

  fNameMap.SetSecond(name);
}
