// @(#):$Name$:$Id$
// Author: Andrei Gheata 10/07/2003

#include "TObjString.h"
#include "TFluka.h"
//#include "TVirtualMCApplication.h"
#include "TFlukaMCGeometry.h"
#include "TGeoManager.h" 
#include "TGeoVolume.h" 

#include "TCallf77.h"

#ifndef WIN32 
# define idnrwr idnrwr_
# define g1wr   g1wr_
# define g1rtwr g1rtwr_
# define conhwr conhwr_
# define inihwr inihwr_
# define jomiwr jomiwr_
# define lkdbwr lkdbwr_
# define lkfxwr lkfxwr_
# define lkmgwr lkmgwr_
# define lkwr lkwr_
# define magfld magfld_
# define nrmlwr nrmlwr_
# define rgrpwr rgrpwr_
# define isvhwr isvhwr_

#else

# define idnrwr IDNRWR
# define g1wr   G1WR
# define g1rtwr G1RTWR
# define conhwr CONHWR
# define inihwr INIHWR
# define jomiwr JOMIWR
# define lkdbwr LKDBWR
# define lkfxwr LKFXWR
# define lkmgwr LKMGWR
# define lkwr   LKWR
# define magfld MAGFLD
# define nrmlwr NRMLWR
# define rgrpwr RGRPWR
# define isvhwr ISVHWR

#endif

//____________________________________________________________________________ 
extern "C" 
{
   //
   // Prototypes for FLUKA navigation methods
   //
   Int_t type_of_call idnrwr(const Int_t & /*nreg*/, const Int_t & /*mlat*/);
   void  type_of_call   g1wr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/, 
                             Double_t * /*pV*/,  Int_t & /*oldReg*/ , const Int_t & /*oldLttc*/, Double_t & /*propStep*/,
                             Int_t & /*nascFlag*/, Double_t & /*retStep*/, Int_t & /*newReg*/,
	                          Double_t & /*saf*/, Int_t & /*newLttc*/, Int_t & /*LttcFlag*/,
                             Double_t *s /*Lt*/, Int_t * /*jrLt*/);
   
   void  type_of_call g1rtwr();
   void  type_of_call conhwr(Int_t & /*intHist*/, Int_t * /*incrCount*/); 
   void  type_of_call inihwr(Int_t & /*intHist*/);
   void  type_of_call jomiwr(const Int_t & /*nge*/, const Int_t & /*lin*/, const Int_t & /*lou*/,
                             Int_t & /*flukaReg*/);
   void  type_of_call lkdbwr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/,
                             Double_t * /*pV*/, const Int_t & /*oldReg*/, const Int_t & /*oldLttc*/,
                             Int_t & /*newReg*/, Int_t & /*flagErr*/, Int_t & /*newLttc*/);
   void  type_of_call lkfxwr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/,
                             Double_t * /*pV*/, const Int_t & /*oldReg*/, const Int_t & /*oldLttc*/,
                             Int_t & /*newReg*/, Int_t & /*flagErr*/, Int_t & /*newLttc*/);
   void  type_of_call lkmgwr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/,
                             Double_t * /*pV*/, const Int_t & /*oldReg*/, const Int_t & /*oldLttc*/,
		                       Int_t & /*flagErr*/, Int_t & /*newReg*/, Int_t & /*newLttc*/);
   void  type_of_call   lkwr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/,
                             Double_t * /*pV*/, const Int_t & /*oldReg*/, const Int_t & /*oldLttc*/,
	                          Int_t & /*newReg*/, Int_t & /*flagErr*/, Int_t & /*newLttc*/);
//   void  type_of_call magfld(const Double_t & /*pX*/, const Double_t & /*pY*/, const Double_t & /*pZ*/,
//                             Double_t & /*cosBx*/, Double_t & /*cosBy*/, Double_t & /*cosBz*/, 
//                             Double_t & /*Bmag*/, Int_t & /*reg*/, Int_t & /*idiscflag*/);	    
   void  type_of_call nrmlwr(Double_t & /*pSx*/, Double_t & /*pSy*/, Double_t & /*pSz*/,
                             Double_t & /*pVx*/, Double_t & /*pVy*/, Double_t & /*pVz*/,
	                          Double_t * /*norml*/, const Int_t & /*oldReg*/, 
	                          const Int_t & /*newReg*/, Int_t & /*flagErr*/);
   void  type_of_call rgrpwr(const Int_t & /*flukaReg*/, const Int_t & /*ptrLttc*/, Int_t & /*g4Reg*/,
                             Int_t * /*indMother*/, Int_t * /*repMother*/, Int_t & /*depthFluka*/);
   Int_t type_of_call isvhwr(const Int_t & /*fCheck*/, const Int_t & /*intHist*/);
};
   
// TFluka global pointer
TFluka *fluka = 0;
TFlukaMCGeometry *mcgeom = 0;

ClassImp(TFlukaMCGeometry)

TFlukaMCGeometry* TFlukaMCGeometry::fgInstance=0;

//_____________________________________________________________________________
TFlukaMCGeometry::TFlukaMCGeometry(const char *name, const char *title) 
  : TVirtualMCGeometry(name, title)
{
  //
  // Standard constructor
  //
  fDebug        = kFALSE;
  fLastMaterial = 0;
  fNextRegion   = 0;
  fNextLattice  = 0;
  fRegionList   = 0;
  fluka = (TFluka*)gMC;
  mcgeom = this;
}

//_____________________________________________________________________________
TFlukaMCGeometry::TFlukaMCGeometry()
  : TVirtualMCGeometry()
{    
  //
  // Default constructor
  //
  fDebug        = kFALSE;
  fLastMaterial = 0;
  fNextRegion   = 0;
  fNextLattice  = 0;
  fRegionList   = 0;
  fluka = (TFluka*)gMC;
  mcgeom = this;
}

//_____________________________________________________________________________
TFlukaMCGeometry::~TFlukaMCGeometry() 
{
  //
  // Destructor
  //
  fgInstance=0;
  if (fRegionList) delete [] fRegionList;
  if (gGeoManager) delete gGeoManager;
}

//
// private methods
//
//_____________________________________________________________________________
TFlukaMCGeometry::TFlukaMCGeometry(const TFlukaMCGeometry &)
  : TVirtualMCGeometry()
{    
  //
  // Copy constructor
  //
}

//_____________________________________________________________________________
Double_t* TFlukaMCGeometry::CreateDoubleArray(Float_t* array, Int_t size) const
{
// Converts Float_t* array to Double_t*,
// !! The new array has to be deleted by user.
// ---

  Double_t* doubleArray;
  if (size>0) {
    doubleArray = new Double_t[size]; 
    for (Int_t i=0; i<size; i++) doubleArray[i] = array[i];
  }
  else {
    //doubleArray = 0; 
    doubleArray = new Double_t[1]; 
  }  
  return doubleArray;
}
//
// public methods
//_____________________________________________________________________________
void TFlukaMCGeometry::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		       Float_t &dens, Float_t &radl, Float_t &absl,
		       Float_t* /*ubuf*/, Int_t& /*nbuf*/)
{
   if (fDebug) printf("Gfmate %i\n", imat);
   TGeoMaterial *mat;
   TIter next (gGeoManager->GetListOfMaterials());
   while ((mat = (TGeoMaterial*)next())) {
     if (mat->GetUniqueID() == (UInt_t)imat) break;
   }
   if (!mat) {
      Error("Gfmate", "no material with index %i found", imat);
      return;
   }
   sprintf(name, "%s", mat->GetName());
   a = mat->GetA();
   z = mat->GetZ();
   dens = mat->GetDensity();
   radl = mat->GetRadLen();
   absl = mat->GetIntLen();
   if (fDebug) printf("   ->material found : %s a=%g, z=%g, dens=%g, radl=%g, absl=%g\n", name, a,z,dens,radl,absl);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
		       Double_t &dens, Double_t &radl, Double_t &absl,
		       Double_t* /*ubuf*/, Int_t& /*nbuf*/)
{
   if (fDebug) printf("Gfmate %i\n", imat);
   TGeoMaterial *mat;
   TIter next (gGeoManager->GetListOfMaterials());
   while ((mat = (TGeoMaterial*)next())) {
     if (mat->GetUniqueID() == (UInt_t)imat) break;
   }
   if (!mat) {
      Error("Gfmate", "no material with index %i found", imat);
      return;
   }
   sprintf(name, "%s", mat->GetName());
   a = mat->GetA();
   z = mat->GetZ();
   dens = mat->GetDensity();
   radl = mat->GetRadLen();
   absl = mat->GetIntLen();
   if (fDebug) printf("   ->material found : %s a=%g, z=%g, dens=%g, radl=%g, absl=%g\n", name, a,z,dens,radl,absl);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
		       Double_t dens, Double_t radl, Double_t absl, Float_t* buf,
		       Int_t nwbuf)
{
  //
  // Defines a Material
  // 
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //
  
  Double_t* dbuf = CreateDoubleArray(buf, nwbuf);  
  Material(kmat, name, a, z, dens, radl, absl, dbuf, nwbuf);
  delete [] dbuf;
}  

//_____________________________________________________________________________
void TFlukaMCGeometry::Material(Int_t& kmat, const char* name, Double_t a, Double_t z,
		       Double_t dens, Double_t radl, Double_t absl, Double_t* /*buf*/,
		       Int_t /*nwbuf*/)
{
  //
  // Defines a Material
  // 
  //  kmat               number assigned to the material
  //  name               material name
  //  a                  atomic mass in au
  //  z                  atomic number
  //  dens               density in g/cm3
  //  absl               absorbtion length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -absl is taken
  //  radl               radiation length in cm
  //                     if >=0 it is ignored and the program 
  //                     calculates it, if <0. -radl is taken
  //  buf                pointer to an array of user words
  //  nbuf               number of user words
  //

  kmat = gGeoManager->GetListOfMaterials()->GetSize();
  gGeoManager->Material(name, a, z, dens, kmat, radl, absl);
  if (fDebug) printf("Material %s: kmat=%i, a=%g, z=%g, dens=%g\n", name, kmat, a, z, dens);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Mixture(Int_t& kmat, const char* name, Float_t* a, Float_t* z, 
		      Double_t dens, Int_t nlmat, Float_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by 
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture. 
  // 
  // If nlmat < 0 then WMAT contains the number of atoms 
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //
  
  Double_t* da = CreateDoubleArray(a, TMath::Abs(nlmat));  
  Double_t* dz = CreateDoubleArray(z, TMath::Abs(nlmat));  
  Double_t* dwmat = CreateDoubleArray(wmat, TMath::Abs(nlmat));  

  Mixture(kmat, name, da, dz, dens, nlmat, dwmat);
  for (Int_t i=0; i<nlmat; i++) {
    a[i] = da[i]; z[i] = dz[i]; wmat[i] = dwmat[i];
  }  

  delete [] da;
  delete [] dz;
  delete [] dwmat;
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Mixture(Int_t& kmat, const char* name, Double_t* a, Double_t* z, 
		      Double_t dens, Int_t nlmat, Double_t* wmat)
{
  //
  // Defines mixture OR COMPOUND IMAT as composed by 
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture. 
  // 
  // If nlmat < 0 then WMAT contains the number of atoms 
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //

  if (nlmat < 0) {
     nlmat = - nlmat;
     Double_t amol = 0;
     Int_t i;
     for (i=0;i<nlmat;i++) {
        amol += a[i]*wmat[i];
     }
     for (i=0;i<nlmat;i++) {
        wmat[i] *= a[i]/amol;
     }
  }
  kmat = gGeoManager->GetListOfMaterials()->GetSize();
  if (fDebug) {
     printf("Mixture %s with %i elem: kmat=%i, dens=%g\n", name, nlmat, kmat, dens);
     for (Int_t j=0; j<nlmat; j++) printf("  Elem %i: z=%g  a=%g  w=%g\n",j,z[j],a[j],wmat[j]);
  }   
  gGeoManager->Mixture(name, a, z, dens, nlmat, wmat, kmat);
}
//_____________________________________________________________________________
Int_t TFlukaMCGeometry::GetMedium() const
{
// Get current medium number
   Int_t imed = 0;
   TGeoNode *node = gGeoManager->GetCurrentNode();
   if (!node) imed = gGeoManager->GetTopNode()->GetVolume()->GetMedium()->GetId();
   else       imed = node->GetVolume()->GetMedium()->GetId();
   if (fDebug) printf("GetMedium=%i\n", imed);
   return imed;
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::GetFlukaMaterial(Int_t imed) const
{
// Returns FLUKA material index for medium IMED
   TGeoMedium *med = (TGeoMedium*)gGeoManager->GetListOfMedia()->At(imed-1);
   if (!med) {
      Error("GetFlukaMaterial", "MEDIUM %i nor found", imed);
      return -1;
   }
   Int_t imatfl = med->GetMaterial()->GetIndex();
   return imatfl;   
}

//_____________________________________________________________________________
Int_t *TFlukaMCGeometry::GetRegionList(Int_t imed, Int_t &nreg)
{
// Get an ordered list of regions matching a given medium number
   nreg = 0;
   if (!fRegionList) fRegionList = new Int_t[NofVolumes()+1];
   TIter next(gGeoManager->GetListOfUVolumes());
   TGeoVolume *vol;
   Int_t imedium, ireg;
   while ((vol = (TGeoVolume*)next())) {
      imedium = vol->GetMedium()->GetId();
      if (imedium == imed) {
         ireg = vol->GetNumber();
         fRegionList[nreg++] = ireg;
      }
   }
   return fRegionList;
}         

//_____________________________________________________________________________
Int_t *TFlukaMCGeometry::GetMaterialList(Int_t imat, Int_t &nreg)
{
// Get an ordered list of regions matching a given medium number
   nreg = 0;
   if (!fRegionList) fRegionList = new Int_t[NofVolumes()+1];
   TIter next(gGeoManager->GetListOfUVolumes());
   TGeoVolume *vol;
   Int_t imaterial, ireg;
   while ((vol = (TGeoVolume*)next())) {
      imaterial = vol->GetMedium()->GetMaterial()->GetIndex();
      if (imaterial == imat) {
         ireg = vol->GetNumber();
         fRegionList[nreg++] = ireg;
      }
   }
   return fRegionList;
}         
//_____________________________________________________________________________
void TFlukaMCGeometry::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
		     Int_t ifield, Double_t fieldm, Double_t tmaxfd,
		     Double_t stemax, Double_t deemax, Double_t epsil,
		     Double_t stmin, Float_t* ubuf, Int_t nbuf)
{
  //
  //  kmed      tracking medium number assigned
  //  name      tracking medium name
  //  nmat      material number
  //  isvol     sensitive volume flag
  //  ifield    magnetic field
  //  fieldm    max. field value (kilogauss)
  //  tmaxfd    max. angle due to field (deg/step)
  //  stemax    max. step allowed
  //  deemax    max. fraction of energy lost in a step
  //  epsil     tracking precision (cm)
  //  stmin     min. step due to continuous processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
  //  performed with g3helix; ifield = 3 if tracking performed with g3helx3.
  //  

  //printf("Creating mediuma: %s, numed=%d, nmat=%d\n",name,kmed,nmat);
  Double_t* dubuf = CreateDoubleArray(ubuf, nbuf);  
  Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, epsil,
         stmin, dubuf, nbuf);
  delete [] dubuf;	 
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Medium(Int_t& kmed, const char* name, Int_t nmat, Int_t isvol,
		     Int_t ifield, Double_t fieldm, Double_t tmaxfd,
		     Double_t stemax, Double_t deemax, Double_t epsil,
		     Double_t stmin, Double_t* /*ubuf*/, Int_t /*nbuf*/)
{
  //
  //  kmed      tracking medium number assigned
  //  name      tracking medium name
  //  nmat      material number
  //  isvol     sensitive volume flag
  //  ifield    magnetic field
  //  fieldm    max. field value (kilogauss)
  //  tmaxfd    max. angle due to field (deg/step)
  //  stemax    max. step allowed
  //  deemax    max. fraction of energy lost in a step
  //  epsil     tracking precision (cm)
  //  stmin     min. step due to continuos processes (cm)
  //
  //  ifield = 0 if no magnetic field; ifield = -1 if user decision in guswim;
  //  ifield = 1 if tracking performed with g3rkuta; ifield = 2 if tracking
  //  performed with g3helix; ifield = 3 if tracking performed with g3helx3.
  //  

  kmed = gGeoManager->GetListOfMedia()->GetSize()+1;
  gGeoManager->Medium(name,kmed,nmat, isvol, ifield, fieldm, tmaxfd, stemax,deemax, epsil, stmin);
  if (fDebug) printf("Medium %s: kmed=%i, nmat=%i, isvol=%i\n", name, kmed, nmat,isvol);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Matrix(Int_t& krot, Double_t thex, Double_t phix, Double_t they,
		     Double_t phiy, Double_t thez, Double_t phiz)
{
  //
  //  krot     rotation matrix number assigned
  //  theta1   polar angle for axis i
  //  phi1     azimuthal angle for axis i
  //  theta2   polar angle for axis ii
  //  phi2     azimuthal angle for axis ii
  //  theta3   polar angle for axis iii
  //  phi3     azimuthal angle for axis iii
  //
  //  it defines the rotation matrix number irot.
  //  

  krot = gGeoManager->GetListOfMatrices()->GetEntriesFast();
  gGeoManager->Matrix(krot, thex, phix, they, phiy, thez, phiz);  
  if (fDebug) printf("Rotation %i defined\n", krot);
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		      Float_t *upar, Int_t npar) 
{ 
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //  

  Double_t* dupar = CreateDoubleArray(upar, npar);
  Int_t id = Gsvolu(name, shape, nmed, dupar, npar);
  delete [] dupar;  
  return id;
} 

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		      Double_t *upar, Int_t npar) 
{ 
  //
  //  NAME   Volume name
  //  SHAPE  Volume type
  //  NUMED  Tracking medium number
  //  NPAR   Number of shape parameters
  //  UPAR   Vector containing shape parameters
  //
  //  It creates a new volume in the JVOLUM data structure.
  //  
  char vname[5];
  Vname(name,vname);
  char vshape[5];
  Vname(shape,vshape);

  TGeoVolume* vol = gGeoManager->Volume(vname, shape, nmed, upar, npar); 
  if (fDebug) printf("Volume %s: id=%i shape=%s, nmed=%i\n", vname, vol->GetNumber(), shape, nmed);
  return vol->GetNumber();
} 
 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsdvn(const char *name, const char *mother, Int_t ndiv,
		     Int_t iaxis) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //  NAME   Volume name
  //  MOTHER Mother volume name
  //  NDIV   Number of divisions
  //  IAXIS  Axis value
  //
  //  X,Y,Z of CAXIS will be translated to 1,2,3 for IAXIS.
  //  It divides a previously defined volume.
  //  
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  gGeoManager->Division(vname, vmother, iaxis, ndiv, 0, 0, 0, "n");
  if (fDebug) printf("Division %s: mother=%s iaxis=%i ndiv=%i\n", vname, vmother, iaxis, ndiv);
} 
 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsdvn2(const char *name, const char *mother, Int_t ndiv,
		      Int_t iaxis, Double_t c0i, Int_t numed) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  // Divides mother into ndiv divisions called name
  // along axis iaxis starting at coordinate value c0.
  // the new volume created will be medium number numed.
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  gGeoManager->Division(vname, vmother, iaxis, ndiv, c0i, 0, numed, "nx");
} 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsdvt(const char *name, const char *mother, Double_t step,
		     Int_t iaxis, Int_t numed, Int_t /*ndvmx*/) 
{ 
  //
  // Create a new volume by dividing an existing one
  // 
  //       Divides MOTHER into divisions called NAME along
  //       axis IAXIS in steps of STEP. If not exactly divisible 
  //       will make as many as possible and will centre them 
  //       with respect to the mother. Divisions will have medium 
  //       number NUMED. If NUMED is 0, NUMED of MOTHER is taken.
  //       NDVMX is the expected maximum number of divisions
  //          (If 0, no protection tests are performed) 
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  
  gGeoManager->Division(vname, vmother, iaxis, 0, 0, step, numed, "s");
} 

//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsdvt2(const char *name, const char *mother, Double_t step,
		      Int_t iaxis, Double_t c0, Int_t numed, Int_t /*ndvmx*/) 
{ 
  //
  // Create a new volume by dividing an existing one
  //                                                                    
  //           Divides MOTHER into divisions called NAME along          
  //            axis IAXIS starting at coordinate value C0 with step    
  //            size STEP.                                              
  //           The new volume created will have medium number NUMED.    
  //           If NUMED is 0, NUMED of mother is taken.                 
  //           NDVMX is the expected maximum number of divisions        
  //             (If 0, no protection tests are performed)              
  //
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  
  gGeoManager->Division(vname, vmother, iaxis, 0, c0, step, numed, "sx");
} 

//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsord(const char * /*name*/, Int_t /*iax*/) 
{ 
  //
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
  //

  // TBC - keep this function
  // nothing to be done for TGeo  //xx
} 
 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gspos(const char *name, Int_t nr, const char *mother, Double_t x,
		     Double_t y, Double_t z, Int_t irot, const char *konly) 
{ 
  //
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
  //  
    
  TString only = konly;
  only.ToLower();
  Bool_t isOnly = kFALSE;
  if (only.Contains("only")) isOnly = kTRUE;
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);
  
  Double_t *upar=0;
  gGeoManager->Node(vname, nr, vmother, x, y, z, irot, isOnly, upar);
  if (fDebug) printf("Adding daughter %s to %s: cpy=%i irot=%i only=%s\n", vname,vmother,nr,irot,only.Data());
} 
 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsposp(const char *name, Int_t nr, const char *mother,  
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Float_t *upar, Int_t np ) 
{ 
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

  Double_t* dupar = CreateDoubleArray(upar, np);
  Gsposp(name, nr, mother, x, y, z, irot, konly, dupar, np); 
  delete [] dupar;
} 
 
//_____________________________________________________________________________
void  TFlukaMCGeometry::Gsposp(const char *name, Int_t nr, const char *mother,  
		      Double_t x, Double_t y, Double_t z, Int_t irot,
		      const char *konly, Double_t *upar, Int_t np ) 
{ 
  //
  //      Place a copy of generic volume NAME with user number
  //      NR inside MOTHER, with its parameters UPAR(1..NP)
  //

  TString only = konly;
  only.ToLower();
  Bool_t isOnly = kFALSE;
  if (only.Contains("only")) isOnly = kTRUE;
  char vname[5];
  Vname(name,vname);
  char vmother[5];
  Vname(mother,vmother);

  gGeoManager->Node(vname,nr,vmother, x,y,z,irot,isOnly,upar,np);
  if (fDebug) printf("Adding daughter(s) %s to %s: cpy=%i irot=%i only=%s\n", vname,vmother,nr,irot,only.Data());
} 
 
//_____________________________________________________________________________
Int_t TFlukaMCGeometry::VolId(const Text_t *name) const
{
  //
  // Return the unique numeric identifier for volume name
  //

  Int_t uid = gGeoManager->GetUID(name);
  if (uid<0) {
     printf("VolId: Volume %s not found\n",name);
     return 0;
  }
  if (fDebug) printf("VolId for %s: %i\n", name, uid);
  return uid;     
}

//_____________________________________________________________________________
const char* TFlukaMCGeometry::VolName(Int_t id) const
{
  //
  // Return the volume name given the volume identifier
  //
  TGeoVolume *volume = gGeoManager->GetVolume(id);
  if (!volume) {
     Error("VolName","volume with id=%d does not exist",id);
     return "NULL";
  }
  if (fDebug) printf("VolName for id=%i: %s\n", id, volume->GetName());
  return volume->GetName();
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::NofVolumes() const 
{
  //
  // Return total number of volumes in the geometry
  //

  return gGeoManager->GetListOfUVolumes()->GetEntriesFast()-1;
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::VolId2Mate(Int_t id) const 
{
  //
  // Return material number for a given volume id
  //
  TGeoVolume *volume = gGeoManager->GetVolume(id);
  if (!volume) {
     Error("VolId2Mate","volume with id=%d does not exist",id);
     return 0;
  }
  TGeoMedium *med = volume->GetMedium();
  if (!med) return 0;
  if (fDebug) printf("VolId2Mate id=%i: idmed=%i\n", id, med->GetId());
  return med->GetId();
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::CurrentVolID(Int_t& copyNo) const
{
  // Returns the current volume ID and copy number
  if (gGeoManager->IsOutside()) return 0;
  TGeoNode *node = gGeoManager->GetCurrentNode();
  copyNo = node->GetNumber();
  Int_t id = node->GetVolume()->GetNumber();
  if (fDebug) printf("CurrentVolId(cpy=%i) = %i\n", copyNo, id); 
  return id;
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
  // Return the current volume "off" upward in the geometrical tree 
  // ID and copy number
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolID(copyNo);
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  copyNo = node->GetNumber();
  if (fDebug) printf("CurrentVolOffId(off=%i,cpy=%i) = %i\n", off,copyNo,node->GetVolume()->GetNumber() ); 
  return node->GetVolume()->GetNumber();
}
// FLUKA specific

//_____________________________________________________________________________
const char* TFlukaMCGeometry::CurrentVolName() const
{
  //
  // Returns the current volume name
  //
  if (gGeoManager->IsOutside()) return 0;
  if (fDebug) printf("CurrentVolName : %s\n", gGeoManager->GetCurrentVolume()->GetName()); 
  return gGeoManager->GetCurrentVolume()->GetName();
}
//_____________________________________________________________________________
const char* TFlukaMCGeometry::CurrentVolOffName(Int_t off) const
{
  //
  // Return the current volume "off" upward in the geometrical tree 
  // ID, name and copy number
  // if name=0 no name is returned
  //
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolName();
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  if (fDebug) printf("CurrentVolOffName(off=%i) : %s\n", off,node->GetVolume()->GetName()); 
  return node->GetVolume()->GetName();
}
  
//_____________________________________________________________________________
void TFlukaMCGeometry::Gsatt(const char *name, const char *att, Int_t val)
{ 
  //
  //  NAME   Volume name
  //  IOPT   Name of the attribute to be set
  //  IVAL   Value to which the attribute is to be set
  // see: TFluka::Gsatt
  char vname[5];
  Vname(name,vname);
  char vatt[5];
  Vname(att,vatt);
  gGeoManager->SetVolumeAttribute(vname, vatt, val);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Gdtom(Float_t *xd, Float_t *xm, Int_t iflag) 
{ 
  //
  //  Computes coordinates XM (Master Reference System
  //  knowing the coordinates XD (Detector Ref System)
  //  The local reference system can be initialized by
  //    - the tracking routines and GDTOM used in GUSTEP
  //    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
  //        (inverse routine is GMTOD)
  // 
  //   If IFLAG=1  convert coordinates
  //      IFLAG=2  convert direction cosinus
  //
   Double_t XM[3], XD[3];
   Int_t i;
   for (i=0;i<3;i++) XD[i] = xd[i];
   if (iflag == 1) gGeoManager->LocalToMaster(XD,XM);
   else            gGeoManager->LocalToMasterVect(XD,XM);
   for (i=0;i<3;i++) xm[i]=XM[i];
}   

//_____________________________________________________________________________
void TFlukaMCGeometry::Gdtom(Double_t *xd, Double_t *xm, Int_t iflag) 
{ 
   if (iflag == 1) gGeoManager->LocalToMaster(xd,xm);
   else            gGeoManager->LocalToMasterVect(xd,xm);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::Gmtod(Float_t *xm, Float_t *xd, Int_t iflag) 
{ 
  //
  //       Computes coordinates XD (in DRS) 
  //       from known coordinates XM in MRS 
  //       The local reference system can be initialized by
  //         - the tracking routines and GMTOD used in GUSTEP
  //         - a call to GMEDIA(XM,NUMED,CHECK)
  //         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
  //             (inverse routine is GDTOM) 
  //
  //        If IFLAG=1  convert coordinates 
  //           IFLAG=2  convert direction cosinus
  //
   Double_t XM[3], XD[3];
   Int_t i;
   for (i=0;i<3;i++) XM[i]=xm[i];
   if (iflag == 1) gGeoManager->MasterToLocal(XM,XD);
   else            gGeoManager->MasterToLocalVect(XM,XD);
   for (i=0;i<3;i++) xd[i] = XD[i];
}
  
//_____________________________________________________________________________
void TFlukaMCGeometry::Gmtod(Double_t *xm, Double_t *xd, Int_t iflag) 
{ 
   if (iflag == 1) gGeoManager->MasterToLocal(xm,xd);
   else            gGeoManager->MasterToLocalVect(xm,xd);
}
   
//_____________________________________________________________________________
void TFlukaMCGeometry::CreateFlukaMatFile(const char *fname)
{
  // ==== from FLUGG ====
  // NAMES OF ELEMENTS AND COMPOUNDS: the names must be written in upper case,
  // according to the fluka standard. In addition,. they must be equal to the
  // names of the fluka materials - see fluka manual - in order that the 
  // program load the right cross sections, and equal to the names included in
  // the .pemf. Otherwise the user must define the LOW-MAT CARDS, and make his
  // own .pemf, in order to get the right cross sections loaded in memory.
  
   Int_t zelem[128];
   static char elNames[220] = {
   //  1 ============================= 5 ==================================== 10 ===================================== 15 ===
      'H','_','H','E','L','I','B','E','B','_','C','_','N','_','O','_','F','_','N','E','N','A','M','G','A','L','S','I','P','_',
      'S','_','C','L','A','R','K','_','C','A','S','C','T','I','V','_','C','R','M','N','F','E','C','O','N','I','C','U','Z','N',
      'G','A','G','E','A','S','S','E','B','R','K','R','R','B','S','R','Y','_','Z','R','N','B','M','O','T','C','R','U','R','H',
      'P','D','A','G','C','D','I','N','S','N','S','B','T','E','I','_','X','E','C','S','B','A','L','A','C','E','P','R','N','D',
      'P','M','S','M','E','U','G','D','T','B','D','Y','H','O','E','R','T','M','Y','B','L','U','H','F','T','A','W','_','R','E',
      'O','S','I','R','P','T','A','U','H','G','T','L','P','B','B','I','P','O','A','T','R','N','F','R','R','A','A','C','T','H',
      'P','A','U','_','N','P','P','U','A','M','C','M','B','K','C','F','E','S','F','M','M','D','N','O','L','R','R','F','D','B',
      'S','G','B','H','H','S','M','T','D','S'};
   memset(zelem, 0, 128*sizeof(Int_t));
   TString sname;
   gGeoManager->Export("flgeom.root");
   if (fname) sname = fname;
   else       sname = "flukaMat.inp";
   ofstream out;
   out.open(sname.Data(), ios::out);
   if (!out.good()) {
      Fatal("CreateFlukaMatFile", "could not open file %s for writing", sname.Data());
      return;
   }
   PrintHeader(out, "MATERIALS AND COMPOUNDS");
   PrintHeader(out, "MATERIALS");   
   TList *matlist = gGeoManager->GetListOfMaterials();
   TIter next(matlist);
   Int_t nmater = matlist->GetSize();
   Int_t nfmater = 0;
   TObjArray *listfluka = new TObjArray(nmater+50);
   TObjArray *listflukanames = new TObjArray(nmater+50);
   TGeoMaterial *mat, *matorig;
   TGeoMixture *mix = 0;
   TString matname;
   TObjString *objstr, *objstrother;
   Int_t i,j,k,idmat;
   Bool_t done;
   Int_t nelem, nidmat;
   Double_t amat,zmat,rhomat;
   Double_t  zel, ael, wel, rho;
   char elname[8] = {' ',' ','_', 'E','L','E','M','\0'}; 
   char digit[3];
   Bool_t found = kFALSE;
   
   if (fDebug) printf("Creating materials and compounds\n");
   for (i=0; i<nmater; i++) {
      mat = (TGeoMaterial*)matlist->At(i);
      if (mat->GetZ()<1E-1) {
         mat->SetIndex(2); // vacuum, built-in inside FLUKA
         continue;
      }     
//      printf("material: %s index=%i: Z=%f A=%f rho=%f\n", mat->GetName(), mat->GetIndex(),mat->GetZ(),mat->GetA(),mat->GetDensity());
      matorig = gGeoManager->FindDuplicateMaterial(mat);
      if (matorig) {
         idmat = matorig->GetIndex();
         mat->SetIndex(idmat);
//         printf(" -> found a duplicate: %s with index %i\n", matorig->GetName(), idmat);
         matorig = 0;
      } else  {
//         printf(" Adding to temp list with index %i\n", nfmater+3);
         listfluka->Add(mat);
         mat->SetIndex(nfmater+3);
         matorig = mat;
         objstr = new TObjString(mat->GetName());
         listflukanames->Add(objstr);
         nfmater++;
         // look if name is existing
         nidmat = 0;
         matname = objstr->GetString();
         ToFlukaString(matname);
         objstr->SetString(matname.Data());
         done = kFALSE;
         while (!done) {
            if (nfmater == 1) break;
            for (j=0; j<nfmater-1; j++) {
               objstrother = (TObjString*)listflukanames->At(j);
               if (objstr->IsEqual(objstrother)) {
                  // we have to change the name
                  if (nidmat>98) {
                     Error("CreateFlukaMatFile", "too many materials having same name");
                     return;
                  }
                  nidmat++;
                  k = matname.Index(" ");
                  if (k<0 || k>6) k=6;
                  if (nidmat>9) {
                     sprintf(digit, "%d", nidmat);
                  } else {
                     digit[0] = '0';
                     sprintf(&digit[1], "%d", nidmat);
                  }
                  matname.Insert(k,digit);
                  matname.Remove(8);
                  objstr->SetString(matname.Data());
                  break;
               }
               if (j == nfmater-2) {
                  done = kTRUE;
                  break;
               }    
            }     
         } 
//         printf(" newmat name: %s\n", matname.Data());                               
      }
      // now we have unique materials with unique names in the lists
         
      if (matorig && matorig->IsMixture()) {
      // create dummy materials for elements
         rho = 0.999;
         mix = (TGeoMixture*)matorig;
         nelem = mix->GetNelements();
//         printf(" material is a MIXTURE with %i elements:\n", nelem);
         for (j=0; j<nelem; j++) {
            found = kFALSE;
            zel = (mix->GetZmixt())[j];
            ael = (mix->GetAmixt())[j];
//            printf("   Zelem[%i] = %g\n",j,zel);
            if ((zel-Int_t(zel))>0.01) {
               TGeoMaterial *mat1;
               for (Int_t imat=0; imat<nfmater; imat++) {
                  mat1 = (TGeoMaterial*)listfluka->At(imat);
                  if (TMath::Abs(mat1->GetZ()-zel)>1E-4) continue;
                  if (TMath::Abs(mat1->GetA()-ael)>1E-4) continue;
                  found = kTRUE;
                  break;
               }      
               if (!found) Warning("CreateFlukaMatFile", "element with Z=%f\n", zel);
            }   
            if (!zelem[Int_t(zel)] && !found) {
               // write fluka element
               memcpy(elname, &elNames[2*Int_t(zel-1)], 2);
               zelem[Int_t(zel)] = 1;
               mat = new TGeoMaterial(elname, ael, zel, rho);
               mat->SetIndex(nfmater+3);
//               printf("  element not in list: new material %s at index=%i, Z=%g, A=%g, dummyrho=%g\n",
//                       elname,nfmater+3,zel,ael,rho);
               listfluka->Add(mat);
               objstr = new TObjString(elname);
               listflukanames->Add(objstr);
               nfmater++;
            }   
         }
      }      
   }
   // now dump materials in the file   
//   printf("DUMPING %i materials\n", nfmater);
   for (i=0; i<nfmater; i++) {
      mat = (TGeoMaterial*)listfluka->At(i);
      out << setw(10) << "MATERIAL  ";
      out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
//      matname = mat->GetName();
      objstr = (TObjString*)listflukanames->At(i);
      matname = objstr->GetString();
      ToFlukaString(matname);
      zmat = mat->GetZ();
      if (zmat-Int_t(zmat)>0.01) {
         if (zmat-Int_t(zmat)>0.5) zmat = Int_t(zmat)+1.;
         else zmat = Int_t(zmat);
      }   
      amat = mat->GetA();
      rhomat = mat->GetDensity();
      // write material card
      if (mat->IsMixture()) {
         out << setw(10) << " ";
         out << setw(10) << " ";
         mix = (TGeoMixture*)mat;
      } else {   
         out << setw(10) << setiosflags(ios::fixed) << setprecision(1) << zmat;
         out << setw(10) << setprecision(3) << amat;
      }
      out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
      out << setw(10) << setiosflags(ios::scientific) << setprecision(3) << rhomat;
      out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
      out << setw(10) << setiosflags(ios::fixed) << setprecision(1) << Double_t(i+3);   
      out << setw(10) << " ";
      out << setw(10) << " ";
      out << setw(8) << matname.Data() << endl;
   } 
   // write mixture header           
   PrintHeader(out, "COMPOUNDS");   
   Int_t counttothree;
   TGeoMaterial *element;
   for (i=0; i<nfmater; i++) {
      mat = (TGeoMaterial*)listfluka->At(i);
      if (!mat->IsMixture()) continue;
      mix = (TGeoMixture*)mat;
      counttothree = 0;
      out << setw(10) << "COMPOUND  ";
      nelem = mix->GetNelements();
      objstr = (TObjString*)listflukanames->At(i);
      matname = objstr->GetString();
//      printf("MIXTURE %s with index %i having %i elements\n", matname.Data(), mat->GetIndex(),nelem);
      for (j=0; j<nelem; j++) {
         // dump mixture cards
//         printf(" #elem %i: Z=%g, A=%g, W=%g\n", j, (mix->GetZmixt())[j], 
//                (mix->GetAmixt())[j],(mix->GetWmixt())[j]); 
         wel = (mix->GetWmixt())[j];
         zel = (mix->GetZmixt())[j];       
         ael = (mix->GetAmixt())[j];
         if (zel-Int_t(zel)>0.01) {
            // loop the temporary list
            element = 0;
            TGeoMaterial *mat1;
            for (Int_t imat=0; imat<i; imat++) {
               mat1 = (TGeoMaterial*)listfluka->At(imat);
               if (TMath::Abs(mat1->GetZ()-zel)>1E-4) continue;
               if (TMath::Abs(mat1->GetA()-ael)>1E-4) continue;
               element = mat1;
               break;
            }      
         } else {
            memcpy(elname, &elNames[2*Int_t(zel-1)], 2);
            element = (TGeoMaterial*)listfluka->FindObject(elname);
         }   
         if (!element) {
            Error("CreateFlukaMatFile", "Element Z=%g %s not found", zel, elname);
            return;
         }
         idmat = element->GetIndex();
//         printf("element %s , index=%i\n", element->GetName(), idmat);
         out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
         out << setw(10) << setiosflags(ios::fixed) << setprecision(6) << -wel;   
         out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
         out << setw(10) << setiosflags(ios::fixed) << setprecision(1) << Double_t(idmat);
         counttothree++;
         if (counttothree == 3) {
            out << matname.Data();
            out << endl;
            if ( (j+1) != nelem) out << setw(10) << "COMPOUND  ";
            counttothree = 0;
         }             
      }
      //Unless we have 3, 6, 9... submaterials we need to put some empty
      //space and the compound name
      if (nelem%3) {
         for (j=0; j<(3-(nelem%3)); j++)
            out << setw(10) << " " << setw(10) << " ";
         out << matname.Data();
         out << endl;
      }   
   }      
   
   // Now print the list of regions (volumes in TGeo)
   Int_t nvols = gGeoManager->GetListOfUVolumes()->GetEntriesFast()-1;
   TGeoVolume *vol;
/*
   PrintHeader(out, "TGEO VOLUMES");
   for (i=1; i<=nvols; i++) {
      vol = gGeoManager->GetVolume(i);
      out.setf(std::ios::left, std::ios::adjustfield);
      out << setw(10) << i;
      out << setw(20) << vol->GetName() << endl;
   }   
*/   
   // Now print the material assignments
   Double_t flagfield;
   PrintHeader(out, "TGEO MATERIAL ASSIGNMENTS");   
   for (i=1; i<=nvols; i++) {
      vol = gGeoManager->GetVolume(i);
      mat = vol->GetMedium()->GetMaterial();
      idmat = mat->GetIndex();
//      flagfield = (vol->GetField())?1.:0.;
      flagfield = 1.;
      out << setw(10) << "ASSIGNMAT ";
      out.setf(static_cast<std::ios::fmtflags>(0),std::ios::floatfield);
      out << setw(10) << setiosflags(ios::fixed) << Double_t(idmat);
      out << setw(10) << setiosflags(ios::fixed) << Double_t(i);
      out << setw(10) << "0.0";
      out << setw(10) << "0.0";
      out << setw(10) << setiosflags(ios::fixed) << flagfield;
      out << setw(10) << "0.0";
      out << endl;
   }
   delete listfluka;
   listflukanames->Delete();
   delete listflukanames;   
   out.close();
   fLastMaterial = nfmater+2;
}

//_____________________________________________________________________________
void TFlukaMCGeometry::PrintHeader(ofstream &out, const char *text) const
{
// Print a FLUKA header.
  out << "*\n" << "*\n" << "*\n";
  out << "*********************  " << text << " *********************\n"
     << "*\n";
  out << "*...+....1....+....2....+....3....+....4....+....5....+....6....+....7..."
     << endl;
  out << "*" << endl;
}

//_____________________________________________________________________________
Int_t TFlukaMCGeometry::RegionId() const
{
// Returns current region id <-> TGeo node id
   if (gGeoManager->IsOutside()) return 0;
   return gGeoManager->GetCurrentNode()->GetUniqueID();
}
//_____________________________________________________________________________
void TFlukaMCGeometry::SetMreg(Int_t mreg)
{
// Update if needed next history;
   Int_t curreg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();
   if (mreg==curreg) return;
   if (mreg==fNextRegion) {
      if (fNextLattice!=999999999) gGeoManager->CdNode(fNextLattice-1);
      return;
   }   
   if (fDebug) printf("ERROR: mreg=%i neither current nor next region\n", mreg);
}

//_____________________________________________________________________________
void TFlukaMCGeometry::SetNextRegion(Int_t mreg, Int_t latt)
{
// Set index/history for next entered region
   fNextRegion = mreg;
   fNextLattice = latt;
}   

//_____________________________________________________________________________
void TFlukaMCGeometry::ToFlukaString(TString &str) const
{
// ToFlukaString converts an string to something usefull in FLUKA:
// * Capital letters
// * Only 8 letters
// * Replace ' ' by '_'
   if (str.Length()<8) {
      str += "        ";
   }   
   str.Remove(8);
   Int_t ilast;
   for (ilast=7; ilast>0; ilast--) if (str(ilast)!=' ') break;
   str.ToUpper();
   for (Int_t pos=0; pos<ilast; pos++)
      if (str(pos)==' ') str.Replace(pos,1,"_",1);
   return;
}   
//______________________________________________________________________________
void TFlukaMCGeometry::Vname(const char *name, char *vname) const
{
  //
  //  convert name to upper case. Make vname at least 4 chars
  //
  Int_t l = strlen(name);
  Int_t i;
  l = l < 4 ? l : 4;
  for (i=0;i<l;i++) vname[i] = toupper(name[i]);
  for (i=l;i<4;i++) vname[i] = ' ';
  vname[4] = 0;      
}


// FLUKA GEOMETRY WRAPPERS - to replace FLUGG wrappers

//_____________________________________________________________________________
Int_t idnrwr(const Int_t & /*nreg*/, const Int_t & /*mlat*/)
{
//   from FLUGG:
// Wrapper for setting DNEAR option on fluka side. Must return 0 
// if user doesn't want Fluka to use DNEAR to compute the 
// step (the same effect is obtained with the GLOBAL (WHAT(3)=-1)
// card in fluka input), returns 1 if user wants Fluka always to 
// use DNEAR (in this case, be sure that GEANT4 DNEAR is unique, 
// coming from all directions!!!)
   if (mcgeom->IsDebugging()) printf("========== Dummy IDNRWR\n");
   return 0;
}

//_____________________________________________________________________________
void g1wr(Double_t &pSx, Double_t &pSy, Double_t &pSz, 
          Double_t *pV,  Int_t &oldReg , const Int_t &oldLttc, Double_t & propStep,
          Int_t & /*nascFlag*/, Double_t &retStep, Int_t &newReg,
	       Double_t &saf, Int_t &newLttc, Int_t &lttcFlag,
          Double_t *sLt, Int_t *jrLt)
{
//   from FLUGG:
// Wrapper for geometry tracking: returns approved step of 
// particle and all variables that fluka G1 computes.

   // Initialize current point/direction
   if (mcgeom->IsDebugging()) {
      printf("========== Inside G1WR\n");
      printf("   point/dir:(%14.9f, %14.9f, %14.9f, %g, %g, %g)\n", pSx,pSy,pSz,pV[0],pV[1],pV[2]);
   }   
   gGeoManager->SetCurrentPoint(pSx, pSy, pSz);
   gGeoManager->SetCurrentDirection(pV);
   if (mcgeom->IsDebugging()) printf("   oldReg=%i  oldLttc=%i  pstep=%f\n",oldReg, oldLttc, propStep);
   if (oldLttc==999999999) printf("WOOPS - wrong old lattice\n");
   if (gGeoManager->IsOutside()) {
      gGeoManager->SetOutside(kFALSE);
      gGeoManager->CdTop();
   }   
   Int_t curLttc = gGeoManager->GetCurrentNodeId()+1;
   Int_t curreg = gGeoManager->GetCurrentVolume()->GetNumber();
   if (mcgeom->IsDebugging()) printf("   curReg=%i  curLttc=%i curPath=%s\n", curreg, curLttc, gGeoManager->GetPath());
   Bool_t regsame = (curreg==oldReg)?kTRUE:kFALSE;
   if (!regsame && mcgeom->IsDebugging()) printf("   REGIONS DOES NOT MATCH\n");
   if (oldLttc != curLttc) {
      if (mcgeom->IsDebugging()) printf("   HISTORIES DOES NOT MATCH\n");
      gGeoManager->CdNode(oldLttc-1);
      curLttc = gGeoManager->GetCurrentNodeId()+1;
      curreg  = gGeoManager->GetCurrentVolume()->GetNumber();
      if (mcgeom->IsDebugging()) printf("   re-initialized point: curReg=%i  curLttc=%i curPath=%s\n", curreg, curLttc, gGeoManager->GetPath());
   }         
   lttcFlag = 0;
   sLt[lttcFlag] = 0.;   
   jrLt[lttcFlag] = curLttc;
   // now 'oldregion' contains the real region, matching or not the old history
   
   // Compute geometry step/safety within physical step limit
//   newReg = oldregion;
   Double_t *point = gGeoManager->GetCurrentPoint();
   Double_t *dir = gGeoManager->GetCurrentDirection();
   Double_t steptot = 0.;
   Double_t snext = 0.;
   Int_t istep = 0;
   Bool_t done = kFALSE;
   Double_t pst;
   Int_t i;
   while (!done) {
      gGeoManager->FindNextBoundary(-propStep);
      snext = gGeoManager->GetStep();
      if (mcgeom->IsDebugging()) printf("   FindNextBoundary(%g) snext=%g\n", propStep, snext);
      if (steptot == 0) {
         saf = gGeoManager->GetSafeDistance();
         if (mcgeom->IsDebugging()) printf("   Safety: %g\n", saf);
      }   
      sLt[lttcFlag] = propStep;
      jrLt[lttcFlag] = gGeoManager->GetCurrentNodeId()+1;     
      lttcFlag++; //1
      sLt[lttcFlag] = 0.;
      jrLt[lttcFlag] = -1;     
      newReg = curreg;
      newLttc = oldLttc;
      if (snext<propStep) {
         // There is a boundary on the way.
         // Make a step=snext+1E-6 to force boundary crossing
         lttcFlag--; // 0
         steptot += snext;
         sLt[lttcFlag] = snext;
         retStep = snext;
//         lttcFlag++;
         // make the step to get into the next region
         for (i=0;i<3;i++) point[i]+=(snext+1E-6)*dir[i];
         gGeoManager->FindNode();
         istep = 0;
         if (mcgeom->IsDebugging()) printf("   boundary: step made %g\n", snext);
         while (gGeoManager->IsSameLocation() && steptot<propStep) {
            if (istep>1E3) {
               printf("Geometry error: could not cross boundary after extra 10 microns\n");
               return;
            }   
            for (i=0;i<3;i++) point[i]+=1E-6*dir[i];
            gGeoManager->FindNode();
            sLt[lttcFlag] += 1E-6;
            retStep = sLt[lttcFlag];
            steptot += 1E-6;
            istep++;
         }            
         if (steptot>propStep) {printf("Error\n");return;}
         // we managed to cross the boundary -> in which region
         newReg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();
         lttcFlag++; //1                
         newLttc = (gGeoManager->IsOutside())?999999999:gGeoManager->GetCurrentNodeId()+1;
         sLt[lttcFlag] = snext; // at 1
         jrLt[lttcFlag] = newLttc;
         sLt[lttcFlag+1] = 0.;
         jrLt[lttcFlag+1] = -1;
         // !!!!!!!!!!

         while (newReg==oldReg && steptot<propStep) {
            if (mcgeom->IsDebugging()) printf("   Entered SAME region... continue\n");
            pst = propStep-steptot;
            gGeoManager->FindNextBoundary(-pst);
            snext = gGeoManager->GetStep();
            steptot += snext;
            if (snext<pst) {
               printf("Found new boundary\n");
               sLt[lttcFlag] = snext;
               retStep = steptot; // ???
               for (i=0;i<3;i++) point[i]+=(snext+1E-6)*dir[i];
               steptot += 1E-6;
               gGeoManager->FindNode();
               if (gGeoManager->IsSameLocation()) {
                  printf("Cannot cross boundary\n");
                  break;
               }
               newReg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();  
               newLttc = (gGeoManager->IsOutside())?999999999:gGeoManager->GetCurrentNodeId()+1;  
               if (mcgeom->IsDebugging()) printf("Found newreg=%i, newLttc=%i, lttFlag is: %i\n", newReg, newLttc, lttcFlag);
               sLt[lttcFlag-1] += snext; // correct step in old region
               sLt[lttcFlag] = propStep-snext;
               jrLt[lttcFlag] = newLttc;
               sLt[lttcFlag+1] = 0.;
               jrLt[lttcFlag+1] = -1;
               if (newReg != oldReg) break; // lttcFlag=1
               lttcFlag++;
            } else {
                if (mcgeom->IsDebugging()) printf("Not crossing next\n");
                lttcFlag--; //0
                retStep=steptot;
                sLt[lttcFlag] = retStep;
                sLt[lttcFlag+1] = 0.;
                jrLt[lttcFlag+1] = -1;
                done = kTRUE;  
            }  
         }
            
         lttcFlag++; //2
         if (mcgeom->IsDebugging()) {
            if (!gGeoManager->IsOutside()) {
               printf("   ENTERED region %i, newLttc=%i in: %s\n", newReg,newLttc,gGeoManager->GetPath());
            } else printf("   EXIT GEOMETRY: BLKHOLE reg=%i\n", newReg);
         }   
      } 
      // no boundary within proposed step
      lttcFlag--;
      done = kTRUE;
   }   
   if (mcgeom->IsDebugging()) printf("=> newReg=%i newLttc=%i lttcFlag=%i\n", newReg, newLttc, lttcFlag);
   mcgeom->SetNextRegion(newReg, newLttc);
   if (mcgeom->IsDebugging()) {
      printf("=> snext=%g safe=%g\n", snext, saf);
      for (Int_t i=0; i<lttcFlag+1; i++) printf("   jrLt[%i]=%i  sLt[%i]=%g\n", i,jrLt[i],i,sLt[i]);
   }   
   if (newLttc!=oldLttc) {
      if (gGeoManager->IsOutside()) {
         gGeoManager->SetOutside(kFALSE);
         gGeoManager->CdTop();
      }   
      gGeoManager->CdNode(oldLttc-1);
   }   
   if (mcgeom->IsDebugging()) printf("<= G1WR (in: %s)\n", gGeoManager->GetPath());
}

//_____________________________________________________________________________
void g1rtwr()
{
   if (mcgeom->IsDebugging()) printf("========== Dummy G1RTWR\n");
} 

//_____________________________________________________________________________
void conhwr(Int_t & /*intHist*/, Int_t * /*incrCount*/)
{
   if (mcgeom->IsDebugging()) printf("========== Dummy CONHWR\n");
}

//_____________________________________________________________________________
void inihwr(Int_t &intHist)
{
   if (mcgeom->IsDebugging()) printf("========== Inside INIHWR -> reinitializing history: %i\n", intHist);
   if (gGeoManager->IsOutside()) gGeoManager->CdTop();
   if (intHist<=0) {
//      printf("=== wrong history number\n");
      return;
   }
   if (intHist==0) gGeoManager->CdTop();
   else gGeoManager->CdNode(intHist-1);
   if (mcgeom->IsDebugging()) {
      printf(" --- current path: %s\n", gGeoManager->GetPath());
      printf("<= INIHWR\n");
   }   
}

//_____________________________________________________________________________
void  jomiwr(const Int_t & /*nge*/, const Int_t & /*lin*/, const Int_t & /*lou*/,
             Int_t &flukaReg)
{
// Geometry initialization wrapper called by FLUKAM. Provides to FLUKA the
// number of regions (volumes in TGeo)
   // build application geometry
   if (mcgeom->IsDebugging()) printf("========== Inside JOMIWR\n");
   flukaReg = gGeoManager->GetListOfUVolumes()->GetEntriesFast();
   if (mcgeom->IsDebugging()) printf("<= JOMIWR: last region=%i\n", flukaReg);
}   

//_____________________________________________________________________________
void lkdbwr(Double_t &pSx, Double_t &pSy, Double_t &pSz,
            Double_t * /*pV*/, const Int_t &oldReg, const Int_t &oldLttc,
            Int_t &newReg, Int_t &flagErr, Int_t &newLttc)             
{
   if (mcgeom->IsDebugging()) {
      printf("========== Inside LKDBWR (%f, %f, %f)\n",pSx, pSy, pSz);
//      printf("   in: pV=(%f, %f, %f)\n", pV[0], pV[1], pV[2]);
      printf("   in: oldReg=%i oldLttc=%i\n", oldReg, oldLttc);
   }   
   TGeoNode *node = gGeoManager->FindNode(pSx, pSy, pSz);
   if (gGeoManager->IsOutside()) {
      newReg = mcgeom->NofVolumes()+1;
//      newLttc = gGeoManager->GetCurrentNodeId();
      newLttc = 999999999;
      if (mcgeom->IsDebugging()) {
         printf("OUTSIDE\n");
         printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
         printf("<= LKMGWR\n");
      }   
      flagErr = newReg;
      return;
   } 
   newReg = node->GetVolume()->GetNumber();
   newLttc = gGeoManager->GetCurrentNodeId()+1; 
   flagErr = newReg;
   if (mcgeom->IsDebugging()) {
      printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
      printf("<= LKDBWR\n");
   }   
}

//_____________________________________________________________________________
void lkfxwr(Double_t &pSx, Double_t &pSy, Double_t &pSz,
            Double_t * /*pV*/, const Int_t &oldReg, const Int_t &oldLttc,
            Int_t &newReg, Int_t &flagErr, Int_t &newLttc)
{
   if (mcgeom->IsDebugging()) {
      printf("========== Inside LKFXWR (%f, %f, %f)\n",pSx, pSy, pSz);
//      printf("   in: pV=(%f, %f, %f)\n", pV[0], pV[1], pV[2]);
      printf("   in: oldReg=%i oldLttc=%i\n", oldReg, oldLttc);
   }   
   TGeoNode *node = gGeoManager->FindNode(pSx, pSy, pSz);
   if (gGeoManager->IsOutside()) {
      newReg = mcgeom->NofVolumes()+1;
//      newLttc = gGeoManager->GetCurrentNodeId();
      newLttc = 999999999;
      if (mcgeom->IsDebugging()) {
         printf("OUTSIDE\n");
         printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
         printf("<= LKMGWR\n");
      }   
      flagErr = newReg;
      return;
   } 
   newReg = node->GetVolume()->GetNumber();
   newLttc = gGeoManager->GetCurrentNodeId()+1; 
   flagErr = newReg;
   if (mcgeom->IsDebugging()) {
      printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
      printf("<= LKFXWR\n");
   }   
}

//_____________________________________________________________________________
void lkmgwr(Double_t &pSx, Double_t &pSy, Double_t &pSz,
            Double_t * /*pV*/, const Int_t &oldReg, const Int_t &oldLttc,
		      Int_t &flagErr, Int_t &newReg, Int_t &newLttc)
{
   if (mcgeom->IsDebugging()) {
      printf("========== Inside LKMGWR (%f, %f, %f)\n",pSx, pSy, pSz);
//      printf("   in: pV=(%f, %f, %f)\n", pV[0], pV[1], pV[2]);
      printf("   in: oldReg=%i oldLttc=%i\n", oldReg, oldLttc);
   }   
   TGeoNode *node = gGeoManager->FindNode(pSx, pSy, pSz);
   if (gGeoManager->IsOutside()) {
      newReg = mcgeom->NofVolumes()+1;
//      newLttc = gGeoManager->GetCurrentNodeId();
      newLttc = 999999999;
      if (mcgeom->IsDebugging()) {
         printf("OUTSIDE\n");
         printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
         printf("<= LKMGWR\n");
      }   
      flagErr = newReg;
      return;
   } 
   newReg = node->GetVolume()->GetNumber();
   newLttc = gGeoManager->GetCurrentNodeId()+1; 
   flagErr = newReg;
   if (mcgeom->IsDebugging()) {
      printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
      printf("<= LKMGWR\n");
   }   
}

//_____________________________________________________________________________
void lkwr(Double_t &pSx, Double_t &pSy, Double_t &pSz,
          Double_t * /*pV*/, const Int_t &oldReg, const Int_t &oldLttc,
	       Int_t &newReg, Int_t &flagErr, Int_t &newLttc)
{
   if (mcgeom->IsDebugging()) {
      printf("========== Inside LKWR (%f, %f, %f)\n",pSx, pSy, pSz);
//      printf("   in: pV=(%f, %f, %f)\n", pV[0], pV[1], pV[2]);
      printf("   in: oldReg=%i oldLttc=%i\n", oldReg, oldLttc);
   }   
   TGeoNode *node = gGeoManager->FindNode(pSx, pSy, pSz);
   if (gGeoManager->IsOutside()) {
      newReg = mcgeom->NofVolumes()+1;
//      newLttc = gGeoManager->GetCurrentNodeId();
      newLttc = 999999999;
      if (mcgeom->IsDebugging()) {
         printf("OUTSIDE\n");
         printf("  out: newReg=%i newLttc=%i\n", newReg, newLttc);
         printf("<= LKMGWR\n");
      }   
      flagErr = newReg;
      return;
   } 
   newReg = node->GetVolume()->GetNumber();
   newLttc = gGeoManager->GetCurrentNodeId()+1; 
   flagErr = newReg;
   if (mcgeom->IsDebugging()) {
      printf("  out: newReg=%i newLttc=%i in %s\n", newReg, newLttc, gGeoManager->GetPath());
      printf("<= LKWR\n");
   }   
}

//_____________________________________________________________________________
void nrmlwr(Double_t &pSx, Double_t &pSy, Double_t &pSz,
            Double_t &pVx, Double_t &pVy, Double_t &pVz,
	         Double_t *norml, const Int_t &oldReg, 
	         const Int_t &newReg, Int_t &flagErr)
{
   if (mcgeom->IsDebugging()) {
      printf("========== Inside NRMLWR (%g, %g, %g, %g, %g, %g)\n", pSx,pSy,pSz,pVx,pVy,pVz);
      printf("   oldReg=%i, newReg=%i\n", oldReg,newReg);
   }   
   Int_t curreg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();
   Int_t curLttc = gGeoManager->GetCurrentNodeId()+1;
   if (mcgeom->IsDebugging()) printf("   curReg=%i, curLttc=%i in: %s\n", curreg, curLttc, gGeoManager->GetPath());
   Bool_t regsame = (curreg==oldReg)?kTRUE:kFALSE;
   gGeoManager->SetCurrentPoint(pSx, pSy, pSz);
   gGeoManager->SetCurrentDirection(pVx,pVy,pVz);
   if (!regsame) {
      if (mcgeom->IsDebugging()) printf("   REGIONS DOEN NOT MATCH\n");
      gGeoManager->FindNode();
      curreg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();
      curLttc = gGeoManager->GetCurrentNodeId()+1;
      if (mcgeom->IsDebugging()) printf("   re-initialized point: curReg=%i  curLttc=%i curPath=%s\n", curreg, curLttc, gGeoManager->GetPath());
   }
   Double_t *dnorm = gGeoManager->FindNormalFast();
   flagErr = 0;
   if (!dnorm) {
      printf("   ERROR: Cannot compute fast normal\n");
      flagErr = 1;
      norml[0] = -pVx;   
      norml[1] = -pVy;   
      norml[2] = -pVz; 
   }
   norml[0] = -dnorm[0];   
   norml[1] = -dnorm[1];   
   norml[2] = -dnorm[2]; 
   if (mcgeom->IsDebugging()) printf("   normal to boundary: (%g, %g, %g)\n", norml[0], norml[1], norml[2]);  
   curreg = (gGeoManager->IsOutside())?(mcgeom->NofVolumes()+1):gGeoManager->GetCurrentVolume()->GetNumber();
   curLttc = gGeoManager->GetCurrentNodeId()+1;
   if (mcgeom->IsDebugging()) {
      printf("   final location: curReg=%i, curLttc=%i in %s\n", curreg,curLttc,gGeoManager->GetPath());
      printf("<= NRMLWR\n");
   }   
}

//_____________________________________________________________________________
void rgrpwr(const Int_t & /*flukaReg*/, const Int_t & /*ptrLttc*/, Int_t & /*g4Reg*/,
            Int_t * /*indMother*/, Int_t * /*repMother*/, Int_t & /*depthFluka*/)
{
   if (mcgeom->IsDebugging()) printf("=> Dummy RGRPWR\n");
}

//_____________________________________________________________________________
Int_t isvhwr(const Int_t &check, const Int_t & intHist)
{
//   from FLUGG:
// Wrapper for saving current navigation history (fCheck=default) 
// and returning its pointer. If fCheck=-1 copy of history pointed 
// by intHist is made in NavHistWithCount object, and its pointer 
// is returned. fCheck=1 and fCheck=2 cases are only in debugging 
// version: an array is created by means of FGeometryInit functions
// (but could be a static int * ptrArray = new int[10000] with 
// file scope as well) that stores a flag for deleted/undeleted 
// histories and at the end of event is checked to verify that 
// all saved history objects have been deleted.

// For TGeo, just return the current node ID. No copy need to be made.

   if (mcgeom->IsDebugging()) printf("=> Inside ISVHWR\n");
   if (check<0) return intHist;
   Int_t histInt = gGeoManager->GetCurrentNodeId()+1;
   if (mcgeom->IsDebugging()) printf("<= ISVHWR: history is: %i in: %s\n", histInt, gGeoManager->GetPath());
   return histInt;
}



   
