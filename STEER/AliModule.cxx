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
Revision 1.10  2000/07/11 18:24:59  fca
Coding convention corrections + few minor bug fixes

Revision 1.9  2000/05/16 08:45:08  fca
Correct dtor, thanks to J.Belikov

Revision 1.8  2000/02/23 16:25:22  fca
AliVMC and AliGeant3 classes introduced
ReadEuclid moved from AliRun to AliModule

Revision 1.7  1999/09/29 09:24:29  fca
Introduction of the Copyright and cvs Log

*/

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
  fIdtmed     = 0;
  fIdmate     = 0;
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
  if (id>=0) {
    // Module already added !
     Warning("Ctor","Module: %s already present at %d\n",name,id);
     return;
  }
  //
  // Add this Module to the list of Modules
  gAlice->Modules()->Add(this);
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
AliModule::AliModule(const AliModule &mod)
{
  //
  // Copy constructor
  //
  mod.Copy(*this);
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
  if(fNodes) {
    fNodes->Clear();
    delete fNodes;
  }
  //
  // Delete TArray objects
  delete fIdtmed;
  delete fIdmate;
}
 
//_____________________________________________________________________________
void AliModule::Copy(AliModule & /* mod */) const
{
  //
  // Copy *this onto mod, not implemented for AliModule
  //
  Fatal("Copy","Not implemented!\n");
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
  gMC->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
  (*fIdmate)[imat]=kmat;
}
  
//_____________________________________________________________________________
void AliModule::AliGetMaterial(Int_t imat, char* name, Float_t &a, 
			      Float_t &z, Float_t &dens, Float_t &radl,
			      Float_t &absl)
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

  Float_t buf[10];
  Int_t nwbuf, kmat;
  kmat=(*fIdmate)[imat];
  gMC->Gfmate(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
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
  gMC->Mixture(kmat, name, a, z, dens, nlmat, wmat);
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
  // numed       the medium number is stored into (*fIdtmed)[numed]
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
  gMC->Medium(kmed,name, (*fIdmate)[nmat], isvol, ifield, fieldm,
			 tmaxfd, stemax, deemax, epsil,	stmin, ubuf, nbuf); 
  (*fIdtmed)[numed]=kmed;
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
  gMC->Matrix(nmat, theta1, phi1, theta2, phi2, theta3, phi3); 
} 

//_____________________________________________________________________________
AliModule& AliModule::operator=(const AliModule &mod)
{
  mod.Copy(*this);
  return (*this);
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
void AliModule::ReadEuclid(const char* filnam, char* topvol)
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
  const Int_t kMaxRot=5000;
  Int_t idrot[kMaxRot],istop[7000];
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
  TArrayI &idtmed = *fIdtmed;
  for(i=1; i<kMaxRot; ++i) idrot[i]=-99;
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
    if( irot<=0 || irot>=kMaxRot ) {
      Error("ReadEuclid","ROTM rotation matrix number %d illegal\n",irot);
      exit(1);
    }
    AliMatrix(idrot[irot],teta1,phi1,teta2,phi2,teta3,phi3);
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
    if( irot<0 || irot>=kMaxRot ) {
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
    if( irot<0 || irot>=kMaxRot ) {
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
void AliModule::ReadEuclidMedia(const char* filnam)
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
  Int_t globField=gAlice->Field()->Integ();
  Float_t globMaxField=gAlice->Field()->Max();
  //  TArrayI &idtmed = *fIdtmed;
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
    AliMaterial(imate,namate,a,z,dens,radl,absl,ubuf,nwbuf);
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
    AliMedium(itmed,natmed,nmat,isvol,globField,globMaxField,tmaxfd,
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
 
