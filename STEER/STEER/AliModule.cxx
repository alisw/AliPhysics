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

/* $Id$ */

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

#include <TObjArray.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TSystem.h>
#include <TDirectory.h>
#include <TVirtualMC.h>
#include <TGeoManager.h>
#include <TString.h>

#include "AliLog.h"
#include "AliConfig.h"
#include "AliLoader.h"
#include "AliMagF.h"
#include "AliModule.h"
#include "AliRun.h"
#include "AliTrackReference.h"
#include "AliMC.h"
#include "AliSimulation.h"
#include "AliRawDataHeader.h"
#include "AliDigitizationInput.h"

#include "AliDAQ.h"

ClassImp(AliModule)
 
Float_t AliModule::fgDensityFactor = 1.0;
 
//_______________________________________________________________________
AliModule::AliModule():
  fIdtmed(0),
  fIdmate(0),
  fLoMedium(0),
  fHiMedium(0),
  fActive(0),
  fEnable(1),
  fMaxIterTrackRef(0),
  fCurrentIterTrackRef(0),
  fRunLoader(0),
  fDigInput(0)
{
  //
  // Default constructor for the AliModule class
  //
}
 
//_______________________________________________________________________
AliModule::AliModule(const char* name,const char *title):
  TNamed(name,title),
  fIdtmed(new TArrayI(100)),
  fIdmate(new TArrayI(100)),
  fLoMedium(65536),
  fHiMedium(0),
  fActive(0),
  fEnable(1),
  fMaxIterTrackRef(0),
  fCurrentIterTrackRef(0),
  fRunLoader(0),
  fDigInput(0)
{
  //
  // Normal constructor invoked by all Modules.
  // Create the list for Module specific histograms
  // Add this Module to the global list of Modules in Run.
  //
  // Get the Module numeric ID

  Int_t id = gAlice->GetModuleID(name);
  if (id>=0) {
    // Module already added !
     AliWarning(Form("Module: %s already present at %d",name,id));
     return;
  }
  //
  // Add this Module to the list of Modules

  gAlice->AddModule(this);

  //PH  SetMarkerColor(3);
  //
  // Clear space for tracking media and material indexes

  for(Int_t i=0;i<100;i++) (*fIdmate)[i]=(*fIdtmed)[i]=0;
}
 
//_______________________________________________________________________
AliModule::~AliModule()
{
  //
  // Destructor
  //

  // Remove this Module from the list of Modules
  if (gAlice) {
    TObjArray * modules = gAlice->Modules();
    if (modules) modules->Remove(this);
  }

  // Delete TArray objects
  delete fIdtmed;
  delete fIdmate;

} 

//_______________________________________________________________________
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
  //Build the string uniquename as "DET_materialname"
  TString uniquename = GetName();
  uniquename.Append("_");
  uniquename.Append(name);
  //if geometry loaded from file only fill fIdmate, else create material too
  if(AliSimulation::Instance()->IsGeometryFromFile()){
    TGeoMaterial *mat = gGeoManager->GetMaterial(uniquename.Data());
    kmat = mat->GetUniqueID();
    (*fIdmate)[imat]=kmat;
  }else{
    if (fgDensityFactor != 1.0)
      AliWarning(Form("Material density multiplied by %.2f!", fgDensityFactor));
    gMC->Material(kmat, uniquename.Data(), a, z, dens * fgDensityFactor, radl, absl, buf, nwbuf);
    (*fIdmate)[imat]=kmat;
  }
}
  
//_______________________________________________________________________
void AliModule::AliGetMaterial(Int_t imat, char* name, Float_t &a, 
                               Float_t &z, Float_t &dens, Float_t &radl,
                               Float_t &absl) const
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
  

//_______________________________________________________________________
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
  //Build the string uniquename as "DET_mixturename"
  TString uniquename = GetName();
  uniquename.Append("_");
  uniquename.Append(name);
  //if geometry loaded from file only fill fIdmate, else create mixture too
  if(AliSimulation::Instance()->IsGeometryFromFile()){
    TGeoMaterial *mat = gGeoManager->GetMaterial(uniquename.Data());
    kmat = mat->GetUniqueID();
    (*fIdmate)[imat]=kmat;
  }else{
    if (fgDensityFactor != 1.0)
      AliWarning(Form("Material density multiplied by %.2f!", fgDensityFactor));
    gMC->Mixture(kmat, uniquename.Data(), a, z, dens * fgDensityFactor, nlmat, wmat);
    (*fIdmate)[imat]=kmat;
  }
} 
 
//_______________________________________________________________________
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
  //Build the string uniquename as "DET_mediumname"
  TString uniquename = GetName();
  uniquename.Append("_");
  uniquename.Append(name);
  //if geometry loaded from file only fill fIdtmed, else create medium too
  if(AliSimulation::Instance()->IsGeometryFromFile()){
    TGeoMedium *med = gGeoManager->GetMedium(uniquename.Data());
    kmed = med->GetId();
    (*fIdtmed)[numed]=kmed;
  }else{
    gMC->Medium(kmed, uniquename.Data(), (*fIdmate)[nmat], isvol, ifield,
                fieldm, tmaxfd, stemax, deemax, epsil, stmin, ubuf, nbuf);
    (*fIdtmed)[numed]=kmed;
  }
} 
 
//_______________________________________________________________________
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

//_______________________________________________________________________
Float_t AliModule::ZMin() const
{
  return -500;
}

//_______________________________________________________________________
Float_t AliModule::ZMax() const
{
  return 500;
}

//_______________________________________________________________________
void AliModule::AddAlignableVolumes() const
{
  // 
  if (IsActive())
    AliWarning(Form(" %s still has to implement the AddAlignableVolumes method!",GetName()));
}

//_______________________________________________________________________

AliLoader*  AliModule::MakeLoader(const char* /*topfoldername*/)
{
  return 0x0;
}
 

//_____________________________________________________________________________
AliTrackReference*  AliModule::AddTrackReference(Int_t label, Int_t id){
  //
  // add a trackrefernce to the list
    return (gAlice->GetMCApp()->AddTrackReference(label, id));
}

//_____________________________________________________________________________
TTree* AliModule::TreeTR()
{
  //
  // Return TR tree pointer
  //
  if ( fRunLoader == 0x0)
   {
     AliError("Can not get the run loader");
     return 0x0;
   }

  TTree* tree = fRunLoader->TreeTR();
  return tree;
}


//_____________________________________________________________________________
void AliModule::Digits2Raw()
{
// This is a dummy version that just copies the digits file contents
// to a raw data file.

  AliWarning(Form("Dummy version called for %s", GetName()));

  Int_t nDDLs = AliDAQ::NumberOfDdls(GetName());

  if (!GetLoader()) return;
  fstream digitsFile(GetLoader()->GetDigitsFileName(), ios::in);
  if (!digitsFile) return;

  digitsFile.seekg(0, ios::end);
  UInt_t size = digitsFile.tellg();
  UInt_t ddlSize = 4 * (size / (4*nDDLs));
  Char_t* buffer = new Char_t[ddlSize+1];

  for (Int_t iDDL = 0; iDDL < nDDLs; iDDL++) {
    char fileName[256]="";
    strncpy(fileName,AliDAQ::DdlFileName(GetName(),iDDL),255);
    fstream rawFile(fileName, ios::out);
    if (!rawFile) break;

    AliRawDataHeader header;
    header.fSize = ddlSize + sizeof(header);
    rawFile.write((char*) &header, sizeof(header));

    digitsFile.read(buffer, ddlSize);
    rawFile.write(buffer, ddlSize);
    rawFile.close();
  }

  digitsFile.close();
  delete[] buffer;
}
