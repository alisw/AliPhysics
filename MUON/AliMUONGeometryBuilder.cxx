/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *      SigmaEffect_thetadegrees                                                                  *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// $Id$
//
// Class AliMUONGeometryBuilder
// ----------------------------
// MUON manager class for geometry construction,
// separated form AliMUONv1
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TClonesArray.h>
#include <TGeoMatrix.h>
#include <TVirtualMC.h>

#include "AliMUONGeometryBuilder.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONChamberGeometry.h"	
#include "AliMUONGeometryEnvelope.h"	
#include "AliMUONGeometryEnvelopeStore.h"
#include "AliMUONGeometryConstituent.h"	
#include "AliMagF.h"
#include "AliRun.h"

ClassImp(AliMUONGeometryBuilder)
 
//______________________________________________________________________________//___________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder() 
  : TObject(),
    fMUON(0),
    fAlign(false),
    fGlobalTransformation(0),
    fGeometryBuilders(0)
{
// Default constructor
} 

//______________________________________________________________________________//___________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(AliMUON* muon)
  : TObject(),
    fMUON(muon),
    fAlign(false),
    fGlobalTransformation(0), 
    fGeometryBuilders(0)
{
// Standars constructor

  // Define the global transformation:
  // Transformation from the old ALICE coordinate system to a new one:
  // x->-x, z->-z 
  TGeoRotation* rotGlobal 
    = new TGeoRotation("rotGlobal", 90., 180., 90., 90., 180., 0.);
  fGlobalTransformation = new TGeoCombiTrans(0., 0., 0., rotGlobal);

  fGeometryBuilders = new TObjArray(AliMUONConstants::NCh());
}

//______________________________________________________________________________
AliMUONGeometryBuilder::AliMUONGeometryBuilder(const AliMUONGeometryBuilder& right) 
  : TObject(right) 
{  
  // copy constructor (not implemented)

  Fatal("AliMUONGeometryBuilder", "Copy constructor not provided.");
}

//______________________________________________________________________________
AliMUONGeometryBuilder::~AliMUONGeometryBuilder()
{
// Destructor

  delete fGlobalTransformation;

  if (fGeometryBuilders){
    fGeometryBuilders->Delete();
    delete fGeometryBuilders;
  }
}

//______________________________________________________________________________
AliMUONGeometryBuilder& 
AliMUONGeometryBuilder::operator=(const AliMUONGeometryBuilder& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//
// private functions
//

//______________________________________________________________________________
void AliMUONGeometryBuilder::PlaceVolume(const TString& name, const TString& mName, 
                            Int_t copyNo, const TGeoHMatrix& matrix, 
			    Int_t npar, Double_t* param, const char* only) const
{
// Place the volume specified by name with the given transformation matrix
// ---

  // Do not apply global transformation 
  // if mother volume == DDIP
  // (as it is applied on this volume)
  TGeoHMatrix transform(matrix);
  if (mName == TString("DDIP")) {
    transform = (*fGlobalTransformation) * transform;
               // To be changed to (*fGlobalTransformation).inverse()
	       // when available in TGeo
	       // To make this correct also for a general case when
	       // (*fGlobalTransformation) * *fGlobalTransformation) != 1
  }	       
     
  // Decompose transformation
  const Double_t* xyz = transform.GetTranslation();
  const Double_t* rm = transform.GetRotationMatrix();
	
  //cout << "Got translation: "
  //     << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	
  //cout << "Got rotation: "
  //     << rm[0] << " " << rm[1] << " " << rm[2] << endl
  //     << rm[3] << " " << rm[4] << " " << rm[5] << endl
  //     << rm[6] << " " << rm[7] << " " << rm[8] << endl;

  // Check for presence of rotation
  // (will be nice to be available in TGeo)
  const Double_t kTolerance = 1e-04;
  Bool_t isRotation = true; 
  if (TMath::Abs(rm[0] - 1.) < kTolerance &&
      TMath::Abs(rm[1] - 0.) < kTolerance &&
      TMath::Abs(rm[2] - 0.) < kTolerance &&
      TMath::Abs(rm[3] - 0.) < kTolerance &&
      TMath::Abs(rm[4] - 1.) < kTolerance &&
      TMath::Abs(rm[5] - 0.) < kTolerance &&
      TMath::Abs(rm[6] - 0.) < kTolerance &&
      TMath::Abs(rm[7] - 0.) < kTolerance &&
      TMath::Abs(rm[8] - 1.) < kTolerance) isRotation = false; 

  Int_t krot = 0;
  if (isRotation) {
    TGeoRotation rot;
    rot.SetMatrix(const_cast<Double_t*>(transform.GetRotationMatrix()));
    Double_t theta1, phi1, theta2, phi2, theta3, phi3;
    rot.GetAngles(theta1, phi1, theta2, phi2, theta3, phi3);
	
    //cout << "angles: " 
    //     << theta1 << " " << phi1 << " "
    //     << theta2 << " " << phi2 << " "
    //     << theta3 << " " << phi3 << endl;
	
    fMUON->AliMatrix(krot, theta1, phi1, theta2, phi2, theta3, phi3);
  }	
	
  // Place the volume in ALIC
  if (npar == 0)
    gMC->Gspos(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only);
  else 
    gMC->Gsposp(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only,
                param, npar);

} 

//
// public functions
//

//______________________________________________________________________________
void AliMUONGeometryBuilder::CreateGeometry()
{
//
// Construct geometry using geometry builders.
//

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create geometry with each builder
    if (builder) {
      if (fAlign) {
        builder->ReadTransformations();
        builder->CreateGeometry();
      }
      else {  
        builder->CreateGeometry();
        builder->SetTransformations();
      } 
    }
  }

  for (Int_t j=0; j<AliMUONConstants::NCh(); j++) {

    AliMUONChamberGeometry* geometry = fMUON->Chamber(j).GetGeometry();

    if (!geometry) continue;
          // Skip chambers with not defined geometry  
	  
    // Loop over envelopes
    const TObjArray* kEnvelopes = geometry->GetEnvelopeStore()->GetEnvelopes();
    for (Int_t k=0; k<kEnvelopes->GetEntriesFast(); k++) {

      // Get envelope
      AliMUONGeometryEnvelope* env = (AliMUONGeometryEnvelope*)kEnvelopes->At(k);
      const TGeoCombiTrans* kEnvTrans = env->GetTransformation();
      const char* only = "ONLY";
      if (env->IsMANY()) only = "MANY";

      if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() == 0 ) {
        // virtual envelope + nof constituents = 0 
        //         => not allowed;
        //            empty virtual envelope has no sense 
        Fatal("CreateGeometry", "Virtual envelope must have constituents.");
        return;
      }

      if (!env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
        // non virtual envelope + nof constituents > 0 
        //        => not allowed;
        //           use VMC to place constituents
        Fatal("CreateGeometry", "Non virtual envelope cannot have constituents.");
        return;
      }

      if (!env->IsVirtual() && env->GetConstituents()->GetEntriesFast() == 0 ) {
        // non virtual envelope + nof constituents = 0 
        //        => place envelope in ALICE by composed transformation:
        //           Tglobal * Tch * Tenv

        // Compound chamber transformation with the envelope one
        TGeoHMatrix total 
	  = (*fGlobalTransformation) * 
	    (*geometry->GetTransformation()) * 
	    (*kEnvTrans);
        PlaceVolume(env->GetName(), geometry->GetMotherVolume(),
	            env->GetCopyNo(), total, 0, 0, only);
      }

      if (env->IsVirtual() && env->GetConstituents()->GetEntriesFast() > 0 ) {
        // virtual envelope + nof constituents > 0 
        //         => do not place envelope and place constituents
        //            in ALICE by composed transformation:
        //            Tglobal * Tch * Tenv * Tconst   

        for  (Int_t l=0; l<env->GetConstituents()->GetEntriesFast(); l++) {
          AliMUONGeometryConstituent* constituent
            = (AliMUONGeometryConstituent*)env->GetConstituents()->At(l);

          // Compound chamber transformation with the envelope one + the constituent one
          TGeoHMatrix total 
	    = (*fGlobalTransformation) *
	      (*geometry->GetTransformation()) * 
	      (*kEnvTrans) * 
	      (*constituent->GetTransformation());

          PlaceVolume(constituent->GetName(), geometry->GetMotherVolume(),
	              constituent->GetCopyNo(), total,
                      constituent->GetNpar(), constituent->GetParam(), only);
        }
      }
    } 
  }
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::CreateMaterials()
{
  // Definition of common materials
  // --

  //
  //     Ar-CO2 gas (80%+20%)
  Float_t ag1[3]   = { 39.95,12.01,16. };
  Float_t zg1[3]   = { 18.,6.,8. };
  Float_t wg1[3]   = { .8,.0667,.13333 };
  Float_t dg1      = .001821;
  //
  //     Ar-buthane-freon gas -- trigger chambers 
  Float_t atr1[4]  = { 39.95,12.01,1.01,19. };
  Float_t ztr1[4]  = { 18.,6.,1.,9. };
  Float_t wtr1[4]  = { .56,.1262857,.2857143,.028 };
  Float_t dtr1     = .002599;
  //
  //     Ar-CO2 gas 
  Float_t agas[3]  = { 39.95,12.01,16. };
  Float_t zgas[3]  = { 18.,6.,8. };
  Float_t wgas[3]  = { .74,.086684,.173316 };
  Float_t dgas     = .0018327;
  //
  //     Ar-Isobutane gas (80%+20%) -- tracking 
  Float_t ag[3]    = { 39.95,12.01,1.01 };
  Float_t zg[3]    = { 18.,6.,1. };
  Float_t wg[3]    = { .8,.057,.143 };
  Float_t dg       = .0019596;
  //
  //     Ar-Isobutane-Forane-SF6 gas (49%+7%+40%+4%) -- trigger 
  Float_t atrig[5] = { 39.95,12.01,1.01,19.,32.066 };
  Float_t ztrig[5] = { 18.,6.,1.,9.,16. };
  Float_t wtrig[5] = { .49,1.08,1.5,1.84,0.04 };
  Float_t dtrig    = .0031463;
  //
  //     bakelite: C6 H6 O
  Float_t abak[3] = {12.01 , 1.01 , 16.};
  Float_t zbak[3] = {6.     , 1.   , 8.};
  Float_t wbak[3] = {6.     , 6.   , 1.}; 
  Float_t dbak = 1.4;

  Int_t iSXFLD   = gAlice->Field()->Integ();
  Float_t sXMGMX = gAlice->Field()->Max();
  //
  // --- Define the various materials for GEANT --- 
  fMUON->AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  fMUON->AliMaterial(10, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
  // Air
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  fMUON->AliMixture(15, "AIR$      ", aAir,  zAir, dAir,4, wAir);
  //    fMUON->AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
  fMUON->AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
  fMUON->AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
  fMUON->AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
  fMUON->AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
  fMUON->AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
  fMUON->AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);

  // materials for slat: 
  //     Sensitive area: gas (already defined) 
  //     PCB: copper 
  //     insulating material: vetronite -> replacing by G10 Ch. Finck
  //     spacer: noryl Ch. Finck
  //     panel sandwich: carbon, nomex, carbon replacing rohacell by nomex Ch. Finck

  // G10: SiO2(60%) + C8H14O4(40%)
  Float_t aglass[5] = {12.01, 28.09, 16., 1.01,  16.};
  Float_t zglass[5] = { 6.,   14.,    8., 1.,    8.};
  Float_t wglass[5] = { 0.22, 0.28, 0.32, 0.03,  0.15};
  Float_t dglass    = 1.7;

  // rohacell: C9 H13 N1 O2
  Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
  Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
  Float_t drohac    = 0.03;

  // Nomex: C22 H10 N2 O5
  Float_t aNomex[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zNomex[4] = { 6.,    1.,    7.,     8.};
  Float_t wNomex[4] = { 22.,   10.,   2.,     5.};
  Float_t dNomex    = 0.024; //honey comb
  Float_t dNomex2   = 1.43;  //bulk material


  // Noryl: C8 H8 O polyphenylene oxyde (di-methyl not sure)
  Float_t aNoryl[3] = {12.01,  1.01, 16.};
  Float_t zNoryl[3] = { 6.,    1.,    8.};
  Float_t wNoryl[3] = { 8.,    8.,    1.};
  Float_t dNoryl    = 1.06;

  fMUON->AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,   1.4, 0.);
  fMUON->AliMixture( 32, "G10$",      aglass, zglass, dglass, -5, wglass);
  fMUON->AliMaterial(33, "Carbon$",   12.01,     6.,  2.265,  18.8, 49.9);
  fMUON->AliMixture( 34, "Rohacell$", arohac, zrohac, drohac, -4, wrohac); 
  fMUON->AliMixture( 35, "Nomex$",    aNomex, zNomex, dNomex, -4, wNomex); 
  fMUON->AliMixture( 36, "Noryl$",    aNoryl, zNoryl, dNoryl, -3, wNoryl); 
  fMUON->AliMixture( 37, "Nomex_bulk$",aNomex, zNomex, dNomex2, -4, wNomex); 

  Float_t  epsil  = .001; // Tracking precision, 
  Float_t  stemax = -1.;  // Maximum displacement for multiple scat 
  Float_t  tmaxfd = -20.; // Maximum angle due to field deflection 
  Float_t  deemax = -.3;  // Maximum fractional energy loss, DLS 
  Float_t  stmin  = -.8;
  Float_t  maxDestepAlu = fMUON->GetMaxDestepAlu();
  Float_t  maxDestepGas = fMUON->GetMaxDestepGas();
  Float_t  maxStepAlu = fMUON->GetMaxStepAlu();
  Float_t  maxStepGas = fMUON->GetMaxStepGas();

  //
  //    Air 
  fMUON->AliMedium(1, "AIR_CH_US         ", 15, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
 
  //
  //    Aluminum 
  fMUON->AliMedium(4, "ALU_CH_US          ", 9, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(5, "ALU_CH_US          ", 10, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  //
  //    Ar-isoC4H10 gas 
  fMUON->AliMedium(6, "AR_CH_US          ", 20, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas, 
		   maxDestepGas, epsil, stmin);
  //
  //    Ar-Isobuthane-Forane-SF6 gas 
  fMUON->AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

  fMUON->AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, iSXFLD, sXMGMX, tmaxfd, maxStepAlu, 
		   maxDestepAlu, epsil, stmin);
  //
  // slat medium
  fMUON->AliMedium(9, "ARG_CO2   ", 22, 1, iSXFLD, sXMGMX, tmaxfd, maxStepGas, 
		   maxDestepAlu, epsil, stmin);
  //
  // tracking media for slats: check the parameters!! 
  fMUON->AliMedium(11, "PCB_COPPER        ", 31, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(12, "G10               ", 32, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(13, "CARBON            ", 33, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(14, "Rohacell          ", 34, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(15, "Nomex             ", 35, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(16, "Noryl             ", 36, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  fMUON->AliMedium(17, "Nomex bulk        ", 37, 0, iSXFLD, sXMGMX, tmaxfd, 
		   maxStepAlu, maxDestepAlu, epsil, stmin);
  //.Materials specific to stations
  // created via builders
  
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Create materials with each builder
    if (builder) builder->CreateMaterials();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::InitGeometry()
{
 // Initialize geometry
 // ---

  //
  // Set the chamber (sensitive region) GEANT identifier
  //
  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Set sesitive volumes with each builder
    builder->SetSensitiveVolumes();

    // Read sensitive volume map from a file
    builder->ReadSVMap();
    if (!fAlign)  builder->FillTransformations();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteTransformations()
{
 // Writes transformations into files per builder
 // ---

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Write transformations
    builder->WriteTransformations();
  }
}

//______________________________________________________________________________
void AliMUONGeometryBuilder::WriteSVMaps(Bool_t rebuild)
{
 // Writes sensitive volume maps into files per builder
 // ---

  for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Write transformations      
    builder->WriteSVMap(rebuild);
  }
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::AddBuilder(AliMUONVGeometryBuilder* geomBuilder)
{
// Adds the geometry builder to the list
// ---

  fGeometryBuilders->Add(geomBuilder);
}

//_____________________________________________________________________________
void AliMUONGeometryBuilder::SetAlign(Bool_t align)
{ 
// Sets the option for alignement
// ---

  fAlign = align; 

  for (Int_t j=0; j<AliMUONConstants::NCh(); j++) {

    AliMUONChamberGeometry* geometry = fMUON->Chamber(j).GetGeometry();

    if (!geometry) continue;
          // Skip chambers with not defined geometry  
	  
    geometry->SetAlign(align);	  
  }	  
}


