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

/* $Id$ */

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:MUON version 1    //
/////////////////////////////////////////////////////////

#include <TRandom.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TRandom.h> 
#include <TGeoMatrix.h>
#include <TVirtualMC.h>

#include "AliMUONv1.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONFactory.h"
#include "AliMUONHit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONVGeometryBuilder.h"	
#include "AliMUONChamberGeometry.h"	
#include "AliMUONGeometryEnvelope.h"	
#include "AliMUONGeometryConstituent.h"	
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliMUONv1)
 
//___________________________________________
AliMUONv1::AliMUONv1() 
  : AliMUON(),
    fTrackMomentum(), fTrackPosition(),fGlobalTransformation(0) 
{
// Constructor
    fChambers   = 0;
    fStepManagerVersionOld  = kFALSE;
    fAngleEffect = kTRUE;
    fStepMaxInActiveGas     = 0.6;
    fStepSum    =  0x0;
    fDestepSum  =  0x0;
    fElossRatio =  0x0;
    fAngleEffect10   = 0x0;
    fAngleEffectNorma= 0x0;
} 
//___________________________________________
AliMUONv1::AliMUONv1(const char *name, const char *title)
  : AliMUON(name,title), fTrackMomentum(), fTrackPosition()
{
// Constructor
    // By default include all stations
    AliMUONFactory factory;
    factory.Build(this, title);

    fStepManagerVersionOld = kFALSE;
    fAngleEffect = kTRUE;
    fStepMaxInActiveGas = 0.6;

    fStepSum   = new Float_t [AliMUONConstants::NCh()];
    fDestepSum = new Float_t [AliMUONConstants::NCh()];
    for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {
      fStepSum[i] =0.0;
      fDestepSum[i]=0.0;
    }
    // Ratio of particle mean eloss with respect MIP's Khalil Boudjemline, sep 2003, PhD.Thesis and Particle Data Book
    fElossRatio = new TF1("ElossRatio","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",0.5,5.); 
    fElossRatio->SetParameter(0,1.02138);
    fElossRatio->SetParameter(1,-9.54149e-02);
    fElossRatio->SetParameter(2,+7.83433e-02); 
    fElossRatio->SetParameter(3,-9.98208e-03);
    fElossRatio->SetParameter(4,+3.83279e-04);

    // Angle effect in tracking chambers at theta =10 degres as a function of ElossRatio (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis) (in micrometers)
    fAngleEffect10 = new TF1("AngleEffect10","[0]+[1]*x+[2]*x*x",0.5,3.0);
    fAngleEffect10->SetParameter(0, 1.90691e+02);
    fAngleEffect10->SetParameter(1,-6.62258e+01);
    fAngleEffect10->SetParameter(2,+1.28247e+01);
    // Angle effect: Normalisation form theta=10 degres to theta between 0 and 10 (Khalil BOUDJEMLINE sep 2003 Ph.D Thesis)  
    // Angle with respect to the wires assuming that chambers are perpendicular to the z axis.
    fAngleEffectNorma = new TF1("AngleEffectNorma","[0]+[1]*x+[2]*x*x+[3]*x*x*x",0.0,10.0);
    fAngleEffectNorma->SetParameter(0,4.148);
    fAngleEffectNorma->SetParameter(1,-6.809e-01);
    fAngleEffectNorma->SetParameter(2,5.151e-02);
    fAngleEffectNorma->SetParameter(3,-1.490e-03);

    // Define the global transformation:
    // Transformation from the old ALICE coordinate system to a new one:
    // x->-x, z->-z 
    TGeoRotation* rotGlobal 
      = new TGeoRotation("rotGlobal", 90., 180., 90., 90., 180., 0.);
    fGlobalTransformation = new TGeoCombiTrans(0., 0., 0., rotGlobal);
}

//_____________________________________________________________________________
AliMUONv1::AliMUONv1(const AliMUONv1& right) 
  : AliMUON(right) 
{  
  // copy constructor (not implemented)

  Fatal("AliMUONv1", "Copy constructor not provided.");
}

//___________________________________________
AliMUONv1::~AliMUONv1()
{
// Destructor

  delete fGlobalTransformation;
}

//_____________________________________________________________________________
AliMUONv1& AliMUONv1::operator=(const AliMUONv1& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  Fatal("operator =", "Assignement operator not provided.");
    
  return *this;  
}    

//__________________________________________________
void AliMUONv1::CreateGeometry()
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
      builder->CreateGeometry();
      builder->SetTransformations();
    }
  }

  for (Int_t j=0; j<AliMUONConstants::NCh(); j++) {

    AliMUONChamberGeometry* geometry = Chamber(j).GetGeometry();

    if (!geometry) continue;
          // Skip chambers with not defined geometry  
	  
    // Loop over envelopes
    const TObjArray* kEnvelopes = geometry->GetEnvelopes();
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
//__________________________________________________________________
Int_t  AliMUONv1::GetChamberId(Int_t volId) const
{
// Check if the volume with specified  volId is a sensitive volume (gas) 
// of some chamber and returns the chamber number;
// if not sensitive volume - return 0.
// ---

/*
  for (Int_t i = 1; i <= AliMUONConstants::NCh(); i++)
    if (volId==((AliMUONChamber*)(*fChambers)[i-1])->GetGid()) return i;
*/
  for (Int_t i = 1; i <= AliMUONConstants::NCh(); i++)
    if ( ((AliMUONChamber*)(*fChambers)[i-1])->IsSensId(volId) ) return i;

  return 0;
}
//________________________________________________________________
void AliMUONv1::CreateMaterials()
{

  // *** DEFINITION OF AVAILABLE MUON MATERIALS *** 
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
    //     bakelite 

    Float_t abak[3] = {12.01 , 1.01 , 16.};
    Float_t zbak[3] = {6.     , 1.   , 8.};
    Float_t wbak[3] = {6.     , 6.   , 1.}; 
    Float_t dbak = 1.4;

    Float_t epsil, stmin, deemax, tmaxfd, stemax;

    Int_t iSXFLD   = gAlice->Field()->Integ();
    Float_t sXMGMX = gAlice->Field()->Max();
    //
    // --- Define the various materials for GEANT --- 
    AliMaterial(9, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial(10, "ALUMINIUM$", 26.98, 13., 2.7, 8.9, 37.2);
    AliMaterial(15, "AIR$      ", 14.61, 7.3, .001205, 30423.24, 67500);
    AliMixture(19, "Bakelite$", abak, zbak, dbak, -3, wbak);
    AliMixture(20, "ArC4H10 GAS$", ag, zg, dg, 3, wg);
    AliMixture(21, "TRIG GAS$", atrig, ztrig, dtrig, -5, wtrig);
    AliMixture(22, "ArCO2 80%$", ag1, zg1, dg1, 3, wg1);
    AliMixture(23, "Ar-freon $", atr1, ztr1, dtr1, 4, wtr1);
    AliMixture(24, "ArCO2 GAS$", agas, zgas, dgas, 3, wgas);
    // materials for slat: 
    //     Sensitive area: gas (already defined) 
    //     PCB: copper 
    //     insulating material and frame: vetronite
    //     walls: carbon, rohacell, carbon 
  Float_t aglass[5]={12.01, 28.09, 16.,   10.8,  23.};
  Float_t zglass[5]={ 6.,   14.,    8.,    5.,   11.};
  Float_t wglass[5]={ 0.5,  0.105, 0.355, 0.03,  0.01};
  Float_t dglass=1.74;

  // rohacell: C9 H13 N1 O2
  Float_t arohac[4] = {12.01,  1.01, 14.010, 16.};
  Float_t zrohac[4] = { 6.,    1.,    7.,     8.};
  Float_t wrohac[4] = { 9.,   13.,    1.,     2.};
  Float_t drohac    = 0.03;

  AliMaterial(31, "COPPER$",   63.54,    29.,   8.96,  1.4, 0.);
  AliMixture(32, "Vetronite$",aglass, zglass, dglass,    5, wglass);
  AliMaterial(33, "Carbon$",   12.01,     6.,  2.265, 18.8, 49.9);
  AliMixture(34, "Rohacell$", arohac, zrohac, drohac,   -4, wrohac); 


    epsil  = .001; // Tracking precision, 
    stemax = -1.;  // Maximum displacement for multiple scat 
    tmaxfd = -20.; // Maximum angle due to field deflection 
    deemax = -.3;  // Maximum fractional energy loss, DLS 
    stmin  = -.8;
    //
    //    Air 
    AliMedium(1, "AIR_CH_US         ", 15, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);
    //
    //    Aluminum 

    AliMedium(4, "ALU_CH_US          ", 9, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
    AliMedium(5, "ALU_CH_US          ", 10, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);
    //
    //    Ar-isoC4H10 gas 

    AliMedium(6, "AR_CH_US          ", 20, 1, iSXFLD, sXMGMX, tmaxfd, fMaxStepGas, 
	    fMaxDestepGas, epsil, stmin);
//
    //    Ar-Isobuthane-Forane-SF6 gas 

    AliMedium(7, "GAS_CH_TRIGGER    ", 21, 1, iSXFLD, sXMGMX, tmaxfd, stemax, deemax, epsil, stmin);

    AliMedium(8, "BAKE_CH_TRIGGER   ", 19, 0, iSXFLD, sXMGMX, tmaxfd, fMaxStepAlu, 
	    fMaxDestepAlu, epsil, stmin);

    AliMedium(9, "ARG_CO2   ", 22, 1, iSXFLD, sXMGMX, tmaxfd, fMaxStepGas, 
	    fMaxDestepAlu, epsil, stmin);
    // tracking media for slats: check the parameters!! 
    AliMedium(11, "PCB_COPPER        ", 31, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(12, "VETRONITE         ", 32, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(13, "CARBON            ", 33, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);
    AliMedium(14, "Rohacell          ", 34, 0, iSXFLD, sXMGMX, tmaxfd, 
	      fMaxStepAlu, fMaxDestepAlu, epsil, stmin);



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
void AliMUONv1::PlaceVolume(const TString& name, const TString& mName, 
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
	
    AliMatrix(krot, theta1, phi1, theta2, phi2, theta3, phi3);
  }	
	
  // Place the volume in ALIC
  if (npar == 0)
    gMC->Gspos(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only);
  else 
    gMC->Gsposp(name, copyNo, mName, xyz[0], xyz[1], xyz[2] , krot, only,
                param, npar);

} 

//___________________________________________
void AliMUONv1::Init()
{
   // 
   // Initialize Tracking Chambers
   //

   if(fDebug) printf("\n%s: Start Init for version 1 - CPC chamber type\n\n",ClassName());
   Int_t i;
   for (i=0; i<AliMUONConstants::NCh(); i++) {
       ( (AliMUONChamber*) (*fChambers)[i])->Init();
   }
   
   //
   // Set the chamber (sensitive region) GEANT identifier
   //
   for (Int_t i=0; i<fGeometryBuilders->GetEntriesFast(); i++) {

    // Get the builder
    AliMUONVGeometryBuilder* builder
      = (AliMUONVGeometryBuilder*)fGeometryBuilders->At(i);

    // Set sesitive volumes with each builder
    if (builder) builder->SetSensitiveVolumes();
  }

/*
   //
   // Set the chamber (sensitive region) GEANT identifier
   ((AliMUONChamber*)(*fChambers)[0])->SetGid(gMC->VolId("S01G"));
   ((AliMUONChamber*)(*fChambers)[1])->SetGid(gMC->VolId("S02G"));

   ((AliMUONChamber*)(*fChambers)[2])->SetGid(gMC->VolId("S03G"));
   ((AliMUONChamber*)(*fChambers)[3])->SetGid(gMC->VolId("S04G"));

   ((AliMUONChamber*)(*fChambers)[4])->SetGid(gMC->VolId("S05G"));
   ((AliMUONChamber*)(*fChambers)[5])->SetGid(gMC->VolId("S06G"));

   ((AliMUONChamber*)(*fChambers)[6])->SetGid(gMC->VolId("S07G"));
   ((AliMUONChamber*)(*fChambers)[7])->SetGid(gMC->VolId("S08G"));

   ((AliMUONChamber*)(*fChambers)[8])->SetGid(gMC->VolId("S09G"));
   ((AliMUONChamber*)(*fChambers)[9])->SetGid(gMC->VolId("S10G"));

   ((AliMUONChamber*)(*fChambers)[10])->SetGid(gMC->VolId("SG1A"));
   ((AliMUONChamber*)(*fChambers)[11])->SetGid(gMC->VolId("SG2A"));
   ((AliMUONChamber*)(*fChambers)[12])->SetGid(gMC->VolId("SG3A"));
   ((AliMUONChamber*)(*fChambers)[13])->SetGid(gMC->VolId("SG4A"));
*/
   if(fDebug) printf("\n%s: Finished Init for version 1 - CPC chamber type\n",ClassName());

   //cp 
   if(fDebug) printf("\n%s: Start Init for Trigger Circuits\n",ClassName());
   for (i=0; i<AliMUONConstants::NTriggerCircuit(); i++) {
     ( (AliMUONTriggerCircuit*) (*fTriggerCircuits)[i])->Init(i);
   }
   if(fDebug) printf("%s: Finished Init for Trigger Circuits\n",ClassName());
   //cp
}

//_______________________________________________________________________________
void AliMUONv1::StepManager()
{
  // Stepmanager for the chambers

 if (fStepManagerVersionOld) {
    StepManagerOld();
    return;
  }

  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  // Only charged tracks
  
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  static Int_t   idvol=-1;
  Int_t   iChamber=0;
  Int_t   id=0;
  Int_t   copy;
  const  Float_t kBig = 1.e10;


  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  id=gMC->CurrentVolID(copy);
  iChamber = GetChamberId(id);
  idvol = iChamber -1;

  if (idvol == -1) return;

  // Filling TrackRefs file for MUON. Our Track references are the active volume of the chambers
  if ( (gMC->IsTrackEntering() || gMC->IsTrackExiting() ) )     
    AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
  
   if( gMC->IsTrackEntering() ) {
     Float_t theta = fTrackMomentum.Theta();
     if ((TMath::Pi()-theta)*kRaddeg>=15.) gMC->SetMaxStep(fStepMaxInActiveGas); // We use Pi-theta because z is negative
  }

//  if (GetDebug()) {
//     Float_t z = ( (AliMUONChamber*)(*fChambers)[idvol])->Z() ;
//      Info("StepManager Step","Active volume found %d chamber %d Z chamber is %f ",idvol,iChamber, z);
//   }  
  // Particule id and mass, 
  Int_t     ipart = gMC->TrackPid();
  Float_t   mass  = gMC->TrackMass();

  fDestepSum[idvol]+=gMC->Edep();
  // Get current particle id (ipart), track position (pos)  and momentum (mom)
  if ( fStepSum[idvol]==0.0 )  gMC->TrackMomentum(fTrackMomentum);
  fStepSum[idvol]+=gMC->TrackStep();
  
//   if (GetDebug()) {
//     Info("StepManager Step","iChamber %d, Particle %d, theta %f phi %f mass %f StepSum %f eloss %g",
//       iChamber,ipart, fTrackMomentum.Theta()*kRaddeg, fTrackMomentum.Phi()*kRaddeg, mass, fStepSum[idvol], gMC->Edep());
//     Info("StepManager Step","Track Momentum %f %f %f", fTrackMomentum.X(), fTrackMomentum.Y(), fTrackMomentum.Z()) ;
//     gMC->TrackPosition(fTrackPosition);
//     Info("StepManager Step","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
//   }

  // Track left chamber or StepSum larger than fStepMaxInActiveGas
  if ( gMC->IsTrackExiting() || 
       gMC->IsTrackStop() || 
       gMC->IsTrackDisappeared()||
       (fStepSum[idvol]>fStepMaxInActiveGas) ) {
    
    if   ( gMC->IsTrackExiting() || 
           gMC->IsTrackStop() || 
           gMC->IsTrackDisappeared() ) gMC->SetMaxStep(kBig);

    gMC->TrackPosition(fTrackPosition);
    Float_t theta = fTrackMomentum.Theta();
    Float_t phi   = fTrackMomentum.Phi();
    
    TLorentzVector backToWire( fStepSum[idvol]/2.*sin(theta)*cos(phi),
                               fStepSum[idvol]/2.*sin(theta)*sin(phi),
                               fStepSum[idvol]/2.*cos(theta),0.0       );
    //     if (GetDebug()) 
    //       Info("StepManager Exit","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
    //     if (GetDebug()) 
    //        Info("StepManager Exit ","Track backToWire %f %f %f",backToWire.X(),backToWire.Y(),backToWire.Z()) ;
    fTrackPosition-=backToWire;
    
    //-------------- Angle effect 
    // Ratio between energy loss of particle and Mip as a function of BetaGamma of particle (Energy/Mass)
    
    Float_t betaxGamma    = fTrackMomentum.P()/mass;//  pc/mc2
    Float_t sigmaEffect10degrees;
    Float_t sigmaEffectThetadegrees;
    Float_t eLossParticleELossMip;
    Float_t yAngleEffect=0.;
    Float_t thetawires      =  TMath::Abs( TMath::ASin( TMath::Sin(TMath::Pi()-theta) * TMath::Sin(phi) ) );// We use Pi-theta because z is negative


    if (fAngleEffect){
    if ( (betaxGamma >3.2)   &&  (thetawires*kRaddeg<=15.) ) {
      betaxGamma=TMath::Log(betaxGamma);
      eLossParticleELossMip = fElossRatio->Eval(betaxGamma);
      // 10 degrees is a reference for a model (arbitrary)
      sigmaEffect10degrees=fAngleEffect10->Eval(eLossParticleELossMip);// in micrometers
      // Angle with respect to the wires assuming that chambers are perpendicular to the z axis.
      sigmaEffectThetadegrees =  sigmaEffect10degrees/fAngleEffectNorma->Eval(thetawires*kRaddeg);  // For 5mm gap  
      if ( (iChamber==1)  ||  (iChamber==2) )  
        sigmaEffectThetadegrees/=(1.09833e+00+1.70000e-02*(thetawires*kRaddeg)); // The gap is different (4mm)
      yAngleEffect=1.e-04*gRandom->Gaus(0,sigmaEffectThetadegrees); // Error due to the angle effect in cm
    }
    }
    
    // One hit per chamber
    GetMUONData()->AddHit(fIshunt, gAlice->GetMCApp()->GetCurrentTrackNumber(), iChamber, ipart, 
                          fTrackPosition.X(), fTrackPosition.Y()+yAngleEffect, fTrackPosition.Z(), 0.0, 
                          fTrackMomentum.P(),theta, phi, fStepSum[idvol], fDestepSum[idvol],
                          fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z());
//     if (GetDebug()){
//       Info("StepManager Exit","Particle exiting from chamber %d",iChamber);
//       Info("StepManager Exit","StepSum %f eloss geant %g ",fStepSum[idvol],fDestepSum[idvol]);
//       Info("StepManager Exit","Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;
//     }
    fStepSum[idvol]  =0; // Reset for the next event
    fDestepSum[idvol]=0; // Reset for the next event
  }
}

//__________________________________________
void AliMUONv1::StepManagerOld()
{
  // Old Stepmanager for the chambers
  Int_t          copy, id;
  static Int_t   idvol =-1;
  static Int_t   vol[2];
  Int_t          ipart;
  TLorentzVector pos;
  TLorentzVector mom;
  Float_t        theta,phi;
  Float_t        destep, step;
  
  static Float_t sstep;
  static Float_t eloss, eloss2, xhit, yhit, zhit, tof, tlength;
  const  Float_t kBig = 1.e10;
  static Float_t hits[15];

  TClonesArray &lhits = *fHits;

  //
  //
  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  id=gMC->CurrentVolID(copy);
  vol[0] = GetChamberId(id);
  idvol = vol[0] -1;

  if (idvol == -1) return;

  //
  // Get current particle id (ipart), track position (pos)  and momentum (mom) 
  gMC->TrackPosition(pos);
  gMC->TrackMomentum(mom);

  ipart  = gMC->TrackPid();

  //
  // momentum loss and steplength in last step
  destep = gMC->Edep();
  step   = gMC->TrackStep();
  // cout<<"------------"<<step<<endl;
  //
  // record hits when track enters ...
  if( gMC->IsTrackEntering()) {

      gMC->SetMaxStep(fMaxStepGas);
      Double_t tc = mom[0]*mom[0]+mom[1]*mom[1];
      Double_t rt = TMath::Sqrt(tc);
      Double_t pmom = TMath::Sqrt(tc+mom[2]*mom[2]);
      Double_t tx = mom[0]/pmom;
      Double_t ty = mom[1]/pmom;
      Double_t tz = mom[2]/pmom;
      Double_t s  = ((AliMUONChamber*)(*fChambers)[idvol])
          ->ResponseModel()
          ->Pitch()/tz;
      theta   = Float_t(TMath::ATan2(rt,Double_t(mom[2])))*kRaddeg;
      phi     = Float_t(TMath::ATan2(Double_t(mom[1]),Double_t(mom[0])))*kRaddeg;
      hits[0] = Float_t(ipart);         // Geant3 particle type
      hits[1] = pos[0]+s*tx;            // X-position for hit
      hits[2] = pos[1]+s*ty;            // Y-position for hit
      hits[3] = pos[2]+s*tz;            // Z-position for hit
      hits[4] = theta;                  // theta angle of incidence
      hits[5] = phi;                    // phi angle of incidence 
      hits[8] = 0;//PadHits does not exist anymore  (Float_t) fNPadHits;    // first padhit
      hits[9] = -1;                     // last pad hit
      hits[10] = mom[3];                // hit momentum P
      hits[11] = mom[0];                // Px
      hits[12] = mom[1];                // Py
      hits[13] = mom[2];                // Pz
      tof=gMC->TrackTime();
      hits[14] = tof;                   // Time of flight
      tlength  = 0;
      eloss    = 0;
      eloss2   = 0;
      sstep=0;
      xhit     = pos[0];
      yhit     = pos[1];      
      zhit     = pos[2];      
      Chamber(idvol).ChargeCorrelationInit();
      // Only if not trigger chamber

//       printf("---------------------------\n");
//       printf(">>>> Y =  %f \n",hits[2]);
//       printf("---------------------------\n");
    
      

     //  if(idvol < AliMUONConstants::NTrackingCh()) {
//        //
//        //  Initialize hit position (cursor) in the segmentation model 
//        ((AliMUONChamber*) (*fChambers)[idvol])
//            ->SigGenInit(pos[0], pos[1], pos[2]);
//       } else {
//        //geant3->Gpcxyz();
//        //printf("In the Trigger Chamber #%d\n",idvol-9);
//       }
  }
  eloss2+=destep;
  sstep+=step;

  // cout<<sstep<<endl;

  // 
  // Calculate the charge induced on a pad (disintegration) in case 
  //
  // Mip left chamber ...
  if( gMC->IsTrackExiting() || gMC->IsTrackStop() || gMC->IsTrackDisappeared()){
      gMC->SetMaxStep(kBig);
      eloss   += destep;
      tlength += step;
      
      Float_t x0,y0,z0;
      Float_t localPos[3];
      Float_t globalPos[3] = {pos[0], pos[1], pos[2]};
      gMC->Gmtod(globalPos,localPos,1); 

      if(idvol < AliMUONConstants::NTrackingCh()) {
// tracking chambers
          x0 = 0.5*(xhit+pos[0]);
          y0 = 0.5*(yhit+pos[1]);
          z0 = 0.5*(zhit+pos[2]);
      } else {
// trigger chambers
          x0 = xhit;
          y0 = yhit;
          z0 = 0.;
      }
      

      //      if (eloss >0)  MakePadHits(x0,y0,z0,eloss,tof,idvol);
      
          
      hits[6] = tlength;   // track length
      hits[7] = eloss2;    // de/dx energy loss


      //      if (fNPadHits > (Int_t)hits[8]) {
      //          hits[8] = hits[8]+1;
      //          hits[9] = 0: // PadHits does not exist anymore (Float_t) fNPadHits;
      //}
//
//    new hit 
      
      new(lhits[fNhits++]) 
          AliMUONHit(fIshunt, gAlice->GetMCApp()->GetCurrentTrackNumber(), vol,hits);
      eloss = 0; 
      //
      // Check additional signal generation conditions 
      // defined by the segmentation
      // model (boundary crossing conditions)
      // only for tracking chambers
  } else if 
      ((idvol < AliMUONConstants::NTrackingCh()) &&
       ((AliMUONChamber*) (*fChambers)[idvol])->SigGenCond(pos[0], pos[1], pos[2]))
  {
      ((AliMUONChamber*) (*fChambers)[idvol])
          ->SigGenInit(pos[0], pos[1], pos[2]);
      
      Float_t localPos[3];
      Float_t globalPos[3] = {pos[0], pos[1], pos[2]};
      gMC->Gmtod(globalPos,localPos,1); 

      eloss    += destep;

      // if (eloss > 0 && idvol < AliMUONConstants::NTrackingCh())
      //        MakePadHits(0.5*(xhit+pos[0]),0.5*(yhit+pos[1]),pos[2],eloss,tof,idvol);
      xhit     = pos[0];
      yhit     = pos[1]; 
      zhit     = pos[2];
      eloss = 0;
      tlength += step ;
      //
      // nothing special  happened, add up energy loss
  } else {        
      eloss   += destep;
      tlength += step ;
  }
}
