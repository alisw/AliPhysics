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
#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include "AliMUONv1.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONFactoryV3.h"
#include "AliMUONFactoryV2.h"
#include "AliMUONHit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryBuilder.h"	
#include "AliMUONGeometryModule.h"	
#include "AliMUONGeometrySVMap.h"	
#include "AliMUONGeometryDetElement.h"	
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"

#include <string>

ClassImp(AliMUONv1)
 
//___________________________________________
AliMUONv1::AliMUONv1() 
  : AliMUON(),
    fStepManagerVersionOld(kFALSE),
    fStepManagerVersionDE(kFALSE),
    fAngleEffect(kTRUE),
    fStepMaxInActiveGas(0.6),
    fStepSum(0x0),
    fDestepSum(0x0),
    fTrackMomentum(), 
    fTrackPosition(),
    fElossRatio(0x0),
    fAngleEffect10(0x0),
    fAngleEffectNorma(0x0),
    fFactory(0x0)
{
// Default constructor
	AliDebug(1,Form("default (empty) ctor this=%p",this));
} 

//___________________________________________
AliMUONv1::AliMUONv1(const char *name, const char *title)
  : AliMUON(name,title), 
    fStepManagerVersionOld(kFALSE),
    fStepManagerVersionDE(kFALSE),
    fAngleEffect(kTRUE),
    fStepMaxInActiveGas(0.6),
    fStepSum(0x0),
    fDestepSum(0x0),
    fTrackMomentum(), 
    fTrackPosition(),
    fElossRatio(0x0),
    fAngleEffect10(0x0),
    fAngleEffectNorma(0x0),
    fFactory(0x0)
{
// Standard onstructor

	AliDebug(1,Form("ctor this=%p",this));	
	
    // By default include all stations

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
}

//_____________________________________________________________________________
AliMUONv1::AliMUONv1(const AliMUONv1& right) 
  : AliMUON(right) 
{  
  // copy constructor (not implemented)

  AliFatal("Copy constructor not provided.");
}

//___________________________________________
AliMUONv1::~AliMUONv1()
{
// Destructor
	AliDebug(1,Form("dtor this=%p",this));
  delete [] fStepSum;
  delete [] fDestepSum;
  delete fElossRatio;
  delete fAngleEffect10;
  delete fAngleEffectNorma; 
  delete fFactory;
}

//_____________________________________________________________________________
AliMUONv1& AliMUONv1::operator=(const AliMUONv1& right)
{
  // assignement operator (not implemented)

  // check assignement to self
  if (this == &right) return *this;

  AliFatal("Assignement operator not provided.");
    
  return *this;  
}    

//__________________________________________________
void AliMUONv1::CreateGeometry()
{
//
// Construct geometry using geometry builder
//

  fGeometryBuilder->CreateGeometry();
}

//________________________________________________________________
void AliMUONv1::CreateMaterials()
{
//
// Construct materials using geometry builder
//

  fGeometryBuilder->CreateMaterials();
}

//___________________________________________
void AliMUONv1::Init()
{ 
  AliDebug(1,"Start Init for version 1 - CPC chamber type");
  Int_t i;
   
  //
  // Initialize geometry
  //
  fGeometryBuilder->InitGeometry();
  AliDebug(1,"Finished Init for version 1 - CPC chamber type");   
  
  std::string ftype(GetTitle());
  
  if ( ftype == "default" || ftype == "AliMUONFactoryV2") 
    {
      fFactory = new AliMUONFactoryV2("New MUON Factory");
      (static_cast<AliMUONFactoryV2*>(fFactory))->Build(this,"default");
    }
  else if ( ftype == "AliMUONFactoryV3" )
    {
      fFactory = new AliMUONFactoryV3("New MUON Factory");
      (static_cast<AliMUONFactoryV3*>(fFactory))->Build(this,"default");
    }
  else
    {
      AliFatal(Form("Wrong factory type : %s",ftype.c_str()));      
    }
      
  //
  // Initialize segmentation
  //
  
  for (i=0; i<AliMUONConstants::NCh(); i++) 
    ( (AliMUONChamber*) (*fChambers)[i])->Init(2);// new segmentation
  
  
  // trigger circuit
  // cp 
  for (i=0; i<AliMUONConstants::NTriggerCircuit(); i++) 
    ( (AliMUONTriggerCircuit*) (*fTriggerCircuits)[i])->Init(i);
  
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

//_______________________________________________________________________________
TString  AliMUONv1::CurrentVolumePath() const
{
// Returns current volume path
// (Could be removed when this function is available via gMC)
// ---	     

  TString path = "";
  TString name;
  Int_t copyNo;
  Int_t imother = 0;
  do {
    name = gMC->CurrentVolOffName(imother);
    gMC->CurrentVolOffID(imother++, copyNo);
    TString add = "/";
    add += name;
    add += ".";
    add += copyNo;
    path.Insert(0,add); 
  }
  while ( name != TString("ALIC") );
  
  return path;  
}

//_______________________________________________________________________________
void AliMUONv1::StepManager()
{
  // Stepmanager for the chambers
  // TBR

 if (fStepManagerVersionDE) {
    StepManager2();
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

   //  AliDebug(1,
   //    Form("Active volume found %d chamber %d Z chamber is %f ",idvol,iChamber,
   //    ( (AliMUONChamber*)(*fChambers)[idvol])->Z()));
   // Particule id and mass, 
  Int_t     ipart = gMC->TrackPid();
  Float_t   mass  = gMC->TrackMass();

  fDestepSum[idvol]+=gMC->Edep();
  // Get current particle id (ipart), track position (pos)  and momentum (mom)
  if ( fStepSum[idvol]==0.0 )  gMC->TrackMomentum(fTrackMomentum);
  fStepSum[idvol]+=gMC->TrackStep();
  
  //  if(AliDebugLevel()) {
  //  AliDebug(1,
  //	     Form("iChamber %d, Particle %d, theta %f phi %f mass %f StepSum %f eloss %g",
  //		  iChamber,ipart, fTrackMomentum.Theta()*kRaddeg, fTrackMomentum.Phi()*kRaddeg, 
  //		  mass, fStepSum[idvol], gMC->Edep()));
  //  AliDebug(1,
  //	     Form("Track Momentum %f %f %f", fTrackMomentum.X(), fTrackMomentum.Y(), 
  //		  fTrackMomentum.Z()));
  //  gMC->TrackPosition(fTrackPosition);
  //  AliDebug(1,
  //	     Form("Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),
  //	  fTrackPosition.Z())) ;
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
    //     AliDebug(1,Form("Exit: Track Position %f %f %f",fTrackPosition.X(),
    //                     fTrackPosition.Y(),fTrackPosition.Z())) ;
    //     AliDebug(1,Form("Exit: Track backToWire %f %f %f",backToWire.X(),
    //                     backToWire.Y(),backToWire.Z()) ;
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
    
    // Detection elements ids
    AliMUONGeometryModule* geometry
      = Chamber(iChamber-1).GetGeometry();

    AliMUONGeometryDetElement* detElement
      = geometry->FindBySensitiveVolume(CurrentVolumePath());

    Int_t detElemId = 0;
    if (detElement) detElemId = detElement->GetUniqueID(); 
 
    if (!detElemId) {
      cerr << "Chamber id: "
           << setw(3) << iChamber << "  "
           << "Current SV: " 
           <<  CurrentVolumePath() 
 	   << "  detElemId: "
           << setw(5) << detElemId 
	   << endl;
      Double_t x, y, z;
      gMC->TrackPosition(x, y, z);	   
      cerr << "   global position: "
           << x << ", " << y << ", " << z
           << endl;
      AliWarning("DetElemId not identified.");
    }  
    
    // One hit per chamber
    GetMUONData()->AddHit(fIshunt, 
			  gAlice->GetMCApp()->GetCurrentTrackNumber(), 
			  iChamber, ipart,
			  fTrackPosition.X(), 
			  fTrackPosition.Y()+yAngleEffect, 
			  fTrackPosition.Z(), 
			  gMC->TrackTime(),
			  fTrackMomentum.P(),
			  theta, 
			  phi, 
			  fStepSum[idvol], 
			  fDestepSum[idvol],                        
			  fTrackPosition.X(),
			  fTrackPosition.Y(),
			  fTrackPosition.Z());

//     if (AliDebugLevel()){
//       AliDebug(1,Form("Exit: Particle exiting from chamber %d",iChamber));
//       AliDebug(1,Form("Exit: StepSum %f eloss geant %g ",fStepSum[idvol],fDestepSum[idvol]));
//       AliDebug(1,Form("Exit: Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z())) ;
//     }
    fStepSum[idvol]  =0; // Reset for the next event
    fDestepSum[idvol]=0; // Reset for the next event
  }
}

//_______________________________________________________________________________
void AliMUONv1::StepManager2()
{
  // Stepmanager for the chambers

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

   //   AliDebug(1,
   //	    Form("Active volume found %d chamber %d Z chamber is %f ",idvol,iChamber,
   //		 ( (AliMUONChamber*)(*fChambers)[idvol])->Z())) ;
  // Particule id and mass, 
  Int_t     ipart = gMC->TrackPid();
  Float_t   mass  = gMC->TrackMass();

  fDestepSum[idvol]+=gMC->Edep();
  // Get current particle id (ipart), track position (pos)  and momentum (mom)
  if ( fStepSum[idvol]==0.0 )  gMC->TrackMomentum(fTrackMomentum);
  fStepSum[idvol]+=gMC->TrackStep();
  
  //  if (AliDebugLevel()) {
  //   AliDebug(1,Form("Step, iChamber %d, Particle %d, theta %f phi %f mass %f StepSum %f eloss %g",
  //		     iChamber,ipart, fTrackMomentum.Theta()*kRaddeg, fTrackMomentum.Phi()*kRaddeg,
  //	     mass, fStepSum[idvol], gMC->Edep()));
  // AliDebug(1,Form("Step:Track Momentum %f %f %f", fTrackMomentum.X(), fTrackMomentum.Y(), 
  //	     fTrackMomentum.Z()));
  // gMC->TrackPosition(fTrackPosition);
  // AliDebug(1,Form("Step: Track Position %f %f %f",fTrackPosition.X(),
  //	     fTrackPosition.Y(),fTrackPosition.Z())) ;
  //}

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
    //    AliDebug(1,
    //	     Form("Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()));
    // AliDebug(1,
    //	     Form("Exit: Track backToWire %f %f %f",backToWire.X(),backToWire.Y(),backToWire.Z())) ;
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
    
    // Detection elements ids
    AliMUONGeometryModule* geometry
      = Chamber(iChamber-1).GetGeometry();

    AliMUONGeometryDetElement* detElement
      = geometry->FindBySensitiveVolume(CurrentVolumePath());

    Int_t detElemId = 0;
    if (detElement) detElemId = detElement->GetUniqueID(); 
 
    if (!detElemId) {
      cerr << "Chamber id: "
           << setw(3) << iChamber << "  "
           << "Current SV: " 
           <<  CurrentVolumePath() 
 	   << "  detElemId: "
           << setw(5) << detElemId 
	   << endl;
      Double_t x, y, z;
      gMC->TrackPosition(x, y, z);	   
      cerr << "   global position: "
           << x << ", " << y << ", " << z
           << endl;
      AliError("DetElemId not identified.");
    }  
    
    // One hit per chamber
    GetMUONData()->AddHit2(fIshunt, 
			  gAlice->GetMCApp()->GetCurrentTrackNumber(), 
			  detElemId, ipart,
			  fTrackPosition.X(), 
			  fTrackPosition.Y()+yAngleEffect, 
			  fTrackPosition.Z(), 
			  gMC->TrackTime(),
			  fTrackMomentum.P(),
			  theta, 
			  phi, 
			  fStepSum[idvol], 
			  fDestepSum[idvol],                        
			  fTrackPosition.X(),
			  fTrackPosition.Y(),
			  fTrackPosition.Z());

    //       AliDebug(1,Form("Exit: Particle exiting from chamber %d",iChamber));
    //       AliDebug(1,Form("Exit: StepSum %f eloss geant %g ",fStepSum[idvol],fDestepSum[idvol]));
    //       AliDebug(1,Form("Exit: Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()) ;

    fStepSum[idvol]  =0; // Reset for the next event
    fDestepSum[idvol]=0; // Reset for the next event
  }
}
