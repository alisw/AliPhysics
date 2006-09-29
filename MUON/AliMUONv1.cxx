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
 * about the suitability of this software for any purpeateose. It is      *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

// --------------------
// Class AliMUONv1
// --------------------
// AliDetector class for MUON subsystem which implements
// functions for simulation 

#include "AliMUONv1.h"
#include "AliMUONConstants.h"
#include "AliMUONSegFactory.h"
#include "AliMUONResponseFactory.h"
#include "AliMUONSegmentation.h"
#include "AliMUONHit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONTriggerCircuitNew.h"
#include "AliMUONTriggerCrateStore.h"
#include "AliMUONGeometryBuilder.h"	
#include "AliMUONGeometry.h"	
#include "AliMUONGeometryTransformer.h"	
#include "AliMUONGeometryModule.h"	
#include "AliMUONStringIntMap.h"	
#include "AliMUONGeometryDetElement.h"	

#include "AliConst.h" 
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliTrackReference.h"
#include "AliLog.h"

#include <TRandom.h>
#include <TF1.h>
#include <TClonesArray.h>
#include <TRandom.h> 
#include <TVirtualMC.h>
#include <TGeoMatrix.h>

#include <string>

/// \cond CLASSIMP
ClassImp(AliMUONv1)
/// \endcond
 
//___________________________________________
AliMUONv1::AliMUONv1() 
  : AliMUON(),
    fAngleEffect(kTRUE),
    fStepMaxInActiveGas(0.6),
    fStepSum(0x0),
    fDestepSum(0x0),
    fTrackMomentum(), 
    fTrackPosition(),
    fElossRatio(0x0),
    fAngleEffect10(0x0),
    fAngleEffectNorma(0x0)
{
/// Default constructor
  
  AliDebug(1,Form("default (empty) ctor this=%p",this));
} 

//___________________________________________
AliMUONv1::AliMUONv1(const char *name, const char *title,
                     const char* sDigitizerClassName,
                     const char* digitizerClassName)
: AliMUON(name,title,sDigitizerClassName,digitizerClassName), 
    fAngleEffect(kTRUE),
    fStepMaxInActiveGas(0.6),
    fStepSum(0x0),
    fDestepSum(0x0),
    fTrackMomentum(), 
    fTrackPosition(),
    fElossRatio(0x0),
    fAngleEffect10(0x0),
    fAngleEffectNorma(0x0)
{
/// Standard onstructor

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

//___________________________________________
AliMUONv1::~AliMUONv1()
{
/// Destructor
  
  AliDebug(1,Form("dtor this=%p",this));
  delete [] fStepSum;
  delete [] fDestepSum;
  delete fElossRatio;
  delete fAngleEffect10;
  delete fAngleEffectNorma; 
}

//__________________________________________________
void AliMUONv1::CreateGeometry()
{
/// Construct geometry using geometry builder

  fGeometryBuilder->CreateGeometry();
}

//________________________________________________________________
void AliMUONv1::CreateMaterials()
{
/// Construct materials using geometry builder

  fGeometryBuilder->CreateMaterials();
}

//________________________________________________________________
void AliMUONv1::AddAlignableVolumes() const
{
/// Construct materials using geometry builder

  GetGeometryTransformer()->AddAlignableVolumes();
}


//___________________________________________
void AliMUONv1::Init()
{ 
/// Initialize geometry

  AliDebug(1,"Start Init for version 1 - CPC chamber type");
   
  fGeometryBuilder->InitGeometry();
  AliDebug(1,"Finished Init for version 1 - CPC chamber type");   
 

  std::string ftype(GetTitle());

  // Build segmentation
  // using geometry parametrisation
  //
  AliMUONSegFactory segFactory(GetGeometryTransformer());
  fSegmentation = segFactory.CreateSegmentation(ftype);

  if (!fSegmentation) {
    AliFatal(Form("Wrong factory type : %s",ftype.c_str()));
  }        

  // Build response
  //
  AliMUONResponseFactory respFactory("default");
  respFactory.Build(this);
  

  // Initialize segmentation
  //
  fSegmentation->Init();

  // Initialize trigger circuits
  //
  for (Int_t i=0; i<AliMUONConstants::NTriggerCircuit(); i++)  {
    AliMUONTriggerCircuit* c = (AliMUONTriggerCircuit*)(fTriggerCircuits->At(i));
    c->Init(i);
  }
  
  AliMUONTriggerCrateStore store;
  store.ReadFromFile();
  for (Int_t i=0; i<AliMUONConstants::NTriggerCircuit(); i++)  
  {
    AliMUONTriggerCircuitNew* c = (AliMUONTriggerCircuitNew*)(fTriggerCircuitsNew->At(i));
    c->Init(i,store);
  }
  
}

//__________________________________________________________________
Int_t  AliMUONv1::GetChamberId(Int_t volId) const
{
/// Check if the volume with specified  volId is a sensitive volume (gas) 
/// of some chamber and return the chamber number;
/// if not sensitive volume - return 0.

  for (Int_t i = 0; i < AliMUONConstants::NCh(); i++) {
    if ( GetGeometry()->GetModule(i)->IsSensitiveVolume(volId) )
      return i+1;
  }    

  return 0;
}

//_______________________________________________________________________________
TString  AliMUONv1::CurrentVolumePath() const
{
/// Return current volume path
/// (Could be removed when this function is available via gMC)

  TString path = "";
  TString name;
  Int_t copyNo;
  Int_t imother = 0;
  do {
    name = gMC->CurrentVolOffName(imother);
    gMC->CurrentVolOffID(imother++, copyNo);
    TString add = "/";
    add += name;
    add += "_";
    add += copyNo;
    path.Insert(0,add); 
  }
  while ( name != TString("ALIC") );
  
  return path;  
}

//_______________________________________________________________________________
void AliMUONv1::StepManager()
{
/// Step manager for the chambers

  // Only charged tracks
  if( !(gMC->TrackCharge()) ) return; 
  // Only charged tracks
  
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  static Int_t   idvol=-1, iEnter = 0;
  Int_t   iChamber=0;
  Int_t   id=0;
  Int_t   copy;
  const  Float_t kBig = 1.e10;
  static Double_t xyzEnter[3];

  //
  // Only gas gap inside chamber
  // Tag chambers and record hits when track enters 
  id=gMC->CurrentVolID(copy);
  iChamber = GetChamberId(id);
  idvol = iChamber -1;

  if (idvol == -1) return;
  
  // Detection elements id
  const AliMUONGeometryModule* kGeometryModule
    = GetGeometry()->GetModule(iChamber-1);

  AliMUONGeometryDetElement* detElement
    = kGeometryModule->FindBySensitiveVolume(CurrentVolumePath());

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
    cerr << "	global position: "
  	 << x << ", " << y << ", " << z
  	 << endl;
    AliError("DetElemId not identified.");
  }  
    
  // Filling TrackRefs file for MUON. Our Track references are the active volume of the chambers
  if ( (gMC->IsTrackEntering() || gMC->IsTrackExiting() ) ) {
    AliTrackReference* trackReference    
      = AddTrackReference(gAlice->GetMCApp()->GetCurrentTrackNumber());
    trackReference->SetUserId(detElemId);
  }  
  
  if( gMC->IsTrackEntering() ) {
     Float_t theta = fTrackMomentum.Theta();
     if ((TMath::Pi()-theta)*kRaddeg>=15.) gMC->SetMaxStep(fStepMaxInActiveGas); // We use Pi-theta because z is negative
     iEnter = 1;
     gMC->TrackPosition(xyzEnter[0], xyzEnter[1], xyzEnter[2]); // save coordinates of entrance point
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
    if (fDestepSum[idvol] == 0) {
      // AZ - no energy release
      fStepSum[idvol] = 0; // Reset for the next event
      iEnter = 0;
      return; 
    }

    gMC->TrackPosition(fTrackPosition);
    Float_t theta = fTrackMomentum.Theta();
    Float_t phi   = fTrackMomentum.Phi();
    
    Int_t merge = 0;
    Double_t xyz0[3]={0}, xyz1[3]={0}, tmp[3]={0};
    if (gMC->IsTrackExiting() && iEnter != 0) {
      // AZ - this code is to avoid artificial hit splitting at the
      // "fake" boundary inside the same chamber. It will still produce 
      // 2 hits but with the same coordinates (at the wire) to allow 
      // their merging at the digitization level.

      // Only for a track going from the entrance to the exit from the volume
      // Get local coordinates
      gMC->Gmtod(xyzEnter, xyz0, 1); // local coord. at the entrance

      fTrackPosition.Vect().GetXYZ(tmp);
      gMC->Gmtod(tmp, xyz1, 1); // local coord. at the exit
      Float_t dx = xyz0[0] - xyz1[0];
      Float_t dy = xyz0[1] - xyz1[1];
      Float_t thLoc = TMath::ATan2 (TMath::Sqrt(dx*dx+dy*dy), TMath::Abs(xyz0[2]-xyz1[2]));
      if (thLoc * TMath::RadToDeg() < 15) merge = 1; 
    }

    if (merge) {
      Double_t dz = -0.5;
      if (xyz1[2] != xyz0[2]) dz = xyz0[2] / (xyz1[2] - xyz0[2]);
      tmp[0] = xyz0[0] - (xyz1[0] - xyz0[0]) * dz; // local coord. at the wire
      tmp[1] = xyz0[1] - (xyz1[1] - xyz0[1]) * dz;
      tmp[2] = xyz0[2] - (xyz1[2] - xyz0[2]) * dz;
      gMC->Gdtom(tmp, xyz1, 1); // global coord. at the wire
      fTrackPosition.SetXYZT(xyz1[0], xyz1[1], xyz1[2], fTrackPosition.T());
    } else {
      TLorentzVector backToWire( fStepSum[idvol]/2.*sin(theta)*cos(phi),
				 fStepSum[idvol]/2.*sin(theta)*sin(phi),
				 fStepSum[idvol]/2.*cos(theta),0.0       );
      fTrackPosition-=backToWire;
      //printf(" %d %d %d %f %d \n", gMC->IsTrackExiting(), gMC->IsTrackStop(), gMC->IsTrackDisappeared(), fStepSum[idvol], iEnter);
      //    AliDebug(1,
      //     Form("Track Position %f %f %f",fTrackPosition.X(),fTrackPosition.Y(),fTrackPosition.Z()));
      // AliDebug(1,
      //     Form("Exit: Track backToWire %f %f %f",backToWire.X(),backToWire.Y(),backToWire.Z())) ;
    }
    
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
    iEnter = 0;
  }
}
