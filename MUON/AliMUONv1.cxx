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

#include "AliMUONv1.h"
#include "AliConst.h" 
#include "AliMUONChamber.h"
#include "AliMUONConstants.h"
#include "AliMUONFactory.h"
#include "AliMUONHit.h"
#include "AliMUONTriggerCircuit.h"
#include "AliMUONGeometryBuilder.h"	
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"

ClassImp(AliMUONv1)
 
//___________________________________________
AliMUONv1::AliMUONv1() 
  : AliMUON(),
    fStepManagerVersionOld(kFALSE),
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
// Default constructor
} 

//___________________________________________
AliMUONv1::AliMUONv1(const char *name, const char *title)
  : AliMUON(name,title), 
    fStepManagerVersionOld(kFALSE),
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
// Standard onstructor

    // By default include all stations
    AliMUONFactory factory;
    factory.Build(this, title);

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

  Fatal("AliMUONv1", "Copy constructor not provided.");
}

//___________________________________________
AliMUONv1::~AliMUONv1()
{
// Destructor

  delete [] fStepSum;
  delete [] fDestepSum;
  delete fElossRatio;
  delete fAngleEffect10;
  delete fAngleEffectNorma;  
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
   // 
   // Initialize Tracking Chambers
   //

   if(fDebug) printf("\n%s: Start Init for version 1 - CPC chamber type\n\n",ClassName());
   Int_t i;
   for (i=0; i<AliMUONConstants::NCh(); i++) {
       ( (AliMUONChamber*) (*fChambers)[i])->Init();
   }
   
   //
   // Initialize geometry
   //
   fGeometryBuilder->InitGeometry();
   if(fDebug) printf("\n%s: Finished Init for version 1 - CPC chamber type\n",ClassName());

   //cp 
   if(fDebug) printf("\n%s: Start Init for Trigger Circuits\n",ClassName());
   for (i=0; i<AliMUONConstants::NTriggerCircuit(); i++) {
     ( (AliMUONTriggerCircuit*) (*fTriggerCircuits)[i])->Init(i);
   }
   if(fDebug) printf("%s: Finished Init for Trigger Circuits\n",ClassName());
   //cp
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
