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
*/

/////////////////////////////////////////////////////////
//  Manager and hits classes for set:CPV version 0     //
//  Coarse geometry                                    //
//                                                     //
//  Author: Yuri Kharlov, IHEP, Protvino               //
//  e-mail: Yuri.Kharlov@cern.ch                       //
//  Last modified: 17 September 1999                   //
/////////////////////////////////////////////////////////
 
// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TBRIK.h"
#include "TNode.h"

// --- galice header files ---
#include "AliCPVv0.h"
#include "AliRun.h"
#include "AliMC.h" 

ClassImp(AliCPVv0)

//==============================================================================
//                            AliCPVv0
//==============================================================================

AliCPVv0::AliCPVv0()
{
}
 
//______________________________________________________________________________

AliCPVv0::AliCPVv0(const char *name, const char *title)
          : AliCPV(name, title)
{
}
 
//______________________________________________________________________________

void AliCPVv0::CreateGeometry()
{

  AliCPV *CPV_tmp = (AliCPV*)gAlice->GetModule("CPV");
  if( NULL==CPV_tmp )
  {
    printf("There isn't CPV detector!\n");
    return;
  }

  Int_t    rotation_matrix_number=0;
  Float_t  par[3],
           x,y,z;

  // CPV creation

  par[0] = GetPadZSize()   / 2 * GetNz();
  par[1] = GetPadPhiSize() / 2 * GetNphi();
  par[2] = GetThickness()  / 2;
  gMC->Gsvolu("CPV","BOX ",GetCPVIdtmed(),par,3);

  for( Int_t i=0; i<GetNofCradles(); i++ )
  {
    Float_t cradle_angle_pos = -90+(i-(GetNofCradles()-1)/2.) * GetAngle();

    // Cradles are numerated in clock reversed order. (general way of angle increment)

    Float_t r = GetRadius() + GetThickness()/2;
    x = r*cos(cradle_angle_pos*kPI/180);
    y = r*sin(cradle_angle_pos*kPI/180);
    z = 0;
    AliMatrix(rotation_matrix_number, 0,0 , 90,90+cradle_angle_pos , 90,180+cradle_angle_pos);
    gMC->Gspos("CPV",i+1,"ALIC",x,y,z,rotation_matrix_number,"ONLY");
  }
  AddCPVCradles();

}

//______________________________________________________________________________

void AliCPVv0::StepManager()
{

//    if( gMC->IsTrackEntering() ) {
//      const char *VolumeName = gMC->CurrentVolName();
//      cout << "AliCPVv0::StepManager() entered to CPV to the volume " 
//           << VolumeName << "!\n";
//    }
  
  if( strcmp(gMC->CurrentVolName(),"CPV ")==0 && gMC->IsTrackEntering() )
  {
    // GEANT particle just have entered into CPV detector.

    AliCPV &CPV = *(AliCPV*)gAlice->GetModule("CPV");

    Int_t CradleNumber;
    gMC->CurrentVolOffID(0,CradleNumber);
    CradleNumber--;
    AliCPVCradle  &cradle = CPV.GetCradle(CradleNumber);

    // Current position of the hit in the cradle ref. system

    TLorentzVector xyzt;
    gMC -> TrackPosition(xyzt);
    Float_t xyzm[3], xyzd[3], xyd[2];
    for (Int_t i=0; i<3; i++) xyzm[i] = xyzt[i];
    gMC -> Gmtod (xyzm, xyzd, 1);
    for (Int_t i=0; i<2; i++) xyd[i]  = xyzd[i];

    // Current momentum of the hit's track in the MARS ref. system
    
    TLorentzVector  pmom;
    gMC -> TrackMomentum(pmom);
    Int_t ipart = gMC->TrackPid();

    // Store current particle in the list of Cradle particles.
    cradle.AddHit(pmom,xyd,ipart);

//      if (gMC->TrackCharge()!=0.) CPV.Reconstruction(0,0);

//      CPV.Print("p");
//      CPV.Print("r");

  }
}
