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
//  Manager and hits classes for set:PHOS version 1    //
/////////////////////////////////////////////////////////
 
// --- ROOT system ---
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

// --- galice header files ---
#include "AliPHOSv1.h"
#include "AliRun.h"

ClassImp(AliPHOSv1)

//______________________________________________________________________________


AliPHOSv1::AliPHOSv1()
{
}
 
//______________________________________________________________________________

AliPHOSv1::AliPHOSv1(const char *name, const char *title)
          : AliPHOS(name, title)
{
}
 
//___________________________________________
void AliPHOSv1::CreateGeometry()
{

  AliPHOS *PHOS_tmp = (AliPHOS*)gAlice->GetModule("PHOS");
  if( NULL==PHOS_tmp )
  {
    printf("There isn't PHOS detector!\n");
    return;
  }
//  AliPHOS &PHOS = *PHOS_tmp;

  //////////////////////////////////////////////////////////////////////////////

  Int_t                 rotation_matrix_number=0;
  Float_t               par[11],
                        x,y,z;

  const float           cell_length             = GetCrystalLength()+GetAirThickness()+GetWrapThickness()+GetPIN_Length(),
                        cell_side_size          = GetCrystalSideSize()+2*GetAirThickness()+2*GetWrapThickness(),
                        cradle_thikness         = cell_length;

  //////////////////////////////////////////////////////////////////////////////
  // CELL volume and subvolumes creation
  //////////////////////////////////////////////////////////////////////////////

  par[0] = GetCrystalSideSize()/2 + GetWrapThickness();
  par[1] = GetCrystalSideSize()/2 + GetWrapThickness();
  par[2] = GetCrystalLength()  /2 + GetWrapThickness()/2;
  gMC->Gsvolu("WRAP","BOX ",GetPHOS_IDTMED_Tyvek(),par,3);

  par[0] = GetCrystalSideSize()/2;
  par[1] = GetCrystalSideSize()/2;
  par[2] = GetCrystalLength()/2;
  gMC->Gsvolu("CRST","BOX ",GetPHOS_IDTMED_PbWO4(),par,3);

  // PIN
  par[0] = GetPIN_SideSize()/2;
  par[1] = GetPIN_SideSize()/2;
  par[2] = GetPIN_Length()/2;
  gMC->Gsvolu("PIN ","BOX ",GetPHOS_IDTMED_PIN(),par,3);

  //////////////////////////////////////////////////////////////////////////////
  // CRADLE creation.
  //////////////////////////////////////////////////////////////////////////////

  par[0] = cell_side_size/2 * GetNz();
  par[1] = cell_side_size/2 * GetNphi();
  par[2] = cradle_thikness/2;
  gMC->Gsvolu("PHOS","BOX ",GetPHOS_IDTMED_AIR(),par,3);


  par[0] = cell_side_size/2 * GetNz();
  par[1] = cell_side_size/2 * GetNphi();
  par[2] = cell_length/2;
  gMC->Gsvolu("CRS0","BOX ",GetPHOS_IDTMED_AIR(),par,3);

  x = 0;
  y = 0;
  z = -(cradle_thikness-cell_length)/2;
  gMC->Gspos("CRS0",1,"PHOS",x,y,z,0,"ONLY");

  gMC->Gsdvn("CRS1","CRS0",GetNphi(),2);
  gMC->Gsdvn("CELL","CRS1",GetNz()  ,1);

  //////////////////////////////////////////////////////////////////////////////
  // CELL creation
  //////////////////////////////////////////////////////////////////////////////

  x = 0;
  y = 0;
  z = -GetWrapThickness()/2;
  gMC->Gspos("CRST",1,"WRAP",x,y,z,0,"ONLY");

  x = 0;
  y = 0;
  z = GetPIN_Length()/2;
  gMC->Gspos("WRAP",1,"CELL",x,y,z,0,"ONLY");

  x = 0;
  y = 0;
  z = -GetCrystalLength()/2-GetWrapThickness()/2;
  gMC->Gspos("PIN ",1,"CELL",x,y,z,0,"ONLY");

  //////////////////////////////////////////////////////////////////////////////
  // CELL has been created.
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // End of CRADLE creation.
  //////////////////////////////////////////////////////////////////////////////


  //////////////////////////////////////////////////////////////////////////////
  // PHOS creation
  //////////////////////////////////////////////////////////////////////////////

  for( int i=0; i<GetCradlesAmount(); i++ )
  {
    Float_t cradle_angle     = 27.,
            cradle_angle_pos = -90+(i-(GetCradlesAmount()-1)/2.) *
                               (cradle_angle+GetAngleBetweenCradles());
    // Cradles are numerated in clock reversed order. (general way of angle increment)

    Float_t r = GetRadius() + cradle_thikness/2;
    x = r*cos(cradle_angle_pos*kPI/180);
    y = r*sin(cradle_angle_pos*kPI/180);
    z = 0;
    AliMatrix(rotation_matrix_number, 0,0 , 90,90+cradle_angle_pos , 90,180+cradle_angle_pos);
    gMC->Gspos("PHOS",i+1,"ALIC",x,y,z,rotation_matrix_number,"ONLY");

    GetCradleAngle(i) = cradle_angle_pos;
  }
  AddPHOSCradles();

  //////////////////////////////////////////////////////////////////////////////
  // All is done.
  // Print some information.
  //////////////////////////////////////////////////////////////////////////////
}

void AliPHOSv1::StepManager()
{
  static Bool_t inwold=0;   // Status of previous ctrak->inwvol
  Int_t copy;

  int cradle_number, cell_Z, cell_Phi;  // Variables that describe cell position.

  if( gMC->GetMedium() == GetPHOS_IDTMED_PIN() && (gMC->IsTrackInside() || gMC->IsTrackExiting()==2) && inwold && gMC->TrackCharge()!=0 )
  {
    // GEANT particle just have entered into PIN diode.

    AliPHOS &PHOS = *(AliPHOS*)gAlice->GetModule("PHOS");

    gMC->CurrentVolOffID(4,copy);
    cradle_number  = copy-1;
    gMC->CurrentVolOffID(1,copy);
    cell_Z         = copy-1;
    gMC->CurrentVolOffID(2,copy);
    cell_Phi       = copy-1;
/*
        cradle_number  = cvolu->number[cvolu->nlevel-5]-1;
        cell_Z         = cvolu->number[cvolu->nlevel-2]-1;
        cell_Phi       = cvolu->number[cvolu->nlevel-3]-1;
*/

    TH2S &h = PHOS.GetCradle(cradle_number).fChargedTracksInPIN;
    h.AddBinContent(h.GetBin(cell_Z,cell_Phi));
  }

  //////////////////////////////////////////////////////////////////////////////

  if( gMC->GetMedium() == GetPHOS_IDTMED_PbWO4() )
  {
    // GEANT particle into crystal.

    AliPHOS &PHOS = *(AliPHOS*)gAlice->GetModule("PHOS");

    gMC->CurrentVolOffID(5,copy);
    cradle_number  = copy-1;
    gMC->CurrentVolOffID(2,copy);
    cell_Z         = copy-1;
    gMC->CurrentVolOffID(3,copy);
    cell_Phi       = copy-1;
/*
        cradle_number  = cvolu->number[cvolu->nlevel-6]-1;
        cell_Z         = cvolu->number[cvolu->nlevel-3]-1;
        cell_Phi       = cvolu->number[cvolu->nlevel-4]-1;
*/
    TH2F &h = PHOS.GetCradle(cradle_number).fCellEnergy;
    h.AddBinContent(h.GetBin(cell_Z,cell_Phi),gMC->Edep());
  }

  //////////////////////////////////////////////////////////////////////////////


  inwold=gMC->IsTrackEntering();         // Save current status of GEANT variable.
}

