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

//
// Class for the display of ITS
// Version 11
// with the new geometry
// A proper description of this class
// will be written shortly
//


/*
$Id$
*/

#include <Riostream.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TControlBar.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoNode.h>
#include <TGeoVolume.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TView.h>
class TGeoTube;

#include "AliITS.h"
#include "AliITSv11.h"
#include "DisplayITSv11.h"

ClassImp(DisplayITSv11)
//______________________________________________________________________
DisplayITSv11::DisplayITSv11()
{
  // Default Constructor.
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    A Default constructed DispalyITSv11 class/task.
  
  fmgr         = 0;
  fits         = 0;
  fALICE       = 0;
  fITS         = 0;
  fClip        = 0;
  fITSdebug    = 0;
  fNsegments   = 80;
  fCut         = 0;
  fAxis        = 1;
  fPerspective = 0;
  fSolid       = 0;
  fPhimincut   = 0.0;
  fPhimaxcut   = 180.0;
  fRmin[0] = fRmin[1] = fRmin[2] = fRmax[0] = fRmax[1] = fRmax[2] = 0.0;
  fLongitude   = 90.0;
  fLatitude    = 0.0;
  fPsi         = 0.0;
}

//______________________________________________________________________
DisplayITSv11::~DisplayITSv11()
{
  // Destructor.
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Return:
  //    none..
  
  if(fits) delete fits;
  if(fClip) delete fClip;
  fITS   = 0;
  fALICE = 0;
  fmgr   = 0;
}

//______________________________________________________________________
void DisplayITSv11::Exec(Option_t* opt)
{
  // Main display routine
  // Inputs:
  //      Option_t * opt   Not presently used.
  // Outputs:
  //      none.
  // Return:
  //      none.
  Int_t i;
  
  if(strstr(opt,"debug")) fits->SetDebug();
  gSystem->Load("libGeom");
  //
  if(gGeoManager) delete gGeoManager;
  fmgr = gGeoManager = new TGeoManager("ITSGeometry",
				       " ITS Simulation Geometry Manager");
  //
  TGeoMaterial *vacmat = new TGeoMaterial("Vacume",0,0,0);
  TGeoMedium   *vacmed = new TGeoMedium("Vacume_med",1,vacmat);
  fRmax[0] = 100.;
  fRmax[1] = 100.;
  fRmax[2] = 200.;
  for(i=0;i<3;i++) fRmin[i] = -fRmax[i];
  fALICE = fmgr->MakeBox("ALIC",vacmed,fRmax[0],fRmax[1],fRmax[2]);
  fmgr->SetTopVolume(fALICE);
  //
  AliITSv11 *fits = new AliITSv11();
  fits->SetDebug(GetDebugITS());
  fits->CreateMaterials();
  fits->CreateGeometry();
  //
  fmgr->CloseGeometry();
  fITS = fALICE->FindNode("ITSV_1")->GetVolume();
  //
  TControlBar *bar=new TControlBar("vertical","ITS Geometry Display",10,10);
  bar->AddButton("Set ITS Debug level 1","SetITSdebugOn()","Debug on");
  bar->AddButton("Set ITS Debug level 0","SetITSdebugOff()","Debug off");
  bar->AddButton("Set ITS Theta,Phi cut on","SetCylindericalClipVolume()",
		 "Cut on");
  bar->AddButton("Set ITS Theta,Phi cut off","SetCylindericalCutOff()",
		 "Cut off");
  bar->AddButton("Set axis on","SetAxisOn()","Show Axis on");
  bar->AddButton("Set axis off","SetAxisOff()","Show Axis off");
  bar->AddButton("Set perspective on","SetPerspectiveOn()",
		 "Perspective on");
  bar->AddButton("Set perspective off","SetPerspectiveOff()",
		 "Perspective off");
  bar->AddButton("Draw volumes as Solid","SetSolid()","Solid Volumes");
  bar->AddButton("Draw volumes as wire","SetWire()","Wire Volumes");
  bar->AddButton("Set circle/80","SetCircleSegments(80)",
		 "circles ~ by 80 lines");
  bar->AddButton("Check Overlaps","fmgr->CheckOverlaps()",
		 "Check for overlaps");
  bar->AddButton("Display Geometry","Displayit()","Run Displayit");
  bar->AddButton("Display SPD Thermal Sheald","EngineeringSPDThS()",
		 "Run EngineeringSPDThS");
  bar->AddButton("Display SDD Cone","EngineeringSDDCone()",
		 "Run EngineeringSDDCone");
  bar->AddButton("Display SDD Centeral Cylinder","EngineeringSDDCylinder()",
		 "Run EngineeringSDDCylinder");
  bar->AddButton("Display SUP RB24 side","EngineeringSupRB24()",
		 "Run EngineeringSDDCylinder");
  bar->AddButton("Display SUP RB26 side","EngineeringSupRB26()",
		 "Run EngineeringSupRB26");
  bar->AddButton("Quit/Exit",".q","Exit");
  bar->Show();
  gROOT->SaveContext();
  //Displayit();
}

//______________________________________________________________________
void DisplayITSv11::Displaying(TGeoVolume *v,TCanvas *c,Int_t ipad)
{
  // Display Volume according to existing values
  // Inputs:
  //    TGeoVolume *v volume to be drawn.
  //    TCanvas    *p Pad where drawing is to be done
  //    Int_t      ipad subpad to draw on
  // Outputs:
  //    none.
  // Return:
  //    none.
  Int_t irr;
  TVirtualPad *p = 0;
  TView *view=0;
  
  c->cd(ipad);
  p = c->GetPad(ipad);
  if(!p) return;
  view = p->GetView();
  if(!view){view = new TView(fRmin,fRmax,1);p->SetView(view);}
  view->SetRange(fRmin,fRmax);
  view->SetView(fLongitude,fLatitude,fPsi,irr);
  if(fPerspective) view->SetPerspective();
  else view->SetParralel();
  if(fAxis) view->ShowAxis();
  //
  if(fClip) fmgr->SetClippingShape(fClip);
  v->Draw();
  if(fSolid) v->Raytrace();
  p->Modified(1);
  p->Update();
}
//----------------------------------------------------------------------
void DisplayITSv11::DisplayITS(){
    // Display AliITSv11 Geometry
    // Inputs:
    //    const char* filename output file with the display in it
    // Outputs:
    //    none.
    // Retrurn:
    //    none.
    Double_t lon,lat,psi;
    //
    TCanvas *c1;
    if(!(c1 = (TCanvas*)gROOT->FindObject("C1")))
        c1 = new TCanvas("C1","ITS Simulation Geometry",900,900);
    c1->Divide(2,2);
    //
    fmgr->SetNsegments(fNsegments);
    //
    fmgr->SetVisLevel(6);
    fmgr->SetVisOption(0);
    // Perspective
    Displaying(fALICE,c1,2);
    lon = fLongitude;
    lat = fLatitude;
    psi = fPsi;
    fLongitude = 270.0; // Front
    fLatitude  =  90.0;
    fPsi       =   0.0;
    Displaying(fALICE,c1,1);
    fLongitude = 270.0; // Top
    fLatitude  =   0.0;
    fPsi       =   0.0;
    Displaying(fALICE,c1,3);
    fLongitude =   0.0; // Side
    fLatitude  =  90.0;
    fPsi       =   0.0;
    Displaying(fALICE,c1,4);
    fLongitude = lon;
    fLatitude  = lat;
    fPsi       = psi;
    //
}
//----------------------------------------------------------------------
void DisplayITSv11::EngineeringSPDThS()
{
  // Display SPD Thermal Sheald Geometry
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Retrurn:
  //    none.
  Double_t lon,lat,psi;
  //
  TCanvas *c4;
  TGeoNode *node;
  if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
    c4 = new TCanvas("C4","ITS SDD Cylinder Geometry",900,450);
  c4->Divide(2,1);
  TGeoVolume *sPDThS=0;
  //TArrow *arrow=new TArrow();
  //
  node = fITS->FindNode("ITSspdThermalSheald_1");
  sPDThS = node->GetVolume();
  //
  fmgr->SetNsegments(fNsegments);
  //
  fmgr->SetVisLevel(6);
  fmgr->SetVisOption(0);
  //
  lon = fLongitude;
  lat = fLatitude;
  psi = fPsi;
  fLongitude = 270.0; // Front
  fLatitude  =  90.0;
  fPsi       =   0.0;
  Displaying(sPDThS,c4,1);
  //
  fLongitude = 270.0; // Top
  fLatitude  =   0.0;
  fPsi       =   0.0;
  Displaying(sPDThS,c4,2);
  //
  fLongitude = lon;
  fLatitude  = lat;
  fPsi       = psi;
}

//----------------------------------------------------------------------
void DisplayITSv11::EngineeringSDDCone()
{
  // Display SDD Cone Geometry
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Retrurn:
  //    none.
  Double_t lon,lat,psi;
  //
  TCanvas *c2;
  if(!(c2 = (TCanvas*)gROOT->FindObject("C2")))
    c2 = new TCanvas("C2","ITS SDD Cone Geometry",900,450);
  c2->Divide(2,1);
  TGeoVolume *sDD=0;
  TGeoNode *node;
  //
  node = fITS->FindNode("ITSsddConeL_1");
  sDD = node->GetVolume();
  //
  fmgr->SetNsegments(fNsegments);
  //
  fmgr->SetVisLevel(6);
  fmgr->SetVisOption(0);
  //
  lon = fLongitude;
  lat = fLatitude;
  psi = fPsi;
  fLongitude = 270.0; // Front
  fLatitude  =  90.0;
  fPsi       =   0.0;
  Displaying(sDD,c2,1);
  //
  fLongitude = 270.0; // Top
  fLatitude  =   0.0;
  fPsi       =   0.0;
  Displaying(sDD,c2,2);
  //
  fLongitude = lon;
  fLatitude  = lat;
  fPsi       = psi;
}

//----------------------------------------------------------------------
void DisplayITSv11::EngineeringSDDCylinder()
{
  // Display SDD Cylinder Geometry
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Retrurn:
  //    none.
  Double_t lon,lat,psi;
  //
  TCanvas *c3;
  if(!(c3 = (TCanvas*)gROOT->FindObject("C3")))
    c3 = new TCanvas("C3","ITS SDD Cylinder Geometry",900,450);
  c3->Divide(2,1);
  TGeoVolume *sDD=0;
  TGeoNode *node;
  TArrow *arrow=new TArrow();
  //
  node = fITS->FindNode("ITSsddCentCylCF_1");
  sDD = node->GetVolume();
  Double_t rmin = ((TGeoTube*)(sDD->GetShape()))->GetRmin();
  Double_t rmax = ((TGeoTube*)(sDD->GetShape()))->GetRmax();
  Double_t dz   = ((TGeoTube*)(sDD->GetShape()))->GetDz();
  //
  fmgr->SetNsegments(fNsegments);
  //
  fmgr->SetVisLevel(6);
  fmgr->SetVisOption(0);
  //
  lon = fLongitude;
  lat = fLatitude;
  psi = fPsi;
  fLongitude = 270.0; // Front
  fLatitude  =  90.0;
  fPsi       =   0.0;
  Displaying(sDD,c3,1);
  arrow->DrawArrow(1.01*rmax,-dz,1.01*rmax,+dz);
  //
  fLongitude = 270.0; // Top
  fLatitude  =   0.0;
  fPsi       =   0.0;
  Displaying(sDD,c3,2);
  arrow->DrawArrow(rmax,0.0,rmax,0.0);
  Double_t s = TMath::Sin(0.7),c = TMath::Cos(0.7);
  arrow->DrawArrow(-rmin*c,-rmin*s,rmin*c,rmin*s);
  //
  fLongitude = lon;
  fLatitude  = lat;
  fPsi       = psi;
}

//----------------------------------------------------------------------
void DisplayITSv11::EngineeringSupRB24()
{
  // Display SDD Cylinder Geometry
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Retrurn:
  //    none.
  Double_t lon,lat,psi;
  //
  TCanvas *c4;
  if(!(c4 = (TCanvas*)gROOT->FindObject("C4")))
    c4 = new TCanvas("C4","ITS SDD Cylinder Geometry",900,450);
  c4->Divide(2,1);
  TGeoVolume *sUPRB24=0;
  TGeoNode *node;
  //TArrow *arrow=new TArrow();
  //
  node = fITS->FindNode("ITSsupFrameM24_1");
  sUPRB24 = node->GetVolume();
  //
  fmgr->SetNsegments(fNsegments);
  //
  fmgr->SetVisLevel(6);
  fmgr->SetVisOption(0);
  //
  lon = fLongitude;
  lat = fLatitude;
  psi = fPsi;
  fLongitude = 270.0; // Front
  fLatitude  =  90.0;
  fPsi       =   0.0;
  Displaying(sUPRB24,c4,1);
  //
  fLongitude = 270.0; // Top
  fLatitude  =   0.0;
  fPsi       =   0.0;
  Displaying(sUPRB24,c4,2);
  //
  fLongitude = lon;
  fLatitude  = lat;
  fPsi       = psi;
}

//----------------------------------------------------------------------
void DisplayITSv11::EngineeringSupRB26()
{
  // Display SDD Cylinder Geometry
  // Inputs:
  //    none.
  // Outputs:
  //    none.
  // Retrurn:
  //    none.
  Double_t lon,lat,psi;
  //
  TCanvas *c5;
  if(!(c5 = (TCanvas*)gROOT->FindObject("C5")))
    c5 = new TCanvas("C5","ITS SDD Cylinder Geometry",900,450);
  c5->Divide(2,1);
  TGeoVolume *sUPRB26=0;
  TGeoNode *node;
  //TArrow *arrow=new TArrow();
  //
  node = fITS->FindNode("ITSsupFrameM26_1");
  sUPRB26 = node->GetVolume();
  //
  fmgr->SetNsegments(fNsegments);
  //
  fmgr->SetVisLevel(6);
  fmgr->SetVisOption(0);
  //
  lon = fLongitude;
  lat = fLatitude;
  psi = fPsi;
  fLongitude = 270.0; // Front
  fLatitude  =  90.0;
  fPsi       =   0.0;
  Displaying(sUPRB26,c5,1);
  //
  fLongitude = 270.0; // Top
  fLatitude  =   0.0;
  fPsi       =   0.0;
  Displaying(sUPRB26,c5,2);
  //
  fLongitude = lon;
  fLatitude  = lat;
  fPsi       = psi;
}
