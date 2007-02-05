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
//  Transition Radiation Detector                                            //
//  This class contains the basic functions for the Transition Radiation     //
//  Detector.                                                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <Riostream.h>

#include <TClonesArray.h>
#include <TFile.h>
#include <TGeometry.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TNode.h>
#include <TPGON.h> 
#include <TParticle.h>
#include <TROOT.h>
#include <TTree.h>
#include <TVirtualMC.h>
 
#include "AliConst.h"
#include "AliDigit.h"
#include "AliLoader.h"
#include "AliLog.h"
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliTrackReference.h"

#include "AliTRD.h"
#include "AliTRDdigit.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDhit.h"
#include "AliTRDpoints.h"
#include "AliTRDrawData.h"
#include "AliTRDSimParam.h"
#include "AliTRDRecParam.h"
#include "AliTRDCommonParam.h"
#include "AliTRDcalibDB.h"

ClassImp(AliTRD)
 
//_____________________________________________________________________________
AliTRD::AliTRD()
  :AliDetector()
  ,fGeometry(0)
  ,fGasDensity(0)
  ,fFoilDensity(0)
  ,fDrawTR(0)
  ,fDisplayType(0)
{
  //
  // Default constructor
  //
 
}
 
//_____________________________________________________________________________
AliTRD::AliTRD(const char *name, const char *title)
  :AliDetector(name,title)
  ,fGeometry(0)
  ,fGasDensity(0)
  ,fFoilDensity(0)
  ,fDrawTR(0)
  ,fDisplayType(0)
{
  //
  // Standard constructor for the TRD
  //

  // Check that FRAME is there otherwise we have no place where to put TRD
  AliModule *frame = gAlice->GetModule("FRAME");
  if (!frame) {
    AliError("TRD needs FRAME to be present\n");
    exit(1);
  } 

  // Define the TRD geometry
  if ((frame->IsVersion() == 0) ||
      (frame->IsVersion() == 1)) {
    fGeometry = new AliTRDgeometry();
  }
  else {
    AliError("Could not find valid FRAME version\n");
    exit(1);
  }

  // Save the geometry
  TDirectory *saveDir = gDirectory;
  gAlice->GetRunLoader()->CdGAFile();
  fGeometry->Write("TRDgeometry");
  saveDir->cd();

  // Allocate the hit array
  fHits = new TClonesArray("AliTRDhit",405);
  gAlice->GetMCApp()->AddHitList(fHits);

  //PH SetMarkerColor(kWhite);   

}

//_____________________________________________________________________________
AliTRD::~AliTRD()
{
  //
  // TRD destructor
  //

  if (fGeometry) {
    delete fGeometry;
    fGeometry = 0;
  }

  if (fHits) {
    delete fHits;
    fHits     = 0;
  }

}

//_____________________________________________________________________________
void AliTRD::Hits2Digits()
{
  //
  // Create digits
  //

  AliTRDdigitizer digitizer("TRDdigitizer","TRD digitizer class");
  AliLog::SetClassDebugLevel("TRDdigitizer",AliDebugLevel());
  
  // Initialization
  digitizer.InitDetector();
    
  if (!fLoader->TreeH()) {
    fLoader->LoadHits("read");
  }
  fLoader->LoadDigits("recreate");

  AliRunLoader *runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    digitizer.Open(runLoader,iEvent);
    digitizer.MakeDigits();
    digitizer.WriteDigits();
  }

  fLoader->UnloadHits();
  fLoader->UnloadDigits();

}

//_____________________________________________________________________________
void AliTRD::Hits2SDigits()
{
  //
  // Create summable digits
  //

  AliTRDdigitizer digitizer("TRDdigitizer","TRD digitizer class");
  // For the summable digits
  digitizer.SetSDigits(kTRUE);
  AliLog::SetClassDebugLevel("TRDdigitizer",AliDebugLevel());

  // Initialization
  digitizer.InitDetector();
    
  if (!fLoader->TreeH()) {
    fLoader->LoadHits("read");
  }
  fLoader->LoadSDigits("recreate");

  AliRunLoader *runLoader = fLoader->GetRunLoader(); 

  for (Int_t iEvent = 0; iEvent < runLoader->GetNumberOfEvents(); iEvent++) {
    runLoader->GetEvent(iEvent);
    digitizer.Open(runLoader,iEvent);
    digitizer.MakeDigits();
    digitizer.WriteDigits();
  }

  fLoader->UnloadHits();
  fLoader->UnloadSDigits();
  
}

//_____________________________________________________________________________
AliDigitizer *AliTRD::CreateDigitizer(AliRunDigitizer *manager) const
{
  //
  // Creates a new digitizer object
  //

  return new AliTRDdigitizer(manager);

}

//_____________________________________________________________________________
void AliTRD::SDigits2Digits()
{
  //
  // Create final digits from summable digits
  //

  // Create the TRD digitizer
  AliTRDdigitizer digitizer("TRDdigitizer","TRD digitizer class");  
  AliLog::SetClassDebugLevel("TRDdigitizer",AliDebugLevel());

  // Set the parameter
  digitizer.SetEvent(gAlice->GetEvNumber());

  // Initialization
  digitizer.InitDetector();

  // Read the s-digits via digits manager
  AliTRDdigitsManager sdigitsManager;
 
  AliLog::SetClassDebugLevel("TRDdigitisManager",AliDebugLevel());
  sdigitsManager.SetSDigits(kTRUE);
  sdigitsManager.CreateArrays();
  
  if (!fLoader->TreeS()) { 
    if (fLoader->LoadSDigits("read")) {
      return;
    }
  }
  if (!fLoader->TreeS()) {
    AliError(Form("Error while reading SDigits for event %d",gAlice->GetEvNumber()));
    return;
  }
  
  sdigitsManager.ReadDigits(fLoader->TreeS());

  // Add the s-digits to the input list 
  digitizer.AddSDigitsManager(&sdigitsManager);

  // Convert the s-digits to normal digits
  digitizer.SDigits2Digits();

  // Store the digits
  if (!fLoader->TreeD()) {
    fLoader->MakeTree("D");
  }
  if (digitizer.MakeBranch(fLoader->TreeD())){
    digitizer.WriteDigits();
  }

}

//_____________________________________________________________________________
void AliTRD::Digits2Raw() 
{
  //
  // Convert digits of the current event to raw data
  //

  fLoader->LoadDigits();
  TTree *digits = fLoader->TreeD();
  if (!digits) {
    AliError("No digits tree");
    return;
  }

  AliTRDrawData rawWriter;
  if (!rawWriter.Digits2Raw(digits)) {
    AliError("The raw writer could not load the digits tree");
  }

  fLoader->UnloadDigits();

}

//_____________________________________________________________________________
void AliTRD::AddHit(Int_t track, Int_t det, Float_t *hits, Int_t q
                  , Bool_t inDrift)
{
  //
  // Add a hit for the TRD
  // 

  TClonesArray &lhits = *fHits;
  AliTRDhit *hit = new(lhits[fNhits++]) AliTRDhit(fIshunt,track,det,hits,q);

  if (inDrift) {
    hit->SetDrift();
  }
  else {
    hit->SetAmplification();
  }

  if (q < 0) {
    hit->SetTRphoton();
  }

}

//_____________________________________________________________________________
void AliTRD::BuildGeometry()
{
  //
  // Create the ROOT TNode geometry for the TRD
  //

  TNode *node;
  TNode *top;
  TPGON *pgon;

  // The dimensions of the TRD super module
  const Float_t kRmin  = 291.0;
  const Float_t kRmax  = 369.0;
  const Float_t kZmax1 = 378.35;
  const Float_t kZmax2 = 302.0;

  Float_t rmin;
  Float_t rmax;
  Float_t zmax1;
  Float_t zmax2;

  Int_t   iPlan;
 
  const Int_t kColorTRD = 46;
  
  // Find the top node alice
  top = gAlice->GetGeometry()->GetNode("alice");
  
  if      (fDisplayType == 0) {

    pgon = new TPGON("S_TRD","TRD","void",0,360,AliTRDgeometry::Nsect(),4);
    rmin = kRmin;
    rmax = kRmax;
    pgon->DefineSection(0,-kZmax1,rmax,rmax);
    pgon->DefineSection(1,-kZmax2,rmin,rmax);
    pgon->DefineSection(2, kZmax2,rmin,rmax);
    pgon->DefineSection(3, kZmax1,rmax,rmax);
    top->cd();
    node = new TNode("TRD","TRD","S_TRD",0,0,0,"");
    node->SetLineColor(kColorTRD);
    fNodes->Add(node);

  }
  else if (fDisplayType == 1) {

    Char_t name[7];

    Float_t slope = (kZmax1 - kZmax2) / (kRmax  - kRmin);

    rmin  = kRmin + AliTRDgeometry::CraHght();
    rmax  = rmin  + AliTRDgeometry::CdrHght();

    Float_t thickness = rmin - kRmin;
    zmax2 = kZmax2 + slope * thickness;
    zmax1 = zmax2 + slope * AliTRDgeometry::DrThick();

    for (iPlan = 0; iPlan < AliTRDgeometry::Nplan(); iPlan++) {

      sprintf(name,"S_TR1%d",iPlan);
      pgon  = new TPGON(name,"TRD","void",0,360,AliTRDgeometry::Nsect(),4);
      pgon->DefineSection(0,-zmax1,rmax,rmax);
      pgon->DefineSection(1,-zmax2,rmin,rmax);
      pgon->DefineSection(2, zmax2,rmin,rmax);
      pgon->DefineSection(3, zmax1,rmax,rmax);
      top->cd();
      node = new TNode("TRD","TRD",name,0,0,0,"");
      node->SetLineColor(kColorTRD);
      fNodes->Add(node);

      Float_t height = AliTRDgeometry::Cheight() + AliTRDgeometry::Cspace(); 
      rmin  = rmin  + height;
      rmax  = rmax  + height;
      zmax1 = zmax1 + slope * height;
      zmax2 = zmax2 + slope * height;

    }

    thickness += AliTRDgeometry::DrThick();
    rmin       = kRmin  + thickness;
    rmax       = rmin   + AliTRDgeometry::AmThick();
    zmax2      = kZmax2 + slope * thickness;
    zmax1      = zmax2  + slope * AliTRDgeometry::AmThick();

    for (iPlan = 0; iPlan < AliTRDgeometry::Nplan(); iPlan++) {

      sprintf(name,"S_TR2%d",iPlan);
      pgon  = new TPGON(name,"TRD","void",0,360,AliTRDgeometry::Nsect(),4);
      pgon->DefineSection(0,-zmax1,rmax,rmax);
      pgon->DefineSection(1,-zmax2,rmin,rmax);
      pgon->DefineSection(2, zmax2,rmin,rmax);
      pgon->DefineSection(3, zmax1,rmax,rmax);
      top->cd();
      node = new TNode("TRD","TRD",name,0,0,0,"");
      node->SetLineColor(kColorTRD);
      fNodes->Add(node);

      Float_t height = AliTRDgeometry::Cheight() + AliTRDgeometry::Cspace(); 
      rmin  = rmin  + height;
      rmax  = rmax  + height;
      zmax1 = zmax1 + slope * height;
      zmax2 = zmax2 + slope * height;

    }

  }

}
 
//_____________________________________________________________________________
void AliTRD::CreateGeometry()
{
  //
  // Creates the volumes for the TRD chambers
  //

  // Check that FRAME is there otherwise we have no place where to put the TRD
  AliModule *frame = gAlice->GetModule("FRAME");
  if (!frame) {
    AliFatal("The TRD needs the FRAME to be defined first");
  }

  fGeometry->CreateGeometry(fIdtmed->GetArray() - 1299);

}
 
//_____________________________________________________________________________
void AliTRD::CreateMaterials()
{
  //
  // Create the materials for the TRD
  //

  Int_t   isxfld = gAlice->Field()->Integ();
  Float_t sxmgmx = gAlice->Field()->Max();
  
  // For polyethilene (CH2) 
  Float_t ape[2]    = { 12.011 ,  1.0079 };
  Float_t zpe[2]    = {  6.0   ,  1.0    };
  Float_t wpe[2]    = {  1.0   ,  2.0    };
  Float_t dpe       = 0.95;

  // For CO2 
  Float_t aco[2]    = { 12.011 , 15.9994 };
  Float_t zco[2]    = {  6.0   ,  8.0    };
  Float_t wco[2]    = {  1.0   ,  2.0    };
  Float_t dco       = 0.00186;

  // For water
  Float_t awa[2]    = {  1.0079, 15.9994 };
  Float_t zwa[2]    = {  1.0   ,  8.0    };
  Float_t wwa[2]    = {  2.0   ,  1.0    };
  Float_t dwa       = 1.0;

  // For isobutane (C4H10)
  Float_t ais[2]    = { 12.011 ,  1.0079 };
  Float_t zis[2]    = {  6.0   ,  1.0    };
  Float_t wis[2]    = {  4.0   , 10.0    };
  Float_t dis       = 0.00267;

  // For plexiglas (C5H8O2)
  Float_t apg[3]    = { 12.011 ,  1.0079, 15.9994 };
  Float_t zpg[3]    = {  6.0   ,  1.0   ,  8.0    };
  Float_t wpg[3]    = {  5.0   ,  8.0   ,  2.0    };
  Float_t dpg       = 1.18; 
  
  // For epoxy (C18H19O3)
  Float_t aEpoxy[3] = { 15.9994,  1.0079, 12.011  }; 
  Float_t zEpoxy[3] = {  8.0   ,  1.0   ,  6.0    }; 
  Float_t wEpoxy[3] = {  3.0   , 19.0   , 18.0    }; 
  Float_t dEpoxy    = 1.8 ; 

  // For Araldite, low density epoxy (C18H19O3)
  Float_t aAral[3]  = { 15.9994,  1.0079, 12.011  }; 
  Float_t zAral[3]  = {  8.0   ,  1.0   ,  6.0    }; 
  Float_t wAral[3]  = {  3.0   , 19.0   , 18.0    }; 
  Float_t dAral     = 1.05; 

  // For air  
  Float_t aAir[4]   = { 12.011   , 14.0     , 15.9994  , 36.0      };
  Float_t zAir[4]   = {  6.0     ,  7.0     ,  8.0     , 18.0      };
  Float_t wAir[4]   = {  0.000124,  0.755267,  0.231781,  0.012827 };
  Float_t dAir      = 1.20479e-03;

  // For G10
  Float_t aG10[4]   = {  1.0079  , 12.011   , 15.9994  , 28.086    };
  Float_t zG10[4]   = {  1.0     ,  6.0     ,  8.0     , 14.0      };
  Float_t wG10[4]   = {  0.15201 ,  0.10641 ,  0.49444 ,  0.24714  };
  Float_t dG10      = 1.7;

  // For Xe/CO2-gas-mixture 
  Float_t aXeCO2[3] = { 131.29   ,  12.0107 ,  15.9994  };
  Float_t zXeCO2[3] = {  54.0    ,   6.0    ,   8.0     };
  Float_t wXeCO2[3] = {   8.5    ,   1.5    ,   3.0     }; 
  // Xe-content of the Xe/CO2-mixture (85% / 15%) 
  Float_t fxc       = 0.85;
  Float_t dxe       = 0.00549;
  Float_t dgm       = fxc * dxe + (1.0 - fxc) * dco;
  
  // General tracking parameter
  Float_t tmaxfd    = -10.0;
  Float_t stemax    = -1.0e10;
  Float_t deemax    = -0.1;
  Float_t epsil     =  1.0e-4;
  Float_t stmin     = -0.001;
  
  //////////////////////////////////////////////////////////////////////////
  //     Define Materials 
  //////////////////////////////////////////////////////////////////////////

  AliMaterial( 1, "Al"   ,  26.98, 13.0, 2.7     ,     8.9 ,    37.2);
  AliMaterial( 4, "Xe"   , 131.29, 54.0, dxe     ,  1546.16,     0.0);
  AliMaterial( 5, "Cu"   ,  63.54, 29.0, 8.96    ,     1.43,    14.8);
  AliMaterial( 6, "C"    ,  12.01,  6.0, 2.265   ,    18.8 ,    74.4);
  AliMaterial(15, "Sn"   , 118.71, 50.0, 7.31    ,     1.21,    14.8);
  AliMaterial(16, "Si"   ,  28.09, 14.0, 2.33    ,     9.36,    37.2);
  AliMaterial(18, "Fe"   ,  55.85, 26.0, 7.87    ,     1.76,    14.8);

  // Mixtures 
  AliMixture(2, "Air"         , aAir,   zAir,   dAir,    4, wAir  );
  AliMixture(3, "Polyethilene", ape,    zpe,    dpe,    -2, wpe   );
  AliMixture(8, "CO2",          aco,    zco,    dco,    -2, wco   );
  AliMixture(9, "Isobutane",    ais,    zis,    dis,    -2, wis   );
  AliMixture(10,"Gas mixture",  aXeCO2, zXeCO2, dgm,    -3, wXeCO2);
  AliMixture(12,"G10",          aG10,   zG10,   dG10,    4, wG10  );
  AliMixture(13,"Water",        awa,    zwa,    dwa,    -2, wwa   );
  AliMixture(14,"Plexiglas",    apg,    zpg,    dpg,    -3, wpg   );
  AliMixture(17,"Epoxy",        aEpoxy, zEpoxy, dEpoxy, -3, wEpoxy);
  AliMixture(19,"Araldite",     aAral,  zAral,  dAral,  -3, wAral );

  //////////////////////////////////////////////////////////////////////////
  //     Tracking Media Parameters 
  //////////////////////////////////////////////////////////////////////////

  // Al Frame 
  AliMedium( 1,"Al Frame"   , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Air 
  AliMedium( 2,"Air"        , 2,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Wires
  AliMedium( 3,"Wires"      , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // All other ROB materials (caps, etc.)
  AliMedium( 4,"ROB Other"  , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu pads 
  AliMedium( 5,"Padplane"   , 5,1,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Fee + cables 
  AliMedium( 6,"Readout"    , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // C frame 
  AliMedium( 7,"C Frame"    , 6,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // INOX of cooling bus bars
  AliMedium( 8,"Cooling bus",18,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Gas-mixture (Xe/CO2) 
  AliMedium( 9,"Gas-mix"    ,10,1,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Nomex-honeycomb
  AliMedium(10,"Nomex"      ,12,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Araldite glue
  AliMedium(11,"Glue"       ,19,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // G10-plates
  AliMedium(13,"G10-plates" ,12,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cooling water
  AliMedium(14,"Water"      ,13,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Rohacell (plexiglas) for the radiator
  AliMedium(15,"Rohacell"   ,14,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Al layer in MCMs
  AliMedium(16,"MCM-Al"     , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Sn layer in MCMs
  AliMedium(17,"MCM-Sn"     ,15,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu layer in MCMs
  AliMedium(18,"MCM-Cu"     , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // G10 layer in MCMs
  AliMedium(19,"MCM-G10"    ,12,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Si in readout chips
  AliMedium(20,"Chip-Si"    ,16,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Epoxy in readout chips
  AliMedium(21,"Chip-Ep"    ,17,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // PE in connectors
  AliMedium(22,"Conn-PE"    , 3,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu in connectors
  AliMedium(23,"Chip-Cu"    , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Al of cooling pipes
  AliMedium(24,"Cooling"    , 1,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);
  // Cu in services
  AliMedium(25,"Serv-Cu"    , 5,0,isxfld,sxmgmx
              ,tmaxfd,stemax,deemax,epsil,stmin);

  // Save the density values for the TRD absorbtion
  Float_t dmy  = 1.39;
  fFoilDensity = dmy;
  fGasDensity  = dgm;

}

//_____________________________________________________________________________
void AliTRD::DrawModule() const
{
  //
  // Draw a shaded view of the Transition Radiation Detector version 0
  //

  // Set everything unseen
  gMC->Gsatt("*"   ,"SEEN",-1);
  
  // Set ALIC mother transparent
  gMC->Gsatt("ALIC","SEEN", 0);
  
  // Set the volumes visible
  if (fGeometry->IsVersion() == 0) {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("UTR1","SEEN", 0);
    gMC->Gsatt("UTR2","SEEN", 0);
    gMC->Gsatt("UTR3","SEEN", 0);
  }
  else {
    gMC->Gsatt("B071","SEEN", 0);
    gMC->Gsatt("B074","SEEN", 0);
    gMC->Gsatt("B075","SEEN", 0);
    gMC->Gsatt("B077","SEEN", 0);
    gMC->Gsatt("BTR1","SEEN", 0);
    gMC->Gsatt("BTR2","SEEN", 0);
    gMC->Gsatt("BTR3","SEEN", 0);
    gMC->Gsatt("UTR1","SEEN", 0);
  }
  
  gMC->Gdopt("hide", "on");
  gMC->Gdopt("shad", "on");
  gMC->Gsatt("*", "fill", 7);
  gMC->SetClipBox(".");
  gMC->SetClipBox("*", 0, 2000, -2000, 2000, -2000, 2000);
  gMC->DefaultRange();
  gMC->Gdraw("alic", 40, 30, 0, 12, 9.4, .021, .021);
  gMC->Gdhead(1111, "Transition Radiation Detector");
  gMC->Gdman(18, 4, "MAN");

}

//_____________________________________________________________________________
Int_t AliTRD::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Distance between the mouse and the TRD detector on the screen
  // Dummy routine
  //
  
  return 9999;

}
 
//_____________________________________________________________________________
void AliTRD::Init()
{
  //
  // Initialize the TRD detector after the geometry has been created
  //

  AliDebug(1,"++++++++++++++++++++++++++++++++++++++++++++++");

  if (fGeometry->IsVersion() != 1) {
    AliError("Not a valid geometry");
  }
  // Special tracking options for charged particles for XeCO2
  gMC->Gstpar((* fIdtmed)[9],"DRAY",1.0);
  gMC->Gstpar((* fIdtmed)[9],"STRA",1.0); 
  
}

//_____________________________________________________________________________
void AliTRD::LoadPoints(Int_t )
{
  //
  // Store x, y, z of all hits in memory.
  // Hit originating from TR photons are given a different color
  //

  if (fHits == 0) {
    return;
  }

  Int_t nhits  = fHits->GetEntriesFast();
  if (nhits == 0) {
    return;
  }

  Int_t tracks = gAlice->GetMCApp()->GetNtrack();
  if (fPoints == 0) {
    fPoints = new TObjArray(tracks);
  }

  AliTRDhit *ahit;
  
  Int_t    *ntrkE = new Int_t[tracks];
  Int_t    *ntrkT = new Int_t[tracks];
  Int_t    *limiE = new Int_t[tracks];
  Int_t    *limiT = new Int_t[tracks];
  Float_t **coorE = new Float_t*[tracks];
  Float_t **coorT = new Float_t*[tracks];
  for(Int_t i = 0; i < tracks; i++) {
    ntrkE[i] = 0;
    ntrkT[i] = 0;
    coorE[i] = 0;
    coorT[i] = 0;
    limiE[i] = 0;
    limiT[i] = 0;
  }
  
  AliTRDpoints *points = 0;
  Float_t      *fp     = 0;
  Int_t         trk;
  Int_t         chunk  = nhits / 4 + 1;

  // Loop over all the hits and store their position
  ahit = (AliTRDhit *) FirstHit(-1);
  while (ahit) {

    // dEdx hits
    if (ahit->GetCharge() >= 0) {

      trk = ahit->GetTrack();
      if (ntrkE[trk] == limiE[trk]) {
        // Initialise a new track
        fp = new Float_t[3*(limiE[trk]+chunk)];
        if (coorE[trk]) {
          memcpy(fp,coorE[trk],sizeof(Float_t)*3*limiE[trk]);
          delete [] coorE[trk];
        }
        limiE[trk] += chunk;
        coorE[trk]  = fp;
      } 
      else {
        fp = coorE[trk];
      }
      fp[3*ntrkE[trk]  ] = ahit->X();
      fp[3*ntrkE[trk]+1] = ahit->Y();
      fp[3*ntrkE[trk]+2] = ahit->Z();
      ntrkE[trk]++;

    }
    // TR photon hits
    else if ((ahit->GetCharge() < 0) && 
             (fDrawTR)) {

      trk = ahit->GetTrack();
      if (ntrkT[trk] == limiT[trk]) {
        // Initialise a new track
        fp = new Float_t[3*(limiT[trk]+chunk)];
        if (coorT[trk]) {
          memcpy(fp,coorT[trk],sizeof(Float_t)*3*limiT[trk]);
          delete [] coorT[trk];
        }
        limiT[trk] += chunk;
        coorT[trk]  = fp;
      } 
      else {
        fp = coorT[trk];
      }
      fp[3*ntrkT[trk]  ] = ahit->X();
      fp[3*ntrkT[trk]+1] = ahit->Y();
      fp[3*ntrkT[trk]+2] = ahit->Z();
      ntrkT[trk]++;

    }

    ahit = (AliTRDhit *) NextHit();

  }

  for (trk = 0; trk < tracks; ++trk) {

    if (ntrkE[trk] || ntrkT[trk]) {

      points = new AliTRDpoints();
      points->SetDetector(this);
      points->SetParticle(trk);

      // Set the dEdx points
      if (ntrkE[trk]) {
        points->SetMarkerColor(kWhite); //PH This is the default color in TRD
        points->SetMarkerSize(1); //PH Default size=1
        points->SetPolyMarker(ntrkE[trk],coorE[trk],1); //PH Default style=1
        delete [] coorE[trk];
        coorE[trk] = 0;
      }

      // Set the TR photon points
      if (ntrkT[trk]) {
        points->SetTRpoints(ntrkT[trk],coorT[trk]);
        delete [] coorT[trk];
        coorT[trk] = 0;
      }

      fPoints->AddAt(points,trk);

    }

  }

  delete [] coorE;
  delete [] coorT;
  delete [] ntrkE;
  delete [] ntrkT;
  delete [] limiE;
  delete [] limiT;

}

//_____________________________________________________________________________
void AliTRD::MakeBranch(Option_t *option)
{
  //
  // Create Tree branches for the TRD digits.
  //

  Int_t  buffersize = 4000;
  Char_t branchname[15];
  sprintf(branchname,"%s",GetName());

  const Char_t *cD = strstr(option,"D");

  AliDetector::MakeBranch(option);

  if (fDigits         && 
      gAlice->TreeD() && 
      cD) {
    MakeBranchInTree(gAlice->TreeD(),branchname,&fDigits,buffersize,0);
  }

}

//_____________________________________________________________________________
void AliTRD::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  //

  fNdigits = 0;

  if (fDigits) {
    fDigits->Clear();
  }

}

//_____________________________________________________________________________
void AliTRD::SetTreeAddress()
{
  //
  // Set the branch addresses for the trees.
  //

  if (fLoader->TreeH() && 
      (fHits == 0x0)) {
    fHits = new TClonesArray("AliTRDhit",405);
  }
  AliDetector::SetTreeAddress();

}

//_____________________________________________________________________________
AliTRD &AliTRD::operator=(const AliTRD &trd)
{
  //
  // Assignment operator
  //

  if (this != &trd) {
    ((AliTRD &) trd).Copy(*this);
  }

  return *this;

} 















































































































































































































































































































































