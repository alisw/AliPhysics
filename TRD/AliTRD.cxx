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

#include <TVirtualMC.h>
 
#include "AliMC.h"
#include "AliMagF.h"
#include "AliRun.h"

#include "AliTRD.h"
#include "AliTRDdigitizer.h"
#include "AliTRDdigitsManager.h"
#include "AliTRDgeometry.h"
#include "AliTRDhit.h"
#include "AliTRDrawData.h"
#include "AliTRDSimParam.h"

ClassImp(AliTRD)
 
//_____________________________________________________________________________
AliTRD::AliTRD()
  :AliDetector()
  ,fGeometry(0)
  ,fGasDensity(0)
  ,fFoilDensity(0)
  ,fGasNobleFraction(0)
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
  ,fGasNobleFraction(0)
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

  // Allocate the hit array
  fHits = new TClonesArray("AliTRDhit",405);
  gAlice->GetMCApp()->AddHitList(fHits);

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
                  , Float_t time, Bool_t inDrift)
{
  //
  // Add a hit for the TRD
  // 

  TClonesArray &lhits = *fHits;
  AliTRDhit *hit = new(lhits[fNhits++]) AliTRDhit(fIshunt
                                                 ,track
                                                 ,det
                                                 ,hits
                                                 ,q
                                                 ,time);

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
  Float_t dco       = 0.00186; // at 20C

  // For water
  Float_t awa[2]    = {  1.0079, 15.9994 };
  Float_t zwa[2]    = {  1.0   ,  8.0    };
  Float_t wwa[2]    = {  2.0   ,  1.0    };
  Float_t dwa       = 1.0;

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
  Float_t dxe       = 0.00549; // at 20C
  Float_t dgmXe     = fxc * dxe + (1.0 - fxc) * dco;

  // For Ar/CO2-gas-mixture
  Float_t aArCO2[3] = {  39.948  ,  12.0107 ,  15.9994  };
  Float_t zArCO2[3] = {  18.0    ,   6.0    ,   8.0     };
  Float_t wArCO2[3] = {   8.2    ,   1.8    ,   3.6     }; 
  // Ar-content of the Ar/CO2-mixture (82% / 18%) 
  Float_t fac       = 0.82;
  Float_t dar       = 0.00166; // at 20C
  Float_t dgmAr     = fac * dar + (1.0 - fac) * dco;
  
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
  if      (AliTRDSimParam::Instance()->IsXenon()) {
    AliMixture(10,"XeCO2",        aXeCO2, zXeCO2, dgmXe,  -3, wXeCO2);
  }
  else if (AliTRDSimParam::Instance()->IsArgon()) {
    AliInfo("Gas mixture: Ar C02 (80/20)");
    AliMixture(10,"ArCO2",        aArCO2, zArCO2, dgmAr,  -3, wArCO2);
  }
  else {
    AliFatal("Wrong gas mixture");
    exit(1);
  }
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
  if      (AliTRDSimParam::Instance()->IsXenon()) {
    fGasDensity       = dgmXe;
    fGasNobleFraction = fxc;
  }
  else if (AliTRDSimParam::Instance()->IsArgon()) {
    fGasDensity       = dgmAr;
    fGasNobleFraction = fac;
  }

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
  gMC->Gstpar((* fIdtmed)[9],"DRAY"    , 1.0);
  gMC->Gstpar((* fIdtmed)[9],"STRA"    , 1.0); 
  gMC->Gstpar((* fIdtmed)[9],"LOSS"    ,13.0);      // Specific energy loss
  gMC->Gstpar((* fIdtmed)[9],"PRIMIO_E",23.53);     // 1st ionisation potential
  gMC->Gstpar((* fIdtmed)[9],"PRIMIO_N",19.344431); // Number of primaries

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
Bool_t AliTRD::Raw2SDigits(AliRawReader *rawReader)
{
  //
  // Converts RAW data to SDigits
  //

  AliLoader *loader = fRunLoader->GetLoader("TRDLoader");
  if (!loader) {
    AliError("Can not get TRD loader from Run Loader");
    return kFALSE;
  }
    
  TTree *tree = 0;
  tree = loader->TreeS();
  if (!tree) {
    loader->MakeTree("S");
    tree = loader->TreeS();
  }

  AliTRDrawData       *rawdata        = new AliTRDrawData();
  AliTRDdigitsManager *sdigitsManager = rawdata->Raw2Digits(rawReader);
  if (sdigitsManager) {
    sdigitsManager->SetSDigits(kTRUE);
    sdigitsManager->MakeBranch(tree);
    sdigitsManager->WriteDigits();
    return kTRUE;
  } 
  else {
    return kFALSE;
  }

}

//_____________________________________________________________________________
AliLoader* AliTRD::MakeLoader(const char* topfoldername)
{
 fLoader = new AliLoader(GetName(),topfoldername);

 AliInfo("Adding Tracklets-loader");
 AliDataLoader *dl = new AliDataLoader("TRD.Tracklets.root","tracklets", "tracklets");
 fLoader->AddDataLoader(dl);

 return fLoader;
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
