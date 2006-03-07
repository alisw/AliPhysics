/* $Id$*/

// Script to create alignment parameters and store them into CDB
// Three sets of alignment parameters can be created:
// 0) 1 PHOS module with ideal geometry
// 1) 5 PHOS modules with small disalignments
// 2) 5 PHOS modules with ideal geometry

#if !defined(__CINT__)
#include "TControlBar.h"
#include "TString.h"
#include "TRandom.h"

#include "AliRun.h"
#include "AliPHOSAlignData.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#endif


void AliPHOSSetAlignment()
{
  TControlBar *menu = new TControlBar("vertical","PHOS alignment control");
  menu->AddButton("Help to run PHOS alignment control","Help()",
		  "Explains how to use PHOS alignment control menus");

  menu->AddButton("PHOS 2007","SetAlignment(0)",
		  "Set PHOS alignment for the LHC run 2007");
  menu->AddButton("Full misaligned PHOS","SetAlignment(1)",
		  "Set all 5 modules with random displacement");
  menu->AddButton("Full ideal PHOS","SetAlignment(2)",
		  "Set all 5 modules with random displacement");

  menu->AddButton("Read PHOS 2007","GetAlignment(0)",
		  "Read PHOS geometry for the LHC run 2007");
  menu->AddButton("Read full misaligned PHOS","GetAlignment(1)",
		  "Read full PHOS geometry with random displacements");
  menu->AddButton("Read full ideal PHOS","GetAlignment(2)",
		  "Read full PHOS geometry with random displacements");

  menu->Show();
}

//------------------------------------------------------------------------
void Help()
{
  char *string =
    "\nSet PHOS alignment parameters and write them into ALICE CDB.
Press button \"PHOS 2007\" to create PHOS geometry for the first LHC run in 2007.
Press button \"Full PHOS\" to create full PHOS geometry with radnomly displaced modules\n";
  printf(string);
}

//------------------------------------------------------------------------
void SetAlignment(Int_t flag=0)
{
  // Write alignment parameters into CDB
  // Arguments:
  //   flag=0: ideal geometry with one module
  //   flag=1: disligned geometry with 5 modules
  // Author: Yuri Kharlov

  TString DBFolder;
  Int_t firstRun   =  0;
  Int_t lastRun    = 10;
  Int_t beamPeriod =  1;
  char* objFormat  = "";
  gRandom->SetSeed(0);

  AliPHOSAlignData *alignda=new AliPHOSAlignData("PHOS");
  
  switch (flag) {
  case 0:
    DBFolder  ="local://InitAlignDB";
    firstRun  =  0;
    lastRun   =  0;
    objFormat = "PHOS ideal geometry with 1 module";

    alignda->SetNModules(1);
    
    alignda->SetModuleCenter(0,0,   0.);
    alignda->SetModuleCenter(0,1,-460.);
    alignda->SetModuleCenter(0,2,   0.);
    
    alignda->SetModuleAngle(0,0,0,  90.);
    alignda->SetModuleAngle(0,0,1,   0.);
    alignda->SetModuleAngle(0,1,0,   0.);
    alignda->SetModuleAngle(0,1,1,   0.);
    alignda->SetModuleAngle(0,2,0,  90.);
    alignda->SetModuleAngle(0,2,1, 270.);
    break;
  case 1:
    DBFolder  ="local://DisAlignDB";
    firstRun  =  0;
    lastRun   = 10;
    objFormat = "PHOS disaligned geometry with 5 modules";

    Int_t nModules = 5;
    alignda->SetNModules(nModules);
    
    Float_t dAngle= 20;
    Float_t r0    = 460.;
    Float_t theta, phi;
    for (Int_t iModule=0; iModule<nModules; iModule++) {
      Float_t r = r0 + gRandom->Uniform(0.,40.);
      Float_t angle = dAngle * ( iModule - nModules / 2.0 + 0.5 ) ;
      angle += gRandom->Uniform(-0.1,+0.1);
      Float_t x = r * TMath::Sin(angle * TMath::DegToRad() );
      Float_t y =-r * TMath::Cos(angle * TMath::DegToRad() );
      Float_t z = gRandom->Uniform(-0.05,0.05);

      alignda->SetModuleCenter(iModule,0,x);
      alignda->SetModuleCenter(iModule,1,y);
      alignda->SetModuleCenter(iModule,2,z);
    
      theta = 90 + gRandom->Uniform(-1.,1.);
      phi   = angle;
      alignda->SetModuleAngle(iModule,0,0,theta);
      alignda->SetModuleAngle(iModule,0,1,phi);
      theta = 0;
      phi   = 0;
      alignda->SetModuleAngle(iModule,1,0,theta);
      alignda->SetModuleAngle(iModule,1,1,phi);
      theta = 90 + gRandom->Uniform(-1.,1.);
      phi   = 270+angle;
      alignda->SetModuleAngle(iModule,2,0,theta);
      alignda->SetModuleAngle(iModule,2,1,phi);
    }
  case 2:
    DBFolder  ="local://AlignDB";
    firstRun  =  0;
    lastRun   = 10;
    objFormat = "PHOS disaligned geometry with 5 modules";

    Int_t nModules = 5;
    alignda->SetNModules(nModules);
    
    Float_t dAngle= 20;
    Float_t r0    = 460.;
    Float_t theta, phi;
    for (Int_t iModule=0; iModule<nModules; iModule++) {
      Float_t r = r0;
      Float_t angle = dAngle * ( iModule - nModules / 2.0 + 0.5 ) ;
      Float_t x = r * TMath::Sin(angle * TMath::DegToRad() );
      Float_t y =-r * TMath::Cos(angle * TMath::DegToRad() );
      Float_t z = 0.;

      alignda->SetModuleCenter(iModule,0,x);
      alignda->SetModuleCenter(iModule,1,y);
      alignda->SetModuleCenter(iModule,2,z);
    
      theta = 90;
      phi   = angle;
      alignda->SetModuleAngle(iModule,0,0,theta);
      alignda->SetModuleAngle(iModule,0,1,phi);
      theta = 0;
      phi   = 0;
      alignda->SetModuleAngle(iModule,1,0,theta);
      alignda->SetModuleAngle(iModule,1,1,phi);
      theta = 90;
      phi   = 270+angle;
      alignda->SetModuleAngle(iModule,2,0,theta);
      alignda->SetModuleAngle(iModule,2,1,phi);
    }
    break;
  default:
    printf("Unknown flag %d, can be 0 or 1 only\n",flag);
    return;
  }
  

  //Store calibration data into database
  
  AliCDBMetaData md;
  md.SetComment(objFormat);
  md.SetBeamPeriod(beamPeriod);
  md.SetResponsible("Yuri Kharlov");
  
  AliCDBId id("PHOS/Alignment/Geometry",firstRun,lastRun);

  AliCDBManager* man = AliCDBManager::Instance();  
  AliCDBStorage* loc = man->GetStorage(DBFolder.Data());
  loc->Put(alignda, id, &md);

}

//------------------------------------------------------------------------
void GetAlignment(Int_t flag=0)
{
  // Read alignment parameters into CDB
  // Arguments:
  //   flag=0: ideal geometry with one module
  //   flag=1: disligned geometry with 5 modules
  // Author: Yuri Kharlov

  TString DBFolder;
  Int_t run   =  0;
  AliPHOSAlignData *alignda=new AliPHOSAlignData("PHOS");
  
  switch (flag) {
  case 0:
    DBFolder  ="local://InitAlignDB";
    break;
  case 1:
    DBFolder  ="local://DisAlignDB";
    break;
  case 2:
    DBFolder  ="local://AlignDB";
    break;
  default:
    printf("Unknown flag %d, can be 0 or 1 only\n",flag);
    return;
  }
  AliPHOSAlignData* alignda  = (AliPHOSAlignData*)
    (AliCDBManager::Instance()
     ->GetStorage(DBFolder.Data())
     ->Get("PHOS/Alignment/Geometry",run)->GetObject());
  alignda->Print();
}
