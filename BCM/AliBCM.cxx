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
//  Beam Condition Monitor BCM                                               //
//                                                                           //
//  andreas.morsch@cern.ch                                                   // 
///////////////////////////////////////////////////////////////////////////////
 
#include <TClonesArray.h>
#include <TGeoCompositeShape.h>
#include <TGeoGlobalMagField.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoPgon.h>
#include <TGeoVolume.h>
#include <TGeoXtru.h>
#include <TVirtualMC.h>

#include "AliBCM.h"
#include "AliBCMHit.h"
#include "AliBCMLoader.h"
#include "AliMagF.h"
#include "AliRun.h"
#include "AliMC.h"
#include "AliLog.h"
 
ClassImp(AliBCM)

 
//_____________________________________________________________________________
AliBCM::AliBCM():
    AliDetector(),
    fVolId(0)
{
  //
  // Default constructor for BCM
  //
}
 
//_____________________________________________________________________________
AliBCM::AliBCM(const char *name, const char *title): 
    AliDetector(name,title),
    fVolId(0)  
{
//
//  Constructor
    fHits  = new TClonesArray("AliBCMHit");
    fNhits = 0;    
    gAlice->GetMCApp()->AddHitList(fHits);
}
AliBCM::~AliBCM()
{  
    // Destructor
    if(fHits)      delete fHits;
}

void AliBCM::StepManager()
{
//
// Step Manager for ALICE Beam Condition Monitor
//    

    static Float_t edepT;    
    static Double_t xh[4] = {0., 0., 0., 0.};
    Float_t edep = 0.;
    Int_t   copy = -1; 
    
    
    if (fMC->TrackCharge() &&
    fMC->CurrentVolID(copy) == fVolId) {
	// Charged particle inside sensitive volume
	//
	// Entering
    if (fMC->IsTrackEntering()) {
	    edepT = 0.;
        fMC->TrackPosition(xh[0],xh[1],xh[2]);
        xh[3] = fMC->TrackTime();
	}
	
	//
	// Any step
    if ((edep = fMC->Edep()) > 0.) {
	    Double_t x[3];   
        fMC->TrackPosition(x[0],x[1],x[2]);
	    edepT += edep;
	}
	//
	// Exiting 
    if(fMC->IsTrackExiting()||fMC->IsTrackStop()||fMC->IsTrackDisappeared())
	{
	    Int_t track = gAlice->GetMCApp()->GetCurrentTrackNumber();
	    TClonesArray &lhits = *fHits;
	    Int_t ic = copy + 10;
	    if (xh[2] < 0.) ic+=10;
	    new(lhits[fNhits++]) AliBCMHit(1, track, xh, ic, edepT);
	}
    }

    
}

//_____________________________________________________________________________
void AliBCM::CreateGeometry()
{
    //
    // Create geometry for BCM
    //
    
    //
    // Top volume 
    TGeoVolume* top      = gGeoManager->GetVolume("ALIC");
    // Media 
    TGeoMedium* medPCD   = gGeoManager->GetMedium("BCM_PCD");
    // Rotations
    TGeoRotation* rotxz   = new TGeoRotation("rotxz" ,  270.,   0., 90.,  90., 180., 0.);
    TGeoRotation* rot000  = new TGeoRotation("rot000",   90.,   0., 90.,  90.,   0., 0.);
    TGeoRotation* rot090  = new TGeoRotation("rot090",   90.,  90., 90., 180.,   0., 0.);
    TGeoRotation* rot180  = new TGeoRotation("rot180",   90., 180., 90., 270.,   0., 0.);
    TGeoRotation* rot270  = new TGeoRotation("rot270",   90., 270., 90.,   0.,   0., 0.);
    //
    const Float_t kWidth     = 1.00;
    const Float_t kThickness = 0.05;
    const Float_t rBCM       = 7.;
    
    
    TGeoBBox*   shBCMpcd = new TGeoBBox(kWidth/2., kWidth/2., kThickness/2.);
    TGeoVolume* voBCMpcd = new TGeoVolume("BCMpcd", shBCMpcd, medPCD);
    TGeoVolumeAssembly* voBCM = new TGeoVolumeAssembly("BCM");
    
    voBCM->AddNode(voBCMpcd, 1, new TGeoCombiTrans(+rBCM, 0.    , 0., rot000));
    voBCM->AddNode(voBCMpcd, 2, new TGeoCombiTrans(0.   , +rBCM , 0., rot090));
    voBCM->AddNode(voBCMpcd, 3, new TGeoCombiTrans(-rBCM, 0.,     0., rot180));
    voBCM->AddNode(voBCMpcd, 4, new TGeoCombiTrans(0.   , -rBCM , 0., rot270));
    
    top->AddNode(voBCM, 1, new TGeoTranslation(0., 0., 1561.));
    top->AddNode(voBCM, 2, new TGeoCombiTrans(0., 0., -1908., rotxz));
    
}

//_____________________________________________________________________________
void AliBCM::CreateMaterials()
{
  //
  // Create materials for BCM
  //
  // Polycristalline Diamond
    Float_t rho  = 3.53;
    Float_t absl = 86.3 / rho;
    Float_t radl = 42.7 / rho;
    //
    Float_t epsil  = .001;   // Tracking precision, 
    Float_t stemax = -1.;    // Maximum displacement for multiple scat 
    Float_t tmaxfd = -20. ;  // Maximum angle due to field deflection 
    Float_t deemax = -.01;   // Maximum fractional energy loss, DLS 
    Float_t stmin  = -.8;
    Int_t   isxfld = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Integ();
    Float_t sxmgmx = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->Max();

    AliMaterial(1, "PCD", 12.011, 6., rho, radl, absl);
    //
    AliMedium(1, "PCD", 1, 0, isxfld, sxmgmx, tmaxfd, stemax, deemax, epsil, stmin);
}

//_____________________________________________________________________________
void AliBCM::Init()
{
    //
    // Initialise BCM magnet after it has been built
    //
  Int_t i;
  //
  if(AliLog::GetGlobalDebugLevel()>0) {
    printf("\n%s: ",ClassName());
    for(i=0;i<35;i++) printf("*");
    printf(" BCM_INIT ");
    for(i=0;i<35;i++) printf("*");
    printf("\n%s: ",ClassName());
    for(i=0;i<80;i++) printf("*");
    printf("\n");
  }

    //
    // Here the BCM initialisation code (if any!)
    fVolId =  TVirtualMC::GetMC()->VolId("BCMpcd");
}

void AliBCM::MakeBranch(Option_t* option)
{
//Create Tree branches for the BCM
    
  const Int_t kBufSize = 4000;
  const char *cH = strstr(option,"H");
//  const char *cD = strstr(option,"D");
//  const char *cR = strstr(option,"R");
//  const char *cS = strstr(option,"S");

  if(cH && fLoader->TreeH() && (fHits == 0x0)){
      fHits  = new TClonesArray("AliBCMHit");
      fNhits = 0;    
      MakeBranchInTree(fLoader->TreeH(), "BCM" ,&fHits ,kBufSize, 0);
  }
  AliDetector::MakeBranch(option);
}

void AliBCM::SetTreeAddress()
{
  // Set branch address

    if (fLoader->TreeH() && fHits==0x0)
	fHits   = new TClonesArray("AliBCMHit",  4000);
    AliDetector::SetTreeAddress();
}

//_____________________________________________________________________________
AliLoader* AliBCM::MakeLoader(const char* topfoldername)
{ 
  //
  // Builds BCM getter (AliLoader type)
  AliDebug(1,Form("Creating AliBCMLoader, Top folder is %s ",topfoldername));
  fLoader = new AliBCMLoader(GetName(),topfoldername);
  return fLoader;
}
