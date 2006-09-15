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

// $Id$
//
// -----------------
// Class AliMUONTest
// -----------------
// Class with functions for testing segmentations
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONTest.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONSegFactory.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModuleTransformer.h"
#include "AliMUONGeometryDetElement.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"

#include "AliMpDEIterator.h"
#include "AliMpDEManager.h"
#include "AliMpVSegmentation.h"
#include "AliMpPad.h"

#include "AliRun.h"
#include "AliLog.h"

#include <TStopwatch.h>
#include <Riostream.h>
#include <TH2F.h>
#include <TPave.h>
#include <TCanvas.h>
#include <TGeoMatrix.h>

ClassImp(AliMUONTest)

//_____________________________________________________________________________
  AliMUONTest::AliMUONTest(const TString& option)
  : TObject(),
    fkTransformer(0),
    fSegmentation(0),
    fCanvas(0)
{
/// Standard Constructor

  if ( option != "default" && 
       option != "FactoryV2" &&
       option != "FactoryV3" && 
       option != "FactoryV4" && 
       option != "new" )  
  {
    BuildWithMUON(option);
  }
  else     
    BuildWithoutMUON(option);

  // Create canvas
  fCanvas = new TCanvas("c1","c1", 0, 0, 800, 800);
  fCanvas->Range(-400,-400, 400, 400);
  fCanvas->cd();
}

//_____________________________________________________________________________
AliMUONTest::AliMUONTest()
  : TObject(),
    fkTransformer(0),
    fSegmentation(0),
    fCanvas(0)
{
/// Default Constructor
}

//_____________________________________________________________________________
AliMUONTest::~AliMUONTest()
{
/// Destructor

  delete fCanvas;
}

//
// private methods
//
//_____________________________________________________________________________
void AliMUONTest::BuildWithMUON(const TString& configMacro)
{
/// Build segmentation via AliMUON initialisation

  gAlice->Init(configMacro.Data());
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  
  fkTransformer = muon->GetGeometryTransformer();
  fSegmentation   = muon->GetSegmentation();
}    

//_____________________________________________________________________________
void AliMUONTest::BuildWithoutMUON(const TString& option)
{
/// Fill geometry from transform*.dat files and build segmentation via 
/// SegFactory

  AliMUONSegFactory segFactory("volpaths.dat", "transform.dat");
  fSegmentation = segFactory.CreateSegmentation(option);
}  


//
// public methods
//

//_____________________________________________________________________________
void AliMUONTest::PrintPadsForAll() const
{
/// Print pads for all chambers and first cathod

  TStopwatch timer;
  timer.Start();  

  // Loop over chambers
  for (Int_t iChamber=0; iChamber<AliMUONConstants::NCh(); iChamber++) {

    // Loop over cathods
    //for (Int_t cath=0; cath<2; cath++) {
    for (Int_t cath=0; cath<1; cath++) {

      cout << setw(6) << "Pads in chamber " << iChamber 
           << " cathod " << cath << endl;
      cout << "===================================" << endl;  

      PrintPadsForSegmentation(iChamber, cath);
    }  
  }     
  timer.Stop();
  timer.Print();
}    

//_____________________________________________________________________________
void AliMUONTest::PrintPadsForSegmentation(Int_t moduleId, Int_t cath) const
{
/// Print pads for the given segmentation
  
  TStopwatch timer;
  timer.Start();  

  // Loop over detection elements
  //
  AliMpDEIterator it;
  for ( it.First(moduleId); ! it.IsDone(); it.Next()) {
       
    Int_t detElemId = it.CurrentDE();       
    cout << "Detection element id: " << detElemId << endl;
    
    PrintPadsForDetElement(detElemId, cath);
  }  

  timer.Stop();
  timer.Print();
} 
   
//_____________________________________________________________________________
void AliMUONTest::PrintPadsForDetElement(Int_t detElemId, Int_t cath) const
{
/// Print global pad positions for the given detection element
  
  
  // Get geometry segmentation
  //
  AliMUONGeometrySegmentation* segmentation 
    = fSegmentation->GetModuleSegmentationByDEId(detElemId, cath);
  if (!segmentation) {
    AliErrorStream()
      << "Segmentation for detElemId = " << detElemId << ", cath = " << cath
      << " not defined." << endl;
    return;
  }  
      	
  Int_t counter = 0;

  // Loop over pads in a detection element
  //
  for(Int_t ix=0; ix<=segmentation->Npx(detElemId); ix++)
    for(Int_t iy=0; iy<=segmentation->Npy(detElemId); iy++) 
       PrintPad(counter, detElemId, ix, iy, segmentation);
} 
   
//_____________________________________________________________________________
void AliMUONTest::PrintPad(Int_t& counter,
                           Int_t detElemId, Int_t ix, Int_t iy,
                           AliMUONGeometrySegmentation* segmentation) const
{
/// Print global pad position for the given pad
  
  Float_t x, y, z;
  Bool_t success
    = segmentation->GetPadC(detElemId, ix, iy, x, y, z);
  
  cout << setw(6) << "counter " << counter++ << "   ";
  cout << "Pad indices:  ( " << detElemId << "; " << ix << ", " << iy << " )  " ;

  if (success) {
    cout << "Pad position: ( " << x << ", " << y << ", " << z << " );  ";
    Int_t sector = segmentation->Sector(detElemId, ix, iy);
    Float_t dpx = segmentation->Dpx(detElemId, sector);
    Float_t dpy = segmentation->Dpy(detElemId, sector);
    cout << " dimensions: ( " << dpx << ", " << dpy << " )" << endl;
  }  
  else  {
    counter--; 
    cout << "... no pad " << endl; 
  }  
} 
   
//_____________________________________________________________________________
void AliMUONTest::DrawPadsForAll() const
{
/// Draw pad for all chambers and first cathod

  TStopwatch timer;
  timer.Start();  

  // Loop over chambers
  for (Int_t iChamber=0; iChamber<AliMUONConstants::NCh(); iChamber++) {

    // Loop over cathods
    //for (Int_t cath=0; cath<2; cath++) {
    for (Int_t cath=0; cath<1; cath++) {

      cout << setw(6) << "Pads in chamber " << iChamber 
           << " cathod " << cath << endl;
      cout << "===================================" << endl;  

      DrawPadsForSegmentation(iChamber, cath);
    }  
  }     
  timer.Stop();
  timer.Print();
}    

//_____________________________________________________________________________
void AliMUONTest::DrawPadsForSegmentation(Int_t moduleId, Int_t cath) const
{
/// Draw pads for the given segmentation
  
  TStopwatch timer;
  timer.Start();  

  // Loop over detection elements
  //
  AliMpDEIterator it;
  for ( it.First(moduleId); ! it.IsDone(); it.Next()) {
       
    Int_t detElemId = it.CurrentDE();       
    cout << "Detection element id: " << detElemId << endl;
    
    DrawPadsForDetElement(detElemId, cath);
  }  

  fCanvas->Update();
  cout << "Print any key + enter to continue ..." << endl;
  char c;
  cin >> c;
  fCanvas->Clear();

  timer.Stop();
  timer.Print();
} 
   
//_____________________________________________________________________________
void AliMUONTest::DrawPadsForDetElement(Int_t detElemId, Int_t cath) const
{
/// Draw pads for the given detection element
  
  // Get geometry segmentation
  //
  AliMUONGeometrySegmentation* segmentation 
    = fSegmentation->GetModuleSegmentationByDEId(detElemId, cath);
  if (!segmentation) {
    AliErrorStream()
      << "Segmentation for detElemId = " << detElemId << ", cath = " << cath
      << " not defined." << endl;
    return;
  }  
      	
  Int_t counter = 0;

  // Loop over pads in a detection element
  //
  for(Int_t ix=0; ix<=segmentation->Npx(detElemId); ix++)
    for(Int_t iy=0; iy<=segmentation->Npy(detElemId); iy++) 
      DrawPad(counter, detElemId, ix, iy, segmentation);
} 
   
//_____________________________________________________________________________
void AliMUONTest::DrawPad(Int_t& counter,
                          Int_t detElemId, Int_t ix, Int_t iy,
                          AliMUONGeometrySegmentation* segmentation) const
{
/// Draw the given pad
  
  Float_t x, y, z;
  Bool_t success
    = segmentation->GetPadC(detElemId, ix, iy, x, y, z);

  if (!success) return;  
  
  // PrintPad(counter,detElemId, ix, iy, segmentation); 

  counter++;
  
  Int_t sector = segmentation->Sector(detElemId, ix, iy);
  Float_t dpx = segmentation->Dpx(detElemId, sector);
  Float_t dpy = segmentation->Dpy(detElemId, sector);

  //printf(" ***** Pad position is ix: %d iy: %d x: %f y: %f sector: %d dpx: %f dpy: %f \n",
  //       ix, iy, x, y, sector, dpx, dpy);

  fCanvas->cd();
  TPave* pave = new TPave(x-dpx/2., y-dpy/2., x+dpx/2., y+dpy/2., 1);
  pave->Draw();
} 

//_____________________________________________________________________________
void  AliMUONTest::DetElemTransforms() const
{
/// Print detection elements transformations 

  // Loop over chambers
  for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {

    const AliMUONGeometryModuleTransformer* kModuleTransformer 
      = fkTransformer->GetModuleTransformer(i);
      
    // Loop over detection elements
    AliMpDEIterator it;
    for ( it.First(i); ! it.IsDone(); it.Next()) {
       
      Int_t detElemId = it.CurrentDE();       
      cout << "Detection element id: " << detElemId << endl;
	
      Double_t x, y, z;
      kModuleTransformer->Local2Global(detElemId, 0., 0., 0., x, y, z);
      cout << "  Global DE position:            " 
	   <<  x << ",  " << y << ",  " << z << endl; 

      Double_t x2, y2, z2;
      kModuleTransformer->Global2Local(detElemId, 0., 0., 0., x2, y2, z2);
      cout << "  ALIC center in the local frame: " 
	   <<  x2 << ",  " << y2 << ",  " << z2 << endl; 
	     
      Double_t x3, y3, z3;
      kModuleTransformer->Global2Local(detElemId, x, y, z, x3, y3, z3);
      cout << "  Back in the local frame: " 
           <<  x3 << ",  " << y3 << ",  " << z3 << endl;        
      cout << endl;	     

      AliMUONGeometryDetElement* detElem =  
        kModuleTransformer->GetDetElement(detElemId);
      detElem->PrintGlobalTransform();	
    }
  }
}    	 

