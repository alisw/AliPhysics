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
//
// Author: Ivana Hrivnacova, IPN Orsay

#include "AliMUONTest.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryBuilder.h"
#include "AliMUONSt1GeometryBuilderV2.h"
#include "AliMUONSt2GeometryBuilderV2.h"
#include "AliMUONSlatGeometryBuilder.h"
#include "AliMUONTriggerGeometryBuilder.h"
#include "AliMUONSegFactory.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometryTransformer.h"
#include "AliMUONSegmentation.h"
#include "AliMUONGeometrySegmentation.h"

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
// Standard Constructor
//

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
}

//_____________________________________________________________________________
AliMUONTest::AliMUONTest()
  : TObject(),
    fkTransformer(0),
    fSegmentation(0),
    fCanvas(0)
{
// Default Constructor
//
}

//_____________________________________________________________________________
AliMUONTest::AliMUONTest(const AliMUONTest& rhs)
 : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//_____________________________________________________________________________
AliMUONTest::~AliMUONTest()
{
// Destructor

  delete fCanvas;
}

//_____________________________________________________________________________
AliMUONTest& AliMUONTest::operator = (const AliMUONTest& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//
// private methods
//
//_____________________________________________________________________________
void AliMUONTest::BuildWithMUON(const TString& configMacro)
{
// Build segmentation via AliMUON initialisation

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
// Fill geometry from transform*.dat files and build segmentation via 
// SegFactory

  AliMUONSegFactory segFactory("volpaths.dat", "transform.dat");
  fSegmentation = segFactory.CreateSegmentation(option);
}  


//
// public methods
//

//_____________________________________________________________________________
AliMUONGeometrySegmentation* 
AliMUONTest::GetSegmentation(Int_t chamberId, Int_t cath)
{
// Create geometry segmentation for the specified chamber and cathod

  return fSegmentation->GetModuleSegmentation(chamberId, cath);
}		

#include "AliMUONGeometryDetElement.h"
//_____________________________________________________________________________
void  AliMUONTest::DetElemTransforms()
{
// 
  // Loop over chambers
  for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {

    const AliMUONGeometryModuleTransformer* kModuleTransformer 
      = fkTransformer->GetModuleTransformer(i);
      
    AliMUONGeometryStore* detElements 
      = kModuleTransformer->GetDetElementStore();
    
    // Loop over detection elements
    for (Int_t j=0; j<detElements->GetNofEntries(); j++) {
       
      //Int_t detElemId = kModuleTransformer->GetDetElemId(j);       
      Int_t detElemId = detElements->GetEntry(j)->GetUniqueID();       
      cout << "Detection element Id: " << detElemId << endl;
	
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
        (AliMUONGeometryDetElement*)detElements->GetEntry(j);
      detElem->PrintGlobalTransform();	
    }
  }
}    	 

//_____________________________________________________________________________
void AliMUONTest::ForWhole(AliMUONTests testCase)
{
// Perform test for all chambers and first cathod

  TStopwatch timer;
  timer.Start();  

  // Loop over chambers
  for (Int_t iChamber=0; iChamber<AliMUONConstants::NCh(); iChamber++) {

    // Loop over cathods
    //for (Int_t cath=0; cath<2; cath++) {
    for (Int_t cath=0; cath<1; cath++) {

      AliMUONGeometrySegmentation* segmentation 
        = GetSegmentation(iChamber, cath);
	
      if (!segmentation) continue;
      	
      cout << setw(6) << "Pads in chamber " << iChamber 
           << " cathod " << cath << endl;
      cout << "===================================" << endl;  

      ForSegmentation(testCase, segmentation);
    }  
  }     
  timer.Stop();
  timer.Print();
}    

//_____________________________________________________________________________
void AliMUONTest::ForSegmentation(AliMUONTests testCase,
                                  AliMUONGeometrySegmentation *segmentation)
{
// Perform test for a given segmentation
  
  TStopwatch timer;
  timer.Start();  

  Before(testCase);

  // Loop over detection elements
  //
  AliMUONGeometryStore* detElements 
    = segmentation->GetTransformer()->GetDetElementStore();
    
  for (Int_t j=0; j<detElements->GetNofEntries(); j++) {
       
    Int_t detElemId = detElements->GetEntry(j)->GetUniqueID();       
    cout << "Detection element id: " << detElemId << endl;
    
    ForDetElement(testCase, detElemId, segmentation);
  }  

  After(testCase);

  timer.Stop();
  timer.Print();
} 
   
//_____________________________________________________________________________
void AliMUONTest::ForDetElement(AliMUONTests testCase,
                                Int_t detElemId,
                                AliMUONGeometrySegmentation *segmentation)
{
// Prints global pad positions for given detection element
// in a given geometry segmentation
  
  
  Int_t counter = 0;

  // Loop over pads in a detection element
  //

  for(Int_t ix=1; ix<=segmentation->Npx(detElemId); ix++)
    for(Int_t iy=1; iy<=segmentation->Npy(detElemId); iy++) 
    {
       switch (testCase) {
     
         case kPrintPads:
           PrintPad(counter, detElemId, ix, iy, segmentation);
	   break;
     
         case kDrawPads:
           DrawPad(counter, detElemId, ix, iy, segmentation);
	   break;
      }
    }    
} 
   
//_____________________________________________________________________________
void AliMUONTest::Before(AliMUONTests testCase)
{
// Do some initialization if necessary

  switch (testCase) {
  
    case kPrintPads:
      break;

    case kDrawPads:
      if (!fCanvas) {
        fCanvas = new TCanvas("c1","c1", 0, 0, 600, 600);
        fCanvas->Range(-300,-300, 300, 300);
	fCanvas->cd();
      }  
      break;
  }        
}

//_____________________________________________________________________________
void AliMUONTest::After(AliMUONTests testCase)
{
// Do some cleanup if necessary

  switch (testCase) {
  
    case kPrintPads:
      break;

    case kDrawPads:
      fCanvas->Update();
      cout << "Print any key + enter to continue ..." << endl;
      char c;
      cin >> c;
      fCanvas->Clear();
      break;
  }        
}

//_____________________________________________________________________________
void AliMUONTest::PrintPad(Int_t& counter,
                           Int_t detElemId, Int_t ix, Int_t iy,
                           AliMUONGeometrySegmentation* segmentation)
{
// Prints global pad positions for the given pad
  
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
void AliMUONTest::DrawPad(Int_t& counter,
                          Int_t detElemId, Int_t ix, Int_t iy,
                          AliMUONGeometrySegmentation* segmentation)
{
// Prints global pad positions for the given pad
  
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

  if (!fCanvas) Before(kDrawPads);

  fCanvas->cd();
  TPave* pave = new TPave(x-dpx/2., y-dpy/2., x+dpx/2., y+dpy/2., 1);
  pave->Draw();
} 
   
//_____________________________________________________________________________
void AliMUONTest::DrawSegmentation(AliMUONGeometrySegmentation *seg)
{
// TBR

  // Drawing slat504
  Int_t ix, iy, deId;
  Float_t x, y, z;
  Float_t dpx, dpy;
//   TH2F * frame = new TH2F(" "," ",10,-10.,245.,10, -5., 45.);
//   TH2F * frame = new TH2F(" "," ",10,-300.,300.,10, -300., 300.);
  TH2F * frame = new TH2F(" "," ",10,-350.,350.,10, -350., 350.);
  frame->Draw();
//   (new TPave(  0.,  0., 40., 40.,2))->Draw();
//   (new TPave( 40.,  0., 80., 40.,2))->Draw();
//   (new TPave( 80.,  0.,120., 40.,2))->Draw();
//   (new TPave(120.,  0.,160., 40.,2))->Draw();
//   (new TPave(160.,  0.,200., 40.,2))->Draw();
//   (new TPave(200.,  0.,240., 40.,2))->Draw();
  
  // Loop over detection elements
  //
  AliMUONGeometryStore* detElements 
    = seg->GetTransformer()->GetDetElementStore();
    
  for (Int_t iDE=0; iDE<detElements->GetNofEntries(); iDE++) {
    
    deId = detElements->GetEntry(iDE)->GetUniqueID();       
	
    cout << "Detection element id: " << deId << endl;
    
    
    //   for ( seg->FirstPad(detElementId,  0., 0., 0., 100., 100.);
    // 	seg->MorePads(detElementId); 
    // 	seg->NextPad(detElementId) ) {
    for(ix=1; ix<=seg->Npx(deId); ix++) {
      for(iy=1; iy<=seg->Npy(deId); iy++) {
	
	seg->GetPadC(deId, ix, iy, x, y, z);
	Int_t sector = seg->Sector(deId, ix, iy);
	dpx = seg->Dpx(deId,sector);
	dpy = seg->Dpy(deId,sector);
	
	//printf(" ***** Pad position is ix: %d iy: %d x: %f y: %f sector: %d dpx: %f dpy: %f \n",ix, iy, x, y, sector, dpx, dpy);
	(new TPave(x-dpx/2., y-dpy/2., x+dpx/2., y+dpy/2., 1))->Draw();
      }
      
    }
  }
}



