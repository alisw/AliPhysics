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
// Class AliMUONTest
// -----------------
// Class with functions for testing
//
// Author: Ivana Hrivnacova, IPN Orsay

#include <TStopwatch.h>
#include <Riostream.h>

#include "AliRun.h"
#include "AliSegmentation.h"
#include "AliLog.h"

#include "AliMUONTest.h"
#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONSt12QuadrantSegmentation.h"

ClassImp(AliMUONTest)

//__________________________________________________________________
AliMUONTest::AliMUONTest(const TString& configMacro)
  : TObject()
{
// Standard Constructor
//
  // Initialize AliRoot
  gAlice->Init(configMacro.Data());
}

//__________________________________________________________________
AliMUONTest::AliMUONTest()
  : TObject()
{
// Default Constructor
//
}

//____________________________________________________________________
AliMUONTest::AliMUONTest(const AliMUONTest& rhs)
 : TObject(rhs)
{
// Protected copy constructor

  AliFatal("Not implemented.");
}

//____________________________________________________________________
AliMUONTest::~AliMUONTest()
{
// Destructor
}

//________________________________________________________________________
AliMUONTest& AliMUONTest::operator = (const AliMUONTest& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}

//
// public methods
//

//______________________________________________________________________________
void  AliMUONTest::DetElemTransforms()
{
// 
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  
  
  // Loop over chambers
  for (Int_t i=0; i<AliMUONConstants::NCh(); i++) {

    AliMUONGeometryModule* geometry = muon->Chamber(i).GetGeometry();
    AliMUONGeometryStore* detElements = geometry->GetDetElementStore();
    
    // Loop over detection elements
    for (Int_t j=0; j<detElements->GetNofEntries(); j++) {
       
      //Int_t detElemId = geometry->GetDetElemId(j);       
      Int_t detElemId = detElements->GetEntry(j)->GetUniqueID();       
      cout << "Detection element Id: " << detElemId << endl;
	
      Double_t x, y, z;
      geometry->Local2Global(detElemId, 0., 0., 0., x, y, z);
      cout << "  Global DE position:            " 
	   <<  x << ",  " << y << ",  " << z << endl; 

      Double_t x2, y2, z2;
      geometry->Global2Local(detElemId, 0., 0., 0., x2, y2, z2);
      cout << "  ALIC center in the local frame: " 
	   <<  x2 << ",  " << y2 << ",  " << z2 << endl; 
	     
      Double_t x3, y3, z3;
      geometry->Global2Local(detElemId, x, y, z, x3, y3, z3);
      cout << "  Back in the local frame: " 
           <<  x3 << ",  " << y3 << ",  " << z3 << endl;        
      cout << endl;	     
    }
  }
}    	 

//______________________________________________________________________________
void AliMUONTest::PrintPadPositions1()
{
// Build new segmentations (based on the  detection element local 
// segmentations), iterate over all segmentations and print global 
// pad positions

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  

  // Loop over chambers
  for (Int_t i=0; i<1; i++) {

    // Create chamber segmentations
    AliMUONGeometrySegmentation* seg[2];
    seg[0] = new AliMUONGeometrySegmentation(muon->Chamber(i).GetGeometry());
    seg[1] = new AliMUONGeometrySegmentation(muon->Chamber(i).GetGeometry());
    
    // Quadrant segmentations:
    AliMUONSt12QuadrantSegmentation* bendSt1
      = new AliMUONSt12QuadrantSegmentation(kStation1, kBendingPlane);
    AliMUONSt12QuadrantSegmentation* nonbendSt1
      = new AliMUONSt12QuadrantSegmentation(kStation1, kNonBendingPlane);

    // The same configuration for both chambers of Station 1
    Int_t id0 = (i+1)*100;

    // Configure  St1 chamber segmentations
    seg[0]->Add(id0,      bendSt1);
    seg[0]->Add(id0 +  1, nonbendSt1);
    seg[0]->Add(id0 + 50, nonbendSt1);
    seg[0]->Add(id0 + 51, bendSt1);

    seg[1]->Add(id0,      nonbendSt1);
    seg[1]->Add(id0 +  1, bendSt1);
    seg[1]->Add(id0 + 50, bendSt1);
    seg[1]->Add(id0 + 51, nonbendSt1);

    // Iterate over the whole plane and return pad indices and 
    // global/local positions
    cout << "Go to loop over pads" << endl;
    for (Int_t cath=0; cath<2; cath++) {
      
      cout << setw(6) << "Pads in chamber " << i << " cathod " << cath << endl;
      cout << "===================================" << endl;  
      TStopwatch timer;
      timer.Start();  

      Int_t counter = 0;
      for ( seg[cath]->FirstPad(100,  70., 70., 0., 80., 80.);
            seg[cath]->MorePads(100); 
            seg[cath]->NextPad(100) )
      {
        cout << setw(6) << "counter " << counter++ << "   ";
  
        Int_t ix = seg[cath]->Ix();
        Int_t iy = seg[cath]->Iy();
        Int_t deId = seg[cath]->DetElemId();
        cout << "Pad indices:  ( " << deId << "; " << ix << ", " << iy << " )  " ;

        Float_t x, y, z;
        seg[cath]->GetPadC(deId, ix, iy, x, y, z);
        cout << "Pad position: ( " << x << ", " << y << ", " << z << " )" << endl;
      }
      timer.Stop();
      timer.Print();
    }  
  }  
}    

//________________________________________________________________________
void AliMUONTest::PrintPadPositions2()
{
// Iterate over all chamber segmentations and prints
// global pad positions

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  

  // Loop over chambers
  for (Int_t i=0; i<1; i++) {

    // Create chamber segmentations
    AliSegmentation* seg[2];
    seg[0] = muon->Chamber(i).SegmentationModel(0);
    seg[1] = muon->Chamber(i).SegmentationModel(1);

    // Iterate over the whole plane and return pad indices and 
    // global/local positions
    cout << "Go to loop over pads" << endl;
    for (Int_t cath=0; cath<2; cath++) {
      
      cout << setw(6) << "Pads in chamber " << i << " cathod " << cath << endl;
      cout << "===================================" << endl;  
      TStopwatch timer;
      timer.Start();  

      Int_t counter = 0;
      for ( seg[cath]->FirstPad(70., 70., 0., 80., 80.);
            seg[cath]->MorePads(); 
            seg[cath]->NextPad() )
      {
        cout << setw(6) << "counter " << counter++ << "   ";
  
        Int_t ix = seg[cath]->Ix();
        Int_t iy = seg[cath]->Iy();
        cout << "Pad indices:  ( " << ix << ", " << iy << " )  " ;

        Float_t x, y, z;
        seg[cath]->GetPadC(ix, iy, x, y, z);
        cout << "Pad position: ( " << x << ", " << y << ", " << z << " )" << endl;
      }
      timer.Stop();
      timer.Print();
    }  
  }  
}
 
