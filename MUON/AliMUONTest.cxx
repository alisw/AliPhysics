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
#include <TH2F.h>
#include <TPave.h>

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
#include "AliMUONSt345SlatSegmentation.h"

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
//     cout << "Go to loop over pads" << endl;
//     for (Int_t cath=0; cath<2; cath++) {
      
//       cout << setw(6) << "Pads in chamber " << i << " cathod " << cath << endl;
//       cout << "===================================" << endl;  
//       TStopwatch timer;
//       timer.Start();  

//       Int_t counter = 0;
//       for ( seg[cath]->FirstPad(100,  70., 70., 0., 80., 80.);
//             seg[cath]->MorePads(100); 
//             seg[cath]->NextPad(100) )
//       {
//         cout << setw(6) << "counter " << counter++ << "   ";
  
//         Int_t ix = seg[cath]->Ix();
//         Int_t iy = seg[cath]->Iy();
//         Int_t deId = seg[cath]->DetElemId();
//         cout << "Pad indices:  ( " << deId << "; " << ix << ", " << iy << " )  " ;

//         Float_t x, y, z;
//         seg[cath]->GetPadC(deId, ix, iy, x, y, z);
//         cout << "Pad position: ( " << x << ", " << y << ", " << z << " )" << endl;
//       }
//       timer.Stop();
//       timer.Print();
//     }
    if (i == 0) DrawSegmentation(i,seg[0]);
  } 
}    
//______________________________________________________________________________
void AliMUONTest::DrawSegmentation(Int_t chamber, AliMUONGeometrySegmentation *seg)
{
  // Drawing slat504
  Int_t ix, iy, deId;
  Float_t x, y, z;
  Float_t dpx, dpy;
//   TH2F * frame = new TH2F(" "," ",10,-10.,245.,10, -5., 45.);
  TH2F * frame = new TH2F(" "," ",10,-300.,300.,10, -300., 300.);
//   TH2F * frame = new TH2F(" "," ",10,-300.,300.,10, -25., 25.);
  frame->Draw();
//   (new TPave(  0.,  0., 40., 40.,2))->Draw();
//   (new TPave( 40.,  0., 80., 40.,2))->Draw();
//   (new TPave( 80.,  0.,120., 40.,2))->Draw();
//   (new TPave(120.,  0.,160., 40.,2))->Draw();
//   (new TPave(160.,  0.,200., 40.,2))->Draw();
//   (new TPave(200.,  0.,240., 40.,2))->Draw();
  
  for (Int_t iDE = 0; iDE < 13; iDE++) {
    
    for (Int_t side = 0; side < 2; side++) {
      if (side == 0)
	deId = (chamber+1)*100+iDE;
      else 
	deId = (chamber+1)*100+50+iDE;
	
      
      //   for ( seg->FirstPad(detElementId,  0., 0., 0., 100., 100.);
      // 	seg->MorePads(detElementId); 
      // 	seg->NextPad(detElementId) ) {
      for(ix=1; ix<=seg->Npx(deId); ix++) {
	for(iy=1; iy<=seg->Npy(deId); iy++) {
	  seg->GetPadC(deId, ix, iy, x, y, z);
	  Int_t sector = seg->Sector(deId, ix, iy);
	  dpx = seg->Dpx(deId,sector);
	  dpy = seg->Dpy(deId,sector);
	  //       printf(" ***** Pad position is ix: %d iy: %d x: %f y: %f sector: %d dpx: %f dpy: %f \n",ix, iy, x, y, sector, dpx, dpy);
	  (new TPave(x-dpx/2., y-dpy/2., x+dpx/2., y+dpy/2., 1))->Draw();
	}
      }
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
//_____________________________________________________________________________
void AliMUONTest::St3SlatSegmentation()
{

// Build new segmentations (based on the  detection element local 
// segmentations), iterate over all segmentations and print global 
// pad positions

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  

   //Slats Segmentations
  AliMUONSt345SlatSegmentation *slatseg[20]; // Types of segmentation for St3 in this framework
  Int_t ndiv[4] ={ 4, 4, 2, 1};  // densities zones 
  for(Int_t i=0; i<20;i++) {
    slatseg[i] = new AliMUONSt345SlatSegmentation();
    slatseg[i]->SetPadSize(10.,0.5);
    slatseg[i]->SetPadDivision(ndiv);
    slatseg[i]->SetId(1);
    slatseg[i]->SetDAnod(0.25);
  }
  
  //******************************************************************************************
  // Station 3
  //******************************************************************************************

  // Type 112200 for 500, 501, 508, 509, 510, 517 in Ch5 (similar for Ch6) for the futur official numbering
  // Type 112200 for 503, 504, 505, 555, 554, 553 in Ch5 (similar for Ch6) actual numbering in the code to be changed in jan05
  Int_t n0[4] = { 0, 2, 2, 0 };
  slatseg[0]->SetPcbBoards(n0);
  slatseg[0]->Init(0);
  
  // Type 122200 for 502, 507, 511, 516 (similar in Ch6) for future official numbering of ALICE
  // Type 122200 for 502, 506, 556, 552 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
  Int_t n1[4] = { 0, 1, 3, 0 }; 
  slatseg[1]->SetPcbBoards(n1);
  slatseg[1]->Init(0); 
  
  // Type 222000 for 503, 506, 512, 515 (similar in Ch6) for future official numbering of ALICE
  // Type 222000 for 501, 507, 557, 551 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
  Int_t n2[4] = { 0, 0, 3, 0 };
  slatseg[2]->SetPcbBoards(n2);
  slatseg[2]->Init(0);
  
  // Type 220000 for 504, 505, 513, 514 (similar in Ch6) for future official numbering of ALICE
  // Type 220000 for 500, 508, 558, 550 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
  Int_t n3[4] = { 0, 0, 2, 0 };
  slatseg[3]->SetPcbBoards(n3);
  slatseg[3]->Init(0); 

  //******************************************************************************************
  // Station 4 & 5
  //****************************************************************************************
  
  // Type 122330 for 700, 713 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 122330 for 706, 756 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  // Type 122330 for 900, 913 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 122330 for 906, 956 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
  
  Int_t n4[4] = { 0, 1, 2, 2 };
  slatseg[4]->SetPcbBoards(n4);
  slatseg[4]->Init(0); // 0 detection element id

  // Type 112233 for 701, 712, 714, 725 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 112233 for 705, 707, 755, 757 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  // Type 112233 for 901, 902, 911, 912, 914, 915, 924, 925 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 112233 for 904, 905, 907, 908, 954, 955, 957, 958 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n5[4] = { 0, 2, 2, 2 };
  slatseg[5]->SetPcbBoards(n5);
  slatseg[5]->Init(0); // 0 detection element id

  // Type 112230 for 702, 711, 715, 724 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 112230 for 704, 708, 754, 758 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  Int_t n6[4] = { 0, 2, 2, 1 };
  slatseg[6]->SetPcbBoards(n6);
  slatseg[6]->Init(0); // 0 detection element id

  // Type 222330 for 703, 710, 716, 723 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 222330 for 703, 709, 753, 759 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  Int_t n7[4] = { 0, 0, 3, 2 };
  slatseg[7]->SetPcbBoards(n7);
  slatseg[7]->Init(0); // 0 detection element id

  // Type 223300 for 704, 709, 717, 722 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 223300 for 702, 710, 752, 760 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  Int_t n8[4] = { 0, 0, 2, 2 };
  slatseg[8]->SetPcbBoards(n8);
  slatseg[8]->Init(0); // 0 detection element id

  // Type 333000 for 705, 708, 718, 721 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 333000 for 701, 711, 751, 761 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  // Type 333000 for 906, 907, 919, 920 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 333000 for 900, 912, 950, 962 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
   Int_t n9[4] = { 0, 0, 0, 3 };
  slatseg[9]->SetPcbBoards(n9);
  slatseg[9]->Init(0); // 0 detection element id

  // Type 330000 for 706, 707, 719, 720 in Ch7 (similar for Ch8) for the futur official numbering
  // Type 330000 for 700, 712, 750, 762 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
  Int_t n10[4] = { 0, 0, 0, 2 };
  slatseg[10]->SetPcbBoards(n10);
  slatseg[10]->Init(0); // 0 detection element id

  // Type 222333 for 903, 910, 916, 923 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 222333 for 903, 909, 953, 959 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n11[4] = { 0, 0, 3, 3 };
  slatseg[11]->SetPcbBoards(n11);
  slatseg[11]->Init(0); // 0 detection element id
         
  // Type 223330 for 904, 909, 917, 922 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 223330 for 902, 910, 952, 960 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n12[4] = { 0, 0, 2, 3 };
  slatseg[12]->SetPcbBoards(n12);
  slatseg[12]->Init(0); // 0 detection element id

  // Type 333300 for 905, 908, 918, 921 in Ch9 (similar for Ch10) for the futur official numbering
  // Type 333300 for 901, 911, 951, 961 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n13[4] = { 0, 0, 0, 4 };
  slatseg[13]->SetPcbBoards(n13);
  slatseg[13]->Init(0); // 0 detection element id
  
 
  //Loop for St3 (only segmentation for bending plane)
  for (Int_t i=0; i<2; i++) {
    // Create chamber segmentations  
    Int_t chamber = i+4;
    AliMUONGeometrySegmentation* st3seg = new AliMUONGeometrySegmentation(muon->Chamber(chamber).GetGeometry());
    
    Int_t id0=(chamber+1)*100;
    // type 220000
    st3seg->Add(id0, slatseg[3]);
    st3seg->Add(id0+ 8, slatseg[3]);  
    st3seg->Add(id0+50, slatseg[3]);  
    st3seg->Add(id0+58, slatseg[3]);
    // type 222000
    st3seg->Add(id0+ 1, slatseg[2]);
    st3seg->Add(id0+ 7, slatseg[2]);  
    st3seg->Add(id0+51, slatseg[2]);  
    st3seg->Add(id0+57, slatseg[2]);
    // type 122200
    st3seg->Add(id0+ 2, slatseg[1]);
    st3seg->Add(id0+ 6, slatseg[1]);  
    st3seg->Add(id0+52, slatseg[1]);  
    st3seg->Add(id0+56, slatseg[1]);
    // type 112200
    st3seg->Add(id0+ 3, slatseg[0]);
    st3seg->Add(id0+ 4, slatseg[0]);  
    st3seg->Add(id0+ 5, slatseg[0]);  
    st3seg->Add(id0+53, slatseg[0]);
    st3seg->Add(id0+54, slatseg[0]);     
    st3seg->Add(id0+55, slatseg[0]);
    
//     if (i == 0) DrawSegmentation(chamber,st3seg);
  }
  
  //Loop for St4 (only segmentation for bending plane)
  for (Int_t i=0; i<2; i++) {
    // Create chamber segmentations  
    Int_t chamber = i+6;
    AliMUONGeometrySegmentation* st4seg = new AliMUONGeometrySegmentation(muon->Chamber(chamber).GetGeometry());
    
    Int_t id0=(chamber+1)*100;
    // type 122330
    st4seg->Add(id0+ 6, slatseg[4]);
    st4seg->Add(id0+56, slatseg[4]);
    // type 112233
    st4seg->Add(id0+ 5, slatseg[5]);
    st4seg->Add(id0+ 7, slatseg[5]);  
    st4seg->Add(id0+55, slatseg[5]);  
    st4seg->Add(id0+57, slatseg[5]);
    // type 112230
    st4seg->Add(id0+ 4, slatseg[6]);
    st4seg->Add(id0+ 8, slatseg[6]);  
    st4seg->Add(id0+54, slatseg[6]);  
    st4seg->Add(id0+58, slatseg[6]);
    // type 222330 
    st4seg->Add(id0+ 3, slatseg[7]);
    st4seg->Add(id0+ 9, slatseg[7]);  
    st4seg->Add(id0+53, slatseg[7]);
    st4seg->Add(id0+59, slatseg[7]);
    // type 223300 
    st4seg->Add(id0+ 2, slatseg[8]);
    st4seg->Add(id0+10, slatseg[8]);  
    st4seg->Add(id0+52, slatseg[8]);
    st4seg->Add(id0+60, slatseg[8]);
    // type 333000 
    st4seg->Add(id0+ 1, slatseg[9]);
    st4seg->Add(id0+11, slatseg[9]);  
    st4seg->Add(id0+51, slatseg[9]);
    st4seg->Add(id0+61, slatseg[9]);
    // type 330000 
    st4seg->Add(id0   , slatseg[10]);
    st4seg->Add(id0+12, slatseg[10]);  
    st4seg->Add(id0+50, slatseg[10]);
    st4seg->Add(id0+62, slatseg[10]);
    
//     if (i == 0) DrawSegmentation(chamber,st4seg);
  }

  //Loop for St5 (only segmentation for bending plane)
  for (Int_t i=0; i<2; i++) {
    // Create chamber segmentations  
    Int_t chamber = i+8;
    AliMUONGeometrySegmentation* st5seg = new AliMUONGeometrySegmentation(muon->Chamber(chamber).GetGeometry());
    
    Int_t id0=(chamber+1)*100;
    // type 122330
    st5seg->Add(id0+ 6, slatseg[4]);
    st5seg->Add(id0+56, slatseg[4]);
    // type 112233
    st5seg->Add(id0+ 4, slatseg[5]);
    st5seg->Add(id0+ 5, slatseg[5]);
    st5seg->Add(id0+ 7, slatseg[5]);  
    st5seg->Add(id0+ 8, slatseg[5]);  
    st5seg->Add(id0+54, slatseg[5]);  
    st5seg->Add(id0+55, slatseg[5]);  
    st5seg->Add(id0+57, slatseg[5]);
    st5seg->Add(id0+58, slatseg[5]);
    // type 222333 
    st5seg->Add(id0+ 3, slatseg[11]);
    st5seg->Add(id0+ 9, slatseg[11]);  
    st5seg->Add(id0+53, slatseg[11]);
    st5seg->Add(id0+59, slatseg[11]);
    // type 223330 
    st5seg->Add(id0+ 2, slatseg[12]);
    st5seg->Add(id0+10, slatseg[12]);  
    st5seg->Add(id0+52, slatseg[12]);
    st5seg->Add(id0+60, slatseg[12]);
    // type 333300 
    st5seg->Add(id0+ 1, slatseg[13]);
    st5seg->Add(id0+11, slatseg[13]);  
    st5seg->Add(id0+51, slatseg[13]);
    st5seg->Add(id0+61, slatseg[13]);
    // type 333000 
    st5seg->Add(id0   , slatseg[9]);
    st5seg->Add(id0+12, slatseg[9]);  
    st5seg->Add(id0+50, slatseg[9]);
    st5seg->Add(id0+62, slatseg[9]);
    
    if (i == 0) DrawSegmentation(chamber,st5seg);
  }

} 
