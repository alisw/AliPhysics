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
#include <TCanvas.h>

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
#include "AliMUONTriggerSegmentation.h"
#include "AliMUONTriggerConstants.h"

ClassImp(AliMUONTest)

//__________________________________________________________________
  AliMUONTest::AliMUONTest(const TString& configMacro)
  : TObject(),
    fCanvas(0)
{
// Standard Constructor
//
  // Initialize AliRoot
  gAlice->Init(configMacro.Data());
}

//__________________________________________________________________
AliMUONTest::AliMUONTest()
  : TObject(),
    fCanvas(0)
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

  delete fCanvas;
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
// private methods
//

//________________________________________________________________________
AliMUONGeometrySegmentation*     
AliMUONTest::CreateSt1Segmentation(Int_t chamberId, Int_t cathod)
{
// Create St1 geometry segmentation for given chamber and cathod

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return 0;   
  }  

  AliMUONGeometrySegmentation* segmentation
    = new AliMUONGeometrySegmentation(muon->Chamber(chamberId).GetGeometry());
    
  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kNonBendingPlane);

  // The same configuration for both chambers of Station 1
  Int_t id0 = (chamberId+1)*100;

  // Configure  St1 chamber segmentations
  if (cathod == 0) {
    segmentation->Add(id0,      bendSt1);
    segmentation->Add(id0 +  1, nonbendSt1);
    segmentation->Add(id0 + 50, bendSt1);
    segmentation->Add(id0 + 51, nonbendSt1);
  }
  else if (cathod == 1) {
    segmentation->Add(id0,      nonbendSt1);
    segmentation->Add(id0 +  1, bendSt1);
    segmentation->Add(id0 + 50, nonbendSt1);
    segmentation->Add(id0 + 51, bendSt1);
  }
  else {
    AliError("Wrong cathod number");
    return 0;
  }
  
  return segmentation;
}      

//_____________________________________________________________________________
AliMUONGeometrySegmentation*     
AliMUONTest::CreateSlatSegmentation(Int_t chamberId, Int_t cathod)
{
// Create St1 geometry segmentation for given chamber and cathod

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return 0;
  }  

  AliMUONGeometrySegmentation *chamberSeg = new AliMUONGeometrySegmentation(muon->Chamber(chamberId).GetGeometry());

   //Slats Segmentations
  AliMUONSt345SlatSegmentation *slatsegB[14]; // Types of segmentation for St3 in this framework
  AliMUONSt345SlatSegmentation *slatsegNB[14]; // Types of segmentation for St3 in this framework
                                               // Bending
  Int_t ndiv[4] ={ 4, 4, 2, 1};  // densities zones 
  for(Int_t i=0; i<14; i++) {
    slatsegB[i] = new AliMUONSt345SlatSegmentation(1);
    slatsegB[i]->SetPadSize(10.,0.5);
    slatsegB[i]->SetPadDivision(ndiv);
    slatsegB[i]->SetId(1);
    slatsegB[i]->SetDAnod(0.25);
    slatsegNB[i] = new AliMUONSt345SlatSegmentation(0);
    slatsegNB[i]->SetPadSize(1.,10.); // must be 0.713 !!!
    slatsegNB[i]->SetPadDivision(ndiv);
    slatsegNB[i]->SetId(1);
    slatsegNB[i]->SetDAnod(0.25);
  }
  
  
  //******************************************************************************************
  // Station 3
  //******************************************************************************************

  if (chamberId == 4 || chamberId == 5) {
    // Type 112200 for 500, 501, 508, 509, 510, 517 in Ch5 (similar for Ch6) for the futur official numbering
    // Type 112200 for 503, 504, 505, 555, 554, 553 in Ch5 (similar for Ch6) actual numbering in the code to be changed in jan05
    Int_t n0[4] = { 0, 2, 2, 0 };
    slatsegB[0]->SetPcbBoards(n0);
    slatsegB[0]->Init(0);
    slatsegNB[0]->SetPcbBoards(n0);
    slatsegNB[0]->Init(0);
    
    // Type 122200 for 502, 507, 511, 516 (similar in Ch6) for future official numbering of ALICE
    // Type 122200 for 502, 506, 556, 552 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
    Int_t n1[4] = { 0, 1, 3, 0 }; 
    slatsegB[1]->SetPcbBoards(n1);
    slatsegB[1]->Init(0); 
    slatsegNB[1]->SetPcbBoards(n1);
    slatsegNB[1]->Init(0); 
    
    // Type 222000 for 503, 506, 512, 515 (similar in Ch6) for future official numbering of ALICE
    // Type 222000 for 501, 507, 557, 551 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
    Int_t n2[4] = { 0, 0, 3, 0 };
    slatsegB[2]->SetPcbBoards(n2);
    slatsegB[2]->Init(0);
    slatsegNB[2]->SetPcbBoards(n2);
    slatsegNB[2]->Init(0);
    
    // Type 220000 for 504, 505, 513, 514 (similar in Ch6) for future official numbering of ALICE
    // Type 220000 for 500, 508, 558, 550 (similiarin Ch6) for actual numbering in muon code to be changed in jan05
    Int_t n3[4] = { 0, 0, 2, 0 };
    slatsegB[3]->SetPcbBoards(n3);
    slatsegB[3]->Init(0); 
    slatsegNB[3]->SetPcbBoards(n3);
    slatsegNB[3]->Init(0); 
  }
  
  
  //***************************************************************************************
  // Station 4 & 5
  //****************************************************************************************
  
  if (chamberId >= 6 && chamberId <= 9) {
    // Type 122330 for 700, 713 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 122330 for 706, 756 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    // Type 122330 for 900, 913 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 122330 for 906, 956 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    
    Int_t n4[4] = { 0, 1, 2, 2 };
    slatsegB[4]->SetPcbBoards(n4);
    slatsegB[4]->Init(0); // 0 detection element id
    slatsegNB[4]->SetPcbBoards(n4);
    slatsegNB[4]->Init(0); // 0 detection element id
    
    // Type 112233 for 701, 712, 714, 725 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 112233 for 705, 707, 755, 757 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    // Type 112233 for 901, 902, 911, 912, 914, 915, 924, 925 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 112233 for 904, 905, 907, 908, 954, 955, 957, 958 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    Int_t n5[4] = { 0, 2, 2, 2 };
    slatsegB[5]->SetPcbBoards(n5);
    slatsegB[5]->Init(0); // 0 detection element id
    slatsegNB[5]->SetPcbBoards(n5);
    slatsegNB[5]->Init(0); // 0 detection element id
    
    // Type 112230 for 702, 711, 715, 724 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 112230 for 704, 708, 754, 758 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    Int_t n6[4] = { 0, 2, 2, 1 };
    slatsegB[6]->SetPcbBoards(n6);
    slatsegB[6]->Init(0); // 0 detection element id
    slatsegNB[6]->SetPcbBoards(n6);
    slatsegNB[6]->Init(0); // 0 detection element id
    
    // Type 222330 for 703, 710, 716, 723 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 222330 for 703, 709, 753, 759 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    Int_t n7[4] = { 0, 0, 3, 2 };
    slatsegB[7]->SetPcbBoards(n7);
    slatsegB[7]->Init(0); // 0 detection element id
    slatsegNB[7]->SetPcbBoards(n7);
    slatsegNB[7]->Init(0); // 0 detection element id
    
    // Type 223300 for 704, 709, 717, 722 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 223300 for 702, 710, 752, 760 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    Int_t n8[4] = { 0, 0, 2, 2 };
    slatsegB[8]->SetPcbBoards(n8);
    slatsegB[8]->Init(0); // 0 detection element id
    slatsegNB[8]->SetPcbBoards(n8);
    slatsegNB[8]->Init(0); // 0 detection element id
    
    // Type 333000 for 705, 708, 718, 721 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 333000 for 701, 711, 751, 761 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    // Type 333000 for 906, 907, 919, 920 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 333000 for 900, 912, 950, 962 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    Int_t n9[4] = { 0, 0, 0, 3 };
    slatsegB[9]->SetPcbBoards(n9);
    slatsegB[9]->Init(0); // 0 detection element id
    slatsegNB[9]->SetPcbBoards(n9);
    slatsegNB[9]->Init(0); // 0 detection element id
    
    // Type 330000 for 706, 707, 719, 720 in Ch7 (similar for Ch8) for the futur official numbering
    // Type 330000 for 700, 712, 750, 762 in Ch7 (similar for Ch8) actual numbering in the code to be changed in jan05
    Int_t n10[4] = { 0, 0, 0, 2 };
    slatsegB[10]->SetPcbBoards(n10);
    slatsegB[10]->Init(0); // 0 detection element id
    slatsegNB[10]->SetPcbBoards(n10);
    slatsegNB[10]->Init(0); // 0 detection element id
    
    // Type 222333 for 903, 910, 916, 923 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 222333 for 903, 909, 953, 959 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    Int_t n11[4] = { 0, 0, 3, 3 };
    slatsegB[11]->SetPcbBoards(n11);
    slatsegB[11]->Init(0); // 0 detection element id
    slatsegNB[11]->SetPcbBoards(n11);
    slatsegNB[11]->Init(0); // 0 detection element id
  
    // Type 223330 for 904, 909, 917, 922 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 223330 for 902, 910, 952, 960 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    Int_t n12[4] = { 0, 0, 2, 3 };
    slatsegB[12]->SetPcbBoards(n12);
    slatsegB[12]->Init(0); // 0 detection element id
    slatsegNB[12]->SetPcbBoards(n12);
    slatsegNB[12]->Init(0); // 0 detection element id
    
    // Type 333300 for 905, 908, 918, 921 in Ch9 (similar for Ch10) for the futur official numbering
    // Type 333300 for 901, 911, 951, 961 in Ch9 (similar for Ch10) actual numbering in the code to be changed in jan05
    Int_t n13[4] = { 0, 0, 0, 4 };
    slatsegB[13]->SetPcbBoards(n13);
    slatsegB[13]->Init(0); // 0 detection element id
    slatsegNB[13]->SetPcbBoards(n13);
    slatsegNB[13]->Init(0); // 0 detection element id
    
  }

  Int_t id0 = 0;

  // For St3 
  if (chamberId == 4 || chamberId == 5) {
    // Create chamber segmentations  

    id0=(chamberId+1)*100;
    // type 220000
    if (cathod == 0) {
      chamberSeg->Add(id0, slatsegB[3]);
      chamberSeg->Add(id0+ 8, slatsegB[3]);  
      chamberSeg->Add(id0+50, slatsegB[3]);  
      chamberSeg->Add(id0+58, slatsegB[3]);
    } else {
      chamberSeg->Add(id0, slatsegNB[3]);
      chamberSeg->Add(id0+ 8, slatsegNB[3]);  
      chamberSeg->Add(id0+50, slatsegNB[3]);  
      chamberSeg->Add(id0+58, slatsegNB[3]);
    }
    // type 222000
    if (cathod == 0) {
      chamberSeg->Add(id0+ 1, slatsegB[2]);
      chamberSeg->Add(id0+ 7, slatsegB[2]);  
      chamberSeg->Add(id0+51, slatsegB[2]);  
      chamberSeg->Add(id0+57, slatsegB[2]);
    } else {
      chamberSeg->Add(id0+ 1, slatsegNB[2]);
      chamberSeg->Add(id0+ 7, slatsegNB[2]);  
      chamberSeg->Add(id0+51, slatsegNB[2]);  
      chamberSeg->Add(id0+57, slatsegNB[2]);
    }
    // type 122200
    if (cathod == 0) {
      chamberSeg->Add(id0+ 2, slatsegB[1]);
      chamberSeg->Add(id0+ 6, slatsegB[1]);  
      chamberSeg->Add(id0+52, slatsegB[1]);  
      chamberSeg->Add(id0+56, slatsegB[1]);
    } else {
      chamberSeg->Add(id0+ 2, slatsegNB[1]);
      chamberSeg->Add(id0+ 6, slatsegNB[1]);  
      chamberSeg->Add(id0+52, slatsegNB[1]);  
      chamberSeg->Add(id0+56, slatsegNB[1]);
    }
    // type 112200
    if (cathod == 0) {
      chamberSeg->Add(id0+ 3, slatsegB[0]);
      chamberSeg->Add(id0+ 4, slatsegB[0]);  
      chamberSeg->Add(id0+ 5, slatsegB[0]);  
      chamberSeg->Add(id0+53, slatsegB[0]);
      chamberSeg->Add(id0+54, slatsegB[0]);     
      chamberSeg->Add(id0+55, slatsegB[0]);
    } else {
      chamberSeg->Add(id0+ 3, slatsegNB[0]);
      chamberSeg->Add(id0+ 4, slatsegNB[0]);  
      chamberSeg->Add(id0+ 5, slatsegNB[0]);  
      chamberSeg->Add(id0+53, slatsegNB[0]);
      chamberSeg->Add(id0+54, slatsegNB[0]);     
      chamberSeg->Add(id0+55, slatsegNB[0]);
    }

  }
  
  // For St4 
  if (chamberId == 6 || chamberId == 7) {
    // Create chamber segmentations  
    id0=(chamberId+1)*100;
    // type 122330
    if (cathod == 0) {
      chamberSeg->Add(id0+ 6, slatsegB[4]);
      chamberSeg->Add(id0+56, slatsegB[4]);
    } else {
      chamberSeg->Add(id0+ 6, slatsegNB[4]);
      chamberSeg->Add(id0+56, slatsegNB[4]);
    }
    // type 112233
    if (cathod == 0) {
      chamberSeg->Add(id0+ 5, slatsegB[5]);
      chamberSeg->Add(id0+ 7, slatsegB[5]);  
      chamberSeg->Add(id0+55, slatsegB[5]);  
      chamberSeg->Add(id0+57, slatsegB[5]);
    } else {
      chamberSeg->Add(id0+ 5, slatsegNB[5]);
      chamberSeg->Add(id0+ 7, slatsegNB[5]);  
      chamberSeg->Add(id0+55, slatsegNB[5]);  
      chamberSeg->Add(id0+57, slatsegNB[5]);
    }
    // type 112230
    if (cathod == 0) {
      chamberSeg->Add(id0+ 4, slatsegB[6]);
      chamberSeg->Add(id0+ 8, slatsegB[6]);  
      chamberSeg->Add(id0+54, slatsegB[6]);  
      chamberSeg->Add(id0+58, slatsegB[6]);
    } else {
      chamberSeg->Add(id0+ 4, slatsegNB[6]);
      chamberSeg->Add(id0+ 8, slatsegNB[6]);  
      chamberSeg->Add(id0+54, slatsegNB[6]);  
      chamberSeg->Add(id0+58, slatsegNB[6]);
    }
    // type 222330 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 3, slatsegB[7]);
      chamberSeg->Add(id0+ 9, slatsegB[7]);  
      chamberSeg->Add(id0+53, slatsegB[7]);
      chamberSeg->Add(id0+59, slatsegB[7]);
    } else {
      chamberSeg->Add(id0+ 3, slatsegNB[7]);
      chamberSeg->Add(id0+ 9, slatsegNB[7]);  
      chamberSeg->Add(id0+53, slatsegNB[7]);
      chamberSeg->Add(id0+59, slatsegNB[7]);
    }
    // type 223300 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 2, slatsegB[8]);
      chamberSeg->Add(id0+10, slatsegB[8]);  
      chamberSeg->Add(id0+52, slatsegB[8]);
      chamberSeg->Add(id0+60, slatsegB[8]);
    } else {
      chamberSeg->Add(id0+ 2, slatsegNB[8]);
      chamberSeg->Add(id0+10, slatsegNB[8]);  
      chamberSeg->Add(id0+52, slatsegNB[8]);
      chamberSeg->Add(id0+60, slatsegNB[8]);
    }
    // type 333000 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 1, slatsegB[9]);
      chamberSeg->Add(id0+11, slatsegB[9]);  
      chamberSeg->Add(id0+51, slatsegB[9]);
      chamberSeg->Add(id0+61, slatsegB[9]);
    } else {
      chamberSeg->Add(id0+ 1, slatsegNB[9]);
      chamberSeg->Add(id0+11, slatsegNB[9]);  
      chamberSeg->Add(id0+51, slatsegNB[9]);
      chamberSeg->Add(id0+61, slatsegNB[9]);
    }
    // type 330000 
    if (cathod == 0) {
      chamberSeg->Add(id0   , slatsegB[10]);
      chamberSeg->Add(id0+12, slatsegB[10]);  
      chamberSeg->Add(id0+50, slatsegB[10]);
      chamberSeg->Add(id0+62, slatsegB[10]);
    } else {
      chamberSeg->Add(id0   , slatsegNB[10]);
      chamberSeg->Add(id0+12, slatsegNB[10]);  
      chamberSeg->Add(id0+50, slatsegNB[10]);
      chamberSeg->Add(id0+62, slatsegNB[10]);
    }
  }

  // For St5 
  if (chamberId == 8 || chamberId == 9) {
    // Create chamber segmentations      
    id0=(chamberId+1)*100;
    // type 122330
    if (cathod == 0) {
      chamberSeg->Add(id0+ 6, slatsegB[4]);
      chamberSeg->Add(id0+56, slatsegB[4]);
    } else {
      chamberSeg->Add(id0+ 6, slatsegNB[4]);
      chamberSeg->Add(id0+56, slatsegNB[4]);
    }
    // type 112233
    if (cathod == 0) {
      chamberSeg->Add(id0+ 4, slatsegB[5]);
      chamberSeg->Add(id0+ 5, slatsegB[5]);
      chamberSeg->Add(id0+ 7, slatsegB[5]);  
      chamberSeg->Add(id0+ 8, slatsegB[5]);  
      chamberSeg->Add(id0+54, slatsegB[5]);  
      chamberSeg->Add(id0+55, slatsegB[5]);  
      chamberSeg->Add(id0+57, slatsegB[5]);
      chamberSeg->Add(id0+58, slatsegB[5]);
    } else {
      chamberSeg->Add(id0+ 4, slatsegNB[5]);
      chamberSeg->Add(id0+ 5, slatsegNB[5]);
      chamberSeg->Add(id0+ 7, slatsegNB[5]);  
      chamberSeg->Add(id0+ 8, slatsegNB[5]);  
      chamberSeg->Add(id0+54, slatsegNB[5]);  
      chamberSeg->Add(id0+55, slatsegNB[5]);  
      chamberSeg->Add(id0+57, slatsegNB[5]);
      chamberSeg->Add(id0+58, slatsegNB[5]);
    }
    // type 222333 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 3, slatsegB[11]);
      chamberSeg->Add(id0+ 9, slatsegB[11]);  
      chamberSeg->Add(id0+53, slatsegB[11]);
      chamberSeg->Add(id0+59, slatsegB[11]);
    } else {
      chamberSeg->Add(id0+ 3, slatsegNB[11]);
      chamberSeg->Add(id0+ 9, slatsegNB[11]);  
      chamberSeg->Add(id0+53, slatsegNB[11]);
      chamberSeg->Add(id0+59, slatsegNB[11]);
    }
    // type 223330 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 2, slatsegB[12]);
      chamberSeg->Add(id0+10, slatsegB[12]);  
      chamberSeg->Add(id0+52, slatsegB[12]);
      chamberSeg->Add(id0+60, slatsegB[12]);
    } else {
      chamberSeg->Add(id0+ 2, slatsegNB[12]);
      chamberSeg->Add(id0+10, slatsegNB[12]);  
      chamberSeg->Add(id0+52, slatsegNB[12]);
      chamberSeg->Add(id0+60, slatsegNB[12]);
    }
    // type 333300 
    if (cathod == 0) {
      chamberSeg->Add(id0+ 1, slatsegB[13]);
      chamberSeg->Add(id0+11, slatsegB[13]);  
      chamberSeg->Add(id0+51, slatsegB[13]);
      chamberSeg->Add(id0+61, slatsegB[13]);
    } else {
      chamberSeg->Add(id0+ 1, slatsegNB[13]);
      chamberSeg->Add(id0+11, slatsegNB[13]);  
      chamberSeg->Add(id0+51, slatsegNB[13]);
      chamberSeg->Add(id0+61, slatsegNB[13]);
    }
    // type 333000 
    if (cathod == 0) {
      chamberSeg->Add(id0   , slatsegB[9]);
      chamberSeg->Add(id0+12, slatsegB[9]);  
      chamberSeg->Add(id0+50, slatsegB[9]);
      chamberSeg->Add(id0+62, slatsegB[9]);
    } else {
      chamberSeg->Add(id0   , slatsegNB[9]);
      chamberSeg->Add(id0+12, slatsegNB[9]);  
      chamberSeg->Add(id0+50, slatsegNB[9]);
      chamberSeg->Add(id0+62, slatsegNB[9]);
    }
  }  
  
 
  if (!id0) {
    AliWarning(Form("Segmentation for chamber %d , cathod %d is not yet defined",
		    chamberId, cathod));
    return 0;
    
  }
  
  DrawSegmentation(chamberSeg);
  return chamberSeg;
  
} 
//_____________________________________________________________________________
AliMUONGeometrySegmentation*     
AliMUONTest::CreateTriggerSegmentation(Int_t chamberId, Int_t cathod)
{
// Create Trigger geometry segmentation for given chamber and cathod

    printf("in CreateTriggerSegmentation chamber=%d cathode=%d\n",
	   chamberId,cathod);
  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return 0;
  }  

  AliMUONGeometrySegmentation *chamberSeg = new AliMUONGeometrySegmentation(muon->Chamber(chamberId).GetGeometry());

   //Trigger Segmentations
  AliMUONTriggerSegmentation *trigSegX[9]; 
  AliMUONTriggerSegmentation *trigSegY[9]; 
  for(Int_t i=0; i<9; i++) {
    trigSegX[i] = new AliMUONTriggerSegmentation(0);
    trigSegY[i] = new AliMUONTriggerSegmentation(1);
  }

  AliMUONChamber *iChamber, *iChamber1;
  iChamber1 = &muon->Chamber(10);
  iChamber  = &muon->Chamber(chamberId);
  Float_t zpos1= - iChamber1->Z();  
  Float_t zpos = - iChamber->Z();	     
  Float_t zRatio = zpos / zpos1;

// init
  Int_t nStrip[7]={0,0,0,0,0,0,0};	  
  Float_t stripYsize[7]={0.,0.,0.,0.,0.,0.,0.};
  Float_t stripXsize[7]={0.,0.,0.,0.,0.,0.,0.};

// chamber 8 cathode 0
  for (Int_t i=0; i<7; i++) nStrip[i]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio; 
  trigSegX[8]->Init(0,1,nStrip,stripYsize,stripXsize,0.);  
 
// chamber 8 cathode 1
  for (Int_t i=0; i<6; i++) nStrip[i]=8;
  nStrip[6]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<7; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[8]->Init(0,1,nStrip,stripYsize,stripXsize,0.);  
 
// chamber 7 cathode 0
  for (Int_t i=0; i<6; i++) nStrip[i]=32;
  nStrip[7]=16;  
  for (Int_t i=0; i<6; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[7]->Init(0,2,nStrip,stripYsize,stripXsize,0.);  

// chamber 7 cathode 1
  for (Int_t i=0; i<6; i++) nStrip[i]=8;
  nStrip[6]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<7; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[7]->Init(0,2,nStrip,stripYsize,stripXsize,0.);  
 
// chamber 6 cathode 0
  for (Int_t i=0; i<6; i++) nStrip[i]=32;
  nStrip[7]=16;  
  for (Int_t i=0; i<6; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[6]->Init(0,3,nStrip,stripYsize,stripXsize,0.);  

// chamber 6 cathode 1
  for (Int_t i=0; i<5; i++) nStrip[i]=16;
  nStrip[5]=8;
  nStrip[6]=16;  
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;  
  for (Int_t i=0; i<5; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripXsize[5]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[6]->Init(0,3,nStrip,stripYsize,stripXsize,0.);  

// chamber 5 cathode 0
  nStrip[0]=48;
  for (Int_t i=1; i<3; i++) nStrip[i]=64;  
  for (Int_t i=3; i<6; i++) nStrip[i]=32;
  nStrip[6]=16;
  for (Int_t i=0; i<3; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(0)*zRatio;
  for (Int_t i=3; i<6; i++) 
  stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[5]->Init(0,4,nStrip,stripYsize,stripXsize,AliMUONTriggerConstants::StripLength(0)*zRatio);  

// chamber 5 cathode 1
  for (Int_t i=0; i<5; i++) nStrip[i]=16;
  nStrip[5]=8;  
  nStrip[6]=16;
  stripYsize[0]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  for (Int_t i=1; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<5; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripXsize[5]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[5]->Init(0,4,nStrip,stripYsize,stripXsize,AliMUONTriggerConstants::StripLength(0)*zRatio);  

// chamber 4 cathode 0
  nStrip[0]=0;
  for (Int_t i=1; i<3; i++) nStrip[i]=64;  
  for (Int_t i=3; i<6; i++) nStrip[i]=32;  
  nStrip[6]=16;
  stripYsize[0]=0.;
  for (Int_t i=1; i<3; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(0)*zRatio;
  for (Int_t i=3; i<6; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[0]=0;  
  stripXsize[1]=AliMUONTriggerConstants::StripLength(0)*zRatio;
  for (Int_t i=2; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[4]->Init(0,5,nStrip,stripYsize,stripXsize,0.);  

// chamber 4 cathode 1
  nStrip[0]=0;  
  nStrip[1]=8;  
  for (Int_t i=2; i<5; i++) nStrip[i]=16;
  nStrip[5]=8;
  nStrip[6]=16;  
  stripYsize[0]=0.;  
  for (Int_t i=1; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  stripXsize[0]=0.;
  for (Int_t i=1; i<5; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripXsize[5]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[4]->Init(0,5,nStrip,stripYsize,stripXsize,0.);

// chamber 3 cathode 0
  nStrip[0]=48;
  for (Int_t i=1; i<3; i++) nStrip[i]=64;  
  for (Int_t i=3; i<6; i++) nStrip[i]=32;
  nStrip[6]=16;
  for (Int_t i=0; i<3; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(0)*zRatio;
  for (Int_t i=3; i<6; i++) 
  stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[3]->Init(0,6,nStrip,stripYsize,stripXsize,0.);  

// chamber 3 cathode 1
  for (Int_t i=0; i<5; i++) nStrip[i]=16;
  nStrip[5]=8;  
  nStrip[6]=16;
  stripYsize[0]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  for (Int_t i=1; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<5; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripXsize[5]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[3]->Init(0,6,nStrip,stripYsize,stripXsize,0.);  

// chamber 2 cathode 0
  for (Int_t i=0; i<6; i++) nStrip[i]=32;
  nStrip[7]=16;  
  for (Int_t i=0; i<6; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[2]->Init(0,7,nStrip,stripYsize,stripXsize,0.);  

// chamber 2 cathode 1
  for (Int_t i=0; i<5; i++) nStrip[i]=16;
  nStrip[5]=8;
  nStrip[6]=16;  
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;  
  for (Int_t i=0; i<5; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripXsize[5]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[2]->Init(0,7,nStrip,stripYsize,stripXsize,0.);  

// chamber 1 cathode 0
  for (Int_t i=0; i<6; i++) nStrip[i]=32;
  nStrip[7]=16;  
  for (Int_t i=0; i<6; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(1)*zRatio;
  stripYsize[6]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio;
  trigSegX[1]->Init(0,8,nStrip,stripYsize,stripXsize,0.);  

// chamber 1 cathode 1
  for (Int_t i=0; i<6; i++) nStrip[i]=8;
  nStrip[6]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<7; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[1]->Init(0,8,nStrip,stripYsize,stripXsize,0.);  

// chamber 0 cathode 0
  for (Int_t i=0; i<7; i++) nStrip[i]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  for (Int_t i=0; i<6; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripLength(1)*zRatio;
  stripXsize[6]=AliMUONTriggerConstants::StripLength(2)*zRatio; 
  trigSegX[0]->Init(0,9,nStrip,stripYsize,stripXsize,0.);  
 
// chamber 0 cathode 1
  for (Int_t i=0; i<6; i++) nStrip[i]=8;
  nStrip[6]=16;
  for (Int_t i=0; i<7; i++) 
      stripYsize[i]=AliMUONTriggerConstants::StripLength(3)*zRatio;
  for (Int_t i=0; i<7; i++) 
      stripXsize[i]=AliMUONTriggerConstants::StripWidth(2)*zRatio;
  trigSegY[0]->Init(0,9,nStrip,stripYsize,stripXsize,0.);

  Int_t icount=chamberId-10;  // chamber counter (0 1 2 3)
  Int_t id0=(10+icount+1)*100;

  printf("in CreateTriggerSegmentation here 0 id0=%i \n",id0);  

  for (Int_t i=0; i<9; i++) {      
      if (cathod==0) {	  
	  chamberSeg->Add(id0+i,     trigSegX[i]);
	  chamberSeg->Add(id0+50+i,  trigSegX[i]);
      } else if (cathod==1) {       
	  chamberSeg->Add(id0+i,     trigSegY[i]);
	  chamberSeg->Add(id0+50+i,  trigSegY[i]);
      }
  }

  printf("in CreateTriggerSegmentation here 1\n");  

  if (!id0) {
      AliWarning(Form("Segmentation for chamber %d , cathod %d is not yet defined",chamberId, cathod));
      return 0;      
  }

  DrawSegmentation(chamberSeg);  
  return chamberSeg;
}


//
// public methods
//

//______________________________________________________________________________
AliMUONGeometrySegmentation* 
AliMUONTest::CreateSegmentation(Int_t chamberId, Int_t cath)
{
// Create geometry segmentation for the specified chamber and cathod

  switch (chamberId) {

    // Station1
    case 0: 
    case 1:
        return CreateSt1Segmentation(chamberId, cath);
	break;

    // Station2
    case 2: 
    case 3:
        AliWarning("Not yet implemented");
	return 0;
	break;

    // Slat stations
    case 4: 
    case 5: 
    case 6: 
    case 7: 
    case 8: 
    case 9:
        return CreateSlatSegmentation(chamberId, cath);
	break;
		
    // Trigger stations
    case 10: 
    case 11: 
    case 12: 
    case 13:
        return CreateTriggerSegmentation(chamberId, cath);
	break;

    default:
        AliWarning("Wrong chamber Id");
	return 0;
	break;
  }	
}		
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

//________________________________________________________________________
void AliMUONTest::PrintPadPositionsOld()
{
// Iterate over all old chamber segmentations and prints
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

//______________________________________________________________________________
void AliMUONTest::ForWhole(AliMUONTests testCase)
{
// Perform test for all chambers and first cathod

  AliMUON* muon = (AliMUON*)gAlice->GetModule("MUON");
  if (!muon) {
    AliFatal("MUON detector not defined.");
    return;
  }  

  TStopwatch timer;
  timer.Start();  

  // Loop over chambers
//  for (Int_t iChamber=0; iChamber<AliMUONConstants::NCh(); iChamber++) {
  for (Int_t iChamber=10; iChamber<14; iChamber++) {

    // Loop over cathods
    //for (Int_t cath=0; cath<2; cath++) {
    for (Int_t cath=0; cath<1; cath++) {

      AliMUONGeometrySegmentation* segmentation 
        = CreateSegmentation(iChamber, cath);
	
      if (!segmentation) continue;
      	
      cout << setw(6) << "Pads in chamber " << iChamber 
           << " cathod " << cath << endl;
      cout << "===================================" << endl;  

      ForSegmentation(testCase, segmentation);
           
      //if (testCase == kDrawPads) {
      //}	
    }  
  }     
  timer.Stop();
  timer.Print();
}    

//______________________________________________________________________________
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
    = segmentation->GetGeometry()->GetDetElementStore();
    
  for (Int_t j=0; j<detElements->GetNofEntries(); j++) {
       
    Int_t detElemId = detElements->GetEntry(j)->GetUniqueID();       
    cout << "Detection element id: " << detElemId << endl;
    
    ForDetElement(testCase, detElemId, segmentation);
  }  

  After(testCase);

  timer.Stop();
  timer.Print();
} 
   
//______________________________________________________________________________
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
   
//______________________________________________________________________________
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

//______________________________________________________________________________
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

//______________________________________________________________________________
void AliMUONTest::PrintPad(Int_t& counter,
                           Int_t detElemId, Int_t ix, Int_t iy,
                           AliMUONGeometrySegmentation* segmentation)
{
// Prints global pad positions for the given pad
  
  Float_t x, y, z;
  Bool_t success
    = segmentation->GetPadC(detElemId, ix, iy, x, y, z);
  
  if (!success) return;  

  cout << setw(6) << "counter " << counter++ << "   ";
  cout << "Pad indices:  ( " << detElemId << "; " << ix << ", " << iy << " )  " ;
  cout << "Pad position: ( " << x << ", " << y << ", " << z << " )" << endl;
} 
   
//______________________________________________________________________________
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
   
//______________________________________________________________________________
void AliMUONTest::DrawSegmentation(AliMUONGeometrySegmentation *seg)
{
// TBR

  // Drawing slat504
  Int_t ix, iy, deId;
  Float_t x, y, z;
  Float_t dpx, dpy;
//   TH2F * frame = new TH2F(" "," ",10,-10.,245.,10, -5., 45.);
//   TH2F * frame = new TH2F(" "," ",10,-300.,300.,10, -300., 300.);
  TH2F * frame = new TH2F(" "," ",10,-200.,200.,10, -200., 200.);
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
    = seg->GetGeometry()->GetDetElementStore();
    
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



