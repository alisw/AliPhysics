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

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////

/* $Id$ */

#include "AliMUONFactoryV2.h"
#include "AliRun.h"
#include "AliLog.h"

#include "AliMUON.h"
#include "AliMUONConstants.h"
#include "AliMUONTriggerConstants.h"
#include "AliMUONChamber.h"
#include "AliMUONResponseV0.h"
#include "AliMUONGeometryModule.h"
#include "AliMUONGeometryStore.h"
#include "AliMUONGeometrySegmentation.h"
#include "AliMUONVGeometryDEIndexing.h"
#include "AliMUONSt12QuadrantSegmentation.h"
#include "AliMUONSt345SlatSegmentation.h"
#include "AliMUONTriggerSegmentation.h"
#include "AliMUONResponseTrigger.h"

ClassImp(AliMUONFactoryV2)

//__________________________________________________________________________
  AliMUONFactoryV2::AliMUONFactoryV2(const char* name)
    : TNamed(name, ""),
      fMUON(0),
      fResponse0(0),
      fDESegmentations(0)
{
  // FactoryV2 inherite from Factory for switching in AliMUONv1::Init()
  // to be changed when old segmentation will be removed.

  fDESegmentations = new TObjArray();
  fDESegmentations->SetOwner(kTRUE);
}

//__________________________________________________________________________
  AliMUONFactoryV2::AliMUONFactoryV2()
    : TNamed(),
      fMUON(0),
      fResponse0(0),
      fDESegmentations(0)
{
// Default constructor
}

//__________________________________________________________________________
AliMUONFactoryV2::AliMUONFactoryV2(const AliMUONFactoryV2& rhs)
 : TNamed(rhs)
{
  // Protected copy constructor

  AliFatal("Not implemented.");
}

//__________________________________________________________________________

AliMUONFactoryV2::~AliMUONFactoryV2()
{
// Destructor

  delete fDESegmentations;
}

//__________________________________________________________________________
AliMUONFactoryV2&  AliMUONFactoryV2::operator=(const AliMUONFactoryV2& rhs)
{
  // Protected assignement operator

  if (this == &rhs) return *this;

  AliFatal("Not implemented.");
    
  return *this;  
}    
          
//__________________________________________________________________________
Bool_t AliMUONFactoryV2::IsGeometryDefined(Int_t ichamber)
{
// Return true, if det elements for the chamber with the given ichamber Id
// are defined in geometry (the geometry builder for this chamber was activated)

  if ( ! fMUON ||
       ! fMUON->Chamber(ichamber).GetGeometry() ||
       ! fMUON->Chamber(ichamber).GetGeometry()->GetDEIndexing() ||
       ! fMUON->Chamber(ichamber).GetGeometry()->GetDEIndexing()->GetNofDetElements() )
       
    return kFALSE;
  
  return kTRUE;
}  

//__________________________________________________________________________
void AliMUONFactoryV2::BuildCommon() 
{
  //
  // Construct the default response.
  //

  // Default response: 5 mm of gas
  fResponse0 = new AliMUONResponseV0;
  fResponse0->SetSqrtKx3AndDeriveKx2Kx4(0.7131); // sqrt(0.5085)
  fResponse0->SetSqrtKy3AndDeriveKy2Ky4(0.7642); // sqrt(0.5840)
  fResponse0->SetPitch(AliMUONConstants::Pitch()); // anode-cathode distance
  fResponse0->SetSigmaIntegration(10.);
  fResponse0->SetChargeSlope(10);
  fResponse0->SetChargeSpread(0.18, 0.18);
  fResponse0->SetMaxAdc(4096);
  fResponse0->SetSaturation(3000);
  fResponse0->SetZeroSuppression(6);
}       
        
//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation1() 
{
  //--------------------------------------------------------
  // Configuration for Chamber TC1/2  (Station 1) ----------           
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  // Response for 4 mm of gas (station 1)
  // automatic consistency with width of sensitive medium in CreateGeometry ????
  AliMUONResponseV0* responseSt1 = new AliMUONResponseV0;
  // Mathieson parameters from L.Kharmandarian's thesis, page 190
  responseSt1->SetSqrtKx3AndDeriveKx2Kx4(0.7000); // sqrt(0.4900)
  responseSt1->SetSqrtKy3AndDeriveKy2Ky4(0.7550); // sqrt(0.5700)
  responseSt1->SetPitch(AliMUONConstants::PitchSt1()); // anode-cathode distance
  responseSt1->SetSigmaIntegration(10.);
  // ChargeSlope larger to compensate for the smaller anode-cathode distance
  // and keep the same most probable ADC channel for mip's
  responseSt1->SetChargeSlope(62.5); 
  // assumed proportionality to anode-cathode distance for ChargeSpread
  responseSt1->SetChargeSpread(0.144, 0.144);
  responseSt1->SetMaxAdc(4096);
  responseSt1->SetSaturation(3000);
  responseSt1->SetZeroSuppression(6);

  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt1
    = new AliMUONSt12QuadrantSegmentation(kStation1, kNonBendingPlane);
  
  // Add in the array (for safe deleting)  
  fDESegmentations->Add(bendSt1);  
  fDESegmentations->Add(nonbendSt1);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 0; chamber < 2; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC1/2  (Station 1) ----------           

    // cathode 0
    segmentation[0]->Add(id0,      bendSt1);
    segmentation[0]->Add(id0 +  3, nonbendSt1);
    segmentation[0]->Add(id0 +  2, bendSt1);
    segmentation[0]->Add(id0 +  1, nonbendSt1); 
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      nonbendSt1);
    segmentation[1]->Add(id0 +  3, bendSt1);
    segmentation[1]->Add(id0 +  2, nonbendSt1);
    segmentation[1]->Add(id0 +  1, bendSt1);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);
        
    fMUON->SetResponseModel(chamber, responseSt1); // special response      
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
        
  }
}

//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation2() 
{
  //
  //--------------------------------------------------------
  // Configuration for Chamber TC3/4 (Station 2) -----------
  ///^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  // Quadrant segmentations:
  AliMUONSt12QuadrantSegmentation* bendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kBendingPlane);
  AliMUONSt12QuadrantSegmentation* nonbendSt2
    = new AliMUONSt12QuadrantSegmentation(kStation2, kNonBendingPlane);

  // Add in the array (for safe deleting)  
  fDESegmentations->Add(bendSt2);  
  fDESegmentations->Add(nonbendSt2);  

  AliMUONGeometrySegmentation* segmentation[2];

  for (Int_t chamber = 2; chamber < 4; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC3/4  (Station 2) ----------           

    // cathode 0
    segmentation[0]->Add(id0,      bendSt2);
    segmentation[0]->Add(id0 +  3, nonbendSt2);
    segmentation[0]->Add(id0 +  2, bendSt2);
    segmentation[0]->Add(id0 +  1, nonbendSt2); 
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    segmentation[1]->Add(id0,      nonbendSt2);
    segmentation[1]->Add(id0 +  3, bendSt2);
    segmentation[1]->Add(id0 +  2, nonbendSt2);
    segmentation[1]->Add(id0 +  1, bendSt2);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);
        
    fMUON->SetResponseModel(chamber, fResponse0); // normal response        
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
        
  }
}       
        
        
//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation3() 
{
  //--------------------------------------------------------
  // Configuration for Chamber TC5/6  (Station 3) ----------          
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  AliMUONGeometrySegmentation* segmentation[2];

  //Slats Segmentations
  AliMUONSt345SlatSegmentation *slatsegB[4]; // Types of segmentation for St3
  AliMUONSt345SlatSegmentation *slatsegNB[4]; 
  // Bending

  Int_t ndiv[4] ={ 4, 4, 2, 1};  // densities zones 
  for(Int_t i=0; i<4; i++) {
    slatsegB[i] = new AliMUONSt345SlatSegmentation(1);
    fDESegmentations->Add(slatsegB[i]);  
    slatsegB[i]->SetPadSize(10.,0.5);
    slatsegB[i]->SetPadDivision(ndiv);
    slatsegB[i]->SetId(1); // Id elt ????
    slatsegB[i]->SetDAnod(AliMUONConstants::Pitch());
    slatsegNB[i] = new AliMUONSt345SlatSegmentation(0);
    fDESegmentations->Add(slatsegNB[i]);  
    slatsegNB[i]->SetPadSize(1./1.4,10.); // Nbending
    slatsegNB[i]->SetPadDivision(ndiv);
    slatsegNB[i]->SetId(1);
    slatsegNB[i]->SetDAnod(AliMUONConstants::Pitch());
  }

  // Type 112200 for 500, 501, 508, 509, 510, 517 
  // in Ch5 (similar for Ch6) for the futur official numbering
  // Type 112200 for 503, 504, 505, 555, 554, 553 
  // in Ch5 (similar for Ch6) actual numbering in the code to be changed in jan05
  Int_t n0[4] = { 0, 2, 2, 0 };
  slatsegB[0]->SetPcbBoards(n0);
  slatsegB[0]->Init(0);
  slatsegNB[0]->SetPcbBoards(n0);
  slatsegNB[0]->Init(0);
    
  // Type 122200 for 502, 507, 511, 516 (similar in Ch6) 
  // for future official numbering of ALICE
  // Type 122200 for 502, 506, 556, 552 (similiarin Ch6) 
  // for actual numbering in muon code to be changed in jan05
  Int_t n1[4] = { 0, 1, 3, 0 }; 
  slatsegB[1]->SetPcbBoards(n1);
  slatsegB[1]->Init(0); 
  slatsegNB[1]->SetPcbBoards(n1);
  slatsegNB[1]->Init(0); 
    
  // Type 222000 for 503, 506, 512, 515 (similar in Ch6) 
  // for future official numbering of ALICE
  // Type 222000 for 501, 507, 557, 551 (similiarin Ch6) 
  // for actual numbering in muon code to be changed in jan05
  Int_t n2[4] = { 0, 0, 3, 0 };
  slatsegB[2]->SetPcbBoards(n2);
  slatsegB[2]->Init(0);
  slatsegNB[2]->SetPcbBoards(n2);
  slatsegNB[2]->Init(0);
    
  // Type 220000 for 504, 505, 513, 514 (similar in Ch6) 
  // for future official numbering of ALICE
  // Type 220000 for 500, 508, 558, 550 (similiarin Ch6) 
  // for actual numbering in muon code to be changed in jan05
  Int_t n3[4] = { 0, 0, 2, 0 };
  slatsegB[3]->SetPcbBoards(n3);
  slatsegB[3]->Init(0); 
  slatsegNB[3]->SetPcbBoards(n3);
  slatsegNB[3]->Init(0); 

  for (Int_t chamber = 4; chamber < 6; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    // cathode 0
    // type 220000
    segmentation[0]->Add(id0+14, slatsegB[3]);
    segmentation[0]->Add(id0+ 4, slatsegB[3]);  
    segmentation[0]->Add(id0+13, slatsegB[3]);  
    segmentation[0]->Add(id0+ 5, slatsegB[3]);
    // type 222000
    segmentation[0]->Add(id0+15, slatsegB[2]);
    segmentation[0]->Add(id0+ 3, slatsegB[2]);  
    segmentation[0]->Add(id0+12, slatsegB[2]);  
    segmentation[0]->Add(id0+ 6, slatsegB[2]);
    // type 122200
    segmentation[0]->Add(id0+16, slatsegB[1]);
    segmentation[0]->Add(id0+ 2, slatsegB[1]);  
    segmentation[0]->Add(id0+11, slatsegB[1]);  
    segmentation[0]->Add(id0+ 7, slatsegB[1]);
    // type 112200
    segmentation[0]->Add(id0+17, slatsegB[0]);
    segmentation[0]->Add(id0, slatsegB[0]);  
    segmentation[0]->Add(id0+ 1, slatsegB[0]);  
    segmentation[0]->Add(id0+10, slatsegB[0]);
    segmentation[0]->Add(id0+ 9, slatsegB[0]);     
    segmentation[0]->Add(id0+ 8, slatsegB[0]);
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    // type 220000
    segmentation[1]->Add(id0+14, slatsegNB[3]);
    segmentation[1]->Add(id0+ 4, slatsegNB[3]);  
    segmentation[1]->Add(id0+13, slatsegNB[3]);  
    segmentation[1]->Add(id0+ 5, slatsegNB[3]);
    // type 222000
    segmentation[1]->Add(id0+15, slatsegNB[2]);
    segmentation[1]->Add(id0+ 3, slatsegNB[2]);  
    segmentation[1]->Add(id0+12, slatsegNB[2]);  
    segmentation[1]->Add(id0+ 6, slatsegNB[2]);
    // type 122200
    segmentation[1]->Add(id0+16, slatsegNB[1]);
    segmentation[1]->Add(id0+ 2, slatsegNB[1]);  
    segmentation[1]->Add(id0+11, slatsegNB[1]);  
    segmentation[1]->Add(id0+ 7, slatsegNB[1]);
    // type 112200
    segmentation[1]->Add(id0+17, slatsegNB[0]);
    segmentation[1]->Add(id0, slatsegNB[0]);  
    segmentation[1]->Add(id0+ 1, slatsegNB[0]);  
    segmentation[1]->Add(id0+10, slatsegNB[0]);
    segmentation[1]->Add(id0+ 9, slatsegNB[0]);     
    segmentation[1]->Add(id0+ 8, slatsegNB[0]);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);

    fMUON->SetResponseModel(chamber, fResponse0);      
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}
        
//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation4() 
{
  //--------------------------------------------------------
  // Configuration for Chamber TC7/8  (Station 4) ----------           
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  AliMUONGeometrySegmentation* segmentation[2];

  //Slats Segmentations
  AliMUONSt345SlatSegmentation *slatsegB[7]; // Types of segmentation for St4
  AliMUONSt345SlatSegmentation *slatsegNB[7]; 
  // Bending

  Int_t ndiv[4] ={ 4, 4, 2, 1};  // densities zones 
  for(Int_t i = 0; i < 7; i++) {
    slatsegB[i] = new AliMUONSt345SlatSegmentation(1);
    fDESegmentations->Add(slatsegB[i]);  
    slatsegB[i]->SetPadSize(10.,0.5);
    slatsegB[i]->SetPadDivision(ndiv);
    slatsegB[i]->SetId(1);
    slatsegB[i]->SetDAnod(AliMUONConstants::Pitch());
    slatsegNB[i] = new AliMUONSt345SlatSegmentation(0);
    fDESegmentations->Add(slatsegNB[i]);  
    slatsegNB[i]->SetPadSize(1./1.4,10.); 
    slatsegNB[i]->SetPadDivision(ndiv);
    slatsegNB[i]->SetId(1);
    slatsegNB[i]->SetDAnod(AliMUONConstants::Pitch());
  }

  Int_t n4[4] = { 0, 1, 2, 2 };
  slatsegB[0]->SetPcbBoards(n4);
  slatsegB[0]->Init(0); // 0 detection element id
  slatsegNB[0]->SetPcbBoards(n4);
  slatsegNB[0]->Init(0); // 0 detection element id
    
  // Type 112233 for 701, 712, 714, 725 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 112233 for 705, 707, 755, 757 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  // Type 112233 for 901, 902, 911, 912, 914, 915, 924, 925 in Ch9 
  // (similar for Ch10) for the futur official numbering
  // Type 112233 for 904, 905, 907, 908, 954, 955, 957, 958 in Ch9 
  // (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n5[4] = { 0, 2, 2, 2 };
  slatsegB[1]->SetPcbBoards(n5);
  slatsegB[1]->Init(0); // 0 detection element id
  slatsegNB[1]->SetPcbBoards(n5);
  slatsegNB[1]->Init(0); // 0 detection element id
    
  // Type 112230 for 702, 711, 715, 724 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 112230 for 704, 708, 754, 758 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  Int_t n6[4] = { 0, 2, 2, 1 };
  slatsegB[2]->SetPcbBoards(n6);
  slatsegB[2]->Init(0); // 0 detection element id
  slatsegNB[2]->SetPcbBoards(n6);
  slatsegNB[2]->Init(0); // 0 detection element id
    
  // Type 222330 for 703, 710, 716, 723 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 222330 for 703, 709, 753, 759 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  Int_t n7[4] = { 0, 0, 3, 2 };
  slatsegB[3]->SetPcbBoards(n7);
  slatsegB[3]->Init(0); // 0 detection element id
  slatsegNB[3]->SetPcbBoards(n7);
  slatsegNB[3]->Init(0); // 0 detection element id
    
  // Type 223300 for 704, 709, 717, 722 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 223300 for 702, 710, 752, 760 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  Int_t n8[4] = { 0, 0, 2, 2 };
  slatsegB[4]->SetPcbBoards(n8);
  slatsegB[4]->Init(0); // 0 detection element id
  slatsegNB[4]->SetPcbBoards(n8);
  slatsegNB[4]->Init(0); // 0 detection element id
    
  // Type 333000 for 705, 708, 718, 721 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 333000 for 701, 711, 751, 761 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  // Type 333000 for 906, 907, 919, 920 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 333000 for 900, 912, 950, 962 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n9[4] = { 0, 0, 0, 3 };
  slatsegB[5]->SetPcbBoards(n9);
  slatsegB[5]->Init(0); // 0 detection element id
  slatsegNB[5]->SetPcbBoards(n9);
  slatsegNB[5]->Init(0); // 0 detection element id
    
  // Type 330000 for 706, 707, 719, 720 in Ch7 (similar for Ch8) 
  // for the futur official numbering
  // Type 330000 for 700, 712, 750, 762 in Ch7 (similar for Ch8) 
  // actual numbering in the code to be changed in jan05
  Int_t n10[4] = { 0, 0, 0, 2 };
  slatsegB[6]->SetPcbBoards(n10);
  slatsegB[6]->Init(0); // 0 detection element id
  slatsegNB[6]->SetPcbBoards(n10);
  slatsegNB[6]->Init(0); // 0 detection element id


  for (Int_t chamber = 6; chamber < 8; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
        
    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC6/7  (Station 4) ----------           

    // cathode 0
    // type 122330
    segmentation[0]->Add(id0+13, slatsegB[0]);
    segmentation[0]->Add(id0   , slatsegB[0]);
  
    // type 112233
    segmentation[0]->Add(id0+14, slatsegB[1]);
    segmentation[0]->Add(id0+12, slatsegB[1]);  
    segmentation[0]->Add(id0+25, slatsegB[1]);  
    segmentation[0]->Add(id0+ 1, slatsegB[1]);
   
    // type 112230
    segmentation[0]->Add(id0+15, slatsegB[2]);
    segmentation[0]->Add(id0+11, slatsegB[2]);  
    segmentation[0]->Add(id0+24, slatsegB[2]);  
    segmentation[0]->Add(id0+ 2, slatsegB[2]);

    // type 222330 
    segmentation[0]->Add(id0+16, slatsegB[3]);
    segmentation[0]->Add(id0+10, slatsegB[3]);  
    segmentation[0]->Add(id0+23, slatsegB[3]);
    segmentation[0]->Add(id0+ 3, slatsegB[3]);

    // type 223300 
    segmentation[0]->Add(id0+17, slatsegB[4]);
    segmentation[0]->Add(id0+ 9, slatsegB[4]);  
    segmentation[0]->Add(id0+22, slatsegB[4]);
    segmentation[0]->Add(id0+ 4, slatsegB[4]);

    // type 333000 
    segmentation[0]->Add(id0+18, slatsegB[5]);
    segmentation[0]->Add(id0+ 8, slatsegB[5]);  
    segmentation[0]->Add(id0+21, slatsegB[5]);
    segmentation[0]->Add(id0+ 5, slatsegB[5]);

    // type 330000 
    segmentation[0]->Add(id0+19, slatsegB[6]);
    segmentation[0]->Add(id0+ 7, slatsegB[6]);  
    segmentation[0]->Add(id0+20, slatsegB[6]);
    segmentation[0]->Add(id0+ 6, slatsegB[6]);
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);   

    // cathode 1
    // type 122330
    segmentation[1]->Add(id0+13, slatsegNB[0]);
    segmentation[1]->Add(id0   , slatsegNB[0]);

    // type 112233
    segmentation[1]->Add(id0+14, slatsegNB[1]);
    segmentation[1]->Add(id0+12, slatsegNB[1]);  
    segmentation[1]->Add(id0+25, slatsegNB[1]);  
    segmentation[1]->Add(id0+ 1, slatsegNB[1]);
  
    // type 112230
    segmentation[1]->Add(id0+15, slatsegNB[2]);
    segmentation[1]->Add(id0+11, slatsegNB[2]);  
    segmentation[1]->Add(id0+24, slatsegNB[2]);  
    segmentation[1]->Add(id0+ 2, slatsegNB[2]);

    // type 222330 
    segmentation[1]->Add(id0+16, slatsegNB[3]);
    segmentation[1]->Add(id0+10, slatsegNB[3]);  
    segmentation[1]->Add(id0+23, slatsegNB[3]);
    segmentation[1]->Add(id0+ 3, slatsegNB[3]);

    // type 223300 
    segmentation[1]->Add(id0+17, slatsegNB[4]);
    segmentation[1]->Add(id0+ 9, slatsegNB[4]);  
    segmentation[1]->Add(id0+22, slatsegNB[4]);
    segmentation[1]->Add(id0+ 4, slatsegNB[4]);

    // type 333000 
    segmentation[1]->Add(id0+18, slatsegNB[5]);
    segmentation[1]->Add(id0+ 8, slatsegNB[5]);  
    segmentation[1]->Add(id0+21, slatsegNB[5]);
    segmentation[1]->Add(id0+ 5, slatsegNB[5]);

    // type 330000 
    segmentation[1]->Add(id0+19, slatsegNB[6]);
    segmentation[1]->Add(id0+ 7, slatsegNB[6]);  
    segmentation[1]->Add(id0+20, slatsegNB[6]);
    segmentation[1]->Add(id0+ 6, slatsegNB[6]);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);

    fMUON->SetResponseModel(chamber, fResponse0);      
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}

//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation5() 
{       
  //--------------------------------------------------------
  // Configuration for Chamber TC9/10  (Station 5) ---------           
  //^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  AliMUONGeometrySegmentation* segmentation[2];

  //Slats Segmentations
  AliMUONSt345SlatSegmentation *slatsegB[6]; // Types of segmentation for St5
  AliMUONSt345SlatSegmentation *slatsegNB[6]; 
  // Bending

  Int_t ndiv[4] ={ 4, 4, 2, 1};  // densities zones 
  for(Int_t i = 0; i < 6; i++) {
    slatsegB[i] = new AliMUONSt345SlatSegmentation(1);
    fDESegmentations->Add(slatsegB[i]);  
    slatsegB[i]->SetPadSize(10.,0.5);
    slatsegB[i]->SetPadDivision(ndiv);
    slatsegB[i]->SetId(1);
    slatsegB[i]->SetDAnod(AliMUONConstants::Pitch());
    slatsegNB[i] = new AliMUONSt345SlatSegmentation(0);
    fDESegmentations->Add(slatsegNB[i]);  
    slatsegNB[i]->SetPadSize(1./1.4,10.); 
    slatsegNB[i]->SetPadDivision(ndiv);
    slatsegNB[i]->SetId(1);
    slatsegNB[i]->SetDAnod(AliMUONConstants::Pitch());
  }

  // Type 122330 for 900, 913 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 122330 for 906, 956 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n4[4] = { 0, 1, 2, 2 };
  slatsegB[0]->SetPcbBoards(n4);
  slatsegB[0]->Init(0); // 0 detection element id
  slatsegNB[0]->SetPcbBoards(n4);
  slatsegNB[0]->Init(0); // 0 detection element id
    
  // Type 112233 for 901, 902, 911, 912, 914, 915, 924, 925 in Ch9 
  // (similar for Ch10) for the futur official numbering
  // Type 112233 for 904, 905, 907, 908, 954, 955, 957, 958 in Ch9 
  // (similar for Ch10) actual numbering in the code to be changed in jan05
  Int_t n5[4] = { 0, 2, 2, 2 };
  slatsegB[1]->SetPcbBoards(n5);
  slatsegB[1]->Init(0); // 0 detection element id
  slatsegNB[1]->SetPcbBoards(n5);
  slatsegNB[1]->Init(0); // 0 detection element id

  // Type 333000 for 906, 907, 919, 920 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 333000 for 900, 912, 950, 962 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n9[4] = { 0, 0, 0, 3 };
  slatsegB[2]->SetPcbBoards(n9);
  slatsegB[2]->Init(0); // 0 detection element id
  slatsegNB[2]->SetPcbBoards(n9);
  slatsegNB[2]->Init(0); // 0 detection element id

  // Type 222333 for 903, 910, 916, 923 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 222333 for 903, 909, 953, 959 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n11[4] = { 0, 0, 3, 3 };
  slatsegB[3]->SetPcbBoards(n11);
  slatsegB[3]->Init(0); // 0 detection element id
  slatsegNB[3]->SetPcbBoards(n11);
  slatsegNB[3]->Init(0); // 0 detection element id
  
  // Type 223330 for 904, 909, 917, 922 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 223330 for 902, 910, 952, 960 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n12[4] = { 0, 0, 2, 3 };
  slatsegB[4]->SetPcbBoards(n12);
  slatsegB[4]->Init(0); // 0 detection element id
  slatsegNB[4]->SetPcbBoards(n12);
  slatsegNB[4]->Init(0); // 0 detection element id
    
  // Type 333300 for 905, 908, 918, 921 in Ch9 (similar for Ch10) 
  // for the futur official numbering
  // Type 333300 for 901, 911, 951, 961 in Ch9 (similar for Ch10) 
  // actual numbering in the code to be changed in jan05
  Int_t n13[4] = { 0, 0, 0, 4 };
  slatsegB[5]->SetPcbBoards(n13);
  slatsegB[5]->Init(0); // 0 detection element id
  slatsegNB[5]->SetPcbBoards(n13);
  slatsegNB[5]->Init(0); // 0 detection element id

  for (Int_t chamber = 8; chamber < 10; chamber++) {

    segmentation[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
    segmentation[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());

    // id detection elt for chamber 1
    Int_t id0 = (chamber+1)*100;

    //--------------------------------------------------------
    // Configuration for Chamber TC8/9  (Station 5) ----------           

    // cathode 0
    // type 122330
    segmentation[0]->Add(id0+13, slatsegB[0]);
    segmentation[0]->Add(id0   , slatsegB[0]);
  
    // type 112233
    segmentation[0]->Add(id0+15, slatsegB[1]);
    segmentation[0]->Add(id0+14, slatsegB[1]);
    segmentation[0]->Add(id0+12, slatsegB[1]);  
    segmentation[0]->Add(id0+11, slatsegB[1]);  
    segmentation[0]->Add(id0+24, slatsegB[1]);  
    segmentation[0]->Add(id0+25, slatsegB[1]);  
    segmentation[0]->Add(id0+ 1, slatsegB[1]);
    segmentation[0]->Add(id0+ 2, slatsegB[1]);

    // type 333000 
    segmentation[0]->Add(id0+19, slatsegB[2]);
    segmentation[0]->Add(id0+ 7, slatsegB[2]);  
    segmentation[0]->Add(id0+20, slatsegB[2]);
    segmentation[0]->Add(id0+ 6, slatsegB[2]);
 
    // type 222333 
    segmentation[0]->Add(id0+16, slatsegB[3]);
    segmentation[0]->Add(id0+10, slatsegB[3]);  
    segmentation[0]->Add(id0+23, slatsegB[3]);
    segmentation[0]->Add(id0+ 3, slatsegB[3]);
 
    // type 223330 
    segmentation[0]->Add(id0+17, slatsegB[4]);
    segmentation[0]->Add(id0+ 9, slatsegB[4]);  
    segmentation[0]->Add(id0+22, slatsegB[4]);
    segmentation[0]->Add(id0+ 4, slatsegB[4]);
  
    // type 333300 
    segmentation[0]->Add(id0+18, slatsegB[5]);
    segmentation[0]->Add(id0+ 8, slatsegB[5]);  
    segmentation[0]->Add(id0+21, slatsegB[5]);
    segmentation[0]->Add(id0+ 5, slatsegB[5]);
    fMUON->SetSegmentationModel(chamber, 1, segmentation[0]);

    // cathode 1
    // type 122330
    segmentation[1]->Add(id0+13, slatsegNB[0]);
    segmentation[1]->Add(id0   , slatsegNB[0]);
  
    // type 112233
    segmentation[1]->Add(id0+15, slatsegNB[1]);
    segmentation[1]->Add(id0+14, slatsegNB[1]);
    segmentation[1]->Add(id0+12, slatsegNB[1]);  
    segmentation[1]->Add(id0+11, slatsegNB[1]);  
    segmentation[1]->Add(id0+24, slatsegNB[1]);  
    segmentation[1]->Add(id0+25, slatsegNB[1]);  
    segmentation[1]->Add(id0+ 1, slatsegNB[1]);
    segmentation[1]->Add(id0+ 2, slatsegNB[1]);

    // type 333000 
    segmentation[1]->Add(id0+19 , slatsegNB[2]);
    segmentation[1]->Add(id0+ 7, slatsegNB[2]);  
    segmentation[1]->Add(id0+20, slatsegNB[2]);
    segmentation[1]->Add(id0+ 6, slatsegNB[2]);
 
    // type 222333 
    segmentation[1]->Add(id0+16, slatsegNB[3]);
    segmentation[1]->Add(id0+10, slatsegNB[3]);  
    segmentation[1]->Add(id0+23, slatsegNB[3]);
    segmentation[1]->Add(id0+ 3, slatsegNB[3]);
 
    // type 223330 
    segmentation[1]->Add(id0+17, slatsegNB[4]);
    segmentation[1]->Add(id0+ 9, slatsegNB[4]);  
    segmentation[1]->Add(id0+22, slatsegNB[4]);
    segmentation[1]->Add(id0+ 4, slatsegNB[4]);
  
    // type 333300 
    segmentation[1]->Add(id0+18, slatsegNB[5]);
    segmentation[1]->Add(id0+ 8, slatsegNB[5]);  
    segmentation[1]->Add(id0+21, slatsegNB[5]);
    segmentation[1]->Add(id0+ 5, slatsegNB[5]);
    fMUON->SetSegmentationModel(chamber, 2, segmentation[1]);
        
    fMUON->SetResponseModel(chamber, fResponse0);           
    fMUON->Chamber(chamber).SetChargeCorrel(0.11); // 11% charge spread
  }
}

//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation6() 
{       
 // Create Trigger geometry segmentation for given chamber and cathod

 
    AliMUONGeometrySegmentation *chamberSeg[2];
// Cluster-size off
        AliMUONResponseTrigger* responseTrigger0 =  new AliMUONResponseTrigger;
// Cluster-size on  
//  AliMUONResponseTriggerV1* responseTrigger0 =  new AliMUONResponseTriggerV1;

    for (Int_t chamber = 10; chamber < 14; chamber++) {

      //Trigger Segmentation
      AliMUONTriggerSegmentation *trigSegX[9]; 
      AliMUONTriggerSegmentation *trigSegY[9]; 
      for(Int_t i=0; i<9; i++) {
        trigSegX[i] = new AliMUONTriggerSegmentation(1);
        trigSegY[i] = new AliMUONTriggerSegmentation(0);
        fDESegmentations->Add(trigSegX[i]);  
        fDESegmentations->Add(trigSegY[i]);  
        trigSegX[i]->SetLineNumber(9-i);    
        trigSegY[i]->SetLineNumber(9-i);    
      }

      AliMUONChamber *iChamber, *iChamber1;
      iChamber1 = &fMUON->Chamber(10);
      iChamber  = &fMUON->Chamber(chamber);
//cp      Float_t zpos1= - iChamber1->Z();  
      Float_t zpos1= iChamber1->Z();  
//cp      Float_t zpos = - iChamber->Z();        
      Float_t zpos = iChamber->Z();        
      Float_t zRatio = zpos / zpos1;

      // init
      Float_t stripWidth[3]={0.,0.,0.};     // 1.0625 2.125 4.25
      Float_t stripLength[4]={0.,0.,0.,0.}; // 17. 34. 51. 68.
      for (Int_t i=0; i<3; i++) 
        stripWidth[i]=AliMUONTriggerConstants::StripWidth(i)*zRatio;
      for (Int_t i=0; i<4; i++) 
        stripLength[i]=AliMUONTriggerConstants::StripLength(i)*zRatio;
      Int_t nStrip[7]={0,0,0,0,0,0,0};    
      Float_t stripYsize[7]={0.,0.,0.,0.,0.,0.,0.};
      Float_t stripXsize[7]={0.,0.,0.,0.,0.,0.,0.};

      // chamber 8 0 cathode 0
      for (Int_t i=0; i<7; i++) nStrip[i]=16;
      for (Int_t i=0; i<7; i++) stripYsize[i]=stripWidth[2];
      for (Int_t i=0; i<6; i++) stripXsize[i]=stripLength[1];
      stripXsize[6]=stripLength[2];
      trigSegX[8]->Init(0,nStrip,stripYsize,stripXsize,0.); 
      trigSegX[0]->Init(0,nStrip,stripYsize,stripXsize,0.); 

      // chamber 8 7 1 0 cathode 1
      for (Int_t i=0; i<6; i++) nStrip[i]=8;
      nStrip[6]=16;
      for (Int_t i=0; i<7; i++) stripYsize[i]=stripLength[3];  
      for (Int_t i=0; i<7; i++) stripXsize[i]=stripWidth[2];
      trigSegY[8]->Init(0,nStrip,stripYsize,stripXsize,0.);  
      trigSegY[7]->Init(0,nStrip,stripYsize,stripXsize,0.);
      trigSegY[1]->Init(0,nStrip,stripYsize,stripXsize,0.);
      trigSegY[0]->Init(0,nStrip,stripYsize,stripXsize,0.);
 
      // chamber 7 6 2 1 cathode 0
      for (Int_t i=0; i<6; i++) nStrip[i]=32;
      nStrip[6]=16;  
      for (Int_t i=0; i<6; i++) stripYsize[i]=stripWidth[1];
      stripYsize[6]=stripWidth[2];
      for (Int_t i=0; i<6; i++) stripXsize[i]=stripLength[1];
      stripXsize[6]=stripLength[2];
      trigSegX[7]->Init(0,nStrip,stripYsize,stripXsize,0.);  
      trigSegX[6]->Init(0,nStrip,stripYsize,stripXsize,0.);
      trigSegX[2]->Init(0,nStrip,stripYsize,stripXsize,0.);  
      trigSegX[1]->Init(0,nStrip,stripYsize,stripXsize,0.);

      // chamber 6 2 cathode 1
      for (Int_t i=0; i<5; i++) nStrip[i]=16;
      for (Int_t i=5; i<6; i++) nStrip[i]=8;
      nStrip[6]=16;
      for (Int_t i=0; i<7; i++) stripYsize[i]=stripLength[3];
      for (Int_t i=0; i<5; i++) stripXsize[i]=stripWidth[1];
      for (Int_t i=5; i<7; i++) stripXsize[i]=stripWidth[2];
      trigSegY[6]->Init(0,nStrip,stripYsize,stripXsize,0.);  
      trigSegY[2]->Init(0,nStrip,stripYsize,stripXsize,0.);  

      // chamber 5 3 cathode 0
      nStrip[0]=48;
      for (Int_t i=1; i<3; i++) nStrip[i]=64;
      for (Int_t i=3; i<6; i++) nStrip[i]=32;
      nStrip[6]=16;  
      for (Int_t i=0; i<3; i++) stripYsize[i]=stripWidth[0];
      for (Int_t i=3; i<6; i++) stripYsize[i]=stripWidth[1];
      stripYsize[6]=stripWidth[2];
      for (Int_t i=0; i<6; i++) stripXsize[i]=stripLength[1];
      stripXsize[6]=stripLength[2];
      trigSegX[5]->Init(0,nStrip,stripYsize,stripXsize,stripLength[0]);  
      trigSegX[3]->Init(0,nStrip,stripYsize,stripXsize,0.);

      // chamber 5 3 cathode 1
      for (Int_t i=0; i<5; i++) nStrip[i]=16;
      for (Int_t i=5; i<6; i++) nStrip[5]=8;  
      nStrip[6]=16;  
      stripYsize[0]=stripLength[2];
      for (Int_t i=1; i<7; i++) stripYsize[i]=stripLength[3];
      for (Int_t i=0; i<5; i++) stripXsize[i]=stripWidth[1];
      for (Int_t i=5; i<7; i++) stripXsize[i]=stripWidth[2];
      trigSegY[5]->Init(0,nStrip,stripYsize,stripXsize,stripLength[0]);  
      trigSegY[3]->Init(0,nStrip,stripYsize,stripXsize,0.);

      // chamber 4 cathode 0
      nStrip[0]=0;
      for (Int_t i=1; i<3; i++) nStrip[i]=64;  
      for (Int_t i=3; i<6; i++) nStrip[i]=32;  
      nStrip[6]=16;
      stripYsize[0]=0.;
      for (Int_t i=1; i<3; i++) stripYsize[i]=stripWidth[0];
      for (Int_t i=3; i<6; i++) stripYsize[i]=stripWidth[1];
      stripYsize[6]=stripWidth[2];
      stripXsize[0]=0;  
      stripXsize[1]=stripLength[0];  
      for (Int_t i=2; i<6; i++) stripXsize[i]=stripLength[1];
      stripXsize[6]=stripLength[2];
      trigSegX[4]->Init(0,nStrip,stripYsize,stripXsize,0.);  

      // chamber 4 cathode 1
      nStrip[0]=0;  
      nStrip[1]=8;  
      for (Int_t i=2; i<5; i++) nStrip[i]=16;
      for (Int_t i=5; i<6; i++) nStrip[i]=8;
      nStrip[6]=16;
      stripYsize[0]=0.;  
      for (Int_t i=1; i<7; i++) stripYsize[i]=stripLength[3];
      stripXsize[0]=0.;
      for (Int_t i=1; i<5; i++) stripXsize[i]=stripWidth[1];
      for (Int_t i=5; i<7; i++) stripXsize[i]=stripWidth[2];
      trigSegY[4]->Init(0,nStrip,stripYsize,stripXsize,0.);

      chamberSeg[0] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());
      chamberSeg[1] = new AliMUONGeometrySegmentation(fMUON->Chamber(chamber).GetGeometry());

      Int_t icount=chamber-10;  // chamber counter (0 1 2 3)
      Int_t id0=(10+icount+1)*100;

      //  printf("in CreateTriggerSegmentation here 0 id0=%i \n",id0);  

      chamberSeg[0]->Add(id0+0,      trigSegX[4]);
      chamberSeg[0]->Add(id0+1,      trigSegX[5]);
      chamberSeg[0]->Add(id0+2,      trigSegX[6]);
      chamberSeg[0]->Add(id0+3,      trigSegX[7]);
      chamberSeg[0]->Add(id0+4,      trigSegX[8]);
      chamberSeg[0]->Add(id0+5,      trigSegX[8]);
      chamberSeg[0]->Add(id0+6,      trigSegX[7]);
      chamberSeg[0]->Add(id0+7,      trigSegX[6]);
      chamberSeg[0]->Add(id0+8,      trigSegX[5]);
      chamberSeg[0]->Add(id0+9,      trigSegX[4]);
      chamberSeg[0]->Add(id0+10,     trigSegX[3]);
      chamberSeg[0]->Add(id0+11,     trigSegX[2]);
      chamberSeg[0]->Add(id0+12,     trigSegX[1]);
      chamberSeg[0]->Add(id0+13,     trigSegX[0]);
      chamberSeg[0]->Add(id0+14,     trigSegX[0]);
      chamberSeg[0]->Add(id0+15,     trigSegX[1]);
      chamberSeg[0]->Add(id0+16,     trigSegX[2]);
      chamberSeg[0]->Add(id0+17,     trigSegX[3]);

      chamberSeg[1]->Add(id0+0,      trigSegY[4]);
      chamberSeg[1]->Add(id0+1,      trigSegY[5]);
      chamberSeg[1]->Add(id0+2,      trigSegY[6]);
      chamberSeg[1]->Add(id0+3,      trigSegY[7]);
      chamberSeg[1]->Add(id0+4,      trigSegY[8]);
      chamberSeg[1]->Add(id0+5,      trigSegY[8]);
      chamberSeg[1]->Add(id0+6,      trigSegY[7]);
      chamberSeg[1]->Add(id0+7,      trigSegY[6]);
      chamberSeg[1]->Add(id0+8,      trigSegY[5]);
      chamberSeg[1]->Add(id0+9,      trigSegY[4]);
      chamberSeg[1]->Add(id0+10,     trigSegY[3]);
      chamberSeg[1]->Add(id0+11,     trigSegY[2]);
      chamberSeg[1]->Add(id0+12,     trigSegY[1]);
      chamberSeg[1]->Add(id0+13,     trigSegY[0]);
      chamberSeg[1]->Add(id0+14,     trigSegY[0]);
      chamberSeg[1]->Add(id0+15,     trigSegY[1]);
      chamberSeg[1]->Add(id0+16,     trigSegY[2]);
      chamberSeg[1]->Add(id0+17,     trigSegY[3]);

      fMUON->SetSegmentationModel(chamber, 1, chamberSeg[0]);
      fMUON->SetSegmentationModel(chamber, 2, chamberSeg[1]);

      fMUON->SetResponseModel(chamber, responseTrigger0);      
      fMUON->Chamber(chamber).SetChargeCorrel(0); // same charge on cathodes

  
      //  printf("in CreateTriggerSegmentation here 1\n");  

      if (!id0) {
        AliWarning(Form("Segmentation for chamber %d is not yet defined",chamber));
        return ;      
      }
    }
}       
//__________________________________________________________________________
void AliMUONFactoryV2::Build(AliMUON* where, const char* what) 
{
  //
  // Construct MUON from chambers, segmentation and responses
  //

  fMUON = where;
  char tmp[20];
  strcpy(tmp, what);

  if (strcmp(tmp, "default")==0) {
    // Set default parameters
    fMUON->SetIshunt(0);
    fMUON->SetMaxStepGas(0.1);
    fMUON->SetMaxStepAlu(0.1);

    // Build stations
    BuildCommon();
    if (IsGeometryDefined(0))  BuildStation1();
    if (IsGeometryDefined(2))  BuildStation2();
    if (IsGeometryDefined(4))  BuildStation3();
    if (IsGeometryDefined(6))  BuildStation4();
    if (IsGeometryDefined(8))  BuildStation5();
    if (IsGeometryDefined(10)) BuildStation6();
  } 
  else
    AliDebug(0,"Non default version of MUON selected. You have to construct yourself the MUON elements !!");
}

//__________________________________________________________________________
void AliMUONFactoryV2::BuildStation(AliMUON* where, Int_t stationNumber) 
{
  //
  // Construct MUON from chambers, segmentation and responses
  //
  // Version 0
  //
  // Then chose for each chamber (chamber plane) the segmentation
  // and response model.
  // They should be equal for the two chambers of each station. In a future
  // version this will be enforced.
  //

  fMUON = where;
  if (!fResponse0) BuildCommon(); 
    
  switch (stationNumber) {    
  case 1:  BuildStation1(); break;
  case 2:  BuildStation2(); break;
  case 3:  BuildStation3(); break;
  case 4:  BuildStation4(); break;
  case 5:  BuildStation5(); break;
  case 6:  BuildStation6(); break;
    
  default: AliFatal("Wrong station number");
  }  
}         
