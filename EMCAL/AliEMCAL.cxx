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

//_________________________________________________________________________
// Base Class for EMCAL description:
// This class contains material definitions    
// for the EMCAL - It does not place the detector in Alice
//*-- Author: Yves Schutz (SUBATECH) 
//
//*-- Additional Contributions: Sahal Yacoob (LBNL/UCT)
//
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
class TFile;
#include <TFolder.h> 
#include <TTree.h>
#include <TVirtualMC.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliEMCAL.h"
#include "AliEMCALGetter.h"
#include "AliRun.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALDigitizer.h"
#include "AliAltroBuffer.h"

ClassImp(AliEMCAL)
//____________________________________________________________________________
AliEMCAL::AliEMCAL():AliDetector()
{
  // Default ctor 
  fName="EMCAL";
  //  fGeom = 0 ; 
  fRan    = new TRandom(123456) ; 
}

//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title): AliDetector(name,title)
{
  //   ctor : title is used to identify the layout
  // fGeom = 0;
  fRan    = new TRandom(123456) ; 
}

//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{

}

//____________________________________________________________________________
void AliEMCAL::Copy(AliEMCAL & emcal)
{
  TObject::Copy(emcal) ; 
}

//____________________________________________________________________________
AliDigitizer* AliEMCAL::CreateDigitizer(AliRunDigitizer* manager) const
{
  return new AliEMCALDigitizer(manager);
}

//____________________________________________________________________________
void AliEMCAL::CreateMaterials()
{
  // Definitions of materials to build EMCAL and associated tracking media.
  // media number in idtmed are 1599 to 1698.

  // --- Air ---               
  Float_t aAir[4]={12.0107,14.0067,15.9994,39.948};
  Float_t zAir[4]={6.,7.,8.,18.};
  Float_t wAir[4]={0.000124,0.755267,0.231781,0.012827};
  Float_t dAir = 1.20479E-3;
  AliMixture(0, "Air$", aAir, zAir, dAir, 4, wAir) ;

  // --- Lead ---                                                                     
  AliMaterial(1, "Pb$", 207.2, 82, 11.35, 0.56, 0., 0, 0) ;


  // --- The polysterene scintillator (CH) ---
  Float_t aP[2] = {12.011, 1.00794} ;
  Float_t zP[2] = {6.0, 1.0} ;
  Float_t wP[2] = {1.0, 1.0} ;
  Float_t dP = 1.032 ;

  AliMixture(2, "Polystyrene$", aP, zP, dP, -2, wP) ;

  // --- Aluminium ---
  AliMaterial(3, "Al$", 26.98, 13., 2.7, 8.9, 999., 0, 0) ;
  // ---         Absorption length is ignored ^

  // DEFINITION OF THE TRACKING MEDIA

  // for EMCAL: idtmed[1599->1698] equivalent to fIdtmed[0->100]
  Int_t * idtmed = fIdtmed->GetArray() - 1599 ; 
  Int_t   isxfld = gAlice->Field()->Integ() ;
  Float_t sxmgmx = gAlice->Field()->Max() ;

  // Air                                                                         -> idtmed[1599]
 AliMedium(0, "Air          $", 0, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // The Lead                                                                      -> idtmed[1600]
 
  AliMedium(1, "Lead      $", 1, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[1601]
  AliMedium(2, "CPV scint.   $", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, 0.1, 0.001, 0.001, 0, 0) ;

  // Various Aluminium parts made of Al                                            -> idtmed[1602]
  AliMedium(3, "Al parts     $", 3, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;


// --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in Lead 
  gMC->Gstpar(idtmed[1600],"CUTGAM",0.00008) ;
  gMC->Gstpar(idtmed[1600],"CUTELE",0.001) ;
  gMC->Gstpar(idtmed[1600],"BCUTE",0.0001) ;

  // --- Generate explicitly delta rays in Lead ---
  gMC->Gstpar(idtmed[1600], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1600], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1600], "DCUTE",0.00001) ;
  gMC->Gstpar(idtmed[1600], "DCUTM",0.00001) ;

// --- in aluminium parts ---
  gMC->Gstpar(idtmed[1602], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1602], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1602], "DCUTE",0.00001) ;
  gMC->Gstpar(idtmed[1602], "DCUTM",0.00001) ;

// --- and finally thresholds for photons and electrons in the scintillator ---
  gMC->Gstpar(idtmed[1601],"CUTGAM",0.00008) ;
  gMC->Gstpar(idtmed[1601],"CUTELE",0.001) ;
  gMC->Gstpar(idtmed[1601],"BCUTE",0.0001) ;

  //set constants for Birk's Law implentation
  fBirkC0 =  1;
  fBirkC1 =  0.013/dP;
  fBirkC2 =  9.6e-6/(dP * dP);


}
      
//____________________________________________________________________________
void AliEMCAL::Digits2Raw()
{
// convert digits of the current event to raw data

  // get the digits
  AliEMCALGetter * gime = AliEMCALGetter::Instance(AliRunLoader::GetGAliceName()) ; 
  if (!gime) {
    Error("Digits2Raw", "EMCAL Getter not instantiated") ;
    return ; 
  }
  gime->Event(gime->EventNumber(), "D") ; 
  TClonesArray* digits = gime->Digits() ;

  if (!digits) {
    Error("Digits2Raw", "no digits found !");
    return;
  }

  // get the geometry
  AliEMCALGeometry* geom = gime->EMCALGeometry();
  if (!geom) {
    Error("Digits2Raw", "no geometry found !");
    return;
  }

  // some digitization constants
  const Int_t    kDDLOffset = 0x800;
  const Double_t kTimeMax = 1.28E-5;
  const Int_t    kTimeBins = 256;
  const Double_t kTimePeak = 2.0E-6;
  const Double_t kTimeRes = 1.5E-6;
  const Int_t    kThreshold = 3;
  const Int_t    kHighGainFactor = 40;
  const Int_t    kHighGainOffset = 0x200;
  // PHOS has 4 DDL per module; I assume therefore that kChannelsperDDL=896+1 EMCAL channel go to one DDL
  const Int_t    kChannelsperDDL = 897 ; 

  AliAltroBuffer* buffer = NULL;
  Int_t prevDDL = -1;
  Int_t adcValuesLow[kTimeBins];
  Int_t adcValuesHigh[kTimeBins];

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliEMCALDigit* digit = gime->Digit(iDigit);
    if (digit->GetAmp() < kThreshold) 
      continue;
    Int_t iDDL = digit->GetId() / kChannelsperDDL ;
    // for each DDL id is numbered from 1 to  kChannelsperDDL -1 
    Int_t idDDL = digit->GetId() - iDDL * ( kChannelsperDDL - 1 ) ;  
    // new DDL
    if (iDDL != prevDDL) {
      // write real header and close previous file
      if (buffer) {
	buffer->Flush();
	buffer->WriteDataHeader(kFALSE, kFALSE);
	delete buffer;
      }

      // open new file and write dummy header
      TString fileName("EMCAL_") ;
      fileName += (iDDL + kDDLOffset) ; 
      fileName += ".ddl" ; 
      buffer = new AliAltroBuffer(fileName.Data(), 1);
      buffer->WriteDataHeader(kTRUE, kFALSE);  //Dummy;

      prevDDL = iDDL;
    }

    // out of time range signal (?)
    if (digit->GetTime() > kTimeMax) {
      buffer->FillBuffer(digit->GetAmp());
      buffer->FillBuffer(kTimeBins);  // time bin
      buffer->FillBuffer(3);          // bunch length
      buffer->WriteTrailer(3, idDDL, 0, 0);  // trailer
      
    // simulate linear rise and gaussian decay of signal
    } else {
      Bool_t highGain = kFALSE;

      // fill time bin values :
      // 1. the signal starts at the time given by the digit
      // 2. the rise is linear and the maximum is reached kTimePeak after start
      // 3. the decay is gaussian with a sigma of kTimeRes
      // 4. the signal is binned into kTimeBins bins 
      for (Int_t iTime = 0; iTime < kTimeBins; iTime++) {
	Double_t time = iTime * kTimeMax/kTimeBins;
	Int_t signal = 0;
	if (time < digit->GetTime() + kTimePeak) {	// signal is rising
	  signal = static_cast<Int_t>((fRan->Rndm() + digit->GetAmp()) * 
				      (time - digit->GetTime() / kTimePeak) + 0.5);
	} else {                                        // signal is decaying
	  signal = static_cast<Int_t>((fRan->Rndm() + digit->GetAmp()) * 
				      TMath::Gaus(time, digit->GetTime() + kTimePeak, kTimeRes) + 0.5);
	}
	if (signal < 0) 
	  signal = 0;
	adcValuesLow[iTime] = signal;
	if (signal > 0x3FF) // larger than 10 bits 
	  adcValuesLow[iTime] = 0x3FF;
	adcValuesHigh[iTime] = static_cast<Int_t>(0.5 + (signal / kHighGainFactor));
	if (adcValuesHigh[iTime] > 0) 
	  highGain = kTRUE;
      }

      // write low and eventually high gain channel
      buffer->WriteChannel(idDDL, 0, 0, 
			   kTimeBins, adcValuesLow, kThreshold);
      if (highGain) {
	buffer->WriteChannel(idDDL, 0, kHighGainOffset, 
			     kTimeBins, adcValuesHigh, 1);
      }
    }
  }

  // write real header and close last file
  if (buffer) {
    buffer->Flush();
    buffer->WriteDataHeader(kFALSE, kFALSE);
    delete buffer;
  }

  gime->EmcalLoader()->UnloadDigits();
}

//____________________________________________________________________________
void AliEMCAL::Hits2SDigits()  
{ 
// create summable digits

  AliEMCALSDigitizer* emcalDigitizer = 
    new AliEMCALSDigitizer(fLoader->GetRunLoader()->GetFileName().Data()) ;
  emcalDigitizer->SetEventRange(0, -1) ; // do all the events
  emcalDigitizer->ExecuteTask() ;
}

//____________________________________________________________________________
AliLoader* AliEMCAL::MakeLoader(const char* topfoldername)
{
//different behaviour than standard (singleton getter)
// --> to be discussed and made eventually coherent
 fLoader = new AliEMCALLoader(GetName(),topfoldername);
 return fLoader;
}

//____________________________________________________________________________
void AliEMCAL::SetTreeAddress()
{ 
  // Linking Hits in Tree to Hits array
  TBranch *branch;
  char branchname[20];
  sprintf(branchname,"%s",GetName());
  
  // Branch address for hit tree
  TTree *treeH = TreeH();
  if (treeH) {
    branch = treeH->GetBranch(branchname);
    if (branch) 
      { 
	if (fHits == 0x0) 
	  fHits= new TClonesArray("AliEMCALHit",1000);
	branch->SetAddress(&fHits);
      }
    else
      {
	Warning("SetTreeAddress","<%s> Failed",GetName());
      }
  }
}





