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
#include <TH1F.h> 
#include <TF1.h> 
#include <TRandom.h> 

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
Double_t AliEMCAL::fgCapa        = 1.;        // 1pF 
Int_t    AliEMCAL::fgOrder       = 2 ;
Double_t AliEMCAL::fgTimeMax     = 2.56E-5 ;  // each sample is over 100 ns fTimeMax/fTimeBins
Double_t AliEMCAL::fgTimePeak    = 4.1E-6 ;   // 4 micro seconds
Double_t AliEMCAL::fgTimeTrigger = 100E-9 ;      // 100ns, just for a reference
 
//____________________________________________________________________________
AliEMCAL::AliEMCAL():AliDetector()
{
  // Default ctor 
  fName = "EMCAL" ;
}

//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title): AliDetector(name,title)
{
  //   ctor : title is used to identify the layout

  fHighCharge        = 8.2 ;          // adjusted for a high gain range of 5.12 GeV (10 bits)
  fHighGain          = 6.64 ; 
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits) 
  fLowGainOffset     = 1 ;            // offset added to the module id to distinguish high and low gain data
}

//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{

}

//____________________________________________________________________________
void AliEMCAL::Copy(AliEMCAL & emcal)
{
  TObject::Copy(emcal) ; 
  emcal.fHighCharge        = fHighCharge ;
  emcal.fHighGain          = fHighGain ; 
  emcal.fHighLowGainFactor = fHighLowGainFactor ;  
  emcal.fLowGainOffset     = fLowGainOffset;   
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
  AliEMCALLoader * loader = dynamic_cast<AliEMCALLoader*>(fLoader) ; 

  // get the digits
  loader->LoadDigits();
  TClonesArray* digits = loader->Digits() ;

  if (!digits) {
    Error("Digits2Raw", "no digits found !");
    return;
  }

  // get the digitizer 
  loader->LoadDigitizer();
  AliEMCALDigitizer * digitizer = dynamic_cast<AliEMCALDigitizer *>(loader->Digitizer())  ; 
  
  // get the geometry
  AliEMCALGeometry* geom = GetGeometry();
  if (!geom) {
    Error("Digits2Raw", "no geometry found !");
    return;
  }

  // some digitization constants
  const Int_t    kDDLOffset = 0x800;
  const Int_t    kThreshold = 1;
  const Int_t    kChannelsperDDL = 897 ; 
  AliAltroBuffer* buffer = NULL;
  Int_t prevDDL = -1;
  Int_t adcValuesLow[fkTimeBins];
  Int_t adcValuesHigh[fkTimeBins];
  
  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliEMCALDigit* digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit)) ;
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
    if (digit->GetTimeR() > GetRawFormatTimeMax() ) {
      buffer->FillBuffer(digit->GetAmp());
      buffer->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffer->FillBuffer(3);          // bunch length
      buffer->WriteTrailer(3, idDDL, 0, 0);  // trailer

      // calculate the time response function
    } else {
      Double_t energy = 0 ;  
      energy = digit->GetAmp() * digitizer->GetECAchannel() + digitizer->GetECApedestal() ; 
      
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), energy, adcValuesHigh, adcValuesLow) ; 
      
      if (lowgain) 
	buffer->WriteChannel(iDDL, 0, fLowGainOffset, 
			     GetRawFormatTimeBins(), adcValuesLow, kThreshold);
      else 
	buffer->WriteChannel(iDDL, 0, 0, 
			     GetRawFormatTimeBins(), adcValuesHigh, kThreshold);
      
    }
  }
  
  // write real header and close last file
  if (buffer) {
    buffer->Flush();
    buffer->WriteDataHeader(kFALSE, kFALSE);
    delete buffer;
  }

  loader->UnloadDigits();
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

//__________________________________________________________________
Double_t AliEMCAL::RawResponseFunction(Double_t *x, Double_t *par)
{
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  // v(t) = n**n * Q * A**n / C *(t/tp)**n * exp(-n * t/tp) with 
  // tp : peaking time par[0]
  // n  : order of the function
  // C  : integrating capacitor in the preamplifier
  // A  : open loop gain of the preamplifier
  // Q  : the total APD charge to be measured Q = C * energy
  
  Double_t signal ;
  Double_t xx = x[0] - ( fgTimeTrigger + par[3] ) ; 
  
  if (xx < 0 || xx > fgTimeMax) 
    signal = 0. ;  
  else { 
    Double_t fac = par[0] * TMath::Power(fgOrder, fgOrder) * TMath::Power(par[1], fgOrder) / fgCapa ; 
    signal = fac * par[2] * TMath::Power(xx / fgTimePeak, fgOrder) * TMath::Exp(-fgOrder * (xx / fgTimePeak)) ; 
  }
  return signal ;  
}

//__________________________________________________________________
Double_t AliEMCAL::RawResponseFunctionMax(Double_t charge, Double_t gain) 
{
  return ( charge * TMath::Power(fgOrder, fgOrder) * TMath::Power(gain, fgOrder) 
     / ( fgCapa * TMath::Exp(fgOrder) ) );  

}
//__________________________________________________________________
Bool_t AliEMCAL::RawSampledResponse(
const Double_t dtime, const Double_t damp, Int_t * adcH, Int_t * adcL) const 
{
  // for a start time dtime and an amplitude damp given by digit, 
  // calculates the raw sampled response AliEMCAL::RawResponseFunction

  const Int_t kRawSignalOverflow = 0x3FF ; 
  Bool_t lowGain = kFALSE ; 

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    signalF.SetParameter(0, GetRawFormatHighCharge() ) ; 
    signalF.SetParameter(1, GetRawFormatHighGain() ) ; 
    signalF.SetParameter(2, damp) ; 
    signalF.SetParameter(3, dtime) ; 
    Double_t time = iTime * GetRawFormatTimeMax() / GetRawFormatTimeBins() ;
    Double_t signal = signalF.Eval(time) ;     
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow ){  // larger than 10 bits 
      signal = kRawSignalOverflow ;
      lowGain = kTRUE ; 
    }
    adcH[iTime] =  static_cast<Int_t>(signal + 0.5) ;

    signalF.SetParameter(0, GetRawFormatLowCharge() ) ;     
    signalF.SetParameter(1, GetRawFormatLowGain() ) ; 
    signal = signalF.Eval(time) ;  
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow)  // larger than 10 bits 
      signal = kRawSignalOverflow ;
    adcL[iTime] = static_cast<Int_t>(0.5 + signal ) ; 

  }
  return lowGain ; 
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





