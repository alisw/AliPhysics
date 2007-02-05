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
/* History of cvs commits:
 *
 * $Log$
 * Revision 1.49  2007/01/22 17:29:12  pavlinov
 * EMCAL geometry can be created independently form anything now
 *
 * Revision 1.48  2006/12/19 02:34:13  pavlinov
 * clean up the EMCAL name scheme : super module -> module -> tower (or cell)
 *
 * Revision 1.47  2006/12/05 17:12:03  gustavo
 * Updated AliEMCAL::Digits2Raw, reads first provisional RCU mapping files to make Raw data with new AliCaloAltroMapping and AliCaloRawStream
 *
 *
 */
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
#include <TGraph.h> 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliMagF.h"
#include "AliEMCAL.h"
#include "AliRun.h"
#include "AliEMCALLoader.h"
#include "AliEMCALSDigitizer.h"
#include "AliEMCALDigitizer.h"
#include "AliEMCALDigit.h"
//#include "AliAltroMapping.h"
#include "AliCaloAltroMapping.h"
#include "AliAltroBuffer.h"
#include "AliRawReader.h"
#include "AliCaloRawStream.h"
#include "AliDAQ.h"

ClassImp(AliEMCAL)
Double_t AliEMCAL::fgCapa        = 1.;        // 1pF 
Int_t    AliEMCAL::fgOrder       = 2 ;
Double_t AliEMCAL::fgTimeMax     = 2.56E-5 ;  // each sample is over 100 ns fTimeMax/fTimeBins
Double_t AliEMCAL::fgTimePeak    = 4.1E-6 ;   // 4 micro seconds
Double_t AliEMCAL::fgTimeTrigger = 100E-9 ;      // 100ns, just for a reference
// some digitization constants
Int_t    AliEMCAL::fgThreshold = 1;
Int_t    AliEMCAL::fgDDLPerSuperModule = 2;  // 2 ddls per SuperModule

//____________________________________________________________________________
AliEMCAL::AliEMCAL()
  : AliDetector(),
    fBirkC0(0),
    fBirkC1(0.),
    fBirkC2(0.),
    fHighCharge(0.),
    fHighGain(0.),
    fHighLowGainFactor(0.),
    fLowGainOffset(0)
{
  // Default ctor 
  fName = "EMCAL" ;
  InitConstants();

}

//____________________________________________________________________________
AliEMCAL::AliEMCAL(const char* name, const char* title)
  : AliDetector(name,title),
    fBirkC0(0),
    fBirkC1(0.),
    fBirkC2(0.),
    fHighCharge(0.),
    fHighGain(0.),
    fHighLowGainFactor(0.),
    fLowGainOffset(0)
{
  //   ctor : title is used to identify the layout
  InitConstants();

}

//____________________________________________________________________________
AliEMCAL::~AliEMCAL()
{
  //dtor
}

//____________________________________________________________________________
void AliEMCAL::InitConstants()
{
  //initialize EMCAL values
  fBirkC0 = 1;
  fBirkC1 = 0.013/1.032;
  fBirkC2 = 9.6e-6/(1.032 * 1.032);
  
  fHighCharge        = 8.2 ;          // adjusted for a high gain range of 5.12 GeV (10 bits)
  fHighGain          = 6.64 ; 
  fHighLowGainFactor = 16. ;          // adjusted for a low gain range of 82 GeV (10 bits) 
  fLowGainOffset     = 1 ;            // offset added to the module id to distinguish high and low gain data
}

//____________________________________________________________________________
void AliEMCAL::Init()
{
  //EMCAL cuts
  Int_t * idtmed = fIdtmed->GetArray() - 1599 ; 
// --- Set decent energy thresholds for gamma and electron tracking

  // Tracking threshold for photons and electrons in Lead 
  Float_t cutgam=10.e-5; // 100 kev;
  Float_t cutele=10.e-5; // 100 kev;
  TString ntmp(GetTitle()); 
  ntmp.ToUpper();
  if(ntmp.Contains("10KEV")) {
    cutele = cutgam = 1.e-5;
  } else if(ntmp.Contains("50KEV")) {
    cutele = cutgam = 5.e-5;
  } else if(ntmp.Contains("100KEV")) {
    cutele = cutgam = 1.e-4;
  } else if(ntmp.Contains("200KEV")) {
    cutele = cutgam = 2.e-4;
  } else if(ntmp.Contains("500KEV")) {
    cutele = cutgam = 5.e-4;
  }

  gMC->Gstpar(idtmed[1600],"CUTGAM", cutgam);
  gMC->Gstpar(idtmed[1600],"CUTELE", cutele); // 1MEV -> 0.1MEV; 15-aug-05
  gMC->Gstpar(idtmed[1600],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1600],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  // --- Generate explicitly delta rays in Lead ---
  gMC->Gstpar(idtmed[1600], "LOSS", 3) ;
  gMC->Gstpar(idtmed[1600], "DRAY", 1) ;
  gMC->Gstpar(idtmed[1600], "DCUTE", cutele) ;
  gMC->Gstpar(idtmed[1600], "DCUTM", cutele) ;

// --- in aluminium parts ---
  gMC->Gstpar(idtmed[1602],"CUTGAM", cutgam) ;
  gMC->Gstpar(idtmed[1602],"CUTELE", cutele) ;
  gMC->Gstpar(idtmed[1602],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1602],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1602], "LOSS",3.) ;
  gMC->Gstpar(idtmed[1602], "DRAY",1.) ;
  gMC->Gstpar(idtmed[1602], "DCUTE", cutele) ;
  gMC->Gstpar(idtmed[1602], "DCUTM", cutele) ;

// --- and finally thresholds for photons and electrons in the scintillator ---
  gMC->Gstpar(idtmed[1601],"CUTGAM", cutgam) ;
  gMC->Gstpar(idtmed[1601],"CUTELE", cutele) ;// 1MEV -> 0.1MEV; 15-aug-05
  gMC->Gstpar(idtmed[1601],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1601],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1601], "LOSS",3) ; // generate delta rays 
  gMC->Gstpar(idtmed[1601], "DRAY",1) ;
  gMC->Gstpar(idtmed[1601], "DCUTE", cutele) ;
  gMC->Gstpar(idtmed[1601], "DCUTM", cutele) ;

  // S steel - 
  gMC->Gstpar(idtmed[1603],"CUTGAM", cutgam);
  gMC->Gstpar(idtmed[1603],"CUTELE", cutele);
  gMC->Gstpar(idtmed[1603],"BCUTE",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  gMC->Gstpar(idtmed[1603],"BCUTM",  cutgam);  // BCUTE and BCUTM start from GUTGUM
  // --- Generate explicitly delta rays 
  gMC->Gstpar(idtmed[1603], "LOSS",3);
  gMC->Gstpar(idtmed[1603], "DRAY",1);
  gMC->Gstpar(idtmed[1603], "DCUTE", cutele) ;
  gMC->Gstpar(idtmed[1603], "DCUTM", cutele) ;

  AliEMCALGeometry* geom = GetGeometry();
  if(geom->GetILOSS()>=0) {
    for(int i=1600; i<=1603; i++) gMC->Gstpar(idtmed[i], "LOSS", geom->GetILOSS()) ; 
  } 
  if(geom->GetIHADR()>=0) {
    for(int i=1600; i<=1603; i++) gMC->Gstpar(idtmed[i], "HADR", geom->GetIHADR()) ; 
  }
}

//____________________________________________________________________________
AliDigitizer* AliEMCAL::CreateDigitizer(AliRunDigitizer* manager) const
{
  //create and return the digitizer
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

  // 25-aug-04 by PAI - see  PMD/AliPMDv0.cxx for STEEL definition
  Float_t asteel[4] = { 55.847,51.9961,58.6934,28.0855 };
  Float_t zsteel[4] = { 26.,24.,28.,14. };
  Float_t wsteel[4] = { .715,.18,.1,.005 };
  AliMixture(4, "STAINLESS STEEL$", asteel, zsteel, 7.88, 4, wsteel);

  // DEFINITION OF THE TRACKING MEDIA

  // for EMCAL: idtmed[1599->1698] equivalent to fIdtmed[0->100]
  Int_t   isxfld = gAlice->Field()->Integ() ;
  Float_t sxmgmx = gAlice->Field()->Max() ;

  // Air                                                                         -> idtmed[1599]
 AliMedium(0, "Air$", 0, 0,
	     isxfld, sxmgmx, 10.0, 1.0, 0.1, 0.1, 10.0, 0, 0) ;

  // The Lead                                                                      -> idtmed[1600]
 
  AliMedium(1, "Lead$", 1, 0,
	     isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.1, 0.1, 0, 0) ;

 // The scintillator of the CPV made of Polystyrene scintillator                   -> idtmed[1601]
  float deemax = 0.1; // maximum fractional energy loss in one step (0 < DEEMAX â¤ 1);i
  AliMedium(2, "Scintillator$", 2, 1,
            isxfld, sxmgmx, 10.0, 0.001, deemax, 0.001, 0.001, 0, 0) ;

  // Various Aluminium parts made of Al                                            -> idtmed[1602]
  AliMedium(3, "Al$", 3, 0,
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;

  // 25-aug-04 by PAI : see  PMD/AliPMDv0.cxx for STEEL definition                 -> idtmed[1603]
  AliMedium(4, "S steel$", 4, 0, 
             isxfld, sxmgmx, 10.0, 0.1, 0.1, 0.001, 0.001, 0, 0) ;


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
  loader->LoadDigits("EMCAL");
  loader->GetEvent();
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
    AliError(Form("No geometry found !"));
    return;
  }
  
  AliAltroBuffer* buffer = NULL;
  Int_t prevDDL = -1;
  Int_t adcValuesLow[fgkTimeBins];
  Int_t adcValuesHigh[fgkTimeBins];

  //Load Mapping RCU files once
  TString path = gSystem->Getenv("ALICE_ROOT");
  path += "/EMCAL/mapping/RCU";
  TString path0 = path+"0.data";//This file will change in future
  TString path1 = path+"1.data";//This file will change in future
  AliAltroMapping * mapping[2] ; // For the moment only 2
  mapping[0] = new AliCaloAltroMapping(path0.Data());
  mapping[1] = new AliCaloAltroMapping(path1.Data());

  // loop over digits (assume ordered digits)
  for (Int_t iDigit = 0; iDigit < digits->GetEntries(); iDigit++) {
    AliEMCALDigit* digit = dynamic_cast<AliEMCALDigit *>(digits->At(iDigit)) ;
    if (digit->GetAmp() < fgThreshold) 
      continue;

    //get cell indeces
    Int_t nSM = 0;
    Int_t nIphi = 0;
    Int_t nIeta = 0;
    Int_t iphi = 0;
    Int_t ieta = 0;
    Int_t nModule = 0;
    geom->GetCellIndex(digit->GetId(), nSM, nModule, nIphi, nIeta);
    geom->GetCellPhiEtaIndexInSModule(nSM, nModule, nIphi, nIeta,iphi, ieta) ;
    
    //Check which is the RCU of the cell.
    Int_t iRCU = -111;
    //RCU0
    if (0<=iphi&&iphi<8) iRCU=0; // first cable row
    else if (8<=iphi&&iphi<16 && 0<=ieta&&ieta<24) iRCU=0; // first half; 
    //second cable row
    //RCU1
    else if(8<=iphi&&iphi<16 && 24<=ieta&&ieta<48) iRCU=1; // second half; 
    //second cable row
    else if(16<=iphi&&iphi<24) iRCU=1; // third cable row
    
    //Which DDL?
    Int_t iDDL = fgDDLPerSuperModule* nSM + iRCU;
    
    // new DDL
    if (iDDL != prevDDL) {
      // write real header and close previous file
       
      if (buffer) {
	buffer->Flush();
	buffer->WriteDataHeader(kFALSE, kFALSE);
	delete buffer;
      }
      
      // open new file and write dummy header
      TString fileName = AliDAQ::DdlFileName("EMCAL",iDDL);
      buffer = new AliAltroBuffer(fileName.Data(),mapping[iRCU]);
      buffer->WriteDataHeader(kTRUE, kFALSE);  //Dummy;
      prevDDL = iDDL;
    }
    
    // out of time range signal (?)
    if (digit->GetTimeR() > GetRawFormatTimeMax() ) {
      AliInfo("Signal is out of time range.\n");
      buffer->FillBuffer((Int_t)digit->GetAmp());
      buffer->FillBuffer(GetRawFormatTimeBins() );  // time bin
      buffer->FillBuffer(3);          // bunch length      
      buffer->WriteTrailer(3, ieta, iphi, nSM);  // trailer
      // calculate the time response function
    } else {

      Double_t energy = digit->GetAmp() * digitizer->GetECAchannel() + digitizer->GetECApedestal() ; 
      
      Bool_t lowgain = RawSampledResponse(digit->GetTimeR(), energy, adcValuesHigh, adcValuesLow) ; 
      
      if (lowgain) 
	buffer->WriteChannel(ieta, iphi, 0, GetRawFormatTimeBins(), adcValuesLow, fgThreshold);
      else 
	buffer->WriteChannel(ieta,iphi, 1, GetRawFormatTimeBins(), adcValuesHigh, fgThreshold);
           
    }
  }
  
  // write real header and close last file
  if (buffer) {
    buffer->Flush();
    buffer->WriteDataHeader(kFALSE, kFALSE);
    delete buffer;
  }
  mapping[0]->Delete();
  mapping[1]->Delete();
  loader->UnloadDigits();
}

//____________________________________________________________________________
void AliEMCAL::Raw2Digits(AliRawReader* reader)
{
  // convert raw data of the current event to digits
  AliEMCALGeometry * geom = GetGeometry();
  AliEMCALLoader * loader = dynamic_cast<AliEMCALLoader*>(fLoader) ; 

  // get the digits
  loader->CleanDigits(); // start from scratch
  loader->LoadDigits("EMCAL");
  TClonesArray* digits = loader->Digits() ;
  digits->Clear(); // yes, this is perhaps somewhat paranoid.. [clearing an extra time]

  if (!digits) {
    Error("Raw2Digits", "no digits found !");
    return;
  }
  if (!reader) {
    Error("Raw2Digits", "no raw reader found !");
    return;
  }

  // and get the digitizer too 
  loader->LoadDigitizer();
  AliEMCALDigitizer * digitizer = dynamic_cast<AliEMCALDigitizer *>(loader->Digitizer())  ; 

  // Use AliAltroRawStream to read the ALTRO format.  No need to
  // reinvent the wheel :-) 
  AliCaloRawStream in(reader,"EMCAL");
  // Select EMCAL DDL's;
  reader->Select("EMCAL");

  // reading is from previously existing AliEMCALGetter.cxx
  // ReadRaw method
  Bool_t first = kTRUE ;
 
  TF1 * signalF = new TF1("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);
  signalF->SetParNames("Charge", "Gain", "Amplitude", "TimeZero"); 
  
  Bool_t lowGainFlag = kFALSE ; 
  Int_t id =  -1;
  Int_t idigit = 0 ; 
  Int_t amp = 0 ; 
  Double_t time = 0. ; 
  Double_t energy = 0. ; 

  TGraph * gLowGain = new TGraph(GetRawFormatTimeBins()) ; 
  TGraph * gHighGain= new TGraph(GetRawFormatTimeBins()) ;  

  while ( in.Next() ) { // EMCAL entries loop 
    if ( in.IsNewRow() ) {//phi
      if ( in.IsNewColumn() ) {//eta
	id =  geom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
	if (!first) {
	  FitRaw(lowGainFlag, gLowGain, gHighGain, signalF, energy, time) ; 
	  
	  if (time == 0. && energy == 0.) { 
	    amp = 0 ; 
	  }
	  else {
	    amp = static_cast<Int_t>( (energy - digitizer->GetECApedestal()) / digitizer->GetECAchannel() + 0.5 ) ; 
	  }
	  
	  if (amp > 0) {
	    new((*digits)[idigit]) AliEMCALDigit( -1, -1, id, amp, time) ;	
	    idigit++ ; 
	  }
	  Int_t index ; 
	  for (index = 0; index < GetRawFormatTimeBins(); index++) {
	    gLowGain->SetPoint(index, index * GetRawFormatTimeMax() / GetRawFormatTimeBins(), 0) ;  
	    gHighGain->SetPoint(index, index * GetRawFormatTimeMax() / GetRawFormatTimeBins(), 0) ; 
	  } 
	} // not first  
	first = kFALSE ; 
	id =  geom->GetAbsCellIdFromCellIndexes(in.GetModule(), in.GetRow(), in.GetColumn()) ;
	if (in.GetModule() == GetRawFormatLowGainOffset() ) {
	  lowGainFlag = kTRUE ; 
	}
	else { 
	  lowGainFlag = kFALSE ; 
	}
      } //new column-eta
    }// new row-phi
    if (lowGainFlag) {
      gLowGain->SetPoint(in.GetTime(), 
			 in.GetTime()* GetRawFormatTimeMax() / GetRawFormatTimeBins(), 
			 in.GetSignal()) ;
    }
    else { 
      gHighGain->SetPoint(in.GetTime(), 
			  in.GetTime() * GetRawFormatTimeMax() / GetRawFormatTimeBins(), 
			  in.GetSignal() ) ;
    }
  } // EMCAL entries loop
  digits->Sort() ; 
  
  delete signalF ; 
  delete gLowGain;
  delete gHighGain ; 
  
  return ; 
}

//____________________________________________________________________________ 
void AliEMCAL::FitRaw(Bool_t lowGainFlag, TGraph * gLowGain, TGraph * gHighGain, TF1* signalF, Double_t & energy, Double_t & time)
{
  // Fits the raw signal time distribution; from AliEMCALGetter 

  const Int_t kNoiseThreshold = 0 ;
  Double_t timezero1 = 0., timezero2 = 0., timemax = 0. ;
  Double_t signal = 0., signalmax = 0. ;       
  energy = time = 0. ; 

  if (lowGainFlag) {
    timezero1 = timezero2 = signalmax = timemax = 0. ;
    signalF->FixParameter(0, GetRawFormatLowCharge()) ; 
    signalF->FixParameter(1, GetRawFormatLowGain()) ; 
    Int_t index ; 
    for (index = 0; index < GetRawFormatTimeBins(); index++) {
      gLowGain->GetPoint(index, time, signal) ; 
      if (signal > kNoiseThreshold && timezero1 == 0.) 
	timezero1 = time ;
      if (signal <= kNoiseThreshold && timezero1 > 0. && timezero2 == 0.)
	timezero2 = time ; 
      if (signal > signalmax) {
	signalmax = signal ; 
	timemax   = time ; 
      }
    }
    signalmax /= RawResponseFunctionMax(GetRawFormatLowCharge(), 
						GetRawFormatLowGain()) ;
    if ( timezero1 + GetRawFormatTimePeak() < GetRawFormatTimeMax() * 0.4 ) { // else its noise 
      signalF->SetParameter(2, signalmax) ; 
      signalF->SetParameter(3, timezero1) ;    	    
      gLowGain->Fit(signalF, "QRON", "", 0., timezero2); //, "QRON") ; 
      energy = signalF->GetParameter(2) ; 
      time   = signalF->GetMaximumX() - GetRawFormatTimePeak() - GetRawFormatTimeTrigger() ;
    }
  } else {
    timezero1 = timezero2 = signalmax = timemax = 0. ;
    signalF->FixParameter(0, GetRawFormatHighCharge()) ; 
    signalF->FixParameter(1, GetRawFormatHighGain()) ; 
    Int_t index ; 
    for (index = 0; index < GetRawFormatTimeBins(); index++) {
      gHighGain->GetPoint(index, time, signal) ;               
      if (signal > kNoiseThreshold && timezero1 == 0.) 
	timezero1 = time ;
      if (signal <= kNoiseThreshold && timezero1 > 0. && timezero2 == 0.)
	timezero2 = time ; 
      if (signal > signalmax) {
	signalmax = signal ;   
	timemax   = time ; 
      }
    }
    signalmax /= RawResponseFunctionMax(GetRawFormatHighCharge(), 
						GetRawFormatHighGain()) ;;
    if ( timezero1 + GetRawFormatTimePeak() < GetRawFormatTimeMax() * 0.4 ) { // else its noise  
      signalF->SetParameter(2, signalmax) ; 
      signalF->SetParameter(3, timezero1) ;               
      gHighGain->Fit(signalF, "QRON", "", 0., timezero2) ; 
      energy = signalF->GetParameter(2) ; 
      time   = signalF->GetMaximumX() - GetRawFormatTimePeak() - GetRawFormatTimeTrigger() ;
    }
  }
  
  return;
}

//____________________________________________________________________________
void AliEMCAL::Hits2SDigits()  
{ 
// create summable digits

  GetGeometry();
  AliEMCALSDigitizer emcalDigitizer(fLoader->GetRunLoader()->GetFileName().Data()) ;
  emcalDigitizer.SetEventRange(0, -1) ; // do all the events
  emcalDigitizer.ExecuteTask() ;
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
  //compute the maximum of the raw response function and return
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
