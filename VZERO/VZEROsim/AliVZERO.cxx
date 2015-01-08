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

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//                          V-Zero   Detector                            //
//  This class contains the base procedures for the VZERO  detector      //
//  Default geometry of November 2003 :   V0R box is 4.4 cm thick        //
//                                  scintillators are 2 cm thick         //
//  All comments should be sent to Brigitte CHEYNIS :                    //
//                                 b.cheynis@ipnl.in2p3.fr               //
//                                                                       //
//                                                                       //
///////////////////////////////////////////////////////////////////////////


// --- Standard libraries ---
#include <Riostream.h>
#include <stdlib.h>

// --- ROOT libraries ---
#include <TNamed.h>
#include "TROOT.h"
#include "TFile.h"
#include "TNetFile.h"
#include "TRandom.h"
#include "TTree.h"
#include "TBranch.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
#include "TParameter.h"
#include "TF1.h"

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliMC.h"
#include "AliVZERO.h"
#include "AliVZEROLoader.h"
#include "AliVZERODigitizer.h"
#include "AliVZEROBuffer.h"
#include "AliDigitizationInput.h"
#include "AliVZEROdigit.h"
#include "AliVZEROSDigit.h"
#include "AliDAQ.h"
#include "AliRawReader.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliVZERORawStream.h"
#include "AliVZEROCalibData.h"
#include "AliVZERORecoParam.h"
#include "AliVZEROReconstructor.h"

ClassImp(AliVZERO)
 //__________________________________________________________________
AliVZERO::AliVZERO(): AliDetector(),
          fIdSens1(0),
          fThickness(0.),
	  fThickness1(0.),
	  fMaxStepQua(0.),
	  fMaxStepAlu(0.),
	  fMaxDestepQua(0.),
          fMaxDestepAlu(0.),
          fCalibData(NULL),
          fTimeSlewing(NULL),
          fSignalShape(NULL),
          fRecoParam(NULL)
{
/// Default Constructor
    
    AliDebug(1,Form("default (empty) ctor this=%p",this));
    fIshunt          = 0;

    for(Int_t i = 0 ; i < 64; ++i) {
      fNBins[i] = 0;
      fBinSize[i] = 0;
    }
}

//_____________________________________________________________________________
AliVZERO::AliVZERO(const char *name, const char *title)
       : AliDetector(name,title),
         fIdSens1(0),
         fThickness(4.4),
	 fThickness1(2.0),
	 fMaxStepQua(0.05),
	 fMaxStepAlu(0.01),
	 fMaxDestepQua(-1.0),
	 fMaxDestepAlu(-1.0),
	 fCalibData(NULL),
	 fTimeSlewing(NULL),
         fSignalShape(NULL),
	 fRecoParam(NULL)
{
  
  // Standard constructor for VZERO Detector
  
  AliDebug(1,Form("ctor this=%p",this));
  
  //  fIshunt       =  1;  // All hits are associated with primary particles  
   
  fHits         =  new TClonesArray("AliVZEROhit", 400);
  fDigits       =  new TClonesArray("AliVZEROdigit",400); 
   
  gAlice->GetMCApp()->AddHitList(fHits);

//   fThickness    =  4.4;   // total thickness of the V0R box in cm
//   fThickness1   =  2.0;   // thickness of scintillating cells in cm
//   
//   fMaxStepQua   =  0.05; 
//   fMaxStepAlu   =  0.01; 
//   
//   fMaxDestepQua =  -1.0;
//   fMaxDestepAlu =  -1.0;
  
  for(Int_t i = 0 ; i < 64; ++i) {
    fNBins[i] = 0;
    fBinSize[i] = 0;
  }
}

//_____________________________________________________________________________
AliVZERO::~AliVZERO()
{
  //
  // Default destructor for VZERO Detector
  //
  
    if (fHits) {
        fHits->Delete();
        delete fHits;
	fHits=0; }
    
    if (fDigits) {
        fDigits->Delete();
        delete fDigits;
        fDigits=0; }
    if (fSignalShape) {
      delete fSignalShape;
      fSignalShape = NULL;
    }
    if (fRecoParam) {
      delete fRecoParam;
      fRecoParam = NULL;
    }
}

//_____________________________________________________________________________
void AliVZERO::CreateGeometry()
{
  //
  // Builds simple Geant3 geometry 
  //
}
//_____________________________________________________________________________
void AliVZERO::CreateMaterials()
{
  //
  // Creates materials used for Geant3 geometry 
  //
}

//_____________________________________________________________________________
void AliVZERO::Init()
{
  //
  // Initialises the VZERO  class after it has been built
  //
}


//_____________________________________________________________________________
void AliVZERO::SetMaxStepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials
  //
     fMaxStepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxStepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum foils (not used in 
  // version v2)
  //
    fMaxStepAlu = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepQua(Float_t p1)
{
  //
  // Possible parametrisation of steps in active materials (quartz)
  //
    fMaxDestepQua = p1;
}

//_____________________________________________________________________________
void AliVZERO::SetMaxDestepAlu(Float_t p1)
{
  //
  // Possible parametrisation of steps in Aluminum (not used in 
  // version v2)
  //
    fMaxDestepAlu = p1;
}

//_____________________________________________________________________________
AliLoader* AliVZERO::MakeLoader(const char* topfoldername)
{ 
  //
  // Builds VZEROgetter (AliLoader type)
  // if detector wants to use customized getter, it must overload this method
  //
//  Info("MakeLoader","Creating AliVZEROLoader. Top folder is %s.",topfoldername); 
 
  AliDebug(1,Form("Creating AliVZEROLoader, Top folder is %s ",topfoldername));
  fLoader = new AliVZEROLoader(GetName(),topfoldername);
  return fLoader;
}

//_____________________________________________________________________________
void AliVZERO::SetTreeAddress()
{
  //
  // Sets tree address for hits.
  //
  if (fLoader->TreeH() && (fHits == 0x0))
    fHits = new  TClonesArray("AliVZEROhit", 400);

  AliDetector::SetTreeAddress();
}

//_____________________________________________________________________________
AliDigitizer* AliVZERO::CreateDigitizer(AliDigitizationInput* digInput) const
{
  //
  // Creates a digitizer for VZERO
  //
  return new AliVZERODigitizer(digInput);
}

//_____________________________________________________________________________
void AliVZERO::Hits2Digits(){
  //
  // Converts hits to digits
  //
  // Creates the VZERO digitizer 
  AliVZERODigitizer* dig = new AliVZERODigitizer(this,AliVZERODigitizer::kHits2Digits);

  // Creates the digits
  dig->Digitize("");

  // deletes the digitizer
  delete dig;
}

//_____________________________________________________________________________
void AliVZERO::Hits2SDigits(){
  //
  // Converts hits to summable digits
  //
  // Creates the VZERO digitizer 
  AliVZERODigitizer* dig = new AliVZERODigitizer(this,AliVZERODigitizer::kHits2SDigits);

  // Creates the sdigits
  dig->Digitize("");

  // deletes the digitizer
  delete dig;
}

//_____________________________________________________________________________
void AliVZERO::Digits2Raw()
{
   //
   //  Converts digits of the current event to raw data
   //
  
   AliVZERO *fVZERO = (AliVZERO*)gAlice->GetDetector("VZERO");
   fLoader->LoadDigits();
   TTree* digits = fLoader->TreeD();
   if (!digits) {
      Error("Digits2Raw", "no digits tree");
      return;
   }
   TClonesArray * VZEROdigits = new TClonesArray("AliVZEROdigit",1000);
   fVZERO->SetTreeAddress();  		
   digits->GetBranch("VZERODigit")->SetAddress(&VZEROdigits); 
  
   const char *fileName    = AliDAQ::DdlFileName("VZERO",0);
   AliVZEROBuffer* buffer  = new AliVZEROBuffer(fileName);
  
   // Get Trigger information first
   // Read trigger inputs from trigger-detector object
   AliDataLoader * dataLoader = fLoader->GetDigitsDataLoader();
   if( !dataLoader->IsFileOpen() ) 
        dataLoader->OpenFile( "READ" );
   AliTriggerDetector* trgdet = (AliTriggerDetector*)dataLoader->GetDirectory()->Get( "Trigger" );
   UInt_t triggerInfo = 0;
   if(trgdet) {
      triggerInfo = trgdet->GetMask() & 0xffff;
   }
   else {
      AliError(Form("There is no trigger object for %s",fLoader->GetName()));
   }

   buffer->WriteTriggerInfo((UInt_t)triggerInfo); 
   buffer->WriteTriggerScalers(); 
   buffer->WriteBunchNumbers(); 
  
   Int_t aBBflagsV0A = 0;
   Int_t aBBflagsV0C = 0;
   Int_t aBGflagsV0A = 0;
   Int_t aBGflagsV0C = 0;

   if (digits->GetUserInfo()->FindObject("BBflagsV0A")) {
     aBBflagsV0A = ((TParameter<int>*)digits->GetUserInfo()->FindObject("BBflagsV0A"))->GetVal();
   }
   else
     AliWarning("V0A beam-beam flags not found in digits tree UserInfo! The flags will not be written to the raw-data stream!");

   if (digits->GetUserInfo()->FindObject("BBflagsV0C")) {
     aBBflagsV0C = ((TParameter<int>*)digits->GetUserInfo()->FindObject("BBflagsV0C"))->GetVal();
   }
   else
     AliWarning("V0C beam-beam flags not found in digits tree UserInfo! The flags will not be written to the raw-data stream!");

   if (digits->GetUserInfo()->FindObject("BGflagsV0A")) {
     aBGflagsV0A = ((TParameter<int>*)digits->GetUserInfo()->FindObject("BGflagsV0A"))->GetVal();
   }
   else
     AliWarning("V0A beam-gas flags not found in digits tree UserInfo! The flags will not be written to the raw-data stream!");

   if (digits->GetUserInfo()->FindObject("BGflagsV0C")) {
     aBGflagsV0C = ((TParameter<int>*)digits->GetUserInfo()->FindObject("BGflagsV0C"))->GetVal();
   }
   else
     AliWarning("V0C beam-gas flags not found in digits tree UserInfo! The flags will not be written to the raw-data stream!");

   // Now retrieve the channel information: charge smaples + time and 
   // dump it into ADC and Time arrays
   Int_t nEntries = Int_t(digits->GetEntries());
   Short_t aADC[64][AliVZEROdigit::kNClocks];
   Float_t aTime[64];
   Float_t aWidth[64];
   Bool_t  aIntegrator[64];
   Bool_t  aBBflag[64];
   Bool_t  aBGflag[64];
  
   for (Int_t i = 0; i < nEntries; i++) {
     fVZERO->ResetDigits();
     digits->GetEvent(i);
     Int_t ndig = VZEROdigits->GetEntriesFast(); 
   
     if(ndig == 0) continue;
     for(Int_t k=0; k<ndig; k++){
         AliVZEROdigit* fVZERODigit = (AliVZEROdigit*) VZEROdigits->At(k);
         // Convert aliroot channel k into FEE channel iChannel before writing data
	 Int_t iChannel       = AliVZEROCalibData::GetBoardNumber(fVZERODigit->PMNumber()) * 8 +
	   AliVZEROCalibData::GetFEEChannelNumber(fVZERODigit->PMNumber());
	 for(Int_t iClock = 0; iClock < AliVZEROdigit::kNClocks; ++iClock) aADC[iChannel][iClock] = fVZERODigit->ChargeADC(iClock);
	 aTime[iChannel]      = fVZERODigit->Time();
	 aWidth[iChannel]     = fVZERODigit->Width();
	 aIntegrator[iChannel]= fVZERODigit->Integrator();
	 if(fVZERODigit->PMNumber() < 32) {
	   aBBflag[iChannel] = (aBBflagsV0C >> fVZERODigit->PMNumber()) & 0x1;
	   aBGflag[iChannel] = (aBGflagsV0C >> fVZERODigit->PMNumber()) & 0x1;
	 }
	 else {
	   aBBflag[iChannel] = (aBBflagsV0A >> (fVZERODigit->PMNumber()-32)) & 0x1;
	   aBGflag[iChannel] = (aBGflagsV0A >> (fVZERODigit->PMNumber()-32)) & 0x1;
	 }
         AliDebug(1,Form("DDL: %s\tdigit number: %d\tPM number: %d\tADC: %d\tTime: %f",
			 fileName,k,fVZERODigit->PMNumber(),aADC[iChannel][AliVZEROdigit::kNClocks/2],aTime[iChannel])); 
     }        
   }

   // Now fill raw data
          
   for (Int_t  iCIU = 0; iCIU < 8; iCIU++) { 
 
   // decoding of one Channel Interface Unit numbered iCIU - there are 8 channels per CIU (and 8 CIUs) :
  
      for(Int_t iChannel_Offset = iCIU*8; iChannel_Offset < (iCIU*8)+8; iChannel_Offset=iChannel_Offset+4) { 
         for(Int_t iChannel = iChannel_Offset; iChannel < iChannel_Offset+4; iChannel++) {
             buffer->WriteChannel(iChannel, aADC[iChannel], aIntegrator[iChannel]);       
         }
         buffer->WriteBeamFlags(&aBBflag[iChannel_Offset],&aBGflag[iChannel_Offset]); 
         buffer->WriteMBInfo(); 
         buffer->WriteMBFlags();   
         buffer->WriteBeamScalers(); 
      } 

      for(Int_t iChannel = iCIU*8 + 7; iChannel >= iCIU*8; iChannel--) {
          buffer->WriteTiming(aTime[iChannel], aWidth[iChannel]); 
      }

    // End of decoding of one CIU card
    
  } // end of decoding the eight CIUs
     
  delete buffer;
  fLoader->UnloadDigits();  
}

//_____________________________________________________________________________
Bool_t AliVZERO::Raw2SDigits(AliRawReader* rawReader){
  // Converts the VZERO raw data into digits
  // The method is used for merging simulated and
  // real data events
  TStopwatch timer;
  timer.Start();

  if(!fLoader) {
    AliError("no VZERO loader found");
    return kFALSE;
  }
  fLoader->LoadSDigits("UPDATE");

  if (!fLoader->TreeS()) fLoader->MakeTree("S");
  fLoader->MakeSDigitsContainer();
  TTree* treeS  = fLoader->TreeS();

  TClonesArray *sdigits = new TClonesArray("AliVZEROSDigit", 64);
  treeS->Branch("VZEROSDigit", &sdigits); 

  {
    rawReader->Reset();
    AliVZERORawStream rawStream(rawReader);    
     
    if (!rawStream.Next()) return kFALSE; // No VZERO data found

    GetCalibData();

    Int_t nSDigits = 0;
    Float_t *charges = NULL;
    Int_t nbins = 0;
    for(Int_t iChannel=0; iChannel < 64; ++iChannel) {
      Int_t offlineCh = rawStream.GetOfflineChannel(iChannel);
      Short_t chargeADC[AliVZEROdigit::kNClocks];
      for(Int_t iClock=0; iClock < AliVZEROdigit::kNClocks; ++iClock) {
	chargeADC[iClock] = rawStream.GetPedestal(iChannel,iClock);
      }
      // Integrator flag
      Bool_t integrator = rawStream.GetIntegratorFlag(iChannel,AliVZEROdigit::kNClocks/2);
      // HPTDC data (leading time and width)
      Int_t board = AliVZEROCalibData::GetBoardNumber(offlineCh);
      Float_t time = rawStream.GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      //      Float_t width = rawStream.GetWidth(iChannel)*fCalibData->GetWidthResolution(board);
      Float_t adc = 0;

      // Pedestal retrieval and suppression
      Float_t maxadc = 0;
      Int_t imax = -1;
      Float_t adcPedSub[AliVZEROdigit::kNClocks];
      Float_t integral = fSignalShape->Integral(0,200);
      for(Int_t iClock=0; iClock < AliVZEROdigit::kNClocks; ++iClock) {
	Bool_t iIntegrator = (iClock%2 == 0) ? integrator : !integrator;
	Int_t k = offlineCh + 64*iIntegrator;
	adcPedSub[iClock] = (Float_t)chargeADC[iClock] - fCalibData->GetPedestal(k);
	if(adcPedSub[iClock] <= fRecoParam->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	  adcPedSub[iClock] = 0;
	  continue;
	}
	if(iClock < fRecoParam->GetStartClock() || iClock > fRecoParam->GetEndClock()) continue;
	if(adcPedSub[iClock] > maxadc) {
	  maxadc = adcPedSub[iClock];
	  imax   = iClock;
	}
      }
      if (imax != -1) {
	Int_t start = imax - fRecoParam->GetNPreClocks();
	if (start < 0) start = 0;
	Int_t end = imax + fRecoParam->GetNPostClocks();
	if (end > 20) end = 20;
	for(Int_t iClock = start; iClock <= end; iClock++) {
	  adc += adcPedSub[iClock];
	}
      }
      Float_t correctedTime = CorrectLeadingTime(offlineCh,time,adc);

      if (!charges) {
	nbins = fNBins[offlineCh];
	charges = new Float_t[nbins];
      }
      else if (nbins != fNBins[offlineCh]) {
	delete [] charges;
	nbins = fNBins[offlineCh];
	charges = new Float_t[nbins];
      }
      memset(charges,0,nbins*sizeof(Float_t));

      // Now lets produce SDigit
      if ((correctedTime > (AliVZEROReconstructor::kInvalidTime + 1e-6)) &&
	  (adc > 1e-6)) {
	for(Int_t iBin = 0; iBin < nbins; ++iBin) {
	  Float_t t = fBinSize[offlineCh]*Float_t(iBin);
	  if ((t < correctedTime) ||
	      (t > (correctedTime+200.))) continue;
	  charges[iBin] = kChargePerADC*adc*(fSignalShape->Eval(t-correctedTime)*fBinSize[offlineCh]/integral);
	}
      }

      TClonesArray &sdigitsref = *sdigits;  
      new (sdigitsref[nSDigits++]) AliVZEROSDigit(offlineCh,fNBins[offlineCh],charges);
    }
    if (charges) delete [] charges;
  }
  
  treeS->Fill();
  fLoader->WriteSDigits("OVERWRITE");
  fLoader->UnloadSDigits();	
	
  timer.Stop();
  timer.Print();
  return kTRUE;
}

//_____________________________________________________________________________
void AliVZERO::GetCalibData()
{
  // Gets calibration object for VZERO set
  // Do nothing in case it is already loaded
  if (fCalibData) return;

  AliCDBEntry *entry = AliCDBManager::Instance()->Get("VZERO/Calib/Data");
  if (entry) fCalibData = (AliVZEROCalibData*) entry->GetObject();
  if (!fCalibData)  AliFatal("No calibration data from calibration database !");

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("VZERO/Calib/TimeSlewing");
  if (!entry2) AliFatal("VZERO time slewing function is not found in OCDB !");
  fTimeSlewing = (TF1*)entry2->GetObject();

  for(Int_t i = 0 ; i < 64; ++i) {
    Int_t board = AliVZEROCalibData::GetBoardNumber(i);
    fNBins[i] = TMath::Nint(((Float_t)(fCalibData->GetMatchWindow(board)+1)*25.0+
			     (Float_t)kMaxTDCWidth*fCalibData->GetWidthResolution(board))/
			    fCalibData->GetTimeResolution(board));
    fBinSize[i] = fCalibData->GetTimeResolution(board);
  }

  fSignalShape = new TF1("VZEROSDigitSignalShape",this,&AliVZERO::SignalShape,0,200,6,"AliVZERO","SignalShape");
  fSignalShape->SetParameters(0,1.57345e1,-4.25603e-1,
			      2.9,6.40982,3.69339e-01);

  fRecoParam = new AliVZERORecoParam;

  return;
}

Float_t AliVZERO::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6) return -1024;

  // In case of pathological signals
  if (adc < 1e-6) return time;

  // Slewing correction
  Float_t thr = fCalibData->GetCalibDiscriThr(i,kTRUE);
  time -= fTimeSlewing->Eval(adc/thr);

  return time;
}

double AliVZERO::SignalShape(double *x, double *par)
{
  // this function simulates the signal
  // shape used in Raw->SDigits method

  Double_t xx = x[0];
  if (xx <= par[0]) return 0;
  Double_t a = 1./TMath::Power((xx-par[0])/par[1],1./par[2]);
  if (xx <= par[3]) return a;
  Double_t b = 1./TMath::Power((xx-par[3])/par[4],1./par[5]);
  Double_t f = a*b/(a+b);
  AliDebug(100,Form("x=%f func=%f",xx,f));
  return f;
}
