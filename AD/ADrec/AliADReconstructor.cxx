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

/* $Id: AliADReconstructor.cxx 20956 2007-09-26 14:22:18Z mrodrigu $ */
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  Class for AD reconstruction                                         //
//////////////////////////////////////////////////////////////////////////////
#include <TParameter.h>
#include <TGeoManager.h>
#include <TGeoMatrix.h>
#include <TGeoPhysicalNode.h>
#include <AliGeomManager.h>

#include "AliRawReader.h"
#include "AliGRPObject.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliESDEvent.h"
#include "AliRunInfo.h"
#include "AliCTPTimeParams.h"
#include "AliLHCClockPhase.h"
#include "TSpline.h"

#include "AliADReconstructor.h"
#include "AliADdigit.h"
#include "AliESDAD.h"
#include "AliESDADfriend.h"
#include "AliADConst.h"
#include "AliADCalibData.h"
#include "AliADRawStream.h"
#include "AliADDecision.h"

ClassImp(AliADReconstructor)
//_____________________________________________________________________________
AliADReconstructor:: AliADReconstructor():
  AliReconstructor(),
  fESDAD(0x0),
  fESD(0x0),
  fESDADfriend(0x0),
  fCalibData(NULL),
  fDigitsArray(0),
  fCollisionMode(0),
  fBeamEnergy(0.)

{
  
  // Default constructor  
  // Get calibration data
  fCalibData = GetCalibData();
  //Get time slewing
  GetTimeSlewingSplines();
  
  // Now get the CTP L0->L1 delay
  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/CTP/CTPtiming");
  if (!entry) AliFatal("CTP timing parameters are not found in OCDB !");
  AliCTPTimeParams *ctpParams = (AliCTPTimeParams*)entry->GetObject();
  Float_t l1Delay = (Float_t)ctpParams->GetDelayL1L0()*25.0;

  AliCDBEntry *entry1 = AliCDBManager::Instance()->Get("GRP/CTP/TimeAlign");
  if (!entry1) AliFatal("CTP time-alignment is not found in OCDB !");
  AliCTPTimeParams *ctpTimeAlign = (AliCTPTimeParams*)entry1->GetObject();
  l1Delay += ((Float_t)ctpTimeAlign->GetDelayL1L0()*25.0);

  AliCDBEntry *entry2 = AliCDBManager::Instance()->Get("AD/Calib/TimeDelays");
  if (!entry2) AliFatal("AD time delays are not found in OCDB !");
  TH1F *TimeDelays = (TH1F*)entry2->GetObject();

  AliCDBEntry *entry3 = AliCDBManager::Instance()->Get("GRP/Calib/LHCClockPhase");
  if (!entry3) AliFatal("LHC clock-phase shift is not found in OCDB !");
  AliLHCClockPhase *phase = (AliLHCClockPhase*)entry3->GetObject();

  for(Int_t i = 0 ; i < 16; ++i) {
    Int_t board = AliADCalibData::GetBoardNumber(i); 
    fHptdcOffset[i] = (((Float_t)fCalibData->GetRollOver(board)-
			(Float_t)fCalibData->GetTriggerCountOffset(board))*25.0
		       +fCalibData->GetTimeOffset(i)
		       -l1Delay
		       -phase->GetMeanPhase()
		       -TimeDelays->GetBinContent(i+1)
		       -kADOffset);
   }
   
  Float_t zADC1 = TMath::Abs(GetZPosition("AD/ADC1"));//outer 
  Float_t zADC2 = TMath::Abs(GetZPosition("AD/ADC2"));//inner
  Float_t zADA1 = TMath::Abs(GetZPosition("AD/ADA1"));//inner
  Float_t zADA2 = TMath::Abs(GetZPosition("AD/ADA2"));//outer
  // distance in time units from nominal vertex to AD
  fLayerDist[0] = zADC1/TMath::Ccgs()*1e9;
  fLayerDist[1] = zADC2/TMath::Ccgs()*1e9;
  fLayerDist[2] = zADA1/TMath::Ccgs()*1e9;
  fLayerDist[3] = zADA2/TMath::Ccgs()*1e9; 

}
//________________________________________________________________________________
Double_t AliADReconstructor::GetZPosition(const char* symname)
{
// Get the global z coordinate of the given AD alignable volume
//
  Double_t *tr;
  TGeoPNEntry *pne = gGeoManager->GetAlignableEntry(symname);
  if (!pne) {
    AliFatalClass(Form("TGeoPNEntry with symbolic name %s does not exist!",symname));
    return 0;
  }

  TGeoPhysicalNode *pnode = pne->GetPhysicalNode();
  if(pnode){
          TGeoHMatrix* hm = pnode->GetMatrix();
           tr = hm->GetTranslation();
  }else{
          const char* path = pne->GetTitle();
          if(!gGeoManager->cd(path)){
                  AliFatalClass(Form("Volume path %s not valid!",path));
                  return 0;
          }
         tr = gGeoManager->GetCurrentMatrix()->GetTranslation();
  }
  return tr[2];

}


//_____________________________________________________________________________
AliADReconstructor& AliADReconstructor::operator = 
  (const AliADReconstructor& /*reconstructor*/)
{
// assignment operator

  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliADReconstructor::~AliADReconstructor()
{
// destructor
  delete fESDAD;
  delete fESDADfriend;
  delete fDigitsArray;
}

//_____________________________________________________________________________
void AliADReconstructor::Init()
{
// initializer
    fESDAD  = new AliESDAD;
    fESDADfriend  = new AliESDADfriend;
}

//_____________________________________________________________________________
void AliADReconstructor::ConvertDigits(AliRawReader* rawReader, TTree* digitsTree) const
{
// converts RAW to digits 

  if (!digitsTree) {
    AliError("No digits tree!");
    return;
  }

  if (!fDigitsArray) fDigitsArray = new TClonesArray("AliADdigit", 16);
  digitsTree->Branch("ADDigit", &fDigitsArray);

  rawReader->Reset();
  AliADRawStream rawStream(rawReader);
  if (rawStream.Next()) {
   
    Bool_t aBBflag[16];
    Bool_t aBGflag[16];
    
    for(Int_t iChannel=0; iChannel < 16; ++iChannel) {
      Int_t offlineCh = kOfflineChannel[iChannel];
      // ADC charge samples
      Short_t chargeADC[kADNClocks];
      aBBflag[offlineCh] = kFALSE;
      aBGflag[offlineCh] = kFALSE;
      for(Int_t iClock=0; iClock < kADNClocks; ++iClock) {
	chargeADC[iClock] = rawStream.GetPedestal(iChannel,iClock);
	if(rawStream.GetBBFlag(iChannel,iClock)) aBBflag[offlineCh]=kTRUE;
	if(rawStream.GetBGFlag(iChannel,iClock)) aBGflag[offlineCh]=kTRUE;
      }
      // Integrator flag
      Bool_t integrator = rawStream.GetIntegratorFlag(iChannel,kADNClocks/2);
      Bool_t BBflag = rawStream.GetBBFlag(iChannel,kADNClocks/2); 
      Bool_t BGflag = rawStream.GetBGFlag(iChannel,kADNClocks/2);
   
      // HPTDC data (leading time and width)
      Int_t board = AliADCalibData::GetBoardNumber(offlineCh);
      Float_t time = rawStream.GetTime(iChannel)*fCalibData->GetTimeResolution(board);
      Float_t width = rawStream.GetWidth(iChannel)*fCalibData->GetWidthResolution(board);
      // Add a digit
      //if(!fCalibData->IsChannelDead(offlineCh)){ Off for the moment
	  new ((*fDigitsArray)[fDigitsArray->GetEntriesFast()]) AliADdigit(offlineCh, time, width,integrator, chargeADC, BBflag, BGflag);
      //}
      
      fESDADfriend->SetBBScalers(offlineCh,rawStream.GetBBScalers(iChannel));
      fESDADfriend->SetBGScalers(offlineCh,rawStream.GetBGScalers(iChannel));
      for (Int_t iEv = 0; iEv < kADNClocks; iEv++) {
	  fESDADfriend->SetBBFlag(offlineCh,iEv,rawStream.GetBBFlag(iChannel,iEv));
	  fESDADfriend->SetBGFlag(offlineCh,iEv,rawStream.GetBGFlag(iChannel,iEv));
      }
    }
    //BC Unmasked triggers
    Int_t pBBmulADA = 0;
    Int_t pBBmulADC = 0;
    Int_t pBGmulADA = 0;
    Int_t pBGmulADC = 0;
  
    for(Int_t iChannel=0; iChannel<4; iChannel++) {//Loop over pairs of pads
    //Enable time is used to turn off the coincidence 
      if((!fCalibData->GetEnableTiming(iChannel) || aBBflag[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || aBBflag[iChannel+4])) pBBmulADC++;
      if((!fCalibData->GetEnableTiming(iChannel) || aBGflag[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || aBGflag[iChannel+4])) pBGmulADC++;

      if((!fCalibData->GetEnableTiming(iChannel+8) || aBBflag[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBBflag[iChannel+12])) pBBmulADA++;
      if((!fCalibData->GetEnableTiming(iChannel+8) || aBGflag[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBGflag[iChannel+12])) pBGmulADA++;
	}
	
    Bool_t UBA = kFALSE;
    Bool_t UBC = kFALSE;
    Bool_t UGA = kFALSE;
    Bool_t UGC = kFALSE;

    if(pBBmulADA>=fCalibData->GetBBAThreshold()) UBA = kTRUE;
    if(pBBmulADC>=fCalibData->GetBBCThreshold()) UBC = kTRUE;
    if(pBGmulADA>=fCalibData->GetBGAThreshold()) UGA = kTRUE;
    if(pBGmulADC>=fCalibData->GetBGCThreshold()) UGC = kTRUE;
  
    UShort_t fTrigger = 0;
    if(UBA && UBC) fTrigger |= 1;
    if(UBA || UBC) fTrigger |= (1<<1);
    if(UGA && UBC) fTrigger |= (1<<2);
    if(UGA) fTrigger |= (1<<3);
    if(UGC && UBA) fTrigger |= (1<<4);
    if(UGC) fTrigger |= (1<<5);
    if(UBA) fTrigger |= (1<<12);
    if(UBC) fTrigger |= (1<<13); 
    if(UGA || UGC) fTrigger |= (1<<14);
    if((UGA && UBC) || (UGC && UBA)) fTrigger |= (1<<15);
   
    fESDADfriend->SetTriggerInputs(fTrigger);
    fESDADfriend->SetTriggerInputsMask(rawStream.GetTriggerInputsMask());
  
    for(Int_t iScaler = 0; iScaler < AliESDADfriend::kNScalers; iScaler++)
      fESDADfriend->SetTriggerScalers(iScaler,rawStream.GetTriggerScalers(iScaler));

    digitsTree->Fill();
  }

  fDigitsArray->Clear();

}

//_____________________________________________________________________________
void AliADReconstructor::FillESD(TTree* digitsTree, TTree* /*clustersTree*/,AliESDEvent* esd) const
{

  printf("Running AD Reconstruction \n");

  if (!digitsTree) {
      AliError("No digits tree!");
      return;
  }

  TBranch* digitBranch = digitsTree->GetBranch("ADDigit");
  digitBranch->SetAddress(&fDigitsArray);

  Float_t   mult[16];  
  Float_t    adc[16]; 
  Float_t   time[16]; 
  Float_t  width[16];
  Bool_t aBBflag[16];
  Bool_t aBGflag[16];
   
  for (Int_t i=0; i<16; i++){
       adc[i]    = 0.0;
       mult[i]   = 0.0;
       time[i]   = kInvalidTime;
       width[i]  = 0.0;
       aBBflag[i] = kFALSE;
       aBGflag[i] = kFALSE;
  }

  Int_t nEntries = (Int_t)digitsTree->GetEntries();
  for (Int_t e=0; e<nEntries; e++) {
    digitsTree->GetEvent(e);

    Int_t nDigits = fDigitsArray->GetEntriesFast();
    
    for (Int_t d=0; d<nDigits; d++) {    
        AliADdigit* digit = (AliADdigit*) fDigitsArray->At(d);      
        Int_t  pmNumber = digit->PMNumber();

        // Pedestal retrieval and suppression
	Bool_t integrator = digit->Integrator();
        Float_t maxadc = 0;
        Int_t imax = -1;
        Float_t adcPedSub[kADNClocks];
        for(Int_t iClock=0; iClock < kADNClocks; ++iClock) {
	  Short_t charge = digit->ChargeADC(iClock);
	  Bool_t iIntegrator = (iClock%2 == 0) ? integrator : !integrator;
	  Int_t k = pmNumber + 16*iIntegrator;
	  adcPedSub[iClock] = (Float_t)charge - fCalibData->GetPedestal(k);
	  if(adcPedSub[iClock] <= GetRecoParam()->GetNSigmaPed()*fCalibData->GetSigma(k)) {
	    adcPedSub[iClock] = 0;
	    continue;
	  }
	  
	  if(iClock < GetRecoParam()->GetStartClock() || iClock > GetRecoParam()->GetEndClock()) continue;
	  if(adcPedSub[iClock] > maxadc) {
	    maxadc = adcPedSub[iClock];
	    imax   = iClock;
	  }
	}

	if (imax != -1) {
	  Int_t start = imax - GetRecoParam()->GetNPreClocks();
	  if (start < 0) start = 0;
	  Int_t end = imax + GetRecoParam()->GetNPostClocks();
	  if (end > 20) end = 20;
	  for(Int_t iClock = start; iClock <= end; iClock++) {
	    adc[pmNumber] += adcPedSub[iClock];
	  }
	}

	// HPTDC leading time and width
	// Correction for slewing and various time delays
        time[pmNumber]  =  CorrectLeadingTime(pmNumber,digit->Time(),adc[pmNumber]);
	width[pmNumber] =  digit->Width();
	aBBflag[pmNumber] = digit->GetBBflag();
	aBGflag[pmNumber] = digit->GetBGflag();
	
	//MIP multiplicity
	if(fCalibData->GetADCperMIP(pmNumber) != 0) mult[pmNumber] = adc[pmNumber]/fCalibData->GetADCperMIP(pmNumber);

	if (adc[pmNumber] > 0) {
	  AliDebug(1,Form("PM = %d ADC = %.2f (%.2f) TDC %.2f (%.2f)   Int %d (%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d)    %.2f %.2f   %.2f %.2f    %d %d",pmNumber, adc[pmNumber],
		       digit->ChargeADC(11)+digit->ChargeADC(10)+digit->ChargeADC(9)+digit->ChargeADC(8)+
		       digit->ChargeADC(7)+digit->ChargeADC(6)+digit->ChargeADC(5)+digit->ChargeADC(4)-
		       4.*fCalibData->GetPedestal(pmNumber)-4.*fCalibData->GetPedestal(pmNumber+16),
			  digit->Time(),time[pmNumber],
			  integrator,
		          digit->ChargeADC(0),digit->ChargeADC(1),digit->ChargeADC(2),digit->ChargeADC(3),digit->ChargeADC(4),digit->ChargeADC(5),digit->ChargeADC(6),digit->ChargeADC(7),
			  digit->ChargeADC(8),digit->ChargeADC(9),digit->ChargeADC(10),
			  digit->ChargeADC(11),digit->ChargeADC(12),
		          digit->ChargeADC(13),digit->ChargeADC(14),digit->ChargeADC(15),digit->ChargeADC(16),digit->ChargeADC(17),digit->ChargeADC(18),digit->ChargeADC(19),digit->ChargeADC(20),
			  fCalibData->GetPedestal(pmNumber),fCalibData->GetSigma(pmNumber),
			  fCalibData->GetPedestal(pmNumber+16),fCalibData->GetSigma(pmNumber+16),
			  aBBflag[pmNumber],aBGflag[pmNumber]));
	    };

	// Fill ESD friend object
	for (Int_t iEv = 0; iEv < kADNClocks; iEv++) {
	  fESDADfriend->SetPedestal(pmNumber,iEv,(Float_t)digit->ChargeADC(iEv));
	  fESDADfriend->SetIntegratorFlag(pmNumber,iEv,(iEv%2 == 0) ? integrator : !integrator);
	}
	fESDADfriend->SetTime(pmNumber,digit->Time());
	fESDADfriend->SetWidth(pmNumber,digit->Width());

    } // end of loop over digits
  } // end of loop over events in digits tree
         
  fESDAD->SetBit(AliESDAD::kCorrectedLeadingTime,kTRUE);
  fESDAD->SetBit(AliESDAD::kCorrectedForSaturation,kFALSE);
  fESDAD->SetMultiplicity(mult);
  fESDAD->SetADC(adc);
  fESDAD->SetTime(time);
  fESDAD->SetWidth(width);
  fESDAD->SetBBFlag(aBBflag);
  fESDAD->SetBGFlag(aBGflag);
  
  // Fill BB and BG flags for all channel in 21 clocks (called past-future flags)
  for(Int_t i = 0; i < 16; ++i) {
    for (Int_t iClock=0; iClock < kADNClocks; ++iClock) {
      fESDAD->SetPFBBFlag(i,iClock,fESDADfriend->GetBBFlag(i,iClock));
      fESDAD->SetPFBGFlag(i,iClock,fESDADfriend->GetBGFlag(i,iClock));
    }
  }
  fESDAD->SetBit(AliESDAD::kPastFutureFlagsFilled,kTRUE);
   
  Int_t	pBBmulADA = 0;
  Int_t	pBBmulADC = 0;
  Int_t	pBGmulADA = 0;
  Int_t	pBGmulADC = 0;
  
  for(Int_t iChannel=0; iChannel<4; iChannel++) {//Loop over pairs of pads
    //Enable time is used to turn off the coincidence 
    if((!fCalibData->GetEnableTiming(iChannel) || aBBflag[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || aBBflag[iChannel+4])) pBBmulADC++;
    if((!fCalibData->GetEnableTiming(iChannel) || aBGflag[iChannel]) && (!fCalibData->GetEnableTiming(iChannel+4) || aBGflag[iChannel+4])) pBGmulADC++;

    if((!fCalibData->GetEnableTiming(iChannel+8) || aBBflag[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBBflag[iChannel+12])) pBBmulADA++;
    if((!fCalibData->GetEnableTiming(iChannel+8) || aBGflag[iChannel+8]) && (!fCalibData->GetEnableTiming(iChannel+12) || aBGflag[iChannel+12])) pBGmulADA++;
	}
	
  Bool_t UBA = kFALSE;
  Bool_t UBC = kFALSE;
  Bool_t UGA = kFALSE;
  Bool_t UGC = kFALSE;

  if(pBBmulADA>=fCalibData->GetBBAThreshold()) UBA = kTRUE;
  if(pBBmulADC>=fCalibData->GetBBCThreshold()) UBC = kTRUE;
  if(pBGmulADA>=fCalibData->GetBGAThreshold()) UGA = kTRUE;
  if(pBGmulADC>=fCalibData->GetBGCThreshold()) UGC = kTRUE;
  
  UShort_t fTrigger = 0;
  if(UBA && UBC) fTrigger |= 1;
  if(UBA || UBC) fTrigger |= (1<<1);
  if(UGA && UBC) fTrigger |= (1<<2);
  if(UGA) fTrigger |= (1<<3);
  if(UGC && UBA) fTrigger |= (1<<4);
  if(UGC) fTrigger |= (1<<5);
  
  //CTA1 and CTC1
  //CTA1 or CTC1
  //CTA2 and CTC2
  //CTA2 or CTC2 
  //MTA and MTC
  //MTA or MTC
  
  if(UBA) fTrigger |= (1<<12);
  if(UBC) fTrigger |= (1<<13); 
  if(UGA || UGC) fTrigger |= (1<<14);
  if((UGA && UBC) || (UGC && UBA)) fTrigger |= (1<<15);

  fESDAD->SetTriggerBits(fTrigger);
  fESDAD->SetBit(AliESDAD::kOnlineBitsFilled,kTRUE);

  // now fill the AD decision
  AliADDecision offlineDecision;
  offlineDecision.SetRecoParam(GetRecoParam());
  offlineDecision.FillDecisions(fESDAD);
  
  if (esd) { 
     AliDebug(1, Form("Writing AD data to ESD tree"));
     esd->SetADData(fESDAD);
 
     AliESDfriend *fr = (AliESDfriend*)esd->FindListObject("AliESDfriend");
     if (fr) {
        AliDebug(1, Form("Writing AD friend data to ESD tree"));
        fr->SetADfriend(fESDADfriend);
    }
  }

  fDigitsArray->Clear();

}

//_____________________________________________________________________________
AliADCalibData* AliADReconstructor::GetCalibData() const
{
  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("AD/Calib/Data");

  AliADCalibData *calibdata = 0;

  if (entry) calibdata = (AliADCalibData*) entry->GetObject();
  if (!calibdata)  AliFatal("No calibration data from calibration database !");

  return calibdata;
}

//_____________________________________________________________________________
void AliADReconstructor::GetTimeSlewingSplines()
{

  AliCDBManager *man = AliCDBManager::Instance();

  AliCDBEntry *entry=0;

  entry = man->Get("AD/Calib/TimeSlewing");
  
  TList *fListSplines = 0;

  if (entry) fListSplines = (TList*) entry->GetObject();
  if (!fListSplines)  AliFatal("No time slewing correction from calibration database !");
  
  for(Int_t i=0; i<16; i++) fTimeSlewingSpline[i] = (TSpline3*)(fListSplines->At(i));
  

}

//_____________________________________________________________________________
AliCDBStorage* AliADReconstructor::SetStorage(const char *uri) 
{
// Sets the storage  

  Bool_t deleteManager = kFALSE;
  
  AliCDBManager *manager = AliCDBManager::Instance();
  AliCDBStorage *defstorage = manager->GetDefaultStorage();
  
  if(!defstorage || !(defstorage->Contains("AD"))){ 
     AliWarning("No default storage set or default storage doesn't contain AD!");
     manager->SetDefaultStorage(uri);
     deleteManager = kTRUE;
  }
 
  AliCDBStorage *storage = manager->GetDefaultStorage();

  if(deleteManager){
     AliCDBManager::Instance()->UnsetDefaultStorage();
     defstorage = 0;   // the storage is killed by AliCDBManager::Instance()->Destroy()
  }

  return storage; 
}

//____________________________________________________________________________
void AliADReconstructor::GetCollisionMode()
{
  // Retrieval of collision mode 

  TString beamType = GetRunInfo()->GetBeamType();
  if(beamType==AliGRPObject::GetInvalidString()){
     AliError("AD cannot retrieve beam type");
     return;
  }

  if( (beamType.CompareTo("P-P") ==0)  || (beamType.CompareTo("p-p") ==0) ){
    fCollisionMode=0;
  }
  else if( (beamType.CompareTo("Pb-Pb") ==0)  || (beamType.CompareTo("A-A") ==0) ){
    fCollisionMode=1;
  }
    
  fBeamEnergy = GetRunInfo()->GetBeamEnergy();
  if(fBeamEnergy==AliGRPObject::GetInvalidFloat()) {
     AliError("Missing value for the beam energy ! Using 0");
     fBeamEnergy = 0.;
  }
  
  AliDebug(1,Form("\n ++++++ Beam type and collision mode retrieved as %s %d @ %1.3f GeV ++++++\n\n",beamType.Data(), fCollisionMode, fBeamEnergy));

}
//____________________________________________________________________________
Float_t AliADReconstructor::CorrectLeadingTime(Int_t i, Float_t time, Float_t adc) const
{
  // Correct the leading time
  // for slewing effect and
  // misalignment of the channels
  if (time < 1e-6) return kInvalidTime;

  // In case of pathological signals
  if (adc < 1) return time;

  // Slewing and offset correction
  Int_t board = AliADCalibData::GetBoardNumber(i);
  time -= fTimeSlewingSpline[i]->Eval(TMath::Log10(1/adc))*fCalibData->GetTimeResolution(board);
  time += fLayerDist[i/4];
  
  // Channel alignment and general offset subtraction
  //time -= fHptdcOffset[i];

  return time;
}
