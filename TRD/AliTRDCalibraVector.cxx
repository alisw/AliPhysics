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

////////////////////////////////////////////////////////////////////////////
//                                                                        //   
// AliTRDCalibraVector                                                    //       
//                                                                        //   
// This class is for the vector method of the TRD calibration.            //
//                                                                        //
// Author:                                                                //
//   R. Bailhache (R.Bailhache@gsi.de)                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TGraphErrors.h>
#include <TH1F.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TMath.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TFile.h>
#include <TString.h>

#include "AliLog.h"

#include "AliTRDCalibraVector.h"
#include "AliTRDCommonParam.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDPhInfo.h"
#include "AliTRDEntriesInfo.h"
#include "AliTRDPrfInfo.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDCalibraVector)

//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector()
  :TObject()
  ,fModeCH(0)
  ,fModePH(0)
  ,fModePRF(0)
  ,fNbGroupPRF(0)
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fStopFillCH(kFALSE)
  ,fHisto(0x0)
  ,fGraph(0x0)
  ,fCalVector(0x0)
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fTimeMax(0)
  ,fPRFRange(1.5)
{
  //
  // Default constructor
  //

  for (Int_t idet = 0; idet < 540; idet++){
    
    fPHEntries[idet]= 0x0;
    fPHMean[idet]= 0x0;
    fPHSquares[idet]= 0x0;

    fPRFEntries[idet]= 0x0;
    fPRFMean[idet]= 0x0;
    fPRFSquares[idet]= 0x0;


    fCHEntries[idet]= 0x0;
    
  }
  
  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = 0;
    fDetCha2[k] = 0;
  }
 
}
//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector(const AliTRDCalibraVector &c)
  :TObject(c)
  ,fModeCH(c.fModeCH)
  ,fModePH(c.fModePH)
  ,fModePRF(c.fModePRF)
  ,fNbGroupPRF(c.fNbGroupPRF)
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fStopFillCH(kFALSE)
  ,fHisto(0x0)
  ,fGraph(0x0)
  ,fCalVector(0x0)
  ,fNumberBinCharge(c.fNumberBinCharge)
  ,fNumberBinPRF(c.fNumberBinPRF)
  ,fTimeMax(c.fTimeMax)
  ,fPRFRange(c.fPRFRange)
{
  //
  // Copy constructor
  //
  
  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = c.fDetCha0[k];
    fDetCha2[k] = c.fDetCha2[k];
  }

  for (Int_t idet = 0; idet < 540; idet++){
    
    const AliTRDEntriesInfo *phEntries  = (AliTRDEntriesInfo*)c.fPHEntries[idet];
    const AliTRDPhInfo      *phMean     = (AliTRDPhInfo *)c.fPHMean[idet];
    const AliTRDPhInfo      *phSquares  = (AliTRDPhInfo *)c.fPHSquares[idet];

    const AliTRDEntriesInfo *prfEntries = (AliTRDEntriesInfo*)c.fPRFEntries[idet];
    const AliTRDPrfInfo     *prfMean     = (AliTRDPrfInfo *)c.fPRFMean[idet];
    const AliTRDPrfInfo     *prfSquares  = (AliTRDPrfInfo *)c.fPRFSquares[idet];

    const AliTRDEntriesInfo *chEntries  = (AliTRDEntriesInfo*)c.fCHEntries[idet];

    if ( chEntries != 0x0 ) fCHEntries[idet] = new AliTRDEntriesInfo(*chEntries);
    
    if ( phEntries != 0x0 ) {
      fPHMean[idet]    = new AliTRDPhInfo(*phMean);
      fPHSquares[idet] = new AliTRDPhInfo(*phSquares);
      fPHEntries[idet] = new AliTRDEntriesInfo(*phEntries);
    }

    if ( prfEntries != 0x0 ) {
      fPRFEntries[idet] = new AliTRDEntriesInfo(*prfEntries);
      fPRFMean[idet] = new AliTRDPrfInfo(*prfMean);
      fPRFSquares[idet] = new AliTRDPrfInfo(*prfSquares);
    }
    
  }
   
}
//_____________________________________________________________________
AliTRDCalibraVector& AliTRDCalibraVector::operator = (const  AliTRDCalibraVector &source)
{
  //
  // assignment operator
  //
  if (&source == this) return *this;
  new (this) AliTRDCalibraVector(source);

  return *this;
}
//____________________________________________________________________________________
AliTRDCalibraVector::~AliTRDCalibraVector()
{
  //
  // AliTRDCalibraVector destructor
  //

  if(fPHEntries) delete fPHEntries;
  if(fPHMean) delete fPHMean;
  if(fPHSquares) delete fPHSquares;
  if(fPRFEntries) delete fPRFEntries;
  if(fPRFMean) delete fPRFMean;
  if(fPRFSquares) delete fPRFSquares;
  if(fCHEntries) delete fCHEntries;
  if(fHisto) delete fHisto;
  if(fGraph) delete fGraph;
  if(fCalVector) delete fCalVector;

}
//_____________________________________________________________________________
Long64_t AliTRDCalibraVector::Merge(const TCollection* list) 
{
  // Merge list of objects (needed by PROOF)

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;
  
  TIterator* iter = list->MakeIterator();
  TObject* obj = 0;
  
  // collection of generated histograms
  Int_t count=0;
  while((obj = iter->Next()) != 0) 
    {
      AliTRDCalibraVector* entry = dynamic_cast<AliTRDCalibraVector*>(obj);
      if (entry == 0) continue; 
      
      if(this->Add(entry)) count++;
    
    }
  
  return count;
}
//_____________________________________________________________________________
void AliTRDCalibraVector::TestInit(Int_t i, Int_t detmax)
{
  //
  // Init to see the size
  //

  for(Int_t det = 0; det < detmax; det++){

    if(i==2) {
      
      fPRFEntries[det] = ((AliTRDEntriesInfo *)GetPRFEntries(det,kTRUE));
      fPRFMean[det]    = ((AliTRDPrfInfo  *)GetPRFMean(det,kTRUE));
      fPRFSquares[det] = ((AliTRDPrfInfo  *)GetPRFSquares(det,kTRUE));

    }

    if(i==1) {

      fPHEntries[det] = ((AliTRDEntriesInfo *)GetPHEntries(det,kTRUE));
      fPHMean[det]    = ((AliTRDPhInfo  *)GetPHMean(det,kTRUE));
      fPHSquares[det] = ((AliTRDPhInfo  *)GetPHSquares(det,kTRUE));
      
    }
    
    if(i==0) fCHEntries[det] = ((AliTRDEntriesInfo *)GetCHEntries(det,kTRUE));

  }

}
//_____________________________________________________________________________
Int_t AliTRDCalibraVector::SearchBin(Float_t value, Int_t i) const
{
  //
  // Search the bin
  //

  Int_t reponse        = 0;
  Float_t fbinmin      = 0;
  Float_t fbinmax      = value;
  Int_t fNumberOfBin   = -1;

  switch(i)
    {
    case 0:
      fbinmax      = 300.0;
      fbinmin      = 0.0;
      fNumberOfBin = fNumberBinCharge;
      break;

    case 2:
      fbinmax      =   TMath::Abs(fPRFRange);
      fbinmin      =  -TMath::Abs(fPRFRange);
      fNumberOfBin = fNumberBinPRF;
      break;
      
    default: 
      return -1;
    }
  
  // Return -1 if out
  if ((value >= fbinmax) || 
      (value <  fbinmin)) {
    return -1;
  }
  else {
    reponse = (Int_t) ((fNumberOfBin*(value-fbinmin)) / (fbinmax-fbinmin));
  }

  return reponse;

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorCH(Int_t det, Int_t group, Float_t value)
{
  //
  // Fill the vector CH   
  //

  // Search bin
  Int_t bin = SearchBin(value,0);
  // Out
  if (bin == -1) {
    return kFALSE; 
  }



  if(fDetectorCH != det){
    fCHEntries[det] = ((AliTRDEntriesInfo *)GetCHEntries(det,kTRUE));
  }

  Int_t entries  = ((AliTRDEntriesInfo *)fCHEntries[det])->At(group*fNumberBinCharge+bin);
  
  Int_t entriesn = entries+1;

  if(entriesn > 65535) {
    fStopFillCH = kTRUE;
    return kTRUE;
  }

  ((AliTRDEntriesInfo *)fCHEntries[det])->AddAt(entriesn,group*fNumberBinCharge+bin);
  
  fDetectorCH = det;

 
  return kTRUE;

}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorPRF(Int_t det, Int_t group, Float_t x, Float_t y)
{
  //
  // Fill the vector PRF
  //

  // Search bin
  Int_t bin = SearchBin(x,2);
  // Out
  if (bin == -1) {
    return kFALSE; 
  }

  
  if(fDetectorPRF != det){
    fPRFEntries[det] = ((AliTRDEntriesInfo *)GetPRFEntries(det,kTRUE));
    fPRFMean[det]    = ((AliTRDPrfInfo  *)GetPRFMean(det,kTRUE));
    fPRFSquares[det] = ((AliTRDPrfInfo  *)GetPRFSquares(det,kTRUE));
  }

  Int_t entries  = ((AliTRDEntriesInfo *)fPRFEntries[det])->At((Int_t)group*fNumberBinPRF+bin);
  Float_t mean   = ((AliTRDPrfInfo  *)fPRFMean[det])->At((Int_t)group*fNumberBinPRF+bin);
  Float_t square = ((AliTRDPrfInfo  *)fPRFSquares[det])->At((Int_t)group*fNumberBinPRF+bin);
  
  Int_t entriesn = entries+1;

  if(entriesn > 65535) return kTRUE;
  
  ((AliTRDEntriesInfo *)fPRFEntries[det])->AddAt(entriesn,(Int_t)group*fNumberBinPRF+bin);
  Float_t meann = (mean*((Float_t)entries)+y)/((Float_t)entriesn);
  ((AliTRDPrfInfo *)fPRFMean[det])->AddAt(meann,(Int_t)group*fNumberBinPRF+bin);
  Float_t squaren = ((square*((Float_t)entries))+(y*y))/((Float_t)entriesn);
  ((AliTRDPrfInfo *)fPRFSquares[det])->AddAt(squaren,(Int_t)group*fNumberBinPRF+bin);

  
  fDetectorPRF = det;

  return kTRUE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorPH(Int_t det, Int_t group, Int_t time, Float_t value)
{
  //
  // Fill the vector PH  
  //

  // Search bin
  Int_t bin = time;
  // Out
  if ((bin <         0) || 
      (bin >= fTimeMax)) {
    return kFALSE; 
  }


  if(fDetectorPH != det){
    fPHEntries[det] = ((AliTRDEntriesInfo *)GetPHEntries(det,kTRUE));
    fPHMean[det]    = ((AliTRDPhInfo  *)GetPHMean(det,kTRUE));
    fPHSquares[det] = ((AliTRDPhInfo  *)GetPHSquares(det,kTRUE));
  }

  Int_t entries  = ((AliTRDEntriesInfo *)fPHEntries[det])->At(group*fTimeMax+bin);
  Float_t mean   = ((AliTRDPhInfo  *)fPHMean[det])->At(group*fTimeMax+bin);
  Float_t square = ((AliTRDPhInfo  *)fPHSquares[det])->AtS(group*fTimeMax+bin);
  
  Int_t entriesn = entries+1;
  Float_t meann = (mean*((Float_t)entries)+value)/((Float_t)entriesn);
  Float_t squaren = ((square*((Float_t)entries))+(value*value))/((Float_t)entriesn);
  
  if(entriesn > 65535) return kTRUE;
  //printf("meann %f, squaren %f\n",meann,squaren);
  if((meann > 3000.0) || (meann < 0.0) || (squaren > (3000.0*3000.0)) || (squaren < 0.0)) return kFALSE;

  ((AliTRDEntriesInfo *)fPHEntries[det])->AddAt(entriesn,group*fTimeMax+bin);
  ((AliTRDPhInfo *)fPHMean[det])->AddAt(meann,group*fTimeMax+bin);
  ((AliTRDPhInfo *)fPHSquares[det])->AddAtS(squaren,group*fTimeMax+bin);
  
  fDetectorPH = det;

  return kTRUE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::FillVectorCH(Int_t det, Int_t group, Int_t bin, Int_t entries)
{
  //
  // Fill the vector CH   
  //

  if(entries > 65535) return kFALSE;

  if(fDetectorCH != det){
    fCHEntries[det] = ((AliTRDEntriesInfo *)GetCHEntries(det,kTRUE));
  }

  ((AliTRDEntriesInfo *)fCHEntries[det])->AddAt(entries,group*fNumberBinCharge+bin);
  
  fDetectorCH = det;
  
  
  return kTRUE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::FillVectorPRF(Int_t det, Int_t group, Int_t bin, Int_t entries, Float_t mean, Float_t square)
{
  //
  // Fill the vector PRF
  //

  if((entries > 65535) || (mean > 1.0) || (mean < 0.0) || (square > 1.0) || (square < 0.0)) return kFALSE;

  if(fDetectorPRF != det){
    fPRFEntries[det] = ((AliTRDEntriesInfo *)GetPRFEntries(det,kTRUE));
    fPRFMean[det]    = ((AliTRDPrfInfo  *)GetPRFMean(det,kTRUE));
    fPRFSquares[det] = ((AliTRDPrfInfo  *)GetPRFSquares(det,kTRUE));
  }

  ((AliTRDEntriesInfo *)fPRFEntries[det])->AddAt(entries,(Int_t)group*fNumberBinPRF+bin);
  ((AliTRDPrfInfo  *)fPRFMean[det])->AddAt(mean,(Int_t)group*fNumberBinPRF+bin);
  ((AliTRDPrfInfo  *)fPRFSquares[det])->AddAt(square,(Int_t)group*fNumberBinPRF+bin);

  
  fDetectorPRF = det;

  return kTRUE;
  
}
//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::FillVectorPH(Int_t det, Int_t group, Int_t bin, Int_t entries, Float_t mean, Float_t square)
{
  //
  // Fill the vector PH  
  //

  if((entries > 65535) || (mean > 3000.0) || (mean < 0.0) || (square > (3000.0*3000.0)) || (square < 0.0)) return kFALSE;


  if(fDetectorPH != det){
    fPHEntries[det] = ((AliTRDEntriesInfo *)GetPHEntries(det,kTRUE));
    fPHMean[det]    = ((AliTRDPhInfo  *)GetPHMean(det,kTRUE));
    fPHSquares[det] = ((AliTRDPhInfo  *)GetPHSquares(det,kTRUE));
  }

  ((AliTRDEntriesInfo *)fPHEntries[det])->AddAt(entries,group*fTimeMax+bin);
  ((AliTRDPhInfo  *)fPHMean[det])->AddAt(mean,group*fTimeMax+bin);
  ((AliTRDPhInfo  *)fPHSquares[det])->AddAtS(square,group*fTimeMax+bin);
  
  fDetectorPH = det;

  return kTRUE;
  
}
//__________________________________________________________________________________
Bool_t AliTRDCalibraVector::Add(AliTRDCalibraVector *calvect)
{
  //
  // Add a other AliTRCalibraVector to this one
  //

  Bool_t result = kTRUE;

  // Check compatibility
  if(fNumberBinCharge != calvect->GetNumberBinCharge()) return kFALSE;
  if(fNumberBinPRF    != calvect->GetNumberBinPRF()) return kFALSE;
  if(fPRFRange        != calvect->GetPRFRange()) return kFALSE;
  if(fTimeMax         != calvect->GetTimeMax()) return kFALSE;
  for(Int_t k = 0; k < 3; k++){
    if(fDetCha0[k] != calvect->GetDetCha0(k)) return kFALSE;
    if(fDetCha2[k] != calvect->GetDetCha2(k)) return kFALSE;
  }

  //printf("All ok for variables before adding!\n"); 

  // Add
  for (Int_t idet = 0; idet < 540; idet++){

    //printf("Detector %d\n",idet);
    
    const AliTRDEntriesInfo *phEntriesvect    = (AliTRDEntriesInfo *) calvect->GetPHEntries(idet);
    const AliTRDPhInfo      *phMeanvect       = (AliTRDPhInfo *) calvect->GetPHMean(idet);
    const AliTRDPhInfo      *phSquaresvect    = (AliTRDPhInfo *) calvect->GetPHSquares(idet);
    
    const AliTRDEntriesInfo *prfEntriesvect   = (AliTRDEntriesInfo *) calvect->GetPRFEntries(idet);
    const AliTRDPrfInfo    *prfMeanvect       = (AliTRDPrfInfo *) calvect->GetPRFMean(idet);
    const AliTRDPrfInfo    *prfSquaresvect    = (AliTRDPrfInfo *) calvect->GetPRFSquares(idet);
    
    const AliTRDEntriesInfo *chEntriesvect    = (AliTRDEntriesInfo *) calvect->GetCHEntries(idet);

    if ( phEntriesvect != 0x0 ){
      //Take the stuff
      fPHEntries[idet] = ((AliTRDEntriesInfo *)GetPHEntries(idet,kTRUE));
      fPHMean[idet]    = ((AliTRDPhInfo *)GetPHMean(idet,kTRUE));
      fPHSquares[idet] = ((AliTRDPhInfo *)GetPHSquares(idet,kTRUE));

      Int_t total = fPHEntries[idet]->GetSize();
      //printf("Total size PH %d\n",total);
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entriesv  = ((AliTRDEntriesInfo *)phEntriesvect)->At(k);
	Float_t meanv   = ((AliTRDPhInfo  *)phMeanvect)->At(k);
	Float_t squarev = ((AliTRDPhInfo  *)phSquaresvect)->AtS(k);
	
	Int_t entries  = ((AliTRDEntriesInfo *)fPHEntries[idet])->At(k);
	Float_t mean   = ((AliTRDPhInfo  *)fPHMean[idet])->At(k);
	Float_t square = ((AliTRDPhInfo  *)fPHSquares[idet])->AtS(k);
  
	Int_t entriesn = entries+entriesv;
	Float_t meann = (mean*((Float_t)entries)+meanv*((Float_t)entriesv))/((Float_t)entriesn);
	Float_t squaren = ((square*((Float_t)entries))+(squarev*((Float_t)entriesv)))/((Float_t)entriesn);
	
	if((entriesn > 0) && (entriesn <= 65535) && (meann >= 0) && (meann < 3000.0) && (squaren >= 0.0) && (squaren < (3000.0*3000.0))) {
	
	  ((AliTRDEntriesInfo *)fPHEntries[idet])->AddAt(entriesn,k);
	  ((AliTRDPhInfo *)fPHMean[idet])->AddAt(meann,k);
	  ((AliTRDPhInfo *)fPHSquares[idet])->AddAtS(squaren,k);
      
	}
      }
    }     

    if ( prfEntriesvect != 0x0 ){
      //Take the stuff
      fPRFEntries[idet] = ((AliTRDEntriesInfo *)GetPRFEntries(idet,kTRUE));
      fPRFMean[idet]    = ((AliTRDPrfInfo  *)GetPRFMean(idet,kTRUE));
      fPRFSquares[idet] = ((AliTRDPrfInfo  *)GetPRFSquares(idet,kTRUE));
      
      Int_t total = fPRFEntries[idet]->GetSize(); 
      //Int_t total0 = fPRFMean[idet]->GetSize(); 
      //Int_t total1 = fPRFSquares[idet]->GetSize(); 
      //printf("Total size PRF %d\n",total);     
      //printf("Total0 size PRF %d\n",total0);     
      //printf("Total1 size PRF %d\n",total1);     
      // Add
      for(Int_t k = 0; k < total; k++){


	Int_t entries  = ((AliTRDEntriesInfo *)fPRFEntries[idet])->At(k);
	Float_t mean   = ((AliTRDPrfInfo  *)fPRFMean[idet])->At(k);
	Float_t square = ((AliTRDPrfInfo  *)fPRFSquares[idet])->At(k);

	Int_t entriesv  = ((AliTRDEntriesInfo *)prfEntriesvect)->At(k);
	Float_t meanv   = ((AliTRDPrfInfo  *)prfMeanvect)->At(k);
	Float_t squarev = ((AliTRDPrfInfo  *)prfSquaresvect)->At(k);

	Int_t entriesn = entries + entriesv;
	
	if((entriesn > 0) && (entriesn <= 65535)) {
	  
	  ((AliTRDEntriesInfo *)fPRFEntries[idet])->AddAt(entriesn,k);
	  
	  Float_t meann = (mean*((Float_t)entries)+meanv*((Float_t)entriesv))/((Float_t)entriesn);
	  //printf("test0\n");
	  ((AliTRDPrfInfo *)fPRFMean[idet])->AddAt(meann,k);
	  //printf("test1\n");
	  
	  Float_t squaren = ((square*((Float_t)entries))+(squarev*((Float_t)entriesv)))/((Float_t)entriesn);
	  //printf("test2\n");
	  ((AliTRDPrfInfo *)fPRFSquares[idet])->AddAt(squaren,k);
	  //printf("test3\n");
		  
	}
      }
    }

    if ( chEntriesvect != 0x0 ){
      //Take the stuff
      fCHEntries[idet] = ((AliTRDEntriesInfo *)GetCHEntries(idet,kTRUE));
      //printf("TestAdd\n");
      fStopFillCH = ((AliTRDEntriesInfo *)fCHEntries[idet])->TestAdd((AliTRDEntriesInfo *)chEntriesvect);
      //
      if(!fStopFillCH) {
	fStopFillCH = kTRUE;
	result = kFALSE;
      }
      else {
	
	((AliTRDEntriesInfo *)fCHEntries[idet])->Add(chEntriesvect);
	//printf("Add\n");
      }
    }           
  }
  
  return result;
}
//_____________________________________________________________________________________________________________________
AliTRDCalibraVector *AliTRDCalibraVector::AddStatsPerDetectorCH()
{
  //
  // Create a AliTRDCalibraVector detector wise
  //

  // Use a AliTRDCalibraMode to navigate in the calibration groups
  AliTRDCalibraMode calibMode = AliTRDCalibraMode();
  calibMode.SetNz(0,GetNz(0));
  calibMode.SetNrphi(0,GetNrphi(0));
  if(((calibMode.GetNz(0) == 100) && (calibMode.GetNrphi(0) == 100)) || ((calibMode.GetNz(0) == 10) && (calibMode.GetNrphi(0) == 10))) return 0x0;
  
  Int_t nybins = 6*4*18*fDetCha0[0]+6*18*fDetCha2[0];
  Int_t nxbins = fNumberBinCharge;
  
  // Check 
  Int_t perChamber2 = 0;
  Int_t perChamber0 = 0;
  calibMode.ModePadCalibration(2,0);
  calibMode.ModePadFragmentation(0,2,0,0);
  calibMode.SetDetChamb2(0);
  perChamber2 = (Int_t) calibMode.GetDetChamb2(0);
  calibMode.ModePadCalibration(0,0);
  calibMode.ModePadFragmentation(0,0,0,0);
  calibMode.SetDetChamb0(0);
  perChamber0 = (Int_t) calibMode.GetDetChamb0(0);
  if(nybins != (6*18*perChamber2+6*4*18*perChamber0)) return 0x0;
    
  // Create calvector 
  if(!fCalVector) fCalVector = new AliTRDCalibraVector();
  else{ 
    fCalVector->~AliTRDCalibraVector();
    new(fCalVector) AliTRDCalibraVector();
  }
  fCalVector->SetNumberBinCharge(nxbins);
  fCalVector->SetDetCha0(0,1);
  fCalVector->SetDetCha2(0,1);
  fCalVector->SetNzNrphi(0,0,0);
 
  
  for(Int_t det = 0; det < 540; det++){
    
    // Take
    AliTRDEntriesInfo *entriesch = (AliTRDEntriesInfo *)GetCHEntries(det,kFALSE);
    if(!entriesch) continue;  
  
    // Number of groups
    Int_t numberofgroup = 0;
    if(AliTRDgeometry::GetStack(det) == 2) numberofgroup = perChamber2;
    else numberofgroup = perChamber0;
  
    // Check if one can merge calibration groups for this chamber
    // entries is the number of entries in each bin after adding the different the calibration group in the detector
    fStopFillCH = kFALSE;
    Int_t firstnumberofgroup = -1;
    Int_t entries[500];
    for(Int_t k = 0; k < nxbins; k++){
      entries[k] = 0;
    }
    // Loop over group in the detector
    for(Int_t k = 0; k < numberofgroup; k++){
      // Loop over bins
      for(Int_t nx = 0; nx < nxbins; nx++) {
	Int_t binnumber = k*nxbins+nx;
	entries[nx] += entriesch->At(binnumber);
	// as soon as one bin is over threshold stop 
	if(!fStopFillCH) {
	  if(entries[nx] > 65535) {
	    firstnumberofgroup = k;
	    fStopFillCH = kTRUE;
	  }
	}
	else continue;
      }
    }
    if(fStopFillCH && (firstnumberofgroup == 0)) return 0x0;
    if(!fStopFillCH) firstnumberofgroup = numberofgroup;
    
    // Now add up to possible 
    for(Int_t k = 0; k < nxbins; k++){
      entries[k] = 0;
    }
    for(Int_t k = 0; k < firstnumberofgroup; k++){
      for(Int_t nx = 0; nx < nxbins; nx++) {
	Int_t binnumber = k*nxbins+nx;
	entries[nx] += entriesch->At(binnumber);
      }
    }

    // Finally fill
    for(Int_t nx = 0; nx < nxbins; nx++){
      fCalVector->FillVectorCH(det,0,nx,(Int_t)entries[nx]);  
    }
  }
  
  return fCalVector;
} 
//_____________________________________________________________________________________________________________________
AliTRDCalibraVector *AliTRDCalibraVector::AddStatsPerDetectorPH()
{
  //
  // Create a AliTRDCalibraVector detector wise
  //

  AliTRDCalibraMode calibMode = AliTRDCalibraMode();
  calibMode.SetNz(1,GetNz(1));
  calibMode.SetNrphi(1,GetNrphi(1));
  if(((calibMode.GetNz(1) == 100) && (calibMode.GetNrphi(1) == 100)) || ((calibMode.GetNz(1) == 10) && (calibMode.GetNrphi(1) == 10))) return 0x0;
  

  // Check
  Int_t nybins = 6*4*18*fDetCha0[1]+6*18*fDetCha2[1];
  Int_t nxbins = fTimeMax;
 
  Int_t perChamber2 = 0;
  Int_t perChamber0 = 0;
  calibMode.ModePadCalibration(2,1);
  calibMode.ModePadFragmentation(0,2,0,1);
  calibMode.SetDetChamb2(1);
  perChamber2 = (Int_t) calibMode.GetDetChamb2(1);
  calibMode.ModePadCalibration(0,1);
  calibMode.ModePadFragmentation(0,0,0,1);
  calibMode.SetDetChamb0(1);
  perChamber0 = (Int_t) calibMode.GetDetChamb0(1);
  
  if(nybins != (6*18*perChamber2+6*4*18*perChamber0)) return 0x0;
  
  // Create calvector 
  if(!fCalVector) fCalVector = new AliTRDCalibraVector();
  else{ 
    fCalVector->~AliTRDCalibraVector();
    new(fCalVector) AliTRDCalibraVector();
  }
  fCalVector->SetTimeMax(nxbins);
  fCalVector->SetDetCha0(1,1);
  fCalVector->SetDetCha2(1,1);
  fCalVector->SetNzNrphi(1,0,0);
  
 
  for(Int_t det = 0; det < 540; det++){
    
    // Take
    AliTRDEntriesInfo *entriesph = (AliTRDEntriesInfo *)GetPHEntries(det,kFALSE);
    if(!entriesph) continue;
    AliTRDPhInfo      *meanph    = (AliTRDPhInfo *)GetPHMean(det,kFALSE);
    AliTRDPhInfo      *squaresph = (AliTRDPhInfo *)GetPHSquares(det,kFALSE);

    // Number of groups
    Int_t numberofgroup = 0;
    if(AliTRDgeometry::GetStack(det) == 2) numberofgroup = perChamber2;
    else numberofgroup = perChamber0;
    
    // PH
    for(Int_t nx = 0; nx < nxbins; nx++) {
      
      Double_t entries = 0.0;
      Double_t sumw2   = 0.0;
      Double_t sumw    = 0.0;

      // Sum the contributions of the different calibration group in the detector      
      for(Int_t k = 0; k < numberofgroup; k++){
	  
	Int_t binnumber = k*nxbins+nx;	  
	
	Int_t entriesv  = ((AliTRDEntriesInfo *)entriesph)->At(binnumber);
       	Float_t sumw2v  = ((AliTRDPhInfo *)squaresph)->AtS(binnumber)*entriesv;
	Float_t sumwv   = ((AliTRDPhInfo *)meanph)->At(binnumber)*entriesv;
	
	
	if(((entries+entriesv) > 65535) || ((entries+entriesv) <= 0)) continue;

	entries += entriesv;
	sumw2   += sumw2v;
	sumw    += sumwv;
      
      }

      if(entries > 0) {
	sumw2 = sumw2/((Float_t)entries);
	sumw  = sumw/((Float_t)entries);
      }
      
      fCalVector->FillVectorPH(det,0,nx,(Int_t)entries,(Float_t)sumw,(Float_t)sumw2);
    }
  }

  return fCalVector;
  
} 
//_____________________________________________________________________________________________________________________
AliTRDCalibraVector *AliTRDCalibraVector::AddStatsPerDetectorPRF()
{
  //
  // Create a AliTRDCalibraVector detector wise
  //

  AliTRDCalibraMode calibMode = AliTRDCalibraMode();
  calibMode.SetNz(2,GetNz(2));
  calibMode.SetNrphi(2,GetNrphi(2));
  if(((calibMode.GetNz(2) == 100) && (calibMode.GetNrphi(2) == 100)) || ((calibMode.GetNz(2) == 10) && (calibMode.GetNrphi(2) == 10))) return 0x0;

  // Check  
  Int_t nybins =  6*4*18*fDetCha0[2]+ 6*18*fDetCha2[2];
  Int_t nxbins = fNumberBinPRF;
  
  Int_t perChamber2 = 0;
  Int_t perChamber0 = 0;
  calibMode.ModePadCalibration(2,2);
  calibMode.ModePadFragmentation(0,2,0,2);
  calibMode.SetDetChamb2(2);
  perChamber2 = (Int_t) calibMode.GetDetChamb2(2);
  calibMode.ModePadCalibration(0,2);
  calibMode.ModePadFragmentation(0,0,0,2);
  calibMode.SetDetChamb0(2);
  perChamber0 = (Int_t) calibMode.GetDetChamb0(2);
  
  if(nybins != (6*18*perChamber2+6*4*18*perChamber0)) return 0x0;
    
  // Create calvector 
  if(!fCalVector) fCalVector = new AliTRDCalibraVector();
  else{ 
    fCalVector->~AliTRDCalibraVector();
    new(fCalVector) AliTRDCalibraVector();
  }
  fCalVector->SetNumberBinPRF(nxbins);
  fCalVector->SetDetCha0(2,1);
  fCalVector->SetDetCha2(2,1);
  fCalVector->SetNzNrphi(2,0,0);
  fCalVector->SetNbGroupPRF(fNbGroupPRF);

  
  for(Int_t det = 0; det < 540; det++){
    
    // Take
    AliTRDEntriesInfo *entriesprf = (AliTRDEntriesInfo *) GetPRFEntries(det,kFALSE);
    if(!entriesprf) continue;
    AliTRDPrfInfo     *meanprf    = (AliTRDPrfInfo *) GetPRFMean(det,kFALSE);
    AliTRDPrfInfo     *squaresprf = (AliTRDPrfInfo *) GetPRFSquares(det,kFALSE);

    // Number of groups
    Int_t numberofgroup = 0;
    if(AliTRDgeometry::GetStack(det) == 2) numberofgroup = perChamber2;
    else numberofgroup = perChamber0;
    
    for(Int_t nx = 0; nx < nxbins; nx++) {
      
      Double_t entries = 0.0;
      Double_t sumw2   = 0.0;
      Double_t sumw    = 0.0;
      
      // Sum the contributions of the different groups in the detector for one bin
      for(Int_t k = 0; k < numberofgroup; k++){
	  
	Int_t binnumber = k*nxbins+nx;	  

	Int_t entriesv  = ((AliTRDEntriesInfo *)entriesprf)->At(binnumber);
       	Float_t sumw2v  = ((AliTRDPrfInfo *)squaresprf)->At(binnumber)*entriesv;
	Float_t sumwv   = ((AliTRDPrfInfo *)meanprf)->At(binnumber)*entriesv;
	
	if(((entries+entriesv) > 65535) || ((entries+entriesv) <= 0)) continue;

	entries += entriesv;
	sumw2   += sumw2v;
	sumw    += sumwv;
      
      }

      if(entries > 0) {
	sumw2 = sumw2/((Float_t)entries);
	sumw  = sumw/((Float_t)entries);
      }
      
      fCalVector->FillVectorPRF(det,0,nx,(Int_t)entries,(Float_t)sumw,(Float_t)sumw2);
      
    }
  }

  return fCalVector;
}
//_______________________________________________________________________________
Bool_t AliTRDCalibraVector::FindTheMaxEntries(Int_t i, Int_t &detectormax, Int_t &groupmax)
{
  //
  // Find detectormax and groupmax with the biggest number of entries
  //

  Int_t numberofTB = 0;
  if(i==0) numberofTB = (Int_t) GetNumberBinCharge();
  if(i==1) numberofTB = GetTimeMax();
  if(i==2) numberofTB = GetNumberBinPRF();
  if((i!=0) && (i!=1) && (i!=2)) AliInfo("Didn't understand i");


  // Init
  Double_t entries [540];
  for(Int_t idet = 0; idet < 540; idet++){
    entries[idet] = 0.0;
  }

  AliTRDEntriesInfo *entriesd = 0x0;
  // Take the number of entries per detector
  for(Int_t idet = 0; idet < 540; idet++){
 
    if(i==0) entriesd = (AliTRDEntriesInfo *) GetCHEntries(idet,kFALSE);
    if(i==1) entriesd = (AliTRDEntriesInfo *) GetPHEntries(idet,kFALSE);
    if(i==2) entriesd = (AliTRDEntriesInfo *) GetPRFEntries(idet,kFALSE);
    if(!entriesd) continue;

    entries[idet] = entriesd->GetSum();
    
  }

  // Search detector max
  Double_t max = -10;
  detectormax = -1;
  for(Int_t idet = 0; idet < 540; idet++){
    if(entries[idet] > max) {
      max = entries[idet];
      detectormax = idet;
    }
  }
  if((max == 0.0) || (detectormax <0.0)) return kFALSE;

  // Search group max
  if(i==0) entriesd = (AliTRDEntriesInfo *) GetCHEntries(detectormax,kFALSE);
  if(i==1) entriesd = (AliTRDEntriesInfo *) GetPHEntries(detectormax,kFALSE);
  if(i==2) entriesd = (AliTRDEntriesInfo *) GetPRFEntries(detectormax,kFALSE);  
  if(!entriesd) return kFALSE;
  // Number of groups
  Int_t numberofgroup = 0;
  if(AliTRDgeometry::GetStack(detectormax) == 2) numberofgroup = fDetCha2[i];
  else numberofgroup = fDetCha0[i];
  // Init
  Double_t nbgroup [2304];
  for(Int_t k = 0; k < 2304; k++){
    nbgroup[k] = 0.0;
  }
  Int_t nxbins = 0;
  if(i==0) nxbins = fNumberBinCharge;
  if(i==1) nxbins = fTimeMax;
  if(i==2) nxbins = fNumberBinPRF;
  // Compute the number of entries per group
  for(Int_t k = 0; k < numberofgroup; k++){
    for(Int_t nx = 0; nx < nxbins; nx++) {
      Int_t binnumber = k*nxbins+nx;	  
      nbgroup[k] += ((AliTRDEntriesInfo  *)entriesd)->At(binnumber);
    }
  }
  max = -10.0;
  groupmax = -1;
  for(Int_t io = 0; io < numberofgroup; io++){
    if(nbgroup[io] > max){
      max = nbgroup[io];
      groupmax = io;
    }
  }
  if((max == 0.0) || (groupmax < 0.0) || (groupmax >= numberofgroup)) return kFALSE;

  return kTRUE;

}
//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPHTGraphErrors(Int_t det, Int_t group , const Char_t * name)
{
  //
  // Convert the fVectorPHMean, fVectorPHSquares and fVectorPHEntries in TGraphErrors
  //

  // Take the info
  fPHEntries[det] = ((AliTRDEntriesInfo  *)GetPHEntries(det,kTRUE));
  fPHMean[det]    = ((AliTRDPhInfo *)GetPHMean(det,kTRUE));
  fPHSquares[det] = ((AliTRDPhInfo *)GetPHSquares(det,kTRUE));
  

  // Axis
  Float_t sf = 10.0;
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam, take the default 10MHz");
  }
  sf = parCom->GetSamplingFrequency();
  // Axis
  Double_t x[35];  // Xaxis
  Double_t y[35];  // Sum/entries
  Double_t ex[35]; // Nentries
  Double_t ey[35]; // Sum of square/nentries
  Double_t step = 0.0;
  Double_t min  = 0.0;
  if(sf > 0.0) step = 1.0 / sf;
  min  = 0.0;
  Int_t offset = group*fTimeMax;
  
  // Fill histo
  for (Int_t k = 0; k < fTimeMax; k++) {
    x[k]  = min + k*step;
    y[k]  = 0.0;
    ex[k] = 0.0;
    ey[k] = 0.0;   
    Int_t bin = offset+k;
    // Fill only if there is more than 0 something
    if (((AliTRDEntriesInfo *)fPHEntries[det])->At(bin) > 0) {
      ex[k] = ((AliTRDEntriesInfo *)fPHEntries[det])->At(bin);
      y[k]  = ((AliTRDPhInfo *)fPHMean[det])->At(bin);
      ey[k] = ((AliTRDPhInfo *)fPHSquares[det])->AtS(bin);
    }
  }

  // Define the TGraphErrors
  if(!fGraph) fGraph = new TGraphErrors(fTimeMax,&x[0],&y[0],&ex[0],&ey[0]);
  else{ 
    fGraph->~TGraphErrors();
    new(fGraph) TGraphErrors(fTimeMax,&x[0],&y[0],&ex[0],&ey[0]);
  } 
  fGraph->SetTitle(name); 

  return fGraph;

}
//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPRFTGraphErrors(Int_t det, Int_t group , const Char_t * name)
{
  //
  // Convert the fVectorPRFMean, fVectorPRFSquares and fVectorPRFEntries in TGraphErrors
  //

  // Take the info
  fPRFEntries[det] = ((AliTRDEntriesInfo *)GetPRFEntries(det,kTRUE));
  fPRFMean[det]    = ((AliTRDPrfInfo     *)GetPRFMean(det,kTRUE));
  fPRFSquares[det] = ((AliTRDPrfInfo     *)GetPRFSquares(det,kTRUE));
  

  // Axis
  Double_t x[200];  // Xaxis
  Double_t y[200];  // Sum/entries
  Double_t ex[200]; //Nentries
  Double_t ey[200]; // Sum of square/nentries
  Double_t step = 0.0;
  Double_t min  = 0.0;
  if(fNumberBinPRF) step = (2*TMath::Abs(fPRFRange)) / fNumberBinPRF;
  min  = -TMath::Abs(fPRFRange) + step / 2.0;
  Int_t offset = group*fNumberBinPRF;
  //printf("number of total: %d\n",fNumberBinPRF);
  // Fill histo
  for (Int_t k = 0; k < fNumberBinPRF; k++) {
    x[k]  = min + k*step;
    y[k]  = 0.0;
    ex[k] = 0.0;
    ey[k] = 0.0;
    Int_t bin = offset+k;
    // Fill only if there is more than 0 something
    if (((AliTRDEntriesInfo *)fPRFEntries[det])->At(bin) > 0) {
      ex[k] = ((AliTRDEntriesInfo *)fPRFEntries[det])->At(bin);
      y[k]  = ((AliTRDPrfInfo *)fPRFMean[det])->At(bin);
      ey[k] = ((AliTRDPrfInfo *)fPRFSquares[det])->At(bin);
    }
    //printf("Number of entries %f for %d\n",ex[k],k);
  }

  // Define the TGraphErrors
  if(!fGraph) fGraph = new TGraphErrors(fNumberBinPRF,&x[0],&y[0],&ex[0],&ey[0]);
  else{ 
    fGraph->~TGraphErrors();
    new(fGraph) TGraphErrors(fNumberBinPRF,&x[0],&y[0],&ex[0],&ey[0]);
  }
  fGraph->SetTitle(name); 

  return fGraph;



}
//_____________________________________________________________________________
TH1F *AliTRDCalibraVector::CorrectTheError(const TGraphErrors *hist, Int_t &nbEntries)
{
  //
  // In the case of the vectors method the trees contains TGraphErrors for PH and PRF
  // to be able to add them after
  // We convert it to a TH1F to be able to applied the same fit function method
  // After having called this function you can not add the statistics anymore
  //

  Int_t nbins       = hist->GetN();
  Double_t *x       = hist->GetX();
  Double_t *entries = hist->GetEX();
  Double_t *mean    = hist->GetY();
  Double_t *square  = hist->GetEY();
  nbEntries   = 0;

  if (nbins < 2) {
    return 0x0; 
  }

  Double_t step     = x[1] - x[0]; 
  Double_t minvalue = x[0] - step/2;
  Double_t maxvalue = x[(nbins-1)] + step/2;

  if(!fHisto) fHisto = new TH1F("projcorrecterror","",nbins,minvalue,maxvalue);
  else{ 
    fHisto->~TH1F();
    new(fHisto) TH1F("projcorrecterror","",nbins,minvalue,maxvalue);
  }

  for (Int_t k = 0; k < nbins; k++) {
    fHisto->SetBinContent(k+1,mean[k]);
    if (entries[k] > 0.0) {
      nbEntries += (Int_t) entries[k];
      Double_t d = TMath::Abs(square[k] - (mean[k]*mean[k]));
      fHisto->SetBinError(k+1,TMath::Sqrt(d/entries[k]));
    }
    else {
      fHisto->SetBinError(k+1,0.0);
    }
  }

  return fHisto;
 
}  
//_____________________________________________________________________________
TH1F *AliTRDCalibraVector::ConvertVectorCHHisto(Int_t det, Int_t group, const Char_t * name)
{
  //
  // Convert the fVectorCHEntries in TH1F
  //

  // Take the info
  fCHEntries[det] = ((AliTRDEntriesInfo *)GetCHEntries(det,kTRUE));
  
  // Init the stuff
  if(!fHisto) fHisto = new TH1F(name,name,fNumberBinCharge,0,300);
  else{ 
    fHisto->~TH1F();
    new(fHisto) TH1F(name,name,fNumberBinCharge,0,300);
  }
  fHisto->Sumw2();
  Int_t offset = group*fNumberBinCharge;
  // Fill histo
  for (Int_t k = 0; k < fNumberBinCharge; k++) {
    Int_t bin = offset+k;
    fHisto->SetBinContent(k+1,((AliTRDEntriesInfo *)fCHEntries[det])->At(bin));
    fHisto->SetBinError(k+1,TMath::Sqrt(TMath::Abs(((AliTRDEntriesInfo *)fCHEntries[det])->At(bin))));
  }
  
  return fHisto;

} 
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    AliTRDEntriesInfo**arr = &fPHEntries[0];
    return (TObject *) GetEntriesPH(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPRFEntries(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    AliTRDEntriesInfo**arr = &fPRFEntries[0];
    return (TObject *) GetEntriesPRF(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetCHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    AliTRDEntriesInfo**arr = &fCHEntries[0];
    return (TObject *) GetEntriesCH(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPHMean(Int_t det
                                           , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new array
    //
    AliTRDPhInfo**arr = &fPHMean[0];
    return (TObject *) GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPHSquares(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new array
    //
    AliTRDPhInfo**arr = &fPHSquares[0];
    return (TObject *)  GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPRFMean(Int_t det
                                            , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new array
    //
    AliTRDPrfInfo**arr = &fPRFMean[0];
    return (TObject *) GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
TObject* AliTRDCalibraVector::GetPRFSquares(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to ROC Calibration
    // if force is true create a new array
    //
    AliTRDPrfInfo**arr = &fPRFSquares[0];
    return (TObject *) GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
AliTRDEntriesInfo* AliTRDCalibraVector::GetEntriesCH(Int_t det
                                               , AliTRDEntriesInfo** arr
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to UShort_t array Entries
    // if force is true create a new UShort_t array if it doesn't exist allready
    //
  if ( (!force) || (arr[det]))
    return (AliTRDEntriesInfo*)arr[det];

  // if we are forced and TArrayI doesn't yes exist create it
  Int_t ngroup = GetTotalNumberOfBinsInDetector(det,0,fNumberBinCharge); 
  // init
  arr[det] = new AliTRDEntriesInfo(ngroup);
  
  return (AliTRDEntriesInfo*)arr[det];

}
//_____________________________________________________________________
AliTRDEntriesInfo* AliTRDCalibraVector::GetEntriesPRF(Int_t det
                                               , AliTRDEntriesInfo** arr
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to UShort_t array Entries
    // if force is true create a new UShort_t array if it doesn't exist allready
    //
  if ( (!force) || (arr[det]))
    return (AliTRDEntriesInfo*)arr[det];

  // if we are forced and TArrayI doesn't yes exist create it
  Int_t ngroup = GetTotalNumberOfBinsInDetector(det,2,fNumberBinPRF); 
  // init
  arr[det] = new AliTRDEntriesInfo(ngroup);
  
  return (AliTRDEntriesInfo*)arr[det];

}
//_____________________________________________________________________
AliTRDEntriesInfo *AliTRDCalibraVector::GetEntriesPH(Int_t det
                                              , AliTRDEntriesInfo ** arr
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to UShort_t array Entries
    // if force is true create a new UShort_t array if it doesn't exist allready
    //
    if ( (!force) || (arr[det]))
	return (AliTRDEntriesInfo *)arr[det];

    // if we are forced and UShort_t doesn't yet exist create it
    Int_t ngroup = GetTotalNumberOfBinsInDetector(det,1,fTimeMax); 
    // init
    arr[det] = new AliTRDEntriesInfo(ngroup);
    
    return (AliTRDEntriesInfo*)arr[det];
   
}
//_____________________________________________________________________
AliTRDPhInfo* AliTRDCalibraVector::GetMeanSquaresPH(Int_t det
                                                  , AliTRDPhInfo** arr
                                                  , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Float_t array Mean or Squares
    // if force is true create a new Float_t array if it doesn't exist allready
    //
    if ( (!force) || (arr[det]))
	return (AliTRDPhInfo*)arr[det];

    // if we are forced and Float_t array doesn't yes exist create it
    Int_t ngroup = GetTotalNumberOfBinsInDetector(det,1,fTimeMax); 
    // init
    arr[det] = new AliTRDPhInfo(ngroup);
    
    return ((AliTRDPhInfo *)arr[det]);
}
//_____________________________________________________________________
AliTRDPrfInfo* AliTRDCalibraVector::GetMeanSquaresPRF(Int_t det
                                                   , AliTRDPrfInfo** arr
                                                   , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Float_t array Mean or Squares
    // if force is true create a new array if it doesn't exist allready
    //
  if ( (!force) || (arr[det]))
    return arr[det];
  
  // if we are forced and the array doesn't yet exist create it
  Int_t ngroup = GetTotalNumberOfBinsInDetector(det,2,fNumberBinPRF); 
  // init
  arr[det] = new AliTRDPrfInfo(ngroup);
  
  return (AliTRDPrfInfo*)arr[det];
   
}
//_____________________________________________________________________________
Int_t AliTRDCalibraVector::GetTotalNumberOfBinsInDetector(Int_t det, Int_t i, Int_t nbBin) const
{
  //
  // Get the total number of bins (Nb of bins*Nb of groups) in the detector det for the group i
  //

  Int_t ngroup = 0;
  Int_t stack  = AliTRDgeometry::GetStack(det);
  if(stack == 2) ngroup = fDetCha2[i]*nbBin;
  else ngroup = fDetCha0[i]*nbBin;

  return ngroup;
 
}
//____________________________________________________________________________
Int_t AliTRDCalibraVector::GetNz(Int_t i) const
{
  //
  // Get Nz the granularity in row
  //

  Int_t nz = 0;
  if(i==0) nz = (Int_t)(fModeCH>>4);
  if(i==1) nz = (Int_t)(fModePH>>4);
  if(i==2) nz = (Int_t)(fModePRF>>4);
  
  return nz;

}
//____________________________________________________________________________
Int_t AliTRDCalibraVector::GetNrphi(Int_t i) const
{
  //
  // Get Nrphi the granularity in col
  //

  Int_t nrphi = 0;
  if(i==0) nrphi = (Int_t)(fModeCH&15);
  if(i==1) nrphi = (Int_t)(fModePH&15);
  if(i==2) nrphi = (Int_t)(fModePRF&15);
  
  return nrphi;

}
//_________________________________________________________________________________
TString AliTRDCalibraVector::GetNamePH() const
{
  //
  // Get the name of PH to know the granularity
  //
  
  Int_t nz = GetNz(1);
  Int_t nrphi = GetNrphi(1);

  TString name("Nz");
  name += nz;
  name += "Nrphi";
  name += nrphi;
  
  return name;

}   
//_________________________________________________________________________________
TString AliTRDCalibraVector::GetNameCH() const
{
  //
  // Get the name of CH to know the granularity
  //
  
  Int_t nz = GetNz(0);
  Int_t nrphi = GetNrphi(0);

  TString name("Nz");
  name += nz;
  name += "Nrphi";
  name += nrphi;
  
  return name;

}   
//_________________________________________________________________________________
TString AliTRDCalibraVector::GetNamePRF() const
{
  //
  // Get the name of PRF to know the granularity
  //

  Int_t nz = GetNz(2);
  Int_t nrphi = GetNrphi(2);
  
  TString name("Nz");
  name += nz;
  name += "Nrphi";
  name += nrphi;
  name += "Ngp";
  name += fNbGroupPRF;
  
  return name;

}
//____________________________________________________________________________
void AliTRDCalibraVector::SetNzNrphi(Int_t i, Int_t nz, Int_t nrphi) 
{
  //
  // Set NzNrphi for the granularity
  //
  
  if(i==0) {
    fModeCH = nz;
    fModeCH = fModeCH << 4;
    fModeCH |= nrphi;
  }
  if(i==1) {
    fModePH = nz;
    fModePH = fModePH << 4;
    fModePH |= nrphi;
  }
  if(i==2) {
    fModePRF = nz;
    fModePRF = fModePRF << 4;
    fModePRF |= nrphi;
  }
  
}

