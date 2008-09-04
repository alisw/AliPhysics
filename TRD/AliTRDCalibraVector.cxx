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
#include <TMath.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TFile.h>

#include "AliLog.h"

#include "AliTRDCalibraVector.h"
#include "AliTRDCommonParam.h"
#include "AliTRDarrayF.h"
#include "AliTRDarrayI.h"

ClassImp(AliTRDCalibraVector)

//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector()
  :TObject()
  ,fPHEntries(0x0)
  ,fPHMean(0x0)
  ,fPHSquares(0x0)
  ,fPRFEntries(0x0)
  ,fPRFMean(0x0)
  ,fPRFSquares(0x0)
  ,fCHEntries(0x0)
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fVectorPHMean(540)
  ,fVectorPHSquares(540)
  ,fVectorPHEntries(540)
  ,fVectorCHEntries(540)
  ,fVectorPRFMean(540)
  ,fVectorPRFSquares(540)
  ,fVectorPRFEntries(540)
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fTimeMax(0)
  ,fPRFRange(1.5)
{
  //
  // Default constructor
  //
  
  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = 0;
    fDetCha2[k] = 0;
  }
 
}
//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector(const AliTRDCalibraVector &c)
  :TObject(c)
  ,fPHEntries(0x0)
  ,fPHMean(0x0)
  ,fPHSquares(0x0)
  ,fPRFEntries(0x0)
  ,fPRFMean(0x0)
  ,fPRFSquares(0x0)
  ,fCHEntries(0x0)
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fVectorPHMean(540)
  ,fVectorPHSquares(540)
  ,fVectorPHEntries(540)
  ,fVectorCHEntries(540)
  ,fVectorPRFMean(540)
  ,fVectorPRFSquares(540)
  ,fVectorPRFEntries(540)
  ,fNumberBinCharge(c.fNumberBinCharge)
  ,fNumberBinPRF(c.fNumberBinPRF)
  ,fTimeMax(c.fTimeMax)
  ,fPRFRange(c.fPRFRange)
{
  //
  // Copy constructor
  //
  for (Int_t idet = 0; idet < 540; idet++){
    
    const AliTRDarrayI *phEntries  = (AliTRDarrayI*)c.fVectorPHEntries.UncheckedAt(idet);
    const AliTRDarrayF *phMean     = (AliTRDarrayF*)c.fVectorPHMean.UncheckedAt(idet);
    const AliTRDarrayF *phSquares  = (AliTRDarrayF*)c.fVectorPHSquares.UncheckedAt(idet);

    const AliTRDarrayI *prfEntries  = (AliTRDarrayI*)c.fVectorPRFEntries.UncheckedAt(idet);
    const AliTRDarrayF *prfMean     = (AliTRDarrayF*)c.fVectorPRFMean.UncheckedAt(idet);
    const AliTRDarrayF *prfSquares  = (AliTRDarrayF*)c.fVectorPRFSquares.UncheckedAt(idet);

    const AliTRDarrayI *chEntries  = (AliTRDarrayI*)c.fVectorCHEntries.UncheckedAt(idet);

    if ( phEntries != 0x0 )  fVectorPHEntries.AddAt(new AliTRDarrayI(*phEntries), idet);
    if ( phMean != 0x0 )     fVectorPHMean.AddAt(new AliTRDarrayF(*phMean), idet);
    if ( phSquares != 0x0 )  fVectorPHSquares.AddAt(new AliTRDarrayF(*phSquares), idet);

    if ( prfEntries != 0x0 )  fVectorPRFEntries.AddAt(new AliTRDarrayI(*prfEntries), idet);
    if ( prfMean != 0x0 )     fVectorPRFMean.AddAt(new AliTRDarrayF(*prfMean), idet);
    if ( prfSquares != 0x0 )  fVectorPRFSquares.AddAt(new AliTRDarrayF(*prfSquares), idet);


    if ( chEntries != 0x0 )   fVectorCHEntries.AddAt(new AliTRDarrayI(*chEntries), idet);

  }

  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = c.fDetCha0[k];
    fDetCha2[k] = c.fDetCha2[k];
  }

  fVectorPHEntries.SetName(c.fVectorPHEntries.GetName());
  fVectorCHEntries.SetName(c.fVectorCHEntries.GetName());
  fVectorPRFEntries.SetName(c.fVectorPRFEntries.GetName());
  
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
  fVectorPHMean.Delete();
  fVectorPHSquares.Delete();
  fVectorPHEntries.Delete();
  fVectorCHEntries.Delete();
  fVectorPRFMean.Delete();
  fVectorPRFSquares.Delete();
  fVectorPRFEntries.Delete();

  if(fPHEntries) delete fPHEntries;
  if(fPHMean) delete fPHMean;
  if(fPHSquares) delete fPHSquares;
  if(fPRFEntries) delete fPRFEntries;
  if(fPRFMean) delete fPRFMean;
  if(fPRFSquares) delete fPRFSquares;
  if(fCHEntries) delete fCHEntries;

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
    fCHEntries = ((AliTRDarrayI *)GetCHEntries(det,kTRUE));
  }

  Int_t entries  = fCHEntries->At(group*fNumberBinCharge+bin);
  
  Int_t entriesn = entries+1;
  fCHEntries->AddAt(entriesn,group*fNumberBinCharge+bin);
  
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
    fPRFEntries = ((AliTRDarrayI *)GetPRFEntries(det,kTRUE));
    fPRFMean    = ((AliTRDarrayF *)GetPRFMean(det,kTRUE));
    fPRFSquares = ((AliTRDarrayF *)GetPRFSquares(det,kTRUE));
  }

  Int_t entries  = fPRFEntries->At(group*fNumberBinPRF+bin);
  Float_t mean   = fPRFMean->At(group*fNumberBinPRF+bin);
  Float_t square = fPRFSquares->At(group*fNumberBinPRF+bin);
  
  Int_t entriesn = entries+1;
  fPRFEntries->AddAt(entriesn,group*fNumberBinPRF+bin);
  Float_t meann = (mean*((Float_t)entries)+y)/((Float_t)entriesn);
  fPRFMean->AddAt(meann,group*fNumberBinPRF+bin);
  Float_t squaren = ((square*((Float_t)entries))+(y*y))/((Float_t)entriesn);
  fPRFSquares->AddAt(squaren,group*fNumberBinPRF+bin);

  
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
    fPHEntries = ((AliTRDarrayI *)GetPHEntries(det,kTRUE));
    fPHMean    = ((AliTRDarrayF *)GetPHMean(det,kTRUE));
    fPHSquares = ((AliTRDarrayF *)GetPHSquares(det,kTRUE));
  }

  Int_t entries  = fPHEntries->At(group*fTimeMax+bin);
  Float_t mean   = fPHMean->At(group*fTimeMax+bin);
  Float_t square = fPHSquares->At(group*fTimeMax+bin);
  
  Int_t entriesn = entries+1;
  fPHEntries->AddAt(entriesn,group*fTimeMax+bin);
  Float_t meann = (mean*((Float_t)entries)+value)/((Float_t)entriesn);
  fPHMean->AddAt(meann,group*fTimeMax+bin);
  Float_t squaren = ((square*((Float_t)entries))+(value*value))/((Float_t)entriesn);
  fPHSquares->AddAt(squaren,group*fTimeMax+bin);

  
  fDetectorPH = det;

  return kTRUE;
  
}
//__________________________________________________________________________________
Bool_t AliTRDCalibraVector::Add(AliTRDCalibraVector *calvect)
{
  //
  // Add a other AliTRCalibraVector to this one
  //

  // Check compatibility
  if(fNumberBinCharge != calvect->GetNumberBinCharge()) return kFALSE;
  if(fNumberBinPRF    != calvect->GetNumberBinPRF()) return kFALSE;
  if(fPRFRange        != calvect->GetPRFRange()) return kFALSE;
  if(fTimeMax         != calvect->GetTimeMax()) return kFALSE;
  for(Int_t k = 0; k < 3; k++){
    if(fDetCha0[k] != calvect->GetDetCha0(k)) return kFALSE;
    if(fDetCha2[k] != calvect->GetDetCha2(k)) return kFALSE;
  }

  //printf("pass0!\n");

  // Add
  for (Int_t idet = 0; idet < 540; idet++){
    
    const AliTRDarrayI *phEntriesvect   = calvect->GetPHEntries(idet);
    const AliTRDarrayF *phMeanvect      = calvect->GetPHMean(idet);
    const AliTRDarrayF *phSquaresvect   = calvect->GetPHSquares(idet);
    
    const AliTRDarrayI *prfEntriesvect  = calvect->GetPRFEntries(idet);
    const AliTRDarrayF *prfMeanvect     = calvect->GetPRFMean(idet);
    const AliTRDarrayF *prfSquaresvect  = calvect->GetPRFSquares(idet);
    
    const AliTRDarrayI *chEntriesvect   = calvect->GetCHEntries(idet);

    //printf("idet %d!\n",idet);

    //printf("phEntriesvect %d\n",(Bool_t) phEntriesvect);
    //printf("phMeanvect %d\n",(Bool_t) phMeanvect);
    //printf("phSquaresvect %d\n",(Bool_t) phSquaresvect);

    //printf("prfEntriesvect %d\n",(Bool_t) prfEntriesvect);
    //printf("prfMeanvect %d\n",(Bool_t) prfMeanvect);
    //printf("prfSquaresvect %d\n",(Bool_t) prfSquaresvect);

    //printf("chEntriesvect %d\n",(Bool_t) chEntriesvect);

    if ( phEntriesvect != 0x0 ){
      //Take the stuff
      fPHEntries = ((AliTRDarrayI *)GetPHEntries(idet,kTRUE));
      fPHMean    = ((AliTRDarrayF *)GetPHMean(idet,kTRUE));
      fPHSquares = ((AliTRDarrayF *)GetPHSquares(idet,kTRUE));
      Int_t total = fPHEntries->GetSize();
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = fPHEntries->At(k);
	Float_t mean   = fPHMean->At(k);
	Float_t square = fPHSquares->At(k);
  
	Int_t entriesn = entries+phEntriesvect->At(k);
	if(entriesn <= 0) continue;
	fPHEntries->AddAt(entriesn,k);
	Float_t meann = (mean*((Float_t)entries)+phMeanvect->At(k)*((Float_t)phEntriesvect->At(k)))/((Float_t)entriesn);
	fPHMean->AddAt(meann,k);
	Float_t sq      = phSquaresvect->At(k)*((Float_t)phEntriesvect->At(k));
	Float_t squaren = ((square*((Float_t)entries))+sq)/((Float_t)entriesn);
	fPHSquares->AddAt(squaren,k);
	//printf("test ph!\n");
      }
    }     

    if ( prfEntriesvect != 0x0 ){
      //Take the stuff
      fPRFEntries = ((AliTRDarrayI *)GetPRFEntries(idet,kTRUE));
      fPRFMean    = ((AliTRDarrayF *)GetPRFMean(idet,kTRUE));
      fPRFSquares = ((AliTRDarrayF *)GetPRFSquares(idet,kTRUE));
      Int_t total = fPRFEntries->GetSize();
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = fPRFEntries->At(k);
	Float_t mean   = fPRFMean->At(k);
	Float_t square = fPRFSquares->At(k);

	//printf("entries0 %d\n",entries);
	//printf("mean0 %f\n",mean);
	//printf("square0 %f\n",square);

	//printf("entries1 %d\n",prfEntriesvect->At(k));
	//printf("mean1 %f\n",prfMeanvect->At(k));
	//printf("square1 %f\n",prfSquaresvect->At(k));
  

	//printf("entries0 size %d\n",fPRFEntries->GetSize());
	//printf("mean0 size %d\n",fPRFMean->GetSize());
	//printf("squares0 size %d\n",fPRFSquares->GetSize());

	//printf("entries1 size %d\n",prfEntriesvect->GetSize());
	//printf("mean1 size %d\n",prfMeanvect->GetSize());
	//printf("squares1 size %d\n",prfSquaresvect->GetSize());


	Int_t entriesn = entries+prfEntriesvect->At(k);
	if(entriesn <= 0) continue;
	fPRFEntries->AddAt(entriesn,k);
	Float_t meann = (mean*((Float_t)entries)+prfMeanvect->At(k)*((Float_t)prfEntriesvect->At(k)))/((Float_t)entriesn);
	fPRFMean->AddAt(meann,k);
	Float_t sq      = prfSquaresvect->At(k)*((Float_t)prfEntriesvect->At(k));
	Float_t squaren = ((square*((Float_t)entries))+sq)/((Float_t)entriesn);
	fPRFSquares->AddAt(squaren,k);
	//printf("test prf!\n");
      }
    }

    if ( chEntriesvect != 0x0 ){
      //Take the stuff
      fCHEntries = ((AliTRDarrayI *)GetCHEntries(idet,kTRUE));
      Int_t total = fCHEntries->GetSize();
      //if(idet == 180) printf("total %d\n",total);
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = fCHEntries->At(k);
	Int_t entriesn = entries+chEntriesvect->At(k);
	//if((idet == 180) && ((entries != 0) || (entriesn != 0))) printf("for k %d we have entries %d and entriesn %d\n",k,entries,entriesn);
	if(entriesn <= 0) continue;
	fCHEntries->AddAt(entriesn,k);
      }
      //printf("test ch!\n");
    }           
  }
  
  return kTRUE;
} 
//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPHTGraphErrors(Int_t det, Int_t group
                                        , const Char_t * name)
{
  //
  // Convert the fVectorPHMean, fVectorPHSquares and fVectorPHEntries in TGraphErrors
  //

  // Take the info
  fPHEntries = ((AliTRDarrayI *)GetPHEntries(det,kTRUE));
  fPHMean    = ((AliTRDarrayF *)GetPHMean(det,kTRUE));
  fPHSquares = ((AliTRDarrayF *)GetPHSquares(det,kTRUE));
  

  // Init the stuff
  TGraphErrors *histo;
  Float_t sf = 10.0;
  AliTRDCommonParam *parCom = AliTRDCommonParam::Instance();
  if (!parCom) {
    AliInfo("Could not get CommonParam, take the default 10MHz");
  }
  sf = parCom->GetSamplingFrequency();
  // Axis
  Double_t *x;
  Double_t *y;
  Double_t *ex;
  Double_t *ey;
  Double_t step = 0.0;
  Double_t min  = 0.0;
  x  = new Double_t[fTimeMax]; // Xaxis
  y  = new Double_t[fTimeMax]; // Sum/entries
  ex = new Double_t[fTimeMax]; // Nentries
  ey = new Double_t[fTimeMax]; // Sum of square/nentries
  step = 1.0 / sf;
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
    if (fPHEntries->At(bin) > 0) {
      ex[k] = fPHEntries->At(bin);
      y[k]  = fPHMean->At(bin);
      ey[k] = fPHSquares->At(bin);
    }
  }

  // Define the TGraphErrors
  histo = new TGraphErrors(fTimeMax,x,y,ex,ey);
  histo->SetTitle(name); 
  return histo;

}
//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPRFTGraphErrors(Int_t det, Int_t group
                                        , const Char_t * name)
{
  //
  // Convert the fVectorPRFMean, fVectorPRFSquares and fVectorPRFEntries in TGraphErrors
  //

  // Take the info
  fPRFEntries = ((AliTRDarrayI *)GetPRFEntries(det,kTRUE));
  fPRFMean    = ((AliTRDarrayF *)GetPRFMean(det,kTRUE));
  fPRFSquares = ((AliTRDarrayF *)GetPRFSquares(det,kTRUE));
  

  // Init the stuff
  TGraphErrors *histo;
  // Axis
  Double_t *x;
  Double_t *y;
  Double_t *ex;
  Double_t *ey;
  Double_t step = 0.0;
  Double_t min  = 0.0;
  x  = new Double_t[fNumberBinPRF]; // Xaxis
  y  = new Double_t[fNumberBinPRF]; // Sum/entries
  ex = new Double_t[fNumberBinPRF]; // Nentries
  ey = new Double_t[fNumberBinPRF]; // Sum of square/nentries
  step = (2*TMath::Abs(fPRFRange)) / fNumberBinPRF;
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
    if (fPRFEntries->At(bin) > 0) {
      ex[k] = fPRFEntries->At(bin);
      y[k]  = fPRFMean->At(bin);
      ey[k] = fPRFSquares->At(bin);
    }
    //printf("Number of entries %f for %d\n",ex[k],k);
  }

  // Define the TGraphErrors
  histo = new TGraphErrors(fNumberBinPRF,x,y,ex,ey);
  histo->SetTitle(name); 
  return histo;

}  
//_____________________________________________________________________________
TH1F *AliTRDCalibraVector::ConvertVectorCHHisto(Int_t det, Int_t group
                                              , const Char_t * name)
{
  //
  // Convert the fVectorCHEntries in TH1F
  //

  // Take the info
  fCHEntries = ((AliTRDarrayI *)GetCHEntries(det,kTRUE));
  
  // Init the stuff
  TH1F *histo = new TH1F(name,name,fNumberBinCharge,0,300);
  histo->Sumw2();
  Int_t offset = group*fNumberBinCharge;
  // Fill histo
  for (Int_t k = 0; k < fNumberBinCharge; k++) {
    Int_t bin = offset+k;
    histo->SetBinContent(k+1,fCHEntries->At(bin));
    histo->SetBinError(k+1,TMath::Sqrt(TMath::Abs(fCHEntries->At(bin))));
  }
  
  return histo;

} 
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetPHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPHEntries;
    return GetEntriesPH(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetPRFEntries(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPRFEntries;
    return GetEntriesPRF(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetCHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorCHEntries;
    return GetEntriesCH(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetPHMean(Int_t det
                                           , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPHMean;
    return GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetPHSquares(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPHSquares;
    return GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetPRFMean(Int_t det
                                            , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPRFMean;
    return GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetPRFSquares(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TObjArray *arr = &fVectorPRFSquares;
    return GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetEntriesCH(Int_t det
                                              , TObjArray* arr
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to AliTRDarrayI Entries
    // if force is true create a new AliTRDarrayI if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDarrayI*)arr->UncheckedAt(det);

    // if we are forced and AliTRDarrayI doesn't yes exist create it
    AliTRDarrayI *croc = new AliTRDarrayI();
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[0]*fNumberBinCharge;
    else ngroup = fDetCha0[0]*fNumberBinCharge;
    // init
    croc->Expand(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      croc->AddAt(0,k);
    }
    arr->AddAt(croc,det);
    return croc;
}
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetEntriesPRF(Int_t det
                                               , TObjArray* arr
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to AliTRDarrayI Entries
    // if force is true create a new AliTRDarrayI if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDarrayI*)arr->UncheckedAt(det);

    // if we are forced and AliTRDarrayI doesn't yes exist create it
    AliTRDarrayI *croc = new AliTRDarrayI();
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[2]*fNumberBinPRF;
    else ngroup = fDetCha0[2]*fNumberBinPRF;
    // init
    croc->Expand(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      croc->AddAt(0,k);
    }
    arr->AddAt(croc,det);
    return croc;

}
//_____________________________________________________________________
AliTRDarrayI* AliTRDCalibraVector::GetEntriesPH(Int_t det
                                              , TObjArray* arr
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to AliTRDarrayI Entries
    // if force is true create a new AliTRDarrayI if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDarrayI*)arr->UncheckedAt(det);

    // if we are forced and AliTRDarrayI doesn't yes exist create it
    AliTRDarrayI *croc = new AliTRDarrayI();
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[1]*fTimeMax;
    else ngroup = fDetCha0[1]*fTimeMax;
    // init
    croc->Expand(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      croc->AddAt(0,k);
    }
    arr->AddAt(croc,det);
    return croc;

}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetMeanSquaresPH(Int_t det
                                                  , TObjArray* arr
                                                  , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to AliTRDarrayF Mean or Squares
    // if force is true create a new AliTRDarrayF if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDarrayF*)arr->UncheckedAt(det);

    // if we are forced and AliTRDarrayF doesn't yes exist create it
    AliTRDarrayF *croc = new AliTRDarrayF();
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[1]*fTimeMax;
    else ngroup = fDetCha0[1]*fTimeMax;
    // init
    croc->Expand(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      croc->AddAt(0.0,k);
    }
    arr->AddAt(croc,det);
    return croc;
}
//_____________________________________________________________________
AliTRDarrayF* AliTRDCalibraVector::GetMeanSquaresPRF(Int_t det
                                                   , TObjArray* arr
                                                   , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to AliTRDarrayF Mean or Squares
    // if force is true create a new AliTRDarrayF if it doesn't exist allready
    //
    if ( !force || arr->UncheckedAt(det) )
	return (AliTRDarrayF*)arr->UncheckedAt(det);

    // if we are forced and AliTRDarrayF doesn't yes exist create it
    AliTRDarrayF *croc = new AliTRDarrayF();
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[2]*fNumberBinPRF;
    else ngroup = fDetCha0[2]*fNumberBinPRF;
    // init
    croc->Expand(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      croc->AddAt(0.0,k);
    }
    arr->AddAt(croc,det);
    return croc;
}
//_____________________________________________________________________________
Int_t AliTRDCalibraVector::GetStack(Int_t d) const
{
  //
  // Reconstruct the stack number from the detector number
  //

  return ((Int_t) (d % 30) / 6);

}

