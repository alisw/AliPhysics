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
#include "TArrayF.h"
#include "TArrayI.h"

ClassImp(AliTRDCalibraVector)

//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector()
  :TObject()
  ,fNameCH("CH2d")
  ,fNamePH("PH2d")
  ,fNamePRF("PRF2d")
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fTimeMax(0)
  ,fPRFRange(1.5)
{
  //
  // Default constructor
  //

  for (Int_t idet = 0; idet < 540; idet++){
    
    fPHEntries[idet]=new TArrayI();
    fPHMean[idet]=new TArrayF();
    fPHSquares[idet]=new TArrayF();

    fPRFEntries[idet]=new TArrayI();
    fPRFMean[idet]=new TArrayF();
    fPRFSquares[idet]=new TArrayF();


    fCHEntries[idet]=new TArrayI();
    
  }
  
  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = 0;
    fDetCha2[k] = 0;
  }
 
}
//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector(const AliTRDCalibraVector &c)
  :TObject(c)
  ,fNameCH("CH2d")
  ,fNamePH("PH2d")
  ,fNamePRF("PRF2d")
  ,fDetectorPH(-1)
  ,fDetectorCH(-1)
  ,fDetectorPRF(-1)
  ,fNumberBinCharge(c.fNumberBinCharge)
  ,fNumberBinPRF(c.fNumberBinPRF)
  ,fTimeMax(c.fTimeMax)
  ,fPRFRange(c.fPRFRange)
{
  //
  // Copy constructor
  //
  for (Int_t idet = 0; idet < 540; idet++){
    
    const TArrayI *phEntries  = (TArrayI*)c.fPHEntries[idet];
    const TArrayF *phMean     = (TArrayF*)c.fPHMean[idet];
    const TArrayF *phSquares  = (TArrayF*)c.fPHSquares[idet];

    const TArrayI *prfEntries  = (TArrayI*)c.fPRFEntries[idet];
    const TArrayF *prfMean     = (TArrayF*)c.fPRFMean[idet];
    const TArrayF *prfSquares  = (TArrayF*)c.fPRFSquares[idet];

    const TArrayI *chEntries  = (TArrayI*)c.fCHEntries[idet];

    if ( phEntries != 0x0 )  fPHEntries[idet]=new TArrayI(*phEntries);
    if ( phMean != 0x0 )     fPHMean[idet]=new TArrayF(*phMean);
    if ( phSquares != 0x0 )  fPHSquares[idet]=new TArrayF(*phSquares);

    if ( prfEntries != 0x0 )  fPRFEntries[idet]=new TArrayI(*prfEntries);
    if ( prfMean != 0x0 )     fPRFMean[idet]=new TArrayF(*prfMean);
    if ( prfSquares != 0x0 )  fPRFSquares[idet]=new TArrayF(*prfSquares);


    if ( chEntries != 0x0 )   fCHEntries[idet]=new TArrayI(*chEntries);

  }

  for(Int_t k = 0; k < 3; k++){
    fDetCha0[k] = c.fDetCha0[k];
    fDetCha2[k] = c.fDetCha2[k];
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
    fCHEntries[det] = ((TArrayI *)GetCHEntries(det,kTRUE));
  }

  Int_t entries  = ((TArrayI *)fCHEntries[det])->At(group*fNumberBinCharge+bin);
  
  Int_t entriesn = entries+1;
  ((TArrayI *)fCHEntries[det])->AddAt(entriesn,group*fNumberBinCharge+bin);
  
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
    fPRFEntries[det] = ((TArrayI *)GetPRFEntries(det,kTRUE));
    fPRFMean[det]    = ((TArrayF *)GetPRFMean(det,kTRUE));
    fPRFSquares[det] = ((TArrayF *)GetPRFSquares(det,kTRUE));
  }

  Int_t entries  = ((TArrayI *)fPRFEntries[det])->At((Int_t)group*fNumberBinPRF+bin);
  Float_t mean   = ((TArrayI *)fPRFMean[det])->At((Int_t)group*fNumberBinPRF+bin);
  Float_t square = ((TArrayI *)fPRFSquares[det])->At((Int_t)group*fNumberBinPRF+bin);
  
  Int_t entriesn = entries+1;
  ((TArrayI *)fPRFEntries[det])->AddAt(entriesn,(Int_t)group*fNumberBinPRF+bin);
  Float_t meann = (mean*((Float_t)entries)+y)/((Float_t)entriesn);
  ((TArrayF *)fPRFMean[det])->AddAt(meann,(Int_t)group*fNumberBinPRF+bin);
  Float_t squaren = ((square*((Float_t)entries))+(y*y))/((Float_t)entriesn);
  ((TArrayF *)fPRFSquares[det])->AddAt(squaren,(Int_t)group*fNumberBinPRF+bin);

  
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
    fPHEntries[det] = ((TArrayI *)GetPHEntries(det,kTRUE));
    fPHMean[det]    = ((TArrayF *)GetPHMean(det,kTRUE));
    fPHSquares[det] = ((TArrayF *)GetPHSquares(det,kTRUE));
  }

  Int_t entries  = ((TArrayI *)fPHEntries[det])->At(group*fTimeMax+bin);
  Float_t mean   = ((TArrayI *)fPHMean[det])->At(group*fTimeMax+bin);
  Float_t square = ((TArrayI *)fPHSquares[det])->At(group*fTimeMax+bin);
  
  Int_t entriesn = entries+1;
  ((TArrayI *)fPHEntries[det])->AddAt(entriesn,group*fTimeMax+bin);
  Float_t meann = (mean*((Float_t)entries)+value)/((Float_t)entriesn);
  ((TArrayF *)fPHMean[det])->AddAt(meann,group*fTimeMax+bin);
  Float_t squaren = ((square*((Float_t)entries))+(value*value))/((Float_t)entriesn);
  ((TArrayF *)fPHSquares[det])->AddAt(squaren,group*fTimeMax+bin);

  
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
    
    const TArrayI *phEntriesvect   = calvect->GetPHEntries(idet);
    const TArrayF *phMeanvect      = calvect->GetPHMean(idet);
    const TArrayF *phSquaresvect   = calvect->GetPHSquares(idet);
    
    const TArrayI *prfEntriesvect  = calvect->GetPRFEntries(idet);
    const TArrayF *prfMeanvect     = calvect->GetPRFMean(idet);
    const TArrayF *prfSquaresvect  = calvect->GetPRFSquares(idet);
    
    const TArrayI *chEntriesvect   = calvect->GetCHEntries(idet);

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
      fPHEntries[idet] = ((TArrayI *)GetPHEntries(idet,kTRUE));
      fPHMean[idet]    = ((TArrayF *)GetPHMean(idet,kTRUE));
      fPHSquares[idet] = ((TArrayF *)GetPHSquares(idet,kTRUE));
      Int_t total = ((TArrayI *)fPHEntries[idet])->GetSize();
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = ((TArrayI *)fPHEntries[idet])->At(k);
	Float_t mean   = ((TArrayF *)fPHMean[idet])->At(k);
	Float_t square = ((TArrayF *)fPHSquares[idet])->At(k);
  
	Int_t entriesn = entries+((TArrayI *)phEntriesvect)->At(k);
	if(entriesn <= 0) continue;
	((TArrayI *)fPHEntries[idet])->AddAt(entriesn,k);
	Float_t meann = (mean*((Float_t)entries)+((TArrayF *)phMeanvect)->At(k)*((Float_t)((TArrayI *)phEntriesvect)->At(k)))/((Float_t)entriesn);
	((TArrayF *)fPHMean[idet])->AddAt(meann,k);
	Float_t sq      = ((TArrayF *)phSquaresvect)->At(k)*((Float_t)((TArrayI *)phEntriesvect)->At(k));
	Float_t squaren = ((square*((Float_t)entries))+sq)/((Float_t)entriesn);
	((TArrayF *)fPHSquares[idet])->AddAt(squaren,k);
	//printf("test ph!\n");
      }
    }     

    if ( prfEntriesvect != 0x0 ){
      //Take the stuff
      fPRFEntries[idet] = ((TArrayI *)GetPRFEntries(idet,kTRUE));
      fPRFMean[idet]    = ((TArrayF *)GetPRFMean(idet,kTRUE));
      fPRFSquares[idet] = ((TArrayF *)GetPRFSquares(idet,kTRUE));
      Int_t total = fPRFEntries[idet]->GetSize();
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = ((TArrayI *)fPRFEntries[idet])->At(k);
	Float_t mean   = ((TArrayF *)fPRFMean[idet])->At(k);
	Float_t square = ((TArrayF *)fPRFSquares[idet])->At(k);

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


	Int_t entriesn = entries+((TArrayI *)prfEntriesvect)->At(k);
	if(entriesn <= 0) continue;
	((TArrayI *)fPRFEntries[idet])->AddAt(entriesn,k);
	Float_t meann = (mean*((Float_t)entries)+((TArrayF *)prfMeanvect)->At(k)*((Float_t)((TArrayI *)prfEntriesvect)->At(k)))/((Float_t)entriesn);
	((TArrayF *)fPRFMean[idet])->AddAt(meann,k);
	Float_t sq      = ((TArrayF *)prfSquaresvect)->At(k)*((Float_t)((TArrayI *)prfEntriesvect)->At(k));
	Float_t squaren = ((square*((Float_t)entries))+sq)/((Float_t)entriesn);
	((TArrayF *)fPRFSquares[idet])->AddAt(squaren,k);
	//printf("test prf!\n");
      }
    }

    if ( chEntriesvect != 0x0 ){
      //Take the stuff
      fCHEntries[idet] = ((TArrayI *)GetCHEntries(idet,kTRUE));
      Int_t total = fCHEntries[idet]->GetSize();
      //if(idet == 180) printf("total %d\n",total);
      // Add
      for(Int_t k = 0; k < total; k++){
	Int_t entries  = ((TArrayI *)fCHEntries[idet])->At(k);
	Int_t entriesn = entries+((TArrayI *)chEntriesvect)->At(k);
	//if((idet == 180) && ((entries != 0) || (entriesn != 0))) printf("for k %d we have entries %d and entriesn %d\n",k,entries,entriesn);
	if(entriesn <= 0) continue;
	((TArrayI *)fCHEntries[idet])->AddAt(entriesn,k);
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
  fPHEntries[det] = ((TArrayI *)GetPHEntries(det,kTRUE));
  fPHMean[det]    = ((TArrayF *)GetPHMean(det,kTRUE));
  fPHSquares[det] = ((TArrayF *)GetPHSquares(det,kTRUE));
  

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
    if (fPHEntries[det]->At(bin) > 0) {
      ex[k] = ((TArrayI *)fPHEntries[det])->At(bin);
      y[k]  = ((TArrayF *)fPHMean[det])->At(bin);
      ey[k] = ((TArrayF *)fPHSquares[det])->At(bin);
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
  fPRFEntries[det] = ((TArrayI *)GetPRFEntries(det,kTRUE));
  fPRFMean[det]    = ((TArrayF *)GetPRFMean(det,kTRUE));
  fPRFSquares[det] = ((TArrayF *)GetPRFSquares(det,kTRUE));
  

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
    if (fPRFEntries[det]->At(bin) > 0) {
      ex[k] = ((TArrayF *)fPRFEntries[det])->At(bin);
      y[k]  = ((TArrayF *)fPRFMean[det])->At(bin);
      ey[k] = ((TArrayF *)fPRFSquares[det])->At(bin);
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
  fCHEntries[det] = ((TArrayI *)GetCHEntries(det,kTRUE));
  
  // Init the stuff
  TH1F *histo = new TH1F(name,name,fNumberBinCharge,0,300);
  histo->Sumw2();
  Int_t offset = group*fNumberBinCharge;
  // Fill histo
  for (Int_t k = 0; k < fNumberBinCharge; k++) {
    Int_t bin = offset+k;
    histo->SetBinContent(k+1,((TArrayI *)fCHEntries[det])->At(bin));
    histo->SetBinError(k+1,TMath::Sqrt(TMath::Abs(((TArrayI *)fCHEntries[det])->At(bin))));
  }
  
  return histo;

} 
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetPHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayI**arr = &fPHEntries[0];
    return GetEntriesPH(det, arr, force);
}
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetPRFEntries(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayI**arr = &fPRFEntries[0];
    return GetEntriesPRF(det, arr, force);
}
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetCHEntries(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayI**arr = &fCHEntries[0];
    return GetEntriesCH(det, arr, force);
}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetPHMean(Int_t det
                                           , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayF**arr = &fPHMean[0];
    return GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetPHSquares(Int_t det
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayF**arr = &fPHSquares[0];
    return GetMeanSquaresPH(det, arr, force);
}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetPRFMean(Int_t det
                                            , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayF**arr = &fPRFMean[0];
    return GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetPRFSquares(Int_t det
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to Carge ROC Calibration
    // if force is true create a new histogram if it doesn't exist allready
    //
    TArrayF**arr = &fPRFSquares[0];
    return GetMeanSquaresPRF(det, arr, force);
}
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetEntriesCH(Int_t det
                                              , TArrayI** arr
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to TArrayI Entries
    // if force is true create a new TArrayI if it doesn't exist allready
    //
  if ( !force || (((TArrayI *)arr[det])->GetSize()>0))
	return (TArrayI*)arr[det];

    // if we are forced and TArrayI doesn't yes exist create it
  Int_t stack  = GetStack(det);
  Int_t ngroup = 0;
  if(stack == 2) ngroup = fDetCha2[0]*fNumberBinCharge;
  else ngroup = fDetCha0[0]*fNumberBinCharge;
  // init
  ((TArrayI *)arr[det])->Set(ngroup);
  for(Int_t k = 0; k < ngroup; k++){
    ((TArrayI *)arr[det])->AddAt(0,k);
  }
  return (TArrayI*)arr[det];
}
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetEntriesPRF(Int_t det
                                               , TArrayI** arr
                                               , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to TArrayI Entries
    // if force is true create a new TArrayI if it doesn't exist allready
    //
  if ( !force || (((TArrayI *)arr[det])->GetSize()>0))
	return (TArrayI*)arr[det];

  // if we are forced and TArrayI doesn't yes exist create it
  Int_t stack  = GetStack(det);
  Int_t ngroup = 0;
  if(stack == 2) ngroup = fDetCha2[2]*fNumberBinPRF;
  else ngroup = fDetCha0[2]*fNumberBinPRF;
  // init
  ((TArrayI *)arr[det])->Set(ngroup);
  for(Int_t k = 0; k < ngroup; k++){
    ((TArrayI *)arr[det])->AddAt(0,k);
  }
  return (TArrayI*)arr[det];

}
//_____________________________________________________________________
TArrayI* AliTRDCalibraVector::GetEntriesPH(Int_t det
                                              , TArrayI** arr
                                              , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to TArrayI Entries
    // if force is true create a new TArrayI if it doesn't exist allready
    //
    if ( !force || (((TArrayI *)arr[det])->GetSize()>0))
	return (TArrayI*)arr[det];

    // if we are forced and TArrayI doesn't yes exist create it
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[1]*fTimeMax;
    else ngroup = fDetCha0[1]*fTimeMax;
    // init
    ((TArrayI *)arr[det])->Set(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      ((TArrayI *)arr[det])->AddAt(0,k);
    }
    return (TArrayI*)arr[det];

}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetMeanSquaresPH(Int_t det
                                                  , TArrayF** arr
                                                  , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to TArrayF Mean or Squares
    // if force is true create a new TArrayF if it doesn't exist allready
    //
    if ( !force || (((TArrayF *)arr[det])->GetSize()>0))
	return (TArrayF*)arr[det];

    // if we are forced and TArrayF doesn't yes exist create it
    Int_t stack  = GetStack(det);
    Int_t ngroup = 0;
    if(stack == 2) ngroup = fDetCha2[1]*fTimeMax;
    else ngroup = fDetCha0[1]*fTimeMax;
    // init
    ((TArrayF *)arr[det])->Set(ngroup);
    for(Int_t k = 0; k < ngroup; k++){
      ((TArrayF *)arr[det])->AddAt(0.0,k);
    }
    return ((TArrayF *)arr[det]);
}
//_____________________________________________________________________
TArrayF* AliTRDCalibraVector::GetMeanSquaresPRF(Int_t det
                                                   , TArrayF** arr
                                                   , Bool_t force) /*FOLD00*/
{
    //
    // return pointer to TArrayF Mean or Squares
    // if force is true create a new TArrayF if it doesn't exist allready
    //
  if ( !force || (((TArrayF *)arr[det])->GetSize()>0))
    return arr[det];
  
  // if we are forced and TArrayF doesn't yes exist create it
  Int_t stack  = GetStack(det);
  Int_t ngroup = 0;
  if(stack == 2) ngroup = fDetCha2[2]*fNumberBinPRF;
  else ngroup = fDetCha0[2]*fNumberBinPRF;
  // init
  ((TArrayF *)arr[det])->Set(ngroup);
  for(Int_t k = 0; k < ngroup; k++){
    ((TArrayF *)arr[det])->AddAt(0.0,k);
  }
  return ((TArrayF *)arr[det]);
}
//_____________________________________________________________________________
Int_t AliTRDCalibraVector::GetStack(Int_t d) const
{
  //
  // Reconstruct the stack number from the detector number
  //

  return ((Int_t) (d % 30) / 6);

}

