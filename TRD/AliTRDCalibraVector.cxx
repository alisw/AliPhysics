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

/////////////////////////////////////////////////////////////////////////////////
//                                                                             
// AliTRDCalibraVector                                                             
//                                                                             
// This class is for the vector methode of the TRD calibration.        
//                            
// Author:
//   R. Bailhache (R.Bailhache@gsi.de)
//                            
//////////////////////////////////////////////////////////////////////////////////////

#include <TGraphErrors.h>
#include <TH1F.h>
#include <TTree.h>
#include <TObjArray.h>
#include <TMath.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TROOT.h>
#include <TFile.h>

#include "AliLog.h"

#include "AliTRDCalibraVector.h"
#include "AliTRDCommonParam.h"

ClassImp(AliTRDCalibraVector)

//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector()
  :TObject()
  ,fVectorPH(new TObjArray())
  ,fPlaPH(new TObjArray())
  ,fVectorCH(new TObjArray())
  ,fPlaCH(new TObjArray())
  ,fVectorPRF(new TObjArray())
  ,fPlaPRF(new TObjArray())
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fTimeMax(0)
{
  //
  // Default constructor
  //
   
}

//______________________________________________________________________________________
AliTRDCalibraVector::AliTRDCalibraVector(const AliTRDCalibraVector &c)
  :TObject(c)
  ,fVectorPH(new TObjArray())
  ,fPlaPH(new TObjArray())
  ,fVectorCH(new TObjArray())
  ,fPlaCH(new TObjArray())
  ,fVectorPRF(new TObjArray())
  ,fPlaPRF(new TObjArray())
  ,fNumberBinCharge(0)
  ,fNumberBinPRF(0)
  ,fTimeMax(0)
{
  //
  // Copy constructor
  //

}

//____________________________________________________________________________________
AliTRDCalibraVector::~AliTRDCalibraVector()
{
  //
  // AliTRDCalibraVector destructor
  //

}

//_____________________________________________________________________________
Int_t AliTRDCalibraVector::SearchBin(Float_t value, Int_t i) const
{
  //
  // Search the bin
  //

  Int_t reponse      = 0;
  Int_t fbinmin      = 0;
  Int_t fbinmax      = (Int_t) value;
  Int_t fNumberOfBin = -1;

  // Charge
  if (i == 0) {
    fbinmax      = 300;
    fbinmin      = 0;
    fNumberOfBin = fNumberBinCharge;
  }

  // PRF
  if (i == 2) {
    fbinmax      =   1;
    fbinmin      =  -1;
    fNumberOfBin = fNumberBinPRF;
  }

  // Return -1 if out
  if ((value >= fbinmax) || 
      (value <  fbinmin)) {
    return -1;
  }
  // Sinon
  else {
    reponse = (Int_t) ((fNumberOfBin*(value-fbinmin)) / (fbinmax-fbinmin));
  }

  return reponse;

}

//_____________________________________________________________________________
Int_t AliTRDCalibraVector::SearchInVector(Int_t group, Int_t i) const
{
  //
  // Search if the calibration group "group" has already been
  // initialised by a previous track in the vector
  //

  if (i == 0) {
    for (Int_t k = 0; k < (Int_t) fPlaCH->GetEntriesFast(); k++) {
      if (((AliTRDPlace *) fPlaCH->At(k))->GetPlace() == group) {
        return k;
      }
    }
    return -1;
  }

  if (i == 1) {
    for (Int_t k = 0; k < (Int_t) fPlaPH->GetEntriesFast(); k++) {
      if (((AliTRDPlace *) fPlaPH->At(k))->GetPlace() == group) {
        return k;
      }
    }
    return -1;
  }

  if (i == 2) {
    for (Int_t k = 0; k < (Int_t) fPlaPRF->GetEntriesFast(); k++) {
      if (((AliTRDPlace *) fPlaPRF->At(k))->GetPlace() == group) {
        return k;
      }
    }
    return -1;
  }

  return -1;

}

//_____________________________________________________________________________
Int_t AliTRDCalibraVector::SearchInTreeVector(TObjArray *vectorplace, Int_t group) const
{
  //
  // Search if the calibration group "group" is present in the tree
  //

  for (Int_t k = 0; k < (Int_t) vectorplace->GetEntriesFast(); k++) {
    if (((AliTRDPlace *) vectorplace->At(k))->GetPlace() == group) {
      return k;
    }
  }

  return -1;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorCH(Int_t group, Float_t value)
{
  //
  // Fill the vector if a new calibration group "group" or update the
  // values of the calibration group "group" if already here  
  //

  // Search bin
  Int_t bin = SearchBin(value,0);
  // Out
  if ((bin < 0) || (bin >= fNumberBinCharge)) {
    return kFALSE; 
  }

  // Search place
  Int_t place = SearchInVector(group,0);

  // New group
  if (place == -1) {
    AliTRDPlace *placegroup = new AliTRDPlace();
    placegroup->SetPlace(group);
    fPlaCH->Add((TObject *) placegroup);
    // Variable
    AliTRDCTInfo *fCHInfo = new AliTRDCTInfo();
    UShort_t *entries = new UShort_t[fNumberBinCharge];
    // Initialise first
    for(Int_t k = 0; k < fNumberBinCharge; k++) {
      entries[k] = 0;
    }
    // Add the value
    entries[bin]= 1;
    // Set
    fCHInfo->SetEntries(entries);
    // Set in the vector
    fVectorCH->Add((TObject *) fCHInfo);
  }
  // Group already exits
  else {
    // Variable
    AliTRDCTInfo *fCHInfo = new AliTRDCTInfo();
    // Retrieve
    fCHInfo = ((AliTRDCTInfo *) fVectorCH->At(place));
    UShort_t *entries = fCHInfo->GetEntries();
    // Add
    entries[bin]++;
    // Set
    fCHInfo->SetEntries(entries);
    // Update the vector
    fVectorCH->AddAt((TObject *) fCHInfo,place);
  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorPRF(Int_t group, Float_t x, Float_t y)
{
  //
  // Fill the vector if a new calibration group "group" or update the
  // values of the calibration group "group" if already here  
  //

  // Search bin
  Int_t bin = SearchBin(x,2);
  // Out
  if ((bin < 0) || (bin >= fNumberBinPRF)) {
    return kFALSE; 
  }

  // Search place
  Int_t place = SearchInVector(group,2);

  // New group
  if (place == -1) {

    AliTRDPlace *placegroup = new AliTRDPlace();
    placegroup->SetPlace(group);
    fPlaPRF->Add((TObject *) placegroup);
    AliTRDPInfo *fPRFInfo = new AliTRDPInfo();

    Float_t  *sum       = new Float_t[fNumberBinPRF];
    Float_t  *sumsquare = new Float_t[fNumberBinPRF];
    UShort_t *entries   = new UShort_t[fNumberBinPRF];

    // Initialise first
    for (Int_t k = 0; k < fNumberBinPRF; k++) {
      sum[k]       = 0.0;
      sumsquare[k] = 0.0;
      entries[k]   = 0;
    }

    // Add the value
    sum[bin]       += y;
    sumsquare[bin] += y*y;
    entries[bin]++;

    // Set
    fPRFInfo->SetSum(sum);
    fPRFInfo->SetSumSquare(sumsquare);
    fPRFInfo->SetEntries(entries);

    // Set in the vector
    fVectorPRF->Add((TObject *) fPRFInfo);
        
  }
  // Group already exits
  else {

    AliTRDPInfo *fPRFInfo = new AliTRDPInfo();
    // Retrieve
    fPRFInfo = (AliTRDPInfo *) fVectorPRF->At(place);

    Float_t  *sum       = fPRFInfo->GetSum();
    Float_t  *sumsquare = fPRFInfo->GetSumSquare();
    UShort_t *entries   = fPRFInfo->GetEntries();

    // Add
    Double_t calcul       = (((Double_t) fPRFInfo->GetEntries()[bin])
                           * ((Double_t) fPRFInfo->GetSum()[bin]) + (Double_t) y)
                          / (((Double_t) fPRFInfo->GetEntries()[bin]) + 1);
    sum[bin]       = (Float_t) calcul;
    Double_t calculsquare = (((Double_t) fPRFInfo->GetSumSquare()[bin])
                           * ((Double_t) fPRFInfo->GetEntries()[bin]) + ((Double_t) y)*((Double_t) y))
                          / (((Double_t) fPRFInfo->GetEntries()[bin]) + 1);
    sumsquare[bin] = (Float_t) calculsquare;
    entries[bin]++;

    // Set
    fPRFInfo->SetSum(sum);
    fPRFInfo->SetSumSquare(sumsquare);
    fPRFInfo->SetEntries(entries);
 
    // Update the vector
    fVectorPRF->AddAt((TObject *) fPRFInfo,place);

  }

  return kTRUE;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::UpdateVectorPH(Int_t group, Int_t time, Float_t value)
{
  //
  // Fill the vector if a new calibration group "group" or update
  // the values of the calibration group "group" if already here  
  //

  // Search bin
  Int_t bin = time;
  // Out
  if ((bin <         0) || 
      (bin >= fTimeMax)) {
    return kFALSE; 
  }

  // Search place
  Int_t place = SearchInVector(group,1);

  // New group
  if(place == -1){

    AliTRDPlace *placegroup = new AliTRDPlace();
    placegroup->SetPlace(group);
    fPlaPH->Add((TObject *) placegroup);
    AliTRDPInfo *fPHInfo = new AliTRDPInfo();

    Float_t  *sum       = new Float_t[fTimeMax];
    Float_t  *sumsquare = new Float_t[fTimeMax];
    UShort_t *entries   = new UShort_t[fTimeMax];

    // Initialise first
    for (Int_t k = 0; k < fTimeMax; k++) {
      sum[k]       = 0.0;
      sumsquare[k] = 0.0;
      entries[k]   = 0;
    }

    // Add the value
    sum[bin]       += value;
    sumsquare[bin] += value*value;
    entries[bin]++;

    // Set
    fPHInfo->SetSum(sum);
    fPHInfo->SetSumSquare(sumsquare);
    fPHInfo->SetEntries(entries);

    // Set in the vector
    fVectorPH->Add((TObject *) fPHInfo);

  }
  // Group already exits
  else {

    AliTRDPInfo *fPHInfo = new AliTRDPInfo();
    // Retrieve
    fPHInfo = (AliTRDPInfo *) fVectorPH->At(place);

    Float_t  *sum       = fPHInfo->GetSum();
    Float_t  *sumsquare = fPHInfo->GetSumSquare();
    UShort_t *entries   = fPHInfo->GetEntries();

    // Add
    Double_t calcul       = (((Double_t) fPHInfo->GetEntries()[bin])
                           * ((Double_t) fPHInfo->GetSum()[bin]) + (Double_t) value)
                          / (((Double_t) fPHInfo->GetEntries()[bin]) + 1);
    sum[bin]       = (Float_t) calcul;
    Double_t calculsquare = ((((Double_t) fPHInfo->GetSumSquare()[bin])
                            * ((Double_t) fPHInfo->GetEntries()[bin])) 
                          + (((Double_t) value) * ((Double_t)value))) 
                          / (((Double_t) fPHInfo->GetEntries()[bin]) + 1);
    sumsquare[bin] = (Float_t) calculsquare;
    entries[bin]++;

    // Set
    fPHInfo->SetSum(sum);
    fPHInfo->SetSumSquare(sumsquare);
    fPHInfo->SetEntries(entries);

    // Update the vector
    fVectorPH->AddAt((TObject *) fPHInfo,place);

  }

  return kTRUE;

}  

//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPHisto(Int_t place
                                        , const Char_t * name) const
{
  //
  // Convert the PInfo in a TGraphErrors
  //

  AliTRDPInfo *pInfo = new AliTRDPInfo();
  // Retrieve
  pInfo = ((AliTRDPInfo *) fVectorPH->At(place));
  TGraphErrors *histo;
  histo = ConvertVectorPHistoI((AliTRDPInfo *)pInfo, name);
  
  return histo;

} 

//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::ConvertVectorPHistoI(AliTRDPInfo *pInfo
                                               , const Char_t *name) const
{
  //
  // Convert the PInfo in a 1D grapherror, name must contains "PRF"
  // if PRF calibration and not "PRF" for Vdrift calibration
  //
 
  TGraphErrors *histo;
  const Char_t *pattern1 = "PRF";

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

  // Ntimes
  Int_t ntimes = 0;
  if (strstr(name,pattern1)) {
    ntimes = fNumberBinPRF;
  }
  else {
    ntimes = fTimeMax;
  }
  x  = new Double_t[ntimes]; // Xaxis
  y  = new Double_t[ntimes]; // Mean
  ex = new Double_t[ntimes]; // Nentries
  ey = new Double_t[ntimes]; // Sum of square/nentries

  // Init histo
  if (!strstr(name,pattern1)) {
    step = 1.0 / sf;
    min  = 0.0;
  }
  else {
    step = (1.0 - (-1.0)) / fNumberBinPRF;
    min  = -1.0 + step / 2.0;
  }

  // Fill histo
  for (Int_t k = 0; k < ntimes; k++) {
    x[k]  = min + k*step;
    y[k]  = 0.0;
    ex[k] = 0.0;
    ey[k] = 0.0;
    // Fill only if there is more than 0 something
    if (pInfo->GetEntries()[k] > 0) {
      ex[k] = pInfo->GetEntries()[k];
      y[k]  = pInfo->GetSum()[k];
      ey[k] =  pInfo->GetSumSquare()[k];
    }
  }

  // Define the TGraphErrors
  histo = new TGraphErrors(ntimes,x,y,ex,ey);
  histo->SetTitle(name); 
  return histo;

}

//_____________________________________________________________________________
TH1F *AliTRDCalibraVector::ConvertVectorCTHisto(Int_t place
                                              , const Char_t * name) const
{
  //
  // Convert the CTInfo in a 1D histo
  //

  AliTRDCTInfo *cHInfo = new AliTRDCTInfo();
  // Retrieve
  cHInfo = ((AliTRDCTInfo *) fVectorCH->At(place));
  TH1F *histo;
  histo = ConvertVectorCTHistoI((AliTRDCTInfo *)cHInfo,(const Char_t *) name);
  
  return histo;

} 

//_____________________________________________________________________________
TH1F *AliTRDCalibraVector::ConvertVectorCTHistoI(AliTRDCTInfo *cTInfo
                                               , const Char_t * name) const
{
  //
  // Convert the CTInfo in a 1D histo
  //

  TH1F *histo;
  
  Int_t     ntimes  = fNumberBinCharge;
  UShort_t *entries = cTInfo->GetEntries();
  
  // Init histo
  histo = new TH1F(name,name,fNumberBinCharge,0,300);
  histo->Sumw2();
  // Fill histo
  for (Int_t k = 0; k < ntimes; k++) {
    histo->SetBinContent(k+1,entries[k]);
    histo->SetBinError(k+1,TMath::Sqrt(TMath::Abs(entries[k])));
  }
  
  return histo;

}

//_____________________________________________________________________________
TTree *AliTRDCalibraVector::ConvertVectorCTTreeHisto(TObjArray *vVectorCT
                                                   , TObjArray *pPlaCT
                                                   , const Char_t *name
                                                   , const Char_t *nametitle) const
{
  //
  // Convert the vector in a tree with two branchs: the group number
  // and the TH1F histo reconstructed from the vector
  //

  // Size of the things
  Int_t ntotal = (Int_t) pPlaCT->GetEntriesFast();
  if (ntotal == 0) {
    AliInfo("nothing to write!");
    TTree *treeCT = new TTree(name,nametitle);
    return treeCT;
  }
  
  // Variable of the tree
  Int_t groupnumber = -1; // Group calibration
  TH1F      *histo = 0x0;
  TObjArray  vectorCT = *vVectorCT;
  TObjArray  plaCT    = *pPlaCT;

  // Init the tree
  TTree *treeCT = new TTree(name,nametitle);
  treeCT->Branch("groupnumber",&groupnumber,"groupnumber/I");
  treeCT->Branch("histo","TH1F",&histo,32000,0);

  // Fill
  Int_t k = 0;
  while (k < ntotal) {
    TString nome(name);
    groupnumber  = ((AliTRDPlace *) plaCT.At(0))->GetPlace();
    nome        += groupnumber;
    histo        = ConvertVectorCTHistoI(((AliTRDCTInfo *) vectorCT.At(0)),nome);
    treeCT->Fill();
    vectorCT.RemoveAt(0);
    vectorCT.Compress();
    plaCT.RemoveAt(0);
    plaCT.Compress();
    k++;
  } 

  return treeCT;

}

//_____________________________________________________________________________
TTree *AliTRDCalibraVector::ConvertVectorPTreeHisto(TObjArray *vVectorP
                                                  , TObjArray *pPlaP
                                                  , const Char_t *name
                                                  , const Char_t *nametitle) const
{
  //
  // Convert the vector in a tree with two branchs: the group number
  // and the TGraphErrors histo reconstructed from the vector.
  // The name must contain "PRF" for PRF calibration and not "PRF"
  // for Vdrift calibration
  //

  // Size of the things
  Int_t ntotal = (Int_t) pPlaP->GetEntriesFast();
  if (ntotal == 0) {
    AliInfo("nothing to write!");
    TTree *treeP = new TTree(name,nametitle);
    return treeP;
  }

  // Variable of the tree
  Int_t groupnumber = -1; // Group calibration
  TGraphErrors *histo   = 0x0;
  TObjArray     vectorP = *vVectorP;
  TObjArray     plaP    = *pPlaP;

  // Init the tree
  TTree *treeP = new TTree(name,nametitle);
  treeP->Branch("groupnumber",&groupnumber,"groupnumber/I");
  treeP->Branch("histo","TGraphErrors",&histo,32000,0);

  // Fill
  Int_t k = 0;
  while (k < ntotal) {
    TString nome(name);
    groupnumber = ((AliTRDPlace *) plaP.At(0))->GetPlace();
    nome       += groupnumber;
    histo       = ConvertVectorPHistoI((AliTRDPInfo *) vectorP.At(0),nome);
    treeP->Fill();
    vectorP.RemoveAt(0);
    vectorP.Compress();
    plaP.RemoveAt(0);
    plaP.Compress();
    k++;
  } 

  return treeP;

}

//_____________________________________________________________________________
TObjArray *AliTRDCalibraVector::ConvertTreeVector(TTree *tree) const
{
  //
  // Convert the branch groupnumber of the tree taken from
  // TRD.calibration.root in case of vector method in a std::vector 
  // to be faster
  //

  // Initialise
  TObjArray *vectorplace = new TObjArray();
  
  // Variable of the tree
  Int_t groupnumber = -1; // Group calibration

  // Set the branch
  tree->SetBranchAddress("groupnumber",&groupnumber);
    
  // Fill
  Int_t ntotal = tree->GetEntries();
  for (Int_t k = 0; k < ntotal; k++) {
    tree->GetEntry(k);
    AliTRDPlace *placegroupnumber = new AliTRDPlace();
    placegroupnumber->SetPlace(groupnumber);
    vectorplace->Add((TObject *) placegroupnumber);
  }
  
  return vectorplace;

}

//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::MergeVectorCT(TObjArray *vVectorCT2, TObjArray *pPlaCT2)
{
  //
  // Add the two vectors and place the result in the first
  //

  if (((Int_t) pPlaCT2->GetEntriesFast()) != ((Int_t) vVectorCT2->GetEntriesFast())){
    AliInfo("VectorCT2 doesn't correspond to PlaCT2!");
    return kFALSE;
  }
 
  // CH case
  for (Int_t k = 0; k < (Int_t) fPlaCH->GetEntriesFast(); k++) {
    
    // Look if PlaCT1[k] it is also in the second vector
    Int_t place = -1;
    for (Int_t j = 0; j < (Int_t) pPlaCT2->GetEntriesFast(); j++) {
      if (((AliTRDPlace *) pPlaCT2->At(j))->GetPlace() == 
            ((AliTRDPlace *) fPlaCH->At(k))->GetPlace()) {
	place = j;
	break;
      }
    }
    
    // If not in the second vector nothing to do

    // If in the second vector
    if (place != -1) {
      
      AliTRDCTInfo *fCTInfo = new AliTRDCTInfo();
      UShort_t *entries = new UShort_t[fNumberBinCharge];
      
      for (Int_t nu = 0; nu < fNumberBinCharge; nu++) {
	entries[nu] = ((AliTRDCTInfo *)  fVectorCH->At(((AliTRDPlace *) fPlaCH->At(k))->GetPlace()))->GetEntries()[nu]
                    + ((AliTRDCTInfo *) vVectorCT2->At(((AliTRDPlace *) fPlaCH->At(k))->GetPlace()))->GetEntries()[nu];
      }
      
      // Set
      fCTInfo->SetEntries(entries);

      // Nothing to do on PlaCT1
      
      // Update the vector 
      fVectorCH->AddAt((TObject *) fCTInfo,((AliTRDPlace *) fPlaCH->At(k))->GetPlace());

    }
    
  } 
 
  // And at the end the vector in CT2 but not in CH1
  for (Int_t k = 0; k < (Int_t) pPlaCT2->GetEntriesFast(); k++) {
    
    // Look if pPlaCT2[k] it is also in the second vector
    Int_t place = -1;
    for (Int_t j = 0; j < (Int_t) fPlaCH->GetEntriesFast(); j++) {
      if (((AliTRDPlace *) fPlaCH->At(j))->GetPlace() == ((AliTRDPlace *) pPlaCT2->At(k))->GetPlace()) {
	place = j;
	break;
      }
    }

    // If not in the first vector
    if (place == -1) {
      
      AliTRDCTInfo *fCTInfo = new AliTRDCTInfo();     
      fCTInfo = ((AliTRDCTInfo *) vVectorCT2->At(((AliTRDPlace *) pPlaCT2->At(k))->GetPlace()));
      
      // Add at the end 
      fPlaCH->Add((TObject *) (pPlaCT2->At(k)));
      fVectorCH->Add((TObject *) fCTInfo);

    }
    
  }
  
  return kTRUE;
  
}

//_____________________________________________________________________________
Bool_t AliTRDCalibraVector::MergeVectorP(TObjArray *vVectorP2
                                       , TObjArray *pPlaP2
                                       , Int_t i)
{
  //
  // Add the two vectors and place the result in the first
  //

  if (((Int_t) pPlaP2->GetEntriesFast()) != ((Int_t) vVectorP2->GetEntriesFast())) {
    AliInfo("VectorP2 doesn't correspond to PlaP2!");
    return kFALSE;
  }

  // PH case
  if (i == 1) {

     for (Int_t k = 0; k < (Int_t) fPlaPH->GetEntriesFast(); k++) {
       
       // Look if fPlaPH[k] it is also in the second vector
       Int_t place = -1;
       for (Int_t j = 0; j < (Int_t) pPlaP2->GetEntriesFast(); j++) {
	 if (((AliTRDPlace *) pPlaP2->At(j))->GetPlace() == ((AliTRDPlace *) fPlaPH->At(k))->GetPlace()) {
	   place = j;
	   break;
	 }
       }
       
       // If not in the second vector nothing to do

       // If in the second vector
       if (place != -1) {

	 AliTRDPInfo *fPInfo = new AliTRDPInfo();
	 UShort_t *entries   = new UShort_t[fTimeMax];
	 Float_t  *sum       = new Float_t[fTimeMax];
	 Float_t  *sumsquare = new Float_t[fTimeMax];

	 for (Int_t nu = 0; nu < fTimeMax; nu++) {
	   
	   entries[nu]   = ((AliTRDPInfo *) fVectorPH->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu]
                         + ((AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu];
	   
	   Double_t calcul       = ((((Double_t) ((AliTRDPInfo *) fVectorPH->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetSum()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) fVectorPH->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu]))
                                  + (((Double_t) ((AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetSum()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu])))
                                 / ((Double_t) fPInfo->GetEntries()[nu]);

	   sum[nu]       = (Float_t) calcul;
	   
	   Double_t calculsquare = ((((Double_t) ((AliTRDPInfo *) fVectorPH->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetSumSquare()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) fVectorPH->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu]))
                                  + (((Double_t) ((AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetSumSquare()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) fPlaPH->At(k))->GetPlace()))->GetEntries()[nu])))
                                 / ((Double_t) fPInfo->GetEntries()[nu]);

	   sumsquare[nu] = calculsquare;

	 }

	 // Set
	 fPInfo->SetSum(sum);
	 fPInfo->SetSumSquare(sumsquare);
	 fPInfo->SetEntries(entries);
	 
	 // Nothing to do on PlaCT1
	 
	 // Update the vector VectorCT1
	 fVectorPH->AddAt((TObject *) fPInfo,((AliTRDPlace *) fPlaPH->At(k))->GetPlace());
	 
       }

     }

     // And at the end the vector in P2 but not in CH1
     for (Int_t k = 0; k < (Int_t) pPlaP2->GetEntriesFast(); k++) {
       
       // Look if PlaCT2[k] it is also in the second vector
       Int_t place = -1;
       for (Int_t j = 0; j < (Int_t) fPlaPH->GetEntriesFast(); j++) {
	 if (((AliTRDPlace *) fPlaPH->At(j))->GetPlace() == ((AliTRDPlace *) pPlaP2->At(k))->GetPlace()) {
	   place = j;
	   break;
	 }
       }
       
       // If not in the first vector
       if (place == -1) {
	 	 
	 AliTRDPInfo *fPInfo = new AliTRDPInfo();
	 fPInfo = (AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) pPlaP2->At(k))->GetPlace());
	 
	 // Add at the end of CH1
	 fPlaPH->Add(((TObject *) pPlaP2->At(k)));
	 fVectorPH->Add((TObject *) fPInfo);

       }

     }

   }

   // PRF case
   if (i == 1) {

     for (Int_t k = 0; k < (Int_t) fPlaPRF->GetEntriesFast(); k++) {

       // Look if fPlaPRF[k] it is also in the second vector
       Int_t place = -1;
       for (Int_t j = 0; j < (Int_t) pPlaP2->GetEntriesFast(); j++) {
	 if (((AliTRDPlace *) pPlaP2->At(j))->GetPlace() == ((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()) {
	   place = j;
	   break;
	 }
       }

       // If not in the second vector nothing to do

       // If in the second vector
       if (place != -1) {
	
	 AliTRDPInfo *fPInfo = new AliTRDPInfo();
	 UShort_t *entries   = new UShort_t[fNumberBinPRF];
	 Float_t  *sum       = new Float_t[fNumberBinPRF];
	 Float_t  *sumsquare = new Float_t[fNumberBinPRF];

	 for (Int_t nu = 0; nu < fNumberBinPRF; nu++) {
	   
	   entries[nu]           = ((AliTRDPInfo *) fVectorPRF->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu]
                                 + ((AliTRDPInfo *)  vVectorP2->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu];
	   
	   Double_t calcul       = ((((Double_t) ((AliTRDPInfo *) fVectorPRF->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetSum()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) fVectorPRF->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu]))
                                  + (((Double_t) ((AliTRDPInfo *)  vVectorP2->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetSum()[nu])
                                   * ((Double_t) ((AliTRDPInfo *)  vVectorP2->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu])))
                                 / ((Double_t) fPInfo->GetEntries()[nu]);

	   sum[nu]               = (Float_t) calcul;
	   
	   Double_t calculsquare = ((((Double_t) ((AliTRDPInfo *) fVectorPRF->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetSumSquare()[nu])
                                   * ((Double_t) ((AliTRDPInfo *) fVectorPRF->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu]))
                                  + (((Double_t) ((AliTRDPInfo *)  vVectorP2->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetSumSquare()[nu])
                                   * ((Double_t) ((AliTRDPInfo *)  vVectorP2->At(((AliTRDPlace *) fPlaPRF->At(k))->GetPlace()))->GetEntries()[nu])))
                                 / ((Double_t) fPInfo->GetEntries()[nu]);

	   sumsquare[nu]         = calculsquare;

	 }

	 // Set
	 fPInfo->SetSum(sum);
	 fPInfo->SetSumSquare(sumsquare);
	 fPInfo->SetEntries(entries);

	 // Nothing to do on PlaCT1
	 
	 // Update the vector VectorCT1
	 fVectorPRF->AddAt((TObject *) fPInfo,((AliTRDPlace *) fPlaPRF->At(k))->GetPlace());
	 
       }

     }

     // And at the end the vector in P2 but not in CH1
     for (Int_t k = 0; k < (Int_t) pPlaP2->GetEntriesFast(); k++) {
       
       // Look if PlaCT2[k] it is also in the second vector
       Int_t place = -1;
       for (Int_t j = 0; j < (Int_t) fPlaPRF->GetEntriesFast(); j++) {
	 if (((AliTRDPlace *) fPlaPRF->At(j))->GetPlace() == ((AliTRDPlace *) pPlaP2->At(k))->GetPlace()) {
	   place = j;
	   break;
	 }
       }

       // If not in the first vector
       if (place == -1) {

	 AliTRDPInfo *fPInfo = new AliTRDPInfo();
	 fPInfo = (AliTRDPInfo *) vVectorP2->At(((AliTRDPlace *) pPlaP2->At(k))->GetPlace());

	 // Add at the end of CH1
	 fPlaPRF->Add(((TObject *) pPlaP2->At(k)));
	 fVectorPRF->Add((TObject *) fPInfo);

       }
       
     }

   } 

   return kTRUE;

}   

//_____________________________________________________________________________
TTree *AliTRDCalibraVector::Sum2Trees(const Char_t *filename1
                                    , const Char_t *filename2
                                    , const Char_t *variablecali)
{
  //
  // It returns the sum of two trees with the name variablecali
  // in the files filenam1 and filename2 equivalent of merging two 2D histos
  // The name of the resulting tree is the same as the two input trees
  // variablecali can be treeCH2d, treePH2d or treePRF2d 
  //

  // Variables
  TChain    *treeChain   = new TChain(variablecali);
  TObjArray *vectorplace = new TObjArray();
  TObjArray *where       = new TObjArray();
  
  // First tree
  // Take the tree
  TFile *file1 = new TFile(filename1,"READ");
  TTree *tree1 = (TTree *) file1->Get(variablecali);

  gDirectory = gROOT;

  // Take the places
  vectorplace = ConvertTreeVector(tree1);

  // Say where it is in tree 1
  for (Int_t jui = 0; jui < (Int_t) vectorplace->GetEntriesFast(); jui++) {
    AliTRDPlace *placejui = new AliTRDPlace();
    placejui->SetPlace(jui);
    TObjArray *chainplace = new TObjArray();
    chainplace->Add((TObject *) placejui);
    where->Add((TObject *) chainplace);
  }

  // Add to the chain
  treeChain->Add(filename1);
  delete file1;

  // Second tree
  // Take the tree
  TFile *file2 = new TFile(filename2,"READ");
  TTree *tree2 = (TTree *) file2->Get(variablecali);

  gDirectory = gROOT;

  // Take the places
  TObjArray *vector2 = ConvertTreeVector(tree2);
  Int_t j = treeChain->GetEntries();

  for (Int_t jui = 0; jui < (Int_t) vector2->GetEntriesFast(); jui++) {
    // Search if already found
    Int_t place = SearchInTreeVector(vectorplace,((AliTRDPlace *) vector2->At(jui))->GetPlace());
    // Create a new element in the two std vectors
    if (place == -1) {
      AliTRDPlace *placejjui  = new AliTRDPlace();
      placejjui->SetPlace((j+jui));
      TObjArray   *chainplace = new TObjArray();
      chainplace->Add((TObject *) placejjui);
      vectorplace->Add((TObject *) (vector2->At(jui)));
      where->Add((TObject *) chainplace);
    }
    // Update the element at the place "place" in the std vector whereinthechain
    else {
      AliTRDPlace *placejjui  = new AliTRDPlace();
      placejjui->SetPlace((j+jui));
      TObjArray   *chainplace = ((TObjArray *) where->At(place));
      chainplace->Add((TObject *) placejjui);
      where->AddAt((TObject *) chainplace,place);
    }
  }

  // Add to the Chain
  treeChain->Add(filename2);
  delete file2; 

  // Take care of the profile
  const Char_t *pattern = "P";
  TTree *tree = 0x0;

  if (!strstr(variablecali,pattern)) {

    // Ready to read the chain
    TH1F *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);

    // Initialise the final tree
    Int_t group   = -1;
    TH1F *histsum = 0x0;
   
    tree = new TTree(variablecali,variablecali);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TH1F",&histsum,32000,0);

    // Init histsum
    if (treeChain->GetEntries() < 1) {
      return tree1; 
    }
    
    for (Int_t h = 0; h < (Int_t) vectorplace->GetEntriesFast(); h++) {
      group = ((AliTRDPlace *) vectorplace->At(h))->GetPlace();
      TObjArray *chainplace = ((TObjArray *) where->At(h));
      treeChain->GetEntry(((AliTRDPlace *) chainplace->At(0))->GetPlace());
      //Init for the first time
      if (h == 0)  {
	histsum = new TH1F("","",his->GetXaxis()->GetNbins()
                                ,his->GetXaxis()->GetBinLowEdge(1)
                                ,his->GetXaxis()->GetBinUpEdge(his->GetXaxis()->GetNbins()));
	histsum->Sumw2();
      }
      // Reset for each new group
      histsum->SetEntries(0.0);
      for (Int_t l = 0; l <= histsum->GetXaxis()->GetNbins(); l++) {
	histsum->SetBinContent(l,0.0);
	histsum->SetBinError(l,0.0);
      }
      histsum->Add(his,1);
      if ((Int_t) chainplace->GetEntriesFast() > 1) {
	for (Int_t s = 1; s < (Int_t) chainplace->GetEntriesFast(); s++) {
	  treeChain->GetEntry(((AliTRDPlace *) chainplace->At(s))->GetPlace());
	  histsum->Add(his,1);
	}
      }
      tree->Fill();
    }

  }
  else {

    // Ready to read the chain
    TGraphErrors *his = 0x0;
    treeChain->SetBranchAddress("histo",&his);
    
    // Initialise the final tree
    Int_t         group   = -1;
    TGraphErrors *histsum = 0x0;
    Double_t     *xref    = 0x0;
  
    tree = new TTree(variablecali,variablecali);
    tree->Branch("groupnumber",&group,"groupnumber/I");
    tree->Branch("histo","TGraphErrors",&histsum,32000,0);

    // Init histsum
    if (treeChain->GetEntries() < 1) {
      return tree1; 
    }

    for (Int_t h = 0; h < (Int_t) vectorplace->GetEntriesFast(); h++) {

      group = ((AliTRDPlace *) vectorplace->At(h))->GetPlace();
      TObjArray *chainplace = ((TObjArray *) where->At(h));
      treeChain->GetEntry(((AliTRDPlace *) chainplace->At(0))->GetPlace());
      //Init or reset for a new group
      Int_t nbins = his->GetN();
      Double_t *x;
      x    = new Double_t[nbins];
      xref = his->GetX();
      Double_t *ex;
      ex   = new Double_t[nbins];
      Double_t *y;
      y    = new Double_t[nbins];
      Double_t *ey;
      ey   = new Double_t[nbins];
     
      for (Int_t lo = 0; lo < nbins; lo++) {
	x[lo]  = xref[lo];
	ex[lo] = 0.0;
	y[lo]  = 0.0;
	ey[lo] = 0.0;
      }
      delete histsum;
      histsum = new TGraphErrors(nbins,x,y,ex,ey);

      // Add the first
      histsum = AddProfiles(his,histsum);
      if ((Int_t) chainplace->GetEntriesFast() > 1) {
	for (Int_t s = 1; s < (Int_t) chainplace->GetEntriesFast(); s++) {
	  treeChain->GetEntry(((AliTRDPlace *) chainplace->At(s))->GetPlace());
	  histsum = AddProfiles(his,histsum);
	}
      }

      tree->Fill();

    }

  }
    
  return tree;

}

//_____________________________________________________________________________
TGraphErrors *AliTRDCalibraVector::AddProfiles(TGraphErrors *hist1
                                             , TGraphErrors *hist2) const
{
  //
  // In the case of the vectors method we use TGraphErrors for PH and PRF
  // to be able to add the them after
  // Here we add the TGraphErrors  
  //

  // First TGraphErrors
  Int_t     nbins1 = hist1->GetN();
  Double_t *x1     = hist1->GetX();
  Double_t *ex1    = hist1->GetEX();
  Double_t *y1     = hist1->GetY();
  Double_t *ey1    = hist1->GetEY();

  TGraphErrors *rehist = new TGraphErrors(nbins1);

  // Second TGraphErrors
  Double_t *ex2    = hist2->GetEX();
  Double_t *y2     = hist2->GetY();
  Double_t *ey2    = hist2->GetEY();

  // Define the Variables for the new TGraphErrors
  Double_t x;
  Double_t ex;
  Double_t y;
  Double_t ey;
  
  for (Int_t k = 0; k < nbins1; k++) {
    Double_t nentries = 0.0;
    x  = x1[k];
    y  = 0.0;
    ey = 0.0;
    ex = 0.0;
    if ((ex2[k] == 0.0) && 
        (ex1[k] == 0.0)) {
      nentries = 0.0;
    }
    if ((ex2[k] == 0.0) && 
        (ex1[k]  > 0.0)) {
      nentries = ex1[k];
      y  = y1[k];
      ey = ey1[k];
      ex = ex1[k];
    }
    if ((ex2[k]  > 0.0) && 
        (ex1[k] == 0.0)) {
      nentries = ex2[k];
      y  = y2[k];
      ey = ey2[k];
      ex = ex2[k];
    }
    if ((ex2[k] > 0.0) && 
        (ex1[k] > 0.0)) { 
     nentries = ex1[k] + ex2[k];
     y  = ( y1[k]*ex1[k]+ y2[k]*ex2[k]) / nentries;
     ey = (ey1[k]*ex1[k]+ey2[k]*ex2[k]) / nentries;
     ex = nentries;
   }
   rehist->SetPoint(k,x,y);
   rehist->SetPointError(k,ex,ey);
 }

 return rehist;

}
