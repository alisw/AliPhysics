/* $Id$ */
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
//--------------------------------------------------------------------//
//                                                                    //
// AliCFEffGrid Class                                                 //
// Class to handle efficiency grids                                   // 
//                                                                    //
// -- Author : S.Arcelli                                              //
//                                                                    //
//                                                                    //
//                                                                    //
//--------------------------------------------------------------------//
//
//
#include "TMath.h"
#include "AliLog.h"
#include "AliCFEffGrid.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

//____________________________________________________________________
ClassImp(AliCFEffGrid)

//____________________________________________________________________
AliCFEffGrid::AliCFEffGrid() : 
  AliCFGridSparse(),
  fContainer(0x0),
  fSelNum(-1),
  fSelDen(-1)
{
  //
  // default constructor
  //
}

//____________________________________________________________________
AliCFEffGrid::AliCFEffGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn) :
  AliCFGridSparse(name,title,nVarIn,nBinIn),
  fContainer(0x0),
  fSelNum(-1),
  fSelDen(-1)
{
  //
  // ctor
  //
  SumW2();
}
//____________________________________________________________________
AliCFEffGrid::AliCFEffGrid(const Char_t* name, const Char_t* title, const AliCFContainer &c) :  
  AliCFGridSparse(name,title,c.GetNVar(),c.GetNBins()),
  fContainer(NULL),
  fSelNum(-1),
  fSelDen(-1)
{
  //
  // main constructor
  //
  SumW2();
  //assign the container;
  fContainer=&c;
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) {
    Int_t nbins = c.GetNBins(iVar);
    Double_t* array=new Double_t[nbins+1] ;
    c.GetBinLimits(iVar,array);
    SetBinLimits(iVar,array);
    delete [] array ;
  }
  for (Int_t iVar=0; iVar<GetNVar(); iVar++) SetVarTitle(iVar,c.GetVarTitle(iVar));
}
//____________________________________________________________________
AliCFEffGrid::AliCFEffGrid(const AliCFEffGrid& eff) : 
  AliCFGridSparse(eff),
  fContainer(0x0),
  fSelNum(-1),
  fSelDen(-1)
{
  //
  // copy constructor
  //
  ((AliCFEffGrid &)eff).Copy(*this);
}

//____________________________________________________________________
AliCFEffGrid::~AliCFEffGrid()
{
  //
  // destructor
  //
}

//____________________________________________________________________
AliCFEffGrid &AliCFEffGrid::operator=(const AliCFEffGrid &eff)
{
  //
  // assigment operator
  //
  if (this != &eff) eff.Copy(*this);
  return *this;
} 

//____________________________________________________________________
void AliCFEffGrid::CalculateEfficiency(Int_t istep1,Int_t istep2, Option_t *option)
{
  //
  // Calculate the efficiency matrix and its error between selection
  // Steps istep1 and istep2
  //
  // 'option' is used as an argument for THnSparse::Divide
  // default is "B" : binomial error calculation
  //

  fSelNum=istep1;
  fSelDen=istep2;
  AliCFGridSparse *num=GetNum();
  AliCFGridSparse *den=GetDen();
  num->SumW2();
  den->SumW2();
  this->SumW2();
  this->Divide(num,den,1.,1.,option);
  SetTitle(Form("Efficiency: %s / %s",fContainer->GetStepTitle(istep1),fContainer->GetStepTitle(istep2)));

  AliInfo(Form("Efficiency calculated for steps %i and %i.",fSelNum,fSelDen));
} 
//_____________________________________________________________________
Double_t AliCFEffGrid::GetAverage() const 
{
  //
  // Get the average efficiency 
  //

  Double_t val=0;
  Double_t valnum=0;
  Double_t valden=0;

  THnSparse* num = ((AliCFGridSparse*)GetNum())->GetGrid() ;
  THnSparse* den = ((AliCFGridSparse*)GetDen())->GetGrid() ;

  for (Long_t iBin=0; iBin<num->GetNbins(); iBin++) valnum+=num->GetBinContent(iBin);
  for (Long_t iBin=0; iBin<den->GetNbins(); iBin++) valden+=den->GetBinContent(iBin);
  if (valden>0) val=valnum/valden;
  AliInfo(Form(" The Average Efficiency = %f ",val)); 
  return val;
} 
//_____________________________________________________________________
// Double_t AliCFEffGrid::GetAverage(const Double_t *varMin, const Double_t* varMax ) const 
// {
//   //
//   // Get ave efficiency in a range
//   // (may not work properly, should be modified)

//   Double_t val=0;
//   Int_t *indexMin = new Int_t[GetNVar()];
//   Int_t *indexMax = new Int_t[GetNVar()];
//   Int_t *index    = new Int_t[GetNVar()];
  
//   //Find out the min and max bins
  
//   for(Int_t i=0;i<GetNVar();i++){
//     Double_t xmin=varMin[i]; // the min values  
//     Double_t xmax=varMax[i]; // the max values  
//     Int_t nbins = GetNBins(i)+1;
// //     Double_t *bins = new Double_t[nbins];
// //     for (Int_t ibin =0; ibin<nbins; ibin++) {
// //       bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
// //     }
//     bins = GetBinLimits(i);
//     indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
//     indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
//     if(xmax>=bins[nbins-1]){
//       indexMax[i]=indexMax[i]-1;
//     }  
//     delete [] bins;
//   }
  
//   Double_t valnum=0;
//   Double_t valden=0;
//   for (Int_t i=0; i<fNDim; i++) {
//     for (Int_t j=0; j<GetNVar(); j++) index[j]=GetBinIndex(j,i);
//     Bool_t isIn=kTRUE;
//     for (Int_t j=0;j<GetNVar();j++) {
//       if(!(index[j]>=indexMin[j] && index[j]<=indexMax[j]))isIn=kFALSE;   
//     }
//     if(isIn){
//       valnum+=GetNum()->GetElement(i);
//       valden+=GetDen()->GetElement(i);
//     }
//   } 
//   delete [] index;
//   delete [] indexMin;
//   delete [] indexMax;
//   if(valden>0)val=valnum/valden;
//   AliInfo(Form(" the Average Efficiency = %f ",val)); 
//   return val;
// } 
//____________________________________________________________________
// void AliCFEffGrid::Copy(TObject& eff) const
// {
//   //
//   // copy function
//   //
//   Copy(eff);
//   AliCFEffGrid& target = (AliCFEffGrid &) eff;
  
//   target.fSelNum=fSelNum; 
//   target.fSelDen=fSelDen; 
//   if(fContainer)
//     target.fContainer=fContainer;
// }
//___________________________________________________________________
TH1D *AliCFEffGrid::Project(Int_t ivar) const
{
  //
  // Make a 1D projection along variable ivar 
  //
  
  if (fSelNum<0 || fSelDen<0) {
    AliError("You must call CalculateEfficiency() first !");
    return 0x0;
  }
  const Int_t nDim = 1 ;
  Int_t dim[nDim] = {ivar} ;
  THnSparse* hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim,dim);
  THnSparse* hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim,dim);
  THnSparse* ratio = (THnSparse*)hNum->Clone();
  ratio->Divide(hNum,hDen,1.,1.,"B");
  delete hNum; delete hDen;
  TH1D* h = ratio->Projection(0);
  h->SetXTitle(GetVarTitle(ivar));
  h->SetName(Form("%s_proj-%s",GetName(),GetVarTitle(ivar)));
  h->SetTitle(Form("%s projected on %s",GetTitle(),GetVarTitle(ivar)));
  return h ;
} 
//___________________________________________________________________
TH2D *AliCFEffGrid::Project(Int_t ivar1,Int_t ivar2) const
{
  //
  // Make a 2D projection along variable ivar1,ivar2 
  //
  
  if (fSelNum<0 || fSelDen<0) {
    AliError("You must call CalculateEfficiency() first !");
    return 0x0;
  }
  const Int_t nDim = 2 ;
  Int_t dim[nDim] = {ivar1,ivar2} ;
  THnSparse* hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim,dim);
  THnSparse* hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim,dim);
  THnSparse* ratio = (THnSparse*)hNum->Clone();
  ratio->Divide(hNum,hDen,1.,1.,"B");
  delete hNum; delete hDen;
  TH2D* h = ratio->Projection(1,0);
  h->SetXTitle(GetVarTitle(ivar1));
  h->SetYTitle(GetVarTitle(ivar2));
  h->SetName(Form("%s_proj-%s,%s",GetName(),GetVarTitle(ivar1),GetVarTitle(ivar2)));
  h->SetTitle(Form("%s projected on %s-%s",GetTitle(),GetVarTitle(ivar1),GetVarTitle(ivar2)));
  return h;
} 
//___________________________________________________________________
TH3D *AliCFEffGrid::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
{
  //
  // Make a 3D projection along variable ivar1,ivar2,ivar3 
  //

  if (fSelNum<0 || fSelDen<0) {
    AliError("You must call CalculateEfficiency() first !");
    return 0x0;
  }
  const Int_t nDim = 3 ;
  Int_t dim[nDim] = {ivar1,ivar2,ivar3} ;
  THnSparse* hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim,dim);
  THnSparse* hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim,dim);
  THnSparse* ratio = (THnSparse*)hNum->Clone();
  ratio->Divide(hNum,hDen,1.,1.,"B");
  delete hNum; delete hDen;
  TH3D* h = ratio->Projection(0,1,2);
  h->SetXTitle(GetVarTitle(ivar1));
  h->SetYTitle(GetVarTitle(ivar2));
  h->SetZTitle(GetVarTitle(ivar3));
  h->SetName(Form("%s_proj-%s,%s,%s",GetName(),GetVarTitle(ivar1),GetVarTitle(ivar2),GetVarTitle(ivar3)));
  h->SetTitle(Form("%s projected on %s-%s-%s",GetTitle(),GetVarTitle(ivar1),GetVarTitle(ivar2),GetVarTitle(ivar3)));
  return h;
} 
//___________________________________________________________________
AliCFEffGrid* AliCFEffGrid::MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Int_t numStep, Int_t denStep, Bool_t useBins) const {
  //
  // Makes a slice along the "nVars" variables defined in the array "vars[nVars]" for all the container steps.
  // The ranges of ALL the container variables must be defined in the array varMin[fNVar] and varMax[fNVar].
  // This function returns the efficiency relative to this new 'sliced' container, between steps defined in numStep and denStep
  // If useBins=true, varMin and varMax are taken as bin numbers
  //
  
  AliCFContainer* cont = fContainer->MakeSlice(nVars,vars,varMin,varMax,useBins);
  AliCFEffGrid  * eff  = new AliCFEffGrid(Form("%s_sliced",GetName()), Form("%s_sliced",GetTitle()), *cont);
  eff->CalculateEfficiency(numStep,denStep);
  return eff;
}
