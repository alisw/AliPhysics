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
//___________________________________________________________________
TH1* AliCFEffGrid::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
{
  //
  // Make a projection along variable ivar1 (and ivar2 (and ivar3))
  //

  if (fSelNum<0 || fSelDen<0) {
    AliError("You must call CalculateEfficiency() first !");
    return 0x0;
  }
  const Int_t nDim = 3 ;
  Int_t dim[nDim] = {ivar1,ivar2,ivar3} ;
  
  THnSparse *hNum, *hDen, *ratio;
  TH1* h ;

  if (ivar3<0) {
    if (ivar2<0) {
      hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim-2,dim);
      hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim-2,dim);
      ratio = (THnSparse*)hNum->Clone();
      ratio->Divide(hNum,hDen,1.,1.,"B");
      h = ratio->Projection(0);
    }
    else{
      hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim-1,dim);
      hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim-1,dim);
      ratio = (THnSparse*)hNum->Clone();
      ratio->Divide(hNum,hDen,1.,1.,"B");
      h = ratio->Projection(1,0);
    }
  }
  else {
    hNum = ((AliCFGridSparse*)GetNum())->GetGrid()->Projection(nDim,dim);
    hDen = ((AliCFGridSparse*)GetDen())->GetGrid()->Projection(nDim,dim);
    ratio = (THnSparse*)hNum->Clone();
    ratio->Divide(hNum,hDen,1.,1.,"B");
    h = ratio->Projection(0,1,2);
  }

  delete hNum; delete hDen; delete ratio;
  return h ;
} 
//___________________________________________________________________
AliCFEffGrid* AliCFEffGrid::MakeSlice(Int_t nVars, const Int_t* vars, const Double_t* varMin, const Double_t* varMax, Bool_t useBins) const {
  //
  // returns a slice of the efficiency grid (slice is actually done on the container, and efficiency recomputed)
  //
  
  AliCFContainer* cont = fContainer->MakeSlice(nVars,vars,varMin,varMax,useBins);
  AliCFEffGrid  * eff  = new AliCFEffGrid(Form("%s_sliced",GetName()), Form("%s_sliced",GetTitle()), *cont);
  eff->CalculateEfficiency(fSelNum,fSelDen);
  return eff;
}
