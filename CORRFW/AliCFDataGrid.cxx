/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFDataGrid Class                                        //
// Class to handle observed data and correct them                     // 
//                                                                    //
// -- Author : S.Arcelli                                              //
//                                                                    //
//                                                                    //
//                                                                    //
//--------------------------------------------------------------------//
//
//
#include <TROOT.h>
#include <TMath.h>
#include <TFile.h>
#include <AliLog.h>
#include "AliCFDataGrid.h"

//____________________________________________________________________
ClassImp(AliCFDataGrid)

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid() : 
  AliCFGrid(),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // default constructor
  //
  SumW2(); //errors saved
}

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const Char_t* name,const Char_t* title) : 
  AliCFGrid(name,title),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // default constructor
  //
  SumW2(); //errors saved
}

//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Float_t *binLimitsIn) :  
  AliCFGrid(name,title,nVarIn,nBinIn,binLimitsIn),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // main constructor
  //
  SumW2();// errors saved
}
//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const Char_t* name, const Char_t* title, const AliCFContainer &c) :  
  AliCFGrid(name,title,c.GetNVar(),c.GetNBins(),c.GetBinLimits()),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // main constructor
  //
  SumW2();
   //assign the container;
  fContainer=&c;

}
//____________________________________________________________________
AliCFDataGrid::AliCFDataGrid(const AliCFDataGrid& data) :   AliCFGrid(),
  fSelData(-1),
  fContainer(0x0)
{
  //
  // copy constructor
  //
  ((AliCFDataGrid &)data).Copy(*this);
}

//____________________________________________________________________
AliCFDataGrid::~AliCFDataGrid()
{
  //
  // destructor
  //
}
//____________________________________________________________________
AliCFDataGrid &AliCFDataGrid::operator=(const AliCFDataGrid &c)
{
  //
  // assigment operator
  //
  if (this != &c)
    ((AliCFDataGrid &) c).Copy(*this);
  return *this;
} 
//____________________________________________________________________

void AliCFDataGrid::SetMeasured(Int_t istep)
{
  //
  // Deposit observed data over the grid
  //
  Int_t nEmptyBins=0;
  fSelData=istep;
  //Initially, set the corrected data to the measured data
  for(Int_t i=0;i<fNDim;i++){
    Float_t meas=fContainer->GetGrid(fSelData)->GetElement(i);
    Float_t dmeas=fContainer->GetGrid(fSelData)->GetElementError(i);
    SetElement(i,meas);
    SetElementError(i,dmeas);
    if(meas <=0)nEmptyBins++;
  }
  fNentriesTot=fNDim;
  AliInfo(Form("retrieving measured data from Container %s at selection step %i: %i empty bins were found.",fContainer->GetName(),fSelData,nEmptyBins));
} 
//____________________________________________________________________
void AliCFDataGrid::ApplyEffCorrection(const AliCFEffGrid &c)
{

  //
  // Apply the efficiency correction
  //
  if(c.GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot apply correction");
    return;
  }
  if(c.GetNDim()!=fNDim){
    AliInfo("Different number of dimension, cannot apply correction");
    return;
  }

  //Get the data
  Int_t ncorr=0;    
  Int_t nnocorr=0;    
  Float_t eff,deff,unc,dunc,corr,dcorr;
  //Apply the correction
  for(Int_t i=0;i<fNDim;i++){
    eff =c.GetElement(i);    
    deff =TMath::Sqrt(c.GetElementError(i));    
    unc =GetElement(i);    
    dunc =TMath::Sqrt(GetElementError(i));    
    if(eff>0 && unc>0){      
      ncorr++;
      corr=unc/eff;
      dcorr=TMath::Sqrt(dunc*dunc/unc/unc+deff*deff/eff/eff)*corr;
      SetElement(i,corr);
      SetElementError(i,dcorr*dcorr);
      
    } else{
      if(unc>0)nnocorr++;
      SetElement(i,0);
      SetElementError(i,0);
    }
  }
  AliInfo(Form("correction applied for %i cells in correction matrix of Container %s, having entries in Data Container %s.",ncorr,c.GetName(),GetName()));
  AliInfo(Form("No correction applied for %i empty bins in correction matrix of Container %s, having entries in Data Container %s. Their content in the corrected data container was set to zero",nnocorr,c.GetName(),GetName()));
}
//____________________________________________________________________
void AliCFDataGrid::ApplyBGCorrection(const AliCFDataGrid &c)
{

  //
  // Apply correction for background
  //
  if(c.GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot apply correction");
    return;
  }
  if(c.GetNDim()!=fNDim){
    AliInfo("Different number of dimension, cannot apply correction");
    return;
  }

  //Get the data
  Float_t bkg,dbkg,unc,dunc,corr,dcorr;

  //Apply the correction

  for(Int_t i=0;i<fNDim;i++){
    bkg =c.GetElement(i);    
    dbkg =TMath::Sqrt(c.GetElementError(i));    
    unc =GetElement(i);    
    dunc =TMath::Sqrt(GetElementError(i));    
    corr=unc-bkg;
    dcorr=TMath::Sqrt(unc+bkg); //stats only, check this
    SetElement(i,corr);
    SetElementError(i,dcorr*dcorr);
      
  }
}

//____________________________________________________________________
void AliCFDataGrid::Copy(TObject& eff) const
{
  // copy function

  Copy(eff);
  AliCFDataGrid& target = (AliCFDataGrid &) eff;
  target.fContainer=fContainer;
  target.fSelData=fSelData;

}
