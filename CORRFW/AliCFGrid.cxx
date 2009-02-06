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
//---------------------------------------------------------------------//
//                                                                     //
// AliCFGrid Class                                                     //
// Class to accumulate data on an N-dimensional grid, to be used       //
// as input to get corrections for Reconstruction & Trigger efficiency // 
// The class uses a one-dimensional array of floats to store the grid  //     
// --Author : S.Arcelli                                                //
//   Still to be done:                                                 //
// --Implement methods to merge cells                                  //
// --Interpolate among bins in a range                                 // 
// This implementation will be eventually replaced by AliCFGridSparse  //
//---------------------------------------------------------------------//
//
//
#include "AliLog.h"
#include "AliCFGrid.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

//____________________________________________________________________
ClassImp(AliCFGrid)

//____________________________________________________________________
AliCFGrid::AliCFGrid() : 
  AliCFVGrid(),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  // default constructor
}
//____________________________________________________________________
AliCFGrid::AliCFGrid(const Char_t* name, const Char_t* title) : 
  AliCFVGrid(name,title),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  // default constructor
}

//____________________________________________________________________
AliCFGrid::AliCFGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Double_t *binLimitsIn) :  
  AliCFVGrid(name,title,nVarIn,nBinIn,binLimitsIn),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  //
  // main constructor
  //
 
  //The over/underflows
  fNunfl=new Float_t[fNVar];
  fNovfl= new Float_t[fNVar];


  //Initialization
  fNentriesTot =0;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j] =0;
    fNovfl[j] =0;
  }


  // the grid
 
  fData = new Float_t[fNDim]; 

  //Initialization
  for(Int_t j=0;j<fNDim;j++){
    fData[j] =0;
  }
 
}

//____________________________________________________________________
AliCFGrid::AliCFGrid(const AliCFGrid& c) : 
  AliCFVGrid(c),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  //
  // copy constructor
  //
  ((AliCFGrid &)c).Copy(*this);
}

//____________________________________________________________________
AliCFGrid::~AliCFGrid()
{
  //
  // destructor
  //
  if(fNunfl)delete fNunfl;
  if(fNovfl)delete fNovfl;
  if(fData)delete fData;
  if(fSumW2)delete fErr2;

}
//____________________________________________________________________
AliCFGrid &AliCFGrid::operator=(const AliCFGrid &c)
{
  //
  // assigment operator
  //
  if (this != &c)
    ((AliCFGrid &) c).Copy(*this);
  return *this;
} 
//____________________________________________________________________
Float_t AliCFGrid::GetElement(Int_t bin) const
{
  //
  // Returns grid element bin
  //
  if(bin>=fNDim){
    AliInfo(Form(" element index outside the grid, return -1"));
    return -1.;
  }
  return fData[bin];
}
//____________________________________________________________________
Float_t AliCFGrid::GetElement(Int_t *bin) const
{
 //
  // Get the content in a bin corresponding to a set of bin indexes
  //
  //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
  for(Int_t i=0;i<fNVar; i++)fIndex[i]=bin[i]-1;
  Int_t ind =GetBinIndex(fIndex);
  return GetElement(ind);
}  
//____________________________________________________________________
Float_t AliCFGrid::GetElement(Double_t *var) const
{
  //
  // Get the content in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);

    //underflows

    if(var[i] < bins[0]){
      unfl=1;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl=1;
    }
    delete [] bins;
  }


  //move to the TH/THnSparse convention in N-dim bin numbering 
  for(Int_t i=0;i<fNVar; i++)fIndex[i]+=1;

  if(!(ovfl==1 || unfl==1)){
    return GetElement(fIndex);
  }
  else{
    AliInfo(Form(" input variables outside the grid, return -1"));
    return -1.;
  }
} 
//____________________________________________________________________
Float_t AliCFGrid::GetElementError(Int_t iel) const
{
  //
  // Return the error on grid element iel
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, return -1"));
    return -1.;
  }
  if(fSumW2)return TMath::Sqrt(fErr2[iel]);
  return TMath::Sqrt(fData[iel]);
}
//____________________________________________________________________
Float_t AliCFGrid::GetElementError(Int_t *bin) const
{
  //
  // Get the error in a bin corresponding to a set of bin indeces
  //
  //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
  for(Int_t i=0;i<fNVar; i++)fIndex[i]=bin[i]-1;
  Int_t ind =GetBinIndex(fIndex);
  return GetElementError(ind);

}
//____________________________________________________________________
Float_t AliCFGrid::GetElementError(Double_t *var) const
{
  //
  // Get the error in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);
    //underflows

    if(var[i] < bins[0]){
      unfl=1;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl=1;
    }
    delete [] bins;
  }

  //move to the TH/THnSparse convention in N-dim bin numbering 
  for(Int_t i=0;i<fNVar; i++)fIndex[i]+=1;
  if(!(ovfl==1 || unfl==1)){
    return GetElementError(fIndex);
  }
  else{
    AliInfo(Form(" input variables outside the grid, return -1"));
    return -1.;
  }
} 
//____________________________________________________________________
void AliCFGrid::SetElement(Int_t iel, Float_t val)
{
  //
  // Sets grid element iel to val
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, no value set"));
  }else {
    fData[iel]=val;
  }
}
//____________________________________________________________________
void AliCFGrid::SetElement(Int_t *bin, Float_t val)
{
  //
  // Sets grid element of bin indeces bin to val
  //
  //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
  for(Int_t i=0;i<fNVar; i++)fIndex[i]=bin[i]-1;
  Int_t ind =GetBinIndex(fIndex);
  SetElement(ind,val);
}
//____________________________________________________________________
void AliCFGrid::SetElement(Double_t *var, Float_t val) 
{
  //
  // Set the content in a bin to value val corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);
    //underflows

    if(var[i] < bins[0]){
      unfl=1;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl=1;
    }
    delete [] bins;
  }

  //move to the TH/THnSparse convention in N-dim bin numbering 
  for(Int_t i=0;i<fNVar; i++)fIndex[i]+=1;
  if(!(ovfl==1 || unfl==1)){
    SetElement(fIndex,val);
  }
  else{
    AliInfo(Form(" input variables outside the grid, no value set"));
  }
} 
//____________________________________________________________________
void AliCFGrid::SetElementError(Int_t iel, Float_t val) 
{
  //
  // Set squared error on grid element iel to val*val
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, no value set"));
    return;
  }
  if(!fErr2)SumW2();
  fErr2[iel]=val*val;   
}
//____________________________________________________________________
void AliCFGrid::SetElementError(Int_t *bin, Float_t val) 
{
  //
  // Set squared error to val on grid element of bin indeces bin
  //
  //-1 is to move from TH/ThnSparse N-dim bin convention to one in AliCFFrame
  for(Int_t i=0;i<fNVar; i++)fIndex[i]=bin[i]-1;
  Int_t ind =GetBinIndex(fIndex);
  SetElementError(ind,val);
}
//____________________________________________________________________
void AliCFGrid::SetElementError(Double_t *var, Float_t val)
{
  //
  // Set squared error to val in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);
    //underflows

    if(var[i] < bins[0]){
      unfl=1;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl=1;
    }
    delete [] bins;
  }

  //move to the TH/THnSparse convention in N-dim bin numbering 
  for(Int_t i=0;i<fNVar; i++)fIndex[i]+=1;

  if(!(ovfl==1 || unfl==1)){
    SetElementError(fIndex,val);
  }
  else{
    AliInfo(Form(" input variables outside the grid, no value set"));
  }
} 
//____________________________________________________________________
void AliCFGrid::Fill(Double_t *var, Double_t weight)
{

  //
  // Fill the grid,
  // given a set of values of the input variable, 
  // with weight (by default w=1)
  //


  Int_t isunfl=0;  
  Int_t isovfl=0;  
  Int_t *unfl=new Int_t[fNVar];  
  Int_t *ovfl=new Int_t[fNVar];  

  for(Int_t i=0;i<fNVar;i++){
    unfl[i]=0;
    ovfl[i]=0;
  }

  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Double_t *bins=new Double_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);
    //underflows

    if(var[i] < bins[0]){
      unfl[i]=1;  
      isunfl=1;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl[i]=1;  
      isovfl=1;
    }
    delete [] bins;
  }

  //exclusive under/overflows

  for(Int_t i=0;i<fNVar;i++){
    Bool_t add=kTRUE;
    for(Int_t j=0;j<fNVar;j++){
      if(i==j)continue;
      if(!(unfl[j]==0 && ovfl[j]==0))add=kFALSE;
    }
    if(add && unfl[i]==1)fNunfl[i]++;
    if(add && ovfl[i]==1)fNovfl[i]++;
  }

  delete [] unfl;
  delete [] ovfl;

  // Total number of entries, overflows and underflows

  fNentriesTot++;

  //if not ovfl/unfl, fill the element  

  if(!(isovfl==1 || isunfl==1)){
    Int_t ind =GetBinIndex(fIndex);
    fData[ind]+=weight;
    if(fSumW2)fErr2[ind]+=(weight*weight);
  }
} 
//____________________________________________________________________
Float_t AliCFGrid::GetOverFlows( Int_t ivar) const {
  //
  // Get overflows in variable var 
  //
  return fNovfl[ivar];
} 
//____________________________________________________________________
Float_t AliCFGrid::GetUnderFlows( Int_t ivar) const {
  //
  // Get overflows in variable var 
  //
  return fNunfl[ivar];
} 
//____________________________________________________________________
Float_t AliCFGrid::GetEntries( ) const {
  //
  // Get total entries (in grid + over/underflows) 
  //
  return fNentriesTot;
} 
//___________________________________________________________________
TH1D *AliCFGrid::Project(Int_t ivar) const
{
  //
  // Make a 1D projection along variable ivar 


  Int_t nbins =fNVarBins[ivar];
  Float_t *bins = new Float_t[nbins+1];    
  for (Int_t i=0;i<=fNVar;i++){
  }
  for(Int_t ibin =0;ibin<nbins+1;ibin++){
    bins[ibin] = fVarBinLimits[ibin+fOffset[ivar]];
  }

  char pname[40];
  sprintf(pname,"%s%s_%i",GetName(),"_proj1D_var", ivar);
  char htitle[40];
  sprintf(htitle,"%s%s_%i",GetName(),"_proj1D_var", ivar);

  TH1D *proj1D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH1D")) {
    proj1D = (TH1D*)obj;
    proj1D->Reset();
  }

  if(!proj1D){
    proj1D =new TH1D(pname,htitle, nbins, bins);
  }  

  delete [] bins;
  Float_t sum=0;
  Float_t *data= new Float_t[nbins];
  Float_t *err= new Float_t[nbins];
  
  for(Int_t ibin=0;ibin<nbins;ibin++)data[ibin]=0;
  for(Int_t ibin=0;ibin<nbins;ibin++)err[ibin]=0;
  for(Int_t iel=0;iel<fNDim;iel++){
      data[GetBinIndex(ivar,iel)]+=fData[iel];
      if(fSumW2)err[GetBinIndex(ivar,iel)]+=fErr2[iel];
  }

  for(Int_t ibin =0;ibin<nbins;ibin++){
    proj1D->SetBinContent(ibin+1,data[ibin]);
    proj1D->SetBinError(ibin+1,TMath::Sqrt(data[ibin]));
    if(fSumW2)proj1D->SetBinError(ibin+1,TMath::Sqrt(err[ibin]));
    sum+=data[ibin];
  }

  delete [] data;
  delete [] err;
  proj1D->SetBinContent(nbins+1,GetOverFlows(ivar));
  proj1D->SetBinContent(0,GetUnderFlows(ivar));
  proj1D->SetEntries(fNentriesTot);
  return proj1D;
} 

//___________________________________________________________________
TH2D *AliCFGrid::Project(Int_t ivar1, Int_t ivar2) const
{
  //
  // Make a 2D projection along variable ivar 

  Int_t nbins1 =fNVarBins[ivar1];
  Int_t nbins2 =fNVarBins[ivar2];

  Float_t *bins1 = new Float_t[nbins1+1];         
  Float_t *bins2 = new Float_t[nbins2+1];     

  for(Int_t ibin =0;ibin<nbins1+1;ibin++){
    bins1[ibin] = fVarBinLimits[ibin+fOffset[ivar1]];
  }
  for(Int_t ibin =0;ibin<nbins2+1;ibin++){
    bins2[ibin] = fVarBinLimits[ibin+fOffset[ivar2]];
  }

  char pname[40];
  sprintf(pname,"%s%s_%i_%i",GetName(),"_proj2D_var",ivar1,ivar2);
  char htitle[40];
  sprintf(htitle,"%s%s_%i_%i",GetName(),"_proj2D_var",ivar1,ivar2);

  TH2D *proj2D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH2D")) {
    proj2D = (TH2D*)obj;
    proj2D->Reset();
  }

  if(!proj2D){
    proj2D =new TH2D(pname,htitle, nbins1, bins1,nbins2,bins2);
  }  

  delete [] bins1;
  delete [] bins2;


  Float_t sum=0;
  Float_t **data=new Float_t*[nbins1];
  Float_t *data2=new Float_t[nbins1*nbins2];
  Float_t **err=new Float_t*[nbins1];
  Float_t *err2=new Float_t[nbins1*nbins2];
  for(Int_t i=0;i<nbins1;i++)data[i] = data2+i*nbins2;
  for(Int_t i=0;i<nbins1;i++)err[i] = err2+i*nbins2;

  for(Int_t ibin1 =0;ibin1<nbins1;ibin1++){
    for(Int_t ibin2 =0;ibin2<nbins2;ibin2++){
      data[ibin1][ibin2]=0;
      err[ibin1][ibin2]=0;
    }
  }

  for(Int_t iel=0;iel<fNDim;iel++){
    data[GetBinIndex(ivar1,iel)][GetBinIndex(ivar2,iel)]+=fData[iel];
    if(fSumW2)err[GetBinIndex(ivar1,iel)][GetBinIndex(ivar2,iel)]+=fErr2[iel];
  }

  for(Int_t ibin1 =0;ibin1<nbins1;ibin1++){
    for(Int_t ibin2 =0;ibin2<nbins2;ibin2++){
      proj2D->SetBinContent(ibin1+1,ibin2+1,data[ibin1][ibin2]);
      proj2D->SetBinError(ibin1+1,ibin2+1,TMath::Sqrt(data[ibin1][ibin2]));
      if(fSumW2)proj2D->SetBinError(ibin1+1,ibin2+1,TMath::Sqrt(err[ibin1][ibin2]));
      sum+=data[ibin1][ibin2];
    }

  } 
  delete data;
  delete data2;
  delete err;
  delete err2;
  proj2D->SetBinContent(0,nbins2/2,GetUnderFlows(ivar1));
  proj2D->SetBinContent(nbins1+1,nbins2/2,GetOverFlows(ivar1));
  proj2D->SetBinContent(nbins1/2,0,GetUnderFlows(ivar2));
  proj2D->SetBinContent(nbins1/2,nbins2+1,GetOverFlows(ivar2));
  proj2D->SetEntries(fNentriesTot);
  return proj2D;
} 
//___________________________________________________________________
TH3D *AliCFGrid::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
{
  //
  // Make a 3D projection along variable ivar 

  Int_t nbins1 =fNVarBins[ivar1];
  Int_t nbins2 =fNVarBins[ivar2];
  Int_t nbins3 =fNVarBins[ivar3];

  Float_t *bins1 = new Float_t[nbins1+1];         
  Float_t *bins2 = new Float_t[nbins2+1];     
  Float_t *bins3 = new Float_t[nbins3+1];     

  for(Int_t ibin =0;ibin<nbins1+1;ibin++){
    bins1[ibin] = fVarBinLimits[ibin+fOffset[ivar1]];
  }
  for(Int_t ibin =0;ibin<nbins2+1;ibin++){
    bins2[ibin] = fVarBinLimits[ibin+fOffset[ivar2]];
  }
  for(Int_t ibin =0;ibin<nbins3+1;ibin++){
    bins3[ibin] = fVarBinLimits[ibin+fOffset[ivar3]];
  }

  char pname[40];
  sprintf(pname,"%s%s_%i_%i_%i",GetName(),"_proj3D_var",ivar1,ivar2,ivar3);
  char htitle[40];
  sprintf(htitle,"%s%s_%i_%i_%i",GetName(),"_proj3D_var",ivar1,ivar2,ivar3);

  TH3D *proj3D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH3D")) {
    proj3D = (TH3D*)obj;
    proj3D->Reset();
  }

  if(!proj3D){
    proj3D =new TH3D(pname,htitle, nbins1,bins1,nbins2,bins2,nbins3,bins3);
  }  

  delete [] bins1;
  delete [] bins2;
  delete [] bins3;


  Float_t sum=0;
  Float_t ***data=new Float_t**[nbins1];
  Float_t **data2=new Float_t*[nbins1*nbins2];
  Float_t *data3=new Float_t[nbins1*nbins2*nbins3];
  Float_t ***err=new Float_t**[nbins1];
  Float_t **err2=new Float_t*[nbins1*nbins2];
  Float_t *err3=new Float_t[nbins1*nbins2*nbins3];
  for(Int_t i=0;i<nbins1;i++)data[i] = data2+i*nbins2;
  for(Int_t i=0;i<nbins1;i++)err[i] = err2+i*nbins2;
  for(Int_t i=0;i<nbins1;i++){
    for(Int_t j=0;j<nbins2;j++){
      data[i][j] = data3+i*nbins2*nbins3+j*nbins3;
      err[i][j] = err3+i*nbins2*nbins3+j*nbins3;
    }
  }
  for(Int_t ibin1 =0;ibin1<nbins1;ibin1++){
    for(Int_t ibin2 =0;ibin2<nbins2;ibin2++){
      for(Int_t ibin3 =0;ibin3<nbins3;ibin3++){
	data[ibin1][ibin2][ibin3]=0;
	err[ibin1][ibin2][ibin3]=0;
      }
    }
  }
  
  for(Int_t iel=0;iel<fNDim;iel++){
    data[GetBinIndex(ivar1,iel)][GetBinIndex(ivar2,iel)][GetBinIndex(ivar3,iel)]+=fData[iel];
    if(fSumW2)err[GetBinIndex(ivar1,iel)][GetBinIndex(ivar2,iel)][GetBinIndex(ivar3,iel)]+=fErr2[iel];
  }

  for(Int_t ibin1 =0;ibin1<nbins1;ibin1++){
    for(Int_t ibin2 =0;ibin2<nbins2;ibin2++){
      for(Int_t ibin3 =0;ibin3<nbins3;ibin3++){
	proj3D->SetBinContent(ibin1+1,ibin2+1,ibin3+1,data[ibin1][ibin2][ibin3]);
	proj3D->SetBinError(ibin1+1,ibin2+1,ibin3+1,TMath::Sqrt(data[ibin1][ibin2][ibin3]));
	if(fSumW2)proj3D->SetBinError(ibin1+1,ibin2+1,ibin3+1,TMath::Sqrt(err[ibin1][ibin2][ibin3]));
	sum+=data[ibin1][ibin2][ibin3];
      }
    }
  } 

  delete data;
  delete data2;
  delete data3;
  delete err;
  delete err2;
  delete err3;

  proj3D->SetEntries(fNentriesTot);
  return proj3D;
} 

//___________________________________________________________________
TH1D *AliCFGrid::Slice(Int_t ivar, Double_t *varMin, Double_t* varMax) const
{
  //
  // Make a slice along variable ivar in range [varMin,varMax]


  Int_t nbins =fNVarBins[ivar];
  Float_t *bins = new Float_t[nbins+1];    
  for (Int_t i=0;i<=fNVar;i++){
  }
  for(Int_t ibin =0;ibin<nbins+1;ibin++){
    bins[ibin] = fVarBinLimits[ibin+fOffset[ivar]];
  }

  char pname[40];
  sprintf(pname,"%s%s_%i",GetName(),"_proj1D_var", ivar);
  char htitle[40];
  sprintf(htitle,"%s%s_%i",GetName(),"_proj1D_var", ivar);

  TH1D *proj1D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH1D")) {
    proj1D = (TH1D*)obj;
    proj1D->Reset();
  }

  if(!proj1D){
    proj1D =new TH1D(pname,htitle, nbins, bins);
  }  

  delete [] bins;


  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];


  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Double_t xmin=varMin[i]; // the min values  
    Double_t xmax=varMax[i]; // the max values  
    Int_t nBins=fNVarBins[i]+1;
    Double_t *Bins=new Double_t[nBins];
    for(Int_t ibin =0;ibin<nBins;ibin++){
     Bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nBins,Bins,xmin);
    indexMax[i] = TMath::BinarySearch(nBins,Bins,xmax);
    delete [] Bins;
  }

  Float_t sum=0;
  Float_t *data= new Float_t[nbins];
  for(Int_t ibin=0;ibin<nbins;ibin++)data[ibin]=0;

  Int_t *index= new Int_t[fNVar];
  Int_t ielmin=GetBinIndex(indexMin);
  Int_t ielmax=GetBinIndex(indexMax);


  for(Int_t iel=ielmin;iel<=ielmax;iel++){
    GetBinIndex(iel,index);
    Bool_t isIn=kTRUE;
    for (Int_t j=0;j<fNVar;j++){
      if(!(index[j]>=indexMin[j] && index[j]<=indexMax[j]))isIn=kFALSE;   
      break;
    }
    if(isIn)data[GetBinIndex(ivar,iel)]+=fData[iel];
  }

  delete [] index;
  delete [] indexMin;
  delete [] indexMax;


  for(Int_t ibin =0;ibin<nbins;ibin++){
    proj1D->SetBinContent(ibin+1,data[ibin]);
    proj1D->SetBinError(ibin+1,TMath::Sqrt(data[ibin]));
    sum+=data[ibin];
  }

  delete [] data;

  proj1D->SetEntries(sum);
  return proj1D;
} 


//____________________________________________________________________
void AliCFGrid::Add(AliCFVGrid* aGrid, Double_t c)
{
  //
  //add aGrid to the current one
  //

  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot add the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot add the grids!");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();

  for(Int_t iel=0;iel<fNDim;iel++){
    fData[iel]+=(c*aGrid->GetElement(iel));
    if(fSumW2){
      Float_t err=aGrid->GetElementError(iel);  
      fErr2[iel]+=c*c*err*err;
    }
  }

  //Add entries, overflows and underflows

  fNentriesTot+= c*aGrid->GetEntries();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]+= c*aGrid->GetUnderFlows(j);
    fNovfl[j]+= c*aGrid->GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Add(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1,Double_t c2)
{
  //
  //add aGrid1 and aGrid2
  //

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot add the grids");
    return;
  } 
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot add the grids!");
    return;
  } 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();

  Float_t cont1,cont2,err1,err2;  

  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=aGrid1->GetElement(iel);
    cont2=aGrid2->GetElement(iel);
    SetElement(iel,c1*cont1+c2*cont2);
    if(fSumW2){
      err1=aGrid1->GetElementError(iel);
      err2=aGrid2->GetElementError(iel);
      SetElementError(iel,TMath::Sqrt(c1*c1*err1*err1+c2*c2*err2*err2));
    }
  }

  //Add entries, overflows and underflows

  fNentriesTot= c1*aGrid1->GetEntries()+c2*aGrid2->GetEntries();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= c1*aGrid1->GetUnderFlows(j)+c2*aGrid2->GetUnderFlows(j);
    fNovfl[j]= c1*aGrid1->GetOverFlows(j)+c2*aGrid2->GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Multiply(AliCFVGrid* aGrid, Double_t c)
{
  //
  //multiply grid aGrid by the current one
  //

  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot multiply the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot multiply the grids!");
    return;
  } 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();
  
  Float_t cont1,cont2,err1,err2;  

  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=GetElement(iel);  
    cont2=c*aGrid->GetElement(iel);  
    SetElement(iel,cont1*cont2);
    if(fSumW2){
      err1=GetElementError(iel);  
      err2=aGrid->GetElementError(iel);  
      SetElementError(iel,TMath::Sqrt(c*c*(cont2*cont2*err1*err1+cont1*cont1*err2*err2)));
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Multiply(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1, Double_t c2)
{
  //
  //multiply grids aGrid1 and aGrid2
  //

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot multiply the grids");
    return;
  } 
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot multiply the grids!");
    return;
  } 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();

  Float_t cont1,cont2,err1,err2;  
  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=c1*aGrid1->GetElement(iel);  
    cont2=c2*aGrid2->GetElement(iel);  
    SetElement(iel,cont1*cont2);
    if(fSumW2){
      err1=aGrid1->GetElementError(iel);  
      err2=aGrid2->GetElementError(iel);  
      SetElementError(iel,TMath::Sqrt(c1*c1*c2*c2*(cont2*cont2*err1*err1+cont1*cont1*err2*err2)));
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Divide(AliCFVGrid* aGrid, Double_t c)
{
  //
  //divide current grid by grid aGrid
  //

  if(aGrid->GetNVar()!=fNVar){
    AliInfo("Different number of variables, cannot divide the grids");
    return;
  } 
  if(aGrid->GetNDim()!=fNDim){
    AliInfo("Different number of dimensions, cannot divide the grids!");
    return;
  } 
 if(!c){AliInfo(Form("c is %f, cannot divide!",c)); return;} 
  
  if(!fSumW2  && aGrid->GetSumW2())SumW2();

  Float_t cont1,cont2,err1,err2,den;
  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=GetElement(iel);  
    cont2=aGrid->GetElement(iel);
    if(cont2)SetElement(iel,cont1/(c*cont2));
    else SetElement(iel,0);
    if(fSumW2){
      err1=GetElementError(iel);  
      err2=aGrid->GetElementError(iel);  
      if(!cont2){SetElementError(iel,0.); continue;}
      den=cont2*cont2*cont2*c*c;
      SetElementError(iel,TMath::Sqrt((cont2*cont2*err1*err1+cont1*cont1*err2*err2)/den));
    }
  }
  
  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Divide(AliCFVGrid* aGrid1, AliCFVGrid* aGrid2, Double_t c1,Double_t c2, Option_t *option)
{
  //
  //divide grids aGrid1,aGrid2
  //

  TString opt = option;
  opt.ToUpper();

  if(fNVar!=aGrid1->GetNVar()|| fNVar!=aGrid2->GetNVar()){
    AliInfo("Different number of variables, cannot divide the grids");
    return;
  }
  if(fNDim!=aGrid1->GetNDim()|| fNDim!=aGrid2->GetNDim()){
    AliInfo("Different number of dimensions, cannot divide the grids!");
    return;
  }
  if(!c2){AliInfo(Form("c2 is %f, cannot divide!",c2)); return;} 
  
  if(!fSumW2  && (aGrid1->GetSumW2() || aGrid2->GetSumW2()))SumW2();

  Float_t cont1,cont2,err1,err2,r,den;
 
  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=aGrid1->GetElement(iel);  
    cont2=aGrid2->GetElement(iel);  
    if(cont2)SetElement(iel,c1*cont1/(c2*cont2));
    else SetElement(iel,0);
     if(fSumW2){
      err1=aGrid1->GetElementError(iel);  
      err2=aGrid2->GetElementError(iel);  
      if(!cont2){SetElementError(iel,0.); continue;}
      if (opt.Contains("B")){
	if(cont1!=cont2){
	  r=cont1/cont2;	    
	  SetElementError(iel,TMath::Sqrt(TMath::Abs(((1.-2.*r)*err1*err1+r*r*err2*err2)/(cont2*cont2))));
	}else{
	  SetElementError(iel,0.);
	}
      }else{
        den=cont2*cont2*cont2*cont2*c2*c2;
	SetElementError(iel,TMath::Sqrt(c1*c1*(cont2*cont2*err1*err1+cont1*cont1*err2*err2)/den));
      }
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetOverFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Rebin(const Int_t* group)
{
  //
  // Not yet implemented
  //
  for(Int_t i=0;i<fNVar;i++){
    if(group[i]!=1)AliInfo(Form(" merging bins along dimension %i in groups of %i bins", i,group[i]));
  }
  AliInfo(Form("This method was so far not implemented for AliCFGrid, but it is available for AliCFGridSparse"));

}
//____________________________________________________________________
void AliCFGrid::SumW2()
{
  //
  //set calculation of the squared sum of the weighted entries
  //
  if(!fSumW2){
    fErr2=new Float_t [fNDim];
    //init....
    for(Int_t iel=0;iel<fNDim;iel++){
      fErr2[iel]=fData[iel];
    }
  }

  fSumW2=kTRUE;
}
//____________________________________________________________________
void AliCFGrid::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFGrid& target = (AliCFGrid &) c;

  target.fNentriesTot = fNentriesTot;
  if (fNunfl)
    target.fNunfl = fNunfl;
  if (fNovfl)
    target.fNovfl = fNovfl;
  if (fData)
    target.fData = fData;
  if (fErr2)
    target.fErr2 = fErr2;
  
}
