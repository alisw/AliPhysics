/* $Id$ */

//--------------------------------------------------------------------//
//                                                                    //
// AliCFGrid Class                                                 //
// Class to accumulate data on an N-dimensional grid, to be used      //
// as input to get corrections for Reconstruction & Trigger efficiency// 
//                                                                    //
// -- Author : S.Arcelli                                              //
// Still to be done:                                                  //
// --Implement methods to merge cells                                 //
// --Interpolate among bins in a range                                // 
//--------------------------------------------------------------------//
//
//
#include <AliLog.h>
#include "AliCFGrid.h"
#include "TMath.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"

//____________________________________________________________________
ClassImp(AliCFGrid)

//____________________________________________________________________
AliCFGrid::AliCFGrid() : 
  AliCFFrame(),
  fSumW2(kFALSE),
  fNunflTot(0),
  fNovflTot(0),
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
  AliCFFrame(name,title),
  fSumW2(kFALSE),
  fNunflTot(0),
  fNovflTot(0),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  // default constructor
}

//____________________________________________________________________
AliCFGrid::AliCFGrid(const Char_t* name, const Char_t* title, const Int_t nVarIn, const Int_t * nBinIn, const Float_t *binLimitsIn) :  
  AliCFFrame(name,title,nVarIn,nBinIn,binLimitsIn),
  fSumW2(kFALSE),
  fNunflTot(0),
  fNovflTot(0),
  fNentriesTot(0),
  fNunfl(0x0),
  fNovfl(0x0),
  fData(0x0),
  fErr2(0x0)
{
  //
  // main constructor
  //


  //The underflows
  fNunfl=new Float_t[fNVar];
  fNovfl= new Float_t[fNVar];


  // the grid
 
  fData = new Float_t[fNDim]; //num

  //Initialization
 
  for(Int_t i = 0;i<fNDim;i++)fData[i]=0;

  fNunflTot =0;
  fNovflTot =0;
  fNentriesTot =0;
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j] =0;
    fNovfl[j] =0;
  }
}

//____________________________________________________________________
AliCFGrid::AliCFGrid(const AliCFGrid& c) : 
  AliCFFrame(),
  fSumW2(kFALSE),
  fNunflTot(0),
  fNovflTot(0),
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
  delete fData;
  if(fSumW2)delete fErr2;
  delete fNunfl;
  delete fNovfl;

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
    Int_t ind =GetBinIndex(bin);
    return GetElement(ind);
}  
//____________________________________________________________________
Float_t AliCFGrid::GetElement(Float_t *var) const
{
  //
  // Get the content in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
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
  // Return the squared error on grid element iel
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, return -1"));
    return -1.;
  }
  if(fSumW2)return fErr2[iel];
  return fData[iel];
}
//____________________________________________________________________
Float_t AliCFGrid::GetElementError(Int_t *bin) const
{
  //
  // Get the squared error in a bin corresponding to a set of bin indeces
  //
    Int_t ind =GetBinIndex(bin);
    return GetElementError(ind);

}
//____________________________________________________________________
Float_t AliCFGrid::GetElementError(Float_t *var) const
{
  //
  // Get the squared error in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
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
    Int_t ind =GetBinIndex(bin);
    SetElement(ind,val);
}
//____________________________________________________________________
void AliCFGrid::SetElement(Float_t *var, Float_t val) 
{
  //
  // Set the content in a bin to value val corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
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
  // Set squared error to val on grid element iel
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, no value set"));
    return;
  }
  if(!fErr2)SumW2();
  fErr2[iel]=val;   
}
//____________________________________________________________________
void AliCFGrid::SetElementError(Int_t *bin, Float_t val) 
{
  //
  // Set squared error to val on grid element of bin indeces bin
  //
    Int_t ind =GetBinIndex(bin);
    SetElementError(ind,val);
}
//____________________________________________________________________
void AliCFGrid::SetElementError(Float_t *var, Float_t val)
{
  //
  // Set squared error to val in a bin corresponding to a set of input variables
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
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

  if(!(ovfl==1 || unfl==1)){
    SetElementError(fIndex,val);
  }
  else{
    AliInfo(Form(" input variables outside the grid, no value set"));
  }
} 
//____________________________________________________________________
void AliCFGrid::Scale(Int_t iel, Float_t *fact)
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  if(iel>=fNDim){
    AliInfo(Form(" element index outside the grid, no scaling"));
    return;
  }
  Float_t el,del,elsc,delsc;  
  if(GetElement(iel)>0 && fact[0]>0){
    el=GetElement(iel);
    del=TMath::Sqrt(GetElementError(iel));
    elsc=el*fact[0];
    delsc=TMath::Sqrt(del*del/el/el
		    +fact[1]*fact[1]/fact[0]/fact[0])
      *elsc;
    SetElement(iel,elsc);
    if(fSumW2)SetElementError(iel,delsc*elsc);
  }
}
//____________________________________________________________________
void AliCFGrid::Scale(Int_t *bin, Float_t *fact)
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  Int_t iel=GetBinIndex(bin);
  Scale(iel,fact);
}
//____________________________________________________________________
void AliCFGrid::Scale(Float_t *var, Float_t *fact) 
{
  //
  //scale content of a certain cell by (positive) fact (with error)
  //
  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
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

  if(!(ovfl==1 || unfl==1)){
    Int_t iel=GetBinIndex(fIndex);
    Scale(iel,fact);
  }
  else{
    AliInfo(Form(" input variables outside the grid, no scaling done"));
  }
}
//____________________________________________________________________
void AliCFGrid::Scale( Float_t *fact) 
{
  //
  //scale contents of the whole grid by fact
  //
  for(Int_t iel=0;iel<fNDim;iel++){
    Scale(iel,fact);
  }
}
//____________________________________________________________________
void AliCFGrid::Fill(Float_t *var, Float_t weight)
{

  //
  // Fill the grid,
  // given a set of values of the input variable, 
  // with weight (by default w=1)
  //

  Int_t unfl=0;  
  Int_t ovfl=0;  
  for(Int_t i=0;i<fNVar;i++){
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }

    fIndex[i] = TMath::BinarySearch(nbins,bins,var[i]);
    //underflows

    if(var[i] < bins[0]){
      unfl=1;
      fNunfl[i]++;
    }

    //overflows

    if(var[i] > bins[nbins-1]){
      ovfl=1;
      fNovfl[i]++;
    }
    delete [] bins;
  }

  // Total number of entries, overflows and underflows

  fNentriesTot++;
  if(unfl)fNunflTot++;
  if(ovfl)fNovflTot++;

  //if not ovfl/unfl, fill the element  
  if(!(ovfl==1 || unfl==1)){
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
Float_t AliCFGrid::GetOverFlows() const {
  //
  // Get overflows 
  //
  return fNovflTot;
} 
//____________________________________________________________________
Float_t AliCFGrid::GetUnderFlows( Int_t ivar) const {
  //
  // Get overflows in variable var 
  //
  return fNunfl[ivar];
} 
//____________________________________________________________________
Float_t AliCFGrid::GetUnderFlows() const {
  //
  // Get overflows 
  //
  return fNunflTot;
} 
//____________________________________________________________________
Float_t AliCFGrid::GetEntries( ) const {
  //
  // Get total entries (in grid + over/underflows) 
  //
  return fNentriesTot;
} 
//____________________________________________________________________
Int_t AliCFGrid::GetEmptyBins() const {
  //
  // Get empty bins 
  //
  Int_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    if(fData[i]<=0)val++;     
  }
  return val;
} 
//_____________________________________________________________________
Int_t AliCFGrid::GetEmptyBins( Float_t *varMin, Float_t* varMax ) const 
{
  //
  // Get empty bins in a range
  //

  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];

  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Float_t xmin=varMin[i]; // the min values  
    Float_t xmax=varMax[i]; // the min values  
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
    indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
    if(xmax>=bins[nbins-1]){
      indexMax[i]=indexMax[i]-1;
    }  
    delete [] bins;
  }

  Int_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    for (Int_t j=0;j<fNVar;j++)fIndex[j]=GetBinIndex(j,i);
    Bool_t isIn=kTRUE;
    for (Int_t j=0;j<fNVar;j++){
      if(!(fIndex[j]>=indexMin[j] && fIndex[j]<=indexMax[j]))isIn=kFALSE;   
    }
    if(isIn && fData[i]<=0)val++;     
  }
  AliInfo(Form(" the empty bins = %i ",val)); 

  delete [] indexMin;
  delete [] indexMax;
  return val;
} 
//____________________________________________________________________
Int_t AliCFGrid::CheckEfficiencyStats(Float_t thr) const
{
  //
  // Count the cells below a certain threshold
  //
  Int_t ncellsLow=0;
  for(Int_t i=0;i<fNDim;i++){
    if(GetElement(i)<thr)ncellsLow++;
  }
  return ncellsLow;
}
//_____________________________________________________________________
Float_t AliCFGrid::GetIntegral() const 
{
  //
  // Get full Integral
  //
  Float_t val=0;
  for(Int_t i=0;i<fNDim;i++){
    val+=fData[i];     
  }
  return val;  
} 
//_____________________________________________________________________
Float_t AliCFGrid::GetIntegral(Int_t *binMin, Int_t* binMax ) const 
{
  //
  // Get Integral in a range of bin indeces (extremes included)
  //

  Float_t val=0;
  for(Int_t i=0;i<fNVar;i++){
    if((binMin[i]<0) || (binMax[i]>=fNVarBins[i]) || (binMin[i]>binMax[i])){
      AliInfo(Form(" Bin indeces in variable %i outside allowed range or in reverse order, please check!", i));
      return val;
    }
  }
  val=GetSum(0,binMin,binMax);
  return val;
} 
//_____________________________________________________________________
Float_t AliCFGrid::GetIntegral(Float_t *varMin, Float_t* varMax ) const 
{
  //
  // Get Integral in a range (extremes included)
  //

  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];

  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Float_t xmin=varMin[i]; // the min values  
    Float_t xmax=varMax[i]; // the min values  
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
    indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
    if(xmax>=bins[nbins-1]){
      indexMax[i]=indexMax[i]-1;
    }  
    delete [] bins;
  }

  Float_t val=GetIntegral(indexMin,indexMax);

  delete [] indexMin;
  delete [] indexMax;

  return val;
} 
//___________________________________________________________________
TH1F *AliCFGrid::Project(Int_t ivar) const
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

  TH1F *proj1D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH1F")) {
    proj1D = (TH1F*)obj;
    proj1D->Reset();
  }

  if(!proj1D){
    proj1D =new TH1F(pname,htitle, nbins, bins);
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
  proj1D->SetEntries(GetEntries());
  return proj1D;
} 

//___________________________________________________________________
TH2F *AliCFGrid::Project(Int_t ivar1, Int_t ivar2) const
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

  TH2F *proj2D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH2F")) {
    proj2D = (TH2F*)obj;
    proj2D->Reset();
  }

  if(!proj2D){
    proj2D =new TH2F(pname,htitle, nbins1, bins1,nbins2,bins2);
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
  proj2D->SetEntries(GetEntries());
  return proj2D;
} 
//___________________________________________________________________
TH3F *AliCFGrid::Project(Int_t ivar1, Int_t ivar2, Int_t ivar3) const
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

  TH3F *proj3D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH3F")) {
    proj3D = (TH3F*)obj;
    proj3D->Reset();
  }

  if(!proj3D){
    proj3D =new TH3F(pname,htitle, nbins1,bins1,nbins2,bins2,nbins3,bins3);
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
  proj3D->SetEntries(GetEntries());
  return proj3D;
} 

//___________________________________________________________________
TH1F *AliCFGrid::Slice(Int_t ivar, Float_t *varMin, Float_t* varMax) const
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

  TH1F *proj1D=0;

  //check if a projection with identical name exist
  TObject *obj = gROOT->FindObject(pname);
  if (obj && obj->InheritsFrom("TH1F")) {
    proj1D = (TH1F*)obj;
    proj1D->Reset();
  }

  if(!proj1D){
    proj1D =new TH1F(pname,htitle, nbins, bins);
  }  

  delete [] bins;


  Int_t *indexMin=new Int_t[fNVar];
  Int_t *indexMax=new Int_t[fNVar];


  //Find out the min and max bins

  for(Int_t i=0;i<fNVar;i++){
    Float_t xmin=varMin[i]; // the min values  
    Float_t xmax=varMax[i]; // the min values  
    Int_t nbins=fNVarBins[i]+1;
    Float_t *bins=new Float_t[nbins];
    for(Int_t ibin =0;ibin<nbins;ibin++){
     bins[ibin] = fVarBinLimits[ibin+fOffset[i]];
    }
    indexMin[i] = TMath::BinarySearch(nbins,bins,xmin);
    indexMax[i] = TMath::BinarySearch(nbins,bins,xmax);
    if(xmax>=bins[nbins-1]){
      indexMax[i]=indexMax[i]-1;
    }  
    delete [] bins;
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
Long64_t AliCFGrid::Merge(TCollection* list)
{
  // Merge a list of AliCorrection objects with this (needed for
  // PROOF). 
  // Returns the number of merged objects (including this).

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  TIterator* iter = list->MakeIterator();
  TObject* obj;
  
  Int_t count = 0;
  while ((obj = iter->Next())) {
    AliCFGrid* entry = dynamic_cast<AliCFGrid*> (obj);
    if (entry == 0) 
      continue;
    this->Add(entry);
    count++;
  }

  return count+1;
}

//____________________________________________________________________
void AliCFGrid::Add(AliCFGrid* aGrid, Float_t c)
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
      Float_t err=TMath::Sqrt(aGrid->GetElementError(iel));  
      fErr2[iel]+=c*c*err*err;
    }
  }

  //Add entries, overflows and underflows

  fNentriesTot+= c*aGrid->GetEntries();
  fNunflTot+= c*aGrid->GetUnderFlows();
  fNovflTot+= c*aGrid->GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]+= c*aGrid->GetUnderFlows(j);
    fNovfl[j]+= c*aGrid->GetUnderFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Add(AliCFGrid* aGrid1, AliCFGrid* aGrid2, Float_t c1,Float_t c2)
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
      err1=TMath::Sqrt(aGrid1->GetElementError(iel));
      err2=TMath::Sqrt(aGrid2->GetElementError(iel));
      SetElementError(iel,c1*c1*err1*err1+c2*c2*err2*err2);
    }
  }

  //Add entries, overflows and underflows

  fNentriesTot= c1*aGrid1->GetEntries()+c2*aGrid2->GetEntries();
  fNunflTot= c1*aGrid1->GetUnderFlows()+c2*aGrid2->GetUnderFlows();
  fNovflTot= c1*aGrid1->GetOverFlows()+c2*aGrid2->GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= c1*aGrid1->GetUnderFlows(j)+c2*aGrid2->GetUnderFlows(j);
    fNovfl[j]= c1*aGrid1->GetUnderFlows(j)+c2*aGrid2->GetUnderFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Multiply(AliCFGrid* aGrid, Float_t c)
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
      err1=TMath::Sqrt(GetElementError(iel));  
      err2=TMath::Sqrt(aGrid->GetElementError(iel));  
      SetElementError(iel,c*c*(cont2*cont2*err1*err1+cont1*cont1*err2*err2));
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  fNunflTot=GetUnderFlows();
  fNovflTot=GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetUnderFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Multiply(AliCFGrid* aGrid1,AliCFGrid* aGrid2, Float_t c1,Float_t c2)
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
      err1=TMath::Sqrt(aGrid1->GetElementError(iel));  
      err2=TMath::Sqrt(aGrid2->GetElementError(iel));  
      SetElementError(iel,c1*c1*c2*c2*(cont2*cont2*err1*err1+cont1*cont1*err2*err2));
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  fNunflTot=GetUnderFlows();
  fNovflTot=GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetUnderFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Divide(AliCFGrid* aGrid, Float_t c, Option_t *option)
{
  //
  //divide current grid by grid aGrid
  //

  TString opt = option;
  opt.ToUpper();

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

  Float_t cont1,cont2,err1,err2,r,den;
  for(Int_t iel=0;iel<fNDim;iel++){
    cont1=GetElement(iel);  
    cont2=aGrid->GetElement(iel);
    if(cont2)SetElement(iel,cont1/(c*cont2));
    else SetElement(iel,0);
    if(fSumW2){
      err1=TMath::Sqrt(GetElementError(iel));  
      err2=TMath::Sqrt(aGrid->GetElementError(iel));  
      if(!cont2){SetElementError(iel,0.); continue;}
      if (opt.Contains("B")){
	if(cont1!=cont2){
	  r=cont1/cont2;	    
	  SetElementError(iel,TMath::Abs(((1-2.*r)*err1*err1+r*r*err2*err2)/(cont2*cont2)));
	}else{
	  SetElementError(iel,0.);
	}
      }else{
        den=cont2*cont2*cont2*c*c;
	SetElementError(iel,(cont2*cont2*err1*err1+cont1*cont1*err2*err2)/den);
      }
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  fNunflTot=GetUnderFlows();
  fNovflTot=GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetUnderFlows(j);
  }
}
//____________________________________________________________________
void AliCFGrid::Divide(AliCFGrid* aGrid1, AliCFGrid* aGrid2, Float_t c1,Float_t c2, Option_t *option)
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
      err1=TMath::Sqrt(aGrid1->GetElementError(iel));  
      err2=TMath::Sqrt(aGrid2->GetElementError(iel));  
      if(!cont2){SetElementError(iel,0.); continue;}
      if (opt.Contains("B")){
	if(cont1!=cont2){
	  r=cont1/cont2;	    
	  SetElementError(iel,TMath::Abs(((1.-2.*r)*err1*err1+r*r*err2*err2)/(cont2*cont2)));
	}else{
	  SetElementError(iel,0.);
	}
      }else{
        den=cont2*cont2*cont2*cont2*c2*c2;
	SetElementError(iel,c1*c1*(cont2*cont2*err1*err1+cont1*cont1*err2*err2)/den);
      }
    }
  }

  //Set entries to the number of bins, preserve original overflows and underflows

  fNentriesTot=fNDim;
  fNunflTot=GetUnderFlows();
  fNovflTot=GetOverFlows();
  for(Int_t j=0;j<fNVar;j++){
    fNunfl[j]= GetUnderFlows(j);
    fNovfl[j]= GetUnderFlows(j);
  }
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
//_____________________________________________________________________
Float_t AliCFGrid::GetSum(Int_t ivar, Int_t *binMin, Int_t* binMax) const 
{
  //
  // recursively add over nested loops.... 
  //
  static Float_t val;
  if(ivar==0)val=0.;
  for(Int_t ibin=binMin[ivar];ibin<=binMax[ivar];ibin++){
    fIndex[ivar]=ibin;
    if(ivar<fNVar-1) {
      val=GetSum(ivar+1,binMin,binMax);
    }
    else {
      Int_t iel=GetBinIndex(fIndex);
      val+=fData[iel];
    }
  }

  return val;
}
//____________________________________________________________________
void AliCFGrid::Copy(TObject& c) const
{
  //
  // copy function
  //
  AliCFGrid& target = (AliCFGrid &) c;

  target.fNVar=fNVar;
  target.fNDim=fNDim;
  target.fSumW2=fSumW2;
  target.fNVarBinLimits=fNVarBinLimits;
  target.fNunflTot = fNunflTot;
  target.fNovflTot = fNovflTot;
  target.fNentriesTot = fNentriesTot;
  if (fNVarBins)
    target.fNVarBins = fNVarBins;
  if (fVarBinLimits)
    target.fVarBinLimits = fVarBinLimits;
  if (fNunfl)
    target.fNunfl = fNunfl;
  if (fNunfl)
    target.fNunfl = fNunfl;
  if (fNovfl)
    target.fNovfl = fNovfl;
  if (fProduct)
    target.fProduct = fProduct;
  if (fOffset)
    target.fOffset = fOffset;
  if (fData)
    target.fData = fData;
  if (fErr2)
    target.fErr2 = fErr2;
  
}
