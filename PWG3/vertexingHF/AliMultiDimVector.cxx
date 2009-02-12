/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id: $ */

///////////////////////////////////////////////////////////////////
//                                                               //
// Implementation of the class to store the number of signal     //
// and background events in bins of the cut values               //
// Origin:       Elena Bruna (bruna@to.infn.it)                  //
// Updated:      Sergey Senyukov (senyukov@to.infn.it)           //
// Last updated: Francesco Prino (prino@to.infn.it)              //
//                                                               //
///////////////////////////////////////////////////////////////////


#include "TH2.h"
#include "AliMultiDimVector.h"
#include "AliLog.h"
#include "TMath.h"
#include "TString.h"

ClassImp(AliMultiDimVector)
//___________________________________________________________________________
AliMultiDimVector::AliMultiDimVector():TNamed("AliMultiDimVector","default"),
fNVariables(0),
fNPtBins(0),
fVett(0),
fNTotCells(0),
fIsIntegrated(0)
{
  // default constructor
}
//___________________________________________________________________________
AliMultiDimVector::AliMultiDimVector(const char *name,const char *title, const Int_t npars, Int_t nptbins, Int_t *nofcells, Float_t *loosecuts, Float_t *tightcuts, TString *axisTitles):TNamed(name,title),
fNVariables(npars),
fNPtBins(nptbins),
fVett(0),
fNTotCells(0),
fIsIntegrated(0){
// standard constructor
  ULong64_t ntot=1; 
  for(Int_t i=0;i<fNVariables;i++){
    ntot*=nofcells[i];
    fNCutSteps[i]=nofcells[i];
    if(loosecuts[i] <= tightcuts[i]){
      fMinLimits[i]=loosecuts[i];
      fMaxLimits[i]=tightcuts[i];
      fGreaterThan[i]=kTRUE;
    }else{
      fMinLimits[i]=tightcuts[i];
      fMaxLimits[i]=loosecuts[i];
      fGreaterThan[i]=kFALSE;
    }
    fAxisTitles[i]=axisTitles[i].Data();
  }
  fNTotCells=ntot*fNPtBins;
  fVett.Set(fNTotCells); 
  for (ULong64_t j=0;j<fNTotCells;j++) fVett.AddAt(0,j);
}
//___________________________________________________________________________
AliMultiDimVector::AliMultiDimVector(const AliMultiDimVector &mv):TNamed(mv.GetName(),mv.GetTitle()),
fNVariables(mv.fNVariables),
fNPtBins(mv.fNPtBins),
fVett(0),
fNTotCells(mv.fNTotCells),
fIsIntegrated(mv.fIsIntegrated)
{
  // copy constructor
  for(Int_t i=0;i<fNVariables;i++){
    fNCutSteps[i]=mv.GetNCutSteps(i);
    fMinLimits[i]=mv.GetMinLimit(i);
    fMaxLimits[i]=mv.GetMaxLimit(i);
    fGreaterThan[i]=mv.GetGreaterThan(i);
    fAxisTitles[i]=mv.GetAxisTitle(i);
  }
  fVett.Set(fNTotCells); 
  for(ULong64_t i=0;i<fNTotCells;i++) fVett[i]=mv.GetElement(i);
}
//___________________________________________________________________________
void AliMultiDimVector::CopyStructure(const AliMultiDimVector* mv){
// Sets dimensions and limit from mv
  fNVariables=mv->GetNVariables();
  fNPtBins=mv->GetNPtBins();
  fNTotCells=mv->GetNTotCells();
  fIsIntegrated=mv->IsIntegrated();
  for(Int_t i=0;i<fNVariables;i++){
    fNCutSteps[i]=mv->GetNCutSteps(i);
    fMinLimits[i]=mv->GetMinLimit(i);
    fMaxLimits[i]=mv->GetMaxLimit(i);
    fGreaterThan[i]=mv->GetGreaterThan(i);
    fAxisTitles[i]=mv->GetAxisTitle(i);
  }
  fVett.Set(fNTotCells);
  
}
//______________________________________________________________________
Bool_t AliMultiDimVector::GetIndicesFromGlobalAddress(ULong64_t globadd, Int_t *ind, Int_t &ptbin) const {
// returns matrix element indices and Pt bin from global index
  if(globadd>=fNTotCells) return kFALSE;
  ULong64_t r=globadd;
  Int_t prod=1;
  Int_t nOfCellsPlusLevel[fgkMaxNVariables+1];
  for(Int_t k=0;k<fNVariables;k++) nOfCellsPlusLevel[k]=fNCutSteps[k];
  nOfCellsPlusLevel[fNVariables]=fNPtBins;
	
  for(Int_t i=0;i<fNVariables+1;i++) prod*=nOfCellsPlusLevel[i];
  for(Int_t i=0;i<fNVariables+1;i++){
    prod/=nOfCellsPlusLevel[i];
    if(i<fNVariables) ind[i]=r/prod;
    else ptbin=r/prod;
    r=globadd%prod;
  }
  return kTRUE;
}
//______________________________________________________________________
Bool_t AliMultiDimVector::GetCutValuesFromGlobalAddress(ULong64_t globadd, Float_t *cuts, Int_t &ptbin) const {
  Int_t ind[fgkMaxNVariables];
  Bool_t retcode=GetIndicesFromGlobalAddress(globadd,ind,ptbin);
  if(!retcode) return kFALSE;
  for(Int_t i=0;i<fNVariables;i++) cuts[i]=GetCutValue(i,ind[i]);
  return kTRUE;
}
//______________________________________________________________________
ULong64_t AliMultiDimVector::GetGlobalAddressFromIndices(Int_t *ind, Int_t ptbin) const {
  // Returns the global index of the cell in the matrix
  Int_t prod=1;
  ULong64_t elem=0;
  Int_t indexPlusLevel[fgkMaxNVariables+1];
  Int_t nOfCellsPlusLevel[fgkMaxNVariables+1];
  for(Int_t i=0;i<fNVariables;i++){
    indexPlusLevel[i]=ind[i];
    nOfCellsPlusLevel[i]=fNCutSteps[i];
  }
  indexPlusLevel[fNVariables]=ptbin;
  nOfCellsPlusLevel[fNVariables]=fNPtBins;
	
  for(Int_t i=0;i<fNVariables+1;i++){
    prod=indexPlusLevel[i];
    if(i<fNVariables){
      for(Int_t j=i+1;j<fNVariables+1;j++){
	prod*=nOfCellsPlusLevel[j];
      }
    }
    elem+=prod;
  }
  return elem;
}
//______________________________________________________________________
Bool_t AliMultiDimVector::GetIndicesFromValues(Float_t *values, Int_t *ind) const {
  // Fills the array of matrix indices strating from variable values
  for(Int_t i=0;i<fNVariables;i++){
    if(fGreaterThan[i]){ 
      if(values[i]<GetMinLimit(i)) return kFALSE;
      ind[i]=(Int_t)((values[i]-fMinLimits[i])/GetCutStep(i));
      if(ind[i]>=GetNCutSteps(i)) ind[i]=GetNCutSteps(i)-1;
    }else{
      if(values[i]>GetMaxLimit(i)) return kFALSE;
      ind[i]=(Int_t)((fMaxLimits[i]-values[i])/GetCutStep(i));
      if(ind[i]>=GetNCutSteps(i)) ind[i]=GetNCutSteps(i)-1;
    }
  }
  return kTRUE;
}
//______________________________________________________________________
ULong64_t AliMultiDimVector::GetGlobalAddressFromValues(Float_t *values, Int_t ptbin) const {
  // Returns the global index of the cell in the matrix
   Int_t ind[fgkMaxNVariables];
   Bool_t retcode=GetIndicesFromValues(values,ind);
   if(retcode) return GetGlobalAddressFromIndices(ind,ptbin);
   else{
     AliError("Values out of range");
     return fNTotCells+999;
   }
}
//_____________________________________________________________________________
void AliMultiDimVector::MultiplyBy(Float_t factor){
  // multiply the AliMultiDimVector by a constant factor
  for(ULong64_t i=0;i<fNTotCells;i++){
    if(fVett.At(i)!=-1)
      fVett.AddAt(fVett.At(i)*factor,i);
    else fVett.AddAt(-1,i);
  }
  
}
//_____________________________________________________________________________
void AliMultiDimVector::Multiply(AliMultiDimVector* mv,Float_t factor){
  //  Sets AliMultiDimVector=mv*constant factor
  for(ULong64_t i=0;i<fNTotCells;i++){
    if(mv->GetElement(i)!=-1)
      fVett.AddAt(mv->GetElement(i)*factor,i);
    else fVett.AddAt(-1,i); 
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::Multiply(AliMultiDimVector* mv1,AliMultiDimVector* mv2){
  //  Sets AliMultiDimVector=mv1*mv2
  for(ULong64_t i=0;i<fNTotCells;i++){
    if(mv1->GetElement(i)!=-1 && mv2->GetElement(i)!=-1)
      fVett.AddAt(mv1->GetElement(i)*mv2->GetElement(i),i);
    else fVett.AddAt(-1,i); 
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::Add(AliMultiDimVector* mv){
  // Sums contents of mv to AliMultiDimVector
  if (mv->GetNTotCells()!=fNTotCells){ 
    AliError("Different dimension of the vectors!!");
  }else{
    for(ULong64_t i=0;i<fNTotCells;i++) 
      if(mv->GetElement(i)!=-1 && fVett.At(i)!=-1)
	fVett.AddAt(fVett.At(i)+mv->GetElement(i),i);
      else fVett.AddAt(-1,i); 
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::Sum(AliMultiDimVector* mv1,AliMultiDimVector* mv2){
  // Sets AliMultiDimVector=mv1+mv2
  if (fNTotCells!=mv1->GetNTotCells()&&mv1->GetNTotCells()!=mv2->GetNTotCells()) {
    AliError("Different dimension of the vectors!!");
  }
  else{
    for(ULong64_t i=0;i<mv1->GetNTotCells();i++) {
      if(mv1->GetElement(i)!=-1 && mv2->GetElement(i)!=-1)
	fVett.AddAt(mv1->GetElement(i)+mv2->GetElement(i),i); 
      else fVett.AddAt(-1,i); 
    }
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::LinearComb(AliMultiDimVector* mv1, Float_t norm1, AliMultiDimVector* mv2, Float_t norm2){
  // Sets AliMultiDimVector=n1*mv1+n2*mv2
  if (fNTotCells!=mv1->GetNTotCells()&&mv1->GetNTotCells()!=mv2->GetNTotCells()) {
    AliError("Different dimension of the vectors!!");
  }
  else{
    for(ULong64_t i=0;i<mv1->GetNTotCells();i++) {
      if(mv1->GetElement(i)!=-1 && mv2->GetElement(i)!=-1)
	fVett.AddAt(norm1*mv1->GetElement(i)+norm2*mv2->GetElement(i),i); 
      else fVett.AddAt(-1,i); 
    }
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::DivideBy(AliMultiDimVector* mv){
  // Divide AliMulivector by mv
  if (mv->GetNTotCells()!=fNTotCells) {
    AliError("Different dimension of the vectors!!");
  }
  else{
    for(ULong64_t i=0;i<fNTotCells;i++) 
      if(mv->GetElement(i)!=0 &&mv->GetElement(i)!=-1 && fVett.At(i)!=-1)
	fVett.AddAt(fVett.At(i)/mv->GetElement(i),i);
      else fVett.AddAt(-1,i);
  }

}
//_____________________________________________________________________________
void AliMultiDimVector::Divide(AliMultiDimVector* mv1,AliMultiDimVector* mv2){
  // Sets AliMultiDimVector=mv1/mv2
  if (fNTotCells!=mv1->GetNTotCells()&&mv1->GetNTotCells()!=mv2->GetNTotCells()) {
    AliError("Different dimension of the vectors!!");
  }
  else{
    for(ULong64_t i=0;i<mv1->GetNTotCells();i++) 
      if(mv2->GetElement(i)!=0&& mv2->GetElement(i)!=-1&& mv1->GetElement(i)!=-1)
	{
	  fVett.AddAt(mv1->GetElement(i)/mv2->GetElement(i),i);
	}
      else fVett.AddAt(-1,i);
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::Sqrt(){
  // Sqrt of elements of AliMultiDimVector
  for(ULong64_t i=0;i<fNTotCells;i++) {
    if(fVett.At(i)>=0) fVett.AddAt(TMath::Sqrt(fVett.At(i)),i);
    else {
      fVett.AddAt(-1,i);
    }
  }
}
//_____________________________________________________________________________
void AliMultiDimVector::Sqrt(AliMultiDimVector* mv){
  // Sets AliMultiDimVector=sqrt(mv)
  for(ULong64_t i=0;i<fNTotCells;i++) 
    if(mv->GetElement(i)>=0) fVett.AddAt(TMath::Sqrt(mv->GetElement(i)),i);
    else fVett.AddAt(-1,i);
}
//_____________________________________________________________________________
void AliMultiDimVector::FindMaximum(Float_t& maxValue, Int_t *ind , Int_t ptbin){
  // finds the element with maximum contents
  const ULong64_t nelem=fNTotCells/fNPtBins;
  TArrayF vett;
  vett.Set(nelem);
  ULong64_t runningAddress;
  for(ULong64_t i=0;i<nelem;i++){
    runningAddress=ptbin+i*fNPtBins;
    vett.AddAt(fVett[runningAddress],i);
  }
  maxValue=TMath::MaxElement(nelem,vett.GetArray());
  ULong64_t maxAddress=TMath::LocMax(nelem,vett.GetArray());
  ULong64_t maxGlobalAddress=ptbin+maxAddress*fNPtBins;
  Int_t checkedptbin;
  GetIndicesFromGlobalAddress(maxGlobalAddress,ind,checkedptbin);
}

//_____________________________________________________________________________
TH2F*  AliMultiDimVector::Project(Int_t firstVar, Int_t secondVar, Int_t* fixedVars, Int_t ptbin, Float_t norm){
  // Project the AliMultiDimVector on a 2D histogram

  TString hisName=Form("hproj%s%dv%d",GetName(),secondVar,firstVar);
  TString hisTit=Form("%s vs. %s",fAxisTitles[secondVar].Data(),fAxisTitles[firstVar].Data());
  TH2F* h2=new TH2F(hisName.Data(),hisTit.Data(),fNCutSteps[firstVar],fMinLimits[firstVar],fMaxLimits[firstVar],fNCutSteps[secondVar],fMinLimits[secondVar],fMaxLimits[secondVar]);
	
  Int_t index[fgkMaxNVariables];
  for(Int_t i=0;i<fNVariables;i++){
    index[i]=fixedVars[i];
  }
  
  for(Int_t i=0;i<fNCutSteps[firstVar];i++){
    for(Int_t j=0;j<fNCutSteps[secondVar];j++){
      index[firstVar]=i;
      index[secondVar]=j;
      Float_t cont=GetElement(index,ptbin)/norm;
      Int_t bin1=i+1;
      if(!fGreaterThan[firstVar]) bin1=fNCutSteps[firstVar]-i;
      Int_t bin2=j+1;
      if(!fGreaterThan[secondVar]) bin2=fNCutSteps[secondVar]-j;
      h2->SetBinContent(bin1,bin2,cont);
    }
  }
  return h2;
}
//_____________________________________________________________________________ 
void  AliMultiDimVector::GetIntegrationLimits(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const {
  // computes bin limits for integrating the AliMultiDimVector
  minbin=0;
  maxbin=0;
  if(iVar<fNVariables){
    minbin=iCell;
    maxbin=fNCutSteps[iVar]-1;
  }
}
//_____________________________________________________________________________ 
void  AliMultiDimVector::GetFillRange(Int_t iVar, Int_t iCell, Int_t& minbin, Int_t& maxbin) const {
  // computes range of cells passing the cuts for FillAndIntegrate
  minbin=0;
  maxbin=0;
  if(iVar<fNVariables){
    minbin=0; // bin 0 corresponds to loose cuts
    maxbin=iCell;
  }
}
//_____________________________________________________________________________ 
void AliMultiDimVector::Integrate(){
  // integrates the matrix
  if(fIsIntegrated){
    AliError("MultiDimVector already integrated");
    return;
  }
  TArrayF integral(fNTotCells);
  for(ULong64_t i=0;i<fNTotCells;i++) integral[i]=CountsAboveCell(i);
  for(ULong64_t i=0;i<fNTotCells;i++) fVett[i]= integral[i];
  fIsIntegrated=kTRUE;
}
//_____________________________________________________________________________ 
Float_t AliMultiDimVector::CountsAboveCell(ULong64_t globadd) const{
  // integrates the counts of cells above cell with address globadd
  Int_t ind[fgkMaxNVariables];
  Int_t ptbin;
  GetIndicesFromGlobalAddress(globadd,ind,ptbin);
  for(Int_t i=fNVariables; i<fgkMaxNVariables; i++) ind[i]=0;
  Int_t mink[fgkMaxNVariables];
  Int_t maxk[fgkMaxNVariables];
  for(Int_t i=0;i<fgkMaxNVariables;i++){
    GetIntegrationLimits(i,ind[i],mink[i],maxk[i]);
  }
  Float_t sumcont=0.;
  for(Int_t k0=mink[0]; k0<=maxk[0]; k0++){
    for(Int_t k1=mink[1]; k1<=maxk[1]; k1++){
      for(Int_t k2=mink[2]; k2<=maxk[2]; k2++){
	for(Int_t k3=mink[3]; k3<=maxk[3]; k3++){
	  for(Int_t k4=mink[4]; k4<=maxk[4]; k4++){
	    for(Int_t k5=mink[5]; k5<=maxk[5]; k5++){
	      for(Int_t k6=mink[6]; k6<=maxk[6]; k6++){
		for(Int_t k7=mink[7]; k7<=maxk[7]; k7++){
		  for(Int_t k8=mink[8]; k8<=maxk[8]; k8++){
		    for(Int_t k9=mink[9]; k9<=maxk[9]; k9++){
		      Int_t currentBin[fgkMaxNVariables]={k0,k1,k2,k3,k4,k5,k6,k7,k8,k9};
		      sumcont+=GetElement(currentBin,ptbin);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return sumcont;
}
//_____________________________________________________________________________ 
void AliMultiDimVector::Fill(Float_t* values, Int_t ptbin){
  // fills the cells of AliMultiDimVector corresponding to values
  if(fIsIntegrated){
    AliError("MultiDimVector already integrated -- Use FillAndIntegrate");
    return;
  }
  Int_t ind[fgkMaxNVariables];
  Bool_t retcode=GetIndicesFromValues(values,ind);
  for(Int_t i=fNVariables; i<fgkMaxNVariables; i++) ind[i]=0;
  if(retcode) IncrementElement(ind,ptbin);
}
//_____________________________________________________________________________ 
void AliMultiDimVector::FillAndIntegrate(Float_t* values, Int_t ptbin){
  // fills the cells of AliMultiDimVector passing the cuts
  // The number of nested loops must match fgkMaxNVariables!!!!
  fIsIntegrated=kTRUE;
  Int_t ind[fgkMaxNVariables];
  Bool_t retcode=GetIndicesFromValues(values,ind);
  if(!retcode) return;
  for(Int_t i=fNVariables; i<fgkMaxNVariables; i++) ind[i]=0;
  Int_t mink[fgkMaxNVariables];
  Int_t maxk[fgkMaxNVariables];
  for(Int_t i=0;i<fgkMaxNVariables;i++){
    GetFillRange(i,ind[i],mink[i],maxk[i]);
  }
  for(Int_t k0=mink[0]; k0<=maxk[0]; k0++){
    for(Int_t k1=mink[1]; k1<=maxk[1]; k1++){
      for(Int_t k2=mink[2]; k2<=maxk[2]; k2++){
	for(Int_t k3=mink[3]; k3<=maxk[3]; k3++){
	  for(Int_t k4=mink[4]; k4<=maxk[4]; k4++){
	    for(Int_t k5=mink[5]; k5<=maxk[5]; k5++){
	      for(Int_t k6=mink[6]; k6<=maxk[6]; k6++){
		for(Int_t k7=mink[7]; k7<=maxk[7]; k7++){
		  for(Int_t k8=mink[8]; k8<=maxk[8]; k8++){
		    for(Int_t k9=mink[9]; k9<=maxk[9]; k9++){
		      Int_t currentBin[fgkMaxNVariables]={k0,k1,k2,k3,k4,k5,k6,k7,k8,k9};
		      IncrementElement(currentBin,ptbin);
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

}
//_____________________________________________________________________________ 
void AliMultiDimVector::SuppressZeroBKGEffect(AliMultiDimVector* mvBKG){
  // Sets to zero elements for which mvBKG=0
  for(ULong64_t i=0;i<fNTotCells;i++)
    if(mvBKG->GetElement(i)==0) fVett.AddAt(0,i);
}
//_____________________________________________________________________________ 
AliMultiDimVector* AliMultiDimVector:: ShrinkPtBins(Int_t firstBin, Int_t lastBin){
  // sums the elements of pt bins between firstBin and lastBin
  if(firstBin<0 || lastBin>=fNPtBins || firstBin>=lastBin){
    AliError("Bad numbers of Pt bins to be shrinked");
    return 0;
  }
  Int_t nofcells[fgkMaxNVariables];
  Float_t loosecuts[fgkMaxNVariables];
  Float_t tightcuts[fgkMaxNVariables];
  TString axisTitles[fgkMaxNVariables];
  for(Int_t i=0;i<fNVariables;i++){
    nofcells[i]=fNCutSteps[i];
    if(fGreaterThan[i]){
      loosecuts[i]=fMinLimits[i];
      tightcuts[i]=fMaxLimits[i];
    }else{
      loosecuts[i]=fMaxLimits[i];
      tightcuts[i]=fMinLimits[i];
    }
    axisTitles[i]=fAxisTitles[i];
  }
  Int_t newNptbins=fNPtBins-(lastBin-firstBin);
  AliMultiDimVector* shrinkedMV=new AliMultiDimVector(GetName(),GetTitle(),fNVariables,newNptbins,nofcells,loosecuts,tightcuts,axisTitles);
  
  ULong64_t nOfPointsPerPtbin=fNTotCells/fNPtBins;
  ULong64_t addressOld,addressNew;
  Int_t npb,opb;
  for(npb=0;npb<firstBin;npb++){
    opb=npb;
    for(ULong64_t k=0;k<nOfPointsPerPtbin;k++){
      addressOld=opb+k*fNPtBins;
      addressNew=npb+k*newNptbins;
      shrinkedMV->SetElement(addressNew,fVett[addressOld]);
    }
  }
  npb=firstBin;
  for(ULong64_t k=0;k<nOfPointsPerPtbin;k++){
    Float_t summedValue=0.;
    for(opb=firstBin;opb<=lastBin;opb++){
      addressOld=opb+k*fNPtBins;
      summedValue+=fVett[addressOld];
    }
    addressNew=npb+k*newNptbins;
    shrinkedMV->SetElement(addressNew,summedValue);
  }
  for(npb=firstBin+1;npb<newNptbins;npb++){
    opb=npb+(lastBin-firstBin);
    for(ULong64_t k=0;k<nOfPointsPerPtbin;k++){
      addressOld=opb+k*fNPtBins;
      addressNew=npb+k*newNptbins;
      shrinkedMV->SetElement(addressNew,fVett[addressOld]);
    }
  }
  return shrinkedMV;
}
