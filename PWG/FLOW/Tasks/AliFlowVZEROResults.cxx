#include<stdio.h>
#include "AliFlowVZEROResults.h"
#include "TList.h"

ClassImp(AliFlowVZEROResults);

AliFlowVZEROResults::AliFlowVZEROResults(const char *name,const Int_t nvar,const Int_t* binVar) :
  TNamed(name,name),
  fNbinVar(new TArrayI(nvar)),
  fXmin(new TArrayF(nvar)),
  fXmax(new TArrayF(nvar)),
  fNameVar(new TClonesArray("TNamed")),
  fV2(new TClonesArray("TProfile"))
{

  for(Int_t i=0;i < GetNvar();i++){
    (*fNbinVar)[i] = binVar[i];
  }
  
  for(Int_t i=0; i<GetNvar();i++){
    new((*fNameVar)[i]) TNamed("","");
  }
}

AliFlowVZEROResults::AliFlowVZEROResults() :
  TNamed("v2","v2"),
  fNbinVar(new TArrayI(0)),
  fXmin(new TArrayF(0)),
  fXmax(new TArrayF(0)),
  fNameVar(new TClonesArray("TNamed")),
  fV2(new TClonesArray("TProfile"))
{
  
}

AliFlowVZEROResults::~AliFlowVZEROResults(){
  for(Int_t i=fNameVar->GetEntries();i>0;i--){
    delete fNameVar->At(i-1);
    fNameVar->RemoveAt(i-1);
  }

  for(Int_t i=fV2->GetEntries();i>0;i--){
    delete fV2->At(i-1);
    fV2->RemoveAt(i-1);
  }

  delete fNbinVar;
  delete fXmin;
  delete fXmax;
}

void AliFlowVZEROResults::Reset(){
  for(Int_t i=fNameVar->GetEntries();i>0;i--){
    delete fNameVar->At(i-1);
    fNameVar->RemoveAt(i-1);
  }

  for(Int_t i=fV2->GetEntries();i>0;i--){
    delete fV2->At(i-1);
    fV2->RemoveAt(i-1);
  }

  delete fNbinVar;
  delete fXmin;
  delete fXmax;

  fNbinVar = new TArrayI(0);
  fXmin = new TArrayF(0);
  fXmax = new TArrayF(0);

}

AliFlowVZEROResults::AliFlowVZEROResults(const AliFlowVZEROResults &old) :
  TNamed(old),
  fNbinVar(NULL),
  fXmin(NULL),
  fXmax(NULL),
  fNameVar(new TClonesArray("TNamed")),
  fV2(new TClonesArray("TProfile"))
{

  fNbinVar = new TArrayI(old.GetNvar());
  fXmin = new TArrayF(old.GetNvar());
  fXmax = new TArrayF(old.GetNvar());

  for(Int_t i=0; i<old.GetNhistos();i++){
    new((*fV2)[i]) TProfile(*((TProfile *) old.GetV2(i)));
  }

  for(Int_t i=0; i<old.GetNvar();i++){
    new((*fNameVar)[i]) TNamed(old.GetVarName(i),old.GetVarName(i));
  }

  for(Int_t i=0;i < GetNvar();i++){
    (*fNbinVar)[i] = (*old.fNbinVar)[i];
    (*fXmin)[i] = (*old.fXmin)[i];
    (*fXmax)[i] = (*old.fXmax)[i];
   }

}

AliFlowVZEROResults& AliFlowVZEROResults::operator=(const AliFlowVZEROResults &old){

  if(this != &old){
     printf("different\n");
   }

   for(Int_t i=0; i<old.GetNhistos();i++){
     new((*fV2)[i]) TProfile(*((TProfile *) old.GetV2(i)));
   }
   
   fNbinVar = new TArrayI(old.GetNvar());
   fXmin = new TArrayF(old.GetNvar());
   fXmax = new TArrayF(old.GetNvar());
   
   for(Int_t i=0;i < old.GetNvar();i++){
     (*fNbinVar)[i] = (*old.fNbinVar)[i];
     (*fXmin)[i] = (*old.fXmin)[i];
     (*fXmax)[i] = (*old.fXmax)[i];
   }
   
   fNameVar = new TClonesArray("TNamed");
   for(Int_t i=0; i<old.GetNvar();i++){
     new((*fNameVar)[i]) TNamed(old.GetVarName(i),old.GetVarName(i));
   }

  return *this;
}
  
Int_t AliFlowVZEROResults::GetNspecies() const{
  Int_t n = fV2->GetEntries();

  for(Int_t i=0;i < GetNvar();i++){
    n /= (*fNbinVar)[i];
  }

  return n;
}

void AliFlowVZEROResults::AddSpecies(const char *name,Int_t nXbin,const Double_t *bin){
  
  Bool_t kErr = kFALSE;
  for(Int_t i=0;i < GetNvar();i++){ // check the var ranges are set properly    
    if((*fNbinVar)[i] < 1 || (*fXmin)[i] >= (*fXmax)[i]){
    printf("var ranges are not set properly for variable %i please chek it before to define the species\n",i);
	kErr = kTRUE;
    }
  }
  if(kErr){
    printf("AddSpecies: NOTHING DONE\n");
    return;
  }

  Int_t ncomb = 1;
  for(Int_t i=0;i < GetNvar();i++){
    ncomb *= (*fNbinVar)[i];
  }

  char nameHisto[200];
  char title[300];
  for(Int_t i=0; i < ncomb;i++){
    snprintf(nameHisto,200,"%s_%s_%i",GetName(),name,i);
    snprintf(title,300,"%s",name);
    Int_t ncombTemp = i;
    for(Int_t j=0;j < GetNvar();j++){
      Int_t ibin = ncombTemp%(*fNbinVar)[j];
      snprintf(title,300,"%s_%04.1f<%s<%04.1f",title,(*fXmin)[j] + ((*fXmax)[j]-(*fXmin)[j])/(*fNbinVar)[j]*ibin,fNameVar->At(j)->GetName(),(*fXmin)[j] + ((*fXmax)[j]-(*fXmin)[j])/(*fNbinVar)[j]*(ibin+1));
      ncombTemp /= (*fNbinVar)[j];
    }

    new((*fV2)[GetNhistos()]) TProfile(nameHisto,title,nXbin,bin);
    ((TProfile *) GetV2(GetNhistos()-1))->GetXaxis()->SetTitle("p_{t} (GeV/c)");
    ((TProfile *) GetV2(GetNhistos()-1))->GetYaxis()->SetTitle("v_{2}");
  }
}

Int_t AliFlowVZEROResults::Add(const AliFlowVZEROResults *oth){
  if(GetNhistos() == oth->GetNhistos()){
    for(Int_t i=0;i < GetNhistos();i++){
      GetV2(i)->Add(oth->GetV2(i));
    }
    return 0;
  }
  else{
    printf("ADD error: number of objects is different (%i != %i)\n",GetNhistos(),oth->GetNhistos());
    return 1;
  }
}

void AliFlowVZEROResults::SetVarRange(Int_t ivar,Float_t xMin,Float_t xMax){
  if(!GetNhistos()){
    (*fXmin)[ivar]=xMin;
    (*fXmax)[ivar]=xMax;
  }
  else{ // to avoid different range among the histos
    printf("Ranges should be set before to define the species\nNOTHING DONE\n");
  }

}

void AliFlowVZEROResults::Fill(Int_t species,Float_t pt,Float_t v2,Float_t x[]){
  Int_t ncomb = 1;
  Int_t histo = 0;
 
  for(Int_t i=0;i < GetNvar();i++){
    Int_t ibin = GetBin(i,x[i]);
    if(ibin < 0 || ibin >= (*fNbinVar)[i]){
      printf("%i) %i not good w.r.t. %i (%f) (%f,%f)\n",i,ibin,(*fNbinVar)[i],x[i],(*fXmin)[i],(*fXmax)[i]);
      return;
    }
    histo += ncomb * ibin;
    ncomb *= (*fNbinVar)[i];
  }
  histo += species*ncomb;
  DirectFill(histo,pt,v2);
 };

TProfile *AliFlowVZEROResults::GetV2(Int_t species,Float_t x[]) const{
  Int_t ncomb = 1;
  Int_t histo = 0;
 
  for(Int_t i=0;i < GetNvar();i++){
    Int_t ibin = GetBin(i,x[i]);
    if(ibin < 0 || ibin >= (*fNbinVar)[i]){
      printf("%i) %i not good w.r.t. %i (%f) (%f,%f)\n",i,ibin,(*fNbinVar)[i],x[i],(*fXmin)[i],(*fXmax)[i]);
      return NULL;
    }
    histo += ncomb * ibin;
    ncomb *= (*fNbinVar)[i];
  }
  histo += species*ncomb;


  return GetV2(histo);

}

TProfile *AliFlowVZEROResults::GetV2(Int_t species,Float_t xMin[],Float_t xMax[]) const{
  if(GetNvar()){
    char title[300];
    char title2[300];
    Int_t ncomb = 1;
    for(Int_t i=0;i < GetNvar();i++){
      ncomb *= (*fNbinVar)[i];
    }

    TProfile *htemplate = GetV2(species*ncomb);
    TProfile *temp = new TProfile(*htemplate);
    temp->SetName("histo");
    temp->Reset();
    snprintf(title,300,"%i",species);
    for(Int_t i=0;i < GetNvar();i++){
      Int_t imin = GetBin(i,xMin[i]);
      if(imin < 0) imin = 0;
      else if(imin >= (*fNbinVar)[i]) imin = (*fNbinVar)[i]-1;
      Int_t imax = GetBin(i,xMax[i]);
      if(imax < imin) imax = imin;
      else if(imax >= (*fNbinVar)[i]) imax = (*fNbinVar)[i]-1;
      snprintf(title2,300,"%s",title);
      snprintf(title,300,"%s_%04.1f<%s<%04.1f",title2,
                (*fXmin)[i] + ((*fXmax)[i]-(*fXmin)[i])/(*fNbinVar)[i]*imin,
                fNameVar->At(i)->GetName(),
                (*fXmin)[i] + ((*fXmax)[i]-(*fXmin)[i])/(*fNbinVar)[i]*(imax+1));
    }
    temp->SetTitle(title);

    for(Int_t i=species*ncomb;i < (species+1)*ncomb;i++){
      Bool_t kGood = kTRUE;

      Int_t ncombTemp = i;
      for(Int_t j=0;j < GetNvar();j++){
	Int_t imin = GetBin(j,xMin[j]);
	if(imin < 0) imin = 0;
	else if(imin >= (*fNbinVar)[j]) imin = (*fNbinVar)[j]-1;
	Int_t imax = GetBin(j,xMax[j]);
	if(imax < imin) imax = imin;
	else if(imax >= (*fNbinVar)[j]) imax = (*fNbinVar)[j]-1;
	
	Int_t ibin = ncombTemp%(*fNbinVar)[j];
	ncombTemp /= (*fNbinVar)[j];

	if(ibin < imin || ibin > imax){
	  kGood = kFALSE;
	  j = GetNvar();
	}
      }

      if(kGood) temp->Add(GetV2(i));
    }
    return temp;
  }

  return GetV2(species);
}

Long64_t AliFlowVZEROResults::Merge(TCollection* list){
  Long64_t res=0;
  if (!list) return 0;
  if (list->IsEmpty()) return 0;

  TList *listObj = new TList();
  listObj->AddAll(list);

  for(Int_t i=0;i < listObj->GetEntries();i++){
    AliFlowVZEROResults *obj = (AliFlowVZEROResults *) listObj->At(i);
    Add(obj);
    res++;
  }
  return res;
}
