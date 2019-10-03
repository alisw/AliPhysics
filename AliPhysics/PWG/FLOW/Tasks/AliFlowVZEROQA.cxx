#include<stdio.h>
#include "AliFlowVZEROQA.h"
#include "TList.h"

ClassImp(AliFlowVZEROQA);

AliFlowVZEROQA::AliFlowVZEROQA(const char *name,const Int_t nvar,const Int_t* binVar) :
  TNamed(name,name),
  fNbinVar(new TArrayI(nvar)),
  fXmin(new TArrayF(nvar)),
  fXmax(new TArrayF(nvar)),
  fNameVar(new TClonesArray("TNamed")),
  fQA(new TClonesArray("TH2F"))
{

  for(Int_t i=0;i < GetNvar();i++){
    (*fNbinVar)[i] = binVar[i];
  }
  
  for(Int_t i=0; i<GetNvar();i++){
    new((*fNameVar)[i]) TNamed("","");
  }
}

AliFlowVZEROQA::AliFlowVZEROQA() :
  TNamed("qa","qa"),
  fNbinVar(new TArrayI(0)),
  fXmin(new TArrayF(0)),
  fXmax(new TArrayF(0)),
  fNameVar(new TClonesArray("TNamed")),
  fQA(new TClonesArray("TH2F"))
{
  
}

AliFlowVZEROQA::~AliFlowVZEROQA(){
  for(Int_t i=fNameVar->GetEntries();i>0;i--){
    delete fNameVar->At(i-1);
    fNameVar->RemoveAt(i-1);
  }

  for(Int_t i=fQA->GetEntries();i>0;i--){
    delete fQA->At(i-1);
    fQA->RemoveAt(i-1);
  }

  delete fNbinVar;
  delete fXmin;
  delete fXmax;
}

void AliFlowVZEROQA::Reset(){
  for(Int_t i=fNameVar->GetEntries();i>0;i--){
    delete fNameVar->At(i-1);
    fNameVar->RemoveAt(i-1);
  }

  for(Int_t i=fQA->GetEntries();i>0;i--){
    delete fQA->At(i-1);
    fQA->RemoveAt(i-1);
  }

  delete fNbinVar;
  delete fXmin;
  delete fXmax;

  fNbinVar = new TArrayI(0);
  fXmin = new TArrayF(0);
  fXmax = new TArrayF(0);

}

AliFlowVZEROQA::AliFlowVZEROQA(const AliFlowVZEROQA &old) :
  TNamed(old),
  fNbinVar(NULL),
  fXmin(NULL),
  fXmax(NULL),
  fNameVar(new TClonesArray("TNamed")),
  fQA(new TClonesArray("TH2F"))
{

  fNbinVar = new TArrayI(old.GetNvar());
  fXmin = new TArrayF(old.GetNvar());
  fXmax = new TArrayF(old.GetNvar());

  for(Int_t i=0; i<old.GetNhistos();i++){
    new((*fQA)[i]) TH2F(*((TH2F *) old.GetQA(i)));
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

AliFlowVZEROQA& AliFlowVZEROQA::operator=(const AliFlowVZEROQA &old){

  if(this != &old){
     printf("different\n");
   }

   for(Int_t i=0; i<old.GetNhistos();i++){
     new((*fQA)[i]) TH2F(*((TH2F *) old.GetQA(i)));
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
  
Int_t AliFlowVZEROQA::GetNspecies() const{
  Int_t n = fQA->GetEntries();

  for(Int_t i=0;i < GetNvar();i++){
    n /= (*fNbinVar)[i];
  }

  return n;
}

void AliFlowVZEROQA::AddSpecies(const char *name,Int_t nXbin,const Double_t *xbin,Int_t nYbin,const Double_t *ybin){
  
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
  char title2[300];
  for(Int_t i=0; i < ncomb;i++){
    snprintf(nameHisto,200,"%s_%s_%i",GetName(),name,i);
    snprintf(title,300,"%s",name);
    Int_t ncombTemp = i;
    for(Int_t j=0;j < GetNvar();j++){
      Int_t ibin = ncombTemp%(*fNbinVar)[j];
      snprintf(title2,300,"%s",title);
      snprintf(title,300,"%s_%04.1f<%s<%04.1f",title2,(*fXmin)[j] + ((*fXmax)[j]-(*fXmin)[j])/(*fNbinVar)[j]*ibin,fNameVar->At(j)->GetName(),(*fXmin)[j] + ((*fXmax)[j]-(*fXmin)[j])/(*fNbinVar)[j]*(ibin+1));
      ncombTemp /= (*fNbinVar)[j];
    }

    new((*fQA)[GetNhistos()]) TH2F(nameHisto,title,nXbin,xbin,nYbin,ybin);
    ((TH2F *) GetQA(GetNhistos()-1))->GetXaxis()->SetTitle("N_{#sigma}^{TPC}");
    ((TH2F *) GetQA(GetNhistos()-1))->GetYaxis()->SetTitle("N_{#sigma}^{TOF}");
  }
}

Int_t AliFlowVZEROQA::Add(const AliFlowVZEROQA *oth){
  if(GetNhistos() == oth->GetNhistos()){
    for(Int_t i=0;i < GetNhistos();i++){
      GetQA(i)->Add(oth->GetQA(i));
    }
    return 0;
  }
  else{
    printf("ADD error: number of objects is different (%i != %i)\n",GetNhistos(),oth->GetNhistos());
    return 1;
  }
}

void AliFlowVZEROQA::SetVarRange(Int_t ivar,Float_t xMin,Float_t xMax){
  if(!GetNhistos()){
    (*fXmin)[ivar]=xMin;
    (*fXmax)[ivar]=xMax;
  }
  else{ // to avoid different range among the histos
    printf("Ranges should be set before to define the species\nNOTHING DONE\n");
  }

}

void AliFlowVZEROQA::Fill(Int_t species,Float_t var1,Float_t var2,Float_t x[]){
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
  DirectFill(histo,var1,var2);
 };

TH2F *AliFlowVZEROQA::GetQA(Int_t species,Float_t x[]) const{
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


  return GetQA(histo);

}

TH2F *AliFlowVZEROQA::GetQA(Int_t species,Float_t xMin[],Float_t xMax[]) const{
  if(GetNvar()){
    char title[300];
    char title2[300];
    Int_t ncomb = 1;
    for(Int_t i=0;i < GetNvar();i++){
      ncomb *= (*fNbinVar)[i];
    }

    TH2F *htemplate = GetQA(species*ncomb);
    TH2F *temp = new TH2F(*htemplate);
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

      if(kGood) temp->Add(GetQA(i));
    }
    return temp;
  }

  return GetQA(species);
}

Long64_t AliFlowVZEROQA::Merge(TCollection* list){
  Long64_t res=0;
  if (!list) return 0;
  if (list->IsEmpty()) return 0;

  TList *listObj = new TList();
  listObj->AddAll(list);

  for(Int_t i=0;i < listObj->GetEntries();i++){
    AliFlowVZEROQA *obj = (AliFlowVZEROQA *) listObj->At(i);
    Add(obj);
    res++;
  }
  return res;
}
