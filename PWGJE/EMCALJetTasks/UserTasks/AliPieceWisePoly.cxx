#include "TF1.h"
#include "TString.h"
#include "TFile.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TGrid.h"

#include <cmath>

#include "AliPID.h"

#include "AliPieceWisePoly.h"

#include <string.h>
using namespace std;

AliPieceWisePoly::AliPieceWisePoly(Int_t parts, Double_t* cutxvalues, Int_t* polys, Double_t xmin, Double_t xmax,  Double_t* params, Int_t smooth)
  : doSmoothing(2)
  , nParts(parts)
  , cuts(0x0)
  , nOfFreeParams(0)
  , polyParameters(0x0)
  , piecewisepolynom(0x0)
  {
  nParts = parts;
  if (parts > 1) {
    cuts = new Double_t[nParts-1];
    
    for (Int_t i=0;i<nParts-1;i++) 
      cuts[i] = cutxvalues[i];
  }
  
  polyParameters = new Int_t[nParts];
  for (Int_t i=0;i<nParts;++i)
    polyParameters[i]=polys[i];
  
  doSmoothing = smooth;
  
  TString functionString = "";
  
  Int_t remParts = nParts;
  Int_t previousNOfParameters = 0;
  while (remParts > 1) {
    functionString += TString::Format("x < %f ? pol%d(%d) : ",cuts[nParts-remParts],polyParameters[nParts-remParts]-1,previousNOfParameters);
    previousNOfParameters += polyParameters[nParts-remParts];
    remParts--;
  }
  functionString += TString::Format("pol%d(%d)",polyParameters[nParts-remParts]-1,previousNOfParameters);
  piecewisepolynom = new TF1("piecewisepolynom",functionString.Data(),xmin,xmax);
  nOfFreeParams = previousNOfParameters += polyParameters[nParts-remParts];
  
  if (doSmoothing==1)
    nOfFreeParams -= nParts-1;
  else if (doSmoothing==2)
    nOfFreeParams -= 2 * (nParts-1);
  
  SetParam(params);
}

AliPieceWisePoly::~AliPieceWisePoly() {
  delete [] cuts;
  cuts = 0x0;
  delete[] piecewisepolynom;
  piecewisepolynom = 0x0;
}

void AliPieceWisePoly::SetParam(Double_t* params) {
  if (!params)
    return;
  
  Bool_t changed = kFALSE;
  
  for (Int_t j=0;j<polyParameters[0];++j) {
    if (changed) {
      piecewisepolynom->SetParameter(j,params[j]);
    }
    else if (piecewisepolynom->GetParameter(j) != params[j]) {
        changed = kTRUE;
        piecewisepolynom->SetParameter(j,params[j]);
    }
  }

  Int_t parNumber = polyParameters[0];
  Int_t freeparNumber = polyParameters[0];
  
  for (Int_t i=1;i<nParts;++i) {
    for (Int_t j=doSmoothing;j<polyParameters[i];++j) {
      Int_t internalparNumber = parNumber + j;
      if (changed) {
        piecewisepolynom->SetParameter(internalparNumber,params[freeparNumber]);
      }
      else if (piecewisepolynom->GetParameter(internalparNumber) != params[freeparNumber]) {
          changed = kTRUE;
          piecewisepolynom->SetParameter(internalparNumber,params[freeparNumber]);
      }
      freeparNumber++;
    }
    if (changed) {
      for (Int_t j=doSmoothing-1;j>=0;--j) 
        piecewisepolynom->SetParameter(parNumber+j,SumUp(cuts[i-1],parNumber-polyParameters[i-1],parNumber-1,j,j)-SumUp(cuts[i-1],parNumber,parNumber+polyParameters[i]-1,j,j+1));
    }
    
    parNumber += polyParameters[i];
  }
  return;
}

Double_t AliPieceWisePoly::SumUp(Double_t constant, Int_t start, Int_t end, Int_t derivative, Int_t startn) {
  Double_t sum = 0.0;
  if (derivative == 1) {
    for (Int_t i=startn;i<=end-start;++i) {
      Double_t potence = 1.0;
      for (Int_t j=1;j<i;++j)
        potence *= constant;
      
      sum += i * piecewisepolynom->GetParameter(start+i) * potence;
    }
  }
  else {
    for (Int_t i=startn;i<=end-start;++i) {
      Double_t potence = 1.0;
      for (Int_t j=1;j<=i;++j)
        potence *= constant;
      
      sum += piecewisepolynom->GetParameter(start+i) * potence;
    }
  }
  return sum;
}

TF1* AliPieceWisePoly::GetPartFunction(Int_t i) {
  if (i>nParts)
    return 0x0;
  
  TString functionString = TString::Format("pol%d(0)",polyParameters[i]-1);
  Int_t previousNOfParameters = 0;
  for (Int_t j=0;j<i;++j) 
    previousNOfParameters += polyParameters[j];
  
  TF1* func = new TF1(TString::Format("func_%d",i).Data(),functionString.Data(),0,1);
  
  for (Int_t j=0;j<polyParameters[i];++j)
    func->SetParameter(j,piecewisepolynom->GetParameter(previousNOfParameters+j));
  
  return func;
}


double AliPieceWisePoly::Eval (double x, double* p) {
  SetParam(p);
  return piecewisepolynom->Eval(x);
}

TString AliPieceWisePoly::ReadFSParameters(TString parameterFile, TF1** effFunctions) 
{
  if(parameterFile.Contains("alien://") && !gGrid)
    TGrid::Connect("alien://");
  
  TString joinedString = "";
  char separator = ':';
  
  TFile* f = TFile::Open(parameterFile.Data(), "READ");
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    for (Int_t charge=0;charge<=1;++charge) {     
      TString name = TString::Format("fastSimulationParameters_%s_%s", AliPID::ParticleShortName(species), charge ? "pos" : "neg");
			
      TNamed* cont = (TNamed*)f->FindObjectAny(name.Data());

      if (!name)
        cout << "Could not find " << name << "." << endl;
      
      TString parString = TString(cont->GetTitle());
      joinedString += parString + TString(separator);
      
      TString nameFunction = TString::Format("fastSimulationFunction_%s_%s", AliPID::ParticleShortName(species), charge ? "pos" : "neg");
      TF1* func = GetFSFunctionFromString(parString, nameFunction);
      
      effFunctions[2*species + charge] = func;
    }
  }
  f->Close();
  delete f;
  joinedString.Remove(TString::kTrailing, separator);
  return joinedString;
}

void AliPieceWisePoly::ReadFSParametersFromString(TString parameterString, TF1** effFunctions)
{
  TObjArray* parameters = parameterString.Tokenize(":");
  
  Int_t position = 0;
  for (Int_t species=0;species<AliPID::kSPECIES;++species) {
    for (Int_t charge=0;charge<=1;++charge) {     
      TString parString = ((TObjString*)(parameters->At(position)))->GetString();
      
      TString nameFunction = TString::Format("fastSimulationFunction_%s_%s", AliPID::ParticleShortName(species), charge ? "pos" : "neg");
      TF1* func = GetFSFunctionFromString(parString, nameFunction);
      
      effFunctions[2*species + charge] = func;
      position++;
    }
  }
}


TF1* AliPieceWisePoly::GetFSFunctionFromString(TString parameters, TString nameFunction) 
{
  TObjArray* arrPar = parameters.Tokenize(";");
  Int_t nOfParts = (((TObjString*)(arrPar->At(0)))->GetString()).Atoi();
  Double_t* cuts = new Double_t[nOfParts-1];
  Int_t* nparameters = new Int_t[nOfParts]; 
  
  for (Int_t part = 0;part<nOfParts - 1;++part) {
    nparameters[part] = (((TObjString*)(arrPar->At(nOfParts + part)))->GetString()).Atoi();
    cuts[part] = (((TObjString*)(arrPar->At(1 + part)))->GetString()).Atof();
  }
  nparameters[nOfParts - 1] = (((TObjString*)(arrPar->At(2*nOfParts -1)))->GetString()).Atoi();
  
  AliPieceWisePoly* pwp = new AliPieceWisePoly(nOfParts,cuts,nparameters,0,50,0x0,2);
  TF1* function = new TF1(nameFunction.Data(),pwp,0,50,pwp->GetNOfParam());
  for (Int_t param=0;param<pwp->GetNOfParam();++param) {
    function->SetParameter(param,(((TObjString*)(arrPar->At(2*nOfParts + param)))->GetString()).Atof());
  }
  delete [] nparameters;
  delete [] cuts;
  delete arrPar;
  return function;
}

  
