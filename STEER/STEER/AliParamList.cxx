#include "AliParamList.h"
#include <TString.h>
#include <AliLog.h>

ClassImp(AliParamList)

//_____________________________________________________________________
AliParamList::AliParamList(Int_t n, const Double_t *parVal)
: fID(0)
  , fNPar(0)
  , fNames(0)
  , fParams(0)
{
  // def. c-tor
  //
  if (n>0) {
    SetNParams(n);
    if (parVal) SetParameters(parVal);
  }
}

//_____________________________________________________________________
AliParamList::AliParamList(const AliParamList& src)
  : TNamed(src)
  , fID(src.fID)
  , fNPar(src.fNPar)
  , fNames(0)
  , fParams(0)
{
  // copy c-tor
  if (fNPar>0) {
    fParams = new Double_t[fNPar];
    for (int i=fNPar;i--;) {
      fNames[i] = src.fNames[i];
      fParams[i] = src.fParams[i];
    }
  }
}

//_____________________________________________________________________
AliParamList& AliParamList::operator=(const AliParamList& src)
{
  // copy op.
  if (this != &src) {
    this->~AliParamList();
    new(this) AliParamList(src);
  }
  return *this;
  //
}

//_____________________________________________________________________
AliParamList::~AliParamList()
{
  // d-tor
  delete[] fNames;
  delete[] fParams;
}

//_____________________________________________________________________
void AliParamList::SetNParams(Int_t n)
{
  // init params structure
  if (fNPar) AliFatal(Form("N params was already set to %d",fNPar));
  fNPar = n;
  fParams = new Double_t[fNPar];
  fNames = new TString[fNPar];
  for (int i=fNPar;i--;) fParams[i] = 0.;
  //
}

//_____________________________________________________________________
void AliParamList::SetParName(Int_t i,const char* nm)
{
  // assign param name
  if (i<0||i>=fNPar) AliFatal(Form("Param %d accessed while the range is %d : %d",i,0,fNPar));
  fNames[i] = nm;
}

//_____________________________________________________________________
void AliParamList::SetParameter(Int_t i, Double_t v, const char* nm)
{
  // assign param value and optionally name
  if (i<0||i>=fNPar) AliFatal(Form("Param %d accessed while the range is %d : %d",i,0,fNPar));
  fParams[i] = v;
  fNames[i] = nm;
}

//_____________________________________________________________________
void AliParamList::Print(Option_t *) const
{
  // print itself
  printf("ParamList#%d/%d %s %s\n",fID,GetUniqueID(),GetName(),GetTitle());
  for (int i=0;i<fNPar;i++) printf("#%2d\t%20s\t%e\n",i,GetParName(i),GetParameter(i));
}

//_____________________________________________________________________
const Char_t* AliParamList::GetParName(Int_t i) const 
{
  // get par name
  return (fNames[i].IsNull()) ? Form("par%d",i) : fNames[i].Data();;
}
