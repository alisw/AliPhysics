#include <TString.h>
#include <TMath.h>
#include "AliITSUSeed.h"
using namespace TMath;

ClassImp(AliITSUSeed)

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed() 
: fHitsPattern(0)
  ,fClID(0)
  ,fChi2Glo(0)
  ,fChi2Cl(0)
  ,fParent(0)
{
  // def c-tor
}

//_________________________________________________________________________
AliITSUSeed::~AliITSUSeed()
{
  // d-rot
}

//_________________________________________________________________________
AliITSUSeed::AliITSUSeed(const AliITSUSeed& src) 
  :AliExternalTrackParam(src)
  ,fHitsPattern(src.fHitsPattern)
  ,fClID(src.fClID)
  ,fChi2Glo(src.fChi2Glo)
  ,fChi2Cl(src.fChi2Cl)
  ,fParent(src.fParent) 
{
  // def c-tor
}

//_________________________________________________________________________
AliITSUSeed &AliITSUSeed::operator=(const AliITSUSeed& src) 
{
  // def c-tor
  if (this == &src) return *this;
  fClID        = src.fClID;
  fHitsPattern = src.fHitsPattern;
  fChi2Glo     = src.fChi2Glo;
  fChi2Cl      = src.fChi2Cl;
  fParent      = src.fParent;
  AliExternalTrackParam::operator=(src);
  return *this;
}

//_________________________________________________________________________
void AliITSUSeed::Print(Option_t* opt) const
{
  // print seed info
  int lr,cl = GetLrCluster(lr);
  printf("%cLr%d Cl:%4d Chi2Glo:%7.2f(%7.2f) Chi2Cl:",IsKilled() ? '-':' ',
	 lr,cl,GetChi2Glo(),GetChi2GloNrm());
  cl<0 ? printf("   NA  ") : printf("%7.2f",GetChi2Cl());
  printf(" |"); 
  for (int i=0;i<=12;i++) printf("%c",HasClusterOnLayer(i) ? '+':'-'); printf("|\n");
  TString opts = opt; opts.ToLower();
  if (opts.Contains("etp")) AliExternalTrackParam::Print();
  if (opts.Contains("parent") && GetParent()) GetParent()->Print(opt);
}

//______________________________________________________________________________
Float_t AliITSUSeed::GetChi2GloNrm() const
{
  int ndf = 2*GetNLayersHit() - 5;
  return ndf>0 ? fChi2Glo/ndf : fChi2Glo;
}


//______________________________________________________________________________
Int_t AliITSUSeed::Compare(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUSeed* sd = (const AliITSUSeed*)obj;
  const Float_t kTol = 1e-5;
  if (!IsKilled() && sd->IsKilled()) return -1;
  if ( IsKilled() &&!sd->IsKilled()) return  1;
  //
  if      (GetChi2Glo()+kTol<sd->GetChi2Glo()) return -1;
  else if (GetChi2Glo()-kTol>sd->GetChi2Glo()) return  1;
  return 0;
}

//______________________________________________________________________________
Bool_t AliITSUSeed::IsEqual(const TObject* obj)  const
{
  // compare clusters accodring to specific mode
  const AliITSUSeed* sd = (const AliITSUSeed*)obj;
  const Float_t kTol = 1e-5;
  if (IsKilled() != sd->IsKilled()) return kFALSE;
  return Abs(GetChi2Glo() - sd->GetChi2Glo())<kTol;
}
