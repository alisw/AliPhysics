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
  printf("Lr%d Cl:%4d Chi2Glo:%6.3f Chi2Cl:",lr,cl,GetChi2Glo());
  cl<0 ? printf("  NA  ") : printf("%6.3f",GetChi2Cl());
  printf(" |"); 
  for (int i=0;i<=lr;i++) printf("%c",HasClusterOnLayer(i) ? '+':'-');
  TString opts = opt; opts.ToLower();
  if (opts.Contains("etp")) AliExternalTrackParam::Print();
  if (opts.Contains("parent") && GetParent()) GetParent()->Print(opt);
}
