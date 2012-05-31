#include "AliHLTGenerator.h"


ClassImp(AliHLTGenerator)

AliHLTGenerator::AliHLTGenerator()
{

}

AliHLTGenerator::~AliHLTGenerator()
{

}

void AliHLTGenerator::Generate()
{
  Float_t pt[3];
  for(Int_t i=0; i<fNpart; i++)
    {
      TParticle *p = part[i];
      gAlice->SetTrack(fTrackIt,-1,p->GetPdgCode(),fP,fOrigin.GetArray(),polar,0,kPPrimary,nt);
      
    }

}
