#include "AliHBTMonSplitting.h"

ClassImp(AliHBTMonSplittingQosl)

AliHBTMonSplittingQosl::AliHBTMonSplittingQosl(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                               Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  Rename("splitosl","Q_{out}-Q_{side}-Q_{long} Splitting Monitoring Function");
}

void   AliHBTMonSplittingQosl::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  AliVAODParticle* p1 = partpair->Particle1();
  AliVAODParticle* p2 = partpair->Particle2();
  
  if (p1->Px() != p2->Px()) return;
  if (p1->Py() != p2->Py()) return;
  if (p1->Pz() != p2->Pz()) return;
  
  Double_t out = trackpair->GetQOutLCMS();
  Double_t side = trackpair->GetQSideLCMS();
  Double_t lon = trackpair->GetQLongLCMS();
    
  fNumerator->Fill(out,side,lon);//here we fill in q's corresponding to track pair 
                                          //weight calculated for the simulated one
}

ClassImp(AliHBTMonSplittingDptDthetaDphi)

AliHBTMonSplittingDptDthetaDphi::AliHBTMonSplittingDptDthetaDphi(Int_t nXbins, Double_t maxXval, Double_t minXval,
                                               Int_t nYbins, Double_t maxYval, Double_t minYval,
                                                Int_t nZbins, Double_t maxZval, Double_t minZval):
 AliHBTTwoPairFctn3D(nXbins,maxXval,minXval,nYbins,maxYval,minYval,nZbins,maxZval,minZval)
{
  //ctor
  Rename("splitdpdthedphi","\\Deltap_{t}-\\Delta\\theta-\\Delta\\phi Splitting Monitoring Function");
}

void   AliHBTMonSplittingDptDthetaDphi::ProcessSameEventParticles(AliHBTPair* trackpair, AliHBTPair* partpair)
{
  AliVAODParticle* p1 = partpair->Particle1();
  AliVAODParticle* p2 = partpair->Particle2();
  
  if (p1->Px() != p2->Px()) return;
  if (p1->Py() != p2->Py()) return;
  if (p1->Pz() != p2->Pz()) return;
  
  Double_t dpt = trackpair->GetDeltaPt();
  Double_t dphi = trackpair->GetDeltaPhi();
  Double_t dtheta = trackpair->GetDeltaTheta();
    
  fNumerator->Fill(dpt, dphi, dtheta);//here we fill in q's corresponding to track pair 
                                          //weight calculated for the simulated one
}
