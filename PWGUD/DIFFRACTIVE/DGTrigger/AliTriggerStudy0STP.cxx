// -*- C++ -*-
// $Id$

#include "AliLog.h"
#include "AliTriggerStudy0STP.h"

ClassImp(AliTriggerStudy0STP);

Bool_t AliTriggerStudy0STP::CheckForTrackletPair(const TBits* mult, std::bitset<40> &phi) const
{
  phi.reset();
  
  const std::bitset<20> phi_L0 = ExtractPhi_L0(mult);
  const std::bitset<40> phi_L1 = ExtractPhi_L1(mult);
  
  const Int_t nOuter = mult->CountBits(400);
  const Int_t nInner = mult->CountBits()-nOuter;
   if (nOuter > 7 || nInner > 7)
     return false;
  
  for (std::vector<TablePhi>::const_iterator i=fLookupTablePhi.begin(), iend=fLookupTablePhi.end(); i!=iend; ++i) {
    if ((i->p.first  & phi_L0).any() &&
	(i->p.second & phi_L1).any()) 
      phi[i->phi] = 1;
  }
  AliDebug(3, Form("CheckForTrackletPair: L1  %s", phi_L1.to_string().c_str()));
  AliDebug(3, Form("CheckForTrackletPair: L0  %s", phi_L0.to_string().c_str()));
  AliDebug(3, Form("CheckForTrackletPair: phi %s", phi.to_string().c_str()));
  for (std::vector<TableDeltaPhi>::const_iterator i=fLookupTableDeltaPhi.begin(), iend=fLookupTableDeltaPhi.end(); i!=iend; ++i) {
    //     Printf(".....ForTrackletPair:     %s", i->p.to_string().c_str());
    if ((phi & i->p) == i->p)
      return true;
  }
  
  return false;
}

std::bitset<20> AliTriggerStudy0STP::CheckForVertex(const TBits* mult, const std::bitset<40> &phi, std::vector<int>& v) const
{
  std::bitset<20> vtx;
  for (std::vector<TableDeltaPhi>::const_iterator i=fLookupTableDeltaPhi.begin(), iend=fLookupTableDeltaPhi.end(); i!=iend; ++i) {
    if ((phi & i->p) == i->p) {
      const std::bitset<20> z1_L0 = ExtractZ_L0(mult, i->phi1);
      const std::bitset<20> z1_L1 = ExtractZ_L1(mult, i->phi1);
      
      const std::bitset<20> z2_L0 = ExtractZ_L0(mult, i->phi2);
      const std::bitset<20> z2_L1 = ExtractZ_L1(mult, i->phi2);
      
      std::bitset<20> v1, v2;
      for (std::vector<TableZ>::const_iterator j=fLookupTableZ.begin(), jend=fLookupTableZ.end(); j!=jend; ++j) {
	if ((j->p.first & z1_L0).any() && (j->p.second & z1_L1).any()) 
	  v1 |= j->vtx;
	if ((j->p.first & z2_L0).any() && (j->p.second & z2_L1).any()) 
	  v2 |= j->vtx;	
      }
      AliDebug(3, Form("VTX_L1: %s", z1_L1.to_string().c_str()));
      AliDebug(3, Form("VTX_L0: %s", z1_L0.to_string().c_str()));
      AliDebug(3, Form("VTX1:   %s", v1.to_string().c_str()));
      AliDebug(3, Form("VTX2:   %s", v2.to_string().c_str()));
      AliDebug(3, Form("VTX_L0: %s", z2_L0.to_string().c_str()));
      AliDebug(3, Form("VTX_L1: %s", z2_L1.to_string().c_str()));
      vtx |= (v1 & v2);
      
      const std::bitset<20> v12 = (v1 & v2);
      for (Int_t j=0; j<20; ++j)
	if (v12[j])
	  ++v[j];
    }
  }
  return vtx;
}


void AliTriggerStudy0STP::MakeTableZ()
{
  for (Int_t i=4; i<16; ++i) {
    for (Int_t j=0; j<20; ++j) {
      if (std::abs(i-j) > 10) continue;
      TableZ t;
      t.vtx[i] = 1;
      t.p.first[Int_t(i+3.9/7.6*(j-i)+0.5)] = 1;
      for (Int_t k=j-0; k<j+1; ++k) {
	if (k<0 || k>19) continue;
	t.p.second[k] = 1;
      }
      AliDebug(3, Form("i=%2d\n\t %s\n\t %s\n\t %s", i, t.p.second.to_string().c_str(), t.p.first.to_string().c_str(), t.vtx.to_string().c_str()));
      fLookupTableZ.push_back(t);
    }
  }
  AliDebug(3, Form("size: %ld", fLookupTableZ.size()));
  //   for (std::vector<TableZ>::const_iterator i=fLookupTableZ.begin(), end=fLookupTableZ.end(); i!=end; ++i) {
  //     Printf("vtx=%2d %s %s", i->vtx, i->p.first.to_string().c_str(), i->p.second.to_string().c_str());
  //   }
}

void AliTriggerStudy0STP::MakeTablePhi()
{
  for (Int_t k=0; k<40; ++k) {
    const Int_t i=k/2;
    TablePhi t;
    t.phi          = k;
    t.p.first[ i   %20]  = 1;
    t.p.first[(i+1)%20]  = 1;
    t.p.second[ k   %40] = 1;
    t.p.second[(k+1)%40] = 1;
    t.p.second[(k+2)%40] = 1;
    AliDebug(3, Form("i,k= %d %d", i,k));
    TString s0 = "TablePhi: L0  ";
    for (Int_t l=0; l<20; ++l)
      s0 += TString::Format("%2d", Int_t(t.p.first[19-l]));
    AliDebug(3, s0.Data());
    AliDebug(3, Form("TablePhi: L1  %s", t.p.second.to_string().c_str()));
    fLookupTablePhi.push_back(t);
  }
}

void AliTriggerStudy0STP::MakeTableDeltaPhi(Int_t deltaPhiMin, Int_t deltaPhiMax)
{
  for (Int_t k=0; k<40; ++k) {
    for (Int_t j=deltaPhiMin; j<deltaPhiMax+1; ++j) {
      TableDeltaPhi t;
      t.phi1 = k;
      t.phi2 = ((k+j)%40);
      t.p[t.phi1] = 1;
      t.p[t.phi2] = 1;
      AliDebug(3, Form("TableDeltaPhi: %s", t.p.to_string().c_str()));
      fLookupTableDeltaPhi.push_back(t);
    }
  }
}

std::bitset<20> AliTriggerStudy0STP::ExtractPhi_L0(const TBits *mult)
{
  std::bitset<20> b;
  for (Int_t i=0; i<400; ++i)
    if (mult->TestBitNumber(i))
      b[i/20] = 1;
  return b;
}
std::bitset<40> AliTriggerStudy0STP::ExtractPhi_L1(const TBits *mult)
{
  std::bitset<40> b;
  for (Int_t i=400; i<1200; ++i)
    if (mult->TestBitNumber(i))
      b[(i-400)/20] = 1;
  return b;
}
std::bitset<20> AliTriggerStudy0STP::ExtractZ_L0(const TBits *mult, Int_t k)
{
  const Int_t i = k/2;
  std::bitset<20> b;
  for (Int_t l=0; l<20; ++l)
    if (mult->TestBitNumber(20*  i       +l) ||
	mult->TestBitNumber(20*((i+1)%20)+l))
      b[l] = 1;  
  return b;
}
std::bitset<20> AliTriggerStudy0STP::ExtractZ_L1(const TBits *mult, Int_t k)
{
  std::bitset<20> b;
  for (Int_t l=0; l<20; ++l)
    if (mult->TestBitNumber(400+20*  k       +l) ||
	mult->TestBitNumber(400+20*((k+1)%40)+l) ||
	mult->TestBitNumber(400+20*((k+2)%40)+l))
      b[l] = 1;
  return b;
}
