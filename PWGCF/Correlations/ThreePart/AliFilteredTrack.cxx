/*************************************************************************
* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
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
#include "AliFilteredTrack.h"
#include "TString.h"
#include "TObjArray.h"
#include "TMCParticle.h"
#include "AliMCParticle.h"
#include <iostream>
#include <memory>

using namespace std;

ClassImp(AliFilteredTrack)

AliFilteredTrack::AliFilteredTrack()
  : AliVTrack()
  , fPt(0.0)
  , fPhi(0.0)
  , fTheta(0.0)
  , feff(-1.0)
  , fP()
  , fOneOverPt(0.0)
  , fPtot(0.0)
  , fEta(-20.0)
{
  // default constructor
  // 
  //

  memset(fP, 0, sizeof(fP));
  Calculate();
  SetBit(kSignBit); // default 'plus'
  SetBit(kLeading,false);//is not leading
}

// AliFilteredTrack::AliFilteredTrack(float pt, float phi, float theta, int charge)
//   : AliVTrack()
//   , fPt(pt)
//   , fPhi(phi)
//   , fTheta(theta)
//   , fP()
//   , fOneOverPt(0.0)
//   , fPtot(0.0)
//   , fEta(-20.0)
// {
//   // constructor
//   Calculate;
//   if (charge==-1) SetBit(kSignBit);
//   else ResetBit(kSignBit);
// }

AliFilteredTrack::AliFilteredTrack(const AliFilteredTrack& other)
  : AliVTrack(other)
  , fPt(0.0)
  , fPhi(0.0)
  , fTheta(0.0)
  , feff(other.feff)
  , fP()
  , fOneOverPt(0.0)
  , fPtot(0.0)
  , fEta(-20.0)
{
  // copy constructor
  Set(other);
  Calculate();
  SetBit(kLeading,false);//is not leading
  
}

AliFilteredTrack::AliFilteredTrack(const AliVParticle& other)
  : AliVTrack()
  , fPt(0.0)
  , fPhi(0.0)
  , fTheta(0.0)
  , feff(-1.0)
  , fP()
  , fOneOverPt(0.0)
  , fPtot(0.0)
  , fEta(-20.0)
{
  // copy constructor
  Set(other);
  SetBit(kLeading,false);//is not leading
}

AliFilteredTrack::AliFilteredTrack(const TMCParticle& other)
  : AliVTrack()
  , fPt(0.0)
  , fPhi(0.0)
  , fTheta(0.0)
  , feff(-1.0)
  , fP()
  , fOneOverPt(0.0)
  , fPtot(0.0)
  , fEta(-20.0)
{
  // copy constructor
  float p[3]={other.GetPx(), other.GetPy(), other.GetPz()};
  Set(p);
  if (other.GetKF()<0) SetBit(kSignBit);
  else ResetBit(kSignBit);
  SetBit(kLeading,false);//is not leading
}

AliFilteredTrack::AliFilteredTrack(const AliAODTrack& other)
  : AliVTrack()
  , fPt(0.0)
  , fPhi(0.0)
  , fTheta(0.0)
  , feff(-1.0)
  , fP()
  , fOneOverPt(0.0)
  , fPtot(0.0)
  , fEta(-20.0)
{
  // copy constructor
  Set(other);
  SetBit(kLeading,false);//is not leading
}

AliFilteredTrack& AliFilteredTrack::operator=(const AliFilteredTrack& other)
{
  // assignment operator
  // forward to copy constructor, but should be implemented if the class
  // members get more complex
  if (this==&other) return *this;
  this->~AliFilteredTrack();
  new (this) AliFilteredTrack(other);
  return *this;
}

AliFilteredTrack& AliFilteredTrack::operator=(const AliVParticle& other)
{
  // assignment operator
  // forward to copy constructor, but should be implemented if the class
  // members get more complex
  if (this==&other) return *this;
  this->~AliFilteredTrack();
  new (this) AliFilteredTrack(other);
  return *this;
}

AliFilteredTrack::~AliFilteredTrack()
{
  // destructor
  //
  //

}

void AliFilteredTrack::Calculate(bool bCalculateMommentumComponents)
{
  // calculate internal members once from the momentum vector
  if (bCalculateMommentumComponents) {
    fP[kX]=fPt * TMath::Cos(fPhi);
    fP[kY]=fPt * TMath::Sin(fPhi);
    if (TMath::Tan(fTheta)>0.) fP[kZ]=fPt / TMath::Tan(fTheta);
    else if (fPt>0.) fP[kZ]=-999;
    else fP[kZ]=0.;
  }
  if (fPt>0) fOneOverPt=1/fPt;
  else fOneOverPt=-999.;
  fPtot = TMath::Sqrt(fPt*fPt + fP[kZ]*fP[kZ]);
  const double dtheta=9.07998594625855556e-05; // resp rapidity +/-10
  if (fTheta>dtheta && fTheta<(TMath::Pi()-dtheta)) {
    fEta = -TMath::Log(TMath::Tan(fTheta/2.));
  } else {
    fEta = -TMath::Log(TMath::Tan(dtheta/2.));
  }
}

void AliFilteredTrack::Set(const AliVParticle& track)
{
  // calculate internal members once from the momentum vector
  fPt = track.Pt();
  fPhi = track.Phi();
  fTheta = track.Theta();
  SetCharge(track);
  
  Calculate();
}

void AliFilteredTrack::Set(const AliAODTrack& track)
{
  // calculate internal members once from the momentum vector
  fPt = track.Pt();
  fPhi = track.Phi();
  fTheta = track.Theta();
  SetCharge(track);
  SetAODFilterBits(&track);
  
  Calculate();
}

void AliFilteredTrack::SetCharge(const AliVParticle& track)
{
  // charge is stored in one bit of TObject
  if (track.Charge()==-1) SetBit(kSignBit);
  else if (track.Charge()==-3&&dynamic_cast<const AliMCParticle*>(&track)) SetBit(kSignBit);
  else ResetBit(kSignBit);
}

void AliFilteredTrack::SetAODFilterBits(const AliAODTrack* t)
{
  // The filter bits are stored in bits of TObject
  if(t->IsHybridGlobalConstrainedGlobal()){SetBit(kGlobalHybrid);}
  else {ResetBit(kGlobalHybrid);}
  if(t->TestFilterBit(BIT(4)))SetBit(kBIT4);
  else ResetBit(kBIT4);
  if(t->TestFilterBit(BIT(5)))SetBit(kBIT5);
  else ResetBit(kBIT5);
  if(t->TestFilterBit(BIT(6)))SetBit(kBIT6);
  else ResetBit(kBIT6);  
  ResetBit(kMC);
}


void AliFilteredTrack::Clear(Option_t * /*option*/)
{
  /// overloaded from TObject: cleanup
}

void AliFilteredTrack::Print(Option_t *option) const
{
  /// overloaded from TObject: print info
  TString strOption(option);
  const char* fullIndent="          ";
  int maxIndent=strlen(fullIndent);
  const char* indent=fullIndent+maxIndent;
  std::unique_ptr<TObjArray> tokens(strOption.Tokenize(" "));
  if (!tokens.get()) return;
  for (int i=0; i<tokens->GetEntriesFast(); i++) {
    if (!tokens->At(i)) continue;
    const char* key="";
    TString arg=tokens->At(i)->GetName();

    key="indent=";
    if (arg.BeginsWith(key)) {
      arg.ReplaceAll(key, "");
      int argf=arg.Atoi();
      indent=fullIndent;
      if (argf<maxIndent)
	indent+=maxIndent-argf;
    }
  }

  cout << indent << "Pt:    " << Pt() << endl;
  cout << indent << "Phi:   " << Phi() << endl;
  cout << indent << "Theta: " << Theta() << endl;
  cout << indent << "Charge:" << Charge  () << endl;
  cout << indent << "px:    " << Px() << endl;
  cout << indent << "py:    " << Py() << endl;
  cout << indent << "pz:    " << Pz() << endl;
}

// void AliFilteredTrack::Streamer(TBuffer &buffer)
// {
//   // custom streamer to calculate members after reading the object 
//   if (buffer.IsReading()) {
//     buffer.ReadClassBuffer(AliFilteredTrack::Class(),this);
//     Calculate();
//   } else {
//     buffer.WriteClassBuffer(AliFilteredTrack::Class(),this);
//   }
// }
