#ifndef ALIJETPARTICLESREADERKINE_H
#define ALIJETPARTICLESREADERKINE_H

/* $Id$ */

//_______________________________________________________________________
/////////////////////////////////////////////////////////////////////////
//
// class AliJetParticlesReaderKine
//
// Reader for Kinematics
//
// loizides@ikf.uni-frankfurt.de
//
/////////////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include <TString.h>
#include <TParticle.h>
#include "AliJetParticlesReader.h"
class AliRunLoader;

class AliJetParticlesReaderKine: public AliJetParticlesReader
{
  public:
  AliJetParticlesReaderKine();
  AliJetParticlesReaderKine(TString&);
  AliJetParticlesReaderKine(TObjArray*, const Char_t *filename="galice.root");
  virtual ~AliJetParticlesReaderKine();

  void Rewind();
  void SetNeutral(Bool_t b){fNeutral=b;}
  void SetCharged(Bool_t b){fCharged=b;}
  void SetEM(Bool_t b){fEM=b;}
  void SetUseTracks(Bool_t b){fUseTracks=b;}

  protected:
  Int_t  ReadNext();
  Int_t  OpenFile(Int_t n);
  Bool_t IsAcceptedParticle(TParticle *p) const;
   
  TString       fFileName;  //file name of galice 
  AliRunLoader* fRunLoader; //!pointer to loader

  Bool_t  fNeutral; //neutral cut
  Bool_t  fCharged; //charged cut
  Bool_t  fEM;      //em (e+,e-,gamma) cut
  Bool_t  fUseTracks; // use ntracks instead of nprimaries
  ClassDef(AliJetParticlesReaderKine,1)
};

inline Bool_t AliJetParticlesReaderKine::IsAcceptedParticle(TParticle *p) const
{
  Int_t pcode=p->GetPdgCode();  
  if ((pcode==11)||(pcode==-11)||(pcode==22)) {
    if(!fEM) return kFALSE;
  }  else {
    TParticlePDG *pdg=p->GetPDG();
    Int_t ch=(Int_t)pdg->Charge(); 
    if((!fCharged)&&(ch)) return kFALSE;
    if((!fNeutral)&&(!ch)) return kFALSE;
  }

  Float_t eta=0.;//p->Eta();
  Float_t pz=p->Pz();
  //if(TMath::Abs(pz)>150.) return kFALSE;
  Float_t pabs=p->P();
  if(pabs-TMath::Abs(pz)>1.e-8) eta=0.5*TMath::Log((pabs+pz)/(pabs-pz));
  else return kFALSE;

  if((eta<fEtaMin)||(eta>fEtaMax)) return kFALSE;

  Float_t phi=p->Phi();
  if((phi<fPhiMin)||(phi>fPhiMax)) return kFALSE;

  Float_t pt=p->Pt();
  if((pt<fPtMin)||(pt>fPtMax)) return kFALSE;

  return kTRUE;
}
#endif
