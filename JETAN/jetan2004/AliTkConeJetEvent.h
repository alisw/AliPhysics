// $Id$

#ifndef ALITKCONEJETEVENT_H
#define ALITKCONEJETEVENT_H

#include <Riostream.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TString.h>
#include "AliTkConeJet.h"
#ifdef ALICEINTERFACE
#include <AliJetEventParticles.h>
#endif

class AliTkConeJetEvent : public TObject {
 public:
  AliTkConeJetEvent();
  ~AliTkConeJetEvent();
  
  void addJet(AliTkConeJet *jet);
  void sortJets(){if(fNJets>1) fJets->Sort(fNJets);}
  void Clear(Option_t *option="");

  Int_t getNJets()        const { return fNJets; }
  TClonesArray *getJets() const { return fJets; }
#ifdef ALICEINTERFACE
  void setJetParticles(const AliJetEventParticles* p) { 
    if(fParticles) delete fParticles;
    if(!p) {
      fParticles=new AliJetEventParticles(0);
      fParticles->SetHeader(*new TString("--- no event here, use aliev instead ---"));
    } else
      fParticles=new AliJetEventParticles(*p);
  }
  AliJetEventParticles* getJetParticles() const { return fParticles; }
#endif

  void Print(Option_t *) const {
    cout << "AliTkConeJetEvent " << fNJets << endl;
    for(Int_t i=0;i<fNJets;i++)
      cout << i <<": " << *(AliTkConeJet*)fJets->At(i) << endl;
  }

  void setRadius(Float_t r) {fRadius=r;}
  void setPtCut(Float_t p)  {fPtCut=p;}
  void setEtCut(Float_t p)  {fEtCut=p;}
  void setDesc(TString &s)  {fDesc=s;}

  Float_t getRadius() const {return fRadius;}
  Float_t getPtCut()  const {return fPtCut;}
  Float_t getEtCut()  const {return fEtCut;}
  TString getDesc()   const {return fDesc;}

 private:
  Int_t fNJets;
  TClonesArray *fJets; //->
#ifdef ALICEINTERFACE
  AliJetEventParticles *fParticles; //->
#endif

  TString fDesc;   //description to remember event
  Float_t fRadius; //radius used in the finder
  Float_t fPtCut;  //pT cut used on original event
  Float_t fEtCut;  //jet Et cut used in the finder

  ClassDef(AliTkConeJetEvent,3)
};
#endif
