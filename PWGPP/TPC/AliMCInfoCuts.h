#ifndef ALIMCINFOCUTS_H
#define ALIMCINFOCUTS_H

//------------------------------------------------------------------------------
// Class to keep selection cuts for MC tracks. 
// 
// Author: J.Otwinowski 04/02/2008 
//------------------------------------------------------------------------------

#include <TPDGCode.h>
#include "AliAnalysisCuts.h"

class TArrayI;

class AliMCInfoCuts : public AliAnalysisCuts
{
public:
  AliMCInfoCuts(const Char_t* name ="AliMCInfoCuts", const Char_t *title ="");
  virtual ~AliMCInfoCuts(); 
 
  // setters 
  void SetMinRowsWithDigits(const Int_t min=0) {fMinRowsWithDigits = min;}
  void SetMaxR(const Float_t max=1e99)         {fMaxR = max;}
  void SetMaxVz(const Float_t max=1e99)        {fMaxVz = max;}
  void SetRangeTPCSignal(const Float_t min=0, const Float_t max=1e99)  {fMinTPCSignal = min; fMaxTPCSignal = max;}
  void SetMinTrackLength(const Int_t min=0)    {fMinTrackLength = min;}

  // getters 
  Int_t GetMinRowsWithDigits() const {return fMinRowsWithDigits;}
  Float_t GetMaxR()            const {return fMaxR;}
  Float_t GetMaxVz()           const {return fMaxVz;}
  Float_t GetMinTPCSignal()    const {return fMinTPCSignal;}
  Float_t GetMaxTPCSignal()    const {return fMaxTPCSignal;}
  Float_t GetMinTrackLength()    const {return fMinTrackLength;}

  Float_t GetEP()  const       {return ep;}
  Float_t GetEM()  const       {return em;}
  Float_t GetMuP() const       {return mup;}
  Float_t GetMuM() const       {return mum;}
  Float_t GetPiP() const       {return pip;}
  Float_t GetPiM() const       {return pim;}
  Float_t GetKP()  const       {return kp;}
  Float_t GetKM()  const       {return km;}
  Float_t GetProt() const      {return prot;}
  Float_t GetProtBar() const   {return protbar;}
 
  // cuts init function
  void InitME();

  // check MC tracks
  virtual Bool_t IsSelected(TObject *) {return kTRUE;}
  virtual Bool_t IsSelected(TList *) {return kTRUE;}

  // add particle to array 
  void AddPdgParticle(Int_t idx=-1, Int_t pdgcode=0) const;

  // check particle in array 
  Bool_t IsPdgParticle(Int_t pdgcode=0) const;

  // check particle in array 
  Bool_t IsPosPdgParticle(Int_t pdgcode=0) const;

  // Merge output objects (needed by PROOF) 
  virtual Long64_t Merge(TCollection* list);

private:
  Int_t fMinRowsWithDigits; // min. number of TPC raws with digits  
  Float_t fMaxR;            // max. R distance from MC vertex 
  Float_t fMaxVz;           // max. Z distance from MC vertex
  Float_t fMinTPCSignal;    // min. TPC Signal calculated from Bethe Bloch formula
  Float_t fMaxTPCSignal;    // max. TPC Signal calculated from Bethe Bloch formula
  Float_t fMinTrackLength;  // min. TPC minimum track length
 
  TArrayI* aTrackParticles; // array of tracked particles 

  // PDG tracked particles (later added to aTrackParticles)
  enum enumData {
    kNParticles = 10, // number of particles below
    ep = kPositron,
    em = kElectron,
    mup = kMuonPlus,
    mum = kMuonMinus,
    pip = kPiPlus,
    pim = kPiMinus,
    kp = kKPlus,
    km = kKMinus,
    prot = kProton,
    protbar = kProtonBar
  };

  AliMCInfoCuts(const AliMCInfoCuts&); // not implemented
  AliMCInfoCuts& operator=(const AliMCInfoCuts&); // not implemented

  ClassDef(AliMCInfoCuts, 1)
};

#endif // ALIMCINFOCUTS_H
