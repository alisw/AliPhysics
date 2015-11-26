#ifndef ALIAODTZERO_H
#define ALIAODTZERO_H

//-------------------------------------------------------------------------
//     Container class for AOD TZERO data
//     Author: Filip Krizek
//     filip.krizek@cern.ch 23/02/2012
//-------------------------------------------------------------------------

#include <TObject.h>
#include <TBits.h>

class AliAODTZERO : public TObject 
{
public:
  AliAODTZERO();
  AliAODTZERO(const AliAODTZERO& source);
  AliAODTZERO &operator=(const AliAODTZERO& source);

  virtual ~AliAODTZERO() {};

  // Getters 
  //1st
  Double32_t GetT0TOF(Int_t i)  const {return fT0TOF[i];}
  const Double32_t * GetT0TOF() const {return fT0TOF;}
  //best
  Double32_t GetT0TOFbest(Int_t i)  const {return fT0TOFbest[i];}
  const Double32_t * GetT0TOFbest() const {return fT0TOFbest;}
 
  Bool_t GetBackgroundFlag() const {return fBackground;}
  Bool_t GetPileupFlag()     const {return fPileup;}
  Bool_t GetSatellite()      const {return fSattelite;}
  
  Float_t GetT0VertexRaw()      const {return fT0VertexRaw;}
  Double32_t GetT0zVertex()      const {return fT0zVertex;}

  Float_t GetAmp(Int_t pmt)  const {return fT0Amp[pmt];}
  
  //Setters
  void SetT0TOF(Int_t icase, Double32_t time) { fT0TOF[icase] = time;}
  void SetT0TOFbest(Int_t icase, Double32_t time) { fT0TOFbest[icase] = time;}
   
  void SetBackgroundFlag(Bool_t back = false) {fBackground = back;}
  void SetPileupFlag(Bool_t back = false) {fPileup  = back;}
  void SetSatelliteFlag(Bool_t sat = false) { fSattelite = sat;}
  
  void SetT0VertexRaw(Float_t vtx) { fT0VertexRaw = vtx;}
  void SetT0zVertex(Double32_t z) {fT0zVertex = z;}
  void SetAmp(Int_t pmt, Float_t amp) {fT0Amp[pmt]=amp;}
  //pile up bits
  void SetPileupBits(TBits pileup) {fPileupBits=pileup;}
  TBits GetT0PileupBits() const {return fPileupBits;}
     
  
protected:
  Double32_t   fT0TOF[3];    // interaction time in ps with 1st time( A&C, A, C)
  Bool_t       fPileup;      // pile-up flag
  Bool_t       fSattelite;   // sattelite flag
  Bool_t       fBackground;  // sattelite flag
  Double32_t   fT0TOFbest[3];// interaction time in ps ( A&C, A, C) with best time
  Float_t      fT0VertexRaw; // raw T0 vertex without any cuts 
  Double32_t   fT0zVertex;    // reconstructed T0 vertex
  Float_t fT0Amp[26];          //amplitude on PMTs and MPD
  TBits fPileupBits;     //BC number

  ClassDef(AliAODTZERO,5)
};

#endif
