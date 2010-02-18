//-*- Mode: C++ -*-
#ifndef ALIGAMMACONVERSIONBGHANDLER_H
#define ALIGAMMACONVERSIONBGHANDLER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

////////////////////////////////////////////////
//--------------------------------------------- 
// Class  for handling of background calculation
//---------------------------------------------
////////////////////////////////////////////////

#include <vector>


// --- ROOT system ---
#include <TObject.h> 
#include "AliKFParticle.h"
#include "TClonesArray.h"

#if __GNUC__ >= 3
using namespace std;
#endif

typedef vector<AliKFParticle*> AliGammaConversionKFVector;

class AliGammaConversionBGHandler : public TObject {

 public: 
  
  typedef vector<AliGammaConversionKFVector> AliGammaConversionBGEventVector;
  typedef vector<AliGammaConversionBGEventVector> AliGammaConversionMultipicityVector;
  typedef vector<AliGammaConversionMultipicityVector> AliGammaConversionBGVector;

  AliGammaConversionBGHandler();                                                        //constructor
  AliGammaConversionBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t fNEvents);    //constructor
  AliGammaConversionBGHandler(const AliGammaConversionBGHandler & g);                   //copy constructor
  AliGammaConversionBGHandler & operator = (const AliGammaConversionBGHandler & g);     //assignment operator
  virtual ~AliGammaConversionBGHandler();                                               //virtual destructor

  void Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray);

  Int_t GetZBinIndex(Double_t z) const;

  Int_t GetMultiplicityBinIndex(Int_t mult) const;

  void AddEvent(TClonesArray * const eventGammas, Double_t zvalue, Int_t multiplicity);

  Int_t GetNBGEvents()const {return fNEvents;}

  AliGammaConversionKFVector* GetBGGoodV0s(Int_t event, Double_t zvalue, Int_t multiplicity);

  void PrintBGArray();

 private:

  Int_t fNEvents; // number of events
  Int_t ** fBGEventCounter; // bg counter
  
  Int_t fNBinsZ; //n z bins
  Int_t fNBinsMultiplicity; //n bins multiplicity
  Double_t *fBinLimitsArrayZ; //bin limits z array
  Double_t *fBinLimitsArrayMultiplicity; //bin limit multiplicity array
  AliGammaConversionBGVector fBGEvents; //gackground events

  ClassDef(AliGammaConversionBGHandler,0)
};
#endif
