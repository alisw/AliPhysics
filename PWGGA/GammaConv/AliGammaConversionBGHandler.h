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
#include "AliESDVertex.h"

#if __GNUC__ >= 3
using namespace std;
#endif

typedef vector<AliKFParticle*> AliGammaConversionKFVector;

class AliGammaConversionBGHandler : public TObject {

 public: 
  struct GammaConversionVertex
  {
    Double_t fX;
    Double_t fY;
    Double_t fZ;
  };
  typedef struct GammaConversionVertex GammaConversionVertex; //!

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

  void AddEvent(TClonesArray * const eventGammas, Double_t xvalue,Double_t yvalue,Double_t zvalue, Int_t multiplicity);
  void AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity);

  Int_t GetNBGEvents()const {return fNEvents;}

  AliGammaConversionKFVector* GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event);
  AliGammaConversionKFVector* GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity);
  void PrintBGArray();

  GammaConversionVertex * GetBGEventVertex(Int_t zbin, Int_t mbin, Int_t event){return &fBGEventVertex[zbin][mbin][event];}

  Double_t GetBGProb(Int_t z, Int_t m){return fBGProbability[z][m];}

 private:

  Int_t fNEvents; // number of events
  Int_t ** fBGEventCounter; //! bg counter
  Int_t ** fBGEventENegCounter;//! bg electron counter
  Double_t ** fBGProbability; //! prob per bin
  GammaConversionVertex *** fBGEventVertex;//! array of event vertex
  Int_t fNBinsZ; //n z bins
  Int_t fNBinsMultiplicity; //n bins multiplicity
  Double_t *fBinLimitsArrayZ;//! bin limits z array
  Double_t *fBinLimitsArrayMultiplicity;//! bin limit multiplicity array
  AliGammaConversionBGVector fBGEvents; //gackground events
  AliGammaConversionBGVector fBGEventsENeg; //background electron events
  ClassDef(AliGammaConversionBGHandler,2)
};
#endif
