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
struct AliGammaConversionBGEvent
{
  AliGammaConversionBGEvent() :
    fReconstructedGammas(),
    fChargedTrackMultiplicity(0),
    fZVertexPosition(0.)
  {
  }
  AliGammaConversionKFVector fReconstructedGammas;
  UInt_t fChargedTrackMultiplicity;
  Double_t fZVertexPosition;
}; typedef struct AliGammaConversionBGEvent AliGammaConversionBGEvent; //!	

  
  typedef vector<AliGammaConversionBGEvent> AliGammaConversionBGEventVector;
  //  typedef vector<Int_t> AliGammaConversionBGEventVector;
  typedef vector<AliGammaConversionBGEventVector> AliGammaConversionMultipicityVector;
  typedef vector<AliGammaConversionMultipicityVector> AliGammaConversionBGVector;
  
    

  AliGammaConversionBGHandler();                                                        //constructor
  AliGammaConversionBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t fNEvents);    //constructor
  AliGammaConversionBGHandler(const AliGammaConversionBGHandler & g);                   //copy constructor
  AliGammaConversionBGHandler & operator = (const AliGammaConversionBGHandler & g);     //assignment operator
  virtual ~AliGammaConversionBGHandler();                                               //virtual destructor

  void Initialize(Double_t *zBinLimitsArray, Double_t *multiplicityBinLimitsArray);

  Int_t GetZBinIndex(Double_t z);

  Int_t GetMultiplicityBinIndex(Int_t mult);

  void AddEvent(TClonesArray * eventGammas, Double_t zvalue, Int_t multiplicity);

  Int_t GetNBGEvents(){return fNEvents;}

  AliGammaConversionKFVector* GetBGGoodV0s(Int_t event, Double_t zvalue, Int_t multiplicity);

  void PrintBGArray();

 private:

  Int_t fNEvents;
  Int_t ** fBGEventCounter;
  
  Int_t fNBinsZ;
  Int_t fNBinsMultiplicity;
  Double_t *fBinLimitsArrayZ;
  Double_t *fBinLimitsArrayMultiplicity;
  AliGammaConversionBGVector fBGEvents;

  ClassDef(AliGammaConversionBGHandler,0)
};
#endif
