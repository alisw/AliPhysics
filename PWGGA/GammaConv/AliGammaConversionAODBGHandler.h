//-*- Mode: C++ -*-
#ifndef ALIGAMMACONVERSIONAODBGHANDLER_H
#define ALIGAMMACONVERSIONAODBGHANDLER_H
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
#include "AliAODConversionPhoton.h"
#include "AliAODConversionMother.h"
#include "TClonesArray.h"
#include "AliESDVertex.h"

#if __GNUC__ >= 3
using namespace std;
#endif

typedef vector<AliAODConversionPhoton*> AliGammaConversionAODVector;

class AliGammaConversionAODBGHandler : public TObject {

public: 
   struct GammaConversionVertex
   {
      Double_t fX;
      Double_t fY;
      Double_t fZ;
   };
   typedef struct GammaConversionVertex GammaConversionVertex; //!

   typedef vector<AliGammaConversionAODVector> AliGammaConversionBGEventVector;
   typedef vector<AliGammaConversionBGEventVector> AliGammaConversionMultipicityVector;
   typedef vector<AliGammaConversionMultipicityVector> AliGammaConversionBGVector;

   AliGammaConversionAODBGHandler();                                                        //constructor
   AliGammaConversionAODBGHandler(UInt_t binsZ,UInt_t binsMultiplicity,UInt_t nEvents); // constructor
   AliGammaConversionAODBGHandler(UInt_t collisionSystem,UInt_t centMin,UInt_t centMax,UInt_t nEvents, Bool_t useTrackMult);
   AliGammaConversionAODBGHandler(const AliGammaConversionAODBGHandler & g);                   //copy constructor
   AliGammaConversionAODBGHandler & operator = (const AliGammaConversionAODBGHandler & g);     //assignment operator
   virtual ~AliGammaConversionAODBGHandler();                                               //virtual destructor

   void Initialize(Double_t * const zBinLimitsArray, Double_t * const multiplicityBinLimitsArray);

   Int_t GetZBinIndex(Double_t z) const;

   Int_t GetMultiplicityBinIndex(Int_t mult) const;

   void AddEvent(TList* const eventGammas, Double_t xvalue,Double_t yvalue,Double_t zvalue, Int_t multiplicity);
   void AddElectronEvent(TClonesArray* const eventENeg, Double_t zvalue, Int_t multiplicity);

   Int_t GetNBGEvents()const {return fNEvents;}

   AliGammaConversionAODVector* GetBGGoodV0s(Int_t zbin, Int_t mbin, Int_t event);
   AliGammaConversionAODVector* GetBGGoodENeg(Int_t event, Double_t zvalue, Int_t multiplicity);
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
   AliGammaConversionBGVector fBGEvents; //background events
   AliGammaConversionBGVector fBGEventsENeg; //background electron events
   ClassDef(AliGammaConversionAODBGHandler,2)
};
#endif
