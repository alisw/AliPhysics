#ifndef ALIRICHPATREC_H
#define ALIRICHPATREC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////
//   Pattern Recognition classes for set:RICH version 0       //
////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TMath.h>

#include "AliRICH.h"
#include "AliHitMap.h"


class AliRICHPatRec;

class AliRICHPatRec : public TObject {
    
 public:
  AliRICHPatRec();
  AliRICHPatRec(const char *name, const char *title);
  virtual       ~AliRICHPatRec() {}
  void   PatRec();
  Int_t   TrackParam(Int_t itr, Int_t &ich);
  //Old CERENK
  Float_t EstimationAtLimits(Float_t lim, Float_t radius, Float_t phiphot);  
  //Old REC_ETAPHOT
  Float_t PhotonCerenkovAngle();
  //Old GIME_EMISSPOINT
  void EmissionPoint();
  //Old ITER_CUT
  void PhotonSelection(Int_t track, Int_t &nphot, Float_t &thetamean);
  //Old BKG_SUBTRACT
  void BackgroundEstimation();
  //Old FLAG_PHOTONS
  void FlagPhotons(Int_t track, Float_t theta);
  //Old NEWINBAND
  Int_t PhotonInBand();
  //Old RADII
  Float_t DistanceFromMip(Float_t nf,Float_t nq,
				    Float_t Em,Float_t th, Float_t ph, Float_t pointsOnCathode[3]);
  //Old GIME_PHI
  Float_t PhiPad();
  //Old THREECOORD
  //void CoordSphere(Float_t r, Float_t theta, Float_t phi, Float_t *x);
  //Old ANGT
  Float_t SnellAngle(Float_t n1, Float_t n2, Float_t theta1);
  //Old TETCER
  Float_t CherenkovAngle(Float_t n, Float_t beta);
  // Old hough_filtering
  void HoughFiltering(float HCS[]);
  // Old hough_analysis
  void HoughResponse();
  //new
  Float_t BetaCerenkov(Float_t n, Float_t theta);  
  //new
  //void CerenkovRingDrawing();

 public:
  Int_t   fNumEtaPhotons;                 // Number of photons
  Float_t fEtaPhotons[1000];              // Cerenkov angle each photon
  Float_t fWeightPhotons[1000];           // weight for each photon
  Float_t fThetaCerenkov;                 // Theta angle for photon
  Float_t fThetaPeakPos;                  // Peak position

  Float_t fDTheta;                        //Step for sliding window
  Float_t fWindowWidth;                   //Hough width of sliding window

private:

  Float_t fRw,fQw,fTgap;                  //  

  Float_t fTrackLoc[3];                   // Position of track
  Float_t fTrackTheta;                    // Theta angle for track (dip angle)
  Float_t fTrackPhi;                      // Phi angle
  Float_t fTrackMom;                      // Track momentum
  Float_t fXpad;                          // x-position of pad hit
  Float_t fYpad;                          // y-position of pad hit
  Int_t   fQpad;                          // Charge deposited

  Float_t fXshift,fYshift;                // Shifts in x and y
  Float_t fEmissPoint;                    // Emission point for Cerenkov photons
  Float_t fCerenkovAnglePad;              // Cerenkov angle of single pad
  Float_t fOnCathode[3];                 // Point of exit from photocathode 


  ClassDef(AliRICHPatRec,1)                //Pat Rec module for :RICH version 0
      
};


	
	
#endif
