#ifndef ALICDFJETFINDER_H
#define ALICDFJETFINDER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliJetFinder.h"
#include <math.h>

//
//  Definition of constants, structures and functions
//

using namespace std ;

const Double_t pi = TMath::Pi();

struct varContainer // container for Particle and Jets
  { // variables of container struct
  Double_t  pt; Double_t eta; Double_t phi;
  Int_t njet; // if jets are stored in varContainer njet is multiplicity of jet
              // if particles are stored , njet is index number of jet
  } ;



class AliCdfJetHeader;

class AliCdfJetFinder : public AliJetFinder
  {
  public:

    AliCdfJetFinder();
    virtual ~AliCdfJetFinder();

    void           CreateOutputObjects(TList *histos);
    void           FindJets();
    void           InitData();
    void           FindCones();
    void           ComputeConesWeight();
    void           WriteJets();
    void           AnalizeJets();
    void           Clean();

    virtual void   FinishRun();

    Double_t DeltaPhiNorm (Double_t dphi)
      {
      if ( dphi < - pi ) { dphi = -dphi - 2.0 * pi ; }
      if ( dphi >   pi ) { dphi = -dphi + 2.0 * pi ; }
      return dphi;
      }

    inline long double Distance (long double x, long double y)
      { return TMath::Sqrt ( pow(x,2) + pow(y,2) ); }



  protected:
    AliCdfJetFinder ( const AliCdfJetFinder& jf );
    AliCdfJetFinder& operator = ( const AliCdfJetFinder& jf );

    TList         *fHistos;    // List of histograms

    Bool_t        fDebug;   //  enables debugging

    Bool_t fFromAod;
    Bool_t fAODwrite;   // write jets to AOD
    Bool_t fAODtracksWrite;  // write jet tracks to AOD
    TRefArray *fRefArr ; // pointer to references array of tracks from AOD


    Int_t         fNJets;     // counter of number of jets
    Int_t         fNPart;     // number of particles in event

    Double_t      fRadius;

    Int_t fMinJetParticles;
    Double_t fJetPtCut;

    varContainer **fVectParticle; // container for Particles
    varContainer **fVectJet;      // container for Jets

    Double_t *fPtArray;  // momentum array
    Int_t   *fIdxArray;  // index array of sorted pts






    ClassDef ( AliCdfJetFinder, 1 )
  };//
#endif

