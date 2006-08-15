#ifndef ALITRDSEED_H
#define ALITRDSEED_H   

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */ 

/* $Id$ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  The TRD track seed                                                    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TObject.h" 

class AliTRDcluster;

class AliTRDseed : public TObject {

  friend class AliTRDtracker;

 public:

  AliTRDseed(); 
  AliTRDseed(const AliTRDseed &s);
  ~AliTRDseed() {};                 

  AliTRDseed      &operator=(const AliTRDseed &/*s*/)       { return *this;   } 

  static  void     EvaluateUni(Int_t nvectors, Double_t *data, Double_t &mean, Double_t &sigma, Int_t hh);
  static  Float_t  FitRiemanTilt(AliTRDseed *seed, Bool_t error);
          void     UseClusters();
          void     Update();
          void     CookLabels();
          void     UpdateUsed();
          void     Reset();
          Bool_t   IsOK() const                             { return fN2 > 8; }

 private:

          Float_t  fTilt;               //  Tilting angle
          Float_t  fPadLength;          //  Pad length
          Float_t  fX0;                 //  X0 position
          Float_t  fX[25];              //! X position
          Float_t  fY[25];              //! Y position
          Float_t  fZ[25];              //! Z position
          Int_t    fIndexes[25];        //! Indexes
          AliTRDcluster *fClusters[25]; //! Clusters
          Bool_t   fUsable[25];         //! Indication  - usable cluster
          Float_t  fYref[2];            //  Reference y
          Float_t  fZref[2];            //  Reference z
          Float_t  fYfit[2];            //  Y fit position +derivation
          Float_t  fYfitR[2];           //  Y fit position +derivation
          Float_t  fZfit[2];            //  Z fit position
          Float_t  fZfitR[2];           //  Z fit position
          Float_t  fSigmaY;             //  "Robust" sigma in Y - constant fit
          Float_t  fSigmaY2;            //  "Robust" sigma in Y - line fit
          Float_t  fMeanz;              //  Mean vaue of z
          Float_t  fZProb;              //  Max probbable z
          Int_t    fLabels[2];          //  Labels
          Int_t    fN;                  //  Number of associated clusters
          Int_t    fN2;                 //  Number of not crossed
          Int_t    fNUsed;              //  Number of used clusters
          Int_t    fFreq;               //  Frequency
          Int_t    fNChange;            //  Change z counter
          Float_t  fMPads;              //  Mean number of pads per cluster

          Float_t  fC;                  //  Curvature
          Float_t  fCC;                 //  Curvature with constrain
          Float_t  fChi2;               //  Global chi2
          Float_t  fChi2Z;              //  Global chi2

  ClassDef(AliTRDseed,1)                //  Seed for a local TRD track

};

#endif 
