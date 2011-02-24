#ifndef ALIDIELECTRONTRACKROTATOR_H
#define ALIDIELECTRONTRACKROTATOR_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//#############################################################
//#                                                           # 
//#         Class AliDielectronTrackRotator                   #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni HD / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <TNamed.h>

class TObjArray;
class AliVTrack;

class AliDielectronTrackRotator : public TNamed {
public:
  enum ERotationType {kRotatePositive, kRotateNegative, kRotateBothRandom};
  
  AliDielectronTrackRotator();
  AliDielectronTrackRotator(const char*name, const char* title);

  virtual ~AliDielectronTrackRotator();

  void SetTrackArrays(const TObjArray * const arrP, const TObjArray * const arrN) {fkArrTracksP=arrP;fkArrTracksN=arrN;}
  void Reset();
  Bool_t NextCombination();

  //Setters
  void SetIterations(UInt_t niter)         { fIterations=niter;  }
  void SetRotationType(ERotationType type) { fRotationType=type; }
  void SetStartAnglePhi(Double_t phi)      { fStartAnglePhi=phi; }
  void SetConeAnglePhi(Double_t phi)       { fConeAnglePhi=phi;  }

  //Getters
  Int_t GetIterations() const           { return fIterations;    }
  ERotationType GetRotationType() const { return fRotationType;  }
  Double_t GetStartAnglePhi() const     { return fStartAnglePhi; }
  Double_t GetConeAnglePhi() const      { return fConeAnglePhi;  }


  AliVTrack* GetTrackP() const {return fTrackP;}
  AliVTrack* GetTrackN() const {return fTrackN;}

private:
  UInt_t   fIterations;             // number of iterations

  ERotationType fRotationType;      // which track to rotate
  
  Double_t fStartAnglePhi;          // starting angle for rotation
  Double_t fConeAnglePhi;           // opening angle in phi for multiple rotation

  const TObjArray *fkArrTracksP;           //! array of positive tracks
  const TObjArray *fkArrTracksN;           //! array of negative tracks

  UInt_t   fCurrentIteration;       //! current iteration step
  Int_t    fCurrentTackP;           //! current positive track in array
  Int_t    fCurrentTackN;           //! current negative track in array

  AliVTrack *fTrackP;               //! Positive track
  AliVTrack *fTrackN;               //! Negative track

  Bool_t RotateTracks();
  
  AliDielectronTrackRotator(const AliDielectronTrackRotator &c);
  AliDielectronTrackRotator &operator=(const AliDielectronTrackRotator &c);

  
  ClassDef(AliDielectronTrackRotator,1)         // Dielectron TrackRotator
};



#endif
