#ifndef ALIDIELECTRONTRACKROTATOR_H
#define ALIDIELECTRONTRACKROTATOR_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

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

#include <AliKFParticle.h>

class TObjArray;
class AliVTrack;
class AliVEvent;

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

  void SetEvent(AliVEvent * const ev)   { fEvent = ev;           }
  void SetPdgLegs(Int_t pdfLeg1, Int_t pdfLeg2) { fPdgLeg1=pdfLeg1; fPdgLeg2=pdfLeg2; }

  const AliKFParticle& GetKFTrackP() const {return fTrackP;}
  const AliKFParticle& GetKFTrackN() const {return fTrackN;}

  AliVTrack* GetVTrackP() const {return fVTrackP;}
  AliVTrack* GetVTrackN() const {return fVTrackN;}
  
private:
  UInt_t   fIterations;             // number of iterations

  ERotationType fRotationType;      // which track to rotate
  
  Double_t fStartAnglePhi;          // starting angle for rotation
  Double_t fConeAnglePhi;           // opening angle in phi for multiple rotation

  const TObjArray *fkArrTracksP;    //! array of positive tracks
  const TObjArray *fkArrTracksN;    //! array of negative tracks

  UInt_t   fCurrentIteration;       //! current iteration step
  Int_t    fCurrentTackP;           //! current positive track in array
  Int_t    fCurrentTackN;           //! current negative track in array

  AliVEvent *fEvent;                //! current event
  
  AliKFParticle fTrackP;            //! Positive track
  AliKFParticle fTrackN;            //! Negative track
  
  AliVTrack *fVTrackP;              //! Positive track
  AliVTrack *fVTrackN;              //! Negative track
  
  Int_t fPdgLeg1;                   //! pdg code leg1
  Int_t fPdgLeg2;                   //! pdg code leg2
  

  Bool_t RotateTracks();
  
  AliDielectronTrackRotator(const AliDielectronTrackRotator &c);
  AliDielectronTrackRotator &operator=(const AliDielectronTrackRotator &c);

  
  ClassDef(AliDielectronTrackRotator,1)         // Dielectron TrackRotator
};



#endif
