#ifndef ALITRDQAJPSI_H
#define ALITRDQAJPSI_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//
// This class is a part of a package of high level QA monitoring for TRD.
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliAnalysisTask.h"  

class TTree; 
class AliESDEvent; 
class AliESDtrack;
class TH1D; 
class TH2D;
class AliExternalTrackParam;
class AliKFParticle;
class TLorentzVector;

class AliTRDqaJPsi : public AliAnalysisTask {

 public:

  AliTRDqaJPsi();
  AliTRDqaJPsi(const char *name);
  AliTRDqaJPsi(AliTRDqaJPsi& trd);
  AliTRDqaJPsi& operator = (const AliTRDqaJPsi& /*g*/) { return *this; };
  virtual ~AliTRDqaJPsi() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

 private:
 
  TTree        * fChain;             //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD;               //! Declaration of leave types

  TObjArray * fOutputContainer; //! output data container
  
  // histograms

  static const Int_t knSteps = 5;   // number of analysis steps (arbitrary)

  TH1D *fStatus[knSteps];      // track status
  TH1D *fnTracks[2*knSteps];   // number of tracks
  TH1D *fPt[2*knSteps];        // transverse momentum
  TH1D *fPID[2*knSteps];       // PID LQ
  TH1D *fAngleSM[knSteps];     // difference in SM ID
  
  //TH2D *fnGoodTracks;          // correlation of the final number of Pos and Neg tracks
  TH1D *fInvMass[knSteps];       // invariant mass using different cuts
  TH1D *fInvMassVec[knSteps];    //
  TH1D *fInvMassDiff[knSteps];   //

  TH2D *fPtAngle[knSteps];

  // tracks
  AliKFParticle *tracks[1000];
  TLorentzVector *vec[1000];
  Int_t fInSample[1000][knSteps];
  Int_t fSM[1000];
  Int_t fnKFtracks; //[2];
  
  // helper functions
  void FillHist(AliESDtrack *track, Int_t step);
  TLorentzVector *CreateVector(AliESDtrack *track);


  ClassDef(AliTRDqaJPsi, 0); // a TRD analysis task 
};
#endif // ALITRDQAJPSI_H
