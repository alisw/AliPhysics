#ifndef ALITRDQAESDFRIENDS_H
#define ALITRDQAESDFRIENDS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

/* $Id: AliTRDqaESDFriends.h  $ */

//
// This class is a part of a package of high level QA monitoring for TRD.
// The residuals of cluster with respect to tracklets are analyzed 
// in this class. This class needs ESDfriends.root
//
// S. Radomski
// radomski@physi.uni-heidelberg.de
// March 2008
//

#include "AliAnalysisTask.h"  

class TTree; 
class AliESDEvent; 
class TH1D; 
class TH2D;
class AliExternalTrackParam;

class AliTRDqaESDFriends : public AliAnalysisTask {

public:
  AliTRDqaESDFriends();
  AliTRDqaESDFriends(const char *name);
  AliTRDqaESDFriends(AliTRDqaESDFriends& trd);
  AliTRDqaESDFriends& operator = (const AliTRDqaESDFriends& /*g*/) { return *this; };
  virtual ~AliTRDqaESDFriends() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

protected:
 
  TTree        * fChain;        //!pointer to the analyzed TTree or TChain
  AliESDEvent  * fESD;          //! Declaration of leave types

  TObjArray * fOutputContainer; //! output data container
  
  // histograms
  TH1D *fResiduals;             // residuals distribution
  TH2D *fResidualsAngle;        // diferential resisuals distribution
  //TH2D *fResidualsAngleChamber[540];   // per chamber

  ClassDef(AliTRDqaESDFriends, 0); // a TRD analysis task 
};
#endif // ALITRDQAESDFRIENDS_H
