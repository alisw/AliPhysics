#ifndef ALITRDQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the TRD data in simulated data
//
//*-- Sylwester Radomski
//////////////////////////////////////////////////////////////////////////////

#include <TTree.h> 
#include "AliAnalysisTask.h"  

class AliESD; 
class TH1D; 
class TH2D;

class AliTRDQATask : public AliAnalysisTask {

public:
  AliTRDQATask(const char *name);
  virtual ~AliTRDQATask() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void Init(Option_t * opt = ""); 
  virtual void Terminate(Option_t * opt = "");

private:

  int  GetSector(double alpha);
  int  CheckSector(int sector);
  void CalculateEff();
  void DrawESD() ; 
  void DrawGeoESD() ; 
  void DrawConvESD() ; 
  void DrawPidESD() ; 

  TTree   * fChain;             //!pointer to the analyzed TTree or TChain
  AliESD  * fESD;               //! Declaration of leave types

  TObjArray * fOutputContainer; //! output data container

  // options
  int fConfSM;
  
  // Histograms
  TH1D *fNTracks;
  TH1D *fEventSize;
  TH1D *fTrackStatus;

  TH1D *fParIn;
  TH1D *fParOut;
  TH1D *fKinkIndex;
   
  // TPC clusters histograms
  //TH1D *fTpcNCls;
  //TH1D *fTpcFCls;
  //TH1D *fTpcRCls; 
  
  // last measurement X plane
  TH1D *fXIn;
  TH1D *fXOut;
  
  // sector
  TH1D *fAlpha[4];
  TH1D *fSectorTRD;
 
  //static const int knbits = 5;
  
  // track parameters
  TH1D *fPt[6];
  TH1D *fTheta[6];
  TH1D *fSigmaY[6]; 
  TH1D *fChi2[6];
  TH2D *fPlaneYZ[6];

  TH1D *fEffPt[4];

  // track features
  TH1D *fClustersTRD[3];

  // for good refitted tracks only
  TH1D *fTime;
  TH1D *fBudget;
  TH1D *fQuality;
  TH1D *fSignal;

  // PID for TPC and TRD  
  TH2D *fTrdSigMom;
  TH2D *fTpcSigMom;
  
  TH1D *fTrdPID[6];
  TH2D *fTrdSigMomPID[6];
  
  TH1D *fTpcPID[6];
  TH2D *fTpcSigMomPID[6];
      
  
  ClassDef(AliTRDQATask, 0); // a TRD analysis task 
};
#endif // ALITRDQATASK_H
