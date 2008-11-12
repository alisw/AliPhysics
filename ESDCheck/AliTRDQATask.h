#ifndef ALITRDQATASK_H
#define ALITRDQATASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// An analysis task to check the TRD data in simulated data
// Starting from ESD
// Producing Histograms and plots
// Part of an analysis Train
//*-- Sylwester Radomski
//////////////////////////////////////////////////////////////////////////////

#include "AliAnalysisTask.h"  

class TTree; 
class AliESD; 
class TH1D; 
class TH2D;

class AliTRDQATask : public AliAnalysisTask {

public:
  AliTRDQATask(const char *name);
  virtual ~AliTRDQATask() {}
   
  virtual void Exec(Option_t * opt = "");
  virtual void ConnectInputData(Option_t *);
  virtual void CreateOutputObjects();
  virtual void Terminate(Option_t * opt = "");

private:

  int  GetSector(const double alpha) const;
  int  CheckSector(const int sector) const;
  void CalculateEff();
  void DrawESD() const ; 
  void DrawGeoESD() const ; 
  void DrawConvESD() const  ; 
  void DrawPidESD() const ; 

  TTree   * fChain;             //!pointer to the analyzed TTree or TChain
  AliESD  * fESD;               //! Declaration of leave types

  TObjArray * fOutputContainer; //! output data container

  // options
  int fConfSM; //!Super Module Configuration
  
  // Histograms
  TH1D *fNTracks;     // Number of tracks
  TH1D *fEventSize;   // Event size
  TH1D *fTrackStatus; // Status of tracks

  TH1D *fParIn;       // Par In
  TH1D *fParOut;      // Par out
  TH1D *fKinkIndex;   // Kink Index
   
  // TPC clusters histograms
  //TH1D *fTpcNCls;
  //TH1D *fTpcFCls;
  //TH1D *fTpcRCls; 
  
  // last measurement X plane
  TH1D *fXIn;        // input Xplane
  TH1D *fXOut;       // output Xplane
  
  // sector
  TH1D *fAlpha[4];   // alpha sectors
  TH1D *fSectorTRD;  // TRD sectors
 
  //static const int knbits = 5;
  
  // track parameters
  TH1D *fPt[6];      // Transverse momentum
  TH1D *fTheta[6];   // Theta distribution
  TH1D *fSigmaY[6];  // Sigma Y
  TH1D *fChi2[6];    // Chi 2
  TH2D *fPlaneYZ[6]; // YZ Plane

  TH1D *fEffPt[4];   // Eff transverse momentum

  // track features
  TH1D *fClustersTRD[3]; // Clusters

  // for good refitted tracks only
  TH1D *fTime;           // time
  TH1D *fBudget;         // Budget
  TH1D *fQuality;        // Quality
  TH1D *fSignal;         // Signal 

  // PID for TPC and TRD  
  TH2D *fTrdSigMom;      // Sig TRD
  TH2D *fTpcSigMom;      // Sig TPC
  
  TH1D *fTrdPID[6];      // Pid TRD
  TH2D *fTrdSigMomPID[6];// Pid TRD
  
  TH1D *fTpcPID[6];      // Pid TPC
  TH2D *fTpcSigMomPID[6];// Pid TPC
      
  AliTRDQATask(const AliTRDQATask&); // Not implemented
  AliTRDQATask& operator=(const AliTRDQATask&); // Not implemented
  
  ClassDef(AliTRDQATask, 0); // a TRD analysis task 
};
#endif // ALITRDQATASK_H
