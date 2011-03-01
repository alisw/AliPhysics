#ifndef ALITRACKCOMAPRISONESD_H
#define ALITRACKCOMAPRISONESD_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////
////
////

#include "AliAnalysisTask.h"
#include "TObjArray.h"
#include "TRefArray.h"
class TGeoHMatrix;
class AliEMCALGeometry;
class AliESDEvent;
class AliESDfriend;

class AliESDtrack;
class AliESDfriendTrack;
class AliESDtrackCuts;

class AliESDCaloCells;

class AliTrackComparison;

class AliTrackComparisonESD:public AliAnalysisTask {
public:
  AliTrackComparisonESD();
  AliTrackComparisonESD(const char *name);
  virtual ~AliTrackComparisonESD();

  virtual void ConnectInputData(Option_t *option);
  virtual void CreateOutputObjects();
  virtual void Exec(Option_t *option);
  virtual void Terminate(Option_t *option);
  virtual void FinishTaskOutput();
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}
  //void         SetOcdbPath(const char *path){fOcdbPath=path;}
  //  TString      GetOcdbPath(){return fOcdbPath;}

  Bool_t SetupEvent();
  void ProcessTOF(AliESDtrack *track, AliESDfriendTrack *friendTrack, Double_t *vPos);
  void ProcessEMCAL(AliESDtrack *track, AliESDfriendTrack *friendTrack, TRefArray *clusters, Double_t *vPos);
  void ProcessHMPID(AliESDtrack *track, AliESDfriendTrack *friendTrack,Double_t *vPos);

  TObjArray *GetComparisonOutput() {return fOutput;}

  void RecalClusterPos(TRefArray *clusters, AliESDCaloCells *cells);
  void  SetCutWindow(Double_t cutX, Double_t cutY, Double_t cutZ)
  {fCutX=cutX;fCutY=cutY;fCutZ=cutZ;}

protected:
  virtual Long64_t Merge(TCollection *li);
  virtual void     Analyze();
  void             RegisterDebugOutput();
private:
  AliESDEvent *fESD;         //! current esd
  AliESDtrackCuts *fESDCuts; //! esd track cuts
  AliESDfriend *fESDfriend;  //! current esd friend
  Int_t fCurrentRun;
  TString      fDebugOutputPath; // debug output path
  //  TString      fOcdbPath;
  
  TObjArray    *fOutput;    //Output array for fEMCAL,fHMPID,fTOF
  AliTrackComparison *fEMCAL;            // EMCAL track comparison
  AliTrackComparison *fHMPID;             //HMPID track comparison
  AliTrackComparison *fTOF;              // TOF track comparison
  //
  AliEMCALGeometry *fGeom; //EMCAL geometry for position calculation
  Double_t  fCutX; 
  Double_t  fCutY;
  Double_t  fCutZ;
  TGeoHMatrix *fTransMatrix[4];

  AliTrackComparisonESD(const AliTrackComparisonESD&);
  AliTrackComparisonESD& operator=(const AliTrackComparisonESD&);


  ClassDef(AliTrackComparisonESD,1)
};

#endif
