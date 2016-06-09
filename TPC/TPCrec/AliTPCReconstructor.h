#ifndef ALITPCRECONSTRUCTOR_H
#define ALITPCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"
#include "AliTPCRecoParam.h"
#include "TVector.h"
#include <TString.h>

class AliTPCParam;
class AliTPCclusterer;
class AliTPCtracker;
class AliTPCAltroEmulator;
class TObjArray;
class TTreeSRedirector;

class AliTPCReconstructor: public AliReconstructor {
public:
  AliTPCReconstructor();
  virtual ~AliTPCReconstructor();
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const;
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const;

  virtual AliTracker*  CreateTracker() const;

  virtual void         FillESD(TTree* /*digitsTree*/, TTree* /*clustersTree*/, 
			       AliESDEvent* esd) const;
  virtual void         FillESD(AliRawReader* /*rawReader*/, TTree* /*clustersTree*/, 
			       AliESDEvent* esd) const
  {FillESD((TTree*)NULL,(TTree*)NULL,esd);}

  static const AliTPCRecoParam* GetRecoParam() { return dynamic_cast<const AliTPCRecoParam*>(AliReconstructor::GetRecoParam(1)); }
  virtual void                 GetPidSettings(AliESDpid *esdPID);
  //
  static void        SetPIDRespnonsePath(const char* pth) {fgPIDRespnonsePath = pth;}
  static const char* GetPIDRespnonsePath() {return fgPIDRespnonsePath.Data();}  
  //
  static Double_t GetCtgRange()     { return GetRecoParam()->GetCtgRange();}
  static Double_t GetMaxSnpTracker(){ return GetRecoParam()->GetMaxSnpTracker();}
  static Double_t GetMaxSnpTrack()  { return GetRecoParam()->GetMaxSnpTrack();}

  static Int_t StreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}
  static void  SetAltroEmulator(AliTPCAltroEmulator *altro) { fAltroEmulator=altro;}
  static AliTPCAltroEmulator *  GetAltroEmulator() { return fAltroEmulator;}
  static TTreeSRedirector    *GetDebugStreamer(){return fgDebugStreamer;}
  static void SetDebugStreamer(TTreeSRedirector    *debugStreamer){fgDebugStreamer=debugStreamer;}
  void ParseOptions(AliTPCtracker* tracker) const;
  static  const Double_t * GetSystematicError()  { return (fSystematicErrors)? fSystematicErrors->GetMatrixArray():0;}
  static  const Double_t * GetSystematicErrorCluster() { return (fSystematicErrorClusters) ? fSystematicErrorClusters->GetMatrixArray():0;}
  static  const Double_t * GetExtendedRoads()  { return (fgExtendedRoads)? fgExtendedRoads->GetMatrixArray():0; }
  static  const Double_t * GetPrimaryDCACut()  { return (fgPrimaryDCACut)? fgPrimaryDCACut->GetMatrixArray():0; }
  static        Double_t   GetPrimaryZ2XCut()  { return fgPrimaryZ2XCut; }
  static        Double_t   GetZOutSectorCut()  { return fgZOutSectorCut; }
  static        Bool_t     GetCompactClusters() { return fgCompactClusters;}

  static  void SetSystematicError( TVectorD *vec)  { fSystematicErrors=vec;}
  static  void SetSystematicErrorCluster( TVectorD *vec );
  static  void SetExtendedRoads( TVectorD *extendedRoads ) { fgExtendedRoads=extendedRoads;}
  static  void SetPrimaryDCACut( TVectorD *dcacut );
  static  void SetPrimaryZ2XCut( double v )                { fgPrimaryZ2XCut=v; }
  static  void SetZOutSectorCut( double v )                { fgZOutSectorCut=v; }
  static  void SetCompactClusters(Bool_t v)                { fgCompactClusters=v;}

private:
  AliTPCReconstructor(const AliTPCReconstructor&); //Not implemented
  AliTPCReconstructor& operator=(const AliTPCReconstructor&); //Not implemented
  AliTPCParam*         GetTPCParam() const;
  static Int_t               fgStreamLevel; // flag for streaming      - for TPC reconstruction
  static TTreeSRedirector    *fgDebugStreamer; // pointer to the streamer
  AliTPCclusterer*           fClusterer;   // TPC clusterer
  static AliTPCAltroEmulator * fAltroEmulator;    // ALTRO emulator
  static TString             fgPIDRespnonsePath;           // path to PIDResponse
  //
  // varaibles which overwrite content of the TPCRecoParam in case of custom recosntrcution (e.g CPass0 with imperfect calibration)
  static TVectorD            * fSystematicErrors;    // systematic errors for the TPC tracks
  static TVectorD            * fSystematicErrorClusters;    // systematic errors for the TPC tracks
  static TVectorD            * fgExtendedRoads;       // extended roads for clusters
  static TVectorD            * fgPrimaryDCACut;       // only primaries passing DCAYZ cut are reconstructed
  static Double_t              fgPrimaryZ2XCut;       // cut on Z2X for fast primaries reco 
  static Double_t              fgZOutSectorCut;       // cut on Z going on other side of CE 
  static Bool_t                fgCompactClusters;     // if true, cluster coordinates will be set to 0 in clusterizer
  TObjArray *fArrSplines;                  // array of pid splines

  void SetSplinesFromOADB(const char* tmplt, AliESDpid *esdPID);
  
  ClassDef(AliTPCReconstructor, 0)   // class for the TPC reconstruction
};

#endif
