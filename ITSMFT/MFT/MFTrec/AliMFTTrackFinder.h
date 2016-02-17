/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//====================================================================================================================================================
//
//      Class for finding Tracks from Cluster of the ALICE Muon Forward Tracker
//
//      Contact author: raphael.tieulent@cern.ch
//
//====================================================================================================================================================

#ifndef AliMFTTrackFinder_H
#define AliMFTTrackFinder_H

#include "TObject.h"
#include "TClonesArray.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"

#include "AliMFTConstants.h"
#include "AliLog.h"

class AliMFTCATrack;
class AliMFTCACell;
class AliMFTCAHit;
class AliMFTCALayer;
class AliMFTCARoad;
class TF1;

//============================================================================================

class AliMFTTrackFinder : public TObject {
  
private:
  static const Int_t fNDetMax = 10;
public:
  
  AliMFTTrackFinder();
  ~AliMFTTrackFinder();
  
  void Init(Char_t *parfile);
  void ReadParam(Char_t *parfile = "param.txt");
  void SetDebug(Int_t debug) { fDebug = debug; }
  
  virtual void Clear(Option_t *);
  
  AliMFTCATrack *AddTrack(Int_t gid);
  AliMFTCATrack *AddTrack(Int_t gid, const AliMFTCATrack& track);
  
  const Int_t GetNlayers() { return fNlayers; }
  const Int_t GetNtracks() { return fNtracks; }
  const Int_t GetNcells()  { return fCellGID; }
  
  void CreateCellsOld(Bool_t calcVertex = kFALSE);
  void CreateCells(Bool_t calcVertex = kFALSE);
  void CreateGapCells();
  void ClearCells();
  void ResetCells();
  
  AliMFTCACell *GetCellByGID(Int_t gid);
  
  void RunForward();
  void RunBackward();
  void LoadClusters( TClonesArray *clusterArrayFront[AliMFTConstants::fNMaxPlanes],  TClonesArray *clusterArrayBack[AliMFTConstants::fNMaxPlanes] );
  
  Bool_t RuleSelect(AliMFTCACell *cellL, AliMFTCACell *cellR);
  Bool_t RuleSelectCell(AliMFTCACell *cell);
  Bool_t RuleSelectCell(Double_t *h1, Double_t *h2, Int_t iL1, TF1 *f = 0, Bool_t acalc = kFALSE);
  Bool_t RuleSelect2LayersGap(Int_t iL1, Int_t iL2, Double_t *hit1, Double_t *hit2);
  
  AliMFTCALayer *GetLayer(Int_t nl) { return fLayers[nl]; }
  AliMFTCATrack *GetTrack(Int_t nt) { return (AliMFTCATrack*)fTracks->At(nt); }
  
  TList * GetHistograms(){return fHistList;};
  
  void UpdateCellStatus();
  
  void PrintAll();
  void DrawHisto();
  void DrawTracks(Double_t *pTot, Double_t *Theta);
  void FilterTracks();
  void PrintTrack(Int_t id);
  
  void CalculateVertex();
  Double_t GetZVertCalc() { return fZVertCalc; }
  void SetZVertRange(Double_t *zvr, Double_t zvd) {
    fZVertRange[0] = zvr[0]; fZVertRange[1] = zvr[1];
    fZVertDet = zvd;
  }
  
  // parameters
  const Double_t GetThick()              { return fThick; }
  const Double_t GetPixelNoise()         { return fPixelNoise; }
  const Double_t GetPlaneDetEff(Int_t i) { return fPlaneDetEff[i]; }
  const Double_t GetMBRate()             { return fMBRate; }
  const Double_t GetCMOSIntTime()        { return fCMOSIntTime; }
  const Bool_t   AddNoise()              { return fAddNoise; }
  const Bool_t   UseTF()                 { return fUseTF; }
  const Bool_t   AddQED()                { return fAddQED; }
  const Bool_t   CalcVertex()            { return fCalcVertex; }
  const Bool_t   ReadGeom()              { return fReadGeom; }
  const Char_t  *GetGeomName()           { return (fGeomName.Append(".root")).Data(); }
  
  void PrintParam();
  
  void SetCAtime(Double_t cpu, Double_t real) { fCPUTime += cpu; fRealTime += real; }
  
  const Double_t GetCPUTime() { return fCPUTime; }
  const Double_t GetRealTime() { return fRealTime; }
  
  const Int_t GetNDet() { return fNlayers; }
  
  void FindTracks();   // simple track finding
  
  void RemoveLastTrack() {
    AliInfo("Removing Last Track");
    fTracks->RemoveAt(fTracks->GetLast());
    fNtracks--;
  }
  
  void AnalyzeCells();
  const Int_t GetNDifTracks() { return fNDifTracks; }
  
  void SetPlanesZ(Double_t *z, Int_t n) {
    for (Int_t i = 0; i < n; i++) fPlanesZ[i] = z[i];
  }
  
  void BuildRoads();
  AliMFTCARoad *AddRoad();
  AliMFTCARoad *GetRoad(Int_t nr) { return (AliMFTCARoad*)fRoads->At(nr); }
  const Int_t GetNRoads() { return fNRoads; }
  
  // CA using roads
  void CreateCellsR(AliMFTCARoad *road);
  void RunForwardR(AliMFTCARoad *road, Int_t &trackGID);
  void RunBackwardR(AliMFTCARoad *road, Int_t &trackGID);
  void UpdateCellStatusR();
  Bool_t LinFit(Int_t nDet, Double_t *xcl,
                Double_t *ycl, Double_t *yerr,
                Double_t &a, Double_t &ae, Double_t &b, Double_t &be,
                Int_t skip = -1);
  Double_t GetErrX() { return fErrX; }
  Double_t GetErrY() { return fErrY; }

private:

  Float_t fXCut;                   ///< Cut in x difference; RuleSelect
  Float_t fYCut;                   ///< Cut in y difference; RuleSelect
  Float_t fACutV[fNDetMax];        ///< Cut in angle difference: for cell vertex compatibility
  Float_t fACutN[fNDetMax];        ///< Cut in angle difference: for neighbor cells compatibility
  Float_t fMaxSegAngle;            ///< Max cut of the Theta angle of segments [deg]
  
  Int_t fCellGID;                  //!<! Cell global identifier
  Int_t fMaxCellStatus;            //!<! Maximum value of a cell status after RunForward
  Int_t fNlayers;                  ///< Number of detection planes
  Int_t fNtracks;                  //!<! Number of tracks
  
  AliMFTCALayer *fLayers[fNDetMax]; //!<! Array of layers
  
  TClonesArray *fTracks;           //!<! Array of tracks
  TList * fHistList;               //!<! List of histograms
  
  TH1F *hDA[fNDetMax];             //!<! Histogram with angle between cells
  TH1F *hDAv[fNDetMax];            //!<! Histogram with angle with respect to the vertex
  TH2F *hDXY[fNDetMax];            //!<! Histogram with X,Y distance between end of cells
  TH1F *hNGoodCell;                //!<! Histogram showing numbers of good cells in the track
  TH1F *hTrackType;                //!<! Histogram track Type: 0 = good track ; 1 = hits from other track ; 2 = hits from noisy pixels
  TH1F *hAngleCells;               //!<! Angle between adjacent cells
  TH1F *hThetaCells;               //!<! Theta of the cells
  
  Bool_t   fCalcVertex;            // Calculate vertex z from cells
  Double_t fZVertCalc;             //!<! Calculated vertez z
  Double_t fZVertDet;              //!<! Vertex z given by ext detector
  Double_t fZVertRange[2];         //!<! Limits of vertex z accepted range
  Int_t    fMinTrackLength;        // Minimum requested track length in hits
  Double_t fPlaneDetEff[fNDetMax]; // Single plane efficiency
  Bool_t   fAddNoise;              // Add pixel noise
  Double_t fPixelNoise;            // Probability of pixel noise
  Bool_t   fAddQED;                // Add hits for QED processes
  Double_t fMBRate;                // Hadronic MB rate [kHz]
  Double_t fCMOSIntTime;           // CMOS Integration Time [micro seconds]
  Double_t fThick;                 // Single plane thickness in X0
  Bool_t   fReadGeom;              // Read geometry from ROOT file
  TString  fGeomName;              // Geometry file name
  Bool_t   fUseTF;                 // Use the TrackFinder
  
  Int_t fDebug;                    // Debugging level
  
  Double_t fCPUTime;               // CA CPU time
  Double_t fRealTime;              // CA real time
  
  Double_t GetCellAngleDif(AliMFTCACell *cell1, AliMFTCACell *cell2);
  Double_t GetCellInterceptX(AliMFTCACell *cell, Double_t z);
  Double_t GetCellInterceptY(AliMFTCACell *cell, Double_t z);
  
  TH2F *hHitDifXY;                 //!<!
  TH1F *hAngDifAll;                //!<!
  TH1F *hAngDifDup;                //!<!
  TH2F *hIntDifXYAll;              //!<!
  TH2F *hIntDifXYDup;              //!<!
  TH1F * hVertZ;
  TH1F * hVertZa;

  Int_t fNDifTracks;
  
  Double_t fPlanesZ[fNDetMax];     // median z values of the two faces
  Double_t fZGap;                  // distance between the two faces
  
  Int_t fNRoads;       //!<! number of built roads
  
  TClonesArray *fRoads;      //!<! Array of roads
  
  Double_t fErrX;            //!<! Error in X
  Double_t fErrY;            //!<! Error in Y

  ClassDef(AliMFTTrackFinder,2);
  
};



#endif
