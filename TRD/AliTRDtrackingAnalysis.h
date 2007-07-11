#ifndef AliTRDtrackingAnalysis_H
#define AliTRDtrackingAnalysis_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////////////
//                                                                //
// Fills a set of QA histograms to check the correctness of       //
// the TRD reconstruction                                         // 
//                                                                //
////////////////////////////////////////////////////////////////////

#include "TObject.h"

class TH1D;
class TH2D;
class TTree;
class TObjArray;
class TGraphErrors;

class AliRunLoader;
class AliTRDgeometry;
class AliESDEvent;
class AliTRDcluster;
class AliTRDtracker;

class AliTRDtrackingAnalysis : public TObject {

 public:
  
  AliTRDtrackingAnalysis();
  AliTRDtrackingAnalysis(const AliTRDtrackingAnalysis &t);
  virtual ~AliTRDtrackingAnalysis() {}
  AliTRDtrackingAnalysis  &operator=(const AliTRDtrackingAnalysis &/*t*/) { return *this; }
    
  void SetPath(const char *path) {fPath = path;}

  void DrawResolutionPt(int startEvent, int stopEvent);  
  void DrawRecPointResolution(int startEvent, int stopEvent);
  //void DrawTrackletResolution(int startEvent, int stopEvent);

 protected:

  void CheckFiles();
  void LoadRecPointsFile();
  void LoadRefs();
  Int_t GetReference(Int_t label); 
  Int_t GetMCPosition(Int_t label, Double_t x, Double_t &Y, Double_t &Z, Double_t &tgphi);

  Int_t GetPhiBin(Double_t phi) const;
  Double_t GetPhi(Int_t bin) const;

  const char *fPath;              // Path to data directory
  
  TObjArray *fRefTPC;             // TPC track references
  TObjArray *fRefTRD;             // TRD track references
  Int_t fLabels[100000];          // Track lables

  AliRunLoader *fLoader;          // AliRunLoader instance
  TTree  *fEsdTree;               // ESD tree
  AliESDEvent *fESD;                   // ESD

  AliTRDtracker *fTracker;        // TRD tracker instance

  // histograms 
  TH1D *fDeltaPt;                 // Histogram for the pt resolution
  TH1D *fDeltaZ;                  // Histogram for the cluster z deviation
  TH1D *fDeltaX;                  // Histogram for the cluster x deviation
  TH1D *fDeltaYPos;               // Histogram for the cluster y deviation (positives)
  TH1D *fDeltaYNeg;               // Histogram for the cluster y deviation (negatives)

  TH1D *fNPoints;                 // Histogram for the number of points
  TH1D *fNGood;                   // Histogram for the number of good points
  
  TH2D *fRefSpace;                // Histogram for reference space 

  AliTRDgeometry *fGeo;           // TRD geometry

  TH1D *fClY2;                    // Histogram for cluster studies Y
  TH1D *fClY3;                    // Histogram for cluster studies Y

  TH1D *fTgPhi;                   // Histogram for tangens(phi)
  TH1D *fClYTgPhi[12];            // Histogram cluster Y tangen phi

  TGraphErrors *fGrResTgPhi;      // Graph resolution tangens phi
  TGraphErrors *fGrMeanTgPhi;     // Graph mean tangens phi

  //TH1D *fPullY2;
  //TH1D *fPullY3;

  TH1D *fTrklY;                   // QA histogram
  TH1D *fTrklZ;                   // QA histogram

  TH1D *fClZ;                     // QA histogram
  TH2D *fClZZ;                    // QA histogram
  TH2D *fClYY;                    // QA histogram
  TH2D *fClYX;                    // QA histogram
  TH1D *fNLabels;                 // QA histogram
  TH1D *fTestBits;                // QA histogram
  TH1D *fRefDx;                   // QA histogram

  TH2D *fClZXref;                 // QA histogram
  TH2D *fClZXcl;                  // QA histogram

  TH2D *fClPos;                   // QA histogram
 
  ClassDef(AliTRDtrackingAnalysis,1)            // qa for Digits
};

#endif


