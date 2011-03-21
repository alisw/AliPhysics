#ifndef ALITPCCALIBCOSMIC_H
#define ALITPCCALIBCOSMIC_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliTPCcalibBase.h"
class TH2F;
class TH1F;
class TList;
class AliESDEvent;
class AliESDtrack;
class THnSparse;

class AliTPCcalibCosmic:public AliTPCcalibBase {
public:
  AliTPCcalibCosmic(); 
  AliTPCcalibCosmic(const Text_t *name, const Text_t *title);
  virtual ~AliTPCcalibCosmic();
  
  virtual void      Process(AliESDEvent *event);
  virtual Long64_t  Merge(TCollection *const li);
  void              Add(const AliTPCcalibCosmic* cosmic);
  //
  //
  void              Init();
  void              FindPairs(const AliESDEvent *event);
  void              FindCosmicPairs(const AliESDEvent * event);

  Bool_t            IsPair(AliExternalTrackParam *tr0, AliExternalTrackParam *tr1) const;
  static void       CalculateBetheParams(TH2F *hist, Double_t * initialParam);
  static Double_t   CalculateMIPvalue(TH1F * hist);
  AliExternalTrackParam *MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1);
  AliExternalTrackParam *MakeCombinedTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1);

  void UpdateTrack(AliExternalTrackParam &track0, const AliExternalTrackParam &track1);
  //
  void FillHistoPerformance(const AliExternalTrackParam *par0, const AliExternalTrackParam *par1, const AliExternalTrackParam *inner0, const AliExternalTrackParam *inner1, AliTPCseed *seed0,  AliTPCseed *seed1, const AliExternalTrackParam *param0Combined, Int_t cross);
  static void MakeFitTree(TTree * treeInput, TTreeSRedirector *pcstream, const TObjArray * corrArray, Int_t step, Int_t run);
  TTree * GetCosmicTree() const {return fCosmicTree;}
  //
  TH1F   *          GetHistNTracks() const {return fHistNTracks;};
  TH1F   *          GetHistClusters() const {return fClusters;};
  TH2F   *          GetHistAcorde()const {return fModules;};
  TH1F   *          GetHistPt() const {return fHistPt;};
  TH2F   *          GetHistDeDx() const {return fDeDx;};
  TH1F   *          GetHistMIP() const {return fDeDxMIP;};
  //
  Double_t          GetMIPvalue()const {return fMIPvalue;};
  //
  static void       BinLogX(TH1 *const h);   // method for correct histogram binning
  static void       BinLogX(THnSparse *const h, Int_t axisDim);   // method for correct histogram binning

  void     Process(AliESDtrack *const track, Int_t runNo=-1) {AliTPCcalibBase::Process(track,runNo);};
  void     Process(AliTPCseed *const track)  {return AliTPCcalibBase::Process(track);}
  virtual void  Terminate();
  static Double_t GetDeltaTime(Double_t rmin0, Double_t rmax0, Double_t rmin1, Double_t rmax1, Double_t tmin0, Double_t tmax0, Double_t tmin1, Double_t tmax1, Double_t dcaR, TVectorD& vectorDT);
public:  
  //
  // Performance histograms
  //
  THnSparse   *fHistoDelta[6];  // histograms of tracking performance delta
  THnSparse   *fHistoPull[6];   // histograms of tracking performance pull
  THnSparse   *fHistodEdxMax[4];   // histograms of dEdx perfomance - max charge
  THnSparse   *fHistodEdxTot[4];   // histograms of dEdx perfomance - tot charge
  static void AddTree(TTree* treeOutput, TTree * treeInput);
private:
  
  void              FillAcordeHist(AliESDtrack *upperTrack);

  

  TH1F  *fHistNTracks;            //  histogram showing number of ESD tracks per event
  TH1F  *fClusters;               //  histogram showing the number of clusters per track
  TH2F  *fModules;                //  2d histogram of tracks which are propagated to the ACORDE scintillator array
  TH1F  *fHistPt;                 //  Pt histogram of reconstructed tracks
  TH2F  *fDeDx;                   //  dEdx spectrum showing the different particle types
  TH1F  *fDeDxMIP;                //  TPC signal close to the MIP region of muons 0.4 < p < 0.45 GeV

  Double_t fMIPvalue;             //  MIP value calculated via a fit to fDeDxMIP
  //
  // cuts
  //
  Float_t fCutMaxD;     // maximal distance in rfi ditection
  Float_t fCutMaxDz;     // maximal distance in z ditection
  Float_t fCutTheta;    // maximal distance in theta ditection
  Float_t fCutMinDir;   // direction vector products

  TTree  *fCosmicTree;  // tree with the cosmic tracks
  AliTPCcalibCosmic(const AliTPCcalibCosmic&); 
  AliTPCcalibCosmic& operator=(const AliTPCcalibCosmic&); 

  ClassDef(AliTPCcalibCosmic, 3); 
};

#endif

