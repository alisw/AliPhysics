#ifndef ALIESDMUONTRACKCUTS_H
#define ALIESDMUONTRACKCUTS_H

/* $Id$ */ 

//
//  Class for handling of ESD Muon track cuts 
//  (based on ANALYSIS/AliESDtrackCuts).
//
//  The class manages some kinematic cuts. Two methods
//  can be used to figure out if an ESD Muon track survives the cuts:
//  AcceptTrack which takes a single AliESDMuonTrack as argument and
//  returns kTRUE/kFALSE or GetAcceptedTracks which takes an AliESD
//  object and returns an TObjArray (of AliESDMuonTracks) with the tracks
//  in the ESD that survived the cuts.
//

#include <TF1.h>
#include <TH2.h>
#include "AliAnalysisCuts.h"

class AliESD;
class AliESDEvent;
class AliESDMuonTrack;
class AliLog;
class TTree;

class AliESDMuonTrackCuts : public AliAnalysisCuts
{
public:
  AliESDMuonTrackCuts(const Char_t* name = "AliESDMuonTrackCuts", const Char_t* title = "");
  virtual ~AliESDMuonTrackCuts();
  Bool_t IsSelected(TObject* obj)
       {return AcceptTrack((AliESDMuonTrack*)obj);}
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
  Bool_t AcceptTrack(AliESDMuonTrack* esdMuTrack);
  TObjArray* GetAcceptedTracks(AliESD* esd);
  Int_t CountAcceptedTracks(AliESD* esd);
  TObjArray* GetAcceptedTracks(AliESDEvent* esd);
  Int_t CountAcceptedTracks(AliESDEvent* esd);

  virtual Long64_t Merge(TCollection* list);
  virtual void Copy(TObject &c) const;
  AliESDMuonTrackCuts(const AliESDMuonTrackCuts& pd);  // Copy Constructor
  AliESDMuonTrackCuts &operator=(const AliESDMuonTrackCuts &c);

  //######################################################

  // track kinematic cut setters
  void SetPRange(Float_t r1=0, Float_t r2=1e10)       {fPMin=r1;   fPMax=r2;}
  void SetPtRange(Float_t r1=0, Float_t r2=1e10)      {fPtMin=r1;  fPtMax=r2;}
  void SetPxRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPxMin=r1;  fPxMax=r2;}
  void SetPyRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPyMin=r1;  fPyMax=r2;}
  void SetPzRange(Float_t r1=-1e10, Float_t r2=1e10)  {fPzMin=r1;  fPzMax=r2;}
  void SetEtaRange(Float_t r1=-1e10, Float_t r2=1e10) {fEtaMin=r1; fEtaMax=r2;}
  void SetRapRange(Float_t r1=-1e10, Float_t r2=1e10) {fRapMin=r1; fRapMax=r2;}

  //######################################################
  
  void SetHistogramsOn(Bool_t b=kFALSE) {fHistogramsOn = b;}
  void DefineHistograms(Int_t color=1);
  virtual Bool_t LoadHistograms(const Char_t* dir = 0);
  void SaveHistograms(const Char_t* dir = 0);
  void DrawHistograms();

  static void EnableNeededBranches(TTree* tree);

protected:
  void Init(); // sets everything to 0

  enum { kNCuts = 7 };

  //######################################################
  static const Char_t* fgkCutNames[kNCuts]; //! names of cuts (for internal use)

  // esd kinematics cuts
  Float_t fPMin,   fPMax;             // definition of the range of the P
  Float_t fPtMin,  fPtMax;            // definition of the range of the Pt
  Float_t fPxMin,  fPxMax;            // definition of the range of the Px
  Float_t fPyMin,  fPyMax;            // definition of the range of the Py
  Float_t fPzMin,  fPzMax;            // definition of the range of the Pz
  Float_t fEtaMin, fEtaMax;           // definition of the range of the eta
  Float_t fRapMin, fRapMax;           // definition of the range of the y

  //######################################################
  // diagnostics histograms
  Bool_t fHistogramsOn;               // histograms on/off
  
  TH1F* fhPt[2];                      //-> pt of esd tracks
  TH1F* fhEta[2];                     //-> eta of esd tracks

  TH1F*  fhCutStatistics;             //-> statistics of what cuts the tracks did not survive
  TH2F*  fhCutCorrelation;            //-> 2d statistics plot

  ClassDef(AliESDMuonTrackCuts, 1)
};


#endif
