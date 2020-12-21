#ifndef ALIDIELECTRONCLUSTERCUTS_H
#define ALIDIELECTRONCLUSTERCUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           #
//#         Class AliDielectronClusterCuts                    #
//#                                                           #
//#  Authors:                                                 #
//#   Lucas Altenkamper, UiB / lucas.altenkamper@cern.ch      #
//#                                                           #
//#############################################################

#include <AliPID.h>
#include <AliAnalysisCuts.h>

class AliDielectronClusterCuts : public AliAnalysisCuts {
  
public:
  
  enum Detector {
    kEMCal=0, // calorimeter = EMCal/DCal
    kPHOS,    // calorimeter = PHOS
    kAny
  };
  AliDielectronClusterCuts();
  AliDielectronClusterCuts(const char*name, const char* title);
  virtual ~AliDielectronClusterCuts();

  // setters
  void SetCaloType(Short_t type) { fCaloType=type; }
  void SetNCellsMinCut(Short_t val) { fMinNCells=val; }
  void SetRejectExotics() { fRejectExotics=kTRUE; }
  void SetRequireTrackMatch(Bool_t req) { fRequireTrackMatch=req; }
  void SetM02Cut(Float_t min, Float_t max) { fM02Min=min; fM02Max=max; }
  void SetM20Cut(Float_t min, Float_t max) { fM20Min=min; fM20Max=max; }
  void SetTrackDxCut(Float_t min, Float_t max) { fTrackDxMin=min; fTrackDxMax=max; }
  void SetTrackDzCut(Float_t min, Float_t max) { fTrackDzMin=min; fTrackDzMax=max; }

  // getters
  Short_t GetCaloType() const { return fCaloType; }
  Short_t GetNCellsMinCut() const { return fMinNCells; }
  Bool_t  GetRejectExotics() const { return fRejectExotics; }
  Bool_t  GetRequireTrackMatch() const { return fRequireTrackMatch; }
  void    GetM02Cut(Float_t &min, Float_t &max) const { min=fM02Min; max=fM02Max; }
  void    GetM20Cut(Float_t &min, Float_t &max) const { min=fM20Min; max=fM20Max; }
  void    GetTrackDxCut(Float_t &min, Float_t &max) const { min=fTrackDxMin; max=fTrackDxMax; }
  void    GetTrackDzCut(Float_t &min, Float_t &max) const { min=fTrackDzMin; max=fTrackDzMax; }

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* cluster);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}

private:
  
  AliDielectronClusterCuts(const AliDielectronClusterCuts &c);
  AliDielectronClusterCuts &operator=(const AliDielectronClusterCuts &c);

  Short_t fCaloType;          // calorimeter type (EMCal/DCal, PHOS, any)
  Short_t fMinNCells;         // min number of cells per cluster
  Bool_t  fRejectExotics;     // reject exotic cluster
  Bool_t  fRequireTrackMatch; // require at least one track matched to cluster
  Float_t fM02Min;            // min cluster M02
  Float_t fM02Max;            // max cluster M02
  Float_t fM20Min;            // min cluster M20
  Float_t fM20Max;            // max cluster M20
  Float_t fTrackDxMin;        // min track Dx
  Float_t fTrackDxMax;        // max track Dx
  Float_t fTrackDzMin;        // min track Dz
  Float_t fTrackDzMax;        // max track Dz
  
  ClassDef(AliDielectronClusterCuts,2);
};

#endif
