#ifndef ALITPCRECONSTRUCTOR_H
#define ALITPCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"
#include "AliTPCRecoParam.h"

class AliTPCParam;


class AliTPCReconstructor: public AliReconstructor {
public:
  AliTPCReconstructor();
  virtual ~AliTPCReconstructor() {if (fgkRecoParam) delete fgkRecoParam;};

  virtual void         Reconstruct(AliRunLoader* runLoader) const;
  virtual void         Reconstruct(AliRunLoader* runLoader,
				   AliRawReader* rawReader) const;
  virtual void         Reconstruct(TTree* digitsTree, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(digitsTree,clustersTree);
  }
  virtual void         Reconstruct(AliRawReader* rawReader, TTree* clustersTree) const {
    AliReconstructor::Reconstruct(rawReader,clustersTree);
  }
  virtual AliTracker*  CreateTracker(AliRunLoader* runLoader) const;
  virtual void         FillESD(AliRunLoader* runLoader, AliESD* esd) const;
  virtual void         FillESD(TTree* digitsTree, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(digitsTree,clustersTree,esd);
  }
  virtual void         FillESD(AliRawReader* rawReader, TTree* clustersTree, 
			       AliESD* esd) const {
    AliReconstructor::FillESD(rawReader,clustersTree,esd);
  }
  virtual void         FillESD(AliRunLoader* runLoader, 
			       AliRawReader* rawReader, AliESD* esd) const {
    AliReconstructor::FillESD(runLoader,rawReader,esd);
  }

  void SetRecoParam(AliTPCRecoParam * param){ fgkRecoParam = param;}
  static const AliTPCRecoParam* GetRecoParam(){ return fgkRecoParam;}
  //
  static Double_t GetCtgRange()     { return fgkRecoParam->GetCtgRange();}
  static Double_t GetMaxSnpTracker(){ return fgkRecoParam->GetMaxSnpTracker();}
  static Double_t GetMaxSnpTrack()  { return fgkRecoParam->GetMaxSnpTrack();}

  static Int_t StreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}

private:
  AliTPCParam*         GetTPCParam(AliRunLoader* runLoader) const;
  static AliTPCRecoParam *   fgkRecoParam; // reconstruction parameters
  static Int_t               fgStreamLevel; // flag for streaming      - for TPC reconstruction

  ClassDef(AliTPCReconstructor, 0)   // class for the TPC reconstruction
};

#endif
