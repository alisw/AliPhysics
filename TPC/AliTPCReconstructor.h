#ifndef ALITPCRECONSTRUCTOR_H
#define ALITPCRECONSTRUCTOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliReconstructor.h"
#include "AliTPCRecoParam.h"

class AliTPCParam;
class AliTPCclustererMI;

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

  void SetRecoParam(AliTPCRecoParam * param){ fgkRecoParam = param;}
  static const AliTPCRecoParam* GetRecoParam(){ return fgkRecoParam;}
  //
  static Double_t GetCtgRange()     { return fgkRecoParam->GetCtgRange();}
  static Double_t GetMaxSnpTracker(){ return fgkRecoParam->GetMaxSnpTracker();}
  static Double_t GetMaxSnpTrack()  { return fgkRecoParam->GetMaxSnpTrack();}

  static Int_t StreamLevel()               { return fgStreamLevel;}
  static void  SetStreamLevel(Int_t level) { fgStreamLevel = level;}

private:
  AliTPCParam*         GetTPCParam() const;
  static AliTPCRecoParam *   fgkRecoParam; // reconstruction parameters
  static Int_t               fgStreamLevel; // flag for streaming      - for TPC reconstruction
  AliTPCclustererMI*         fClusterer;   // TPC clusterer

  ClassDef(AliTPCReconstructor, 0)   // class for the TPC reconstruction
};

#endif
