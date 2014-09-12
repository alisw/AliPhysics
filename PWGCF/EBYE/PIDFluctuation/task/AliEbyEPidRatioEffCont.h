#ifndef ALIEBYEPIDRSTIOEFFCONT_H
#define ALIEBYEPIDRSTIOEFFCONT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//=========================================================================//

class AliVTrack;

#include "THnSparse.h"

#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioEffCont: public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioEffCont();
  virtual ~AliEbyEPidRatioEffCont();
  virtual void Process();
  THnSparseF* GetHnEff()  {return fHnEff;}
  THnSparseF* GetHnCont() {return fHnCont;}

 private:

  AliEbyEPidRatioEffCont(const AliEbyEPidRatioEffCont&); // not implemented
  AliEbyEPidRatioEffCont& operator=(const AliEbyEPidRatioEffCont&); // not implemented

  
  virtual void Init();
  virtual void CreateHistograms();
  virtual void Reset();
  virtual Int_t Setup();
  void FillMCLabels(); 
  void FillMCEffHist();
  void CheckContTrack(AliVTrack* track, Int_t iPid, Int_t gPdgCode);
    
  Int_t            ***fLabelsRec;             //! 2x nTracks large array with labels for MC particles
  THnSparseF         *fHnEff;                 //  THnSparseF efficiency 
  THnSparseF         *fHnCont;                //  THnSparseF contamination
  
  ClassDef(AliEbyEPidRatioEffCont, 1);
};

#endif
