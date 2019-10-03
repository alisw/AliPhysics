#ifndef ALIEBYEPIDRSTIOEFFCONTEXTRA_H
#define ALIEBYEPIDRSTIOEFFCONTEXTRA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */
 
//=========================================================================//
//             AliEbyE Analysis for Particle Ratio Fluctuation             //
//                   Deepika Rathee  | Satyajit Jena                       //
//                   drathee@cern.ch | sjena@cern.ch                       //
//                  Date: Wed Jul  9 18:38:30 CEST 2014                    //
//          New approch to find particle ratio to reduce memory            //
//                             (Test Only)                                 //
//        Copied from NetParticle Classes
//        Origin: Authors: Jochen Thaeder <jochen@thaeder.de>
//                         Michael Weber <m.weber@cern.ch>
//=========================================================================//


class AliVTrack;

#include "THnSparse.h"

#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioEffContExtra: public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioEffContExtra();
  virtual ~AliEbyEPidRatioEffContExtra();
  
  virtual void Process();

  THnSparseF* GetHnEff(Int_t i, Int_t type)  {
    if (i == 0) {
      if (type == 0)      return  fHnNchEMc;
      else if (type == 1) return  fHnNchERec;
      else if (type == 2) return  fHnNchCMc;
      else if (type == 3) return  fHnNchCRec;
    } else if (i == 1) {
      if (type == 0)      return  fHnNpiEMc;
      else if (type == 1) return  fHnNpiERec;
      else if (type == 2) return  fHnNpiCMc;
      else if (type == 3) return  fHnNpiCRec;
    } else if (i == 2) {
      if (type == 0)      return  fHnNkaEMc;
      else if (type == 1) return  fHnNkaERec;
      else if (type == 2) return  fHnNkaCMc;
      else if (type == 3) return  fHnNkaCRec;
    } else if (i == 3) {
      if (type == 0)      return  fHnNprEMc;
      else if (type == 1) return  fHnNprERec;
      else if (type == 2) return  fHnNprCMc;
      else if (type == 3) return  fHnNprCRec;
    }
    return 0; 
  }

 private:

  AliEbyEPidRatioEffContExtra(const AliEbyEPidRatioEffContExtra&); // not implemented
  AliEbyEPidRatioEffContExtra& operator=(const AliEbyEPidRatioEffContExtra&); // not implemented

  virtual void Init();

  virtual void CreateHistograms();
  virtual void Reset();
  virtual Int_t Setup();
  void FillMCLabels(Int_t ipid); 
  void FillMCEffHist(Int_t ipid);
  void CheckContTrack(AliVTrack* track, Int_t ipid);
      
 
  Int_t             ***fLabelsRec;     //! 2x nTracks large array with labels for MC particles
 
 
  THnSparseF         *fHnNchEMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNchERec;    //  THnSparseF efficiency 

  THnSparseF         *fHnNpiEMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNpiERec;    //  THnSparseF efficiency 
  
  THnSparseF         *fHnNkaEMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNkaERec;    //  THnSparseF efficiency 
  
  THnSparseF         *fHnNprEMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNprERec;    //  THnSparseF efficiency 
 
  THnSparseF         *fHnNchCMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNchCRec;    //  THnSparseF efficiency 

  THnSparseF         *fHnNpiCMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNpiCRec;    //  THnSparseF efficiency 
  
  THnSparseF         *fHnNkaCMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNkaCRec;    //  THnSparseF efficiency 
  
  THnSparseF         *fHnNprCMc;     //  THnSparseF efficiency 
  THnSparseF         *fHnNprCRec;    //  THnSparseF efficiency 

  // -----------------------------------------------------------------------

  ClassDef(AliEbyEPidRatioEffContExtra, 1);
};

#endif
