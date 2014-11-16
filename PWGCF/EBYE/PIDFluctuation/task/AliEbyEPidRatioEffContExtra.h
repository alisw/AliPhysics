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
//=========================================================================//


class AliVTrack;

#include "THnSparse.h"

#include "AliEbyEPidRatioBase.h"

class AliEbyEPidRatioEffContExtra: public AliEbyEPidRatioBase {

 public:

  AliEbyEPidRatioEffContExtra();
  virtual ~AliEbyEPidRatioEffContExtra();
  
  virtual void Process();

  THnSparseF* GetHnEff(Int_t type, Int_t i)  {
    if (i == 0) {
      if (type == 0)      return  fHnNchEMc;
      else if (type == 1) return  fHnNchERec;
      else if (type == 2) return  fHnNchCMc;
      else if (type == 4) return  fHnNchCRec;
    } else if (i == 1) {
      if (type == 0)      return  fHnNpiEMc;
      else if (type == 1) return  fHnNpiERec;
      else if (type == 2) return  fHnNpiCMc;
      else if (type == 4) return  fHnNpiCRec;
    } else if (i == 2) {
      if (type == 0)      return  fHnNkaEMc;
      else if (type == 1) return  fHnNkaERec;
      else if (type == 2) return  fHnNkaCMc;
      else if (type == 4) return  fHnNkaCRec;
    } else if (i == 3) {
      if (type == 0)      return  fHnNprEMc;
      else if (type == 1) return  fHnNprERec;
      else if (type == 2) return  fHnNprCMc;
      else if (type == 4) return  fHnNprCRec;
    }
    return 0; 
  }

 private:

  AliEbyEPidRatioEffContExtra(const AliEbyEPidRatioEffContExtra&); // not implemented
  AliEbyEPidRatioEffContExtra& operator=(const AliEbyEPidRatioEffContExtra&); // not implemented

  virtual void Init();

  /** Create the efficiency / contamination THnSparse */
  virtual void CreateHistograms();

  /** Event-wise Reset - Can be implemented by every class */
  virtual void Reset();

  /** Event-wise Setup - Can be implemented by every class */
  virtual Int_t Setup();

  // -----------------------------------------------------------------------

  /** Fill MC labels */
  void FillMCLabels(Int_t ipid); 

  /** Fill efficiency THnSparse */
  void FillMCEffHist(Int_t ipid);

  /** Check if particle is contamination */
  void CheckContTrack(AliVTrack* track, Int_t ipid);
      
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  // =======================================================================
  Int_t             ***fLabelsRec;     //! 2x nTracks large array with labels for MC particles
   // =======================================================================
 
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
