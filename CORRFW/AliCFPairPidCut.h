/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


////////////////////////////////////////////////
// Class to define PID cuts on a pair of tracks
// The track pair object to use is AliCFPair
//
// author : renaud.vernet@cern.ch
////////////////////////////////////////////////

#ifndef ALICFPAIRPIDCUT_H
#define ALICFPAIRPIDCUT_H

#include "AliCFCutBase.h"
#include "AliCFTrackCutPid.h"

//__________________________________________________________________________________
// CUT ON TRACK PID FOR V0 DAUGHTERS
//__________________________________________________________________________________

class AliCFPairPidCut : public AliCFCutBase
{
  public :
  AliCFPairPidCut() ;
  AliCFPairPidCut(const Char_t* name, const Char_t* title) ;
  AliCFPairPidCut(const AliCFPairPidCut& c) ;
  AliCFPairPidCut& operator=(const AliCFPairPidCut& c) ;
  virtual ~AliCFPairPidCut() {delete fCutNeg; delete fCutPos; };

  virtual AliCFTrackCutPid* GetNegCut() const {return fCutNeg;}
  virtual AliCFTrackCutPid* GetPosCut() const {return fCutPos;}

  Bool_t IsSelected(TObject *obj); //boolean for detectors
  Bool_t IsSelected(TList* /*list*/) {return kTRUE;}
 private:
  AliCFTrackCutPid* fCutNeg; // PID cut on negative daughter
  AliCFTrackCutPid* fCutPos; // PID cut on positive daughter

  ClassDef(AliCFPairPidCut,1);
};
#endif
    
