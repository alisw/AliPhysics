#ifndef ALIHFECUTSTEP_H
#define ALIHFECUTSTEP_H

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

/* $Id$ */ 

//
// Cut step class
// Select all tracks surviving cuts in one special cut step
// Used in AliHFEtrackFilter
// 
#include <TNamed.h>

class TObjArray;
class AliAnalysisCuts;
class AliMCEvent;
class AliVEvent;

class AliHFEcutStep : public TNamed{
    public:
      AliHFEcutStep(const Char_t *name);
      AliHFEcutStep(const AliHFEcutStep &o);
      AliHFEcutStep &operator=(const AliHFEcutStep &o);
      virtual void Copy(TObject &o) const;
      ~AliHFEcutStep();
      
      void AddCut(AliAnalysisCuts *cut);
      AliAnalysisCuts *GetCut(const Char_t *name);
      Bool_t IsSelected(TObject *o);

      void SetMC(AliMCEvent *mc);
      void SetRecEvent(AliVEvent *mc);

    private:
      TObjArray *fCuts; // List of cuts in one cut step
      
      ClassDef(AliHFEcutStep, 1)
};
#endif

