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
//
// Detector PID Class for Monte Carlo PID
// Returns PDG code of a particle
// For more information see the impelmentation file
//
#ifndef ALIHFEPIDMC_H
#define ALIHFEPIDMC_H

 #ifndef ALIHFEPIDBASE_H
 #include "AliHFEpidBase.h"
 #endif

class AliHFEpidQAmanager;

class AliHFEpidMC : public AliHFEpidBase{
  public:
    AliHFEpidMC();
    AliHFEpidMC(const Char_t *name);
    virtual ~AliHFEpidMC(){};
    
    virtual Bool_t InitializePID(Int_t /*run*/);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;

  private:

  ClassDef(AliHFEpidMC, 1)    // MC electron ID class
};

#endif
