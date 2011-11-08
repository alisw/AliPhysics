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
// Steering class for electron identification
// For more information please check the implementation file
// ...
//
#ifndef ALIHFEPIDESD_H
#define ALIHFEPIDESD_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class AliVParticle;

class AliHFEpidESD : public AliHFEpidBase{
  public:
    AliHFEpidESD();
    AliHFEpidESD(const AliHFEpidESD &ref);
    AliHFEpidESD &operator=(const AliHFEpidESD &ref);
    virtual ~AliHFEpidESD();  

    virtual void InitializePID();
    virtual Int_t IsSelected(const AliVParticle *track);
    void SetRequireTOFRange(Double_t pmin, Double_t pmax);
    void SetRequireTRDRange(Double_t pmin, Double_t pmax);
    void SetRequireMinTRDtracklets();
  private:
    Double_t fPminTRD;		// Min. Momentum where TRD PID is required
    Double_t fPmaxTRD;		// Max. Momentum where TRD PID is required
    Double_t fPminTOF;		// Min. Momentum where TOF PID is required
    Double_t fPmaxTOF;		// Max. Momentum where TOF PID is required
    Int_t fMinTrackletsTRD;	// Min. Number of TRD tracklets in region where TRD PID is required
  ClassDef(AliHFEpidESD) // ESD PID class
};
#endif
