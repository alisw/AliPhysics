#ifndef ALIHFEV0INF0_H
#define ALIHFEV0INFO_H

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
// Stores aditional information about the V0 candidates
// author: M.Fasel@gsi.de
//
#ifndef ROOT_TObject
#include <TObject.h>
#endif

class AliVParticle;

class AliHFEV0info : public TObject{
  public:
    AliHFEV0info();
    AliHFEV0info(AliVParticle *track, Int_t idPartnerTrack, Int_t v0id);
    ~AliHFEV0info() {};

    AliVParticle *GetTrack() const { return fTrack; }
    Int_t GetPartnerID() const { return fIDpartnerTrack; }
    Int_t GetV0ID() const { return fIDV0; }
  private:
    AliHFEV0info(const AliHFEV0info &);
    AliHFEV0info &operator=(const AliHFEV0info &);
    AliVParticle *fTrack;       // The track itself
    Int_t fIDpartnerTrack;      // The ID of the parter
    Int_t fIDV0;                // the V0 ID

    ClassDef(AliHFEV0info, 1)
};
#endif
