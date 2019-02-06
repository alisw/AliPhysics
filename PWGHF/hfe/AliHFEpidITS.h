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
// PID development class for ITS
// does proton rejection via dE/dx
// For more information see implementation file
#ifndef ALIHFEPIDITS_H
#define ALIHFEPIDITS_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class AliVParticle;
class AliVTrack;
class AliPID;

class AliHFEpidQAmanager;

class AliHFEpidITS : public AliHFEpidBase{
  public:
    AliHFEpidITS();
    AliHFEpidITS(const Char_t *name);
    AliHFEpidITS(const AliHFEpidITS &ref);
    AliHFEpidITS& operator=(const AliHFEpidITS &ref);
    virtual ~AliHFEpidITS();

    void SetITSnSigma(Float_t nSigmalow, Float_t nSigmahigh);
    void SetMeanShift(Double_t meanshift) { fMeanShift = meanshift; }
    virtual Bool_t InitializePID(Int_t /*run*/);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;

    Double_t GetITSNsigmaCorrected(const AliVTrack *track) const;
  protected:
    void Copy(TObject &o) const;
  private:
    enum{
      kITSsigV1 = 0,
      kITSsigV2 = 1
    };
    Float_t    fNsigmaITSlow;          // ITS sigma band
    Float_t    fNsigmaITShigh;          // ITS sigma band
    Double_t   fMeanShift;          // Correction for possible shift of the electron band
    ClassDef(AliHFEpidITS, 2)  // PID class for ITS
};
#endif

