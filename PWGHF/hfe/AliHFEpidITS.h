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
class AliHFEpidQAmanager;

class AliHFEpidITS : public AliHFEpidBase{
  public:
    AliHFEpidITS(const Char_t *name);
    AliHFEpidITS(const AliHFEpidITS &ref);
    AliHFEpidITS& operator=(const AliHFEpidITS &ref);
    virtual ~AliHFEpidITS();

    virtual Bool_t InitializePID(Int_t /*run*/);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;

    Double_t GetITSSignalV1(AliVParticle *track);
    Double_t GetITSSignalV2(AliVParticle *track);
  protected:
    void Copy(TObject &o) const;
  private:
    enum{
      kITSsigV1 = 0,
      kITSsigV2 = 1
    };

    ClassDef(AliHFEpidITS, 0)  // PID class for ITS
};
#endif

