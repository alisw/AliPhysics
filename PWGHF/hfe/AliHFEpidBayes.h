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
// Class for Bayes PID
// electron selection via combined PID
// For more information please check the implementation file
//
#ifndef ALIHFEPIDBAYES_H
#define ALIHFEPIDBAYES_H

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

//#ifndef ALIPID_H
//#include "AliPID.h"
//#endif

class TList;
class AliAODTrack;
class AliAODMCParticle;
class AliESDtrack;
class AliMCParticle;
class AliVParticle;
class AliHFEcollection;
class AliHFEpidQAmanager;
class AliPIDCombined;
class AliPID;

class AliHFEpidBayes : public AliHFEpidBase{
  public:
    AliHFEpidBayes();
    AliHFEpidBayes(const Char_t *name);
    AliHFEpidBayes(const AliHFEpidBayes &ref);
    AliHFEpidBayes &operator=(const AliHFEpidBayes &ref);
    void Copy(TObject &ref) const;
    virtual ~AliHFEpidBayes();
    
    virtual Bool_t InitializePID(Int_t /*run*/);
    virtual Int_t IsSelected(const AliHFEpidObject *track, AliHFEpidQAmanager *pidqa) const;
    void CalcCombProb(const AliHFEpidObject *track,const AliPIDResponse *fkPIDResponse, Double_t* probTPCTOF) const;
    void SetBayesDetectorMask(Int_t detmask) { fDetMask = detmask; };
    void SetBayesPIDThreshold(Float_t pidthres) {fpidthres=pidthres;};

protected:
  
  private:
  AliPIDCombined       *fPIDCombined;            //! combined PID object
  Int_t                fDetMask;                 // detector mask
  Int_t                fDetMaskDefault;          // default detector mask
  Int_t                fDetMaskDefaultandTRD;    // default detector mask +TRD
  Float_t              fpidthres;                // pid threshold
  ClassDef(AliHFEpidBayes, 1)                    // Bayes Electron ID class
};

 
#endif
