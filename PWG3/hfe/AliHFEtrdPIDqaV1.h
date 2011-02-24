#ifndef ALIHFETRDPIDQAV1_H
#define ALIHFETRDPIDQAV1_H

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
// Class AliHFEtrdPIDqaV1
// Monitoring TRD PID in the HFE PID montioring framework. 
// More information can be found inside the implementation file
//
#ifndef ALIHFEDETPIDQAV1_H
#include "AliHFEdetPIDqa.h"
#endif

class AliESDtrack;
class AliAODTrack;
class AliHFEcollection;
class AliHFEpidObject;
class TBrowser;
class TCollection;

class AliHFEtrdPIDqaV1 : public AliHFEdetPIDqa{
  public:
    AliHFEtrdPIDqaV1();
    AliHFEtrdPIDqaV1(const Char_t *name);
    AliHFEtrdPIDqaV1(const AliHFEtrdPIDqaV1 &c);
    AliHFEtrdPIDqaV1 &operator=(const AliHFEtrdPIDqaV1 &o);
    ~AliHFEtrdPIDqaV1(){}
    virtual Long64_t Merge(TCollection *coll);
    virtual void Browse(TBrowser *b);
    virtual Bool_t IsFolder() const { return kTRUE; };

    virtual void Initialize();
    virtual void ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step);

    TH2 *MakeTPCspectrumNsigma(AliHFEdetPIDqa::EStep_t step, Int_t species = -1);
    TH2 *MakeTRDspectrumTM(AliHFEdetPIDqa::EStep_t step, Int_t species = -1);
    TH2 *MakeTRDlikelihoodDistribution(AliHFEdetPIDqa::EStep_t step, Int_t species = -1);
    TH2 *MakeTRDchargeDistribution(AliHFEdetPIDqa::EStep_t step, Int_t species = -1);
  protected:
    void ProcessESDtrack(const AliESDtrack *esdtrack, AliHFEdetPIDqa::EStep_t step, Int_t species);
    void ProcessAODtrack(const AliAODTrack *aodtrack, AliHFEdetPIDqa::EStep_t step, Int_t species);
    AliHFEcollection *fHistos; // QA histograms

    ClassDef(AliHFEtrdPIDqaV1, 1)     // Base class for detector PID QA
};

#endif
