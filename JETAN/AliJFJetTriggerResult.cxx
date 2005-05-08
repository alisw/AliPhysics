// $Id$

#include <Riostream.h>
#include "AliJFJetTriggerResult.h"

ClassImp(AliJFJetTriggerResult)

AliJFJetTriggerResult::AliJFJetTriggerResult() : TObject(),
                                           fPtMin(0),fPtMax(0),
                                           fEtaMin(0),fEtaMax(0),
                                           fPhiMin(0),fPhiMax(0),
				           fNeutral(kTRUE),fCharged(kTRUE),fEM(kTRUE)
{
}

AliJFJetTriggerResult::AliJFJetTriggerResult(Float_t PtMin,Float_t PtMax,Float_t EtaMin,Float_t EtaMax,Float_t PhiMin,Float_t PhiMax,Bool_t Neutral,Bool_t EM,Bool_t Charged) : TObject(),
                              fPtMin(PtMin),fPtMax(PtMax),
			      fEtaMin(EtaMin),fEtaMax(EtaMax),
                              fPhiMin(PhiMin),fPhiMax(PhiMax),
			      fNeutral(Neutral),fCharged(Charged),fEM(EM)

{
}


AliJFJetTriggerResult::~AliJFJetTriggerResult()
{
}
