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
// Class AliHFEemcalPIDqa
// Monitoring EMCAL PID in the HFE PID montioring framework
// More information can be found inside the implementation file
//
#ifndef ALIHFEEMCALPIDQA_H
#define ALIHFEEMCALPIDQA_H

#ifndef ALIHFEDETPIDQA_H
#include "AliHFEdetPIDqa.h"
#endif

class TH1;
class TH2;
class AliHFEcollection;
class AliHFEpidObject;
class AliESDtrack;
class AliAODTrack;

class AliHFEemcalPIDqa : public AliHFEdetPIDqa{
  public:
    AliHFEemcalPIDqa();
    AliHFEemcalPIDqa(const char*name);
    AliHFEemcalPIDqa(const AliHFEemcalPIDqa &o);
    AliHFEemcalPIDqa &operator=(const AliHFEemcalPIDqa &o);
    ~AliHFEemcalPIDqa();
    void Copy(TObject &o) const;
    virtual Long64_t Merge(TCollection *col);
  
    virtual void Initialize();
    virtual void ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step);

    //void HistEnergyMomMatch(double pT, double eop);
    void HistEnergyMomMatch(const AliHFEpidObject *track, double eop);
    TH1 *GetHistogram(const char *name); 
    AliHFEcollection *GetHistoCollection() const { return fHistos; }

  protected:
    void ProcessESDtrack(const AliESDtrack *track, AliHFEdetPIDqa::EStep_t step, Int_t species);
    void ProcessAODtrack(const AliAODTrack *track, AliHFEdetPIDqa::EStep_t step, Int_t species);
  private:
    AliHFEcollection *fHistos;        // Container for Histograms

    ClassDef(AliHFEemcalPIDqa, 1);
};
#endif
