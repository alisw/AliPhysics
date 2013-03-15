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
// Class AliHFEtpcPIDqa
// Monitoring TPC PID in the HFE PID montioring framework
// More information can be found inside the implementation file
//
#ifndef ALIHFETPCPIDQA_H
#define ALIHFETPCPIDQA_H

#ifndef ALIHFEDETPIDQA_H
#include "AliHFEdetPIDqa.h"
#endif

#ifndef ALIHFEPIDBASE_H
#include "AliHFEpidBase.h"
#endif

class TBrowser;
class TH2;
class AliHFEcollection;
class AliVParticle;

class AliHFEtpcPIDqa : public AliHFEdetPIDqa{
  public:
    AliHFEtpcPIDqa();
    AliHFEtpcPIDqa(const char*name);
    AliHFEtpcPIDqa(const AliHFEtpcPIDqa &o);
    AliHFEtpcPIDqa &operator=(const AliHFEtpcPIDqa &o);
    ~AliHFEtpcPIDqa();
    void Copy(TObject &o) const;
    virtual Long64_t Merge(TCollection *col);
    virtual void Browse(TBrowser *b);
    virtual Bool_t IsFolder() const { return kTRUE; };
  
    virtual void Initialize();
    virtual void ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step);

    void SetBrowseCentrality(Int_t browseCentrality) { browseCentrality < 11  && browseCentrality >= -1 ? fBrowseCentrality = browseCentrality : -1;} // *MENU*

    AliHFEcollection *GetHistograms() const { return fHistos; }
    TH2 *MakeSpectrumdEdx(AliHFEdetPIDqa::EStep_t step, Int_t species = -1, Int_t centralityClass = -1);
    TH2 *MakeSpectrumNSigma(AliHFEdetPIDqa::EStep_t step, Int_t species = -1, Int_t centralityClass = -1);

  protected:
    Double_t GetTPCsignal(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype);
    Double_t GetEta(const AliVParticle *track, AliHFEpidObject::AnalysisType_t anatype);

  private:
    AliHFEcollection *fHistos;        // Container for Histograms
    Int_t fBrowseCentrality;          // Centrality Class for Browser

    ClassDef(AliHFEtpcPIDqa, 1);
};
#endif
