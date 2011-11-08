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
// Class AliHFEbayesPIDqa
// Monitoring Bayes PID in the HFE PID montioring framework
// More information can be found inside the implementation file
//
#ifndef ALIHFEBAYESPIDQA_H
#define ALIHFEBAYESPIDQA_H

#ifndef ALIHFEDETPIDQA_H
#include "AliHFEdetPIDqa.h"
#endif

//#ifndef ALIHFEPIDBASE_H
//#include "AliHFEpidBase.h"
//#endif

class TBrowser;
class TH2;
class AliHFEcollection;
class AliVParticle;
class AliHFEpidBayes;
class AliHFEpidBase;

class AliHFEbayesPIDqa : public AliHFEdetPIDqa{
  public:
    AliHFEbayesPIDqa();
    AliHFEbayesPIDqa(const char*name);
    AliHFEbayesPIDqa(const AliHFEbayesPIDqa &c);
    AliHFEbayesPIDqa &operator=(const AliHFEbayesPIDqa &o);
    ~AliHFEbayesPIDqa();
  
    virtual void Initialize();
    virtual void ProcessTrack(const AliHFEpidObject *track, AliHFEdetPIDqa::EStep_t step);


    AliHFEcollection *GetHistograms() const { return fHistos; }

  protected:
      void CreateProbabilityHistograms();
      void CreateDetectorSignalHistograms();
      Double_t CalcTOFMass(const AliHFEpidObject *track);
private:
    static const Int_t       fgkNBinsProb[5];         // number of bins
    static const Double_t    fgkMinBinsProb[5];       // bin minimum
    static const Double_t    fgkMaxBinsProb[5];       // bin maximum
//    AliHFEpidBayes *fBAYESpid;        // HFE PID for TRD
    AliHFEcollection *fHistos;        // Container for Histograms

    ClassDef(AliHFEbayesPIDqa, 1);
};
#endif
