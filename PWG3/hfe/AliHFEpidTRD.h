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
#ifndef __ALIHFEPIDTRD_H__
#define __ALIHFEPIDTRD_H__

 #ifndef __ALIHFEPIDBASE_H__
 #include "AliHFEpidBase.h"
 #endif

class TAxis;
class TH1F;
class TMap;
class TString;
class AliVParticle;

class AliHFEpidTRD : public AliHFEpidBase{
  public:
    typedef enum{
      kLQ = 0,
      kNN = 1
    } PIDMethodTRD_t;
    AliHFEpidTRD(const Char_t *name);
    AliHFEpidTRD(const AliHFEpidTRD &ref);
    AliHFEpidTRD& operator=(const AliHFEpidTRD &ref);
    virtual ~AliHFEpidTRD();
    
    virtual Bool_t InitializePID();
    virtual Int_t IsSelected(AliVParticle *track);
    virtual Bool_t HasQAhistos() const { return kFALSE; };

    void LoadTRDthresholds();

    void SetPIDMethod(PIDMethodTRD_t method) { fPIDMethod = method; };
    void SetThresholdFile(Char_t *thresholdFile) { fThresholdFile = thresholdFile; };

  protected:
    void Copy(TObject &ref) const;
    TH1F *GetTRDthresholds(Float_t electronEff);
    void SetTRDthresholds(TH1F *thresholds, Float_t electronEff);

  private:
    TString fThresholdFile;                                 // Threshold file name
    PIDMethodTRD_t fPIDMethod;                              // PID Method: 2D Likelihood or Neural Network
    TMap *fTRDthresholds;                                   //! TRD Thresholds
    TAxis *fTRDelectronEfficiencies;                        //! Electron Efficiencies corresponding to reference histos

  ClassDef(AliHFEpidTRD, 1)     // TRD electron ID class
};

#endif
