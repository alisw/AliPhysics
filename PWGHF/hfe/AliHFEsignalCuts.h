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
// Signal cuts
// Checks whether a particle (reconstructed or MC) is coming from MC Signal
// For more information see implementation file
//
#ifndef ALIHFESIGNALCUTS_H
#define ALIHFESIGNALCUTS_H

#ifndef ALIANALYSISCUTS_H
#include "AliAnalysisCuts.h"
#endif

class TList;
class AliMCEvent;
class AliVParticle;
class AliHFEmcQA;
class TClonesArray;

class AliHFEsignalCuts : public AliAnalysisCuts{
  public: 
    enum ESignalSource_t{
       kEleCharm = 0,
       kEleBeauty = 1,
       kEleGamma = 2,
       kEleNonHFE = 3,
       kEleJPsi = 4,
       kEleBtoJPsi = 5,
       kEleKe3 =6,
       kOther = 7
    };
    AliHFEsignalCuts();
    AliHFEsignalCuts(const Char_t *name, const Char_t *title);
    AliHFEsignalCuts(const AliHFEsignalCuts &ref);
    AliHFEsignalCuts &operator=(const AliHFEsignalCuts &ref);
    virtual ~AliHFEsignalCuts();

    virtual Bool_t IsSelected(TObject *o);
    virtual Bool_t IsSelected(TList * /*l*/) { return kTRUE; };

    ESignalSource_t GetSignalSource(const TObject *const o) const;

    Bool_t IsCharmElectron(const TObject * const o) const { return GetSignalSource(o) == kEleCharm; }
    Bool_t IsBeautyElectron(const TObject * const o) const { return GetSignalSource(o) == kEleBeauty; }
    Bool_t IsGammaElectron(const TObject * const o) const { return GetSignalSource(o) == kEleGamma; }
    Bool_t IsNonHFElectron(const TObject * const o) const { return GetSignalSource(o) == kEleNonHFE; }
    Bool_t IsJpsiElectron(const TObject * const o) const { return GetSignalSource(o) == kEleJPsi; }
    Bool_t IsB2JpsiElectron(const TObject * const o) const { return GetSignalSource(o) == kEleBtoJPsi; }
    Bool_t IsKe3Electron(const TObject * const o) const { return GetSignalSource(o) == kEleKe3; }

    /*********************************************
     *           Old legacy code                 *
     *********************************************/
    Bool_t IsCharmElectronOld(const TObject * const o) const;
    Bool_t IsBeautyElectronOld(const TObject * const o) const;
    Bool_t IsGammaElectronOld(const TObject * const o) const;
    Bool_t IsNonHFElectronOld(const TObject * const o) const;
    Bool_t IsJpsiElectronOld(const TObject * const o) const;
    Bool_t IsB2JpsiElectronOld(const TObject * const o) const;
    Bool_t IsKe3ElectronOld(const TObject * const o) const;

    //void SetMCEvent(AliMCEvent *mc) { fMC = mc; }
    void SetMCEvent(AliMCEvent *mc);
    void SetMCAODInfo(TClonesArray *mcarray);
    const AliHFEmcQA *GetMCQAObject() const { return fMCQA; }
  
  protected:
    Int_t GetMotherPDG(const AliVParticle * const track) const;
    Int_t GetTrackPDG(const AliVParticle * const track) const;
    Int_t GetElecSource(const AliVParticle * const track) const ;

  private:
    AliMCEvent *fMC;                   //! MC event
    TClonesArray *fAODArrayMCInfo;     //! MC info particle AOD
    AliHFEmcQA *fMCQA;                 //! MC QA

    ClassDef(AliHFEsignalCuts, 3)
};
#endif

