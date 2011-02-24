#ifndef ALIHFESIGNALCUTS_H
#define ALIHFESIGNALCUTS_H

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
// Signal cuts
// Checks whether a particle (reconstructed or MC) is coming from MC Signal
// For more information see implementation file
//
#ifndef ALIANALYSISCUTS_H
#include "AliAnalysisCuts.h"
#endif

class TList;
class AliMCEvent;
class AliVParticle;
class AliHFEmcQA;

class AliHFEsignalCuts : public AliAnalysisCuts{
  public: 
    AliHFEsignalCuts();
    AliHFEsignalCuts(const Char_t *name, const Char_t *title);
    AliHFEsignalCuts(const AliHFEsignalCuts &ref);
    AliHFEsignalCuts &operator=(const AliHFEsignalCuts &ref);
    virtual ~AliHFEsignalCuts();

    virtual Bool_t IsSelected(TObject *o);
    virtual Bool_t IsSelected(TList * /*l*/) { return kTRUE; };

    Bool_t IsCharmElectron(const TObject * const o) const;
    Bool_t IsBeautyElectron(const TObject * const o) const;
    Bool_t IsGammaElectron(const TObject * const o) const;

    //void SetMCEvent(AliMCEvent *mc) { fMC = mc; }
    void SetMCEvent(AliMCEvent *mc);
  
  protected:
    Int_t GetMotherPDG(const AliVParticle * const track) const;
    Int_t GetTrackPDG(const AliVParticle * const track) const;
    Int_t GetElecSource(const AliVParticle * const track) const ;

  private:
    AliMCEvent *fMC;   //! MC event
    AliHFEmcQA *fMCQA; //! MC QA

    ClassDef(AliHFEsignalCuts, 2)
};
#endif

