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
// Toolkit containing various usefull things
// Usable everywhere in the hfe software package
// For more information see the cxx file
//
#ifndef ALIHFETOOLS_H
#define ALIHFETOOLS_H

#include <TObject.h>

class TParticle;
class AliAODMCParticle;

class AliHFEtools : public TObject{
  public:
    AliHFEtools();
    ~AliHFEtools() {};

    static Double_t *MakeLinearBinning(Int_t nBins, Double_t ymin, Double_t ymax);
    static Double_t *MakeLogarithmicBinning(Int_t nBins, Double_t ymin, Double_t ymax);
    Bool_t    BinLogAxis(TObject *o, Int_t dim);
    static Float_t GetRapidity(TParticle *part);
    static Float_t GetRapidity(AliAODMCParticle *part); // return rapidity

    ClassDef(AliHFEtools, 0)
};
#endif
