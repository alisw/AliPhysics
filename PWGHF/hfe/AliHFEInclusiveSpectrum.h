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
// Class for spectrum correction
// Subtraction of hadronic background, Unfolding of the data and
// Renormalization done here
// For more information see the implementation file
//
#ifndef ALIHFEINCLUSIVESPECTRUM_H
#define ALIHFEINCLUSIVESPECTRUM_H

#include "AliHFECorrectSpectrumBase.h"


class TGraphErrors;
class TObject;
class TH1;
class TF1;
class TList;
class TObjArray;
class AliCFContainer;
class AliHFEcontainer;
class AliCFDataGrid;
class AliCFEffGrid;
class AliHFEInclusiveSpectrumQA;

class AliHFEInclusiveSpectrum : public AliHFECorrectSpectrumBase{
  public:
  
  AliHFEInclusiveSpectrum(const char* name);
    ~AliHFEInclusiveSpectrum();
    

    virtual Bool_t Init(const AliHFEcontainer *datahfecontainer, const AliHFEcontainer *mchfecontainer, const AliHFEcontainer */*bghfecontainer*/=0x0, const AliHFEcontainer *v0hfecontainer=0x0,AliCFContainer *photoniccontainerD=0x0);
    virtual Bool_t Correct(Bool_t subtractcontamination=kTRUE,  Bool_t subtractphotonic=kFALSE);
   
    AliCFDataGrid *SubtractBackground();
    AliCFDataGrid *SubtractPhotonicBackground();
    AliCFDataGrid *CorrectV0Efficiency(AliCFDataGrid* const bgsubpectrum = 0x0);
    AliCFDataGrid *CorrectParametrizedEfficiency(AliCFDataGrid* const bgsubpectrum = 0x0);
    THnSparse *Unfold(AliCFDataGrid* const bgsubpectrum = 0x0);
    AliCFDataGrid *CorrectForEfficiency(AliCFDataGrid* const bgsubpectrum = 0x0);

    void WriteResults(const char *filename);
   
 private:
    AliHFEInclusiveSpectrum(const AliHFEInclusiveSpectrum &ref);
    AliHFEInclusiveSpectrum &operator=(const AliHFEInclusiveSpectrum &ref);
    virtual void Copy(TObject &o) const;
 
    AliHFEInclusiveSpectrumQA *fQA; // QA
   
    ClassDef(AliHFEInclusiveSpectrum, 1) 
};
#endif

