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
#ifndef ALIHFEINCLUSIVESPECTRUMQA_H
#define ALIHFEINCLUSIVESPECTRUMQA_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif


class TObjArray;
class TGraphErrors;
class TObject;


class AliHFEInclusiveSpectrumQA : public TNamed{
 public:

  enum Results_t{
    kDataProjection = 0,
    kCMProjection = 1,
    kBeforeSC = 2,
    kAfterSC = 3,
    kBeforeV0 = 4,
    kAfterV0 = 5,
    kV0Efficiency = 6,
    kBeforePE = 7,
    kAfterPE = 8,
    kPEfficiency = 9,
    kBeforeMCE = 10,
    kAfterMCE = 11,
    kMCEfficiency = 12,
    kBeforeU = 13,
    kAfterU = 14,
    kResidualU = 15,
    kUEfficiency = 16,
    kFinalResultUnfolded = 17,
    kFinalResultDirectEfficiency = 18,
    kBeforeSPB = 19,
    kAfterSPB = 20,
    kNResults = 21
  };


  enum EfficiencyCorrection_t{
    kV0 = 0,
    kMC = 1,
    kParametrized = 2,
    kNTypeEfficiency = 3
  };
   
  void AddResultAt(TObject *obj,Int_t index);
  TObject *GetResult(Int_t index);
  
  void DrawProjections() const;
  void DrawSubtractContamination() const;
  void DrawSubtractPhotonicBackground() const;
  void DrawCorrectWithEfficiency(Int_t typeeff) const;
  void DrawUnfolding() const;
  void DrawResult();
  
  void SetStyle() const;
  void SetWriteToFile(Bool_t writetofile) {fWriteToFile=writetofile; };
  void SetPtMax(Double_t ptmax) {fPtMax = ptmax; };


  TH1D *DivideSpectra(TGraphErrors *ga, TGraphErrors *gb);
  

  AliHFEInclusiveSpectrumQA();
  AliHFEInclusiveSpectrumQA(const char* name);
  ~AliHFEInclusiveSpectrumQA();
  
  protected:
  
  private:

  static const Char_t* fgkNameCanvas[kNTypeEfficiency];     // Name of canvas

    AliHFEInclusiveSpectrumQA(const AliHFEInclusiveSpectrumQA &);
    AliHFEInclusiveSpectrumQA &operator=(const AliHFEInclusiveSpectrumQA &);
 
    Double_t   fPtMax;             // Pt max to plot
    TObjArray *fListOfResult;     // ListOfResults
    Bool_t fWriteToFile;           // Write plots to eps files
   
    ClassDef(AliHFEInclusiveSpectrumQA, 1) 
};
#endif

