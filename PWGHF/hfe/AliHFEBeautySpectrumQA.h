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
#ifndef ALIHFEBEAUTYSPECTRUMQA_H
#define ALIHFEBEAUTYSPECTRUMQA_H

#ifndef ROOT_TNamed
#include <TNamed.h>
#endif


class TObjArray;
class TGraphErrors;
class TObject;


class AliHFEBeautySpectrumQA : public TNamed{
 public:

  enum Results_t{
    kDataProjection = 0,
    kCMProjection = 1,
    kBeforeSC = 2,
    kAfterSC = 3,
    kBeforePE = 4,
    kAfterPE = 5,
    kPEfficiency = 6,
    kBeforeMCE = 7,
    kAfterMCE = 8,
    kMCEfficiency = 9,
    kBeforeU = 10,
    kAfterU = 11,
    kResidualU = 12,
    kUEfficiency = 13,
    kFinalResultUnfolded = 14,
    kFinalResultDirectEfficiency = 15,
    kFinalResultUnfSparse = 16,
    kFinalResultDirectEffSparse = 17,
    kMeasBG = 18,
    kNResults = 19
  };


  enum EfficiencyCorrection_t{
    kMC = 0,
    kParametrized = 1,
    kNTypeEfficiency = 2
  };
   
  void AddResultAt(TObject *obj,Int_t index);
  TObject *GetResult(Int_t index);
  
  void DrawProjections() const;
  void DrawSubtractContamination() const;
  void DrawCorrectWithEfficiency(Int_t typeeff) const;
  void DrawUnfolding() const;
  void DrawResult();
  
  void SetStyle() const;
  void SetWriteToFile(Bool_t writetofile) {fWriteToFile=writetofile; };
  void SetPtMax(Double_t ptmax) {fPtMax = ptmax; };


  TH1D *DivideSpectra(TGraphErrors *ga, TGraphErrors *gb);
  

  AliHFEBeautySpectrumQA();
  AliHFEBeautySpectrumQA(const char* name);
  ~AliHFEBeautySpectrumQA();
  
  protected:
  
  private:

  static const Char_t* fgkNameCanvas[kNTypeEfficiency];     // Name of canvas

    AliHFEBeautySpectrumQA(const AliHFEBeautySpectrumQA &);
    AliHFEBeautySpectrumQA &operator=(const AliHFEBeautySpectrumQA &);
 
    Double_t   fPtMax;             // Pt max to plot
    TObjArray *fListOfResult;     // ListOfResults
    Bool_t fWriteToFile;           // Write plots to eps files
   
    ClassDef(AliHFEBeautySpectrumQA, 1) 
};
#endif

