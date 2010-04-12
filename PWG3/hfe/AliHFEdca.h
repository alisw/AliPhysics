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
// Class for checking impact parameter (DCA) study 
// Study DCA in rphi (xy) and z
// resolution and pull
// 

#ifndef ALIHFEDCA_H
#define ALIHFEDCA_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

class TChain;
class TTree;
class TFile;

class TString;
class TList;

class TObjArray;
class AliStack;
class AliMCEvent;

class AliESDEvent;
class AliESDtrack;
class AliESDVertex;

class AliHFEdca : public TObject{

 public:  
  enum{
    kPDGelectron = 11,
    kPDGmuon = 13,
    kPDGpion = 211,
    kPDGkaon = 321,
    kPDGproton = 2212
  };
 
  enum{
    kNParticles = 12,
    kNPtBins = 43,   
    kNDcaVar = 2, 
    kNPullVar = 2
  };

  AliHFEdca();
  AliHFEdca(const AliHFEdca &p); // copy constructor
  AliHFEdca &operator=(const AliHFEdca &); // assignment operator

  virtual ~AliHFEdca();

  void Initialize();
  void CreateHistogramsPull(TList *pullList);  
  void CreateHistogramsResidual(TList *residualList);  
  void InitAnalysis();  
  void FillHistograms(AliESDEvent *esdEvent,  AliESDtrack *track,  AliMCEvent *mcEvent);
  void PostAnalysis() const;


 private:   
  static const Char_t *fgkParticles[kNParticles];  // particle names
  static const Int_t fgkColorPart[kNParticles]; // colors for particles

  static const Float_t fgkPtIntv[kNPtBins+1];  // pt intervals

  static const Char_t* fgkDcaVar[kNDcaVar];  // dca variables
  static const Char_t* fgkDcaVarTitle[kNDcaVar]; // titles for dca variables

  static const Char_t* fgkPullDcaVar[kNPullVar];  // pull variables
  static const Char_t* fgkPullDcaVarTitle[kNPullVar]; // titles for pull variables

  TH1F* fHistDcaXYRes[kNParticles][kNPtBins];  //! residuals in XY
  TH1F* fHistDcaZRes[kNParticles][kNPtBins];   //! residuals in Z
  TH1F* fHistDcaXYPull[kNParticles][kNPtBins]; //! pulls XY
  TH1F* fHistDcaZPull[kNParticles][kNPtBins];  //! pulls Z

  TList *fResidualList;   //! collection of histograms of residual
  TList *fPullList;       //! collection of histograms of pull
  
  ClassDef(AliHFEdca, 1);
};

#endif
