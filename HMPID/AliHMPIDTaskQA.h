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

//==============================================================================
// AliHMPIDTaskQA - Class representing a quality check tool of HMPID
// A set of histograms is created.
//==============================================================================

#ifndef ALIHMPIDTASKQA_H
#define ALIHMPIDTASKQA_H

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class TH1;
class TParticle;
class TFile;
class AliESDtrack;
class AliESDEvent;

class AliHMPIDTaskQA : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDTaskQA();
  AliHMPIDTaskQA(const Char_t* name);
  AliHMPIDTaskQA& operator= (const AliHMPIDTaskQA& c);
  AliHMPIDTaskQA(const AliHMPIDTaskQA& c);
  virtual ~AliHMPIDTaskQA();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

          void   SetUseMC(Bool_t useMC) { fUseMC = useMC; }
          Bool_t Equal(Double_t x, Double_t y, Double_t tolerance);

 protected:
     
 private:     
 
  AliESDEvent *fESD;               //! ESD object
  AliMCEvent  *fMC;                //! MC event

  Bool_t       fUseMC;             // decide whether use or not the MC information

  TList         *fHmpHistList ;    // list of histograms

  TH1F          *fHmpMipTrkDistX;  // Xtrk - Xmip
  TH1F          *fHmpMipTrkDistY;  // Ytrk - Ymip
  TH1F          *fHmpMipCharge3cm; // Mip charge with 3 cm distance cut
  TH1F          *fHmpTrkFlags;     // track flags

  Int_t          fN1;              // number of points for pi and K
  Int_t          fN2;              // number of point for p
  TH1F          *fPionEff;         // identified pions
  TH1F          *fKaonEff;         // identified kaons
  TH1F          *fProtEff;         // identified protons
  TH1I          *fPionTot;         // total pions
  TH1I          *fKaonTot;         // total kaons
  TH1I          *fProtTot;         // total protons
  TH1F          *fPionNot;         // non-pion tracks
  TH1F          *fKaonNot;         // non-kaon tracks
  TH1F          *fProtNot;         // non-proton tracks
  TH1I          *fPionCon;         // tracks identified as pions
  TH1I          *fKaonCon;         // tracks identified as kaons
  TH1I          *fProtCon;         // tracks identified as protons

  TTree         *fTree;            // tree with useful data for subsequent analysis
  Float_t        fVar[28];         // array of data to fill the tree

  ClassDef(AliHMPIDTaskQA,1);
};

#endif