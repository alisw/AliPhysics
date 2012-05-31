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
// AliHMPIDAnalysysTask - Class representing a basic analysis tool of HMPID data at  
// level of ESD.
// A set of histograms is created.
//==============================================================================

#ifndef ALIHMPIDANALYSISTASK_H
#define ALIHMPIDANALYSISTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class TH1;
class TParticle;
class TFile;
class AliESDtrack;
class AliESDEvent;

class AliHMPIDAnalysisTask : public AliAnalysisTaskSE {
  public:

  enum {kChamber = 7};

  AliHMPIDAnalysisTask();
  AliHMPIDAnalysisTask(const Char_t* name);
  AliHMPIDAnalysisTask& operator= (const AliHMPIDAnalysisTask& c);
  AliHMPIDAnalysisTask(const AliHMPIDAnalysisTask& c);
  virtual ~AliHMPIDAnalysisTask();
  
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

  TH2F          *fHmpPesdPhmp;     // HMP momentum vs ESD momentum
  TH2F          *fHmpCkovPesd;     // Ckov angle vs ESD momentum
  TH2F          *fHmpCkovPhmp;     // Ckov angle vs HMP momenutm

  TH1F          *fHmpMipTrkDist;   // Track-Mip distance distribution
  TH1F          *fHmpMipTrkDistX;  // Xtrk - Xmip
  TH1F          *fHmpMipTrkDistY;  // Ytrk - Ymip
  TH1F          *fHmpMipCharge3cm; // Mip charge with 3 cm distance cut
  TH1F          *fHmpMipCharge1cm; // Mip charge with 1 cm distance cut
  TH1F          *fHmpNumPhots;     // Number of reconstructed photo-electrons
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
  TH2F          *fThetavsPiFromK;  // theta chkov of pis from Ks
  TH2F          *fThetapivsPesd;   // theta chkov of pions vs Pesd
  TH2F          *fThetaKvsPesd;    // theta chkov of kaons vs Pesd
  TH2F          *fThetaPvsPesd;    // theta chkov of protons vs Pesd
  TH1I          *fProtGen;         // generated protons from MC
  TH1I          *fPbarGen;         // generated antiprotons from MC
  TH1I          *fProtHmp;         // protons into the HMPID
  TH1I          *fPbarHmp;         // antiprotons into the HMPID

  TTree         *fTree;            // tree with useful data for subsequent analysis
  Float_t        fVar[28];         // array of data to fill the tree

  ClassDef(AliHMPIDAnalysisTask,4);
};

#endif
