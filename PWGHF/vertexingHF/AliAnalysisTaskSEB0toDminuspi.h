#ifndef ALIANALYSISTASKSEB0toDminuspi_H
#define ALIANALYSISTASKSEB0toDminuspi_H
/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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



#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisTaskSE.h"

class AliRDHFCutsDplustoKpipi;
class AliNormalizationCounter;

class AliAnalysisTaskSEB0toDminuspi : public AliAnalysisTaskSE 
{
  
 public:
  
  AliAnalysisTaskSEB0toDminuspi();
  AliAnalysisTaskSEB0toDminuspi(const Char_t* name, AliRDHFCutsDplustoKpipi* dpluscutsana);
  virtual ~AliAnalysisTaskSEB0toDminuspi();

  AliAODVertex* RecalculateVertex(const AliVVertex *old,TObjArray *tracks,Double_t bField);

  void SetCutsDistr(Bool_t cutsDistr=kTRUE){fCutsDistr=cutsDistr;}

  // Implementation of interface methods  
  virtual void UserCreateOutputObjects();
  virtual void Init();
  virtual void LocalInit() {Init();}
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
 

    // histos
  void FillSpectrum(Int_t nB0, Int_t* arrayMClabel, Double_t Bfield);
  void     DefineHistograms();
  
  
  void SetUseBit(Bool_t dols=kTRUE){fUseBit=dols;}
  void SetUseOnlyPositiveEta(){fEtaSelection=1;}
  void SetUseOnlyNegativeEta(){fEtaSelection=-1;}
  void SetUseFullEta(){fEtaSelection=0;}
  void SetReadMC(Bool_t readMC=kTRUE){fReadMC=readMC;}
  void SetUseQuarkLevelTag(Bool_t opt){fUseQuarkTagInKine=opt;}
  
  void SetTopomaticCutOnDDaughters(Double_t TopoCut){fInvertedTopomaticCutOnDDaughters=TopoCut;}
  void SetDpTCut(Double_t pTD){fpTD=pTD;} 
  void SetpipTCut(Double_t pTpi){fpTpi=pTpi;} 
  void Setprodd0Cut(Double_t prodd0){fprodd0=prodd0;} 
  void SetCosBCut(Double_t cosB){fcosB=cosB;} 
  void SetCosXYBCut(Double_t cosBXY){fcosBXY=cosBXY;} 
  void SetdlBCut(Double_t dlB){fdlB=dlB;} 
  void SetNdlXYBCut(Double_t NdlBXY){fNdlBXY=NdlBXY;} 
  void SetTopomaticD(Double_t TopomaticD){fTopomaticD=TopomaticD;} 
  void SetTopomaticpi(Double_t Topomaticpi){fTopomaticpi=Topomaticpi;} 
  void SetcosoaDpi(Double_t cosoaDpi){fcosoaDpi=cosoaDpi;}
  
  //to define the range of the D plus mass candidates to form a B meson

  void SetDplusMassLimit(Double_t DplusMassLowLimit,Double_t DplusMassUpLimit){
    fDplusMassUpLimit=DplusMassUpLimit;
    fDplusMassLowLimit=DplusMassLowLimit;}

  //range of the B meson invariant mass distributions
  void SetBMassLimits(Float_t lowlimit, Float_t uplimit);
  void SetBinWidth(Float_t w);


  Float_t GetUpperBMassLimit(){return fBMassUpLimit;}
  Float_t GetLowerBMassLimit(){return fBMassLowLimit;}
  Float_t GetBinWidth(){return fBinWidth;}
  Int_t GetNBinsHistos();

  
  //   void SetBMassLimitForHisto(Double_t BMassLowLimit,Double_t BMassUpLimit){
  //  fBMassLowLimit=BMassLowLimit;fBMassUpLimit=BMassUpLimit;}


  

 private:
  
  AliAnalysisTaskSEB0toDminuspi(const AliAnalysisTaskSEB0toDminuspi &source);
  AliAnalysisTaskSEB0toDminuspi& operator=(const AliAnalysisTaskSEB0toDminuspi& source); 
  
  enum {kMaxPtBins=10};
 
  
  
  TList *fOutput;                //!<!  User output
  TList *fOutputBins;             //!<!  User output for different pt bins
  TList *fListCuts; //!<! list of cuts
  
  AliRDHFCutsDplustoKpipi *fRDCutsAnalysis; /// Cuts for Analysis
  std::vector<AliAODRecoDecayHF2Prong*> fArrayB0; //!<! vector to save B candidate 
 
  // define the histograms
  TH1F *fCEvents;                            //!<!hist. for No. of events
  TH1F *fB0InvMass;                          //!<!hist. for B0 mass 
  TH1F *fDplusInvMass;                       //!<!hist. for D+ mass 
  TH1F* fB0InvMass_PassCut1;                 //!<!hist. for B0 mass after cut
  TH1F* fB0InvMassMinusD_PassCut1;           //!<!hist. for B0minusD mass after cut   
  TH1F* fDplusInvMass_PassCut1;              //!<!hist. for D+ mass after cut 
  TH1F* fDplusInvMass_DplusMCmatch;          //!<!hist. for D+ mass matched to B0 in MC 
  TH1F* fDplusInvMass_MCmatch;               //!<!hist. for  (after cut) matched to B0 in MC 
  TH1F* fBInvMass_MCmatch;                   //!<!hist. for B0 mass matched to B0 in MC 
  TH1F* fBInvMass_MCmatch_PassCut1;          //!<!hist. for B0 mass after cut matched to B0 in MC 
  TH1F* fB0InvMassMinusD_MCmatch_PassCut1;   //!<!hist. for B0minusD mass after cut matched to B0 in MC 
  TH1F* fpiPt;                //!<!hist. for pion Pt (after cut)
  TH1F* fpiEta;               //!<!hist. for pion eta (after cut)
  TH1F* fpiPhi;               //!<!hist. for pion phi (after cut)
  TH1F* fDplusPt;             //!<!hist. for D Pt (after cut)
  TH1F* fDplusEta;            //!<!hist. for D eta (after cut)
  TH1F* fDplusPhi;            //!<!hist. for D phi (after cut)
  TH1F* fDplusPointing;       //!<!hist. for D cos pointing (after cut)
  TH1F* fDplusDecayLength;    //!<!hist. for D decay length (after cut)
  TH1F* fDplusNormalizedDecayLength; //!<!hist. for D norm decay lengthXY (after cut)
  TH1F* fDplusSigmaVert;      //!<!hist. for D sigma vertex (after cut)
  TH1F* fBPt;                 //!<!hist. for B Pt (after cut)
  TH1F* fBEta;                //!<!hist. for B eta (after cut)
  TH1F* fBPhi;                //!<!hist. for B phi  (after cut)
  TH1F* fBPointing;           //!<!hist. for B cos poining (after cut)
  TH1F* fBPointingXY;         //!<!hist. for B cos poining XY (after cut)
  TH1F* fBDecayLength;        //!<!hist. for B decay length (after cut)
  TH1F* fBNormalizedDecayLength; //!<!hist. for B norm decay lengthXY  (after cut)
  TH1F* fDplusd0;             //!<!hist. for D d0 (after cut)
  TH1F* fpid0;                //!<!hist. for pion d0  (after cut)
  TH1F* fproductd0;           //!<!hist. for d0 product (after cut)
  TH1F* fpiPt_MCmatch;        //!<!hist. for pi Pt (after cut) matched to B0 in MC
  TH1F* fpiEta_MCmatch;       //!<!hist. for pi eta (after cut) matched to B0 in MC
  TH1F* fpiPhi_MCmatch;       //!<!hist. for pi phi (after cut) matched to B0 in MC
  TH1F* fDplusPt_MCmatch;     //!<!hist. for D pT (after cut) matched to B0 in MC
  TH1F* fDplusEta_MCmatch;    //!<!hist. for D eta (after cut) matched to B0 in MC
  TH1F* fDplusPhi_MCmatch;    //!<!hist. for D phi (after cut) matched to B0 in MC
  TH1F* fDplusPointing_MCmatch; //!<!hist. for D cos pointing (after cut) matched to B0 in MC
  TH1F* fDplusDecayLength_MCmatch; //!<!hist. for D decay length (after cut) matched to B0 in MC
  TH1F* fDplusNormalizedDecayLength_MCmatch;             //!<!hist. for D norm dec length XY(after cut) matched to B0 in MC
  TH1F* fDplusSigmaVert_MCmatch;  //!<!hist. for D sigma vertex (after cut) matched to B0 in MC
  TH1F* fBPt_MCmatch;             //!<!hist. for B pt (after cut) matched to B0 in MC
  TH1F* fBEta_MCmatch;            //!<!hist. for B eta (after cut) matched to B0 in MC
  TH1F* fBPhi_MCmatch;            //!<!hist. for B phi (after cut) matched to B0 in MC
  TH1F* fBPointing_MCmatch;       //!<!hist. for B pointing (after cut) matched to B0 in MC
  TH1F* fBPointingXY_MCmatch;     //!<!hist. for B pointing XY (after cut) matched to B0 in MC
  TH1F* fBDecayLength_MCmatch;    //!<!hist. for B decay length (after cut) matched to B0 in MC
  TH1F* fBNormalizedDecayLength_MCmatch; //!<!hist. for B norm dec length XY (after cut) matched to B0 in MC
  TH1F* fDplusd0_MCmatch;         //!<!hist. for D d0 (after cut) matched to B0 in MC
  TH1F* fpid0_MCmatch;            //!<!hist. for pion d0 (after cut) matched to B0 in MC 
  TH1F* fproductd0_MCmatch;       //!<!hist. for d0 product (after cut) matched to B0 in MC 


  TH1F* fB0InvMassBin[kMaxPtBins];                        //!<!hist. for B0 mass in pT bins
  TH1F* fDplusInvMassBin[kMaxPtBins];                     //!<!hist. for D+ mass in pT bins
  TH1F* fB0InvMassBin_MCmatch[kMaxPtBins];                //!<!hist. for B0 mass matched to B0 in MC in pT bins 
  TH1F* fDplusInvMassBin_MCmatch[kMaxPtBins];             //!<!hist. for D+ mass matched to B0 in MC in pT bins  
  TH1F* fB0InvMassBin_PassCut1[kMaxPtBins];               //!<!hist. for B0 mass after cut in pT bins   
  TH1F* fDplusInvMassBin_PassCut1[kMaxPtBins];            //!<!hist. for D mass after cut matched to B0 in MC in pT bins  
  TH1F* fB0InvMassBinMinusD_PassCut1[kMaxPtBins];         //!<!hist. for B0minusD mass after cut in MC in pT bins  
  TH1F* fB0InvMassBin_MCmatch_PassCut1[kMaxPtBins];       //!<!hist. for B0 mass after cut matched to B0 in MC in pT bins 
  TH1F* fB0InvMassBinMinusD_MCmatch_PassCut1[kMaxPtBins]; //!<!hist. for B0minusD mass after cut matched to B0 in MC in pT bins 
  TH1F* fBDecayLengthBin[kMaxPtBins];                     //!<!hist. for B decay length (after cut) in pT bins
  TH1F* fBNormalizedDecayLengthBin[kMaxPtBins];           //!<!hist. for B norm decay length XY (after cut) in pT bins
  TH1F* fBPointingBin[kMaxPtBins];                        //!<!hist. for B cos pointing (after cut) in pT bins
  TH1F* fBPointingXYBin[kMaxPtBins];                       //!<!hist. for B cos pointing XY (after cut) in pT bins
  TH1F* fDplusd0Bin[kMaxPtBins];                          //!<!hist. for D d0  (after cut) in pT bins
  TH1F* fBpTBin[kMaxPtBins];                              //!<!hist. for B pT (after cut) in pT bins
  TH1F* fDpluspTBin[kMaxPtBins];                          //!<!hist. for D pT (after cut) in pT bins
  TH1F* fproductd0Bin[kMaxPtBins];                        //!<!hist. for d0 product (after cut) in pT bins
  TH1F* fDplusPointingBin[kMaxPtBins];                    //!<!hist. for D cos pointing (after cut) in pT bins
  TH1F* fDplusDecayLengthBin[kMaxPtBins];                 //!<!hist. for D decay length (after cut) in pT bins
  TH1F* fDplusNormalizedDecayLengthBin[kMaxPtBins];       //!<!hist. for D norm decay length XY (after cut) in pT bins
  TH1F* fDplusSigmaVertBin[kMaxPtBins];                   //!<!hist. for D sigma vertex (after cut) in pT bins
  TH1F* fBDecayLengthBin_MCmatch[kMaxPtBins];             //!<!hist. for B decay length (after cut) in pT bins
  TH1F* fBNormalizedDecayLengthBin_MCmatch[kMaxPtBins];   //!<!hist. for B norm decay length XY (after cut) in pT bins
  TH1F* fBPointingBin_MCmatch[kMaxPtBins];                //!<!hist. for B cos pointnig (after cut) matched to B0 in MC in pT bins
  TH1F* fBPointingXYBin_MCmatch[kMaxPtBins];              //!<!hist. for B cos pointnig XY (after cut) matched to B0 in MC in pT bins
  TH1F* fDplusd0Bin_MCmatch[kMaxPtBins];                  //!<!hist. for D d0 (after cut) matched to B0 in MC in pT bins
  TH1F* fBpTBin_MCmatch[kMaxPtBins];                      //!<!hist. for B pT (after cut) matched to B0 in MC in pT bins
  TH1F* fDpluspTBin_MCmatch[kMaxPtBins];                  //!<!hist. for D pT (after cut) matched to B0 in MC in pT bins
  TH1F* fproductd0Bin_MCmatch[kMaxPtBins];                //!<!hist. for product d0 (after cut) matched to B0 in MC in pT bins
  TH1F* fDplusPointingBin_MCmatch[kMaxPtBins];            //!<!hist. for D pointing (after cut) matched to B0 in MC in pT bins
  TH1F* fDplusDecayLengthBin_MCmatch[kMaxPtBins];         //!<!hist. for D decay length (after cut) matched to B0 in MC in pT bins
  TH1F* fDplusNormalizedDecayLengthBin_MCmatch[kMaxPtBins];//!<!hist. for D norm dec length XY (after cut) matched to B0 in MC in pT bins
  TH1F* fDplusSigmaVertBin_MCmatch[kMaxPtBins];      //!<!hist. for D sigma vertex (after cut) matched to B0 in MC in pT bins

  TH1F* fd0MMExpD;                                  //!<!hist. for D topomatic (after cut)
  TH1F* fd0MMExpDMCmatch;                           //!<!hist. for D topomatic (after cut) matched to B0 in MC 
  TH1F* fd0MMExppi;                                 //!<!hist. for pion topomatic (after cut)
  TH1F* fd0MMExppiMCmatch;                          //!<!hist. for pion topomatic (after cut) matched to B0 in MC 
  TH1F* fd0MMExpDBin[kMaxPtBins];                   //!<!hist. for D topomatic (after cut) in pT bins
  TH1F* fd0MMExpDBin_MCmatch[kMaxPtBins];           //!<!hist. for D topomatic (after cut) matched to B0 in MC  in pT bins
  TH1F* fd0MMExppiBin[kMaxPtBins];                  //!<!hist. for pion topomatic (after cut) in pT bins
  TH1F* fd0MMExppiBin_MCmatch[kMaxPtBins];          //!<!hist. for pion topomatic (after cut) matched to B0 in MC  in pT bins
  TH1F* fd0MMExpDDaughters;                         //!<!hist. for D daughters topomatic (after cut)
  TH1F* fd0MMExpDDaughters_MCmatch;                 //!<!hist. for D daughters topomatic (after cut) matched to B0 in MC
  TH1F* fd0MMExpDDaughtersBin[kMaxPtBins];          //!<!hist. for D daughters topomatic (after cut) in pT bins
  TH1F* fd0MMExpDDaughtersBin_MCmatch[kMaxPtBins];  //!<!hist. for D daughters topomatic (after cut) matched to B0 in MC in pT bins
  TH1F* fcosoa;                                     //!<!hist. for cos opening angle between D and pion (after cut)
  TH1F* fcosoa_MCmatch;                             //!<!hist. for cos opening angle between D and pion (after cut) matched to B0 in MC
 
  TH2F* fBDplusPt;           //!<!hist of B pT vs D pT 
  TH2F* fBDplusPt_MCmatch;   //!<!hist of B pT vs D pT matched to B0 in MC
  TH2F* fBpiPt;              //!<!hist of B pT vs pion pT 
  TH2F* fBpiPt_MCmatch;      //!<!hist of B pT vs pion pT matched to B0 in MC
  TH2F* fBDplusPtBin[kMaxPtBins];           //!<!hist of B pT vs D pT in pT bins
  TH2F* fBDplusPtBin_MCmatch[kMaxPtBins];   //!<!hist of B pT vs D pT matched to B0 in MC in pT bins
  TH2F* fBpiPtBin[kMaxPtBins];              //!<!hist of B pT vs pion pT in pT bins
  TH2F* fBpiPtBin_MCmatch[kMaxPtBins];      //!<!hist of B pT vs pion pT matched to B0 in MC in pT bins

  AliNormalizationCounter *fCounter; //!Counter for normalization slot 


  Int_t  fEvents;       ///  n. of events
  Bool_t fUseBit;       /// flag to use bitmask  
  Int_t fEtaSelection;  /// eta region to accept D+ 0=all, -1 = negative, 1 = positive
  Bool_t fReadMC;       /// flag for access to MC 
  Bool_t kbins = kTRUE; /// flag to add histo in pT bins
  Bool_t fUseQuarkTagInKine; /// flag for quark/hadron level identification of prompt and feeddown
  Bool_t fCutsDistr;    /// flag to activate cuts distr histos
  
  Double_t fDplusMassLowLimit; ///DplusMassLowLimit to select D
  Double_t fDplusMassUpLimit;  ///DplusMassUpLimit to select D
  Double_t fBMassLowLimit;     ///BMassLowLimit for histo 
  Double_t fBMassUpLimit;      ///BMassUpLimit for histo 
  Float_t fBinWidth;           /// width of one bin in B mass output histos
 
 

  Double_t fInvertedTopomaticCutOnDDaughters;///only candidates with a minimum topomatic of fUseInvertedTopomaticCutOnDDaughters are considered
  Double_t fpTD;   ///cut to applied on B selection: pT D
  Double_t fpTpi;  ///cut to applied on B selection: pT pion
  Double_t fprodd0;///cut to applied on B selection: product d0
  Double_t fcosB;  ///cut to applied on B selection: B cos pointing 
  Double_t fcosBXY;  ///cut to applied on B selection: B cos pointing XY
  Double_t fdlB;   ///cut to applied on B selection: B decay length
  Double_t fNdlBXY; ///cut to applied on B selection: B norm decay length XY
  Double_t fTopomaticD; ///cut to applied on B selection: D topomatic
  Double_t fTopomaticpi;///cut to applied on B selection: pion topomatic
  Double_t fcosoaDpi;   ///cut to applied on B selection: opening angke
  
  

  ClassDef(AliAnalysisTaskSEB0toDminuspi,2); /// class for B spectra
};

#endif

 
