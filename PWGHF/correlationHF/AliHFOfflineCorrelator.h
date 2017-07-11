#ifndef AliHFOfflineCorrelator_H
#define AliHFOfflineCorrelator_H

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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
//  Base class to perform D-h correlations offline (starting from D-meson and hadron TTrees)
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
//-----------------------------------------------------------------------

#include <iostream>
#include "TObject.h"
#include "TMath.h"
#include "TFile.h"
#include "TDirectoryFile.h"
#include "TList.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TDatabasePDG.h"
#include "AliHFAssociatedTrackCuts.h"

using std::vector;

class AliHFCorrelationBranchD : public TObject
{
  public:
    AliHFCorrelationBranchD(); // default constructor
    Float_t  phi_D;
    Float_t  eta_D;
    Float_t  pT_D;
    Float_t  mult_D;
    Float_t  zVtx_D;
    Float_t  cent_D;
    Float_t  invMass_D;
    UInt_t   period_D;
    UInt_t   orbit_D;
    UShort_t BC_D;
    Short_t  IDtrig_D;
    Short_t  sel_D;
    Float_t  pXdaug1_D;
    Float_t  pXdaug2_D;
    Float_t  pYdaug1_D;
    Float_t  pYdaug2_D;
    Float_t  pZdaug1_D;
    Float_t  pZdaug2_D;
    UShort_t hyp_D; //1 or 2 (no 3, because double hypotheses are filled as separate entries)

    ClassDef(AliHFCorrelationBranchD,2);
};

class AliHFCorrelationBranchTr : public TObject
{
  public:
    AliHFCorrelationBranchTr(); // default constructor
    Float_t  phi_Tr;
    Float_t  eta_Tr;
    Float_t  pT_Tr;
    Float_t  mult_Tr;
    Float_t  zVtx_Tr;
    Float_t  cent_Tr;
    UInt_t   period_Tr;
    UInt_t   orbit_Tr;
    UShort_t BC_Tr;
    Short_t  IDtrig_Tr;
    Short_t  IDtrig2_Tr;
    Short_t  IDtrig3_Tr;
    Short_t  IDtrig4_Tr;
    Short_t  sel_Tr;

    ClassDef(AliHFCorrelationBranchTr,3);
};

class AliHFOfflineCorrelator : public TObject
{

public:
    
    enum DMesonSpecies {kD0toKpi, kDplusKpipi, kDStarD0pi};
    enum AnalysisType {kSE, kME};

    AliHFOfflineCorrelator(); // default constructor
    AliHFOfflineCorrelator(const AliHFOfflineCorrelator &source);
    AliHFOfflineCorrelator& operator=(const AliHFOfflineCorrelator& source); 
    virtual ~AliHFOfflineCorrelator();

    Bool_t SetDmesonSpecie(DMesonSpecies k);
    void SetAnalysisType(AnalysisType k) {fAnType=k;}
    void SetUseEfficiency(Bool_t eff) {fUseEff=eff;}
    void AddInputFile(TString str) {fFileList.push_back(str); fNinputFiles++;}
    void SetDirName(TString dirname) {fDirName=dirname;}
    void SetTreeNames(TString nameD, TString nameTr) {fNameTreeD=nameD; fNameTreeTr=nameTr;}
    void SetEffMapNames(TString cutObj, TString nameD, TString nameTr) {fNameCutObj=cutObj; fNameMapD=nameD; fNameMapTr=nameTr;}
    void SetDPtBins(Int_t nBins, Double_t* ptDarray);  
    void AddAssocPtRange(Double_t lowPt, Double_t upPt) {fPtBinsTrLow.push_back(lowPt); fPtBinsTrUp.push_back(upPt);}
    void MakeDetaDphiPlots(Double_t* MsigL, Double_t* MsigR, Double_t* MSB1L, Double_t* MSB1R, Double_t* MSB2L=0x0, Double_t* MSB2R=0x0); //makes 2D plots for sign.region and SB
    void SetOutputFileName(TString namefile) {fOutputFileName=namefile;}
    void SetPoolBins(Int_t nMultPools, Double_t* multBins, Int_t nzVtxPools, Double_t* zVtxBins);
    void SetMaxTracksToCorrelate(Int_t maxTracks=-1) {fMaxTracks=maxTracks;}
    void SetDLoopRange(Int_t minD=-1, Int_t maxD=-1) {fMinD=minD; fMaxD=maxD;}
    void SetWeightPeriods(Bool_t wgh) {fWeightPeriods=wgh;}
    void SetFirstBinNum(Int_t num=0) {fFirstBinNum=num;}
    void SetNumSelD(Int_t sel) {fNumSelD=sel;}
    void SetNumSelTr(Int_t sel) {fNumSelTr=sel;}
    void SetCentralitySelection(Double_t min, Double_t max) {fMinCent=min; fMaxCent=max;} //activated only if both values are != 0
    void SetRejectSoftPion(Bool_t store) {fRejectSoftPi=store;}
    void SetDebugLevel(Int_t deb=0) {fDebug=deb;}

    Bool_t Correlate();

    void DefineOutputObjects();
    void PrintCfg() const;
    Bool_t CorrelateSingleFile(Int_t iFile);
    void GetCorrelationsValue(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr, Double_t &deltaPhi, Double_t &deltaEta);
    Double_t GetEfficiencyWeight(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr);
    Double_t GetEfficiencyWeightDOnly(AliHFCorrelationBranchD *brD);
    Bool_t IsSoftPionFromDstar(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr);
    Int_t PtBin(Double_t pt) const;
    Int_t GetPoolBin(Double_t mult, Double_t zVtx) const;
    Bool_t DefinePeriodWeights();
    void NormalizeMEPlots();
    void SaveOutputPlots();

private:
    
    std::vector<TString>  fFileList;    //container of input filenames
    Int_t fNinputFiles;			//number of input files

    TFile *fFile;        		//file containing the analysis output

    TTree *fTreeD;    			//for D Tree allocation
    TTree *fTreeTr;    			//for associated track Tree allocation
    
    TList *fOutputDistr;		//allocates the outputcorrelation plots
    TList *fOutputMass;			//allocates the mass plots

    Int_t fNBinsPt;			//number of D-meson pT bins
    Int_t fnPools;			//number of ME pools (total)
    Int_t fnMultPools;			//number of ME pools (multiplicity)
    Int_t fnzVtxPools;			//number of ME pools (z vertex)
    Int_t fFirstBinNum;			//number of first bin for the name of the output plots
    Int_t fNumSelD;			//number of selection for D meson (0=default selection; 1,2,3,... = alternate selections)
    Int_t fNumSelTr;			//number of selection for assoc tracks (0=default selection; 1,2,3,... = alternate selections)
    Int_t fDebug;

    Int_t fMinD;			//start of loop on D mesons
    Int_t fMaxD;			//end of loop on D mesons
    Int_t fMaxTracks;			//maximum number of tracks (of the tree, i.e. integrated over everything) to be used in ME correlations
    Double_t fMinCent;			//minimum centrality 
    Double_t fMaxCent;			//maximum centrality 

    std::vector<Double_t>  fPtBinsDLow;	//lower edges of D-meson pT bins
    std::vector<Double_t>  fPtBinsDUp; 	//upper edges of D-meson pT bins
    std::vector<Double_t>  fPtBinsTrLow;//lower edges of assoc. track pT bins
    std::vector<Double_t>  fPtBinsTrUp;	//upper edges of assoc. track pT bins
    std::vector<Double_t>  fMassSignL;  //lower mass edge of signal region
    std::vector<Double_t>  fMassSignR;  //lower mass edge of signal region
    std::vector<Double_t>  fMassSB1L;   //lower mass edge of signal region
    std::vector<Double_t>  fMassSB1R;   //lower mass edge of signal region
    std::vector<Double_t>  fMassSB2L;   //lower mass edge of signal region
    std::vector<Double_t>  fMassSB2R; 	//lower mass edge of signal region
    std::vector<Double_t>  fMultBins; 	//bin edges for ME pools (multiplicity)
    std::vector<Double_t>  fzVtxBins; 	//bin edges for ME pools (z vertex)
    std::vector<Double_t>  fPrdWeights;	//period weights (if different number of tracks is used, in ME analysis)

    TH3F* fMapEffTr;
    TH2F* fMapEffD;

    DMesonSpecies fDmesonSpecies;
    AnalysisType fAnType;

    TString fDmesonLabel;
    TString fDirName; 
    TString fNameTreeTr;
    TString fNameTreeD;
    TString fNameMapTr;
    TString fNameMapD;
    TString fOutputFileName;
    TString fNameCutObj;

    Bool_t fUseEff;			    //flag to use D-meson and track efficiency
    Bool_t fMake2DPlots; 		//flag to produce 2D plots for sign.region and SB
    Bool_t fWeightPeriods;		//flag to weight periods in ME analysis with max number of tracks used
    Bool_t fRejectSoftPi;	     //flag to remove soft pions in SE and ME analysis for D0 meson (ME rejection is done in extraction code)

    ClassDef(AliHFOfflineCorrelator,4); // class for plotting HF correlations

};

#endif

