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

#include "AliHFEcuts.h"
#include <AliCFContainer.h>
#include <AliCFDataGrid.h>

class TArrayD;
class TParticle;
class AliAODMCParticle;
class AliPIDResponse;
class AliVParticle;
class TGraph;
class TGraphErrors;
class TGraphAsymmErrors;
class TH1D;
class TH1;
class TF1;
class TString;

class AliHFEtools : public TObject{
  public:
    AliHFEtools();
    ~AliHFEtools() {};

    static Double_t *MakeLinearBinning(Int_t nBins, Double_t ymin, Double_t ymax);
    static Double_t *MakeLogarithmicBinning(Int_t nBins, Double_t ymin, Double_t ymax);
    static void FillLinearBinning(TArrayD &bins, Int_t nBins, Double_t ymin, Double_t ymax);
    static void FillLogarithmicBinning(TArrayD &bins, Int_t nBins, Double_t ymin, Double_t ymax);
    Bool_t    BinLogAxis(TObject *o, Int_t dim);
    static Float_t GetRapidity(const TParticle *part);
    static Float_t GetRapidity(const AliAODMCParticle *part); // return rapidity
    static Int_t GetPdg(const AliVParticle *track);
    static Int_t PDG2AliPID(Int_t pdg);
    static AliPIDResponse *GetDefaultPID(Bool_t isMC = kTRUE, Bool_t isESD = kTRUE);
    static void DestroyDefaultPID();
    static void SetLogLevel(Int_t loglevel) { fgLogLevel = loglevel ;}
    static TH1D* GraphErrorsToHist(TGraphErrors* g = 0, Double_t firstBinWidth = -1, Bool_t exchange=kFALSE, Int_t markerstyle=8, Int_t markercolor=2, Float_t markersize=0.7);
    static TH1D* GraphToHist(TGraph* g = 0, Double_t firstBinWidth = -1, Bool_t exchange=kFALSE, Int_t markerstyle=8, Int_t markercolor=2, Float_t markersize=0.7);
    static Bool_t ExchangeXYGraph(TGraph* g = 0);
    static Bool_t ExchangeXYGraphErrors(TGraphErrors* g = 0);
    static void BinParameterisation(const TF1 &fun, const TArrayD &xbins, TArrayD &bincontent);
    static TList *GetHFEResultList(const TString str);
    static TList *GetHFEResultListGRID(const TString str, const TString dir);

    static TList *GetHFEQAList(const TString str);
    static TList *GetHFEQAListGRID(const TString str, const TString dir);
    static void NormaliseBinWidth(TH1 *histo);
    static void NormaliseBinWdith(TGraphErrors *graph);
    static void NormaliseBinWdithAsymm(TGraphAsymmErrors *graph);
    static TGraphErrors *Normalise(const TH1 * const input, Int_t fNEvents, Double_t fEtaRange, Int_t fNCharges = 2);
    static Int_t GetNumberOfEvents(const TString filename, Double_t centmin = -1., Double_t centmax = 11.);

  private:
      AliHFEtools(const AliHFEtools &);
      AliHFEtools &operator=(const AliHFEtools &);
      static AliPIDResponse *fgDefaultPID;      // Default PID object
      static Int_t fgLogLevel;                  // Log Level

    ClassDef(AliHFEtools, 0)
};
#endif
