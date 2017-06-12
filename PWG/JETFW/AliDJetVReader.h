/**
 * \file AliDJetVReader.h
 * \brief Declaration of class AliDJetVReader
 *
 * In this header file the class AliDJetVReader is declared.
 * Class to read the invariant mass histograms used to extract the raw yield.
 *
 * \author Fabio Colamaria <fabio.colamaria@cern.ch>, INFN Bari
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Mar 7, 2017
 */

#ifndef ALIDJETVREADER_H
#define ALIDJETVREADER_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>

/**
 * \class AliDJetVReader
 * \brief Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 *
 * Implementation of an abstract class to read the invariant mass histograms used to extract the raw yield.
 */
class AliDJetVReader : public TObject {

public:

  AliDJetVReader();
  AliDJetVReader(const AliDJetVReader &source);
  virtual ~AliDJetVReader();

  void SetPtBinEdgesForMassPlot(Double_t ptmin, Double_t ptmax)    { fpTmin                  = ptmin ; fpTmax                  = ptmax ; }
  void SetZedges(Double_t zmin, Double_t zmax)                     { fzmin                   = zmin  ; fzmax                   = zmax  ; }
  void SetMassRebin(Int_t r)                                       { fMassRebin              = r > 0 ? r : 1; }

  void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetDmesonEfficiency(Double_t* effvalues=0x0);

  virtual Bool_t ExtractInputMassPlotEffScale() = 0;
  virtual Bool_t ExtractInputMassPlotSideband() = 0;

  TH1D* GetMassPlot()         { return fMassPlot         ; }
  TH2D* GetMassVsJetPtPlot()  { return fMassVsJetPtPlot  ; }

protected:
  Double_t           fpTmin                      ; ///< pT lower edge of mass plot to evaluate variations of yields
  Double_t           fpTmax                      ; ///< pT upper edge of mass plot to evaluate variations of yields
  Double_t           fzmin                       ; ///< z minimum value to extract jet pT spectrum
  Double_t           fzmax                       ; ///< z maximum value to extract jet pT spectrum
  Int_t              fnDbins                     ; ///< Number of D-meson pT bins (for eff scaling)
  Double_t          *fDbinpTedges                ; ///< D-meson pt bin edges values
  Int_t              fnJetbins                   ; ///< Number of pT-bins to be used for spectrum
  Double_t          *fJetbinpTedges              ; ///< Jet pT bin edges to be used for spectrum
  Double_t          *fDEffValues                 ; ///< D-meson efficiency values
  UInt_t             fMassRebin                  ; ///< Rebin the mass histogram axis

  TH1D              *fMassPlot                   ; //!<!Mass spectra to be fitted
  TH2D              *fMassVsJetPtPlot            ; //!<!Mass vs jet pt (SB method)

private:
  ClassDef(AliDJetVReader,2);
};

#endif
