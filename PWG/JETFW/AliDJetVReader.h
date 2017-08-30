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
  void SetMassRebin(UInt_t r)                                      { fMassRebin              = r > 0 ? r : 1; }

  UInt_t GetMassRebin() const { return fMassRebin; }

  void SetDmesonPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetJetPtBins(Int_t nbins=0, Double_t* ptedges=0x0);
  void SetJetzBins(Int_t nbins=0, Double_t* zedges=0x0);
  void SetDmesonEfficiency(Double_t* effvalues=0x0);

  virtual Bool_t ExtractInputMassPlotEffScale() = 0;
  virtual Bool_t ExtractInputMassPlotSideband() = 0;

  TH1D* GetMassPlot()         { return fMassPlot         ; }
  TH2D* GetMassVsJetPtPlot()  { return fMassVsJetPtPlot  ; }
  TH2D* GetMassVsJetzPlot()   { return fMassVsJetzPlot   ; }

protected:
  Double_t           fpTmin                      ; ///< pT lower edge of mass plot to evaluate variations of yields
  Double_t           fpTmax                      ; ///< pT upper edge of mass plot to evaluate variations of yields
  Int_t              fnDbins                     ; ///< Number of D-meson pT bins (for eff scaling)
  Double_t          *fDbinpTedges                ; ///< D-meson pt bin edges values
  Int_t              fnJetPtbins                 ; ///< Number of jet pT bins to be used for spectrum
  Double_t          *fJetPtBinEdges              ; ///< Jet pT bin edges to be used for spectrum
  Int_t              fnJetzbins                  ; ///< Number of jet z bins to be used for spectrum
  Double_t          *fJetzBinEdges               ; ///< Jet z bin edges to be used for spectrum
  Double_t          *fDEffValues                 ; ///< D-meson efficiency values
  UInt_t             fMassRebin                  ; ///< Rebin the mass histogram axis

  TH1D              *fMassPlot                   ; //!<!Mass spectra to be fitted
  TH2D              *fMassVsJetPtPlot            ; //!<!Mass vs jet pt (SB method)
  TH2D              *fMassVsJetzPlot             ; //!<!Mass vs jet z (SB method)

private:
  ClassDef(AliDJetVReader,3);
};

#endif
