#ifndef ALIMUONGAIN_H
#define ALIMUONGAIN_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONGain
/// \brief Implementation of the pedestal and gain computing
/// 
//  Author: Alberto Baldisseri, JL Charvet (05/05/2009)

#include "AliMUONPedestal.h"

class AliMUONGain : public AliMUONPedestal
{
  public:
    AliMUONGain();
    virtual ~AliMUONGain();

    TString WriteGainData(Int_t bp, Int_t manu, Int_t ch, Double_t p1, Double_t p2, Int_t threshold, Int_t q);
    TString WriteGainHeader(Int_t nInit, Int_t nEntries, Int_t nbpf2, Int_t *numrun, Double_t *injCharge);
    /// Store Pedmean and sigma in pedestal-like ascii file
    void MakePedStoreForGain(TString flatfile);
    /// Computes gain parameters and store in an ascii file
    void MakeGainStore(TString flatfile); 
    /// Set filename of root file containing pedmean and sigma values for each DAC value
  void SetAliRootDataFileName(const char* name="MUONTRKGAINda_data.root") { fRootDataFileName=name; }
    /// Get filename of root file containing pedmean and sigma values
    const char* GetRootDataFileName() const {return fRootDataFileName.Data();}
    /// Write Header in Dummy ascii file 
    TString WriteDummyHeader();
    /// Set InjCharge = DAC value
    void SetAliInjCharge(Int_t charge) {fInjCharge = charge;}
    /// Set PrintLevel
    void SetAliPrintLevel(Int_t pri) {fPrintLevel = pri;}
    /// Set fnInit (=1: first DAC=0 removed to compute the fit)
    void SetAliInit(Int_t ini) {fnInit = ini;}
    /// Set nEntries = Nb of DAC values
    void SetAliEntries(Int_t ent) {fnEntries = ent;}
    /// Set Nbpf1 = nb of DAC values for the linear fit
    void SetAliNbpf1(Int_t nf1) {fnbpf1 = nf1;}
    /// Set PlotLevel
    void SetAliPlotLevel(Int_t plo) {fPlotLevel = plo;}

  private:
    Int_t fInjCharge; ///< DAC value
    TString fRootDataFileName; ///< Root data file name
    Int_t fnInit; ///< fnInit (for expert)
    Int_t fnEntries; ///< Nb of DAC values
    Int_t fnbpf1; ///< nb  of DAC values for linear fit (for expert)
    Int_t fPrintLevel; ///< Print level
    Int_t fPlotLevel; ///< Plot level
	
  ClassDef(AliMUONGain,2) // 
};

#endif
