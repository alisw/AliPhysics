#ifndef ALIMUONPEDESTAL_H
#define ALIMUONPEDESTAL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

// $Id$

/// \ingroup calib
/// \class AliMUONPedestal
/// \brief Implementation of the pedestal computing
/// 
//  Author Alberto Baldisseri, JL Charvet 

#ifndef ROOT_TObject
#  include <TObject.h>
#endif
#ifndef ROOT_TString
#  include <TString.h>
#endif

class AliMUONVStore;

class TTimeStamp;

class AliMUONPedestal : public TObject
{
  public:
    AliMUONPedestal();
    virtual ~AliMUONPedestal();
    
    /// return the number of current events
    void SetAliNCurrentEvents(Int_t events) {fNCurrentEvents = events;}
    /// return the number of events
    void SetAliNEvents(Int_t events) {fNEvents = events;}
    /// return the Run number
    void SetAliRunNumber(Int_t run) {fRunNumber = run;}
    /// return the total number of channels (or pads)
    void SetAliNChannel(Int_t nch) {fNChannel = nch;}
    /// output .log file of DAs
    void SetAlifilcout(ofstream* stream) {fFilcout = stream;}
    /// return date and time
    TTimeStamp* GetDate() const {return fDate;}
    /// Count parity errors per Buspatch
    AliMUONVStore* GetErrorBuspatchTable() const {return fErrorBuspatchTable;}
    /// return the name of DAPedestal .root file
    const char* GetHistoFileName() const;
    /// load MuonTrk configuration from ascii dbfile
    void LoadConfig(const char* dbfile);
    /// sum pedestal values for mean and sigma determination
    void MakePed(Int_t bp,Int_t manu,Int_t ch,Int_t charge);

    /// set config flag
    void SetconfigDA(Int_t ind) {fConfig = ind;}
    /// set Nb of evt threshold to calculate pedestal
    void SetnEvthreshold(Int_t ind) {fNEvthreshold = ind;}
    /// set specific  DA prefixname
    void SetprefixDA(const char* folder) { fPrefixDA=folder;}
    /// set the index of calibration runs
    void SetAliIndex(Int_t ind) {fIndex = ind;}
    /// Compute the pedestal data (mean, sigma)
    void Finalize();
    /// Create String to be put into file or AMORE DB
    void MakeASCIIoutput(ostream& out) const;
    /// Fill Histograms
    void MakeControlHistos();

  Int_t ADCMax() const { return 4095; }

protected:
    //    Int_t fN; ///<
    Int_t fNCurrentEvents; ///< Number of current events
    Int_t fNEvthreshold; ///< Nbevt threshold (pedestal calculation)
    Int_t fNEvents; ///< Number of events
    Int_t fRunNumber; ///< run number
    Int_t fNChannel; ///< Nb of channels (pads)
    Int_t fNManu; ///<  Nb of Manu
    Int_t fNManuConfig; ///<  Nb of Manu in the current detector configuration
    Int_t fConfig; ///< flag 1(0) for reading(or not) configuration ascii file
    AliMUONVStore* fErrorBuspatchTable; ///< Table for buspatches with parity errors 
    AliMUONVStore* fManuBuspatchTable; ///< Occupancy rate for each (buspatch, manu)
    AliMUONVStore* fManuBPoutofconfigTable; ///< (buspatch, manu) out of config
 
    TTimeStamp* fDate; ///< date
    ofstream* fFilcout; ///< .log output file
    TString fHistoFileName; ///< .root histo file
    AliMUONVStore* fPedestalStore; ///< data container:  (Pedmean,sigma) values for each (BP,manuId)
    Int_t fIndex; ///< calibration run index
    TString fPrefixDA; ///< specific DA prefixname

  static const Int_t fgkADCMax; ///< max channel count
  
  private:
    /// Not implemented
    AliMUONPedestal(const AliMUONPedestal& rhs);
    /// Not implemented
    AliMUONPedestal& operator = (const AliMUONPedestal& rhs);

  ClassDef(AliMUONPedestal,5) // Pedestal computing for DA 
};

#endif
