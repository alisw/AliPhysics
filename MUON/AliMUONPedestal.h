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

#include <TObject.h>

#include "AliMpConstants.h"

class AliMUONVStore;

class TTimeStamp;

// global variables
const Int_t kNChannels = AliMpConstants::ManuNofChannels();
const Int_t kADCMax    = 4095;

class AliMUONPedestal : public TObject
{
  public:
    AliMUONPedestal();
    virtual ~AliMUONPedestal();
    
    /// return the number of events
    void SetAliNEvents(Int_t events) {fNEvents = events;}
    /// return the Run number
    void SetAliRunNumber(Int_t run) {fRunNumber = run;}
    /// return the total number of channels (or pads)
    void SetAliNChannel(Int_t nch) {fNChannel = nch;}
    /// output .log file of DAs
    void SetAlifilcout(ofstream* stream) {fFilcout = stream;}
    /// return date and time
    TTimeStamp* GetDate() {return fDate;}
    /// Count parity errors per Buspatch
    AliMUONVStore* GetErrorBuspatchTable() {return fErrorBuspatchTable;}
    /// return the name of DAPedestal .root file
    Char_t* GetHistoFileName() {return fHistoFileName;}
    void MakePed(Int_t bp,Int_t manu,Int_t ch,Int_t charge);

    /// set specific  DA prefixname
    void SetprefixDA(char* folder) {sprintf(fprefixDA,"%s",folder);}
    /// set the index of calibration runs
    void SetAliIndex(Int_t ind) {fIndex = ind;}
    /// Compute the pedestal data (mean, sigma)
    void Finalize();
    /// Create String to be put into file or AMORE DB
    void MakeASCIIoutput(ostream& out) const;
    /// Fill Histograms
    void MakeControlHistos();

  protected:
    Int_t fN; ///<
    Int_t fNEvents; ///< Number of events
    Int_t fRunNumber; ///< run number
    Int_t fNChannel; ///< Nb of channels (pads)
    Int_t fNManu; ///<  Nb of Manu
    AliMUONVStore* fErrorBuspatchTable; ///< Table for buspatches with parity errors 
    AliMUONVStore* fManuBuspatchTable; ///< Occupancy rate for each (buspatch, manu)
 
    TTimeStamp* fDate; ///< date
    ofstream* fFilcout; ///< .log output file
    Char_t fHistoFileName[256]; ///< .root histo file
    AliMUONVStore* fPedestalStore; ///< 
    Int_t fIndex; ///< calibration run index
    Char_t fprefixDA[256]; ///< specific DA prefixname

  private:
    /// Not implemented
    AliMUONPedestal(const AliMUONPedestal& rhs);
    /// Not implemented
    AliMUONPedestal& operator = (const AliMUONPedestal& rhs);

  ClassDef(AliMUONPedestal,2) // 
};

#endif
