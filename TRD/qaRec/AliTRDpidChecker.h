#ifndef ALITRDPIDCHECKER_H
#define ALITRDPIDCHECKER_H

//////////////////////////////////////////////////////
//
// Task to check PID performance of the TRD
//
// Author : Alex Wilk <wilka@uni-muenster.de>
//
///////////////////////////////////////////////////////

#include "AliPID.h"
#include "../Cal/AliTRDCalPID.h"

#ifndef ALITRDRECOTASK_H
#include "AliTRDrecoTask.h"
#endif

class TObjArray;
class TList;
class TClonesArray;
class TTreeSRedirector;
class AliTRDReconstructor;
class AliTRDpidChecker : public AliTRDrecoTask 
{

  enum{
    kLQlikelihood    = 0                                           // place for 2-dim LQ electron likelihood distributions
    ,kNNlikelihood = 1 * AliTRDCalPID::kNMom * AliPID::kSPECIES  // place for NN electron likelihood distributions
    ,kdEdx         = 2 * AliTRDCalPID::kNMom * AliPID::kSPECIES  // place for the dE/dx spectra
    ,kPH           = 3 * AliTRDCalPID::kNMom * AliPID::kSPECIES  // place for pulse height spectra
    ,kMomentum     = 4 * AliTRDCalPID::kNMom * AliPID::kSPECIES  // place for the momentum distribution
    ,kMomentumBin  = kMomentum +1                                // place for the momentum distribution
    ,kGraphLQ      = kMomentumBin +1                             // place for the 2-dim LQ pion efficiencies
    ,kGraphLQerr   = kGraphLQ +1                                 // place for the 2-dim LQ pion efficiency errors
    ,kGraphNN      = kGraphLQerr +1                              // place for the NN pion efficiencies
    ,kGraphNNerr   = kGraphNN +1                                 // place for the NN pion efficiency errors
  };

  enum{
    kGraphStart = kGraphLQ
  };

public:
  AliTRDpidChecker();
  virtual ~AliTRDpidChecker();
  
  void    CreateOutputObjects();
  void    Exec(Option_t *option);
  void    GetRefFigure(Int_t ifig, Int_t &first, Int_t &last, Option_t *opt);  
  Bool_t  PostProcess();
  void    Terminate(Option_t *);


private:
  AliTRDpidChecker(const AliTRDpidChecker&);               // not implemented
  AliTRDpidChecker& operator=(const AliTRDpidChecker&);    // not implemented

  Double_t GetPionEfficiency(Int_t Index1, Int_t Index2);  // calculates the pion efficiency
  Double_t GetError(Int_t Index1, Int_t Index2);           // calculates the error
  

  AliTRDReconstructor *fReconstructor;     //! reconstructor needed for recalculation the PID

  enum{
    kBins = 12001                // binning of the likelihood histograms
  };

  ClassDef(AliTRDpidChecker, 1); // TRD PID checker
};

#endif
