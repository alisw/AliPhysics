#ifndef ALIHFCUTVARFDSUBMINIMISER_H
#define ALIHFCUTVARFDSUBMINIMISER_H
/// \class AliHFCutVarFDsubMinimiser
/// \brief Minimiser for the cut variation feed down method analysis
///
///
///
///
/// \author Felix Reidt <felix.reidt@cern.ch>, CERN
/// \author Fabrizio Grosa <grosa@to.infn.it>, INFN Torino
/// \date Aug 17, 2015

#include "TObject.h"
#include "THnSparse.h"
#include "TH2F.h"

class TList;
class TH1F;

class AliHFCutVarFDsubMinimiser : public TObject {
protected:
  TH1F* fhRawYields; //!<!
  TH1F* fhEffPrompt; //!<!
  TH1F* fhEffFD; //!<!

  UInt_t fNiterations; ///
  Bool_t fUseWeights; ///
  Double_t fRelSystEffErr; ///

  Double_t fPromptYield; ///
  Double_t fPromptYieldErr; ///
  Double_t fFDYield; ///
  Double_t fFDYieldErr; ///

  TH1F* fhResiduals; ///
  TH1F* fhPulls; ///
  TH1F* fhfPrompt; ///
  TH1F* fhfPromptRaw; ///

  Int_t fnSim; ///
  TH2F* fhIncDistError; //!<!

  Bool_t fMinimised; ///

  Bool_t MinimiseDefault(); ///
  Bool_t MinimiseInCentre(); ///

  Bool_t InCentre(Double_t *effPrompt, Double_t *effFD, Double_t *rawYield, Double_t *Ncorr, Double_t* fprompt, Double_t* fpromptraw,Double_t* residuals); /// calculates the incenter

  AliHFCutVarFDsubMinimiser(const AliHFCutVarFDsubMinimiser& m); /// Copy constructor
  AliHFCutVarFDsubMinimiser operator=(const AliHFCutVarFDsubMinimiser& m); /// Assignment operator

public:
  AliHFCutVarFDsubMinimiser(); /// Default constructor
  AliHFCutVarFDsubMinimiser(TH1F* hRawYields, TH1F* hEffPrompt, TH1F* hEffFD,
                            UInt_t method = 0, UInt_t nIterations=10, Bool_t useWeights=kTRUE,
                            Double_t relSystEffErr=0., Int_t nSim=1000); /// Constructor
  ~AliHFCutVarFDsubMinimiser(); /// Destructor


  Bool_t GetStatus() { return fMinimised; } /// Get the minimisation status
  Double_t GetPromptYield()    { return fPromptYield;    }
  Double_t GetPromptYieldErr() { return fPromptYieldErr; }
  Double_t GetFDYield()        { return fFDYield;        }
  Double_t GetFDYieldErr()     { return fFDYieldErr;     }
  TH1F* GetResiduals()  { return fhResiduals;  }
  TH1F* GetPulls()      { return fhPulls;      }
  TH1F* GetFprompt()    { return fhfPrompt;    }
  TH1F* GetFpromptRaw() { return fhfPromptRaw; }
  TH2F* GetIncDistError() { return fhIncDistError; }

  /// \cond CLASSDEF
  ClassDef(AliHFCutVarFDsubMinimiser, 1);
  /// \endcond
};
#endif // ALIHFCUTVARFDSUBMINIMISER_H
