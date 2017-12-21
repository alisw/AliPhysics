#ifndef ALIMTRPARAMETERIZEDRESPONSE_H
#define ALIMTRPARAMETERIZEDRESPONSE_H

/// \class AliMTRParameterizedResponse
/// \brief Class to handle the Lpt/Apt
///
/// The class provides tools to handle the trigger response function
/// for a better/easier inclusion in simulations
///
/// \author Diego Stocco <dstocco@cern.ch>, Subatech
/// \date June 3, 2016

#include "TObject.h"

class TH2;
class TObjArray;
class TAxis;
class TF1;

class AliMTRParameterizedResponse : public TObject {
 public:
  AliMTRParameterizedResponse();
  virtual ~AliMTRParameterizedResponse();

  AliMTRParameterizedResponse ( const AliMTRParameterizedResponse& );
  AliMTRParameterizedResponse& operator = ( const AliMTRParameterizedResponse& );

  enum {
    kAptOverAll,
    kLptOverApt,
    kHptOverLpt
  };

  Bool_t FitResponses ( Bool_t buildDataAptOverAllFromGraph = kTRUE, Bool_t fitLowPtIncrease = kFALSE );

  Double_t WeightPerBoard ( Double_t pt, Int_t iboard, Int_t itype, Bool_t isMC, Bool_t useFit ) const;
  Double_t WeightPerEta ( Double_t pt, Double_t eta, Int_t itype, Bool_t isMC, Bool_t useFit ) const;

  Bool_t SetAptOverAllMC ( TH2* histoApt, TH2* histoAll );
  Bool_t SetHptOverLpt ( TH2* histoHpt, TH2* histoLpt, TH2* histoMCHpt = 0x0, TH2* histoMCLpt = 0x0 );
  Bool_t SetLptOverApt ( TH2* histoLpt, TH2* histoApt, TH2* histoMCLpt = 0x0, TH2* histoMCApt = 0x0 );

  Bool_t SetFromMTRResponseTaskOutput ( Bool_t perBoard, const char* filenameData,
                                        const char* filenameMC = "",
                                        const char* filenameMCApt = "",
                                        Int_t rebin = 3 );
//                                       const char* outputNameData = "MTRResponseOut", const char* outputNameMC = "MTRResponseOut" const char* identifierData = "out", const char* identifierMC = "out");

  Bool_t ShowResponses ( Int_t itype, Bool_t isMC, Bool_t perBoard ) const;
  Bool_t CompareResponses ( Int_t itype, Bool_t perBoard ) const;

  static void ZoomPad();

 private:

  Double_t GetWeight ( Double_t pt, Int_t ibin, Int_t itype, Bool_t isMC, Bool_t useFit, Bool_t perBoard ) const;
  Bool_t SetRatio ( TH2* histoNum, TH2* histoDen, Int_t itype, Bool_t isMC, Int_t rebin = 1 );
  TObjArray* GetResponse ( Int_t itype, Bool_t isMC, Bool_t perBoard, Bool_t warn = kTRUE ) const;
  Bool_t AddResponse ( Int_t itype, Bool_t isMC, Bool_t perBoard, TObjArray* resp );

  TString GetStdName ( Int_t itype, Bool_t isMC, Bool_t perBoard ) const;

  Double_t FitFunctionErf ( Double_t* xVal, Double_t *par );
  Double_t FitFunctionErfApt ( Double_t* xVal, Double_t *par );

  void InitFunctionParams ( TF1* func, Bool_t fitLowPtIncrease, Int_t itype ) const;

//  TH1* GetHisto ( Int_t itype, Int_t isMC, Bool_t perBoard, Double_t val );

  TObjArray* fResponseBoard; /// Response per board
  TObjArray* fResponseEta; /// Response per eta

  TAxis* fEtaBinning; /// Eta binning
  TAxis* fPtBinning; /// Pt binning

  /// \cond CLASSIMP
  ClassDef(AliMTRParameterizedResponse, 1); // MTR trigger response
  /// \endcond
};

#endif
