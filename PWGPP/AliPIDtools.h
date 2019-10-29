#ifndef ALIPIDTOOLS_H
#define ALIPIDTOOLS_H

/// \ingroup PWGPP
/// \class AliPIDtools
/// \brief Wrapper for the AliPID classes - enable to use function in TFormula, TTreeFormula and Python

/// #### Example 1: register PID for run 246751 and draw dEdx as function of BetaGamma
/// \code
/// Int_t hash= AliPIDtools::LoadPID(246751,1,"pass1",0);
///  TF1 * fbb = new TF1("fbb",Form("AliPIDtools::BetheBlochAleph(%d,x)",hash),1,100);
///  fbb->SetNpx(1000); fbb->Draw();
/// \endcode
/// #### Example 2: Draw dEdx ratios pi/proton resp. pi/Kaon
/// \code
/// TF1 *fPionToProton1P=new TF1("fPionToProton",Form("AliPIDtools::GetExpectedTPCSignal(%d,1/x,2)/AliPIDtools::GetExpectedTPCSignal(%d,1/x*%f/%f,4)",hash,hash,massProton,massPi),0.,5);
///  TF1 *fPionToKaon1P=new TF1("fPionToKaon",Form("AliPIDtools::GetExpectedTPCSignal(%d,1/x,2)/AliPIDtools::GetExpectedTPCSignal(%d,1/x*%f/%f,3)",hash,hash,massKaon,massPi),0.,5);
///  fPionToKaon->SetLineColor(4);
///  fPionToProton1P->Draw(); fPionToKaon1P->Draw("same");
/// \endcode


#include "map"
#include  "AliESDtrack.h"
class AliPIDResponse;
class AliTPCPIDResponse;
class AliITSPIDResponse;
class AliTOFPIDResponse;



class AliPIDtools {
public:
  static Int_t GetHash(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC);
  static Int_t LoadPID(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC);
  static AliPIDResponse *GetPID(Int_t hash);
  static AliTPCPIDResponse &GetTPCPID(Int_t hash);
  static AliTOFPIDResponse &GetTOFPID(Int_t hash);
  static AliITSPIDResponse &GetITSPID(Int_t hash);
  //
  static Double_t BetheBlochAleph(Int_t hash, Double_t bg);
  static Double_t BetheBlochITS(Int_t hash, Double_t p, Double_t mass);
  static Double_t GetExpectedTPCSignal(Int_t hash, Double_t p, Int_t  particle);
  static Double_t GetExpectedITSSignal(Int_t hash, Double_t p, Int_t  particle);
  static Double_t GetExpectedTOFSigma(Int_t hash, Float_t mom, Int_t type);
  static Double_t GetExpectedTOFSignal(Int_t hash, const AliVTrack *track, Int_t  type);
  //
  static std::map<Int_t, AliTPCPIDResponse *> pidTPC;     /// we should use better hash map
  static std::map<Int_t, AliPIDResponse *> pidAll;        /// we should use better hash map
private:
  static AliESDtrack  dummyTrack;/// dummy value to save CPU - unfortunately PID object use AliVtrack - for the moment create global varaible t avoid object constructions
};

#endif //ALIPIDTOOLS_H
