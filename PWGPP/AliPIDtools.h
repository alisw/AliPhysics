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
/// #### Example 3: Draw Expected dEdx
/// AliPIDtools::SetFilteredTreeV0(treeV0)
/// treeV0->Draw("log(track0.fTPCsignal/(AliPIDtools::GetExpectedTPCSignalV0(pidHash,0,0x1,0)))","type==1&&abs(log(track1.fTPCsignal/(AliPIDtools::GetExpectedTPCSignalV0(pidHash,0,0x1,1))))<0.1","colz",20000)

#include "map"
#include  "AliESDtrack.h"
class AliPIDResponse;
class AliTPCPIDResponse;
class AliITSPIDResponse;
class AliTOFPIDResponse;



class AliPIDtools {
  enum  { kEtaCorr=0x1, kMultCorr=0x2, kPileUpCorr=0x4 };
public:
  static Int_t GetHash(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC);
  static Int_t LoadPID(Int_t run, Int_t passNumber, TString recoPass, Bool_t isMC);
  static AliPIDResponse *GetPID(Int_t hash);
  static AliTPCPIDResponse &GetTPCPID(Int_t hash);
  static AliTOFPIDResponse &GetTOFPID(Int_t hash);
  static AliITSPIDResponse &GetITSPID(Int_t hash);
  //
  static Double_t BetheBlochAleph(Int_t hash, Double_t bg);
  static Double_t BetheBlochAleph(Int_t hash, Double_t p, Int_t type);
  static Double_t BetheBlochITS(Int_t hash, Double_t p, Double_t mass);
  static Double_t GetExpectedTPCSignal(Int_t hash, Double_t p, Int_t  particle);
  static Double_t GetExpectedITSSignal(Int_t hash, Double_t p, Int_t  particle);
  static Double_t GetExpectedTOFSigma(Int_t hash, Float_t mom, Int_t type);
  static Double_t GetExpectedTOFSignal(Int_t hash, const AliVTrack *track, Int_t  type);
  // TTree interface
  static AliESDtrack* GetCurrentTrack();
  static AliESDtrack* GetCurrentTrackV0(Int_t index);
  static TVectorD*    GetTOFInfo(Int_t infoType);
  static TVectorD*    GetTOFInfoV0(Int_t source, Int_t infoType);
  static Float_t      GetTOFInfoAt(Int_t infoType, Int_t index) {return (*GetTOFInfo(infoType))[index];}
  static Float_t      GetTOFInfoV0At(Int_t source, Int_t infoType,Int_t index) {return (*GetTOFInfoV0(source, infoType))[index];}
  static Bool_t       SetTPCEventInfo(Int_t hash,  Int_t corrMask);
  static Bool_t       SetTPCEventInfoV0(Int_t hash, Int_t corrMask);
  //
  static Double_t GetExpectedTPCSignal(Int_t hash, Int_t particleType, Int_t corrMask, Int_t returnType);
  static Double_t GetExpectedTPCSignalV0(Int_t hash, Int_t particleType, Int_t corrMask, Int_t index, Int_t returnType);
  static Double_t GetITSPID(Int_t hash, Int_t particleType, Int_t valueType, Float_t resol=0);
  //
  static Float_t NumberOfSigmas(Int_t hash, Int_t detCode, Int_t particleType, Int_t source=-1, Int_t corrMask=-1);
  static Float_t GetSignalDelta(Int_t hash, Int_t detCode, Int_t particleType, Int_t source=-1, Int_t corrMask=-1);
  static Float_t ComputePIDProbability(Int_t hash, Int_t detCode, Int_t particleType, Int_t source=-1, Int_t corrMask=-1,Int_t norm=1, Float_t fakeProb=0.01, Float_t* pidVector=0);
  static Float_t ComputePIDProbabilityCombined(Int_t hash, Int_t detMask, Int_t particleType, Int_t source=-1, Int_t corrMask=-1,Int_t norm=1, Float_t fakeProb=0.01);
  //
  static Bool_t    RegisterPIDAliases(Int_t pidHash, TString fakeRate="0.1", Int_t suffix=-1);
  //
  //
  static std::map<Int_t, AliTPCPIDResponse *> pidTPC;     /// we should use better hash map
  static std::map<Int_t, AliPIDResponse *> pidAll;        /// we should use better hash map
  //  TTree interface for interactive queries
  static Bool_t SetFilteredTree(TTree * filteredTree); /// set variable address in filtered trees
  static Bool_t SetFilteredTreeV0(TTree * filteredTreeV0); /// set variable address in filtered trees
  //static Bool_t SetPileUpProperties(const TVectorF * tpcVertex, const TVectorF *itsMult, AliTPCPIDResponse *pidTPC);
  static Bool_t SetPileUpProperties(const TVectorF & tpcVertexInfo, const TVectorF &itsClustersPerLayer, Int_t primMult, AliTPCPIDResponse *pidTPC);
  //
  static TTree *       fFilteredTree;  /// pointer to filteredTree
  static TTree *       fFilteredTreeV0;  /// pointer to filteredTree V0
  static void UnitTest();                       /// unit test of invariants
private:
  static AliESDtrack  dummyTrack;     /// dummy value to save CPU - unfortunately PID object use AliVtrack - for the moment create global varaible t avoid object constructions

};

#endif //ALIPIDTOOLS_H
