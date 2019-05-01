#ifndef AliAnTOFtrack_H
#define AliAnTOFtrack_H

#include "AliUtilTOFParams.h"
#include "TMath.h"
#include "iostream"
using std::cout;
using std::endl;

using namespace AliUtilTOFParams;

///////////////////////////////////////////////////////////////////////////////
///                                                                          //
///               Class container for TOF analysis of Pi/K/p.                //
///                                                                          //
///                                                                          //
/// Authors:                                                                 //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                           //
///////////////////////////////////////////////////////////////////////////////

class AliAnTOFtrack {
  public:
  //Constructors and destructor
  AliAnTOFtrack();
  virtual ~AliAnTOFtrack();

  //***************************
  ////////Data Members/////////
  //***************************

  UShort_t fTrkMask;                 //Mask for track information
  UChar_t fTPCPIDMask;               //Mask for TPC PID information
  UShort_t fTrkCutMask;              //Mask for Track cuts information
  Double32_t fDCAXY;                 //[-3,3,12]  Binned version of the XY impact parameters
  Double32_t fDCAZ;                  //[-3,3,12]  Binned version of the Z impact parameters
  Float_t fLength;                   //Track length
  Float_t fLengthRatio;              //Ratio of length of the track matched and the one from a distant cluster
  Float_t fTOFTime;                  //Time measured by TOF
  Float_t fTOFMismatchTime;          //Time measured by TOF from a distant cluster
  Float_t fTOFExpTime[kExpSpecies];  //Expected time for all species
  Float_t fTOFExpSigma[kExpSpecies]; //Expected sigmas for all species
  Float_t fT0TrkTime;                //T0 best for the event
  Double32_t fT0TrkSigma;            //[0,1048.576,20]  T0 best resolution for the event
  Double32_t fTOFchan;               //[-1.5,262142.5,18]  Channel of the matched track
  Float_t fEta;                      //Eta distribution of the track
  Double32_t fPhi;                   //[0,6.5536,16]  Phi distribution of the track
  Float_t fPt;                       //Transverse momentum
  Double32_t fPTPC;                  //[0,40,19]  Momentum in the TPC
  Double32_t fNTOFClusters;          //[-1.5,1022.5,10]  Number of clusters matchable to the one matched
  Double32_t fTPCSignal;             //[0,1000,20]  Signal in the TPC

  //******************************
  ////////Utility methods/////////
  //******************************

  //DCA binning

  ///
  /// Method to put the DCAxy or DCAz of the class into bins
  void ComputeDCABin(const Double_t dca, const Bool_t xy);

  ///
  /// Method to put the DCAxy and DCAz of the class into bins
  void ComputeDCABin(const Double_t dcaxy, const Double_t dcaz);

  ///
  /// Method to get the value of the track DCA
  Double_t GetDCA(const Bool_t xy) const { return xy ? fDCAXY : fDCAZ; }

  //TOF utilities

  ///
  /// Method to get the diffenece between the track time and the expected one
  Float_t GetDeltaT(const UInt_t id) const;

  ///
  /// Method to get the diffenece between the track time and the expected one in Number of sigmas
  Float_t GetDeltaSigma(const UInt_t id, const UInt_t hypo) const;

  ///
  /// Method to get the resolution on T0 based on the expected sigma of electrons
  Double_t GetT0Resolution(const Double_t TOFsigma = 80) const;

  ///
  /// Method to get the track beta
  Double_t GetBeta() const { return fLength / ((fTOFTime - fT0TrkTime) * CSPEED); }

  ///
  /// Method to get the track gamma * beta
  Double_t GetGammaBeta() const;

  ///
  /// Method to get the track mass from the momentum at vertex
  Double_t GetTOFMass() const;

  //T0 Methods

  ///
  /// Method to get if the T0 time is T0 TOF
  Bool_t IsT0TOF(const Bool_t exclusive = kFALSE) const;

  ///
  /// Method to get if the T0 time is T0 T0A
  Bool_t IsT0A(const Bool_t exclusive = kFALSE) const;

  ///
  /// Method to get if the T0 time is T0 T0C
  Bool_t IsT0C(const Bool_t exclusive = kFALSE) const;

  ///
  /// Method to get if the T0 time is T0 Fill
  Bool_t IsT0Fill() const;

  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0A
  Bool_t IsT0TOF_T0A() const;

  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A
  Bool_t IsOrT0TOF_T0A() const;

  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0C
  Bool_t IsT0TOF_T0C() const;

  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0C
  Bool_t IsOrT0TOF_T0C() const;

  ///
  /// Method to get if the T0 time is T0 T0A and T0 T0C
  Bool_t IsT0A_T0C() const;

  ///
  /// Method to get if the T0 time is T0 T0A or T0 T0C
  Bool_t IsOrT0A_T0C() const;

  ///
  /// Method to get if the T0 time is T0 TOF, T0 T0A and T0 T0C
  Bool_t IsT0TOF_T0A_T0C() const;

  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A or T0 T0C
  Bool_t IsOrT0TOF_T0A_T0C() const;

  // Cut flags

  ///
  /// Method to check the standard cuts
  Bool_t PassStdCut() const;

  ///
  /// Method to check all the cut variations
  Bool_t PassCut(const Int_t cut = -1) const;

  // TPC PID

  ///
  /// TPC PID for electrons
  Bool_t IsTPCElectron() const;

  ///
  /// Method to check if it is TPC Pi K P
  Bool_t IsTPCPiKP(const UInt_t i) const;

  ///
  /// Method to check that in TPC the PID is for Pi K P
  Bool_t IsTPCPiKP() const;

  ///
  /// Method to check consistency between TOF and TPC for Pi K P to remove mismatch
  Bool_t ConsistentTPCTOF() const;

  ///
  /// Method to get the particle charge
  Bool_t IsNegative() const;

  ///
  /// Method to get the particle theta
  Double_t GetTheta() const { return 2. * TMath::ATan(TMath::Exp(-fEta)); }

  ///
  /// Method to get the particle momentum
  Double_t GetMomentum() const { return fPt / TMath::Sin(GetTheta()); }

  ///
  /// Method to get the particle momentum
  void Print() const;
};

#endif
