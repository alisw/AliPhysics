#ifndef AliAnTOFtrack_H
#define AliAnTOFtrack_H

#include "AliUtilTOFParams.h"
#include "TMath.h"

using namespace AliUtilTOFParams;



///////////////////////////////////////////////////////////////////////////////
///                                                                          //
///               Class container for TOF analysis of Pi/K/p.                //
///                                                                          //
///                                                                          //
/// Authors:                                                                 //
/// N. Jacazio,  nicolo.jacazio[AROBASe]bo.infn.it                           //
///////////////////////////////////////////////////////////////////////////////

class TObject;

class AliAnTOFtrack : public TObject{
public:
  //Constructors and destructor
  AliAnTOFtrack();
  virtual ~AliAnTOFtrack();
  
  //Masks
  UShort_t        fTrkMask;                 //Mask for track information
  UChar_t         fTPCPIDMask;              //Mask for TPC PID information
  UShort_t        fTrkCutMask;              //Mask for Track cuts information
  UShort_t        fDCAXYIndex;              //Binned version of the XY impact parameters
  UShort_t        fDCAZIndex;               //Binned version of the Z impact parameters
  Double_t        fLength;                  //Track length
  Float_t         fLengthRatio;             //Ratio of length of the track matched and the one from a distant cluster
  Double_t        fTOFTime;                 //Time measured by TOF
  Float_t         fTOFMismatchTime;         //Time measured by TOF from a distant cluster
  Float_t         fTOFExpTime[kExpSpecies]; //Expected time for all species
  Float_t         fTOFExpSigma[kExpSpecies];//Expected sigmas for all species
  Float_t         fT0TrkTime;               //T0 best for the event
  Int_t           fTOFchan;                 //Channel of the matched track
  Float_t         fEta;                     //Eta distribution of the track
  Float_t         fPhi;                     //Phi distribution of the track
  Double_t        fPt;                      //Transverse momentum
  Double_t        fPTPC;                    //Momentum in the TPC
  Int_t           fNTOFClusters;            //Number of clusters matchable to the one matched
  Float_t         fTPCSignal;               //Signal in the TPC
  
  
  //******************************
  ////////Utility methods/////////
  //******************************
  
  
  //DCA binning
  ///
  /// Method to tests the DCA binning
  void TestDCAXYBinning();
  
  ///
  /// Method to put the DCAxy or DCAz of the class into bins
  void ComputeDCABin(const Double_t dca, const Bool_t xy);
  
  ///
  /// Method to put the DCAxy and DCAz of the class into bins
  void ComputeDCABin(const Double_t dcaxy, const Double_t dcaz);
  
  ///
  /// Method to convert compute the bin limits of the DCA binning
  void GetBinnedDCA(Double_t &down, Double_t &up, const Bool_t xy);
  
  ///
  /// Method to get the value of the track DCA
  Double_t GetDCA(const Bool_t xy);
  
  //TOF utilities
  
  ///
  /// Method to get the diffenece between the track time and the expected one
  Float_t GetDeltaT(const UInt_t id);
  
  ///
  /// Method to get the diffenece between the track time and the expected one in Number of sigmas
  Float_t GetDeltaSigma(const UInt_t id, const UInt_t hypo);
  
  ///
  /// Method to get the resolution on T0 based on the expected sigma of electrons
  Double_t GetT0Resolution(const Double_t TOFsigma = 80) const;
  
  //T0 Methods
  
  ///
  /// Method to get if the T0 time is T0 TOF
  Bool_t IsT0TOF(const Bool_t exclusive = kFALSE);
  
  ///
  /// Method to get if the T0 time is T0 T0A
  Bool_t IsT0A(const Bool_t exclusive = kFALSE);
  
  ///
  /// Method to get if the T0 time is T0 T0C
  Bool_t IsT0C(const Bool_t exclusive = kFALSE);
  
  ///
  /// Method to get if the T0 time is T0 Fill
  Bool_t IsT0Fill();
  
  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0A
  Bool_t IsT0TOF_T0A();
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A
  Bool_t IsOrT0TOF_T0A();
  
  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0C
  Bool_t IsT0TOF_T0C();
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0C
  Bool_t IsOrT0TOF_T0C();
  
  ///
  /// Method to get if the T0 time is T0 T0A and T0 T0C
  Bool_t IsT0A_T0C();
  
  ///
  /// Method to get if the T0 time is T0 T0A or T0 T0C
  Bool_t IsOrT0A_T0C();
  
  ///
  /// Method to get if the T0 time is T0 TOF, T0 T0A and T0 T0C
  Bool_t IsT0TOF_T0A_T0C();
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A or T0 T0C
  Bool_t IsOrT0TOF_T0A_T0C();
  
  // Cut flags 
  
  ///
  /// Method to check the standard cuts
  Bool_t PassStdCut();
  
  ///
  /// Method to check all the cut variations
  Bool_t PassCut(const Int_t cut = -1);
  
  // TPC PID
  
  ///
  /// TPC PID for electrons
  Bool_t IsTPCElectron();
  
  ///
  /// Method to check if it is TPC Pi K P
  Bool_t IsTPCPiKP(const UInt_t i);
  
  ///
  /// Method to check that in TPC the PID is for Pi K P
  Bool_t IsTPCPiKP();
  
  ///
  /// Method to check consistency between TOF and TPC for Pi K P to remove mismatch
  Bool_t ConsistentTPCTOF();
  
  ///
  /// Method to get the particle charge
  Bool_t IsNegative();
  
  ///
  /// Method to get the particle theta
  Double_t GetTheta() const { return 2.*TMath::ATan(TMath::Exp(-fEta)); }
  
  ///
  /// Method to get the particle momentum
  Double_t GetMomentum() const { return fPt/TMath::Sin(GetTheta()); }
  
  ClassDef(AliAnTOFtrack, 7);//AliAnTOFtrack : TOF analysis container class
};

#endif
