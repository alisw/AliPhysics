#ifndef AliAnTOFtrack_H
#define AliAnTOFtrack_H

#include "AliUtilTOFParams.h"
#include "TMath.h"

using namespace AliUtilTOFParams;



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//               Class container for TOF analysis of Pi/K/p.                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

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
  void TestDCAXYBinning(){
    Int_t dim = static_cast<Int_t>(TMath::Power(2., 8.*sizeof(fDCAXYIndex)));
    Double_t range = (2.*fDCAXYRange);
    std::cout<<"fDCAXYIndex has "<<dim<<" possibilities"<<std::endl;
    std::cout<<"DCAXY range is ["<<-fDCAXYRange<<","<<fDCAXYRange<<"] -> "<<range<<" bin width = "<<range/dim<<std::endl;
    std::cout<<"Variable dimension is "<<dim<<" while actual bin required are "<<kDCAXYBins<<std::endl;
    
    dim = static_cast<Int_t>(TMath::Power(2., 8.*sizeof(fDCAZIndex)));
    range = (2.*fDCAZRange);
    std::cout<<"fDCAZIndex has "<<dim<<" possibilities"<<std::endl;
    std::cout<<"DCAZ range is ["<<-fDCAZRange<<","<<fDCAZRange<<"] -> "<<range<<" bin width = "<<range/dim<<std::endl;
    std::cout<<"Variable dimension is "<<dim<<" while actual bin required are "<<kDCAZBins<<std::endl;
    
    for (Int_t var = 0; var <=20; ++var) {
      //       Double_t rnd = 2.*fDCAXYRange*gRandom->Rndm()-1;
      Double_t rnd = 2.*fDCAXYRange*(Double_t)var/20.-fDCAXYRange;
      ComputeDCABin(rnd, kTRUE);
      
      Double_t binlow, binup;
      GetBinnedDCA(binlow, binup, 1);
      
      std::cout<<rnd<<" Index: "<<fDCAXYIndex<<"   ["<<binlow<<" ; "<<binup<<"] --> Diff --> ["<< binlow - rnd <<" ; "<<binup - rnd<<"]"<<std::endl;
      //       cout<<rnd<<"  "<<TMath::Abs((rnd+fDCAXYRange)*kDCAXYBins/(2*fDCAXYRange))<<" Index: "<<fDCAXYIndex<<"   ["<<binlow<<" ; "<<binup<<"] -->  ["<< binlow - rnd <<" ; "<<binup - rnd<<"]"<<endl;
    }
    
    //     for (Int_t c = 0 ; c <= kDCAXYBins; ++c) cout<<c<<"  "<<-fDCAXYRange+2.*fDCAXYRange*c/kDCAXYBins<<endl;
    
  }
  
  ///
  /// Method to put the DCAxy or DCAz of the class into bins
  void ComputeDCABin(const Double_t dca, const Bool_t xy){
    if(xy){
      if(dca > fDCAXYRange){//Overflow
        fDCAXYIndex = kDCAXYBins+2;
        return;
      }
      else if(dca <= -fDCAXYRange){//Underflow
        fDCAXYIndex = kDCAXYBins+1;
        return;
      }
      fDCAXYIndex = static_cast<UShort_t>(TMath::Ceil((dca+fDCAXYRange)/((2*fDCAXYRange)/kDCAXYBins)) -1);
    }
    else{
      if(dca > fDCAZRange){//Overflow
        fDCAZIndex = kDCAZBins+2;
        return;
      }
      else if(dca <= -fDCAZRange){//Underflow
        fDCAZIndex = kDCAZBins+1;
        return;
      }
      fDCAZIndex = static_cast<UShort_t>(TMath::Ceil((dca+fDCAZRange)/((2*fDCAZRange)/kDCAZBins)) -1);
    }
    
  };
  
  ///
  /// Method to put the DCAxy and DCAz of the class into bins
  void ComputeDCABin(const Double_t dcaxy, const Double_t dcaz){
    
    ComputeDCABin(dcaxy, (Bool_t) kTRUE);
    ComputeDCABin(dcaz, (Bool_t) kFALSE);
  };
  
  ///
  /// Method to convert compute the bin limits of the DCA binning
  void GetBinnedDCA(Double_t &down, Double_t &up, const Bool_t xy){
    if(xy){
      if(fDCAXYIndex == kDCAXYBins+2){//Overflow
        down = fDCAXYRange;
        up = -999;
        return;
      }
      else if(fDCAXYIndex == kDCAXYBins+1){//Underflow
        down = -999;
        up = -fDCAXYRange;
        return;
      }
      
      down = (Double_t) (2.*fDCAXYRange)*fDCAXYIndex/kDCAXYBins - fDCAXYRange;
      up = (Double_t) (2.*fDCAXYRange)*(fDCAXYIndex+1)/kDCAXYBins - fDCAXYRange;
    }
    else{
      if(fDCAZIndex == kDCAZBins+2){//Overflow
        down = fDCAZRange;
        up = -999;
        return;
      }
      else if(fDCAZIndex == kDCAZBins+1){//Underflow
        down = -999;
        up = -fDCAZRange;
        return;
      }
      
      down = (Double_t) (2.*fDCAZRange)*fDCAZIndex/kDCAZBins - fDCAZRange;
      up = (Double_t) (2.*fDCAZRange)*(fDCAZIndex+1)/kDCAZBins - fDCAZRange;
    }
  }
  
  ///
  /// Method to get the value of the track DCA
  Double_t GetDCA(const Bool_t xy){
    Double_t u, d;
    GetBinnedDCA(d, u, xy);
    return d + (u-d)/2;
  }
  
  //TOF utilities
  
  ///
  /// Method to get the diffenece between the track time and the expected one
  Float_t GetDeltaT(const UInt_t id);
  
  ///
  /// Method to get the diffenece between the track time and the expected one in Number of sigmas
  Float_t GetDeltaSigma(const UInt_t id, const UInt_t hypo);
  
  ///
  /// Method to get the resolution on T0 based on the expected sigma of electrons
  Double_t GetT0Resolution(const Double_t TOFsigma = 80) const {
    return TMath::Sqrt(TMath::Power(fTOFExpSigma[0], 2) - TMath::Power(TOFsigma, 2) - TMath::Power(15./GetMomentum(), 2));
  }
  
  //T0 Methods
  
  ///
  /// Method to get if the T0 time is T0 TOF
  Bool_t IsT0TOF(){
    if(GetMaskBit(fTrkMask, kT0_0)) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 T0A
  Bool_t IsT0A(){
    if(GetMaskBit(fTrkMask, kT0_1)) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 T0C
  Bool_t IsT0C(){
    if(GetMaskBit(fTrkMask, kT0_2)) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 Fill
  Bool_t IsT0Fill(){
    if(!GetMaskBit(fTrkMask, kT0_0) && !GetMaskBit(fTrkMask, kT0_1) && !GetMaskBit(fTrkMask, kT0_2)) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0A
  Bool_t IsT0TOF_T0A(){
    if(IsT0TOF() && IsT0A()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A
  Bool_t IsOrT0TOF_T0A(){
    if(IsT0TOF() || IsT0A()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF and T0 T0C
  Bool_t IsT0TOF_T0C(){
    if(IsT0TOF() && IsT0C()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0C
  Bool_t IsOrT0TOF_T0C(){
    if(IsT0TOF() || IsT0C()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 T0A and T0 T0C
  Bool_t IsT0A_T0C(){
    if(IsT0A() && IsT0C()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 T0A or T0 T0C
  Bool_t IsOrT0A_T0C(){
    if(IsT0A() || IsT0C()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF, T0 T0A and T0 T0C
  Bool_t IsT0TOF_T0A_T0C(){
    if(IsT0TOF() && IsT0A() && IsT0C()) return kTRUE;
    return kFALSE;
  }
  
  ///
  /// Method to get if the T0 time is T0 TOF or T0 T0A or T0 T0C
  Bool_t IsOrT0TOF_T0A_T0C(){
    if(IsT0TOF() || IsT0A() || IsT0C()) return kTRUE;
    return kFALSE;
  }
  
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
  
  ClassDef(AliAnTOFtrack, 6);
};

#endif
