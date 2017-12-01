#ifndef AliFlatTPCdEdxInfo_H
#define AliFlatTPCdEdxInfo_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               *
 * Primary Authors : Sergey Gorbunov, Jochen Thaeder, Chiara Zampolli     */

/**
 * >> Flat structure representingTPCdEdxInfo <<
 */

#include "Rtypes.h"
#include "AliVMisc.h"
#include "AliTPCdEdxInfo.h"

class AliFlatTPCdEdxInfo
{
 public:

  // -- Constructor / Destructors
 
  AliFlatTPCdEdxInfo();
  ~AliFlatTPCdEdxInfo(){}
 
  // constructor and method for reinitialisation of virtual table
  AliFlatTPCdEdxInfo( AliVConstructorReinitialisationFlag );
  void Reinitialize() { new (this) AliFlatTPCdEdxInfo( AliVReinitialize ); }
 
  //--
 
  static size_t GetSize(){ return sizeof(AliFlatTPCdEdxInfo); }

  // -- Set / Get methods
 
  void SetFromTPCdEdxInfo( const AliTPCdEdxInfo *info );
  void GetTPCdEdxInfo( AliTPCdEdxInfo *info ) const;

  void SetSignalTot( Float_t signalTot[4] );
  void SetSignalMax( Float_t signalMax[4] );
  void SetNumberOfClusters( Char_t ncl[3] );
  void SetNumberOfCrossedRows( Char_t nrows[3] );

  void SetTPCsignal( Float_t signal, Float_t rms, UShort_t ncl );
  void GetTPCsignal( Float_t &signal, Float_t &rms, UShort_t &ncl ) const;

  // -- implementation of AliTPCdEdxInfo methods

  Double_t GetWeightedMean(Int_t qType, Int_t wType, Double_t w0, Double_t w1, Double_t w2) const;
  Double_t GetFractionOfClusters(Int_t iregion) const;
  
  Double_t GetSignalTot(Int_t index){ return fTPCsignalRegion[index];}
  Double_t GetSignalMax(Int_t index){ return fTPCsignalRegionQmax[index];}
  Double_t GetNumberOfClusters(Int_t index) {return fTPCsignalNRegion[index%3];}
  Double_t GetNumberOfCrossedRows(Int_t index) {return fTPCsignalNRowRegion[index%3];}
  //
  Double_t GetTPCsignalShortPad()      const {return fTPCsignalRegion[0];}
  Double_t GetTPCsignalMediumPad()     const {return fTPCsignalRegion[1];}
  Double_t GetTPCsignalLongPad()       const {return fTPCsignalRegion[2];}
  Double_t GetTPCsignalOROC()          const {return fTPCsignalRegion[3];}
  
  Double_t GetTPCsignalShortPadQmax()  const {return fTPCsignalRegionQmax[0];}
  Double_t GetTPCsignalMediumPadQmax() const {return fTPCsignalRegionQmax[1];}
  Double_t GetTPCsignalLongPadQmax()   const {return fTPCsignalRegionQmax[2];}
  Double_t GetTPCsignalOROCQmax()      const {return fTPCsignalRegionQmax[3];}

  
 private:

  // === content of AliTPCdEdxInfo class ==
  
  Float_t  fTPCsignalRegion[4]; //[0.,0.,10] TPC dEdx signal in 4 different regions - 0 - IROC, 1- OROC medium, 2 - OROC long, 3- OROC all, (default truncation used)  - for qTot
  Float_t  fTPCsignalRegionQmax[4]; //[0.,0.,10] TPC dEdx signal in 4 different regions - 0 - IROC, 1- OROC medium, 2 - OROC long, 3- OROC all, (default truncation used) - for qMax
  Char_t   fTPCsignalNRegion[3]; // number of clusters above threshold used in the dEdx calculation
  Char_t   fTPCsignalNRowRegion[3]; // number of crosed rows used in the dEdx calculation - signal below threshold included

  // === content of dEdX part of AliESDtrack ===

  Float_t  fTPCsignal;     // [0.,0.,10] detector's PID signal
  Float_t  fTPCsignalRMS;  // [0.,0.,10] RMS of dEdx measurement
  UShort_t fTPCsignalN;    // number of points used for dEdx
  
};

#pragma GCC diagnostic ignored "-Weffc++" 
inline AliFlatTPCdEdxInfo::AliFlatTPCdEdxInfo( AliVConstructorReinitialisationFlag ) 
{
  // no virtual table:  do nothing 
} 
#pragma GCC diagnostic warning "-Weffc++" 

inline AliFlatTPCdEdxInfo::AliFlatTPCdEdxInfo()
  :
  fTPCsignal(0.),
  fTPCsignalRMS(0.),
  fTPCsignalN(0.)			  
{
  for( int i=0; i<4; i++){
    fTPCsignalRegion[i] = 0.;
    fTPCsignalRegionQmax[i] = 0.;
  }
  for( int i=0; i<3; i++){
    fTPCsignalNRegion[i] = 0;
    fTPCsignalNRowRegion[i] = 0;
  }
}


inline void AliFlatTPCdEdxInfo::SetSignalTot( Float_t signalTot[4] )
{
  for( int i=0; i<4; i++)  fTPCsignalRegion[i] = signalTot[i];
}

inline void AliFlatTPCdEdxInfo::SetSignalMax( Float_t signalMax[4] )
{
  for( int i=0; i<4; i++)  fTPCsignalRegionQmax[i] = signalMax[i];
}

inline void AliFlatTPCdEdxInfo::SetNumberOfClusters( Char_t ncl[3] )
{
  for( int i=0; i<3; i++)  fTPCsignalNRegion[i] = ncl[i];
}

inline void AliFlatTPCdEdxInfo::SetNumberOfCrossedRows( Char_t nrows[3] )
{
  for( int i=0; i<3; i++) fTPCsignalNRowRegion[i] = nrows[i];
}

inline void AliFlatTPCdEdxInfo::SetTPCsignal( Float_t signal, Float_t rms, UShort_t ncl )
{ 
  fTPCsignal = signal;  
  fTPCsignalRMS = rms;
  fTPCsignalN = ncl;
}

inline void AliFlatTPCdEdxInfo::GetTPCsignal( Float_t &signal, Float_t &rms, UShort_t &ncl ) const
{ 
  signal = fTPCsignal;
  rms = fTPCsignalRMS;
  ncl = fTPCsignalN;
}

inline void AliFlatTPCdEdxInfo::SetFromTPCdEdxInfo( const AliTPCdEdxInfo *info )
{
  if( !info ) return;
  for( int i=0; i<4; i++){
    fTPCsignalRegion[i] = info->GetSignalTot(i);
    fTPCsignalRegionQmax[i] = info->GetSignalMax(i);
  }
  for( int i=0; i<3; i++){
    fTPCsignalNRegion[i] = info->GetNumberOfClusters(i);
    fTPCsignalNRowRegion[i] = info->GetNumberOfCrossedRows(i);
  } 
}

inline void AliFlatTPCdEdxInfo::GetTPCdEdxInfo( AliTPCdEdxInfo *info ) const
{
  if( !info ) return;
  Double_t signalTot[4], signalMax[4];
  for( int i=0; i<4; i++ ){ // convert to double to satisfy AliTPCdEdxInfo interface
    signalTot[i] = fTPCsignalRegion[i];
    signalMax[i] = fTPCsignalRegionQmax[i];
  }
  info->SetTPCSignalRegionInfoQmax(signalTot , fTPCsignalNRegion, fTPCsignalNRowRegion );
  info->SetTPCSignalsQmax(signalMax );
}

inline Double_t AliFlatTPCdEdxInfo::GetWeightedMean(Int_t qType, Int_t wType, Double_t w0, Double_t w1, Double_t w2) const
{
  //
  // Get weighted mean of the dEdx information
  //
  const Float_t *info = (qType==0)? fTPCsignalRegion :  fTPCsignalRegionQmax;
  const Char_t *ninfo = (wType==0)? fTPCsignalNRegion:  fTPCsignalNRowRegion;
  Double_t weight[3]={w0,w1,w2};
  Double_t sum=0;
  Double_t sumw=0;
  for (Int_t i=0; i<3; i++){
    sum+= info[i]*Double_t(ninfo[i])*weight[i];
    sumw+= ninfo[i]*weight[i];
  }
  Double_t result = (sumw>0) ? sum/sumw:0;
  return result;
}

inline  Double_t AliFlatTPCdEdxInfo::GetFractionOfClusters(Int_t iregion) const
{
  return fTPCsignalNRowRegion[iregion]>0 ? Double_t(fTPCsignalNRegion[iregion])/Double_t(fTPCsignalNRowRegion[iregion]):0.;
}

#endif
