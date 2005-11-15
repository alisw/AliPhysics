#ifndef ALIANALYSISTRACKCUTS_H
#define ALIANALYSISTRACKCUTS_H
/*  See cxx source for full Copyright notice */


/* $Id$ */

//-------------------------------------------------------------------------
//                       Class AliAnalysisTrackCuts
//   This is the class for the cuts in event & track level
//
//    Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

#include <TObject.h>

class AliESD;
class AliESDtrack;

class TPaveText;

class AliAnalysisTrackCuts : public TObject
{
 public:
  AliAnalysisTrackCuts();
  
  ~AliAnalysisTrackCuts();

  void Reset();
  
  void SetPRange(Float_t r1, Float_t r2);
  void SetPtRange(Float_t r1, Float_t r2);
  void SetPxRange(Float_t r1, Float_t r2);
  void SetPyRange(Float_t r1, Float_t r2);
  void SetPzRange(Float_t r1, Float_t r2);
  void SetBrRange(Float_t r1, Float_t r2);
  void SetBzRange(Float_t r1, Float_t r2);
  void SetEtaRange(Float_t r1, Float_t r2);
  void SetRapRange(Float_t r1, Float_t r2);
  
  Bool_t IsAccepted(AliESD *esd,AliESDtrack *esdtrack);

  TPaveText *GetTrackCuts();
  void PrintTrackCuts();
  void GetTrackStats();
  void GetPStats();
  void GetPxStats();
  void GetPyStats();
  void GetPzStats();
  void GetPtStats();
  void GetEtaStats();
  void GetRapStats();
  void GetBrStats();
  void GetBzStats();

 private:
  Float_t fPMin, fPMax;  //Definition of the range of the P
  Float_t fPtMin, fPtMax;  //Definition of the range of the Pt
  Float_t fPxMin, fPxMax;  //Definition of the range of the Px
  Float_t fPyMin, fPyMax;  //Definition of the range of the Py
  Float_t fPzMin, fPzMax;  //Definition of the range of the Pz
  Float_t fEtaMin, fEtaMax;  //Definition of the range of the eta
  Float_t fRapMin, fRapMax;  //Definition of the range of the y
  Float_t fBrMin, fBrMax;  //Definition of the range of the br
  Float_t fBzMin, fBzMax;  //Definition of the range of the bz

  Int_t fP;  //Number of events rejected due to P cut
  Int_t fPt;  //Number of events rejected due to Pt cut
  Int_t fPx;  //Number of events rejected due to Px cut
  Int_t fPy;  //Number of events rejected due to Py cut
  Int_t fPz;  //Number of events rejected due to Pz cut
  Int_t fEta;  //Number of events rejected due to eta cut
  Int_t fRap;  //Number of events rejected due to y cut
  Int_t fbr;  //Number of events rejected due to br cut
  Int_t fbz;  //Number of events rejected due to bz cut
  Int_t fTotalTracks;  //Total number of tracks
  Int_t fAcceptedTracks;  //Total number of accepted tracks
   
  Int_t fFlagP;  //Flag that shows if the P cut was imposed
  Int_t fFlagPt;  //Flag that shows if the Pt cut was imposed
  Int_t fFlagPx;  //Flag that shows if the Px cut was imposed
  Int_t fFlagPy;  //Flag that shows if the Py cut was imposed
  Int_t fFlagPz;  //Flag that shows if the Pz cut was imposed
  Int_t fFlagEta;  //Flag that shows if the eta cut was imposed
  Int_t fFlagRap;  //Flag that shows if the y cut was imposed
  Int_t fFlagbr;  //Flag that shows if the br cut was imposed
  Int_t fFlagbz;  //Flag that shows if the bz cut was imposed
 
  
  ClassDef(AliAnalysisTrackCuts, 1)
} ;

#endif
