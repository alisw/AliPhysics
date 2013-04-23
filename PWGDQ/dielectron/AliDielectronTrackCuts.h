#ifndef ALIDIELECTRONTRACKCUTS_H
#define ALIDIELECTRONTRACKCUTS_H

/* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//#############################################################
//#                                                           # 
//#         Class AliDielectronTrackCuts                     #
//#                                                           #
//#  Authors:                                                 #
//#   Anton     Andronic, GSI / A.Andronic@gsi.de             #
//#   Ionut C.  Arsene,   GSI / I.C.Arsene@gsi.de             #
//#   Julian    Book,     Uni Ffm / Julian.Book@cern.ch       #
//#   Frederick Kramer,   Uni Ffm, / Frederick.Kramer@cern.ch #
//#   Magnus    Mager,    CERN / Magnus.Mager@cern.ch         #
//#   WooJin J. Park,     GSI / W.J.Park@gsi.de               #
//#   Jens      Wiechula, Uni Tü / Jens.Wiechula@cern.ch      #
//#                                                           #
//#############################################################

#include <AliPID.h>
#include <AliAnalysisCuts.h>

class AliDielectronTrackCuts : public AliAnalysisCuts {
public:
  enum ITSClusterRequirement { kOff = 0, kNone, kAny, kFirst, kOnlyFirst, kSecond, kOnlySecond, kBoth };
  enum Detector { kSPD = 0, kSDD, kSSD };
  enum ITScluster { kSPD0=0x0001, kSPD1=0x0002,
                    kSDD0=0x0004, kSDD1=0x0008,
                    kSSD0=0x0010, kSSD1=0x0020};
  enum ITSclusterCutType { kOneOf=0, kAtLeast, kExact };
  enum EFilterBit  { kSwitchOff=0, kTPCqual=1, kTPCqualSPDany=4, kTPCqualSPDanyPIDele=8 };

  AliDielectronTrackCuts();
  AliDielectronTrackCuts(const char*name, const char* title);

  virtual ~AliDielectronTrackCuts();

  void SetV0DaughterCut(AliPID::EParticleType type, Bool_t negate=kFALSE);
  void SetClusterRequirementITS(Detector det, ITSClusterRequirement req = kOff) { fCutClusterRequirementITS[det] = req; }
  
  void SetRequireITSRefit(Bool_t req) { fRequireITSRefit=req; }
  void SetRequireTPCRefit(Bool_t req) { fRequireTPCRefit=req; }

  void SetTPCNclFRobust(Int_t cut) { fTPCNclRobustCut=cut; }
  void SetMinNCrossedRowsOverFindable(Double_t CrossedOverFindable) { fTPCcrossedOverFindable = CrossedOverFindable; }
  
  Int_t GetV0DaughterCut() const { return fV0DaughterCut; }
  ITSClusterRequirement GetClusterRequirementITS(Detector det) const { return fCutClusterRequirementITS[det]; }

  void SetITSclusterCut(ITSclusterCutType type, UChar_t map) { fITSclusterBitMap=map; fITSclusterCutType=type; }


  void SetAODFilterBit(EFilterBit type) { fAODFilterBit = type; }
  void SetMaxWaivedITSNcls(Int_t max) { fWaiveITSNcls = max; }

  //
  //Analysis cuts interface
  //
  virtual Bool_t IsSelected(TObject* track);
  virtual Bool_t IsSelected(TList*   /* list */ ) {return kFALSE;}
  

private:

  AliDielectronTrackCuts(const AliDielectronTrackCuts &c);
  AliDielectronTrackCuts &operator=(const AliDielectronTrackCuts &c);

  Int_t fV0DaughterCut;                                // Bit for V0 daughter selection
  Bool_t fNegateV0DauterCut;                           // If to negate the V0 daughter cut
  ITSClusterRequirement fCutClusterRequirementITS[3];  // detailed ITS cluster requirements for (SPD, SDD, SSD)

  UChar_t fITSclusterBitMap;                           // map of requested ITS clusters
  ITSclusterCutType fITSclusterCutType;                // logic of requested ITS clusters
  
  Bool_t fRequireITSRefit;                             // require ITS refit
  Bool_t fRequireTPCRefit;                             // require TPC refit

  Int_t fTPCNclRobustCut;                              // TPC Ncl cut, Robust, corresponds to 'crossed Rows' in ESDTrackCuts
  Double_t fTPCcrossedOverFindable;			           // TPC Crossed Rows / Findable Clusters Cut, analogous to ESDTrackCuts

  Int_t fAODFilterBit;                                 // Filter bit for AOD analysis
  Int_t fWaiveITSNcls;                                 // max number of waived ITS clusters after first hit 

  Bool_t CheckITSClusterRequirement(ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2) const;
  Bool_t CheckITSClusterCut(UChar_t itsBits) const;

  ClassDef(AliDielectronTrackCuts,3)         // Dielectron TrackCuts
};



#endif
