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
  enum EFilterBit  { kSwitchOff=0, kTPCqual=1, kITSonly=2, kTPCqualSPDany=4, kTPCqualSPDanyPIDele=8, kGlobalNoDCA=16, kTPConly=128 };

  AliDielectronTrackCuts();
  AliDielectronTrackCuts(const char*name, const char* title);

  virtual ~AliDielectronTrackCuts();

  // setters
  void SetV0DaughterCut(AliPID::EParticleType type, Bool_t negate=kFALSE);
  void SetClusterRequirementITS(Detector det, ITSClusterRequirement req = kOff) { fCutClusterRequirementITS[det] = req; }

  void SetRequireITSRefit(Bool_t req) { fRequireITSRefit=req; }
  void SetRequireTPCRefit(Bool_t req) { fRequireTPCRefit=req; }

  void SetTPCNclFRobust(Int_t cut) { fTPCNclRobustCut=cut; }
  void SetMinNCrossedRowsOverFindable(Double_t CrossedOverFindable) { fTPCcrossedOverFindable = CrossedOverFindable; }

  void SetITSclusterCut(ITSclusterCutType type, UChar_t map) { fITSclusterBitMap=map; fITSclusterCutType=type; }

  void SetGlobalTracksOnly(Bool_t setter = kTRUE) {fSelectGlobalTrack = setter;}
  void SetAODFilterBit(EFilterBit type) { fAODFilterBit = type; }
  void SetMaxWaivedITSNcls(Int_t max) { fWaiveITSNcls = max; }

  void SetRequireTRDUpdate(Bool_t req) { fRequireTRDUpdate=req; }

  void SetRequireCaloClusterMatch(Bool_t req, Short_t caloType) { fRequireCaloClusterMatch=req; fClusterMatchCaloType=caloType; }

  // getters
  Int_t GetV0DaughterCut() const { return fV0DaughterCut; }
  ITSClusterRequirement GetClusterRequirementITS(Detector det) const { return fCutClusterRequirementITS[det]; }

  Bool_t GetRequireITSRefit() const { return fRequireITSRefit; }
  Bool_t GetRequireTPCRefit() const { return fRequireTPCRefit; }

  Int_t GetTPCNclFRobust() const { return fTPCNclRobustCut; }
  Double_t GetMinNCrossedRowsOverFindable() const { return fTPCcrossedOverFindable; }

  UChar_t GetITSclusterCutMap() const { return fITSclusterBitMap; }
  ITSclusterCutType GetITSclusterCutType() const { return fITSclusterCutType; }

  Bool_t GetGlobalTracksOnly() const { return fSelectGlobalTrack; }
  Int_t GetAODFilterBit() const { return fAODFilterBit; }
  Int_t GetMaxWaivedITSNcls() const { return fWaiveITSNcls; }

  Bool_t GetRequireTRDUpdate() const { return fRequireTRDUpdate; }

  Bool_t  GetRequireCaloClusterMatch() const { return fRequireCaloClusterMatch; }
  Short_t GetCaloClusterMatchDetector() const { return fClusterMatchCaloType; }

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

  Bool_t fSelectGlobalTrack;                           // Cut for AODs to select global tracks (should in principle be covered by the filter bits, but if one wants to check efficiencies and cut only on acceptance, this selection is needed!)

  Bool_t fRequireITSRefit;                             // require ITS refit
  Bool_t fRequireTPCRefit;                             // require TPC refit

  Int_t fTPCNclRobustCut;                              // TPC Ncl cut, Robust, corresponds to 'crossed Rows' in ESDTrackCuts
  Double_t fTPCcrossedOverFindable;			           // TPC Crossed Rows / Findable Clusters Cut, analogous to ESDTrackCuts

  Int_t fAODFilterBit;                                 // Filter bit for AOD analysis
  Int_t fWaiveITSNcls;                                 // max number of waived ITS clusters after first hit

  Bool_t fRequireTRDUpdate;                            // require TRD update

  Bool_t  fRequireCaloClusterMatch;                    // require calo cluster matched to track
  Short_t fClusterMatchCaloType;                       // calo type for track match: AliDielectronClusterCuts::Detector

  Bool_t CheckITSClusterRequirement(ITSClusterRequirement req, Bool_t clusterL1, Bool_t clusterL2) const;
  Bool_t CheckITSClusterCut(UChar_t itsBits) const;

  ClassDef(AliDielectronTrackCuts,7)         // Dielectron TrackCuts
};



#endif
