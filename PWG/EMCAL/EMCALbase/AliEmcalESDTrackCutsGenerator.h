/**
 * \file AliEmcalESDTrackCutsGenerator.h
 * \brief Declaration of class AliEmcalESDTrackCutsGenerator
 *
 * In this header file the class AliEmcalESDTrackCutsGenerator, which handles the
 * generation of ESD track cuts, is declared.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Jan 21, 2016
 */
#ifndef ALIEMCALESDTRACKCUTSGENERATOR_H
#define ALIEMCALESDTRACKCUTSGENERATOR_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class AliESDtrackCuts;
class AliEmcalTrackSelection;
class TString;

class AliEmcalESDTrackCutsGenerator {
public:
  enum EDataSet_t {
    kUnknown   = 0,
    kLHC10bcde = 1,
    kLHC10h    = 2,
    kLHC11a    = 3,
    kLHC11c    = 4,
    kLHC11d    = 5,
    kLHC11h    = 6
  };

  enum EStdCutMode_t {
    kRAA2011 = 1000,
    kGlobalTracksNCls90NoSPD = 1001,
    kGlobalTracksNCls80NoSPD = 1002,
    kGlobalTracks2010NCrossRows120 = 1003,
    kGlobalTracksNCls70NoSPD = 1004,
    kGlobalTracksNCls70NoSPDNoPtCut = 1005,
    kGlobalTracksNClsPtDepNoSPDNoPtCut = 1006,
    kGlobalTracks2011 = 1007,
    kGlobalTracks2011NoSPD = 1008,
    kGlobalTracksNCls90NoITS = 2000,
    kTPCOnlyTracksNCls70 = 2001,
    kTPCOnlyTracksNCrossRows120 = 2002
  };

  enum EAddCutMode_t {
    kSPDAny = 1000,
    kSPDNone = 1001,
    kNoITSChi2 = 1002,
    kNoMinTPCCls = 1003,
    kNoITSRefit = 1004,
    kSPDOff = 1005
  };

  static const Int_t fgkAddCutFactor;

  static AliESDtrackCuts* CreateTrackCutsPWGJE(Int_t cutMode);
  static AliESDtrackCuts* CreateTrackCutsPWGJE(Int_t stdCutMode, Int_t addCutMode);
  static AliESDtrackCuts* CreateTrackCutsPWGJE(Int_t stdCutMode, Int_t addCutMode1, Int_t addCutMode2);
  static TString SetStandardCuts(AliESDtrackCuts*& trackCuts, Int_t stdCutMode);
  static TString SetAdditionalCuts(AliESDtrackCuts*& trackCuts, Int_t addCutMode);
  static EDataSet_t SteerDataSetFromString(TString period);
  static void AddHybridTrackCuts(AliEmcalTrackSelection* trkSel, TString period) { AddHybridTrackCuts(trkSel, SteerDataSetFromString(period)); }
  static void AddHybridTrackCuts(AliEmcalTrackSelection* trkSel, EDataSet_t period);
  static void AddTPCOnlyTrackCuts(AliEmcalTrackSelection* trkSel, TString period) { AddTPCOnlyTrackCuts(trkSel, SteerDataSetFromString(period)); }
  static void AddTPCOnlyTrackCuts(AliEmcalTrackSelection* trkSel, EDataSet_t period);
};

#endif /* ALIEMCALESDTRACKCUTSGENERATOR_H */
