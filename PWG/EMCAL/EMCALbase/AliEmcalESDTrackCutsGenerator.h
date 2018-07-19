/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
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

class AliESDtrackCuts;
class TString;
class AliEmcalTrackSelection;

namespace PWG {

namespace EMCAL{


class AliEmcalESDTrackCutsGenerator {
public:
  enum EDataSet_t {
    kUnknown   = 0,
    kLHC10bcde = 1,
    kLHC10h    = 2,
    kLHC11a    = 3,
    kLHC11c    = 4,
    kLHC11d    = 5,
    kLHC11h    = 6,
    kLHC17o_TRD = 7
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

}

}

#endif /* ALIEMCALESDTRACKCUTSGENERATOR_H */
