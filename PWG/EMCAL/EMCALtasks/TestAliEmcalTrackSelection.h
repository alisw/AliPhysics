// clang-format off
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
// clang-format on
#ifndef TESTALIEMCALTRACKSELECTION_H
#define TESTALIEMCALTRACKSELECTION_H

#include "AliAnalysisTaskEmcalLight.h"
#include <TString.h>

class AliAODEvent;
class AliAODTrack;
class AliEmcalTrackSelection;
class THistManager;
class TObjArray;

namespace PWG {

namespace EMCAL {

class TestImplAliEmcalTrackSelection : public TNamed {
 public:
  TestImplAliEmcalTrackSelection();
  TestImplAliEmcalTrackSelection(const char* name);
  virtual ~TestImplAliEmcalTrackSelection();

  bool RunTest(const AliAODEvent* const ev);

 protected:
  virtual bool IsTrueTrack(const AliAODTrack* const trk) const = 0;

  AliEmcalTrackSelection* fTrackSelection; ///< Object to be tested

 private:
  TestImplAliEmcalTrackSelection(const TestImplAliEmcalTrackSelection&);
  TestImplAliEmcalTrackSelection* operator=(const TestImplAliEmcalTrackSelection&);

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalTrackSelection, 1);
  /// \endcond
};

class TestImplAliEmcalTrackSelectionITSpure : public TestImplAliEmcalTrackSelection {
 public:
  TestImplAliEmcalTrackSelectionITSpure();
  TestImplAliEmcalTrackSelectionITSpure(const char* name, const char* period);
  virtual ~TestImplAliEmcalTrackSelectionITSpure();

 protected:
  virtual bool IsTrueTrack(const AliAODTrack* const trk) const;

 private:
  AliESDtrackCuts* fRefCuts; ///< Reference cuts

  TestImplAliEmcalTrackSelectionITSpure(const TestImplAliEmcalTrackSelectionITSpure&);
  TestImplAliEmcalTrackSelectionITSpure& operator=(const TestImplAliEmcalTrackSelectionITSpure&);

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalTrackSelectionITSpure, 1);
  /// \endcond
};

class TestImplAliEmcalTrackSelectionHybrid : public TestImplAliEmcalTrackSelection {
 public:
  TestImplAliEmcalTrackSelectionHybrid() : TestImplAliEmcalTrackSelection() {}
  TestImplAliEmcalTrackSelectionHybrid(const char* name, const char* period);
  virtual ~TestImplAliEmcalTrackSelectionHybrid() {}

 protected:
  virtual bool IsTrueTrack(const AliAODTrack* trk) const;

 private:
  TestImplAliEmcalTrackSelectionHybrid(const TestImplAliEmcalTrackSelectionHybrid&);
  TestImplAliEmcalTrackSelectionHybrid& operator=(const TestImplAliEmcalTrackSelectionHybrid&);

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalTrackSelectionHybrid, 1);
  /// \endcond
};

class TestImplAliEmcalTrackSelectionTPConly : public TestImplAliEmcalTrackSelection {
 public:
  TestImplAliEmcalTrackSelectionTPConly() : TestImplAliEmcalTrackSelection() {}
  TestImplAliEmcalTrackSelectionTPConly(const char* name, const char* period);
  virtual ~TestImplAliEmcalTrackSelectionTPConly() {}

 protected:
  virtual bool IsTrueTrack(const AliAODTrack* const trk) const;

 private:
  TestImplAliEmcalTrackSelectionTPConly(const TestImplAliEmcalTrackSelectionTPConly&);
  TestImplAliEmcalTrackSelectionTPConly& operator=(const TestImplAliEmcalTrackSelectionTPConly&);

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalTrackSelectionTPConly, 1);
  /// \endcond
};

class TestAliEmcalTrackSelection : public AliAnalysisTaskEmcalLight {
 public:
  TestAliEmcalTrackSelection();
  TestAliEmcalTrackSelection(const char* name);
  virtual ~TestAliEmcalTrackSelection();

  void SetPeriod(const char* period) { fPeriod = period; }
  void AddTestImpl(TestImplAliEmcalTrackSelection* test);
  void GenerateTestSuite();
  static TestAliEmcalTrackSelection* AddTestAliEmcalTrackSelection(const char* name);

 protected:
  bool EvaluateTest(TestImplAliEmcalTrackSelection* test);
  virtual void UserCreateOutputObjects();
  virtual bool Run();

 private:
  TString fPeriod;            ///< Period
  TObjArray* fTestSuite;      ///< Test suite
  THistManager* fTestResults; ///< Test result

  TestAliEmcalTrackSelection(const TestAliEmcalTrackSelection&);
  TestAliEmcalTrackSelection& operator=(const TestAliEmcalTrackSelection&);

  /// \cond CLASSIMP
  ClassDef(TestAliEmcalTrackSelection, 1);
  /// \endcond
};
}
}
#endif