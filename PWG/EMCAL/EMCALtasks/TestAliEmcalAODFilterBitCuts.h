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
#ifndef PWG_EMCAL_TESTALIEMCALAODFILTERBITCUTS_H_
#define PWG_EMCAL_TESTALIEMCALAODFILTERBITCUTS_H_

#include "AliAnalysisTaskEmcalLight.h"
#include <TNamed.h>

class THistManager;
class TObjArray;
class AliAODEvent;
class AliAODTrack;

namespace PWG {
namespace EMCAL {

class AliEmcalAODFilterBitCuts;

class TestImplAliEmcalAODFilterBitCuts : public TNamed {
public:
  TestImplAliEmcalAODFilterBitCuts();
  TestImplAliEmcalAODFilterBitCuts(const char *name);
  virtual ~TestImplAliEmcalAODFilterBitCuts();

  bool RunTest(const AliAODEvent *const ev);

protected:
  virtual bool IsTrueTrack(const AliAODTrack *const trk) const = 0;

  AliEmcalAODFilterBitCuts          *fTestObject;   ///< Object to be tested

private:
  TestImplAliEmcalAODFilterBitCuts(const TestImplAliEmcalAODFilterBitCuts &);
  TestImplAliEmcalAODFilterBitCuts &operator=(const TestImplAliEmcalAODFilterBitCuts &);

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalAODFilterBitCuts, 1);
  /// \endcond
};

class TestImplAliEmcalAODFilterBitCutsTPCconstrained : public TestImplAliEmcalAODFilterBitCuts {
public:
  TestImplAliEmcalAODFilterBitCutsTPCconstrained() : TestImplAliEmcalAODFilterBitCuts() {}
  TestImplAliEmcalAODFilterBitCutsTPCconstrained(const char *name);
  virtual ~TestImplAliEmcalAODFilterBitCutsTPCconstrained(){}

protected:
  virtual bool IsTrueTrack(const AliAODTrack *const trk) const;

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalAODFilterBitCutsTPCconstrained, 1);
  /// \endcond
};

class TestImplAliEmcalAODFilterBitCutsHybrid : public TestImplAliEmcalAODFilterBitCuts {
public:
  TestImplAliEmcalAODFilterBitCutsHybrid() : TestImplAliEmcalAODFilterBitCuts() {}
  TestImplAliEmcalAODFilterBitCutsHybrid(const char *name);
  virtual ~TestImplAliEmcalAODFilterBitCutsHybrid(){}

protected:
  virtual bool IsTrueTrack(const AliAODTrack *const trk) const;

  /// \cond CLASSIMP
  ClassDef(TestImplAliEmcalAODFilterBitCutsHybrid, 1);
  /// \endcond
};

class TestAliEmcalAODFilterBitCuts : public AliAnalysisTaskEmcalLight{
public:
  TestAliEmcalAODFilterBitCuts();
  TestAliEmcalAODFilterBitCuts(const char *testname);
  virtual ~TestAliEmcalAODFilterBitCuts();

  void AddTestImpl(TestImplAliEmcalAODFilterBitCuts *test);
  void GenerateTestSuite();

  static TestAliEmcalAODFilterBitCuts *AddTestAliEmcalAODFilterBitCuts(const char *testname);

protected:
  virtual void UserCreateOutputObjects();
  virtual bool Run();

private:
  void EvaluateTest(TestImplAliEmcalAODFilterBitCuts *test, const AliAODEvent *const ev);

  TObjArray                             *fTestSuite;    ///< Test suite to be executed
  THistManager                          *fTestStatus;   ///< Histograms with test results

  TestAliEmcalAODFilterBitCuts(const TestAliEmcalAODFilterBitCuts &);
  TestAliEmcalAODFilterBitCuts &operator=(const TestAliEmcalAODFilterBitCuts &);

  /// \cond CLASSIMP
  ClassDef(TestAliEmcalAODFilterBitCuts, 1);
  /// \endcond
};

} /* namespace EMCAL */
} /* namespace PWG */

#endif /* PWG_EMCAL_TESTALIEMCALAODFILTERBITCUTS_H_ */
