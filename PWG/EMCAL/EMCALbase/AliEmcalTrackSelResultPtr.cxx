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
#include <iostream>
#include <TNamed.h> // for user object
#include <TClass.h>
#include "AliEmcalTrackSelResultPtr.h"
#include "AliVTrack.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalTrackSelResultPtr)
ClassImp(PWG::EMCAL::AliEmcalTrackSelResultUserPtr)
ClassImp(PWG::EMCAL::AliEmcalTrackSelResultUserStorage)
ClassImp(PWG::EMCAL::TestAliEmcalTrackSelResultPtr)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr() :
  TObject(),
  fTrack(nullptr),
  fSelectionResult(false),
  fUserInfo()
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(AliVTrack *trk, Bool_t selectionStatus, TObject * userinfo) :
  TObject(),
  fTrack(trk),
  fSelectionResult(selectionStatus),
  fUserInfo(userinfo)
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(const AliEmcalTrackSelResultPtr &ref) :
  TObject(ref),
  fTrack(ref.fTrack),
  fSelectionResult(ref.fSelectionResult),
  fUserInfo(ref.fUserInfo)
{

}

AliEmcalTrackSelResultPtr::AliEmcalTrackSelResultPtr(AliEmcalTrackSelResultPtr &&ref) :
  TObject(ref),
  fTrack(ref.fTrack),
  fSelectionResult(ref.fSelectionResult),
  fUserInfo()
{
  ref.fTrack = nullptr;
  fUserInfo = std::move(ref.fUserInfo);
}

AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultPtr::operator =(const AliEmcalTrackSelResultPtr &ref){
  TObject::operator =(ref);
  if(this != &ref){
    fTrack = ref.fTrack;
    fSelectionResult = ref.fSelectionResult;
    fUserInfo = ref.fUserInfo;
  }
  return *this;
}

AliEmcalTrackSelResultPtr &AliEmcalTrackSelResultPtr::operator =(AliEmcalTrackSelResultPtr &&ref){
  TObject::operator =(ref);
  if(this != &ref){
    fTrack = ref.fTrack;
    fSelectionResult = ref.fSelectionResult;
    fUserInfo = std::move(ref.fUserInfo);

    ref.fTrack = nullptr;
  }
  return *this;
}

Bool_t AliEmcalTrackSelResultPtr::operator ==(const AliEmcalTrackSelResultPtr &other) const {
  return fTrack == other.fTrack;
}

Bool_t AliEmcalTrackSelResultPtr::operator <(const AliEmcalTrackSelResultPtr &other) const {
  return fTrack < other.fTrack;
}

Bool_t AliEmcalTrackSelResultPtr::IsEqual(const TObject *o) const {
  const AliEmcalTrackSelResultPtr *otherobj = static_cast<const AliEmcalTrackSelResultPtr *>(o);
  if(!otherobj) return false;
  return *this == *otherobj;
}

Int_t AliEmcalTrackSelResultPtr::Compare(const TObject *o) const {
  const AliEmcalTrackSelResultPtr *otherobj = static_cast<const AliEmcalTrackSelResultPtr *>(o);
  if(!otherobj) return 1;
  if (*this == *otherobj) return 0;
  if (*this < *otherobj) return -1;
  return 1;
}

AliVTrack *AliEmcalTrackSelResultPtr::operator *() const {
  return fTrack;
}

AliVTrack *AliEmcalTrackSelResultPtr::operator ->() const {
  return fTrack;
}

void AliEmcalTrackSelResultPtr::PrintStream(std::ostream &stream) const {
  stream  << "Track selection result for track with address " << fTrack
          << ": Selection status: " << (fSelectionResult ? "true" : "false")
          << ", HasUserInfo: " << (fUserInfo.GetData() ? "Yes" : "No");
}

std::ostream &PWG::EMCAL::operator<<(std::ostream &stream, const PWG::EMCAL::AliEmcalTrackSelResultPtr &o){
  o.PrintStream(stream);
  return stream;
}

/*************************************************************
 *                                                           *
 * Content of class AliEmcalTrackSelResultUserStorage        *
 *                                                           *
 *************************************************************/

AliEmcalTrackSelResultUserStorage::AliEmcalTrackSelResultUserStorage():
  TObject(),
  fData(nullptr),
  fReferenceCount(0)
{

}

AliEmcalTrackSelResultUserStorage::AliEmcalTrackSelResultUserStorage(TObject *data):
  TObject(),
  fData(data),
  fReferenceCount(1)
{
}

AliEmcalTrackSelResultUserStorage::~AliEmcalTrackSelResultUserStorage(){
  //std::cout << "Destructing user storage with reference count 0" << std::endl;
  if(fData) delete fData;
}

void AliEmcalTrackSelResultUserStorage::Connect() {
  fReferenceCount++; 
}

void AliEmcalTrackSelResultUserStorage::Release(){
  fReferenceCount--;
  if(!fReferenceCount) delete this;
}

/*************************************************************
 *                                                           *
 * Content of class AliEmcalTrackSelResultUserPtr            *
 *                                                           *
 *************************************************************/

AliEmcalTrackSelResultUserPtr::AliEmcalTrackSelResultUserPtr():
  TObject(),
  fUserStorage(nullptr)
{

}

AliEmcalTrackSelResultUserPtr::AliEmcalTrackSelResultUserPtr(TObject *data):
  TObject(),
  fUserStorage(nullptr)
{
  if(data) fUserStorage = new AliEmcalTrackSelResultUserStorage(data);
}

AliEmcalTrackSelResultUserPtr::AliEmcalTrackSelResultUserPtr(const AliEmcalTrackSelResultUserPtr &ref):
  TObject(ref),
  fUserStorage(nullptr)
{
  if(ref.fUserStorage) {
    fUserStorage = ref.fUserStorage;
    fUserStorage->Connect();
  }
}

AliEmcalTrackSelResultUserPtr::AliEmcalTrackSelResultUserPtr(AliEmcalTrackSelResultUserPtr &&ref):
  TObject(ref),
  fUserStorage(nullptr)
{
  if(ref.fUserStorage) {
    // reference count does not change in case of move semantics
    fUserStorage = ref.fUserStorage;
    ref.fUserStorage = nullptr;
  }
}

AliEmcalTrackSelResultUserPtr &AliEmcalTrackSelResultUserPtr::operator=(const AliEmcalTrackSelResultUserPtr &ref){
  TObject::operator=(ref);

  if(&ref != this){
    if(fUserStorage) {
      fUserStorage->Release();
     // if(fUserStorage->GetReferenceCount() < 1) delete fUserStorage;
    } 

    if(ref.fUserStorage) {
      fUserStorage = ref.fUserStorage;
      fUserStorage->Connect();
    }
  }
  return *this;
}

AliEmcalTrackSelResultUserPtr &AliEmcalTrackSelResultUserPtr::operator=(AliEmcalTrackSelResultUserPtr &&ref){
  TObject::operator=(ref);

  if(&ref != this){
    if(fUserStorage) {
      fUserStorage->Release();
      //if(fUserStorage->GetReferenceCount() < 1) delete fUserStorage;
    } 

    if(ref.fUserStorage) {
      // reference count does not change in case of move semantics
      fUserStorage = ref.fUserStorage;
      ref.fUserStorage = nullptr;
    }
  }
  return *this;
}

AliEmcalTrackSelResultUserPtr::~AliEmcalTrackSelResultUserPtr(){
  if(fUserStorage) {
    fUserStorage->Release();
    // The last pointer connected to the storage deletes it
    //if(fUserStorage->GetReferenceCount() < 1) delete fUserStorage;
  }
}

bool TestAliEmcalTrackSelResultPtr::RunAllTests() const {
  return TestOperatorBool() && TestUserInfo() && TestCopyConstructor() && TestOperatorAssign();
}

bool TestAliEmcalTrackSelResultPtr::TestOperatorBool() const {
  AliEmcalTrackSelResultPtr testtrue(nullptr, true), testfalse(nullptr, false);

  bool testresult(true);
  if(!(testtrue == true)) testresult = false;
  if(testfalse == true) testresult = false;
  return testresult;
}

bool TestAliEmcalTrackSelResultPtr::TestCopyConstructor() const {
  int failure(0);
  for(int i = 0; i < 10; i++) {
    TNamed *payloadTrue = new TNamed("truewith", "true, with object"), 
           *payloadFalse = new TNamed("falsewith", "false, with user object");
    AliEmcalTrackSelResultPtr truewith(nullptr, true, payloadTrue),
                              truewithout(nullptr, true),
                              falsewith(nullptr, false, payloadFalse),
                              falsewithout(nullptr, false);
    // make copies
    AliEmcalTrackSelResultPtr cpytruewith(truewith), cpyfalsewith(falsewith), cpytruewithout(truewithout), cpyfalsewithout(falsewithout);
    if(!(AssertBool(cpytruewith, true) && AssertBool(cpytruewithout, true) && AssertBool(cpyfalsewith, false) && AssertBool(cpyfalsewithout, false))) failure++;
    if(!(AssertPayload(cpytruewith, payloadTrue) && AssertPayload(cpytruewithout, nullptr) && AssertPayload(cpyfalsewith, payloadFalse) && AssertPayload(cpyfalsewithout, nullptr))) failure++;
  }
  return failure == 0;
}

bool TestAliEmcalTrackSelResultPtr::TestOperatorAssign() const {
  int failure(0);
  for(int i = 0; i < 10; i++) {
    TNamed *payloadTrue = new TNamed("truewith", "true, with object"), 
           *payloadFalse = new TNamed("falsewith", "false, with user object");
    AliEmcalTrackSelResultPtr truewith(nullptr, true, payloadTrue),
                              truewithout(nullptr, true),
                              falsewith(nullptr, false, payloadFalse),
                              falsewithout(nullptr, false);
    // make assignments
    AliEmcalTrackSelResultPtr asgtruewith = truewith, asgfalsewith = falsewith, asgtruewithout = truewithout, asgfalsewithout = falsewithout;
    if(!(AssertBool(asgtruewith, true) && AssertBool(asgtruewithout, true) && AssertBool(asgfalsewith, false) && AssertBool(asgfalsewithout, false))) failure++;
    if(!(AssertPayload(asgtruewith, payloadTrue) && AssertPayload(asgtruewithout, nullptr) && AssertPayload(asgfalsewith, payloadFalse) && AssertPayload(asgfalsewithout, nullptr))) failure++;
  }
  return failure == 0;
}

bool TestAliEmcalTrackSelResultPtr::TestUserInfo() const {
  int failure(0);
  for(int i = 0; i < 10; i++) {
    AliEmcalTrackSelResultPtr testwith(nullptr, true, new TNamed("testuserobject", "Test userobject"));
    if(!testwith.GetUserInfo()) failure++;

    AliEmcalTrackSelResultPtr testwithout(nullptr, true);
    if(testwithout.GetUserInfo()) failure++;
  }
  return failure == 0;
}

bool TestAliEmcalTrackSelResultPtr::AssertBool(const AliEmcalTrackSelResultPtr &test, bool assertval) const {
  return test.GetSelectionResult() == assertval;
}

bool TestAliEmcalTrackSelResultPtr::AssertPayload(const AliEmcalTrackSelResultPtr &test, void *payload) const {
  return test.GetUserInfo() == reinterpret_cast<const TObject *>(payload);
}