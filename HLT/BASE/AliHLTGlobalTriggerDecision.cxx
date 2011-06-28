// $Id$
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTGlobalTriggerDecision.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Implementation of the AliHLTGlobalTriggerDecision class.
/// 
/// The global trigger decision class stores the global HLT decision.

#include "AliHLTGlobalTriggerDecision.h"
#include "Riostream.h"
#include "TClass.h"
#include "AliHLTMisc.h"
#include "AliHLTLogging.h"

ClassImp(AliHLTGlobalTriggerDecision)


namespace {
const char* kgNULLObjMessage =
  "The global decision object contains NULL pointers in the fContributingTriggers array."
  " This is unexpected and probably indicates a serious problem.";
const char* kgSplitModeMessage =
  "The global decision object contains NULL pointers in the fInputObjects array."
  " This means that it was written to a TTree branch with the split-mode > 0."
  " This causes the custom streamer to be skipped and prevents the"
  " fInputObjects array to be setup properly. In addition this can cause memory leaks.";
};


AliHLTGlobalTriggerDecision::AliHLTGlobalTriggerDecision() :
AliHLTTriggerDecision(0, "HLTGlobalTrigger"),
  fContributingTriggers(AliHLTTriggerDecision::Class()),
  fInputObjects(),
  fCounters()
{
  // Default constructor.
  
  // We set the ownership to false since the objects are deleted manually by this
  // class in DeleteInputObjects().
  fInputObjects.SetOwner(kFALSE);
}


AliHLTGlobalTriggerDecision::AliHLTGlobalTriggerDecision(
    bool result, const AliHLTTriggerDomain& triggerDomain, const char* description
  ) :
  AliHLTTriggerDecision(result, "HLTGlobalTrigger", triggerDomain, description),
  fContributingTriggers(AliHLTTriggerDecision::Class()),
  fInputObjects(),
  fCounters()
{
  // Constructor specifying multiple information fields.
  
  Result(result);
  // We set the ownership to false since the objects are deleted manually by this
  // class in DeleteInputObjects().
  fInputObjects.SetOwner(kFALSE);
}


AliHLTGlobalTriggerDecision::~AliHLTGlobalTriggerDecision()
{
  // Default destructor.
  
  DeleteInputObjects();
}


void AliHLTGlobalTriggerDecision::Print(Option_t* option) const
{
  // Prints the contents of the trigger decision.
  
  TString opt(option);
  if (opt.Contains("compact"))
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print("");
  }
  else if (opt.Contains("short"))
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print(option);
    cout << "#################### Input trigger decisions ####################" << endl;
    for (Int_t i = 0; i < NumberOfTriggerInputs(); i++)
    {
      if (TriggerInput(i) != NULL)
      {
        TriggerInput(i)->Print(option);
      }
      else
      {
        AliHLTLogging log;
        log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgNULLObjMessage);
      }
    }
    if (NumberOfTriggerInputs() == 0)
    {
      cout << "(none)" << endl;
    }
  }
  else if (opt.Contains("counters"))
  {
    cout << "Counter\tValue" << endl;
    for (Int_t i = 0; i < fCounters.GetSize(); i++)
    {
      cout << i << "\t" << fCounters[i] << endl;
    }
    if (fCounters.GetSize() == 0)
    {
      cout << "(none)" << endl;
    }
  }
  else
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print(option);
    cout << "#################### Input trigger decisions ####################" << endl;
    for (Int_t i = 0; i < NumberOfTriggerInputs(); i++)
    {
      cout << "-------------------- Input trigger decision " << i << " --------------------" << endl;
      if (TriggerInput(i) != NULL)
      {
        TriggerInput(i)->Print(option);
      }
      else
      {
        AliHLTLogging log;
        log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgNULLObjMessage);
      }
    }
    if (NumberOfTriggerInputs() == 0)
    {
      cout << "(none)" << endl;
    }
    cout << "###################### Other input objects ######################" << endl;
    for (Int_t i = 0; i < NumberOfInputObjects(); i++)
    {
      cout << "------------------------ Input object " << i << " ------------------------" << endl;
      if (InputObject(i) != NULL)
      {
        InputObject(i)->Print(option);
      }
      else
      {
        AliHLTLogging log;
        log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgSplitModeMessage);
      }
    }
    if (NumberOfInputObjects() == 0)
    {
      cout << "(none)" << endl;
    }
    cout << "#################### Event class counters ####################" << endl;
    cout << "Counter\tValue" << endl;
    for (Int_t i = 0; i < fCounters.GetSize(); i++)
    {
      cout << i << "\t" << fCounters[i] << endl;
    }
    if (fCounters.GetSize() == 0)
    {
      cout << "(none)" << endl;
    }
  }
}

void AliHLTGlobalTriggerDecision::Copy(TObject &object) const
{
  // copy this to the specified object

  if (object.IsA() == AliHLTMisc::Instance().IsAliESDHLTDecision()) {
    AliHLTMisc::Instance().Copy(this, &object);
    return;
  }

  AliHLTGlobalTriggerDecision* pDecision=dynamic_cast<AliHLTGlobalTriggerDecision*>(&object);
  if (pDecision)
  {
    // copy members if target is a AliHLTGlobalTriggerDecision
    *pDecision=*this;
  }
  // copy the base class
  AliHLTTriggerDecision::Copy(object);
}

TObject *AliHLTGlobalTriggerDecision::Clone(const char */*newname*/) const
{
  // create a new clone, classname is ignored

  return new AliHLTGlobalTriggerDecision(*this);
}

AliHLTGlobalTriggerDecision::AliHLTGlobalTriggerDecision(const AliHLTGlobalTriggerDecision& src) :
  AliHLTTriggerDecision(src),
  fContributingTriggers(AliHLTTriggerDecision::Class()),
  fInputObjects(),
  fCounters()
{
  // Copy constructor performs a deep copy.
  
  // We set the ownership to false since the objects are deleted manually by this
  // class in DeleteInputObjects().
  fInputObjects.SetOwner(kFALSE);
  
  *this=src;
}

AliHLTGlobalTriggerDecision& AliHLTGlobalTriggerDecision::operator=(const AliHLTGlobalTriggerDecision& src)
{
  // assignment operator performs a deep copy.

  fContributingTriggers.Delete();
  for (int triggerInput=0; triggerInput<src.NumberOfTriggerInputs(); triggerInput++) {
    const AliHLTTriggerDecision* pTriggerObject=src.TriggerInput(triggerInput);
    if (pTriggerObject != NULL)
    {
      // the AddTriggerInput function uses the copy constructor and
      // makes a new object from the reference
      AddTriggerInput(*pTriggerObject);
    }
    else
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgNULLObjMessage);
    }
  }

  DeleteInputObjects();
  for (int inputObject=0; inputObject<src.NumberOfInputObjects(); inputObject++) {
    const TObject* pInputObject=src.InputObject(inputObject);
    if (pInputObject != NULL)
    {
      // the AddInputObject function uses Clone() and
      // makes a new object from the reference
      AddInputObject(pInputObject);
    }
    else
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgSplitModeMessage);
    }
  }
  
  SetCounters(src.Counters());

  return *this;
}


void AliHLTGlobalTriggerDecision::AddInputObject(const TObject* object)
{
  // Adds an object to the list of input objects considered in the global trigger.
  
  if (object == NULL) return;
  TObject* obj = object->Clone();
  obj->SetBit(kCanDelete);
  fInputObjects.Add(obj);
}


void AliHLTGlobalTriggerDecision::AddInputObjectRef(TObject* object, bool own)
{
  // Adds an object to the list of input objects considered in the global trigger.
  
  if (object == NULL) return;
  if (own)
  {
    object->SetBit(kCanDelete);
  }
  else
  {
    object->ResetBit(kCanDelete);
  }
  fInputObjects.Add(object);
}


void AliHLTGlobalTriggerDecision::SetCounters(const TArrayL64& counters, Long64_t eventCount)
{
  // Sets the counter array.
  // If the number of events is specified, an additional counter is added at the end.
  fCounters = counters;
  if (eventCount>=0) {
    int size=fCounters.GetSize();
    fCounters.Set(size+1);
    fCounters[size]=eventCount;
  }
}


void AliHLTGlobalTriggerDecision::Clear(Option_t* option)
{
  // Clears the trigger domain and resets the decision result.
  
  AliHLTTriggerDecision::Clear(option);
  // because of TClonesArray members in AliHLTTriggerDecision it is not
  // enough to call Clear. Delete will also invoke the destructor of the
  // elements which is necessary to do the proper internal cleanup
  fContributingTriggers.Delete();
  DeleteInputObjects();
  fCounters.Set(0);
}


TObject* AliHLTGlobalTriggerDecision::FindObject(const char* name) const
{
  // Finds the first object in fContributingTriggers or fInputObjects that has the given name.
  
  TObject* result = fContributingTriggers.FindObject(name);
  if (result != NULL) return result;
  return fInputObjects.FindObject(name);
}


TObject* AliHLTGlobalTriggerDecision::FindObject(const TObject* obj) const
{
  // Finds the first object in fContributingTriggers or fInputObjects that matches
  // based on a IsEqual() comparison.
  
  TObject* result = fContributingTriggers.FindObject(obj);
  if (result != NULL) return result;
  return fInputObjects.FindObject(obj);
}


void AliHLTGlobalTriggerDecision::DeleteInputObjects()
{
  // Deletes the objects marked with kCanDelete in fInputObjects and clears the array.
  
  for (Int_t i = 0; i < NumberOfInputObjects(); i++)
  {
    TObject* obj = fInputObjects.UncheckedAt(i);
    if (obj == NULL)
    {
      AliHLTLogging log;
      log.LoggingVarargs(kHLTLogError, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__, kgSplitModeMessage);
      continue;
    }
    if (obj->TestBit(kCanDelete)) delete obj;
  }
  fInputObjects.Clear();
}


void AliHLTGlobalTriggerDecision::MarkInputObjectsAsOwned()
{
  // Marks all input objects as owned.

  // We must mark all the objects that were read into fInputObjects as owned.
  // Otherwise we will have a memory leak in DeleteInputObjects.
  bool loggedWarning = false;
  for (Int_t i = 0; i < fInputObjects.GetEntriesFast(); ++i)
  {
    TObject* obj = fInputObjects.UncheckedAt(i);
    // We must check if the object pointer is NULL. This could happen because the
    // class dictionary has not been loaded, so the ReadClassBuffer streamer just
    // silently skips the object but fills the fInputObjects array with a NULL pointer.
    if (obj == NULL)
    {
      fInputObjects.RemoveAt(i);
      if (not loggedWarning)
      {
        AliHLTLogging log;
        log.LoggingVarargs(kHLTLogWarning, this->ClassName(), FUNCTIONNAME(), __FILE__, __LINE__,
          "The global trigger decision contains NULL pointers in the input object array."
          " This is probably due to the fact that some class dictionaries have not been loaded."
          " Will just remove the NULL pointer and continue."
        );
        loggedWarning = true;  // Prevent multiple warnings, one is enough.
      }
    }
    else
    {
      obj->SetBit(kCanDelete);
    }
  }
  // Compress the input object array to prevent any seg-faults due to access of
  // NULL pointers if the objects were not loaded due to missing dictionaries.
  fInputObjects.Compress();
}

#if ROOT_VERSION_CODE < ROOT_VERSION(5,26,0)
void AliHLTGlobalTriggerDecision::Streamer(TBuffer &b)
{
   // Stream an object of class AliHLTGlobalTriggerDecision.

   if (b.IsReading())
   {
     b.ReadClassBuffer(AliHLTGlobalTriggerDecision::Class(), this);
     MarkInputObjectsAsOwned();
   }
   else
   {
     b.WriteClassBuffer(AliHLTGlobalTriggerDecision::Class(), this);
   }
}
#endif // ROOT version check.
