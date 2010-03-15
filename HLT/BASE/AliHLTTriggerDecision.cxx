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

/// @file   AliHLTTriggerDecision.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   21 Nov 2008
/// @brief  Implementation of the AliHLTTriggerDecision class.
/// 
/// The trigger decision class stores the HLT decision from an AliHLTTrigger component.

#include "AliHLTTriggerDecision.h"
#include "Riostream.h"

ClassImp(AliHLTTriggerDecision)


AliHLTTriggerDecision::AliHLTTriggerDecision() :
  TObject(),
  fName(),
  fDescription(),
  fTriggerDomain()
{
  // Default constructor.
}


AliHLTTriggerDecision::AliHLTTriggerDecision(const AliHLTTriggerDecision& obj) :
  TObject(obj),
  fName(obj.fName),
  fDescription(obj.fDescription),
  fTriggerDomain(obj.fTriggerDomain)
{
  // Copy constructor performs a deep copy.
  
  // The following is a backward compatibility fix to be able to read trigger
  // decisions recorded before the fix in AliRoot trunk rev. 35998 correctly.
  if (obj.TestBits(15) == 15)
  {
    ResetBit(15);  // We can clear 'this' objects bit because we performed a deep copy.
    SetBit(BIT(15));
  }
}


AliHLTTriggerDecision::AliHLTTriggerDecision(bool result, const char* name) :
  TObject(),
  fName(name),
  fDescription(),
  fTriggerDomain()
{
  // Constructor specifying the name and result of the trigger decision.
  
  Result(result);
}


AliHLTTriggerDecision::AliHLTTriggerDecision(
    bool result, const char* name,
    const AliHLTTriggerDomain& triggerDomain,
    const char* description
  ) :
  TObject(),
  fName(name),
  fDescription(description),
  fTriggerDomain(triggerDomain)
{
  // Constructor specifying all information fields.
  
  Result(result);
}


AliHLTTriggerDecision::~AliHLTTriggerDecision()
{
  // Default destructor.
}


bool AliHLTTriggerDecision::Result() const
{
  // Returns the result of the trigger decision.
  
  // The following is a backward compatibility fix to be able to read trigger
  // decisions recorded before the fix in AliRoot trunk rev. 35998 correctly.
  if (TestBits(15) == 15) return true;
  
  return TestBit(BIT(15)) == 1;
}


void AliHLTTriggerDecision::Result(bool value)
{
  // Sets the result of the trigger decision.
  SetBit(BIT(15), value);
  
  // The following is a backward compatibility fix to be able to read trigger
  // decisions recorded before the fix in AliRoot trunk rev. 35998 correctly.
  // It looks like bit 1 and 2 of fBits are not used in the case of the
  // AliHLTTriggerDecision class, so reset those to prevent "TestBits(15) == 15"
  // from succeeding in the "Result() const" method above.
  // We do not touch the other two bits because they could affect memory handling
  // and cleanup.
  if (TestBits(15) == 15) ResetBit(6);
}


void AliHLTTriggerDecision::ReadoutList(const AliHLTReadoutList& value)
{
  // Replaces the readout list in the trigger domain with the new value.
  
  AliHLTReadoutList fullReadout = ~ AliHLTReadoutList(0x0);
  fTriggerDomain.Remove(fullReadout);
  fTriggerDomain.Add(value);
}


void AliHLTTriggerDecision::Print(Option_t* option) const
{
  // Prints the contents of the trigger decision.
  
  cout << "Trigger (" << fName.Data() << ") result = " << Result() << endl;
  TString opt(option);
  if (opt.Contains("short")) return;
  cout << "Description = \"" << fDescription.Data() << "\"" << endl;
  fTriggerDomain.Print();
}

void AliHLTTriggerDecision::Copy(TObject &object) const
{
  // copy this to the specified object

  AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(&object);
  if (pDecision) {
    // copy members if target is a AliHLTTriggerDecision
    *pDecision=*this;
  }

  // copy the base class
  TObject::Copy(object);
}

TObject *AliHLTTriggerDecision::Clone(const char */*newname*/) const
{
  // create a new clone, classname is ignored

  return new AliHLTTriggerDecision(*this);
}

Option_t *AliHLTTriggerDecision::GetOption() const
{
  // Return the result of the trigger.
  // "0" or "1"
  if (Result()) return "1";
  return "0";
}


AliHLTTriggerDecision& AliHLTTriggerDecision::operator = (const AliHLTTriggerDecision& obj)
{
  // Assignment operator performs a deep copy.
  
  if (this == &obj) return *this;
  
  TObject::operator = (obj);
  // The following is a backward compatibility fix to be able to read trigger
  // decisions recorded before the fix in AliRoot trunk rev. 35998 correctly.
  if (obj.TestBits(15) == 15)
  {
    ResetBit(15);  // We can clear 'this' objects bit because we performed a deep copy.
    SetBit(BIT(15));
  }
  
  fName = obj.fName;
  fDescription = obj.fDescription;
  fTriggerDomain = obj.fTriggerDomain;
  return *this;
}


void AliHLTTriggerDecision::Clear(Option_t* option)
{
  // Clears the trigger domain and resets the decision result.
  
  Result(false);
  fTriggerDomain.Clear(option);
}
