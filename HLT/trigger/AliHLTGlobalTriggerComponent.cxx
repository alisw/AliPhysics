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

/// @file   AliHLTGlobalTriggerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Implementation of the AliHLTGlobalTriggerComponent component class.
///
/// The AliHLTGlobalTriggerComponentComponent class applies the global HLT trigger to all
/// trigger information produced by components deriving from AliHLTTrigger.

#include "AliHLTGlobalTriggerComponent.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTGlobalTrigger.h"
#include "TUUID.h"
#include "TROOT.h"
#include "TSystem.h"
#include <fstream>
#include <cerrno>

ClassImp(AliHLTGlobalTriggerComponent)


AliHLTGlobalTriggerComponent::AliHLTGlobalTriggerComponent() :
	AliHLTTrigger(),
	fTrigger(NULL)
{
  // Default constructor.
}


AliHLTGlobalTriggerComponent::~AliHLTGlobalTriggerComponent()
{
  // Default destructor.
  
  if (fTrigger != NULL) delete fTrigger;
}


void AliHLTGlobalTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // Returns the output data size estimate.

  constBase = sizeof(AliHLTGlobalTriggerDecision);
  inputMultiplier = 1;
}


Int_t AliHLTGlobalTriggerComponent::DoInit(int /*argc*/, const char** /*argv*/)
{
  // Initialises the global trigger component.
  
  AliHLTTriggerMenu* menu = NULL;
  TString classname;
  int result = GenerateTrigger(menu, classname);
  if (result != 0) return result;
  
  fTrigger = AliHLTGlobalTrigger::CreateNew(classname.Data());
  if (fTrigger == NULL)
  {
    HLTError("Could not create a new instance of '%s'.", classname.Data());
    return -EIO;
  }
  
  return 0;
}


Int_t AliHLTGlobalTriggerComponent::DoDeinit()
{
  // Cleans up the global trigger component.
  
  if (fTrigger != NULL)
  {
    delete fTrigger;
    fTrigger = NULL;
  }
  
  return 0;
}


AliHLTComponent* AliHLTGlobalTriggerComponent::Spawn()
{
  // Creates a new object instance.
  
  return new AliHLTGlobalTriggerComponent;
}


int AliHLTGlobalTriggerComponent::DoTrigger()
{
  // This method will apply the global trigger decision.

  if (fTrigger == NULL)
  {
    HLTFatal("Global trigger implementation object is NULL!");
    return -EIO;
  }
  
  fTrigger->NewEvent();
  
  // Fill in the input data.
  const TObject* obj = GetFirstInputObject();
  while (obj != NULL)
  {
    if (obj->IsA() == AliHLTTriggerDecision::Class())
    {
      const AliHLTTriggerDecision* decision = static_cast<const AliHLTTriggerDecision*>(obj);
      fTrigger->Add(decision);
    }
    else
    {
      fTrigger->Add(obj, GetDataType(), GetSpecification());
    }
    obj = GetNextInputObject();
  }

  // Apply the trigger.
  TriggerEvent(fTrigger->CalculateTriggerDecision());
  return 0;
}


int AliHLTGlobalTriggerComponent::GenerateTrigger(const AliHLTTriggerMenu* /*menu*/, TString& name)
{
  // Generates the global trigger class that will implement the specified trigger menu.
  
  // Create a new UUID and replace the '-' characters with '_' to make it a valid
  // C++ symbol name.
  TUUID uuid;
  TString uuidstr = uuid.AsString();
  for (Int_t i = 0; i < uuidstr.Length(); i++)
  {
    if (uuidstr[i] == '-') uuidstr[i] = '_';
  }
  
  // Create the name of the new class.
  name = "AliHLTGlobalTriggerImpl_";
  name += uuidstr;
  TString filename = name + ".cxx";
  
  fstream code(filename.Data(), ios_base::out | ios_base::trunc);
  if (not code.good())
  {
    HLTError("Could not open file '%s' for writing.", filename.Data());
    return -EIO;
  }
  
  code << "#include \"AliHLTGlobalTrigger.h\"" << endl;
  code << "#include \"AliHLTGlobalTriggerDecision.h\"" << endl;
  code << "class " << name << " : public AliHLTGlobalTrigger" << endl;
  code << "{" << endl;
  code << "public:" << endl;
  code << "  " << name << "() : AliHLTGlobalTrigger(), fDecision() {" << endl;
  code << "  }" << endl;
  code << "  virtual ~" << name << "() {" << endl;
  code << "  }" << endl;
  code << "  virtual void NewEvent() {" << endl;
  //code << "    ;" << endl;
  code << "  }" << endl;
  code << "  virtual void Add(const AliHLTTriggerDecision* decision) {" << endl;
  //code << "    ;" << endl;
  code << "  }" << endl;
  code << "  virtual void Add(const TObject* object, const AliHLTComponentDataType& type, AliHLTUInt32_t spec) {" << endl;
  //code << "    ;" << endl;
  code << "  }" << endl;
  code << "  virtual AliHLTGlobalTriggerDecision* CalculateTriggerDecision() {" << endl;
  code << "    return &fDecision;" << endl;
  code << "  }" << endl;
  code << "  class FactoryImpl : public AliHLTGlobalTrigger::Factory" << endl;
  code << "  {" << endl;
  code << "  public:" << endl;
  code << "    virtual const char* ClassName() const {" << endl;
  code << "      return \"" << name << "\";" << endl;
  code << "    }" << endl;
  code << "    virtual AliHLTGlobalTrigger* New() const {" << endl;
  code << "      return new " << name << "();" << endl;
  code << "    }" << endl;
  code << "  private:" << endl;
  code << "    static FactoryImpl fFactoryImpl; // for registration only." << endl;
  code << "  };" << endl;
  code << "private:" << endl;
  code << "  AliHLTGlobalTriggerDecision fDecision;" << endl;
  code << "};" << endl;
  code << name << "::FactoryImpl " << name << "::FactoryImpl::fFactoryImpl;" << endl;
  
  TString includePath = "-I${ALICE_ROOT}/include -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/trigger";
  gSystem->SetIncludePath(includePath);
  
  TString cmd = ".L ";
  cmd += filename;
  cmd += "++";
  gROOT->ProcessLine(cmd);
  
  code.close();
  
  return 0;
}

