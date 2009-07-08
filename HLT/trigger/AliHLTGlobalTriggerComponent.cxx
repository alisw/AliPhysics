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
#include "AliHLTGlobalTriggerConfig.h"
#include "AliHLTTriggerMenu.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "TUUID.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TRegexp.h"
#include "TClonesArray.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TInterpreter.h"
#include <fstream>
#include <cerrno>

ClassImp(AliHLTGlobalTriggerComponent)

const char* AliHLTGlobalTriggerComponent::fgkTriggerMenuCDBPath = "HLT/ConfigHLT/HLTGlobalTrigger";


AliHLTGlobalTriggerComponent::AliHLTGlobalTriggerComponent() :
	AliHLTTrigger(),
	fTrigger(NULL),
	fDebugMode(false),
	fCodeFileName()
{
  // Default constructor.
  
  ClearInfoForNewEvent(false);
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


Int_t AliHLTGlobalTriggerComponent::DoInit(int argc, const char** argv)
{
  // Initialises the global trigger component.
  
  fDebugMode = false;
  const char* configFileName = NULL;
  const char* codeFileName = NULL;
  TString classname;
  TClonesArray includePaths(TObjString::Class());
  TClonesArray includeFiles(TObjString::Class());
  
  for (int i = 0; i < argc; i++)
  {
    if (strcmp(argv[i], "-config") == 0)
    {
      if (configFileName != NULL)
      {
        HLTWarning("Trigger configuration macro was already specified."
                   " Will replace previous value given by -config."
        );
      }
      if (argc <= i+1)
      {
        HLTError("The trigger configuration macro filename was not specified." );
        return -EINVAL;
      }
      configFileName = argv[i+1];
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-includepath") == 0)
    {
      if (argc <= i+1)
      {
        HLTError("The include path was not specified." );
        return -EINVAL;
      }
      new (includePaths[includePaths.GetEntriesFast()]) TObjString(argv[i+1]);
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-include") == 0)
    {
      if (argc <= i+1)
      {
        HLTError("The include file name was not specified." );
        return -EINVAL;
      }
      new (includeFiles[includeFiles.GetEntriesFast()]) TObjString(argv[i+1]);
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-debug") == 0)
    {
      if (fDebugMode == true)
      {
        HLTWarning("The debug flag was already specified. Ignoring this instance.");
      }
      fDebugMode = true;
      i++;
      continue;
    }
    
    if (strcmp(argv[i], "-usecode") == 0)
    {
      if (codeFileName != NULL)
      {
        HLTWarning("Custom trigger code file was already specified."
                   " Will replace previous value given by -usecode."
        );
      }
      if (argc <= i+1)
      {
        HLTError("The custom trigger code filename was not specified." );
        return -EINVAL;
      }
      codeFileName = argv[i+1];
      if (argc <= i+2)
      {
        HLTError("The custom trigger class name was not specified." );
        return -EINVAL;
      }
      classname = argv[i+2];
      i += 2;
      continue;
    }
    
    HLTError("Unknown option '%s'.", argv[i]);
    return -EINVAL;
  } // for loop
  
  const AliHLTTriggerMenu* menu = NULL;
  if (configFileName != NULL)
  {
    TString cmd = ".x ";
    cmd += configFileName;
    gROOT->ProcessLine(cmd);
    menu = AliHLTGlobalTriggerConfig::Menu();
  }
  
  // Try load the trigger menu from the CDB if it is not loaded yet with the
  // -config option
  if (menu == NULL and AliCDBManager::Instance() != NULL)
  {
    AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
    if (store == NULL)
    {
      HLTError("Could not get the the default storage for the CDB.");
      return -EIO;
    }
    Int_t version = store->GetLatestVersion(fgkTriggerMenuCDBPath, GetRunNo());
    Int_t subVersion = store->GetLatestSubVersion(fgkTriggerMenuCDBPath, GetRunNo(), version);
    AliCDBEntry* entry = AliCDBManager::Instance()->Get(fgkTriggerMenuCDBPath, GetRunNo(), version, subVersion);
    if (entry == NULL)
    {
      HLTError("Could not get the CDB entry for \"%s\".", fgkTriggerMenuCDBPath);
      return -EIO;
    }
    TObject* obj = entry->GetObject();
    if (obj == NULL)
    {
      HLTError("Configuration object for \"%s\" is missing.", fgkTriggerMenuCDBPath);
      return -ENOENT;
    }
    if (obj->IsA() != AliHLTTriggerMenu::Class())
    {
      HLTError("Wrong type for configuration object in \"%s\". Found a %s but we expect a AliHLTTriggerMenu.",
               fgkTriggerMenuCDBPath, obj->ClassName()
      );
      return -EPROTO;
    }
    menu = dynamic_cast<AliHLTTriggerMenu*>(obj);
  }
  
  if (menu == NULL)
  {
    HLTError("No trigger menu configuration found or specified.");
    return -ENOENT;
  }
  
  int result = 0;
  if (codeFileName == NULL)
  {
    HLTDebug("Generating custom HLT trigger class.");
    result = GenerateTrigger(menu, classname, includePaths, includeFiles);
  }
  else
  {
    HLTDebug("Loading HLT trigger class from file '%s'.", codeFileName);
    result = LoadTriggerClass(codeFileName, includePaths);
  }
  if (result != 0) return result;
  
  fTrigger = AliHLTGlobalTrigger::CreateNew(classname.Data());
  if (fTrigger == NULL)
  {
    HLTError("Could not create a new instance of '%s'.", classname.Data());
    return -EIO;
  }
  
  fTrigger->FillFromMenu(*menu);
  fTrigger->ResetCounters(menu->NumberOfItems());
  
  // Set the default values from the trigger menu.
  SetDescription(menu->DefaultDescription());
  SetTriggerDomain(menu->DefaultTriggerDomain());
  
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

  if (!fCodeFileName.IsNull() && gSystem->AccessPathName(fCodeFileName)==0) {
    TString command="rm "; command+=fCodeFileName;
    gSystem->Exec(command);
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
    fTrigger->Add(obj, GetDataType(), GetSpecification());
    obj = GetNextInputObject();
  }

  // Calculate the global trigger result and trigger domain, then create and push
  // back the new global trigger decision object.
  TString description;
  AliHLTTriggerDomain triggerDomain;
  bool triggerResult = fTrigger->CalculateTriggerDecision(triggerDomain, description);
  
  AliHLTGlobalTriggerDecision decision(
      triggerResult,
      // The following will cause the decision to be generated with default values
      // (set in fTriggerDomain and fDescription) if the trigger result is false.
      (triggerResult == true) ? triggerDomain : GetTriggerDomain(),
      (triggerResult == true) ? description.Data() : GetDescription()
    );
  decision.SetCounters(fTrigger->Counters());
  
  // Add the input objects used to the global decision.
  obj = GetFirstInputObject();
  while (obj != NULL)
  {
    if (obj->IsA() == AliHLTTriggerDecision::Class())
    {
      decision.AddTriggerInput( *static_cast<const AliHLTTriggerDecision*>(obj) );
    }
    else
    {
      decision.AddInputObject(obj);
    }
    obj = GetNextInputObject();
  }
  
  TriggerEvent(&decision);
  return 0;
}


int AliHLTGlobalTriggerComponent::GenerateTrigger(
    const AliHLTTriggerMenu* menu, TString& name,
    const TClonesArray& includePaths, const TClonesArray& includeFiles
  )
{
  // Generates the global trigger class that will implement the specified trigger menu.
  // See header for more details.
  
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
  fCodeFileName = name + ".cxx";
  
  // Open a text file to write the code and generate the new class.
  fstream code(fCodeFileName.Data(), ios_base::out | ios_base::trunc);
  if (not code.good())
  {
    HLTError("Could not open file '%s' for writing.", fCodeFileName.Data());
    return -EIO;
  }
  
  TClonesArray symbols(AliHLTTriggerMenuSymbol::Class());
  int result = BuildSymbolList(menu, symbols);
  if (result != 0) return result;
  
  code << "#include <cstring>" << endl;
  code << "#include \"TString.h\"" << endl;
  code << "#include \"TClonesArray.h\"" << endl;
  code << "#include \"AliHLTGlobalTrigger.h\"" << endl;
  code << "#include \"AliHLTGlobalTriggerDecision.h\"" << endl;
  code << "#include \"AliHLTDomainEntry.h\"" << endl;
  code << "#include \"AliHLTTriggerDomain.h\"" << endl;
  code << "#include \"AliHLTReadoutList.h\"" << endl;
  code << "#include \"AliHLTTriggerMenu.h\"" << endl;
  code << "#include \"AliHLTTriggerMenuItem.h\"" << endl;
  code << "#include \"AliHLTTriggerMenuSymbol.h\"" << endl;
  
  // Add any include files that were specified on the command line.
  for (Int_t i = 0; i < includeFiles.GetEntriesFast(); i++)
  {
    TString file = static_cast<const TObjString*>(includeFiles.UncheckedAt(i))->String();
    code << "#include \"" << file.Data() << "\"" << endl;
  }
  
  code << "class " << name << " : public AliHLTGlobalTrigger" << endl;
  code << "{" << endl;
  code << "public:" << endl;
  
  code << "  " << name << "() : AliHLTGlobalTrigger()";
  // Write the symbols in the trigger menu in the initialisation list.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    code << "," << endl;
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "    " << symbol->Name() << "()," << endl;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "    " << symbol->Name() << "TriggerDomain()," << endl;
    }
    code << "    " << symbol->Name() << "DomainEntry(kAliHLTAnyDataType)";
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "," << endl << "    fMenuItemDescription" << i << "()";
  }
  code << endl << "  {" << endl;
  if (fDebugMode)
  {
    code << "    SetLocalLoggingLevel(kHLTLogAll);" << endl;
    code << "    HLTInfo(\"Creating new instance at %p.\", this);" << endl;
  }
  code << "  }" << endl;
  
  code << "  virtual ~" << name << "() {" << endl;
  if (fDebugMode)
  {
    code << "    HLTInfo(\"Deleting instance at %p.\", this);" << endl;
  }
  code << "  }" << endl;
  
  // Generate the FillFromMenu method.
  code << "  virtual void FillFromMenu(const AliHLTTriggerMenu& menu) {" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Filling description entries from trigger menu for global trigger %p.\", this);" << endl;
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "    fMenuItemDescription" << i << " = (menu.Item(" << i
         << ") != NULL) ? menu.Item(" << i << ")->Description() : \"\";" << endl;
  }
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Finished filling description entries from trigger menu.\");" << endl;
    code << "    HLTDebug(\"Filling domain entries from trigger menu symbols for global trigger %p.\", this);" << endl;
  }
  code << "    for (Int_t i = 0; i < menu.SymbolArray().GetEntriesFast(); i++) {" << endl;
  code << "      const AliHLTTriggerMenuSymbol* symbol = dynamic_cast<const"
           " AliHLTTriggerMenuSymbol*>(menu.SymbolArray().UncheckedAt(i));" << endl;
  code << "      if (symbol == NULL) continue;" << endl;
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "      if (strcmp(symbol->Name(), \"" << symbol->Name() << "\") == 0) {" << endl;
    if (fDebugMode)
    {
      code << "        HLTDebug(\"Assinging domain entry value to match for symbol '%s' to '%s'.\","
              " symbol->Name(), symbol->BlockType().AsString().Data());" << endl;
    }
    code << "        " << symbol->Name() << "DomainEntry = symbol->BlockType();" << endl;
    code << "        continue;" << endl;
    code << "      }" << endl;
  }
  code << "    }" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Finished filling domain entries from trigger menu symbols.\");" << endl;
  }
  code << "  }" << endl;
  
  // Generate the NewEvent method.
  code << "  virtual void NewEvent() {" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(\"New event for global trigger object %p, initialising variables to default values.\", this);" << endl;
  }
  // Write code to initialise the symbols in the trigger menu to their default values.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "    " << symbol->Name() << " = " << symbol->DefaultValue() << ";" << endl;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "    " << symbol->Name() << "TriggerDomain.Clear();" << endl;
    }
  }
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Finished initialising variables.\");" << endl;
  }
  code << "  }" << endl;
  
  // Generate the Add method.
  code << "  virtual void Add(const TObject* _object_, const AliHLTComponentDataType& _type_, AliHLTUInt32_t _spec_) {" << endl;
  code << "    AliHLTDomainEntry _type_spec_(_type_, _spec_);" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Adding TObject %p, with class name '%s' from data block"
            " '%s', to global trigger object %p\", _object_, _object_->ClassName(),"
            " _type_spec_.AsString().Data(), this);" << endl;
    code << "    _object_->Print();" << endl;
    code << "    bool _object_assigned_ = false;" << endl;
  }
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    // Write code to check if the block type, specification and class name is correct.
    // Then write code to assign the variable from the input object.
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    TString expr = symbol->AssignExpression();
    if (expr == "") continue; // Skip entries that have no assignment expression.
    bool isTrigDecision = strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0;
    if (fDebugMode)
    {
      if (isTrigDecision)
      {
        code << "    HLTDebug(\"Trying to match input object to class '"
             << symbol->ObjectClass() << "', trigger name '" << symbol->Name()
             << "' and block type '%s'\", " << symbol->Name()
             << "DomainEntry.AsString().Data());" << endl;
      }
      else
      {
        code << "    HLTDebug(\"Trying to match input object to class '"
             << symbol->ObjectClass() << "' and block type '%s'\", "
             << symbol->Name() << "DomainEntry.AsString().Data());" << endl;
      }
    }
    code << "    const " << symbol->ObjectClass() << "* " << symbol->Name()
         << "_object_ = dynamic_cast<const " << symbol->ObjectClass()
         << "*>(_object_);" << endl;
    code << "    if (" << symbol->Name() << "_object_ != NULL and ";
    if (isTrigDecision)
    {
      code << "strcmp(" << symbol->Name() << "_object_->Name(), \""
           << symbol->Name() << "\") == 0 and ";
    }
    code << symbol->Name() << "DomainEntry == _type_spec_) {" << endl;
    TString fullname = symbol->Name();
    fullname += "_object_";
    expr.ReplaceAll("this", fullname);
    code << "      this->" << symbol->Name() << " = " << expr.Data() << ";" << endl;
    if (isTrigDecision)
    {
      code << "      this->" << symbol->Name() << "TriggerDomain = "
           << fullname.Data() << "->TriggerDomain();" << endl;
    }
    if (fDebugMode)
    {
      code << "      HLTDebug(\"Added TObject %p with class name '%s' to variable "
           << symbol->Name() << "\", _object_, _object_->ClassName());" << endl;
      code << "      _object_assigned_ = true;" << endl;
    }
    code << "    }" << endl;
  }
  if (fDebugMode)
  {
    code << "    if (not _object_assigned_) HLTDebug(\"Did not assign TObject %p"
            " with class name '%s' to any variable.\", _object_, _object_->ClassName());"
         << endl;
  }
  code << "  }" << endl;
  
  // Generate the CalculateTriggerDecision method.
  code << "  virtual bool CalculateTriggerDecision(AliHLTTriggerDomain& _domain_, TString& _description_) {" << endl;
  if (fDebugMode)
  {
    code << "    HLTDebug(\"Calculating global HLT trigger result with trigger object at %p.\", this);" << endl;
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    const AliHLTTriggerMenuItem* item = menu->Item(i);
    TString mergeExpr = item->MergeExpression();
    for (Int_t j = 0; j < symbols.GetEntriesFast(); j++)
    {
      AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(j) );
      if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") != 0) continue;
      TString newname = symbol->Name();
      newname += "TriggerDomain";
      mergeExpr.ReplaceAll(symbol->Name(), newname);
    }
    if (fDebugMode)
    {
      code << "    HLTDebug(\"Trying trigger condition " << i
           << " (Description = '%s').\", fMenuItemDescription" << i << ".Data());"
           << endl;
    }
    code << "    if (" << item->TriggerCondision() << ") {" << endl;
    code << "      IncrementCounter(" << i << ");" << endl;
    const char* indentation = "";
    if (item->PreScalar() != 0)
    {
      indentation = "  ";
      code << "      if ((GetCounter(" << i << ") % " << item->PreScalar() << ") == 1) {" << endl;
    }
    code << indentation << "      _domain_ = " << mergeExpr.Data() << ";" << endl;
    code << indentation << "      _description_ = fMenuItemDescription" << i << ";" << endl;
    if (fDebugMode)
    {
      code << indentation << "      HLTDebug(\"Matched trigger condition " << i
           << " (Description = '%s').\", fMenuItemDescription" << i << ".Data());" << endl;
    }
    code << indentation << "      return true;" << endl;
    if (item->PreScalar() != 0)
    {
      code << "      }" << endl;
    }
    code << "    }" << endl;
  }
  code << "    return false;" << endl;
  code << "  }" << endl;
  
  // Generate the custom Factory class.
  code << "  class FactoryImpl : public AliHLTGlobalTrigger::Factory" << endl;
  code << "  {" << endl;
  code << "  public:" << endl;
  code << "    virtual const char* ClassName() const {" << endl;
  code << "      return \"" << name << "\";" << endl;
  code << "    }" << endl;
  code << "    virtual AliHLTGlobalTrigger* New() const {" << endl;
  code << "      return new " << name << "();" << endl;
  code << "    }" << endl;
  code << "  };" << endl;
  
  code << "private:" << endl;
  // Add the symbols in the trigger menu to the list of private variables.
  for (Int_t i = 0; i < symbols.GetEntriesFast(); i++)
  {
    AliHLTTriggerMenuSymbol* symbol = static_cast<AliHLTTriggerMenuSymbol*>( symbols.UncheckedAt(i) );
    code << "  " << symbol->Type() << " " << symbol->Name() << ";" << endl;
    if (strcmp(symbol->ObjectClass(), "AliHLTTriggerDecision") == 0)
    {
      code << "  AliHLTTriggerDomain " << symbol->Name() << "TriggerDomain;" << endl;
    }
    code << "  AliHLTDomainEntry " << symbol->Name() << "DomainEntry;" << endl;
  }
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    code << "  TString fMenuItemDescription" << i << ";" << endl;
  }
  code << "};" << endl;
  
  // Write a global object for the Factory class for automatic registration.
  code << "namespace {" << endl;
  code << "  const " << name << "::FactoryImpl gkFactoryImpl;" << endl;
  code << "};" << endl;
  
  code.close();
  
  // Now we need to compile and load the new class.
  result = LoadTriggerClass(fCodeFileName, includePaths);
  return result;
}


int AliHLTGlobalTriggerComponent::LoadTriggerClass(
    const char* filename, const TClonesArray& includePaths
  )
{
  // Loads the code for a custom global trigger class implementation on the fly.
  
  TString compiler = gSystem->GetBuildCompilerVersion();
  if (compiler.Contains("gcc") or compiler.Contains("icc"))
  {
    TString includePath = "-I${ALICE_ROOT}/include -I${ALICE_ROOT}/HLT/BASE -I${ALICE_ROOT}/HLT/trigger";
    // Add any include paths that were specified on the command line.
    for (Int_t i = 0; i < includePaths.GetEntriesFast(); i++)
    {
      TString path = static_cast<const TObjString*>(includePaths.UncheckedAt(i))->String();
      includePath += " ";
      includePath += path;
    }
    gSystem->SetIncludePath(includePath);
    gSystem->SetFlagsOpt("-O3 -DNDEBUG");
    gSystem->SetFlagsDebug("-g3 -DDEBUG -D__DEBUG");
    
    int result = kTRUE;
    if (fDebugMode)
    {
      result = gSystem->CompileMacro(filename, "g");
    }
    else
    {
      result = gSystem->CompileMacro(filename, "O");
    }
    if (result != kTRUE)
    {
      HLTFatal("Could not compile and load global trigger menu implementation.");
      return -ENOENT;
    }
  }
  else
  {
    // If we do not support the compiler then try interpret the class instead.
    TString cmd = ".L ";
    cmd += filename;
    Int_t errorcode = TInterpreter::kNoError;
    gROOT->ProcessLine(cmd, &errorcode);
    if (errorcode != TInterpreter::kNoError)
    {
      HLTFatal("Could not load interpreted global trigger menu implementation"
               " (Interpreter error code = %d).",
               errorcode
      );
      return -ENOENT;
    }
  }
  
  return 0;
}


int AliHLTGlobalTriggerComponent::FindSymbol(const char* name, const TClonesArray& list)
{
  // Searches for the named symbol in the given list.
  // See header for more details.
  
  for (int i = 0; i < list.GetEntriesFast(); i++)
  {
    const AliHLTTriggerMenuSymbol* symbol = dynamic_cast<const AliHLTTriggerMenuSymbol*>( list.UncheckedAt(i) );
    if (symbol == NULL) continue;
    if (strcmp(symbol->Name(), name) == 0) return i;
  }
  return -1;
}


int AliHLTGlobalTriggerComponent::BuildSymbolList(const AliHLTTriggerMenu* menu, TClonesArray& list)
{
  // Builds the list of symbols to use in the custom global trigger menu
  // implementation class.
  // See header for more details.
  
  for (UInt_t i = 0; i < menu->NumberOfSymbols(); i++)
  {
    const AliHLTTriggerMenuSymbol* symbol = menu->Symbol(i);
    if (FindSymbol(symbol->Name(), list) != -1)
    {
      HLTError("Multiple symbols with the name '%s' defined in the trigger menu.", symbol->Name());
      return -EIO;
    }
    new (list[list.GetEntriesFast()]) AliHLTTriggerMenuSymbol(*symbol);
  }
  
  TRegexp exp("[_a-zA-Z][_a-zA-Z0-9]*");
  TRegexp hexexp("x[a-fA-F0-9]+");
  for (UInt_t i = 0; i < menu->NumberOfItems(); i++)
  {
    const AliHLTTriggerMenuItem* item = menu->Item(i);
    TString str = item->TriggerCondision();
    Ssiz_t start = 0;
    do
    {
      Ssiz_t length = 0;
      Ssiz_t pos = exp.Index(str, &length, start);
      if (pos == kNPOS) break;
      start = pos+length;
      
      // Check if there is a numerical character before the found
      // regular expression. If so, then the symbol is not a valid one
      // and should be skipped.
      if (pos > 0)
      {
        bool notValid = false;
        switch (str[pos-1])
        {
        case '0': case '1': case '2': case '3': case '4':
        case '5': case '6': case '7': case '8': case '9':
          notValid = true;
          break;
        default:
          notValid = false;
          break;
        }
        if (notValid) continue;
      }
      TString s = str(pos, length);
      
      if (s == "and" or s == "and_eq" or s == "bitand" or s == "bitor" or
          s == "compl" or s == "not" or s == "not_eq" or s == "or" or
          s == "or_eq" or s == "xor" or s == "xor_eq" or s == "true" or
          s == "false"
         )
      {
        // Ignore iso646.h and other keywords.
        continue;
      }

      if (FindSymbol(s.Data(), list) == -1)
      {
        AliHLTTriggerMenuSymbol newSymbol;
        newSymbol.Name(s.Data());
        newSymbol.Type("bool");
        newSymbol.ObjectClass("AliHLTTriggerDecision");
        newSymbol.AssignExpression("this->Result()");
        newSymbol.DefaultValue("false");
        new (list[list.GetEntriesFast()]) AliHLTTriggerMenuSymbol(newSymbol);
      }
    }
    while (start < str.Length());
  }
  
  return 0;
}

