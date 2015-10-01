#include "AliHLTTriggerCTPRE.h"
#include "TPRegexp.h"
#include "AliHLTCTPData.h"
#include "TObjString.h"

ClassImp(AliHLTTriggerCTPRE)

//_____________________________________________________________________________
AliHLTTriggerCTPRE::AliHLTTriggerCTPRE()
  : AliHLTTrigger()
  , fName()
  , fGlob()
  , fRegexp(NULL)
{
  //ctor
}

const char* AliHLTTriggerCTPRE::fgkDefaultOCDBEntry="HLT/ConfigHLT/CTPREtrigger";

//_____________________________________________________________________________
AliHLTTriggerCTPRE::~AliHLTTriggerCTPRE()
{
  //dtor
  delete fRegexp;
}

//_____________________________________________________________________________
const char* AliHLTTriggerCTPRE::GetTriggerName() const
{
  if (!fName.IsNull())
    return fName.Data();
  else
    return "CTPREtrigger";
}

//_____________________________________________________________________________
AliHLTComponent* AliHLTTriggerCTPRE::Spawn()
{
  return new AliHLTTriggerCTPRE;
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::DoTrigger()
{
  int retCode = 0;
  if (!IsDataEvent()) {
    IgnoreEvent();  // dont generate any trigger decision.
  }

  Bool_t decision = kFALSE;
  const AliHLTCTPData* ctpData = CTPData();
  
  //take the current trigger, zero out anything that should be zero
  AliHLTTriggerMask_t currentTrigger = ctpData->Triggers();
  currentTrigger &= ctpData->Mask();

  //check every fired trigger agains the expressions
  for (int i=0; i<NCTPTRIGGERCLASSES; i++)
  {
    if (!currentTrigger.test(i)) continue;

    const char* triggerName = ctpData->Name(i);
    
    if (!fGlob.IsNull()) 
    {
      if (Globncmp(triggerName, fGlob.Data(), strnlen(triggerName,100), fGlob.Length()))
      {
        decision = kTRUE;
        break;
      }
    }
    if (!(fRegexp->GetPattern()).IsNull())
    {
      if ((fRegexp->Match(triggerName)>0))
      {
        decision = kTRUE;
        break;
      }
    }
  }
  
  TriggerEvent(decision);

  return retCode;
}

//_____________________________________________________________________________
Bool_t AliHLTTriggerCTPRE::Globncmp(const char* triggerName, const char* glob, int triggerNameSize, int globSize )
{
  if (globSize == 0) return kFALSE;
  for (int i=0; i<((triggerNameSize<globSize)?triggerNameSize:globSize); i++)
  {
    if (!(glob[i]=='*' || triggerName[i]==glob[i])) {return kFALSE;}
  }
  return kTRUE;
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::DoInit(int argc, const char** argv)
{
  int retCode=0;

  // check if the -triggername argument is used
  // the name of the trigger determines the following initialization
  vector<const char*> remainingArgs;
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i], "-triggername")==0) {
      if (++i<argc) fName=argv[i];
      else {
	HLTError("invalid parameter for argument '-triggername', string expected");
	return -EINVAL;
      }
      continue;
    }
    remainingArgs.push_back(argv[i]);
  }

  //Init the CTP data
  if (SetupCTPData() == -ENOMEM) 
  {
    HLTError("could not SetupCTPData(); ENOMEM");
    return -ENOMEM;
  }

  // get path from triggername, use default object otherwise
  TString cdbPath;
  if (!fName.IsNull()) {
    cdbPath="HLT/ConfigHLT/";
    cdbPath+=fName;
  } else {
    cdbPath=fgkDefaultOCDBEntry;
  }

  // -- Check if CDB object is AliHLTESDTrackCuts or TObjString 
  //    and configure from it. Replace "-" by "_._" if needed in the cdbPath
  retCode = ConfigureFromCDBObject(cdbPath);

  // -- Configure from the command line parameters if specified
  if (retCode>=0 && argc>0)
    retCode=ConfigureFromArgumentString(remainingArgs.size(), &(remainingArgs[0]));

  //make sure the regexp pointer is not null, even if the regexp is empty
  if (!fRegexp) fRegexp=new TPRegexp();

  return retCode;
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::DoDeinit()
{
  // see header file for class documentation

  return 0;
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::Reconfigure(const char* cdbEntry, const char* /*chainId*/)
{
  // configure from the specified antry or the default one
  TString cdbPath;
  if (!cdbEntry || cdbEntry[0]==0) {
    if (!fName.IsNull()) {
      cdbPath="HLT/ConfigHLT/";
      cdbPath+=fName;
    } else {
      cdbPath=fgkDefaultOCDBEntry;
    }
  } else {
    cdbPath=cdbEntry;
  }

  return ConfigureFromCDBObject(cdbPath);
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::ReadPreprocessorValues(const char* /*modules*/)
{
  return 0;
}

//_____________________________________________________________________________
Int_t AliHLTTriggerCTPRE::ConfigureFromCDBObject(TString cdbPath)
{
  Int_t retCode = 0;
  TString arguments;

  // -- check for "-" and replace by "_._" in the path name
  cdbPath.ReplaceAll("-",1,"_._",3);

  TObject* pCDBObject = LoadAndExtractOCDBObject(cdbPath);
  if (pCDBObject) {
      TObjString* pString = dynamic_cast<TObjString*>(pCDBObject);
      if (pString) {
	HLTInfo("Received configuration object string: \'%s\'", pString->GetString().Data());
	arguments+=pString->GetString().Data();
      } 
      else {
	HLTError("Configuration object \"%s\" has wrong type, required AliHLTESDTrackCuts or TObjString", cdbPath.Data());
	retCode=-EINVAL;
      }
    }
  else {
    HLTError("Can not fetch object \"%s\" from CDB", cdbPath.Data());
    retCode=-ENOENT;
  }
  
  if ( retCode>=0 && !arguments.IsNull() ) {
    const Char_t* array = arguments.Data();
    retCode = ConfigureFromArgumentString(1, &array);
  }

  return retCode;
}

//_____________________________________________________________________________
int AliHLTTriggerCTPRE::ScanConfigurationArgument(int argc, const char** argv)
{
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -glob
  // selection glob (fast)
  if (argument.CompareTo("-glob")==0) {
    if (++i >= argc) return -EPROTO;
    argument = argv[i];
    fGlob = argument;
    return 2;
  }

  if (argument.CompareTo("-regex")==0) {
    if (++i >= argc) return -EPROTO;
    argument = argv[i];
    delete fRegexp;
    fRegexp = new TPRegexp(argument);
    if (!fRegexp->IsValid()) {
      HLTError("regexp %s not valid!\n", fRegexp->GetPattern().Data());
      delete fRegexp;
      fRegexp = new TPRegexp();
    }
    return 2;
  }

  // unknown argument
  return -EINVAL;
}
