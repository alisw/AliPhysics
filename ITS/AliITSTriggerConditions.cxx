////////////////////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                                         //
//                                                                                //
// Implementation of conditions data from Pixel Trigger (PIT)                     //
//                                                                                //
// The information is propagated from pixel trigger system to DCS file exchange   //
// server (text file format). The ReadFromTextFile method will populate this      //
// object with the values from the text file. Via a Preprocessor, this object     //
// can be stored in OCDB.                                                         //
//                                                                                //
// One can also manually create conditions data that may be interesting for       //
// simulation.                                                                    //
//                                                                                //
////////////////////////////////////////////////////////////////////////////////////

#include "AliITSTriggerConditions.h"
#include "AliITSTriggerAlgorithmConditions.h"
#include <TError.h>
#include <fstream>

ClassImp(AliITSTriggerConditions)

//__________________________________________________________________________________
AliITSTriggerConditions::AliITSTriggerConditions() :
TObject(),
fRunNumber(0),
fFirmWareVersion(0),
fGlobalDescription("n/a"),
fVersionRegister(0),
fInputConditionsVersion(0),
fParametersVersion(0),
fInActiveChips(1200),
fNumAlgo(0),
fAlgoList(TObjArray(10))
{
  // default constructor
  fAlgoList.SetOwner(kTRUE);
}
//__________________________________________________________________________________
AliITSTriggerConditions::AliITSTriggerConditions(const AliITSTriggerConditions& cond) :
TObject(),
fRunNumber(cond.fRunNumber),
fFirmWareVersion(cond.fFirmWareVersion),
fGlobalDescription(cond.fGlobalDescription),
fVersionRegister(cond.fVersionRegister),
fInputConditionsVersion(cond.fInputConditionsVersion),
fParametersVersion(cond.fParametersVersion),
fInActiveChips(cond.fInActiveChips),
fNumAlgo(cond.fNumAlgo),
fAlgoList(cond.fAlgoList)
{
  // copy constructor
  fAlgoList.SetOwner(kTRUE);
}
//__________________________________________________________________________________
AliITSTriggerConditions::~AliITSTriggerConditions() 
{
  // destructor
  ClearAlgorithms();
}
//______________________________________________________________________
AliITSTriggerConditions& AliITSTriggerConditions::operator=(const AliITSTriggerConditions& cond) {
  // assignment operator
  if (this!=&cond) {
    fRunNumber = cond.fRunNumber;
    fFirmWareVersion = cond.fFirmWareVersion;
    fGlobalDescription = cond.fGlobalDescription;
    fVersionRegister = cond.fVersionRegister;
    fInputConditionsVersion = cond.fInputConditionsVersion;
    fParametersVersion = cond.fParametersVersion;
    fInActiveChips = cond.fInActiveChips;
    fNumAlgo = cond.fNumAlgo;
    fAlgoList = cond.fAlgoList;
  }
  return *this;
}
//__________________________________________________________________________________
void AliITSTriggerConditions::DumpAll() const {
  // Dumps all conditions data

  printf("[Header]\n");
  printf("RUN_NUMBER = %d\n",fRunNumber);
  printf("PROCESSING_FIRMWARE_VERSION = %d\n",fFirmWareVersion);
  printf("GLOBAL_DESCRIPTION = %s\n",fGlobalDescription.Data());
  printf("VERSION_REGISTER_VALUE = %d\n",fVersionRegister);
  printf("INPUT_CONDITIONS_VERSION = %d\n",fInputConditionsVersion);
  printf("PARAMETERS_VERSION = %d\n",fParametersVersion);
  printf("\n");

  printf("[Outputs]\n");
  for (UInt_t i=0; i<fNumAlgo; i++) {
    printf("%d = '%s', '%s'\n", GetAlgoIDI(i), GetAlgoLabelI(i), GetAlgoDescriptionI(i));
  }
  printf("\n");

  printf("[Output_parameters]\n");
  for (UInt_t i=0; i<fNumAlgo; i++) {
    printf("%d =", GetAlgoIDI(i));
    for (Short_t p=0; p<GetNumAlgoParamI(i); p++) {
      printf(" '%s', %d;", GetAlgoParamNameII(i,p), GetAlgoParamValueII(i,p));
    }
    printf("\n");
  }
  printf("\n");

  printf("[Active_chips]\n");
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      UInt_t nActiveOnHs = 0;
      TString inactiveStr = "";
      for (UInt_t chip=0; chip<10; chip++) {
	Bool_t isChipActive = IsChipActive(eq,hs,chip);
	inactiveStr.Append(Form("%d",isChipActive));
	nActiveOnHs+=isChipActive;
      }
      if (nActiveOnHs<10) {
	printf("%d,%c,%d=%s\n", eq%10, eq<10 ? 'A' : 'C', hs, inactiveStr.Data());
      }
    }
  }

}
//__________________________________________________________________________________
void AliITSTriggerConditions::ResetAll() {
  // clear all data, and put default values
  fRunNumber=0;
  fFirmWareVersion=0;
  fGlobalDescription="n/a";
  fVersionRegister=0;
  fInputConditionsVersion=0;
  fParametersVersion=0;
  ResetInActiveChips();
  ClearAlgorithms();
}
//__________________________________________________________________________________
void AliITSTriggerConditions::ClearAlgorithms() {
  // clears the list of algorithms
  fAlgoList.Clear();
  fNumAlgo=0;
}
//__________________________________________________________________________________
void AliITSTriggerConditions::ClearAlgoParamsI(UShort_t aIndex) {
  // clears the list of parameters for algorithm with index aIndex  
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::ClearAlgoParamsI", "index %d out of range", aIndex);
    return;
  }
  ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->ClearParams();
}
//__________________________________________________________________________________
void AliITSTriggerConditions::ClearAlgoParamsL(const Char_t* aLabel) {
  // clears the list of parameters for algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(),aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::ClearAlgoParamsL", "label %s not found", aLabel);
  }
  ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->ClearParams();
}
//__________________________________________________________________________________
void AliITSTriggerConditions::AddAlgo(UShort_t id, const Char_t* aLabel, const Char_t* aDescr) {
  // adds a new algorithm with id 'id', label aLabel, and description aDescr
  // if the id or label is already used in the list of algorithms, the old entry will be over-written
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0 || 
	((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetID() == id ) {
      findIndex = i;
      break;
    }
  }
  if (findIndex<fNumAlgo) {
    delete fAlgoList.At(findIndex);
    fAlgoList.AddAt(new AliITSTriggerAlgorithmConditions(id,aLabel,aDescr),findIndex);
  }
  else {
    fAlgoList.AddAtAndExpand(new AliITSTriggerAlgorithmConditions(id,aLabel,aDescr),fNumAlgo);
    fNumAlgo++;
  }
}
//__________________________________________________________________________________
void AliITSTriggerConditions::AddAlgoParam(UShort_t id, const Char_t* name, Int_t value) {
  // adds a new parameter with name 'name' and value 'value', for the algorithm with id 'id'
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetID() == id) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::AddAlgoParam", "id %d not found", id);
    return;
  }
  ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->AddParam(name, value);
}
//__________________________________________________________________________________
Short_t AliITSTriggerConditions::GetAlgoIndexL(const Char_t* aLabel) const {
  // returns the index for the algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoIndexL", "label %s not found", aLabel);
    return -1;
  }
  return findIndex;
}
//__________________________________________________________________________________
Short_t AliITSTriggerConditions::GetAlgoIDI(UShort_t aIndex) const {
  // returns the ID for the algorithm with index aIndex (in real life, could be 1-10)
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoIDI", "index %d out of range", aIndex);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetID();
}
//__________________________________________________________________________________
const Char_t* AliITSTriggerConditions::GetAlgoLabelI(UShort_t aIndex) const {
  // returns the label for the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoLabelI", "index %d out of range", aIndex);
    return "";
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetLabel();
}
//__________________________________________________________________________________
const Char_t* AliITSTriggerConditions::GetAlgoDescriptionI(UShort_t aIndex) const {
  // returns the description for the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoDescriptionI", "index %d out of range", aIndex);
    return "";
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetDescription();
}
//__________________________________________________________________________________
Short_t AliITSTriggerConditions::GetNumAlgoParamI(UShort_t aIndex) const {
  // returns the number of parameters, corresponding to the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetNumAlgoParamI", "index %d out of range", aIndex);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetNumParam();
}
//__________________________________________________________________________________
const Char_t* AliITSTriggerConditions::GetAlgoParamNameII(UShort_t aIndex, UShort_t pIndex) const {
  // returns the parameter name for the parameter with index pIndex, corresponding to the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamNameII", "index %d out of range", aIndex);
    return "";
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetParamNameI(pIndex);
}
//__________________________________________________________________________________
Int_t AliITSTriggerConditions::GetAlgoParamValueII(UShort_t aIndex, UShort_t pIndex) const {
  // returns the parameter value for the parameter with index pIndex, corresponding to the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamValueII", "index %d out of range", aIndex);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetParamValueI(pIndex);
}
//__________________________________________________________________________________
Int_t AliITSTriggerConditions::GetAlgoParamValueIN(UShort_t aIndex, const Char_t* pName) const {
  // returns parameter value for the parameter named pName, corresponding to the algorithm with index aIndex
  if (aIndex>=fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamValueIN", "index %d out of range", aIndex);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(aIndex))->GetParamValueN(pName);
}
//__________________________________________________________________________________
Short_t AliITSTriggerConditions::GetNumAlgoParamL(const Char_t* aLabel) const {
  // returns the number of parameters, corresponding to the algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::GetNumAlgoParamL", "label %s not found", aLabel);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->GetNumParam();
}
//__________________________________________________________________________________
const Char_t* AliITSTriggerConditions::GetAlgoParamNameLI(const Char_t* aLabel, UShort_t pIndex) const {
  // returns parameter name for the parameter with index pIndex, corresponding to the algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamNameLI", "label %s not found", aLabel);
    return "";
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->GetParamNameI(pIndex);
}
//__________________________________________________________________________________
Int_t AliITSTriggerConditions::GetAlgoParamValueLI(const Char_t* aLabel, UShort_t pIndex) const {
  // returns parameter value for the parameter with index pIndex, corresponding to the algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamValueLI", "label %s not found", aLabel);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->GetParamValueI(pIndex);
}
//__________________________________________________________________________________
Int_t AliITSTriggerConditions::GetAlgoParamValueLN(const Char_t* aLabel, const Char_t* pName) const {
  // returns parameter value for the parameter named pName, corresponding to the algorithm with label aLabel
  UShort_t findIndex=fNumAlgo;
  for (UInt_t i=0; i<fNumAlgo; i++) {
    if (strcmp(((AliITSTriggerAlgorithmConditions*)fAlgoList.At(i))->GetLabel(), aLabel) == 0) {
      findIndex = i;
      break;
    }
  }
  if (findIndex==fNumAlgo) {
    Error("AliITSTriggerConditions::GetAlgoParamValueLN", "label %s not found", aLabel);
    return -1;
  }
  return ((AliITSTriggerAlgorithmConditions*)fAlgoList.At(findIndex))->GetParamValueN(pName);
}
//__________________________________________________________________________________
UInt_t AliITSTriggerConditions::GetChipKey(Int_t eq, Int_t hs, Int_t chip) const {
  // translates eq,hs,chip numbers into one integer key (0-1199)
  if (eq<0 || eq>=20 || hs<0 || hs>=6 || chip<0 || chip>=10) {
    Error("AliITSTriggerConditions::GetChipKey", "eq,hs,chip = %d,%d,%d out of range",eq,hs,chip);
    return 0;
  }
  return eq*60 + hs*10 + chip;
}
//__________________________________________________________________________________
void AliITSTriggerConditions::GetChipFromKey(UInt_t key, Int_t& eq, Int_t& hs, Int_t& chip) const {
  // translates a chip key back into eq,hs,chip numbers
  if (key>=1200) {
    Error("AliITSTriggerConditions::GetChipFromKey", "key = %d out of range", key);
    return;
  }
  eq   = key/60;
  hs   = (key%60)/10;
  chip = key%10;
}
//__________________________________________________________________________________
Bool_t AliITSTriggerConditions::GetNextInActiveChip(Int_t& eq, Int_t& hs, Int_t& chip) const {
  // Returns true if an in-active chip was found (start looking after the bit number
  // corresponding to the input parameters eq,hs,chip).
  // If either of eq,hs,chip < 0 , start from beginning of TBits array.
  // See example of usage in AliITSDetTypeRec::RemoveFastOrFiredInActive.
  UInt_t searchIndex;
  if (eq<0 || hs<0 || chip<0) searchIndex = 0;
  else searchIndex = GetChipKey(eq, hs, chip) + 1;
  UInt_t nextIndex = fInActiveChips.FirstSetBit(searchIndex);
  if (nextIndex==1200) return kFALSE;
  GetChipFromKey(nextIndex, eq, hs, chip);
  return kTRUE;
}
//__________________________________________________________________________________
void AliITSTriggerConditions::DumpInActiveChips() const {
  // Prints a list of all inactive chips
  UInt_t startBit=0;
  UInt_t occ=0;
  while (startBit<1200) {
    startBit = fInActiveChips.FirstSetBit(startBit);
    if (startBit<1200) {
      occ++;
      Int_t eq,hs,chip;
      GetChipFromKey(startBit,eq,hs,chip);
      printf("%3d: %d,%d,%d\n",occ,eq,hs,chip);
      startBit++;
    }
  }
}
//__________________________________________________________________________________
void AliITSTriggerConditions::ReadFromTextFile(const Char_t* fileName) {
  // Reads conditions from text file (file format is used online by PIT system)

  ResetAll();

  const Int_t maxS = 500;
  enum headers {HEAD, OUTPUT, PARAM, ACTIVECHIP};

  ifstream file;
  file.open(fileName, ifstream::in);
  if (file.fail()) {
    Error("AliITSTriggerConditions::ReadFromTextFile","No file (%s) present.",fileName);
    return;
  }

  Int_t headType = -1; // no header read from start
  UInt_t nl = 0;
  Char_t cline[maxS];
  while(!file.eof()) {
    // *** get line
    nl++;
    file.getline(cline,maxS);
    TString line(cline);

    // *** remove comments from line
    Int_t skipPos = line.First('#');
    if (skipPos>=0) {
      line.Remove(skipPos,maxS);
    }

    // *** check what type of information the line has (header or not...)
    Int_t brackPos1 = line.First('[');
    Int_t brackPos2 = line.First(']');
    if (brackPos1==0 && brackPos2-1>brackPos1) {
      // *** parse header line (header has to come first on the line)
      TString headword = line(brackPos1+1,brackPos2-1-brackPos1);
      if      (headword.CompareTo("Header",           TString::kIgnoreCase) == 0) headType = HEAD;
      else if (headword.CompareTo("Outputs",          TString::kIgnoreCase) == 0) headType = OUTPUT;
      else if (headword.CompareTo("Output_parameters",TString::kIgnoreCase) == 0) headType = PARAM;
      else if (headword.CompareTo("Active_chips",     TString::kIgnoreCase) == 0) headType = ACTIVECHIP;
    }
    else {
      // *** parse non-header line

      // HEAD data
      if (headType==HEAD) {
	TString descrWord, valueWord;
	if (! SplitStringIn2(line,descrWord,valueWord,'=')) continue;
	descrWord.ReplaceAll(" ","");
	valueWord.Remove(TString::kBoth,' ');

	if (descrWord.CompareTo("RUN_NUMBER",TString::kIgnoreCase) == 0 && valueWord.IsDigit()) {
	  SetRunNumber(valueWord.Atoi());
	}
	else if (descrWord.CompareTo("PROCESSING_FIRMWARE_VERSION",TString::kIgnoreCase) == 0 && valueWord.IsDigit()) {
	  SetFirmWareVersion(valueWord.Atoi());
	}
	else if (descrWord.CompareTo("GLOBAL_DESCRIPTION",TString::kIgnoreCase) == 0) {
	  SetGlobalDescription(valueWord.Data());
	}
	else if (descrWord.CompareTo("VERSION_REGISTER_VALUE",TString::kIgnoreCase) == 0 && valueWord.IsDigit()) {
	  SetVersionRegister(valueWord.Atoi());
	}
	else if (descrWord.CompareTo("INPUT_CONDITIONS_VERSION",TString::kIgnoreCase) == 0 && valueWord.IsDigit()) {
	  SetInputConditionsVersion(valueWord.Atoi());
	}
	else if (descrWord.CompareTo("PARAMETERS_VERSION",TString::kIgnoreCase) == 0 && valueWord.IsDigit()) {
	  SetParametersVersion(valueWord.Atoi());
	}
      }
      // OUTPUT data
      else if (headType==OUTPUT) {
	TString idWord, labelWord, descrWord, restWord;
	if (! SplitStringIn2(line,idWord,restWord,'=')) continue;
	if (! idWord.IsDigit()) continue;
	if (! SplitStringIn2(restWord,labelWord,descrWord,',')) continue;
	labelWord = GetStringBetween(labelWord,'\'','\'');
	descrWord = GetStringBetween(descrWord,'\'','\'');
	if (labelWord.Length()==0 || descrWord.Length()==0) continue;
	//	printf("id: %d , label '%s' , descr '%s'\n", idWord.Atoi(),labelWord.Data(),descrWord.Data());
	AddAlgo(idWord.Atoi(),labelWord.Data(),descrWord.Data());
      }
      // PARAM data
      else if (headType==PARAM) {
	TString idWord, restWord;
	if (! SplitStringIn2(line,idWord,restWord,'=')) continue;
	if (! idWord.IsDigit()) continue;
	while (restWord.Length()>0) {
	  TString parWord, nameWord, valWord;
	  SplitStringIn2(restWord,parWord,restWord,';');
	  if (! SplitStringIn2(parWord,nameWord,valWord,',')) break;
	  nameWord = GetStringBetween(nameWord,'\'','\'');
	  if (nameWord.Length()==0 || valWord.Length()==0 || ! valWord.IsDigit()) break;
	  //	  printf("id %d , param %s , value %d\n",idWord.Atoi(),nameWord.Data(),valWord.Atoi());
	  AddAlgoParam(idWord.Atoi(),nameWord.Data(),valWord.Atoi());
	}
      }
      // ACTIVECHIP data
      if (headType==ACTIVECHIP) {
	TString eqWord, sideWord, hsWord, chipWord, restWord;
	if (! SplitStringIn2(line,eqWord,restWord,',')) continue;
	if (! eqWord.IsDigit()) continue;
	UInt_t eq = eqWord.Atoi();
	if (eq>=20) continue;
	if (! SplitStringIn2(restWord,sideWord,restWord,',')) continue;
	sideWord.ReplaceAll(" ","");
	if (sideWord.CompareTo("A",TString::kIgnoreCase) == 0) {}
	else if (sideWord.CompareTo("C",TString::kIgnoreCase) == 0) eq+=10;
	else continue;
	if (! SplitStringIn2(restWord,hsWord,chipWord,'=')) continue;
	if (! hsWord.IsDigit()) continue;
	UInt_t hs = hsWord.Atoi();
	if (hs>=6) continue;
	chipWord.ReplaceAll(" ","");
	if (chipWord.Length()!=10) continue;
	for (UInt_t chip=0; chip<10; chip++) {
	  if (chipWord[9-chip]=='0') {
	    //	    printf("Chip %d,%d,%d inactive\n",eq,hs,chip);
	    SetInActiveChip(eq,hs,chip);
	  }
	}
      }

    }
  }
  file.close();
}
//__________________________________________________________________________________
Bool_t AliITSTriggerConditions::SplitStringIn2(TString orig, TString& word1, TString& word2, Char_t sep) {
  // splits a string in two parts (one before separator character and one after)
  Int_t sepPos = orig.First(sep);
  if (sepPos<0) sepPos = orig.Length();
  word1 = orig(0,sepPos);
  word2 = orig(sepPos+1,orig.Length());
  return (word1.Length()>0 && word2.Length()>0);
}
//__________________________________________________________________________________
TString AliITSTriggerConditions::GetStringBetween(TString orig, Char_t sep1, Char_t sep2) {
  // returns string between separator character 1 and separator character 2
  Int_t pos1 = orig.First(sep1);
  if (pos1<0) return "";
  TString ret = orig(pos1+1,orig.Length());
  Int_t pos2 = ret.First(sep2);
  if (pos2<0) return "";
  ret = ret(0,pos2);
  return ret;
}
//__________________________________________________________________________________
Bool_t AliITSTriggerConditions::IsEqualTo(AliITSTriggerConditions *cond) const {
  // checks if this object contains the same information as the input cond object
  if (fRunNumber != cond->GetRunNumber()) return kFALSE;
  if (fFirmWareVersion != cond->GetFirmWareVersion()) return kFALSE;
  if (fGlobalDescription.CompareTo(cond->GetGlobalDescription()) !=0) return kFALSE;
  if (fVersionRegister != cond->GetVersionRegister()) return kFALSE;
  if (fInputConditionsVersion != cond->GetInputConditionsVersion()) return kFALSE;
  if (fParametersVersion != cond->GetParametersVersion()) return kFALSE;

  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	if (IsChipActive(eq,hs,chip) != cond->IsChipActive(eq,hs,chip)) return kFALSE;
      }
    }
  }

  if (fNumAlgo != cond->GetNumAlgo()) return kFALSE;
  for (Short_t alg1=0; alg1<fNumAlgo; alg1++) {
    Short_t alg2 = cond->GetAlgoIndexL(GetAlgoLabelI(alg1));
    if (alg2<0) return kFALSE;
    if (GetAlgoIDI(alg1) != cond->GetAlgoIDI(alg2)) return kFALSE;
    if (strcmp(GetAlgoDescriptionI(alg1), cond->GetAlgoDescriptionI(alg2)) != 0) return kFALSE;
    if (GetNumAlgoParamI(alg1) != cond->GetNumAlgoParamI(alg2)) return kFALSE;
    for (Short_t par1=0; par1<GetNumAlgoParamI(alg1); par1++) {      
      const Char_t* paramName = GetAlgoParamNameII(alg1,par1);
      if (GetAlgoParamValueIN(alg1,paramName) != cond->GetAlgoParamValueIN(alg2,paramName)) return kFALSE;
    }
  }

  return kTRUE;
}
