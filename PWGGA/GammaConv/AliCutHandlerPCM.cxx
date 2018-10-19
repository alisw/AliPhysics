#include "AliCutHandlerPCM.h"

using namespace std;

/// \cond CLASSIMP
ClassImp(AliCutHandlerPCM)
/// \endcond

//________________________________________________________________________
AliCutHandlerPCM::AliCutHandlerPCM() :
fMode(0),
fNCuts(0),
fNMaxCuts(10),
fValidCuts(kTRUE),
fValidCutsEvent(kTRUE),
fValidCutsPCM(kTRUE),
fValidCutsCalo(kTRUE),
fValidCutsMergedCalo(kTRUE),
fValidCutsMeson(kTRUE),
fEventCutArray(0),
fPhotonCutArray(0),
fMesonCutArray(0),
fClusterCutArray(0),
fMergedClusterCutArray(0)
{
  fNCuts                  = 0;
  fNMaxCuts               = 10;
  fEventCutArray          = new TString[fNMaxCuts];
  fPhotonCutArray         = new TString[fNMaxCuts];
  fMesonCutArray          = new TString[fNMaxCuts];
  fClusterCutArray        = new TString[fNMaxCuts];
  fMergedClusterCutArray  = new TString[fNMaxCuts];
  for(Int_t i=0; i<fNMaxCuts; i++) {
    fEventCutArray[i]           = "";
    fPhotonCutArray[i]          = "";
    fMesonCutArray[i]           = "";
    fClusterCutArray[i]         = "";
    fMergedClusterCutArray[i]   = "";
  }
}

//________________________________________________________________________
AliCutHandlerPCM::AliCutHandlerPCM(Int_t nMax=10) :
  fMode(0),
  fNCuts(0),
  fNMaxCuts(10),
  fValidCuts(kTRUE),
  fValidCutsEvent(kTRUE),
  fValidCutsPCM(kTRUE),
  fValidCutsCalo(kTRUE),
  fValidCutsMergedCalo(kTRUE),
  fValidCutsMeson(kTRUE),
  fEventCutArray(0),
  fPhotonCutArray(0),
  fMesonCutArray(0),
  fClusterCutArray(0),
  fMergedClusterCutArray(0)
{
  fNCuts                  = 0;
  fNMaxCuts               = nMax;
  fEventCutArray          = new TString[fNMaxCuts];
  fPhotonCutArray         = new TString[fNMaxCuts];
  fMesonCutArray          = new TString[fNMaxCuts];
  fClusterCutArray        = new TString[fNMaxCuts];
  fMergedClusterCutArray  = new TString[fNMaxCuts];
  for(Int_t i=0; i<fNMaxCuts; i++) {
    fEventCutArray[i]           = "";
    fPhotonCutArray[i]          = "";
    fMesonCutArray[i]           = "";
    fClusterCutArray[i]         = "";
    fMergedClusterCutArray[i]   = "";
  }
}

void AliCutHandlerPCM::AddCutPCM(TString eventCut, TString photonCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsPCM   = kFALSE;
    fValidCutsMeson = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 0;
  fEventCutArray[fNCuts]    = eventCut;
  fPhotonCutArray[fNCuts]   = photonCut;
  fMesonCutArray[fNCuts]    = mesonCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutPCM(TString eventCut, TString photonCut, TString mesonCut, TString clusterCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || mesonCut.Length()!=16 || clusterCut.Length()!=19 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsPCM   = kFALSE;
    fValidCutsCalo  = kFALSE;
    fValidCutsMeson = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 0;
  fEventCutArray[fNCuts]    = eventCut;
  fPhotonCutArray[fNCuts]   = photonCut;
  fMesonCutArray[fNCuts]    = mesonCut;
  fClusterCutArray[fNCuts]  = clusterCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutCalo(TString eventCut, TString clusterCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsCalo  = kFALSE;
    fValidCutsMeson = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 2;
  fEventCutArray[fNCuts]    = eventCut;
  fMesonCutArray[fNCuts]    = mesonCut;
  fClusterCutArray[fNCuts]  = clusterCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutMergedCalo(TString eventCut, TString clusterCut, TString clusterMergedCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || clusterCut.Length()!=19 || mesonCut.Length()!=16 || clusterMergedCut.Length()!=19 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsCalo  = kFALSE;
    fValidCutsMeson = kFALSE;
    fValidCutsMergedCalo = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                           = 3;
  fEventCutArray[fNCuts]          = eventCut;
  fMesonCutArray[fNCuts]          = mesonCut;
  fClusterCutArray[fNCuts]        = clusterCut;
  fMergedClusterCutArray[fNCuts]  = clusterMergedCut;
  fNCuts++;
  return;
}


void AliCutHandlerPCM::AddCutPCMCalo(TString eventCut, TString photonCut, TString clusterCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || mesonCut.Length()!=16 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsPCM   = kFALSE;
    fValidCutsCalo  = kFALSE;
    fValidCutsMeson = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 1;
  fEventCutArray[fNCuts]    = eventCut;
  fPhotonCutArray[fNCuts]   = photonCut;
  fMesonCutArray[fNCuts]    = mesonCut;
  fClusterCutArray[fNCuts]  = clusterCut;
  fNCuts++;
  return;
}


Int_t AliCutHandlerPCM::GetNCuts(){
  if(fValidCuts) return fNCuts;
  else return 0;
}

TString AliCutHandlerPCM::GetEventCut(Int_t i){
  if(fValidCutsEvent&&i<fNMaxCuts&&i>=0)
      return fEventCutArray[i];
  else{
    cout << "ERROR in AliCutHandlerPCM: GetEventCut wrong index i" << endl;
    return "";
  }
}

TString AliCutHandlerPCM::GetPhotonCut(Int_t i){
  if(fValidCutsPCM&&i<fNMaxCuts&&i>=0)
    return fPhotonCutArray[i];
  else {
    cout << "ERROR in AliCutHandlerPCM: GetPhotonCut wrong index i" << endl;
    return "";
  }
}

TString AliCutHandlerPCM::GetClusterCut(Int_t i){
  if(fValidCutsCalo&&i<fNMaxCuts&&i>=0)
    return fClusterCutArray[i];
  else {
    cout << "ERROR in AliCutHandlerPCM: GetClusterCut wrong index i" << endl;
    return "";
  }
}


TString AliCutHandlerPCM::GetClusterMergedCut(Int_t i){
  if(fValidCutsMergedCalo&&i<fNMaxCuts&&i>=0)
    return fMergedClusterCutArray[i];
  else {
    cout << "ERROR in AliCutHandlerPCM: GetClusterMergedCut wrong index i" << endl;
    return "";
  }
}

TString AliCutHandlerPCM::GetMesonCut(Int_t i){
  if(fValidCutsMeson&&i<fNMaxCuts&&i>=0)
    return fMesonCutArray[i];
  else {
    cout << "ERROR in AliCutHandlerPCM: GetMesonCut wrong index i" << endl;
    return "";
  }
}


TString AliCutHandlerPCM::GetSpecialFileNameFromString (TString fileNameExternalInputs = "", TString configString = ""){
  TObjArray *rfileNameExternalInputs = fileNameExternalInputs.Tokenize(";");
  if(rfileNameExternalInputs->GetEntries()<1){
    cout << "WARNING: Empty string in GetSpecialFileNameFromString during parsing of fileNameExternalInputs '" << fileNameExternalInputs.Data() << "'" << endl;
    return "";
  }
  for(Int_t i = 0; i<rfileNameExternalInputs->GetEntries() ; i++){
    TObjString* temp = (TObjString*) rfileNameExternalInputs->At(i);
    TString tempStr = temp->GetString();
    if(tempStr.BeginsWith(configString.Data())){
      cout << "INFO: Found special file " << tempStr.Data() <<"!" << endl;
      tempStr.Replace(0,5,"");
      return tempStr;
    }
  }
  return "";
}


TString AliCutHandlerPCM::GetSpecialSettingFromAddConfig (TString additionalTrainConfig = "", TString configString = "", TString fileNameMatBudWeights = ""){
  TObjArray *rAddConfigArr = additionalTrainConfig.Tokenize("_");
  if(rAddConfigArr->GetEntries()<1){
    cout << "WARNING: Empty string in GetSpecialSettingFromAddConfig during parsing of additionalTrainConfig String '" << additionalTrainConfig.Data() << "'" << endl;
    return "";
  }
  for(Int_t i = 0; i<rAddConfigArr->GetEntries() ; i++){
    if(!configString.CompareTo("")){
      TObjString* temp = (TObjString*)rAddConfigArr->At(0);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;
      return tempStr;
    } else {
      TObjString* temp = (TObjString*) rAddConfigArr->At(i);
      TString tempStr = temp->GetString();
      cout<< tempStr.Data()<<endl;
      if(tempStr.Contains("MaterialBudgetWeights") && !configString.CompareTo("MaterialBudgetWeights")){
        TObjArray *fileNameMatBudWeightsArr = fileNameMatBudWeights.Tokenize("/");
        if(fileNameMatBudWeightsArr->GetEntries()<1 ){
          cout<<"ERROR: Empty string in GetSpecialSettingFromAddConfig when reading material budget weights file name" << fileNameMatBudWeights.Data()<< "'" << endl;
          return "";
        }
        TObjString * oldMatObjStr = (TObjString*)fileNameMatBudWeightsArr->At( fileNameMatBudWeightsArr->GetEntries()-1);
        TString  oldfileName  = oldMatObjStr->GetString();
        TString  newFileName  = Form("MCInputFile%s.root",tempStr.Data());
        cout<<newFileName.Data()<<endl;
        if( oldfileName.EqualTo(newFileName.Data()) == 0 ){
          fileNameMatBudWeights.ReplaceAll(oldfileName.Data(),newFileName.Data());
          cout << "INFO: GetSpecialSettingFromAddConfig the material budget weights file has been changed to " <<fileNameMatBudWeights.Data()<<"'"<< endl;
          return fileNameMatBudWeights;
        }
      } else if(tempStr.BeginsWith("CF") && !configString.CompareTo("CF")){
        cout << "INFO: GetSpecialSettingFromAddConfig will use custom branch from Correction Framework!" << endl;
        tempStr.Replace(0,2,"");
        return tempStr;
      } else if(tempStr.BeginsWith("TM") && !configString.CompareTo("TM")){
        tempStr.Replace(0,2,"");
        cout << Form("INFO: GetSpecialSettingFromAddConfig will use running mode '%i' for the TrackMatcher!",tempStr.Atoi()) << endl;
        return tempStr;
      }else if(tempStr.CompareTo("EPCLUSTree") == 0&& !configString.CompareTo("EPCLUSTree")){
        cout << "INFO: AddTask_GammaCalo_pp activating 'EPCLUSTree'" << endl;
        return "1";
      }else if(tempStr.CompareTo("INVMASSCLUSTree") == 0&& !configString.CompareTo("INVMASSCLUSTree")){
        cout << "INFO: AddTask_GammaCalo_pp activating 'INVMASSCLUSTree'" << endl;
        return "1";
      }else if(tempStr.BeginsWith("MODIFYACC")&& !configString.CompareTo("MODIFYACC")){
        return tempStr;
      }else if(tempStr.BeginsWith("LOCALDEBUGFLAG")&& !configString.CompareTo("LOCALDEBUGFLAG")){
        cout << "INFO: AddTask_GammaCalo_pp activating 'LOCALDEBUGFLAG'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,14,"");
        return tempType;
      }
    }
  }
  return "";
}
