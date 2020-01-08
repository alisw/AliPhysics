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
fValidCutsElectron(kTRUE),
fValidCutsNDM(kTRUE),
fValidCutsChargedPion(kTRUE),
fValidCutsChargedKaon(kTRUE),
fValidCutsProton(kTRUE),
fValidCutsDeuteron(kTRUE),
fEventCutArray(0),
fPhotonCutArray(0),
fMesonCutArray(0),
fClusterCutArray(0),
fClusterCutArray2(0),
fMergedClusterCutArray(0),
fElectronCutArray(0),
fNeutralDecayMesonCutArray(0),
fChargedPionCutArray(0),
fChargedKaonCutArray(0),
fProtonCutArray(0),
fDeuteronCutArray(0)
{
  fNCuts                     = 0;
  fNMaxCuts                  = 10;
  fEventCutArray             = new TString[fNMaxCuts];
  fPhotonCutArray            = new TString[fNMaxCuts];
  fMesonCutArray             = new TString[fNMaxCuts];
  fClusterCutArray           = new TString[fNMaxCuts];
  fClusterCutArray2           = new TString[fNMaxCuts];
  fMergedClusterCutArray     = new TString[fNMaxCuts];
  fElectronCutArray          = new TString[fNMaxCuts];
  fNeutralDecayMesonCutArray = new TString[fNMaxCuts];
  fChargedPionCutArray       = new TString[fNMaxCuts];
  fChargedKaonCutArray       = new TString[fNMaxCuts];
  fProtonCutArray            = new TString[fNMaxCuts];
  fDeuteronCutArray          = new TString[fNMaxCuts];


  for(Int_t i=0; i<fNMaxCuts; i++) {
    fEventCutArray[i]             = "";
    fPhotonCutArray[i]            = "";
    fMesonCutArray[i]             = "";
    fClusterCutArray[i]           = "";
    fClusterCutArray2[i]           = "";
    fMergedClusterCutArray[i]     = "";
    fElectronCutArray[i]          = "";
    fNeutralDecayMesonCutArray[i] = "";
    fChargedPionCutArray[i]       = "";
    fChargedKaonCutArray[i]       = "";
    fProtonCutArray[i]            = "";
    fDeuteronCutArray[i]          = "";
  }
}

//________________________________________________________________________
AliCutHandlerPCM::AliCutHandlerPCM(Int_t nMax) :
  fMode(0),
  fNCuts(0),
  fNMaxCuts(1),
  fValidCuts(kTRUE),
  fValidCutsEvent(kTRUE),
  fValidCutsPCM(kTRUE),
  fValidCutsCalo(kTRUE),
  fValidCutsMergedCalo(kTRUE),
  fValidCutsMeson(kTRUE),
  fValidCutsElectron(kTRUE),
  fValidCutsNDM(kTRUE),
  fValidCutsChargedPion(kTRUE),
  fValidCutsChargedKaon(kTRUE),
  fValidCutsProton(kTRUE),
  fValidCutsDeuteron(kTRUE),
  fEventCutArray(0),
  fPhotonCutArray(0),
  fMesonCutArray(0),
  fClusterCutArray(0),
  fClusterCutArray2(0),
  fMergedClusterCutArray(0),
  fElectronCutArray(0),
  fNeutralDecayMesonCutArray(0),
  fChargedPionCutArray(0),
  fChargedKaonCutArray(0),
  fProtonCutArray(0),
  fDeuteronCutArray(0)
{
  fNCuts                     = 0;
  fNMaxCuts                  = nMax;
  fEventCutArray             = new TString[fNMaxCuts];
  fPhotonCutArray            = new TString[fNMaxCuts];
  fMesonCutArray             = new TString[fNMaxCuts];
  fClusterCutArray           = new TString[fNMaxCuts];
  fClusterCutArray2           = new TString[fNMaxCuts];
  fMergedClusterCutArray     = new TString[fNMaxCuts];
  fElectronCutArray          = new TString[fNMaxCuts];
  fNeutralDecayMesonCutArray = new TString[fNMaxCuts];
  fChargedPionCutArray       = new TString[fNMaxCuts];
  fChargedKaonCutArray       = new TString[fNMaxCuts];
  fProtonCutArray            = new TString[fNMaxCuts];
  fDeuteronCutArray          = new TString[fNMaxCuts];


  for(Int_t i=0; i<fNMaxCuts; i++) {
    fEventCutArray[i]             = "";
    fPhotonCutArray[i]            = "";
    fMesonCutArray[i]             = "";
    fClusterCutArray[i]           = "";
    fClusterCutArray2[i]           = "";
    fMergedClusterCutArray[i]     = "";
    fElectronCutArray[i]          = "";
    fNeutralDecayMesonCutArray[i] = "";
    fChargedPionCutArray[i]       = "";
    fChargedKaonCutArray[i]       = "";
    fProtonCutArray[i]            = "";
    fDeuteronCutArray[i]          = "";
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

void AliCutHandlerPCM::AddCutCaloCalo(TString eventCut, TString clusterCut1, TString clusterCut2, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || clusterCut1.Length()!=19 || clusterCut2.Length()!=19 || mesonCut.Length()!=16 ) {
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
  fClusterCutArray[fNCuts]  = clusterCut1;
  fClusterCutArray2[fNCuts]  = clusterCut2;
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


void AliCutHandlerPCM::AddCutPCMDalitz(TString eventCut, TString photonCut, TString electronCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length() !=8 || photonCut.Length() !=26 || electronCut.Length()!=20 || mesonCut.Length()!=16 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent     = kFALSE;
    fValidCutsPCM       = kFALSE;
    fValidCutsElectron  = kFALSE;
    fValidCutsMeson     = kFALSE;
    fValidCuts          = false;
    return;
  }
  fMode                      = 1;
  fEventCutArray[fNCuts]     = eventCut;
  fPhotonCutArray[fNCuts]    = photonCut;
  fMesonCutArray[fNCuts]     = mesonCut;
  fElectronCutArray[fNCuts]  = electronCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutHeavyMesonPCM(TString eventCut, TString photonCut, TString pionCut, TString ndmCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || pionCut.Length()!=9 || ndmCut.Length()!=16 || mesonCut.Length()!=16 ){
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent       = kFALSE;
    fValidCutsPCM         = kFALSE;
    fValidCutsChargedPion = kFALSE;
    fValidCutsNDM         = kFALSE;
    fValidCutsMeson       = kFALSE;
    fValidCuts            = false;
  }
  fMode                              = 0;
  fEventCutArray[fNCuts]             = eventCut;
  fPhotonCutArray[fNCuts]            = photonCut;
  fChargedPionCutArray[fNCuts]       = pionCut;
  fNeutralDecayMesonCutArray[fNCuts] = ndmCut;
  fMesonCutArray[fNCuts]             = mesonCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutHeavyMesonCalo(TString eventCut, TString clusterCut, TString pionCut, TString ndmCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || clusterCut.Length()!=19 || pionCut.Length()!=9 || ndmCut.Length()!=16 || mesonCut.Length()!=16 ){
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent       = kFALSE;
    fValidCutsCalo        = kFALSE;
    fValidCutsChargedPion = kFALSE;
    fValidCutsNDM         = kFALSE;
    fValidCutsMeson       = kFALSE;
    fValidCuts            = false;
  }
  fMode                              = 0;
  fEventCutArray[fNCuts]             = eventCut;
  fClusterCutArray[fNCuts]           = clusterCut;
  fChargedPionCutArray[fNCuts]       = pionCut;
  fNeutralDecayMesonCutArray[fNCuts] = ndmCut;
  fMesonCutArray[fNCuts]             = mesonCut;
  fNCuts++;
  return;
}
void AliCutHandlerPCM::AddCutHeavyMesonPCMCalo(TString eventCut,TString photonCut, TString clusterCut, TString pionCut, TString ndmCut, TString mesonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 || clusterCut.Length()!=19 || pionCut.Length()!=9 || ndmCut.Length()!=16 || mesonCut.Length()!=16 ){
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent       = kFALSE;
    fValidCutsPCM         = kFALSE;
    fValidCutsCalo        = kFALSE;
    fValidCutsChargedPion = kFALSE;
    fValidCutsNDM         = kFALSE;
    fValidCutsMeson       = kFALSE;
    fValidCuts            = false;
  }
  fMode                              = 0;
  fEventCutArray[fNCuts]             = eventCut;
  fPhotonCutArray[fNCuts]            = photonCut;
  fClusterCutArray[fNCuts]           = clusterCut;
  fChargedPionCutArray[fNCuts]       = pionCut;
  fNeutralDecayMesonCutArray[fNCuts] = ndmCut;
  fMesonCutArray[fNCuts]             = mesonCut;
  fNCuts++;
  return;
}

void AliCutHandlerPCM::AddCutPCMMaterial(TString eventCut, TString photonCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || photonCut.Length()!=26 ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsPCM   = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 0;
  fEventCutArray[fNCuts]    = eventCut;
  fPhotonCutArray[fNCuts]   = photonCut;
  fNCuts++;
  return;
}
void AliCutHandlerPCM::AddCutTrackQA(TString eventCut, TString pionCut, TString kaonCut, TString protonCut, TString deuteronCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || pionCut.Length()!=10 || kaonCut.Length()!=10  || protonCut.Length()!=10  || deuteronCut.Length()!=10    ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsChargedPion = kFALSE;
    fValidCutsChargedKaon = kFALSE;
    fValidCutsProton = kFALSE;
    fValidCutsDeuteron = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 0;
  fEventCutArray[fNCuts]    = eventCut;
  fChargedPionCutArray[fNCuts] = pionCut;
  fChargedKaonCutArray[fNCuts] = kaonCut;
  fProtonCutArray[fNCuts] = protonCut;
  fDeuteronCutArray[fNCuts] = deuteronCut;
  fNCuts++;
  return;
}
void AliCutHandlerPCM::AddCutTrackQAPion(TString eventCut, TString pionCut){
  if(fNCuts>=fNMaxCuts) {
    cout << "ERROR in AliCutHandlerPCM: Exceeded maximum number of cuts!" << endl;
    fValidCuts = false;
    return;
  }
  if( eventCut.Length()!=8 || pionCut.Length()!=10    ) {
    cout << "ERROR in AliCutHandlerPCM: Incorrect length of cut string!" << endl;
    fValidCutsEvent = kFALSE;
    fValidCutsChargedPion   = kFALSE;
    fValidCuts      = false;
    return;
  }
  fMode                     = 0;
  fEventCutArray[fNCuts]    = eventCut;
  fChargedPionCutArray[fNCuts] = pionCut;
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
TString AliCutHandlerPCM::GetClusterCut2(Int_t i){
  if(fValidCutsCalo&&i<fNMaxCuts&&i>=0)
    return fClusterCutArray2[i];
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

TString AliCutHandlerPCM::GetElectronCut(Int_t i){
  if(fValidCutsElectron&&i<fNMaxCuts&&i>=0)
    return fElectronCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetElectronCut wrong index i "<<endl;
    return "";
  }
}

TString AliCutHandlerPCM::GetNDMCut(Int_t i){ // Neutral Decay Meson Cut
  if(fValidCutsNDM&&i<fNMaxCuts&&i>=0)
    return fNeutralDecayMesonCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetNDMCut wrong index i "<<endl;
    return "";
  }
}

TString AliCutHandlerPCM::GetPionCut(Int_t i){ // Get charged pion cut
  if(fValidCutsChargedPion&&i<fNMaxCuts&&i>=0)
    return fChargedPionCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetPionCut wrong index i "<<endl;
    return "";
  }
}
TString AliCutHandlerPCM::GetKaonCut(Int_t i){ // Get charged kaon cut
  if(fValidCutsChargedKaon&&i<fNMaxCuts&&i>=0)
    return fChargedKaonCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetKaonCut wrong index i "<<endl;
    return "";
  }
}
TString AliCutHandlerPCM::GetProtonCut(Int_t i){ // Get  proton cut
  if(fValidCutsProton&&i<fNMaxCuts&&i>=0)
    return fProtonCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetProtonCut wrong index i "<<endl;
    return "";
  }
}
TString AliCutHandlerPCM::GetDeuteronCut(Int_t i){ // Get deuteron cut
  if(fValidCutsDeuteron&&i<fNMaxCuts&&i>=0)
    return fDeuteronCutArray[i];
  else {
    cout<<" ERROR in AliCutHandlerPCM: GetDeuteronCut wrong index i "<<endl;
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


TString AliCutHandlerPCM::GetSpecialSettingFromAddConfig (
  TString additionalTrainConfig   = "",
  TString configString            = "",
  TString fileNameMatBudWeights   = "",
  TString addTaskName             = "AddTask_GammaCalo_pp"
){

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
        cout << "INFO: "<< addTaskName.Data() << " activating 'EPCLUSTree'" << endl;
        return "1";
      }else if(tempStr.CompareTo("INVMASSCLUSTree") == 0&& !configString.CompareTo("INVMASSCLUSTree")){
        cout << "INFO: "<< addTaskName.Data() << " activating 'INVMASSCLUSTree'" << endl;
        return "1";
      }else if(tempStr.BeginsWith("MODIFYACC")&& !configString.CompareTo("MODIFYACC")){
        cout << "INFO: "<< addTaskName.Data() << " activating 'MODIFYACC'" << endl;
        return tempStr;
      }else if(tempStr.BeginsWith("LOCALDEBUGFLAG")&& !configString.CompareTo("LOCALDEBUGFLAG")){
        cout << "INFO: "<< addTaskName.Data() << " activating 'LOCALDEBUGFLAG'" << endl;
        TString tempType = tempStr;
        tempType.Replace(0,14,"");
        return tempType;
      }
    }
  }
  return "";
}
