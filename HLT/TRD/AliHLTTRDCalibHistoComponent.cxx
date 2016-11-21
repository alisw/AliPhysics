// $Id: AliHLTTRDCalibHistoComponent.cxx 40282 2010-04-09 13:29:10Z richterm $

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors:                                                               *
 *          for The ALICE HLT Project.                                    *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

//  @file   AliHLTTRDCalibHistoComponent.cxx
//  @author Theodor Rascanu
//  @date   25.04.2010
//  @brief  A TRDCalibration histogramming component for the HLT. 
// 

#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"
#include "TH2I.h"
#include "TH2.h"
#include "TProfile2D.h"

#include "AliHLTReadoutList.h"

#include "AliHLTTRDCalibHistoComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliCDBEntry.h"
#include "AliRawReaderMemory.h"

#include "AliTRDCalPad.h"
#include "AliTRDCalDet.h"

#include "AliTRDCalibraFillHisto.h"
#include "AliTRDtrackV1.h"

#include "AliTRDCalibraFit.h"
#include "AliTRDCalibraMode.h"
#include "AliTRDCalibraVector.h"
#include "AliTRDCalibraVdriftLinearFit.h"
#include "AliTRDReconstructor.h"
#include "AliTRDrecoParam.h"

#include <cstdlib>
#include <cerrno>
#include <string>

using namespace std;

ClassImp(AliHLTTRDCalibHistoComponent);

AliHLTTRDCalibHistoComponent::AliHLTTRDCalibHistoComponent()
  : AliHLTProcessor(),
    fOutputSize(500000),
    fSpec(0),
    fTracksArray(NULL),
    fOutArray(NULL),
    fTRDCalibraFillHisto(NULL),
    fSavedTimeBins(kFALSE),
    fTrgStrings(NULL),
    fAccRejTrg(0),
    fMinClusters(0),
    fMinTracklets(0),
    fTakeAllEvents(kFALSE)
{
  // Default constructor
}

AliHLTTRDCalibHistoComponent::~AliHLTTRDCalibHistoComponent()
{
  // Destructor
}

const char* AliHLTTRDCalibHistoComponent::GetComponentID()
{
  // Return the component ID const char *
  return "TRDCalibHisto"; // The ID of this component
}

void AliHLTTRDCalibHistoComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // Get the list of input data
  list.clear(); // We do not have any requirements for our input data type(s).
  list.push_back(AliHLTTRDDefinitions::fgkTracksV1DataType);
}

AliHLTComponentDataType AliHLTTRDCalibHistoComponent::GetOutputDataType()
{
  // Get the output data type
  return kAliHLTMultipleDataType;
  //  return AliHLTTRDDefinitions::fgkCalibrationDataType;
 
}

int AliHLTTRDCalibHistoComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // Get the output data type
  tgtList.clear();
  tgtList.push_back(AliHLTTRDDefinitions::fgkCalibrationDataType);
  return tgtList.size();
}

void AliHLTTRDCalibHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // Get the output data size
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDCalibHistoComponent::Spawn()
{
  // Spawn function, return new instance of this class
  return new AliHLTTRDCalibHistoComponent;
};

int AliHLTTRDCalibHistoComponent::DoInit( int argc, const char** argv )
{
  int iResult=0;
  if(fTrgStrings)
    delete fTrgStrings;
  fTrgStrings = new TObjArray();
  
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }

  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }
  if(iResult>=0){
    iResult=SetParams();
  }
  return iResult;
}

int AliHLTTRDCalibHistoComponent::Configure(const char* arguments){
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("output_size")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting output size to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fOutputSize=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-minClusters")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting minCusters to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fMinClusters=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-minTracklets")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting minTracklets to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fMinTracklets=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else if (argument.CompareTo("-TrgStr")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Select TrgStr: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fTrgStrings->Add(new TObjString(((TObjString*)pTokens->At(i))->GetString().Data()));
	continue;
      } 
      else if (argument.CompareTo("-acceptTrgStr")==0) {
	fAccRejTrg=1;
	HLTInfo("Accept selected Trigger Strings only");
	continue;
      }
      else if (argument.CompareTo("-rejectTrgStr")==0) {
	fAccRejTrg=-1;
	HLTInfo("Reject all selected Trigger Strings");
	continue;
      }
      else if (argument.CompareTo("-takeAllEvents")==0) {
	fAccRejTrg=0;
	fTakeAllEvents = kTRUE;
	HLTInfo("Take all events independently of the trigger strings");
	continue;
      }
      
      else {
	HLTError("unknown argument: %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTRDCalibHistoComponent::SetParams()
{

  if(!fTrgStrings)
    fTrgStrings = new TObjArray();

  if(!AliCDBManager::Instance()->IsDefaultStorageSet()){
    HLTError("DefaultStorage is not set in CDBManager");
    return -EINVAL;
  }
  if(AliCDBManager::Instance()->GetRun()<0){
    HLTError("Run Number is not set in CDBManager");
    return -EINVAL;
  }
  HLTInfo("CDB default storage: %s; RunNo: %i", (AliCDBManager::Instance()->GetDefaultStorage()->GetBaseFolder()).Data(), AliCDBManager::Instance()->GetRun());

  if(fTrgStrings->GetEntriesFast()>0 && !fAccRejTrg){
    HLTError("Trigger string(s) given, but acceptTrgStr or rejectTrgStr not selected");
    return -EINVAL;
  }

  fTRDCalibraFillHisto = AliTRDCalibraFillHisto::Instance();
  fTRDCalibraFillHisto->SetIsHLT(kTRUE);
  fTRDCalibraFillHisto->SetHisto2d(); // choose to use histograms
  fTRDCalibraFillHisto->SetCH2dOn();  // choose to calibrate the gain
  fTRDCalibraFillHisto->SetPH2dOn();  // choose to calibrate the drift velocity
  fTRDCalibraFillHisto->SetPRF2dOn(); // choose to look at the PRF
  fTRDCalibraFillHisto->SetIsHLT(); // per detector
  //fTRDCalibraFillHisto->SetDebugLevel(1);// debug
  fTRDCalibraFillHisto->SetFillWithZero(kTRUE);
  fTRDCalibraFillHisto->SetLinearFitterOn(kTRUE);
  fTRDCalibraFillHisto->SetNumberBinCharge(100);
  
  if(!fTracksArray) fTracksArray = new TClonesArray("AliTRDtrackV1");
  if(!fOutArray)fOutArray = new TObjArray(4);

  HLTDebug("run SetupCTPData");
  SetupCTPData();

  return 0;
}

int AliHLTTRDCalibHistoComponent::DoDeinit()
{
  
  // Deinitialization of the component
  
  HLTDebug("DeinitCalibration");
  delete fTracksArray; fTracksArray=0;
  fTRDCalibraFillHisto->DestroyDebugStreamer();
  //fTRDCalibraFillHisto->Destroy();
  //fOutArray->Delete();
  delete fOutArray; fOutArray=0;
  fTrgStrings->Delete();
  delete fTrgStrings; fTrgStrings=0;
  return 0;
}

Int_t AliHLTTRDCalibHistoComponent::DoEvent(const AliHLTComponent_EventData& /*evtData*/,
					    const AliHLTComponent_BlockData* /*blocks*/,
					    AliHLTComponent_TriggerData& /*trigData*/,
					    AliHLTUInt8_t* /*outputPtr*/,
					    AliHLTUInt32_t& /*size*/,
					    vector<AliHLTComponent_BlockData>& /*outputBlocks*/)
{
  // Process an event
 
  TClonesArray* TCAarray[18] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Int_t usedEntries = 0;
  Int_t blockOrObject = 0;
  Int_t nTimeBins = -1;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(AliHLTTRDDefinitions::fgkTracksV1DataType); pBlock; pBlock=GetNextInputBlock()) 
    {
      TCAarray[0] = fTracksArray;
      AliHLTTRDUtils::ReadTracks(TCAarray[0], pBlock->fPtr, pBlock->fSize, &nTimeBins);
      fSpec |= pBlock->fSpecification;
      usedEntries = 1;
      blockOrObject = -1;
    }  

  for(const TObject *iter = GetFirstInputObject(AliHLTTRDDefinitions::fgkHiLvlTracksDataType); iter; iter = GetNextInputObject()) 
    {
      if(blockOrObject<0){
	HLTError("You may not mix high level and low level!");
	return -1;
      }

      TCAarray[usedEntries] = dynamic_cast<TClonesArray*>(const_cast<TObject*>(iter));
      if(TCAarray[usedEntries])continue;
      TObjString* strg = dynamic_cast<TObjString*>(const_cast<TObject*>(GetNextInputObject()));
      if(!strg)continue;

      nTimeBins = strg->String().Atoi();
      fSpec |= GetSpecification(iter);
      usedEntries++;
      blockOrObject = 1;
    }

  if(!blockOrObject)
    return 0;

  if(!fSavedTimeBins){
    if(nTimeBins<0){
      HLTFatal("Number of timebins is negative!");
      return -1;
    }
    HLTDebug("Saving number of time bins which was read from input block. Value is: %d", nTimeBins);
    fTRDCalibraFillHisto->Init2Dhistos(nTimeBins); // initialise the histos
    fTRDCalibraFillHisto->SetNumberClusters(fMinClusters); // At least fMinClusters clusters
    fTRDCalibraFillHisto->SetNumberClustersf(nTimeBins); // Not more than %d  clusters
    fSavedTimeBins=kTRUE;
  }

  Bool_t bTriggerPassed = fTakeAllEvents;

  if(fAccRejTrg){
    if(fAccRejTrg>0){
      bTriggerPassed=kFALSE;
      for(int i = 0; i < fTrgStrings->GetEntriesFast(); i++){
	const TObjString *const obString=(TObjString*)fTrgStrings->At(i);
	const TString tString=obString->GetString();
	if(CheckCTPTrigger(tString.Data())>0){bTriggerPassed=kTRUE; break;}
      }
    }
    else{
      bTriggerPassed=kTRUE;
      for(int i = 0; i < fTrgStrings->GetEntriesFast(); i++){
	const TObjString *const obString=(TObjString*)fTrgStrings->At(i);
	const TString tString=obString->GetString();
	if(CheckCTPTrigger(tString.Data())>0){bTriggerPassed=kFALSE; break;}
      }
    }
  }
  
  fTRDCalibraFillHisto->SetCH2dOn(bTriggerPassed);
  fTRDCalibraFillHisto->SetPH2dOn(bTriggerPassed);
  for(int i=0; i<usedEntries; i++){
    const TClonesArray *const inArr = TCAarray[i];
    Int_t nbEntries = inArr->GetEntries();
    HLTDebug(" %i TRDtracks in tracksArray", nbEntries);
    AliTRDtrackV1* trdTrack = 0x0;
    for (Int_t ii = 0; ii < nbEntries; ii++){
      HLTDebug("%i/%i: ", ii+1, nbEntries);
      trdTrack = (AliTRDtrackV1*)inArr->At(ii);
      if(trdTrack->GetNumberOfTracklets()<fMinTracklets)continue;
      fTRDCalibraFillHisto->UpdateHistogramsV1(trdTrack);
      // for(int i3=0; i3<7; i3++)
      //   if(trdTrack->GetTracklet(i3))trdTrack->GetTracklet(i3)->Bootstrap(fReconstructor);
    }
  }

  if(!fOutArray->At(0))FormOutput();
  PushBack(fOutArray, AliHLTTRDDefinitions::fgkCalibrationDataType, fSpec);

  if(blockOrObject<0){
    TCAarray[0]->Delete();
  }

  return 0;
}

/**
 * Form output array of histrograms
 */
//============================================================================
void AliHLTTRDCalibHistoComponent::FormOutput()
{
  // gain histo
  TH2I *hCH2d = fTRDCalibraFillHisto->GetCH2d();
  fOutArray->Add(hCH2d);
  
  // drift velocity histo
  TProfile2D *hPH2d = fTRDCalibraFillHisto->GetPH2d();
  fOutArray->Add(hPH2d);
  
  // PRF histo
  TProfile2D *hPRF2d = fTRDCalibraFillHisto->GetPRF2d();
  fOutArray->Add(hPRF2d);
  
  // Vdrift Linear Fit
  AliTRDCalibraVdriftLinearFit *hVdriftLinearFitOne=(AliTRDCalibraVdriftLinearFit *)fTRDCalibraFillHisto->GetVdriftLinearFit();
  fOutArray->Add(hVdriftLinearFitOne);
  
  HLTDebug("GetCH2d = 0x%x; NEntries = %i; size = %i", hCH2d, hCH2d->GetEntries(), sizeof(*hCH2d));
  hCH2d->Print();
  HLTDebug("GetPH2d = 0x%x; NEntries = %i; size = %i", hPH2d, hPH2d->GetEntries(), sizeof(*hPH2d));
  hPH2d->Print();
  HLTDebug("GetPRF2d = 0x%x; NEntries = %i; size = %i", hPRF2d, hPRF2d->GetEntries(), sizeof(*hPRF2d));
  hPRF2d->Print();
  HLTDebug("GetVdriftLinearFit = 0x%x; size = %i", hVdriftLinearFitOne, sizeof(hVdriftLinearFitOne)); 
  
  HLTDebug("output Array: pointer = 0x%x; NEntries = %i; size = %i", fOutArray, fOutArray->GetEntries(), sizeof(fOutArray));
   
}

int AliHLTTRDCalibHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{
  // see header file for class documentation

  int iResult=0;
  const char* path="HLT/ConfigTRD/CalibHistoComponent";
  const char* defaultNotify="";
  if (cdbEntry) {
    path=cdbEntry;
    defaultNotify=" (default)";
  }
  if (path) {
    HLTInfo("reconfigure from entry %s%s, chain id %s", path, defaultNotify,(chainId!=NULL && chainId[0]!=0)?chainId:"<none>");
    AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
    if (pEntry) {
      TObjString* pString=dynamic_cast<TObjString*>(pEntry->GetObject());
      if (pString) {
  	HLTInfo("received configuration object string: \'%s\'", pString->GetString().Data());
  	iResult=Configure(pString->GetString().Data());
      } else {
  	HLTError("configuration object \"%s\" has wrong type, required TObjString", path);
      }
    } else {
      HLTError("cannot fetch object \"%s\" from CDB", path);
    }
  }

  return iResult;
}
