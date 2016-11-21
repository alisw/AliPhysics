//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sylwester Radomski radomski@physi.uni-heidelberg.de    *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTRDTrackHistoComponent.cxx
    @author Raphaelle and Theodor
    @brief  Component for ploting charge in clusters
*/

#include <time.h>

#include "AliHLTTRDTrackHistoComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TClonesArray.h"
#include "TTimeStamp.h"
#include "AliHLTTRDUtils.h"
#include "TH1F.h"
#include "AliTRDcluster.h"
#include "AliTRDtrackV1.h"
#include "AliTRDseedV1.h"

//#include "AliHLTTRD.h"
//#include <stdlib.h>
//#include <cerrno>

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDTrackHistoComponent)

AliHLTTRDTrackHistoComponent::AliHLTTRDTrackHistoComponent()
: AliHLTProcessor(),
  fOutputSize(100000),
  fSpec(0),
  fTracksArray(NULL),
  fClPerTrkl(NULL),
  fTrklPerTrk(NULL),
  fEvSize(NULL),
  fEtaDistrib(NULL),
  fPhiDistrib(NULL),
  fPtDistrib(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTTRDTrackHistoComponent::~AliHLTTRDTrackHistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTRDTrackHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "TRDTrackHisto";
}

void AliHLTTRDTrackHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTRDDefinitions::fgkTracksV1DataType );
}

AliHLTComponentDataType AliHLTTRDTrackHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram  | kAliHLTDataOriginTRD;

}

void AliHLTTRDTrackHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDTrackHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTRDTrackHistoComponent;
}

int AliHLTTRDTrackHistoComponent::DoInit(int argc, const char** argv)
{
  // Initialize histograms
  int iResult=0;
  
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }

  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } 

  fTracksArray = new TClonesArray("AliTRDtrackV1");

  fClPerTrkl = new TH1F("TrdClPerTrkl","Clusters per Tracklet", AliTRDseedV1::kNtb, -0.5, AliTRDseedV1::kNtb - 0.5);
  fTrklPerTrk = new TH1F("TrdTrklPerTrk","Tracklets per Track", 7, -0.5, 6.5);
  fEvSize = new TH1F("TrdTrEvSize", "Tracks size per event per ddl in kbyte", 512, 0, 512);
  fEtaDistrib = new TH1F("TrdTrEtaDistrib", "Eta distribution of tracks", 51, -1, 1);
  fPhiDistrib = new TH1F("TrdTrPhiDistrib", "Phi distribution of tracks", 63, 0, 6.3);
  fPtDistrib = new TH1F("TrdTrPtDistrib", "Pt distribution of tracks", 101, 0, 10);
  return 0;
}
  
int AliHLTTRDTrackHistoComponent::DoDeinit()
{
  // see header file for class documentation

  fTracksArray->Delete();
  delete fTracksArray;

  // delete histograms
  if (fClPerTrkl) delete fClPerTrkl;
  if (fTrklPerTrk) delete fTrklPerTrk;
  
  return 0;
}

int AliHLTTRDTrackHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					    AliHLTComponentTriggerData& /*trigData*/)
{

  // if (GetFirstInputBlock(kAliHLTDataTypeSOR)) return 0;
  // else if (GetFirstInputBlock(kAliHLTDataTypeEOR))
  //   {
  //     TString fileName="/tmp/TracksHistoDump_run";
  //     fileName+=AliCDBManager::Instance()->GetRun();
  //     fileName+=".root";
  //     HLTInfo("Dumping Histogram file to %s",fileName.Data());
  //     TFile* file = TFile::Open(fileName, "RECREATE");
  //     fClPerTrkl->Write();
  //     fTrklPerTrk->Write();
  //     file->Close();
  //     HLTInfo("Histogram file dumped");
  //     return 0;
  //   }

  if(!IsDataEvent())return 0;

  const AliHLTComponentBlockData* iter = NULL;
  Bool_t gotData = kFALSE;
  
  for(iter = GetFirstInputBlock(AliHLTTRDDefinitions::fgkTracksV1DataType); 
	iter != NULL; iter = GetNextInputBlock() ) {
    
    fEvSize->Fill((iter->fSize+0.5f)/1024);
    AliHLTTRDUtils::ReadTracks(fTracksArray, iter->fPtr, iter->fSize);
    HLTDebug("TClonesArray of tracks: nbEntries = %i", fTracksArray->GetEntriesFast());
    gotData=kTRUE;
    fSpec |= iter->fSpecification;
  }
  
  if(!gotData) return 0;
  
  AliTRDtrackV1 *trk;
  
  // loop over tracks
  for(int i=0;i<fTracksArray->GetEntriesFast();i++) {
    trk=(AliTRDtrackV1*)fTracksArray->At(i);
    fEtaDistrib->Fill(trk->Eta());
    fPhiDistrib->Fill(trk->Phi());
    fPtDistrib->Fill(trk->Pt());
    Int_t nrOfTrkls=0;
    for(int seedNr=0; seedNr<6; seedNr++){
      AliTRDseedV1* seed = trk->GetTracklet(seedNr);
      if(!seed)continue;
      nrOfTrkls++;
      Int_t nrOfCls=0;
      for(int clsNr=0; clsNr<AliTRDseedV1::kNtb; clsNr++)
	if(seed->GetClusters(clsNr))nrOfCls++;
      fClPerTrkl->Fill(nrOfCls);
    }
    fTrklPerTrk->Fill(nrOfTrkls);
  }
  
  fTracksArray->Delete();
  
  PushBack((TObject*)fClPerTrkl, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);   
  PushBack((TObject*)fTrklPerTrk, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);  
  PushBack((TObject*)fEvSize, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fEtaDistrib, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);   
  PushBack((TObject*)fPhiDistrib, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);  
  PushBack((TObject*)fPtDistrib, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  
  return 0;
}

int AliHLTTRDTrackHistoComponent::Configure(const char* arguments){
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
      if (argument.CompareTo("-everyNevent")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Option -everyNevent depreceated");
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
