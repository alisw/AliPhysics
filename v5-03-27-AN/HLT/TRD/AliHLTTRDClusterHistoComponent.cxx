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

/** @file   AliHLTTRDClusterHistoComponent.cxx
    @author Sylwester Radomski
    @brief  Component for ploting charge in clusters
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TClonesArray.h"
#include "TH1F.h"

#include "AliHLTTRDClusterHistoComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliTRDcluster.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTTRDUtils.h"

//#include "AliHLTTRD.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDClusterHistoComponent)

AliHLTTRDClusterHistoComponent::AliHLTTRDClusterHistoComponent()
: AliHLTProcessor(),
  fOutputSize(100000),
  fSpec(0),
  fClusterArray(NULL),
  fNClsDet(NULL),
  fClsAmp(NULL),
  fClsAmpDrift(NULL),
  fClsTB(NULL),
  fClsAmpDriftDet(),
  fClsAmpDist(NULL),
  fSClsDist(NULL),
  fNScls(NULL),
  fEvSize(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  memset(fClsAmpDriftDet, 0, sizeof(fClsAmpDriftDet));
}

AliHLTTRDClusterHistoComponent::~AliHLTTRDClusterHistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTRDClusterHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "TRDClusterHisto";
}

void AliHLTTRDClusterHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTRDDefinitions::fgkClusterDataType );
}

AliHLTComponentDataType AliHLTTRDClusterHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram  | kAliHLTDataOriginTRD;

}

void AliHLTTRDClusterHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDClusterHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTRDClusterHistoComponent;
}

int AliHLTTRDClusterHistoComponent::DoInit(int argc, const char** argv)
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

  fClusterArray = new TClonesArray("AliTRDcluster");

  fNClsDet = new TH1F("trdClsDet", ";detector", 540, -0.5, 539.5);
  fClsAmp  = new TH1F("trdClsAmp", ";amplitude", 200, -0.5, 1999.5);
  fClsAmpDrift = new TH1F("trdClsAmpDrift", ";amplitude", 200, -0.5, 199.5) ;
  fClsTB = new TH1F("trdClsTB", ";time bin", 35, -0.5, 34.5);
  fClsAmpDist = new TH1F("trdClsAmpDist", "mean amplitude", 200, 0, 1000);
  fSClsDist = new TH1F("sclsdist", "Super cluster spectrum", 200, 0, 8000);
  fNScls = new TH1F("nscls", "No. of Kr clusters per event", 540, 0, 540);
  fEvSize = new TH1F("TrdClEvSize", "Clusters size per event per ddl in kbyte", 512, 0, 512);

  for(int i=0; i<540; i++) {
    fClsAmpDriftDet[i] = new TH1F(Form("trdClsDriftDet_%d",i), "", 200, -0.5, 199.5);
  }
  
  return 0;
}
  
int AliHLTTRDClusterHistoComponent::DoDeinit()
{
  // see header file for class documentation

  fClusterArray->Delete();
  delete fClusterArray;

  // delete histograms
  if (fNClsDet) delete fNClsDet;
  if (fClsAmp) delete fClsAmp;
  if (fClsAmpDrift) delete fClsAmpDrift;
  if (fClsTB) delete fClsTB;
  if (fClsAmpDist) delete fClsAmpDist;
  if (fSClsDist) delete fSClsDist;
  if (fNScls) delete fNScls;
  if (fEvSize) delete fEvSize;

  for(int i=0; i<540; i++){
    if (fClsAmpDriftDet[i]) delete fClsAmpDriftDet[i];
  }

  return 0;
}

int AliHLTTRDClusterHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					    AliHLTComponentTriggerData& /*trigData*/)
{

  // if (GetFirstInputBlock(kAliHLTDataTypeSOR)) return 0;
  // else if (GetFirstInputBlock(kAliHLTDataTypeEOR))
  //   {
  //     TString fileName="/tmp/ClusterHistoDump_run";
  //     fileName+=AliCDBManager::Instance()->GetRun();
  //     fileName+=".root";
  //     HLTInfo("Dumping Histogram file to %s",fileName.Data());
  //     TFile* file = TFile::Open(fileName, "RECREATE");
  //     fNClsDet->Write();
  //     fClsAmp->Write();
  //     fClsAmpDrift->Write();
  //     fClsTB->Write();
  //     fClsAmpDist->Write(); 
  //     fSClsDist->Write();
  //     fNScls->Write();
  //     file->Close();
  //     HLTInfo("Histogram file dumped");
  //     return 0;
  //   }

  if(!IsDataEvent())return 0;

  const AliHLTComponentBlockData* iter = NULL;
  Bool_t gotData = kFALSE;

  for ( iter = GetFirstInputBlock(AliHLTTRDDefinitions::fgkClusterDataType); 
	iter != NULL; iter = GetNextInputBlock() ) {

    fEvSize->Fill((iter->fSize+0.5f)/1024);
    AliHLTTRDUtils::ReadClusters(fClusterArray, iter->fPtr, iter->fSize);
    HLTDebug("TClonesArray of clusters: nbEntries = %i", fClusterArray->GetEntriesFast());
    gotData = kTRUE;
    fSpec |= iter->fSpecification;
  }

  if(!gotData) return 0;

  Float_t sClusterCharge[540] = { 0 };
    AliTRDcluster *cls;

    // loop over clusters 
    for(int i=0;i<fClusterArray->GetEntriesFast();i++) {

      cls=(AliTRDcluster*)fClusterArray->At(i);
      
      fNClsDet->Fill(cls->GetDetector());
      fClsAmp->Fill(cls->GetQ());
      
      int tb = cls->GetPadTime();
      fClsTB->Fill(tb);
    if (tb > 5 && tb <25){
	fClsAmpDrift->Fill(cls->GetQ()); 
    }
      
    //fClsAmpDriftDet[cls->GetDetector()]->Fill(cls->GetQ());

      Int_t det = cls->GetDetector();
      sClusterCharge[det] += cls->GetQ();

    }
    
    fClusterArray->Delete();
    
  //fClsAmpDist->Reset();
  //Int_t nSClusters = 0;
  for(int det=0; det<540; det++) {
    // if (fClsAmpDriftDet[det]->GetSum() > 0) 
    //   fClsAmpDist->Fill(fClsAmpDriftDet[det]->GetMean());
    if(sClusterCharge[det])
    fSClsDist->Fill(sClusterCharge[det]);
  }

  //fNScls->Fill(nSClusters);

  PushBack((TObject*)fNClsDet, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fClsAmp, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fClsAmpDrift, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fClsTB, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  //PushBack((TObject*)fClsAmpDist, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  //PushBack((TObject*)fNScls, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fSClsDist, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);
  PushBack((TObject*)fEvSize, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, fSpec);

  return 0;
}

int AliHLTTRDClusterHistoComponent::Configure(const char* arguments){
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
