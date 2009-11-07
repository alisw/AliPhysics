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

#include <time.h>

#include "AliHLTTRDClusterHistoComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDCluster.h"
#include "AliTRDcluster.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TClonesArray.h"
#include "TTimeStamp.h"
#include "AliHLTTRDUtils.h"
#include "TH2F.h"

//#include "AliHLTTRD.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDClusterHistoComponent)

AliHLTTRDClusterHistoComponent::AliHLTTRDClusterHistoComponent()
: AliHLTProcessor(),
  fOutputSize(100000),
  fClusterArray(NULL),
  fNClsDet(NULL),
  fClsAmp(NULL),
  fClsAmpDrift(NULL),
  fClsTB(NULL),
  fClsAmpDist(NULL),
  fSClsDist(NULL),
  fNScls(NULL),
  fClusterDist(NULL),
  fClusterCandCharge(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

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

  fNClsDet = new TH1D("trdClsDet", ";detector", 540, -0.5, 539.5);
  fClsAmp  = new TH1D("trdClsAmp", ";amplitude", 200, -0.5, 1999.5);
  fClsAmpDrift = new TH1D("trdClsAmpDrift", ";amplitude", 200, -0.5, 199.5) ;
  fClsTB = new TH1D("trdClsTB", ";time bin", 35, -0.5, 34.5);
  fClsAmpDist = new TH1D("trdClsAmpDist", "mean amplitude", 200, 0, 1000);
  fSClsDist = new TH1D("sclsdist", "Super cluster spectrum", 200, 0, 8000);
  fNScls = new TH1D("nscls", "No. of Kr clusters per event", 540, 0, 540);
  fClusterDist = new TH2F("cldist", "Cluster distribution;padrow;padcol", 16*5, -0.5, 79.5, 8*6*18, -0.5, 863.5);
  fClusterCandCharge = new TH1D("qClsCand", "Cluster candidate charge;charge (ADC counts);counts", 200, 0, 8000);

  for(int i=0; i<540; i++) {
    fClsAmpDriftDet[i] = new TH1D(Form("trdClsDriftDet_%d",i), "", 200, -0.5, 199.5);
    fSlidingWindow[i].SetBins(9, -0.5, 8.5, 17, -0.5, 16.5);
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
  if (fClusterDist) delete fClusterDist;
  if (fClusterCandCharge) delete fClusterCandCharge;
  if (fNScls) delete fNScls;

  for(int i=0; i<540; i++)
    if (fClsAmpDriftDet[i]) delete fClsAmpDriftDet[i];

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
  //     fClusterDist->Write();
  //     fClusterCandCharge->Write();
  //     file->Close();
  //     HLTInfo("Histogram file dumped");
  //     return 0;
  //   }

  if (GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR)) return 0;

  const AliHLTComponentBlockData* iter = NULL;
  
  Float_t sClusterCharge[540] = { 0 };
  for (Int_t iDet = 0; iDet < 540; iDet++) {
    fSlidingWindow[iDet].Reset();
  }

  for ( iter = GetFirstInputBlock(AliHLTTRDDefinitions::fgkClusterDataType); 
	iter != NULL; iter = GetNextInputBlock() ) {

    HLTDebug("We get the right data type: Block Ptr: 0x%x; Block Size: %i",
	     iter->fPtr, iter->fSize);

    AliHLTTRDUtils::ReadClusters(fClusterArray, iter->fPtr, iter->fSize);
    HLTDebug("TClonesArray of clusters: nbEntries = %i", fClusterArray->GetEntriesFast());

    AliTRDcluster *cls;

    // loop over clusters 
    for(int i=0;i<fClusterArray->GetEntriesFast();i++) {

      cls=(AliTRDcluster*)fClusterArray->At(i);
      
      fNClsDet->Fill(cls->GetDetector());
      fClsAmp->Fill(cls->GetQ());
      
      int tb = cls->GetPadTime();
      fClsTB->Fill(tb);
      if (tb > 5 && tb <25)
	fClsAmpDrift->Fill(cls->GetQ()); 
      
      fClsAmpDriftDet[cls->GetDetector()]->Fill(cls->GetQ());

      Int_t det = cls->GetDetector();
      sClusterCharge[det] += cls->GetQ();

      fSlidingWindow[det].Fill((cls->GetPadCol() / 18), cls->GetPadRow(), cls->GetQ());
      fSlidingWindow[det].Fill((cls->GetPadCol() / 18), cls->GetPadRow() + 1, cls->GetQ());
      fSlidingWindow[det].Fill((cls->GetPadCol() / 18) + 1, cls->GetPadRow(), cls->GetQ());
      fSlidingWindow[det].Fill((cls->GetPadCol() / 18) + 1, cls->GetPadRow() + 1, cls->GetQ());
    }
    
    fClusterArray->Delete();
    
  }
   
  fClsAmpDist->Reset();
  Int_t nSClusters = 0;
  for(int det=0; det<540; det++) {
    if (fClsAmpDriftDet[det]->GetSum() > 0) 
      fClsAmpDist->Fill(fClsAmpDriftDet[det]->GetMean());
    fSClsDist->Fill(sClusterCharge[det]);

    Int_t xmax;
    Int_t ymax;
    Int_t zmax;
    fSlidingWindow[det].GetMaximumBin(xmax, ymax, zmax);
    Float_t charge = fSlidingWindow[det].GetBinContent(xmax, ymax);
    fClusterCandCharge->Fill(charge);
    if (charge > 2000. && charge < 6500.) {
      nSClusters++;
      fClusterDist->Fill((ymax-1) + 16 * ((det % 30) / 6), (det / 30) * 48 + 8 * (det % 6) + (xmax-1));
    }
  }

  fNScls->Fill(nSClusters);

  PushBack((TObject*)fNClsDet, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);   
  PushBack((TObject*)fClsAmp, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);  
  //  PushBack((TObject*)fClsAmpDrift, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);   
  //  PushBack((TObject*)fClsTB, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);  
  //  PushBack((TObject*)fClsAmpDist, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);  
  PushBack((TObject*)fNScls, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);  
  PushBack((TObject*)fSClsDist, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);  
  PushBack((TObject*)fClusterDist, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);
  PushBack((TObject*)fClusterCandCharge, kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);

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
