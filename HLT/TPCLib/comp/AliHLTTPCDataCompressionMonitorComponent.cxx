// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/// @file   AliHLTTPCDataCompressionMonitorComponent.cxx
/// @author Matthias Richter
/// @date   2011-09-12
/// @brief  TPC component for monitoring of data compression
///

#include "AliHLTTPCDataCompressionMonitorComponent.h"
#include "AliHLTTPCDataCompressionComponent.h"
#include "AliHLTTPCDataCompressionDecoder.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCHWCFData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrackGeometry.h"
#include "AliHLTTPCHWCFSpacePointContainer.h"
#include "AliHLTErrorGuard.h"
#include "AliRawDataHeader.h"
#include "AliTPCclusterMI.h"
#include "AliTPCROC.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TH3I.h"
#include "TH3F.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TList.h"
#include <memory>

ClassImp(AliHLTTPCDataCompressionMonitorComponent)

AliHLTTPCDataCompressionMonitorComponent::AliHLTTPCDataCompressionMonitorComponent()
  : AliHLTProcessor()
  , fpHWClusterDecoder(NULL)
  , fHistoHWCFDataSize(NULL)
  , fHistoHWCFReductionFactor(NULL)
  , fHistoNofClusters(NULL)
  , fHistoNofClustersReductionFactor(NULL)
  , fHistogramFile()
  , fMonitoringContainer(NULL)
  , fVerbosity(0)
  , fFlags(0)
  , fPublishingMode(kPublishSeparate)
{
}

AliHLTTPCDataCompressionMonitorComponent::~AliHLTTPCDataCompressionMonitorComponent()
{
  /// destructor
}


const char* AliHLTTPCDataCompressionMonitorComponent::GetComponentID()
{
  /// inherited from AliHLTComponent: id of the component
  return "TPCDataCompressorMonitor";
}


void AliHLTTPCDataCompressionMonitorComponent::GetInputDataTypes( AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: list of data types in the vector reference
  tgtList.clear();
  tgtList.push_back(AliHLTTPCDefinitions::fgkHWClustersDataType);
  tgtList.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);
  tgtList.push_back(AliHLTTPCDefinitions::fgkRawClustersDataType);
  tgtList.push_back(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
  tgtList.push_back(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());  
}

AliHLTComponentDataType AliHLTTPCDataCompressionMonitorComponent::GetOutputDataType()
{
  /// inherited from AliHLTComponent: output data type of the component.
  return kAliHLTMultipleDataType;
}

int AliHLTTPCDataCompressionMonitorComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  /// inherited from AliHLTComponent: multiple output data types of the component.
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
  return tgtList.size();
}

void AliHLTTPCDataCompressionMonitorComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  /// inherited from AliHLTComponent: output data size estimator
  constBase=100000;
  inputMultiplier=1.0;
}

AliHLTComponent* AliHLTTPCDataCompressionMonitorComponent::Spawn()
{
  /// inherited from AliHLTComponent: spawn function.
  return new AliHLTTPCDataCompressionMonitorComponent;
}

int AliHLTTPCDataCompressionMonitorComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, 
						       const AliHLTComponentBlockData* /*inputBlocks*/, 
						       AliHLTComponentTriggerData& /*trigData*/,
						       AliHLTUInt8_t* /*outputPtr*/,
						       AliHLTUInt32_t& /*size*/,
						       AliHLTComponentBlockDataList& /*outputBlocks*/ )
{
  /// inherited from AliHLTProcessor: data processing
  int iResult=0;

  if (!IsDataEvent()) return 0;

  const AliHLTComponentBlockData* pDesc=NULL;
  unsigned rawDataSize=0;
  unsigned hwclustersDataSize=0;
  unsigned nofClusters=0;
  unsigned compDataSize=0; 
  
  // check size of TPC raw data
  for (pDesc=GetFirstInputBlock(kAliHLTDataTypeDDLRaw | kAliHLTDataOriginTPC);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    fFlags|=kHaveRawData;
    rawDataSize+=pDesc->fSize;
  }

  // check size of HWCF data and add to the MonitoringContainer
  for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::fgkHWClustersDataType);
       pDesc!=NULL; pDesc=GetNextInputBlock()) {
    fFlags|=kHaveHWClusters;
    // FIXME: the decoding can now be handled via the data container
    if (pDesc->fSize<=sizeof(AliRawDataHeader)) continue;
    if (fpHWClusterDecoder) {
      hwclustersDataSize+=pDesc->fSize;
      AliHLTUInt8_t* pData=reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr);
      pData+=sizeof(AliRawDataHeader);
      if (fpHWClusterDecoder->Init(pData, pDesc->fSize-sizeof(AliRawDataHeader))<0 ||
	  (fpHWClusterDecoder->CheckVersion()<0 && (int)pDesc->fSize>fpHWClusterDecoder->GetRCUTrailerSize())) {
	HLTError("data block of type %s corrupted: can not decode format",
		 AliHLTComponent::DataType2Text(pDesc->fDataType).c_str());
      } else {
	nofClusters+=fpHWClusterDecoder->GetNumberOfClusters();
      }
    }
    if (fMonitoringContainer) {
      fMonitoringContainer->AddRawData(pDesc);
    }
  }

  if (fMonitoringContainer) {
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClusterIdsDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=fMonitoringContainer->AddClusterIds(pDesc);
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterIdTracksDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=fMonitoringContainer->AddClusterIds(pDesc);
    }

    // read data
    AliHLTTPCDataCompressionDecoder decoder;
    bool bHaveRawClusters=false;
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RawClustersDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      // Note: until r51411 and v5-01-Rev-03 the compressed cluster format was sent with data
      // type {CLUSTRAW,TPC }, the version member indicated the actual type of data
      // These data do not include the 0.5 shift in pad position, that's wht it has
      // to be added in the unpacking. This is a very special case, this data type and
      // data version==1 only occured in the early TPC data compression test runs with
      // v5-01-Rev-01
      if (pDesc->fSize<sizeof(AliHLTTPCRawClusterData)) continue;
      const AliHLTTPCRawClusterData* clusterData = reinterpret_cast<const AliHLTTPCRawClusterData*>(pDesc->fPtr);
      if (!clusterData) continue;
      if (clusterData->fVersion==1) {
	// compressed clusters without the pad shift
	// data type {CLUSTRAW,TPC } with version==1
	decoder.SetPadShift(0.5);
      } else {
	decoder.SetPadShift(0.0);
      }
      bHaveRawClusters=true;
      iResult=decoder.ReadClustersPartition(fMonitoringContainer->BeginRemainingClusterBlock(0, pDesc->fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					    pDesc->fSize,
					    pDesc->fSpecification);
      if (iResult<0) {
	HLTError("reading of partition clusters failed with error %d", iResult);
      }
    }

    decoder.SetPadShift(0.0);

    if (!bHaveRawClusters) {
    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=decoder.ReadClustersPartition(fMonitoringContainer->BeginRemainingClusterBlock(0, pDesc->fSpecification),
					    reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					    pDesc->fSize,
					    pDesc->fSpecification);
      compDataSize+=pDesc->fSize;
    }

    for (pDesc=GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType());
	 pDesc!=NULL; pDesc=GetNextInputBlock()) {
      iResult=decoder.ReadTrackModelClustersCompressed(fMonitoringContainer->BeginTrackModelClusterBlock(0),
					       reinterpret_cast<AliHLTUInt8_t*>(pDesc->fPtr),
					       pDesc->fSize,
					       pDesc->fSpecification);
      compDataSize+=pDesc->fSize;
    }
    } else {
      if (GetFirstInputBlock(AliHLTTPCDefinitions::RemainingClustersCompressedDataType()) ||
	  GetFirstInputBlock(AliHLTTPCDefinitions::ClusterTracksCompressedDataType())) {
	ALIHLTERRORGUARD(5, "conflicting data blocks, monitoring histograms already filled from raw cluster data, ignoring blocks of compressed partition and track clusters");
      }		     
    }

    fMonitoringContainer->Clear();
  }

  float ratio=0;
  if (compDataSize) {ratio=(float)hwclustersDataSize; ratio/=compDataSize;}
  if (fHistoHWCFDataSize)        fHistoHWCFDataSize       ->Fill(rawDataSize/1024, compDataSize/1024);
  if (fHistoHWCFReductionFactor) fHistoHWCFReductionFactor->Fill(rawDataSize/1024, ratio);
  if (fHistoNofClusters)         fHistoNofClusters        ->Fill(rawDataSize/1024, nofClusters);
  if (fHistoNofClustersReductionFactor) fHistoNofClustersReductionFactor ->Fill(nofClusters, ratio);
  HLTInfo("raw data %d, hwcf data %d, comp data %d, ratio %f, %d clusters\n", rawDataSize, hwclustersDataSize, compDataSize, ratio, nofClusters);

  if (iResult>=0 && fPublishingMode!=kPublishOff) {
    iResult=Publish(fPublishingMode);
  }

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::Publish(int mode)
{
  /// publish to output
  // additional histograms derived from the main ones to publish
  TObjArray *derivedHistos = new TObjArray();
  derivedHistos->SetOwner(kTRUE);

  // FIXME: code needs to be optimized, maybe a bit to much new and delete for the
  // moment, the data type might need adjustment
  int iResult=0;
  TObjArray* pArray=mode==kPublishArray?(new TObjArray):NULL;
  TList* pList=mode==kPublishList?(new TList):NULL;
  if (mode==kPublishSeparate) {
    if (fHistoHWCFDataSize)        PushBack(fHistoHWCFDataSize       , kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
    if (fHistoHWCFReductionFactor) PushBack(fHistoHWCFReductionFactor, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
    if (fHistoNofClusters)         PushBack(fHistoNofClusters        , kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
    if (fHistoNofClustersReductionFactor) PushBack(fHistoNofClustersReductionFactor, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
  } else if (pList) {
    if (fHistoHWCFDataSize)        pList->Add(fHistoHWCFDataSize->Clone());
    if (fHistoHWCFReductionFactor) pList->Add(fHistoHWCFReductionFactor->Clone());
    if (fHistoNofClusters)         pList->Add(fHistoNofClusters->Clone());
    if (fHistoNofClustersReductionFactor) pList->Add(fHistoNofClustersReductionFactor->Clone());
  } else if (pArray) {
    if (fHistoHWCFDataSize)        pArray->Add(fHistoHWCFDataSize->Clone());
    if (fHistoHWCFReductionFactor) pArray->Add(fHistoHWCFReductionFactor->Clone());
    if (fHistoNofClusters)         pArray->Add(fHistoNofClusters->Clone());
    if (fHistoNofClustersReductionFactor) pArray->Add(fHistoNofClustersReductionFactor->Clone());
  }


  if (fMonitoringContainer) {
    static const char* searchIds[] = {"fHistograms", "fHistograms2D", "fHistograms3D", NULL};
    const char** searchId=searchIds;
    while (*searchId && iResult>=0) {
      const TObject* o=fMonitoringContainer->FindObject(*searchId);
      if (o) {
	const TObjArray* histograms=dynamic_cast<const TObjArray*>(o);
	if (histograms) {
	  for (int i=0; i<histograms->GetEntriesFast() && iResult>=0; i++) {
	    if (!histograms->At(i)) continue;
	    ///
	    TString name=histograms->At(i)->GetName();
	    if( (name.CompareTo(fgkHistogramDefinitions2D[kHistogramQMaxSector].fName)==0) ||
		(name.CompareTo(fgkHistogramDefinitions2D[kHistogramSigmaY2Sector].fName)==0) ||
		(name.CompareTo(fgkHistogramDefinitions2D[kHistogramSigmaZ2Sector].fName)==0) ){
	      TH2F *h1=(TH2F*)histograms->At(i);
	      TProfile *h2 = (TProfile*)(h1->ProfileX());
	      derivedHistos->Add(h2);
	    }
	    if( name.CompareTo(fgkHistogramDefinitions3D[kHistogramPadrowPadSector].fName)==0) {
	      TH3F *h1=(TH3F*)histograms->At(i);
	      for (int j=1; j<=72; j++) {
	      h1->GetXaxis()->SetRange(j,j);
	      TString histoname = Form("zy_%d",j);
	      TH2F *h2 = (TH2F*)h1->Project3D(histoname.Data());
	      derivedHistos->Add(h2);
	      }
	    }
	    ///
	    if (mode==kPublishSeparate) {
	      iResult=PushBack(histograms->At(i), kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
	    } else if (pList) {
	      pList->Add(histograms->At(i)->Clone());
	    } else if (pArray) {
	      pArray->Add(histograms->At(i)->Clone());
	    }
	  }
	  for (int i=0; i<derivedHistos->GetEntriesFast() && iResult>=0; i++) {
	    if (mode==kPublishSeparate) {
	      iResult=PushBack(derivedHistos->At(i), kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
	    } else if (pList) {
	      pList->Add(derivedHistos->At(i)->Clone());
	    } else if (pArray) {
	      pArray->Add(derivedHistos->At(i)->Clone());
	    }	    
	  }
	}
      } else {
	HLTError("failed to find object \"%s\"", *searchId);
      }
      searchId++;
    }
  }

  if (pArray) {
    iResult=PushBack(pArray, kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC);
    pArray->SetOwner(kTRUE);
    delete pArray;
    pArray=NULL;
  }
  if (pList) {
    iResult=PushBack(pList, kAliHLTDataTypeTObject|kAliHLTDataOriginTPC);
    pList->SetOwner(kTRUE);
    delete pList;
    pList=NULL;
  }
  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::DoInit( int argc, const char** argv )
{
  /// inherited from AliHLTComponent: component initialisation and argument scan.
  int iResult=0;

  // component configuration
  //Stage 1: default initialization.
  //Default values.
  fFlags=0;

  //Stage 2: OCDB.
  TString cdbPath("HLT/ConfigTPC/");
  cdbPath += GetComponentID();
  //
  // iResult = ConfigureFromCDBTObjString(cdbPath);
  // if (iResult < 0) 
  //   return iResult;

  //Stage 3: command line arguments.
  if (argc && (iResult = ConfigureFromArgumentString(argc, argv)) < 0)
    return iResult;

  std::auto_ptr<AliHLTTPCHWCFData> hwClusterDecoder(new AliHLTTPCHWCFData);
  std::auto_ptr<AliDataContainer> dataContainer(new AliDataContainer);

  std::auto_ptr<TH2I> histoHWCFDataSize(new TH2I("HWCFDataSize",
						 "HW ClusterFinder Size",
						 100, 0., 80000., 100, 0., 80000.));
  if (histoHWCFDataSize.get()) {
    TAxis* xaxis=histoHWCFDataSize->GetXaxis();
    if (xaxis) xaxis->SetTitle("raw data size [kB]");
    TAxis* yaxis=histoHWCFDataSize->GetYaxis();
    if (yaxis) yaxis->SetTitle("compressed data size [kb]");
  }

  std::auto_ptr<TH2I> histoHWCFReductionFactor(new TH2I("HWCFReductionFactor",
							"Data reduction HW ClusterFinder",
							100, 0., 80000., 100, 0., 10.));
  if (histoHWCFReductionFactor.get()) {
    TAxis* xaxis=histoHWCFReductionFactor->GetXaxis();
    if (xaxis) xaxis->SetTitle("raw data size [kB]");
    TAxis* yaxis=histoHWCFReductionFactor->GetYaxis();
    if (yaxis) yaxis->SetTitle("reduction factor");
  }

  std::auto_ptr<TH2I> histoNofClusters(new TH2I("NofClusters",
					       "Number of HLT TPC clusters",
					       100, 0., 80000., 500, 0., 3000000.));
  if (histoNofClusters.get()) {
    TAxis* xaxis=histoNofClusters->GetXaxis();
    if (xaxis) xaxis->SetTitle("raw data size [kB]");
    TAxis* yaxis=histoNofClusters->GetYaxis();
    if (yaxis) yaxis->SetTitle("N. of clusters");
  }

  std::auto_ptr<TH2I> histoNofClustersReductionFactor(new TH2I("ReductionFactorVsNofClusters",
							       "Reduction Factor vs. Number of HLT TPC clusters",
							       500, 0., 3000000., 100, 0., 10.));
  if (histoNofClustersReductionFactor.get()) {
    TAxis* xaxis=histoNofClustersReductionFactor->GetXaxis();
    if (xaxis) xaxis->SetTitle("N. of clusters");
    TAxis* yaxis=histoNofClustersReductionFactor->GetYaxis();
    if (yaxis) yaxis->SetTitle("reduction factor");
  }

  fHistoHWCFDataSize=histoHWCFDataSize.release();
  fHistoHWCFReductionFactor=histoHWCFReductionFactor.release();
  fHistoNofClusters=histoNofClusters.release();
  fHistoNofClustersReductionFactor=histoNofClustersReductionFactor.release();

  fpHWClusterDecoder=hwClusterDecoder.release();
  fMonitoringContainer=dataContainer.release();

  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::DoDeinit()
{
  /// inherited from AliHLTComponent: component cleanup
  int iResult=0;

  if (fpHWClusterDecoder) delete fpHWClusterDecoder;
  fpHWClusterDecoder=NULL;

  if (!fHistogramFile.IsNull()) {
    TFile out(fHistogramFile, "RECREATE");
    if (!out.IsZombie()) {
      out.cd();
      if (fHistoHWCFDataSize) fHistoHWCFDataSize->Write();
      if (fHistoHWCFReductionFactor) fHistoHWCFReductionFactor->Write();
      if (fHistoNofClusters) fHistoNofClusters->Write();
      if (fHistoNofClustersReductionFactor) fHistoNofClustersReductionFactor->Write();
      if (fMonitoringContainer) {
	const TObject* o1=fMonitoringContainer->FindObject("fHistograms");
	const TObject* o2=fMonitoringContainer->FindObject("fHistograms2D");
	const TObject* o3=fMonitoringContainer->FindObject("fHistograms3D");
	if (o1) o1->Write();
	if (o2) o2->Write();
	if (o3) o3->Write();
      }
      out.Close();
    }
  }
  if (fHistoHWCFDataSize) delete fHistoHWCFDataSize;
  fHistoHWCFDataSize=NULL;
  if (fHistoHWCFReductionFactor) delete fHistoHWCFReductionFactor;
  fHistoHWCFReductionFactor=NULL;
  if (fHistoNofClusters) delete fHistoNofClusters;
  fHistoNofClusters=NULL;
  if (fHistoNofClustersReductionFactor) delete fHistoNofClustersReductionFactor;
  fHistoNofClustersReductionFactor=NULL;
  if (fMonitoringContainer) {
    fMonitoringContainer->Clear();
    delete fMonitoringContainer;
  }
  fMonitoringContainer=NULL;


  return iResult;
}

int AliHLTTPCDataCompressionMonitorComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  /// inherited from AliHLTComponent: argument scan
  int iResult=0;
  if (argc<1) return 0;
  int bMissingParam=0;
  int i=0;
  TString argument=argv[i];

  do {
    // -histogram-file
    if (argument.CompareTo("-histogram-file")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      fHistogramFile=argv[i++];
      return i;
    }
    // -publishing-mode
    if (argument.CompareTo("-publishing-mode")==0) {
      if ((bMissingParam=(++i>=argc))) break;
      TString option=argv[i++];
      if (option.CompareTo("off")==0)           fPublishingMode=kPublishOff     ;
      else if (option.CompareTo("separate")==0) fPublishingMode=kPublishSeparate;
      else if (option.CompareTo("list")==0)     fPublishingMode=kPublishList    ;
      else if (option.CompareTo("array")==0)    fPublishingMode=kPublishArray   ;
      else {
	HLTError("invalid option \"%s\" for argument \"%s\", expecting 'off', 'separate', 'list', or 'array'", option.Data(), argument.Data());
	return -EPROTO;
      }
      return i;
    }
  } while (0); // using do-while only to have break available

  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EPROTO;
  }

  return iResult;
}

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AliDataContainer()
  : fHistograms(new TObjArray)  
  , fHistograms2D(new TObjArray)  
  , fHistograms3D(new TObjArray)    
  , fHistogramPointers()
  , fHistogram2DPointers()
  , fHistogram3DPointers()
  , fRemainingClusterIds()
  , fTrackModelClusterIds()
  , fCurrentClusterIds(NULL)
  , fRawData(NULL)
  , fCurrentCluster()
  , fSector(-1)
  , fBegin()
{
  /// constructor
  memset(&fCurrentCluster, 0, sizeof(AliHLTTPCRawCluster));
  if (fHistograms) {
    fHistograms->SetOwner(kTRUE);
    fHistogramPointers.resize(kNumberOfHistograms, NULL);
    for (const AliHistogramDefinition* definition=fgkHistogramDefinitions;
	 definition->fName!=NULL; definition++) {
      fHistogramPointers[definition->fId]=new TH1D(definition->fName,
						  definition->fTitle,
						  definition->fBins,
						  definition->fLowerBound,
						  definition->fUpperBound
						  );
      fHistograms->AddAt(fHistogramPointers[definition->fId], definition->fId);
    }
  }
  ///
  if (fHistograms2D) {
    fHistograms2D->SetOwner(kTRUE);
    fHistogram2DPointers.resize(kNumberOfHistograms2D, NULL);
    for (const AliHistogramDefinition2D* definition=fgkHistogramDefinitions2D;
  	 definition->fName!=NULL; definition++) {
      fHistogram2DPointers[definition->fId]=new TH2D(definition->fName,
						     definition->fTitle,
						     definition->fBinsX,
						     definition->fLowerBoundX,
						     definition->fUpperBoundX,
						     definition->fBinsY,
						     definition->fLowerBoundY,
						     definition->fUpperBoundY
						     );
      fHistograms2D->AddAt(fHistogram2DPointers[definition->fId], definition->fId);
    }
  }
  ///
  if (fHistograms3D) {
    fHistograms3D->SetOwner(kTRUE);
    fHistogram3DPointers.resize(kNumberOfHistograms3D, NULL);
    for (const AliHistogramDefinition3D* definition=fgkHistogramDefinitions3D;
  	 definition->fName!=NULL; definition++) {
      fHistogram3DPointers[definition->fId]=new TH3D(definition->fName,
						     definition->fTitle,
						     definition->fBinsX,
						     definition->fLowerBoundX,
						     definition->fUpperBoundX,
						     definition->fBinsY,
						     definition->fLowerBoundY,
						     definition->fUpperBoundY,
						     definition->fBinsZ,
						     definition->fLowerBoundZ,
						     definition->fUpperBoundZ
						     );
      fHistograms3D->AddAt(fHistogram3DPointers[definition->fId], definition->fId);
    }
  }
  
}

const AliHLTTPCDataCompressionMonitorComponent::AliHistogramDefinition AliHLTTPCDataCompressionMonitorComponent::fgkHistogramDefinitions[] = {
  {kHistogramPadrow,        "padrow"   , "padrow; padrow; counts"                  ,  159,   0.,   159.},
  {kHistogramHWCFPad,       "hwcfpad"  , "hwcfpad; pad; counts"                    ,  280,   0.,   140.},
  {kHistogramPad,           "pad"      , "pad; pad; counts"                        ,  280,   0.,   140.},
  {kHistogramTime,          "timebin"  , "timebin; time; counts"                   , 1024,   0.,  1024.},
  {kHistogramSigmaY2,       "sigmaY2"  , "sigmaY2; #sigma_{Y}^{2}; counts"         ,  100,   0.,     1.},
  {kHistogramSigmaZ2,       "sigmaZ2"  , "sigmaZ2; #sigma_{Z}^{2}; counts"         ,  100,   0.,     1.},
  {kHistogramCharge,        "charge"   , "charge; charge; counts"                  , 1024,   0., 65536.},
  {kHistogramQMax,          "qmax"     , "qmax; Q_{max}; counts"                   ,  128,   0.,  1024.},
  {kHistogramDeltaPadrow,   "d_padrow" , "d_padrow; #Delta padrow; counts"         , 1000,  -1.,     1.},
  {kHistogramDeltaPad,      "d_pad"    , "d_pad; #Delta pad; counts"               , 1000,  -1.,     1.},
  {kHistogramDeltaTime,     "d_time"   , "d_time; #Delta time; counts"             , 1000,  -1.,     1.},
  {kHistogramDeltaSigmaY2,  "d_sigmaY2", "d_sigmaY2; #Delta #sigma_{Y}^{2}; counts", 1000,  -1.,     1.},
  {kHistogramDeltaSigmaZ2,  "d_sigmaZ2", "d_sigmaZ2; #Delta #sigma_{Z}^{2}; counts", 1000,  -1.,     1.},
  {kHistogramDeltaCharge,   "d_charge" , "d_charge; #Delta charge"                 , 1000,  -1.,     1.},
  {kHistogramDeltaQMax,     "d_qmax"   , "d_qmax; #Delta Q_{max}"                  , 1000,  -1.,     1.},
  {kHistogramOutOfRange,    "ResError" , "Residual Error; padrow; counts"          ,  159,   0.,   159.},
  {kNumberOfHistograms, NULL, NULL, 0,0.,0.}
};

 const AliHLTTPCDataCompressionMonitorComponent::AliHistogramDefinition2D AliHLTTPCDataCompressionMonitorComponent::fgkHistogramDefinitions2D[] = {
   {kHistogramQMaxSector,    "qmaxsector"   , "qmaxsector; sector; Q_{max}"           , 72,0.,72., 1024,0.,1024.},
   {kHistogramSigmaY2Sector, "sigmaY2sector", "sigmaY2sector; sector; #sigma_{Y}^{2}" , 72,0.,72., 100,0.,1.},
   {kHistogramSigmaZ2Sector, "sigmaZ2sector", "sigmaZ2sector; sector; #sigma_{Z}^{2}" , 72,0.,72., 100,0.,1.},
   {kHistogramXY,            "XY", "XY; X[cm]; Y[cm]" , 100,-300.,300., 100,-300.,300.},
   {kNumberOfHistograms2D, NULL, NULL, 0,0.,0., 0,0.,0.}
 };

 const AliHLTTPCDataCompressionMonitorComponent::AliHistogramDefinition3D AliHLTTPCDataCompressionMonitorComponent::fgkHistogramDefinitions3D[] = {
   {kHistogramPadrowPadSector,"padrowpadsector","padrowpadsector; sector; pad;padrow", 72,0.,72., 140,0.,140., 159,0.,159.},
   {kNumberOfHistograms3D, NULL, NULL, 0,0.,0., 0,0.,0., 0,0.,0.}
 };

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::~AliDataContainer()
{
  /// dectructor
  if (fRawData) delete fRawData;
  if (fHistograms) delete fHistograms;
  if (fHistograms2D) delete fHistograms2D;
  if (fHistograms3D) delete fHistograms3D;
}

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::iterator& AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::BeginRemainingClusterBlock(int /*count*/, AliHLTUInt32_t specification)
{
  /// iterator of remaining clusters block of specification
  AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(specification);
  AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(specification);
  unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
  if (index<fRemainingClusterIds.size())
    fCurrentClusterIds=&fRemainingClusterIds[index];
  else
    fCurrentClusterIds=NULL;
  fBegin=iterator(this);
  return fBegin;
}

AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::iterator& AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::BeginTrackModelClusterBlock(int /*count*/)
{
  /// iterator of track model clusters
  if (fTrackModelClusterIds.fIds && fTrackModelClusterIds.fSize>0)
    fCurrentClusterIds=&fTrackModelClusterIds;
  else
    fCurrentClusterIds=NULL;
  fBegin=iterator(this);
  return fBegin;
}

int AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AddRawData(const AliHLTComponentBlockData* pDesc)
{
  /// add raw data bloack
  int iResult=0;
  if (pDesc->fDataType==AliHLTTPCDefinitions::HWClustersDataType()) {
    if (!fRawData) fRawData=new AliHLTTPCHWCFSpacePointContainer(AliHLTTPCHWCFSpacePointContainer::kModeCreateMap);
    if (!fRawData) return -ENOMEM;
    if ((iResult=fRawData->AddInputBlock(pDesc))<0) return iResult;
    AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pSpacePointGrid=AliHLTTPCHWCFSpacePointContainer::AllocateIndexGrid();
    if (pSpacePointGrid) {
      fRawData->PopulateAccessGrid(pSpacePointGrid, pDesc->fSpecification);
      fRawData->SetSpacePointPropertyGrid(pDesc->fSpecification, pSpacePointGrid);
    }
    return 0;
  }
  return -ENODATA;  
}

int AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::AddClusterIds(const AliHLTComponentBlockData* pDesc)
{
  /// add cluster id block for remaining or track model clusters
  if (!pDesc) return -EINVAL;
  if (pDesc->fDataType==AliHLTTPCDefinitions::ClusterIdTracksDataType()) {
    fTrackModelClusterIds.fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fTrackModelClusterIds.fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  if (pDesc->fDataType==AliHLTTPCDefinitions::RemainingClusterIdsDataType()) {
    AliHLTUInt8_t slice=AliHLTTPCDefinitions::GetMinSliceNr(pDesc->fSpecification);
    AliHLTUInt8_t partition=AliHLTTPCDefinitions::GetMinPatchNr(pDesc->fSpecification);
    unsigned index=slice*AliHLTTPCTransform::GetNumberOfPatches()+partition;
    if (fRemainingClusterIds.size()<=index) {
      if ((int)fRemainingClusterIds.size()<AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches()) {
	fRemainingClusterIds.resize(AliHLTTPCTransform::GetNSlice()*AliHLTTPCTransform::GetNumberOfPatches());
      } else {
	fRemainingClusterIds.resize(index+1);
      }
    }
    fRemainingClusterIds[index].fIds=reinterpret_cast<AliHLTUInt32_t*>(pDesc->fPtr);
    fRemainingClusterIds[index].fSize=pDesc->fSize/sizeof(AliHLTUInt32_t);
    return 0;
  }
  return -ENODATA;
}

AliHLTUInt32_t AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::GetClusterId(int clusterNo) const
{
  /// get the cluster id from the current cluster id block (optional)
  if (!fCurrentClusterIds ||
      clusterNo<0 ||
      (int)fCurrentClusterIds->fSize<=clusterNo)
    return kAliHLTVoidDataSpec;
  return fCurrentClusterIds->fIds[clusterNo];
}

AliHLTUInt32_t AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FindNearestCluster(int slice, int partition, const AliHLTTPCRawCluster& cluster) const
{
  /// get the cluster id of the nearest original cluster
  if (!fRawData) return kAliHLTVoidDataSpec;
  AliHLTUInt32_t key=AliHLTTPCDefinitions::EncodeDataSpecification(slice, slice, partition, partition);
  // FIXME: AliHLTIndexGrid::Index is not declared const
  AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid* pGrid=const_cast<AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid*>(fRawData->GetSpacePointPropertyGrid(key));
  if (!pGrid) return kAliHLTVoidDataSpec;
  AliHLTUInt32_t clusterId=kAliHLTVoidDataSpec;
  // search a 4x4 matrix out of the 9x9 matrix around the cell addressed by
  // pad and time
  float padrow=(float)cluster.GetPadRow()-AliHLTTPCTransform::GetFirstRow(partition);
  float pad=cluster.GetPad();
  float time=cluster.GetTime();
  float minr2=-1.;
  const float padpitch=AliHLTTPCTransform::GetPadPitchWidth(partition);
  const float zwidth=AliHLTTPCTransform::GetZWidth();
  float maxDeltaPad=AliHLTTPCDefinitions::GetMaxClusterDeltaPad();
  float maxDeltaTime=AliHLTTPCDefinitions::GetMaxClusterDeltaTime();
  int rowindex=pGrid->GetXIndex(padrow);
  int padstartindex=pGrid->GetYIndex(pad);
  int timestartindex=pGrid->GetZIndex(time);
  int cellindex=pGrid->Index(rowindex, padstartindex, timestartindex);
  float centerpad=pGrid->GetCenterY(cellindex);
  float centertime=pGrid->GetCenterZ(cellindex);
  if ((TMath::Abs(centerpad-pad)>maxDeltaPad && pad>0.) ||
      (TMath::Abs(centertime-time)>maxDeltaTime && time>0.)) {
    ALIHLTERRORGUARD(20, "invalid pad center calculation, please check dimensions if dimensions of index grid match the maximum possible deviation");
  }

  int paddirection=1;
  int timedirection=1;
  if (centerpad>pad) paddirection=-1;
  if (centertime>time) timedirection=-1;
  for (int padcount=0, padindex=padstartindex; padcount<2; padcount++, padindex+=paddirection) {
    if (padindex<0) continue;
    if (padindex>=pGrid->GetDimensionY()) break;
    for (int timecount=0, timeindex=timestartindex; timecount<2; timecount++, timeindex+=timedirection) {
      if (timeindex<0) continue;
      if (timeindex>=pGrid->GetDimensionZ()) break;
      cellindex=pGrid->Index(rowindex, padindex, timeindex);
      float cellpad=pGrid->GetCenterY(cellindex);
      float celltime=pGrid->GetCenterZ(cellindex);
      for (AliHLTSpacePointContainer::AliHLTSpacePointPropertyGrid::iterator& cl=pGrid->begin((float)padrow, cellpad, celltime);
	   cl!=pGrid->end(); cl++) {
	if (cl.Data().fTrackId>=0) continue;
	if (fRawData->GetCharge(cl.Data().fId)!=cluster.GetCharge() ||
	    fRawData->GetQMax(cl.Data().fId)!=cluster.GetQMax()) continue;
	if (TMath::Abs(padrow-fRawData->GetX(cl.Data().fId))>=1.) {
	  HLTError("slice %d, partition %d, cluster 0x%08x: mismatch on padrow: %f  vs. cluster %f", slice, partition, cl.Data().fId, padrow, fRawData->GetX(cl.Data().fId));
	  continue;
	}
	float clusterpad=fRawData->GetY(cl.Data().fId);
	float clustertime=fRawData->GetZ(cl.Data().fId);
	clusterpad-=pad;
	clusterpad*=padpitch;
	clustertime-=time;
	clustertime*=zwidth;
	float r2=clusterpad*clusterpad+clustertime*clustertime;
	if (minr2<0. || r2<minr2) {
	  clusterId=cl.Data().fId;
	  cl.Data().fTrackId=1;
	  minr2=r2;
	}
      }
    }
  }
  return clusterId;
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillPadRow(int row, int slice, AliHLTUInt32_t /*clusterId*/)
{
  /// fill padrow histogram
  unsigned index=kHistogramPadrow;
  fCurrentCluster.SetPadRow(row);
  // the inner sectors consist of readout partitions 0 and 1, if the row
  // is smaller than first row of readout partition 2, its an inner sector
  if (row<AliHLTTPCTransform::GetFirstRow(2)) {
    fSector = slice;
  } else {
    fSector = slice+AliHLTTPCTransform::GetNSlice();
  }
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(row);
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillPad(float pad, AliHLTUInt32_t /*clusterId*/)
{
  /// fill pad histogram
  fCurrentCluster.SetPad(pad);
  int currentRow=fCurrentCluster.GetPadRow();
  unsigned index=kHistogramPad;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(pad);

  index=kHistogramPadrowPadSector;
  if (index<fHistogram3DPointers.size() && fHistogram3DPointers[index]!=NULL)
    fHistogram3DPointers[index]->Fill(fSector,pad,currentRow);
  
  AliTPCROC *roc=AliTPCROC::Instance();
  Float_t pos[2]={0};
  roc->GetPositionGlobal(fSector, fSector>35?currentRow-63:currentRow, (int)pad, pos); 
  index=kHistogramXY;
  if (index<fHistogram2DPointers.size() && fHistogram2DPointers[index]!=NULL)
    fHistogram2DPointers[index]->Fill(pos[0],pos[1]);
  
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillTime(float time, AliHLTUInt32_t /*clusterId*/)
{
  /// fill pad histogram
  fCurrentCluster.SetTime(time);
  unsigned index=kHistogramTime;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(time);
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillSigmaY2(float sigmaY2, AliHLTUInt32_t /*clusterId*/, int partition)
{
  /// fill sigmaY2 histogram
  fCurrentCluster.SetSigmaY2(sigmaY2);
  unsigned index=kHistogramSigmaY2;
  /// take account for different pad widths
  float weight=AliHLTTPCTransform::GetPadPitchWidth(partition);
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(sigmaY2*weight*weight);

  index=kHistogramSigmaY2Sector;
  if (index<fHistogram2DPointers.size() && fHistogram2DPointers[index]!=NULL)
    fHistogram2DPointers[index]->Fill(fSector,sigmaY2*weight*weight);

}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillSigmaZ2(float sigmaZ2, AliHLTUInt32_t /*clusterId*/)
{
  /// fill sigmaZ2 histogram
  fCurrentCluster.SetSigmaZ2(sigmaZ2);
  unsigned index=kHistogramSigmaZ2;
  // FIXME: this is just a fixed value, to be correct the values from the global
  // parameter block has to be used
  float weight=AliHLTTPCTransform::GetZWidth();
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(sigmaZ2*weight*weight);

  index=kHistogramSigmaZ2Sector;
  if (index<fHistogram2DPointers.size() && fHistogram2DPointers[index]!=NULL)
    fHistogram2DPointers[index]->Fill(fSector,sigmaZ2*weight*weight);

}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillCharge(unsigned charge, AliHLTUInt32_t /*clusterId*/)
{
  /// fill charge histogram
  fCurrentCluster.SetCharge(charge);
  unsigned index=kHistogramCharge;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(charge);
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FillQMax(unsigned qmax, AliHLTUInt32_t /*clusterId*/)
{
  /// fill qmax histogram
  fCurrentCluster.SetQMax(qmax);
  unsigned index=kHistogramQMax;
  if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL)
    fHistogramPointers[index]->Fill(qmax);

  index=kHistogramQMaxSector;
  if (index<fHistogram2DPointers.size() && fHistogram2DPointers[index]!=NULL)
    fHistogram2DPointers[index]->Fill(fSector,qmax);
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::Fill(int slice, int partition, AliHLTUInt32_t clusterId)
{
  /// fill cluster histograms requiring the full cluster information
  
  // TODO: the complete filling of histograms can be moved to this function
  // and the cluster struct be filled in the iterator
  // The delta histograms are filled here either by using the specified
  // cluster, or the nearest cluster on the padrow with identical charge
  // and qmax is searched for comparison.
  if (clusterId==kAliHLTVoidDataSpec) {
    clusterId=FindNearestCluster(slice, partition, fCurrentCluster);
  }
  if (clusterId==kAliHLTVoidDataSpec) return;
  bool bResidualError=false;
  int currentRow=fCurrentCluster.GetPadRow();

    unsigned index=kHistogramDeltaPadrow;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(fCurrentCluster.GetPadRow()-fRawData->GetX(clusterId));
      }
    }

    index=kHistogramDeltaPad;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	float dPad=fCurrentCluster.GetPad()-fRawData->GetY(clusterId);
	fHistogramPointers[index]->Fill(dPad);
	static const float maxdPad=0.015; // better 100um for 4 and 6mm pad width
	if (TMath::Abs(dPad)>maxdPad) {
	  //HLTError("cluster 0x%08x slice %d partition %d: pad difference %f - max %f", clusterId, slice, partition, dPad, maxdPad);
	  bResidualError=true;
	}
      }
    }

    index=kHistogramDeltaTime;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	float dTime=fCurrentCluster.GetTime()-fRawData->GetZ(clusterId);
	fHistogramPointers[index]->Fill(dTime);
	static const float maxdTime=0.04; // corresponds to 100um
	if (TMath::Abs(dTime)>maxdTime) {
	  //HLTError("cluster 0x%08x slice %d partition %d: time difference %f - max %f", clusterId, slice, partition, dTime, maxdTime);
	  bResidualError=true;
	}
      }
    }

    index=kHistogramDeltaSigmaY2;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(fCurrentCluster.GetSigmaY2()-fRawData->GetYWidth(clusterId));
      }
    }

    index=kHistogramDeltaSigmaZ2;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(fCurrentCluster.GetSigmaZ2()-fRawData->GetZWidth(clusterId));
      }
    }

    index=kHistogramDeltaCharge;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(fCurrentCluster.GetCharge()-fRawData->GetCharge(clusterId));
      }
    }

    index=kHistogramDeltaQMax;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      if (fRawData->Check(clusterId)) {
	fHistogramPointers[index]->Fill(fCurrentCluster.GetQMax()-fRawData->GetQMax(clusterId));
      }
    }

    if (bResidualError) {
    index=kHistogramOutOfRange;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData) {
      fHistogramPointers[index]->Fill(currentRow>=0?currentRow:0);
    }
    }

    index=kHistogramHWCFPad;
    if (index<fHistogramPointers.size() && fHistogramPointers[index]!=NULL && fRawData)
      fHistogramPointers[index]->Fill(fRawData->GetY(clusterId));
}

void AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::Clear(Option_t * option)
{
  /// internal cleanup
  if (fRawData) fRawData->Clear(option);
}

TObject* AliHLTTPCDataCompressionMonitorComponent::AliDataContainer::FindObject(const char *name) const
{
  /// get histogram object  
  if (!name) return NULL;
  if ( strcmp(name,"fHistograms")   == 0 )
    return fHistograms;
  if ( strcmp(name,"fHistograms2D") == 0 )
    return fHistograms2D;
  if ( strcmp(name,"fHistograms3D") == 0 )
    return fHistograms3D;

  return NULL;
}
