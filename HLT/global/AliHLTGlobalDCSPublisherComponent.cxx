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

/** @file   AliHLTGlobalDCSPublisherComponent.cxx
    @author Matthias Richter
    @date   20010-03-10
    @brief  DIM publisher component for global HLT data
*/

#include "AliHLTGlobalDCSPublisherComponent.h"
#include "AliHLTDimServer.h"
#include "AliESDVertex.h"
#include "TDatime.h"
#include "TMath.h"
#include "TSystem.h"
#include <algorithm>
#include <numeric>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTGlobalDCSPublisherComponent)

AliHLTGlobalDCSPublisherComponent::AliHLTGlobalDCSPublisherComponent()
  : AliHLTDataSink()
  , fServerName()
  , fDimdns()
  , fpServer(NULL)
  , fPosition(0)
  , fEventBufferSize(1000)
  , fLastUpdate(0)
  , fUpdatePeriod(10)
  , fEventBuffers()
  , fDimServices()
  , fDimServiceEventCount(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  for (int i=0; i<kLastService; i++) fDimServices[i]=NULL;
}

const char* AliHLTGlobalDCSPublisherComponent::fgkServiceNames[kLastService]= {
  "Vertex_X",
  "Vertex_Y",
  "Vertex_Z",
  "RmsVertex_X",
  "RmsVertex_Y",
  "RmsVertex_Z",
};

AliHLTGlobalDCSPublisherComponent::~AliHLTGlobalDCSPublisherComponent()
{
  // see header file for class documentation

  // file list and file name list are owner of their objects and
  // delete all the objects
}

const char* AliHLTGlobalDCSPublisherComponent::GetComponentID()
{
  // see header file for class documentation
  return "DCSPublisher";
}

void AliHLTGlobalDCSPublisherComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponent* AliHLTGlobalDCSPublisherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalDCSPublisherComponent;
}

int AliHLTGlobalDCSPublisherComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  iResult=ConfigureFromArgumentString(argc, argv);

  fpServer=new AliHLTDimServer(fServerName.c_str());
  if (!fpServer) return -ENOMEM;
  if ((iResult=fpServer->Init(fDimdns.c_str()))>=0) {
    // add services
    for (int service=0; service<kLastService; service++) {
      fDimServices[service]=fpServer->CreateService(AliHLTDimServer::kDataTypeFloat, fgkServiceNames[service]);
    }
    fDimServiceEventCount=fpServer->CreateService(AliHLTDimServer::kDataTypeInt, "EventCount");

    // now start the server
    iResult=fpServer->Start();
  }

  fPosition=0;
  for (int i=0; i<kDimensions; i++)
    Reset(fEventBuffers[i], fEventBufferSize);

  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation

  // -servername
  if (argc==0) return 0;
  int i=0;
  TString argument=argv[0];
  if (argument.CompareTo("-servername")==0) {
    if (++i>=argc) return -EPROTO;
    fServerName=argv[i];
    return 2;
  }
  
  // --dimdns
  if (argument.CompareTo("-dimdns")==0) {
    if (++i>=argc) return -EPROTO;
    fDimdns=argv[i];
    return 2;
  }

  // TODO: further options for tracklet cut and event buffer size

  return -EINVAL;
}

int AliHLTGlobalDCSPublisherComponent::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  Publish(true);
  // just wait a few seconds to get the updated services through before
  // terminating the server
  gSystem->Sleep(fUpdatePeriod<30?fUpdatePeriod:30);

  if (!fpServer) return -ENODEV;
  fpServer->Stop();
  delete fpServer;
  fpServer=NULL;

  // TODO: proper cleanup of services objects to be synchronized with
  // the server

  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::DumpEvent( const AliHLTComponentEventData& /*evtData*/,
						  AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  if (!IsDataEvent()) return 0;

  int iResult=0;
  for (const TObject* pObject=GetFirstInputObject(kAliHLTDataTypeESDVertex);
       pObject!=NULL;
       pObject=GetNextInputObject()) {
    const AliESDVertex* pVertex=dynamic_cast<const AliESDVertex*>(pObject);
    if (pVertex && pVertex->GetNContributors()>=5) {
      (fEventBuffers[kVertexX])[fPosition%fEventBuffers[kVertexX].size()]=pVertex->GetX();
      (fEventBuffers[kVertexY])[fPosition%fEventBuffers[kVertexY].size()]=pVertex->GetY();
      (fEventBuffers[kVertexZ])[fPosition%fEventBuffers[kVertexZ].size()]=pVertex->GetZ();
      fPosition++;
    }
  }
  
  TDatime time;
  if ((time.Get()-fLastUpdate>fUpdatePeriod) && fPosition>0) {
    fLastUpdate=time.Get();
    Publish();
  }
  return iResult;
}

int AliHLTGlobalDCSPublisherComponent::Publish(bool reset)
{
  // publish the values to DCS
  float mean[kDimensions];
  float rms[kDimensions];
  for (int dimension=0; dimension<kDimensions; dimension++) {
    mean[dimension]= Mean(fEventBuffers[dimension], fPosition);
    rms[dimension] = Rms(fEventBuffers[dimension], mean[dimension], fPosition);
    ((AliHLTDimServer::AliHLTDimServiceFloat*)fDimServices[dimension])->Update(mean[dimension]);
    ((AliHLTDimServer::AliHLTDimServiceFloat*)fDimServices[dimension+kDimensions])->Update(rms[dimension]);
  }
  ((AliHLTDimServer::AliHLTDimServiceInt*)fDimServiceEventCount)->Update(fPosition<fEventBufferSize?fPosition:fEventBufferSize);

  HLTInfo("Vertex from %d samples: X:Y:Z Mean %f:%f:%f  RMS %f:%f:%f", 
	  fPosition<fEventBufferSize?fPosition:fEventBufferSize, 
	  mean[kVertexX], mean[kVertexY], mean[kVertexZ],
	  rms[kVertexX], rms[kVertexY], rms[kVertexZ]
	  );

  if (reset) {
    for (int i=0; i<kDimensions; i++)
      Reset(fEventBuffers[i], fEventBufferSize);
    fPosition=0;
  }

  return 0;
}

void AliHLTGlobalDCSPublisherComponent::Reset(vector<float>& sample, int size) const
{
  // reset the event buffer
  if (size>0) sample.resize(size);
  fill(sample.begin(), sample.end(), 0.0);
}

template <class T> T AliHLTGlobalDCSPublisherComponent::Mean(const vector<T>& sample, int count) const
{
  // calculate mean of the event buffer
  int samplesize=count<(int)sample.size()?count:sample.size();
  if (samplesize==0) return 0.0;
  T sum=std::accumulate(sample.begin(), sample.begin()+samplesize, (T)0, std::plus<T>());
  return sum/samplesize;
}

template <class T> T AliHLTGlobalDCSPublisherComponent::Rms(const vector<T>& sample, T mean, int count) const
{
  // calculate sigma of the event buffer
  int samplesize=count<(int)sample.size()?count:sample.size();
  if (samplesize==0) return 0.0;
  T msd=std::accumulate(sample.begin(), sample.begin()+samplesize, (T)0, MeanSqtOp<T>(mean))/samplesize;
  if (msd<0.00000001) return 0.0;
  return TMath::Sqrt(msd);
}
