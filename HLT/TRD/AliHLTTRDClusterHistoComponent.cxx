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

#include "AliHLTTRDClusterHistoComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDCluster.h"
#include "AliTRDcluster.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include <TFile.h>
#include <TString.h>
#include "TObjString.h"
#include "TObjArray.h"

//#include "AliHLTTRD.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDClusterHistoComponent)

AliHLTTRDClusterHistoComponent::AliHLTTRDClusterHistoComponent()
: fNClsDet(0),
  fClsAmp(0),
  fClsAmpDrift(0),
  fClsTB(0),
  fClsAmpDist(0)
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
  return kAliHLTDataTypeHistogram;

}

void AliHLTTRDClusterHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTRDClusterHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTRDClusterHistoComponent;
}

int AliHLTTRDClusterHistoComponent::DoInit( int argc, const char** argv )
{
  // Initialize histograms

  fNClsDet = new TH1D("trdClsDet", ";detector", 540, -0.5, 539.5);
  fClsAmp  = new TH1D("trdClsAmp", ";amplitude", 200, -0.5, 199.5);
  fClsAmpDrift = new TH1D("trdClsAmpDrift", ";amplitude", 200, -0.5, 199.5) ;
  fClsTB = new TH1D("trdClsTB", ";time bin", 35, -0.5, 34.5);
  fClsAmpDist = new TH1D("trdClsAmpDist", "mean amplitude", 200, 0, 100);

  for(int i=0; i<540; i++)
    fClsAmpDriftDet[i] = new TH1D(Form("trdClsDriftDet_%d",i), "", 200, -0.5, 199.5);

  /*
  // configure

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

  return iResult; 
    */
    
    return 0;
}
  
int AliHLTTRDClusterHistoComponent::DoDeinit()
{
  // see header file for class documentation

  // delete histograms
  if (fNClsDet) delete fNClsDet;
  if (fClsAmp) delete fClsAmp;
  if (fClsAmpDrift) delete fClsAmpDrift;
  if (fClsTB) delete fClsTB;
  if (fClsAmpDist) delete fClsAmpDist;

  for(int i=0; i<540; i++)
    if (fClsAmpDriftDet[i]) delete fClsAmpDriftDet[i];

  return 0;
}

int AliHLTTRDClusterHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, 
					    AliHLTComponentTriggerData& /*trigData*/)
{

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  
  const AliHLTComponentBlockData* iter = NULL;
  AliTRDcluster *cls = new AliTRDcluster();
  
  for ( iter = GetFirstInputBlock(AliHLTTRDDefinitions::fgkClusterDataType); 
	iter != NULL; iter = GetNextInputBlock() ) {
    
    AliHLTTRDCluster* clusterData = (AliHLTTRDCluster*) iter->fPtr;
    Int_t nSpacepoint = 0; 
        
    // loop over clusters 
    for(int i=0;i<nSpacepoint;i++) {
      
      clusterData[i].ExportTRDCluster(cls);
      
      fNClsDet->Fill(cls->GetDetector());
      fClsAmp->Fill(cls->GetQ());
      
      int tb = cls->GetPadTime();
      fClsTB->Fill(tb);
      if (tb > 5 && tb <25)
	fClsAmpDrift->Fill(cls->GetQ()); 
      
      fClsAmpDriftDet[cls->GetDetector()]->Fill(cls->GetQ());
    }
  }
    
  delete cls;
   
  fClsAmpDist->Reset();
  for(int det=0; det<540; det++)
    if (fClsAmpDriftDet[det]->GetSum() > 0) 
      fClsAmpDist->Fill(fClsAmpDriftDet[det]->GetMean());


  //delete til dodeinit
  // if(fPlotChargeOROCAll){
  //  AliHLTUInt32_t fSpecification = AliHLTTRDDefinitions::EncodeDataSpecification(0,35,2,5);
  //  PushBack( (TObject*) fTotalClusterChargeOROCAll,kAliHLTDataTypeHistogram,fSpecification);
  //}
  
  
  return 0;
}

int AliHLTTRDClusterHistoComponent::Configure(const char* arguments)
{
 
  return 0;
}

int AliHLTTRDClusterHistoComponent::Reconfigure(const char* cdbEntry, const char* chainId)
{


  return 0;
}
