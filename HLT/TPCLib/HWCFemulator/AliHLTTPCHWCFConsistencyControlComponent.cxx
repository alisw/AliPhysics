//****************************************************************************
//* This file is property of and copyright by the ALICE HLT Project          * 
//* ALICE Experiment at CERN, All rights reserved.                           *
//*                                                                          *
//* Primary Authors: Sergey Gorbunov, Torsten Alt                            *
//* Developers:      Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de> *
//*                  Torsten Alt <talt@cern.ch>                              *
//*                  for The ALICE HLT Project.                              *
//*                                                                          *
//* Permission to use, copy, modify and distribute this software and its     *
//* documentation strictly for non-commercial purposes is hereby granted     *
//* without fee, provided that the above copyright notice appears in all     *
//* copies and that both the copyright notice and this permission notice     *
//* appear in the supporting documentation. The authors make no claims       *
//* about the suitability of this software for any purpose. It is            *
//* provided "as is" without express or implied warranty.                    *
//****************************************************************************

//  @file   AliHLTTPCHWCFConsistencyControlComponent.cxx
//  @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
//  @author Torsten Alt <talt@cern.ch> 
//  @brief  Comparison of TPC clusters produced by FPGA clusterfinder and by FPGA Emulator 
//  @brief  ( see AliHLTTPCHWCFEmulator class )
//  @note


#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCHWCFConsistencyControlComponent.h"

#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCHWCFDataTypes.h"
#include "AliHLTTPCClusterMCData.h"

#include "AliGRPObject.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliRawDataHeader.h"
#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1F.h"

#include <sys/time.h>
#include "TFile.h"

AliHLTTPCHWCFConsistencyControlComponent::AliHLTTPCHWCFConsistencyControlComponent()
  :
  AliHLTProcessor(),
  fNDismatch(0),
  fNBlocks(0),
  fBenchmark("TPCHWConsistencyControl"), 
  fHistHeaderAll(0),
  fHistHeaderGood(0),
  fHistClusterAll(0),
  fHistClusterGood(0),
  fProfHeader(0),
  fProfCluster(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}


AliHLTTPCHWCFConsistencyControlComponent::AliHLTTPCHWCFConsistencyControlComponent(const AliHLTTPCHWCFConsistencyControlComponent&)
  :
  AliHLTProcessor(),
  fNDismatch(0),
  fNBlocks(0),
  fBenchmark("TPCHWConsistencyControl"),
  fHistHeaderAll(0),
  fHistHeaderGood(0),
  fHistClusterAll(0),
  fHistClusterGood(0),
  fProfHeader(0),
  fProfCluster(0)
{
  // dummy
}

AliHLTTPCHWCFConsistencyControlComponent& AliHLTTPCHWCFConsistencyControlComponent::operator=(const AliHLTTPCHWCFConsistencyControlComponent&)
{
  // dummy
  return *this;
}

AliHLTTPCHWCFConsistencyControlComponent::~AliHLTTPCHWCFConsistencyControlComponent()
{
  // see header file for class documentation
  DoDeinit();
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCHWCFConsistencyControlComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCHWCFConsistenyControl";
}

void AliHLTTPCHWCFConsistencyControlComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkHWClustersDataType| kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCHWCFConsistencyControlComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram  | kAliHLTDataOriginTPC;
}


void AliHLTTPCHWCFConsistencyControlComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 1000;
  inputMultiplier = 0;
}


AliHLTComponent* AliHLTTPCHWCFConsistencyControlComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCHWCFConsistencyControlComponent();
}

void AliHLTTPCHWCFConsistencyControlComponent::GetOCDBObjectDescription( TMap* const /*targetMap*/){
// Get a list of OCDB object description needed for the particular component
  /*
  if (!targetMap) return;
  
  // OCDB entries for component arguments
  targetMap->Add(new TObjString("HLT/ConfigTPC/TPCHWClusterFinder"), new TObjString("component arguments, empty at the moment"));
  */
}


int AliHLTTPCHWCFConsistencyControlComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation

  fHistHeaderAll = new TH1F("hHWCFHeaderAll", "fHistHeaderAll",6,0.,6.);
  fHistHeaderGood = new TH1F("hHWCFHeaderGood", "fHistHeaderGood",6,0.,6.);
  fHistClusterAll = new TH1F("hHWCFClusterAll", "fHistClusterAll",8,0.,8.);
  fHistClusterGood = new TH1F("hHWCFClusterGood", "fHistClusterGood",8,0.,8.);
  fProfHeader = new TH1F("pHWCFHeader", "HWCF: Consistency of header data", 6, 0., 6.);
  fProfCluster = new TH1F("pHWCFClusters", "HWCF: Consisteny of cluster data", 8, 0., 8.);

  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }

  return Configure( NULL, NULL, arguments.Data()  );
}

int AliHLTTPCHWCFConsistencyControlComponent::Reconfigure( const char* cdbEntry, const char* chainId )
{
  // Reconfigure the component from OCDB

  return Configure( cdbEntry, chainId, NULL );
}

int AliHLTTPCHWCFConsistencyControlComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  TString arguments = "";
  for ( int i = 0; i < argc; i++ ) {
    if ( !arguments.IsNull() ) arguments += " ";
    arguments += argv[i];
  }
  return ReadConfigurationString(arguments);
}



void AliHLTTPCHWCFConsistencyControlComponent::SetDefaultConfiguration()
{
  // Set default configuration for the FPGA ClusterFinder Emulator component
  // Some parameters can be later overwritten from the OCDB

  fNDismatch = 0;
  fNBlocks = 0;
  fBenchmark.Reset();
  fBenchmark.SetTimer(0,"total");
  fBenchmark.SetTimer(1,"reco");    
}

int AliHLTTPCHWCFConsistencyControlComponent::ReadConfigurationString(  const char* arguments )
{
  // Set configuration parameters for the FPGA ClusterFinder Emulator component
  // from the string

  int iResult = 0;
  if ( !arguments ) return iResult;

  TString allArgs = arguments;
  TString argument;
  int bMissingParam = 0;

  TObjArray* pTokens = allArgs.Tokenize( " " );

  int nArgs =  pTokens ? pTokens->GetEntries() : 0;

  for ( int i = 0; i < nArgs; i++ ) {
    argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
    if ( argument.IsNull() ) continue;

    //if ( argument.CompareTo( "-deconvolute-time" ) == 0 ) {
    //if ( ( bMissingParam = ( ++i >= pTokens->GetEntries() ) ) ) break;
    //fDoDeconvTime  = ( ( TObjString* )pTokens->At( i ) )->GetString().Atoi();
    //HLTInfo( "Time deconvolution is set to: %d", fDoDeconvTime );
    //continue;
    //}

    HLTError( "Unknown option \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }
  delete pTokens;

  if ( bMissingParam ) {
    HLTError( "Specifier missed for parameter \"%s\"", argument.Data() );
    iResult = -EINVAL;
  }

  return iResult;
}


int AliHLTTPCHWCFConsistencyControlComponent::ReadCDBEntry( const char* /*cdbEntry*/, const char* /*chainId*/ )
{
  // Read configuration from OCDB
  return 0;
  /*
  const char* defaultNotify = "";

  if ( !cdbEntry ){
    cdbEntry = "HLT/ConfigTPC/TPCHWClusterFinder";
    defaultNotify = " (default)";
    chainId = 0;
  }

  HLTInfo( "configure from entry \"%s\"%s, chain id %s", cdbEntry, defaultNotify, ( chainId != NULL && chainId[0] != 0 ) ? chainId : "<none>" );
  AliCDBEntry *pEntry = AliCDBManager::Instance()->Get( cdbEntry );//,GetRunNo());
  
  if ( !pEntry ) {
    HLTError( "cannot fetch object \"%s\" from CDB", cdbEntry );
    return -EINVAL;
  }

  TObjString* pString = dynamic_cast<TObjString*>( pEntry->GetObject() );

  if ( !pString ) {
    HLTError( "configuration object \"%s\" has wrong type, required TObjString", cdbEntry );
    return -EINVAL;
  }

  HLTInfo( "received configuration object string: \"%s\"", pString->GetString().Data() );

  return  ReadConfigurationString( pString->GetString().Data() );
  */
}


int AliHLTTPCHWCFConsistencyControlComponent::Configure( const char* cdbEntry, const char* chainId, const char *commandLine )
{
  // Configure the component
  // There are few levels of configuration,
  // parameters which are set on one step can be overwritten on the next step

  //* read hard-coded values

  SetDefaultConfiguration();

  //* read the default CDB entry

  int iResult1 = ReadCDBEntry( NULL, chainId );

  //* read the actual CDB entry if required

  int iResult2 = ( cdbEntry ) ? ReadCDBEntry( cdbEntry, chainId ) : 0;

  //* read extra parameters from input (if they are)

  int iResult3 = 0;

  if ( commandLine && commandLine[0] != '\0' ) {
    HLTInfo( "received configuration string from HLT framework: \"%s\"", commandLine );
    iResult3 = ReadConfigurationString( commandLine );
  }

  return iResult1 ? iResult1 : ( iResult2 ? iResult2 : iResult3 );
}


int AliHLTTPCHWCFConsistencyControlComponent::DoDeinit()
{
  // see header file for class documentation 
  if( fHistHeaderAll ) delete fHistHeaderAll;
  if( fHistHeaderGood ) delete fHistHeaderGood;
  if( fHistClusterAll ) delete fHistClusterAll;
  if( fHistClusterGood ) delete fHistClusterGood;
  if( fProfHeader ) delete fProfHeader;
  if( fProfCluster ) delete fProfCluster;

  fHistHeaderAll = 0;
  fHistHeaderGood = 0;
  fHistClusterAll = 0;
  fHistClusterGood = 0;
  fProfHeader = 0;
  fProfCluster = 0;
  return 0;
}


int AliHLTTPCHWCFConsistencyControlComponent::DoEvent( const AliHLTComponentEventData& evtData, 
						       const AliHLTComponentBlockData* blocks, 
						       AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, 
						       AliHLTUInt32_t& size, 
						       vector<AliHLTComponentBlockData>& /*outputBlocks*/ )
{
  // see header file for class documentation

  int iResult=0;
  //AliHLTUInt32_t maxSize = size;
  size = 0;
  
  if(!IsDataEvent()){
    return 0;
  }

  fBenchmark.StartNewEvent();
  fBenchmark.Start(0);
  bool *checkedBlocks = new bool[evtData.fBlockCnt];
  for( unsigned long i=0; i<evtData.fBlockCnt; i++ ) checkedBlocks[i] = 0;
  
  for ( unsigned long ndx1 = 0; ndx1 < evtData.fBlockCnt; ndx1++ ){
    const AliHLTComponentBlockData* iter1 = blocks+ndx1;
    fBenchmark.AddInput(iter1->fSize);
    
    if (  iter1->fDataType != (AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC) 
	  ) continue;
 
    if( checkedBlocks[ndx1] ) continue;// block already checked

    checkedBlocks[ndx1] = 1;

    int slice1 = AliHLTTPCDefinitions::GetMinSliceNr( *iter1 );
    int patch1 = AliHLTTPCDefinitions::GetMinPatchNr( *iter1 );
   
    bool header1OK = 1;
    if( iter1->fSize>0 ){
      header1OK = ( iter1->fSize % sizeof(AliHLTUInt32_t) == 0 ) 
	&& ( iter1->fSize >= sizeof(AliRawDataHeader) ) 
	&& (iter1->fPtr != NULL );      
    }
    
    int  nMatchedBlocks=0;
        
    bool sameHLTHd = 1;     
    bool sameCDH = 1;     
    bool sameRCU = 1;      
    bool sameSize = 1;
   
    bool sameFlag = 1;
    bool sameQMax = 1;
    bool sameCharge = 1;
    bool sameRow = 1;
    bool sameFloat[4] = {1,1,1,1};
    
    for ( unsigned long ndx2 = ndx1+1; ndx2 < evtData.fBlockCnt; ndx2++ ){      
      const AliHLTComponentBlockData* iter2 = blocks+ndx2;
      if (  iter2->fDataType != (AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC)
	    ) continue;
            
      int slice2 = AliHLTTPCDefinitions::GetMinSliceNr( *iter2 );
      int patch2 = AliHLTTPCDefinitions::GetMinPatchNr( *iter2 );
      if( slice1!=slice2 || patch1!=patch2 ) continue;
      
      if( checkedBlocks[ndx2] ) continue;
      checkedBlocks[ndx2] = 1;

      nMatchedBlocks++;

      bool header2OK = 1;
      if( iter2->fSize>0 ){
	header2OK = ( iter2->fSize % sizeof(AliHLTUInt32_t) == 0 ) 
	  && ( iter2->fSize >= sizeof(AliRawDataHeader) ) 
	  && (iter2->fPtr != NULL );
      }
 
      fHistHeaderAll->Fill(2);
      sameHLTHd = ( header1OK == header2OK );
      if( sameHLTHd ) fHistHeaderGood->Fill(2);
      if( !header1OK || !header2OK ) continue;
      
      int nWordsHeader =  sizeof(AliRawDataHeader)/sizeof(AliHLTUInt32_t);
      int nWords1 = iter1->fSize/sizeof(AliHLTUInt32_t);
      int nWords2 = iter2->fSize/sizeof(AliHLTUInt32_t);
      const AliHLTUInt32_t *p1 = (const AliHLTUInt32_t *) iter1->fPtr;
      const AliHLTUInt32_t *p2 = (const AliHLTUInt32_t *) iter2->fPtr;

      // compare CDH headers
      
       for(AliHLTInt32_t i=0; i<nWordsHeader; i++ ){
	if( p1[i] != p2[i] ) sameCDH = 0;
      }

      fHistHeaderAll->Fill(3);
      if( sameCDH ) fHistHeaderGood->Fill(3);

      // find rcu headers
      
      int startRCU1 = nWordsHeader;
      int startRCU2 = nWordsHeader;

      for( ; startRCU1<nWords1; startRCU1+=6 ){
	if( p1[startRCU1]>>30 == 0x2 ) break;	  
      }

      for( ; startRCU2<nWords2; startRCU2+=6 ){
	if( p2[startRCU2]>>30 == 0x2 ) break;	  
      } 
      

      // compare RCU headers

      if( startRCU1 < nWords1 || startRCU2 < nWords2 ){
	if( nWords1 - startRCU1 != nWords2 - startRCU2 ) sameRCU = 0;
	else{
	  for( AliHLTInt32_t i1=startRCU1, i2=startRCU2; (i1<nWords1) && (i2<nWords2); i1++,i2++ ){
	    if( p1[i1]!=p2[i2] ) sameRCU = 0;
	  }
	}
	fHistHeaderAll->Fill(4);
	if( sameRCU ) fHistHeaderGood->Fill(4);
      }
      
      sameSize = ( startRCU1 == startRCU2 );

      fHistHeaderAll->Fill(5);
      if( sameSize ) fHistHeaderGood->Fill(5);

      // compare clusters

      if( startRCU1 == startRCU2 ){
 	for( AliHLTInt32_t i=nWordsHeader; i<startRCU1; i+=6){
	  AliHLTUInt32_t header1 = p1[i];
	  AliHLTUInt32_t charge1 = p1[i+1];
	  AliHLTUInt32_t header2 = p2[i];
	  AliHLTUInt32_t charge2 = p2[i+1];


	  AliHLTUInt32_t flag1 = header1 >> 30;
	  AliHLTUInt32_t qmax1 = ( header1 & 0xFFFFFF ) >> 6;
	  AliHLTUInt32_t row1 = (header1>>24) & 0x3f;
	  AliHLTUInt32_t flag2 = header2 >> 30;
	  AliHLTUInt32_t qmax2 = ( header2 & 0xFFFFFF ) >> 6;
	  AliHLTUInt32_t row2 = (header2>>24) & 0x3f;
	  
	  fHistClusterAll->Fill(0);
	  if( flag1 == flag2 ) fHistClusterGood->Fill(0);
	  else sameFlag = 0;
	  
	  if( flag1!=0x3 || flag2!=0x3 ) continue;
	  
	  // compare cluster qMax
	  
	  fHistClusterAll->Fill(1);
	  if( qmax1 == qmax2 ) fHistClusterGood->Fill(1);
	  else sameQMax = 0;
	  
	  // compare cluster row index
	  
	  fHistClusterAll->Fill(2);
	  if( row1 == row2 ) fHistClusterGood->Fill(2);
	  else sameRow = 0;
	  
	  // compare cluster charge
	  
	  fHistClusterAll->Fill(3);
	  if( charge1 == charge2 ) fHistClusterGood->Fill(3);
	  else sameCharge = 0;

	  // compare floating point data

	  for( int j=0; j<4; j++ ){
	    AliHLTFloat64_t f1     = *((AliHLTFloat32_t*)&p1[i+2+j]);
	    AliHLTFloat64_t f2     = *((AliHLTFloat32_t*)&p2[i+2+j]);
	    double w = (fabs(f1) + fabs(f2))/2;
	    if( w>1.e-20 ){
	      fHistClusterAll->Fill(4+j);
	      if( fabs(f1 - f2) < 1.e-6*w ) fHistClusterGood->Fill(4+j);
	      else sameFloat[j] = 0;
	    }
	  }
	}
      }
    }
        
    fNBlocks++;
      
    bool err = 0;
    TString warn;
    warn.Form("HWCF consistency check for slice %d patch %d. Values  not matched:", slice1, patch1 );
    
    if( nMatchedBlocks!=1 ){ 
      err=1; 
      TString x;
      x.Form(" NMatchedBlocks(%d)",nMatchedBlocks);
      warn+=x;
    }
    if( !sameHLTHd ){ err=1; warn+=", HLT header";}
    if( !sameCDH ){ err=1; warn+=", CDH header";}
    if( !sameRCU ){ err=1; warn+=", RCU header";}
    if( !sameSize ){ err=1; warn+=", N Clusters";}
    if( !sameFlag ){ err=1; warn+=", Cluster Header";}
    if( !sameQMax ){ err=1; warn+=", Cluster Qmax";}
    if( !sameCharge ){ err=1; warn+=", Cluster Charge";}
    if( !sameRow ){ err=1; warn+=", Cluster Row";}
    if( !sameFloat[0] ){ err=1; warn+=", Cluster Pad value";}
    if( !sameFloat[1] ){ err=1; warn+=", Cluster Time value";}
    if( !sameFloat[2] ){ err=1; warn+=", Cluster PadErr value";}
    if( !sameFloat[3] ){ err=1; warn+=", Cluster TimeErr value";}
      
    if( err ){
      fNDismatch++;		
      HLTWarning(warn.Data());
    } else {
      //warn+=" NO ";
      //HLTWarning(warn.Data());     
    }
       
    fHistHeaderAll->Fill(0);
    if( nMatchedBlocks >=1 ) fHistHeaderGood->Fill(0);
    fHistHeaderAll->Fill(1);
    if( nMatchedBlocks <=1 ) fHistHeaderGood->Fill(1);

  } // first block
    
  
  delete[] checkedBlocks;
   
  HLTInfo("HWCF consistency check: %.10f %s of %ld data blocks are OK",
	     (double)(fNBlocks-fNDismatch)/(double)fNBlocks*100.,"%",fNBlocks);
  
  //if( fNDismatch>0 ){
  //HLTWarning("HWCF inconsistency: %ld of %ld data blocks are not OK",
  //fNDismatch,fNBlocks);      
  //}  

  fProfHeader->Divide(fHistHeaderGood, fHistHeaderAll,1,1.,"b");
  fProfCluster->Divide(fHistClusterGood, fHistClusterAll, 1, 1, "b");

  PushBack( (TObject*) fProfHeader, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC,0);
  fBenchmark.AddOutput(GetLastObjectSize());
  PushBack( (TObject*) fProfCluster, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC,0);
  fBenchmark.AddOutput(GetLastObjectSize());
  

  fBenchmark.Stop(0);  
  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}


