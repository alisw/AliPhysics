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


#include <sys/time.h>
#include "TFile.h"

AliHLTTPCHWCFConsistencyControlComponent::AliHLTTPCHWCFConsistencyControlComponent()
  :
  AliHLTProcessor(),
  fNDismatch(0),
  fNBlocks(0),
  fBenchmark("TPCHWConsistencyControl")
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
  fBenchmark("TPCHWConsistencyControl")
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
  return kAliHLTVoidDataType;//kMultipleDataType;
}

int AliHLTTPCHWCFConsistencyControlComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList)
{
  // see header file for class documentation
  tgtList.clear();
  return tgtList.size();
}


void AliHLTTPCHWCFConsistencyControlComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  inputMultiplier = 6;
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
    
    int nSecondBlocksFound=0;    
    bool sizeInconsistency = 0;	
    AliHLTInt64_t intDiff = 0;
    AliHLTFloat64_t floatDiff = 0;
    
    for ( unsigned long ndx2 = ndx1+1; ndx2 < evtData.fBlockCnt; ndx2++ ){
      
      const AliHLTComponentBlockData* iter2 = blocks+ndx2;      
      if (  iter2->fDataType != (AliHLTTPCDefinitions::fgkHWClustersDataType | kAliHLTDataOriginTPC)
	    ) continue;
            
      int slice2 = AliHLTTPCDefinitions::GetMinSliceNr( *iter2 );
      int patch2 = AliHLTTPCDefinitions::GetMinPatchNr( *iter2 );
      if( slice1!=slice2 || patch1!=patch2 ) continue;

      if( checkedBlocks[ndx2] ) continue;

      checkedBlocks[ndx2] = 1;

      nSecondBlocksFound++;
      
      int nWordsHeader =  sizeof(AliRawDataHeader)/4;
      int nWords = iter1->fSize/4;
 
      const AliHLTUInt32_t *p1 = (const AliHLTUInt32_t *) iter1->fPtr;
      const AliHLTUInt32_t *p2 = (const AliHLTUInt32_t *) iter2->fPtr;

      // compare size 

      bool sizeInc = 
	( iter1->fSize != iter2->fSize )
	|| ( iter1->fSize % 4 != 0 ) || ( iter2->fSize % 4 !=0 )
	|| iter1->fSize < sizeof(AliRawDataHeader)      
	|| iter2->fSize < sizeof(AliRawDataHeader)
	|| (nWords>0 && (!iter1->fPtr || !iter2->fPtr) );
      
      sizeInconsistency = sizeInconsistency || sizeInc;

      if( !sizeInc ){
	bool RCUtrailer = 0;
 	for( AliHLTInt32_t i=0; i<nWords; i++){
	  if( RCUtrailer || i<nWordsHeader || ( (i-nWordsHeader)%5==0) ){
	    // integer data
	    if( p1[i]>>30 == 0x2 ) RCUtrailer = 1;
	    long int d = (long int) p1[i] - (long int) p2[i];
	    if( abs(d)>intDiff ) intDiff = d;
	  } else {
	    // float data
	    AliHLTFloat64_t f1 = *(AliHLTFloat32_t*)&p1[i];
	    AliHLTFloat64_t f2 = *(AliHLTFloat32_t*)&p2[i];
	    double w = fabs(f1 + f2)/2;
	    if( w>1.e-20 ){
	      AliHLTFloat64_t d = fabs(f1 - f2)/w;
	      if( d > floatDiff ) floatDiff = d;	      
	    }
	  }
	}
      }         
    }

    HLTInfo("HWCF consistency check for slice %d patch %d: wrong NBlocks: %d, wrong size: %d, intDiff: %d, floatDiff: %.10f %s", slice1, patch1, 
	    (nSecondBlocksFound!=1),sizeInconsistency,intDiff, floatDiff*100.,"%"); 
    if( (nSecondBlocksFound!=1) || sizeInconsistency || intDiff!=0 || floatDiff>1.e-6 ){
      fNDismatch++;
      HLTWarning("HWCF consistency check for slice %d patch %d: wrong NBlocks: %d, wrong size: %d, intDiff: %d, floatDiff: %.10f %s", slice1, patch1, 
		 (nSecondBlocksFound!=1),sizeInconsistency,intDiff, floatDiff*100.,"%"); 
    }
    fNBlocks++;
  }
  
  delete[] checkedBlocks;
   
  HLTInfo("HWCF consistency check: %.10f %s of %ld data blocks are OK",
	  (double)(fNBlocks-fNDismatch)/(double)fNBlocks*100.,"%",fNBlocks);
  
  if( fNDismatch>0 ){
    HLTWarning("HWCF inconsistency: %ld of %ld data blocks are not OK",
	       fNDismatch,fNBlocks);      
  }  
  
  fBenchmark.Stop(0);  
  HLTInfo(fBenchmark.GetStatistics());
  return iResult;
}


