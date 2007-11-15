// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCClusterFinderComponent.cxx
    @author Timm Steinbeck, Matthias Richter, Jochen Thaeder, Kenneth Aamodt
    @date   
    @brief  The TPC cluster finder processing component
*/

// see header file for class documentation                                   //
// or                                                                        //
// refer to README to build package                                          //
// or                                                                        //
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt                          //

#if __GNUC__>= 3
using namespace std;
#endif
#include "AliHLTTPCClusterFinderComponent.h"
#include "AliHLTTPCDigitReaderPacked.h"
#include "AliHLTTPCDigitReaderUnpacked.h"
#include "AliHLTTPCDigitReaderRaw.h"
#include "AliHLTTPCClusterFinder.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusters.h"
#include "AliHLTTPCDefinitions.h"
#include <cstdlib>
#include <cerrno>
#include "TString.h"
#include <sys/time.h>

// this is a global object used for automatic component registration, do not use this
// use fPackedSwitch = true for packed inputtype "gkDDLPackedRawDataType"
// use fPackedSwitch = false for unpacked inputtype "gkUnpackedRawDataType"
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponentPacked(true);
AliHLTTPCClusterFinderComponent gAliHLTTPCClusterFinderComponentUnpacked(false);

ClassImp(AliHLTTPCClusterFinderComponent)

AliHLTTPCClusterFinderComponent::AliHLTTPCClusterFinderComponent(bool packed)
  :
  fClusterFinder(NULL),
  fReader(NULL),
  fClusterDeconv(true),
  fXYClusterError(-1),
  fZClusterError(-1),
  fPackedSwitch(packed),
  fUnsorted(0),
  fPatch(0),
  fPadArray(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCClusterFinderComponent::~AliHLTTPCClusterFinderComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCClusterFinderComponent::GetComponentID()
{
  // see header file for class documentation
  if (fPackedSwitch) return "TPCClusterFinderPacked";
  else return "TPCClusterFinderUnpacked";
}

void AliHLTTPCClusterFinderComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.clear(); 
  if (fPackedSwitch) list.push_back( AliHLTTPCDefinitions::fgkDDLPackedRawDataType );
  else list.push_back( AliHLTTPCDefinitions::fgkUnpackedRawDataType );
   
}

AliHLTComponentDataType AliHLTTPCClusterFinderComponent::GetOutputDataType()
{
  // see header file for class documentation
  return AliHLTTPCDefinitions::fgkClustersDataType;
}

void AliHLTTPCClusterFinderComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.  
  constBase = 0;
  if (fPackedSwitch)  inputMultiplier = (6 * 0.4);
  else  inputMultiplier = 0.4;
}

AliHLTComponent* AliHLTTPCClusterFinderComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCClusterFinderComponent(fPackedSwitch);
}
	
int AliHLTTPCClusterFinderComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  if ( fClusterFinder )
    return EINPROGRESS;

  fClusterFinder = new AliHLTTPCClusterFinder();

  Int_t rawreadermode =  -1;
  Int_t sigthresh = -1;
  Float_t occulimit = 1.0;
  Int_t oldRCUFormat=0;
  // Data Format version numbers:
  // 0: RCU Data format as delivered during TPC commissioning, pads/padrows are sorted, RCU trailer is one 32 bit word.
  // 1: As 0, but pads/padrows are delivered "as is", without sorting
  // 2: As 0, but RCU trailer is 3 32 bit words.
  // 3: As 1, but RCU trailer is 3 32 bit words.
  // -1: use offline raw reader

  Int_t i = 0;
  Char_t* cpErr;

  while ( i < argc ) {      

    // -- raw reader mode option
    if ( !strcmp( argv[i], "rawreadermode" ) ) {
      if ( argc <= i+1 ) {
	Logging( kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Missing rawreadermode", "Raw Reader Mode not specified" );
	return ENOTSUP;
      }

      // Decodes the rawreader mode: either number or string and returns the rawreadermode
      // -1 on failure, -2 for offline
      rawreadermode = AliHLTTPCDigitReaderRaw::DecodeMode( argv[i+1] );

      if (rawreadermode == -1 ) {
	Logging( kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Missing rawreadermode", "Cannot convert rawreadermode specifier '%s'.", argv[i+1] );
	return EINVAL;
      }

      i += 2;
      continue;
    }

    // -- pp run option
    if ( !strcmp( argv[i], "pp-run" ) ) {
      fClusterDeconv = false;
      i++;
      continue;
    }

    // -- zero suppression threshold
    if ( !strcmp( argv[i], "adc-threshold" ) ) {
      sigthresh = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ) {
	HLTError("Cannot convert threshold specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- pad occupancy limit
    if ( !strcmp( argv[i], "occupancy-limit" ) ) {
      occulimit = strtod( argv[i+1], &cpErr);
      if ( *cpErr ) {
	HLTError("Cannot convert occupancy specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- number of timebins (default 1024)
    if ( !strcmp( argv[i], "timebins" ) ) {
      TString parameter(argv[i+1]);
      parameter.Remove(TString::kLeading, ' '); // remove all blanks
      if (parameter.IsDigit()) {
	AliHLTTPCTransform::SetNTimeBins(parameter.Atoi());
	HLTInfo("number of timebins set to %d, zbin=%f", AliHLTTPCTransform::GetNTimeBins(), AliHLTTPCTransform::GetZWidth());
      } else {
	HLTError("Cannot timebin specifier '%s'.", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    // -- checking for rcu format
    if ( !strcmp( argv[i], "oldrcuformat" ) ) {
      oldRCUFormat = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert oldrcuformat specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
      
    // -- checking for unsorted clusterfinding
    if ( !strcmp( argv[i], "unsorted" ) ) {
      fUnsorted = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert unsorted specifier '%s'. Should  be 0(off) or 1(on), must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }
      
    // -- checking for unsorted clusterfinding
    if ( !strcmp( argv[i], "patch" ) ) {
      fPatch = strtoul( argv[i+1], &cpErr ,0);
      if ( *cpErr ){
	HLTError("Cannot convert patch specifier '%s'. Should  be between 0 and 5, must be integer", argv[i+1]);
	return EINVAL;
      }
      i+=2;
      continue;
    }

    Logging(kHLTLogError, "HLT::TPCClusterFinder::DoInit", "Unknown Option", "Unknown option '%s'", argv[i] );
    return EINVAL;

  }

  // Choose reader

  if (fPackedSwitch) { 
    if (rawreadermode == -2) {
#if defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
      fReader = new AliHLTTPCDigitReaderPacked();
      if(oldRCUFormat==1){
	fReader->SetOldRCUFormat(kTRUE);
      }
      else if(oldRCUFormat!=0){
	HLTWarning("Wrong oldrcuformat specifier %d; oldrcuformat set to default(kFALSE)",oldRCUFormat);
      }
      if(fUnsorted){
	fReader->SetUnsorted(kTRUE);
      }
      fClusterFinder->SetReader(fReader);
#else // ! defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
      HLTFatal("DigitReaderPacked not available - check your build");
      return -ENODEV;
#endif //  defined(HAVE_ALIRAWDATA) && defined(HAVE_ALITPCRAWSTREAM_H)
    } else {
#if defined(HAVE_TPC_MAPPING)
      fReader = new AliHLTTPCDigitReaderRaw(rawreadermode);
      fClusterFinder->SetReader(fReader);
#else //! defined(HAVE_TPC_MAPPING)
      HLTFatal("DigitReaderRaw not available - check your build");
      return -ENODEV;
#endif //defined(HAVE_TPC_MAPPING)
    }
  }
  else {
    fReader = new AliHLTTPCDigitReaderUnpacked();
    fClusterFinder->SetReader(fReader);
  }

  // if pp-run use occupancy limit else set to 1. ==> use all 
  if ( !fClusterDeconv )
    fClusterFinder->SetOccupancyLimit(occulimit);
  else 
    fClusterFinder->SetOccupancyLimit(1.0);
      
  // Variables to setup the Clusterfinder
  // TODO: this sounds strange and has to be verified; is the cluster finder not working when
  // fClusterDeconv = false ?
  fClusterDeconv = true;
  fXYClusterError = -1;
  fZClusterError = -1;

 
  fClusterFinder->SetDeconv( fClusterDeconv );
  fClusterFinder->SetXYError( fXYClusterError );
  fClusterFinder->SetZError( fZClusterError );
  if ( (fXYClusterError>0) && (fZClusterError>0) )
    fClusterFinder->SetCalcErr( false );
  fClusterFinder->SetSignalThreshold(sigthresh);
    
  if(fUnsorted&&fPatch>-1&&fPatch<6){
    fPadArray = new AliHLTTPCPadArray(fPatch);
    fPadArray->InitializeVector();
  }

  return 0;
}

int AliHLTTPCClusterFinderComponent::DoDeinit()
{
  // see header file for class documentation

  if ( fClusterFinder )
    delete fClusterFinder;
  fClusterFinder = NULL;
 
  if ( fReader )
    delete fReader;
  fReader = NULL;
    
  return 0;
}

int AliHLTTPCClusterFinderComponent::DoEvent( const AliHLTComponentEventData& evtData, 
					      const AliHLTComponentBlockData* blocks, 
					      AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					      AliHLTUInt32_t& size, 
					      vector<AliHLTComponentBlockData>& outputBlocks )
{
  // see header file for class documentation

  //  == init iter (pointer to datablock)
  const AliHLTComponentBlockData* iter = NULL;
  unsigned long ndx;

  //  == OUTdatatype pointer
  AliHLTTPCClusterData* outPtr;

  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, nSize, tSize = 0;

  outBPtr = outputPtr;
  outPtr = (AliHLTTPCClusterData*)outBPtr;

  Int_t slice, patch, row[2];
  unsigned long maxPoints, realPoints = 0;

  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;


      if (fPackedSwitch) {	
	char tmp1[14], tmp2[14];
	DataType2Text( iter->fDataType, tmp1 );
	DataType2Text( AliHLTTPCDefinitions::fgkDDLPackedRawDataType, tmp2 );
	Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Event received", 
		 "Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
		 evtData.fEventID, evtData.fEventID, tmp1, tmp2 );

	if ( iter->fDataType != AliHLTTPCDefinitions::fgkDDLPackedRawDataType ) continue;

      }
      else {
	char tmp1[14], tmp2[14];
	DataType2Text( iter->fDataType, tmp1 );
	DataType2Text( AliHLTTPCDefinitions::fgkUnpackedRawDataType, tmp2 );
	Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Event received", 
		 "Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",
		 evtData.fEventID, evtData.fEventID, tmp1, tmp2 );

	if ( iter->fDataType != AliHLTTPCDefinitions::fgkUnpackedRawDataType ) continue;

      }
    	
      slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
      patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
      row[0] = AliHLTTPCTransform::GetFirstRow( patch );
      row[1] = AliHLTTPCTransform::GetLastRow( patch );
	
      outPtr = (AliHLTTPCClusterData*)outBPtr;

#ifndef KENNETH
      maxPoints = (size-tSize-sizeof(AliHLTTPCClusterData))/sizeof(AliHLTTPCSpacePointData);
#else
      maxPoints = (size-tSize-sizeof(AliHLTTPCClusters))/sizeof(AliHLTTPCSpacePointData);
#endif 

      fClusterFinder->InitSlice( slice, patch, row[0], row[1], maxPoints );
      fClusterFinder->SetOutputArray( (AliHLTTPCSpacePointData*)outPtr->fSpacePoints );
	
      if(fUnsorted){


	fClusterFinder->SetPadArray(fPadArray);
	  
	double totalT=0;
	struct timeval startT, endT;
	gettimeofday( &startT, NULL );

	fClusterFinder->ReadDataUnsorted(iter->fPtr, iter->fSize );

	gettimeofday( &endT, NULL );
	unsigned long long dt;
	dt = endT.tv_sec-startT.tv_sec;
	dt *= 1000000ULL;
	dt += endT.tv_usec-startT.tv_usec;
	double dtd = ((double)dt);
	totalT += dtd;
	//	  dtd = dtd / (double)eventIterations;
	//	  if ( iterations<=1 )
	cout<<endl;
	printf( "Time needed to read data: %f microsec. / %f millisec. / %f s\n", 
		dtd, dtd/1000.0, dtd/1000000.0 );
	  
	cout<<endl;
	fClusterFinder->FindClusters();
      }
      else{
	fClusterFinder->Read(iter->fPtr, iter->fSize );
	fClusterFinder->ProcessDigits();
      }
      realPoints = fClusterFinder->GetNumberOfClusters();
	
      outPtr->fSpacePointCnt = realPoints;
      nSize = sizeof(AliHLTTPCSpacePointData)*realPoints;
#ifndef KENNETH
      mysize += nSize+sizeof(AliHLTTPCClusterData);
#else
      mysize += nSize+sizeof(AliHLTTPCClusters);
#endif 

      Logging( kHLTLogDebug, "HLT::TPCClusterFinder::DoEvent", "Spacepoints", 
	       "Number of spacepoints: %lu Slice/Patch/RowMin/RowMax: %d/%d/%d/%d.",
	       realPoints, slice, patch, row[0], row[1] );
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fSpecification = iter->fSpecification;
      //AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
      outputBlocks.push_back( bd );
	
      tSize += mysize;
      outBPtr += mysize;
      outPtr = (AliHLTTPCClusterData*)outBPtr;
	
      if ( tSize > size )
	{
	  Logging( kHLTLogFatal, "HLT::TPCClusterFinder::DoEvent", "Too much data", 
		   "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",
		   tSize, size );
	  return EMSGSIZE;
	}
    }
    
  size = tSize;

  return 0;
}
