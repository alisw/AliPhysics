// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTCOMPHuffmanAltroComponent.cxx
    @author Jenny Wagner
    @date   29-08-2007
    @brief  The Huffman compressor component.
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTCOMPHuffmanAltroComponent.h"
#include "AliHLTCOMPHuffmanAltro.h"
#include "AliHLTCOMPHuffmanData.h"
#include "AliHLTCompDefinitions.h"
#include "AliHLTStdIncludes.h"
#include "TFile.h"

ClassImp(AliHLTCOMPHuffmanAltroComponent)

/* constructur with arguments */
AliHLTCOMPHuffmanAltroComponent::AliHLTCOMPHuffmanAltroComponent(bool compression)
  :
  fHuffmanCompressor(NULL),
  fCompressionSwitch(compression),
  fTrainingMode(kFALSE),
  fOrigin(kAliHLTVoidDataOrigin),
  fRunNumber(0),
  fDataSpec(0),
  fTablePath(),
  fNrcuTrailerwords(0),
  fHuffmanData(NULL) 
{
   // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCOMPHuffmanAltroComponent::~AliHLTCOMPHuffmanAltroComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTCOMPHuffmanAltroComponent::GetComponentID()
{
  // see header file for class documentation
  if(fCompressionSwitch)
    {
      return "COMPHuffmanCompressor";
    }
  else
    {
      return "COMPHuffmanDecompressor";
    }
}

void AliHLTCOMPHuffmanAltroComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  // initialise list
  list.clear(); 

  // if compression is to be done, input data is packed raw data
  // else (decompression): input data is special entropy encoded raw data
  if (fCompressionSwitch) list.push_back( kAliHLTDataTypeDDLRaw );
  else list.push_back( AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType);
  
}

AliHLTComponentDataType AliHLTCOMPHuffmanAltroComponent::GetOutputDataType()
{
  // see header file for class documentation
  // if compression is to be one, output data is special entropy encoded raw data
  // else (decompression): output data is packed raw data
  AliHLTComponentDataType dt=kAliHLTDataTypeDDLRaw;
  if(fCompressionSwitch)
    dt=AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType;
  if (!fOrigin.IsNull()) SetDataType(dt, NULL, fOrigin.Data());
  return dt;
}

void AliHLTCOMPHuffmanAltroComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // reserved outputsize = inputside * inputMultiplier
  constBase = 0;
  if (fCompressionSwitch == kFALSE)
    {
      // for decompression: compressed * 4 = (enough space for) decompressed
      inputMultiplier = 4.0 ;
    }
  else
    {
      // for compression: original * 1 = (enough space for) compressed
      inputMultiplier = 1.0;
    }
}

AliHLTComponent* AliHLTCOMPHuffmanAltroComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTCOMPHuffmanAltroComponent(fCompressionSwitch);
}


int AliHLTCOMPHuffmanAltroComponent::DoInit( int argc, const char** argv )
{
  
  // see header file for class documentation
  
  if ( fHuffmanCompressor )
    return EINPROGRESS;
  
  
  Int_t i = 0;
  //Char_t* cpErr;
  
  while ( i < argc ) 
    {     
      // -- training mode wrongly called here
      if ( !strcmp( argv[i], "-training" ) ) 
	{
	  fTrainingMode = kTRUE;
	  
	  HLTInfo("HuffmanCompressor called in training mode, please call HuffmanCalibration instead.");

	  return EINVAL;
	}
      
      // mode option: compress or decompress
      //  if ( !strcmp( argv[i], "-compress" ) ) 
      //	{
      //	  fCompressionSwitch = kTRUE;

      //  ++i;
      //  continue;
	  //}
    
      //if(!strcmp( argv[i], "-decompress" ) ) 
      //	{
      //	  fCompressionSwitch = kFALSE;
	  
      //	  ++i;
      //  continue;
      //	}
   
      // -- argument to load correct code table for respective data specifications (origin, runnumber, specification)
      if ( !strcmp( argv[i], "-origin" ) ) 
	{
	 
	  if ( argc <= i+1 ) 
	    {
	      HLTError("Missing data origin specification");
	      return ENOTSUP;
	    }

	  // get data origin (TPC, PHOS, ...)
	  fOrigin=argv[i+1];

	  while(fOrigin.Length() <  kAliHLTComponentDataTypefOriginSize)
	    {
	      fOrigin.Append(" ");
	    }
     
	  HLTDebug("Origin is set to %s.", fOrigin.Data());
         
	  i += 2;
	  continue;
	}

      // -- get run number specification
      if ( !strcmp( argv[i], "-runnumber" ) ) 
	{
	  if ( argc <= i+1 ) 
	    {
	      HLTError("Missing run number specification");
	      return ENOTSUP;
	    }

	  // get run number
	  const char* runnumber = argv[i+1];
	  
	  fRunNumber =  atoi(runnumber);

	  HLTDebug("Run number set to %d (Dec) = %X (Hex).", fRunNumber, fRunNumber);
	  
	  // validation check of run number?!

	  i += 2;
	  continue;
	}


      // -- get data specification (e.g. TPC: slice and patch information contained in "dataspec")
      if ( !strcmp( argv[i], "-dataspec" ) ) 
	{
	  if ( argc <= i+1 ) 
	    {
	      HLTError("Missing data specification");
	      return ENOTSUP;
	    }

	  // get data spec
	  fDataSpec = strtoul(argv[i+1], NULL, 16);

	  HLTDebug("Dataspecification set to %d (Dec) = %08X (Hex).", fDataSpec, fDataSpec);
	  
	  // validation check of specification?!

	  i += 2;
	  continue;
	}

      // -- get tablepathname (e.g. ../HLT-data/)
      if ( !strcmp( argv[i], "-tablepath" ) ) 
	{
	  if ( argc <= i+1 ) 
	    {
	      HLTDebug("Missing table path argument.");
	    }

	  // get data spec
	  fTablePath=argv[i+1];

	  HLTDebug("Table path set to %s.", fTablePath.Data());
	  
	  // validation check of specification?!

	  i += 2;
	  continue;
	}

      // -- number of trailerwords: from 1 to 3
      if ( !strcmp( argv[i], "-trailerwords" ) ) 
	{
	  if ( argc <= i+1 ) 
	    {
	      HLTError("Missing trailerword specification");
	      return ENOTSUP;
	    }
	  
	  if ( !strcmp( argv[i+1], "1" ) ) 
	      fNrcuTrailerwords = 1;
	  else if ( !strcmp( argv[i+1], "2" ) ) 
	    fNrcuTrailerwords = 2; 
	  else if ( !strcmp( argv[i+1], "3" ) ) 
	    fNrcuTrailerwords = 3; 
	  else
	    {
	      HLTError("Missing number of trailerwords, cannot accept argument '%s'.", argv[i+1] );
	      
	      return EINVAL;
	    }	
	  
	  i += 2;
	  continue;
	}
            
      HLTError("Unknown Option '%s'", argv[i] );
      return EINVAL;
    } // end while-loop
  
  // load HuffmanData from root-file to acquire translation table
  fHuffmanData = new AliHLTCOMPHuffmanData();

  TString rootfilename;
  if(fTablePath.IsNull())
    {
      // if no table path is explicity set, take current path as table path
      rootfilename.Form("huffmanData_%s_%08X_%08X.root", fOrigin.Data(), (unsigned)fRunNumber, (unsigned)fDataSpec);   
    }
  else
    {
      rootfilename.Form("%shuffmanData_%s_%08X_%08X.root", fTablePath.Data(), fOrigin.Data(), (unsigned)fRunNumber, (unsigned)fDataSpec);
    }

  TFile* huffmancodefile = new TFile(rootfilename, "READ");
  
  if(huffmancodefile->IsZombie())
  { 
    HLTFatal("No Huffman code table available for %s.", rootfilename.Data());
      return EINVAL;
    }    
  fHuffmanData = (AliHLTCOMPHuffmanData*) huffmancodefile->Get("HuffmanData");

  // create a new Huffman compressor 
  fHuffmanCompressor = new AliHLTCOMPHuffmanAltro(fCompressionSwitch, kFALSE, NULL, fNrcuTrailerwords);
  
  // get translation table for pure encoding and decoding from HuffmanData
  fHuffmanCompressor->GetTranslationTable(fHuffmanData);
    
  return 0;
}

int AliHLTCOMPHuffmanAltroComponent::DoDeinit()
{

  // see header file for class documentation
  if (fHuffmanCompressor) 
    delete fHuffmanCompressor;
  fHuffmanCompressor = NULL;

  if ( fHuffmanData )
    delete fHuffmanData; 
  fHuffmanData = NULL;
  
  return 0;
}

int AliHLTCOMPHuffmanAltroComponent::DoEvent( const AliHLTComponentEventData& evtData, 
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
  //    AliHLTTPCClusterData* outPtr;
  
  AliHLTUInt8_t* outBPtr;
  UInt_t offset, mysize, tSize = 0;
  
  outBPtr = outputPtr;
  //    outPtr = (AliHLTTPCClusterData*)outBPtr;
  
  
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      mysize = 0;
      offset = tSize;
      
      if(fCompressionSwitch) // show selected mode
	{
	  HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",evtData.fEventID, evtData.fEventID, 
		   DataType2Text(iter->fDataType).c_str(), DataType2Text(kAliHLTDataTypeDDLRaw|fOrigin.Data()).c_str());

	  // check if current block has correct data format
	  // if not, take next block
	  if ( iter->fDataType != (kAliHLTDataTypeDDLRaw|fOrigin.Data()) ) continue;
	}
      else
	{
	  HLTDebug("Event 0x%08LX (%Lu) received datatype: %s - required datatype: %s",evtData.fEventID, evtData.fEventID,
		   DataType2Text(iter->fDataType).c_str(), DataType2Text(AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType|fOrigin.Data()).c_str());

	  // check if current block has correct data format
	  // if not, take next block
	  if ( iter->fDataType != (AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType|fOrigin.Data()) ) continue;
	}

      // HLTDebug("HLT::HuffmanCompressor::DoEvent", "Event received", "Starting to process data");

      // FIXME: set ddl no
      fHuffmanCompressor->AddInputData(reinterpret_cast<UChar_t*>(iter->fPtr), iter->fSize, 768);

      // validation test
      // HLTDebug("input data pointer (HEX) = %x ", iter->fPtr);
      // HLTDebug("input data size (bytes) = %i ", iter->fSize);
      
      fHuffmanCompressor->SetOutputData(outBPtr, size);
      
      // validation test
      // HLTDebug("output data pointer (HEX) = %x ", outBPtr);
      // HLTDebug("reserved output data size (bytes) = %i ", size);
      
      fHuffmanCompressor->ProcessData();
      
      //	outPtr = (AliHLTTPCClusterData*)outBPtr;
      
      mysize = fHuffmanCompressor->GetOutputDataSize();
      
      if(mysize != 0)
	{
	  AliHLTComponentBlockData bd;
	  FillBlockData( bd );
	  bd.fOffset = offset;
	  bd.fSize = mysize;
	  bd.fSpecification = iter->fSpecification;
	  //AliHLTSubEventDescriptor::FillBlockAttributes( bd.fAttributes );
	  outputBlocks.push_back( bd );
	    
	  tSize += mysize;
	    outBPtr += mysize;
	    //outPtr = (AliHLTTPCClusterData*)outBPtr;
	    
	    if ( tSize > size )
	      {
		HLTFatal("HLT::TPCHuffmanCompressor::DoEvent: Too much data, data written over allowed buffer. Amount written: %lu, allowed amount: %lu.",tSize, size );
		return EMSGSIZE;
	      }
	    
	} // end of output-block-generation
      
    }
  
  size = tSize;

  return 0;   
}
   
