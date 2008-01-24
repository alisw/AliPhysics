// $Id: AliHLTTPCCompModelInflater.cxx,v 1.2 2006/08/10 09:46:51 richterm Exp $

/**************************************************************************
 * TPCCompModelInflaterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCCompModelInflater.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelInflater.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include "AliHLTDataTypes.h"
#include <cerrno>

AliHLTTPCCompModelInflater::AliHLTTPCCompModelInflater():
    fBitDataCurrentWord(0),
    fBitDataCurrentPosInWord(0),
    fBitDataCurrentInput(0),
    fBitDataCurrentInputStart(0),
    fBitDataCurrentInputEnd(0)
    {
      // see header file for class documentation
    }

AliHLTTPCCompModelInflater::~AliHLTTPCCompModelInflater()
    {
      // see header file for class documentation
    }

int AliHLTTPCCompModelInflater::DecompressTracks( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize )
    {
      // see header file for class documentation
      AliHLTUInt8_t* inputPtr = inData;
      //AliHLTUInt8_t* inputEndPtr = inData+inputSize;
      AliHLTUInt8_t* outputPtr = output;
      AliHLTUInt8_t* outputEndPtr = output+outputSize;
      
      printf( "outuptSize: %lu\n", (unsigned long)outputSize );
      
      InitBitDataInput( inputPtr, inputSize );
      AliHLTUInt8_t version;
      if ( !InputBits( version, 4 ) ) // Version information
	{
	  HLTError( "Corrupt input data. Cannot read data version number at position %u",
		    (unsigned)(inputPtr-inData) );
	  return EIO;
	}
      if ( version != 0 )
	{
	  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
	  return EIO;
	}
      AliHLTUInt8_t readShape;
      if ( !InputBit( readShape ) ) // Data format flag
	{
	  HLTError( "Corrupt input data. Cannot read shape flag at position %u",
		    (unsigned)(inputPtr-inData) );
	  return EIO;
	}
      
      if ( outputPtr+sizeof(AliHLTUInt32_t)>outputEndPtr )
	{
	  HLTError( "Not enough space to write decompressed data. %lu already written",
		    (unsigned)(outputPtr-output) );
	  return ENOBUFS;
	}
      *(AliHLTUInt32_t*)outputPtr = 0; // Write format version number
      outputPtr += sizeof(AliHLTUInt32_t);
      
      Pad8Bits();
      
      AliHLTTPCClusterModel *cluster;
      
      Int_t timeo,pado,chargeo,padshapeo,timeshapeo;
      timeo=pado=chargeo=padshapeo=timeshapeo=0;
      unsigned trackCnt=0;
      while( !EndOfBitInput() )
	{
	  if ( outputPtr+sizeof(AliHLTTPCTrackModel)>outputEndPtr )
	    {
	      HLTError( "Not enough space to write decompressed data. %lu already written",
			(unsigned)(outputPtr-output) );
	      return ENOBUFS;
	    }
	  if ( !InputBytes( outputPtr, sizeof(AliHLTTPCTrackModel) ) )
	    {
	      HLTError( "Corrupt input data. Cannot read track model data at position %u",
			(unsigned)(inputPtr-inData) );
	      return EIO;
	    }
	  outputPtr += sizeof(AliHLTTPCTrackModel);
	  
	  Int_t clustercount=0;
	  AliHLTUInt32_t slice;
	  for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(); i++)
	    {
	      if ( outputPtr+sizeof(AliHLTTPCClusterModel)>outputEndPtr )
		{
		  HLTError( "Not enough space to write decompressed data. %lu already written",
			    (unsigned)(outputPtr-output) );
		  return ENOBUFS;
		}
	      cluster = (AliHLTTPCClusterModel*)outputPtr;
	      outputPtr += sizeof(AliHLTTPCClusterModel);
	      
	      // Read present flag:
	      AliHLTUInt8_t present;
	      if ( !InputBit( present ) )
		{
		  HLTError( "Corrupt input data. Cannot read cluster presence bit at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      HLTDebug( "Cluster for row %d %s", i, (present ? "present" : "not present") );
	      cluster->fPresent = present;
	      if ( !present )
		continue;
	      
	      
	      //Read slice number of first point
	      if ( clustercount==0 )
		{
		  if ( !InputBits( slice,6 ) ) //Need 6 bits to encode slice number
		    {
		      HLTError( "Corrupt input data. Cannot read cluster slice number at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		}
	      else
		{
		  AliHLTUInt8_t sliceChange;
		  if ( !InputBit( sliceChange ) )
		    {
		      HLTError( "Corrupt input data. Cannot read cluster slice change bit at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  if ( sliceChange )
		    {  //Change of slice
		      if ( !InputBits( slice, 6 ) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster slice number at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		    }
		}
	      HLTDebug( "Slice: %d", slice );
	      cluster->fSlice = slice;
	      if ( cluster->fSlice<0 || cluster->fSlice>35 )
		{
		  HLTError( "Inconsistent slice number %u (track %u, cluster %d)", cluster->fSlice, trackCnt, i );
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  return EINVAL;
		}
	      
	      AliHLTUInt8_t signBit;
	      Int_t sign;
	      AliHLTUInt64_t temp;
	      Int_t val;
	      //Read time information:
	      if ( !InputBit( signBit ) )
		{
		  HLTError( "Corrupt input data. Cannot read DTime sign bit at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      sign = signBit;
	      sign = -1+sign*2;
	      if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNTimeBits()-1 ) )
		{
		  HLTError( "Corrupt input data. Cannot read DTime data at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      val = (Int_t)temp;
	      cluster->fDTime = val*sign;
	      
	      //Read pad information:
	      if ( !InputBit( signBit ) )
		{
		  HLTError( "Corrupt input data. Cannot read DPad sign bit at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      sign = signBit;
	      sign = -1+sign*2;
	      if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNPadBits()-1 ) )
		{
		  HLTError( "Corrupt input data. Cannot read DPad data at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      
	      val = (Int_t)temp;
	      cluster->fDPad = val*sign;
	      
	      // Read charge information:
	      if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNChargeBits() ) )
		{
		  HLTError( "Corrupt input data. Cannot read charge data at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      cluster->fDCharge = temp;
	      
	      if ( readShape )
		{
		  // Read shape information:
		  if ( !InputBit( signBit ) )
		    {
		      HLTError( "Corrupt input data. Cannot read DSigmaY sign bit at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  sign = signBit;
		  sign = -1+sign*2;
		  if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1 ) )
		    {
		      HLTError( "Corrupt input data. Cannot read DSigmaY data at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  
		  val = (Int_t)temp;
		  cluster->fDSigmaY = val*sign;
		  //HLTInfo("DSigmaY: %f", cluster->fDSigmaY);
		  
		  
		  if ( !InputBit( signBit ) )
		    {
		      HLTError( "Corrupt input data. Cannot read DSigmaZ sign bit at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  sign = signBit;
		  sign = -1+sign*2;
		  if ( !InputBits( temp, AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1 ) )
		    {
		      HLTError( "Corrupt input data. Cannot read DSigmaZ data at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  val = (Int_t)temp;
		  cluster->fDSigmaZ = val*sign;
		}
	      
	      clustercount++;
	    }
	  Pad8Bits();
	  HLTDebug( "Track %u: %d clusters", trackCnt, clustercount );
	}
      
      outputSize = (UInt_t)( outputPtr - output );
      return 0;
    }

int AliHLTTPCCompModelInflater::DecompressRemainingClusters( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize )
    {
      // see header file for class documentation
      AliHLTUInt8_t* inputPtr = inData;
      AliHLTUInt8_t* outputPtr = output;
      AliHLTUInt8_t* outputEndPtr = output+outputSize;
      
      InitBitDataInput( inputPtr, inputSize );
      AliHLTUInt8_t version;
      if ( !InputBits( version, 4 ) ) // Version information
	{
	  HLTError( "Corrupt input data. Cannot read data version number at position %u",
		    (unsigned)(inputPtr-inData) );
	  return EIO;
	}
      if ( version != 0 )
	{
	  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
	}
      Pad8Bits();
      
      if ( outputPtr+sizeof(AliHLTUInt32_t)>outputEndPtr )
	{
	  HLTError( "Not enough space to write uncompressed data. %lu already written",
		    (unsigned long)(outputPtr-output) );
	  outputSize = (unsigned long)(outputPtr-output);
	  return ENOBUFS;
	}
      
      *(AliHLTUInt32_t*)outputPtr = 0; // Write format version
      outputPtr += sizeof(AliHLTUInt32_t);
      
      //Read the remaining clusters in a compressed format.
      
      for(Int_t slice=0; slice<=35; slice++)
	{
	  for(Int_t patch=0; patch < 6; patch++)
	    {
	      UInt_t i;
	      HLTDebug( "slice %u patch %u: %u padrows",
			(unsigned)slice, (unsigned)patch, (unsigned)*inputPtr );
	      //Write number of padrows with clusters
	      if ( outputPtr>=outputEndPtr )
		{
		  HLTError( "Not enough space to write uncompressed data. %lu already written",
			    (unsigned long)(outputPtr-output) );
		  outputSize = (unsigned long)(outputPtr-output);
		  return ENOBUFS;
		}
	      if ( !InputBits( *outputPtr,8 ) )
		{
		  HLTError( "Corrupt input data. Cannot read padrow count at position %u",
			    (unsigned)(inputPtr-inData) );
		  return EIO;
		}
	      if ( !*outputPtr )
		{
		  outputPtr++;
		  continue;
		}
	      UInt_t nRows=(UInt_t)*outputPtr;
	      outputPtr++;
	      if ( outputPtr>=outputEndPtr )
		{
		  HLTError( "Not enough space to write uncompressed data. %lu already written",
			    (unsigned long)(outputPtr-output) );
		  outputSize = (unsigned long)(outputPtr-output);
		  return ENOBUFS;
		}
	      
	      for ( UInt_t jj=0; jj<nRows; jj++ )
		{
		  
		  AliHLTTPCRemainingRow *thisRow = (AliHLTTPCRemainingRow*)outputPtr;
		  if ( outputPtr+sizeof(AliHLTTPCRemainingRow)>outputEndPtr )
		    {
		      HLTError( "Not enough space to write uncompressed data. %lu already written",
				(unsigned long)(outputPtr-output) );
		      outputSize = (unsigned long)(outputPtr-output);
		      return ENOBUFS;
		    }
		  AliHLTTPCRemainingCluster *cl = thisRow->fClusters;
		  if ( !InputBits(thisRow->fPadRow,8) ) //Read padrow #
		    {
		      HLTError( "Corrupt input data. Cannot read padrow number at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  if ( !InputBits(thisRow->fNClusters,10) )//Read number of clusters on this padrow
		    {
		      HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				(unsigned)(inputPtr-inData) );
		      return EIO;
		    }
		  if ( outputPtr+sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster)>outputEndPtr )
		    {
		      HLTError( "Not enough space to write uncompressed data. %lu already written",
				(unsigned long)(outputPtr-output) );
		      outputSize = (unsigned long)(outputPtr-output);
		      return ENOBUFS;
		    }
		  for ( i=0; i<thisRow->fNClusters; i++ )
		    {
		      //Read pad
		      Int_t buff;
		      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNPadBitsRemaining()) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		      cl[i].fPad = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor() );
		      
		      //Read time
		      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNTimeBitsRemaining()) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		      cl[i].fTime = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor() );
		      
		      //Read widths
		      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		      
		      //HLTInfo("fDSgimaY = %d",buff);
		      
		      Float_t padw = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor() );
		      cl[i].fSigmaY2 = padw*padw;
		      //HLTInfo("sigmaY2: %f", cl[i].fSigmaY2);
		      
		      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		      Float_t timew = (Float_t)( ((double)buff) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor() );
		      cl[i].fSigmaZ2 = timew*timew;
		      
		      //Read charge 
		      if ( !InputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNChargeBits()) )
			{
			  HLTError( "Corrupt input data. Cannot read cluster count at position %u",
				    (unsigned)(inputPtr-inData) );
			  return EIO;
			}
		      cl[i].fCharge = buff;
		    }
		  outputPtr += sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster);
		}
	    }
	  
	}
      outputSize = (UInt_t)( outputPtr - output );
      return 0;
    }
