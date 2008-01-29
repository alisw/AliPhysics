// $Id$

/**************************************************************************
 * TPCCompModelDeflaterright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

/** @file   AliHLTTPCCompModelDeflater.cxx
    @author Timm Steinbeck
    @date   
    @brief  A copy processing component for the HLT. */

#if __GNUC__ >= 3
using namespace std;
#endif

#include "AliHLTTPCCompModelDeflater.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCTrack.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCCompDataCompressorHelper.h"
#include "AliHLTDataTypes.h"
#include <cerrno>

AliHLTTPCCompModelDeflater::AliHLTTPCCompModelDeflater():
    fWriteShape(true),
    fBitDataCurrentWord(0),
    fBitDataCurrentPosInWord(0),
    fBitDataCurrentOutput(0),
    fBitDataCurrentOutputStart(0),
    fBitDataCurrentOutputEnd(0)
    {
      // see header file for class documentation
    }

AliHLTTPCCompModelDeflater::~AliHLTTPCCompModelDeflater()
    {
      // see header file for class documentation
    }

int AliHLTTPCCompModelDeflater::CompressTracks( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize )
    {
      // see header file for class documentation
      AliHLTUInt8_t* inputPtr = inData;
      AliHLTUInt8_t* inputEndPtr = inData+inputSize;
      
      if ( inputPtr+sizeof(AliHLTUInt32_t)>inputEndPtr )
	{
	  HLTError( "Cannot read input data version number" );
	  return EIO;
	}
      if ( *(AliHLTUInt32_t*)inputPtr != 0 )
	{
	  HLTError( "Only input data format version 0 is supported. Found version: %u",
		    (unsigned)*(AliHLTUInt32_t*)inputPtr );
	  return EINVAL;
	}
      inputPtr += sizeof(AliHLTUInt32_t);
      
      printf( "outuptSize: %lu\n", (unsigned long)outputSize );
      
      InitBitDataOutput( output, outputSize );
      HLTDebug( "Output: Position: %lu / %u (0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte() );
      OutputBits( 0, 4 ); // Version information
      HLTDebug( "Output: Position: %lu / %u (0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte() );
      OutputBit( fWriteShape ? 1 : 0 ); // Data format flag
      HLTDebug( "Output: Position: %lu / %u (0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte() );
      Pad8Bits();
      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
      
      AliHLTTPCClusterModel *cluster;
      Int_t temp;
      Int_t power;
      
      Int_t timeo,pado,chargeo,padshapeo,timeshapeo;
      timeo=pado=chargeo=padshapeo=timeshapeo=0;
      unsigned trackCnt=0;
      while( inputPtr<inputEndPtr )
	{
	  if ( !OutputBytes( inputPtr, sizeof(AliHLTTPCTrackModel) ) )
	    {
	      HLTError( "Not enough space to write compressed data. %lu already written",
			(unsigned long)GetBitDataOutputSizeBytes() );
	      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
	      outputSize = GetBitDataOutputSizeBytes();
	      return ENOBUFS;
	    }
	  HLTDebug( "sizeof(AliHLTTPCTrackModel): %d", sizeof(AliHLTTPCTrackModel) );
	  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
	  inputPtr += sizeof(AliHLTTPCTrackModel);
	  
	  Int_t origslice=-1,slice,clustercount=0;
	  for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(); i++)
	    {
	      cluster = (AliHLTTPCClusterModel*)inputPtr;
	      inputPtr += sizeof(AliHLTTPCClusterModel);
	      
	      //Write empty flag:
	      if ( !OutputBit( cluster->fPresent ? 1 : 0 ) )
		{
		  HLTError( "Not enough space to write compressed data. %lu already written",
			    (unsigned long)GetBitDataOutputSizeBytes() );
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  outputSize = GetBitDataOutputSizeBytes();
		  return ENOBUFS;
		}
	      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
	      if ( !cluster->fPresent )
		continue;
	    
	      if ( cluster->fSlice<0 || cluster->fSlice>35 )
		{
		  HLTError( "Inconsistent slice number %u (track %u, cluster %d)", cluster->fSlice, trackCnt, i );
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  return EINVAL;
		}
	      
	      //Write slice number of first point
	      if ( clustercount==0 )
		{
		  origslice = cluster->fSlice;
		  if ( !OutputBits( origslice,6 ) ) //Need 6 bits to encode slice number
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      else
		{
		  slice = cluster->fSlice;
		  if( slice == origslice )
		    {
		      if ( !OutputBit( 0 ) )  //No change of slice
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  return ENOBUFS;
			}
		      HLTDebug( "No slice change (%d/%d)", (int)origslice, (int)slice );
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		    }
		  else
		    {
		      if ( !OutputBit( 1 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  return ENOBUFS;
			}
		      HLTDebug( "Slice change (%d/%d)", (int)origslice, (int)slice );
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		      if ( !OutputBits( slice, 6 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  return ENOBUFS;
			}
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		      origslice=slice;
		    }
		}
	      
	      //Write time information:
	      temp = (Int_t)rint(cluster->fDTime);
	      if( temp<0 )
		{
		if ( !OutputBit( 0 ) )
		  {
		    HLTError( "Not enough space to write compressed data. %lu already written",
			      (unsigned long)GetBitDataOutputSizeBytes() );
		    outputSize = GetBitDataOutputSizeBytes();
		    printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		    return ENOBUFS;
		  }
		HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      else
		{
		  if ( !OutputBit( 1 ) )
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      power = 1<<(AliHLTTPCCompDataCompressorHelper::GetNTimeBits()-1);
	      if ( abs(temp)>=power )
		{
		  //cout<<abs(temp)<<" "<<power<<endl;
		  timeo++;
		  temp=power - 1;
		}
	      temp = abs(temp);
	      if ( !OutputBits(temp,(AliHLTTPCCompDataCompressorHelper::GetNTimeBits()-1)) )
		{
		  HLTError( "Not enough space to write compressed data. %lu already written",
			    (unsigned long)GetBitDataOutputSizeBytes() );
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  outputSize = GetBitDataOutputSizeBytes();
		  return ENOBUFS;
		}
	      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-2), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
	      
	      //Write pad information:
	      temp = (Int_t)rint(cluster->fDPad);
	      HLTDebug( "cluster->fDPad (%d): %f - temp: %d", clustercount, cluster->fDPad, temp );
	      if ( temp<0 )
		{
		  if ( !OutputBit( 0 ) )
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      outputSize = GetBitDataOutputSizeBytes();
		      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		      return ENOBUFS;
		    }
		  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      else
		{
		if ( !OutputBit( 1 ) )
		  {
		    HLTError( "Not enough space to write compressed data. %lu already written",
			      (unsigned long)GetBitDataOutputSizeBytes() );
		    outputSize = GetBitDataOutputSizeBytes();
		    printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		    return ENOBUFS;
		  }
		HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      power = 1<<(AliHLTTPCCompDataCompressorHelper::GetNPadBits()-1);
	      if ( abs(temp)>=power )
		{
		  pado++;
		  temp=power - 1;
		}
	      temp = abs(temp);
	      if ( !OutputBits(temp,(AliHLTTPCCompDataCompressorHelper::GetNPadBits()-1)) )
		{
		  HLTError( "Not enough space to write compressed data. %lu already written",
			    (unsigned long)GetBitDataOutputSizeBytes() );
		  outputSize = GetBitDataOutputSizeBytes();
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  return ENOBUFS;
		}
	      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-2), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
	      
	      //Write charge information:
	      temp = (Int_t)cluster->fDCharge;
	      power = 1<<(AliHLTTPCCompDataCompressorHelper::GetNChargeBits());
	      if ( abs(temp)>=power )
		{
		  chargeo++;
		  temp=power - 1;
		}
	      temp = abs(temp);
	      if ( !OutputBits(temp,(AliHLTTPCCompDataCompressorHelper::GetNChargeBits())) )
		{
		  HLTError( "Not enough space to write compressed data. %lu already written",
			    (unsigned long)GetBitDataOutputSizeBytes() );
		  outputSize = GetBitDataOutputSizeBytes();
		  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		  return ENOBUFS;
		}
	      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-2), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
	      
	      if ( fWriteShape )
		{
		  //Write shape information:
		  // HLTInfo("DSigmaY %f", cluster->fDSigmaY);
		  temp = (Int_t)rint(cluster->fDSigmaY);
		  // HLTInfo("temp %d", temp);
		  if( temp<0 )
		    {
		      if ( !OutputBit( 0 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		    }
		  else
		    {
		      if ( !OutputBit( 1 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  return ENOBUFS;
			}
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		    }
		  power = 1<<(AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1);
		  if ( abs(temp) >= power )
		    {
		      padshapeo++;
		      temp = power - 1;
		    }
		  temp = abs(temp);
		  if ( !OutputBits(temp,(AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1)) )
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-2), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		  
		  temp = (Int_t)rint(cluster->fDSigmaZ);
		  if ( temp<0 )
		    {
		      if ( !OutputBit( 0 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		    }
		  else
		    {
		      if ( !OutputBit( 1 ) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		    }
		  power = 1<<(AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1);
		  if ( abs(temp) >= power )
		    {
		      timeshapeo++;
		      temp=power - 1;
		    }
		  temp = abs(temp);
		  if ( !OutputBits(temp,(AliHLTTPCCompDataCompressorHelper::GetNShapeBits()-1)) )
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      printf( "TRACE: %s:%d\n", __FILE__, __LINE__ );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-2), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
		}
	      
	      clustercount++;
	    }
	  trackCnt++;
	}
      
      CloseBitDataOutput();
      HLTDebug( "Output: Position: %lu / %u (0x%02X / 0x%02X)", GetCurrentByteOutputPosition(), GetCurrentBitOutputPosition(), (unsigned)GetCurrentOutputByte(-1), (unsigned)GetCurrentOutputByte() );
      outputSize = GetBitDataOutputSizeBytes();
      return 0;
    }

int AliHLTTPCCompModelDeflater::CompressRemainingClusters( AliHLTUInt8_t* inData, UInt_t const& inputSize, AliHLTUInt8_t* output, UInt_t& outputSize )
    {
      // see header file for class documentation
      AliHLTUInt8_t* inputPtr = inData;
      AliHLTUInt8_t* inputEndPtr = inData+inputSize;
      
      AliHLTUInt32_t version = *(AliHLTUInt32_t*)inputPtr;
      inputPtr += sizeof(AliHLTUInt32_t);
      if ( version != 0 )
	{
	  HLTError( "Unsupported version %hu. Only version 0 supported currently.", version );
	  return EIO;
	}
      
      InitBitDataOutput( output, outputSize );
      OutputBits( 0, 4 ); // Version information
      //OutputBit( fWriteShape ); // Data format flag
      Pad8Bits();
      
      //Write the remaining clusters in a compressed format.
      
      for(Int_t slice=0; slice<=35; slice++)
	{
	  for(Int_t patch=0; patch < 6; patch++)
	    {
	      UInt_t i;
	      HLTDebug( "slice %u patch %u: %u padrows",
			(unsigned)slice, (unsigned)patch, (unsigned)*inputPtr );
	      //Write number of padrows with clusters
	      if ( inputPtr>=inputEndPtr )
		{
		  HLTError( "Corrupt input data, cannot read row counter for slice %u, partition %u", (unsigned)slice, (unsigned)patch );
		  return EIO;
		}
	      if ( !OutputBits( *inputPtr,8 ) )
		{
		  HLTError( "Not enough space to write compressed data. %lu already written",
			    (unsigned long)GetBitDataOutputSizeBytes() );
		  outputSize = GetBitDataOutputSizeBytes();
		  return ENOBUFS;
		}
	      if ( !*inputPtr )
		{
		  inputPtr++;
		  continue;
		}
	      UInt_t nRows=(UInt_t)*inputPtr;
	      inputPtr++;
	      if ( inputPtr>=inputEndPtr )
		{
		  HLTError( "Corrupt input data, unexpected end of data after row counter for slice %u, partition %u", (unsigned)slice, (unsigned)patch );
		  return EIO;
		}
	      
	      for ( UInt_t jj=0; jj<nRows; jj++ )
		{
		  
		  AliHLTTPCRemainingRow *thisRow = (AliHLTTPCRemainingRow*)inputPtr;
		  if ( inputPtr+sizeof(AliHLTTPCRemainingRow)>inputEndPtr )
		    {
		      HLTError( "Corrupt input data, cannot read row data for row %u of slice %u, partition %u", (unsigned)jj, (unsigned)slice, (unsigned)patch );
		      return EIO;
		    }
		  AliHLTTPCRemainingCluster *cl = thisRow->fClusters;
		  HLTDebug( "Row %u: %u clusters", (unsigned)thisRow->fPadRow, (unsigned)thisRow->fNClusters );
		  if ( inputPtr+sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster)>inputEndPtr )
		    {
		      HLTError( "Corrupt input data, unable to read clusters for row %u, slice %u, partition %u", (unsigned)jj, (unsigned)slice, (unsigned)patch );
		      return EIO;
		    }
		  Int_t padrow = thisRow->fPadRow;
		  if ( !OutputBits(padrow,8) ) //Write padrow #
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  if( thisRow->fNClusters >= 1<<10)
		    {
		      HLTError( "Too many remaining clusters (%u)", (unsigned)thisRow->fNClusters );
		      return ERANGE;
		    }
		  if ( !OutputBits(thisRow->fNClusters,10) )//Write number of clusters on this padrow
		    {
		      HLTError( "Not enough space to write compressed data. %lu already written",
				(unsigned long)GetBitDataOutputSizeBytes() );
		      outputSize = GetBitDataOutputSizeBytes();
		      return ENOBUFS;
		    }
		  for ( i=0; i<thisRow->fNClusters; i++ )
		    {
		      
		      Float_t padw = sqrt(cl[i].fSigmaY2);
		      //HLTInfo( "padw0: %f", padw );
		      Float_t timew = sqrt( cl[i].fSigmaZ2 );
		      
		      //Check for saturation in the widths.
		      //Basically only store a certain number of decimals here, and cut the widths which is higher:
		      if(padw >= (1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor())
			padw = (1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor() - 1/AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor();
		      //HLTInfo( "padw1: %f", padw );
		      if(timew >= (1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor())
			timew = (1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor() - 1/AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor();;
		      
		      //Write pad
		      Int_t buff;
		      buff = (Int_t)rint(cl[i].fPad*AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor());
		      if(buff<0)
			{
			  HLTError( "Wrong pad value %d (%f, %f, row %u, i: %u)",buff, cl[i].fPad, AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor(), (unsigned)thisRow->fNClusters, (unsigned)i );
			  return EINVAL;
			}
		      if ( !OutputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNPadBitsRemaining()) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			outputSize = GetBitDataOutputSizeBytes();
			return ENOBUFS;
			}
		      
		    
		      //Write time
		      buff = (Int_t)rint(cl[i].fTime*AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor());
		      if(buff<0)
			{
			  HLTError( "Wrong time value %d",buff);
			return EINVAL;
			}
		      if ( !OutputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNTimeBitsRemaining()) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		      
		      //Write widths
		      buff = (Int_t)rint(padw*AliHLTTPCCompDataCompressorHelper::GetPadPrecisionFactor());
		      HLTDebug( "padw/buff: %d (%d / 0x%08X)", buff, 
				(buff & ((1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining())-1)),
				(buff & ((1<<AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining())-1)) );
		      
		      if ( !OutputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
			{
			HLTError( "Not enough space to write compressed data. %lu already written",
				  (unsigned long)GetBitDataOutputSizeBytes() );
			outputSize = GetBitDataOutputSizeBytes();
			return ENOBUFS;
			}
		      buff = (Int_t)rint(timew*AliHLTTPCCompDataCompressorHelper::GetTimePrecisionFactor());
		      if ( !OutputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNShapeBitsRemaining()) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		      
		      //Write charge 
		      buff = cl[i].fCharge;
		      if(buff >= 1<<(AliHLTTPCCompDataCompressorHelper::GetNChargeBits()))
			buff = (1<<(AliHLTTPCCompDataCompressorHelper::GetNChargeBits()))-1;
		      if ( !OutputBits(buff,AliHLTTPCCompDataCompressorHelper::GetNChargeBits()) )
			{
			  HLTError( "Not enough space to write compressed data. %lu already written",
				    (unsigned long)GetBitDataOutputSizeBytes() );
			  outputSize = GetBitDataOutputSizeBytes();
			  return ENOBUFS;
			}
		    }
		  inputPtr += sizeof(AliHLTTPCRemainingRow)+thisRow->fNClusters*sizeof(AliHLTTPCRemainingCluster);
		}
	    }
	  
	}
      CloseBitDataOutput();
      outputSize = GetBitDataOutputSizeBytes();
      return 0;
    }
