// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#if __GNUC__ == 3
using namespace std;
#endif

#include "AliL3StandardIncludes.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3Compress.h"
#include "AliL3DataCompressorHelper.h"
#include "AliL3CompressAC.h"


//_____________________________________________________________
//
//  AliL3CompressAC
//
// Compression class which performs Arithmetic Coding of the quantized residuals.
// The implemented algorithm is inspired by the examples in The Data Compression Book 
// by Nelson & Gailly.

ClassImp(AliL3CompressAC)

AliL3CompressAC::AliL3CompressAC()
{
  fCount=0;
  fTotals=0;
  fMax=0;
  fRange=0;
  fLow=0;
  fHigh=0;
  fUnderflowBits=0;
  fCode=0;
}

AliL3CompressAC::AliL3CompressAC(Int_t slice,Int_t patch,Char_t *path,Bool_t writeshape,Int_t event) :
  AliL3Compress(slice,patch,path,writeshape,event)
{
  fCount=0;
  fTotals=0;
  fMax=0;
  fRange=0;
  fLow=0;
  fHigh=0;
  fUnderflowBits=0;
  fCode=0;
}

AliL3CompressAC::~AliL3CompressAC()
{
  ClearArrays();
}

void AliL3CompressAC::ClearArrays()
{
  fMax=0;
  if(fCount)
    delete [] fCount;
  if(fTotals)
    delete [] fTotals;
}

void AliL3CompressAC::BuildModel(BIT_FILE *output)
{
  //Build the model from the input data, i.e. probability distributions of the quantized residuals.
  
  ClearArrays();
  ReadFile('m');
  
  UInt_t nmax=10000,qres;
  UInt_t * temp = new UInt_t[nmax];
  memset(&temp[0],0,nmax*sizeof(UInt_t));

  AliL3TrackArray *tracks = GetTracks();
  for(Int_t t=0; t<tracks->GetNTracks(); t++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(t);
      if(!track) continue;
      for(Int_t padrow=0; padrow<AliL3Transform::GetNRows(); padrow++)
	{
	  if(!track->IsPresent(padrow)) continue;
	  qres = abs((Int_t)rint(track->GetClusterModel(padrow)->fDPad));
	  if(qres >= nmax)
	    {
	      cerr<<"AliL3CompressAC::BuildModel() : Residual values seems way too big!"<<endl;
	      continue;
	    }
	  if(qres > fMax)
	    fMax=qres;
	  temp[qres]++;
	  qres = abs((Int_t)rint(track->GetClusterModel(padrow)->fDTime));
	  if(qres > fMax)
	    fMax = qres;
	  temp[qres]++;

	}
    }
  
  fCount = new UChar_t[fMax+1];
  
  //Find the highest counts in order to do scaling:
  UInt_t i,max_count=0;
  for(i=0; i<=fMax; i++)
    {
      if(temp[i] > max_count)
	max_count=temp[i];
    }

  //Perform the scaling
  UInt_t scale,total=1;
  scale = max_count / 256 + 1;
  for(i=0; i<=fMax; i++)
    {
      fCount[i] = (UChar_t)(temp[i]/scale);
      if(fCount[i]==0 && temp[i]!=0)
	fCount[i] = 1;
      total += (UInt_t)fCount[i];
    }
  if(total > (32767 - 256))
    scale=4;
  else if(total > 16383)
    scale=2;
  else
    scale=1;
  if(scale > 1)
    for(i=0; i<=fMax; i++)
      fCount[i] /= scale;

  cout<<"Writing "<<sizeof(UChar_t)*fMax+1<<" bytes with model information to compressed file"<<endl;
  fwrite(&fMax,sizeof(UShort_t),1,output->file);
  fwrite(fCount,sizeof(UChar_t),fMax+1,output->file);
  
  FillTotals();
  delete [] temp;
}

void AliL3CompressAC::RebuildModel(BIT_FILE *input)
{
  //Rebuild the model from the counts written to the beginning of the compressed file.
  
  ClearArrays();
  fread(&fMax,sizeof(UShort_t),1,input->file);
  fCount = new UChar_t[fMax+1];
  fread(fCount,sizeof(UChar_t),fMax+1,input->file);
  FillTotals();
}

void AliL3CompressAC::FillTotals()
{
  //Fill the array of totals, which is actually the model being used during encoding/decoding.
  if(fMax == 0)
    cerr<<"AliL3CompressAC::FillTotals : max value is zero!"<<endl;

  fTotals = new UInt_t[fMax+3];//up to max, and one reserved for endofstream symbol
  
  UInt_t i;
  fTotals[0]=0;
  for(i=0; i<=fMax; i++)
    {
      fTotals[i+1] = fTotals[i] + fCount[i];
    }
  fTotals[fMax+2] = fTotals[fMax+1]+1;//Used for the scale
}

void AliL3CompressAC::PrintTotals()
{
  cout<<"Totals:"<<endl;
  for(UInt_t i=0; i<=fMax; i++)
    {
      cout<<"Totals "<<i<<" "<<fTotals[i]<<" count "<<(int)fCount[i]<<endl;
    }
}

void AliL3CompressAC::InitEncoder()
{
  fLow = 0;
  fHigh = 0xffff;
  fUnderflowBits=0;
}

void AliL3CompressAC::InitDecoder(BIT_FILE *input)
{
  fCode=0;
  for(Int_t i=0; i<16; i++)
    {
      fCode <<= 1;
      fCode += InputBit(input);
    }
  fLow = 0;
  fHigh = 0xffff;
}

void AliL3CompressAC::ConvertIntToSymbol(Int_t value)
{
  UInt_t range = fHigh - fLow + 1;
  fHigh = fLow + (UShort_t)((range*fTotals[value+1])/fTotals[fMax+2] - 1);
  fLow = fLow + (UShort_t)((range*fTotals[value])/fTotals[fMax+2]);
}

UInt_t AliL3CompressAC::ConvertSymbolToInt()
{
  UInt_t range = (UInt_t)(fHigh-fLow) + 1;
  UShort_t count = (UShort_t)((((UInt_t)(fCode-fLow)+1)*fTotals[fMax+2] - 1)/range);
  UInt_t j=fMax+1;
  while(count < fTotals[j])
    j--;
  
  return j;
}

void AliL3CompressAC::EncodeSymbol(BIT_FILE *output)
{
  while(1)
    {
      if( (fHigh & 0x8000) == (fLow & 0x8000) )
	{
	  OutputBit(output,fHigh & 0x8000);
	  while(fUnderflowBits > 0)
	    {
	      OutputBit(output,~fHigh & 0x8000);
	      fUnderflowBits--;
	    }
	}
      else if( (fLow & 0x4000) && !(fHigh & 0x4000) )
	{
	  fUnderflowBits++;
	  fLow &= 0x3fff;
	  fHigh |= 0x4000;
	}
      else
	return;
      fLow <<= 1;
      fHigh <<= 1;
      fHigh |= 1;
    }
}

void AliL3CompressAC::RemoveSymbolFromStream(BIT_FILE *input,Int_t j)
{
  UInt_t range = (UInt_t)(fHigh-fLow)+1;
  fHigh = fLow + (UShort_t)((range*fTotals[j+1])/fTotals[fMax+2]-1);
  fLow = fLow + (UShort_t)((range*fTotals[j])/fTotals[fMax+2]);
  while(1)
    {
      if( (fHigh & 0x8000)==(fLow & 0x8000) )
	{}
      else if((fLow & 0x4000) == 0x4000 && (fHigh & 0x4000)==0)
	{
	  fCode ^= 0x4000;
	  fLow &= 0x3fff;
	  fHigh |= 0x4000;
	}
      else
	return;
      fLow <<= 1;
      fHigh <<= 1;
      fHigh |= 1;
      fCode <<= 1;
      fCode += InputBit(input);
    }
}

void AliL3CompressAC::FlushEncoder(BIT_FILE *output)
{
  //Flush the encoder:
  OutputBit(output,fLow & 0x4000);
  fUnderflowBits++;
  while(fUnderflowBits-- > 0)
    OutputBit(output,~fLow & 0x4000);
  
}


Bool_t AliL3CompressAC::CompressFile()
{
  Char_t fname[100];
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_ac_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_ac_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  BIT_FILE *output = OpenOutputBitFile(fname);
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_m_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  
  FILE *input = fopen(fname,"r");
  if(!input)
    {
      cerr<<"AliL3CompressAC::CompressFileAC() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }
  
  BuildModel(output);

  AliL3TrackModel track;
  Int_t temp,power,i,j;
  
  fseek(input,0,SEEK_END);
  UInt_t size = ftell(input);
  rewind(input);
  Int_t trackcount = size/(sizeof(AliL3TrackModel) + sizeof(AliL3ClusterModel)*AliL3Transform::GetNRows(fPatch));
  
  //Write the number of tracks in the beginning of stream.
  fwrite(&trackcount,sizeof(Int_t),1,output->file);
  
  AliL3ClusterModel **clusters = new AliL3ClusterModel*[trackcount];
  Int_t *clustercount = new Int_t[trackcount];
  i=0;
  
  //Read all the tracks from input file, and write them all to the outputfile.
  //Store the clusters in memory for later encoding and storing.
  while(!feof(input))
    {
      if(fread(&track,sizeof(AliL3TrackModel),1,input)!=1) break;
      fwrite(&track,sizeof(AliL3TrackModel),1,output->file);
      
      clusters[i] = new AliL3ClusterModel[AliL3Transform::GetNRows()];
      clustercount[i]=0;
      
      //Read in the clusters:
      fread(clusters[i],sizeof(AliL3ClusterModel),AliL3Transform::GetNRows(fPatch),input);
      i++;
    }
  if(i != trackcount)
    {
      cerr<<"AliL3CompressAC::CompressFile : Mismatching file size and trackcount "<<i<<" "<<trackcount<<endl;
      exit(5);
    }
  fclose(input);
  
  //Write all the fixed size variables of the clusters:
  for(i=0; i<trackcount; i++)
    {
      Int_t origslice=-1,slice;
      for(j=0; j<AliL3Transform::GetNRows(fPatch); j++)
	{
	  temp = (Int_t)clusters[i][j].fPresent;
	  OutputBit(output,temp);
	  if(!temp) continue;
	  
	  if(clusters[i][j].fSlice<0 || clusters[i][j].fSlice>35)
	    {
	      cerr<<"AliL3DataCompress::CompressFile : Fucked up slice number :"<<clusters[i][j].fSlice<<endl;
	      exit(5);
	    }
	  
	  //Write slice number of first point
	  if(clustercount[i]==0)
	    {
	      origslice = clusters[i][j].fSlice;
	      OutputBits(output,origslice,6); //Need 6 bits to encode slice number
	    }
	  else
	    {
	      slice = clusters[i][j].fSlice;
	      if(slice == origslice)
		{
		  OutputBit(output,0);  //No change of slice
		}
	      else
		{
		  OutputBit(output,1);
		  if(abs(slice - origslice)==1)
		    {
		      if(slice > origslice)
			OutputBit(output,1);
		      else
			OutputBit(output,0);
		    }
		  else
		    {
		      if( (slice == 0 && origslice == 17) || (slice == 18 && origslice == 35) )
			OutputBit(output,1);
		      else if( (slice == 17 && origslice == 0) || (slice == 35 && origslice == 18) )
			OutputBit(output,0);
		    }
		  origslice=slice;
		}
	    }
	  
	  //Write charge information:
	  temp = (Int_t)clusters[i][j].fDCharge;
	  power = 1<<(AliL3DataCompressorHelper::GetNChargeBits());
	  if(abs(temp)>=power)
	    {
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliL3DataCompressorHelper::GetNChargeBits()));

	  //Write sign information of the residuals:
	  temp = (Int_t)rint(clusters[i][j].fDTime);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  temp = (Int_t)rint(clusters[i][j].fDPad);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);

	  //Write shape information if requested:
	  if(fWriteShape)
	    {
	      temp = (Int_t)rint(clusters[i][j].fDSigmaY);
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliL3DataCompressorHelper::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  temp = power - 1;
		}
	      OutputBits(output,abs(temp),(AliL3DataCompressorHelper::GetNShapeBits()-1));
	      
	      temp = (Int_t)rint(clusters[i][j].fDSigmaZ);
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliL3DataCompressorHelper::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  temp=power - 1;
		}
	      OutputBits(output,abs(temp),(AliL3DataCompressorHelper::GetNShapeBits()-1));
	    }
	  clustercount[i]++;
	}
    }

  InitEncoder();
  
  //Start the arithmetic coding of the residuals.
  //All the residuals (both pad and time) are coded in one go,
  //i.e. for all tracks and clusters in the input file.
  for(i=0; i<trackcount; i++)
    {
      Int_t counter=0;
      
      for(j=0; j<AliL3Transform::GetNRows(fPatch); j++)
	{
	  if(!clusters[i][j].fPresent) continue;
	  temp = abs((Int_t)rint(clusters[i][j].fDTime));
	  ConvertIntToSymbol(temp);
	  EncodeSymbol(output);

	  temp = abs((Int_t)rint(clusters[i][j].fDPad));
	  ConvertIntToSymbol(temp);
 	  EncodeSymbol(output);
	  counter++;
	}
      if(counter != clustercount[i])
	{
	  cerr<<"AliL3CompressAC::CompressFile : Mismatching clustercount "<<counter<<" "<<clustercount[i]<<endl;
	  exit(5);
	}

    }
  
  ConvertIntToSymbol(fMax+1);//End of stream symbol
  EncodeSymbol(output);
  FlushEncoder(output);
  OutputBits(output,0,16);

  for(i=0; i<trackcount; i++)
    delete [] clusters[i];
  delete [] clusters;
  delete [] clustercount;
  

  CloseOutputBitFile(output);
  
  return kTRUE;
}

Bool_t AliL3CompressAC::ExpandFile()
{
  Char_t fname[100];
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_ac_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_ac_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  BIT_FILE *input = OpenInputBitFile(fname);
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_u_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_u_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  FILE *output = fopen(fname,"w");
  if(!output)
    {
      cerr<<"AliL3Compress::ExpandFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }
  
  RebuildModel(input);

  int trackcount,i,j;
  fread(&trackcount,sizeof(Int_t),1,input->file);

  AliL3TrackModel *trackmodels = new AliL3TrackModel[trackcount];
  AliL3ClusterModel **clusters = new AliL3ClusterModel*[trackcount];
  Int_t *clustercount = new Int_t[trackcount];

  fread(trackmodels,sizeof(AliL3TrackModel),trackcount,input->file);
  
  for(i=0; i<trackcount; i++)
    {
      clusters[i] = new AliL3ClusterModel[AliL3Transform::GetNRows(fPatch)];
      clustercount[i]=0;

      //Read the fixed size variables:
      Int_t origslice=-1;
      Int_t temp,sign;
      for(j=0; j<AliL3Transform::GetNRows(fPatch); j++)
	{
	  //Read empty flag:
	  temp = InputBit(input);
	  if(!temp) 
	    {
	      clusters[i][j].fPresent=kFALSE;
	      continue;
	    }
	  clusters[i][j].fPresent=kTRUE;
	  
	  //Read slice information
	  if(clustercount[i]==0)
	    {
	      temp = InputBits(input,6);
	      clusters[i][j].fSlice = temp;
	      origslice = temp;
	    }
	  else
	    {
	      temp = InputBit(input);
	      if(!temp)//no change
		clusters[i][j].fSlice = origslice;
	      else
		{
		  temp = InputBit(input);
		  if(temp)
		    {
		      if(origslice == 17)
			origslice = 0;
		      else if(origslice == 35)
			origslice = 18;
		      else
			origslice++;
		    }
		  else
		    {
		      if(origslice == 0)
			origslice = 17;
		      else if(origslice == 18)
			origslice = 35;
		      else
			origslice--;
		    }
		  if(origslice < 0 || origslice > 35)
		    {
		      cerr<<"AliL3CompressAC::ExpandFile : Bad slice number "<<temp<<endl;
		      exit(5);
		    }
		  clusters[i][j].fSlice = origslice;
		}
	    }
	  
	  //Read charge information:
	  temp=InputBits(input,(AliL3DataCompressorHelper::GetNChargeBits()));
	  clusters[i][j].fDCharge = temp;
	  
	  //Read sign information of the residuals:
	  sign=InputBit(input);
	  if(!sign)
	    clusters[i][j].fDTime = -1;
	  else
	    clusters[i][j].fDTime = 1;
	  sign=InputBit(input);
	  if(!sign)
	    clusters[i][j].fDPad = -1;
	  else
	    clusters[i][j].fDPad = 1;
	  
	  //Read shape information if requested
	  if(fWriteShape)
	    {
	      sign = InputBit(input);
	      temp = InputBits(input,(AliL3DataCompressorHelper::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i][j].fDSigmaY = temp;
	      
	      sign = InputBit(input);
	      temp = InputBits(input,(AliL3DataCompressorHelper::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i][j].fDSigmaZ = temp;
	    }
	  clustercount[i]++;
	}
    }
  
  InitDecoder(input);
  
  Int_t temp;
  for(i=0; i<trackcount; i++)
    {
      Int_t count=0;
      for(j=0; j<AliL3Transform::GetNRows(fPatch); j++)
	{
	  if(!clusters[i][j].fPresent) continue;
	    
	  temp = ConvertSymbolToInt();
	  RemoveSymbolFromStream(input,temp);
	  clusters[i][j].fDTime *= temp;
	  
	  temp = ConvertSymbolToInt();
	  RemoveSymbolFromStream(input,temp);
	  clusters[i][j].fDPad *= temp;
	  count++;
	}

      if(count != clustercount[i])
	{
	  cerr<<"AliL3CompressAC::ExpandFile : Mismatching clustercount "<<count<<" "<<clustercount[i]<<endl;
	  exit(5);
	}
    }
  
  //Now there should be a endofstream indicator, if not something went wrong during encoding/decoding.
  temp = ConvertSymbolToInt();
  if((UShort_t)temp != fMax + 1)
    cerr<<"AliL3CompressAC::ExpandFile : Missing the endofstream indicator!"<<endl;
  
  CloseInputBitFile(input);

  //Write everything to the uncompressed outfile:
  for(i=0; i<trackcount; i++)
    {
      fwrite(&trackmodels[i],sizeof(AliL3TrackModel),1,output);
      fwrite(clusters[i],AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel),1,output);
    }
  
  fclose(output);
  
  for(i=0; i<trackcount; i++)
    delete [] clusters[i];
  delete [] clusters;
  delete [] trackmodels;
  delete [] clustercount;
  
  return kTRUE;
}

void AliL3CompressAC::PrintCompRatio(ofstream *outfile)
{
  AliL3MemHandler *mem = new AliL3MemHandler();
  Char_t fname[1024];
  UInt_t remain_size=0,digit_size=0;
  for(Int_t i=0; i<36; i++)
    {
      if(fEvent<0)
	sprintf(fname,"%s/comp/remains_%d_%d.raw",fPath,i,-1);
      else
	sprintf(fname,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,i,-1);
      mem->SetBinaryInput(fname);
      remain_size += mem->GetFileSize();
      mem->CloseBinaryInput();

      sprintf(fname,"%s/binaries/digits_c8_%d_%d_%d.raw",fPath,fEvent,i,-1);
      mem->SetBinaryInput(fname);
      digit_size += mem->GetFileSize();
      mem->CloseBinaryInput();
    }
  
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_ac_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_ac_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);

  mem->SetBinaryInput(fname);
  UInt_t compress_size = mem->GetFileSize();
  mem->CloseBinaryInput();
  
  if(digit_size==0)
    {
      cerr<<"AliL3Compress::PrintCompRatio : Zero digit size, not able to obtain comp. ratios!"<<endl;
      return;
    }
  
  Float_t compratio = (Float_t)(compress_size + remain_size)/(Float_t)digit_size;
  Float_t entropy[3];
  Int_t track_size = GetEntropy(entropy[0],entropy[1],entropy[2])*sizeof(AliL3TrackModel);
  if(outfile)
    {
      ofstream &out = *outfile;
      out<<compress_size<<' '<<remain_size<<' '<<digit_size<<' '<<track_size<<' '<<entropy[0]<<' '<<entropy[1]<<endl;
    }
  
  cout<<"=========================================="<<endl;
  cout<<"Original digits size : "<<digit_size/1000<<" kByte ( 100 % )"<<endl;
  cout<<"Compressed file size : "<<compress_size/1000<<" kByte ( "<<(Float_t)compress_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"Remaining file size  : "<<remain_size/1000<<" kByte ( "<<(Float_t)remain_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"Relative track size  : "<<track_size/1000<<" kByte ( "<<(Float_t)track_size*100/(Float_t)compress_size<<" % )"<<endl;
  cout<<"Relative cluster size: "<<(compress_size-track_size)/1000<<" kByte ( "<<(Float_t)(compress_size-track_size)*100/(Float_t)compress_size<<" % )"<<endl;
  cout<<"---------------------- "<<endl;
  cout<<"Compression ratio    : "<<compratio*100<<" %"<<endl;
  cout<<"=========================================="<<endl;
  cout<<"Entropy of residual and charge : "<<entropy[0]<<" "<<entropy[1]<<" "<<entropy[2]<<endl;
}

