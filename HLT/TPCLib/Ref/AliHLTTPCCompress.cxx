// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTTPCStandardIncludes.h"

#include "bitio.h"
#include "AliHLTTPCRootTypes.h"
#include "AliHLTTPCModels.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCLogging.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCModelTrack.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDataCompressor.h"

#if 0
#ifdef use_root
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#endif
#ifdef use_aliroot
#include "AliHLTTPCFileHandler.h"
#endif
#endif

#include "AliHLTTPCCompress.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliHLTTPCCompress
//
// Class for compressing and uncompressing data.

ClassImp(AliHLTTPCCompress)

AliHLTTPCCompress::AliHLTTPCCompress()
{
  fTracks=0;
  fSlice =0;
  fPatch=0;
  fWriteShape=kFALSE;
  fEvent=-1;
}

AliHLTTPCCompress::AliHLTTPCCompress(Int_t slice,Int_t patch,Char_t *path,Bool_t writeshape,Int_t event)
{
  fEvent=event;
  fSlice=slice;
  fPatch=patch;
  fTracks=0;
  sprintf(fPath,"%s",path);
  fWriteShape=writeshape;
}

AliHLTTPCCompress::~AliHLTTPCCompress()
{
  if(fTracks)
    delete fTracks;
}

Bool_t AliHLTTPCCompress::WriteFile(AliHLTTPCTrackArray *tracks,Char_t *filename)
{
  Char_t fname[1024];
  if(filename)
    sprintf(fname,"%s/comp/%s",fPath,filename);
  else if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_m_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);

  FILE *file = fopen(fname,"w");
  if(!file)
    {
      cerr<<"AliHLTTPCCompress::WriteFile : Error opening file "<<fname<<endl;
      return kFALSE;
    }
  Short_t ntracks = tracks->GetNTracks();
    
  Int_t count=0;
  AliHLTTPCClusterModel *clusters=0;
  AliHLTTPCTrackModel *model=0;
  for(Int_t i=0; i<ntracks; i++)
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;

      //Do not save useless tracks or clusters:
      //if(track->GetNPresentClusters() == 0)
      //continue;
      
      track->FillModel();
      model = track->GetModel();
      //if(model->fNClusters==0) continue;
      clusters = track->GetClusters();
      if(fwrite(model,sizeof(AliHLTTPCTrackModel),1,file)!=1) break;
      //if(fwrite(clusters,model->fNClusters*sizeof(AliHLTTPCClusterModel),1,file)!=1) break;
      if(fwrite(clusters,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(AliHLTTPCClusterModel),1,file)!=1) break;
      count++;
      
    }
  fclose(file);
  return kTRUE;
}

Bool_t AliHLTTPCCompress::ReadFile(Char_t which,Char_t *filename)
{
  //Read the trackfile.

  Char_t fname[1024];
  if(filename)
    sprintf(fname,"%s/comp/%s",fPath,filename);
  else
    {
      if(which == 'm')
	{
	  if(fEvent<0)
	    sprintf(fname,"%s/comp/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
	  else
	    sprintf(fname,"%s/comp/tracks_m_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
	}
      else if(which == 'u')
	{
	  if(fEvent<0)
	    sprintf(fname,"%s/comp/tracks_u_%d_%d.raw",fPath,fSlice,fPatch);
	  else
	    sprintf(fname,"%s/comp/tracks_u_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
	}
      else
	{
	  cerr<<"AliHLTTPCCompress::ReadFile() : Wrong option"<<endl;
	  return kFALSE;
	}
    }

  FILE *file = fopen(fname,"r");
  if(!file)
    {
      cerr<<"AliHLTTPCCompress::ReadFile : Cannot open file "<<fname<<endl;
      return kFALSE;
    }

  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTPCTrackArray("AliHLTTPCModelTrack");
  
  while(!feof(file))
    {
      AliHLTTPCModelTrack *track = (AliHLTTPCModelTrack*)fTracks->NextTrack();
      track->Init(fSlice,fPatch);
      AliHLTTPCTrackModel *model = track->GetModel();
      AliHLTTPCClusterModel *clusters = track->GetClusters();
      if(fread(model,sizeof(AliHLTTPCTrackModel),1,file)!=1) break;
      //if(fread(clusters,model->fNClusters*sizeof(AliHLTTPCClusterModel),1,file)!=1) break;
      if(fread(clusters,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(AliHLTTPCClusterModel),1,file)!=1) break;
      track->FillTrack();
    }

  fTracks->RemoveLast();
  fclose(file);
  return kTRUE;
}

Bool_t AliHLTTPCCompress::CompressFile()
{
  Char_t fname[100];
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_c_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  BIT_FILE *output = OpenOutputBitFile(fname);
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_m_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  
  FILE *input = fopen(fname,"r");
  if(!input)
    {
      cerr<<"AliHLTTPCCompress::CompressFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliHLTTPCTrackModel track;
  AliHLTTPCClusterModel cluster;
  Int_t temp;
  Int_t power;
  
  Int_t timeo,pado,chargeo,shapeo;
  timeo=pado=chargeo=shapeo=0;
  while(!feof(input))
    {
      if(fread(&track,sizeof(AliHLTTPCTrackModel),1,input)!=1) break;
      
      if(output->mask != 0x80) //Write the current byte to file.
	{
	  //cerr<<"\nAliHLTTPCCompress::CompressFile() : Writing overhead bits!!!"<<endl;
	  if(putc(output->rack,output->file )!=output->rack)
	    cerr<<"AliHLTTPCCompress::ComressFile : Error writing to bitfile"<<endl;
	  output->mask=0x80;
	  output->rack=0;
	}
      
      //Write track parameters:
      fwrite(&track,sizeof(AliHLTTPCTrackModel),1,output->file);
      
      Int_t origslice=-1,slice,clustercount=0;
      for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(fPatch); i++)
	{
	  if(fread(&cluster,sizeof(AliHLTTPCClusterModel),1,input)!=1) break;
	  
	  //Write empty flag:
	  temp = (Int_t)cluster.fPresent;
	  OutputBit(output,temp);
	  if(!temp) continue;
	  
	  if(cluster.fSlice<0 || cluster.fSlice>35)
	    {
	      cerr<<"AliHLTTPCDataCompress::CompressFile : Fucked up slice number :"<<cluster.fSlice<<endl;
	      exit(5);
	    }
	  
	  //Write slice number of first point
	  if(clustercount==0)
	    {
	      origslice = cluster.fSlice;
	      OutputBits(output,origslice,6); //Need 6 bits to encode slice number
	    }
	  else
	    {
	      slice = cluster.fSlice;
	      if(slice == origslice)
		OutputBit(output,0);
	      else
		{
		  OutputBit(output,1);
		  OutputBits(output,slice,6);
		  origslice=slice;
		}
	    }
	  
	  //Write time information:
	  temp = (Int_t)rint(cluster.fDTime);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<(AliHLTTPCDataCompressor::GetNTimeBits()-1);
	  if(abs(temp)>=power)
	    {
	      timeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTTPCDataCompressor::GetNTimeBits()-1));
	  
	  //Write pad information:
	  temp = (Int_t)rint(cluster.fDPad);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<(AliHLTTPCDataCompressor::GetNPadBits()-1);
	  if(abs(temp)>=power)
	    {
	      pado++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTTPCDataCompressor::GetNPadBits()-1));
	  
	  //Write charge information:
	  temp = (Int_t)cluster.fDCharge;
	  power = 1<<(AliHLTTPCDataCompressor::GetNChargeBits());
	  if(abs(temp)>=power)
	    {
	      chargeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTTPCDataCompressor::GetNChargeBits()));
	  
	  if(fWriteShape)
	    {
	      //Write shape information:
	      temp = (Int_t)cluster.fDSigmaY2;
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliHLTTPCDataCompressor::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  shapeo++;
		  temp = power - 1;
		}
	      OutputBits(output,abs(temp),(AliHLTTPCDataCompressor::GetNShapeBits()-1));
	      
	      temp = (Int_t)cluster.fDSigmaZ2;
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliHLTTPCDataCompressor::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  shapeo++;
		  temp=power - 1;
		}
	      OutputBits(output,abs(temp),(AliHLTTPCDataCompressor::GetNShapeBits()-1));
	    }
	  
	  clustercount++;
	}
    }
  
  fclose(input);
  CloseOutputBitFile(output);
  if(pado || timeo || chargeo || shapeo)
    {
#if 0
      cout<<endl<<"Saturations: "<<endl
	  <<"Pad "<<pado<<endl
	  <<"Time "<<timeo<<endl
	  <<"Charge "<<chargeo<<endl
	  <<"Shape "<<shapeo<<endl<<endl;
#endif
    }
  return kTRUE;
}

Bool_t AliHLTTPCCompress::ExpandFile()
{
  Char_t fname[100];
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_c_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  BIT_FILE *input = OpenInputBitFile(fname);
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_u_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_u_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);
  FILE *output = fopen(fname,"w");
  if(!output)
    {
      cerr<<"AliHLTTPCCompress::ExpandFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliHLTTPCTrackModel trackmodel;
  AliHLTTPCClusterModel *clusters=0;
  Int_t count=0;
  
  clusters = new AliHLTTPCClusterModel[(AliHLTTPCTransform::GetNRows(fPatch))];
  while(!feof(input->file))
    {
      input->mask=0x80;//make sure we read a new byte from file.
      
      //Read and write track:
      if(fread(&trackmodel,sizeof(AliHLTTPCTrackModel),1,input->file)!=1) break;
      fwrite(&trackmodel,sizeof(AliHLTTPCTrackModel),1,output);

      memset(clusters,0,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(AliHLTTPCClusterModel));
      Int_t origslice=-1,clustercount=0;
      for(Int_t i=0; i<AliHLTTPCTransform::GetNRows(fPatch); i++)
	{
	  Int_t temp,sign;
	  
	  //Read empty flag:
	  temp = InputBit(input);
	  if(!temp) 
	    {
	      clusters[i].fPresent=kFALSE;
	      continue;
	    }
	  clusters[i].fPresent=kTRUE;
	  
	  //Read slice information
	  if(clustercount==0)
	    {
	      temp = InputBits(input,6);
	      clusters[i].fSlice = temp;
	      origslice = temp;
	    }
	  else
	    {
	      temp = InputBit(input);
	      if(!temp)//no change
		clusters[i].fSlice = origslice;
	      else
		{
		  temp = InputBits(input,6);//read new slice
		  clusters[i].fSlice = temp;
		  origslice = temp;//store new slice
		}
	    }
	  
	  //Read time information:
	  sign=InputBit(input);
	  temp = InputBits(input,(AliHLTTPCDataCompressor::GetNTimeBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDTime = temp;
	  
	  //Read pad information:
	  sign=InputBit(input);
	  temp = InputBits(input,(AliHLTTPCDataCompressor::GetNPadBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDPad = temp;
	  
	  //Read charge information:
	  temp=InputBits(input,(AliHLTTPCDataCompressor::GetNChargeBits()));
	  clusters[i].fDCharge = temp;
	  
	  if(fWriteShape)
	    {
	      //Read shape information:
	      sign = InputBit(input);
	      temp = InputBits(input,(AliHLTTPCDataCompressor::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaY2 = temp;
	      
	      sign = InputBit(input);
	      temp = InputBits(input,(AliHLTTPCDataCompressor::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaZ2 = temp;
	    }
	  clustercount++;
	}
      count++;
      //fwrite(clusters,(trackmodel.fNClusters)*sizeof(AliHLTTPCClusterModel),1,output);
      fwrite(clusters,AliHLTTPCTransform::GetNRows(fPatch)*sizeof(AliHLTTPCClusterModel),1,output);
    }
  
  delete [] clusters;
  fclose(output);
  CloseInputBitFile(input);
  return kTRUE;
}

void AliHLTTPCCompress::PrintCompRatio(ofstream *outfile)
{
  AliHLTTPCMemHandler *mem = new AliHLTTPCMemHandler();
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
    sprintf(fname,"%s/comp/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_c_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);

  mem->SetBinaryInput(fname);
  UInt_t compress_size = mem->GetFileSize();
  mem->CloseBinaryInput();
  
  if(digit_size==0)
    {
      cerr<<"AliHLTTPCCompress::PrintCompRatio : Zero digit size, not able to obtain comp. ratios!"<<endl;
      return;
    }
  
  if(outfile)
    {
      ofstream &out = *outfile;
      out<<compress_size<<' '<<remain_size<<' '<<digit_size<<endl;
    }

#if 0
  Float_t compratio = (Float_t)(compress_size + remain_size)/(Float_t)digit_size;
  cout<<"=========================================="<<endl;
  cout<<"Original digits size : "<<digit_size/1000<<" kByte ( 100 % )"<<endl;
  cout<<"Compressed file size : "<<compress_size/1000<<" kByte ( "<<(Float_t)compress_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"Remainig file size   : "<<remain_size/1000<<" kByte ( "<<(Float_t)remain_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"---------------------- "<<endl;
  cout<<"Compression ratio    : "<<compratio*100<<" %"<<endl;
  cout<<"=========================================="<<endl;
#endif
}
