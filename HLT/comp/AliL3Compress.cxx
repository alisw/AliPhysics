// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliL3StandardIncludes.h"

#include "bitio.h"
#include "AliL3RootTypes.h"
#include "AliL3Models.h"
#include "AliL3DigitData.h"
#include "AliL3Logging.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3DataCompressor.h"

#if 0
#ifdef use_root
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#endif
#ifdef use_aliroot
#include "AliL3FileHandler.h"
#endif
#endif

#include "AliL3Compress.h"

#if GCCVERSION == 3
using namespace std;
#endif

//_____________________________________________________________
//
//  AliL3Compress
//
// Class for compressing and uncompressing data.

ClassImp(AliL3Compress)

AliL3Compress::AliL3Compress()
{
  fTracks=0;
  fSlice =0;
  fPatch=0;
  fWriteShape=kFALSE;
  fEvent=-1;
}

AliL3Compress::AliL3Compress(Int_t slice,Int_t patch,Char_t *path,Bool_t writeshape,Int_t event)
{
  fEvent=event;
  fSlice=slice;
  fPatch=patch;
  fTracks=0;
  sprintf(fPath,"%s",path);
  fWriteShape=writeshape;
}

AliL3Compress::~AliL3Compress()
{
  if(fTracks)
    delete fTracks;
}

Bool_t AliL3Compress::WriteFile(AliL3TrackArray *tracks,Char_t *filename)
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
      cerr<<"AliL3Compress::WriteFile : Error opening file "<<fname<<endl;
      return kFALSE;
    }
  Short_t ntracks = tracks->GetNTracks();
    
  Int_t count=0;
  AliL3ClusterModel *clusters=0;
  AliL3TrackModel *model=0;
  for(Int_t i=0; i<ntracks; i++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;

      //Do not save useless tracks or clusters:
      //if(track->GetNPresentClusters() == 0)
      //continue;
      
      track->FillModel();
      model = track->GetModel();
      //if(model->fNClusters==0) continue;
      clusters = track->GetClusters();
      if(fwrite(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      //if(fwrite(clusters,model->fNClusters*sizeof(AliL3ClusterModel),1,file)!=1) break;
      if(fwrite(clusters,AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel),1,file)!=1) break;
      count++;
      
    }
  fclose(file);
  return kTRUE;
}

Bool_t AliL3Compress::ReadFile(Char_t which,Char_t *filename)
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
	  cerr<<"AliL3Compress::ReadFile() : Wrong option"<<endl;
	  return kFALSE;
	}
    }

  FILE *file = fopen(fname,"r");
  if(!file)
    {
      cerr<<"AliL3Compress::ReadFile : Cannot open file "<<fname<<endl;
      return kFALSE;
    }

  if(fTracks)
    delete fTracks;
  fTracks = new AliL3TrackArray("AliL3ModelTrack");
  
  while(!feof(file))
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->NextTrack();
      track->Init(fSlice,fPatch);
      AliL3TrackModel *model = track->GetModel();
      AliL3ClusterModel *clusters = track->GetClusters();
      if(fread(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      //if(fread(clusters,model->fNClusters*sizeof(AliL3ClusterModel),1,file)!=1) break;
      if(fread(clusters,AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel),1,file)!=1) break;
      track->FillTrack();
    }

  fTracks->RemoveLast();
  fclose(file);
  return kTRUE;
}

Bool_t AliL3Compress::CompressFile()
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
      cerr<<"AliL3Compress::CompressFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliL3TrackModel track;
  AliL3ClusterModel cluster;
  Int_t temp;
  Int_t power;
  
  Int_t timeo,pado,chargeo,shapeo;
  timeo=pado=chargeo=shapeo=0;
  while(!feof(input))
    {
      if(fread(&track,sizeof(AliL3TrackModel),1,input)!=1) break;
      
      if(output->mask != 0x80) //Write the current byte to file.
	{
	  //cerr<<"\nAliL3Compress::CompressFile() : Writing overhead bits!!!"<<endl;
	  if(putc(output->rack,output->file )!=output->rack)
	    cerr<<"AliL3Compress::ComressFile : Error writing to bitfile"<<endl;
	  output->mask=0x80;
	  output->rack=0;
	}
      
      //Write track parameters:
      fwrite(&track,sizeof(AliL3TrackModel),1,output->file);
      
      Int_t origslice=-1,slice,clustercount=0;
      for(Int_t i=0; i<AliL3Transform::GetNRows(fPatch); i++)
	{
	  if(fread(&cluster,sizeof(AliL3ClusterModel),1,input)!=1) break;
	  
	  //Write empty flag:
	  temp = (Int_t)cluster.fPresent;
	  OutputBit(output,temp);
	  if(!temp) continue;
	  
	  if(cluster.fSlice<0 || cluster.fSlice>35)
	    {
	      cerr<<"AliL3DataCompress::CompressFile : Fucked up slice number :"<<cluster.fSlice<<endl;
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
	  power = 1<<(AliL3DataCompressor::GetNTimeBits()-1);
	  if(abs(temp)>=power)
	    {
	      timeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliL3DataCompressor::GetNTimeBits()-1));
	  
	  //Write pad information:
	  temp = (Int_t)rint(cluster.fDPad);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<(AliL3DataCompressor::GetNPadBits()-1);
	  if(abs(temp)>=power)
	    {
	      pado++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliL3DataCompressor::GetNPadBits()-1));
	  
	  //Write charge information:
	  temp = (Int_t)cluster.fDCharge;
	  power = 1<<(AliL3DataCompressor::GetNChargeBits());
	  if(abs(temp)>=power)
	    {
	      chargeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliL3DataCompressor::GetNChargeBits()));
	  
	  if(fWriteShape)
	    {
	      //Write shape information:
	      temp = (Int_t)cluster.fDSigmaY2;
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliL3DataCompressor::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  shapeo++;
		  temp = power - 1;
		}
	      OutputBits(output,abs(temp),(AliL3DataCompressor::GetNShapeBits()-1));
	      
	      temp = (Int_t)cluster.fDSigmaZ2;
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliL3DataCompressor::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  shapeo++;
		  temp=power - 1;
		}
	      OutputBits(output,abs(temp),(AliL3DataCompressor::GetNShapeBits()-1));
	    }
	  
	  clustercount++;
	}
    }
  
  fclose(input);
  CloseOutputBitFile(output);
  if(pado || timeo || chargeo || shapeo)
    {
      cout<<endl<<"Saturations: "<<endl
	  <<"Pad "<<pado<<endl
	  <<"Time "<<timeo<<endl
	  <<"Charge "<<chargeo<<endl
	  <<"Shape "<<shapeo<<endl<<endl;
    }
  return kTRUE;
}

Bool_t AliL3Compress::ExpandFile()
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
      cerr<<"AliL3Compress::ExpandFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliL3TrackModel trackmodel;
  AliL3ClusterModel *clusters=0;
  Int_t count=0;
  
  clusters = new AliL3ClusterModel[(AliL3Transform::GetNRows(fPatch))];
  while(!feof(input->file))
    {
      input->mask=0x80;//make sure we read a new byte from file.
      
      //Read and write track:
      if(fread(&trackmodel,sizeof(AliL3TrackModel),1,input->file)!=1) break;
      fwrite(&trackmodel,sizeof(AliL3TrackModel),1,output);

      memset(clusters,0,AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel));
      Int_t origslice=-1,clustercount=0;
      for(Int_t i=0; i<AliL3Transform::GetNRows(fPatch); i++)
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
	  temp = InputBits(input,(AliL3DataCompressor::GetNTimeBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDTime = temp;
	  
	  //Read pad information:
	  sign=InputBit(input);
	  temp = InputBits(input,(AliL3DataCompressor::GetNPadBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDPad = temp;
	  
	  //Read charge information:
	  temp=InputBits(input,(AliL3DataCompressor::GetNChargeBits()));
	  clusters[i].fDCharge = temp;
	  
	  if(fWriteShape)
	    {
	      //Read shape information:
	      sign = InputBit(input);
	      temp = InputBits(input,(AliL3DataCompressor::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaY2 = temp;
	      
	      sign = InputBit(input);
	      temp = InputBits(input,(AliL3DataCompressor::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaZ2 = temp;
	    }
	  clustercount++;
	}
      count++;
      //fwrite(clusters,(trackmodel.fNClusters)*sizeof(AliL3ClusterModel),1,output);
      fwrite(clusters,AliL3Transform::GetNRows(fPatch)*sizeof(AliL3ClusterModel),1,output);
    }
  
  delete [] clusters;
  fclose(output);
  CloseInputBitFile(input);
  return kTRUE;
}

void AliL3Compress::PrintCompRatio(ofstream *outfile)
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
    sprintf(fname,"%s/comp/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_c_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);

  mem->SetBinaryInput(fname);
  UInt_t compress_size = mem->GetFileSize();
  mem->CloseBinaryInput();
  
  if(digit_size==0)
    {
      cerr<<"AliL3Compress::PrintCompRatio : Zero digit size, not able to obtain comp. ratios!"<<endl;
      return;
    }
  
  Float_t compratio = (Float_t)(compress_size + remain_size)/(Float_t)digit_size;
  if(outfile)
    {
      ofstream &out = *outfile;
      out<<compress_size<<' '<<remain_size<<' '<<digit_size<<endl;
    }

  cout<<"=========================================="<<endl;
  cout<<"Original digits size : "<<digit_size/1000<<" kByte ( 100 % )"<<endl;
  cout<<"Compressed file size : "<<compress_size/1000<<" kByte ( "<<(Float_t)compress_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"Remainig file size   : "<<remain_size/1000<<" kByte ( "<<(Float_t)remain_size*100/(Float_t)digit_size<<" % )"<<endl;
  cout<<"---------------------- "<<endl;
  cout<<"Compression ratio    : "<<compratio*100<<" %"<<endl;
  cout<<"=========================================="<<endl;
}
