//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stdio.h>
#include <stream.h>
#include <stdlib.h>

#include "AliL3Compress.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "bitio.h"

ClassImp(AliL3Compress)

AliL3Compress::AliL3Compress()
{
  
}

AliL3Compress::~AliL3Compress()
{
  
}

void AliL3Compress::Write2File(AliL3TrackArray *tracks)
{
  FILE *file = fopen("data.raw","w");
  Short_t ntracks = tracks->GetNTracks();
  cout<<"Writing "<<ntracks<<" tracks to file"<<endl;
  //Write the number of tracks at the beginning:
  //fwrite(&ntracks,sizeof(Short_t),1,file);
  
  for(Int_t i=0; i<ntracks; i++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      AliL3TrackModel *model = track->GetModel();
      cout<<"Writing "<<model->fFirstPointX<<endl;
      if(fwrite(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      for(Int_t j=0; j<model->fNClusters; j++)
	{
	  AliL3ClusterModel *cl = track->GetClusterModel(j);
	  fwrite(cl,sizeof(AliL3ClusterModel),1,file);
	}
    }
  fclose(file);
}

void AliL3Compress::ReadFile()
{
  FILE *file = fopen("data.raw","r");
  
  AliL3TrackArray *tracks = new AliL3TrackArray("AliL3ModelTrack");
  Int_t ntracks=0;
  
  while(!feof(file))
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->NextTrack();
      track->Init(0,0);
      AliL3TrackModel *model = track->GetModel();
      AliL3ClusterModel *clusters = track->GetClusters();
      if(fread(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      cout<<"Read model "<<model->fFirstPointX<<endl;
      if(fread(clusters,(model->fNClusters)*sizeof(AliL3ClusterModel),1,file)!=1) break;
      ntracks++;
    }
  
  delete tracks;
  cout<<"Read "<<ntracks<<" tracks from file"<<endl;
  fclose(file);
}

void AliL3Compress::CompressFile()
{
  
  BIT_FILE *output = OpenOutputBitFile("test.raw");
  FILE *input = fopen("data.raw","r");

  AliL3TrackModel track;
  AliL3ClusterModel cluster;
  Int_t temp;

  while(!feof(input))
    {
      if(fread(&track,sizeof(AliL3TrackModel),1,input)!=1) break;
      cout<<"Writing "<<sizeof(AliL3TrackModel)<<endl;
      fwrite(&track,sizeof(AliL3TrackModel),1,output->file);
      Int_t bitcount=0;
      for(Int_t i=0; i<track.fNClusters; i++)
	{
	  if(fread(&cluster,sizeof(AliL3ClusterModel),1,input)!=1) break;
	  Int_t flag = cluster.fPresent;
	  OutputBit(output,flag);
	  bitcount++;
	  if(!flag) continue;
	  temp = (Int_t)cluster.fDTime;
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  OutputBits(output,abs(temp),8);
	  temp = (Int_t)cluster.fDPad;
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  OutputBits(output,abs(temp),8);
	  bitcount+=8;
	  temp = (Int_t)cluster.fDCharge;
	  OutputBits(output,temp,10);
	  bitcount+=10;
	  
	  //Short_t temp=(Short_t)cluster.fDTime;
	  cout<<"flag "<<(int)flag<<" dtime "<<(int)cluster.fDTime<<" dpad "<<(int)cluster.fDPad<<" charge "<<cluster.fDCharge<<endl;
	}
      
    }

  fclose(input);
  CloseOutputBitFile(output);
}

void AliL3Compress::ExpandFile()
{
  BIT_FILE *input = OpenInputBitFile("test.raw");
  
  AliL3TrackModel track;
  AliL3ClusterModel cluster;
  
  fread(&track,sizeof(AliL3TrackModel),1,input->file);
  for(Int_t i=0; i<track.fNClusters; i++)
    {
      Int_t temp,sign;
      temp = InputBit(input);
      if(!temp) break;
      sign=InputBit(input);
      temp = InputBits(input,8);
      if(!sign)
	temp*=-1;
      cout<<"Dtime "<<temp;
      sign=InputBit(input);
      temp = InputBits(input,8);
      if(!sign)
	temp*=-1;
      cout<<" DPad "<<temp;
      temp=InputBits(input,10);
      cout<<" Charge "<<temp<<endl;
    }
  CloseInputBitFile(input);
}
