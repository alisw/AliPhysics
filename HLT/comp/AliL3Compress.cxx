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
  fTracks=0;
}

AliL3Compress::~AliL3Compress()
{
  if(fTracks)
    delete fTracks;
}

void AliL3Compress::WriteFile(AliL3TrackArray *tracks,Char_t *filename)
{
  FILE *file = fopen(filename,"w");
  Short_t ntracks = tracks->GetNTracks();
  //cout<<"Writing "<<ntracks<<" tracks to file"<<endl;
    
  Int_t count=0;
  AliL3ClusterModel *clusters=0;
  AliL3TrackModel *model=0;
  for(Int_t i=0; i<ntracks; i++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
      model = track->GetModel();
      if(model->fNClusters==0) continue;
      clusters = track->GetClusters();
      //cout<<"Writing "<<(int)model->fNClusters<<" clusters"<<endl;
      if(fwrite(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      //cout<<"Writing "<<(int)model->fNClusters<<" clusters to file"<<endl;
      if(fwrite(clusters,model->fNClusters*sizeof(AliL3ClusterModel),1,file)!=1) break;
      count++;
    }
  cout<<"Wrote "<<count<<" tracks "<<endl;
  fclose(file);
}

void AliL3Compress::ReadFile(Char_t *filename)
{
  FILE *file = fopen(filename,"r");
  
  if(fTracks)
    delete fTracks;
  fTracks = new AliL3TrackArray("AliL3ModelTrack");
  
  while(!feof(file))
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->NextTrack();
      track->Init(0,0);
      AliL3TrackModel *model = track->GetModel();
      AliL3ClusterModel *clusters = track->GetClusters();
      if(fread(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      if(fread(clusters,(model->fNClusters)*sizeof(AliL3ClusterModel),1,file)!=1) break;
    }
  
  cout<<"Read "<<fTracks->GetNTracks()<<" tracks from file"<<endl;
  fclose(file);
}

void AliL3Compress::CompressFile(Char_t *infile,Char_t *outfile)
{
  
  BIT_FILE *output = OpenOutputBitFile(outfile);
  FILE *input = fopen(infile,"r");

  AliL3TrackModel track;
  AliL3ClusterModel cluster;
  Int_t temp;

  while(!feof(input))
    {
      if(fread(&track,sizeof(AliL3TrackModel),1,input)!=1) break;
      
      if(output->mask != 0x80) //Write the current byte to file.
	{
	  if(putc(output->rack,output->file )!=output->rack)
	    cerr<<"AliL3Compress::ComressFile : Error writing to bitfile"<<endl;
	  output->mask=0x80;
	}
      fwrite(&track,sizeof(AliL3TrackModel),1,output->file);
      for(Int_t i=0; i<track.fNClusters; i++)
	{
	  if(fread(&cluster,sizeof(AliL3ClusterModel),1,input)!=1) break;
	  Int_t flag = (Int_t)cluster.fPresent;
	  OutputBit(output,flag);
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
	  temp = (Int_t)cluster.fDCharge;
	  OutputBits(output,temp,10);
	  
	  //Short_t temp=(Short_t)cluster.fDTime;
	  // cout<<"flag "<<(int)flag<<" dtime "<<(int)cluster.fDTime<<" dpad "<<(int)cluster.fDPad<<" charge "<<cluster.fDCharge<<endl;
	}
      
    }
  
  fclose(input);
  CloseOutputBitFile(output);
}

void AliL3Compress::ExpandFile(Char_t *infile,Char_t *outfile)
{
  BIT_FILE *input = OpenInputBitFile(infile);
  FILE *output = fopen(outfile,"w");

  AliL3TrackModel trackmodel;
  AliL3ClusterModel *clusters=0;
  
  while(!feof(input->file))
    {
      
      if(fread(&trackmodel,sizeof(AliL3TrackModel),1,input->file)!=1) break;
      fwrite(&trackmodel,sizeof(AliL3TrackModel),1,output);
      input->mask=0x80;//make sure we read a new byte from file.
      clusters = new AliL3ClusterModel[(trackmodel.fNClusters)];
      for(Int_t i=0; i<trackmodel.fNClusters; i++)
	{
	  Int_t temp,sign;
	  temp = InputBit(input);
	  if(!temp) 
	    {
	      clusters[i].fPresent=kFALSE;
	      continue;
	    }
	  clusters[i].fPresent=kTRUE;
	  sign=InputBit(input);
	  temp = InputBits(input,8);
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDTime = temp;
	  sign=InputBit(input);
	  temp = InputBits(input,8);
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDPad = temp;
	  temp=InputBits(input,10);
	  clusters[i].fDCharge = temp;
	}
      fwrite(clusters,(trackmodel.fNClusters)*sizeof(AliL3ClusterModel),1,output);
      delete [] clusters;
    }
  fclose(output);
  CloseInputBitFile(input);
}
