// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ALICE HLT Group
//_____________________________________________________________
//
//  AliHLTCompress
//
// Class for compressing and uncompressing data.

#include "AliHLTStandardIncludes.h"

#include "bitio.h"
#include "AliHLTRootTypes.h"
#include "AliHLTModels.h"
#include "AliHLTDigitData.h"
#include "AliHLTLogging.h"
#include "AliHLTTrackArray.h"
#include "AliHLTModelTrack.h"
#include "AliHLTTransform.h"
#include "AliHLTMemHandler.h"
#include "AliHLTDataCompressorHelper.h"
#include "AliHLTDataCompressor.h"
#include "AliHLTSpacePointData.h"

#if 0
#ifdef use_root
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>
#endif
#ifdef use_aliroot
#include "AliHLTFileHandler.h"
#endif
#endif

#include "AliHLTCompress.h"

#if __GNUC__ >= 3
using namespace std;
#endif


ClassImp(AliHLTCompress)

AliHLTCompress::AliHLTCompress()
{
  // default constructor
  fTracks=0;
  fSlice =0;
  fPatch=0;
  fWriteShape=kFALSE;
  fEvent=-1;
}

AliHLTCompress::AliHLTCompress(Int_t slice,Int_t patch,Char_t *path,Bool_t writeshape,Int_t event)
{
  // constructor
  fEvent=event;
  fSlice=slice;
  fPatch=patch;
  fTracks=0;
  sprintf(fPath,"%s",path);
  fWriteShape=writeshape;
}

AliHLTCompress::~AliHLTCompress()
{
  // destructor
  if(fTracks)
    delete fTracks;
}

Bool_t AliHLTCompress::WriteFile(AliHLTTrackArray *tracks,Char_t *filename)
{
  // writes file
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
      cerr<<"AliHLTCompress::WriteFile : Error opening file "<<fname<<endl;
      return kFALSE;
    }
  Short_t ntracks = tracks->GetNTracks();
    
  Int_t count=0;
  AliHLTClusterModel *clusters=0;
  AliHLTTrackModel *model=0;
  for(Int_t i=0; i<ntracks; i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;

      //Do not save useless tracks or clusters:
      //if(track->GetNPresentClusters() == 0)
      //continue;
      
      track->FillModel();
      model = track->GetModel();
      //if(model->fNClusters==0) continue;
      clusters = track->GetClusters();
      if(fwrite(model,sizeof(AliHLTTrackModel),1,file)!=1) break;
      //if(fwrite(clusters,model->fNClusters*sizeof(AliHLTClusterModel),1,file)!=1) break;
      if(fwrite(clusters,AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTClusterModel),1,file)!=1) break;
      count++;
      
    }
  fclose(file);
  return kTRUE;
}

Bool_t AliHLTCompress::ReadFile(Char_t which,Char_t *filename)
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
	  cerr<<"AliHLTCompress::ReadFile() : Wrong option"<<endl;
	  return kFALSE;
	}
    }

  FILE *file = fopen(fname,"r");
  if(!file)
    {
      cerr<<"AliHLTCompress::ReadFile : Cannot open file "<<fname<<endl;
      return kFALSE;
    }

  if(fTracks)
    delete fTracks;
  fTracks = new AliHLTTrackArray("AliHLTModelTrack");
  
  while(!feof(file))
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->NextTrack();
      track->Init(fSlice,fPatch);
      AliHLTTrackModel *model = track->GetModel();
      AliHLTClusterModel *clusters = track->GetClusters();
      if(fread(model,sizeof(AliHLTTrackModel),1,file)!=1) break;
      //if(fread(clusters,model->fNClusters*sizeof(AliHLTClusterModel),1,file)!=1) break;
      if(fread(clusters,AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTClusterModel),1,file)!=1) break;
      track->FillTrack();
    }

  fTracks->RemoveLast();
  fclose(file);
  return kTRUE;
}

Bool_t AliHLTCompress::CompressFile()
{
  // compresses file
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
      cerr<<"AliHLTCompress::CompressFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliHLTTrackModel track;
  AliHLTClusterModel cluster;
  Int_t temp;
  Int_t power;
  
  Int_t timeo,pado,chargeo,padshapeo,timeshapeo;
  timeo=pado=chargeo=padshapeo=timeshapeo=0;
  while(!feof(input))
    {
      if(fread(&track,sizeof(AliHLTTrackModel),1,input)!=1) break;
      
      if(output->mask != 0x80) //Write the current byte to file.
	{
	  //cerr<<"\nAliHLTCompress::CompressFile() : Writing overhead bits!!!"<<endl;
	  if(putc(output->rack,output->file )!=output->rack)
	    cerr<<"AliHLTCompress::ComressFile : Error writing to bitfile"<<endl;
	  output->mask=0x80;
	  output->rack=0;
	}
      
      //Write track parameters:
      fwrite(&track,sizeof(AliHLTTrackModel),1,output->file);
      
      Int_t origslice=-1,slice,clustercount=0;
      for(Int_t i=0; i<AliHLTTransform::GetNRows(fPatch); i++)
	{
	  if(fread(&cluster,sizeof(AliHLTClusterModel),1,input)!=1) break;
	  
	  //Write empty flag:
	  temp = (Int_t)cluster.fPresent;
	  OutputBit(output,temp);
	  if(!temp) continue;
	  
	  if(cluster.fSlice<0 || cluster.fSlice>35)
	    {
	      cerr<<"AliHLTDataCompress::CompressFile : Fucked up slice number :"<<cluster.fSlice<<endl;
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
		OutputBit(output,0);  //No change of slice
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
	  power = 1<<(AliHLTDataCompressorHelper::GetNTimeBits()-1);
	  if(abs(temp)>=power)
	    {
	      //cout<<abs(temp)<<" "<<power<<endl;
	      timeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTDataCompressorHelper::GetNTimeBits()-1));
	  
	  //Write pad information:
	  temp = (Int_t)rint(cluster.fDPad);
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<(AliHLTDataCompressorHelper::GetNPadBits()-1);
	  if(abs(temp)>=power)
	    {
	      pado++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTDataCompressorHelper::GetNPadBits()-1));
	  
	  //Write charge information:
	  temp = (Int_t)cluster.fDCharge;
	  power = 1<<(AliHLTDataCompressorHelper::GetNChargeBits());
	  if(abs(temp)>=power)
	    {
	      chargeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),(AliHLTDataCompressorHelper::GetNChargeBits()));
	  
	  if(fWriteShape)
	    {
	      //Write shape information:
	      temp = (Int_t)rint(cluster.fDSigmaY);
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliHLTDataCompressorHelper::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  padshapeo++;
		  temp = power - 1;
		}
	      OutputBits(output,abs(temp),(AliHLTDataCompressorHelper::GetNShapeBits()-1));
	      
	      temp = (Int_t)rint(cluster.fDSigmaZ);
	      if(temp<0)
		OutputBit(output,0);
	      else
		OutputBit(output,1);
	      power = 1<<(AliHLTDataCompressorHelper::GetNShapeBits()-1);
	      if(abs(temp) >= power)
		{
		  timeshapeo++;
		  temp=power - 1;
		}
	      OutputBits(output,abs(temp),(AliHLTDataCompressorHelper::GetNShapeBits()-1));
	    }
	  
	  clustercount++;
	}
    }
  
  fclose(input);
  CloseOutputBitFile(output);
  if(pado || timeo || chargeo || padshapeo || timeshapeo)
    {
      cout<<endl<<"Saturations: "<<endl
	  <<"Pad "<<pado<<endl
	  <<"Time "<<timeo<<endl
	  <<"Charge "<<chargeo<<endl
	  <<"Padshape "<<padshapeo<<endl
	  <<"Timeshape "<<timeshapeo<<endl<<endl;
    }
  return kTRUE;
}

Bool_t AliHLTCompress::ExpandFile()
{
  // expands file
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
      cerr<<"AliHLTCompress::ExpandFile() : Error opening file: "<<fname<<endl;
      return kFALSE;
    }

  AliHLTTrackModel trackmodel;
  AliHLTClusterModel *clusters=0;
  Int_t count=0;
  
  clusters = new AliHLTClusterModel[(AliHLTTransform::GetNRows(fPatch))];
  while(!feof(input->file))
    {
      input->mask=0x80;//make sure we read a new byte from file.
      
      //Read and write track:
      if(fread(&trackmodel,sizeof(AliHLTTrackModel),1,input->file)!=1) break;
      fwrite(&trackmodel,sizeof(AliHLTTrackModel),1,output);

      memset(clusters,0,AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTClusterModel));
      Int_t origslice=-1,clustercount=0;
      for(Int_t i=0; i<AliHLTTransform::GetNRows(fPatch); i++)
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
	  temp = InputBits(input,(AliHLTDataCompressorHelper::GetNTimeBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDTime = temp;
	  
	  //Read pad information:
	  sign=InputBit(input);
	  temp = InputBits(input,(AliHLTDataCompressorHelper::GetNPadBits()-1));
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDPad = temp;
	  
	  //Read charge information:
	  temp=InputBits(input,(AliHLTDataCompressorHelper::GetNChargeBits()));
	  clusters[i].fDCharge = temp;
	  
	  if(fWriteShape)
	    {
	      //Read shape information:
	      sign = InputBit(input);
	      temp = InputBits(input,(AliHLTDataCompressorHelper::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaY = temp;
	      
	      sign = InputBit(input);
	      temp = InputBits(input,(AliHLTDataCompressorHelper::GetNShapeBits()-1));
	      if(!sign)
		temp*=-1;
	      clusters[i].fDSigmaZ = temp;
	    }
	  clustercount++;
	}
      count++;
      //fwrite(clusters,(trackmodel.fNClusters)*sizeof(AliHLTClusterModel),1,output);
      fwrite(clusters,AliHLTTransform::GetNRows(fPatch)*sizeof(AliHLTClusterModel),1,output);
    }
  
  delete [] clusters;
  fclose(output);
  CloseInputBitFile(input);
  return kTRUE;
}

void AliHLTCompress::CompressRemaining(AliHLTSpacePointData *clusters[36][6],UInt_t nclusters[36][6])
{
  //Write the remaining clusters in a compressed format.

  Char_t filename[1024];
  Int_t nrows = AliHLTTransform::GetNRows();
  Int_t *npoints = new Int_t[nrows];
  for(Int_t slice=0; slice<=35; slice++)
    {
      for(Int_t patch=0; patch < 1; patch++)
	{
	  sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  BIT_FILE *output = OpenOutputBitFile(filename);
	  if(!output)
	    {
	      cerr<<"AliHLTCompress::CompressRemaining : Cannot open file "<<filename<<endl;
	      exit(5);
	    }
	  AliHLTSpacePointData *cl = clusters[slice][patch];
	  memset(npoints,0,nrows*sizeof(Int_t));
	  
	  UInt_t i;
	  for(i=0; i<nclusters[slice][patch]; i++)
	    {
	      if(cl[i].fCharge == 0) continue; //has been used
	      npoints[cl[i].fPadRow]++;
	    }
	  Int_t rowspresent=0;
	  for(Int_t j=0; j<nrows; j++)
	    {
	      if(!npoints[j]) continue;
	      rowspresent++;
	    }
	  
	  //Write number of padrows with clusters
	  OutputBits(output,rowspresent,8);
	  
	  Int_t lastPadrow=-1;
	  for(i=0; i<nclusters[slice][patch]; i++)
	    {
	      if(cl[i].fCharge == 0) continue; //has been used
	      Int_t padrow = cl[i].fPadRow;
	      if(padrow != lastPadrow)
		{
		  OutputBits(output,padrow,8);//Write padrow #
		  if(npoints[padrow] >= 1<<10)
		    {
		      cerr<<"AliHLTCompress::CompressRemaining : Too many remaining clusters "<<npoints[padrow]<<endl;
		      exit(5);
		    }
		  OutputBits(output,npoints[padrow],10);//Write number of clusters on this padrow
		  lastPadrow = padrow;
		}
	      
	      Float_t xyz[3] = {cl[i].fX,cl[i].fY,cl[i].fZ};
	      Int_t sector,row,buff;
	      AliHLTTransform::Slice2Sector(slice,padrow,sector,row);
	      AliHLTTransform::Global2Raw(xyz,sector,row);
	      
	      Float_t padw = sqrt(cl[i].fSigmaY2) / AliHLTTransform::GetPadPitchWidth(AliHLTTransform::GetPatch(padrow));
	      Float_t timew = sqrt(cl[i].fSigmaZ2) / AliHLTTransform::GetZWidth();
	      
	      //Check for saturation in the widths.
	      //Basically only store a certain number of decimals here, and cut the widths which is higher:
	      if(padw >= (1<<AliHLTDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTDataCompressorHelper::GetPadPrecisionFactor())
		padw = (1<<AliHLTDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTDataCompressorHelper::GetPadPrecisionFactor() - 1/AliHLTDataCompressorHelper::GetPadPrecisionFactor();
	      if(timew >= (1<<AliHLTDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTDataCompressorHelper::GetTimePrecisionFactor())
		timew = (1<<AliHLTDataCompressorHelper::GetNShapeBitsRemaining()) / AliHLTDataCompressorHelper::GetTimePrecisionFactor() - 1/AliHLTDataCompressorHelper::GetTimePrecisionFactor();;
	      
	      //Write pad
	      buff = (Int_t)rint(xyz[1]*AliHLTDataCompressorHelper::GetPadPrecisionFactor());
	      if(buff<0)
		{
		  cerr<<"AliHLTCompress:CompressRemaining : Wrong pad value "<<buff<<endl;
		  exit(5);
		}
	      OutputBits(output,buff,AliHLTDataCompressorHelper::GetNPadBitsRemaining());
	      
	      //Write time
	      buff = (Int_t)rint(xyz[2]*AliHLTDataCompressorHelper::GetTimePrecisionFactor());
	      if(buff<0)
		{
		  cerr<<"AliHLTCompress:CompressRemaining : Wrong time value "<<buff<<endl;
		  exit(5);
		}
	      OutputBits(output,buff,AliHLTDataCompressorHelper::GetNTimeBitsRemaining());
	      
	      //Write widths
	      buff = (Int_t)rint(padw*AliHLTDataCompressorHelper::GetPadPrecisionFactor());
	      OutputBits(output,buff,AliHLTDataCompressorHelper::GetNShapeBitsRemaining());
	      buff = (Int_t)rint(timew*AliHLTDataCompressorHelper::GetTimePrecisionFactor());
	      OutputBits(output,buff,AliHLTDataCompressorHelper::GetNShapeBitsRemaining());
	      
	      //Write charge 
	      buff = cl[i].fCharge;
	      if(buff >= 1<<(AliHLTDataCompressorHelper::GetNChargeBits()))
		buff = (1<<(AliHLTDataCompressorHelper::GetNChargeBits()))-1;
	      OutputBits(output,buff,AliHLTDataCompressorHelper::GetNChargeBits());
	    }
	  
	  CloseOutputBitFile(output);
	}
      
    }
  delete [] npoints;
}

void AliHLTCompress::ExpandRemaining(TempCluster **clusters,Int_t *ncl, Int_t /*maxpoints*/)
{
  //Expand the remaining clusters stored using function CompressRemaining
  
  Char_t filename[1024];
  Int_t buff;
  for(Int_t slice=0; slice<=35; slice++)
    {
      for(Int_t p=0; p<1; p++)
	{
	  sprintf(filename,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,slice,-1);
	  BIT_FILE *input = OpenInputBitFile(filename);
	  
	  //Read number of padrows 
	  buff = InputBits(input,8);
	  Int_t nrows = buff;
	  
	  for(Int_t i=0; i<nrows; i++)
	    {
	      //Read padrow
	      buff = InputBits(input,8);
	      Int_t padrow = buff;
	      
	      //Read nclusters;
	      buff = InputBits(input,10);
	      Int_t npoints = buff;
	      
	      for(Int_t i=0; i<npoints; i++)
		{
		  clusters[slice][ncl[slice]].fPadrow = padrow;

		  //Read pad
		  buff = InputBits(input,AliHLTDataCompressorHelper::GetNPadBitsRemaining());
		  clusters[slice][ncl[slice]].fPad = (Float_t)buff/AliHLTDataCompressorHelper::GetPadPrecisionFactor();
		  
		  //Read time
		  buff = InputBits(input,AliHLTDataCompressorHelper::GetNTimeBitsRemaining());
		  clusters[slice][ncl[slice]].fTime = (Float_t)buff/AliHLTDataCompressorHelper::GetTimePrecisionFactor();
		  
		  //Read widths 
		  buff = InputBits(input,AliHLTDataCompressorHelper::GetNShapeBitsRemaining());
		  clusters[slice][ncl[slice]].fSigmaY2 = pow((Float_t)buff/AliHLTDataCompressorHelper::GetPadPrecisionFactor(),2);
		  buff = InputBits(input,AliHLTDataCompressorHelper::GetNShapeBitsRemaining());
		  clusters[slice][ncl[slice]].fSigmaZ2 = pow((Float_t)buff/AliHLTDataCompressorHelper::GetPadPrecisionFactor(),2);
		  
		  //Read charge
		  buff = InputBits(input,AliHLTDataCompressorHelper::GetNChargeBits());
		  clusters[slice][ncl[slice]].fCharge = buff;
		  
		  ncl[slice]++;
		}
	      
	    }
	  CloseInputBitFile(input);
	}
    }
}

void AliHLTCompress::PrintCompRatio(ofstream *outfile)
{
  // prints the compression ratio
  AliHLTMemHandler *mem = new AliHLTMemHandler();
  Char_t fname[1024];
  UInt_t remainSize=0,digitSize=0;
  for(Int_t i=0; i<36; i++)
    {
      if(fEvent<0)
	sprintf(fname,"%s/comp/remains_%d_%d.raw",fPath,i,-1);
      else
	sprintf(fname,"%s/comp/remains_%d_%d_%d.raw",fPath,fEvent,i,-1);
      mem->SetBinaryInput(fname);
      remainSize += mem->GetFileSize();
      mem->CloseBinaryInput();

      sprintf(fname,"%s/binaries/digits_c8_%d_%d_%d.raw",fPath,fEvent,i,-1);
      mem->SetBinaryInput(fname);
      digitSize += mem->GetFileSize();
      mem->CloseBinaryInput();
    }
  
  
  if(fEvent<0)
    sprintf(fname,"%s/comp/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  else
    sprintf(fname,"%s/comp/tracks_c_%d_%d_%d.raw",fPath,fEvent,fSlice,fPatch);

  mem->SetBinaryInput(fname);
  UInt_t compressSize = mem->GetFileSize();
  mem->CloseBinaryInput();
  
  if(digitSize==0)
    {
      cerr<<"AliHLTCompress::PrintCompRatio : Zero digit size, not able to obtain comp. ratios!"<<endl;
      return;
    }
  
  Float_t compratio = (Float_t)(compressSize + remainSize)/(Float_t)digitSize;
  Float_t entropy[3];
  Int_t trackSize = GetEntropy(entropy[0],entropy[1],entropy[2])*sizeof(AliHLTTrackModel);
  if(outfile)
    {
      ofstream &out = *outfile;
      out<<compressSize<<' '<<remainSize<<' '<<digitSize<<' '<<trackSize<<' '<<entropy[0]<<' '<<entropy[1]<<endl;
    }
  
  cout<<"=========================================="<<endl;
  cout<<"Original digits size : "<<digitSize/1000<<" kByte ( 100 % )"<<endl;
  cout<<"Compressed file size : "<<compressSize/1000<<" kByte ( "<<(Float_t)compressSize*100/(Float_t)digitSize<<" % )"<<endl;
  cout<<"Remaining file size  : "<<remainSize/1000<<" kByte ( "<<(Float_t)remainSize*100/(Float_t)digitSize<<" % )"<<endl;
  cout<<"Relative track size  : "<<trackSize/1000<<" kByte ( "<<(Float_t)trackSize*100/(Float_t)digitSize<<" % )"<<endl;
  cout<<"Relative cluster size: "<<(compressSize-trackSize)/1000<<" kByte ( "<<(Float_t)(compressSize-trackSize)*100/(Float_t)digitSize<<" % )"<<endl;
  cout<<"---------------------- "<<endl;
  cout<<"Compression ratio    : "<<compratio*100<<" %"<<endl;
  cout<<"=========================================="<<endl;
  cout<<"Entropy of residuals : "<<entropy[0]<<" "<<entropy[1]<<endl;
}

Int_t AliHLTCompress::GetEntropy(Float_t &padEntropy,Float_t &timeEntropy,Float_t &chargeEntropy)
{
  //Calculate the entropy of the quantized residuals in both directions
  
  if(!ReadFile('m'))
    return 0;
  const Int_t knmax=100000;
  Float_t pads[knmax];
  Float_t times[knmax];
  Float_t charge[knmax];
  memset(&pads[0],0,knmax*sizeof(Float_t));
  memset(&times[0],0,knmax*sizeof(Float_t));
  memset(&charge[0],0,knmax*sizeof(Float_t));
  Float_t counter=0;

  for(Int_t i=0; i<fTracks->GetNTracks(); i++)
    {
      AliHLTModelTrack *track = (AliHLTModelTrack*)fTracks->GetCheckedTrack(i);
      if(!track) continue;
      for(Int_t padrow=0; padrow<AliHLTTransform::GetNRows(); padrow++)
	{
	  if(!track->IsPresent(padrow)) continue;
	  Int_t dpad = abs((Int_t)rint(track->GetClusterModel(padrow)->fDPad));
	  Int_t dtime = abs((Int_t)rint(track->GetClusterModel(padrow)->fDTime));
	  Int_t dcharge = (Int_t)track->GetClusterModel(padrow)->fDCharge;
	  if(dpad >= knmax || dtime >= knmax || dcharge >= knmax)
	    {
	      cerr<<"AliHLTCompress::GetEntropy : Quantization out of range: "<<dpad<<" "<<dtime<<" "<<dcharge<<endl;
	      break;
	    }
	  pads[dpad]++;
	  times[dtime]++;
	  charge[dcharge]++;
	  counter++;
	}
    }
  padEntropy=timeEntropy=chargeEntropy=0;
  for(Int_t i=0; i<knmax; i++)
    {
      if(pads[i]>0)
	padEntropy += (pads[i]/counter)*(log(pads[i]/counter)/log(2.0));
      if(times[i]>0)
	timeEntropy += (times[i]/counter)*(log(times[i]/counter)/log(2.0));
      if(charge[i]>0)
	chargeEntropy += (charge[i]/counter)*(log(charge[i]/counter)/log(2.0));
    }
  
  padEntropy*=-1;
  timeEntropy*=-1;
  chargeEntropy*=-1;
  return fTracks->GetNTracks();
}
