//$Id$

// Author: Anders Vestbo <mailto:vestbo$fi.uib.no>
//*-- Copyright &copy ASV

#include <stdio.h>
#include <stream.h>
#include <stdlib.h>
#include <TH1.h>
#include <TH2.h>
#include <TRandom.h>

#include "AliL3Compress.h"
#include "AliL3TrackArray.h"
#include "AliL3ModelTrack.h"
#include "AliL3Transform.h"
#include "AliL3MemHandler.h"
#include "AliL3FileHandler.h"
#include "bitio.h"

//_____________________________________________________________
//
//  AliL3Compress
//
// Class for compressing and uncompressing data.

ClassImp(AliL3Compress)

AliL3Compress::AliL3Compress()
{
  fTracks=0;
  SetBitNumbers(0,0,0,0);
  fSlice =0;
  fPatch=0;
  fDigits=0;
  fDPt=0;
}

AliL3Compress::AliL3Compress(Int_t slice,Int_t patch,Char_t *path)
{
  fSlice=slice;
  fPatch=patch;
  SetBitNumbers(0,0,0,0);
  fTracks=0;
  fDigits=0;
  fDPt=0;
  sprintf(fPath,"%s",path);
}

AliL3Compress::~AliL3Compress()
{
  if(fTracks)
    delete fTracks;
  if(fDigits)
    delete [] fDigits;
  if(fDPt)
    delete [] fDPt;
}

void AliL3Compress::SetBitNumbers(Int_t pad,Int_t time,Int_t charge,Int_t shape)
{
  fNumPadBits=pad;
  fNumTimeBits=time;
  fNumChargeBits=charge;
  fNumShapeBits=shape;
}

void AliL3Compress::WriteFile(AliL3TrackArray *tracks)
{
  Char_t fname[100];
  sprintf(fname,"%s/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  FILE *file = fopen(fname,"w");
  Short_t ntracks = tracks->GetNTracks();
  //cout<<"Writing "<<ntracks<<" tracks to file"<<endl;
    
  Int_t count=0;
  AliL3ClusterModel *clusters=0;
  AliL3TrackModel *model=0;
  for(Int_t i=0; i<ntracks; i++)
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)tracks->GetCheckedTrack(i);
      if(!track) continue;
            
      track->FillModel();
      model = track->GetModel();
      if(model->fNClusters==0) continue;
      clusters = track->GetClusters();
      //cout<<"Writing "<<(int)model->fNClusters<<" clusters"<<endl;
      if(fwrite(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      //cout<<"Writing "<<(int)model->fNClusters<<" clusters to file"<<endl;
      if(fwrite(clusters,model->fNClusters*sizeof(AliL3ClusterModel),1,file)!=1) break;
      //track->Print();
      count++;
      
    }
  cout<<"Wrote "<<count<<" tracks "<<endl;
  fclose(file);
}

void AliL3Compress::ReadFile(Char_t which)
{
  //Read the trackfile.

  Char_t fname[100];
  if(which == 'm')
    sprintf(fname,"%s/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  else if(which == 'u')
    sprintf(fname,"%s/tracks_u_%d_%d.raw",fPath,fSlice,fPatch);
  else
    {
      cerr<<"AliL3Compress::ReadFile() : Wrong option"<<endl;
      return;
    }

  FILE *file = fopen(fname,"r");
  if(!file)
    {
      cerr<<"Cannot open file "<<fname<<endl;
      return;
    }

  if(fTracks)
    delete fTracks;
  fTracks = new AliL3TrackArray("AliL3ModelTrack");
  
  cout<<"Reading file "<<fname<<endl;
  while(!feof(file))
    {
      AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->NextTrack();
      track->Init(fSlice,fPatch);
      AliL3TrackModel *model = track->GetModel();
      AliL3ClusterModel *clusters = track->GetClusters();
      //cout<<"Reading model "<<(int)model<<endl;
      if(fread(model,sizeof(AliL3TrackModel),1,file)!=1) break;
      //cout<<"Reading clusters "<<(int)clusters<<endl;
      if(fread(clusters,(model->fNClusters)*sizeof(AliL3ClusterModel),1,file)!=1) break;
      //cout<<"Filling track"<<endl;
      track->FillTrack();
      //track->Print();
    }

  fTracks->RemoveLast();
  cout<<"Read "<<fTracks->GetNTracks()<<" tracks from file"<<endl;
  fclose(file);
}

void AliL3Compress::CompressFile()
{
  if(fNumTimeBits==0)
    {
      cerr<<"AliL3Compress::CompressFile() : Bitnumbers not set"<<endl;
      return;
    }
  
  Char_t fname[100];
  sprintf(fname,"%s/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  BIT_FILE *output = OpenOutputBitFile(fname);
  
  sprintf(fname,"%s/tracks_m_%d_%d.raw",fPath,fSlice,fPatch);
  FILE *input = fopen(fname,"r");
  if(!input)
    {
      cerr<<"AliL3Compress::CompressFile() : Error opening file: "<<fname<<endl;
      return;
    }

  AliL3TrackModel track;
  AliL3ClusterModel cluster;
  Int_t temp;
  Short_t power;
  
  Int_t timeo,pado,chargeo,shapeo;
  timeo=pado=chargeo=shapeo=0;
  while(!feof(input))
    {
      if(fread(&track,sizeof(AliL3TrackModel),1,input)!=1) break;
      
      if(output->mask != 0x80) //Write the current byte to file.
	{
	  if(putc(output->rack,output->file )!=output->rack)
	    cerr<<"AliL3Compress::ComressFile : Error writing to bitfile"<<endl;
	  output->mask=0x80;
	  output->rack=0;
	}
      
      //Write track parameters:
      fwrite(&track,sizeof(AliL3TrackModel),1,output->file);
      for(Int_t i=0; i<track.fNClusters; i++)
	{
	  if(fread(&cluster,sizeof(AliL3ClusterModel),1,input)!=1) break;
	  
	  //Write empty flag:
	  temp = (Int_t)cluster.fPresent;
	  OutputBit(output,temp);
	  if(!temp) continue;
	  
	  //Write time information:
	  temp = (Int_t)cluster.fDTime;
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<fNumTimeBits;
	  if(abs(temp)>=power)
	    {
	      timeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),fNumTimeBits);
	  
	  //Write pad information:
	  temp = (Int_t)cluster.fDPad;
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<fNumPadBits;
	  if(abs(temp)>=power)
	    {
	      pado++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),fNumPadBits);
	  
	  //Write charge information:
	  temp = (Int_t)cluster.fDCharge;
	  if(temp<0)
	    OutputBit(output,0);
	  else
	    OutputBit(output,1);
	  power = 1<<fNumChargeBits;
	  if(abs(temp)>=power)
	    {
	      chargeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),fNumChargeBits);
	  
	  //Write shape information:
	  temp = (Int_t)cluster.fDSigmaY2;
	  power = 1<<fNumShapeBits;
	  if(abs(temp) >= power)
	    {
	      shapeo++;
	      temp = power - 1;
	    }
	  OutputBits(output,abs(temp),fNumShapeBits);
	  
	  temp = (Int_t)cluster.fDSigmaZ2;
	  if(abs(temp) >= power)
	    {
	      shapeo++;
	      temp=power - 1;
	    }
	  OutputBits(output,abs(temp),fNumShapeBits);
	}
    }
  
  fclose(input);
  CloseOutputBitFile(output);

  cout<<endl<<"There was following number of overflows: "<<endl
      <<"Pad "<<pado<<endl
      <<"Time "<<timeo<<endl
      <<"Charge "<<chargeo<<endl
      <<"Shape "<<shapeo<<endl;
}

void AliL3Compress::ExpandFile()
{
  if(fNumTimeBits==0)
    {
      cerr<<"AliL3Compress::ExpandFile() : Bitnumbers not set"<<endl;
      return;
    }
  
  Char_t fname[100];
  sprintf(fname,"%s/tracks_c_%d_%d.raw",fPath,fSlice,fPatch);
  BIT_FILE *input = OpenInputBitFile(fname);
  
  sprintf(fname,"%s/tracks_u_%d_%d.raw",fPath,fSlice,fPatch);
  FILE *output = fopen(fname,"w");
  if(!output)
    {
      cerr<<"AliL3Compress::ExpandFile() : Error opening file: "<<fname<<endl;
      return;
    }

  AliL3TrackModel trackmodel;
  AliL3ClusterModel *clusters=0;
  Int_t count=0;
  
  clusters = new AliL3ClusterModel[(NumRows[fPatch])];
  while(!feof(input->file))
    {
      input->mask=0x80;//make sure we read a new byte from file.
      
      //Read and write track:
      if(fread(&trackmodel,sizeof(AliL3TrackModel),1,input->file)!=1) break;
      fwrite(&trackmodel,sizeof(AliL3TrackModel),1,output);
      
      for(Int_t i=0; i<NumRows[fPatch]; i++)
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
	  
	  //Read time information:
	  sign=InputBit(input);
	  temp = InputBits(input,fNumTimeBits);
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDTime = temp;
	  
	  //Read pad information:
	  sign=InputBit(input);
	  temp = InputBits(input,fNumPadBits);
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDPad = temp;
	  
	  //Read charge information:
	  sign = InputBit(input);
	  temp=InputBits(input,fNumChargeBits);
	  if(!sign)
	    temp*=-1;
	  clusters[i].fDCharge = temp;
	  
	  //Read shape information:
	  temp = InputBits(input,fNumShapeBits);
	  clusters[i].fDSigmaY2 = temp;
	  
	  temp = InputBits(input,fNumShapeBits);
	  clusters[i].fDSigmaZ2 = temp;
	}
      

      count++;
      fwrite(clusters,(trackmodel.fNClusters)*sizeof(AliL3ClusterModel),1,output);
      
    }
  
  delete [] clusters;
  fclose(output);
  CloseInputBitFile(input);
}

void AliL3Compress::CreateDigitArray(Int_t maxnumber)
{
  fNUsed=0;
  fNDigits = 0;
  fMaxDigits=maxnumber;
  if(fDigits) delete [] fDigits;
  fDigits = new AliL3RandomDigitData[maxnumber];
  if(fDPt) delete [] fDPt;
  fDPt = new AliL3RandomDigitData*[maxnumber];
}

void AliL3Compress::RestoreData()
{
  
  //Read the uncompressed file:
  ReadFile('u');
  
  CreateDigitArray(100000);
  
  Float_t pad,time,sigmaY2,sigmaZ2;
  Int_t charge;
  for(Int_t j=NRows[fPatch][0]; j<=NRows[fPatch][1]; j++)
    {
      for(Int_t i=0; i<fTracks->GetNTracks(); i++)
	{
	  AliL3ModelTrack *track = (AliL3ModelTrack*)fTracks->GetCheckedTrack(i);
	  if(!track) continue;
	  if(!track->GetPad(j,pad) || 
	     !track->GetTime(j,time) || 
	     !track->GetClusterCharge(j,charge) ||
	     !track->GetXYWidth(j,sigmaY2) || 
	     !track->GetZWidth(j,sigmaZ2))
	    continue;
	  
	  CreateDigits(j,pad,time,charge,sigmaY2,sigmaZ2);
	}
    }
  
  QSort(fDPt,0,fNDigits);
}

void AliL3Compress::PrintDigits()
{
  Int_t pad,time,charge,row;
  for(Int_t i=0; i<fNDigits; i++)
    {
      row = fDPt[i]->fRow;
      pad = fDPt[i]->fPad;
      time = fDPt[i]->fTime;
      charge = fDPt[i]->fCharge;
      if(i>0 && row != fDPt[i-1]->fRow)
	cout<<"---Padrow "<<row<<"---"<<endl;
      cout<<"Pad "<<pad<<" time "<<time<<" charge "<<charge<<endl;
    }
}

void AliL3Compress::WriteRestoredData()
{
  Char_t fname[100];
  
  //Get the remaining raw data array:
  AliL3MemHandler *mem = new AliL3MemHandler();
  sprintf(fname,"%s/remains_%d_%d.raw",fPath,fSlice,fPatch);
  mem->SetBinaryInput(fname);
  UInt_t numdigits;
  AliL3DigitRowData *origRow = mem->CompBinary2Memory(numdigits);
  mem->CloseBinaryInput();
  
  //Allocate memory for the merged data:
  UInt_t size = mem->GetAllocatedSize() + fNDigits*sizeof(AliL3DigitData);
  cout<<"Allocating "<<size<<" bytes for merged data array "<<endl;
  Byte_t *data = new Byte_t[size];
  memset(data,0,size);
  AliL3DigitRowData *tempRow = (AliL3DigitRowData*)data;

  Int_t ndigits,action,charge;
  UShort_t pad,time;
      
  UInt_t digit_counter;
  Int_t row_counter=0;
  for(Int_t i=NRows[fPatch][0]; i<=NRows[fPatch][1]; i++)
    {
      tempRow->fRow = i;
      ndigits=0;
      AliL3DigitData *origDig = origRow->fDigitData;
      AliL3DigitData *tempDig = tempRow->fDigitData;
      if((Int_t)origRow->fRow != i)
	cerr<<"AliL3Compress::WriteRestoredData() : Mismatching row numbering "<<(Int_t)origRow->fRow<<" "<<i<<endl;

      //cout<<"Writing row "<<i<<" with "<<(Int_t)origRow->fNDigit<<" old digits"<<endl;
      digit_counter=0;
      
      while(1)
	{
	  while(digit_counter < origRow->fNDigit)
	    {
	      pad = origDig[digit_counter].fPad;
	      time = origDig[digit_counter].fTime;
	      charge = origDig[digit_counter].fCharge;
	      digit_counter++;
	      while((action=ComparePoints(i,pad,time)) == 1)
		{
		  tempDig[ndigits].fPad = fDPt[fNUsed]->fPad;
		  tempDig[ndigits].fTime = fDPt[fNUsed]->fTime;
		  tempDig[ndigits].fCharge = fDPt[fNUsed]->fCharge;
		  ndigits++;
		  fNUsed++;

		}
	      if(action == 0)
		{
		  tempDig[ndigits].fPad = pad;
		  tempDig[ndigits].fTime = time;
		  tempDig[ndigits].fCharge = charge;
		  ndigits++;
		}
	    }
	  
	  if(fNUsed >= fNDigits) 
	    break;
	  if(fDPt[fNUsed]->fRow != i) //we are on a new row
	    break;
	  tempDig[ndigits].fPad = fDPt[fNUsed]->fPad;
	  tempDig[ndigits].fTime = fDPt[fNUsed]->fTime;
	  tempDig[ndigits].fCharge = fDPt[fNUsed]->fCharge;
	  ndigits++;
	  fNUsed++;
	}
      //cout<<"Writing "<<ndigits<<" digits on row "<<i<<endl;
      if(ndigits < 4)
	{
	  row_counter++;
	  cout<<"Few digits on row "<<i<<endl;
	}
      tempRow->fNDigit = ndigits;
      Int_t size = sizeof(AliL3DigitData)*tempRow->fNDigit + sizeof(AliL3DigitRowData);
      Byte_t *byte_pt = (Byte_t*)tempRow;
      byte_pt += size;
      tempRow = (AliL3DigitRowData*)byte_pt;
      mem->UpdateRowPointer(origRow);
    }
  
  if(row_counter != NumRows[fPatch])
    cerr<<"AliL3Compress::WriteRestoredData() : Written rows: "<<row_counter<<" total rows "<<NumRows[fPatch]<<endl;
  
  mem->Free();  
  sprintf(fname,"%s/restored_%d_%d.raw",fPath,fSlice,fPatch);
  mem->SetBinaryOutput(fname);
  mem->Memory2CompBinary((UInt_t)NumRows[fPatch],(AliL3DigitRowData*)data);
  mem->CloseBinaryOutput();
  
  delete [] data;
  delete mem;
  
}

void AliL3Compress::CreateDigits(Int_t row,Float_t pad,Float_t time,Int_t charge,Float_t sigmaY2,Float_t sigmaZ2)
{
  //Create raw data out of the cluster.
  
  AliL3Transform *tr = new AliL3Transform();
  TRandom *random = new TRandom();
  
  Int_t entries=10000;
  TH1F *hist1 = new TH1F("hist1","",tr->GetNPads(row),0,tr->GetNPads(row)-1);
  TH1F *hist2 = new TH1F("hist2","",tr->GetNTimeBins(),0,tr->GetNTimeBins()-1);
  TH2F *hist3 = new TH2F("hist3","",tr->GetNPads(row),0,tr->GetNPads(row)-1,tr->GetNTimeBins(),0,tr->GetNTimeBins()-1);
  
  //Convert back the sigmas:
  Float_t padw,timew;
  if(fPatch < 3)
    padw = tr->GetPadPitchWidthLow();
  else
    padw = tr->GetPadPitchWidthUp();
  timew = tr->GetZWidth();

  if(fPatch < 3)
    sigmaY2 = sigmaY2/2.07;
  sigmaY2 = sigmaY2/0.108;
  sigmaY2 = sigmaY2/(padw*padw);
  sigmaY2 = sigmaY2 - 1./12;
  
  if(fPatch < 3)
    sigmaZ2 = sigmaZ2/1.77;
  sigmaZ2 = sigmaZ2/0.169;
  sigmaZ2 = sigmaZ2/(timew*timew);
  sigmaZ2 = sigmaZ2 - 1./12;
  
  if(sigmaY2 <= 0 || sigmaZ2 <= 0)
    {
      cerr<<"AliL3Compress::CreateDigits() : Wrong sigmas : "<<sigmaY2<<" "<<sigmaZ2<<endl;
      return;
    }
  
  //Create the distributions in pad and time:
  for(Int_t i=0; i<entries; i++)
    {
      hist1->Fill(random->Gaus(pad,sqrt(sigmaY2)));
      hist2->Fill(random->Gaus(time,sqrt(sigmaZ2)));
    }
    
  //Create the cluster:
  Int_t bin1,bin2;
  Double_t content1,content2,dpad,dtime;
  for(Int_t i=0; i<hist1->GetEntries(); i++)
    {
      bin1 = hist1->GetBin(i);
      content1 = hist1->GetBinContent(bin1);
      if((Int_t)content1==0) continue;
      content1 = charge*content1/entries;
      dpad = hist1->GetBinCenter(bin1);
      for(Int_t j=0; j<hist2->GetEntries(); j++)
	{
	  bin2 = hist2->GetBin(j);
	  content2 = hist2->GetBinContent(bin2);
	  if((Int_t)content2==0) continue;
	  content2 = content1*content2/entries;
	  dtime = hist2->GetBinCenter(bin2);
	  hist3->Fill(dpad,dtime,content2);
	}
    }
  
  //Fill it into the digit array:
  for(Int_t i=0; i<hist3->GetNbinsX(); i++)
    {
      for(Int_t j=0; j<hist3->GetNbinsY(); j++)
	{
	  bin1 = hist3->GetBin(i,j);
	  content1 = hist3->GetBinContent(bin1);
	  if((Int_t)content1 < 3) continue;
	  if(content1 >= 1024)
	    content1 = 1023;
	  if(fNDigits >= fMaxDigits)
	    {
	      cerr<<"AliL3Compress::CreateDigits() : Array index out of range : "<<fNDigits<<endl;
	      return;
	    }
	  fDigits[fNDigits].fCharge=(Int_t)content1;
	  fDigits[fNDigits].fRow = row;
	  fDigits[fNDigits].fPad = (Int_t)hist3->GetXaxis()->GetBinCenter(i);
	  fDigits[fNDigits].fTime = (Int_t)hist3->GetYaxis()->GetBinCenter(j);
	  fDPt[fNDigits] = &fDigits[fNDigits];
	  fNDigits++;
	}
    }
  
  delete random;
  delete hist1;
  delete hist2;
  delete hist3;
  delete tr;
}

void AliL3Compress::QSort(AliL3RandomDigitData **a, Int_t first, Int_t last)
{
  
  // Sort array of AliL3RandomDigitData pointers using a quicksort algorithm.
  // Uses CompareDigits() to compare objects.
  // Thanks to Root!
  
  static AliL3RandomDigitData *tmp;
  static int i;           // "static" to save stack space
  int j;
  
  while (last - first > 1) {
    i = first;
    j = last;
    for (;;) {
      while (++i < last && CompareDigits(a[i], a[first]) < 0)
	;
      while (--j > first && CompareDigits(a[j], a[first]) > 0)
	;
      if (i >= j)
	break;
      
      tmp  = a[i];
      a[i] = a[j];
      a[j] = tmp;
    }
    if (j == first) {
      ++first;
      continue;
    }
    tmp = a[first];
    a[first] = a[j];
    a[j] = tmp;
    if (j - first < last - (j + 1)) {
      QSort(a, first, j);
      first = j + 1;   // QSort(j + 1, last);
    } else {
      QSort(a, j + 1, last);
      last = j;        // QSort(first, j);
    }
  }
}

void AliL3Compress::WriteRootFile(Char_t *digitsfile,Char_t *rootfile)
{
  Char_t fname[100];
  AliL3MemHandler *mem = new AliL3MemHandler();
  sprintf(fname,"%s/restored_%d_%d.raw",fPath,fSlice,fPatch);
  mem->SetBinaryInput(fname);
  UInt_t ndigits;
  AliL3DigitRowData *rowPt = (AliL3DigitRowData*)mem->CompBinary2Memory(ndigits);
  mem->CloseBinaryInput();

  AliL3FileHandler *file = new AliL3FileHandler();
  if(!file->SetAliInput(digitsfile))
    {
      cerr<<"AliL3Compress::WriteRootFile() : Error opening file: "<<digitsfile<<endl;
      return;
    }
  file->Init(fSlice,fPatch,NRows[fPatch]);
  file->AliDigits2RootFile(rowPt,rootfile);
  file->CloseAliInput();

  delete mem;
  delete file;
}
