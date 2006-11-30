// @(#) $Id$

// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>, Anders Vestbo <mailto:vestbo$fi.uib.no>, Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group 

/** \class AliHLTMemHandler 
<pre>
//_____________________________________________________________
// AliHLTMemHandler
//
// The HLT Binary File handler 
//
//  This class does all the memory I/O handling of HLT binary files.
//  
//  Examples:
//  ---------
//
//  1) Reading a binary file:
//  
//  AliHLTMemHandler file;
//  file.SetBinaryInput(filename);
//  file.Init(slice,patch);
//
//  UInt_t nrowss;
//  AliHLTDigitRowData *data = file.CompBinary2Memory(nrows);
//  
//  for(int i=0; i<nrows; i++) 
//    {
//    
//    AliHLTDigitData *dataPt = (AliHLTDigitData*)data->fDigitData;
//    for(int j=0; j<data->fNDigit; j++) 
//      {
//        pad = dataPt[j].fPad;
//        time = dataPt[j].fTime;
//        charge = dataPt[j].fCharge;
//      }
//     
//    file.UpdateRowPointer(data);
//  
//    }
//  file.CloseBinaryInput();
//  ________________________
//  
//  2) Writing a binary file:
//  
//  //First of all you need to store the data in memory,
//  //and have a pointer to it of type AliHLTDigitRowData.
//  //E.g. if you just want to write the data you read in example 1)
//  //into a new file, you can do the following:
//  
//  AliHLTMemHandler newfile;
//  newfile.Init(slice,patch);
//  newfile.SetBinaryOutput(newfilename);
//  newfile.Memory2CompBinary((UInt_t)NumberOfRowsInPatch,(AliHLTDigitRowData*)data);
//  newfile.CloseBinaryOutput();
//
//
// Compressed file format:
// -----------------------
//
// The data is RLE encoded and currently using _10_ bit range for the ADC-values.
</pre>
*/  

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"
#include "AliHLTDigitData.h"
#include "AliHLTLogging.h"
#include "AliHLTTransform.h"
#include "AliHLTTrackSegmentData.h"
#include "AliHLTSpacePointData.h"
#include "AliHLTTrackArray.h"
#include "AliHLTMemHandler.h"

#if __GNUC__ >= 3
using namespace std;
#endif
  
ClassImp(AliHLTMemHandler)
  
AliHLTMemHandler::AliHLTMemHandler()
{ 
  //Constructor
  fPt = 0;
  fSize =0;
  fInBinary = 0;
  fOutBinary = 0;
  fNRandom = 0;
  Init(0,0);
  fIsRandom = kFALSE;
  fRandomDigits = 0;
  fDPt =0;
  fNGenerate = 0;
  fNUsed = 0;
  fNDigits = 0;
  ResetROI();
}


AliHLTMemHandler::~AliHLTMemHandler()
{
  //Destructor
  if(fPt) delete[] fPt;
  if(fRandomDigits) delete [] fRandomDigits;
  if(fDPt) delete [] fDPt;
}

void AliHLTMemHandler::Init(Int_t s,Int_t p, Int_t *r)
{
  //init handler
  fSlice=s;fPatch=p;
  if(r) {
    fRowMin=r[0];
    fRowMax=r[1];
  }else{
     fRowMin=AliHLTTransform::GetFirstRow(p);
     fRowMax=AliHLTTransform::GetLastRow(p); 
  }
  ResetROI();
}

void AliHLTMemHandler::ResetROI()
{
  //Resets the Look-up table for Region of Interest mode.
  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      fEtaMinTimeBin[i] = 0;
      fEtaMaxTimeBin[i] = AliHLTTransform::GetNTimeBins()-1;
    }
}

void AliHLTMemHandler::SetROI(Float_t *eta,Int_t */*slice*/)
{
  // Init the Look-up table for the Region of Interest mode.
  //   Here you can specify a certain etaregion, - all data
  //   outside this region will be discarded:
  //   eta[0] = mimium eta
  //   eta[1] = maximum eta
  //   slice[0] = mimumum slice
  //   slice[1] = maximum slice


  if(eta[1]==0)
    {
      LOG(AliHLTLog::kWarning,"AliHLTMemHandler::SetROI","Eta Values")
	<<"Bad ROI parameters. IDIOT! "<<ENDLOG;
      for(Int_t i=fRowMin; i<=fRowMax; i++)
	{
	  fEtaMinTimeBin[i]=0;
	  fEtaMaxTimeBin[i]=0;
	}
      return;
    }
  
  for(Int_t i=fRowMin; i<=fRowMax; i++)
    {
      Int_t sector,row;
      Float_t xyz[3];
      
      Float_t thetamax = 2*atan(exp(-1.*eta[1]));
      
      xyz[0] = AliHLTTransform::Row2X(i);
      xyz[1]=0;
      xyz[2] = xyz[0]/tan(thetamax);
      AliHLTTransform::Slice2Sector(fSlice,i,sector,row);
      AliHLTTransform::Local2Raw(xyz,sector,row);
      
      fEtaMinTimeBin[i] = (Int_t)xyz[2];
      
      if(eta[0]==0)
	fEtaMaxTimeBin[i] = 445;
      else
	{
	  Float_t thetamin = 2*atan(exp(-1.*eta[0]));
	  xyz[0] = AliHLTTransform::Row2X(i);
	  xyz[1] = AliHLTTransform::GetMaxY(i);
	  Float_t radii = sqrt(pow(xyz[0],2) + pow(xyz[1],2));
	  xyz[2] = radii/tan(thetamin);
	  AliHLTTransform::Local2Raw(xyz,sector,row);
	  fEtaMaxTimeBin[i] = (Int_t)xyz[2];
	}
    }
  
}

Bool_t AliHLTMemHandler::SetBinaryInput(char *name)
{
  //Set the input binary file.
  fInBinary = fopen(name,"r");
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::SetBinaryInput","File Open")
      <<"Error opening file "<<name<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTMemHandler::SetBinaryInput(FILE *file)
{
  //Set the input binary file.
  fInBinary = file;
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::SetBinaryInput","File Open")
    <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliHLTMemHandler::CloseBinaryInput()
{
  //Close the input file.
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CloseBinaryInput","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fInBinary);
  fInBinary =0;
}

Bool_t AliHLTMemHandler::SetBinaryOutput(char *name)
{
  //Set the binary output file.
    fOutBinary = fopen(name,"w");
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::SetBinaryOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

Bool_t AliHLTMemHandler::SetBinaryOutput(FILE *file)
{
  //Set the binary output file.
    fOutBinary = file;
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::SetBinaryOutput","File Open")
      <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}

void AliHLTMemHandler::CloseBinaryOutput()
{
  //close binary  
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CloseBinaryOutPut","File Close")
      <<"Nothing to Close"<<ENDLOG;
    return;
  }
  fclose(fOutBinary);
  fOutBinary =0;
}

UInt_t AliHLTMemHandler::GetFileSize()
{
  //Returns the file size in bytes of the input file.
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::GetFileSize","File")
      <<"No Input File"<<ENDLOG;
    return 0;
  }
  fseek(fInBinary,0,SEEK_END);
  UInt_t size = (UInt_t) ftell(fInBinary);
  rewind(fInBinary);
  return size; 
}

Byte_t *AliHLTMemHandler::Allocate()
{
  //Allocate
  return Allocate(GetFileSize()); 
}

Byte_t *AliHLTMemHandler::Allocate(AliHLTTrackArray *array)
{
  //Allocate memory for tracks in memory. Used by TrackArray2Binary()
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Allocate","Memory")
      <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return 0;
  }
  return Allocate(array->GetOutSize()); 
}

Byte_t *AliHLTMemHandler::Allocate(UInt_t size)
{
  //Allocate memory of size in bytes.
  if(fPt){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Allocate","Memory")
      <<"Delete Memory"<<ENDLOG;
    Free();
  } 
  fPt = new Byte_t[size];
  fSize = size;
  memset(fPt,0,fSize);
  LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Allocate","Memory")
  <<AliHLTLog::kDec<<"Allocate "<<size<<" Bytes of Memory"<<ENDLOG;
  return fPt;
}

void AliHLTMemHandler::Free()
{
  //Clear the memory, if allocated.
  if(!fPt){
    //    LOG(AliHLTLog::kInformational,"AliHLTMemHandler::Free","Memory")
    //      <<"No Memory allocated - can't Free"<<ENDLOG;
    return;
  }  
  delete[] fPt;
  fPt = 0;
  fSize =0;
}

///////////////////////////////////////// Random
void AliHLTMemHandler::SetRandomSeed()
{
  //If you are adding random data to the original data.
  time_t *tp=0;
  SetRandomSeed(time(tp));
}

void AliHLTMemHandler::SetRandomCluster(Int_t maxnumber)
{
  //If you are adding random data to the original data.
  
  fIsRandom = kTRUE;
  fNRandom = maxnumber;
  fNDigits = 0;
  if(fRandomDigits) delete [] fRandomDigits;
  fRandomDigits = new AliHLTRandomDigitData[fNRandom*9];
  if(fDPt) delete [] fDPt;
  fDPt = new AliHLTRandomDigitData *[fNRandom*9];
}

void AliHLTMemHandler::QSort(AliHLTRandomDigitData **a, Int_t first, Int_t last)
{

   // Sort array of AliHLTRandomDigitData pointers using a quicksort algorithm.
   // Uses CompareDigits() to compare objects.
   // Thanks to Root!

   static AliHLTRandomDigitData *tmp;
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

UInt_t AliHLTMemHandler::GetRandomSize() const
{
  //get random size
  Int_t nrandom = 0;
  for(Int_t r=fRowMin;r<=fRowMax;r++){
    Int_t npad=AliHLTTransform::GetNPads(r);
    nrandom  += Int_t (fNGenerate * ((Double_t) npad/141.));
  }
  return 9 * nrandom * sizeof(AliHLTDigitData);
}

void AliHLTMemHandler::Generate(Int_t row)
{
  //Generate random data on row, if you didn't 
  //ask for this, nothing happens here.
  
  if(!fIsRandom) return;
  ResetRandom();
  fNDigits = 0;
  Int_t npad=AliHLTTransform::GetNPads(row);
  Int_t ntime = fEtaMaxTimeBin[row] - fEtaMinTimeBin[row];
  Int_t nrandom  = Int_t (fNGenerate * ((Double_t) npad/141.) * 
			  (Double_t) ntime/(Double_t) AliHLTTransform::GetNTimeBins() );
  
  for(Int_t n=0;n<nrandom;n++){
    Int_t pad = (int)((float)rand()/RAND_MAX*npad);
    Int_t time =(int)((float)rand()/RAND_MAX*ntime+fEtaMinTimeBin[row] );
    Int_t charge = (int)((float)rand()/RAND_MAX*AliHLTTransform::GetADCSat());
    DigitizePoint(row,pad,time,charge);
  }
  QSort(fDPt,0,fNDigits);
}


void AliHLTMemHandler::DigitizePoint(Int_t row, Int_t pad, 
				    Int_t time,Int_t charge)
{
  //Making one single random cluster.
  for(Int_t j=-1;j<2;j++){
    for(Int_t k=-1;k<2;k++){
      Int_t dcharge = charge;
      if(j) dcharge /=2;
      if(k) dcharge /=2;
      if(dcharge<10) continue;
      Int_t dpad  = j + pad;
      Int_t dtime = k + time;
      
      if(dpad<0||dpad>=AliHLTTransform::GetNPads(row))  continue;
      if(dtime<0||dtime>=AliHLTTransform::GetNTimeBins()) continue;
      
      fRandomDigits[fNDigits].fCharge = dcharge;
      fRandomDigits[fNDigits].fRow = row;
      fRandomDigits[fNDigits].fPad = dpad;
      fRandomDigits[fNDigits].fTime = dtime;
      fDPt[fNDigits] = &fRandomDigits[fNDigits];
      fNDigits++;
    }
  }
}

///////////////////////////////////////// Digit IO  
Bool_t AliHLTMemHandler::Memory2Binary(UInt_t nrow,AliHLTDigitRowData *data)
{
  //Write data to the outputfile as is. No run-length encoding is done.
  
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","File")
      <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","Memory")
      <<"Pointer to AliHLTDigitRowData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  
  AliHLTDigitRowData *rowPt = data; 
  Int_t outsize = 0;
  for(UInt_t i=0;i<nrow;i++){
    Int_t size = sizeof(AliHLTDigitData) * rowPt->fNDigit 
      + sizeof(AliHLTDigitRowData);
    outsize += size;
    fwrite(rowPt,size,1,fOutBinary);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTDigitRowData *) bytePt;
  }
  LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Memory2Binary","Memory")
    <<AliHLTLog::kDec<<"Wrote "<<outsize<<" Bytes to Memory ("
    <<nrow<<" Rows)"<<ENDLOG;
  return kTRUE;
}

Bool_t AliHLTMemHandler::Binary2Memory(UInt_t & nrow,AliHLTDigitRowData *data)
{
  //Read inputfile into memory as is, and store it in data. 
  // No run-length encoding is assumed.

  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","File")
      <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","Memory")
      <<"Pointer to AliHLTDigitRowData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  rewind(fInBinary);
  AliHLTDigitRowData *rowPt = data;
  UInt_t rowcount = 0;
  Int_t outsize =0;
  while(!feof(fInBinary)){
    Byte_t  *bytePt =(Byte_t *) rowPt;

    if(fread(rowPt,sizeof(AliHLTDigitRowData),1,fInBinary)!=1) break;

    bytePt += sizeof(AliHLTDigitRowData);
    outsize += sizeof(AliHLTDigitRowData);

    Int_t size = sizeof(AliHLTDigitData) * rowPt->fNDigit;

    //if(fread(bytePt,size,1,fInBinary)!=1) break;
    fread(bytePt,size,1,fInBinary);
    bytePt += size;
    outsize += size;
    rowPt = (AliHLTDigitRowData *) bytePt;
    rowcount++;
  }  
  nrow= rowcount;
    LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Binary2Memory","Memory")
    <<AliHLTLog::kDec<<"Wrote "<<outsize<<" Bytes to Memory ("
    <<rowcount<<" Rows)"<<ENDLOG;
  return kTRUE;
}

void AliHLTMemHandler::AddData(AliHLTDigitData *data,UInt_t & ndata,
			      UInt_t /*row*/,UShort_t pad,UShort_t time,UShort_t charge) const
{
  //add some data
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  ndata++;
}

void AliHLTMemHandler::AddRandom(AliHLTDigitData *data, UInt_t & ndata)
{
  //add some random data
  data[ndata].fPad = fDPt[fNUsed]->fPad;
  data[ndata].fTime = fDPt[fNUsed]->fTime;
  data[ndata].fCharge = fDPt[fNUsed]->fCharge;
  ndata++;
  fNUsed++;
}

void AliHLTMemHandler::MergeDataRandom(AliHLTDigitData *data, UInt_t & ndata,
				      UInt_t row, UShort_t pad, UShort_t time, UShort_t charge)
{
  //merge random data
  data[ndata].fPad = pad;
  data[ndata].fTime = time;
  data[ndata].fCharge = charge;
  while(ComparePoints(row,pad,time)==0){
    Int_t ch = data[ndata].fCharge + fDPt[fNUsed]->fCharge;
    if(charge>=AliHLTTransform::GetADCSat()) ch = AliHLTTransform::GetADCSat();
    data[ndata].fCharge = ch;
    fNUsed++;
  }
  ndata++;
}

void AliHLTMemHandler::AddDataRandom(AliHLTDigitData *data, UInt_t & ndata,
                   UInt_t row, UShort_t pad, UShort_t time, UShort_t charge)
{
  //add data random
  Int_t action;
  while((action=ComparePoints(row,pad,time))==1){
    AddRandom(data,ndata);
  }
  if(action==0){
    MergeDataRandom(data,ndata,row,pad,time,charge);
  }
  if(action<0){
    AddData(data,ndata,row,pad,time,charge);
  }  
}

void AliHLTMemHandler::Write(UInt_t *comp, UInt_t & index, 
			    UInt_t & subindex, UShort_t value) const
{
  //write compressed data
  UInt_t shift[3] = {0,10,20};
  if(subindex==0) comp[index] =0; //clean up memory
  comp[index] |= (value&0x03ff)<<shift[subindex];
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
}

UShort_t AliHLTMemHandler::Read(UInt_t *comp, UInt_t & index, UInt_t & subindex) const
{ 
  //read compressed data
  UInt_t shift[3] = {0,10,20};
  UShort_t value = (comp[index]>>shift[subindex])&0x03ff;
  if(subindex == 2){
    subindex = 0;
    index++;
  }
  else subindex++;
  
  return value;
}

UShort_t AliHLTMemHandler::Test(UInt_t *comp, 
			       UInt_t index, UInt_t  subindex) const
{
  //supi dupi test
  UInt_t shift[3] = {0,10,20};
  return (comp[index]>>shift[subindex])&0x03ff;
}

Int_t AliHLTMemHandler::Memory2CompMemory(UInt_t nrow,
					 AliHLTDigitRowData *data,UInt_t *comp)
{
  //Performs run-length encoding on data stored in memory pointed to by data.
  //The compressed data is written to comp.
  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2CompMemory","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2CompMemory","Memory")
      <<"Pointer to AliHLTDigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliHLTDigitRowData *rowPt = data;
  UInt_t index=0;
  UInt_t subindex=0;
  
  for(UInt_t i=0;i<nrow;i++){
    UShort_t value = rowPt->fRow;
    Write(comp,index,subindex,value);
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<rowPt->fNDigit;dig++){
      if(rowPt->fDigitData[dig].fPad <200){ 
        ddd[rowPt->fDigitData[dig].fPad]++;
      }
    }
    for(Int_t d=0;d<200;d++){ 
      if(ddd[d]){
        npad++;
        maxpad =d;
      }
    }
    Write(comp,index,subindex,npad);
    UInt_t digit=0;
    for(UShort_t pad=0;pad <= maxpad;pad++){
      if(digit>=rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
        continue;
      Write(comp,index,subindex,pad);
//    write zero if time != 0
      if(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        if(rowPt->fDigitData[digit].fTime>0){
          Write(comp,index,subindex,0);
          Write(comp,index,subindex,rowPt->fDigitData[digit].fTime);
        }
      }
      while(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        UShort_t charge = rowPt->fDigitData[digit].fCharge;
        if(charge>=1023){
          charge=1023;
        }
        Write(comp,index,subindex,charge);
        if(digit+1<rowPt->fNDigit&&rowPt->fDigitData[digit+1].fPad == pad){
          if(rowPt->fDigitData[digit].fTime +1 !=
                     rowPt->fDigitData[digit+1].fTime){
            Write(comp,index,subindex,0);
            UShort_t nzero = rowPt->fDigitData[digit+1].fTime - 
                             (rowPt->fDigitData[digit].fTime +1);
            Write(comp,index,subindex,nzero);
          }  
        }
        digit++;
      }
      Write(comp,index,subindex,0);
      Write(comp,index,subindex,0);
    }
    
    Int_t size = sizeof(AliHLTDigitData) * rowPt->fNDigit+
                                            sizeof(AliHLTDigitRowData);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTDigitRowData *) bytePt;
  }
  while(subindex)
    Write(comp,index,subindex,0);
  return index * sizeof(UInt_t);
}

Int_t AliHLTMemHandler::CompMemory2Memory(UInt_t  nrow,
					 AliHLTDigitRowData *data,UInt_t *comp)
{
  //Uncompress the run-length encoded data in memory pointed to by comp, and
  //  store it in data.

  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompMemory2Memory","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompMemory2Memory","Memory")
      <<"Pointer to AliHLTDigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  Int_t outsize=0;
  
  AliHLTDigitRowData *rowPt = data;
  UInt_t index=0;
  UInt_t subindex=0;
  
  for(UInt_t i=0;i<nrow;i++){
    UInt_t ndigit=0;
    UInt_t row =Read(comp,index,subindex);
    rowPt->fRow=row;
    Generate(row);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      UShort_t charge;
      UShort_t time =0;
      UShort_t pad = Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if( (time = Read(comp,index,subindex)) == 0 ){
          continue;
        }
      }
      for(;;){
        while( (charge=Read(comp,index,subindex)) != 0){
          if(time>=fEtaMinTimeBin[row]&&time<=fEtaMaxTimeBin[row])
	    //AddData(rowPt->fDigitData,ndigit,row,pad,time,charge);
	    //seems we are using this function... but dont know why
            AddDataRandom(rowPt->fDigitData,ndigit,row,pad,time,charge);
          time++;
        }
        UShort_t tshift = Read(comp,index,subindex);
        if(tshift == 0) break;
        time += tshift;
      }
    }
    rowPt->fNDigit = ndigit;
    Int_t size = sizeof(AliHLTDigitData) * rowPt->fNDigit+
      sizeof(AliHLTDigitRowData);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    outsize += size;
    rowPt = (AliHLTDigitRowData *) bytePt;
  }
  return outsize;
}

UInt_t AliHLTMemHandler::GetCompMemorySize(UInt_t nrow,AliHLTDigitRowData *data) const
{
  //Return the size of RLE data, after compressing data.
  
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::GetCompMemorySize","Memory")
      <<"Pointer to AliHLTDigitRowData = 0x0 "<<ENDLOG;
    return 0;
  }
  AliHLTDigitRowData *rowPt = data;
  UInt_t index=0;
  
  for(UInt_t i=0;i<nrow;i++){
    index++;
    UShort_t maxpad=0; 
    UShort_t npad=0;
    Int_t ddd[1000];
    for(Int_t d=0;d<200;d++) ddd[d]=0;
    for(UInt_t dig=0;dig<rowPt->fNDigit;dig++){
      if(rowPt->fDigitData[dig].fPad <200){ 
        ddd[rowPt->fDigitData[dig].fPad]++;
      }
    }
    for(Int_t d=0;d<200;d++){ 
      if(ddd[d]){
        npad++;
        maxpad =d;
      }
    }
    index++;
    UInt_t digit=0;
    for(UShort_t pad=0;pad <= maxpad;pad++){
      if(digit>=rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
        continue;
      index++;
      //    write zero if time != 0
      if(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        if(rowPt->fDigitData[digit].fTime>0){
          index++;
          index++;
        }
      }
      while(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad){
        index++;
        if(digit+1<rowPt->fNDigit&&rowPt->fDigitData[digit+1].fPad == pad){
          if(rowPt->fDigitData[digit].fTime +1 !=
                     rowPt->fDigitData[digit+1].fTime){
            index++;
            index++;
          }  
        }
        digit++;
      }
      index++;
      index++;
    }

    Int_t size = sizeof(AliHLTDigitData) * rowPt->fNDigit+
                                            sizeof(AliHLTDigitRowData);
    Byte_t  *bytePt =(Byte_t *) rowPt;
    bytePt += size;
    rowPt = (AliHLTDigitRowData *) bytePt;
  }
  while(index%3)
    index++;
  return (index/3) * sizeof(UInt_t);
}

UInt_t AliHLTMemHandler::GetMemorySize(UInt_t nrow,UInt_t *comp) const
{
  //get memory size
  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::GetMemorySize","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  Int_t outsize=0;

  UInt_t index=0;
  UInt_t subindex=0;

  for(UInt_t i=0;i<nrow;i++){
    UInt_t ndigit=0;
    Read(comp,index,subindex);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if(Read(comp,index,subindex)== 0) continue;
      }
      for(;;){
        while(Read(comp,index,subindex)!=0) ndigit++;
        if(Read(comp,index,subindex)==0) break;
      }
    }
    Int_t size = sizeof(AliHLTDigitData) * ndigit+
                                        sizeof(AliHLTDigitRowData);
    outsize += size;
  }
   
  return outsize;
}

UInt_t AliHLTMemHandler::GetNRow(UInt_t *comp,UInt_t size)
{
  //get number of rows
  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::GetNRow","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return 0;
  }
  size = size /4;
  UInt_t nrow=0;
  UInt_t index=0;
  UInt_t subindex=0;
  while(index<size-1){ //don't start with last word
    nrow++;
    UInt_t ndigit=0;
    Read(comp,index,subindex);
    UShort_t npad = Read(comp,index,subindex);
    for(UShort_t p=0;p<npad;p++){
      Read(comp,index,subindex);
      if(Test(comp,index,subindex)==0){
        Read(comp,index,subindex);
        if(Read(comp,index,subindex)==0)continue;
      }
      for(;;){
        while(Read(comp,index,subindex)!=0) ndigit++;
        if(Read(comp,index,subindex)==0) break;
      }
    }
  }
  if(index==size-1){  //last word
    if(subindex<2){
      if(Read(comp,index,subindex)!=0) nrow++;
    }
  }
  return nrow;
}

Bool_t AliHLTMemHandler::CompMemory2CompBinary(UInt_t nrow,UInt_t *comp,
					      UInt_t size)
{
  //Write the RLE data in comp to the output file.
  
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompMemory2CompBinary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompMemory2CompBinary","Memory")
    <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(size==0)
    size=GetMemorySize(nrow,comp);
  if(!size){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompMemory2CompBinary","Memory")
    <<"Memory size = 0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t length = size/sizeof(UInt_t);
  fwrite(&length,sizeof(UInt_t),1,fOutBinary);  
  fwrite(comp,size,1,fOutBinary);
  return kTRUE;
}

Bool_t AliHLTMemHandler::CompBinary2CompMemory(UInt_t & nrow,UInt_t *comp)
{
  //Read the RLE data from file, and store it in comp. No unpacking yet.

  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompBinary2CompMemory","File")
      <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!comp){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::CompBinary2CompMemory","Memory")
      <<"Pointer to compressed data = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  rewind(fInBinary);
  UInt_t length;
  if(fread(&length,sizeof(UInt_t),1,fInBinary)!=1) return kFALSE;
  UInt_t size = length*sizeof(UInt_t);
  if(fread(comp,size,1,fInBinary)!=1) return kFALSE;
  // now find the number of dig
  nrow =  GetNRow(comp,size);
  return kTRUE;
}

AliHLTDigitRowData *AliHLTMemHandler::CompBinary2Memory(UInt_t & nrow)
{
  // Read the RLE inputfile, unpack it and return the pointer to it.
  AliHLTMemHandler * handler = new AliHLTMemHandler();
  handler->SetBinaryInput(fInBinary);
  UInt_t *comp =(UInt_t *)handler->Allocate();
  handler->CompBinary2CompMemory(nrow,comp);
  UInt_t size = GetMemorySize(nrow,comp);
  AliHLTDigitRowData *data = (AliHLTDigitRowData *)Allocate(size);
  CompMemory2Memory(nrow,data,comp);
  handler->Free();
  delete handler;
  return data;  
}

Bool_t AliHLTMemHandler::Memory2CompBinary(UInt_t nrow,AliHLTDigitRowData *data)
{
  //Perform RLE on the data, and write it to the output file.
  Bool_t out = kTRUE;
  AliHLTMemHandler * handler = new AliHLTMemHandler();
  UInt_t size = GetCompMemorySize(nrow,data);
  UInt_t *comp =(UInt_t *)handler->Allocate(size);
  Memory2CompMemory(nrow,data,comp);
  CompMemory2CompBinary(nrow,comp,size);
  handler->Free();
  delete handler;
  return out;
}


///////////////////////////////////////// Point IO  
Bool_t AliHLTMemHandler::Memory2Binary(UInt_t npoint,AliHLTSpacePointData *data)
{
  //Writing spacepoints stored in data to the outputfile.
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","File")
      <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","Memory")
      <<"Pointer to AliHLTSpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  UInt_t size = npoint*sizeof(AliHLTSpacePointData);
  fwrite(data,size,1,fOutBinary);
  
  return kTRUE;
}

Bool_t AliHLTMemHandler::Transform(UInt_t npoint,AliHLTSpacePointData *data,Int_t slice)
{
  //Transform the space points in data, to global coordinates in slice.
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Transform","Memory")
    <<"Pointer to AliHLTSpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  
  for(UInt_t i=0;i<npoint;i++){
    Float_t xyz[3];
    xyz[0] = data[i].fX;
    xyz[1] = data[i].fY;
    xyz[2] = data[i].fZ;
    AliHLTTransform::Local2Global(xyz,slice);
    data[i].fX = xyz[0];
    data[i].fY = xyz[1];
    data[i].fZ = xyz[2];
  }
  return kTRUE;
}

Bool_t AliHLTMemHandler::Binary2Memory(UInt_t & npoint,AliHLTSpacePointData *data)
{
  //Read the space points in inputfile, and store it in data.
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","Memory")
    <<"Pointer to AliHLTSpacePointData = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  Int_t size = GetFileSize(); 
  npoint = size/sizeof(AliHLTSpacePointData);
  if(size==0) {
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","File")
    <<"File Size == 0"<<ENDLOG;
    return kFALSE;
  }

  if(fread(data,size,1,fInBinary)!=1){
    LOG(AliHLTLog::kFatal,"AliHLTMemHandler::Binary2Memory","File")
    <<"File Read Error "<<ENDLOG;
    return kFALSE;
  }
  if(size%sizeof(AliHLTSpacePointData)){
    LOG(AliHLTLog::kFatal,"AliHLTMemHandler::Binary2Memory","File Size")
    <<"File Size wrong "<<ENDLOG;
    return kFALSE; 
  }
  LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Binary2Memory","File")
  <<AliHLTLog::kDec<<"Wrote  "<<size<<" Bytes to Memory"<<ENDLOG;
  return kTRUE;
}

///////////////////////////////////////// Track IO  
Bool_t AliHLTMemHandler::Memory2Binary(UInt_t ntrack,AliHLTTrackSegmentData *data)
{
  //Write the tracks stored in data, to outputfile.
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2Binary","Memory")
    <<"Pointer to AliHLTTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliHLTTrackSegmentData *trackPt = data;
  for(UInt_t i=0;i<ntrack;i++){
    Int_t size=sizeof(AliHLTTrackSegmentData)+trackPt->fNPoints*sizeof(UInt_t); 
    fwrite(trackPt,size,1,fOutBinary);
    Byte_t *bytePt = (Byte_t*) trackPt;
    bytePt += size; 
    trackPt = (AliHLTTrackSegmentData*) bytePt;
  }
  LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Memory2Binary","File")
  <<AliHLTLog::kDec<<"Wrote  "<<ntrack<<" Tracks to File"<<ENDLOG;
  
  return kTRUE;
}

Bool_t AliHLTMemHandler::Binary2Memory(UInt_t & ntrack,AliHLTTrackSegmentData *data)
{
  //Read the tracks in inputfile, and store it in data.
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2Memory","Memory")
    <<"Pointer to AliHLTTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  ntrack=0;
  AliHLTTrackSegmentData *trackPt = data;
  rewind(fInBinary);

  while(!feof(fInBinary)){
    if(fread(trackPt,sizeof(AliHLTTrackSegmentData),1,fInBinary)!=1) break;
    Int_t size=trackPt->fNPoints*sizeof(UInt_t);
    if(fread(trackPt->fPointIDs,size,1,fInBinary)!=1) break;
    Byte_t *bytePt = (Byte_t*) trackPt;
    bytePt += sizeof(AliHLTTrackSegmentData)+size;
    trackPt = (AliHLTTrackSegmentData*) bytePt;
    ntrack++; 
  }
  LOG(AliHLTLog::kDebug,"AliHLTMemHandler::Binary2Memory","File")
  <<AliHLTLog::kDec<<"Wrote  "<<ntrack<<" Tracks to Memory"<<ENDLOG;
  return kTRUE;
}

Bool_t AliHLTMemHandler::TrackArray2Binary(AliHLTTrackArray *array)
{
  //Write the trackarray to the outputfile.
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::TrackArray2Binary","File")
    <<"No Output File"<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::TrackArray2Binary","Memory")
    <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliHLTTrackSegmentData *data = (AliHLTTrackSegmentData *)Allocate(array);

  UInt_t ntrack;
  TrackArray2Memory(ntrack,data,array);
  Memory2Binary(ntrack,data);
  Free();
  return kTRUE;
}

Bool_t AliHLTMemHandler::Binary2TrackArray(AliHLTTrackArray *array)
{
  //Read the tracks in inputfile, and fill it in trackarray. 
  //array should already be constructed.
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2TrackArray","File")
    <<"No Input File"<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Binary2TrackArray","Memory")
    <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  AliHLTTrackSegmentData *data = (AliHLTTrackSegmentData *)Allocate();
  UInt_t ntrack;
  Binary2Memory(ntrack,data);
  Memory2TrackArray(ntrack,data,array);  
  Free();
  return kTRUE;
}

Bool_t AliHLTMemHandler::TrackArray2Memory(UInt_t & ntrack,AliHLTTrackSegmentData *data,AliHLTTrackArray *array) const
{
  //Fill the trackarray into the AliTrackSegmentData structures before writing to outputfile.
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliHLTTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::TrackArray2Memory","Memory")
    <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }

  array->WriteTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliHLTMemHandler::Memory2TrackArray(UInt_t ntrack,AliHLTTrackSegmentData *data,AliHLTTrackArray *array) const
{
  //Fill the tracks in data into trackarray.
  
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data);
  return kTRUE;
}

Bool_t AliHLTMemHandler::Memory2TrackArray(UInt_t ntrack,AliHLTTrackSegmentData *data,AliHLTTrackArray *array,Int_t slice) const
{
  //Fill the tracks in data into trackarray, and rotate the tracks to global coordinates.
    
  if(!data){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTrackSegmentData = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  if(!array){
    LOG(AliHLTLog::kWarning,"AliHLTMemHandler::Memory2TrackArray","Memory")
    <<"Pointer to AliHLTTrackArray = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  array->FillTracks(ntrack,data,slice);
  return kTRUE;
}

void AliHLTMemHandler::UpdateRowPointer(AliHLTDigitRowData *&tempPt)
{
  //Update the data pointer to the next padrow in memory.
  
  Byte_t *tmp = (Byte_t*)tempPt;
  Int_t size = sizeof(AliHLTDigitRowData) + tempPt->fNDigit*sizeof(AliHLTDigitData);
  tmp += size;
  tempPt = (AliHLTDigitRowData*)tmp;
}

Int_t  AliHLTMemHandler::ComparePoints(UInt_t /*row*/,UShort_t pad,UShort_t time) const
{
  //compare two points
  if(fNUsed>=fNDigits) return -2;

  if(pad==fDPt[fNUsed]->fPad&&time==fDPt[fNUsed]->fTime) return 0;

  if(pad<fDPt[fNUsed]->fPad) return -1;
  if(pad==fDPt[fNUsed]->fPad&&time<fDPt[fNUsed]->fTime)  return -1;

  return 1;
}

Int_t AliHLTMemHandler::CompareDigits(AliHLTRandomDigitData *a,AliHLTRandomDigitData *b) const
{
  //compare two digits
  if(a->fPad==b->fPad && a->fTime == b->fTime) return 0;

  if(a->fPad<b->fPad) return -1;
  if(a->fPad==b->fPad && a->fTime<b->fTime) return -1;
  
  return 1;
}
