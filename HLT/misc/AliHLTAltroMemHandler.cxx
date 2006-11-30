// @(#) $Id$

// Author: Constantin Loizides <mailto: loizides@ikf.uni-frankfurt.de>
// *-- Copyright & copy ALICE HLT Group
/** \class AliHLTAltroMemHandler
<pre>
//--------------------------------------------------------------------
// AliHLTAltroMemHandler
//
// Converts digits in memory into a backlinked ALTRO like data format.
// Its output file is used as input to the various VHDL testbenches.
// The file misc/read.cxx shows how to use this class.
</pre>
*/

#include "AliHLTStandardIncludes.h"

#include "AliHLTAltroMemHandler.h"
#include "AliHLTLogging.h"

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTAltroMemHandler)

AliHLTAltroMemHandler::AliHLTAltroMemHandler(){
  // default constructor
  Clear();
  ClearRead();
};

void AliHLTAltroMemHandler::Clear(){
  // clear everything
  memset(fAltroMem,0,ALTRO_SIZE);
  memset(fTimesPerPad,0,1024);
  memset(fChargesPerPad,0,1024);
  fCounter=ALTRO_SIZE;
  fTCounter=0;
  fLPad=0;
  fLRow=0;
  fFlag=kFALSE;
};

void AliHLTAltroMemHandler::ClearRead(){
  // clears the reading
  fRCounter=0;
  fSCounter=0;
  fRPad=0;
  fRRow=0;
  fRTime=0;
}

void AliHLTAltroMemHandler::Write(UShort_t row, UChar_t pad, UShort_t time, UShort_t charge)
{
  // writes packets
  if(fTCounter==0){
    fLRow=row;
    fLPad=pad;
  } else if((fLRow!=row) || (fLPad!=pad)){
    MakeAltroPackets(); //make packets
    Write();            //write packets
    Clear();            //clear up for next pad

    fLRow=row;
    fLPad=pad;
  }

  Add(charge,time);
}

void AliHLTAltroMemHandler::Add(UShort_t charge, UShort_t time)
{
  // adds new time and charge
  fTimesPerPad[fTCounter]=time;
  fChargesPerPad[fTCounter]=charge;
  fTCounter++;
}

void AliHLTAltroMemHandler::MakeAltroPackets()
{
  // makes Altro packets
  UShort_t i=0,j=0;
  UShort_t t=0,seqlength;
  UShort_t htime,ltime;

  while(t<fTCounter){
    fPCounter=0;
    fAltroMem[--fCounter]=0xFFFF; //mark return

    //make packets
    while((fPCounter<ALTRO_SIZE) && (t<fTCounter)){

      //find sequence
      i=t;
      ltime=fTimesPerPad[t];
      j=0;
      while((i+1<fTCounter)&&(fTimesPerPad[i+1]==ltime+j+1)){
	i++;
	j++;
      }
      seqlength=j+1; //number of consecutive times   
      htime=fTimesPerPad[i]; //abs. time for sequence

      //don't store sequence if it doesn't fit into packet
      if(fPCounter+seqlength>=ALTRO_SIZE) break;

      //store charges of sequence
      for(UShort_t k=0;k<seqlength;k++){
	fAltroMem[--fCounter]=fChargesPerPad[t];
	fPCounter++;
	t++;
      }

      fAltroMem[--fCounter]=htime;       //store abs. time of sequence
      fPCounter++;
      fAltroMem[--fCounter]=seqlength+2; //store length of sequence
      fPCounter++;
    }

    AddTrailer();
  }
}

void AliHLTAltroMemHandler::AddTrailer()
{
  // adds data trailer
  UShort_t savepcounter=fPCounter;

  while(fPCounter%4!=0){
    fAltroMem[--fCounter]=0x2AA;
    fPCounter++;
  } 
  fAltroMem[--fCounter]=0x2AA;
  fAltroMem[--fCounter]=savepcounter;
  fAltroMem[--fCounter]=fLPad;
  fAltroMem[--fCounter]=fLRow;
}

void AliHLTAltroMemHandler::Write()
{
  // default form of Write
  if(fCounter==ALTRO_SIZE) return;

  //save packets (reversed) to file
  UShort_t *ptr=fAltroMem+fCounter;
  for (int i=fCounter;i<ALTRO_SIZE;i++) {

    if(fFlag==kTRUE){
      if (*ptr==0xFFFF) continue; //no return for end of packet
      fwrite(ptr,sizeof(UShort_t),1,fOutBinary);
    } else {
      if (*ptr==0xFFFF) fprintf(fOutBinary,"\n");
      else fprintf(fOutBinary,"%X ",*ptr);
    }

    ptr++;
  }
}

void AliHLTAltroMemHandler::WriteFinal()
{
  // makes Altro packets and writes them
  if(fTCounter>0){
    MakeAltroPackets();
    Write();
  }
}

Bool_t AliHLTAltroMemHandler::Read(UShort_t &row, UChar_t &pad, UShort_t &time, UShort_t &charge)
{
  // reads the packets
  if(fFlag==kTRUE) {
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::Read","File Open")<<"Binary not supported!"<<ENDLOG;
    return kFALSE;
  }
  
  if(feof(fInBinary)){
    LOG(AliHLTLog::kDebug,"AliHLTAltroMemHandler::Read","File Open")<<"End of File!"<<ENDLOG;
    return kFALSE;
  }

  unsigned int dummy,dummy1;

  if(fRCounter==0){
    fscanf(fInBinary,"%x",&dummy);
    fRRow=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    fRPad=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    fRCounter=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    if(dummy!=682){
      if(feof(fInBinary)){
	ClearRead();
	return kFALSE;
      } else {
	LOG(AliHLTLog::kError,"AliHLTAltroMemHandler::Read","Trailer not found!")<<ENDLOG;
	return kFALSE;
      }
    }
    dummy1 = fRCounter;
    while(dummy1 % 4 != 0) {
      fscanf(fInBinary,"%x",&dummy);
      dummy1++;
    } 
  } 
  if(fSCounter==0){
    fscanf(fInBinary,"%x",&dummy);
    fSCounter=(UShort_t)dummy-2;
    fscanf(fInBinary,"%x",&dummy);
    fRTime=(UShort_t)dummy;
    fRCounter-=2;
  }
  fscanf(fInBinary,"%x",&dummy);

  row=fRRow;
  pad=fRPad;
  time=fRTime;
  charge=(UShort_t)dummy;

  fSCounter--;
  fRCounter--;
  fRTime--;
  return kTRUE;
}

Bool_t AliHLTAltroMemHandler::ReadSequence(UShort_t &row, UChar_t &pad, UShort_t &time, UChar_t &n, UShort_t **charges)
{
  // reads sequence
  if(fFlag==kTRUE) {
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::ReadSequence","File Open")<<"Binary not supported!"<<ENDLOG;
    return kFALSE;
  }
  
  if(feof(fInBinary)){
    LOG(AliHLTLog::kDebug,"AliHLTAltroMemHandler::ReadSequence","File Open")<<"End of File!"<<ENDLOG;
    return kFALSE;
  }

  unsigned int dummy,dummy1,dummy2;

  if(fRCounter==0){
    fscanf(fInBinary,"%x",&dummy);
    fRRow=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    fRPad=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    fRCounter=(UShort_t)dummy;
    fscanf(fInBinary,"%x",&dummy);
    if(dummy!=682){
      if(feof(fInBinary)){
	ClearRead();
	return kFALSE;
      } else {
	LOG(AliHLTLog::kError,"AliHLTAltroMemHandler::ReadSequence","Format") <<"Trailer not found!"<<ENDLOG;
	return kFALSE;
      }
    }
    dummy1 = fRCounter;
    while(dummy1 % 4 != 0) {
      fscanf(fInBinary,"%x",&dummy);
      dummy1++;
    } 
  } 

  fscanf(fInBinary,"%x",&dummy);
  fSCounter=(UShort_t)dummy-2;
  fscanf(fInBinary,"%x",&dummy);
  fRTime=(UShort_t)dummy;
  fRCounter-=2;

  if(n<fSCounter){
    if(*charges) delete[] *charges;
    *charges=new UShort_t[fSCounter];
  }

  dummy2=0;
  while(dummy2<fSCounter){
    fscanf(fInBinary,"%x",&dummy);
    (*charges)[dummy2]=(UShort_t)dummy;
    dummy2++;
    fRCounter--;
  }

  row=fRRow;
  pad=fRPad;
  time=fRTime;
  n=fSCounter;
  return kTRUE;
}

Bool_t AliHLTAltroMemHandler::SetBinaryInput(FILE *file)
{
  // sets binary input
  fInBinary = file;
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::SetBinaryInput","File Open")<<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  ClearRead();
  fFlag=kTRUE;
  return kTRUE;
}

Bool_t AliHLTAltroMemHandler::SetASCIIInput(FILE *file)
{
  // sets ASCII input
  fInBinary = file;
  if(!fInBinary){
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::SetASCIIInput","File Open")<<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  ClearRead();
  fFlag=kFALSE;
  return kTRUE;
}

Bool_t AliHLTAltroMemHandler::SetBinaryOutput(FILE *file){
  // sets binary output
  fOutBinary = file;
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::SetBinaryOutput","File Open") <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  Clear();
  fFlag=kTRUE;
  return kTRUE;
}

Bool_t AliHLTAltroMemHandler::SetASCIIOutput(FILE *file){
  // sets ASCII output
  fOutBinary = file;
  if(!fOutBinary){
    LOG(AliHLTLog::kWarning,"AliHLTAltroMemHandler::SetASCIIOutput","File Open") <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  Clear();
  fFlag=kFALSE;
  return kTRUE;
}

