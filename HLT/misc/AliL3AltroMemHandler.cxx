/* $Id$
Author: Constantin Loizides <mailto: loizides@ikf.physik.uni-frankfurt.de>
*/

#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include "AliL3AltroMemHandler.h"

/** \class AliL3AltroMemHandler
// Converts digits in memory into a backlinked ALTRO like data format.
// Its output file is used as input to the various VHDL testbenches.
// The file misc/read.cxx shows how to use this class.
*/

ClassImp(AliL3AltroMemHandler)

AliL3AltroMemHandler::AliL3AltroMemHandler(){
  Clear();
};

void AliL3AltroMemHandler::Clear(){
  memset(altromem,0,ALTRO_SIZE);
  memset(times_per_pad,0,1024);
  memset(charges_per_pad,0,1024);
  counter=ALTRO_SIZE;
  tcounter=0;
  lpad=0;
  lrow=0;
  flag=kFALSE;
};

void AliL3AltroMemHandler::Write(UShort_t row, UChar_t pad, UShort_t charge, UShort_t time)
{
  if(tcounter==0){
    lrow=row;
    lpad=pad;
  } else if((lrow!=row) || (lpad!=pad)){
    MakeAltroPackets(); //make packets
    Write();       //write packets
    Clear();            //clear up for next pad

    lrow=row;
    lpad=pad;
  }

  Add(charge,time);
}

void AliL3AltroMemHandler::Add(UShort_t charge, UShort_t time)
{
  times_per_pad[tcounter]=time;
  charges_per_pad[tcounter]=charge;
  tcounter++;
}

void AliL3AltroMemHandler::MakeAltroPackets()
{
  UShort_t i=0,j=0;
  UShort_t t=0,seqlength;
  UShort_t htime,ltime;
  int ddd=0;
  while(t<tcounter){
    pcounter=0;
    altromem[--counter]=0xFFFF; //mark return

    //make packets
    while((pcounter<ALTRO_SIZE) && (t<tcounter)){

      //find sequence
      i=t;
      ltime=times_per_pad[t];
      j=0;
      while((i+1<tcounter)&&(times_per_pad[i+1]==ltime+j+1)){
	i++;
	j++;
      }
      seqlength=j+1; //number of consecutive times   
      htime=times_per_pad[i]; //abs. time for sequence
/*
  cout << " counter " << counter << endl;
  if(htime!=ltime+j){
  cerr << "Error: " << ltime << " - " << htime << endl;
  exit(1);
  }
*/
      //store charges of sequence
      for(UShort_t k=0;k<seqlength;k++){
	altromem[--counter]=charges_per_pad[t];
	pcounter++;
	t++;
      }

      altromem[--counter]=htime;       //store abs. time of sequence
      pcounter++;
      altromem[--counter]=seqlength+2; //store length of sequence
      pcounter++;
    }

    AddTrailer();
  }
}

void AliL3AltroMemHandler::AddTrailer()
{
  UShort_t savepcounter=pcounter;

  while(pcounter%4!=0){
    altromem[--counter]=0x2AA;
    pcounter++;
  } 
  altromem[--counter]=0x2AA;
  altromem[--counter]=savepcounter;
  altromem[--counter]=lpad;
  altromem[--counter]=lrow;
}

void AliL3AltroMemHandler::Write()
{
  if(counter==ALTRO_SIZE) return;

  //save packets (reversed) to file
  UShort_t *ptr=altromem+counter;
  for (int i=counter;i<ALTRO_SIZE;i++) {

    if(flag==kTRUE){
      if (*ptr==0xFFFF) continue; //no return for end of packet
      fwrite(ptr,sizeof(UShort_t),1,fOutBinary);
    } else {
      if (*ptr==0xFFFF) fprintf(fOutBinary,"\n");
      else fprintf(fOutBinary,"%X ",*ptr);
    }

    ptr++;
  }
}

void AliL3AltroMemHandler::WriteFinal()
{
  if(tcounter>0){
    MakeAltroPackets();
    Write();
  }
}

/*
Bool_t AliL3AltroMemHandler::SetBinaryInput(FILE *file){
  fInBinary = file;
  if(!fInBinary){
    //LOG(AliL3Log::kWarning,"AliL3AltroMem::SetBinaryInput","File Open")<<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  return kTRUE;
}
*/

Bool_t AliL3AltroMemHandler::SetBinaryOutput(FILE *file){
  fOutBinary = file;
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3AltroMemHandler::SetBinaryOutput","File Open") <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  flag=kTRUE;
  return kTRUE;
}

Bool_t AliL3AltroMemHandler::SetASCIIOutput(FILE *file){
  fOutBinary = file;
  if(!fOutBinary){
    LOG(AliL3Log::kWarning,"AliL3AltroMemHandler::SetASCIIOutput","File Open") <<"Pointer to File = 0x0 "<<ENDLOG;
    return kFALSE;
  }
  flag=kFALSE;
  return kTRUE;
}
