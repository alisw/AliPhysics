/* $Id$
Author: Constantin Loizides <mailto: loizides@ikf.physik.uni-frankfurt.de>
*/

#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include "AliL3AltroMemHandler.h"

/** \class AliL3AltroMemHandler
// Converts digits in memory into a backlinked ALTRO like data format.
// The file misc/read.cxx shows how to use this class.
// As soon as you have more than 35 digits per pad, this class needs to
// improved to properly handle the backlinked list.
*/

ClassImp(AliL3AltroMemHandler)

AliL3AltroMemHandler::AliL3AltroMemHandler(){
  Clear();
};

Bool_t AliL3AltroMemHandler::Write(UShort_t row, UChar_t pad, UShort_t charge, UShort_t time)
{
  Bool_t redo=kFALSE;

  if((counter==0)||(tcounter==0)){ //new package or new sequence
    ltime=time;
    lpad=pad;
    lrow=row;
    altromem[counter]=charge;
    counter++;
    tcounter++;
  } else if ((row==lrow)&&(pad==lpad)) {
    if(time-ltime==1){ //add charge to sequence
      ltime=time;
      altromem[counter]=charge;
      counter++;
      tcounter++;

      if(counter==MAX_VALS){
	WriteFinal();
	Clear();
      }    
    } else { //finish old sequence
      altromem[counter]=ltime;
      counter++;
      altromem[counter]=tcounter+2;
      counter++;

      //start new sequence
      if(counter==MAX_VALS){ 
	WriteTrailer();
	Clear();
      } else tcounter=0;
      
      redo=kTRUE; 
    }
  } else { //finish old package
    WriteFinal();
    Clear();

    //start new package
    redo=kTRUE;
  }

  if(redo==kTRUE) Write(row,pad,charge,time);
  return kTRUE;
};

void AliL3AltroMemHandler::Clear(){
  memset(altromem,0,ALTRO_SIZE);
  counter=0;
  tcounter=0;
  lpad=0;
  ltime=0;
  lrow=0;
  flag=kFALSE;
};

Bool_t AliL3AltroMemHandler::WriteFinal(){
 if(counter>0){
   altromem[counter]=ltime;
   counter++;
   altromem[counter]=tcounter+2;
   counter++;
   WriteTrailer();
 }
 return kTRUE;
}

void AliL3AltroMemHandler::WriteTrailer(){

  UShort_t savecounter=counter;
  while(counter%4!=0){
    altromem[counter]=0x2AA;
    counter++;
  } 
  altromem[counter++]=0x2AA;
  altromem[counter++]=savecounter;
  altromem[counter++]=lpad;
  altromem[counter++]=lrow;

  UShort_t *ptr=altromem+counter;
  for (int i=0;i<counter;i++) {
    ptr--;
    if(flag==kTRUE) fwrite(ptr,sizeof(UShort_t),1,fOutBinary);
    else fprintf(fOutBinary,"%X ",*ptr);
  }
  if(flag==kFALSE) fprintf(fOutBinary,"\n");

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
  return kTRUE;
}
