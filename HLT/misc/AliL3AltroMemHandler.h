#ifndef ALIL3ALTROMEMHANDLER_H
#define ALIL3ALTROMEMHANDLER_H

#include "AliL3RootTypes.h"
#include "AliL3Logging.h"

//Maximum Size of Altro Packet equals 1000 Bit
#define ALTRO_SIZE 125
//Maximum 10Bit data to be stored in one Packet
#define MAX_VALS 94

class AliL3AltroMemHandler {

  public: 
   AliL3AltroMemHandler();
   Bool_t Write(UShort_t row, UChar_t pad, UShort_t charge, UShort_t time);
   Bool_t WriteFinal();
   Bool_t SetBinaryOutput(FILE *file);
   //Bool_t SetBinaryInput(FILE *file);
   void SetBinary(Bool_t flag_=kTRUE){flag=flag_;};

  private:
   UShort_t altromem[ALTRO_SIZE];
   //FILE *fInBinary;
   FILE *fOutBinary;
   UShort_t lrow;
   UChar_t lpad;
   UShort_t ltime;
   UShort_t counter;
   UShort_t tcounter;
   Bool_t flag;
   void WriteTrailer();
   void Clear();
   ClassDef(AliL3AltroMemHandler,1)
};

#endif






