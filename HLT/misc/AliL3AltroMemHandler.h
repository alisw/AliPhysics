// @(#) $Id$

#ifndef ALIL3ALTROMEMHANDLER_H
#define ALIL3ALTROMEMHANDLER_H

#include "AliL3RootTypes.h"

//Maximum size of Altro packet equals 1000 bit
#define ALTRO_PACKET_SIZE 125

//Maximum 10 bit data to be stored in one packet
#define MAX_VALS 94

//Maximum size of array to store whole pad
#define ALTRO_SIZE (100*ALTRO_PACKET_SIZE)

class AliL3AltroMemHandler {

  public: 

   AliL3AltroMemHandler();
   void Write(UShort_t row, UChar_t pad, UShort_t time, UShort_t charge);
   Bool_t Read(UShort_t &row, UChar_t &pad, UShort_t &time, UShort_t &charge);
   Bool_t ReadSequence(UShort_t &row, UChar_t &pad, UShort_t &time, UChar_t &n, UShort_t **charges);
   void WriteFinal();

   Bool_t SetASCIIOutput(FILE *file);
   Bool_t SetBinaryOutput(FILE *file);
   Bool_t SetASCIIInput(FILE *file);
   Bool_t SetBinaryInput(FILE *file);

  private:

   UShort_t altromem[ALTRO_SIZE];
   UShort_t times_per_pad[1024];
   UShort_t charges_per_pad[1024];

   FILE *fInBinary;  //!
   FILE *fOutBinary; //!
   UShort_t lrow;
   UChar_t lpad;
   UShort_t rrow;  //read row
   UChar_t rpad;   //read pad
   UShort_t rtime; //read time
   UShort_t counter;  //total counter
   UShort_t tcounter; //time counter
   UShort_t pcounter; //packet counter
   UShort_t rcounter; //read counter
   UShort_t scounter; //sequence counter
   Bool_t flag; //Binary File?

   void Clear();
   void ClearRead();
   void Add(UShort_t charge, UShort_t time);
   void MakeAltroPackets();
   void AddTrailer();
   void Write();

   ClassDef(AliL3AltroMemHandler,1)
};

#endif






