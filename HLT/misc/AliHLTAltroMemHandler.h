// @(#) $Id$

#ifndef ALIL3ALTROMEMHANDLER_H
#define ALIL3ALTROMEMHANDLER_H

#include "AliHLTRootTypes.h"

//Maximum size of Altro packet equals 1000 bit
#define ALTRO_PACKET_SIZE 125

//Maximum 10 bit data to be stored in one packet
#define MAX_VALS 94

//Maximum size of array to store whole pad
#define ALTRO_SIZE (100*ALTRO_PACKET_SIZE)

class AliHLTAltroMemHandler {

  public: 

   AliHLTAltroMemHandler();
   virtual ~AliHLTAltroMemHandler() {}
   void Write(UShort_t row, UChar_t pad, UShort_t time, UShort_t charge);
   Bool_t Read(UShort_t &row, UChar_t &pad, UShort_t &time, UShort_t &charge);
   Bool_t ReadSequence(UShort_t &row, UChar_t &pad, UShort_t &time, UChar_t &n, UShort_t **charges);
   void WriteFinal();

   Bool_t SetASCIIOutput(FILE *file);
   Bool_t SetBinaryOutput(FILE *file);
   Bool_t SetASCIIInput(FILE *file);
   Bool_t SetBinaryInput(FILE *file);

  private:

   UShort_t fAltroMem[ALTRO_SIZE]; // Altro memoru
   UShort_t fTimesPerPad[1024]; // time samples per pad
   UShort_t fChargesPerPad[1024]; // charges per pad

   FILE *fInBinary;  //! // binary input file
   FILE *fOutBinary; //! // binary output file
   UShort_t fLRow;  //read row
   UChar_t fLPad;   //read pad
   UShort_t fRRow;  //read row
   UChar_t fRPad;   //read pad
   UShort_t fRTime; //read time
   UShort_t fCounter;  //total counter
   UShort_t fTCounter; //time counter
   UShort_t fPCounter; //packet counter
   UShort_t fRCounter; //read counter
   UShort_t fSCounter; //sequence counter
   Bool_t fFlag; //flag to indicate no return at the end of packet (binary)

   void Clear();
   void ClearRead();
   void Add(UShort_t charge, UShort_t time);
   void MakeAltroPackets();
   void AddTrailer();
   void Write();

   ClassDef(AliHLTAltroMemHandler,1)
};

typedef AliHLTAltroMemHandler AliL3AltroMemHandler; // for backward compatibility

#endif






