/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
// Class for Compression and Decompression          //
//////////////////////////////////////////////////////

 
#ifndef AliTPCCOMPRESSION_H
#define AliTPCCOMPRESSION_H

#ifdef __CINT__
class fstream;
#else
#include "Riostream.h"
#endif

class AliTPCHNode;
class AliTPCHTable;

class AliTPCCompression:public TObject{
 public:
  AliTPCCompression();
  virtual ~AliTPCCompression(){;}
  AliTPCCompression(const AliTPCCompression &source); // copy constructor
  AliTPCCompression& operator=(const AliTPCCompression &source); // ass. op.
  Int_t  CompressDataOptTables(Int_t NumTable,const char* fSource,const char* fDest);
  //This method is used to compress an Altro file using a set of general table previously calculated  and
  //stored as a sequence of txt file. 
  Int_t  DecompressDataOptTables(Int_t NumTables,const char* fname,const char* fDest="SourceDecompressed.dat");
  //This method is used to decompress a file compressed using the CompressDataOptTable method
  //It expects a set of table used for compressing the file in the same directory of the compressed file
  Int_t  Decompress(AliTPCHNode *RootNode[],Int_t NumTables,char* PointBuffer,UInt_t BufferSize,UShort_t out[],UInt_t &dim);
  //This method is used to decompress data stored in a char* buffer
  Int_t  DecompressWithLUTs(AliTPCHNode *RootNode[],UInt_t *LUTDimension[],AliTPCHNode **LUTNode[],Int_t /*NumTables*/,char* PointBuffer,UInt_t BufferSize,UShort_t out[],UInt_t &dim);
  Int_t  FillTables(const char* fSource,AliTPCHTable* table[],Int_t NumTables);
  //This method is used to compute the frequencies of the symbols in the source file
  Int_t  CreateTables(const char* fSource,Int_t NumTables);
  //This method is used to create and store the tables 
  Int_t  CreateTableFormula(Double_t beta,UInt_t  M,Int_t dim,Int_t Type);
  //This method is used to create and store the Bunch length or Time Gap Table using a formula
  void   SetVerbose(Int_t val){fVerbose=val;}
  //This method is used to set up the verbose level
  //   0 ==> No output messages are displayed
  //   1 ==> Some output messages are displayed during the running
  //   2 ==> A complete output is displayed
  void   ReadAltroFormat(char* fileOut,char* fileIn)const;
  //This method is used to read an Altro file and generate a text file containing the same information
  //It's is useful for debugging
  Int_t  CreateTablesFromTxtFiles(Int_t NumTable);
  //This method creates a set of binary tables starting from a set of txt tables
  Int_t  DestroyTables(AliTPCHNode *RootNode[],Int_t NumTables);
  //This methods is use to deallocate the memory used by Huffman trees
  Int_t  CreateTreesFromFile(AliTPCHNode *RootNode[],Int_t NumTables);
  //This method is like the previous one but the tables are stored in binary files
  //It is used in the decompression phase (DecompressDataOptTables())
  Int_t  CreateLUTsFromTrees(AliTPCHNode *RootNode[],Int_t NumTables,UInt_t *LUTDimension[],AliTPCHNode **LUTNode[]);
 private:
  Int_t   StoreTables(AliTPCHTable* table[],const Int_t NumTable);
  //This method is used to store an array of tables in a sequence of binary files
  //Each file contains the Size of the table (number of words) and for each word contains the correspondent 
  //codeword and codelength
  Int_t   RetrieveTables(AliTPCHTable* table[],Int_t NumTable);
  //This method is used to retrieve an array of tables from a sequence of binary files created using 
  //the previous method 
  void    DeleteHuffmanTree(AliTPCHNode* node);
  //This method is used to delete an Huffman tree
  void    VisitHuffmanTree(AliTPCHNode* node);
  //This method realizes an in order visit of a binary tree
  void    NextTable(Int_t Val,Int_t &NextTableType,Int_t &BunchLen,Int_t &Count)const;
  //This method is used to deduce which is the next table that as to be used to interpret the next value
  //reading the Altro format
  void    StoreValue(UInt_t val,UChar_t len);
  //This method is used to store a value in the compressed file 
  UInt_t  Mirror(UInt_t val,UChar_t len)const;
  //This method is used to get the specular value of a given value
  //for instance the specular value of 12345 is 54321
  void    Flush();
  //This method is used to complete and store the buffer in the output file when it isn't completely full 
  UInt_t  ReadWord(UInt_t NumberOfBit);
  //this method is used to read a specified number of bits from the compressed file
  UInt_t  ReadWordBuffer(UInt_t NumberOfBit);
  //this method is used to read a specified number of bits from the compressed memory buffer
  inline void   AdjustWordBuffer(UInt_t NumberOfBit);
  inline UInt_t ReadWordBufferWithLUTs();
  inline UInt_t ReadBitFromWordBuffer();
  inline UInt_t ReadBitFromLUTBuffer(UInt_t *buffer);
  void    ReadTrailer(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumberr,Bool_t Memory);
  //This method is used to read the trailer 
  inline UInt_t GetDecodedWordBuffer(AliTPCHNode* root);
  //This method is used to get a decoded word from the compressed file
  inline UInt_t GetDecodedWord(AliTPCHNode* root);
  //This method is used to get a decoded word from the compressed file
  inline UInt_t GetDecodedWordBufferWithLUTs(UInt_t* LUTDimension,AliTPCHNode** LUTNode);
  inline UInt_t GetDecodedLUTBuffer(AliTPCHNode** node,UInt_t *buffer);

  fstream f;                  // f is the logical name for the compressed and uncompressed file
  ofstream fStat;             // Logical name for the Statistics file
  UInt_t  fBuffer;            // buffer 
  Int_t   fDimBuffer;         // buffer dimension (32 bit)
  Int_t   fFreeBitsBuffer;    // number of free bits inside the buffer
  UInt_t  fBitsToRead;        // number of bits to be read
  UInt_t  fPos;               // current file position
  Int_t   fVerbose;           // verbose level (0 silent, !=0 output messages)
  UInt_t  fFillWords;         // Number of hexadecimally words (2AA pattern) inside a pad data block 
  UInt_t* fPointBuffer;       //pointer to the compressed raw data
  static const UInt_t   fgkTableDimension = 7; // size of LUTs for fast decompression

  ClassDef(AliTPCCompression,1)
};
#endif
