/* Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////
// Utility Class for Compression and Decompression  //
//////////////////////////////////////////////////////

 
#ifndef AliTPCCOMPRESSION_H
#define AliTPCCOMPRESSION_H

class AliTPCHNode;
class AliTPCHTable;

class AliTPCCompression:public TObject{
 public:
  AliTPCCompression();
  virtual ~AliTPCCompression(){;}
  AliTPCCompression(const AliTPCCompression &source); // copy constructor
  AliTPCCompression& operator=(const AliTPCCompression &source); // ass. op.
  //This method is used to compress the data store in the altro format file using specific tables
  //calculate on a particular file that has to be compressed
  //The tables are stored at the beginning of the compressed file
  Int_t  CompressData(AliTPCHTable* table[],Int_t NumTable,const char* fSource,const char* fDest);
  //This methos is used to compress an Altro file using a set of general table previously calculated  and
  //stored as a sequence of txt file. In this case the tables are not stored in the compressed file
  Int_t  CompressDataOptTables(Int_t NumTable,const char* fSource,const char* fDest);
  //This method is used tho decompress a file compressed using the CompressData method
  Int_t  DecompressData(Int_t NumTables,const char* fname,char* fDest="SourceDecompressed.dat");
  //This methos is used yo decompress a file compressed using the CompressDataOptTable method
  //It expects a set of table used for compressing the file in the same direcotory of the compressed file
  Int_t  DecompressDataOptTables(Int_t NumTables,const char* fname,char* fDest="SourceDecompressed.dat");
  //This method is used to compute the frequencies of the symbols in the source file
  Int_t  FillTables(const char* fSource,AliTPCHTable* table[],const Int_t NumTables);
  //This method is used to create and store the tables 
  Int_t  CreateTables(const char* fSource,const Int_t NumTables);
  //This method is used to set up the verbose level
  //   0 ==> No output messages are displayed
  //   1 ==> Some output messages are displayed during the running
  //   2 ==> A complete output is displayed
  void   SetVerbose(Int_t val){fVerbose=val;}
  //This method is used to read an Altro file and generate a text file containing the same information
  //It's is useful for debugging
  void   ReadAltroFormat(char* fileOut,char* fileIn);
 private:
  //This method is used to store an array of tables in a sequence of binary files
  //Each file contains the Size of the table (number of words) and for each word contains the corrispondent 
  //codeword and codelength
  Int_t   StoreTables(AliTPCHTable* table[],const Int_t NumTable);
  //This method is used to retrieve an array of tables from a sequence of binaruy files created using 
  //the previous method 
  Int_t   RetrieveTables(AliTPCHTable* table[],Int_t NumTable);
  //This method is used to delete an Huffamn tree
  void    DeleteHuffmanTree(AliTPCHNode* node);
  //This method realizes an in order visit of a binary tree
  void    VisitHuffmanTree(AliTPCHNode* node);
  //This methos is used to create one or more Huffman tree strarting from one or more tables 
  //It is used in the decompression phase (DecompressData())
  void    CreateTrees(AliTPCHNode *RootNode[],const Int_t NumTables);
  //This method is like the previous one but the tables are stored in binary files
  //It is used in the decompression phase (DecompressDataOptTables())
  void    CreateTreesFromFile(AliTPCHNode *RootNode[],const Int_t NumTables);
  //This method is used to deduce which is the next table that as to be used to interpret the next value
  //reading the Altro format
  void    NextTable(Int_t Val,Int_t &NextTableType,Int_t &BunchLen,Int_t &Count);
  //This method is used to store a value in the compressed file 
  void    StoreValue(ULong_t val,UChar_t len);
  //This methos is used to get the specular value of a given value
  //for istance the specular value of 12345 is 54321
  ULong_t Mirror(ULong_t val,UChar_t len);
  //This method is used to complete and store the buffer in the output file when it isn't completely full 
  void    Flush();
  //this method is used to read a specified number of bits from the compressed file
  ULong_t ReadWord(Int_t NumberOfBit);
  //This method is used to read the trailer 
  void    ReadTrailer(Int_t &WordsNumber,Int_t &PadNumber,Int_t &RowNumber,Int_t &SecNumber);
  //This method is used to get a decoded word from the compressed file
  ULong_t GetDecodedWord(AliTPCHNode* root);

  fstream f;                  // f is the logical name for the compressed and uncompressed file
  ofstream stat;              // Statistics 
  ULong_t fBuffer;            // buffer 
  Int_t   fDimBuffer;         // buffer dimension (32 bit)
  Int_t   fFreeBitsBuffer;    // number of free bits inside the buffer
  Int_t   fReadBits;          // number of bit read
  ULong_t fPos;               // current file position
  Int_t   fVerbose;           // verbose level
  ULong_t fFillWords;
  ClassDef(AliTPCCompression,1)
};
#endif
