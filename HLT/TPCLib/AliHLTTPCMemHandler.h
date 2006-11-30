// @(#) $Id$
// Original: AliHLTMemHandler.h,v 1.30 2004/10/06 08:51:20 cvetan 
#ifndef ALIHLTTPC_MEMHANDLER_H
#define ALIHLTTPC_MEMHANDLER_H

//_____________________________________________________________
// AliHLTTPCMemHandler
//
// The HLT Binary File handler 
//
//  This class does all the memory I/O handling of HLT binary files.
//  
// Author: Uli Frankenfeld <mailto:franken@fi.uib.no>, 
//         Anders Vestbo <mailto:vestbo$fi.uib.no>, 
//         Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
// *-- Copyright &copy ALICE HLT Group 

class AliHLTTPCDigitData;
class AliHLTTPCSpacePointData;
class AliHLTTPCDigitRowData;
class AliHLTTPCTrackSegmentData;
class AliHLTTPCTrackArray;
class AliHLTTPCRandomPointData;
class AliHLTTPCRandomDigitData;

#ifdef use_newio
class AliRunLoader;
class AliRawEvent;
#endif
class AliTPCRawStream;

class AliHLTTPCMemHandler { 

 public:
  AliHLTTPCMemHandler();
  virtual ~AliHLTTPCMemHandler();
  AliHLTTPCMemHandler(const AliHLTTPCMemHandler& /*m*/){};
  AliHLTTPCMemHandler& operator=(const AliHLTTPCMemHandler& /*&m*/)
    {return (*this);}
   
  void Reset(){CloseBinaryInput();CloseBinaryOutput();Free();}  
  void Init(Int_t s,Int_t p, Int_t *r=0);

  Bool_t SetBinaryInput(char *name);
  Bool_t SetBinaryInput(FILE *file);
  void CloseBinaryInput();
  
  Bool_t SetBinaryOutput(char *name);
  Bool_t SetBinaryOutput(FILE *file);
  void CloseBinaryOutput();

  //Random cluster
  void SetRandomCluster(Int_t maxnumber);
  void SetRandomSeed(UInt_t seed){srand(seed);}
  void SetRandomSeed();

  void ResetRandom(){fNDigits = 0; fNUsed = 0;}
  void Generate(Int_t row);
  void SetNGenerate(Int_t number){(number>fNRandom)?fNGenerate=fNRandom:fNGenerate = number;}

  void SetROI(Float_t *eta,Int_t *slice);
  void ResetROI();

  //Digit IO
  Bool_t Memory2Binary(UInt_t nrow,AliHLTTPCDigitRowData *data);
  Bool_t Binary2Memory(UInt_t & nrow,AliHLTTPCDigitRowData *data, UInt_t& sz);
  Bool_t Binary2Memory(UInt_t & nrow,AliHLTTPCDigitRowData *data){UInt_t tmp;return Binary2Memory(nrow,data,tmp);};

  Int_t Memory2CompMemory(UInt_t nrow,AliHLTTPCDigitRowData *data,UInt_t *comp);
  Int_t CompMemory2Memory(UInt_t nrow,AliHLTTPCDigitRowData *data,UInt_t *comp,UInt_t& sz);
  Int_t CompMemory2Memory(UInt_t nrow,AliHLTTPCDigitRowData *data,UInt_t *comp) {UInt_t tmp;return CompMemory2Memory(nrow,data,comp,tmp);};
  Bool_t CompMemory2CompBinary(UInt_t nrow,UInt_t *comp, UInt_t size=0);
  Bool_t CompBinary2CompMemory(UInt_t & nrow,UInt_t *comp);

  virtual AliHLTTPCDigitRowData *CompBinary2Memory(UInt_t & nrow, UInt_t& sz);
  virtual AliHLTTPCDigitRowData *CompBinary2Memory(UInt_t & nrow){UInt_t tmp;return CompBinary2Memory(nrow,tmp);};
  virtual Bool_t Memory2CompBinary(UInt_t nrow,AliHLTTPCDigitRowData *data);
  
  UInt_t GetNRow(UInt_t *comp,UInt_t size);

  //Point IO
  Bool_t Memory2Binary(UInt_t npoint,AliHLTTPCSpacePointData *data);
  Bool_t Binary2Memory(UInt_t & npoint,AliHLTTPCSpacePointData *data, UInt_t& sz);
  Bool_t Binary2Memory(UInt_t & npoint,AliHLTTPCSpacePointData *data) {UInt_t tmp; return Binary2Memory(npoint,data,tmp);};
  Bool_t Transform(UInt_t npoint,AliHLTTPCSpacePointData *data,Int_t slice);
  static void UpdateRowPointer(AliHLTTPCDigitRowData *&tempPt);
  
  //Track IO
  Bool_t Memory2Binary(UInt_t ntrack,AliHLTTPCTrackSegmentData *data);
  Bool_t Binary2Memory(UInt_t & ntrack,AliHLTTPCTrackSegmentData *data);
  Bool_t TrackArray2Binary(AliHLTTPCTrackArray *array);
  Bool_t Binary2TrackArray(AliHLTTPCTrackArray *array);
  Bool_t TrackArray2Memory(UInt_t & ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array) const;
  Bool_t Memory2TrackArray(UInt_t ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array) const;
  Bool_t Memory2TrackArray(UInt_t ntrack,AliHLTTPCTrackSegmentData *data,AliHLTTPCTrackArray *array,Int_t slice) const;
    
  //Memory Allocation
  UInt_t GetAllocatedSize() const {return fSize;}  
  UInt_t GetFileSize();
  UInt_t GetMemorySize(UInt_t nrow,UInt_t *comp) const;
  UInt_t GetCompMemorySize(UInt_t nrow,AliHLTTPCDigitRowData *data) const;
  UInt_t GetRandomSize() const;

  Byte_t *Allocate(UInt_t size);
  Byte_t *Allocate();  // allocate size of Binary Input File
  Byte_t *Allocate(AliHLTTPCTrackArray *array);
  Byte_t *GetDataPointer(UInt_t &size) {size = fSize; return fPt;}
  FILE *GetFilePointer() {return fInBinary;}
  void   Free();
  
  //Getters:
  Int_t GetRowMin() const {return fRowMin;}
  Int_t GetRowMax() const {return fRowMax;}
  Int_t GetSlice() const {return fSlice;}
  Int_t GetPatch() const {return fPatch;}
  
  //virtual functions:
  virtual void FreeDigitsTree() {fDummy=0; return;}
  virtual Bool_t SetAliInput(char */*name*/){fDummy=0; return 0;}
#ifdef use_newio
  virtual Bool_t SetAliInput(AliRunLoader */*runloader*/){fDummy=0; return 0;}
#endif
  virtual void CloseAliInput(){fDummy=0; return;} 
  virtual Bool_t IsDigit(Int_t /*i*/=0){fDummy=0; return 0;}
  virtual Bool_t SetMCOutput(char */*name*/){fDummy=0; return 0;}
  virtual Bool_t SetMCOutput(FILE */*file*/){fDummy=0; return 0;}
  virtual void CloseMCOutput(){fDummy=0; return;}
  virtual Bool_t AliDigits2Binary(Int_t /*event*/=0,Bool_t /*altro*/=kFALSE){fDummy=0; return 0;}
  virtual Bool_t AliDigits2CompBinary(Int_t /*event*/=0,Bool_t /*altro*/=kFALSE){fDummy=0; return 0;}  
  virtual AliHLTTPCDigitRowData *AliDigits2Memory(UInt_t & /*nrow*/,Int_t /*event*/=0){fDummy=0; return 0;}
  virtual AliHLTTPCDigitRowData *AliAltroDigits2Memory(UInt_t & /*nrow*/,Int_t /*event*/=0,Bool_t /*eventmerge*/=kFALSE){fDummy=0; return 0;}
  virtual void AliDigits2RootFile(AliHLTTPCDigitRowData */*rowPt*/,Char_t */*new_digitsfile*/){fDummy=0; return;}
  virtual Bool_t AliPoints2Binary(Int_t /*eventn*/=0){fDummy=0; return 0;}
  virtual AliHLTTPCSpacePointData *AliPoints2Memory(UInt_t & /*npoint*/,Int_t /*eventn*/=0){fDummy=0; return 0;}

  //AliHLTTPCRawDataFileHandler
  virtual Bool_t SetRawInput(Char_t */*name*/){fDummy=0; return 0;}
  virtual Bool_t SetRawInput(ifstream */*file*/){fDummy=0; return 0;}
  virtual void CloseRawInput(){} 
  virtual Int_t ReadRawInput(){fDummy=0; return 0;}
  virtual Short_t** GetRawData(Int_t &/*channels*/, Int_t & /*timebins*/){fDummy=0; return 0;}

  virtual Bool_t SetRawOutput(Char_t */*name*/){fDummy=0; return 0;}
  virtual Bool_t SetRawOutput(ofstream */*file*/){fDummy=0; return 0;}
  virtual void CloseRawOutput(){} 
  virtual Bool_t SaveRawOutput(){fDummy=0; return 0;}

  virtual Bool_t SetMappingFile(Char_t */*name*/){fDummy=0; return 0;}
  virtual Bool_t SetMappingFile(FILE */*file*/){fDummy=0; return 0;}
  virtual void CloseMappingFile(){} 
  virtual Int_t ReadMappingFile(){fDummy=0; return 0;}
  
  virtual Bool_t SetRawPedestalsInput(Char_t */*name*/){fDummy=0; return 0;}
  virtual Bool_t SetRawPedestalsInput(ifstream */*file*/){fDummy=0; return 0;}
  virtual void CloseRawPedestalsInput(){} 
  virtual Int_t ReadRawPedestalsInput(){fDummy=0; return 0;}

  virtual AliHLTTPCDigitRowData* RawData2Memory(UInt_t &/*nrow*/,Int_t /*event*/=-1){fDummy=0; return 0;}
  virtual Bool_t RawData2CompMemory(Int_t /*event*/=-1){fDummy=0; return 0;}

  //AliHLTTPCDDLDataFileHandler
#ifdef use_newio
  virtual Bool_t SetReaderInput(AliRawEvent */*rawevent*/){fDummy=0; return 0;}
  virtual Bool_t SetReaderInput(Char_t */*name*/,Int_t /*event*/=0){fDummy=0; return 0;}
#else
  virtual Bool_t SetReaderInput(Char_t */*name*/,Bool_t /*add*/=kTRUE){fDummy=0; return 0;}
#endif
  virtual void CloseReaderInput(){};

  virtual AliHLTTPCDigitRowData* DDLData2Memory(UInt_t &/*nrow*/,Int_t /*event*/=-1){fDummy=0; return 0;}
  virtual Bool_t DDLData2CompBinary(Int_t /*event*/=-1){fDummy=0; return 0;}

  virtual AliTPCRawStream* GetTPCRawStream(){fDummy=0; return 0;}

 protected:
  Int_t fRowMin; //min row
  Int_t fRowMax; //max row
  Int_t fSlice;  //slice
  Int_t fPatch;  //patch

  Int_t fEtaMinTimeBin[159]; //for ROI in eta only
  Int_t fEtaMaxTimeBin[159]; //for ROI in eta only
  
  FILE *fInBinary;//!
  FILE *fOutBinary;//!

 private:
  
  Byte_t *fPt;//!
  UInt_t fSize; //size of allocated data structure

  Bool_t fIsRandom; //random data generated
  Int_t fNRandom;   //count random digits 
  Int_t fNGenerate; //count generated digits
  Int_t fNUsed;     //count used digits
  Int_t fNDigits;   //count digits from digitstree

  AliHLTTPCRandomDigitData **fDPt;//!
  AliHLTTPCRandomDigitData *fRandomDigits;//!

  Int_t fDummy; // to fool the virtual const problem 
                // of the coding conventions tool

  void Write(UInt_t *comp, UInt_t & index, UInt_t & subindex, UShort_t value) const;
  UShort_t Read(UInt_t *comp, UInt_t & index, UInt_t & subindex) const;
  UShort_t Test(UInt_t *comp, UInt_t index, UInt_t subindex) const; 
  
  void DigitizePoint(Int_t row,Int_t pad, Int_t time,Int_t charge);
  void QSort(AliHLTTPCRandomDigitData **a, Int_t first, Int_t last);
  Int_t ComparePoints(UInt_t row,UShort_t pad,UShort_t time) const ;
  Int_t CompareDigits(AliHLTTPCRandomDigitData *a,AliHLTTPCRandomDigitData *b) const;
  void AddData(AliHLTTPCDigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge) const;
  void AddRandom(AliHLTTPCDigitData *data,UInt_t & ndata);
  void MergeDataRandom(AliHLTTPCDigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge);
  void AddDataRandom(AliHLTTPCDigitData *data,UInt_t & ndata,
                      UInt_t row,UShort_t pad,UShort_t time,UShort_t charge);


  ClassDef(AliHLTTPCMemHandler,1) // Memory handler class
};
#endif
