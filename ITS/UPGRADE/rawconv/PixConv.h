#ifndef PIXCONV_H
#define PIXCONV_H

#include <stdio.h>
#include <vector>

#ifndef _BIG_ENDIAN_
#define _BIG_ENDIAN_  // store data in big endian format
#endif

//#define _DEBUG_PIX_CONV_ // uncomment for debug mode

class PixConv
{
public:
  struct HitsRecord {       // single records for hits (i.e. DATASHORT or DATALONG)
  HitsRecord() : region(0),dcolumn(0),address(0),hitmap(0) {}
  HitsRecord(unsigned char r,unsigned char dc, unsigned short adr, unsigned char hmap) : 
    region(r),dcolumn(dc),address(adr),hitmap(hmap) {}
    unsigned char  region;  // region ID
    unsigned char  dcolumn; // double column ID
    unsigned short address; // address in double colimn
    unsigned char  hitmap;  // hitmap for extra hits
  };
  typedef struct HitsRecord HitsRecord_t;
  enum { // masks for expected input from the stream
    kExpLinkHeader     = 0x1<<0
    ,kExpLinkTrailer   = 0x1<<1
    ,kExpChipHeader    = 0x1<<2
    ,kExpChipTrailer   = 0x1<<3
    ,kExpChipEmpty     = 0x1<<4
    ,kExpRegion        = 0x1<<5
    ,kExpData          = 0x1<<6
  };
  enum {kNRows=512, kNCols=1024, kNRegions=32,kNDColInReg=kNCols/kNRegions/2,kHitMapSize=7};
  enum {
    // masks for records components
    kMaskEncoder = 0x3c00   // encoder (double column) ID takes 4 bit max (0:15)
    ,kMaskPixID  = 0x3ff    // pixel ID within encoder (double column) takes 10 bit max (0:1023)
    ,kMaskDColID = kMaskEncoder|kMaskPixID // mask for encoder + dcolumn combination
    ,kMaskRegion = 0x1f     // region ID takes 5 bits max (0:31)
    ,kMaskChipID = 0x0f     // chip id in module takes 4 bit max
    ,kMaskROFlags= 0x0f     // RO flags in chip header takes 4 bit max
    ,kMaskFrameStartData = 0xff // Frame start data takes 8 bit max
    ,kMaskReserved       = 0xff // mask for reserved byte
    ,kMaskHitMap = 0x7f     // mask for hit map: at most 7 hits in bits (0:6) 
    ,kMaskLinkID = 0xffff   // THIS IS IMPROVISATION
    ,kMaskCycleID= 0xffff   // THIS IS IMPROVISATION
    //
    // record sizes in bytes
    ,kSizeRegion     = 1    // size of region marker in bytes
    ,kSizeChipHeader = 2    // size of chip header in bytes
    ,kSizeChipEmpty  = 3    // size of empty chip header in bytes
    ,kSizeLinkData   = 5    // size of the link headed/trailer
    //
    // flags for data records
    ,kREGION     = 0xc0     // flag for region
    ,kCHIPHEADER = 0xa0     // flag for chip header
    ,kCHIPTRAILER= 0xb0     // flag for chip trailer
    ,kCHIPEMPTY  = 0xe0     // flag for empty chip
    ,kDATALONG   = 0x0000   // flag for DATALONG
    ,kDATASHORT  = 0x4000   // flag for DATASHORT
    //
    ,kLINKHEADER = 0x01      // THIS IS IMPROVISATION
    ,kLINKTRAILER= 0x81      // THIS IS IMPROVISATION
    //
    ,kError = -1          // flag for decoding error
    ,kEOF=-100             // flag for EOF in reading
  };
  PixConv(const char* inp = 0);
  virtual ~PixConv();
  // 
  // IO
  bool OpenOutput(const char* name);
  void OpenInput(const char* name);
  void CloseIO();
  bool FlushBuffer();
  int  LoadInBuffer(int chunk=1000);

  // building chip hit map
  void AddPixel(short row, short col);
  void ResetMap();
  //
  // converting hitmap to raw data
  int  ProcDoubleCol(short reg, short dcol);
  int  ProcRegion(short reg);
  int  ProcChip(short chipInModule, short framestartdata, short roflags);
  //
  // reading raw data
  int ReadChipData(std::vector<PixConv::HitsRecord_t> &hits 
		   ,unsigned short &chip   // updated on change, don't modify returned value
		   ,unsigned short &link   // updated on change, don't modify returned value
		   ,unsigned short &cycle);// updated on change, don't modify returned value
  
  int UnexpectedEOF(const char* message) const;
  //
  void Print() const;
  void Reset();
  //
  unsigned short MakeChipHeader(short chipId, short framestartdata);
  unsigned short MakeChipTrailer(short roflags, short reserved=0);
  unsigned int   MakeChipEmpty(short chipId, short framestartdata, short reserved=0);
  //  
  unsigned char  MakeRegion(short reg);
  unsigned short MakeDataShort(short encoder, short address);
  unsigned short MakeDataLong(short encoder, short address);
  //
  // THIS IS IMPROVISATION
  unsigned long long   MakeLinkHeader(unsigned short link, unsigned short cycle);
  unsigned long long   MakeLinkTrailer(unsigned short link, unsigned short cycle);
  //
protected:
  struct PixLink {
    PixLink() {col=0; nextInRow=0;}
    short col;
    PixLink* nextInRow;
  };
  typedef struct PixLink PixLink_t;
  //
  void ExpandBuffer(int add=1000);
  void AddToBuffer(unsigned short v);
  void AddToBuffer(unsigned char v);
  void AddToBuffer(unsigned int v, int nbytes);
  void AddToBuffer(unsigned long long v, int nbytes);
  void EraseInBuffer(int nbytes);
  //
  bool GetFromBuffer(unsigned char &v);
  bool GetFromBuffer(unsigned short &v);
  bool GetFromBuffer(unsigned int &v, int nbytes);
  bool GetFromBuffer(unsigned long long &v, int nbytes);
  void StepBackInBuffer();
  //
  short ReadNextDCol(short* dest, short& dcol, short& region, short& chip, short& evid);
  //
protected:
  FILE*                  fIOFile;   //! handler for output
  //
  // cluster map
  unsigned int           fNPixels;  //! number of pixels seen
  int                    fCurrRow;  //! last row processed
  PixLink_t*             fCurrLink; //! last link added
  std::vector<PixLink*>  fRows;     //! column of 1st links in a row
  std::vector<short>     fRowIDs;   //! row ID's
  std::vector<PixLink_t*> fLinkPool; //! pool of links
  //
  unsigned char*         fWrBuffer;  //! write buffer
  unsigned char*         fBufferPointer;  //! current pointer in reading 
  unsigned char*         fBufferEnd;  //! end of filled buffer + 1
  int                    fWrBufferSize; //! buffer size
  int                    fWrBufferFill; //! entries in the buffer
  //
  unsigned int           fExpectInp;    //! type of input expected by reader
  //
  ClassDef(PixConv,1)
};

//_____________________________________
inline unsigned long long PixConv::MakeLinkHeader(unsigned short link, unsigned short cycle)
{
  // prepare link header: 00000001<linkID[31:16]><cycle[15:0]>
  unsigned long long v = kLINKHEADER;
  v = (v<<32) | ((kMaskLinkID&link)<<16) | (cycle&kMaskCycleID);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeLinkHeader: Link:%d Cycle:%d -> 0x%llx\n",link,cycle,v);
#endif
  return v;
}

//_____________________________________
inline unsigned long long PixConv::MakeLinkTrailer(unsigned short link, unsigned short cycle)
{
  // prepare link trailer: 10000001<linkID[31:16]><cycle[15:0]>
  unsigned long long v = kLINKTRAILER;
  v = (v<<32) | ((kMaskLinkID&link)<<16) | (cycle&kMaskCycleID);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeLinkTrailer: Link:%d Cycle:%d -> 0x%llx\n",link,cycle,v);
#endif
  return v;
}

//_____________________________________
inline unsigned short PixConv::MakeChipHeader(short chipID, short framestartdata) {
  // prepare chip header: 1010<chip id[3:0]><frame start data[7:0]>
  unsigned short v = kCHIPHEADER | (kMaskChipID&chipID);
  v = (v<<8) | (framestartdata&kMaskFrameStartData);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeChipHeader: chip:%d framdata:%d -> 0x%x\n",chipID,framestartdata,v);
#endif
  return v;
}
  
//_____________________________________
inline unsigned short PixConv::MakeChipTrailer(short roflags, short reserved) {
  // prepare chip trailer: 1011<readout flags[3:0]><reserved[7:0]>
  unsigned short v = kCHIPTRAILER | (kMaskROFlags&roflags);
  v = (v<<8) | (reserved&kMaskReserved);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeChipTrailer: ROflags:%d framdata:%d -> 0x%x\n",roflags,reserved,v);
#endif
  return v;
}

//_____________________________________
inline unsigned PixConv::MakeChipEmpty(short chipID, short framestartdata, short reserved) {
  // prepare chip empty marker: 1110<chip id[3:0]><frame start data[7:0] ><reserved[7:0]>
  unsigned int v = kCHIPEMPTY | (kMaskChipID&chipID);
  v = (((v<<8) | (framestartdata&kMaskFrameStartData))<<8) | (reserved&kMaskReserved);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeChipEmpty: chip:%d framdata:%d -> 0x%x\n",chipID,framestartdata,v);
#endif
  return v;
}

//_____________________________________
inline unsigned char PixConv::MakeRegion(short reg) {
  // packs the address of region
  unsigned char v = kREGION|(reg&kMaskRegion);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeRegion: region:%d -> 0x%x\n",reg,v);
#endif
  return v;
}

//_____________________________________
inline unsigned short PixConv::MakeDataShort(short encoder, short address) {
  // packs the address for data short
  unsigned short v = kDATASHORT | (kMaskEncoder&(encoder<<10)) | (address&kMaskPixID);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeDataShort: DCol:%d address:%d -> 0x%x\n",encoder,address,v);
#endif
  return v;
}

//_____________________________________
inline unsigned short PixConv::MakeDataLong(short encoder, short address) {
  // packs the address for data short
  unsigned short v = kDATALONG | (kMaskEncoder&(encoder<<10)) | (address&kMaskPixID);
#ifdef _DEBUG_PIX_CONV_
  printf("MakeDataLong: DCol:%d address:%d -> 0x%x\n",encoder,address,v);
#endif
  return v;
}

//_____________________________________
inline void PixConv::AddToBuffer(unsigned char v)
{
  // add character value to buffer 
  if (fWrBufferFill>=fWrBufferSize-1) ExpandBuffer();
  //
#ifdef _DEBUG_PIX_CONV_
  printf("AddToBuffer:C 0x%x\n",v);
#endif  
  //
  fWrBuffer[fWrBufferFill++] = v;
}

//_____________________________________
inline void PixConv::AddToBuffer(unsigned short v)
{
  // add short value to buffer 
  if (fWrBufferFill>=fWrBufferSize-2) ExpandBuffer();
  //
#ifdef _DEBUG_PIX_CONV_
  printf("AddToBuffer:S 0x%x\n",v);
#endif
  //
#ifdef _BIG_ENDIAN_
  fWrBuffer[fWrBufferFill++] = (v>>8)&0xff;
  fWrBuffer[fWrBufferFill++] = v&0xff;
#else 
  unsigned short* bfs = reinterpret_cast<unsigned short*>(fWrBuffer+fWrBufferFill);
  *bfs = v;
  fWrBufferFill+=2;
#endif
}

//_____________________________________
inline void PixConv::AddToBuffer(unsigned int v, int nbytes)
{
  // add 1-4 bytes to buffer 
  if (fWrBufferFill>=fWrBufferSize-nbytes) ExpandBuffer();
#ifdef _DEBUG_PIX_CONV_
  printf("AddToBuffer:I%d 0x%x\n",nbytes,v);
#endif  
  //
#ifdef _BIG_ENDIAN_
  for (int ib=nbytes;ib--;) {
    fWrBuffer[fWrBufferFill+ib] = (unsigned char)(v&0xff);
    v >>= 8;
  }
  fWrBufferFill += nbytes;
#else
  for (int ib=0;ib<nbytes;ib++) {
    fWrBuffer[fWrBufferFill++] = (unsigned char)(v&0xff);
    v >>= 8;
  }
#endif
  //
}

//_____________________________________
inline void PixConv::AddToBuffer(unsigned long long v, int nbytes)
{
  // add 1-8 bytes to buffer 
  if (fWrBufferFill>=fWrBufferSize-nbytes) ExpandBuffer();
  //
#ifdef _DEBUG_PIX_CONV_
  printf("AddToBuffer:LL%d 0x%llx\n",nbytes,v);
#endif  
  //
#ifdef _BIG_ENDIAN_
  for (int ib=nbytes;ib--;) {
    fWrBuffer[fWrBufferFill+ib] = (unsigned char)(v&0xff);
    v >>= 8;
  }
  fWrBufferFill += nbytes;
#else
  for (int ib=0;ib<nbytes;ib++) {
    fWrBuffer[fWrBufferFill++] = (unsigned char)(v&0xff);
    v >>= 8;
  }
#endif
}

//_____________________________________
inline void PixConv::EraseInBuffer(int nbytes)
{
  // erase last nbytes of buffer
  fWrBufferFill -= nbytes;
#ifdef _DEBUG_PIX_CONV_
  printf("EraseInBuffer: %d\n",nbytes);
#endif
}

//_____________________________________
inline bool PixConv::GetFromBuffer(unsigned char &v)
{
  // read character value from buffer
  if (fBufferPointer>=fBufferEnd && !LoadInBuffer()) return false; // upload or finish
  v = *fBufferPointer++;
  return true;
}

//_____________________________________
inline bool PixConv::GetFromBuffer(unsigned short &v)
{
  // read short value from buffer
  if (fBufferPointer>=fBufferEnd-(sizeof(short)-1) && !LoadInBuffer()) return false; // upload or finish
#ifdef _BIG_ENDIAN_
  v = (*fBufferPointer++)<<8;
  v |= (*fBufferPointer++);
#else
  v = (*fBufferPointer++);
  v |= (*fBufferPointer++)<<8;
  //  v = * reinterpret_cast<unsigned short*>(fBufferPointer);
  //  fBufferPointer+=sizeof(short);
#endif
  return true;
}

//_____________________________________
inline bool PixConv::GetFromBuffer(unsigned int &v, int nbytes)
{
  // read nbytes characters to int from buffer (no check for nbytes>4)
  if (fBufferPointer>=fBufferEnd-(nbytes-1) && !LoadInBuffer()) return false; // upload or finish
  v = 0;
#ifdef _BIG_ENDIAN_
  for (int ib=nbytes;ib--;) v |= (unsigned int)((*fBufferPointer++)<<(8<<ib));
#else
  for (int ib=0;ib<nbytes;ib++) v |= (unsigned int)((*fBufferPointer++)<<(8<<ib));
#endif
  return true;
}

//_____________________________________
inline bool PixConv::GetFromBuffer(unsigned long long &v, int nbytes)
{
  // read nbytes characters to int from buffer (no check for nbytes>8)
  if (fBufferPointer>=fBufferEnd-(nbytes-1) && !LoadInBuffer()) return false; // upload or finish
  v = 0;
#ifdef _BIG_ENDIAN_
  for (int ib=nbytes;ib--;) v |= (unsigned int)((*fBufferPointer++)<<(8<<ib));
#else
  for (int ib=0;ib<nbytes;ib++) v |= (unsigned int)((*fBufferPointer++)<<(8<<ib));
#endif
  return true;
}

//_____________________________________
inline void PixConv::StepBackInBuffer()
{
  // step back by 1 byte
  fBufferPointer--;
}


#endif
