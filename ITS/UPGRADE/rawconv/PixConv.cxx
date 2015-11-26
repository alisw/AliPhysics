#include "PixConv.h"

ClassImp(PixConv)

//_____________________________________
PixConv::PixConv(const char* inp) 
: fIOFile(0),
  fNPixels(0)
  ,fCurrRow(-1)
  ,fCurrLink(0)
  ,fRows()
  ,fRowIDs()
  ,fLinkPool()
  ,fWrBuffer(0)
  ,fBufferPointer(0)
  ,fBufferEnd(0)
  ,fWrBufferSize(0)
  ,fWrBufferFill(0)
  ,fExpectInp(kExpLinkHeader)
{
  if (inp && inp[0]!=0) OpenInput(inp);
}

//_____________________________________
PixConv::~PixConv()
{
  for (int i=fLinkPool.size();i--;) delete fLinkPool[i];
  CloseIO();
}

//_____________________________________
void PixConv::Print() const
{
  int nr = fRowIDs.size();
  for (int ir=0;ir<nr;ir++) {
    PixLink_t* link = fRows[ir];
    printf("r%3d |",fRowIDs[ir]);
    do printf(" %3d",link->col); while(link=link->nextInRow);
    printf("\n");
  }
}

//_____________________________________
void PixConv::Reset()
{
  // reset before processing next chip
  fCurrLink = 0;
  fNPixels = 0;
  fCurrRow = -1;
  fRows.clear();
  fRowIDs.clear();
}

//_____________________________________
void PixConv::AddPixel(short row, short col)
{
  // add pixed to compressed matrix
  // the data must be sorted in row/col, no check is done
  if (fLinkPool.size()<fNPixels+1) fLinkPool.push_back(new PixLink_t());
  PixLink_t* link = fLinkPool[fNPixels++];
  link->col = col;
  link->nextInRow = 0;
  if (row==fCurrRow) { // extend current row
    fCurrLink->nextInRow = link;
    fCurrLink = link;
  }
  else { // create new row
    fCurrLink = link;
    fCurrRow = row;
    fRows.push_back(fCurrLink);
    fRowIDs.push_back(row);
  }
}

//_____________________________________
int PixConv::ProcChip(short chipInModule, short framestartdata, short roflags)
{
  //Print();
  // process chip data
  int nfound = 0;
  // chip header
  AddToBuffer(MakeChipHeader(chipInModule,framestartdata));
  for (int ir=0;ir<kNRegions;ir++) nfound += ProcRegion(ir);
  if (nfound) AddToBuffer(MakeChipTrailer(roflags,0));
  else {
    EraseInBuffer(kSizeChipHeader);
    AddToBuffer(MakeChipEmpty(chipInModule,framestartdata,0),kSizeChipEmpty);
  }
  //
  ResetMap();
  return nfound;
}

//_____________________________________
int PixConv::ProcRegion(short reg)
{
  // process region (16 double columns)
  AddToBuffer(MakeRegion(reg));
  int nfound = 0;
  for (int idc=0;idc<kNDColInReg;idc++) nfound += ProcDoubleCol(reg,idc);
  if (!nfound) EraseInBuffer(kSizeRegion);
  return nfound;
}

//_____________________________________
int PixConv::ProcDoubleCol(short reg, short dcol)
{
  // process double column
  short hits[2*kNRows];
  int nHits=0,nData=0;
  //
  int nr = fRows.size();
  short col0=((reg*kNDColInReg+dcol)<<1), col1=col0+1; // 1st,2nd column of double column 
  int prevRow = -1;
  for (int ir=0; ir<nr; ir++) {
    PixLink_t* lnk0 = fRows[ir];
    if (!lnk0 || lnk0->col>col1) continue; // no pixels left on this row or higher column IDs  
    //
    // process fired pixels
    bool left=0,right=0;
    if (lnk0->col==col0) { // pixel in left column
      left = 1;
      lnk0 = lnk0->nextInRow; // unlink processed pixel
    }
    if (lnk0 && lnk0->col==col1) { // pixel in right column
      right = 1;
      lnk0 = lnk0->nextInRow; // unlink processed pixel
    }
    int rowID = fRowIDs[ir];
    int addr0 = rowID<<1;
    if (rowID&0x1) { // odd rows: right to left numbering
      if (right) hits[nHits++] = addr0;
      if (left)  hits[nHits++] = addr0+1;
    }
    else { // even rows: left to right numbering
      if (left)  hits[nHits++] = addr0;
      if (right) hits[nHits++] = addr0+1;
    }
    if (!lnk0) {fRows[ir] = 0; continue;} // this row is finished
    fRows[ir] = lnk0;  // link remaining hit pixels to row
  }
  //
  int ih = 0;
  while(ih<nHits) {
    short addrE,addrW = hits[ih++];  // address of the reference hit
    unsigned char mask=0,npx=0;
    short addrLim = addrW + kHitMapSize+1; // 1+address of furthest hit can be put in the map
    while (ih<nHits && (addrE=hits[ih])<addrLim) {
      mask |= 0x1<<(addrE-addrW-1);
      ih++;
    }
    if (mask) { // flag DATALONG
      AddToBuffer(MakeDataLong(dcol,addrW));
      AddToBuffer(mask);
    }
    else AddToBuffer(MakeDataShort(dcol,addrW));
    nData++;
  }
  //
  return nData;
}

//_____________________________________
void PixConv::ExpandBuffer(int add)
{
#ifdef _DEBUG_PIX_CONV_
  printf("ExpandBuffer: %d -> %d\n",fWrBufferSize,fWrBufferSize+add);
#endif
  unsigned char* bfcopy = new unsigned char[(fWrBufferSize+=add)];
  if (fWrBufferFill) memcpy(bfcopy,fWrBuffer,fWrBufferFill*sizeof(unsigned char));
  delete[] fWrBuffer;
  fWrBuffer = bfcopy;
}

//_____________________________________
bool PixConv::OpenOutput(const char* name)
{
  // open output for raw data
  printf("Opening raw data output file: %s\n",name);
  fIOFile = fopen(name, "wb");
  return fIOFile!=0;
}

//_____________________________________
void PixConv::OpenInput(const char* name)
{
  // open raw data input 
  printf("Opening raw data input file: %s\n",name);
  if (!(fIOFile=fopen(name, "rb"))) {
    printf("Failed to open input file\n");
    exit(1);
  }
  fWrBufferFill = 0;
  fBufferPointer = fWrBuffer;
  fBufferEnd = fWrBuffer;
  fExpectInp = kExpLinkHeader;
  LoadInBuffer();
}

//_____________________________________
void PixConv::CloseIO()
{
  // close io
  if (fIOFile) fclose(fIOFile);
}

//_____________________________________
bool PixConv::FlushBuffer()
{
  // flush current content of buffer
  if (!fWrBufferFill) return false;
  if (!fIOFile) {
    printf("Output handler is not created");
    exit(1);
  }
  fwrite(fWrBuffer,sizeof(char),fWrBufferFill, fIOFile);
  fWrBufferFill = 0; // reset the counter
  return true;
}

//_____________________________________
int PixConv::LoadInBuffer(int chunk)
{
  // upload next chunk of data to buffer
  int save = fBufferEnd - fBufferPointer;
  // move unprocessed part to the beginning of the buffer
  if (save>0) memmove(fWrBuffer,fBufferPointer, save);
  else save = 0;
  //
  if (fWrBufferSize<(chunk+save)) ExpandBuffer(chunk+save+1000);
  fBufferPointer = fWrBuffer;
#ifdef _DEBUG_PIX_CONV_
  printf("LoadInBuffer: %d bytes placed with offset %d\n",chunk,save);
#endif
  int nc = (int)fread(fWrBuffer+save, sizeof(char), chunk, fIOFile);
  fWrBufferFill = save + nc;
  fBufferEnd = fWrBuffer + fWrBufferFill;
  return fWrBufferFill;
}

//_____________________________________
void PixConv::ResetMap()
{
  // reset map of hits for current chip
  fNPixels = 0;
  fCurrRow = -1;
  fCurrLink = 0;
  fRows.clear();
  fRowIDs.clear();
  fLinkPool.clear();
}

//_____________________________________
int PixConv::ReadChipData(std::vector<PixConv::HitsRecord_t> &hits 
			  ,unsigned short &chip   // updated on change, don't modify returned value
			  ,unsigned short &link   // updated on change, don't modify returned value
			  ,unsigned short &cycle  // updated on change, don't modify returned value
			  ) 
{
  // read record for single non-empty chip, updating on change link and cycle.
  // return number of records filled (>0), kEOF or kError  
  //
  hits.clear();
  unsigned char dataC=0;
  unsigned char region=0;
  unsigned char framestartdata=0;
  unsigned short dataS=0;
  unsigned int dataI=0;
  //
  int nHitsRec = 0;
  //
  while(1) {
    //
    if (!GetFromBuffer(dataC)) return kEOF;
    //
    if ( fExpectInp&kExpLinkHeader && dataC==kLINKHEADER) { // new link header
      if (!GetFromBuffer(link) || !GetFromBuffer(cycle)) return UnexpectedEOF("LINK_HEADER");
      fExpectInp = kExpLinkTrailer|kExpChipHeader|kExpChipEmpty;
      continue;
    }
    //
    if ( fExpectInp&kExpLinkTrailer && dataC==kLINKTRAILER) { // link trailer
      unsigned short linkT,cycleT;
      if (!GetFromBuffer(linkT) || !GetFromBuffer(cycleT)) return UnexpectedEOF("LINK_HEADER");
      if (linkT!=link || cycleT!=cycle) {
	printf("Error: expected link trailer for link%d/cycle%d, got for link%d/cycle%d\n",link,cycle,linkT,cycleT);
	return kError;
      }
      fExpectInp = kExpLinkHeader;
      continue;
    }
    // ---------- chip info ?
    unsigned char dataCM = dataC&(~kMaskChipID);
    //
    if ( (fExpectInp&kExpChipHeader) && dataCM==kCHIPHEADER) { // chip header was expected
      chip = dataC & kMaskChipID; 
      if (!GetFromBuffer(framestartdata)) return UnexpectedEOF("CHIP_HEADER");
      fExpectInp = kExpRegion;                                 // now expect region info
      continue;
    }
    //
    if ( (fExpectInp&kExpChipEmpty) && dataCM==kCHIPEMPTY) { // chip trailer was expected
      chip = dataC & kMaskChipID;
      if (!GetFromBuffer(framestartdata)) return UnexpectedEOF("CHIP_EMPTY:FrameStartData");
      if (!GetFromBuffer(dataC))          return UnexpectedEOF("CHIP_EMPTY:ReservedWord");
      fExpectInp = kExpLinkTrailer|kExpChipHeader|kExpChipEmpty;
      continue;
    }
    //
    if ( (fExpectInp&kExpChipTrailer) && dataCM==kCHIPTRAILER) { // chip trailer was expected
      if (!GetFromBuffer(framestartdata)) return UnexpectedEOF("CHIP_TRAILER:FrameStartData");
      fExpectInp = kExpLinkTrailer|kExpChipHeader|kExpChipEmpty;
      return hits.size();
    }
    // region info ? 
    if ( (fExpectInp&kExpRegion) && (dataC&kREGION)==kREGION ) { // chip header was seen, or hit data read
      region = dataC & kMaskRegion;
      fExpectInp = kExpData;
      continue;
    }
    // hit info ? 
    if ( (fExpectInp&kExpData) ) { // region header was seen, expect data
      StepBackInBuffer(); // need to reinterpred as short
      if (!GetFromBuffer(dataS)) return UnexpectedEOF("CHIPDATA");
      unsigned short dataSM = dataS & (~kMaskDColID); // check hit data mask
      if (dataSM==kDATASHORT) { // single hit
	unsigned char  dColID = (dataS & kMaskEncoder)>>10;
	unsigned short pixID  = dataS & kMaskPixID;
	hits.push_back( HitsRecord(region,dColID,pixID,0) );
	fExpectInp = kExpData|kExpRegion|kExpChipTrailer;
	continue;
      }
      else if (dataSM==kDATALONG) { // multiple hits
	unsigned char  dColID = (dataS & kMaskEncoder)>>10;
	unsigned short pixID  = dataS & kMaskPixID;
	unsigned char  hitsPattern=0;
	if (!GetFromBuffer(hitsPattern)) return UnexpectedEOF("CHIP_DATA_LONG:Pattern");
	hits.push_back( HitsRecord(region,dColID,pixID,hitsPattern) );
	fExpectInp = kExpData|kExpRegion|kExpChipTrailer;
	continue;
      }
      else {
	printf("Is this an error?\n");
	return kError;
      }
    }
    
  }
}


//_____________________________________
int PixConv::UnexpectedEOF(const char* message) const
{
  // error message on unexpected EOF
  printf("Error: Unexpected EOF on %s\n",message);
  return kError;
}
