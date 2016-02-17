#if !defined(__CINT__) || defined(__MAKECINT__)
#include "PixConv.h"
#include <vector>
#endif

int PrintNative(std::vector<PixConv::HitsRecord_t> &hitrec);
int PrintColRow(std::vector<PixConv::HitsRecord_t> &hitrec);

void read(const char* inp, bool conv2ColRow=1)
{
  PixConv converter(inp);
  unsigned short chipID,linkID,cycleID;
  //
  std::vector<PixConv::HitsRecord_t> hitsRecord;
  //
  int res = 0;
  int nhitsTot = 0;
  int nCycles = 0, cycleOld=-1;
  while(1) {
    int nhitsChip = 0;
    int res = converter.ReadChipData(hitsRecord,chipID,linkID,cycleID);
    if (cycleOld!=int(cycleID)) {
      printf("Reading new cycle %d\n",cycleID);
      cycleOld = cycleID;
      nCycles++;
    }
    if (!res) continue; // empty
    if (res<0) break; // EOF of Error
    printf("\nCycle %d | chip#%d of link#%d\n",cycleID,chipID,linkID);
    if (conv2ColRow) nhitsTot += PrintColRow(hitsRecord);
    else             nhitsTot += PrintNative(hitsRecord);
  }
  //
  printf("Read %d hits in %d cycles\n",nhitsTot,nCycles);
  //
  if (!res || res==PixConv::kEOF) printf("EOF reached successfully\n");
  else                            printf("ERROR was produced: %d\n",res);
  //
}

//_______________________________________________________
int PrintNative(std::vector<PixConv::HitsRecord_t> &hitrec)
{
  // print in native region, double column, address format
  printf("%4s %4s %4s %s\n","Reg.","DCol","Pix","ExtraHits");
  int nrec = hitrec.size(), nhitsChip=0;
  for (int ih=0;ih<nrec;ih++) {
    PixConv::HitsRecord_t& hitr = hitrec[ih];
    printf("%4d %4d %4d 0x%02x",hitr.region,hitr.dcolumn,hitr.address,hitr.hitmap);
    nhitsChip++;
    if (hitr.hitmap) { // extra hits
      printf("=[ ");
      for (int ip=0;ip<PixConv::kHitMapSize;ip++) {
	if (hitr.hitmap&(0x1<<ip)) {
	  printf("%4d ",hitr.address+ip+1);
	  nhitsChip++;
	}
      }
      printf("]");
    }
    printf("\n");
  }
  printf("read %d hits in %d records\n",nhitsChip,nrec);
  return nhitsChip;
}

//_______________________________________________________
int PrintColRow(std::vector<PixConv::HitsRecord_t> &hitrec)
{
  // print in row/column format
  printf(" col/row  | [extra]\n");
  int nrec = hitrec.size(), nhitsChip=0;
  for (int ih=0;ih<nrec;ih++) {
    PixConv::HitsRecord_t& hitr = hitrec[ih];
    int row = hitr.address>>1;
    int cold = (int(hitr.region)*PixConv::kNDColInReg + int(hitr.dcolumn))<<1;
    int col = cold + ((row&0x1) ? 1-(hitr.address&0x1) : (hitr.address&0x1));
    printf("%04d/%03d ",col,row);
    nhitsChip++;
    if (hitr.hitmap) { // extra hits
      printf(" [ ");
      for (int ip=0;ip<PixConv::kHitMapSize;ip++) {
	if (hitr.hitmap&(0x1<<ip)) {
	  int addr = hitr.address + ip + 1;
	  int rowE = addr>>1;
	  int colE = cold + ((rowE&0x1) ? 1-(addr&0x1) : (addr&0x1));
	  printf("%04d/%03d ",colE,rowE);
	  nhitsChip++;
	}
      }
      printf("]");
    }
    printf("\n");
  }
  printf("read %d hits in %d records\n",nhitsChip,nrec);
  return nhitsChip;
}
