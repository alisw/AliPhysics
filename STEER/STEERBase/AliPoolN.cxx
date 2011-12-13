#include "AliPoolN.h"
#include "AliLog.h"
#include <TMath.h>


ClassImp(AliPoolN)

//_______________________________________________
AliPoolN& AliPoolN::operator=(const AliPoolN& src)
{
  // assign
  if (this!=&src) {
    if (fPool) delete[] fPool;
    fSize = src.fSize;
    fPool = new Int_t[fSize];
    memcpy(fPool,src.fPool,sizeof(Int_t)*fSize);
    fPtrSlot = src.fPtrSlot;
    fFreeSlot = src.fFreeSlot;
    fNBooked  = src.fNBooked;
    fNFreed   = src.fNFreed;
  }
  return *this;
}


//_______________________________________________
void AliPoolN::Clear(Option_t*) 
{
  // reset 
  fNBooked = fFreeSlot = fNFreed = 0;
  SetUniqueID(0);
}

//_______________________________________________
int* AliPoolN::BookSlots(void** addr, int n, int wsize)
{
  // book n slots of size wsize (aligning to sizeof(int))
  const int kIntW = sizeof(Int_t);
  const int kDivW = TMath::Log2(sizeof(Int_t));
  const int kPtrW = sizeof(void*);
  const int kPtrIW = kPtrW/kIntW; // pointer length in int words
  static int eCounter = 0;
  //
  if (n<1) return 0;
  int nbtot = n*wsize; // total bytes to book (apart from addr)
  if (nbtot%kIntW) nbtot = (((nbtot>>kDivW)+1)<<kDivW); // make sure it is in Int words (normally 4 bites)
  int nw = kPtrIW + (nbtot>>kDivW); // number of int words to book (including space for pointer)
  UInt_t newFreeSlot = fFreeSlot+nw;
  //
  if (UInt_t(fPtrSlot.GetSize())<=fNBooked) fPtrSlot.Set( 2*(fNBooked+10) );
  int* ptrS = fPtrSlot.GetArray();
  if (fSize<=newFreeSlot) {   // do we need to expand the array ?
    //
    int sizeN = int((newFreeSlot+10000)*1.5);
    Int_t* oldArr = fPool;
    fPool = new Int_t[sizeN];
    if (fSize) memcpy(fPool,oldArr, fSize*sizeof(Int_t));
    memset(&fPool[fSize],0,(sizeN-fSize)*sizeof(Int_t));
    //
    for (UInt_t i=0;i<fNBooked;i++) {       // fix old addresses
      void* oldP  = (void*)&oldArr[ ptrS[i]+kPtrIW ];
      void** slotP = (void**)&fPool[ ptrS[i] ];
      void** ppp = (void**)slotP[0];
      //      printf("1>>%d %p %p old:%p\n",i,ppp, ppp ? *ppp : 0,oldP);
      if (!ppp || *ppp!=oldP) continue; // array was discarded      
      *ppp = &fPool[ ptrS[i]+kPtrIW ];
    }
    fSize = sizeN;
    delete[] oldArr;
    AliInfo(Form("Expansion %d: Size:%d, Object Booked:%d, Freed:%d, Next Free Slot:%d\n",
		 ++eCounter,fSize,fNBooked,fNFreed,fFreeSlot));
  }
  ptrS[fNBooked++] = fFreeSlot;
  void** slotP = (void**)&fPool[ fFreeSlot ];
  slotP[0] = addr;
  SetUniqueID(fFreeSlot+kPtrIW);
  fFreeSlot = newFreeSlot;
  //
  return &fPool[GetUniqueID()];
}

//_______________________________________________
void AliPoolN::PrintSummary(Option_t*) const
{
  // print summary
  printf("Pool: Num.Arrays: Size:%d, Object Booked:%d, Freed:%d, Next Free Slot:%d\n",
	 fSize,fNBooked,fNFreed,fFreeSlot);
}
 

//_______________________________________________
AliPoolN::AliPoolN(Int_t nini) : 
  fSize(nini>100 ? nini:10000),fPool(new Int_t[fSize]),fPtrSlot(fSize/10),fNBooked(0),
  fNFreed(0),fFreeSlot(0) 
{
  // c-tor
}

//_______________________________________________
AliPoolN::AliPoolN(const AliPoolN& src) : 
  TObject(src), fSize(src.fSize),fPool(new Int_t[fSize]),fPtrSlot(src.fPtrSlot), 
  fNBooked(src.fNBooked), fNFreed(src.fNFreed),fFreeSlot(src.fFreeSlot) 
{
  // c-tor
}

//_______________________________________________
void AliPoolN::FreeSlot(void* adr) 
{
  char** adrc = (char**)adr; 
  if (!(*adrc)) return; //{ printf("ATTENTION: 0 pointer is freed\n"); return; }
  *adrc = NULL; fNFreed++;
}

//_______________________________________________
Int_t AliPoolN::GetUsedSize()
{
  // estimate which fraction of array is actually used
  const int kIntW = sizeof(Int_t);
  const int kPtrW = sizeof(void*);
  const int kPtrIW = kPtrW/kIntW; // pointer length in int words
  //
  if (fFreeSlot<1) return 0.; // empty
  int nFreeW=0, freeStart=-1,freeEnd=-1;
  int* ptrS = fPtrSlot.GetArray();
  //
  for (UInt_t i=0;i<fNBooked;i++) {
    void* arrP  = (void*)&fPool[ ptrS[i]+kPtrIW ]; // address of the booked array
    void** slotP = (void**)&fPool[ ptrS[i] ];      // address of the external pointer to which this array was attached
    void** ppp = (void**)slotP[0];
    if (!ppp || *ppp!=arrP) { // the array was discarded
      if (freeStart<0) {freeStart = ptrS[i]; freeEnd = -1;}
    }
    else if (freeStart!=-1) {
      freeEnd = ptrS[i]; // end of free blocks
      nFreeW += freeEnd-freeStart;
      freeStart = -1;
    }
  }
  if (freeStart!=-1 && freeEnd<0) { // last blocks were freed
    freeEnd = fFreeSlot;
    nFreeW += freeEnd-freeStart;
  }
  return fFreeSlot - nFreeW;
  //
}

//_______________________________________________
void AliPoolN::Defragment()
{
  // deframent the pool
  const int kIntW = sizeof(Int_t);
  const int kPtrW = sizeof(void*);
  const int kPtrIW = kPtrW/kIntW; // pointer length in int words
  //
  if (fFreeSlot<1) return; // empty
  int freeStart=-1,useStart=-1,useFirst=-1;
  int* ptrS = fPtrSlot.GetArray();
  int cnt = 0, shift = 0;
  
  //
  for (UInt_t i=0;i<fNBooked;i++) {
    void* arrP  = (void*)&fPool[ ptrS[i]+kPtrIW ]; // address of the booked array
    void** slotP = (void**)&fPool[ ptrS[i] ];      // address of the external pointer to which this array was attached
    void** ppp = (void**)slotP[0];
    if (!ppp || *ppp!=arrP) { // the array was discarded
      if (useStart!=-1) { // and it is also the end of the used block
	if (freeStart!=-1) { // and there was a free block before
	  int moved = ptrS[i]-useStart;
	  shift = useStart-freeStart;
	  memmove(&fPool[freeStart], &fPool[useStart], kIntW*moved);
	  // update pointers for moved block
	  for (UInt_t j=useFirst;j<i;j++) { // used block covered these arrays
	    int slt = ptrS[cnt++] = ptrS[j]-shift; // updated index of the slot corresponding to array in the used block
	    slotP = (void**)&fPool[ slt ];
	    ppp = (void**)slotP[0];
	    *ppp = &fPool[ slt+kPtrIW ];           // assign external array address the new pointer on the begginning of the array
	  }
	  freeStart = ptrS[i]-shift; // mark beginning of the new free block
	  useStart = -1; // forget already treated used block 
	}
	else {
	  useStart = -1;       // forget this used block, it is already in defragmented part
	  freeStart = ptrS[i]; // no free block before, just mark beginning of new free block
	}
      }
      else if (freeStart==-1) freeStart = ptrS[i]; // previous slot was not used: mark beginning of new free block
    }
    else { // the array is used
      if (useStart==-1) { // and it is in the beginning of the used block
	useStart = ptrS[i];
	useFirst = i;	
      }
      if (freeStart==-1) cnt++; // no free block before, this array will not be moved
    }
  }
  //
  if (useStart!=-1 && freeStart!=-1) { // is there a last used block separated from defragmented part by empty block
    int moved = fFreeSlot-useStart;
    shift = useStart-freeStart;
    memmove(&fPool[freeStart], &fPool[useStart], kIntW*moved);
    // update pointers for moved block
    for (UInt_t j=useFirst;j<fNBooked;j++) { // used block covered these arrays
      int slt = ptrS[cnt++] = ptrS[j]-shift; // updated index of the slot corresponding to array in the used block
      void** slotP = (void**)&fPool[ slt ];
      void** ppp = (void**)slotP[0];       // address of the pointer to which the array was assigned
      *ppp = &fPool[ slt+kPtrIW ];         // assign to it the new begginning of the array
    }
  }
  fFreeSlot -= shift;
  fNBooked = cnt;
  fNFreed  = 0;
  AliInfo(Form("Freed %d words",shift));
  //
}
