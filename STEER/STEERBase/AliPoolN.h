#ifndef ALIPOOL_N
#define ALIPOOL_N

#include <TObject.h>
#include <TArrayI.h>

class AliPoolN : public TObject
{
 public:
  //
 AliPoolN() : fSize(0),fPool(0),fPtrSlot(),fNBooked(0),fNFreed(0),fFreeSlot(0) {}
  AliPoolN(Int_t nini);
  AliPoolN(const AliPoolN& src);
  AliPoolN& operator=(const AliPoolN& src);
  //
  virtual  ~AliPoolN() {delete[] fPool;}
  virtual void Clear(Option_t* option="");
  void    PrintSummary(Option_t* option="") const;
  void    Reset()      {Clear();}
  //
  Int_t*    BookI( int n, Int_t         **addr) { return      (Int_t*)BookSlots((void**)addr,n,sizeof(Int_t));}
  UInt_t*   BookUI(int n, UInt_t        **addr) { return     (UInt_t*)BookSlots((void**)addr,n,sizeof(UInt_t));}
  Short_t*  BookS( int n, Short_t       **addr) { return    (Short_t*)BookSlots((void**)addr,n,sizeof(Short_t));}
  UShort_t* BookUS(int n, UShort_t      **addr) { return   (UShort_t*)BookSlots((void**)addr,n,sizeof(UShort_t));}
  Char_t*   BookC( int n, Char_t        **addr) { return     (Char_t*)BookSlots((void**)addr,n,sizeof(Char_t));}
  Bool_t*   BookB( int n, Bool_t        **addr) { return     (Bool_t*)BookSlots((void**)addr,n,sizeof(Bool_t));}
  Float_t*  BookF( int n, Float_t       **addr) { return    (Float_t*)BookSlots((void**)addr,n,sizeof(Float_t));}
  Double_t* BookD( int n, Double_t      **addr) { return   (Double_t*)BookSlots((void**)addr,n,sizeof(Double_t));}
  Long64_t* BookL64(int n, Long64_t     **addr) { return   (Long64_t*)BookSlots((void**)addr,n,sizeof(Long64_t));}
  Double32_t* BookD32(int n, Double32_t **addr) { return (Double32_t*)BookSlots((void**)addr,n,sizeof(Double32_t));}
  //
  void      FreeSlot(void* adr);/* {char** adrc = (char**)adr; 
    if (*adrc) {
      printf("ATTENTION: 0 pointed is freed\n");
    }
    *adrc = NULL; fNFreed++;}*/
  //
  Bool_t    IsReset() const {return GetUniqueID()==0;}
  Int_t*    GetArr() {return (Int_t*)fPool;}
  UInt_t    GetSize() {return fSize;}
  Int_t     GetUsedSize();
  void      Defragment();
  //
 protected:
  int*      BookSlots(void** addr, int n, int wsize);
  //
 protected:
  UInt_t  fSize;       // pool length
  Int_t*  fPool;       // container for values
  TArrayI fPtrSlot;    // container for object's addresses slots
  UInt_t  fNBooked;    // number of booked arrays
  UInt_t  fNFreed;     // number of freed arrays
  UInt_t  fFreeSlot;   // first free slot
  //
  ClassDef(AliPoolN,1) // pool for numerical values
};

#endif
