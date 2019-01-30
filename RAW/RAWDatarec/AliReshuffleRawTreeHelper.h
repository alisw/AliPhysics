#ifndef __ALIRESHUFFLERAWTREEHELPER_H__
#define __ALIRESHUFFLERAWTREEHELPER_H__

#include "TObject.h"
#include "TString.h"
#include "TTree.h"

class AliReshuffleRawTreeHelper : public TObject {

 public:
  AliReshuffleRawTreeHelper();
  AliReshuffleRawTreeHelper(Int_t size);
  ~AliReshuffleRawTreeHelper();

  static AliReshuffleRawTreeHelper * GetInstance() ;

  void SetTree(TTree *t) {

    Char_t  ChunkName[30] ; 
    Int_t   eventNb;
    ULong64_t globID;


    t->SetBranchAddress("ChunkName" , &ChunkName );
    t->SetBranchAddress("eventNb"   , &eventNb   );
    t->SetBranchAddress("globID"    , &globID    );

    Long64_t NentriestreeID = t->GetEntries();

    if(fSIZE != NentriestreeID)
      {
	if(fChunkName) delete [] fChunkName;
	if(fEventNb)   delete [] fEventNb;
	if(fGlobID)    delete [] fGlobID;

	fChunkName = new TString[fSIZE] ;
	fEventNb   = new Int_t[fSIZE];
	fGlobID    = new ULong64_t[fSIZE];
      }


    for(Long64_t i=0; i<NentriestreeID; i++)
      {
	t->GetEntry(i);

	fChunkName[i]= ChunkName; 
	fEventNb[i] = eventNb;
	fGlobID[i] = globID;
      }

    //	fTree = t;

  }

  Int_t      GetNoffsets()                const {return fNoffsets;}
  Int_t      GetOffset(Int_t i)           const {return fOffsets[i];}
  Int_t      GetSMGroupIdex(Int_t i)      const {return fSMgroupInd[i];}
  Int_t   GetASize()                   const {return fSIZE;}
  Int_t      GetEventNumber(Long64_t i)   const {return fEventNb[i];}
  ULong64_t  GetGlobID(Int_t i)           const {return fGlobID[i];}
  TString    GetChunkName(Int_t i)        const  {return fChunkName[i];}
  TString    GetChunkPath ()              const { return fChunkPath;}
  TString    GetChunkFullName(Long64_t i) const  {return fChunkPath+fChunkName[i];}
  Int_t      GetGlobIDIndex(ULong64_t globid);
  Int_t      GetGlobIDIndex(UInt_t period, UInt_t orbit, UInt_t bc )
  {
    ULong64_t glid = (ULong64_t)bc+ (ULong64_t)orbit*3564+ (ULong64_t)period*16777215*3564;
    return GetGlobIDIndex(glid);
  }



  void       SetChunkPath (TString s)  {fChunkPath = s;}
  void       SetGroupIndices(Int_t *o) {for(Int_t i=0; i<20; i++) fSMgroupInd[i] = o[i];}
  void       SetGroupOffsets(Int_t n, Int_t *o) { fNoffsets=n; for(Int_t i=0; i<n; i++) fOffsets[i]=o[i];}
  void       SetEvent(Int_t i, Int_t nev, ULong64_t globid, TString st)
  { fEventNb[i] = nev; fGlobID[i] = globid; fChunkName[i] = st;}
  void       SetSize(Int_t size)
  {
    if(fSIZE != size)
      {
	if(fChunkName) delete [] fChunkName;
	if(fEventNb)   delete [] fEventNb;
	if(fGlobID)    delete [] fGlobID;

	fChunkName = new TString[fSIZE] ;
	fEventNb   = new Int_t[fSIZE];
	fGlobID    = new ULong64_t[fSIZE];
      }

  }

  void      PrintInfo(Int_t first=0, Int_t last=-1) const;


  AliReshuffleRawTreeHelper(const AliReshuffleRawTreeHelper & hl) ; 
  AliReshuffleRawTreeHelper & operator = (const AliReshuffleRawTreeHelper & hl);




 private:

  static AliReshuffleRawTreeHelper * fgInstance ; // pointer to the unique instance of the class

  Int_t  fSIZE;// =500000;
  Int_t fNoffsets;    // n SM offsets
  Int_t fSMgroupInd[20];// SM offset group
  Int_t fOffsets[20];// SM offsets
  TString   fChunkPath;
  TString   *fChunkName;//[fSIZE] 
  Int_t     *fEventNb;//[fSIZE]
  ULong64_t *fGlobID;//[fSIZE]


  Int_t  BinarySearch(Int_t l, Int_t r, ULong64_t globid);

  ClassDef(AliReshuffleRawTreeHelper,1);   // tbd

};
#endif
