//    AliReshuffleRawTreeHelper.cxx
//
// map needed for reshuffling EMCAL LHC18r data
// Martin Poghosyan
// 03.01.2019
// v.1.3

#include "TMath.h"
#include "AliReshuffleRawTreeHelper.h"
#include "AliLog.h"


ClassImp(AliReshuffleRawTreeHelper);



AliReshuffleRawTreeHelper  * AliReshuffleRawTreeHelper::fgInstance = 0 ;


//____________________________________________________________________________
AliReshuffleRawTreeHelper::AliReshuffleRawTreeHelper():TObject(),
						       fSIZE(0),
						       fNoffsets(0),
						       fChunkPath(0),
						       fChunkName(0x0),
						       fEventNb(0x0),
						       fGlobID(0x0)
{
  for(Int_t i=0; i<20; i++) 
    {
      fOffsets[i]=0;
      fSMgroupInd[i]=0;
    }
  fgInstance=this;
  //printf("AliReshuffleRawTreeHelper::AliReshuffleRawTreeHelper() called\n");
}



//____________________________________________________________________________
AliReshuffleRawTreeHelper::AliReshuffleRawTreeHelper(Int_t size):TObject(),
								 fSIZE(size),
								 fNoffsets(0),
								 fChunkPath(0),
								 fChunkName(0x0),
								 fEventNb(0x0),
								 fGlobID(0x0)
{

  for(Int_t i=0; i<20; i++) 
    {
      fOffsets[i]=0;
      fSMgroupInd[i]=0;
    }

  fChunkName = new TString[fSIZE] ;
  fEventNb   = new Int_t[fSIZE];
  fGlobID    = new ULong64_t[fSIZE];

  fgInstance=this;
}


//____________________________________________________________________________
AliReshuffleRawTreeHelper::~AliReshuffleRawTreeHelper() 
{

  if(fChunkName) delete [] fChunkName;
  if(fEventNb)   delete [] fEventNb;
  if(fGlobID)    delete [] fGlobID;

  fSIZE=0;
  fChunkName=0x0;
  fEventNb=0x0;
  fGlobID=0x0;
  
}



//____________________________________________________________________________
AliReshuffleRawTreeHelper::AliReshuffleRawTreeHelper(const AliReshuffleRawTreeHelper & hl): TObject(),
											    fSIZE(0),
											    fChunkName(0x0),
											    fEventNb(0x0),
											    fGlobID(0x0)
{
  

  for(Int_t i=0; i<20; i++) 
    {
      fOffsets[i]=hl.fOffsets[i];
      fSMgroupInd[i]=hl.fSMgroupInd[i];
    }

  fNoffsets = hl.fNoffsets;
  fChunkPath = hl.fChunkPath;


  if(fChunkName) delete [] fChunkName;
  if(fEventNb)   delete [] fEventNb;
  if(fGlobID)    delete [] fGlobID;

  fSIZE =hl.fSIZE;
  fChunkName = new TString[fSIZE] ;
  fEventNb   = new Int_t[fSIZE];
  fGlobID    = new ULong64_t[fSIZE];


  for(Int_t i=0; i<fSIZE; i++)
    {
      fChunkName[i]=hl.fChunkName[i]; 
      fEventNb[i]=hl.fEventNb[i];
      fGlobID[i]=hl.fGlobID[i];
    }

  fgInstance = this;
}  



//____________________________________________________________________________
AliReshuffleRawTreeHelper & AliReshuffleRawTreeHelper::operator = (const AliReshuffleRawTreeHelper & hl)
{
  if(this!= &hl)
    /*
      return *this ;
      for(Int_t i=0; i<20; i++)
      fSMoffsets[i]=hl.fSMoffsets[i] ;  
      if(fChunkName) delete [] fChunkName;
      if(fEventNb)   delete [] fEventNb;
      if(fGlobID)    delete [] fGlobID;

      fSIZE =hl.fSIZE;
      fChunkName = new TString[fSIZE] ;
      fEventNb   = new Int_t[fSIZE];
      fGlobID    = new ULong64_t[fSIZE];


      for(Int_t i=0; i<fSIZE; i++)
      {
      fChunkName[i]=hl.fChunkName[i]; 
      fEventNb[i]=hl.fEventNb[i];
      fGlobID[i]=hl.fGlobID[i];
      }
    */
    AliError("Should not use operator= for singleton\n") ;

  return *this;
}



AliReshuffleRawTreeHelper * AliReshuffleRawTreeHelper::GetInstance()
{
  if(!fgInstance)
    fgInstance = new AliReshuffleRawTreeHelper() ;

  return fgInstance ;
}


Int_t  AliReshuffleRawTreeHelper::GetGlobIDIndex(ULong64_t globid)
{
  if(globid < fGlobID[0])
    return -1;

  if(globid > fGlobID[fSIZE-1 ])
    return -1;

  return BinarySearch(0, fSIZE-1, globid); 

}


Int_t AliReshuffleRawTreeHelper::BinarySearch(Int_t l, Int_t r, ULong64_t globid) 
{ 
  if (r>=l) 
    { 
      Int_t mid = l + (r - l)/2; 

      if ( fGlobID[mid] == globid) 
	return mid; 
      if (fGlobID[mid] > globid) 
	return BinarySearch(l, mid-1, globid); 
      return BinarySearch(mid+1, r, globid); 
    } 
  return -1; 
} 


void  AliReshuffleRawTreeHelper::PrintInfo(Int_t first, Int_t last) const
{

  printf("offsets ( %d groups)\n",fNoffsets);
  for(Int_t i=0; i<fNoffsets; i++)
    printf("group%d = %d,  ", i, fOffsets[i]);
  printf("\n");


  for(Int_t i=0; i<20; i++) printf(" SM%02d |",i);
  printf("\n");
  for(Int_t i=0; i<20; i++) printf(" %3d  |", fSMgroupInd[i]);
  printf("\n");

  printf("Chunk Path: %s\n",fChunkPath.Data());

  if(last> fSIZE) last=fSIZE-1;

  for(Int_t i=first; i<=last; i++)
    {
      printf("%d: %s  %d  %llu\n", i, fChunkName[i].Data(), fEventNb[i], fGlobID[i]);
    }

}

