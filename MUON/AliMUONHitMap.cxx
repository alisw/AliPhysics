#include "AliMUONHitMap.h"
ClassImp(AliMUONHitMap)
ClassImp(AliMUONHitMapA1)

AliMUONHitMapA1::AliMUONHitMapA1(AliMUONsegmentation *seg, TObjArray *dig)
{
    fSegmentation = seg;
    fNpx  = fSegmentation->Npx();
    fNpy  = fSegmentation->Npy();
    fMaxIndex=2*(fNpx+1)*2*(fNpy+1)+2*fNpy;
    
    fHitMap = new Int_t[fMaxIndex];
    fDigits =  dig;
    fNdigits = fDigits->GetEntriesFast();
    Clear();
}


AliMUONHitMapA1::~AliMUONHitMapA1()
{
//    if (fDigits) delete   fDigits;
    if (fHitMap) delete[] fHitMap;
}

void AliMUONHitMapA1::Clear()
{
    memset(fHitMap,0,sizeof(int)*fMaxIndex);
}

Int_t AliMUONHitMapA1::CheckedIndex(Int_t ix, Int_t iy)
{
    Int_t index=2*fNpy*(ix+fNpx)+(iy+fNpy);
    if (index > fMaxIndex) {
	printf("\n \n \n Try to read/write outside array !!!! \n \n %d %d %d %d %d %d",ix,iy, fMaxIndex, index, fNpx, fNpy);
	return  fMaxIndex-1;
    } else {
	return index;
    }
}

	
void  AliMUONHitMapA1::FillHits()
{
    Int_t ndigits = fDigits->GetEntriesFast();
    //printf("\n Filling hits into HitMap\n");
    //printf("FindRawClusters -- ndigits %d \n",ndigits);
    if (!ndigits) return;
    AliMUONdigit *dig;
    for (Int_t ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliMUONdigit*)fDigits->UncheckedAt(ndig);
	SetHit(dig->fPadX,dig->fPadY,ndig);
    }
}


void  AliMUONHitMapA1::SetHit(Int_t ix, Int_t iy, Int_t idigit)
{
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=idigit+1;
    fHitMap[CheckedIndex(ix, iy)]=idigit+1;
}

void AliMUONHitMapA1::DeleteHit(Int_t ix, Int_t iy)
{
//    fHitMap[kMaxNpady*(ix+fNpx)+(iy+fNpy)]=0;
    fHitMap[CheckedIndex(ix, iy)]=0;
}

void AliMUONHitMapA1::FlagHit(Int_t ix, Int_t iy)
{
    fHitMap[CheckedIndex(ix, iy)]=
	-TMath::Abs(fHitMap[CheckedIndex(ix, iy)]);
}

Int_t AliMUONHitMapA1::GetHitIndex(Int_t ix, Int_t iy)
{
    return TMath::Abs(fHitMap[CheckedIndex(ix, iy)])-1;
}

TObject* AliMUONHitMapA1::GetHit(Int_t ix, Int_t iy)
{
    Int_t index=GetHitIndex(ix,iy);
    // Force crash if index does not exist ! (Manu)
    return (index <0) ? 0 : fDigits->UncheckedAt(GetHitIndex(ix,iy));
}

Flag_t AliMUONHitMapA1::TestHit(Int_t ix, Int_t iy)
{
    Int_t inf=fHitMap[CheckedIndex(ix, iy)];
    if (inf < 0) {
	return used;
    } else if (inf == 0) {
	return empty;
    } else {
	return unused;
    }
}


