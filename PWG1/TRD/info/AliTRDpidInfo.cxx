#include "AliTRDpidInfo.h"
#include "AliTRDgeometry.h"

#include "string.h"

ClassImp(AliTRDpidInfo)
ClassImp(AliTRDpidInfo::AliTRDpidData)

//________________________________________________________________________
AliTRDpidInfo::AliTRDpidData::AliTRDpidData() : fPLbin(0xff) 
{
  memset(fdEdx, 0, 8*sizeof(Float_t));
}

//________________________________________________________________________
AliTRDpidInfo::AliTRDpidInfo()
  :TObject()
  ,fPID(-1)
  ,fNtracklets(0)
  ,fData(0x0)
{
  // Constructor of data array
  fData = new AliTRDpidData[AliTRDgeometry::kNlayer];
}

//________________________________________________________________________
AliTRDpidInfo::AliTRDpidInfo(Int_t idx)
  :TObject()
  ,fPID(idx)
  ,fNtracklets(0)
  ,fData(0x0)
{
  // Constructor of data array
  fData = new AliTRDpidData[AliTRDgeometry::kNlayer];
}

//________________________________________________________________________
AliTRDpidInfo::~AliTRDpidInfo()
{
  // Destructor
  if(fData) delete [] fData;
}

//________________________________________________________________________
AliTRDpidInfo::AliTRDpidData const* AliTRDpidInfo::GetDataInLayer(Int_t ily) const
{
  if(!fData) return NULL;
  if(ily<0 || ily>=AliTRDgeometry::kNlayer) return NULL;

  AliTRDpidData *data(NULL);
  for(Int_t itrklt=fNtracklets; itrklt--;){
    if(!(data=&fData[itrklt])) return NULL;
    if(data->Layer()==ily) return data;
  }
  return NULL;
}

//________________________________________________________________________
void AliTRDpidInfo::PushBack(Int_t ly, Int_t p, const Float_t *dedx)
{
// Add PID data to the end of the array 
  fData[fNtracklets].fPLbin= (ly<<4) | (p&0xf);
  memcpy(fData[fNtracklets].fdEdx, dedx, 8*sizeof(Float_t));
  fNtracklets++;
}

//________________________________________________________________________
void AliTRDpidInfo::Reset()
{
// Reset content

  if(!fNtracklets) return;
  while(fNtracklets--){
    fData[fNtracklets].fPLbin = 0xff;
    memset(fData[fNtracklets].fdEdx, 0, 8*sizeof(Float_t));
  }
  fNtracklets=0;
  fPID=-1;
}

