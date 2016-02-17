#include "AliITSURecoSens.h"
#include "AliITSUGeomTGeo.h"
#include "AliExternalTrackParam.h"
#include "AliITSMFTAux.h"
#include "AliLog.h"

using namespace AliITSMFTAux;
using namespace TMath;

ClassImp(AliITSURecoSens)

//______________________________________________________
AliITSURecoSens::AliITSURecoSens(Int_t id)
:  fNClusters(0)
  ,fFirstClusterId(-1)
  ,fXTF(0)
  ,fPhiTF(0)
  ,fPhiMin(0)
  ,fPhiMax(0)
  ,fZMin(0)
  ,fZMax(0)
{
  // def. c-tor
  SetID(id);
}

//______________________________________________________
AliITSURecoSens::AliITSURecoSens(const AliITSURecoSens &source)
  :TObject(source)
  ,fNClusters(source.fNClusters)
  ,fFirstClusterId(source.fFirstClusterId)
  ,fXTF(source.fXTF)
  ,fPhiTF(source.fPhiTF)
  ,fPhiMin(source.fPhiMin)
  ,fPhiMax(source.fPhiMax)
  ,fZMin(source.fZMin)
  ,fZMax(source.fZMax)
{
  // copy c-tor
}

//______________________________________________________
AliITSURecoSens& AliITSURecoSens::operator=(const AliITSURecoSens &source)
{
  // = operator
  if (&source==this) return *this;
  TObject::operator=(source);
  fNClusters = source.fNClusters;
  fFirstClusterId = source.fFirstClusterId;
  fXTF = source.fXTF;
  fPhiTF = source.fPhiTF;
  fPhiMin = source.fPhiMin;
  fPhiMax = source.fPhiMax;
  fZMin   = source.fZMin;
  fZMax   = source.fZMax;
  //
  return *this;
}

//______________________________________________________
void AliITSURecoSens::SetBoundaries(double phiMn,double phiMx, double zMn, double zMx)
{
  // set phi,z limits 
  fPhiMin = phiMn;
  fPhiMax = phiMx;
  fZMin = zMn;
  fZMax = zMx;
}

//______________________________________________________
void AliITSURecoSens::Print(Option_t*) const			      
{
  //print 
  printf("Sensor%4d xTF=%+.3e phiTF=%+.3e | Phi:[%5.3f:%5.3f] Z:[%+7.3f:%+7.3f]\n",
	 GetID(),GetXTF(),GetPhiTF(), fPhiMin,fPhiMax, fZMin,fZMax);
}

//______________________________________________________
void AliITSURecoSens::ResetClusters()
{
  // discard old clusters
  fNClusters = 0;
  fFirstClusterId = -1;
}

//______________________________________________________
void AliITSURecoSens::ProcessClusters(Int_t)
{
  // create structures for fast finding
  //
  // to do
}

//______________________________________________________________________________
Int_t AliITSURecoSens::Compare(const TObject* obj)  const
{
  // compare sensor positions
  AliITSURecoSens* copy = (AliITSURecoSens*)obj;
  double phi  = MeanPhiSmall(fPhiMin,fPhiMax);
  double phiC = MeanPhiSmall(copy->fPhiMin,copy->fPhiMax);
  double span = DeltaPhiSmall(fPhiMin,fPhiMax)/2;
  double dPhi = DeltaPhiSmall(phi,phiC);
  //
  // special case to well define 1st raw (closest to 0 from above): wrap around 0/2pi
  if (dPhi>span) return phi<phiC ? -1 : 1;
  //
  double phiT = phi+span;
  BringTo02Pi(phiT);
  //  if (phiT<phiC && OKforPhiMin(phiT,phiC)) return -1;
  if (OKforPhiMin(phiT,phiC)) return -1;
  phiT = phi-span;
  BringTo02Pi(phiT);
  //if (phiT>phiC && OKforPhiMax( phiT,phiC)) return 1;
  if (OKforPhiMax( phiT,phiC)) return 1;
  //
  // sane phi range, check Z
  double dz = (fZMax-fZMin)/2;
  if (fZMax+dz < copy->fZMax) return -1;
  if (fZMin-dz > copy->fZMin) return 1;
  AliError(Form("Same chip compared? %d %d",GetID(),copy->GetID()));
  Print();
  copy->Print();
  return 0;
  //
}

//______________________________________________________________________________
Int_t AliITSURecoSens::CheckCoverage(double phi, double z) const
{
  // check if the sensor contains the impact point phi, z
  // if not, tell in which direction to move. 
  // kLeft, kRight are for smaller/larger angles, kUp,kDown for larger/smaller Z
  //
  int res = 0;
  if      (z<fZMin) res |= kDown;
  else if (z>fZMax) res |= kUp;
  //
  if      (!OKforPhiMin(fPhiMin,phi)) res |= kLeft;
  else if (!OKforPhiMax(fPhiMax,phi)) res |= kRight;
  return res;
}
