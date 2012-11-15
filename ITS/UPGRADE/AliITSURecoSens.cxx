#include "AliITSURecoSens.h"
#include "AliITSUGeomTGeo.h"
#include "AliITSsegmentation.h"
#include "AliExternalTrackParam.h"

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
  for (int i=kNNeighbors;i--;) fNeighbors[i] = -1;
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
  for (int i=kNNeighbors;i--;) fNeighbors[i] = source.fNeighbors[i];
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
  for (int i=kNNeighbors;i--;) fNeighbors[i] = source.fNeighbors[i];
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
  printf("Sensor%4d xTF=%+.3e phiTF=%+.3e | Phi:[%5.3f:%5.3f] Z:[%+7.3f:%+7.3f]| Neighb.:",
	 GetID(),GetXTF(),GetPhiTF(), fPhiMin,fPhiMax, fZMin,fZMax);
  for (int i=0;i<kNNeighbors;i++) printf(" %4d",fNeighbors[i]); printf("\n");
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
