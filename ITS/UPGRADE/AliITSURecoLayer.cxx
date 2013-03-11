#include <TClonesArray.h>
#include "AliITSURecoLayer.h"
#include "AliITSsegmentation.h"
#include "AliITSUAux.h"
#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliLog.h"

using namespace AliITSUAux;
using namespace TMath;

ClassImp(AliITSURecoLayer)


//______________________________________________________
AliITSURecoLayer::AliITSURecoLayer(const char* name)
  :fActiveID(-1)
  ,fNSensors(0)
  ,fNSensInLadder(0)
  ,fNLadders(0)
  ,fR(0)
  ,fRMax(0)
  ,fRMin(0)
  ,fZMax(0)
  ,fZMin(0)
  ,fPhiLadMax(0)
  ,fPhiLadMin(0)
  ,fPhiOffs(0)
  ,fSensDZInv(0)
  ,fDPhiLadInv(0)
  ,fMaxStep(0.5)
  ,fSensors(0)
  ,fITSGeom(0)
  ,fClusters(0)
{
  // def. c-tor
  SetNameTitle(name,name);
}

//______________________________________________________
AliITSURecoLayer::AliITSURecoLayer(const char* name, Int_t activeID, AliITSUGeomTGeo* gm)
  :fActiveID(activeID)
  ,fNSensors(0)
  ,fNSensInLadder(0)
  ,fNLadders(0)
  ,fR(0)
  ,fRMax(0)
  ,fRMin(0)
  ,fZMax(0)
  ,fZMin(0)
  ,fPhiLadMax(0)
  ,fPhiLadMin(0)
  ,fPhiOffs(0)
  ,fSensDZInv(0)
  ,fDPhiLadInv(0)
  ,fMaxStep(0.5)
  ,fSensors(0)
  ,fITSGeom(gm)
  ,fClusters(0)
{
  // def. c-tor
  SetNameTitle(name,name);
  Build();
}

//______________________________________________________
AliITSURecoLayer::~AliITSURecoLayer()
{
  // def. d-tor
  delete[] fSensors;
  delete[] fPhiLadMax;
  delete[] fPhiLadMin;
  if (GetOwnsClusterArray()) delete fClusters;
}

//______________________________________________________
void AliITSURecoLayer::Print(Option_t* opt) const			      
{
  //print 
  printf("Lr %-15s %d (act:%+d), NSens: %4d | MaxStep:%.2f ",GetName(),GetID(),GetActiveID(),GetNSensors(),fMaxStep);
  printf("%6.3f<R<%6.3f | %+8.3f<Z<%+8.3f dZ:%6.3f\n",fRMin,fRMax,fZMin,fZMax,fSensDZInv>0 ? 1/fSensDZInv : 0);
  TString opts = opt; opts.ToLower();
  if (opts.Contains("sn")) for (int i=0;i<GetNSensors();i++) GetSensor(i)->Print(opt);
}

//______________________________________________________
void AliITSURecoLayer::Build()
{
  // build internal structures
  if (fActiveID<0) return;
  fNLadders = fITSGeom->GetNLadders(fActiveID);
  fNSensInLadder = fITSGeom->GetNDetectors(fActiveID);
  fNSensors = fNLadders*fNSensInLadder;
  fSensors = new AliITSURecoSens*[fNSensors];
  const AliITSsegmentation* kSegm = fITSGeom->GetSegmentation(fActiveID);
  //
  // name layer according its active id, detector type and segmentation tyoe
  TGeoHMatrix mmod;
  const TGeoHMatrix* mt2l;
  fRMin=fZMin=1e9;
  fRMax=fZMax=-1e9;
  double phiTF,rTF, loc[3]={0,0,0},glo[3];
  fNSensors = 0;
  fPhiLadMin = new Double_t[fNLadders];
  fPhiLadMax = new Double_t[fNLadders];
  fSensDZInv = 0;
  fDPhiLadInv = fNLadders/TwoPi();
  //
  for (int ild=0;ild<fNLadders;ild++) {
    fPhiLadMin[ild] = 1e9;
    fPhiLadMax[ild] = -1e9;
    //
    for (int idt=0;idt<fNSensInLadder;idt++) {
      AliITSURecoSens* sens = new AliITSURecoSens(fNSensors++);
      fSensors[ild*fNSensInLadder+idt] = sens;
      //
      double phiMin=1e9,phiMax=-1e9,zMin=1e9,zMax=-1e9;
      mmod = *fITSGeom->GetMatrix(fActiveID,ild,idt);
      for (int ix=0;ix<2;ix++) {
	loc[0] = (ix-0.5)*kSegm->Dx(); // +-DX/2
	for (int iy=0;iy<2;iy++) {
	  loc[1] = (iy-0.5)*kSegm->Dy(); // +-DY/2
	  for (int iz=0;iz<2;iz++) {
	    loc[2] = (iz-0.5)*kSegm->Dz(); // +-DZ/2
	    //
	    mmod.LocalToMaster(loc,glo);
	    double phi = ATan2(glo[1],glo[0]);
	    double r   = glo[0]*glo[0] + glo[1]*glo[1];
	    if (fRMin>r) fRMin = r;
	    if (fRMax<r) fRMax = r;
	    BringTo02Pi(phi);
	    if      (phiMin>1e8) phiMin=phi; 
	    else if (!OKforPhiMin(phiMin,phi)) phiMin=phi;
	    if      (phiMax<-1e8) phiMax=phi;
	    else if (!OKforPhiMax(phiMax,phi)) phiMax=phi;	      
	    if (glo[2]>zMax) zMax=glo[2];
	    if (glo[2]<zMin) zMin=glo[2];
	  }
	}
      }
      sens->SetBoundaries(phiMin,phiMax,zMin,zMax);
      mt2l = fITSGeom->GetMatrixT2L(fActiveID,ild,idt);
      mmod.Multiply(mt2l);	
      loc[0]=loc[1]=loc[2]=0;
      mmod.LocalToMaster(loc,glo);
      rTF   = Sqrt(glo[0]*glo[0] + glo[1]*glo[1]);  //  tracking params (misaligned)
      phiTF = ATan2(glo[1],glo[0]);
      BringTo02Pi(phiTF);
      //
      sens->SetXTF(rTF);
      sens->SetPhiTF(phiTF);
      //
      if      (fPhiLadMin[ild]>1e8)  fPhiLadMin[ild] = phiMin;
      else if (!OKforPhiMin(fPhiLadMin[ild],phiMin)) fPhiLadMin[ild] = phiMin;
      if      (fPhiLadMax[ild]<-1e8) fPhiLadMax[ild] = phiMax;
      else if (!OKforPhiMax(fPhiLadMax[ild],phiMax)) fPhiLadMax[ild] = phiMax;
      if (fZMin>zMin) fZMin = zMin;
      if (fZMax<zMax) fZMax = zMax;
      //
      if (idt>0) fSensDZInv += zMax - GetSensor(ild,idt-1)->GetZMax(); // z interval to previous
    }
  }
  //
  fRMin = Sqrt(fRMin);
  fRMax = Sqrt(fRMax);
  fR = 0.5*(fRMin+fRMax);
  double dz = fNSensInLadder>0 ? fSensDZInv/(fNSensInLadder-1)/fNLadders : fZMax-fZMin;
  fSensDZInv = 1./dz;

  const int kNBId[3][3] = { 
    {AliITSURecoSens::kNghbBL,AliITSURecoSens::kNghbB,AliITSURecoSens::kNghbBR},
    {AliITSURecoSens::kNghbL,          -1            ,AliITSURecoSens::kNghbR },
    {AliITSURecoSens::kNghbTL,AliITSURecoSens::kNghbT,AliITSURecoSens::kNghbTR}
  };

  // add neighbours info
  double zTol = 0.45*dz, phiTol = 0.45*TwoPi()/fNLadders;
  for (int ild=0;ild<fNLadders;ild++) {
    for (int idt=0;idt<fNSensInLadder;idt++) {
      AliITSURecoSens* sens = GetSensor(ild,idt);
      //
      for (int ils=-1;ils<=1;ils++) {
	int ildN = ild+ils;  // ladders of neighbouring sensors
	if (ildN<0) ildN = fNLadders-1; else if (ildN==fNLadders) ildN = 0;
	for (int ids=-1;ids<=1;ids++) {
	  int idtN = idt+ids;
	  if (idtN<0 || idtN==fNSensInLadder || (ids==0&&ils==0)) continue;
	  AliITSURecoSens* sensN = GetSensor(ildN,idtN); // potential neighbor
	  int neighbID = ildN*fNSensInLadder+idtN;
	  //	  
	  int zType = 1;  // side
	  if (sens->GetZMin()-zTol  > sensN->GetZMax()) continue; // too large distance
	  if (sensN->GetZMin()-zTol > sens->GetZMax() ) continue; // too large distance
	  if      (sens->GetZMin()-zTol>sensN->GetZMin()) zType =  0;     // bottom
	  else if (sensN->GetZMin()-zTol>sens->GetZMin()) zType =  2;     // top
	  //
	  int phiType = 1;

	  double phiTstMn = sensN->GetPhiMin()-phiTol;
	  BringTo02Pi(phiTstMn);
	  if (!OKforPhiMax(sens->GetPhiMax(),phiTstMn)) continue; // too large angle	  
	  double phiTstMx = sensN->GetPhiMax()+phiTol;	  
	  BringTo02Pi(phiTstMx);
	  if (!OKforPhiMin(sens->GetPhiMin(),phiTstMx)) continue; // too large angle
	  //
	  phiTstMn = sensN->GetPhiMin()+phiTol;
	  BringTo02Pi(phiTstMn);
	  phiTstMx = sensN->GetPhiMax()-phiTol;	  
	  BringTo02Pi(phiTstMx);
	  if      (!OKforPhiMax(sens->GetPhiMax(),phiTstMx)) phiType = 0; // left
	  else if (!OKforPhiMin(sens->GetPhiMin(),phiTstMn)) phiType = 2; // right
	  //
	  sens->SetNeighborID(kNBId[zType][phiType], neighbID);
	} // phi scan
      } // z scan
    } // sensors
  } // ladders
  //
}

//______________________________________________________
Int_t AliITSURecoLayer::FindSensors(const double* impPar, AliITSURecoSens *sensors[AliITSURecoSens::kNNeighbors])
{
  // find sensors having intersection with track
  // impPar contains: lab phi of track, dphi, labZ, dz
  double z = impPar[2];
  if (z>fZMax+impPar[3]) return 0; // outside of Z coverage
  z -= fZMin;
  if (z<-impPar[3]) return 0; // outside of Z coverage
  int sensInLad = int(z*fSensDZInv);
  if      (sensInLad<0) sensInLad = 0;
  else if (sensInLad>=fNSensInLadder) sensInLad = fNSensInLadder-1;
  //
  double phi = impPar[0] - fPhiOffs;
  BringTo02Pi(phi);
  int ladID = int(phi*fDPhiLadInv);  // ladder id
  int nsens = 0;
  //
  AliITSURecoSens* sensN,*sens = GetSensor(ladID*fNSensInLadder+sensInLad);
  sensors[nsens++] = sens;
  // check neighbours
  double zMn=impPar[2]-impPar[3], zMx=impPar[2]+impPar[3], phiMn=impPar[0]-impPar[1], phiMx=impPar[0]+impPar[1];
  BringTo02Pi(phiMn);
  BringTo02Pi(phiMx);
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbR)); // neighbor on the right (smaller phi)
  if (sensN && OKforPhiMin(phiMn,sensN->GetPhiMax())) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbTR)); // neighbor on the top right (smaller phi, larger Z)
  if (sensN && OKforPhiMin(phiMn,sensN->GetPhiMax()) && sensN->GetZMin()<zMx) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbT)); // neighbor on the top (larger Z)
  if (sensN && sensN->GetZMin()<zMx) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbTL)); // neighbor on the top left (larger Z, larger phi)
  if (sensN && OKforPhiMax(phiMx,sensN->GetPhiMin()) && sensN->GetZMin()<zMx) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbL)); // neighbor on the left (larger phi)
  if (sensN && OKforPhiMax(phiMx,sensN->GetPhiMin())) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbBL)); // neighbor on the bottom left (smaller Z, larger phi)
  if (sensN && OKforPhiMax(phiMx,sensN->GetPhiMin()) && sensN->GetZMax()>zMn) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbB));  // neighbor on the bottom (smaller Z)
  if (sensN && sensN->GetZMax()>zMn) sensors[nsens++] = sensN;
  //
  sensN = GetSensor(sens->GetNeighborID(AliITSURecoSens::kNghbBR)); // neighbor on the bottom right (smaller Z, smaller phi)
  if (sensN && OKforPhiMin(phiMn,sensN->GetPhiMax()) && sensN->GetZMax()>zMn) sensors[nsens++] = sensN;
  //
  return nsens;
}

//______________________________________________________
void AliITSURecoLayer::ProcessClusters(Int_t mode)
{
  // register in each sensor of the layer its cluster.
  // the clusters of the layer must be sorted per sensor
  int ncl = fClusters->GetEntriesFast();
  int curSensID = -1;
  for (int i=fNSensors;i--;) GetSensor(i)->SetNClusters(0);
  AliITSURecoSens* curSens = 0;
  for (int icl=0;icl<ncl;icl++) {
    AliITSUClusterPix* cl = (AliITSUClusterPix*) fClusters->UncheckedAt(icl);
    cl->GoToFrameTrk();
    int vID = cl->GetVolumeId();
    if (vID<curSensID) {AliFatal("Clusters are not sorted in increasing sensorID");}
    if (vID>curSensID) {
      if (curSens) curSens->ProcessClusters(mode);    // prepare clusters for reconstruction
      curSens   = GetSensor(vID - fITSGeom->GetFirstModIndex(fActiveID));
      curSensID = vID;
      curSens->SetFirstClusterId(icl);
    }
    curSens->IncNClusters();
  }
  if (curSens) curSens->ProcessClusters(mode); // last sensor was not processed yet
  //
}

//______________________________________________________
Bool_t AliITSURecoLayer::IsEqual(const TObject* obj) const
{
  // check if layers are equal in R
  const AliITSURecoLayer* lr = (const AliITSURecoLayer*)obj;
  return Abs(lr->GetR()-GetR())<1e-6 ? kTRUE : kFALSE;
}

//______________________________________________________
Int_t  AliITSURecoLayer::Compare(const TObject* obj) const
{
  // compare two layers
  const AliITSURecoLayer* lr = (const AliITSURecoLayer*)obj;
  double dr = GetR() - lr->GetR();
  if (Abs(dr)<1e-6) return 0;
  return dr>0 ? 1:-1;
  //      
}

//_________________________________________________________________
AliITSURecoSens* AliITSURecoLayer::GetSensorFromID(Int_t i) const 
{
  // get sensor from its global id
  i -= fITSGeom->GetFirstModIndex(fActiveID);
  if (i<0||i>=fNSensors) AliFatal(Form("Sensor with id=%d is not in layer %d",i+fITSGeom->GetFirstModIndex(fActiveID),fActiveID));
  return (AliITSURecoSens*)fSensors[i];
}
