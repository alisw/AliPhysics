#include <TClonesArray.h>
#include "AliITSURecoLayer.h"
#include "AliITSMFTSegmentationPix.h"
#include "AliITSMFTAux.h"
#include "AliITSUClusterPix.h"
#include "AliITSUGeomTGeo.h"
#include "AliLog.h"

using namespace AliITSMFTAux;
using namespace TMath;

ClassImp(AliITSURecoLayer)


//______________________________________________________
AliITSURecoLayer::AliITSURecoLayer(const char* name)
 :fActiveID(-1)
  ,fNSensors(0)
  ,fNSensorRows(0)
  ,fNSensorsPerRow(0)
  ,fSensVIDtoMatrixID(0)
  ,fR(0)
  ,fRMax(0)
  ,fRMin(0)
  ,fZMax(0)
  ,fZMin(0)
  ,fPhiOffs(0)
  ,fSensDZInv(0)
  ,fSensDPhiInv(0)
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
  ,fNSensorRows(0)
  ,fNSensorsPerRow(0)
  ,fSensVIDtoMatrixID(0)
  ,fR(0)
  ,fRMax(0)
  ,fRMin(0)
  ,fZMax(0)
  ,fZMin(0)
  ,fPhiOffs(0)
  ,fSensDZInv(0)
  ,fSensDPhiInv(0)
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
  delete fSensors;
  delete[] fSensVIDtoMatrixID;
  if (GetOwnsClusterArray()) delete fClusters;
}

//______________________________________________________
void AliITSURecoLayer::Print(Option_t* opt) const			      
{
  //print 
  printf("Lr %-15s %d (act:%+d), NSens: %4d in %3d rows| MaxStep:%.2f ",
	 GetName(),GetID(),GetActiveID(),GetNSensors(),GetNSensorRows(),fMaxStep);
  printf("%6.3f<R<%6.3f | %+8.3f<Z<%+8.3f dZ:%6.3f dPhi:%6.3f\n",fRMin,fRMax,fZMin,fZMax,
	 fSensDZInv>0 ? 1/fSensDZInv : 0, fSensDPhiInv>0 ? 1/fSensDPhiInv : 0);
  TString opts = opt; opts.ToLower();
  if (opts.Contains("sn")) for (int i=0;i<GetNSensors();i++) GetSensor(i)->Print(opt);
}

//______________________________________________________
void AliITSURecoLayer::Build()
{
  // build internal structures
  const double kSafeR = 0.05; // safety margin for Rmin,Rmax of the layer
  if (fActiveID<0) return;
  //
  int nStaves=fITSGeom->GetNStaves(fActiveID);
  // determine number of sensor rows (sensors aligned at same phi and spanning the Z range of the layer)
  fNSensorRows = nStaves;
  //
  // if the stave has susbtaves, each substave can have multiple rows of sensors (but just 1 row of modules)
  if (fITSGeom->GetNHalfStaves(fActiveID)>0) fNSensorRows *= fITSGeom->GetNHalfStaves(fActiveID);
  //
  // if there are modules defined, the module may have multiple rows of sensors (though not spanning full Z)
  if (fITSGeom->GetNModules(fActiveID)>0) fNSensorRows *= fITSGeom->GetNChipRowsPerModule(fActiveID);
  //
  fNSensors = fITSGeom->GetNChipsPerLayer(fActiveID);
  fNSensorsPerRow = fNSensors/fNSensorRows;
  //
  fSensors = new TObjArray(fNSensors);
  fSensVIDtoMatrixID = new Int_t[fNSensors];
  const AliITSMFTSegmentationPix* kSegm = fITSGeom->GetSegmentation(fActiveID);
  //
  TGeoHMatrix mmod;
  const TGeoHMatrix* mt2l;
  double phiTF,rTF, loc[3]={0,0,0},glo[3];
  //
  int nSensPerStave = fITSGeom->GetNChipsPerStave(fActiveID);
  for (int staveI=0;staveI<nStaves;staveI++) {
    for (int sensI=0;sensI<nSensPerStave;sensI++) {
      int sID = fITSGeom->GetChipIndex(fActiveID,staveI,sensI);
      AliITSURecoSens* sens = new AliITSURecoSens( sID );
      fSensors->AddLast(sens);
      double phiMin=1e9,phiMax=-1e9,zMin=1e9,zMax=-1e9;
      // this is NOT the sensor matrix, just the ideal chip matrix to get neighbors correct
      fITSGeom->GetOrigMatrix(sID,mmod); 
      //
      for (int ix=0;ix<2;ix++) {       // determine sensor boundaries (ideal)
	loc[0] = (ix-0.5)*kSegm->Dx(); // +-DX/2
	for (int iy=0;iy<2;iy++) {
	  loc[1] = (iy-0.5)*kSegm->Dy(); // +-DY/2
	  for (int iz=0;iz<2;iz++) {
	    loc[2] = (iz-0.5)*kSegm->Dz(); // +-DZ/2
	    mmod.LocalToMaster(loc,glo);
	    double phi = ATan2(glo[1],glo[0]);
	    BringTo02Pi(phi);
	    if      (phiMin>1e8)  phiMin=phi;
	    else if (!OKforPhiMin(phiMin,phi)) phiMin=phi;
	    if      (phiMax<-1e8) phiMax=phi; 
	    else if (!OKforPhiMax(phiMax,phi)) phiMax=phi;	      
	    if (glo[2]>zMax) zMax=glo[2];
	    if (glo[2]<zMin) zMin=glo[2];
	  }
	}
      }
      sens->SetBoundaries(phiMin,phiMax,zMin,zMax);
    }
  }
  fSensors->Sort(); // sort sensors to get the neighborhood correct
  //
  // now fill real sensor angles, Z's, accounting for misalignment
  fRMin=fZMin=1e9;
  fRMax=fZMax=-1e9;
  //
  fPhiOffs = 0;
  int firstSensID = fITSGeom->GetFirstChipIndex(fActiveID);
  for (int sensI=0;sensI<fNSensors;sensI++) {
    AliITSURecoSens* sens = GetSensor(sensI);
    mmod = *fITSGeom->GetMatrixSens(sens->GetID());
    fSensVIDtoMatrixID[sens->GetID() - firstSensID] = sensI;
    double phiMin=1e9,phiMax=-1e9,zMin=1e9,zMax=-1e9;
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
    mt2l = fITSGeom->GetMatrixT2L( sens->GetID() );
    mmod.Multiply(mt2l);	
    loc[0]=loc[1]=loc[2]=0;
    mmod.LocalToMaster(loc,glo);
    rTF   = Sqrt(glo[0]*glo[0] + glo[1]*glo[1]);  //  tracking params (misaligned)
    phiTF = ATan2(glo[1],glo[0]);
    BringTo02Pi(phiTF);
    //
    sens->SetXTF(rTF);
    sens->SetPhiTF(phiTF);
    sens->SetBoundaries(phiMin,phiMax,zMin,zMax);
    if (fZMin>zMin) fZMin = zMin;
    if (fZMax<zMax) fZMax = zMax;
    //
    if (sensI<fNSensorsPerRow) fPhiOffs += MeanPhiSmall(phiMax,phiMin);
  }
  //
  fPhiOffs /= fNSensorsPerRow; // average phi of the 1st row
  fSensDZInv = fNSensorsPerRow/(fZMax-fZMin);
  fSensDPhiInv = fNSensorRows/(2*Pi());
  //
  fRMin = Sqrt(fRMin);
  fRMax = Sqrt(fRMax);
  fR = 0.5*(fRMin+fRMax);
  fRMin -= kSafeR;
  fRMax += kSafeR;
  //
}

//______________________________________________________
Int_t AliITSURecoLayer::FindSensors(const double* impPar, AliITSURecoSens *sensors[kMaxSensMatching])
{
  // find sensors having intersection with track
  // impPar contains: lab phi of track, dphi, labZ, dz
  //
  double zMn=impPar[2]-impPar[3], zMx=impPar[2]+impPar[3]; 
  if (zMn>fZMax) return 0;
  if (zMx<fZMin) return 0;
  //
  int zCenID = int((impPar[2]-fZMin)*fSensDZInv);
  if      (zCenID<0) zCenID = 0;
  else if (zCenID>=fNSensorsPerRow) zCenID = fNSensorsPerRow-1;
  double phiCn = impPar[0] - fPhiOffs;
  //  BringTo02Pi(phiCn); 
  int rowCenID = int(phiCn*fSensDPhiInv);
  //
  // due to the misalignments the actual sensorID's might be shifted
  int res = 0;
  AliITSURecoSens* sensPrev=0, *sens = GetSensor(rowCenID,zCenID);
  //
  //  printf("Guess: Primary Sensor: phiID: %d zID: %d ->",rowCenID,zCenID); sens->Print();
  //
  while ( (res=sens->CheckCoverage(impPar[0], impPar[2])) ) {
    if      (res&AliITSURecoSens::kRight) {if (++rowCenID==fNSensorRows) rowCenID=0;} // neighbor on the right (larger phi)
    else if (res&AliITSURecoSens::kLeft)  {if (--rowCenID<0) rowCenID = fNSensorRows-1;}      // neighbor on the left (smaller phi)
    if      (res&AliITSURecoSens::kUp)    {if (++zCenID==fNSensorsPerRow) zCenID = fNSensorsPerRow-1;}    // neighbor at larger Z (if any)
    else if (res&AliITSURecoSens::kDown)  {if (--zCenID<0) zCenID = 0;}               // neighbor at smaller Z (if any)
    //
    AliITSURecoSens* sensAlt = GetSensor(rowCenID,zCenID);
    if (sensAlt==sens || sensAlt==sensPrev) break;  // there is no better neighbor (z edge) or the point falls in dead area
    //
    sensPrev = sens;
    sens = sensAlt;
  }
  //  printf("Found: Primary Sensor: phiID: %d zID: %d ->",rowCenID,zCenID); sens->Print();
  //
  int nFnd = 0;
  sensors[nFnd++] = sens;
  //
  double phiMn = impPar[0]-impPar[1], phiMx = impPar[0]+impPar[1];
  BringTo02Pi(phiMn);
  BringTo02Pi(phiMx);
  //
  const int kNNeighb = 8;
  const int kCheckNeighb[2][kNNeighb] = { // phi and Z neighbours to check
    { 1, 1, 0,-1,-1,-1, 0, 1},
    { 0, 1, 1, 1, 0,-1,-1,-1}
  };
  //  
  //  printf("Search: %+.4f %+.4f | %+.4f %+.4f\n",phiMn,phiMx, zMn,zMx);
  for (int inb=kNNeighb;inb--;) {
    int idz = kCheckNeighb[1][inb];
    int iz   = zCenID   + idz;
    //    printf("#%d  dp:%+d dz:%+d IZ: %d\n",inb, kCheckNeighb[0][inb], kCheckNeighb[1][inb], iz);
    //
    if (iz<0 || iz>=fNSensorsPerRow) continue;
    int idphi = kCheckNeighb[0][inb];
    int iphi = rowCenID + idphi;
    if      (iphi<0) iphi += fNSensorRows;
    else if (iphi>=fNSensorRows) iphi -= fNSensorRows;
    sens = GetSensor(iphi,iz);
    //
    if      (idz>0) {if (zMx<sens->GetZMin()) continue;}
    else if (idz<0) {if (zMn>sens->GetZMax()) continue;}
    //
    // Z range matches
    if      (idphi>0) {if (!OKforPhiMin(sens->GetPhiMin(),phiMx)) continue;}
    else if (idphi<0) {if (!OKforPhiMax(sens->GetPhiMax(),phiMn)) continue;}
    //
    //    printf("Add %d\n",nFnd);
    sensors[nFnd++] = sens;
    if (nFnd==kMaxSensMatching) break;
  }
  return nFnd;
}

//*/
/*
Int_t AliITSURecoLayer::FindSensors(const double* impPar, AliITSURecoSens *sensors[AliITSURecoSens::kNNeighbors],int mcLab)
{
  // find sensors having intersection with track
  // impPar contains: lab phi of track, dphi, labZ, dz

  //tmp>>>
  int nFnd = 0;
  int fndSens[50];
  if (mcLab>=0) { // find correct sensors from MC info
    int ncl = GetNClusters();
    for (int icl=ncl;icl--;) {
      AliCluster* cl = GetCluster(icl);
      for (int ilb=0;ilb<3;ilb++) {
	if (cl->GetLabel(ilb)<0) break;
	if (cl->GetLabel(ilb)==mcLab) {fndSens[nFnd++] = cl->GetVolumeId(); break;}
      }
    }
    if (nFnd>0) {
      Int_t layS,staS,sensS;
      for (int is=0;is<nFnd;is++) {
	fITSGeom->GetChipId(fndSens[is],layS,staS,sensS);
	printf("SNMC#%d(%d): %d %d %d | ",is,mcLab,layS,staS,sensS); GetSensorFromID(fndSens[is])->Print();
      }
    }
  }

  //tmp<<<

  double z = impPar[2];
  if (z>fZMax+impPar[3]) {
    if (nFnd>0) printf("MissedSens!!!\n");
    return 0; // outside of Z coverage
  }
  z -= fZMin;
  if (z<-impPar[3]) {
    if (nFnd>0) printf("MissedSens!!!\n");
    return 0; // outside of Z coverage
  }
  int sensInSta = int(z*fSensDZInv);
  if      (sensInSta<0) sensInSta = 0;
  else if (sensInSta>=fNSensInStave) sensInSta = fNSensInStave-1;
  //
  double phi = impPar[0] - fPhiOffs;
  BringTo02Pi(phi);
  int staID = int(phi*fSensDPhiInv);  // stave id
  int nsens = 0;
  //
  AliITSURecoSens* sensN,*sens = GetSensor(staID*fNSensInStave+sensInSta);
  sensors[nsens++] = sens;
  //
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
  if (mcLab>=0) {
    Int_t layS,staS,sensS;
    printf("Found %d sensors for phi %.3f : %.3f | Z %.4f %.4f\n", nsens,phiMn,phiMx,zMn,zMx); 
    for (int is=0;is<nsens;is++) {
      fITSGeom->GetChipId(sensors[is]->GetID()+fITSGeom->GetFirstModIndex(fActiveID),layS,staS,sensS);
      printf("*SNF#%d: %d %d %d | ",is,layS,staS,sensS); sensors[is]->Print();
    }
    for (int ism=0;ism<nFnd;ism++) {
      AliITSURecoSens* snMC = GetSensorFromID(fndSens[ism]);
      Bool_t ok=kFALSE;
      for (int isf=0;isf<nsens;isf++) {
	if (snMC==sensors[isf]) {ok=kTRUE;break;}
      }
      if (!ok) printf("MissedSens %d!!!\n",ism);
    }
  }
  return nsens;
}
*/
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
    cl->SetRecoInfo(0);
    cl->GoToFrameTrk();
    int vID = cl->GetVolumeId();
    if (vID<curSensID) {AliFatal("Clusters are not sorted in increasing sensorID");}
    if (vID>curSensID) {
      if (curSens) curSens->ProcessClusters(mode);    // prepare clusters for reconstruction
      curSens   = GetSensorFromID(vID); //GetSensor(vID - fITSGeom->GetFirstChipIndex(fActiveID));
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
  i -= fITSGeom->GetFirstChipIndex(fActiveID);
  if (i<0||i>=fNSensors) AliFatal(Form("Sensor with id=%d is not in layer %d",i+fITSGeom->GetFirstChipIndex(fActiveID),fActiveID));
  return GetSensor(SensVIDtoMatrixID(i));
}
