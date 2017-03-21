// $Id$
#ifndef ALIHLTTOFRAWCLUSTER_H
#define ALIHLTTOFRAWCLUSTER_H

#include "Rtypes.h"
#include "AliTOFcluster.h"
#include "AliTOFClusterFinder.h"

/**
 * @struct AliHLTTOFRawCluster
 * Primitive data of a TOF cluster in raw coordinates. The plan is to store the
 * data in a compressed format by limiting the resolution of the float values.
 * @ingroup alihlt_tof_datastructs
 */
struct AliHLTTOFRawCluster {

AliHLTTOFRawCluster()
  :   fX(0)
    , fY(0)
    , fZ(0)
    , fSigmaY2(0)
    , fSigmaZ2(0) 
    , fSigmaYZ(0)
    , fVolumeId(0)
    , fIsMisaligned(kFALSE)
    , fSigmaX2(0)
    , fSigmaXY(0)
    , fSigmaXZ(0)
    , fIdx(-1)
    , fQuality(-100)
    , fR(0)
    , fPhi(0)
    , fTDC(0.)
    , fToT(0)
    , fADC(0)
    , fTdcND(0)
    , fTdcRAW(0)
    , fStatus(kTRUE)
    , fDeltaBC(0)
    , fL0L1Latency(0)
    , fESDID(-1)
  {  
  for (int ii=0; ii<3; ii++) fTracks[ii] = -3141593;
  for (int ii=0; ii<5; ii++) fdetIndex[ii] = -1;
  }

AliHLTTOFRawCluster(short volId, 
		    float x,   float y,   float z,
		    float sx2, float sxy, float sxz,
		    float sy2, float syz, 
		    float sz2, int *lab,
		    int * const ind,
		    int *par, bool status, int idx) {

  // Constructor
  
  AliTOFcluster* temp = new AliTOFcluster(volId, x, y, z, sx2, sxy, sxz, sy2, syz, sz2, lab, ind, par, status, idx);
  
  fX = temp->GetX();
  fY = temp->GetY();
  fZ = temp->GetZ();
  fSigmaY2 = temp->GetSigmaY2();
  fSigmaZ2 = temp->GetSigmaZ2();
  fSigmaYZ = temp->GetSigmaYZ();
  fVolumeId = temp->GetVolumeId();
  fIsMisaligned = temp->Misalign();
  // AliCluster3D
  fSigmaX2 = temp->GetSigmaX2();
  fSigmaXY = temp->GetSigmaXY();
  fSigmaXZ = temp->GetSigmaXZ();
  // AliTOFcluster   
  fIdx = temp->GetIndex();
  fQuality = temp->GetQuality();
  fR = temp->GetR();
  fPhi = GetPhi();
  fTDC = GetTDC();
  fToT = GetToT();
  fADC = GetADC();
  fTdcND = GetTDCND(); 
  fTdcRAW = GetTDCRAW();
  fStatus = GetStatus(); 
  fDeltaBC = GetDeltaBC();
  fL0L1Latency = GetL0L1Latency(); 
  fESDID = GetESDID(); 
  
  if (lab) {
    fTracks[0] = lab[0];
    fTracks[1] = lab[1];
    fTracks[2] = lab[2];
  }
  else
    fTracks[0]=fTracks[1]=fTracks[2]=-3141593;
  
}

  
  AliHLTTOFRawCluster(const AliHLTTOFRawCluster& other)
      : fX(other.fX)
      , fY(other.fY)
      , fZ(other.fZ)
      , fSigmaY2(other.fSigmaY2)
      , fSigmaZ2(other.fSigmaZ2)
      , fSigmaYZ(other.fSigmaXZ)
      , fVolumeId(other.fVolumeId)
      , fIsMisaligned(other.fIsMisaligned)
      // AliCluster3D
      , fSigmaX2(other.fSigmaX2)
      , fSigmaXY(other.fSigmaXY)
      , fSigmaXZ(other.fSigmaXZ)
      // AliTOFcluster   
      , fIdx(other.fIdx)
      , fQuality(other.fQuality) 
      , fR(other.fR)
      , fPhi(other.fPhi)
      , fTDC(other.fTDC)
      , fToT(other.fToT)
      , fADC(other.fADC)
      , fTdcND(other.fTdcND)
      , fTdcRAW(other.fTdcRAW)
      , fStatus(other.fStatus)
      , fDeltaBC(other.fDeltaBC)
      , fL0L1Latency(other.fL0L1Latency)
      , fESDID(other.fESDID)  {

	for (int ii=0; ii<5; ii++) fdetIndex[ii] = other.fdetIndex[ii];
      }

  AliHLTTOFRawCluster& operator=(const AliHLTTOFRawCluster& other) {
    if (this==&other) return *this;
    this->~AliHLTTOFRawCluster();
    new (this) AliHLTTOFRawCluster(other);
    return *this;
  }

  void Clear() {
    this->~AliHLTTOFRawCluster();
    new (this) AliHLTTOFRawCluster;
  }

  // AliCluster
  int    fTracks[3];     // MC labels
  float  fX;           // X of the cluster in the tracking c.s.
  float  fY;        // Y of the cluster in the tracking c.s.
  float  fZ;        // Z of the cluster in the tracking c.s.
  float  fSigmaY2;  // Sigma Y square of cluster
  float  fSigmaZ2;  // Sigma Z square of cluster
  float  fSigmaYZ;  // Non-diagonal element of cov.matrix [y, z]
  short  fVolumeId; // Volume ID of the detector element
  bool   fIsMisaligned; // Cluster was misagned or not?

  // AliCluster3D
  float  fSigmaX2;  // Sigma X square of cluster
  float  fSigmaXY;  // Non-diagonal element of cov.matrix [x, y]
  float  fSigmaXZ;  // Non-diagonal element of cov.matrix

  // AliTOFcluster
  int fIdx;         // index of the digit related to this cluster
  int fdetIndex[5]; // Cluster detector indices
                      // (sector, plate, strip, padz, padx)
  // Cluster Quality
  double fQuality;  // quality of the best track 
  // Cluster Global Position
  double fR;        // r-coordinate
  double fPhi;      // phi-coordinate
  // TOF Signal parameters
  int  fTDC;        // TDC count
  int  fToT;        // ToT
  int  fADC;        // ADC count
  int  fTdcND;      // TDC count
  int  fTdcRAW;     // RAW TDC count
  bool fStatus;     // cluster online status 
  int  fDeltaBC; // deltaBC
  int  fL0L1Latency; // L0L1 latency
  //
  int fESDID;      //! id in ESD clusters list (if stored)

  // Getters and Setters
  // AliCluster
  int    GetLabel(int i) const {return fTracks[i];}
  float  GetX()            const {return fX;}
  float  GetY()            const {return fY;}
  float  GetZ()            const {return fZ;}
  float  GetSigmaY2()      const {return fSigmaY2;}
  float  GetSigmaZ2()      const {return fSigmaZ2;}
  float  GetSigmaYZ()      const {return fSigmaYZ;}
  short  GetVolumeId()     const {return fVolumeId;}

  // AliCluster3D
  float GetSigmaX2() const {return fSigmaX2;}
  float GetSigmaXY() const {return fSigmaXY;}
  float GetSigmaXZ() const {return fSigmaXZ;}

  // AliTOFcluster 
  double GetR() const   {return fR;}   // Cluster Radius
  double GetPhi() const {return fPhi;} // Cluster Phi

  double GetQuality() const {return fQuality;} // Cluster quality getter
  bool   GetStatus()  const {return fStatus;}  // Cluster status getter
  int    GetToT()    const {return fToT;}    // Cluster Charge getter
  int    GetTDC()    const {return fTDC;}    // Cluster ToF getter
  int    GetTDCND()  const {return fTdcND;}  // Cluster ToF getter
  int    GetTDCRAW() const {return fTdcRAW;} // Cluster Raw time getter
  int    GetADC()    const {return TMath::Abs(fADC);} // Cluster Charge getter
  int    IsUsed()    const {return (fADC<0) ? 1 : 0;} // Flagging
  int    GetDetInd(int n) const {return fdetIndex[n];} // Cluster Detector Indices getter
  int    GetIndex()  const {return fIdx;} // Digit Index getter
  int    GetDeltaBC() const {return fDeltaBC;}; // deltaBC
  int    GetL0L1Latency() const {return fL0L1Latency;}; // L0L1 latency
  int GetESDID()  const  {return fESDID;}  
  
  // Setters
  // AliCluster
  void     SetLabel(int lab, int i)
  { if (i>=0 && i<3) fTracks[i] = lab;}
  void     SetX(float x) {fX = x;}
  void     SetY(float y) {fY = y;}
  void     SetZ(float z) {fZ = z;}
  void     SetSigmaY2(float sigy2) {fSigmaY2 = sigy2;}
  void     SetSigmaZ2(float sigz2) {fSigmaZ2 = sigz2;}
  void     SetSigmaYZ(float sigyz) {fSigmaYZ = sigyz;};
  void     SetVolumeId(short id)  {fVolumeId = id;}

  // AliCluster3D
  void     SetSigmaX2(float sigx2) {fSigmaX2 = sigx2;}
  void     SetSigmaXY(float sigxy) {fSigmaXY = sigxy;}
  void     SetSigmaXZ(float sigxz) {fSigmaXZ = sigxz;};

  // AliTOFcluster
  void  Use(int = 0) {fADC=-fADC;}        //  setter
  void  SetQuality(double quality) {fQuality = quality;} // Cluster quality setter
  void  SetStatus(bool status) {fStatus = status;}       // Cluster status setter
  void  SetToT(int ToT) {fToT = ToT;}       // Cluster ToT setter
  void  SetTDC(int Tdc) {fTDC = Tdc;}       // Cluster ToF setter
  void  SetTDCND(int Tdc) {fTdcND = Tdc;}   // Cluster ToFnd setter
  void  SetTDCRAW(int Tdc) {fTdcRAW = Tdc;} // Cluster ToF-raw setter
  void  SetDeltaBC(int value) {fDeltaBC = value;}; // deltaBC
  void  SetL0L1Latency(int value) {fL0L1Latency = value;}; // L0-L1 latency
  void  SetDetInd(int* ind){for (int ii = 0; ii<5; ii++) fdetIndex[ii] = ind[ii];}; // detector indices
  void  SetDetInd(int ii, int index){fdetIndex[ii] = index;}; // detector index
  
  //
  void  SetESDID(int id) {fESDID = id;}


};

typedef struct AliHLTTOFRawCluster AliHLTTOFRawCluster;

struct AliHLTTOFRawClusterData
{
  UInt_t fVersion; // version number
  UInt_t fCount;   // number of clusters
#if defined(__HP_aCC) || defined(__DECCXX) || defined(__SUNPRO_CC)
  AliHLTTOFRawCluster  fClusters[1]; // array of clusters  
#else
  AliHLTTOFRawCluster  fClusters[0]; // array of clusters 
#endif
};
typedef struct AliHLTTOFRawClusterData AliHLTTOFRawClusterData;

#endif
