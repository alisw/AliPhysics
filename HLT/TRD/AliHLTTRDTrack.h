#ifndef ALIHLTTRDTRACK_H
#define ALIHLTTRDTRACK_H

#include "AliESDtrack.h"
#include "AliTRDtrackV1.h"
#include "AliHLTLogging.h"
class AliHLTTRDTracklet;


class AliHLTTRDTrack
{
 public:
  AliHLTTRDTrack();
  ~AliHLTTRDTrack();
  AliHLTTRDTrack( AliTRDtrackV1* inTrack);

  void AddTracklets();
  void CopyDataMembers();
  void ExportTRDTrack(AliTRDtrackV1* outTrack);
  AliHLTUInt8_t *GetEndPointer() // Returns pointer to the end of the track
    { return ((AliHLTUInt8_t *) this + fSize); };
  AliHLTUInt32_t GetSize(){return fSize;};
  void Print(Bool_t printTracklets = kTRUE);
  void ReadTrackletsFromMemory(void* );
  
 private:
  AliHLTTRDTrack(const AliHLTTRDTrack& inTrack);
  AliHLTTRDTrack& operator=(const AliHLTTRDTrack& inTrack);
  void InitArrays();
  
  AliHLTUInt32_t fSize; // Size of the track with tracklets and clusters in the memory
  AliTRDtrackV1* fTRDtrack;
  
  /* ======== From AliTRDtrackV1 ======== */
    enum { kNdet      = 540
        , kNstacks   =  90
        , kNplane    =   AliESDtrack::kTRDnPlanes
        , kNcham     =   5
        , kNsect     =  18
        , kNslice    =   3
        , kNMLPslice =   8 };

  /* Defenitely need */
  UChar_t      fPIDquality;           //  No of planes used for PID calculation	
  AliHLTTRDTracklet *fTracklet[kNplane];   //  Tracklets array defining the track
  
  /* Probably need */
  Double32_t   fPID[AliPID::kSPECIES];//  PID probabilities
  Double32_t   fBudget[3];            //  Integrated material budget
  Double32_t   fDE;                   //  Integrated delta energy

  /* ======== From AliKalmanTrack ======== */
  
  /* Defenitely need */
  Double32_t fFakeRatio;  // fake ratio
  Double32_t fChi2;       // total chi2 value for this track
  Double32_t fMass;       // mass hypothesis
  Int_t fLab;             // track label

  /* Probably need */
  Int_t fN;               // number of associated clusters
  Double32_t fIntegratedLength;        // integrated length  // variables for time integration (S.Radomski@gsi.de)

  /* ======= From AliExternalTrackParam ======== */

  /* Defenitely need */
  Double32_t           fX;     // X coordinate for the point of parametrisation
  Double32_t           fAlpha; // Local <-->global coor.system rotation angle
  Double32_t           fP[5];  // The track parameters
  Double32_t           fC[15]; // The track parameter covariance matrix

  /* Probably need */
  //  static Double32_t    fgMostProbablePt; // "Most probable" pt (to be used if Bz=0)

};

#endif
