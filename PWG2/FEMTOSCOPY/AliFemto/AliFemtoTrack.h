///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoTrack: main class holding all the necessary information       //
// about a track (before the identification) that is required during     //
// femtoscopic analysis. This class is filled with information from the  //
// input stream by the reader. A particle has a link back to the Track   //
// it was created from, so we do not copy the information.               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef ALIFEMTOTRACK_H
#define ALIFEMTOTRACK_H

#include "AliFemtoTypes.h"
#include "AliFmPhysicalHelixD.h"
#include "TBits.h"
/* Th stuff */
#include "AliFemtoHiddenInfo.h"
/***/

#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
class StEvent;
class StTrack;
class StMuDst;
class StMuTrack;
#endif

class AliFemtoTrack{
public:
  AliFemtoTrack();
  AliFemtoTrack(const AliFemtoTrack& aTrack);// copy constructor
  AliFemtoTrack& operator=(const AliFemtoTrack& aTrack);
#ifndef __NO_STAR_DEPENDENCE_ALLOWED__
#ifdef __ROOT__
 //AliFemtoTrack(const StTrack*, AliFemtoThreeVector);   // c-tor from StTrack of STAR DSTs
  //AliFemtoTrack(const StMuDst* dst, const StMuTrack* t);
#endif
  //AliFemtoTrack(const StEvent*, const StTrack*);
#endif

  ~AliFemtoTrack();
//    ~AliFemtoTrack(){/* no-op*/};

  short Charge() const;
  float PidProbElectron() const;
  float PidProbPion() const;
  float PidProbKaon() const;
  float PidProbProton() const;
  float PidProbMuon() const;
  
  AliFemtoThreeVector P() const;
  float Pt() const;
  const AliFmPhysicalHelixD& Helix() const;
  short TrackId() const;
  long int Flags() const;
  int Label()const;
  float ImpactD()const;
  float ImpactZ()const;
  float Cdd() const;
  float Cdz() const;
  float Czz() const; 
  
  float ITSchi2() const;    
  int   ITSncls() const;     
  float TPCchi2() const;       
  int   TPCncls() const;       
  short TPCnclsF() const;      
  float TPCsignal() const;
  short TPCsignalN() const;    
  float TPCsignalS() const;   

  const TBits& TPCclusters() const;
  const TBits& TPCsharing()  const;
  
  void SetCharge(const short& s);
  void SetPidProbElectron(const float& x);
  void SetPidProbPion(const float& x);
  void SetPidProbKaon(const float& x);
  void SetPidProbProton(const float& x);
  void SetPidProbMuon(const float& x);
   
  void SetP(const AliFemtoThreeVector& p);
  void SetPt(const float& x);
  void SetHelix(const AliFmPhysicalHelixD& h);
  void SetTrackId(const short& s);
  void SetFlags(const long int& i);
  void SetLabel(const int& i);
  void SetImpactD(const float& x);
  void SetImpactZ(const float& x);
  void SetCdd(const float& x);
  void SetCdz(const float& x);
  void SetCzz(const float& x);
  
  void SetITSchi2(const float& x);    
  void SetITSncls(const int& i);     
  void SetTPCchi2(const float& x);       
  void SetTPCncls(const int& i);       
  void SetTPCnclsF(const short& s);      
  void SetTPCsignal(const float& s);
  void SetTPCsignalN(const short& s);    
  void SetTPCsignalS(const float& x);   

  void SetTPCcluster(const short& aNBit, const Bool_t& aValue);
  void SetTPCshared(const short& aNBit, const Bool_t& aValue);

  void SetTPCClusterMap(const TBits& aBits);
  void SetTPCSharedMap(const TBits& aBits);

  void SetKinkIndexes(int points[3]);
  int  KinkIndex(int aIndex) const;

  /* Th stuff */
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  bool ValidHiddenInfo() const;
  // Fab private : (official : const AliFemtoHiddenInfo* HiddenInfo() const;
  AliFemtoHiddenInfo* GetHiddenInfo() const;
  /***/
  
  const AliFemtoThreeVector& NominalTpcExitPoint() const;
  const AliFemtoThreeVector& NominalTpcEntrancePoint() const;
    
  void SetNominalTPCEntrancePoint(const AliFemtoThreeVector& aXTPC);
  void SetNominalTPCEntrancePoint(double *aXTPC);

  void SetNominalTPCExitPoint(const AliFemtoThreeVector& aXTPC);
  void SetNominalTPCExitPoint(double *aXTPC);
  void SetSigmaToVertex(const float& Sigma);
  float SigmaToVertex() const;
 
  //Alice stuff
  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kRICHpid=0x20000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000
  };
  
 private:
  char fCharge;           // track charge
  float fPidProbElectron; // electron pid
  float fPidProbPion;     // pion pid 
  float fPidProbKaon;     // kaon pid 
  float fPidProbProton;   // proton pid
  float fPidProbMuon;     // muon pid
  unsigned int fTrackId;  // track unique id


  AliFemtoThreeVector fP; // track momentum
  float fPt;              // transverse momenta
  AliFmPhysicalHelixD fHelix; // track helix
  //alice stuff
  long int fFlags; //Reconsruction status flags
  int fLabel; //Track label  
  float fImpactD; //impact parameter in xy plane
  float fImpactZ;//impacct parameter in z
  float fCdd,fCdz,fCzz;//covariance matrix of the impact parameters
  // ITS related track information
  float fITSchi2;        // chi2 in the ITS
  int   fITSncls;        // number of clusters assigned in the ITS
  // TPC related track information
  float  fTPCchi2;       // chi2 in the TPC
  int    fTPCncls;       // number of clusters assigned in the TPC
  short fTPCnclsF;       // number of findable clusters in the TPC
  float fTPCsignal;      // dEdx TPC value
  short fTPCsignalN;     // number of points used for dEdx
  float fTPCsignalS;     // RMS of dEdx measurement
  float fSigmaToVertex;  // Distance from track to vertex in sigmas
  TBits fClusters;       // Cluster per padrow map
  TBits fShared;         // Sharing per padrow map
  AliFemtoThreeVector fNominalTpcEntrancePoint; // Nominal track entrance point into TPC
  AliFemtoThreeVector fNominalTpcExitPoint;     // Nominal track exit point from TPC

  int   fKinkIndexes[3]; // Kink Index list
  /* Th stuff */
  // Fab private : add mutable
  //  mutable 
  AliFemtoHiddenInfo* fHiddenInfo; //! hidden info
  /***/
};

//inline const float* AliFemtoTrack::NSigma() const 
//{return &mNSigmaElectron;} // Fab private 
inline float AliFemtoTrack::PidProbElectron() const {return fPidProbElectron;}
inline float AliFemtoTrack::PidProbPion() const {return fPidProbPion;}
inline float AliFemtoTrack::PidProbKaon() const {return fPidProbKaon;}
inline float AliFemtoTrack::PidProbProton() const {return fPidProbProton;}
inline float AliFemtoTrack::PidProbMuon() const {return fPidProbMuon;}

#endif
