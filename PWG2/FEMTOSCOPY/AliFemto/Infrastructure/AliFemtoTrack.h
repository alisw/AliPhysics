/***************************************************************************
 *
 * $Id$
 *
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *    Intermediate format for particle.  This is built from the
 *    input particle format (e.g. StTrack of StEvent) and presented to
 *    the Analyses for ParticleCuts.
 *
 ***************************************************************************
 *  ... Lots of stuff deleted here ...
 **************************************************************************/

#ifndef AliFemtoTrack_hh
#define AliFemtoTrack_hh

#include "Infrastructure/AliFemtoTypes.h"
#include "AliFmPhysicalHelixD.h"
#include "TBits.h"
/* Th stuff */
#include "Base/AliFemtoHiddenInfo.h"
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
  AliFemtoTrack(const AliFemtoTrack&);// copy constructor
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
  short TPCsignalN() const;    
  float TPCsignalS() const;   

  const TBits& TPCclusters() const;
  const TBits& TPCsharing()  const;
  
  void SetCharge(const short&);
  void SetPidProbElectron(const float&);
  void SetPidProbPion(const float&);
  void SetPidProbKaon(const float&);
  void SetPidProbProton(const float&);
  void SetPidProbMuon(const float&);
   
  void SetP(const AliFemtoThreeVector&);
  void SetPt(const float&);
  void SetHelix(const AliFmPhysicalHelixD&);
  void SetTrackId(const short&);
  void SetFlags(const long int&);
  void SetLabel(const int&);
  void SetImpactD(const float&);
  void SetImpactZ(const float&);
  void SetCdd(const float&);
  void SetCdz(const float&);
  void SetCzz(const float&);
  
  void SetITSchi2(const float&);    
  void SetITSncls(const int&);     
  void SetTPCchi2(const float&);       
  void SetTPCncls(const int&);       
  void SetTPCnclsF(const short&);      
  void SetTPCsignalN(const short&);    
  void SetTPCsignalS(const float&);   

  void SetTPCcluster(const short& aNBit, const Bool_t& aValue);
  void SetTPCshared(const short& aNBit, const Bool_t& aValue);

  /* Th stuff */
  void SetHiddenInfo(AliFemtoHiddenInfo* aHiddenInfo);
  bool ValidHiddenInfo() const;
  // Fab private : (official : const AliFemtoHiddenInfo* HiddenInfo() const;
  AliFemtoHiddenInfo* getHiddenInfo() const;
  /***/
  
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
  char fCharge;
  float fPidProbElectron; // new
  float fPidProbPion; // new
  float fPidProbKaon; // new
  float fPidProbProton; // new
  float fPidProbMuon; //new 
  unsigned int fTrackId;


  AliFemtoThreeVector fP;
  float fPt;
  AliFmPhysicalHelixD fHelix;

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
  short fTPCnclsF;      // number of findable clusters in the TPC
  short fTPCsignalN;    // number of points used for dEdx
  float fTPCsignalS;    // RMS of dEdx measurement
  TBits fClusters;
  TBits fShared;
  
  /* Th stuff */
  // Fab private : add mutable
  //  mutable 
  AliFemtoHiddenInfo* fHiddenInfo; //!
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
