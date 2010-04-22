// $Id: AliJTrack.h,v 1.3 2008/01/21 11:56:39 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliJTrack.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.1 $
  \date $Date: 2008/05/02 11:56:39 $
*/
////////////////////////////////////////////////////

#ifndef ALIJTRACK_H
#define ALIJTRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


#include "AliPhJBaseTrack.h"
#include "JConst.h"

//class TObject;

class AliJTrack : public AliPhJBaseTrack {

public:

   enum AliJTrkPID {
    kElectronAli     =  0,
    kMuonAli         =  1,
    kPionAli         =  2,
    kKaonAli         =  3,
    kProtonAli       =  4,
    kDeuteronAli     =  5,
    kTritonAli       =  6,
    kHelium3Ali      =  7,
    kAlphaAli        =  8,
    kUnknownAli      =  9,
    kMostProbableAli = -1
   };
  enum {
    kITSin=0x0001,kITSout=0x0002,kITSrefit=0x0004,kITSpid=0x0008,
    kTPCin=0x0010,kTPCout=0x0020,kTPCrefit=0x0040,kTPCpid=0x0080,
    kTRDin=0x0100,kTRDout=0x0200,kTRDrefit=0x0400,kTRDpid=0x0800,
    kTOFin=0x1000,kTOFout=0x2000,kTOFrefit=0x4000,kTOFpid=0x8000,
    kHMPIDout=0x10000,kHMPIDpid=0x20000,
    kEMCALmatch=0x40000,
    kTRDbackup=0x80000,
    kTRDStop=0x20000000,
    kESDpid=0x40000000,
    kTIME=0x80000000,
    kGlobalMerge=0x08000000
  };

  AliJTrack();		         // default constructor
  AliJTrack(const AliJTrack& a); // copy constructor
  ~AliJTrack(){;}		//destructor
	
  void GetPID(Double_t *pid) const {
    for(Int_t i=0; i<10; ++i) pid[i]=ftrkPID[i];
  }
    
  void SetPID(const Double_t *pid); 
  
  void        ConvertAliPID();
  Double_t    GetChi2perNDF() const {return fChi2perNDF;}
  Double_t    GetChi2Trig()   const {return fChi2Trig;}

  void     SetChi2perNDF(Double_t chi2) {fChi2perNDF = chi2;}
  void     SetChi2Trig(Double_t chi2) {fChi2Trig = chi2;}
  
  ULong_t  GetRecFlags() const { return fRecFlags; }
  void     SetRecFlags(ULong_t flags) { fRecFlags = flags; }
  
  Double_t GetTPCdEdx() const {return fTPCdEdx;}
  void     SetTPCdEdx(Double_t dedx) {fTPCdEdx = dedx;}

  int      GetTPCnClust() const {return fTPCnClust;}
  void     SetTPCnClust(int ival) {fTPCnClust = ival;}

  Double_t GetTPCDCAXY() const {return fTPCDCAXY;}
  void     SetTPCDCAXY(Double_t ival) {fTPCDCAXY = ival;}

  Double_t GetTPCDCAZ() const {return fTPCDCAZ;}
  void     SetTPCDCAZ(Double_t ival) {fTPCDCAZ = ival;}

  Double_t GetTPCClustPerFindClust() const {return fTPCClustPerFindClust;}
  void     SetTPCClustPerFindClust(Double_t ival) {fTPCClustPerFindClust = ival;}

  Double_t GetTPCChi2PerClust() const {return fTPCChi2PerClust;}
  void     SetTPCChi2PerClust(Double_t ival) {fTPCChi2PerClust = ival;}

  Double_t GetImapactXY() const {return fImapactXY;}
  void     SetImapactXY(Double_t ival) {fImapactXY = ival;}

  Double_t GetImapactZ() const {return fImapactZ;}
  void     SetImapactZ(Double_t ival) {fImapactZ = ival;}

  Int_t GetKinkIndex() const {return fKinkIndex;}
  void     SetKinkIndex(Int_t ival) {fKinkIndex = ival;}

  UInt_t GetStatus() const {return fstatus;}
  void     SetStatus(UInt_t ival) {fstatus = ival;}

  void GetExternalDiaCovariance(Double_t *ecov) const {
    for(Int_t i=0; i<5; i++) ecov[i]=fextDiaCov[i];
  }
    
  void SetExternalDiaCovariance(const Double_t *ecov);

  AliJTrack& operator=(const AliJTrack& trk);
  
private:

  Double_t   ftrkPID[10];     // [0.,1.,8] pointer to PID object
  Double_t   fChi2perNDF;     // chi2/NDF of mometum fit !! see details 
  Double_t   fChi2Trig;       // chi2 of trigger/track matching   !! see details 
 
  ULong_t    fRecFlags;       // reconstruction status flags 
  
  //TPC 
  Double_t   fTPCdEdx;        // track TPC dE/dx
  int        fTPCnClust;      // track TPC nclusters 
  Double_t   fImapactXY;     // distance of a track to the event vertex in xy plane
  Double_t   fImapactZ;      // distance of a track to the event vertex in z direction
  Double_t   fTPCDCAXY;      // track impact parameter in XY
  Double_t   fTPCDCAZ;        // track impact parameter in Z
  Double_t   fTPCClustPerFindClust; // tpc
  Double_t   fTPCChi2PerClust;  //tpc chi2 per cluster
  //ESD track cuts
  Int_t fKinkIndex;  //kink index ... indication of  kink daughters
  UInt_t fstatus;    // reconstruction flag status
  Double_t fextDiaCov[5];//track parameters covariance matrix
  

  ClassDef(AliJTrack,1)
};

#endif
