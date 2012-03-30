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

#include "AliJBaseTrack.h"

class AliJTrack : public AliJBaseTrack {

public:
   enum AliJTrkPID {
    kMostProbableAliJ = -1, 
    kElectronAliJ, 
    kMuonAliJ,
    kPionAliJ, 
    kKaonAliJ, 
    kProtonAliJ, 
    kNAliJTrkPID
   };

   enum AliJTrkPIDmethod {
    kTOF=0, kTPC=1, kTPCTOF=2, kNAliJTrkPIDmethod
   };

  AliJTrack();		         // default constructor
  AliJTrack(const AliJTrack& a); // copy constructor
  ~AliJTrack(){;}		//destructor

  Double32_t  GetPID(AliJTrkPID p, AliJTrkPIDmethod m) const { return fTrkPID[p][m]; }
  UInt_t      GetFilterMap() const { return fFilterMap; }
  Bool_t      IsFiltered( int i ) const { return TESTBIT( fFilterMap, i); }
  Bool_t      IsFilteredMask( UInt_t mask ) const { return ((Bool_t)(((fFilterMap) & mask) != 0)); };
  int         GetTPCnClust() const {return fTPCnClust;}
  Float_t     GetTPCdEdx()  const {return fTPCdEdx; }
  Float_t  GetTOFbeta() const {return fTOFbeta;}
  Float_t  GetExpectedTOFbeta( AliJTrkPID p ) const {return fExpTOFbeta[p];}
  Double_t  GetExpectedTPCdEdx( AliJTrkPID p ) const {return fExpTPCdEdx[p];}
  Double_t  GetTPCsigma( AliJTrkPID p) const {return fTPCsigma[p];}
  Double_t  GetTOFsigma( AliJTrkPID p) const {return fTOFsigma[p];}

  void SetPID(AliJTrkPID p, double pro, AliJTrkPIDmethod m){ fTrkPID[p][m] = pro; }
  void SetFilterMap( UInt_t map ){ fFilterMap = map; }
  void SetFilterMap( int i, bool t ){ t?SETBIT(fFilterMap,i):CLRBIT(fFilterMap,i); } 
  void SetTPCnClust(int ival) {fTPCnClust = ival;}
  void SetTPCdEdx( Float_t dedx ){ fTPCdEdx = dedx; }
  void SetTOFbeta(Float_t beta) {fTOFbeta = beta;}
  void SetExpectedTOFbeta(AliJTrkPID p, float exbeta) {fExpTOFbeta[p] = exbeta;}
  void SetExpectedTPCdEdx(AliJTrkPID p, Double_t exdEdx) {fExpTPCdEdx[p] = exdEdx;}
  void SetTPCsigma(AliJTrkPID p, Double_t TPCsigma) {fTPCsigma[p] = TPCsigma;}
  void SetTOFsigma(AliJTrkPID p, Double_t TOFsigma) {fTOFsigma[p] = TOFsigma;} // Return the expected sigma of the PID signal for the specified  particle type.

  void SetTPCTrack( Float_t px,  Float_t py,  Float_t pz ){ fTPCTrack[0]=px;fTPCTrack[1]=py;fTPCTrack[2]=pz;}
  void SetUseTPCTrack(){ SetPxPyPzE( fTPCTrack[0], fTPCTrack[1], fTPCTrack[2], 0 ); } 
  Float_t* GetTPCTrack(){ return fTPCTrack; }

  AliJTrack& operator=(const AliJTrack& trk);

private:
  Double32_t  fTrkPID[kNAliJTrkPID][kNAliJTrkPIDmethod];    //[0.,1.,8] Array for PID. 
  UInt_t      fFilterMap;               // bit serious of cuts
  Short_t     fTPCnClust;               // track TPC nclusters 
  Float_t     fTPCdEdx;
  Float_t    fTOFbeta;                    //!
  Float_t    fExpTOFbeta[kNAliJTrkPID];   //!
  Double_t   fExpTPCdEdx[kNAliJTrkPID];   //!
  Double_t   fTPCsigma[kNAliJTrkPID];     //!
  Double_t   fTOFsigma[kNAliJTrkPID];     //!
  Float_t    fTPCTrack[3];              // px, py, pz for TPCTrack;

  ClassDef(AliJTrack,1)
};

#endif
