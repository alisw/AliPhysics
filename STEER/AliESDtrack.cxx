/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-----------------------------------------------------------------
//           Implementation of the ESD track class
//   ESD = Event Summary Data
//   This is the class to deal with during the phisics analysis of data
//      Origin: Iouri Belikov, CERN
//      e-mail: Jouri.Belikov@cern.ch
//-----------------------------------------------------------------

#include "TMath.h"

#include "AliESDtrack.h"
#include "AliKalmanTrack.h"
#include "AliLog.h"

ClassImp(AliESDtrack)

//_______________________________________________________________________
AliESDtrack::AliESDtrack() : 
  TObject(),
  fFlags(0),
  fLabel(0),
  fID(0),
  fTrackLength(0),
  fD(0),
  fZ(0),
  fStopVertex(0),
  fRalpha(0),
  fRx(0),
  fCalpha(0),
  fCx(0),
  fCchi2(1e10),
  fIalpha(0),
  fIx(0),
  fTalpha(0),
  fTx(0),
  fITSchi2(0),
  fITSncls(0),
  fITSsignal(0),
  fITSLabel(0),
  fITSFakeRatio(0),
  fITStrack(0),
  fTPCchi2(0),
  fTPCncls(0),
  fTPCClusterMap(159),//number of padrows
  fTPCsignal(0),
  fTPCLabel(0),
  fTRDchi2(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDsignal(0),
  fTRDLabel(0),
  fTRDQuality(0),
  fTRDtrack(0),
  fTOFchi2(0),
  fTOFindex(0),
  fTOFsignal(-1),
  fPHOSsignal(-1),
  fEMCALsignal(-1),
  fRICHchi2(1e10),
  fRICHncls(0),
  fRICHindex(0),
  fRICHsignal(-1),
  fRICHtheta(0),
  fRICHphi(0),
  fRICHdx(0),
  fRICHdy(0)
{
  //
  // The default ESD constructor 
  //
  for (Int_t i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=1.;
    fITSr[i]=1.;
    fTPCr[i]=1.;
    fTRDr[i]=1.;
    fTOFr[i]=1.;
    fRICHr[i]=1.;
  }
  
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) {
    fPHOSr[i]  = 1.;
    fEMCALr[i] = 1.;
  }

 
  fPHOSpos[0]=fPHOSpos[1]=fPHOSpos[2]=0.;
  fEMCALpos[0]=fEMCALpos[1]=fEMCALpos[2]=0.;
  Int_t i;
  for (i=0; i<5; i++)  { 
    fRp[i]=fCp[i]=fIp[i]=fTp[i]=0.;
  }
  for (i=0; i<15; i++) { 
    fRc[i]=fCc[i]=fIc[i]=fTc[i]=0.;  
  }
  for (i=0; i<6; i++)  { fITSindex[i]=0; }
  for (i=0; i<180; i++){ fTPCindex[i]=0; }
  for (i=0; i<3;i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3;i++)   { fV0Indexes[i]=-1;}
  for (i=0; i<130; i++) { fTRDindex[i]=0; }
  for (i=0;i<kNPlane;i++) {fTRDsignals[i]=0.; fTRDTimBin[i]=-1;}
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=-1;}
  for (Int_t i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (Int_t i=0;i<10;i++) {fTOFInfo[i]=-1;}
  fTPCLabel = 0;
  fTRDLabel = 0;
  fTRDQuality =0;
  fITSLabel = 0;
  fITStrack = 0;
  fTRDtrack = 0;  
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack(const AliESDtrack& track):
  TObject(track),
  fFlags(track.fFlags),
  fLabel(track.fLabel),
  fID(track.fID),
  fTrackLength(track.fTrackLength),
  fD(track.fD),
  fZ(track.fZ),
  fStopVertex(track.fStopVertex),
  fRalpha(track.fRalpha),
  fRx(track.fRx),
  fCalpha(track.fCalpha),
  fCx(track.fCx),
  fCchi2(track.fCchi2),
  fIalpha(track.fIalpha),
  fIx(track.fIx),
  fTalpha(track.fTalpha),
  fTx(track.fTx),
  fITSchi2(track.fITSchi2),
  fITSncls(track.fITSncls),
  fITSsignal(track.fITSsignal),
  fITSLabel(track.fITSLabel),
  fITSFakeRatio(track.fITSFakeRatio),
  fITStrack(0),    //coping separatelly - in user code
  fTPCchi2(track.fTPCchi2),
  fTPCncls(track.fTPCncls),
  fTPCClusterMap(track.fTPCClusterMap),
  fTPCsignal(track.fTPCsignal),
  fTPCLabel(track.fTPCLabel),
  fTRDchi2(track.fTRDchi2),
  fTRDncls(track.fTRDncls),
  fTRDncls0(track.fTRDncls0),
  fTRDsignal(track.fTRDsignal),
  fTRDLabel(track.fTRDLabel),
  fTRDQuality(track.fTRDQuality),
  fTRDtrack(0),
  fTOFchi2(track.fTOFchi2),
  fTOFindex(track.fTOFindex),
  fTOFsignal(track.fTOFsignal),
  fPHOSsignal(track.fPHOSsignal),
  fEMCALsignal(track.fEMCALsignal),
  fRICHchi2(track.fRICHchi2),
  fRICHncls(track.fRICHncls),
  fRICHindex(track.fRICHindex),
  fRICHsignal(track.fRICHsignal),
  fRICHtheta(track.fRICHtheta),
  fRICHphi(track.fRICHphi),
  fRICHdx(track.fRICHdx),
  fRICHdy(track.fRICHdy)
{
  //
  //copy constructor
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] =track.fTrackTime[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++)  fR[i] =track.fR[i];
  //
  for (Int_t i=0;i<5;i++) fRp[i] =track.fRp[i];
  for (Int_t i=0;i<15;i++) fRc[i] =track.fRc[i];
  //
  for (Int_t i=0;i<5;i++) fCp[i] =track.fCp[i];
  for (Int_t i=0;i<15;i++)  fCc[i] =track.fCc[i];
  //
  for (Int_t i=0;i<5;i++) fIp[i] =track.fIp[i];
  for (Int_t i=0;i<15;i++)  fIc[i] =track.fIc[i];
  //
  for (Int_t i=0;i<5;i++) fTp[i] =track.fTp[i];
  for (Int_t i=0;i<15;i++)  fTc[i] =track.fTc[i];
  //
  for (Int_t i=0;i<12;i++) fITSchi2MIP[i] =track.fITSchi2MIP[i];
  for (Int_t i=0;i<6;i++) fITSindex[i]=track.fITSindex[i];    
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]=track.fITSr[i]; 
  //
  for (Int_t i=0;i<180;i++) fTPCindex[i]=track.fTPCindex[i];  
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=track.fTPCr[i]; 
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=track.fTPCPoints[i];}
  for (Int_t i=0; i<3;i++)   { fKinkIndexes[i]=track.fKinkIndexes[i];}
  for (Int_t i=0; i<3;i++)   { fV0Indexes[i]=track.fV0Indexes[i];}
  //
  for (Int_t i=0;i<130;i++) fTRDindex[i]=track.fTRDindex[i];   
  for (Int_t i=0;i<kNPlane;i++) {
      fTRDsignals[i]=track.fTRDsignals[i]; 
      fTRDTimBin[i]=track.fTRDTimBin[i];
  }
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i]=track.fTRDr[i]; 
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i]=track.fTOFr[i];
  for (Int_t i=0;i<3;i++) fTOFLabel[i]=track.fTOFLabel[i];
  for (Int_t i=0;i<10;i++) fTOFInfo[i]=track.fTOFInfo[i];
  //
  for (Int_t i=0;i<3;i++) fPHOSpos[i]=track.fPHOSpos[i]; 
  for (Int_t i=0;i<AliPID::kSPECIESN;i++) fPHOSr[i]=track.fPHOSr[i]; 
  //
  for (Int_t i=0;i<3;i++) fEMCALpos[i]=track.fEMCALpos[i]; 
  for (Int_t i=0;i<AliPID::kSPECIESN;i++) fEMCALr[i]=track.fEMCALr[i]; 
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fRICHr[i]=track.fRICHr[i];
}
//_______________________________________________________________________
AliESDtrack::~AliESDtrack(){ 
  //
  // This is destructor according Coding Conventrions 
  //
  //printf("Delete track\n");
  delete fITStrack;
  delete fTRDtrack;  
}

//_______________________________________________________________________
void AliESDtrack::MakeMiniESDtrack(){
  // Resets everything except
  // fFlags: Reconstruction status flags 
  // fLabel: Track label
  // fID:  Unique ID of the track
  // fD: Impact parameter in XY-plane
  // fZ: Impact parameter in Z 
  // fR[AliPID::kSPECIES]: combined "detector response probability"
  // Running track parameters
  // fRalpha: track rotation angle
  // fRx: X-coordinate of the track reference plane 
  // fRp[5]: external track parameters  
  // fRc[15]: external cov. matrix of the track parameters
  
  fTrackLength = 0;
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = 0;
  fStopVertex = 0;

  // Reset track parameters constrained to the primary vertex
  fCalpha = 0;
  fCx = 0;
  for (Int_t i=0;i<5;i++) fCp[i] = 0;
  for (Int_t i=0;i<15;i++)  fCc[i] = 0;
  fCchi2 = 0;

  // Reset track parameters at the inner wall of TPC
  fIalpha = 0;
  fIx = 0;
  for (Int_t i=0;i<5;i++) fIp[i] = 0;
  for (Int_t i=0;i<15;i++)  fIc[i] = 0;

  // Reset track parameters at the inner wall of the TRD
  fTalpha = 0;
  fTx = 0;
  for (Int_t i=0;i<5;i++) fTp[i] = 0;
  for (Int_t i=0;i<15;i++)  fTc[i] = 0;

  // Reset ITS track related information
  fITSchi2 = 0;
  for (Int_t i=0;i<12;i++) fITSchi2MIP[i] = 0;
  fITSncls = 0;       
  for (Int_t i=0;i<6;i++) fITSindex[i]= 0;    
  fITSsignal = 0;     
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]= 0; 
  fITSLabel = 0;       
  fITSFakeRatio = 0;   
  fITStrack =0;

  // Reset TPC related track information
  fTPCchi2 = 0;       
  fTPCncls = 0;       
  for (Int_t i=0;i<180;i++) fTPCindex[i] = 0;  
  fTPCClusterMap = 0;  
  fTPCsignal= 0;      
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=0; 
  fTPCLabel=0;       
  for (Int_t i=0;i<4;i++) fTPCPoints[i] = 0;
  for (Int_t i=0; i<3;i++)   fKinkIndexes[i] = 0;
  for (Int_t i=0; i<3;i++)   fV0Indexes[i] = 0;

  // Reset TRD related track information
  fTRDchi2 = 0;        
  fTRDncls = 0;       
  fTRDncls0 = 0;       
  for (Int_t i=0;i<130;i++) fTRDindex[i] = 0;   
  fTRDsignal = 0;      
  for (Int_t i=0;i<kNPlane;i++) {
      fTRDsignals[i] = 0; 
      fTRDTimBin[i]  = 0;
  }
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i] = 0; 
  fTRDLabel = 0;       
  fTRDtrack = 0; 
  fTRDQuality  = 0;

  // Reset TOF related track information
  fTOFchi2 = 0;        
  fTOFindex = 0;       
  fTOFsignal = 0;      
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i] = 0;
  for (Int_t i=0;i<3;i++) fTOFLabel[i] = 0;
  for (Int_t i=0;i<10;i++) fTOFInfo[i] = 0;

  // Reset PHOS related track information
  for (Int_t i=0;i<3;i++) fPHOSpos[i] = 0; 
  fPHOSsignal = 0; 
  for (Int_t i=0;i<AliPID::kSPECIESN;i++) fPHOSr[i] = 0;
 
  // Reset EMCAL related track information
  for (Int_t i=0;i<3;i++) fEMCALpos[i] = 0; 
  fEMCALsignal = 0; 
  for (Int_t i=0;i<AliPID::kSPECIESN;i++) fEMCALr[i] = 0;
 
  // Reset RICH related track information
  fRICHchi2 = 0;     
  fRICHncls = 0;     
  fRICHindex = 0;     
  fRICHsignal = 0;     
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fRICHr[i] = 0;
  fRICHtheta = 0;     
  fRICHphi = 0;      
  fRICHdx = 0;     
  fRICHdy = 0;      

} 
//_______________________________________________________________________
Double_t AliESDtrack::GetMass() const {
  // Returns the mass of the most probable particle type
  Float_t max=0.;
  Int_t k=-1;
  for (Int_t i=0; i<AliPID::kSPECIES; i++) {
    if (fR[i]>max) {k=i; max=fR[i];}
  }
  if (k==0) { // dE/dx "crossing points" in the TPC
     Double_t p=GetP();
     if ((p>0.38)&&(p<0.48))
        if (fR[0]<fR[3]*10.) return AliPID::ParticleMass(AliPID::kKaon);
     if ((p>0.75)&&(p<0.85))
        if (fR[0]<fR[4]*10.) return AliPID::ParticleMass(AliPID::kProton);
     return 0.00051;
  }
  if (k==1) return AliPID::ParticleMass(AliPID::kMuon); 
  if (k==2||k==-1) return AliPID::ParticleMass(AliPID::kPion);
  if (k==3) return AliPID::ParticleMass(AliPID::kKaon);
  if (k==4) return AliPID::ParticleMass(AliPID::kProton);
  AliWarning("Undefined mass !");
  return AliPID::ParticleMass(AliPID::kPion);
}

//_______________________________________________________________________
Bool_t AliESDtrack::UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags) {
  //
  // This function updates track's running parameters 
  //
  Bool_t rc=kTRUE;

  SetStatus(flags);
  fLabel=t->GetLabel();

  if (t->IsStartedTimeIntegral()) {
    SetStatus(kTIME);
    Double_t times[10];t->GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t->GetIntegratedLength());
  }

  fRalpha=t->GetAlpha();
  t->GetExternalParameters(fRx,fRp);
  t->GetExternalCovariance(fRc);

  switch (flags) {
    
  case kITSin: case kITSout: case kITSrefit:
    fITSncls=t->GetNumberOfClusters();
    fITSchi2=t->GetChi2();
    for (Int_t i=0;i<fITSncls;i++) fITSindex[i]=t->GetClusterIndex(i);
    fITSsignal=t->GetPIDsignal();
    fITSLabel = t->GetLabel();
    fITSFakeRatio = t->GetFakeRatio();
    break;
    
  case kTPCin: case kTPCrefit:
    fTPCLabel = t->GetLabel();
    fIalpha=fRalpha;
    fIx=fRx;    
    {
      Int_t i;
      for (i=0; i<5; i++) fIp[i]=fRp[i];
      for (i=0; i<15;i++) fIc[i]=fRc[i];
    }
  case kTPCout:
  
    fTPCncls=t->GetNumberOfClusters();
    fTPCchi2=t->GetChi2();
    
     {//prevrow must be declared in separate namespace, otherwise compiler cries:
      //"jump to case label crosses initialization of `Int_t prevrow'"
       Int_t prevrow = -1;
       //       for (Int_t i=0;i<fTPCncls;i++) 
       for (Int_t i=0;i<160;i++) 
        {
          fTPCindex[i]=t->GetClusterIndex(i);

          // Piotr's Cluster Map for HBT  
          // ### please change accordingly if cluster array is changing 
          // to "New TPC Tracking" style (with gaps in array) 
          Int_t idx = fTPCindex[i];
          Int_t sect = (idx&0xff000000)>>24;
          Int_t row = (idx&0x00ff0000)>>16;
          if (sect > 18) row +=63; //if it is outer sector, add number of inner sectors

          fTPCClusterMap.SetBitNumber(row,kTRUE);

          //Fill the gap between previous row and this row with 0 bits
          //In case  ###  pleas change it as well - just set bit 0 in case there 
          //is no associated clusters for current "i"
          if (prevrow < 0) 
           {
             prevrow = row;//if previous bit was not assigned yet == this is the first one
           }
          else
           { //we don't know the order (inner to outer or reverse)
             //just to be save in case it is going to change
             Int_t n = 0, m = 0;
             if (prevrow < row)
              {
                n = prevrow;
                m = row;
              }
             else
              {
                n = row;
                m = prevrow;
              }

             for (Int_t j = n+1; j < m; j++)
              {
                fTPCClusterMap.SetBitNumber(j,kFALSE);
              }
             prevrow = row; 
           }
          // End Of Piotr's Cluster Map for HBT
        }
     }
    fTPCsignal=t->GetPIDsignal();
    {Double_t mass=t->GetMass();    // preliminary mass setting 
    if (mass>0.5) fR[4]=1.;         //        used by
    else if (mass<0.4) fR[2]=1.;    // the ITS reconstruction
    else fR[3]=1.;}
                     //
    break;

  case kTRDout: case kTRDin: case kTRDrefit:
    fTRDLabel = t->GetLabel(); 
    fTRDncls=t->GetNumberOfClusters();
    fTRDchi2=t->GetChi2();
    for (Int_t i=0;i<fTRDncls;i++) fTRDindex[i]=t->GetClusterIndex(i);
    fTRDsignal=t->GetPIDsignal();
    break;
  case kTRDbackup:
    t->GetExternalParameters(fTx,fTp);
    t->GetExternalCovariance(fTc);
    fTalpha = t->GetAlpha();
    fTRDncls0 = t->GetNumberOfClusters(); 
    break;
  case kTOFin: 
    break;
  case kTOFout: 
    break;
  case kTRDStop:
    break;
  default: 
    AliError("Wrong flag !");
    return kFALSE;
  }

  return rc;
}

//_______________________________________________________________________
void 
AliESDtrack::SetConstrainedTrackParams(const AliKalmanTrack *t, Double_t chi2) {
  //
  // This function sets the constrained track parameters 
  //
  Int_t i;
  Double_t x,buf[15];
  fCalpha=t->GetAlpha();
  t->GetExternalParameters(x,buf); fCx=x;
  for (i=0; i<5; i++) fCp[i]=buf[i];
  t->GetExternalCovariance(buf);
  for (i=0; i<15; i++) fCc[i]=buf[i];
  fCchi2=chi2;
}


//_______________________________________________________________________
void AliESDtrack::GetExternalParameters(Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  //---------------------------------------------------------------------
  x=fRx;
  for (Int_t i=0; i<5; i++) p[i]=fRp[i];
}

//_______________________________________________________________________
Bool_t AliESDtrack::GetExternalParametersAt(Double_t x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  // at the position given by the first argument 
  //---------------------------------------------------------------------
  Double_t dx=x-fRx;
  Double_t f1=fRp[2], f2=f1 + dx*fRp[4]/AliKalmanTrack::GetConvConst();

  if (TMath::Abs(f2) >= 0.9999) return kFALSE;
  
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);
  p[0] = fRp[0] + dx*(f1+f2)/(r1+r2);
  p[1] = fRp[1] + dx*(f1+f2)/(f1*r2 + f2*r1)*fRp[3];
  p[2] = f2;
  p[3] = fRp[3];
  p[4] = fRp[4];

  return kTRUE;
}

//_______________________________________________________________________
void AliESDtrack::GetExternalCovariance(Double_t cov[15]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the cov. matrix
  //---------------------------------------------------------------------
  for (Int_t i=0; i<15; i++) cov[i]=fRc[i];
}


//_______________________________________________________________________
void 
AliESDtrack::GetConstrainedExternalParameters(Double_t &x, Double_t p[5])const{
  //---------------------------------------------------------------------
  // This function returns the constrained external track parameters
  //---------------------------------------------------------------------
  x=fCx;
  for (Int_t i=0; i<5; i++) p[i]=fCp[i];
}
//_______________________________________________________________________
void 
AliESDtrack::GetConstrainedExternalCovariance(Double_t c[15]) const {
  //---------------------------------------------------------------------
  // This function returns the constrained external cov. matrix
  //---------------------------------------------------------------------
  for (Int_t i=0; i<15; i++) c[i]=fCc[i];
}


Double_t AliESDtrack::GetP() const {
  //---------------------------------------------------------------------
  // This function returns the track momentum
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fRp[4])<=0) return 0;
  Double_t pt=1./TMath::Abs(fRp[4]);
  return pt*TMath::Sqrt(1.+ fRp[3]*fRp[3]);
}

Bool_t Local2GlobalMomentum(Double_t p[3],Double_t alpha) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track momentum.
  // When called, the arguments are:
  //    p[0] = 1/pt of the track;
  //    p[1] = sine of local azim. angle of the track momentum;
  //    p[2] = tangent of the track momentum dip angle;
  //   alpha - rotation angle. 
  // The result is returned as:
  //    p[0] = px
  //    p[1] = py
  //    p[2] = pz
  // Results for (nearly) straight tracks are meaningless !
  //----------------------------------------------------------------
  if (TMath::Abs(p[0])<=0)        return kFALSE;
  if (TMath::Abs(p[1])> 0.999999) return kFALSE;

  Double_t pt=1./TMath::Abs(p[0]);
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha);
  Double_t r=TMath::Sqrt(1 - p[1]*p[1]);
  p[0]=pt*(r*cs - p[1]*sn); p[1]=pt*(p[1]*cs + r*sn); p[2]=pt*p[2];

  return kTRUE;
}

Bool_t Local2GlobalPosition(Double_t r[3],Double_t alpha) {
  //----------------------------------------------------------------
  // This function performs local->global transformation of the
  // track position.
  // When called, the arguments are:
  //    r[0] = local x
  //    r[1] = local y
  //    r[2] = local z
  //   alpha - rotation angle. 
  // The result is returned as:
  //    r[0] = global x
  //    r[1] = global y
  //    r[2] = global z
  //----------------------------------------------------------------
  Double_t cs=TMath::Cos(alpha), sn=TMath::Sin(alpha), x=r[0];
  r[0]=x*cs - r[1]*sn; r[1]=x*sn + r[1]*cs;

  return kTRUE;
}

Bool_t AliESDtrack::GetConstrainedPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the constrained global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  p[0]=fCp[4]; p[1]=fCp[2]; p[2]=fCp[3];
  return Local2GlobalMomentum(p,fCalpha);
}  

Bool_t AliESDtrack::GetConstrainedXYZ(Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the constrained global track position
  //---------------------------------------------------------------------
  r[0]=fCx; r[1]=fCp[0]; r[2]=fCp[1];
  return Local2GlobalPosition(r,fCalpha);
}

Bool_t AliESDtrack::GetPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  p[0]=fRp[4]; p[1]=fRp[2]; p[2]=fRp[3];
  return Local2GlobalMomentum(p,fRalpha);
}

Bool_t AliESDtrack::GetXYZ(Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  r[0]=fRx; r[1]=fRp[0]; r[2]=fRp[1];
  return Local2GlobalPosition(r,fRalpha);
}

void AliESDtrack::GetCovariance(Double_t cv[21]) const {
  //---------------------------------------------------------------------
  // This function returns the global covariance matrix of the track params
  // 
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,x)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,x)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
  //
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fRp[4])<=0) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return;
  }
  if (TMath::Abs(fRp[2]) > 0.999999) {
     for (Int_t i=0; i<21; i++) cv[i]=0.;
     return;
  }
  Double_t pt=1./TMath::Abs(fRp[4]);
  Double_t cs=TMath::Cos(fRalpha), sn=TMath::Sin(fRalpha);
  Double_t r=TMath::Sqrt(1-fRp[2]*fRp[2]);

  Double_t m00=-sn, m10=cs;
  Double_t m23=-pt*(sn + fRp[2]*cs/r), m43=-pt*pt*(r*cs - fRp[2]*sn);
  Double_t m24= pt*(cs - fRp[2]*sn/r), m44=-pt*pt*(r*sn + fRp[2]*cs);
  Double_t m35=pt, m45=-pt*pt*fRp[3];

  cv[0]=fRc[0]*m00*m00;
  cv[1]=fRc[0]*m00*m10; 
  cv[2]=fRc[0]*m10*m10;
  cv[3]=fRc[1]*m00; 
  cv[4]=fRc[1]*m10; 
  cv[5]=fRc[2];
  cv[6]=m00*(fRc[3]*m23+fRc[10]*m43); 
  cv[7]=m10*(fRc[3]*m23+fRc[10]*m43); 
  cv[8]=fRc[4]*m23+fRc[11]*m43; 
  cv[9]=m23*(fRc[5]*m23+fRc[12]*m43)+m43*(fRc[12]*m23+fRc[14]*m43);
  cv[10]=m00*(fRc[3]*m24+fRc[10]*m44); 
  cv[11]=m10*(fRc[3]*m24+fRc[10]*m44); 
  cv[12]=fRc[4]*m24+fRc[11]*m44; 
  cv[13]=m23*(fRc[5]*m24+fRc[12]*m44)+m43*(fRc[12]*m24+fRc[14]*m44);
  cv[14]=m24*(fRc[5]*m24+fRc[12]*m44)+m44*(fRc[12]*m24+fRc[14]*m44);
  cv[15]=m00*(fRc[6]*m35+fRc[10]*m45); 
  cv[16]=m10*(fRc[6]*m35+fRc[10]*m45); 
  cv[17]=fRc[7]*m35+fRc[11]*m45; 
  cv[18]=m23*(fRc[8]*m35+fRc[12]*m45)+m43*(fRc[13]*m35+fRc[14]*m45);
  cv[19]=m24*(fRc[8]*m35+fRc[12]*m45)+m44*(fRc[13]*m35+fRc[14]*m45); 
  cv[20]=m35*(fRc[9]*m35+fRc[13]*m45)+m45*(fRc[13]*m35+fRc[14]*m45);
}

Bool_t AliESDtrack::GetInnerPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // af the entrance of the TPC
  //---------------------------------------------------------------------
  p[0]=fIp[4]; p[1]=fIp[2]; p[2]=fIp[3];
  return Local2GlobalMomentum(p,fIalpha);
}

Bool_t AliESDtrack::GetInnerXYZ(Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  // af the entrance of the TPC
  //---------------------------------------------------------------------
  if (fIx==0) return kFALSE;
  r[0]=fIx; r[1]=fIp[0]; r[2]=fIp[1];
  return Local2GlobalPosition(r,fIalpha);
}

void AliESDtrack::GetInnerExternalParameters(Double_t &x, Double_t p[5]) const 
{
  //skowron
 //---------------------------------------------------------------------
  // This function returns external representation of the track parameters at Inner Layer of TPC
  //---------------------------------------------------------------------
  x=fIx;
  for (Int_t i=0; i<5; i++) p[i]=fIp[i];
}
void AliESDtrack::GetInnerExternalCovariance(Double_t cov[15]) const
{
 //skowron
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix at Inner Layer of TPC
 //---------------------------------------------------------------------
 for (Int_t i=0; i<15; i++) cov[i]=fIc[i];
 
}

void  AliESDtrack::GetTRDExternalParameters(Double_t &x, Double_t&alpha, Double_t p[5], Double_t cov[15]) const
{
  //
  //this function returns TRD parameters
  //
  x=fTx;
  alpha = fTalpha; 
  for (Int_t i=0; i<5; i++) p[i]=fTp[i];
  for (Int_t i=0; i<15; i++) cov[i]=fTc[i];
}

Bool_t AliESDtrack::GetPxPyPzAt(Double_t x,Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // at the position "x" using the helix track approximation
  //---------------------------------------------------------------------
  p[0]=fRp[4]; 
  p[1]=fRp[2]+(x-fRx)*fRp[4]/AliKalmanTrack::GetConvConst(); 
  p[2]=fRp[3];
  return Local2GlobalMomentum(p,fRalpha);
}

Bool_t AliESDtrack::GetXYZAt(Double_t x, Double_t *r) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  // af the radius "x" using the helix track approximation
  //---------------------------------------------------------------------
  Double_t dx=x-fRx;
  Double_t f1=fRp[2], f2=f1 + dx*fRp[4]/AliKalmanTrack::GetConvConst();

  if (TMath::Abs(f2) >= 0.9999) return kFALSE;
  
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);
  r[0] = x;
  r[1] = fRp[0] + dx*(f1+f2)/(r1+r2);
  r[2] = fRp[1] + dx*(f1+f2)/(f1*r2 + f2*r1)*fRp[3];
  return Local2GlobalPosition(r,fRalpha);
}

//_______________________________________________________________________
void AliESDtrack::GetIntegratedTimes(Double_t *times) const {
  // Returns the array with integrated times for each particle hypothesis
  for (Int_t i=0; i<AliPID::kSPECIES; i++) times[i]=fTrackTime[i];
}

//_______________________________________________________________________
void AliESDtrack::SetIntegratedTimes(const Double_t *times) {
  // Sets the array with integrated times for each particle hypotesis
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fTrackTime[i]=times[i];
}

//_______________________________________________________________________
void AliESDtrack::SetITSpid(const Double_t *p) {
  // Sets values for the probability of each particle type (in ITS)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fITSr[i]=p[i];
  SetStatus(AliESDtrack::kITSpid);
}

void AliESDtrack::SetITSChi2MIP(const Float_t *chi2mip){
  for (Int_t i=0; i<12; i++) fITSchi2MIP[i]=chi2mip[i];
}
//_______________________________________________________________________
void AliESDtrack::GetITSpid(Double_t *p) const {
  // Gets the probability of each particle type (in ITS)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fITSr[i];
}

//_______________________________________________________________________
Int_t AliESDtrack::GetITSclusters(UInt_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  for (Int_t i=0; i<fITSncls; i++) idx[i]=fITSindex[i];
  return fITSncls;
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTPCclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  if (idx!=0)
    for (Int_t i=0; i<180; i++) idx[i]=fTPCindex[i];  // MI I prefer some constant
  return fTPCncls;
}

//_______________________________________________________________________
void AliESDtrack::SetTPCpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TPC)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fTPCr[i]=p[i];
  SetStatus(AliESDtrack::kTPCpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTPCpid(Double_t *p) const {
  // Gets the probability of each particle type (in TPC)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTPCr[i];
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTRDclusters(UInt_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined TRD clusters 
  //---------------------------------------------------------------------
  if (idx!=0)
    for (Int_t i=0; i<130; i++) idx[i]=fTRDindex[i];  // MI I prefer some constant
  return fTRDncls;
}

//_______________________________________________________________________
void AliESDtrack::SetTRDpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TRD)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fTRDr[i]=p[i];
  SetStatus(AliESDtrack::kTRDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTRDpid(Double_t *p) const {
  // Gets the probability of each particle type (in TRD)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTRDr[i];
}

//_______________________________________________________________________
void    AliESDtrack::SetTRDpid(Int_t iSpecies, Float_t p)
{
  // Sets the probability of particle type iSpecies to p (in TRD)
  fTRDr[iSpecies] = p;
}

Float_t AliESDtrack::GetTRDpid(Int_t iSpecies) const
{
  // Returns the probability of particle type iSpecies (in TRD)
  return fTRDr[iSpecies];
}

//_______________________________________________________________________
void AliESDtrack::SetTOFpid(const Double_t *p) {  
  // Sets the probability of each particle type (in TOF)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fTOFr[i]=p[i];
  SetStatus(AliESDtrack::kTOFpid);
}

//_______________________________________________________________________
void AliESDtrack::SetTOFLabel(const Int_t *p) {  
  // Sets  (in TOF)
  for (Int_t i=0; i<3; i++) fTOFLabel[i]=p[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFpid(Double_t *p) const {
  // Gets probabilities of each particle type (in TOF)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTOFr[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFLabel(Int_t *p) const {
  // Gets (in TOF)
  for (Int_t i=0; i<3; i++) p[i]=fTOFLabel[i];
}

//_______________________________________________________________________
void AliESDtrack::GetTOFInfo(Float_t *info) const {
  // Gets (in TOF)
  for (Int_t i=0; i<10; i++) info[i]=fTOFInfo[i];
}

//_______________________________________________________________________
void AliESDtrack::SetTOFInfo(Float_t*info) {
  // Gets (in TOF)
  for (Int_t i=0; i<10; i++) fTOFInfo[i]=info[i];
}



//_______________________________________________________________________
void AliESDtrack::SetPHOSpid(const Double_t *p) {  
  // Sets the probability of each particle type (in PHOS)
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) fPHOSr[i]=p[i];
  SetStatus(AliESDtrack::kPHOSpid);
}

//_______________________________________________________________________
void AliESDtrack::GetPHOSpid(Double_t *p) const {
  // Gets probabilities of each particle type (in PHOS)
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) p[i]=fPHOSr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetEMCALpid(const Double_t *p) {  
  // Sets the probability of each particle type (in EMCAL)
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) fEMCALr[i]=p[i];
  SetStatus(AliESDtrack::kEMCALpid);
}

//_______________________________________________________________________
void AliESDtrack::GetEMCALpid(Double_t *p) const {
  // Gets probabilities of each particle type (in EMCAL)
  for (Int_t i=0; i<AliPID::kSPECIESN; i++) p[i]=fEMCALr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetRICHpid(const Double_t *p) {  
  // Sets the probability of each particle type (in RICH)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fRICHr[i]=p[i];
  SetStatus(AliESDtrack::kRICHpid);
}

//_______________________________________________________________________
void AliESDtrack::GetRICHpid(Double_t *p) const {
  // Gets probabilities of each particle type (in RICH)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fRICHr[i];
}



//_______________________________________________________________________
void AliESDtrack::SetESDpid(const Double_t *p) {  
  // Sets the probability of each particle type for the ESD track
  for (Int_t i=0; i<AliPID::kSPECIES; i++) fR[i]=p[i];
  SetStatus(AliESDtrack::kESDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetESDpid(Double_t *p) const {
  // Gets probability of each particle type for the ESD track
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fR[i];
}

//_______________________________________________________________________
void AliESDtrack::Print(Option_t *) const {
  // Prints info on the track
  
  printf("ESD track info\n") ; 
  Double_t p[AliPID::kSPECIESN] ; 
  Int_t index = 0 ; 
  if( IsOn(kITSpid) ){
    printf("From ITS: ") ; 
    GetITSpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetITSsignal()) ;
  } 
  if( IsOn(kTPCpid) ){
    printf("From TPC: ") ; 
    GetTPCpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTPCsignal()) ;
  }
  if( IsOn(kTRDpid) ){
    printf("From TRD: ") ; 
    GetTRDpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTRDsignal()) ;
  }
  if( IsOn(kTOFpid) ){
    printf("From TOF: ") ; 
    GetTOFpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTOFsignal()) ;
  }
  if( IsOn(kRICHpid) ){
    printf("From TOF: ") ; 
    GetRICHpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetRICHsignal()) ;
  }
  if( IsOn(kPHOSpid) ){
    printf("From PHOS: ") ; 
    GetPHOSpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIESN; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetPHOSsignal()) ;
  }
  if( IsOn(kEMCALpid) ){
    printf("From EMCAL: ") ; 
    GetEMCALpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIESN; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetEMCALsignal()) ;
  }
} 
