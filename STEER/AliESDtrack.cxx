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
fFlags(0),
fLabel(0),
fTrackLength(0),
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
fOalpha(0),
fOx(0),
fITSchi2(0),
fITSncls(0),
fITSsignal(0),
fTPCchi2(0),
fTPCncls(0),
fTPCClusterMap(159),//number of padrows
fTPCsignal(0),
fTRDchi2(0),
fTRDncls(0),
fTRDncls0(0),
fTRDsignal(0),
fTOFchi2(0),
fTOFindex(0),
fTOFsignal(-1),
fPHOSsignal(-1),
fEMCALsignal(-1),
fRICHsignal(-1)
{
  //
  // The default ESD constructor 
  //
  fID =0;
  for (Int_t i=0; i<kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=1.;
    fITSr[i]=1.;
    fTPCr[i]=1.;
    fTRDr[i]=1.;
    fTOFr[i]=1.;
    fRICHr[i]=1.;
  }
  
  for (Int_t i=0; i<kSPECIESN; i++) {
    fPHOSr[i]  = 1.;
    fEMCALr[i] = 1.;
  }

 
  fPHOSpos[0]=fPHOSpos[1]=fPHOSpos[2]=0.;
  fEMCALpos[0]=fEMCALpos[1]=fEMCALpos[2]=0.;
  Int_t i;
  for (i=0; i<5; i++)  { 
    fRp[i]=fCp[i]=fIp[i]=fOp[i]=fXp[i]=fTp[i]=0.;
  }
  for (i=0; i<15; i++) { 
    fRc[i]=fCc[i]=fIc[i]=fOc[i]=fXc[i]=fTc[i]=0.;  
  }
  for (i=0; i<6; i++)  { fITSindex[i]=0; }
  for (i=0; i<180; i++){ fTPCindex[i]=0; }
  for (i=0; i<3;i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3;i++)   { fV0Indexes[i]=-1;}
  for (i=0; i<130; i++) { fTRDindex[i]=0; }
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=-1;}
  for (Int_t i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (Int_t i=0;i<10;i++) {fTOFInfo[i]=-1;}
  fTPCLabel = 0;
  fTRDLabel = 0;
  fITSLabel = 0;
  fITStrack = 0;
  fTRDtrack = 0;  
}

//_______________________________________________________________________

AliESDtrack::AliESDtrack(const AliESDtrack& track):TObject(track){
  //
  //copy constructor
  //
  fID = track.fID;
  fFlags = track.fFlags;
  fLabel =track.fLabel;
  fTrackLength =track.fTrackLength;
  for (Int_t i=0;i<kSPECIES;i++) fTrackTime[i] =track.fTrackTime[i];
  for (Int_t i=0;i<kSPECIES;i++)  fR[i] =track.fR[i];
  fStopVertex =track.fStopVertex;
  //
  fRalpha =track.fRalpha;
  fRx =track.fRx;
  for (Int_t i=0;i<5;i++) fRp[i] =track.fRp[i];
  for (Int_t i=0;i<15;i++) fRc[i] =track.fRc[i];
  //
  fCalpha =track.fCalpha;
  fCx =track.fCx;
  for (Int_t i=0;i<5;i++) fCp[i] =track.fCp[i];
  for (Int_t i=0;i<15;i++)  fCc[i] =track.fCc[i];
  fCchi2 =track.fCchi2;
  //
  fIalpha =track.fIalpha;
  fIx =track.fIx;
  for (Int_t i=0;i<5;i++) fIp[i] =track.fIp[i];
  for (Int_t i=0;i<15;i++)  fIc[i] =track.fIc[i];
  //
  fTalpha =track.fTalpha;
  fTx =track.fTx;
  for (Int_t i=0;i<5;i++) fTp[i] =track.fTp[i];
  for (Int_t i=0;i<15;i++)  fTc[i] =track.fTc[i];
  //
  fOalpha =track.fOalpha;
  fOx =track.fOx;
  for (Int_t i=0;i<5;i++) fOp[i] =track.fOp[i];
  for (Int_t i=0;i<15;i++)  fOc[i] =track.fOc[i];
  //
  fXalpha =track.fXalpha;
  fXx =track.fXx;
  for (Int_t i=0;i<5;i++) fXp[i] =track.fXp[i];
  for (Int_t i=0;i<15;i++) fXc[i] =track.fXc[i];
  //
  fITSchi2 =track.fITSchi2;
  for (Int_t i=0;i<12;i++) fITSchi2MIP[i] =track.fITSchi2MIP[i];
  fITSncls =track.fITSncls;       
  for (Int_t i=0;i<6;i++) fITSindex[i]=track.fITSindex[i];    
  fITSsignal =track.fITSsignal;     
  for (Int_t i=0;i<kSPECIES;i++) fITSr[i]=track.fITSr[i]; 
  fITSLabel =track.fITSLabel;       
  fITSFakeRatio =track.fITSFakeRatio;   
  fITStrack =0;  //coping separatelly - in user code
  //
  fTPCchi2 =track.fTPCchi2;       
  fTPCncls =track.fTPCncls;       
  for (Int_t i=0;i<180;i++) fTPCindex[i]=track.fTPCindex[i];  
  fTPCClusterMap=track.fTPCClusterMap;  
  fTPCsignal=track.fTPCsignal;      
  for (Int_t i=0;i<kSPECIES;i++) fTPCr[i]=track.fTPCr[i]; 
  fTPCLabel=track.fTPCLabel;       
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=track.fTPCPoints[i];}
  for (Int_t i=0; i<3;i++)   { fKinkIndexes[i]=track.fKinkIndexes[i];}
  for (Int_t i=0; i<3;i++)   { fV0Indexes[i]=track.fV0Indexes[i];}
  //
  fTRDchi2=track.fTRDchi2;        
  fTRDncls=track.fTRDncls;       
  fTRDncls0=track.fTRDncls0;       
  for (Int_t i=0;i<130;i++) fTRDindex[i]=track.fTRDindex[i];   
  fTRDsignal=track.fTRDsignal;      
  for (Int_t i=0;i<kSPECIES;i++) fTRDr[i]=track.fTRDr[i]; 
  fTRDLabel=track.fTRDLabel;       
  fTRDtrack=0; 
  //
  fTOFchi2=track.fTOFchi2;        
  fTOFindex=track.fTOFindex;       
  fTOFsignal=track.fTOFsignal;      
  for (Int_t i=0;i<kSPECIES;i++) fTOFr[i]=track.fTOFr[i];
  for (Int_t i=0;i<3;i++) fTOFLabel[i]=track.fTOFLabel[i];
  for (Int_t i=0;i<10;i++) fTOFInfo[i]=track.fTOFInfo[i];
  //
  for (Int_t i=0;i<3;i++) fPHOSpos[i]=track.fPHOSpos[i]; 
  fPHOSsignal=track.fPHOSsignal; 
  for (Int_t i=0;i<kSPECIESN;i++) fPHOSr[i]=track.fPHOSr[i]; 
  //
  for (Int_t i=0;i<3;i++) fEMCALpos[i]=track.fEMCALpos[i]; 
  fEMCALsignal=track.fEMCALsignal; 
  for (Int_t i=0;i<kSPECIESN;i++) fEMCALr[i]=track.fEMCALr[i]; 
  //
  fRICHsignal=track.fRICHsignal;     
  for (Int_t i=0;i<kSPECIES;i++) fRICHr[i]=track.fRICHr[i];
  
  
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
Double_t AliESDtrack::GetMass() const {
  // Returns the mass of the most probable particle type
  Float_t max=0.;
  Int_t k=-1;
  for (Int_t i=0; i<kSPECIES; i++) {
    if (fR[i]>max) {k=i; max=fR[i];}
  }
  if (k==0) { // dE/dx "crossing points" in the TPC
     Double_t p=GetP();
     if ((p>0.38)&&(p<0.48))
        if (fR[0]<fR[3]*10.) return 0.49368;
     if ((p>0.75)&&(p<0.85))
        if (fR[0]<fR[4]*10.) return 0.93827;
     return 0.00051;
  }
  if (k==1) return 0.10566; 
  if (k==2||k==-1) return 0.13957;
  if (k==3) return 0.49368;
  if (k==4) return 0.93827;
  AliWarning("Undefined mass !");
  return 0.13957;
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

  case kTRDout:
    //requested by the PHOS/EMCAL  ("temporary solution")
    if (GetExternalParametersAt(460.,fOp)) {
       fOalpha=t->GetAlpha();
       fOx=460.;
       t->GetExternalCovariance(fOc); //can be done better
    }
    if (GetExternalParametersAt(450.,fXp)) {
       fXalpha=t->GetAlpha();
       fXx=450.;
       t->GetExternalCovariance(fXc); //can be done better
    }
  case kTRDin: case kTRDrefit:
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

void AliESDtrack::GetConstrainedPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the constrained global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fCp[4])<=0) {
    p[0]=p[1]=p[2]=0;
    return;
  }
  if (TMath::Abs(fCp[2]) > 0.999999) {
     p[0]=p[1]=p[2]=0;
     return;
  }
  Double_t pt=1./TMath::Abs(fCp[4]);
  Double_t cs=TMath::Cos(fCalpha), sn=TMath::Sin(fCalpha);
  Double_t r=TMath::Sqrt(1-fCp[2]*fCp[2]);
  p[0]=pt*(r*cs - fCp[2]*sn); p[1]=pt*(fCp[2]*cs + r*sn); p[2]=pt*fCp[3];
}

void AliESDtrack::GetConstrainedXYZ(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fCalpha), sn=TMath::Sin(fCalpha);
  xyz[0]=fCx*cs - fCp[0]*sn; xyz[1]=fCx*sn + fCp[0]*cs; xyz[2]=fCp[1];
}

void AliESDtrack::GetPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // Results for (nearly) straight tracks are meaningless !
  //---------------------------------------------------------------------
  if (TMath::Abs(fRp[4])<=0) {
     p[0]=p[1]=p[2]=0;
     return;
  }
  if (TMath::Abs(fRp[2]) > 0.999999) {
     p[0]=p[1]=p[2]=0;
     return;
  }
  Double_t pt=1./TMath::Abs(fRp[4]);
  Double_t cs=TMath::Cos(fRalpha), sn=TMath::Sin(fRalpha);
  Double_t r=TMath::Sqrt(1-fRp[2]*fRp[2]);
  p[0]=pt*(r*cs - fRp[2]*sn); p[1]=pt*(fRp[2]*cs + r*sn); p[2]=pt*fRp[3];
}

void AliESDtrack::GetXYZ(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  //---------------------------------------------------------------------
  Double_t cs=TMath::Cos(fRalpha), sn=TMath::Sin(fRalpha);
  xyz[0]=fRx*cs - fRp[0]*sn; xyz[1]=fRx*sn + fRp[0]*cs; xyz[2]=fRp[1];
}

void AliESDtrack::GetCovariance(Double_t cv[21]) const {
  //---------------------------------------------------------------------
  // This function returns the global covariance matrix of the track params
  // 
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,y)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,z)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
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

void AliESDtrack::GetInnerPxPyPz(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // af the entrance of the TPC
  //---------------------------------------------------------------------
  if (fIx==0) {p[0]=p[1]=p[2]=0.; return;}
  Double_t phi=TMath::ASin(fIp[2]) + fIalpha;
  Double_t pt=1./TMath::Abs(fIp[4]);
  p[0]=pt*TMath::Cos(phi); p[1]=pt*TMath::Sin(phi); p[2]=pt*fIp[3]; 
}

void AliESDtrack::GetInnerXYZ(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  // af the entrance of the TPC
  //---------------------------------------------------------------------
  if (fIx==0) {xyz[0]=xyz[1]=xyz[2]=0.; return;}
  Double_t phi=TMath::ATan2(fIp[0],fIx) + fIalpha;
  Double_t r=TMath::Sqrt(fIx*fIx + fIp[0]*fIp[0]);
  xyz[0]=r*TMath::Cos(phi); xyz[1]=r*TMath::Sin(phi); xyz[2]=fIp[1]; 
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

void AliESDtrack::GetOuterPxPyPzPHOS(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // af the radius of the PHOS
  //---------------------------------------------------------------------
  p[0]=p[1]=p[2]=0. ; 
  if (fOx==0) 
    return;
  Double_t phi=TMath::ASin(fOp[2]) + fOalpha;
  Double_t pt=1./TMath::Abs(fOp[4]);
  p[0]=pt*TMath::Cos(phi); 
  p[1]=pt*TMath::Sin(phi); 
  p[2]=pt*fOp[3];
} 
void AliESDtrack::GetOuterPxPyPzEMCAL(Double_t *p) const {
  //---------------------------------------------------------------------
  // This function returns the global track momentum components
  // af the radius of the EMCAL
  //---------------------------------------------------------------------
  if (fXx==0)
    return;
  Double_t phi=TMath::ASin(fXp[2]) + fXalpha;
  Double_t pt=1./TMath::Abs(fXp[4]);
  p[0]=pt*TMath::Cos(phi); 
  p[1]=pt*TMath::Sin(phi); 
  p[2]=pt*fXp[3];
}

void AliESDtrack::GetOuterXYZPHOS(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  // af the radius of the PHOS
  //---------------------------------------------------------------------
  xyz[0]=xyz[1]=xyz[2]=0.;
  if (fOx==0) 
    return;
  Double_t phi=TMath::ATan2(fOp[0],fOx) + fOalpha;
  Double_t r=TMath::Sqrt(fOx*fOx + fOp[0]*fOp[0]);
  xyz[0]=r*TMath::Cos(phi); xyz[1]=r*TMath::Sin(phi); xyz[2]=fOp[1]; 
} 
void AliESDtrack::GetOuterXYZEMCAL(Double_t *xyz) const {
  //---------------------------------------------------------------------
  // This function returns the global track position
  // af the radius of the EMCAL
  //---------------------------------------------------------------------
  if (fXx==0) 
    return;
  Double_t phi=TMath::ATan2(fXp[0],fOx) + fXalpha;
  Double_t r=TMath::Sqrt(fXx*fXx + fXp[0]*fXp[0]);
  xyz[0]=r*TMath::Cos(phi); 
  xyz[1]=r*TMath::Sin(phi); 
  xyz[2]=fXp[1]; 
} 

//_______________________________________________________________________
void AliESDtrack::GetIntegratedTimes(Double_t *times) const {
  // Returns the array with integrated times for each particle hypothesis
  for (Int_t i=0; i<kSPECIES; i++) times[i]=fTrackTime[i];
}

//_______________________________________________________________________
void AliESDtrack::SetIntegratedTimes(const Double_t *times) {
  // Sets the array with integrated times for each particle hypotesis
  for (Int_t i=0; i<kSPECIES; i++) fTrackTime[i]=times[i];
}

//_______________________________________________________________________
void AliESDtrack::SetITSpid(const Double_t *p) {
  // Sets values for the probability of each particle type (in ITS)
  for (Int_t i=0; i<kSPECIES; i++) fITSr[i]=p[i];
  SetStatus(AliESDtrack::kITSpid);
}

void AliESDtrack::SetITSChi2MIP(const Float_t *chi2mip){
  for (Int_t i=0; i<12; i++) fITSchi2MIP[i]=chi2mip[i];
}
//_______________________________________________________________________
void AliESDtrack::GetITSpid(Double_t *p) const {
  // Gets the probability of each particle type (in ITS)
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fITSr[i];
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
  for (Int_t i=0; i<kSPECIES; i++) fTPCr[i]=p[i];
  SetStatus(AliESDtrack::kTPCpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTPCpid(Double_t *p) const {
  // Gets the probability of each particle type (in TPC)
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTPCr[i];
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
  for (Int_t i=0; i<kSPECIES; i++) fTRDr[i]=p[i];
  SetStatus(AliESDtrack::kTRDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTRDpid(Double_t *p) const {
  // Gets the probability of each particle type (in TRD)
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTRDr[i];
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
  for (Int_t i=0; i<kSPECIES; i++) fTOFr[i]=p[i];
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
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fTOFr[i];
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
  for (Int_t i=0; i<kSPECIESN; i++) fPHOSr[i]=p[i];
  SetStatus(AliESDtrack::kPHOSpid);
}

//_______________________________________________________________________
void AliESDtrack::GetPHOSpid(Double_t *p) const {
  // Gets probabilities of each particle type (in PHOS)
  for (Int_t i=0; i<kSPECIESN; i++) p[i]=fPHOSr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetEMCALpid(const Double_t *p) {  
  // Sets the probability of each particle type (in EMCAL)
  for (Int_t i=0; i<kSPECIESN; i++) fEMCALr[i]=p[i];
  SetStatus(AliESDtrack::kEMCALpid);
}

//_______________________________________________________________________
void AliESDtrack::GetEMCALpid(Double_t *p) const {
  // Gets probabilities of each particle type (in EMCAL)
  for (Int_t i=0; i<kSPECIESN; i++) p[i]=fEMCALr[i];
}

//_______________________________________________________________________
void AliESDtrack::SetRICHpid(const Double_t *p) {  
  // Sets the probability of each particle type (in RICH)
  for (Int_t i=0; i<kSPECIES; i++) fRICHr[i]=p[i];
  SetStatus(AliESDtrack::kRICHpid);
}

//_______________________________________________________________________
void AliESDtrack::GetRICHpid(Double_t *p) const {
  // Gets probabilities of each particle type (in RICH)
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fRICHr[i];
}



//_______________________________________________________________________
void AliESDtrack::SetESDpid(const Double_t *p) {  
  // Sets the probability of each particle type for the ESD track
  for (Int_t i=0; i<kSPECIES; i++) fR[i]=p[i];
  SetStatus(AliESDtrack::kESDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetESDpid(Double_t *p) const {
  // Gets probability of each particle type for the ESD track
  for (Int_t i=0; i<kSPECIES; i++) p[i]=fR[i];
}

//_______________________________________________________________________
void AliESDtrack::Print(Option_t *) const {
  // Prints info on the track
  
  printf("ESD track info\n") ; 
  Double_t p[kSPECIESN] ; 
  Int_t index = 0 ; 
  if( IsOn(kITSpid) ){
    printf("From ITS: ") ; 
    GetITSpid(p) ; 
    for(index = 0 ; index < kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetITSsignal()) ;
  } 
  if( IsOn(kTPCpid) ){
    printf("From TPC: ") ; 
    GetTPCpid(p) ; 
    for(index = 0 ; index < kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTPCsignal()) ;
  }
  if( IsOn(kTRDpid) ){
    printf("From TRD: ") ; 
    GetTRDpid(p) ; 
    for(index = 0 ; index < kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTRDsignal()) ;
  }
  if( IsOn(kTOFpid) ){
    printf("From TOF: ") ; 
    GetTOFpid(p) ; 
    for(index = 0 ; index < kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetTOFsignal()) ;
  }
  if( IsOn(kRICHpid) ){
    printf("From TOF: ") ; 
    GetRICHpid(p) ; 
    for(index = 0 ; index < kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetRICHsignal()) ;
  }
  if( IsOn(kPHOSpid) ){
    printf("From PHOS: ") ; 
    GetPHOSpid(p) ; 
    for(index = 0 ; index < kSPECIESN; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetPHOSsignal()) ;
  }
  if( IsOn(kEMCALpid) ){
    printf("From EMCAL: ") ; 
    GetEMCALpid(p) ; 
    for(index = 0 ; index < kSPECIESN; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetEMCALsignal()) ;
  }
} 
