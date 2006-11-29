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

#include <TMath.h>
#include <TParticle.h>

#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliKalmanTrack.h"
#include "AliLog.h"
#include "AliTrackPointArray.h"

ClassImp(AliESDtrack)

void SetPIDValues(Float_t * dest, const Double_t * src, Int_t n) {
  // This function copies "n" PID weights from "scr" to "dest"
  // and normalizes their sum to 1 thus producing conditional probabilities.
  // The negative weights are set to 0.
  // In case all the weights are non-positive they are replaced by
  // uniform probabilities

  if (n<=0) return;

  Float_t uniform = 1./(Float_t)n;

  Float_t sum = 0;
  for (Int_t i=0; i<n; i++) 
    if (src[i]>=0) {
      sum+=src[i];
      dest[i] = src[i];
    }
    else {
      dest[i] = 0;
    }

  if(sum>0)
    for (Int_t i=0; i<n; i++) dest[i] /= sum;
  else
    for (Int_t i=0; i<n; i++) dest[i] = uniform;
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack() : 
  AliExternalTrackParam(),
  fFlags(0),
  fLabel(0),
  fID(0),
  fTrackLength(0),
  fD(0),fZ(0),
  fCdd(0),fCdz(0),fCzz(0),
  fStopVertex(0),
  fCp(0),
  fCchi2(1e10),
  fIp(0),
  fOp(0),
  fITSchi2(0),
  fITSncls(0),
  fITSClusterMap(0),
  fITSsignal(0),
  fITSLabel(0),
  fTPCchi2(0),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCClusterMap(159),//number of padrows
  fTPCsignal(0),
  fTPCsignalN(0),
  fTPCsignalS(0),
  fTPCLabel(0),
  fTRDchi2(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDsignal(0),
  fTRDLabel(0),
  fTRDQuality(0),
  fTRDBudget(0),
  fTOFchi2(0),
  fTOFindex(0),
  fTOFCalChannel(-1),
  fTOFsignal(-1),
  fTOFsignalToT(0),
  fHMPIDchi2(1e10),
  fHMPIDqn(-1),
  fHMPIDcluIdx(-1),
  fHMPIDsignal(-1),
  fHMPIDtrkTheta(-1),
  fHMPIDtrkPhi(-1),
  fHMPIDtrkX(-1),
  fHMPIDtrkY(-1),
  fHMPIDmipX(-1),
  fHMPIDmipY(-1),
  fEMCALindex(kEMCALNoMatch),
  fFriendTrack(new AliESDfriendTrack())
{
  //
  // The default ESD constructor 
  //
  Int_t i, j;
  for (i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=1.;
    fITSr[i]=1.;
    fTPCr[i]=1.;
    fTRDr[i]=1.;
    fTOFr[i]=1.;
    fHMPIDr[i]=1.;
  }
  
  for (i=0; i<3; i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3; i++)   { fV0Indexes[i]=-1;}
  for (i=0;i<kNPlane;i++) {
    for (j=0;j<kNSlice;j++) {
      fTRDsignals[i][j]=0.; 
    }
    fTRDTimBin[i]=-1;
  }
  for (i=0;i<4;i++) {fTPCPoints[i]=-1;}
  for (i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (i=0;i<10;i++) {fTOFInfo[i]=-1;}
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack(const AliESDtrack& track):
  AliExternalTrackParam(track),
  fFlags(track.fFlags),
  fLabel(track.fLabel),
  fID(track.fID),
  fTrackLength(track.fTrackLength),
  fD(track.fD),fZ(track.fZ),
  fCdd(track.fCdd),fCdz(track.fCdz),fCzz(track.fCzz),
  fStopVertex(track.fStopVertex),
  fCp(0),
  fCchi2(track.fCchi2),
  fIp(0),
  fOp(0),
  fITSchi2(track.fITSchi2),
  fITSncls(track.fITSncls),
  fITSClusterMap(track.fITSClusterMap),
  fITSsignal(track.fITSsignal),
  fITSLabel(track.fITSLabel),
  fTPCchi2(track.fTPCchi2),
  fTPCncls(track.fTPCncls),
  fTPCnclsF(track.fTPCnclsF),
  fTPCClusterMap(track.fTPCClusterMap),
  fTPCsignal(track.fTPCsignal),
  fTPCsignalN(track.fTPCsignalN),
  fTPCsignalS(track.fTPCsignalS),
  fTPCLabel(track.fTPCLabel),
  fTRDchi2(track.fTRDchi2),
  fTRDncls(track.fTRDncls),
  fTRDncls0(track.fTRDncls0),
  fTRDsignal(track.fTRDsignal),
  fTRDLabel(track.fTRDLabel),
  fTRDQuality(track.fTRDQuality),
  fTRDBudget(track.fTRDBudget),
  fTOFchi2(track.fTOFchi2),
  fTOFindex(track.fTOFindex),
  fTOFCalChannel(track.fTOFCalChannel),
  fTOFsignal(track.fTOFsignal),
  fTOFsignalToT(track.fTOFsignalToT),
  fHMPIDchi2(track.fHMPIDchi2),
  fHMPIDqn(track.fHMPIDqn),
  fHMPIDcluIdx(track.fHMPIDcluIdx),
  fHMPIDsignal(track.fHMPIDsignal),
  fHMPIDtrkTheta(track.fHMPIDtrkTheta),
  fHMPIDtrkPhi(track.fHMPIDtrkPhi),
  fHMPIDtrkX(track.fHMPIDtrkX),
  fHMPIDtrkY(track.fHMPIDtrkY),
  fHMPIDmipX(track.fHMPIDmipX),
  fHMPIDmipY(track.fHMPIDmipY),
  fEMCALindex(track.fEMCALindex),
  fFriendTrack(0)
{
  //
  //copy constructor
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i]=track.fTrackTime[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++)  fR[i]=track.fR[i];
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]=track.fITSr[i]; 
  //
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=track.fTPCr[i]; 
  for (Int_t i=0;i<4;i++) {fTPCPoints[i]=track.fTPCPoints[i];}
  for (Int_t i=0; i<3;i++)   { fKinkIndexes[i]=track.fKinkIndexes[i];}
  for (Int_t i=0; i<3;i++)   { fV0Indexes[i]=track.fV0Indexes[i];}
  //
  for (Int_t i=0;i<kNPlane;i++) {
    for (Int_t j=0;j<kNSlice;j++) {
      fTRDsignals[i][j]=track.fTRDsignals[i][j]; 
    }
    fTRDTimBin[i]=track.fTRDTimBin[i];
  }
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i]=track.fTRDr[i]; 
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i]=track.fTOFr[i];
  for (Int_t i=0;i<3;i++) fTOFLabel[i]=track.fTOFLabel[i];
  for (Int_t i=0;i<10;i++) fTOFInfo[i]=track.fTOFInfo[i];
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fHMPIDr[i]=track.fHMPIDr[i];

  if (track.fCp) fCp=new AliExternalTrackParam(*track.fCp);
  if (track.fIp) fIp=new AliExternalTrackParam(*track.fIp);
  if (track.fOp) fOp=new AliExternalTrackParam(*track.fOp);

  if (track.fFriendTrack) fFriendTrack=new AliESDfriendTrack(*(track.fFriendTrack));
}

//_______________________________________________________________________
AliESDtrack::AliESDtrack(TParticle * part) : 
  AliExternalTrackParam(),
  fFlags(0),
  fLabel(0),
  fID(0),
  fTrackLength(0),
  fD(0),fZ(0),
  fCdd(0),fCdz(0),fCzz(0),
  fStopVertex(0),
  fCp(0),
  fCchi2(1e10),
  fIp(0),
  fOp(0),
  fITSchi2(0),
  fITSncls(0),
  fITSClusterMap(0),
  fITSsignal(0),
  fITSLabel(0),
  fTPCchi2(0),
  fTPCncls(0),
  fTPCnclsF(0),
  fTPCClusterMap(159),//number of padrows
  fTPCsignal(0),
  fTPCsignalN(0),
  fTPCsignalS(0),
  fTPCLabel(0),
  fTRDchi2(0),
  fTRDncls(0),
  fTRDncls0(0),
  fTRDsignal(0),
  fTRDLabel(0),
  fTRDQuality(0),
  fTRDBudget(0),
  fTOFchi2(0),
  fTOFindex(0),
  fTOFCalChannel(-1),
  fTOFsignal(-1),
  fTOFsignalToT(0),
  fHMPIDchi2(1e10),
  fHMPIDqn(-1),
  fHMPIDcluIdx(-1),
  fHMPIDsignal(-1),
  fHMPIDtrkTheta(-1),
  fHMPIDtrkPhi(-1),
  fHMPIDtrkX(-1),
  fHMPIDtrkY(-1),
  fHMPIDmipX(-1),
  fHMPIDmipY(-1),
  fEMCALindex(kEMCALNoMatch),
  fFriendTrack(0)
{
  //
  // ESD track from TParticle
  //

  // Reset all the arrays
  Int_t i, j;
  for (i=0; i<AliPID::kSPECIES; i++) {
    fTrackTime[i]=0.;
    fR[i]=0.;
    fITSr[i]=0.;
    fTPCr[i]=0.;
    fTRDr[i]=0.;
    fTOFr[i]=0.;
    fHMPIDr[i]=0.;
  }
  
  for (i=0; i<3; i++)   { fKinkIndexes[i]=0;}
  for (i=0; i<3; i++)   { fV0Indexes[i]=-1;}
  for (i=0;i<kNPlane;i++) {
    for (j=0;j<kNSlice;j++) {
      fTRDsignals[i][j]=0.; 
    }
    fTRDTimBin[i]=-1;
  }
  for (i=0;i<4;i++) {fTPCPoints[i]=-1;}
  for (i=0;i<3;i++) {fTOFLabel[i]=-1;}
  for (i=0;i<10;i++) {fTOFInfo[i]=-1;}

  // Calculate the AliExternalTrackParam content

  Double_t xref;
  Double_t alpha;
  Double_t param[5];
  Double_t covar[15];

  // Calculate alpha: the rotation angle of the corresponding local system (TPC sector)
  alpha = part->Phi()*180./TMath::Pi();
  if (alpha<0) alpha+= 360.;
  if (alpha>360) alpha -= 360.;

  Int_t sector = (Int_t)(alpha/20.);
  alpha = 10. + 20.*sector;
  alpha /= 180;
  alpha *= TMath::Pi();

  // Covariance matrix: no errors, the parameters are exact
  for (Int_t i=0; i<15; i++) covar[i]=0.;

  // Get the vertex of origin and the momentum
  TVector3 ver(part->Vx(),part->Vy(),part->Vz());
  TVector3 mom(part->Px(),part->Py(),part->Pz());

  // Rotate to the local coordinate system (TPC sector)
  ver.RotateZ(-alpha);
  mom.RotateZ(-alpha);

  // X of the referense plane
  xref = ver.X();

  Int_t pdgCode = part->GetPdgCode();

  Double_t charge = 
    TDatabasePDG::Instance()->GetParticle(pdgCode)->Charge();

  param[0] = ver.Y();
  param[1] = ver.Z();
  param[2] = TMath::Sin(mom.Phi());
  param[3] = mom.Pz()/mom.Pt();
  param[4] = TMath::Sign(1/mom.Pt(),charge);

  // Set AliExternalTrackParam
  Set(xref, alpha, param, covar);

  // Set the PID
  Int_t indexPID = 99;

  switch (TMath::Abs(pdgCode)) {

  case  11: // electron
    indexPID = 0;
    break;

  case 13: // muon
    indexPID = 1;
    break;

  case 211: // pion
    indexPID = 2;
    break;

  case 321: // kaon
    indexPID = 3;
    break;

  case 2212: // proton
    indexPID = 4;
    break;

  default:
    break;
  }

  // If the particle is not e,mu,pi,K or p the PID probabilities are set to 0
  if (indexPID < AliPID::kSPECIES) {
    fR[indexPID]=1.;
    fITSr[indexPID]=1.;
    fTPCr[indexPID]=1.;
    fTRDr[indexPID]=1.;
    fTOFr[indexPID]=1.;
    fHMPIDr[indexPID]=1.;

  }
  // AliESD track label
  SetLabel(part->GetUniqueID());

}

//_______________________________________________________________________
AliESDtrack::~AliESDtrack(){ 
  //
  // This is destructor according Coding Conventrions 
  //
  //printf("Delete track\n");
  delete fIp; 
  delete fOp;
  delete fCp; 
  delete fFriendTrack;
}

void AliESDtrack::AddCalibObject(TObject * object){
  //
  // add calib object to the list
  //
  if (!fFriendTrack) fFriendTrack  = new AliESDfriendTrack;
  fFriendTrack->AddCalibObject(object);
}

TObject *  AliESDtrack::GetCalibObject(Int_t index){
  //
  // return calib objct at given position
  //
  if (!fFriendTrack) return 0;
  return fFriendTrack->GetCalibObject(index);
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
  // Running track parameters in the base class (AliExternalTrackParam)
  
  fTrackLength = 0;
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTrackTime[i] = 0;
  fStopVertex = 0;

  // Reset track parameters constrained to the primary vertex
  fCp = 0;
  fCchi2 = 0;

  // Reset track parameters at the inner wall of TPC
  fIp = 0;

  // Reset track parameters at the inner wall of the TRD
  fOp = 0;

  // Reset ITS track related information
  fITSchi2 = 0;
  fITSncls = 0;       
  fITSClusterMap=0;
  fITSsignal = 0;     
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fITSr[i]=0; 
  fITSLabel = 0;       

  // Reset TPC related track information
  fTPCchi2 = 0;       
  fTPCncls = 0;       
  fTPCnclsF = 0;       
  fTPCClusterMap = 0;  
  fTPCsignal= 0;      
  fTPCsignalS= 0;      
  fTPCsignalN= 0;      
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTPCr[i]=0; 
  fTPCLabel=0;       
  for (Int_t i=0;i<4;i++) fTPCPoints[i] = 0;
  for (Int_t i=0; i<3;i++)   fKinkIndexes[i] = 0;
  for (Int_t i=0; i<3;i++)   fV0Indexes[i] = 0;

  // Reset TRD related track information
  fTRDchi2 = 0;        
  fTRDncls = 0;       
  fTRDncls0 = 0;       
  fTRDsignal = 0;      
  for (Int_t i=0;i<kNPlane;i++) {
    for (Int_t j=0;j<kNSlice;j++) {
      fTRDsignals[i][j] = 0; 
    }
    fTRDTimBin[i]  = 0;
  }
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTRDr[i] = 0; 
  fTRDLabel = 0;       
  fTRDQuality  = 0;
  fTRDBudget  = 0;

  // Reset TOF related track information
  fTOFchi2 = 0;        
  fTOFindex = 0;       
  fTOFsignal = 0;      
  fTOFCalChannel = -1;
  fTOFsignalToT = 0;
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fTOFr[i] = 0;
  for (Int_t i=0;i<3;i++) fTOFLabel[i] = 0;
  for (Int_t i=0;i<10;i++) fTOFInfo[i] = 0;

  // Reset HMPID related track information
  fHMPIDchi2 = 0;     
  fHMPIDqn = -1;     
  fHMPIDcluIdx = -1;     
  fHMPIDsignal = 0;     
  for (Int_t i=0;i<AliPID::kSPECIES;i++) fHMPIDr[i] = 0;
  fHMPIDtrkTheta = -1;     
  fHMPIDtrkPhi = -1;      
  fHMPIDtrkX = -1;     
  fHMPIDtrkY = -1;      
  fHMPIDmipX = -1;
  fHMPIDmipY = -1;
  fEMCALindex = kEMCALNoMatch;

  delete fFriendTrack; fFriendTrack = 0;
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
Bool_t AliESDtrack::UpdateTrackParams(const AliKalmanTrack *t, ULong_t flags){
  //
  // This function updates track's running parameters 
  //
  Int_t *index=0;
  Bool_t rc=kTRUE;

  SetStatus(flags);
  fLabel=t->GetLabel();

  if (t->IsStartedTimeIntegral()) {
    SetStatus(kTIME);
    Double_t times[10];t->GetIntegratedTimes(times); SetIntegratedTimes(times);
    SetIntegratedLength(t->GetIntegratedLength());
  }

  Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  
  switch (flags) {
    
  case kITSin: case kITSout: case kITSrefit:
    fITSncls=t->GetNumberOfClusters();
    index=fFriendTrack->GetITSindices(); 
    for (Int_t i=0;i<AliESDfriendTrack::kMaxITScluster;i++) {
        index[i]=t->GetClusterIndex(i);
	if (i<fITSncls) {
           Int_t l=(index[i] & 0xf0000000) >> 28;
           SETBIT(fITSClusterMap,l);                 
        }
    }
    fITSchi2=t->GetChi2();
    fITSsignal=t->GetPIDsignal();
    fITSLabel = t->GetLabel();
    break;
    
  case kTPCin: case kTPCrefit:
    fTPCLabel = t->GetLabel();
    if (!fIp) fIp=new AliExternalTrackParam(*t);
    else 
      fIp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
  case kTPCout:
    index=fFriendTrack->GetTPCindices(); 
    if (flags & kTPCout){
      if (!fOp) fOp=new AliExternalTrackParam(*t);
      else 
        fOp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
    }
    fTPCncls=t->GetNumberOfClusters();    
    fTPCchi2=t->GetChi2();
    
     {//prevrow must be declared in separate namespace, otherwise compiler cries:
      //"jump to case label crosses initialization of `Int_t prevrow'"
       Int_t prevrow = -1;
       //       for (Int_t i=0;i<fTPCncls;i++) 
       for (Int_t i=0;i<AliESDfriendTrack::kMaxTPCcluster;i++) 
        {
          index[i]=t->GetClusterIndex(i);
          Int_t idx = index[i];

	  if (idx<0) continue; 

          // Piotr's Cluster Map for HBT  
          // ### please change accordingly if cluster array is changing 
          // to "New TPC Tracking" style (with gaps in array) 
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
    break;

  case kTRDout: case kTRDin: case kTRDrefit:
    index=fFriendTrack->GetTRDindices();
    fTRDLabel = t->GetLabel(); 
    fTRDncls=t->GetNumberOfClusters();
    fTRDchi2=t->GetChi2();
    for (Int_t i=0;i<fTRDncls;i++) index[i]=t->GetClusterIndex(i);
    fTRDsignal=t->GetPIDsignal();
    break;
  case kTRDbackup:
    if (!fOp) fOp=new AliExternalTrackParam(*t);
    else 
      fOp->Set(t->GetX(),t->GetAlpha(),t->GetParameter(),t->GetCovariance());
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
void AliESDtrack::GetExternalParameters(Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters
  //---------------------------------------------------------------------
  x=GetX();
  for (Int_t i=0; i<5; i++) p[i]=GetParameter()[i];
}

//_______________________________________________________________________
void AliESDtrack::GetExternalCovariance(Double_t cov[15]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the cov. matrix
  //---------------------------------------------------------------------
  for (Int_t i=0; i<15; i++) cov[i]=AliExternalTrackParam::GetCovariance()[i];
}

//_______________________________________________________________________
Bool_t AliESDtrack::GetConstrainedExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns the constrained external track parameters
  //---------------------------------------------------------------------
  if (!fCp) return kFALSE;
  alpha=fCp->GetAlpha();
  x=fCp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fCp->GetParameter()[i];
  return kTRUE;
}

//_______________________________________________________________________
Bool_t 
AliESDtrack::GetConstrainedExternalCovariance(Double_t c[15]) const {
  //---------------------------------------------------------------------
  // This function returns the constrained external cov. matrix
  //---------------------------------------------------------------------
  if (!fCp) return kFALSE;
  for (Int_t i=0; i<15; i++) c[i]=fCp->GetCovariance()[i];
  return kTRUE;
}

Bool_t
AliESDtrack::GetInnerExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters 
  // at the inner layer of TPC
  //---------------------------------------------------------------------
  if (!fIp) return kFALSE;
  alpha=fIp->GetAlpha();
  x=fIp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fIp->GetParameter()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetInnerExternalCovariance(Double_t cov[15]) const {
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix 
 // at the inner layer of TPC
 //---------------------------------------------------------------------
  if (!fIp) return kFALSE;
  for (Int_t i=0; i<15; i++) cov[i]=fIp->GetCovariance()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetOuterExternalParameters
                 (Double_t &alpha, Double_t &x, Double_t p[5]) const {
  //---------------------------------------------------------------------
  // This function returns external representation of the track parameters 
  // at the inner layer of TRD
  //---------------------------------------------------------------------
  if (!fOp) return kFALSE;
  alpha=fOp->GetAlpha();
  x=fOp->GetX();
  for (Int_t i=0; i<5; i++) p[i]=fOp->GetParameter()[i];
  return kTRUE;
}

Bool_t 
AliESDtrack::GetOuterExternalCovariance(Double_t cov[15]) const {
 //---------------------------------------------------------------------
 // This function returns external representation of the cov. matrix 
 // at the inner layer of TRD
 //---------------------------------------------------------------------
  if (!fOp) return kFALSE;
  for (Int_t i=0; i<15; i++) cov[i]=fOp->GetCovariance()[i];
  return kTRUE;
}

Int_t AliESDtrack::GetNcls(Int_t idet) const
{
  // Get number of clusters by subdetector index
  //
  Int_t ncls = 0;
  switch(idet){
  case 0:
    ncls = fITSncls;
    break;
  case 1:
    ncls = fTPCncls;
    break;
  case 2:
    ncls = fTRDncls;
    break;
  case 3:
    if (fTOFindex != 0)
      ncls = 1;
    break;
  default:
    break;
  }
  return ncls;
}

Int_t AliESDtrack::GetClusters(Int_t idet, Int_t *idx) const
{
  // Get cluster index array by subdetector index
  //
  Int_t ncls = 0;
  switch(idet){
  case 0:
    ncls = GetITSclusters(idx);
    break;
  case 1:
    ncls = GetTPCclusters(idx);
    break;
  case 2:
    ncls = GetTRDclusters(idx);
    break;
  case 3:
    if (fTOFindex != 0) {
      idx[0] = GetTOFcluster();
      ncls = 1;
    }
    break;
  default:
    break;
  }
  return ncls;
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
  SetPIDValues(fITSr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kITSpid);
}

//_______________________________________________________________________
void AliESDtrack::GetITSpid(Double_t *p) const {
  // Gets the probability of each particle type (in ITS)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fITSr[i];
}

//_______________________________________________________________________
Int_t AliESDtrack::GetITSclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  if (idx!=0) {
     Int_t *index=fFriendTrack->GetITSindices();
     for (Int_t i=0; i<AliESDfriendTrack::kMaxITScluster; i++) idx[i]=index[i];
  }
  return fITSncls;
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTPCclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined ITS clusters 
  //---------------------------------------------------------------------
  if (idx!=0) {
    Int_t *index=fFriendTrack->GetTPCindices();
    for (Int_t i=0; i<AliESDfriendTrack::kMaxTPCcluster; i++) idx[i]=index[i];
  }
  return fTPCncls;
}

Float_t AliESDtrack::GetTPCdensity(Int_t row0, Int_t row1) const{
  //
  // GetDensity of the clusters on given region between row0 and row1
  // Dead zone effect takin into acoount
  //
  Int_t good  = 0;
  Int_t found = 0;
  //  
  Int_t *index=fFriendTrack->GetTPCindices();
  for (Int_t i=row0;i<=row1;i++){     
    Int_t idx = index[i];
    if (idx!=-1)  good++;             // track outside of dead zone
    if (idx>0)    found++;
  }
  Float_t density=0.5;
  if (good>(row1-row0)*0.5) density = Float_t(found)/Float_t(good);
  return density;
}

//_______________________________________________________________________
void AliESDtrack::SetTPCpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TPC)
  SetPIDValues(fTPCr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kTPCpid);
}

//_______________________________________________________________________
void AliESDtrack::GetTPCpid(Double_t *p) const {
  // Gets the probability of each particle type (in TPC)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fTPCr[i];
}

//_______________________________________________________________________
Int_t AliESDtrack::GetTRDclusters(Int_t *idx) const {
  //---------------------------------------------------------------------
  // This function returns indices of the assgined TRD clusters 
  //---------------------------------------------------------------------
  if (idx!=0) {
     Int_t *index=fFriendTrack->GetTRDindices();
     for (Int_t i=0; i<AliESDfriendTrack::kMaxTRDcluster; i++) idx[i]=index[i];
  }
  return fTRDncls;
}

//_______________________________________________________________________
void AliESDtrack::SetTRDpid(const Double_t *p) {  
  // Sets values for the probability of each particle type (in TRD)
  SetPIDValues(fTRDr,p,AliPID::kSPECIES);
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
  SetPIDValues(fTOFr,p,AliPID::kSPECIES);
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
void AliESDtrack::SetHMPIDpid(const Double_t *p) {  
  // Sets the probability of each particle type (in HMPID)
  SetPIDValues(fHMPIDr,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kHMPIDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetHMPIDpid(Double_t *p) const {
  // Gets probabilities of each particle type (in HMPID)
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fHMPIDr[i];
}



//_______________________________________________________________________
void AliESDtrack::SetESDpid(const Double_t *p) {  
  // Sets the probability of each particle type for the ESD track
  SetPIDValues(fR,p,AliPID::kSPECIES);
  SetStatus(AliESDtrack::kESDpid);
}

//_______________________________________________________________________
void AliESDtrack::GetESDpid(Double_t *p) const {
  // Gets probability of each particle type for the ESD track
  for (Int_t i=0; i<AliPID::kSPECIES; i++) p[i]=fR[i];
}

//_______________________________________________________________________
Bool_t AliESDtrack::RelateToVertex
(const AliESDVertex *vtx, Double_t b, Double_t maxd) {
  //
  // Try to relate this track to the vertex "vtx", 
  // if the (rough) transverse impact parameter is not bigger then "maxd". 
  //            Magnetic field is "b" (kG).
  //
  // a) The track gets extapolated to the DCA to the vertex.
  // b) The impact parameters and their covariance matrix are calculated.
  // c) An attempt to constrain this track to the vertex is done.
  //
  //    In the case of success, the returned value is kTRUE
  //    (otherwise, it's kFALSE)
  //  

  if (!vtx) return kFALSE;

  Double_t alpha=GetAlpha();
  Double_t sn=TMath::Sin(alpha), cs=TMath::Cos(alpha);
  Double_t x=GetX(), y=GetParameter()[0], snp=GetParameter()[2];
  Double_t xv= vtx->GetXv()*cs + vtx->GetYv()*sn;
  Double_t yv=-vtx->GetXv()*sn + vtx->GetYv()*cs, zv=vtx->GetZv();
  x-=xv; y-=yv;

  //Estimate the impact parameter neglecting the track curvature
  Double_t d=TMath::Abs(x*snp - y*TMath::Sqrt(1.- snp*snp));
  if (d > maxd) return kFALSE; 

  //Propagate to the DCA
  Double_t crv=kB2C*b*GetParameter()[4];
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;

  Double_t tgfv=-(crv*x - snp)/(crv*y + TMath::Sqrt(1.-snp*snp));
  sn=tgfv/TMath::Sqrt(1.+ tgfv*tgfv);
  if (TMath::Abs(tgfv)>0.) cs = sn/tgfv;
  else cs=1.;

  x = xv*cs + yv*sn;
  yv=-xv*sn + yv*cs; xv=x;

  if (!Propagate(alpha+TMath::ASin(sn),xv,b)) return kFALSE;

  fD = GetParameter()[0] - yv;
  fZ = GetParameter()[1] - zv;
  
  Double_t cov[6]; vtx->GetCovMatrix(cov);

  //***** Improvements by A.Dainese
  alpha=GetAlpha(); sn=TMath::Sin(alpha); cs=TMath::Cos(alpha);
  Double_t s2ylocvtx = cov[0]*sn*sn + cov[2]*cs*cs - 2.*cov[1]*cs*sn;
  fCdd = GetCovariance()[0] + s2ylocvtx;   // neglecting correlations
  fCdz = GetCovariance()[1];               // between (x,y) and z
  fCzz = GetCovariance()[2] + cov[5];      // in vertex's covariance matrix
  //*****

  {//Try to constrain 
    Double_t p[2]={yv,zv}, c[3]={cov[2],0.,cov[5]};
    Double_t chi2=GetPredictedChi2(p,c);

    if (chi2>77.) return kFALSE;

    AliExternalTrackParam tmp(*this);
    if (!tmp.Update(p,c)) return kFALSE;

    fCchi2=chi2;
    if (!fCp) fCp=new AliExternalTrackParam();
    new (fCp) AliExternalTrackParam(tmp);
  }

  return kTRUE;
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
  if( IsOn(kHMPIDpid) ){
    printf("From HMPID: ") ; 
    GetHMPIDpid(p) ; 
    for(index = 0 ; index < AliPID::kSPECIES; index++) 
      printf("%f, ", p[index]) ;
    printf("\n           signal = %f\n", GetHMPIDsignal()) ;
  }
} 

Bool_t AliESDtrack::PropagateTo(Double_t xToGo, Double_t b, Double_t mass, 
Double_t maxStep, Bool_t rotateTo, Double_t maxSnp){
  //----------------------------------------------------------------
  //
  // MI's function
  //
  // Propagates this track to the plane X=xk (cm) 
  // in the magnetic field "b" (kG),
  // the correction for the material is included
  //
  // mass     - mass used in propagation - used for energy loss correction
  // maxStep  - maximal step for propagation
  //----------------------------------------------------------------
  const Double_t kEpsilon = 0.00001;
  Double_t xpos     = GetX();
  Double_t dir      = (xpos<xToGo) ? 1.:-1.;
  //
  while ( (xToGo-xpos)*dir > kEpsilon){
    Double_t step = dir*TMath::Min(TMath::Abs(xToGo-xpos), maxStep);
    Double_t x    = xpos+step;
    Double_t xyz0[3],xyz1[3],param[7];
    GetXYZ(xyz0);   //starting global position
    if (!GetXYZAt(x,b,xyz1)) return kFALSE;   // no prolongation
    xyz1[2]+=kEpsilon; // waiting for bug correction in geo
    AliKalmanTrack::MeanMaterialBudget(xyz0,xyz1,param);	
    if (TMath::Abs(GetSnpAt(x,b)) >= maxSnp) return kFALSE;
    if (!AliExternalTrackParam::PropagateTo(x,b))  return kFALSE;

    Double_t rho=param[0],x0=param[1],distance=param[4];
    Double_t d=distance*rho/x0;

    if (!CorrectForMaterial(d,x0,mass)) return kFALSE;
    if (rotateTo){
      if (TMath::Abs(GetSnp()) >= maxSnp) return kFALSE;
      GetXYZ(xyz0);   // global position
      Double_t alphan = TMath::ATan2(xyz0[1], xyz0[0]); 
      //
      Double_t ca=TMath::Cos(alphan-GetAlpha()), 
               sa=TMath::Sin(alphan-GetAlpha());
      Double_t sf=GetSnp(), cf=TMath::Sqrt(1.- sf*sf);
      Double_t sinNew =  sf*ca - cf*sa;
      if (TMath::Abs(sinNew) >= maxSnp) return kFALSE;
      if (!Rotate(alphan)) return kFALSE;
    }
    xpos = GetX();
  }
  return kTRUE;
}
