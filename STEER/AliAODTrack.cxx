/**************************************************************************
 * Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//-------------------------------------------------------------------------
//     AOD track implementation of AliVParticle
//     Author: Markus Oldenburg, CERN
//     Markus.Oldenburg@cern.ch
//-------------------------------------------------------------------------

#include "AliLog.h"
#include "AliAODTrack.h"

ClassImp(AliAODTrack)

//______________________________________________________________________________
AliAODTrack::AliAODTrack() : 
  AliVParticle(),
  fChi2perNDF(-999.),
  fChi2MatchTrigger(0.),
  fFlags(0),
  fLabel(-999),
  fITSMuonClusterMap(0),
  fFilterMap(0),
  fID(-999),
  fCharge(-99),
  fType(kUndef),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fProdVertex(NULL)
{
  // default constructor

  SetP();
  SetPosition((Float_t*)NULL);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  SetPID((Float_t*)NULL);
}

//______________________________________________________________________________
AliAODTrack::AliAODTrack(Short_t id,
			 Int_t label, 
			 Double_t p[3],
			 Bool_t cartesian,
			 Double_t x[3],
			 Bool_t isDCA,
			 Double_t covMatrix[21],
			 Short_t charge,
			 UChar_t itsClusMap,
			 Double_t pid[10],
			 AliAODVertex *prodVertex,
			 Bool_t usedForVtxFit,
			 Bool_t usedForPrimVtxFit,
			 AODTrk_t ttype,
			 UInt_t selectInfo,
			 Float_t chi2perNDF) :
  AliVParticle(),
  fChi2perNDF(chi2perNDF),
  fChi2MatchTrigger(0.),
  fFlags(0),
  fLabel(label),
  fITSMuonClusterMap(0),
  fFilterMap(selectInfo),
  fID(id),
  fCharge(charge),
  fType(ttype),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fProdVertex(prodVertex)
{
  // constructor
 
  SetP(p, cartesian);
  SetPosition(x, isDCA);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  SetUsedForVtxFit(usedForVtxFit);
  SetUsedForPrimVtxFit(usedForPrimVtxFit);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);
  SetITSClusterMap(itsClusMap);
}

//______________________________________________________________________________
AliAODTrack::AliAODTrack(Short_t id,
			 Int_t label, 
			 Float_t p[3],
			 Bool_t cartesian,
			 Float_t x[3],
			 Bool_t isDCA,
			 Float_t covMatrix[21],
			 Short_t charge,
			 UChar_t itsClusMap,
			 Float_t pid[10],
			 AliAODVertex *prodVertex,
			 Bool_t usedForVtxFit,
			 Bool_t usedForPrimVtxFit,
			 AODTrk_t ttype,
			 UInt_t selectInfo,
			 Float_t chi2perNDF) :
  AliVParticle(),
  fChi2perNDF(chi2perNDF),
  fChi2MatchTrigger(0.),
  fFlags(0),
  fLabel(label),
  fITSMuonClusterMap(0),
  fFilterMap(selectInfo),
  fID(id),
  fCharge(charge),
  fType(ttype),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fProdVertex(prodVertex)
{
  // constructor
 
  SetP(p, cartesian);
  SetPosition(x, isDCA);
  SetXYAtDCA(-999., -999.);
  SetPxPyPzAtDCA(-999., -999., -999.);
  SetUsedForVtxFit(usedForVtxFit);
  SetUsedForPrimVtxFit(usedForPrimVtxFit);
  if(covMatrix) SetCovMatrix(covMatrix);
  SetPID(pid);
  SetITSClusterMap(itsClusMap);
}

//______________________________________________________________________________
AliAODTrack::~AliAODTrack() 
{
  // destructor
  delete fCovMatrix;
  delete fDetPid;
}


//______________________________________________________________________________
AliAODTrack::AliAODTrack(const AliAODTrack& trk) :
  AliVParticle(trk),
  fChi2perNDF(trk.fChi2perNDF),
  fChi2MatchTrigger(trk.fChi2MatchTrigger),
  fFlags(trk.fFlags),
  fLabel(trk.fLabel),
  fITSMuonClusterMap(trk.fITSMuonClusterMap),
  fFilterMap(trk.fFilterMap),
  fID(trk.fID),
  fCharge(trk.fCharge),
  fType(trk.fType),
  fCovMatrix(NULL),
  fDetPid(NULL),
  fProdVertex(trk.fProdVertex)
{
  // Copy constructor

  trk.GetP(fMomentum);
  trk.GetPosition(fPosition);
  SetXYAtDCA(trk.XAtDCA(), trk.YAtDCA());
  SetPxPyPzAtDCA(trk.PxAtDCA(), trk.PyAtDCA(), trk.PzAtDCA());
  SetUsedForVtxFit(trk.GetUsedForVtxFit());
  SetUsedForPrimVtxFit(trk.GetUsedForPrimVtxFit());
  if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<6>(*trk.fCovMatrix);
  if(trk.fDetPid) fDetPid=new AliAODPid(*trk.fDetPid);
  SetPID(trk.fPID);
}

//______________________________________________________________________________
AliAODTrack& AliAODTrack::operator=(const AliAODTrack& trk)
{
  // Assignment operator
  if(this!=&trk) {

    AliVParticle::operator=(trk);

    trk.GetP(fMomentum);
    trk.GetPosition(fPosition);
    trk.GetPID(fPID);

    SetXYAtDCA(trk.XAtDCA(), trk.YAtDCA());
    SetPxPyPzAtDCA(trk.PxAtDCA(), trk.PyAtDCA(), trk.PzAtDCA());
    
    fChi2perNDF = trk.fChi2perNDF;
    fChi2MatchTrigger = trk.fChi2MatchTrigger;

    fFlags = trk.fFlags;
    fLabel = trk.fLabel;    
    
    fITSMuonClusterMap = trk.fITSMuonClusterMap;
    fFilterMap = trk.fFilterMap;

    fID = trk.fID;

    fCharge = trk.fCharge;
    fType = trk.fType;

    delete fCovMatrix;
    if(trk.fCovMatrix) fCovMatrix=new AliAODRedCov<6>(*trk.fCovMatrix);
    else fCovMatrix=NULL;
    fProdVertex = trk.fProdVertex;

    SetUsedForVtxFit(trk.GetUsedForVtxFit());
    SetUsedForPrimVtxFit(trk.GetUsedForPrimVtxFit());

    delete fDetPid;
    if(trk.fDetPid) fDetPid=new AliAODPid(*trk.fDetPid);
    else fDetPid=NULL;
  }

  return *this;
}

//______________________________________________________________________________
Double_t AliAODTrack::M(AODTrkPID_t pid) const
{
  // Returns the mass.
  // Masses for nuclei don't exist in the PDG tables, therefore they were put by hand.

  switch (pid) {

  case kElectron :
    return 0.000510999; //TDatabasePDG::Instance()->GetParticle(11/*::kElectron*/)->Mass();
    break;

  case kMuon :
    return 0.1056584; //TDatabasePDG::Instance()->GetParticle(13/*::kMuonMinus*/)->Mass();
    break;

  case kPion :
    return 0.13957; //TDatabasePDG::Instance()->GetParticle(211/*::kPiPlus*/)->Mass();
    break;

  case kKaon :
    return 0.4937; //TDatabasePDG::Instance()->GetParticle(321/*::kKPlus*/)->Mass();
    break;

  case kProton :
    return 0.9382720; //TDatabasePDG::Instance()->GetParticle(2212/*::kProton*/)->Mass();
    break;

  case kDeuteron :
    return 1.8756; //TDatabasePDG::Instance()->GetParticle(1000010020)->Mass();
    break;

  case kTriton :
    return 2.8089; //TDatabasePDG::Instance()->GetParticle(1000010030)->Mass();
    break;

  case kHelium3 :
    return 2.8084; //TDatabasePDG::Instance()->GetParticle(1000020030)->Mass();
    break;

  case kAlpha :
    return 3.7274; //TDatabasePDG::Instance()->GetParticle(1000020040)->Mass();
    break;

  case kUnknown :
    return -999.;
    break;

  default :
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::E(AODTrkPID_t pid) const
{
  // Returns the energy of the particle of a given pid.
  
  if (pid != kUnknown) { // particle was identified
    Double_t m = M(pid);
    return TMath::Sqrt(P()*P() + m*m);
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::Y(AODTrkPID_t pid) const
{
  // Returns the rapidity of a particle of a given pid.
  
  if (pid != kUnknown) { // particle was identified
    Double_t e = E(pid);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
Double_t AliAODTrack::Y(Double_t m) const
{
  // Returns the rapidity of a particle of a given mass.
  
  if (m >= 0.) { // mass makes sense
    Double_t e = E(m);
    Double_t pz = Pz();
    if (e>=0 && e!=pz) { // energy was positive (e.g. not -999.) and not equal to pz
      return 0.5*TMath::Log((e+pz)/(e-pz));
    } else { // energy not known or equal to pz
      return -999.;
    }
  } else { // pid unknown
    return -999.;
  }
}

//______________________________________________________________________________
AliAODTrack::AODTrkPID_t AliAODTrack::GetMostProbablePID() const 
{
  // Returns the most probable PID array element.
  
  Int_t nPID = 10;
  if (fPID) {
    AODTrkPID_t loc = kUnknown;
    Double_t max = 0.;
    Bool_t allTheSame = kTRUE;
    
    for (Int_t iPID = 0; iPID < nPID; iPID++) {
      if (fPID[iPID] >= max) {
	if (fPID[iPID] > max) {
	  allTheSame = kFALSE;
	  max = fPID[iPID];
	  loc = (AODTrkPID_t)iPID;
	} else {
	  allTheSame = kTRUE;
	}
      }
    }
    
    return allTheSame ? kUnknown : loc;
  } else {
    return kUnknown;
  }
}

//______________________________________________________________________________
void AliAODTrack::ConvertAliPIDtoAODPID()
{
  // Converts AliPID array.
  // The numbering scheme is the same for electrons, muons, pions, kaons, and protons.
  // Everything else has to be set to zero.

  fPID[kDeuteron] = 0.;
  fPID[kTriton]   = 0.;
  fPID[kHelium3]  = 0.;
  fPID[kAlpha]    = 0.;
  fPID[kUnknown]  = 0.;
  
  return;
}


//______________________________________________________________________________
template <class T> void AliAODTrack::SetP(const T *p, const Bool_t cartesian) 
{
  // Set the momentum

  if (p) {
    if (cartesian) {
      Double_t pt2 = p[0]*p[0] + p[1]*p[1];
      Double_t pp  = TMath::Sqrt(pt2 + p[2]*p[2]);
      
      fMomentum[0] = TMath::Sqrt(pt2); // pt
      fMomentum[1] = (pt2 != 0.) ? TMath::Pi()+TMath::ATan2(-p[1], -p[0]) : -999; // phi
      fMomentum[2] = (pp != 0.) ? TMath::ACos(p[2] / pp) : -999.; // theta
    } else {
      fMomentum[0] = p[0];  // pt
      fMomentum[1] = p[1];  // phi
      fMomentum[2] = p[2];  // theta
    }
  } else {
    fMomentum[0] = -999.;
    fMomentum[1] = -999.;
    fMomentum[2] = -999.;
  }
}

//______________________________________________________________________________
template <class T> void AliAODTrack::SetPosition(const T *x, const Bool_t dca) 
{
  // set the position

  if (x) {
    if (!dca) {
      ResetBit(kIsDCA);

      fPosition[0] = x[0];
      fPosition[1] = x[1];
      fPosition[2] = x[2];
    } else {
      SetBit(kIsDCA);
      // don't know any better yet
      fPosition[0] = -999.;
      fPosition[1] = -999.;
      fPosition[2] = -999.;
    }
  } else {
    ResetBit(kIsDCA);

    fPosition[0] = -999.;
    fPosition[1] = -999.;
    fPosition[2] = -999.;
  }
}

//______________________________________________________________________________
void AliAODTrack::SetDCA(Double_t d, Double_t z) 
{
  // set the dca
  fPosition[0] = d;
  fPosition[1] = z;
  fPosition[2] = 0.;
  SetBit(kIsDCA);
}

//______________________________________________________________________________
void AliAODTrack::Print(Option_t* /* option */) const
{
  // prints information about AliAODTrack

  printf("Object name: %s   Track type: %s\n", GetName(), GetTitle()); 
  printf("        px = %f\n", Px());
  printf("        py = %f\n", Py());
  printf("        pz = %f\n", Pz());
  printf("        pt = %f\n", Pt());
  printf("      1/pt = %f\n", OneOverPt());
  printf("     theta = %f\n", Theta());
  printf("       phi = %f\n", Phi());
  printf("  chi2/NDF = %f\n", Chi2perNDF());
  printf("    charge = %d\n", Charge());
}

void AliAODTrack::SetMatchTrigger(Int_t matchTrig){
//
// Set the MUON trigger information
  switch(matchTrig){
    case 0: // 0 track does not match trigger
      fITSMuonClusterMap=fITSMuonClusterMap&0x3fffffff;
      break;
    case 1: // 1 track match but does not pass pt cut
      fITSMuonClusterMap=(fITSMuonClusterMap&0x3fffffff)|0x40000000;
      break;
    case 2: // 2 track match Low pt cut
      fITSMuonClusterMap=(fITSMuonClusterMap&0x3fffffff)|0x80000000;
      break;
    case 3: // 3 track match High pt cut
      fITSMuonClusterMap=fITSMuonClusterMap|0xc0000000;
      break;
    default:
      fITSMuonClusterMap=fITSMuonClusterMap&0x3fffffff;
      AliWarning(Form("unknown case for matchTrig: %d\n",matchTrig));
  }
}

Int_t AliAODTrack::HitsMT(Int_t istation, Int_t iplane, Char_t *cathode){
//
// Retrieve hit information for MUON identified by  (station, plane, cathode)
  if(cathode){
    if(cathode[0]=='x'||cathode[0]=='X'){
      if(istation==1){
        if(iplane==1)
	  return (fITSMuonClusterMap&0x8000)?1:0;
	else if(iplane==2)
	  return (fITSMuonClusterMap&0x4000)?1:0;
	else
	  return 0;
      }else if(istation==2){
        if(iplane==1)
	  return (fITSMuonClusterMap&0x2000)?1:0;
	else if(iplane==2)
	  return (fITSMuonClusterMap&0x1000)?1:0;
	else
	  return 0;
      }else{
        return 0;
      }
    }else if(cathode[0]=='y'||cathode[0]=='Y'){
      if(istation==1){
        if(iplane==1)
	  return (fITSMuonClusterMap&0x0800)?1:0;
	else if(iplane==2)
	  return (fITSMuonClusterMap&0x0400)?1:0;
	else
	  return 0;
      }else if(istation==2){
        if(iplane==1)
	  return (fITSMuonClusterMap&0x0200)?1:0;
	else if(iplane==2)
	  return (fITSMuonClusterMap&0x0100)?1:0;
	else
	  return 0;
      }else{
        return 0;
      }
    }else{
      return 0;
    }
  }else{
    if(istation==1){
      if(iplane==1)
	return (HitsMT(1,1,"X")||HitsMT(1,1,"Y"))?1:0;
      else if(iplane==2)
	return (HitsMT(1,2,"X")||HitsMT(1,2,"Y"))?1:0;
      else
	return 0;
    }else if(istation==2){
      if(iplane==1)
	return (HitsMT(2,1,"X")||HitsMT(2,1,"Y"))?1:0;
      else if(iplane==2)
	return (HitsMT(2,2,"X")||HitsMT(2,2,"Y"))?1:0;
      else
	return 0;
    }else{
      return 0;
    }
  }
}

Int_t AliAODTrack::HitsMuonChamber(Int_t MuonChamber){
  //
  // Retrieve hit information for MUON Tracker/Trigger Chamber
  // WARNING: chamber number start from 1 instead of 0
  
  if (MuonChamber > 0 && MuonChamber < 11) {
    return ((GetMUONClusterMap() & BIT(MuonChamber-1)) != 0) ? 1 : 0;
  } else {
    switch(MuonChamber){
      case 11:
	return HitsMT(1,1);
      case 12:
	return HitsMT(1,2);
      case 13:
	return HitsMT(2,1);
      case 14:
	return HitsMT(2,2);
      default:
	printf("Unknown MUON chamber: %d\n",MuonChamber);
	return 0;
    }
  }
}

Bool_t AliAODTrack::PropagateTo(Double_t xk, Double_t b) {
  //----------------------------------------------------------------
  // Propagate this track to the plane X=xk (cm) in the field "b" (kG)
  // This is in local coordinates!!!
  //----------------------------------------------------------------

  Double_t alpha = 0.;
  Double_t localP[3] = {Px(), Py(), Pz()}; // set global (sic!) p
  Global2LocalMomentum(localP, Charge(), alpha); // convert global to local momentum

  AliAODVertex *origin = (AliAODVertex*)fProdVertex.GetObject();
  Double_t localX[3] = {origin->GetX(), origin->GetY(), origin->GetZ()}; // set global (sic!) location of first track point
  Global2LocalPosition(localX, alpha); // convert global to local position

  Double_t &fX = localX[0];

  Double_t dx=xk-fX;
  if (TMath::Abs(dx)<=kAlmost0)  return kTRUE;

  Double_t crv=localP[0]*b*kB2C;
  if (TMath::Abs(b) < kAlmost0Field) crv=0.;

  Double_t f1=localP[1], f2=f1 + crv*dx;
  if (TMath::Abs(f1) >= kAlmost1) return kFALSE;
  if (TMath::Abs(f2) >= kAlmost1) return kFALSE;

  Double_t &fP0=localX[1], &fP1=localX[2], &fP2=localP[0], &fP3=localP[1], &fP4=localP[2];
  /* covariance matrix to be fixed! 
  Double_t 
  &fC00=fC[0],
  &fC10=fC[1],   &fC11=fC[2],  
  &fC20=fC[3],   &fC21=fC[4],   &fC22=fC[5],
  &fC30=fC[6],   &fC31=fC[7],   &fC32=fC[8],   &fC33=fC[9],  
  &fC40=fC[10],  &fC41=fC[11],  &fC42=fC[12],  &fC43=fC[13], &fC44=fC[14];
  */
  Double_t r1=TMath::Sqrt(1.- f1*f1), r2=TMath::Sqrt(1.- f2*f2);

  fX=xk;
  fP0 += dx*(f1+f2)/(r1+r2);
  fP1 += dx*(r2 + f2*(f1+f2)/(r1+r2))*fP3;
  fP2 += dx*crv;

  //f = F - 1
   
  //Double_t f02=    dx/(r1*r1*r1);            
  Double_t cc=crv/fP4;
  Double_t f04=0.5*dx*dx/(r1*r1*r1);         f04*=cc;
  //Double_t f12=    dx*fP3*f1/(r1*r1*r1);
  Double_t f14=0.5*dx*dx*fP3*f1/(r1*r1*r1);  f14*=cc;
  //Double_t f13=    dx/r1;
  Double_t f24=    dx;                       f24*=cc;
  
  /* covariance matrix to be fixed!
  //b = C*ft
  Double_t b00=f02*fC20 + f04*fC40, b01=f12*fC20 + f14*fC40 + f13*fC30;
  Double_t b02=f24*fC40;
  Double_t b10=f02*fC21 + f04*fC41, b11=f12*fC21 + f14*fC41 + f13*fC31;
  Double_t b12=f24*fC41;
  Double_t b20=f02*fC22 + f04*fC42, b21=f12*fC22 + f14*fC42 + f13*fC32;
  Double_t b22=f24*fC42;
  Double_t b40=f02*fC42 + f04*fC44, b41=f12*fC42 + f14*fC44 + f13*fC43;
  Double_t b42=f24*fC44;
  Double_t b30=f02*fC32 + f04*fC43, b31=f12*fC32 + f14*fC43 + f13*fC33;
  Double_t b32=f24*fC43;
  
  //a = f*b = f*C*ft
  Double_t a00=f02*b20+f04*b40,a01=f02*b21+f04*b41,a02=f02*b22+f04*b42;
  Double_t a11=f12*b21+f14*b41+f13*b31,a12=f12*b22+f14*b42+f13*b32;
  Double_t a22=f24*b42;

  //F*C*Ft = C + (b + bt + a)
  fC00 += b00 + b00 + a00;
  fC10 += b10 + b01 + a01; 
  fC20 += b20 + b02 + a02;
  fC30 += b30;
  fC40 += b40;
  fC11 += b11 + b11 + a11;
  fC21 += b21 + b12 + a12;
  fC31 += b31; 
  fC41 += b41;
  fC22 += b22 + b22 + a22;
  fC32 += b32;
  fC42 += b42;
  */
  
  Local2GlobalMomentum(localP, alpha); // convert local to global momentum
  SetP(localP);

  return kTRUE;
}

