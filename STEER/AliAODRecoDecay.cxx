/**************************************************************************
 * Copyright(c) 1998-2006, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//
// Base class for AOD reconstructed decay
//
// Author: A.Dainese, andrea.dainese@lnl.infn.it
/////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TClonesArray.h>
#include "AliLog.h"
#include "AliVTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODRecoDecay.h"

ClassImp(AliAODRecoDecay)

//--------------------------------------------------------------------------
AliAODRecoDecay::AliAODRecoDecay() :
  AliVTrack(),
  fSecondaryVtx(0x0),
  fOwnSecondaryVtx(0x0),
  fCharge(0),
  fNProngs(0), fNDCA(0), fNPID(0),
  fPx(0x0), fPy(0x0), fPz(0x0),
  fd0(0x0),
  fDCA(0x0),
  fPID(0x0), 
  fEventNumber(-1),fRunNumber(-1)
{
  //
  // Default Constructor
  //
}
//--------------------------------------------------------------------------
AliAODRecoDecay::AliAODRecoDecay(AliAODVertex *vtx2,Int_t nprongs,
				 Short_t charge,
				 Double_t *px,Double_t *py,Double_t *pz,
				 Double_t *d0) :
  AliVTrack(),
  fSecondaryVtx(vtx2),
  fOwnSecondaryVtx(0x0),
  fCharge(charge),
  fNProngs(nprongs), fNDCA(0), fNPID(0),
  fPx(0x0), fPy(0x0), fPz(0x0),
  fd0(0x0),
  fDCA(0x0),
  fPID(0x0), 
  fEventNumber(-1),fRunNumber(-1)
{
  //
  // Constructor with AliAODVertex for decay vertex
  //

  fPx = new Double_t[GetNProngs()];
  fPy = new Double_t[GetNProngs()];
  fPz = new Double_t[GetNProngs()];
  fd0 = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) {
    fPx[i] = px[i];
    fPy[i] = py[i];
    fPz[i] = pz[i];
    fd0[i] = d0[i];
  }
}
//--------------------------------------------------------------------------
AliAODRecoDecay::AliAODRecoDecay(AliAODVertex *vtx2,Int_t nprongs,
				 Short_t charge,
				 Double_t *d0) :
  AliVTrack(),
  fSecondaryVtx(vtx2),
  fOwnSecondaryVtx(0x0),
  fCharge(charge),
  fNProngs(nprongs), fNDCA(0), fNPID(0),
  fPx(0x0), fPy(0x0), fPz(0x0),
  fd0(0x0),
  fDCA(0x0),
  fPID(0x0), 
  fEventNumber(-1),fRunNumber(-1)
{
  //
  // Constructor with AliAODVertex for decay vertex and without prongs momenta
  //

  fd0 = new Double_t[GetNProngs()];
  for(Int_t i=0; i<GetNProngs(); i++) fd0[i] = d0[i];
}
//--------------------------------------------------------------------------
AliAODRecoDecay::AliAODRecoDecay(const AliAODRecoDecay &source) :
  AliVTrack(source),
  fSecondaryVtx(source.fSecondaryVtx),
  fOwnSecondaryVtx(source.fOwnSecondaryVtx),
  fCharge(source.fCharge),
  fNProngs(source.fNProngs), fNDCA(source.fNDCA), fNPID(source.fNPID),
  fPx(0x0), fPy(0x0), fPz(0x0),
  fd0(0x0), 
  fDCA(0x0),
  fPID(0x0), 
  fEventNumber(source.fEventNumber),fRunNumber(source.fRunNumber)
{
  //
  // Copy constructor
  //
  if(source.GetNProngs()>0) {
    fd0 = new Double32_t[GetNProngs()];
    memcpy(fd0,source.fd0,GetNProngs()*sizeof(Double32_t));
    if(source.fPx) {
      fPx = new Double32_t[GetNProngs()];
      fPy = new Double32_t[GetNProngs()];
      fPz = new Double32_t[GetNProngs()];
      memcpy(fPx,source.fPx,GetNProngs()*sizeof(Double32_t));
      memcpy(fPy,source.fPy,GetNProngs()*sizeof(Double32_t));
      memcpy(fPz,source.fPz,GetNProngs()*sizeof(Double32_t));
    }
    if(source.fPID) {
      fPID = new Double32_t[fNPID];
      memcpy(fPID,source.fPID,fNPID*sizeof(Double32_t));
    }
    if(source.fDCA) {
      fDCA = new Double32_t[fNDCA];
      memcpy(fDCA,source.fDCA,fNDCA*sizeof(Double32_t));
    }
  }
}
//--------------------------------------------------------------------------
AliAODRecoDecay &AliAODRecoDecay::operator=(const AliAODRecoDecay &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;
  fSecondaryVtx = source.fSecondaryVtx;
  fOwnSecondaryVtx = source.fOwnSecondaryVtx;
  fCharge = source.fCharge;
  fNProngs = source.fNProngs;
  fNDCA = source.fNDCA;
  fNPID = source.fNPID;
  fEventNumber = source.fEventNumber;
  fRunNumber = source.fRunNumber;
  if(source.GetNProngs()>0) {
    if(fd0)delete [] fd0; 
    fd0 = new Double32_t[GetNProngs()];
    memcpy(fd0,source.fd0,GetNProngs()*sizeof(Double32_t));
    if(source.fPx) {
      if(fPx) delete [] fPx; 
      fPx = new Double32_t[GetNProngs()];
      if(fPy) delete [] fPy; 
      fPy = new Double32_t[GetNProngs()];
      if(fPz) delete [] fPz; 
      fPz = new Double32_t[GetNProngs()];
      memcpy(fPx,source.fPx,GetNProngs()*sizeof(Double32_t));
      memcpy(fPy,source.fPy,GetNProngs()*sizeof(Double32_t));
      memcpy(fPz,source.fPz,GetNProngs()*sizeof(Double32_t));
    }
    if(source.fPID) {
      if(fPID) delete [] fPID; 
      fPID = new Double32_t[fNPID];
      memcpy(fPID,source.fPID,fNPID*sizeof(Double32_t));
    }
    if(source.fDCA) {
      if(fDCA) delete [] fDCA; 
      fDCA = new Double32_t[fNDCA];
      memcpy(fDCA,source.fDCA,fNDCA*sizeof(Double32_t));
    }
  }
  return *this;
}
//--------------------------------------------------------------------------
AliAODRecoDecay::~AliAODRecoDecay() {
  //  
  // Default Destructor
  //
  if(fPx) { delete [] fPx; fPx=NULL; } 
  if(fPy) { delete [] fPy; fPy=NULL; }
  if(fPz) { delete [] fPz; fPz=NULL; }
  if(fd0) { delete [] fd0; fd0=NULL; }
  if(fPID) { delete [] fPID; fPID=NULL; }
  if(fDCA) { delete [] fDCA; fDCA=NULL; }
  if(fOwnSecondaryVtx) { delete fOwnSecondaryVtx; fOwnSecondaryVtx=NULL; }
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::Alpha() const 
{
  //
  // Armenteros-Podolanski alpha for 2-prong decays
  //
  if(GetNProngs()!=2) {
    printf("Can be called only for 2-prong decays");
    return (Double_t)-99999.;
  }
  return 1.-2./(1.+QlProng(0)/QlProng(1));
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::DecayLength(Double_t point[3]) const 
{
  //
  // Decay length assuming it is produced at "point" [cm]
  //
  return TMath::Sqrt((point[0]-GetSecVtxX())
		    *(point[0]-GetSecVtxX())
		    +(point[1]-GetSecVtxY())
		    *(point[1]-GetSecVtxY())
		    +(point[2]-GetSecVtxZ())
		    *(point[2]-GetSecVtxZ()));  
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::DecayLengthXY(Double_t point[3]) const 
{
  //
  // Decay length in XY assuming it is produced at "point" [cm]
  //
  return TMath::Sqrt((point[0]-GetSecVtxX())
		    *(point[0]-GetSecVtxX())
		    +(point[1]-GetSecVtxY())
		    *(point[1]-GetSecVtxY()));  
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::CosPointingAngle(Double_t point[3]) const 
{
  //
  // Cosine of pointing angle in space assuming it is produced at "point"
  //
  TVector3 mom(Px(),Py(),Pz());
  TVector3 fline(GetSecVtxX()-point[0],
		 GetSecVtxY()-point[1],
		 GetSecVtxZ()-point[2]);

  Double_t pta = mom.Angle(fline);

  return TMath::Cos(pta); 
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::CosPointingAngleXY(Double_t point[3]) const 
{
  //
  // Cosine of pointing angle in transverse plane assuming it is produced 
  // at "point"
  //
  TVector3 momXY(Px(),Py(),0.);
  TVector3 flineXY(GetSecVtxX()-point[0],
		   GetSecVtxY()-point[1],
		   0.);

  Double_t ptaXY = momXY.Angle(flineXY);

  return TMath::Cos(ptaXY); 
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::CosThetaStar(Int_t ip,UInt_t pdgvtx,UInt_t pdgprong0,UInt_t pdgprong1) const 
{
  //
  // Only for 2-prong decays: 
  // Cosine of decay angle (theta*) in the rest frame of the mother particle
  // for prong ip (0 or 1) with mass hypotheses pdgvtx for mother particle,
  // pdgprong0 for prong 0 and pdgprong1 for prong1
  //
  if(GetNProngs()!=2) {
    printf("Can be called only for 2-prong decays");
    return (Double_t)-99999.;
  }
  Double_t massvtx = TDatabasePDG::Instance()->GetParticle(pdgvtx)->Mass();
  Double_t massp[2];
  massp[0] = TDatabasePDG::Instance()->GetParticle(pdgprong0)->Mass();
  massp[1] = TDatabasePDG::Instance()->GetParticle(pdgprong1)->Mass();

  Double_t pStar = TMath::Sqrt(TMath::Power(massvtx*massvtx-massp[0]*massp[0]-massp[1]*massp[1],2.)-4.*massp[0]*massp[0]*massp[1]*massp[1])/(2.*massvtx);

  Double_t beta = P()/E(pdgvtx);
  Double_t gamma = E(pdgvtx)/massvtx;

  Double_t cts = (QlProng(ip)/gamma-beta*TMath::Sqrt(pStar*pStar+massp[ip]*massp[ip]))/pStar;

  return cts;
}
//---------------------------------------------------------------------------
Double_t AliAODRecoDecay::Ct(UInt_t pdg,Double_t point[3]) const
{
  //
  // Decay time * c assuming it is produced at "point" [cm]
  //
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  return DecayLength(point)*mass/P();
}
//---------------------------------------------------------------------------
Double_t AliAODRecoDecay::E(UInt_t pdg) const 
{
  //
  // Energy
  //
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  return TMath::Sqrt(mass*mass+P()*P());
}
//---------------------------------------------------------------------------
Double_t AliAODRecoDecay::EProng(Int_t ip,UInt_t pdg) const 
{
  //
  // Energy of ip-th prong 
  //
  Double_t mass = TDatabasePDG::Instance()->GetParticle(pdg)->Mass();
  return TMath::Sqrt(mass*mass+PProng(ip)*PProng(ip));
}
//--------------------------------------------------------------------------
Bool_t AliAODRecoDecay::GetCovarianceXYZPxPyPz(Double_t cv[21]) const {
  //
  // This function returns the global covariance matrix of the track params
  // 
  // Cov(x,x) ... :   cv[0]
  // Cov(y,x) ... :   cv[1]  cv[2]
  // Cov(z,x) ... :   cv[3]  cv[4]  cv[5]
  // Cov(px,x)... :   cv[6]  cv[7]  cv[8]  cv[9]
  // Cov(py,x)... :   cv[10] cv[11] cv[12] cv[13] cv[14]
  // Cov(pz,x)... :   cv[15] cv[16] cv[17] cv[18] cv[19] cv[20]
  //
  // For XYZ we take the cov of the vertex, for PxPyPz we take the 
  // sum of the covs of PxPyPz from the daughters, for the moment 
  // we set the cov between position and momentum as the sum of 
  // the same cov from the daughters.
  //

  Int_t j;
  for(j=0;j<21;j++) cv[j]=0.;

  if(!GetNDaughters()) {
    AliError("No daughters available");
    return kFALSE;
  }

  Double_t v[6];
  AliAODVertex *secv=GetSecondaryVtx();
  if(!secv) {
    AliError("Vertex covariance matrix not available");
    return kFALSE;
  }
  if(!secv->GetCovMatrix(v)) {
    AliError("Vertex covariance matrix not available");
    return kFALSE;
  }

  Double_t p[21]; for(j=0;j<21;j++) p[j]=0.;
  Bool_t error=kFALSE;
  for(Int_t i=1; i<GetNDaughters(); i++) {
    AliVTrack *daugh = (AliVTrack*)GetDaughter(i);
    Double_t dcov[21];
    if(!daugh->GetCovarianceXYZPxPyPz(dcov)) error=kTRUE;
    for(j=0;j<21;j++) p[j] += dcov[j];
  }
  if(error) {
    AliError("No covariance for at least one daughter")
    return kFALSE;
  }

  for(j=0; j<21; j++) {
    if(j<6) {
      cv[j] = v[j];
    } else {
      cv[j] = p[j];
    }
  }

  return kTRUE;
}
//----------------------------------------------------------------------------
UChar_t  AliAODRecoDecay::GetITSClusterMap() const {
  //
  // We take the logical AND of the daughters cluster maps 
  // (only if all daughters have the bit for given layer, we set the bit)
  //
  UChar_t map=0;

  if(!GetNDaughters()) {
    AliError("No daughters available");
    return map;
  }

  for(Int_t l=0; l<12; l++) { // loop on ITS layers (from here we cannot know how many they are; let's put 12 to be conservative)
    Int_t bit = 1;
    for(Int_t i=0; i<GetNDaughters(); i++) {
      AliVTrack *daugh = (AliVTrack*)GetDaughter(i);
      if(!TESTBIT(daugh->GetITSClusterMap(),l)) bit=0; 
    }
    if(bit) SETBIT(map,l);
  }

  return map;
}
//--------------------------------------------------------------------------
ULong_t AliAODRecoDecay::GetStatus() const {
  // 
  // Same as for ITSClusterMap
  //
  ULong_t status=0;

  if(!GetNDaughters()) {
    AliError("No daughters available");
    return status;
  }

  AliVTrack *daugh0 = (AliVTrack*)GetDaughter(0);
  status = status&(daugh0->GetStatus());

  for(Int_t i=1; i<GetNDaughters(); i++) {
    AliVTrack *daugh = (AliVTrack*)GetDaughter(i);
    status = status&(daugh->GetStatus());
  }

  return status;
}
//--------------------------------------------------------------------------
Double_t AliAODRecoDecay::ImpParXY(Double_t point[3]) const 
{
  //
  // Impact parameter in the bending plane of the particle 
  // w.r.t. to "point"
  //
  Double_t k = -(GetSecVtxX()-point[0])*Px()-(GetSecVtxY()-point[1])*Py();
  k /= Pt()*Pt();
  Double_t dx = GetSecVtxX()-point[0]+k*Px();
  Double_t dy = GetSecVtxY()-point[1]+k*Py();
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy);
  TVector3 mom(Px(),Py(),Pz());
  TVector3 fline(GetSecVtxX()-point[0],
		 GetSecVtxY()-point[1],
		 GetSecVtxZ()-point[2]);
  TVector3 cross = mom.Cross(fline);
  return (cross.Z()>0. ? absImpPar : -absImpPar);
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::InvMass(Int_t npdg,UInt_t *pdg) const 
{
  //
  // Invariant mass for prongs mass hypotheses in pdg array
  //
  if(GetNProngs()!=npdg) {
    printf("npdg != GetNProngs()");
    return (Double_t)-99999.;
  }
  Double_t energysum = 0.;

  for(Int_t i=0; i<GetNProngs(); i++) {
    energysum += EProng(i,pdg[i]);
  }

  Double_t mass = TMath::Sqrt(energysum*energysum-P()*P());

  return mass;
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::InvMass2Prongs(Int_t ip1,Int_t ip2,
				      UInt_t pdg1,UInt_t pdg2) const 
{
  //
  // 2-prong(ip1,ip2) invariant mass for prongs mass hypotheses in pdg1,2
  //
  Double_t energysum = EProng(ip1,pdg1) + EProng(ip2,pdg2);
  Double_t psum2 = (PxProng(ip1)+PxProng(ip2))*(PxProng(ip1)+PxProng(ip2))
                  +(PyProng(ip1)+PyProng(ip2))*(PyProng(ip1)+PyProng(ip2))
                  +(PzProng(ip1)+PzProng(ip2))*(PzProng(ip1)+PzProng(ip2));
  Double_t mass = TMath::Sqrt(energysum*energysum-psum2);

  return mass;
}
//----------------------------------------------------------------------------
Int_t AliAODRecoDecay::MatchToMC(Int_t pdgabs,TClonesArray *mcArray) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle
  // 

  if(!GetNDaughters()) {
    AliError("No daughters available");
    return -1;
  }
  
  Int_t *dgLabels = new Int_t[GetNDaughters()];

  // loop on daughters and write the labels
  for(Int_t i=0; i<GetNDaughters(); i++) {
    AliVTrack *trk = (AliVTrack*)GetDaughter(i);
    dgLabels[i] = trk->GetLabel();
  }

  Int_t labMother = MatchToMC(pdgabs,mcArray,dgLabels);

  delete [] dgLabels; dgLabels=NULL;

  return labMother;
}
//----------------------------------------------------------------------------
Int_t AliAODRecoDecay::MatchToMC(Int_t pdgabs,TClonesArray *mcArray,
				 Int_t *dgLabels) const
{
  //
  // Check if this candidate is matched to a MC signal
  // If no, return -1
  // If yes, return label (>=0) of the AliAODMCParticle
  // 

  if(!GetNDaughters()) {
    AliError("No daughters available");
    return -1;
  }
  
  Int_t *labMom = new Int_t[GetNDaughters()];
  Int_t i,lab,labMother,pdgMother;
  AliAODMCParticle *part=0;

  // loop on daughter labels
  for(i=0; i<GetNDaughters(); i++) {
    labMom[i]=-1;
    lab = dgLabels[i];
    if(lab<0) {
      printf("daughter with negative label\n");
      continue;
    }
    part = (AliAODMCParticle*)mcArray->At(lab);
    if(!part) { 
      printf("no MC particle\n");
      continue;
    }
    while(part->GetMother()>=0) {
      labMother=part->GetMother();
      part = (AliAODMCParticle*)mcArray->At(labMother);
      if(!part) {
	printf("no MC mother particle\n");
	break;
      }
      pdgMother = TMath::Abs(part->GetPdgCode());
      if(pdgMother==pdgabs) {
	labMom[i]=labMother;
	break;
      }
    }
  }

  // check if the candidate is signal
  Bool_t isSignal=kTRUE;
  labMother=labMom[0];
  // all labels have to be the same and !=-1
  for(i=0; i<GetNDaughters(); i++) {
    if(labMom[i]==-1 || labMom[i]!=labMother) isSignal=kFALSE;
  }

  delete [] labMom; labMom=NULL;

  if(!isSignal) return -1;

  // check that the mother decayed in <GetNDaughters()> prongs
  part = (AliAODMCParticle*)mcArray->At(labMother);
  Int_t ndg = TMath::Abs(part->GetDaughter(1)-part->GetDaughter(0))+1;

  if(ndg!=GetNDaughters()) return -1;

  return labMother;
}
//---------------------------------------------------------------------------
void AliAODRecoDecay::Print(Option_t* /*option*/) const 
{
  //
  // Print some information
  //
  printf("AliAODRecoDecay with %d prongs\n",GetNProngs());
  printf("Secondary Vertex: (%f, %f, %f)\n",GetSecVtxX(),GetSecVtxY(),GetSecVtxZ());

  return;
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::ProngsRelAngle(Int_t ip1,Int_t ip2) const 
{
  //
  // Relative angle between two prongs
  //
  TVector3 momA(PxProng(ip1),PyProng(ip1),PzProng(ip1));
  TVector3 momB(PxProng(ip2),PyProng(ip2),PzProng(ip2));

  Double_t angle = momA.Angle(momB);

  return angle; 
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::QlProng(Int_t ip) const 
{
  //
  // Longitudinal momentum of prong w.r.t. to total momentum
  //
  TVector3 mom(PxProng(ip),PyProng(ip),PzProng(ip));
  TVector3 momTot(Px(),Py(),Pz());

  return mom.Dot(momTot)/momTot.Mag();
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::QtProng(Int_t ip) const 
{
  //
  // Transverse momentum of prong w.r.t. to total momentum  
  //
  TVector3 mom(PxProng(ip),PyProng(ip),PzProng(ip));
  TVector3 momTot(Px(),Py(),Pz());

  return mom.Perp(momTot);
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::QlProngFlightLine(Int_t ip,Double_t point[3]) const 
{
  //
  // Longitudinal momentum of prong w.r.t. to flight line between "point"
  // and fSecondaryVtx
  //
  TVector3 mom(PxProng(ip),PyProng(ip),PzProng(ip));
  TVector3 fline(GetSecVtxX()-point[0],
		 GetSecVtxY()-point[1],
		 GetSecVtxZ()-point[2]);

  return mom.Dot(fline)/fline.Mag();
}
//----------------------------------------------------------------------------
Double_t AliAODRecoDecay::QtProngFlightLine(Int_t ip,Double_t point[3]) const 
{
  //
  // Transverse momentum of prong w.r.t. to flight line between "point" and 
  // fSecondaryVtx 
  //
  TVector3 mom(PxProng(ip),PyProng(ip),PzProng(ip));
  TVector3 fline(GetSecVtxX()-point[0],
		 GetSecVtxY()-point[1],
		 GetSecVtxZ()-point[2]);

  return mom.Perp(fline);
}
//--------------------------------------------------------------------------
