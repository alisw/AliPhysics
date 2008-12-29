/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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
// Class for AOD reconstructed heavy-flavour cascades
//
// Author: X-M. Zhang, zhangxm@ccnu.iop.edu.cn
/////////////////////////////////////////////////////////////

#include <TVector3.h>
#include <TDatabasePDG.h>
#include "AliAODRecoDecay.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoCascadeHF.h"

ClassImp(AliAODRecoCascadeHF)
//-----------------------------------------------------------------------------

AliAODRecoCascadeHF::AliAODRecoCascadeHF() :
  AliAODRecoDecayHF2Prong(),
  f2Prong()
{
  //
  // Default Constructor
  //
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
					 Double_t *px, Double_t *py, Double_t *pz,
					 Double_t *d0, Double_t *d0err, Double_t dca) :
  AliAODRecoDecayHF2Prong(vtx2, px, py, pz, d0, d0err, dca),
  f2Prong()
{
  //
  //  Constructor with AliAODVertex for decay vertex
  //
  SetCharge(charge);
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(AliAODVertex *vtx2, Short_t charge,
					 Double_t *d0, Double_t *d0err, Double_t dca) :
  AliAODRecoDecayHF2Prong(vtx2, d0, d0err, dca),
  f2Prong()
{
  //
  //  Constructor with decay vertex and without prongs momenta
  //
  SetCharge(charge);
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::AliAODRecoCascadeHF(const AliAODRecoCascadeHF &source) :
  AliAODRecoDecayHF2Prong(source),
  f2Prong()
{
  //
  // Copy constructor
  //
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF &AliAODRecoCascadeHF::operator=(const AliAODRecoCascadeHF &source)
{
  //
  // assignment operator
  //
  if(&source == this) return *this;

  AliAODRecoDecayHF2Prong::operator=(source);

  f2Prong = source.f2Prong;

  return *this;
}
//-----------------------------------------------------------------------------
AliAODRecoCascadeHF::~AliAODRecoCascadeHF()
{
  //
  // Default Destructor
  //
}
//-----------------------------------------------------------------------------
Double_t AliAODRecoCascadeHF::InvMassDstarKpipi() const 
{
  //
  // 3 prong invariant mass of the D0 daughters and the soft pion
  //

  Double_t px[3],py[3],pz[3];
  UInt_t pdg[3]={321,211,211};
  pdg[0] = (Charge()>0 ? 211 : 321); // positive daughter of D0
  px[0] = Get2Prong()->PxProng(0);
  py[0] = Get2Prong()->PyProng(0);
  pz[0] = Get2Prong()->PzProng(0);
  pdg[1] = (Charge()>0 ? 321 : 211); // negative daughter of D0
  px[1] = Get2Prong()->PxProng(1);
  py[1] = Get2Prong()->PyProng(1);
  pz[1] = Get2Prong()->PzProng(1);
  pdg[2] = 211; // soft pion
  px[2] = PxProng(0);
  py[2] = PyProng(0);
  pz[2] = PzProng(0);
  Short_t dummycharge=0;
  Double_t dummyd0[3]={0,0,0};
  AliAODRecoDecay *rd = new AliAODRecoDecay(0x0,3,dummycharge,px,py,pz,dummyd0);

  Double_t minv = rd->InvMass(3,pdg);

  delete rd; rd=NULL;

  return minv;
}
//-----------------------------------------------------------------------------
Bool_t AliAODRecoCascadeHF::SelectDstar(const Double_t *cutsDstar,
					const Double_t *cutsD0,
					Bool_t testD0) const
{
  //
  // cutsDstar[0] = inv. mass half width of D* [GeV]
  // cutsDstar[1] = half width of (M_Kpipi-M_D0) [GeV]
  // cutsDstar[2] = PtMin of pi_s [GeV/c]
  // cutsDstar[3] = PtMax of pi_s [GeV/c]
  // cutsDstar[4] = theta, angle between the pi_s and decay plane of the D0 [rad]
  //
  // cutsD0[0] = inv. mass half width [GeV]   
  // cutsD0[1] = dca [cm]
  // cutsD0[2] = cosThetaStar 
  // cutsD0[3] = pTK [GeV/c]
  // cutsD0[4] = pTPi [GeV/c]
  // cutsD0[5] = d0K [cm]   upper limit!
  // cutsD0[6] = d0Pi [cm]  upper limit!
  // cutsD0[7] = d0d0 [cm^2]
  // cutsD0[8] = cosThetaPoint


  // check that the D0 passes the cuts
  // (if we have a D*+, it has to pass as D0, 
  //  if we have a D*-, it has to pass as D0bar)

  if(testD0) {
    Int_t okD0=0,okD0bar=0;
    Get2Prong()->SelectD0(cutsD0,okD0,okD0bar);
    if((Charge()==+1 && !okD0) || (Charge()==-1 && !okD0bar)) return kFALSE; 
  }
 
  if( (PtProng(0)<cutsDstar[2]) || (PtProng(0)>cutsDstar[3]) ) return kFALSE;

  Double_t mDstar = TDatabasePDG::Instance()->GetParticle(413)->Mass();
  Double_t invmDstar = InvMassDstarKpipi();
  if(TMath::Abs(mDstar-invmDstar)>cutsDstar[0]) return kFALSE;

  Double_t mD0 = TDatabasePDG::Instance()->GetParticle(421)->Mass();
  if(TMath::Abs((mDstar-mD0)-DeltaInvMass())>cutsDstar[1]) return kFALSE;

  TVector3 p3Trk0(Get2Prong()->PxProng(0),Get2Prong()->PyProng(0),Get2Prong()->PzProng(0)); // from D0
  TVector3 p3Trk1(Get2Prong()->PxProng(1),Get2Prong()->PyProng(1),Get2Prong()->PzProng(1)); // from D0
  TVector3 p3Trk2(PxProng(0),PyProng(0),PzProng(0)); // pi_s

  TVector3 perp = p3Trk0.Cross(p3Trk1);
  Double_t theta = p3Trk2.Angle(perp);
  if(theta>(TMath::Pi()-theta)) theta = TMath::Pi() - theta;
  theta = TMath::Pi()/2. - theta;

  if(theta>cutsDstar[4]) return kFALSE;

  Double_t alpha = p3Trk0.Angle(p3Trk2);
  Double_t belta = p3Trk1.Angle(p3Trk2);

  Double_t cosphi01 = TMath::Cos(alpha) / TMath::Cos(theta);
  Double_t cosphi02 = TMath::Cos(belta) / TMath::Cos(theta);

  Double_t phi01 = TMath::ACos(cosphi01);
  Double_t phi02 = TMath::ACos(cosphi02);
  Double_t phi00 = p3Trk0.Angle(p3Trk1);

  if((phi01>phi00) || (phi02>phi00)) return kFALSE;
  
  return kTRUE;
}
//-----------------------------------------------------------------------------
