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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Redo primary vertex on the fly, using the diamond constraint              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include "AliESDEvent.h"
#include "AliTender.h"
#include "AliVParticle.h"
#include "AliLog.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"
#include "AliTrackFixTenderSupply.h"
#include "AliOADBTrackFix.h"

ClassImp(AliTrackFixTenderSupply)

AliTrackFixTenderSupply::AliTrackFixTenderSupply() :
fDebug(0),
  fBz(0),
  fParams(0),
  fOADBObjPath("$OADB/PWGPP/data/CorrPTInv.root"),
  fOADBObjName("CorrPTInv"),
  fOADBCont(0)
{
  // default ctor
}

//_____________________________________________________
AliTrackFixTenderSupply::AliTrackFixTenderSupply(const char *name, const AliTender *tender) :
  AliTenderSupply(name,tender),
  fDebug(0),
  fBz(0),
  fParams(0),
  fOADBObjPath("$OADB/PWGPP/data/CorrPTInv.root"),
  fOADBObjName("CorrPTInv"),
  fOADBCont(0)
{
  // named ctor
  //
}

//_____________________________________________________
AliTrackFixTenderSupply::~AliTrackFixTenderSupply()
{
  // d-tor
  delete fOADBCont;
}


//_____________________________________________________
void AliTrackFixTenderSupply::ProcessEvent()
{
  //
  // Fix track kinematics
  //
  AliESDEvent *event=fTender->GetEvent();
  if (!event) return;
  //
  if (fTender->RunChanged() && !GetRunCorrections(fTender->GetRun())) return;
  //
  int nTracks = event->GetNumberOfTracks();
  fBz = event->GetMagneticField();
  if (TMath::Abs(fBz) < kAlmost0Field) return;
  //
  const AliESDVertex *vtx = 0, *vtxTPC = 0;
  //
  vtx = event->GetPrimaryVertexTracks(); // vertex to be used for update via RelateToVertex
  if (!vtx || vtx->GetStatus()<1) {
    vtx = event->GetPrimaryVertexSPD();
    if (vtx && vtx->GetStatus()<1) vtx = 0;
  }
  vtxTPC = event->GetPrimaryVertexTPC(); // vertex to be used for update via RelateToVertexTPC
  if (vtxTPC && vtxTPC->GetStatus()<1) vtxTPC = 0;
  //  
  AliExternalTrackParam* extPar = 0;
  double xOrig = 0;
  double xyzTPCInner[3] = {0,0,0};
  for (int itr=0;itr<nTracks;itr++) {
    //
    AliESDtrack* trc = event->GetTrack(itr);
    if (!trc->IsOn(AliESDtrack::kTPCin)) continue;
    //
    double sideAfraction = GetSideAFraction(trc);
    // correct the main parameterization
    int cormode = trc->IsOn(AliESDtrack::kITSin) ? AliOADBTrackFix::kCorModeGlob : AliOADBTrackFix::kCorModeTPCInner;
    xOrig = trc->GetX();
    double xIniCor = fParams->GetXIniPtInvCorr(cormode);
    const AliExternalTrackParam* parInner = trc->GetInnerParam();
    if (!parInner) {
      AliError("Failed to extract inner param");
      continue;
    }
    parInner->GetXYZ(xyzTPCInner);
    double phi = TMath::ATan2(xyzTPCInner[1],xyzTPCInner[0]);
    if (phi<0) phi += 2*TMath::Pi();
    //
    if (fDebug>1) {
      AliInfo(Form("Tr:%4d kITSin:%d Phi=%+5.2f at X=%+7.2f | SideA fraction: %.3f",itr,trc->IsOn(AliESDtrack::kITSin),phi,parInner->GetX(),sideAfraction));
      AliInfo(Form("Main Param before corr. in mode %s, xIni:%.1f",cormode== AliOADBTrackFix::kCorModeGlob ?  "Glo":"TPC",xIniCor));
      trc->AliExternalTrackParam::Print();
    }
    //
    if (xIniCor>0) trc->PropagateTo(xIniCor,fBz);
    CorrectTrackPtInv(trc, cormode, sideAfraction, phi);
    if (xIniCor>0) {                             // full update is requested
      if (vtx) trc->RelateToVertex(vtx, fBz, kVeryBig); // redo DCA if vtx is available
      else     trc->PropagateTo(xOrig, fBz);            // otherwise bring to original point
    }
    // 
    if (fDebug>1) {
      AliInfo("Main Param after corr.");
      trc->AliExternalTrackParam::Print();
    }
    // correct TPCinner param
    if ( (extPar=(AliExternalTrackParam*)trc->GetTPCInnerParam()) ) {
      cormode = AliOADBTrackFix::kCorModeTPCInner;
      xOrig = extPar->GetX();
      xIniCor = fParams->GetXIniPtInvCorr(cormode);
      if (fDebug>1) {
	AliInfo(Form("TPCinner Param before corr. in mode %s, xIni:%.1f",cormode== AliOADBTrackFix::kCorModeGlob ?  "Glo":"TPC",xIniCor));
	extPar->AliExternalTrackParam::Print();
      }
      //
      if (xIniCor>0) extPar->PropagateTo(xIniCor,fBz);
      CorrectTrackPtInv(extPar,cormode,sideAfraction, phi);
      if (xIniCor>0) {                              // full update is requested
	if (vtxTPC) trc->RelateToVertexTPC(vtxTPC, fBz, kVeryBig);  // redo DCA if vtx is available
	else        extPar->PropagateTo(xOrig, fBz);                // otherwise bring to original point
      }
      //
      if (fDebug>1) {
	AliInfo("TPCinner Param after corr.");
	extPar->AliExternalTrackParam::Print();
      }      
    }
    //
  }
  //
}

//_____________________________________________________
void AliTrackFixTenderSupply::CorrectTrackPtInv(AliExternalTrackParam* trc, Int_t mode, double sideAfraction, double phi) const
{
  // fix track kinematics
  if (!trc) return;
  double *param = (double*)trc->GetParameter();
  param[4] += fParams->GetPtInvCorr(mode,sideAfraction,phi);
  //
}

//_____________________________________________________
Bool_t AliTrackFixTenderSupply::LoadOADBObjects()
{
  // Load OADB parameters
  TString fileName = fOADBObjPath;
  if (fileName.BeginsWith("$OADB")) fileName.ReplaceAll("$OADB",Form("%s/",AliAnalysisManager::GetOADBPath()));
  gSystem->ExpandPathName(fileName);
  AliInfo(Form("Loading correction parameters %s from %s",fOADBObjName.Data(),fileName.Data()));
  //
  fOADBCont = new AliOADBContainer("OADB");
  if (fOADBCont->InitFromFile(fileName.Data(),fOADBObjName.Data())) {
    AliError("Failed to load OADB Container");
    delete fOADBCont;
    fOADBCont = 0;
    return kFALSE;
  }
  //
  return kTRUE;
}

//_____________________________________________________
Bool_t AliTrackFixTenderSupply::GetRunCorrections(int run)
{
  // extract corrections for given run
  fParams = 0;
  if (!fOADBCont) if (!LoadOADBObjects()) return kFALSE;
  fParams = dynamic_cast<AliOADBTrackFix*>(fOADBCont->GetObject(run,"default"));
  if (!fParams) {AliError(Form("No correction parameters for found for run %d",run)); return kFALSE;}
  AliInfo(Form("Loaded correction parameters for run %d",run));
  //
  return kTRUE;
}

//_____________________________________________________
Double_t AliTrackFixTenderSupply::GetSideAFraction(const AliESDtrack* track) const
{
  // decide A-side fraction of trajectory length in TPC 
  // for the moment just a placeholder...
  //
  const AliExternalTrackParam *trIn = track->GetInnerParam();
  if (!trIn) {AliError("Failed to onbtain InnerParam"); return 0.5;}
  //
  // easy cases:
  double sideAfrac = 0.5;
  if      (trIn->GetZ()>0) { sideAfrac = 1.; if (trIn->GetTgl()>-kAlmost0) return sideAfrac;}
  else if (trIn->GetZ()<0) { sideAfrac = 0.; if (trIn->GetTgl()< kAlmost0) return sideAfrac;}
  //
  // track may cross the membrane
  // convert to usual helix param: 
  // x(t)=(dca+R)*cos(phi0) -R*cos(t+phi0);
  // y(t)=(dca+R)*sin(phi0) -R*sin(t+phi0);
  // z(t)=z0 - R*t*tgl
  //
  int sqb = fBz<0 ? -trIn->Charge() : trIn->Charge();
  double radq = TMath::Abs(trIn->GetC(fBz));
  radq = radq>kAlmost0 ? 1./radq : 1e4; // track radius
  if (sqb<0) radq = -radq;
  double cstCphi0 =  trIn->GetSnp(); // tC is current t
  double sntCphi0 = -TMath::Sqrt((1.-cstCphi0)*(1+cstCphi0));
  double xc0 = trIn->GetX() + radq*cstCphi0; // circle center
  double yc0 = trIn->GetY() + radq*sntCphi0;
  double phi0 = TMath::ATan2( yc0, xc0);
  if (sqb<0) phi0 += TMath::Pi();
  if      (phi0 > TMath::Pi()) phi0 -= 2.*TMath::Pi();
  else if (phi0 <-TMath::Pi()) phi0 += 2.*TMath::Pi();
  //
  double csphi0 = TMath::Cos(phi0);
  double snphi0 = TMath::Sin(phi0);
  double dca    = xc0*csphi0 + yc0*snphi0 - radq;
  double tC     = TMath::ATan2(sntCphi0,cstCphi0) - phi0;
  double tglR   = radq*trIn->GetTgl();
  double z0     = trIn->GetZ() + tglR*tC;
  Double_t tCross = z0/tglR;  // Z=0 is crossed at this t
  // t is changed withing tCross:2pi range for the track moving forward
  if (tCross<tC || tCross>TMath::Pi()) return sideAfrac;  // on the way forward does not cross the membrane
  //
  // find max t reached by the track in the TPC
  const double kROutTPC2 = 250.*250.;
  // find t at which the outer radius is crossed
  // r^2 = (dca+R)^2 + R^2 - 2R(R+dca)*cos(t)
  double dcaR = dca+radq;
  double rmax2 = dcaR*dcaR + radq*radq + 2*radq*dcaR; // max rad at t=pi
  double tmax = TMath::Pi();
  if (rmax2>kROutTPC2) { // traverses full TPC
    double ct = 0.5*(dcaR/radq + radq/dcaR - kROutTPC2/radq/dcaR);    
    if (TMath::Abs(ct)<1-kAlmost0) tmax = TMath::ACos(ct);
  }
  if (tmax<tCross) return sideAfrac;  // does not cross z=0
  sideAfrac = sideAfrac>0.5 ? (tCross-tC)/(tmax-tC) : 1.-(tCross-tC)/(tmax-tC);
  return sideAfrac;
}
