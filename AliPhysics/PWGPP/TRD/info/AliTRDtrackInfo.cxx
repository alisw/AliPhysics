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

/* $Id: AliTRDtrackInfo.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Reconstruction QA                                                     //
//                                                                        //
//  Authors:                                                              //
//    Alex Bercuci <A.Bercuci@gsi.de>                                     //
//    Markus Fasel <M.Fasel@gsi.de>                                       //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <TDatabasePDG.h>
#include <TPDGCode.h>
#include <TVectorT.h>
#include <TGeoMatrix.h>

#include <AliTrackReference.h>
#include <AliTrackPointArray.h>
#include <AliExternalTrackParam.h>

#include <AliTRDseedV1.h>
#include <AliTRDtrackV1.h>
#include <AliTRDgeometry.h>
#include <AliTRDtrackerV1.h>

#include "AliTRDtrackInfo.h"

ClassImp(AliTRDtrackInfo)
ClassImp(AliTRDtrackInfo::AliMCinfo)
ClassImp(AliTRDtrackInfo::AliESDinfo)
Double_t AliTRDtrackInfo::AliMCinfo::fgKalmanStep = 2.;
Bool_t AliTRDtrackInfo::AliMCinfo::fgKalmanUpdate = kTRUE;

//___________________________________________________
AliTRDtrackInfo::AliTRDtrackInfo():
  TObject()
  ,fNClusters(0)
  ,fTRDtrack(NULL)
  ,fMC(NULL)
  ,fESD()
{
  //
  // Default constructor
  //
}


//___________________________________________________
AliTRDtrackInfo::AliTRDtrackInfo(const AliTRDtrackInfo &trdInfo):
  TObject((const TObject&)trdInfo)  
  ,fNClusters(trdInfo.fNClusters)
  ,fTRDtrack(NULL)
  ,fMC(NULL)
  ,fESD(trdInfo.fESD)
{
  //
  // copy Entries
  //

  if(trdInfo.fMC) fMC = new AliMCinfo(*trdInfo.fMC);
  SetTrack(trdInfo.fTRDtrack);
}

//___________________________________________________
AliTRDtrackInfo::AliMCinfo::AliMCinfo()
  :fLabel(0)
  ,fTRDlabel(0)
  ,fPDG(0)
  ,fNTrackRefs(0)
  ,fEta(-999.)
  ,fPhi(-999.)
  ,fPt(-1.)
{
  // Set 0-Pointers
  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
}

//___________________________________________________
AliTRDtrackInfo::AliMCinfo::AliMCinfo(const AliMCinfo &mc)
  :fLabel(mc.fLabel)
  ,fTRDlabel(mc.fTRDlabel)
  ,fPDG(mc.fPDG)
  ,fNTrackRefs(mc.fNTrackRefs)
  ,fEta(mc.fEta)
  ,fPhi(mc.fPhi)
  ,fPt(mc.fPt)
{
  //
  // Constructor
  //

  memset(fTrackRefs, 0, sizeof(AliTrackReference *) * 12);
  for(Int_t ien = 0; ien < 12; ien++){
    if(mc.fTrackRefs[ien])
      fTrackRefs[ien] = new AliTrackReference(*(mc.fTrackRefs[ien]));
  }
}

//________________________________________________________
Int_t AliTRDtrackInfo::AliMCinfo::GetPID() const
{
// Translate pdg code to PID index

  switch(fPDG){
  case kElectron: 
  case kPositron: return AliPID::kElectron;  
  case kMuonPlus:
  case kMuonMinus: return AliPID::kMuon;  
  case kPiPlus: 
  case kPiMinus: return AliPID::kPion;  
  case kKPlus: 
  case kKMinus: return AliPID::kKaon;
  case kProton: 
  case kProtonBar: return AliPID::kProton;
  } 
  return -1;
}

//___________________________________________________
AliTRDtrackInfo::AliESDinfo::AliESDinfo()
  :fSteer(0)
  ,fId(-1)
  ,fStatus(0)
  ,fKinkIndex(0)
  ,fTPCncls(0)
  ,fTPCdedx(0.)
  ,fTOFbeta(0.)
  ,fTOFbc(0)
  ,fTRDpidQuality(0)
  ,fTRDnSlices(0)
  ,fPt(0.)
  ,fPhi(-999.)
  ,fEta(-999.)
  ,fTRDslices(NULL)
  ,fOP(NULL)
  ,fTPCout(NULL)
  ,fITSout(NULL)
  ,fTPArray(NULL)
{
  //
  // Constructor
  //

  memset(fTRDr, 0, AliPID::kSPECIES*sizeof(Double32_t));
  memset(fTRDv0pid, 0, AliPID::kSPECIES*sizeof(Int_t));
}

//___________________________________________________
AliTRDtrackInfo::AliESDinfo::AliESDinfo(const AliESDinfo &esd)
  :fSteer(esd.fSteer)
  ,fId(esd.fId)
  ,fStatus(esd.fStatus)
  ,fKinkIndex(esd.fKinkIndex)
  ,fTPCncls(esd.fTPCncls)
  ,fTPCdedx(esd.fTPCdedx)
  ,fTOFbeta(esd.fTOFbeta)
  ,fTOFbc(esd.fTOFbc)
  ,fTRDpidQuality(esd.fTRDpidQuality)
  ,fTRDnSlices(esd.fTRDnSlices)
  ,fPt(esd.fPt)
  ,fPhi(esd.fPhi)
  ,fEta(esd.fEta)
  ,fTRDslices(NULL)
  ,fOP(NULL)
  ,fTPCout(NULL)
  ,fITSout(NULL)
  ,fTPArray(NULL)
{
  //
  // Constructor
  //

  memcpy(fTRDr, esd.fTRDr, AliPID::kSPECIES*sizeof(Double32_t));
  memcpy(fTRDv0pid, esd.fTRDv0pid, AliPID::kSPECIES*sizeof(Int_t));

  if(fTRDnSlices){
    fTRDslices = new Double32_t[fTRDnSlices];
    memcpy(fTRDslices, esd.fTRDslices, fTRDnSlices*sizeof(Double32_t));
  }
  SetOuterParam(esd.fOP);
  SetTPCoutParam(esd.fTPCout);
  SetITSoutParam(esd.fITSout);
  SetTrackPointArray(esd.fTPArray);
}


//___________________________________________________
AliTRDtrackInfo::~AliTRDtrackInfo()
{
  //
  // Destructor
  //

  if(fMC) delete fMC;
  if(fTRDtrack) delete fTRDtrack;
}

//___________________________________________________
AliTRDtrackInfo::AliMCinfo::~AliMCinfo()
{
  //
  // Destructor
  //

  fNTrackRefs = 0;
  for(Int_t ien = 0; ien < 12; ien++){
    if(fTrackRefs[ien]) delete fTrackRefs[ien];
    fTrackRefs[ien] = NULL;
  }
}

//___________________________________________________
AliTRDtrackInfo::AliESDinfo::~AliESDinfo()
{
  //
  // Destructor
  //

  if(fTRDnSlices){
    delete [] fTRDslices;
    fTRDslices = NULL;
    fTRDnSlices = 0;
  }
  if(fOP) delete fOP; fOP = NULL;
  if(fTPCout) delete fTPCout; fTPCout = NULL;
  if(fITSout) delete fITSout; fITSout = NULL;
  if(fTPArray) delete fTPArray;
}

//___________________________________________________
void AliTRDtrackInfo::AliESDinfo::Delete(const Option_t *){
  //
  // Delete Pointer members 
  // 
  if(fTRDnSlices){
    delete [] fTRDslices;
    fTRDslices = NULL;
    fTRDnSlices = 0;
  }
  if(fOP) delete fOP; fOP = NULL;
  if(fTPCout) delete fTPCout; fTPCout = NULL;
  if(fITSout) delete fITSout; fITSout = NULL;
  if(fTPArray) delete fTPArray; fTPArray = NULL;
}


//___________________________________________________
AliTRDtrackInfo& AliTRDtrackInfo::operator=(const AliTRDtrackInfo &trdInfo)
{
  //
  // = Operator
  //
  if(this == &trdInfo) return *this;

  fNClusters  = trdInfo.fNClusters;
  fESD        = trdInfo.fESD;

  if(trdInfo.fMC){
    if(!fMC) fMC = new AliMCinfo(*trdInfo.fMC);
    else{
      fMC->~AliMCinfo();
      new(fMC) AliMCinfo(*trdInfo.fMC);
    }
  } else {if(fMC) delete fMC; fMC = NULL;}

  SetTrack(trdInfo.fTRDtrack);

  return *this;
}

//___________________________________________________
AliTRDtrackInfo::AliMCinfo& AliTRDtrackInfo::AliMCinfo::operator=(const AliMCinfo &mc)
{
  //
  // Assignment operator
  //

  if(this == &mc) return *this;
  fLabel      = mc.fLabel;
  fTRDlabel   = mc.fTRDlabel;
  fPDG        = mc.fPDG;
  fNTrackRefs = mc.fNTrackRefs;

  AliTrackReference **itr = &fTrackRefs[0];
  AliTrackReference* const *jtr = &mc.fTrackRefs[0];
  for(Int_t ien = 0; ien < 2*AliTRDgeometry::kNlayer; ien++, itr++, jtr++){
    if((*jtr)){
      if(!(*itr)) (*itr) = new AliTrackReference(*(*jtr));
      else{
        (*itr)->~AliTrackReference();
        new(&(*itr)) AliTrackReference(*(*jtr));
      }
    } else {if((*itr)) delete (*itr); (*itr) = NULL;}
  }
  return *this;
}

//___________________________________________________
AliTRDtrackInfo::AliESDinfo& AliTRDtrackInfo::AliESDinfo::operator=(const AliESDinfo &esd)
{
  //
  // Assignment operator
  //

  if(this == &esd) return *this;
  fSteer       = esd.fSteer;
  fId          = esd.fId;
  fStatus      = esd.fStatus;
  fKinkIndex   = esd.fKinkIndex;
  fTPCncls     = esd.fTPCncls;
  fTPCdedx     = esd.fTPCdedx;
  fTOFbeta     = esd.fTOFbeta;
  fTOFbc       = esd.fTOFbc;
  fTRDpidQuality= esd.fTRDpidQuality;
  fTRDnSlices  = esd.fTRDnSlices;
  fPt          = esd.fPt;
  fPhi         = esd.fPhi;
  fEta         = esd.fEta;
  
  memcpy(fTRDr, esd.fTRDr, AliPID::kSPECIES*sizeof(Double32_t));
  memcpy(fTRDv0pid, esd.fTRDv0pid, AliPID::kSPECIES*sizeof(Int_t));

  if(fTRDnSlices){
    if(!fTRDslices) fTRDslices = new Double32_t[fTRDnSlices];
    memcpy(fTRDslices, esd.fTRDslices, fTRDnSlices*sizeof(Double32_t));
  }
  SetOuterParam(esd.fOP);
  SetTPCoutParam(esd.fTPCout);
  SetITSoutParam(esd.fITSout);
  SetTrackPointArray(esd.fTPArray);

  return *this;
}

//___________________________________________________
void AliTRDtrackInfo::Delete(const Option_t *)
{
  //
  // Delete
  //

  AliDebug(2, Form("track[%p] mc[%p]", (void*)fTRDtrack, (void*)fMC));
  fNClusters  = 0;
  if(fMC) delete fMC; fMC = NULL;
  fESD.Delete(NULL);
  if(fTRDtrack) delete fTRDtrack; fTRDtrack = NULL;
}

//___________________________________________________
void AliTRDtrackInfo::SetTrack(const AliTRDtrackV1 *track)
{
  //
  // Set the TRD track
  //

  if(track){
    if(!fTRDtrack) fTRDtrack = new AliTRDtrackV1(*track);
    else{
      fTRDtrack->~AliTRDtrackV1();
      new(fTRDtrack) AliTRDtrackV1(*track);
    }
    if(track->IsOwner()) fTRDtrack->SetOwner();
  } else {
    if(fTRDtrack) delete fTRDtrack; fTRDtrack = NULL;
  }
}

//___________________________________________________
void AliTRDtrackInfo::AliESDinfo::SetTrackPointArray(const AliTrackPointArray *tps)
{
  //
  // Set the track point array for alignment task
  //


  if(tps){
    if(!fTPArray) fTPArray = new AliTrackPointArray(*tps);
    else{
      fTPArray->~AliTrackPointArray();
      new(fTPArray) AliTrackPointArray(*tps);
    }
  } else {
    if(fTPArray) delete fTPArray; fTPArray = NULL;
  }
}

//___________________________________________________
void AliTRDtrackInfo::AddTrackRef(const AliTrackReference *tref)
{
  //
  // Add track reference
  //

  if(fMC->fNTrackRefs >= 2*AliTRDgeometry::kNlayer){ 
    SetCurved();
    return;
  }
  // Make a copy for the object in order to avoid ownership problems
  fMC->fTrackRefs[fMC->fNTrackRefs++] = new AliTrackReference(*tref);
}

//___________________________________________________
AliTrackReference* AliTRDtrackInfo::GetTrackRef(Int_t idx) const
{
//
// Returns a track reference
//
  if(!fMC) return NULL;
  return (idx>=0 && idx < 12) ? fMC->fTrackRefs[idx] : NULL;
}

//___________________________________________________
AliTrackReference* AliTRDtrackInfo::GetTrackRef(const AliTRDseedV1* const tracklet) const
{
//
// Returns a track reference
//
  if(!fMC) return NULL;
  Double_t cw = AliTRDgeometry::CamHght() + AliTRDgeometry::CdrHght();
  AliTrackReference * const* jtr = &(fMC->fTrackRefs[0]);
  for(Int_t itr = 0; itr < fMC->fNTrackRefs; itr++, ++jtr){
    if(!(*jtr)) break;   
    if(TMath::Abs(tracklet->GetX0() - (*jtr)->LocalX()) < cw) return (*jtr);
  }
  return NULL;
}

//___________________________________________________
Int_t AliTRDtrackInfo::GetNumberOfClusters() const
{
  //
  // Returns the number of clusters
  //

  Int_t n = 0;
  if(!fTRDtrack) return 0;
  if(fTRDtrack->GetNumberOfTracklets() == 0) return n;
  AliTRDseedV1 *tracklet = NULL;
  for(Int_t ip=0; ip<AliTRDgeometry::kNlayer; ip++){
    if(!(tracklet = const_cast<AliTRDseedV1 *>(fTRDtrack->GetTracklet(ip)))) continue;
    n+=tracklet->GetN();
  }
  return n;
}


//___________________________________________________
void  AliTRDtrackInfo::AliESDinfo::SetOuterParam(const AliExternalTrackParam *op)
{
  //
  // Set outer track parameters
  //

  if(op){
    if(fOP){
      fOP->~AliExternalTrackParam();
      // RS: Constructor from VTrack was used instead of Constructor from AliExternalTrackParam
      new(fOP) AliExternalTrackParam(*op);
    } else fOP = new AliExternalTrackParam(*op);
  } else {
    if(fOP) delete fOP; fOP = NULL;
  }
}

//___________________________________________________
void  AliTRDtrackInfo::AliESDinfo::SetITSoutParam(const AliExternalTrackParam *op)
{
  //
  // Set TPCout track parameters
  //

  if(op){
    if(fITSout){
      fITSout->~AliExternalTrackParam();
      // RS: Constructor from VTrack was used instead of Constructor from AliExternalTrackParam
      new(fITSout) AliExternalTrackParam(*op);
    } else fITSout = new AliExternalTrackParam(*op);
  } else {
    if(fITSout) delete fITSout; fITSout = NULL;
  }
}

//___________________________________________________
void  AliTRDtrackInfo::AliESDinfo::SetTPCoutParam(const AliExternalTrackParam *op)
{
  //
  // Set TPCout track parameters
  //

  if(op){
    if(fTPCout){
      fTPCout->~AliExternalTrackParam();
      // RS: Constructor from VTrack was used instead of Constructor from AliExternalTrackParam
      new(fTPCout) AliExternalTrackParam(*op);
    } else fTPCout = new AliExternalTrackParam(*op);
  } else {
    if(fTPCout) delete fTPCout; fTPCout = NULL;
  }
}

//___________________________________________________
Int_t AliTRDtrackInfo::GetNTracklets() const
{
  //
  // Return the number of tracklets
  //

  if(!fTRDtrack) return 0;
  return fTRDtrack->GetNumberOfTracklets();
}

//___________________________________________________
void AliTRDtrackInfo::SetSlices(Int_t n, Double32_t *s)
{
  //
  // Set the slices
  //
  if(fESD.fTRDnSlices != n){
    fESD.fTRDnSlices = 0;
    delete [] fESD.fTRDslices;
    fESD.fTRDslices = NULL;
  }

  if(!fESD.fTRDnSlices){
    fESD.fTRDnSlices = n;
    fESD.fTRDslices = new Double32_t[fESD.fTRDnSlices];
  }

  memcpy(fESD.fTRDslices, s, n*sizeof(Double32_t));
}
 
//___________________________________________________
Bool_t AliTRDtrackInfo::AliMCinfo::GetDirections(Float_t x0, Int_t chg, TGeoHMatrix* matrix, Double_t dir[10], UChar_t &status) const
{
// Check for 2 track ref for the tracklet defined bythe radial position x0
// The "status" is a bit map and gives a more informative output in case of failure:
//   - 0 : everything is OK
//   - BIT(0) : 0 track refs found
//   - BIT(1) : 1 track reference found
//   - BIT(2) : dx <= 0 between track references
//   - BIT(3) : dx > 0 && dx < 3.7 - tangent tracks 
//
//   The return array contains the MC info in the following order
//   [0] - pt (local chamber)
//   [1] - p (local chamber)
//   [2] - eta (local chamber)
//   [3] - phi (local chamber)
//   [4] - dydx (local chamber - without alignment chamber tilts - between in and out)
//   [5] - dzdx (local chamber - without alignment chamber tilts - between in and out)
//   [6] - local y (r-phi position at Anode wire in local chamber coordinates corrected for track curvature)
//   [7] - local z (z position at Anode wire in local chamber coordinates corrected for track curvature)
//   [8] - trk y (r-phi position at Anode wire in tracking coordinates corrected for track curvature and radial shift from chamber tilt)
//   [9] - trk z (z position at Anode wire in local tracking coordinates corrected for track curvature and radial shift from chamber tilt)
  
  status = 0; memset(dir, 0, 10*sizeof(Double_t));
  Double_t cw = AliTRDgeometry::CamHght() + AliTRDgeometry::CdrHght();
  Int_t nFound = 0;
  AliTrackReference *tr[2] = {NULL, NULL};
  AliTrackReference * const* jtr = &fTrackRefs[0];
  for(Int_t itr = 0; itr < fNTrackRefs; itr++, ++jtr){
    if(!(*jtr)) break;
/*
    if(fDebugLevel>=5) printf("\t\tref[%2d] x[%6.3f]\n", itr, (*jtr)->LocalX());*/
    if(TMath::Abs(x0 - (*jtr)->LocalX()) > cw) continue;
    tr[nFound++] = (*jtr);
    if(nFound == 2) break;
  } 
  if(nFound < 2){ 
    //AliDebug(1, Form("Missing track ref x0[%6.3f] nref[%d]", x0, nFound));
    if(!nFound) SETBIT(status, 0);
    else SETBIT(status, 1);
    return kFALSE;
  }
  dir[0] = tr[1]->Pt(); Double_t pt(dir[0]); // alias
  dir[1] = tr[1]->P();
  if(dir[0] < 1.e-3) return kFALSE;
  dir[2] =  -TMath::Log(TMath::Tan(0.5 * tr[1]->Theta()));
  dir[3] =  TMath::ATan2(tr[1]->Y(), tr[1]->X());

  // CORRECT FOR TRACK CURVATURE  
  // get local chamber coordinates of TRs
  Double_t  trk_tr0[] = {tr[0]->LocalX(), tr[0]->LocalY(), tr[0]->Z()}, 
            trk_tr1[] = {tr[1]->LocalX(), tr[1]->LocalY(), tr[1]->Z()}, 
            loc_tr0[]={1.,1.,1.}, loc_tr1[]={1.,1.,1.};
  matrix->MasterToLocal(trk_tr0, loc_tr0); matrix->MasterToLocal(trk_tr1, loc_tr1);
  Double_t dx = loc_tr1[0] - loc_tr0[0];
  if(dx <= 0.){
    AliWarningGeneral("AliTRDtrackInfo::AliMCinfo::GetDirections()", Form("Track ref with wrong radial distances refX0[%6.3f] refX1[%6.3f]", tr[0]->LocalX(), tr[1]->LocalX()));
    SETBIT(status, 2);
    return kFALSE;
  }
  if(TMath::Abs(dx-AliTRDgeometry::CamHght()-AliTRDgeometry::CdrHght())>1.E-3) SETBIT(status, 3); 
  Double_t dydx = (loc_tr1[1] - loc_tr0[1]) / dx,
           dzdx = (loc_tr1[2] - loc_tr0[2]) / dx;
  dir[4] = dydx; dir[5] = dzdx; 
  
  // find center of track curvature by solving the circle equation for the 2 TRs
  Double_t  R  = pt/-0.299792458e-3/AliTracker::GetBz(), // local track Radius
            tgp= (trk_tr1[1] - trk_tr0[1])/(trk_tr1[0] - trk_tr0[0]),
            a  = .5*(trk_tr1[0] + trk_tr0[0] +tgp*(trk_tr1[1] + trk_tr0[1])),
            ap = trk_tr1[0]-a,
            bp = 1./(1.+tgp*tgp),
            b  = (trk_tr1[1]-ap*tgp)*bp, b2(b*b),
            c  = (ap*ap+trk_tr1[1]*trk_tr1[1]-R*R)*bp,
            d  = (b2>c?TMath::Sqrt(b2-c):0.),
            yc0= b+d, yc1= b-d, yc(0.), xc(0.), ycn(yc0), ycp(yc1);
  if(yc0>0){ycn=yc1; ycp=yc0;}
  if(AliTracker::GetBz()>0) yc=chg<0?ycp:ycn; 
  else yc=chg<0?ycn:ycp; 
  xc= a-yc*tgp;
  // find position on the circle for the anode wire
  Double_t w = loc_tr1[0] - AliTRDgeometry::AnodePos(),
	      wp = w*(trk_tr1[0] - trk_tr0[0])/dx,
	      trk_trAn[3] = {trk_tr1[0] - wp, 
			     0., 
			     trk_tr1[2]-wp*(trk_tr1[2] - trk_tr0[2])/(trk_tr1[0] - trk_tr0[0])},
	      dxAnTrk= trk_trAn[0] - xc,
	      dyAnTrk= 0,
	      yAnTrk0= yc,
	      yAnTrk1= yc;
	      dyAnTrk = R*R - dxAnTrk*dxAnTrk; 
	      if (dyAnTrk>1e-16) {
		dyAnTrk = TMath::Sqrt(dyAnTrk);
		yAnTrk0 -= dyAnTrk; 
		yAnTrk1 += dyAnTrk;
	      }
	      else dyAnTrk = 0.;
	      trk_trAn[1]=TMath::Abs(yAnTrk0-trk_tr1[1])<TMath::Abs(yAnTrk1-trk_tr1[1])?yAnTrk0:yAnTrk1;         
  // go back to local coordinates
  matrix->MasterToLocal(trk_trAn, loc_tr1);
  Double_t dxTRD  = AliTRDgeometry::AnodePos() - loc_tr1[0];
  // compute radial position in tracking coordinates corresponding to local anode wire 
  Double_t dxTrk(dxTRD*TMath::Sqrt(1. + dydx*dydx)/TMath::Sqrt(1. + tgp*tgp));
  trk_trAn[0] += dxTrk;          
  trk_trAn[1] += dxTrk*tgp; 
  trk_trAn[2] += dxTrk*(trk_tr1[2] - trk_tr0[2])/(trk_tr1[0] - trk_tr0[0]); 
  // get coordinates also in local coordinates
  matrix->MasterToLocal(trk_trAn, loc_tr1);
  dir[6] = loc_tr1[1];  // y local @ anode
  dir[7] = loc_tr1[2];  // z local @ anode
  dir[8] = trk_trAn[1]; // y tracking @ anode
  dir[9] = trk_trAn[2];  // z tracking @ anode
  return kTRUE;
}

//___________________________________________________
Bool_t AliTRDtrackInfo::AliMCinfo::PropagateKalman(
      TVectorD *x, TVectorD *y, TVectorD *z,
      TVectorD *dx, TVectorD *dy, TVectorD *dz,
      TVectorD *pt, TVectorD *dpt, TVectorD *budget, TVectorD *c, Double_t mass) const
{
// Propagate Kalman from the first TRD track reference to 
// last one and save residuals in the y, z and pt.
// 
// This is to calibrate the dEdx and MS corrections

  if(!fNTrackRefs) return kFALSE;
  for(Int_t itr=kNTrackRefs; itr--;){
    (*x)[itr] = 0.;(*y)[itr] = 0.;(*z)[itr] = 0.;
    (*dx)[itr] = -1.; (*dy)[itr] = 100.; (*dz)[itr] = 100.; (*dpt)[itr] = 100.;
  }

  // Initialize TRD track to the first track reference
  AliTrackReference *tr(NULL);
  Int_t itr(0); while(!(tr = fTrackRefs[itr])) itr++;
  if(tr->Pt()<1.e-3) return kFALSE;

  AliTRDtrackV1 tt;
  Double_t xyz[3]={tr->X(),tr->Y(),tr->Z()};
  Double_t pxyz[3]={tr->Px(),tr->Py(),tr->Pz()};
  Double_t var[6] = {1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4};
  Double_t cov[21]={
    var[0],  0.,  0.,  0.,  0.,  0.,
         var[1],  0.,  0.,  0.,  0.,
              var[2],  0.,  0.,  0.,
                   var[3],  0.,  0.,
                        var[4],  0.,
                             var[5]
  };
  TDatabasePDG db;
  const TParticlePDG *pdg=db.GetParticle(fPDG);
  if(!pdg){
    AliWarningGeneral("AliTRDtrackInfo::AliMCinfo::PropagateKalman()", Form("PDG entry missing for code %d. References for track %d", fPDG, fNTrackRefs));
    return kFALSE;
  }
  tt.Set(xyz, pxyz, cov, Short_t(pdg->Charge()));
  if(mass<0){ // mass 0 use PDG
    tt.SetMass(pdg->Mass());
  } else { // use rec value
    tt.SetMass(mass);
  }

//  Double_t bg(tr->P()/pdg->Mass());
//  printf("\n\nNEW track PDG[%d] bg[%f] x[%f]\n", fPDG, bg, TMath::Log(bg));
  Double_t x0(tr->LocalX());
  const Double_t *cc(NULL);
  for(Int_t ip=0; itr<fNTrackRefs; itr++){
    if(!(tr = fTrackRefs[itr])) continue;
//    printf("ip[%d] det[%d]\n", ip, tr->DetectorId());
    if(!AliTRDtrackerV1::PropagateToX(tt, tr->LocalX(), fgKalmanStep)) continue;

    //if(update) ...
    tt.GetXYZ(xyz);
    (*x)[ip]   = xyz[0];
    (*y)[ip]   = xyz[1];
    (*z)[ip]   = xyz[2];
    (*dx)[ip]  = tt.GetX() - x0;
    (*dy)[ip]  = tt.GetY() - tr->LocalY();
    (*dz)[ip]  = tt.GetZ() - tr->Z();
    (*pt)[ip]  = tr->Pt();
    (*dpt)[ip] = tt.Pt()- tr->Pt();
//    printf("pt : kalman[%e] MC[%e]\n", tt.Pt(), tr->Pt());
    (*budget)[ip] = tt.GetBudget(0);
    cc = tt.GetCovariance();
//    printf("dx[%5.2f] sy[%f]\n", (*dx)[ip], TMath::Sqrt(cc[0]));
    for(Int_t ic(0), jp(ip*15); ic<15; ic++, jp++) (*c)[jp]=cc[ic];
    ip++;
  }
  return kTRUE;
}
