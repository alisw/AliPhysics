/**************************************************************************
 * Author: Panos Christakoglou.                                           *
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

/* $Id: AliProtonAnalysisBase.cxx 31056 2009-02-16 14:31:41Z pchrist $ */

//-----------------------------------------------------------------
//                 AliProtonAnalysisBase class
//   This is the class to deal with the proton analysis
//   Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------
#include <Riostream.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TF1.h>

#include <AliExternalTrackParam.h>
#include <AliESDEvent.h>
#include <AliPID.h>
#include <AliVertexerTracks.h>

class AliLog;
class AliESDVertex;

#include "AliProtonAnalysisBase.h"

ClassImp(AliProtonAnalysisBase)

//____________________________________________________________________//
AliProtonAnalysisBase::AliProtonAnalysisBase() : 
  TObject(),  fProtonAnalysisLevel("ESD"),
  fTriggerMode(kMB2), fProtonAnalysisMode(kTPC), fProtonPIDMode(kBayesian),
  fAnalysisEtaMode(kFALSE),
  fVxMax(100.), fVyMax(100.), fVzMax(100.),
  fNBinsX(0), fMinX(0), fMaxX(0),
  fNBinsY(0), fMinY(0), fMaxY(0),
  fMinTPCClusters(0), fMinITSClusters(0),
  fMaxChi2PerTPCCluster(0), fMaxChi2PerITSCluster(0),
  fMaxCov11(0), fMaxCov22(0), fMaxCov33(0), fMaxCov44(0), fMaxCov55(0),
  fMaxSigmaToVertex(0), fMaxSigmaToVertexTPC(0),
  fMaxDCAXY(0), fMaxDCAXYTPC(0),
  fMaxDCAZ(0), fMaxDCAZTPC(0),
  fMaxDCA3D(0), fMaxDCA3DTPC(0),
  fMaxConstrainChi2(0), fMinTPCdEdxPoints(0),
  fMinTPCClustersFlag(kFALSE), fMinITSClustersFlag(kFALSE),
  fMaxChi2PerTPCClusterFlag(kFALSE), fMaxChi2PerITSClusterFlag(kFALSE),
  fMaxCov11Flag(kFALSE), fMaxCov22Flag(kFALSE), 
  fMaxCov33Flag(kFALSE), fMaxCov44Flag(kFALSE), fMaxCov55Flag(kFALSE),
  fMaxSigmaToVertexFlag(kFALSE), fMaxSigmaToVertexTPCFlag(kFALSE),
  fMaxDCAXYFlag(kFALSE), fMaxDCAXYTPCFlag(kFALSE),
  fMaxDCAZFlag(kFALSE), fMaxDCAZTPCFlag(kFALSE),
  fMaxDCA3DFlag(kFALSE), fMaxDCA3DTPCFlag(kFALSE),
  fMaxConstrainChi2Flag(kFALSE),
  fITSRefitFlag(kFALSE), fTPCRefitFlag(kFALSE),
  fESDpidFlag(kFALSE), fTPCpidFlag(kFALSE),
  fPointOnITSLayer1Flag(0), fPointOnITSLayer2Flag(0),
  fPointOnITSLayer3Flag(0), fPointOnITSLayer4Flag(0),
  fPointOnITSLayer5Flag(0), fPointOnITSLayer6Flag(0),
  fMinTPCdEdxPointsFlag(kFALSE),
  fFunctionProbabilityFlag(kFALSE), 
  fNSigma(0),
  fElectronFunction(0), fMuonFunction(0),
  fPionFunction(0), fKaonFunction(0), fProtonFunction(0),
  fDebugMode(kFALSE) {
  //Default constructor
  for(Int_t i = 0; i < 5; i++) fPartFrac[i] = 0.0;
  for(Int_t i = 0; i < 5; i++) {
    fdEdxMean[i] = 0.0;
    fdEdxSigma[i] = 0.0;
  }
}

//____________________________________________________________________//
AliProtonAnalysisBase::~AliProtonAnalysisBase() {
  //Default destructor
  if(fElectronFunction) delete fElectronFunction;
  if(fMuonFunction) delete fMuonFunction;
  if(fPionFunction) delete fPionFunction;
  if(fKaonFunction) delete fKaonFunction;
  if(fProtonFunction) delete fProtonFunction;
}

//____________________________________________________________________//
Double_t AliProtonAnalysisBase::GetParticleFraction(Int_t i, Double_t p) {
  //Return the a priori probs
  Double_t partFrac=0;
  if(fFunctionProbabilityFlag) {
    if(i == 0) partFrac = fElectronFunction->Eval(p);
    if(i == 1) partFrac = fMuonFunction->Eval(p);
    if(i == 2) partFrac = fPionFunction->Eval(p);
    if(i == 3) partFrac = fKaonFunction->Eval(p);
    if(i == 4) partFrac = fProtonFunction->Eval(p);
  }
  else partFrac = fPartFrac[i];

  return partFrac;
}

//____________________________________________________________________//
Bool_t AliProtonAnalysisBase::IsInPhaseSpace(AliESDtrack* const track) {
  // Checks if the track is outside the analyzed y-Pt phase space
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t eta = 0.0;

  if(kTPC) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0; eta = -10.0;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
      eta = tpcTrack->Eta();
    }
  }
  else {
    gPt = track->Pt();
    gPx = track->Px();
    gPy = track->Py();
    gPz = track->Pz();
    eta = track->Eta();
  }
  
  if((gPt < fMinY) || (gPt > fMaxY)) {
      if(fDebugMode)
	Printf("IsInPhaseSpace: Track rejected because it has a Pt value of %lf (accepted interval: %lf - %lf)",gPt,fMinY,fMaxY);
      return kFALSE;
  }
  if(fAnalysisEtaMode) {
    if((eta < fMinX) || (eta > fMaxX)) {
      if(fDebugMode)
	Printf("IsInPhaseSpace: Track rejected because it has an eta value of %lf (accepted interval: %lf - %lf)",eta,fMinX,fMaxX);
      return kFALSE;
    }
  }
  else {
    if((Rapidity(gPx,gPy,gPz) < fMinX) || (Rapidity(gPx,gPy,gPz) > fMaxX)) {
      if(fDebugMode)
	Printf("IsInPhaseSpace: Track rejected because it has a y value of %lf (accepted interval: %lf - %lf)",Rapidity(gPx,gPy,gPz),fMinX,fMaxX);
      return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________//
Bool_t AliProtonAnalysisBase::IsAccepted(AliESDEvent *esd,
					 const AliESDVertex *vertex, 
					 AliESDtrack* track) {
  // Checks if the track is excluded from the cuts
  Double_t gPt = 0.0, gPx = 0.0, gPy = 0.0, gPz = 0.0;
  Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};  //The impact parameters and their covariance.
  Double_t dca3D = 0.0;
  
  if((kTPC)&&(!kHybrid)) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.; dca3D = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
      tpcTrack->PropagateToDCA(vertex,
			       esd->GetMagneticField(),
			       100.,dca,cov);
    }
  }
  else if(kHybrid) {
     AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(!tpcTrack) {
      gPt = 0.0; gPx = 0.0; gPy = 0.0; gPz = 0.0;
      dca[0] = -100.; dca[1] = -100.;
      cov[0] = -100.; cov[1] = -100.; cov[2] = -100.;
    }
    else {
      gPt = tpcTrack->Pt();
      gPx = tpcTrack->Px();
      gPy = tpcTrack->Py();
      gPz = tpcTrack->Pz();
      tpcTrack->PropagateToDCA(vertex,
			       esd->GetMagneticField(),
			       100.,dca,cov);
    }
  }
  else{
    gPt = track->Pt();
    gPx = track->Px();
    gPy = track->Py();
    gPz = track->Pz();
    track->PropagateToDCA(vertex,
			  esd->GetMagneticField(),
			  100.,dca,cov);
  }
  dca3D = TMath::Sqrt(TMath::Power(dca[0],2) +
		      TMath::Power(dca[1],2));
     
  Int_t  fIdxInt[200];
  Int_t nClustersITS = track->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = track->GetTPCclusters(fIdxInt);

  Float_t chi2PerClusterITS = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = track->GetITSchi2()/Float_t(nClustersITS);
  Float_t chi2PerClusterTPC = -1;
  if (nClustersTPC!=0)
    chi2PerClusterTPC = track->GetTPCchi2()/Float_t(nClustersTPC);

  Double_t extCov[15];
  track->GetExternalCovariance(extCov);

  if(fPointOnITSLayer1Flag) {
    if(!track->HasPointOnITSLayer(0)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 1st ITS layer");
      return kFALSE;
    }
  }
  if(fPointOnITSLayer2Flag) {
    if(!track->HasPointOnITSLayer(1)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 2nd ITS layer");
      return kFALSE;
    }
  }
  if(fPointOnITSLayer3Flag) {
    if(!track->HasPointOnITSLayer(2)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 3rd ITS layer");
      return kFALSE;
    }
  }
  if(fPointOnITSLayer4Flag) {
    if(!track->HasPointOnITSLayer(3)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 4th ITS layer");
      return kFALSE;
    }
  }
  if(fPointOnITSLayer5Flag) {
    if(!track->HasPointOnITSLayer(4)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 5th ITS layer");
      return kFALSE;
    }
  }
  if(fPointOnITSLayer6Flag) {
    if(!track->HasPointOnITSLayer(5)) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it doesn't have a point on the 6th ITS layer");
      return kFALSE;
    }
  }
  if(fMinITSClustersFlag) {
    if(nClustersITS < fMinITSClusters) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has %d ITS points (min. requested: %d)",nClustersITS,fMinITSClusters);
      return kFALSE;
    }
  }
  if(fMaxChi2PerITSClusterFlag) {
    if(chi2PerClusterITS > fMaxChi2PerITSCluster) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a chi2 per ITS cluster %lf (max. requested: %lf)",chi2PerClusterITS,fMaxChi2PerITSCluster);
      return kFALSE; 
    }
  }
  if(fMinTPCClustersFlag) {
    if(nClustersTPC < fMinTPCClusters) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has %d TPC clusters (min. requested: %d)",nClustersTPC,fMinTPCClusters);
      return kFALSE;
    }
  }
  if(fMaxChi2PerTPCClusterFlag) {
    if(chi2PerClusterTPC > fMaxChi2PerTPCCluster) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a chi2 per TPC cluster %lf (max. requested: %lf)",chi2PerClusterTPC,fMaxChi2PerTPCCluster);
      return kFALSE; 
    }
  }
  if(fMaxCov11Flag) {
    if(extCov[0] > fMaxCov11) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a cov11 value of %lf (max. requested: %lf)",extCov[0],fMaxCov11);
      return kFALSE;
    }
  }
  if(fMaxCov22Flag) {
    if(extCov[2] > fMaxCov22) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a cov22 value of %lf (max. requested: %lf)",extCov[2],fMaxCov22);
      return kFALSE;
    }
  }
  if(fMaxCov33Flag) {
    if(extCov[5] > fMaxCov33) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a cov33 value of %lf (max. requested: %lf)",extCov[5],fMaxCov33);
      return kFALSE;
    }
  }
  if(fMaxCov44Flag) {
    if(extCov[9] > fMaxCov44) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a cov44 value of %lf (max. requested: %lf)",extCov[9],fMaxCov44);
      return kFALSE;
    }
  }
  if(fMaxCov55Flag) {
    if(extCov[14] > fMaxCov55) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a cov55 value of %lf (max. requested: %lf)",extCov[14],fMaxCov55);
      return kFALSE;
    }
  }
  if(fMaxSigmaToVertexFlag) {
    if(GetSigmaToVertex(track) > fMaxSigmaToVertex) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a %lf sigmas to vertex (max. requested: %lf)",GetSigmaToVertex(track),fMaxSigmaToVertex);
      return kFALSE;
    }
  }
  if(fMaxSigmaToVertexTPCFlag) {
    if(GetSigmaToVertex(track) > fMaxSigmaToVertexTPC) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a %lf sigmas to vertex TPC (max. requested: %lf)",GetSigmaToVertex(track),fMaxSigmaToVertexTPC);
      return kFALSE;
    }
  }
  if(fMaxDCAXYFlag) { 
    if(TMath::Abs(dca[0]) > fMaxDCAXY) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(xy) of %lf (max. requested: %lf)",TMath::Abs(dca[0]),fMaxDCAXY);
      return kFALSE;
    }
  }
  if(fMaxDCAXYTPCFlag) { 
    if(TMath::Abs(dca[0]) > fMaxDCAXYTPC) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(xy) (TPC) of %lf (max. requested: %lf)",TMath::Abs(dca[0]),fMaxDCAXYTPC);
      return kFALSE;
    }
  }
  if(fMaxDCAZFlag) { 
    if(TMath::Abs(dca[1]) > fMaxDCAZ) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(z) of %lf (max. requested: %lf)",TMath::Abs(dca[1]),fMaxDCAZ);
      return kFALSE;
    }
  }
  if(fMaxDCAZTPCFlag) { 
    if(TMath::Abs(dca[1]) > fMaxDCAZTPC) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(z) (TPC) of %lf (max. requested: %lf)",TMath::Abs(dca[1]),fMaxDCAZTPC);
      return kFALSE;
    }
  }
  if(fMaxDCA3DFlag) { 
    if(TMath::Abs(dca3D) > fMaxDCA3D) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(3D) of %lf (max. requested: %lf)",TMath::Abs(dca3D),fMaxDCA3D);
      return kFALSE;
    }
  }
  if(fMaxDCA3DTPCFlag) { 
    if(TMath::Abs(dca3D) > fMaxDCA3DTPC)  {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of dca(3D) (TPC) of %lf (max. requested: %lf)",TMath::Abs(dca3D),fMaxDCA3DTPC);
      return kFALSE;
    }
  }
  if(fMaxConstrainChi2Flag) {
    if(track->GetConstrainedChi2() > 0) 
      if(TMath::Log(track->GetConstrainedChi2()) > fMaxConstrainChi2)  {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has a value of the constrained chi2 to the vertex of %lf (max. requested: %lf)",TMath::Log(track->GetConstrainedChi2()),fMaxConstrainChi2);
      return kFALSE;
      }
  }
  if(fMinTPCdEdxPointsFlag) {
    if(track->GetTPCsignalN() < fMinTPCdEdxPoints) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has %d TPC points for the calculation of the energy loss (min. requested: %d)",track->GetTPCsignalN(),fMinTPCdEdxPoints);
      return kFALSE;
    }
  }
  if(fITSRefitFlag) {
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == 0) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has no ITS refit flag");
      return kFALSE;
    }
  }
  if(fTPCRefitFlag) {
    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == 0) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has no TPC refit flag");
      return kFALSE;
    }
  }
  if(fESDpidFlag) {
    if ((track->GetStatus() & AliESDtrack::kESDpid) == 0) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has no ESD pid flag");
      return kFALSE;
    }
  }
  if(fTPCpidFlag) {
    if ((track->GetStatus() & AliESDtrack::kTPCpid) == 0) {
      if(fDebugMode)
	Printf("IsAccepted: Track rejected because it has no TPC pid flag");
      return kFALSE;
    }
  }

  return kTRUE;
}

//____________________________________________________________________//
Float_t AliProtonAnalysisBase::GetSigmaToVertex(AliESDtrack* esdTrack) const {
  // Calculates the number of sigma to the vertex.
  
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  if((kTPC)&&(!kHybrid))
    esdTrack->GetImpactParametersTPC(b,bCov);
  else
    esdTrack->GetImpactParameters(b,bCov);
  
  if (bCov[0]<=0 || bCov[2]<=0) {
    //AliDebug(1, "Estimated b resolution lower or equal zero!");
    bCov[0]=0; bCov[2]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);
  
  if (bRes[0] == 0 || bRes[1] ==0) return -1;
  
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));
  
  if (TMath::Exp(-d * d / 2) < 1e-10) return 1000;
  
  d = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  
  return d;
}

//____________________________________________________________________//
Double_t AliProtonAnalysisBase::Rapidity(Double_t gPx, 
					 Double_t gPy, 
					 Double_t gPz) const {
  //returns the rapidity of the proton - to be removed
  Double_t fMass = 9.38270000000000048e-01;
  
  Double_t gP = TMath::Sqrt(TMath::Power(gPx,2) + 
                           TMath::Power(gPy,2) + 
			   TMath::Power(gPz,2));
  Double_t energy = TMath::Sqrt(gP*gP + fMass*fMass);
  Double_t y = -999;
  if(energy != gPz) 
    y = 0.5*TMath::Log((energy + gPz)/(energy - gPz));

  return y;
}

//________________________________________________________________________
const AliESDVertex* AliProtonAnalysisBase::GetVertex(AliESDEvent* esd, 
						     AnalysisMode mode,
						     Double_t gVxMax,
						     Double_t gVyMax,
						     Double_t gVzMax) {
  // Get the vertex from the ESD and returns it if the vertex is valid
  // Second argument decides which vertex is used (this selects
  // also the quality criteria that are applied)
  const AliESDVertex* vertex = 0;
  if(mode == kHybrid)
    vertex = esd->GetPrimaryVertexSPD();
  else if(mode == kTPC){
    Double_t kBz = esd->GetMagneticField();
    AliVertexerTracks vertexer(kBz);
    vertexer.SetTPCMode();
    AliESDVertex *vTPC = vertexer.FindPrimaryVertex(esd);
    esd->SetPrimaryVertexTPC(vTPC);
    for (Int_t i=0; i<esd->GetNumberOfTracks(); i++) {
      AliESDtrack *t = esd->GetTrack(i);
      t->RelateToVertexTPC(vTPC, kBz, kVeryBig);
    }
    delete vTPC;
    vertex = esd->GetPrimaryVertexTPC();
  }
  else if(mode == kGlobal)
    vertex = esd->GetPrimaryVertex();
  else
    Printf("GetVertex: ERROR: Invalid second argument %d", mode);
  
  if(!vertex) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because there is no valid vertex object");
    return 0;
  }
  
  // check Ncontributors
  if(vertex->GetNContributors() <= 0) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because the number of contributors for the vertex determination is <= 0");
    return 0;
  }
  
  // check resolution
  Double_t zRes = vertex->GetZRes();
  if(zRes == 0) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because the value of the vertex resolution in z is 0");
    return 0;
  }
  
  //check position
  if(TMath::Abs(vertex->GetXv()) > gVxMax) {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vx value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetXv()),gVxMax,gVxMax);
    return 0;
  }
  if(TMath::Abs(vertex->GetYv()) > gVyMax)  {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vy value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetYv()),gVyMax,gVyMax);
    return 0;
  }
  if(TMath::Abs(vertex->GetZv()) > gVzMax)  {
    if(fDebugMode)
      Printf("GetVertex: Event rejected because it has a Vz value of %lf cm (accepted interval: -%lf - %lf)",TMath::Abs(vertex->GetZv()),gVzMax,gVzMax);
    return 0;
  }
  
  return vertex;
}

//________________________________________________________________________
Bool_t AliProtonAnalysisBase::IsEventTriggered(const AliESDEvent *esd, 
					       TriggerMode trigger) {
  // check if the event was triggered
  ULong64_t triggerMask = esd->GetTriggerMask();
  
  // definitions from p-p.cfg
  ULong64_t spdFO = (1 << 14);
  ULong64_t v0left = (1 << 11);
  ULong64_t v0right = (1 << 12);

  switch (trigger) {
  case kMB1: {
    if (triggerMask & spdFO || ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  }
  case kMB2: {
    if (triggerMask & spdFO && ((triggerMask & v0left) || (triggerMask & v0right)))
      return kTRUE;
    break;
  }
  case kSPDFASTOR: {
    if (triggerMask & spdFO)
      return kTRUE;
    break;
  }
  }//switch

  return kFALSE;
}

//________________________________________________________________________
TCanvas *AliProtonAnalysisBase::GetListOfCuts() {
  // return the list of cuts and their cut values for reference
  TLatex l;
  l.SetTextAlign(12);
  l.SetTextSize(0.04);

  TString listOfCuts;

  TCanvas *c = new TCanvas("cListOfCuts","List of cuts",0,0,900,600);
  c->SetFillColor(10); c->SetHighLightColor(41);
  c->Divide(3,2);

  c->cd(1)->SetFillColor(10);
  l.DrawLatex(0.3,0.9,"Analysis details: List of cuts\n\n");

  listOfCuts = "Analysis level: "; listOfCuts += fProtonAnalysisLevel;
  l.DrawLatex(0.1,0.82,listOfCuts.Data());
  listOfCuts = "Analysis mode: ";
  if(fProtonAnalysisMode == kTPC) listOfCuts += "TPC standalone"; 
  if(fProtonAnalysisMode == kHybrid) listOfCuts += "Hybrid TPC"; 
  if(fProtonAnalysisMode == kGlobal) listOfCuts += "Global tracking"; 
  l.DrawLatex(0.1,0.74,listOfCuts.Data());
  listOfCuts = "Trigger mode: "; 
  if(fTriggerMode == kMB1) listOfCuts += "Minimum bias 1"; 
  if(fTriggerMode == kMB2) listOfCuts += "Minimum bias 2"; 
  if(fTriggerMode == kSPDFASTOR) listOfCuts += "FastOR (SPD)"; 
  l.DrawLatex(0.1,0.66,listOfCuts.Data());
  listOfCuts = "PID mode: "; 
  if(fProtonPIDMode == kBayesian) listOfCuts += "Bayesian PID";
  if(fProtonPIDMode == kRatio) listOfCuts += "Z = ln((dE/dx)_{exp.}/(dE/dx)_{theor.})"; 
  if(fProtonPIDMode == kSigma1) {
    listOfCuts += "N_{#sigma}(1) area: "; listOfCuts += fNSigma;
    listOfCuts += " #sigma";
  }
  if(fProtonPIDMode == kSigma2) {
    listOfCuts += "N_{#sigma}(2) area: "; listOfCuts += fNSigma;
    listOfCuts += " #sigma";
  }
  l.DrawLatex(0.1,0.58,listOfCuts.Data());
  listOfCuts = "Accepted vertex diamond: "; 
  l.DrawLatex(0.1,0.5,listOfCuts.Data());
  listOfCuts = "|V_{x}| < "; listOfCuts += fVxMax; listOfCuts += " [cm]";
  l.DrawLatex(0.6,0.5,listOfCuts.Data());
  listOfCuts = "|V_{y}| < "; listOfCuts += fVyMax; listOfCuts += " [cm]";
  l.DrawLatex(0.6,0.4,listOfCuts.Data());
  listOfCuts = "|V_{z}| < "; listOfCuts += fVzMax; listOfCuts += " [cm]";
  l.DrawLatex(0.6,0.3,listOfCuts.Data());
  listOfCuts = "Phase space: "; 
  l.DrawLatex(0.1,0.2,listOfCuts.Data());
  if(fAnalysisEtaMode) listOfCuts = "|#eta| < ";  
  else listOfCuts = "|y| < ";
  listOfCuts += TMath::Abs(fMinX); 
  l.DrawLatex(0.35,0.2,listOfCuts.Data());
  listOfCuts = "N_{bins} = ";
  listOfCuts += fNBinsX; listOfCuts += " (binning: ";
  listOfCuts += (fMaxX-fMinX)/fNBinsX; listOfCuts += ")";
  l.DrawLatex(0.35,0.15,listOfCuts.Data());
  listOfCuts = ""; listOfCuts += fMinY; listOfCuts += " < P_{T} < ";
  listOfCuts += fMaxY; listOfCuts += "GeV/c"; 
  l.DrawLatex(0.35,0.1,listOfCuts.Data());
  listOfCuts = "N_{bins} = ";
  listOfCuts += fNBinsY; listOfCuts += " (binning: ";
  listOfCuts += (fMaxY-fMinY)/fNBinsY; listOfCuts += ")";
  l.DrawLatex(0.35,0.05,listOfCuts.Data());

  c->cd(2)->SetFillColor(2);
  l.DrawLatex(0.3,0.9,"ITS related cuts");
  listOfCuts = "Request a cluster on SPD1: "; 
  listOfCuts += fPointOnITSLayer1Flag;
  l.DrawLatex(0.1,0.8,listOfCuts.Data());
  listOfCuts = "Request a cluster on SPD2: "; 
  listOfCuts += fPointOnITSLayer2Flag;
  l.DrawLatex(0.1,0.7,listOfCuts.Data());
  listOfCuts = "Request a cluster on SDD1: "; 
  listOfCuts += fPointOnITSLayer3Flag;
  l.DrawLatex(0.1,0.6,listOfCuts.Data());
  listOfCuts = "Request a cluster on SDD2: "; 
  listOfCuts += fPointOnITSLayer4Flag;
  l.DrawLatex(0.1,0.5,listOfCuts.Data());
  listOfCuts = "Request a cluster on SSD1: "; 
  listOfCuts += fPointOnITSLayer5Flag;
  l.DrawLatex(0.1,0.4,listOfCuts.Data());
  listOfCuts = "Request a cluster on SSD2: "; 
  listOfCuts += fPointOnITSLayer6Flag; 
  l.DrawLatex(0.1,0.3,listOfCuts.Data());  
  listOfCuts = "Minimum number of ITS clusters: ";
  if(fMinITSClustersFlag) listOfCuts += fMinITSClusters;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.2,listOfCuts.Data());
  listOfCuts = "Maximum #chi^{2} per ITS cluster: ";
  if(fMaxChi2PerITSClusterFlag) listOfCuts += fMaxChi2PerITSCluster; 
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.1,listOfCuts.Data());

  c->cd(3)->SetFillColor(3);
  l.DrawLatex(0.3,0.9,"TPC related cuts");
  listOfCuts = "Minimum number of TPC clusters: ";
  if(fMinTPCClustersFlag) listOfCuts += fMinTPCClusters;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.7,listOfCuts.Data());
  listOfCuts = "Maximum #chi^{2} per TPC cluster: ";
  if(fMaxChi2PerTPCClusterFlag) listOfCuts += fMaxChi2PerTPCCluster;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.5,listOfCuts.Data());
  listOfCuts = "Minimum number of TPC points for the dE/dx: ";
  if(fMinTPCdEdxPointsFlag) listOfCuts += fMinTPCdEdxPoints;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.3,listOfCuts.Data());

  c->cd(4)->SetFillColor(4);
  l.DrawLatex(0.3,0.9,"Tracking related cuts");
  listOfCuts = "Maximum cov11: ";
  if(fMaxCov11Flag) listOfCuts += fMaxCov11;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.75,listOfCuts.Data());
  listOfCuts = "Maximum cov22: ";
  if(fMaxCov22Flag) listOfCuts += fMaxCov22;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.6,listOfCuts.Data());
  listOfCuts = "Maximum cov33: ";
  if(fMaxCov33Flag) listOfCuts += fMaxCov33;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.45,listOfCuts.Data());
  listOfCuts = "Maximum cov44: ";
  if(fMaxCov44Flag) listOfCuts += fMaxCov44;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.3,listOfCuts.Data());
  listOfCuts = "Maximum cov55: ";
  if(fMaxCov55Flag) listOfCuts += fMaxCov55;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.15,listOfCuts.Data());

  c->cd(5)->SetFillColor(5);
  l.DrawLatex(0.3,0.9,"DCA related cuts");
  listOfCuts = "Maximum sigma to vertex: ";
  if(fMaxSigmaToVertexFlag) listOfCuts += fMaxSigmaToVertex;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.8,listOfCuts.Data());
  listOfCuts = "Maximum sigma to vertex (TPC): ";
  if(fMaxSigmaToVertexTPCFlag) listOfCuts += fMaxSigmaToVertexTPC;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.72,listOfCuts.Data());
  listOfCuts = "Maximum DCA in xy: ";
  if(fMaxDCAXYFlag) listOfCuts += fMaxDCAXY;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.64,listOfCuts.Data());
  listOfCuts = "Maximum DCA in z: ";
  if(fMaxDCAZFlag) listOfCuts += fMaxDCAZ;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.56,listOfCuts.Data());
  listOfCuts = "Maximum DCA in 3D: ";
  if(fMaxDCA3DFlag) listOfCuts += fMaxDCA3D;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.48,listOfCuts.Data());
  listOfCuts = "Maximum DCA in xy (TPC): ";
  if(fMaxDCAXYFlag) listOfCuts += fMaxDCAXYTPC;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.4,listOfCuts.Data());
  listOfCuts = "Maximum DCA in z (TPC): ";
  if(fMaxDCAZFlag) listOfCuts += fMaxDCAZTPC;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.32,listOfCuts.Data());
  listOfCuts = "Maximum DCA in 3D (TPC): ";
  if(fMaxDCA3DFlag) listOfCuts += fMaxDCA3DTPC;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.24,listOfCuts.Data());
  listOfCuts = "Maximum constrained #chi^{2} (vertex): ";
  if(fMaxConstrainChi2Flag) listOfCuts += fMaxConstrainChi2;
  else listOfCuts += "Not used";
  l.DrawLatex(0.1,0.16,listOfCuts.Data());

  c->cd(6)->SetFillColor(6);
  l.DrawLatex(0.3,0.9,"Tracking bits related cuts");
  listOfCuts = "Request the ITS refit bit: "; listOfCuts += fITSRefitFlag;
  l.DrawLatex(0.1,0.7,listOfCuts.Data());
  listOfCuts = "Request the TPC refit bit: "; listOfCuts += fTPCRefitFlag; 
  l.DrawLatex(0.1,0.5,listOfCuts.Data());
  listOfCuts = "Request the TPC pid bit: "; listOfCuts += fTPCpidFlag;
  l.DrawLatex(0.1,0.3,listOfCuts.Data());
  listOfCuts = "Request the ESD pid bit: "; listOfCuts += fESDpidFlag;
  l.DrawLatex(0.1,0.1,listOfCuts.Data());

  return c;
}

//________________________________________________________________________
Bool_t AliProtonAnalysisBase::IsProton(AliESDtrack *track) {
  //Function that checks if a track is a proton
  Double_t probability[5];
  Double_t gPt = 0.0, gP = 0.0, gEta = 0.0;
  Long64_t fParticleType = 0;
 
  //Bayesian approach for the PID
  if(fProtonPIDMode == kBayesian) {
    if((fProtonAnalysisMode == kTPC)||(fProtonAnalysisMode == kHybrid)) {
      AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
      if(tpcTrack) {
	gPt = tpcTrack->Pt();
	gP = tpcTrack->P();
	track->GetTPCpid(probability);
      }
    }//TPC standalone or Hybrid TPC
    else if(fProtonAnalysisMode == kGlobal) {
      gPt = track->Pt();
      gP = track->P();
      track->GetESDpid(probability);
    }//Global tracking    
    
    Double_t rcc = 0.0;
    for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
      rcc += probability[i]*GetParticleFraction(i,gP);
    if(rcc != 0.0) {
      Double_t w[5];
      for(Int_t i = 0; i < AliPID::kSPECIES; i++) 
	w[i] = probability[i]*GetParticleFraction(i,gP)/rcc;
      fParticleType = TMath::LocMax(AliPID::kSPECIES,w);
    }
    if(fParticleType == 4)
      return kTRUE;
  }
  //Ratio of the measured over the theoretical dE/dx a la STAR
  else if(fProtonPIDMode == kRatio) {
    Printf("The kRatio mode is not implemented yet!!!");
    return kFALSE;
  }
  //Definition of an N-sigma area around the dE/dx vs P band
  else if(fProtonPIDMode == kSigma1) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = tpcTrack->P();
      gEta = tpcTrack->Eta();
    }
    //We start the P slices at >0.3GeV/c with a bining of 50MeV/c ==> Int_t(0.3001/0.05) = 6
    Int_t nbinP = Int_t(gP/0.05) - 6;
    Double_t tpcSignal = track->GetTPCsignal();
    Double_t dEdxMean = fdEdxMean[nbinP];
    Double_t dEdxSigma = fdEdxSigma[nbinP];
    if((tpcSignal <= dEdxMean + fNSigma*dEdxSigma)&&
       (tpcSignal <= dEdxMean + fNSigma*dEdxSigma))
      return kTRUE;
  }//kSigma1 PID method
  //Another definition of an N-sigma area around the dE/dx vs P band
  else if(fProtonPIDMode == kSigma2) {
    AliExternalTrackParam *tpcTrack = (AliExternalTrackParam *)track->GetTPCInnerParam();
    if(tpcTrack) {
      gPt = tpcTrack->Pt();
      gP = tpcTrack->P();
      gEta = tpcTrack->Eta();
    }
    //We start the P slices at >0.3GeV/c with a bining of 50MeV/c ==> Int_t(0.3001/0.05) = 6
    Int_t nbinP = Int_t(gP/0.05) - 6;
    Double_t tpcSignal = track->GetTPCsignal();
    Double_t dEdxTheory = Bethe(gP/9.38270000000000048e-01);
    Double_t dEdxSigma = fdEdxSigma[nbinP];
    Double_t nsigma = TMath::Abs(tpcSignal - dEdxTheory)/(tpcSignal*(dEdxSigma/TMath::Sqrt(track->GetTPCsignalN())));
    if(nsigma <= fNSigma) 
      return kTRUE;
  }

  return kFALSE;
}

//________________________________________________________________________
void AliProtonAnalysisBase::SetdEdxBandInfo(const char *filename) {
  // This function is used in case the kSigma1 or kSigma2 PID mode is selected
  // It takes as an argument the name of the ascii file (for the time being) 
  // that is generated as a prior process.
  // This ascii file has three columns: The min. P value (bins of 50MeV/c) 
  // the mean and the sigma of the dE/dx distributions for protons coming 
  // from a gaussian fit.
  ifstream in;
  in.open(filename);

  Double_t gPtMin = 0.0;
  Int_t iCounter = 0;
  while(in.good()) {
    in >> gPtMin >> fdEdxMean[iCounter] >> fdEdxSigma[iCounter];
    if(fDebugMode)
      Printf("Momentum bin: %d - Min momentum: %lf - mean(dE/dx): %lf - sigma(dE/dx): %lf",iCounter+1,gPtMin,fdEdxMean[iCounter],fdEdxSigma[iCounter]);
    iCounter += 1;
  }
}

//________________________________________________________________________
Double_t AliProtonAnalysisBase::Bethe(Double_t bg) const {
  // This is the Bethe-Bloch function normalised to 1 at the minimum
  // We renormalize it based on the MC information
  // WARNING: To be revised soon!!!
  // This is just a temporary fix
  Double_t normalization = 49.2;
  Double_t bg2=bg*bg;
  Double_t beta2 = bg2/(1.+ bg2);
  Double_t bb = 8.62702e-2*(9.14550 - beta2 - TMath::Log(3.51000e-5 + 1./bg2))/beta2;
  //
  const Float_t kmeanCorrection =0.1;
  Double_t meanCorrection =(1+(bb-1)*kmeanCorrection);
  bb *= meanCorrection;

  return normalization*bb; 
}


