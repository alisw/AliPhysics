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

// The class AliCFTrackIsPrimaryCut is designed to select reconstructed tracks
// with a small impact parameter and tracks which are (not) daughters of kink
// decays and to provide corresponding QA histograms.
// This class inherits from the Analysis' Framework abstract base class
// AliAnalysisCuts and is a part of the Correction Framework.
// This class acts on single, reconstructed tracks, it is applicable on
// ESD and AOD data.
// It mainly consists of a IsSelected function that returns a boolean.
// This function checks whether the considered track passes a set of cuts:
// - min. and max. distance to main vertex in transverse plane (xy)
// - min. and max. longitudinal distance to main vertex (z)
// - min. and max. distance to main vertex as ellpise in xy - z plane
// - all above cuts on absolute values or in units of sigma (resolution)
// - min. and max. distance to main vertex in units of sigma (resolution)
// - max. transverse (xy) and longitudinal (z) impact parameter resolution
// - require that the dca calculation doesn't fail
// - accept or not accept daughter tracks of kink decays
//
// By default, the distance to 'vertex calculated from tracks' is used.
// Optionally the SPD (tracklet based) or TPC (TPC only tracks based) vertex
// can be used.
// Note: the distance to the TPC-vertex is already stored in the ESD,
// the distance to the SPD-vertex has to be re-calculated by propagating each
// track while executing this cut.
//
// The cut values for these cuts are set with the corresponding set functions.
// All cut classes provided by the correction framework are supposed to be
// added in the Analysis Framwork's class AliAnalysisFilter and applied by
// the filter via a loop.
//
// author: I. Kraus (Ingrid.Kraus@cern.ch)
// idea taken form
// AliESDtrackCuts writte by Jan Fiete Grosse-Oetringhaus and
// AliRsnDaughterCut class written by A. Pulvirenti.

#include <TCanvas.h>
#include <TDirectory.h>
#include <TH2.h>
#include <TBits.h>

#include <AliESDtrack.h>
#include <AliAODTrack.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliLog.h>
#include "AliCFTrackIsPrimaryCuts.h"

ClassImp(AliCFTrackIsPrimaryCuts)

//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts() :
  AliCFCutBase(),
  fEvt(0x0),
  fUseSPDvertex(0),
  fUseTPCvertex(0),
  fMinDCAToVertexXY(0),
  fMinDCAToVertexZ(0),
  fMaxDCAToVertexXY(0),
  fMaxDCAToVertexZ(0),
  fDCAToVertex2D(0),
  fAbsDCAToVertex(0),
  fNSigmaToVertexMin(0),
  fNSigmaToVertexMax(0),
  fSigmaDCAxy(0),
  fSigmaDCAz(0),
  fRequireSigmaToVertex(0),
  fAODType(AliAODTrack::kUndef),
  fAcceptKinkDaughters(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsNSigma(0),
  fhNBinsRequireSigma(0),
  fhNBinsAcceptKink(0),
  fhNBinsDcaXY(0),
  fhNBinsDcaZ(0),
  fhNBinsDcaXYnorm(0),
  fhNBinsDcaZnorm(0),
  fhNBinsSigmaDcaXY(0),
  fhNBinsSigmaDcaZ(0),
  fhBinLimNSigma(0x0),
  fhBinLimRequireSigma(0x0),
  fhBinLimAcceptKink(0x0),
  fhBinLimDcaXY(0x0),
  fhBinLimDcaZ(0x0),
  fhBinLimDcaXYnorm(0x0),
  fhBinLimDcaZnorm(0x0),
  fhBinLimSigmaDcaXY(0x0),
  fhBinLimSigmaDcaZ(0x0)
{
  //
  // Default constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts(Char_t* name, Char_t* title) :
  AliCFCutBase(name,title),
  fEvt(0x0),
  fUseSPDvertex(0),
  fUseTPCvertex(0),
  fMinDCAToVertexXY(0),
  fMinDCAToVertexZ(0),
  fMaxDCAToVertexXY(0),
  fMaxDCAToVertexZ(0),
  fDCAToVertex2D(0),
  fAbsDCAToVertex(0),
  fNSigmaToVertexMin(0),
  fNSigmaToVertexMax(0),
  fSigmaDCAxy(0),
  fSigmaDCAz(0),
  fRequireSigmaToVertex(0),
  fAODType(AliAODTrack::kUndef),
  fAcceptKinkDaughters(0),
  fhCutStatistics(0),
  fhCutCorrelation(0),
  fBitmap(0x0),
  fhNBinsNSigma(0),
  fhNBinsRequireSigma(0),
  fhNBinsAcceptKink(0),
  fhNBinsDcaXY(0),
  fhNBinsDcaZ(0),
  fhNBinsDcaXYnorm(0),
  fhNBinsDcaZnorm(0),
  fhNBinsSigmaDcaXY(0),
  fhNBinsSigmaDcaZ(0),
  fhBinLimNSigma(0x0),
  fhBinLimRequireSigma(0x0),
  fhBinLimAcceptKink(0x0),
  fhBinLimDcaXY(0x0),
  fhBinLimDcaZ(0x0),
  fhBinLimDcaXYnorm(0x0),
  fhBinLimDcaZnorm(0x0),
  fhBinLimSigmaDcaXY(0x0),
  fhBinLimSigmaDcaZ(0x0)
{
  //
  // Constructor
  //
  Initialise();
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::AliCFTrackIsPrimaryCuts(const AliCFTrackIsPrimaryCuts& c) :
  AliCFCutBase(c),
  fEvt(c.fEvt),
  fUseSPDvertex(c.fUseSPDvertex),
  fUseTPCvertex(c.fUseTPCvertex),
  fMinDCAToVertexXY(c.fMinDCAToVertexXY),
  fMinDCAToVertexZ(c.fMinDCAToVertexZ),
  fMaxDCAToVertexXY(c.fMaxDCAToVertexXY),
  fMaxDCAToVertexZ(c.fMaxDCAToVertexZ),
  fDCAToVertex2D(c.fDCAToVertex2D),
  fAbsDCAToVertex(c.fAbsDCAToVertex),
  fNSigmaToVertexMin(c.fNSigmaToVertexMin),
  fNSigmaToVertexMax(c.fNSigmaToVertexMax),
  fSigmaDCAxy(c.fSigmaDCAxy),
  fSigmaDCAz(c.fSigmaDCAz),
  fRequireSigmaToVertex(c.fRequireSigmaToVertex),
  fAODType(c.fAODType),
  fAcceptKinkDaughters(c.fAcceptKinkDaughters),
  fhCutStatistics(c.fhCutStatistics),
  fhCutCorrelation(c.fhCutCorrelation),
  fBitmap(c.fBitmap),
  fhNBinsNSigma(c.fhNBinsNSigma),
  fhNBinsRequireSigma(c.fhNBinsRequireSigma),
  fhNBinsAcceptKink(c.fhNBinsAcceptKink),
  fhNBinsDcaXY(c.fhNBinsDcaXY),
  fhNBinsDcaZ(c.fhNBinsDcaZ),
  fhNBinsDcaXYnorm(c.fhNBinsDcaXYnorm),
  fhNBinsDcaZnorm(c.fhNBinsDcaZnorm),
  fhNBinsSigmaDcaXY(c.fhNBinsSigmaDcaXY),
  fhNBinsSigmaDcaZ(c.fhNBinsSigmaDcaZ),
  fhBinLimNSigma(c.fhBinLimNSigma),
  fhBinLimRequireSigma(c.fhBinLimRequireSigma),
  fhBinLimAcceptKink(c.fhBinLimAcceptKink),
  fhBinLimDcaXY(c.fhBinLimDcaXY),
  fhBinLimDcaZ(c.fhBinLimDcaZ),
  fhBinLimDcaXYnorm(c.fhBinLimDcaXYnorm),
  fhBinLimDcaZnorm(c.fhBinLimDcaZnorm),
  fhBinLimSigmaDcaXY(c.fhBinLimSigmaDcaXY),
  fhBinLimSigmaDcaZ(c.fhBinLimSigmaDcaZ)
{
  //
  // copy constructor
  //
  ((AliCFTrackIsPrimaryCuts &) c).Copy(*this);
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts& AliCFTrackIsPrimaryCuts::operator=(const AliCFTrackIsPrimaryCuts& c)
{
  //
  // Assignment operator
  //
  if (this != &c) {
    AliCFCutBase::operator=(c) ;
    fEvt = c.fEvt;
    fUseSPDvertex = c.fUseSPDvertex;
    fUseTPCvertex = c.fUseTPCvertex;
    fMinDCAToVertexXY = c.fMinDCAToVertexXY;
    fMinDCAToVertexZ = c.fMinDCAToVertexZ;
    fMaxDCAToVertexXY = c.fMaxDCAToVertexXY;
    fMaxDCAToVertexZ = c.fMaxDCAToVertexZ;
    fDCAToVertex2D = c.fDCAToVertex2D;
    fAbsDCAToVertex = c.fAbsDCAToVertex;
    fNSigmaToVertexMin = c.fNSigmaToVertexMin ;
    fNSigmaToVertexMax = c.fNSigmaToVertexMax ;
    fSigmaDCAxy = c.fSigmaDCAxy ;
    fSigmaDCAz = c.fSigmaDCAz ;
    fRequireSigmaToVertex = c.fRequireSigmaToVertex ;
    fAODType = c.fAODType ;
    fAcceptKinkDaughters = c.fAcceptKinkDaughters ;
    fhCutStatistics = c.fhCutStatistics ;
    fhCutCorrelation = c.fhCutCorrelation ;
    fBitmap =  c.fBitmap;
    fhNBinsNSigma = c.fhNBinsNSigma;
    fhNBinsRequireSigma = c.fhNBinsRequireSigma;
    fhNBinsAcceptKink = c.fhNBinsAcceptKink;
    fhNBinsDcaXY = c.fhNBinsDcaXY;
    fhNBinsDcaZ = c.fhNBinsDcaZ;
    fhNBinsDcaXYnorm = c.fhNBinsDcaXYnorm;
    fhNBinsDcaZnorm = c.fhNBinsDcaZnorm;
    fhNBinsSigmaDcaXY = c.fhNBinsSigmaDcaXY;
    fhNBinsSigmaDcaZ = c.fhNBinsSigmaDcaZ;
    fhBinLimNSigma = c.fhBinLimNSigma;
    fhBinLimRequireSigma = c.fhBinLimRequireSigma;
    fhBinLimAcceptKink = c.fhBinLimAcceptKink;
    fhBinLimDcaXY = c.fhBinLimDcaXY;
    fhBinLimDcaZ = c.fhBinLimDcaZ;
    fhBinLimDcaXYnorm = c.fhBinLimDcaXYnorm;
    fhBinLimDcaZnorm = c.fhBinLimDcaZnorm;
    fhBinLimSigmaDcaXY = c.fhBinLimSigmaDcaXY;
    fhBinLimSigmaDcaZ = c.fhBinLimSigmaDcaZ;

    for (Int_t j=0; j<6; j++) fDCA[j] = c.fDCA[j];
    for (Int_t j=0; j<c.kNStepQA; j++){
      if(c.fhDcaXYvsDcaZ[j]) fhDcaXYvsDcaZ[j] = (TH2F*)c.fhDcaXYvsDcaZ[j]->Clone();
      for (Int_t i=0; i<c.kNHist; i++){
	if(c.fhQA[i][j]) fhQA[i][j] = (TH1F*)c.fhQA[i][j]->Clone();
      }
    }
    ((AliCFTrackIsPrimaryCuts &) c).Copy(*this);
 }
  return *this;
}
//__________________________________________________________________________________
AliCFTrackIsPrimaryCuts::~AliCFTrackIsPrimaryCuts()
{
  //
  // destructor
  //
  if (fhCutStatistics)			delete fhCutStatistics;
  if (fhCutCorrelation)			delete fhCutCorrelation;

  for (Int_t j=0; j<kNStepQA; j++){
    if(fhDcaXYvsDcaZ[j])	delete fhDcaXYvsDcaZ[j];
    for (Int_t i=0; i<kNHist; i++)
      if(fhQA[i][j]) 		delete fhQA[i][j];
  }
  if(fEvt) 			delete fEvt;
  if(fBitmap) 			delete fBitmap;
  if(fhBinLimNSigma) 		delete fhBinLimNSigma;
  if(fhBinLimRequireSigma) 	delete fhBinLimRequireSigma;
  if(fhBinLimAcceptKink) 	delete fhBinLimAcceptKink;
  if(fhBinLimDcaXY) 		delete fhBinLimDcaXY;
  if(fhBinLimDcaZ) 		delete fhBinLimDcaZ;
  if(fhBinLimDcaXYnorm) 	delete fhBinLimDcaXYnorm;
  if(fhBinLimDcaZnorm) 		delete fhBinLimDcaZnorm;
  if(fhBinLimSigmaDcaXY) 	delete fhBinLimSigmaDcaXY;
  if(fhBinLimSigmaDcaZ) 	delete fhBinLimSigmaDcaZ;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::Initialise()
{
  //
  // sets everything to zero
  //
  fUseSPDvertex = 0;
  fUseTPCvertex = 0;
  fMinDCAToVertexXY = 0;
  fMinDCAToVertexZ = 0;
  fMaxDCAToVertexXY = 0;
  fMaxDCAToVertexZ = 0;
  fDCAToVertex2D = 0;
  fAbsDCAToVertex = 0;
  fNSigmaToVertexMin = 0;
  fNSigmaToVertexMax = 0;
  fSigmaDCAxy = 0;
  fSigmaDCAz = 0;
  fRequireSigmaToVertex = 0;
  fAcceptKinkDaughters = 0;
  fAODType = AliAODTrack::kUndef;

  SetMinDCAToVertexXY();
  SetMinDCAToVertexZ();
  SetMaxDCAToVertexXY();
  SetMaxDCAToVertexZ();
  SetDCAToVertex2D();
  SetAbsDCAToVertex();
  SetMinNSigmaToVertex();
  SetMaxNSigmaToVertex();
  SetMaxSigmaDCAxy();
  SetMaxSigmaDCAz();
  SetRequireSigmaToVertex();
  SetAcceptKinkDaughters();
  SetAODType();

  for (Int_t j=0; j<6; j++) fDCA[j] = 0.;
  for (Int_t j=0; j<kNStepQA; j++)  {
    fhDcaXYvsDcaZ[j] = 0x0;
    for (Int_t i=0; i<kNHist; i++)
      fhQA[i][j] = 0x0;
  }
  fhCutStatistics = 0;
  fhCutCorrelation = 0;
  fBitmap=new TBits(0);

  //set default bining for QA histograms
  SetHistogramBins(kCutNSigmaToVertex,100,0.,10.);
  SetHistogramBins(kCutRequireSigmaToVertex,5,-0.75,1.75);
  SetHistogramBins(kCutAcceptKinkDaughters,5,-0.75,1.75);
  SetHistogramBins(kDcaXY,500,-10.,10.);
  SetHistogramBins(kDcaZ,500,-10.,10.);
  SetHistogramBins(kDcaXYnorm,500,-10.,10.);
  SetHistogramBins(kDcaZnorm,500,-10.,10.);
  SetHistogramBins(kSigmaDcaXY,500,-0.1,0.9);
  SetHistogramBins(kSigmaDcaZ,500,-0.1,0.9);
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::Copy(TObject &c) const
{
  //
  // Copy function
  //
  AliCFTrackIsPrimaryCuts& target = (AliCFTrackIsPrimaryCuts &) c;

  target.Initialise();

  if (fhCutStatistics)  target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  for (Int_t j=0; j<6; j++) target.fDCA[j] = fDCA[j];
  for (Int_t j=0; j<kNStepQA; j++){
    if(fhDcaXYvsDcaZ[j]) target.fhDcaXYvsDcaZ[j] = (TH2F*)fhDcaXYvsDcaZ[j]->Clone();
    for (Int_t i=0; i<kNHist; i++)
      if(fhQA[i][j]) target.fhQA[i][j] = (TH1F*)fhQA[i][j]->Clone();
  }
  TNamed::Copy(c);
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SetRecEventInfo(const TObject* evt) {
  //
  // Sets pointer to event information (AliESDEvent or AliAODEvent)
  //
  if (!evt) {
    AliError("Pointer to AliVEvent !");
    return;
  }
  TString className(evt->ClassName());
  if (! (className.CompareTo("AliESDEvent")==0 || className.CompareTo("AliAODEvent")==0)) {
    AliError("argument must point to an AliESDEvent or AliAODEvent !");
    return ;
  }
  fEvt = (AliVEvent*) evt;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::UseSPDvertex(Bool_t b) {
  fUseSPDvertex = b;
  if(fUseTPCvertex && fUseSPDvertex) {
	fUseSPDvertex = kFALSE;
	AliError("SPD and TPC vertex chosen. TPC vertex is preferred.");
  }
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::UseTPCvertex(Bool_t b) {
  fUseTPCvertex = b;
  if(fUseTPCvertex) fUseSPDvertex = kFALSE;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::GetDCA(AliESDtrack* esdTrack)
{
  if (!esdTrack) return;

  Float_t b[2] = {0.,0.};
  Float_t bCov[3] = {0.,0.,0.};
  if(!fUseSPDvertex && !fUseTPCvertex)	esdTrack->GetImpactParameters(b,bCov);
  if( fUseTPCvertex)			esdTrack->GetImpactParametersTPC(b,bCov);

  if( fUseSPDvertex) {
	if (!fEvt) return;
	AliESDEvent * evt = 0x0 ; 
	evt = dynamic_cast<AliESDEvent*>(fEvt);
	if (!evt) {
	  AliError("event not found");
	  return;
	}
	const AliESDVertex *vtx = evt->GetVertex();
	const Double_t Bz = evt->GetMagneticField();
	AliExternalTrackParam *cParam = 0;
	Bool_t success = esdTrack->RelateToVertex(vtx, Bz, kVeryBig, cParam);
	if (success) esdTrack->GetImpactParameters(b,bCov);
  }

  if (bCov[0]<=0 || bCov[2]<=0) {
      bCov[0]=0; bCov[2]=0;
  }
  fDCA[0] = b[0]; // impact parameter xy
  fDCA[1] = b[1]; // impact parameter z
  fDCA[2] = TMath::Sqrt(bCov[0]); // resolution xy
  fDCA[3] = TMath::Sqrt(bCov[2]); // resolution z

  if (!fAbsDCAToVertex) {
	if (fDCA[2] > 0) fDCA[0] = fDCA[0]/fDCA[2]; // normalised impact parameter xy
	if (fDCA[3] > 0) fDCA[1] = fDCA[1]/fDCA[3]; // normalised impact parameter z
  }

  // get n_sigma
  if(!fUseSPDvertex && !fUseTPCvertex)
	fDCA[5] = AliESDtrackCuts::GetSigmaToVertex(esdTrack);

  if(fUseTPCvertex) {
	fDCA[5] = -1;
	if (fDCA[2]==0 || fDCA[3]==0)
	  return;
	fDCA[5] = 1000.;
	Float_t d = TMath::Sqrt(TMath::Power(b[0]/fDCA[2],2) + TMath::Power(b[1]/fDCA[3],2));
	if (TMath::Exp(-d * d / 2) < 1e-15)
	  return;
	fDCA[5] = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);
  }
  return;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::GetDCA(AliAODTrack* aodTrack)
{
  if (!aodTrack) return;

  Double_t p[3] = {0.,0.,0.};
  Double_t cov[21] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  aodTrack->XYZAtDCA(p); // position at DCA
  aodTrack->GetCovarianceXYZPxPyPz(cov);


  fDCA[5] = -1; // n_sigma
  if (p[0]==-999. || p[1]==-999. ||  p[2]==-999.) {
    AliError("dca info not available !");
    fDCA[0] = -999.; // impact parameter xy
    fDCA[1] = -999.; // impact parameter z
    return;
  }

  AliAODEvent * evt = 0x0;
  evt = dynamic_cast<AliAODEvent*>(fEvt);
  if (!evt) return;

  // primary vertex is the "best": tracks, SPD or TPC vertex
  AliAODVertex * primaryVertex = evt->GetVertex(0);
  // dca = track postion - primary vertex position
  p[0] = p[0] - primaryVertex->GetX();
  p[1] = p[1] - primaryVertex->GetY();
  p[2] = p[2] - primaryVertex->GetZ();
  // calculate dca in transverse plane
  Float_t b[2] = {0.,0.};
  b[0] = TMath::Sqrt(p[0]*p[0]+p[1]*p[1]);
  b[1] = p[2];
  // resolution
  Double_t bCov[3] = {0.,0.,0.};
  // how to calculate the resoultion in the transverse plane ?
  bCov[0] = 0.; // to do: calculate cov in transverse plane
  bCov[2] = 0.; // from cov in x and y, need to know correlation

  if (bCov[0]<=0 || bCov[2]<=0) {
      bCov[0]=0; bCov[2]=0;
  }
  fDCA[0] = b[0]; // impact parameter xy
  fDCA[1] = b[1]; // impact parameter z
  fDCA[2] = TMath::Sqrt(bCov[0]); // resolution xy
  fDCA[3] = TMath::Sqrt(bCov[2]); // resolution z

  if (!fAbsDCAToVertex) {
	AliError("resolution of impact parameter not available, use absolute dca cut instead !");
	if (fDCA[2] > 0) fDCA[0] = fDCA[0]/fDCA[2]; // normalised impact parameter xy
	if (fDCA[3] > 0) fDCA[1] = fDCA[1]/fDCA[3]; // normalised impact parameter z
  }

  // get n_sigma
  if (fDCA[2]==0 || fDCA[3]==0)
	return;
  fDCA[5] = 1000.;
  Float_t d = TMath::Sqrt(TMath::Power(b[0]/fDCA[2],2) + TMath::Power(b[1]/fDCA[3],2));
  if (TMath::Exp(-d * d / 2) < 1e-15)
	return;
  fDCA[5] = TMath::ErfInverse(1 - TMath::Exp(-d * d / 2)) * TMath::Sqrt(2);

  return;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SelectionBitMap(TObject* obj)
{
  //
  // test if the track passes the single cuts
  // and store the information in a bitmap
  //

  // bitmap stores the decision of each single cut
  for(Int_t i=0; i<kNCuts; i++)fBitmap->SetBitNumber(i,kFALSE);

  // check TObject and cast into ESDtrack
  if (!obj) return;
  if (!obj->InheritsFrom("AliVParticle")) {
    AliError("object must derived from AliVParticle !");
    return;
  }
  
  AliESDtrack * esdTrack = dynamic_cast<AliESDtrack*>(obj);
  AliAODTrack * aodTrack = dynamic_cast<AliAODTrack*>(obj);

  if (!(esdTrack || aodTrack)) {
    AliError("object must be an ESDtrack or an AODtrack !");
    return;
  }

  Bool_t isESDTrack = kFALSE;
  Bool_t isAODTrack = kFALSE;

  if (esdTrack) isESDTrack = strcmp(obj->ClassName(),"AliESDtrack") == 0 ? kTRUE : kFALSE ;
  if (aodTrack) isAODTrack = strcmp(obj->ClassName(),"AliAODTrack") == 0 ? kTRUE : kFALSE ;

  // get the track to vertex parameter for ESD track
  if (isESDTrack) GetDCA(esdTrack);
  if (isAODTrack) GetDCA(aodTrack);

  // check whether dca info is filled
  Bool_t dcaInfo = 0;
  if (fDCA[0]>-990. && fDCA[1]>-990.) dcaInfo = 1;

  Float_t bxy = 0, bz = 0;
  bxy = TMath::Abs(fDCA[0]);
  bz  = TMath::Abs(fDCA[1]);

  Float_t b2Dmin = 0, b2Dmax = 0;
  if (fMinDCAToVertexXY>0 && fMinDCAToVertexZ>0)
    b2Dmin = fDCA[0]*fDCA[0]/fMinDCAToVertexXY/fMinDCAToVertexXY + fDCA[1]*fDCA[1]/fMinDCAToVertexZ/fMinDCAToVertexZ;
  if (fMaxDCAToVertexXY>0 && fMaxDCAToVertexZ>0) 
    b2Dmax = fDCA[0]*fDCA[0]/fMaxDCAToVertexXY/fMaxDCAToVertexXY + fDCA[1]*fDCA[1]/fMaxDCAToVertexZ/fMaxDCAToVertexZ;


  // fill the bitmap
  Int_t iCutBit = 0;

  if (!dcaInfo || fDCAToVertex2D || (!fDCAToVertex2D && bxy >= fMinDCAToVertexXY && bxy <= fMaxDCAToVertexXY))
	fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || fDCAToVertex2D || (!fDCAToVertex2D && bz  >= fMinDCAToVertexZ && bz  <= fMaxDCAToVertexZ))
	fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || !fDCAToVertex2D || (fDCAToVertex2D && TMath::Sqrt(b2Dmin) > 1  && TMath::Sqrt(b2Dmax) < 1))
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || (fDCA[5] >= fNSigmaToVertexMin && fDCA[5] <= fNSigmaToVertexMax))
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || fDCA[2] < fSigmaDCAxy)
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || fDCA[3] < fSigmaDCAz)
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || !fRequireSigmaToVertex || (fDCA[5]>=0 && fRequireSigmaToVertex))
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (!dcaInfo || fAcceptKinkDaughters || (!fAcceptKinkDaughters && esdTrack->GetKinkIndex(0)<=0))
      fBitmap->SetBitNumber(iCutBit,kTRUE);
  iCutBit++;

  if (isAODTrack) {
    if (fAODType==AliAODTrack::kUndef || fAODType == aodTrack->GetType()) {
      fBitmap->SetBitNumber(iCutBit,kTRUE);
    }
  }
  else fBitmap->SetBitNumber(iCutBit,kTRUE);

  return;
}
//__________________________________________________________________________________
Bool_t AliCFTrackIsPrimaryCuts::IsSelected(TObject* obj) {
  //
  // loops over decisions of single cuts and returns if the track is accepted
  //
  SelectionBitMap(obj);

  if (fIsQAOn) FillHistograms(obj,0);
  Bool_t isSelected = kTRUE;

  for (UInt_t icut=0; icut<fBitmap->GetNbits();icut++) {
    if(!fBitmap->TestBitNumber(icut)) {
	isSelected = kFALSE;
	break;
    }
  }

  if (!isSelected) return kFALSE ;
  if (fIsQAOn) FillHistograms(obj,1);
  return kTRUE;
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t *bins)
{
  //
  // variable bin size
  //
  if(!fIsQAOn) return;

  switch(index){
  case kCutNSigmaToVertex:
    fhNBinsNSigma=nbins+1;
    fhBinLimNSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNSigma[i]=bins[i];
    break;

  case kCutRequireSigmaToVertex:
    fhNBinsRequireSigma=nbins+1;
    fhBinLimRequireSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRequireSigma[i]=bins[i];
    break;

  case kCutAcceptKinkDaughters:
    fhNBinsAcceptKink=nbins+1;
    fhBinLimAcceptKink=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimAcceptKink[i]=bins[i];
    break;

  case kDcaXY:
    fhNBinsDcaXY=nbins+1;
    fhBinLimDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXY[i]=bins[i];
    break;

  case kDcaZ:
    fhNBinsDcaZ=nbins+1;
    fhBinLimDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZ[i]=bins[i];
    break;

  case kDcaXYnorm:
    fhNBinsDcaXYnorm=nbins+1;
    fhBinLimDcaXYnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXYnorm[i]=bins[i];
    break;

  case kDcaZnorm:
    fhNBinsDcaZnorm=nbins+1;
    fhBinLimDcaZnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZnorm[i]=bins[i];
    break;

  case kSigmaDcaXY:
    fhNBinsSigmaDcaXY=nbins+1;
    fhBinLimSigmaDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimSigmaDcaXY[i]=bins[i];
    break;

  case kSigmaDcaZ:
    fhNBinsSigmaDcaZ=nbins+1;
    fhBinLimSigmaDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimSigmaDcaZ[i]=bins[i];
    break;
  }
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SetHistogramBins(Int_t index, Int_t nbins, Double_t xmin, Double_t xmax)
{
  //
  // fixed bin size
  //
  switch(index){
  case kCutNSigmaToVertex:
    fhNBinsNSigma=nbins+1;
    fhBinLimNSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimNSigma[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutRequireSigmaToVertex:
    fhNBinsRequireSigma=nbins+1;
    fhBinLimRequireSigma=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimRequireSigma[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kCutAcceptKinkDaughters:
    fhNBinsAcceptKink=nbins+1;
    fhBinLimAcceptKink=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimAcceptKink[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaXY:
    fhNBinsDcaXY=nbins+1;
    fhBinLimDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXY[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaZ:
    fhNBinsDcaZ=nbins+1;
    fhBinLimDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZ[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaXYnorm:
    fhNBinsDcaXYnorm=nbins+1;
    fhBinLimDcaXYnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaXYnorm[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kDcaZnorm:
    fhNBinsDcaZnorm=nbins+1;
    fhBinLimDcaZnorm=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimDcaZnorm[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kSigmaDcaXY:
    fhNBinsSigmaDcaXY=nbins+1;
    fhBinLimSigmaDcaXY=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimSigmaDcaXY[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;

  case kSigmaDcaZ:
    fhNBinsSigmaDcaZ=nbins+1;
    fhBinLimSigmaDcaZ=new Double_t[nbins+1];
    for(Int_t i=0;i<nbins+1;i++)fhBinLimSigmaDcaZ[i]=xmin+i*(xmax-xmin)/Double_t(nbins);
    break;
  }
}
//__________________________________________________________________________________
 void AliCFTrackIsPrimaryCuts::DefineHistograms() {
  //
  // histograms for cut variables, cut statistics and cut correlations
  //
  Int_t color = 2;

  // book cut statistics and cut correlation histograms
  fhCutStatistics = new TH1F(Form("%s_cut_statistics",GetName()), Form("%s cut statistics",GetName()), kNCuts,0.5,kNCuts+0.5);
  fhCutStatistics->SetLineWidth(2);
  fhCutStatistics->GetXaxis()->SetBinLabel(1,"dca xy");
  fhCutStatistics->GetXaxis()->SetBinLabel(2,"dca z");
  fhCutStatistics->GetXaxis()->SetBinLabel(3,"dca ellipse");
  fhCutStatistics->GetXaxis()->SetBinLabel(4,"n dca");
  fhCutStatistics->GetXaxis()->SetBinLabel(5,"sigma dca xy");
  fhCutStatistics->GetXaxis()->SetBinLabel(6,"sigma dca z");
  fhCutStatistics->GetXaxis()->SetBinLabel(7,"require dca");
  fhCutStatistics->GetXaxis()->SetBinLabel(8,"kink daughter");
  fhCutStatistics->GetXaxis()->SetBinLabel(9,"AOD type");

  fhCutCorrelation = new TH2F(Form("%s_cut_correlation",GetName()), Form("%s cut  correlation",GetName()), kNCuts,0.5,kNCuts+0.5,kNCuts,0.5,kNCuts+0.5);
  fhCutCorrelation->SetLineWidth(2);
  fhCutCorrelation->GetXaxis()->SetBinLabel(1,"dca xy");
  fhCutCorrelation->GetXaxis()->SetBinLabel(2,"dca z");
  fhCutCorrelation->GetXaxis()->SetBinLabel(3,"dca ellipse");
  fhCutCorrelation->GetXaxis()->SetBinLabel(4,"n dca");
  fhCutCorrelation->GetXaxis()->SetBinLabel(5,"sigma dca xy");
  fhCutCorrelation->GetXaxis()->SetBinLabel(6,"sigma dca z");
  fhCutCorrelation->GetXaxis()->SetBinLabel(7,"require dca");
  fhCutCorrelation->GetXaxis()->SetBinLabel(8,"kink daughter");
  fhCutCorrelation->GetXaxis()->SetBinLabel(9,"AOD type");

  fhCutCorrelation->GetYaxis()->SetBinLabel(1,"dca xy");
  fhCutCorrelation->GetYaxis()->SetBinLabel(2,"dca z");
  fhCutCorrelation->GetYaxis()->SetBinLabel(3,"dca ellipse");
  fhCutCorrelation->GetYaxis()->SetBinLabel(4,"n dca");
  fhCutCorrelation->GetYaxis()->SetBinLabel(5,"sigma dca xy");
  fhCutCorrelation->GetYaxis()->SetBinLabel(6,"sigma dca z");
  fhCutCorrelation->GetYaxis()->SetBinLabel(7,"require dca");
  fhCutCorrelation->GetYaxis()->SetBinLabel(8,"kink daughter");
  fhCutCorrelation->GetYaxis()->SetBinLabel(9,"AOD type");

  // book QA histograms
  Char_t str[5];
  for (Int_t i=0; i<kNStepQA; i++) {
    if (i==0) snprintf(str,5," ");
    else snprintf(str,5,"_cut");

    fhDcaXYvsDcaZ[i]            = new  TH2F(Form("%s_dcaXYvsDcaZ%s",GetName(),str),"",200,-10,10,200,-10,10);
    fhQA[kCutNSigmaToVertex][i]	= new TH1F(Form("%s_nSigmaToVertex%s",GetName(),str),"",fhNBinsNSigma-1,fhBinLimNSigma);
    fhQA[kCutRequireSigmaToVertex][i] = new TH1F(Form("%s_requireSigmaToVertex%s",GetName(),str),"",fhNBinsRequireSigma-1,fhBinLimRequireSigma);
    fhQA[kCutAcceptKinkDaughters][i] = new TH1F(Form("%s_acceptKinkDaughters%s",GetName(),str),"",fhNBinsAcceptKink-1,fhBinLimAcceptKink);
    fhQA[kDcaXY][i]		= new TH1F(Form("%s_dcaXY%s",GetName(),str),"",fhNBinsDcaXY-1,fhBinLimDcaXY);
    fhQA[kDcaZ][i]		= new TH1F(Form("%s_dcaZ%s",GetName(),str),"",fhNBinsDcaZ-1,fhBinLimDcaZ);
    fhQA[kDcaXYnorm][i]		= new TH1F(Form("%s_dcaXYnorm%s",GetName(),str),"",fhNBinsDcaXYnorm-1,fhBinLimDcaXYnorm);
    fhQA[kDcaZnorm][i]		= new TH1F(Form("%s_dcaZnorm%s",GetName(),str),"",fhNBinsDcaZnorm-1,fhBinLimDcaZnorm);
    fhQA[kSigmaDcaXY][i]	= new TH1F(Form("%s_sigmaDcaXY%s",GetName(),str),"",fhNBinsSigmaDcaXY-1,fhBinLimSigmaDcaXY);
    fhQA[kSigmaDcaZ][i]		= new TH1F(Form("%s_sigmaDcaZ%s",GetName(),str),"",fhNBinsSigmaDcaZ-1,fhBinLimSigmaDcaZ);

    fhDcaXYvsDcaZ[i]->SetXTitle("impact par. d_{z}");
    fhDcaXYvsDcaZ[i]->SetYTitle("impact par. d_{xy}");

    fhQA[kCutNSigmaToVertex][i]->SetXTitle("n #sigma to vertex");
    fhQA[kCutRequireSigmaToVertex][i]->SetXTitle("require #sigma to vertex");
    fhQA[kCutAcceptKinkDaughters][i]->SetXTitle("accept kink daughters");
    fhQA[kDcaXY][i]->SetXTitle("impact par. d_{xy}");
    fhQA[kDcaZ][i]->SetXTitle("impact par. d_{z}");
    fhQA[kDcaXYnorm][i]->SetXTitle("norm. impact par. d_{xy} / #sigma_{xy}");
    fhQA[kDcaZnorm][i]->SetXTitle("norm. impact par. d_{z} / #sigma_{z}");
    fhQA[kSigmaDcaXY][i]->SetXTitle("impact par. resolution #sigma_{xy}");
    fhQA[kSigmaDcaZ][i]->SetXTitle("impact par. resolution #sigma_{z}");
  }
  for(Int_t i=0; i<kNHist; i++) fhQA[i][1]->SetLineColor(color);
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::FillHistograms(TObject* obj, Bool_t f)
{
  //
  // fill the QA histograms
  //

  if (!obj) return;
  Bool_t isESDTrack = strcmp(obj->ClassName(),"AliESDtrack") == 0 ? kTRUE : kFALSE ;
  Bool_t isAODTrack = strcmp(obj->ClassName(),"AliAODTrack") == 0 ? kTRUE : kFALSE ;

  AliESDtrack * esdTrack = 0x0 ;
  AliAODTrack * aodTrack = 0x0 ; 
  if (isESDTrack) esdTrack = dynamic_cast<AliESDtrack*>(obj);
  if (isAODTrack) aodTrack = dynamic_cast<AliAODTrack*>(obj);

  // f = 0: fill histograms before cuts
  // f = 1: fill histograms after cuts

  // get the track to vertex parameter for ESD track
  if (esdTrack) {

    fhQA[kDcaZ][f]->Fill(fDCA[1]);
    fhQA[kDcaXY][f]->Fill(fDCA[0]);
    fhDcaXYvsDcaZ[f]->Fill(fDCA[1],fDCA[0]);
    fhQA[kSigmaDcaXY][f]->Fill(fDCA[2]);
    fhQA[kSigmaDcaZ][f]->Fill(fDCA[3]);
// // // // // // //  delete histograms
      fhQA[kDcaZnorm][f]->Fill(fDCA[1]);
      fhQA[kDcaXYnorm][f]->Fill(fDCA[0]);

    fhQA[kCutNSigmaToVertex][f]->Fill(fDCA[5]);
    if (fDCA[5]<0 && fRequireSigmaToVertex) 
      fhQA[kCutRequireSigmaToVertex][f]->Fill(0.);
    else 
      fhQA[kCutRequireSigmaToVertex][f]->Fill(1.);

    if (!fAcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0)
      fhQA[kCutAcceptKinkDaughters][f]->Fill(0.);
    else 
      fhQA[kCutAcceptKinkDaughters][f]->Fill(1.);
  }

  // fill cut statistics and cut correlation histograms with information from the bitmap
  if (f) return;

  // Get the bitmap of the single cuts
  SelectionBitMap(obj);

  // Number of single cuts in this class
  UInt_t ncuts = fBitmap->GetNbits();
  for(UInt_t bit=0; bit<ncuts;bit++) {
    if (!fBitmap->TestBitNumber(bit)) {
	fhCutStatistics->Fill(bit+1);
	for (UInt_t bit2=bit; bit2<ncuts;bit2++) {
	  if (!fBitmap->TestBitNumber(bit2)) 
	    fhCutCorrelation->Fill(bit+1,bit2+1);
	}
    }
  }
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::SaveHistograms(const Char_t* dir) {
  //
  // saves the histograms in a directory (dir)
  //
  if(!fIsQAOn) return;

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");

  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t j=0; j<kNStepQA; j++) {
    if (j==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");

    fhDcaXYvsDcaZ[j]    ->Write();

    for(Int_t i=0; i<kNHist; i++) fhQA[i][j]->Write();

    gDirectory->cd("../");
  }
  gDirectory->cd("../");
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::DrawHistograms()
{
  //
  // draws some histograms
  //
  if(!fIsQAOn) return;

  // pad margins
  Float_t right = 0.03;
  Float_t left = 0.175;
  Float_t top = 0.03;
  Float_t bottom = 0.175;

  TCanvas* canvas1 = new TCanvas("Track_QA_Primary_1", "Track QA Primary 1", 800, 500);
  canvas1->Divide(2, 1);

  canvas1->cd(1);
  fhCutStatistics->SetStats(kFALSE);
  fhCutStatistics->LabelsOption("v");
  gPad->SetLeftMargin(left);
  gPad->SetBottomMargin(0.25);
  gPad->SetRightMargin(right);
  gPad->SetTopMargin(0.1);
  fhCutStatistics->Draw();

  canvas1->cd(2);
  fhCutCorrelation->SetStats(kFALSE);
  fhCutCorrelation->LabelsOption("v");
  gPad->SetLeftMargin(0.30);
  gPad->SetRightMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.25);
  fhCutCorrelation->Draw("COLZ");

  canvas1->SaveAs(Form("%s.eps", canvas1->GetName()));
  canvas1->SaveAs(Form("%s.ps", canvas1->GetName()));

  // -----

  TCanvas* canvas2 = new TCanvas("Track_QA_Primary_2", "Track QA Primary 2", 800, 800);
  canvas2->Divide(2, 2);

  canvas2->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaXY][0]->SetStats(kFALSE);
  fhQA[kDcaXY][0]->Draw();
  fhQA[kDcaXY][1]->Draw("same");

  canvas2->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaZ][0]->SetStats(kFALSE);
  fhQA[kDcaZ][0]->Draw();
  fhQA[kDcaZ][1]->Draw("same");

  canvas2->cd(3);
//   fhDXYvsDZ[0]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZ[0]->Draw("COLZ");

  canvas2->cd(4);
//   fhDXYvsDZ[1]->SetStats(kFALSE);
  gPad->SetLogz();
  gPad->SetLeftMargin(bottom);
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(bottom);
  gPad->SetRightMargin(0.2);
  fhDcaXYvsDcaZ[1]->Draw("COLZ");

  canvas2->SaveAs(Form("%s.eps", canvas2->GetName()));
  canvas2->SaveAs(Form("%s.ps", canvas2->GetName()));

  // -----

  TCanvas* canvas3 = new TCanvas("Track_QA_Primary_3", "Track QA Primary 3", 800, 400);
  canvas3->Divide(2, 1);

  canvas3->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaXYnorm][0]->SetStats(kFALSE);
  fhQA[kDcaXYnorm][0]->Draw();
  fhQA[kDcaXYnorm][1]->Draw("same");

  canvas3->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kDcaZnorm][0]->SetStats(kFALSE);
  fhQA[kDcaZnorm][0]->Draw();
  fhQA[kDcaZnorm][1]->Draw("same");

  canvas3->SaveAs(Form("%s.eps", canvas3->GetName()));
  canvas3->SaveAs(Form("%s.ps", canvas3->GetName()));

  // -----

  TCanvas* canvas4 = new TCanvas("Track_QA_Primary_4", "Track QA Primary 4", 1200, 800);
  canvas4->Divide(3, 2);

  canvas4->cd(1);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kSigmaDcaXY][0]->SetStats(kFALSE);
  fhQA[kSigmaDcaXY][0]->Draw();
  fhQA[kSigmaDcaXY][1]->Draw("same");

  canvas4->cd(2);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kSigmaDcaZ][0]->SetStats(kFALSE);
  fhQA[kSigmaDcaZ][0]->Draw();
  fhQA[kSigmaDcaZ][1]->Draw("same");

  canvas4->cd(4);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutNSigmaToVertex][0]->SetStats(kFALSE);
  fhQA[kCutNSigmaToVertex][0]->Draw();
  fhQA[kCutNSigmaToVertex][1]->Draw("same");

  canvas4->cd(5);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutRequireSigmaToVertex][0]->SetStats(kFALSE);
  fhQA[kCutRequireSigmaToVertex][0]->Draw();
  fhQA[kCutRequireSigmaToVertex][1]->Draw("same");

  canvas4->cd(6);
  gPad->SetRightMargin(right);
  gPad->SetLeftMargin(left);
  gPad->SetTopMargin(top);
  gPad->SetBottomMargin(bottom);
  fhQA[kCutAcceptKinkDaughters][0]->SetStats(kFALSE);
  fhQA[kCutAcceptKinkDaughters][0]->Draw();
  fhQA[kCutAcceptKinkDaughters][1]->Draw("same");

  canvas4->SaveAs(Form("%s.eps", canvas4->GetName()));
  canvas4->SaveAs(Form("%s.ps", canvas4->GetName()));
}
//__________________________________________________________________________________
void AliCFTrackIsPrimaryCuts::AddQAHistograms(TList *qaList) {
  //
  // saves the histograms in a TList
  //
  DefineHistograms();

  qaList->Add(fhCutStatistics);
  qaList->Add(fhCutCorrelation);

  for (Int_t j=0; j<kNStepQA; j++) {
    qaList->Add(fhDcaXYvsDcaZ[j]);
    for(Int_t i=0; i<kNHist; i++)
	qaList->Add(fhQA[i][j]);
  }
}
