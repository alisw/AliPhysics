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

/*****************************************************************
  AliFlowEventSimple: A simple event
  for flow analysis

  origin: Naomi van der Kolk (kolk@nikhef.nl)
          Ante Bilandzic     (anteb@nikhef.nl)
          Raimond Snellings  (Raimond.Snellings@nikhef.nl)
  mods:   Mikolaj Krzewicki  (mikolaj.krzewicki@cern.ch)
          Redmer A. Bertens  (rbertens@cern.ch)
*****************************************************************/

#include "Riostream.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TList.h"
#include "TTree.h"
#include "TParticle.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TProfile.h"
#include "TParameter.h"
#include "TBrowser.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowTrackSimpleCuts.h"
#include "AliFlowEventSimple.h"
#include "TRandom.h"
#include <random>

using std::cout;
using std::endl;
using std::cerr;
ClassImp(AliFlowEventSimple)

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple():
  fTrackCollection(NULL),
  fReferenceMultiplicity(0),
  fNumberOfTracks(0),
  fUseGlauberMCSymmetryPlanes(kFALSE),
  fUseExternalSymmetryPlanes(kFALSE),
  fPsi1(0.),
  fPsi2(0.),
  fPsi3(0.),
  fPsi4(0.),
  fPsi5(0.),
  fPsi1Psi3(0x0),
  fPsi2Psi4(0x0),
  fPsi3Psi5(0x0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fUserModified(kFALSE),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fNumberOfPOIsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL),
  fShuffledIndexes(NULL),
  fShuffleTracks(kFALSE),
  fMothersCollection(NULL),
  fCentrality(-1.),
  fCentralityCL1(-1.),
  fNITSCL1(-1.),
  fCentralityTRK(-1.),
  fRun(-1),
  fZNCQ0(0.),
  fZNAQ0(0.),
  fZNCM(0.),
  fZNAM(0.),
  fZPCM(0.),
  fZPAM(0.),
  fAbsOrbit(0),
  fNumberOfPOItypes(2),
  fNumberOfPOIs(NULL)
{
  fZNCQ = AliFlowVector();
  fZNAQ = AliFlowVector();
  for(Int_t i(0); i < 3; i++) {
    fVtxPos[i] = 0.;
  }
  for(Int_t i=0; i < 4; i++) {
    fV0C[i] = AliFlowVector();
    fV0A[i] = AliFlowVector();
  }
  cout << "AliFlowEventSimple: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple( Int_t n,
                                        ConstructionMethod method,
                                        TF1* ptDist,
                                        Double_t phiMin,
                                        Double_t phiMax,
                                        Double_t etaMin,
                                        Double_t etaMax):
  fTrackCollection(new TObjArray(n)),
  fReferenceMultiplicity(0),
  fNumberOfTracks(0),
  fUseGlauberMCSymmetryPlanes(kFALSE),
  fUseExternalSymmetryPlanes(kFALSE),
  fPsi1(0.),
  fPsi2(0.),
  fPsi3(0.),
  fPsi4(0.),
  fPsi5(0.),
  fPsi1Psi3(0x0),
  fPsi2Psi4(0x0),
  fPsi3Psi5(0x0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fUserModified(kFALSE),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fNumberOfPOIsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL),
  fShuffledIndexes(NULL),
  fShuffleTracks(kFALSE),
  fMothersCollection(new TObjArray()),
  fCentrality(-1.),
  fCentralityCL1(-1.),
  fNITSCL1(-1.),
  fCentralityTRK(-1.),
  fRun(-1),
  fZNCQ0(0.),
  fZNAQ0(0.),
  fZNCM(0.),
  fZNAM(0.),
  fZPCM(0.),
  fZPAM(0.),
  fAbsOrbit(0),
  fNumberOfPOItypes(2),
  fNumberOfPOIs(new Int_t[fNumberOfPOItypes])
{
  //ctor
  // if second argument is set to AliFlowEventSimple::kGenerate
  // it generates n random tracks with given Pt distribution
  // (a sane default is provided), phi and eta are uniform

  if (method==kGenerate)
    Generate(n,ptDist,phiMin,phiMax,etaMin,etaMax);

  fZNCQ = AliFlowVector();
  fZNAQ = AliFlowVector();
  for(Int_t i(0); i < 3; i++) {
    fVtxPos[i] = 0.;
  }
  for(Int_t i=0; i < 4; i++) {
    fV0C[i] = AliFlowVector();
    fV0A[i] = AliFlowVector();
  }
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple(const AliFlowEventSimple& anEvent):
  TObject(anEvent),
  fTrackCollection((TObjArray*)(anEvent.fTrackCollection)->Clone()),
  fReferenceMultiplicity(anEvent.fReferenceMultiplicity),
  fNumberOfTracks(anEvent.fNumberOfTracks),
  fUseGlauberMCSymmetryPlanes(anEvent.fUseGlauberMCSymmetryPlanes),
  fUseExternalSymmetryPlanes(anEvent.fUseExternalSymmetryPlanes),
  fPsi1(anEvent.fPsi1),
  fPsi2(anEvent.fPsi2),
  fPsi3(anEvent.fPsi3),
  fPsi4(anEvent.fPsi4),
  fPsi5(anEvent.fPsi5),
  fPsi1Psi3(anEvent.fPsi1Psi3),
  fPsi2Psi4(anEvent.fPsi2Psi4),
  fPsi3Psi5(anEvent.fPsi3Psi5),
  fMCReactionPlaneAngle(anEvent.fMCReactionPlaneAngle),
  fMCReactionPlaneAngleIsSet(anEvent.fMCReactionPlaneAngleIsSet),
  fAfterBurnerPrecision(anEvent.fAfterBurnerPrecision),
  fUserModified(anEvent.fUserModified),
  fNumberOfTracksWrap(anEvent.fNumberOfTracksWrap),
  fNumberOfRPsWrap(anEvent.fNumberOfRPsWrap),
  fNumberOfPOIsWrap(anEvent.fNumberOfPOIsWrap),
  fMCReactionPlaneAngleWrap(anEvent.fMCReactionPlaneAngleWrap),
  fShuffledIndexes(NULL),
  fShuffleTracks(anEvent.fShuffleTracks),
  fMothersCollection(new TObjArray()),
  fCentrality(anEvent.fCentrality),
  fCentralityCL1(anEvent.fCentralityCL1),
  fNITSCL1(anEvent.fNITSCL1),
  fCentralityTRK(anEvent.fCentralityTRK),
  fRun(anEvent.fRun),
  fZNCQ(anEvent.fZNCQ),
  fZNAQ(anEvent.fZNAQ),
  fZNCQ0(anEvent.fZNCQ0),
  fZNAQ0(anEvent.fZNAQ0),
  fZNCM(anEvent.fZNCM),
  fZNAM(anEvent.fZNAM),
  fZPCM(anEvent.fZPCM),
  fZPAM(anEvent.fZPAM),
  fAbsOrbit(anEvent.fAbsOrbit),
  fNumberOfPOItypes(anEvent.fNumberOfPOItypes),
  fNumberOfPOIs(new Int_t[fNumberOfPOItypes])
{
  //copy constructor
  memcpy(fNumberOfPOIs,anEvent.fNumberOfPOIs,fNumberOfPOItypes*sizeof(Int_t));
  for(Int_t i(0); i < 3; i++) {
    fVtxPos[i] = anEvent.fVtxPos[i];
  }
  for(Int_t i=0; i < 4; i++) {
    fV0C[i] = anEvent.fV0C[i];
    fV0A[i] = anEvent.fV0A[i];
  }
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::SetNumberOfPOIs( Int_t numberOfPOIs, Int_t poiType)
{
  //set the number of poi classes, resize the array if larger is needed
  //never shrink the array
  //never decrease the stored number
  if (poiType>=fNumberOfPOItypes)
  {
    Int_t n = poiType+1;
    Int_t* tmp = new Int_t[n];
    for (Int_t j=0; j<n; j++) { tmp[j]=0; }
    memcpy(tmp,fNumberOfPOIs,fNumberOfPOItypes*sizeof(Int_t));
    delete [] fNumberOfPOIs;
    fNumberOfPOIs = tmp;
    fNumberOfPOItypes = n;
  }

  fNumberOfPOIs[poiType] = numberOfPOIs;
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::IncrementNumberOfPOIs(Int_t poiType)
{

  if (poiType>=fNumberOfPOItypes) SetNumberOfPOIs(0,poiType);
  fNumberOfPOIs[poiType]++;
}

//-----------------------------------------------------------------------
AliFlowEventSimple& AliFlowEventSimple::operator=(const AliFlowEventSimple& anEvent)
{
  //assignment operator
  if (&anEvent==this) return *this; //check self-assignment
  if (fTrackCollection) fTrackCollection->Delete();
  delete fTrackCollection;
  fTrackCollection = (TObjArray*)(anEvent.fTrackCollection)->Clone(); //deep copy
  fReferenceMultiplicity = anEvent.fReferenceMultiplicity;
  fNumberOfTracks = anEvent.fNumberOfTracks;
  fNumberOfPOItypes = anEvent.fNumberOfPOItypes;
  delete [] fNumberOfPOIs;
  fNumberOfPOIs=new Int_t[fNumberOfPOItypes];
  memcpy(fNumberOfPOIs,anEvent.fNumberOfPOIs,fNumberOfPOItypes*sizeof(Int_t));
  fUseGlauberMCSymmetryPlanes = anEvent.fUseGlauberMCSymmetryPlanes;
  fUseExternalSymmetryPlanes = anEvent.fUseExternalSymmetryPlanes;
  fPsi1 = anEvent.fPsi1;
  fPsi2 = anEvent.fPsi2;
  fPsi3 = anEvent.fPsi3;
  fPsi4 = anEvent.fPsi4;
  fPsi5 = anEvent.fPsi5;
  fPsi1Psi3 = anEvent.fPsi1Psi3;
  fPsi2Psi4 = anEvent.fPsi2Psi4;
  fPsi3Psi5 = anEvent.fPsi3Psi5;
  fMCReactionPlaneAngle = anEvent.fMCReactionPlaneAngle;
  fMCReactionPlaneAngleIsSet = anEvent.fMCReactionPlaneAngleIsSet;
  fAfterBurnerPrecision = anEvent.fAfterBurnerPrecision;
  fUserModified=anEvent.fUserModified;
  fNumberOfTracksWrap = anEvent.fNumberOfTracksWrap;
  fNumberOfRPsWrap = anEvent.fNumberOfRPsWrap;
  fNumberOfPOIsWrap = anEvent.fNumberOfPOIsWrap;
  fMCReactionPlaneAngleWrap = anEvent.fMCReactionPlaneAngleWrap;
  fShuffleTracks = anEvent.fShuffleTracks;
  fCentrality = anEvent.fCentrality;
  fCentralityCL1 = anEvent.fCentralityCL1;
  fNITSCL1 = anEvent.fNITSCL1;
  fCentralityTRK = anEvent.fCentralityTRK;
  fRun = anEvent.fRun;
  fZNCQ = anEvent.fZNCQ;
  fZNAQ = anEvent.fZNAQ;
  fZNCQ0 = anEvent.fZNCQ0;
  fZNAQ0 = anEvent.fZNAQ0;
  fZNCM = anEvent.fZNCM;
  fZNAM = anEvent.fZNAM;
  fZPCM = anEvent.fZPCM;
  fZPAM = anEvent.fZPAM;
  fAbsOrbit = anEvent.fAbsOrbit;
  for(Int_t i(0); i < 3; i++) {
    fVtxPos[i] = anEvent.fVtxPos[i];
  }
  for(Int_t i=0; i < 4; i++) {
    fV0C[i] = anEvent.fV0C[i];
    fV0A[i] = anEvent.fV0A[i];
  }
  delete [] fShuffledIndexes;
  return *this;
}

//-----------------------------------------------------------------------
AliFlowEventSimple::~AliFlowEventSimple()
{
  //destructor
  if (fTrackCollection) fTrackCollection->Delete();
  delete fTrackCollection;
  delete fNumberOfTracksWrap;
  delete fNumberOfRPsWrap;
  delete fNumberOfPOIsWrap;
  delete fMCReactionPlaneAngleWrap;
  delete fShuffledIndexes;
  delete fMothersCollection;
  delete [] fNumberOfPOIs;
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::SetUseExternalSymmetryPlanes(TF1 *gPsi1Psi3,
						      TF1 *gPsi2Psi4,
						      TF1 *gPsi3Psi5) {
  //Use symmetry planes, setup correlations between different Psi_n
  fUseExternalSymmetryPlanes = kTRUE;

  //Correlations between Psi_1 and Psi_3
  if(gPsi1Psi3) fPsi1Psi3 = gPsi1Psi3;
  else {
    fPsi1Psi3 = new TF1("fPsi1Psi3","[0]*x+[1]",0.,2.*TMath::Pi());
    fPsi1Psi3->SetParameter(0,1.);
    fPsi1Psi3->SetParameter(1,0.);
  }

  //Correlations between Psi_2 and Psi_4
  if(gPsi2Psi4) fPsi2Psi4 = gPsi2Psi4;
  else {
    fPsi2Psi4 = new TF1("fPsi2Psi4","[0]*x+[1]",0.,2.*TMath::Pi());
    fPsi2Psi4->SetParameter(0,1.);
    fPsi2Psi4->SetParameter(1,0.);
  }

  //Correlations between Psi_3 and Psi_5
  if(gPsi3Psi5) fPsi3Psi5 = gPsi3Psi5;
  else {
    fPsi3Psi5 = new TF1("fPsi3Psi5","[0]*x+[1]",0.,2.*TMath::Pi());
    fPsi3Psi5->SetParameter(0,1.);
    fPsi3Psi5->SetParameter(1,0.);
  }
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::Generate(Int_t nParticles,
                                  TF1* ptDist,
                                  Double_t phiMin,
                                  Double_t phiMax,
                                  Double_t etaMin,
                                  Double_t etaMax)
{
  //generate nParticles random tracks uniform in phi and eta
  //according to the specified pt distribution
  if (!ptDist)
  {
    static TF1 ptdistribution("ptSpectra","x*TMath::Exp(-pow(0.13957*0.13957+x*x,0.5)/0.4)",0.1,10.);
    ptDist=&ptdistribution;
  }

  for (Int_t i=0; i<nParticles; i++)
  {
    AliFlowTrackSimple* track = new AliFlowTrackSimple();
    track->SetPhi( gRandom->Uniform(phiMin,phiMax) );
    track->SetEta( gRandom->Uniform(etaMin,etaMax) );
    track->SetPt( ptDist->GetRandom() );
    track->SetCharge( (gRandom->Uniform()-0.5<0)?-1:1 );
    AddTrack(track);
  }
  if(fUseExternalSymmetryPlanes) {
    Double_t betaParameter = gRandom->Gaus(0.,1.3);
    fPsi1Psi3->SetParameter(1,betaParameter);

    betaParameter = gRandom->Gaus(0.,0.9);
    fPsi2Psi4->SetParameter(1,betaParameter);

    betaParameter = gRandom->Gaus(0.,1.5);
    fPsi3Psi5->SetParameter(1,betaParameter);

    fPsi1 = gRandom->Uniform(2.*TMath::Pi());
    fPsi2 = gRandom->Uniform(2.*TMath::Pi());
    fPsi3 = fPsi1Psi3->Eval(fPsi1);
    if(fPsi3 < 0) fPsi3 += 2.*TMath::Pi();
    else if(fPsi3 > 2.*TMath::Pi()) fPsi3 -= 2.*TMath::Pi();
    fPsi4 = fPsi2Psi4->Eval(fPsi2);
    if(fPsi4 < 0) fPsi4 += 2.*TMath::Pi();
    else if(fPsi4 > 2.*TMath::Pi()) fPsi4 -= 2.*TMath::Pi();
    fPsi5 = fPsi3Psi5->Eval(fPsi3);
    if(fPsi5 < 0) fPsi5 += 2.*TMath::Pi();
    else if(fPsi5 > 2.*TMath::Pi()) fPsi5 -= 2.*TMath::Pi();

    fMCReactionPlaneAngle=fPsi2;
    fMCReactionPlaneAngleIsSet=kTRUE;
  }
  else {
    fMCReactionPlaneAngle=gRandom->Uniform(0.0,TMath::TwoPi());
    fMCReactionPlaneAngleIsSet=kTRUE;
  }
  SetUserModified();
}

//-----------------------------------------------------------------------
AliFlowTrackSimple* AliFlowEventSimple::GetTrack(Int_t i)
{
  //get track i from collection
  if (i>=fNumberOfTracks) return NULL;
  Int_t trackIndex=i;
  //if asked use the shuffled index
  if (fShuffleTracks)
  {
    if (!fShuffledIndexes) ShuffleTracks();
    trackIndex=fShuffledIndexes[i];
  }
  AliFlowTrackSimple* pTrack = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(trackIndex)) ;
  return pTrack;
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::ShuffleTracks()
{
  //shuffle track indexes
  if (!fShuffledIndexes)
  {
    //initialize the table with shuffled indexes
    fShuffledIndexes = new Int_t[fNumberOfTracks];
    for (Int_t j=0; j<fNumberOfTracks; j++) { fShuffledIndexes[j]=j; }
  }
  //shuffle
  std::random_device rd;
  std::default_random_engine engine{rd()};
  std::shuffle(&fShuffledIndexes[0], &fShuffledIndexes[fNumberOfTracks],engine);
  Printf("Tracks shuffled! tracks: %i",fNumberOfTracks);
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::AddTrack( AliFlowTrackSimple* track )
{
  //add a track, delete the old one if necessary
  if (fNumberOfTracks < fTrackCollection->GetEntriesFast())
  {
    TObject* o = fTrackCollection->At(fNumberOfTracks);
    delete o;
  }
  fTrackCollection->AddAtAndExpand(track,fNumberOfTracks);
  if (track->GetNDaughters()>0)
  {
    //if there track has daughters cache in the collection of mothers
    fMothersCollection->Add(track);
  }
  TrackAdded();
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::TrackAdded()
{
  //book keeping after a new track has been added
  fNumberOfTracks++;
  if (fShuffledIndexes)
  {
    delete [] fShuffledIndexes;
    fShuffledIndexes=NULL;
  }
}

//-----------------------------------------------------------------------
AliFlowTrackSimple* AliFlowEventSimple::MakeNewTrack()
{
   AliFlowTrackSimple *t=dynamic_cast<AliFlowTrackSimple *>(fTrackCollection->RemoveAt(fNumberOfTracks));
   if( !t ) {  // If there was no track at the end of the list then create a new track
      t=new AliFlowTrackSimple();
   }

   return t;
}

//-----------------------------------------------------------------------
AliFlowVector AliFlowEventSimple::GetQ( Int_t n,
                                        TList *weightsList,
                                        Bool_t usePhiWeights,
                                        Bool_t usePtWeights,
                                        Bool_t useEtaWeights )
{
  // calculate Q-vector in harmonic n without weights (default harmonic n=2)
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  AliFlowVector vQ;
  vQ.Set(0.,0.);

  Int_t iOrder = n;
  Double_t sumOfWeights = 0.;
  Double_t dPhi = 0.;
  Double_t dPt = 0.;
  Double_t dEta = 0.;
  Double_t dWeight = 1.;

  AliFlowTrackSimple* pTrack = NULL;

  Int_t nBinsPhi = 0;
  Double_t dBinWidthPt = 0.;
  Double_t dPtMin = 0.;
  Double_t dBinWidthEta = 0.;
  Double_t dEtaMin = 0.;

  Double_t wPhi = 1.; // weight Phi
  Double_t wPt = 1.;  // weight Pt
  Double_t wEta = 1.; // weight Eta

  TH1F *phiWeights = NULL;
  TH1D *ptWeights  = NULL;
  TH1D *etaWeights = NULL;

  if(weightsList)
  {
    if(usePhiWeights)
    {
      phiWeights = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights"));
      if(phiWeights) nBinsPhi = phiWeights->GetNbinsX();
    }
    if(usePtWeights)
    {
      ptWeights = dynamic_cast<TH1D *>(weightsList->FindObject("pt_weights"));
      if(ptWeights)
      {
        dBinWidthPt = ptWeights->GetBinWidth(1); // assuming that all bins have the same width
        dPtMin = (ptWeights->GetXaxis())->GetXmin();
      }
    }
    if(useEtaWeights)
    {
      etaWeights = dynamic_cast<TH1D *>(weightsList->FindObject("eta_weights"));
      if(etaWeights)
      {
        dBinWidthEta = etaWeights->GetBinWidth(1); // assuming that all bins have the same width
        dEtaMin = (etaWeights->GetXaxis())->GetXmin();
      }
    }
  } // end of if(weightsList)

  // loop over tracks
  for(Int_t i=0; i<fNumberOfTracks; i++)
  {
    pTrack = (AliFlowTrackSimple*)fTrackCollection->At(i);
    if(pTrack)
    {
      if(pTrack->InRPSelection())
      {
        dPhi = pTrack->Phi();
        dPt  = pTrack->Pt();
        dEta = pTrack->Eta();
	      dWeight = pTrack->Weight();

        // determine Phi weight: (to be improved, I should here only access it + the treatment of gaps in the if statement)
        if(phiWeights && nBinsPhi)
        {
          wPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
        }
        // determine v'(pt) weight:
        if(ptWeights && dBinWidthPt)
        {
          wPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt)));
        }
        // determine v'(eta) weight:
        if(etaWeights && dBinWidthEta)
        {
          wEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta)));
        }

        // building up the weighted Q-vector:
        dQX += dWeight*wPhi*wPt*wEta*TMath::Cos(iOrder*dPhi);
        dQY += dWeight*wPhi*wPt*wEta*TMath::Sin(iOrder*dPhi);

        // weighted multiplicity:
        sumOfWeights += dWeight*wPhi*wPt*wEta;

      } // end of if (pTrack->InRPSelection())
    } // end of if (pTrack)
    else
    {
      cerr << "no particle!!!"<<endl;
    }
  } // loop over particles

  vQ.Set(dQX,dQY);
  vQ.SetMult(sumOfWeights);
  vQ.SetHarmonic(iOrder);
  vQ.SetPOItype(AliFlowTrackSimple::kRP);
  vQ.SetSubeventNumber(-1);

  return vQ;

}

//-----------------------------------------------------------------------
void AliFlowEventSimple::Get2Qsub( AliFlowVector* Qarray,
                                   Int_t n,
                                   TList *weightsList,
                                   Bool_t usePhiWeights,
                                   Bool_t usePtWeights,
                                   Bool_t useEtaWeights )
{

  // calculate Q-vector in harmonic n without weights (default harmonic n=2)
  Double_t dQX = 0.;
  Double_t dQY = 0.;

  Int_t iOrder = n;
  Double_t sumOfWeights = 0.;
  Double_t dPhi = 0.;
  Double_t dPt  = 0.;
  Double_t dEta = 0.;
  Double_t dWeight = 1.;

  AliFlowTrackSimple* pTrack = NULL;

  Int_t    iNbinsPhiSub0 = 0;
  Int_t    iNbinsPhiSub1 = 0;
  Double_t dBinWidthPt = 0.;
  Double_t dPtMin      = 0.;
  Double_t dBinWidthEta= 0.;
  Double_t dEtaMin     = 0.;

  Double_t dWphi = 1.;  // weight Phi
  Double_t dWpt  = 1.;  // weight Pt
  Double_t dWeta = 1.;  // weight Eta

  TH1F* phiWeightsSub0 = NULL;
  TH1F* phiWeightsSub1 = NULL;
  TH1D* ptWeights  = NULL;
  TH1D* etaWeights = NULL;

  if(weightsList)
  {
    if(usePhiWeights)
    {
      phiWeightsSub0 = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights_sub0"));
      if(phiWeightsSub0) {
	iNbinsPhiSub0 = phiWeightsSub0->GetNbinsX();
      }
      phiWeightsSub1 = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights_sub1"));
      if(phiWeightsSub1) {
	iNbinsPhiSub1 = phiWeightsSub1->GetNbinsX();
      }
    }
    if(usePtWeights)
    {
      ptWeights = dynamic_cast<TH1D *>(weightsList->FindObject("pt_weights"));
      if(ptWeights)
      {
        dBinWidthPt = ptWeights->GetBinWidth(1); // assuming that all bins have the same width
        dPtMin = (ptWeights->GetXaxis())->GetXmin();
      }
    }
    if(useEtaWeights)
    {
      etaWeights = dynamic_cast<TH1D *>(weightsList->FindObject("eta_weights"));
      if(etaWeights)
      {
        dBinWidthEta = etaWeights->GetBinWidth(1); // assuming that all bins have the same width
        dEtaMin = (etaWeights->GetXaxis())->GetXmin();
      }
    }
  } // end of if(weightsList)

  //loop over the two subevents
  for (Int_t s=0; s<2; s++)
  {
    // loop over tracks
    for(Int_t i=0; i<fNumberOfTracks; i++)
    {
      pTrack = (AliFlowTrackSimple*)fTrackCollection->At(i);
      if(!pTrack)
      {
        cerr << "no particle!!!"<<endl;
        continue;
      }
      if(pTrack->InRPSelection() && (pTrack->InSubevent(s)))
      {
        dPhi    = pTrack->Phi();
        dPt     = pTrack->Pt();
        dEta    = pTrack->Eta();
        dWeight = pTrack->Weight();

        // determine Phi weight: (to be improved, I should here only access it + the treatment of gaps in the if statement)
        //subevent 0
        if(s == 0)  {
          if(phiWeightsSub0 && iNbinsPhiSub0)  {
            Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNbinsPhiSub0/TMath::TwoPi()));
            //use the phi value at the center of the bin
            dPhi  = phiWeightsSub0->GetBinCenter(phiBin);
            dWphi = phiWeightsSub0->GetBinContent(phiBin);
          }
        }
        //subevent 1
        else if (s == 1) {
          if(phiWeightsSub1 && iNbinsPhiSub1) {
            Int_t phiBin = 1+(Int_t)(TMath::Floor(dPhi*iNbinsPhiSub1/TMath::TwoPi()));
            //use the phi value at the center of the bin
            dPhi  = phiWeightsSub1->GetBinCenter(phiBin);
            dWphi = phiWeightsSub1->GetBinContent(phiBin);
          }
        }

        // determine v'(pt) weight:
        if(ptWeights && dBinWidthPt)
        {
          dWpt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt)));
        }

        // determine v'(eta) weight:
        if(etaWeights && dBinWidthEta)
        {
          dWeta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta)));
        }

        // building up the weighted Q-vector:
        dQX += dWeight*dWphi*dWpt*dWeta*TMath::Cos(iOrder*dPhi);
        dQY += dWeight*dWphi*dWpt*dWeta*TMath::Sin(iOrder*dPhi);

        // weighted multiplicity:
        sumOfWeights+=dWeight*dWphi*dWpt*dWeta;

      } // end of if (pTrack->InRPSelection())
    } // loop over particles

    Qarray[s].Set(dQX,dQY);
    Qarray[s].SetMult(sumOfWeights);
    Qarray[s].SetHarmonic(iOrder);
    Qarray[s].SetPOItype(AliFlowTrackSimple::kRP);
    Qarray[s].SetSubeventNumber(s);

    //reset
    sumOfWeights = 0.;
    dQX = 0.;
    dQY = 0.;
  }

}

//------------------------------------------------------------------------------

void AliFlowEventSimple::GetZDC2Qsub(AliFlowVector* Qarray)
{
  Qarray[0] = fZNCQ;
  Qarray[1] = fZNAQ;
}

//------------------------------------------------------------------------------

void AliFlowEventSimple::SetZDC2Qsub(Double_t* QVC, Double_t MC, Double_t* QVA, Double_t MA)
{
  fZNCQ = AliFlowVector(QVC[0],QVC[1],MC,1);
  fZNAQ = AliFlowVector(QVA[0],QVA[1],MA,1);
}

//------------------------------------------------------------------------------

void AliFlowEventSimple::GetV02Qsub(AliFlowVector* Qarray, Int_t har)
{
  if(har>0 && har<5) {
    Qarray[0] = fV0C[har-1];
    Qarray[1] = fV0A[har-1];
  } else {
    printf("WARNING: harmonic %d not available \n",har);
  }
}

//------------------------------------------------------------------------------

void AliFlowEventSimple::SetV02Qsub(Double_t QVCx, Double_t QVCy, Double_t MC, Double_t QVAx, Double_t QVAy, Double_t MA, Int_t har)
{
  if(har>0 && har<5) {
    fV0C[har-1] = AliFlowVector(QVCx,QVCy,MC,har);
    fV0A[har-1] = AliFlowVector(QVAx,QVAy,MA,har);
  }
}

//------------------------------------------------------------------------------

void AliFlowEventSimple::GetVertexPosition(Double_t* pos)
{
  pos[0] = fVtxPos[0];
  pos[1] = fVtxPos[1];
  pos[2] = fVtxPos[2];
}

//------------------------------------------------------------------------------

void AliFlowEventSimple::SetVertexPosition(Double_t* pos)
{
  fVtxPos[0] = pos[0];
  fVtxPos[1] = pos[1];
  fVtxPos[2] = pos[2];
}
//------------------------------------------------------------------------------

void AliFlowEventSimple::Print(Option_t *option) const
{
  //   -*-*-*-*-*Print some global quantities for this histogram collection class *-*-*-*-*-*-*-*
  //             ===============================================
  //   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
  printf( "Class.Print Name = %s, #tracks= %d, Number of RPs= %d, Number of POIs = %d, MC EventPlaneAngle= %f\n",
          GetName(),fNumberOfTracks, fNumberOfPOIs[0], fNumberOfPOIs[1], fMCReactionPlaneAngle );

  TString optionstr(option);
  if (!optionstr.Contains("all")) return;
  if (fTrackCollection)
  {
    fTrackCollection->Print(option);
  }
  else
  {
    printf( "Empty track collection \n");
  }
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::Browse(TBrowser *b)
{
  //browse in TBrowser
  if (!b) return;
  if (!fNumberOfTracksWrap)
  {
    fNumberOfTracksWrap = new TParameter<int>("fNumberOfTracks", fNumberOfTracks);
    b->Add(fNumberOfTracksWrap);
  }
  if (!fNumberOfRPsWrap)
  {
    fNumberOfRPsWrap = new TParameter<int>("fNumberOfRPs", GetNumberOfRPs());
    b->Add(fNumberOfRPsWrap);
  }
  if (!fNumberOfPOIsWrap)
  {
    fNumberOfPOIsWrap = new TParameter<int>("fNumberOfPOIs", GetNumberOfPOIs());
    b->Add(fNumberOfPOIsWrap);
  }
  if (!fMCReactionPlaneAngleWrap)
  {
    fMCReactionPlaneAngleWrap = new TParameter<double>(" fMCReactionPlaneAngle",  fMCReactionPlaneAngle);
    b->Add( fMCReactionPlaneAngleWrap);
  }
  if (fTrackCollection) b->Add(fTrackCollection,"AliFlowTracksSimple");
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple( TTree* inputTree,
                                        const AliFlowTrackSimpleCuts* rpCuts,
                                        const AliFlowTrackSimpleCuts* poiCuts):
  fTrackCollection(NULL),
  fReferenceMultiplicity(0),
  fNumberOfTracks(0),
  fUseGlauberMCSymmetryPlanes(kFALSE),
  fUseExternalSymmetryPlanes(kFALSE),
  fPsi1(0.),
  fPsi2(0.),
  fPsi3(0.),
  fPsi4(0.),
  fPsi5(0.),
  fPsi1Psi3(0x0),
  fPsi2Psi4(0x0),
  fPsi3Psi5(0x0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fUserModified(kFALSE),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fNumberOfPOIsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL),
  fShuffledIndexes(NULL),
  fShuffleTracks(kFALSE),
  fMothersCollection(new TObjArray()),
  fCentrality(-1.),
  fCentralityCL1(-1.),
  fNITSCL1(-1.),
  fCentralityTRK(-1.),
  fRun(-1),
  fZNCQ0(0.),
  fZNAQ0(0.),
  fZNCM(0.),
  fZNAM(0.),
  fZPCM(0.),
  fZPAM(0.),
  fAbsOrbit(0),
  fNumberOfPOItypes(2),
  fNumberOfPOIs(new Int_t[fNumberOfPOItypes])
{
  //constructor, fills the event from a TTree of kinematic.root files
  //applies RP and POI cuts, tags the tracks

  Int_t numberOfInputTracks = inputTree->GetEntries() ;
  fTrackCollection = new TObjArray(numberOfInputTracks/2);

  TParticle* pParticle = new TParticle();
  inputTree->SetBranchAddress("Particles",&pParticle);

  for (Int_t i=0; i<numberOfInputTracks; i++)
  {
    inputTree->GetEntry(i);   //get input particle

    if (!pParticle) continue; //no particle
    if (!pParticle->IsPrimary()) continue;

    Bool_t rpOK = (rpCuts->PassesCuts(pParticle)>0);
    Bool_t poiOK = poiCuts->PassesCuts(pParticle);
    Int_t poiType = poiCuts->GetPOItype();

    if (rpOK || poiOK)
    {
      AliFlowTrackSimple* pTrack = new AliFlowTrackSimple(pParticle);

      //marking the particles used for int. flow:
      if(rpOK)
      {
        pTrack->TagRP(kTRUE);
        IncrementNumberOfPOIs(0);
        cout<<"numberOfRPs = "<<fNumberOfPOIs[0]<<endl;
      }
      //marking the particles used for diff. flow:
      if(poiOK)
      {
        pTrack->Tag(poiType);
        IncrementNumberOfPOIs(poiType);
        printf("fNumberOfPOIs[%i] = %i",poiType,fNumberOfPOIs[poiType]);
      }
      //adding a particles which were used either for int. or diff. flow to the list
      AddTrack(pTrack);
    }
  }//for i
  delete pParticle;

  fZNCQ = AliFlowVector();
  fZNAQ = AliFlowVector();
  for(Int_t i(0); i < 3; i++) {
    fVtxPos[i] = 0.;
  }
  for(Int_t i=0; i < 4; i++) {
    fV0C[i] = AliFlowVector();
    fV0A[i] = AliFlowVector();
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::CloneTracks(Int_t n)
{
  //clone every track n times to add non-flow
  if (n<=0) return; //no use to clone stuff zero or less times
  Int_t ntracks = fNumberOfTracks;
  fTrackCollection->Expand((n+1)*fNumberOfTracks);
  for (Int_t i=0; i<n; i++)
  {
    for (Int_t itrack=0; itrack<ntracks; itrack++)
    {
      AliFlowTrackSimple* track = dynamic_cast<AliFlowTrackSimple*>(fTrackCollection->At(itrack));
      if (!track) continue;
      AddTrack(static_cast<AliFlowTrackSimple*>(track->Clone()));
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::ResolutionPt(Double_t res)
{
  //smear pt of all tracks by gaussian with sigma=res
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->ResolutionPt(res);
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagSubeventsInEta( Double_t etaMinA,
                                            Double_t etaMaxA,
                                            Double_t etaMinB,
                                            Double_t etaMaxB )
{
  //Flag two subevents in given eta ranges
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    track->ResetSubEventTags();
    Double_t eta=track->Eta();
    if (eta >= etaMinA && eta <= etaMaxA) track->SetForSubevent(0);
    if (eta >= etaMinB && eta <= etaMaxB) track->SetForSubevent(1);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagSubeventsByCharge()
{
  //Flag two subevents in given eta ranges
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    track->ResetSubEventTags();
    Int_t charge=track->Charge();
    if (charge<0) track->SetForSubevent(0);
    if (charge>0) track->SetForSubevent(1);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV1( Double_t v1 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) {
      if((fUseExternalSymmetryPlanes)||(fUseGlauberMCSymmetryPlanes))
	track->AddV1(v1, fPsi1, fAfterBurnerPrecision);
      else
	track->AddV1(v1, fMCReactionPlaneAngle, fAfterBurnerPrecision);
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV2( Double_t v2 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) {
      if((fUseExternalSymmetryPlanes)||(fUseGlauberMCSymmetryPlanes))
	track->AddV2(v2, fPsi2, fAfterBurnerPrecision);
      else
	track->AddV2(v2, fMCReactionPlaneAngle, fAfterBurnerPrecision);
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV3( Double_t v3 )
{
  //add v3 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if(track) {
      if((fUseExternalSymmetryPlanes)||(fUseGlauberMCSymmetryPlanes))
	track->AddV3(v3, fPsi3, fAfterBurnerPrecision);
      else
	track->AddV3(v3, fMCReactionPlaneAngle, fAfterBurnerPrecision);
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV4( Double_t v4 )
{
  //add v4 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if(track) {
      if((fUseExternalSymmetryPlanes)||(fUseGlauberMCSymmetryPlanes))
	track->AddV4(v4, fPsi4, fAfterBurnerPrecision);
      else
	track->AddV4(v4, fMCReactionPlaneAngle, fAfterBurnerPrecision);
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV5( Double_t v5 )
{
  //add v4 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if(track) {
      if((fUseExternalSymmetryPlanes)||(fUseGlauberMCSymmetryPlanes))
	track->AddV5(v5, fPsi5, fAfterBurnerPrecision);
      else
	track->AddV5(v5, fMCReactionPlaneAngle, fAfterBurnerPrecision);
    }
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddFlow( Double_t v1, Double_t v2, Double_t v3, Double_t v4, Double_t v5,
                                  Double_t rp1, Double_t rp2, Double_t rp3, Double_t rp4, Double_t rp5 )
{
  //add flow to all tracks wrt the reaction plane angle, for all harmonic separate angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddFlow(v1,v2,v3,v4,v5,rp1,rp2,rp3,rp4,rp5,fAfterBurnerPrecision);
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddFlow( Double_t v1, Double_t v2, Double_t v3, Double_t v4, Double_t v5 )
{
  //add flow to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddFlow(v1,v2,v3,v4,v5,fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV2( TF1* ptDepV2 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    Double_t v2 = ptDepV2->Eval(track->Pt());
    track->AddV2(v2, fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV2( TF2* ptEtaDepV2 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    Double_t v2 = ptEtaDepV2->Eval(track->Pt(), track->Eta());
    track->AddV2(v2, fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
  SetUserModified();
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagRP( const AliFlowTrackSimpleCuts* cuts )
{
  //tag tracks as reference particles (RPs)
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    Bool_t pass=cuts->PassesCuts(track);
    Bool_t rpTrack=track->InRPSelection();
    if (pass)
    {
      if (!rpTrack) fNumberOfPOIs[0]++; //only increase if not already tagged
    }
    else
    {
      if (rpTrack) fNumberOfPOIs[0]--; //only decrease if detagging
    }
    track->SetForRPSelection(pass);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagPOI( const AliFlowTrackSimpleCuts* cuts, Int_t poiType )
{
  //tag tracks as particles of interest (POIs)
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    Bool_t pass=cuts->PassesCuts(track);
    Bool_t poiTrack=track->InPOISelection();
    if (pass)
    {
      if (!poiTrack) fNumberOfPOIs[poiType]++; //only increase if not already tagged
    }
    else
    {
      if (poiTrack) fNumberOfPOIs[poiType]--; //only decrease if detagging
    }
    track->Tag(poiType,pass);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagTracks( const AliFlowTrackSimpleCuts* cutsRP, const AliFlowTrackSimpleCuts* cutsPOI)
{
    // simple interface to tagging poi's and rp's
    TagPOI(cutsRP, 0);
    TagPOI(cutsPOI, 1);
}
//_____________________________________________________________________________
void AliFlowEventSimple::DefineDeadZone( Double_t etaMin,
                                         Double_t etaMax,
                                         Double_t phiMin,
                                         Double_t phiMax )
{
  //mark tracks in given eta-phi region as dead
  //by resetting the flow bits
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    Double_t eta = track->Eta();
    Double_t phi = track->Phi();
    if (eta>etaMin && eta<etaMax && phi>phiMin && phi<phiMax)
    {
      if (track->InRPSelection()) {fNumberOfPOIs[0]--;}
      for (Int_t j=1; j<fNumberOfPOItypes; j++)
      {
        if (track->CheckTag(j)) {fNumberOfPOIs[j]--;}
      }
      track->ResetPOItype();
    }
  }
}

//_____________________________________________________________________________
Int_t AliFlowEventSimple::CleanUpDeadTracks()
{
  //remove tracks that have no flow tags set and cleanup the container
  //returns number of cleaned tracks
  Int_t ncleaned=0;
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    if (track->IsDead()) {delete track;track=NULL;ncleaned++;}
  }
  fTrackCollection->Compress(); //clean up empty slots
  fNumberOfTracks-=ncleaned; //update number of tracks
  delete [] fShuffledIndexes; fShuffledIndexes=NULL;
  return ncleaned;
}

//_____________________________________________________________________________
TF1* AliFlowEventSimple::SimplePtDepV2()
{
  //return a standard pt dependent v2 formula, user has to clean up!
  return new TF1("StandardPtDepV2","((x<1.0)*(0.05/1.0)*x+(x>=1.0)*0.05)");
}

//_____________________________________________________________________________
TF2* AliFlowEventSimple::SimplePtEtaDepV2()
{
    //returna standard pt and eta dependent v2 formula, user has to clean up!
    return new TF2 ("f","((x<1)*.1*x+(x>=1)*.1)*(1-0.2*TMath::Abs(y))");
}

//_____________________________________________________________________________
TF1* AliFlowEventSimple::SimplePtSpectrum()
{
  //return a standard pt spectrum, user has to clean up!
  return new TF1("StandardPtSpectrum","x*TMath::Exp(-pow(0.13957*0.13957+x*x,0.5)/0.4)",0.1,10.);
}

//_____________________________________________________________________________
void AliFlowEventSimple::ClearFast()
{
  //clear the counters without deleting allocated objects so they can be reused
  fReferenceMultiplicity = 0;
  fNumberOfTracks = 0;
  for (Int_t i=0; i<fNumberOfPOItypes; i++)
  {
    fNumberOfPOIs[i] = 0;
  }
  fMCReactionPlaneAngle = 0.0;
  fMCReactionPlaneAngleIsSet = kFALSE;
  fAfterBurnerPrecision = 0.001;
  fUserModified = kFALSE;
  delete [] fShuffledIndexes; fShuffledIndexes=NULL;
}
