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
#include "TProfile.h"
#include "TParameter.h"
#include "TBrowser.h"
#include "AliFlowVector.h"
#include "AliFlowTrackSimple.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimpleCuts.h"
#include "TRandom.h"

ClassImp(AliFlowEventSimple)

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple():
  fTrackCollection(NULL),
  fReferenceMultiplicity(0),
  fNumberOfTracks(0),
  fNumberOfRPs(0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL)
{
  cout << "AliFlowEventSimple: Default constructor to be used only by root for io" << endl;
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple(Int_t aLength):
  fTrackCollection(new TObjArray(aLength)),
  fReferenceMultiplicity(0),
  fNumberOfTracks(0),
  fNumberOfRPs(0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL)
{
  //constructor
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple( Int_t nParticles,
                                        TF1* ptDist,
                                        Double_t phiMin,
                                        Double_t phiMax,
                                        Double_t etaMin,
                                        Double_t etaMax):
  fTrackCollection(new TObjArray(nParticles)),
  fReferenceMultiplicity(nParticles),
  fNumberOfTracks(0),
  fNumberOfRPs(0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL)
{
  //ctor. generates nParticles random tracks with given Pt distribution
  //phi and eta are uniform
  Generate(nParticles,ptDist,phiMin,phiMax,etaMin,etaMax);
}

//-----------------------------------------------------------------------
AliFlowEventSimple::AliFlowEventSimple(const AliFlowEventSimple& anEvent):
  TObject(anEvent),
  fTrackCollection((TObjArray*)(anEvent.fTrackCollection)->Clone()),
  fReferenceMultiplicity(anEvent.fReferenceMultiplicity),
  fNumberOfTracks(anEvent.fNumberOfTracks),
  fNumberOfRPs(anEvent.fNumberOfRPs),
  fMCReactionPlaneAngle(anEvent.fMCReactionPlaneAngle),
  fMCReactionPlaneAngleIsSet(anEvent.fMCReactionPlaneAngleIsSet),
  fAfterBurnerPrecision(anEvent.fAfterBurnerPrecision),
  fNumberOfTracksWrap(anEvent.fNumberOfTracksWrap),
  fNumberOfRPsWrap(anEvent.fNumberOfRPsWrap),
  fMCReactionPlaneAngleWrap(anEvent.fMCReactionPlaneAngleWrap)
{
  //copy constructor
}

//-----------------------------------------------------------------------
AliFlowEventSimple& AliFlowEventSimple::operator=(const AliFlowEventSimple& anEvent)
{
  //assignment operator
  delete fTrackCollection;
  fTrackCollection = (TObjArray*)(anEvent.fTrackCollection)->Clone(); //deep copy
  fReferenceMultiplicity = anEvent.fReferenceMultiplicity;
  fNumberOfTracks = anEvent.fNumberOfTracks;
  fNumberOfRPs = anEvent.fNumberOfRPs;
  fMCReactionPlaneAngle = anEvent.fMCReactionPlaneAngle;
  fMCReactionPlaneAngleIsSet = anEvent.fMCReactionPlaneAngleIsSet;
  fAfterBurnerPrecision = anEvent.fAfterBurnerPrecision;
  fNumberOfTracksWrap = anEvent.fNumberOfTracksWrap;
  fNumberOfRPsWrap = anEvent.fNumberOfRPsWrap;
  fMCReactionPlaneAngleWrap=anEvent.fMCReactionPlaneAngleWrap;
  return *this;
}

//-----------------------------------------------------------------------
AliFlowEventSimple::~AliFlowEventSimple()
{
  //destructor
  if (fTrackCollection) fTrackCollection->Delete();
  delete fTrackCollection;
  if (fNumberOfTracksWrap) delete fNumberOfTracksWrap;
  if (fNumberOfRPsWrap) delete fNumberOfRPsWrap;
  if (fMCReactionPlaneAngleWrap) delete fMCReactionPlaneAngleWrap;
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
  for (Int_t i=0; i<nParticles; i++)
  {
    AddTrack(new AliFlowTrackSimple( gRandom->Uniform(phiMin,phiMax),
                                     gRandom->Uniform(etaMin,etaMax),
                                     ptDist->GetRandom()));
  }
}

//-----------------------------------------------------------------------
AliFlowTrackSimple* AliFlowEventSimple::GetTrack(Int_t i)
{
  //get track i from collection
  if (i>=fNumberOfTracks) return NULL;
  AliFlowTrackSimple* pTrack = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i)) ;
  return pTrack;
}

//-----------------------------------------------------------------------
void AliFlowEventSimple::AddTrack( AliFlowTrackSimple* track )
{
  //add a track
  fTrackCollection->AddLast(track);
  fNumberOfTracks++;
}

//-----------------------------------------------------------------------
AliFlowVector AliFlowEventSimple::GetQ(Int_t n, TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights)
{
  // calculate Q-vector in harmonic n without weights (default harmonic n=2)
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  AliFlowVector vQ;
  vQ.Set(0.,0.);

  Int_t iOrder = n;
  Double_t iUsedTracks = 0;
  Double_t dPhi=0.;
  Double_t dPt=0.;
  Double_t dEta=0.;

  AliFlowTrackSimple* pTrack = NULL;

  Int_t nBinsPhi=0;
  Double_t dBinWidthPt=0.;
  Double_t dPtMin=0.;
  Double_t dBinWidthEta=0.;
  Double_t dEtaMin=0.;

  Double_t wPhi=1.; // weight Phi
  Double_t wPt=1.;  // weight Pt
  Double_t wEta=1.; // weight Eta

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
        dQX += wPhi*wPt*wEta*TMath::Cos(iOrder*dPhi);
        dQY += wPhi*wPt*wEta*TMath::Sin(iOrder*dPhi);

        // weighted multiplicity:
        iUsedTracks+=wPhi*wPt*wEta;

      } // end of if (pTrack->InRPSelection())
    } // end of if (pTrack)
    else
    {
      cerr << "no particle!!!"<<endl;
    }
  } // loop over particles

  vQ.Set(dQX,dQY);
  vQ.SetMult(iUsedTracks);

  return vQ;

}

//-----------------------------------------------------------------------
void AliFlowEventSimple::Get2Qsub(AliFlowVector* Qarray, Int_t n, TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights)
{

  // calculate Q-vector in harmonic n without weights (default harmonic n=2)
  Double_t dQX = 0.;
  Double_t dQY = 0.;

  Int_t iOrder = n;
  Double_t iUsedTracks = 0;
  Double_t dPhi = 0.;
  Double_t dPt  = 0.;
  Double_t dEta = 0.;

  AliFlowTrackSimple* pTrack = NULL;

  Int_t    iNbinsPhi   = 0;
  Double_t dBinWidthPt = 0.;
  Double_t dPtMin      = 0.;
  Double_t dBinWidthEta= 0.;
  Double_t dEtaMin     = 0.;

  Double_t dWphi = 1.;  // weight Phi
  Double_t dWpt  = 1.;  // weight Pt
  Double_t dWeta = 1.;  // weight Eta

  TH1F* phiWeights = NULL;
  TH1D* ptWeights  = NULL;
  TH1D* etaWeights = NULL;

  if(weightsList)
  {
    if(usePhiWeights)
    {
      phiWeights = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights"));
      if(phiWeights)
      {
        iNbinsPhi = phiWeights->GetNbinsX();
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
      if(pTrack)
      {
        if(pTrack->InRPSelection())
        {
          if (pTrack->InSubevent(s))
          {
            dPhi = pTrack->Phi();
            dPt  = pTrack->Pt();
            dEta = pTrack->Eta();

            // determine Phi weight: (to be improved, I should here only access it + the treatment of gaps in the if statement)
            if(phiWeights && iNbinsPhi)
            {
              dWphi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*iNbinsPhi/TMath::TwoPi())));
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
            dQX += dWphi*dWpt*dWeta*TMath::Cos(iOrder*dPhi);
            dQY += dWphi*dWpt*dWeta*TMath::Sin(iOrder*dPhi);

            // weighted multiplicity:
            iUsedTracks+=dWphi*dWpt*dWeta;

          } // end of subevent
        } // end of if (pTrack->InRPSelection())
      } // end of if (pTrack)
      else
      {
        cerr << "no particle!!!"<<endl;
      }
    } // loop over particles
    Qarray[s].Set(dQX,dQY);
    Qarray[s].SetMult(iUsedTracks);
    //reset
    iUsedTracks = 0;
    dQX = 0.;
    dQY = 0.;
  }

}


//-----------------------------------------------------------------------
void AliFlowEventSimple::Print(Option_t *option) const
{
  //   -*-*-*-*-*Print some global quantities for this histogram collection class *-*-*-*-*-*-*-*
  //             ===============================================
  //   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
  printf( "Class.Print Name = %s, #tracks= %d, Number of RPs= %d, MC EventPlaneAngle= %f\n",
          GetName(),fNumberOfTracks, fNumberOfRPs, fMCReactionPlaneAngle );

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
  if (!b) return;
  if (!fNumberOfTracksWrap)
  {
    fNumberOfTracksWrap = new TParameter<int>("fNumberOfTracks", fNumberOfTracks);
    b->Add(fNumberOfTracksWrap);
  }
  if (!fNumberOfRPsWrap)
  {
    fNumberOfRPsWrap = new TParameter<int>("fNumberOfRPs", fNumberOfRPs);
    b->Add(fNumberOfRPsWrap);
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
  fNumberOfRPs(0),
  fMCReactionPlaneAngle(0.),
  fMCReactionPlaneAngleIsSet(kFALSE),
  fAfterBurnerPrecision(0.001),
  fNumberOfTracksWrap(NULL),
  fNumberOfRPsWrap(NULL),
  fMCReactionPlaneAngleWrap(NULL)
{
  //constructor, fills the event from a TTree of kinematic.root files
  //applies RP and POI cuts, tags the tracks

  Int_t numberOfInputTracks = inputTree->GetEntries() ;
  fTrackCollection = new TObjArray(numberOfInputTracks/2);

  TParticle* pParticle = new TParticle();
  inputTree->SetBranchAddress("Particles",&pParticle);

  Int_t iSelParticlesPOI = 0;

  for (Int_t i=0; i<numberOfInputTracks; i++)
  {
    inputTree->GetEntry(i);   //get input particle
    
    if (!pParticle) continue; //no particle
    if (!pParticle->IsPrimary()) continue;

    Bool_t rpOK = rpCuts->PassesCuts(pParticle);
    Bool_t poiOK = poiCuts->PassesCuts(pParticle);
    
    if (rpOK || poiOK)
    {
      AliFlowTrackSimple* pTrack = new AliFlowTrackSimple(pParticle);

      //marking the particles used for int. flow:
      if(rpOK)
      {
        pTrack->SetForRPSelection(kTRUE);
        fNumberOfRPs++;
      }
      //marking the particles used for diff. flow:
      if(poiOK)
      {
        pTrack->SetForPOISelection(kTRUE);
        iSelParticlesPOI++;
      }
      //adding a particles which were used either for int. or diff. flow to the list
      AddTrack(pTrack);
    }
  }//for i
  delete pParticle;
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
    Double_t eta=track->Eta();
    if (eta >= etaMinA && eta <= etaMaxA) track->SetForSubevent(0);
    if (eta >= etaMinB && eta <= etaMaxB) track->SetForSubevent(1);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV1( Double_t v1 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddV1(v1, fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV2( Double_t v2 )
{
  //add v2 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddV2(v2, fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddV4( Double_t v4 )
{
  //add v4 to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddV4(v4, fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::AddFlow( Double_t v1, Double_t v2, Double_t v4 )
{
  //add flow to all tracks wrt the reaction plane angle
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (track) track->AddFlow(v1,v2,v4,fMCReactionPlaneAngle, fAfterBurnerPrecision);
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagRP( AliFlowTrackSimpleCuts* cuts )
{
  //tag tracks as reference particles (RPs)
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    if (cuts->PassesCuts(track)) 
    {
      track->SetForRPSelection();
      fNumberOfRPs++;
    }
  }
}

//_____________________________________________________________________________
void AliFlowEventSimple::TagPOI( AliFlowTrackSimpleCuts* cuts )
{
  //tag tracks as particles of interest (POIs)
  for (Int_t i=0; i<fNumberOfTracks; i++)
  {
    AliFlowTrackSimple* track = static_cast<AliFlowTrackSimple*>(fTrackCollection->At(i));
    if (!track) continue;
    if (cuts->PassesCuts(track)) track->SetForPOISelection();
  }
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
      if (track->InRPSelection()) fNumberOfRPs--;
      track->ResetFlowTags();
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
    if (track->IsDead()) {delete track;ncleaned++;}
  }
  fTrackCollection->Compress(); //clean up empty slots
  return ncleaned;
}
