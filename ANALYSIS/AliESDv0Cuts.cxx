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

/* $Id: AliESDv0Cuts.cxx 24534 2008-03-16 22:22:11Z fca $ */

#include "AliESDv0Cuts.h"

#include <AliESDVertex.h>
#include <AliESDtrack.h>
#include <AliESDv0.h>
#include <AliESD.h>
#include <AliESDEvent.h>
#include <AliLog.h>

#include <TTree.h>
#include <TCanvas.h>
#include <TDirectory.h>

//____________________________________________________________________
ClassImp(AliESDv0Cuts)

// Cut names
const Char_t* AliESDv0Cuts::fgkCutNames[kNCuts] = {
 "dca positive to pvtx",
 "dca negative to pvtx",
 "#Chi^{2}",
 "dca v0 daughters",
 "min decay radius",
 "max decay radius",
 "cosine pointing angle",
 "on-the-fly status",
 "dca v0 to pvtx",
 "p",
 "p_{T}",
 "p_{x}",
 "p_{y}",
 "p_{z}"
};

//____________________________________________________________________
AliESDv0Cuts::AliESDv0Cuts(const Char_t* name, const Char_t* title) : AliAnalysisCuts(name,title),
  fCutMinDcaPosToVertex(0),
  fCutMinDcaNegToVertex(0),
  fCutMaxChi2(40),
  fCutMaxDcaV0Daughters(0),
  fCutMinRadius(0),
  fCutMaxRadius(0),
  fCutMinCosinePointingAngle(0),
  fCutRequireOnFlyStatus(0),
  fCutMaxDcaV0ToVertex(0),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fHistogramsOn(0),
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // constructor
  //

  Init();

  //##############################################################################
  // setting default cuts
  SetMinDcaPosToVertex();
  SetMinDcaNegToVertex();
  SetMaxChi2();
  SetMaxDcaV0Daughters();
  SetMinRadius();
  SetMaxRadius();
  SetMinCosinePointingAngle();
  SetRequireOnFlyStatus();
  SetMaxDcaV0ToVertex();
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();

  SetHistogramsOn();
}

//_____________________________________________________________________________
AliESDv0Cuts::AliESDv0Cuts(const AliESDv0Cuts &c) : AliAnalysisCuts(c),
  fCutMinDcaPosToVertex(0),
  fCutMinDcaNegToVertex(0),
  fCutMaxChi2(0),
  fCutMaxDcaV0Daughters(0),
  fCutMinRadius(0),
  fCutMaxRadius(0),
  fCutMinCosinePointingAngle(0),
  fCutRequireOnFlyStatus(0),
  fCutMaxDcaV0ToVertex(0),
  fPMin(0),
  fPMax(0),
  fPtMin(0),
  fPtMax(0),
  fPxMin(0),
  fPxMax(0),
  fPyMin(0),
  fPyMax(0),
  fPzMin(0),
  fPzMax(0),
  fHistogramsOn(0),
  fhCutStatistics(0),         
  fhCutCorrelation(0)
{
  //
  // copy constructor
  //

  ((AliESDv0Cuts &) c).Copy(*this);
}

AliESDv0Cuts::~AliESDv0Cuts()
{
  //
  // destructor
  //

  for (Int_t i=0; i<2; i++) {
    
    if (fhDcaPosToVertex[i])
      delete fhDcaPosToVertex[i];
    if (fhDcaNegToVertex[i])
      delete fhDcaNegToVertex[i];
    if (fhChi2[i])
      delete fhChi2[i];
    if (fhDcaV0Daughters[i])
      delete fhDcaV0Daughters[i];
    if (fhRadius[i])
      delete fhRadius[i];             
    if (fhCosinePointingAngle[i])
      delete fhCosinePointingAngle[i];
    if (fhOnFlyStatus[i])
    delete fhOnFlyStatus[i];
    if (fhDcaV0ToVertex[i])
      delete fhDcaV0ToVertex[i];
    if (fhPt[i])
      delete fhPt[i];
  }

  if (fhCutStatistics)
    delete fhCutStatistics;             
  if (fhCutCorrelation)
    delete fhCutCorrelation;            
}

void AliESDv0Cuts::Init()
{
  //
  // sets everything to zero
  //
  fCutMinDcaPosToVertex      = 0;
  fCutMinDcaNegToVertex      = 0;
  fCutMaxChi2                = 0;
  fCutMaxDcaV0Daughters      = 0;
  fCutMinRadius              = 0;
  fCutMaxRadius              = 0;
  fCutMinCosinePointingAngle = 0;
  fCutRequireOnFlyStatus     = 0;
  fCutMaxDcaV0ToVertex       = 0;

  fPMin  = 0;
  fPMax  = 0;
  fPtMin = 0;
  fPtMax = 0;
  fPxMin = 0;
  fPxMax = 0;
  fPyMin = 0;
  fPyMax = 0;
  fPzMin = 0;
  fPzMax = 0;

  fHistogramsOn = kFALSE;

  for (Int_t i=0; i<2; ++i)
  {
    fhDcaPosToVertex[i]      = 0;
    fhDcaNegToVertex[i]      = 0;
    fhChi2[i]                = 0;
    fhDcaV0Daughters[i]      = 0;
    fhRadius[i]              = 0;
    fhCosinePointingAngle[i] = 0;
    fhOnFlyStatus[i]         = 0;
    fhDcaV0ToVertex[i]       = 0;
    
    fhPt[i]                  = 0;
  }
  fhCutStatistics = 0;
  fhCutCorrelation = 0;
}

//_____________________________________________________________________________
AliESDv0Cuts &AliESDv0Cuts::operator=(const AliESDv0Cuts &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliESDv0Cuts &) c).Copy(*this);
  return *this;
}

//_____________________________________________________________________________
void AliESDv0Cuts::Copy(TObject &c) const
{
  //
  // Copy function
  //

  AliESDv0Cuts& target = (AliESDv0Cuts &) c;

  target.Init();

  target.fCutMinDcaPosToVertex      = fCutMinDcaPosToVertex;
  target.fCutMinDcaNegToVertex      = fCutMinDcaNegToVertex;
  target.fCutMaxChi2                = fCutMaxChi2;
  target.fCutMaxDcaV0Daughters      = fCutMaxDcaV0Daughters;
  target.fCutMinRadius              = fCutMinRadius;
  target.fCutMaxRadius              = fCutMaxRadius;
  target.fCutMinCosinePointingAngle = fCutMinCosinePointingAngle;
  target.fCutRequireOnFlyStatus     = fCutRequireOnFlyStatus;
  target.fCutMaxDcaV0ToVertex       = fCutMaxDcaV0ToVertex;

  target.fPMin  = fPMin;
  target.fPMax  = fPMax;
  target.fPtMin = fPtMin;
  target.fPtMax = fPtMax;
  target.fPxMin = fPxMin;
  target.fPxMax = fPxMax;
  target.fPyMin = fPyMin;
  target.fPyMax = fPyMax;
  target.fPzMin = fPzMin;
  target.fPzMax = fPzMax;

  target.fHistogramsOn = fHistogramsOn;

  for (Int_t i=0; i<2; ++i)
  {
    if (fhDcaPosToVertex[i]) target.fhDcaPosToVertex[i] = (TH1F*) fhDcaPosToVertex[i]->Clone();
    if (fhDcaNegToVertex[i]) target.fhDcaNegToVertex[i] = (TH1F*) fhDcaNegToVertex[i]->Clone();
    if (fhChi2[i]) target.fhChi2[i] = (TH1F*) fhChi2[i]->Clone();
    if (fhDcaV0Daughters[i]) target.fhDcaV0Daughters[i] = (TH1F*) fhDcaV0Daughters[i]->Clone();
    if (fhRadius[i]) target.fhRadius[i] = (TH1F*) fhRadius[i]->Clone();
    if (fhCosinePointingAngle[i]) target.fhCosinePointingAngle[i] = (TH1F*) fhCosinePointingAngle[i]->Clone();
    if (fhOnFlyStatus[i]) target.fhOnFlyStatus[i] = (TH1F*) fhOnFlyStatus[i]->Clone();
    if (fhDcaV0ToVertex[i]) target.fhDcaV0ToVertex[i] = (TH1F*) fhDcaV0ToVertex[i]->Clone();
    
    if (fhPt[i]) target.fhPt[i] = (TH1F*) fhPt[i]->Clone();
  }
  if (fhCutStatistics) target.fhCutStatistics = (TH1F*) fhCutStatistics->Clone();
  if (fhCutCorrelation) target.fhCutCorrelation = (TH2F*) fhCutCorrelation->Clone();

  TNamed::Copy(c);
}

//_____________________________________________________________________________
Long64_t AliESDv0Cuts::Merge(TCollection* list) {
  // Merge a list of AliESDv0Cuts objects with this (needed for PROOF)
  // Returns the number of merged objects (including this)

  if (!list)
    return 0;
  
  if (list->IsEmpty())
    return 1;

  if (!fHistogramsOn)
    return 0;

  TIterator* iter = list->MakeIterator();
  TObject* obj;


  // collection of measured and generated histograms
  Int_t count = 0;
  while ((obj = iter->Next())) {

    AliESDv0Cuts* entry = dynamic_cast<AliESDv0Cuts*>(obj);
    if (entry == 0)
      continue;

    if (!entry->fHistogramsOn)
      continue;
    
    for (Int_t i=0; i<2; i++) {
      
      fhDcaPosToVertex[i]     ->Add(entry->fhDcaPosToVertex[i]     );
      fhDcaNegToVertex[i]     ->Add(entry->fhDcaNegToVertex[i]     );
      fhChi2[i]               ->Add(entry->fhChi2[i]               );
      fhDcaV0Daughters[i]     ->Add(entry->fhDcaV0Daughters[i]     );
      fhRadius[i]             ->Add(entry->fhRadius[i]             );
      fhCosinePointingAngle[i]->Add(entry->fhCosinePointingAngle[i]);
      fhOnFlyStatus[i]        ->Add(entry->fhOnFlyStatus[i]        );
      fhDcaV0ToVertex[i]      ->Add(entry->fhDcaV0ToVertex[i]      );

      fhPt[i]                ->Add(entry->fhPt[i]);
    }      
    fhCutStatistics  ->Add(entry->fhCutStatistics);        
    fhCutCorrelation ->Add(entry->fhCutCorrelation);      

    count++;
  }

  return count+1;
}

void AliESDv0Cuts::EnableNeededBranches(TTree* tree)
{
  // enables the branches needed by AcceptV0, for a list see comment of AcceptV0

  tree->SetBranchStatus("fV0s.fDcaV0Daughters", 1);
  tree->SetBranchStatus("fV0s.fChi2V0", 1);
  tree->SetBranchStatus("fV0s.fPos*", 1);
  tree->SetBranchStatus("fV0s.fNmom*", 1);
  tree->SetBranchStatus("fV0s.fPmom*", 1);
  tree->SetBranchStatus("fV0s.fRr", 1);
  tree->SetBranchStatus("fV0s.fPointAngle*", 1);
  tree->SetBranchStatus("fV0s.fOnFlyStatus", 1);
}

//____________________________________________________________________
Bool_t
AliESDv0Cuts::IsSelected(TList* listObj) {
// Selection cuts
  if(listObj->GetSize()!=4) return kFALSE;
  AliESDv0           *esdV0     = (AliESDv0*)listObj->At(0);
  AliESDtrack        *trackPos  = (AliESDtrack*)listObj->At(1);
  AliESDtrack        *trackNeg  = (AliESDtrack*)listObj->At(2);
  const AliESDVertex *esdVertex = (AliESDVertex*)listObj->At(3);
  return AcceptV0(esdV0,trackPos,trackNeg,esdVertex);
}

//____________________________________________________________________
Bool_t
AliESDv0Cuts::AcceptV0(AliESDv0* const esdV0, AliESDtrack* const trackPos, AliESDtrack* const trackNeg, const AliESDVertex* esdVertex) {
  // 
  // figure out if the v0s survives all the v0 cuts defined
  //
  // the different quality parameter and kinematic values are first
  // retrieved from the v0. then it is found out what cuts the
  // v0 did not survive and finally the cuts are imposed.

  // this function needs the following branches... but this is not enough
  // fV0s.fDcaV0Daughters
  // fV0s.fChi2V0
  // fV0s.fPos*
  // fV0s.fNmom*
  // fV0s.fPmom*
  // fV0s.fRr
  // fV0s.fPointAngle
  // fV0s.fOnFlyStatus

  Float_t  dcaPosToVertex = 0, dcaNegToVertex = 0;
  Float_t  tdcaPosToVertex[2]={999,999};
  Float_t  tdcaNegToVertex[2]={999,999};

  if (trackPos) trackPos->GetImpactParameters(tdcaPosToVertex[0],tdcaPosToVertex[1]);
  if (trackNeg) trackNeg->GetImpactParameters(tdcaNegToVertex[0],tdcaNegToVertex[1]);

  dcaPosToVertex = TMath::Sqrt(tdcaPosToVertex[0]*tdcaPosToVertex[0]+tdcaPosToVertex[1]*tdcaPosToVertex[1]);
  dcaNegToVertex = TMath::Sqrt(tdcaNegToVertex[0]*tdcaNegToVertex[0]+tdcaNegToVertex[1]*tdcaNegToVertex[1]);

  UInt_t  status = esdV0->GetOnFlyStatus();
  Float_t chi2  = esdV0->GetChi2V0();

  Double_t dcaV0Daughters = esdV0->GetDcaV0Daughters();

  Double_t vtxPosition[3]; esdVertex->GetXYZ(vtxPosition);
  Double_t dcaV0ToVertex  = esdV0->GetD(vtxPosition[0],vtxPosition[1],vtxPosition[2]);

  Double_t v0Position[3];
  esdV0->GetXYZ(v0Position[0],v0Position[1],v0Position[2]);
  Double_t radius = TMath::Sqrt(TMath::Power(v0Position[0],2) + TMath::Power(v0Position[1],2));
  Double_t v0CosinePointingAngle = esdV0->GetV0CosineOfPointingAngle();

  // getting the kinematic variables of the v0
  Double_t p[3];
  esdV0->GetPxPyPz(p[0],p[1],p[2]);
  Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
  Float_t pt       = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2));

  //########################################################################
  // cut the v0?
  
  Bool_t cuts[kNCuts];
  for (Int_t i=0; i<kNCuts; i++) cuts[i]=kFALSE;
  
  // v0 quality cuts
  if (dcaPosToVertex < fCutMinDcaPosToVertex)
    cuts[0]=kTRUE;
  if (dcaNegToVertex < fCutMinDcaNegToVertex) 
    cuts[1]=kTRUE;
  if (chi2 > fCutMaxChi2) 
    cuts[2]=kTRUE; 
  if (dcaV0Daughters > fCutMaxDcaV0Daughters) 
    cuts[3]=kTRUE;  
  if (radius  < fCutMinRadius) 
    cuts[4]=kTRUE;  
  if (radius  > fCutMaxRadius) 
    cuts[5]=kTRUE;  
  if (v0CosinePointingAngle < fCutMinCosinePointingAngle)
    cuts[6]=kTRUE;  
  if (fCutRequireOnFlyStatus && !status)
    cuts[7]=kTRUE;  
  if (dcaV0ToVertex > fCutMaxDcaV0ToVertex)
    cuts[8] = kTRUE;

  // v0 kinematics cut
  if((momentum < fPMin) || (momentum > fPMax)) 
    cuts[9]=kTRUE;
  if((pt < fPtMin) || (pt > fPtMax)) 
    cuts[10] = kTRUE;
  if((p[0] < fPxMin) || (p[0] > fPxMax)) 
    cuts[11] = kTRUE;
  if((p[1] < fPyMin) || (p[1] > fPyMax)) 
    cuts[12] = kTRUE;
  if((p[2] < fPzMin) || (p[2] > fPzMax))
    cuts[13] = kTRUE;

  Bool_t cut=kFALSE;
  for (Int_t i=0; i<kNCuts; i++) 
    if (cuts[i]) cut = kTRUE;
  
  //########################################################################
  // filling histograms
  if (fHistogramsOn) {
    fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n v0s")));
    
    if (cut)
      fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n cut v0s")));
    
    for (Int_t i=0; i<kNCuts; i++) {
      if (cuts[i])
 	fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin(fgkCutNames[i])));
      
      for (Int_t j=i; j<kNCuts; j++) {
 	if (cuts[i] && cuts[j]) {
 	  Float_t x = fhCutCorrelation->GetXaxis()->GetBinCenter(fhCutCorrelation->GetXaxis()->FindBin(fgkCutNames[i]));
 	  Float_t y = fhCutCorrelation->GetYaxis()->GetBinCenter(fhCutCorrelation->GetYaxis()->FindBin(fgkCutNames[j]));
 	  fhCutCorrelation->Fill(x,y);
 	}
      }
    }
    
    fhDcaPosToVertex[0]->Fill(dcaPosToVertex);
    fhDcaNegToVertex[0]->Fill(dcaNegToVertex);
    fhChi2[0]->Fill(chi2);
    fhDcaV0Daughters[0]->Fill(dcaV0Daughters);
    fhRadius[0]->Fill(radius);
    fhCosinePointingAngle[0]->Fill(v0CosinePointingAngle);
    fhOnFlyStatus[0]->Fill(status);
    fhDcaV0ToVertex[0]->Fill(dcaV0ToVertex);
    
    fhPt[0]->Fill(pt);
  }

  //########################################################################
  // cut the v0!
  if (cut) return kFALSE;

  //########################################################################
  // filling histograms after cut
  if (fHistogramsOn) {
    fhDcaPosToVertex[1]->Fill(dcaPosToVertex);
    fhDcaNegToVertex[1]->Fill(dcaNegToVertex);
    fhChi2[1]->Fill(chi2);
    fhDcaV0Daughters[1]->Fill(dcaV0Daughters);
    fhRadius[1]->Fill(radius);
    fhCosinePointingAngle[1]->Fill(v0CosinePointingAngle);
    fhOnFlyStatus[1]->Fill(status);
    fhDcaV0ToVertex[1]->Fill(dcaV0ToVertex);
    
    fhPt[1]->Fill(pt);
  }

  return kTRUE;
}

//____________________________________________________________________
TObjArray* AliESDv0Cuts::GetAcceptedV0s(const AliESD* esd)
{
  //
  // returns an array of all v0s that pass the cuts
  //

  TObjArray* acceptedV0s = new TObjArray();
  //  const AliESDVertex *spdVertex = esd->GetVertex();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;

  // loop over esd v0s
  for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = esd->GetV0(iV0);

    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *trackPos = esd->GetTrack(lIndexTrackPos);
    AliESDtrack *trackNeg = esd->GetTrack(lIndexTrackNeg);

    if (AcceptV0(v0,trackPos,trackNeg,primaryVertex))
      acceptedV0s->Add(v0);
  }

  return acceptedV0s;
}

//____________________________________________________________________
Int_t AliESDv0Cuts::CountAcceptedV0s(const AliESD* esd)
{
  //
  // returns an the number of v0s that pass the cuts
  //

  Int_t count = 0;
  //  const AliESDVertex *spdVertex = esd->GetVertex();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;

  // loop over esd v0s
  for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = esd->GetV0(iV0);

    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *trackPos = esd->GetTrack(lIndexTrackPos);
    AliESDtrack *trackNeg = esd->GetTrack(lIndexTrackNeg);

    if (AcceptV0(v0,trackPos,trackNeg,primaryVertex))
      count++;
  }

  return count;
}

//____________________________________________________________________
TObjArray* AliESDv0Cuts::GetAcceptedV0s(const AliESDEvent* esd)
{
  //
  // returns an array of all v0s that pass the cuts
  //

  TObjArray* acceptedV0s = new TObjArray();
  //  const AliESDVertex *spdVertex = esd->GetVertex();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;

  // loop over esd v0s
  for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = esd->GetV0(iV0);

    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *trackPos = esd->GetTrack(lIndexTrackPos);
    AliESDtrack *trackNeg = esd->GetTrack(lIndexTrackNeg);

    if (AcceptV0(v0,trackPos,trackNeg,primaryVertex))
      acceptedV0s->Add(v0);
  }

  return acceptedV0s;
}

//____________________________________________________________________
Int_t AliESDv0Cuts::CountAcceptedV0s(const AliESDEvent* esd)
{
  //
  // returns an the number of v0s that pass the cuts
  //

  Int_t count = 0;
  //  const AliESDVertex *spdVertex = esd->GetVertex();
  const AliESDVertex *primaryVertex = esd->GetPrimaryVertex();
  Int_t    lIndexTrackPos       = 0, lIndexTrackNeg       = 0;

  // loop over esd v0s
  for (Int_t iV0 = 0; iV0 < esd->GetNumberOfV0s(); iV0++) {
    AliESDv0* v0 = esd->GetV0(iV0);

    lIndexTrackPos = TMath::Abs(v0->GetPindex());
    lIndexTrackNeg = TMath::Abs(v0->GetNindex());
    AliESDtrack *trackPos = esd->GetTrack(lIndexTrackPos);
    AliESDtrack *trackNeg = esd->GetTrack(lIndexTrackNeg);

    if (AcceptV0(v0,trackPos,trackNeg,primaryVertex))
      count++;
  }

  return count;
}

//____________________________________________________________________
 void AliESDv0Cuts::DefineHistograms(Int_t color) {
   // 
   // diagnostics histograms are defined
   // 

   fHistogramsOn=kTRUE;
   
   Bool_t oldStatus = TH1::AddDirectoryStatus();
   TH1::AddDirectory(kFALSE);
   
   //###################################################################################
   // defining histograms

   fhCutStatistics = new TH1F("cut_statistics","cut statistics",kNCuts+4,-0.5,kNCuts+3.5);

   fhCutStatistics->GetXaxis()->SetBinLabel(1,"n v0s");
   fhCutStatistics->GetXaxis()->SetBinLabel(2,"n cut v0s");

   fhCutCorrelation = new TH2F("cut_correlation","cut correlation",kNCuts,-0.5,kNCuts-0.5,kNCuts,-0.5,kNCuts-0.5);;
  
   for (Int_t i=0; i<kNCuts; i++) {
     fhCutStatistics->GetXaxis()->SetBinLabel(i+4,fgkCutNames[i]);
     fhCutCorrelation->GetXaxis()->SetBinLabel(i+1,fgkCutNames[i]);
     fhCutCorrelation->GetYaxis()->SetBinLabel(i+1,fgkCutNames[i]);
   } 

  fhCutStatistics  ->SetLineColor(color);
  fhCutCorrelation ->SetLineColor(color);
  fhCutStatistics  ->SetLineWidth(2);
  fhCutCorrelation ->SetLineWidth(2);

  Char_t str[256];
  for (Int_t i=0; i<2; i++) {
    if (i==0) snprintf(str,256, " ");
    else snprintf(str,256, "_cut");

    fhDcaPosToVertex[i]      = new TH1F(Form("dcaPosToVertex%s",str),"",120,0,3);
    fhDcaNegToVertex[i]      = new TH1F(Form("dcaNegToVertex%s",str),"",120,0,3);
    fhChi2[i]                = new TH1F(Form("chi2%s",str),"",50,0,50);
    fhDcaV0Daughters[i]      = new  TH1F(Form("dcaV0Daughters%s",str),"",200,0,5);
    fhRadius[i]              = new  TH1F(Form("decayRadius%s",str),"",300,0,150);
    fhCosinePointingAngle[i] = new  TH1F(Form("cosinePointingAngle%s",str),"",100,-1,1);
    fhOnFlyStatus[i]         = new  TH1F(Form("onflystatus%s",str),"",5,0,5);
    fhDcaV0ToVertex[i]       = new  TH1F(Form("dcaV0ToVertex%s",str),"",100,0,5);

    fhPt[i]                  = new TH1F(Form("pt%s",str)     ,"p_{T} distribution;p_{T} (GeV/c)",500,0.0,100.0);
    
    fhDcaPosToVertex[i]->SetTitle("Dca of positive daughter to parent vertex");
    fhDcaNegToVertex[i]->SetTitle("Dca of negative daughter to parent vertex");
    fhChi2[i]->SetTitle("#Chi^{2} for v0");
    fhDcaV0Daughters[i]->SetTitle("Dca between daughter tracks");
    fhRadius[i]->SetTitle("Decay radius of the v0");
    fhCosinePointingAngle[i]->SetTitle("Cosine of the Pointing Angle");
    fhOnFlyStatus[i]->SetTitle("On-the-Fly Status");
    fhDcaV0ToVertex[i]->SetTitle("Dca of v0 to parent vertex");

    fhDcaPosToVertex[i]->SetLineColor(color);      fhDcaPosToVertex[i]->SetLineWidth(2);
    fhDcaNegToVertex[i]->SetLineColor(color);      fhDcaNegToVertex[i]->SetLineWidth(2);
    fhChi2[i]->SetLineColor(color);                fhChi2[i]->SetLineWidth(2);
    fhDcaV0Daughters[i]->SetLineColor(color);      fhDcaV0Daughters[i]->SetLineWidth(2);
    fhRadius[i]->SetLineColor(color);              fhRadius[i]->SetLineWidth(2);
    fhCosinePointingAngle[i]->SetLineColor(color); fhCosinePointingAngle[i]->SetLineWidth(2);
    fhOnFlyStatus[i]->SetLineColor(color);         fhOnFlyStatus[i]->SetLineWidth(2);
    fhDcaV0ToVertex[i]->SetLineColor(color);       fhDcaV0ToVertex[i]->SetLineWidth(2); 
  }

  TH1::AddDirectory(oldStatus);
}

//____________________________________________________________________
Bool_t AliESDv0Cuts::LoadHistograms(const Char_t* dir)
{
  //
  // loads the histograms from a file
  // if dir is empty a directory with the name of this object is taken (like in SaveHistogram)
  //

  if (!dir)
    dir = GetName();

  if (!gDirectory->cd(dir))
    return kFALSE;

  fhCutStatistics = dynamic_cast<TH1F*> (gDirectory->Get("cut_statistics"));
  fhCutCorrelation = dynamic_cast<TH2F*> (gDirectory->Get("cut_correlation"));

  Char_t str[5];
  for (Int_t i=0; i<2; i++) {
    if (i==0)
    {
      gDirectory->cd("before_cuts");
      str[0] = 0;
    }
    else
    {
      gDirectory->cd("after_cuts");
      snprintf(str,5, "_cut");
    }

    fhDcaPosToVertex[i]      = dynamic_cast<TH1F*> (gDirectory->Get(Form("dcaPosToVertex%s",str)     ));
    fhDcaNegToVertex[i]      = dynamic_cast<TH1F*> (gDirectory->Get(Form("dcaNegToVertex%s",str)     ));
    fhChi2[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("chi2%s",str)));
    fhDcaV0Daughters[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("dcaV0Daughters%s",str)));
    fhRadius[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("decayRadius%s",str)));
    fhCosinePointingAngle[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("cosinepointingangle%s",str)));
    fhOnFlyStatus[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("onflystatus%s",str)));
    fhDcaV0ToVertex[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("dcaV0ToVertex%s",str)));

    fhPt[i] = dynamic_cast<TH1F*> (gDirectory->Get(Form("pt%s",str)));

    gDirectory->cd("../");
  }

  gDirectory->cd("..");

  return kTRUE;
}

//____________________________________________________________________
void AliESDv0Cuts::SaveHistograms(const Char_t* dir) {
  //
  // saves the histograms in a directory (dir)
  //

  if (!fHistogramsOn) {
    AliDebug(0, "Histograms not on - cannot save histograms!!!");
    return;
  }

  if (!dir)
    dir = GetName();

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");
 
  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t i=0; i<2; i++) {
    if (i==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");

    fhDcaPosToVertex[i]      ->Write();
    fhDcaNegToVertex[i]      ->Write();
    fhChi2[i]                ->Write();
    fhDcaV0Daughters[i]      ->Write();
    fhRadius[i]              ->Write();
    fhCosinePointingAngle[i] ->Write();
    fhOnFlyStatus[i]         ->Write();
    fhDcaV0ToVertex[i]       ->Write();

    fhPt[i]                  ->Write();
    
    gDirectory->cd("../");
  }

  gDirectory->cd("../");
}

//____________________________________________________________________
void AliESDv0Cuts::DrawHistograms()
{
  // draws some histograms

  TCanvas* canvas1 = new TCanvas(Form("%s_1", GetName()), "V0 Quality Results1", 800, 800);
  canvas1->Divide(2, 2);

  canvas1->cd(1);
  fhDcaPosToVertex[0]->SetStats(kFALSE);
  fhDcaPosToVertex[0]->Draw();

  canvas1->cd(2);
  fhChi2[0]->SetStats(kFALSE);
  fhChi2[0]->Draw();

  canvas1->cd(3);
  fhDcaV0ToVertex[0]->SetStats(kFALSE);
  fhDcaV0ToVertex[0]->Draw();

  canvas1->SaveAs(Form("%s_%s.gif", GetName(), canvas1->GetName()));

  TCanvas* canvas2 = new TCanvas(Form("%s_2", GetName()), "V0 Quality Results2", 1200, 800);
  canvas2->Divide(3, 2);

  canvas2->cd(1);
  fhDcaV0Daughters[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhDcaV0Daughters[0]->Draw();

  canvas2->cd(2);
  fhRadius[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhRadius[0]->Draw();


  canvas2->cd(4);
  fhCosinePointingAngle[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhCosinePointingAngle[0]->Draw();

  canvas2->cd(5);
  fhOnFlyStatus[0]->SetStats(kFALSE);
  gPad->SetLogy();
  fhOnFlyStatus[0]->Draw();

  canvas2->SaveAs(Form("%s_%s.gif", GetName(), canvas2->GetName()));

  TCanvas* canvas3 = new TCanvas(Form("%s_4", GetName()), "V0 Quality Results3", 800, 500);
  canvas3->Divide(2, 1);

  canvas3->cd(1);
  fhCutStatistics->SetStats(kFALSE);
  fhCutStatistics->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  fhCutStatistics->Draw();

  canvas3->cd(2);
  fhCutCorrelation->SetStats(kFALSE);
  fhCutCorrelation->LabelsOption("v");
  gPad->SetBottomMargin(0.3);
  gPad->SetLeftMargin(0.3);
  fhCutCorrelation->Draw("COLZ");

  canvas3->SaveAs(Form("%s_%s.gif", GetName(), canvas3->GetName()));
}

