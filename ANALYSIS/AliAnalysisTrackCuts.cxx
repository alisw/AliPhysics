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

/* $Id$ */

//-----------------------------------------------------------------
//           AliAnalysisTrackCuts class
//   This is the class to deal with the event and track level cuts
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------



//ROOT
#include <TROOT.h>
#include <TObject.h>
#include <TSystem.h>
#include <TObject.h>
#include <TPaveText.h>
#include <TText.h>
#include <TLine.h>
#include <TCanvas.h>

#include "AliLog.h"
#include "AliESDtrack.h"
#include "AliESD.h"

#include "AliAnalysisTrackCuts.h"

ClassImp(AliAnalysisTrackCuts)

//----------------------------------------//
AliAnalysisTrackCuts::AliAnalysisTrackCuts()
{
  Reset();
}

//----------------------------------------//
AliAnalysisTrackCuts::~AliAnalysisTrackCuts()
{
}

//----------------------------------------//
void AliAnalysisTrackCuts::Reset()
{
  fPxMin = -1000.0;
  fPxMax = 1000.0; 
  fPyMin = -1000.0;
  fPyMax = 1000.0;  
  fPzMin = -1000.0;
  fPzMax = 1000.0;
  fPtMin = 0.0;
  fPtMax = 1000.0;
  fPMin = 0.0;
  fPMax = 1000.0;
  fBrMin = 0.0;
  fBrMax = 1000.0;
  fBzMin = 0.0;
  fBzMax = 1000.0;
  fEtaMin = -100.0;
  fEtaMax = 100.0;
  fRapMin = -100.0;
  fRapMax = 100.0;
  
  fP = 0;  
  fPt = 0; 
  fPx = 0;  
  fPy = 0; 
  fPz = 0; 
  fbr = 0; 
  fbz = 0; 
  fEta = 0;  
  fRap = 0; 
  fTotalTracks = 0; 
  fAcceptedTracks = 0; 

  fFlagP = 0;  
  fFlagPt = 0;  
  fFlagPx = 0;  
  fFlagPy = 0;  
  fFlagPz = 0;  
  fFlagEta = 0;  
  fFlagRap = 0;  
  fFlagbr = 0;  
  fFlagbz = 0;  
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetPxRange(Float_t r1, Float_t r2)
{
  fPxMin = r1;
  fPxMax = r2;
  fFlagPx = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetPyRange(Float_t r1, Float_t r2)
{
  fPyMin = r1;
  fPyMax = r2; 
  fFlagPy = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetPzRange(Float_t r1, Float_t r2)
{
  fPzMin = r1;
  fPzMax = r2; 
  fFlagPy = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetPRange(Float_t r1, Float_t r2)
{
  fPMin = r1;
  fPMax = r2; 
  fFlagPz = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetPtRange(Float_t r1, Float_t r2)
{
  fPtMin = r1;
  fPtMax = r2; 
  fFlagPt = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetBrRange(Float_t r1, Float_t r2)
{
  fBrMin = r1;
  fBrMax = r2; 
  fFlagbr = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetBzRange(Float_t r1, Float_t r2)
{
  fBzMin = r1;
  fBzMax = r2; 
  fFlagbz = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetEtaRange(Float_t r1, Float_t r2)
{
  fEtaMin = r1;
  fEtaMax = r2; 
  fFlagEta = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::SetRapRange(Float_t r1, Float_t r2)
{
  fRapMin = r1;
  fRapMax = r2; 
  fFlagRap = 1;
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetTrackStats()
{
  AliInfo(Form("Total number of tracks: %d",fTotalTracks));
  AliInfo(Form("Total number of accepted tracks: %d",fAcceptedTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetPStats()
{
  AliInfo(Form("P range: [%f,%f]",fPMin,fPMax));
  if(fTotalTracks != 0)
    AliInfo(Form("Tracks rejected: %f",100.0*fP/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetPtStats()
{
  AliInfo(Form("Pt range: [%f,%f]",fPtMin,fPtMax));
  if(fTotalTracks != 0) 
    AliInfo(Form("Tracks rejected: %f",100.0*fPt/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetPxStats()
{
  AliInfo(Form("Px range: [%f,%f]",fPxMin,fPxMax));
  if(fTotalTracks != 0) 
    AliInfo(Form("Tracks rejected: %f",100.0*fPx/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetPyStats()
{
  AliInfo(Form("Py range: [%f,%f]",fPyMin,fPyMax));
  if(fTotalTracks != 0) 
    AliInfo(Form("Tracks rejected: %f",100.0*fPy/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetPzStats()
{
  AliInfo(Form("Pz range: [%f,%f]",fPzMin,fPzMax));
  if(fTotalTracks != 0) 
    AliInfo(Form("Tracks rejected: %f",100.0*fPz/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetEtaStats()
{
  AliInfo(Form("eta range: [%f,%f]",fEtaMin,fEtaMax));
  if(fTotalTracks != 0)
    AliInfo(Form("Tracks rejected: %f",100.0*fEta/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetRapStats()
{
  AliInfo(Form("y range: [%f,%f]",fRapMin,fRapMax));
  if(fTotalTracks != 0)
    AliInfo(Form("Tracks rejected: %f",100.0*fRap/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetBrStats()
{
  AliInfo(Form("br range: [%f,%f]",fBrMin,fBrMax));
  if(fTotalTracks != 0) 
    AliInfo(Form("Tracks rejected: %f",100.0*fbr/fTotalTracks)); 
}

//----------------------------------------//
void AliAnalysisTrackCuts::GetBzStats()
{
  AliInfo(Form("bz range: [%f,%f]",fBzMin,fBzMax));
  if(fTotalTracks != 0)
    AliInfo(Form("Tracks rejected: %f",100.0*fbz/fTotalTracks)); 
}


//----------------------------------------//
Bool_t AliAnalysisTrackCuts::IsAccepted(AliESD *esd ,AliESDtrack *esdtrack)
{
  fTotalTracks++;
  
  //momentum related calculations
  Double_t p[3];
  esdtrack->GetPxPyPz(p);
  Float_t P = sqrt(pow(p[0],2) + pow(p[1],2) + pow(p[2],2));
  Float_t Pt = sqrt(pow(p[0],2) + pow(p[1],2));
  Float_t E = sqrt(pow(esdtrack->GetMass(),2) + pow(P,2));

  //y-eta related calculations
  Float_t eta = 0.5*log((P - p[2])/(P + p[2]));
  Float_t y = 0.5*log((E - p[2])/(E + p[2]));
 
  //impact parameter related calculations
  Double_t TrackPosition[3];
  esdtrack->GetXYZ(TrackPosition);
  const AliESDVertex * VertexIn = esd->GetVertex();
  Double_t VertexPosition[3];
  VertexIn->GetXYZ(VertexPosition);		    
  for (Int_t ii=0; ii<3; ii++) TrackPosition[ii] -= VertexPosition[ii];
		    
  Float_t br = Float_t(TMath::Sqrt(pow(TrackPosition[0],2) + pow(TrackPosition[1],2)));
  Float_t bz = Float_t(TMath::Abs(TrackPosition[2]));
 
  if((P < fPMin) || (P > fPMax))
    {
      fP++;
      return kFALSE;
    }
  if((Pt < fPtMin) || (Pt > fPtMax))
    {
      fPt++;
      return kFALSE;
    }
  if((p[0] < fPxMin) || (p[0] > fPxMax))
    {
      fPx++;
      return kFALSE;
    }
  if((p[1] < fPyMin) || (p[1] > fPyMax))
    {
      fPy++;
      return kFALSE;
    }
  if((p[2] < fPzMin) || (p[2] > fPzMax))
    {
      fPz++;
      return kFALSE;
    } 
  if((br < fBrMin) || (br > fBrMax))
    {
      fbr++;
      return kFALSE;
    }
  if((bz < fBzMin) || (bz > fBzMax))
    {
      fbz++;
      return kFALSE;
    }
  if((eta < fEtaMin) || (eta > fEtaMax))
    {
      fEta++;
      return kFALSE;
    }
  if((y < fRapMin) || (y > fRapMax))
    {
      fRap++;
      return kFALSE;
    }
  
  fAcceptedTracks++;
  
  return kTRUE;
}


//----------------------------------------//
TPaveText *AliAnalysisTrackCuts::GetTrackCuts()
{
  TCanvas *ccuts2 = new TCanvas("ccuts2","Track cuts",410,10,400,400);
  ccuts2->SetFillColor(10);
  ccuts2->SetHighLightColor(10);

  TPaveText *pave = new TPaveText(0.01,0.01,0.98,0.98);
  pave->SetFillColor(3);
  Char_t CutName[256];
 
  TLine *l1 = pave->AddLine(0,0.89,1,0.89);
  l1->SetLineWidth(2);
  TLine *l2 = pave->AddLine(0,0.79,1,0.79);
  l2->SetLineWidth(2);
  TLine *l3 = pave->AddLine(0,0.69,1,0.69);
  l3->SetLineWidth(2);
  TLine *l4 = pave->AddLine(0,0.59,1,0.59);
  l4->SetLineWidth(2);
  TLine *l5 = pave->AddLine(0,0.49,1,0.49);
  l5->SetLineWidth(2);
  TLine *l6 = pave->AddLine(0,0.39,1,0.39);
  l6->SetLineWidth(2);
  TLine *l7 = pave->AddLine(0,0.29,1,0.29);
  l7->SetLineWidth(2);
  TLine *l8 = pave->AddLine(0,0.19,1,0.19);
  l8->SetLineWidth(2);
  TLine *l9 = pave->AddLine(0,0.09,1,0.09);
  l9->SetLineWidth(2);

  sprintf(CutName,"Total number of tracks: %d",fTotalTracks);
  TText *t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"Total number of accepted tracks: %d",fAcceptedTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"P range: [%f,%f]",fPMin,fPMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fP/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"Pt range: [%f,%f]",fPtMin,fPtMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fPt/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);

  sprintf(CutName,"Px range: [%f,%f]",fPxMin,fPxMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fPx/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"Py range: [%f,%f]",fPyMin,fPyMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fPy/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"Pz range: [%f,%f]",fPzMin,fPzMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fPz/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  
  sprintf(CutName,"br range: [%f,%f]",fBrMin,fBrMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fbr/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
 
  sprintf(CutName,"bz range: [%f,%f]",fBzMin,fBzMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fbz/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);

  sprintf(CutName,"eta range: [%f,%f]",fEtaMin,fEtaMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fEta/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);

  sprintf(CutName,"y range: [%f,%f]",fRapMin,fRapMax);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);
  sprintf(CutName,"Tracks rejected: %f",100.0*fRap/fTotalTracks);
  t1 = pave->AddText(CutName);
  t1->SetTextColor(1);
  t1->SetTextSize(0.04);
  t1->SetTextAlign(11);

  return pave;
}

//----------------------------------------//
void AliAnalysisTrackCuts::PrintTrackCuts()
{
  //GetTrackCuts()->Draw();

  AliInfo(Form("**************TRACK CUTS**************"));
  GetTrackStats();
  if(fFlagP)
    GetPStats();
  if(fFlagPt)
    GetPtStats();
  if(fFlagPx)
    GetPxStats();
  if(fFlagPy)
    GetPyStats();
  if(fFlagPz)
    GetPzStats();
  if(fFlagEta)
    GetEtaStats();
  if(fFlagRap)
    GetRapStats();
  if(fFlagbr)
    GetBrStats();
  if(fFlagbz)
    GetBzStats(); 
  AliInfo(Form("**************************************"));
}
