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

#include "Riostream.h"              //needed as include
#include "AliFlowCommonConstants.h" //needed as include
#include "AliFlowCommonHist.h"
#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"

#include "TString.h" 
#include "TProfile.h"
#include "TMath.h"   //needed as include
#include "AliFlowVector.h"

class TH1F;
class TH1D;

// AliFlowCommonHist:
//
// Description: Class to organise common histograms for Flow Analysis

// authors: N. van der Kolk (kolk@nikhef.nl) and A. Bilandzic (anteb@nikhef.nl)


ClassImp(AliFlowCommonHist)

//-----------------------------------------------------------------------

  AliFlowCommonHist::AliFlowCommonHist(TString input):
   fEvent(0),
   fTrack(0), 
   fHistMultOrig(0),
   fHistMultInt(0),
   fHistMultDiff(0),
   fHistPtInt(0),
   fHistPtDiff(0),
   fHistPhiInt(0),
   fHistPhiDiff(0),
   fHistEtaInt(0),
   fHistEtaDiff(0),
   fHistProMeanPtperBin(0),
   fHistQ(0)
 {
  //constructor creating histograms 
  Int_t fNbinsMult = AliFlowCommonConstants::GetNbinsMult();
  Int_t fNbinsPt = AliFlowCommonConstants::GetNbinsPt();
  Int_t fNbinsPhi = AliFlowCommonConstants::GetNbinsPhi();
  Int_t fNbinsEta = AliFlowCommonConstants::GetNbinsEta();
  Int_t fNbinsQ = AliFlowCommonConstants::GetNbinsQ();
  TString name;

  Double_t  fMultMin = AliFlowCommonConstants::GetMultMin();            
  Double_t  fMultMax = AliFlowCommonConstants::GetMultMax();
  Double_t  fPtMin = AliFlowCommonConstants::GetPtMin();	     
  Double_t  fPtMax = AliFlowCommonConstants::GetPtMax();
  Double_t  fPhiMin = AliFlowCommonConstants::GetPhiMin();	     
  Double_t  fPhiMax = AliFlowCommonConstants::GetPhiMax();
  Double_t  fEtaMin = AliFlowCommonConstants::GetEtaMin();	     
  Double_t  fEtaMax = AliFlowCommonConstants::GetEtaMax();	     
  Double_t  fQMin = AliFlowCommonConstants::GetQMin();	     
  Double_t  fQMax = AliFlowCommonConstants::GetQMax();	
  
  cout<<"The settings for the common histograms are as follows:"<<endl;
  cout<<"Multiplicity: "<<fNbinsMult<<" bins between "<<fMultMin<<" and "<<fMultMax<<endl;
  cout<<"Pt: "<<fNbinsPt<<" bins between "<<fPtMin<<" and "<<fPtMax<<endl;
  cout<<"Phi: "<<fNbinsPhi<<" bins between "<<fPhiMin<<" and "<<fPhiMax<<endl;
  cout<<"Eta: "<<fNbinsEta<<" bins between "<<fEtaMin<<" and "<<fEtaMax<<endl;
  cout<<"Q: "<<fNbinsQ<<" bins between "<<fQMin<<" and "<<fQMax<<endl;

  //Multiplicity
  name = "Control_Flow_OrigMult_";
  name +=input;
  fHistMultOrig = new TH1F(name.Data(), name.Data(),fNbinsMult, fMultMin, fMultMax);
  fHistMultOrig ->SetXTitle("Original Multiplicity");
  fHistMultOrig ->SetYTitle("Counts");

  name = "Control_Flow_MultInt_";
  name +=input;
  fHistMultInt = new TH1F(name.Data(), name.Data(),fNbinsMult, fMultMin, fMultMax);
  fHistMultInt ->SetXTitle("Multiplicity for integrated flow");
  fHistMultInt ->SetYTitle("Counts");

  name = "Control_Flow_MultDiff_";
  name +=input;
  fHistMultDiff = new TH1F(name.Data(), name.Data(),fNbinsMult, fMultMin, fMultMax);
  fHistMultDiff ->SetXTitle("Multiplicity for differential flow");
  fHistMultDiff ->SetYTitle("Counts");

  //Pt
  name = "Control_Flow_PtInt_";
  name +=input;
  fHistPtInt = new TH1F(name.Data(), name.Data(),fNbinsPt, fPtMin, fPtMax); 
  fHistPtInt ->SetXTitle("Pt (GeV/c) for integrated flow");
  fHistPtInt ->SetYTitle("Counts");

  name = "Control_Flow_PtDiff_";
  name +=input;
  fHistPtDiff = new TH1F(name.Data(), name.Data(),fNbinsPt, fPtMin, fPtMax); 
  //binning has to be the same as for fHistProVPt! use to get Nprime!
  fHistPtDiff ->SetXTitle("Pt (GeV/c) for differential flow");
  fHistPtDiff ->SetYTitle("Counts");

  //Phi
  name = "Control_Flow_PhiInt_";
  name +=input;
  fHistPhiInt = new TH1F(name.Data(), name.Data(),fNbinsPhi, fPhiMin, fPhiMax);
  fHistPhiInt ->SetXTitle("Phi for integrated flow");
  fHistPhiInt ->SetYTitle("Counts");

  name = "Control_Flow_PhiDiff_";
  name +=input;
  fHistPhiDiff = new TH1F(name.Data(), name.Data(),fNbinsPhi, fPhiMin, fPhiMax);
  fHistPhiDiff ->SetXTitle("Phi for differential flow");
  fHistPhiDiff ->SetYTitle("Counts");

  //Eta
  name = "Control_Flow_EtaInt_";
  name +=input;
  fHistEtaInt = new TH1F(name.Data(), name.Data(),fNbinsEta, fEtaMin, fEtaMax);
  fHistEtaInt ->SetXTitle("Eta for integrated flow");
  fHistEtaInt ->SetYTitle("Counts");

  name = "Control_Flow_EtaDiff_";
  name +=input;
  fHistEtaDiff = new TH1F(name.Data(), name.Data(),fNbinsEta, fEtaMin, fEtaMax);
  fHistEtaDiff ->SetXTitle("Eta for differential flow");
  fHistEtaDiff ->SetYTitle("Counts");

  //Mean Pt per pt bin 
  name = "Control_FlowPro_MeanPtperBin_";
  name +=input;
  fHistProMeanPtperBin = new TProfile(name.Data(), name.Data(),fNbinsPt,fPtMin,fPtMax);
  fHistProMeanPtperBin ->SetXTitle("Pt");
  fHistProMeanPtperBin ->SetYTitle("<Pt>");

  //Q vector
  name = "Control_Flow_Q_";
  name +=input;
  fHistQ = new TH1F(name.Data(), name.Data(),fNbinsQ, fQMin, fQMax);
  fHistQ ->SetXTitle("Qvector/Mult");
  fHistQ ->SetYTitle("Counts");  
}


//----------------------------------------------------------------------- 

AliFlowCommonHist::~AliFlowCommonHist()
{
  //deletes histograms
  delete fHistMultOrig;      
  delete fHistMultInt;       
  delete fHistMultDiff;      
  delete fHistPtInt;         
  delete fHistPtDiff;       
  delete fHistPhiInt;        
  delete fHistPhiDiff;       
  delete fHistEtaInt;        
  delete fHistEtaDiff;
  delete fHistProMeanPtperBin;
  delete fHistQ;
}


//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHist::FillControlHistograms(AliFlowEventSimple* fEvent)
{
  //Fills the control histograms
  if (!fEvent){
    cout<<"##### FillControlHistograms: FlowEvent pointer null"<<endl;
    return kFALSE;
  }

  Double_t fPt, fPhi, fEta;


  //fill the histograms
  Int_t fNumberOfTracks = fEvent->NumberOfTracks();
  fHistMultOrig->Fill(fNumberOfTracks);

  AliFlowVector fQ = fEvent->GetQ(); 
  //weight by the Multiplicity
  Double_t fQX = fQ.X()/fQ.GetMult();
  Double_t fQY = fQ.Y()/fQ.GetMult();
  fQ.Set(fQX,fQY);
  fHistQ->Fill(fQ.Mod());

  Int_t fMultInt = 0;
  Int_t fMultDiff = 0;
  for (Int_t i=0;i<fNumberOfTracks;i++) {
    fTrack = fEvent->GetTrack(i);
    if (fTrack ) {
      if (fTrack->UseForIntegratedFlow()){
	fPt = fTrack->Pt();
	fHistPtInt->Fill(fPt);
	fPhi = fTrack->Phi();
	if (fPhi<0.) fPhi+=2*TMath::Pi();
	fHistPhiInt->Fill(fPhi);
	fEta = fTrack->Eta();
	fHistEtaInt->Fill(fEta);
	fMultInt++;
      }
      if (fTrack->UseForDifferentialFlow()){
	fPt = fTrack->Pt();
	fHistPtDiff->Fill(fPt);
	fPhi = fTrack->Phi();
	if (fPhi<0.) fPhi+=2*TMath::Pi();
	fHistPhiDiff->Fill(fPhi);
	fEta = fTrack->Eta();
	fHistEtaDiff->Fill(fEta);
	fHistProMeanPtperBin->Fill(fPt,fPt);
	fMultDiff++;
      }
    } //track
  } //loop over tracks
  
  fHistMultInt->Fill(fMultInt);
  fHistMultDiff->Fill(fMultDiff);

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBin(Int_t fBin)
{
  //get entries in bin fBin from fHistPtDiff
  Double_t fEntries = fHistPtDiff->GetBinContent(fBin);

  return fEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetMeanPt(Int_t fBin)
{  
  //Get entry from bin fBin from fHistProMeanPtperBin
  Double_t fMeanPt = fHistProMeanPtperBin->GetBinContent(fBin);

  return fMeanPt;
  
}


