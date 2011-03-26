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
//           Balance Function class
//   This is the class to deal with the Balance Function analysis
//   Origin: Panos Christakoglou, UOA-CERN, Panos.Christakoglou@cern.ch
//-----------------------------------------------------------------


//ROOT
#include <Riostream.h>
#include <TMath.h>
#include <TAxis.h>
#include <TLorentzVector.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TH1F.h>

#include "AliVParticle.h"
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"

#include "AliBalance.h"

ClassImp(AliBalance)

//____________________________________________________________________//
AliBalance::AliBalance() :
  TObject(), 
  fAnalysisLevel("ESD"), fNumberOfBins(0),
  fAnalysisType(0), fAnalyzedEvents(0), fP2Start(0),
  fP2Stop(0), fP2Step(0), fNn(0), fNp(0),
  fHistfNnn(new TH1F("fHistfNnn","(--) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)),
  fHistfNpp(new TH1F("fHistfNpp","(++) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)),
  fHistfNpn(new TH1F("fHistfNpn","(+-) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)) {
  // Default constructor
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 

  switch(fAnalysisType) {
  case 0:
    fHistfNnn->GetXaxis()->SetTitle("#Delta y");
    fHistfNpp->GetXaxis()->SetTitle("#Delta y");
    fHistfNpn->GetXaxis()->SetTitle("#Delta y");
    break;
  case 1:
    fHistfNnn->GetXaxis()->SetTitle("#Delta #eta");
    fHistfNpp->GetXaxis()->SetTitle("#Delta #eta");
    fHistfNpn->GetXaxis()->SetTitle("#Delta #eta");
    break;
  case 2:
    fHistfNnn->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    break;
  case 3:
    fHistfNnn->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    break;
  case 4:
    fHistfNnn->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    break;
  case 5:
    fHistfNnn->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    break;
  case 6:
    fHistfNnn->GetXaxis()->SetTitle("#Delta #phi");
    fHistfNpp->GetXaxis()->SetTitle("#Delta #phi");
    fHistfNpn->GetXaxis()->SetTitle("#Delta #phi");
    break;
  default:
    break;
  }
}

//____________________________________________________________________//
AliBalance::AliBalance(Double_t p2Start, Double_t p2Stop, Int_t p2Bins) :
  TObject(), fAnalysisLevel("ESD"),
  fNumberOfBins(p2Bins), fAnalysisType(0), 
  fAnalyzedEvents(0), fP2Start(p2Start), fP2Stop(p2Stop), 
  fP2Step(TMath::Abs(fP2Start - fP2Stop) / (Double_t)fNumberOfBins), 
  fNn(0), fNp(0),
  fHistfNnn(new TH1F("fHistfNnn","(--) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)),
  fHistfNpp(new TH1F("fHistfNpp","(++) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)),
  fHistfNpn(new TH1F("fHistfNpn","(+-) component;;Entries",
		     fNumberOfBins,fP2Start,fP2Stop)) {
  // Constructor
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 

  switch(fAnalysisType) {
  case 0:
    fHistfNnn->GetXaxis()->SetTitle("#Delta y");
    fHistfNpp->GetXaxis()->SetTitle("#Delta y");
    fHistfNpn->GetXaxis()->SetTitle("#Delta y");
    break;
  case 1:
    fHistfNnn->GetXaxis()->SetTitle("#Delta #eta");
    fHistfNpp->GetXaxis()->SetTitle("#Delta #eta");
    fHistfNpn->GetXaxis()->SetTitle("#Delta #eta");
    break;
  case 2:
    fHistfNnn->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    break;
  case 3:
    fHistfNnn->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    break;
  case 4:
    fHistfNnn->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    break;
  case 5:
    fHistfNnn->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    fHistfNpp->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    fHistfNpn->GetXaxis()->SetTitle("q_{inv.} (GeV/c)");
    break;
  case 6:
    fHistfNnn->GetXaxis()->SetTitle("#Delta #phi");
    fHistfNpp->GetXaxis()->SetTitle("#Delta #phi");
    fHistfNpn->GetXaxis()->SetTitle("#Delta #phi");
    break;
  default:
    break;
  }
}

//____________________________________________________________________//
AliBalance::AliBalance(const AliBalance& balance):
  TObject(balance), fAnalysisLevel(balance.fAnalysisLevel),
  fNumberOfBins(balance.fNumberOfBins),
  fAnalysisType(balance.fAnalysisType),
  fAnalyzedEvents(balance.fAnalyzedEvents),
  fP2Start(balance.fP2Start),
  fP2Stop(balance.fP2Stop),
  fP2Step(balance.fP2Step),
  fNn(balance.fNn),
  fNp(balance.fNp),
  fHistfNnn(balance.fHistfNnn), 
  fHistfNpp(balance.fHistfNpp), 
  fHistfNpn(balance.fHistfNpn) {
  //copy constructor
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 
}

//____________________________________________________________________//
AliBalance::~AliBalance() {
  // Destructor
  if(fHistfNnn) delete fHistfNnn;
  if(fHistfNpp) delete fHistfNpp;
  if(fHistfNpn) delete fHistfNpn;
}

//____________________________________________________________________//
void AliBalance::SetNnn(Double_t *nn) {
  // Setter of the Nnn term
  for(Int_t i = 0; i < fNumberOfBins; i++) fNnn[i] = nn[i];
}

//____________________________________________________________________//
void AliBalance::SetNpp(Double_t *pp) {
  // Setter of the Npp term
  for(Int_t i = 0; i < fNumberOfBins; i++) fNpp[i] = pp[i];
}

//____________________________________________________________________//
void AliBalance::SetNpn(Double_t *pn) {
  // Setter of the Npn term
  for(Int_t i = 0; i < fNumberOfBins; i++) fNpn[i] = pn[i];
}

//____________________________________________________________________//
void AliBalance::SetNumberOfBins(Int_t ibins) {
  // Sets the number of bins for the analyzed interval
  fNumberOfBins = ibins;
}

//____________________________________________________________________//
void AliBalance::SetInterval(Double_t p2Start, Double_t p2Stop) {
  // Sets the analyzed interval. 
  // The analysis variable is set by SetAnalysisType
  fP2Start = p2Start;
  fP2Stop = p2Stop;
  fP2Step = TMath::Abs(p2Start - p2Stop) / (Double_t)fNumberOfBins;
}

//____________________________________________________________________//
void AliBalance::SetAnalysisType(Int_t iType) {
  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  this->fAnalysisType = iType; 
  if(fAnalysisType==0) {
    cout<<" ====================== "<<endl;
    cout<<"||Analysis selected: y||"<<endl;
    cout<<" ====================== "<<endl;
  } 
  else if(fAnalysisType==1) {
    cout<<" ======================== "<<endl;
    cout<<"||Analysis selected: eta||"<<endl;
    cout<<" ======================== "<<endl;
  }
  else if(fAnalysisType==2) {
    cout<<" ========================== "<<endl;
    cout<<"||Analysis selected: Qlong||"<<endl;
    cout<<" ========================== "<<endl;
  }
  else if(fAnalysisType==3) {
    cout<<" ========================= "<<endl;
    cout<<"||Analysis selected: Qout||"<<endl;
    cout<<" ========================= "<<endl;
  }
  else if(fAnalysisType==4) {
    cout<<" ========================== "<<endl;
    cout<<"||Analysis selected: Qside||"<<endl;
    cout<<" ========================== "<<endl;
  }
  else if(fAnalysisType==5) {
    cout<<" ========================= "<<endl;
    cout<<"||Analysis selected: Qinv||"<<endl;
    cout<<" ========================= "<<endl;
  }
  else if(fAnalysisType==6) {
    cout<<" ======================== "<<endl;
    cout<<"||Analysis selected: phi||"<<endl;
    cout<<" ======================== "<<endl;
  }
  else {
    cout<<"Selection of analysis mode failed!!!"<<endl;
    cout<<"Choices are: 0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi"<<endl;
    abort();
  }
}

//____________________________________________________________________//
void AliBalance::PrintAnalysisSettings() {
  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  TString analysisType;
  switch(fAnalysisType) {
  case 0:
    analysisType = "Rapidity"; 
    break;
  case 1:
    analysisType = "Pseudo-rapidity"; 
    break;
  case 2:
    analysisType = "Qlong"; 
    break;
  case 3:
    analysisType = "Qout"; 
    break;
  case 4:
    analysisType = "Qside"; 
    break;
  case 5:
    analysisType = "Qinv"; 
    break;
  case 6:
    analysisType = "Phi"; 
    break;
  default:
    break;
  }
  
  Printf("======================================");
  Printf("Analysis level: %s",fAnalysisLevel.Data());
  Printf("Analysis type: %s",analysisType.Data());
  Printf("Analyzed interval (min.): %lf",fP2Start);
  Printf("Analyzed interval (max.): %lf",fP2Stop);
  Printf("Number of bins: %d",fNumberOfBins);
  Printf("Step: %lf",fP2Step);
  Printf("======================================");
}

//____________________________________________________________________//
void AliBalance::CalculateBalance(TObjArray *gTrackArray) {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t ibin = 0;
  
  AliVParticle* track = 0;
  AliVParticle* track1 = 0;
  AliVParticle* track2 = 0;
    
  //Printf("(AliBalance) Number of tracks: %d",gTrackArray->GetEntries());
  Int_t gNtrack = gTrackArray->GetEntries();
  for(i = 0; i < gNtrack; i++) {
    if(fAnalysisLevel == "ESD")
      track = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "AOD")
      track = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
    else if(fAnalysisLevel == "MC")
      track = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

    if(track) {
      Short_t charge = track->Charge();
      if(charge > 0) fNp += 1.;
      if(charge < 0) fNn += 1.;
    }
    else continue;
  }
  //Printf("Np: %lf - Nn: %lf",fNp,fNn);

  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  if(fAnalysisType==0) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pZ1 = 0., energy1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pZ1 = track1->Pz();
	energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
			      TMath::Power(track1->M(),2));
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pZ2 = track2->Pz();
	  Double_t energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
					 TMath::Power(track2->M(),2));
	  
	  Double_t rap1 = 0.5*log((energy1 + pZ1)/(energy1 - pZ1)); 
	  Double_t rap2 = 0.5*log((energy2 + pZ2)/(energy2 - pZ2)); 
	  Double_t dy = TMath::Abs(rap1 - rap2);
	  ibin = Int_t(dy/fP2Step);
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 0
  if(fAnalysisType==1) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pZ1 = 0., p1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pZ1 = track1->Pz();
	p1 = track1->P();
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pZ2 = track2->Pz();
	  Double_t p2 = track2->P();
	  Double_t eta1 = 0.5*log((p1 + pZ1)/(p1 - pZ1)); 
	  Double_t eta2 = 0.5*log((p2 + pZ2)/(p2 - pZ2)); 
	  Double_t deta = TMath::Abs(eta1 - eta2);
	  ibin = Int_t(deta/fP2Step);
	  
	  if((charge1 > 0.)&&(charge2 > 0.)) fNpp[ibin] += 1.;
	  if((charge1 < 0.)&&(charge2 < 0.)) fNnn[ibin] += 1.;
	  if((charge1 > 0.)&&(charge2 < 0.)) fNpn[ibin] += 1.;
	  if((charge1 < 0.)&&(charge2 > 0.)) fNpn[ibin] += 1.;
	  //Printf("charge1: %d - eta1: %lf - charge2: %d - eta2: %lf - deta: %lf - ibin: %d - fNpp: %lf - fNnn: %lf - fNpn: %lf",charge1,eta1,charge2,eta2,deta,ibin,fNpp[ibin],fNnn[ibin],fNpn[ibin]);      
	}
	else continue;
      }
    }
  }//case 1
  if(fAnalysisType==2) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pX1 = 0., pY1 = 0., pZ1 = 0., energy1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pX1 = track1->Px();
	pY1 = track1->Py();
	pZ1 = track1->Pz();
	energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
				       TMath::Power(track1->M(),2));
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pX2 = track2->Px();
	  Double_t pY2 = track2->Py();
	  Double_t pZ2 = track2->Pz();
	  Double_t energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
					 TMath::Power(track2->M(),2));
	  Double_t eTot = energy1 + energy2;
	  Double_t pxTot = pX1 + pX2;
	  Double_t pyTot = pY1 + pY2;
	  Double_t pzTot = pZ1 + pZ2;
	  Double_t q0Tot = energy1 - energy2;
	  Double_t qzTot = pZ1 - pZ2;
	  Double_t snn = TMath::Power(eTot,2) - TMath::Power(pxTot,2) - TMath::Power(pyTot,2) - TMath::Power(pzTot,2);
	  Double_t ptTot = TMath::Sqrt( TMath::Power(pxTot,2) + TMath::Power(pyTot,2));
	  Double_t qLong = TMath::Abs(eTot*qzTot - pzTot*q0Tot)/TMath::Sqrt(snn + TMath::Power(ptTot,2));
	  ibin = Int_t(qLong/fP2Step);
	  //cout<<ibin<<endl;
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 2
  if(fAnalysisType==3) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pX1 = 0., pY1 = 0., pZ1 = 0., energy1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pX1 = track1->Px();
	pY1 = track1->Py();
	pZ1 = track1->Pz();
	energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
			      TMath::Power(track1->M(),2));
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pX2 = track2->Px();
	  Double_t pY2 = track2->Py();
	  Double_t pZ2 = track2->Pz();
	  Double_t energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
					 TMath::Power(track2->M(),2));
	  Double_t eTot = energy1 + energy2;
	  Double_t pxTot = pX1 + pX2;
	  Double_t pyTot = pY1 + pY2;
	  Double_t pzTot = pZ1 + pZ2;
	  Double_t qxTot = pX1 - pX2;
	  Double_t qyTot = pY1 - pY2;
	  Double_t snn = TMath::Power(eTot,2) - TMath::Power(pxTot,2) - TMath::Power(pyTot,2) - TMath::Power(pzTot,2);
	  Double_t ptTot = TMath::Sqrt( TMath::Power(pxTot,2) + TMath::Power(pyTot,2));
	  Double_t qOut = TMath::Sqrt(snn/(snn + TMath::Power(ptTot,2))) * TMath::Abs(pxTot*qxTot + pyTot*qyTot)/ptTot;
	  ibin = Int_t(qOut/fP2Step);
	  //cout<<ibin<<endl;
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 3
  if(fAnalysisType==4) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pX1 = 0., pY1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pX1 = track1->Px();
	pY1 = track1->Py();
	//Double_t energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
	//TMath::Power(track1->M(),2));
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pX2 = track2->Px();
	  Double_t pY2 = track2->Py();
	  //Double_t energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
	  //TMath::Power(track2->M(),2));
	  //Double_t eTot = energy1 + energy2;
	  Double_t pxTot = pX1 + pX2;
	  Double_t pyTot = pY1 + pY2;
	  Double_t qxTot = pX1 - pX2;
	  Double_t qyTot = pY1 - pY2;
	  Double_t ptTot = TMath::Sqrt( TMath::Power(pxTot,2) + TMath::Power(pyTot,2));
	  Double_t qSide = TMath::Abs(pxTot*qyTot - pyTot*qxTot)/ptTot;
	  ibin = Int_t(qSide/fP2Step);
	  //cout<<ibin<<endl;
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 4
  if(fAnalysisType==5) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t pX1 = 0., pY1 = 0., pZ1 = 0., energy1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	pX1 = track1->Px();
	pY1 = track1->Py();
	pZ1 = track1->Pz();
	energy1 = TMath::Sqrt(TMath::Power(track1->P(),2) +
			      TMath::Power(track1->M(),2));
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t pX2 = track2->Px();
	  Double_t pY2 = track2->Py();
	  Double_t pZ2 = track2->Pz();
	  Double_t energy2 = TMath::Sqrt(TMath::Power(track2->P(),2) +
					 TMath::Power(track2->M(),2));
	  Double_t q0Tot = energy1 - energy2;
	  Double_t qxTot = pX1 - pX2;
	  Double_t qyTot = pY1 - pY2;
	  Double_t qzTot = pZ1 - pZ2;
	  Double_t qInv = TMath::Sqrt(TMath::Abs(-TMath::Power(q0Tot,2) +TMath::Power(qxTot,2) +TMath::Power(qyTot,2) +TMath::Power(qzTot,2)));
	  ibin = Int_t(qInv/fP2Step);
	  //cout<<ibin<<endl;
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 5	
  if(fAnalysisType==6) {
    for(i = 1; i < gNtrack; i++) {
      if(fAnalysisLevel == "ESD")
	track1 = dynamic_cast<AliESDtrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "AOD")
	track1 = dynamic_cast<AliAODTrack *>(gTrackArray->At(i));
      else if(fAnalysisLevel == "MC")
	track1 = dynamic_cast<AliMCParticle *>(gTrackArray->At(i));

      Short_t charge1 = 0;
      Double_t phi1 = 0.;
      if(track1) {
	charge1 = track1->Charge();
	phi1 = track1->Phi();
      }
      else continue;

      for(j = 0; j < i; j++) {
	if(fAnalysisLevel == "ESD")
	  track2 = dynamic_cast<AliESDtrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "AOD")
	  track2 = dynamic_cast<AliAODTrack *>(gTrackArray->At(j));
	else if(fAnalysisLevel == "MC")
	  track2 = dynamic_cast<AliMCParticle *>(gTrackArray->At(j));

	if(track2) {
	  Short_t charge2 = track2->Charge();
	  Double_t phi2 = track2->Phi();
	  Double_t dphi = TMath::Abs(phi1 - phi2);
	  ibin = Int_t(dphi/fP2Step);
	  if((charge1 > 0)&&(charge2 > 0)) fNpp[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 < 0)) fNnn[ibin] += 1.;
	  if((charge1 > 0)&&(charge2 < 0)) fNpn[ibin] += 1.;
	  if((charge1 < 0)&&(charge2 > 0)) fNpn[ibin] += 1.;
	}
	else continue;
      }
    }
  }//case 6

  /*for(Int_t i = 0; i < fNumberOfBins; i++) 
    Printf("bin: %d - Npp: %lf - Nnn: %lf - Nnp: %lf - Npn: %lf",i,fNpp[i],fNnn[i],fNpn[i],fNpn[i]);*/
}

//____________________________________________________________________//
TH1F *AliBalance::GetHistNnn() {
  //Return the control histogram of the -- component
  for(Int_t iBin = 1; iBin <= fNumberOfBins; iBin++)
    fHistfNnn->SetBinContent(iBin,GetNnn(iBin-1));

  return fHistfNnn;
}

//____________________________________________________________________//
TH1F *AliBalance::GetHistNpp() {
  //Return the control histogram of the ++ component
  for(Int_t iBin = 1; iBin <= fNumberOfBins; iBin++)
    fHistfNpp->SetBinContent(iBin,GetNpp(iBin-1));

  return fHistfNpp;
}

//____________________________________________________________________//
TH1F *AliBalance::GetHistNpn() {
  //Return the control histogram of the +- component
  for(Int_t iBin = 1; iBin <= fNumberOfBins; iBin++)
    fHistfNpn->SetBinContent(iBin,GetNpn(iBin-1));

  return fHistfNpn;
}

//____________________________________________________________________//
Double_t AliBalance::GetBalance(Int_t p2) {
  // Returns the value of the balance function in bin p2
  fB[p2] = 0.5*(((fNpn[p2] - 2.0*fNnn[p2])/fNn) + ((fNpn[p2] - 2.0*fNpp[p2])/fNp))/fP2Step;
  
  return fB[p2];
}
    
//____________________________________________________________________//
Double_t AliBalance::GetError(Int_t p2) {
  // Returns the error on the BF value for bin p2
  ferror[p2] = TMath::Sqrt( Double_t(fNpp[p2])/(Double_t(fNp)*Double_t(fNp)) + Double_t(fNnn[p2])/(Double_t(fNn)*Double_t(fNn)) + Double_t(fNpn[p2])*TMath::Power((0.5/Double_t(fNp) + 0.5/Double_t(fNn)),2))/fP2Step;

  return ferror[p2];
}

//____________________________________________________________________//
void AliBalance::PrintResults() {
  // Prints the results
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t fSumXi = 0.0, fSumBi = 0.0, fSumBiXi = 0.0;
  Double_t fSumBiXi2 = 0.0, fSumBi2Xi2 = 0.0;
  Double_t fSumDeltaBi2 = 0.0, fSumXi2DeltaBi2 = 0.0;
  Double_t deltaBalP2 = 0.0, integral = 0.0;
  Double_t deltaErrorNew = 0.0;

  cout<<"=================================================="<<endl;
  for(Int_t i = 0; i < fNumberOfBins; i++) { 
    //x[i] = fP2Start + fP2Step*i + fP2Step/2;
    x[i] = fP2Step*i + fP2Step/2;
    cout<<"B: "<<fB[i]<<"\t Error: "<<ferror[i]<<"\t bin: "<<x[i]<<endl;
  } 
  cout<<"=================================================="<<endl;
  for(Int_t i = 1; i < fNumberOfBins; i++) {
    fSumXi += x[i];
    fSumBi += fB[i];
    fSumBiXi += fB[i]*x[i];
    fSumBiXi2 += fB[i]*TMath::Power(x[i],2);
    fSumBi2Xi2 += TMath::Power(fB[i],2)*TMath::Power(x[i],2);
    fSumDeltaBi2 +=  TMath::Power(ferror[i],2);
    fSumXi2DeltaBi2 += TMath::Power(x[i],2) * TMath::Power(ferror[i],2);
    
    deltaBalP2 += fP2Step*TMath::Power(ferror[i],2);
    integral += fP2Step*fB[i];
  }
  for(Int_t i = 1; i < fNumberOfBins; i++) deltaErrorNew += ferror[i]*(x[i]*fSumBi - fSumBiXi)/TMath::Power(fSumBi,2);
   
  Double_t integralError = TMath::Sqrt(deltaBalP2);
  
  Double_t delta = fSumBiXi / fSumBi;
  Double_t deltaError = (fSumBiXi / fSumBi) * TMath::Sqrt(TMath::Power((TMath::Sqrt(fSumXi2DeltaBi2)/fSumBiXi),2) + TMath::Power((fSumDeltaBi2/fSumBi),2) );
 
  cout<<"Analyzed events: "<<fAnalyzedEvents<<endl;
  cout<<"Width: "<<delta<<"\t Error: "<<deltaError<<endl;
  cout<<"New error: "<<deltaErrorNew<<endl;
  cout<<"Interval: "<<integral<<"\t Error: "<<integralError<<endl;
  cout<<"=================================================="<<endl;
}
  
//____________________________________________________________________//
TGraphErrors *AliBalance::DrawBalance() {
  // Draws the BF
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t xer[MAXIMUM_NUMBER_OF_STEPS];
  Double_t b[MAXIMUM_NUMBER_OF_STEPS];
  Double_t ber[MAXIMUM_NUMBER_OF_STEPS];

  if((fNp == 0)||(fNn == 0)) {
    cout<<"Couldn't find any particles in the analyzed interval!!!"<<endl;
    cout<<"Aborting....."<<endl;
    abort();
  }
  
  for(Int_t i = 0; i < fNumberOfBins; i++) {
    b[i] = GetBalance(i);
    ber[i] = GetError(i);
    //x[i] = fP2Start + fP2Step*i + fP2Step/2;
    x[i] = fP2Step*i + fP2Step/2;
    xer[i] = 0.0;
  }
  
  TGraphErrors *gr = new TGraphErrors(fNumberOfBins,x,b,xer,ber);
  gr->SetMarkerStyle(25);
  gr->GetXaxis()->SetTitleColor(1);
  if(fAnalysisType==0) {
    gr->GetXaxis()->SetTitle("#Delta y");
    gr->GetYaxis()->SetTitle("B(#Delta y)");
  }
  if(fAnalysisType==1) {
    gr->GetXaxis()->SetTitle("#Delta #eta");
    gr->GetYaxis()->SetTitle("B(#Delta #eta)");
  }
  if(fAnalysisType==2) {
    gr->GetXaxis()->SetTitle("q_{long} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{long}) [(GeV/c)^{-1}]");
  }
  if(fAnalysisType==3) {
    gr->GetXaxis()->SetTitle("q_{out} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{out}) [(GeV/c)^{-1}]");
  }
  if(fAnalysisType==4) {
    gr->GetXaxis()->SetTitle("q_{side} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{side}) [(GeV/c)^{-1}]");
  }
  if(fAnalysisType==5) {
    gr->GetXaxis()->SetTitle("q_{inv} (GeV/c)");
    gr->GetYaxis()->SetTitle("B(q_{inv}) [(GeV/c)^{-1}]");
  }
  if(fAnalysisType==6) {
    gr->GetXaxis()->SetTitle("#Delta #phi");
    gr->GetYaxis()->SetTitle("B(#Delta #phi)");
  }

  return gr;
}

//____________________________________________________________________//
void AliBalance::Merge(AliBalance *b) {
  //Merging function to be used for proof and grid
  fNp += b->GetNp();
  fNn += b->GetNn();
  for(Int_t i = 0; i < fNumberOfBins; i++) {
    fNnn[i] += b->GetNnn(i);
    fNpp[i] += b->GetNpp(i);
    fNpn[i] += b->GetNpn(i);
  }
}
