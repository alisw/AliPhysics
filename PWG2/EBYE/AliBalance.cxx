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
#include <TLorentzVector.h>
#include <TGraphErrors.h>

#include "AliBalance.h"

ClassImp(AliBalance)

//____________________________________________________________________//
AliBalance::AliBalance(const char* name, const char* title) :TNamed(name,title) {
  // Default constructor
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 
  fNp = 0.0;
  fNn = 0.0;
  fP2Start = 0.0;
  fP2Stop = 0.0;
  fP2Step = 0.0; 
  fAnalysisType = 0;
  fNumberOfBins = 0;
  fNtrack = 0;
  fAnalyzedEvents = 0;
}

//____________________________________________________________________//
AliBalance::AliBalance(const char* name, const char* title, Double_t p2Start, Double_t p2Stop, Int_t p2Bins) :TNamed(name,title) {
  // Constructor
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 
  fNp = 0.0;
  fNn = 0.0;
  fAnalysisType = 0;
  fNtrack = 0;
  fAnalyzedEvents = 0;

  fP2Start = p2Start;
  fP2Stop = p2Stop;
  fNumberOfBins = p2Bins;
  fP2Step = TMath::Abs(fP2Start - fP2Stop) / (Double_t)fNumberOfBins;
}

//____________________________________________________________________//
AliBalance::AliBalance(const AliBalance& balance):
  TNamed(balance),
  fNtrack(balance.fNtrack),
  fNumberOfBins(balance.fNumberOfBins),
  fAnalysisType(balance.fAnalysisType),
  fAnalyzedEvents(balance.fAnalyzedEvents),
  fP2Start(balance.fP2Start),
  fP2Stop(balance.fP2Stop),
  fP2Step(balance.fP2Step),
  fNn(balance.fNn),
  fNp(balance.fNp),
  fCharge(0),
  fV(0)
{
  //copy constructor

  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++) {
    fNpp[i] = .0;
    fNnn[i] = .0;
    fNpn[i] = .0;
    fB[i] = 0.0;
    ferror[i] = 0.0;
  } 

  if (balance.fV) fV = new TLorentzVector(*(balance.fV));
  if (balance.fCharge) fCharge = new Double_t(*(balance.fCharge));
}

//____________________________________________________________________//
AliBalance::~AliBalance() {
  // Destructor
  delete fV;
  delete fCharge;
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
void AliBalance::SetParticles(TLorentzVector *P, Double_t *charge, Int_t dim) {
  // Sets a new particle with given 4-momentum and charge.
  // dim is the size of the array of charges and corresponds
  // to the number of selected tracks.
  this->fV = P;
  this->fCharge = charge;
  fNtrack = dim;
}


//____________________________________________________________________//
void AliBalance::CalculateBalance() {
  // Calculates the balance function
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t ibin = 0;

  for(i = 0; i < fNtrack; i++) {
    if(fCharge[i] > 0) fNp += 1.;
    if(fCharge[i] < 0) fNn += 1.;
  }

  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  if(fAnalysisType==0) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t rap1 = 0.5*log((fV[i].E() + fV[i].Pz())/(fV[i].E() - fV[i].Pz())); 
	Double_t rap2 = 0.5*log((fV[j].E() + fV[j].Pz())/(fV[j].E() - fV[j].Pz())); 
	Double_t dy = TMath::Abs(rap1 - rap2);
	ibin = Int_t(dy/fP2Step);
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 0
  if(fAnalysisType==1) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t p1 = sqrt(pow(fV[i].Px(),2) + pow(fV[i].Py(),2) + pow(fV[i].Pz(),2)); 
	Double_t p2 = sqrt(pow(fV[j].Px(),2) + pow(fV[j].Py(),2) + pow(fV[j].Pz(),2));
	Double_t eta1 = 0.5*log((p1 + fV[i].Pz())/(p1 - fV[i].Pz())); 
	Double_t eta2 = 0.5*log((p2 + fV[j].Pz())/(p2 - fV[j].Pz())); 
	Double_t deta = TMath::Abs(eta1 - eta2);
	ibin = Int_t(deta/fP2Step);
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 1
  if(fAnalysisType==2) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t eTot = fV[i].E() + fV[j].E();
	Double_t pxTot = fV[i].Px() + fV[j].Px();
	Double_t pyTot = fV[i].Py() + fV[j].Py();
	Double_t pzTot = fV[i].Pz() + fV[j].Pz();
	Double_t q0Tot = fV[i].E() - fV[j].E();
	Double_t qzTot = fV[i].Pz() - fV[j].Pz();
	Double_t snn = pow(eTot,2) - pow(pxTot,2) - pow(pyTot,2) - pow(pzTot,2);
	Double_t ptTot = sqrt( pow(pxTot,2) + pow(pyTot,2));
	Double_t qLong = TMath::Abs(eTot*qzTot - pzTot*q0Tot)/sqrt(snn + pow(ptTot,2));
	ibin = Int_t(qLong/fP2Step);
	//cout<<ibin<<endl;
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 2
  if(fAnalysisType==3) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t eTot = fV[i].E() + fV[j].E();
	Double_t pxTot = fV[i].Px() + fV[j].Px();
	Double_t pyTot = fV[i].Py() + fV[j].Py();
	Double_t pzTot = fV[i].Pz() + fV[j].Pz();
	Double_t qxTot = fV[i].Px() - fV[j].Px();
	Double_t qyTot = fV[i].Py() - fV[j].Py();
	Double_t snn = pow(eTot,2) - pow(pxTot,2) - pow(pyTot,2) - pow(pzTot,2);
	Double_t ptTot = sqrt( pow(pxTot,2) + pow(pyTot,2));
	Double_t qOut = sqrt(snn/(snn + pow(ptTot,2))) * TMath::Abs(pxTot*qxTot + pyTot*qyTot)/ptTot;
	ibin = Int_t(qOut/fP2Step);
	//cout<<ibin<<endl;
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 3
  if(fAnalysisType==4) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t pxTot = fV[i].Px() + fV[j].Px();
	Double_t pyTot = fV[i].Py() + fV[j].Py();
	Double_t qxTot = fV[i].Px() - fV[j].Px();
	Double_t qyTot = fV[i].Py() - fV[j].Py();
	Double_t ptTot = sqrt( pow(pxTot,2) + pow(pyTot,2));
	Double_t qSide = TMath::Abs(pxTot*qyTot - pyTot*qxTot)/ptTot;
	ibin = Int_t(qSide/fP2Step);
	//cout<<ibin<<endl;
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 4
  if(fAnalysisType==5) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t q0Tot = fV[i].E() - fV[j].E();
	Double_t qxTot = fV[i].Px() - fV[j].Px();
	Double_t qyTot = fV[i].Py() - fV[j].Py();
	Double_t qzTot = fV[i].Pz() - fV[j].Pz();
	Double_t qInv = sqrt(TMath::Abs(-pow(q0Tot,2) +pow(qxTot,2) +pow(qyTot,2) +pow(qzTot,2)));
	ibin = Int_t(qInv/fP2Step);
	//cout<<ibin<<endl;
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 5	
  if(fAnalysisType==6) {
    for(i = 1; i < fNtrack; i++) {
      for(j = 0; j < i; j++) {
	Double_t phi1 = TMath::ATan(fV[i].Py()/fV[i].Px())*180.0/TMath::Pi();
	Double_t phi2 = TMath::ATan(fV[j].Py()/fV[j].Px())*180.0/TMath::Pi();
	Double_t dphi = TMath::Abs(phi1 - phi2);
	ibin = Int_t(dphi/fP2Step);
	if((fCharge[i] > 0)&&(fCharge[j] > 0)) fNpp[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] < 0)) fNnn[ibin] += 1.;
	if((fCharge[i] > 0)&&(fCharge[j] < 0)) fNpn[ibin] += 1.;
	if((fCharge[i] < 0)&&(fCharge[j] > 0)) fNpn[ibin] += 1.;
      }
    }
  }//case 6
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
  ferror[p2] = sqrt( Double_t(fNpp[p2])/(Double_t(fNp)*Double_t(fNp)) + Double_t(fNnn[p2])/(Double_t(fNn)*Double_t(fNn)) + Double_t(fNpn[p2])*pow((0.5/Double_t(fNp) + 0.5/Double_t(fNn)),2))/fP2Step;

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
    x[i] = fP2Start + fP2Step*i + fP2Step/2;
    cout<<"B: "<<fB[i]<<"\t Error: "<<ferror[i]<<"\t bin: "<<x[i]<<endl;
  } 
  cout<<"=================================================="<<endl;
  for(Int_t i = 1; i < fNumberOfBins; i++) {
    fSumXi += x[i];
    fSumBi += fB[i];
    fSumBiXi += fB[i]*x[i];
    fSumBiXi2 += fB[i]*pow(x[i],2);
    fSumBi2Xi2 += pow(fB[i],2)*pow(x[i],2);
    fSumDeltaBi2 +=  pow(ferror[i],2);
    fSumXi2DeltaBi2 += pow(x[i],2) * pow(ferror[i],2);
    
    deltaBalP2 += fP2Step*pow(ferror[i],2);
    integral += fP2Step*fB[i];
  }
  for(Int_t i = 1; i < fNumberOfBins; i++) deltaErrorNew += ferror[i]*(x[i]*fSumBi - fSumBiXi)/pow(fSumBi,2);
   
  Double_t integralError = sqrt(deltaBalP2);
  
  Double_t delta = fSumBiXi / fSumBi;
  Double_t deltaError = (fSumBiXi / fSumBi) * sqrt(pow((sqrt(fSumXi2DeltaBi2)/fSumBiXi),2) + pow((fSumDeltaBi2/fSumBi),2) );
 
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
    x[i] = fP2Start + fP2Step*i + fP2Step/2;
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
    gr->GetXaxis()->SetTitle("Q_{long} [GeV]");
    gr->GetYaxis()->SetTitle("B(Q_{long})");
  }
  if(fAnalysisType==3) {
    gr->GetXaxis()->SetTitle("Q_{out} [GeV]");
    gr->GetYaxis()->SetTitle("B(Q_{out})");
  }
  if(fAnalysisType==4) {
    gr->GetXaxis()->SetTitle("Q_{side} [GeV]");
    gr->GetYaxis()->SetTitle("B(Q_{side})");
  }
  if(fAnalysisType==5) {
    gr->GetXaxis()->SetTitle("Q_{inv} [GeV]");
    gr->GetYaxis()->SetTitle("B(Q_{inv})");
  }
  if(fAnalysisType==6) {
    gr->GetXaxis()->SetTitle("#Delta #phi");
    gr->GetYaxis()->SetTitle("B(#Delta #phi)");
  }

  return gr;
}
