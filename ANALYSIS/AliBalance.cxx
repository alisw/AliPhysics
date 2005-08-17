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


#include <stdlib.h>

//ROOT
#include <Riostream.h>
#include <TROOT.h>
#include <TObject.h>
#include <TSystem.h>
#include <TObject.h>
#include <TMath.h>
#include <TLorentzVector.h>

#include "AliBalance.h"

ClassImp(AliBalance)

//----------------------------------------//
AliBalance::AliBalance()
{
  for(Int_t i = 0; i < MAXIMUM_NUMBER_OF_STEPS; i++)
    {
      fNpp[i] = .0;
      fNnn[i] = .0;
      fNpn[i] = .0;
      fB[i] = 0.0;
      ferror[i] = 0.0;
    } 
  fAnalyzedEvents = 0;
}

//----------------------------------------//
AliBalance::AliBalance(Double_t P2_Start, Double_t P2_Stop, Int_t P2_bins)
{
  this->fP2_Start = P2_Start;
  this->fP2_Stop = P2_Stop;
  this->fNumberOfBins = P2_bins;
  this->fP2_Step = fabs(fP2_Start - fP2_Stop) / (Double_t)fNumberOfBins;
}

//----------------------------------------//
AliBalance::~AliBalance()
{
}

//----------------------------------------//
void AliBalance::SetNumberOfBins(Int_t ibins)
{
  this->fNumberOfBins = ibins ;
}

//----------------------------------------//
void AliBalance::SetInterval(Double_t P2_Start, Double_t P2_Stop)
{
  this->fP2_Start = P2_Start;
  this->fP2_Stop = P2_Stop;
  this->fP2_Step = fabs(P2_Start - P2_Stop) / (Double_t)fNumberOfBins;
}

//----------------------------------------//
void AliBalance::SetAnalysisType(Int_t iType)
{
  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
  this->fAnalysisType = iType; 
  if(fAnalysisType==0)
    {
      cout<<" ====================== "<<endl;
      cout<<"||Analysis selected: y||"<<endl;
      cout<<" ====================== "<<endl;
    } 
  else if(fAnalysisType==1)
    {
      cout<<" ======================== "<<endl;
      cout<<"||Analysis selected: eta||"<<endl;
      cout<<" ======================== "<<endl;
    }
  else if(fAnalysisType==2)
    {
      cout<<" ========================== "<<endl;
      cout<<"||Analysis selected: Qlong||"<<endl;
      cout<<" ========================== "<<endl;
    }
  else if(fAnalysisType==3)
    {
      cout<<" ========================= "<<endl;
      cout<<"||Analysis selected: Qout||"<<endl;
      cout<<" ========================= "<<endl;
    }
  else if(fAnalysisType==4)
    {
      cout<<" ========================== "<<endl;
      cout<<"||Analysis selected: Qside||"<<endl;
      cout<<" ========================== "<<endl;
    }
  else if(fAnalysisType==5)
    {
      cout<<" ========================= "<<endl;
      cout<<"||Analysis selected: Qinv||"<<endl;
      cout<<" ========================= "<<endl;
    }
  else if(fAnalysisType==6)
    {
      cout<<" ======================== "<<endl;
      cout<<"||Analysis selected: phi||"<<endl;
      cout<<" ======================== "<<endl;
    }
  else
    {
      cout<<"Selection of analysis mode failed!!!"<<endl;
      cout<<"Choices are: 0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi"<<endl;
      abort();
    }
}

//----------------------------------------//
void AliBalance::SetParticles(TLorentzVector *P, Double_t *charge, Int_t dim)
{
  this->fV = P;
  this->fCharge = charge;
  fNtrack = dim;
}


//----------------------------------------//
void AliBalance::CalculateBalance()
{
  fAnalyzedEvents++;
  Int_t i = 0 , j = 0;
  Int_t ibin = 0;

  for(i = 0; i < fNtrack; i++)
    {
      if(fCharge[i] > 0)
	fNp += 1.;
      if(fCharge[i] < 0)
	fNn += 1.;
    }

  //0:y - 1:eta - 2:Qlong - 3:Qout - 4:Qside - 5:Qinv - 6:phi
   if(fAnalysisType==0)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t rap1 = 0.5*log((fV[i].E() - fV[i].Pz())/(fV[i].E() + fV[i].Pz())); 
	      Double_t rap2 = 0.5*log((fV[j].E() - fV[j].Pz())/(fV[j].E() + fV[j].Pz())); 
	      Double_t dy = TMath::Abs(rap1 - rap2);
	      ibin = Int_t(dy/fP2_Step);
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 0
   if(fAnalysisType==1)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t P1 = sqrt(pow(fV[i].Px(),2) + pow(fV[i].Py(),2) + pow(fV[i].Pz(),2)); 
	      Double_t P2 = sqrt(pow(fV[j].Px(),2) + pow(fV[j].Py(),2) + pow(fV[j].Pz(),2));
	      Double_t eta1 = 0.5*log((P1 - fV[i].Pz())/(P1 + fV[i].Pz())); 
	      Double_t eta2 = 0.5*log((P2 - fV[j].Pz())/(P2 + fV[j].Pz())); 
	      Double_t deta = TMath::Abs(eta1 - eta2);
	      ibin = Int_t(deta/fP2_Step);
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 1
   if(fAnalysisType==2)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t ETot = fV[i].E() + fV[j].E();
	      Double_t PxTot = fV[i].Px() + fV[j].Px();
	      Double_t PyTot = fV[i].Py() + fV[j].Py();
	      Double_t PzTot = fV[i].Pz() + fV[j].Pz();
	      Double_t Q0Tot = fV[i].E() - fV[j].E();
	      Double_t QzTot = fV[i].Pz() - fV[j].Pz();
	      Double_t Snn = pow(ETot,2) - pow(PxTot,2) - pow(PyTot,2) - pow(PzTot,2);
	      Double_t PtTot = sqrt( pow(PxTot,2) + pow(PyTot,2));
	      Double_t Qlong = TMath::Abs(ETot*QzTot - PzTot*Q0Tot)/sqrt(Snn + pow(PtTot,2));
	      ibin = Int_t(Qlong/fP2_Step);
	      //cout<<ibin<<endl;
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 2
  if(fAnalysisType==3)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t ETot = fV[i].E() + fV[j].E();
	      Double_t PxTot = fV[i].Px() + fV[j].Px();
	      Double_t PyTot = fV[i].Py() + fV[j].Py();
	      Double_t PzTot = fV[i].Pz() + fV[j].Pz();
	      Double_t QxTot = fV[i].Px() - fV[j].Px();
	      Double_t QyTot = fV[i].Py() - fV[j].Py();
	      Double_t Snn = pow(ETot,2) - pow(PxTot,2) - pow(PyTot,2) - pow(PzTot,2);
	      Double_t PtTot = sqrt( pow(PxTot,2) + pow(PyTot,2));
	      Double_t Qout = sqrt(Snn/(Snn + pow(PtTot,2))) * TMath::Abs(PxTot*QxTot + PyTot*QyTot)/PtTot;
	      ibin = Int_t(Qout/fP2_Step);
	      //cout<<ibin<<endl;
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 3
  if(fAnalysisType==4)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t PxTot = fV[i].Px() + fV[j].Px();
	      Double_t PyTot = fV[i].Py() + fV[j].Py();
	      Double_t QxTot = fV[i].Px() - fV[j].Px();
	      Double_t QyTot = fV[i].Py() - fV[j].Py();
	      Double_t PtTot = sqrt( pow(PxTot,2) + pow(PyTot,2));
	      Double_t Qside = TMath::Abs(PxTot*QyTot - PyTot*QxTot)/PtTot;
	      ibin = Int_t(Qside/fP2_Step);
	      //cout<<ibin<<endl;
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 4
   if(fAnalysisType==5)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t Q0Tot = fV[i].E() - fV[j].E();
	      Double_t QxTot = fV[i].Px() - fV[j].Px();
	      Double_t QyTot = fV[i].Py() - fV[j].Py();
	      Double_t QzTot = fV[i].Pz() - fV[j].Pz();
	      Double_t Qinv = sqrt(TMath::Abs(-pow(Q0Tot,2) +pow(QxTot,2) +pow(QyTot,2) +pow(QzTot,2)));
	      ibin = Int_t(Qinv/fP2_Step);
	      //cout<<ibin<<endl;
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 5	
  if(fAnalysisType==6)
    {
      for(i = 1; i < fNtrack; i++)
	{
	  for(j = 0; j < i; j++)
	    {
	      Double_t phi1 = TMath::ATan(fV[i].Py()/fV[i].Px())*180.0/TMath::Pi();
	      Double_t phi2 = TMath::ATan(fV[j].Py()/fV[j].Px())*180.0/TMath::Pi();
	      Double_t dphi = TMath::Abs(phi1 - phi2);
	      ibin = Int_t(dphi/fP2_Step);
	      if((fCharge[i] > 0)&&(fCharge[j] > 0))
		fNpp[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] < 0))
		fNnn[ibin] += 1.;
	      if((fCharge[i] > 0)&&(fCharge[j] < 0))
		fNpn[ibin] += 1.;
	      if((fCharge[i] < 0)&&(fCharge[j] > 0))
		fNpn[ibin] += 1.;
	    }
	}
    }//case 6
}

//----------------------------------------//
Double_t AliBalance::GetBalance(Int_t p2)
{
  fB[p2] = 0.5*(((fNpn[p2] - 2.0*fNnn[p2])/fNn) + ((fNpn[p2] - 2.0*fNpp[p2])/fNp))/fP2_Step;

  return fB[p2];
}
    
//----------------------------------------//
Double_t AliBalance::GetError(Int_t p2)
{
  ferror[p2] = sqrt( Double_t(fNpp[p2])/(Double_t(fNp)*Double_t(fNp)) + Double_t(fNnn[p2])/(Double_t(fNn)*Double_t(fNn)) + Double_t(fNpn[p2])*pow((0.5/Double_t(fNp) + 0.5/Double_t(fNn)),2))/fP2_Step;

  return ferror[p2];
}

//----------------------------------------//
void AliBalance::PrintResults()
{
  Double_t x[MAXIMUM_NUMBER_OF_STEPS];
  Double_t fSumXi = 0.0, fSumBi = 0.0, fSumBiXi = 0.0;
  Double_t fSumBiXi2 = 0.0, fSumBi2Xi2 = 0.0;
  Double_t fSumDeltaBi2 = 0.0, fSumXi2DeltaBi2 = 0.0;
  Double_t delta_bal_P2 = 0.0, Integral = 0.0;
  Double_t DeltaErrorNew = 0.0;

  cout<<"=================================================="<<endl;
  for(Int_t i = 0; i < fNumberOfBins; i++)
    { 
      Double_t x = fP2_Start + fP2_Step*i + fP2_Step/2 ;
      cout<<"B: "<<fB[i]<<"\t Error: "<<ferror[i]<<"\t bin: "<<x<<endl;
    } 
  cout<<"=================================================="<<endl;
  for(Int_t i = 1; i < fNumberOfBins; i++)
    {
      fSumXi += x[i];
      fSumBi += fB[i];
      fSumBiXi += fB[i]*x[i];
      fSumBiXi2 += fB[i]*pow(x[i],2);
      fSumBi2Xi2 += pow(fB[i],2)*pow(x[i],2);
      fSumDeltaBi2 +=  pow(ferror[i],2) ;
      fSumXi2DeltaBi2 += pow(x[i],2) * pow(ferror[i],2) ;
      
      delta_bal_P2 += fP2_Step*pow(ferror[i],2) ;
      Integral += fP2_Step*fB[i] ;
    }
  for(Int_t i = 1; i < fNumberOfBins; i++)
    {
      DeltaErrorNew += ferror[i]*(x[i]*fSumBi - fSumBiXi)/pow(fSumBi,2);
    }
  Double_t IntegralError = sqrt(delta_bal_P2) ;
  
  Double_t Delta = fSumBiXi / fSumBi ;
  Double_t DeltaError = (fSumBiXi / fSumBi) * sqrt(pow((sqrt(fSumXi2DeltaBi2)/fSumBiXi),2) + pow((fSumDeltaBi2/fSumBi),2) ) ;
 
  cout<<"Analyzed events: "<<fAnalyzedEvents<<endl;
  cout<<"Width: "<<Delta<<"\t Error: "<<DeltaError<<endl;
  cout<<"New error: "<<DeltaErrorNew<<endl;
  cout<<"Interval: "<<Integral<<"\t Error: "<<IntegralError<<endl;
  cout<<"=================================================="<<endl;
}
  
