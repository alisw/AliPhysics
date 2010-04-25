//-*- Mode: C++ -*-
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Svein Lindal                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTCaloHistoInvMass
 * @author Svein Lindal <slindal@fys.uio.no>
 * @date 
 * @brief  Produces plots of invariant mass of two clusters. 
 */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloHistoInvMass.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "TObjArray.h"
#include "AliESDEvent.h"
#include "TRefArray.h"
#include "TH1F.h"
#include "TString.h"
#include "AliESDCaloCluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"

AliHLTCaloHistoInvMass::AliHLTCaloHistoInvMass(TString det) :
  fHistTwoClusterInvMass0(NULL),
  fHistTwoClusterInvMass1(NULL),
  fHistTwoClusterInvMass2(NULL),
  fHistTwoClusterInvMass3(NULL),
  fHistTwoClusterInvMass4(NULL)
{
  // See header file for documentation
  fHistTwoClusterInvMass0 = new TH1F(Form("%s fHistTwoClusterInvMass0", det.Data()), Form("Invariant mass of two clusters in %s, 0.8 GeV < E < 1.2", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass0->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass0->GetYaxis()->SetTitle("Counts");
  fHistTwoClusterInvMass0->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass0);

  fHistTwoClusterInvMass1 = new TH1F(Form("%s fHistTwoClusterInvMass1", det.Data()), Form("Invariant mass of two clusters in %s, 1.2 GeV < E < 1.6", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass1->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass1->GetYaxis()->SetTitle("Counts");
  fHistTwoClusterInvMass1->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass1);

  fHistTwoClusterInvMass2 = new TH1F(Form("%s fHistTwoClusterInvMass2", det.Data()), Form("Invariant mass of two clusters in %s, 1.6 GeV < E < 2.0", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass2->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass2->GetYaxis()->SetTitle("Counts");
  fHistTwoClusterInvMass2->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass2);

  fHistTwoClusterInvMass3 = new TH1F(Form("%s fHistTwoClusterInvMass3", det.Data()), Form("Invariant mass of two clusters in %s, 2.0 GeV < E < 4.0", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass3->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass3->GetYaxis()->SetTitle("Counts");
  fHistTwoClusterInvMass3->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass3);

  fHistTwoClusterInvMass4 = new TH1F(Form("%s fHistTwoClusterInvMass4", det.Data()), Form("Invariant mass of two clusters in %s E > 4.0 GeV", det.Data()), 200, 0, 1);
  fHistTwoClusterInvMass4->GetXaxis()->SetTitle("m_{#gamma#gamma} GeV");
  fHistTwoClusterInvMass4->GetYaxis()->SetTitle("Counts");
  fHistTwoClusterInvMass4->SetMarkerStyle(21);
  fHistArray->AddLast(fHistTwoClusterInvMass4);


}

AliHLTCaloHistoInvMass::~AliHLTCaloHistoInvMass()
{

  if(fHistTwoClusterInvMass0)
    delete fHistTwoClusterInvMass0;
  fHistTwoClusterInvMass0 = NULL;

  if(fHistTwoClusterInvMass1)
    delete fHistTwoClusterInvMass1;
  fHistTwoClusterInvMass1 = NULL;

  if(fHistTwoClusterInvMass2)
    delete fHistTwoClusterInvMass2;
  fHistTwoClusterInvMass2 = NULL;

  if(fHistTwoClusterInvMass3)
    delete fHistTwoClusterInvMass3;
  fHistTwoClusterInvMass3 = NULL;

  if(fHistTwoClusterInvMass4)
    delete fHistTwoClusterInvMass4;
  fHistTwoClusterInvMass4 = NULL;


}


Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, vector<AliHLTCaloClusterDataStruct*> &cVec) {
  //See header file for documentation
  
  Float_t cPos[nc][3];
  Float_t cEnergy[nc];

  for(int ic = 0; ic < nc; ic++) {
    AliHLTCaloClusterDataStruct * cluster = cVec.at(ic);
    cluster->GetPosition(cPos[ic]);
    cEnergy[ic] = cluster->E();
  }

  return FillInvariantMassHistograms(nc, cPos, cEnergy);
}



Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, TRefArray * clusterArray) {
  //See header file for documentation
  
  Float_t cPos[nc][3];
  Float_t cEnergy[nc];

  for(int ic = 0; ic < nc; ic++) {
    AliESDCaloCluster * cluster = static_cast<AliESDCaloCluster*>(clusterArray->At(ic));
    cluster->GetPosition(cPos[ic]);
    cEnergy[ic] = cluster->E();
  }

  return FillInvariantMassHistograms(nc, cPos, cEnergy);

}


Int_t AliHLTCaloHistoInvMass::FillInvariantMassHistograms(Int_t nc, Float_t cPos[][3], Float_t  cEnergy[]){ 

  Int_t iResult = 0;

  for(Int_t ic = 0; ic<(nc-1); ic++) { 
     
    //BALLE hardcoded variable
    if(cEnergy[ic] < 0.4)
      continue;

    //Get momentum vector
    TVector3 iVec(cPos[ic]);
    iVec = iVec.Unit();
    iVec = cEnergy[ic] * iVec;
    
    
    for(Int_t jc = ic+1; jc<nc; jc++) { 
     
    //BALLE hardcoded variable
      if(cEnergy[jc] < 0.4)
	continue;

      

      
      //Get second momentum vector
      TVector3 jVec(cPos[jc]);
      jVec = jVec.Unit();
      jVec = cEnergy[jc] * jVec;
      
      //Calculate inv mass
      Double_t m = TMath::Sqrt( 2 *(cEnergy[ic]* cEnergy[jc] - iVec.Dot(jVec) ) );
      
      //Fill histograms
      
      Float_t sum = cEnergy[ic]+cEnergy[ic];
      if(sum > 1.2)
      {
	 if(sum > 1.6)
	 {
	    if(sum > 2.0)
	    {
	       if(sum > 4.0)
	       {
		  fHistTwoClusterInvMass4->Fill(m);
	       }
	       else
	       {
		  fHistTwoClusterInvMass3->Fill(m);
	       }
	    }
	    else
	    {
	       fHistTwoClusterInvMass2->Fill(m);
	    }
	 }
	 else
	 {
	    fHistTwoClusterInvMass1->Fill(m);
	 }
      }
      else
      {
	 fHistTwoClusterInvMass0->Fill(m);
      }
    }
  }

  return iResult;

}

// template <class T>
// Int_t AliHLTCaloHistoInvMass::FillHistograms(Int_t nc, vector<T*> clusterVec ) {
//   Float_t cPos[nc][3];
//   Float_t cEnergy[nc];
  
//   for(int ic = 0; ic < nc; ic++) {
//     T * cluster = cVec.at(ic);
//     cluster->GetPosition(cPos[ic]);
//     cEnergy[ic] = cluster->E();
//   }
  

//   for(Int_t ic = 0; ic<(nc-1); ic++) { 

//     //Get momentum vector
//     TVector3 iVec(cPos[ic]);
//     iVec = iVec.Unit();
//     iVec = cEnergy[ic] * iVec;

    
//     for(Int_t jc = ic+1; jc<nc; jc++) { 

//       //Get second momentum vector
//       TVector3 jVec(cPos[jc]);
//       jVec = jVec.Unit();
//       jVec = cEnergy[jc] * jVec;

//       //Calculate inv mass
//       Double_t m = TMath::Sqrt( 2 *(cEnergy[ic]* cEnergy[jc] - iVec.Dot(jVec) ) );

//       //Fill histogram
//       fHistTwoClusterInvMass->Fill(m);
//     }
//   }
// }
