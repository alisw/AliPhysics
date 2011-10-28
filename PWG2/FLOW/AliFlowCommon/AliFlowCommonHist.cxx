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
#include "TList.h"
#include "TH2F.h"
#include "AliFlowVector.h"
#include "TBrowser.h"

class TH1F;
class TH1D;

// AliFlowCommonHist:
//
// Description: Class to organise common histograms for Flow Analysis
//
// authors: N. van der Kolk (kolk@nikhef.nl), A. Bilandzic (anteb@nikhef.nl), RS


ClassImp(AliFlowCommonHist)

//-----------------------------------------------------------------------

AliFlowCommonHist::AliFlowCommonHist():
  TNamed(),
  fBookOnlyBasic(kFALSE),
  fHistMultRP(NULL),
  fHistMultPOI(NULL),
  fHistMultPOIvsRP(NULL),
  fHistPtRP(NULL),
  fHistPtPOI(NULL),
  fHistPtSub0(NULL),
  fHistPtSub1(NULL),
  fHistPhiRP(NULL),
  fHistPhiPOI(NULL),
  fHistPhiSub0(NULL),
  fHistPhiSub1(NULL),
  fHistEtaRP(NULL),
  fHistEtaPOI(NULL),
  fHistEtaSub0(NULL),
  fHistEtaSub1(NULL),
  fHistPhiEtaRP(NULL),
  fHistPhiEtaPOI(NULL),
  fHistProMeanPtperBin(NULL),
  fHistWeightvsPhi(NULL),
  fHistQ(NULL),
  fHistAngleQ(NULL),
  fHistAngleQSub0(NULL),
  fHistAngleQSub1(NULL), 
  fHarmonic(NULL),
  fRefMultVsNoOfRPs(NULL),
  fHistRefMult(NULL),
  fHistList(NULL)
{
  
  //default constructor
  
}

AliFlowCommonHist::AliFlowCommonHist(const AliFlowCommonHist& a):
  TNamed(),
  fBookOnlyBasic(a.fBookOnlyBasic),
  fHistMultRP(new TH1F(*a.fHistMultRP)),
  fHistMultPOI(new TH1F(*a.fHistMultPOI)),
  fHistMultPOIvsRP(new TH2F(*a.fHistMultPOIvsRP)),
  fHistPtRP(new TH1F(*a.fHistPtRP)),
  fHistPtPOI(new TH1F(*a.fHistPtPOI)),
  fHistPtSub0(new TH1F(*a.fHistPtSub0)),
  fHistPtSub1(new TH1F(*a.fHistPtSub1)),
  fHistPhiRP(new TH1F(*a.fHistPhiRP)),
  fHistPhiPOI(new TH1F(*a.fHistPhiPOI)),
  fHistPhiSub0(new TH1F(*a.fHistPhiSub0)),
  fHistPhiSub1(new TH1F(*a.fHistPhiSub1)),
  fHistEtaRP(new TH1F(*a.fHistEtaRP)),
  fHistEtaPOI(new TH1F(*a.fHistEtaPOI)),
  fHistEtaSub0(new TH1F(*a.fHistEtaSub0)),
  fHistEtaSub1(new TH1F(*a.fHistEtaSub1)),
  fHistPhiEtaRP(new TH2F(*a.fHistPhiEtaRP)),
  fHistPhiEtaPOI(new TH2F(*a.fHistPhiEtaPOI)),
  fHistProMeanPtperBin(new TProfile(*a.fHistProMeanPtperBin)),
  fHistWeightvsPhi(new TH2F(*a.fHistWeightvsPhi)),
  fHistQ(new TH1F(*a.fHistQ)),
  fHistAngleQ(new TH1F(*a.fHistAngleQ)),
  fHistAngleQSub0(new TH1F(*a.fHistAngleQSub0)),
  fHistAngleQSub1(new TH1F(*a.fHistAngleQSub1)), 
  fHarmonic(new TProfile(*a.fHarmonic)),
  fRefMultVsNoOfRPs(new TProfile(*a.fRefMultVsNoOfRPs)),
  fHistRefMult(new TH1F(*a.fHistRefMult)),  
  fHistList(NULL)
{
  // copy constructor

  fHistList = new TList();
  fHistList-> Add(fHistMultRP);        
  fHistList-> Add(fHistMultPOI);
  if(!fBookOnlyBasic){fHistList-> Add(fHistMultPOIvsRP);}
  fHistList-> Add(fHistPtRP);          
  fHistList-> Add(fHistPtPOI);
  if(!fBookOnlyBasic){fHistList-> Add(fHistPtSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPtSub1);}
  fHistList-> Add(fHistPhiRP);          
  fHistList-> Add(fHistPhiPOI);
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiSub1);}    
  fHistList-> Add(fHistEtaRP);          
  fHistList-> Add(fHistEtaPOI); 
  if(!fBookOnlyBasic){fHistList-> Add(fHistEtaSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistEtaSub1);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiEtaRP);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiEtaPOI);}
  fHistList-> Add(fHistProMeanPtperBin);
  if(!fBookOnlyBasic){fHistList-> Add(fHistWeightvsPhi);}
  fHistList-> Add(fHarmonic);  
  fHistList-> Add(fRefMultVsNoOfRPs);
  fHistList-> Add(fHistRefMult); 
  if(!fBookOnlyBasic){fHistList-> Add(fHistQ);} 
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQ);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQSub1);}
  //  TListIter next = TListIter(a.fHistList);

}


//-----------------------------------------------------------------------

  AliFlowCommonHist::AliFlowCommonHist(const char *anInput, const char *title, Bool_t bookOnlyBasic):
    TNamed(anInput,title),
    fBookOnlyBasic(bookOnlyBasic),
    fHistMultRP(NULL),
    fHistMultPOI(NULL),
    fHistMultPOIvsRP(NULL),
    fHistPtRP(NULL),
    fHistPtPOI(NULL),
    fHistPtSub0(NULL),
    fHistPtSub1(NULL),
    fHistPhiRP(NULL),
    fHistPhiPOI(NULL),
    fHistPhiSub0(NULL),
    fHistPhiSub1(NULL),
    fHistEtaRP(NULL),
    fHistEtaPOI(NULL),
    fHistEtaSub0(NULL),
    fHistEtaSub1(NULL),
    fHistPhiEtaRP(NULL),
    fHistPhiEtaPOI(NULL),
    fHistProMeanPtperBin(NULL),
    fHistWeightvsPhi(NULL),
    fHistQ(NULL),
    fHistAngleQ(NULL),
    fHistAngleQSub0(NULL),
    fHistAngleQSub1(NULL), 
    fHarmonic(NULL),
    fRefMultVsNoOfRPs(NULL),
    fHistRefMult(NULL),
    fHistList(NULL)
{

  //constructor creating histograms 
  Int_t iNbinsMult = AliFlowCommonConstants::GetMaster()->GetNbinsMult();
  Int_t iNbinsPt = AliFlowCommonConstants::GetMaster()->GetNbinsPt();
  Int_t iNbinsPhi = AliFlowCommonConstants::GetMaster()->GetNbinsPhi();
  Int_t iNbinsEta = AliFlowCommonConstants::GetMaster()->GetNbinsEta();
  Int_t iNbinsQ = AliFlowCommonConstants::GetMaster()->GetNbinsQ();
  TString sName;

  Double_t  dMultMin = AliFlowCommonConstants::GetMaster()->GetMultMin();            
  Double_t  dMultMax = AliFlowCommonConstants::GetMaster()->GetMultMax();
  Double_t  dPtMin = AliFlowCommonConstants::GetMaster()->GetPtMin();	     
  Double_t  dPtMax = AliFlowCommonConstants::GetMaster()->GetPtMax();
  Double_t  dPhiMin = AliFlowCommonConstants::GetMaster()->GetPhiMin();	     
  Double_t  dPhiMax = AliFlowCommonConstants::GetMaster()->GetPhiMax();
  Double_t  dEtaMin = AliFlowCommonConstants::GetMaster()->GetEtaMin();	     
  Double_t  dEtaMax = AliFlowCommonConstants::GetMaster()->GetEtaMax();	     
  Double_t  dQMin = AliFlowCommonConstants::GetMaster()->GetQMin();	     
  Double_t  dQMax = AliFlowCommonConstants::GetMaster()->GetQMax();	
  Double_t  dHistWeightvsPhiMin = AliFlowCommonConstants::GetMaster()->GetHistWeightvsPhiMin();
  Double_t  dHistWeightvsPhiMax = AliFlowCommonConstants::GetMaster()->GetHistWeightvsPhiMax();
  
  cout<<"The settings for the common histograms are as follows:"<<endl;
  cout<<"Multiplicity: "<<iNbinsMult<<" bins between "<<dMultMin<<" and "<<dMultMax<<endl;
  cout<<"Pt: "<<iNbinsPt<<" bins between "<<dPtMin<<" and "<<dPtMax<<endl;
  cout<<"Phi: "<<iNbinsPhi<<" bins between "<<dPhiMin<<" and "<<dPhiMax<<endl;
  cout<<"Eta: "<<iNbinsEta<<" bins between "<<dEtaMin<<" and "<<dEtaMax<<endl;
  cout<<"Q: "<<iNbinsQ<<" bins between "<<dQMin<<" and "<<dQMax<<endl;

  //Multiplicity
  sName = "Control_Flow_MultRP_";
  sName +=anInput;
  fHistMultRP = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultRP ->SetXTitle("Multiplicity for RP selection");
  fHistMultRP ->SetYTitle("Counts");

  sName = "Control_Flow_MultPOI_";
  sName +=anInput;
  fHistMultPOI = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistMultPOI ->SetXTitle("Multiplicity for POI selection");
  fHistMultPOI ->SetYTitle("Counts");

  if(!fBookOnlyBasic){
  sName = "Control_Flow_MultPOIvsRP_";
  sName +=anInput;
  fHistMultPOIvsRP = new TH2F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax,100, dMultMin, dMultMax);
  fHistMultPOIvsRP->SetXTitle("Multiplicity for RP selection");
  fHistMultPOIvsRP->SetYTitle("Multiplicity for POI selection");
  }
  
  //Pt
  sName = "Control_Flow_PtRP_";
  sName +=anInput;
  fHistPtRP = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  fHistPtRP ->SetXTitle("P_{t} (GeV/c) for RP selection");
  fHistPtRP ->SetYTitle("Counts");

  sName = "Control_Flow_PtPOI_";
  sName +=anInput;
  fHistPtPOI = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  //binning has to be the same as for fHistProVPt! use to get Nprime!
  fHistPtPOI ->SetXTitle("P_{t} (GeV/c) for POI selection");
  fHistPtPOI ->SetYTitle("Counts");

  if(!fBookOnlyBasic){
  sName = "Control_Flow_PtSub0_";
  sName +=anInput;
  fHistPtSub0 = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  fHistPtSub0 ->SetXTitle("P_{t} (GeV/c) for Subevent 0 selection");
  fHistPtSub0 ->SetYTitle("Counts");
  }
  
  if(!fBookOnlyBasic){
  sName = "Control_Flow_PtSub1_";
  sName +=anInput;
  fHistPtSub1 = new TH1F(sName.Data(), sName.Data(),iNbinsPt, dPtMin, dPtMax); 
  fHistPtSub1 ->SetXTitle("P_{t} (GeV/c) for Subevent 1 selection");
  fHistPtSub1 ->SetYTitle("Counts");
  }
  
  //Phi
  sName = "Control_Flow_PhiRP_";
  sName +=anInput;
  fHistPhiRP = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiRP ->SetXTitle("#phi for RP selection");
  fHistPhiRP ->SetYTitle("Counts");

  sName = "Control_Flow_PhiPOI_";
  sName +=anInput;
  fHistPhiPOI = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiPOI ->SetXTitle("#phi for POI selection");
  fHistPhiPOI ->SetYTitle("Counts");

  if(!fBookOnlyBasic){
  sName = "Control_Flow_PhiSub0_";
  sName +=anInput;
  fHistPhiSub0 = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiSub0 ->SetXTitle("#phi for Subevent 0 selection");
  fHistPhiSub0 ->SetYTitle("Counts");
  }
   
  if(!fBookOnlyBasic){
  sName = "Control_Flow_PhiSub1_";
  sName +=anInput;
  fHistPhiSub1 = new TH1F(sName.Data(), sName.Data(),iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiSub1 ->SetXTitle("#phi for Subevent 1 selection");
  fHistPhiSub1 ->SetYTitle("Counts");
  }
  
  //Eta
  sName = "Control_Flow_EtaRP_";
  sName +=anInput;
  fHistEtaRP = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaRP ->SetXTitle("#eta for RP selection");
  fHistEtaRP ->SetYTitle("Counts");

  sName = "Control_Flow_EtaPOI_";
  sName +=anInput;
  fHistEtaPOI = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaPOI ->SetXTitle("#eta for POI selection");
  fHistEtaPOI ->SetYTitle("Counts");

  if(!fBookOnlyBasic){
  sName = "Control_Flow_EtaSub0_";
  sName +=anInput;
  fHistEtaSub0 = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaSub0 ->SetXTitle("#eta for Subevent 0 selection");
  fHistEtaSub0 ->SetYTitle("Counts");
  }
  
  if(!fBookOnlyBasic){
  sName = "Control_Flow_EtaSub1_";
  sName +=anInput;
  fHistEtaSub1 = new TH1F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax);
  fHistEtaSub1 ->SetXTitle("#eta for Subevent 1 selection");
  fHistEtaSub1 ->SetYTitle("Counts");
  }

  if(!fBookOnlyBasic){
  //Phi vs Eta
  sName = "Control_Flow_PhiEtaRP_";
  sName +=anInput;
  fHistPhiEtaRP = new TH2F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax, iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiEtaRP ->SetXTitle("#eta");
  fHistPhiEtaRP ->SetYTitle("#phi");
  }
  
  if(!fBookOnlyBasic){
  sName = "Control_Flow_PhiEtaPOI_";
  sName +=anInput;
  fHistPhiEtaPOI = new TH2F(sName.Data(), sName.Data(),iNbinsEta, dEtaMin, dEtaMax, iNbinsPhi, dPhiMin, dPhiMax);
  fHistPhiEtaPOI ->SetXTitle("#eta");
  fHistPhiEtaPOI ->SetYTitle("#phi");
  }
  
  //Mean Pt per pt bin 
  sName = "Control_FlowPro_MeanPtperBin_";
  sName +=anInput;
  fHistProMeanPtperBin = new TProfile(sName.Data(), sName.Data(),iNbinsPt,dPtMin,dPtMax);
  fHistProMeanPtperBin ->SetXTitle("P_{t} (GeV/c) ");
  fHistProMeanPtperBin ->SetYTitle("<P_{t}> (GeV/c) ");

  
  if(!fBookOnlyBasic){
  //Particle weight
  sName = "Control_Flow_WeightvsPhi_";
  sName +=anInput;
  fHistWeightvsPhi = new TH2F(sName.Data(), sName.Data(), iNbinsPhi, dPhiMin, dPhiMax, 300, dHistWeightvsPhiMin, dHistWeightvsPhiMax); 
  fHistWeightvsPhi ->SetXTitle("#phi");
  fHistWeightvsPhi ->SetYTitle("weight");
  }
  
  if(!fBookOnlyBasic){
  //Q vector
  sName = "Control_Flow_Q_";
  sName +=anInput;
  fHistQ = new TH1F(sName.Data(), sName.Data(),iNbinsQ, dQMin, dQMax);
  fHistQ ->SetXTitle("Q_{vector}/Mult");
  fHistQ ->SetYTitle("Counts");  
  }
    
  if(!fBookOnlyBasic){
  //Angle of Q vector
  sName = "Control_Flow_AngleQ_";
  sName +=anInput;
  fHistAngleQ = new TH1F(sName.Data(), sName.Data(),72, 0., TMath::Pi());
  fHistAngleQ ->SetXTitle("Angle of Q_{vector}");
  fHistAngleQ ->SetYTitle("Counts"); 
  }
   
  if(!fBookOnlyBasic){
  sName = "Control_Flow_AngleQSub0_";
  sName +=anInput;
  fHistAngleQSub0 = new TH1F(sName.Data(), sName.Data(),72, 0., TMath::Pi());
  fHistAngleQSub0 ->SetXTitle("Angle of Q_{vector} for Subevent 0");
  fHistAngleQSub0 ->SetYTitle("Counts"); 
  }
  
  if(!fBookOnlyBasic){
  sName = "Control_Flow_AngleQSub1_";
  sName +=anInput;
  fHistAngleQSub1 = new TH1F(sName.Data(), sName.Data(),72, 0., TMath::Pi());
  fHistAngleQSub1 ->SetXTitle("Angle of Q_{vector} for Subevent 1");
  fHistAngleQSub1 ->SetYTitle("Counts"); 
  }
  
  //harmonic
  sName = "Control_Flow_Harmonic_";
  sName +=anInput;
  fHarmonic = new TProfile(sName.Data(),sName.Data(),1,0,1);
  fHarmonic ->SetYTitle("harmonic");
  
  //<reference multiplicity> versus # of RPs
  sName = "Reference_Multiplicity_Vs_Number_Of_RPs_";
  sName +=anInput;
  fRefMultVsNoOfRPs = new TProfile(sName.Data(),sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fRefMultVsNoOfRPs->SetYTitle("<reference multiplicity>");
  fRefMultVsNoOfRPs->SetXTitle("# of RPs");

  //reference multiplicity
  sName = "Control_Flow_Ref_Mult_";
  sName +=anInput;
  fHistRefMult = new TH1F(sName.Data(), sName.Data(),iNbinsMult, dMultMin, dMultMax);
  fHistRefMult->SetXTitle("Reference multiplicity");
  fHistRefMult->SetYTitle("Counts");

  //list of histograms if added here also add in copy constructor
  fHistList = new TList();
  fHistList-> Add(fHistMultRP);        
  fHistList-> Add(fHistMultPOI); 
  if(!fBookOnlyBasic){fHistList-> Add(fHistMultPOIvsRP);}
  fHistList-> Add(fHistPtRP);          
  fHistList-> Add(fHistPtPOI); 
  if(!fBookOnlyBasic){fHistList-> Add(fHistPtSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPtSub1);}
  fHistList-> Add(fHistPhiRP);          
  fHistList-> Add(fHistPhiPOI);
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiSub1);}
  fHistList-> Add(fHistEtaRP);          
  fHistList-> Add(fHistEtaPOI); 
  if(!fBookOnlyBasic){fHistList-> Add(fHistEtaSub0);} 
  if(!fBookOnlyBasic){fHistList-> Add(fHistEtaSub1);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiEtaRP);}  
  if(!fBookOnlyBasic){fHistList-> Add(fHistPhiEtaPOI);}
  fHistList-> Add(fHistProMeanPtperBin);
  if(!fBookOnlyBasic){fHistList-> Add(fHistWeightvsPhi);}
  fHistList-> Add(fHarmonic); 
  fHistList-> Add(fRefMultVsNoOfRPs); 
  fHistList-> Add(fHistRefMult);   
  if(!fBookOnlyBasic){fHistList-> Add(fHistQ);}           
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQ);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQSub0);}
  if(!fBookOnlyBasic){fHistList-> Add(fHistAngleQSub1);} 

}


//----------------------------------------------------------------------- 

AliFlowCommonHist::~AliFlowCommonHist()
{
  //deletes histograms
  delete fHistMultRP;       
  delete fHistMultPOI; 
  if(!fBookOnlyBasic){delete fHistMultPOIvsRP;}
  delete fHistPtRP;         
  delete fHistPtPOI;
  if(!fBookOnlyBasic){delete fHistPtSub0;}
  if(!fBookOnlyBasic){delete fHistPtSub1;}
  delete fHistPhiRP;        
  delete fHistPhiPOI;
  if(!fBookOnlyBasic){delete fHistPhiSub0;}
  if(!fBookOnlyBasic){delete fHistPhiSub1;}
  delete fHistEtaRP;        
  delete fHistEtaPOI;
  if(!fBookOnlyBasic){delete fHistEtaSub0;}
  if(!fBookOnlyBasic){delete fHistEtaSub1;}
  delete fHistPhiEtaRP;
  delete fHistPhiEtaPOI;
  delete fHistProMeanPtperBin;
  if(!fBookOnlyBasic){delete fHistWeightvsPhi;}
  if(!fBookOnlyBasic){delete fHistQ;}
  if(!fBookOnlyBasic){delete fHistAngleQ;}
  if(!fBookOnlyBasic){delete fHistAngleQSub0;}
  if(!fBookOnlyBasic){delete fHistAngleQSub1;}
  delete fHarmonic;
  delete fRefMultVsNoOfRPs;
  delete fHistRefMult;  
  delete fHistList;
}


//----------------------------------------------------------------------- 

Bool_t AliFlowCommonHist::FillControlHistograms(AliFlowEventSimple* anEvent,TList *weightsList, Bool_t usePhiWeights, Bool_t usePtWeights, Bool_t useEtaWeights)
{
  //Fills the control histograms
  if (!anEvent){
    cout<<"##### FillControlHistograms: FlowEvent pointer null"<<endl;
    return kFALSE;
  }

  //track datamembers
  Double_t dPt     = 0.;
  Double_t dPhi    = 0.;
  Double_t dEta    = 0.;
  Double_t dWeight = 1.;

  //weights used for corrections
  Double_t dWPhi = 1.;
  Double_t dWPt  = 1.;
  Double_t dWEta = 1.;

  TH1F *phiWeights     = NULL;
  TH1F *phiWeightsSub0 = NULL;
  TH1F *phiWeightsSub1 = NULL;
  TH1D *ptWeights      = NULL;
  TH1D *etaWeights     = NULL;

  Int_t nBinsPhi     = 0;
  Int_t nBinsPhiSub0 = 0;
  Int_t nBinsPhiSub1 = 0;
  Double_t dBinWidthPt  = 0.;
  Double_t dPtMin       = 0.;
  Double_t dBinWidthEta = 0.;
  Double_t dEtaMin      = 0.;

  if(weightsList)
    {
      if(usePhiWeights)
	{
	  phiWeights = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights"));
	  if(phiWeights) nBinsPhi = phiWeights->GetNbinsX();
	  phiWeightsSub0 = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights_sub0"));
	  if(phiWeightsSub0) nBinsPhiSub0 = phiWeightsSub0->GetNbinsX();
	  phiWeightsSub1 = dynamic_cast<TH1F *>(weightsList->FindObject("phi_weights_sub1"));
	  if(phiWeightsSub1) nBinsPhiSub1 = phiWeightsSub1->GetNbinsX();
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


  
  //fill the histograms
  AliFlowVector vQ = anEvent->GetQ(2, weightsList, usePhiWeights, usePtWeights, useEtaWeights); 
  //weight by the Multiplicity
  Double_t dQX = 0.;
  Double_t dQY = 0.;
  if (vQ.GetMult()!=0) {
    dQX = vQ.X()/vQ.GetMult();
    dQY = vQ.Y()/vQ.GetMult();
  }
  vQ.Set(dQX,dQY);
  if(!fBookOnlyBasic){fHistQ->Fill(vQ.Mod());}
  if(!fBookOnlyBasic){fHistAngleQ->Fill(vQ.Phi()/2);}

  AliFlowVector* vQSub = new AliFlowVector[2];
  anEvent->Get2Qsub(vQSub, 2, weightsList, usePhiWeights, usePtWeights, useEtaWeights);
  AliFlowVector vQa = vQSub[0];
  AliFlowVector vQb = vQSub[1];
  if(!fBookOnlyBasic){fHistAngleQSub0->Fill(vQa.Phi()/2);}
  if(!fBookOnlyBasic){fHistAngleQSub1->Fill(vQb.Phi()/2);}

  Double_t dMultRP = 0.;
  Double_t dMultPOI = 0.;
  
  Int_t iNumberOfTracks = anEvent->NumberOfTracks();
  AliFlowTrackSimple* pTrack = NULL;     

  for (Int_t i=0;i<iNumberOfTracks;i++) {
    pTrack = anEvent->GetTrack(i);
    if (pTrack ) {
      dWeight = pTrack->Weight();
      dPt = pTrack->Pt();
      dPhi = pTrack->Phi();
      if (dPhi<0.) dPhi+=2*TMath::Pi();
      dEta = pTrack->Eta();

      //weights are only used for the RP selection
      if (pTrack->InRPSelection()){
	// determine Phi weight:
	if(phiWeights && nBinsPhi) {
	  dWPhi = phiWeights->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhi/TMath::TwoPi())));
	}
	// determine v'(pt) weight:
	if(ptWeights && dBinWidthPt) {
	  dWPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt)));
	}
	// determine v'(eta) weight:
	if(etaWeights && dBinWidthEta)  {
	  dWEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta)));
	}
	
	//the total weight is the product
	Double_t dW = dWeight*dWPhi*dWPt*dWEta; 

	//pt
	fHistPtRP->Fill(dPt,dW);
	//phi
	fHistPhiRP->Fill(dPhi,dW);
	//eta
	fHistEtaRP->Fill(dEta,dW);
	//eta vs phi
	if(!fBookOnlyBasic){fHistPhiEtaRP->Fill(dEta,dPhi,dW);}
	//weight vs phi
	if(!fBookOnlyBasic){fHistWeightvsPhi->Fill(dPhi,dW);}
	//count
	dMultRP += dW;
      }
      if (pTrack->InRPSelection() && pTrack->InSubevent(0)) {
	// determine Phi weight:
	if(phiWeightsSub0 && nBinsPhiSub0){
	  dWPhi = phiWeightsSub0->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhiSub0/TMath::TwoPi())));
	}
	// determine v'(pt) weight:
	if(ptWeights && dBinWidthPt) {
	  dWPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt)));
	}
	// determine v'(eta) weight:
	if(etaWeights && dBinWidthEta)  {
	  dWEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta)));
	}
	
	//the total weight is the product
	Double_t dW = dWeight*dWPhi*dWPt*dWEta;  

	//pt
	if(!fBookOnlyBasic){fHistPtSub0 ->Fill(dPt,dW);}
	//phi
	if(!fBookOnlyBasic){fHistPhiSub0 ->Fill(dPhi,dW);}
	//eta
	if(!fBookOnlyBasic){fHistEtaSub0 ->Fill(dEta,dW);}
      }
      if (pTrack->InRPSelection() && pTrack->InSubevent(1)) {
	// determine Phi weight:
	if(phiWeightsSub1 && nBinsPhiSub1){
	  dWPhi = phiWeightsSub1->GetBinContent(1+(Int_t)(TMath::Floor(dPhi*nBinsPhiSub1/TMath::TwoPi())));
	}
	// determine v'(pt) weight:
	if(ptWeights && dBinWidthPt) {
	  dWPt=ptWeights->GetBinContent(1+(Int_t)(TMath::Floor((dPt-dPtMin)/dBinWidthPt)));
	}
	// determine v'(eta) weight:
	if(etaWeights && dBinWidthEta)  {
	  dWEta=etaWeights->GetBinContent(1+(Int_t)(TMath::Floor((dEta-dEtaMin)/dBinWidthEta)));
	}
	
	//the total weight is the product
	Double_t dW = dWeight*dWPhi*dWPt*dWEta;  

	//pt
	if(!fBookOnlyBasic){fHistPtSub1 -> Fill(dPt,dW);}
	//phi
	if(!fBookOnlyBasic){fHistPhiSub1 -> Fill(dPhi,dW);}
	//eta
	if(!fBookOnlyBasic){fHistEtaSub1 -> Fill(dEta,dW);}
      }
      if (pTrack->InPOISelection()){

	Double_t dW = dWeight; //no pt, phi or eta weights

	//pt
	fHistPtPOI ->Fill(dPt,dW);
	//phi
	fHistPhiPOI ->Fill(dPhi,dW);
	//eta
	fHistEtaPOI ->Fill(dEta,dW);
	//eta vs phi
	if(!fBookOnlyBasic){fHistPhiEtaPOI ->Fill(dEta,dPhi,dW);}
	//mean pt
	fHistProMeanPtperBin ->Fill(dPt,dPt,dW);
	//count
	dMultPOI += dW;
      }
    } //track
  } //loop over tracks
  
  fHistMultRP->Fill(dMultRP);
  fHistMultPOI->Fill(dMultPOI);
  if(!fBookOnlyBasic){fHistMultPOIvsRP->Fill(dMultRP,dMultPOI);}
  
  //<reference multiplicity> versus # of RPs:
  fRefMultVsNoOfRPs->Fill(dMultRP+0.5,anEvent->GetReferenceMultiplicity(),1.);
  
  //reference multiplicity:
  fHistRefMult->Fill(anEvent->GetReferenceMultiplicity());

  return kTRUE; 
}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtRP
  Double_t dEntries = fHistPtRP->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInPtBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtPOI
  Double_t dEntries = fHistPtPOI->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinRP(Int_t aBin)
{
  //get entries in bin aBin from fHistPtRP
  Double_t dEntries = fHistEtaRP->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetEntriesInEtaBinPOI(Int_t aBin)
{
  //get entries in bin aBin from fHistPtPOI
  Double_t dEntries = fHistEtaPOI->GetBinContent(aBin);

  return dEntries;

}

//----------------------------------------------------------------------- 

Double_t AliFlowCommonHist::GetMeanPt(Int_t aBin)
{  
  //Get entry from bin aBin from fHistProMeanPtperBin
  Double_t dMeanPt = fHistProMeanPtperBin->GetBinContent(aBin);

  return dMeanPt;
  
}


//----------------------------------------------------------------------- 
 Double_t AliFlowCommonHist::Merge(TCollection *aList)
{
  //merge fuction
  //cout<<"entering merge function"<<endl;
  if (!aList) return 0;
  if (aList->IsEmpty()) return 0; //no merging is needed

  Int_t iCount = 0;
  TIter next(aList); // list is supposed to contain only objects of the same type as this
  AliFlowCommonHist *toMerge;
  // make a temporary list
  TList *pTemp = new TList();
  while ((toMerge=(AliFlowCommonHist*)next())) {
    pTemp->Add(toMerge->GetHistList()); 
    iCount++;
  }
  // Now call merge for fHistList providing temp list
  fHistList->Merge(pTemp);
  // Cleanup
  delete pTemp;
    
  //cout<<"Merged"<<endl;
  return (double)iCount;
    
}

void AliFlowCommonHist::Print(Option_t *option) const
{
  //   -*-*-*-*-*Print some global quantities for this histogram collection class *-*-*-*-*-*-*-*
  //             ===============================================
  //   printf( "TH1.Print Name  = %s, Entries= %d, Total sum= %g\n",GetName(),Int_t(fEntries),GetSumOfWeights());
  printf( "Class.Print Name = %s, Histogram list:\n",GetName());

  if (fHistList) {  
    fHistList->Print(option);
  }
  else
    {
      printf( "Empty histogram list \n");
    }
}

//----------------------------------------------------------------------- 
 void AliFlowCommonHist::Browse(TBrowser *b)
{

  if (!b) return;
  if (fHistList) b->Add(fHistList,"AliFlowCommonHistList");
}




