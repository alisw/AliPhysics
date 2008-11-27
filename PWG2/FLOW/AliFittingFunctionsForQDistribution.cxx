/*************************************************************************
* Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/******************************** 
 * integrated flow estimate by  *
 *   fitting q-distribution     * 
 *                              *
 * author: Ante Bilandzic       * 
 *          (anteb@nikhef.nl)   *
 *                              *
 * based on the macro written   *
 *     by Sergei Voloshin       *                        
 *******************************/ 

#define AliFittingFunctionsForQDistribution_cxx

#include "Riostream.h"
#include "AliFlowCommonConstants.h"
#include "TChain.h"
#include "TFile.h"
#include "TList.h"
#include "TParticle.h"

#include "TProfile.h"
#include "TF1.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"
#include "TF1.h"
#include "TMath.h"

#include "AliFlowEventSimple.h"
#include "AliFlowTrackSimple.h"
#include "AliFittingQDistribution.h"
#include "AliFlowCommonConstants.h"
#include "AliFittingFunctionsForQDistribution.h"

ClassImp(AliFittingFunctionsForQDistribution)

//================================================================================================================_

AliFittingFunctionsForQDistribution::AliFittingFunctionsForQDistribution():  
 fAvMultFQD(NULL),
 fQDistributionFQD(NULL), 
 fIntFlowResFQD(NULL)
{
 //default constructor 
}

AliFittingFunctionsForQDistribution::~AliFittingFunctionsForQDistribution()
{
 //destructor
}

AliFittingFunctionsForQDistribution::AliFittingFunctionsForQDistribution(TProfile *AvMult, TH1D *QDistribution, TH1D *intFlowRes):
 fAvMultFQD(AvMult),
 fQDistributionFQD(QDistribution),
 fIntFlowResFQD(intFlowRes)
{
 //custom constructor 
}
  
//================================================================================================================

void AliFittingFunctionsForQDistribution::Calculate()
{
 //fitting q-distribution
 
 Int_t n=2;//harmonic (to be improved)
 
 Double_t AvM = fAvMultFQD->GetBinContent(1);
 
 Int_t nEvts = (Int_t)(fAvMultFQD->GetBinEntries(1));

 Double_t qmin=(fQDistributionFQD->GetXaxis())->GetXmin();  
 Double_t qmax=(fQDistributionFQD->GetXaxis())->GetXmax(); 
 Double_t bin=fQDistributionFQD->GetBinWidth(4);//assuming that all bins have the same width 
 Double_t ent=fQDistributionFQD->GetEntries();
 Double_t norm=bin*ent;//assuming that all bins have the same width

 TF1 *fittingFun = new TF1("fittingFun","[2]*(x/[1])*exp(-(x*x+[0]*[0])/(2.*[1]))*TMath::BesselI0(x*[0]/[1])",qmin,qmax); 
 
 fittingFun->SetParNames("V","sigma","norm");
 fittingFun->SetParameters(0.1*pow(AvM,0.5),0.5,norm);
 
 fittingFun->SetParLimits(0,0.,10.);//to be improved (limits)
 fittingFun->SetParLimits(1,0.0,5.5);//to be improved (limits)
 fittingFun->FixParameter(2,norm); 

 fQDistributionFQD->Fit("fittingFun","NQ","",qmin,qmax);
 
 Double_t v=0.,errorv=0.,sigma2=0.,errorsigma2=0.;
 if(AvM)
 { 
  v = fittingFun->GetParameter(0)/pow(AvM,0.5);
  errorv = fittingFun->GetParError(0)/pow(AvM,0.5);
 }
 
 sigma2 = fittingFun->GetParameter(1);
 errorsigma2 = fittingFun->GetParError(1);
 
 cout<<" "<<endl;
 cout<<"************************************"<<endl;
 cout<<"************************************"<<endl;
 cout<<"    integrated flow by fitting"<<endl;
 cout<<"         q-distribution:      "<<endl;
 cout<<""<<endl;
 cout<<" v_"<<n<<"{FQD} = "<<v<<" +/- "<<errorv<<endl;
 cout<<" sigma^2  = "<<sigma2<<" +/- "<<errorsigma2<<endl; 
 //cout<<"vm       = "<<v*pow(AvM,0.5)<<endl;
 cout<<" "<<endl;
 cout<<"    nEvts = "<<nEvts<<", AvM = "<<AvM<<endl; 
 cout<<"************************************"<<endl;
 cout<<"************************************"<<endl;
 fIntFlowResFQD->SetBinContent(1,v);
 fIntFlowResFQD->SetBinError(1,errorv);
 cout<<" "<<endl;

 //fQDistributionFQD->Draw("");
 //fittingFun->Draw("SAME");

}//end of Calculate()

  













