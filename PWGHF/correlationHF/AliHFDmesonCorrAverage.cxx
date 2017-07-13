/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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


/* $Id: $ */

/////////////////////////////////////////////////////////////
// class to average D meson -hadron correlations
//
// Author: A. Rossi, andrea.rossi@cern.ch
/////////////////////////////////////////////////////////////

#include "AliHFDmesonCorrAverage.h"
#include "AliHFDhadronCorrSystUnc.h"
#include <TMath.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>


//using std::cout;
//using std::endl;

ClassImp(AliHFDmesonCorrAverage)

AliHFDmesonCorrAverage::AliHFDmesonCorrAverage() : 
TNamed(),
  fSystDzero(0),
  fSystDstar(0),
  fSystDplus(0),
  fSystDaverage(0),
  fincludeDzero(kFALSE),
  fincludeDstar(kFALSE),                    
  fincludeDplus(kFALSE),                   
  fmethod(),
  fptminD(),
  fptmaxD(),
  fptminAsso(),
  fptmaxAsso(),
  fhDzero(0x0),
  fhDstar(0x0),
  fhDplus(0x0),
  fhDaverage(0x0),
  fgrTotSystAverage(0x0),
  fgrFDSystAverage(0x0),
  fgrNonFDSystAverage(0x0),
  fweightsDzeroStat(0x0),
  fweightsDstarStat(0x0),
  fweightsDplusStat(0x0),
  fweightsDzeroSystYield(0x0),
  fweightsDstarSystYield(0x0),
  fweightsDplusSystYield(0x0),
  fweightsDzeroSystBkg(0x0),
  fweightsDstarSystBkg(0x0),
  fweightsDplusSystBkg(0x0),
  fnbinsphi(),  
  fsys(-2),
  fyear(-2),
  fSystAlreadySet(kFALSE),
  fArithmeticAverage(kFALSE),
  fhUsedWeightsDzero(0x0),
  fhUsedWeightsDstar(0x0),
  fhUsedWeightsDplus(0x0)
{// standard constructor 
  
}


AliHFDmesonCorrAverage::AliHFDmesonCorrAverage(const char* name) : 
  TNamed(name,name),
  fSystDzero(0),
  fSystDstar(0),
  fSystDplus(0),
  fSystDaverage(0),
  fincludeDzero(kFALSE),
  fincludeDstar(kFALSE),                    
  fincludeDplus(kFALSE),                   
  fmethod(),
  fptminD(),
  fptmaxD(),
  fptminAsso(),
  fptmaxAsso(),
  fhDzero(0x0),
  fhDstar(0x0),
  fhDplus(0x0),
  fhDaverage(0x0),
  fgrTotSystAverage(0x0),
  fgrFDSystAverage(0x0),
  fgrNonFDSystAverage(0x0),
  fweightsDzeroStat(0x0),
  fweightsDstarStat(0x0),
  fweightsDplusStat(0x0),
  fweightsDzeroSystYield(0x0),
  fweightsDstarSystYield(0x0),
  fweightsDplusSystYield(0x0),
  fweightsDzeroSystBkg(0x0),
  fweightsDstarSystBkg(0x0),
  fweightsDplusSystBkg(0x0),
  fnbinsphi(),
  fsys(-2),
  fyear(-2),
  fSystAlreadySet(kFALSE),
  fArithmeticAverage(kFALSE),
  fhUsedWeightsDzero(0x0),
  fhUsedWeightsDstar(0x0),
  fhUsedWeightsDplus(0x0)
{// default constructor 

}

AliHFDmesonCorrAverage::~AliHFDmesonCorrAverage(){
  
  delete fSystDzero;
  delete fSystDstar;
  delete fSystDplus;
  delete fSystDaverage;
  delete fhDzero;
  delete fhDstar; 
  delete fhDplus;
  delete fhDaverage;
  delete fgrTotSystAverage;
  delete fgrFDSystAverage;
  delete fgrNonFDSystAverage;
  delete fweightsDzeroStat;
  delete fweightsDstarStat;
  delete fweightsDplusStat;
  delete fweightsDzeroSystYield;
  delete fweightsDstarSystYield;
  delete fweightsDplusSystYield;
  delete fweightsDzeroSystBkg;
  delete fweightsDstarSystBkg;
  delete fweightsDplusSystBkg;
  delete  fhUsedWeightsDzero;
  delete fhUsedWeightsDstar;
  delete fhUsedWeightsDplus;

}



Bool_t AliHFDmesonCorrAverage::InitSystematicUncertainty(Int_t system,Int_t year){
  if(!fSystAlreadySet){
    if(system!=-1){
      if(system==0){ //pp
	if(fincludeDzero){
	  if(year==2010){
	    fSystDzero=new AliHFDhadronCorrSystUnc("fSystDzero");
	    fSystDzero->InitStandardUncertaintiesPP2010(0,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDzero->BuildTotalUncHisto();
	    fSystDzero->BuildTotalNonFDUncHisto();
	    fSystDzero->BuildGraphsUnc(fhDzero);
	    fSystDzero->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
	
	if(fincludeDstar){
	  if(year==2010){
	    fSystDstar=new AliHFDhadronCorrSystUnc("fSystDstar");
	    fSystDstar->InitStandardUncertaintiesPP2010(1,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDstar->BuildTotalUncHisto();
	    fSystDstar->BuildTotalNonFDUncHisto();
	    fSystDstar->BuildGraphsUnc(fhDstar);
	    fSystDstar->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
	
	if(fincludeDplus){
	  if(year==2010){
	    fSystDplus=new AliHFDhadronCorrSystUnc("fSystDplus");
	    fSystDplus->InitStandardUncertaintiesPP2010(2,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDplus->BuildTotalUncHisto();
	    fSystDplus->BuildTotalNonFDUncHisto();
	    fSystDplus->BuildGraphsUnc(fhDplus);
	    fSystDplus->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
      }
      else if(system==1){ //p-Pb
	if(fincludeDzero){
	  if(year==2013){
	    fSystDzero=new AliHFDhadronCorrSystUnc("fSystDzero");
	    fSystDzero->InitStandardUncertaintiesPPb2013(0,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDzero->BuildTotalUncHisto();
	    fSystDzero->BuildTotalNonFDUncHisto();
	    fSystDzero->BuildGraphsUnc(fhDzero);
	    fSystDzero->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
	if(fincludeDstar){
	  if(year==2013){
	    fSystDstar=new AliHFDhadronCorrSystUnc("fSystDstar");
	    fSystDstar->InitStandardUncertaintiesPPb2013(1,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDstar->BuildTotalUncHisto();
	    fSystDstar->BuildTotalNonFDUncHisto();
	    fSystDstar->BuildGraphsUnc(fhDstar);
	    fSystDstar->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
	
	if(fincludeDplus){
	  if(year==2013){
	    fSystDplus=new AliHFDhadronCorrSystUnc("fSystDplus");
	    fSystDplus->InitStandardUncertaintiesPPb2013(2,(fptminD+fptmaxD)*0.5,fptminAsso);  
	    fSystDplus->BuildTotalUncHisto();
	    fSystDplus->BuildTotalNonFDUncHisto();
	    fSystDplus->BuildGraphsUnc(fhDplus);
	    fSystDplus->BuildGraphsRelUnc();
	  }
	  else {
	    Printf("No values for syst unc foreseen for this dataset");
	  }
	}
      }
      else if(system==2){ //p-Pb 2016
  if(fincludeDzero){
    if(year==2016){
      fSystDzero=new AliHFDhadronCorrSystUnc("fSystDzero");
      fSystDzero->InitStandardUncertaintiesPPb2016(0,(fptminD+fptmaxD)*0.5,fptminAsso,fptmaxAsso);  
      fSystDzero->BuildTotalUncHisto();
      fSystDzero->BuildTotalNonFDUncHisto();
      fSystDzero->BuildGraphsUnc(fhDzero);
      fSystDzero->BuildGraphsRelUnc();
    }
    else {
      Printf("No values for syst unc foreseen for this dataset");
    }
  }
  if(fincludeDstar){
    if(year==2016){
      fSystDstar=new AliHFDhadronCorrSystUnc("fSystDstar");
      fSystDstar->InitStandardUncertaintiesPPb2016(1,(fptminD+fptmaxD)*0.5,fptminAsso,fptmaxAsso);  
      fSystDstar->BuildTotalUncHisto();
      fSystDstar->BuildTotalNonFDUncHisto();
      fSystDstar->BuildGraphsUnc(fhDstar);
      fSystDstar->BuildGraphsRelUnc();
    }
    else {
      Printf("No values for syst unc foreseen for this dataset");
    }
  }
  
  if(fincludeDplus){
    if(year==2016){
      fSystDplus=new AliHFDhadronCorrSystUnc("fSystDplus");
      fSystDplus->InitStandardUncertaintiesPPb2016(2,(fptminD+fptmaxD)*0.5,fptminAsso,fptmaxAsso);  
      fSystDplus->BuildTotalUncHisto();
      fSystDplus->BuildTotalNonFDUncHisto();
      fSystDplus->BuildGraphsUnc(fhDplus);
      fSystDplus->BuildGraphsRelUnc();
    }
    else {
      Printf("No values for syst unc foreseen for this dataset");
    }
  }
      }
      else {
	Printf("Cannot initiate syst uncertainties: wrong system selected");
	return kFALSE;
      }
    }
  }
  
  fSystDaverage=new AliHFDhadronCorrSystUnc("fSystDaverage");
  if(fincludeDzero)fSystDaverage->SetHistoTemplate(fSystDzero->GetHistoTemplate(),"fhDeltaPhiTemplateDaverage");
  else if(fincludeDstar)fSystDaverage->SetHistoTemplate(fSystDstar->GetHistoTemplate(),"fhDeltaPhiTemplateDaverage");
  else if(fincludeDplus)fSystDaverage->SetHistoTemplate(fSystDplus->GetHistoTemplate(),"fhDeltaPhiTemplateDaverage");
  
  fSystDaverage->InitEmptyHistosFromTemplate();
  
  return kTRUE;
}


void AliHFDmesonCorrAverage::InitAverageHisto(TH1D *h){
  fhDaverage=(TH1D*)h->Clone("fhDaverage");  
  fhDaverage->SetXTitle("#Delta#varphi = #varphi_{assoc} - #varphi_{D}");
  fhDaverage->SetYTitle("#frac{1}{N_{D}}#frac{dN^{h}}{d#Delta#varphi} (rad^{-1})");
  fhDaverage->SetTitle("");
  fhDaverage->SetMinimum(0);

  fhDaverage->Reset(0);
  fhDaverage->SetLineWidth(2);
  fhDaverage->SetLineColor(kBlack);
  fhDaverage->SetMarkerStyle(20);
  fhDaverage->SetMarkerSize(1.2);
  fhDaverage->SetMarkerColor(kBlack);

  // The following histos are created here to use the final binning
  fhUsedWeightsDzero=(TH1D*)h->Clone("fhUsedWeightsDzero");
  fhUsedWeightsDzero->SetTitle("Dzero final weights used");
  fhUsedWeightsDzero->SetXTitle("#Delta#varphi = #varphi_{assoc} - #varphi_{D}");
  fhUsedWeightsDzero->SetYTitle("weight");

  fhUsedWeightsDplus=(TH1D*)h->Clone("fhUsedWeightsDplus");
  fhUsedWeightsDplus->SetTitle("Dplus final weights used");
  fhUsedWeightsDplus->SetXTitle("#Delta#varphi = #varphi_{assoc} - #varphi_{D}");
  fhUsedWeightsDplus->SetYTitle("weight");

  fhUsedWeightsDstar=(TH1D*)h->Clone("fhUsedWeightsDstar");
  fhUsedWeightsDstar->SetTitle("Dstar final weights used");
  fhUsedWeightsDstar->SetXTitle("#Delta#varphi = #varphi_{assoc} - #varphi_{D}");
  fhUsedWeightsDstar->SetYTitle("weight");



}

void AliHFDmesonCorrAverage::CalculateAverage(){
  if(fmethod!=10&&fmethod!=20){
    Printf("This option for the calculation of the average has not been implemented yet");
    return;
  }
  // check that histos were set
  if(fincludeDzero&&(!fhDzero)){
    Printf("Dzero histo not set");
    return;
  }
  if(fincludeDstar&&(!fhDstar)){
    Printf("Dstar histo not set");
    return;
  }
  if(fincludeDplus&&(!fhDplus)){
    Printf("Dplus histo not set");
    return;
  }
  

  // check that they have the same binning
  if(fincludeDzero&&fincludeDstar){
    if(fhDzero->GetNbinsX()!=fhDstar->GetNbinsX()){
      Printf("Dzero and Dstar histos have different binning");    
      return ;
    }
    fnbinsphi=fhDzero->GetNbinsX();
    for(Int_t j=1;j<=fnbinsphi;j++){
      if(TMath::Abs(fhDzero->GetBinLowEdge(j)-fhDstar->GetBinLowEdge(j))>0.001){
	Printf("Dzero and Dstar histos have different binning");    
	return;
      }
    }
    InitAverageHisto(fhDzero);
  }
  if(fincludeDzero&&fincludeDplus){
    if(fhDzero->GetNbinsX()!=fhDplus->GetNbinsX()){
    Printf("Dzero and Dplus histos have different binning");    
    return ;
    }
    fnbinsphi=fhDzero->GetNbinsX();
    for(Int_t j=1;j<=fnbinsphi;j++){
      if(TMath::Abs(fhDzero->GetBinLowEdge(j)-fhDplus->GetBinLowEdge(j))>0.001){
	Printf("Dzero and Dplus histos have different binning");    
	return;
      }
    }
    if(!fhDaverage)InitAverageHisto(fhDzero);
  }
  if(fincludeDstar&&fincludeDplus){
    if(fhDstar->GetNbinsX()!=fhDplus->GetNbinsX()){
      Printf("Dstar and Dplus histos have different binning");    
      return ;
    }
    fnbinsphi=fhDstar->GetNbinsX();
    for(Int_t j=1;j<=fnbinsphi;j++){
      if(TMath::Abs(fhDstar->GetBinLowEdge(j)-fhDplus->GetBinLowEdge(j))>0.001){
	Printf("Dplus and Dstar histos have different binning");    
	return;
      }
    }
    if(!fhDaverage)InitAverageHisto(fhDstar);
  }
  
  if(!(fincludeDstar||fincludeDplus||fincludeDzero)){
    Printf("No mesons choosen for average");
    return ;
  }
  Int_t nmeson=0;
  if(fincludeDzero)nmeson++;
  if(fincludeDstar)nmeson++;
  if(fincludeDplus)nmeson++;
  
  printf("Settin syst \n");

  if(!InitSystematicUncertainty(fsys,fyear)){
    Printf("Cannot init syst uncertainties \n");
    return;
  }
  
  SetWeights();
  printf("Weights set \n");
  
  for(Int_t j=1;j<=fnbinsphi;j++){
    Double_t value=0.,errStatValue=0.,errSystYieldValue=0.,errSystBkgValue=0.,weight=0.;
    Double_t systMCclosureMin=0.,systMCcorrectionsMin=0.,systFDmin=0.,systMCDefficiencyMin=0.,systSecContaminationMin=0.;
    Double_t systMCclosureMax=0.,systMCcorrectionsMax=0.,systFDmax=0.,systMCDefficiencyMax=0.,systSecContaminationMax=0.;
    Double_t sumweights=0;
    Double_t errsyst;

    if(fincludeDzero){
      // stat error + yield unc + bkg subtr
      if(fArithmeticAverage){
	weight=1./nmeson;
      }
      else {
	weight=1./(1./fweightsDzeroStat[j-1]+1./fweightsDzeroSystYield[j-1]+1./fweightsDzeroSystBkg[j-1]);// need to do this way since we stored separately the stat and syst weight (=1/variance)
      }
      fhUsedWeightsDzero->SetBinContent(j,weight);

      value+=fhDzero->GetBinContent(j)*weight;
      errStatValue+=1./fweightsDzeroStat[j-1]*weight*weight;
      errSystYieldValue+=1./fweightsDzeroSystYield[j-1]*weight*weight;
      errSystBkgValue+=1./fweightsDzeroSystBkg[j-1]*weight*weight;
      sumweights+=weight;
      
      printf("Dzero the value is: %f, weight: %f \n",value, weight);
      // MCcorrections  : associated tracks, correlated
      errsyst=TMath::Abs(fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCcorrectionsMin()->GetBinContent(j));
      systMCcorrectionsMin+=errsyst*weight;
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCcorrectionsMax()->GetBinContent(j);
      systMCcorrectionsMax+=errsyst*weight;

      printf(" Dzero trackeff the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoMCcorrectionsMin()->GetBinContent(j), systMCcorrectionsMin);
      printf(" Dzero trackeff the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoMCcorrectionsMax()->GetBinContent(j), systMCcorrectionsMax);

      // MCDefficiency  : uncorrelated
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCDefficiencyMin()->GetBinContent(j);
      systMCDefficiencyMin+=errsyst*errsyst*weight*weight;
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCDefficiencyMax()->GetBinContent(j);
      systMCDefficiencyMax+=errsyst*errsyst*weight*weight;

      // SecContamination  : correlated
      errsyst=TMath::Abs(fhDzero->GetBinContent(j)*fSystDzero->GetHistoSecContaminationMin()->GetBinContent(j));
      systSecContaminationMin+=errsyst*weight;
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoSecContaminationMax()->GetBinContent(j);
      systSecContaminationMax+=errsyst*weight;

      // MC closure: fully correlated
      errsyst=TMath::Abs(fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCclosureTestMin()->GetBinContent(j));
      systMCclosureMin+=errsyst*weight;
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoMCclosureTestMax()->GetBinContent(j);
      systMCclosureMax+=errsyst*weight;

      printf(" Dzero Mcclosure the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoMCclosureTestMin()->GetBinContent(j), systMCclosureMin);
      printf(" Dzero Mcclosure the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoMCclosureTestMax()->GetBinContent(j), systMCclosureMax);

      // FD (assume full correlations)
      errsyst=TMath::Abs(fhDzero->GetBinContent(j)*fSystDzero->GetHistoFDmin()->GetBinContent(j));
      systFDmin+=errsyst*weight;
      errsyst=fhDzero->GetBinContent(j)*fSystDzero->GetHistoFDmax()->GetBinContent(j);
      systFDmax+=errsyst*weight;
      printf(" Dzero feeddown the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoFDmin()->GetBinContent(j), systFDmin);
      printf(" Dzero feeddown the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDzero->GetHistoFDmax()->GetBinContent(j), systFDmax);

      printf("Dzero the value is: %f, weight: %f \n",value, weight);
    } 
    if(fincludeDstar){
      if(fArithmeticAverage){
	weight=1./nmeson;
      }
      else{
	// stat error + yield unc + bkg subtr
	weight=1./(1./fweightsDstarStat[j-1]+1./fweightsDstarSystYield[j-1]+1./fweightsDstarSystBkg[j-1]);// need to do this way since we stored separately the stat and syst weight (=1/variance)
      }
      fhUsedWeightsDstar->SetBinContent(j,weight);
      value+=fhDstar->GetBinContent(j)*weight;
      errStatValue+=1./fweightsDstarStat[j-1]*weight*weight;
      errSystYieldValue+=1./fweightsDstarSystYield[j-1]*weight*weight;
      errSystBkgValue+=1./fweightsDstarSystBkg[j-1]*weight*weight;
      sumweights+=weight;

      printf("Dstar the value is: %f, weight: %f \n",value, weight);

      // MCcorrections: associated tracks, correlated
      errsyst=TMath::Abs(fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCcorrectionsMin()->GetBinContent(j));
      systMCcorrectionsMin+=errsyst*weight;
      errsyst=TMath::Abs(fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCcorrectionsMax()->GetBinContent(j));
      systMCcorrectionsMax+=errsyst*weight;

      printf(" Dstar trackeff the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoMCcorrectionsMin()->GetBinContent(j), systMCcorrectionsMin);
      printf(" Dstar trackeff the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoMCcorrectionsMax()->GetBinContent(j), systMCcorrectionsMax);

      // MCDefficiency  : uncorrelated
      errsyst=fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCDefficiencyMin()->GetBinContent(j);
      systMCDefficiencyMin+=errsyst*errsyst*weight*weight;
      errsyst=fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCDefficiencyMax()->GetBinContent(j);
      systMCDefficiencyMax+=errsyst*errsyst*weight*weight;

      // SecContamination  : correlated
      errsyst=TMath::Abs(fhDstar->GetBinContent(j)*fSystDstar->GetHistoSecContaminationMin()->GetBinContent(j));
      systSecContaminationMin+=errsyst*weight;
      errsyst=fhDstar->GetBinContent(j)*fSystDstar->GetHistoSecContaminationMax()->GetBinContent(j);
      systSecContaminationMax+=errsyst*weight;

      // MC closure
      errsyst=TMath::Abs(fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCclosureTestMin()->GetBinContent(j));
      systMCclosureMin+=errsyst*weight;
      errsyst=fhDstar->GetBinContent(j)*fSystDstar->GetHistoMCclosureTestMax()->GetBinContent(j);
      systMCclosureMax+=errsyst*weight;

      printf(" Dstar Mcclosure the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoMCclosureTestMin()->GetBinContent(j), systMCclosureMin);
      printf(" Dstar Mcclosure the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoMCclosureTestMax()->GetBinContent(j), systMCclosureMax);


      // FD (assume full correlations)
      errsyst=TMath::Abs(fhDstar->GetBinContent(j)*fSystDstar->GetHistoFDmin()->GetBinContent(j));
      systFDmin+=errsyst*weight;
      errsyst=fhDstar->GetBinContent(j)*fSystDstar->GetHistoFDmax()->GetBinContent(j);
      systFDmax+=errsyst*weight;

      printf(" Dstar feeddown the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoFDmin()->GetBinContent(j), systFDmin);
      printf(" Dstar feeddown the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDstar->GetHistoFDmax()->GetBinContent(j), systFDmax);

      printf("Dstar the value is: %f, weight: %f \n",value, weight);
    }
    if(fincludeDplus){
      if(fArithmeticAverage){
	weight=1./nmeson;
      }
      else{
	// stat error + yield unc + bkg subtr
	weight=1./(1./fweightsDplusStat[j-1]+1./fweightsDplusSystYield[j-1]+1./fweightsDplusSystBkg[j-1]);// need to do this way since we stored separately the stat and syst weight (=1/variance)
      }
      fhUsedWeightsDplus->SetBinContent(j,weight);
      value+=fhDplus->GetBinContent(j)*weight;
      errStatValue+=1./fweightsDplusStat[j-1]*weight*weight;
      errSystYieldValue+=1./fweightsDplusSystYield[j-1]*weight*weight;      
      errSystBkgValue+=1./fweightsDplusSystBkg[j-1]*weight*weight;
      sumweights+=weight;

      printf("Dplus the value is: %f, weight: %f \n",value, weight);

      // MCcorrections  : associated tracsk, correlated
      errsyst=TMath::Abs(fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCcorrectionsMin()->GetBinContent(j));
      systMCcorrectionsMin+=errsyst*weight;
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCcorrectionsMax()->GetBinContent(j);
      systMCcorrectionsMax+=errsyst*weight;

      printf(" Dplus trackeff the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCcorrectionsMin()->GetBinContent(j), systMCcorrectionsMin);
      printf(" Dplus trackeff the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCcorrectionsMax()->GetBinContent(j), systMCcorrectionsMax);

      // MCDefficiency  : uncorrelated
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCDefficiencyMin()->GetBinContent(j);
      systMCDefficiencyMin+=errsyst*errsyst*weight*weight;
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCDefficiencyMax()->GetBinContent(j);
      systMCDefficiencyMax+=errsyst*errsyst*weight*weight;

      printf(" Dplus cutvariat the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCDefficiencyMin()->GetBinContent(j), systMCDefficiencyMin);
      printf(" Dplus cutvariat the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCDefficiencyMax()->GetBinContent(j), systMCDefficiencyMax);

      // SecContamination  : correlated
      errsyst=TMath::Abs(fhDplus->GetBinContent(j)*fSystDplus->GetHistoSecContaminationMin()->GetBinContent(j));
      systSecContaminationMin+=errsyst*weight;
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoSecContaminationMax()->GetBinContent(j);
      systSecContaminationMax+=errsyst*weight;

      // MC closure
      errsyst=TMath::Abs(fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCclosureTestMin()->GetBinContent(j));
      systMCclosureMin+=errsyst*weight;
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoMCclosureTestMax()->GetBinContent(j);
      systMCclosureMax+=errsyst*weight;

      printf(" Dplus Mcclosure the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCclosureTestMin()->GetBinContent(j), systMCclosureMin);
      printf(" Dplus Mcclosure the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoMCclosureTestMax()->GetBinContent(j), systMCclosureMax);

      // FD (assume full correlations)
      errsyst=TMath::Abs(fhDplus->GetBinContent(j)*fSystDplus->GetHistoFDmin()->GetBinContent(j));
      systFDmin+=errsyst*weight;
      errsyst=fhDplus->GetBinContent(j)*fSystDplus->GetHistoFDmax()->GetBinContent(j);
      systFDmax+=errsyst*weight;

      printf(" Dplus feeddown the min syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoFDmin()->GetBinContent(j), systMCclosureMin);
      printf(" Dplus feeddown the max syst value is: %f (rel %f), sum: %f \n",errsyst,fSystDplus->GetHistoFDmax()->GetBinContent(j), systMCclosureMax);

      printf("Dplus the value is: %f, weight: %f \n",value, weight);

    }
    
    value/=sumweights;
    errStatValue/=sumweights*sumweights;// was: ((Double_t)nmeson*sumweights);
    errSystYieldValue/=sumweights*sumweights;// was: ((Double_t)nmeson*sumweights);
    errSystBkgValue/=sumweights*sumweights; // was: ((Double_t)nmeson*sumweights);
    fhDaverage->SetBinContent(j,value);
    fhDaverage->SetBinError(j,TMath::Sqrt(errStatValue));

    // Settting unc on yield and back sub
    fSystDaverage->GetHistoYieldUnc()->SetBinContent(j,TMath::Sqrt(errSystYieldValue)/TMath::Abs(value));
    fSystDaverage->GetHistoBackSubUncMin()->SetBinContent(j,-TMath::Sqrt(errSystBkgValue)/TMath::Abs(value));
    fSystDaverage->GetHistoBackSubUncMax()->SetBinContent(j,TMath::Sqrt(errSystBkgValue)/TMath::Abs(value));



    // MCcorrections, associated tracks  
    if(systMCcorrectionsMin>=0.){
      systMCcorrectionsMin=systMCcorrectionsMin/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }
    
    if(systMCcorrectionsMax>=0.){
      systMCcorrectionsMax=systMCcorrectionsMax/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }


    // MCDefficiency  
    if(systMCDefficiencyMin>=0.){
      systMCDefficiencyMin=TMath::Sqrt(systMCDefficiencyMin)/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }
    
    if(systMCDefficiencyMax>=0.){
      systMCDefficiencyMax=TMath::Sqrt(systMCDefficiencyMax)/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }

    // SecContamination  
    if(systSecContaminationMin>=0.){
      systSecContaminationMin=systSecContaminationMin/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }
    
    if(systSecContaminationMax>=0.){
      systSecContaminationMax=systSecContaminationMax/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }

    // MCclosureTest  
    if(systMCclosureMin>=0.){
      systMCclosureMin=systMCclosureMin/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }

    if(systMCclosureMax>=0.){
      systMCclosureMax=systMCclosureMax/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }



    // Feed down  
    if(systFDmin>=0.){
      systFDmin=systFDmin/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }
    
    if(systFDmax>=0.){
      systFDmax=systFDmax/sumweights;
    }
    else {
      Printf("Negative sum in error calculation \n");
      return;
    }

    value=TMath::Abs(value);// protection to avoid flipping of min syst histo sign

    fSystDaverage->GetHistoMCcorrectionsMin()->SetBinContent(j,-systMCcorrectionsMin/value);
    fSystDaverage->GetHistoMCcorrectionsMax()->SetBinContent(j,systMCcorrectionsMax/value);

    fSystDaverage->GetHistoMCDefficiencyMin()->SetBinContent(j,-systMCDefficiencyMin/value);
    fSystDaverage->GetHistoMCDefficiencyMax()->SetBinContent(j,systMCDefficiencyMax/value);

    fSystDaverage->GetHistoSecContaminationMin()->SetBinContent(j,-systSecContaminationMin/value);
    fSystDaverage->GetHistoSecContaminationMax()->SetBinContent(j,systSecContaminationMax/value);

    fSystDaverage->GetHistoMCclosureTestMin()->SetBinContent(j,-systMCclosureMin/value);
    fSystDaverage->GetHistoMCclosureTestMax()->SetBinContent(j,systMCclosureMax/value);
    
    fSystDaverage->GetHistoFDmin()->SetBinContent(j,-systFDmin/value);
    fSystDaverage->GetHistoFDmax()->SetBinContent(j,systFDmax/value);

  }
  fSystDaverage->BuildTotalUncHisto();
  fSystDaverage->BuildTotalNonFDUncHisto();
  fSystDaverage->BuildGraphsUnc(fhDaverage);
  fSystDaverage->BuildGraphsRelUnc();
  
  return ;
  
}


TH1D* AliHFDmesonCorrAverage::ReflectHisto(TH1D *h){
  
  TH1D *h2=new TH1D(Form("%sReflected",h->GetName()),Form("%sReflected",h->GetName()),h->GetNbinsX()/2.,0.,TMath::Pi());
  for(Int_t j=1;j<=h->GetNbinsX();j++){
    Double_t x=h->GetBinCenter(j);
    Double_t y0=h->GetBinContent(j);
    Double_t ey0=h->GetBinError(j);
    Int_t j2;
    if(x>0&&x<TMath::Pi())j2=h2->FindBin(x);
    else if(x<0)j2=h2->FindBin(-1.*x);
    else if(x>TMath::Pi())j2=h2->FindBin(2.*TMath::Pi()-x);
    else {
      printf("Point %d excluded \n",j);
      continue;
    }
    Double_t y=h2->GetBinContent(j2);
    Double_t ey=h2->GetBinError(j2);
    h2->SetBinContent(j2,y+y0);
    h2->SetBinError(j2,TMath::Sqrt(ey0*ey0+ey*ey));
  }
  
  return h2;
}

void AliHFDmesonCorrAverage::SetWeights(){

  
  if(fincludeDzero){    
    TH1D *hYieldUnc=fSystDzero->GetHistoYieldUnc();
    TH1D *hBkgUncMin=fSystDzero->GetHistoBackSubUncMin();
    TH1D *hBkgUncMax=fSystDzero->GetHistoBackSubUncMax();
    if(fweightsDzeroStat)delete fweightsDzeroStat;
    if(fweightsDzeroSystYield)delete fweightsDzeroSystYield;
    if(fweightsDzeroSystBkg)delete fweightsDzeroSystBkg;
    fweightsDzeroStat=new Double_t[fnbinsphi];    
    fweightsDzeroSystYield=new Double_t[fnbinsphi];    
    fweightsDzeroSystBkg=new Double_t[fnbinsphi];    
    //    fhGlobalWeightDzero=new TH1F("fhGlobalWeightDzero","fhGlobalWeightDzero",fnbinsphi
    for(Int_t j=0;j<fnbinsphi;j++){
      if(fmethod==10){
	fweightsDzeroStat[j]=1./(fhDzero->GetBinError(j+1)*fhDzero->GetBinError(j+1));
	fweightsDzeroSystYield[j]=1./(hYieldUnc->GetBinContent(j+1)*fhDzero->GetBinContent(j+1)*hYieldUnc->GetBinContent(j+1)*fhDzero->GetBinContent(j+1)); 
	//for asymmetric uncertainties...
	if(TMath::Abs(hBkgUncMin->GetBinContent(j+1)) > TMath::Abs(hBkgUncMax->GetBinContent(j+1))) fweightsDzeroSystBkg[j]=1./(hBkgUncMin->GetBinContent(j+1)*fhDzero->GetBinContent(j+1)*hBkgUncMin->GetBinContent(j+1)*fhDzero->GetBinContent(j+1));
	else fweightsDzeroSystBkg[j]=1./(hBkgUncMax->GetBinContent(j+1)*fhDzero->GetBinContent(j+1)*hBkgUncMax->GetBinContent(j+1)*fhDzero->GetBinContent(j+1));
      }
      else{
	Printf("This option for the calculation of the average has not been implemented yet");
	return;
      }
    }
  }

  if(fincludeDstar){    
    TH1D *hYieldUnc=fSystDstar->GetHistoYieldUnc();
    TH1D *hBkgUncMin=fSystDstar->GetHistoBackSubUncMin();
    TH1D *hBkgUncMax=fSystDstar->GetHistoBackSubUncMax();
    if(fweightsDstarStat)delete fweightsDstarStat;
    if(fweightsDstarSystYield)delete fweightsDstarSystYield;
    if(fweightsDstarSystBkg)delete fweightsDstarSystBkg;
    fweightsDstarStat=new Double_t[fnbinsphi];
    fweightsDstarSystYield=new Double_t[fnbinsphi];    
    fweightsDstarSystBkg=new Double_t[fnbinsphi];    
    for(Int_t j=0;j<fnbinsphi;j++){
      if(fmethod==10){
	fweightsDstarStat[j]=1./(fhDstar->GetBinError(j+1)*fhDstar->GetBinError(j+1));
	fweightsDstarSystYield[j]=1./(hYieldUnc->GetBinContent(j+1)*fhDstar->GetBinContent(j+1)*hYieldUnc->GetBinContent(j+1)*fhDstar->GetBinContent(j+1)); 
	//for asymmetric uncertainties...
	if(TMath::Abs(hBkgUncMin->GetBinContent(j+1)) > TMath::Abs(hBkgUncMax->GetBinContent(j+1))) fweightsDstarSystBkg[j]=1./(hBkgUncMin->GetBinContent(j+1)*fhDstar->GetBinContent(j+1)*hBkgUncMin->GetBinContent(j+1)*fhDstar->GetBinContent(j+1));
	else fweightsDstarSystBkg[j]=1./(hBkgUncMax->GetBinContent(j+1)*fhDstar->GetBinContent(j+1)*hBkgUncMax->GetBinContent(j+1)*fhDstar->GetBinContent(j+1));
      }
      else{
	Printf("This option for the calculation of the average has not been implemented yet");
	return;
      }
    }
  }
  
  if(fincludeDplus){    
    TH1D *hYieldUnc=fSystDplus->GetHistoYieldUnc();
    TH1D *hBkgUncMin=fSystDplus->GetHistoBackSubUncMin();
    TH1D *hBkgUncMax=fSystDplus->GetHistoBackSubUncMax();
    if(fweightsDplusStat)delete fweightsDplusStat;
    if(fweightsDplusSystYield)delete fweightsDplusSystYield;
    if(fweightsDplusSystBkg)delete fweightsDplusSystBkg;
    fweightsDplusStat=new Double_t[fnbinsphi];    
    fweightsDplusSystYield=new Double_t[fnbinsphi];    
    fweightsDplusSystBkg=new Double_t[fnbinsphi];    
    for(Int_t j=0;j<fnbinsphi;j++){
      if(fmethod==10){
	fweightsDplusStat[j]=1./(fhDplus->GetBinError(j+1)*fhDplus->GetBinError(j+1));
	fweightsDplusSystYield[j]=1./(hYieldUnc->GetBinContent(j+1)*fhDplus->GetBinContent(j+1)*hYieldUnc->GetBinContent(j+1)*fhDplus->GetBinContent(j+1)); 
	//for asymmetric uncertainties...
	if(TMath::Abs(hBkgUncMin->GetBinContent(j+1)) > TMath::Abs(hBkgUncMax->GetBinContent(j+1))) fweightsDplusSystBkg[j]=1./(hBkgUncMin->GetBinContent(j+1)*fhDplus->GetBinContent(j+1)*hBkgUncMin->GetBinContent(j+1)*fhDplus->GetBinContent(j+1));
	else fweightsDplusSystBkg[j]=1./(hBkgUncMax->GetBinContent(j+1)*fhDplus->GetBinContent(j+1)*hBkgUncMax->GetBinContent(j+1)*fhDplus->GetBinContent(j+1));
      }
      else{
	Printf("This option for the calculation of the average has not been implemented yet");
	return;
      }
    }
  }
  
  
  
  
}


