/*************************************************************************
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

/* $Id$ */

//
// Class for impact parameter (DCA) of charged particles
// Study resolution and pull: prepare for beauty study
//
// Authors:
//   Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//   Carlo Bombonati <carlo.bombonati@cern.ch>
//

#include "TMath.h"
#include "TH1F.h"
#include "TList.h"
#include <TParticle.h>
#include "AliMCParticle.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliESDVertex.h"

#include "AliPID.h"

#include "AliHFEdca.h"


ClassImp(AliHFEdca)

//________________________________________________________________________
const Char_t* AliHFEdca::fgkParticles[12] = {
 // particles name
 "electron", "muonMinus","pionMinus", "kaonMinus", "protonMinus", 
 "positron", "muonPlus", "pionPlus", "kaonPlus", "protonPlus",
 "allNegative", "allPositive"
};

const Int_t AliHFEdca::fgkPdgParticle[10] = { 
//   11, 13, -211, -233, -2122,  
//   -11, -13, 211, 233, 2122};
 kPDGelectron, kPDGmuon, -kPDGpion, -kPDGkaon, -kPDGproton, 
 -kPDGelectron, -kPDGmuon, kPDGpion, kPDGkaon, kPDGproton};

//________________________________________________________________________
const Int_t AliHFEdca::fgkColorPart[12] = { 
 // colors assigned to particles
 kRed, kBlue, kGreen+2, kYellow+2, kMagenta, 
 kRed+2, kBlue+2, kGreen+4, kYellow+4, kMagenta+2,
 kBlack, kGray+1
};

//________________________________________________________________________
const Float_t AliHFEdca::fgkPtIntv[51] = {
 // define pT bins
 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 
 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 
 1.70, 1.90, 2.10, 2.30, 2.50, 2.70, 2.90, 3.10, 3.30, 3.50, 
 3.80, 4.10, 4.40, 4.70, 5.00, 5.30, 5.60, 5.90, 6.20, 6.50, 
 7.00, 7.50, 8.00, 9.00, 10.0, 11.0, 12.0, 13.0, 15.0, 18.0, 
 20.00};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkDcaVar[2] = {
 "DcaXY",  "DcaZ"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkDcaVarTitle[2] ={
 ";dca_{xy} [#mum];counts", ";dca_{z} [#mum];counts"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkVertexVar[3] = {
 "VertexX", "VertexY", "VertexZ"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkVertexVarTitle[3] ={
 ";vertex_{x} [#mum];counts", ";vertex_{y} [#mum];counts", ";vertex_{z} [#mum];counts"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkResDcaVar[2] = {
 "deltaDcaXY",  "deltaDcaZ"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkResDcaVarTitle[2] ={
 ";residual #Delta(d_{xy}) [#mum];counts", ";residual #Delta(d_{z}) [#mum];counts"};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkPullDcaVar[2] = {
 "pullDcaXY", "pullDcaZ"
};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkPullDcaVarTitle[2] = {
 ";residual dca_{xy}/(error dca_{xy});counts", 
 ";residual dca_{z}/(error dca_{z});counts"
};

//________________________________________________________________________
const Char_t *AliHFEdca::fgkPullDataDcaVarTitle[2] = {
 ";dca_{xy}^{data}/error dca_{xy};counts", 
 ";dca_{z}^{data}/error dca_{z};counts"
};

//________________________________________________________________________
AliHFEdca::AliHFEdca():
  TObject(),
  fStat(NULL)
{
 // default constructor

 for(Int_t j=0; j<kNParticles; j++){
   fHistMcPid[j] = new TH1F();
   fHistEsdPid[j] = new TH1F();
   fHistDataEsdPid[j] = new TH1F();
 }

 for(Int_t i=0; i<3; i++){
   fHistMCvertex[i] = new TH1F();
   fHistESDvertex[i] = new TH1F();
   fHistDatavertex[i] = new TH1F();   
 }
 
 for(Int_t iEle=0; iEle<2; iEle++){
   fHistDataHfePid[iEle] = new TH1F();
 }

}

//________________________________________________________________________
AliHFEdca::AliHFEdca(const AliHFEdca &ref):
 TObject(ref),
 fStat(ref.fStat)
{
 // copy constructor

 for(Int_t j=0; j<kNParticles; j++){
   fHistMcPid[j] = ref.fHistMcPid[j];
   fHistEsdPid[j] = ref.fHistEsdPid[j];
   fHistDataEsdPid[j] = ref.fHistDataEsdPid[j];
 }

 for(Int_t i=0; i<3; i++){
   fHistMCvertex[i] = ref.fHistMCvertex[i];
   fHistESDvertex[i] = ref.fHistESDvertex[i];
   fHistDatavertex[i] = ref.fHistDatavertex[i];
 }
 
 for(Int_t iEle=0; iEle<2; iEle++){
   fHistDataHfePid[iEle] = ref.fHistDataHfePid[iEle];
 }

}
//_______________________________________________________________________________________________
AliHFEdca&AliHFEdca::operator=(const AliHFEdca &ref)
{
 //
 // Assignment operator
 //

 if(this == &ref) return *this;
 TObject::operator=(ref);
 return *this;

}

//________________________________________________________________________
AliHFEdca::~AliHFEdca()
{
 // default destructor

 for(Int_t j=0; j<kNParticles; j++){
   for(Int_t i=0; i<kNPtBins; i++){
     if(fHistDcaXYRes[j][i]) delete fHistDcaXYRes[j][i];
     if(fHistDcaZRes[j][i]) delete fHistDcaZRes[j][i];

     if(fHistDcaXYPull[j][i]) delete fHistDcaXYPull[j][i];
     if(fHistDcaZPull[j][i]) delete fHistDcaZPull[j][i];

     if(fHistDcaXY[j][i]) delete fHistDcaXY[j][i];
     if(fHistDcaZ[j][i]) delete fHistDcaZ[j][i];

     if(j<(kNParticles-2)){
	if(fHistEPDcaXYRes[j][i]) delete fHistEPDcaXYRes[j][i];
	if(fHistEPDcaZRes[j][i]) delete fHistEPDcaZRes[j][i];
	
	if(fHistEPDcaXYPull[j][i]) delete fHistEPDcaXYPull[j][i];
	if(fHistEPDcaZPull[j][i]) delete fHistEPDcaZPull[j][i];
	
	if(fHistEPDcaXY[j][i]) delete fHistEPDcaXY[j][i];
	if(fHistEPDcaZ[j][i]) delete fHistEPDcaZ[j][i];
     }

     if(fHistKFDcaXY[j][i]) delete fHistKFDcaXY[j][i];
     if(fHistKFDcaZ[j][i]) delete fHistKFDcaZ[j][i];

     if(fHistDataDcaXY[j][i]) delete fHistDataDcaXY[j][i];
     if(fHistDataDcaZ[j][i]) delete fHistDataDcaZ[j][i];
     if(fHistDataWoDcaXY[j][i]) delete fHistDataWoDcaXY[j][i];
     if(fHistDataWoDcaZ[j][i]) delete fHistDataWoDcaZ[j][i];

     if(fHistDataDcaXYPull[j][i]) delete fHistDataDcaXYPull[j][i];
     if(fHistDataDcaZPull[j][i]) delete fHistDataDcaZPull[j][i];
     if(fHistDataWoDcaXYPull[j][i]) delete fHistDataWoDcaXYPull[j][i];
     if(fHistDataWoDcaZPull[j][i]) delete fHistDataWoDcaZPull[j][i];
   }

   if(fHistMcPid[j]) delete fHistMcPid[j];
   if(fHistEsdPid[j]) delete fHistEsdPid[j];
   if(fHistDataEsdPid[j]) delete fHistDataEsdPid[j];
 }

 for(Int_t i=0; i<3; i++){
   if(fHistMCvertex[i]) delete fHistMCvertex[i];
   if(fHistESDvertex[i]) delete fHistESDvertex[i];
   if(fHistDatavertex[i]) delete fHistDatavertex[i];    
 }

 // for the HFEpid
 for(Int_t iEle=0; iEle<2; iEle++){
   for(Int_t iPt=0; iPt<kNPtBins; iPt++){
     if(fHistHPDcaXYRes[iEle][iPt]) delete fHistHPDcaXYRes[iEle][iPt]; 
     if(fHistHPDcaZRes[iEle][iPt]) delete fHistHPDcaZRes[iEle][iPt];   
     if(fHistHPDcaXYPull[iEle][iPt]) delete fHistHPDcaXYPull[iEle][iPt]; 
     if(fHistHPDcaZPull[iEle][iPt]) delete fHistHPDcaZPull[iEle][iPt];  
     if(fHistHPDcaXY[iEle][iPt]) delete fHistHPDcaXY[iEle][iPt];     
     if(fHistHPDcaZ[iEle][iPt]) delete fHistHPDcaZ[iEle][iPt];      
     
     
     // Data
     if(fHistHPDataDcaXY[iEle][iPt]) delete fHistHPDataDcaXY[iEle][iPt];   
     if(fHistHPDataDcaZ[iEle][iPt]) delete fHistHPDataDcaZ[iEle][iPt];   
     if(fHistHPDataDcaXYPull[iEle][iPt]) delete fHistHPDataDcaXYPull[iEle][iPt];   
     if(fHistHPDataDcaZPull[iEle][iPt]) delete fHistHPDataDcaZPull[iEle][iPt];    

   }
   for(Int_t i=0; i<2; i++)
     if(fHistHfePid[iEle][i]) delete fHistHfePid[iEle][i];

   if(fHistDataHfePid[iEle]) delete fHistDataHfePid[iEle];

 }

 if(fStat) delete fStat;

 //Printf("analysis done\n");

}

//________________________________________________________________________
void AliHFEdca::InitAnalysis(){

 //Printf("initialize analysis\n");

}


//________________________________________________________________________
void AliHFEdca::PostAnalysis() const
{
 // do fit
 // moved to dcaPostAnalysis.C

}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsResidual(TList *residualList){
 // define histogram
 // 1. residual

 // for residuals
  // fHistDcaXYRes[kNParticles][kNPtBins]=0x0;
  // fHistDcaZRes[kNParticles][kNPtBins]=0x0;
  // fHistEPDcaXYRes[kNParticles-2][kNPtBins]=0x0;
 // fHistEPDcaZRes[kNParticles-2][kNPtBins]=0x0;

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 1000.;
 const Float_t maxZBin = 1000.;


 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkResDcaVarTitle[k]);

   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);
	histName += Form("_MCpid_%s_pT-%.2f-%.2f", (const char*)fgkResDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	TString histEPName((const char*)fgkParticles[j]);
	histEPName += Form("_ESDpid_%s_pT-%.2f-%.2f", (const char*)fgkResDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDcaXYRes[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistDcaXYRes[j][i]->SetLineColor((int)fgkColorPart[j]);
	  if(j<(kNParticles-2)){
	    fHistEPDcaXYRes[j][i] = new TH1F((const char*)histEPName, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	    fHistEPDcaXYRes[j][i]->SetLineColor((int)fgkColorPart[j]);}
	}	    
	if(k==1){
	  fHistDcaZRes[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistDcaZRes[j][i]->SetLineColor((int)fgkColorPart[j]);
	  if(j<(kNParticles-2)){
	    fHistEPDcaZRes[j][i] = new TH1F((const char*)histEPName, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	    fHistEPDcaZRes[j][i]->SetLineColor((int)fgkColorPart[j]); }
	}   
     } // 50 pt bins
   } //12 nparticles
 } // 2 dca var

 //  TList *fResidualList = 0;
 residualList->SetOwner();
 residualList->SetName("residual");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     residualList->Add(fHistDcaXYRes[iPart][iPtBin]);  
     residualList->Add(fHistDcaZRes[iPart][iPtBin]);  
     if(iPart<(kNParticles-2)){
	residualList->Add(fHistEPDcaXYRes[iPart][iPtBin]);  
	residualList->Add(fHistEPDcaZRes[iPart][iPtBin]);  
     }
   } // loop over pt bins
 }  // loop over particles (pos, neg)



}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsPull(TList *pullList){
 // define histogram
 // 2. pull

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 20.;
 const Float_t maxZBin = 20.;


 // for pull -----------------------------------------------------------------------
 // fHistDcaXYPull[kNParticles][kNPtBins]=0x0;
 // fHistDcaZPull[kNParticles][kNPtBins]=0x0;
 // fHistEPDcaXYPull[kNParticles-2][kNPtBins]=0x0;
 // fHistEPDcaZPull[kNParticles-2][kNPtBins]=0x0;


 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkPullDcaVarTitle[k]);

   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);	
	histName += Form("_MCpid_%s_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	TString histEPName((const char*)fgkParticles[j]);	
	histEPName += Form("_ESDpid_%s_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDcaXYPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxXYBin, 1+maxXYBin);
	  fHistDcaXYPull[j][i]->SetLineColor((int)fgkColorPart[j]);
	  if(j<(kNParticles-2))    {
	    fHistEPDcaXYPull[j][i] = new TH1F((const char*)histEPName, (const char*)histTitle, nBins, 1-maxXYBin, 1+maxXYBin);
	    fHistEPDcaXYPull[j][i]->SetLineColor((int)fgkColorPart[j]);}
	}	    
	if(k==1){
	  fHistDcaZPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxZBin, 1+maxZBin);
	  fHistDcaZPull[j][i]->SetLineColor((int)fgkColorPart[j]);
	  if(j<(kNParticles-2))    {
	    fHistEPDcaZPull[j][i] = new TH1F((const char*)histEPName, (const char*)histTitle, nBins, 1-maxZBin, 1+maxZBin);
	    fHistEPDcaZPull[j][i]->SetLineColor((int)fgkColorPart[j]);}
	}   
     } // 50 pt bins
   } //6 nparticles
 } // 2 dca var

 //  TList *fPullList = 0;
 pullList->SetOwner();
 pullList->SetName("pull");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     pullList->Add(fHistDcaXYPull[iPart][iPtBin]);  
     pullList->Add(fHistDcaZPull[iPart][iPtBin]);  
     if(iPart<(kNParticles-2)){
	pullList->Add(fHistDcaXYPull[iPart][iPtBin]);  
	pullList->Add(fHistDcaZPull[iPart][iPtBin]); }
   } // loop over pt bins
 }  // loop over particles (pos, neg)

}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsDca(TList *dcaList){
 // 
 // define histograms: MC dca
 //

 // statistics
 fStat = 0x0;
 fStat = new TH1I("fStatistics", "allStatistics;ID;counts", 7, -3.5, 3.5);
 fStat->SetMarkerStyle(20); 
 fStat->SetMarkerColor(3); 
 fStat->SetMarkerSize(1); 

 // for dca
 // fHistDcaXY[kNParticles][kNPtBins]=0x0;
 // fHistDcaZ[kNParticles][kNPtBins]=0x0;
 // fHistEPDcaXY[kNParticles-2][kNPtBins]=0x0;
 // fHistEPDcaZ[kNParticles-2][kNPtBins]=0x0;

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 1000.;
 const Float_t maxZBin = 1000.;


 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkDcaVarTitle[k]);

   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);	
	histName += Form("_MCpid_%s_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	TString histNameEP((const char*)fgkParticles[j]);	
	histNameEP += Form("_ESDpid_%s_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);

	if(k==0){
	  fHistDcaXY[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistDcaXY[j][i]->SetLineColor((int)fgkColorPart[j]);
	  
	  if(j<(kNParticles-2)){
	    fHistEPDcaXY[j][i] = new TH1F((const char*)histNameEP, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	    fHistEPDcaXY[j][i]->SetLineColor((int)fgkColorPart[j]);}
	}	    
	if(k==1){
	  fHistDcaZ[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistDcaZ[j][i]->SetLineColor((int)fgkColorPart[j]);
	  if(j<(kNParticles-2)){
	    fHistEPDcaZ[j][i] = new TH1F((const char*)histNameEP, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	    fHistEPDcaZ[j][i]->SetLineColor((int)fgkColorPart[j]);}
	}   
     } // 50 pt bins
   } //12 nparticles
 } // 2 dca var

   //  TList *fDcaList = 0;
 dcaList->SetOwner();
 dcaList->SetName("mcDcaDistr");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     dcaList->Add(fHistDcaXY[iPart][iPtBin]);  
     dcaList->Add(fHistDcaZ[iPart][iPtBin]);  
     if(iPart<(kNParticles-2)) {
	dcaList->Add(fHistEPDcaXY[iPart][iPtBin]);  
	dcaList->Add(fHistEPDcaZ[iPart][iPtBin]); }
   } // loop over pt bins
 }  // loop over particles (pos, neg)

 dcaList->Add(fStat);

}

//________________________________________________________________________
void AliHFEdca::CreateHistogramsKfDca(TList *kfDcaList){
 // 
 // define histograms: MC dca
 //

 // statistics
 fStat = 0x0;
 fStat = new TH1I("fStatistics", "allStatistics;ID;counts", 7, -3.5, 3.5);
 fStat->SetMarkerStyle(20); 
 fStat->SetMarkerColor(3); 
 fStat->SetMarkerSize(1); 

 // for kf dca
 // fHistKFDcaXY[kNParticles][kNPtBins]=0x0;
 // fHistKFDcaZ[kNParticles][kNPtBins]=0x0;

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 1000.;
 const Float_t maxZBin = 1000.;


 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkDcaVarTitle[k]);

   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	TString histNameKF((const char*)fgkParticles[j]);	
	histNameKF += Form("_MCpid_KF%s_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistKFDcaXY[j][i] = new TH1F((const char*)histNameKF, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistKFDcaXY[j][i]->SetLineColor((int)fgkColorPart[j]);
	}	    
	if(k==1){
	  fHistKFDcaZ[j][i] = new TH1F((const char*)histNameKF, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistKFDcaZ[j][i]->SetLineColor((int)fgkColorPart[j]);
	}   
     } // 50 pt bins
   } //12 nparticles
 } // 2 dca var

 kfDcaList->SetOwner();
 kfDcaList->SetName("mcKfDcaDistr");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     kfDcaList->Add(fHistKFDcaXY[iPart][iPtBin]);  
     kfDcaList->Add(fHistKFDcaZ[iPart][iPtBin]);  
   } // loop over pt bins
 }  // loop over particles (pos, neg)

 kfDcaList->Add(fStat);

}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsDataDca(TList *dataDcaList){
 //
 // define histograms: real Data
 //

 // for dca
//  fHistDataDcaXY[kNParticles][kNPtBins]=0x0;
//  fHistDataDcaZ[kNParticles][kNPtBins]=0x0;
//  fHistDataWoDcaXY[kNParticles][kNPtBins]=0x0;
//  fHistDataWoDcaZ[kNParticles][kNPtBins]=0x0;

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 1000.;
 const Float_t maxZBin = 1000.;

 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkDcaVarTitle[k]);    
   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);	
	histName += Form("_%s_Data_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);

	TString histNameWo((const char*)fgkParticles[j]);	
	histNameWo += Form("_%s_Data_wo_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDataDcaXY[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistDataDcaXY[j][i]->SetLineColor((int)fgkColorPart[j]);

	  fHistDataWoDcaXY[j][i] = new TH1F((const char*)histNameWo, (const char*)histTitle, nBins, -maxXYBin, maxXYBin);
	  fHistDataWoDcaXY[j][i]->SetLineColor((int)fgkColorPart[j]);
	}	    
	if(k==1){
	  fHistDataDcaZ[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistDataDcaZ[j][i]->SetLineColor((int)fgkColorPart[j]);

	  fHistDataWoDcaZ[j][i] = new TH1F((const char*)histNameWo, (const char*)histTitle, nBins, -maxZBin, maxZBin);
	  fHistDataWoDcaZ[j][i]->SetLineColor((int)fgkColorPart[j]);
	}   
     } // 50 pt bins
   } //12 nparticles
 } // 2 dca var

   //  TList *fDcaList = 0;
 dataDcaList->SetOwner();
 dataDcaList->SetName("dataDcaDistr");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     dataDcaList->Add(fHistDataDcaXY[iPart][iPtBin]);  
     dataDcaList->Add(fHistDataDcaZ[iPart][iPtBin]);  

     dataDcaList->Add(fHistDataWoDcaXY[iPart][iPtBin]);  
     dataDcaList->Add(fHistDataWoDcaZ[iPart][iPtBin]);  
   } // loop over pt bins
 }  // loop over particles (pos, neg)


}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsDataPull(TList *dataPullList){
 // define histogram
 // 2. pull

 const Int_t nBins = 1000;
 const Float_t maxXYBin = 20.;
 const Float_t maxZBin = 20.;

 // for pull -----------------------------------------------------------------------
//  fHistDataDcaXYPull[kNParticles][kNPtBins]=0x0;
//  fHistDataDcaZPull[kNParticles][kNPtBins]=0x0;

//  fHistDataWoDcaXYPull[kNParticles][kNPtBins]=0x0;
//  fHistDataWoDcaZPull[kNParticles][kNPtBins]=0x0;


 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitle((const char*)fgkPullDataDcaVarTitle[k]);

   for(Int_t j=0; j<kNParticles; j++){
     for(Int_t i=0; i<kNPtBins; i++){
	
	TString histName((const char*)fgkParticles[j]);	
	histName += Form("_%s_Data_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);

	TString histNameWo((const char*)fgkParticles[j]);	
	histNameWo += Form("_%s_Data_wo_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistDataDcaXYPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxXYBin, 1+maxXYBin);
	  fHistDataDcaXYPull[j][i]->SetLineColor((int)fgkColorPart[j]);

	  fHistDataWoDcaXYPull[j][i] = new TH1F((const char*)histNameWo, (const char*)histTitle, nBins, 1-maxXYBin, 1+maxXYBin);
	  fHistDataWoDcaXYPull[j][i]->SetLineColor((int)fgkColorPart[j]);
	}	    
	if(k==1){
	  fHistDataDcaZPull[j][i] = new TH1F((const char*)histName, (const char*)histTitle, nBins, 1-maxZBin, 1+maxZBin);
	  fHistDataDcaZPull[j][i]->SetLineColor((int)fgkColorPart[j]);

	  fHistDataWoDcaZPull[j][i] = new TH1F((const char*)histNameWo, (const char*)histTitle, nBins, 1-maxZBin, 1+maxZBin);
	  fHistDataWoDcaZPull[j][i]->SetLineColor((int)fgkColorPart[j]);
	}   
     } // 50 pt bins
   } //6 nparticles
 } // 2 dca var

 //  TList *fDataPullList = 0;
 dataPullList->SetOwner();
 dataPullList->SetName("dataPull");
 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     dataPullList->Add(fHistDataDcaXYPull[iPart][iPtBin]);  
     dataPullList->Add(fHistDataDcaZPull[iPart][iPtBin]);  

     dataPullList->Add(fHistDataWoDcaXYPull[iPart][iPtBin]);  
     dataPullList->Add(fHistDataWoDcaZPull[iPart][iPtBin]);  
   } // loop over pt bins
 }  // loop over particles (pos, neg)

}

//________________________________________________________________________
void AliHFEdca::CreateHistogramsVertex(TList *vertexList){
 //
 // define histograms: vertex
 //
 // for  vertex

//  fHistMCvertex[kNVertexVar]=0x0;
//  fHistESDvertex[kNVertexVar]=0x0;

 const Int_t nBins = 1000;
 const Float_t minXBin = -0.2e4;
 const Float_t maxXBin = 0.2e4;
 const Float_t minYBin = -0.5e4;
 const Float_t maxYBin = 0.5e4;
 const Float_t minZBin = -1.5e5;
 const Float_t maxZBin = 1.5e5;

 const Float_t minBin[kNVertexVar] = {minXBin, minYBin, minZBin};
 const Float_t maxBin[kNVertexVar] = {maxXBin, maxYBin, maxZBin};

 for(Int_t k=0; k<kNVertexVar; k++){
   TString histTitle((const char*)fgkVertexVarTitle[k]);    
   TString histNameMC((const char*)fgkVertexVar[k]);
   histNameMC += Form("_MC");
   TString histNameESD((const char*)fgkVertexVar[k]);
   histNameESD += Form("_ESD");

   fHistMCvertex[k] = new TH1F((const char*)histNameMC, (const char*)histTitle, nBins, minBin[k], maxBin[k]);
   fHistMCvertex[k]->SetLineColor(k+2);

   fHistESDvertex[k] = new TH1F((const char*)histNameESD, (const char*)histTitle, nBins, minBin[k], maxBin[k]);
   fHistESDvertex[k]->SetLineColor(k+2);
 } // 3 vertex var

 vertexList->SetOwner();
 vertexList->SetName("vertexDistr");

 for(Int_t k=0; k<kNVertexVar; k++){
   vertexList->Add(fHistMCvertex[k]);  
   vertexList->Add(fHistESDvertex[k]);  
 }

}



//________________________________________________________________________
void AliHFEdca::CreateHistogramsDataVertex(TList *dataVertexList){
 //
 // define histograms: vertex
 //
 // for data vertex

//  fHistDatavertex[kNVertexVar]=0x0;

 const Int_t nBins = 1000;
 const Float_t minXBin = -0.2e4;
 const Float_t maxXBin = 0.2e4;
 const Float_t minYBin = -0.5e4;
 const Float_t maxYBin = 0.5e4;
 const Float_t minZBin = -1.5e5;
 const Float_t maxZBin = 1.5e5;

 const Float_t minBin[kNVertexVar] = {minXBin, minYBin, minZBin};
 const Float_t maxBin[kNVertexVar] = {maxXBin, maxYBin, maxZBin};

 for(Int_t k=0; k<kNVertexVar; k++){
   TString histTitle((const char*)fgkVertexVarTitle[k]);    
   TString histNameDataESD((const char*)fgkVertexVar[k]);
   histNameDataESD += Form("_data");

   fHistDatavertex[k] = new TH1F((const char*)histNameDataESD, (const char*)histTitle, nBins, minBin[k], maxBin[k]);
   fHistDatavertex[k]->SetLineColor(k+2);
 } // 3 vertex var

   //  TList *fVDaraVertexList = 0;
 dataVertexList->SetOwner();
 dataVertexList->SetName("dataVertexDistr");

 for(Int_t k=0; k<kNVertexVar; k++){
   dataVertexList->Add(fHistDatavertex[k]);  
 }

}

//_______________________________________________________________________________________________
void AliHFEdca::CreateHistogramsPid(TList *mcPidList){
 //
 // define histograms which fills combined PID
 //

 const Char_t *mcOResd[2]={"mcPt", "esdPt"};

 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   TString histTitleMc((const char*)fgkParticles[iPart]);
   TString histTitleEsd((const char*)fgkParticles[iPart]);
   histTitleMc += Form("_McPid_%s;p_{T} [GeV/c];counts", mcOResd[0]);
   histTitleEsd += Form("_EsdPid_%s;p_{T} [GeV/c];counts", mcOResd[1]);

   TString histNameMc((const char*)fgkParticles[iPart]);
   TString histNameEsd((const char*)fgkParticles[iPart]);
   histNameMc+=Form("_McPid_%s", mcOResd[0]);     
   histNameEsd+=Form("_EsdPid_%s", mcOResd[1]);     

   fHistMcPid[iPart] = new TH1F(histNameMc, histTitleMc, kNPtBins, fgkPtIntv);
   fHistEsdPid[iPart] = new TH1F(histNameEsd, histTitleEsd, kNPtBins, fgkPtIntv);
 }


 mcPidList->SetOwner();
 mcPidList->SetName("combinedPid");

 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   mcPidList->Add(fHistMcPid[iPart]);
   mcPidList->Add(fHistEsdPid[iPart]);
 }
}


//_______________________________________________________________________________________________
void AliHFEdca::CreateHistogramsDataPid(TList *pidList){
 //
 // define histograms which fills combined PID: data
 //



 for(Int_t iPart=0; iPart<kNParticles; iPart++){
     TString histTitleEsd((const char*)fgkParticles[iPart]);
     histTitleEsd+=Form("_DataEsdPid_esdPt;p_{T} [GeV/c];counts");
     TString histNameEsd((const char*)fgkParticles[iPart]);
     histNameEsd+=Form("_DataEsdPid");

     fHistDataEsdPid[iPart] = new TH1F(histNameEsd, histTitleEsd, kNPtBins, fgkPtIntv);
 }


 pidList->SetOwner();
 pidList->SetName("dataCombinedPid");

 for(Int_t iPart=0; iPart<kNParticles; iPart++)
   pidList->Add(fHistDataEsdPid[iPart]);


}



//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsDca(AliESDEvent * const esdEvent, AliESDtrack * const track, AliMCEvent * const mcEvent)
{
 // the kDca plugin
 // MC vertex
 AliMCVertex *mcPrimVtx = (AliMCVertex *)mcEvent->GetPrimaryVertex();      
 Double_t mcPrimV[3];
 mcPrimV[0] = mcPrimVtx->GetX();
 mcPrimV[1] = mcPrimVtx->GetY();
 mcPrimV[2] = mcPrimVtx->GetZ();

 Double_t mcVtxXY = TMath::Abs(mcPrimV[0]*mcPrimV[0] + mcPrimV[1]*mcPrimV[1]);

// filling historgams track by track
// obtaining reconstructed dca ------------------------------------------------------------------

 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  

// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();

 Float_t magneticField = 0;  // initialized as 5kG
 magneticField = esdEvent->GetMagneticField();  // in kG

 Double_t beampiperadius=3.;
 Double_t dz[2];   // error of dca in cm
 Double_t covardz[3];

 if(!track->PropagateToDCA(primVtx,magneticField, beampiperadius, dz, covardz)) return;  // protection

 AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel()))); 
 if(!mctrack) return;
 TParticle *part = mctrack->Particle();

 Float_t vx = part->Vx();  // in cm
 Float_t vy = part->Vy();  // in cm
 Float_t vz = part->Vz();   // in cm

 Float_t vxy = TMath::Sqrt(vx*vx+vy*vy);

 Float_t mcpx = part->Px();
 Float_t mcpy = part->Py();
 Float_t mcpt = TMath::Sqrt(mcpx*mcpx+mcpy*mcpy);

 Int_t pdg = part->GetPdgCode();
 Int_t esdPid = GetCombinedPid(track);

 Int_t charge = 1;
 if(pdg==kPDGelectron || pdg==kPDGmuon 
    || pdg==-kPDGpion || pdg==-kPDGkaon || pdg==-kPDGproton) charge = -1;  

 // calculate mcDca ------------------------------------------------------------------ 
 const Float_t conv[2] = {1.783/1.6, 2.99792458};
 Float_t radiusMc = mcpt/(TMath::Abs(magneticField)/10.)*conv[0]*conv[1]; // pt in GeV/c, magnetic field in Tesla, radius in meter

 Float_t nx = esdpx/mcpt;
 Float_t ny = esdpy/mcpt;

 Float_t radius;
 radius = TMath::Abs(radiusMc);

 Double_t dxy = vxy - mcVtxXY;   // in cm
 Double_t dvx = vx - mcPrimV[0]; // in cm
 Double_t dvy = vy - mcPrimV[1]; // in cm

 Float_t mcDcaXY = (radius - TMath::Sqrt(dxy*dxy/100./100. + radius*radius + 2*radius*charge*(dvx*ny-dvy*nx)/100.)) ;  // in meters

 Double_t mcDca[2] = {mcDcaXY*100, vz};  // in cm
 Double_t residual[2] = {0, 0};
 Double_t pull[2] = {0, 0};
 Double_t error[2] ={TMath::Sqrt(covardz[0]), TMath::Sqrt(covardz[2])};
 for(Int_t i=0; i<2; i++){
   residual[i] = dz[i] - mcDca[i]; // in centimeters       
   if(error[i]!=0)pull[i] = residual[i]/error[i];   // unitless
 }


 for(Int_t iPart=0; iPart<(kNParticles-2); iPart++){    
   // identified ones
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
	if(pdg==fgkPdgParticle[iPart]) {
	  fHistDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1.0e4);  
	  fHistDcaZRes[iPart][iPtBin]->Fill(residual[1]*1.0e4);   
	  fHistDcaXYPull[iPart][iPtBin]->Fill(pull[0]);  
	  fHistDcaZPull[iPart][iPtBin]->Fill(pull[1]); 
	  fHistDcaXY[iPart][iPtBin]->Fill(dz[0]*1.0e4); 
	  fHistDcaZ[iPart][iPtBin]->Fill(dz[1]*1.0e4);  
	}  // mc pdg  	  
	
	if(esdPid==fgkPdgParticle[iPart]) {
	  fHistEPDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1.0e4);  
	  fHistEPDcaZRes[iPart][iPtBin]->Fill(residual[1]*1.0e4);   
	  fHistEPDcaXYPull[iPart][iPtBin]->Fill(pull[0]);  
	  fHistEPDcaZPull[iPart][iPtBin]->Fill(pull[1]); 
	  fHistEPDcaXY[iPart][iPtBin]->Fill(dz[0]*1.0e4); 
	  fHistEPDcaZ[iPart][iPtBin]->Fill(dz[1]*1.0e4);  
	}  // esd pid
	
     } // pt range

     else
	continue;
   }  // pt loop
 } // particle id loop

 // for charged particles: no pid
 for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
   if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
     Int_t iPart = 10;
     if(charge>0) iPart = 11;
     fHistDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1e4);
     fHistDcaZRes[iPart][iPtBin]->Fill(residual[1]*1e4);
     fHistDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
     fHistDcaZPull[iPart][iPtBin]->Fill(pull[1]);
     fHistDcaXY[iPart][iPtBin]->Fill(dz[0]*1e4);
     fHistDcaZ[iPart][iPtBin]->Fill(dz[1]*1e4);      
   }
   else
     continue;
 } // pt  

}

//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsKfDca(AliESDEvent * const esdEvent, AliESDtrack * const track, const AliMCEvent * const mcEvent)
 {
 // the kKfDca plugin

// filling historgams track by track

// obtaining reconstructed dca ------------------------------------------------------------------  
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
   Float_t magneticField = 0;  // initialized as 5kG
 magneticField = esdEvent->GetMagneticField();  // in kG

 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  

 Int_t charge = track->Charge();

 Double_t beampiperadius=3.;  
 Double_t dz[2];   // error of dca in cm
 Double_t covardz[3];
 if(!track->PropagateToDCA(primVtx,magneticField, beampiperadius, dz, covardz)) return; // protection 

 AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));  
 if(!mctrack) return;
 TParticle *part = mctrack->Particle();    
 Int_t pdg = part->GetPdgCode();

 // calculate dca using AliKFParticle class------------------------------------------------------------------
 Double_t kfdz[3] = {0, 0, 0};  
 Double_t kfdzwith[3] = {0, 0, 0};  

 Int_t trkID = track->GetID();

 AliKFParticle::SetField(magneticField);
 AliKFParticle kfParticle(*track, pdg);

 // prepare kfprimary vertex
 AliKFVertex kfESDprimary;
 // Reconstruct Primary Vertex (with ESD tracks)
 Int_t n=primVtx->GetNIndices();
 if (n>0 && primVtx->GetStatus()){
   kfESDprimary = AliKFVertex(*primVtx);

   Double_t dcaXYWithTrk = kfParticle.GetDistanceFromVertexXY(kfESDprimary);
   Double_t dcaWithTrk = kfParticle.GetDistanceFromVertex(kfESDprimary);
   Double_t dcaZWithTrk = 0;
   if(TMath::Abs(dcaWithTrk)>=TMath::Abs(dcaXYWithTrk)) 
     dcaZWithTrk =TMath::Sqrt(dcaWithTrk*dcaWithTrk-dcaXYWithTrk*dcaXYWithTrk)*((dz[1]*-1<=0)?1:-1);
   kfdzwith[0] = dcaXYWithTrk;  
   kfdzwith[1] = dcaZWithTrk; 
   kfdzwith[2] = dcaWithTrk;  // with current track

   Double_t dcaXYWoTrk = 0;
   Double_t dcaZWoTrk = 0;
   Double_t dcaWoTrk = 0;

   UShort_t *priIndex = primVtx->GetIndices();

   for (Int_t i=0;i<n;i++){

     Int_t idx = Int_t(priIndex[i]);
     if (idx == trkID){
	kfESDprimary -= kfParticle;
	dcaXYWoTrk = kfParticle.GetDistanceFromVertexXY(kfESDprimary);
	dcaWoTrk = kfParticle.GetDistanceFromVertex(kfESDprimary);
	if((dcaWoTrk-dcaXYWoTrk)>=0)
	  dcaZWoTrk = TMath::Abs(dcaWoTrk*dcaWoTrk - dcaXYWoTrk*dcaXYWoTrk)*((dz[1]*-1<=0)?1:-1);
     }  // remove current track from this calculation
   }  // loop over all primary vertex contributors


   kfdz[0] = dcaXYWoTrk;  
   kfdz[1] = dcaZWoTrk;
   kfdz[2] = dcaWoTrk;  

 }  // only if n contributor > 0 and primVtx constructed

 fStat->Fill(0);

 if(dz[0]!=0 && dz[0]*kfdzwith[0]>0 && TMath::Abs(kfdzwith[0]/dz[0])>0.9999 && TMath::Abs(kfdzwith[0]/dz[0])<1.0001)fStat->Fill(1);; // same 
 if(dz[0]!=0 && dz[0]*kfdzwith[0]<0 && TMath::Abs(kfdzwith[0]/dz[0])>0.9999 && TMath::Abs(kfdzwith[0]/dz[0])<1.0001) fStat->Fill(2); // swapped sign
 if(kfdzwith[0]==0 && dz[0]!=0) fStat->Fill(3);  // 0 from KF particle (with current track)

 if(dz[0]!=0 && dz[0]*kfdz[0]>0 && TMath::Abs(kfdz[0]/dz[0])>0.8 && TMath::Abs(kfdz[0]/dz[0])<1.2) fStat->Fill(-1);; // same 
 if(dz[0]!=0 && dz[0]*kfdz[0]<0 && TMath::Abs(kfdz[0]/dz[0])>0.8 && TMath::Abs(kfdz[0]/dz[0]) <1.2) fStat->Fill(-2); // swapped sign
 if(kfdz[0]==0 && dz[0]!=0) fStat->Fill(-3);  // 0 from KF particle (without current track)

 for(Int_t iPart=0; iPart<(kNParticles-2); iPart++){    
   // identified ones
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
	if(pdg==fgkPdgParticle[iPart]) {
	  fHistKFDcaXY[iPart][iPtBin]->Fill(kfdzwith[0]*1.0e4); 
	  fHistKFDcaZ[iPart][iPtBin]->Fill(kfdzwith[1]*1.0e4);  
	}  // mc pdg  	  
     } // pt range

     else
	continue;
   }  // pt loop
 } // particle id loop

 // for charged particles: no pid
 for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
   if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
     Int_t iPart = 10;
     if(charge>0) iPart = 11;
     fHistKFDcaXY[iPart][iPtBin]->Fill(kfdzwith[0]*1e4);
     fHistKFDcaZ[iPart][iPtBin]->Fill(kfdzwith[1]*1e4);

   }
   else
     continue;
 } // pt  

}  // KF dca


//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsVtx(AliESDEvent *const esdEvent, AliMCEvent *const mcEvent)
{

 // MC vertex
 AliMCVertex *mcPrimVtx = (AliMCVertex *)mcEvent->GetPrimaryVertex();      
 Double_t mcPrimV[3];
 mcPrimV[0] = mcPrimVtx->GetX();
 mcPrimV[1] = mcPrimVtx->GetY();
 mcPrimV[2] = mcPrimVtx->GetZ();

// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();

 for(Int_t i=0; i<kNVertexVar; i++){
   fHistMCvertex[i]->Fill(mcPrimV[i]*1.0e4); 
   fHistESDvertex[i]->Fill(primV[i]*1.0e4); 
 }

}

//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsPid(AliESDtrack * const track, const AliMCEvent * const mcEvent)
{


// filling historgams track by track
 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);    

 AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));  
 if(!mctrack) return;
 TParticle *part = mctrack->Particle();

 Float_t mcpx = part->Px();
 Float_t mcpy = part->Py();
 Float_t mcpt = TMath::Sqrt(mcpx*mcpx+mcpy*mcpy);

 Int_t pdg = part->GetPdgCode();
 Int_t esdPid = GetCombinedPid(track);


 Double_t ptMom[2] = {mcpt, esdpt};
 // for combined PID
 for(Int_t iPart=0; iPart<(kNParticles-2); iPart++){
   if(pdg==fgkPdgParticle[iPart])                // pid all by MC 
     fHistMcPid[iPart]->Fill(ptMom[0]);

   if(esdPid==fgkPdgParticle[iPart])             // pid all by combined pid
     fHistEsdPid[iPart]->Fill(ptMom[1]);    
 } // loop over particles

 // for charged
 if(pdg==kPDGelectron || pdg==kPDGmuon || pdg==-kPDGpion || pdg==-kPDGkaon || pdg==-kPDGproton)
   fHistMcPid[10]->Fill(ptMom[0]);
 if(pdg==-kPDGelectron || pdg==-kPDGmuon || pdg==kPDGpion || pdg==kPDGkaon || pdg==kPDGproton)
   fHistMcPid[11]->Fill(ptMom[0]);
 if(esdPid==kPDGelectron || esdPid==kPDGmuon || esdPid==-kPDGpion || esdPid==-kPDGkaon || esdPid==-kPDGproton)
   fHistEsdPid[10]->Fill(ptMom[1]);  
 if(esdPid==-kPDGelectron || esdPid==-kPDGmuon || esdPid==kPDGpion || esdPid==kPDGkaon || esdPid==kPDGproton)
   fHistEsdPid[11]->Fill(ptMom[1]);  
}


////_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsDataDca(AliESDEvent * const esdEvent, AliESDtrack * const track, AliESDVertex * const vtxESDSkip)
{
// filling historgams track by track
// obtaining reconstructed dca --------------------------------------------------------------

 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  
 Int_t charge = track->Charge();

// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();


 Float_t magneticField = 0;  // initialized as 5kG
 magneticField = esdEvent->GetMagneticField();  // in kG

 Double_t beampiperadius=3.;
 Double_t dz[2];   // error of dca in cm
 Double_t covardz[3];

 if(!track->PropagateToDCA(primVtx,magneticField, beampiperadius, dz, covardz)) return;  // protection


 Double_t pull[2] = {0, 0};
 Double_t error[2] ={TMath::Sqrt(covardz[0]), TMath::Sqrt(covardz[2])};
 for(Int_t i=0; i<2; i++){
   if(error[i]!=0)pull[i] = dz[i]/error[i];   // unitless                                                
 }

 // get dca when current track is not included

 Double_t dzwo[2], covardzwo[3];
 Double_t pullwo[2] = {0, 0};
 if(!track->PropagateToDCA(vtxESDSkip, magneticField, beampiperadius, dzwo, covardzwo)) return;   // protection

 Double_t errorwo[2] ={TMath::Sqrt(TMath::Abs(covardzwo[0])), TMath::Sqrt(TMath::Abs(covardzwo[2]))};
 for(Int_t i=0; i<2; i++){
   if(errorwo[i]!=0) pullwo[i] = dzwo[i]/errorwo[i];   // unitless                                                
 }

 // do pid 
 Int_t esdPid = GetCombinedPid(track);

 for(Int_t iPart=0; iPart<(kNParticles-2); iPart++){    
   // identified ones
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     if(esdPid==fgkPdgParticle[iPart] && (esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1])) {
	fHistDataDcaXY[iPart][iPtBin]->Fill(dz[0]*1e4);
	fHistDataDcaZ[iPart][iPtBin]->Fill(dz[1]*1e4);      
	fHistDataDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
	fHistDataDcaZPull[iPart][iPtBin]->Fill(pull[1]);	
	// w/o current track
	fHistDataWoDcaXY[iPart][iPtBin]->Fill(dzwo[0]*1e4);
	fHistDataWoDcaZ[iPart][iPtBin]->Fill(dzwo[1]*1e4);      
	fHistDataWoDcaXYPull[iPart][iPtBin]->Fill(pullwo[0]);
	fHistDataWoDcaZPull[iPart][iPtBin]->Fill(pullwo[1]);	
     }
     else
	continue;
   }
 } 

 // for charged particles
 for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
   if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
     Int_t iPart = 10;
     if(charge>0) iPart = 11;
     fHistDataDcaXY[iPart][iPtBin]->Fill(dz[0]*1e4);
     fHistDataDcaZ[iPart][iPtBin]->Fill(dz[1]*1e4);      
     fHistDataDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
     fHistDataDcaZPull[iPart][iPtBin]->Fill(pull[1]);
     // without current track
     fHistDataWoDcaXY[iPart][iPtBin]->Fill(dzwo[0]*1e4);
     fHistDataWoDcaZ[iPart][iPtBin]->Fill(dzwo[1]*1e4);      
     fHistDataWoDcaXYPull[iPart][iPtBin]->Fill(pullwo[0]);
     fHistDataWoDcaZPull[iPart][iPtBin]->Fill(pullwo[1]);	

   }
   else
     continue;
 } 

}

//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsDataVtx(AliESDEvent * const esdEvent)
{


// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();

 // require events with at least 3 contributors for primary vertex construction
 Int_t nEsdPrimVtxCtrb = primVtx->GetNContributors(); 
 if(nEsdPrimVtxCtrb<1) return; // for pass 1, no diomond constrain, each event has at least 1 contributor to Vtx
 for(Int_t i=0; i<kNVertexVar; i++)
   fHistDatavertex[i]->Fill(primV[i]*1e4);

}


////_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsDataPid(AliESDtrack * const track)
{
// filling historgams track by track
// obtaining reconstructed dca --------------------------------------------------------------

 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  
 Int_t charge = track->Charge();

 Int_t esdPid = GetCombinedPid(track);

 for(Int_t iPart=0; iPart<kNParticles; iPart++){
   if(iPart<(kNParticles-2)){
     if(esdPid==fgkPdgParticle[iPart])  fHistDataEsdPid[iPart]->Fill(esdpt);
   }     // for identified

   else {
     if(charge<0) fHistDataEsdPid[10]->Fill(esdpt);
     if(charge>0) fHistDataEsdPid[11]->Fill(esdpt);
   }
 }
}  

//_________________________________________________________________________________________________
void AliHFEdca::ApplyExtraCuts(AliESDEvent * const esdEvent, Int_t nMinPrimVtxContributor)
{ 

 //
 // only one extra cut, number of contributors to each primary vertex
 //

 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Int_t nEsdPrimVtxCtrb = primVtx->GetNContributors(); 
 if(nEsdPrimVtxCtrb<nMinPrimVtxContributor) return; 
 // for pass 1, no diomond constrain, each event has at least 1 contributor to Vtx

}

//_____________________________________________________
Int_t AliHFEdca::GetCombinedPid(AliESDtrack *const track) 
{

 // combined detector pid             
 Double_t prob[AliPID::kSPECIES];
 track->GetESDpid(prob);

 // setting priors!
 Double_t priors[AliPID::kSPECIESN];
 priors[0] = 0.01;
 priors[1] = 0.01;
 priors[2] = 0.85;
 priors[3] = 0.10;
 priors[4] = 0.05;

 Int_t charge = track->Charge();
 Int_t esdPid = -1; 

 AliPID pid;
 pid.SetPriors(priors);
 pid.SetProbabilities(prob);

 // identify particle as the most probable 

 Double_t pelectron = pid.GetProbability(AliPID::kElectron);
 if(pelectron > pid.GetProbability(AliPID::kMuon) && 
    pelectron > pid.GetProbability(AliPID::kPion) && 
    pelectron > pid.GetProbability(AliPID::kKaon) && 
    pelectron > pid.GetProbability(AliPID::kProton) )  esdPid = -kPDGelectron; 

 Double_t pmuon = pid.GetProbability(AliPID::kMuon);
 if(pmuon > pid.GetProbability(AliPID::kElectron) && 
    pmuon > pid.GetProbability(AliPID::kPion) && 
    pmuon > pid.GetProbability(AliPID::kKaon) && 
    pmuon > pid.GetProbability(AliPID::kProton) )  esdPid = -kPDGmuon; 

 Double_t ppion = pid.GetProbability(AliPID::kPion);
 if(ppion > pid.GetProbability(AliPID::kElectron) && 
    ppion > pid.GetProbability(AliPID::kMuon) && 
    ppion > pid.GetProbability(AliPID::kKaon) && 
    ppion > pid.GetProbability(AliPID::kProton) )  esdPid = kPDGpion; 

 Double_t pkaon = pid.GetProbability(AliPID::kKaon);
 if(pkaon > pid.GetProbability(AliPID::kElectron) && 
    pkaon > pid.GetProbability(AliPID::kMuon) && 
    pkaon > pid.GetProbability(AliPID::kPion) && 
    pkaon > pid.GetProbability(AliPID::kProton) )  esdPid = kPDGkaon; 

 Double_t pproton = pid.GetProbability(AliPID::kProton);
 if(pproton > pid.GetProbability(AliPID::kElectron) && 
    pproton > pid.GetProbability(AliPID::kMuon) && 
    pproton > pid.GetProbability(AliPID::kPion) && 
    pproton > pid.GetProbability(AliPID::kKaon) )  esdPid = kPDGproton; 


 return charge*esdPid;

}


// for the HFE pid

//________________________________________________________________________
void AliHFEdca::CreateHistogramsHfeDca(TList *hfeDcaList){
 //
 // define histograms: hfe pid electrons in MC
 //

 const Int_t nBinsDca = 1000;
 const Float_t maxXYBinDca = 1000.;
 const Float_t maxZBinDca = 1000.;

 const Int_t nBinsPull = 1000;
 const Float_t maxXYBinPull = 20.;
 const Float_t maxZBinPull = 20.;

 const Char_t *mcOResd[2]={"mcPt", "esdPt"};

//  fHistHPDcaXY[2][kNPtBins]=0x0;
//  fHistHPDcaZ[2][kNPtBins]=0x0;
//  fHistHPDcaXYRes[2][kNPtBins]=0x0;
//  fHistHPDcaZRes[2][kNPtBins]=0x0;
//  fHistHPDcaXYPull[2][kNPtBins]=0x0;
//  fHistHPDcaZPull[2][kNPtBins]=0x0;


 for(Int_t k=0; k<kNDcaVar; k++){

   TString histTitleDca((const char*)fgkDcaVarTitle[k]);
   TString histTitleRes((const char*)fgkPullDcaVarTitle[k]);
   TString histTitlePull((const char*)fgkResDcaVarTitle[k]);

   for(Int_t iPart=0; iPart<2; iPart++){
     for(Int_t i=0; i<kNPtBins; i++){		
	TString histHPName((const char*)fgkParticles[iPart*5]);
	histHPName += Form("_HFEpid_%s_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	TString histHPNameRes((const char*)fgkParticles[iPart*5]);
	histHPNameRes += Form("_HFEpid_%s_pT-%.2f-%.2f", (const char*)fgkResDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	TString histHPNamePull((const char*)fgkParticles[iPart*5]);
	histHPNamePull += Form("_HFEpid_%s_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistHPDcaXY[iPart][i] = new TH1F((const char*)histHPName, (const char*)histTitleDca, nBinsDca, 1-maxXYBinDca, 1+maxXYBinDca);
	  fHistHPDcaXY[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDcaXYRes[iPart][i] = new TH1F((const char*)histHPNameRes, (const char*)histTitleRes, nBinsDca, 1-maxXYBinDca, 1+maxXYBinDca);
	  fHistHPDcaXYRes[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDcaXYPull[iPart][i] = new TH1F((const char*)histHPNamePull, (const char*)histTitlePull, nBinsPull, 1-maxXYBinPull, 1+maxXYBinPull);
	  fHistHPDcaXYPull[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	}
	
	if(k==1){
	  fHistHPDcaZ[iPart][i] = new TH1F((const char*)histHPName, (const char*)histTitleDca, nBinsDca, 1-maxZBinDca, 1+maxZBinDca);
	  fHistHPDcaZ[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDcaZRes[iPart][i] = new TH1F((const char*)histHPNameRes, (const char*)histTitleRes, nBinsDca, 1-maxZBinDca, 1+maxZBinDca);
	  fHistHPDcaZRes[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDcaZPull[iPart][i] = new TH1F((const char*)histHPNamePull, (const char*)histTitlePull, nBinsPull, 1-maxZBinPull, 1+maxZBinPull);
	  fHistHPDcaZPull[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	}
     } // 50 pt bins
   } //2 nparticles
 } // 2 dca var

 // fHistHfePid[2][2] = 0x0; //!  HFE pid  
 for(Int_t id=0; id<2; id++){
   for(Int_t iPart=0; iPart<2; iPart++){
     TString histTitleHfe((const char*)fgkParticles[iPart*5]);
     histTitleHfe+=Form("_MC_HfePid_esdPt;p_{T} [GeV/c];counts");
     TString histNameHfe((const char*)fgkParticles[iPart*5]);
     histNameHfe+=Form("_MC_HfePid_%s", mcOResd[id]);
     fHistHfePid[id][iPart] = new TH1F(histNameHfe, histTitleHfe, kNPtBins, fgkPtIntv);
   }
 }  

 hfeDcaList->SetOwner();
 hfeDcaList->SetName("hfeDca");
 for(Int_t iPart=0; iPart<2; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     hfeDcaList->Add(fHistHPDcaXY[iPart][iPtBin]);  
     hfeDcaList->Add(fHistHPDcaZ[iPart][iPtBin]); 
     hfeDcaList->Add(fHistHPDcaXYRes[iPart][iPtBin]);  
     hfeDcaList->Add(fHistHPDcaZRes[iPart][iPtBin]); 
     hfeDcaList->Add(fHistHPDcaXYPull[iPart][iPtBin]);  
     hfeDcaList->Add(fHistHPDcaZPull[iPart][iPtBin]);
   }
   for(Int_t id=0; id<2; id++)
     hfeDcaList->Add(fHistHfePid[id][iPart]);
 }

}


//________________________________________________________________________
void AliHFEdca::CreateHistogramsHfeDataDca(TList *hfeDataDcaList){
 //
 // define histograms: hfe pid electrons in data
 //

 const Int_t nBinsDca = 1000;
 const Float_t maxXYBinDca = 1000.;
 const Float_t maxZBinDca = 1000.;

 const Int_t nBinsPull = 1000;
 const Float_t maxXYBinPull = 20.;
 const Float_t maxZBinPull = 20.;


//  fHistHPDataDcaXY[2][kNPtBins]=0x0;
//  fHistHPDataDcaZ[2][kNPtBins]=0x0;
//  fHistHPDataDcaXYPull[2][kNPtBins]=0x0;
//  fHistHPDataDcaZPull[2][kNPtBins]=0x0;

 for(Int_t k=0; k<kNDcaVar; k++){
   TString histTitleDca((const char*)fgkDcaVarTitle[k]);
   TString histTitlePull((const char*)fgkPullDcaVarTitle[k]);
   for(Int_t iPart=0; iPart<2; iPart++){
     for(Int_t i=0; i<kNPtBins; i++){		
	TString histHPName((const char*)fgkParticles[iPart*5]);
	histHPName += Form("_HFEpid_%s_pT-%.2f-%.2f", (const char*)fgkDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	TString histHPNamePull((const char*)fgkParticles[iPart*5]);
	histHPNamePull += Form("_HFEpid_%s_pT-%.2f-%.2f", (const char*)fgkPullDcaVar[k], fgkPtIntv[i], fgkPtIntv[i+1]);
	
	if(k==0){
	  fHistHPDataDcaXY[iPart][i] = new TH1F((const char*)histHPName, (const char*)histTitleDca, nBinsDca, 1-maxXYBinDca, 1+maxXYBinDca);
	  fHistHPDataDcaXY[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDataDcaXYPull[iPart][i] = new TH1F((const char*)histHPNamePull, (const char*)histTitlePull, nBinsPull, 1-maxXYBinPull, 1+maxXYBinPull);
	  fHistHPDataDcaXYPull[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	}
	
	if(k==1){
	  fHistHPDataDcaZ[iPart][i] = new TH1F((const char*)histHPName, (const char*)histTitleDca, nBinsDca, 1-maxZBinDca, 1+maxZBinDca);
	  fHistHPDataDcaZ[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	  fHistHPDataDcaZPull[iPart][i] = new TH1F((const char*)histHPNamePull, (const char*)histTitlePull, nBinsDca, 1-maxZBinPull, 1+maxZBinPull);
	  fHistHPDataDcaZPull[iPart][i]->SetLineColor((int)fgkColorPart[iPart*5]);
	}
	
     } // 50 pt bins
   } // 2 particle type
 } // 2 dca var

 //fHistDataHfePid[2] = 0x0; //!  HFE pid  
 for(Int_t iPart=0; iPart<2; iPart++){
   TString histTitleHfe((const char*)fgkParticles[iPart*5]);
   histTitleHfe+=Form("_Data_HfePid_esdPt;p_{T} [GeV/c];counts");
   TString histNameHfe((const char*)fgkParticles[iPart*5]);
   histNameHfe+=Form("_Data_HfePid");
   fHistDataHfePid[iPart] = new TH1F(histNameHfe, histTitleHfe, kNPtBins, fgkPtIntv);
 }


 hfeDataDcaList->SetOwner();
 hfeDataDcaList->SetName("hfeDataDca");
 for(Int_t iPart=0; iPart<2; iPart++){
   for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
     hfeDataDcaList->Add(fHistHPDataDcaXY[iPart][iPtBin]);  
     hfeDataDcaList->Add(fHistHPDataDcaZ[iPart][iPtBin]); 
     hfeDataDcaList->Add(fHistHPDataDcaXYPull[iPart][iPtBin]);  
     hfeDataDcaList->Add(fHistHPDcaZPull[iPart][iPtBin]); 

     hfeDataDcaList->Add(fHistDataHfePid[iPart]);
   }
 }


}


//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsHfeDca(AliESDEvent * const esdEvent, AliESDtrack * const track, AliMCEvent * const mcEvent)
{
 // the kHFEpid plugin

 AliMCVertex *mcPrimVtx = (AliMCVertex *)mcEvent->GetPrimaryVertex();      
 Double_t mcPrimV[3];
 mcPrimV[0] = mcPrimVtx->GetX();
 mcPrimV[1] = mcPrimVtx->GetY();
 mcPrimV[2] = mcPrimVtx->GetZ();

 Double_t mcVtxXY = TMath::Abs(mcPrimV[0]*mcPrimV[0] + mcPrimV[1]*mcPrimV[1]);

// filling historgams track by track
// obtaining reconstructed dca ------------------------------------------------------------------
 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  

// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex();      
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();

 Float_t magneticField = 0;  // initialized as 5kG
 magneticField = esdEvent->GetMagneticField();  // in kG

 Double_t beampiperadius=3.;
 Double_t dz[2];   // error of dca in cm
 Double_t covardz[3];
 if(!track->PropagateToDCA(primVtx,magneticField, beampiperadius, dz, covardz)) return; // protection

 AliMCParticle *mctrack = dynamic_cast<AliMCParticle *>(mcEvent->GetTrack(TMath::Abs(track->GetLabel())));  
 if(!mctrack) return;
 TParticle *part = mctrack->Particle();

 Float_t vx = part->Vx();  // in cm
 Float_t vy = part->Vy();  // in cm
 Float_t vz = part->Vz();   // in cm

 Float_t vxy = TMath::Sqrt(vx*vx+vy*vy);

 Float_t mcpx = part->Px();
 Float_t mcpy = part->Py();
 Float_t mcpt = TMath::Sqrt(mcpx*mcpx+mcpy*mcpy);

 Int_t pdg = part->GetPdgCode();

 Int_t charge = 1;
 if(pdg==kPDGelectron || pdg==kPDGmuon 
    || pdg==-kPDGpion || pdg==-kPDGkaon || pdg==-kPDGproton) charge = -1;

 // calculate mcDca ------------------------------------------------------------------ 
 const Float_t conv[2] = {1.783/1.6, 2.99792458};
 Float_t radiusMc = mcpt/(TMath::Abs(magneticField)/10.)*conv[0]*conv[1]; // pt in GeV/c, magnetic field in Tesla, radius in meter

 Float_t nx = esdpx/mcpt;
 Float_t ny = esdpy/mcpt;

 Float_t radius;
 radius = TMath::Abs(radiusMc);

 Double_t dxy = vxy - mcVtxXY;   // in cm
 Double_t dvx = vx - mcPrimV[0]; // in cm
 Double_t dvy = vy - mcPrimV[1]; // in cm

 Float_t mcDcaXY = (radius - TMath::Sqrt(dxy*dxy/100./100. + radius*radius + 2*radius*charge*(dvx*ny-dvy*nx)/100.)) ;  // in meters

 Double_t mcDca[2] = {mcDcaXY*100, vz};  // in cm
 Double_t residual[2] = {0, 0};
 Double_t pull[2] = {0, 0};
 Double_t error[2] ={TMath::Sqrt(covardz[0]), TMath::Sqrt(covardz[2])};
 for(Int_t i=0; i<2; i++){
   residual[i] = dz[i] - mcDca[i]; // in centimeters       
   if(error[i]!=0)pull[i] = residual[i]/error[i];   // unitless
 }

 Int_t iPart = -1;
 if(track->Charge()<0) iPart = 0;  // electron
 if(track->Charge()>0) iPart = 1;  // positron
 if(track->Charge()==0) {
   printf("this is not an electron! Check HFEpid method");
   return;
 }
 for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
   if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]){
     fHistHPDcaXYRes[iPart][iPtBin]->Fill(residual[0]*1.0e4);  
     fHistHPDcaZRes[iPart][iPtBin]->Fill(residual[1]*1.0e4);   
     fHistHPDcaXYPull[iPart][iPtBin]->Fill(pull[0]);  
     fHistHPDcaZPull[iPart][iPtBin]->Fill(pull[1]); 
     fHistHPDcaXY[iPart][iPtBin]->Fill(dz[0]*1.0e4); 
     fHistHPDcaZ[iPart][iPtBin]->Fill(dz[1]*1.0e4);  

   } // pt range

   else
     continue;
 }  // pt loop

 fHistHfePid[iPart][0]->Fill(esdpt);
 fHistHfePid[iPart][1]->Fill(mcpt);

}


//_______________________________________________________________________________________________
void AliHFEdca::FillHistogramsHfeDataDca(AliESDEvent * const esdEvent, AliESDtrack * const track, AliESDVertex * const vtxESDSkip)
{
// filling historgams track by track
// obtaining reconstructed dca --------------------------------------------------------------

 Float_t esdpx = track->Px();
 Float_t esdpy = track->Py();
 Float_t esdpt = TMath::Sqrt(esdpx*esdpx+esdpy*esdpy);  
 Int_t charge = track->Charge();

// obtaining errors of dca ------------------------------------------------------------------
 const AliESDVertex *primVtx = esdEvent->GetPrimaryVertex(); // UNUSED!     
 Double_t primV[3];
 primV[0] = primVtx->GetXv();
 primV[1] = primVtx->GetYv();
 primV[2] = primVtx->GetZv();

 Float_t magneticField = 0;  // initialized as 5kG
 magneticField = esdEvent->GetMagneticField();  // in kG
 Double_t beampiperadius=3.; 

 Double_t dz[2];   // error of dca in cm
 Double_t covardz[3];

 if(!track->PropagateToDCA(vtxESDSkip,magneticField, beampiperadius, dz, covardz)) return; // protection

 Double_t pull[2] = {0, 0};
 Double_t error[2] ={TMath::Sqrt(covardz[0]), TMath::Sqrt(covardz[2])};
 for(Int_t i=0; i<2; i++){
   if(error[i]!=0) pull[i] = dz[i]/error[i];   // unitless    
 }

 Int_t iPart = -1;
 if(charge<0) iPart = 0;  // electron
 if(charge>0) iPart = 1;  // positron
 if(charge==0) {
   printf("this is not an electron! Check HFEpid method\n");
   return;
 }

 for(Int_t iPtBin=0; iPtBin<kNPtBins; iPtBin++){
   if(esdpt>fgkPtIntv[iPtBin] && esdpt<=fgkPtIntv[iPtBin+1]) {
     fHistHPDataDcaXY[iPart][iPtBin]->Fill(dz[0]*1e4);
     fHistHPDataDcaZ[iPart][iPtBin]->Fill(dz[1]*1e4);      
     fHistHPDataDcaXYPull[iPart][iPtBin]->Fill(pull[0]);
     fHistHPDataDcaZPull[iPart][iPtBin]->Fill(pull[1]);		      

   }
   else continue;
 }

 fHistDataHfePid[iPart]->Fill(esdpt);

}

