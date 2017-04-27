
/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

//
//   Base class to perform D-h correlations offline (starting from D-meson and hadron TTrees)
//
//-----------------------------------------------------------------------
//  Author F.Colamaria
//  INFN Bari
//  fabio.colamaria@cern.ch
//-----------------------------------------------------------------------

#include "AliHFOfflineCorrelator.h"

//___________________________________________________________________________________________
AliHFCorrelationBranchD::AliHFCorrelationBranchD():
// default constructor
phi_D(0),
eta_D(0),
pT_D(0),
mult_D(0),
zVtx_D(0),
cent_D(0),
invMass_D(0),
period_D(0),
orbit_D(0),
BC_D(0),
IDtrig_D(0),
sel_D(0)
{

}

//___________________________________________________________________________________________
AliHFCorrelationBranchTr::AliHFCorrelationBranchTr():
// default constructor
phi_Tr(0),
eta_Tr(0),
pT_Tr(0),
mult_Tr(0),
zVtx_Tr(0),
cent_Tr(0),
period_Tr(0),
orbit_Tr(0),
BC_Tr(0),
IDtrig_Tr(0),
IDtrig2_Tr(0),
IDtrig3_Tr(0),
IDtrig4_Tr(0),
sel_Tr(0)
{

}

//___________________________________________________________________________________________
AliHFOfflineCorrelator::AliHFOfflineCorrelator():
// default constructor
fFileList(0),
fNinputFiles(0),
fFile(0x0),
fTreeD(0x0),
fTreeTr(0x0),
fOutputDistr(0x0),
fNBinsPt(0),
fnPools(0),
fnMultPools(0),
fnzVtxPools(0),
fFirstBinNum(0),
fMaxTracks(-1),
fMinCent(0.),
fMaxCent(0.),
fNumSelD(-1),
fNumSelTr(-1),
fDebug(0),
fMinD(-1),
fMaxD(-1),
fPtBinsDLow(0),
fPtBinsDUp(0),
fPtBinsTrLow(0),
fPtBinsTrUp(0),
fMassSignL(0),
fMassSignR(0),
fMassSB1L(0),
fMassSB1R(0),
fMassSB2L(0),
fMassSB2R(0),
fMultBins(0),
fzVtxBins(0),
fPrdWeights(0),
fMapEffTr(0x0),
fMapEffD(0x0),
fDmesonSpecies(kD0toKpi),
fAnType(kSE),
fDmesonLabel("Dzero"),
fDirName(0),
fNameTreeTr(0),
fNameTreeD(0),
fNameMapTr(0),
fNameMapD(0),
fOutputFileName(0),
fNameCutObj(0),
fUseEff(0),
fMake2DPlots(kFALSE),
fWeightPeriods(kTRUE),
fRemoveSoftPiInME(kTRUE)
{

}

//___________________________________________________________________________________________
AliHFOfflineCorrelator::AliHFOfflineCorrelator(const AliHFOfflineCorrelator &source):
// copy constructor
TObject(source),
fFileList(source.fFileList),
fNinputFiles(source.fNinputFiles),
fFile(source.fFile),
fTreeD(source.fTreeD),
fTreeTr(source.fTreeTr),
fOutputDistr(source.fOutputDistr),
fNBinsPt(source.fNBinsPt),
fnPools(source.fnPools),
fnMultPools(source.fnMultPools),
fnzVtxPools(source.fnzVtxPools),
fFirstBinNum(source.fFirstBinNum),
fMaxTracks(source.fMaxTracks),
fMinCent(source.fMinCent),
fMaxCent(source.fMaxCent),
fNumSelD(source.fNumSelD),
fNumSelTr(source.fNumSelTr),
fDebug(source.fDebug),
fMinD(source.fMinD),
fMaxD(source.fMaxD),
fPtBinsDLow(source.fPtBinsDLow),
fPtBinsDUp(source.fPtBinsDUp),
fPtBinsTrLow(source.fPtBinsTrLow),
fPtBinsTrUp(source.fPtBinsTrUp),
fMassSignL(source.fMassSignL),
fMassSignR(source.fMassSignR),
fMassSB1L(source.fMassSB1L),
fMassSB1R(source.fMassSB1R),
fMassSB2L(source.fMassSB2L),
fMassSB2R(source.fMassSB2R),
fMultBins(source.fMultBins),
fzVtxBins(source.fzVtxBins),
fPrdWeights(source.fPrdWeights),
fMapEffTr(source.fMapEffTr),
fMapEffD(source.fMapEffD),
fDmesonSpecies(source.fDmesonSpecies),
fAnType(source.fAnType),
fDmesonLabel(source.fDmesonLabel),
fDirName(source.fDirName),
fNameTreeTr(source.fNameTreeTr),
fNameTreeD(source.fNameTreeD),
fNameMapTr(source.fNameMapTr),
fNameMapD(source.fNameMapD),
fOutputFileName(source.fOutputFileName),
fNameCutObj(source.fNameCutObj),
fUseEff(source.fUseEff),
fMake2DPlots(source.fMake2DPlots),
fWeightPeriods(source.fWeightPeriods),
fRemoveSoftPiInME(source.fRemoveSoftPiInME)
{

}

//___________________________________________________________________________________________
AliHFOfflineCorrelator& AliHFOfflineCorrelator::operator=(const AliHFOfflineCorrelator& orig)
{
// Assignment
if (&orig == this) return *this; //if address is the same (same object), returns itself

TObject::operator=(orig); //Uses the TObject operator to assign the inherited part of the class
fFileList = orig.fFileList;
fNinputFiles = orig.fNinputFiles;
fFile = orig.fFile;
fTreeD = orig.fTreeD;
fTreeTr = orig.fTreeTr;
fOutputDistr = orig.fOutputDistr;
fNBinsPt = orig.fNBinsPt;
fnPools = orig.fnPools;
fnMultPools = orig.fnMultPools;
fnzVtxPools = orig.fnzVtxPools;
fFirstBinNum = orig.fFirstBinNum;
fMaxTracks = orig.fMaxTracks;
fMinCent = orig.fMinCent;
fMaxCent = orig.fMaxCent;
fNumSelD = orig.fNumSelD;
fNumSelTr = orig.fNumSelTr;
fDebug = orig.fDebug;
fMinD = orig.fMinD;
fMaxD = orig.fMaxD;
fPtBinsDLow = orig.fPtBinsDLow;
fPtBinsDUp = orig.fPtBinsDUp;
fPtBinsTrLow = orig.fPtBinsTrLow;
fPtBinsTrUp = orig.fPtBinsTrUp;
fMassSignL = orig.fMassSignL;
fMassSignR = orig.fMassSignR;
fMassSB1L = orig.fMassSB1L;
fMassSB1R = orig.fMassSB1R;
fMassSB2L = orig.fMassSB2L;
fMassSB2R = orig.fMassSB2R;
fMultBins = orig.fMultBins;
fzVtxBins = orig.fzVtxBins;
fPrdWeights = orig.fPrdWeights;
fMapEffTr = orig.fMapEffTr;
fMapEffD = orig.fMapEffD;
fDmesonSpecies = orig.fDmesonSpecies;
fAnType = orig.fAnType;
fDmesonLabel = orig.fDmesonLabel;
fDirName = orig.fDirName;
fNameTreeTr = orig.fNameTreeTr;
fNameTreeD = orig.fNameTreeD;
fNameMapTr = orig.fNameMapTr;
fNameMapD = orig.fNameMapD;
fOutputFileName = orig.fOutputFileName;
fNameCutObj = orig.fNameCutObj;
fUseEff = orig.fUseEff;
fMake2DPlots = orig.fMake2DPlots;
fWeightPeriods = orig.fWeightPeriods;
fRemoveSoftPiInME = orig.fRemoveSoftPiInME;

return *this; //returns pointer of the class
}

//___________________________________________________________________________________________
AliHFOfflineCorrelator::~AliHFOfflineCorrelator() {
//destructor

}

//___________________________________________________________________________________________
Bool_t AliHFOfflineCorrelator::SetDmesonSpecie(DMesonSpecies k) {

  if(k<0 || k>2) {
    printf("Error! D meson specie not correctly set!\n");
    return kFALSE;
  } else if(k==0) fDmesonLabel="Dzero";
  else if(k==1) fDmesonLabel="Dplus";  
  else fDmesonLabel="Dstar";

  fDmesonSpecies=k;
  return kTRUE;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::SetDPtBins(Int_t nBins, Double_t* ptDarray) {

  fNBinsPt=nBins;

  for(Int_t i=0;i<nBins;i++) {
    fPtBinsDLow.push_back(ptDarray[i]);
    fPtBinsDUp.push_back(ptDarray[i+1]);
  }    

  return;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::SetPoolBins(Int_t nMultPools, Double_t* multBins, Int_t nzVtxPools, Double_t* zVtxBins) {

  fnPools = nMultPools*nzVtxPools;
  fnMultPools = nMultPools;
  fnzVtxPools = nzVtxPools;

  for(Int_t i=0;i<=nMultPools;i++) {
    fMultBins.push_back(multBins[i]);
  }    
  for(Int_t i=0;i<=nzVtxPools;i++) {
    fzVtxBins.push_back(zVtxBins[i]);
  }
  
  return;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::DefineOutputObjects(){

  fOutputDistr = new TList();

  Int_t massBins;
  Double_t massLow, massUp;
  TString namePlot;

  if(fDmesonSpecies==kD0toKpi) { //D0
    massBins = 150;
    massLow = 1.5848;
    massUp = 2.1848;
  } else if(fDmesonSpecies==kDplusKpipi) { //D+
    massBins = 200;
    massLow = 1.665;
    massUp = 2.065;  
  } else { //D*
    massBins = 400;
    massLow = 0.13;
    massUp = 0.19;  
  }

  //Mass plots for signal yield extraction and sideband normalization
  fOutputMass = new TList();
  fOutputMass->SetOwner();
  fOutputMass->SetName("MassPlots");
  
  for(Int_t iBin=0; iBin<fNBinsPt; iBin++) {
    //Defines 3D (DPhi,DEta,Minv) plots
    namePlot = Form("histMass_%d",fFirstBinNum+iBin);
    TH1F* h1D = new TH1F(namePlot.Data(),"Mass plots for triggers",150,1.5648,2.1648);
    h1D->Sumw2();
    fOutputMass->Add(h1D);
    
    if(fUseEff) {
	  namePlot = Form("histMass_WeigD0Eff_%d",fFirstBinNum+iBin);
      TH1F* h1Dw = new TH1F(namePlot.Data(),"Mass plots for triggers",150,1.5648,2.1648);
      h1Dw->Sumw2();
      fOutputMass->Add(h1Dw);
    }
  }

  for(Int_t iBin=0; iBin<fNBinsPt; iBin++) {
    for(Int_t iAssPt=0; iAssPt<(int)fPtBinsTrLow.size(); iAssPt++) {
      for(Int_t iPool=0; iPool<fnPools; iPool++) {
       //Defines 3D (DPhi,DEta,Minv) plots
        namePlot = Form("h3DCorrelations_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
        TH3F* h3D = new TH3F(namePlot.Data(),"3D Correlation distribution",32,-TMath::Pi()/2.,3.*TMath::Pi()/2.,16,-1.6,1.6,massBins,massLow,massUp);
        h3D->Sumw2();
        fOutputDistr->Add(h3D);
      
        if(fDebug) { //Defines Eta distribution for all D triggers and tracks in a given pool
          namePlot = Form("hEtaD_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
          TH1F* hEtaD = new TH1F(namePlot.Data(),"D-meson eta distribution",32,-0.8,0.8);
          fOutputDistr->Add(hEtaD);

          namePlot = Form("hEtaTr_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
          TH1F* hEtaTr = new TH1F(namePlot.Data(),"Track eta distribution",32,-0.8,0.8);
          fOutputDistr->Add(hEtaTr);
        }
      }
    }
  }
      
  if(fMake2DPlots) { 
    for(Int_t iBin=0; iBin<fNBinsPt; iBin++) {
      for(Int_t iAssPt=0; iAssPt<(int)fPtBinsTrLow.size(); iAssPt++) {
        for(Int_t iPool=0; iPool<fnPools; iPool++) {
          //Defines 2D (DPhi,DEta,Minv) plots
          namePlot = Form("h2DCorrelations_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
          TH2F* h2D_Sign = new TH2F(namePlot.Data(),"2D Correlation distribution - Signal Region",32,-TMath::Pi()/2.,3.*TMath::Pi()/2.,16,-1.6,1.6);
          h2D_Sign->Sumw2();
          fOutputDistr->Add(h2D_Sign);

          namePlot = Form("h2DCorrelations_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
          TH2F* h2D_SB = new TH2F(namePlot.Data(),"2D Correlation distribution - Sidebands",32,-TMath::Pi()/2.,3.*TMath::Pi()/2.,16,-1.6,1.6);
          h2D_SB->Sumw2();
          fOutputDistr->Add(h2D_SB);

          if(fDebug) { //Defines Eta distribution for all D triggers and tracks in a given pool (only when D is in signal region/SB)
            namePlot = Form("hEtaD_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
            TH1F* hEtaD_Sign = new TH1F(namePlot.Data(),"D-meson eta distribution (signal region)",32,-0.8,0.8);
            fOutputDistr->Add(hEtaD_Sign);

            namePlot = Form("hEtaTr_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
            TH1F* hEtaTr_Sign = new TH1F(namePlot.Data(),"Track eta distribution (signal region)",32,-0.8,0.8);
            fOutputDistr->Add(hEtaTr_Sign);

            namePlot = Form("hEtaD_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
            TH1F* hEtaD_SB = new TH1F(namePlot.Data(),"D-meson eta distribution (sidebands)",32,-0.8,0.8);
            fOutputDistr->Add(hEtaD_SB);

            namePlot = Form("hEtaTr_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
            TH1F* hEtaTr_SB = new TH1F(namePlot.Data(),"Track eta distribution (sidebands)",32,-0.8,0.8);
            fOutputDistr->Add(hEtaTr_SB);
          }
        }
      }
    }
  }
  
  return;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::MakeDetaDphiPlots(Double_t* MsigL, Double_t* MsigR, Double_t* MSB1L, Double_t* MSB1R, Double_t* MSB2L, Double_t* MSB2R) {

  fMake2DPlots=kTRUE;

  //Sets mass limits
  for(Int_t i=0;i<fNBinsPt;i++) {
    fMassSignL.push_back(MsigL[i]);
    fMassSignR.push_back(MsigR[i]);
    fMassSB1L.push_back(MSB1L[i]);
    fMassSB1R.push_back(MSB1R[i]);
    if(fDmesonSpecies!=kDStarD0pi) {
      fMassSB2L.push_back(MSB2L[i]);
      fMassSB2R.push_back(MSB2R[i]);
    }
  }

  return;
}

//___________________________________________________________________________________________
Bool_t AliHFOfflineCorrelator::Correlate() {

  PrintCfg(); 

  //Define weights for periods, in ME analysis
  if(fWeightPeriods && fAnType==kME) {
    if(!DefinePeriodWeights()) return kFALSE;
  }

  //Creates the output plots looking at the configuration of pT bins
  DefineOutputObjects();

  //Safety checks
  if(fnPools<1) {
    std::cout << "ERROR! Pools are not correctly set! Exiting..." << std::endl;
    return kFALSE;
  }
  if(fMake2DPlots) {
    if(fMassSignL.size()<1 || fMassSignR.size()<1 || fMassSB1L.size()<1 || fMassSB1R.size()<1) {
      std::cout << "ERROR! You are requiring (dPhi,dEta) plots, but mass lmits of signal/SB regions are not correctly set! Exiting..." << std::endl;
      return kFALSE;
    }
    if(fDmesonSpecies!=kDStarD0pi && (fMassSB2L.size()<1 || fMassSB2R.size()<1)) {
      std::cout << "ERROR! You are requiring (dPhi,dEta) plots, but mass lmits of signal/SB regions are not correctly set! Exiting..." << std::endl;
      return kFALSE;
    }
  }

  for(Int_t iFile=0; iFile<(int)fFileList.size(); iFile++) {
    Bool_t success = CorrelateSingleFile(iFile);
    if(!success) {
      std::cout << "Error in the evaluation of correlations for file #" << iFile << ". Exiting..." << std::endl;
      return kFALSE;
    }
  }

  //Normalize 2D plots, if they are produced
  if(fMake2DPlots) NormalizeMEPlots();

  //Save in the output file the distributions
  SaveOutputPlots();

  return kTRUE; //if everything was fine, returns kTRUE!
}

//___________________________________________________________________________________________
Bool_t AliHFOfflineCorrelator::CorrelateSingleFile(Int_t iFile) {

  std::cout << "Opening file: " << fFileList.at(iFile) << std::endl;

  fFile = TFile::Open((TString)(fFileList.at(iFile)).Data());
  if(!fFile){
    std::cout << "File " << fFileList.at(iFile) << " cannot be opened! check your file path!" << std::endl;
    return kFALSE;
  }

  TDirectoryFile *dir = (TDirectoryFile*)fFile->Get(fDirName.Data());
  if(!dir){
    std::cout << "Directory " << fDirName << " is missing! Check its spelling/the file content" << std::endl;
    fFile->ls();
    return kFALSE;
  }  

  fTreeD = (TTree*)dir->Get(fNameTreeD.Data());
  fTreeTr = (TTree*)dir->Get(fNameTreeTr.Data());
  if(!fTreeD || !fTreeTr){
    std::cout << "TTrees not found! Check its spelling/the directory content" << std::endl;
    dir->ls();
    return kFALSE;
  }  
  
  if(fUseEff) {
    AliHFAssociatedTrackCuts *cutObj = (AliHFAssociatedTrackCuts*)dir->Get(fNameCutObj.Data());
    if(!cutObj){
      std::cout << "Wrong cut file name, or missing cut file! (you chose: " << fNameCutObj << ")" << std::endl;
      fFile->ls();
      return kFALSE;
    } 
    fMapEffD = (TH2F*)cutObj->GetTrigEfficiencyWeight();
    fMapEffTr = (TH3F*)cutObj->GetEfficiencyWeight();
    if(!fMapEffD || !fMapEffTr){
      std::cout << "Efficiency maps missing! Check the spelling (you chose: " << fNameMapD << "/" << fNameMapTr << ") or the file content content" << std::endl;
      fFile->ls();
      return kFALSE;
    }  
  }

  AliHFCorrelationBranchD *brD = 0;
  AliHFCorrelationBranchTr *brTr = 0;

  fTreeD->SetBranchAddress("branchD",&brD);
  fTreeTr->SetBranchAddress("branchTr",&brTr);

  std::cout << "File contains a total of " << fTreeD->GetEntries() << " D mesons and of " << fTreeTr->GetEntries() << " associated tracks" << std::endl;
  std::cout << "Correlating..." << std::endl;

  TString namePlot = "";
  Int_t poolD = 0, poolTr = 0;
  Int_t minDLoop = 0, maxDLoop = fTreeD->GetEntries();
  Int_t minTrackLoop = 0, maxTrackLoop = fTreeTr->GetEntries();

  if(fMinD>=0) minDLoop=fMinD;
  if(fMaxD>=0) maxDLoop=fMaxD;
  if(fMinD>fMaxD) {printf("Warning! Wrong settings of D-meson loop edges! Exiting...\n"); return kFALSE;}
  if(fMinD>fTreeD->GetEntries()) {printf("Warning! The lower edge of D meson loop exceeds the number of D in the TTree! No loop will be done\n"); return kTRUE;}
  if(fMaxD>fTreeD->GetEntries()) {printf("Warning! The upper edge of D meson loop exceeds the number of D in the TTree!\n"); maxDLoop = fTreeD->GetEntries();}

  TRandom3 *tRnd = new TRandom3();
  tRnd->SetSeed(1);

  TStopwatch *tim = new TStopwatch();
  tim->Start();

  for(Int_t iD=minDLoop; iD<maxDLoop; iD++) {  //loop on D-mesons in tree   

    //time monitoring
    if(iD%10==0) {
      tim->Stop();
      std::cout << "--- D-meson " << iD << std::endl;
      tim->Print();
      tim->Continue();
    }

    fTreeD->GetEntry(iD); 
    Int_t ptBinD = PtBin(brD->pT_D);
    if(ptBinD<0) continue;  
    if(fNumSelD>=0 && (brD->sel_D>>fNumSelD)%2!=1) continue; //important in case of multiple selection (default selection is 0)
    if(fMinCent!=0 && fMaxCent!=0) {if(brD->cent_D < fMinCent || brD->cent_D > fMaxCent) continue;} //skip triggers outside centrality range
    
    poolD = GetPoolBin(brD->mult_D,brD->zVtx_D); 

    if(fMaxTracks>0) { //select random range of 'fMaxTracks' tracks in the TTree of tracks (the range changes for each D meson to use all the sample)
      if(fMaxTracks>=fTreeTr->GetEntries()) printf("Warning! Requested to loop on more tracks than the available number! Standard loop being done\n");
      else {
        minTrackLoop = tRnd->Rndm()*(fTreeTr->GetEntries()-fMaxTracks);
        maxTrackLoop = fMaxTracks+minTrackLoop;
      }
    }

    Int_t fillOnce[(int)fPtBinsTrLow.size()]; for(int ii=0;ii<(int)fPtBinsTrLow.size();ii++) fillOnce[ii]=0;

    //Fill mass plots
    ((TH1F*)(fOutputMass->FindObject(Form("histMass_%d",fFirstBinNum+ptBinD))))->Fill(brD->invMass_D);
    if(fUseEff) ((TH1F*)(fOutputMass->FindObject(Form("histMass_WeigD0Eff_%d",fFirstBinNum+ptBinD))))->Fill(brD->invMass_D,GetEfficiencyWeightDOnly(brD));
    
    //Correlation plots!
    for(Int_t iTr=minTrackLoop; iTr<maxTrackLoop; iTr++) {  //loop on associated tracks in tree

      fTreeTr->GetEntry(iTr); 
      if(fAnType==kSE && (brD->period_D!=brTr->period_Tr || brD->orbit_D!=brTr->orbit_Tr || brD->BC_D!=brTr->BC_Tr)) continue; //skips D and tracks from different events in ME 
      if(fAnType==kME && (brD->period_D==brTr->period_Tr && brD->orbit_D==brTr->orbit_Tr && brD->BC_D==brTr->BC_Tr)) continue; //skips D and tracks from same event in SE 

      if(fAnType==kSE && brD->IDtrig_D==brTr->IDtrig_Tr) continue; //skips D0 daughter association with their own trigger (or own soft-pion, for the D0)
      if(fAnType==kSE && brD->IDtrig_D==brTr->IDtrig2_Tr) continue; //skips D0 daughter association with their own trigger (or own soft-pion, for the D0)
      if(fAnType==kSE && brD->IDtrig_D==brTr->IDtrig3_Tr) continue; //skips D0 daughter association with their own trigger (or own soft-pion, for the D0)
      if(fAnType==kSE && brD->IDtrig_D==brTr->IDtrig4_Tr) continue; //skips D0 daughter association with their own trigger (or own soft-pion, for the D0)

      if(fNumSelTr>=0 && (brTr->sel_Tr>>fNumSelTr)%2!=1) continue; //important in case of multiple selection (default selection is 0)
	  if(fMinCent!=0 && fMaxCent!=0) {if(brTr->cent_Tr < fMinCent || brTr->cent_Tr > fMaxCent) continue;} //skip tracks outside centrality range

      poolTr = GetPoolBin(brTr->mult_Tr,brTr->zVtx_Tr);
      if(poolD<0 || poolTr<0 || poolD!=poolTr) continue;  //skips if pools of D and tracks do not match, or if pool number is wrong

      Double_t weight = 1.;
      if(fUseEff) weight = GetEfficiencyWeight(brD,brTr); //efficiency weighting
      if(fWeightPeriods && fAnType==kME) weight*=fPrdWeights.at(iFile); //period-by-period weighting
      Double_t deltaPhi, deltaEta;
      GetCorrelationsValue(brD,brTr,deltaPhi,deltaEta);

      if(fRemoveSoftPiInME && fDmesonSpecies==kD0toKpi && fAnType==kME && deltaPhi > -0.2 && deltaPhi < 0.2 && deltaEta > -0.2 && deltaEta < 0.2) { //ME fake soft pi cut
		Bool_t reject = IsSoftPionFromDstar(brD,brTr);
		if(reject) continue;
	  }

      for(Int_t iRng=0; iRng<(int)fPtBinsTrLow.size(); iRng++) {  //loop on associated track ranges

        //fill 3D and 2D correlation plots
        if(brTr->pT_Tr < fPtBinsTrLow.at(iRng) || brTr->pT_Tr > fPtBinsTrUp.at(iRng)) continue; //skip cases where associated track pT is out of range
        namePlot = Form("h3DCorrelations_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
        ((TH3F*)(fOutputDistr->FindObject(namePlot)))->Fill(deltaPhi,deltaEta,brD->invMass_D,weight);

        if(fMake2DPlots) {
	      if(brD->invMass_D > fMassSignL.at(ptBinD) && brD->invMass_D < fMassSignR.at(ptBinD)) {
            namePlot = Form("h2DCorrelations_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
            ((TH2F*)(fOutputDistr->FindObject(namePlot)))->Fill(deltaPhi,deltaEta,weight);
          }     
          if(brD->invMass_D > fMassSB1L.at(ptBinD) && brD->invMass_D < fMassSB1R.at(ptBinD)) {
            namePlot = Form("h2DCorrelations_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
            ((TH2F*)(fOutputDistr->FindObject(namePlot)))->Fill(deltaPhi,deltaEta,weight);
          }
          if(fDmesonSpecies!=kDStarD0pi && (brD->invMass_D > fMassSB2L.at(ptBinD) && brD->invMass_D < fMassSB2R.at(ptBinD))) {
            namePlot = Form("h2DCorrelations_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
            ((TH2F*)(fOutputDistr->FindObject(namePlot)))->Fill(deltaPhi,deltaEta,weight);
          }
        } //end if 2D plots

        //***fill debug plots***
        if(fDebug) {
          if(brTr->pT_Tr < fPtBinsTrLow.at(iRng) || brTr->pT_Tr > fPtBinsTrUp.at(iRng)) continue; //skip cases where associated track pT is out of range
          namePlot = Form("hEtaD_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
          if(fillOnce[iRng]==0) ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brD->eta_D);  //in the track loop, fill only once for D-meson!
          namePlot = Form("hEtaTr_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
          ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brTr->eta_Tr);  //fill at each track iteration for the tracks!
          if(fMake2DPlots) {
	        if(brD->invMass_D > fMassSignL.at(ptBinD) && brD->invMass_D < fMassSignR.at(ptBinD)) {
              namePlot = Form("hEtaD_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
              if(fillOnce[iRng]==0) ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brD->eta_D);
 	          namePlot = Form("hEtaTr_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
   	          ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brTr->eta_Tr);
            }     
            if(brD->invMass_D > fMassSB1L.at(ptBinD) && brD->invMass_D < fMassSB1R.at(ptBinD)) {
              namePlot = Form("hEtaD_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
              if(fillOnce[iRng]==0) ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brD->eta_D);
 	          namePlot = Form("hEtaTr_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
   	          ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brTr->eta_Tr);
            }
            if(fDmesonSpecies!=kDStarD0pi && (brD->invMass_D > fMassSB2L.at(ptBinD) && brD->invMass_D < fMassSB2R.at(ptBinD))) {
              namePlot = Form("hEtaD_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
              if(fillOnce[iRng]==0) ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brD->eta_D);
 	          namePlot = Form("hEtaTr_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+ptBinD,fPtBinsTrLow.at(iRng),fPtBinsTrUp.at(iRng),poolD);
   	          ((TH1F*)(fOutputDistr->FindObject(namePlot)))->Fill(brTr->eta_Tr);
            }
          } //end if 2D plots (for debug plots)
          fillOnce[iRng]++; //to avoid re-filling of D-meson debug plots with further tracks for the same meson
        } //***end fill debug plots***

      } //end ass track ranges
    } //end ass track loop
  } //end D-meson loop

  std::cout << "Done! Closing file." << std::endl;

  fFile->TFile::Close();

  return kTRUE;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::GetCorrelationsValue(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr, Double_t &deltaPhi, Double_t &deltaEta) {

  deltaPhi = brD->phi_D - brTr->phi_Tr;
  if(deltaPhi < -TMath::Pi()/2.)   deltaPhi = deltaPhi + 2*TMath::Pi();
  if(deltaPhi > 3.*TMath::Pi()/2.) deltaPhi = deltaPhi - 2*TMath::Pi();

  deltaEta =  brD->eta_D - brTr->eta_Tr;

  return;
}

//___________________________________________________________________________________________
Double_t AliHFOfflineCorrelator::GetEfficiencyWeight(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr) {

  Double_t effD = 1, effTr = 1;
   
  Int_t binD=fMapEffD->FindBin(brD->pT_D,brD->mult_D);
  if(fMapEffD->IsBinUnderflow(binD)||fMapEffD->IsBinOverflow(binD))return 1.;
  effD = fMapEffD->GetBinContent(binD);

  Int_t binTr=fMapEffTr->FindBin(brTr->pT_Tr,brTr->eta_Tr,brTr->zVtx_Tr);
  if(fMapEffTr->IsBinUnderflow(binTr)||fMapEffTr->IsBinOverflow(binTr))return 1.;
  effTr = fMapEffTr->GetBinContent(binTr);

  return 1./(effD*effTr);
}

//___________________________________________________________________________________________
Double_t AliHFOfflineCorrelator::GetEfficiencyWeightDOnly(AliHFCorrelationBranchD *brD) {

  Double_t effD = 1;
   
  Int_t binD=fMapEffD->FindBin(brD->pT_D,brD->mult_D);
  if(fMapEffD->IsBinUnderflow(binD)||fMapEffD->IsBinOverflow(binD))return 1.;
  effD = fMapEffD->GetBinContent(binD);

  return 1./(effD);
}

//___________________________________________________________________________________________
Int_t AliHFOfflineCorrelator::PtBin(Double_t pt) const {

  Int_t ptbin=-1;
  if(pt<fPtBinsDLow.at(0)) return ptbin;
  for (Int_t i=0;i<fNBinsPt;i++) {
    if(pt<fPtBinsDUp.at(i)) {
      ptbin=i;
      break;
    }
  }
  if(pt>fPtBinsDUp.at(fNBinsPt-1)) ptbin=-1;  
  return ptbin;
}

//--------------------------------------------------------------------------
Int_t AliHFOfflineCorrelator::GetPoolBin(Double_t mult, Double_t zVtx) const {

    Int_t poolbin = -1;
    Int_t centbin = -1;
    Int_t zvtxbin = -1;
        
    if(mult < fMultBins.at(0)) return poolbin;
    if(zVtx < fzVtxBins.at(0)) return poolbin;
  

    for (Int_t i=0;i<fnMultPools;i++){
        if(mult<fMultBins.at(i+1)) {
            centbin=i;
            break;
        }
    }

    for (Int_t i=0;i<fnzVtxPools;i++){
        if(zVtx<fzVtxBins.at(i+1)) {
            zvtxbin=i;
            break;
        }
    }

    if(mult>fMultBins.at(fnMultPools) || zVtx>fzVtxBins.at(fnzVtxPools)) return -1;

    poolbin = centbin  + zvtxbin*fnMultPools;

    return poolbin;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::NormalizeMEPlots() {

  TString namePlot, nameOrigPlot;
  Double_t factorAdd = 0;
  Double_t bin0phi, bin0eta;

  for(Int_t iBin=0; iBin<fNBinsPt; iBin++) {
    for(Int_t iAssPt=0; iAssPt<(int)fPtBinsTrLow.size(); iAssPt++) {
      for(Int_t iPool=0; iPool<fnPools; iPool++) {
      
        //Defines 3D (DPhi,DEta,Minv) plots
        nameOrigPlot = Form("h2DCorrelations_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
        namePlot = Form("Norm_h2DCorrelations_Sign_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
        TH2F* h2D_Sign_Norm = (TH2F*)(fOutputDistr->FindObject(nameOrigPlot))->Clone(namePlot);

        if(iBin==0 && iAssPt==0) {
          bin0phi = h2D_Sign_Norm->GetXaxis()->FindBin(0.);
          bin0eta = h2D_Sign_Norm->GetYaxis()->FindBin(0.);
        }

        factorAdd = 0;
        for(int in=-1;in<=0;in++) factorAdd+=h2D_Sign_Norm->GetBinContent(bin0phi+in,bin0eta);
        for(int in=-1;in<=0;in++) factorAdd+=h2D_Sign_Norm->GetBinContent(bin0phi+in,bin0eta-1);
        factorAdd/=4.;
        if(factorAdd) h2D_Sign_Norm->Scale(1./factorAdd);

        fOutputDistr->Add(h2D_Sign_Norm);

        nameOrigPlot = Form("h2DCorrelations_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
        namePlot = Form("Norm_h2DCorrelations_SB_Bin%d_%1.1fto%1.1f_p%d",fFirstBinNum+iBin,fPtBinsTrLow.at(iAssPt),fPtBinsTrUp.at(iAssPt),iPool);
        TH2F* h2D_SB_Norm = (TH2F*)(fOutputDistr->FindObject(nameOrigPlot))->Clone(namePlot);

        factorAdd = 0;
        for(int in=-1;in<=0;in++) factorAdd+=h2D_SB_Norm->GetBinContent(bin0phi+in,bin0eta);
        for(int in=-1;in<=0;in++) factorAdd+=h2D_SB_Norm->GetBinContent(bin0phi+in,bin0eta-1);
        factorAdd/=4.;
        if(factorAdd) h2D_SB_Norm->Scale(1./factorAdd);

        fOutputDistr->Add(h2D_SB_Norm);
      }
    }
  }

}

//___________________________________________________________________________________________
Bool_t AliHFOfflineCorrelator::DefinePeriodWeights() {

  Double_t maxTracks = 0;
  std::vector<Double_t> periodTracks;

  for(Int_t iFile=0; iFile<(int)fFileList.size(); iFile++) {
    std::cout << "Defining weights for the different periods" << std::endl;

    fFile = TFile::Open((TString)(fFileList.at(iFile)).Data());
    if(!fFile){
      std::cout << "File " << fFileList.at(iFile) << " cannot be opened! check your file path!" << std::endl;
      return kFALSE;
    }
    
    TDirectoryFile *dir = (TDirectoryFile*)fFile->Get(fDirName.Data());
    if(!dir){
      std::cout << "Directory " << fDirName << " is missing! Check its spelling/the file content" << std::endl;
      fFile->ls();
      return kFALSE;
    }  

    fTreeTr = (TTree*)dir->Get(fNameTreeTr.Data());
    if(!fTreeTr){
      std::cout << "TTrees not found! Check its spelling/the directory content" << std::endl;
      dir->ls();
      return kFALSE;
    }  

    Double_t treeEntries = fTreeTr->GetEntries();
    if(fMaxTracks>0) { //max number of tracks used
      if(maxTracks<TMath::Min((Double_t)treeEntries,(Double_t)fMaxTracks)) maxTracks = (TMath::Min((Double_t)treeEntries,(Double_t)fMaxTracks));
      periodTracks.push_back(TMath::Min((Double_t)treeEntries,(Double_t)fMaxTracks));
    } else { //no limit in track usage
      if(maxTracks<treeEntries) maxTracks = treeEntries;
      periodTracks.push_back(treeEntries);
    }
  }

  for(Int_t iFile=0; iFile<(int)fFileList.size(); iFile++) {
    fPrdWeights.push_back((Double_t)maxTracks/(Double_t)periodTracks.at(iFile));
    std::cout << "--> File #" << iFile << " weight: " << (Double_t)fPrdWeights.at(iFile) << std::endl;
  }

  return kTRUE;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::SaveOutputPlots() {

  TFile *f = new TFile(fOutputFileName.Data(),"RECREATE");
  fOutputDistr->Write();
  fOutputMass->Write();
  f->Close();
  return;
}

//___________________________________________________________________________________________
Bool_t AliHFOfflineCorrelator::IsSoftPionFromDstar(AliHFCorrelationBranchD *brD, AliHFCorrelationBranchTr *brTr) {
	//
	// Calculates invmass of track+D0 and rejects if compatible with D*
	// (to remove fake pions from D* in ME events - in SE the cut is applied in the main task)
	// 
	Double_t nsigma = 3.;
	
	Double_t mPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
	Double_t mK = TDatabasePDG::Instance()->GetParticle(321)->Mass();
	Double_t pxD = brD->pT_D*TMath::Cos(brD->phi_D);
	Double_t pyD = brD->pT_D*TMath::Sin(brD->phi_D);
	Double_t pzD = brD->pT_D*TMath::SinH(brD->eta_D);
	Double_t pxTr = brTr->pT_Tr*TMath::Cos(brTr->phi_Tr);
	Double_t pyTr = brTr->pT_Tr*TMath::Sin(brTr->phi_Tr);
	Double_t pzTr = brTr->pT_Tr*TMath::SinH(brTr->eta_Tr);
	Double_t invmassDstar1 = 0, invmassDstar2 = 0; 
	
	//hyp 1 (pi,K) - D0
	Double_t e1Pi = TMath::Sqrt(mPi*mPi + TMath::Power(brD->pXdaug1_D,2) + TMath::Power(brD->pYdaug1_D,2) + TMath::Power(brD->pZdaug1_D,2));
	Double_t e2K = TMath::Sqrt(mK*mK + TMath::Power(brD->pXdaug2_D,2) + TMath::Power(brD->pYdaug2_D,2) + TMath::Power(brD->pZdaug2_D,2));
	//hyp 2 (K,pi) - D0bar
	Double_t e1K = TMath::Sqrt(mK*mK + TMath::Power(brD->pXdaug1_D,2) + TMath::Power(brD->pYdaug1_D,2) + TMath::Power(brD->pZdaug1_D,2));
	Double_t e2Pi = TMath::Sqrt(mPi*mPi + TMath::Power(brD->pXdaug2_D,2) + TMath::Power(brD->pYdaug2_D,2) + TMath::Power(brD->pZdaug2_D,2));
		
	Double_t psum2 = TMath::Power(pxD+pxTr,2) + TMath::Power(pyD+pyTr,2) + TMath::Power(pzD+pzTr,2);

	switch(brD->hyp_D) { //there's no 3, since each mass hypothesis is filled separately (and 3's are fileld two times, one as 1 and one as 2)
		case 1:
			invmassDstar1 = TMath::Sqrt(TMath::Power(e1Pi+e2K+TMath::Sqrt(mPi*mPi + pxTr*pxTr + pyTr*pyTr + pzTr*pzTr),2.)-psum2);
			if ((TMath::Abs(invmassDstar1-brD->invMass_D)-0.14543) < nsigma*800.*pow(10.,-6.)) return kTRUE;
			break;
		case 2:
			invmassDstar2 = TMath::Sqrt(TMath::Power(e2Pi+e1K+TMath::Sqrt(mPi*mPi + pxTr*pxTr + pyTr*pyTr + pzTr*pzTr),2.)-psum2);
			if ((TMath::Abs(invmassDstar2-brD->invMass_D)-0.14543) < nsigma*800.*pow(10.,-6.)) return kTRUE;
			break;
	}
	
	return kFALSE;
}

//___________________________________________________________________________________________
void AliHFOfflineCorrelator::PrintCfg() const {

  std::cout << "*************** Configuration ***************\n";
  std::cout << " N# PtBins = " << fNBinsPt << "\n";
  std::cout << "--------------- PtBin limits ----------------\n";
  for (int i=0; i<fNBinsPt; i++) {
    std::cout << "Bin "<<i<<" = "<<fPtBinsDLow.at(i)<<" - "<<fPtBinsDUp.at(i)<<"\n";
  }
  std::cout << "-------------- Assoc pT ranges --------------\n";
  for (int i=0; i<(int)fPtBinsTrLow.size(); i++) {
    std::cout << "Range "<<i<<" = "<<fPtBinsTrLow.at(i)<<" - "<<fPtBinsTrUp.at(i)<<"\n";
  }
  std::cout << "------------ Pool ranges (mult.) -------------\n";
  for (int i=0; i<fnMultPools; i++) {
    std::cout << "Pool "<<i<<" = "<<fMultBins.at(i)<<" - "<<fMultBins.at(i+1)<<"\n";
  }
  std::cout << "------------ Pool ranges (z Vtx) -------------\n";
  for (int i=0; i<fnzVtxPools; i++) {
    std::cout << "Pool "<<i<<" = "<<fzVtxBins.at(i)<<" - "<<fzVtxBins.at(i+1)<<"\n";
  }
  std::cout << "----------------------------------------------\n";
  std::cout << "Total number of pools = " << fnPools << "\n";
  std::cout << "Max #tracks used for ME correlations = " << fMaxTracks << "\n";    

  if(fMake2DPlots) {
    std::cout << "---------------- Mass ranges ----------------\n";
    for (int i=0; i<fNBinsPt; i++) {
      if(fDmesonSpecies!=kDStarD0pi) std::cout << "Bin "<<i<<": Signal region = "<<fMassSignL.at(i)<<" - "<<fMassSignR.at(i)<<": SB1 = "<<fMassSB1L.at(i)<<" - "<<fMassSB1R.at(i)<<": SB2 = "<<fMassSB2L.at(i)<<" - "<<fMassSB2R.at(i)<<"\n";
      else std::cout << "Bin "<<i<<": Signal region = "<<fMassSignL.at(i)<<" - "<<fMassSignR.at(i)<<": SB1 = "<<fMassSB1L.at(i)<<" - "<<fMassSB1R.at(i)<<"\n";
    }
    std::cout << "----------------------------------------------\n";
  }

  std::cout << "----------------------------------------------\n";
  std::cout << " D-meson specie = "<<fDmesonLabel<<std::endl;
  std::cout << "----------------------------------------------\n";
  std::cout << " Analysis Type (SE=0, ME=1)= "<<fAnType<<"\n";
  std::cout << "----------------------------------------------\n";
  std::cout << " Efficiency weights = "<<fUseEff<<"\n";
  std::cout << "----------------------------------------------\n";
  std::cout << " Produce 2D plots = "<<fMake2DPlots<<"\n";
  std::cout << "----------------------------------------------\n";
  std::cout << " Weight periods (if ME) = "<<fWeightPeriods<<"\n";
  std::cout << "----------------------------------------------\n";
  std::cout << " Fake soft pi rejection in Me (D0) = "<<fRemoveSoftPiInME<<"\n";
  std::cout << "----------------------------------------------\n";
}

