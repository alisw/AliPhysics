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

#include <Riostream.h>

#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <TProfile.h>
#include <TList.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

#include "AliAnalysisTaskFragmentationFunction.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODHandler.h"
#include "AliAODTrack.h"
#include "AliJetHeader.h"
#include "AliAODEvent.h"
#include "AliAODJet.h"
#include "AliAODDiJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliAnalysisHelperJetTasks.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskFragmentationFunction)

//#######################################################################
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction():
  AliAnalysisTaskSE(),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fBranchRec(""),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseAODJetInput(kFALSE),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fFilterMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef),  
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1.),
  fExternalWeight(1.),    
  fRecEtaWindow(0.5),
  fR(0.),
  fdRdNdxi(0.7),
  fPartPtCut(0.),
  fEfactor(1.),
  fNff(5),
  fNim(5),
  fList(0x0),
  fGlobVar(1),
  // Number of energy bins
  fnEBin(6),         
  fEmin(10.),
  fEmax(70.),
  fnEInterval(6),
  // Number of radius bins
  fnRBin(10),        
  fRmin(0.1),
  fRmax(1.),	
  fnRInterval(9),
  // eta histograms
  fnEtaHBin(50),
  fEtaHBinMin(-0.9),
  fEtaHBinMax(0.9),
  // phi histograms
  fnPhiHBin(60),
  fPhiHBinMin(0.),
  fPhiHBinMax(2*TMath::Pi()),
  // pt histograms
  fnPtHBin(300),
  fPtHBinMin(0.),
  fPtHBinMax(300.),
  // E histograms
  fnEHBin(300),
  fEHBinMin(0.),
  fEHBinMax(300.),
  // Xi histograms
  fnXiHBin(27),
  fXiHBinMin(0.),
  fXiHBinMax(9.),
  // Pthad histograms
  fnPthadHBin(60),
  fPthadHBinMin(0.),
  fPthadHBinMax(30.),
  // z histograms
  fnZHBin(30), 
  fZHBinMin(0.),
  fZHBinMax(1.),
  // theta histograms
  fnThetaHBin(200),
  fThetaHBinMin(-0.5),
  fThetaHBinMax(1.5),
  fnCosThetaHBin(100),
  fcosThetaHBinMin(0.),
  fcosThetaHBinMax(1.),
  // kT histograms
  fnkTHBin(25), 
  fkTHBinMin(0.),
  fkTHBinMax(5.),
  // kT histograms
  fnRHBin(10),
  fRHBinMin(0),
  fRHBinMax(1),
  // pt trig
  fnPtTrigBin(10),
  //Histograms
  fEtaMonoJet1H(0x0),
  fPhiMonoJet1H(0x0),
  fPtMonoJet1H(0x0),
  fEMonoJet1H(0x0),
  fdNdXiMonoJet1H(0x0),
  fdNdPtMonoJet1H(0x0),
  fdNdZMonoJet1H(0x0),
  fdNdThetaMonoJet1H(0x0),
  fdNdcosThetaMonoJet1H(0x0),
  fdNdkTMonoJet1H(0x0),
  fdNdpTvsZMonoJet1H(0x0),
  fShapeMonoJet1H(0x0),
  fNMonoJet1sH(0x0),
  fThetaPtPartMonoJet1H(0x0),
  fcosThetaPtPartMonoJet1H(0x0),
  fkTPtPartMonoJet1H(0x0),
  fThetaPtJetMonoJet1H(0x0),
  fcosThetaPtJetMonoJet1H(0x0),
  fkTPtJetMonoJet1H(0x0),
  fpTPtJetMonoJet1H(0x0),
  farrayEmin(0x0),
  farrayEmax(0x0),
  farrayRadii(0x0),
  farrayPtTrigmin(0x0),
  farrayPtTrigmax(0x0),
  // Track control plots
  fptAllTracks(0x0),
  fetaAllTracks(0x0),
  fphiAllTracks(0x0),
  fetaphiptAllTracks(0x0),
  fetaphiAllTracks(0x0),
  fptAllTracksCut(0x0),
  fetaAllTracksCut(0x0),
  fphiAllTracksCut(0x0),
  fetaphiptAllTracksCut(0x0),
  fetaphiAllTracksCut(0x0),
  fptTracks(0x0),
  fetaTracks(0x0),
  fphiTracks(0x0),
  fdetaTracks(0x0),
  fdphiTracks(0x0),
  fetaphiptTracks(0x0),
  fetaphiTracks(0x0),
  fdetadphiTracks(0x0),
  fptTracksCut(0x0),
  fetaTracksCut(0x0),
  fphiTracksCut(0x0),
  fdetaTracksCut(0x0),
  fdphiTracksCut(0x0),
  fetaphiptTracksCut(0x0),
  fetaphiTracksCut(0x0),
  fdetadphiTracksCut(0x0),
  fNPtTrig(0x0),
  fNPtTrigCut(0x0),
  fvertexXY(0x0),
  fvertexZ(0x0),
  fEvtMult(0x0),
  fEvtMultvsJetPt(0x0),
  fPtvsEtaJet(0x0),
  fNpvsEtaJet(0x0),
  fNpevtvsEtaJet(0x0),
  fPtvsPtJet(0x0),
  fNpvsPtJet(0x0),
  fNpevtvsPtJet(0x0),
  fPtvsPtJet1D(0x0),
  fNpvsPtJet1D(0x0),
  fNpevtvsPtJet1D(0x0),
  fptLeadingJet(0x0),
  fetaLeadingJet(0x0),
  fphiLeadingJet(0x0),
  fptJet(0x0),
  fetaJet(0x0),
  fphiJet(0x0),
  fHistList(0x0),
  fNBadRuns(0),
  fNBadRunsH(0x0)
 {
  //
  // Default constructor
  //
  /*
  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  }
  */
//   for(int ij  = 0;ij<kMaxJets;++ij){
//     fh1E[ij] = fh1PtRecIn[ij] = fh1PtRecOut[ij] = fh1PtGenIn[ij] = fh1PtGenOut[ij] = 0;
//     fh1Eta[ij] = fh1Phi[ij] = 0;
//   }
}


//#######################################################################
AliAnalysisTaskFragmentationFunction::AliAnalysisTaskFragmentationFunction(const char* name):
  AliAnalysisTaskSE(name),
  fJetHeaderRec(0x0),
  fJetHeaderGen(0x0),
  fAOD(0x0),
  fBranchRec(""),
  fBranchGen(""),
  fUseAODInput(kFALSE),
  fUseAODJetInput(kFALSE),
  fUseAODTrackInput(kFALSE),
  fUseAODMCInput(kFALSE),
  fUseGlobalSelection(kFALSE),
  fUseExternalWeightOnly(kFALSE),
  fLimitGenJetEta(kFALSE),
  fFilterMask(0),
  fAnalysisType(0),
  fTrackTypeRec(kTrackUndef), 
  fTrackTypeGen(kTrackUndef),
  fAvgTrials(1.),
  fExternalWeight(1.),    
  fRecEtaWindow(0.5),
  fR(0.),
  fdRdNdxi(0.7),
  fPartPtCut(0.),
  fEfactor(1.),
  fNff(5),
  fNim(5),
  fList(0x0),
  fGlobVar(1),
  fCDFCut(1),
  // Number of energy bins
  fnEBin(6),         
  fEmin(10.),
  fEmax(70.),
  fnEInterval(6),
  // Number of radius bins
  fnRBin(10),        
  fRmin(0.1),
  fRmax(1.),	
  fnRInterval(9),
  // eta histograms
  fnEtaHBin(50),
  fEtaHBinMin(-0.9),
  fEtaHBinMax(0.9),
  // phi histograms
  fnPhiHBin(60),
  fPhiHBinMin(0.),
  fPhiHBinMax(2*TMath::Pi()),
  // pt histograms
  fnPtHBin(300),
  fPtHBinMin(0.),
  fPtHBinMax(300.),
  // E histograms
  fnEHBin(300),
  fEHBinMin(0.),
  fEHBinMax(300.),
  // Xi histograms
  fnXiHBin(27),
  fXiHBinMin(0.),
  fXiHBinMax(9.),
  // Pthad histograms
  fnPthadHBin(60),
  fPthadHBinMin(0.),
  fPthadHBinMax(30.),
  // z histograms
  fnZHBin(30),
  fZHBinMin(0.),
  fZHBinMax(1.),
  fnPtTrigBin(10),
  // theta histograms
  fnThetaHBin(200),
  fThetaHBinMin(-0.5),
  fThetaHBinMax(1.5),
  fnCosThetaHBin(100),
  fcosThetaHBinMin(0.),
  fcosThetaHBinMax(1.),
  // kT histograms
  fnkTHBin(25),
  fkTHBinMin(0.),
  fkTHBinMax(5.),
  // kT histograms
  fnRHBin(10),
  fRHBinMin(0.),
  fRHBinMax(1.),
  //Histograms
  fEtaMonoJet1H(0x0),
  fPhiMonoJet1H(0x0),
  fPtMonoJet1H(0x0),
  fEMonoJet1H(0x0),
  fdNdXiMonoJet1H(0x0),
  fdNdPtMonoJet1H(0x0),
  fdNdZMonoJet1H(0x0),
  fdNdThetaMonoJet1H(0x0),
  fdNdcosThetaMonoJet1H(0x0),
  fdNdkTMonoJet1H(0x0),
  fdNdpTvsZMonoJet1H(0x0),
  fShapeMonoJet1H(0x0),
  fNMonoJet1sH(0x0),
  fThetaPtPartMonoJet1H(0x0),
  fcosThetaPtPartMonoJet1H(0x0),
  fkTPtPartMonoJet1H(0x0),
  fThetaPtJetMonoJet1H(0x0),
  fcosThetaPtJetMonoJet1H(0x0),
  fkTPtJetMonoJet1H(0x0),
  fpTPtJetMonoJet1H(0x0),
  farrayEmin(0x0),
  farrayEmax(0x0),
  farrayRadii(0x0),
  farrayPtTrigmin(0x0),
  farrayPtTrigmax(0x0),
  // Track control plots
  fptAllTracks(0x0),
  fetaAllTracks(0x0),
  fphiAllTracks(0x0),
  fetaphiptAllTracks(0x0),
  fetaphiAllTracks(0x0),
  fptAllTracksCut(0x0),
  fetaAllTracksCut(0x0),
  fphiAllTracksCut(0x0),
  fetaphiptAllTracksCut(0x0),
  fetaphiAllTracksCut(0x0),
  fptTracks(0x0),
  fetaTracks(0x0),
  fphiTracks(0x0),
  fdetaTracks(0x0),
  fdphiTracks(0x0),
  fetaphiptTracks(0x0),
  fetaphiTracks(0x0),
  fdetadphiTracks(0x0),
  fptTracksCut(0x0),
  fetaTracksCut(0x0),
  fphiTracksCut(0x0),
  fdetaTracksCut(0x0),
  fdphiTracksCut(0x0),
  fetaphiptTracksCut(0x0),
  fetaphiTracksCut(0x0),
  fdetadphiTracksCut(0x0),
  fNPtTrig(0x0),
  fNPtTrigCut(0x0),
  fvertexXY(0x0),
  fvertexZ(0x0),
  fEvtMult(0x0),
  fEvtMultvsJetPt(0x0),
  fPtvsEtaJet(0x0),
  fNpvsEtaJet(0x0),
  fNpevtvsEtaJet(0x0),
  fPtvsPtJet(0x0),
  fNpvsPtJet(0x0),
  fNpevtvsPtJet(0x0),
  fPtvsPtJet1D(0x0),
  fNpvsPtJet1D(0x0),
  fNpevtvsPtJet1D(0x0),
  fptLeadingJet(0x0),
  fetaLeadingJet(0x0),
  fphiLeadingJet(0x0),
  fptJet(0x0),
  fetaJet(0x0),
  fphiJet(0x0),
  fHistList(0x0),
  fNBadRuns(0),
  fNBadRunsH(0x0)
{
  //
  // Default constructor
  //
  /*
  for(int i = 0;i < kMaxStep*2;++i){
    fhnJetContainer[i] = 0;
  } 
  */
//   for(int ij  = 0;ij<kMaxJets;++ij){
//     fh1E[ij] = fh1PtRecIn[ij] = fh1PtRecOut[ij] = fh1PtGenIn[ij] = fh1PtGenOut[ij] = 0;
//     fh1Eta[ij] = fh1Phi[ij] = 0;
//   }
  
  DefineOutput(1, TList::Class());  
}

//////////////////////////////////////////////////////////////////////////////

Bool_t AliAnalysisTaskFragmentationFunction::Notify()
{
  //
  // Implemented Notify() to read the cross sections
  // and number of trials from pyxsec.root
  // 

//   TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
//   UInt_t ntrials  = 0;
//   Float_t ftrials  = 0;
//   if(tree){
//     TFile *curfile = tree->GetCurrentFile();
//     if (!curfile) {
//       Error("Notify","No current file");
//       return kFALSE;
//     }

//     if(!fh1Xsec||!fh1Trials){
//       Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
//       return kFALSE;
//     }

//     TString fileName(curfile->GetName());
//     if(fileName.Contains("AliESDs.root")){
//         fileName.ReplaceAll("AliESDs.root", "pyxsec.root");
//     }
//     else if(fileName.Contains("AliAOD.root")){
//         fileName.ReplaceAll("AliAOD.root", "pyxsec.root");
//     }
//     else if(fileName.Contains("AliAODs.root")){
//         fileName.ReplaceAll("AliAODs.root", "");
//     }
//     else if(fileName.Contains("galice.root")){
//         // for running with galice and kinematics alone...                      
//         fileName.ReplaceAll("galice.root", "pyxsec.root");
//     }
//     TFile *fxsec = TFile::Open(fileName.Data());
//     if(!fxsec){
//       Printf("%s:%d %s not found in the Input",(char*)__FILE__,__LINE__,fileName.Data());
//       // no a severe condition
//       return kTRUE;
//     }
//     TTree *xtree = (TTree*)fxsec->Get("Xsection");
//     if(!xtree){
//       Printf("%s:%d tree not found in the pyxsec.root",(char*)__FILE__,__LINE__);
//     }
//     xtree->SetBranchAddress("xsection",&fXsection);
//     xtree->SetBranchAddress("ntrials",&ntrials);
//     ftrials = ntrials;
//     xtree->GetEntry(0);

//     fh1Xsec->Fill("<#sigma>",fXsection);
//     fh1Trials->Fill("#sum{ntrials}",ftrials);

//   }

  return kTRUE;

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::UserCreateOutputObjects()
{
  //
  // Create the output container
  //

  //**** Connect the AOD
  if(fUseAODInput) // Use AODs as input not ESDs
  { 
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD)
    {
      Printf("%s:%d AODEvent not found in Input Manager %d",(char*)__FILE__,__LINE__,fUseAODInput);
      return;
    }

    // fetch the header
    fJetHeaderRec = dynamic_cast<AliJetHeader*>(fInputHandler->GetTree()->GetUserInfo()->FindObject(Form("AliJetHeader_%s",fBranchRec.Data())));
    if(!fJetHeaderRec)
    {
      Printf("%s:%d Jet Header not found in the Input",(char*)__FILE__,__LINE__);
    }
  }


  else // Use the AOD on the flight 
  {
    //  assume that the AOD is in the general output...
    fAOD  = AODEvent();
    if(!fAOD)
    {
      Printf("%s:%d AODEvent not found in the Output",(char*)__FILE__,__LINE__);
      return;
    }

    //    ((TList*)OutputTree()->GetUserInfo())->Dump();
    fJetHeaderRec = dynamic_cast<AliJetHeader*>(OutputTree()->GetUserInfo()->FindObject(Form("AliJetHeader_%s",fBranchRec.Data())));    
    if(!fJetHeaderRec)
    {
      Printf("%s:%d Jet Header not found in the Output",(char*)__FILE__,__LINE__);
    }
    else
    {
      if(fDebug>10)fJetHeaderRec->Dump();
    }
  }

////////////////////////////

  if (fDebug > 1) printf("AnalysisTaskJetSpectrum::UserCreateOutputObjects() \n");

  OpenFile(1);
  if(!fHistList)fHistList = new TList();
  fHistList->SetOwner(kTRUE);
  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

//////////////////////////////////////////////////
//////////// HISTOGRAM DECLARATION ///////////////
//////////////////////////////////////////////////

  DefineJetH();

//////////////////////////////////////////////////
////////////// HISTOGRAM SAVING //////////////////
//////////////////////////////////////////////////

  for (Int_t i3 = 0; i3 < fnEBin; i3++)
  {
    fHistList->Add(fEtaMonoJet1H[i3]);
    fHistList->Add(fPhiMonoJet1H[i3]);
    fHistList->Add(fPtMonoJet1H[i3]);
    fHistList->Add(fEMonoJet1H[i3]);
   
    for(Int_t i4 = 0; i4 < fnRBin; i4++)
    {
      fHistList->Add(fdNdXiMonoJet1H[i3][i4]);
      fHistList->Add(fdNdPtMonoJet1H[i3][i4]);
      fHistList->Add(fdNdZMonoJet1H[i3][i4]);
      fHistList->Add(fNMonoJet1sH[i3][i4]);
    }
  }
 
  // Theta, kT particles/jet
  for (Int_t i3 = 0; i3 < fnEBin; i3++)
  {
    for(Int_t i4 = 0; i4 < fnRBin; i4++)
    {
      fHistList->Add(fdNdThetaMonoJet1H[i3][i4]);
      fHistList->Add(fdNdcosThetaMonoJet1H[i3][i4]);
      fHistList->Add(fdNdkTMonoJet1H[i3][i4]);
      fHistList->Add(fdNdpTvsZMonoJet1H[i3][i4]);
      fHistList->Add(fShapeMonoJet1H[i3][i4]);
	  
      fHistList->Add(fThetaPtPartMonoJet1H[i3][i4]);
      fHistList->Add(fcosThetaPtPartMonoJet1H[i3][i4]);
      fHistList->Add(fkTPtPartMonoJet1H[i3][i4]);
      fHistList->Add(fThetaPtJetMonoJet1H[i3][i4]);
      fHistList->Add(fcosThetaPtJetMonoJet1H[i3][i4]);
      fHistList->Add(fkTPtJetMonoJet1H[i3][i4]);
      fHistList->Add(fpTPtJetMonoJet1H[i3][i4]);
    }
  }
  
  // Track QA - Correlations
  for (Int_t iPtBin=0; iPtBin<fnPtTrigBin; iPtBin++)
    {
      fHistList->Add(fptTracks[iPtBin]);
      fHistList->Add(fetaTracks[iPtBin]);
      fHistList->Add(fphiTracks[iPtBin]);
      fHistList->Add(fdetaTracks[iPtBin]);
      fHistList->Add(fdphiTracks[iPtBin]);
      fHistList->Add(fetaphiptTracks[iPtBin]);
      fHistList->Add(fetaphiTracks[iPtBin]);
      fHistList->Add(fdetadphiTracks[iPtBin]);
      fHistList->Add(fptTracksCut[iPtBin]);
      fHistList->Add(fetaTracksCut[iPtBin]);
      fHistList->Add(fphiTracksCut[iPtBin]);
      fHistList->Add(fdetaTracksCut[iPtBin]);
      fHistList->Add(fdphiTracksCut[iPtBin]);
      fHistList->Add(fetaphiptTracksCut[iPtBin]);
      fHistList->Add(fetaphiTracksCut[iPtBin]);
      fHistList->Add(fdetadphiTracksCut[iPtBin]);
      fHistList->Add(fNPtTrig[iPtBin]);
      fHistList->Add(fNPtTrigCut[iPtBin]);
    }

  // Track QA 
  fHistList->Add(fptAllTracks);
  fHistList->Add(fetaAllTracks);
  fHistList->Add(fphiAllTracks);
  fHistList->Add(fetaphiptAllTracks);
  fHistList->Add(fetaphiAllTracks);
  fHistList->Add(fptAllTracksCut);
  fHistList->Add(fetaAllTracksCut);
  fHistList->Add(fphiAllTracksCut);
  fHistList->Add(fetaphiptAllTracksCut);
  fHistList->Add(fetaphiAllTracksCut);

  // Event caracterisation QA
  fHistList->Add(fvertexXY);
  fHistList->Add(fvertexZ);
  fHistList->Add(fEvtMult);
  fHistList->Add(fEvtMultvsJetPt);
  fHistList->Add(fPtvsEtaJet);
  fHistList->Add(fNpvsEtaJet);
  fHistList->Add(fNpevtvsEtaJet);
  fHistList->Add(fPtvsPtJet);
  fHistList->Add(fNpvsPtJet);
  fHistList->Add(fNpevtvsPtJet);
  fHistList->Add(fPtvsPtJet1D);
  fHistList->Add(fNpvsPtJet1D);
  fHistList->Add(fNpevtvsPtJet1D);
  fHistList->Add(fptLeadingJet);
  fHistList->Add(fetaLeadingJet);
  fHistList->Add(fphiLeadingJet);
  fHistList->Add(fptJet);
  fHistList->Add(fetaJet);
  fHistList->Add(fphiJet);
  fHistList->Add(fNBadRunsH);

//////////////////////////////////////////////////
///////// END OF HISTOGRAM DECLARATION ///////////
//////////////////////////////////////////////////

  // =========== Switch on Sumw2 for all histos ===========
  for (Int_t i=0; i<fHistList->GetEntries(); ++i) 
  {
    TH1 *h1 = dynamic_cast<TH1*>(fHistList->At(i));
    if (h1)
    {
      // Printf("%s ",h1->GetName());
      h1->Sumw2();
      continue;
    }
  }
  
  TH1::AddDirectory(oldStatus);
  
  /*
    if (fDebug > 1) printf("AnalysisTaskDiJets::CreateOutPutData() \n");
    fDiJets = new TClonesArray("AliAODDiJet", 0);
    fDiJets->SetName("Dijets");
    AddAODBranch("TClonesArray", &fDiJets);
  */
  
  
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::Init()
{
  //
  // Initialization
  //

  Printf(">>> AnalysisTaskFragmentationFunction::Init() debug level %d\n",fDebug);
  if (fDebug > 1) printf("AnalysisTaskDiJets::Init() \n");
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::UserExec(Option_t */*option*/) 
{
  //
  // Execute analysis for current event
  //

  //****
  //**** Check of input data
  //****

  printf("Analysing event # %5d\n", (Int_t) fEntry);
  if (fDebug > 1)printf("Analysing event # %5d\n", (Int_t) fEntry);
  
  AliAODHandler *aodH = dynamic_cast<AliAODHandler*>(AliAnalysisManager::GetAnalysisManager()->GetOutputEventHandler());
  if(!aodH)
  {
    Printf("%s:%d no output aodHandler found Jet",(char*)__FILE__,__LINE__);
    return;
  }
  
  TClonesArray *aodRecJets = dynamic_cast<TClonesArray*>(fAOD->FindListObject(fBranchRec.Data()));
  if(!aodRecJets)
  {
    Printf("%s:%d no reconstructed Jet array with name %s in AOD",(char*)__FILE__,__LINE__,fBranchRec.Data());
    return;
  }
  if (fDebug > 10) Printf("%s:%d",(char*)__FILE__,__LINE__);
  
  //****
  //**** Check of primary vertex
  //****
  AliAODVertex * pvtx = dynamic_cast<AliAODVertex*>(fAOD->GetPrimaryVertex());
  if( (!pvtx) ||
      (pvtx->GetZ()<-10. || pvtx->GetZ()>10.) ||
      (pvtx->GetNContributors()<0) )
    { 
	fNBadRuns++;
	fNBadRunsH->Fill(0.5);
 	return;
    }

  //****
  //**** Check number of tracks 
  //****
  TClonesArray* tracks = dynamic_cast<TClonesArray*>(fAOD->GetTracks());

  //****
  //**** Declaration of arrays and variables 
  //****
  // We use static array, not to fragment the memory
  AliAODJet recJets[kMaxJets];
  Int_t nRecJets       = 0;
  Int_t nTracks = 0;

//////////////////////////////////////////////////
///////// Get the reconstructed jets /////////////
//////////////////////////////////////////////////

  nRecJets = aodRecJets->GetEntries();
  nRecJets = TMath::Min(nRecJets,kMaxJets);
  nTracks  = fAOD->GetNumberOfTracks();
    
  for(int ir = 0;ir < nRecJets;++ir)
  {
    AliAODJet *tmp = dynamic_cast<AliAODJet*>(aodRecJets->At(ir));
    if(!tmp)continue;
    recJets[ir] = *tmp;
    cout << "recJets[" << ir << "].Eta(): " << recJets[ir].Eta() << ", recJets[" << ir <<"].Phi(): " << recJets[ir].Phi() << ", recJets[" << ir << "].E(): " << recJets[ir].E() << endl;
  }
  if(nRecJets>1) {
    Float_t detaJ = recJets[0].Eta() - recJets[1].Eta();
    Float_t dphiJ = recJets[0].Phi() - recJets[1].Phi();
    Float_t detJ = recJets[0].Pt() - recJets[1].Pt();
    cout << "detaJ: " << detaJ << ", dphiJ: " << dphiJ << ", detJ: " << detJ << endl;
  }

  // Get vertex information
  fvertexXY->Fill(pvtx->GetX(),pvtx->GetY());
  fvertexZ->Fill(pvtx->GetZ());

  //////////////////////////////////////////////////
  ////////////// TRACK QUALITY ASSURANCE ///////////           TO BE OPTIMISED!!
  //////////////////////////////////////////////////
  Int_t evtMult = 0;
  for(Int_t it=0; it<nTracks; it++)
  {
    AliAODTrack* aodTrack = (AliAODTrack*)tracks->At(it);
    Float_t etaT = aodTrack->Eta();
    Float_t phiT = aodTrack->Phi();
    Float_t ptT = aodTrack->Pt(); 
    phiT = ((phiT < 0) ? phiT + 2 * TMath::Pi() : phiT);
    //    cout << "etaT: " << etaT << ", phiT: " << phiT << endl;
    fptAllTracks->Fill(ptT);
    fetaAllTracks->Fill(etaT);
    fphiAllTracks->Fill(phiT);
    fetaphiptAllTracks->Fill(etaT,phiT,ptT);
    fetaphiAllTracks->Fill(etaT,phiT,1);
    UInt_t status = aodTrack->GetStatus();

    for(Int_t i=0; i<fnPtTrigBin; i++)
    {
      if(ptT>=farrayPtTrigmin[i] && ptT<farrayPtTrigmax[i])
      {
	fptTracks[i]->Fill(ptT);
	fetaTracks[i]->Fill(etaT);
	fphiTracks[i]->Fill(phiT);
	fNPtTrig[i]->Fill(0.5);
	// Compute deta/dphi
	Float_t etaT2 = 0.; 
	Float_t phiT2 = 0.;
	Float_t ptT2 = 0.;
	Float_t detaT = 0.;
	Float_t dphiT = 0.;
	for(Int_t it2 = 0; it2< nTracks; it2++)
	{
	   AliAODTrack* aodTrack2 = (AliAODTrack*)tracks->At(it2);
	   etaT2 = aodTrack2->Eta(); phiT2 = aodTrack2->Phi(); ptT2 = aodTrack2->Pt();
	   phiT2 = ((phiT2 < 0) ? phiT2 + 2 * TMath::Pi() : phiT2);
	   //	      cout << "etaT2: " << etaT2 << ", phiT2: " << phiT2 << endl;
	   if(ptT2 > 2.*i+4.) continue;
	   if(it2==it) continue;
	   detaT = etaT - etaT2; 
	   dphiT = phiT - phiT2; 
	   if (dphiT  >   TMath::Pi()) dphiT = (-TMath::Pi() +TMath::Abs(dphiT - TMath::Pi()));
	   if (dphiT  < -1.0*TMath::Pi())  dphiT = (TMath::Pi() -  TMath::Abs(dphiT + TMath::Pi()));
	      
	   fdetaTracks[i]->Fill(detaT);
	   fdphiTracks[i]->Fill(dphiT);
	   fdetadphiTracks[i]->Fill(detaT,dphiT,1);
	}
        fetaphiptTracks[i]->Fill(etaT,phiT,ptT);
	fetaphiTracks[i]->Fill(etaT,phiT,1);
      }
    } // End loop over trigger ranges 

    if (status == 0) continue;
    if((fFilterMask>0)&&!(aodTrack->TestFilterBit(fFilterMask))) continue;
    fptAllTracksCut->Fill(ptT);
    fetaAllTracksCut->Fill(etaT);
    fphiAllTracksCut->Fill(phiT);
    fetaphiptAllTracksCut->Fill(etaT,phiT,ptT);
    fetaphiAllTracksCut->Fill(etaT,phiT,1);
    if(ptT > 0.150 && TMath::Abs(etaT) < 0.9) evtMult++;
  } // end loop over tracks
  fEvtMult->Fill(evtMult);   

//////////////////////////////////////////////////
///////////////// MONOJET PART ///////////////////
//////////////////////////////////////////////////

  if (nRecJets == 0) return;

  Double_t jetEnergy = recJets[0].E();
  Int_t    goodBin   = 999;

  for (Int_t i1 = 0; i1 < fnEBin; i1++)
  {
    if (jetEnergy < farrayEmax[i1] && jetEnergy >= farrayEmin[i1]) 
    {
      goodBin = i1;
      continue;
    }
  }

  fptLeadingJet->Fill(recJets[0].Pt());
  fetaLeadingJet->Fill(recJets[0].Eta());
  fphiLeadingJet->Fill(recJets[0].Phi());

  for(Int_t ij=0; ij<nRecJets; ij++)
  {  
    fptJet->Fill(recJets[ij].Pt());
    fetaJet->Fill(recJets[ij].Eta());
    fphiJet->Fill(recJets[ij].Phi());
  }

  // Get track ref 
  TRefArray* ref = recJets[0].GetRefTracks();
  for(Int_t it=0; it<ref->GetEntries(); it++)
  {
    Float_t ptTrack = ((AliVTrack*)ref->At(it))->Pt();
    fPtvsEtaJet->Fill(recJets[0].Eta(),ptTrack);
    fNpvsEtaJet->Fill(recJets[0].Eta(),ref->GetEntries());
    fNpevtvsEtaJet->Fill(recJets[0].Eta(),evtMult);
    fPtvsPtJet->Fill(recJets[0].Pt(),ptTrack);
    fNpvsPtJet->Fill(recJets[0].Pt(),ref->GetEntries());
    fNpevtvsPtJet->Fill(recJets[0].Pt(),evtMult);
    fPtvsPtJet1D->Fill(recJets[0].Pt(),ptTrack);
    fNpvsPtJet1D->Fill(recJets[0].Pt(),ref->GetEntries());
    fNpevtvsPtJet1D->Fill(recJets[0].Pt(),evtMult);
  }

  FillMonoJetH(goodBin, recJets, tracks);

  //////////////////////////////////////////////////
  ////////////////// DIJET PART ////////////////////       
  //////////////////////////////////////////////////

  // UNDER CONSTRUCTION

  PostData(1, fHistList);
}

//#######################################################################
void AliAnalysisTaskFragmentationFunction::Terminate(Option_t */*option*/)
{
// Terminate analysis
//
    if (fDebug > 1) printf("AnalysisDiJets: Terminate() \n");
    printf("Number of events with vertex out of bound: %d", fNBadRuns);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::DefineJetH()
{

  /////////////////////////////////////// HISTOGRAMS FIRST JET
  fEtaMonoJet1H   = new TH1F*[fnEBin+1];
  fPhiMonoJet1H   = new TH1F*[fnEBin+1];
  fPtMonoJet1H    = new TH1F*[fnEBin+1];
  fEMonoJet1H     = new TH1F*[fnEBin+1];

  fdNdXiMonoJet1H = new TH1F**[fnEBin+1];
  fdNdPtMonoJet1H = new TH1F**[fnEBin+1];
  fdNdZMonoJet1H  = new TH1F**[fnEBin+1];
  fdNdThetaMonoJet1H    = new TH1F**[fnEBin+1];
  fdNdcosThetaMonoJet1H    = new TH1F**[fnEBin+1];
  fdNdkTMonoJet1H       = new TH1F**[fnEBin+1];
  fdNdpTvsZMonoJet1H    = new TH1F**[fnEBin+1];
  fShapeMonoJet1H       = new TH1F**[fnEBin+1];
  fNMonoJet1sH          = new TH1F**[fnEBin+1];

  fThetaPtPartMonoJet1H = new TH2F**[fnEBin+1];
  fcosThetaPtPartMonoJet1H = new TH2F**[fnEBin+1];
  fkTPtPartMonoJet1H    = new TH2F**[fnEBin+1];
  fThetaPtJetMonoJet1H = new TH2F**[fnEBin+1];
  fcosThetaPtJetMonoJet1H = new TH2F**[fnEBin+1];
  fkTPtJetMonoJet1H    = new TH2F**[fnEBin+1];
  fpTPtJetMonoJet1H    = new TH2F**[fnEBin+1];

  for (Int_t iEbin=0;iEbin<fnEBin+1;iEbin++)
  {
    fdNdXiMonoJet1H[iEbin]       = new TH1F*[fnRBin+1];
    fdNdPtMonoJet1H[iEbin]       = new TH1F*[fnRBin+1];
    fdNdZMonoJet1H[iEbin]        = new TH1F*[fnRBin+1];
    fdNdThetaMonoJet1H[iEbin]    = new TH1F*[fnRBin+1];
    fdNdcosThetaMonoJet1H[iEbin] = new TH1F*[fnRBin+1];
    fdNdkTMonoJet1H[iEbin]       = new TH1F*[fnRBin+1];
    fdNdpTvsZMonoJet1H[iEbin]    = new TH1F*[fnRBin+1];
    fShapeMonoJet1H[iEbin]       = new TH1F*[fnRBin+1];
    fNMonoJet1sH[iEbin]          = new TH1F*[fnRBin+1];

    fThetaPtPartMonoJet1H[iEbin]    = new TH2F*[fnRBin+1];
    fcosThetaPtPartMonoJet1H[iEbin] = new TH2F*[fnRBin+1];
    fkTPtPartMonoJet1H[iEbin]       = new TH2F*[fnRBin+1];
    fThetaPtJetMonoJet1H[iEbin]     = new TH2F*[fnRBin+1];
    fcosThetaPtJetMonoJet1H[iEbin]  = new TH2F*[fnRBin+1];
    fkTPtJetMonoJet1H[iEbin]        = new TH2F*[fnRBin+1];
    fpTPtJetMonoJet1H[iEbin]        = new TH2F*[fnRBin+1];
  }

  for (Int_t iEbin=0;iEbin<fnEBin+1;iEbin++)
  {
    fEtaMonoJet1H[iEbin]    = 0;
    fPhiMonoJet1H[iEbin]    = 0;
    fPtMonoJet1H[iEbin]     = 0;
    fEMonoJet1H[iEbin]      = 0;

    for (Int_t iRbin=0;iRbin<fnRBin+1;iRbin++)
    {
      fdNdXiMonoJet1H[iEbin][iRbin]       = 0;
      fdNdPtMonoJet1H[iEbin][iRbin]       = 0;
      fdNdZMonoJet1H[iEbin][iRbin]        = 0;
      fdNdThetaMonoJet1H[iEbin][iRbin]    = 0;
      fdNdcosThetaMonoJet1H[iEbin][iRbin] = 0;
      fdNdkTMonoJet1H[iEbin][iRbin]       = 0;
      fdNdpTvsZMonoJet1H[iEbin][iRbin]    = 0;
      fShapeMonoJet1H[iEbin][iRbin]       = 0;
      fNMonoJet1sH[iEbin][iRbin]          = 0;

      fThetaPtPartMonoJet1H[iEbin][iRbin]    = 0;
      fcosThetaPtPartMonoJet1H[iEbin][iRbin] = 0;
      fkTPtPartMonoJet1H[iEbin][iRbin]       = 0;
      fThetaPtJetMonoJet1H[iEbin][iRbin]     = 0;
      fcosThetaPtJetMonoJet1H[iEbin][iRbin]  = 0;
      fkTPtJetMonoJet1H[iEbin][iRbin]        = 0;
      fpTPtJetMonoJet1H[iEbin][iRbin]        = 0;
    }
  }

  fptTracks          = new TH1F*[fnPtTrigBin+1];
  fetaTracks         = new TH1F*[fnPtTrigBin+1];
  fphiTracks         = new TH1F*[fnPtTrigBin+1];
  fdetaTracks        = new TH1F*[fnPtTrigBin+1];
  fdphiTracks        = new TH1F*[fnPtTrigBin+1];
  fetaphiptTracks    = new TH2F*[fnPtTrigBin+1];
  fetaphiTracks      = new TH2F*[fnPtTrigBin+1];
  fdetadphiTracks    = new TH2F*[fnPtTrigBin+1];
  fptTracksCut       = new TH1F*[fnPtTrigBin+1];
  fetaTracksCut      = new TH1F*[fnPtTrigBin+1];
  fphiTracksCut      = new TH1F*[fnPtTrigBin+1];
  fdetaTracksCut     = new TH1F*[fnPtTrigBin+1];
  fdphiTracksCut     = new TH1F*[fnPtTrigBin+1];
  fetaphiptTracksCut = new TH2F*[fnPtTrigBin+1];
  fetaphiTracksCut   = new TH2F*[fnPtTrigBin+1];
  fdetadphiTracksCut = new TH2F*[fnPtTrigBin+1];
  fNPtTrig           = new TH1F*[fnPtTrigBin+1];
  fNPtTrigCut        = new TH1F*[fnPtTrigBin+1];

  for(Int_t iPtTrigBin=0; iPtTrigBin<fnPtTrigBin; iPtTrigBin++)
  {
    fptTracks[iPtTrigBin]          = 0;
    fetaTracks[iPtTrigBin]         = 0;
    fphiTracks[iPtTrigBin]         = 0;
    fdetaTracks[iPtTrigBin]        = 0;
    fdphiTracks[iPtTrigBin]        = 0;
    fetaphiptTracks[iPtTrigBin]    = 0;
    fetaphiTracks[iPtTrigBin]      = 0;
    fdetadphiTracks[iPtTrigBin]    = 0;
    fptTracksCut[iPtTrigBin]       = 0;
    fetaTracksCut[iPtTrigBin]      = 0;
    fphiTracksCut[iPtTrigBin]      = 0;
    fdetaTracksCut[iPtTrigBin]     = 0;
    fdphiTracksCut[iPtTrigBin]     = 0;
    fetaphiptTracksCut[iPtTrigBin] = 0;
    fetaphiTracksCut[iPtTrigBin]   = 0;
    fdetadphiTracksCut[iPtTrigBin] = 0;
    fNPtTrig[iPtTrigBin]           = 0;
    fNPtTrigCut[iPtTrigBin]        = 0;
  }

  farrayEmin = new Double_t[fnEBin];
  farrayEmax = new Double_t[fnEBin];

  farrayPtTrigmin = new Double_t[fnPtTrigBin];
  farrayPtTrigmax = new Double_t[fnPtTrigBin];

  farrayRadii = new Double_t[fnRBin];

  Double_t eMin    = 0;
  Double_t eMax    = 0;
  Double_t pasE    = (Double_t)((fEmax-fEmin)/fnEInterval);
  TString energy;
  TString energy2;

  Double_t r    = 0.;
  Double_t pasR = (Double_t)((fRmax-fRmin)/fnRInterval);
  TString radius;

  for (Int_t i = 0; i < fnEBin; i++)
  {
    eMin       = 0;
    eMax       = 0;

    if (i==0) eMin = fEmin;
    if (i!=0) eMin = fEmin + pasE*i;
    eMax = eMin+pasE;
    energy2 = "E_{jet1} : ";
    energy = "E_{jet} : ";
    energy += eMin;
    energy2 += eMin;
    energy += "-";
    energy2 += "-";
    energy +=eMax;
    energy2 += eMax;
    energy += "GeV";
    energy2 += "GeV";

    farrayEmin[i]      = eMin;
    farrayEmax[i]      = eMax;

    for (Int_t j = 0; j < fnRBin; j++)
    {
 
      if (j==0) r = fRmin;
      if (j!=0) r = fRmin + pasR*j;
      radius = ", R = ";
      radius += r;    

      farrayRadii[j] = r;

      fEtaMonoJet1H[i]   = new TH1F("fEtaMonoJet1H,"+energy, "#eta_{jet1},"+energy, fnEtaHBin, fEtaHBinMin, fEtaHBinMax);
      fPhiMonoJet1H[i]   = new TH1F("fPhiMonoJet1H,"+energy, "#phi_{jet1},"+energy, fnPhiHBin, fPhiHBinMin, fPhiHBinMax);
      fPtMonoJet1H[i]    = new TH1F("fPtMonoJet1H,"+energy, "pT_{jet1},"+energy, fnPtHBin, fPtHBinMin, fPtHBinMax);
      fEMonoJet1H[i]     = new TH1F("fEMonoJet1H,"+energy, "E_{jet1},"+energy, fnEHBin, fEHBinMin, fEHBinMax);

      fdNdXiMonoJet1H[i][j] = new TH1F("fdNdXiMonoJet1H,"+energy+radius, "dN_{ch}/d#xi,"+energy+radius, fnXiHBin, fXiHBinMin, fXiHBinMax);
      fdNdPtMonoJet1H[i][j] = new TH1F("fdNdPtMonoJet1H,"+energy+radius, "dN_{ch}/dPt_{had},"+energy+radius, fnPthadHBin, fPthadHBinMin, fPthadHBinMax);
      fdNdZMonoJet1H[i][j]  = new TH1F("fdNdZMonoJet1H,"+energy+radius, "dN_{ch}/dz,"+energy+radius, fnZHBin, fZHBinMin, fZHBinMax);

      fdNdThetaMonoJet1H[i][j]    = new TH1F("fdNdThetaMonoJet1H,"+energy+radius, "dN_{ch}/d#Theta,"+energy+radius, fnThetaHBin, fThetaHBinMin, fThetaHBinMax);
      fdNdcosThetaMonoJet1H[i][j] = new TH1F("fdNdcosThetaMonoJet1H,"+energy+radius, "dN_{ch}/dcos(#Theta),"+energy+radius, fnCosThetaHBin, fcosThetaHBinMin, fcosThetaHBinMax);
      fdNdkTMonoJet1H[i][j]   = new TH1F("fdNdkTMonoJet1H,"+energy+radius, "dN_{ch}/dk_{T},"+energy+radius, fnkTHBin, fkTHBinMin, fkTHBinMax);
      fdNdpTvsZMonoJet1H[i][j]= new TH1F("fdNdpTvsZMonoJet1H,"+energy+radius, "dN_{ch}/dk_{T} vs Z,"+energy+radius, fnZHBin, fZHBinMin, fZHBinMax);
      fShapeMonoJet1H[i][j]   = new TH1F("fShapeMonoJet1H,"+energy+radius, "E(R=x)/E(R=1)"+energy+radius, fnRHBin, fRHBinMin, fRHBinMax);

      fThetaPtPartMonoJet1H[i][j] = new TH2F("fThetaPtPartMonoJet1H,"+energy+radius, "#Theta vs Pt particle,"+energy+radius, fnPthadHBin, fPthadHBinMin, fPthadHBinMax, fnThetaHBin, fThetaHBinMin, fThetaHBinMax);
      fcosThetaPtPartMonoJet1H[i][j] = new TH2F("fcosThetaPtPartMonoJet1H,"+energy+radius, "cos(#Theta) vs Pt particle,"+energy+radius, fnPthadHBin, fPthadHBinMin, fPthadHBinMax, fnCosThetaHBin, fcosThetaHBinMin, fcosThetaHBinMax);
      fkTPtPartMonoJet1H[i][j] = new TH2F("fkTPtPartMonoJet1H,"+energy+radius, "kT vs Pt particle,"+energy+radius, fnPthadHBin, fPthadHBinMin, fPthadHBinMax, fnkTHBin, fkTHBinMin, fkTHBinMax);
      fThetaPtJetMonoJet1H[i][j] = new TH2F("fThetaPtJetMonoJet1H,"+energy+radius, "#Theta vs Pt jet,"+energy+radius, fnPtHBin, fPtHBinMin, fPtHBinMax, fnThetaHBin, fThetaHBinMin, fThetaHBinMax);
      fcosThetaPtJetMonoJet1H[i][j] = new TH2F("fcosThetaPtJetMonoJet1H,"+energy+radius, "cos(#Theta) vs Pt jet,"+energy+radius, fnPtHBin, fPtHBinMin, fPtHBinMax, fnCosThetaHBin, fcosThetaHBinMin, fcosThetaHBinMax);
      fkTPtJetMonoJet1H[i][j] = new TH2F("fkTPtJetMonoJet1H,"+energy+radius, "kT vs Pt jet,"+energy+radius, fnPtHBin, fPtHBinMin, fPtHBinMax, fnkTHBin, fkTHBinMin, fkTHBinMax);
      fpTPtJetMonoJet1H[i][j] = new TH2F("fpTPtJetMonoJet1H,"+energy+radius, "pT vs Pt jet,"+energy+radius, fnPtHBin, fPtHBinMin, fPtHBinMax, fnkTHBin, fkTHBinMin, fkTHBinMax);

      fNMonoJet1sH[i][j]    = new TH1F("fNMonoJet1sH,"+energy+radius, "N_{jets1},"+energy+radius, 1, 0., 1.);
      fNBadRunsH = new TH1F("fNBadRunsH","Number of events with Z vertex out of range", 1, 0., 1.);

      SetProperties(fEtaMonoJet1H[i], "#eta_{jet1}", "Entries");
      SetProperties(fPhiMonoJet1H[i], "#phi_{jet1}", "Entries");
      SetProperties(fPtMonoJet1H[i], "p_{Tjet1} (GeV/c)", "Entries");
      SetProperties(fEMonoJet1H[i], "E_{jet1} (GeV)", "Entries");

      SetProperties(fdNdXiMonoJet1H[i][j], "#xi = ln(E_{jet1}/p_{Thad})", "dN_{had}/d#xi");
      SetProperties(fdNdPtMonoJet1H[i][j], "p_{Thad} (GeV/c)", "dN_{had}/dp_{Thad}");
      SetProperties(fdNdZMonoJet1H[i][j], "z = (p_{Thad}/E_{jet1})", "dN_{had}/dz");
      SetProperties(fdNdThetaMonoJet1H[i][j], "#Theta", "dN_{had}/d#Theta");
      SetProperties(fdNdcosThetaMonoJet1H[i][j], "cos(#Theta)", "dN_{had}/dcos(#Theta)");
      SetProperties(fdNdkTMonoJet1H[i][j], "k_{Thad}", "dN_{had}/dk_{Thad}");
      SetProperties(fdNdpTvsZMonoJet1H[i][j], "z = (p_{Thad}/E_{jet1})", "dN_{had}/dp_{T}");
      SetProperties(fShapeMonoJet1H[i][j], "R", "#Psi(R)");
      SetProperties(fNMonoJet1sH[i][j], "Bin", "N_{jets1}");

      SetProperties(fThetaPtPartMonoJet1H[i][j],"p_{Thad} [GeV/c]","#Theta");
      SetProperties(fcosThetaPtPartMonoJet1H[i][j],"p_{Thad} [GeV/c]","cos(#Theta)");
      SetProperties(fkTPtPartMonoJet1H[i][j],"p_{Thad} [GeV/c]","k_{Thad}");
      SetProperties(fThetaPtJetMonoJet1H[i][j], "p_{Tjet} [GeV/c]", "#Theta");
      SetProperties(fcosThetaPtJetMonoJet1H[i][j], "p_{Tjet} [GeV/c]", "cos(#Theta)");
      SetProperties(fkTPtJetMonoJet1H[i][j], "p_{Tjet} [GeV/c]", "k_{Thad} [GeV/c]");
      SetProperties(fpTPtJetMonoJet1H[i][j], "p_{Tjet} [GeV/c]", "p_{Thad} [GeV/c]");
    }
  }

  for(Int_t i=0; i<fnPtTrigBin; i++)
  {
    if(i==0) farrayPtTrigmin[i] = 1.; 
    else farrayPtTrigmin[i] = i*5.;
    farrayPtTrigmax[i] = i*5+5.;
    
    TString ptTrigRange;
    ptTrigRange = "; p_{T} trig range: ";
    ptTrigRange += farrayPtTrigmin[i];
    ptTrigRange += "-";
    ptTrigRange += farrayPtTrigmax[i];
    ptTrigRange += " [GeV]";

    fptTracks[i] = new TH1F("fptTracks"+ptTrigRange, "Track transverse momentum [GeV]"+ptTrigRange,300,0.,150.);
    fetaTracks[i] = new TH1F("fetaTracks"+ptTrigRange, "#eta tracks"+ptTrigRange, 36, -0.9, 0.9);
    fphiTracks[i] = new TH1F("fphiTracks"+ptTrigRange, "#phi tracks"+ptTrigRange, 60, 0., 2*TMath::Pi());
    fdetaTracks[i] = new TH1F("fdetaTracks"+ptTrigRange, "#Delta #eta tracks"+ptTrigRange,80, -2., 2.);
    fdphiTracks[i] = new TH1F("fdphiTracks"+ptTrigRange, "#Delta #phi tracks"+ptTrigRange, 120, -TMath::Pi(), TMath::Pi());
    fetaphiptTracks[i] = new TH2F("fetaphiptTracks"+ptTrigRange,"#eta/#phi track p_{T} mapping"+ptTrigRange,36, -0.9, 0.9,60, 0., 2*TMath::Pi());
    fetaphiTracks[i] = new TH2F("fetaphiTracks"+ptTrigRange,"#eta/#phi track mapping"+ptTrigRange,36, -0.9, 0.9,60, 0., 2*TMath::Pi());
    fdetadphiTracks[i] = new TH2F("fdetadphiTracks"+ptTrigRange,"#Delta #eta/#Delta #phi track mapping"+ptTrigRange,80, -2., 2., 120, -TMath::Pi(), TMath::Pi());
    fptTracksCut[i] = new TH1F("fptTracksCut"+ptTrigRange, "Track transverse momentum after cuts [GeV]"+ptTrigRange,300,0.,150.);
    fetaTracksCut[i] = new TH1F("fetaTracksCut"+ptTrigRange, "#eta tracks after cuts"+ptTrigRange, 36, -0.9, 0.9);
    fphiTracksCut[i] = new TH1F("fphiTracksCuts"+ptTrigRange, "#phi tracks after cuts"+ptTrigRange, 60, 0., 2*TMath::Pi());
    fdetaTracksCut[i] = new TH1F("fdetaTracksCuts"+ptTrigRange, "#Delta #eta tracks after cuts"+ptTrigRange,80, -2., 2.);
    fdphiTracksCut[i] = new TH1F("fdphiTracksCuts"+ptTrigRange, "#Delta #phi tracks after cuts"+ptTrigRange, 120, -TMath::Pi(), TMath::Pi());
    fetaphiptTracksCut[i] = new TH2F("fetaphiptTracksCuts"+ptTrigRange,"#eta/#phi track p_{T} mapping after cuts"+ptTrigRange,36, -0.9, 0.9,60, 0., 2*TMath::Pi());
    fetaphiTracksCut[i] = new TH2F("fetaphiTracksCuts"+ptTrigRange,"#eta/#phi track mapping after cuts"+ptTrigRange,36, -0.9, 0.9,60, 0., 2*TMath::Pi());
    fdetadphiTracksCut[i] = new TH2F("fdetadphiTracksCuts"+ptTrigRange,"#Delta #eta/#Delta #phi track mapping after cuts"+ptTrigRange,80, -2., 2., 120, -TMath::Pi(), TMath::Pi());
    fNPtTrig[i] = new TH1F("fNPtTrig"+ptTrigRange,"Number of triggers"+ptTrigRange,1,0.,1.);
    fNPtTrigCut[i] = new TH1F("fNPtTrigCut"+ptTrigRange,"Number of triggers after cut"+ptTrigRange,1,0.,1.);

    SetProperties(fptTracks[i], "Track p_{T} [GeV]", "dN/dp_{T}"); 
    SetProperties(fetaTracks[i], "Track #eta", "dN/d#eta"); 
    SetProperties(fphiTracks[i], "Track #phi", "dN/d#phi");
    SetProperties(fdetaTracks[i], "#eta_{track} - #eta_{trig}", "dN/d#Delta#eta");
    SetProperties(fdphiTracks[i], "#phi_{track} - #phi_{trig}", "dN/d#Delta#phi");
    SetProperties(fetaphiptTracks[i], "#eta_{track}", "#phi_{track}"); 
    SetProperties(fetaphiTracks[i], "#eta_{track}", "#phi_{track}");
    SetProperties(fdetadphiTracks[i], "#Delta #eta_{track}", "#Delta #phi_{track}");
    SetProperties(fptTracksCut[i], "p_{T}track [GeV]", "dN/dp_{T}");
    SetProperties(fetaTracksCut[i], "#eta_{track}", "dN/d#eta"); 
    SetProperties(fphiTracksCut[i], "#phi_{track}", "dN/d#phi"); 
    SetProperties(fdetaTracksCut[i], "#eta_{track} - #eta_{trig}", "dN/d#Delta#eta");
    SetProperties(fdphiTracksCut[i], "#phi_{track} - #phi_{trig}", "dN/d#Delta#phi");
    SetProperties(fetaphiptTracksCut[i], "#eta_{track}", "#phi_{track}");
    SetProperties(fetaphiTracksCut[i], "#eta_{track}", "#phi_{track}");
    SetProperties(fdetadphiTracksCut[i], "#Delta #eta_{track}", "#Delta #phi_{track}");
    SetProperties(fNPtTrig[i], "", "Number of triggers");
    SetProperties(fNPtTrigCut[i], "", "Number of triggers");

  }

  fptAllTracks = new TH1F("fptAllTracks", "Track transverse momentum [GeV]",300,0.,150.);
  fetaAllTracks = new TH1F("fetaAllTracks", "#eta tracks", 36, -0.9, 0.9);
  fphiAllTracks = new TH1F("fphiAllTracks", "#phi tracks", 60, 0., 2*TMath::Pi());
  fetaphiptAllTracks = new TH2F("fetaphiptAllTracks","#eta/#phi track p_{T} mapping",36, -0.9, 0.9,60, 0., 2*TMath::Pi());
  fetaphiAllTracks = new TH2F("fetaphiAllTracks","#eta/#phi track mapping",36, -0.9, 0.9,60, 0., 2*TMath::Pi());
  fptAllTracksCut = new TH1F("fptAllTracksCut", "Track transverse momentum after cuts [GeV]",300,0.,150.);
  fetaAllTracksCut = new TH1F("fetaAllTracksCut", "#eta tracks after cuts", 36, -0.9, 0.9);
  fphiAllTracksCut = new TH1F("fphiAllTracksCuts", "#phi tracks after cuts", 60, 0., 2*TMath::Pi());
  fetaphiptAllTracksCut = new TH2F("fetaphiptAllTracksCuts","#eta/#phi track p_{T} mapping after cuts",36, -0.9, 0.9,60, 0., 2*TMath::Pi());
  fetaphiAllTracksCut = new TH2F("fetaphiAllTracksCuts","#eta/#phi track mapping after cuts",36, -0.9, 0.9,60, 0., 2*TMath::Pi());

  SetProperties(fptAllTracks, "Track p_{T} [GeV]", "dN/dp_{T}"); 
  SetProperties(fetaAllTracks, "Track #eta", "dN/d#eta");
  SetProperties(fphiAllTracks, "Track #phi", "dN/d#phi");
  SetProperties(fetaphiptAllTracks, "#eta_{track}", "#phi_{track}");
  SetProperties(fetaphiAllTracks, "#eta_{track}", "#phi_{track}");
  SetProperties(fptAllTracksCut, "p_{T}track [GeV]", "dN/dp_{T}"); 
  SetProperties(fetaAllTracksCut, "#eta_{track}", "dN/d#eta"); 
  SetProperties(fphiAllTracksCut, "#phi_{track}", "dN/d#phi"); 
  SetProperties(fetaphiptAllTracksCut, "#eta_{track}", "#phi_{track}");
  SetProperties(fetaphiAllTracksCut, "#eta_{track}", "#phi_{track}");

  fvertexXY = new TH2F("fvertexXY","X-Y vertex position",30,0.,10.,30,0.,10.);
  fvertexZ = new TH1F("fvertexZ","Z vertex position",60,-30.,30.);
  fEvtMult = new TH1F("fEvtMult","Event multiplicity, track pT cut > 150 MeV/c, |#eta| < 0.9",100,0.,100.);
  fEvtMultvsJetPt = new TH2F("fEvtMultvsJetPt","Event multiplicity vs pT_{jet}",60,0.,60.,100,0.,100.);
  fPtvsEtaJet = new TH2F("fPtvsEtaJet","Pt vs #eta_{jet}",20,-1.,1.,60,0.,60.);
  fNpvsEtaJet = new TH2F("fNpvsEtaJet","N_{part} inside jet vs #eta_{jet}",20,-1.,1.,20,0,20); 
  fNpevtvsEtaJet = new TH2F("fNpevtvsEtaJet","N_{part} in evt vs #eta_{jet}",20,-1.,1.,90,0,90); 
  fPtvsPtJet = new TH2F("fPtvsPtJet","Pt vs #p_{Tjet}",60,0.,60.,60,0.,60.);
  fNpvsPtJet = new TH2F("fNpvsPtJet","N_{part} inside jet vs #pt_{Tjet}",60,0.,60.,20,0,20); 
  fNpevtvsPtJet = new TH2F("fNpevtvsPtJet","N_{part} in evt vs #pt_{Tjet}",60,0.,60.,90,0,90); 
  fPtvsPtJet1D = new TH1F("fPtvsPtJet1D","Pt vs #p_{Tjet}",60,0.,60.);
  fNpvsPtJet1D = new TH1F("fNpvsPtJet1D","N_{part} inside jet vs #pt_{Tjet}",60,0.,60.); 
  fNpevtvsPtJet1D = new TH1F("fNpevtvsPtJet1D","N_{part} in evt vs #pt_{Tjet}",60,0.,60.); 
  fptLeadingJet = new TH1F("fptLeadingJet","Pt leading Jet [GeV/c]",60,0.,60.); 
  fetaLeadingJet = new TH1F("fetaLeadingJet","#eta leading jet",20,-1.,1.);
  fphiLeadingJet = new TH1F("fphiLeadingJet","#phi leading jet",12,0.,2*TMath::Pi());
  fptJet = new TH1F("fptJet","Pt Jets [GeV/c]",60,0.,60.); 
  fetaJet = new TH1F("fetaJet","#eta jet",20,-1.,1.);
  fphiJet = new TH1F("fphiJet","#phi jet",12,0.,2*TMath::Pi());
  fNBadRunsH = new TH1F("fNBadRunsH","Number of events with Z vertex out of range", 1, 0., 1.);

  SetProperties(fvertexXY, "vtx X", "vtx Y");
  SetProperties(fvertexZ, "vtx Z", "Count");
  SetProperties(fEvtMult, "N_{part} / event", "Count");
  SetProperties(fEvtMultvsJetPt, "p_{T jet}", "Event multiplicity");
  SetProperties(fPtvsEtaJet, "#eta_{leading jet}", "p_{T} part [GeV/c]");
  SetProperties(fNpvsEtaJet, "#eta_{leading jet}", "N_{part} in leading jet");
  SetProperties(fNpevtvsEtaJet, "#eta_{leading jet}", "N_{part} in event");
  SetProperties(fPtvsPtJet, "#p_{T leading jet}", "p_{T} part [GeV/c]");
  SetProperties(fNpvsPtJet, "#p_{T leading jet}", "N_{part} in leading jet");
  SetProperties(fNpevtvsPtJet, "#p_{T leading jet}", "N_{part} in event");
  SetProperties(fPtvsPtJet1D, "#p_{T leading jet}", "<p_{T}> part [GeV/c]");
  SetProperties(fNpvsPtJet1D, "#p_{T leading jet}", "<N_{part}> in leading jet");
  SetProperties(fNpevtvsPtJet1D, "#p_{T leading jet}", "<N_{part}> in event");
  SetProperties(fptLeadingJet, "p_{T} leading jet", "dN/dp_{T} leading jet");
  SetProperties(fetaLeadingJet, "#eta leading jet", "dN/d#eta leading jet");
  SetProperties(fphiLeadingJet, "#phi leading jet", "dN/d#phi leading jet");
  SetProperties(fptJet, "p_{T} jet [GeV/c]", "dN/dp_{T}");
  SetProperties(fetaJet, "#eta jet", "dN/d#eta");
  SetProperties(fphiJet, "#phi jet", "dN/d#phi");

}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::FillMonoJetH(Int_t goodBin, AliAODJet* recJets, TClonesArray* tracks)
{
  if (goodBin == 999) return;

  Int_t nTracks = tracks->GetEntries();
  Float_t xi,t1,ene,dr2,deta2,dphi2, z, cosTheta, theta, kt;
  Bool_t jetOk1 = 0;
  Bool_t jetOk2 = 0;
  Bool_t jetOk3 = 0;
  Bool_t jetOk4 = 0;
  Bool_t jetOk5 = 0;
  Bool_t jetOk6 = 0;
  Bool_t jetOk7 = 0;
  Bool_t jetOk8 = 0;
  Bool_t jetOk9 = 0;
  Bool_t jetOk10 = 0;

  fEtaMonoJet1H[goodBin]->Fill(recJets[0].Eta());
  fPhiMonoJet1H[goodBin]->Fill(recJets[0].Phi());
  fPtMonoJet1H[goodBin]->Fill(recJets[0].Pt());
  fEMonoJet1H[goodBin]->Fill(recJets[0].E());

  Int_t mult = 0;
  for(Int_t it=0; it<nTracks; it++)
  {
    AliAODTrack* aodTrack = (AliAODTrack*)tracks->At(it);
    // Apply track cuts
    UInt_t status = aodTrack->GetStatus();
    if (status == 0) continue;
    if((fFilterMask>0)&&!(aodTrack->TestFilterBit(fFilterMask)))continue;
    mult++;
    Float_t etaT = aodTrack->Eta();
    Float_t phiT = aodTrack->Phi();
    Float_t ptT = aodTrack->Pt();
    // For Theta distribution
    Float_t pxT = aodTrack->Px();
    Float_t pyT = aodTrack->Py();
    Float_t pzT = aodTrack->Pz();
    Float_t pT = aodTrack->P();
    // Compute theta
    cosTheta = (pxT*recJets[0].Px()+pyT*recJets[0].Py()+pzT*recJets[0].Pz())/(pT*recJets[0].P());
    theta = TMath::ACos(cosTheta);
    // Compute xi
    deta2 = etaT - recJets[0].Eta();
    dphi2 = phiT - recJets[0].Phi();
    if (dphi2 < -TMath::Pi()) dphi2= -dphi2 - 2.0 * TMath::Pi();
    if (dphi2 > TMath::Pi()) dphi2 = 2.0 * TMath::Pi() - dphi2;
    dr2 = TMath::Sqrt(deta2 * deta2 + dphi2 * dphi2);
    t1  = TMath::Tan(2.0*TMath::ATan(TMath::Exp(-etaT)));
    ene = ptT*TMath::Sqrt(1.+1./(t1*t1));
    xi  = (Float_t) TMath::Log(recJets[0].E()/ptT);
    // Compute z
    z   = (Double_t)(ptT/recJets[0].E());
    // Compute kT/jT
    TVector3 partP; TVector3 jetP;
    jetP[0] = recJets[0].Px();
    jetP[1] = recJets[0].Py();
    jetP[2] = recJets[0].Pz();
    partP.SetPtEtaPhi(ptT,etaT,phiT);
    kt = TMath::Sin(partP.Angle(jetP))*partP.Mag();
    // Compute Jet shape

    for(Int_t i2 = 0; i2 < fnRBin; i2++)
    {
      if ((dr2<farrayRadii[i2]) && (ptT > fPartPtCut)) 
      {
        if (i2 == 0) jetOk1 = 1;
        if (i2 == 1) jetOk2 = 1;
        if (i2 == 2) jetOk3 = 1;
        if (i2 == 3) jetOk4 = 1;
        if (i2 == 4) jetOk5 = 1;
        if (i2 == 5) jetOk6 = 1;
        if (i2 == 6) jetOk7 = 1;
        if (i2 == 7) jetOk8 = 1;
        if (i2 == 8) jetOk9 = 1;
        if (i2 == 9) jetOk10 = 1;

        fdNdXiMonoJet1H[goodBin][i2]->Fill(xi);
        fdNdPtMonoJet1H[goodBin][i2]->Fill(ptT);
        fdNdZMonoJet1H[goodBin][i2]->Fill(z);
	fdNdThetaMonoJet1H[goodBin][i2]->Fill(theta);
	fdNdcosThetaMonoJet1H[goodBin][i2]->Fill(cosTheta);
	fdNdkTMonoJet1H[goodBin][i2]->Fill(kt);
	fdNdpTvsZMonoJet1H[goodBin][i2]->Fill(z,1/((fPthadHBinMax-fPthadHBinMin)/fnPthadHBin));

	fThetaPtPartMonoJet1H[goodBin][i2]->Fill(ptT,theta);
	fcosThetaPtPartMonoJet1H[goodBin][i2]->Fill(ptT,cosTheta);
	fkTPtPartMonoJet1H[goodBin][i2]->Fill(ptT,kt);
	fThetaPtJetMonoJet1H[goodBin][i2]->Fill(recJets[0].Pt(),theta);
	fcosThetaPtJetMonoJet1H[goodBin][i2]->Fill(recJets[0].Pt(),cosTheta);
	fkTPtJetMonoJet1H[goodBin][i2]->Fill(recJets[0].Pt(),kt);
	fpTPtJetMonoJet1H[goodBin][i2]->Fill(recJets[0].Pt(),ptT);
      }
    }
  }
  fEvtMultvsJetPt->Fill(recJets[0].Pt(),mult);

  if (jetOk1)  fNMonoJet1sH[goodBin][0]->Fill(0.5);
  if (jetOk2)  fNMonoJet1sH[goodBin][1]->Fill(0.5);
  if (jetOk3)  fNMonoJet1sH[goodBin][2]->Fill(0.5);
  if (jetOk4)  fNMonoJet1sH[goodBin][3]->Fill(0.5);
  if (jetOk5)  fNMonoJet1sH[goodBin][4]->Fill(0.5);
  if (jetOk6)  fNMonoJet1sH[goodBin][5]->Fill(0.5);
  if (jetOk7)  fNMonoJet1sH[goodBin][6]->Fill(0.5);
  if (jetOk8)  fNMonoJet1sH[goodBin][7]->Fill(0.5);
  if (jetOk9)  fNMonoJet1sH[goodBin][8]->Fill(0.5);
  if (jetOk10) fNMonoJet1sH[goodBin][9]->Fill(0.5);

}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

void AliAnalysisTaskFragmentationFunction::SetProperties(TH1* h,const char* x, const char* y)
{
  //Set properties of histos (x and y title and error propagation)
  h->SetXTitle(x);
  h->SetYTitle(y);
  h->GetXaxis()->SetTitleColor(1);
  h->GetYaxis()->SetTitleColor(1);
  h->Sumw2();
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
