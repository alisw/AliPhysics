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

///////////////////////////////////////////////////////////////////////////
/// AliAnalysisTaskSE for the extraction of the various histograms to
/// study the pt spectra of identified hadrons vs multiplicity:
/// - multiplicity estimated with Reference Multiplicity and "V0M" percentiles
/// - log(dEdx)-log(dEdxBB) distributions for pions, kaons and protons in pt bins
/// - Pt distributions of pions, kaons and protons with nSigma PID
/// Authors: 
/// \author E. Biolcati, biolcati@to.infn.it
/// \author L. Milano,   milano@to.infn.it
/// \author F. Prino,    prino@to.infn.it
/// \author N. Jacazio,  jacazio@to.infn.it
///////////////////////////////////////////////////////////////////////////

#define LOG_NO_INFO
// #define LOG_NO_DEBUG

#include <TH1F.h>
#include <TF1.h>
#include <TRandom3.h>
#include <TH2F.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TParticle.h>
#include <Rtypes.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisDataContainer.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliPhysicsSelection.h"
#include "AliAnalysisTaskSEITSsaSpectraMultiplicity.h"
#include "AliESDtrackCuts.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliESD.h"
#include "AliITSPIDResponse.h"
#include "THnSparse.h"
#include "AliPPVsMultUtils.h"
#include "AliESDUtils.h"
#include "AliAnalysisUtils.h"
#include "TStopwatch.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "AliLog.h"

ClassImp(AliAnalysisTaskSEITSsaSpectraMultiplicity)
/* $Id$ */

//________________________________________________________________________
AliAnalysisTaskSEITSsaSpectraMultiplicity::AliAnalysisTaskSEITSsaSpectraMultiplicity():
AliAnalysisTaskSE("taskITSsaSpectra"),
fESDtrackCuts(0),
fPPVsMultUtils(0),
fUtils(0),
fESD(0),
fOutput(0),
fListCuts(0),
fHistNEvents(0),
fHistMCSampSel(0),
fHistMult(0),
fHistMultAftEvSel(0),
fHistEventMultiplicity(0),
fHistBefPileUp(0),
fHistTaggedPileUp(0),
fHistSparse(0),
fHistSparseBefEvSel(0),
fHistPtCorr(0x0),
fHistCen(0),
fHistNTracks(0),
fHistNTracksPos(0),
fHistNTracksNeg(0),
fHistDEDX(0),
fHistDEDXdouble(0),
fHistDEDXMulti(0),
fHistDEDXdoubleMulti(0),
fHistBeforeEvSel(0),
fHistAfterEvSel(0),
fHistITSchi2ncls(0),
fDCAxyCutFunc(0x0),
fDCAzCutFunc(0x0),
fITSPIDResponse(0),
fMinSPDPts(1),
fMinNdEdxSamples(3),
fMindEdx(0.),
fMinNSigma(1.5),
fMaxY(0.5),
fMaxChi2Clu(2.5),
fNSigmaDCAxy(7.),
fNSigmaDCAz(7.),
fEtaRange(0.8),
fvZRange(10.),
fLowMult(-1),
fUpMult(-1),
fLowCentrality(-1.0),
fUpCentrality(-1.0),
fMultEstimator(0),
fHImode(0),
fPileupContributors(3),
fPileupDistance(0.8),
fYear(2010),
fMC(kFALSE), 
fSmearMC(kFALSE),
fSmearP(0.),
fSmeardEdx(0.),
fRandGener(0),
fFillTree(kFALSE),
fLowEnergypp(kFALSE),
fuseV0(kFALSE),
fPileupRej(kTRUE),
fUseThnSparse(kFALSE),
fTreeNSigma(0x0),
fTreeMC(0x0),
fMult(0),
fV0Mult(0),
fdedx(0),
fP(0),
fTrkSign(0),
fClumap(0),
fEta(0),
fevSelmask(0),
fImpactXY(0),
fImpactZ(0),
fParticleMap(0),
fptMC(0),
fChi2(0)
{
  AliLog::SetClassDebugLevel("AliAnalysisTaskSEITSsaSpectraMultiplicity",2);
  infomsg("Begin of AliAnalysisTaskSEITSsaSpectraMultiplicity");
  
  // Constructor
  Double_t xbins[kNbins+1]={0.08,0.10,0.12,0.14,0.16,0.18,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0};
  for(Int_t iBin=0; iBin<kNbins+1; iBin++) fPtBinLimits[iBin]=xbins[iBin];
  fRandGener=new TRandom3(0);
  
  for(Int_t j=0; j<3; j++){
    fHistPosNSigmaSep[j]=0x0;
    fHistNegNSigmaSep[j]=0x0;
    
    fHistPrimMCpos[j]=0x0;
    fHistPrimMCneg[j]=0x0;
    fHistPrimMCposEtaY[j]=0x0;
    fHistPrimMCnegEtaY[j]=0x0;
    fHistSecStrMCpos[j]=0x0;
    fHistSecStrMCneg[j]=0x0;
    fHistSecMatMCpos[j]=0x0;
    fHistSecMatMCneg[j]=0x0;
    fHistPrimMCposBefEvSel[j]=0x0;
    fHistPrimMCposBefEvSelEtaY[j]=0x0;
    fHistPrimMCposBefEvSelEta[j]=0x0;
    fHistPrimMCnegBefEvSel[j]=0x0;
    fHistPrimMCnegBefEvSelEtaY[j]=0x0;
    fHistPrimMCnegBefEvSelEta[j]=0x0;
    fHistSecStrMCposBefEvSel[j]=0x0;
    fHistSecStrMCnegBefEvSel[j]=0x0;
    fHistSecMatMCposBefEvSel[j]=0x0;
    fHistSecMatMCnegBefEvSel[j]=0x0;
    fHistPrimMCposReco[j]=0x0;
    fHistPrimMCposRecoEtaY[j]=0x0;
    fHistPrimMCnegReco[j]=0x0;
    fHistPrimMCnegRecoEtaY[j]=0x0;
    fHistSecStrMCposReco[j]=0x0;
    fHistSecStrMCnegReco[j]=0x0;
    fHistSecMatMCposReco[j]=0x0;
    fHistSecMatMCnegReco[j]=0x0;
  }
  
  for(Int_t j=0; j<4; j++) fHistCharge[j]=0x0;
  
  for(Int_t j=0; j<kNbins; j++){
    fHistPosPi[j]=0x0;
    fHistPosK[j]=0x0;
    fHistPosP[j]=0x0;
    fHistNegPi[j]=0x0;
    fHistNegK[j]=0x0;
    fHistNegP[j]=0x0;
    fHistDCAPosPi[j]=0x0;
    fHistDCAPosK[j]=0x0;
    fHistDCAPosP[j]=0x0;
    fHistDCANegPi[j]=0x0;
    fHistDCANegK[j]=0x0;
    fHistDCANegP[j]=0x0;
    fHistMCPrimDCAPosPi[j]=0x0;
    fHistMCPrimDCAPosK[j]=0x0;
    fHistMCPrimDCAPosP[j]=0x0;
    fHistMCPrimDCANegPi[j]=0x0;
    fHistMCPrimDCANegK[j]=0x0;
    fHistMCPrimDCANegP[j]=0x0;
    fHistMCSecStDCAPosPi[j]=0x0;
    fHistMCSecStDCAPosK[j]=0x0;
    fHistMCSecStDCAPosP[j]=0x0;
    fHistMCSecStDCANegPi[j]=0x0;
    fHistMCSecStDCANegK[j]=0x0;
    fHistMCSecStDCANegP[j]=0x0;
    fHistMCSecMatDCAPosPi[j]=0x0;
    fHistMCSecMatDCAPosK[j]=0x0;
    fHistMCSecMatDCAPosP[j]=0x0;
    fHistMCSecMatDCANegPi[j]=0x0;
    fHistMCSecMatDCANegK[j]=0x0;
    fHistMCSecMatDCANegP[j]=0x0;
    fHistMCPosOtherHypPion[j]=0x0;
    fHistMCPosOtherHypKaon[j]=0x0;
    fHistMCPosOtherHypProton[j]=0x0;
    fHistMCPosElHypPion[j]=0x0;
    fHistMCPosElHypKaon[j]=0x0;
    fHistMCPosElHypProton[j]=0x0;
    fHistMCPosPiHypPion[j]=0x0;
    fHistMCPosPiHypKaon[j]=0x0;
    fHistMCPosPiHypProton[j]=0x0;
    fHistMCPosKHypPion[j]=0x0;
    fHistMCPosKHypKaon[j]=0x0;
    fHistMCPosKHypProton[j]=0x0;
    fHistMCPosPHypPion[j]=0x0;
    fHistMCPosPHypKaon[j]=0x0;
    fHistMCPosPHypProton[j]=0x0;
    fHistMCNegOtherHypPion[j]=0x0;
    fHistMCNegOtherHypKaon[j]=0x0;
    fHistMCNegOtherHypProton[j]=0x0;
    fHistMCNegElHypPion[j]=0x0;
    fHistMCNegElHypKaon[j]=0x0;
    fHistMCNegElHypProton[j]=0x0;
    fHistMCNegPiHypPion[j]=0x0;
    fHistMCNegPiHypKaon[j]=0x0;
    fHistMCNegPiHypProton[j]=0x0;
    fHistMCNegKHypPion[j]=0x0;
    fHistMCNegKHypKaon[j]=0x0;
    fHistMCNegKHypProton[j]=0x0;
    fHistMCNegPHypPion[j]=0x0;
    fHistMCNegPHypKaon[j]=0x0;
    fHistMCNegPHypProton[j]=0x0;
  }
  
  for(Int_t j=0; j<3; j++){
    fHistPosNSigmaMean[j]=0x0;
    fHistPosNSigmaMCMean[j]=0x0;
    fHistPosNSigmaPrimMean[j]=0x0;
    fHistPosNSigmaPrimMCMean[j]=0x0;
    fHistNegNSigmaMean[j]=0x0;
    fHistNegNSigmaMCMean[j]=0x0;
    fHistNegNSigmaPrimMean[j]=0x0;
    fHistNegNSigmaPrimMCMean[j]=0x0;
    
    fHistPosNSigma[j]=0x0;
    fHistPosNSigmaMC[j]=0x0;
    fHistPosNSigmaPrim[j]=0x0;
    fHistPosNSigmaPrimMC[j]=0x0;
    fHistNegNSigma[j]=0x0;
    fHistNegNSigmaMC[j]=0x0;
    fHistNegNSigmaPrim[j]=0x0;
    fHistNegNSigmaPrimMC[j]=0x0;
  }
  
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  infomsg("End of AliAnalysisTaskSEITSsaSpectraMultiplicity");
}

//___________________________________________________________________________
AliAnalysisTaskSEITSsaSpectraMultiplicity::~AliAnalysisTaskSEITSsaSpectraMultiplicity(){
  // Destructor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if (fListCuts) {
    delete fListCuts;
    fListCuts = 0;
  }
  
  if(fRandGener) delete fRandGener;
  if(fITSPIDResponse) delete fITSPIDResponse;
  
  if (fPPVsMultUtils){
    delete fPPVsMultUtils;
    fPPVsMultUtils = 0x0;
  }
  if (fUtils){
    delete fUtils;
    fUtils = 0x0;
  }
  
  delete fDCAxyCutFunc;
  delete fDCAzCutFunc;
}

//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectraMultiplicity::CookdEdx(Double_t *s) const {
  // truncated mean for the dEdx
  Int_t nc=0; 
  Double_t dedx[4]={0.,0.,0.,0.};
  for (Int_t il=0; il<4; il++) { // count good (>0) dE/dx values
    if(s[il]>fMindEdx){
      dedx[nc]= s[il];
      nc++;
    }
  }
  if(nc<fMinNdEdxSamples) return -1.;
  
  Double_t tmp;
  Int_t swap; // sort in ascending order 
  do {
    swap=0;
    for (Int_t i=0; i<nc-1; i++) {
      if (dedx[i]<=dedx[i+1]) continue;
      tmp=dedx[i];
      dedx[i]=dedx[i+1];
      dedx[i+1]=tmp;
      swap++;
    } 
  } while (swap);
  
  Double_t sumamp=0,sumweight=0;
  Double_t weight[4]={1.,1.,0.,0.};
  if(nc==3) weight[1]=0.5;
  else if(nc<3) weight[1]=0.;
  for (Int_t i=0; i<nc; i++) {
    sumamp+= dedx[i]*weight[i];
    sumweight+=weight[i];
  }
  return sumamp/sumweight;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectraMultiplicity::DCAcut(Double_t impactXY, Double_t impactZ, Double_t pt) const {
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  
  return DCAcutXY(impactXY, pt)*DCAcutZ(impactZ, pt);
}


//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectraMultiplicity::DCAcutXY(Double_t impactXY, Double_t pt) const {
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  
  Double_t xyMax = fDCAxyCutFunc->Eval(pt); //in micron
  AliDebugF(3,"Max value for the DCAxy Cut is:%f Measured value for DCAxy is:%f cut to %.0f sigmas\n",xyMax,TMath::Abs(impactXY)*10000,fDCAxyCutFunc->GetParameter(3));
  if((TMath::Abs(impactXY)*10000)>xyMax) return kFALSE;  
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSEITSsaSpectraMultiplicity::DCAcutZ(Double_t impactZ, Double_t pt) const {
  // cut on transverse impact parameter updated on 20-5-2010
  // from the study of L. Milano, F. Prino on the ITS standalone tracks
  // using the common binning of the TPC tracks
  
  Double_t zMax = fDCAzCutFunc->Eval(pt); //in micron
  AliDebugF(3,"Max value for the DCAz Cut is:%f Measured value for DCAz is:%f cut to %.0f sigmas\n",zMax,TMath::Abs(impactZ)*10000,fDCAzCutFunc->GetParameter(3));
  if((TMath::Abs(impactZ)*10000)>zMax) return kFALSE;
  
  return kTRUE;
}


//________________________________________________________________________
Double_t AliAnalysisTaskSEITSsaSpectraMultiplicity::Eta2y(Double_t pt, Double_t m, Double_t eta) const {
  // convert eta to y
  Double_t mt = TMath::Sqrt(m*m + pt*pt);
  return TMath::ASinH(pt/mt*TMath::SinH(eta));
}

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectraMultiplicity::Init(){
  //
  // Initialization
  //
  infomsg("Begin of Init");
  
  fListCuts=new TList();
  fListCuts->SetOwner();
  Double_t xyP[3];
  Double_t zP[3];
  if(fMC){
    if(fYear==2009){
      xyP[0]=88.63; //MC LHC10a12
      xyP[1]=19.57;
      xyP[2]=1.65;
      zP[0]=140.98;
      zP[1]=62.33;
      zP[2]=1.15;
    }else{
      xyP[0]=36.; //MC LHC10d1
      xyP[1]=43.9;
      xyP[2]=1.3;
      zP[0]=111.9;
      zP[1]=59.8;
      zP[2]=1.2;
    }
  }else{
    if(fYear==2009){
      xyP[0]=85.28;//DATA 900 GeV pass6
      xyP[1]=25.78;
      xyP[2]=1.55;
      zP[0]=146.80;
      zP[1]=70.07;
      zP[2]=1.11;
    }else{
      xyP[0]=32.7;//DATA 7 TeV pass2
      xyP[1]=44.8;
      xyP[2]=1.3;
      zP[0]=117.3;
      zP[1]=66.8;
      zP[2]=1.2;
    }
  }
  fDCAxyCutFunc = new TF1("fDCAxyCutFunc","[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))",0.05,10.);
  for(Int_t ipar=0; ipar<3; ipar++) fDCAxyCutFunc->SetParameter(ipar,xyP[ipar]);
  infomsg(Form("Setting DCA cut value to %.0f",fNSigmaDCAxy));
  fDCAxyCutFunc->SetParameter(3,fNSigmaDCAxy);
  
  fDCAzCutFunc = new TF1("fDCAzCutFunc","[3]*([0]+[1]/TMath::Power(TMath::Abs(x),[2]))",0.05,10.);
  for(Int_t ipar=0; ipar<3; ipar++) fDCAzCutFunc->SetParameter(ipar,zP[ipar]);
  infomsg(Form("Setting DCA cut value to %.0f",fNSigmaDCAz));
  fDCAzCutFunc->SetParameter(3,fNSigmaDCAz);
  
  fListCuts->Add(fDCAxyCutFunc);
  fListCuts->Add(fDCAzCutFunc);
  
  PostData(2,fListCuts);
  infomsg("End of Init()");
}

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectraMultiplicity::UserCreateOutputObjects(){
  // Create a TList with histograms and a TNtuple
  // Called once
  infomsg("Begin of CreateOutputObjects");
  
  //General Utilities  
  // Multiplicity
  if(! fESDtrackCuts ){
    fESDtrackCuts = new AliESDtrackCuts();
  }
  //Helper
  if(! fPPVsMultUtils ){
    fPPVsMultUtils = new AliPPVsMultUtils();
  }
  //Analysis Utils
  if(! fUtils ){
    fUtils = new AliAnalysisUtils();
  }
  //PID Response
  if(!fITSPIDResponse){
    fITSPIDResponse = new AliITSPIDResponse(fMC); 
  }
  
  fOutput = new TList();
  fOutput->SetOwner();
  fOutput->SetName("Spiderman");
  
  fHistNEvents = new TH1F("fHistNEvents", "Number of processed events;Ev. Sel. Step;Counts",11,-1.5,9.5);
  fHistNEvents->Sumw2();
  fHistNEvents->SetMinimum(0);
  fHistNEvents->GetXaxis()->SetBinLabel(1,"Read from ESD");
  fHistNEvents->GetXaxis()->SetBinLabel(2,"Pass Phys. Sel. + Trig");
  fHistNEvents->GetXaxis()->SetBinLabel(3,"SDD read out");
  fHistNEvents->GetXaxis()->SetBinLabel(4,"In mult. range");
  fHistNEvents->GetXaxis()->SetBinLabel(5,"IsMinimumBias");
  fHistNEvents->GetXaxis()->SetBinLabel(6,"HasGoodVertex");
  fHistNEvents->GetXaxis()->SetBinLabel(7,"IsAcceptedVertexPosition");
  fHistNEvents->GetXaxis()->SetBinLabel(8,"IsINELgtZERO");
  fHistNEvents->GetXaxis()->SetBinLabel(9,"IsNotPileupSPDInMultBins");
  fHistNEvents->GetXaxis()->SetBinLabel(10,"HasNoInconsistentSPDandTrackVertices");
  fHistNEvents->GetXaxis()->SetBinLabel(11,"IsEventSelected");
  fOutput->Add(fHistNEvents);
  
  fHistMCSampSel = new TH1F("fHistMCSampSel", "MC Sample",3,-1.5,1.5);
  fHistMCSampSel->Sumw2();
  fHistMCSampSel->SetMinimum(0);
  fHistMCSampSel->GetXaxis()->SetBinLabel(1,"Read from ESD");
  fHistMCSampSel->GetXaxis()->SetBinLabel(2,"In mult. range");
  fHistMCSampSel->GetXaxis()->SetBinLabel(3,"Pass Samp. Sel.");
  fOutput->Add(fHistMCSampSel);
  
  //Two cases: if it's run on DATA 5 axis, if it's run on MC 7 axis (MCpt and MCPID)
  UInt_t dimsparse;
  if(!fMC) dimsparse=5;
  else dimsparse=7;
  Int_t binsparse[dimsparse];
  Double_t minsparse[dimsparse];
  Double_t maxsparse[dimsparse];
  
  const UInt_t numberofv0multibin=12;
  Double_t v0multbin[numberofv0multibin+1] = {0.,0.01,0.1,1.,5.,10.,15.,20.,30.,40.,50.,70.,100.};
  
  if(!fMC){//Five axis: multiplicity (Nch or V0), pt, hipotesis particle asym, hipotesis particle sym, DCA 
    
    if(!fuseV0){
      binsparse[0]=3000;//Multiplicity
      minsparse[0]=-0.5;
      maxsparse[0]=2999.5;
    }
    else{
      binsparse[0]=numberofv0multibin;//Multiplicity V0
      minsparse[0]=0.;
      maxsparse[0]=100;
    }
    
    binsparse[1]=kNbins;//Pt
    minsparse[1]=0;
    maxsparse[1]=kNbins;
    
    binsparse[2]=12;//Asymmetric Particle Hipotesis plus information from the DCA cut: [0,1] not DCA approved [1,2] DCA approved et ansi de suite
    minsparse[2]=-6;
    maxsparse[2]=6;
    
    binsparse[3]=6;//Symmetric Particle Hipotesis Pi-K-P
    minsparse[3]=-3;
    maxsparse[3]=3;
    
    binsparse[4]=2000;//DCA
    minsparse[4]=-1;
    maxsparse[4]=1;
  }
  else {//Seven axis: multiplicity (Nch or V0), pt, hipotesis particle asym, hipotesis particle sym, DCA , MC Hipotesis Prim;SecSt;SetMat on pt, MCpt
    
    if(!fuseV0){
      binsparse[0]=3000;//Multiplicity
      minsparse[0]=-0.5;
      maxsparse[0]=2999.5;
    }
    else{
      binsparse[0]=numberofv0multibin;//Multiplicity V0
      minsparse[0]=0.;
      maxsparse[0]=100;
    }
    
    binsparse[1]=kNbins;//Pt
    minsparse[1]=0;
    maxsparse[1]=kNbins;
    
    binsparse[2]=12;//Particle Hipotesis asym
    minsparse[2]=-6;
    maxsparse[2]=6;
    
    binsparse[3]=6;//Particle Hipotesis sym
    minsparse[3]=-3;
    maxsparse[3]=3;
    
    binsparse[4]=2000;//DCA
    minsparse[4]=-1;
    maxsparse[4]=1;
    
    binsparse[5]=21;//To identificate the type of particle there are 3 bins for each one (pos and neg) plus the partial information from primaries in this order: Prim-SecSt-SecMat
    minsparse[5]=-9;
    maxsparse[5]=12;
    
    binsparse[6]=kNbins;//MC Pt
    minsparse[6]=0;
    maxsparse[6]=kNbins;
  }
  
  fHistSparse = new THnSparseD("fHistSparse","Combination of data",dimsparse,binsparse,minsparse,maxsparse);
  if(!fuseV0) fHistSparse->GetAxis(0)->SetTitle("Event Multiplicity Nch");
  else{
    fHistSparse->GetAxis(0)->Set(numberofv0multibin,v0multbin);//Rebin V0 multiplicity axis
    fHistSparse->GetAxis(0)->SetTitle("Event Multiplicity V0");
  }
  fHistSparse->GetAxis(1)->Set(kNbins,fPtBinLimits);//Rebin Pt axis
  
  fHistSparse->GetAxis(1)->SetTitle("Particle transverse momentum");
  fHistSparse->GetAxis(2)->SetTitle("Identified Particle aSym");
  fHistSparse->GetAxis(2)->SetBinLabel(1,"Negative Proton DCA passed");
  fHistSparse->GetAxis(2)->SetBinLabel(2,"Negative Proton");
  fHistSparse->GetAxis(2)->SetBinLabel(3,"Negative Kaon DCA passed");
  fHistSparse->GetAxis(2)->SetBinLabel(4,"Negative Kaon");
  fHistSparse->GetAxis(2)->SetBinLabel(5,"Negative Pion DCA passed");
  fHistSparse->GetAxis(2)->SetBinLabel(6,"Negative Pion");
  fHistSparse->GetAxis(2)->SetBinLabel(7,"Positive Pion");
  fHistSparse->GetAxis(2)->SetBinLabel(8,"Positive Pion DCA passed");
  fHistSparse->GetAxis(2)->SetBinLabel(9,"Positive Kaon");
  fHistSparse->GetAxis(2)->SetBinLabel(10,"Positive Kaon DCA passed");
  fHistSparse->GetAxis(2)->SetBinLabel(11,"Positive Proton");
  fHistSparse->GetAxis(2)->SetBinLabel(12,"Positive Proton DCA passed");
  fHistSparse->GetAxis(3)->SetTitle("Identified Particle Sym");
  fHistSparse->GetAxis(3)->SetBinLabel(1,"Negative Proton");
  fHistSparse->GetAxis(3)->SetBinLabel(2,"Negative Kaon");
  fHistSparse->GetAxis(3)->SetBinLabel(3,"Negative Pion");
  fHistSparse->GetAxis(3)->SetBinLabel(4,"Positive Pion");
  fHistSparse->GetAxis(3)->SetBinLabel(5,"Positive Kaon");
  fHistSparse->GetAxis(3)->SetBinLabel(6,"Positive Proton");
  fHistSparse->GetAxis(4)->SetTitle("DCAxy");
  
  if(fMC){
    fHistSparse->GetAxis(5)->SetTitle("MC Identified Particle and Production");
    fHistSparse->GetAxis(5)->SetBinLabel(1,"MC Negative Proton From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(2,"MC Negative Kaon From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(3,"MC Negative Pion From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(4,"MC Negative Proton From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(5,"MC Negative Kaon From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(6,"MC Negative Pion From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(7,"MC Negative Proton From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(8,"MC Negative Kaon From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(9,"MC Negative Pion From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(10,"MC Positive Pion From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(11,"MC Positive Kaon From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(12,"MC Positive Proton From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(13,"MC Positive Pion From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(14,"MC Positive Kaon From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(15,"MC Positive Proton From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(16,"MC Positive Pion From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(17,"MC Positive Kaon From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(18,"MC Positive Proton From Material");
    fHistSparse->GetAxis(5)->SetBinLabel(19,"MC Positive Not identified From Primary");
    fHistSparse->GetAxis(5)->SetBinLabel(20,"MC Positive Not identified From Strangeness");
    fHistSparse->GetAxis(5)->SetBinLabel(21,"MC Positive Not identified From Material");
    fHistSparse->GetAxis(6)->Set(kNbins,fPtBinLimits);//Rebin
    fHistSparse->GetAxis(6)->SetTitle("Particle transverse MC momentum");
  }
  
  fHistSparse->Sumw2();
  if(fUseThnSparse) fOutput->Add(fHistSparse);
  
  if(fMC){//four axis: multiplicity (Nch or V0), MCpt, hipotesis particle (with primary Y and eta cut information) and event selection
    UInt_t dimsparseMC=4;
    Int_t binsparseMC[dimsparseMC];
    Double_t minsparseMC[dimsparseMC];
    Double_t maxsparseMC[dimsparseMC];
    if(!fuseV0){
      binsparseMC[0]=3000;//Multiplicity Reference
      minsparseMC[0]=-0.5;
      maxsparseMC[0]=2999.5;
    }
    else{
      binsparseMC[0]=numberofv0multibin;//Multiplicity V0
      minsparseMC[0]=0.;
      maxsparseMC[0]=100;
    }
    
    binsparseMC[1]=kNbins;//MCPt
    minsparseMC[1]=0;
    maxsparseMC[1]=kNbins;
    
    binsparseMC[2]=30;//Particle Hipotesis 18 bins for Primaries (with eta & y cut information) and 12 bins for Secondaries (strangeness & material)
    minsparseMC[2]=-15;
    maxsparseMC[2]=15;
    
    binsparseMC[3]=4;//EventSelection
    minsparseMC[3]=-0.5;
    maxsparseMC[3]=3.5;
    
    fHistSparseBefEvSel = new THnSparseD("fHistSparseBefEvSel","Combination of data",dimsparseMC,binsparseMC,minsparseMC,maxsparseMC);
    if(!fuseV0)fHistSparseBefEvSel->GetAxis(0)->SetTitle("Event Multiplicity Nch");
    else{
      fHistSparseBefEvSel->GetAxis(0)->Set(numberofv0multibin,v0multbin);//Rebin V0 multiplicity axis
      fHistSparseBefEvSel->GetAxis(0)->SetTitle("Event Multiplicity V0");
    }
    fHistSparseBefEvSel->GetAxis(1)->Set(kNbins,fPtBinLimits);//Rebin Pt axis
    
    fHistSparseBefEvSel->GetAxis(1)->SetTitle("Particle MC momentum");
    fHistSparseBefEvSel->GetAxis(2)->SetTitle("Hipotesis Particle");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(1,"MC Negative Proton From Material");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(2,"MC Negative Kaon From Material");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(3,"MC Negative Pion From Material");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(4,"MC Negative Proton From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(5,"MC Negative Kaon From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(6,"MC Negative Pion From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(7,"MC Negative Proton From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(8,"MC Negative Proton From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(9,"MC Negative Proton From Primary Eta only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(10,"MC Negative Kaon From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(11,"MC Negative Kaon From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(12,"MC Negative Kaon From Primary Eta only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(13,"MC Negative Pion From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(14,"MC Negative Pion From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(15,"MC Negative Pion From Primary Eta only");
    
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(16,"MC Positive Pion From Primary Eta only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(17,"MC Positive Pion From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(18,"MC Positive Pion From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(19,"MC Positive Kaon From Primary Eta only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(20,"MC Positive Kaon From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(21,"MC Positive Kaon From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(22,"MC Positive Proton From Primary Eta only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(24,"MC Positive Proton From Primary Eta and Y");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(23,"MC Positive Proton From Primary Y Only");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(25,"MC Positive Pion From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(26,"MC Positive Kaon From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(27,"MC Positive Proton From Strangeness");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(28,"MC Positive Pion From Material");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(29,"MC Positive Kaon From Material");
    fHistSparseBefEvSel->GetAxis(2)->SetBinLabel(30,"MC Positive Proton From Material");
    fHistSparseBefEvSel->GetAxis(3)->SetTitle("Event selection");
    fHistSparseBefEvSel->GetAxis(3)->SetBinLabel(1,"Event Selection Not Passed - Sample Selection Not Passed");
    fHistSparseBefEvSel->GetAxis(3)->SetBinLabel(2,"Event Selection Not Passed - Sample Selection Passed");
    fHistSparseBefEvSel->GetAxis(3)->SetBinLabel(3,"Event Selection Passed - Sample Selection Not Passed");
    fHistSparseBefEvSel->GetAxis(3)->SetBinLabel(4,"Event Selection Passed - Sample Selection Passed");
    
    fHistSparseBefEvSel->Sumw2();
    if(fUseThnSparse) fOutput->Add(fHistSparseBefEvSel);
  }
  
  //Pt correlation histo
  Int_t ptbins[4] = {kNbins,kNbins,3,2};
  Double_t minptbins[4] = {0.08,0.08,0.,-1.};
  Double_t maxptbins[4] = {1.,1.,3.,1.};
  
  fHistPtCorr = new THnSparseD("fHistPtCorr","Pt correlation",4,ptbins,minptbins,maxptbins);
  fHistPtCorr->GetAxis(0)->Set(kNbins,fPtBinLimits);//Rebin Pt axis
  fHistPtCorr->GetAxis(1)->Set(kNbins,fPtBinLimits);//Rebin Pt axis
  fOutput->Add(fHistPtCorr);
  
  //Multiplicity Histos
  Int_t nv0multibin=102;
  Double_t v0bin[nv0multibin+1];
  Double_t x=0.;
  for(Int_t c=0; c<=nv0multibin; c++){
    v0bin[c]=x;
    if(c==0) x=0.01;
    else if(c==1) x=0.1;
    else if(c==2) x=1.;
    else if(c>2) x+=1.;
  }
  
  if(!fuseV0){
    fHistMult = new TH1F("fHistMult", "Event Multiplicity",3000,-0.5,2999.5);
    fHistMult->Sumw2();
    fHistMult->SetMinimum(0);
    if(fMultEstimator==0) fHistMult->GetXaxis()->SetTitle("Multiplicity |#eta|<0.8");
    else if(fMultEstimator==1) fHistMult->GetXaxis()->SetTitle("Tracklets |#eta|<0.8");
    else if(fMultEstimator==2) fHistMult->GetXaxis()->SetTitle("Clusters on SPD1");
    fOutput->Add(fHistMult);
    
    fHistMultAftEvSel = new TH1F("fHistMultAftEvSel", "Event Multiplicity afte event selection",3000,-0.5,2999.5);
    if(fMultEstimator==0) fHistMultAftEvSel->GetXaxis()->SetTitle("Multiplicity |#eta|<0.8");
    else if(fMultEstimator==1) fHistMultAftEvSel->GetXaxis()->SetTitle("Tracklets |#eta|<0.8");
    else if(fMultEstimator==2) fHistMultAftEvSel->GetXaxis()->SetTitle("Clusters on SPD1");
    fOutput->Add(fHistMultAftEvSel);
  }
  else{
    fHistMult = new TH1F("fHistMult", "Event V0 Multiplicity Percentile; %V0M; Counts",nv0multibin,v0bin);
    fHistMult->Sumw2();
    fHistMult->SetMinimum(0);
    fOutput->Add(fHistMult); 
    
    fHistMultAftEvSel = new TH1F("fHistMultAftEvSel", "Event V0 Multiplicity Percentile; %V0M; Counts",nv0multibin,v0bin);
    fHistMultAftEvSel->Sumw2();
    fHistMultAftEvSel->SetMinimum(0);
    fOutput->Add(fHistMultAftEvSel); 
  }
  
  UInt_t dimMultEst=3;//Reference Multiplicity-V0M Percentile-Event Selection
  if(fMC) dimMultEst=4;//Reference Multiplicity-V0M Percentile-Event Selection-Sample Selection
  Int_t binMultEst[dimMultEst];
  Double_t minMultEst[dimMultEst];
  Double_t maxMultEst[dimMultEst];
  binMultEst[0]=3000;//ReferenceMultiplicity
  minMultEst[0]=-0.5;
  maxMultEst[0]=2999.5;
  
  binMultEst[1]=nv0multibin;//V0MPercentile
  minMultEst[1]=0.;
  maxMultEst[1]=100;
  
  binMultEst[2]=2;//Event Selection
  minMultEst[2]=-.5;
  maxMultEst[2]=1.5; 
  
  if(fMC){
    binMultEst[3]=2;//Event Selection
    minMultEst[3]=-.5;
    maxMultEst[3]=1.5; 
  }
  
  fHistEventMultiplicity = new THnSparseD("fHistEventMultiplicity","Reference Multiplicity-V0M Percentile- EvSel-sampSel",dimMultEst,binMultEst,minMultEst,maxMultEst);
  fHistEventMultiplicity->GetAxis(0)->SetTitle("Reference Multiplicity");
  fHistEventMultiplicity->GetAxis(1)->SetTitle("V0M Percentile");
  fHistEventMultiplicity->GetAxis(1)->Set(nv0multibin,v0bin);//Rebin V0MPercentile axis
  fHistEventMultiplicity->GetAxis(2)->SetTitle("Event selection");
  if(fMC) fHistEventMultiplicity->GetAxis(3)->SetTitle("Sample selection");
  fOutput->Add(fHistEventMultiplicity);
  
  fHistBefPileUp = new TH1F("fHistBefPileUp","Difference between the primary vertex and secondaries vertices;Z_{prim}-Z_{sec};Counts",80000,-40,40);
  fOutput->Add(fHistBefPileUp);
  
  fHistTaggedPileUp = new TH1F("fHistTaggedPileUp","Difference between the primary vertex and secondaries vertices in case the event is tagged as pile up;Z_{prim}-Z_{sec};Counts",80000,-40,40);
  fOutput->Add(fHistTaggedPileUp);
  
  fHistCen = new TH1F("fHistCen", "Event Centrality",101,-0.5,100.5);
  fHistCen->Sumw2();
  fHistCen->SetMinimum(0);
  fOutput->Add(fHistCen);
  
  
  //Histo with track cuts
  fHistNTracks = new TH1F("fHistNTracks", "Number of ITSsa tracks",20,0.5,20.5);
  fHistNTracks->Sumw2();
  fHistNTracks->SetMinimum(0);
  
  TString label;
  
  label="no selection";//1
  fHistNTracks->GetXaxis()->SetBinLabel(kHasNoSelection,label.Data());
  label="ITSsa";//2
  fHistNTracks->GetXaxis()->SetBinLabel(kIsITSsa,label.Data());
  label="ITSrefit";//3
  fHistNTracks->GetXaxis()->SetBinLabel(kIsITSrefit,label.Data());
  label="neutral particle";//4
  fHistNTracks->GetXaxis()->SetBinLabel(kIsNotNeutralParticle,label.Data());
  label="OneSPDcl";//5
  fHistNTracks->GetXaxis()->SetBinLabel(kHasOneSPD,label.Data());
  label="TwoSPDcl";//6
  fHistNTracks->GetXaxis()->SetBinLabel(kHasTwoSPD,label.Data());
  label="SPDcls";//7
  fHistNTracks->GetXaxis()->SetBinLabel(kPassSPD,label.Data());
  label="3 SDD+SSD cls";//8
  fHistNTracks->GetXaxis()->SetBinLabel(kHas3PIDcls,label.Data());
  label="4 SDD+SSD cls";//9
  fHistNTracks->GetXaxis()->SetBinLabel(kHas4PIDcls,label.Data());
  label="SDD+SSD cls";//10
  fHistNTracks->GetXaxis()->SetBinLabel(kPassPIDcls,label.Data());
  label="chi2/ncls";//11
  fHistNTracks->GetXaxis()->SetBinLabel(kPassChi2Ncls,label.Data());
  label="eta";//12
  fHistNTracks->GetXaxis()->SetBinLabel(kIsInEta,label.Data());
  label="de/dx<0";//13
  fHistNTracks->GetXaxis()->SetBinLabel(kPassdEdx,label.Data());
  label="DCAz";//14
  fHistNTracks->GetXaxis()->SetBinLabel(kPassDCAzcut,label.Data());
  label="DCAxy";//15
  fHistNTracks->GetXaxis()->SetBinLabel(kPassDCAxycut,label.Data());
  
  fOutput->Add(fHistNTracks);
  
  //Histo with track cuts for positive particles
  fHistNTracksPos = new TH1F("fHistNTracksPos", "Number of positive ITSsa tracks",20,0.5,20.5);
  fHistNTracksPos->Sumw2();
  fHistNTracksPos->SetMinimum(0);
  
  label="no selection";//1
  fHistNTracksPos->GetXaxis()->SetBinLabel(kHasNoSelection,label.Data());
  label="ITSsa";//2
  fHistNTracksPos->GetXaxis()->SetBinLabel(kIsITSsa,label.Data());
  label="ITSrefit";//3
  fHistNTracksPos->GetXaxis()->SetBinLabel(kIsITSrefit,label.Data());
  label="neutral particle";//4
  fHistNTracksPos->GetXaxis()->SetBinLabel(kIsNotNeutralParticle,label.Data());
  label="OneSPDcl";//5
  fHistNTracksPos->GetXaxis()->SetBinLabel(kHasOneSPD,label.Data());
  label="TwoSPDcl";//6
  fHistNTracksPos->GetXaxis()->SetBinLabel(kHasTwoSPD,label.Data());
  label="SPDcls";//7
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassSPD,label.Data());
  label="3 SDD+SSD cls";//8
  fHistNTracksPos->GetXaxis()->SetBinLabel(kHas3PIDcls,label.Data());
  label="4 SDD+SSD cls";//9
  fHistNTracksPos->GetXaxis()->SetBinLabel(kHas4PIDcls,label.Data());
  label="SDD+SSD cls";//10
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassPIDcls,label.Data());
  label="chi2/ncls";//11
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassChi2Ncls,label.Data());
  label="eta";//12
  fHistNTracksPos->GetXaxis()->SetBinLabel(kIsInEta,label.Data());
  label="de/dx<0";//13
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassdEdx,label.Data());
  label="DCAz";//14
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassDCAzcut,label.Data());
  label="DCAxy";//15
  fHistNTracksPos->GetXaxis()->SetBinLabel(kPassDCAxycut,label.Data());
  
  fOutput->Add(fHistNTracksPos);
  
  //Histo with track cuts for negative particles
  fHistNTracksNeg = new TH1F("fHistNTracksNeg", "Number of negative ITSsa tracks",20,0.5,20.5);
  fHistNTracksNeg->Sumw2();
  fHistNTracksNeg->SetMinimum(0);
  
  label="no selection";//1
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kHasNoSelection,label.Data());
  label="ITSsa";//2
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kIsITSsa,label.Data());
  label="ITSrefit";//3
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kIsITSrefit,label.Data());
  label="neutral particle";//4
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kIsNotNeutralParticle,label.Data());
  label="OneSPDcl";//5
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kHasOneSPD,label.Data());
  label="TwoSPDcl";//6
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kHasTwoSPD,label.Data());
  label="SPDcls";//7
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassSPD,label.Data());
  label="3 SDD+SSD cls";//8
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kHas3PIDcls,label.Data());
  label="4 SDD+SSD cls";//9
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kHas4PIDcls,label.Data());
  label="SDD+SSD cls";//10
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassPIDcls,label.Data());
  label="chi2/ncls";//11
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassChi2Ncls,label.Data());
  label="eta";//12
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kIsInEta,label.Data());
  label="de/dx<0";//13
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassdEdx,label.Data());
  label="DCAz";//14
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassDCAzcut,label.Data());
  label="DCAxy";//15
  fHistNTracksNeg->GetXaxis()->SetBinLabel(kPassDCAxycut,label.Data());
  
  
  fOutput->Add(fHistNTracksNeg);
  
  //binning for the histogram
  const Int_t hnbins=400;
  Double_t hxmin = 0.01;
  Double_t hxmax = 10;
  Double_t hlogxmin = TMath::Log10(hxmin);
  Double_t hlogxmax = TMath::Log10(hxmax);
  Double_t hbinwidth = (hlogxmax-hlogxmin)/hnbins;
  Double_t hxbins[hnbins+1];
  hxbins[0] = 0.01; 
  for (Int_t i=1;i<=hnbins;i++) {
    hxbins[i] = hxmin + TMath::Power(10,hlogxmin+i*hbinwidth);
  }
  
  fHistDEDX = new TH2F("fHistDEDX","",hnbins,hxbins,900,0,1000);
  fOutput->Add(fHistDEDX);
  
  fHistDEDXdouble = new TH2F("fHistDEDXdouble","",500,-5,5,900,0,1000);
  fOutput->Add(fHistDEDXdouble);
  
  Int_t numbindedx[4];
  Double_t minbindedx[4];
  Double_t maxbindedx[4];
  numbindedx[0] = hnbins;//P
  minbindedx[0] = hxbins[0];
  maxbindedx[0] = hxbins[hnbins];
  
  numbindedx[1] = 900;//dE/dx
  minbindedx[1] = 0;
  maxbindedx[1] = 1000;
  
  if(fuseV0){
    numbindedx[2] = nv0multibin;//V0M Percentile
    minbindedx[2] = 0;
    maxbindedx[2] = 100;
  }
  else{
    numbindedx[2] = 3000;//Reference Multiplicity
    minbindedx[2] = -0.5;
    maxbindedx[2] = 2999.5;
  }
  
  numbindedx[3] = 3;//Particle
  minbindedx[3] = 0;
  maxbindedx[3] = 3;
  
  fHistDEDXMulti = new THnSparseD("fHistDEDXMulti","dE/dx vs p",4,numbindedx,minbindedx,maxbindedx);
  fHistDEDXMulti->GetAxis(0)->Set(hnbins,hxbins);
  fHistDEDXMulti->GetAxis(0)->SetTitle("P (GeV/c)");
  fHistDEDXMulti->GetAxis(1)->SetTitle("dE/dx");
  if(fuseV0){
    fHistDEDXMulti->GetAxis(2)->SetTitle("V0M %");
    fHistDEDXMulti->GetAxis(2)->Set(nv0multibin,v0bin);//Rebin V0MPercentile axis
  }
  else fHistDEDXMulti->GetAxis(2)->SetTitle("Nch");
  fHistDEDXMulti->GetAxis(3)->SetTitle("Particle");
  if(fLowMult==-1 && fUpMult ==-1) fOutput->Add(fHistDEDXMulti);
  
  numbindedx[0] = 500;//P
  minbindedx[0] = -5;
  maxbindedx[0] = 5;
  
  numbindedx[1] = 900;//dE/dx
  minbindedx[1] = 0;
  maxbindedx[1] = 1000;
  
  fHistDEDXdoubleMulti = new THnSparseD("fHistDEDXdoubleMulti","dE/dx vs signed p",4,numbindedx,minbindedx,maxbindedx);
  fHistDEDXdoubleMulti->GetAxis(0)->SetTitle("P (GeV/c)");
  fHistDEDXdoubleMulti->GetAxis(1)->SetTitle("dE/dx");
  fHistDEDXdoubleMulti->GetAxis(2)->SetTitle("Nch");
  fHistDEDXdoubleMulti->GetAxis(3)->SetTitle("Particle");
  if(fLowMult==-1 && fUpMult ==-1) fOutput->Add(fHistDEDXdoubleMulti);
  
  fHistBeforeEvSel = new TH1F("fHistBeforeEvSel","fHistBeforeEvSel",kNbins,fPtBinLimits);
  fOutput->Add(fHistBeforeEvSel);
  fHistAfterEvSel = new TH1F("fHistAfterEvSel","fHistAfterEvSel",kNbins,fPtBinLimits);
  fOutput->Add(fHistAfterEvSel);
  
  fHistITSchi2ncls = new TH1F("fHistITSchi2ncls","Chi2/ncls for the tracks",2000,0,20);  
  fOutput->Add(fHistITSchi2ncls);
  
  for(Int_t j=0;j<3;j++){
    
    fHistPosNSigmaSep[j] = new TH2F(Form("fHistPosNSigmaSep%d",j),"",hnbins,hxbins,1000,-10,10);
    fOutput->Add(fHistPosNSigmaSep[j]);
    fHistNegNSigmaSep[j] = new TH2F(Form("fHistNegNSigmaSep%d",j),"",hnbins,hxbins,1000,-10,10);
    fOutput->Add(fHistNegNSigmaSep[j]);
    
    fHistPrimMCpos[j] = new TH1F(Form("fHistPrimMCpos%d",j),Form("fHistPrimMCpos%d",j),kNbins,fPtBinLimits);
    fHistPrimMCneg[j] = new TH1F(Form("fHistPrimMCneg%d",j),Form("fHistPrimMCneg%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCpos[j]);
    fOutput->Add(fHistPrimMCneg[j]);
    fHistPrimMCposEtaY[j] = new TH1F(Form("fHistPrimMCposEtaY%d",j),Form("fHistPrimMCposEtaY%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegEtaY[j] = new TH1F(Form("fHistPrimMCnegEtaY%d",j),Form("fHistPrimMCnegEtaY%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCposEtaY[j]);
    fOutput->Add(fHistPrimMCnegEtaY[j]);
    fHistSecStrMCpos[j] = new TH1F(Form("fHistSecStrMCpos%d",j),Form("fHistSecStrMCpos%d",j),kNbins,fPtBinLimits);
    fHistSecStrMCneg[j] = new TH1F(Form("fHistSecStrMCneg%d",j),Form("fHistSecStrMCneg%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecStrMCneg[j]);
    fOutput->Add(fHistSecStrMCpos[j]);
    fHistSecMatMCpos[j] = new TH1F(Form("fHistSecMatMCpos%d",j),Form("fHistSecMatMCpos%d",j),kNbins,fPtBinLimits);
    fHistSecMatMCneg[j] = new TH1F(Form("fHistSecMatMCneg%d",j),Form("fHistSecMatMCneg%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecMatMCneg[j]);
    fOutput->Add(fHistSecMatMCpos[j]);
    //
    
    
    fHistPrimMCposBefEvSel[j] = new TH1F(Form("fHistPrimMCposBefEvSel%d",j),Form("fHistPrimMCposBefEvSel%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegBefEvSel[j] = new TH1F(Form("fHistPrimMCnegBefEvSel%d",j),Form("fHistPrimMCnegBefEvSel%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCnegBefEvSel[j]);
    fOutput->Add(fHistPrimMCposBefEvSel[j]);
    fHistPrimMCposBefEvSelEta[j] = new TH1F(Form("fHistPrimMCposBefEvSelEta%d",j),Form("fHistPrimMCposBefEvSelEta%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegBefEvSelEta[j] = new TH1F(Form("fHistPrimMCnegBefEvSelEta%d",j),Form("fHistPrimMCnegBefEvSelEta%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCposBefEvSelEta[j]);
    fOutput->Add(fHistPrimMCnegBefEvSelEta[j]);    
    fHistPrimMCposBefEvSelEtaY[j] = new TH1F(Form("fHistPrimMCposBefEvSelEtaY%d",j),Form("fHistPrimMCposBefEvSelEtaY%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegBefEvSelEtaY[j] = new TH1F(Form("fHistPrimMCnegBefEvSelEtaY%d",j),Form("fHistPrimMCnegBefEvSelEtaY%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCposBefEvSelEtaY[j]);
    fOutput->Add(fHistPrimMCnegBefEvSelEtaY[j]);
    fHistSecStrMCposBefEvSel[j] = new TH1F(Form("fHistSecStrMCposBefEvSel%d",j),Form("fHistSecStrMCposBefEvSel%d",j),kNbins,fPtBinLimits);
    fHistSecStrMCnegBefEvSel[j] = new TH1F(Form("fHistSecStrMCnegBefEvSel%d",j),Form("fHistSecStrMCnegBefEvSel%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecStrMCnegBefEvSel[j]);
    fOutput->Add(fHistSecStrMCposBefEvSel[j]);
    fHistSecMatMCposBefEvSel[j] = new TH1F(Form("fHistSecMatMCposBefEvSel%d",j),Form("fHistSecMatMCposBefEvSel%d",j),kNbins,fPtBinLimits);
    fHistSecMatMCnegBefEvSel[j] = new TH1F(Form("fHistSecMatMCnegBefEvSel%d",j),Form("fHistSecMatMCnegBefEvSel%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecMatMCnegBefEvSel[j]);
    fOutput->Add(fHistSecMatMCposBefEvSel[j]);
    //
    fHistPrimMCposReco[j] = new TH1F(Form("fHistPrimMCposReco%d",j),Form("fHistPrimMCposReco%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegReco[j] = new TH1F(Form("fHistPrimMCnegReco%d",j),Form("fHistPrimMCnegReco%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCnegReco[j]);
    fOutput->Add(fHistPrimMCposReco[j]);
    fHistPrimMCposRecoEtaY[j] = new TH1F(Form("fHistPrimMCposRecoEtaY%d",j),Form("fHistPrimMCposRecoEtaY%d",j),kNbins,fPtBinLimits);
    fHistPrimMCnegRecoEtaY[j] = new TH1F(Form("fHistPrimMCnegRecoEtaY%d",j),Form("fHistPrimMCnegRecoEtaY%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPrimMCposRecoEtaY[j]);
    fOutput->Add(fHistPrimMCnegRecoEtaY[j]);
    fHistSecStrMCposReco[j] = new TH1F(Form("fHistSecStrMCposReco%d",j),Form("fHistSecStrMCposReco%d",j),kNbins,fPtBinLimits);
    fHistSecStrMCnegReco[j] = new TH1F(Form("fHistSecStrMCnegReco%d",j),Form("fHistSecStrMCnegReco%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecStrMCnegReco[j]);
    fOutput->Add(fHistSecStrMCposReco[j]);
    fHistSecMatMCposReco[j] = new TH1F(Form("fHistSecMatMCposReco%d",j),Form("fHistSecMatMCposReco%d",j),kNbins,fPtBinLimits);
    fHistSecMatMCnegReco[j] = new TH1F(Form("fHistSecMatMCnegReco%d",j),Form("fHistSecMatMCnegReco%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistSecMatMCnegReco[j]);
    fOutput->Add(fHistSecMatMCposReco[j]);
    
  }
  
  for(Int_t i=0; i<4; i++){
    fHistCharge[i] = new TH1F(Form("fHistChargeLay%d",i),Form("fHistChargeLay%d",i),100,0,300);
    fOutput->Add(fHistCharge[i]);
  }
  
  for(Int_t i=0; i<kNbins; i++){
    fHistPosPi[i] = new TH1F(Form("fHistPosPi%d",i),Form("fHistPosPi%d",i),175,-3.5,3.5);	
    fHistPosK[i]  = new TH1F(Form("fHistPosK%d",i),Form("fHistPosK%d",i),175,-3.5,3.5);	
    fHistPosP[i]  = new TH1F(Form("fHistPosP%d",i),Form("fHistPosP%d",i),175,-3.5,3.5);	
    fHistNegPi[i] = new TH1F(Form("fHistNegPi%d",i),Form("fHistNegPi%d",i),175,-3.5,3.5);	
    fHistNegK[i]  = new TH1F(Form("fHistNegK%d",i),Form("fHistNegK%d",i),175,-3.5,3.5);	
    fHistNegP[i]  = new TH1F(Form("fHistNegP%d",i),Form("fHistNegP%d",i),175,-3.5,3.5);	
    
    
    fHistDCAPosPi[i] = new TH1F(Form("fHistDCAPosPi%d",i),Form("fHistDCAPosPi%d",i),2000,-1,1);  //DCA distr. with NSigma PID
    fHistDCAPosK[i]  = new TH1F(Form("fHistDCAPosK%d",i),Form("fHistDCAPosK%d",i),2000,-1,1);	
    fHistDCAPosP[i]  = new TH1F(Form("fHistDCAPosP%d",i),Form("fHistDCAPosP%d",i),2000,-1,1);	
    fHistDCANegPi[i] = new TH1F(Form("fHistDCANegPi%d",i),Form("fHistDCANegPi%d",i),2000,-1,1);	
    fHistDCANegK[i]  = new TH1F(Form("fHistDCANegK%d",i),Form("fHistDCANegK%d",i),2000,-1,1);	
    fHistDCANegP[i]  = new TH1F(Form("fHistDCANegP%d",i),Form("fHistDCANegP%d",i),2000,-1,1);	
    
    fHistMCPrimDCAPosPi[i] = new TH1F(Form("fHistMCPrimDCAPosPi%d",i),Form("fHistMCPrimDCAPosPi%d",i),2000,-1,1);  //DCA distr. with MC truth
    fHistMCPrimDCAPosK[i]  = new TH1F(Form("fHistMCPrimDCAPosK%d",i),Form("fHistMCPrimDCAPosK%d",i),2000,-1,1);	
    fHistMCPrimDCAPosP[i]  = new TH1F(Form("fHistMCPrimDCAPosP%d",i),Form("fHistMCPrimDCAPosP%d",i),2000,-1,1);	
    fHistMCPrimDCANegPi[i] = new TH1F(Form("fHistMCPrimDCANegPi%d",i),Form("fHistMCPrimDCANegPi%d",i),2000,-1,1);	
    fHistMCPrimDCANegK[i]  = new TH1F(Form("fHistMCPrimDCANegK%d",i),Form("fHistMCPrimDCANegK%d",i),2000,-1,1);	
    fHistMCPrimDCANegP[i]  = new TH1F(Form("fHistMCPrimDCANegP%d",i),Form("fHistMCPrimDCANegP%d",i),2000,-1,1);	
    
    fHistMCSecStDCAPosPi[i] = new TH1F(Form("fHistMCSecStDCAPosPi%d",i),Form("fHistMCSecStDCAPosPi%d",i),2000,-1,1);  //DCA distr. with MC truth
    fHistMCSecStDCAPosK[i]  = new TH1F(Form("fHistMCSecStDCAPosK%d",i),Form("fHistMCSecStDCAPosK%d",i),2000,-1,1);	
    fHistMCSecStDCAPosP[i]  = new TH1F(Form("fHistMCSecStDCAPosP%d",i),Form("fHistMCSecStDCAPosP%d",i),2000,-1,1);	
    fHistMCSecStDCANegPi[i] = new TH1F(Form("fHistMCSecStDCANegPi%d",i),Form("fHistMCSecStDCANegPi%d",i),2000,-1,1);	
    fHistMCSecStDCANegK[i]  = new TH1F(Form("fHistMCSecStDCANegK%d",i),Form("fHistMCSecStDCANegK%d",i),2000,-1,1);	
    fHistMCSecStDCANegP[i]  = new TH1F(Form("fHistMCSecStDCANegP%d",i),Form("fHistMCSecStDCANegP%d",i),2000,-1,1);	
    
    fHistMCSecMatDCAPosPi[i] = new TH1F(Form("fHistMCSecMatDCAPosPi%d",i),Form("fHistMCSecMatDCAPosPi%d",i),2000,-1,1);  //DCA distr. with MC truth
    fHistMCSecMatDCAPosK[i]  = new TH1F(Form("fHistMCSecMatDCAPosK%d",i),Form("fHistMCSecMatDCAPosK%d",i),2000,-1,1);	
    fHistMCSecMatDCAPosP[i]  = new TH1F(Form("fHistMCSecMatDCAPosP%d",i),Form("fHistMCSecMatDCAPosP%d",i),2000,-1,1);	
    fHistMCSecMatDCANegPi[i] = new TH1F(Form("fHistMCSecMatDCANegPi%d",i),Form("fHistMCSecMatDCANegPi%d",i),2000,-1,1);	
    fHistMCSecMatDCANegK[i]  = new TH1F(Form("fHistMCSecMatDCANegK%d",i),Form("fHistMCSecMatDCANegK%d",i),2000,-1,1);	
    fHistMCSecMatDCANegP[i]  = new TH1F(Form("fHistMCSecMatDCANegP%d",i),Form("fHistMCSecMatDCANegP%d",i),2000,-1,1);	
    
    fHistMCPosOtherHypPion[i] = new TH1F(Form("fHistMCPosOtherHypPion%d",i),Form("fHistMCPosOtherHypPion%d",i),175,-3.5,3.5);	//MC truth
    fHistMCPosOtherHypKaon[i] = new TH1F(Form("fHistMCPosOtherHypKaon%d",i),Form("fHistMCPosOtherHypKaon%d",i),175,-3.5,3.5);
    fHistMCPosOtherHypProton[i] = new TH1F(Form("fHistMCPosOtherHypProton%d",i),Form("fHistMCPosOtherHypProton%d",i),175,-3.5,3.5);
    fHistMCPosElHypPion[i] = new TH1F(Form("fHistMCPosElHypPion%d",i),Form("fHistMCPosElHypPion%d",i),175,-3.5,3.5);	
    fHistMCPosElHypKaon[i] = new TH1F(Form("fHistMCPosElHypKaon%d",i),Form("fHistMCPosElHypKaon%d",i),175,-3.5,3.5);
    fHistMCPosElHypProton[i] = new TH1F(Form("fHistMCPosElHypProton%d",i),Form("fHistMCPosElHypProton%d",i),175,-3.5,3.5);
    fHistMCPosPiHypPion[i] = new TH1F(Form("fHistMCPosPiHypPion%d",i),Form("fHistMCPosPiHypPion%d",i),175,-3.5,3.5);	
    fHistMCPosPiHypKaon[i] = new TH1F(Form("fHistMCPosPiHypKaon%d",i),Form("fHistMCPosPiHypKaon%d",i),175,-3.5,3.5);
    fHistMCPosPiHypProton[i] = new TH1F(Form("fHistMCPosPiHypProton%d",i),Form("fHistMCPosPiHypProton%d",i),175,-3.5,3.5);
    fHistMCPosKHypPion[i] = new TH1F(Form("fHistMCPosKHypPion%d",i),Form("fHistMCPosKHypPion%d",i),175,-3.5,3.5);	
    fHistMCPosKHypKaon[i]  = new TH1F(Form("fHistMCPosKHypKaon%d",i),Form("fHistMCPosKHypKaon%d",i),175,-3.5,3.5);	
    fHistMCPosKHypProton[i] = new TH1F(Form("fHistMCPosKHypProton%d",i),Form("fHistMCPosKHypProton%d",i),175,-3.5,3.5);
    fHistMCPosPHypPion[i] = new TH1F(Form("fHistMCPosPHypPion%d",i),Form("fHistMCPosPHypPion%d",i),175,-3.5,3.5);	
    fHistMCPosPHypKaon[i] = new TH1F(Form("fHistMCPosPHypKaon%d",i),Form("fHistMCPosPHypKaon%d",i),175,-3.5,3.5);
    fHistMCPosPHypProton[i]  = new TH1F(Form("fHistMCPosPHypProton%d",i),Form("fHistMCPosPHypProton%d",i),175,-3.5,3.5);	
    
    fHistMCNegOtherHypPion[i] = new TH1F(Form("fHistMCNegOtherHypPion%d",i),Form("fHistMCNegOtherHypPion%d",i),175,-3.5,3.5);	//MC truth
    fHistMCNegOtherHypKaon[i] = new TH1F(Form("fHistMCNegOtherHypKaon%d",i),Form("fHistMCNegOtherHypKaon%d",i),175,-3.5,3.5);
    fHistMCNegOtherHypProton[i] = new TH1F(Form("fHistMCNegOtherHypProton%d",i),Form("fHistMCNegOtherHypProton%d",i),175,-3.5,3.5);
    fHistMCNegElHypPion[i] = new TH1F(Form("fHistMCNegElHypPion%d",i),Form("fHistMCNegElHypPion%d",i),175,-3.5,3.5);
    fHistMCNegElHypKaon[i] = new TH1F(Form("fHistMCNegElHypKaon%d",i),Form("fHistMCNegElHypKaon%d",i),175,-3.5,3.5);
    fHistMCNegElHypProton[i] = new TH1F(Form("fHistMCNegElHypProton%d",i),Form("fHistMCNegElHypProton%d",i),175,-3.5,3.5);
    fHistMCNegPiHypPion[i] = new TH1F(Form("fHistMCNegPiHypPion%d",i),Form("fHistMCNegPiHypPion%d",i),175,-3.5,3.5);
    fHistMCNegPiHypKaon[i] = new TH1F(Form("fHistMCNegPiHypKaon%d",i),Form("fHistMCNegPiHypKaon%d",i),175,-3.5,3.5);
    fHistMCNegPiHypProton[i] = new TH1F(Form("fHistMCNegPiHypProton%d",i),Form("fHistMCNegPiHypProton%d",i),175,-3.5,3.5);
    fHistMCNegKHypPion[i] = new TH1F(Form("fHistMCNegKHypPion%d",i),Form("fHistMCNegKHypPion%d",i),175,-3.5,3.5);	
    fHistMCNegKHypKaon[i]  = new TH1F(Form("fHistMCNegKHypKaon%d",i),Form("fHistMCNegKHypKaon%d",i),175,-3.5,3.5);	
    fHistMCNegKHypProton[i] = new TH1F(Form("fHistMCNegKHypProton%d",i),Form("fHistMCNegKHypProton%d",i),175,-3.5,3.5);
    fHistMCNegPHypPion[i] = new TH1F(Form("fHistMCNegPHypPion%d",i),Form("fHistMCNegPHypPion%d",i),175,-3.5,3.5);	
    fHistMCNegPHypKaon[i] = new TH1F(Form("fHistMCNegPHypKaon%d",i),Form("fHistMCNegPHypKaon%d",i),175,-3.5,3.5);
    fHistMCNegPHypProton[i]  = new TH1F(Form("fHistMCNegPHypProton%d",i),Form("fHistMCNegPHypProton%d",i),175,-3.5,3.5);	
    
    
    fOutput->Add(fHistPosPi[i]);
    fOutput->Add(fHistPosK[i]);
    fOutput->Add(fHistPosP[i]);
    fOutput->Add(fHistNegPi[i]);
    fOutput->Add(fHistNegK[i]);
    fOutput->Add(fHistNegP[i]);
    
    
    fOutput->Add(fHistDCAPosPi[i]); //DCA distr
    fOutput->Add(fHistDCAPosK[i]);
    fOutput->Add(fHistDCAPosP[i]);
    fOutput->Add(fHistDCANegPi[i]);
    fOutput->Add(fHistDCANegK[i]);
    fOutput->Add(fHistDCANegP[i]);
    
    fOutput->Add(fHistMCPrimDCAPosPi[i]);//DCA distr.
    fOutput->Add(fHistMCPrimDCAPosK[i]);
    fOutput->Add(fHistMCPrimDCAPosP[i]);
    fOutput->Add(fHistMCPrimDCANegPi[i]);
    fOutput->Add(fHistMCPrimDCANegK[i]);
    fOutput->Add(fHistMCPrimDCANegP[i]);
    
    fOutput->Add(fHistMCSecStDCAPosPi[i]);//DCA distr.
    fOutput->Add(fHistMCSecStDCAPosK[i]);
    fOutput->Add(fHistMCSecStDCAPosP[i]);
    fOutput->Add(fHistMCSecStDCANegPi[i]);
    fOutput->Add(fHistMCSecStDCANegK[i]);
    fOutput->Add(fHistMCSecStDCANegP[i]);
    
    fOutput->Add(fHistMCSecMatDCAPosPi[i]);//DCA distr.
    fOutput->Add(fHistMCSecMatDCAPosK[i]);
    fOutput->Add(fHistMCSecMatDCAPosP[i]);
    fOutput->Add(fHistMCSecMatDCANegPi[i]);
    fOutput->Add(fHistMCSecMatDCANegK[i]);
    fOutput->Add(fHistMCSecMatDCANegP[i]);
    
    fOutput->Add(fHistMCPosOtherHypPion[i]);//MC truth
    fOutput->Add(fHistMCPosOtherHypKaon[i]);
    fOutput->Add(fHistMCPosOtherHypProton[i]);
    fOutput->Add(fHistMCPosElHypPion[i]);
    fOutput->Add(fHistMCPosElHypKaon[i]);
    fOutput->Add(fHistMCPosElHypProton[i]);
    fOutput->Add(fHistMCPosPiHypPion[i]);
    fOutput->Add(fHistMCPosPiHypKaon[i]);
    fOutput->Add(fHistMCPosPiHypProton[i]);
    fOutput->Add(fHistMCPosKHypPion[i]);
    fOutput->Add(fHistMCPosKHypKaon[i]);
    fOutput->Add(fHistMCPosKHypProton[i]);
    fOutput->Add(fHistMCPosPHypPion[i]);
    fOutput->Add(fHistMCPosPHypKaon[i]);
    fOutput->Add(fHistMCPosPHypProton[i]);
    
    fOutput->Add(fHistMCNegOtherHypPion[i]);//MC truth
    fOutput->Add(fHistMCNegOtherHypKaon[i]);
    fOutput->Add(fHistMCNegOtherHypProton[i]);
    fOutput->Add(fHistMCNegElHypPion[i]);
    fOutput->Add(fHistMCNegElHypKaon[i]);
    fOutput->Add(fHistMCNegElHypProton[i]);
    fOutput->Add(fHistMCNegPiHypPion[i]);
    fOutput->Add(fHistMCNegPiHypKaon[i]);
    fOutput->Add(fHistMCNegPiHypProton[i]);
    fOutput->Add(fHistMCNegKHypPion[i]);
    fOutput->Add(fHistMCNegKHypKaon[i]);
    fOutput->Add(fHistMCNegKHypProton[i]);
    fOutput->Add(fHistMCNegPHypPion[i]);
    fOutput->Add(fHistMCNegPHypKaon[i]);
    fOutput->Add(fHistMCNegPHypProton[i]);
    
  }
  
  //NSigma Histos
  for(Int_t j=0;j<3;j++){
    
    fHistPosNSigmaMean[j] = new TH1F(Form("hHistPosNSigmaMean%d",j),Form("hHistPosNSigmaMean%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaMean[j] = new TH1F(Form("hHistNegNSigmaMean%d",j),Form("hHistNegNSigmaMean%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaMCMean[j] = new TH1F(Form("hHistPosNSigmaMCMean%d",j),Form("hHistPosNSigmaMCMean%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaMCMean[j] = new TH1F(Form("hHistNegNSigmaMCMean%d",j),Form("hHistNegNSigmaMCMean%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrimMean[j] = new TH1F(Form("hHistPosNSigmaPrimMean%d",j),Form("hHistPosNSigmaPrimMean%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrimMean[j] = new TH1F(Form("hHistNegNSigmaPrimMean%d",j),Form("hHistNegNSigmaPrimMean%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrimMCMean[j] = new TH1F(Form("hHistPosNSigmaPrimMCMean%d",j),Form("hHistPosNSigmaPrimMCMean%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrimMCMean[j] = new TH1F(Form("hHistNegNSigmaPrimMCMean%d",j),Form("hHistNegNSigmaPrimMCMean%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPosNSigmaMean[j]);
    fOutput->Add(fHistNegNSigmaMean[j]);
    fOutput->Add(fHistPosNSigmaMCMean[j]);
    fOutput->Add(fHistNegNSigmaMCMean[j]);
    fOutput->Add(fHistPosNSigmaPrimMean[j]);
    fOutput->Add(fHistNegNSigmaPrimMean[j]);
    fOutput->Add(fHistPosNSigmaPrimMCMean[j]);
    fOutput->Add(fHistNegNSigmaPrimMCMean[j]);
    
    fHistPosNSigma[j] = new TH1F(Form("hHistPosNSigma%d",j),Form("hHistPosNSigma%d",j),kNbins,fPtBinLimits);
    fHistNegNSigma[j] = new TH1F(Form("hHistNegNSigma%d",j),Form("hHistNegNSigma%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaMC[j] = new TH1F(Form("hHistPosNSigmaMC%d",j),Form("hHistPosNSigmaMC%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaMC[j] = new TH1F(Form("hHistNegNSigmaMC%d",j),Form("hHistNegNSigmaMC%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrim[j] = new TH1F(Form("hHistPosNSigmaPrim%d",j),Form("hHistPosNSigmaPrim%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrim[j] = new TH1F(Form("hHistNegNSigmaPrim%d",j),Form("hHistNegNSigmaPrim%d",j),kNbins,fPtBinLimits);
    fHistPosNSigmaPrimMC[j] = new TH1F(Form("hHistPosNSigmaPrimMC%d",j),Form("hHistPosNSigmaPrimMC%d",j),kNbins,fPtBinLimits);
    fHistNegNSigmaPrimMC[j] = new TH1F(Form("hHistNegNSigmaPrimMC%d",j),Form("hHistNegNSigmaPrimMC%d",j),kNbins,fPtBinLimits);
    fOutput->Add(fHistPosNSigma[j]);
    fOutput->Add(fHistNegNSigma[j]);
    fOutput->Add(fHistPosNSigmaMC[j]);
    fOutput->Add(fHistNegNSigmaMC[j]);
    fOutput->Add(fHistPosNSigmaPrim[j]);
    fOutput->Add(fHistNegNSigmaPrim[j]);
    fOutput->Add(fHistPosNSigmaPrimMC[j]);
    fOutput->Add(fHistNegNSigmaPrimMC[j]);
    
  }
  
  if(fFillTree){
    OpenFile(2);
    fTreeNSigma = new TTree("fTreeNSigma","NSigma candidates");
    fTreeNSigma->Branch("fMult",&fMult,"fMult/I");
    fTreeNSigma->Branch("fV0Mult",&fV0Mult,"fV0Mult/F");
    fTreeNSigma->Branch("dedx",&fdedx,"fdedx/F");
    fTreeNSigma->Branch("p",&fP,"fP/F");
    fTreeNSigma->Branch("sign",&fTrkSign,"fTrkSign/B");
    fTreeNSigma->Branch("clumap",&fClumap,"fClumap/b");
    fTreeNSigma->Branch("eta",&fEta,"fEta/F");
    fTreeNSigma->Branch("evselmask",&fevSelmask,"fevSelmask/b");
    fTreeNSigma->Branch("impactxy",&fImpactXY,"fImpactXY/F");
    fTreeNSigma->Branch("impactz",&fImpactZ,"fImpactZ/F");
    fTreeNSigma->Branch("partmap",&fParticleMap,"fParticleMap/b");
    fTreeNSigma->Branch("ptmc",&fptMC,"fptMC/F");
    fTreeNSigma->Branch("chi2",&fChi2,"fChi2/F");
    
    fTreeNSigma->SetAutoSave(100000000);
    PostData(3,fTreeNSigma);
    
    if(fMC){
      OpenFile(3);
      fTreeMC = new TTree("fTreeMC","MC candidates");
      
      fTreeMC->Branch("fMult",&fMult,"fMult/I");
      fTreeMC->Branch("fV0Mult",&fV0Mult,"fV0Mult/F");
      fTreeMC->Branch("p",&fP,"fP/F");
      fTreeMC->Branch("sign",&fTrkSign,"fTrkSign/B");
      fTreeMC->Branch("eta",&fEta,"fEta/F");
      fTreeMC->Branch("evselmask",&fevSelmask,"fevSelmask/b");
      fTreeMC->Branch("partmap",&fParticleMap,"fParticleMap/b");
      fTreeMC->Branch("ptmc",&fptMC,"fptMC/F");
      
      fTreeMC->SetAutoSave(100000000);
      PostData(4,fTreeMC);
    }    
  } 
  else{
    fTreeNSigma = new TTree();
    fTreeMC = new TTree();
  }
  
  
  //////////////////////////////////
  
  // Post output data.
  PostData(1,fOutput);  
  
  
  infomsg("End of CreateOutputObjects");
}

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectraMultiplicity::UserExec(Option_t *){
  // Main loop
  // Called for each event
  ///////////////////////////////////////
  //variables
  Float_t pdgmass[4]={0.13957,0.493677,0.938272,1.8756}; //mass for pi, K, P, deuton (Gev/c^2)
  Int_t listcode[3]={211,321,2212};//code for pi, K, P (Gev/c^2)
  Double_t dedxLay[4];
  Float_t ptMC=-999;
  Int_t code=-999, signMC=-999,isph=-999,mfl=-999,uniqueID=-999;
  fImpactXY=-999;
  fImpactZ=-999;
  Int_t evSel=1;
  for(UInt_t c=0; c<kbitEv; c++) fevSelmask = fevSelmask & 0<<c;//Set all bit to 0
  Int_t sampSel=1;
  AliESDtrack *track;
  ULong_t status; 
  AliStack *stack=0;
  TParticle *part=0;
  TParticlePDG *pdgPart=0;
  ///////////////////////////////////////
  
  fESD=(AliESDEvent*)InputEvent(); //UserExec is called for every event fetched from the handler
  if(!fESD) {
    printf("AliAnalysisTaskSDDRP::Exec(): bad ESD\n");
    PostData(1,fOutput);
    PostData(2,fListCuts);
    return;
  } 
  fHistNEvents->Fill(-1);//Number of events opened -->Read from ESD
  if(fMC) fHistMCSampSel->Fill(-1);
  
  
  UInt_t maskPhysSel = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  TString firedTriggerClasses=fESD->GetFiredTriggerClasses();
  //  if(!firedTriggerClasses.Contains("CINT1B")) return;
  if((maskPhysSel & AliVEvent::kMB)==0){
    PostData(1,fOutput);
    PostData(2,fListCuts);
    return;
  }
  fHistNEvents->Fill(0); //-->Pass Phys. Sel. + Trig
  
  if(fLowEnergypp && !fMC){ // remove events without SDD in pp 2.76 TeV
    if(!firedTriggerClasses.Contains("ALL")){
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
  }
  fHistNEvents->Fill(1); //-->SDD read out
  
  if(fMC){//MC info and sample selection
    AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if (!eventHandler) {
      Printf("ERROR: Could not retrieve MC event handler");
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    AliMCEvent* mcEvent = eventHandler->MCEvent();
    if (!mcEvent) {
      Printf("ERROR: Could not retrieve MC event");
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    stack = mcEvent->Stack();
    if (!stack) {
      printf("ERROR: stack not available\n");
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    //Selection of the MC sample 
    const  AliVVertex *mcvtx = mcEvent->GetPrimaryVertex();
    if (!mcvtx) {
      printf("ERROR: mcvtx not available\n");
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    
    if(TMath::Abs(mcvtx->GetZ())>fvZRange) sampSel=0; //Position on Z of the vertex 
    else fevSelmask = fevSelmask | 1<<kVtxInRange;
    
  }
  
  //flags for MC
  Int_t nTrackMC=0; 
  if(stack) nTrackMC = stack->GetNtrack();
  
  //Reconstructed Vertex of the event
  const AliESDVertex *vtx =  fESD->GetPrimaryVertexSPD();
  
  ///////////selection of the centrality or multiplicity bin
  
  //selection on the event centrality
  if(fHImode){
    if(!(fLowCentrality<0.0)&&fUpCentrality>0.0)
    {
      AliCentrality *centrality = fESD->GetCentrality();
      if(!centrality->IsEventInCentralityClass(fLowCentrality,fUpCentrality,"V0M")){
	PostData(1,fOutput);
	PostData(2,fListCuts);
	return;
      }
      Printf("Centrality of the event: %.1f",centrality->GetCentralityPercentile("V0M"));
      Printf("Centrality cut: %.1f to %.1f",fLowCentrality,fUpCentrality);
      fHistCen->Fill(centrality->GetCentralityPercentile("V0M"));
    }
  }
  
  UInt_t dimsparse;//Definition of the THnSparse dimension and Initialization of the xsparse array
  if(!fMC) dimsparse=5;//Multiplicity (Nch or V0) - pt - PID asym - PID sym - DCA
  else dimsparse=7;//Multiplicity (Nch or V0) - pt - PID asym - PID sym - DCA - MCPID - MCpt
  Double_t xsparse[dimsparse];
  for(UInt_t counter = 0; counter < dimsparse ; counter++) xsparse[counter]=-999;//Default position for THnSparse 
  
  
  UInt_t dimsparseMC=0;//Definition of the THnSparse dimension and Initialization of the xsparse array for MC
  if(fMC) dimsparseMC=4;//Multiplicity (Nch or V0) - MCpt - MCPID - Event Selection
  Double_t xsparseMC[dimsparseMC];
  for(UInt_t counter = 0; counter < dimsparseMC ; counter++) xsparseMC[counter]=-999;//Default position for THnSparse 
  
  //selection on the event multiplicity based on global tracks  
  fMult = -1;  
  if(fMultEstimator==0){
    // tracks+tracklets
    //     fMult = fESDtrackCuts->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTrackletsITSTPC, 0.8); 
    fMult = AliPPVsMultUtils::GetStandardReferenceMultiplicity(fESD); 
  }else if(fMultEstimator==1){
    // tracklets
    fMult = fESDtrackCuts->GetReferenceMultiplicity(fESD, AliESDtrackCuts::kTracklets, 0.8);    
  }else if(fMultEstimator==2){
    // clusters in SPD1
    const AliMultiplicity *mult = fESD->GetMultiplicity();
    Float_t nClu1 = (Float_t)mult->GetNumberOfITSClusters(1);
    fMult = (Int_t)(AliESDUtils::GetCorrSPD2(nClu1,vtx->GetZ())+0.5);
  }
  
  //selection on the signal in the V0
  fV0Mult =-1;
  fV0Mult = fPPVsMultUtils->GetMultiplicityPercentile(fESD,"V0M",kFALSE);//No event selection here
  
  if(!fuseV0){//Use reference multiplicity
    if(fLowMult>-1 && fMult<fLowMult){//check on multiplicity lower bin
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    if(fUpMult>-1 && fMult>fUpMult){//check on multiplicity higher bin
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    fHistMult->Fill(fMult);
    fHistNEvents->Fill(2);//-->In mult. range
    if(fMC) fHistMCSampSel->Fill(0);
    xsparse[0]=fMult;
    if(fMC) xsparseMC[0]=fMult;
  }
  else{//Use V0M estimator
    if(fLowMult>-1 && fV0Mult<fLowMult){//check on multiplicity lower bin
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    if(fUpMult>-1 && fV0Mult>=fUpMult){//check on multiplicity higher bin
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    fHistMult->Fill(fV0Mult);
    fHistNEvents->Fill(2);//-->In mult. range
    if(fMC) fHistMCSampSel->Fill(0);
    xsparse[0]=fV0Mult;
    if(fMC) xsparseMC[0]=fV0Mult;
  }
  
  //sample selection
  if(fMC && sampSel==1) fHistMCSampSel->Fill(1);
  
  if(!vtx->GetStatus() && AliPPVsMultUtils::IsINELgtZERO( fESD )) errormsg("UserExec: Check on the consistency in the event selection returned a bad value!");
  
  //------------------------------------------------
  // Selection Investigation with AliPPVsMultUtils
  //------------------------------------------------
  
  //------------------------------------------------
  //Check if event is selected 
  //------------------------------------------------
  if( AliPPVsMultUtils::IsEventSelected (fESD) ) fHistNEvents->Fill(9);
  else evSel = 0;
  
  //------------------------------------------------
  //Step 1: Check for Min-Bias Trigger
  //------------------------------------------------
  if( AliPPVsMultUtils::IsMinimumBias( fESD ) ){
    fHistNEvents->Fill(3);
    //------------------------------------------------
    //Step 2: Check for INEL>0
    //------------------------------------------------
    if( vtx->GetStatus()){
      fHistNEvents->Fill(4);
      if( AliPPVsMultUtils::IsAcceptedVertexPosition( fESD ) ){
	fHistNEvents->Fill(5);
	//------------------------------------------------
	//Step 3: Check for Vertex-Z position
	//------------------------------------------------
	if( AliPPVsMultUtils::IsINELgtZERO( fESD ) ){
	  fHistNEvents->Fill(6);
	  //------------------------------------------------
	  //Step 4: Check for SPD Pileup
	  //------------------------------------------------
	  if( AliPPVsMultUtils::IsNotPileupSPDInMultBins( fESD ) ){
	    fHistNEvents->Fill(7);
	    //------------------------------------------------
	    //Step 5: Check for SPD / track vertex consistency
	    //------------------------------------------------
	    if( AliPPVsMultUtils::HasNoInconsistentSPDandTrackVertices( fESD ) ){
	      fHistNEvents->Fill(8);
	      
	    }
	  }
	}
      }
    }
  }
  
  
  xsparseMC[3]=2*evSel+sampSel;
  
  if(evSel==1){
    if(!fuseV0) fHistMultAftEvSel->Fill(fMult);
    else fHistMultAftEvSel->Fill(fV0Mult);
  }
  
  UInt_t dimsparseMult = 3;
  if(fMC) dimsparseMult = 4;
  Double_t multicontainer[dimsparseMult];
  multicontainer[0] = fMult;
  multicontainer[1] = fV0Mult;
  multicontainer[2] = evSel;
  if(fMC) multicontainer[3] = sampSel;
  
  if(fLowMult ==-1 && fUpMult ==-1)fHistEventMultiplicity->Fill(multicontainer,1);
  
  //first loop on stack, before event selection, filling MC ntuple
  for(Int_t imc=0; imc<nTrackMC; imc++){
    if(imc == 1) fevSelmask = fevSelmask & 0<<kIsNewEvent;
    
    for(UInt_t counter = 1; counter < dimsparseMC-1; counter++) xsparseMC[counter]=-999;//Default position for THnSparse Multiplicity and event selection are defined for the whole event
    
    part = stack->Particle(imc);
    isph=1;    
    if(!stack->IsPhysicalPrimary(imc)) isph=0;
    pdgPart = part->GetPDG();
    if(!pdgPart)continue;
    if(pdgPart->Charge()==0) continue; //no neutral particles
    Float_t yMC=-999.;
    if(part->Energy() != TMath::Abs(part->Pz())) yMC = 0.5*TMath::Log((part->Energy()+part->Pz())/(part->Energy()-part->Pz()));
    if(pdgPart->Charge()>0) signMC=1;
    else signMC=-1;
    ptMC=part->Pt();
    xsparseMC[1]=ptMC;
    code=pdgPart->PdgCode();
    Int_t jpart=-1;
    for(Int_t j=0; j<3; j++){
      if(TMath::Abs(code)==listcode[j]){
	jpart=j;
	break;
      }
    }
    
    if(jpart>=0 && isph && sampSel && TMath::Abs(part->Eta())< fEtaRange){
      if(signMC>0){
	fHistPrimMCposBefEvSelEta[jpart]->Fill(ptMC);
	xsparseMC[2]=jpart*3+0.5;
      }
      else{
	fHistPrimMCnegBefEvSelEta[jpart]->Fill(ptMC);
	xsparseMC[2]=-jpart*3-0.5;	    
      }
    }
    
    //     if(TMath::Abs(yMC) > fMaxY) continue; //rapidity cut
    if(TMath::Abs(yMC) <= fMaxY){ //rapidity cut
      
      Int_t indexMoth=part->GetFirstMother();//Information of the provenance of the particle
      if(indexMoth>=0){
	TParticle* moth = stack->Particle(indexMoth);
	Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
      }
      uniqueID = part->GetUniqueID();
      
      //filling MC Ntuple
      if((TMath::Abs(code)==211 || TMath::Abs(code)==321 || TMath::Abs(code)==2212) && fFillTree){
	fptMC=ptMC;
	fTreeMC->Fill();
      }
      
      if(jpart>=0){
	if(isph){//MC Primaries 
	  if(sampSel==1){//Selection of the MC sample before event selection
	    if(signMC>0){
	      fHistPrimMCposBefEvSel[jpart]->Fill(ptMC);
	      if(TMath::Abs(part->Eta())<fEtaRange) fHistPrimMCposBefEvSelEtaY[jpart]->Fill(ptMC);
	    }
	    else{
	      fHistPrimMCnegBefEvSel[jpart]->Fill(ptMC);
	      if(TMath::Abs(part->Eta())<fEtaRange) fHistPrimMCnegBefEvSelEtaY[jpart]->Fill(ptMC);
	    }
	  }
	  if(evSel==1){//Event selection 
	    if(signMC>0){
	      fHistPrimMCpos[jpart]->Fill(ptMC);
	      if(TMath::Abs(part->Eta())<fEtaRange) fHistPrimMCposEtaY[jpart]->Fill(ptMC);
	    }
	    else{
	      fHistPrimMCneg[jpart]->Fill(ptMC);
	      if(TMath::Abs(part->Eta())<fEtaRange) fHistPrimMCnegEtaY[jpart]->Fill(ptMC);  
	    }
	  }
	  
	  //With THnSparse Information of the sample and event selection is already in memory!
	  if(signMC>0){
	    if(xsparseMC[2]>-999) xsparseMC[2]= jpart*3+1.5;
	    else xsparseMC[2]= jpart*3+2.5;
	  }
	  else{
	    if(xsparseMC[2]>-999) xsparseMC[2]= -jpart*3-1.5;
	    else xsparseMC[2]= -jpart*3-2.5;
	  }  
	}
	else{//MC Secondaries
	  if(mfl==3 && uniqueID == kPDecay){ // If a particle is not a physical primary, check if it comes from weak decay
	    if(sampSel==1){//Selection of the MC sample before event selection
	      if(signMC>0) fHistSecStrMCposBefEvSel[jpart]->Fill(ptMC);
	      else fHistSecStrMCnegBefEvSel[jpart]->Fill(ptMC);		
	    }
	    if(evSel==1){//Event selection
	      if(signMC>0) fHistSecStrMCpos[jpart]->Fill(ptMC);
	      else  fHistSecStrMCneg[jpart]->Fill(ptMC);	    
	    }
	    //With THnSparse Information of the sample and event selection is already in memory!
	    if(signMC>0)xsparseMC[2]=jpart+9.5;
	    else xsparseMC[2]=-jpart-9.5;
	  }
	  else{//From material
	    if(sampSel==1){//Selection of the MC sample before event selection
	      if(signMC>0) fHistSecMatMCposBefEvSel[jpart]->Fill(ptMC);
	      else fHistSecMatMCnegBefEvSel[jpart]->Fill(ptMC);	      
	    }
	    if(evSel==1){//Event selection
	      if(signMC>0) fHistSecMatMCpos[jpart]->Fill(ptMC);
	      else  fHistSecMatMCneg[jpart]->Fill(ptMC);	    
	    }
	    //With THnSparse Information of the sample and event selection is already in memory!
	    if(signMC>0) xsparseMC[2]=jpart+12.5;
	    else xsparseMC[2]=-jpart-12.5;
	  }
	}
      }
    }
    //Filling Time!
    if(xsparseMC[0]>-999 && xsparseMC[1]>-999 && xsparseMC[2]>-999 && xsparseMC[3]>-999 && fUseThnSparse) fHistSparseBefEvSel->Fill(xsparseMC,1);
  }//Loop on MC stack
  
  //event selection
  if(evSel==0 && fFillTree==kFALSE){
    PostData(1,fOutput);
    PostData(2,fListCuts);
    return;  
  }
  
  //loop on tracks
  for (Int_t iTrack=0; iTrack<fESD->GetNumberOfTracks(); iTrack++) {  
    if(iTrack==0) fevSelmask = fevSelmask | 1<<kIsNewEvent;
    else if(iTrack ==1) fevSelmask = fevSelmask & 0<<kIsNewEvent;
    
    
    isph=-999;
    code=-999;
    mfl=-999;
    uniqueID=-999;
    
    for(UInt_t counter = 1; counter < dimsparse ; counter++){//Reset of array variables, Multiplicity is set for the entire event
      xsparse[counter]=-999;//Default position for THnSparse 
    } 
    
    track = (AliESDtrack*)fESD->GetTrack(iTrack);      
    if (!track) continue;
    
    //track selection
    if(track->GetSign()>0) fTrkSign = 1;
    else if(track->GetSign()) fTrkSign = -1;
    
    TString label;
    
    label="no selection";
    fHistNTracks->Fill(kHasNoSelection);
    if(fTrkSign>0)fHistNTracksPos->Fill(kHasNoSelection);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kHasNoSelection);
    
    status=track->GetStatus();
    if((status&AliESDtrack::kITSpureSA)==0) continue; //its standalone
    
    label="ITSsa";
    fHistNTracks->Fill(kIsITSsa);
    if(fTrkSign>0)fHistNTracksPos->Fill(kIsITSsa);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kIsITSsa);
    
    if((status&AliESDtrack::kITSrefit)==0) continue; //its refit
    
    label="ITSrefit";
    fHistNTracks->Fill(kIsITSrefit);
    if(fTrkSign>0)fHistNTracksPos->Fill(kIsITSrefit);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kIsITSrefit);
    
    if(TMath::Abs(track->GetSign())<0.0001) continue; //no neutral particles
    
    label="neutral particle";
    fHistNTracks->Fill(kIsNotNeutralParticle);
    if(fTrkSign>0)fHistNTracksPos->Fill(kIsNotNeutralParticle);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kIsNotNeutralParticle);
    
    //cluster in ITS
    fClumap = track->GetITSClusterMap();
    Int_t nSPD=0;
    for(Int_t il=0; il<2; il++) if(TESTBIT(fClumap,il)) nSPD++;
    
    if(nSPD == 1){
      label="OneSPDcl";
      fHistNTracks->Fill(kHasOneSPD);
      if(fTrkSign>0)fHistNTracksPos->Fill(kHasOneSPD);
      if(fTrkSign<0)fHistNTracksNeg->Fill(kHasOneSPD);
    }
    
    if(nSPD == 2){
      label="TwoSPDcl";
      fHistNTracks->Fill(kHasTwoSPD);
      if(fTrkSign>0)fHistNTracksPos->Fill(kHasTwoSPD);
      if(fTrkSign<0)fHistNTracksNeg->Fill(kHasTwoSPD);
    }
    
    if(nSPD<fMinSPDPts) continue;//At least one point in the SPD
    
    label="SPDcls";
    fHistNTracks->Fill(kPassSPD);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassSPD);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassSPD);
    
    Int_t nPtsForPid=0;
    for(Int_t j=2;j<6;j++) if(TESTBIT(fClumap,j)) nPtsForPid++;
    
    if(nPtsForPid == 3){
      label="3 SDD+SSD cls";
      fHistNTracks->Fill(kHas3PIDcls);
      if(fTrkSign>0)fHistNTracksPos->Fill(kHas3PIDcls);
      if(fTrkSign<0)fHistNTracksNeg->Fill(kHas3PIDcls); 
    }
    else if(nPtsForPid == 4){
      label="4 SDD+SSD cls";
      fHistNTracks->Fill(kHas4PIDcls);
      if(fTrkSign>0)fHistNTracksPos->Fill(kHas4PIDcls);
      if(fTrkSign<0)fHistNTracksNeg->Fill(kHas4PIDcls);
    }
    
    if(nPtsForPid<fMinNdEdxSamples) continue; //at least 3 points on SSD/SDD
    
    label="SDD+SSD cls";
    fHistNTracks->Fill(kPassPIDcls);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassPIDcls);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassPIDcls);
    
    //chisquare/nclusters	
    Int_t nclu=nSPD+nPtsForPid;
    fChi2 = track->GetITSchi2();
    Double_t chi2ncls = fChi2/nclu;
    fHistITSchi2ncls->Fill(chi2ncls);
    if(chi2ncls > fMaxChi2Clu) continue; 
    
    label="chi2/ncls";
    fHistNTracks->Fill(kPassChi2Ncls);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassChi2Ncls);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassChi2Ncls);
    
    //pseudorapidity
    fEta = track->Eta();
    if(TMath::Abs(fEta) > fEtaRange) continue;
    
    label="eta";
    fHistNTracks->Fill(kIsInEta);
    if(fTrkSign>0)fHistNTracksPos->Fill(kIsInEta);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kIsInEta);
    
    //truncated mean
    track->GetITSdEdxSamples(dedxLay);
    fdedx = CookdEdx(dedxLay);
    if(fdedx<0) continue;
    
    label="de/dx<0";
    fHistNTracks->Fill(kPassdEdx);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassdEdx);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassdEdx);
    
    Float_t pt = track->Pt();
    xsparse[1]=TMath::Abs(pt);
    Int_t theBin=-1;
    for(Int_t m=0; m<kNbins; m++){
      if(TMath::Abs(pt) >= fPtBinLimits[m] && TMath::Abs(pt) < fPtBinLimits[m+1]){//Choice of the pt bin 
	theBin=m;
	break;
      }
    }
    track->GetImpactParameters(fImpactXY, fImpactZ);
    xsparse[4]=fImpactXY;
    
    //DCAz Cut
    //Because we fit the DCA distributions, we need them after the DCAz cut, otherwise the templates and distributions are modified by this cut (if it is done with the one of the DCAxy)
    if(!DCAcutZ(fImpactZ,pt)) continue;
    label="DCAz";   
    
    fHistNTracks->Fill(kPassDCAzcut);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassDCAzcut);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassDCAzcut);
    
    //DCA Templates
    //information from the MC kinematics
    if(fMC){
      Int_t trklabel = TMath::Abs(track->GetLabel());
      if(trklabel<0)isph=-1;
      if(trklabel>=0){
	part = (TParticle*)stack->Particle(trklabel);
	pdgPart = part->GetPDG();
	code = pdgPart->PdgCode();
	ptMC = part->Pt();
	if(pdgPart->Charge()>0) signMC=1;
	else signMC=-1;
	if(stack->IsPhysicalPrimary(trklabel)) isph=1;
	else{//If it is not from php check the origin
	  isph=0;
	  Int_t indexMoth=part->GetFirstMother();
	  if(indexMoth>=0){
	    TParticle* moth = stack->Particle(indexMoth);
	    Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	    mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	  }
	  uniqueID = part->GetUniqueID();
	}
	
	//Filling DCA distribution with MC truth
	if(theBin>=0 && theBin<kNbins){
	  if(isph==1){//primaries in MC
	    if(fTrkSign>0){//Positive from Primary
	      if(TMath::Abs(code)==listcode[0]){ 
		fHistMCPrimDCAPosPi[theBin]->Fill(fImpactXY);
	      }
	      if(TMath::Abs(code)==listcode[1]){
		fHistMCPrimDCAPosK[theBin]->Fill(fImpactXY);
	      }
	      if(TMath::Abs(code)==listcode[2]){
		fHistMCPrimDCAPosP[theBin]->Fill(fImpactXY);
	      }
	    }else{//Negative from primary
	      if(TMath::Abs(code)==listcode[0]){
		fHistMCPrimDCANegPi[theBin]->Fill(fImpactXY);
	      }
	      if(TMath::Abs(code)==listcode[1]){
		fHistMCPrimDCANegK[theBin]->Fill(fImpactXY);
	      }
	      if(TMath::Abs(code)==listcode[2]){
		fHistMCPrimDCANegP[theBin]->Fill(fImpactXY);
	      }
	    }
	  }  
	  else if(isph==0){//Secondaries in MC
	    if(mfl==3 && uniqueID == kPDecay){//Secondaries from strangeness
	      if(fTrkSign>0){//Positive from strangeness
		if(TMath::Abs(code)==listcode[0]) {
		  fHistMCSecStDCAPosPi[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[1]) {
		  fHistMCSecStDCAPosK[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[2]) {
		  fHistMCSecStDCAPosP[theBin]->Fill(fImpactXY);
		}
	      }else{//Negative from strangeness
		if(TMath::Abs(code)==listcode[0]) {
		  fHistMCSecStDCANegPi[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[1]) {
		  fHistMCSecStDCANegK[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[2]) {
		  fHistMCSecStDCANegP[theBin]->Fill(fImpactXY);
		}
	      }
	    }else{//Secondaries from material
	      if(fTrkSign>0){//Positive from material
		if(TMath::Abs(code)==listcode[0]) {
		  fHistMCSecMatDCAPosPi[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[1]) {
		  fHistMCSecMatDCAPosK[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[2]) {
		  fHistMCSecMatDCAPosP[theBin]->Fill(fImpactXY);
		}
	      }else{//Negative from material
		if(TMath::Abs(code)==listcode[0]) {
		  fHistMCSecMatDCANegPi[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[1]) {
		  fHistMCSecMatDCANegK[theBin]->Fill(fImpactXY);
		}
		if(TMath::Abs(code)==listcode[2]) {
		  fHistMCSecMatDCANegP[theBin]->Fill(fImpactXY);
		}
	      }
	    }
	  }
	}
	
	if(1){//with THnSparse
	  if(isph==1){//primaries in MC
	    xsparse[5]=9.5;//Primaries without selection
	    if(fTrkSign>0){//Positive from Primary
	      if(TMath::Abs(code)==listcode[0]){ 
		xsparse[5]=0.5;//Positive Pion Primary Hipotesis
	      }
	      if(TMath::Abs(code)==listcode[1]){
		xsparse[5]=1.5;//Positive Kaon Primary Hipotesis
	      }
	      if(TMath::Abs(code)==listcode[2]){
		xsparse[5]=2.5;//Positive Proton Primary Hipotesis
	      }
	    }else{//Negative from primary
	      if(TMath::Abs(code)==listcode[0]){
		xsparse[5]=-0.5;//Negative Pion Primary Hipotesis
	      }
	      if(TMath::Abs(code)==listcode[1]){
		xsparse[5]=-1.5;//Negative Kaon Primary Hipotesis
	      }
	      if(TMath::Abs(code)==listcode[2]){
		xsparse[5]=-2.5;//Negative Proton Primary Hipotesis
	      }
	    }
	  }
	  if(isph==0){//Secondaries in MC
	    if(mfl==3 && uniqueID == kPDecay){//Secondaries from strangeness
	      xsparse[5]=10.5;//Secondaries from strangeness without selection
	      if(fTrkSign>0){//Positive from strangeness
		if(TMath::Abs(code)==listcode[0]) {
		  xsparse[5]=0.5+3;//Positive Pion Secondary from strangeness Hipotesis
		}
		if(TMath::Abs(code)==listcode[1]) {
		  xsparse[5]=1.5+3;//Positive Kaon Secondary from strangeness Hipotesis
		}
		if(TMath::Abs(code)==listcode[2]) {
		  xsparse[5]=2.5+3;//Positive Proton Secondary from strangeness Hipotesis
		}
	      }else{//Negative from strangeness
		if(TMath::Abs(code)==listcode[0]) {
		  xsparse[5]=-0.5-3;//Negative Pion Secondary from strangeness Hipotesis
		}
		if(TMath::Abs(code)==listcode[1]) {
		  xsparse[5]=-1.5-3;//Negative Kaon Secondary from strangeness Hipotesis
		}
		if(TMath::Abs(code)==listcode[2]) {
		  xsparse[5]=-2.5-3;//Negative Proton Secondary from strangeness Hipotesis
		}
	      }
	    }else{//Secondaries from material
	      xsparse[5]=11.5;//Secondaries from material without selection	      
	      if(fTrkSign>0){//Positive from material
		if(TMath::Abs(code)==listcode[0]) {
		  xsparse[5]=0.5+6;//Positive Pion Secondary from material Hipotesis
		}
		if(TMath::Abs(code)==listcode[1]) {
		  xsparse[5]=1.5+6;//Positive Pion Secondary from material Hipotesis
		}
		if(TMath::Abs(code)==listcode[2]) {
		  xsparse[5]=2.5+6;//Positive Pion Secondary from material Hipotesis
		}
	      }else{//Negative from material
		if(TMath::Abs(code)==listcode[0]) {
		  xsparse[5]=-0.5-6;//Positive Pion Secondary from material Hipotesis
		}
		if(TMath::Abs(code)==listcode[1]) {
		  xsparse[5]=-1.5-6;//Positive Pion Secondary from material Hipotesis
		}
		if(TMath::Abs(code)==listcode[2]) {
		  xsparse[5]=-2.5-6;//Positive Pion Secondary from material Hipotesis
		}
	      }
	    }
	  }
	}
	
	for(UInt_t c=0; c<8; c++) fParticleMap = fParticleMap& 0<<c;//Set all bit to 0
	if(signMC>0) fParticleMap = fParticleMap | 1<<kIsPositive;
	if(TMath::Abs(code)==listcode[0]) fParticleMap = fParticleMap | 1<<kIsPion;
	if(TMath::Abs(code)==listcode[1]) fParticleMap = fParticleMap | 1<<kIsKaon;
	if(TMath::Abs(code)==listcode[2]) fParticleMap = fParticleMap | 1<<kIsProton;
	if(isph) fParticleMap = fParticleMap | 1<<kIsPhysicalPrimary;
	else if(mfl==3 && uniqueID == kPDecay) fParticleMap = fParticleMap | 1<<kIsFromStrangeness;
	else fParticleMap = fParticleMap | 1<<kIsFromMaterial;
	
      }
    }
    
    //Filling Ntuple
    if(fFillTree){
      fP=track->GetP();
      
      if(fMC) fptMC = ptMC;
      else fptMC = -999;
      fTreeNSigma->Fill();
      
      PostData(3,fTreeNSigma);
      if(fMC) PostData(4,fTreeMC);    
      
      PostData(1,fOutput);
      PostData(2,fListCuts);
      return;
    }
    
    
    //Compute y and bb
    Double_t y[4],bbtheo[4],logdiff[4];
    Float_t p=track->GetP();
    if(fMC && fSmearMC){
      fdedx=fRandGener->Gaus(fdedx,fSmeardEdx*fdedx);
      p=fRandGener->Gaus(p,fSmearP*p);     
    }
    
    //Nsigma Method
    Float_t resodedx[4];
    for(Int_t ires=0;ires<4;ires++){//calculation of resolution on dedx
      resodedx[ires]=fITSPIDResponse->GetResolution(1,ires+1,kTRUE);
    }
    
    for(Int_t i=0;i<4;i++){//calculation of y and bb and log difference
      y[i] = Eta2y(pt,pdgmass[i],track->Eta());
      //bbtheo[i]=fITSPIDResponse->Bethe(p,pdgmass[i],kTRUE); //Pure PHOBOS BB
      bbtheo[i]=fITSPIDResponse->BetheITSsaHybrid(p,pdgmass[i]); //PHOBOS + polinomial at low pt (below beta*gamma = 0.76)
      logdiff[i]=TMath::Log(fdedx) - TMath::Log(bbtheo[i]);
    }
    
    Int_t resocls=(Int_t)nPtsForPid-1;
    
    //NSigma Method, with asymmetric bands->Mean is for asymmetric
    Int_t minPosMean=-1;
    for(Int_t isp=0; isp<3; isp++){
      if(fdedx<bbtheo[0])continue;
      Double_t bb=bbtheo[isp];
      if(fdedx<bb){
	Double_t bbdistance=TMath::Abs((bbtheo[isp]-bbtheo[isp-1])/2);
	if(bbdistance==0) continue;
	Double_t nsigma=TMath::Abs((fdedx-bb)/bbdistance);
	if(nsigma<1.)minPosMean=isp;
      }
      else{
	Double_t bbdistance=TMath::Abs((bbtheo[isp]-bbtheo[isp+1])/2);
	if(bbdistance==0) continue;
	Double_t nsigma=TMath::Abs((fdedx-bb)/bbdistance);
	if(nsigma<1.)minPosMean=isp;
      }
    }
    if(fdedx<bbtheo[0] && TMath::Abs((fdedx-bbtheo[0])/(resodedx[resocls]*bbtheo[0]))<2)minPosMean=0;
    
    //NSigma method with simmetric bands-> no Mean is for symmetric 
    Double_t nsigmas[3];
    Double_t min=999999.;
    Int_t minPos=-1;
    for(Int_t isp=0; isp<3; isp++){
      Double_t bb=bbtheo[isp];
      nsigmas[isp]=TMath::Abs((fdedx-bb)/(resodedx[resocls]*bb));
      if(nsigmas[isp]<min){
	min=nsigmas[isp];
	minPos=isp;
      }
      //Filling histos with nsigma separation
      if(fTrkSign>0)fHistPosNSigmaSep[isp]->Fill(track->GetP(),((fdedx-bb)/(resodedx[resocls]*bb)));
      else fHistNegNSigmaSep[isp]->Fill(track->GetP(),((fdedx-bb)/(resodedx[resocls]*bb)));
    }
    
    // y calculation
    Double_t yPartMean=-999;
    Double_t yPart=-999;
    if(minPosMean >-1)yPartMean =y[minPosMean];
    if(minPos >-1)yPart=y[minPos];
    

    if(TMath::Abs(yPartMean)<fMaxY){//DCA distributions, before the DCAxy cut, based on asymmetric nsigma approach
      if(xsparse[1]>=0){
	if(fTrkSign>0){//with THnSparse
	  if(minPosMean==0){ 
	    xsparse[2]=0.5;
	  }
	  else if(minPosMean==1){
	    xsparse[2]=2.5;
	  }
	  else if(minPosMean==2){
	    xsparse[2]=4.5;
	  }
	}
	else{
	  if(minPosMean==0){
	    xsparse[2]=-0.5;
	  }else if(minPosMean==1){
	    xsparse[2]=-2.5;
	  }else if(minPosMean==2){
	    xsparse[2]=-4.5;
	  }
	}
      }
      
      if(theBin>=0 && theBin<kNbins){
	if(fTrkSign>0){
	  if(minPosMean==0) fHistDCAPosPi[theBin]->Fill(fImpactXY);
	  else if(minPosMean==1) fHistDCAPosK[theBin]->Fill(fImpactXY);
	  else if(minPosMean==2) fHistDCAPosP[theBin]->Fill(fImpactXY); 
	}else{
	  if(minPosMean==0)fHistDCANegPi[theBin]->Fill(fImpactXY);
	  else if(minPosMean==1)fHistDCANegK[theBin]->Fill(fImpactXY);
	  else if(minPosMean==2)fHistDCANegP[theBin]->Fill(fImpactXY);
	}
      }
    } 
    
    /////////////////////
    //DCA cut on xy and z
    /////////////////////
    if(!DCAcutXY(fImpactXY,pt)){
      //cout<<"DCA cut "<<endl;
      if(xsparse[0]>-999 && fUseThnSparse) fHistSparse->Fill(xsparse,1);
      continue;
    }
    
    fHistNTracks->Fill(kPassDCAxycut);
    if(fTrkSign>0)fHistNTracksPos->Fill(kPassDCAxycut);
    if(fTrkSign<0)fHistNTracksNeg->Fill(kPassDCAxycut);
    
    //Information that the DCA cut has been passed
    if(xsparse[2]>0) xsparse[2] = xsparse[2]+1;
    else if(xsparse[2]<0 && xsparse[2]>-999) xsparse[2] = xsparse[2]-1;
    
    Int_t jpart=-1;
    
    //Filling Histos for Reco Efficiency
    //information from the MC kinematics
    if(fMC){
      Int_t trklabel = TMath::Abs(track->GetLabel());
      if(trklabel<0)isph=-1;
      if(trklabel>=0){
	part = (TParticle*)stack->Particle(trklabel);
	pdgPart = part->GetPDG();
	code = pdgPart->PdgCode();
	for(Int_t j=0; j<3; j++){//jpart is 0 for pion, 1 for kaon, 2 for proton
	  if(TMath::Abs(code)==listcode[j]){
	    jpart=j;
	    break;
	  }
	}
	if(jpart>=0){
	  if(pdgPart->Charge()>0) signMC=1;
	  else signMC=-1;
	  ptMC=part->Pt();
	  xsparse[6]=ptMC;
	  if(stack->IsPhysicalPrimary(trklabel)){
	    if(signMC>0){
	      fHistPrimMCposReco[jpart]->Fill(ptMC);
	      if(TMath::Abs(y[jpart])<fMaxY) fHistPrimMCposRecoEtaY[jpart]->Fill(ptMC);
	    }else{
	      fHistPrimMCnegReco[jpart]->Fill(ptMC);
	      if(TMath::Abs(y[jpart])<fMaxY) fHistPrimMCnegRecoEtaY[jpart]->Fill(ptMC);
	    }
	  }
	  else{ 
	    Int_t indexMoth=part->GetFirstMother();
	    if(indexMoth>=0){
	      TParticle* moth = stack->Particle(indexMoth);
	      Float_t codemoth = TMath::Abs(moth->GetPdgCode());
	      mfl = Int_t (codemoth/ TMath::Power(10, Int_t(TMath::Log10(codemoth))));
	    }
	    uniqueID = part->GetUniqueID();
	    if(mfl==3 && uniqueID == kPDecay){ // strangeness
	      if(signMC>0) fHistSecStrMCposReco[jpart]->Fill(ptMC);
	      else fHistSecStrMCnegReco[jpart]->Fill(ptMC);
	      
	    }else{//material
	      if(signMC>0) fHistSecMatMCposReco[jpart]->Fill(ptMC);
	      else fHistSecMatMCnegReco[jpart]->Fill(ptMC);      
	    }
	  }
	  
	  Double_t ptcontainer[4] = {ptMC,pt,jpart + 0.5,0.5*signMC};
	  fHistPtCorr->Fill(ptcontainer);
	}
      }
    }
    
    //Nsigma histos with eventually the MC truth
    
    //asymmetric approach
    if(TMath::Abs(yPartMean)<fMaxY && minPosMean>-1){
      //nsigma histos
      if(fTrkSign>0) fHistPosNSigmaMean[minPosMean]->Fill(pt);
      else fHistNegNSigmaMean[minPosMean]->Fill(pt);
      
      if(fMC){
	//nsigma histos with MC truth on PID
	if(TMath::Abs(code)==listcode[minPosMean]){
	  if(fTrkSign>0) fHistPosNSigmaMCMean[minPosMean]->Fill(pt);
	  else fHistNegNSigmaMCMean[minPosMean]->Fill(pt);
	}
	//nsigma histos with MC truth on IsPhysicalPrimary
	if(isph==1){
	  if(fTrkSign>0){
	    fHistPosNSigmaPrimMean[minPosMean]->Fill(pt);
	    
	  }
	  else fHistNegNSigmaPrimMean[minPosMean]->Fill(pt);
	  //nsigma histos with MC truth on IsPhysicalPrimary and PID
	  if(TMath::Abs(code)==listcode[minPosMean]){
	    if(fTrkSign>0) fHistPosNSigmaPrimMCMean[minPosMean]->Fill(pt);
	    else fHistNegNSigmaPrimMCMean[minPosMean]->Fill(pt);
	  }
	}
      }
    }
    
    //symmetric bands
    if(min<fMinNSigma && TMath::Abs(yPart)<fMaxY && minPos>-1){
      //nsigma histos
      if(fTrkSign>0){
	xsparse[3]=minPos+0.5;
	fHistPosNSigma[minPos]->Fill(pt);
      }
      else{
	xsparse[3]=-minPos-0.5;
	fHistNegNSigma[minPos]->Fill(pt);
      }
      
      if(fMC){
	//nsigma histos with MC truth on PID
	if(TMath::Abs(code)==listcode[minPos]){
	  if(fTrkSign>0) fHistPosNSigmaMC[minPos]->Fill(pt);
	  else fHistNegNSigmaMC[minPos]->Fill(pt);
	}
	//nsigma histos with MC truth on IsPhysicalPrimary
	if(isph==1){
	  if(fTrkSign>0) fHistPosNSigmaPrim[minPos]->Fill(pt);
	  else fHistNegNSigmaPrim[minPos]->Fill(pt);
	  //nsigma histos with MC truth on IsPhysicalPrimary and PID
	  if(TMath::Abs(code)==listcode[minPos]){
	    if(fTrkSign>0) fHistPosNSigmaPrimMC[minPos]->Fill(pt);
	    else fHistNegNSigmaPrimMC[minPos]->Fill(pt);
	  }
	}
      }
    }
    
    //integral approach histograms
    if(theBin>=0 && theBin<kNbins){
      if(fTrkSign>0){
	if(TMath::Abs(y[0]) < fMaxY)fHistPosPi[theBin]->Fill(logdiff[0]);
	if(TMath::Abs(y[1]) < fMaxY)fHistPosK[theBin]->Fill(logdiff[1]);
	if(TMath::Abs(y[2]) < fMaxY)fHistPosP[theBin]->Fill(logdiff[2]);
	if(fMC){
	  if(TMath::Abs(y[0])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCPosOtherHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==11)fHistMCPosElHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==211)fHistMCPosPiHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==321)fHistMCPosKHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==2212)fHistMCPosPHypPion[theBin]->Fill(logdiff[0]);
	  }
	  if(TMath::Abs(y[1])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCPosOtherHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==11)fHistMCPosElHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==211)fHistMCPosPiHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==321)fHistMCPosKHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==2212)fHistMCPosPHypKaon[theBin]->Fill(logdiff[1]);
	  }
	  if(TMath::Abs(y[2])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCPosOtherHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==11)fHistMCPosElHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==211)fHistMCPosPiHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==321)fHistMCPosKHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==2212)fHistMCPosPHypProton[theBin]->Fill(logdiff[2]);
	  }
	}
      }else{
	if(TMath::Abs(y[0]) < fMaxY)fHistNegPi[theBin]->Fill(logdiff[0]);
	if(TMath::Abs(y[1]) < fMaxY)fHistNegK[theBin]->Fill(logdiff[1]);
	if(TMath::Abs(y[2]) < fMaxY)fHistNegP[theBin]->Fill(logdiff[2]);
	if(fMC){
	  if(TMath::Abs(y[0])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCNegOtherHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==11)fHistMCNegElHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==211)fHistMCNegPiHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==321)fHistMCNegKHypPion[theBin]->Fill(logdiff[0]);
	    if(TMath::Abs(code)==2212)fHistMCNegPHypPion[theBin]->Fill(logdiff[0]);
	  }
	  if(TMath::Abs(y[1])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCNegOtherHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==11)fHistMCNegElHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==211)fHistMCNegPiHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==321)fHistMCNegKHypKaon[theBin]->Fill(logdiff[1]);
	    if(TMath::Abs(code)==2212)fHistMCNegPHypKaon[theBin]->Fill(logdiff[1]);
	  }
	  if(TMath::Abs(y[2])<fMaxY){
	    if(TMath::Abs(code)!=11 && jpart<0)fHistMCNegOtherHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==11)fHistMCNegElHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==211)fHistMCNegPiHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==321)fHistMCNegKHypProton[theBin]->Fill(logdiff[2]);
	    if(TMath::Abs(code)==2212)fHistMCNegPHypProton[theBin]->Fill(logdiff[2]);
	  }
	}
      }
    }
    
    
    //Filling Time!
    if(xsparse[0]>-999 && fUseThnSparse) fHistSparse->Fill(xsparse,1);
    
    //fill propaganda plot with dedx
    fHistDEDX->Fill(track->GetP(),fdedx);
    fHistDEDXdouble->Fill(track->GetP()*fTrkSign,fdedx);
    if(fLowMult==-1 && fUpMult==-1){
      Double_t xdedx[4];
      xdedx[0] = track->GetP();
      xdedx[1] = fdedx;
      xdedx[2] = fMult;
      xdedx[3] = minPosMean+0.5;
      fHistDEDXMulti->Fill(xdedx);
      xdedx[0] = track->GetP()*fTrkSign;
      fHistDEDXdoubleMulti->Fill(xdedx);
    }
    //fill charge distribution histo to check the calibration
    for(Int_t j=0;j<4;j++){
      if(dedxLay[j]<5) continue;
      fHistCharge[j]->Fill(dedxLay[j]);
    }
  }//Ends loop on tracks
  
  
  // Post output data.
  PostData(1,fOutput);
  PostData(2,fListCuts);
  
  infomsg("End of UserExec");
}      

//________________________________________________________________________
void AliAnalysisTaskSEITSsaSpectraMultiplicity::Terminate(Option_t *) {
  // Merge output
  // Called once at the end of the query
  
  fOutput = dynamic_cast<TList*>(GetOutputData(1));
  if (!fOutput) {
    printf("ERROR: fOutput not available\n");
    return;
  }
  if(fHistNEvents){
    fHistNEvents = dynamic_cast<TH1F*>(fOutput->FindObject("fHistNEvents"));
    printf("Number of Analyzed Events = %f\n",fHistNEvents->GetBinContent(1));
  }else{
    printf("ERROR: fHistNEvents not available\n");
  }
  
  infomsg("End of Terminate");
  return;
}
