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
//
// QA class of Heavy Flavor quark and fragmeted/decayed particles
// -Check kinematics of Heavy Quarks/hadrons, and decayed leptons
//    pT, rapidity
//    decay lepton kinematics w/wo acceptance
//    heavy hadron decay length, electron pT fraction carried from decay
// -Check yield of Heavy Quarks/hadrons
//    Number of produced heavy quark
//    Number of produced hadron of given pdg code
//
//
// Authors:
//   MinJung Kweon <minjung@physi.uni-heidelberg.de>
//

#include <TH2F.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TParticle.h>
#include "TTreeStream.h"

#include <AliLog.h>
#include <AliMCEvent.h>
#include <AliGenEventHeader.h>
#include <AliAODMCParticle.h>
#include <AliStack.h>

#include "AliHFEmcQA.h"
#include "AliHFEtools.h"
#include "AliHFEcollection.h"

ClassImp(AliHFEmcQA)

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA() :
fMCEvent(NULL)
  ,fMCHeader(NULL)
  ,fMCArray(NULL)
  ,fQAhistos(NULL)
  ,fMCQACollection(NULL)
  ,fNparents(0)
  ,fCentrality(0)
  ,fPerCentrality(-1)
  ,fIsPbPb(kFALSE)
  ,fIsppMultiBin(kFALSE)
  ,fContainerStep(0)
  ,fIsDebugStreamerON(kFALSE)
  ,fRecPt(-999)
  ,fRecEta(-999)
  ,fRecPhi(-999)
  ,fLyrhit(0)
  ,fLyrstat(0)
  ,fHfeImpactR(-999)
  ,fHfeImpactnsigmaR(-999)
  ,fTreeStream(NULL)
  ,fGetWeightHist(kFALSE)
{
  // Default constructor
  for(Int_t mom = 0; mom < 9; mom++){
    fhD[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 50; mom++){
    fHeavyQuark[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 2; mom++){
    fIsHeavy[mom] = 0;
  }
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
}

//_______________________________________________________________________________________________
AliHFEmcQA::AliHFEmcQA(const AliHFEmcQA&p):
  TObject(p)
  ,fMCEvent(NULL)
  ,fMCHeader(NULL)
  ,fMCArray(NULL)
  ,fQAhistos(p.fQAhistos)
  ,fMCQACollection(p.fMCQACollection)
  ,fNparents(p.fNparents)
  ,fCentrality(0)
  ,fPerCentrality(-1)
  ,fIsPbPb(kFALSE)
  ,fIsppMultiBin(kFALSE)
  ,fContainerStep(0)
  ,fIsDebugStreamerON(kFALSE)
  ,fRecPt(-999)
  ,fRecEta(-999)
  ,fRecPhi(-999)
  ,fLyrhit(0)
  ,fLyrstat(0)
  ,fHfeImpactR(0)
  ,fHfeImpactnsigmaR(0)
  ,fTreeStream(NULL)
  ,fGetWeightHist(kFALSE)
{
  // Copy constructor
  for(Int_t mom = 0; mom < 9; mom++){
    fhD[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 50; mom++){
    fHeavyQuark[mom] = NULL;
  }
  for(Int_t mom = 0; mom < 2; mom++){
    fIsHeavy[mom] = 0;
  }
  memset(fElecBackgroundFactor, 0, sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memset(fBinLimit, 0, sizeof(Double_t) * (kBgPtBins+1));
}
//_______________________________________________________________________________________________
AliHFEmcQA&
AliHFEmcQA::operator=(const AliHFEmcQA &)
{
  // Assignment operator
  
  AliInfo("Not yet implemented.");
  return *this;
}

//_______________________________________________________________________________________________
AliHFEmcQA::~AliHFEmcQA()
{
  // Destructor
  
  if(fTreeStream && fIsDebugStreamerON) delete fTreeStream;
  AliInfo("Analysis Done.");
}
//_______________________________________________________________________________________________
void AliHFEmcQA::PostAnalyze() const
{
  //
  // Post analysis
  //
}
//_______________________________________________________________________________________________
void AliHFEmcQA::SetBackgroundWeightFactor(Double_t *elecBackgroundFactor, Double_t *binLimit)
{
  //
  // copy background weighting factors into data member
  //
  
  memcpy(fElecBackgroundFactor,elecBackgroundFactor,sizeof(Double_t) * kElecBgSpecies * kBgPtBins * kCentBins * kBgLevels);
  memcpy(fBinLimit,binLimit,sizeof(Double_t) * (kBgPtBins+1));
}
//__________________________________________
void AliHFEmcQA::CreatDefaultHistograms(TList * const qaList)
{      
  //
  // make default histograms
  //
  
  if(!qaList) return;
  
  fQAhistos = qaList;
  fQAhistos->SetName("MCqa");
  
  CreateHistograms(AliHFEmcQA::kCharm);               // create histograms for charm
  CreateHistograms(AliHFEmcQA::kBeauty);              // create histograms for beauty
  CreateHistograms(AliHFEmcQA::kOthers);              // create histograms for beauty
  
  // prepare 2D(pt vs Y) histogram for D spectra, we consider following 9 particles
  const Int_t nbspecies = 9;
  TString kDspecies[nbspecies];
  kDspecies[0]="411";   //D+
  kDspecies[1]="421";   //D0
  kDspecies[2]="431";   //Ds+
  kDspecies[3]="4122";  //Lambdac+
  kDspecies[4]="4132";  //Ksic0
  kDspecies[5]="4232";  //Ksic+
  kDspecies[6]="4332";  //OmegaC0
  kDspecies[7]="413";   //D*(2010)+
  kDspecies[8]="423";   //D*(2007)0

  const Double_t kPtbound[2] = {0.1, 20.}; //bin taken for considering inclusive e analysis binning
  Int_t iBin[2];
  iBin[0] = 44; // bins in pt for log binning
  iBin[1] = 23; // bins in pt for pi0 measurement binning
  //Double_t* binEdges[1];
  //binEdges[0] =  AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);

  // bin size is chosen to consider ALICE D measurement
  const Int_t nptbins = 15;
  const Int_t nybins = 9;
  Double_t xbins[nptbins+1]={0,0.5,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,12,16,24,32,40,50}; //pt binning for the final 7 TeV D measurement 
  Double_t ybins[nybins+1]={-7.5,-1.0,-0.9,-0.8,-0.5,0.5,0.8,0.9,1.0,7.5}; // y binning
  TString hname;
  for (Int_t iDmeson=0; iDmeson<nbspecies; iDmeson++){
     hname = "Dmeson"+kDspecies[iDmeson];
     fhD[iDmeson] = new TH2F(hname,hname+";p_{T} (GeV/c)",nptbins,xbins,nybins,ybins);
     if(fQAhistos) fQAhistos->Add(fhD[iDmeson]);
  }

  const Double_t kPtRange[24] = {0.,0.3,0.4,0.5,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.,3.5,4.,5.,6.,7.,20.,30.}; // to cope with Ana's bin

  Int_t kNcent;
  if(fIsPbPb) kNcent=11;
  else
  {
      if(fIsppMultiBin) kNcent=8;
      else kNcent = 1;
  }

  fMCQACollection = new AliHFEcollection("TaskMCQA", "MC QA histos for meason pt spectra");

  for(Int_t centbin=0; centbin<kNcent; centbin++)
  {
      fMCQACollection->CreateTH1Farray(Form("pionspectra_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("etaspectra_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("omegaspectra_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("phispectra_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("etapspectra_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[1],kPtRange);
      fMCQACollection->CreateTH1Farray(Form("rhospectra_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[1],kPtRange);

      fMCQACollection->CreateTH1F(Form("pionspectraLog_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etaspectraLog_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("omegaspectraLog_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("phispectraLog_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etapspectraLog_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("rhospectraLog_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("kaonspectraLog_centrbin%i",centbin), "kaon yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("k0LspectraLog_centrbin%i",centbin), "k0L yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("k0SspectraLog_centrbin%i",centbin), "k0S yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("lamdaspectraLog_centrbin%i",centbin), "lamda yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("sigmaspectraLog_centrbin%i",centbin), "sigma yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);

      fMCQACollection->CreateTH2F(Form("pionspectraLog2D_centrbin%i",centbin), "pion yields: MC p_{t} ", 42, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH2F(Form("etaspectraLog2D_centrbin%i",centbin), "eta yields: MC p_{t} ", 42, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH2F(Form("omegaspectraLog2D_centrbin%i",centbin), "omega yields: MC p_{t} ", 40, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH2F(Form("phispectraLog2D_centrbin%i",centbin), "phi yields: MC p_{t} ", 42, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH2F(Form("etapspectraLog2D_centrbin%i",centbin), "etap yields: MC p_{t} ", 45, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH2F(Form("rhospectraLog2D_centrbin%i",centbin), "rho yields: MC p_{t} ", 42, -1.5, 40.5, iBin[0],kPtbound[0], kPtbound[1], 1);

      fMCQACollection->CreateTH1F(Form("piondaughters_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etadaughters_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("omegadaughters_centrbin%i",centbin), "omega yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("phidaughters_centrbin%i",centbin), "phi yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etapdaughters_centrbin%i",centbin), "etap yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("rhodaughters_centrbin%i",centbin), "rho yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);

      fMCQACollection->CreateTH1F(Form("pionspectraPrimary_centrbin%i",centbin), "pion yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
      fMCQACollection->CreateTH1F(Form("etaspectraPrimary_centrbin%i",centbin), "eta yields: MC p_{t} ", iBin[0],kPtbound[0], kPtbound[1], 1);
  }

  fQAhistos->Add(fMCQACollection->GetList());

  if(!fTreeStream && fIsDebugStreamerON){
   fTreeStream = new TTreeSRedirector(Form("HFEmcqadebugTree%s.root", GetName()));
  }

}
  
//__________________________________________
void AliHFEmcQA::CreateHistograms(const Int_t kquark) 
{
  // create histograms

  if (!(kquark == kCharm || kquark == kBeauty || kquark == kOthers)) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 

  TString kqTypeLabel[fgkqType];
  if (kquark == kCharm){
    kqTypeLabel[kQuark]="c";
    kqTypeLabel[kantiQuark]="cbar";
    kqTypeLabel[kHadron]="cHadron";
    kqTypeLabel[keHadron]="ceHadron";
    kqTypeLabel[kDeHadron]="nullHadron";
    kqTypeLabel[kElectron]="ce";
    kqTypeLabel[kElectron2nd]="nulle";
  } else if (kquark == kBeauty){
    kqTypeLabel[kQuark]="b";
    kqTypeLabel[kantiQuark]="bbar";
    kqTypeLabel[kHadron]="bHadron";
    kqTypeLabel[keHadron]="beHadron";
    kqTypeLabel[kDeHadron]="bDeHadron";
    kqTypeLabel[kElectron]="be";
    kqTypeLabel[kElectron2nd]="bce";
  } else if (kquark == kOthers){
    kqTypeLabel[kGamma-4]="gammae";
    kqTypeLabel[kPi0-4]="pi0e";
    kqTypeLabel[kElse-4]="elsee";
    kqTypeLabel[kMisID-4]="miside";
  }

  TString kqEtaRangeLabel[fgkEtaRanges];
  kqEtaRangeLabel[0] = "mcqa_";
  kqEtaRangeLabel[1] = "mcqa_barrel_";
  kqEtaRangeLabel[2] = "mcqa_unitY_";

  //const Double_t kPtbound[2] = {0.1, 20.}; //bin taken for considering inclusive e analysis binning
  const Int_t nptbinning1 = 35;
  Int_t iBin[2];
  iBin[0] = 44; // bins in pt
  iBin[1] = nptbinning1; // bins in pt
  //Double_t* binEdges[1];
  //binEdges[0] =  AliHFEtools::MakeLogarithmicBinning(iBin[0], kPtbound[0], kPtbound[1]);

  // new binning for final electron analysis
  const Double_t kPtbinning1[nptbinning1+1] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1., 1.1, 1.2, 1.3, 1.4, 1.5, 1.75, 2., 2.25, 2.5, 2.75, 3., 3.5, 4., 4.5, 5., 5.5, 6., 7., 8., 10., 12., 14., 16., 18., 20.};

  const Int_t ndptbins = 500;
  Double_t xcorrbin[ndptbins+1];
  for (int icorrbin = 0; icorrbin< ndptbins+1; icorrbin++){
    xcorrbin[icorrbin]=icorrbin*0.1;
  }

  Int_t fCentrmax = 0;
  if(!fIsPbPb&&!fIsppMultiBin) fCentrmax=0;
  else fCentrmax = kCentBins;
  TString hname; 
  if(kquark == kOthers){
   for (Int_t icut = 0; icut < fgkEtaRanges; icut++ ){
       for (Int_t iqType = 0; iqType < 4; iqType++ ){
	   for(Int_t icentr = 0; icentr < (fCentrmax+1); icentr++)
	   {
	       hname = kqEtaRangeLabel[icut]+"Pt_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fPt = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;p_{T} (GeV/c)",hname.Data(),icentr),60,0.25,30.25);
	       hname = kqEtaRangeLabel[icut]+"Y_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fY = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;y",hname.Data(),icentr),150,-7.5,7.5);
	       hname = kqEtaRangeLabel[icut]+"Eta_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fEta = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;eta",hname.Data(),icentr),150,-7.5,7.5);
	       // Fill List
	       if(fQAhistos) fHist[iq][iqType][icut][icentr].FillList(fQAhistos);
	   }
       }
   }
   return;
  }
  for (Int_t icut = 0; icut < fgkEtaRanges; icut++ ){
   for (Int_t iqType = 0; iqType < fgkqType; iqType++ ){
       if (iqType < keHadron && icut > 0) continue; // don't duplicate histogram for quark and hadron
        for(Int_t icentr = 0; icentr<(fCentrmax+1); icentr++)
	   {
	       hname = kqEtaRangeLabel[icut]+"PdgCode_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fPdgCode = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;PdgCode",hname.Data(),icentr),20001,-10000.5,10000.5);
	       hname = kqEtaRangeLabel[icut]+"Pt_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fPt = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;p_{T} (GeV/c)",hname.Data(),icentr),iBin[1],kPtbinning1); // new binning
	       hname = kqEtaRangeLabel[icut]+"Y_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fY = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;y",hname.Data(),icentr),150,-7.5,7.5);
	       hname = kqEtaRangeLabel[icut]+"Eta_"+kqTypeLabel[iqType];
	       fHist[iq][iqType][icut][icentr].fEta = new TH1F(Form("%sCentr_%i",hname.Data(),icentr),Form("%sCentr_%i;eta",hname.Data(),icentr),150,-7.5,7.5);
	       // Fill List
	       if(fQAhistos) fHist[iq][iqType][icut][icentr].FillList(fQAhistos);
	   }
   }
  }

  for (Int_t icut = 0; icut < fgkEtaRanges; icut++ ){
    hname = kqEtaRangeLabel[icut]+"PtCorr_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fPtCorr = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); 
    hname = kqEtaRangeLabel[icut]+"PtCorrDp_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fPtCorrDp= new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1);
    hname = kqEtaRangeLabel[icut]+"PtCorrD0_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fPtCorrD0 = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1);
    hname = kqEtaRangeLabel[icut]+"PtCorrDrest_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fPtCorrDrest = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1);

    hname = kqEtaRangeLabel[icut]+"ePtRatio_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",200,0,20,100,0,1);
    hname = kqEtaRangeLabel[icut]+"DePtRatio_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fDePtRatio = new TH2F(hname,hname+";p_{T} (GeV/c);momentum fraction",100,0,20,100,0,1);
    hname = kqEtaRangeLabel[icut]+"eDistance_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].feDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,20,200,0,2);
    hname = kqEtaRangeLabel[icut]+"DeDistance_"+kqTypeLabel[kQuark];
    fHistComm[iq][icut].fDeDistance= new TH2F(hname,hname+";p_{T} (GeV/c);distance (cm)",100,0,20,200,0,2);

    if(icut <1){
      hname = kqEtaRangeLabel[icut]+"PtCorrDinein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDinein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDineout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDineout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDoutein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDoutein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDouteout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDouteout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning

      hname = kqEtaRangeLabel[icut]+"PtCorrDpDinein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDpDinein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDpDineout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDpDineout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDpDoutein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDpDoutein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDpDouteout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDpDouteout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning

      hname = kqEtaRangeLabel[icut]+"PtCorrD0Dinein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrD0Dinein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrD0Dineout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrD0Dineout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrD0Doutein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrD0Doutein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrD0Douteout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrD0Douteout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning

      hname = kqEtaRangeLabel[icut]+"PtCorrDrestDinein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDrestDinein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDrestDineout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDrestDineout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDrestDoutein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDrestDoutein = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrDrestDouteout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrDrestDouteout = new TH2F(hname,hname+";D p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning

      hname = kqEtaRangeLabel[icut]+"fEtaCorrD_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrD = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrDp_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrDp = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrD0_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrD0 = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrDrest_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrDrest = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 

      hname = kqEtaRangeLabel[icut]+"fEtaCorrGD_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrGD = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrGDp_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrGDp = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrGD0_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrGD0 = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrGDrest_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrGDrest = new TH2F(hname,hname+";D Y;e eta",200,-10,10,200,-10,10); 

      hname = kqEtaRangeLabel[icut]+"fEtaCorrB_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrB = new TH2F(hname,hname+";B Y;e eta",200,-10,10,200,-10,10); 
      hname = kqEtaRangeLabel[icut]+"fEtaCorrGB_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fEtaCorrGB = new TH2F(hname,hname+";B Y;e eta",200,-10,10,200,-10,10); 

      hname = kqEtaRangeLabel[icut]+"PtCorrBinein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrBinein = new TH2F(hname,hname+";B p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrBineout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrBineout = new TH2F(hname,hname+";B p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrBoutein_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrBoutein = new TH2F(hname,hname+";B p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
      hname = kqEtaRangeLabel[icut]+"PtCorrBouteout_"+kqTypeLabel[kQuark];
      fHistComm[iq][icut].fPtCorrBouteout = new TH2F(hname,hname+";B p_{T} (GeV/c);e p_{T} (GeV/c)",ndptbins,xcorrbin,iBin[1],kPtbinning1); // new binning
    }
    if(fQAhistos) fHistComm[iq][icut].FillList(fQAhistos);
  }


  hname = kqEtaRangeLabel[0]+"Nq_"+kqTypeLabel[kQuark];
  fHistComm[iq][0].fNq = new TH1F(hname,hname,50,-0.5,49.5);
  hname = kqEtaRangeLabel[0]+"ProcessID_"+kqTypeLabel[kQuark];
  fHistComm[iq][0].fProcessID = new TH1F(hname,hname,21,-10.5,10.5);
  
}

//__________________________________________
void AliHFEmcQA::Init()
{
  // called at begining every event
  
  for (Int_t i=0; i<2; i++){
     fIsHeavy[i] = 0;
  } 

  fNparents = 7;

  fParentSelect[0][0] =  411; //D+  
  fParentSelect[0][1] =  421; //D0
  fParentSelect[0][2] =  431; //Ds+
  fParentSelect[0][3] = 4122; //Lambdac+
  fParentSelect[0][4] = 4132; //Ksic0
  fParentSelect[0][5] = 4232; //Ksic+
  fParentSelect[0][6] = 4332; //OmegaC0

  fParentSelect[1][0] =  511; //B0
  fParentSelect[1][1] =  521; //B+
  fParentSelect[1][2] =  531; //Bs0
  fParentSelect[1][3] = 5122; //Lambdab0
  fParentSelect[1][4] = 5132; //Ksib-
  fParentSelect[1][5] = 5232; //Ksib0
  fParentSelect[1][6] = 5332; //Omegab-


}

//__________________________________________
void AliHFEmcQA::GetMesonKine() 
{
  //
  // get meson pt spectra
  //

  AliVParticle *mctrack2 = NULL;
  AliMCParticle *mctrack0 = NULL;
  AliVParticle *mctrackdaugt= NULL;
  AliMCParticle *mctrackd= NULL;
  Int_t id1=0, id2=0;

 
  if(fCentrality>11) {
      AliWarning(Form("Centrality out of histogram array limits: %d", fCentrality));
      return;
  }
  if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

  for(Int_t imc = 0; imc <fMCEvent->GetNumberOfPrimaries(); imc++){
     if(!(mctrack2 = fMCEvent->GetTrack(imc))) continue;
     TParticle* mcpart0 = fMCEvent->Stack()->Particle(imc);
     if(!mcpart0) continue;
     mctrack0 = dynamic_cast<AliMCParticle *>(mctrack2);
     if(!mctrack0) continue;

//     if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;
     Float_t mcsource = 1;
     if(fGetWeightHist) mcsource = GetElecSource(mctrack0, kFALSE);

     if(TMath::Abs(mctrack0->PdgCode()) == 111) // pi0 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            if(mcpart0->IsPrimary()){
                fMCQACollection->Fill(Form("pionspectraPrimary_centrbin%i",fCentrality),mctrack0->Pt());
            }
            fMCQACollection->Fill(Form("pionspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("pionspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("pionspectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
	}
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==2)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("piondaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 221) // eta 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            if(mcpart0->IsPrimary()){
                fMCQACollection->Fill(Form("etaspectraPrimary_centrbin%i",fCentrality),mctrack0->Pt());
            }
            fMCQACollection->Fill(Form("etaspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etaspectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
          } 
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==2||(id2-id1)==3)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("etadaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 223) // omega
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("omegaspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("omegaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("omegaspectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
          }
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==1||(id2-id1)==2)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("omegadaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 333) // phi 
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("phispectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("phispectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("phispectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
          } 
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==1)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("phidaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 331) // eta prime
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("etapspectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etapspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("etapspectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
          }
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==2||(id2-id1)==3)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("etapdaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 113) // rho
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("rhospectra_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("rhospectraLog_centrbin%i",fCentrality),mctrack0->Pt());
            fMCQACollection->Fill(Form("rhospectraLog2D_centrbin%i",fCentrality),mcsource,mctrack0->Pt());
          }
          id1=mctrack0->GetDaughterFirst();
          id2=mctrack0->GetDaughterLast();
          if(!((id2-id1)==1)) continue;
          for(int idx=id1; idx<=id2; idx++){
            if(!(mctrackdaugt = fMCEvent->GetTrack(idx))) continue;
            if(!(mctrackd = dynamic_cast<AliMCParticle *>(mctrackdaugt))) continue;
            if(TMath::Abs(mctrackd->PdgCode()) == 11 && TMath::Abs(mctrackd->Eta())<0.8)
             fMCQACollection->Fill(Form("rhodaughters_centrbin%i",fCentrality),mctrackd->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 321) // kaon+-
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("kaonspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 130) // k0L
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("k0LspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     
     else if(TMath::Abs(mctrack0->PdgCode()) == 310) // k0S
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("k0SspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 3122) // lamda
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("lamdaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
     else if(TMath::Abs(mctrack0->PdgCode()) == 3222) // sigma
       {
          if(TMath::Abs(AliHFEtools::GetRapidity(mcpart0))<0.8) {
            fMCQACollection->Fill(Form("sigmaspectraLog_centrbin%i",fCentrality),mctrack0->Pt());
          }
       }
    }

}
//__________________________________________
void AliHFEmcQA::GetQuarkKine(TParticle *part, Int_t iTrack, const Int_t kquark) 
{
  // get heavy quark kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 

    if (iTrack < 0 || !part) { 
      AliDebug(1, "Stack label is negative or no mcparticle, return\n");
      return; 
    }

    if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

    AliMCParticle *mctrack = NULL;
    Int_t partPdgcode = TMath::Abs(part->GetPdgCode());

    // select heavy hadron or not fragmented heavy quark 
    if ( int(partPdgcode/100.)==kquark || int(partPdgcode/1000.)==kquark || (partPdgcode==kquark && (part->GetNDaughters()==0 && iTrack>5)) ){ 

      TParticle *partMother;
      Int_t iLabel;

      if (partPdgcode == kquark){ // in case of not fragmented heavy quark  
        partMother = part; 
        iLabel = iTrack;
      } else{ // in case of heavy hadron, start to search for mother heavy parton 
        iLabel = part->GetFirstMother(); 
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
        if (iLabel>-1) { partMother = mctrack->Particle(); }
        else {
          AliDebug(1, "Stack label is negative, return\n");
          return; 
        }
      }

      // heavy parton selection as a mother of heavy hadron 
      // if the heavy particle comes from string which is denoted as particle status 12|12|12...12|11,[PYTHIA p.60]
      // in this case, the mother of heavy particle can be one of the fragmented parton of the string
      // should I make a condition that partMother should be quark or diquark? -> not necessary
      if ( TMath::Abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11) ){
      //if ( TMath::Abs(partMother->GetPdgCode()) == kquark || (partMother->GetStatusCode() == 11 || partMother->GetStatusCode() == 12) ){

        if ( TMath::Abs(partMother->GetPdgCode()) != kquark ){
          // search fragmented partons in the same string
          Bool_t isSameString = kTRUE; 
          for (Int_t i=1; i<fgkMaxIter; i++){
             iLabel = iLabel - 1;
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
             if (iLabel>-1) { partMother = mctrack->Particle(); }
             else {
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }
             if ( TMath::Abs(partMother->GetPdgCode()) == kquark ) break;
             if ( partMother->GetStatusCode() != 12 ) isSameString = kFALSE;
             if (!isSameString) return; 
          }
        }
        AliDebug(1, "Can not find heavy parton of this heavy hadron in the string, return\n");
        if (TMath::Abs(partMother->GetPdgCode()) != kquark) return; 

        if (fIsHeavy[iq] >= 50) return;  
        fHeavyQuark[fIsHeavy[iq]] = partMother;
        fIsHeavy[iq]++;

        // fill kinematics for heavy parton
        if (partMother->GetPdgCode() > 0) { // quark
          fHist[iq][kQuark][0][fCentrality].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kQuark][0][fCentrality].fPt->Fill(partMother->Pt());
          fHist[iq][kQuark][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMother));
          fHist[iq][kQuark][0][fCentrality].fEta->Fill(partMother->Eta());
        } else{ // antiquark
          fHist[iq][kantiQuark][0][fCentrality].fPdgCode->Fill(partMother->GetPdgCode());
          fHist[iq][kantiQuark][0][fCentrality].fPt->Fill(partMother->Pt());
          fHist[iq][kantiQuark][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMother));
          fHist[iq][kantiQuark][0][fCentrality].fEta->Fill(partMother->Eta());
        }

      } // end of heavy parton slection loop 

    } // end of heavy hadron or quark selection

}

//__________________________________________
void AliHFEmcQA::EndOfEventAna(const Int_t kquark)
{
  // end of event analysis

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return; 
  }
  Int_t iq = kquark - kCharm; 


  // # of heavy quark per event
  AliDebug(1,Form("Number of heavy quark in this event = %d \n",fIsHeavy[iq]));
  fHistComm[iq][0].fNq->Fill(fIsHeavy[iq]);

  Int_t motherID[fgkMaxGener];
  Int_t motherType[fgkMaxGener];
  Int_t motherLabel[fgkMaxGener];
  Int_t ancestorPdg[fgkMaxGener];
  Int_t ancestorLabel[fgkMaxGener];

  for (Int_t i = 0; i < fgkMaxGener; i++){ // initialization
     motherID[i] = 0;
     motherType[i] = 0;
     motherLabel[i] = 0;
     ancestorPdg[i] = 0;
     ancestorLabel[i] = 0;
  }


  // check history of found heavy quarks
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){

     if(!fHeavyQuark[i]) return;

     ancestorLabel[0] = i;
     ancestorPdg[0] = fHeavyQuark[i]->GetPdgCode(); 
     ancestorLabel[1] = fHeavyQuark[i]->GetFirstMother(); 

     AliDebug(1,Form("pdg code= %d\n",ancestorPdg[0]));
     AliDebug(1,Form("ancestor label= %d\n",ancestorLabel[1]));

     Int_t ig = 1;
     while (ancestorLabel[ig] != -1){
          // in case there is mother, get mother's pdg code and grandmother's label
          IdentifyMother(ancestorLabel[ig], ancestorPdg[ig], ancestorLabel[ig+1]); 
          // if mother is still heavy, find again mother's ancestor
          if (ancestorPdg[ig-1] == ancestorPdg[ig]) {
            ig++;
            continue; // if it is from same heavy
          }
          // if the heavy's mother is not heavy, check the mother's label to know if it comes from inital or final parton shower
          if (IsFromInitialShower(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          if (IsFromFinalParton(ancestorLabel[ig],motherID[i],motherType[i],motherLabel[i])) break;
          // if it is not the above case, something is strange
          ReportStrangeness(motherID[i],motherType[i],motherLabel[i]);
          break;
     } 
     if (ancestorLabel[ig] == -1){ // from hard scattering
       HardScattering(kquark, motherID[i],motherType[i], motherLabel[i]);
     }

  } // end of found heavy quark loop


  // check process type
  Int_t processID = 0;
  for (Int_t i = 0; i < fIsHeavy[iq]; i++){
     AliDebug(1,Form("Mother ID= %d type= %d label= %d\n",motherID[i],motherType[i],motherLabel[i]));
  }


  Int_t nheavypair = Int_t(fIsHeavy[iq]/2.); 
  for (Int_t ipair = 0; ipair < nheavypair; ipair++){

     Int_t id1 = ipair*2;
     Int_t id2 = ipair*2 + 1;

     if (motherType[id1] == 2 && motherType[id2] == 2){
       if (motherLabel[id1] == motherLabel[id2]) processID = kGluonSplitting; // gluon spliting
       else processID = -9;
     }
     else if (motherType[id1] == -1 && motherType[id2] == -1) {
       if (motherLabel[id1] == -1 && motherLabel[id2] == -1) {
         if (motherID[id1] == fgkGluon) processID = kPairCreationFromg; // gluon fusion
         else processID = kPairCreationFromq; // q-qbar pair creation
       }
       else processID = -8;
     }
     else if (motherType[id1] == -1 || motherType[id2] == -1) {
       if ((motherLabel[id1] == -1 || motherLabel[id2] == -1) && (motherLabel[id1]*motherLabel[id2] == -2 || motherLabel[id1]*motherLabel[id2] == -3)) {
         if(motherID[id1]*motherID[id2] == kquark*fgkGluon) processID = kFlavourExitation; // flavour exitation 
         else processID = kLightQuarkShower;
       }
       else processID = -7;
     }
     else if (motherType[id1] == -2 || motherType[id2] == -2) {
       if (motherLabel[id1] == motherLabel[id2]) processID = kInitialPartonShower; // initial parton shower
       else processID = -6;
       
     }
     else processID = -5;

     if (nheavypair >1) AliDebug(1,Form("Multi pair found : process ID = %d\n",processID));
     else fHistComm[iq][0].fProcessID->Fill(processID);
     AliDebug(1,Form("Process ID = %d\n",processID));
  } // end of # heavy quark pair loop

}

//__________________________________________
void AliHFEmcQA::GetHadronKine(TParticle* mcpart, const Int_t kquark)
{
    // decay electron kinematics

    if (kquark != kCharm && kquark != kBeauty) {
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return;
    }
    Int_t iq = kquark - kCharm;

    if(!mcpart){
      AliDebug(1, "no mc particle, return\n");
      return;
    }

    Int_t iLabel = mcpart->GetFirstMother();
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return;
    }

    if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

    TParticle *partCopy = mcpart;
    Int_t pdgcode = mcpart->GetPdgCode();
    Int_t pdgcodeCopy = pdgcode;

    AliMCParticle *mctrack = NULL;

    // if the mother is charmed hadron  
    Bool_t isDirectCharm = kFALSE;
    if ( int(TMath::Abs(pdgcode)/100.) == kCharm || int(TMath::Abs(pdgcode)/1000.) == kCharm ) {

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = mcpart->GetFirstMother();
             if (jLabel == -1){
               isDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return;
             }
             // if there is an ancester
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return; 
             TParticle* mother = mctrack->Particle();
             Int_t motherPDG = mother->GetPdgCode();
    
             for (Int_t j=0; j<fNparents; j++){
                if (TMath::Abs(motherPDG)==fParentSelect[1][j]) return; // return if this hadron is originated from b
             }

             mcpart = mother;
          } // end of iteration 
    } // end of if
    if((isDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (TMath::Abs(pdgcodeCopy)==fParentSelect[iq][i]){

              // fill hadron kinematics
              fHist[iq][kHadron][0][fCentrality].fPdgCode->Fill(pdgcodeCopy);
              fHist[iq][kHadron][0][fCentrality].fPt->Fill(partCopy->Pt());
              fHist[iq][kHadron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partCopy));
              fHist[iq][kHadron][0][fCentrality].fEta->Fill(partCopy->Eta());

              if(iq==0) {
               fhD[i]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
              }
            }
         }
	 // I also want to store D* info to compare with D* measurement 
	 if (TMath::Abs(pdgcodeCopy)==413 && iq==0) { //D*+
               fhD[7]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
	 }
	 if (TMath::Abs(pdgcodeCopy)==423 && iq==0) { //D*0
               fhD[8]->Fill(partCopy->Pt(),AliHFEtools::GetRapidity(partCopy));
	 }
    } // end of if
}

//__________________________________________
void AliHFEmcQA::GetDecayedKine(TParticle* mcpart, const Int_t kquark, Int_t kdecayed) 
{
    // decay electron kinematics
    
    if (!(kquark == kCharm || kquark == kBeauty || kquark == kOthers)){
      AliDebug(1, "This task is only for heavy quark QA, return\n");
      return; 
    }
    Int_t iq = kquark - kCharm; 
    Bool_t isFinalOpenCharm = kFALSE;

    if(!mcpart){
      AliDebug(1, "no mcparticle, return\n");
      return;
    }

    if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

    Double_t eabsEta = TMath::Abs(mcpart->Eta());
    Double_t eabsY = TMath::Abs(AliHFEtools::GetRapidity(mcpart));

    if(kquark==kOthers){
      Int_t esource = -1;
      if ( TMath::Abs(mcpart->GetPdgCode()) != kdecayed ) esource = kMisID-4;
      else esource =GetSource(mcpart)-4; // return for the cases kGamma=4, kPi0=5, kElse=6
      if(esource==0|| esource==1 || esource==2 || esource==3){
        fHist[iq][esource][0][fCentrality].fPt->Fill(mcpart->Pt());
        fHist[iq][esource][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
        fHist[iq][esource][0][fCentrality].fEta->Fill(mcpart->Eta());
        if(eabsEta<0.9){
          fHist[iq][esource][1][fCentrality].fPt->Fill(mcpart->Pt());
          fHist[iq][esource][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
          fHist[iq][esource][1][fCentrality].fEta->Fill(mcpart->Eta());
        }
        if(eabsY<0.5){
          fHist[iq][esource][2][fCentrality].fPt->Fill(mcpart->Pt());
          fHist[iq][esource][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
          fHist[iq][esource][2][fCentrality].fEta->Fill(mcpart->Eta());
        }
        return; 
      }
      else {
        AliDebug(1, "e source is out of defined ranges, return\n");
        return;
      }
    }

    if ( TMath::Abs(mcpart->GetPdgCode()) != kdecayed ) return;

    Int_t iLabel = mcpart->GetFirstMother(); 
    if (iLabel<0){
      AliDebug(1, "Stack label is negative, return\n");
      return; 
    }

    AliMCParticle *mctrack = NULL;
    if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return; 
    TParticle *partMother = mctrack->Particle();
    TParticle *partMotherCopy = partMother;
    Int_t maPdgcode = partMother->GetPdgCode();
    Int_t maPdgcodeCopy = maPdgcode;

    // get mc primary vertex
    /*
    TArrayF mcPrimVtx;
    if(fMCHeader) fMCHeader->PrimaryVertex(mcPrimVtx);

    // get electron production vertex   
    TLorentzVector ePoint;
    mcpart->ProductionVertex(ePoint);

    // calculated production vertex to primary vertex (in xy plane)
    Float_t decayLxy = TMath::Sqrt((mcPrimVtx[0]-ePoint[0])*(mcPrimVtx[0]-ePoint[0])+(mcPrimVtx[1]-ePoint[1])*(mcPrimVtx[1]-ePoint[1]));
    */ 
    Float_t decayLxy = 0;

    // if the mother is charmed hadron  
    Bool_t isMotherDirectCharm = kFALSE;
    if ( int(TMath::Abs(maPdgcode)/100.) == kCharm || int(TMath::Abs(maPdgcode)/1000.) == kCharm ) { 

         for (Int_t i=0; i<fNparents; i++){
            if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
              isFinalOpenCharm = kTRUE;
            } 
         }  
         if (!isFinalOpenCharm) return ;

          // iterate until you find B hadron as a mother or become top ancester 
          for (Int_t i=1; i<fgkMaxIter; i++){

             Int_t jLabel = partMother->GetFirstMother(); 
             if (jLabel == -1){
               isMotherDirectCharm = kTRUE;
               break; // if there is no ancester
             }
             if (jLabel < 0){ // safety protection
               AliDebug(1, "Stack label is negative, return\n");
               return; 
             }

             // if there is an ancester
             if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return; 
             TParticle* grandMa = mctrack->Particle();
             Int_t grandMaPDG = grandMa->GetPdgCode();

             for (Int_t j=0; j<fNparents; j++){
                if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){

                  if (kquark == kCharm) return;
                  // fill electron kinematics
                  fHist[iq][kElectron2nd][0][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
                  fHist[iq][kElectron2nd][0][fCentrality].fPt->Fill(mcpart->Pt());
                  fHist[iq][kElectron2nd][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                  fHist[iq][kElectron2nd][0][fCentrality].fEta->Fill(mcpart->Eta());

                  // fill mother hadron kinematics
                  fHist[iq][kDeHadron][0][fCentrality].fPdgCode->Fill(grandMaPDG); 
                  fHist[iq][kDeHadron][0][fCentrality].fPt->Fill(grandMa->Pt());
                  fHist[iq][kDeHadron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
                  fHist[iq][kDeHadron][0][fCentrality].fEta->Fill(grandMa->Eta());

                  if(eabsEta<0.9){
                    fHist[iq][kElectron2nd][1][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
                    fHist[iq][kElectron2nd][1][fCentrality].fPt->Fill(mcpart->Pt());
                    fHist[iq][kElectron2nd][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                    fHist[iq][kElectron2nd][1][fCentrality].fEta->Fill(mcpart->Eta());

                    // fill mother hadron kinematics
                    fHist[iq][kDeHadron][1][fCentrality].fPdgCode->Fill(grandMaPDG); 
                    fHist[iq][kDeHadron][1][fCentrality].fPt->Fill(grandMa->Pt());
                    fHist[iq][kDeHadron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
                    fHist[iq][kDeHadron][1][fCentrality].fEta->Fill(grandMa->Eta());
                  }

                  if(eabsY<0.5){
                    fHist[iq][kElectron2nd][2][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
                    fHist[iq][kElectron2nd][2][fCentrality].fPt->Fill(mcpart->Pt());
                    fHist[iq][kElectron2nd][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                    fHist[iq][kElectron2nd][2][fCentrality].fEta->Fill(mcpart->Eta());

                    // fill mother hadron kinematics
                    fHist[iq][kDeHadron][2][fCentrality].fPdgCode->Fill(grandMaPDG); 
                    fHist[iq][kDeHadron][2][fCentrality].fPt->Fill(grandMa->Pt());
                    fHist[iq][kDeHadron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
                    fHist[iq][kDeHadron][2][fCentrality].fEta->Fill(grandMa->Eta());
                  }

                  //mj: to calculate B to e eta correlation to calculate total heavy quark cross section
                  Int_t kLabel0 = grandMa->GetFirstMother();
                  Bool_t isGGrandmaYes = kFALSE;
                  Double_t ggmrapidwstmp=0;
                  if (!(kLabel0 < 0)){ // safety protection
                    if((mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(kLabel0))))){
                      TParticle* ggrandMatmp = mctrack->Particle();
                      Int_t ggrandMaPDGtmp = ggrandMatmp->GetPdgCode();
                      if ( int(TMath::Abs(ggrandMaPDGtmp)/100.) == kBeauty || int(TMath::Abs(ggrandMaPDGtmp)/1000.) == kBeauty) isGGrandmaYes = kTRUE;
                      ggmrapidwstmp = AliHFEtools::GetRapidity(ggrandMatmp);
                    }
                  }

                  Double_t gmrapidwstmp0 = AliHFEtools::GetRapidity(grandMa);
                  Double_t eetawstmp0 = mcpart->Eta();
  
                  Double_t gmrapidtmp0 = TMath::Abs(gmrapidwstmp0);
                  Double_t eetatmp0 = TMath::Abs(eetawstmp0);

                  fHistComm[iq][0].fEtaCorrB->Fill(gmrapidwstmp0,eetawstmp0);
                  if(isGGrandmaYes) fHistComm[iq][0].fEtaCorrGB->Fill(ggmrapidwstmp,eetawstmp0);
                  else fHistComm[iq][0].fEtaCorrGB->Fill(gmrapidwstmp0,eetawstmp0);

                  if(gmrapidtmp0<0.5 && eetatmp0<0.5 ) fHistComm[iq][0].fPtCorrBinein->Fill(grandMa->Pt(),mcpart->Pt());
                  else if(gmrapidtmp0<0.5 && eetatmp0>0.5 ) fHistComm[iq][0].fPtCorrBineout->Fill(grandMa->Pt(),mcpart->Pt());
                  else if(gmrapidtmp0>0.5 && eetatmp0<0.5 ) fHistComm[iq][0].fPtCorrBoutein->Fill(grandMa->Pt(),mcpart->Pt());
                  else if(gmrapidtmp0>0.5 && eetatmp0>0.5 ) fHistComm[iq][0].fPtCorrBouteout->Fill(grandMa->Pt(),mcpart->Pt());
                  //======================================================================================

                  // ratio between pT of electron and pT of mother B hadron 
                  if(grandMa->Pt()) {
                    fHistComm[iq][0].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());
                    if(eabsEta<0.9){
                      fHistComm[iq][1].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());
                    }
                    if(eabsY<0.5){
                      fHistComm[iq][2].fDePtRatio->Fill(grandMa->Pt(),mcpart->Pt()/grandMa->Pt());
                    }
                  }

                  // distance between electron production point and primary vertex
                  fHistComm[iq][0].fDeDistance->Fill(grandMa->Pt(),decayLxy);
                  if(eabsEta<0.9){
                    fHistComm[iq][1].fDeDistance->Fill(grandMa->Pt(),decayLxy);
                  }
                  if(eabsY<0.5){
                    fHistComm[iq][2].fDeDistance->Fill(grandMa->Pt(),decayLxy);
                  }
                  return;
                }
             } 

             partMother = grandMa;
          } // end of iteration 
    } // end of if
    if((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
         for (Int_t i=0; i<fNparents; i++){
            if (TMath::Abs(maPdgcodeCopy)==fParentSelect[iq][i]){

              fHist[iq][kElectron][0][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
              fHist[iq][kElectron][0][fCentrality].fPt->Fill(mcpart->Pt());
              fHist[iq][kElectron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
              fHist[iq][kElectron][0][fCentrality].fEta->Fill(mcpart->Eta());  

              // fill mother hadron kinematics
              fHist[iq][keHadron][0][fCentrality].fPdgCode->Fill(maPdgcodeCopy); 
              fHist[iq][keHadron][0][fCentrality].fPt->Fill(partMotherCopy->Pt());
              fHist[iq][keHadron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
              fHist[iq][keHadron][0][fCentrality].fEta->Fill(partMotherCopy->Eta());

              if(eabsEta<0.9){
                fHist[iq][kElectron][1][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
                fHist[iq][kElectron][1][fCentrality].fPt->Fill(mcpart->Pt());
                fHist[iq][kElectron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                fHist[iq][kElectron][1][fCentrality].fEta->Fill(mcpart->Eta());  

                // fill mother hadron kinematics
                fHist[iq][keHadron][1][fCentrality].fPdgCode->Fill(maPdgcodeCopy); 
                fHist[iq][keHadron][1][fCentrality].fPt->Fill(partMotherCopy->Pt());
                fHist[iq][keHadron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
                fHist[iq][keHadron][1][fCentrality].fEta->Fill(partMotherCopy->Eta());
              }

              if(eabsY<0.5){
                fHist[iq][kElectron][2][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
                fHist[iq][kElectron][2][fCentrality].fPt->Fill(mcpart->Pt());
                fHist[iq][kElectron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
                fHist[iq][kElectron][2][fCentrality].fEta->Fill(mcpart->Eta());  

                // fill mother hadron kinematics
                fHist[iq][keHadron][2][fCentrality].fPdgCode->Fill(maPdgcodeCopy); 
                fHist[iq][keHadron][2][fCentrality].fPt->Fill(partMotherCopy->Pt());
                fHist[iq][keHadron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
                fHist[iq][keHadron][2][fCentrality].fEta->Fill(partMotherCopy->Eta());
              }

              // ratio between pT of electron and pT of mother B or direct D hadron 
              if(partMotherCopy->Pt()) {
                 fHistComm[iq][0].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());
                 if(eabsEta<0.9){
                   fHistComm[iq][1].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());
                 }
                 if(eabsY<0.5){
                   fHistComm[iq][2].fePtRatio->Fill(partMotherCopy->Pt(),mcpart->Pt()/partMotherCopy->Pt());
                 }
              }
              fHistComm[iq][0].fPtCorr->Fill(partMotherCopy->Pt(),mcpart->Pt());
              if(eabsEta<0.9){
                fHistComm[iq][1].fPtCorr->Fill(partMotherCopy->Pt(),mcpart->Pt());
              }
              if(eabsY<0.5){
                fHistComm[iq][2].fPtCorr->Fill(partMotherCopy->Pt(),mcpart->Pt());
              }
              if(TMath::Abs(partMotherCopy->GetPdgCode())==411) {
                fHistComm[iq][0].fPtCorrDp->Fill(partMotherCopy->Pt(),mcpart->Pt());
                if(eabsEta<0.9){
                  fHistComm[iq][1].fPtCorrDp->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
                if(eabsY<0.5){
                  fHistComm[iq][2].fPtCorrDp->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
              }
              else if(TMath::Abs(partMotherCopy->GetPdgCode())==421) {
                fHistComm[iq][0].fPtCorrD0->Fill(partMotherCopy->Pt(),mcpart->Pt());
                if(eabsEta<0.9){
                  fHistComm[iq][1].fPtCorrD0->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
                if(eabsY<0.5){
                  fHistComm[iq][2].fPtCorrD0->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
              }
              else {
                fHistComm[iq][0].fPtCorrDrest->Fill(partMotherCopy->Pt(),mcpart->Pt());
                if(eabsEta<0.9){
                  fHistComm[iq][1].fPtCorrDrest->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
                if(eabsY<0.5){
                  fHistComm[iq][2].fPtCorrDrest->Fill(partMotherCopy->Pt(),mcpart->Pt());
                }
              }

              //mj: to calculate D to e eta correlation to calculate total heavy quark cross section
              Int_t kLabel = partMotherCopy->GetFirstMother();
              Bool_t isGrandmaYes = kFALSE;
              Double_t gmrapidwstmp =0;
              if (!(kLabel < 0)){ // safety protection
                if((mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(kLabel))))){ 
                  TParticle* grandMatmp = mctrack->Particle();
                  Int_t grandMaPDGtmp = grandMatmp->GetPdgCode();
                  if ( int(TMath::Abs(grandMaPDGtmp)/100.) == kCharm || int(TMath::Abs(grandMaPDGtmp)/1000.) == kCharm ) isGrandmaYes = kTRUE;
                  if ( int(TMath::Abs(grandMaPDGtmp)/100.) == kBeauty || int(TMath::Abs(grandMaPDGtmp)/1000.) == kBeauty) isGrandmaYes = kTRUE;
                  gmrapidwstmp = AliHFEtools::GetRapidity(grandMatmp);
                }
              }

              Double_t mrapidwstmp = AliHFEtools::GetRapidity(partMotherCopy);
              Double_t eetawstmp = mcpart->Eta();

              Double_t mrapidtmp = TMath::Abs(mrapidwstmp);
              Double_t eetatmp = TMath::Abs(eetawstmp);

              fHistComm[iq][0].fEtaCorrD->Fill(mrapidwstmp,eetawstmp);
              if(isGrandmaYes) fHistComm[iq][0].fEtaCorrGD->Fill(gmrapidwstmp,eetawstmp);
              else fHistComm[iq][0].fEtaCorrGD->Fill(mrapidwstmp,eetawstmp);

              if(mrapidtmp<0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDinein->Fill(partMotherCopy->Pt(),mcpart->Pt());
              else if(mrapidtmp<0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDineout->Fill(partMotherCopy->Pt(),mcpart->Pt());
              else if(mrapidtmp>0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDoutein->Fill(partMotherCopy->Pt(),mcpart->Pt());
              else if(mrapidtmp>0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDouteout->Fill(partMotherCopy->Pt(),mcpart->Pt());
              if(TMath::Abs(partMotherCopy->GetPdgCode())==411) {
                fHistComm[iq][0].fEtaCorrDp->Fill(mrapidwstmp,eetawstmp);
                if(isGrandmaYes) fHistComm[iq][0].fEtaCorrGDp->Fill(gmrapidwstmp,eetawstmp);
                else fHistComm[iq][0].fEtaCorrGDp->Fill(mrapidwstmp,eetawstmp);
                if(mrapidtmp<0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDpDinein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp<0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDpDineout->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDpDoutein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDpDouteout->Fill(partMotherCopy->Pt(),mcpart->Pt());
              }
              else if(TMath::Abs(partMotherCopy->GetPdgCode())==421) {
                fHistComm[iq][0].fEtaCorrD0->Fill(mrapidwstmp,eetawstmp);
                if(isGrandmaYes) fHistComm[iq][0].fEtaCorrGD0->Fill(gmrapidwstmp,eetawstmp);
                else fHistComm[iq][0].fEtaCorrGD0->Fill(mrapidwstmp,eetawstmp);
                if(mrapidtmp<0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrD0Dinein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp<0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrD0Dineout->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrD0Doutein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrD0Douteout->Fill(partMotherCopy->Pt(),mcpart->Pt());
              }
              else {
                fHistComm[iq][0].fEtaCorrDrest->Fill(mrapidwstmp,eetawstmp);
                if(isGrandmaYes) fHistComm[iq][0].fEtaCorrGDrest->Fill(gmrapidwstmp,eetawstmp);
                else fHistComm[iq][0].fEtaCorrGDrest->Fill(mrapidwstmp,eetawstmp);
                if(mrapidtmp<0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDrestDinein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp<0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDrestDineout->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp<0.5 ) fHistComm[iq][0].fPtCorrDrestDoutein->Fill(partMotherCopy->Pt(),mcpart->Pt());
                else if(mrapidtmp>0.5 && eetatmp>0.5 ) fHistComm[iq][0].fPtCorrDrestDouteout->Fill(partMotherCopy->Pt(),mcpart->Pt());
              }

              // distance between electron production point and primary vertex
              fHistComm[iq][0].feDistance->Fill(partMotherCopy->Pt(),decayLxy);
              if(eabsEta<0.9){
                fHistComm[iq][1].feDistance->Fill(partMotherCopy->Pt(),decayLxy);
              }
              if(eabsY<0.5){
                fHistComm[iq][2].feDistance->Fill(partMotherCopy->Pt(),decayLxy);
              }
            }
         }
    } // end of if
}

//____________________________________________________________________
void  AliHFEmcQA::GetDecayedKine(AliAODMCParticle *mcpart, const Int_t kquark, Int_t kdecayed)
{
  // decay electron kinematics

  if (kquark != kCharm && kquark != kBeauty) {
    AliDebug(1, "This task is only for heavy quark QA, return\n");
    return;
  }

  Int_t iq = kquark - kCharm;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return;
  }

  if ( TMath::Abs(mcpart->GetPdgCode()) != kdecayed ) return;

  Double_t eabsEta = TMath::Abs(mcpart->Eta());
  Double_t eabsY = TMath::Abs(AliHFEtools::GetRapidity(mcpart));

  // mother
  Int_t iLabel = mcpart->GetMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return;
  }

  if(!fIsPbPb&&!fIsppMultiBin) fCentrality=0;

  AliAODMCParticle *partMother = (AliAODMCParticle*)fMCArray->At(iLabel);
  AliAODMCParticle *partMotherCopy = partMother;
  Int_t maPdgcode = partMother->GetPdgCode();
  Int_t maPdgcodeCopy = maPdgcode;

  Bool_t isMotherDirectCharm = kFALSE;
  if ( int(TMath::Abs(maPdgcode)/100.) == kCharm || int(TMath::Abs(maPdgcode)/1000.) == kCharm ) {

    for (Int_t i=0; i<fNparents; i++){
       if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    } 
    if (!isFinalOpenCharm) return;

    for (Int_t i=1; i<fgkMaxIter; i++){

       Int_t jLabel = partMother->GetMother();
       if (jLabel == -1){
         isMotherDirectCharm = kTRUE;
         break; // if there is no ancester
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return;
       }

       // if there is an ancester
       AliAODMCParticle* grandMa = (AliAODMCParticle*)fMCArray->At(jLabel);
       Int_t grandMaPDG = grandMa->GetPdgCode();

       for (Int_t j=0; j<fNparents; j++){
          if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){

            if (kquark == kCharm) return;
            // fill electron kinematics
            fHist[iq][kElectron2nd][0][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
            fHist[iq][kElectron2nd][0][fCentrality].fPt->Fill(mcpart->Pt());
            fHist[iq][kElectron2nd][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
            fHist[iq][kElectron2nd][0][fCentrality].fEta->Fill(mcpart->Eta());

            // fill mother hadron kinematics
            fHist[iq][kDeHadron][0][fCentrality].fPdgCode->Fill(grandMaPDG);
            fHist[iq][kDeHadron][0][fCentrality].fPt->Fill(grandMa->Pt());
            fHist[iq][kDeHadron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
            fHist[iq][kDeHadron][0][fCentrality].fEta->Fill(grandMa->Eta());

            if(eabsEta<0.9){
              // fill electron kinematics
              fHist[iq][kElectron2nd][1][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
              fHist[iq][kElectron2nd][1][fCentrality].fPt->Fill(mcpart->Pt());
              fHist[iq][kElectron2nd][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
              fHist[iq][kElectron2nd][1][fCentrality].fEta->Fill(mcpart->Eta());

              // fill mother hadron kinematics
              fHist[iq][kDeHadron][1][fCentrality].fPdgCode->Fill(grandMaPDG);
              fHist[iq][kDeHadron][1][fCentrality].fPt->Fill(grandMa->Pt());
              fHist[iq][kDeHadron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
              fHist[iq][kDeHadron][1][fCentrality].fEta->Fill(grandMa->Eta());
            }
            if(eabsY<0.5){
              // fill electron kinematics
              fHist[iq][kElectron2nd][2][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
              fHist[iq][kElectron2nd][2][fCentrality].fPt->Fill(mcpart->Pt());
              fHist[iq][kElectron2nd][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
              fHist[iq][kElectron2nd][2][fCentrality].fEta->Fill(mcpart->Eta());

              // fill mother hadron kinematics
              fHist[iq][kDeHadron][2][fCentrality].fPdgCode->Fill(grandMaPDG);
              fHist[iq][kDeHadron][2][fCentrality].fPt->Fill(grandMa->Pt());
              fHist[iq][kDeHadron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(grandMa));
              fHist[iq][kDeHadron][2][fCentrality].fEta->Fill(grandMa->Eta());
            }

            return;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  if ((isMotherDirectCharm == kTRUE && kquark == kCharm) || kquark == kBeauty) {
    for (Int_t i=0; i<fNparents; i++){
       if (TMath::Abs(maPdgcodeCopy)==fParentSelect[iq][i]){

         // fill electron kinematics
         fHist[iq][kElectron][0][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
         fHist[iq][kElectron][0][fCentrality].fPt->Fill(mcpart->Pt());
         fHist[iq][kElectron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
         fHist[iq][kElectron][0][fCentrality].fEta->Fill(mcpart->Eta());

         // fill mother hadron kinematics
         fHist[iq][keHadron][0][fCentrality].fPdgCode->Fill(maPdgcodeCopy);
         fHist[iq][keHadron][0][fCentrality].fPt->Fill(partMotherCopy->Pt());
         fHist[iq][keHadron][0][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
         fHist[iq][keHadron][0][fCentrality].fEta->Fill(partMotherCopy->Eta());

         if(eabsEta<0.9){
           // fill electron kinematics
           fHist[iq][kElectron][1][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
           fHist[iq][kElectron][1][fCentrality].fPt->Fill(mcpart->Pt());
           fHist[iq][kElectron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
           fHist[iq][kElectron][1][fCentrality].fEta->Fill(mcpart->Eta());

           // fill mother hadron kinematics
           fHist[iq][keHadron][1][fCentrality].fPdgCode->Fill(maPdgcodeCopy);
           fHist[iq][keHadron][1][fCentrality].fPt->Fill(partMotherCopy->Pt());
           fHist[iq][keHadron][1][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
           fHist[iq][keHadron][1][fCentrality].fEta->Fill(partMotherCopy->Eta());
         }
         if(eabsY<0.5){
           // fill electron kinematics
           fHist[iq][kElectron][2][fCentrality].fPdgCode->Fill(mcpart->GetPdgCode());
           fHist[iq][kElectron][2][fCentrality].fPt->Fill(mcpart->Pt());
           fHist[iq][kElectron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(mcpart));
           fHist[iq][kElectron][2][fCentrality].fEta->Fill(mcpart->Eta());

           // fill mother hadron kinematics
           fHist[iq][keHadron][2][fCentrality].fPdgCode->Fill(maPdgcodeCopy);
           fHist[iq][keHadron][2][fCentrality].fPt->Fill(partMotherCopy->Pt());
           fHist[iq][keHadron][2][fCentrality].fY->Fill(AliHFEtools::GetRapidity(partMotherCopy));
           fHist[iq][keHadron][2][fCentrality].fEta->Fill(partMotherCopy->Eta());
         }

       }
    }
  } // end of if

}

//__________________________________________
void AliHFEmcQA::IdentifyMother(Int_t motherlabel, Int_t &motherpdg, Int_t &grandmotherlabel)
{
       // find mother pdg code and label 

       if (motherlabel < 0) { 
         AliDebug(1, "Stack label is negative, return\n");
         return; 
       }
       AliMCParticle *mctrack = NULL;
       if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(motherlabel))))) return; 
       TParticle *heavysMother = mctrack->Particle();
       motherpdg = heavysMother->GetPdgCode();
       grandmotherlabel = heavysMother->GetFirstMother();
       AliDebug(1,Form("ancestor pdg code= %d\n",motherpdg));
}

//__________________________________________
void AliHFEmcQA::HardScattering(const Int_t kquark, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -1 means this heavy quark coming from hard vertex

       AliMCParticle *mctrack1 = NULL;
       AliMCParticle *mctrack2 = NULL;
       if(!(mctrack1 = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(4))))) return; 
       if(!(mctrack2 = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(5))))) return; 
       TParticle *afterinitialrad1  = mctrack1->Particle();
       TParticle *afterinitialrad2  = mctrack2->Particle();
           
       motherlabel = -1;

       if (TMath::Abs(afterinitialrad1->GetPdgCode()) == fgkGluon && TMath::Abs(afterinitialrad2->GetPdgCode()) == fgkGluon){
         AliDebug(1,"heavy from gluon gluon pair creation!\n");
         mothertype = -1;
         motherID = fgkGluon;
       }
       else if (TMath::Abs(afterinitialrad1->GetPdgCode()) == kquark || TMath::Abs(afterinitialrad2->GetPdgCode()) == kquark){ // one from Q and the other from g
         AliDebug(1,"heavy from flavor exitation!\n");
         mothertype = -1;
         motherID = kquark;
       }
       else if  (TMath::Abs(afterinitialrad1->GetPdgCode()) == TMath::Abs(afterinitialrad2->GetPdgCode())){
         AliDebug(1,"heavy from q-qbar pair creation!\n");
         mothertype = -1;
         motherID = 1;
       }
       else {
         AliDebug(1,"something strange!\n");
         mothertype = -999;
         motherlabel = -999;
         motherID = -999;
       }
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromInitialShower(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype -2 means this heavy quark coming from initial state 

       AliMCParticle *mctrack = NULL;
       if (inputmotherlabel==2 || inputmotherlabel==3){ // mother exist before initial state radiation
         if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(inputmotherlabel))))) return kFALSE; 
         TParticle *heavysMother = mctrack->Particle();
         motherID = heavysMother->GetPdgCode(); 
         mothertype = -2; // there is mother before initial state radiation
         motherlabel = inputmotherlabel;
         AliDebug(1,"initial parton shower! \n");

         return kTRUE;
       }

       return kFALSE;
}

//__________________________________________
Bool_t AliHFEmcQA::IsFromFinalParton(Int_t inputmotherlabel, Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
       // mothertype 2 means this heavy quark coming from final state 

       AliMCParticle *mctrack = NULL;
       if (inputmotherlabel > 5){ // mother exist after hard scattering
         if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(inputmotherlabel))))) return kFALSE; 
         TParticle *heavysMother = mctrack->Particle();
         motherID = heavysMother->GetPdgCode(); 
         mothertype = 2; // 
         motherlabel = inputmotherlabel;
         AliDebug(1,Form("heavy quark from %d after hard scattering! \n",motherID));

         return kTRUE;
       }
       return kFALSE;
}

//__________________________________________
void AliHFEmcQA::ReportStrangeness(Int_t &motherID, Int_t &mothertype, Int_t &motherlabel)
{
      // mark strange behavior  

       mothertype = -888;
       motherlabel = -888;
       motherID = -888;
       AliDebug(1,"something strange!\n");
}

//__________________________________________
Int_t AliHFEmcQA::GetSource(const AliVParticle* const mcpart) const
{        
  // decay particle's origin 

  //if ( TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;
       
  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "Stack label is negative or no mcparticle, return\n");
    return -1;
  }

  // mother
  // Information not in the base class, cast necessary
  Int_t iLabel = GetMother(mcpart);
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  } 
       
  const AliVParticle *partMother = fMCEvent->GetTrack(iLabel);
  Int_t maPdgcode = partMother->PdgCode();
  
  // if the mother is charmed hadron  
  if ( int(TMath::Abs(maPdgcode)/100.) == kCharm || int(TMath::Abs(maPdgcode)/1000.) == kCharm ) {
    
    for (Int_t i=0; i<fNparents; i++){
       if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
         isFinalOpenCharm = kTRUE;
       }
    }
    if (!isFinalOpenCharm) return -1;

    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<fgkMaxIter; i++){
      
       Int_t jLabel = GetMother(partMother);
       if (jLabel == -1){
         origin = kDirectCharm;
         return origin;
       }
       if (jLabel < 0){ // safety protection
         AliDebug(1, "Stack label is negative, return\n");
         return -1;
       }

       // if there is an ancester
       const AliVParticle* grandMa = fMCEvent->GetTrack(jLabel);
       Int_t grandMaPDG = grandMa->PdgCode();

       for (Int_t j=0; j<fNparents; j++){
          if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
            origin = kBeautyCharm;
            return origin;
          }
       }

       partMother = grandMa;
    } // end of iteration 
  } // end of if
  else if ( int(TMath::Abs(maPdgcode)/100.) == kBeauty || int(TMath::Abs(maPdgcode)/1000.) == kBeauty ) {
    for (Int_t i=0; i<fNparents; i++){
       if (TMath::Abs(maPdgcode)==fParentSelect[1][i]){
         origin = kDirectBeauty;
         return origin;
       }
    }
  } // end of if
  else if ( TMath::Abs(maPdgcode) == 22 ) {
    origin = kGamma;
    return origin;
  } // end of if
  else if ( TMath::Abs(maPdgcode) == 111 ) {
    origin = kPi0;
    return origin;
  } // end of if

  return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetSource(const TParticle * const mcpart) const
{
  // decay particle's origin 

  //if ( TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return -1;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  AliMCParticle *mctrack = NULL;
  if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return -1; 
  TParticle *partMother = mctrack->Particle();
  Int_t maPdgcode = partMother->GetPdgCode();

   // if the mother is charmed hadron  
   if ( int(TMath::Abs(maPdgcode)/100.) == kCharm || int(TMath::Abs(maPdgcode)/1000.) == kCharm ) {

     for (Int_t i=0; i<fNparents; i++){
        if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
          isFinalOpenCharm = kTRUE;
        }
     }
     if (!isFinalOpenCharm) return -1;

     // iterate until you find B hadron as a mother or become top ancester 
     for (Int_t i=1; i<fgkMaxIter; i++){

        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          origin = kDirectCharm;
          return origin;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return -1; 
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t j=0; j<fNparents; j++){
           if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
             origin = kBeautyCharm;
             return origin;
           }
        }

        partMother = grandMa;
     } // end of iteration 
   } // end of if
   else if ( int(TMath::Abs(maPdgcode)/100.) == kBeauty || int(TMath::Abs(maPdgcode)/1000.) == kBeauty ) {
     for (Int_t i=0; i<fNparents; i++){
        if (TMath::Abs(maPdgcode)==fParentSelect[1][i]){
          origin = kDirectBeauty;
          return origin;
        }
     }
   } // end of if
   else if ( TMath::Abs(maPdgcode) == 22 ) {
     origin = kGamma;
     return origin;
   } // end of if
   else if ( TMath::Abs(maPdgcode) == 111 ) {
     origin = kPi0;
     return origin;
   } // end of if
   else origin = kElse;

   return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetElecSource(TParticle * const mcpart, Bool_t isElec) const
{
  Double_t mpt;
  return GetElecSource(mcpart, isElec, mpt);
}

//__________________________________________
Int_t AliHFEmcQA::GetElecSource(TParticle * const mcpart, Bool_t isElec, Double_t &mpt) const
{
  // decay particle's origin 

  if(!mcpart){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  if(isElec) if ( TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return kMisID;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  Int_t iLabel = mcpart->GetFirstMother();
  if (iLabel<0){
    AliDebug(1, "Stack label is negative, return\n");
    return -1;
  }

  AliMCParticle *mctrack = NULL;
  Int_t tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(iLabel))))) return -1; 
  TParticle *partMother = mctrack->Particle();
  TParticle *partMotherCopy = mctrack->Particle();
  Int_t maPdgcode = partMother->GetPdgCode();
  mpt = partMother->Pt();
  Int_t grmaPdgcode = 0;
  Int_t ggrmaPdgcode;
  Double_t gmpt, ggmpt;
   // if the mother is charmed hadron  

   if(TMath::Abs(maPdgcode)==443){ // J/spi
      Int_t jLabel = partMother->GetFirstMother();
      if((jLabel>=0) && (mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))){
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();
        mpt = grandMa->Pt();
        if((int(TMath::Abs(grandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(grandMaPDG)/1000.)%10) == kBeauty) return kB2Jpsi;
      }
      return kJpsi;   
   } 
   else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(maPdgcode)/1000.)%10) == kCharm ) {

     for (Int_t i=0; i<fNparents; i++){
        if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
	  mpt = partMother->Pt();
          isFinalOpenCharm = kTRUE;
        }
     }
     if (!isFinalOpenCharm) return -1;

     // iterate until you find B hadron as a mother or become top ancester 
     for (Int_t i=1; i<fgkMaxIter; i++){

        
        Int_t jLabel = partMother->GetFirstMother();
        if (jLabel == -1){
          return kDirectCharm;
        }
        if (jLabel < 0){ // safety protection
          AliDebug(1, "Stack label is negative, return\n");
          return -1;
        }

        // if there is an ancester
        if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(jLabel))))) return -1; 
        TParticle* grandMa = mctrack->Particle();
        Int_t grandMaPDG = grandMa->GetPdgCode();

        for (Int_t j=0; j<fNparents; j++){
           if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
             mpt = grandMa->Pt();
             return kBeautyCharm;
           }
        }

        partMother = grandMa;
     } // end of iteration 
   } // end of if
   else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(maPdgcode)/1000.)%10) == kBeauty ) {
     for (Int_t i=0; i<fNparents; i++){
        if (TMath::Abs(maPdgcode)==fParentSelect[1][i]){
          mpt = partMotherCopy->Pt();
          return kDirectBeauty;
        }
     }
   } // end of if
   else if ( TMath::Abs(maPdgcode) == 22 ) { //conversion

     tmpMomLabel = partMotherCopy->GetFirstMother();//mother of photon
     mpt = partMotherCopy->Pt();//pT of photon
     if(tmpMomLabel==-1) return kGamma; // no grand mother
     if(tmpMomLabel<0) return -1;
     if(!(mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(tmpMomLabel)))) return -1;
     partMother = mctrack->Particle();//partMother = photon's mom
     partMotherCopy = mctrack->Particle();//partMother = photon's mom
     mpt = partMother->Pt(); //pt of photon's mom
     maPdgcode = partMother->GetPdgCode();//pdg of photon's mom

     // check if the ligth meson is the decay product of heavy mesons
     tmpMomLabel = partMother->GetFirstMother();
     if((tmpMomLabel>=0) && (mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(tmpMomLabel)))) { // grandgrandmother
      partMother = mctrack->Particle();
      grmaPdgcode = partMother->GetPdgCode();
      mpt = partMother->Pt();
      gmpt = partMother->Pt();

      if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) return kGammaB2M;
      if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) return kGammaD2M;

      tmpMomLabel = partMother->GetFirstMother();
      if((tmpMomLabel>=0) && (mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(tmpMomLabel)))) { // grandgrandgrandmother
       partMother = mctrack->Particle();
       ggrmaPdgcode = partMother->GetPdgCode();
       mpt = partMother->Pt();
       ggmpt = partMother->Pt();

       if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) return kGammaB2M;
       if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) return kGammaD2M;
      } // grandgrandgrandmother

      if ( TMath::Abs(maPdgcode) == 111 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
        else if(grmaPdgcode == 310) return kGammaK0s2P;
        else if(grmaPdgcode == 130) return kGammaK0l2P;
        else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
        else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
        else if(grmaPdgcode == 3222) return kGammaSigma2P;
        mpt = partMotherCopy->Pt();
        return kGammaPi0;
      } 
      else if ( TMath::Abs(maPdgcode) == 221 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
        mpt = partMotherCopy->Pt();
        return kGammaEta;
      } 
      else if ( TMath::Abs(maPdgcode) == 223 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
        mpt = partMotherCopy->Pt();
        return kGammaOmega;
      } 
      else if ( TMath::Abs(maPdgcode) == 333 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
        mpt = partMotherCopy->Pt();
        return kGammaPhi;
      }
      else if ( TMath::Abs(maPdgcode) == 331 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
        mpt = partMotherCopy->Pt();
        return kGammaEtaPrime; 
      }
      else if ( TMath::Abs(maPdgcode) == 113 ) {
        mpt = gmpt; 
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
        mpt = partMotherCopy->Pt();
        return kGammaRho0;
      }
      else {
	origin = kElse; // grandgrandmother there but nothing we identify
      }
     } // grandgrandmother
     else { // no grandgrandmother

       if ( TMath::Abs(maPdgcode) == 111 ) {
	 return kGammaPi0;
       } 
       else if ( TMath::Abs(maPdgcode) == 221 ) {
	 return kGammaEta;
       } 
       else if ( TMath::Abs(maPdgcode) == 223 ) {
	 return kGammaOmega;
       } 
       else if ( TMath::Abs(maPdgcode) == 333 ) {
	 return kGammaPhi;
       }
       else if ( TMath::Abs(maPdgcode) == 331 ) {
	 return kGammaEtaPrime; 
       }
       else if ( TMath::Abs(maPdgcode) == 113 ) {
	 return kGammaRho0;
       }
       else origin = kElse;//grandmother is primary but nothing we identify
     }
     return origin;

   } 
   else {

     // check if the ligth meson is the decay product of heavy mesons
     tmpMomLabel = partMotherCopy->GetFirstMother();
     mpt = partMotherCopy->Pt();
     if((tmpMomLabel>=0) && (mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tmpMomLabel))))) {//grandmother
      partMother = mctrack->Particle();
      grmaPdgcode = partMother->GetPdgCode();
      mpt = partMother->Pt();
      gmpt = partMother->Pt();

      if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) return kB2M;
      if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) return kD2M;

      tmpMomLabel = partMother->GetFirstMother();
      if((tmpMomLabel>=0) && (mctrack = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(TMath::Abs(tmpMomLabel))))) { // grandgrandmother
       partMother = mctrack->Particle();
       ggrmaPdgcode = partMother->GetPdgCode();
       mpt = partMother->Pt();
       ggmpt = partMother->Pt();

       if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) return kB2M;
       if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) return kD2M;
      } //grandgrandmother 

      if ( TMath::Abs(maPdgcode) == 111 ) {
        mpt = gmpt;
        if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
        else if(grmaPdgcode == 310) return kK0s2P;
        else if(grmaPdgcode == 130) return kK0l2P;
        else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
        else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
        else if(grmaPdgcode == 3222) return kSigma2P;
        mpt = partMotherCopy->Pt();
        return kPi0;
      } 
      else if ( TMath::Abs(maPdgcode) == 221 ) {
        mpt = gmpt;
        if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
        mpt = partMotherCopy->Pt();
        return kEta;
      } 
      else if ( TMath::Abs(maPdgcode) == 223 ) {
        mpt = gmpt;
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
        mpt = partMotherCopy->Pt();
        return kOmega;
      } 
      else if ( TMath::Abs(maPdgcode) == 333 ) {
        mpt = gmpt;
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
        mpt = partMotherCopy->Pt();
        return kPhi;
      } 
      else if ( TMath::Abs(maPdgcode) == 331 ) {
        mpt = gmpt;
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
        mpt = partMotherCopy->Pt();
        return kEtaPrime;
      } 
      else if ( TMath::Abs(maPdgcode) == 113 ) {
        mpt = gmpt;
        if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
        mpt = partMotherCopy->Pt();
        return kRho0;
      } 
      else if ( TMath::Abs(maPdgcode) == 321 ) {
        mpt = partMotherCopy->Pt();
        return kKe3;
      }
      else if ( TMath::Abs(maPdgcode) == 130 ) {
        mpt = partMotherCopy->Pt();
        return kK0L;
      }
      else origin = kElse; //grandmother exist but nothing we identify
     }
     else {// mother is primary
       
       if ( TMath::Abs(maPdgcode) == 111 ) {
	 return kPi0;
       } 
       else if ( TMath::Abs(maPdgcode) == 221 ) {
	 return kEta;
       } 
       else if ( TMath::Abs(maPdgcode) == 223 ) {
	 return kOmega;
       } 
       else if ( TMath::Abs(maPdgcode) == 333 ) {
	 return kPhi;
       } 
       else if ( TMath::Abs(maPdgcode) == 331 ) {
	 return kEtaPrime;
       } 
       else if ( TMath::Abs(maPdgcode) == 113 ) {
	 return kRho0;
       } 
       else if ( TMath::Abs(maPdgcode) == 321 ) {
	 return kKe3;
       }
       else if ( TMath::Abs(maPdgcode) == 130 ) {
	 return kK0L;
       }else origin = kElse;
     }

   }//mother is something different from J/psi,charm,beauty or gamma
   return origin;
}

//__________________________________________
Int_t AliHFEmcQA::GetElecSource(const AliVParticle * const mctrack, Bool_t isElec) const
{
  Double_t mpt;
  return GetElecSource(mctrack, isElec, mpt);
}

//__________________________________________
Int_t AliHFEmcQA::GetElecSource(const AliVParticle * const mctrack, Bool_t isElec, Double_t &mpt) const
{
  //
  // decay particle's origin 
  //

  if(!mctrack){
    AliDebug(1, "no mcparticle, return\n");
    return -1;
  }

  TClass *type = mctrack->IsA();  

  if(type == AliMCParticle::Class()) {
    const AliMCParticle *esdmc = dynamic_cast<const AliMCParticle *>(mctrack);
    if(esdmc){
      TParticle *mcpart =  esdmc->Particle();
      return GetElecSource(mcpart, isElec, mpt);
    }
    else return -1;
  }
  if(type == AliAODMCParticle::Class()) {
    const AliAODMCParticle *aodmc = dynamic_cast<const AliAODMCParticle *>(mctrack);
    if(aodmc){
      return GetElecSource(aodmc, isElec);
    }
    else return -1;
  }
  return -1;
}
//__________________________________________
Int_t AliHFEmcQA::GetElecSource(const AliAODMCParticle * const mcpart, Bool_t isElec) const
{
  //
  // Function for AliAODMCParticle
  //

  if (!mcpart) return -1;
  if (!fMCArray) return -1;
  
  if(isElec) if ( TMath::Abs(mcpart->GetPdgCode()) != AliHFEmcQA::kElectronPDG ) return kMisID;

  Int_t origin = -1;
  Bool_t isFinalOpenCharm = kFALSE;

  Int_t iLabel = mcpart->GetMother();
  if ((iLabel<0) || (iLabel>=fMCArray->GetEntriesFast())){
    AliDebug(1, "label is out of range, return\n");
    return -1;
  }

  AliAODMCParticle *mctrack = NULL; // will change all the time
  Int_t tmpMomLabel=0;
  if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return -1; 
  AliAODMCParticle *partMother = mctrack; 
  AliAODMCParticle *partMotherCopy = mctrack;
  Int_t maPdgcode = mctrack->GetPdgCode();
  Int_t grmaPdgcode;
  Int_t ggrmaPdgcode;

  // if the mother is charmed hadron  

  if(TMath::Abs(maPdgcode)==443){ 
    //
    // J/spi
    //
    Int_t jLabel = partMother->GetMother();
    if ((jLabel>=0) && (jLabel<fMCArray->GetEntriesFast())){
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(jLabel))))){
	Int_t grandMaPDG = mctrack->GetPdgCode();
	if((int(TMath::Abs(grandMaPDG)/100.)%10) == kBeauty || (int(TMath::Abs(grandMaPDG)/1000.)%10) == kBeauty) {
	  return kB2Jpsi;
	}
      }
    }
    return kJpsi;   
  } 
  else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(maPdgcode)/1000.)%10) == kCharm ) {
    //
    // charm
    //
    for (Int_t i=0; i<fNparents; i++){
      if (TMath::Abs(maPdgcode)==fParentSelect[0][i]){
	isFinalOpenCharm = kTRUE;
      }
    }
    if (!isFinalOpenCharm) {
      return -1;
    }
    
    // iterate until you find B hadron as a mother or become top ancester 
    for (Int_t i=1; i<fgkMaxIter; i++){
      
      Int_t jLabel = partMother->GetMother();
      if (jLabel == -1){
	return kDirectCharm;
      }
      if ((jLabel<0) || (jLabel>=fMCArray->GetEntriesFast())){
	AliDebug(1, "Stack label is negative, return\n");
	return -1;
      }
      
      // if there is an ancester
      if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(jLabel))))) {
	return -1; 
      }
      Int_t grandMaPDG = mctrack->GetPdgCode();
      for (Int_t j=0; j<fNparents; j++){
	if (TMath::Abs(grandMaPDG)==fParentSelect[1][j]){
	  return kBeautyCharm;
	}
      }
      partMother = mctrack;
    } // end of iteration 

  } // end of if
  else if ( (int(TMath::Abs(maPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(maPdgcode)/1000.)%10) == kBeauty ) {
    //
    // beauty
    //
    for (Int_t i=0; i<fNparents; i++){
      if (TMath::Abs(maPdgcode)==fParentSelect[1][i]){
	return kDirectBeauty;
      }
    }
  } // end of if
  else if ( TMath::Abs(maPdgcode) == 22 ) { 
    //
    //conversion
    //
    tmpMomLabel = partMotherCopy->GetMother();
    if(tmpMomLabel==-1) return kGamma;
    if((tmpMomLabel<0) || (tmpMomLabel>=fMCArray->GetEntriesFast())) {
      return -1;
    }
    if(!(mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
      return -1;
    }
    partMother = mctrack;
    maPdgcode = partMother->GetPdgCode();
    
    // check if the ligth meson is the decay product of heavy mesons
    tmpMomLabel = partMother->GetMother();
    if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandmother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
	partMother = mctrack;
	grmaPdgcode = partMother->GetPdgCode();
	
	if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) {	 
	  return kGammaB2M;
	}
	if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) {	 
	  return kGammaD2M;
	}
	
	tmpMomLabel = partMother->GetMother();
	if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandgrandmother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
	    partMother = mctrack;
	    ggrmaPdgcode = partMother->GetPdgCode();
	    
	    if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) {
	      return kGammaB2M;
	    }
	    if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) {
	      return kGammaD2M;
	    }
	  }
	}//grandgrandgrandmother
	
	if ( TMath::Abs(maPdgcode) == 111 ) {
	  if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
	    else if(grmaPdgcode == 310) return kGammaK0s2P;
        else if(grmaPdgcode == 130) return kGammaK0l2P;
        else if(TMath::Abs(grmaPdgcode) == 321) return kGammaK2P;
        else if(TMath::Abs(grmaPdgcode) == 3122) return kGammaLamda2P;
        else if(grmaPdgcode == 3222) return kGammaSigma2P;
	  return kGammaPi0;
	} 
	else if ( TMath::Abs(maPdgcode) == 221 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
	  return kGammaEta;
	} 
	else if ( TMath::Abs(maPdgcode) == 223 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
	  return kGammaOmega;
	} 
	else if ( TMath::Abs(maPdgcode) == 333 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kGammaM2M;
	  return kGammaPhi;
	}
	else if ( TMath::Abs(maPdgcode) == 331 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kGammaM2M;
	  return kGammaEtaPrime; 
	}
	else if ( TMath::Abs(maPdgcode) == 113 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kGammaM2M;
	  return kGammaRho0;
	}
	else origin = kElse;//grandgrandmother but nothing we identify
      }//mctrack grandgrandmother
    }
    else {
      // grandmother is primary
      if ( TMath::Abs(maPdgcode) == 111 ) {
	return kGammaPi0;
      } 
      else if ( TMath::Abs(maPdgcode) == 221 ) {
	return kGammaEta;
      } 
      else if ( TMath::Abs(maPdgcode) == 223 ) {
	return kGammaOmega;
      } 
      else if ( TMath::Abs(maPdgcode) == 333 ) {
	return kGammaPhi;
      }
      else if ( TMath::Abs(maPdgcode) == 331 ) {
	return kGammaEtaPrime; 
      }
      else if ( TMath::Abs(maPdgcode) == 113 ) {
	return kGammaRho0;
      }
      else origin = kElse;//grandmother is primary but nothing we identify
    }
    
    return origin;
    
  } 
  else {
    //
    // check if the ligth meson is the decay product of heavy mesons
    //
    tmpMomLabel = partMotherCopy->GetMother();
    if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandmother
      if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
	partMother = mctrack;
	grmaPdgcode = partMother->GetPdgCode();
	
	if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kBeauty ) {
	  return kB2M;
	}
	if ( (int(TMath::Abs(grmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(grmaPdgcode)/1000.)%10) == kCharm ) {
	  return kD2M;
	}
	
	tmpMomLabel = partMother->GetMother();
	if((tmpMomLabel>=0) && (tmpMomLabel<fMCArray->GetEntriesFast())) {//grandgrandmother
	  if((mctrack = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(tmpMomLabel))))) {
	    partMother = mctrack;
	    ggrmaPdgcode = partMother->GetPdgCode();
	    
	    if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kBeauty || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kBeauty ) {
	      return kB2M;
	    }
	    if ( (int(TMath::Abs(ggrmaPdgcode)/100.)%10) == kCharm || (int(TMath::Abs(ggrmaPdgcode)/1000.)%10) == kCharm ) {
	      return kD2M;
	    }
	  }
	}//grandgrandmother
	
	if ( TMath::Abs(maPdgcode) == 111 ) {
	  if(grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
	   else if(grmaPdgcode == 310) return kK0s2P;
       else if(grmaPdgcode == 130) return kK0l2P;
       else if(TMath::Abs(grmaPdgcode) == 321) return kK2P;
       else if(TMath::Abs(grmaPdgcode) == 3122) return kLamda2P;
       else if(grmaPdgcode == 3222) return kSigma2P;
	  return kPi0;
	} 
	else if ( TMath::Abs(maPdgcode) == 221 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
	  return kEta;
	} 
	else if ( TMath::Abs(maPdgcode) == 223 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 333 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
	  return kOmega;
	} 
	else if ( TMath::Abs(maPdgcode) == 333 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 331 || grmaPdgcode == 113) return kM2M;
	  return kPhi;
	} 
	else if ( TMath::Abs(maPdgcode) == 331 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 113) return kM2M;
	  return kEtaPrime;
	} 
	else if ( TMath::Abs(maPdgcode) == 113 ) {
	  if(grmaPdgcode == 111 || grmaPdgcode == 221 || grmaPdgcode == 223 || grmaPdgcode == 333 || grmaPdgcode == 331) return kM2M;
	  return kRho0;
	} 
	else if ( TMath::Abs(maPdgcode) == 321 ) {
	  return kKe3;
	}
	else if ( TMath::Abs(maPdgcode) == 130 ) {
	  return kK0L;
	}
	else origin = kElse;//grandmother but nothing we identidy
      }//mctrack grandmother
    }
    else {
      // no grandmother
      if ( TMath::Abs(maPdgcode) == 111 ) {
	return kPi0;
      } 
      else if ( TMath::Abs(maPdgcode) == 221 ) {
	return kEta;
      } 
      else if ( TMath::Abs(maPdgcode) == 223 ) {
	return kOmega;
      } 
      else if ( TMath::Abs(maPdgcode) == 333 ) {
	return kPhi;
      } 
      else if ( TMath::Abs(maPdgcode) == 331 ) {
	return kEtaPrime;
      } 
      else if ( TMath::Abs(maPdgcode) == 113 ) {
	return kRho0;
      } 
      else if ( TMath::Abs(maPdgcode) == 321 ) {
	return kKe3;
      }
      else if ( TMath::Abs(maPdgcode) == 130 ) {
	return kK0L;
      }
      else origin = kElse;//mother but nothing we identify
    }
  }//mother is something different from J/psi,charm,beauty or gamma
  return origin;
}
//__________________________________________
Double_t AliHFEmcQA::GetWeightFactor(AliMCParticle *mctrack, const Int_t iBgLevel){
  //
  // Get weighting factor for the realistic background estimation, for three possible background yield levels, indicated by the argument "iLevel": the best estimate (0), the lower uncertainty level (1), and the upper uncertainty level (2)
  //
  AliMCParticle *mctrackmother = NULL;  
  Double_t weightElecBg = 0.;
  Double_t mesonPt = 0.;
  Double_t mesonMotherPt = 0.;
  Double_t bgcategory = 0.;
  //Bool_t condition = kTRUE;
  Int_t mArr = -1;  
  Int_t mesonID = GetElecSource(mctrack->Particle(), kTRUE);
  if(mesonID==kGammaPi0 || mesonID==kPi0) mArr=0;                //pion
  else if(mesonID==kGammaEta || mesonID==kEta) mArr=1;           //eta
  else if(mesonID==kGammaOmega || mesonID==kOmega) mArr=2;       //omega
  else if(mesonID==kGammaPhi || mesonID==kPhi) mArr=3;           //phi
  else if(mesonID==kGammaEtaPrime || mesonID==kEtaPrime) mArr=4; //etaprime
  else if(mesonID==kGammaRho0 || mesonID==kRho0) mArr=5;         //rho
  else if(mesonID==kKe3 || mesonID==kGammaK2P|| mesonID==kK2P) mArr=6; //ke3 or K->pi->e
  else if(mesonID==kK0L || mesonID==kGammaK0s2P|| mesonID==kK0s2P) mArr=7; //K0L->e+X or k0s->pi->e
  else if(mesonID==kGammaLamda2P|| mesonID==kLamda2P) mArr=8;    //lambda->pi->e

  Double_t datamc[30]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999,-999};
  Double_t xr[3]={-999,-999,-999};
  datamc[0] = mesonID;
  datamc[17] = mctrack->Pt(); //electron pt
  datamc[18] = mctrack->Eta(); //electron eta

  mctrack->XvYvZv(xr);
  datamc[9] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
  datamc[10] = xr[2];

  TParticle *mcpart = mctrack->Particle();
  if(mcpart){
    datamc[14] = mcpart->GetUniqueID();
    datamc[19] = (mcpart->IsPrimary()) ? 1 : 0;
  }
  datamc[24] = mctrack->Label();
  if(fMCEvent) datamc[20] = fMCEvent->IsPhysicalPrimary(mctrack->Label()) ? 1 : 0;

  if(!(mArr<0)){
     if((mesonID>=kGammaPi0 && mesonID<=kGammaRho0) || mesonID==kGammaK2P || mesonID==kGammaK0s2P || mesonID==kGammaLamda2P) {  // conversion electron, be careful with the enum odering 
        Int_t glabel=TMath::Abs(mctrack->GetMother()); // gamma label
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's label
          if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
            mesonPt = mctrackmother->Pt(); //meson pt
            bgcategory = 1.;
            datamc[1] = bgcategory;
            datamc[2] = mesonPt;
            datamc[28] = mctrackmother->Eta();
            mctrackmother->XvYvZv(xr);
            datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
            datamc[12] = xr[2];
	   
            mcpart = mctrackmother->Particle();
            if(mcpart){
              datamc[15] = mcpart->GetUniqueID();
              datamc[21] = mcpart->IsPrimary() ? 1 : 0;
	      //if(TMath::Abs(AliHFEtools::GetRapidity(mcpart))>0.8) condition = kFALSE;
            }
            datamc[22] = fMCEvent->IsPhysicalPrimary(glabel) ? 1 : 0;
            datamc[25] = glabel;
            datamc[27] = fMCEvent->GetNumberOfPrimaries();
            if(glabel>fMCEvent->GetNumberOfPrimaries()) {
	      //condition = kFALSE;
              bgcategory = 1.; 
              datamc[1] = 2;
              //printf("I should be gamma meson = %d  mesonlabel= %d  NumberOfPrimaries= %d \n",mctrackmother->PdgCode(),glabel,fMCEvent->GetNumberOfPrimaries()); 
              glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
              datamc[26] = glabel;
              if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                datamc[3]=mctrackmother->PdgCode();
                mesonMotherPt=mctrackmother->Pt();
                datamc[4]=mesonMotherPt; 
                datamc[29]=mctrackmother->Eta();
                if(TMath::Abs(mctrackmother->PdgCode())==310){
                  bgcategory = 1.; // it was 3 in old code
                  datamc[1] = 3;
                  glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother's mother
                  if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                    datamc[5]=mctrackmother->PdgCode();
                    glabel=TMath::Abs(mctrackmother->GetMother()); // gamma's mother's mother
                    if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                       datamc[6]=mctrackmother->PdgCode();
                    }
                  }
                }
              }
            } 
          } 
        }
     }
     else{ // nonHFE except for the conversion electron
        Int_t glabel=TMath::Abs(mctrack->GetMother()); 
        if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
          mesonPt=mctrackmother->Pt(); //meson pt
          if(mesonID==kEta) bgcategory = -1.41; // to consider new branching ratio for the eta Dalitz decay
          else bgcategory = -1.;
          datamc[1] = bgcategory;
          datamc[2] = mesonPt;
          datamc[28] = mctrackmother->Eta();
          datamc[23] = mctrackmother->PdgCode();
          mctrackmother->XvYvZv(xr);
          datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
          datamc[12] = xr[2];

          mcpart = mctrackmother->Particle();
          if(mcpart){
            datamc[15] = mcpart->GetUniqueID();
            datamc[21] = mcpart->IsPrimary() ? 1 : 0;
	    //if(TMath::Abs(AliHFEtools::GetRapidity(mcpart))>0.8) condition = kFALSE;
          }
          datamc[22] = fMCEvent->IsPhysicalPrimary(glabel) ? 1 : 0;
          datamc[25] = glabel;
          datamc[27] = fMCEvent->GetNumberOfPrimaries();
          if(glabel>fMCEvent->GetNumberOfPrimaries()) {
	    //condition = kFALSE;
            if(mesonID==kEta) bgcategory = -1.41; // to consider new branching ratio for the eta Dalitz decay (1.41). It was -2.82 in old code
            else bgcategory = -1.; // it was -2 in old code
            datamc[1] = -2;
            glabel=TMath::Abs(mctrackmother->GetMother()); // mesons's mother
            datamc[26] = glabel;
            if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
              datamc[3]=mctrackmother->PdgCode();
              mesonMotherPt=mctrackmother->Pt();
              datamc[4]=mesonMotherPt;
              datamc[29]=mctrackmother->Eta();
              if(TMath::Abs(mctrackmother->PdgCode())==310){
               bgcategory = -1.; // it was -3 in old code
               datamc[1] = -3;
               glabel=TMath::Abs(mctrackmother->GetMother()); // mesons's mother's mother
               if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                 datamc[5]=mctrackmother->PdgCode();
                 glabel=TMath::Abs(mctrackmother->GetMother()); // mesons's mother's mother's mother 
                 if((mctrackmother = dynamic_cast<AliMCParticle *>(fMCEvent->GetTrack(glabel)))){
                   datamc[6]=mctrackmother->PdgCode();
                 }
               }
              }  
            }
          }
        }
     }

     Int_t centBin = 0;
     if(fIsPbPb){
       centBin = GetWeightCentralityBin(fPerCentrality);
       if(centBin < 0) return 0.;
     }

     weightElecBg=fElecBackgroundFactor[iBgLevel][centBin][mArr][kBgPtBins-1];                         
     if(mArr<=5 || mesonID==kKe3 || mesonID==kK0L){ 
       for(int ii=0; ii<kBgPtBins; ii++){              
	 if((mesonPt > fBinLimit[ii]) && (mesonPt < fBinLimit[ii+1])){
	   weightElecBg = fElecBackgroundFactor[iBgLevel][centBin][mArr][ii];
  	   break;
	 }
       }
     }
     else{
       for(int ii=0; ii<kBgPtBins; ii++){              
	 if((mesonMotherPt > fBinLimit[ii]) && (mesonMotherPt < fBinLimit[ii+1])){
	   weightElecBg = fElecBackgroundFactor[iBgLevel][centBin][mArr][ii];
  	   break;
	 }
       }
     }
  }

  //if(!condition) weightElecBg=0.;

  datamc[13] = weightElecBg;
  datamc[16] = Double_t(fContainerStep);

  datamc[7] = fHfeImpactR;
  datamc[8] = fHfeImpactnsigmaR;


  if(fIsDebugStreamerON && fTreeStream){
   if(!iBgLevel){
    (*fTreeStream)<<"nonhfeQA"<<
        "mesonID="<<datamc[0]<<
        "bgcategory="<<datamc[1]<<
        "mesonPt="<<datamc[2]<<
        "mesonMomPdg="<<datamc[3]<<
        "mesonMomPt="<<datamc[4]<<
        "mesonGMomPdg="<<datamc[5]<<
        "mesonGGMomPdg="<<datamc[6]<<
        "eIPAbs="<<datamc[7]<<
        "eIPSig="<<datamc[8]<<
        "eR="<<datamc[9]<<
        "eZ="<<datamc[10]<<
        "mesonR="<<datamc[11]<<
        "mesonZ="<<datamc[12]<<
        "weightElecBg="<<datamc[13]<< 
        "eUniqID="<<datamc[14]<<
        "mesonUniqID="<<datamc[15]<<
        "containerStep="<<datamc[16]<<
        "emcpt="<<datamc[17]<<
        "emceta="<<datamc[18]<<
        "eprim="<<datamc[19]<<
        "ephysiprim="<<datamc[20]<<
        "mprim="<<datamc[21]<< 
        "mphysiprim="<<datamc[22]<<
        "mesonPdg="<<datamc[23]<<
        "glabel="<<datamc[24]<<
        "motherglabel="<<datamc[25]<<
        "gmotherglabel="<<datamc[26]<<
        "NumberOfPrimaries="<<datamc[27]<<
        "mesonEta="<<datamc[28]<<
        "mesonMomEta="<<datamc[29]
        << "\n";
   }
  }

  Double_t returnval = bgcategory*weightElecBg;

  return returnval;
}

//__________________________________________
Double_t AliHFEmcQA::GetWeightFactor(const AliAODMCParticle * const mcpart, const Int_t iBgLevel){
  //
  // Get weighting factor for the realistic background estimation, for three possible background yield levels, indicated by the argument "iLevel": the best estimate (0), the lower uncertainty level (1), and the upper uncertainty level (2)
  //

  if (!mcpart) return 0;
  if (!fMCArray) return 0;

  Double_t weightElecBg = 0.;
  Double_t mesonPt = 0.;
  Double_t mesonMotherPt = 0.;
  Double_t bgcategory = 0.;
  Bool_t condition = kTRUE;
  Int_t mArr = -1;
  Int_t mesonID = GetElecSource(mcpart, kTRUE);
  if(mesonID==kGammaPi0 || mesonID==kPi0) mArr=0;                //pion
  else if(mesonID==kGammaEta || mesonID==kEta) mArr=1;           //eta
  else if(mesonID==kGammaOmega || mesonID==kOmega) mArr=2;       //omega
  else if(mesonID==kGammaPhi || mesonID==kPhi) mArr=3;           //phi
  else if(mesonID==kGammaEtaPrime || mesonID==kEtaPrime) mArr=4; //etaprime
  else if(mesonID==kGammaRho0 || mesonID==kRho0) mArr=5;         //rho
  else if(mesonID==kKe3 || mesonID==kGammaK2P|| mesonID==kK2P) mArr=6; //ke3 or K->pi->e
  else if(mesonID==kK0L || mesonID==kGammaK0s2P|| mesonID==kK0s2P) mArr=7; //K0L->e+X or k0s->pi->e
  else if(mesonID==kGammaLamda2P|| mesonID==kLamda2P) mArr=8;    //lambda->pi->e

  Double_t datamc[30]={-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999,-999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999, -999,-999};
  Double_t xr[3]={-999,-999,-999};
  datamc[0] = mesonID;
  datamc[17] = mcpart->Pt(); //electron pt
  datamc[18] = mcpart->Eta(); //electron eta

  mcpart->XvYvZv(xr);
  datamc[9] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
  datamc[10] = xr[2];

  datamc[19] = (mcpart->IsPrimary()) ? 1 : 0;
  datamc[20] = (mcpart->IsPhysicalPrimary()) ? 1 : 0;

  datamc[24] = mcpart->Label();

  if(!(mArr<0)){

     AliAODMCParticle *mctrackmother = NULL; // will change all the time

     if((mesonID>=kGammaPi0 && mesonID<=kGammaRho0) || mesonID==kGammaK2P || mesonID==kGammaK0s2P || mesonID==kGammaLamda2P) { // conversion electron
        Int_t iLabel = mcpart->GetMother(); //gamma label
        if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
        iLabel = mctrackmother->GetMother(); //gamma's mother's label
        if(!(mctrackmother= dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
        mesonPt = mctrackmother->Pt(); //meson pt
        datamc[28] = mctrackmother->Eta(); //meson eta
        bgcategory = 1.;
	if(TMath::Abs(AliHFEtools::GetRapidity(mctrackmother))>0.8) condition = kFALSE;
	if(!mctrackmother->IsPrimary()) condition = kFALSE;
        datamc[1] = bgcategory;
        datamc[2] = mesonPt;
        mctrackmother->XvYvZv(xr);
        datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
        datamc[12] = xr[2];
        
        datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
        datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;

        //bgcategory 2, 3 is not defined for AOD
 
        iLabel=mctrackmother->GetMother(); // gamma's mother's mother
        datamc[26] = iLabel;
        if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
          datamc[3]=mctrackmother->PdgCode();
          mesonMotherPt=mctrackmother->Pt();
          datamc[4]=mesonMotherPt;
          datamc[29] = mctrackmother->Eta(); //meson mother's eta
          if(TMath::Abs(mctrackmother->PdgCode())==310){
            datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
            datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;
            iLabel=mctrackmother->GetMother(); // gamma's mother's mother's mother
            if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
              datamc[5]=mctrackmother->PdgCode();
              iLabel=mctrackmother->GetMother(); // gamma's mother's mother's mother
              if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
                datamc[6]=mctrackmother->PdgCode();
              }
            }
          }
        }
     }
     else{ // nonHFE except for the conversion electron
        Int_t iLabel = mcpart->GetMother(); //meson label
        if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
        mesonPt = mctrackmother->Pt(); //meson pt
        datamc[28] = mctrackmother->Eta(); //meson eta
        if(mesonID==kEta) bgcategory = -1.41; // to consider new branching ratio for the eta Dalitz decay
        else bgcategory = -1.;
	if(TMath::Abs(AliHFEtools::GetRapidity(mctrackmother))>0.8) condition = kFALSE;
	if(!mctrackmother->IsPrimary()) condition = kFALSE;
        datamc[1] = bgcategory;
        datamc[2] = mesonPt;
        datamc[23] = mctrackmother->PdgCode();
        mctrackmother->XvYvZv(xr);
        datamc[11] = TMath::Sqrt(xr[0]*xr[0]+xr[1]*xr[1]);
        datamc[12] = xr[2];

        datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
        datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;

        //bgcategory 2, 3 is not defined for AOD

        iLabel=mctrackmother->GetMother(); // mesons' mother
        datamc[26] = iLabel;
        if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
          datamc[3]=mctrackmother->PdgCode();
          mesonMotherPt=mctrackmother->Pt();
          datamc[4]=mesonMotherPt;
          datamc[29] = mctrackmother->Eta(); //meson mother's eta
          if(TMath::Abs(mctrackmother->PdgCode())==310){
            datamc[21] = (mctrackmother->IsPrimary()) ? 1 : 0;
            datamc[22] = (mctrackmother->IsPhysicalPrimary()) ? 1 : 0;
            iLabel=mctrackmother->GetMother(); // mesons' mother's mother
            if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){
              datamc[5]=mctrackmother->PdgCode();
              iLabel=mctrackmother->GetMother(); // meson's mother's mother's mother
              if((mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))){ 
                datamc[6]=mctrackmother->PdgCode();
              }
            }
          }
        }
     }

     Int_t centBin = 0;
     if(fIsPbPb) {
       centBin = GetWeightCentralityBin(fPerCentrality);
       if(centBin < 0) return 0.;
     }

     weightElecBg=fElecBackgroundFactor[iBgLevel][centBin][mArr][kBgPtBins-1];
     if(mArr<=5 || mesonID==kKe3 || mesonID==kK0L){
       for(int ii=0; ii<kBgPtBins; ii++){
         if((mesonPt > fBinLimit[ii]) && (mesonPt < fBinLimit[ii+1])){
           weightElecBg = fElecBackgroundFactor[iBgLevel][centBin][mArr][ii];
           break;
         }
       }
     }
     else{
       for(int ii=0; ii<kBgPtBins; ii++){
         if((mesonMotherPt > fBinLimit[ii]) && (mesonMotherPt < fBinLimit[ii+1])){
           weightElecBg = fElecBackgroundFactor[iBgLevel][centBin][mArr][ii];
           break;
         }
       }
     }
  }

  if(!condition) weightElecBg=0.;

  datamc[13] = weightElecBg;
  datamc[16] = Double_t(fContainerStep);

  datamc[7] = fHfeImpactR;
  datamc[8] = fHfeImpactnsigmaR;


  if(fIsDebugStreamerON && fTreeStream){
   if(!iBgLevel){
    (*fTreeStream)<<"nonhfeQA"<<
        "mesonID="<<datamc[0]<<
        "bgcategory="<<datamc[1]<<
        "mesonPt="<<datamc[2]<<
        "mesonMomPdg="<<datamc[3]<<
        "mesonMomPt="<<datamc[4]<<
        "mesonGMomPdg="<<datamc[5]<<
        "mesonGGMomPdg="<<datamc[6]<<
        "eIPAbs="<<datamc[7]<<
        "eIPSig="<<datamc[8]<<
        "eR="<<datamc[9]<<
        "eZ="<<datamc[10]<<
        "mesonR="<<datamc[11]<<
        "mesonZ="<<datamc[12]<<
        "weightElecBg="<<datamc[13]<<
        "eUniqID="<<datamc[14]<<
        "mesonUniqID="<<datamc[15]<<
        "containerStep="<<datamc[16]<<
        "emcpt="<<datamc[17]<<
        "emceta="<<datamc[18]<<
        "eprim="<<datamc[19]<<
        "ephysiprim="<<datamc[20]<<
        "mprim="<<datamc[21]<<
        "mphysiprim="<<datamc[22]<<
        "mesonPdg="<<datamc[23]<<
        "glabel="<<datamc[24]<<
        "motherglabel="<<datamc[25]<<
        "gmotherglabel="<<datamc[26]<<
        "NumberOfPrimaries="<<datamc[27]<<
        "mesonEta="<<datamc[28]<<
        "mesonMomEta="<<datamc[29]
        << "\n";
   }
  }

  Double_t returnval = bgcategory*weightElecBg;

  return returnval;
}

//__________________________________________
Double_t AliHFEmcQA::GetWeightFactorForPrimaries(const AliAODMCParticle * const mcpart, const Int_t iBgLevel){
  //
  // Get weighting factor for the realistic background estimation, for three possible background yield levels, indicated by the argument "iLevel": the best estimate (0), the lower uncertainty level (1), and the upper uncertainty level (2)
  // weighting will only be non zero for electrons from primary pi0 and eta mesons (via Dalitz or gamma conversions)  
  //

  if (!mcpart) return 0;
  if (!fMCArray) return 0;

  Int_t mArr = -1;
  Double_t mesonPt = 0.;
  AliAODMCParticle *mctrackmother = NULL; // temp pointer
  Int_t mesonID = GetElecSource(mcpart, kTRUE); // get source of electron 
  if(mesonID==kGammaPi0 || mesonID==kPi0) mArr=0;                //pion
  else if(mesonID==kGammaEta || mesonID==kEta) mArr=1;           //eta
  else return 0;

  Int_t iLabel = mcpart->GetMother(); //mother label
  switch (mesonID) {
      case kGammaPi0:
      case kGammaEta:
          if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
          iLabel = mctrackmother->GetMother(); //reset label to mother of gamma
      case kPi0:
      case kEta:
          if(!(mctrackmother = dynamic_cast<AliAODMCParticle *>(fMCArray->At(TMath::Abs(iLabel))))) return 0;
          mesonPt = mctrackmother->Pt(); //pt of (grand)mother 
          //check mother is primary 
          if(!(mctrackmother->IsPrimary())) return 0;
          break;
      default:
          return 0; 
  }
  Double_t weightElecBg = fElecBackgroundFactor[iBgLevel][0][mArr][kBgPtBins-1]; //set weighting for pt > max pt
  for(int ii=0; ii<kBgPtBins; ii++){
      if((mesonPt >= fBinLimit[ii]) && (mesonPt < fBinLimit[ii+1])){
          weightElecBg = fElecBackgroundFactor[iBgLevel][0][mArr][ii];
          break;
      }
  }

  return weightElecBg;
}

//__________________________________________
//__________________________________________
Int_t AliHFEmcQA::GetMother(const AliVParticle * const mcpart) const {
  //
  // Wrapper to get the mother label
  //
  Int_t label = -1; 
  TClass *type = mcpart->IsA();
  if(type == AliMCParticle::Class()){
    // ESD analysis
    const AliMCParticle *emcpart = static_cast<const AliMCParticle *>(mcpart);
    label = emcpart->GetMother();
  } else if(type == AliAODMCParticle::Class()){
    // AOD analysis
    const AliAODMCParticle *amcpart = static_cast<const AliAODMCParticle *>(mcpart);
    label = amcpart->GetMother();
  }
  return label;
}
//__________________________________________
Int_t AliHFEmcQA::GetWeightCentralityBin(const Float_t percentile) const {
  //
  //translate the centrality percentile into the centrality bin of the reference weighting histograms for electron background
  //

  Float_t centralityLimits[12]= {0.,5.,10., 20., 30., 40., 50., 60.,70.,80., 90., 100.};
  Int_t bin = -1;
  for(Int_t ibin = 0; ibin < 11; ibin++){
    if(percentile >= centralityLimits[ibin] && percentile < centralityLimits[ibin+1]){
      bin = ibin;
      break;
    }
  } 
  return bin;
}
//__________________________________________
void AliHFEmcQA::AliHists::FillList(TList *l) const {
  //
  // Fill Histos into a list for output
  //
  if(fPdgCode) l->Add(fPdgCode);
  if(fPt) l->Add(fPt);
  if(fY) l->Add(fY);
  if(fEta) l->Add(fEta);
}

//__________________________________________
void AliHFEmcQA::AliHistsComm::FillList(TList *l) const { 
  //
  // Fill Histos into a list for output
  //
  if(fNq) l->Add(fNq);
  if(fProcessID) l->Add(fProcessID);
  if(fePtRatio) l->Add(fePtRatio);
  if(fPtCorr) l->Add(fPtCorr);
  if(fPtCorrDp) l->Add(fPtCorrDp);
  if(fPtCorrD0) l->Add(fPtCorrD0);
  if(fPtCorrDrest) l->Add(fPtCorrDrest);

  if(fPtCorrDinein) l->Add(fPtCorrDinein);
  if(fPtCorrDineout) l->Add(fPtCorrDineout);
  if(fPtCorrDoutein) l->Add(fPtCorrDoutein);
  if(fPtCorrDouteout) l->Add(fPtCorrDouteout);
  if(fPtCorrDpDinein) l->Add(fPtCorrDpDinein);
  if(fPtCorrDpDineout) l->Add(fPtCorrDpDineout);
  if(fPtCorrDpDoutein) l->Add(fPtCorrDpDoutein);
  if(fPtCorrDpDouteout) l->Add(fPtCorrDpDouteout);
  if(fPtCorrD0Dinein) l->Add(fPtCorrD0Dinein);
  if(fPtCorrD0Dineout) l->Add(fPtCorrD0Dineout);
  if(fPtCorrD0Doutein) l->Add(fPtCorrD0Doutein);
  if(fPtCorrD0Douteout) l->Add(fPtCorrD0Douteout);
  if(fPtCorrDrestDinein) l->Add(fPtCorrDrestDinein);
  if(fPtCorrDrestDineout) l->Add(fPtCorrDrestDineout);
  if(fPtCorrDrestDoutein) l->Add(fPtCorrDrestDoutein);
  if(fPtCorrDrestDouteout) l->Add(fPtCorrDrestDouteout);

  if(fEtaCorrD) l->Add(fEtaCorrD);
  if(fEtaCorrDp) l->Add(fEtaCorrDp);
  if(fEtaCorrD0) l->Add(fEtaCorrD0);
  if(fEtaCorrDrest) l->Add(fEtaCorrDrest);

  if(fEtaCorrGD) l->Add(fEtaCorrGD);
  if(fEtaCorrGDp) l->Add(fEtaCorrGDp);
  if(fEtaCorrGD0) l->Add(fEtaCorrGD0);
  if(fEtaCorrGDrest) l->Add(fEtaCorrGDrest);

  if(fEtaCorrB) l->Add(fEtaCorrB);
  if(fEtaCorrGB) l->Add(fEtaCorrGB);
  if(fPtCorrBinein) l->Add(fPtCorrBinein);
  if(fPtCorrBineout) l->Add(fPtCorrBineout);
  if(fPtCorrBoutein) l->Add(fPtCorrBoutein);
  if(fPtCorrBouteout) l->Add(fPtCorrBouteout);

  if(fDePtRatio) l->Add(fDePtRatio);
  if(feDistance) l->Add(feDistance);
  if(fDeDistance) l->Add(fDeDistance);
}

