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

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "THashList.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0Flow.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "TGeoManager.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include <AliAODCaloCluster.h>
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "TProfile.h"
#include "AliOADBContainer.h"
#include "AliEPFlattener.h"
#include "AliAnalysisUtils.h"

// Analysis task to fill histograms with PHOS ESD or AOD clusters and cells
// Authors : Dmitri Peressounko
// Date    : 28.05.2011
// Modified: 03.08.2012 Henrik Qvigstad
/* $Id: AliAnalysisTaskPi0Flow.cxx 64584 2013-10-17 15:37:28Z kharlov $ */

ClassImp(AliAnalysisTaskPi0Flow);

const Double_t AliAnalysisTaskPi0Flow::kLogWeight         = 4.5 ;
const Double_t AliAnalysisTaskPi0Flow::kAlphaCut          = 0.1 ;
const Bool_t   AliAnalysisTaskPi0Flow::doESDReCalibration = kTRUE;
const Double_t AliAnalysisTaskPi0Flow::kMinClusterEnergy  = 0.3;
const Double_t AliAnalysisTaskPi0Flow::kMinBCDistance     = 2.5;
const Int_t    AliAnalysisTaskPi0Flow::kMinNCells         = 3;
const Double_t AliAnalysisTaskPi0Flow::kMinM02            = 0.2;
const Int_t    AliAnalysisTaskPi0Flow::kNVtxZBins         = 1;
const Double_t AliAnalysisTaskPi0Flow::kCentCutoff        = 90.;

//________________________________________________________________________
Double_t rnlin(Double_t *x, Double_t * /*par*/)
{
  //a = par[0], b = par[1].
  //1+a*exp(-e/b)

  return 0.0241+1.0504*x[0]+0.000249*x[0]*x[0] ;
}

//________________________________________________________________________
AliAnalysisTaskPi0Flow::AliAnalysisTaskPi0Flow(const char *name, Period period)
: AliAnalysisTaskSE(name),
  fCentEdges(10),
  fCentNMixed(10),
  fNEMRPBins(9),
  fPeriod(period),
  fInternalTriggerSelection(kNoSelection),
  fMaxAbsVertexZ(10.),
  fManualV0EPCalc(false),
  fTOFCutEnabled(false),
  fTOFCut(100.e-9),
  fFillWideTOF(false),
  fTrigName(0x0),
  fOutputContainer(0x0),
  fNonLinCorr(0),
  fEvent(0x0),
  fEventESD(0x0),
  fEventAOD(0x0),
  fRunNumber(-999),
  fInternalRunNumber(0),
  fPHOSGeo(0),
  fMultV0(0x0),
  fV0Cpol(0.),fV0Apol(0.),
  fESDtrackCuts(0x0),
  fPHOSCalibData(0x0),
  fEPcalibFileName("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root"), 
  fTPCFlat(0x0),
  fV0AFlat(0x0),
  fV0CFlat(0x0),
  fVertexVector(),
  fVtxBin(0),
  fCentralityEstimator("V0M"),
  fCentrality(0.),
  fCentBin(0),
  fHaveTPCRP(0),
  fRP(0),
  fRPV0A(0),
  fRPV0C(0),
  fEMRPBin(0),
  fCaloPhotonsPHOS(0x0),
  fCaloPhotonsPHOSLists(0x0)
{
  const int nbins = 9;
  Double_t edges[nbins+1] = {0., 5., 10., 20., 30., 40., 50., 60., 70., 80.};
  TArrayD centEdges(nbins+1, edges);
  Int_t nMixed[nbins] = {4,4,6,10,20,30,50,100,100};
  TArrayI centNMixed(nbins, nMixed);
  SetCentralityBinning(centEdges, centNMixed);
  
  for(int mod=1; mod <= kNMod; ++mod)
    fModuleEnabled[mod-1] = kTRUE;

  for(Int_t i=0;i<kNCenBins;i++){
    for(Int_t j=0;j<2; j++)
      for(Int_t k=0; k<2; k++) {
        fMeanQ[i][j][k]=0.;
	fWidthQ[i][j][k]=0.;
      }
  }
  
  fVertex[0]=0; fVertex[1]=0; fVertex[2]=0; 

  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());



  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }


  // Initialize non-linrarity correction
  fNonLinCorr = new TF1("nonlib",rnlin,0.,40.,0);



}
//___________________________________________________________________________
AliAnalysisTaskPi0Flow::~AliAnalysisTaskPi0Flow()
{
  delete fNonLinCorr;
  delete fESDtrackCuts;
  delete fPHOSCalibData;
  delete fCaloPhotonsPHOSLists;
  if(fTPCFlat)delete fTPCFlat;  fTPCFlat=0x0;
  if(fV0AFlat)delete fV0AFlat;  fV0AFlat=0x0;
  if(fV0CFlat)delete fV0CFlat;  fV0CFlat=0x0;
  
}
//________________________________________________________________________
void AliAnalysisTaskPi0Flow::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  const Int_t nRuns=200 ;

  // histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  //========QA histograms=======

  //Event selection
  fOutputContainer->Add(new TH2F("hSelEvents","Event selection", kTotalSelected+1, 0., double(kTotalSelected+1), nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH1F("hTotSelEvents","Event selection", kTotalSelected+1, 0., double(kTotalSelected+1))) ;

  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;

  //Centrality
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu","ClusterMult vs E",200,0.,10.,100,0.,100.)) ;


  //Reaction plane
  fOutputContainer->Add(new TH3F("hPHOSphi","cos" ,10,0.,100.,20,0.,10.,100,-TMath::Pi(),TMath::Pi()));

  fOutputContainer->Add(new TH2F("cos2AC","RP correlation between TPC subs", 100,-1.,1.,20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("cos2V0AC","RP correlation between VO A and C sides", 100,-1.,1.,20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("cos2V0ATPC","RP correlation between TPC and V0A", 100,-1.,1.,20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("cos2V0CTPC","RP correlation between TPC and V0C", 100,-1.,1.,20,0.,100.)) ;

  fOutputContainer->Add(new TH2F("phiRP","RP distribution with TPC", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPflat","RP distribution with TPC flat", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0A","RP distribution with V0A", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0C","RP distribution with V0C", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH3F("phiRPV0AC","RP distribution with V0A and V0C", 100,0.,TMath::Pi(),100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0Aflat","RP distribution with V0 flat", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0Cflat","RP distribution with V0 flat", 100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH3F("phiRPV0ATPC","RP distribution with V0A + TPC", 100,0.,TMath::Pi(),100,0.,TMath::Pi(),20,0.,100.)) ;
  fOutputContainer->Add(new TH3F("phiRPV0CTPC","RP distribution with V0C + TPC", 100,0.,TMath::Pi(),100,0.,TMath::Pi(),20,0.,100.)) ;


  //PHOS QA
  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM4","PHOS cell multiplicity per event, M4",2000,0,2000));

  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM1","Cell energy in module 1",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM2","Cell energy in module 2",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM3","Cell energy in module 3",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM4","Cell energy in module 4",3000,0.,30.));

  fOutputContainer->Add(new TH2F("hCellNXZM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM4","Cell (X,Z), M4" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCellEXZM1","Cell E(X,Z), M1",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM2","Cell E(X,Z), M2",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM3","Cell E(X,Z), M3",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM4","Cell E(X,Z), M4",64,0.5,64.5, 56,0.5,56.5));

  //Bad Map
  fOutputContainer->Add(new TH2F("hCluLowM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluLowM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluLowM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluLowM4","Cell (X,Z), M4" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluHighM1","Cell (X,Z), M1 (E_{clu}>1.5 GeV)" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluHighM2","Cell (X,Z), M2 (E_{clu}>1.5 GeV)" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluHighM3","Cell (X,Z), M3 (E_{clu}>1.5 GeV)" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluHighM4","Cell (X,Z), M4 (E_{clu}>1.5 GeV)" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluVetoM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluVetoM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluVetoM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluVetoM4","Cell (X,Z), M4" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluDispM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluDispM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluDispM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluDispM4","Cell (X,Z), M4" ,64,0.5,64.5, 56,0.5,56.5));


  //Single photon and pi0 spectrum
  const Int_t nPtPhot = 400 ;
  const Double_t ptPhotMax = 40 ;
  const Int_t nM       = 500;
  const Double_t mMin  = 0.0;
  const Double_t mMax  = 1.0;

  //PHOS calibration QA
  fOutputContainer->Add(new TH2F("hPi0M11","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M12","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M13","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M14","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M22","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M23","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M24","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M33","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M34","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M44","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

  // Histograms for different centralities
  const int kNPID = 17;
  const char* pidNames[kNPID] = {"All", "Allcore", "Allwou", "Disp", "Disp2", "Dispcore",  "Disp2core", "Dispwou", "CPV", "CPVcore", "CPV2", "CPV2core", "Both", "Bothcore", "Both2", "Both2core", "WideTOF"};
  char key[55];
  TString name, title;
  for(Int_t cent=0; cent < fCentEdges.GetSize()-1; cent++){
    for(Int_t ipid=0; ipid < kNPID; ipid++){
      if( !fFillWideTOF && TString(pidNames[ipid]).EqualTo("WideTOF") ) continue;

      name = Form("hPhot%s_cen%i", pidNames[ipid], cent );
      title = Form("%s clusters", pidNames[ipid]);
      fOutputContainer->Add(new TH1F(name.Data(), title.Data(), nPtPhot,0.,ptPhotMax));

      name = Form("hPi0%s_cen%i", pidNames[ipid], cent );
      title = Form("%s clusters", pidNames[ipid]);
      fOutputContainer->Add(new TH2F(name.Data(), title.Data(), nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

      name = Form("hSingle%s_cen%i", pidNames[ipid], cent );
      title = Form("%s clusters", pidNames[ipid]);
      fOutputContainer->Add(new TH2F(name.Data(), title.Data(), nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

      name = Form("hMiPi0%s_cen%i", pidNames[ipid], cent );
      title = Form("%s clusters", pidNames[ipid]);
      fOutputContainer->Add(new TH2F(name.Data(), title.Data(), nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

      name = Form("hMiSingle%s_cen%i", pidNames[ipid], cent );
      title = Form("%s clusters", pidNames[ipid]);
      fOutputContainer->Add(new TH2F(name.Data(), title.Data(), nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

      // PhotPhi histograms
      const Int_t nPt      = 20;
      const Double_t xPt[21]={0.6,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.} ;
      const Int_t nPhi=10 ;
      Double_t xPhi[nPhi+1] ;
      for(Int_t i=0;i<=nPhi;i++)
	xPhi[i]=i*TMath::Pi() /nPhi ;
      const Int_t nMm=200 ;
      Double_t xM[nMm+1] ;
      for(Int_t i=0;i<=nMm;i++)
	xM[i]=i*0.5 /nMm;
      const Int_t kNPhiTitles = 3;
      const char* phiTitles[kNPhiTitles] = {"TPC", "V0A", "V0C"};
      for(Int_t iRP=0; iRP<3; iRP++){
	name = Form("hPhotPhi%s%s_cen%i", phiTitles[iRP], pidNames[ipid], cent );
	title = Form("(M,p_{T},d#phi)_{#gamma#gamma}");
	fOutputContainer->Add(new TH2F(name.Data(), title.Data(), nPt,xPt,nPhi,xPhi));
    
	name = Form("hMassPt%s%s_cen%i", phiTitles[iRP], pidNames[ipid], cent );
	title = Form("(M,p_{T},d#phi)_{#gamma#gamma}");
	fOutputContainer->Add(new TH3F(name.Data(), title.Data(), nMm,xM,nPt,xPt,nPhi,xPhi));
    
	name = Form("hMiMassPt%s%s_cen%i", phiTitles[iRP], pidNames[ipid], cent );
	title = Form("(M,p_{T},d#phi)_{#gamma#gamma}");
	fOutputContainer->Add(new TH3F(name.Data(), title.Data(), nMm,xM,nPt,xPt,nPhi,xPhi));
      }
    }
    
    // a07 histograms
    snprintf(key,55,"hPi0All_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0Disp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0CPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0CPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0Both_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

    snprintf(key,55,"hMiPi0All_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Disp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Both_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  }
  
  // Setup photon lists
  Int_t kapacity = kNVtxZBins * GetNumberOfCentralityBins() * fNEMRPBins;
  fCaloPhotonsPHOSLists = new TObjArray(kapacity);
  fCaloPhotonsPHOSLists->SetOwner();
  
  PostData(1, fOutputContainer);
}

void AliAnalysisTaskPi0Flow::ProcessMC()
{
  //empty, MC extensions should override
}


//________________________________________________________________________
void AliAnalysisTaskPi0Flow::UserExec(Option_t *)
{
  // Main loop, called for each event
  // Analyze ESD/AOD


  // Step 0: Event Objects
  fEvent = GetEvent();
  fEventESD = dynamic_cast<AliESDEvent*> (fEvent);
  fEventAOD = dynamic_cast<AliAODEvent*> (fEvent);
  LogProgress(0);


  // Step 1: Run Number, Misalignment Matrix, and Calibration
  // fRunNumber, fInternalRunNumber, fMultV0, fV0Cpol, fV0Apol, fMeanQ, fWidthQ
  if( fRunNumber != fEvent->GetRunNumber()) { // Check run number
    // this should run only at first call of UserExec(),
    // or if task runs over multiple runs, which should not occur in normal use.

    // if run number has changed, set run variables
    fRunNumber = fEvent->GetRunNumber();
    fInternalRunNumber = ConvertToInternalRunNumber(fRunNumber);
    // then set misalignment and V0 calibration
    SetGeometry();
    SetMisalignment();
    SetV0Calibration();
    SetESDTrackCuts();
    SetPHOSCalibData();
    SetFlatteningData();
  }
  LogProgress(1);
  LogSelection(kTotal, fInternalRunNumber);


  // Step 2: Internal Trigger Selection
  if( RejectTriggerMaskSelection() ) {
    PostData(1, fOutputContainer);
    return; // Reject!
  }

  if( RejectFiredTriggerClassSelection() ) {
    PostData(1, fOutputContainer);
    return; // Reject!
  }
  
  // Step 3: Vertex
  // fVertex, fVertexVector, fVtxBin
  SetVertex();
  if( RejectEventVertex() ) {
    PostData(1, fOutputContainer);
    return; // Reject!
  }
  LogProgress(2);

// Step 3:
//   if(event->IsPileupFromSPD()){
//     PostData(1, fOutputContainer);
//     return; // Reject!
//   }
  LogProgress(3);


  // Step 4: Centrality
  // fCentrality, fCentBin
  if (fRunNumber < 224994) SetCentrality();
    else SetCentralityRun2();
    
  if( RejectCentrality() ){
    PostData(1, fOutputContainer);
    return; // Reject! 
  }
  LogProgress(4);


  // Step 5: Reaction Plane
  // fHaveTPCRP, fRP, fRPV0A, fRPV0C, fRPBin
  EvalReactionPlane(); //TODO: uncomment this, or at least deal with it
  EvalV0ReactionPlane(); //TODO: uncomment this, or at least deal with it
  fEMRPBin = GetRPBin(); //TODO: uncomment this, or at least deal with it
  LogProgress(5);

  // Step 6: QA PHOS cells
  FillPHOSCellQAHists();
  LogProgress(6);

  // Step 7: Event Photons (PHOS Clusters) selection
  this->SelectPhotonClusters();
  this->FillSelectedClusterHistograms();
  LogProgress(7);

  // Step 8: MC
  this->ProcessMC() ;
  LogProgress(8);

  if( ! fCaloPhotonsPHOS->GetEntriesFast() )
    return;
  else
    LogSelection(kHasPHOSClusters, fInternalRunNumber);

  LogSelection(kTotalSelected, fInternalRunNumber);

  // Step 9: Consider pi0 (photon/cluster) pairs.
  this->ConsiderPi0s();
  LogProgress(9);

  // Step 10; Mixing
  this->ConsiderPi0sMix();
  LogProgress(10);
  
  // Step 12: Update lists
  UpdateLists();
  LogProgress(11);

  
  // Post output data.
  PostData(1, fOutputContainer);
}

//________________________________________________________________________
// void AliAnalysisTaskPi0Flow::Terminate(Option_t *)
// {
//   // Draw result to the screen
//   // Called once at the end of the query
//   // new TCanvas;
//   // TH1 * hTotSelEvents = dynamic_cast<TH1*>(fOutputContainer->FindObject("hTotSelEvents"));
//   // hTotSelEvents->Draw();
// }
//________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetCentralityBinning(const TArrayD& edges, const TArrayI& nMixed)
{
  // Define centrality bins by their edges
  for(int i=0; i<edges.GetSize()-1; ++i)
    if(edges.At(i) > edges.At(i+1)) AliFatal("edges are not sorted");
  if( edges.GetSize() != nMixed.GetSize()+1) AliFatal("edges and nMixed don't have appropriate relative sizes");
  
  fCentEdges = edges;
  fCentNMixed = nMixed;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetEnablePHOSModule(int module, Bool_t enable)
{
  if( module < 1 || kNMod < module )
    AliFatal(Form("PHOS Module must be between 1 and %i", kNMod));
  else
    fModuleEnabled[module-1] = enable;
}


//________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetPHOSBadMap(Int_t mod, TH2I* badMapHist)
  {
    if(fPHOSBadMap[mod])
      delete fPHOSBadMap[mod];

    fPHOSBadMap[mod]=new TH2I(*badMapHist);
    if(fDebug)
      AliInfo(Form("Setting Bad Map Histogram  %s",fPHOSBadMap[mod]->GetName()));
  }

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d",mod)) ;
      return kTRUE ;
    }
    if(!fPHOSBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d, !fPHOSBadMap[mod]",mod)) ;
      return kTRUE ;
    }
    if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
      return kFALSE ;
    else
      return kTRUE ;
  }
  else{
    AliError(Form("Can not find bad channels for detector %s ",det)) ;
  }

  //Remove 6 noisy channels in run 139036, LHC10h
  if( 139036 == fRunNumber
    && mod==1
    && (ix==9||ix==10||ix==11)
    && (iz==45 || iz==46))
    return kFALSE;

  return kTRUE ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::FillPHOSCellQAHists()
{
  // Fill cell occupancy per module

  AliVCaloCells * cells = fEvent->GetPHOSCells();

  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,fEvent->GetNumberOfTracks()) ;


  Int_t nCellModule[4] = {0,0,0,0};
  for (Int_t iCell=0; iCell<cells->GetNumberOfCells(); iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    Int_t relId[4] = {0,0,0,0};
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    Int_t mod1  = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    Float_t energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if(mod1==1) {
      nCellModule[0]++;
      FillHistogram("hCellEnergyM1",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM1",cellX,cellZ,1.);
      FillHistogram("hCellEXZM1",cellX,cellZ,energy);
    }
    else if (mod1==2) {
      nCellModule[1]++;
      FillHistogram("hCellEnergyM2",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM2",cellX,cellZ,1.);
      FillHistogram("hCellEXZM2",cellX,cellZ,energy);
    }
    else if (mod1==3) {
      nCellModule[2]++;
      FillHistogram("hCellEnergyM3",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM3",cellX,cellZ,1.);
      FillHistogram("hCellEXZM3",cellX,cellZ,energy);
    }
    else if (mod1==4) {
      nCellModule[3]++;
      FillHistogram("hCellEnergyM4",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM4",cellX,cellZ,1.);
      FillHistogram("hCellEXZM4",cellX,cellZ,energy);
    }
  }
  FillHistogram("hCellMultEvent",nCellModule[0]+nCellModule[1]+nCellModule[2]+nCellModule[3]);
  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  FillHistogram("hCellMultEventM4",nCellModule[3]);

}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SelectPhotonClusters()
{
  // clear (or create) array for holding events photons/clusters
  if(fCaloPhotonsPHOS)
    fCaloPhotonsPHOS->Clear();
  else {
    fCaloPhotonsPHOS = new TObjArray(200);
    fCaloPhotonsPHOS->SetOwner();
  }

  
  AliVCaloCells* cells = dynamic_cast<AliVCaloCells*> (fEvent->GetPHOSCells());
  for (Int_t i=0; i<fEvent->GetNumberOfCaloClusters(); i++) {
    AliVCluster *clu = fEvent->GetCaloCluster(i);

    if ( clu->GetType() !=AliVCluster::kPHOSNeutral ) continue; // reject CPV clusters
    if ( !clu->IsPHOS() || clu->E()< kMinClusterEnergy) continue; // reject cluster



    // check if cell/channel is good.
    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    if ( ! fModuleEnabled[mod-1] )
      continue;
    if ( !IsGoodChannel("PHOS",mod,cellX,cellZ) )
      continue ; // reject if not.

    Double_t distBC=clu->GetDistanceToBadChannel();
    if(distBC<kMinBCDistance)
      continue ;
      
    FillHistogram("hCluEvsClu", clu->E(), clu->GetNCells()) ;

    if(clu->GetNCells() < kMinNCells) continue ;
    if(clu->GetM02() < kMinM02)   continue ;

    TLorentzVector lorentzMomentum;
    Double_t ecore;
    ecore = CoreEnergy(clu,cells);

    Double_t origo[3] = {0,0,0}; // don't rely on event vertex, assume (0,0,0)

    AliESDCaloCluster* aodCluster = (AliESDCaloCluster*) (clu);
    aodCluster->GetMomentum(lorentzMomentum ,origo);

    FillHistogram(Form("hCluLowM%d",mod),cellX,cellZ,1.);
    if(lorentzMomentum.E()>1.5){
      FillHistogram(Form("hCluHighM%d",mod),cellX,cellZ,1.);
    }

    fCaloPhotonsPHOS->Add(new  AliCaloPhoton(lorentzMomentum.X(),lorentzMomentum.Py(),lorentzMomentum.Z(),lorentzMomentum.E()) );
    AliCaloPhoton * ph = (AliCaloPhoton*) fCaloPhotonsPHOS->At( fCaloPhotonsPHOS->GetLast() );
    ph->SetCluster(clu);

    ph->SetModule(mod) ;
    lorentzMomentum*= ecore/lorentzMomentum.E() ;
    ph->SetMomV2(&lorentzMomentum) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    //Evaluate CoreDispersion
    Double_t m02=0.,m20=0. ;
    EvalCoreLambdas(clu, cells, m02, m20) ;
    ph->SetDisp2Bit(TestCoreLambda(clu->E(),m20,m02)) ; //Correct order m20,m02
//    ph->SetDisp2Bit(TestCoreLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    if(ph->IsDispOK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
    }

    // Track Matching
    Double_t dx=clu->GetTrackDx() ;
    Double_t dz=clu->GetTrackDz() ;
    Bool_t cpvBit=kTRUE ; //No track matched by default. True means: not from charged, according to veto.
    Bool_t cpvBit2=kTRUE ; //More Strict criterion
    if( fEventESD ) {
      
      TArrayI * itracks = static_cast<AliESDCaloCluster*> (clu)->GetTracksMatched() ;
      if(itracks->GetSize()>0){
	Int_t iTr = itracks->At(0);
	if(iTr>=0 && iTr<fEvent->GetNumberOfTracks()){
	  AliVParticle* track = fEvent->GetTrack(iTr);
	  Double_t pt = track->Pt() ;
	  Short_t charge = track->Charge() ;
	  Double_t r=TestCPV(dx, dz, pt, charge) ;
	  cpvBit=(r>2.) ;
	  cpvBit2=(r>4.) ;
	}
      }
    }
    else if ( fEventAOD ) {
      int nTracksMatched = clu->GetNTracksMatched();
      if(nTracksMatched > 0) {
	AliVTrack* track = dynamic_cast<AliVTrack*> (clu->GetTrackMatched(0));
	if ( track ) {
	  Double_t pt = track->Pt();
	  Short_t charge = track->Charge();
	  Double_t r = TestCPV(dx, dz, pt, charge) ;
	  cpvBit=(r>2.) ;
	  cpvBit2=(r>4.) ;
	}
      }
    }
    ph->SetCPVBit(cpvBit) ;
    ph->SetCPV2Bit(cpvBit2) ;
    if(cpvBit){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
    }
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;
    //    ph->SetLambdas(clu->GetM20(),clu->GetM02()) ;
    ph->SetUnfolded(clu->GetNExMax()<2); // Remember, if it is unfolde

    // Time of Flight (TOF)
    Double_t tof = clu->GetTOF();
    ph->SetTOFBit( TMath::Abs(tof) < fTOFCut );
  }
  FillHistogram("hCenPHOS",fCentrality, fCaloPhotonsPHOS->GetEntriesFast()) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::FillSelectedClusterHistograms()
{
  for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;

    Double_t dphiA=ph1->Phi()-fRPV0A ;
    while(dphiA<0)dphiA+=TMath::Pi() ;
    while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

    Double_t dphiC=ph1->Phi()-fRPV0C ;
    while(dphiC<0)dphiC+=TMath::Pi() ;
    while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

    Double_t dphiT=ph1->Phi()-fRP ;
    while(dphiT<0)dphiT+=TMath::Pi() ;
    while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;
    
    Double_t pt = ph1->Pt() ;
    Double_t ptcore = ph1->GetMomV2()->Pt() ;

    if( fFillWideTOF ) {
      FillHistogram(Form("hPhotWideTOF_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotPhiV0AWideTOF_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CWideTOF_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
	FillHistogram(Form("hPhotPhiTPCWideTOF_cen%d",fCentBin),pt,dphiT) ;
    }
    if(fTOFCutEnabled && !ph1->IsTOFOK() )
      continue;

    FillHistogram(Form("hPhotPhiV0AAll_cen%d",fCentBin),pt,dphiA) ;
    FillHistogram(Form("hPhotPhiV0CAll_cen%d",fCentBin),pt,dphiC) ;
    if(fHaveTPCRP)
      FillHistogram(Form("hPhotPhiTPCAll_cen%d",fCentBin),pt,dphiT) ;
    FillHistogram(Form("hPhotPhiV0AAllcore_cen%d",fCentBin),ptcore,dphiA) ;
    FillHistogram(Form("hPhotPhiV0CAllcore_cen%d",fCentBin),ptcore,dphiC) ;
    if(fHaveTPCRP)
      FillHistogram(Form("hPhotPhiTPCAllcore_cen%d",fCentBin),ptcore,dphiT) ;

    FillHistogram(Form("hPhotAll_cen%d",fCentBin),pt) ;
    FillHistogram(Form("hPhotAllcore_cen%d",fCentBin),ptcore) ;
    if(ph1->IsntUnfolded()){
      FillHistogram(Form("hPhotAllwou_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotPhiV0AAllwou_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CAllwou_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCAllwou_cen%d",fCentBin),pt,dphiT) ;
    }
    if(ph1->IsCPVOK()){
      FillHistogram(Form("hPhotPhiV0ACPV_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CCPV_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV_cen%d",fCentBin),pt,dphiT) ;

      FillHistogram(Form("hPhotPhiV0ACPVcore_cen%d",fCentBin),ptcore,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CCPVcore_cen%d",fCentBin),ptcore,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPVcore_cen%d",fCentBin),ptcore,dphiT) ;

      FillHistogram(Form("hPhotCPV_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCentBin),ptcore) ;
    }
    if(ph1->IsCPV2OK()){
      FillHistogram(Form("hPhotPhiV0ACPV2_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CCPV2_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV2_cen%d",fCentBin),pt,dphiT) ;

      FillHistogram(Form("hPhotPhiV0ACPV2core_cen%d",fCentBin),ptcore,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CCPV2core_cen%d",fCentBin),ptcore,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCCPV2core_cen%d",fCentBin),ptcore,dphiT) ;
      FillHistogram(Form("hPhotCPV2_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotCPV2core_cen%d",fCentBin),ptcore) ;
    }
    if(ph1->IsDispOK()){
      FillHistogram(Form("hPhotPhiV0ADisp_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CDisp_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp_cen%d",fCentBin),pt,dphiT) ;

      FillHistogram(Form("hPhotPhiV0ADispcore_cen%d",fCentBin),ptcore,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CDispcore_cen%d",fCentBin),ptcore,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDispcore_cen%d",fCentBin),ptcore,dphiT) ;

      if(ph1->IsntUnfolded()){
        FillHistogram(Form("hPhotPhiV0ADispwou_cen%d",fCentBin),pt,dphiA) ;
        FillHistogram(Form("hPhotPhiV0CDispwou_cen%d",fCentBin),pt,dphiC) ;
        if(fHaveTPCRP)
          FillHistogram(Form("hPhotPhiTPCDispwou_cen%d",fCentBin),pt,dphiT) ;

      }
      FillHistogram(Form("hPhotDisp_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotDispcore_cen%d",fCentBin),ptcore) ;
      if(ph1->IsntUnfolded()){
        FillHistogram(Form("hPhotDispwou_cen%d",fCentBin),pt) ;
      }
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotPhiV0ABoth_cen%d",fCentBin),pt,dphiA) ;
	FillHistogram(Form("hPhotPhiV0CBoth_cen%d",fCentBin),pt,dphiC) ;
        if(fHaveTPCRP)
  	  FillHistogram(Form("hPhotPhiTPCBoth_cen%d",fCentBin),pt,dphiT) ;

	FillHistogram(Form("hPhotPhiV0ABothcore_cen%d",fCentBin),ptcore,dphiA) ;
	FillHistogram(Form("hPhotPhiV0CBothcore_cen%d",fCentBin),ptcore,dphiC) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hPhotPhiTPCBothcore_cen%d",fCentBin),ptcore,dphiT) ;

	FillHistogram(Form("hPhotBoth_cen%d",fCentBin),pt) ;
	FillHistogram(Form("hPhotBothcore_cen%d",fCentBin),ptcore) ;
      }
    }
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hPhotPhiV0ADisp2_cen%d",fCentBin),pt,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CDisp2_cen%d",fCentBin),pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp2_cen%d",fCentBin),pt,dphiT) ;
      FillHistogram(Form("hPhotPhiV0ADisp2core_cen%d",fCentBin),ptcore,dphiA) ;
      FillHistogram(Form("hPhotPhiV0CDisp2core_cen%d",fCentBin),ptcore,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hPhotPhiTPCDisp2core_cen%d",fCentBin),ptcore,dphiT) ;

      FillHistogram(Form("hPhotDisp2_cen%d",fCentBin),pt) ;
      FillHistogram(Form("hPhotDisp2core_cen%d",fCentBin),ptcore) ;
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotPhiV0ABoth2_cen%d",fCentBin),pt,dphiA) ;
	FillHistogram(Form("hPhotPhiV0CBoth2_cen%d",fCentBin),pt,dphiC) ;
        if(fHaveTPCRP)
  	  FillHistogram(Form("hPhotPhiTPCBoth2_cen%d",fCentBin),pt,dphiT) ;

	FillHistogram(Form("hPhotPhiV0ABoth2core_cen%d",fCentBin),ptcore,dphiA) ;
	FillHistogram(Form("hPhotPhiV0CBoth2core_cen%d",fCentBin),ptcore,dphiC) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hPhotPhiTPCBoth2core_cen%d",fCentBin),ptcore,dphiT) ;

	FillHistogram(Form("hPhotBoth2_cen%d",fCentBin),pt) ;
	FillHistogram(Form("hPhotBoth2core_cen%d",fCentBin),ptcore) ;
      }
    }
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::ConsiderPi0s()
{
  char key[55];
  for (Int_t i1=0; i1 < fCaloPhotonsPHOS->GetEntriesFast()-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    for (Int_t i2=i1+1; i2<fCaloPhotonsPHOS->GetEntriesFast(); i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      FillHistogram("hPHOSphi",fCentrality,p12.Pt(),p12.Phi());
      Double_t dphiA=p12.Phi()-fRPV0A ;
      while(dphiA<0)dphiA+=TMath::Pi() ;
      while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

      Double_t dphiC=p12.Phi()-fRPV0C ;
      while(dphiC<0)dphiC+=TMath::Pi() ;
      while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

      Double_t dphiT=p12.Phi()-fRP ;
      while(dphiT<0)dphiT+=TMath::Pi() ;
      while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;

      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
      Double_t m=p12.M() ;
      Double_t mcore=pv12.M() ;
      Double_t pt=p12.Pt() ;
      Double_t ptcore=pv12.Pt() ;
      Double_t pt1=ph1->Pt() ;
      Double_t pt2=ph2->Pt() ;
      Double_t ptcore1=ph1->GetMomV2()->Pt() ;
      Double_t ptcore2=ph2->GetMomV2()->Pt() ;

      if( fFillWideTOF ) {
	FillHistogram(Form("hPi0WideTOF_cen%d",fCentBin),m,pt) ;
	FillHistogram(Form("hSingleWideTOF_cen%d",fCentBin),m,pt1) ;
	FillHistogram(Form("hSingleWideTOF_cen%d",fCentBin),m,pt2) ;
	if(fHaveTPCRP)
	  FillHistogram(Form("hMassPtTPCWideTOF_cen%d",fCentBin),m,pt,dphiT) ;
      }

      if( fTOFCutEnabled && !(ph1->IsTOFOK() && ph2->IsTOFOK()) )
	continue;

      FillHistogram(Form("hMassPtV0AAll_cen%d",fCentBin),m,pt,dphiA) ;
      FillHistogram(Form("hMassPtV0CAll_cen%d",fCentBin),m,pt,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hMassPtTPCAll_cen%d",fCentBin),m,pt,dphiT) ;

      FillHistogram(Form("hMassPtV0AAllcore_cen%d",fCentBin),mcore,ptcore,dphiA) ;
      FillHistogram(Form("hMassPtV0CAllcore_cen%d",fCentBin),mcore,ptcore,dphiC) ;
      if(fHaveTPCRP)
        FillHistogram(Form("hMassPtTPCAllcore_cen%d",fCentBin),mcore,ptcore,dphiT) ;


      FillHistogram(Form("hPi0All_cen%d",fCentBin),m,pt) ;
      FillHistogram(Form("hPi0Allcore_cen%d",fCentBin),mcore,ptcore) ;
      if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
        FillHistogram(Form("hPi0Allwou_cen%d",fCentBin),m,pt) ;
        FillHistogram(Form("hMassPtV0AAllwou_cen%d",fCentBin),m,pt,dphiA) ;
        FillHistogram(Form("hMassPtV0CAllwou_cen%d",fCentBin),m,pt,dphiC) ;
        if(fHaveTPCRP)
          FillHistogram(Form("hMassPtTPCAllwou_cen%d",fCentBin),m,pt,dphiT) ;
      }

      FillHistogram(Form("hSingleAll_cen%d",fCentBin),m,pt1) ;
      FillHistogram(Form("hSingleAll_cen%d",fCentBin),m,pt2) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCentBin),mcore,ptcore1) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCentBin),mcore,ptcore2) ;
      if(ph1->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCentBin),m,pt1) ;
      if(ph2->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCentBin),m,pt2) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph2->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCentBin),mcore,ptcore2) ;
      }
      if(ph1->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hSingleCPV2core_cen%d",fCentBin),mcore,ptcore2) ;
      }
      if(ph2->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hSingleCPV2core_cen%d",fCentBin),mcore,ptcore2) ;
      }
      if(ph1->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCentBin),m,pt1) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCentBin),m,pt1) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph2->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCentBin),m,pt2) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCentBin),m,pt2) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCentBin),mcore,ptcore2) ;
      }
      if(ph1->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph2->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph1->IsDispOK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph2->IsDispOK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCentBin),mcore,ptcore2) ;
      }
      if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCentBin),mcore,ptcore1) ;
      }
      if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCentBin),mcore,ptcore2) ;
      }


      if(a<kAlphaCut){
        FillHistogram(Form("hPi0All_a07_cen%d",fCentBin),m,pt) ;
      }

      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCentBin) ;
	FillHistogram(Form("hMassPtV0ACPV_cen%d",fCentBin),m,pt,dphiA) ;
	FillHistogram(Form("hMassPtV0CCPV_cen%d",fCentBin),m,pt,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV_cen%d",fCentBin),m,pt,dphiT) ;

	FillHistogram(Form("hMassPtV0ACPVcore_cen%d",fCentBin),mcore,ptcore,dphiA) ;
	FillHistogram(Form("hMassPtV0CCPVcore_cen%d",fCentBin),mcore,ptcore,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPVcore_cen%d",fCentBin),mcore,ptcore,dphiT) ;

	FillHistogram(Form("hPi0CPV_cen%d",fCentBin),m,pt) ;
	FillHistogram(Form("hPi0CPVcore_cen%d",fCentBin),mcore, ptcore) ;

        if(a<kAlphaCut){
          FillHistogram(Form("hPi0CPV_a07_cen%d",fCentBin),m,pt) ;
        }
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hMassPtV0ACPV2_cen%d",fCentBin),m,pt,dphiA) ;
	FillHistogram(Form("hMassPtV0CCPV2_cen%d",fCentBin),m,pt,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV2_cen%d",fCentBin),m,pt,dphiT) ;
	FillHistogram(Form("hMassPtV0ACPV2core_cen%d",fCentBin),mcore,ptcore,dphiA) ;
	FillHistogram(Form("hMassPtV0CCPV2core_cen%d",fCentBin),mcore,ptcore,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCCPV2core_cen%d",fCentBin),mcore,ptcore,dphiT) ;
	
	FillHistogram(Form("hPi0CPV2_cen%d",fCentBin),m,pt) ;
	FillHistogram(Form("hPi0CPV2core_cen%d",fCentBin),mcore, ptcore) ;
        if(a<kAlphaCut){
          FillHistogram(Form("hPi0CPV2_a07_cen%d",fCentBin),m,pt) ;
        }
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCentBin) ;
	FillHistogram(Form("hMassPtV0ADisp_cen%d",fCentBin),m,pt,dphiA) ;
	FillHistogram(Form("hMassPtV0CDisp_cen%d",fCentBin),m,pt,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCDisp_cen%d",fCentBin),m,pt,dphiT) ;
	
	FillHistogram(Form("hMassPtV0ADispcore_cen%d",fCentBin),mcore, ptcore,dphiA) ;
	FillHistogram(Form("hMassPtV0CDispcore_cen%d",fCentBin),mcore, ptcore,dphiC) ;
	if(fHaveTPCRP)
	  FillHistogram(Form("hMassPtTPCDispcore_cen%d",fCentBin),mcore, ptcore,dphiT) ;

	FillHistogram(Form("hPi0Disp_cen%d",fCentBin),m,pt) ;
	FillHistogram(Form("hPi0Dispcore_cen%d",fCentBin),mcore, ptcore) ;
	
	if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	  FillHistogram(Form("hPi0Dispwou_cen%d",fCentBin),m,pt) ;

	  FillHistogram(Form("hMassPtV0ADispwou_cen%d",fCentBin),m,pt,dphiA) ;
 	  FillHistogram(Form("hMassPtV0CDispwou_cen%d",fCentBin),m,pt,dphiC) ;
	  if(fHaveTPCRP)
  	    FillHistogram(Form("hMassPtTPCDispwou_cen%d",fCentBin),m,pt,dphiT) ;
	}

        if(a<kAlphaCut){
          FillHistogram(Form("hPi0Disp_a07_cen%d",fCentBin),m,pt) ;
        }
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtV0ABoth_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMassPtV0CBoth_cen%d",fCentBin),m,pt,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth_cen%d",fCentBin),m,pt,dphiT) ;

	  FillHistogram(Form("hMassPtV0ABothcore_cen%d",fCentBin),mcore,ptcore,dphiA) ;
	  FillHistogram(Form("hMassPtV0CBothcore_cen%d",fCentBin),mcore,ptcore,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBothcore_cen%d",fCentBin),mcore,ptcore,dphiT) ;

	  FillHistogram(Form("hPi0Both_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hPi0Bothcore_cen%d",fCentBin),mcore,ptcore) ;

          if(a<kAlphaCut){
            snprintf(key,55,"hPi0Both_a07_cen%d",fCentBin) ;
            FillHistogram(Form("hPi0Both_a07_cen%d",fCentBin),m,pt) ;
          }
          if(ph1->Module()==1 && ph2->Module()==1)
	    FillHistogram("hPi0M11",m,pt );
          else if(ph1->Module()==2 && ph2->Module()==2)
	    FillHistogram("hPi0M22",m,pt );
          else if(ph1->Module()==3 && ph2->Module()==3)
	    FillHistogram("hPi0M33",m,pt );
          else if(ph1->Module()==4 && ph2->Module()==4)
	    FillHistogram("hPi0M44",m,pt );
          else if(ph1->Module()==1 && ph2->Module()==2)
	    FillHistogram("hPi0M12",m,pt );
          else if(ph1->Module()==1 && ph2->Module()==3)
	    FillHistogram("hPi0M13",m,pt );
          else if(ph1->Module()==1 && ph2->Module()==4)
	    FillHistogram("hPi0M14",m,pt );
          else if(ph1->Module()==2 && ph2->Module()==3)
	    FillHistogram("hPi0M23",m,pt );
          else if(ph1->Module()==2 && ph2->Module()==4)
	    FillHistogram("hPi0M24",m,pt );
          else if(ph1->Module()==3 && ph2->Module()==4)
	    FillHistogram("hPi0M34",m,pt );
        }
	
      }
      
      
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hPi0Disp2_cen%d",fCentBin),m,pt) ;
  	FillHistogram(Form("hPi0Disp2core_cen%d",fCentBin),mcore, ptcore) ;	

	FillHistogram(Form("hMassPtV0ADisp2_cen%d",fCentBin),m,pt,dphiA) ;
	FillHistogram(Form("hMassPtV0CDisp2_cen%d",fCentBin),m,pt,dphiC) ;
	if(fHaveTPCRP)
  	  FillHistogram(Form("hMassPtTPCDisp2_cen%d",fCentBin),m,pt,dphiT) ;

	FillHistogram(Form("hMassPtV0ADisp2core_cen%d",fCentBin),mcore, ptcore,dphiA) ;
	FillHistogram(Form("hMassPtV0CDisp2core_cen%d",fCentBin),mcore, ptcore,dphiC) ;
	if(fHaveTPCRP)
	  FillHistogram(Form("hMassPtTPCDisp2core_cen%d",fCentBin),mcore, ptcore,dphiT) ;
	  
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtV0ABoth2_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMassPtV0CBoth2_cen%d",fCentBin),m,pt,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth2_cen%d",fCentBin),m,pt,dphiT) ;

	  FillHistogram(Form("hMassPtV0ABoth2core_cen%d",fCentBin),mcore,ptcore,dphiA) ;
	  FillHistogram(Form("hMassPtV0CBoth2core_cen%d",fCentBin),mcore,ptcore,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMassPtTPCBoth2core_cen%d",fCentBin),mcore,ptcore,dphiT) ;

	  FillHistogram(Form("hPi0Both2_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hPi0Both2core_cen%d",fCentBin),mcore,ptcore) ;
	}

      }
    } // end of loop i2
  } // end of loop i1
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::ConsiderPi0sMix()
{
  char key[55];

  TList * arrayList = GetCaloPhotonsPHOSList(fVtxBin, fCentBin, fEMRPBin);

  for (Int_t i1=0; i1<fCaloPhotonsPHOS->GetEntriesFast(); i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fCaloPhotonsPHOS->At(i1) ;
    for(Int_t evi=0; evi<arrayList->GetEntries();evi++){
      TObjArray * mixPHOS = static_cast<TObjArray*>(arrayList->At(evi));
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	TLorentzVector p12  = *ph1  + *ph2;
	TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());

	Double_t dphiA=p12.Phi()-fRPV0A ;
	while(dphiA<0)dphiA+=TMath::Pi() ;
	while(dphiA>TMath::Pi())dphiA-=TMath::Pi() ;

	Double_t dphiC=p12.Phi()-fRPV0C ;
	while(dphiC<0)dphiC+=TMath::Pi() ;
	while(dphiC>TMath::Pi())dphiC-=TMath::Pi() ;

	Double_t dphiT=p12.Phi()-fRP ;
	while(dphiT<0)dphiT+=TMath::Pi() ;
	while(dphiT>TMath::Pi())dphiT-=TMath::Pi() ;


        Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
        Double_t m=p12.M() ;
        Double_t mcore=pv12.M() ;
        Double_t pt=p12.Pt() ;
        Double_t ptcore=pv12.Pt() ;
        Double_t pt1=ph1->Pt() ;
        Double_t pt2=ph2->Pt() ;
        Double_t ptcore1=ph1->GetMomV2()->Pt() ;
        Double_t ptcore2=ph2->GetMomV2()->Pt() ;

	snprintf(key,55,"hMiMassPtAll_cen%d",fCentBin) ; // probably not needed, consider removing this line!
	if( fFillWideTOF ) {
	  FillHistogram(Form("hMiPi0WideTOF_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hMiSingleWideTOF_cen%d",fCentBin),m,pt1) ;
	  FillHistogram(Form("hMiSingleWideTOF_cen%d",fCentBin),m,pt2) ;
	  if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCWideTOF_cen%d",fCentBin),m,pt,dphiT) ;
	}

	if( fTOFCutEnabled && !(ph1->IsTOFOK() && ph2->IsTOFOK()) )
	  continue;

	FillHistogram(Form("hMiMassPtV0AAll_cen%d",fCentBin),m,pt,dphiA) ;
	FillHistogram(Form("hMiMassPtV0CAll_cen%d",fCentBin),m,pt,dphiC) ;
	if(fHaveTPCRP)
 	  FillHistogram(Form("hMiMassPtTPCAll_cen%d",fCentBin),m,pt,dphiT) ;

	FillHistogram(Form("hMiMassPtV0AAllcore_cen%d",fCentBin),mcore, ptcore, dphiA) ;
	FillHistogram(Form("hMiMassPtV0CAllcore_cen%d",fCentBin),mcore, ptcore, dphiC) ;
        if(fHaveTPCRP)
	  FillHistogram(Form("hMiMassPtTPCAllcore_cen%d",fCentBin),mcore, ptcore, dphiT) ;

	FillHistogram(Form("hMiPi0All_cen%d",fCentBin),m,pt) ;
	FillHistogram(Form("hMiPi0Allcore_cen%d",fCentBin),mcore,ptcore) ;
	if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	  FillHistogram(Form("hMiPi0Allwou_cen%d",fCentBin),m,pt) ;
          FillHistogram(Form("hMiMassPtV0AAllwou_cen%d",fCentBin),m,pt,dphiA) ;
          FillHistogram(Form("hMiMassPtV0CAllwou_cen%d",fCentBin),m,pt,dphiC) ;
          if(fHaveTPCRP)
            FillHistogram(Form("hMiMassPtTPCAllwou_cen%d",fCentBin),m,pt,dphiT) ;
	}

	FillHistogram(Form("hMiSingleAll_cen%d",fCentBin),m,pt1) ;
        FillHistogram(Form("hMiSingleAll_cen%d",fCentBin),m,pt2) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCentBin),mcore,ptcore1) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCentBin),mcore,ptcore2) ;
        if(ph1->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCentBin),m,pt1) ;
        if(ph2->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCentBin),m,pt2) ;
        if(ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCentBin),m,pt1) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCentBin),mcore,ptcore1) ;
        }
        if(ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCentBin),m,pt2) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCentBin),mcore,ptcore2) ;
        }
        if(ph1->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCentBin),m,pt1) ;
          FillHistogram(Form("hMiSingleCPV2core_cen%d",fCentBin),mcore,ptcore1) ;
        }
        if(ph2->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCentBin),m,pt2) ;
          FillHistogram(Form("hMiSingleCPV2core_cen%d",fCentBin),mcore,ptcore2) ;
        }
        if(ph1->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCentBin),m,pt1) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCentBin),m,pt1) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCentBin),mcore,ptcore1) ;
        }
        if(ph2->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCentBin),m,pt2) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCentBin),m,pt2) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCentBin),mcore,ptcore2) ;
        }
        if(ph1->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCentBin),m,pt1) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCentBin),mcore,ptcore1) ;
        }
        if(ph2->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCentBin),m,pt2) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCentBin),mcore,ptcore2) ;
        }
        if(ph1->IsDispOK() && ph1->IsCPVOK()){
          snprintf(key,55,"hMiSingleBoth_cen%d",fCentBin) ;
          FillHistogram(key,m,pt1) ;
          snprintf(key,55,"hMiSingleBothcore_cen%d",fCentBin) ;
          FillHistogram(key,mcore,ptcore1) ;
        }
        if(ph2->IsDispOK() && ph2->IsCPVOK()){
          snprintf(key,55,"hMiSingleBoth_cen%d",fCentBin) ;
          FillHistogram(key,m,pt2) ;
          snprintf(key,55,"hMiSingleBothcore_cen%d",fCentBin) ;
          FillHistogram(key,mcore,ptcore2) ;
        }
        if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCentBin),m,pt1) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCentBin),mcore,ptcore1) ;
        }
        if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCentBin),m,pt2) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCentBin),mcore,ptcore2) ;
        }



        if(a<kAlphaCut){
          FillHistogram(Form("hMiPi0All_a07_cen%d",fCentBin),m,pt) ;
        }
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMiMassPtV0ACPV_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CCPV_cen%d",fCentBin),m,pt,dphiC) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV_cen%d",fCentBin),m,pt,dphiT) ;

	  FillHistogram(Form("hMiMassPtV0ACPVcore_cen%d",fCentBin),mcore, ptcore,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CCPVcore_cen%d",fCentBin),mcore, ptcore,dphiC) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPVcore_cen%d",fCentBin),mcore, ptcore,dphiT) ;

	  FillHistogram(Form("hMiPi0CPV_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hMiPi0CPVcore_cen%d",fCentBin),mcore, ptcore) ;

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0CPV_a07_cen%d",fCentBin),m,pt) ;
          }
	}
	if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	  FillHistogram(Form("hMiPi0CPV2_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hMiPi0CPV2core_cen%d",fCentBin),mcore, ptcore) ;

	  FillHistogram(Form("hMiMassPtV0ACPV2_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CCPV2_cen%d",fCentBin),m,pt,dphiC) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV2_cen%d",fCentBin),m,pt,dphiT) ;
	  FillHistogram(Form("hMiMassPtV0ACPV2core_cen%d",fCentBin),mcore,ptcore,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CCPV2core_cen%d",fCentBin),mcore,ptcore,dphiC) ;
	  if(fHaveTPCRP)
 	    FillHistogram(Form("hMiMassPtTPCCPV2core_cen%d",fCentBin),mcore,ptcore,dphiT) ;

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0CPV2_a07_cen%d",fCentBin),m,pt) ;
          }
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  FillHistogram(Form("hMiMassPtV0ADisp_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CDisp_cen%d",fCentBin),m,pt,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp_cen%d",fCentBin),m,pt,dphiT) ;

	  FillHistogram(Form("hMiMassPtV0ADispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CDispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDispcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT) ;


	  FillHistogram(Form("hMiPi0Disp_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hMiPi0Dispcore_cen%d",fCentBin),pv12.M(),pv12.Pt()) ;
          if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	    FillHistogram(Form("hMiPi0Dispwou_cen%d",fCentBin),m,pt) ;
	    FillHistogram(Form("hMiMassPtV0ADispwou_cen%d",fCentBin),m,pt,dphiA) ;
	    FillHistogram(Form("hMiMassPtV0CDispwou_cen%d",fCentBin),m,pt,dphiC) ;
            if(fHaveTPCRP)
	      FillHistogram(Form("hMiMassPtTPCDispwou_cen%d",fCentBin),m,pt,dphiT) ;
	  }

	  if(a<kAlphaCut){
            FillHistogram(Form("hMiPi0Disp_a07_cen%d",fCentBin),m,pt) ;
          }
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtV0ABoth_cen%d",fCentBin),m,pt,dphiA) ;
	    FillHistogram(Form("hMiMassPtV0CBoth_cen%d",fCentBin),m,pt,dphiC) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth_cen%d",fCentBin),m,pt,dphiT) ;

	    FillHistogram(Form("hMiMassPtV0ABothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA) ;
	    FillHistogram(Form("hMiMassPtV0CBothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBothcore_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT) ;

	    FillHistogram(Form("hMiPi0Both_cen%d",fCentBin),m,pt) ;
	    FillHistogram(Form("hMiPi0Bothcore_cen%d",fCentBin),pv12.M(),pv12.Pt()) ;

	    if(a<kAlphaCut){
              FillHistogram(Form("hMiPi0Both_a07_cen%d",fCentBin),m,pt) ;
            }
	  }
	}
	
  	if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	  FillHistogram(Form("hMiMassPtV0ADisp2_cen%d",fCentBin),m,pt,dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CDisp2_cen%d",fCentBin),m,pt,dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp2_cen%d",fCentBin),m,pt,dphiT) ;

	  FillHistogram(Form("hMiMassPtV0ADisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA) ;
	  FillHistogram(Form("hMiMassPtV0CDisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC) ;
          if(fHaveTPCRP)
	    FillHistogram(Form("hMiMassPtTPCDisp2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT) ;


	  FillHistogram(Form("hMiPi0Disp2_cen%d",fCentBin),m,pt) ;
	  FillHistogram(Form("hMiPi0Disp2core_cen%d",fCentBin),pv12.M(),pv12.Pt()) ;

	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtV0ABoth2_cen%d",fCentBin),m,pt,dphiA) ;
	    FillHistogram(Form("hMiMassPtV0CBoth2_cen%d",fCentBin),m,pt,dphiC) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth2_cen%d",fCentBin),m,pt,dphiT) ;

	    FillHistogram(Form("hMiMassPtV0ABoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiA) ;
	    FillHistogram(Form("hMiMassPtV0CBoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiC) ;
	    if(fHaveTPCRP)
  	      FillHistogram(Form("hMiMassPtTPCBoth2core_cen%d",fCentBin),pv12.M(),pv12.Pt(),dphiT) ;

	    FillHistogram(Form("hMiPi0Both2_cen%d",fCentBin),m,pt) ;
	    FillHistogram(Form("hMiPi0Both2core_cen%d",fCentBin),pv12.M(),pv12.Pt()) ;

	  }
	}
      } // end of loop i2
    }
  } // end of loop i1
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::UpdateLists()
{
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed

  TList * arrayList = GetCaloPhotonsPHOSList(fVtxBin, fCentBin, fEMRPBin);
  if( fDebug >= 2 )
    AliInfo( Form("fCentBin=%d, fCentNMixed[]=%d",fCentBin,fCentNMixed[fCentBin]) );
  if(fCaloPhotonsPHOS->GetEntriesFast()>0){
    arrayList->AddFirst(fCaloPhotonsPHOS) ;
    fCaloPhotonsPHOS=0;
    if(arrayList->GetEntries() > fCentNMixed[fCentBin]){ // Remove redundant events
      TObjArray * tmp = static_cast<TObjArray*>(arrayList->Last()) ;
      arrayList->RemoveLast() ;
      delete tmp ; // TODO: may conflict with delete done by list being owner.
    }
  }
  else
    fCaloPhotonsPHOS->Clear(); // TODO: redundant???
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
  if(th1)
    th1->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * obj = fOutputContainer->FindObject(key);
  
  TH2 * th2 = dynamic_cast<TH2*> (obj);
  if(th2) {
    th2->Fill(x, y, z) ;
    return;
  }

  TH3 * th3 = dynamic_cast<TH3*> (obj);
  if(th3) {
    th3->Fill(x, y, z) ;
    return;
  }
  
  AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}

void AliAnalysisTaskPi0Flow::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const{
  //Fills 1D histograms with key
  TObject * obj = fOutputContainer->FindObject(key);
  
  TH3 * th3 = dynamic_cast<TH3*> (obj);
  if(th3) {
    th3->Fill(x, y, z, w) ;
    return;
  }
  
  AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}


//_____________________________________________________________________________
AliVEvent* AliAnalysisTaskPi0Flow::GetEvent()
{
  fEvent = InputEvent();
  if( ! fEvent ) {
    AliError("Event could not be retrieved");
    PostData(1, fOutputContainer);
  }
  return fEvent;
}


//___________________________________________________________________________
// AliStack* AliAnalysisTaskPi0Flow::GetMCStack()
// {
//   fMCStack = 0;
//   AliVEventHandler* eventHandler = AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler();
//   if(eventHandler){
//     AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*> (eventHandler);
//     if( mcEventHandler)
//       fMCStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
//   }
//   return fMCStack;
// }

//___________________________________________________________________________
Int_t AliAnalysisTaskPi0Flow::GetCentralityBin(Float_t centralityV0M)
{
 /* fCentBin=1+Int_t(centralityV0M/100. *kNCenBins) ;
  if(centralityV0M < 5. || fCentBin < 0)
   fCentBin=0 ;
  if(fCentBin > kNCenBins-1)
    fCentBin = kNCenBins-1 ;
 */
  int lastBinUpperIndex = fCentEdges.GetSize() -1;
  if( centralityV0M > fCentEdges[lastBinUpperIndex] ) {
    if( fDebug >= 1 )
      AliWarning( Form("centrality (%f) larger then upper edge of last centrality bin (%f)!", centralityV0M, fCentEdges[lastBinUpperIndex]) );
    return lastBinUpperIndex-1;
  }
  if( centralityV0M < fCentEdges[0] ) {
    if( fDebug >= 1 )
      AliWarning( Form("centrality (%f) smaller then lower edge of first bin (%f)!", centralityV0M, fCentEdges[0]) );
    return 0;
  }
  
  fCentBin = TMath::BinarySearch<Double_t> ( GetNumberOfCentralityBins(), fCentEdges.GetArray(), centralityV0M );
  return fCentBin;
}

//___________________________________________________________________________
Int_t AliAnalysisTaskPi0Flow::GetRPBin()
{
  Double_t averageRP;
  if(fHaveTPCRP)
    averageRP = fRP ;// If possible, it is better to have EP bin from TPC
                     // to have similar events for miximng (including jets etc)   (fRPV0A+fRPV0C+fRP) /3.;
  else
    averageRP = (fRPV0A+fRPV0C) /2.;

  fEMRPBin = Int_t(fNEMRPBins*(averageRP)/TMath::Pi());

  if(fEMRPBin> (Int_t) fNEMRPBins-1)
    fEMRPBin=fNEMRPBins-1 ;
  else if(fEMRPBin<0)
    fEMRPBin=0;

  if ( fDebug >= 2 )
    AliInfo(Form("Event Mixing Reaction Plane bin is: %d", fEMRPBin));

  return fEMRPBin;
}


//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::LogProgress(int step)
{
  if(fDebug >= 2) {
    AliInfo(Form("step %d completed", step));
  }
  // the +0.5 is not realy neccisarry, but oh well... -henrik
  //FillHistogram("hSelEvents", step+0.5, internalRunNumber-0.5);
  //FillHistogram("hTotSelEvents", step+0.5);
}

void AliAnalysisTaskPi0Flow::LogSelection(int step, int internalRunNumber)
{
  // if(fDebug > 1) {
  //   AliInfo(Form("step %d completed", step));
  // }
  // the +0.5 is not realy neccisarry, but oh well... -henrik
  FillHistogram("hSelEvents", step+0.5, internalRunNumber-0.5);
  FillHistogram("hTotSelEvents", step+0.5);
}


//___________________________________________________________________________
Int_t AliAnalysisTaskPi0Flow::ConvertToInternalRunNumber(Int_t run){
  if( kLHC11h == fPeriod ) {
    switch(run){
    case  170593 : return 179 ;
    case  170572 : return 178 ;
    case  170556 : return 177 ;
    case  170552 : return 176 ;
    case  170546 : return 175 ;
    case  170390 : return 174 ;
    case  170389 : return 173 ;
    case  170388 : return 172 ;
    case  170387 : return 171 ;
    case  170315 : return 170 ;
    case  170313 : return 169 ;
    case  170312 : return 168 ;
    case  170311 : return 167 ;
    case  170309 : return 166 ;
    case  170308 : return 165 ;
    case  170306 : return 164 ;
    case  170270 : return 163 ;
    case  170269 : return 162 ;
    case  170268 : return 161 ;
    case  170267 : return 160 ;
    case  170264 : return 159 ;
    case  170230 : return 158 ;
    case  170228 : return 157 ;
    case  170208 : return 156 ;
    case  170207 : return 155 ;
    case  170205 : return 154 ;
    case  170204 : return 153 ;
    case  170203 : return 152 ;
    case  170195 : return 151 ;
    case  170193 : return 150 ;
    case  170163 : return 149 ;
    case  170162 : return 148 ;
    case  170159 : return 147 ;
    case  170155 : return 146 ;
    case  170152 : return 145 ;
    case  170091 : return 144 ;
    case  170089 : return 143 ;
    case  170088 : return 142 ;
    case  170085 : return 141 ;
    case  170084 : return 140 ;
    case  170083 : return 139 ;
    case  170081 : return 138 ;
    case  170040 : return 137 ;
    case  170038 : return 136 ;
    case  170036 : return 135 ;
    case  170027 : return 134 ;
    case  169981 : return 133 ;
    case  169975 : return 132 ;
    case  169969 : return 131 ;
    case  169965 : return 130 ;
    case  169961 : return 129 ;
    case  169956 : return 128 ;
    case  169926 : return 127 ;
    case  169924 : return 126 ;
    case  169923 : return 125 ;
    case  169922 : return 124 ;
    case  169919 : return 123 ;
    case  169918 : return 122 ;
    case  169914 : return 121 ;
    case  169859 : return 120 ;
    case  169858 : return 119 ;
    case  169855 : return 118 ;
    case  169846 : return 117 ;
    case  169838 : return 116 ;
    case  169837 : return 115 ;
    case  169835 : return 114 ;
    case  169683 : return 113 ;
    case  169628 : return 112 ;
    case  169591 : return 111 ;
    case  169590 : return 110 ;
    case  169588 : return 109 ;
    case  169587 : return 108 ;
    case  169586 : return 107 ;
    case  169584 : return 106 ;
    case  169557 : return 105 ;
    case  169555 : return 104 ;
    case  169554 : return 103 ;
    case  169553 : return 102 ;
    case  169550 : return 101 ;
    case  169515 : return 100 ;
    case  169512 : return 99 ;
    case  169506 : return 98 ;
    case  169504 : return 97 ;
    case  169498 : return 96 ;
    case  169475 : return 95 ;
    case  169420 : return 94 ;
    case  169419 : return 93 ;
    case  169418 : return 92 ;
    case  169417 : return 91 ;
    case  169415 : return 90 ;
    case  169411 : return 89 ;
    case  169238 : return 88 ;
    case  169236 : return 87 ;
    case  169167 : return 86 ;
    case  169160 : return 85 ;
    case  169156 : return 84 ;
    case  169148 : return 83 ;
    case  169145 : return 82 ;
    case  169144 : return 81 ;
    case  169143 : return 80 ;
    case  169138 : return 79 ;
    case  169099 : return 78 ;
    case  169094 : return 77 ;
    case  169091 : return 76 ;
    case  169045 : return 75 ;
    case  169044 : return 74 ;
    case  169040 : return 73 ;
    case  169035 : return 72 ;
    case  168992 : return 71 ;
    case  168988 : return 70 ;
    case  168984 : return 69 ;
    case  168826 : return 68 ;
    case  168777 : return 67 ;
    case  168514 : return 66 ;
    case  168512 : return 65 ;
    case  168511 : return 64 ;
    case  168467 : return 63 ;
    case  168464 : return 62 ;
    case  168461 : return 61 ;
    case  168460 : return 60 ;
    case  168458 : return 59 ;
    case  168362 : return 58 ;
    case  168361 : return 57 ;
    case  168356 : return 56 ;
    case  168342 : return 55 ;
    case  168341 : return 54 ;
    case  168325 : return 53 ;
    case  168322 : return 52 ;
    case  168318 : return 51 ;
    case  168311 : return 50 ;
    case  168310 : return 49 ;
    case  168213 : return 48 ;
    case  168212 : return 47 ;
    case  168208 : return 46 ;
    case  168207 : return 45 ;
    case  168206 : return 44 ;
    case  168205 : return 43 ;
    case  168204 : return 42 ;
    case  168203 : return 41 ;
    case  168181 : return 40 ;
    case  168177 : return 39 ;
    case  168175 : return 38 ;
    case  168173 : return 37 ;
    case  168172 : return 36 ;
    case  168171 : return 35 ;
    case  168115 : return 34 ;
    case  168108 : return 33 ;
    case  168107 : return 32 ;
    case  168105 : return 31 ;
    case  168104 : return 30 ;
    case  168103 : return 29 ;
    case  168076 : return 28 ;
    case  168069 : return 27 ;
    case  168068 : return 26 ;
    case  168066 : return 25 ;
    case  167988 : return 24 ;
    case  167987 : return 23 ;
    case  167986 : return 22 ;
    case  167985 : return 21 ;
    case  167921 : return 20 ;
    case  167920 : return 19 ;
    case  167915 : return 18 ;
    case  167909 : return 17 ;
    case  167903 : return 16 ;
    case  167902 : return 15 ;
    case  167818 : return 14 ;
    case  167814 : return 13 ;
    case  167813 : return 12 ;
    case  167808 : return 11 ;
    case  167807 : return 10 ;
    case  167806 : return 9 ;
    case  167713 : return 8 ;
    case  167712 : return 7 ;
    case  167711 : return 6 ;
    case  167706 : return 5 ;
    case  167693 : return 4 ;
    case  166532 : return 3 ;
    case  166530 : return 2 ;
    case  166529 : return 1 ;

    default : return 199;
    }
  }
  if( kLHC10h == fPeriod ) {
    switch(run){
    case  139517 : return 137;
    case  139514 : return 136;
    case  139513 : return 135;
    case  139511 : return 134;
    case  139510 : return 133;
    case  139507 : return 132;
    case  139505 : return 131;
    case  139504 : return 130;
    case  139503 : return 129;
    case  139470 : return 128;
    case  139467 : return 127;
    case  139466 : return 126;
    case  139465 : return 125;
    case  139440 : return 124;
    case  139439 : return 123;
    case  139438 : return 122;
    case  139437 : return 121;
    case  139360 : return 120;
    case  139329 : return 119;
    case  139328 : return 118;
    case  139314 : return 117;
    case  139311 : return 116;
    case  139310 : return 115;
    case  139309 : return 114;
    case  139308 : return 113;
    case  139173 : return 112;
    case  139172 : return 111;
    case  139110 : return 110;
    case  139107 : return 109;
    case  139105 : return 108;
    case  139104 : return 107;
    case  139042 : return 106;
    case  139038 : return 105;
    case  139037 : return 104;
    case  139036 : return 103;
    case  139029 : return 102;
    case  139028 : return 101;
    case  138983 : return 100;
    case  138982 : return 99;
    case  138980 : return 98;
    case  138979 : return 97;
    case  138978 : return 96;
    case  138977 : return 95;
    case  138976 : return 94;
    case  138973 : return 93;
    case  138972 : return 92;
    case  138965 : return 91;
    case  138924 : return 90;
    case  138872 : return 89;
    case  138871 : return 88;
    case  138870 : return 87;
    case  138837 : return 86;
    case  138830 : return 85;
    case  138828 : return 84;
    case  138826 : return 83;
    case  138796 : return 82;
    case  138795 : return 81;
    case  138742 : return 80;
    case  138732 : return 79;
    case  138730 : return 78;
    case  138666 : return 77;
    case  138662 : return 76;
    case  138653 : return 75;
    case  138652 : return 74;
    case  138638 : return 73;
    case  138624 : return 72;
    case  138621 : return 71;
    case  138583 : return 70;
    case  138582 : return 69;
    case  138579 : return 68;
    case  138578 : return 67;
    case  138534 : return 66;
    case  138469 : return 65;
    case  138442 : return 64;
    case  138439 : return 63;
    case  138438 : return 62;
    case  138396 : return 61;
    case  138364 : return 60;
    case  138359 : return 59;
    case  138275 : return 58;
    case  138225 : return 57;
    case  138201 : return 56;
    case  138200 : return 55;
    case  138197 : return 54;
    case  138192 : return 53;
    case  138190 : return 52;
    case  138154 : return 51;
    case  138153 : return 50;
    case  138151 : return 49;
    case  138150 : return 48;
    case  138126 : return 47;
    case  138125 : return 46;
    case  137848 : return 45;
    case  137847 : return 44;
    case  137844 : return 43;
    case  137843 : return 42;
    case  137752 : return 41;
    case  137751 : return 40;
    case  137748 : return 39;
    case  137724 : return 38;
    case  137722 : return 37;
    case  137718 : return 36;
    case  137704 : return 35;
    case  137693 : return 34;
    case  137692 : return 33;
    case  137691 : return 32;
    case  137689 : return 31;
    case  137686 : return 30;
    case  137685 : return 29;
    case  137639 : return 28;
    case  137638 : return 27;
    case  137608 : return 26;
    case  137595 : return 25;
    case  137549 : return 24;
    case  137546 : return 23;
    case  137544 : return 22;
    case  137541 : return 21;
    case  137539 : return 20;
    case  137531 : return 19;
    case  137530 : return 18;
    case  137443 : return 17;
    case  137441 : return 16;
    case  137440 : return 15;
    case  137439 : return 14;
    case  137434 : return 13;
    case  137432 : return 12;
    case  137431 : return 11;
    case  137430 : return 10;
    case  137366 : return 9;
    case  137243 : return 8;
    case  137236 : return 7;
    case  137235 : return 6;
    case  137232 : return 5;
    case  137231 : return 4;
    case  137165 : return 3;
    case  137162 : return 2;
    case  137161 : return 1;
    default : return 199;
    }
  }
  if( kLHC13 == fPeriod ) {
    switch(run){
    case  195344 : return 1;
    case  195346 : return 2;
    case  195351 : return 3;
    case  195389 : return 4;
    case  195390 : return 5;
    case  195391 : return 6;
    case  195478 : return 7;
    case  195479 : return 8;
    case  195480 : return 9;
    case  195481 : return 10;
    case  195482 : return 11;
    case  195483 : return 12;
    case  195529 : return 13;
    case  195531 : return 14;
    case  195532 : return 15;
    case  195566 : return 16;
    case  195567 : return 17;
    case  195568 : return 18;
    case  195592 : return 19;
    case  195593 : return 20;
    case  195596 : return 21;
    case  195633 : return 22;
    case  195635 : return 23;
    case  195644 : return 24;
    case  195673 : return 25;
    case  195675 : return 26;
    case  195676 : return 27;
    case  195677 : return 28;
    case  195681 : return 29;
    case  195682 : return 30;
    case  195720 : return 31;
    case  195721 : return 32;
    case  195722 : return 33;
    case  195724 : return 34;
    case  195725 : return 34;
    case  195726 : return 35;
    case  195727 : return 36;
    case  195760 : return 37;
    case  195761 : return 38;
    case  195765 : return 39;
    case  195767 : return 40;
    case  195783 : return 41;
    case  195787 : return 42;
    case  195826 : return 43;
    case  195827 : return 44;
    case  195829 : return 45;
    case  195830 : return 46;
    case  195831 : return 47;
    case  195867 : return 48;
    case  195869 : return 49;
    case  195871 : return 50;
    case  195872 : return 51;
    case  195873 : return 52;
    case  195935 : return 53;
    case  195949 : return 54;
    case  195950 : return 55;
    case  195954 : return 56;
    case  195955 : return 57;
    case  195958 : return 58;
    case  195989 : return 59;
    case  195994 : return 60;
    case  195998 : return 61;
    case  196000 : return 62;
    case  196006 : return 63;
    case  196085 : return 64;
    case  196089 : return 65;
    case  196090 : return 66;
    case  196091 : return 67;
    case  196099 : return 68;
    case  196105 : return 69;
    case  196107 : return 70;
    case  196185 : return 71;
    case  196187 : return 72;
    case  196194 : return 73;
    case  196197 : return 74;
    case  196199 : return 75;
    case  196200 : return 76;
    case  196201 : return 77;
    case  196203 : return 78;
    case  196208 : return 79;
    case  196214 : return 80;
    case  196308 : return 81;
    case  196309 : return 82;
    case  196310 : return 83;
    case  196311 : return 84;
    case  196433 : return 85;
    case  196474 : return 86;
    case  196475 : return 87;
    case  196477 : return 88;
    case  196528 : return 89;
    case  196533 : return 90;
    case  196535 : return 91;
    case  196563 : return 92;
    case  196564 : return 93;
    case  196566 : return 94;
    case  196568 : return 95;
    case  196601 : return 96;
    case  196605 : return 97;
    case  196608 : return 98;
    case  196646 : return 99;
    case  196648 : return 100;
    case  196701 : return 101;
    case  196702 : return 102;
    case  196703 : return 103;
    case  196706 : return 104;
    case  196714 : return 105;
    case  196720 : return 106;
    case  196721 : return 107;
    case  196722 : return 108;
    case  196772 : return 109;
    case  196773 : return 110;
    case  196774 : return 111;
    case  196869 : return 112;
    case  196870 : return 113;
    case  196874 : return 114;
    case  196876 : return 115;
    case  196965 : return 116;
    case  196967 : return 117;
    case  196972 : return 118;
    case  196973 : return 119;
    case  196974 : return 120;
    case  197003 : return 121;
    case  197011 : return 122;
    case  197012 : return 123;
    case  197015 : return 124;
    case  197027 : return 125;
    case  197031 : return 126;
    case  197089 : return 127;
    case  197090 : return 128;
    case  197091 : return 129;
    case  197092 : return 130;
    case  197094 : return 131;
    case  197098 : return 132;
    case  197099 : return 133;
    case  197138 : return 134;
    case  197139 : return 135;
    case  197142 : return 136;
    case  197143 : return 137;
    case  197144 : return 138;
    case  197145 : return 139;
    case  197146 : return 140;
    case  197147 : return 140;
    case  197148 : return 141;
    case  197149 : return 142;
    case  197150 : return 143;
    case  197152 : return 144;
    case  197153 : return 145;
    case  197184 : return 146;
    case  197189 : return 147;
    case  197247 : return 148;
    case  197248 : return 149;
    case  197254 : return 150;
    case  197255 : return 151;
    case  197256 : return 152;
    case  197258 : return 153;
    case  197260 : return 154;
    case  197296 : return 155;
    case  197297 : return 156;
    case  197298 : return 157;
    case  197299 : return 158;
    case  197300 : return 159;
    case  197302 : return 160;
    case  197341 : return 161;
    case  197342 : return 162;
    case  197348 : return 163;
    case  197349 : return 164;
    case  197351 : return 165;
    case  197386 : return 166;
    case  197387 : return 167;
    case  197388 : return 168;
    default : return 199;
    }
  }
  if((fPeriod == kUndefinedPeriod) && (fDebug >= 1) ) {
    AliWarning("Period not defined");
  }
  return 1;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::TestLambda(Double_t pt,Double_t l1,Double_t l2){

  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
    0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
    0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;

}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::TestLambda2(Double_t pt,Double_t l1,Double_t l2){

  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma +
    0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
    0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<1.5*1.5) ;

}
//____________________________________________________________________________
TList* AliAnalysisTaskPi0Flow::GetCaloPhotonsPHOSList(UInt_t vtxBin, UInt_t centBin, UInt_t rpBin)
{
  int offset = vtxBin * GetNumberOfCentralityBins() * fNEMRPBins
	      + centBin * fNEMRPBins
	      + rpBin;
  if( fCaloPhotonsPHOSLists->At(offset) ) { // list exists
    TList* list = dynamic_cast<TList*> (fCaloPhotonsPHOSLists->At(offset));
    if( ! list )
      AliError("object in fCaloPhotonsPHOSLists at %i did not cast");
    return list;
  }
  else {// no list for this bin has been created, yet
    TList* list = new TList();
    list->SetOwner();
    fCaloPhotonsPHOSLists->AddAt(list, offset);
    return list;
  }
}

//____________________________________________________________________________
Double_t AliAnalysisTaskPi0Flow::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_

  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
			 6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  Double_t mf = 0.; //
  if(fEventAOD) mf = fEventAOD->GetMagneticField(); //Positive for ++ and negative for -- 
  else if(fEventESD) mf = fEventESD->GetMagneticField(); //Positive for ++ and negative for --
  

  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ= -0.468318;
    if(charge>0)
      meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
}
//____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetFlatteningData(){
  //Read objects with flattening parameters 
  AliOADBContainer flatContainer("phosFlat");
  flatContainer.InitFromFile(fEPcalibFileName.Data(),"phosFlat");
  TObjArray *maps = (TObjArray*)flatContainer.GetObject(fRunNumber,"phosFlat");
  if(!maps){
      AliError(Form("Can not read Flattening for run %d. \n From file >%s<\n",fRunNumber,fEPcalibFileName.Data())) ;    
  }
  else{
    AliInfo(Form("Setting PHOS flattening with name %s \n",maps->GetName())) ;
    AliEPFlattener * h = (AliEPFlattener*)maps->At(0) ;  
    if(fTPCFlat) delete fTPCFlat ;
    fTPCFlat = new AliEPFlattener() ;
    fTPCFlat = h ;
    h = (AliEPFlattener*)maps->At(1) ;  
    if(fV0AFlat) delete fV0AFlat ;
    fV0AFlat = new AliEPFlattener() ;
    fV0AFlat = h ;
    h = (AliEPFlattener*)maps->At(2) ;  
    if(fV0CFlat) delete fV0CFlat ;
    fV0CFlat = new AliEPFlattener() ;
    fV0CFlat = h ;
  }    
  
}
 //____________________________________________________________________________
Double_t  AliAnalysisTaskPi0Flow::ApplyFlattening(Double_t phi, Double_t c){
  
  if(fTPCFlat)
    return fTPCFlat->MakeFlat(phi,c);
  return phi ;

}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0Flow::ApplyFlatteningV0A(Double_t phi, Double_t c){
  
  if(fV0AFlat)
    return fV0AFlat->MakeFlat(phi,c);
  return phi ;

}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0Flow::ApplyFlatteningV0C(Double_t phi, Double_t c){
  
  if(fV0CFlat)
    return fV0CFlat->MakeFlat(phi,c);
  return phi ;

}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0Flow::CoreEnergy(AliVCluster * clu, AliVCaloCells * cells)
{
  //calculate energy of the cluster in the circle with radius distanceCut around the maximum

  //Can not use already calculated coordinates?
  //They have incidence correction...
  const Double_t distanceCut =3.5 ;
  const Double_t logWeight=4.5 ;

  const Double32_t * elist = clu->GetCellsAmplitudeFraction() ;
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  const Int_t mulDigit=clu->GetNCells() ;
  Double_t xc[mulDigit] ;
  Double_t zc[mulDigit] ;
  Double_t ei[mulDigit] ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    ei[iDigit]=elist[iDigit]*cells->GetCellAmplitude(clu->GetCellsAbsId()[iDigit]);
    if( fDebug >= 3 )
      printf("%f ",ei[iDigit]);
    if (clu->E()>0 && ei[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( ei[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  Double_t coreE=0. ;
  for(Int_t iDigit=0; iDigit < mulDigit; iDigit++) {
    Double_t distance = TMath::Sqrt((xc[iDigit]-x)*(xc[iDigit]-x)+(zc[iDigit]-z)*(zc[iDigit]-z)) ;
    if(distance < distanceCut)
      coreE += ei[iDigit] ;
  }
  //Apply non-linearity correction
  return fNonLinCorr->Eval(coreE) ;
}
//____________________________________________________________________________
Bool_t  AliAnalysisTaskPi0Flow::AreNeibors(Int_t id1,Int_t id2){
  // return true if absId are "Neighbors" (adjacent, including diagornaly,)
  // false if not.

  Int_t relid1[4] ;
  fPHOSGeo->AbsToRelNumbering(id1, relid1) ;

  Int_t relid2[4] ;
  fPHOSGeo->AbsToRelNumbering(id2, relid2) ;

  // if inside the same PHOS module
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) {
    const Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;
    const Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;

    // and if diff in both direction is 1 or less
    if (( coldiff <= 1 )  && ( rowdiff <= 1 ))
      return true; // are neighbors
  }

  // else false
  return false;
}
//____________________________________________________________________________
void  AliAnalysisTaskPi0Flow::Reclusterize(AliVCluster * clu){
  //Re-clusterize to make continues cluster

  const Int_t oldMulDigit=clu->GetNCells() ;
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;
  UShort_t * dlist = clu->GetCellsAbsId();

  Int_t index[oldMulDigit] ;
  Bool_t used[oldMulDigit] ;
  for(Int_t i=0; i<oldMulDigit; i++) used[i]=0 ;
  Int_t inClu=0 ;
  Double_t eMax=0. ;
  //find maximum
  for(Int_t iDigit=0; iDigit<oldMulDigit; iDigit++) {
    if(eMax<elist[iDigit]){
      eMax=elist[iDigit];
      index[0]=iDigit ;
      inClu=1 ;
    }
  }
  if(inClu==0){ //empty cluster
    return ;
  }
  used[index[0]]=kTRUE ; //mark as used
  for(Int_t i=0; i<inClu; i++){
    for(Int_t iDigit=0 ;iDigit<oldMulDigit; iDigit++){
       if(used[iDigit]) //already used
         continue ;
       if(AreNeibors(dlist[index[i]],dlist[iDigit])){
	 index[inClu]= iDigit ;
	 inClu++ ;
	 used[iDigit]=kTRUE ;
       }
    }
  }

  if(inClu==oldMulDigit) //no need to modify
    return ;

  clu->SetNCells(inClu);
  //copy
  UShort_t tmpD[oldMulDigit] ;
  Double_t tmpE[oldMulDigit] ;
  for(Int_t i=0; i<oldMulDigit; i++){
    tmpD[i]=dlist[i] ;
    tmpE[i]=elist[i] ;
  }
  //change order of digits in list so that
  //first inClu cells were true ones
  for(Int_t i=0; i<inClu; i++){
    dlist[i]=tmpD[index[i]] ;
    elist[i]=tmpE[index[i]] ;
  }


}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetMisalignment(){
  // sets the misalignment vertex if ESD
  if( fEventESD ) {
    for(Int_t mod=0; mod<5; mod++) {
      const TGeoHMatrix* modMatrix = fEvent->GetPHOSMatrix(mod);
      if( ! modMatrix) {
	if( fDebug )
	  AliInfo(Form("no PHOS Geometric Misalignment Matrix for module %d", mod));
	continue;
      }
      else {
	fPHOSGeo->SetMisalMatrix(modMatrix, mod);
	if( fDebug )
	  AliInfo(Form("PHOS Geometric Misalignment Matrix set for module %d", mod));
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetV0Calibration(){
    // assigns: fMultV0, fV0Cpol, fV0Apol, fMeanQ, and fWidthQ

    if ( ! fManualV0EPCalc ) {
      if( fDebug >=2 )
	AliInfo("Not setting V0Calibration, only needed for manual V0 EP Calculation");
      return; 
    }

    int runNumber = this->fRunNumber;
    
    TString oadbfilename = "$ALICE_PHYSICS/OADB/PWGCF/VZERO/VZEROcalibEP.root";
    TFile *foadb = TFile::Open(oadbfilename.Data());

    if(!foadb){
	AliError(Form("OADB file %s cannot be opened\n", oadbfilename.Data()));
	AliError("V0 Calibration not set !\n");
	return;
    }

    AliOADBContainer *cont = (AliOADBContainer*) foadb->Get("hMultV0BefCorr");
    if(!cont){
	AliError("OADB object hMultV0BefCorr is not available in the file");
	AliError("V0 Calibration not set!\n");
	return;
    }

    if(!(cont->GetObject(runNumber))){
	AliError(Form("OADB object hMultV0BefCorr is not available for run %i, trying 137366)",runNumber));
	runNumber = 137366;
    }
    if(!(cont->GetObject(runNumber))){
	AliError(Form("OADB object hMultV0BefCorr is not available for run %i ",runNumber));
	AliError("V0 Calibration not set!\n");
	return;
    }

    if( fDebug )  AliInfo("Setting V0 calibration") ;
    fMultV0 = ((TH2F *) cont->GetObject(runNumber))->ProfileX();

    TF1 *fpol0 = new TF1("fpol0","pol0");
    fMultV0->Fit(fpol0,"Q0","",0,31);
    fV0Cpol = fpol0->GetParameter(0);
    fMultV0->Fit(fpol0,"Q0","",32,64);
    fV0Apol = fpol0->GetParameter(0);

    for(Int_t iside=0;iside<2;iside++){
	for(Int_t icoord=0;icoord<2;icoord++){
	    for(Int_t i=0;i  < kNCenBins;i++){
		char namecont[100];
  		if(iside==0 && icoord==0)
		    snprintf(namecont,100,"hQxc2_%i",i);
		else if(iside==1 && icoord==0)
		    snprintf(namecont,100,"hQxa2_%i",i);
		else if(iside==0 && icoord==1)
		    snprintf(namecont,100,"hQyc2_%i",i);
		else if(iside==1 && icoord==1)
		    snprintf(namecont,100,"hQya2_%i",i);

		cont = (AliOADBContainer*) foadb->Get(namecont);
		if(!cont){
		    AliError(Form("OADB object %s is not available in the file %s", namecont, oadbfilename.Data()));
		    AliError("V0 Calibration not fully set!\n");
		    return;
		}

		if(!(cont->GetObject(runNumber))){
		    AliError(Form("OADB object %s is not available for run %i, trying run 137366",namecont,runNumber));
		    runNumber = 137366;
		}
		if(!(cont->GetObject(runNumber))){
		  AliError(Form("OADB object %s is not available for run %i",namecont,runNumber));
		  AliError("V0 Calibration not fully set!\n");
		  return;
		}
		fMeanQ[i][iside][icoord] = ((TH1F *) cont->GetObject(runNumber))->GetMean();
		fWidthQ[i][iside][icoord] = ((TH1F *) cont->GetObject(runNumber))->GetRMS();

		//for v3
// 		if(iside==0 && icoord==0)
// 		    snprintf(namecont,100,"hQxc3_%i",i);
// 		else if(iside==1 && icoord==0)
// 		    snprintf(namecont,100,"hQxa3_%i",i);
// 		else if(iside==0 && icoord==1)
// 		    snprintf(namecont,100,"hQyc3_%i",i);
// 		else if(iside==1 && icoord==1)
// 		    snprintf(namecont,100,"hQya3_%i",i);
// 
// 		cont = (AliOADBContainer*) foadb->Get(namecont);
// 		if(!cont){
// 		    AliError(Form("OADB object %s is not available in the file",namecont));
// 		    AliError("V0 Calibration not fully set!\n");
// 		    return;
// 		}
// 
// 		if(!(cont->GetObject(runNumber))){
// 		    AliError(Form("OADB object %s is not available for run %i, trying run 137366",namecont,runNumber));
// 		    runNumber = 137366;
// 		}
// 		if(!(cont->GetObject(runNumber))){
// 		  AliError(Form("OADB object %s is not available for run %i",namecont,runNumber));
// 		  AliError("V0 Calibration not fully set!\n");
// 		  return;
// 		}
//		fMeanQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetMean();
//		fWidthQv3[i][iside][icoord] = ((TH1F *) cont->GetObject(run))->GetRMS();

     	    }
	}
    }

    delete fpol0; fpol0=0;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetESDTrackCuts()
{
  if( fEventESD ) {
    // Create ESD track cut
    fESDtrackCuts = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts() ;
    //fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
    fESDtrackCuts->SetRequireTPCRefit(kTRUE);
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetGeometry()
{
  // Initialize the PHOS geometry 
  if(!fPHOSGeo){
    
    AliOADBContainer geomContainer("phosGeo");
    geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSGeometry.root","PHOSRotationMatrixes");
    TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(fRunNumber,"PHOSRotationMatrixes");
    
    if (fRunNumber < 224994)
      fPHOSGeo =  AliPHOSGeometry::GetInstance("IHEP") ;
    else
      fPHOSGeo =  AliPHOSGeometry::GetInstance("Run2") ;
    
    for(Int_t mod=0; mod<5; mod++) {
      if(!matrixes->At(mod)) {
	if( fDebug )
	  AliInfo(Form("No PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
	continue;
      }
      else {
	fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
	if( fDebug >1 )
	  AliInfo(Form("Adding PHOS Matrix for mod:%d, geo=%p\n", mod, fPHOSGeo));
      }
    }
  } 
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetPHOSCalibData()
{
  if( fPHOSCalibData )
    delete fPHOSCalibData; 
  fPHOSCalibData = 0;
  
  // Calibration only needed for ESD
  if( fEventESD /*&& */ ) {
    if( kLHC10h == fPeriod && fEventESD ) {
      //We have to apply re-calibration for pass1 LCH10h
      // Initialize decalibration factors in the form of the OCDB object
      AliCDBManager * man = AliCDBManager::Instance();
      man->SetRun(140000) ; //TODO; revise, this should probably not b done.
      man->SetDefaultStorage("local://OCDB");
    }
    fPHOSCalibData = new AliPHOSCalibData();
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::RejectTriggerMaskSelection()
{
  const Bool_t REJECT = true;
  const Bool_t ACCEPT = false;

  // No need to check trigger mask if no selection is done
  if( kNoSelection == fInternalTriggerSelection )
    return ACCEPT;
  
  Bool_t reject = REJECT;
  
  Bool_t isMB = (fEvent->GetTriggerMask() & (ULong64_t(1)<<1));
  Bool_t isCentral = (fEvent->GetTriggerMask() & (ULong64_t(1)<<4));
  Bool_t isSemiCentral = (fEvent->GetTriggerMask() & (ULong64_t(1)<<7));

  if( kCentralInclusive == fInternalTriggerSelection
    && isCentral ) reject = ACCEPT; // accept event.
  else if( kCentralExclusive == fInternalTriggerSelection
    && isCentral && !isSemiCentral && !isMB ) reject = ACCEPT; // accept event.

  else if( kSemiCentralInclusive == fInternalTriggerSelection
    && isSemiCentral ) reject = ACCEPT; // accept event
  else if( kSemiCentralExclusive == fInternalTriggerSelection
    && isSemiCentral && !isCentral && !isMB ) reject = ACCEPT; // accept event.

  else if( kMBInclusive == fInternalTriggerSelection
    && isMB ) reject = ACCEPT; // accept event.
  else if( kMBExclusive == fInternalTriggerSelection
    && isMB && !isCentral && !isSemiCentral ) reject = ACCEPT; // accept event.

  if( REJECT == reject )
    return REJECT;
  else {
    LogSelection(kInternalTriggerMaskSelection, fInternalRunNumber);
    return ACCEPT;
  }
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::RejectFiredTriggerClassSelection()
{
    
  TString trigClasses = fInputEvent->GetFiredTriggerClasses();
  if(fDebug) std::cout<<" Fired trigger classes: "<<trigClasses;
  
  if(!fTrigName)
    return false;
  
  if ( trigClasses.Contains(fTrigName) ) {
    if(fDebug) std::cout << "     ==> Selected for " << fTrigName << std::endl;
    return false; // selected!
  }
  else {
    if(fDebug) std::cout << std::endl;
    return true; // rejected!
  }
  
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetVertex()
{
  const AliVVertex *primaryVertex = fEvent->GetPrimaryVertex();
  if( primaryVertex ) {
    fVertex[0] = primaryVertex->GetX();
    fVertex[1] = primaryVertex->GetY();
    fVertex[2] = primaryVertex->GetZ();
  }
  else {
    AliError("Event has 0x0 Primary Vertex, defaulting to origo");
    fVertex[0] = 0;
    fVertex[1] = 0;
    fVertex[2] = 0;
  }
  fVertexVector = TVector3(fVertex);
  FillHistogram("hZvertex", fVertexVector.z(), fInternalRunNumber-0.5);
  
  if( fDebug >= 2 )
    AliInfo(Form("Vertex is set to (%.1f,%.1f,%.1f)", fVertex[0], fVertex[1], fVertex[2]));

  fVtxBin=0 ;// No support for vtx binning implemented.
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::RejectEventVertex()
{
  if( ! fEvent->GetPrimaryVertex() )
    return true; // reject
  LogSelection(kHasVertex, fInternalRunNumber);

  if ( TMath::Abs(fVertexVector.z()) > fMaxAbsVertexZ )
    return true; // reject
  LogSelection(kHasAbsVertex, fInternalRunNumber);

  if( kLHC13 == fPeriod ) {//pPb vertex and pileup cut
    const bool vertexSelected = GetAnalysisUtils()->IsVertexSelected2013pA(fEvent);
    if(! vertexSelected ) return true;//reject
    const bool pileupSelected = GetAnalysisUtils()->IsPileUpEvent(fEvent);
    if( pileupSelected ) return true;//reject   
  }

  return false; // accept event.
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetCentrality()
{
  AliCentrality *centrality = fEvent->GetCentrality();
  if( centrality )
    fCentrality=centrality->GetCentralityPercentile(fCentralityEstimator);
  else {
    AliError("Event has 0x0 centrality");
    fCentrality = -1.;
  }
  FillHistogram("hCentrality",fCentrality,fInternalRunNumber-0.5) ;

  fCentBin = GetCentralityBin(fCentrality);

  if ( fDebug >= 2 )
    AliInfo(Form("Centrality (bin) is: %f (%d)", fCentrality, fCentBin));
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::SetCentralityRun2()
{
    // Centrality calculation of Run2.
    // See https://twiki.cern.ch/twiki/bin/viewauth/ALICE/CentralityCodeSnippets
    
    fCentrality = 300;
    AliMultSelection *MultSelection = 0x0;
    
    MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
    
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        fCentrality = MultSelection->GetMultiplicityPercentile(fCentralityEstimator);
    }
    
    FillHistogram("hCentrality",fCentrality,fInternalRunNumber-0.5) ;
    fCentBin = GetCentralityBin(fCentrality);
    
    if ( fDebug >= 2 )
        AliInfo(Form("Centrality (bin) is: %f (%d)", fCentrality, fCentBin));
    
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::RejectCentrality()
{
    
    if (fRunNumber < 224994) { // Run1
        if( ! fEvent->GetCentrality() )
            return true; // reject
    }
    
    LogSelection(kHasCentrality, fInternalRunNumber);
    
    //   if( fCentrality <= 0. || fCentrality>80. )
    //     return true; // reject
    
    int lastBinUpperIndex = fCentEdges.GetSize() -1;
    if( fCentrality > fCentEdges[lastBinUpperIndex] ) {
        if( fDebug )
            AliInfo("Rejecting due to centrality outside of binning.");
        return true; // reject
    }
    LogSelection(kCentUnderUpperBinUpperEdge, fInternalRunNumber);
    
    if( fCentrality < fCentEdges[0] ) {
        if( fDebug )
            AliInfo("Rejecting due to centrality outside of binning.");
        return true; // reject
    }
    LogSelection(kCentOverLowerBinLowerEdge, fInternalRunNumber);
    
    return false;
}


//_____________________________________________________________________________
void AliAnalysisTaskPi0Flow::EvalReactionPlane()
{
  // assigns: fHaveTPCRP and fRP
  // also does a few histogram fills

  AliEventplane *eventPlane = fEvent->GetEventplane();
  if( ! eventPlane ) { AliError("Event has no event plane"); return; }
  
  Double_t reactionPlaneQ = eventPlane->GetEventplane("Q");
  FillHistogram("phiRP",reactionPlaneQ,fCentrality) ;

  if(reactionPlaneQ==999 || reactionPlaneQ < 0.){ //reaction plain was not defined
    if( fDebug ) AliInfo(Form("No Q Reaction Plane, value is %f", reactionPlaneQ));
    fHaveTPCRP = kFALSE;
  }
  else{
    if( fDebug >= 2 ) AliInfo(Form("Q Reaction Plane is %f", reactionPlaneQ));
    fHaveTPCRP = kTRUE;
  }

  if(fHaveTPCRP){
    fRP = ApplyFlattening(reactionPlaneQ, fCentrality) ;

    while(fRP<0)  fRP+=TMath::Pi();
    while(fRP>TMath::Pi())  fRP-=TMath::Pi();
    FillHistogram("phiRPflat",fRP,fCentrality) ;
    Double_t dPsi = eventPlane->GetQsubRes() ;
    FillHistogram("cos2AC",TMath::Cos(2.*dPsi),fCentrality) ;
  }
  else
    fRP=0.;
}


//____________________________________________________________________________
void  AliAnalysisTaskPi0Flow::EvalV0ReactionPlane(){
  // set: fRPV0A and fRPV0C

  // Do Manual V0 EP Calculation
  if ( fManualV0EPCalc ) 
    {
      //VZERO data
      AliVVZERO* v0 = fEvent->GetVZEROData();

      //reset Q vector info
      Double_t Qxa2 = 0, Qya2 = 0;
      Double_t Qxc2 = 0, Qyc2 = 0;

      for (Int_t iv0 = 0; iv0 < 64; iv0++) {
	Double_t phiV0 = TMath::PiOver4()*(0.5 + iv0 % 8);
	Float_t multv0 = v0->GetMultiplicity(iv0);
	if (iv0 < 32){ // V0C
	  Qxc2 += TMath::Cos(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	  Qyc2 += TMath::Sin(2*phiV0) * multv0*fV0Cpol/fMultV0->GetBinContent(iv0+1);
	} else {       // V0A
	  Qxa2 += TMath::Cos(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	  Qya2 += TMath::Sin(2*phiV0) * multv0*fV0Apol/fMultV0->GetBinContent(iv0+1);
	}
      }

      Int_t iC = -1;
      // centrality bins
      if(fCentrality < 5) iC = 0;
      else if(fCentrality < 10) iC = 1;
      else if(fCentrality < 20) iC = 2;
      else if(fCentrality < 30) iC = 3;
      else if(fCentrality < 40) iC = 4;
      else if(fCentrality < 50) iC = 5;
      else if(fCentrality < 60) iC = 6;
      else if(fCentrality < 70) iC = 7;
      else iC = 8;

      //grab for each centrality the proper histo with the Qx and Qy to do the recentering
      Double_t Qxamean2 = fMeanQ[iC][1][0];
      Double_t Qxarms2  = fWidthQ[iC][1][0];
      Double_t Qyamean2 = fMeanQ[iC][1][1];
      Double_t Qyarms2  = fWidthQ[iC][1][1];

      Double_t Qxcmean2 = fMeanQ[iC][0][0];
      Double_t Qxcrms2  = fWidthQ[iC][0][0];
      Double_t Qycmean2 = fMeanQ[iC][0][1];
      Double_t Qycrms2  = fWidthQ[iC][0][1];

      Double_t QxaCor2 = (Qxa2 - Qxamean2)/Qxarms2;
      Double_t QyaCor2 = (Qya2 - Qyamean2)/Qyarms2;
      Double_t QxcCor2 = (Qxc2 - Qxcmean2)/Qxcrms2;
      Double_t QycCor2 = (Qyc2 - Qycmean2)/Qycrms2;

      fRPV0A = TMath::ATan2(QyaCor2, QxaCor2)/2.;
      fRPV0C = TMath::ATan2(QycCor2, QxcCor2)/2.;
    }
  else // Use Official V0 EP Calculation. 
    {
      AliEventplane *eventPlane = fEvent->GetEventplane();
      if( ! eventPlane ) { AliError("Event has no event plane"); return; }
      fRPV0A = eventPlane->GetEventplane("V0A", fEvent);
      fRPV0C = eventPlane->GetEventplane("V0C", fEvent);
    }
  
  // Check that the A&C RP are within allowed range.
  if( fDebug >= 3 && (fRPV0A<0 || fRPV0A>TMath::Pi() ) )
    AliInfo(Form("RPV0A outside of permited range [0,pi]: %f, correcting", fRPV0A));
  if( fDebug >= 3 && (fRPV0C<0 || fRPV0C>TMath::Pi() ) )
    AliInfo(Form("RPV0C outside of permited range [0,pi]: %f, correcting", fRPV0C));
  while (fRPV0A<0          ) fRPV0A+=TMath::Pi() ;
  while (fRPV0A>TMath::Pi()) fRPV0A-=TMath::Pi() ;
  while (fRPV0C<0          ) fRPV0C+=TMath::Pi() ;
  while (fRPV0C>TMath::Pi()) fRPV0C-=TMath::Pi() ;

  // Reaction plane histograms before flattening
  if( fDebug >= 2 )
    AliInfo(Form("V0 Reaction Plane before flattening: A side: %f, C side: %f", fRPV0A, fRPV0C));

  FillHistogram("phiRPV0A" ,fRPV0A,fCentrality);
  FillHistogram("phiRPV0C" ,fRPV0C,fCentrality);
  FillHistogram("phiRPV0AC",fRPV0A,fRPV0C,fCentrality) ;

  // Flattening
  fRPV0A=ApplyFlatteningV0A(fRPV0A,fCentrality) ;
  while (fRPV0A<0          ) fRPV0A+=TMath::Pi() ;
  while (fRPV0A>TMath::Pi()) fRPV0A-=TMath::Pi() ;

  fRPV0C=ApplyFlatteningV0C(fRPV0C,fCentrality) ;
  while (fRPV0C<0          ) fRPV0C+=TMath::Pi() ;
  while (fRPV0C>TMath::Pi()) fRPV0C-=TMath::Pi() ;
  
  if( fDebug >= 2 )
    AliInfo(Form("V0 Reaction Plane after  flattening: A side: %f, C side: %f", fRPV0A, fRPV0C));

  FillHistogram("phiRPV0Aflat",fRPV0A,fCentrality) ;
  FillHistogram("cos2V0AC",TMath::Cos(2.*(fRPV0A-fRPV0C)),fCentrality) ;
  if(fHaveTPCRP){
    FillHistogram("phiRPV0ATPC",fRP,fRPV0A,fCentrality) ;
    FillHistogram("cos2V0ATPC",TMath::Cos(2.*(fRP-fRPV0A)),fCentrality) ;
  }

  FillHistogram("phiRPV0Cflat",fRPV0C,fCentrality) ;
  if(fHaveTPCRP){
    FillHistogram("phiRPV0CTPC",fRP,fRPV0C,fCentrality) ;
    FillHistogram("cos2V0CTPC",TMath::Cos(2.*(fRP-fRPV0C)),fCentrality) ;
  }
}
//____________________________________________________________________________
void  AliAnalysisTaskPi0Flow::EvalCoreLambdas(AliVCluster * clu, AliVCaloCells * cells,Double_t &m02, Double_t &m20){ 
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum
    
  const Double_t rCut=4.5 ;  
    
  const Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  const Int_t mulDigit=clu->GetNCells() ;
  Double_t xc[mulDigit] ;
  Double_t zc[mulDigit] ;
  Double_t wi[mulDigit] ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  const Double_t logWeight=4.5 ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi=0. ;
    Float_t zi=0. ;
    Int_t absId = clu->GetCellAbsId(iDigit) ;
    fPHOSGeo->AbsToRelNumbering(absId, relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    Double_t ei = elist[iDigit]*cells->GetCellAmplitude(absId) ;
    wi[iDigit]=0. ;
    if (clu->E()>0 && ei>0) {
      wi[iDigit] = TMath::Max( 0., logWeight + TMath::Log( ei / clu->E() ) ) ;
      Double_t w=wi[iDigit];
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
     
  wtot = 0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t xCut = 0. ;
  Double_t zCut = 0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Double_t w=wi[iDigit];
    if (w>0.) {
        Double_t xi= xc[iDigit] ;
        Double_t zi= zc[iDigit] ;
	if((xi-x)*(xi-x)+(zi-z)*(zi-z) < rCut*rCut){
          xCut += w * xi ;
          zCut += w * zi ; 
          dxx  += w * xi * xi ;
          dzz  += w * zi * zi ;
          dxz  += w * xi * zi ; 
          wtot += w ;
	}
    }
    
  }
  if (wtot>0) {
    xCut/= wtot ;
    zCut/= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xCut * xCut ;
    dzz -= zCut * zCut ;
    dxz -= xCut * zCut ;

    m02 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    m20 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    m20=m02=0.;
  }

}
//____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Flow::TestCoreLambda(Double_t pt,Double_t l1,Double_t l2){
  //Evaluates if lambdas correspond to photon cluster
  //Tuned using pp date 
  //For core radius R=4.5
  Double_t   l1Mean  = 1.150200 + 0.097886/(1.+1.486645*pt+0.000038*pt*pt) ;
  Double_t   l2Mean = 1.574706 + 0.997966*exp(-0.895075*pt)-0.010666*pt ;
  Double_t   l1Sigma = 0.100255 + 0.337177*exp(-0.517684*pt)+0.001170*pt ;
  Double_t   l2Sigma = 0.232580 + 0.573401*exp(-0.735903*pt)-0.002325*pt ;
  Double_t   c = -0.110983 -0.017353/(1.-1.836995*pt+0.934517*pt*pt) ;

  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;  
}


AliAnalysisUtils* AliAnalysisTaskPi0Flow::GetAnalysisUtils()
{
  static AliAnalysisUtils* utils = 0x0;
  if(utils) 
    return utils;

  utils = new AliAnalysisUtils();

  return utils;
}
