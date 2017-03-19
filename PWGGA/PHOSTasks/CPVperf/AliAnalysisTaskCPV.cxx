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
/* $Id$ */
 
// Analysis task for CPV performance via correlation studies
// between CPV clusters and projection of global tracks to CPV
// Author: Sergey Evdokimov

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "THashList.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskCPV.h"
#include "AliPHOSGeometry.h"
#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliVCaloCells.h"
#include "AliVCluster.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliPIDCombined.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliPHOSCPVGeometry.h"
#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliVParticle.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

using namespace std;

ClassImp(AliAnalysisTaskCPV)

//________________________________________________________________________
AliAnalysisTaskCPV::AliAnalysisTaskCPV(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer1(0),
  fOutputContainer2(0),
  fPHOSGeo(0),
  fEventCounter(0),
  fITSTPCTrackCuts(0), //Track cuts
  fPIDResponse(0x0),
  fPIDCombined(0x0),
  fGeometryCPV(0x0),
  frCPV(428.3)
{
  // Constructor
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());
  DefineOutput(2,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("Run2") ;
  fGeometryCPV  = new AliPHOSCPVGeometry() ;

}

//________________________________________________________________________
void AliAnalysisTaskCPV::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // ESD histograms
  if(fOutputContainer1 != NULL){
    delete fOutputContainer1;
  }
  fOutputContainer1 = new THashList();
  fOutputContainer1->SetOwner(kTRUE);
  
  if(fOutputContainer2 != NULL){
    delete fOutputContainer2;
  }
  fOutputContainer2 = new THashList();
  fOutputContainer2->SetOwner(kTRUE);


  fOutputContainer1->Add(new TH1I("hPHSCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer1->Add(new TH1I("hPHSCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  fOutputContainer1->Add(new TH1I("hPHSCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
  fOutputContainer1->Add(new TH1I("hPHSCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));
  fOutputContainer1->Add(new TH1I("hPHSCellMultEventM4","PHOS cell multiplicity per event, M4",2000,0,2000));
  fOutputContainer1->Add(new TH1I("hCPVCellMultEvent"  ,"CPV cell multiplicity per event"     ,2000,0,2000));

  fOutputContainer1->Add(new TH1I("hPHSClusterMult"  ,"PHOS cluster multiplicity"    ,200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSClusterMultM1","PHOS cluster multiplicity, M1",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSClusterMultM2","PHOS cluster multiplicity, M2",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSClusterMultM3","PHOS cluster multiplicity, M3",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSClusterMultM4","PHOS cluster multiplicity, M4",200,0,200));
  fOutputContainer1->Add(new TH1I("hCPVClusterMult"  ,"CPV cluster multiplicity"     ,400,0,400));

  fOutputContainer1->Add(new TH1F("hPHSCellEnergy"  ,"PHOS cell energy"            ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSCellEnergyM1","PHOS cell energy in module 1",500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSCellEnergyM2","PHOS cell energy in module 2",500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSCellEnergyM3","PHOS cell energy in module 3",500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSCellEnergyM4","PHOS cell energy in module 4",500,0.,50.));
  fOutputContainer1->Add(new TH1F("hCPVCellEnergy"  ,"CPV cell energy"             ,500,0.,5000.));

  fOutputContainer1->Add(new TH1F("hPHSClusterEnergy"  ,"PHOS cluster energy"      ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSClusterEnergyM1","PHOS cluster energy, M1"  ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSClusterEnergyM2","PHOS cluster energy, M2"  ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSClusterEnergyM3","PHOS cluster energy, M3"  ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hPHSClusterEnergyM4","PHOS cluster energy, M4"  ,500,0.,50.));
  fOutputContainer1->Add(new TH1F("hCPVClusterEnergy"  ,"CPV cluster energy"       ,1000,0.,20000.));

  fOutputContainer1->Add(new TH1F("hPHSTrackDx"  ,"PHOS track Dx"    ,200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDxM1","PHOS track Dx, M1",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDxM2","PHOS track Dx, M2",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDxM3","PHOS track Dx, M3",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDxM4","PHOS track Dx, M4",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hCPVTrackDx"  ,"CPV track Dx"     ,200,-50.,50.));
  
  fOutputContainer1->Add(new TH1F("hPHSTrackDz"  ,"PHOS track Dz"    ,200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDzM1","PHOS track Dz, M1",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDzM2","PHOS track Dz, M2",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDzM3","PHOS track Dz, M3",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hPHSTrackDzM4","PHOS track Dz, M4",200,-50.,50.));
  fOutputContainer1->Add(new TH1F("hCPVTrackDz"  ,"CPV track Dz"     ,200,-50.,50.));
  
  fOutputContainer1->Add(new TH2F("hPHSClusterEvsN"  ,"PHOS cluster energy vs digit multiplicity"    ,500,0.,50.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hPHSClusterEvsNM1","PHOS cluster energy vs digit multiplicity, M1",500,0.,50.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hPHSClusterEvsNM2","PHOS cluster energy vs digit multiplicity, M2",500,0.,50.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hPHSClusterEvsNM3","PHOS cluster energy vs digit multiplicity, M3",500,0.,50.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hPHSClusterEvsNM4","PHOS cluster energy vs digit multiplicity, M4",500,0.,50.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCPVClusterEvsN"  ,"CPV cluster energy vs digit multiplicity"     ,200,0.,10000.,40,0.,40.));

  fOutputContainer1->Add(new TH1I("hPHSCellMultClu"  ,"PHOS cell multiplicity per cluster"    ,200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSCellMultCluM1","PHOS cell multiplicity per cluster, M1",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSCellMultCluM2","PHOS cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSCellMultCluM3","PHOS cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer1->Add(new TH1I("hPHSCellMultCluM4","PHOS cell multiplicity per cluster, M4",200,0,200));
  fOutputContainer1->Add(new TH1I("hCPVCellMultClu"  ,"CPV cell multiplicity per cluster"     ,200,0,200));

  fOutputContainer1->Add(new TH2F("hPHSCellNXZM1","PHOS cell (X,Z), M1" ,64,0.,64., 56,0.,56.));
  fOutputContainer1->Add(new TH2F("hPHSCellNXZM2","PHOS cell (X,Z), M2" ,64,0.,64., 56,0.,56.));
  fOutputContainer1->Add(new TH2F("hPHSCellNXZM3","PHOS cell (X,Z), M3" ,64,0.,64., 56,0.,56.));
  fOutputContainer1->Add(new TH2F("hPHSCellNXZM4","PHOS cell (X,Z), M4" ,64,0.,64., 56,0.,56.));
  fOutputContainer1->Add(new TH2F("hCPVCellNXZ"  ,"CPV cell (X,Z)"      ,128,0.,128., 60,0.,60.));
  fOutputContainer1->Add(new TH2F("hCPVCellNXZweighted"  ,"CPV cell (X,Z) weighted"      ,128,0.,128., 60,0.,60.));

  fOutputContainer1->Add(new TH2F("hPHSClusNXZM1","PHOS cluster (X,Z), M1" ,75,-75.,75., 70,-70.,70.));
  fOutputContainer1->Add(new TH2F("hPHSClusNXZM2","PHOS cluster (X,Z), M2" ,75,-75.,75., 70,-70.,70.));
  fOutputContainer1->Add(new TH2F("hPHSClusNXZM3","PHOS cluster (X,Z), M3" ,75,-75.,75., 70,-70.,70.));
  fOutputContainer1->Add(new TH2F("hPHSClusNXZM4","PHOS cluster (X,Z), M4" ,75,-75.,75., 70,-70.,70.));
  fOutputContainer1->Add(new TH2F("hCPVClusNXZ"  ,"CPV cluster (X,Z)"      ,75,-75.,75., 70,-70.,70.));

  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrR"  ,"distance R bwt CPV and PHOS nearest clust M3"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrX"  ,"distance X bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrZ"  ,"distance Z bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM3corrXZ"  ,"distance XvsZ bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.,1500,-75.,75.));


  fOutputContainer1->Add(new TH1F("hCPVPHOSM1corrR"  ,"distance R bwt CPV and PHOS nearest clust M1"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM1corrX"  ,"distance X bwt CPV and PHOS nearest clust M1"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM1corrZ"  ,"distance Z bwt CPV and PHOS nearest clust M1"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM1corrXZ"  ,"distance XvsZ bwt CPV and PHOS nearest clust M1"      ,1500,-75.,75.,1500,-75.,75.));


  fOutputContainer1->Add(new TH1F("hCPVPHOSM2corrR"  ,"distance R bwt CPV and PHOS nearest clust M2"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM2corrX"  ,"distance X bwt CPV and PHOS nearest clust M2"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM2corrZ"  ,"distance Z bwt CPV and PHOS nearest clust M2"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM2corrXZ"  ,"distance XvsZ bwt CPV and PHOS nearest clust M2"      ,1500,-75.,75.,1500,-75.,75.));


  fOutputContainer1->Add(new TH1F("hCPVPHOSM4corrR"  ,"distance R bwt CPV and PHOS nearest clust M4"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM4corrX"  ,"distance X bwt CPV and PHOS nearest clust M4"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM4corrZ"  ,"distance Z bwt CPV and PHOS nearest clust M4"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM4corrXZ"  ,"distance XvsZ bwt CPV and PHOS nearest clust M4"      ,1500,-75.,75.,1500,-75.,75.));


  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrREgt03"  ,"distance R bwt CPV and PHOS nearest clust M3"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrXEgt03"  ,"distance X bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrZEgt03"  ,"distance Z bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM3corrXZEgt03"  ,"distance XvsZ bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.,1500,-75.,75.));

  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrRMIP"  ,"distance R bwt CPV and PHOS nearest clust M3"      ,750,0.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrXMIP"  ,"distance X bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCPVPHOSM3corrZMIP"  ,"distance Z bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.));
  fOutputContainer1->Add(new TH2F("hCPVPHOSM3corrXZMIP"  ,"distance XvsZ bwt CPV and PHOS nearest clust M3"      ,1500,-75.,75.,1500,-75.,75.));


  fOutputContainer1->Add(new TH2F("hPHSClusNXZM3Egt03","PHOS cluster (X,Z), M3" ,75,-75.,75., 70,-70.,70.));
  fOutputContainer1->Add(new TH2F("hPHSClusNXZM3MIP","PHOS cluster (X,Z), M3" ,75,-75.,75., 70,-70.,70.));

  for (Int_t i = 0;i<5;i++){
    fOutputContainer1->Add(new TH1F(Form("hCPVPHOSM3corrRE%d-%d",i,i+1)  ,Form("distance R bwt CPV and PHOS nearest clust E(%d,%d)  M3",i,i+1)      ,750,0.,75.));
    fOutputContainer1->Add(new TH1F(Form("hCPVPHOSM3corrXE%d-%d",i,i+1)  ,Form("distance X bwt CPV and PHOS nearest clust E(%d,%d) M3",i,i+1)      ,1500,-75.,75.));
    fOutputContainer1->Add(new TH1F(Form("hCPVPHOSM3corrZE%d-%d",i,i+1)  ,Form("distance Z bwt CPV and PHOS nearest clust E(%d,%d) M3",i,i+1)      ,1500,-75.,75.));
    fOutputContainer1->Add(new TH2F(Form("hCPVPHOSM3corrXZE%d-%d",i,i+1)  ,Form("distance XvsZ bwt CPV and PHOS nearest clust E(%d,%d) M3",i,i+1)      ,1500,-75.,75.,1500,-75.,75.));
    fOutputContainer1->Add(new TH2F(Form("hPHSClusNXZM3E%d-%d",i,i+1),Form("PHOS cluster (X,Z), E(%d,%d) M3",i,i+1) ,75,-75.,75., 70,-70.,70.));
  }

  //ITS-TPC tracks matching 
  fOutputContainer1->Add(new TH1F("hITSTPCtrackN"  ,"Number of ITSTPC tracks"  ,100,0.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackNinAcc"  ,"Number of ITSTPC tracks in CPV acceptance"  ,100,0.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackPinAcc"  ,"Momentum of ITSTPC tracks in CPV acceptance"  ,100,0.,10.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackXinAcc"  ,"X coord of ITSTPC tracks on CPV plane"  ,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackZinAcc"  ,"Z coord of ITSTPC tracks on CPV plane"  ,2000,-100.,100.));
  fOutputContainer1->Add(new TH2F("hITSTPCtrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCtrackInAccXZ"  ,"ITSTPC tracks projection on CPV module"      ,750,-75.,75.,750,-75.,75.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackProbEinAcc"  ,"Probability of ITSTPC track to be electron"  ,100,0.,1.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackProbMuinAcc"  ,"Probability of ITSTPC track to be muon"  ,100,0.,1.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackProbPiinAcc"  ,"Probability of ITSTPC track to be pion"  ,100,0.,1.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackProbKinAcc"  ,"Probability of ITSTPC track to be kaon"  ,100,0.,1.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackProbPinAcc"  ,"Probability of ITSTPC track to be proton"  ,100,0.,1.));





  fOutputContainer1->Add(new TH1F("hITSTPCtrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCtrackMatchedInAccXZ"  ,"Matched ITSTPC tracks projection on CPV module" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackMatchedInAccP"  ,"Matched ITSTPC tracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackMatchedInAccPID"  ,"Matched ITSTPC tracks PID" ,10,0.,10.));

  fOutputContainer1->Add(new TH2F("hITSTPCtrackNonMatchedInAccXZ"  ,"Not Matched ITSTPC tracks projection on CPV module" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackNonMatchedInAccP"  ,"Not Matched ITSTPC tracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH1F("hITSTPCtrackNonMatchedInAccPID"  ,"Not Matched ITSTPC tracks PID" ,10,0.,10.));

  fOutputContainer1->Add(new TH2F("hCpvCluNonMatchedXZ"  ,"Not Matched CPV clusters" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCpvCluNonMatchedEn"  ,"Not Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluNonMatchedN"  ,"Not Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluNonMatchedEvsN"  ,"Not Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedXZ"  ,"Matched CPV clusters" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));


  //electron histograms
  fOutputContainer1->Add(new TH2F("hITSTPCelectrondXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest electron" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCelectronInAccXZ"  ,"ITSTPC electrons projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCelectronDist"  ,"distance b/w CPV cluster and nearest ITSTPC electron in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCelectronMatchedInAccXZ"  ,"Matched ITSTPC electrons projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCelectronMatchedInAccP"  ,"Matched ITSTPC electrons momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCelectronNonMatchedInAccXZ"  ,"Matched ITSTPC electrons projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCelectronNonMatchedInAccP"  ,"Matched ITSTPC electrons momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedElectronXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedElectronEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedElectronN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedElectronEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  //high momentum histograms (P>0.5GeV/c)
  fOutputContainer1->Add(new TH2F("hITSTPCHMtrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest HMtrack" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCHMtrackInAccXZ"  ,"ITSTPC HMtracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCHMtrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC HMtrack in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCHMtrackMatchedInAccXZ"  ,"Matched ITSTPC HMtracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCHMtrackMatchedInAccP"  ,"Matched ITSTPC HMtracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCHMtrackNonMatchedInAccXZ"  ,"Matched ITSTPC HMtracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCHMtrackNonMatchedInAccP"  ,"Matched ITSTPC HMtracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedHMtrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedHMtrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedHMtrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedHMtrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  
  //positive tracks
  fOutputContainer1->Add(new TH2F("hITSTPCPosTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest PosTrack" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosTrackInAccXZ"  ,"ITSTPC PosTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC PosTrack in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosTrackMatchedInAccXZ"  ,"Matched ITSTPC PosTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosTrackMatchedInAccP"  ,"Matched ITSTPC PosTracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosTrackNonMatchedInAccXZ"  ,"Matched ITSTPC PosTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosTrackNonMatchedInAccP"  ,"Matched ITSTPC PosTracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));
  
  //negative tracks
  fOutputContainer1->Add(new TH2F("hITSTPCNegTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest NegTrack" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegTrackInAccXZ"  ,"ITSTPC NegTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC NegTrack in CPV acceptance"  ,250,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegTrackMatchedInAccXZ"  ,"Matched ITSTPC NegTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegTrackMatchedInAccP"  ,"Matched ITSTPC NegTracks momentum" ,1000,0.,100.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegTrackNonMatchedInAccXZ"  ,"Matched ITSTPC NegTracks projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegTrackNonMatchedInAccP"  ,"Matched ITSTPC NegTracks momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  
  //Z>0 tracks
  fOutputContainer1->Add(new TH2F("hITSTPCPosZTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest Z>0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZTrackInAccXZ"  ,"ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC Z>0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZTrackMatchedInAccXZ"  ,"Matched ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZTrackMatchedInAccP"  ,"Matched ITSTPC Z>0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZTrackNonMatchedInAccXZ"  ,"Matched ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZTrackNonMatchedInAccP"  ,"Matched ITSTPC Z>0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosZTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosZTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosZTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosZTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  //Z<0 tracks
  fOutputContainer1->Add(new TH2F("hITSTPCNegZTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest Z<0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZTrackInAccXZ"  ,"ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC Z<0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZTrackMatchedInAccXZ"  ,"Matched ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZTrackMatchedInAccP"  ,"Matched ITSTPC Z<0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZTrackNonMatchedInAccXZ"  ,"Matched ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZTrackNonMatchedInAccP"  ,"Matched ITSTPC Z<0 track momentum" ,1000,0.,100.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegZTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegZTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegZTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegZTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

    //Z>0 HM tracks
  fOutputContainer1->Add(new TH2F("hITSTPCPosZHMTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest Z>0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZHMTrackInAccXZ"  ,"ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZHMTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC Z>0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZHMTrackMatchedInAccXZ"  ,"Matched ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZHMTrackMatchedInAccP"  ,"Matched ITSTPC Z>0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCPosZHMTrackNonMatchedInAccXZ"  ,"Matched ITSTPC Z>0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCPosZHMTrackNonMatchedInAccP"  ,"Matched ITSTPC Z>0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosZHMTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosZHMTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedPosZHMTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedPosZHMTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  //Z<0 HM tracks
  fOutputContainer1->Add(new TH2F("hITSTPCNegZHMTrackdXdZ"  ,"distance dXvsdZ bwt CPV  and TPC nearest Z<0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZHMTrackInAccXZ"  ,"ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZHMTrackDist"  ,"distance b/w CPV cluster and nearest ITSTPC Z<0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZHMTrackMatchedInAccXZ"  ,"Matched ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZHMTrackMatchedInAccP"  ,"Matched ITSTPC Z<0 track momentum" ,100,0.,10.));
  fOutputContainer1->Add(new TH2F("hITSTPCNegZHMTrackNonMatchedInAccXZ"  ,"Matched ITSTPC Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hITSTPCNegZHMTrackNonMatchedInAccP"  ,"Matched ITSTPC Z<0 track momentum" ,1000,0.,100.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegZHMTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegZHMTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer1->Add(new TH1F("hCpvCluMatchedNegZHMTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer1->Add(new TH2F("hCpvCluMatchedNegZHMTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

    //Z>0 Positive tracks
  fOutputContainer2->Add(new TH2F("hPosZPosTrackdXdZ"  ,"dist. dXvsdZ bwt CPV  and TPC  pos. Z>0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer2->Add(new TH2F("hPosZPosTrackInAccXZ"  ," pos. Z>0  track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZPosTrackDist"  ,"dist. b/w CPV cluster and   pos. Z>0  track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer2->Add(new TH2F("hPosZPosTrackMatchedInAccXZ"  ,"Matched  pos. Z>0  track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZPosTrackMatchedInAccP"  ,"Matched  pos. Z>0  track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hPosZPosTrackNonMatchedInAccXZ"  ,"Matched  pos. Z>0  track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZPosTrackNonMatchedInAccP"  ,"Matched  pos. Z>0  track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedPosZPosTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedPosZPosTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedPosZPosTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedPosZPosTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  //Z<0 Negative tracks
  fOutputContainer2->Add(new TH2F("hNegZNegTrackdXdZ"  ,"dist. dXvsdZ bwt CPV  and TPC  neg. Z<0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer2->Add(new TH2F("hNegZNegTrackInAccXZ"  ," neg. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZNegTrackDist"  ,"dist. b/w CPV cluster and   neg. Z<0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer2->Add(new TH2F("hNegZNegTrackMatchedInAccXZ"  ,"Matched  neg. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZNegTrackMatchedInAccP"  ,"Matched  neg. Z<0 track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hNegZNegTrackNonMatchedInAccXZ"  ,"Matched  neg. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZNegTrackNonMatchedInAccP"  ,"Matched  neg. Z<0 track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedNegZNegTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedNegZNegTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedNegZNegTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedNegZNegTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

    //Z>0 negative tracks
  fOutputContainer2->Add(new TH2F("hPosZNegTrackdXdZ"  ,"dist. dXvsdZ bwt CPV  and TPC  neg. Z>0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer2->Add(new TH2F("hPosZNegTrackInAccXZ"  ," neg. Z>0  track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZNegTrackDist"  ,"dist. b/w CPV cluster and   neg. Z>0  track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer2->Add(new TH2F("hPosZNegTrackMatchedInAccXZ"  ,"Matched  neg. Z>0  track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZNegTrackMatchedInAccP"  ,"Matched  neg. Z>0  track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hPosZNegTrackNonMatchedInAccXZ"  ,"Matched  neg. Z>0  track projection500 on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hPosZNegTrackNonMatchedInAccP"  ,"Matched  neg. Z>0  track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedPosZNegTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedPosZNegTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedPosZNegTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedPosZNegTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));

  //Z<0 Positive tracks
  fOutputContainer2->Add(new TH2F("hNegZPosTrackdXdZ"  ,"dist. dXvsdZ bwt CPV  and TPC  pos. Z<0 track" ,500,-25.,25.,500,-25.,25.));
  fOutputContainer2->Add(new TH2F("hNegZPosTrackInAccXZ"  ," pos. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZPosTrackDist"  ,"dist. b/w CPV cluster and   pos. Z<0 track in CPV acceptance"  ,250,0.,25.));
  fOutputContainer2->Add(new TH2F("hNegZPosTrackMatchedInAccXZ"  ,"Matched  pos. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZPosTrackMatchedInAccP"  ,"Matched  pos. Z<0 track momentum" ,100,0.,10.));
  fOutputContainer2->Add(new TH2F("hNegZPosTrackNonMatchedInAccXZ"  ,"Matched  pos. Z<0 track projection on CPV module" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hNegZPosTrackNonMatchedInAccP"  ,"Matched  pos. Z<0 track momentum" ,200,0.,10.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedNegZPosTrackXZ"  ,"Matched CPV clusters" ,2000,-100.,100.,2000,-100.,100.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedNegZPosTrackEn"  ,"Matched CPV clusters energie" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluMatchedNegZPosTrackN"  ,"Matched CPV clusters size" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluMatchedNegZPosTrackEvsN"  ,"Matched CPV clusters energy vs size" ,2000,0.,10000.,40,0.,40.));



  //unmatched clusters from hot zones
  fOutputContainer2->Add(new TH2F("hCpvCluHotNonMatchedXZ"  ,"Matched CPV clusters from hot zones" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer2->Add(new TH1F("hCpvCluHotNonMatchedEn"  ,"Matched CPV clusters energie from hot zones" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluHotNonMatchedN"  ,"Matched CPV clusters size from hot zones" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluHotNonMatchedEvsN"  ,"Matched CPV clusters energy vs size from hot zones" ,2000,0.,10000.,40,0.,40.));

    //unmatched clusters from cold zones
  fOutputContainer2->Add(new TH2F("hCpvCluColdNonMatchedXZ"  ,"Matched CPV clusters from cold zones" ,750,-75.,75.,750,-75.,75.));
  fOutputContainer2->Add(new TH1F("hCpvCluColdNonMatchedEn"  ,"Matched CPV clusters energie from cold zones" ,2000,0.,10000.));
  fOutputContainer2->Add(new TH1F("hCpvCluColdNonMatchedN"  ,"Matched CPV clusters size from cold zones" ,40,0.,40.));
  fOutputContainer2->Add(new TH2F("hCpvCluColdNonMatchedEvsN"  ,"Matched CPV clusters energy vs size from hot zones" ,2000,0.,10000.,40,0.,40.));

  

  PostData(1, fOutputContainer1);
  PostData(2, fOutputContainer2);

  fITSTPCTrackCuts = ITSTPCTrackCuts();
  
  fPIDCombined=new AliPIDCombined;
  fPIDCombined->SetDefaultTPCPriors();
  fPIDCombined->SetDetectorMask(AliPIDResponse::kDetITS|AliPIDResponse::kDetTPC); 

}

//________________________________________________________________________
void AliAnalysisTaskCPV::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD

  AliVEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) {
     AliError("ERROR: Could not retrieve event");
     return;
  }

  TString trigClasses     = event->GetFiredTriggerClasses();
  Int_t eventNumberInFile = event->GetEventNumberInFile();

  //PID responce
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse=inputHandler->GetPIDResponse();

  if(!fPIDResponse) AliError("No PIDResponse found !!! PID from ESD will be used");
  
  //ITSTPC tracks
  Int_t ntr = event->GetNumberOfTracks();
  if(ntr==0) {
    AliError("ERROR: no tracks found in current event!");
    PostData(1, fOutputContainer1);
    PostData(2, fOutputContainer2);
    return;
  }
  //main CPV geometry parameters
  if (frCPV<0) frCPV  = fPHOSGeo->GetIPtoCPVDistance(); //Distance to center of  CPV module
  const Double_t kXmax = fPHOSGeo->GetCPVBoxSize(0); //Size of the CPV module 
  const Double_t kZmax = fPHOSGeo->GetCPVBoxSize(2); //Size of the CPV module
  const Double_t kAlpha0=270./180.*TMath::Pi() ; //angle of the center of CPV module M=3

  // Printf("AliAnalysisTaskCPV::UserExec(): rCPV=%.1f cm\n",frCPV);

//=====================  prepare and count ITSTPCtracks in CPV acceptance  ===================
  AliESDtrack** ITSTPCtracks = new  AliESDtrack*[ntr];
  AliESDtrack** ITSTPCtracksInAcc = new  AliESDtrack*[ntr];
  // AliESDtrack** ITSTPCelectronsInAcc = new  AliESDtrack*[ntr];
  // AliESDtrack** ITSTPCHMtracksInAcc = new  AliESDtrack*[ntr];
  // AliESDtrack** ITSTPCPosTracksInAcc = new  AliESDtrack*[ntr];
  // AliESDtrack** ITSTPCNegTracksInAcc = new  AliESDtrack*[ntr];

  Float_t * trackInAccX=new Float_t[ntr];
  Float_t * trackInAccZ=new Float_t[ntr];
  // Float_t * electronInAccX=new Float_t[ntr];
  // Float_t * electronInAccZ=new Float_t[ntr];
  // Float_t * HMtrackInAccX=new Float_t[ntr];
  // Float_t * HMtrackInAccZ=new Float_t[ntr];
  // Float_t * PosTrackInAccX=new Float_t[ntr];
  // Float_t * PosTrackInAccZ=new Float_t[ntr];
  // Float_t * NegTrackInAccX=new Float_t[ntr];
  // Float_t * NegTrackInAccZ=new Float_t[ntr];

  Float_t * trackInAccProbE=new Float_t[ntr];
  Float_t * trackInAccProbMu=new Float_t[ntr];
  Float_t * trackInAccProbPi=new Float_t[ntr];
  Float_t * trackInAccProbK=new Float_t[ntr];
  Float_t * trackInAccProbP=new Float_t[ntr];

  Int_t nITSTPCtracks=0,nITSTPCtracksInAcc=0,nITSTPCelectronsInAcc=0,nITSTPCHMtracksInAcc=0,nITSTPCPosTracksInAcc=0,nITSTPCNegTracksInAcc=0,nITSTPCPosZTracksInAcc=0,nITSTPCNegZTracksInAcc=0;

  Bool_t * isHMtrack = new Bool_t[ntr];
  Bool_t * isElectron = new Bool_t[ntr];
  Bool_t * isPosTrack = new Bool_t[ntr];
  Bool_t * isNegTrack = new Bool_t[ntr];
  Bool_t * isPosZTrack = new Bool_t[ntr];
  Bool_t * isNegZTrack = new Bool_t[ntr];



  for(Int_t i=0;i<ntr;i++) {
    isHMtrack[i]=0;
    isElectron[i]=0;
    isPosTrack[i]=0;
    isNegTrack[i]=0;
    isPosZTrack[i]=0;
    isNegZTrack[i]=0;

    AliESDtrack *track = (AliESDtrack*)event->GetTrack(i);
    if(fITSTPCTrackCuts->AcceptTrack(track)) {//good ITSTPC track      
      ITSTPCtracks[nITSTPCtracks]=track;
      nITSTPCtracks++;
      //is it in CPV acceptance?
      Double_t bz = AliTracker::GetBz() ; //B-Field for approximate matching
      Double_t b[3]; 
      ULong_t status = track->GetStatus();
      if ((status & AliESDtrack::kTPCout)   == 0) continue;//track doesn't leave TPC
      const AliExternalTrackParam *outerParam = track->GetOuterParam();
      if (!outerParam) continue;
      Double_t z; 
      if(!outerParam->GetZAt(frCPV,bz,z)) continue ;
      if (TMath::Abs(z) > kZmax/2.) continue; // Some tracks miss the CPV module in Z
	
	
      AliExternalTrackParam t(*outerParam);
      if(!t.RotateParamOnly(kAlpha0)) continue ; // Rotate track to have CPV module direction angle = 0 degree.
      Double_t y;                       // Some tracks do not reach the CPV module
      if (!t.GetYAt(frCPV,bz,y)) continue; //    because of the bending
	
      if(TMath::Abs(y) < kXmax/2.){//note rotation here
	  //track is in CPV acceptance
	t.GetBxByBz(b) ;
	t.PropagateParamOnlyBxByBzTo(frCPV,b);        // Propagate to the matching module
	Double_t gposTrack[3] ; 
	t.GetXYZ(gposTrack) ;
	TVector3 globalPositionTr(gposTrack) ;
	TVector3 localPositionTr ;
	fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,3) ;
	FillHistogram("hITSTPCtrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
     	ITSTPCtracksInAcc[nITSTPCtracksInAcc]=track;
	trackInAccX[nITSTPCtracksInAcc]=localPositionTr.X();
	trackInAccZ[nITSTPCtracksInAcc]=localPositionTr.Z();
	Double_t bayesProbabilities[5]={0.};
	fPIDCombined->ComputeProbabilities(track, fPIDResponse, bayesProbabilities);
	trackInAccProbE[nITSTPCtracksInAcc]=bayesProbabilities[0];
	trackInAccProbMu[nITSTPCtracksInAcc]=bayesProbabilities[1];
	trackInAccProbPi[nITSTPCtracksInAcc]=bayesProbabilities[2];
	trackInAccProbK[nITSTPCtracksInAcc]=bayesProbabilities[3];
	trackInAccProbP[nITSTPCtracksInAcc]=bayesProbabilities[4];
	FillHistogram("hITSTPCtrackProbEinAcc",bayesProbabilities[0]);
	FillHistogram("hITSTPCtrackProbMuinAcc",bayesProbabilities[1]);
	FillHistogram("hITSTPCtrackProbPiinAcc",bayesProbabilities[2]);
	FillHistogram("hITSTPCtrackProbKinAcc",bayesProbabilities[3]);
	FillHistogram("hITSTPCtrackProbPinAcc",bayesProbabilities[4]);
        nITSTPCtracksInAcc++;
	//continue;
	if(bayesProbabilities[0]>0.5){//we have electron
	  // electronInAccX[nITSTPCelectronsInAcc]=localPositionTr.X();
	  // electronInAccZ[nITSTPCelectronsInAcc]=localPositionTr.Z();
	  // ITSTPCelectronsInAcc[nITSTPCelectronsInAcc]=track;
	  FillHistogram("hITSTPCelectronInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCelectronsInAcc++;
	  isElectron[nITSTPCtracksInAcc-1]=1;
	}
	if(track->P()>1.){//we have high momentum track
	  // HMtrackInAccX[nITSTPCHMtracksInAcc]=localPositionTr.X();
	  // HMtrackInAccZ[nITSTPCHMtracksInAcc]=localPositionTr.Z();
	  // ITSTPCHMtracksInAcc[nITSTPCHMtracksInAcc]=track;
	  FillHistogram("hITSTPCHMtrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCHMtracksInAcc++;
	  isHMtrack[nITSTPCtracksInAcc-1]=1;
	}
	if(track->Charge()>0.0){//we have positive track
	  // PosTrackInAccX[nITSTPCPosTracksInAcc]=localPositionTr.X();
	  // PosTrackInAccZ[nITSTPCPosTracksInAcc]=localPositionTr.Z();
	  // ITSTPCPosTracksInAcc[nITSTPCPosTracksInAcc]=track;
	  FillHistogram("hITSTPCPosTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCPosTracksInAcc++;
	  isPosTrack[nITSTPCtracksInAcc-1]=1;
	}
	if(track->Charge()<0.0){//we have negative track
	  // NegTrackInAccX[nITSTPCNegTracksInAcc]=localPositionTr.X();
	  // NegTrackInAccZ[nITSTPCNegTracksInAcc]=localPositionTr.Z();
	  // ITSTPCNegTracksInAcc[nITSTPCNegTracksInAcc]=track;
	  FillHistogram("hITSTPCNegTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCNegTracksInAcc++;
	  isNegTrack[nITSTPCtracksInAcc-1]=1;

	}
	if(z<0.){//Z<0 track
	  FillHistogram("hITSTPCNegZTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCNegZTracksInAcc++;
	  isNegZTrack[nITSTPCtracksInAcc-1]=1;
	}
	if(z>=0.){//Z>0 track
	  FillHistogram("hITSTPCPosZTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  nITSTPCPosZTracksInAcc++;
	  isPosZTrack[nITSTPCtracksInAcc-1]=1;
	}
	if(z<0.&&track->P()>1.){//Z<0 HM track
	  FillHistogram("hITSTPCNegZHMTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}
	if(z>=0.&&track->P()>1.){//Z>0 HM track
	  FillHistogram("hITSTPCPosZHMTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}
	if(z<0.&&track->Charge()>0.0){//Z<0 positive tracks
	  FillHistogram("hNegZPosTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}
	if(z<0.&&track->Charge()<0.0){//Z<0 negative tracks
	  FillHistogram("hNegZNegTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}
	if(z>0.&&track->Charge()>0.0){//Z>0 positive tracks
	  FillHistogram("hPosZPosTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}
	if(z>0.&&track->Charge()<0.0){//Z>0 negative tracks
	  FillHistogram("hPosZNegTrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	}


	
      }
    }
  }
  FillHistogram("hITSTPCtrackN",nITSTPCtracks);
  FillHistogram("hITSTPCtrackNinAcc",nITSTPCtracksInAcc);


//================  analyse clusters ======================= 
  //CPV clusters
  AliVCluster *clu1;
  AliVCaloCells *cells      = event->GetPHOSCells();

  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t multCells = cells->GetNumberOfCells();

  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  if(fEventCounter == 0) {
    const TGeoHMatrix *PHOSMatrix;
    for(Int_t mod=0; mod<5; mod++) {
      PHOSMatrix = ((AliESDEvent*)event)->GetPHOSMatrix(mod);
      if(!PHOSMatrix) {
	AliError(Form("PHOS geo matrix for module # %d is missing\n", mod));
	continue;
      }
      fPHOSGeo->SetMisalMatrix(PHOSMatrix,mod) ;
      AliInfo(Form("PHOS geo matrix for module # %d is set to %p\n", mod, PHOSMatrix));
    }
  }

  Float_t  energy, tof, dX, dZ;
  Int_t    mod1, detector, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells
  Int_t nCPVCells=0;
  Int_t nPHSCellModule[5] = {0,0,0,0,0};
  Bool_t isPHS=kFALSE, isCPV=kFALSE;

  for (Int_t iCell=0; iCell<multCells; iCell++) {
    cellAbsId = cells->GetCellNumber(iCell);
    if (cellAbsId > 0) {
      isPHS = kTRUE;
      isCPV = kFALSE;
    }
    if (cellAbsId < 0) {
      isPHS = kFALSE;
      isCPV = kTRUE;
      nCPVCells++;
      cellAbsId = fPHOSGeo->GetNModules() * fPHOSGeo->GetNCristalsInModule() - cellAbsId;
    }

    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1     = relId[0];
    detector = relId[1];
    cellX    = relId[2];
    cellZ    = relId[3] ;
    energy   = cells->GetAmplitude(iCell);

    if (isPHS) {
      nPHSCellModule[0]++;
      FillHistogram("hPHSCellEnergy",energy);
      if      (mod1==1) {
	nPHSCellModule[1]++;
	FillHistogram("hPHSCellEnergyM1",energy);
	FillHistogram("hPHSCellNXZM1",cellX-0.5,cellZ-0.5,1.);
      }
      else if (mod1==2) {
	nPHSCellModule[2]++;
	FillHistogram("hPHSCellEnergyM2",energy);
	FillHistogram("hPHSCellNXZM2",cellX-0.5,cellZ-0.5,1.);
      }
      else if (mod1==3) {
	nPHSCellModule[3]++;
	FillHistogram("hPHSCellEnergyM3",energy);
	FillHistogram("hPHSCellNXZM3",cellX-0.5,cellZ-0.5,1.);
      }
      else if (mod1==4) {
	nPHSCellModule[4]++;
	FillHistogram("hPHSCellEnergyM4",energy);
	FillHistogram("hPHSCellNXZM4",cellX-0.5,cellZ-0.5,1.);
      }
    }
    else if (isCPV) {
      nCPVCells++;
      FillHistogram("hCPVCellEnergy",energy);
      FillHistogram("hCPVCellNXZ",cellX-0.5,cellZ-0.5,1.);
      FillHistogram("hCPVCellNXZweighted",cellX-0.5,cellZ-0.5,energy);
    }
  }
  FillHistogram("hPHSCellMultEvent"  ,nPHSCellModule[0]);
  FillHistogram("hPHSCellMultEventM1",nPHSCellModule[1]);
  FillHistogram("hPHSCellMultEventM2",nPHSCellModule[2]);
  FillHistogram("hPHSCellMultEventM3",nPHSCellModule[3]);
  FillHistogram("hPHSCellMultEventM4",nPHSCellModule[4]);
  FillHistogram("hCPVCellMultEvent"  ,nCPVCells);
  
  // Single loop over clusters fills cluster histograms

  Int_t    module;
  Int_t    digMult;
  Int_t    multPHOSClust[5]  = {0,0,0,0,0};
  Int_t    multPHOSClustM3Egt03 =0, multPHOSClustM3MIP =0;
  Int_t    multCPVClust = 0;
  Float_t  position[3];
  
  Float_t * cpvX = new Float_t[multClust];
  Float_t * cpvZ = new Float_t[multClust];
  Float_t * cpvE = new Float_t[multClust];
  Int_t *   cpvN = new Int_t  [multClust];
  Float_t * cpvTrackDZ = new Float_t[multClust];
  Float_t * cpvTrackDX = new Float_t[multClust];
  Int_t *   cpvMatchedTrack = new Int_t[multClust];//number of ITSTPCtrackInAcc that matched i-th cluster
  for(Int_t i=0;i<multClust;i++)
    cpvMatchedTrack[i]=-1;
  Int_t *   trackMatchedCluster=new Int_t[nITSTPCtracksInAcc];//number of CPV cluster that matched i-th ITSTPCtrackInAcc
  Float_t *   trackMatchedClusterDist=new Float_t[nITSTPCtracksInAcc];
  for(Int_t i=0;i<nITSTPCtracksInAcc;i++){
    trackMatchedClusterDist[i]=999999.;
    trackMatchedCluster[i]=-1;
  }
  Float_t * phsX = new Float_t[multClust];
  Float_t * phsZ = new Float_t[multClust];
  Float_t * phsXM1 = new Float_t[multClust];
  Float_t * phsZM1 = new Float_t[multClust];
  Float_t * phsXM2 = new Float_t[multClust];
  Float_t * phsZM2 = new Float_t[multClust];
  Float_t * phsXM4 = new Float_t[multClust];
  Float_t * phsZM4 = new Float_t[multClust];
  Float_t * phsXEgt03 = new Float_t[multClust];
  Float_t * phsZEgt03 = new Float_t[multClust];
  Float_t * phsEEgt03 = new Float_t[multClust];
  Float_t * phsXMIP = new Float_t[multClust];
  Float_t * phsZMIP = new Float_t[multClust];
  
  Double_t minDistance = 1.e6;

  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() ) continue;
    isPHS = kFALSE;
    isCPV = kFALSE;

    if (clu1->GetType() == AliVCluster::kPHOSNeutral) {
      isPHS = kTRUE;
      isCPV = kFALSE;
    }
    if (clu1->GetType() == AliVCluster::kPHOSCharged) {
      isPHS = kFALSE;
      isCPV = kTRUE;
    }
    
    digMult = clu1->GetNCells();
    //if (digMult < 2) continue;
    clu1->GetPosition(position);
    TVector3 global(position) ;
    TVector3 local;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
    fPHOSGeo->Global2Local(local,global,mod1) ;
    
    cellAbsId = clu1->GetCellAbsId(0);
    if (isCPV){
      cellAbsId = fPHOSGeo->GetNModules() * fPHOSGeo->GetNCristalsInModule() - cellAbsId;
    
//=====================   Track matching   =======================//
      Double_t minDistance = 999999.,dx=999.,dz=999.;
      Double_t bz = AliTracker::GetBz() ; //B-Field for approximate matching
      Double_t b[3]; 

      for (Int_t i=0; i<nITSTPCtracksInAcc; i++) {
	AliESDtrack *track=ITSTPCtracksInAcc[i];
	ULong_t status = track->GetStatus();
	if ((status & AliESDtrack::kTPCout)   == 0) continue;
	const AliExternalTrackParam *outerParam = track->GetOuterParam();
	if (!outerParam) continue;
	Double_t z; 
	if(!outerParam->GetZAt(frCPV,bz,z)) continue ;
	if (TMath::Abs(z) > kZmax/2.) continue; // Some tracks miss the CPV module in Z
	
	
	AliExternalTrackParam t(*outerParam);
	if(!t.RotateParamOnly(kAlpha0)) continue ; // Rotate track to have CPV module direction angle = 0 degree.
	Double_t y;                       // Some tracks do not reach the CPV module
	if (!t.GetYAt(frCPV,bz,y)) continue; //    because of the bending
	
	if(TMath::Abs(y) < kXmax/2.){//note rotation here
	  //track is in CPV acceptance
	  t.GetBxByBz(b) ;
	  t.PropagateParamOnlyBxByBzTo(frCPV,b);        // Propagate to the matching module
	  Double_t gposTrack[3] ; 
	  t.GetXYZ(gposTrack) ;
	  TVector3 globalPositionTr(gposTrack) ;
	  TVector3 localPositionTr ;
	  fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,3) ;
	  //FillHistogram("hITSTPCtrackInAccXZ",localPositionTr.X(),localPositionTr.Z());
	  Double_t ddx = local.X()-localPositionTr.X();
	  Double_t ddz = local.Z()-localPositionTr.Z();
	  Double_t d2 = ddx*ddx + ddz*ddz;
	  if(d2 < minDistance) {
	    dx = ddx ;
	    dz = ddz ;
	    minDistance=d2 ;
	    if((TMath::Abs(dz-3.436)<3.*2.285)&&(TMath::Abs(dx-0.359)<3.*1.736)){//Track matched cluster
	      cpvMatchedTrack[multCPVClust]=i;
	      if(trackMatchedClusterDist[i]>d2){
		trackMatchedCluster[i]=multCPVClust;
		trackMatchedClusterDist[i]=d2;
	      }
	    } 
	  }
	}
      }
      cpvTrackDZ[multCPVClust]=dz;
      cpvTrackDX[multCPVClust]=dx;

    }
    // fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    // mod1   = relId[0];
    energy = clu1->E();
    tof    = clu1->GetTOF();
    dX     = dynamic_cast<AliESDCaloCluster*>(clu1)->GetTrackDx();
    dZ     = dynamic_cast<AliESDCaloCluster*>(clu1)->GetTrackDz();
    
    // Printf("Global=(%.2f,%.2f,%.2f), local=(%.2f,%.2f,%.2f), module=%d",
    // 	   global.X(),global.Y(),global.Z(),
    // 	   local.X(),local.Y(),local.Z(),mod1);
    
    if (isPHS) {
      multPHOSClust[0]++;
      FillHistogram("hPHSClusterEnergy",energy);
      FillHistogram("hPHSClusterEvsN",energy,digMult);
      FillHistogram("hPHSCellMultClu",digMult);
      FillHistogram("hPHSTrackDx"    ,dX);
      FillHistogram("hPHSTrackDz"    ,dZ);
      if      (mod1==1) {
	multPHOSClust[1]++;
	FillHistogram("hPHSClusterEvsNM1",energy,digMult);
	FillHistogram("hPHSCellMultCluM1",digMult);
	FillHistogram("hPHSClusterEnergyM1",energy);
	FillHistogram("hPHSClusNXZM1",local.X(),local.Z(),1.);
	FillHistogram("hPHSTrackDxM1",dX);
	FillHistogram("hPHSTrackDzM1",dZ);
	phsXM1[multPHOSClust[1]-1]=local.X();
	phsZM1[multPHOSClust[1]-1]=local.Z();

      }
      else if (mod1==2) {
	multPHOSClust[2]++;
	FillHistogram("hPHSClusterEvsNM2",energy,digMult);
	FillHistogram("hPHSCellMultCluM2",digMult);
	FillHistogram("hPHSClusterEnergyM2",energy);
	FillHistogram("hPHSClusNXZM2",local.X(),local.Z(),1.);
	FillHistogram("hPHSTrackDxM2",dX);
	FillHistogram("hPHSTrackDzM2",dZ);
	phsXM2[multPHOSClust[2]-1]=local.X();
	phsZM2[multPHOSClust[2]-1]=local.Z();
      }
      else if (mod1==3) {
	multPHOSClust[3]++;
	FillHistogram("hPHSClusterEvsNM3",energy,digMult);
	FillHistogram("hPHSCellMultCluM3",digMult);
	FillHistogram("hPHSClusterEnergyM3",energy);
	FillHistogram("hPHSClusNXZM3",local.X(),local.Z(),1.);
	FillHistogram("hPHSTrackDxM3",dX);
	FillHistogram("hPHSTrackDzM3",dZ);
	phsX[multPHOSClust[3]-1]=local.X();
	phsZ[multPHOSClust[3]-1]=local.Z();
	if(energy>0.3 && digMult>2){
	  multPHOSClustM3Egt03++;
	  FillHistogram("hPHSClusNXZM3Egt03",local.X(),local.Z(),1.);
	  for (Int_t i=0;i<5;i++)
	    if(energy>i && energy <i+1 )
	      FillHistogram(Form("hPHSClusNXZM3E%d-%d",i,i+1),local.X(),local.Z(),1.);
	  phsXEgt03[multPHOSClustM3Egt03-1]=local.X();
	  phsZEgt03[multPHOSClustM3Egt03-1]=local.Z();
	  phsEEgt03[multPHOSClustM3Egt03-1]=energy;
	}
	if(energy<0.3 && energy>0.15 ){
	  multPHOSClustM3MIP++;
	  FillHistogram("hPHSClusNXZM3MIP",local.X(),local.Z(),1.);
	  phsXMIP[multPHOSClustM3MIP-1]=local.X();
	  phsZMIP[multPHOSClustM3MIP-1]=local.Z();
	}
      }
      else if (mod1==4) {
	multPHOSClust[4]++;
	FillHistogram("hPHSClusterEvsNM4",energy,digMult);
	FillHistogram("hPHSCellMultCluM4",digMult);
	FillHistogram("hPHSClusterEnergyM4",energy);
	FillHistogram("hPHSClusNXZM4",local.X(),local.Z(),1.);
	FillHistogram("hPHSTrackDxM4",dX);
	FillHistogram("hPHSTrackDzM4",dZ);
	phsXM4[multPHOSClust[4]-1]=local.X();
	phsZM4[multPHOSClust[4]-1]=local.Z();
      }
    }
    if (isCPV) {
      multCPVClust++;
      FillHistogram("hCPVClusterEnergy",energy);
      FillHistogram("hCPVClusterEvsN",energy,digMult);
      FillHistogram("hCPVCellMultClu",digMult);
      FillHistogram("hCPVClusNXZ",local.X(),local.Z(),1.);
      FillHistogram("hCPVTrackDx",dX);
      FillHistogram("hCPVTrackDz",dZ);
      FillHistogram("hITSTPCtrackdXdZ",cpvTrackDX[multCPVClust-1],cpvTrackDZ[multCPVClust-1]);
      FillHistogram("hITSTPCtrackDist",TMath::Sqrt(pow(cpvTrackDX[multCPVClust-1],2)+pow(cpvTrackDZ[multCPVClust-1],2)));//hITSTPNegTrackdXdZ//
      
      cpvX[multCPVClust-1]=local.X();
      cpvZ[multCPVClust-1]=local.Z();
      cpvE[multCPVClust-1]=energy;
      cpvN[multCPVClust-1]=digMult;
    }
    FillHistogram("hPHSClusterMult"  ,multPHOSClust[0]);
    FillHistogram("hPHSClusterMultM1",multPHOSClust[1]);
    FillHistogram("hPHSClusterMultM2",multPHOSClust[2]);
    FillHistogram("hPHSClusterMultM3",multPHOSClust[3]);
    FillHistogram("hPHSClusterMultM4",multPHOSClust[4]);
    FillHistogram("hCPVClusterMult"  ,multCPVClust);
  }
  //for every ITSTPC track in cpv acceptance and cpv cluster check their correlations and so on
  for(Int_t iClu=0;iClu<multCPVClust;iClu++){
    if(cpvMatchedTrack[iClu]>=0){//this cluster has matched track
      FillHistogram("hCpvCluMatchedXZ",cpvX[iClu],cpvZ[iClu]);
      FillHistogram("hCpvCluMatchedEn",cpvE[iClu]);
      FillHistogram("hCpvCluMatchedN",cpvN[iClu]);
      FillHistogram("hCpvCluMatchedEvsN",cpvE[iClu],cpvN[iClu]);
    }
    else{//this cluster has no matched track
      FillHistogram("hCpvCluNonMatchedXZ",cpvX[iClu],cpvZ[iClu]);
      FillHistogram("hCpvCluNonMatchedEn",cpvE[iClu]);
      FillHistogram("hCpvCluNonMatchedN",cpvN[iClu]);
      FillHistogram("hCpvCluNonMatchedEvsN",cpvE[iClu],cpvN[iClu]);
      //unmatched clusters from hot zones
      if(IsHotZone(cpvX[iClu],cpvZ[iClu])){
	FillHistogram("hCpvCluHotNonMatchedXZ",cpvX[iClu],cpvZ[iClu]);
	FillHistogram("hCpvCluHotNonMatchedEn",cpvE[iClu]);
	FillHistogram("hCpvCluHotNonMatchedN",cpvN[iClu]);
	FillHistogram("hCpvCluHotNonMatchedEvsN",cpvE[iClu],cpvN[iClu]);

      }
      else{//unmatched clusters from cold zones
	FillHistogram("hCpvCluColdNonMatchedXZ",cpvX[iClu],cpvZ[iClu]);
	FillHistogram("hCpvCluColdNonMatchedEn",cpvE[iClu]);
	FillHistogram("hCpvCluColdNonMatchedN",cpvN[iClu]);
	FillHistogram("hCpvCluColdNonMatchedEvsN",cpvE[iClu],cpvN[iClu]);

      }
    }
  }
  for(Int_t iTr=0;iTr<nITSTPCtracksInAcc;iTr++){
    AliESDtrack *track=ITSTPCtracksInAcc[iTr];
    Double_t P = track->P();
    Double_t E = track->E();
    Int_t PID = track->GetPID();
    if(trackMatchedCluster[iTr]>=0){//this track has matched cluster
      FillHistogram("hITSTPCtrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
      FillHistogram("hITSTPCtrackMatchedInAccP",P);
      FillHistogram("hITSTPCtrackMatchedInAccPID",PID);
      if(isElectron[iTr]){//this is electron
	FillHistogram("hITSTPCelectronMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCelectronMatchedInAccP",P);
	FillHistogram("hITSTPCelectrondXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCelectronDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isHMtrack[iTr]){//this is HMtrack
	FillHistogram("hITSTPCHMtrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCHMtrackMatchedInAccP",P);
	FillHistogram("hITSTPCHMtrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCHMtrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isPosTrack[iTr]){//this is PosTrack
	FillHistogram("hITSTPCPosTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosTrackMatchedInAccP",P);
	FillHistogram("hITSTPCPosTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCPosTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isNegTrack[iTr]){//this is NegTrack
	FillHistogram("hITSTPCNegTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegTrackMatchedInAccP",P);
	FillHistogram("hITSTPCNegTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCNegTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isPosZTrack[iTr]){//this is PosZTrack
	FillHistogram("hITSTPCPosZTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosZTrackMatchedInAccP",P);
	FillHistogram("hITSTPCPosZTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCPosZTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isNegZTrack[iTr]){//this is NegZTrack
	FillHistogram("hITSTPCNegZTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegZTrackMatchedInAccP",P);
	FillHistogram("hITSTPCNegZTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCNegZTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isPosZTrack[iTr]&&isHMtrack[iTr]){//this is PosZHMTrack
	FillHistogram("hITSTPCPosZHMTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosZHMTrackMatchedInAccP",P);
	FillHistogram("hITSTPCPosZHMTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCPosZHMTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isNegZTrack[iTr]&&isHMtrack[iTr]){//this is NegZHMTrack
	FillHistogram("hITSTPCNegZHMTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegZHMTrackMatchedInAccP",P);
	FillHistogram("hITSTPCNegZHMTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hITSTPCNegZHMTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isNegZTrack[iTr]&&isNegTrack[iTr]){//this is NegZNegTrack
	FillHistogram("hNegZNegTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hNegZNegTrackMatchedInAccP",P);
	FillHistogram("hNegZNegTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hNegZNegTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isNegZTrack[iTr]&&isPosTrack[iTr]){//this is NegZPosTrack
	FillHistogram("hNegZPosTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hNegZPosTrackMatchedInAccP",P);
	FillHistogram("hNegZPosTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hNegZPosTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isPosZTrack[iTr]&&isPosTrack[iTr]){//this is PosZPosTrack
	FillHistogram("hPosZPosTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hPosZPosTrackMatchedInAccP",P);
	FillHistogram("hPosZPosTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hPosZPosTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }
      if(isPosZTrack[iTr]&&isNegTrack[iTr]){//this is NegZNegTrack
	FillHistogram("hPosZNegTrackMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hPosZNegTrackMatchedInAccP",P);
	FillHistogram("hPosZNegTrackdXdZ",cpvTrackDX[trackMatchedCluster[iTr]],cpvTrackDZ[trackMatchedCluster[iTr]]);
	FillHistogram("hPosZNegTrackDist",TMath::Sqrt(pow(cpvTrackDX[trackMatchedCluster[iTr]],2)+pow(cpvTrackDZ[trackMatchedCluster[iTr]],2)));
      }


      
    }
    else{////this track has no matched cluster
      FillHistogram("hITSTPCtrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
      FillHistogram("hITSTPCtrackNonMatchedInAccP",P);
      FillHistogram("hITSTPCtrackNonMatchedInAccPID",PID);
      if(isElectron[iTr]){//this is electron
	FillHistogram("hITSTPCelectronNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCelectronNonMatchedInAccP",P);
      }
      if(isHMtrack[iTr]){//this is HMtrack
	FillHistogram("hITSTPCHMtrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCHMtrackNonMatchedInAccP",P);
      }
      if(isPosTrack[iTr]){//this is PosTrack
	FillHistogram("hITSTPCPosTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosTrackNonMatchedInAccP",P);
      }
      if(isNegTrack[iTr]){//this is NegTrack
	FillHistogram("hITSTPCNegTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegTrackNonMatchedInAccP",P);
      }
      if(isPosZTrack[iTr]){//this is PosZTrack
	FillHistogram("hITSTPCPosZTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosZTrackNonMatchedInAccP",P);
      }
      if(isNegZTrack[iTr]){//this is NegZTrack
	FillHistogram("hITSTPCNegZTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegZTrackNonMatchedInAccP",P);
      }
      if(isPosZTrack[iTr]&&isHMtrack[iTr]){//this is PosZTrack
	FillHistogram("hITSTPCPosZHMTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCPosZHMTrackNonMatchedInAccP",P);
      }
      if(isNegZTrack[iTr]&&isHMtrack[iTr]){//this is NegZTrack
	FillHistogram("hITSTPCNegZHMTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
	FillHistogram("hITSTPCNegZHMTrackNonMatchedInAccP",P);
      }
      if(isPosZTrack[iTr]&&isPosTrack[iTr]){//this is PosZ positive Track
        FillHistogram("hPosZPosTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
        FillHistogram("hPosZPosTrackNonMatchedInAccP",P);
      }
      if(isNegZTrack[iTr]&&isPosTrack[iTr]){//this is NegZ positive Track
        FillHistogram("hNegZPosTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
        FillHistogram("hNegZPosTrackNonMatchedInAccP",P);
      }
      if(isPosZTrack[iTr]&&isNegTrack[iTr]){//this is PosZ negative Track
        FillHistogram("hPosZNegTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
        FillHistogram("hPosZNegTrackNonMatchedInAccP",P);
      }
      if(isNegZTrack[iTr]&&isNegTrack[iTr]){//this is NegZ negative Track
        FillHistogram("hNegZNegTrackNonMatchedInAccXZ",trackInAccX[iTr],trackInAccZ[iTr]);
        FillHistogram("hNegZNegTrackNonMatchedInAccP",P);
      }


    }
  }


  //for every cpv cluster find nearest phos cluster in M3 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClust[3];iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsX[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZ[iCluPHS],2));
      if(R<minR) {
	minR=R;
	minX=cpvX[iCluCPV]-phsX[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZ[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM3corrR",minR);
    FillHistogram("hCPVPHOSM3corrX",minX);
    FillHistogram("hCPVPHOSM3corrZ",minZ);
    FillHistogram("hCPVPHOSM3corrXZ",minX,minZ,1.);

    //FillHistogram("");
  }

  //for every cpv cluster find nearest phos cluster in M3 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    Int_t iCluPHSmin = -1;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClustM3Egt03;iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsXEgt03[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZEgt03[iCluPHS],2));
      if(R<minR) {
	iCluPHSmin=iCluPHS;
	minR=R;
	minX=cpvX[iCluCPV]-phsXEgt03[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZEgt03[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM3corrREgt03",minR);
    FillHistogram("hCPVPHOSM3corrXEgt03",minX);
    FillHistogram("hCPVPHOSM3corrZEgt03",minZ);
    FillHistogram("hCPVPHOSM3corrXZEgt03",minX,minZ,1.);
    for (Int_t i = 0; i<5;i++)
      if(iCluPHSmin>=0 )
	if(phsEEgt03[iCluPHSmin]>i && phsEEgt03[iCluPHSmin]<i+1){
	  FillHistogram(Form("hCPVPHOSM3corrRE%d-%d",i,i+1),minR);
	  FillHistogram(Form("hCPVPHOSM3corrXE%d-%d",i,i+1),minX);
	  FillHistogram(Form("hCPVPHOSM3corrZE%d-%d",i,i+1),minZ);
	  FillHistogram(Form("hCPVPHOSM3corrXZE%d-%d",i,i+1),minX,minZ,1.);
			
	}
  }

  //for every cpv cluster find nearest phos cluster in M3 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClustM3MIP;iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsXMIP[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZMIP[iCluPHS],2));
      if(R<minR) {
	minR=R;
	minX=cpvX[iCluCPV]-phsXMIP[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZMIP[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM3corrRMIP",minR);
    FillHistogram("hCPVPHOSM3corrXMIP",minX);
    FillHistogram("hCPVPHOSM3corrZMIP",minZ);
    FillHistogram("hCPVPHOSM3corrXZMIP",minX,minZ,1.);
  }

  //for every cpv cluster find nearest phos cluster in M1 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClust[1];iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsXM1[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZM1[iCluPHS],2));
      if(R<minR) {
	minR=R;
	minX=cpvX[iCluCPV]-phsXM1[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZM1[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM1corrR",minR);
    FillHistogram("hCPVPHOSM1corrX",minX);
    FillHistogram("hCPVPHOSM1corrZ",minZ);
    FillHistogram("hCPVPHOSM1corrXZ",minX,minZ,1.);

  }

  //for every cpv cluster find nearest phos cluster in M2 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClust[2];iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsXM2[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZM2[iCluPHS],2));
      if(R<minR) {
	minR=R;
	minX=cpvX[iCluCPV]-phsXM2[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZM2[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM2corrR",minR);
    FillHistogram("hCPVPHOSM2corrX",minX);
    FillHistogram("hCPVPHOSM2corrZ",minZ);
    FillHistogram("hCPVPHOSM2corrXZ",minX,minZ,1.);

  }

  //for every cpv cluster find nearest phos cluster in M4 and check their correlation
  for (Int_t iCluCPV=0;iCluCPV<multCPVClust;iCluCPV++){
    Float_t minR = 10000., minX = 10000.,minZ = 10000.,R;
    for (Int_t iCluPHS=0;iCluPHS<multPHOSClust[4];iCluPHS++){
      R=sqrt(pow(cpvX[iCluCPV]-phsXM4[iCluPHS],2)+pow(cpvZ[iCluCPV]-phsZM4[iCluPHS],2));
      if(R<minR) {
	minR=R;
	minX=cpvX[iCluCPV]-phsXM4[iCluPHS];
	minZ=cpvZ[iCluCPV]-phsZM4[iCluPHS];
      }
    }
    FillHistogram("hCPVPHOSM4corrR",minR);
    FillHistogram("hCPVPHOSM4corrX",minX);
    FillHistogram("hCPVPHOSM4corrZ",minZ);
    FillHistogram("hCPVPHOSM4corrXZ",minX,minZ,1.);

  }

  AliDebug(2,Form("Event %d, trig.class %s, period %d, bc %d, orbit %d, N PHS clusters: %d, N CPV clusters: %d",
		  eventNumberInFile,trigClasses.Data(),event->GetPeriodNumber(),
		  event->GetBunchCrossNumber(),event->GetOrbitNumber(),multPHOSClust[0],multCPVClust));
  
  // Post output data.
  PostData(1, fOutputContainer1);
  PostData(2, fOutputContainer2);
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisTaskCPV::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

}

//_____________________________________________________________________________
void AliAnalysisTaskCPV::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1 * hist1 = dynamic_cast<TH1*>(fOutputContainer1->FindObject(key)) ;
  TH1 * hist2 = dynamic_cast<TH1*>(fOutputContainer2->FindObject(key)) ;
  if(hist1)
    hist1->Fill(x) ;
  else if(hist2)
    hist2->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskCPV::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  TH1 * th11 = dynamic_cast<TH1*> (fOutputContainer1->FindObject(key));
  TH1 * th12 = dynamic_cast<TH1*> (fOutputContainer2->FindObject(key));
  if(th11)
    th11->Fill(x, y) ;
  else if (th12)
    th12->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskCPV::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const
{
  //Fills 1D histograms with key
  TObject * obj1 = fOutputContainer1->FindObject(key);
  TObject * obj2 = fOutputContainer2->FindObject(key);
  TH2 * th2 = dynamic_cast<TH2*> (obj1);
  if(!obj1) th2 = dynamic_cast<TH2*> (obj2);
  if(th2) {
    th2->Fill(x, y, z) ;
    return;
  }
  
  TH3 * th3 = dynamic_cast<TH3*> (obj1);
  if(!th3) th3 = dynamic_cast<TH3*> (obj2);
  if(th3) {
    th3->Fill(x, y, z) ;
    return;
  }
  
  AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}
AliESDtrackCuts* AliAnalysisTaskCPV::ITSTPCTrackCuts()
{
  // cuts for normal tracks (ITS + TPC)
  // i.e. GetStandardITSTPCTrackCuts2010(kTRUE);
  // (same, just typed in full detail...)
  Int_t clusterCut = 1;
  AliESDtrackCuts* esdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE, clusterCut);
  esdTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);

  return esdTrackCuts;

}
Bool_t AliAnalysisTaskCPV::IsHotZone(Double_t x, Double_t z)
{
  Bool_t hotZone = 0;
  if((x>-71.  && x<-62.) && (z> 28.4  && z<  39.2))hotZone=1;
  if((x> 20.8 && x<28. ) && (z> -60.3 && z< -47.4))hotZone=1;
  return hotZone;
}
