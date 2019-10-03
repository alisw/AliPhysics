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
/* $Id: AliAnalysisTaskPi0.cxx 55426 2012-03-30 08:09:28Z kharlov $ */
 
// Analysis task for pi0 and eta meson analysis in pp collisions
// Authors: Yuri Kharlov, Dmitri Peressounko

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
#include "AliAnalysisTaskPi0PP.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"
#include "AliOADBContainer.h"

// Analysis task to fill histograms with PHOS AOD clusters and cells
// Authors: Pooja Pareek
// Date   : 26.01.2016
// Date   : 29.01.2016
ClassImp(AliAnalysisTaskPi0PP)

//________________________________________________________________________
AliAnalysisTaskPi0PP::AliAnalysisTaskPi0PP(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),
  fPHOSEvent(0),
  fnCINT1B(0),
  fnCINT1A(0),
  fnCINT1C(0),
  fnCINT1E(0),
  fBCgap(525e-09),
  fPHOSGeo(0),
  fEventCounter(0),
  fTriggerAnalysis(new AliTriggerAnalysis)
{
  // Constructor
  Int_t nBin=10 ;
  for(Int_t i=0;i<nBin;i++){
    for(Int_t j=0;j<2;j++)
      fPHOSEvents[i][j]=0 ;
  }
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
  AliOADBContainer geomContainer("phosGeo");  
  geomContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSMCGeometry.root","PHOSMCRotationMatrixes");
  Int_t runNumber=117116; //LHC10b 
  TObjArray *matrixes = (TObjArray*)geomContainer.GetObject(runNumber,"PHOSRotationMatrixes");
        for(Int_t mod=0; mod<6; mod++) {
          if(!matrixes->At(mod)) continue;
          fPHOSGeo->SetMisalMatrix(((TGeoHMatrix*)matrixes->At(mod)),mod) ;
          printf(".........Adding Matrix(%d), geo=%p\n",mod,fPHOSGeo) ;
          ((TGeoHMatrix*)matrixes->At(mod))->Print() ;
	}
  
  // Absolute recalibration for LHC11a. Use SetRecalib(mod,recalib) to change it
  fRecalib[0] = 0.9822696;// 0.9942*0.988;
  fRecalib[1] = 0.9861288;//0.9822*1.004;
  fRecalib[2] = 1.0072;

}

//________________________________________________________________________
void AliAnalysisTaskPi0PP::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // AOD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));
  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     ,100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMult"  ,"PHOS cluster multiplicity"    ,100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM1","PHOS cluster multiplicity, M1",100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM2","PHOS cluster multiplicity, M2",100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM3","PHOS cluster multiplicity, M3",100,0,100));
  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,500,0.,50.));
  fOutputContainer->Add(new TH1F("hCellEnergyM1","Cell energy in module 1",500,0.,50.));
  fOutputContainer->Add(new TH1F("hCellEnergyM2","Cell energy in module 2",500,0.,50.));
  fOutputContainer->Add(new TH1F("hCellEnergyM3","Cell energy in module 3",500,0.,50.));
  fOutputContainer->Add(new TH1F("hClusterEnergy"  ,"Cluster energy"      ,500,0.,50.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM1","Cluster energy, M1"  ,500,0.,50.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM2","Cluster energy, M2"  ,500,0.,50.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM3","Cluster energy, M3"  ,500,0.,50.));
  fOutputContainer->Add(new TH2F("hClusterEvsN"  ,"Cluster energy vs digit multiplicity"    ,500,0.,50.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM1","Cluster energy vs digit multiplicity, M1",500,0.,50.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM2","Cluster energy vs digit multiplicity, M2",500,0.,50.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM3","Cluster energy vs digit multiplicity, M3",500,0.,50.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsTM1","Cluster energy vs time, M1", 500,0.,50., 1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hClusterEvsTM2","Cluster energy vs time, M2", 500,0.,50., 1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hClusterEvsTM3","Cluster energy vs time, M3", 500,0.,50., 1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH1I("hCellMultClu"  ,"Cell multiplicity per cluster"    ,200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM1","Cell multiplicity per cluster, M1",200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM2","Cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM3","Cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer->Add(new TH1I("hModule","Module events",5,0.,5.));
  // fOutputContainer->Add(new TH1F("hSelEvents","Selected events",8,-0.5,7.5));
  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",9,-0.5,8.5));
  fOutputContainer->Add(new TH1F("hSelEventspi0","Selected events for pi0",9,-0.5,8.5));//new to check event selection

  fOutputContainer->Add(new TH2F("hCellNXZM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM1","Cell E(X,Z), M1",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM2","Cell E(X,Z), M2",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM3","Cell E(X,Z), M3",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM1_0","Clu (X,Z), M1, E>0.5 GeV"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM2_0","Clu (X,Z), M2, E>0.5 GeV"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM3_0","Clu (X,Z), M3, E>0.5 GeV"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM1_0","Clu E(X,Z), M1, E>0.5 GeV"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM2_0","Clu E(X,Z), M2, E>0.5 GeV"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM3_0","Clu E(X,Z), M3, E>0.5 GeV"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM1_1","Clu (X,Z), M1, E>1 GeV"     ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM2_1","Clu (X,Z), M2, E>1 GeV"     ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM3_1","Clu (X,Z), M3, E>1 GeV"     ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM1_1","Clu E(X,Z), M1, E>1 GeV"    ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM2_1","Clu E(X,Z), M2, E>1 GeV"    ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM3_1","Clu E(X,Z), M3, E>1 GeV"    ,64,0.5,64.5, 56,0.5,56.5));

  Int_t nM       = 750;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.5;
  Int_t nPt      = 400;
  Double_t ptMin = 0;
  Double_t ptMax = 40;
  fOutputContainer->Add(new TH2F("hAsymPtPi0","(A,p_{T})_{#gamma#gamma} #pi^{0}"     ,20,0.,1.,    40,0.,20.));
  fOutputContainer->Add(new TH2F("hAsymPtEta","(A,p_{T})_{#gamma#gamma} #eta"        ,20,0.,1.,    40,0.,20.));

  fOutputContainer->Add(new TH2F("hAsymPtPi0M1" ,"(A,p_{T})_{#gamma#gamma} #pi^{0}. M1"  ,20,0.,1.,    40,0.,20.));
  fOutputContainer->Add(new TH2F("hAsymPtPi0M2" ,"(A,p_{T})_{#gamma#gamma} #pi^{0}. M2"  ,20,0.,1.,    40,0.,20.));
  fOutputContainer->Add(new TH2F("hAsymPtPi0M3" ,"(A,p_{T})_{#gamma#gamma} #pi^{0}. M3"  ,20,0.,1.,    40,0.,20.));
  fOutputContainer->Add(new TH2F("hAsymPtPi0M12","(A,p_{T})_{#gamma#gamma} #pi^{0}. M12" ,20,0.,1.,    40,0.,20.));
  fOutputContainer->Add(new TH2F("hAsymPtPi0M23","(A,p_{T})_{#gamma#gamma} #pi^{0}. M23" ,20,0.,1.,    40,0.,20.));

  fOutputContainer->Add(new TH2F("hMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10BC0" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1=BC2=0",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA10BC1" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1!=BC2" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA10BC2" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1=0"    ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10vtx10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, |Zvtx|<10 cm"     ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07vtx10","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, |Zvtx|<10 cm"     ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, AOD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, AOD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH3F("hMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtM1","(M,p_{T})_{#gamma#gamma}, module 1"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM2","(M,p_{T})_{#gamma#gamma}, module 2"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM3","(M,p_{T})_{#gamma#gamma}, module 3"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM12","(M,p_{T})_{#gamma#gamma}, modules 1,2",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM13","(M,p_{T})_{#gamma#gamma}, modules 1,3",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM23","(M,p_{T})_{#gamma#gamma}, modules 2,3",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPt20cm","(M,p_{T})_{#gamma#gamma}, |z|<20 cm"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPt40cm","(M,p_{T})_{#gamma#gamma}, 20<|z|<40 cm",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPt60cm","(M,p_{T})_{#gamma#gamma}, |z|>40 cm"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiAsymPt","(A,p_{T})_{#gamma#gamma}"                ,50,0.,1.,    nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10BC0" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1=BC2=0",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA10BC1" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1!=BC2" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA10BC2" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, BC1=0"    ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10vtx10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, |Zvtx|<10 cm",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07vtx10","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, |Zvtx|<10 cm",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, AOD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, AOD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH3F("hMiMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMiMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtM1","(M,p_{T})_{#gamma#gamma}, module 1"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM2","(M,p_{T})_{#gamma#gamma}, module 2"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM3","(M,p_{T})_{#gamma#gamma}, module 3"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM12","(M,p_{T})_{#gamma#gamma}, modules 1,2",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM13","(M,p_{T})_{#gamma#gamma}, modules 1,3",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM23","(M,p_{T})_{#gamma#gamma}, modules 2,3",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPt20cm","(M,p_{T})_{#gamma#gamma}, |z|<20 cm"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPt40cm","(M,p_{T})_{#gamma#gamma}, 20<|z|<40 cm",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPt60cm","(M,p_{T})_{#gamma#gamma}, |z|>40 cm"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH1F("hPhotonKappa","#kappa(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPt","p_{T}(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPx","p_{x}(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPy","p_{y}(#gamma)",200,0.,20.));

  fOutputContainer->Add(new TH1F("hTrigClass","Trigger class",5,0.5,5.5));

  fOutputContainer->Add(new TH1F("hNPileupVtx","Number of SPD pileup vertices",10,0.,10.));
//   fOutputContainer->Add(new TH1F("hZPileupVtx","#Delta_{Z} vtx_{0}-vtx_{PU}",200,-50.,+50.));

  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity",150,0.,150.));

  fOutputContainer->Add(new TH1F("hV0Atime","V0A time",1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH1F("hV0Ctime","V0C time",1200,-6.e-6,+6.e-6));
  fOutputContainer->Add(new TH2F("hV0AV0Ctime","V0A time vs V0C time",120,-6.e-6,+6.e-6 ,120,-6.e-6,+6.e-6));

  // Create AOD track cut

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPi0PP::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze AOD

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  //Skip events from fast cluster

  // UInt_t triggerMask = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected());
  // if(triggerMask& AliVEvent::kFastOnly) return; // reject events in the fast cluster only
  // if(!(triggerMask& AliVEvent::kMB)) return; // check the trigger mask as usual

  FillHistogram("hSelEvents",0) ; // All events accepted by PSel

  TString trigClasses = event->GetFiredTriggerClasses();
  if (trigClasses.Contains("FAST")  && !trigClasses.Contains("ALL")) {
    AliWarning(Form("Skip event with triggers %s",trigClasses.Data()));
    return;
  }

  // Event selection flags

  Bool_t eventVtxExist    = kFALSE;
  Bool_t eventVtxZ10cm    = kFALSE;
  Bool_t eventPileup      = kFALSE;
  Bool_t eventV0AND       = kFALSE;

  Int_t eventNumberInFile = event->GetEventNumberInFile();
  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",50) ;

  // Checks if we have a primary vertex
  // Get primary vertices form AOD

  // if(event->GetPrimaryVertexTracks())
  //   if(event->GetPrimaryVertexTracks()->GetNContributors()>0)
  //     eventVtxExist    = kTRUE;

  // if (event->GetPrimaryVertexSPD())
  //   if(event->GetPrimaryVertexSPD()->GetNContributors()>0)
  //     eventVtxExist    = kTRUE;
  
  if (event->GetPrimaryVertex())
    if(event->GetPrimaryVertex()->GetNContributors()>0)
      eventVtxExist    = kTRUE;

  const AliAODVertex *esdVertexBest = event->GetPrimaryVertex();
  const AliAODVertex *esdVertexSPD  = event->GetPrimaryVertexSPD();

  Double_t vtx0[3] = {0,0,0}; // don't rely on AOD vertex, assume (0,0,0)
  Double_t vtxBest[3];
  vtxBest[0] = esdVertexBest->GetX();
  vtxBest[1] = esdVertexBest->GetY();
  vtxBest[2] = esdVertexBest->GetZ();

  FillHistogram("hNvertexTracks",esdVertexBest->GetNContributors());
  
 
  if (event->GetPrimaryVertex() && event->GetPrimaryVertex()->GetNContributors()>0)
    {
      FillHistogram("hZvertex"      ,esdVertexBest->GetZ());
      if (TMath::Abs(esdVertexBest->GetZ()) < 10. )
	eventVtxZ10cm = kTRUE;
    }

  // Check for pileup and fill pileup histograms
  if (event->IsPileupFromSPD()) {
    eventPileup = kTRUE;
    //     TClonesArray *pileupVertices = event->GetPileupVerticesSPD();
    //     Int_t nPileupVertices = pileupVertices->GetEntriesFast();
    FillHistogram("hNPileupVtx",event->GetNumberOfPileupVerticesSPD());
    //     for (Int_t puVtx=0; puVtx<nPileupVertices; puVtx++) {
    //       Double_t dZpileup = esdVertexSPD->GetZ() - event->GetPileupVertexSPD(puVtx)->GetZ();
    //       FillHistogram("hZPileupVtx",dZpileup);
    //     }
  }

  //   eventV0AND = fTriggerAnalysis->IsOfflineTriggerFired(event, AliTriggerAnalysis::kV0AND);
  //   
  //   ULong64_t GetTriggerMask()        const { return fHeader ? fHeader->GetTriggerMask() : 0;}
  //   UChar_t   GetTriggerCluster()     const { return fHeader ? fHeader->GetTriggerCluster() : 0;}
  //   TString   GetFiredTriggerClasses()const { return fHeader->GetFiredTriggerClasses();};
  //   
  //  // Checks if corresponding bit in mask is on
  //   ULong64_t trigmask = aEsd->GetTriggerMask();
  //   return (trigmask & (1ull << (tclass-1)));
  // } 
  
  ULong64_t trigmask = event->GetTriggerMask();
  
  eventV0AND = (trigmask & (1ull << (AliTriggerAnalysis::kV0AND-1)));

  // Fill event statistics for different selection criteria

  FillHistogram("hSelEvents",1) ;
  if (eventVtxExist) 
    FillHistogram("hSelEvents",2) ;
  if (eventVtxExist && eventVtxZ10cm)
    FillHistogram("hSelEvents",3) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND)
    FillHistogram("hSelEvents",4) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND && eventPileup)
    FillHistogram("hSelEvents",5) ;
  if (eventPileup)
    FillHistogram("hSelEvents",6) ;
  if(eventV0AND)
    FillHistogram("hSelEvents",7) ;
  if(eventVtxZ10cm)
    FillHistogram("hSelEvents",8) ;  
  
  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtxBest[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;
  
  if (trigClasses.Contains("CINT1B")) fnCINT1B++;
  if (trigClasses.Contains("CINT1A")) fnCINT1A++;
  if (trigClasses.Contains("CINT1C")) fnCINT1C++;
  if (trigClasses.Contains("CINT1-E")) fnCINT1E++;
  
  //Calculate charged multiplicity
  //   Int_t trackMult = 0;
  //   for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
  //     AliAODtrack *track = new AliAODtrack(*event->GetTrack(i)) ;
  //     if(fAODtrackCuts->AcceptTrack(track) &&  TMath::Abs(track->Eta())< 0.8)
  //       trackMult++;
  //     delete track;
  //   }
  Int_t trackMult =event->GetNumberOfTracks() ;
  FillHistogram("hTrackMult",trackMult+0.5) ;
  
  Float_t tV0A = event->GetVZEROData()->GetV0ATime();
  Float_t tV0C = event->GetVZEROData()->GetV0CTime();
  FillHistogram("hV0Atime",tV0A);
  FillHistogram("hV0Atime",tV0C);
  FillHistogram("hV0AV0Ctime",tV0A,tV0C);

  Int_t centr=0 ;
  //always zero centrality
  if(!fPHOSEvents[zvtx][centr]) fPHOSEvents[zvtx][centr]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][centr] ;
  
  if(trackMult<=2)
    centr=0 ;
  else 
    if(trackMult<=5)
      centr=1 ;
    else
      if(trackMult<=9)
        centr=2 ;
      else
        if(trackMult<=14)
          centr=3 ;
        else
          if(trackMult<=22)
            centr=4 ;
          else
            if(trackMult<=35)
              centr=5 ;
            else
              if(trackMult<=50)
                centr=6 ;
              else
                centr=7 ;
  
  //check what events are selected for pi0 cluster pairs
  if(esdVertexBest->GetNContributors()<1.0 || esdVertexSPD->GetNContributors()<1.0 ) return;
  if (!eventVtxExist) return;
  if (!eventVtxZ10cm) return;
  if (eventPileup) return;
  
  FillHistogram("hSelEventspi0",1) ;
  if (eventVtxExist) 
    FillHistogram("hSelEventspi0",2) ;
  if (eventVtxExist && eventVtxZ10cm)
    FillHistogram("hSelEventspi0",3) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND)
    FillHistogram("hSelEventspi0",4) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND && eventPileup)
    FillHistogram("hSelEventspi0",5) ;
  if (eventPileup)
    FillHistogram("hSelEventspi0",6) ;
  if(eventV0AND)
    FillHistogram("hSelEventspi0",7) ;
  if(eventVtxZ10cm)
    FillHistogram("hSelEventspi0",8) ;  
 


  AliAODCaloCluster *clu1;
  TLorentzVector p1,p2,p12, pv1,pv2,pv12;
  AliAODCaloCells *cells      = event->GetPHOSCells();

  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t multCells = cells->GetNumberOfCells();
  FillHistogram("hClusterMult",multClust);
  FillHistogram("hCellMultEvent",multCells);

  //   Printf("Event %d, trig.class %s, period %d, bc %d, orbit %d",
  //      eventNumberInFile,trigClasses.Data(),event->GetPeriodNumber(),
  //      event->GetBunchCrossNumber(),event->GetOrbitNumber());
  //   Printf("\tthere are %d caloclusters and %d calocells",
  //      multClust,multCells);

  // Get PHOS rotation matrices from AOD and set them to the PHOS geometry
  // if(fEventCounter == 0) {
  //   for(Int_t mod=0; mod<5; mod++) {
  //     if(!event->GetPHOSMatrix(mod)) continue;
  //     fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
  //     Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
  //   }
  // }

  Float_t  energy, tof;
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells

  Int_t nCellModule[3] = {0,0,0};
  for (Int_t iCell=0; iCell<multCells; iCell++) {
    cellAbsId = cells->GetCellNumber(iCell);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1  = relId[0];
    cellX = relId[2];
    cellZ = relId[3] ;
    energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if      (mod1==1) {
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
  }
  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  
  // Single loop over clusters fills cluster histograms
  
  Int_t    digMult;
  Int_t    multPHOSClust[4]  = {0,0,0,0};
  Float_t  position[3];
  
  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() ) continue;

    digMult = clu1->GetNCells();
    clu1->GetPosition(position);
    TVector3 global1(position) ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod1,cellX,cellZ) ) continue ;
    
    cellAbsId = clu1->GetCellAbsId(0);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1   = relId[0];
    energy = clu1->E();
    tof    = clu1->GetTOF();

    // Printf("\tmodule=%d, xyz=(%.3f,%.3f,%.3f) cm, E=%.3f GeV",
    // 	   mod1,position[0],position[1],position[2],energy);
    
    multPHOSClust[0]++;
    FillHistogram("hClusterEnergy",energy);
    FillHistogram("hClusterEvsN",energy,digMult);
    FillHistogram("hCellMultClu",digMult);
    if      (mod1==1) {
      multPHOSClust[1]++;
      FillHistogram("hClusterEvsNM1",energy,digMult);
      FillHistogram("hClusterEvsTM1",energy,tof);
      FillHistogram("hCellMultCluM1",digMult);
      FillHistogram("hClusterEnergyM1",energy);
      if (energy > 0.5) {
	FillHistogram("hCluNXZM1_0",cellX,cellZ,1.);
	FillHistogram("hCluEXZM1_0",cellX,cellZ,energy);
      }
      if (energy > 1.0) {
	FillHistogram("hCluNXZM1_1",cellX,cellZ,1.);
	FillHistogram("hCluEXZM1_1",cellX,cellZ,energy);
      }
    }
    else if (mod1==2) {
      multPHOSClust[2]++;
      FillHistogram("hClusterEvsNM2",energy,digMult);
      FillHistogram("hClusterEvsTM2",energy,tof);
      FillHistogram("hCellMultCluM2",digMult);
      FillHistogram("hClusterEnergyM2",energy);
      if (energy > 0.5) {
	FillHistogram("hCluNXZM2_0",cellX,cellZ,1.);
	FillHistogram("hCluEXZM2_0",cellX,cellZ,energy);
      }
      if (energy > 1.0) {
	FillHistogram("hCluNXZM2_1",cellX,cellZ,1.);
	FillHistogram("hCluEXZM2_1",cellX,cellZ,energy);
      }
    }
    else if (mod1==3) {
      multPHOSClust[3]++;
      FillHistogram("hClusterEvsNM3",energy,digMult);
      FillHistogram("hClusterEvsTM3",energy,tof);
      FillHistogram("hCellMultCluM3",digMult);
      FillHistogram("hClusterEnergyM3",energy);
      if (energy > 0.5) {
	FillHistogram("hCluNXZM3_0",cellX,cellZ,1.);
	FillHistogram("hCluEXZM3_0",cellX,cellZ,energy);
      }
      if (energy > 1.0) {
	FillHistogram("hCluNXZM3_1",cellX,cellZ,1.);
	FillHistogram("hCluEXZM3_1",cellX,cellZ,energy);
      }
    }
    
    if (digMult > 2) {
      clu1 ->GetMomentum(p1 ,vtx0);
      Double_t pAbs = p1.P();
      Double_t pT   = p1.Pt();
      Double_t pX   = p1.Px();
      Double_t pY   = p1.Py();
      if (pAbs<1.e-10) break;
      Double_t kappa = pAbs - TMath::Power(0.135,2)/4./pAbs;
      
      FillHistogram("hPhotonKappa",kappa);
      FillHistogram("hPhotonPt",pT);
      FillHistogram("hPhotonPx",pX);
      FillHistogram("hPhotonPy",pY);
    }
  }
  FillHistogram("hPHOSClusterMult"  ,multPHOSClust[0]);
  FillHistogram("hPHOSClusterMultM1",multPHOSClust[1]);
  FillHistogram("hPHOSClusterMultM2",multPHOSClust[2]);
  FillHistogram("hPHOSClusterMultM3",multPHOSClust[3]);

  //Select photons for inv mass calculation
  Int_t inPHOS=0 ;
  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() || clu1->E()<0.3) continue;

    clu1->GetPosition(position);
    TVector3 global1(position) ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod1,cellX,cellZ) ) continue ;
    
    if (mod1 < 1 || mod1 > 3) {
      AliError(Form("Wrong module number %d",mod1));
      return;
    }

    //.................................Not required for pass4 AOD.................
    // Apply module misalignment
    /*
    Float_t dXmodule[3] = {-2.30, -2.11, -1.53}; // X-shift in local system for module 1,2,3
    Float_t dZmodule[3] = {-0.40, +0.52, +0.80}; // Z-shift in local system for module 1,2,3

    TVector3 globalXYZ(position[0],position[1],position[2]);
    TVector3 localXYZ;
    fPHOSGeo->Global2Local(localXYZ,globalXYZ,mod1) ;
    fPHOSGeo->Local2Global(mod1,localXYZ.X()+dXmodule[mod1-1],localXYZ.Z()+dZmodule[mod1-1],globalXYZ);
    for (Int_t ixyz=0; ixyz<3; ixyz++) position[ixyz]=globalXYZ[ixyz] ;
    clu1->SetPosition(position) ;
    */
    //..................................................

    clu1 ->GetMomentum(p1 ,vtx0);
    clu1 ->GetMomentum(pv1,vtxBest);
    
    p1 *= fRecalib[mod1-1];
    
    digMult   = clu1->GetNCells();
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(p1.X(),p1.Py(),p1.Z(),p1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod1) ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu1->GetNCells());
    ph->SetEMCx(global1.X());
    ph->SetEMCy(global1.Y());
    ph->SetEMCz(global1.Z());
    ph->SetDispBit(TestLambda(clu1->GetM20(),clu1->GetM02())) ;
    ph->SetCPVBit(clu1->GetEmcCpvDistance()>10.) ;
    ph->SetBC(TestBC(clu1->GetTOF()));
    
    inPHOS++ ;
  }

  // Fill Real disribution

  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Bool_t mainBC = (ph1->GetBC()==0 && ph2->GetBC()==0);
      Bool_t mainBC1= (ph1->GetBC()==0 || ph2->GetBC()==0);
      Bool_t diffBC = ((ph1->GetBC()==0 && ph2->GetBC()!=ph1->GetBC()) || 
		       (ph2->GetBC()==0 && ph2->GetBC()!=ph1->GetBC()));
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));
      Double_t ma12 = p12.M();
      Double_t pt12 = p12.Pt();

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) {
        FillHistogram("hMassPtA10",ma12 ,pt12 );
        FillHistogram("hMassPtvA10",pv12.M(),pv12.Pt());
        FillHistogram("hMassPtCA10",ma12 ,pt12, centr+0.5);
        FillHistogram("hMassSingle_all",ma12,ph1->Pt()) ;
        FillHistogram("hMassSingle_all",ma12,ph2->Pt()) ;
	if (mainBC) 
	  FillHistogram("hMassPtA10BC0",ma12 ,pt12 );
	if (diffBC) 
	  FillHistogram("hMassPtA10BC1",ma12 ,pt12 );
	if (mainBC1) 
	  FillHistogram("hMassPtA10BC2",ma12 ,pt12 );

	if(!eventVtxExist)
	  FillHistogram("hMassPtA10nvtx",ma12 ,pt12 );
	if(eventVtxExist)
	  FillHistogram("hMassPtA10vtx"  ,ma12 ,pt12 );
	if(eventVtxExist & eventVtxZ10cm)
	  FillHistogram("hMassPtA10vtx10",ma12 ,pt12 );
	if(eventV0AND)
	  FillHistogram("hMassPtA10V0AND",ma12 ,pt12 );
	if(eventPileup)
	  FillHistogram("hMassPtA10PU"   ,ma12 ,pt12 );

        if(ph1->IsCPVOK())
          FillHistogram("hMassSingle_cpv",ma12,ph1->Pt()) ;
        if(ph2->IsCPVOK())
          FillHistogram("hMassSingle_cpv",ma12,ph2->Pt()) ;
        if(ph1->IsDispOK())
          FillHistogram("hMassSingle_disp",ma12,ph1->Pt()) ;
        if(ph2->IsDispOK())
          FillHistogram("hMassSingle_disp",ma12,ph2->Pt()) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK())
          FillHistogram("hMassSingle_both",ma12,ph1->Pt()) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK())
          FillHistogram("hMassSingle_both",ma12,ph2->Pt()) ;
 

        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMassPtCA10_cpv",ma12 ,pt12, centr+0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK()){
          FillHistogram("hMassPtCA10_disp",ma12 ,pt12, centr+0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA10_both",ma12 ,pt12, centr+0.5);
        }
        if (asym<0.8) {
          FillHistogram("hMassPtA08",ma12,pt12);
        }
        if (asym<0.7) {
          FillHistogram("hMassPtA07",ma12,pt12);
	  FillHistogram("hMassPtvA07",pv12.M(),pv12.Pt());
          FillHistogram("hMassPtCA07",ma12 ,pt12, centr+0.5);
	  if(!eventVtxExist)
	    FillHistogram("hMassPtA07nvtx",ma12 ,pt12 );
	  if(eventVtxExist)
	    FillHistogram("hMassPtA07vtx"  ,ma12 ,pt12 );
	  if(eventVtxExist && eventV0AND)
	    FillHistogram("hMassPtA07V0AND",ma12 ,pt12 );
	  if(eventPileup)
	    FillHistogram("hMassPtA07PU"   ,ma12 ,pt12 );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA07_cpv",ma12 ,pt12, centr+0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK()){
            FillHistogram("hMassPtCA07_disp",ma12 ,pt12, centr+0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMassPtCA07_both",ma12 ,pt12, centr+0.5);
        }

        }
        if (asym<0.1) {
          FillHistogram("hMassPtA01",ma12,pt12);
        }
	if (TMath::Abs(ma12-0.135)<0.03)
	  FillHistogram("hAsymPtPi0",asym   ,pt12);
	if (TMath::Abs(ma12-0.547)<0.09)
	  FillHistogram("hAsymPtEta",asym   ,pt12);

        if (ph1->Module()==1 && ph2->Module()==1) {
	  FillHistogram("hMassPtM1",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M1",asym   ,pt12);
	}
        if (ph1->Module()==2 && ph2->Module()==2) {
	  FillHistogram("hMassPtM2",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M2",asym   ,pt12);
	}
        if (ph1->Module()==3 && ph2->Module()==3) {
	  FillHistogram("hMassPtM3",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M3",asym   ,pt12);
	}
        if ((ph1->Module()==1 && ph2->Module()==2) ||
	    (ph1->Module()==2 && ph2->Module()==1)) {
	  FillHistogram("hMassPtM12",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M12",asym   ,pt12);
	}
        if ((ph1->Module()==2 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==2)) {
	  FillHistogram("hMassPtM23",ma12 ,pt12 );
	  if (TMath::Abs(ma12-0.135)<0.03) FillHistogram("hAsymPtPi0M23",asym   ,pt12);
	}
        if ((ph1->Module()==1 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==1)) FillHistogram("hMassPtM13",ma12 ,pt12 );

	if ( TMath::Abs(ph1->EMCz()) < 20. || TMath::Abs(ph2->EMCz()) < 20.)
	  FillHistogram("hMassPt20cm",ma12 ,pt12 );
	if ((TMath::Abs(ph1->EMCz()) > 20. && TMath::Abs(ph1->EMCz()) < 40.) ||
	    (TMath::Abs(ph2->EMCz()) > 20. && TMath::Abs(ph2->EMCz()) < 40.))
	  FillHistogram("hMassPt40cm",ma12 ,pt12 );
	if ( TMath::Abs(ph1->EMCz()) > 40. || TMath::Abs(ph2->EMCz()) > 40.)
	  FillHistogram("hMassPt60cm",ma12 ,pt12 );

      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) {
        FillHistogram("hMassPtN3",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) {
        FillHistogram("hMassPtN4",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) {
        FillHistogram("hMassPtN5",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) {
        FillHistogram("hMassPtN6",ma12 ,pt12 );
      }

    } // end of loop i2
  } // end of loop i1
  
  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
      AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Bool_t mainBC = (ph1->GetBC()==0 && ph2->GetBC()==0);
      Bool_t mainBC1= (ph1->GetBC()==0 || ph2->GetBC()==0);
      Bool_t diffBC = ((ph1->GetBC()==0 && ph2->GetBC()!=ph1->GetBC()) || 
		       (ph2->GetBC()==0 && ph2->GetBC()!=ph1->GetBC()));
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));
      Double_t ma12 = p12.M();
      Double_t pt12 = p12.Pt();

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) {
        FillHistogram("hMiMassPtA10",ma12 ,pt12 );
        FillHistogram("hMiMassPtvA10",pv12.M(),pv12.Pt());
        FillHistogram("hMiMassPtCA10",ma12 ,pt12, centr+0.5);
        FillHistogram("hMiMassSingle_all",ma12,ph1->Pt()) ;
        FillHistogram("hMiMassSingle_all",ma12,ph2->Pt()) ;
	if (mainBC) 
	  FillHistogram("hMiMassPtA10BC0",ma12 ,pt12 );
	if (diffBC) 
	  FillHistogram("hMiMassPtA10BC1",ma12 ,pt12 );
	if (mainBC1) 
	  FillHistogram("hMiMassPtA10BC2",ma12 ,pt12 );

	if(!eventVtxExist)
	  FillHistogram("hMiMassPtA10nvtx",ma12 ,pt12 );
	if(eventVtxExist)
	  FillHistogram("hMiMassPtA10vtx"  ,ma12 ,pt12 );
	if(eventVtxExist & eventVtxZ10cm)
	  FillHistogram("hMiMassPtA10vtx10",ma12 ,pt12 );
	if(eventV0AND)
	  FillHistogram("hMiMassPtA10V0AND",ma12 ,pt12 );
	if(eventPileup)
	  FillHistogram("hMiMassPtA10PU"   ,ma12 ,pt12 );

        if(ph1->IsCPVOK())
          FillHistogram("hMiMassSingle_cpv",ma12,ph1->Pt()) ;
        if(ph2->IsCPVOK())
          FillHistogram("hMiMassSingle_cpv",ma12,ph2->Pt()) ;
        if(ph1->IsDispOK())
          FillHistogram("hMiMassSingle_disp",ma12,ph1->Pt()) ;
        if(ph2->IsDispOK())
          FillHistogram("hMiMassSingle_disp",ma12,ph2->Pt()) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK())
          FillHistogram("hMiMassSingle_both",ma12,ph1->Pt()) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK())
          FillHistogram("hMiMassSingle_both",ma12,ph2->Pt()) ;


        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMiMassPtCA10_cpv",ma12 ,pt12, centr+0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK()){
          FillHistogram("hMiMassPtCA10_disp",ma12 ,pt12, centr+0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA10_both",ma12 ,pt12, centr+0.5);
        }
        if (asym<0.8) {
          FillHistogram("hMiMassPtA08",ma12,pt12);
        }
        if (asym<0.7) {
          FillHistogram("hMiMassPtA07",ma12,pt12);
 	  FillHistogram("hMiMassPtvA07",pv12.M(),pv12.Pt());
	  FillHistogram("hMiMassPtCA07",ma12 ,pt12, centr+0.5);
	  if(!eventVtxExist)
	    FillHistogram("hMiMassPtA07nvtx",ma12 ,pt12 );
	  if(eventVtxExist)
	    FillHistogram("hMiMassPtA07vtx"  ,ma12 ,pt12 );
	  if(eventVtxExist && eventV0AND)
	    FillHistogram("hMiMassPtA07V0AND",ma12 ,pt12 );
	  if(eventPileup)
	    FillHistogram("hMiMassPtA07PU"   ,ma12 ,pt12 );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA07_cpv",ma12 ,pt12, centr+0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK()){
            FillHistogram("hMiMassPtCA07_disp",ma12 ,pt12, centr+0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMiMassPtCA07_both",ma12 ,pt12, centr+0.5);
          }
        }
        if (asym<0.1) {
          FillHistogram("hMiMassPtA01",ma12,pt12);
        }
        FillHistogram("hMiAsymPt",asym   ,pt12);

        if (ph1->Module()==1 && ph2->Module()==1) FillHistogram("hMiMassPtM1",ma12 ,pt12 );
        if (ph1->Module()==2 && ph2->Module()==2) FillHistogram("hMiMassPtM2",ma12 ,pt12 );
        if (ph1->Module()==3 && ph2->Module()==3) FillHistogram("hMiMassPtM3",ma12 ,pt12 );
        if ((ph1->Module()==1 && ph2->Module()==2) ||
	    (ph1->Module()==2 && ph2->Module()==1)) FillHistogram("hMiMassPtM12",ma12 ,pt12 );
        if ((ph1->Module()==2 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==2)) FillHistogram("hMiMassPtM23",ma12 ,pt12 );
        if ((ph1->Module()==1 && ph2->Module()==3) ||
	    (ph1->Module()==3 && ph2->Module()==1)) FillHistogram("hMiMassPtM13",ma12 ,pt12 );

	if (TMath::Abs(ph1->EMCz()) < 20. || TMath::Abs(ph2->EMCz()) < 20.)
	  FillHistogram("hMiMassPt20cm",ma12 ,pt12 );
	if ((TMath::Abs(ph1->EMCz()) > 20. && TMath::Abs(ph1->EMCz()) < 40.) ||
	    (TMath::Abs(ph2->EMCz()) > 20. && TMath::Abs(ph2->EMCz()) < 40.))
	  FillHistogram("hMiMassPt40cm",ma12 ,pt12 );
	if (TMath::Abs(ph1->EMCz()) > 40. || TMath::Abs(ph2->EMCz()) > 40.)
	  FillHistogram("hMiMassPt60cm",ma12 ,pt12 );

      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) {
        FillHistogram("hMiMassPtN3",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) {
        FillHistogram("hMiMassPtN4",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) {
        FillHistogram("hMiMassPtN5",ma12 ,pt12 );
      }
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) {
        FillHistogram("hMiMassPtN6",ma12 ,pt12 );
      }
      
      } // end of loop i2
    }
  } // end of loop i1
 
  
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(prevPHOS->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }
  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisTaskPi0PP::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  /*
  Int_t nPP = fnCINT1B - fnCINT1A - fnCINT1C + 2*fnCINT1E;
  FillHistogram("hTrigClass",1,fnCINT1B);
  FillHistogram("hTrigClass",2,fnCINT1A);
  FillHistogram("hTrigClass",3,fnCINT1C);
  FillHistogram("hTrigClass",4,fnCINT1E);
  FillHistogram("hTrigClass",5,nPP);
  Printf("fnCINT1B=%d, fnCINT1A=%d ,fnCINT1C=%d ,fnCINT1E=%d, nPP=%d",
     fnCINT1B,fnCINT1A,fnCINT1C,fnCINT1E,nPP);
  */
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0PP::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d ",mod)) ;
      return kTRUE ;
    }
    if(!fPHOSBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d",mod)) ;
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
  return kTRUE ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0PP::FillHistogram(const char * key,Double_t x)const
{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0PP::FillHistogram(const char * key,Double_t x,Double_t y)const
{
  //FillHistogram
  TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
  if(th1)
    th1->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0PP::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const
{
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

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0PP::TestLambda(Double_t l1,Double_t l2){
  Double_t l1Mean=1.22 ;
  Double_t l2Mean=2.0 ;
  Double_t l1Sigma=0.42 ;
  Double_t l2Sigma=0.71 ;
  Double_t c=-0.59 ;
  Double_t R2=(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma+(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma-c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<9.) ;
}
//_____________________________________________________________________________
Int_t AliAnalysisTaskPi0PP::TestBC(Double_t tof){
  Int_t bc = (Int_t)(TMath::Ceil((tof + fBCgap/2)/fBCgap) - 1);
  return bc;
}
