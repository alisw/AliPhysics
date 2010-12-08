#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TList.h"
#include "TCanvas.h"

//#include <iostream.h>
//#include <TVector3.h>
//#include <TGeoHMatrix.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
//#include "AliESDCaloCluster.h"
//#include "AliESDCaloCells.h"
//#include "AliESDHeader.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
#include "AliESDZDC.h"

#include "AliAnalysisTaskPt.h"
//#include "AliEMCALGeoUtils.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) 
  : AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fMyTr(0), jev(0), iev(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPt::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  // fEMCALGeo = new AliEMCALGeoUtils("EMCAL_FIRSTYEAR");

  // float pi=acos(-1.0);
  fOutputList = new TList();
/*
  fHstPt = new TH1F("fHstPt", "Pt distribution", 100, 0.0, 20.0);
  fEMCe = new TH1F("fEMCe", "E distribution", 100, 0.0, 10.0);
  fEMCt = new TH1F("fEMCt", "T distribution", 100, 0.0, 1000.0);
  fEMCn = new TH1F("fEMCn", "N distribution", 100, -0.5, 99.5);
  fEMCm = new TH2F("fEMCm", "phi-zps map", 50, 0, pi, 50, -300.0, 300.0);
  fCelle = new TH1F("fCelle", "Cell E distribution", 100, 0.0, 10.0);
  fCellf = new TH1F("fCellf", "Cell E fraction", 100, -0.1, 1.1);
  fCellt = new TH1F("fCellt", "Cell T distribution", 100, 0.0, 1000.0);
  fCellc = new TH2F("fCellc", "T-ID map", 5000, -200.0, 4800.0, 50, 0.0, 1000.0);
  fCelld = new TH2F("fCelld", "E-ID map", 5000, -200.0, 4800.0, 50, 0.0, 10.0);
  fCellm = new TH2F("fCellm", "phi-eta map", 50, 0, pi, 50, -300.0, 300.0);
  
  fOutputList->Add(fHstPt);
  fOutputList->Add(fEMCe);
  fOutputList->Add(fEMCt);
  fOutputList->Add(fEMCn);
  fOutputList->Add(fEMCm);
  fOutputList->Add(fCelle);
  fOutputList->Add(fCellf);
  fOutputList->Add(fCellt);
  fOutputList->Add(fCellc);
  fOutputList->Add(fCelld);
  fOutputList->Add(fCellm);
  iii=0;
*/
  jev=0;
  iev=0;
  fMyTr = new TTree("fMyTr","beam tree");
  fOutputList->Add(fMyTr);
}

//________________________________________________________________________
void AliAnalysisTaskPt::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

//AliESDHeader *header = fESD->GetHeader();
//fEMCALGeo->SetMisalMatrixes(header->GetEMCALMisalMatrix());
//
//for(Int_t mod=0; mod<(fEMCALGeo->GetEMCGeometry())->GetNumberOfSuperModules(); mod++){
//  if(fESD){
//    const TGeoHMatrix* m=fESD->GetEMCALMatrix(mod) ;
//    fEMCALGeo->SetMisalMatrix(m, mod) ;
//  }
//}
//

  int          ntr, etr;
  float        z1n, z2n, z1p, z2p, zvt;
  float        z1s, z2s, w1s, w2s, t1s, t2s, v1s, v2s;
  float        z1nt[10], z2nt[10], z1pt[10], z2pt[10];
  float        t0amp[24], t0tim[24], v0mul[64];
  ntr = (int)(fESD->GetNumberOfTracks());
  if (ntr>0 && iev==0) {
   fMyTr->Branch("jev",&jev,"jev/I");
   fMyTr->Branch("iev",&iev,"iev/I");
   fMyTr->Branch("ntr",&ntr,"ntr/I");
   fMyTr->Branch("etr",&etr,"etr/I");
   fMyTr->Branch("z1n",&z1n,"z1n/F");
   fMyTr->Branch("z2n",&z2n,"z2n/F");
   fMyTr->Branch("z1p",&z1p,"z1p/F");
   fMyTr->Branch("z2p",&z2p,"z2p/F");
   fMyTr->Branch("zvt",&zvt,"zvt/F");
   fMyTr->Branch("z1s",&z1s,"z1s/F");
   fMyTr->Branch("z2s",&z2s,"z2s/F");
   fMyTr->Branch("w1s",&w1s,"w1s/F");
   fMyTr->Branch("w2s",&w2s,"w2s/F");
   fMyTr->Branch("t1s",&t1s,"t1s/F");
   fMyTr->Branch("t2s",&t2s,"t2s/F");
   fMyTr->Branch("v1s",&v1s,"v1s/F");
   fMyTr->Branch("v2s",&v2s,"v2s/F");
   fMyTr->Branch("z1nt",z1nt,"z1nt[10]/F");
   fMyTr->Branch("z2nt",z2nt,"z2nt[10]/F");
   fMyTr->Branch("z1pt",z1pt,"z1pt[10]/F");
   fMyTr->Branch("z2pt",z2pt,"z2pt[10]/F");
   fMyTr->Branch("t0amp",t0amp,"t0amp[24]/F");
   fMyTr->Branch("t0tim",t0tim,"t0tim[24]/F");
   fMyTr->Branch("v0mul",v0mul,"v0mul[64]/F");
   printf("my tree branches are created\n");
  }
  const AliESDZDC* zdcData = fESD->GetESDZDC();
  const AliESDTZERO* tzrData = fESD->GetESDTZERO();
  const AliESDVZERO* vzrData = fESD->GetVZEROData();
//const AliESDVZERO* vzrData = (const_cast<AliESDEvent*>(fESD))->GetVZEROData();
  const AliESDVertex* vtxData = fESD->GetPrimaryVertex();

  etr = (int)(fESD->GetNumberOfCaloClusters());
  zvt = (float)(vtxData->GetZ());
  z1n = (float)(zdcData->GetZDCN1Energy());
  z2n = (float)(zdcData->GetZDCN2Energy());
  z1p = (float)(zdcData->GetZDCP1Energy());
  z2p = (float)(zdcData->GetZDCP2Energy());
  const Double_t *z1ntTmp,*z2ntTmp,*z1ptTmp,*z2ptTmp;
  const Double_t *z1ntTmq,*z2ntTmq,*z1ptTmq,*z2ptTmq;
  const Double32_t *t0timTmp,*t0ampTmp;
  z1ntTmp = zdcData->GetZN1TowerEnergyLR();
  z2ntTmp = zdcData->GetZN2TowerEnergyLR();
  z1ptTmp = zdcData->GetZP1TowerEnergyLR();
  z2ptTmp = zdcData->GetZP2TowerEnergyLR();
  z1ntTmq = zdcData->GetZN1TowerEnergy();
  z2ntTmq = zdcData->GetZN2TowerEnergy();
  z1ptTmq = zdcData->GetZP1TowerEnergy();
  z2ptTmq = zdcData->GetZP2TowerEnergy();
  t0timTmp = tzrData->GetT0time();
  t0ampTmp = tzrData->GetT0amplitude();
  for (int i=0; i<5; i++) {
   z1nt[i]=(float)(z1ntTmp[i]);  z1nt[i+5]=(float)(z1ntTmq[i]);
   z2nt[i]=(float)(z2ntTmp[i]);  z2nt[i+5]=(float)(z2ntTmq[i]);
   z1pt[i]=(float)(z1ptTmp[i]);  z1pt[i+5]=(float)(z1ptTmq[i]);
   z2pt[i]=(float)(z2ptTmp[i]);  z2pt[i+5]=(float)(z2ptTmq[i]);
  }
  for (int i=0; i<24; i++) {
   t0tim[i] = (float)(t0timTmp[i]);
   t0amp[i] = (float)(t0ampTmp[i]);
  }
  for (int i=0; i<64; i++) {
   v0mul[i] = (float)(vzrData->GetMultiplicity(i));
  }
  z1s=z1nt[1]+z1nt[2]+z1nt[3]+z1nt[4];
  z2s=z2nt[1]+z2nt[2]+z2nt[3]+z2nt[4];
  w1s=z1pt[1]+z1pt[2]+z1pt[3]+z1pt[4];
  w2s=z2pt[1]+z2pt[2]+z2pt[3]+z2pt[4];
  t1s=0; t2s=0; v1s=0; v2s=0;
  for (Int_t n=0; n<24; n++) {
   if (t0tim[n]>100 && t0amp[n]>0) {
    if (n<12) t1s+=t0amp[n];
    else      t2s+=t0amp[n];
   }
  }
  for (Int_t n=0; n<64; n++) {
   if (v0mul[n]>0) {
    if (n<32) v1s+=v0mul[n];
    else      v2s+=v0mul[n];
   }
  }
  jev++;
  if (ntr>0) {
   fMyTr->Fill();
   iev++;
   printf(": %d %d %d %d %f : %f %f %f %f :\n",jev,iev,ntr,etr,zvt,z1s+z2s,w1s+w2s,t1s+t2s,v1s+v2s);
  }

/*
  if(Ntrack==0 || nCluster==0) {
    printf("no trk or clus %d %d\n",Ntrack,nCluster);
    return;
  }

  Double_t bfield = fESD->GetMagneticField();
  printf("bfield %f\n",bfield);

  // Track loop to fill a pT spectrum
  Double_t posI[5],posO[5],posE[5];
  Double_t minR=300.0;
  Double_t maxR=500.0;
  Double_t stpR=10.0;
  Float_t arr[5][22][1000]; 
  int k=0;
  for (Int_t iTracks = 0; iTracks < fESD->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = fESD->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    float pt=track->Pt();
    fHstPt->Fill(pt);
    float px=track->Px();
    float py=track->Py();
    float ph=atan2(py,px);
    track->GetInnerXYZ(posI);
    posI[3] = sqrt(posI[0]*posI[0]+posI[1]*posI[1]);
    posI[4] = atan2(posI[1],posI[0]);
    track->GetOuterXYZ(posO);
    posO[3] = sqrt(posO[0]*posO[0]+posO[1]*posO[1]);
    posO[4] = atan2(posO[1],posO[0]);
    if (pt>1.0 && ph>1.0 && ph<2.5 &&
        posI[3]>80.0  && posI[3]<90.0 && 
        posO[3]>280.0 && posO[3]<300.0) {
      for (int i=0; i<5; i++) arr[i][0][k]=posI[i];
      for (int i=0; i<5; i++) arr[i][1][k]=posO[i];
      int j=0;
      for (Double_t rad=minR; rad<maxR; rad+=stpR) {
        track->GetOuterParam()->GetXYZAt(rad,bfield,posE);
        posE[3] = sqrt(posE[0]*posE[0]+posE[1]*posE[1]);
        posE[4] = atan2(posE[1],posE[0]);
        for (int i=0; i<5; i++) arr[i][2+j][k]=posE[i];
        j++;
      }
      if (k<1000) k++;
      else printf("error max trk 1000\n");
    }
  } //track loop 
  int nTrkEMC=k;

  // EMC loop
    float arr0[10000][7],arr1[10000][7];
    float scl=1.0E+09;
//  TVector3 CaloCellPos; 
    AliESDCaloCells &cells= *(fESD->GetEMCALCells());
    int totalCell = cells.GetNumberOfCells();
    int totCell = 0;
    int totClus = 0;
    for (Int_t iCluster=0; iCluster<fESD->GetNumberOfCaloClusters(); iCluster++) {
       AliESDCaloCluster * clust = fESD->GetCaloCluster(iCluster);
       if (clust->IsEMCAL())
          {
            Int_t nCells = clust->GetNCells();
	    UShort_t * index = clust->GetCellsAbsId() ;
	    Double_t * fraction = clust->GetCellsAmplitudeFraction() ;
            Float_t pos[3]={0,0,0};
            clust->GetPosition(pos);
            float rad=sqrt(pos[0]*pos[0]+pos[1]*pos[1]);
            float phi=atan2(pos[1],pos[0]);
            float zps=pos[2];
            float tof = scl*clust->GetTOF();
            float ene = clust->E();
            fEMCe->Fill(ene);
            fEMCt->Fill(tof);
            fEMCn->Fill((float)nCells);
            fEMCm->Fill(phi,zps);
            arr0[totClus][0]=ene;
            arr0[totClus][1]=tof;
            arr0[totClus][2]=rad;
            arr0[totClus][3]=phi;
            arr0[totClus][4]=zps;
            arr0[totClus][5]=totCell;
            arr0[totClus][6]=totCell+nCells;
            if (totClus<10000) totClus++;
            else printf("error max cls 10000\n");
	    for(Int_t i = 0; i < nCells ; i++)
                 {
	       	   Int_t absId =   index[i]; // or clus->GetCellNumber(i) ;
                   Double_t cpos[3]={0,0,0};
//
                   if (fEMCALGeo) { 
//                   fEMCALGeo->GetGlobal(absId,CaloCellPos); 
                     fEMCALGeo->GetGlobal(absId,cpos); 
                   } else { 
                     printf("Error, fEMCALGeo=0\n");
                   }
//                 Double_t ccphi=CaloCellPos.Phi();
//                 Double_t cceta=CaloCellPos.Eta();
//
                   Double_t ccrad=sqrt(cpos[0]*cpos[0]+cpos[1]*cpos[1]);
                   Double_t ccphi=atan2(cpos[1],cpos[0]);
                   Double_t cczps=cpos[2];

		   Double_t ampFract =  fraction[i];
		   Float_t amp       = cells.GetCellAmplitude(absId) ;
		   Float_t time      = scl*cells.GetCellTime(absId);
                   fCelle->Fill(amp);
                   fCellf->Fill(ampFract);
                   fCellt->Fill(time);
                   fCellc->Fill((float)absId,time);
                   fCelld->Fill((float)absId,amp);
                   fCellm->Fill(ccphi,cczps);
                   arr1[totCell][0]=i;
                   arr1[totCell][1]=absId;
                   arr1[totCell][2]=amp;
                   arr1[totCell][3]=time;
                   arr1[totCell][4]=ccrad;
                   arr1[totCell][5]=ccphi;
                   arr1[totCell][6]=cczps;
                   if (totCell<10000) totCell++;
                   else printf("error max cel 10000\n");
               }  // cell loop

          } // end of cluster loop

//      Printf("+++++ # EMC cluster E = %f", mesE);

      }
  if(nTrkEMC==0 || totClus==0) {
    printf("no emc trk or emc clus %d %d\n",nTrkEMC,totClus);
    return;
  }
  FILE *ofs;
  char opt[4];
  if (iii==0) sprintf(opt,"w");
  else        sprintf(opt,"a");
  if ((ofs = fopen("mydata.dat",opt)) == NULL) {
    printf("error in fopen\n");
  }
  iii++;
  printf("%d %s : %d %d : %d %d : %d %d\n", iii, opt,
  Ntrack, nTrkEMC, nCluster, totClus, totalCell, totCell);
//fprintf(ofs,"%d %d %d %d %d %d\n", 
//Ntrack, nTrkEMC, nCluster, totClus, totalCell, totCell);
  for (int itrk=0; itrk<nTrkEMC; itrk++) {
    for (int ihit=0; ihit<22; ihit++) {
//    fprintf(ofs,"%f %f %f\n",
//    arr[3][ihit][itrk],arr[4][ihit][itrk],arr[2][ihit][itrk]);
    }
  }
  for (int icls=0; icls<totClus; icls++) {
//  for (int ihit=0; ihit<5; ihit++) fprintf(ofs,"%f ",arr0[icls][ihit]);
//  for (int ihit=5; ihit<7; ihit++) fprintf(ofs,"%d ",(int)arr0[icls][ihit]);
//  fprintf(ofs,"\n");
  }
  for (int icel=0; icel<totCell; icel++) {
//  for (int ihit=0; ihit<2; ihit++) fprintf(ofs,"%d ",(int)arr1[icel][ihit]);
//  for (int ihit=2; ihit<7; ihit++) fprintf(ofs,"%f ",arr1[icel][ihit]);
//  fprintf(ofs,"\n");
  }
  fclose(ofs);
*/
  // Post output data.
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

/*
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
  fHstPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHstPt) {
    printf("ERROR: fHstPt not available\n");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHstPt->DrawCopy("E");
*/

}
