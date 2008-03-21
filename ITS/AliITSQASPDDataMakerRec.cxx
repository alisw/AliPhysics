/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
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
/* $Id:$  */
//  *************************************************************
//  Checks the quality assurance 
//  by comparing with reference data
//  contained in a DB
//  -------------------------------------------------------------
//  W. Ferrarese + P. Cerello Feb 2008
//  INFN Torino
//  M. Nicassio D. Elia INFN Bari March 2008
//  maria.nicassio@ba.infn.it
    

// --- ROOT system ---
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
// --- Standard library ---

// --- AliRoot header files ---
#include "AliITSQADataMakerRec.h"
#include "AliITSQASPDDataMakerRec.h"
#include "AliLog.h"
#include "AliQA.h"
#include "AliRawReader.h"
#include "AliITSRecPoint.h"
ClassImp(AliITSQASPDDataMakerRec)

//____________________________________________________________________________
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(AliITSQADataMakerRec *aliITSQADataMakerRec, Bool_t kMode, Short_t ldc) :
TObject(),
fAliITSQADataMakerRec(aliITSQADataMakerRec),
fkOnline(kMode),
fLDC(ldc),
fSPDhRaws(0),
fSPDhRecs(0),
fRawsOffset(0),
fRecsOffset(0)
{
  //ctor used to discriminate OnLine-Offline analysis  
}

//____________________________________________________________________________ 
AliITSQASPDDataMakerRec::AliITSQASPDDataMakerRec(const AliITSQASPDDataMakerRec& qadm) :
TObject(),
fAliITSQADataMakerRec(qadm.fAliITSQADataMakerRec),
fkOnline(qadm.fkOnline),
fLDC(qadm.fLDC),
fSPDhRaws(qadm.fSPDhRaws),
fSPDhRecs(qadm.fSPDhRecs),
fRawsOffset(qadm.fRawsOffset),
fRecsOffset(qadm.fRecsOffset)
{
  //copy ctor 
  fAliITSQADataMakerRec->SetName((const char*)qadm.fAliITSQADataMakerRec->GetName()) ; 
  fAliITSQADataMakerRec->SetTitle((const char*)qadm.fAliITSQADataMakerRec->GetTitle());
  }

//__________________________________________________________________
AliITSQASPDDataMakerRec::~AliITSQASPDDataMakerRec(){
  // destructor

}
//__________________________________________________________________

AliITSQASPDDataMakerRec& AliITSQASPDDataMakerRec::operator = (const AliITSQASPDDataMakerRec& qac )
{
  // Equal operator.
  this->~AliITSQASPDDataMakerRec();
  new(this) AliITSQASPDDataMakerRec(qac);
  return *this;
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::StartOfDetectorCycle()
{
  //Detector specific actions at start of cycle
  AliDebug(1,"AliITSQADM::Start of SPD Cycle\n");
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::EndOfDetectorCycle(AliQA::TASKINDEX_t /*task*/, TObjArray* /*list*/)
{
  // launch the QA checking
  AliDebug(1,"AliITSDM instantiates checker with Run(AliQA::kITS, task, list)\n"); 
  
  //AliQAChecker::Instance()->Run( AliQA::kITS , task, list);
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::InitRaws()
{ 
  // Initialization for RAW data - SPD -
  fRawsOffset = (fAliITSQADataMakerRec->fRawsQAList)->GetEntries();

  // custom code here

  //fSPDhRaws must be incremented by one unit every time a histogram is ADDED to the QA List

  AliDebug(1,Form("%d SPD Raws histograms booked\n",fSPDhRaws));

}


//____________________________________________________________________________
void AliITSQASPDDataMakerRec::MakeRaws(AliRawReader* /*rawReader*/)
{ 
  // Fill QA for RAW - SPD -
}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::InitRecPoints()
{
  // Initialization for RECPOINTS - SPD -
  fRecsOffset = (fAliITSQADataMakerRec->fRecPointsQAList)->GetEntries();


  TH1F* hlayer= new TH1F("LayPattern_SPD","Layer map - SPD",6,0.,6.);
  hlayer->GetXaxis()->SetTitle("Layer number");
  hlayer->GetYaxis()->SetTitle("Entries");
  fAliITSQADataMakerRec->Add2RecPointsList(hlayer, fSPDhRecs+fRecsOffset); 
  fSPDhRecs++;

  TH1F** hmod = new TH1F*[2];
  TH1F** hxl  = new TH1F*[2];
  TH1F** hzl  = new TH1F*[2];
  TH1F** hxg  = new TH1F*[2];
  TH1F** hyg  = new TH1F*[2];
  TH1F** hzg  = new TH1F*[2];
  TH1F** hr   = new TH1F*[2];
  TH1F** hphi = new TH1F*[2];
  TH1F** hMultSPDcl = new TH1F*[2];
  TH2F** hNyNz = new TH2F*[2];  // y and z cluster length
  TH2F** hPhiZ = new TH2F*[2];

  Float_t xlim[2]={4.5,8.};
  Float_t zlim[2]={15.,15.};

  Char_t name[50];
  Char_t title[50];
  for (Int_t iLay=0;iLay<2;iLay++) {
    sprintf(name,"ModPattern_SPD%d",iLay+1);
    sprintf(title,"Module map - SPD Layer %d",iLay+1);
    hmod[iLay]=new TH1F(name,title,fgknSPDmodules,0,fgknSPDmodules);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hmod[iLay], fSPDhRecs +fRecsOffset); 
    fSPDhRecs++;

    sprintf(name,"xLoc_SPD%d",iLay+1);
    sprintf(title,"Local x coordinate - SPD Layer %d",iLay+1);
    hxl[iLay]=new TH1F(name,title,100,-4.,4.);
    hxl[iLay]->GetXaxis()->SetTitle("Local x [cm]");
    hxl[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hxl[iLay], fSPDhRecs +fRecsOffset);
    fSPDhRecs++;

    sprintf(name,"zLoc_SPD%d",iLay+1);
    sprintf(title,"Local z coordinate - SPD Layer %d",iLay+1);
    hzl[iLay]=new TH1F(name,title,100,-4.,4.);
    hzl[iLay]->GetXaxis()->SetTitle("Local z [cm]");
    hzl[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hzl[iLay], fSPDhRecs+fRecsOffset); 
    fSPDhRecs++;

    sprintf(name,"xGlob_SPD%d",iLay+1);
    sprintf(title,"Global x coordinate - SPD Layer %d",iLay+1);
    hxg[iLay]=new TH1F(name,title,100,-xlim[iLay],xlim[iLay]);
    hxg[iLay]->GetXaxis()->SetTitle("Global x [cm]");
    hxg[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hxg[iLay],fSPDhRecs+fRecsOffset);  
    fSPDhRecs++;

    sprintf(name,"yGlob_SPD%d",iLay+1);
    sprintf(title,"Global y coordinate - SPD Layer %d",iLay+1);
    hyg[iLay]=new TH1F(name,title,100,-xlim[iLay],xlim[iLay]);
    hyg[iLay]->GetXaxis()->SetTitle("Global y [cm]");
    hyg[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hyg[iLay], fSPDhRecs+fRecsOffset); 
    fSPDhRecs++;

    sprintf(name,"zGlob_SPD%d",iLay+1);
    sprintf(title,"Global z coordinate - SPD Layer %d",iLay+1);
    hzg[iLay]=new TH1F(name,title,150,-zlim[iLay],zlim[iLay]);
    hzg[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hzg[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hzg[iLay], fSPDhRecs+fRecsOffset); 
    fSPDhRecs++;

    sprintf(name,"r_SPD%d",iLay+1);
    sprintf(title,"Radius - SPD Layer %d",iLay+1);
    hr[iLay]=new TH1F(name,title,100,0.,10.);
    hr[iLay]->GetXaxis()->SetTitle("r [cm]");
    hr[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hr[iLay], fSPDhRecs+fRecsOffset);  
    fSPDhRecs++;

    sprintf(name,"phi_SPD%d",iLay+1);
    sprintf(title,"#varphi - SPD Layer %d",iLay+1);
    hphi[iLay]=new TH1F(name,title,600,0.,2*TMath::Pi());
    hphi[iLay]->GetXaxis()->SetTitle("#varphi [rad]");
    hphi[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hphi[iLay], fSPDhRecs+fRecsOffset);
    fSPDhRecs++;
    
    sprintf(name,"SizeYvsZ_SPD%d",iLay+1);
    sprintf(title,"Cluster dimension - SPD Layer %d",iLay+1);
    hNyNz[iLay]=new TH2F(name,title,100,0.,100.,100,0.,100.);
    hNyNz[iLay]->GetXaxis()->SetTitle("z length");
    hNyNz[iLay]->GetYaxis()->SetTitle("y length");
    fAliITSQADataMakerRec->Add2RecPointsList(hNyNz[iLay], fSPDhRecs+fRecsOffset); 
    fSPDhRecs++;

    sprintf(name,"phi_z_SPD%d",iLay+1);
    sprintf(title,"#varphi vs z - SPD Layer %d",iLay+1);
    hPhiZ[iLay]=new TH2F(name,title,150,-zlim[iLay],zlim[iLay],100,0.,2*TMath::Pi());
    hPhiZ[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hPhiZ[iLay]->GetYaxis()->SetTitle("#varphi [rad]");
    fAliITSQADataMakerRec->Add2RecPointsList(hPhiZ[iLay], fSPDhRecs+fRecsOffset);
    fSPDhRecs++;

  }

  TH2F *hrPhi=new TH2F("r_phi_SPD","#varphi vs r - SPD",100,0.,10.,100,0.,2*TMath::Pi());
  hrPhi->GetXaxis()->SetTitle("r [cm]");
  hrPhi->GetYaxis()->SetTitle("#varphi [rad]");
  fAliITSQADataMakerRec->Add2RecPointsList(hrPhi, fSPDhRecs+fRecsOffset);
  fSPDhRecs++;

  TH2F *hxy=new TH2F("x_y_SPD","Global y vs x - SPD",200,-10.,10.,200,-10.,10.);
  hxy->GetXaxis()->SetTitle("Global x [cm]");
  hxy->GetYaxis()->SetTitle("Global y [cm]");
  fAliITSQADataMakerRec->Add2RecPointsList(hxy, fSPDhRecs+fRecsOffset);
  fSPDhRecs++;

  for (Int_t iLay=0;iLay<2;iLay++) {
    sprintf(name,"Multiplicity_SPD%d",iLay+1);
    sprintf(title,"Cluster multiplicity - SPD Layer %d",iLay+1);
    hMultSPDcl[iLay]=new TH1F(name,title,200,0.,200.);
    hMultSPDcl[iLay]->GetXaxis()->SetTitle("Cluster multiplicity");
    hMultSPDcl[iLay]->GetYaxis()->SetTitle("Entries");
    fAliITSQADataMakerRec->Add2RecPointsList(hMultSPDcl[iLay], fSPDhRecs+fRecsOffset);
    fSPDhRecs++;
  } 

  TH2F *hMultSPDcl2MultSPDcl1 =
            new TH2F("MultCorrelation_SPD","Cluster multiplicity correlation - SPD",200,0.,200.,200,0.,200.);
  hMultSPDcl2MultSPDcl1->GetXaxis()->SetTitle("Clusters multiplicity (Layer 1)");
  hMultSPDcl2MultSPDcl1->GetYaxis()->SetTitle("Clusters multiplicity (Layer 2)"); 
  fAliITSQADataMakerRec->Add2RecPointsList(hMultSPDcl2MultSPDcl1, fSPDhRecs+fRecsOffset);
  fSPDhRecs++;

  AliDebug(1,Form("%d SPD Recs histograms booked\n",fSPDhRecs));

}

//____________________________________________________________________________ 
void AliITSQASPDDataMakerRec::MakeRecPoints(TTree * clusterTree)
{
  // Fill QA for RecPoints - SPD -
      TClonesArray* itsClusters = new TClonesArray("AliITSRecPoint");
      TBranch* itsClusterBranch=clusterTree->GetBranch("ITSRecPoints");
      if (!itsClusterBranch) {
        AliError("can't get the branch with the ITS clusters !");
        return;
      }
      itsClusterBranch->SetAddress(&itsClusters);
      Int_t nItsMods = (Int_t)clusterTree->GetEntries();

      Float_t cluGlo[3] = {0.,0.,0.};
      Int_t nClusters[2] = {0,0};

      for (Int_t iIts=0; iIts < nItsMods; iIts++) {

      if (!clusterTree->GetEvent(iIts))    continue;
      Int_t nCluster = itsClusters->GetEntriesFast();
      // loop over clusters
      while(nCluster--) {
        AliITSRecPoint* cluster = (AliITSRecPoint*)itsClusters->UncheckedAt(nCluster);

        if (cluster->GetLayer()>1)        continue;
        Int_t lay=cluster->GetLayer();
        fAliITSQADataMakerRec->GetRecPointsData(0 +fRecsOffset)->Fill(lay);
        cluster->GetGlobalXYZ(cluGlo);
        Float_t rad=TMath::Sqrt(cluGlo[0]*cluGlo[0]+cluGlo[1]*cluGlo[1]);
        Float_t phi= TMath::Pi() + TMath::ATan2(-cluGlo[1],-cluGlo[0]);
        if (lay==0) {
	  fAliITSQADataMakerRec->GetRecPointsData(1 +fRecsOffset)->Fill(iIts);
	  fAliITSQADataMakerRec->GetRecPointsData(2 +fRecsOffset)->Fill(cluster->GetDetLocalX());
	  fAliITSQADataMakerRec->GetRecPointsData(3 +fRecsOffset)->Fill(cluster->GetDetLocalZ());
	  fAliITSQADataMakerRec->GetRecPointsData(4 +fRecsOffset)->Fill(cluGlo[0]);
	  fAliITSQADataMakerRec->GetRecPointsData(5 +fRecsOffset)->Fill(cluGlo[1]);
	  fAliITSQADataMakerRec->GetRecPointsData(6 +fRecsOffset)->Fill(cluGlo[2]);
	  fAliITSQADataMakerRec->GetRecPointsData(7 +fRecsOffset)->Fill(rad);
	  fAliITSQADataMakerRec->GetRecPointsData(8 +fRecsOffset)->Fill(phi);
	  fAliITSQADataMakerRec->GetRecPointsData(9 +fRecsOffset)->Fill(cluster->GetNz(),cluster->GetNy());
	  fAliITSQADataMakerRec->GetRecPointsData(10 +fRecsOffset)->Fill(cluGlo[2],phi);
 	} else  {
          fAliITSQADataMakerRec->GetRecPointsData(11 +fRecsOffset)->Fill(iIts);
          fAliITSQADataMakerRec->GetRecPointsData(12 +fRecsOffset)->Fill(cluster->GetDetLocalX());
          fAliITSQADataMakerRec->GetRecPointsData(13 +fRecsOffset)->Fill(cluster->GetDetLocalZ());
          fAliITSQADataMakerRec->GetRecPointsData(14 +fRecsOffset)->Fill(cluGlo[0]);
          fAliITSQADataMakerRec->GetRecPointsData(15 +fRecsOffset)->Fill(cluGlo[1]);
          fAliITSQADataMakerRec->GetRecPointsData(16 +fRecsOffset)->Fill(cluGlo[2]);
          fAliITSQADataMakerRec->GetRecPointsData(17 +fRecsOffset)->Fill(rad);
          fAliITSQADataMakerRec->GetRecPointsData(18 +fRecsOffset)->Fill(phi);
          fAliITSQADataMakerRec->GetRecPointsData(19 +fRecsOffset)->Fill(cluster->GetNz(),cluster->GetNy());
          fAliITSQADataMakerRec->GetRecPointsData(20 +fRecsOffset)->Fill(cluGlo[2],phi);
        }
        fAliITSQADataMakerRec->GetRecPointsData(21 +fRecsOffset)->Fill(rad,phi);
        fAliITSQADataMakerRec->GetRecPointsData(22 +fRecsOffset)->Fill(cluGlo[0],cluGlo[1]);

	nClusters[lay]++; 
      } // end of cluster loop
    } // end of its "subdetector" loop

    for (Int_t iLay=0; iLay<2; iLay++)
      fAliITSQADataMakerRec->GetRecPointsData(23+iLay +fRecsOffset)->Fill(nClusters[iLay]);

    fAliITSQADataMakerRec->GetRecPointsData(25 +fRecsOffset)->Fill(nClusters[0],nClusters[1]);

    if (itsClusters) {
      itsClusters->Delete();
      delete itsClusters;
      itsClusters = 0;
    }

}
