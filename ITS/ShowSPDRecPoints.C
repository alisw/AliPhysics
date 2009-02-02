//////////////////////////////////////////////////////////////////
// Macro to check clusters in the 2 SPD layers                  //
// Provides:                                                    //
//  2 canvases with                                             //
//     - cluster loc and glob coordinates  (each layer)         //
//  1 canvas with                                               //
//      - cluster glob coordinates 2D and 3D                    //
//  1 canvas with                                               //
//      - correlations of #clusters for sectors                 //
//      - correlations of #clusters for half-sectors            //
//                                                              //
//  Maria.Nicassio@ba.infn.it                                   //
//  Domenico.Elia@ba.infn.it                                    //
//////////////////////////////////////////////////////////////////

#if !defined(__CINT__) || defined(__MAKECINT__)
                                                                                
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TClassTable.h>
#include <TInterpreter.h> 
#include <TGraph2D.h>
#include <TGraph.h>
#include <TGeoManager.h>

#include "AliCDBManager.h"
#include "AliRunLoader.h"
#include "AliESD.h"
#include "AliRun.h"
#include "AliGeomManager.h"                                                                                
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSLoader.h"
#include <AliITSRecPoint.h>
#include <AliHeader.h>                                                                                
#include <AliITSDetTypeRec.h>                                                                                

#endif

/* $Id$ */

void ShowSPDRecPoints(Int_t RunStart, Int_t RunStop){

  // Set data directory
  Char_t str[256];
  Char_t* dir = "/data/alipix/PhysicsEvents/pp/ppProd/pdc07/oldgeom";        
 
  // Variables for histo booking and filling 
  Int_t modmin=0;
  Int_t modmax=240;
  Int_t totmod=modmax-modmin;
  
  Float_t xlim[2]={4.5,8.};   
  Float_t zlim[2]={15.,15.};

  Int_t nClusters[2];

  Float_t clustersXCoord[2][200];
  Float_t clustersYCoord[2][200];
  Float_t clustersZCoord[2][200];
 
  Float_t cluGlo[3]={0.,0.,0.};  
 
  Int_t nClustersPerLayerPerSector[2][10];
  Int_t nClustersPerLayerPerHalfSector[2][20];

  Int_t iSector=0;  
  Int_t iHalfSector=0;

  // Booking of histograms  
  gStyle->SetPalette(1,0);                                                                                                     
  TH1F* hlayer= new TH1F("hlayer","",6,0.,6.);
  TH1F** hmod = new TH1F*[2];
  TH1F** hxl  = new TH1F*[2];
  TH1F** hzl  = new TH1F*[2];
  TH1F** hxg  = new TH1F*[2];
  TH1F** hyg  = new TH1F*[2];
  TH1F** hzg  = new TH1F*[2];
  TH1F** hr   = new TH1F*[2];
  TH1F** hphi = new TH1F*[2];
  TH1F** hMultSPDcl = new TH1F*[2];
  TH1F** hType = new TH1F*[2];  // cluster type ?

  TH2F** hr_phi   = new TH2F*[2];
  TH2F** hx_y   = new TH2F*[2];
  TH3F** hx_y_z   = new TH3F*[2];

  TH2F** hnCl2_nCl1_Sec = new TH2F*[10];
  TProfile** hnCl2vsnCl1_Sec = new TProfile*[10];
  TH2F** hnCl2_nCl1_HSec = new TH2F*[20];
//  TProfile** hnCl2vsnCl1_HSec = new TProfile*[20];
  TH2F** hNy_Nz = new TH2F*[2];  // y and z length
  TH2F** hPhi_Z = new TH2F*[2];
   
  TH2F *hMultSPDcl2_MultSPDcl1 = 
            new TH2F("nCl2_nCl1","# SPD clusters on layer 2 vs # on layer 1",200,0.,200.,200,0.,200.);
  hMultSPDcl2_MultSPDcl1->GetXaxis()->SetTitle("# clusters (inner layer)");
  hMultSPDcl2_MultSPDcl1->GetYaxis()->SetTitle("# clusters (outer layer)");                                                                                                                         
  Char_t name[10];
  Char_t title[20];
  for (Int_t iLay=0;iLay<2;iLay++) {
    sprintf(name,"hmod%d",iLay+1);
    hmod[iLay]=new TH1F(name,"SPD clusters - Module number",totmod,modmin,modmax);
    hmod[iLay]->GetXaxis()->SetTitle("Module number");
    hmod[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hxloc%d",iLay+1);
    hxl[iLay]=new TH1F(name,"SPD clusters - Local x coordinate",100,-4.,4.);
    hxl[iLay]->GetXaxis()->SetTitle("Local x [cm]");
    hxl[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hzloc%d",iLay+1);
    hzl[iLay]=new TH1F(name,"SPD clusters - Local z coordinate",100,-4.,4.);
    hzl[iLay]->GetXaxis()->SetTitle("Local z [cm]");
    hzl[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hxgl%d",iLay+1);
    hxg[iLay]=new TH1F(name,"SPD clusters - Global x coordinate",100,-xlim[iLay],xlim[iLay]);
    hxg[iLay]->GetXaxis()->SetTitle("Global x [cm]");
    hxg[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hygl%d",iLay+1);
    hyg[iLay]=new TH1F(name,"SPD clusters - Global y coordinate",100,-xlim[iLay],xlim[iLay]);
    hyg[iLay]->GetXaxis()->SetTitle("Global y [cm]");
    hyg[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hzgl%d",iLay+1);
    hzg[iLay]=new TH1F(name,"SPD clusters - Global z coordinate",150,-zlim[iLay],zlim[iLay]);
    hzg[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hzg[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hr%d",iLay+1);
    hr[iLay]=new TH1F(name,"SPD clusters - r",100,0.,50.);
    hr[iLay]->GetXaxis()->SetTitle("r [cm]");
    hr[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hphi%d",iLay+1);
    hphi[iLay]=new TH1F(name,"SPD clusters - #phi",100,0.,2*TMath::Pi()); 
    hphi[iLay]->GetXaxis()->SetTitle("#phi [rad]");
    hphi[iLay]->GetYaxis()->SetTitle("Entries");
    sprintf(name,"hType%d",iLay+1);
    hType[iLay]=new TH1F(name,"SPD clusters - Type",100,0.,100.);
    hType[iLay]->GetXaxis()->SetTitle("Cluster type");
    hType[iLay]->GetYaxis()->SetTitle("Entries");

    sprintf(name,"hMultSPDcl%d",iLay+1);
    hMultSPDcl[iLay]=new TH1F(name,"Cluster multiplicity",200,0.,200.);
    hMultSPDcl[iLay]->GetXaxis()->SetTitle("Cluster multiplicity");
    hMultSPDcl[iLay]->GetYaxis()->SetTitle("Entries");

    sprintf(name,"hNy_Nz%d",iLay+1);
    hNy_Nz[iLay]=new TH2F(name,"SPD clusters - Length",100,0.,100.,100,0.,100.);
    hNy_Nz[iLay]->GetXaxis()->SetTitle("z length");
    hNy_Nz[iLay]->GetYaxis()->SetTitle("y length");
    sprintf(name,"hPhi_Z%d",iLay+1);
    hPhi_Z[iLay]=new TH2F(name,"SPD clusters - #phi vs z",150,-zlim[iLay],zlim[iLay],100,0.,2*TMath::Pi());
    hPhi_Z[iLay]->GetXaxis()->SetTitle("Global z [cm]");
    hPhi_Z[iLay]->GetYaxis()->SetTitle("#phi [rad]");
    sprintf(name,"hr_phi%d",iLay+1);
    hr_phi[iLay]=new TH2F(name,"SPD clusters - #phi_r",500,0.,50.,100,0.,2*TMath::Pi());
    hr_phi[iLay]->GetXaxis()->SetTitle("r [cm]");
    hr_phi[iLay]->GetYaxis()->SetTitle("#phi [rad]");
    sprintf(name,"hx_y%d",iLay+1);
    hx_y[iLay]=new TH2F(name,"SPD clusters - y_x",200,-10.,10.,200,-10.,10.);
    hx_y[iLay]->GetXaxis()->SetTitle("x [cm]");
    hx_y[iLay]->GetYaxis()->SetTitle("y [cm]");
    sprintf(name,"hx_y_z%d",iLay+1);
    hx_y_z[iLay]=new TH3F(name,"SPD clusters - y_x",200,-10.,10.,200,-10.,10.,150,-15.,15.);
    hx_y_z[iLay]->GetXaxis()->SetTitle("z [cm]");
    hx_y_z[iLay]->GetYaxis()->SetTitle("x [cm]");
    hx_y_z[iLay]->GetZaxis()->SetTitle("y [cm]");
  }
  for (Int_t iSec=0; iSec<10; iSec++) {
    sprintf(name,"hnCl2_nCl1_Sector%d",iSec);
    sprintf(title,"Sector %d",iSec+1);
    hnCl2_nCl1_Sec[iSec]=new TH2F(name,title,200,0.,200.,200,0.,200.);
    hnCl2_nCl1_Sec[iSec]->GetXaxis()->SetTitle("# clusters (inner layer)");
    hnCl2_nCl1_Sec[iSec]->GetYaxis()->SetTitle("# clusters (outer layer)");
    sprintf(name,"hnCl2vsnCl1_Sector%d",iSec);
    sprintf(title,"Sector %d",iSec+1);
    hnCl2vsnCl1_Sec[iSec]=new TProfile(name,title,200,0.,200.,0.,200.);
    hnCl2vsnCl1_Sec[iSec]->GetXaxis()->SetTitle("# clusters (inner layer)");
    hnCl2vsnCl1_Sec[iSec]->GetYaxis()->SetTitle("# clusters (outer layer)");
  }
  for (Int_t iHalfSec=0; iHalfSec<20; iHalfSec++) {
    sprintf(name,"hnCl2_nCl1_HalfSector%d",iHalfSec);
    sprintf(title,"Half-Sector %d",iHalfSec+1);
    hnCl2_nCl1_HSec[iHalfSec]=new TH2F(name,title,200,0.,200.,200,0.,200.);
    hnCl2_nCl1_HSec[iHalfSec]->GetXaxis()->SetTitle("# clusters (inner layer)");
    hnCl2_nCl1_HSec[iHalfSec]->GetYaxis()->SetTitle("# clusters (outer layer)");
  }

  // Loop over runs
  for (Int_t run=RunStart; run<RunStop+1; run++) {
                                                                                
    // setup galice and runloader
//    cout << "File nr --> " << run << endl;
 
    if (gClassTable->GetID("AliRun") < 0) {
      gInterpreter->ExecuteMacro("loadlibs.C");
    }
    else { 
      if (gAlice){                        
        delete AliRunLoader::Instance();   
        delete gAlice;                    
        gAlice=0;                        
      }                                  
    }

    // Set OfflineConditionsDataBase if needed
    AliCDBManager* man = AliCDBManager::Instance();
    if (!man->IsDefaultStorageSet()) {
      printf("Setting a local default storage and run number 0\n");
      man->SetDefaultStorage("local://$ALICE_ROOT");
      man->SetRun(0);
    }
    else {
      printf("Using deafult storage \n");
    }
 
    // retrives geometry 
    if (!gGeoManager) {
      sprintf(str,"%s/ppMinBias%04d/geometry.root",dir,run);
      AliGeomManager::LoadGeometry(str);  
    }

    sprintf(str,"%s/ppMinBias%04d/galice.root",dir,run);
    AliRunLoader*  rl = AliRunLoader::Open(str);

    if (rl == 0x0){                                 
      cerr<<"Can not open session RL=NULL"<< endl;  
      return;                                    
    }                                               
    Int_t retval = rl->LoadgAlice();                
    if (retval){
      cerr<<"LoadgAlice returned error"<<endl;
      return;
    }
    gAlice=rl->GetAliRun();                          

    retval = rl->LoadHeader();                       
    if (retval){
      cerr<<"LoadHeader returned error"<<endl;
      return;
    }

    AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");  
    if (!ITSloader){							  
      cerr<<"ITS loader not found"<<endl;					  
      return;								  
    }
    ITSloader->LoadRecPoints("read");                                       

    AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
    ITS->SetTreeAddress();
    AliITSgeom *geom = ITS->GetITSgeom();                                    
    if (!geom) {                                                             
      cout << " Can't get the ITS geometry!" << endl;                      
      return ;                                                             
    }                                                                      

    AliITSDetTypeRec* detTypeRec = new AliITSDetTypeRec(); 

    detTypeRec->SetITSgeom(ITSloader->GetITSgeom());
    detTypeRec->SetDefaults();
 
    Int_t nEvents=rl->GetNumberOfEvents();                      
    printf("Total Number of events = %d\n",nEvents);            

    for (Int_t iev=0;iev<nEvents;iev++){                         
      rl->GetEvent(iev);                                        

      TTree *TR = ITSloader->TreeR();                           
      TClonesArray* ITSClusters;
      ITSClusters  = detTypeRec->RecPoints();  
      TBranch *branch = 0;
      if (TR && ITSClusters) {
        branch = ITSloader->TreeR()->GetBranch("ITSRecPoints");
        if (branch) branch->SetAddress(&ITSClusters);
      }

//      Int_t nparticles = rl->GetHeader()->GetNtrack();
//      cout<<"Event # "<<iev<<"   #Particles="<<nparticles<<endl;

      // Reset cluster counters
      for (Int_t iLay=0; iLay<2; iLay++) {
        nClusters[iLay]=0;
        for (Int_t iSec=0; iSec<10; iSec++) {
          nClustersPerLayerPerSector[iLay][iSec]=0;
        }
        for (Int_t iHSec=0; iHSec<20; iHSec++) {
          nClustersPerLayerPerHalfSector[iLay][iHSec]=0;
        }
      }    
      for (Int_t subd=0;subd<1;subd++) {

        Int_t first = geom->GetStartDet(subd);
        Int_t last = geom->GetLastDet(subd);
        for (Int_t mod=first; mod<=last; mod++) {
          detTypeRec->ResetRecPoints();  
	  branch->GetEvent(mod);
	  Int_t nrecp = ITSClusters->GetEntries();     
          while (nrecp--) {
	  
            AliITSRecPoint *recp = (AliITSRecPoint*)ITSClusters->At(nrecp); 
            
	    Int_t lay=recp->GetLayer();
//            cout<<"lay"<<lay<<endl;
            recp->GetGlobalXYZ(cluGlo);

	    Float_t rad=TMath::Sqrt(cluGlo[0]*cluGlo[0]+cluGlo[1]*cluGlo[1]); 
	    Float_t phi= TMath::Pi() + TMath::ATan2(-cluGlo[1],-cluGlo[0]); 
            hlayer->Fill(lay);
	    hmod[lay]->Fill(mod);
	    hzl[lay]->Fill(recp->GetDetLocalZ());
	    hxl[lay]->Fill(recp->GetDetLocalX());
	    hzg[lay]->Fill(cluGlo[2]);
	    hyg[lay]->Fill(cluGlo[1]);
	    hxg[lay]->Fill(cluGlo[0]);
	    hr[lay]->Fill(rad);
	    hphi[lay]->Fill(phi);
            hPhi_Z[lay]->Fill(cluGlo[2],phi);       
            hr_phi[lay]->Fill(rad,phi);
            hx_y[lay]->Fill(cluGlo[0],cluGlo[1]);
            hx_y_z[lay]->Fill(cluGlo[2],cluGlo[0],cluGlo[1]);
            clustersXCoord[lay][nClusters[lay]]=cluGlo[0]; 
            clustersYCoord[lay][nClusters[lay]]=cluGlo[1];
            clustersZCoord[lay][nClusters[lay]]=cluGlo[2];
            nClusters[lay]++;
  
            hType[lay]->Fill(recp->GetType());
//            cout<<"Clusters type"<<recp->GetType()<<endl;
            hNy_Nz[lay]->Fill(recp->GetNz(),recp->GetNy());
            
            // Set Sector number and increase the counter
            if (lay==0) {
              for (Int_t nRange=0; nRange<10; nRange++) {
                if ((mod>=nRange*8) && (mod<=(nRange*8+7))) {
                  iSector = nRange;
                  nClustersPerLayerPerSector[lay][iSector]++;
                  break;
                }
              }
            }
            if (lay==1) {
              for (Int_t nRange=0; nRange<10; nRange++) {
                if ((mod>=80+nRange*16) && (mod<=(80+nRange*16+15))) {
                  iSector = nRange;
                  nClustersPerLayerPerSector[lay][iSector]++;
                  break;
                }
              }
            }
            // Set HalfSector number and increase the counter
            if (lay==0) {
              for (Int_t nRange=0; nRange<20; nRange++) {
                 if ((mod>=nRange*4) && (mod<=(nRange*4+3))) {
                   iHalfSector = nRange;
                   nClustersPerLayerPerHalfSector[lay][iHalfSector]++;
                   break;
                 }
              }
            } else {
              for (Int_t nRange=0; nRange<20; nRange++) {
                if ((mod>=80+nRange*4) && (mod<=(80+nRange*4+3))) {
                  iHalfSector = nRange;
                  nClustersPerLayerPerHalfSector[lay][iHalfSector]++;
                  break;
                }
              }
            }
          }
	}
      }

      for (Int_t iLay=0; iLay<2; iLay++)   hMultSPDcl[iLay]->Fill(nClusters[iLay]);
      for (Int_t iSec=0; iSec<10; iSec++) {
        hnCl2_nCl1_Sec[iSec]->Fill(nClustersPerLayerPerSector[0][iSec],nClustersPerLayerPerSector[1][iSec]);
        hnCl2vsnCl1_Sec[iSec]->Fill(nClustersPerLayerPerSector[0][iSec],nClustersPerLayerPerSector[1][iSec]);
      }
      for (Int_t iHSec=0; iHSec<20; iHSec++) {
        hnCl2_nCl1_HSec[iHSec]->Fill(nClustersPerLayerPerHalfSector[0][iHSec],nClustersPerLayerPerHalfSector[1][iHSec]);
      }   
      hMultSPDcl2_MultSPDcl1->Fill(nClusters[0],nClusters[1]);
    
    } //end loop over events
    rl->UnloadAll();
    delete rl;
                                                                                
  } //end loop over runs

  // Draw and Write histos

  TFile* fout = new TFile("out_ShowSPDRecPoints.root","RECREATE");
  
  hlayer->Write();
//  cev0->Write();
  
  TCanvas **c=new TCanvas*[2];
  Char_t ctit[30];
  for(Int_t iLay=0;iLay<2;iLay++){
    hNy_Nz[iLay]->Write();
    sprintf(name,"can%d",iLay+1);
    sprintf(ctit,"Layer %d",iLay+1);
    c[iLay]=new TCanvas(name,ctit,1200,900);
    c[iLay]->Divide(3,3);
    c[iLay]->cd(1);
    hmod[iLay]->Draw();
    hmod[iLay]->Write();
    c[iLay]->cd(2);
    hxl[iLay]->Draw();
    hxl[iLay]->Write();
    c[iLay]->cd(3);
    hzl[iLay]->Draw();
    hzl[iLay]->Write();
    c[iLay]->cd(4);
    hxg[iLay]->Draw();
    hxg[iLay]->Write();
    c[iLay]->cd(5);
    hyg[iLay]->Draw();
    hyg[iLay]->Write();
    c[iLay]->cd(6);
    hzg[iLay]->Draw();
    hzg[iLay]->Write();
    c[iLay]->cd(7);
    hr[iLay]->Draw();
    hr[iLay]->Write();
    c[iLay]->cd(8);
    hphi[iLay]->Draw();   
    hphi[iLay]->Write();
    c[iLay]->cd(9);
    hPhi_Z[iLay]->Draw("colz"); 
    hPhi_Z[iLay]->Write();
    hr_phi[iLay]->Write();
    hx_y[iLay]->Write();
    hx_y_z[iLay]->Write();
  }

  TCanvas *cCoord=new TCanvas("cCoord","Cluster coordinates",1200,900);
  cCoord->Divide(2,2);

  for (Int_t iLay=0;iLay<2;iLay++) {
    cCoord->cd(1);
    hx_y[iLay]->SetMarkerStyle(22);
    hx_y[iLay]->SetMarkerSize(0.3);
    hx_y[iLay]->SetMarkerColor(iLay+1);
    if (iLay==0) hx_y[iLay]->Draw("p");
    else hx_y[iLay]->Draw("p,same");
    cCoord->cd(2);
    hx_y_z[iLay]->SetMarkerStyle(23);
    hx_y_z[iLay]->SetMarkerSize(0.3);
    hx_y_z[iLay]->SetMarkerColor(iLay+1);
    hx_y_z[iLay]->Draw("p,same");
    cCoord->cd(3);
    hr_phi[iLay]->SetMarkerColor(iLay+1);
    hr_phi[iLay]->Draw("p,same");
  }
  TCanvas *cCorrelations_Sectors=new TCanvas("cCorrelations_S","SPD cluster correlations (Sectors)",1200,900);
  cCorrelations_Sectors->Divide(3,5);
  cCorrelations_Sectors->cd(1);
  hMultSPDcl2_MultSPDcl1->Draw();
  hMultSPDcl2_MultSPDcl1->Write();
  cCorrelations_Sectors->cd(2);
  hMultSPDcl[0]->Draw();
  hMultSPDcl[0]->Write();
  cCorrelations_Sectors->cd(3);
  hMultSPDcl[1]->Draw();
  hMultSPDcl[1]->Write();
  for (Int_t iS=0; iS<10; ++iS) { 
    cCorrelations_Sectors->cd(iS+4); 
    hnCl2_nCl1_Sec[iS]->Draw();
    hnCl2_nCl1_Sec[iS]->Write();
    hnCl2vsnCl1_Sec[iS]->Write();
  }
 
  TCanvas *cCorrelations_HSectors=new TCanvas("cCorrelations_HS","SPD cluster correlations (Half-Sectors)",1200,900);
  cCorrelations_HSectors->Divide(5,4);

  for (Int_t iHS=0; iHS<20; ++iHS) {
    cCorrelations_HSectors->cd(iHS+1);
    hnCl2_nCl1_HSec[iHS]->Write();
    hnCl2_nCl1_HSec[iHS]->Draw(); 
  }

  fout->Close();

}
