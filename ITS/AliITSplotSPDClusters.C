#if !defined(__CINT__) || defined(__MAKECINT__)

// Standard includes
#include <iostream.h>
// Root includes
#include <TFile.h>
#include <TObjArray.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
// AliRoot includes
#include "STEER/AliRun.h"
#include "ITS/AliITS.h"
#include "ITS/AliITSRawCluster.h"

#endif

void AliITSplotSPDClusters(const char* filename="galice_80.root"){
    //------------------------------------------------------------------
    // This macro will read the TreeC produced during reconstruction and
    // plot the number and type of clusters found in the SPD.
    // root [0] .L AliITSplotSPDClusters.C  //this loads the macro in memory
    // root [1] AliITSplotSPDClusters();    //by default process first event
    // or
    // root [1] AliITSplotSPDClusters("galice2.root"); // process all events
    //                                                    from the file
    //                                                    galice2.root.
    // or
    // root [0] .x AliITSplotSPDClusters.C("galice2.root") // will all of the
    //                                                        events from the
    //                                                        galice2.root
    //------------------------------------------------------------------

    // delete the current gAlice object, the one from input file will be used
    if(gAlice){
	delete gAlice;
	gAlice = 0;
    }else{ // Dynamically link some shared libs
        if(gClassTable->GetID("AliRun") < 0) {
            gROOT->LoadMacro("loadlibs.C");
            loadlibs();
        } // end if
    } // end if gAlice

    // Connect the Root Galice file containing Geometry, Kine and Hits
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if(!file) file = new TFile(filename);
    // Get AliRun object from file or create it if not on file
    if(!gAlice) {
        gAlice = (AliRun*)file->Get("gAlice");
        if(gAlice) cout << "AliRun object found on file" << endl;
        if(!gAlice) gAlice = new AliRun("gAlice","SPD Clulster analysis");
    } // end if !gAlice
    Int_t nevents = gAlice->GetEventsPerRun();

    // Get pointers to the ITS Alice detectors.
    AliITS *ITS = (AliITS*) gAlice->GetDetector("ITS");
    if(ITS==0){
        cout << "ITS not found. Exiting." << endl;
        return;
    } // end if ITS==0
    AliITSresponseSPDdubna *resp0 = ITS->DetType(0)->GetResponseModel();
    Float_t diffCoeff= resp0->DistanceOverVoltage();//Get Value of Diffusion Coefficient parameter d/v.
    diffCoeff *= resp0->Temperature();
    diffCoeff = TMath::Sqrt(diffCoeff*8.6173376e-5);// K_Boltzman/e_Coulumb
    char    aDiffCoeff[100];            //Character array for sprintf

    sprintf(aDiffCoeff,"Number of SPD Clusters in x, layer 1, DiffCoeff=%f #sqrt{cm}",
	    diffCoeff);
    TH1D *hclx1 = new TH1D("hclx1",aDiffCoeff,15,0.5,15.5);
    sprintf(aDiffCoeff,"Number of SPD Clusters in z, layer 1, DiffCoeff=%f #sqrt{cm}",
	    diffCoeff);
    TH1D *hclz1 = new TH1D("hclz1",aDiffCoeff,5,0.5,5.5);
    sprintf(aDiffCoeff,"Number of SPD Clusters in x, layer 2, DiffCoeff=%f #sqrt{cm}",
	    diffCoeff);
    TH1D *hclx2 = new TH1D("hclx2",aDiffCoeff,15,0.5,15.5);
    sprintf(aDiffCoeff,"Number of SPD Clusters in z, layer 2, DiffCoeff=%f #sqrt{cm}",
	    diffCoeff);
    TH1D *hclz2 = new TH1D("hclz2",aDiffCoeff,5,0.5,5.5);
    // Create Arrays with clusters from:  Data, Ba/Sa Model, old version of
    // Dubna
    Float_t dataX1[9] = {0.493, 0.493, 0.0273, 0.00617, 0.00112, 0.000257,
			 0.000257,  0.000257, 0.000257};
                        //clusters from the data in the x-dimension with
                        //standard paramaters
    Float_t baSaX1[11] = {0.942, 0.0180, 0.0112, 0.00389, 0.00203, 0.000257,
			  0.001, 0.000257, 0.000257, 0.000257, 0.001};
                        //same for Ba/Sa model
    Float_t dubnaX1[7] =  {0.889, 0.0789, 0.0151, 0.00389, 0.001, 0.001,
			   0.001};
                        //same for old version of Dubna model
    Float_t dataX2[9] = {0.482, 0.482, 0.0325, 0.00791, 0.00157, 0.000275,
                         0.000275, 0.000275, 0.000275};
                        //clusters from the data in the x-dimension with
                        //optimized paramaters
    Float_t baSaX2[11] = {0.946, 0.0175, 0.0245, 0.00482, 0.001, 0.001,
			  0.000275,  0.001, 0.000275, 0.000275, 0.001};
                        //same for Ba/Sa model
    Float_t dubnaX2[8] = {0.638, 0.325, 0.0275, 0.01, 0.00431, 0.000275,
			  0.001, 0.001};
                        //same for old version of Dubna model
    Float_t dataZ1[4] = {0.889, 0.0624, 0.000257, 0.000257};
                        //clusters from the data in the z-dimension with
                        //standard paramaters
    Float_t baSaZ1[3] = {1., 0.0160, 0.000257}; //same for Ba/Sa model
    Float_t dubnaZ1[3] = {0.889, 0.0180, 0.000257};
                        //same for old version of Dubna model
    Float_t dataZ2[4] = {0.889, 0.0624, 0.000257, 0.000257};
                        //clusters from the data in the z-dimension with
                        //optimized paramaters
    Float_t baSaZ2[4] = {1., 0.0160, 0.00215, 0.000257}; //same for Ba/Sa model
    Float_t dubnaZ2[3] = {0.889, 0.0412, 0.000257};
                        //same for old version of Dubna model
    
    TH1F *hdataX1  = new TH1F("hdataX1","Number of SPD Clusters in x, layer 1",
			      9, 0.5, 9.5);
    hdataX1->SetMarkerColor(2);
    hdataX1->SetMarkerStyle(20);
    hdataX1->SetMarkerSize(0.7);
    TH1F *hbaSaX1  = new TH1F("hbaSaX1", "Ba/Sa", 11, 0.5, 11.5);
    TH1F *hdubnaX1 = new TH1F("hdubnaX1", "Old Dubna", 7, 0.5, 7.5);
    TH1F *hdataX2  = new TH1F("hdataX2","Number of SPD Clusters in x, layer 2",
			      9, 0.5, 9.5);
    hdataX2->SetMarkerColor(2);
    hdataX2->SetMarkerStyle(20);
    hdataX2->SetMarkerSize(0.7);
    TH1F *hbaSaX2  = new TH1F("hbaSaX2", "Ba/Sa", 11, 0.5, 11.5);
    TH1F *hdubnaX2 = new TH1F("hdubnaX2", "Old Dubna", 8, 0.5, 8.5);
    TH1F *hdataZ1  = new TH1F("hdataZ1","Number of SPD Clusters in z, layer 1",
			      4, 0.5, 4.5);
    hdataZ1->SetMarkerColor(2);
    hdataZ1->SetMarkerStyle(20);
    hdataZ1->SetMarkerSize(0.7);
    TH1F *hbaSaZ1  = new TH1F("hbaSaZ1", "Ba/Sa", 3, 0.5, 3.5);
    TH1F *hdubnaZ1 = new TH1F("hdubnaZ1", "Old Dubna", 3, 0.5, 3.5);
    TH1F *hdataZ2  = new TH1F("hdataZ2","Number of SPD Clusters in z, layer 2",
			      4, 0.5, 4.5);
    hdataZ2->SetMarkerColor(2);
    hdataZ2->SetMarkerStyle(20);
    hdataZ2->SetMarkerSize(0.7);
    TH1F *hbaSaZ2  = new TH1F("hbaSaZ2", "Ba/Sa", 4, 0.5, 4.5);
    TH1F *hdubnaZ2 = new TH1F("hdubnaZ2", "Old Dubna", 3, 0.5, 3.5);

    Int_t j = 0; //j will be a counter for the loops to fill the histograms
                 //with the values for clusters for the Data, the Ba/Sa model
                 //and the old Dubna model
    for(j=0; j<9; j++){
	hdataX1->SetBinContent((j+1), (Double_t)dataX1[j]);
	hdataX1->SetBinError((j+1), 0.00001);
    }
    for(j=0; j<11; j++) hbaSaX1->Fill((Float_t)j +0.5, baSaX1[j]);
    for(j=0; j<7; j++)  hdubnaX1->Fill((Float_t)j +0.5, dubnaX1[j]);
    for(j=0; j<9; j++){
	hdataX2->SetBinContent((j+1), (Double_t)dataX2[j]);
	hdataX2->SetBinError((j+1), 0.00001);
    }
    for(j=0; j<11; j++) hbaSaX2->Fill((Float_t)j +0.5, baSaX2[j]);
    for(j=0; j<8; j++)  hdubnaX2->Fill((Float_t)j +0.5, dubnaX2[j]);
    for(j=0; j<4; j++){
	hdataZ1->SetBinContent((j+1), (Double_t)dataZ1[j]);
	hdataZ1->SetBinError((j+1), 0.00001);
    }
    for(j=0; j<3; j++)  hbaSaZ1->Fill((Float_t)j +0.5, baSaZ1[j]);
    for(j=0; j<3; j++)  hdubnaZ1->Fill((Float_t)j +0.5, dubnaZ1[j]);
    for(j=0; j<4; j++){
	hdataZ2->SetBinContent((j+1), (Double_t)dataZ2[j]);
	hdataZ2->SetBinError((j+1), 0.00001);
    }
    for(j=0; j<4; j++)  hbaSaZ2->Fill((Float_t)j +0.5, baSaZ2[j]);
    for(j=0; j<3; j++)  hdubnaZ2->Fill((Float_t)j +0.5, dubnaZ2[j]);

    TTree *tc = 0;
    TBranch *spdBranch = 0;
    TClonesArray *spdClustcA = new TClonesArray("AliITSRawClusterSPD",1000);
    AliITSRawClusterSPD *spdClust = 0;
    char tn[20];
    Int_t evnt,i,n,nc;
    Float_t nclx = 0,nclz = 0;
    for(evnt=0;evnt<nevents;evnt++){
	sprintf(tn,"TreeC%d",evnt);
	tc = (TTree*) file->Get(tn);
	spdBranch = tc->GetBranch("ITSClustersSPD");
	spdBranch->SetAddress(&spdClustcA);
  	n = (Int_t) tc->GetEntries();
	for(i=0;i<n;i++){
  	    tc->GetEvent(i);
	    nc = spdClustcA->GetEntries();
	    if(nc>240) nc = 240;
	    for(j=0;j<nc;j++){
		spdClust = (AliITSRawClusterSPD*) spdClustcA->At(j);
		nclx = spdClust->NclX();
		nclz = spdClust->NclZ();
		if(spdClust->Module()<80){
		    hclx1->Fill(nclx-0.5);
		    hclz1->Fill(nclz-0.5);
		}else if(spdClust->Module()<240){
		    hclx2->Fill(nclx-0.5);
		    hclz2->Fill(nclz-0.5);
		} //endif
	    } // end for j
	} // end for i
	spdClustcA->Clear();
	delete spdBranch; spdBranch = 0;
	delete tc; tc = 0;
    } // end for evnt

    //Begin Process of normalizing histograms
    Double_t integral = (Double_t) hclx1->Integral();
    if(integral>0.0) hclx1->Scale(1./integral);
    integral = hclz1->Integral();
    if(integral>0.0) hclz1->Scale(1./integral);
    integral = hclx2->Integral();
    if(integral>0.0) hclx2->Scale(1./integral);
    integral = hclz2->Integral();
    if(integral>0.0) hclz2->Scale(1./integral);

    hclx1->SetMinimum(0.000257);
    hclx1->SetMaximum(1.01);
    hclz1->SetMinimum(0.000274);
    hclz1->SetMaximum(1.01);
    hclx2->SetMinimum(0.000257);
    hclx2->SetMaximum(1.01);
    hclz2->SetMinimum(0.000257);
    hclz2->SetMaximum(1.01);

    sprintf(aDiffCoeff,"SPD Clusters with Diffusion Coefficent=%f",diffCoeff);
    TCanvas *cSPDclusters = new TCanvas("cSPDclusters",aDiffCoeff,
				       400,10,600,776);
    cSPDclusters->Divide(2, 2);

    cSPDclusters->cd(1);
    cSPDclusters->SetLogy(1);
    hclx1->Draw("hist");
    hdataX1->Draw("same e1p");
    hbaSaX1->SetLineColor(4);
    hbaSaX1->SetLineStyle(2);    
    hbaSaX1->Draw("same");
    hdubnaX1->SetLineColor(3);
    hdubnaX1->SetLineStyle(4);
    hdubnaX1->Draw("same");
    TLegend *l1 = new TLegend(0.55,0.65,0.76,0.82);
    l1->AddEntry(hclx1,"New simulation","l");
    l1->AddEntry(hdataX1,"Test Beam Results","p");
    l1->AddEntry(hbaSaX1,"Bari/Selero simulation","l");
    l1->AddEntry(hdubnaX1,"Dubna simulation","l");
    l1->Draw();

    cSPDclusters->cd(2);
    cSPDclusters->SetLogy(1);
    hclz1->Draw("hist");
    hdataZ1->Draw("same e1p");
    hbaSaZ2->SetLineColor(4);
    hbaSaZ1->SetLineStyle(2);
    hbaSaZ1->Draw("same");
    hdubnaZ1->SetLineColor(3);
    hdubnaZ1->SetLineStyle(4);
    hdubnaZ1->Draw("same");
    TLegend *l2 = new TLegend(0.55,0.65,0.76,0.82);
    l2->AddEntry(hclz1,"New simulation","l");
    l2->AddEntry(hdataZ1,"Test Beam Results","p");
    l2->AddEntry(hbaSaZ1,"Bari/Selero simulation","l");
    l2->AddEntry(hdubnaZ1,"Dubna simulation","l");
    l2->Draw();

    cSPDclusters->cd(3);
    cSPDclusters->SetLogy(1);
    hclx2->Draw("hist");
    TLegend *l3 = new TLegend(0.55,0.65,0.76,0.82);
    l3->AddEntry(hclx2,"New simulation","l");
    l3->Draw();
/*
    hdataX2->Draw("e1p");
    hbaSaX2->SetLineColor(4);
    hbaSaX2->SetLineStyle(2);
    hbaSaX2->Draw("same");
    hdubnaX2->SetLineColor(3);
    hdubnaX2->SetLineStyle(4);
    hdubnaX2->Draw("same");
*/
    cSPDclusters->cd(4);
    cSPDclusters->SetLogy(1);
    hclz2->Draw("hist");
    TLegend *l4 = new TLegend(0.55,0.65,0.76,0.82);
    l4->AddEntry(hclz2,"New simulation","l");
    l4->Draw();
/*
    hdataZ2->Draw("e1p");
    hbaSaZ2->SetLineColor(4);
    hbaSaZ2->SetLineStyle(2);
    hbaSaZ2->Draw("same");
    hdubnaZ2->SetLineColor(3);
    hdubnaZ2->SetLineStyle(4);
    hdubnaZ2->Draw("same");
*/
    cSPDclusters->Update();

    if(gROOT->IsBatch()){
	cSPDclusters->Print("SPDClusters.eps");
    } // end if gROOT->IsBatch()

}
