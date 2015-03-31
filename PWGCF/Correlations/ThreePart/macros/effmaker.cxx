#include <iostream>
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TNamed.h"
#include "TList.h"
#include "TDirectory.h"
#include <X3DDefs.h>
#include "TPaveText.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TObjString.h"
#include "THn.h"


using namespace std;

THnF * Rebin(Int_t Nmult,const Double_t * MultAxis,Int_t nvz, const Double_t * VZAxis,Int_t nphi, Double_t phimin, Double_t phimax,Int_t neta, Double_t etamin,Double_t etamax,Int_t npt, const Double_t * pTAxis,THnF * input, const char* name){
  //Function that takes a THnF, and rebins it to the given axes.
  Int_t    nbins[5] 	= {Nmult,nvz,nphi,neta,npt};
  Double_t xmin[5] 	= {MultAxis[0],VZAxis[0],phimin,etamin,pTAxis[0]};
  Double_t xmax[5] 	= {MultAxis[Nmult],VZAxis[nvz],phimax,etamax,pTAxis[npt]};
  
  THnF * Resulthist = new THnF(name,input->GetTitle(), 5,nbins,xmin,xmax);
  Resulthist->GetAxis(0)->SetTitle(input->GetAxis(0)->GetTitle());
  Resulthist->GetAxis(0)->Set(Nmult,MultAxis);
  Resulthist->GetAxis(1)->SetTitle(input->GetAxis(1)->GetTitle());
  Resulthist->GetAxis(1)->Set(nvz,VZAxis);
  Resulthist->GetAxis(2)->SetTitle(input->GetAxis(2)->GetTitle());
  Resulthist->GetAxis(3)->SetTitle(input->GetAxis(3)->GetTitle());
  Resulthist->GetAxis(4)->SetTitle(input->GetAxis(4)->GetTitle());
  Resulthist->GetAxis(4)->Set(npt,pTAxis);

  for(int mult = 0;mult<=input->GetAxis(0)->GetNbins();mult++){
    //find the corresponding bin:
    double multval = input->GetAxis(0)->GetBinCenter(mult);
    for(int vz = 0;vz<=input->GetAxis(1)->GetNbins();vz++){
      //find the corresponding bin:
      double vzval = input->GetAxis(1)->GetBinCenter(vz);
      for(int phi = 0;phi<=input->GetAxis(2)->GetNbins();phi++){
	//find the corresponding bin:
	double phival = input->GetAxis(2)->GetBinCenter(phi);
	for(int eta = 0;eta<=input->GetAxis(3)->GetNbins();eta++){
	  //find the corresponding bin:
	  double etaval = input->GetAxis(3)->GetBinCenter(eta);
	  for(int pT = 0;pT<=input->GetAxis(4)->GetNbins();pT++){
	    //find the corresponding bin:
	    double pTval = input->GetAxis(4)->GetBinCenter(pT);
	    int bin[5] = {mult,vz,phi,eta,pT};
	    Double_t value[5] = {multval,vzval,phival,etaval,pTval};
	    Resulthist->Fill(value,input->GetBinContent(bin));
	  }	
	}
      }
    }
  }

  return Resulthist;
}

void MakeEffHistspp(){
    //open the input file:
  TFile * infile = TFile::Open("AnalysisResults.root","READ");
  THnF * hnTracksInBins;
  THnF * hnTracksInBinsRecPP;
  THnF * hnTracksInBinsMC;
  TList * dir = dynamic_cast<TList*>(infile->GetDirectory("AddTaskThreePartTrackEfficiencies_0_16")->Get("AddTaskThreePartTrackEfficiencies_0_16Coutput1"));
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBins = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsRecPP")==0)hnTracksInBinsRecPP = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMC = dynamic_cast<THnF*>(dir->At(i));
  }
  //Axis: 0 - Centrality or multiplicity, 1 - Vertex , 2 - phi , 3 - eta , 4 - pT
  TFile* outfile = TFile::Open("eff.root","RECREATE");
  TAxis * multaxis = hnTracksInBins->GetAxis(0);
  multaxis->SetTitle("Centrality [%]");
  hnTracksInBinsRecPP->GetAxis(0)->SetTitle("Centrality [%]");
  hnTracksInBinsMC->GetAxis(0)->SetTitle("Centrality [%]");
  TAxis * vzAxis = hnTracksInBins->GetAxis(1);
  vzAxis->SetTitle("Z-Vertex [cm]");  
  hnTracksInBinsRecPP->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  hnTracksInBinsMC->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  TAxis * phiaxis = hnTracksInBins->GetAxis(2);
  phiaxis->SetTitle("#Phi [radians]");    
  hnTracksInBinsRecPP->GetAxis(2)->SetTitle("#Phi [radians]");
  hnTracksInBinsMC->GetAxis(2)->SetTitle("#Phi [radians]");
  TAxis * Etaaxis = hnTracksInBins->GetAxis(3);
  Etaaxis->SetTitle("#eta []");  
  hnTracksInBinsRecPP->GetAxis(3)->SetTitle("#eta []");
  hnTracksInBinsMC->GetAxis(3)->SetTitle("#eta []");
  TAxis * pTaxis = hnTracksInBins->GetAxis(4);
  pTaxis->SetTitle("pT [GeV/c]");    
  hnTracksInBinsRecPP->GetAxis(4)->SetTitle("pT [GeV/c]");
  hnTracksInBinsMC->GetAxis(4)->SetTitle("pT [GeV/c]");
  
  
  TH1D * multRec = hnTracksInBins->Projection(0);
  TH1D * multRecPP = hnTracksInBinsRecPP->Projection(0);
  TH1D * multMC = hnTracksInBinsMC->Projection(0);
  multRec->Sumw2();
  multRec->Write("Multiplicity_All_Reconstructed");
  multRecPP->Sumw2();
  multRecPP->Write("Multiplicity_Reconstructed_from_Physical_Primary");  
  multMC->Sumw2();
  multMC->Write("Multiplicity_All_Produced");
  
  TH1D * multEffRec = dynamic_cast<TH1D*>(multRec->Clone("Efficiency_Mult_Rec"));
  multEffRec->Divide(multMC);
  multEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRec->Write();
  TH1D * multEffRecPP = dynamic_cast<TH1D*>(multRecPP->Clone("Efficiency_Mult_Rec_pp"));
  multEffRecPP->Divide(multMC);
  multEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  multEffRecPP->Write();

  TH1D * vzRec = hnTracksInBins->Projection(1);
  TH1D * vzRecPP = hnTracksInBinsRecPP->Projection(1);
  TH1D * vzMC = hnTracksInBinsMC->Projection(1);
  vzRec->Sumw2();
  vzRec->Write("VZ_All_Reconstructed");
  vzRecPP->Sumw2();
  vzRecPP->Write("VZ_Reconstructed_from_Physical_Primary");  
  vzMC->Sumw2();
  vzMC->Write("VZ_All_Produced");
  
  TH1D * vzEffRec = dynamic_cast<TH1D*>(vzRec->Clone("Efficiency_VZ_Rec"));
  vzEffRec->Divide(vzMC);
  vzEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRec->Write();
  TH1D * vzEffRecPP = dynamic_cast<TH1D*>(vzRecPP->Clone("Efficiency_VZ_Rec_pp"));
  vzEffRecPP->Divide(vzMC);
  vzEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  vzEffRecPP->Write();
  
  TH1D * phiRec = hnTracksInBins->Projection(2);
  TH1D * phiRecPP = hnTracksInBinsRecPP->Projection(2);
  TH1D * phiMC = hnTracksInBinsMC->Projection(2);
  phiRec->Sumw2();
  phiRec->Write("phi_All_Reconstructed");
  phiRecPP->Sumw2();
  phiRecPP->Write("phi_Reconstructed_from_Physical_Primary");  
  phiMC->Sumw2();
  phiMC->Write("phi_All_Produced");
  
  TH1D * phiEffRec = dynamic_cast<TH1D*>(phiRec->Clone("Efficiency_phi_Rec"));
  phiEffRec->Divide(phiMC);
  phiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiEffRec->Write();
  TH1D * phiEffRecPP = dynamic_cast<TH1D*>(phiRecPP->Clone("Efficiency_phi_Rec_pp"));
  phiEffRecPP->Divide(phiMC);
  phiEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  phiEffRecPP->Write();
  
  TH1D * etaRec = hnTracksInBins->Projection(3);
  TH1D * etaRecPP = hnTracksInBinsRecPP->Projection(3);
  TH1D * etaMC = hnTracksInBinsMC->Projection(3);
  etaRec->Sumw2();
  etaRec->Write("eta_All_Reconstructed");
  etaRecPP->Sumw2();
  etaRecPP->Write("eta_Reconstructed_from_Physical_Primary");  
  etaMC->Sumw2();
  etaMC->Write("eta_All_Produced");
  
  TH1D * etaEffRec = dynamic_cast<TH1D*>(etaRec->Clone("Efficiency_eta_Rec"));
  etaEffRec->Divide(etaMC);
  etaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaEffRec->Write();
  TH1D * etaEffRecPP = dynamic_cast<TH1D*>(etaRecPP->Clone("Efficiency_eta_Rec_pp"));
  etaEffRecPP->Divide(etaMC);
  etaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  etaEffRecPP->Write();

  TH1D * pTRec = hnTracksInBins->Projection(4);
  TH1D * pTRecPP = hnTracksInBinsRecPP->Projection(4);
  TH1D * pTMC = hnTracksInBinsMC->Projection(4);
  pTRec->Sumw2();
  pTRec->Write("pT_All_Reconstructed");
  pTRecPP->Sumw2();
  pTRecPP->Write("pT_Reconstructed_from_Physical_Primary");  
  pTMC->Sumw2();
  pTMC->Write("pT_All_Produced");
  
  TH1D * pTEffRec = dynamic_cast<TH1D*>(pTRec->Clone("Efficiency_pT_Rec"));
  pTEffRec->Divide(pTMC);
  pTEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  pTEffRec->Write();
  TH1D * pTEffRecPP = dynamic_cast<TH1D*>(pTRecPP->Clone("Efficiency_pT_Rec_pp"));
  pTEffRecPP->Divide(pTMC);
  pTEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  pTEffRecPP->Write();
  
  //2d projections in correlated variables:
  TH2D * vzetaRec = hnTracksInBins->Projection(3,1);
  vzetaRec->Sumw2();
  vzetaRec->Write("VZ_Eta_All_Reconstructed");
  TH2D * vzetaRecPP = hnTracksInBinsRecPP->Projection(3,1);
  vzetaRecPP->Sumw2();
  vzetaRecPP->Write("VZ_Eta_All_Reconstructed_from_Physical_Primary");  
  TH2D * vzetaMC = hnTracksInBinsMC->Projection(3,1);
  vzetaMC->Sumw2();
  vzetaMC->Write("VZ_Eta_All_Produced");    
  
  TH2D * vzetaEffRec = dynamic_cast<TH2D*>(vzetaRec->Clone("Efficiency_vz_eta_Rec"));
  vzetaEffRec->Divide(vzetaMC);
  vzetaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzetaEffRec->Write();
  TH2D * vzetaEffRecPP = dynamic_cast<TH2D*>(vzetaRecPP->Clone("Efficiency_vz_eta_Rec_pp"));
  vzetaEffRecPP->Divide(vzetaMC);
  vzetaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  vzetaEffRecPP->Write();
  
  TH2D * multptRec = hnTracksInBins->Projection(4,0);
  multptRec->Sumw2();
  multptRec->Write("mult_pT_All_Reconstructed");
  TH2D * multptRecPP = hnTracksInBinsRecPP->Projection(4,0);
  multptRecPP->Sumw2();
  multptRecPP->Write("mult_pT_All_Reconstructed_from_Physical_Primary");  
  TH2D * multptMC = hnTracksInBinsMC->Projection(4,0);
  multptMC->Sumw2();
  multptMC->Write("mult_pT_All_Produced");    
  
  TH2D * multptEffRec = dynamic_cast<TH2D*>(multptRec->Clone("Efficiency_mult_pT_Rec"));
  multptEffRec->Divide(multptMC);
  multptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multptEffRec->Write();
  TH2D * multptEffRecPP = dynamic_cast<TH2D*>(multptRecPP->Clone("Efficiency_mult_pT_Rec_pp"));
  multptEffRecPP->Divide(multptMC);
  multptEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  multptEffRecPP->Write();
  
  
  
  
  //make more meaningful axes:
  //mult: keep for now, unsure how it actually should look.
  Double_t multaxisArray[6] = {0.0,33,66,99,133,165};
  //vz: 4 cm bins from -6 to 6,6-8,8-10 and the same on the other side:
  Double_t  vzaxisArray[8] = {-10.0,-8.0,-6.0,-2.0,2.0,6.0,8.0,10.0};
  //phi need all the bins I can get... 
//   const Double_t * phibins = hnTracksInBins->GetAxis(2)->GetXbins()->GetArray();
  //eta: same
//   const Double_t * etabins = hnTracksInBins->GetAxis(3)->GetXbins()->GetArray();
  //pt: 0.1GeV/c up to 3 GeV/c = 25 bins, 0.5GeV/c up to 5 GeV/c = 4 bins, then 1 5GeV/c bin and one 6GeV/c bin: total 30 bins
  Double_t  pTaxisArray[31] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,5.0,8.0,16.0};
  
  
  
  THnF * RebinRec = Rebin(5,multaxisArray,7,vzaxisArray,hnTracksInBins->GetAxis(2)->GetNbins(),hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins(),hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBins,"hnRebinRec");
  RebinRec->Write();
  THnF * RebinMC = Rebin(5,multaxisArray,7,vzaxisArray,hnTracksInBins->GetAxis(2)->GetNbins(),hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins(),hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBinsMC,"hnRebinMC");
  RebinMC->Write();  
  
  THnF * weights = dynamic_cast<THnF*>(RebinRec->Clone("hnWeights"));
  THnF * weightsnoerr = dynamic_cast<THnF*>(RebinRec->Clone("hnWeight"));
  weights->Clear();
  weightsnoerr->Clear();

  int nerr = 0;
  int nnoerr = 0;
  for(int mult = 1;   mult<=weights->GetAxis(0)->GetNbins();mult++){
    for(int vz = 1;   vz<=weights->GetAxis(1)->GetNbins();vz++){
      for(int phi = 1;phi<=weights->GetAxis(2)->GetNbins();phi++){
	for(int eta=1;eta<=weights->GetAxis(3)->GetNbins();eta++){
	  for(int pT=1;pT<=weights->GetAxis(4)->GetNbins();pT++){
	    Int_t index[5] = {mult,vz,phi,eta,pT};
	    Double_t BinContentRec = RebinRec->GetBinContent(index) ;
	    Double_t BinContentMC = RebinMC->GetBinContent(index) ;
	    Double_t BinContent = 0.0;
	    Double_t BinError = 0.0;
	    if(!BinContentRec<0.5){
	      BinContent = BinContentMC/BinContentRec;
	      BinError = BinContent*BinContent*(RebinRec->GetBinError2(RebinRec->GetBin(index))/(BinContentRec*BinContentRec) + RebinMC->GetBinError2(RebinMC->GetBin(index))/(BinContentMC*BinContentMC));
	    }
	    weights->SetBinContent(index,BinContent);
	    weightsnoerr->SetBinContent(index,BinContent);
	    if(BinError >1.0E-10)weights->SetBinError2(weights->GetBin(index),BinError);
	    
	    if(BinError*100>BinContent) {nerr+=1;}
	    else nnoerr+=1;
	    
	  }
	}
      }
    }
  }
  cout << nerr<<" "<<nnoerr<<endl;
  weights->Write("hnWeights");
  
  
  //Canvases to compare efficiencies in different binnings:
  TCanvas * plotcanvas = new TCanvas("plotcanvas");
  multEffRec->Draw("E");
  TH1D * multRecrb = RebinRec->Projection(0);
  TH1D * multMCrb = RebinMC->Projection(0);
  multRecrb->Sumw2();
  multMCrb->Sumw2();
  TH1D * multEffRecrb = dynamic_cast<TH1D*>(multRecrb->Clone("Efficiency_Mult_Rec_reb"));
  multEffRecrb->Divide(multMCrb);
  multEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRecrb->Write();
  multEffRecrb->Draw("sameE");
  plotcanvas->Write("MultEffComparison");
  plotcanvas->Clear();
  
  vzEffRec->Draw("E");
  TH1D * vzRecrb = RebinRec->Projection(1);
  TH1D * vzMCrb = RebinRec->Projection(1);
  vzRec->Sumw2();
  vzMC->Sumw2();
  TH1D * vzEffRecrb = dynamic_cast<TH1D*>(vzRecrb->Clone("Efficiency_VZ_Rec_rb"));
  vzEffRecrb->Divide(vzMCrb);
  vzEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRecrb->Write();
  vzEffRecrb->Draw("sameE");
  plotcanvas->Write("VZEffComparison");

  
  
  
  //2d projections in correlated variables: 
  TH2D * vzetaweight = weights->Projection(3,1);
  vzetaweight->Write("Weights_Eta_VZ");
  TH2D * multptweight = weights->Projection(4,0);
  multptweight->Write("weights_pT_All_Reconstructed");
  outfile->Close();
  TFile* outfile2 = TFile::Open("LHC11aWeight.root","RECREATE");
  weightsnoerr->Write("hnWeight");
  outfile2->Close();
  
  infile->Close();
// //   delete multEffRec;delete multEffRecPP;
}

void MakeEffHistsPbPb(){
  //open the input file:
  TFile * infile = TFile::Open("AnalysisResults.root","READ");
  THnF * hnTracksInBins;
  THnF * hnTracksInBinsRecPP;
  THnF * hnTracksInBinsMC;
  TList * dir = dynamic_cast<TList*>(infile->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")->Get("ThreePartTrackEfficienciesPbPb_0_16_0_0Coutput1"));
  for(int i=0;i<dir->GetEntries();i++){
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBins")==0)     hnTracksInBins = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsRecPP")==0)hnTracksInBinsRecPP = dynamic_cast<THnF*>(dir->At(i));
    if(TString(dir->At(i)->GetName()).CompareTo("hnTracksinBinsMC")==0)   hnTracksInBinsMC = dynamic_cast<THnF*>(dir->At(i));
  }
  //Axis: 0 - Centrality or multiplicity, 1 - Vertex , 2 - phi , 3 - eta , 4 - pT
  TFile* outfile = TFile::Open("eff.root","RECREATE");
  TAxis * multaxis = hnTracksInBins->GetAxis(0);
  multaxis->SetTitle("Centrality [%]");
  hnTracksInBinsRecPP->GetAxis(0)->SetTitle("Centrality [%]");
  hnTracksInBinsMC->GetAxis(0)->SetTitle("Centrality [%]");
  TAxis * vzAxis = hnTracksInBins->GetAxis(1);
  vzAxis->SetTitle("Z-Vertex [cm]");  
  hnTracksInBinsRecPP->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  hnTracksInBinsMC->GetAxis(1)->SetTitle("Z-Vertex [cm]");
  TAxis * phiaxis = hnTracksInBins->GetAxis(2);
  phiaxis->SetTitle("#Phi [radians]");    
  hnTracksInBinsRecPP->GetAxis(2)->SetTitle("#Phi [radians]");
  hnTracksInBinsMC->GetAxis(2)->SetTitle("#Phi [radians]");
  TAxis * Etaaxis = hnTracksInBins->GetAxis(3);
  Etaaxis->SetTitle("#eta []");  
  hnTracksInBinsRecPP->GetAxis(3)->SetTitle("#eta []");
  hnTracksInBinsMC->GetAxis(3)->SetTitle("#eta []");
  TAxis * pTaxis = hnTracksInBins->GetAxis(4);
  pTaxis->SetTitle("pT [GeV/c]");    
  hnTracksInBinsRecPP->GetAxis(4)->SetTitle("pT [GeV/c]");
  hnTracksInBinsMC->GetAxis(4)->SetTitle("pT [GeV/c]");
  
  
  TH1D * multRec = hnTracksInBins->Projection(0);
  TH1D * multRecPP = hnTracksInBinsRecPP->Projection(0);
  TH1D * multMC = hnTracksInBinsMC->Projection(0);
  multRec->Sumw2();
  multRec->Write("Multiplicity_All_Reconstructed");
  multRecPP->Sumw2();
  multRecPP->Write("Multiplicity_Reconstructed_from_Physical_Primary");  
  multMC->Sumw2();
  multMC->Write("Multiplicity_All_Produced");
  
  TH1D * multEffRec = dynamic_cast<TH1D*>(multRec->Clone("Efficiency_Mult_Rec"));
  multEffRec->Divide(multMC);
  multEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRec->Write();
  TH1D * multEffRecPP = dynamic_cast<TH1D*>(multRecPP->Clone("Efficiency_Mult_Rec_pp"));
  multEffRecPP->Divide(multMC);
  multEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  multEffRecPP->Write();

  TH1D * vzRec = hnTracksInBins->Projection(1);
  TH1D * vzRecPP = hnTracksInBinsRecPP->Projection(1);
  TH1D * vzMC = hnTracksInBinsMC->Projection(1);
  vzRec->Sumw2();
  vzRec->Write("VZ_All_Reconstructed");
  vzRecPP->Sumw2();
  vzRecPP->Write("VZ_Reconstructed_from_Physical_Primary");  
  vzMC->Sumw2();
  vzMC->Write("VZ_All_Produced");
  
  TH1D * vzEffRec = dynamic_cast<TH1D*>(vzRec->Clone("Efficiency_VZ_Rec"));
  vzEffRec->Divide(vzMC);
  vzEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRec->Write();
  TH1D * vzEffRecPP = dynamic_cast<TH1D*>(vzRecPP->Clone("Efficiency_VZ_Rec_pp"));
  vzEffRecPP->Divide(vzMC);
  vzEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  vzEffRecPP->Write();
  
  TH1D * phiRec = hnTracksInBins->Projection(2);
  TH1D * phiRecPP = hnTracksInBinsRecPP->Projection(2);
  TH1D * phiMC = hnTracksInBinsMC->Projection(2);
  phiRec->Sumw2();
  phiRec->Write("phi_All_Reconstructed");
  phiRecPP->Sumw2();
  phiRecPP->Write("phi_Reconstructed_from_Physical_Primary");  
  phiMC->Sumw2();
  phiMC->Write("phi_All_Produced");
  
  TH1D * phiEffRec = dynamic_cast<TH1D*>(phiRec->Clone("Efficiency_phi_Rec"));
  phiEffRec->Divide(phiMC);
  phiEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  phiEffRec->Write();
  TH1D * phiEffRecPP = dynamic_cast<TH1D*>(phiRecPP->Clone("Efficiency_phi_Rec_pp"));
  phiEffRecPP->Divide(phiMC);
  phiEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  phiEffRecPP->Write();
  
  TH1D * etaRec = hnTracksInBins->Projection(3);
  TH1D * etaRecPP = hnTracksInBinsRecPP->Projection(3);
  TH1D * etaMC = hnTracksInBinsMC->Projection(3);
  etaRec->Sumw2();
  etaRec->Write("eta_All_Reconstructed");
  etaRecPP->Sumw2();
  etaRecPP->Write("eta_Reconstructed_from_Physical_Primary");  
  etaMC->Sumw2();
  etaMC->Write("eta_All_Produced");
  
  TH1D * etaEffRec = dynamic_cast<TH1D*>(etaRec->Clone("Efficiency_eta_Rec"));
  etaEffRec->Divide(etaMC);
  etaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  etaEffRec->Write();
  TH1D * etaEffRecPP = dynamic_cast<TH1D*>(etaRecPP->Clone("Efficiency_eta_Rec_pp"));
  etaEffRecPP->Divide(etaMC);
  etaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  etaEffRecPP->Write();

  TH1D * pTRec = hnTracksInBins->Projection(4);
  TH1D * pTRecPP = hnTracksInBinsRecPP->Projection(4);
  TH1D * pTMC = hnTracksInBinsMC->Projection(4);
  pTRec->Sumw2();
  pTRec->Write("pT_All_Reconstructed");
  pTRecPP->Sumw2();
  pTRecPP->Write("pT_Reconstructed_from_Physical_Primary");  
  pTMC->Sumw2();
  pTMC->Write("pT_All_Produced");
  
  TH1D * pTEffRec = dynamic_cast<TH1D*>(pTRec->Clone("Efficiency_pT_Rec"));
  pTEffRec->Divide(pTMC);
  pTEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  pTEffRec->Write();
  TH1D * pTEffRecPP = dynamic_cast<TH1D*>(pTRecPP->Clone("Efficiency_pT_Rec_pp"));
  pTEffRecPP->Divide(pTMC);
  pTEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  pTEffRecPP->Write();
  
  //2d projections in correlated variables:
  TH2D * vzetaRec = hnTracksInBins->Projection(3,1);
  vzetaRec->Sumw2();
  vzetaRec->Write("VZ_Eta_All_Reconstructed");
  TH2D * vzetaRecPP = hnTracksInBinsRecPP->Projection(3,1);
  vzetaRecPP->Sumw2();
  vzetaRecPP->Write("VZ_Eta_All_Reconstructed_from_Physical_Primary");  
  TH2D * vzetaMC = hnTracksInBinsMC->Projection(3,1);
  vzetaMC->Sumw2();
  vzetaMC->Write("VZ_Eta_All_Produced");    
  
  TH2D * vzetaEffRec = dynamic_cast<TH2D*>(vzetaRec->Clone("Efficiency_vz_eta_Rec"));
  vzetaEffRec->Divide(vzetaMC);
  vzetaEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzetaEffRec->Write();
  TH2D * vzetaEffRecPP = dynamic_cast<TH2D*>(vzetaRecPP->Clone("Efficiency_vz_eta_Rec_pp"));
  vzetaEffRecPP->Divide(vzetaMC);
  vzetaEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  vzetaEffRecPP->Write();
  
  TH2D * multptRec = hnTracksInBins->Projection(4,0);
  multptRec->Sumw2();
  multptRec->Write("mult_pT_All_Reconstructed");
  TH2D * multptRecPP = hnTracksInBinsRecPP->Projection(4,0);
  multptRecPP->Sumw2();
  multptRecPP->Write("mult_pT_All_Reconstructed_from_Physical_Primary");  
  TH2D * multptMC = hnTracksInBinsMC->Projection(4,0);
  multptMC->Sumw2();
  multptMC->Write("mult_pT_All_Produced");    
  
  TH2D * multptEffRec = dynamic_cast<TH2D*>(multptRec->Clone("Efficiency_mult_pT_Rec"));
  multptEffRec->Divide(multptMC);
  multptEffRec->SetTitle("Reconstructed tracks divided by produced MC particles");
  multptEffRec->Write();
  TH2D * multptEffRecPP = dynamic_cast<TH2D*>(multptRecPP->Clone("Efficiency_mult_pT_Rec_pp"));
  multptEffRecPP->Divide(multptMC);
  multptEffRecPP->SetTitle("Reconstructed tracks that come from a Physical Primary divided by produced MC particles");  
  multptEffRecPP->Write();
  
  
  
  
  //make more meaningful axes:
  //mult: keep for now, unsure how it actually should look.
  Double_t multaxisArray[12] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,95.0};
  //vz: 4 cm bins from -6 to 6,6-8,8-10 and the same on the other side:
  Double_t  vzaxisArray[8] = {-10.0,-8.0,-6.0,-2.0,2.0,6.0,8.0,10.0};
  //phi need all the bins I can get... 
//   const Double_t * phibins = hnTracksInBins->GetAxis(2)->GetXbins()->GetArray();
  //eta: same
//   const Double_t * etabins = hnTracksInBins->GetAxis(3)->GetXbins()->GetArray();
  //pt: 0.1GeV/c up to 3 GeV/c = 25 bins, 0.5GeV/c up to 5 GeV/c = 4 bins, then 1 5GeV/c bin and one 6GeV/c bin: total 30 bins
  Double_t  pTaxisArray[31] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,5.0,8.0,16.0};
  
  
  
  THnF * RebinRec = Rebin(11,multaxisArray,7,vzaxisArray,hnTracksInBins->GetAxis(2)->GetNbins(),hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins(),hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBins,"hnRebinRec");
  RebinRec->Write();
  THnF * RebinMC = Rebin(11,multaxisArray,7,vzaxisArray,hnTracksInBins->GetAxis(2)->GetNbins(),hnTracksInBins->GetAxis(2)->GetBinLowEdge(1),hnTracksInBins->GetAxis(2)->GetBinUpEdge(hnTracksInBins->GetAxis(2)->GetNbins()),hnTracksInBins->GetAxis(3)->GetNbins(),hnTracksInBins->GetAxis(3)->GetBinLowEdge(1),hnTracksInBins->GetAxis(3)->GetBinUpEdge(hnTracksInBins->GetAxis(3)->GetNbins()),30,pTaxisArray,hnTracksInBinsMC,"hnRebinMC");
  RebinMC->Write();  
  
  THnF * weights = dynamic_cast<THnF*>(RebinRec->Clone("hnWeights"));
  THnF * weightsnoerr = dynamic_cast<THnF*>(RebinRec->Clone("hnWeight"));
  weights->Clear();
  weightsnoerr->Clear();

  int nerr = 0;
  int nnoerr = 0;
  for(int mult = 1;   mult<=weights->GetAxis(0)->GetNbins();mult++){
    for(int vz = 1;   vz<=weights->GetAxis(1)->GetNbins();vz++){
      for(int phi = 1;phi<=weights->GetAxis(2)->GetNbins();phi++){
	for(int eta=1;eta<=weights->GetAxis(3)->GetNbins();eta++){
	  for(int pT=1;pT<=weights->GetAxis(4)->GetNbins();pT++){
	    Int_t index[5] = {mult,vz,phi,eta,pT};
	    Double_t BinContentRec = RebinRec->GetBinContent(index) ;
	    Double_t BinContentMC = RebinMC->GetBinContent(index) ;
	    Double_t BinContent = 0.0;
	    Double_t BinError = 0.0;
	    if(!BinContentRec<0.5){
	      BinContent = BinContentMC/BinContentRec;
	      BinError = BinContent*BinContent*(RebinRec->GetBinError2(RebinRec->GetBin(index))/(BinContentRec*BinContentRec) + RebinMC->GetBinError2(RebinMC->GetBin(index))/(BinContentMC*BinContentMC));
	    }
	    weights->SetBinContent(index,BinContent);
	    weightsnoerr->SetBinContent(index,BinContent);
	    if(BinError >1.0E-10)weights->SetBinError2(weights->GetBin(index),BinError);
	    
	    if(BinError*100>BinContent) {nerr+=1;}
	    else nnoerr+=1;
	    
	  }
	}
      }
    }
  }
  cout << nerr<<" "<<nnoerr<<endl;
  weights->Write("hnWeights");
  
  
  //Canvases to compare efficiencies in different binnings:
  TCanvas * plotcanvas = new TCanvas("plotcanvas");
  multEffRec->Draw("E");
  TH1D * multRecrb = RebinRec->Projection(0);
  TH1D * multMCrb = RebinMC->Projection(0);
  multRecrb->Sumw2();
  multMCrb->Sumw2();
  TH1D * multEffRecrb = dynamic_cast<TH1D*>(multRecrb->Clone("Efficiency_Mult_Rec_reb"));
  multEffRecrb->Divide(multMCrb);
  multEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  multEffRecrb->Write();
  multEffRecrb->Draw("sameE");
  plotcanvas->Write("MultEffComparison");
  plotcanvas->Clear();
  
  vzEffRec->Draw("E");
  TH1D * vzRecrb = RebinRec->Projection(1);
  TH1D * vzMCrb = RebinRec->Projection(1);
  vzRec->Sumw2();
  vzMC->Sumw2();
  TH1D * vzEffRecrb = dynamic_cast<TH1D*>(vzRecrb->Clone("Efficiency_VZ_Rec_rb"));
  vzEffRecrb->Divide(vzMCrb);
  vzEffRecrb->SetTitle("Reconstructed tracks divided by produced MC particles");
  vzEffRecrb->Write();
  vzEffRecrb->Draw("sameE");
  plotcanvas->Write("VZEffComparison");

  
  
  
  //2d projections in correlated variables: 
  TH2D * vzetaweight = weights->Projection(3,1);
  vzetaweight->Write("Weights_Eta_VZ");
  TH2D * multptweight = weights->Projection(4,0);
  multptweight->Write("weights_pT_All_Reconstructed");
  outfile->Close();
  TFile* outfile2 = TFile::Open("LHC10hWeight.root","RECREATE");
  weightsnoerr->Write("hnWeight");
  outfile2->Close();
  
  infile->Close();
// //   delete multEffRec;delete multEffRecPP;
//   
}