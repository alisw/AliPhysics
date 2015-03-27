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

void MakeEffHists(){
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
  hnTracksInBinsRecPP->GetAxis(3)->SetTitle("pT [GeV/c]");
  hnTracksInBinsMC->GetAxis(3)->SetTitle("pT [GeV/c]");
  
  
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
  //vz: 4 cm bins from -6 to 6,6-8,8-10 and the same on the other side:
  Double_t  vzaxisArray[8] = {-10.0,-8.0,-6.0,-2.0,2.0,6.0,8.0,10.0};
  //phi need all the bins I can get... 
  //eta: rebin by two
  //pt: 0.1GeV/c up to 3 GeV/c = 25 bins, 0.5GeV/c up to 5 GeV/c = 4 bins, then 1 5GeV/c bin and one 6GeV/c bin: total 31 bins
  Double_t  pTaxisArray[32] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0,4.5,5.0,10.0,16.0};
  
  Int_t newaxes[5] = {multaxis->GetNbins(),7,phiaxis->GetNbins(),Etaaxis->GetNbins()/2,31};
  Double_t xmin[5] = {multaxis->GetBinLowEdge(1),vzaxisArray[0],phiaxis->GetBinLowEdge(1),Etaaxis->GetBinLowEdge(1),pTaxisArray[0]};
  Double_t xmax[5] = {multaxis->GetBinUpEdge(multaxis->GetNbins()),vzaxisArray[7],phiaxis->GetBinUpEdge(phiaxis->GetNbins()),Etaaxis->GetBinUpEdge(Etaaxis->GetNbins()),pTaxisArray[31]};
  THnF * weights =  new THnF("hnWeights","Weights for filling.",5,newaxes,xmin,xmax);
  weights->GetAxis(1)->Set(7,vzaxisArray);
  weights->GetAxis(4)->Set(31,pTaxisArray);

  
  THnF * RebinRec =  new THnF("hnRebinRec","Counts in Reconstruction.",5,newaxes,xmin,xmax);
  RebinRec->GetAxis(1)->Set(7,vzaxisArray);
  RebinRec->GetAxis(4)->Set(31,pTaxisArray);  
  THnF * RebinMC =  new THnF("hnRebinMC","Counts produced.",5,newaxes,xmin,xmax);
  RebinMC->GetAxis(1)->Set(7,vzaxisArray);
  RebinMC->GetAxis(4)->Set(31,pTaxisArray);  
  
  for(int mult = 1;mult<=hnTracksInBins->GetAxis(0)->GetNbins();mult++){
    Int_t multindex = weights->GetAxis(0)->FindBin(hnTracksInBins->GetAxis(0)->GetBinCenter(mult));
    for(int vz = 1;vz<=hnTracksInBins->GetAxis(1)->GetNbins();vz++){
      Int_t vzindex = weights->GetAxis(1)->FindBin(hnTracksInBins->GetAxis(1)->GetBinCenter(vz));
      for(int phi = 1;phi<=hnTracksInBins->GetAxis(2)->GetNbins();phi++){
	Int_t phiindex = weights->GetAxis(2)->FindBin(hnTracksInBins->GetAxis(2)->GetBinCenter(phi));
	for(int eta=1;eta<=hnTracksInBins->GetAxis(3)->GetNbins();eta++){
	  Int_t etaindex = weights->GetAxis(3)->FindBin(hnTracksInBins->GetAxis(3)->GetBinCenter(eta));
	  for(int pT=1;pT<=hnTracksInBins->GetAxis(4)->GetNbins();pT++){
	    Int_t ptindex = weights->GetAxis(4)->FindBin(hnTracksInBins->GetAxis(4)->GetBinCenter(pT));
	    Int_t indexold[5] = {mult,vz,phi,eta,pT};
	    Int_t index[5] = {multindex,vzindex,phiindex,etaindex,ptindex};
	    Double_t BinContentRec = hnTracksInBins->GetBinContent(indexold) ;
	    Double_t BinContentMC = hnTracksInBinsMC->GetBinContent(indexold) ;
	    Double_t BinContent = 0.0;
	    Double_t BinError = 0.0;
	    if(!BinContentRec<0.5){
	      BinContent = BinContentRec;
	      BinError = BinContentRec;
	    }
	    RebinRec->SetBinContent(index,BinContent);
	    if(BinError >1.0E-10)RebinRec->SetBinError2(RebinRec->GetBin(index),BinError);
	    BinContent = 0.0;
	    BinError = 0.0;
	    if(!BinContentMC<0.5){
	      BinContent = BinContentMC;
	      BinError = BinContentMC;
	    }
	    RebinMC->SetBinContent(index,BinContent);
	    if(BinError >1.0E-10)RebinMC->SetBinError2(RebinMC->GetBin(index),BinError);
	    
	  }
	}
      }
    }
  }
  for(int mult = 1;mult<=weights->GetAxis(0)->GetNbins();mult++){
    for(int vz = 1;vz<=weights->GetAxis(1)->GetNbins();vz++){
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
	    if(BinError >1.0E-10)weights->SetBinError2(weights->GetBin(index),BinError);
	  }
	}
      }
    }
  }
  RebinRec->Write();
  RebinMC->Write();
  weights->Write();
  
  outfile->Close();
  infile->Close();
//   delete multEffRec;delete multEffRecPP;

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
}