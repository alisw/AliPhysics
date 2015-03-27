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
  
  outfile->Close();
  infile->Close();
//   delete multEffRec;delete multEffRecPP;
}