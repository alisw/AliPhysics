#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TH3.h"
#include "TF1.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "AliPHOSEmcBadChannelsMap.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#endif
void makeBadModESD(Int_t mod=3){
  //Here we use 3D distributions of pedestals and spectra 
  //produced by the macro scanESD.C 
  //Having these distributions we calculate distributions
  //over number of clusters in "Soft" and "Hard" regions and 
  //exlude cells, havig considerable deviations from avrage
  //Result is stored in TH2S in root file and as CDB root file
  //
  //Currect version is valid for one PHOS module (e.g. cosmic run)
  //Author D.Peressounko 
  
  TFile * f = new TFile("scan.root") ;
  TH3D * hSp = (TH3D*)f->Get("hEsd") ;

  //Boundary of "Soft" and "Hard" regions
  //Soft region used to exclude dead or almost dead cells
  //Hard region used to exclude "singers"
  //All energies in GeV
  Double_t aSoftMin=0.22 ;  
  Double_t aSoftMax=0.5 ;
  Double_t aHardMin=1. ;
  Double_t aHardMax=2. ;

  //Final bad map
  TH2D * hBadMap = new TH2D("hBadMap","Bad Modules map",64,0.,64.,56,0.,56.) ;

  TH2D * hSoft = new TH2D("hSoft","Number of soft clusters",64,0.,64.,56,0.,56.) ;
  TH2D * hHard = new TH2D("hHard","Number of hard clusters",64,0.,64.,56,0.,56.) ;

  //The range of these histos depends on statistics, one probaly will have to increase it  
  TH1D * hAllSoft = new TH1D("hallSoft","Summary of Soft",200,0.,200.) ;
  TH1D * hAllHard = new TH1D("hallHard","Summary of Hard",100,0.,100.) ;

  TAxis * ax = hSp->GetZaxis() ;
  Int_t iSoftMin=ax->FindBin(aSoftMin) ;
  Int_t iSoftMax=ax->FindBin(aSoftMax) ;
  Int_t iHardMin=ax->FindBin(aHardMin) ;
  Int_t iHardMax=ax->FindBin(aHardMax) ;

  //Fill histograms with number of counts per cell
  for(Int_t i=1; i<=hSp->GetXaxis()->GetNbins(); i++){
    for(Int_t j=1; j<=hSp->GetYaxis()->GetNbins(); j++){
      Double_t soft=0 ;
      for(Int_t k=iSoftMin; k<=iSoftMax;k++){
         soft+=hSp->GetBinContent(i,j,k) ;
      }
      Double_t hard=0. ;
      for(Int_t k=iHardMin; k<=iHardMax;k++){
         hard+=hSp->GetBinContent(i,j,k) ;
      }
      hSoft->SetBinContent(i,j,soft) ;
      hHard->SetBinContent(i,j,hard) ;
      hAllSoft->Fill(soft) ;
      hAllHard->Fill(hard) ;
    }
  }


  //=====================================================================================
  //
  //Usually the part below needs some manual intervention to properly select good modules region
  //
  //Fit overall distribution with Gauss and calculate limits of "goodness" 
  //usually +-3-4 sigma for soft and ~4-5 sigma for hard (less statistics)
  TF1 * gaus = new TF1("gaus","[0]*exp(-(x-[1])*(x-[1])/2./[2]/[2])",0.,1000.) ;
  gaus->SetParameters(100.,10.,20.) ;
  gaus->FixParameter(1,0.) ;
  hAllSoft->Fit(gaus,"","",10,80) ;
  Double_t m = gaus->GetParameter(1) ;
  Double_t s = TMath::Abs(gaus->GetParameter(2)) ;
  //one can use either sigma or hardwired parameterization
  //of limits of goodness
  Double_t cutSoftMin= 8. ; //TMath::Max(0.,m-4.*s) ;
  Double_t cutSoftMax= 80. ; //m+4.*s ;
  Double_t softH= gaus->GetParameter(0) ;

  //repeat same procedure for hard
  gaus->SetParameters(100.,2.,10.) ;
  gaus->FixParameter(1,0.) ;
  hAllHard->Fit(gaus,"","",1.,20.) ;
  m = gaus->GetParameter(1) ;
  s = TMath::Abs(gaus->GetParameter(2)) ;
  Double_t cutHardMin= TMath::Max(0.,m-5.*s) ;
  Double_t cutHardMax= m+5.*s ;
  Double_t hardH= gaus->GetParameter(0) ;

  //=====================================================================================

  //Draw  results
  TCanvas * cm = new TCanvas("Soft","Soft",15,22,500,300) ;
  hSoft->Draw("colz") ;

  TCanvas * cr = new TCanvas("Hard","Hard",520,22,500,300) ;
  hHard->Draw("colz") ;

  TCanvas * cma = new TCanvas("AllSoft","Soft Multiplicity",15,330,500,300) ;
  cma->SetLogy() ;
  hAllSoft->Draw() ;
  TLine * l1 = new TLine(cutSoftMin,0.1,cutSoftMin,2.*softH) ;
  l1->SetLineColor(2) ;
  l1->Draw() ;
  TLine * l2 = new TLine(cutSoftMax,0.1,cutSoftMax,2.*softH) ;
  l2->SetLineColor(2) ;
  l2->Draw() ;

  TCanvas * cra = new TCanvas("AllHard","Hard Multiplicity",520,330,500,300) ;
  cra->SetLogy() ;
  hAllHard->Draw() ;
  printf("Min=%f, max=%f \n",cutHardMin,cutHardMax) ;
  TLine * l3 = new TLine(cutHardMin,0.1,cutHardMin,2.*hardH) ;
  l3->SetLineColor(2) ;
  l3->Draw() ;
  TLine * l4 = new TLine(cutHardMax,0.1,cutHardMax,2.*hardH) ;
  l4->SetLineColor(2) ;
  l4->Draw() ;


  //now calculate bad map
  //Exclude everething beyond the limits
  AliPHOSEmcBadChannelsMap badMap;
  for(Int_t i=1; i<=hSp->GetXaxis()->GetNbins(); i++){
    for(Int_t j=1; j<=hSp->GetYaxis()->GetNbins(); j++){
      Double_t soft = hSoft->GetBinContent(i,j) ;
      Double_t hard= hHard->GetBinContent(i,j) ;
      if(soft>=cutSoftMin   && soft <=cutSoftMax &&
        hard>=cutHardMin && hard <=cutHardMax){
        hBadMap->SetBinContent(i,j,1.) ;
        hBadMap->SetBinContent(i,j,1.) ;
      }
      else{
        hBadMap->SetBinContent(i,j,0.) ;
        badMap.SetBadChannel(mod,i,j); //module, col,row
      }
    }
  }
  
  TCanvas * cb = new TCanvas("BadMap","Bad map",15,630,500,300) ;
  hBadMap->DrawClone("colz") ;

  TFile * fout = new TFile("BadMap.root","recreate") ;
  hBadMap->Write() ;
  fout->Close() ;

  //put now result to local CDB
  AliCDBManager *CDB = AliCDBManager::Instance();
  //  CDB->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  CDB->SetDefaultStorage("local://./");
  //  CDB->SetSpecificStorage("local://./","PHOS");
  
  AliCDBMetaData *md= new AliCDBMetaData();
  md->SetResponsible("Dmitri Peressounko");
  md->SetComment("Dead channels for run 1234");
  AliCDBId id("PHOS/Calib/EmcBadChannels",0,999999);
  CDB->Put(&badMap,id, md);
 

}
