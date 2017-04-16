/// \file QAplots.C
/// \ingroup CaloTrackCorrMacrosQA
/// \brief Do EMCal QA plots
///
/// Macro to do some EMCAL QA plots (from QA outputs (QAresults.root) 
/// of the QA train (AliAnaCalorimeterQA task) )
///
/// * 1. macro needs AliEMCALGeometry to compute maps
/// * 2. this macro makes the QA plots for
///    * period or single runs
/// * 3. to use it you should have the
///    * QAresults.root files in period/pass/runnumber.root
///    * and if you are running in period mode a runlist.txt
///   file in the same directory: period/pass/runlist.txt
///   were runlist is the list of the runs you want to check
/// * 4. To save the plots you must create prior tu usig the macro in the
///   period/pass/ directory as many subdirectories as number of runs
///   with the name runnumber
///   all gif files will be saved in those subdirectories
/// * 5. The trigger corresponds to the name of the output directory in the root
///   file so "EMC7" for EMC triggers and "default" for AnyInt (Minbias)
///
/// \author Alexis Mas, SUBATECH
/// \author Marie Germain,  SUBATECH


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <Riostream.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TFileMerger.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TGridCollection.h>
#include <TROOT.h>
#include <TGridResult.h>
#include <TClonesArray.h>
#include <TObjString.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include "AliEMCALGeometry.h" 
#endif
using namespace std;

void QAplots(TString fCalorimeter = "EMCAL", TString period = "LHC11h", TString pass = "pass1_HLT", TString trigger= "default")
{
  FILE * pFile;
  TString file = "";
  if (trigger=="default") file = "/scratch/alicehp2/germain/QA/"+period+"/"+ pass + "/runlistMB.txt" ;
  else                    file = "/scratch/alicehp2/germain/QA/"+period+"/"+ pass + "/runlistEMC.txt" ;
  
  pFile = fopen(file.Data(), "r"); //open the text file where include the run list and correct run index
  
  cout << " fcalo: " << fCalorimeter << "; period: " << period << "; pass: " << pass << "  trigger "<<trigger<<  endl; 
    
  Int_t p;
  Int_t q;
  Int_t ncols;
  Int_t nlines = 0 ;
  Int_t RunId[500] ;
  
  Double_t x[500] ;
  Double_t xrun[500] ;
  
  while (1){
    ncols = fscanf(pFile,"%d  %d ",&p,&q);
    if (ncols< 0) break;
    x[nlines]=p;
    RunId[nlines]=q;
    xrun[nlines]=1.*q;
    nlines++;
  }
  fclose(pFile);
  
  const Int_t nRun = nlines ;
  TString base ;
  for(Int_t i = 0 ; i < nRun ; i++) { 
    base = "/scratch/alicehp2/germain/QA/";
    base += period ;
    base += "/";
    base += pass ;
    base += "/";
    base += RunId[i] ;
    TString infile ;
    infile = base + ".root" ;
    TFile *f = TFile::Open(infile);
    DrawOccupancy(RunId[i],period,pass,trigger,f);
    DrawRun(RunId[i],period,pass,trigger,f);
    f->Close();
  }
}

void QAplots(Int_t run, TString period ="LHC11h", TString pass="pass1_HLT", TString trigger= "default")
{
  TString base ;
  base = "/scratch/alicehp2/germain/QA/";
  base += period ;
  base += "/";
  base += pass ;
  base += "/";
  base += run ;
  TString infile ;
  infile = base + ".root" ;
  TFile *f = TFile::Open(infile);
  DrawOccupancy(run,period,pass,trigger,f);
  DrawRun(run,period,pass,trigger,f);
  f->Close();
}

void DrawOccupancy(Int_t run , TString period ="LHC11h", TString pass="pass1_HLT",
                   TString trigger= "default", TFile *f =0x0)
{
  TH2D *hEnergyMap = new TH2D("hEnergyMap","",96,-48,48,120,-0,120);
  TH2D *hOccupancyMap = new TH2D("hOccupancyMap","",96,-48,48,120,-0,120); 
  TH2D *hEnergyMapReal = new TH2D("hEnergyMapReal","",96,-48,48,120,-0,120);
  TH2D *hOccupancyMapReal = new TH2D("hOccupancyMapReal","",96,-48,48,120,-0,120);
  hEnergyMapReal->SetXTitle("eta (bin)");
  hEnergyMapReal->SetYTitle("phi (bin)");
  hEnergyMapReal->SetZTitle("E(GeV)/event");
  hEnergyMapReal->GetZaxis()->SetLabelSize(0.02);
  hOccupancyMapReal->SetXTitle("eta (bin)");
  hOccupancyMapReal->SetZTitle("entries/evt");
  hOccupancyMapReal->GetZaxis()->SetLabelSize(0.02);
  hOccupancyMapReal->SetXTitle("eta (bin)");
  hEnergyMap->SetXTitle("eta (bin)");
  hEnergyMap->SetYTitle("phi (bin)");
  hOccupancyMap->SetXTitle("eta (bin)");
  hOccupancyMap->SetYTitle("phi (bin)");
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  
  AliEMCALGeometry *geom = new AliEMCALGeometry("EMCAL_COMPLETEv1","EMCAL");
  Int_t nSupMod, nModule, nIphi, nIeta, nModulo, iRCU;
  Int_t iphi, ieta,jj,kk;
  Int_t icol, irow;
  Int_t bineta, binphi;
  Int_t realbineta, realbinphi;
  TVector3 vg, gg;
  Double_t eta, phi, glob[3];
  
  
  //LHC11d 
  //Int_t mask[224] = {74, 147, 152, 189, 191, 198, 208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 324, 491, 759, 917, 949, 1376, 1386, 1429, 1487, 1490, 1519, 1967, 2014, 2015, 2112, 2114, 2115, 2116, 2118, 2119, 2120, 2123, 2124, 2125, 2158, 2159, 2326, 2332, 2333, 2350, 2351, 2436, 2448, 2506, 2518, 2534, 2540, 2590, 2640, 2793, 2805, 2815, 2828, 2830, 2869, 2878, 2880, 2881, 2891, 2926, 2985, 3022, 3024, 3070, 3135, 3169, 3263, 3503, 4366, 4560, 4623, 6093, 6331, 6481, 7089, 7113, 7190, 7246, 7425, 7495, 7874,8358, 8811, 9024, 9269, 9302, 9387, 9696, 9697, 9698, 9699, 9701, 9702, 9703, 9704, 9705, 9706, 9707, 9710, 9711, 9748, 9792, 9793, 9794, 9795, 9796, 9797, 9798, 9799, 9800, 9801, 9802, 9803, 9804, 9805, 9806, 9807, 9815, 9819, 9824, 9828, 9829, 9830, 9831, 9832, 9834, 9835, 9836, 9837, 9838, 9839, 9850, 9872, 9874, 9875, 9877, 9878, 9879, 9881, 9882, 9883, 9888, 9890, 9891, 9894, 9896, 9897, 9898, 9899, 9902, 9927, 9938, 9939, 9942, 9943, 9945, 9946, 9947, 9948, 9949, 9950, 9951, 10035, 10073, 10084, 10085, 10086, 10090, 10093, 10112, 10113, 10114, 10115, 10116, 10117, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10126, 10127, 10718, 10723, 10728, 10771, 10831, 11042, 11043, 11090, 11363, 22222222};
  
  
  //LHC11e  
  //Int_t mask[174] = {74, 152, 167, 191, 759, 1059, 1175, 1204, 1288, 1376, 1382, 1386, 1519, 1967, 1968, 2026, 2047, 2112, 2114, 2115, 2116, 2118, 2119, 2120, 2123, 2124, 2125, 2210, 2339, 2350, 2391, 2506, 2540, 2793, 2828, 2869, 2891, 2985, 3135, 3503, 4377, 4817, 5600, 5601, 5602, 5603, 5612, 5613, 5614, 5615, 5648, 5649, 5650, 5651, 5660, 5661, 5662, 5663, 5836, 6104, 6331, 6481, 7089, 7371, 7375, 7425, 7572, 7874, 8358, 9269, 9302, 9389, 9696, 9697, 9698, 9699, 9700, 9701, 9702, 9703, 9705, 9706, 9707, 9708, 9709, 9710, 9711, 9750, 9758, 9792, 9793, 9794, 9795, 9798, 9800, 9801, 9803, 9804, 9815, 9819, 9824, 9825, 9828, 9829, 9830, 9831, 9832, 9833, 9834, 9835, 9836, 9838, 9872, 9874, 9875, 9878, 9882, 9883, 9889, 9890, 9891, 9892, 9893, 9894, 9896, 9897, 9898, 9899, 9900, 9901, 9902, 9903, 9927, 9936, 9937, 9938, 9939, 9940, 9941, 9942, 9943, 9945, 9947, 9948, 9949, 9950, 9951, 10086, 10112, 10113, 10114, 10115, 10116, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10126, 10127, 10134, 10135, 10138, 10143, 10718, 10723, 10771, 11042, 11091, 11363, 2222222};
  
  //LHC11f    
  // Int_t mask[134] = {74, 152, 167, 759, 1204, 1267, 1288, 1376, 1382, 1386, 1424, 1519, 1967, 2026, 2047, 2112, 2114, 2115, 2116, 2118, 2119, 2120, 2123, 2124, 2125, 2506, 2540, 2793, 2828, 2869, 2891, 2985, 3135, 3503, 3785, 4817, 6104, 6331, 6481, 7371, 7375, 7425, 7572, 7874, 8218, 8220, 8222, 9269, 9282, 9302, 9455, 9698, 9699, 9700, 9701, 9702, 9703, 9705, 9706, 9707, 9708, 9709, 9710, 9711, 9748, 9758, 9792, 9793, 9794, 9795, 9796, 9797, 9798, 9799, 9800, 9801, 9803, 9804, 9805,9815, 9828, 9829, 9830, 9831, 9832, 9833, 9834, 9835, 9836, 9838, 9850, 9875, 9891, 9898, 9900, 9927, 9936, 9937, 9938, 9939, 9940, 9941, 9942, 9943, 9944, 9945, 9947, 9948, 9949, 9950, 9951, 10112, 10113, 10114, 10115, 10116, 10118, 10119, 10120, 10121, 10122, 10123, 10124, 10125, 10126, 10127, 10138, 10143, 10363, 10718, 10723, 11091, 11363, 2222222};
  
  //NO MASK
  Int_t mask[1] = {2222222};
  
  TH2F *hCellAmplitude;
  TH1F *hNEvents;
  Int_t Events;
  Int_t n=0;
  TString base = "/scratch/alicehp2/germain/QA/";
  base += period ;
  base += "/";
  base += pass ;
  base += "/";
  base += run ;
  base += "/"; 
  base += trigger;
  // TString infile ;
  //infile = base + ".root" ;
  //*f = TFile::Open(infile);  
  TString direct = "CaloQA_";
  direct += trigger; 
  TDirectoryFile *dir = (TDirectoryFile *)f->Get(direct);
  TList *outputList = (TList*)dir->Get(direct);
  
  hNEvents =(TH1F *)outputList->FindObject("hNEvents");
  Events = hNEvents->GetEntries();
  hCellAmplitude =(TH2F *)outputList->FindObject("EMCAL_hAmpId");
  //hCellAmplitude->Draw();
  //TH2 *hCellAmplitude = (TH2*) gFile->Get("hCellAmplitude");
  Double_t numb =0 ;
  Double_t Eth=0;
  Eth = 5.;
  if ( trigger=="EMC7") Eth=20.;
  
  
  for(Int_t i = 0; i < 11520 ; i++)
  {
    Double_t Esum = 0;
    Double_t Nsum = 0;  
    Double_t EsumH = 0;
    Double_t NsumH = 0;
    Double_t Ratio = 0;
    
    for (Int_t j = 1; j <= hCellAmplitude->GetNbinsX(); j++)
    {
      Double_t E = hCellAmplitude->GetXaxis()->GetBinCenter(j);
      Double_t N = hCellAmplitude->GetBinContent(j, i+1);
      
      if (E < 0.07) continue; 
      //  if (E > 0.07) continue; 
      
      if (E <= Eth) {
        Esum += E*N;
        Nsum += N;
      }
      else { 
        EsumH += E*N;
        NsumH += N;
      }
    }
    
    if(NsumH > 100)  Ratio = Nsum/NsumH ; 
    //  if(Nsum > 20000 && Nsum < 22000 ) cout<<"  "<<i ;
    
    Int_t absId = i;
    if(n!=0) {if(mask[n]<=mask[n-1]) cout<<"not sorted list !!"<<endl;}
    if(i==mask[n]){n++ ; continue; }
    // 	if(Esum/(Double_t) Events > 0.5) cout<<"BAD : "<<i<<endl;
    
    //	hBadCellMap->Fill(1)
    
    
    geom->GetCellIndex(absId,  nSupMod, nModule, nIphi, nIeta);
    geom->GetCellPhiEtaIndexInSModule(nSupMod,nModule,nIphi,nIeta, iphi,ieta);
    
    realbinphi = 120-(nSupMod/2)*24 -iphi -1;
    
    
    if (nSupMod%2==0) realbineta= 48-ieta -1;	
    if (nSupMod%2==1) realbineta= -ieta -1;
    
    // to do as usual (Gustavo) SM0 high left SM9 low right 
    
    binphi = 96 - (nSupMod/2)*24 +iphi;
    if (nSupMod%2==1) bineta=ieta;	 
    /*
     Int_t compt; if(i==0) compt = 1; 
     if(ieta==0){ cout<<compt <<endl ; compt ++;}
     */
    if (nSupMod%2==0) bineta=ieta-48;
    
    hEnergyMapReal->Fill(realbineta,realbinphi,Esum/(Double_t)Events);
    hOccupancyMapReal->Fill(realbineta,realbinphi,Nsum/(Double_t)Events);
    // }
  }
    
  cout<<"N events : "<<Events<<endl;
  
  TString Energy ; Energy = base +  "MapEnergy.pdf";
  TString Entries  ;  Entries = base + "MapEntries.pdf";
  TCanvas *c1 = new TCanvas("c1","energymap",800,800); 
  c1->SetFillColor(0);
  c1->SetGrid();
  c1->cd();  
  c1->SetRightMargin(0.14); 
  TString title = "run ";
  title += run ;
  title += " Summed energy map";
  if(trigger=="EMC7") title += " EMC ";
  hEnergyMapReal->SetTitle(title);
  if(trigger== "EMC7"){
    hEnergyMapReal->SetMinimum(0.001); //FOR Esum EMC triggers
    hEnergyMapReal->SetMaximum(0.05); //FOR Esum
  }else{
    hEnergyMapReal->SetMinimum(0.005); //FOR Esum MB
    hEnergyMapReal->SetMaximum(0.02); //FOR Esum
  }
  hEnergyMapReal->Draw("colz");  
  
  c1->cd();
  c1->SaveAs(Energy);
  
  TCanvas *c2 = new TCanvas("c2","occupancy",800,800); 
  //c2->SetLogz();
  c2->SetFillColor(0);
  c2->SetGrid();
  c2->cd();  
  c2->SetRightMargin(0.14);  
  TString title2 = "run ";
  title2 += run ;
  title2 += " Occupancy map";
  if(trigger=="EMC7") title2 += " EMC ";
  
  hOccupancyMapReal->SetTitle(title2);
  //SCALE USE FOR LHC11h modify fotr pp
  if(trigger== "EMC7"){
    hOccupancyMapReal->SetMinimum(0.01); //FOR Nsum
    hOccupancyMapReal->SetMaximum(0.5);} //FOR Nsum}
  else{
    hOccupancyMapReal->SetMinimum(0.01); //FOR Nsum
    hOccupancyMapReal->SetMaximum(0.1); //FOR Nsum
  } 
  hOccupancyMapReal->Draw("colz");  
  c2->cd();
  c2->SaveAs(Entries);
  c2->SaveAs(Entries2);
}

void  DrawRun(const Int_t run = 167713, TString period ="LHC11h",
              TString pass="pass1_HLT", TString trigger= "default", TFile *f =0x0)
{
  gStyle->SetPalette(1);
  TString base = "/scratch/alicehp2/germain/QA/"; 
  base += period ;
  base += "/";
  base += pass ;
  base += "/";
  TString outfilename ;
  TString outfilename2 ;
  base += run ;
  base += "/" ;
  base += trigger ;
  char legend [100]     ;
  char legend2 [100]     ;
  TString direct = "CaloQA_";
  direct += trigger;
  
  
  TDirectoryFile *dir = (TDirectoryFile *)f->Get(direct);
  TList *outputList = (TList*)dir->Get(direct);
  
  
  if (trigger =="EMC7"){ sprintf(legend,"Run %i EMC ",run);}
  else sprintf(legend,"Run %i ",run);
  
  
  
  hClusterTimeEnergy =(TH2F *)outputList->FindObject("EMCAL_hClusterTimeEnergy");
  
  TString title3 =" Time Vs Energy";
  title3 += legend ;
  hClusterTimeEnergy->SetStats(kFALSE);
  hClusterTimeEnergy->SetTitle(title3);
  
  hClusterVsTrack =(TH2F *)outputList->FindObject("EMCAL_hCaloTrackMNClusters");
  hClusterVsTrack->SetStats(kFALSE);
  TString title3 =" N cluster Vs N track";
  title3 += legend ;
  hClusterVsTrack->SetTitle(title3);
  
  
  hClusterEVsTrack =(TH2F *)outputList->FindObject("EMCAL_hCaloTrackMEClusters");
  hClusterEVsTrack->SetStats(kFALSE);
  hClusterEVsTrack->SetTitle(legend);
  TString title3 =" Sum E cluster Vs N track";
  title3 += legend ;
  hClusterEVsTrack->SetTitle(title3);
  
  
  hClusterEVsV0S =(TH2F *)outputList->FindObject("EMCAL_hCaloV0SEClusters");
  hClusterEVsV0S->SetStats(kFALSE);
  TString title3 =" Sum E cluster Vs V0 signal";
  title3 += legend ;
  hClusterEVsV0S->SetTitle(title3);
  
  
  hNCellsPerClusterMod0 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod0");
  hNCellsPerClusterMod0->SetStats(kFALSE);
  
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 0",run);}
  else {sprintf(legend2,"Run %i Mod 0",run);}
  
  
  hNCellsPerClusterMod0->SetTitle(legend2);
  
  hNCellsPerClusterMod1 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod1");
  hNCellsPerClusterMod1->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 1",run);}
  else {sprintf(legend2,"Run %i Mod 1",run);}
  
  
  hNCellsPerClusterMod1->SetTitle(legend2);
  
  hNCellsPerClusterMod2 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod2");
  hNCellsPerClusterMod2->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 2",run);}
  else {sprintf(legend2,"Run %i Mod 2",run);}
  
  hNCellsPerClusterMod2->SetTitle(legend2);
  
  hNCellsPerClusterMod3 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod3");
  hNCellsPerClusterMod3->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 3",run);}
  else {sprintf(legend2,"Run %i Mod 3",run);}
  
  hNCellsPerClusterMod3->SetTitle(legend2);
  
  
  hNCellsPerClusterMod4 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod4");
  hNCellsPerClusterMod4->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 4",run);}
  else {sprintf(legend2,"Run %i Mod 4",run);}
  
  hNCellsPerClusterMod4->SetTitle(legend2);
  
  hNCellsPerClusterMod5 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod5");
  hNCellsPerClusterMod5->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 5",run);}
  else {sprintf(legend2,"Run %i Mod 5",run);}
  hNCellsPerClusterMod5->SetTitle(legend2);
  
  hNCellsPerClusterMod6 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod6");
  hNCellsPerClusterMod6->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 6",run);}
  else {sprintf(legend2,"Run %i Mod 6",run);}
  hNCellsPerClusterMod6->SetTitle(legend2);
  
  hNCellsPerClusterMod7 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod7");
  hNCellsPerClusterMod7->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 7",run);}
  else {sprintf(legend2,"Run %i Mod 7",run);}
  hNCellsPerClusterMod7->SetTitle(legend2);
  
  hNCellsPerClusterMod8 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod8");
  hNCellsPerClusterMod8->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 8",run);}
  else {sprintf(legend2,"Run %i Mod 8",run);}
  hNCellsPerClusterMod8->SetTitle(legend2);
  
  hNCellsPerClusterMod9 =(TH2F *)outputList->FindObject("EMCAL_hNCellsPerCluster_Mod9");
  hNCellsPerClusterMod9->SetStats(kFALSE);
  if(trigger=="EMC7"){sprintf(legend2,"Run %i EMC Mod 9",run);}
  else {sprintf(legend2,"Run %i Mod 9",run);} 
  hNCellsPerClusterMod9->SetTitle(legend2);
  
  TCanvas  * c31 = new TCanvas("Time Vs E ", "Time Vs E", 600, 600);
  c31->SetLogz();
  c31->SetFillColor(0);
  c31->SetBorderSize(0);
  c31->SetFrameBorderMode(0);     
  // c31->SetOptStat(0); 
  hClusterTimeEnergy->SetYTitle("EMCAL ToF(ns)");
  hClusterTimeEnergy->GetYaxis()->SetTitleOffset(1.2);
  hClusterTimeEnergy->GetYaxis()->SetLabelSize(0.03);
  
  hClusterTimeEnergy->GetZaxis()->SetLabelSize(0.02);
  hClusterTimeEnergy->GetZaxis()->SetTitleSize(0.03);
  hClusterTimeEnergy->GetZaxis()->SetTitleOffset(0.03);
  hClusterTimeEnergy->SetAxisRange(0.,50.);
  hClusterTimeEnergy->Draw("colz");
  c31->cd();
  
  outfilename = base + "TimeRun.pdf" ;
  outfilename2 = base + "TimeRun.gif" ;
  
  c31->SaveAs(outfilename);
  c31->SaveAs(outfilename2);
  
  
  TCanvas  * c32 = new TCanvas(" Cluster Vs Track ","Cluster vs Track", 600, 600);
  c32->SetLogz();
  c32->SetFillColor(0);
  c32->SetBorderSize(0);
  c32->SetFrameBorderMode(0);     
  // c31->SetOptStat(0); 
  
  
  hClusterVsTrack->GetYaxis()->SetTitleOffset(1.2);
  hClusterVsTrack->GetYaxis()->SetLabelSize(0.03);
  
  hClusterVsTrack->GetZaxis()->SetLabelSize(0.02);
  hClusterVsTrack->GetZaxis()->SetTitleSize(0.03);
  hClusterVsTrack->GetZaxis()->SetTitleOffset(0.03);
  hClusterVsTrack->SetAxisRange(0.,1200.);
  hClusterVsTrack->Draw("colz");
  c32->cd();
  
  outfilename = base + "CaloTrackMult.pdf";
  
  outfilename2 = base + "CaloTrackMult.gif";
  
  c32->SaveAs(outfilename);
  
  c32->SaveAs(outfilename2);
  
  
  TCanvas  * c33 = new TCanvas(" Cluster E Vs Track ","Cluster E vs Track", 600, 600);
  c33->SetLogz();
  c33->SetFillColor(0);
  c33->SetBorderSize(0);
  c33->SetFrameBorderMode(0);     
  // c31->SetOptStat(0); 
  
  hClusterEVsTrack->GetYaxis()->SetTitleOffset(1.2);
  hClusterEVsTrack->GetYaxis()->SetLabelSize(0.03);
  
  hClusterEVsTrack->GetZaxis()->SetLabelSize(0.02);
  hClusterEVsTrack->GetZaxis()->SetTitleSize(0.03);
  hClusterEVsTrack->GetZaxis()->SetTitleOffset(0.03);
  hClusterEVsTrack->SetAxisRange(0.,1200.);
  hClusterEVsTrack->Draw("colz");
  c33->cd();
  
  outfilename = base + "CaloETrackMult.pdf";
  outfilename2 = base + "CaloETrackMult.gif";
  
  c33->SaveAs(outfilename);
  
  c33->SaveAs(outfilename2);
  
  TCanvas  * c34 = new TCanvas(" Cluster E Vs V0 ","Cluster E vs V0", 600, 600);
  c34->SetLogz();
  c34->SetFillColor(0);
  c34->SetBorderSize(0);
  c34->SetFrameBorderMode(0);     
  // c31->SetOptStat(0); 
  
  
  hClusterEVsV0S->GetYaxis()->SetTitleOffset(1.2);
  hClusterEVsV0S->GetYaxis()->SetLabelSize(0.03);
  
  hClusterEVsV0S->GetZaxis()->SetLabelSize(0.02);
  hClusterEVsV0S->GetZaxis()->SetTitleSize(0.03);
  hClusterEVsV0S->GetZaxis()->SetTitleOffset(0.03);
  hClusterEVsV0S->GetYaxis()->SetRangeUser(0.,50.);
  hClusterEVsV0S->Draw("colz");
  c34->cd();
  
  
  outfilename = base + "CaloEVsV0s.pdf";
  outfilename2 = base + "CaloEVsV0s.gif";
  
  c34->SaveAs(outfilename);
  c34->SaveAs(outfilename2);
  
  
  
  
  
  TCanvas  * c35 = new TCanvas(" Cells per Cluster 0 ","Cells per Cluster 0", 1000, 1200);
  c35->SetLogz();
  c35->SetFillColor(0);
  c35->SetBorderSize(0);
  c35->SetFrameBorderMode(0);     
  // c31->SetOptStat(0); 
  c35->Divide(3,4);
  
  c35->cd(1);
  gPad->SetLogz();
  hNCellsPerClusterMod0->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod0->GetYaxis()->SetLabelSize(0.03);
  
  hNCellsPerClusterMod0->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod0->GetZaxis()->SetTitleSize(0.4);
  hNCellsPerClusterMod0->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod0->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod0->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod0->Draw("colz");
  c35->cd(2);
  
  gPad->SetLogz();
  hNCellsPerClusterMod1->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod1->GetYaxis()->SetLabelSize(0.03);
  
  hNCellsPerClusterMod1->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod1->GetZaxis()->SetTitleSize(0.2);
  hNCellsPerClusterMod1->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod1->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod1->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod1->Draw("colz");
  c35->cd(3);
  gPad->SetLogz();
  
  
  hNCellsPerClusterMod2->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod2->GetYaxis()->SetLabelSize(0.03);
  
  hNCellsPerClusterMod2->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod2->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod2->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod2->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod2->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod2->Draw("colz");
  c35->cd(4);
  
  gPad->SetLogz();
  hNCellsPerClusterMod3->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod3->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod3->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod3->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod3->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod3->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod3->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod3->Draw("colz");
  
  c35->cd(5);
  gPad->SetLogz();
  
  hNCellsPerClusterMod4->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod4->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod4->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod4->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod4->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod4->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod4->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod4->Draw("colz");
  
  
  
  c35->cd(6);
  gPad->SetLogz();
  
  hNCellsPerClusterMod5->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod5->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod5->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod5->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod5->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod5->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod5->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod5->Draw("colz");
  
  
  c35->cd(7);
  gPad->SetLogz();
  
  hNCellsPerClusterMod6->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod6->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod6->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod6->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod6->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod6->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod6->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod6->Draw("colz");
  
  c35->cd(8);
  gPad->SetLogz();
  
  hNCellsPerClusterMod7->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod7->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod7->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod7->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod7->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod7->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod7->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod7->Draw("colz");
  
  
  c35->cd(9);
  
  gPad->SetLogz();
  hNCellsPerClusterMod8->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod8->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod8->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod8->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod8->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod8->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod8->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod8->Draw("colz");
  
  c35->cd(10);
  gPad->SetLogz();
  
  hNCellsPerClusterMod9->GetYaxis()->SetTitleOffset(1.2);
  hNCellsPerClusterMod9->GetYaxis()->SetLabelSize(0.03);
  hNCellsPerClusterMod9->GetZaxis()->SetLabelSize(0.05);
  hNCellsPerClusterMod9->GetZaxis()->SetTitleSize(0.04);
  hNCellsPerClusterMod9->GetZaxis()->SetTitleOffset(0.03);
  hNCellsPerClusterMod9->GetYaxis()->SetRangeUser(0.,30.);
  hNCellsPerClusterMod9->GetXaxis()->SetRangeUser(0.,20.);
  hNCellsPerClusterMod9->Draw("colz");
  
  outfilename = base + "CellsperCluster.pdf";
  outfilename2 = base + "CellsperCluster.gif";
  
  
  c35->SaveAs(outfilename);
  c35->SaveAs(outfilename2);
}
