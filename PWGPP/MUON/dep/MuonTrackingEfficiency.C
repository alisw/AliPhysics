#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TH1I.h>
#include <TH2I.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TString.h>
#include <Riostream.h>
#include <TFile.h>
#include <TList.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <vector>
#include <TGaxis.h>

/*

This macro is meant to be used with the output of the Efficiency Task. It computes the uncorrelated efficiencies per chammber, station and globally. It then looks for correlated dead arrays and corrects the previous calculated efficiency.

Macro created by Lizardo Valencia Palomo

lizardo.valencia.palomo@cern.ch

*/



TH1D* GraphToHist(Bool_t ToBePlotted, const TGraphAsymmErrors* Graph, Int_t Chamber);
std::vector<Double_t> FillVector(TH1D* Chamber, Bool_t Errors);
std::vector<Double_t> CorrelatedDeadArray(std::vector<Double_t> Chj, std::vector<Double_t> Chi, Int_t NumberOfDEs, Int_t SymmetricGroups);
void Drawing(std::vector<Double_t> Chi, std::vector<Double_t> Chj, std::vector<Double_t> ChiErr, std::vector<Double_t> ChjErr, Int_t Bins, Int_t SymGroups, Int_t Station);
Double_t N00(Double_t Average, Double_t CDA, Int_t SimGroup, Int_t SimGroupsPerCh);


void MuonTrackingEfficiency(TString fileName = "./AnalysisResults.root")
{

  TFile *file = new TFile(fileName, "read");
  TDirectoryFile *Directory = (TDirectoryFile*)((file->GetDirectory("MUON_Efficiency")));
  
  std::vector<TH1D*> CorrectedHistos;
  std::vector<TH1D*> CorrectedChiHisto;
  std::vector<TH1D*> CorrectedChjHisto;
  
  TList *listTT = (TList *)Directory->Get("TotalTracksPerChamber;1");
  TList *listTD = (TList *)Directory->Get("TracksDetectedPerChamber;1");
  TList *listSD = (TList *)Directory->Get("SingleDetectedPerChamber;1");

  TH1D *EfficiencyHisto = 0x0; 
  TH1D *TD = 0x0;
  TH1D *TT = 0x0;
  
  for (Int_t Chamber = 1; Chamber < 12; Chamber++) //Loop over all chambers + 1
    {

      TH2F *TD2Dim = (TH2F *) listTD->At(Chamber - 1);
      TH2F *TT2Dim = (TH2F *) listTT->At(Chamber - 1);
      TH2F *SD2Dim = (TH2F *) listSD->At(Chamber - 1);

      Int_t nDEch = 0;
      nDEch = TD2Dim->GetNbinsX();

      TH1D *Ncluster = 0x0;
      TH1D *Corrected = 0x0;
      TH1D *SD = 0x0;

      TH1D *EfficiencyHistoCh11 = 0x0; 
      TH1D *TDch11 = 0x0;
      TH1D *TTch11 = 0x0;

      if (Chamber == 11)
      {

	TDch11 = new TH1D(Form("TDch11 %d",Chamber), "", nDEch, - 0.5, nDEch - 0.5);
        TTch11 = new TH1D(Form("TTch11 %d",Chamber), "", nDEch, - 0.5, nDEch - 0.5);

      }

      TD = new TH1D(Form("TD %d",Chamber), "", nDEch, - 0.5, nDEch - 0.5);
      TT = new TH1D(Form("TT %d",Chamber), "", nDEch, - 0.5, nDEch - 0.5);
      SD = new TH1D(Form("SD %d",Chamber), "", nDEch, - 0.5, nDEch - 0.5);

      std::vector<Double_t> TDentries (nDEch, 0.0);
      std::vector<Double_t> TTentries (nDEch, 0.0);
      std::vector<Double_t> SDentries (nDEch, 0.0);

      for (Int_t ix = 1; ix <= nDEch; ix++)
        {
          for (Int_t iy = 1; iy <= TD2Dim->GetNbinsY(); iy++)
            {

	      TDentries[ix-1] += TD2Dim->GetBinContent(ix,iy); 
	      TTentries[ix-1] += TT2Dim->GetBinContent(ix,iy); 
	      SDentries[ix-1] += SD2Dim->GetBinContent(ix,iy); 
	    
            }

	  TD->SetBinContent(ix,TDentries[ix-1]);
          TD->SetBinError(ix,sqrt(TDentries[ix-1]));
          TT->SetBinContent(ix,TTentries[ix-1]);
          TT->SetBinError(ix,sqrt(TTentries[ix-1]));
	  SD->SetBinContent(ix,SDentries[ix-1]);
          SD->SetBinError(ix,sqrt(SDentries[ix-1]));

	  if (Chamber != 11) continue;

	  TDch11->SetBinContent(ix,TDentries[ix-1]);
          TDch11->SetBinError(ix,sqrt(TDentries[ix-1]));
          TTch11->SetBinContent(ix,TTentries[ix-1]);
          TTch11->SetBinError(ix,sqrt(TTentries[ix-1]));
         	  
        }
      
      TGraphAsymmErrors *EfficiencyGraph = new TGraphAsymmErrors(TD, TT);

      EfficiencyHisto = GraphToHist(kFALSE,EfficiencyGraph,Chamber);

      Ncluster = new TH1D(Form("Ncluster %d",Chamber), "N_{Cluster}", nDEch, - 0.5, nDEch - 0.5);
      Ncluster->Sumw2();
      Ncluster->Add(TD,SD,1,1);

      Corrected = new TH1D(Form("Corrected %d",Chamber), Form("Corrected N_{Cluster} Ch %d",Chamber), nDEch, - 0.5, nDEch - 0.5);
      Corrected->Sumw2();
      Corrected->Divide(Ncluster,EfficiencyHisto,1,1);

      CorrectedHistos.push_back(Corrected);

      if (Chamber != 11) continue;

      TGraphAsymmErrors *EfficiencyGraphCh11 = new TGraphAsymmErrors(TDch11, TTch11);

      EfficiencyHistoCh11 = GraphToHist(kTRUE,EfficiencyGraphCh11,Chamber);

      TCanvas *CanvasEffPerCh = new TCanvas("Efficiency per chambers", "Efficiency per chambers");
      gStyle->SetOptStat(0);
      CanvasEffPerCh->SetFillColor(0);
      CanvasEffPerCh->SetBorderSize(1);
      CanvasEffPerCh->SetFrameLineWidth(2);

      EfficiencyHistoCh11->SetTitle(kFALSE);
      EfficiencyHistoCh11->SetMarkerStyle(20);
      EfficiencyHistoCh11->SetMarkerColor(1);
      EfficiencyHistoCh11->SetMarkerSize(1.3);
      EfficiencyHistoCh11->GetXaxis()->SetTitle("Chamber");
      EfficiencyHistoCh11->GetXaxis()->SetNdivisions(Chamber);
      EfficiencyHistoCh11->GetXaxis()->SetTitleSize(0.05);
      EfficiencyHistoCh11->GetXaxis()->SetTitleOffset(0.91);
      EfficiencyHistoCh11->GetXaxis()->SetLabelSize(0.05);
      EfficiencyHistoCh11->GetYaxis()->SetTitle("Efficiency");
      EfficiencyHistoCh11->GetYaxis()->CenterTitle(kTRUE);
      EfficiencyHistoCh11->GetYaxis()->SetTitleSize(0.06);
      EfficiencyHistoCh11->GetYaxis()->SetTitleOffset(0.79);
      EfficiencyHistoCh11->GetYaxis()->SetLabelSize(0.05);
      EfficiencyHistoCh11->Draw("e");

    }

  /*************************Correlated Dead Arrays check****************************/


  Double_t MissingTracks [5] = {0.0};//Station one is in slot cero!!!
 
  for (Int_t CDAperSt = 0; CDAperSt < 5; CDAperSt++)
    {

      std::vector<Double_t> FilledVectorChi = FillVector(CorrectedHistos[CDAperSt*2], kFALSE);
      std::vector<Double_t> FilledVectorChj = FillVector(CorrectedHistos[(CDAperSt*2)+1], kFALSE);
      std::vector<Double_t> FilledVectorChiErr = FillVector(CorrectedHistos[CDAperSt*2], kTRUE);
      std::vector<Double_t> FilledVectorChjErr = FillVector(CorrectedHistos[(CDAperSt*2)+1], kTRUE);

      Int_t DEs = FilledVectorChi.size();
      Int_t SimGroups = 0;
      Double_t ToFillMissingTracks = 0.0;
      
      if ( CDAperSt < 2 ) SimGroups = static_cast<Int_t>(DEs/2);
      else SimGroups = static_cast<Int_t>(DEs/2) + 1;

      std::vector<Double_t> CDAs = CorrelatedDeadArray(FilledVectorChi, FilledVectorChj, DEs, SimGroups);
      Drawing(FilledVectorChi, FilledVectorChj, FilledVectorChiErr, FilledVectorChjErr, DEs, SimGroups, CDAperSt);

      for (Int_t DEst = 0; DEst < DEs; DEst++)
        {
          if ( CDAs[DEst] != 0.0 ) 
	    {

	      printf ("Looks like there is a CDA in DE %d of Station %d: N00 = %.2f \n", DEst, CDAperSt+1, CDAs[DEst]);
	      ToFillMissingTracks += CDAs[DEst];

	    }
        }

      MissingTracks[CDAperSt] = ToFillMissingTracks;

    }



  /**********************************Efficiency from Efficiency Task************************************/
  
  
      vector<Double_t> chambersEff;
      vector<Double_t> chambersEffErr;
      vector<Double_t> StEfficiency;
      vector<Double_t> StEfficiencyErr;

      chambersEff.push_back(0.0);
      chambersEffErr.push_back(0.0);
      StEfficiency.push_back(0.0);
      StEfficiencyErr.push_back(0.0);

      for (Int_t EfCh = 0; EfCh < 10; EfCh++)
	{

	  chambersEff.push_back(EfficiencyHisto->GetBinContent(EfCh+1));
	  chambersEffErr.push_back(EfficiencyHisto->GetBinError(EfCh+1));

	  cout << "Efficiency chamber " << EfCh+1 << " = " << chambersEff[EfCh+1] << " +- " << chambersEffErr[EfCh+1] << endl;

	}

      for (Int_t EffSt = 0; EffSt < 4; EffSt++)
	{

	  if (EffSt < 3)
	    {

	      StEfficiency.push_back(1.0 - (1.0-chambersEff[(2*EffSt)+1])*(1.0-chambersEff[(2*EffSt)+2]));
	      StEfficiencyErr.push_back(StEfficiency[EffSt+1]*(chambersEffErr[(2*EffSt)+1]/chambersEff[(2*EffSt)+1] + chambersEffErr[(2*EffSt)+2]/chambersEff[(2*EffSt)+2]));

	  cout << "Efficiency station " << EffSt+1 <<" = " <<StEfficiency[EffSt+1] << " +- " <<StEfficiencyErr[EffSt+1]<<endl;

	    }

	  else
	    {

	      StEfficiency.push_back(chambersEff[(2*EffSt)+1] * chambersEff[(2*EffSt)+2] * chambersEff[(2*EffSt)+3] * chambersEff[(2*EffSt)+4] 
              + (1.0 - chambersEff[(2*EffSt)+1]) * chambersEff[(2*EffSt)+2] * chambersEff[(2*EffSt)+3] * chambersEff[(2*EffSt)+4] 
              + chambersEff[(2*EffSt)+1] * (1.0 - chambersEff[(2*EffSt)+2]) * chambersEff[(2*EffSt)+3] * chambersEff[(2*EffSt)+4] 
              + chambersEff[(2*EffSt)+1] * chambersEff[(2*EffSt)+2] * (1.0 - chambersEff[(2*EffSt)+3]) * chambersEff[(2*EffSt)+4] 
              + chambersEff[(2*EffSt)+1] * chambersEff[(2*EffSt)+2] * chambersEff[(2*EffSt)+3] * (1.0 - chambersEff[(2*EffSt)+4]));

	      StEfficiencyErr.push_back(StEfficiency[EffSt+1]*(chambersEffErr[(2*EffSt)+1]/chambersEff[(2*EffSt)+1] 
              + chambersEffErr[(2*EffSt)+2]/chambersEff[(2*EffSt)+2] + chambersEffErr[(2*EffSt)+3]/chambersEff[(2*EffSt)+3] 
              + chambersEffErr[(2*EffSt)+4]/chambersEff[(2*EffSt)+4]));

	      cout << "Efficiency station 45"<<" = " <<StEfficiency[EffSt+1] << " +- " <<StEfficiencyErr[EffSt+1]<<endl;
	      //Notice efficiency from st 45 is in slot 4 of the vector!!!

	    }

       }
      
      StEfficiency.push_back(StEfficiency[1]*StEfficiency[2]*StEfficiency[3]*StEfficiency[4]);
      StEfficiencyErr.push_back(StEfficiencyErr[1]+StEfficiencyErr[2]+StEfficiencyErr[3]+StEfficiencyErr[4]);

      cout << "Total efficiency = " <<StEfficiency[5]<< " +- " <<StEfficiencyErr[5]<< endl;



 /*************************************Correcting efficiencies********************************************/

  
  TH1D *TTcorrected = new TH1D("TTcorrected", "TDcorrected", 10, -0.5, 9.5);

  for (Int_t St = 0; St < 5; St++)
    {

      Double_t NewBinContentChi = TT->GetBinContent((2*St)+1) + MissingTracks[St];
      TTcorrected->SetBinContent((2*St)+1,NewBinContentChi);
      TTcorrected->SetBinContent((2*St)+2,TT->GetBinContent((2*St)+2));

    }

  TGraphAsymmErrors *EfficiencyGraphCorrected = new TGraphAsymmErrors(TD,TTcorrected);

  TH1D *EfficiencyHistoCorrected = GraphToHist(1,EfficiencyGraphCorrected,11);//In this case doesnt matter which value is used in the first slot.

  vector<Double_t> chambersEffCorrected;
  vector<Double_t> chambersEffErrCorrected;

  chambersEffCorrected.push_back(0.0);
  chambersEffErrCorrected.push_back(0.0);

  for (Int_t EfChCorr = 0; EfChCorr < 10; EfChCorr++)
    {

      chambersEffCorrected.push_back(EfficiencyHistoCorrected->GetBinContent(EfChCorr+1));
      chambersEffErrCorrected.push_back(EfficiencyHistoCorrected->GetBinError(EfChCorr+1));

      if ( (chambersEffCorrected[EfChCorr+1]/chambersEff[EfChCorr+1]) == 1 ) continue;

      cout << "Corrected efficiency chamber " << EfChCorr+1 << " = " << chambersEffCorrected[EfChCorr+1] << " +- " << chambersEffErrCorrected[EfChCorr+1] << endl;

    }

  vector<Double_t> StEfficiencyCorrected;
  vector<Double_t> StEfficiencyErrCorrected;

  StEfficiencyCorrected.push_back(0.0);
  StEfficiencyErrCorrected.push_back(0.0);

  //By default CDA in Station 4 is not corrected.
  
  for (Int_t Station = 0; Station < 4; Station++)
    {

      if (MissingTracks[Station] == 0.0) 
	{

	  StEfficiencyCorrected.push_back(StEfficiency[Station+1]);
	  StEfficiencyErrCorrected.push_back(StEfficiencyErr[Station+1]);

	}

      else
	{

	  if (Station < 3)
	    {

	  StEfficiencyCorrected.push_back(StEfficiency[Station+1]*(chambersEffCorrected[(2*Station)+1]/chambersEff[(2*Station)+1]));
          StEfficiencyErrCorrected.push_back(StEfficiencyCorrected[Station+1]*sqrt((StEfficiencyErr[Station+1]/StEfficiency[Station+1])*(StEfficiencyErr[Station+1]/StEfficiency[Station+1])+(chambersEffErrCorrected[(2*Station)+1]/chambersEffCorrected[(2*Station)+1])*(chambersEffErrCorrected[(2*Station)+1]/chambersEffCorrected[(2*Station)+1]) + (chambersEffErr[(2*Station)+1]/chambersEff[(2*Station)+1])*(chambersEffErr[(2*Station)+1]/chambersEff[(2*Station)+1])));

           cout<<"Corrected efficiency station "<<Station+1<<" = "<<StEfficiencyCorrected[Station+1]<<" +- "<<StEfficiencyErrCorrected[Station+1]<<endl;

	   }

	 if (Station == 3)
	   {

	     StEfficiencyCorrected.push_back(StEfficiency[Station+1]*(chambersEffCorrected[(2*Station)+3]/chambersEff[(2*Station)+3]));
             StEfficiencyErrCorrected.push_back(StEfficiencyCorrected[Station+1]*sqrt((StEfficiencyErr[Station+1]/StEfficiency[Station+1])*(StEfficiencyErr[Station+1]/StEfficiency[Station+1])+(chambersEffErrCorrected[(2*Station)+3]/chambersEffCorrected[(2*Station)+3])*(chambersEffErrCorrected[(2*Station)+3]/chambersEffCorrected[(2*Station)+3]) + (chambersEffErr[(2*Station)+3]/chambersEff[(2*Station)+3])*(chambersEffErr[(2*Station)+3]/chambersEff[(2*Station)+3])));

             cout<<"Corrected efficiency station 45"<<" = "<<StEfficiencyCorrected[Station+1]<<" +- "<<StEfficiencyErrCorrected[Station+1]<<endl;

	   }


	}

    }

  StEfficiencyCorrected.push_back(StEfficiencyCorrected[1]*StEfficiencyCorrected[2]*StEfficiencyCorrected[3]*StEfficiencyCorrected[4]);
  StEfficiencyErrCorrected.push_back(StEfficiencyErrCorrected[1]+StEfficiencyErrCorrected[2]+StEfficiencyErrCorrected[3]+StEfficiencyErrCorrected[4]);

  cout<<"Corrected total efficiency"<<" = "<<StEfficiencyCorrected[5]<< " +- " <<StEfficiencyErrCorrected[5]<< endl;

}


/*************************************Functions************************************/


std::vector<Double_t> FillVector(TH1D* Chamber, Bool_t Errors)
{

  Int_t NumOfDEs = Chamber->GetNbinsX();
  std::vector<Double_t> ChVector (NumOfDEs, 0);
  for (Int_t i = 1; i <= NumOfDEs; i++)
	{
	  if (Errors == kFALSE) ChVector.at(i-1) = Chamber->GetBinContent(i);
	  if (Errors == kTRUE) ChVector.at(i-1) = Chamber->GetBinError(i);
	}

  return ChVector;

}

std::vector<Double_t> CorrelatedDeadArray(std::vector<Double_t> Chi, std::vector<Double_t> Chj, Int_t NumberOfDEs, Int_t SymmetricGroups)
{

  std::vector<Double_t> Average (SymmetricGroups, 0.0);
  std::vector<Double_t> CDA (NumberOfDEs, 0.0);
 
  for (Int_t Sim = 0; Sim < SymmetricGroups/2; Sim++)
    {
       if (Sim == 0)
	{

	  if ( (Sim == 0) && (SymmetricGroups/2 == 1) )
	    {

	      Average[Sim] = (Chi[Sim]+Chi[Sim+1]+Chi[Sim+2]+Chi[Sim+3])/4;
              Average[Sim+(SymmetricGroups/2)] = (Chj[Sim]+Chj[Sim+1]+Chj[Sim+2]+Chj[Sim+3])/4;

	      for (Int_t i = 1; i < 5; i++)
                 {

                  if ( (Chi[i-1] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[i-1] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
	            {
	              CDA[i-1] = N00(Average[Sim], Chi[i-1], Sim, SymmetricGroups/2);
	            }

                 }

	    }

	  else
	    {

	  Average[Sim] = (Chi[Sim] + Chi[Sim+(NumberOfDEs/2)])/2;
	  Average[Sim+(SymmetricGroups/2)] = (Chj[Sim] + Chj[Sim+(NumberOfDEs/2)])/2;


         if ( (Chi[Sim] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[Sim] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
            {
	      CDA[Sim] = N00(Average[Sim], Chi[Sim], Sim, SymmetricGroups/2);
	    }
	 if ( (Chi[Sim+(NumberOfDEs/2)] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[Sim+(NumberOfDEs/2)] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
            {
	      CDA[Sim+(NumberOfDEs/2)] = N00(Average[Sim], Chi[Sim+(NumberOfDEs/2)], Sim, SymmetricGroups/2);
	    }
	    
	   }

         }
     
      else
	{

	  Average[Sim] = (Chi[Sim] + Chi[(NumberOfDEs/2)-Sim] + Chi[Sim+(NumberOfDEs/2)] + Chi[NumberOfDEs-Sim])/4;
	  Average[Sim+(SymmetricGroups/2)] = (Chj[Sim] + Chj[(NumberOfDEs/2)-Sim] + Chj[Sim+(NumberOfDEs/2)] + Chj[NumberOfDEs-Sim])/4;

	  if ( (Chi[Sim] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[Sim] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
	     {
	       CDA[Sim] = N00(Average[Sim], Chi[Sim], Sim, SymmetricGroups/2);
	     }
	   if ( (Chi[(NumberOfDEs/2)-Sim] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[(NumberOfDEs/2)-Sim] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
	     {
	       CDA[(NumberOfDEs/2)-Sim] = N00(Average[Sim], Chi[(NumberOfDEs/2)-Sim], Sim, SymmetricGroups/2);
	     }
	   if ( (Chi[Sim+(NumberOfDEs/2)] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[Sim+(NumberOfDEs/2)] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
	     {
	       CDA[Sim+(NumberOfDEs/2)] = N00(Average[Sim], Chi[Sim+(NumberOfDEs/2)], Sim, SymmetricGroups/2);
	     }

	   if ( (Chi[NumberOfDEs-Sim] < (Average[Sim]-(Average[Sim]*0.2))) && (Chj[NumberOfDEs-Sim] < (Average[Sim+(SymmetricGroups/2)]-(Average[Sim+(SymmetricGroups/2)]*0.2))) )
	     {
	       CDA[NumberOfDEs-Sim] = N00(Average[Sim], Chi[NumberOfDEs-Sim], Sim, SymmetricGroups/2);
	       
	     }  

	}
      
    }

  return CDA;

}


Double_t N00(Double_t Average, Double_t CDA, Int_t SimGroup, Int_t SimGroupsPerCh)
{

  Double_t RecomputedAverage = 0.0;
  Double_t ToBeAdded = 0.0;

  if ( (SimGroup == 0) && (SimGroupsPerCh != 1) ) RecomputedAverage = (2*Average) - CDA;
  else RecomputedAverage = ((4*Average) - CDA)/3;

  ToBeAdded = RecomputedAverage - CDA;

  return ToBeAdded;

}


void Drawing(std::vector<Double_t> Chi, std::vector<Double_t> Chj, std::vector<Double_t> ChiErr, std::vector<Double_t> ChjErr, Int_t Bins, Int_t SymGroups, Int_t Station)
{

  TCanvas *Canvas = new TCanvas(Form("Station %d",Station+1), Form("Station %d",Station+1));
  gStyle->SetOptStat(0);
  Canvas->SetFillColor(0);
  if ( (Station+1) < 3) Canvas->Divide(SymGroups/2);
  if ( (Station+1) == 3) Canvas->Divide(Station+1,Station);
  if ( (Station+1) > 3) Canvas->Divide(3,3);
  Canvas->SetBorderSize(1);
  Canvas->SetFrameLineWidth(2);

  for (Int_t Sim = 0; Sim < SymGroups/2; Sim++)
    {

      TH1D *Chamberi = new TH1D(Form("Chi SimGroup %d St %d",Sim,Station+1), "", Bins, - 0.5, Bins - 0.5);
      TGaxis::SetMaxDigits(2);
      Chamberi->SetTitle(kFALSE);
      Chamberi->SetMarkerStyle(20);
      Chamberi->SetMarkerColor(1);
      Chamberi->SetMarkerSize(1.3);
      Chamberi->GetXaxis()->SetTitle("Detection Element");
      Chamberi->GetXaxis()->SetNdivisions(Bins+1);
      Chamberi->GetXaxis()->SetTitleSize(0.05);
      Chamberi->GetXaxis()->SetTitleOffset(0.91);
      Chamberi->GetXaxis()->SetLabelSize(0.05);
      Chamberi->GetYaxis()->SetTitle("Corrected N_{Cluster}");
      Chamberi->GetYaxis()->CenterTitle(kTRUE);
      Chamberi->GetYaxis()->SetTitleSize(0.06);
      Chamberi->GetYaxis()->SetTitleOffset(0.79);
      Chamberi->GetYaxis()->SetLabelSize(0.05);
      TH1D *Chamberj = new TH1D(Form("Chj SimGroup %d St %d",Sim,Station+1), "", Bins, - 0.5, Bins - 0.5);
      Chamberj->SetTitle(kFALSE);
      Chamberj->SetMarkerStyle(20);
      Chamberj->SetMarkerColor(2);
      Chamberj->SetMarkerSize(1.3);

      if (Sim == 0)
	{

	  if ( (Sim == 0) && (SymGroups/2 == 1) )
	    {

	      for (Int_t i = 1; i < 5; i++)
                 {

	           Chamberi->SetBinContent(i, Chi[i-1]);
		   Chamberi->SetBinError(i, ChiErr[i-1]);
	           Chamberj->SetBinContent(i, Chj[i-1]);
		   Chamberj->SetBinError(i, ChjErr[i-1]);

                 }

	    }

	  else
	    {

	      Chamberi->SetBinContent(Sim+1, Chi[Sim]);
	      Chamberi->SetBinContent((Sim+1)+(Bins/2), Chi[Sim+(Bins/2)]);
	      Chamberi->SetBinError(Sim+1, ChiErr[Sim]);
	      Chamberi->SetBinError((Sim+1)+(Bins/2), ChiErr[Sim+(Bins/2)]);
	      Chamberj->SetBinContent(Sim+1, Chj[Sim]);
	      Chamberj->SetBinContent((Sim+1)+(Bins/2), Chj[Sim+(Bins/2)]);
	      Chamberj->SetBinError(Sim+1, ChjErr[Sim]);
	      Chamberj->SetBinError((Sim+1)+(Bins/2), ChjErr[Sim+(Bins/2)]);

	    }

      }
 
      else
	{

	  Chamberi->SetBinContent(Sim+1, Chi[Sim]);
	  Chamberi->SetBinContent((Bins/2)-(Sim-1), Chi[(Bins/2)-Sim]);
	  Chamberi->SetBinContent((Sim+1)+(Bins/2), Chi[Sim+(Bins/2)]);
	  Chamberi->SetBinContent(Bins-(Sim-1), Chi[Bins-Sim]);

	  Chamberi->SetBinError(Sim+1, ChiErr[Sim]);
	  Chamberi->SetBinError((Bins/2)-(Sim-1), ChiErr[(Bins/2)-Sim]);
	  Chamberi->SetBinError((Sim+1)+(Bins/2), ChiErr[Sim+(Bins/2)]);
	  Chamberi->SetBinError(Bins-(Sim-1), ChiErr[Bins-Sim]);

	  Chamberj->SetBinContent(Sim+1, Chj[Sim]);
	  Chamberj->SetBinContent((Bins/2)-(Sim-1), Chj[(Bins/2)-Sim]);
	  Chamberj->SetBinContent((Sim+1)+(Bins/2), Chj[Sim+(Bins/2)]);
	  Chamberj->SetBinContent(Bins-(Sim-1), Chj[Bins-Sim]);

	  Chamberj->SetBinError(Sim+1, ChjErr[Sim]);
	  Chamberj->SetBinError((Bins/2)-(Sim-1), ChjErr[(Bins/2)-Sim]);
	  Chamberj->SetBinError((Sim+1)+(Bins/2), ChjErr[Sim+(Bins/2)]);
	  Chamberj->SetBinError(Bins-(Sim-1), ChjErr[Bins-Sim]);

	}

      Canvas->cd(Sim+1);
      Chamberi->Draw("e");
      Chamberj->Draw("esame");

      TLegend *Legend = new TLegend (0.75, 0.8, 1.0, 1.0);
      Legend->SetBorderSize(0);
      Legend->SetFillColor(0);
      if ( (Station+1) < 3 ) Legend->SetTextSize(0.05);
      if ( (Station+1) > 2 ) Legend->SetTextSize(0.05);
      Legend->AddEntry(Chamberi,Form("Chamber %d",2*Station+1),"p");
      Legend->AddEntry(Chamberj,Form("Chamber %d",2*Station+2),"p");
      Legend->Draw();
      
    }

}


TH1D* GraphToHist(Bool_t ToBePlotted,const TGraphAsymmErrors* Graph,Int_t Chamber)
{

  Int_t nBins = 0;
  if ( Chamber < 5 ) nBins = 4;
  if ( (Chamber == 5) || (Chamber == 6) ) nBins = 18;
  if ( (Chamber > 6) && (Chamber < 11) ) nBins = 26;
  if ( Chamber == 11 ) nBins = 10;
  TArrayF  bins(nBins+1);
  TArrayF  y(nBins);
  TArrayF  ey(nBins);
  Double_t dx = 1;
  Double_t xmin = 10000;
  Double_t xmax = -10000;
  for (Int_t i = 0; i < nBins; i++) { 
    Double_t x   = Graph->GetX()[i];
    Double_t exl = Graph->GetEXlow()[i];
    Double_t exh = Graph->GetEXhigh()[i];
    xmin             = TMath::Min(x-exl, xmin);
    xmax             = TMath::Max(x+exh, xmax);
    bins.fArray[i]   = x-exl;
    bins.fArray[i+1] = x+exh;
    Double_t dxi = exh+exl;
    if (dxi == 0 && i != 0) dxi = bins.fArray[i]-bins.fArray[i-1];
    if (dx == 0) dx  = dxi;
    else if (dxi != dx) dx = 0;
    
    y.fArray[i]  = Graph->GetY()[i];
    ey.fArray[i] = TMath::Max(Graph->GetEYlow()[i],Graph->GetEYhigh()[i]);

  }
  TString name(Graph->GetName());
  TString title(Graph->GetTitle());
  TH1D* h = 0;
  if (ToBePlotted == kTRUE) {
    h = new TH1D(name.Data(), title.Data(), nBins, bins[0]-dx/2 + 1.5, bins[nBins]-dx/2 + 1.5);
  }
  else {
    if ( Chamber < 5 ) h = new TH1D(name.Data(), title.Data(), 4, -0.5, 3.5);
    if ( (Chamber == 5) || (Chamber == 6) ) h = new TH1D(name.Data(), title.Data(), 18,  -0.5, 17.5);
    if ( (Chamber > 6) && (Chamber < 11) ) h = new TH1D(name.Data(), title.Data(), 26, -0.5, 25.5);
    if ( Chamber == 11 ) h = new TH1D(name.Data(), title.Data(), 10, -0.5, 9.5);
  }
  for (Int_t i = 1; i <= nBins; i++) { 
    h->SetBinContent(i, y.fArray[i-1]);
    h->SetBinError(i, ey.fArray[i-1]);
  }
  h->SetMarkerStyle(Graph->GetMarkerStyle());
  h->SetMarkerColor(Graph->GetMarkerColor());
  h->SetMarkerSize(Graph->GetMarkerSize());
  h->SetDirectory(0);
    
  return h;
}
