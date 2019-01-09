#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TList.h>
#include <TStyle.h>
#include <TPad.h>
void checkMiniV0QA()
{
  
    TFile* original = new TFile("AnalysisResults.root");
    TList* listorigin = (TList*)original->Get("LifetimesFiltering_summary");
    TFile* worked = new TFile("MiniV0QAchecks.root");
    TFile risultati("control.root", "RECREATE");
    gStyle->SetOptStat(0);

    TCanvas cv("cv","cv",1600,500);
    cv.Divide(2);
    cv.Print("control.pdf[");
    for (const auto&& objOriginal : *listorigin) 
    {
       TH1* histWorked = (TH1*)worked->Get(objOriginal->GetName());
       if (!histWorked)
         continue;
       TH1* histOriginal = (TH1*)objOriginal;
       cv.cd(1);
       histOriginal->DrawCopy("PLC PMC");
       histWorked->DrawCopy("PLC PMC same");
       histOriginal->Divide(histWorked);
       cv.cd(2);
       histOriginal->DrawCopy("colz");
       cv.Print("control.pdf");
       risultati.cd();
       histOriginal->Write();
    }
    cv.Print("control.pdf]");

}