#include "BHHelper.h"


TString genmcfile = "LHC161718genMCs";
TString inputDirectory = "Xi1530MB";
std::vector<std::vector<double>> Multibins = {
    {0,1},
    {1,5},
    {5,10},
    {10,15},
    {15,20},
    {20,30},
    {30,40},
    {40,50},
    {50,60},
    {60,70},
    {70,80},
    {80,90},
    {90,100},
    {0,5},
    {0,10},
    {10,30},
    {30,50},
    {50,70},
    {70,100},
    {0,100}
};
TFile* LoadXi1530Results(TString name, TString runnum);
TObject* LoadXi1530ResultList(TFile* fh, TString clistname);
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum = "");
void GetEventCorrections(){

    auto clist_MC_General = LoadXi1530ResultList(
         genmcfile.Data(), inputDirectory.Data());  // From General Purpose MC

    auto hInvMass_MC_General_Trig =
        BHTHnSparseHelper::Load("htriggered_CINT7", clist_MC_General);
    auto htrue_cent =
        hInvMass_MC_General_Trig.GetTH1("true", 1, {1, -1, -1}); // MC True INEL>0
    auto hReco_cent =
        hInvMass_MC_General_Trig.GetTH1("Reco", 1, {2, -1, -1}); // True INEL>0 && Triggered
    auto hVtx_cent =
        hInvMass_MC_General_Trig.GetTH1("Vtx", 1, {3, -1, -1}); // True INEL>0 && triggered && good vtx  

    std::vector<double> TrigEffi  = {};
    std::vector<double> TrigEffie = {};
    std::vector<double> VertEffi  = {};
    std::vector<double> VertEffie = {};
    for(auto iCentBin = 0; iCentBin < Multibins.size(); iCentBin++){
        auto CentArray = Multibins[iCentBin];
        double sumtrue = htrue_cent->Integral(htrue_cent->GetXaxis()->FindBin(CentArray[0]+0.0001),
                htrue_cent->GetXaxis()->FindBin(CentArray[1]-0.0001));
        double sumreco = hReco_cent->Integral(hReco_cent->GetXaxis()->FindBin(CentArray[0]+0.0001),
                hReco_cent->GetXaxis()->FindBin(CentArray[1]-0.0001));
        double sumvtx = hVtx_cent->Integral(hVtx_cent->GetXaxis()->FindBin(CentArray[0]+0.0001),
                hVtx_cent->GetXaxis()->FindBin(CentArray[1]-0.0001));
        double sumtrig = hTrig_cent->Integral(hTrig_cent->GetXaxis()->FindBin(CentArray[0]+0.0001),
                hTrig_cent->GetXaxis()->FindBin(CentArray[1]-0.0001));
        
        auto triggereffi = sumreco / sumtrue;
        auto triggereffi_e = sqrt(triggereffi*(1-triggereffi)/sumtrue);
        auto vertexeffi = sumvtx/sumreco;
        auto vertexeffi_e = sqrt(vertexeffi*(1-vertexeffi)/sumreco);

        TrigEffi.push_back(triggereffi);
        TrigEffie.push_back(triggereffi_e);
        VertEffi.push_back(vertexeffi);
        VertEffie.push_back(vertexeffi_e);
    }
    std::cout << "MultiBin,Trig.Effi.,Trig.Effi.err,Vert.Effi.,Vert.Effi.err" << std::endl;
    for(auto iCentBin = 0; iCentBin < Multibins.size(); iCentBin++){
        std::cout << Multibins[iCentBin][0] << "-" << Multibins[iCentBin][1] << "," << TrigEffi[iCentBin]  << "," << TrigEffie[iCentBin] << "," << VertEffi[iCentBin] << "," << VertEffie[iCentBin]  << std::endl;
    }
}
TFile* LoadXi1530Results(TString name, TString runnum) {
    if (!runnum.IsNull())
        runnum = "_" + runnum;
    //auto fname = "/alice/home/blim/postprocessing/data/AnalysisResults_" + name + runnum + ".root";
    auto fname = "data_evtcor/AnalysisResults_" + name + runnum + ".root";
    return LoadRoot(fname);
}
//__________________________________________________________
TObject* LoadXi1530ResultList(TFile* fh, TString clistname) {
    auto l = fh->Get(clistname);
    if (!l)
        ErrorExit("No list " + clistname);
    return l;
}
//__________________________________________________________
TObject* LoadXi1530ResultList(TString fname,
                              TString clistname,
                              TString runnum) {
    auto f = LoadXi1530Results(fname, runnum);
    return LoadXi1530ResultList(f, clistname);
}