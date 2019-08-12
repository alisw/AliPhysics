#define MiniV0efficiency_cxx

#include <TGraph.h>
#include "MiniV0efficiency.h"
#include <TH2.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TF1.h>
#include <AliPDG.h>
#include "MiniV0.h"
#include "MCparticle.h"
#include "HyperTriton2Body.h"
#include "Riostream.h"
using namespace Lifetimes;
using namespace std;

MiniV0efficiency::MiniV0efficiency(TTree *) :
fChain{nullptr},
fOutputFileName{"MiniV0efficiency_Acceptance.root"} {}




void MiniV0efficiency::Begin(TTree * /*tree*/) {
  // The Begin() function is called at the start of the query.
  // When running with PROOF Begin() is only called on the client.
  // The tree argument is deprecated (on PROOF 0 is passed).
  TString option = GetOption();
}



void MiniV0efficiency::SlaveBegin(TTree * /*tree*/) {
  // The SlaveBegin() function is called after the Begin() function.
  // When running with PROOF SlaveBegin() is called on each slave server.
  // The tree argument is deprecated (on PROOF 0 is passed).
  ////////////////////////////////////////
  TString option = GetOption();
  ////
  ////Kaon
  

  fHistV0ptData[0] =
      new TH1D("fHistV0ptDataK", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 4.);
  fHistV0ptMC[0] =
      new TH1D("fHistV0ptMCK", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 4.);

  fHistV0ctData[0] =
      new TH1D("fHistV0ctDataK", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 20.);
  fHistV0ctMC[0] =
      new TH1D("fHistV0ctMCK", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 20.);

  ctAnalysis[0]=new TH2D("ctAnalysisK","K_ct_Analysis",40,-1.,1.,40,0.,20.);
  ctAnalysis[0]->SetXTitle("ctRec-ctGen(#it{cm})");
  ctAnalysis[0]->SetYTitle("ctGen(#it{cm})");

  ////Lambda
  fHistV0ptData[1] =
      new TH1D("fHistV0ptDataL", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 4.);
  fHistV0ptMC[1] =
      new TH1D("fHistV0ptMCL", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 4.);

  fHistV0ctData[1] =
      new TH1D("fHistV0ctDataL", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 40.);
  fHistV0ctMC[1] =
      new TH1D("fHistV0ctMCL", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 40.);

  ctAnalysis[1]=new TH2D("ctAnalysisL","L_ct_Analysis",40,-1.,1.,40,0.,40.);
  ctAnalysis[1]->SetXTitle("ctRec-ctGen(#it{cm})");
  ctAnalysis[1]->SetYTitle("ctGen(#it{cm})");

  ///Hypertriton
  fHistV0ptData[2] =
      new TH1D("fHistV0ptDataH", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 10.);
  fHistV0ptMC[2] =
      new TH1D("fHistV0ptMCH", ";V0 #it{p}_{T} (GeV/#it{c}); Counts", 40, 0., 10.);
  fHistV0ctData[2] =
      new TH1D("fHistV0ctDataH", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 40.);
  fHistV0ctMC[2] =
      new TH1D("fHistV0ctMCH", ";V0 #it{ct} (#it{cm}); Counts", 40, 0., 40.);

  ptAnalysisH = new TH2D("ptAnalysisH","H_pt_Analysis",40,-1.,1.,40,0.,10.);

  ctAnalysis[2]=new TH2D("ctAnalysisH","H_ct_Analysis",40,-2.,2.,40,0.,40.);
  ctAnalysis[2]->SetXTitle("ctRec-ctGen(#it{cm})");
  ctAnalysis[2]->SetYTitle("ctGen(#it{cm})");

  AliPDG::AddParticlesToPdgDataBase();
  ////
}



Bool_t MiniV0efficiency::Process(Long64_t entry) {
  
  // The Process() function is called for each entry in the tree (or possibly
  // keyed object in the case of PROOF) to be processed. The entry argument
  // specifies which entry in the currently loaded tree is to be processed.
  // When processing keyed objects with PROOF, the object is already loaded
  // and is available via the fObject pointer.
  //
  // This function should contain the \"body\" of the analysis. It can contain
  // simple or elaborate selection criteria, run algorithms on the data
  // of the event and typically fill histograms.
  //
  // The processing can be stopped by calling Abort().
  //
  // Use fStatus to set the return value of TTree::Process().
  //
  // The return value is currently not used.
  fReader.SetEntry(entry);
  int p_vec[3]={310,3122,1010010030};
  for(int i=0;i<(static_cast<int>(MCparticles.GetSize()));i++){
    auto& miniMC=MCparticles[i];
    if(miniMC.IsPrimary()){ 
      int part=miniMC.GetPDGcode();
      int ind=miniMC.GetRecoIndex();
      for (int j=0;j<3;j++){
        if(p_vec[j]==part){
          float MCmass=miniMC.GetMass();
          if(miniMC.GetNBodies()==2 || part!=p_vec[2]){
            fHistV0ptMC[j]->Fill(miniMC.GetPt());
            fHistV0ctMC[j]->Fill(MCmass*(miniMC.GetDistOverP()));
          }
          if(ind>=0){
            if(miniMC.GetNBodies()==2 && part==p_vec[2]){
              auto& minihyper= V0Hyper[ind];  
              if(minihyper.GetCandidateInvMass()!=-1){ 
                fHistV0ptData[2]->Fill(minihyper.GetV0pt());
                fHistV0ctData[2]->Fill(MCmass*(minihyper.GetDistOverP()));
                ctAnalysis[2]->Fill(MCmass*(minihyper.GetDistOverP())-MCmass*(miniMC.GetDistOverP()),MCmass*(miniMC.GetDistOverP()));
                ptAnalysisH->Fill(miniMC.GetPt()-minihyper.GetV0pt(),miniMC.GetPt());            

              }
            }  
            else if(part!=p_vec[2]){              
                auto& minidata= V0s[ind];  
                fHistV0ptData[j]->Fill(minidata.GetV0pt());
                fHistV0ctData[j]->Fill(MCmass*(minidata.GetDistOverP()));
                ctAnalysis[j]->Fill(MCmass*(minidata.GetDistOverP())-MCmass*(miniMC.GetDistOverP()),MCmass*(miniMC.GetDistOverP()));
            }
          }  
        }
 
     }      
  }
}

return kTRUE;


} 



void MiniV0efficiency::SlaveTerminate() {
  // The SlaveTerminate() function is called after all entries or objects
  // have been processed. When running with PROOF SlaveTerminate() is called
  // on each slave server.

}



void MiniV0efficiency::Terminate() {
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
  string vec[]={"K","L","H"};
  for(int j=0;j<3;j++){ 
    GetOutputList()->Add(fHistV0ptData[j]);
    GetOutputList()->Add(fHistV0ptMC[j]);
    GetOutputList()->Add(fHistV0ctData[j]);
    GetOutputList()->Add(fHistV0ctMC[j]);
    EffvsPt[j]=(TH1D*)fHistV0ptData[j]->Clone(Form("EffvsPt%s",vec[j].data()));
    EffvsPt[j]->Divide(fHistV0ptMC[j]);
    EffvsPt[j]->SetTitle(Form("%s",vec[j].data()));
    EffvsPt[j]->SetYTitle("MiniV0efficiency x Acceptance");
    GetOutputList()->Add(EffvsPt[j]); 
    Effvsct[j]=(TH1D*)fHistV0ctData[j]->Clone(Form("Effvsct%s",vec[j].data()));
    Effvsct[j]->SetTitle(Form("%s",vec[j].data()));
    Effvsct[j]->Divide(fHistV0ctMC[j]);
    Effvsct[j]->SetYTitle("MiniV0efficiency x Acceptance"); 
    GetOutputList()->Add(Effvsct[j]); 
    GetOutputList()->Add(ctAnalysis[j]);
   }
  GetOutputList()->Add(ptAnalysisH);
  TFile output(Form("results/%s", fOutputFileName.data()),"RECREATE");
  GetOutputList()->Write();
  output.Close();

}



void MiniV0efficiency::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the reader is initialized.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  fReader.SetTree(tree);

}



Bool_t MiniV0efficiency::Notify() {
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}
