/*
  Used in the Expert QA for the tracking:
  HTML pages:
  
  http://aliqatrk.web.cern.ch/aliqatrk/<data|sim>/<year>/period>/<run>
  e.g:
  http://aliqatrk.web.cern.ch/aliqatrk/data/2010/LHC10d/pass2/outputBneg/RawPowerLawFit_Combined.png

  To run the code ses:
  $ALICE_PHYSICS/../PWGPP/QA/detectorQAscripts/TRK.sh
  $ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/makePlots.C


  To test code  compilation:
  gSystem->SetIncludePath("-I$ROOTSYS/include  -I$ALICE_ROOT/include -I$ALICE_ROOT/STEER  -I$ALICE_ROOT/RAW  -I$ALICE_ROOT/STAT -I$ALICE_ROOT/TPC/TPCBase -I$ALICE_ROOT/TPC/TPCCalib");

  .L  $ALICE_PHYSICS/../src/PWGPP/QA/Tracking/ExpertQA/AliHighPtTreeAnalysis.cxx+
  
  To do add unit test before further modifiing

*/

/*
  Original author:
     tbroeker@ikf.uni-frankfurt.de
      
  Rewritten by: (ongoing - started at 11.05.2015)
    marian.ivanov@cern.ch
    j.gronefeld@cern.ch

*/



#include <inttypes.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TPaveStats.h>
#include <TVectorT.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TH3.h>
#include <TGraphErrors.h>
#include <TObjString.h>
#include <TObject.h>
#include "AliESDVertex.h"
#include <TNamed.h>
#include "AliVertex.h"
#include "AliESDtrack.h"
#include "AliESDv0.h"
#include "AliExternalTrackParam.h"
#include <TBits.h>
#include <TParticle.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TTreeStream.h>
#include "TString.h"
#include "TSystem.h"
#include "TMinuit.h"
#include "TCanvas.h"



#define AliHighPtTreeAnalysis_cxx
#include "AliHighPtTreeAnalysis.h"

ClassImp(AliHighPtTreeAnalysis)


AliHighPtTreeAnalysis::AliHighPtTreeAnalysis( ) :
  TObject(),
  fChain(0),                     //!pointer to the analyzed TTree or TChain
  fV0Chain(0),                   //!pointer to the V0 tree    
  OutTree(0),                    // output tree
  fStreamer(0),                  // streamer for the trending output
  //
  // Declaration of leaf types for the highPt tree
  //
  fChunkName(0),          // highPt - chunk name  -  never used -  TO REMOVE    MI!!!
  runNumber(0),           // highPt - run number as double - should be carefully removed MI!!! 
  fRunNumberInt(0),       // highPt - run number as Ingeter
  fMult(0),
  triggerClass(0),
  Bz(0),
  BzInt(0),
  vtxESD(0),
  esdTrack(0),
  particle(0),
  extTPCInnerC(0),
  extInnerParamC(0),
  extInnerParam(0),
  extInnerParamRef(0),
  chi2TPCInnerC(0),
  chi2InnerC(0),
  //
  fPeriodName(0),                // LHC period name
  fMakePeriod(0),                // switch to run per period ExperQA - who is setting ????      
  fHasMC(0),                     // swith has MC information - automatically got from the tree
  fPtCut(1),                     // pt cut lowe rvalue of pt ranges for some histograms - name is misleading as it is not general
  fBfield(0),                    // cache value for B field -- TO FIX - use abasolute field, currently only sign e.g -0.5 T, 0.5 T
  //
  //   
  fApplyCorrections(0),          // switch to apply q/pt correction
  fMakePlots(kTRUE),                 // switch to run  MakePlot and produce figure 
  fMakeV0s(0),                   // switch to make plots/trending for the V0s  (currently harrdwired kTRUE)
  fMakeFitPerfomancePlots(0),    // swith to make performance plot
  // set apply q/pt correction - here we stroe parameters of corrections - factors com from external files set by ::SetApplyCorrections( const char *correctionFile )
  fCorrectionAside(0),           // delta(q/pt) per sector on A side  - combined tracks
  fCorrectionCside(0),           // delta(q/pt) per sector on C side  - combined tracks
  fCorrectionAsideTPCInner(0),   // delta(q/pt) per sector on A side  - TPC only
  fCorrectionCsideTPCInner(0),   // delta(q/pt) per sector on C side  - TPC only
  fCorrectionAsideTPCInnerC(0),  // delta(q/pt) per sector on A side  - TPC constrained 
  fCorrectionCsideTPCInnerC(0),  // delta(q/pt) per sector on C side  - TPC constrained 
  
  fNtracks_TPCLowPt(0),                    // 
  fNtracks_TPCHighPt(0),
  fNtracks_TPCITSLowPt(0),
  fNtracks_TPCITSHighPt(0),
   
  // Declaration of V0 data members
  v0(0),                  
  v0track0(0),
  v0track1(0),
  // Histos
  hPulldcaRTPConly_vs_eta_1pT(0),          
  hPulldcaRcomb_vs_eta_1pT(0),
  hResdcaRTPConly_vs_eta_1pT(0),           
  hResdcaRcomb_vs_eta_1pT(0),
  hphiPull_vs_eta_1pT(0),                  
  hphiRes_vs_eta_1pT(0), 
  hPulldcaR_vs_eta_pT_Aside(0),            
  hPulldcaR_vs_eta_pT_Cside(0),
  hPulldcaRTPCInner_vs_eta_pT_Aside(0),    
  hPulldcaRTPCInner_vs_eta_pT_Cside(0),
  hResdcaR_vs_eta_pT_Aside(0),             
  hResdcaR_vs_eta_pT_Cside(0),
  hResdcaRTPCInner_vs_eta_pT_Aside(0),     
  hResdcaRTPCInner_vs_eta_pT_Cside(0),
  hphiPull_vs_eta_pT_Aside(0),             
  hphiPull_vs_eta_pT_Cside(0),
  hphiRes_vs_eta_pT_Aside(0),              
  hphiRes_vs_eta_pT_Cside(0),
  hPulldcaR_vs_phi_pT_Aside(0),            
  hPulldcaR_vs_phi_pT_Cside(0),
  hPulldcaRTPCInner_vs_phi_pT_Aside(0),    
  hPulldcaRTPCInner_vs_phi_pT_Cside(0),
  hResdcaR_vs_phi_pT_Aside(0),             
  hResdcaR_vs_phi_pT_Cside(0),
  hResdcaRTPCInner_vs_phi_pT_Aside(0),     
  hResdcaRTPCInner_vs_phi_pT_Cside(0),
  hphiPull_vs_phi_pT_Aside(0),             
  hphiPull_vs_phi_pT_Cside(0),
  hphiRes_vs_phi_pT_Aside(0),              
  hphiRes_vs_phi_pT_Cside(0),
  heta_phi_pT(0),
  hphi_vs_eta_pT_cutTPC(0),                
  hphi_vs_eta_pT_cutTPCITS(0),
  

  // histogram for 1/pt shift calculation
  h1pt_vs_eta_phi(0),
  h1ptRes_vs_phi_pT_Aside(0),        
  h1ptRes_vs_phi_pT_Cside(0),        // 1/pT resolution from cov. matrix
  h1ptRes_vs_mult_pT_Aside(0),       
  h1ptRes_vs_mult_pT_Cside(0),           // 1/pT resolution from cov. matrix vs mult.
  h1ptSigma_vs_phi_pT_Aside(0),      
  h1ptSigma_vs_phi_pT_Cside(0),      // sigma 1/pT from cov. matrix
  h1ptSigma_vs_mult_pT_Aside(0),     
  h1ptSigma_vs_mult_pT_Cside(0),      // sigma 1/pT from cov. matrix vs mult.
  //
  // histogram for 1/pt shift calculation for TPCInnerC
  //
  h1ptTPCInnerC_vs_eta_phi(0),
  h1ptResTPCInnerC_vs_phi_pT_Aside(0),    
  h1ptResTPCInnerC_vs_phi_pT_Cside(0),  // 1/pT resolution from cov. matrix TPCInnerC
  h1ptResTPCInnerC_vs_mult_pT_Aside(0),   
  h1ptResTPCInnerC_vs_mult_pT_Cside(0),  // 1/pT resolution from cov. matrix vs mult. TPCInnerC
  h1ptSigmaTPCInnerC_vs_phi_pT_Aside(0),  
  h1ptSigmaTPCInnerC_vs_phi_pT_Cside(0), // 1/pT sigma from cov. matrix TPCInnerC
  h1ptSigmaTPCInnerC_vs_mult_pT_Aside(0), 
  h1ptSigmaTPCInnerC_vs_mult_pT_Cside(0), // 1/pT sigma from cov. matrix vs mult. TPCInnerC
  //
  // histogram for 1/pt shift calculation for TPCInner
  h1ptTPCInner_vs_eta_phi(0), 
  h1ptResTPCInner_vs_phi_pT_Aside(0),    
  h1ptResTPCInner_vs_phi_pT_Cside(0),     // 1/pT resolution from cov. matrix TPCInner
  h1ptResTPCInner_vs_mult_pT_Aside(0),  
  h1ptResTPCInner_vs_mult_pT_Cside(0),   // 1/pT resolution from cov. matrix vs mult. TPCInner
  h1ptSigmaTPCInner_vs_phi_pT_Aside(0),  
  h1ptSigmaTPCInner_vs_phi_pT_Cside(0),   // 1/pT sigma from cov. matrix TPCInner
  h1ptSigmaTPCInner_vs_mult_pT_Aside(0), 
  h1ptSigmaTPCInner_vs_mult_pT_Cside(0),   // 1/pT sigma from cov. matrix vs mult. TPCInner
  //
  // Histogramm for V0s
  //
  hK0sPull_vs_alpha_1pT_pos(0),
  hK0sRes_vs_alpha_1pT_pos(0),
  hK0sPull_vs_alpha_1pT_neg(0),
  hK0sRes_vs_alpha_1pT_neg(0),
  // MC info
  hptPull_vs_eta_pT(0),
  hptRes_vs_eta_pT(0),
  hptPullTPCInnerC_vs_eta_pT(0),
  hptResTPCInnerC_vs_eta_pT(0),
  hptPullTPCInner_vs_eta_pT(0),
  hptResTPCInner_vs_eta_pT(0),

  hptPull_vs_phi_pT(0),
  hptRes_vs_phi_pT(0),
  hptPullTPCInnerC_vs_phi_pT(0),
  hptResTPCInnerC_vs_phi_pT(0),
  hptPullTPCInner_vs_phi_pT(0),
  hptResTPCInner_vs_phi_pT(0)
{
}


AliHighPtTreeAnalysis::~AliHighPtTreeAnalysis()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}



void AliHighPtTreeAnalysis::InitAnalysis( TString file ) 
{
  //
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  //
  if( file.Length()==0 ){
    TFile * f = TFile::Open("AliHighPtTreeAnalysis.root");
    if (!f){
      ::Error("AliHighPtTreeAnalysis::InitAnalysis","You have to specify the Input file");
      exit(1);
    }
    Int_t nbytes= Read("AliHighPtTreeAnalysis");
    if (nbytes<=0){
      ::Error("AliHighPtTreeAnalysis::InitAnalysis","Wrong input file AliHighPtTreeAnalysis.root");
    }
    return;
  }
  
  TTree *tree = NULL;
  TTree *V0tree = NULL;
  TFile *f = TFile::Open(file);
  if (!f) return;
  f->GetObject("highPt",tree);
  f->GetObject("V0s",V0tree);  
  Bool_t highPtTreeOK=kFALSE;
  Bool_t V0sTreeOK=kFALSE;
  if (V0tree) { if(V0tree->GetEntries() > 1 &&
                   V0tree->GetBranchStatus("v0.") &&
                   V0tree->GetBranchStatus("track0.") &&
                   V0tree->GetBranchStatus("track1.") 
                  )   V0sTreeOK = kTRUE;}

  if (tree)   { if(tree->GetEntries() > 1 &&
                   tree->GetBranchStatus("Bz") &&
                   tree->GetBranchStatus("esdTrack.") &&
                   tree->GetBranchStatus("vtxESD.") &&
                   tree->GetBranchStatus("extTPCInnerC.") &&
                   tree->GetBranchStatus("chi2TPCInnerC") &&
                   tree->GetBranchStatus("mult") &&
                   tree->GetBranchStatus("runNumber") 
                  )   highPtTreeOK = kTRUE;}
  //   fMakePlots = kTRUE;
  if(highPtTreeOK)   
  {
    Init(tree);
  } 
  else printf("tree highPt not found or empty!\n");
  if(V0sTreeOK) 
  {
    InitV0tree(V0tree);
  }
  else printf("tree V0s not found or empty!\n");
}





void AliHighPtTreeAnalysis::Loop(){
  //
  // Main function loop fill standard histograms
  // and store object with histograms in file for later analysis
  //
  if(fChain == 0) return;
  BookHistos();
  //
  Long64_t nentriesHPT = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentriesHPT;jentry++) {  
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if(runNumber == 0) runNumber = (Double_t) fRunNumberInt;
    if(Bz == 0) Bz = (Double_t) BzInt;
    if(Bz > 0) fBfield = 1;
    else fBfield = -1;

    FillHistos();
  }

  if(fMakeV0s){
    Double_t massK0 = .497614000000000001;
    Int_t nentriesV0s = fV0Chain->GetEntriesFast();
    nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentriesV0s;jentry++) {  
      Long64_t ientry = LoadV0Tree(jentry);
      if (ientry < 0) break;
      nb = fV0Chain->GetEntry(jentry);   nbytes += nb;
      if( !(v0->GetRr()<7&&v0->GetCausalityP()[0]<0.3&&v0->GetCausalityP()[1]<0.3&&v0track0->fITSncls>1&&v0track1->fITSncls>1)) continue;
      if( v0->AlphaV0()*v0->GetParamP()->GetSigned1Pt() > 0 ){
        hK0sPull_vs_alpha_1pT_pos->Fill( 1/v0->Pt(), v0->AlphaV0(), (v0->GetKFInfo(2,2,0)-massK0)/v0->GetKFInfo(2,2,1) );
        hK0sRes_vs_alpha_1pT_pos ->Fill( 1/v0->Pt(), v0->AlphaV0(), (v0->GetKFInfo(2,2,0)-massK0)                      );
      }
      if( v0->AlphaV0()*v0->GetParamP()->GetSigned1Pt() < 0 ){
        hK0sPull_vs_alpha_1pT_neg->Fill( 1/v0->Pt(), v0->AlphaV0(), (v0->GetKFInfo(2,2,0)-massK0)/v0->GetKFInfo(2,2,1) );
        hK0sRes_vs_alpha_1pT_neg ->Fill( 1/v0->Pt(), v0->AlphaV0(), (v0->GetKFInfo(2,2,0)-massK0)                      );
      }
    }
  }
  TFile * fout = TFile::Open("AliHighPtTreeAnalysis.root","recreate");
  fout->cd();
  this->Write("AliHighPtTreeAnalysis");
  fout->Close();  
  Terminate();
}

void AliHighPtTreeAnalysis::FillHistos(){
  //
  //
  //


  if(esdTrack->IsOn(0x0040)&&esdTrack->GetTPCclusters(0)>0&&esdTrack->GetTPCClusterInfo(2,1)>120.&&(esdTrack->GetTPCNclsF()>0&&(esdTrack->GetTPCClusterInfo(2,1)/esdTrack->GetTPCNclsF())>0.8)&&(esdTrack->GetTPCchi2()/esdTrack->GetTPCclusters(0)<4.0)&&(esdTrack->GetTPCnclsS()/esdTrack->GetTPCclusters(0)<0.4)&&abs(esdTrack->fdTPC)<3&&abs(esdTrack->fzTPC)<3 ){
    hphi_vs_eta_pT_cutTPC->Fill( esdTrack->Phi(), esdTrack->Eta(), (1./abs(esdTrack->GetSigned1Pt())) );
    if(1./abs(esdTrack->GetSigned1Pt()) < 1.) fNtracks_TPCLowPt++;
    if(1./abs(esdTrack->GetSigned1Pt()) > 4.) fNtracks_TPCHighPt++;

    if(esdTrack->IsOn(0x0004)&&(esdTrack->HasPointOnITSLayer(0)||esdTrack->HasPointOnITSLayer(1))&&esdTrack->fITSncls>0&&sqrt(esdTrack->fITSchi2/esdTrack->fITSncls)<6){
      hphi_vs_eta_pT_cutTPCITS->Fill( esdTrack->Phi(), esdTrack->Eta(), (1./abs(esdTrack->GetSigned1Pt())) );
      if(1./abs(esdTrack->GetSigned1Pt()) < 1.) fNtracks_TPCITSLowPt++;
      if(1./abs(esdTrack->GetSigned1Pt()) > 4.) fNtracks_TPCITSHighPt++;
    }
  }
  if( BaseCut() ){ 
    if(fApplyCorrections){
      h1pt_vs_eta_phi->Fill( esdTrack->GetSigned1Pt() + qoverptCorr(esdTrack->Eta(), esdTrack->Phi(), 0), esdTrack->Eta(), esdTrack->Phi() ); 
      h1ptTPCInner_vs_eta_phi->Fill( esdTrack->GetTPCInnerParam()->GetSigned1Pt() + qoverptCorr(esdTrack->GetTPCInnerParam()->Eta(),esdTrack->GetTPCInnerParam()->Phi(),1), esdTrack->GetTPCInnerParam()->Eta(), esdTrack->GetTPCInnerParam()->Phi() );
      h1ptTPCInnerC_vs_eta_phi->Fill( extTPCInnerC->GetSigned1Pt() + qoverptCorr(extTPCInnerC->Eta(),extTPCInnerC->Phi(),2), extTPCInnerC->Eta(), extTPCInnerC->Phi() );
    }
    h1pt_vs_eta_phi->Fill( esdTrack->GetSigned1Pt(), esdTrack->Eta(), esdTrack->Phi() ); 
    h1ptTPCInner_vs_eta_phi->Fill( esdTrack->GetTPCInnerParam()->GetSigned1Pt(), esdTrack->GetTPCInnerParam()->Eta(), esdTrack->GetTPCInnerParam()->Phi() );
    h1ptTPCInnerC_vs_eta_phi->Fill( extTPCInnerC->GetSigned1Pt(), extTPCInnerC->Eta(), extTPCInnerC->Phi() );


    heta_phi_pT->Fill( esdTrack->Eta(), esdTrack->Phi(), (1./abs(esdTrack->GetSigned1Pt())) );

    Double_t pullDCA, pullDCATPCinner;

    if( esdTrack->fCdd == 0 ) 
      pullDCA = 0;
    else
      pullDCA = (esdTrack->fD/sqrt(esdTrack->fCdd));

    if( esdTrack->fCddTPC == 0 ) 
      pullDCATPCinner = 0;
    else
      pullDCATPCinner = (esdTrack->fdTPC/sqrt(esdTrack->fCddTPC));

    hPulldcaRTPConly_vs_eta_1pT->Fill( abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt()), esdTrack->GetTPCInnerParam()->Eta(), pullDCATPCinner );
    hPulldcaRcomb_vs_eta_1pT   ->Fill( abs(esdTrack->GetSigned1Pt()), esdTrack->Eta(), pullDCA );
    hResdcaRTPConly_vs_eta_1pT ->Fill( abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt()), esdTrack->GetTPCInnerParam()->Eta(), esdTrack->fdTPC );
    hResdcaRcomb_vs_eta_1pT    ->Fill( abs(esdTrack->GetSigned1Pt()), esdTrack->Eta(), esdTrack->fD );

    hphiPull_vs_eta_1pT        ->Fill( abs(esdTrack->GetSigned1Pt()), esdTrack->Eta(), (esdTrack->fP[2]-extTPCInnerC->fP[2])/sqrt(esdTrack->fC[5]+extTPCInnerC->fC[5]) );                  
    hphiRes_vs_eta_1pT         ->Fill( abs(esdTrack->GetSigned1Pt()), esdTrack->Eta(), (esdTrack->fP[2]-extTPCInnerC->fP[2]) );  

    if(esdTrack->Eta() > 0){
      hPulldcaR_vs_eta_pT_Aside          ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Eta(), pullDCA );
      hResdcaR_vs_phi_pT_Aside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), esdTrack->fD );
      hResdcaR_vs_eta_pT_Aside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Eta(), esdTrack->fD );
      hphiPull_vs_phi_pT_Aside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), (esdTrack->fP[2]-extTPCInnerC->fP[2])/sqrt(esdTrack->fC[5]+extTPCInnerC->fC[5]) );
      hphiRes_vs_phi_pT_Aside            ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), (esdTrack->fP[2]-extTPCInnerC->fP[2]) );
      hPulldcaR_vs_phi_pT_Aside          ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), pullDCA );  

      hPulldcaRTPCInner_vs_eta_pT_Aside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Eta(), pullDCATPCinner );
      hResdcaRTPCInner_vs_eta_pT_Aside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Eta(), esdTrack->fdTPC );  
      hPulldcaRTPCInner_vs_phi_pT_Aside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), pullDCATPCinner );
      hResdcaRTPCInner_vs_phi_pT_Aside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), esdTrack->fdTPC );

      h1ptRes_vs_phi_pT_Aside            ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), 1./abs(esdTrack->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptRes_vs_mult_pT_Aside           ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , fMult                                     , 1./abs(esdTrack->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptSigma_vs_phi_pT_Aside          ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , esdTrack->Phi()                          , TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptSigma_vs_mult_pT_Aside         ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , fMult                                     , TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptResTPCInnerC_vs_phi_pT_Aside   ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), extTPCInnerC->Phi()                      , 1./abs(extTPCInnerC->GetSigned1Pt())*TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptResTPCInnerC_vs_mult_pT_Aside  ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), fMult                                     , 1./abs(extTPCInnerC->GetSigned1Pt())*TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptSigmaTPCInnerC_vs_phi_pT_Aside ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), extTPCInnerC->Phi()                      , TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptSigmaTPCInnerC_vs_mult_pT_Aside->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), fMult                                     , TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptResTPCInner_vs_phi_pT_Aside    ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), 1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptResTPCInner_vs_mult_pT_Aside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), fMult                               , 1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptSigmaTPCInner_vs_phi_pT_Aside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptSigmaTPCInner_vs_mult_pT_Aside ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), fMult                               , TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
    }
    if(esdTrack->Eta() < 0){
      hPulldcaR_vs_eta_pT_Cside          ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Eta(), pullDCA );
      hResdcaR_vs_phi_pT_Cside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), esdTrack->fD );
      hResdcaR_vs_eta_pT_Cside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Eta(), esdTrack->fD );
      hphiPull_vs_phi_pT_Cside           ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), (esdTrack->fP[2]-extTPCInnerC->fP[2])/sqrt(esdTrack->fC[5]+extTPCInnerC->fC[5]) );
      hphiRes_vs_phi_pT_Cside            ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), (esdTrack->fP[2]-extTPCInnerC->fP[2]) );
      hPulldcaR_vs_phi_pT_Cside          ->Fill( (1./abs(esdTrack->GetSigned1Pt())), esdTrack->Phi(), pullDCA );

      hPulldcaRTPCInner_vs_eta_pT_Cside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Eta(), pullDCATPCinner );
      hResdcaRTPCInner_vs_eta_pT_Cside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Eta(), esdTrack->fdTPC );  
      hPulldcaRTPCInner_vs_phi_pT_Cside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), pullDCATPCinner );
      hResdcaRTPCInner_vs_phi_pT_Cside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), esdTrack->fdTPC );

      h1ptRes_vs_phi_pT_Cside            ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , esdTrack->Phi()                          , 1./abs(esdTrack->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptRes_vs_mult_pT_Cside           ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , fMult                                     , 1./abs(esdTrack->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptSigma_vs_phi_pT_Cside          ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , esdTrack->Phi()                          , TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptSigma_vs_mult_pT_Cside         ->Fill( (1./abs(esdTrack->GetSigned1Pt()))    , fMult                                     , TMath::Sqrt(esdTrack->GetSigma1Pt2()) );
      h1ptResTPCInnerC_vs_phi_pT_Cside   ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), extTPCInnerC->Phi()                      , 1./abs(extTPCInnerC->GetSigned1Pt())*TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptResTPCInnerC_vs_mult_pT_Cside  ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), fMult                                     , 1./abs(extTPCInnerC->GetSigned1Pt())*TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptSigmaTPCInnerC_vs_phi_pT_Cside ->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), extTPCInnerC->Phi()                      , TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptSigmaTPCInnerC_vs_mult_pT_Cside->Fill( (1./abs(extTPCInnerC->GetSigned1Pt())), fMult                                     , TMath::Sqrt(extTPCInnerC->GetSigma1Pt2()) );
      h1ptResTPCInner_vs_phi_pT_Cside    ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), 1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptResTPCInner_vs_mult_pT_Cside   ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), fMult                               , 1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())*TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptSigmaTPCInner_vs_phi_pT_Cside  ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), esdTrack->GetTPCInnerParam()->Phi(), TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) );
      h1ptSigmaTPCInner_vs_mult_pT_Cside ->Fill( (1./abs(esdTrack->GetTPCInnerParam()->GetSigned1Pt())), fMult                               , TMath::Sqrt(esdTrack->GetTPCInnerParam()->GetSigma1Pt2()) ); 
    }
  }
}

Bool_t AliHighPtTreeAnalysis::BaseCut(){
  if( (esdTrack->GetTPCclusters(0)>0)&&(esdTrack->GetTPCClusterInfo(2,1)>120.)&&(esdTrack->GetTPCNclsF()>0&&(esdTrack->GetTPCClusterInfo(2,1)/esdTrack->GetTPCNclsF())>0.8)&&(esdTrack->GetTPCchi2()/esdTrack->GetTPCclusters(0)<4.0)&&(esdTrack->GetTPCnclsS()/esdTrack->GetTPCclusters(0)<0.4)&&(esdTrack->IsOn(0x0004))&&(esdTrack->HasPointOnITSLayer(0)||esdTrack->HasPointOnITSLayer(1))&&(esdTrack->fITSncls>0&&esdTrack->fITSchi2>0&&sqrt(esdTrack->fITSchi2/esdTrack->fITSncls)<6)&&sqrt(chi2TPCInnerC)<6&&abs(esdTrack->fZ)<2.0&&abs(esdTrack->fD)<(0.018+0.035*abs(esdTrack->fP[4])) )
    return kTRUE;
  else return kFALSE;
}



void AliHighPtTreeAnalysis::BookHistos(){
  // Book the Histograms which are filled in the loop

  // Define binning
  // log pt axis
  Double_t xmin = 0.05; 
  Double_t xmax = 100.; 
  Int_t nbins = 30; 

  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t *xbins =  new Double_t[nbins+1];
  //  Double_t xbins[nbins + 1];
  xbins[0] = xmin;

  for (Int_t i=0;i<=nbins;i++) {
    xbins[i] = xmin + TMath::Power(10,logxmin+i*binwidth);
  }
  // eta bins
  Double_t ymin = -1.; 
  Double_t ymax = 1.; 
  Double_t nybins = 20; 
  Double_t ybins[21];
  for (Int_t i=0;i<=nybins;i++) {
    ybins[i] = ymin + i*(ymax-ymin)/nybins;
  }
  // phi bins
  Double_t minPhi = 0.; 
  Double_t maxPhi = 6.5; 
  Double_t nbinsPhi = 100; 
  Double_t binsPhi[101];
  for (Int_t i=0;i<=nbinsPhi;i++) {
    binsPhi[i] = minPhi + i*(maxPhi-minPhi)/nbinsPhi;
  }
  // pT pull bins
  Double_t zmin = -15.; 
  Double_t zmax = 15.; 
  Double_t nzbins = 175; 
  Double_t zbins[176];
  for (Int_t i=0;i<=nzbins;i++) {
    zbins[i] = zmin + i*(zmax-zmin)/nzbins;
  }
  // pT resol bins
/*
  Double_t zzmin = -1.; 
  Double_t zzmax = 1.; 
  Double_t nzzbins = 200; 
  Double_t zzbins[201];
  for (Int_t i=0;i<=nzzbins;i++) {
    zzbins[i] = zzmin + i*(zzmax-zzmin)/nzzbins;
  }
*/
  // DCAr global resol bins
  Double_t zzzmin = -0.5; 
  Double_t zzzmax = 0.5; 
  Double_t nzzzbins = 200; 
  Double_t zzzbins[201];
  for (Int_t i=0;i<=nzzzbins;i++) {
    zzzbins[i] = zzzmin + i*(zzzmax-zzzmin)/nzzzbins;
  }
  // DCAr TPC resol bins
  Double_t zzzzmin = -5.0; 
  Double_t zzzzmax = 5.0; 
  Double_t nzzzzbins = 200; 
  Double_t zzzzbins[201];
  for (Int_t i=0;i<=nzzzzbins;i++) {
    zzzzbins[i] = zzzzmin + i*(zzzzmax-zzzzmin)/nzzzzbins;
  }
  // Dphi bins
  Double_t minDPhi = -0.05;
  Double_t maxDPhi = 0.05;
  Double_t nbinsDPhi = 100;
  Double_t binsDPhi[101];
  for (Int_t i=0;i<=nbinsDPhi;i++) {
    binsDPhi[i] = minDPhi + i*(maxDPhi-minDPhi)/nbinsDPhi;
  }
  // 1pT resol cov matrix bins
  Double_t min1PtRes = 0.;
  Double_t max1PtRes = 0.2;
  Double_t nbins1PtRes = 200;
  Double_t bins1PtRes[201];
  for (Int_t i=0;i<=nbins1PtRes;i++) {
    bins1PtRes[i] = min1PtRes + i*(max1PtRes-min1PtRes)/nbins1PtRes;
  }
  // mult bins
  Double_t minMult = 0.;
  Double_t maxMult = 4000;
  Double_t nbinsMult = 100;
  Double_t binsMult[101];
  for (Int_t i=0;i<=nbinsMult;i++) {
    binsMult[i] = minMult + i*(maxMult-minMult)/nbinsMult;
  }
  // 1pt sigma bins
  Double_t min1PtSigma = 0.;
  Double_t max1PtSigma = 0.1;
  Double_t nbins1PtSigma = 200;
  Double_t bins1PtSigma[201];
  for (Int_t i=0;i<=nbins1PtSigma;i++) {
    bins1PtSigma[i] = min1PtSigma + i*(max1PtSigma-min1PtSigma)/nbins1PtSigma;
  }


  // Book Histograms
  heta_phi_pT                         = new TH3D("heta_phi_pT", "eta vs phi vs pt", nybins,ybins,nbinsPhi,binsPhi,nbins,xbins);
  hphi_vs_eta_pT_cutTPC               = new TH3D("hphi_vs_eta_pT_cutTPC", "phi vs eta vs pt (TPC cut)",nbinsPhi,binsPhi,nybins,ybins,nbins,xbins);
  hphi_vs_eta_pT_cutTPCITS            = new TH3D("hphi_vs_eta_pT_cutTPCITS", "phi vs eta vs pt (TPC&ITS cut)",nbinsPhi,binsPhi,nybins,ybins,nbins,xbins);
  h1pt_vs_eta_phi                     = new TH3F("h1pt_vs_eta_phi", "1/pt vs eta vs phi",200,-1.,1.,10,-1.,1.,18,0,6.28);
  h1ptTPCInner_vs_eta_phi             = new TH3F("h1ptTPCInner_vs_eta_phi", "1/pt TPCInner vs eta vs phi",200,-1.,1.,10,-1.,1.,18,0,6.28);
  h1ptTPCInnerC_vs_eta_phi            = new TH3F("h1ptTPCInnerC_vs_eta_phi", "1/pt TPCInnerC vs eta vs phi",200,-1.,1.,10,-1.,1.,18,0,6.28);
  hPulldcaR_vs_eta_pT_Aside           = new TH3D("hPulldcaR_vs_eta_pT_Aside","dcaR/sigma vs eta vs pT for global tracks",nbins,xbins,nybins,ybins,nzbins,zbins);
  hPulldcaR_vs_eta_pT_Cside           = new TH3D("hPulldcaR_vs_eta_pT_Cside","dcaR/sigma vs eta vs pT for global tracks",nbins,xbins,nybins,ybins,nzbins,zbins);
  hPulldcaR_vs_phi_pT_Aside           = new TH3D("hPulldcaR_vs_phi_pT_Aside","dcaR/sigma vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hPulldcaR_vs_phi_pT_Cside           = new TH3D("hPulldcaR_vs_phi_pT_Cside","dcaR/sigma vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hPulldcaRTPCInner_vs_eta_pT_Aside   = new TH3D("hPulldcaRTPCInner_vs_eta_pT_Aside","dcaR/sigma vs eta vs pT for tpc inner tracks",nbins,xbins,nybins,ybins,nzbins,zbins);
  hPulldcaRTPCInner_vs_eta_pT_Cside   = new TH3D("hPulldcaRTPCInner_vs_eta_pT_Cside","dcaR/sigma vs eta vs pT for tpc inner tracks",nbins,xbins,nybins,ybins,nzbins,zbins);
  hPulldcaRTPCInner_vs_phi_pT_Aside   = new TH3D("hPulldcaRTPCInner_vs_phi_pT_Aside","dcaR/sigma vs phi vs pT for tpc inner tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hPulldcaRTPCInner_vs_phi_pT_Cside   = new TH3D("hPulldcaRTPCInner_vs_phi_pT_Cside","dcaR/sigma vs phi vs pT for tpc inner tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hResdcaR_vs_phi_pT_Aside            = new TH3D("hResdcaR_vs_phi_pT_Aside","dcaR vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nzzzbins,zzzbins);
  hResdcaR_vs_phi_pT_Cside            = new TH3D("hResdcaR_vs_phi_pT_Cside","dcaR vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nzzzbins,zzzbins);
  hResdcaR_vs_eta_pT_Aside            = new TH3D("hResdcaR_vs_eta_pT_Aside","dcaR vs eta vs pT for global tracks",nbins,xbins,nybins,ybins,nzzzbins,zzzbins);
  hResdcaR_vs_eta_pT_Cside            = new TH3D("hResdcaR_vs_eta_pT_Cside","dcaR vs eta vs pT for global tracks",nbins,xbins,nybins,ybins,nzzzbins,zzzbins);
  hResdcaRTPCInner_vs_eta_pT_Aside    = new TH3D("hResdcaRTPCInner_vs_eta_pT_Aside","dcaR vs eta vs pT for tpc inner tracks",nbins,xbins,nybins,ybins,nzzzzbins,zzzzbins);
  hResdcaRTPCInner_vs_eta_pT_Cside    = new TH3D("hResdcaRTPCInner_vs_eta_pT_Cside","dcaR vs eta vs pT for tpc inner tracks",nbins,xbins,nybins,ybins,nzzzzbins,zzzzbins);
  hResdcaRTPCInner_vs_phi_pT_Aside    = new TH3D("hResdcaRTPCInner_vs_phi_pT_Aside","dcaR vs phi vs pT for tpc inner tracks",nbins,xbins,nbinsPhi,binsPhi,nzzzzbins,zzzzbins);
  hResdcaRTPCInner_vs_phi_pT_Cside    = new TH3D("hResdcaRTPCInner_vs_phi_pT_Cside","dcaR vs phi vs pT for tpc inner tracks",nbins,xbins,nbinsPhi,binsPhi,nzzzzbins,zzzzbins);
  hphiPull_vs_phi_pT_Aside            = new TH3D("hphiPull_vs_phi_pT_Aside","#Delta(phi)/sigma vs phi vs pT between global and TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hphiPull_vs_phi_pT_Cside            = new TH3D("hphiPull_vs_phi_pT_Cside","#Delta(phi)/sigma vs phi vs pT between global and TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nzbins,zbins);
  hphiRes_vs_phi_pT_Aside             = new TH3D("hphiRes_vs_phi_pT_Aside","#Delta(phi) vs phi vs pT between global and TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbinsDPhi,binsDPhi);
  hphiRes_vs_phi_pT_Cside             = new TH3D("hphiRes_vs_phi_pT_Cside","#Delta(phi) vs phi vs pT between global and TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbinsDPhi,binsDPhi);
  h1ptRes_vs_phi_pT_Aside             = new TH3D("h1ptRes_vs_phi_pT_Aside","1/pT resolution from cov. matrix vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptRes_vs_phi_pT_Cside             = new TH3D("h1ptRes_vs_phi_pT_Cside","1/pT resolution from cov. matrix vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptRes_vs_mult_pT_Aside            = new TH3D("h1ptRes_vs_mult_pT_Aside","1/pT resolution from cov. matrix vs mult vs pT for global tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptRes_vs_mult_pT_Cside            = new TH3D("h1ptRes_vs_mult_pT_Cside","1/pT resolution from cov. matrix vs mult vs pT for global tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptSigma_vs_phi_pT_Aside           = new TH3D("h1ptSigma_vs_phi_pT_Aside","1/pT sigma from cov. matrix vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigma_vs_phi_pT_Cside           = new TH3D("h1ptSigma_vs_phi_pT_Cside","1/pT sigma from cov. matrix vs phi vs pT for global tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigma_vs_mult_pT_Aside          = new TH3D("h1ptSigma_vs_mult_pT_Aside","1/pT sigma from cov. matrix vs mult vs pT for global tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);
  h1ptSigma_vs_mult_pT_Cside          = new TH3D("h1ptSigma_vs_mult_pT_Cside","1/pT sigma from cov. matrix vs mult vs pT for global tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);
  h1ptResTPCInnerC_vs_phi_pT_Aside    = new TH3D("h1ptResTPCInnerC_vs_phi_pT_Aside","1/pT resolution from cov. matrix vs phi vs pT for TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptResTPCInnerC_vs_phi_pT_Cside    = new TH3D("h1ptResTPCInnerC_vs_phi_pT_Cside","1/pT resolution from cov. matrix vs phi vs pT for TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptResTPCInnerC_vs_mult_pT_Aside   = new TH3D("h1ptResTPCInnerC_vs_mult_pT_Aside","1/pT resolution from cov. matrix vs mult vs pT for TPC constrained tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptResTPCInnerC_vs_mult_pT_Cside   = new TH3D("h1ptResTPCInnerC_vs_mult_pT_Cside","1/pT resolution from cov. matrix vs mult vs pT for TPC constrained tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptSigmaTPCInnerC_vs_phi_pT_Aside  = new TH3D("h1ptSigmaTPCInnerC_vs_phi_pT_Aside","1/pT sigma from cov. matrix vs phi vs pT for TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInnerC_vs_phi_pT_Cside  = new TH3D("h1ptSigmaTPCInnerC_vs_phi_pT_Cside","1/pT sigma from cov. matrix vs phi vs pT for TPC constrained tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInnerC_vs_mult_pT_Aside = new TH3D("h1ptSigmaTPCInnerC_vs_mult_pT_Aside","1/pT sigma from cov. matrix vs mult vs pT for TPC constrained tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInnerC_vs_mult_pT_Cside = new TH3D("h1ptSigmaTPCInnerC_vs_mult_pT_Cside","1/pT sigma from cov. matrix vs mult vs pT for TPC constrained tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);
  h1ptResTPCInner_vs_phi_pT_Aside     = new TH3D("h1ptResTPCInner_vs_phi_pT_Aside","1/pT resolution from cov. matrix vs phi vs pT for TPC inner tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptResTPCInner_vs_phi_pT_Cside     = new TH3D("h1ptResTPCInner_vs_phi_pT_Cside","1/pT resolution from cov. matrix vs phi vs pT for TPC inner tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtRes,bins1PtRes);
  h1ptResTPCInner_vs_mult_pT_Aside    = new TH3D("h1ptResTPCInner_vs_mult_pT_Aside","1/pT resolution from cov. matrix vs mult vs pT for TPC inner tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptResTPCInner_vs_mult_pT_Cside    = new TH3D("h1ptResTPCInner_vs_mult_pT_Cside","1/pT resolution from cov. matrix vs mult vs pT for TPC inner tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtRes,bins1PtRes);
  h1ptSigmaTPCInner_vs_phi_pT_Aside   = new TH3D("h1ptSigmaTPCInner_vs_phi_pT_Aside","1/pT sigma from cov. matrix vs phi vs pT for TPC inner tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInner_vs_phi_pT_Cside   = new TH3D("h1ptSigmaTPCInner_vs_phi_pT_Cside","1/pT sigma from cov. matrix vs phi vs pT for TPC inner tracks",nbins,xbins,nbinsPhi,binsPhi,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInner_vs_mult_pT_Aside  = new TH3D("h1ptSigmaTPCInner_vs_mult_pT_Aside","1/pT sigma from cov. matrix vs mult vs pT for TPC inner tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);
  h1ptSigmaTPCInner_vs_mult_pT_Cside  = new TH3D("h1ptSigmaTPCInner_vs_mult_pT_Cside","1/pT sigma from cov. matrix vs mult vs pT for TPC inner tracks",nbins,xbins,nbinsMult,binsMult,nbins1PtSigma,bins1PtSigma);

  hK0sPull_vs_alpha_1pT_pos           = new TH3D("hK0sPull_vs_alpha_1pT_pos","(InvMass_K0s-InvMass_K0s_PDG)/#sigma vs alpha vs pT",20,0,1.,10,-0.8,0.8,100,-10.  ,10.  );
  hK0sRes_vs_alpha_1pT_pos            = new TH3D("hK0sRes_vs_alpha_1pT_pos" ,"(InvMass_K0s-InvMass_K0s_PDG) vs alpha vs pT"       ,20,0,1.,10,-0.8,0.8,100,- 0.05, 0.05);
  hK0sPull_vs_alpha_1pT_neg           = new TH3D("hK0sPull_vs_alpha_1pT_neg","(InvMass_K0s-InvMass_K0s_PDG)/#sigma vs alpha vs pT",20,0,1.,10,-0.8,0.8,100,-10.  ,10.  );
  hK0sRes_vs_alpha_1pT_neg            = new TH3D("hK0sRes_vs_alpha_1pT_neg" ,"(InvMass_K0s-InvMass_K0s_PDG) vs alpha vs pT"       ,20,0,1.,10,-0.8,0.8,100,- 0.05, 0.05);

  hPulldcaRTPConly_vs_eta_1pT         = new TH3D("hPulldcaRTPConly_vs_eta_1pT","dcaR/sigma vs eta vs 1pT (tpc only tracking)" ,20,0,1.,10,-0.8,0.8,175,-15.,15.);
  hPulldcaRcomb_vs_eta_1pT            = new TH3D("hPulldcaRcomb_vs_eta_1pT","dcaR/sigma vs eta vs 1pT (combined tracking)"    ,20,0,1.,10,-0.8,0.8,175,-15.,15. );
  hResdcaRTPConly_vs_eta_1pT          = new TH3D("hResdcaRTPConly_vs_eta_1pT","dcaR vs eta vs 1pT (tpc only tracking)"        ,20,0,1.,10,-0.8,0.8,200, -5. ,5. );
  hResdcaRcomb_vs_eta_1pT             = new TH3D("hResdcaRcomb_vs_eta_1pT","dcaR vs eta vs 1pT (combined tracking)"           ,20,0,1.,10,-0.8,0.8,200, -0.5,0.5);

  hphiPull_vs_eta_1pT                 = new TH3D("hphiPull_vs_eta_1pT","#Delta(phi)/sigma vs eta vs 1/pT between global and TPC constrained tracks",20,0,1.,10,-0.8,0.8,175,-15.  ,15.  );
  hphiRes_vs_eta_1pT                  = new TH3D("hphiRes_vs_eta_1pT","#Delta(phi) vs eta vs 1/pT between global and TPC constrained tracks"       ,20,0,1.,10,-0.8,0.8,100,- 0.05, 0.05);
}

void AliHighPtTreeAnalysis::Terminate(){

  TFile *file = 0;
  if(fApplyCorrections){
    //gSystem->Exec(Form("if [ ! -d ./corrected/runs/%d ] ; then mkdir -p corrected/runs/%d ; fi",(Int_t) runNumber,(Int_t) runNumber));
    //gSystem->cd(Form("corrected/runs/%d",(Int_t) runNumber));

    if( fBfield == 1 )
      file = new TFile("genericHistos_Bpos.root", "RECREATE");
    else 
      file = new TFile("genericHistos_Bneg.root", "RECREATE");
  }
  else{
    //gSystem->Exec(Form("if [ ! -d ./runs/%d ] ; then mkdir -p runs/%d ; fi",(Int_t) runNumber,(Int_t) runNumber));
    //gSystem->cd(Form("runs/%d",(Int_t) runNumber));

    if( fBfield == 1 )
      file = new TFile("genericHistos_Bpos.root", "RECREATE");
    else 
      file = new TFile("genericHistos_Bneg.root", "RECREATE");
  }

  if(hPulldcaRTPConly_vs_eta_1pT) hPulldcaRTPConly_vs_eta_1pT   ->Write();      if(hPulldcaRcomb_vs_eta_1pT) hPulldcaRcomb_vs_eta_1pT   ->Write();
  if(hResdcaRTPConly_vs_eta_1pT)  hResdcaRTPConly_vs_eta_1pT    ->Write();      if(hResdcaRcomb_vs_eta_1pT)  hResdcaRcomb_vs_eta_1pT   ->Write();

  if(hphiPull_vs_eta_1pT) hphiPull_vs_eta_1pT   ->Write();
  if(hphiRes_vs_eta_1pT)  hphiRes_vs_eta_1pT    ->Write();

  if(hPulldcaR_vs_eta_pT_Aside)         hPulldcaR_vs_eta_pT_Aside         ->Write();  if(hPulldcaR_vs_eta_pT_Cside)         hPulldcaR_vs_eta_pT_Cside ->Write();
  if(hPulldcaRTPCInner_vs_eta_pT_Aside) hPulldcaRTPCInner_vs_eta_pT_Aside ->Write();  if(hPulldcaRTPCInner_vs_eta_pT_Cside) hPulldcaRTPCInner_vs_eta_pT_Cside ->Write();
  if(hResdcaR_vs_eta_pT_Aside)          hResdcaR_vs_eta_pT_Aside          ->Write();  if(hResdcaR_vs_eta_pT_Cside)          hResdcaR_vs_eta_pT_Cside ->Write();
  if(hResdcaRTPCInner_vs_eta_pT_Aside)  hResdcaRTPCInner_vs_eta_pT_Aside  ->Write();  if(hResdcaRTPCInner_vs_eta_pT_Cside)  hResdcaRTPCInner_vs_eta_pT_Cside ->Write();
  if(hphiPull_vs_eta_pT_Aside)          hphiPull_vs_eta_pT_Aside          ->Write();  if(hphiPull_vs_eta_pT_Cside)          hphiPull_vs_eta_pT_Cside ->Write();
  if(hphiRes_vs_eta_pT_Aside)           hphiRes_vs_eta_pT_Aside           ->Write();  if(hphiRes_vs_eta_pT_Cside)           hphiRes_vs_eta_pT_Cside ->Write();
  if(hPulldcaR_vs_phi_pT_Aside)         hPulldcaR_vs_phi_pT_Aside         ->Write();  if(hPulldcaR_vs_phi_pT_Cside)         hPulldcaR_vs_phi_pT_Cside ->Write();
  if(hPulldcaRTPCInner_vs_phi_pT_Aside) hPulldcaRTPCInner_vs_phi_pT_Aside ->Write();  if(hPulldcaRTPCInner_vs_phi_pT_Cside) hPulldcaRTPCInner_vs_phi_pT_Cside ->Write();
  if(hResdcaR_vs_phi_pT_Aside)          hResdcaR_vs_phi_pT_Aside          ->Write();  if(hResdcaR_vs_phi_pT_Cside)          hResdcaR_vs_phi_pT_Cside ->Write();
  if(hResdcaRTPCInner_vs_phi_pT_Aside)  hResdcaRTPCInner_vs_phi_pT_Aside  ->Write();  if(hResdcaRTPCInner_vs_phi_pT_Cside)  hResdcaRTPCInner_vs_phi_pT_Cside ->Write();
  if(hphiPull_vs_phi_pT_Aside)          hphiPull_vs_phi_pT_Aside          ->Write();  if(hphiPull_vs_phi_pT_Cside)          hphiPull_vs_phi_pT_Cside ->Write();
  if(hphiRes_vs_phi_pT_Aside)           hphiRes_vs_phi_pT_Aside           ->Write();  if(hphiRes_vs_phi_pT_Cside)           hphiRes_vs_phi_pT_Cside ->Write();
  if(heta_phi_pT)                       heta_phi_pT                       ->Write();
  if(hphi_vs_eta_pT_cutTPC)             hphi_vs_eta_pT_cutTPC             ->Write();  if(hphi_vs_eta_pT_cutTPCITS)          hphi_vs_eta_pT_cutTPCITS ->Write();

  // histogram for 1/pt shift calculation
  if(h1pt_vs_eta_phi)            h1pt_vs_eta_phi ->Write();
  if(h1ptRes_vs_phi_pT_Aside)    h1ptRes_vs_phi_pT_Aside    ->Write();  if(h1ptRes_vs_phi_pT_Cside)    h1ptRes_vs_phi_pT_Cside ->Write();        // 1/pT resolution from cov. matrix
  if(h1ptRes_vs_mult_pT_Aside)   h1ptRes_vs_mult_pT_Aside   ->Write();  if(h1ptRes_vs_mult_pT_Cside)   h1ptRes_vs_mult_pT_Cside ->Write();           // 1/pT resolution from cov. matrix vs mult.
  if(h1ptSigma_vs_phi_pT_Aside)  h1ptSigma_vs_phi_pT_Aside  ->Write();  if(h1ptSigma_vs_phi_pT_Cside)  h1ptSigma_vs_phi_pT_Cside ->Write();      // sigma 1/pT from cov. matrix
  if(h1ptSigma_vs_mult_pT_Aside) h1ptSigma_vs_mult_pT_Aside ->Write();  if(h1ptSigma_vs_mult_pT_Cside) h1ptSigma_vs_mult_pT_Cside ->Write();      // sigma 1/pT from cov. matrix vs mult.

  // histogram for 1/pt shift calculation for TPCInnerC
  if(h1ptTPCInnerC_vs_eta_phi)            h1ptTPCInnerC_vs_eta_phi            ->Write();
  if(h1ptResTPCInnerC_vs_phi_pT_Aside)    h1ptResTPCInnerC_vs_phi_pT_Aside    ->Write();  if(h1ptResTPCInnerC_vs_phi_pT_Cside)    h1ptResTPCInnerC_vs_phi_pT_Cside    ->Write();  // 1/pT resolution from cov. matrix TPCInnerC
  if(h1ptResTPCInnerC_vs_mult_pT_Aside)   h1ptResTPCInnerC_vs_mult_pT_Aside   ->Write();  if(h1ptResTPCInnerC_vs_mult_pT_Cside)   h1ptResTPCInnerC_vs_mult_pT_Cside   ->Write();  // 1/pT resolution from cov. matrix vs mult. TPCInnerC
  if(h1ptSigmaTPCInnerC_vs_phi_pT_Aside)  h1ptSigmaTPCInnerC_vs_phi_pT_Aside  ->Write();  if(h1ptSigmaTPCInnerC_vs_phi_pT_Cside)  h1ptSigmaTPCInnerC_vs_phi_pT_Cside  ->Write(); // 1/pT sigma from cov. matrix TPCInnerC
  if(h1ptSigmaTPCInnerC_vs_mult_pT_Aside) h1ptSigmaTPCInnerC_vs_mult_pT_Aside ->Write();  if(h1ptSigmaTPCInnerC_vs_mult_pT_Cside) h1ptSigmaTPCInnerC_vs_mult_pT_Cside ->Write(); // 1/pT sigma from cov. matrix vs mult. TPCInnerC

  // histogram for 1/pt shift calculation for TPCInner
  if(h1ptTPCInner_vs_eta_phi)            h1ptTPCInner_vs_eta_phi            ->Write();
  if(h1ptResTPCInner_vs_phi_pT_Aside)    h1ptResTPCInner_vs_phi_pT_Aside    ->Write();  if(h1ptResTPCInner_vs_phi_pT_Cside)    h1ptResTPCInner_vs_phi_pT_Cside    ->Write();     // 1/pT resolution from cov. matrix TPCInner
  if(h1ptResTPCInner_vs_mult_pT_Aside)   h1ptResTPCInner_vs_mult_pT_Aside   ->Write();  if(h1ptResTPCInner_vs_mult_pT_Cside)   h1ptResTPCInner_vs_mult_pT_Cside   ->Write();   // 1/pT resolution from cov. matrix vs mult. TPCInner
  if(h1ptSigmaTPCInner_vs_phi_pT_Aside)  h1ptSigmaTPCInner_vs_phi_pT_Aside  ->Write();  if(h1ptSigmaTPCInner_vs_phi_pT_Cside)  h1ptSigmaTPCInner_vs_phi_pT_Cside  ->Write();   // 1/pT sigma from cov. matrix TPCInner
  if(h1ptSigmaTPCInner_vs_mult_pT_Aside) h1ptSigmaTPCInner_vs_mult_pT_Aside ->Write();  if(h1ptSigmaTPCInner_vs_mult_pT_Cside) h1ptSigmaTPCInner_vs_mult_pT_Cside ->Write();   // 1/pT sigma from cov. matrix vs mult. TPCInner
  // Histogramm for V0s
  if(hK0sPull_vs_alpha_1pT_pos) hK0sPull_vs_alpha_1pT_pos    ->Write();
  if(hK0sRes_vs_alpha_1pT_pos)  hK0sRes_vs_alpha_1pT_pos     ->Write();
  if(hK0sPull_vs_alpha_1pT_neg) hK0sPull_vs_alpha_1pT_neg    ->Write();
  if(hK0sRes_vs_alpha_1pT_neg)  hK0sRes_vs_alpha_1pT_neg     ->Write();
  // MC info
  if(hptPull_vs_eta_pT)          hptPull_vs_eta_pT           ->Write();
  if(hptRes_vs_eta_pT)           hptRes_vs_eta_pT            ->Write();
  if(hptPullTPCInnerC_vs_eta_pT) hptPullTPCInnerC_vs_eta_pT  ->Write();
  if(hptResTPCInnerC_vs_eta_pT)  hptResTPCInnerC_vs_eta_pT   ->Write();
  if(hptPullTPCInner_vs_eta_pT)  hptPullTPCInner_vs_eta_pT   ->Write();
  if(hptResTPCInner_vs_eta_pT)   hptResTPCInner_vs_eta_pT    ->Write();

  if(hptPull_vs_phi_pT)          hptPull_vs_phi_pT           ->Write();
  if(hptRes_vs_phi_pT)           hptRes_vs_phi_pT            ->Write();
  if(hptPullTPCInnerC_vs_phi_pT) hptPullTPCInnerC_vs_phi_pT  ->Write();
  if(hptResTPCInnerC_vs_phi_pT)  hptResTPCInnerC_vs_phi_pT   ->Write();
  if(hptPullTPCInner_vs_phi_pT)  hptPullTPCInner_vs_phi_pT   ->Write();
  if(hptResTPCInner_vs_phi_pT)   hptResTPCInner_vs_phi_pT    ->Write();


  file->Close();

  TFile *OutTreeFile = new TFile("trending.root","RECREATE");
  OutTree = new TTree("trending","trending");
  OutTree->Fill();
  TBranch *brRun = OutTree->Branch("run",&runNumber);
  brRun->Fill();
  TBranch *brBz = OutTree->Branch("bz",&Bz);
  brBz->Fill();
  MakePowerFit(-1);
  MakeDCArPullFits();
  MakeDCArResFits();
  MakePhiFits();
  Make1pTresCovFits();
  MakeTPCITSMatchingEff();
  MakedcaRTrends();
  MakeDeltaPhiTrends();
  MakeK0trends();
  MakeEfficiencyTrends();

  OutTree->Write();
  OutTreeFile->Close();

  if(fMakePlots) MakeAllPlots();

}

void AliHighPtTreeAnalysis::MakeDCArPullFitsMI(){
  //
  // fit function
  //
  if (fStreamer==NULL) InitStreamer();
  
  TF1 * fFun = new TF1("fFun","[0]*TMath::Gaus(x,[1],[2])",-0.5,0.5); 
  if(!fFun) return;

  TH3D *h1PulldcaRcomb = hPulldcaR_vs_phi_pT_Aside;
  TH3D *h2PulldcaRcomb = hPulldcaR_vs_phi_pT_Cside;
  if(!h1PulldcaRcomb) return;
  if(!h2PulldcaRcomb) return;
  TH3D *h1PulldcaRTPConly = hPulldcaRTPCInner_vs_phi_pT_Aside;
  TH3D *h2PulldcaRTPConly = hPulldcaRTPCInner_vs_phi_pT_Cside;
  if(!h1PulldcaRTPConly) return;
  if(!h2PulldcaRTPConly) return;

  static Int_t    countsPulldcaR_TPCAside[2] = {0,0};
  static Int_t    countsPulldcaR_TPCCside[2] = {0,0};
  static Double_t meanPulldcaR_TPCAside[2]   = {0,0};  
  static Double_t meanPulldcaR_TPCCside[2]   = {0,0};
  static Double_t eMeanPulldcaR_TPCAside[2]  = {0,0};  
  static Double_t eMeanPulldcaR_TPCCside[2]  = {0,0};
  static Double_t rmsPulldcaR_TPCAside[2]    = {0,0};  
  static Double_t rmsPulldcaR_TPCCside[2]    = {0,0};
  static Double_t eRMSPulldcaR_TPCAside[2]   = {0,0};  
  static Double_t eRMSPulldcaR_TPCCside[2]   = {0,0};
  static Double_t shiftPulldcaR_TPCAside[2]  = {0,0};  
  static Double_t shiftPulldcaR_TPCCside[2]  = {0,0};
  static Double_t eShiftPulldcaR_TPCAside[2] = {0,0};  
  static Double_t eShiftPulldcaR_TPCCside[2] = {0,0};
  static Double_t sigmaPulldcaR_TPCAside[2]  = {0,0};  
  static Double_t sigmaPulldcaR_TPCCside[2]  = {0,0};
  static Double_t eSigmaPulldcaR_TPCAside[2] = {0,0};  
  static Double_t eSigmaPulldcaR_TPCCside[2] = {0,0};

  const char * typeName[2]={"TPConly","Combined"};

  for(Int_t itype = 0; itype <= 1; itype++){  //
    TH3D *h1PulldcaR = 0;
    TH3D *h2PulldcaR = 0;
    if(itype == 0){ 
      h1PulldcaR=h1PulldcaRTPConly;
      h2PulldcaR=h2PulldcaRTPConly;
    }
    if(itype == 1){
      h1PulldcaR=h1PulldcaRcomb;
      h2PulldcaR=h2PulldcaRcomb;
    }
    h1PulldcaR->GetXaxis()->SetRangeUser(fPtCut,100.);
    h2PulldcaR->GetXaxis()->SetRangeUser(fPtCut,100.);

    countsPulldcaR_TPCAside[itype]  = h1PulldcaR->Integral();
    countsPulldcaR_TPCCside[itype]  = h2PulldcaR->Integral();

    // calculate average values and RMS of profiles in the given pt range
    TGraphErrors *gr1PulldcaR = Calc2DProfileContent(h1PulldcaR,"zy");
    TGraphErrors *gr2PulldcaR = Calc2DProfileContent(h2PulldcaR,"zy");

    meanPulldcaR_TPCAside[itype]  = gr1PulldcaR->GetX()[0];
    eMeanPulldcaR_TPCAside[itype] = gr1PulldcaR->GetEX()[0];
    rmsPulldcaR_TPCAside[itype]   = gr1PulldcaR->GetY()[0];
    eRMSPulldcaR_TPCAside[itype]  = gr1PulldcaR->GetEY()[0];
    meanPulldcaR_TPCCside[itype]  = gr2PulldcaR->GetX()[0];
    eMeanPulldcaR_TPCCside[itype] = gr2PulldcaR->GetEX()[0];
    rmsPulldcaR_TPCCside[itype]   = gr2PulldcaR->GetY()[0];
    eRMSPulldcaR_TPCCside[itype]  = gr2PulldcaR->GetEY()[0];


    TH1D *h1PulldcaRproj = (TH1D*)h1PulldcaR->Project3D("z");
    TH1D *h2PulldcaRproj = (TH1D*)h2PulldcaR->Project3D("z");

    if(!h1PulldcaRproj) return;
    if(!h2PulldcaRproj) return;

    //  ConfigGausFit(fFun,h1PulldcaRproj);
    fFun->SetParameter(0,h1PulldcaRproj->GetMaximum());
    fFun->SetParameter(1,0);
    fFun->SetParameter(2,h1PulldcaRproj->GetRMS());
    h1PulldcaRproj->Fit(fFun,"Q");
    h1PulldcaRproj->Fit(fFun,"Q"); 
    //h1PulldcaRproj->Draw("e"); 
    shiftPulldcaR_TPCAside[itype]  = fFun->GetParameter(1);
    eShiftPulldcaR_TPCAside[itype] = fFun->GetParError(1);    
    sigmaPulldcaR_TPCAside[itype]  = TMath::Abs(fFun->GetParameter(2));
    eSigmaPulldcaR_TPCAside[itype] = TMath::Abs(fFun->GetParError(2));    

    //  ConfigGausFit(fFun,h2PulldcaRproj);
    fFun->SetParameter(0,h2PulldcaRproj->GetMaximum());
    fFun->SetParameter(1,0);
    fFun->SetParameter(2,h2PulldcaRproj->GetRMS());
    h2PulldcaRproj->Fit(fFun,"Q");
    h2PulldcaRproj->Fit(fFun,"Q");
    //h1PulldcaRproj->Draw("e"); 
    shiftPulldcaR_TPCCside[itype]  = fFun->GetParameter(1);
    eShiftPulldcaR_TPCCside[itype] = fFun->GetParError(1);    
    sigmaPulldcaR_TPCCside[itype]  = TMath::Abs(fFun->GetParameter(2));
    eSigmaPulldcaR_TPCCside[itype] = TMath::Abs(fFun->GetParError(2));    

    if(h1PulldcaRproj) delete h1PulldcaRproj;
    if(h2PulldcaRproj) delete h2PulldcaRproj;
    // set back full range
    h1PulldcaR->GetXaxis()->SetRange(1,h1PulldcaR->GetXaxis()->GetNbins());
    h2PulldcaR->GetXaxis()->SetRange(1,h2PulldcaR->GetXaxis()->GetNbins());
    (*fStreamer)<<"trending"<<
      // TString::Format("grPulldcaR_TPCAside_%s.=",typeName[itype]).Data()<< gr1PulldcaR<<  // mean pulls graphs A side
      // TString::Format("grPulldcaR_TPCCside_%s.=",typeName[itype]).Data()<< gr2PulldcaR<<  // mean pulls graphs C side 
      //
      TString::Format("countsPulldcaR_TPCAside_%s=",typeName[itype]).Data()<<countsPulldcaR_TPCAside[itype]<<  // histogram counter A side GetEntries()
      TString::Format("countsPulldcaR_TPCCside_%s=",typeName[itype]).Data()<<countsPulldcaR_TPCCside[itype]<<  // histogram counter C side GetEntries()
      TString::Format("meanPulldcaR_TPCAside_%s=",typeName[itype]).Data()<<meanPulldcaR_TPCAside[itype]<<      // mean Pull A side      
      TString::Format("meanPulldcaR_TPCCside_%s=",typeName[itype]).Data()<<meanPulldcaR_TPCCside[itype]<<      //
      TString::Format("eMeanPulldcaR_TPCAside_%s",typeName[itype]).Data()<<eMeanPulldcaR_TPCAside[itype]<<
      TString::Format("eMeanPulldcaR_TPCCside_%s",typeName[itype]).Data()<<eMeanPulldcaR_TPCCside[itype]<<
      TString::Format("rmsPulldcaR_TPCAside_%s",typeName[itype]).Data()<<rmsPulldcaR_TPCAside[itype]<<
      TString::Format("rmsPulldcaR_TPCCside_%s",typeName[itype]).Data()<<rmsPulldcaR_TPCCside[itype]<<
      TString::Format("eRMSPulldcaR_TPCAside_%s",typeName[itype]).Data()<<eRMSPulldcaR_TPCAside[itype]<<
      TString::Format("eRMSPulldcaR_TPCCside_%s",typeName[itype]).Data()<<eRMSPulldcaR_TPCCside[itype]<<
      TString::Format("shiftPulldcaR_TPCAside_%s",typeName[itype]).Data()<<shiftPulldcaR_TPCAside[itype]<<      
      TString::Format("shiftPulldcaR_TPCCside_%s",typeName[itype]).Data()<<shiftPulldcaR_TPCCside[itype]<<      
      TString::Format("eShiftPulldcaR_TPCAside_%s",typeName[itype]).Data()<<eShiftPulldcaR_TPCAside[itype]<<
      TString::Format("eShiftPulldcaR_TPCCside_%s",typeName[itype]).Data()<<eShiftPulldcaR_TPCCside[itype]<<
      TString::Format("sigmaPulldcaR_TPCAside_%s",typeName[itype]).Data()<<sigmaPulldcaR_TPCAside[itype]<<
      TString::Format("sigmaPulldcaR_TPCCside_%s",typeName[itype]).Data()<<sigmaPulldcaR_TPCCside[itype]<<
      TString::Format("eSigmaPulldcaR_TPCAside_%s",typeName[itype]).Data()<<eSigmaPulldcaR_TPCAside[itype]<<
      TString::Format("eSigmaPulldcaR_TPCCside_%s",typeName[itype]).Data()<<eSigmaPulldcaR_TPCCside[itype];
  }
}



void AliHighPtTreeAnalysis::MakeDCArPullFits(){

  // fit function
  TF1 * fFun = new TF1("fFun","[0]*TMath::Gaus(x,[1],[2])",-0.5,0.5); 
  if(!fFun) return;

  TH3D *h1PulldcaRcomb = hPulldcaR_vs_phi_pT_Aside;
  TH3D *h2PulldcaRcomb = hPulldcaR_vs_phi_pT_Cside;
  if(!h1PulldcaRcomb) return;
  if(!h2PulldcaRcomb) return;
  TH3D *h1PulldcaRTPConly = hPulldcaRTPCInner_vs_phi_pT_Aside;
  TH3D *h2PulldcaRTPConly = hPulldcaRTPCInner_vs_phi_pT_Cside;
  if(!h1PulldcaRTPConly) return;
  if(!h2PulldcaRTPConly) return;

  Int_t    countsPulldcaR_TPCAside = 0;  Int_t    countsPulldcaR_TPCCside = 0;
  Double_t meanPulldcaR_TPCAside   = 0;  Double_t meanPulldcaR_TPCCside   = 0;
  Double_t eMeanPulldcaR_TPCAside  = 0;  Double_t eMeanPulldcaR_TPCCside  = 0;
  Double_t rmsPulldcaR_TPCAside    = 0;  Double_t rmsPulldcaR_TPCCside    = 0;
  Double_t eRMSPulldcaR_TPCAside   = 0;  Double_t eRMSPulldcaR_TPCCside   = 0;
  Double_t shiftPulldcaR_TPCAside  = 0;  Double_t shiftPulldcaR_TPCCside  = 0;
  Double_t eShiftPulldcaR_TPCAside = 0;  Double_t eShiftPulldcaR_TPCCside = 0;
  Double_t sigmaPulldcaR_TPCAside  = 0;  Double_t sigmaPulldcaR_TPCCside  = 0;
  Double_t eSigmaPulldcaR_TPCAside = 0;  Double_t eSigmaPulldcaR_TPCCside = 0;

  const char * typeName[2]={"TPConly","Combined"};

  for(Int_t itype = 0; itype <= 1; itype++){
    TH3D *h1PulldcaR = 0;
    TH3D *h2PulldcaR = 0;
    if(itype == 0){ 
      h1PulldcaR=h1PulldcaRTPConly;
      h2PulldcaR=h2PulldcaRTPConly;
    }
    if(itype == 1){
      h1PulldcaR=h1PulldcaRcomb;
      h2PulldcaR=h2PulldcaRcomb;
    }
    h1PulldcaR->GetXaxis()->SetRangeUser(fPtCut,100.);
    h2PulldcaR->GetXaxis()->SetRangeUser(fPtCut,100.);

    countsPulldcaR_TPCAside  = h1PulldcaR->Integral();
    countsPulldcaR_TPCCside  = h2PulldcaR->Integral();

    // calculate average values and RMS of profiles in the given pt range
    TGraphErrors *gr1PulldcaR = Calc2DProfileContent(h1PulldcaR,"zy");
    TGraphErrors *gr2PulldcaR = Calc2DProfileContent(h2PulldcaR,"zy");

    meanPulldcaR_TPCAside  = gr1PulldcaR->GetX()[0];
    eMeanPulldcaR_TPCAside = gr1PulldcaR->GetEX()[0];
    rmsPulldcaR_TPCAside   = gr1PulldcaR->GetY()[0];
    eRMSPulldcaR_TPCAside  = gr1PulldcaR->GetEY()[0];
    meanPulldcaR_TPCCside  = gr2PulldcaR->GetX()[0];
    eMeanPulldcaR_TPCCside = gr2PulldcaR->GetEX()[0];
    rmsPulldcaR_TPCCside   = gr2PulldcaR->GetY()[0];
    eRMSPulldcaR_TPCCside  = gr2PulldcaR->GetEY()[0];

    if(!gr1PulldcaR) delete gr1PulldcaR;
    if(!gr2PulldcaR) delete gr2PulldcaR;

    TH1D *h1PulldcaRproj = (TH1D*)h1PulldcaR->Project3D("z");
    TH1D *h2PulldcaRproj = (TH1D*)h2PulldcaR->Project3D("z");

    if(!h1PulldcaRproj) return;
    if(!h2PulldcaRproj) return;

    //  ConfigGausFit(fFun,h1PulldcaRproj);
    fFun->SetParameter(0,h1PulldcaRproj->GetMaximum());
    fFun->SetParameter(1,0);
    fFun->SetParameter(2,h1PulldcaRproj->GetRMS());
    h1PulldcaRproj->Fit(fFun,"Q");
    h1PulldcaRproj->Fit(fFun,"Q"); 
    //h1PulldcaRproj->Draw("e"); 
    shiftPulldcaR_TPCAside  = fFun->GetParameter(1);
    eShiftPulldcaR_TPCAside = fFun->GetParError(1);    
    sigmaPulldcaR_TPCAside  = TMath::Abs(fFun->GetParameter(2));
    eSigmaPulldcaR_TPCAside = TMath::Abs(fFun->GetParError(2));    

    //  ConfigGausFit(fFun,h2PulldcaRproj);
    fFun->SetParameter(0,h2PulldcaRproj->GetMaximum());
    fFun->SetParameter(1,0);
    fFun->SetParameter(2,h2PulldcaRproj->GetRMS());
    h2PulldcaRproj->Fit(fFun,"Q");
    h2PulldcaRproj->Fit(fFun,"Q");
    //h1PulldcaRproj->Draw("e"); 
    shiftPulldcaR_TPCCside  = fFun->GetParameter(1);
    eShiftPulldcaR_TPCCside = fFun->GetParError(1);    
    sigmaPulldcaR_TPCCside  = TMath::Abs(fFun->GetParameter(2));
    eSigmaPulldcaR_TPCCside = TMath::Abs(fFun->GetParError(2));    

    // make trend plots
    //  if(h1PulldcaRproj->GetEntries()) Plot1DProj(h1PulldcaRproj,0,0,"DCAr/#sigma(DCAr)","dcaRPull_HighPt_TPCAside");
    //  if(h2PulldcaRproj->GetEntries()) Plot1DProj(h2PulldcaRproj,0,0,"DCAr/#sigma(DCAr)","dcaRPull_HighPt_TPCCside");
    if(h1PulldcaRproj) delete h1PulldcaRproj;
    if(h2PulldcaRproj) delete h2PulldcaRproj;
    // set back full range
    h1PulldcaR->GetXaxis()->SetRange(1,h1PulldcaR->GetXaxis()->GetNbins());
    h2PulldcaR->GetXaxis()->SetRange(1,h2PulldcaR->GetXaxis()->GetNbins());

    TBranch *br_CountsPulldcaR_Aside   = OutTree->Branch( Form("countsPulldcaR_TPCAside_%s",typeName[itype]) ,&countsPulldcaR_TPCAside);
    br_CountsPulldcaR_Aside->Fill();
    TBranch *br_CountsPulldcaR_Cside   = OutTree->Branch( Form("countsPulldcaR_TPCCside_%s",typeName[itype]) ,&countsPulldcaR_TPCCside);
    br_CountsPulldcaR_Cside->Fill();
    TBranch *br_MeanPulldcaR_Aside   = OutTree->Branch( Form("meanPulldcaR_TPCAside_%s",typeName[itype]) ,&meanPulldcaR_TPCAside);
    br_MeanPulldcaR_Aside->Fill();
    TBranch *br_MeanPulldcaR_Cside   = OutTree->Branch( Form("meanPulldcaR_TPCCside_%s",typeName[itype]) ,&meanPulldcaR_TPCCside);
    br_MeanPulldcaR_Cside->Fill();
    TBranch *br_eMeanPulldcaR_Aside   = OutTree->Branch( Form("eMeanPulldcaR_TPCAside_%s",typeName[itype]) ,&eMeanPulldcaR_TPCAside);
    br_eMeanPulldcaR_Aside->Fill();
    TBranch *br_eMeanPulldcaR_Cside   = OutTree->Branch( Form("eMeanPulldcaR_TPCCside_%s",typeName[itype]) ,&eMeanPulldcaR_TPCCside);
    br_eMeanPulldcaR_Cside->Fill();
    TBranch *br_rmsPulldcaR_Aside   = OutTree->Branch( Form("rmsPulldcaR_TPCAside_%s",typeName[itype]) ,&rmsPulldcaR_TPCAside);
    br_rmsPulldcaR_Aside->Fill();
    TBranch *br_rmsPulldcaR_Cside   = OutTree->Branch( Form("rmsPulldcaR_TPCCside_%s",typeName[itype]) ,&rmsPulldcaR_TPCCside);
    br_rmsPulldcaR_Cside->Fill();      
    TBranch *br_eRMSPulldcaR_Aside   = OutTree->Branch( Form("eRMSPulldcaR_TPCAside_%s",typeName[itype]) ,&eRMSPulldcaR_TPCAside);
    br_eRMSPulldcaR_Aside->Fill();
    TBranch *br_eRMSPulldcaR_Cside   = OutTree->Branch( Form("eRMSPulldcaR_TPCCside_%s",typeName[itype]) ,&eRMSPulldcaR_TPCCside);
    br_eRMSPulldcaR_Cside->Fill();   
    TBranch *br_shiftPulldcaR_Aside   = OutTree->Branch( Form("shiftPulldcaR_TPCAside_%s",typeName[itype]) ,&shiftPulldcaR_TPCAside);
    br_shiftPulldcaR_Aside->Fill();
    TBranch *br_shiftPulldcaR_Cside   = OutTree->Branch( Form("shiftPulldcaR_TPCCside_%s",typeName[itype]) ,&shiftPulldcaR_TPCCside);
    br_shiftPulldcaR_Cside->Fill(); 
    TBranch *br_eShiftPulldcaR_Aside   = OutTree->Branch( Form("eShiftPulldcaR_TPCAside_%s",typeName[itype]) ,&eShiftPulldcaR_TPCAside);
    br_eShiftPulldcaR_Aside->Fill();
    TBranch *br_eShiftPulldcaR_Cside   = OutTree->Branch( Form("eShiftPulldcaR_TPCCside_%s",typeName[itype]) ,&eShiftPulldcaR_TPCCside);
    br_eShiftPulldcaR_Cside->Fill(); 
    TBranch *br_sigmaPulldcaR_Aside   = OutTree->Branch( Form("sigmaPulldcaR_TPCAside_%s",typeName[itype]) ,&sigmaPulldcaR_TPCAside);
    br_sigmaPulldcaR_Aside->Fill();
    TBranch *br_sigmaPulldcaR_Cside   = OutTree->Branch( Form("sigmaPulldcaR_TPCCside_%s",typeName[itype]) ,&sigmaPulldcaR_TPCCside);
    br_sigmaPulldcaR_Cside->Fill(); 
    TBranch *br_eSigmaPulldcaR_Aside   = OutTree->Branch( Form("eSigmaPulldcaR_TPCAside_%s",typeName[itype]) ,&eSigmaPulldcaR_TPCAside);
    br_eSigmaPulldcaR_Aside->Fill();
    TBranch *br_eSigmaPulldcaR_Cside   = OutTree->Branch( Form("eSigmaPulldcaR_TPCCside_%s",typeName[itype]) ,&eSigmaPulldcaR_TPCCside);
    br_eSigmaPulldcaR_Cside->Fill();  
  }
}

void AliHighPtTreeAnalysis::MakeDCArResFits(){

  // fit function
  TF1 * fFun = new TF1("fFun","[0]*TMath::Gaus(x,[1],[2])",-0.5,0.5);
  if(!fFun) return;
  Double_t meanResdcaR_TPCAside   = 0;  Double_t meanResdcaR_TPCCside   = 0;
  Double_t eMeanResdcaR_TPCAside  = 0;  Double_t eMeanResdcaR_TPCCside  = 0;
  Double_t rmsResdcaR_TPCAside    = 0;  Double_t rmsResdcaR_TPCCside    = 0;
  Double_t eRMSResdcaR_TPCAside   = 0;  Double_t eRMSResdcaR_TPCCside   = 0;
  Double_t shiftResdcaR_TPCAside  = 0;  Double_t shiftResdcaR_TPCCside  = 0;
  Double_t eShiftResdcaR_TPCAside = 0;  Double_t eShiftResdcaR_TPCCside = 0;
  Double_t sigmaResdcaR_TPCAside  = 0;  Double_t sigmaResdcaR_TPCCside  = 0;
  Double_t eSigmaResdcaR_TPCAside = 0;  Double_t eSigmaResdcaR_TPCCside = 0;
  Int_t    countsResdcaR_TPCAside = 0;  Int_t    countsResdcaR_TPCCside = 0;
  const char * typeName[2]={"TPConly","Combined"};

  TH3D *h1ResdcaRcomb = hResdcaR_vs_phi_pT_Aside;
  TH3D *h2ResdcaRcomb = hResdcaR_vs_phi_pT_Cside;
  if(!h1ResdcaRcomb) return;
  if(!h2ResdcaRcomb) return;
  TH3D *h1ResdcaRTPConly = hResdcaRTPCInner_vs_phi_pT_Aside;
  TH3D *h2ResdcaRTPConly = hResdcaRTPCInner_vs_phi_pT_Cside;
  if(!h1ResdcaRTPConly) return;
  if(!h2ResdcaRTPConly) return;

  {for(Int_t itype = 0; itype <= 1; itype++){
                                              TH3D *h1ResdcaR = 0;
                                              TH3D *h2ResdcaR = 0;
                                              if(itype == 0){ 
                                                h1ResdcaR=h1ResdcaRTPConly;
                                                h2ResdcaR=h2ResdcaRTPConly;
                                              }
                                              if(itype == 1){
                                                h1ResdcaR=h1ResdcaRcomb;
                                                h2ResdcaR=h2ResdcaRcomb;
                                              }
                                              h1ResdcaR->GetXaxis()->SetRangeUser(fPtCut,100.);
                                              h2ResdcaR->GetXaxis()->SetRangeUser(fPtCut,100.);
                                              h1ResdcaR->GetZaxis()->SetRangeUser(-0.15,0.1499);
                                              h2ResdcaR->GetZaxis()->SetRangeUser(-0.15,0.1499);

                                              countsResdcaR_TPCAside  = h1ResdcaR->Integral();
                                              countsResdcaR_TPCCside  = h2ResdcaR->Integral();

                                              // calculate average values and RMS of profiles in the given pt range
                                              TGraphErrors *gr1ResdcaR = Calc2DProfileContent(h1ResdcaR,"zy");
                                              TGraphErrors *gr2ResdcaR = Calc2DProfileContent(h2ResdcaR,"zy");

                                              meanResdcaR_TPCAside  = gr1ResdcaR->GetX()[0];
                                              eMeanResdcaR_TPCAside = gr1ResdcaR->GetEX()[0];
                                              rmsResdcaR_TPCAside   = gr1ResdcaR->GetY()[0];
                                              eRMSResdcaR_TPCAside  = gr1ResdcaR->GetEY()[0];
                                              meanResdcaR_TPCCside  = gr2ResdcaR->GetX()[0];
                                              eMeanResdcaR_TPCCside = gr2ResdcaR->GetEX()[0];
                                              rmsResdcaR_TPCCside   = gr2ResdcaR->GetY()[0];
                                              eRMSResdcaR_TPCCside  = gr2ResdcaR->GetEY()[0];

                                              if(!gr1ResdcaR) delete gr1ResdcaR;
                                              if(!gr2ResdcaR) delete gr2ResdcaR;
                                              TH1D *h1ResdcaRproj = (TH1D*)h1ResdcaR->Project3D("z");
                                              TH1D *h2ResdcaRproj = (TH1D*)h2ResdcaR->Project3D("z");
                                              if(!h1ResdcaRproj) return;
                                              if(!h2ResdcaRproj) return;

                                              //  ConfigGausFit(fFun,h1ResdcaRproj);
                                              fFun->SetParameter(0,h1ResdcaRproj->GetMaximum());
                                              fFun->SetParameter(1,0);
                                              fFun->SetParameter(2,h1ResdcaRproj->GetRMS());

                                              h1ResdcaRproj->Fit(fFun,"Q");
                                              h1ResdcaRproj->Fit(fFun,"Q");
                                              //h1ResdcaRproj->Draw("e"); 
                                              shiftResdcaR_TPCAside  = fFun->GetParameter(1);
                                              eShiftResdcaR_TPCAside = fFun->GetParError(1);    
                                              sigmaResdcaR_TPCAside  = TMath::Abs(fFun->GetParameter(2));
                                              eSigmaResdcaR_TPCAside = TMath::Abs(fFun->GetParError(2));    
                                              //  ConfigGausFit(fFun,h2ResdcaRproj);
                                              fFun->SetParameter(0,h2ResdcaRproj->GetMaximum());
                                              fFun->SetParameter(1,0);
                                              fFun->SetParameter(2,h2ResdcaRproj->GetRMS());

                                              h2ResdcaRproj->Fit(fFun,"Q"); 
                                              h2ResdcaRproj->Fit(fFun,"Q");
                                              //h1ResdcaRproj->Draw("e"); 
                                              shiftResdcaR_TPCCside  = fFun->GetParameter(1);
                                              eShiftResdcaR_TPCCside = fFun->GetParError(1);    
                                              sigmaResdcaR_TPCCside  = TMath::Abs(fFun->GetParameter(2));
                                              eSigmaResdcaR_TPCCside = TMath::Abs(fFun->GetParError(2));    

                                              if(h1ResdcaRproj) delete h1ResdcaRproj;
                                              if(h2ResdcaRproj) delete h2ResdcaRproj;
                                              // set back full range
                                              h1ResdcaR->GetXaxis()->SetRange(1,h1ResdcaR->GetXaxis()->GetNbins());
                                              h2ResdcaR->GetXaxis()->SetRange(1,h2ResdcaR->GetXaxis()->GetNbins());
                                              h1ResdcaR->GetZaxis()->SetRange(1,h1ResdcaR->GetZaxis()->GetNbins());
                                              h2ResdcaR->GetZaxis()->SetRange(1,h2ResdcaR->GetZaxis()->GetNbins());


                                              TBranch *br_CountsResdcaR_Aside   = OutTree->Branch( Form("countsResdcaR_TPCAside_%s",typeName[itype]) ,&countsResdcaR_TPCAside);
                                              br_CountsResdcaR_Aside->Fill();
                                              TBranch *br_CountsResdcaR_Cside   = OutTree->Branch( Form("countsResdcaR_TPCCside_%s",typeName[itype]) ,&countsResdcaR_TPCCside);
                                              br_CountsResdcaR_Cside->Fill();
                                              TBranch *br_MeanResdcaR_Aside   = OutTree->Branch( Form("meanResdcaR_TPCAside_%s",typeName[itype]) ,&meanResdcaR_TPCAside);
                                              br_MeanResdcaR_Aside->Fill();
                                              TBranch *br_MeanResdcaR_Cside   = OutTree->Branch( Form("meanResdcaR_TPCCside_%s",typeName[itype]) ,&meanResdcaR_TPCCside);
                                              br_MeanResdcaR_Cside->Fill();
                                              TBranch *br_eMeanResdcaR_Aside   = OutTree->Branch( Form("eMeanResdcaR_TPCAside_%s",typeName[itype]) ,&eMeanResdcaR_TPCAside);
                                              br_eMeanResdcaR_Aside->Fill();
                                              TBranch *br_eMeanResdcaR_Cside   = OutTree->Branch( Form("eMeanResdcaR_TPCCside_%s",typeName[itype]) ,&eMeanResdcaR_TPCCside);
                                              br_eMeanResdcaR_Cside->Fill();
                                              TBranch *br_rmsResdcaR_Aside   = OutTree->Branch( Form("rmsResdcaR_TPCAside_%s",typeName[itype]) ,&rmsResdcaR_TPCAside);
                                              br_rmsResdcaR_Aside->Fill();
                                              TBranch *br_rmsResdcaR_Cside   = OutTree->Branch( Form("rmsResdcaR_TPCCside_%s",typeName[itype]) ,&rmsResdcaR_TPCCside);
                                              br_rmsResdcaR_Cside->Fill();      
                                              TBranch *br_eRMSResdcaR_Aside   = OutTree->Branch( Form("eRMSResdcaR_TPCAside_%s",typeName[itype]) ,&eRMSResdcaR_TPCAside);
                                              br_eRMSResdcaR_Aside->Fill();
                                              TBranch *br_eRMSResdcaR_Cside   = OutTree->Branch( Form("eRMSResdcaR_TPCCside_%s",typeName[itype]) ,&eRMSResdcaR_TPCCside);
                                              br_eRMSResdcaR_Cside->Fill();   
                                              TBranch *br_shiftResdcaR_Aside   = OutTree->Branch( Form("shiftResdcaR_TPCAside_%s",typeName[itype]) ,&shiftResdcaR_TPCAside);
                                              br_shiftResdcaR_Aside->Fill();
                                              TBranch *br_shiftResdcaR_Cside   = OutTree->Branch( Form("shiftResdcaR_TPCCside_%s",typeName[itype]) ,&shiftResdcaR_TPCCside);
                                              br_shiftResdcaR_Cside->Fill(); 
                                              TBranch *br_eShiftResdcaR_Aside   = OutTree->Branch( Form("eShiftResdcaR_TPCAside_%s",typeName[itype]) ,&eShiftResdcaR_TPCAside);
                                              br_eShiftResdcaR_Aside->Fill();
                                              TBranch *br_eShiftResdcaR_Cside   = OutTree->Branch( Form("eShiftResdcaR_TPCCside_%s",typeName[itype]) ,&eShiftResdcaR_TPCCside);
                                              br_eShiftResdcaR_Cside->Fill(); 
                                              TBranch *br_sigmaResdcaR_Aside   = OutTree->Branch( Form("sigmaResdcaR_TPCAside_%s",typeName[itype]) ,&sigmaResdcaR_TPCAside);
                                              br_sigmaResdcaR_Aside->Fill();
                                              TBranch *br_sigmaResdcaR_Cside   = OutTree->Branch( Form("sigmaResdcaR_TPCCside_%s",typeName[itype]) ,&sigmaResdcaR_TPCCside);
                                              br_sigmaResdcaR_Cside->Fill(); 
                                              TBranch *br_eSigmaResdcaR_Aside   = OutTree->Branch( Form("eSigmaResdcaR_TPCAside_%s",typeName[itype]) ,&eSigmaResdcaR_TPCAside);
                                              br_eSigmaResdcaR_Aside->Fill();
                                              TBranch *br_eSigmaResdcaR_Cside   = OutTree->Branch( Form("eSigmaResdcaR_TPCCside_%s",typeName[itype]) ,&eSigmaResdcaR_TPCCside);
                                              br_eSigmaResdcaR_Cside->Fill();
                                            }
  }
}

void AliHighPtTreeAnalysis::MakePhiFits(){

  Double_t meanPhiPull_TPCAside   = 0;  Double_t meanPhiPull_TPCCside   = 0;
  Double_t eMeanPhiPull_TPCAside  = 0;  Double_t eMeanPhiPull_TPCCside  = 0;
  Double_t rmsPhiPull_TPCAside    = 0;  Double_t rmsPhiPull_TPCCside    = 0;
  Double_t eRMSPhiPull_TPCAside   = 0;  Double_t eRMSPhiPull_TPCCside   = 0;
  Double_t shiftPhiPull_TPCAside  = 0;  Double_t shiftPhiPull_TPCCside  = 0;
  Double_t eShiftPhiPull_TPCAside = 0;  Double_t eShiftPhiPull_TPCCside = 0;
  Double_t sigmaPhiPull_TPCAside  = 0;  Double_t sigmaPhiPull_TPCCside  = 0;
  Double_t eSigmaPhiPull_TPCAside = 0;  Double_t eSigmaPhiPull_TPCCside = 0;
  Int_t    countsPhiPull_TPCAside = 0;  Int_t    countsPhiPull_TPCCside = 0;  
  Double_t meanPhiRes_TPCAside    = 0;  Double_t meanPhiRes_TPCCside    = 0;
  Double_t eMeanPhiRes_TPCAside   = 0;  Double_t eMeanPhiRes_TPCCside   = 0;
  Double_t rmsPhiRes_TPCAside     = 0;  Double_t rmsPhiRes_TPCCside     = 0;
  Double_t eRMSPhiRes_TPCAside    = 0;  Double_t eRMSPhiRes_TPCCside    = 0;
  Double_t shiftPhiRes_TPCAside   = 0;  Double_t shiftPhiRes_TPCCside   = 0;
  Double_t eShiftPhiRes_TPCAside  = 0;  Double_t eShiftPhiRes_TPCCside  = 0;
  Double_t sigmaPhiRes_TPCAside   = 0;  Double_t sigmaPhiRes_TPCCside   = 0;
  Double_t eSigmaPhiRes_TPCAside  = 0;  Double_t eSigmaPhiRes_TPCCside  = 0;
  Int_t    countsPhiRes_TPCAside  = 0;  Int_t    countsPhiRes_TPCCside  = 0;  
  // fit function
  TF1 * fFun = new TF1("fFun","[0]*TMath::Gaus(x,[1],[2])",-0.5,0.5);
  if(!fFun) return;

  if (!hphiPull_vs_phi_pT_Aside || !hphiPull_vs_phi_pT_Cside)
  {
    printf ("WARNING: !hphiPull_vs_phi_pT_Aside || !hphiPull_vs_phi_pT_Cside\n");
    return;
  }
  TH3D *h1phiPull = (TH3D*) hphiPull_vs_phi_pT_Aside->Clone();
  TH3D *h2phiPull = (TH3D*) hphiPull_vs_phi_pT_Cside->Clone();
  if(!h1phiPull) return;
  if(!h2phiPull) return;
  h1phiPull->GetXaxis()->SetRangeUser(fPtCut,100.);
  h2phiPull->GetXaxis()->SetRangeUser(fPtCut,100.);

  countsPhiPull_TPCAside = h1phiPull->Integral();
  countsPhiPull_TPCCside = h2phiPull->Integral();
  // calculate average values and RMS of profiles in the given pt range
  TGraphErrors *gr1phiPull = Calc2DProfileContent(h1phiPull,"zy");
  TGraphErrors *gr2phiPull = Calc2DProfileContent(h2phiPull,"zy");

  meanPhiPull_TPCAside  = gr1phiPull->GetX()[0];
  eMeanPhiPull_TPCAside = gr1phiPull->GetEX()[0];
  rmsPhiPull_TPCAside   = gr1phiPull->GetY()[0];
  eRMSPhiPull_TPCAside  = gr1phiPull->GetEY()[0];
  meanPhiPull_TPCCside  = gr2phiPull->GetX()[0];
  eMeanPhiPull_TPCCside = gr2phiPull->GetEX()[0];
  rmsPhiPull_TPCCside   = gr2phiPull->GetY()[0];
  eRMSPhiPull_TPCCside  = gr2phiPull->GetEY()[0];

  if(!gr1phiPull) delete gr1phiPull;
  if(!gr2phiPull) delete gr2phiPull;

  TH1D *h1phiPullproj = (TH1D*)h1phiPull->Project3D("z");
  TH1D *h2phiPullproj = (TH1D*)h2phiPull->Project3D("z");
  if(!h1phiPullproj) return;
  if(!h2phiPullproj) return;
  //  ConfigGausFit(fFun,h1phiPullproj);
  fFun->SetParameter(0,h1phiPullproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h1phiPullproj->GetRMS());

  h1phiPullproj->Fit(fFun,"Q");
  h1phiPullproj->Fit(fFun,"Q");
  //h1phiPullproj->Draw("e"); 
  shiftPhiPull_TPCAside  = fFun->GetParameter(1);
  eShiftPhiPull_TPCAside = fFun->GetParError(1);    
  sigmaPhiPull_TPCAside  = TMath::Abs(fFun->GetParameter(2));
  eSigmaPhiPull_TPCAside = TMath::Abs(fFun->GetParError(2));    

  //  ConfigGausFit(fFun,h2phiPullproj);
  fFun->SetParameter(0,h2phiPullproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h2phiPullproj->GetRMS());

  h2phiPullproj->Fit(fFun,"Q"); 
  h2phiPullproj->Fit(fFun,"Q");
  //h1phiPullproj->Draw("e"); 
  shiftPhiPull_TPCCside  = fFun->GetParameter(1);
  eShiftPhiPull_TPCCside = fFun->GetParError(1);    
  sigmaPhiPull_TPCCside  = TMath::Abs(fFun->GetParameter(2));
  eSigmaPhiPull_TPCCside = TMath::Abs(fFun->GetParError(2));    

  if(!h1phiPullproj) delete h1phiPullproj;
  if(!h2phiPullproj) delete h2phiPullproj;
  // set back full range
  h1phiPull->GetXaxis()->SetRange(1,h1phiPull->GetXaxis()->GetNbins());
  h2phiPull->GetXaxis()->SetRange(1,h2phiPull->GetXaxis()->GetNbins());

  // delta phi resol ///////////////////////////

  TH3D *h1phiRes = (TH3D*) hphiRes_vs_phi_pT_Aside->Clone();
  TH3D *h2phiRes = (TH3D*) hphiRes_vs_phi_pT_Cside->Clone();
  if(!h1phiRes) return;
  if(!h2phiRes) return;

  h1phiRes->GetXaxis()->SetRangeUser(fPtCut,100.);
  h2phiRes->GetXaxis()->SetRangeUser(fPtCut,100.);

  countsPhiRes_TPCAside = h1phiRes->Integral();
  countsPhiRes_TPCCside = h2phiRes->Integral();
  // calculate average values and RMS of profiles in the given pt range
  TGraphErrors *gr1phiRes = Calc2DProfileContent(h1phiRes,"zy");
  TGraphErrors *gr2phiRes = Calc2DProfileContent(h2phiRes,"zy");

  meanPhiRes_TPCAside  = gr1phiRes->GetX()[0];
  eMeanPhiRes_TPCAside = gr1phiRes->GetEX()[0];
  rmsPhiRes_TPCAside   = gr1phiRes->GetY()[0];
  eRMSPhiRes_TPCAside  = gr1phiRes->GetEY()[0];
  meanPhiRes_TPCCside  = gr2phiRes->GetX()[0];
  eMeanPhiRes_TPCCside = gr2phiRes->GetEX()[0];
  rmsPhiRes_TPCCside   = gr2phiRes->GetY()[0];
  eRMSPhiRes_TPCCside  = gr2phiRes->GetEY()[0];

  if(!gr1phiRes) delete gr1phiRes;
  if(!gr2phiRes) delete gr2phiRes;

  TH1D *h1phiResproj = (TH1D*) h1phiRes->Project3D("z");
  TH1D *h2phiResproj = (TH1D*) h2phiRes->Project3D("z");
  if(!h1phiResproj) return;
  if(!h2phiResproj) return;
  //  ConfigGausFit(fFun,h1phiResproj);
  fFun->SetParameter(0,h1phiResproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h1phiResproj->GetRMS());

  h1phiResproj->Fit(fFun,"Q");
  h1phiResproj->Fit(fFun,"Q");
  //h1phiResproj->Draw("e"); 
  shiftPhiRes_TPCAside  = fFun->GetParameter(1);
  eShiftPhiRes_TPCAside = fFun->GetParError(1);    
  sigmaPhiRes_TPCAside  = TMath::Abs(fFun->GetParameter(2));
  eSigmaPhiRes_TPCAside = TMath::Abs(fFun->GetParError(2));    
  //  ConfigGausFit(fFun,h2phiResproj);
  fFun->SetParameter(0,h2phiResproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h2phiResproj->GetRMS());

  h2phiResproj->Fit(fFun,"Q");
  h2phiResproj->Fit(fFun,"Q");
  //h1phiResproj->Draw("e"); 
  shiftPhiRes_TPCCside  = fFun->GetParameter(1);
  eShiftPhiRes_TPCCside = fFun->GetParError(1);    
  sigmaPhiRes_TPCCside  = TMath::Abs(fFun->GetParameter(2));
  eSigmaPhiRes_TPCCside = TMath::Abs(fFun->GetParError(2));    

  if(!h1phiResproj) delete h1phiResproj;
  if(!h2phiResproj) delete h2phiResproj;
  // set back full range
  h1phiRes->GetXaxis()->SetRange(1,h1phiRes->GetXaxis()->GetNbins());
  h2phiRes->GetXaxis()->SetRange(1,h2phiRes->GetXaxis()->GetNbins());

  TBranch *br_CountsPhiPull_Aside   = OutTree->Branch("countsPhiPull_TPCAside",&countsPhiPull_TPCAside);
  br_CountsPhiPull_Aside->Fill();
  TBranch *br_CountsPhiPull_Cside   = OutTree->Branch("countsPhiPull_TPCCside",&countsPhiPull_TPCCside);
  br_CountsPhiPull_Cside->Fill();
  TBranch *br_MeanPhiPull_Aside   = OutTree->Branch("meanPhiPull_TPCAside",&meanPhiPull_TPCAside);
  br_MeanPhiPull_Aside->Fill();
  TBranch *br_MeanPhiPull_Cside   = OutTree->Branch("meanPhiPull_TPCCside",&meanPhiPull_TPCCside);
  br_MeanPhiPull_Cside->Fill();
  TBranch *br_eMeanPhiPull_Aside   = OutTree->Branch("eMeanPhiPull_TPCAside",&eMeanPhiPull_TPCAside);
  br_eMeanPhiPull_Aside->Fill();
  TBranch *br_eMeanPhiPull_Cside   = OutTree->Branch("eMeanPhiPull_TPCCside",&eMeanPhiPull_TPCCside);
  br_eMeanPhiPull_Cside->Fill();
  TBranch *br_rmsPhiPull_Aside   = OutTree->Branch("rmsPhiPull_TPCAside",&rmsPhiPull_TPCAside);
  br_rmsPhiPull_Aside->Fill();
  TBranch *br_rmsPhiPull_Cside   = OutTree->Branch("rmsPhiPull_TPCCside",&rmsPhiPull_TPCCside);
  br_rmsPhiPull_Cside->Fill();      
  TBranch *br_eRMSPhiPull_Aside   = OutTree->Branch("eRMSPhiPull_TPCAside",&eRMSPhiPull_TPCAside);
  br_eRMSPhiPull_Aside->Fill();
  TBranch *br_eRMSPhiPull_Cside   = OutTree->Branch("eRMSPhiPull_TPCCside",&eRMSPhiPull_TPCCside);
  br_eRMSPhiPull_Cside->Fill();   
  TBranch *br_shiftPhiPull_Aside   = OutTree->Branch("shiftPhiPull_TPCAside",&shiftPhiPull_TPCAside);
  br_shiftPhiPull_Aside->Fill();
  TBranch *br_shiftPhiPull_Cside   = OutTree->Branch("shiftPhiPull_TPCCside",&shiftPhiPull_TPCCside);
  br_shiftPhiPull_Cside->Fill(); 
  TBranch *br_eShiftPhiPull_Aside   = OutTree->Branch("eShiftPhiPull_TPCAside",&eShiftPhiPull_TPCAside);
  br_eShiftPhiPull_Aside->Fill();
  TBranch *br_eShiftPhiPull_Cside   = OutTree->Branch("eShiftPhiPull_TPCCside",&eShiftPhiPull_TPCCside);
  br_eShiftPhiPull_Cside->Fill(); 
  TBranch *br_sigmaPhiPull_Aside   = OutTree->Branch("sigmaPhiPull_TPCAside",&sigmaPhiPull_TPCAside);
  br_sigmaPhiPull_Aside->Fill();
  TBranch *br_sigmaPhiPull_Cside   = OutTree->Branch("sigmaPhiPull_TPCCside",&sigmaPhiPull_TPCCside);
  br_sigmaPhiPull_Cside->Fill(); 
  TBranch *br_eSigmaPhiPull_Aside   = OutTree->Branch("eSigmaPhiPull_TPCAside",&eSigmaPhiPull_TPCAside);
  br_eSigmaPhiPull_Aside->Fill();
  TBranch *br_eSigmaPhiPull_Cside   = OutTree->Branch("eSigmaPhiPull_TPCCside",&eSigmaPhiPull_TPCCside);
  br_eSigmaPhiPull_Cside->Fill();   

  TBranch *br_CountsPhiRes_Aside   = OutTree->Branch("countsPhiRes_TPCAside",&countsPhiRes_TPCAside);
  br_CountsPhiRes_Aside->Fill();
  TBranch *br_CountsPhiRes_Cside   = OutTree->Branch("countsPhiRes_TPCCside",&countsPhiRes_TPCCside);
  br_CountsPhiRes_Cside->Fill();
  TBranch *br_MeanPhiRes_Aside   = OutTree->Branch("meanPhiRes_TPCAside",&meanPhiRes_TPCAside);
  br_MeanPhiRes_Aside->Fill();
  TBranch *br_MeanPhiRes_Cside   = OutTree->Branch("meanPhiRes_TPCCside",&meanPhiRes_TPCCside);
  br_MeanPhiRes_Cside->Fill();
  TBranch *br_eMeanPhiRes_Aside   = OutTree->Branch("eMeanPhiRes_TPCAside",&eMeanPhiRes_TPCAside);
  br_eMeanPhiRes_Aside->Fill();
  TBranch *br_eMeanPhiRes_Cside   = OutTree->Branch("eMeanPhiRes_TPCCside",&eMeanPhiRes_TPCCside);
  br_eMeanPhiRes_Cside->Fill();
  TBranch *br_rmsPhiRes_Aside   = OutTree->Branch("rmsPhiRes_TPCAside",&rmsPhiRes_TPCAside);
  br_rmsPhiRes_Aside->Fill();
  TBranch *br_rmsPhiRes_Cside   = OutTree->Branch("rmsPhiRes_TPCCside",&rmsPhiRes_TPCCside);
  br_rmsPhiRes_Cside->Fill();      
  TBranch *br_eRMSPhiRes_Aside   = OutTree->Branch("eRMSPhiRes_TPCAside",&eRMSPhiRes_TPCAside);
  br_eRMSPhiRes_Aside->Fill();
  TBranch *br_eRMSPhiRes_Cside   = OutTree->Branch("eRMSPhiRes_TPCCside",&eRMSPhiRes_TPCCside);
  br_eRMSPhiRes_Cside->Fill();   
  TBranch *br_shiftPhiRes_Aside   = OutTree->Branch("shiftPhiRes_TPCAside",&shiftPhiRes_TPCAside);
  br_shiftPhiRes_Aside->Fill();
  TBranch *br_shiftPhiRes_Cside   = OutTree->Branch("shiftPhiRes_TPCCside",&shiftPhiRes_TPCCside);
  br_shiftPhiRes_Cside->Fill(); 
  TBranch *br_eShiftPhiRes_Aside   = OutTree->Branch("eShiftPhiRes_TPCAside",&eShiftPhiRes_TPCAside);
  br_eShiftPhiRes_Aside->Fill();
  TBranch *br_eShiftPhiRes_Cside   = OutTree->Branch("eShiftPhiRes_TPCCside",&eShiftPhiRes_TPCCside);
  br_eShiftPhiRes_Cside->Fill(); 
  TBranch *br_sigmaPhiRes_Aside   = OutTree->Branch("sigmaPhiRes_TPCAside",&sigmaPhiRes_TPCAside);
  br_sigmaPhiRes_Aside->Fill();
  TBranch *br_sigmaPhiRes_Cside   = OutTree->Branch("sigmaPhiRes_TPCCside",&sigmaPhiRes_TPCCside);
  br_sigmaPhiRes_Cside->Fill(); 
  TBranch *br_eSigmaPhiRes_Aside   = OutTree->Branch("eSigmaPhiRes_TPCAside",&eSigmaPhiRes_TPCAside);
  br_eSigmaPhiRes_Aside->Fill();
  TBranch *br_eSigmaPhiRes_Cside   = OutTree->Branch("eSigmaPhiRes_TPCCside",&eSigmaPhiRes_TPCCside);
  br_eSigmaPhiRes_Cside->Fill();
}

void AliHighPtTreeAnalysis::Make1pTresCovFits(){

  Double_t cov1ptRes_TPCAside    = 0;   Double_t cov1ptRes_TPCCside    = 0;
  Double_t eCov1ptRes_TPCAside   = 0;   Double_t eCov1ptRes_TPCCside   = 0;
  Int_t    counts1ptRes_TPCAside = 0;   Int_t    counts1ptRes_TPCCside = 0;

  // fit function
  TF1 * fFun = new TF1("fFun","[0]*TMath::Gaus(x,[1],[2])",-0.5,0.5);
  if(!fFun) return;

  if (!h1ptRes_vs_phi_pT_Aside) return;
  if (!h1ptRes_vs_phi_pT_Cside) return;
  TH3D *h11ptRes = (TH3D*) h1ptRes_vs_phi_pT_Aside->Clone();
  TH3D *h21ptRes = (TH3D*) h1ptRes_vs_phi_pT_Cside->Clone();
  if(!h11ptRes) return;
  if(!h21ptRes) return;

  h11ptRes->GetXaxis()->SetRangeUser(fPtCut,100.);
  h21ptRes->GetXaxis()->SetRangeUser(fPtCut,100.);

  counts1ptRes_TPCAside = h11ptRes->Integral();
  counts1ptRes_TPCCside = h21ptRes->Integral();

  //  TH1D *h11ptResproj = Project1D(h11ptRes,"z");
  TH1D *h11ptResproj = (TH1D*) h11ptRes->Project3D("z");
  //  TH1D *h21ptResproj = Project1D(h21ptRes,"z");
  TH1D *h21ptResproj = (TH1D*) h21ptRes->Project3D("z");
  if(!h11ptResproj) return;
  if(!h21ptResproj) return;

  //  ConfigGausFit(fFun,h11ptResproj);
  fFun->SetParameter(0,h11ptResproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h11ptResproj->GetRMS());
  fFun->SetParameter(1,h11ptResproj->GetMean());
  h11ptResproj->Fit(fFun,"Q");
  h11ptResproj->Fit(fFun,"Q");
  //h11ptResproj->Draw("e"); 
  cov1ptRes_TPCAside  = fFun->GetParameter(1);
  eCov1ptRes_TPCAside = fFun->GetParError(1);
  //  ConfigGausFit(fFun,h21ptResproj);
  fFun->SetParameter(0,h21ptResproj->GetMaximum());
  fFun->SetParameter(1,0);
  fFun->SetParameter(2,h21ptResproj->GetRMS());

  fFun->SetParameter(1,h21ptResproj->GetMean());
  h21ptResproj->Fit(fFun,"Q");
  h21ptResproj->Fit(fFun,"Q");
  //h11ptResproj->Draw("e"); 
  cov1ptRes_TPCCside  = fFun->GetParameter(1);
  eCov1ptRes_TPCCside = fFun->GetParError(1);

  if(!h11ptResproj) delete h11ptResproj;
  if(!h21ptResproj) delete h21ptResproj;

  // set back full range
  h11ptRes->GetXaxis()->SetRange(1,h11ptRes->GetXaxis()->GetNbins());
  h21ptRes->GetXaxis()->SetRange(1,h21ptRes->GetXaxis()->GetNbins());

  TBranch *br_Counts1ptRes_TPCAside   = OutTree->Branch("counts1ptRes_TPCAside",&counts1ptRes_TPCAside);
  br_Counts1ptRes_TPCAside->Fill();
  TBranch *br_Counts1ptRes_TPCCside   = OutTree->Branch("counts1ptRes_TPCCside",&counts1ptRes_TPCCside);
  br_Counts1ptRes_TPCCside->Fill();  
  TBranch *br_cov1ptRes_TPCAside   = OutTree->Branch("cov1ptRes_TPCAside",&cov1ptRes_TPCAside);
  br_cov1ptRes_TPCAside->Fill();
  TBranch *br_cov1ptRes_TPCCside   = OutTree->Branch("cov1ptRes_TPCCside",&cov1ptRes_TPCCside);
  br_cov1ptRes_TPCCside->Fill();  
  TBranch *br_eCov1ptRes_TPCAside   = OutTree->Branch("eCov1ptRes_TPCAside",&eCov1ptRes_TPCAside);
  br_eCov1ptRes_TPCAside->Fill();
  TBranch *br_eCov1ptRes_TPCCside   = OutTree->Branch("eCov1ptRes_TPCCside",&eCov1ptRes_TPCCside);
  br_eCov1ptRes_TPCCside->Fill(); 
}

TGraphErrors* AliHighPtTreeAnalysis::Calc2DProfileContent(TH3D *h1, const char *projAxisName){
  // Make projections
  if(!h1) return 0;
  char name[256];

  sprintf(name,"%s_1",h1->GetName());
  TH3D *h1c = (TH3D*)h1->Clone(name);
  if(!h1c) return 0;
  if(fPtCut>0.)  h1c->GetXaxis()->SetRangeUser(fPtCut,100);

  TH2D *h1proj = (TH2D*)h1c->Project3D(projAxisName);
  if(!h1proj) return 0;

  // Fit slices
  TObjArray *arr1 = new TObjArray();
  h1proj->FitSlicesY(0,0,-1,0,"QNR",arr1);

  if(!arr1->At(1)) return 0;
  if(!arr1->At(2)) return 0;

  // Get histo
  TH1D *h1mean  = (TH1D*)arr1->At(1);
  TH1D *h1width = (TH1D*)arr1->At(2);

  TH1D *htempMean = new TH1D("htempMean","htempMean",100,h1mean->GetMinimum(),h1mean->GetMaximum());
  if(!htempMean) return 0;
  TH1D *htempWidth = new TH1D("htempWidth","htempWidth",100,h1width->GetMinimum(),h1width->GetMaximum());
  if(!htempWidth) return 0;

  for (Int_t i=1; i<h1mean->GetNbinsX(); i++) { htempMean->Fill(h1mean->GetBinContent(i)); }
  for (Int_t i=1; i<h1width->GetNbinsX(); i++) { htempWidth->Fill(h1width->GetBinContent(i)); }

  //
  Double_t x[1] = {0};
  Double_t y[1] = {0};
  Double_t ex[1]= {0};
  Double_t ey[1]= {0};

  x[0] = htempMean->GetMean();
  y[0] = htempWidth->GetMean();
  ex[0] = htempMean->GetRMS();
  ey[0] = htempWidth->GetRMS();

  sprintf(name,"%s_gr_1",h1->GetName());
  TGraphErrors *gr = new TGraphErrors(1,x,y,ex,ey);
  gr->SetName(name);


  delete h1proj;
  delete arr1;
  delete htempMean;
  delete htempWidth;

  return gr;
}

void AliHighPtTreeAnalysis::RunPeriod()
{
  //if(fBfield ==  1)
  //{ 
  //  gSystem->mkdir("periodBpos");
  //  gSystem->cd("periodBpos"); 
  //}
  //if(fBfield == -1)
  //{ 
  //  gSystem->mkdir("periodBneg");
  //  gSystem->cd("periodBneg"); 
  //}

  TFile *OutTreeFile = new TFile("PeriodOutput.root","RECREATE");
  OutTree = new TTree("PeriodTree","PeriodTree");
  OutTree->Fill();
  TBranch *brPeriod = OutTree->Branch("period",&fPeriodName);
  brPeriod->Fill();

  MakePowerFit(-1);
  MakeDCArPullFits();
  MakeDCArResFits();
  MakePhiFits();
  Make1pTresCovFits();
  MakeTPCITSMatchingEff();
  MakedcaRTrends();
  MakeDeltaPhiTrends();
  MakeK0trends();
  MakeEfficiencyTrends();

  OutTree->Write();
  OutTreeFile->Close();

  if(fMakePlots) MakeAllPlots();

}

Double_t AliHighPtTreeAnalysis::qoverptCorr(Double_t trEta, Double_t trPhi, Int_t type){
  // Type 0:  Combined Tracks
  // Type 1:  TPConly Tracks
  // Type 2:  TPCconstrained Tracks

  Double_t tr1ptCorr = 0;
  if(type == 0){
    if(trEta > 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionAside[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionAside[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionAside[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionAside[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionAside[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionAside[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionAside[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionAside[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionAside[17];
    }
    if(trEta < 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionCside[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionCside[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionCside[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionCside[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionCside[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionCside[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionCside[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionCside[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionCside[17];
    }
  }
  if(type == 1){
    if(trEta > 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInner[17];
    }
    if(trEta < 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInner[17];
    }
  }
  if(type == 2){
    if(trEta > 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionAsideTPCInnerC[17];
    }
    if(trEta < 0){
      if( trPhi > 0  * 6.28/18 && trPhi < 1  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 0];
      if( trPhi > 1  * 6.28/18 && trPhi < 2  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 1];
      if( trPhi > 2  * 6.28/18 && trPhi < 3  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 2];
      if( trPhi > 3  * 6.28/18 && trPhi < 4  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 3];
      if( trPhi > 4  * 6.28/18 && trPhi < 5  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 4];
      if( trPhi > 5  * 6.28/18 && trPhi < 6  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 5];
      if( trPhi > 6  * 6.28/18 && trPhi < 7  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 6];
      if( trPhi > 7  * 6.28/18 && trPhi < 8  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 7];
      if( trPhi > 8  * 6.28/18 && trPhi < 9  * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 8];
      if( trPhi > 9  * 6.28/18 && trPhi < 10 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[ 9];
      if( trPhi > 10 * 6.28/18 && trPhi < 11 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[10];
      if( trPhi > 11 * 6.28/18 && trPhi < 12 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[11];
      if( trPhi > 12 * 6.28/18 && trPhi < 13 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[12];
      if( trPhi > 13 * 6.28/18 && trPhi < 14 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[13];
      if( trPhi > 14 * 6.28/18 && trPhi < 15 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[14];
      if( trPhi > 15 * 6.28/18 && trPhi < 16 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[15];
      if( trPhi > 16 * 6.28/18 && trPhi < 17 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[16];
      if( trPhi > 17 * 6.28/18 && trPhi < 18 * 6.28/18 )  tr1ptCorr = fCorrectionCsideTPCInnerC[17];
    }
  }
  return tr1ptCorr;
}

void AliHighPtTreeAnalysis::MakePowerFit(Int_t entries){
  //
  // Make a 1/pt histograms : TPC only, TPConly+constraint, Combined 
  // Fit the PowerLaw model + 1/pt shift 
  // common, A side, C side, per sector

  //  Trending Variables
  Double_t  qptShiftCombinedSinAside = -999.;   Double_t  qptShiftCombinedSinCside = -999.;
  Double_t  qptShiftCombinedCosAside = -999.;   Double_t  qptShiftCombinedCosCside = -999.;
  Double_t eqptShiftCombinedSinAside = -999.;   Double_t eqptShiftCombinedSinCside = -999.;
  Double_t eqptShiftCombinedCosAside = -999.;   Double_t eqptShiftCombinedCosCside = -999.;

  Double_t  qptShiftTPCconstSinAside = -999.;   Double_t  qptShiftTPCconstSinCside = -999.;
  Double_t  qptShiftTPCconstCosAside = -999.;   Double_t  qptShiftTPCconstCosCside = -999.;
  Double_t eqptShiftTPCconstSinAside = -999.;   Double_t eqptShiftTPCconstSinCside = -999.;
  Double_t eqptShiftTPCconstCosAside = -999.;   Double_t eqptShiftTPCconstCosCside = -999.;

  Double_t  qptShiftTPConlySinAside  = -999.;   Double_t  qptShiftTPConlySinCside  = -999.;
  Double_t  qptShiftTPConlyCosAside  = -999.;   Double_t  qptShiftTPConlyCosCside  = -999.;
  Double_t eqptShiftTPConlySinAside  = -999.;   Double_t eqptShiftTPConlySinCside  = -999.;
  Double_t eqptShiftTPConlyCosAside  = -999.;   Double_t eqptShiftTPConlyCosCside  = -999.;

  if (entries<0) entries=1000000000;
  TH3F * phis1ptThetaAlphaComb      =  h1pt_vs_eta_phi;
  TH3F * phis1ptThetaAlphaTPCC      =  h1ptTPCInnerC_vs_eta_phi;
  TH3F * phis1ptThetaAlphaTPC       =  h1ptTPCInner_vs_eta_phi;
  //
  //return;
  /*
     TFile f("highPtPic.root");
     TH3F *phis1ptThetaAlphaTPC=f.Get("his1ptThetaAlphaTPC");
     TH3F *phis1ptThetaAlphaTPCC=f.Get("his1ptThetaAlphaTPCC");
     TH3F *phis1ptThetaAlphaComb=f.Get("his1ptThetaAlphaComb");
     */
  //
  // fit functions
  //

  //  TF1 *fpower = new TF1("fpower","[0]*(abs(1/(x-[2])))/pow(abs(1/(x-[2])),[1])*abs( 1/(x+0.0001)- 1/(x-0.0001))",-1,1);
  //  TF1 *fpower = new TF1("fpower","[0]*(abs(1/(x-[2])))/pow(abs(1/(x-[2])),[1])",-1,1);
  //  TF1 *fpower = new TF1("fpower","[0]*(abs(1/(x-[2])))/pow(abs(1/(x-[2])),max(min([1],8),1))",-1,1);
  TF1 *fpower = new TF1("fpower","[0]*pow( abs(x-[2]) , abs([1]))",-1,1);
  TF1 *fphi  =  new TF1("fphi","[0]+[1]*cos(x)+[2]*sin(x)");
  fpower->SetNpx(300);
  fphi->SetParName(0,"Offset");
  fphi->SetParName(1,"cos(#phi)");
  fphi->SetParName(2,"sin(#phi)");
  fpower->SetParName(0,"Norm.");
  fpower->SetParName(1,"Slope");
  fpower->SetParName(2,"#Delta_{1/pt}");

  TCanvas *workCanvas = new TCanvas;
  TVectorD vecFitAll(4),  vecFitAside(4), vecFitCside(4);            // fit of the 1/pt shift per side
  TVectorD eFitAll(4),  eFitAside(4), eFitCside(4); 
  TVectorD vecFitSecA(18), vecFitSecC(18),sectors(18),phiPos(18);    // fit per sector  
  TVectorD vecFitSlopeSecA(18), vecFitSlopeSecC(18);    // fit per sector  
  TVectorD eFitSecA(18), eFitSecC(18), eFitSlopeSecA(18), eFitSlopeSecC(18); 
  TVectorD vecFitA(3), vecFitC(3),  eFitA(3), eFitC(3);              // global fit    
  //TGraphErrors *grSlopeA[3]={0x0,0x0,0x0};
  //TGraphErrors *grSlopeC[3]={0x0,0x0,0x0};

  TVectorD vecChi2SecA(18), vecChi2SecC(18);
  TVectorD vecNpointsSecA(18), vecNpointsSecC(18);

  Double_t chi2All(0), chi2Aside(0), chi2Cside(0);
  Double_t nPointsAll(0), nPointsAside(0), nPointsCside(0);

  const char * typeTitles[3]={"TPC only", "TPC constrained", "Combined"};
  const char * typeName[3]={"TPConly", "TPCconstrained","Combined"};

  for (Int_t itype=0; itype<3; itype++){
    workCanvas->cd();
    TH3 * phis1ptThetaAlpha=0;
    if (itype==2) phis1ptThetaAlpha=phis1ptThetaAlphaComb;
    if (itype==0) phis1ptThetaAlpha=phis1ptThetaAlphaTPC;
    if (itype==1) phis1ptThetaAlpha=phis1ptThetaAlphaTPCC;
    if (!phis1ptThetaAlpha) 
    {
      printf("WARNING: phis1ptThetaAlpha in NULL, continuing\n");
      continue;
    }
    //
    gStyle->SetOptTitle(1);
    gStyle->SetOptFit(1);


    //    TH1 * hisAll       = phis1ptThetaAlpha->ProjectionX("all",1,7);
    //    TH1 * hisAside     = phis1ptThetaAlpha->ProjectionX("aside",1,3);
    //    TH1 * hisCside     = phis1ptThetaAlpha->ProjectionX("cside",4,7);

    TH1 * hisAll       = phis1ptThetaAlpha->ProjectionX("all",1,10);
    TH1 * hisAside     = phis1ptThetaAlpha->ProjectionX("aside",1,5);
    TH1 * hisCside     = phis1ptThetaAlpha->ProjectionX("cside",6,10);

    //    std::cout << "hisAll   Entries: " << hisAll->Integral(hisAll->FindBin(-0.12),hisAll->FindBin(0.12))       << std::endl;
    //    std::cout << "hisAside Entries: " << hisAside->Integral(hisAside->FindBin(-0.12),hisAside->FindBin(0.12)) << std::endl;
    //    std::cout << "hisCside Entries: " << hisCside->Integral(hisCside->FindBin(-0.12),hisCside->FindBin(0.12)) << std::endl;
    if(hisAside->Integral(hisAside->FindBin(-0.12),hisAside->FindBin(0.12)) >= 100 && hisCside->Integral(hisCside->FindBin(-0.12),hisCside->FindBin(0.12)) >= 100){


      hisAll->Sumw2();
      hisAside->Sumw2();
      hisCside->Sumw2();
      hisAll->Scale(1/hisAll->GetEntries());
      hisAside->Scale(1/hisAside->GetEntries());
      hisCside->Scale(1/hisCside->GetEntries());

      hisAll->SetMarkerStyle(25); 
      hisAside->SetMarkerStyle(25); hisAside->SetMarkerColor(2); 
      hisCside->SetMarkerStyle(25);  hisCside->SetMarkerColor(4);
      hisAll->SetName(Form("%s: Both sides",typeTitles[itype]));   hisAll->SetTitle(Form("%s Both side",typeTitles[itype]));  
      hisAll->GetXaxis()->SetTitle("1/p_{t} (GeV/c)");
      hisAside->SetName(Form("%s:A side",typeTitles[itype]));  hisAside->SetTitle(Form("%s A side",typeTitles[itype]));
      hisAside->GetXaxis()->SetTitle("1/p_{t} (GeV/c)");
      hisCside->SetName(Form("%s: C side",typeTitles[itype]));  hisCside->SetTitle(Form("%s C side",typeTitles[itype]));
      hisCside->GetXaxis()->SetTitle("1/p_{t} (GeV/c)");
      //
      fpower->SetParameters(1,4,0);
      Double_t norm = 1;
      if (fpower->Eval(hisAll->GetBinCenter(hisAll->FindBin(-0.12))) > 0) { 
        norm = hisAll->GetBinContent(hisAll->FindBin(-0.12))/ fpower->Eval(hisAll->GetBinCenter(hisAll->FindBin(-0.12))); 
        fpower->SetParameter(0,norm);
      }
      fpower->FixParameter(2,0);
      hisAll->Fit(fpower,"Q","err",-0.12,0.12);
      if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
      fpower->ReleaseParameter(2);
      fpower->SetLineColor(1);
      hisAll->Fit(fpower,"Q","err",-0.12,0.12);
      if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
      fpower->GetParameters(vecFitAll.GetMatrixArray());
      chi2All = fpower->GetChisquare()/fpower->GetNDF();
      nPointsAll = hisAll->Integral( hisAll->FindBin(-0.12),hisAll->FindBin(0.12) );

      Double_t slopeAll = fpower->GetParameter(1);
      Double_t shiftAll = fpower->GetParameter(2);
      for (Int_t i=0; i<3;i++) eFitAll[i]=fpower->GetParError(i);
      //Aside
      fpower->SetParameters(1,slopeAll,shiftAll);
      if (fpower->Eval(hisAside->GetBinCenter(hisAside->FindBin(-0.12))) > 0) { 
        norm = hisAside->GetBinContent(hisAside->FindBin(-0.12))/ fpower->Eval(hisAside->GetBinCenter(hisAside->FindBin(-0.12)));
        fpower->SetParameter(0,norm);
      }
      fpower->SetLineColor(2);
      hisAside->Fit(fpower,"Q","err",-0.12,0.12);
      if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
      fpower->GetParameters(vecFitAside.GetMatrixArray()); 
      chi2Aside = fpower->GetChisquare()/fpower->GetNDF();
      nPointsAside = hisAside->Integral( hisAside->FindBin(-0.12),hisAside->FindBin(0.12) );
      for (Int_t i=0; i<3;i++) eFitAside[i]=fpower->GetParError(i);
      //Cside
      fpower->SetParameters(1,slopeAll,shiftAll);
      if (fpower->Eval(hisCside->GetBinCenter(hisCside->FindBin(-0.12))) > 0) { 
        norm = hisCside->GetBinContent(hisCside->FindBin(-0.12))/ fpower->Eval(hisCside->GetBinCenter(hisCside->FindBin(-0.12))); 
        fpower->SetParameter(0,norm);
      }
      fpower->SetLineColor(4);
      hisCside->Fit(fpower,"Q","err",-0.12,0.12);
      if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
      fpower->GetParameters(vecFitCside.GetMatrixArray());
      chi2Cside = fpower->GetChisquare()/fpower->GetNDF();
      nPointsCside = hisCside->Integral( hisCside->FindBin(-0.12),hisCside->FindBin(0.12) );
      for (Int_t i=0; i<3;i++) eFitCside[i]=fpower->GetParError(i);

      gStyle->SetOptTitle(1);
      gStyle->SetOptStat(1);
      workCanvas->cd(); workCanvas->Update();
      hisAll->Draw("err"); 
      workCanvas->Update();
      TPaveStats *statAll = (TPaveStats*)hisAll->FindObject("stats");
      hisAside->Draw("err"); 
      workCanvas->Update();
      TPaveStats *statA = (TPaveStats*)hisAside->FindObject("stats");
      hisCside->Draw("err");
      workCanvas->Update(); 
      TPaveStats *statC = (TPaveStats*)hisCside->FindObject("stats");
      statAll->SetBorderSize(4); statA->SetBorderSize(4); statC->SetBorderSize(4);
      statAll->SetTextColor(1); statA->SetTextColor(2); statC->SetTextColor(4);
      statAll->SetOptStat(1000000111); 
      statA->  SetOptStat(1000000111); 
      statC->  SetOptStat(1000000111);
      statAll->SetX1NDC(0.75); statAll->SetX2NDC(0.99);
      statA->SetX1NDC(0.75); statA->SetX2NDC(0.99);
      statC->SetX1NDC(0.75); statC->SetX2NDC(0.99);
      statA->SetY1NDC(0.7); statA->SetY2NDC(0.95);
      statC->SetY1NDC(0.40); statC->SetY2NDC(0.65);
      statAll->SetY1NDC(0.1); statAll->SetY2NDC(0.35);
      gPad->SetLogy();
      gPad->SetTicky(3);gPad->SetTickx(3);
      gPad->SetLeftMargin(0.05);gPad->SetRightMargin(0.25);
      gStyle->SetOptTitle(0); 
      hisAll->Draw("err func same");
      hisCside->Draw("err func same");
      hisAside->Draw("err func same");

      //
      for (Int_t isec=0; isec<18; isec++){
        Int_t isec0=isec, isec1=isec+2;
        if (isec0<1) isec0=1;
        if (isec1>18) isec1=18;
        TH1 *hisAphi = phis1ptThetaAlpha->ProjectionX(Form("Aside%d",isec),1, 5,isec0,isec1);
        TH1 *hisCphi = phis1ptThetaAlpha->ProjectionX(Form("Cside%d",isec),6,10,isec0,isec1);
        fpower->SetParameters(vecFitAll.GetMatrixArray());

        Bool_t fitA = ( hisAphi->Integral(hisAphi->FindBin(-0.12),hisAphi->FindBin(0.12)) > 50 );
        Bool_t fitC = ( hisCphi->Integral(hisCphi->FindBin(-0.12),hisCphi->FindBin(0.12)) > 50 );

        //A side
        if( fitA ){
          fpower->SetParameter(0, vecFitAll[0]/6);
          fpower->FixParameter(1, vecFitAside[1]);
          hisAphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          vecFitSecA[isec]= fpower->GetParameter(2);
          eFitSecA[isec]  = fpower->GetParError(2)*TMath::Sqrt(fpower->GetChisquare()/fpower->GetNDF());
          fpower->FixParameter(2, vecFitAside[2]);
          hisAphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          fpower->ReleaseParameter(1);
          hisAphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          vecFitSlopeSecA[isec]= fpower->GetParameter(1);
          eFitSlopeSecA[isec]  = fpower->GetParError(1)*TMath::Sqrt(fpower->GetChisquare()/fpower->GetNDF());
          fpower->ReleaseParameter(2);
          vecChi2SecA[isec]    = fpower->GetChisquare()/fpower->GetNDF();
          vecNpointsSecA[isec] = hisAphi->Integral( hisAphi->FindBin(-0.12),hisAphi->FindBin(0.12) );
        }
        //        else std::cout << "Skipping Aside of Sector: " << isec << std::endl;
        //C side
        if( fitC ){
          fpower->FixParameter(1, vecFitCside[1]);
          hisCphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          vecFitSecC[isec]= fpower->GetParameter(2);
          eFitSecC[isec]  = fpower->GetParError(2)*TMath::Sqrt(fpower->GetChisquare()/fpower->GetNDF());
          fpower->FixParameter(2, vecFitCside[2]);
          hisCphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          fpower->ReleaseParameter(1);
          hisCphi->Fit(fpower,"Q","err",-0.12,0.12);
          if(gMinuit->fStatus != 0) std::cout << "Minuit Status:   " << gMinuit->fStatus <<std::endl;
          vecFitSlopeSecC[isec]= fpower->GetParameter(1);
          eFitSlopeSecC[isec]  = fpower->GetParError(1)*TMath::Sqrt(fpower->GetChisquare()/fpower->GetNDF());
          fpower->ReleaseParameter(2);
          vecChi2SecC[isec]    = fpower->GetChisquare()/fpower->GetNDF();
          vecNpointsSecC[isec] = hisCphi->Integral( hisCphi->FindBin(-0.12),hisCphi->FindBin(0.12) );
        }
        //       else std::cout << "Skipping Cside of Sector: " << isec << std::endl;

        sectors[isec]=isec;
        phiPos[isec]=phis1ptThetaAlpha->GetZaxis()->GetBinCenter(isec+1); 
        //	
      }

      // use common errors
      Double_t median = 0;
      Double_t medianP = 0;
      median=TMath::Median(18, eFitSecA.GetMatrixArray());
      medianP=TMath::Median(18, vecFitSecA.GetMatrixArray());
      for (Int_t isec=0; isec<18; isec++) eFitSecA[isec]=median;	
      for (Int_t isec=0; isec<18; isec++) if (TMath::Abs(vecFitSecA[isec])>0.004+median) vecFitSecA[isec]=medianP;	
      //
      median=TMath::Median(18, eFitSecC.GetMatrixArray());
      medianP=TMath::Median(18, vecFitSecC.GetMatrixArray());
      for (Int_t isec=0; isec<18; isec++) eFitSecC[isec]=median;	
      for (Int_t isec=0; isec<18; isec++) if (TMath::Abs(vecFitSecC[isec])>0.004+median) vecFitSecC[isec]=medianP;	

      //
      median=TMath::Median(18, eFitSlopeSecA.GetMatrixArray());
      for (Int_t isec=0; isec<18; isec++) eFitSlopeSecA[isec]=median;	
      median=TMath::Median(18, eFitSlopeSecC.GetMatrixArray());
      for (Int_t isec=0; isec<18; isec++) eFitSlopeSecC[isec]=median;	


      gStyle->SetOptTitle(1);

      TGraph* grChi2Aside    = new TGraph(18, phiPos.GetMatrixArray(),vecChi2SecA.GetMatrixArray());
      TGraph* grChi2Cside    = new TGraph(18, phiPos.GetMatrixArray(),vecChi2SecC.GetMatrixArray());
      TGraph* grNpointsAside = new TGraph(18, phiPos.GetMatrixArray(),vecNpointsSecA.GetMatrixArray());
      TGraph* grNpointsCside = new TGraph(18, phiPos.GetMatrixArray(),vecNpointsSecC.GetMatrixArray());

      TGraphErrors* grAside= new TGraphErrors(18, phiPos.GetMatrixArray(),vecFitSecA.GetMatrixArray(), 0, eFitSecA.GetMatrixArray());
      TGraphErrors* grCside= new TGraphErrors(18, phiPos.GetMatrixArray(),vecFitSecC.GetMatrixArray(), 0, eFitSecC.GetMatrixArray());
      TGraphErrors* grAsideSlope= new TGraphErrors(18, phiPos.GetMatrixArray(),vecFitSlopeSecA.GetMatrixArray(), 0, eFitSlopeSecA.GetMatrixArray());
      TGraphErrors* grCsideSlope= new TGraphErrors(18, phiPos.GetMatrixArray(),vecFitSlopeSecC.GetMatrixArray(), 0, eFitSlopeSecC.GetMatrixArray());
      grAside->SetMarkerStyle(25);  grCside->SetMarkerStyle(25);
      grAside->SetMarkerColor(2);   grCside->SetMarkerColor(4);
      grAsideSlope->SetMarkerStyle(25);  grCsideSlope->SetMarkerStyle(25);
      grAsideSlope->SetMarkerColor(2);   grCsideSlope->SetMarkerColor(4);
      fphi->SetLineColor(2); grAside->Fit(fphi,"Q");
      fphi->GetParameters(vecFitA.GetMatrixArray());
      for (Int_t i=0; i<3;i++) eFitA[i]=fphi->GetParError(i);
      fphi->SetLineColor(4); grCside->Fit(fphi,"Q");
      fphi->GetParameters(vecFitC.GetMatrixArray());
      grAside->SetMaximum(0.006);
      grAside->SetMinimum(-0.006);
      grAside->SetName(Form("%s_ASideShift",typeTitles[itype]));
      grCside->SetName(Form("%s_CSideShift",typeTitles[itype]));
      grAside->SetTitle(Form("%s: A Side Shift",typeTitles[itype]));
      grCside->SetTitle(Form("%s: C side Shift",typeTitles[itype]));
      grAside->GetXaxis()->SetTitle("#phi");
      grAside->GetYaxis()->SetTitle("#Delta_{1/pt} (1/GeV)");
      grAsideSlope->SetName(Form("%s_ASideSlope",typeTitles[itype]));
      grCsideSlope->SetName(Form("%s_CSideSlope",typeTitles[itype]));
      grAsideSlope->SetTitle(Form("%s: A Side Slope",typeTitles[itype]));
      grCsideSlope->SetTitle(Form("%s: C side Slope",typeTitles[itype]));
      grAsideSlope->GetXaxis()->SetTitle("#phi");

      grChi2Aside->SetName(Form("%s_AsideChi2",typeTitles[itype]));
      grChi2Cside->SetName(Form("%s_CSideChi2",typeTitles[itype]));
      grChi2Aside->SetTitle(Form("%s: A Side Chi2",typeTitles[itype]));
      grChi2Cside->SetTitle(Form("%s: C side Chi2",typeTitles[itype]));
      grChi2Aside->GetXaxis()->SetTitle("#phi");
      grChi2Aside->GetYaxis()->SetTitle("#Chi^2");
      grChi2Cside->GetXaxis()->SetTitle("#phi");
      grChi2Cside->GetYaxis()->SetTitle("#Chi^2");

      grNpointsAside->SetName(Form("%s_AsideNpoints",typeTitles[itype]));
      grNpointsCside->SetName(Form("%s_CSideNpoints",typeTitles[itype]));
      grNpointsAside->SetTitle(Form("%s: A Side Npoints",typeTitles[itype]));
      grNpointsCside->SetTitle(Form("%s: C side Npoints",typeTitles[itype]));
      grNpointsAside->GetXaxis()->SetTitle("#phi");
      grNpointsAside->GetYaxis()->SetTitle("Number of points");
      grNpointsCside->GetXaxis()->SetTitle("#phi");
      grNpointsCside->GetYaxis()->SetTitle("Number of points");

      //Store Trending variables  {itype=0,1,2}={"TPConly","TPCconstrained","Combined"};
      if(itype == 0){
        qptShiftCombinedSinAside = vecFitA[2];  qptShiftCombinedSinCside = vecFitC[2];  qptShiftCombinedCosAside = vecFitA[1];  qptShiftCombinedCosCside = vecFitC[1];
        eqptShiftCombinedSinAside =   eFitA[2]; eqptShiftCombinedSinCside =   eFitC[2]; eqptShiftCombinedCosAside =   eFitA[1]; eqptShiftCombinedCosCside =   eFitC[1];
      }
      if(itype == 1){
        qptShiftTPCconstSinAside = vecFitA[2];  qptShiftTPCconstSinCside = vecFitC[2];  qptShiftTPCconstCosAside = vecFitA[1];  qptShiftTPCconstCosCside = vecFitC[1];
        eqptShiftTPCconstSinAside =   eFitA[2]; eqptShiftTPCconstSinCside =   eFitC[2]; eqptShiftTPCconstCosAside =   eFitA[1]; eqptShiftTPCconstCosCside =   eFitC[1];   
      }
      if(itype == 2){
        qptShiftTPConlySinAside  = vecFitA[2];  qptShiftTPConlySinCside  = vecFitC[2];  qptShiftTPConlyCosAside  = vecFitA[1];  qptShiftTPConlyCosCside  = vecFitC[1];
        eqptShiftTPConlySinAside  =   eFitA[2]; eqptShiftTPConlySinCside  =   eFitC[2]; eqptShiftTPConlyCosAside  =   eFitA[1]; eqptShiftTPConlyCosCside  =   eFitC[1];
      }


      grAsideSlope->GetYaxis()->SetTitle("#Slope (1/GeV)");;
      for (Int_t i=0; i<3;i++) eFitA[i]=fphi->GetParError(i);
      gPad->SetLogy(0);
      grAside->Draw("ap"); workCanvas->Update();
      TPaveStats *statfA = (TPaveStats*)grAside->FindObject("stats");
      grCside->Draw("ap"); workCanvas->Update();
      TPaveStats *statfC = (TPaveStats*)grCside->FindObject("stats");
      statfA->SetBorderSize(4);statfA->SetTextColor(2);
      statfC->SetBorderSize(4);statfC->SetTextColor(4);
      statfA->SetOptStat(1000000111);
      statfC->SetOptStat(1000000111);
      statfA->SetX1NDC(0.75); statfA->SetX2NDC(0.99);
      statfC->SetX1NDC(0.75); statfC->SetX2NDC(0.99);
      statfA->SetY1NDC(0.7); statfA->SetY2NDC(0.95);
      statfC->SetY1NDC(0.40); statfC->SetY2NDC(0.65);


      gStyle->SetOptTitle(0);
      TCanvas * canvasFit = new TCanvas(Form("canvas%s",  typeTitles[itype]), Form("canvas%s",  typeTitles[itype]),600,700);
      canvasFit->Divide(1,2);
      canvasFit->cd(1)->SetTicks(1);
      gPad->SetLogy(0);  
      gPad->SetTicky(2);gPad->SetTickx(1); 
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.08);gPad->SetRightMargin(0.31);
      hisAll->SetMinimum(0);

      hisAll->  GetXaxis()->SetRangeUser(-0.2,0.2);
      hisAside->GetXaxis()->SetRangeUser(-0.2,0.2);
      hisCside->GetXaxis()->SetRangeUser(-0.2,0.2);

      hisAll->Draw("err func");
      hisCside->Draw("err func same");
      hisAside->Draw("err func same");
      hisAll->Draw("err func same");
      canvasFit->cd(2)->SetTicks(1);
      gPad->SetLogy(0);
      gPad->SetTicky(2);gPad->SetTickx(1);
      gPad->SetTopMargin(0.05);
      gPad->SetLeftMargin(0.08);gPad->SetRightMargin(0.31);
      grAside->Draw("ap");
      grCside->Draw("p");
      canvasFit->cd(2);
      Double_t grBoth[36]={0};
      Double_t grBothE[36]={0};
      for (Int_t i=0; i<18; i++) {grBothE[i]= grAside->GetErrorY(i); grBothE[i+18]= grCside->GetErrorY(i); }
      Double_t medianErr = TMath::Median(36, grBothE);
      for (Int_t i=0; i<18; i++) {
        if (grAside->GetErrorY(i)>4*medianErr){
          grAside->GetY()[i]=0;
          grAside->SetPointError(i,0,0);
        }
        if (grCside->GetErrorY(i)>4*medianErr){
          grCside->GetY()[i]=0;
          grCside->SetPointError(i,0,0);
        }
        grBothE[i]= grAside->GetErrorY(i); 
        grBothE[i+18]= grCside->GetErrorY(i); 
        grBoth[i]= grAside->GetY()[i]; 
        grBoth[i+18]= grCside->GetY()[i];
      }

      Double_t meanB=TMath::Mean(36,  grBoth);
      Double_t meanError=TMath::Mean(36,  grBothE);
      Double_t rmsB=TMath::RMS(36,  grBoth);
      Double_t meanA=TMath::Mean(18,  grAside->GetY());
      Double_t rmsA=TMath::RMS(18,  grAside->GetY());
      Double_t meanC=TMath::Mean(18,  grCside->GetY());
      Double_t rmsC=TMath::RMS(18,  grCside->GetY());
      rmsB=TMath::Sqrt(TMath::Abs(rmsB*rmsB-meanError*meanError));
      rmsA=TMath::Sqrt(TMath::Abs(rmsA*rmsA-meanError*meanError));
      rmsC=TMath::Sqrt(TMath::Abs(rmsC*rmsC-meanError*meanError));
      TLegend * legend = new TLegend(.08,.8, .69, 0.99, Form("%s tracks Run=%d: Bz=%1.1f T",typeTitles[itype],(Int_t)runNumber,Bz/10.));
      legend->AddEntry(hisAll,  Form("Both sides Mean=%1.4f RMS=%1.4f",meanB,rmsB));
      legend->AddEntry(hisAside,Form("A side     Mean=%1.4f RMS=%1.4f",meanA,rmsA));
      legend->AddEntry(hisCside,Form("C side     Mean=%1.4f RMS=%1.4f",meanC,rmsC));
      legend->Draw();

      canvasFit->Write();
      //   canvasFit->SaveAs(Form("RawPowerLawFit_%d.eps", itype));
      canvasFit->SaveAs(Form("RawPowerLawFit_%s.png", typeName[itype]));
      Int_t hentries = hisAll->GetEntries();
      //  pcstream->GetFile()->cd();
      grAside->Write();
      grCside->Write();
      grAsideSlope->Write();
      grCsideSlope->Write();
      hisAll->Write();
      hisAside->Write();
      hisCside->Write();
      //grSlopeA[itype]= (TGraphErrors*)grAsideSlope->Clone();
      //grSlopeC[itype]= (TGraphErrors*)grCsideSlope->Clone();

      grChi2Aside->Write();
      grChi2Cside->Write();
      grNpointsAside->Write();
      grNpointsCside->Write();

      if(!fMakePeriod){      
        //     TBranch *brItype  = OutTree->Branch("itype",&itype);
        //     brItype->Fill();
        TBranch *brfAll   = OutTree->Branch( Form("fAll_%s.",typeName[itype]) ,&vecFitAll);
        brfAll->Fill();
        TBranch *brfAside = OutTree->Branch(Form("fAside_%s.",typeName[itype]),&vecFitAside);
        brfAside->Fill();
        TBranch *brfCside = OutTree->Branch(Form("fCside_%s.",typeName[itype]),&vecFitCside);
        brfCside->Fill();
        TBranch *breAll   = OutTree->Branch(Form("eAll_%s.",typeName[itype]),&vecFitAll);
        breAll->Fill();
        TBranch *breAside = OutTree->Branch(Form("eAside_%s.",typeName[itype]),&vecFitAside);
        breAside->Fill();
        TBranch *breCside = OutTree->Branch(Form("eCside_%s.",typeName[itype]),&vecFitCside);
        breCside->Fill();
        TBranch *brentries = OutTree->Branch(Form("entries_%s",typeName[itype]),&hentries);
        brentries->Fill();
        TBranch *brchi2All = OutTree->Branch(Form("chi2All_%s",typeName[itype]),&chi2All);
        brchi2All->Fill();
        TBranch *brchi2Aside = OutTree->Branch(Form("chi2Aside_%s",typeName[itype]),&chi2Aside);
        brchi2Aside->Fill();
        TBranch *brchi2Cside = OutTree->Branch(Form("chi2Cside_%s",typeName[itype]),&chi2Cside);
        brchi2Cside->Fill();
        TBranch *brnPointsAll = OutTree->Branch(Form("nPointsAll_%s",typeName[itype]),&nPointsAll);
        brnPointsAll->Fill();
        TBranch *brnPointsAside = OutTree->Branch(Form("nPointsAside_%s",typeName[itype]),&nPointsAside);
        brnPointsAside->Fill();
        TBranch *brnPointsCside = OutTree->Branch(Form("nPointsCside_%s",typeName[itype]),&nPointsCside);
        brnPointsCside->Fill();

        TBranch *brgrAside = OutTree->Branch(Form("grAside_%s",typeName[itype]),&grAside);
        brgrAside->Fill();
        TBranch *brgrCside = OutTree->Branch(Form("grCside_%s",typeName[itype]),&grCside);
        brgrCside->Fill();
        TBranch *brgrChi2Aside = OutTree->Branch(Form("grChi2Aside_%s",typeName[itype]),&grChi2Aside);
        brgrChi2Aside->Fill();
        TBranch *brgrChi2Cside = OutTree->Branch(Form("grChi2Cside_%s",typeName[itype]),&grChi2Cside);
        brgrChi2Cside->Fill();
        TBranch *brgrNpointsAside = OutTree->Branch(Form("grNpointsAside_%s",typeName[itype]),&grNpointsAside);
        brgrNpointsAside->Fill();
        TBranch *brgrNpointsCside = OutTree->Branch(Form("grNpointsCside_%s",typeName[itype]),&grNpointsCside);
        brgrNpointsCside->Fill();
        TBranch *brfA = OutTree->Branch(Form("vecFitA_%s.",typeName[itype]),&vecFitA);
        brfA->Fill();
        TBranch *brfC = OutTree->Branch(Form("vecFitC_%s.",typeName[itype]),&vecFitC);
        brfC->Fill();
        TBranch *breA = OutTree->Branch(Form("eFitA_%s.",typeName[itype]),&eFitA);
        breA->Fill();
        TBranch *breC = OutTree->Branch(Form("eFitC_%s.",typeName[itype]),&eFitC);
        breC->Fill();
      }
    }

  }

  //Store Trending
  TBranch *brqptShiftCombinedSinAside  = OutTree->Branch("qptShiftCombinedSinAside",&qptShiftCombinedSinAside);
  brqptShiftCombinedSinAside->Fill();
  TBranch *breqptShiftCombinedSinAside = OutTree->Branch("eqptShiftCombinedSinAside",&eqptShiftCombinedSinAside);
  breqptShiftCombinedSinAside->Fill();
  TBranch *brqptShiftCombinedSinCside  = OutTree->Branch("qptShiftCombinedSinCside",&qptShiftCombinedSinCside);
  brqptShiftCombinedSinCside->Fill();
  TBranch *breqptShiftCombinedSinCside = OutTree->Branch("eqptShiftCombinedSinCside",&eqptShiftCombinedSinCside);
  breqptShiftCombinedSinCside->Fill();
  TBranch *brqptShiftCombinedCosAside  = OutTree->Branch("qptShiftCombinedCosAside",&qptShiftCombinedCosAside);
  brqptShiftCombinedCosAside->Fill();
  TBranch *breqptShiftCombinedCosAside = OutTree->Branch("eqptShiftCombinedCosAside",&eqptShiftCombinedCosAside);
  breqptShiftCombinedCosAside->Fill();
  TBranch *brqptShiftCombinedCosCside  = OutTree->Branch("qptShiftCombinedCosCside",&qptShiftCombinedCosCside);
  brqptShiftCombinedCosCside->Fill();
  TBranch *breqptShiftCombinedCosCside = OutTree->Branch("eqptShiftCombinedCosCside",&eqptShiftCombinedCosCside);
  breqptShiftCombinedCosCside->Fill();
  TBranch *brqptShiftTPCconstSinAside  = OutTree->Branch("qptShiftTPCconstSinAside",&qptShiftTPCconstSinAside);
  brqptShiftTPCconstSinAside->Fill();
  TBranch *breqptShiftTPCconstSinAside = OutTree->Branch("eqptShiftTPCconstSinAside",&eqptShiftTPCconstSinAside);
  breqptShiftTPCconstSinAside->Fill();
  TBranch *brqptShiftTPCconstSinCside  = OutTree->Branch("qptShiftTPCconstSinCside",&qptShiftTPCconstSinCside);
  brqptShiftTPCconstSinCside->Fill();
  TBranch *breqptShiftTPCconstSinCside = OutTree->Branch("eqptShiftTPCconstSinCside",&eqptShiftTPCconstSinCside);
  breqptShiftTPCconstSinCside->Fill();
  TBranch *brqptShiftTPCconstCosAside  = OutTree->Branch("qptShiftTPCconstCosAside",&qptShiftTPCconstCosAside);
  brqptShiftTPCconstCosAside->Fill();
  TBranch *breqptShiftTPCconstCosAside = OutTree->Branch("eqptShiftTPCconstCosAside",&eqptShiftTPCconstCosAside);
  breqptShiftTPCconstCosAside->Fill();
  TBranch *brqptShiftTPCconstCosCside  = OutTree->Branch("qptShiftTPCconstCosCside",&qptShiftTPCconstCosCside);
  brqptShiftTPCconstCosCside->Fill();
  TBranch *breqptShiftTPCconstCosCside = OutTree->Branch("eqptShiftTPCconstCosCside",&eqptShiftTPCconstCosCside);
  breqptShiftTPCconstCosCside->Fill();
  TBranch *brqptShiftTPConlySinAside   = OutTree->Branch("qptShiftTPConlySinAside",&qptShiftTPConlySinAside);
  brqptShiftTPConlySinAside->Fill();
  TBranch *breqptShiftTPConlySinAside  = OutTree->Branch("eqptShiftTPConlySinAside",&eqptShiftTPConlySinAside);
  breqptShiftTPConlySinAside->Fill();
  TBranch *brqptShiftTPConlySinCside   = OutTree->Branch("qptShiftTPConlySinCside",&qptShiftTPConlySinCside);
  brqptShiftTPConlySinCside->Fill();
  TBranch *breqptShiftTPConlySinCside  = OutTree->Branch("eqptShiftTPConlySinCside",&eqptShiftTPConlySinCside);
  breqptShiftTPConlySinCside->Fill();
  TBranch *brqptShiftTPConlyCosAside   = OutTree->Branch("qptShiftTPConlyCosAside",&qptShiftTPConlyCosAside);
  brqptShiftTPConlyCosAside->Fill();
  TBranch *breqptShiftTPConlyCosAside  = OutTree->Branch("eqptShiftTPConlyCosAside",&eqptShiftTPConlyCosAside);
  breqptShiftTPConlyCosAside->Fill();
  TBranch *brqptShiftTPConlyCosCside   = OutTree->Branch("qptShiftTPConlyCosCside",&qptShiftTPConlyCosCside);
  brqptShiftTPConlyCosCside->Fill();
  TBranch *breqptShiftTPConlyCosCside  = OutTree->Branch("eqptShiftTPConlyCosCside",&eqptShiftTPConlyCosCside);
  breqptShiftTPConlyCosCside->Fill();


}

void AliHighPtTreeAnalysis::MakeTPCITSMatchingEff(){
  Double_t e0 = 0.;
  Double_t e1 = 0.;
  Double_t c12 = 0.;

  Double_t effTPCITS_TPCAside = 0.;
  Double_t eEffTPCITS_TPCAside = 0.;
  Double_t countsTPCITS_TPCAside = 0;
  Double_t countsTPC_TPCAside = 0;

  Double_t effTPCITS_TPCCside = 0.;
  Double_t eEffTPCITS_TPCCside = 0.;
  Double_t countsTPCITS_TPCCside = 0;
  Double_t countsTPC_TPCCside = 0;


  TH3D *h1TPC    = hphi_vs_eta_pT_cutTPC;
  TH3D *h1TPCITS = hphi_vs_eta_pT_cutTPCITS;
  if(!h1TPC) return;
  if(!h1TPCITS) return;

  // pT above 3 GeV/c
  h1TPC->GetZaxis()->SetRangeUser(fPtCut,100.);
  h1TPCITS->GetZaxis()->SetRangeUser(fPtCut,100.);

  // TPC A side
  h1TPC->GetYaxis()->SetRangeUser(0.,0.999);
  h1TPCITS->GetYaxis()->SetRangeUser(0.,0.999);

  countsTPCITS_TPCAside = h1TPCITS->Integral();
  countsTPC_TPCAside = h1TPC->Integral();

  if(countsTPC_TPCAside) { 
    e0 = TMath::Sqrt(h1TPCITS->Integral());
    e1 = TMath::Sqrt(h1TPC->Integral());
    c12 = countsTPC_TPCAside*countsTPC_TPCAside;
    effTPCITS_TPCAside = (Double_t)h1TPCITS->Integral() / (Double_t)h1TPC->Integral();
    eEffTPCITS_TPCAside = TMath::Sqrt((e0*e0*countsTPC_TPCAside*countsTPC_TPCAside + e1*e1*countsTPCITS_TPCAside*countsTPCITS_TPCAside)/(c12*c12));
  } 

  // TPC C side
  h1TPC->GetYaxis()->SetRangeUser(-1.0,-0.0001);
  h1TPCITS->GetYaxis()->SetRangeUser(-1.0,-0.0001);

  countsTPCITS_TPCCside = h1TPCITS->Integral();
  countsTPC_TPCCside = h1TPC->Integral();

  if(countsTPC_TPCCside) { 
    e0 = TMath::Sqrt(h1TPCITS->Integral());
    e1 = TMath::Sqrt(h1TPC->Integral());
    c12 = countsTPC_TPCCside*countsTPC_TPCCside;
    effTPCITS_TPCCside = (Double_t)h1TPCITS->Integral() / (Double_t)h1TPC->Integral();
    eEffTPCITS_TPCCside = TMath::Sqrt((e0*e0*countsTPC_TPCCside*countsTPC_TPCCside + e1*e1*countsTPCITS_TPCCside*countsTPCITS_TPCCside)/(c12*c12));
  } 

  // set back full range
  h1TPC->GetZaxis()->SetRange(1,h1TPC->GetZaxis()->GetNbins());
  h1TPCITS->GetZaxis()->SetRange(1,h1TPCITS->GetZaxis()->GetNbins());
  h1TPC->GetYaxis()->SetRange(1,h1TPC->GetYaxis()->GetNbins());
  h1TPCITS->GetYaxis()->SetRange(1,h1TPCITS->GetYaxis()->GetNbins());

  if(!fMakePeriod){      
    //     TBranch *brItype  = OutTree->Branch("itype",&itype);
    //     brItype->Fill();
    TBranch *brCountsTPC_TPCAside = OutTree->Branch("countsTPC_TPCAside",&countsTPC_TPCAside);
    brCountsTPC_TPCAside->Fill();
    TBranch *brCountsTPC_TPCCside = OutTree->Branch("countsTPC_TPCCside",&countsTPC_TPCCside);
    brCountsTPC_TPCCside->Fill();
    TBranch *brCountsTPCITS_TPCAside = OutTree->Branch("countsTPCITS_TPCAside",&countsTPCITS_TPCAside);
    brCountsTPCITS_TPCAside->Fill();
    TBranch *brCountsTPCITS_TPCCside = OutTree->Branch("countsTPCITS_TPCCside",&countsTPCITS_TPCCside);
    brCountsTPCITS_TPCCside->Fill();
    TBranch *breffTPCITS_TPCAside = OutTree->Branch("effTPCITS_TPCAside",&effTPCITS_TPCAside);
    breffTPCITS_TPCAside->Fill();
    TBranch *breffTPCITS_TPCCside = OutTree->Branch("effTPCITS_TPCCside",&effTPCITS_TPCCside);
    breffTPCITS_TPCCside->Fill();
    TBranch *breEffTPCITS_TPCAside = OutTree->Branch("eEffTPCITS_TPCAside",&eEffTPCITS_TPCAside);
    breEffTPCITS_TPCAside->Fill();
    TBranch *breEffTPCITS_TPCCside = OutTree->Branch("eEffTPCITS_TPCCside",&eEffTPCITS_TPCCside);
    breEffTPCITS_TPCCside->Fill();
  }
}

void AliHighPtTreeAnalysis::MakeK0trends(){
  Double_t shiftK0sPullNeg  = 0;  Double_t shiftK0sPullPos  = 0;
  Double_t eShiftK0sPullNeg = 0;  Double_t eShiftK0sPullPos = 0;
  Double_t sigmaK0sPullNeg  = 0;  Double_t sigmaK0sPullPos  = 0;
  Double_t eSigmaK0sPullNeg = 0;  Double_t eSigmaK0sPullPos = 0;

  Double_t shiftK0sResNeg  = 0;   Double_t shiftK0sResPos  = 0;
  Double_t eShiftK0sResNeg = 0;   Double_t eShiftK0sResPos = 0;
  Double_t sigmaK0sResNeg  = 0;   Double_t sigmaK0sResPos  = 0;
  Double_t eSigmaK0sResNeg = 0;   Double_t eSigmaK0sResPos = 0;

  Int_t countsK0sNeg = 0;         Int_t countsK0sPos = 0;

  Double_t  shiftK0sPullNegHigh1pt = -999.; Double_t  shiftK0sPullPosHigh1pt = -999.;
  Double_t eShiftK0sPullNegHigh1pt = -999.; Double_t eShiftK0sPullPosHigh1pt = -999.;
  Double_t  shiftK0sPullNegLow1pt  = -999.; Double_t  shiftK0sPullPosLow1pt  = -999.;
  Double_t eShiftK0sPullNegLow1pt  = -999.; Double_t eShiftK0sPullPosLow1pt  = -999.;

  Double_t  sigmaK0sPullNegHigh1pt = -999.; Double_t  sigmaK0sPullPosHigh1pt = -999.;
  Double_t eSigmaK0sPullNegHigh1pt = -999.; Double_t eSigmaK0sPullPosHigh1pt = -999.;
  Double_t  sigmaK0sPullNegLow1pt  = -999.; Double_t  sigmaK0sPullPosLow1pt  = -999.;
  Double_t eSigmaK0sPullNegLow1pt  = -999.; Double_t eSigmaK0sPullPosLow1pt  = -999.;

  Double_t  shiftK0sResNegHigh1pt = -999.; Double_t  shiftK0sResPosHigh1pt = -999.;
  Double_t eShiftK0sResNegHigh1pt = -999.; Double_t eShiftK0sResPosHigh1pt = -999.;
  Double_t  shiftK0sResNegLow1pt  = -999.; Double_t  shiftK0sResPosLow1pt  = -999.;
  Double_t eShiftK0sResNegLow1pt  = -999.; Double_t eShiftK0sResPosLow1pt  = -999.;

  Double_t  sigmaK0sResNegHigh1pt = -999.; Double_t  sigmaK0sResPosHigh1pt = -999.;
  Double_t eSigmaK0sResNegHigh1pt = -999.; Double_t eSigmaK0sResPosHigh1pt = -999.;
  Double_t  sigmaK0sResNegLow1pt  = -999.; Double_t  sigmaK0sResPosLow1pt  = -999.;
  Double_t eSigmaK0sResNegLow1pt  = -999.; Double_t eSigmaK0sResPosLow1pt  = -999.;

  if(fMakeV0s){
    Double_t low1pt  = 0.;
    Double_t high1pt = 1.;
    TF1 *fLinearFitK0sShift = new TF1("fLinearFitK0sShift","[0] + [1]*x",0.,1.);
    TF1 *fLinearFitK0sSigma = new TF1("fLinearFitK0sSigma","[0] + [1]*x",0.,1.);

    Bool_t bfit(kFALSE);

    fLinearFitK0sShift->SetParameters(0,1);
    fLinearFitK0sSigma->SetParameters(0,1);
    bfit = GetK0TrendFitFunction(fLinearFitK0sShift,fLinearFitK0sSigma,0,-1);
    if(bfit){
      shiftK0sResNegLow1pt   = fLinearFitK0sShift->Eval(low1pt);
      eShiftK0sResNegLow1pt  = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(low1pt*fLinearFitK0sShift->GetParError(1),2) );
      shiftK0sResNegHigh1pt  = fLinearFitK0sShift->Eval(high1pt);
      eShiftK0sResNegHigh1pt = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(high1pt*fLinearFitK0sShift->GetParError(1),2) );

      sigmaK0sResNegLow1pt   = fLinearFitK0sSigma->Eval(low1pt);
      eSigmaK0sResNegLow1pt  = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(low1pt*fLinearFitK0sSigma->GetParError(1),2) );
      sigmaK0sResNegHigh1pt  = fLinearFitK0sSigma->Eval(high1pt);
      eSigmaK0sResNegHigh1pt = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(high1pt*fLinearFitK0sSigma->GetParError(1),2) );  
      bfit = kFALSE;
    }
    fLinearFitK0sShift->SetParameters(0,1);
    fLinearFitK0sSigma->SetParameters(0,1);
    bfit = GetK0TrendFitFunction(fLinearFitK0sShift,fLinearFitK0sSigma,0, 1);
    if(bfit){
      shiftK0sResPosLow1pt   = fLinearFitK0sShift->Eval(low1pt);
      eShiftK0sResPosLow1pt  = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(low1pt*fLinearFitK0sShift->GetParError(1),2) );
      shiftK0sResPosHigh1pt  = fLinearFitK0sShift->Eval(high1pt);
      eShiftK0sResPosHigh1pt = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(high1pt*fLinearFitK0sShift->GetParError(1),2) );

      sigmaK0sResPosLow1pt   = fLinearFitK0sSigma->Eval(low1pt);
      eSigmaK0sResPosLow1pt  = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(low1pt*fLinearFitK0sSigma->GetParError(1),2) );
      sigmaK0sResPosHigh1pt  = fLinearFitK0sSigma->Eval(high1pt);
      eSigmaK0sResPosHigh1pt = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(high1pt*fLinearFitK0sSigma->GetParError(1),2) ); 
      bfit = kFALSE;
    }  
    fLinearFitK0sShift->SetParameters(0,1);
    fLinearFitK0sSigma->SetParameters(0,1);
    bfit = GetK0TrendFitFunction(fLinearFitK0sShift,fLinearFitK0sSigma,1,-1);
    if(bfit){
      shiftK0sPullNegLow1pt   = fLinearFitK0sShift->Eval(low1pt);
      eShiftK0sPullNegLow1pt  = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(low1pt*fLinearFitK0sShift->GetParError(1),2) );
      shiftK0sPullNegHigh1pt  = fLinearFitK0sShift->Eval(high1pt);
      eShiftK0sPullNegHigh1pt = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(high1pt*fLinearFitK0sShift->GetParError(1),2) );

      sigmaK0sPullNegLow1pt   = fLinearFitK0sSigma->Eval(low1pt);
      eSigmaK0sPullNegLow1pt  = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(low1pt*fLinearFitK0sSigma->GetParError(1),2) );
      sigmaK0sPullNegHigh1pt  = fLinearFitK0sSigma->Eval(high1pt);
      eSigmaK0sPullNegHigh1pt = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(high1pt*fLinearFitK0sSigma->GetParError(1),2) );
      bfit = kFALSE;
    }
    fLinearFitK0sShift->SetParameters(0,1);
    fLinearFitK0sSigma->SetParameters(0,1);
    bfit = GetK0TrendFitFunction(fLinearFitK0sShift,fLinearFitK0sSigma,1, 1);
    if(bfit){
      shiftK0sPullPosLow1pt   = fLinearFitK0sShift->Eval(low1pt);
      eShiftK0sPullPosLow1pt  = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(low1pt*fLinearFitK0sShift->GetParError(1),2) );
      shiftK0sPullPosHigh1pt  = fLinearFitK0sShift->Eval(high1pt);
      eShiftK0sPullPosHigh1pt = sqrt( pow(fLinearFitK0sShift->GetParError(0),2) + pow(high1pt*fLinearFitK0sShift->GetParError(1),2) );

      sigmaK0sPullPosLow1pt   = fLinearFitK0sSigma->Eval(low1pt);
      eSigmaK0sPullPosLow1pt  = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(low1pt*fLinearFitK0sSigma->GetParError(1),2) );
      sigmaK0sPullPosHigh1pt  = fLinearFitK0sSigma->Eval(high1pt);
      eSigmaK0sPullPosHigh1pt = sqrt( pow(fLinearFitK0sSigma->GetParError(0),2) + pow(high1pt*fLinearFitK0sSigma->GetParError(1),2) ); 
      bfit = kFALSE;
    }
  }

  if(fMakeV0s){
    // fit function
    TF1 * fmass = new TF1("fmass","[0]+[1]*x+[2]*TMath::Gaus(x,[3],[4])",-0.5,0.5); 

    TH3D *h1K0sPull = hK0sPull_vs_alpha_1pT_neg;
    TH3D *h2K0sPull = hK0sPull_vs_alpha_1pT_pos;
    TH3D *h1K0sRes =  hK0sRes_vs_alpha_1pT_neg;
    TH3D *h2K0sRes =  hK0sRes_vs_alpha_1pT_pos;
    if(!h1K0sRes || !h2K0sRes || !h1K0sPull || !h2K0sPull) return; 

    // pT above 3 GeV/c
    h1K0sPull->GetXaxis()->SetRangeUser(0.0,1./fPtCut);
    h2K0sPull->GetXaxis()->SetRangeUser(0.0,1./fPtCut);
    countsK0sNeg = h1K0sPull->Integral();
    countsK0sPos = h2K0sPull->Integral();

    TH1D *h1K0sPullproj = (TH1D*) h1K0sPull->Project3D("z");
    TH1D *h2K0sPullproj = (TH1D*) h2K0sPull->Project3D("z");
    if(!h1K0sPullproj || !h2K0sPullproj) return;

    Double_t y1 = h1K0sPullproj->GetBinContent(5);
    Double_t y2 = h1K0sPullproj->GetBinContent(95);
    fmass->SetParameter(0,(y1+y2)*0.5);
    fmass->SetParameter(1,(y2-y1)/20.);
    fmass->SetParameter(2,h1K0sPullproj->GetMaximum());
    fmass->SetParameter(3,0);
    fmass->SetParameter(4,h1K0sPullproj->GetRMS());
    h1K0sPullproj->Fit(fmass,"Q"); 
    h1K0sPullproj->Fit(fmass,"Q"); 
    //h1K0sPullproj->Draw("e"); 
    shiftK0sPullNeg  = fmass->GetParameter(3);
    eShiftK0sPullNeg = fmass->GetParError(3);    
    sigmaK0sPullNeg  = TMath::Abs(fmass->GetParameter(4));
    eSigmaK0sPullNeg = TMath::Abs(fmass->GetParError(4));    

    y1 = h2K0sPullproj->GetBinContent(5);
    y2 = h2K0sPullproj->GetBinContent(95);
    fmass->SetParameter(0,(y1+y2)*0.5);
    fmass->SetParameter(1,(y2-y1)/20.);
    fmass->SetParameter(2,h2K0sPullproj->GetMaximum());
    fmass->SetParameter(3,0);
    fmass->SetParameter(4,h2K0sPullproj->GetRMS());
    h2K0sPullproj->Fit(fmass,"Q"); 
    h2K0sPullproj->Fit(fmass,"Q"); 
    //h2K0sPullproj->Draw("e"); 
    shiftK0sPullPos  = fmass->GetParameter(3);
    eShiftK0sPullPos = fmass->GetParError(3);    
    sigmaK0sPullPos  = TMath::Abs(fmass->GetParameter(4));
    eSigmaK0sPullPos = TMath::Abs(fmass->GetParError(4));    
    //
    if(h1K0sPullproj) delete h1K0sPullproj;  
    if(h2K0sPullproj) delete h2K0sPullproj;  
    // set back full range
    h1K0sPull->GetXaxis()->SetRange(1,h1K0sPull->GetXaxis()->GetNbins());
    h2K0sPull->GetXaxis()->SetRange(1,h2K0sPull->GetXaxis()->GetNbins());

    // pT above 3 GeV/c
    h1K0sRes->GetXaxis()->SetRangeUser(0.0,1./fPtCut);
    h2K0sRes->GetXaxis()->SetRangeUser(0.0,1./fPtCut);

    countsK0sNeg = h1K0sRes->Integral();
    countsK0sPos = h2K0sRes->Integral();

    TH1D *h1K0sResproj = (TH1D*) h1K0sRes->Project3D("z");
    TH1D *h2K0sResproj = (TH1D*) h2K0sRes->Project3D("z");
    if(!h1K0sResproj || !h2K0sResproj) return;

    y1 = h1K0sResproj->GetBinContent(5);
    y2 = h1K0sResproj->GetBinContent(95);
    fmass->SetParameter(0,(y1+y2)*0.5);
    fmass->SetParameter(1,(y2-y1)/20.);
    fmass->SetParameter(2,h1K0sResproj->GetMaximum());
    fmass->SetParameter(3,0);
    fmass->SetParameter(4,h1K0sResproj->GetRMS());
    h1K0sResproj->Fit(fmass,"Q"); 
    h1K0sResproj->Fit(fmass,"Q"); 
    //h1K0sResproj->Draw("e"); 
    shiftK0sResNeg  = fmass->GetParameter(3);
    eShiftK0sResNeg = fmass->GetParError(3);    
    sigmaK0sResNeg  = TMath::Abs(fmass->GetParameter(4));
    eSigmaK0sResNeg = TMath::Abs(fmass->GetParError(4));    

    y1 = h2K0sResproj->GetBinContent(5);
    y2 = h2K0sResproj->GetBinContent(95);
    fmass->SetParameter(0,(y1+y2)*0.5);
    fmass->SetParameter(1,(y2-y1)/20.);
    fmass->SetParameter(2,h2K0sResproj->GetMaximum());
    fmass->SetParameter(3,0);
    fmass->SetParameter(4,h2K0sResproj->GetRMS());
    h2K0sResproj->Fit(fmass,"Q"); 
    h2K0sResproj->Fit(fmass,"Q"); 
    //h2K0sResproj->Draw("e"); 
    shiftK0sResPos  = fmass->GetParameter(3);
    eShiftK0sResPos = fmass->GetParError(3);    
    sigmaK0sResPos  = TMath::Abs(fmass->GetParameter(4));
    eSigmaK0sResPos = TMath::Abs(fmass->GetParError(4));    
    if(h1K0sResproj) delete h1K0sResproj;  
    if(h2K0sResproj) delete h2K0sResproj;  

    // set back full range
    h1K0sRes->GetXaxis()->SetRange(1,h1K0sRes->GetXaxis()->GetNbins());
    h2K0sRes->GetXaxis()->SetRange(1,h2K0sRes->GetXaxis()->GetNbins());

  }

  TBranch *brcountsK0sNeg = OutTree->Branch("countsK0sNeg",&countsK0sNeg);
  brcountsK0sNeg->Fill();  
  TBranch *brcountsK0sPos = OutTree->Branch("countsK0sPos",&countsK0sPos);
  brcountsK0sPos->Fill();
  TBranch *brsigmaK0sPullNeg = OutTree->Branch("sigmaK0sPullNeg",&sigmaK0sPullNeg);
  brsigmaK0sPullNeg->Fill();
  TBranch *brsigmaK0sPullPos = OutTree->Branch("sigmaK0sPullPos",&sigmaK0sPullPos);
  brsigmaK0sPullPos->Fill();
  TBranch *breSigmaK0sPullNeg = OutTree->Branch("eSigmaK0sPullNeg",&eSigmaK0sPullNeg);
  breSigmaK0sPullNeg->Fill();
  TBranch *breSigmaK0sPullPos = OutTree->Branch("eSigmaK0sPullPos",&eSigmaK0sPullPos);
  breSigmaK0sPullPos->Fill();
  TBranch *brshiftK0sPullNeg = OutTree->Branch("shiftK0sPullNeg",&shiftK0sPullNeg);
  brshiftK0sPullNeg->Fill();
  TBranch *brshiftK0sPullPos = OutTree->Branch("shiftK0sPullPos",&shiftK0sPullPos);
  brshiftK0sPullPos->Fill();
  TBranch *breShiftK0sPullNeg = OutTree->Branch("eShiftK0sPullNeg",&eShiftK0sPullNeg);
  breShiftK0sPullNeg->Fill();
  TBranch *breShiftK0sPullPos = OutTree->Branch("eShiftK0sPullPos",&eShiftK0sPullPos);
  breShiftK0sPullPos->Fill();
  TBranch *brsigmaK0sResNeg = OutTree->Branch("sigmaK0sResNeg",&sigmaK0sResNeg);
  brsigmaK0sResNeg->Fill();
  TBranch *brsigmaK0sResPos = OutTree->Branch("sigmaK0sResPos",&sigmaK0sResPos);
  brsigmaK0sResPos->Fill();
  TBranch *breSigmaK0sResNeg = OutTree->Branch("eSigmaK0sResNeg",&eSigmaK0sResNeg);
  breSigmaK0sResNeg->Fill();
  TBranch *breSigmaK0sResPos = OutTree->Branch("eSigmaK0sResPos",&eSigmaK0sResPos);
  breSigmaK0sResPos->Fill();
  TBranch *brshiftK0sResNeg = OutTree->Branch("shiftK0sResNeg",&shiftK0sResNeg);
  brshiftK0sResNeg->Fill();
  TBranch *brshiftK0sResPos = OutTree->Branch("shiftK0sResPos",&shiftK0sResPos);
  brshiftK0sResPos->Fill();
  TBranch *breShiftK0sResNeg = OutTree->Branch("eShiftK0sResNeg",&eShiftK0sResNeg);
  breShiftK0sResNeg->Fill();
  TBranch *breShiftK0sResPos = OutTree->Branch("eShiftK0sResPos",&eShiftK0sResPos);
  breShiftK0sResPos->Fill();
  ////
  //// New Trending variables
  ////
  TBranch *brshiftK0sPullNegHigh1pt = OutTree->Branch("shiftK0sPullNegHigh1pt",&shiftK0sPullNegHigh1pt);
  brshiftK0sPullNegHigh1pt->Fill();
  TBranch *brshiftK0sPullPosHigh1pt = OutTree->Branch("shiftK0sPullPosHigh1pt",&shiftK0sPullPosHigh1pt);
  brshiftK0sPullPosHigh1pt->Fill();
  TBranch *breShiftK0sPullNegHigh1pt = OutTree->Branch("eShiftK0sPullNegHigh1pt",&eShiftK0sPullNegHigh1pt);
  breShiftK0sPullNegHigh1pt->Fill();
  TBranch *breShiftK0sPullPosHigh1pt = OutTree->Branch("eShiftK0sPullPosHigh1pt",&eShiftK0sPullPosHigh1pt);
  breShiftK0sPullPosHigh1pt->Fill();
  TBranch *brshiftK0sPullNegLow1pt = OutTree->Branch("shiftK0sPullNegLow1pt",&shiftK0sPullNegLow1pt);
  brshiftK0sPullNegLow1pt->Fill();
  TBranch *brshiftK0sPullPosLow1pt = OutTree->Branch("shiftK0sPullPosLow1pt",&shiftK0sPullPosLow1pt);
  brshiftK0sPullPosLow1pt->Fill();
  TBranch *breShiftK0sPullNegLow1pt = OutTree->Branch("eShiftK0sPullNegLow1pt",&eShiftK0sPullNegLow1pt);
  breShiftK0sPullNegLow1pt->Fill();
  TBranch *breShiftK0sPullPosLow1pt = OutTree->Branch("eShiftK0sPullPosLow1pt",&eShiftK0sPullPosLow1pt);
  breShiftK0sPullPosLow1pt->Fill();
  TBranch *brsigmaK0sPullNegHigh1pt = OutTree->Branch("sigmaK0sPullNegHigh1pt",&sigmaK0sPullNegHigh1pt);
  brsigmaK0sPullNegHigh1pt->Fill();
  TBranch *brsigmaK0sPullPosHigh1pt = OutTree->Branch("sigmaK0sPullPosHigh1pt",&sigmaK0sPullPosHigh1pt);
  brsigmaK0sPullPosHigh1pt->Fill();
  TBranch *breSigmaK0sPullNegHigh1pt = OutTree->Branch("eSigmaK0sPullNegHigh1pt",&eSigmaK0sPullNegHigh1pt);
  breSigmaK0sPullNegHigh1pt->Fill();
  TBranch *breSigmaK0sPullPosHigh1pt = OutTree->Branch("eSigmaK0sPullPosHigh1pt",&eSigmaK0sPullPosHigh1pt);
  breSigmaK0sPullPosHigh1pt->Fill();
  TBranch *brsigmaK0sPullNegLow1pt = OutTree->Branch("sigmaK0sPullNegLow1pt",&sigmaK0sPullNegLow1pt);
  brsigmaK0sPullNegLow1pt->Fill();
  TBranch *brsigmaK0sPullPosLow1pt = OutTree->Branch("sigmaK0sPullPosLow1pt",&sigmaK0sPullPosLow1pt);
  brsigmaK0sPullPosLow1pt->Fill();
  TBranch *breSigmaK0sPullNegLow1pt = OutTree->Branch("eSigmaK0sPullNegLow1pt",&eSigmaK0sPullNegLow1pt);
  breSigmaK0sPullNegLow1pt->Fill();
  TBranch *breSigmaK0sPullPosLow1pt = OutTree->Branch("eSigmaK0sPullPosLow1pt",&eSigmaK0sPullPosLow1pt);
  breSigmaK0sPullPosLow1pt->Fill();
  TBranch *brshiftK0sResNegHigh1pt = OutTree->Branch("shiftK0sResNegHigh1pt",&shiftK0sResNegHigh1pt);
  brshiftK0sResNegHigh1pt->Fill();
  TBranch *brshiftK0sResPosHigh1pt = OutTree->Branch("shiftK0sResPosHigh1pt",&shiftK0sResPosHigh1pt);
  brshiftK0sResPosHigh1pt->Fill();
  TBranch *breShiftK0sResNegHigh1pt = OutTree->Branch("eShiftK0sResNegHigh1pt",&eShiftK0sResNegHigh1pt);
  breShiftK0sResNegHigh1pt->Fill();
  TBranch *breShiftK0sResPosHigh1pt = OutTree->Branch("eShiftK0sResPosHigh1pt",&eShiftK0sResPosHigh1pt);
  breShiftK0sResPosHigh1pt->Fill();
  TBranch *brshiftK0sResNegLow1pt = OutTree->Branch("shiftK0sResNegLow1pt",&shiftK0sResNegLow1pt);
  brshiftK0sResNegLow1pt->Fill();
  TBranch *brshiftK0sResPosLow1pt = OutTree->Branch("shiftK0sResPosLow1pt",&shiftK0sResPosLow1pt);
  brshiftK0sResPosLow1pt->Fill();
  TBranch *breShiftK0sResNegLow1pt = OutTree->Branch("eShiftK0sResNegLow1pt",&eShiftK0sResNegLow1pt);
  breShiftK0sResNegLow1pt->Fill();
  TBranch *breShiftK0sResPosLow1pt = OutTree->Branch("eShiftK0sResPosLow1pt",&eShiftK0sResPosLow1pt);
  breShiftK0sResPosLow1pt->Fill();
  TBranch *brsigmaK0sResNegHigh1pt = OutTree->Branch("sigmaK0sResNegHigh1pt",&sigmaK0sResNegHigh1pt);
  brsigmaK0sResNegHigh1pt->Fill();
  TBranch *brsigmaK0sResPosHigh1pt = OutTree->Branch("sigmaK0sResPosHigh1pt",&sigmaK0sResPosHigh1pt);
  brsigmaK0sResPosHigh1pt->Fill();
  TBranch *breSigmaK0sResNegHigh1pt = OutTree->Branch("eSigmaK0sResNegHigh1pt",&eSigmaK0sResNegHigh1pt);
  breSigmaK0sResNegHigh1pt->Fill();
  TBranch *breSigmaK0sResPosHigh1pt = OutTree->Branch("eSigmaK0sResPosHigh1pt",&eSigmaK0sResPosHigh1pt);
  breSigmaK0sResPosHigh1pt->Fill();
  TBranch *brsigmaK0sResNegLow1pt = OutTree->Branch("sigmaK0sResNegLow1pt",&sigmaK0sResNegLow1pt);
  brsigmaK0sResNegLow1pt->Fill();
  TBranch *brsigmaK0sResPosLow1pt = OutTree->Branch("sigmaK0sResPosLow1pt",&sigmaK0sResPosLow1pt);
  brsigmaK0sResPosLow1pt->Fill();
  TBranch *breSigmaK0sResNegLow1pt = OutTree->Branch("eSigmaK0sResNegLow1pt",&eSigmaK0sResNegLow1pt);
  breSigmaK0sResNegLow1pt->Fill();
  TBranch *breSigmaK0sResPosLow1pt = OutTree->Branch("eSigmaK0sResPosLow1pt",&eSigmaK0sResPosLow1pt);
  breSigmaK0sResPosLow1pt->Fill();

}

void AliHighPtTreeAnalysis::MakedcaRTrends(){

  Double_t  dcaRresCombinedLow1pt  = -999.;   Double_t  dcaRresCombinedHigh1pt  = -999.;
  Double_t edcaRresCombinedLow1pt  = -999.;   Double_t edcaRresCombinedHigh1pt  = -999.;

  Double_t  dcaRresTPCAsideLow1pt  = -999.;   Double_t  dcaRresTPCAsideHigh1pt  = -999.;
  Double_t edcaRresTPCAsideLow1pt  = -999.;   Double_t edcaRresTPCAsideHigh1pt  = -999.;

  Double_t  dcaRresTPCCsideLow1pt  = -999.;   Double_t  dcaRresTPCCsideHigh1pt  = -999.;
  Double_t edcaRresTPCCsideLow1pt  = -999.;   Double_t edcaRresTPCCsideHigh1pt  = -999.;


  Double_t  dcaRpullCombinedLow1pt = -999.;   Double_t  dcaRpullCombinedHigh1pt = -999.;
  Double_t edcaRpullCombinedLow1pt = -999.;   Double_t edcaRpullCombinedHigh1pt = -999.;

  Double_t  dcaRpullTPCAsideLow1pt = -999.;   Double_t  dcaRpullTPCAsideHigh1pt = -999.;
  Double_t edcaRpullTPCAsideLow1pt = -999.;   Double_t edcaRpullTPCAsideHigh1pt = -999.;

  Double_t  dcaRpullTPCCsideLow1pt = -999.;   Double_t  dcaRpullTPCCsideHigh1pt = -999.;
  Double_t edcaRpullTPCCsideLow1pt = -999.;   Double_t edcaRpullTPCCsideHigh1pt = -999.;

  Double_t low1pt  = 0.;
  Double_t high1pt = 1.;

  // Fit models
  //
  TF1 * fRes = new TF1("fRes","sqrt(([0]**2+[1]**2*x^abs([2])))");
  fRes->SetParName(0,"Offset");
  fRes->SetParName(1,"MS  slope");
  fRes->SetParName(2,"MS  exponent");
  TF1 * fPull = new TF1("fPull","[0]*sqrt((1+[1]**2*x^abs([4]))/([2]**2+[3]**2*x^abs([4])))");



  TH3D *hResdcaRTPConly_vs_eta_1pT_cl  = (TH3D*) hResdcaRTPConly_vs_eta_1pT ->Clone( Form("%s_cl",hResdcaRTPConly_vs_eta_1pT->GetName()) );
  TH3D *hPulldcaRTPConly_vs_eta_1pT_cl = (TH3D*) hPulldcaRTPConly_vs_eta_1pT->Clone( Form("%s_cl",hPulldcaRTPConly_vs_eta_1pT->GetName()) );

  TH2D *hResdcaRcombined_1pt       = (TH2D*) hResdcaRcomb_vs_eta_1pT->Project3D("zx");
  hResdcaRTPConly_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(0.,0.8);
  TH2D *hResdcaRTPConly_1pt_Aside  = (TH2D*) hResdcaRTPConly_vs_eta_1pT_cl->Project3D("zx");
  hResdcaRTPConly_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(-0.8,0.);
  TH2D *hResdcaRTPConly_1pt_Cside  = (TH2D*) hResdcaRTPConly_vs_eta_1pT_cl->Project3D("zx");

  TH2D *hPulldcaRcombined_1pt      = (TH2D*) hPulldcaRcomb_vs_eta_1pT->Project3D("zx");
  hPulldcaRTPConly_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(0.,0.8);
  TH2D *hPulldcaRTPConly_1pt_Aside = (TH2D*) hPulldcaRTPConly_vs_eta_1pT_cl->Project3D("zx");
  hPulldcaRTPConly_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(-0.8,0.);
  TH2D *hPulldcaRTPConly_1pt_Cside = (TH2D*) hPulldcaRTPConly_vs_eta_1pT_cl->Project3D("zx");

  TObjArray *arr1 = new TObjArray();

  TH1D *hdcaRresCombined  = 0;
  TH1D *hdcaRresTPCAside  = 0;
  TH1D *hdcaRresTPCCside  = 0;
  TH1D *hdcaRpullCombined = 0;
  TH1D *hdcaRpullTPCAside = 0;
  TH1D *hdcaRpullTPCCside = 0;

  hResdcaRcombined_1pt->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRresCombined = (TH1D*) arr1->At(2)->Clone("hdcaRresCombined");
  arr1->Delete();

  hResdcaRTPConly_1pt_Aside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRresTPCAside = (TH1D*) arr1->At(2)->Clone("hdcaRresTPCAside");
  arr1->Delete();

  hResdcaRTPConly_1pt_Cside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRresTPCCside = (TH1D*) arr1->At(2)->Clone("hdcaRresTPCCside");
  arr1->Delete();

  hPulldcaRcombined_1pt->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRpullCombined = (TH1D*) arr1->At(2)->Clone("hdcaRpullCombined");
  arr1->Delete();

  hPulldcaRTPConly_1pt_Aside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRpullTPCAside = (TH1D*) arr1->At(2)->Clone("hdcaRpullTPCAside");
  arr1->Delete(); 

  hPulldcaRTPConly_1pt_Cside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hdcaRpullTPCCside = (TH1D*) arr1->At(2)->Clone("hdcaRpullTPCCside");
  arr1->Delete();

  if( hdcaRresCombined->GetEntries() > 3 && hdcaRpullCombined->GetEntries() > 4 ){
    //Fit
    fRes->SetParameters(1,0.2,2);
    hdcaRresCombined->Fit(fRes,"Q");
    fPull->SetParameters(1,0.1,1,0.2,2);
    fPull->FixParameter(4,fRes->GetParameter(2));
    hdcaRpullCombined->Fit(fPull,"Qw","w");
    //Store Trending values
    dcaRresCombinedLow1pt  = fRes->Eval(low1pt);
    edcaRresCombinedLow1pt  = sqrt( pow(fRes->GetParError(0),2) + pow(low1pt*fRes->GetParError(1),2) );
    dcaRresCombinedHigh1pt = fRes->Eval(high1pt);
    edcaRresCombinedHigh1pt = sqrt( pow(fRes->GetParError(0),2) + pow(high1pt*fRes->GetParError(1),2) );

    dcaRpullCombinedLow1pt  = fPull->Eval(low1pt);
    edcaRpullCombinedLow1pt  = sqrt( pow(fPull->GetParError(0),2) + pow(low1pt*fPull->GetParError(1),2) );
    dcaRpullCombinedHigh1pt = fPull->Eval(high1pt);
    edcaRpullCombinedHigh1pt = sqrt( pow(fPull->GetParError(0),2) + pow(high1pt*fPull->GetParError(1),2) );
    //check fit performance
    if(fMakeFitPerfomancePlots){ TCanvas *c_tmp = new TCanvas("can","can",1100,550); c_tmp->Divide(2,1); c_tmp->cd(1); hdcaRresCombined->GetYaxis()->SetTitle("resol (unit)"); 
      //hdcaRresCombined->GetYaxis()->SetRangeUser(0.,0.01); 
      hdcaRresCombined->GetXaxis()->SetTitle("1/p_{t} (1/GeV)"); hdcaRresCombined->SetMarkerStyle(21); hdcaRresCombined->Draw(); c_tmp->cd(2);
      hdcaRpullCombined->GetYaxis()->SetTitle("pull (unit)"); //hdcaRpullCombined->GetYaxis()->SetRangeUser(0.,2.); 
      hdcaRpullCombined->GetXaxis()->SetTitle("1/p_{t} (1/GeV)");
      hdcaRpullCombined->SetMarkerStyle(21); hdcaRpullCombined->Draw();
      c_tmp->SaveAs("fit_dcaRcombined.png"); delete c_tmp; }
  }

  if( hdcaRresTPCAside->GetEntries() > 3 && hdcaRpullTPCAside->GetEntries() > 4 ){
    //Fit
    fRes->SetParameters(1,0.2,2);
    hdcaRresTPCAside->Fit(fRes,"Q");
    fPull->SetParameters(1,0.1,1,0.2,2);
    fPull->FixParameter(4,fRes->GetParameter(2));
    hdcaRpullTPCAside->Fit(fPull,"Qw","w");
    //Store Trending values
    dcaRresTPCAsideLow1pt  = fRes->Eval(low1pt);
    edcaRresTPCAsideLow1pt  = sqrt( pow(fRes->GetParError(0),2) + pow(low1pt*fRes->GetParError(1),2) );
    dcaRresTPCAsideHigh1pt = fRes->Eval(high1pt);
    edcaRresTPCAsideHigh1pt = sqrt( pow(fRes->GetParError(0),2) + pow(high1pt*fRes->GetParError(1),2) );

    dcaRpullTPCAsideLow1pt  = fPull->Eval(low1pt);
    edcaRpullTPCAsideLow1pt  = sqrt( pow(fPull->GetParError(0),2) + pow(low1pt*fPull->GetParError(1),2) );
    dcaRpullTPCAsideHigh1pt = fPull->Eval(high1pt);
    edcaRpullTPCAsideHigh1pt = sqrt( pow(fPull->GetParError(0),2) + pow(high1pt*fPull->GetParError(1),2) );
    //check fit performance
    if(fMakeFitPerfomancePlots){ TCanvas *c_tmp = new TCanvas("can","can",1100,550); c_tmp->Divide(2,1); c_tmp->cd(1); hdcaRresTPCAside->GetYaxis()->SetTitle("resol (unit)"); 
      //hdcaRresTPCAside->GetYaxis()->SetRangeUser(0.,0.35); 
      hdcaRresTPCAside->GetXaxis()->SetTitle("1/p_{t} (1/GeV)"); hdcaRresTPCAside->SetMarkerStyle(21); hdcaRresTPCAside->Draw(); c_tmp->cd(2);
      hdcaRpullTPCAside->GetYaxis()->SetTitle("pull (unit)");//hdcaRpullTPCAside->GetYaxis()->SetRangeUser(0.,3.); 
      hdcaRpullTPCAside->GetXaxis()->SetTitle("1/p_{t} (1/GeV)");
      hdcaRpullTPCAside->SetMarkerStyle(21); hdcaRpullTPCAside->Draw();
      c_tmp->SaveAs("fit_dcaRTPCAside.png"); delete c_tmp; }
  }

  if( hdcaRresTPCCside->GetEntries() > 3 && hdcaRpullTPCCside->GetEntries() > 4 ){
    //Fit
    fRes->SetParameters(1,0.2,2);
    hdcaRresTPCCside->Fit(fRes,"Q");
    fPull->SetParameters(1,0.1,1,0.2,2);
    fPull->FixParameter(4,fRes->GetParameter(2));
    hdcaRpullTPCCside->Fit(fPull,"Qw","w");
    //Store Trending values
    dcaRresTPCCsideLow1pt  = fRes->Eval(low1pt);
    edcaRresTPCCsideLow1pt  = sqrt( pow(fRes->GetParError(0),2) + pow(low1pt*fRes->GetParError(1),2) );
    dcaRresTPCCsideHigh1pt = fRes->Eval(high1pt);
    edcaRresTPCCsideHigh1pt = sqrt( pow(fRes->GetParError(0),2) + pow(high1pt*fRes->GetParError(1),2) );

    dcaRpullTPCCsideLow1pt  = fPull->Eval(low1pt);
    edcaRpullTPCCsideLow1pt  = sqrt( pow(fPull->GetParError(0),2) + pow(low1pt*fPull->GetParError(1),2) );
    dcaRpullTPCCsideHigh1pt = fPull->Eval(high1pt);
    edcaRpullTPCCsideHigh1pt = sqrt( pow(fPull->GetParError(0),2) + pow(high1pt*fPull->GetParError(1),2) );
    //check fit performance
    if(fMakeFitPerfomancePlots){ TCanvas *c_tmp = new TCanvas("can","can",1100,550); c_tmp->Divide(2,1); c_tmp->cd(1); hdcaRresTPCCside->GetYaxis()->SetTitle("resol (unit)"); 
      //hdcaRresTPCCside->GetYaxis()->SetRangeUser(0.,0.35);  
      hdcaRresTPCCside->GetXaxis()->SetTitle("1/p_{t} (1/GeV)"); hdcaRresTPCCside->SetMarkerStyle(21); hdcaRresTPCCside->Draw(); c_tmp->cd(2);
      hdcaRpullTPCCside->GetYaxis()->SetTitle("pull (unit)");// hdcaRpullTPCCside->GetYaxis()->SetRangeUser(0.,3.); 
      hdcaRpullTPCCside->GetXaxis()->SetTitle("1/p_{t} (1/GeV)");
      hdcaRpullTPCCside->SetMarkerStyle(21); hdcaRpullTPCCside->Draw();
      c_tmp->SaveAs("fit_dcaRTPCCside.png"); delete c_tmp; }
  }




  TBranch *brdcaRresCombinedLow1pt = OutTree->Branch("dcaRresCombinedLow1pt",&dcaRresCombinedLow1pt);
  brdcaRresCombinedLow1pt->Fill();
  TBranch *bredcaRresCombinedLow1pt = OutTree->Branch("edcaRresCombinedLow1pt",&edcaRresCombinedLow1pt);
  bredcaRresCombinedLow1pt->Fill();
  TBranch *brdcaRresCombinedHigh1pt = OutTree->Branch("dcaRresCombinedHigh1pt",&dcaRresCombinedHigh1pt);
  brdcaRresCombinedHigh1pt->Fill();
  TBranch *bredcaRresCombinedHigh1pt = OutTree->Branch("edcaRresCombinedHigh1pt",&edcaRresCombinedHigh1pt);
  bredcaRresCombinedHigh1pt->Fill();
  TBranch *brdcaRresTPCAsideLow1pt = OutTree->Branch("dcaRresTPCAsideLow1pt",&dcaRresTPCAsideLow1pt);
  brdcaRresTPCAsideLow1pt->Fill();
  TBranch *bredcaRresTPCAsideLow1pt = OutTree->Branch("edcaRresTPCAsideLow1pt",&edcaRresTPCAsideLow1pt);
  bredcaRresTPCAsideLow1pt->Fill();
  TBranch *brdcaRresTPCAsideHigh1pt = OutTree->Branch("dcaRresTPCAsideHigh1pt",&dcaRresTPCAsideHigh1pt);
  brdcaRresTPCAsideHigh1pt->Fill();
  TBranch *bredcaRresTPCAsideHigh1pt = OutTree->Branch("edcaRresTPCAsideHigh1pt",&edcaRresTPCAsideHigh1pt);
  bredcaRresTPCAsideHigh1pt->Fill();
  TBranch *brdcaRresTPCCsideLow1pt = OutTree->Branch("dcaRresTPCCsideLow1pt",&dcaRresTPCCsideLow1pt);
  brdcaRresTPCCsideLow1pt->Fill();
  TBranch *bredcaRresTPCCsideLow1pt = OutTree->Branch("edcaRresTPCCsideLow1pt",&edcaRresTPCCsideLow1pt);
  bredcaRresTPCCsideLow1pt->Fill();
  TBranch *brdcaRresTPCCsideHigh1pt = OutTree->Branch("dcaRresTPCCsideHigh1pt",&dcaRresTPCCsideHigh1pt);
  brdcaRresTPCCsideHigh1pt->Fill();
  TBranch *bredcaRresTPCCsideHigh1pt = OutTree->Branch("edcaRresTPCCsideHigh1pt",&edcaRresTPCCsideHigh1pt);
  bredcaRresTPCCsideHigh1pt->Fill();
  TBranch *brdcaRpullCombinedLow1pt = OutTree->Branch("dcaRpullCombinedLow1pt",&dcaRpullCombinedLow1pt);
  brdcaRpullCombinedLow1pt->Fill();
  TBranch *bredcaRpullCombinedLow1pt = OutTree->Branch("edcaRpullCombinedLow1pt",&edcaRpullCombinedLow1pt);
  bredcaRpullCombinedLow1pt->Fill();
  TBranch *brdcaRpullCombinedHigh1pt = OutTree->Branch("dcaRpullCombinedHigh1pt",&dcaRpullCombinedHigh1pt);
  brdcaRpullCombinedHigh1pt->Fill();
  TBranch *bredcaRpullCombinedHigh1pt = OutTree->Branch("edcaRpullCombinedHigh1pt",&edcaRpullCombinedHigh1pt);
  bredcaRpullCombinedHigh1pt->Fill();
  TBranch *brdcaRpullTPCAsideLow1pt = OutTree->Branch("dcaRpullTPCAsideLow1pt",&dcaRpullTPCAsideLow1pt);
  brdcaRpullTPCAsideLow1pt->Fill();
  TBranch *bredcaRpullTPCAsideLow1pt = OutTree->Branch("edcaRpullTPCAsideLow1pt",&edcaRpullTPCAsideLow1pt);
  bredcaRpullTPCAsideLow1pt->Fill();
  TBranch *brdcaRpullTPCAsideHigh1pt = OutTree->Branch("dcaRpullTPCAsideHigh1pt",&dcaRpullTPCAsideHigh1pt);
  brdcaRpullTPCAsideHigh1pt->Fill();
  TBranch *bredcaRpullTPCAsideHigh1pt = OutTree->Branch("edcaRpullTPCAsideHigh1pt",&edcaRpullTPCAsideHigh1pt);
  bredcaRpullTPCAsideHigh1pt->Fill();
  TBranch *brdcaRpullTPCCsideLow1pt = OutTree->Branch("dcaRpullTPCCsideLow1pt",&dcaRpullTPCCsideLow1pt);
  brdcaRpullTPCCsideLow1pt->Fill();
  TBranch *bredcaRpullTPCCsideLow1pt = OutTree->Branch("edcaRpullTPCCsideLow1pt",&edcaRpullTPCCsideLow1pt);
  bredcaRpullTPCCsideLow1pt->Fill();
  TBranch *brdcaRpullTPCCsideHigh1pt = OutTree->Branch("dcaRpullTPCCsideHigh1pt",&dcaRpullTPCCsideHigh1pt);
  brdcaRpullTPCCsideHigh1pt->Fill();
  TBranch *bredcaRpullTPCCsideHigh1pt = OutTree->Branch("edcaRpullTPCCsideHigh1pt",&edcaRpullTPCCsideHigh1pt);
  bredcaRpullTPCCsideHigh1pt->Fill();



}

void AliHighPtTreeAnalysis::MakeDeltaPhiTrends(){

  Double_t  dPhiResTPCAsideLow1pt = -999.;   Double_t  dPhiResTPCAsideHigh1pt = -999.;
  Double_t edPhiResTPCAsideLow1pt = -999.;   Double_t edPhiResTPCAsideHigh1pt = -999.;

  Double_t  dPhiResTPCCsideLow1pt = -999.;   Double_t  dPhiResTPCCsideHigh1pt = -999.;
  Double_t edPhiResTPCCsideLow1pt = -999.;   Double_t edPhiResTPCCsideHigh1pt = -999.;  

  Double_t  dPhiPullTPCAsideLow1pt = -999.;   Double_t  dPhiPullTPCAsideHigh1pt = -999.;
  Double_t edPhiPullTPCAsideLow1pt = -999.;   Double_t edPhiPullTPCAsideHigh1pt = -999.;

  Double_t  dPhiPullTPCCsideLow1pt = -999.;   Double_t  dPhiPullTPCCsideHigh1pt = -999.;
  Double_t edPhiPullTPCCsideLow1pt = -999.;   Double_t edPhiPullTPCCsideHigh1pt = -999.;


  Double_t low1pt  = 0.;
  Double_t high1pt = 1.;
  TF1 *fLinear = new TF1("fLinear","[0] + [1]*x",0.,1.);

  TH3D *hphiRes_vs_eta_1pT_cl  = (TH3D*) hphiRes_vs_eta_1pT ->Clone(Form( "%s_cl",hphiRes_vs_eta_1pT ->GetName() ));   
  TH3D *hphiPull_vs_eta_1pT_cl = (TH3D*) hphiPull_vs_eta_1pT->Clone(Form( "%s_cl",hphiPull_vs_eta_1pT->GetName() ));

  hphiRes_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(0.,0.8);
  TH2D *hphiRes_1pt_Aside = (TH2D*) hphiRes_vs_eta_1pT_cl->Project3D("zx");
  hphiRes_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(-0.8,0.);
  TH2D *hphiRes_1pt_Cside = (TH2D*) hphiRes_vs_eta_1pT_cl->Project3D("zx");

  hphiPull_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(0.,0.8);
  TH2D *hphiPull_1pt_Aside = (TH2D*) hphiPull_vs_eta_1pT_cl->Project3D("zx");
  hphiPull_vs_eta_1pT_cl->GetYaxis()->SetRangeUser(-0.8,0.);
  TH2D *hphiPull_1pt_Cside = (TH2D*) hphiPull_vs_eta_1pT_cl->Project3D("zx");

  TObjArray *arr1 = new TObjArray();

  TH1D *hphiRes_Aside  = 0;
  TH1D *hphiRes_Cside  = 0;
  TH1D *hphiPull_Aside = 0;
  TH1D *hphiPull_Cside = 0;

  hphiRes_1pt_Aside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hphiRes_Aside  = (TH1D*) arr1->At(2)->Clone("hphiRes_Aside");
  arr1->Delete();

  hphiRes_1pt_Cside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hphiRes_Cside  = (TH1D*) arr1->At(2)->Clone("hphiRes_Cside");
  arr1->Delete();

  hphiPull_1pt_Aside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hphiPull_Aside = (TH1D*) arr1->At(2)->Clone("hphiPull_Aside");
  arr1->Delete();

  hphiPull_1pt_Cside->FitSlicesY(0,0,-1,0,"QNR",arr1);
  if(arr1->At(2)) hphiPull_Cside = (TH1D*) arr1->At(2)->Clone("hphiPull_Cside");
  arr1->Delete();

  if( hphiRes_Aside->GetEntries() > 2 ){
    hphiRes_Aside->Fit(fLinear,"Q");

    dPhiResTPCAsideLow1pt  = fLinear->Eval(low1pt);
    edPhiResTPCAsideLow1pt  = sqrt( pow(fLinear->GetParError(0),2) + pow(low1pt*fLinear->GetParError(1),2) );
    dPhiResTPCAsideHigh1pt = fLinear->Eval(high1pt);
    edPhiResTPCAsideHigh1pt = sqrt( pow(fLinear->GetParError(0),2) + pow(high1pt*fLinear->GetParError(1),2) );

    fLinear->SetParameters(0,1);
  } 
  if( hphiRes_Cside->GetEntries() > 2 ){
    hphiRes_Cside->Fit(fLinear,"Q");

    dPhiResTPCCsideLow1pt  = fLinear->Eval(low1pt);
    edPhiResTPCCsideLow1pt  = sqrt( pow(fLinear->GetParError(0),2) + pow(low1pt*fLinear->GetParError(1),2) );
    dPhiResTPCCsideHigh1pt = fLinear->Eval(high1pt);
    edPhiResTPCCsideHigh1pt = sqrt( pow(fLinear->GetParError(0),2) + pow(high1pt*fLinear->GetParError(1),2) );

    fLinear->SetParameters(0,1);
  }
  if( hphiPull_Aside->GetEntries() > 2 ){
    hphiPull_Aside->Fit(fLinear,"Q");

    dPhiPullTPCAsideLow1pt  = fLinear->Eval(low1pt);
    edPhiPullTPCAsideLow1pt  = sqrt( pow(fLinear->GetParError(0),2) + pow(low1pt*fLinear->GetParError(1),2) );
    dPhiPullTPCAsideHigh1pt = fLinear->Eval(high1pt);
    edPhiPullTPCAsideHigh1pt = sqrt( pow(fLinear->GetParError(0),2) + pow(high1pt*fLinear->GetParError(1),2) );

    fLinear->SetParameters(0,1);
  }
  if( hphiPull_Cside->GetEntries() > 2 ){
    hphiPull_Cside->Fit(fLinear,"Q");

    dPhiPullTPCCsideLow1pt  = fLinear->Eval(low1pt);
    edPhiPullTPCCsideLow1pt  = sqrt( pow(fLinear->GetParError(0),2) + pow(low1pt*fLinear->GetParError(1),2) );
    dPhiPullTPCCsideHigh1pt = fLinear->Eval(high1pt);
    edPhiPullTPCCsideHigh1pt = sqrt( pow(fLinear->GetParError(0),2) + pow(high1pt*fLinear->GetParError(1),2) );

    fLinear->SetParameters(0,1);
  }  

  TBranch *brdPhiResTPCAsideLow1pt = OutTree->Branch("dPhiResTPCAsideLow1pt",&dPhiResTPCAsideLow1pt);
  brdPhiResTPCAsideLow1pt->Fill();
  TBranch *bredPhiResTPCAsideLow1pt = OutTree->Branch("edPhiResTPCAsideLow1pt",&edPhiResTPCAsideLow1pt);
  bredPhiResTPCAsideLow1pt->Fill();
  TBranch *brdPhiResTPCAsideHigh1pt = OutTree->Branch("dPhiResTPCAsideHigh1pt",&dPhiResTPCAsideHigh1pt);
  brdPhiResTPCAsideHigh1pt->Fill();
  TBranch *bredPhiResTPCAsideHigh1pt = OutTree->Branch("edPhiResTPCAsideHigh1pt",&edPhiResTPCAsideHigh1pt);
  bredPhiResTPCAsideHigh1pt->Fill();
  TBranch *brdPhiResTPCCsideLow1pt = OutTree->Branch("dPhiResTPCCsideLow1pt",&dPhiResTPCCsideLow1pt);
  brdPhiResTPCCsideLow1pt->Fill();
  TBranch *bredPhiResTPCCsideLow1pt = OutTree->Branch("edPhiResTPCCsideLow1pt",&edPhiResTPCCsideLow1pt);
  bredPhiResTPCCsideLow1pt->Fill();
  TBranch *brdPhiResTPCCsideHigh1pt = OutTree->Branch("dPhiResTPCCsideHigh1pt",&dPhiResTPCCsideHigh1pt);
  brdPhiResTPCCsideHigh1pt->Fill();
  TBranch *bredPhiResTPCCsideHigh1pt = OutTree->Branch("edPhiResTPCCsideHigh1pt",&edPhiResTPCCsideHigh1pt);
  bredPhiResTPCCsideHigh1pt->Fill();
  TBranch *brdPhiPullTPCAsideLow1pt = OutTree->Branch("dPhiPullTPCAsideLow1pt",&dPhiPullTPCAsideLow1pt);
  brdPhiPullTPCAsideLow1pt->Fill();
  TBranch *bredPhiPullTPCAsideLow1pt = OutTree->Branch("edPhiPullTPCAsideLow1pt",&edPhiPullTPCAsideLow1pt);
  bredPhiPullTPCAsideLow1pt->Fill();
  TBranch *brdPhiPullTPCAsideHigh1pt = OutTree->Branch("dPhiPullTPCAsideHigh1pt",&dPhiPullTPCAsideHigh1pt);
  brdPhiPullTPCAsideHigh1pt->Fill();
  TBranch *bredPhiPullTPCAsideHigh1pt = OutTree->Branch("edPhiPullTPCAsideHigh1pt",&edPhiPullTPCAsideHigh1pt);
  bredPhiPullTPCAsideHigh1pt->Fill();
  TBranch *brdPhiPullTPCCsideLow1pt = OutTree->Branch("dPhiPullTPCCsideLow1pt",&dPhiPullTPCCsideLow1pt);
  brdPhiPullTPCCsideLow1pt->Fill();
  TBranch *bredPhiPullTPCCsideLow1pt = OutTree->Branch("edPhiPullTPCCsideLow1pt",&edPhiPullTPCCsideLow1pt);
  bredPhiPullTPCCsideLow1pt->Fill();
  TBranch *brdPhiPullTPCCsideHigh1pt = OutTree->Branch("dPhiPullTPCCsideHigh1pt",&dPhiPullTPCCsideHigh1pt);
  brdPhiPullTPCCsideHigh1pt->Fill();
  TBranch *bredPhiPullTPCCsideHigh1pt = OutTree->Branch("edPhiPullTPCCsideHigh1pt",&edPhiPullTPCCsideHigh1pt);
  bredPhiPullTPCCsideHigh1pt->Fill();

}

void AliHighPtTreeAnalysis::MakeEfficiencyTrends(){

  Double_t  EfficiencyLowPt    = -999.;  Double_t  EfficiencyHighPt    = -999.;
  Double_t eEfficiencyLowPt    = -999.;  Double_t eEfficiencyHighPt    = -999.;

  if(fNtracks_TPCLowPt > 0 && fNtracks_TPCITSLowPt > 0){
    EfficiencyLowPt  = (Double_t) fNtracks_TPCITSLowPt / (Double_t) fNtracks_TPCLowPt;
    eEfficiencyLowPt  = sqrt( pow(sqrt(fNtracks_TPCITSLowPt)/fNtracks_TPCLowPt,2) + pow(sqrt(fNtracks_TPCLowPt)*fNtracks_TPCITSLowPt/(fNtracks_TPCLowPt*fNtracks_TPCLowPt),2) );
  }
  if(fNtracks_TPCHighPt > 0 && fNtracks_TPCITSHighPt > 0){
    EfficiencyHighPt = (Double_t) fNtracks_TPCITSHighPt / (Double_t) fNtracks_TPCHighPt;
    eEfficiencyHighPt = sqrt( pow(sqrt(fNtracks_TPCITSHighPt)/fNtracks_TPCHighPt,2) + pow(sqrt(fNtracks_TPCHighPt)*fNtracks_TPCITSHighPt/(fNtracks_TPCHighPt*fNtracks_TPCHighPt),2) );
  }

  TBranch *brEfficiencyLowPt = OutTree->Branch("EfficiencyLowPt",&EfficiencyLowPt);
  brEfficiencyLowPt->Fill();
  TBranch *breEfficiencyLowPt = OutTree->Branch("eEfficiencyLowPt",&eEfficiencyLowPt);
  breEfficiencyLowPt->Fill();
  TBranch *brEfficiencyHighPt = OutTree->Branch("EfficiencyHighPt",&EfficiencyHighPt);
  brEfficiencyHighPt->Fill();
  TBranch *breEfficiencyHighPt = OutTree->Branch("eEfficiencyHighPt",&eEfficiencyHighPt);
  breEfficiencyHighPt->Fill();

}

void AliHighPtTreeAnalysis::PlotEff(TH3D *hTPCITS, TH3D *hTPC, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName ,const char *plotName){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);

  if(!hTPCITS || !hTPC) return;
  char name[256];
  sprintf(name,"%s_1",hTPCITS->GetName());
  TH3D *h1c = (TH3D*)hTPCITS->Clone(name);
  sprintf(name,"%s_1",hTPC->GetName());
  TH3D *h2c = (TH3D*)hTPC->Clone(name);
  if(!h1c || !h2c) return;
  if(fPtCut>0.) h1c->GetZaxis()->SetRangeUser(fPtCut,100);
  if(fPtCut>0.) h2c->GetZaxis()->SetRangeUser(fPtCut,100);
  TH1D *h1proj = (TH1D*)h1c->Project3D(projAxisName);
  TH1D *h2proj = (TH1D*)h2c->Project3D(projAxisName);
  if(!h1proj || !h2proj) return;

  h1proj->Sumw2();
  h1proj->Divide(h2proj);
  h1proj->SetMarkerColor(1);
  h1proj->SetLineColor(1);
  h1proj->SetMarkerStyle(20);
  h1proj->GetXaxis()->SetTitle(xaxisName);
  h1proj->GetXaxis()->SetTitleOffset(1.4);
  h1proj->GetYaxis()->SetTitleOffset(1.4);
  h1proj->GetXaxis()->SetLabelSize(0.035);
  h1proj->GetYaxis()->SetLabelSize(0.035);
  h1proj->GetXaxis()->SetLabelFont(62);
  h1proj->GetYaxis()->SetLabelFont(62);
  h1proj->SetTitle(plotName);

  TCanvas *can = new TCanvas("can","can",550,550);
  can->cd(1);
  if(logX) gPad->SetLogx();
  if(histoType==0) gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1proj->Draw("e1hist");

  sprintf(name,"%s.png",plotName);
  can->SaveAs(name);
  //    sprintf(name,"%s.eps",plotName);
  //    can->SaveAs(name);
  //    sprintf(name,"%s.pdf",plotName);
  //    can->SaveAs(name);

  if(fPtCut>0.) h1c->GetZaxis()->SetRangeUser(0.,100);
  if(fPtCut>0.) h2c->GetZaxis()->SetRangeUser(0.,100);

  delete can;
  delete h1proj;
  delete h2proj;
}

void AliHighPtTreeAnalysis::Plot1D(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  // Make projections

  if(!h1) return;

  char name[256];
  sprintf(name,"%s_1",h1->GetName());
  TH3D *h1c = (TH3D*)h1->Clone(name);
  if(!h1c) return;

  if(fPtCut>0.) h1c->GetZaxis()->SetRangeUser(fPtCut,100);

  TH1D *h1proj = (TH1D*)h1c->Project3D(projAxisName);
  if(!h1proj) return;

  // Set properties

  h1proj->SetMarkerColor(1);
  h1proj->SetLineColor(1);
  h1proj->SetMarkerStyle(20);
  h1proj->GetXaxis()->SetTitle(xaxisName);
  h1proj->GetXaxis()->SetTitleOffset(1.4);
  h1proj->GetYaxis()->SetTitleOffset(1.4);
  h1proj->GetXaxis()->SetLabelSize(0.035);
  h1proj->GetYaxis()->SetLabelSize(0.035);
  h1proj->GetXaxis()->SetLabelFont(62);
  h1proj->GetYaxis()->SetLabelFont(62);
  h1proj->SetTitle(plotName);

  // Draw histo
  TCanvas *can = new TCanvas("can","can",550,550);
  can->cd(1);
  if(logX) gPad->SetLogx();
  if(histoType==0) gPad->SetLogy();
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1proj->Draw("e1hist");

  sprintf(name,"%s.png",plotName);
  can->SaveAs(name);
  //    sprintf(name,"%s.eps",plotName);
  //    can->SaveAs(name);
  //    sprintf(name,"%s.pdf",plotName);
  //    can->SaveAs(name);

  if(fPtCut>0.) h1c->GetZaxis()->SetRangeUser(0.,100);

  delete can;
  delete h1proj;
}

void AliHighPtTreeAnalysis::Plot2D(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  // Make projections

  if(!h1) return;

  char name[256];
  sprintf(name,"%s_1",h1->GetName());
  TH3D *h1c = (TH3D*)h1->Clone(name);
  if(!h1c) return;


  if(fPtCut>0.) h1c->GetXaxis()->SetRangeUser(fPtCut,100);

  TH2D *h1proj = (TH2D*)h1c->Project3D(projAxisName);
  if(!h1proj) return;

  // Fit slices
  TObjArray *arr1 = new TObjArray();
  h1proj->FitSlicesY(0,0,-1,0,"QNR",arr1);

  if(!arr1->At(1)) return;
  if(!arr1->At(2)) return;

  // Get histo
  TH1D *h1mean  = (TH1D*)arr1->At(1);
  TH1D *h1width = (TH1D*)arr1->At(2);


  // Set properties

  if (histoType==0) { // pulls dca
    SetHistoProperties(h1mean, 20, 1, -1.0, 1.0);
    SetHistoProperties(h1width, 20, 1, 0, 3);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean DCAr/#sigma(DCAr)");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma DCAr/#sigma(DCAr)");
  }

  if (histoType==1) { // resolution dca
    SetHistoProperties(h1mean, 20, 1, -0.01, 0.01);
    SetHistoProperties(h1width, 20, 1, 0, 0.03);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean DCAr (cm)");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma DCAr (cm)");
  }

  if (histoType==2) { // resolution dca tpc
    SetHistoProperties(h1mean, 20, 1, -0.2, 0.2);
    SetHistoProperties(h1width, 20, 1, 0, 1);
    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean DCAr (cm)");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma DCAr (cm)");
  }

  if (histoType==3) { // pulls pT
    SetHistoProperties(h1mean, 20, 1, -1.0, 1.0);
    SetHistoProperties(h1width, 20, 1, 0, 3);
    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (1/p_{T}-1/p_{T,MC})/#sigma(1/p_{T})");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (1/p_{T}-1/p_{T,MC})/#sigma(1/p_{T})");
  }


  if (histoType==4) { // resolution pt
    SetHistoProperties(h1mean, 20, 1, -0.02, 0.02);
    SetHistoProperties(h1width, 20, 1, 0, 0.2);

    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
  }

  if (histoType==5) { // resolution pt tpc
    SetHistoProperties(h1mean, 20, 1, -0.02, 0.02);
    SetHistoProperties(h1width, 20, 1, 0, 1.0);

    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
  }

  if (histoType==6) { // resolution pt tpc constrain
    SetHistoProperties(h1mean, 20, 1, -0.02, 0.02);
    SetHistoProperties(h1width, 20, 1, 0, 1.0);

    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
  }

  if (histoType==7) { // pulls phi
    SetHistoProperties(h1mean, 20, 1, -1.0, 1.0);
    if(logX == 1) SetHistoProperties(h1width, 20, 1, 0, 3);
    if(logX == 0) SetHistoProperties(h1width, 20, 1, 0, 1.5);
    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    //h1mean->GetYaxis()->SetTitle("mean (#phi_{TPC+ITS}-#phi_{TPCc})/#sqrt{#sigma(#phi_{TPC+ITS}))^2+(#sigma(#phi_{TPCc}))^2}");
    h1mean->GetYaxis()->SetTitle("mean (#phi_{TPC+ITS}-#phi_{TPCc})/#sigma(#Delta#phi)");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (#phi_{TPC+ITS}-#phi_{TPCc})/#sigma(#Delta#phi)");
  }

  if (histoType==8) { // resol phi
    SetHistoProperties(h1mean, 20, 1, -0.002, 0.002);
    if(logX == 1) SetHistoProperties(h1width, 20, 1, 0, 0.01);
    if(logX == 0) SetHistoProperties(h1width, 20, 1, 0, 0.002);
    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (#phi_{TPC+ITS}-#phi_{TPCc}) (rad)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (#phi_{TPC+ITS}-#phi_{TPCc}) (rad)");
  }


  if (histoType==9) { // resolution pt vs phi
    SetHistoProperties(h1mean, 20, 1, -0.02, 0.02);
    SetHistoProperties(h1width, 20, 1, 0, 0.1);

    //h1mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (1/p_{T}-1/p_{T,MC})/(1/p_{T,MC})");
  }

  if (histoType==10) { // pulls k0s vs 1/pT
    SetHistoProperties(h1mean, 20, 1, -1.0, 1.0);
    SetHistoProperties(h1width, 20, 1, 0, 3);
    //h1mean->GetXaxis()->SetTitle("1/p_{T} (GeV/c)");
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (M_K0s-M_K0s_PDG)/#sigma(M_K0s))");
    //h1width->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (M_K0s-M_K0s_PDG)/#sigma(M_K0s))");
  }

  h1mean ->GetYaxis()->SetTitleOffset(1.38);
  h1width->GetYaxis()->SetTitleOffset(1.38);

  h1mean->SetTitle(plotName);
  h1width->SetTitle(plotName);

  // Draw histo
  TCanvas *can = new TCanvas("can","can",550,700);
  can->Divide(1,2);

  can->cd(1);
  if(logX) gPad->SetLogx();
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1mean->Draw();

  can->cd(2);
  if(logX) gPad->SetLogx();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  h1width->Draw();

  sprintf(name,"%s.png",plotName);
  can->SaveAs(name);

  if(fPtCut>0.) h1c->GetXaxis()->SetRangeUser(0.,100.);

  delete can;
  delete h1proj;
  delete arr1;
}

void AliHighPtTreeAnalysis::Plot2DK0s(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName){
  if(!h1) return;
  gStyle->SetOptTitle(1);

  char name[256];
  sprintf(name,"%s_1",h1->GetName());
  TH3D *h1c = (TH3D*)h1->Clone(name);
  if(!h1c) return;

  if(fPtCut>0.) h1c->GetXaxis()->SetRangeUser(0,1./fPtCut);

  TH2D *h1proj = (TH2D*)h1c->Project3D(projAxisName);
  if(!h1proj) return;
  const char *OutXaxis  = "x";
  if(fPtCut>0.) OutXaxis = "y";
  TH1D *h1mean = (TH1D*)h1c->Project3D(OutXaxis);
  if(!h1mean) return;
  h1mean->SetName("h1mean");
  h1mean->Reset();
  TH1D *h1width = (TH1D*)h1c->Project3D(OutXaxis);
  if(!h1width) return;
  h1width->SetName("h1width");
  h1width->Reset();
  TF1 * fK0sFit = new TF1("fK0sFit","[0]+[1]*x+[2]*TMath::Gaus(x,[3],[4])",-0.5,0.5); 
  Int_t nBinsX = h1proj->GetXaxis()->GetNbins();
  for(Int_t i = 1; i<nBinsX; i++){
    h1proj->GetXaxis()->SetRange(i,i);
    if(h1proj->Integral()<50) continue;
    sprintf(name,"h1K0sproj_%s_%d",h1->GetName(),i);
    TH1D *h1K0sproj = (TH1D *)h1proj->ProjectionY(name);
    Double_t y1 = h1K0sproj->GetBinContent(5);
    Double_t y2 = h1K0sproj->GetBinContent(95);
    fK0sFit->SetParameter(0,(y1+y2)*0.5);
    fK0sFit->SetParameter(1,(y2-y1)/20.);
    fK0sFit->SetParameter(2,h1K0sproj->GetMaximum());
    fK0sFit->SetParameter(3,0);
    fK0sFit->SetParameter(4,h1K0sproj->GetRMS());
    h1K0sproj->Fit(fK0sFit,"Q"); 
    h1K0sproj->Fit(fK0sFit,"Q"); 

    h1mean->SetBinContent(i,fK0sFit->GetParameter(3));
    h1mean->SetBinError(i,fK0sFit->GetParError(3));
    h1width->SetBinContent(i,TMath::Abs(fK0sFit->GetParameter(4)));
    h1width->SetBinError(i,fK0sFit->GetParError(4));
  }

  if(histoType == 0){
    SetHistoProperties(h1mean, 20, 1, -1.0, 1.0);
    SetHistoProperties(h1width, 20, 1, 0, 3);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (M_K0s-M_K0s_PDG)/#sigma(M_K0s))");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (M_K0s-M_K0s_PDG)/#sigma(M_K0s))");
  }

  if(histoType == 1){
    SetHistoProperties(h1mean, 20, 1, -0.01, 0.01);
    SetHistoProperties(h1width, 20, 1, 0, 0.02);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("mean (M_K0s-M_K0s_PDG) (GeV/c^{2})");
    h1width->GetXaxis()->SetTitle(xaxisName);
    h1width->GetYaxis()->SetTitle("sigma (M_K0s-M_K0s_PDG) (GeV/c^{2})");
  }

  h1mean->GetXaxis()->SetLabelSize(0.05);
  h1width->GetXaxis()->SetLabelSize(0.05);
  h1mean->GetXaxis()->SetTitleSize(0.05);
  h1width->GetXaxis()->SetTitleSize(0.05);
  h1mean->SetTitle(plotName);
  h1width->SetTitle(plotName);

  h1mean ->GetYaxis()->SetTitleOffset(1.38);
  h1width->GetYaxis()->SetTitleOffset(1.38);

  // Draw histo
  TCanvas *can = new TCanvas("can","can",550,700);
  can->Divide(1,2);
  can->cd(1);
  if(logX) gPad->SetLogx();
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1mean->Draw("e");

  can->cd(2);
  if(logX) gPad->SetLogx();
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  h1width->Draw("e");

  sprintf(name,"%s.png",plotName);
  can->SaveAs(name);

  if(fPtCut>0.) h1c->GetXaxis()->SetRangeUser(0.,1.);

  delete can;
  delete h1proj;
  delete fK0sFit;
}

void AliHighPtTreeAnalysis::Plot1PtRes(TH3D *h1, const char *projAxisName, Int_t histoType, Int_t logX, const char *xaxisName, const char *plotName){

  // Make projections
  gStyle->SetOptTitle(1);

  if(!h1) return;

  char name[256];
  sprintf(name,"%s_1",h1->GetName());
  TH3D *h1c = (TH3D*)h1->Clone(name);
  if(!h1c) return;


  if(fPtCut>0.) {
    h1c->GetXaxis()->SetRangeUser(fPtCut,100);
  }

  TH2D *h1proj = (TH2D*)h1c->Project3D(projAxisName);
  if(!h1proj) return;

  // Fit slices
  TObjArray *arr1 = new TObjArray();
  h1proj->FitSlicesY(0,0,-1,0,"QNR",arr1);

  if(!arr1->At(1)) return;
  if(!arr1->At(2)) return;

  // Get histo
  TH1D *h1mean  = (TH1D*)arr1->At(1);
  //  TH1D *h1width = (TH1D*)arr1->At(2);

  // Set properties

  if (histoType==0) { // resolution pt vs pt
    SetHistoProperties(h1mean, 20, 1, 0, 0.2);
    h1mean->GetYaxis()->SetTitleOffset(2.0);
    h1mean->GetXaxis()->SetTitleSize(0.035);
    h1mean->GetYaxis()->SetTitleSize(0.035);
    h1mean->GetXaxis()->SetLabelSize(0.035);
    h1mean->GetYaxis()->SetLabelSize(0.035);
    h1mean->GetXaxis()->SetLabelFont(62);
    h1mean->GetYaxis()->SetLabelFont(62);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T}");
  }

  if (histoType==1) { // resolution pt vs phi
    SetHistoProperties(h1mean, 20, 1, 0, 0.05);
    h1mean->GetYaxis()->SetTitleOffset(2.0);
    h1mean->GetXaxis()->SetTitleSize(0.035);
    h1mean->GetYaxis()->SetTitleSize(0.035);
    h1mean->GetXaxis()->SetLabelSize(0.035);
    h1mean->GetYaxis()->SetLabelSize(0.035);
    h1mean->GetXaxis()->SetLabelFont(62);
    h1mean->GetYaxis()->SetLabelFont(62);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("#sigma(p_{T})/p_{T}");
  }

  if (histoType==2) { // sigma 1pt vs pt
    SetHistoProperties(h1mean, 20, 1, 0, 0.1);
    h1mean->GetYaxis()->SetTitleOffset(2.0);
    h1mean->GetXaxis()->SetTitleSize(0.035);
    h1mean->GetYaxis()->SetTitleSize(0.035);
    h1mean->GetXaxis()->SetLabelSize(0.035);
    h1mean->GetYaxis()->SetLabelSize(0.035);
    h1mean->GetXaxis()->SetLabelFont(62);
    h1mean->GetYaxis()->SetLabelFont(62);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("#sigma(p_{T}) (1/GeV/c)");
  }

  if (histoType==3) { // sigma 1pt vs phi
    SetHistoProperties(h1mean, 20, 1, 0, 0.007);
    h1mean->GetYaxis()->SetTitleOffset(2.0);
    h1mean->GetXaxis()->SetTitleSize(0.035);
    h1mean->GetYaxis()->SetTitleSize(0.035);
    h1mean->GetXaxis()->SetLabelSize(0.035);
    h1mean->GetYaxis()->SetLabelSize(0.035);
    h1mean->GetXaxis()->SetLabelFont(62);
    h1mean->GetYaxis()->SetLabelFont(62);
    h1mean->GetXaxis()->SetTitle(xaxisName);
    h1mean->GetYaxis()->SetTitle("#sigma(1/p_{T}) (1/GeV/c)");
  }


  h1mean->SetTitle(plotName);


  // Draw histo
  TCanvas *can = new TCanvas("can","can",550,550);
  can->cd();
  if(logX) gPad->SetLogx();
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1mean->Draw();

  sprintf(name,"%s.png",plotName);
  can->SaveAs(name);
  //    sprintf(name,"%s.eps",plotName);
  //    can->SaveAs(name);
  //    sprintf(name,"%s.pdf",plotName);
  //    can->SaveAs(name);

  //
  if(fPtCut>0.) {
    h1c->GetXaxis()->SetRangeUser(0.,100.);
  }

  if(can) delete can;

}

void AliHighPtTreeAnalysis::SetHistoProperties(TH1D *hist, Int_t marker, Int_t color, Double_t yMin, Double_t yMax){
  hist->SetMarkerStyle(marker);
  hist->SetMarkerColor(color);
  hist->SetLineColor(color);
  hist->GetYaxis()->SetRangeUser(yMin,yMax);
}

void AliHighPtTreeAnalysis::MakeAllPlots(){

  //vs Phi
  //Aside
  if(hPulldcaR_vs_phi_pT_Aside->GetEntries())          Plot2D(hPulldcaR_vs_phi_pT_Aside,             "zy",0,0,"#phi (rad)","dcaRPull_phi_HighPt_TPCAside");
  if(hPulldcaRTPCInner_vs_phi_pT_Aside->GetEntries())  Plot2D(hPulldcaRTPCInner_vs_phi_pT_Aside,     "zy",0,0,"#phi (rad)","dcaRPullTPCInner_phi_HighPt_TPCAside");
  if(hResdcaR_vs_phi_pT_Aside->GetEntries())           Plot2D(hResdcaR_vs_phi_pT_Aside,              "zy",1,0,"#phi (rad)","dcaRRes_phi_HighPt_TPCAside");
  if(hResdcaRTPCInner_vs_phi_pT_Aside->GetEntries())   Plot2D(hResdcaRTPCInner_vs_phi_pT_Aside,      "zy",2,0,"#phi (rad)","dcaRResTPCInner_phi_HighPt_TPCAside");
  if(hphiPull_vs_phi_pT_Aside->GetEntries())           Plot2D(hphiPull_vs_phi_pT_Aside,              "zy",7,0,"#phi (rad)","phiPullTPCInner_phi_HighPt_TPCAside");
  if(hphiRes_vs_phi_pT_Aside->GetEntries())            Plot2D(hphiRes_vs_phi_pT_Aside,               "zy",8,0,"#phi (rad)","phiResTPCInner_phi_HighPt_TPCAside");

  if(h1ptRes_vs_phi_pT_Aside->GetEntries())            Plot1PtRes(h1ptRes_vs_phi_pT_Aside,           "zy",1,0,"#phi (rad)","1ptResCov_phi_HighPt_TPCAside");
  if(h1ptResTPCInnerC_vs_phi_pT_Aside->GetEntries())   Plot1PtRes(h1ptResTPCInnerC_vs_phi_pT_Aside,  "zy",1,0,"#phi (rad)","1ptResTPCInnerCCov_phi_HighPt_TPCAside");
  if(h1ptResTPCInner_vs_phi_pT_Aside->GetEntries())    Plot1PtRes(h1ptResTPCInner_vs_phi_pT_Aside,   "zy",1,0,"#phi (rad)","1ptResTPCInnerCov_phi_HighPt_TPCAside");
  if(h1ptSigma_vs_phi_pT_Aside->GetEntries())          Plot1PtRes(h1ptSigma_vs_phi_pT_Aside,         "zy",3,0,"#phi (rad)","1ptSigmaCov_phi_HighPt_TPCAside");
  if(h1ptSigmaTPCInnerC_vs_phi_pT_Aside->GetEntries()) Plot1PtRes(h1ptSigmaTPCInnerC_vs_phi_pT_Aside,"zy",3,0,"#phi (rad)","1ptSigmaTPCInnerCCov_phi_HighPt_TPCAside");
  if(h1ptSigmaTPCInner_vs_phi_pT_Aside->GetEntries())  Plot1PtRes(h1ptSigmaTPCInner_vs_phi_pT_Aside, "zy",3,0,"#phi (rad)","1ptSigmaTPCInnerCov_phi_HighPt_TPCAside");
  //Cside
  if(hPulldcaR_vs_phi_pT_Cside->GetEntries())          Plot2D(hPulldcaR_vs_phi_pT_Cside,             "zy",0,0,"#phi (rad)","dcaRPull_phi_HighPt_TPCCside");
  if(hPulldcaRTPCInner_vs_phi_pT_Cside->GetEntries())  Plot2D(hPulldcaRTPCInner_vs_phi_pT_Cside,     "zy",0,0,"#phi (rad)","dcaRPullTPCInner_phi_HighPt_TPCCside");
  if(hResdcaR_vs_phi_pT_Cside->GetEntries())           Plot2D(hResdcaR_vs_phi_pT_Cside,              "zy",1,0,"#phi (rad)","dcaRRes_phi_HighPt_TPCCside");
  if(hResdcaRTPCInner_vs_phi_pT_Cside->GetEntries())   Plot2D(hResdcaRTPCInner_vs_phi_pT_Cside,      "zy",2,0,"#phi (rad)","dcaRResTPCInner_phi_HighPt_TPCCside");
  if(hphiPull_vs_phi_pT_Cside->GetEntries())           Plot2D(hphiPull_vs_phi_pT_Cside,              "zy",7,0,"#phi (rad)","phiPullTPCInner_phi_HighPt_TPCCside");
  if(hphiRes_vs_phi_pT_Cside->GetEntries())            Plot2D(hphiRes_vs_phi_pT_Cside,               "zy",8,0,"#phi (rad)","phiResTPCInner_phi_HighPt_TPCCside");

  if(h1ptRes_vs_phi_pT_Cside->GetEntries())            Plot1PtRes(h1ptRes_vs_phi_pT_Cside,           "zy",1,0,"#phi (rad)","1ptResCov_phi_HighPt_TPCCside");
  if(h1ptResTPCInnerC_vs_phi_pT_Cside->GetEntries())   Plot1PtRes(h1ptResTPCInnerC_vs_phi_pT_Cside,  "zy",1,0,"#phi (rad)","1ptResTPCInnerCCov_phi_HighPt_TPCCside");
  if(h1ptResTPCInner_vs_phi_pT_Cside->GetEntries())    Plot1PtRes(h1ptResTPCInner_vs_phi_pT_Cside,   "zy",1,0,"#phi (rad)","1ptResTPCInnerCov_phi_HighPt_TPCCside");
  if(h1ptSigma_vs_phi_pT_Cside->GetEntries())          Plot1PtRes(h1ptSigma_vs_phi_pT_Cside,         "zy",3,0,"#phi (rad)","1ptSigmaCov_phi_HighPt_TPCCside");
  if(h1ptSigmaTPCInnerC_vs_phi_pT_Cside->GetEntries()) Plot1PtRes(h1ptSigmaTPCInnerC_vs_phi_pT_Cside,"zy",3,0,"#phi (rad)","1ptSigmaTPCInnerCCov_phi_HighPt_TPCCside");
  if(h1ptSigmaTPCInner_vs_phi_pT_Cside->GetEntries())  Plot1PtRes(h1ptSigmaTPCInner_vs_phi_pT_Cside, "zy",3,0,"#phi (rad)","1ptSigmaTPCInnerCov_phi_HighPt_TPCCside");
  //V0 pos
  if(hK0sPull_vs_alpha_1pT_pos->GetEntries())          Plot2DK0s(hK0sPull_vs_alpha_1pT_pos,          "zy",0,0,"#phi (rad)","pullK0sPos_phi_HighPt");
  if(hK0sRes_vs_alpha_1pT_pos->GetEntries())           Plot2DK0s(hK0sRes_vs_alpha_1pT_pos,           "zy",1,0,"#phi (rad)","resK0sPos_phi_HighPt");
  //V0 neg
  if(hK0sPull_vs_alpha_1pT_neg->GetEntries())          Plot2DK0s(hK0sPull_vs_alpha_1pT_neg,          "zy",0,0,"#phi (rad)","pullK0sNeg_phi_HighPt");
  if(hK0sRes_vs_alpha_1pT_neg->GetEntries())           Plot2DK0s(hK0sRes_vs_alpha_1pT_neg,           "zy",1,0,"#phi (rad)","resK0sNeg_phi_HighPt");

  if(heta_phi_pT->GetEntries())                        Plot1D(heta_phi_pT,"y",1,0,"#phi (rad)","phi_HighPt");
  if(heta_phi_pT->GetEntries())                        Plot1D(heta_phi_pT,"x",1,0,"#eta","eta_HighPt");

  if(hphi_vs_eta_pT_cutTPC->GetEntries()&&hphi_vs_eta_pT_cutTPCITS->GetEntries()) PlotEff(hphi_vs_eta_pT_cutTPCITS,hphi_vs_eta_pT_cutTPC,"x",1,0,"#phi (rad)","TPCITSMatchingEff_phi_HighPt");
  if(hphi_vs_eta_pT_cutTPC->GetEntries()&&hphi_vs_eta_pT_cutTPCITS->GetEntries()) PlotEff(hphi_vs_eta_pT_cutTPCITS,hphi_vs_eta_pT_cutTPC,"y",1,0,"#eta","TPCITSMatchingEff_eta_HighPt");    

  if(heta_phi_pT->GetEntries()){
    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    TH3D *h1c = (TH3D*) heta_phi_pT->Clone(Form("%s_1",heta_phi_pT->GetName()));
    h1c->GetZaxis()->SetRangeUser(fPtCut,100);
    TH2D *heta_vs_phi = (TH2D*) h1c->Project3D("yx");
    heta_vs_phi->GetYaxis()->SetTitle("#phi (rad)");
    heta_vs_phi->GetXaxis()->SetTitle("#eta");
    heta_vs_phi->SetTitle("#eta vs. #phi (pT > 3 GeV)");
    TCanvas *can = new TCanvas("can","can",550,550);
    can->cd();
    //  gPad->SetLogz();
    //   gPad->SetLeftMargin(0.16);
    //   gPad->SetBottomMargin(0.15);
    heta_vs_phi->Draw("COLZ");
    can->SaveAs("eta_vs_phi.png");
    if(can) delete can;
    if(h1c) delete h1c;
    if(heta_vs_phi) delete heta_vs_phi;
  }

  //vs Pt
  fPtCut = 0;
  //Aside
  if(hPulldcaR_vs_phi_pT_Aside->GetEntries())          Plot2D(hPulldcaR_vs_phi_pT_Aside,             "zx",0,1,"p_{T} (GeV/c)","dcaRPull_pT_TPCAside");
  if(hPulldcaRTPCInner_vs_phi_pT_Aside->GetEntries())  Plot2D(hPulldcaRTPCInner_vs_phi_pT_Aside,     "zx",0,1,"p_{T} (GeV/c)","dcaRPullTPCInner_pT_TPCAside");
  if(hResdcaR_vs_phi_pT_Aside->GetEntries())           Plot2D(hResdcaR_vs_phi_pT_Aside,              "zx",1,1,"p_{T} (GeV/c)","dcaRRes_pT_TPCAside");
  if(hResdcaRTPCInner_vs_phi_pT_Aside->GetEntries())   Plot2D(hResdcaRTPCInner_vs_phi_pT_Aside,      "zx",2,1,"p_{T} (GeV/c)","dcaRResTPCInner_pT_TPCAside");
  if(hphiPull_vs_phi_pT_Aside->GetEntries())           Plot2D(hphiPull_vs_phi_pT_Aside,              "zx",7,1,"p_{T} (GeV/c)","phiPullTPCInner_pT_TPCAside");
  if(hphiRes_vs_phi_pT_Aside->GetEntries())            Plot2D(hphiRes_vs_phi_pT_Aside,               "zx",8,1,"p_{T} (GeV/c)","phiResTPCInner_pT_TPCAside");

  if(h1ptRes_vs_phi_pT_Aside->GetEntries())            Plot1PtRes(h1ptRes_vs_phi_pT_Aside,           "zx",0,1,"p_{T} (GeV/c)","1ptResCov_pT_TPCAside");
  if(h1ptResTPCInnerC_vs_phi_pT_Aside->GetEntries())   Plot1PtRes(h1ptResTPCInnerC_vs_phi_pT_Aside,  "zx",0,1,"p_{T} (GeV/c)","1ptResTPCInnerCCov_pT_TPCAside");
  if(h1ptResTPCInner_vs_phi_pT_Aside->GetEntries())    Plot1PtRes(h1ptResTPCInner_vs_phi_pT_Aside,   "zx",0,1,"p_{T} (GeV/c)","1ptResTPCInnerCov_pT_TPCAside");
  if(h1ptSigma_vs_phi_pT_Aside->GetEntries())          Plot1PtRes(h1ptSigma_vs_phi_pT_Aside,         "zx",2,1,"p_{T} (GeV/c)","1ptSigmaCov_pT_TPCAside");
  if(h1ptSigmaTPCInnerC_vs_phi_pT_Aside->GetEntries()) Plot1PtRes(h1ptSigmaTPCInnerC_vs_phi_pT_Aside,"zx",2,1,"p_{T} (GeV/c)","1ptSigmaTPCInnerCCov_pT_TPCAside");
  if(h1ptSigmaTPCInner_vs_phi_pT_Aside->GetEntries())  Plot1PtRes(h1ptSigmaTPCInner_vs_phi_pT_Aside, "zx",2,1,"p_{T} (GeV/c)","1ptSigmaTPCInnerCov_pT_TPCAside");
  //Cside
  if(hPulldcaR_vs_phi_pT_Cside->GetEntries())          Plot2D(hPulldcaR_vs_phi_pT_Cside,             "zx",0,1,"p_{T} (GeV/c)","dcaRPull_pT_TPCCside");
  if(hPulldcaRTPCInner_vs_phi_pT_Cside->GetEntries())  Plot2D(hPulldcaRTPCInner_vs_phi_pT_Cside,     "zx",0,1,"p_{T} (GeV/c)","dcaRPullTPCInner_pT_TPCCside");
  if(hResdcaR_vs_phi_pT_Cside->GetEntries())           Plot2D(hResdcaR_vs_phi_pT_Cside,              "zx",1,1,"p_{T} (GeV/c)","dcaRRes_pT_TPCCside");
  if(hResdcaRTPCInner_vs_phi_pT_Cside->GetEntries())   Plot2D(hResdcaRTPCInner_vs_phi_pT_Cside,      "zx",2,1,"p_{T} (GeV/c)","dcaRResTPCInner_pT_TPCCside");
  if(hphiPull_vs_phi_pT_Cside->GetEntries())           Plot2D(hphiPull_vs_phi_pT_Cside,              "zx",7,1,"p_{T} (GeV/c)","phiPullTPCInner_pT_TPCCside");
  if(hphiRes_vs_phi_pT_Cside->GetEntries())            Plot2D(hphiRes_vs_phi_pT_Cside,               "zx",8,1,"p_{T} (GeV/c)","phiResTPCInner_pT_TPCCside");

  if(h1ptRes_vs_phi_pT_Cside->GetEntries())            Plot1PtRes(h1ptRes_vs_phi_pT_Cside,           "zx",0,1,"p_{T} (GeV/c)","1ptResCov_pT_TPCCside");
  if(h1ptResTPCInnerC_vs_phi_pT_Cside->GetEntries())   Plot1PtRes(h1ptResTPCInnerC_vs_phi_pT_Cside,  "zx",0,1,"p_{T} (GeV/c)","1ptResTPCInnerCCov_pT_TPCCside");
  if(h1ptResTPCInner_vs_phi_pT_Cside->GetEntries())    Plot1PtRes(h1ptResTPCInner_vs_phi_pT_Cside,   "zx",0,1,"p_{T} (GeV/c)","1ptResTPCInnerCov_pT_TPCCside");
  if(h1ptSigma_vs_phi_pT_Cside->GetEntries())          Plot1PtRes(h1ptSigma_vs_phi_pT_Cside,         "zx",2,1,"p_{T} (GeV/c)","1ptSigmaCov_pT_TPCCside");
  if(h1ptSigmaTPCInnerC_vs_phi_pT_Cside->GetEntries()) Plot1PtRes(h1ptSigmaTPCInnerC_vs_phi_pT_Cside,"zx",2,1,"p_{T} (GeV/c)","1ptSigmaTPCInnerCCov_pT_TPCCside");
  if(h1ptSigmaTPCInner_vs_phi_pT_Cside->GetEntries())  Plot1PtRes(h1ptSigmaTPCInner_vs_phi_pT_Cside, "zx",2,1,"p_{T} (GeV/c)","1ptSigmaTPCInnerCov_pT_TPCCside");

  //V0 pos
  if(hK0sPull_vs_alpha_1pT_pos->GetEntries())          Plot2DK0s(hK0sPull_vs_alpha_1pT_pos,          "zx",0,0,"1/p_{T} (1/GeV/c)","pullK0sPos_1pT");
  if(hK0sRes_vs_alpha_1pT_pos->GetEntries())           Plot2DK0s(hK0sRes_vs_alpha_1pT_pos,           "zx",1,0,"1/p_{T} (1/GeV/c)","resK0sPos_1pT");
  //V0 neg
  if(hK0sPull_vs_alpha_1pT_neg->GetEntries())          Plot2DK0s(hK0sPull_vs_alpha_1pT_neg,          "zx",0,0,"1/p_{T} (1/GeV/c)","pullK0sNeg_1pT");
  if(hK0sRes_vs_alpha_1pT_neg->GetEntries())           Plot2DK0s(hK0sRes_vs_alpha_1pT_neg,           "zx",1,0,"1/p_{T} (1/GeV/c)","resK0sNeg_1pT");

  if(heta_phi_pT->GetEntries())                        Plot1D(heta_phi_pT,"z",0,1,"p_{T} (GeV/c)","Pt");
  if(hphi_vs_eta_pT_cutTPC->GetEntries()&&hphi_vs_eta_pT_cutTPCITS->GetEntries()) PlotEff(hphi_vs_eta_pT_cutTPCITS,hphi_vs_eta_pT_cutTPC,"z",1,1,"p_{T} (GeV/c)","TPCITSMatchingEff_pT");
}

Bool_t AliHighPtTreeAnalysis::GetK0TrendFitFunction(TF1 *fLinearFitK0sShift, TF1 *fLinearFitK0sSigma, Int_t Type, Int_t Charge){
  // Type = 0: resolution
  // Type = 1: pull
  // Charge =  1: pos
  // Charge = -1: neg
  gStyle->SetOptStat(0);
  char name[256];
  TH3D *h1c = 0;
  if(Type == 0 && Charge == -1) h1c = (TH3D*) hK0sRes_vs_alpha_1pT_neg->Clone(Form("%s_1",hK0sRes_vs_alpha_1pT_neg->GetName()));
  if(Type == 0 && Charge ==  1) h1c = (TH3D*) hK0sRes_vs_alpha_1pT_pos->Clone(Form("%s_1",hK0sRes_vs_alpha_1pT_pos->GetName()));
  if(Type == 1 && Charge == -1) h1c = (TH3D*) hK0sPull_vs_alpha_1pT_neg->Clone(Form("%s_1",hK0sPull_vs_alpha_1pT_neg->GetName()));
  if(Type == 1 && Charge ==  1) h1c = (TH3D*) hK0sPull_vs_alpha_1pT_pos->Clone(Form("%s_1",hK0sPull_vs_alpha_1pT_pos->GetName()));
  if(!h1c) return kFALSE;

  TH2D *h1proj = (TH2D*)h1c->Project3D("zx");
  if(!h1proj) return kFALSE;
  TH1D *h1mean = (TH1D*)h1c->Project3D("x");
  if(!h1mean) return kFALSE;
  h1mean->SetName("h1mean");
  h1mean->Reset();
  TH1D *h1width = (TH1D*)h1c->Project3D("x");
  if(!h1width) return kFALSE;
  h1width->SetName("h1width");
  h1width->Reset();

  TF1 * fK0sFit = new TF1("fK0sFit","[0]+[1]*x+[2]*TMath::Gaus(x,[3],[4])",-0.5,0.5);
  Int_t nBinsX = h1proj->GetXaxis()->GetNbins();
  Int_t NbinsFitted(0);
  for(Int_t i = 1; i<nBinsX; i++){
    h1proj->GetXaxis()->SetRange(i,i);
    if(h1proj->Integral()<15) continue;
    sprintf(name,"h1K0sproj_%s_%d",h1c->GetName(),i);
    TH1D *h1K0sproj = (TH1D *)h1proj->ProjectionY(name);
    Double_t y1 = h1K0sproj->GetBinContent(5);
    Double_t y2 = h1K0sproj->GetBinContent(95);
    fK0sFit->SetParameter(0,(y1+y2)*0.5);
    fK0sFit->SetParameter(1,(y2-y1)/20.);
    fK0sFit->SetParameter(2,h1K0sproj->GetMaximum());
    fK0sFit->SetParameter(3,0);
    fK0sFit->SetParameter(4,h1K0sproj->GetRMS());
    h1K0sproj->Fit(fK0sFit,"Q"); 
    h1K0sproj->Fit(fK0sFit,"Q"); 

    h1mean->SetBinContent(i,fK0sFit->GetParameter(3));
    h1mean->SetBinError(i,fK0sFit->GetParError(3));
    h1width->SetBinContent(i,TMath::Abs(fK0sFit->GetParameter(4)));
    h1width->SetBinError(i,fK0sFit->GetParError(4));
    NbinsFitted++;
  }

  if(NbinsFitted < 3) return kFALSE;
  h1mean ->Fit(fLinearFitK0sShift,"Q");
  h1width->Fit(fLinearFitK0sSigma,"Q");

  h1mean ->SetMarkerStyle(20);
  h1mean ->SetMarkerColor(1);
  h1mean ->SetLineColor(1);

  h1width->SetMarkerStyle(20);
  h1width->SetMarkerColor(1);
  h1width->SetLineColor(1);

  TCanvas *can = new TCanvas("can","can",550,700);
  can->Divide(1,2);

  can->cd(1);
  gPad->SetLeftMargin(0.16);
  gPad->SetBottomMargin(0.15);
  h1mean->Draw();

  can->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  h1width->Draw();

  sprintf(name,"%s.png",h1c->GetName());
  can->SaveAs(name);

  delete can;

  return kTRUE;
}

Int_t AliHighPtTreeAnalysis::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

Long64_t AliHighPtTreeAnalysis::LoadTree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;

  return centry;
}

Long64_t AliHighPtTreeAnalysis::LoadV0Tree(Long64_t entry)
{
  // Set the environment to read one entry
  if (!fV0Chain) return -5;
  Long64_t centry = fV0Chain->LoadTree(entry);
  if (centry < 0) return centry;

  return centry;
}

void AliHighPtTreeAnalysis::InitV0tree(TTree *tree)
{
  // Set branch addresses and branch pointers

  if(!tree) return;
  fMakeV0s=kTRUE;
  fV0Chain = tree;

  fV0Chain->SetBranchAddress("v0.", &v0);
  fV0Chain->SetBranchAddress("track0.", &v0track0);
  fV0Chain->SetBranchAddress("track1.", &v0track1);

}

void AliHighPtTreeAnalysis::Init(TTree *tree)
{
  // Set branch addresses and branch pointers
  if (!tree) return;

  fChain = tree;

  TString str(fChain->GetBranch("runNumber")->GetTitle());
  if(str[str.Length()-1]=='I')
    fChain->SetBranchAddress("runNumber", &fRunNumberInt);
  if(str[str.Length()-1]=='D')
    fChain->SetBranchAddress("runNumber", &runNumber);
  str = fChain->GetBranch("Bz")->GetTitle();
  if(str[str.Length()-1]=='I')
    fChain->SetBranchAddress("Bz", &BzInt);
  if(str[str.Length()-1]=='D')
    fChain->SetBranchAddress("Bz", &Bz);

  fChain->SetBranchAddress("esdTrack.", &esdTrack);
  fChain->SetBranchAddress("vtxESD.", &vtxESD);
  //   fChain->SetBranchAddress("runNumber", &runNumber);
  fChain->SetBranchAddress("extTPCInnerC.", &extTPCInnerC);
  fChain->SetBranchAddress("chi2TPCInnerC", &chi2TPCInnerC);
  fChain->SetBranchAddress("mult", &fMult);
  //   fChain->SetBranchAddress("Bz", &Bz);
  if( fChain->GetBranchStatus("particle.") ){
    fHasMC = kTRUE;
    fChain->SetBranchAddress("particle.", &particle);
  }
  else
    fHasMC = kFALSE;
}

Bool_t AliHighPtTreeAnalysis::ConnectGenericHistos( const char *genericHistoFile )
{
  TFile *f = TFile::Open(genericHistoFile);
  if(!f) return kFALSE;
  if(strstr(genericHistoFile,"Bpos")) fBfield = 1;
  else if(strstr(genericHistoFile,"Bneg")) fBfield = -1;
  else fBfield = 0;


  heta_phi_pT                        = (TH3D*) f->Get("heta_phi_pT");
  hphi_vs_eta_pT_cutTPC              = (TH3D*) f->Get("hphi_vs_eta_pT_cutTPC");
  hphi_vs_eta_pT_cutTPCITS           = (TH3D*) f->Get("hphi_vs_eta_pT_cutTPCITS");

  h1pt_vs_eta_phi                    = (TH3F*) f->Get("h1pt_vs_eta_phi");
  h1ptTPCInner_vs_eta_phi            = (TH3F*) f->Get("h1ptTPCInner_vs_eta_phi");
  h1ptTPCInnerC_vs_eta_phi           = (TH3F*) f->Get("h1ptTPCInnerC_vs_eta_phi");

  hPulldcaR_vs_eta_pT_Aside          = (TH3D*) f->Get("hPulldcaR_vs_eta_pT_Aside");           hPulldcaR_vs_eta_pT_Cside          = (TH3D*) f->Get("hPulldcaR_vs_eta_pT_Cside");
  hPulldcaR_vs_phi_pT_Aside          = (TH3D*) f->Get("hPulldcaR_vs_phi_pT_Aside");           hPulldcaR_vs_phi_pT_Cside          = (TH3D*) f->Get("hPulldcaR_vs_phi_pT_Cside");
  hPulldcaRTPCInner_vs_eta_pT_Aside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_eta_pT_Aside");   hPulldcaRTPCInner_vs_eta_pT_Cside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_eta_pT_Cside");
  hPulldcaRTPCInner_vs_phi_pT_Aside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_phi_pT_Aside");   hPulldcaRTPCInner_vs_phi_pT_Cside  = (TH3D*) f->Get("hPulldcaRTPCInner_vs_phi_pT_Cside");
  hResdcaR_vs_eta_pT_Aside           = (TH3D*) f->Get("hResdcaR_vs_eta_pT_Aside");            hResdcaR_vs_eta_pT_Cside           = (TH3D*) f->Get("hResdcaR_vs_eta_pT_Cside");
  hResdcaR_vs_phi_pT_Aside           = (TH3D*) f->Get("hResdcaR_vs_phi_pT_Aside");            hResdcaR_vs_phi_pT_Cside           = (TH3D*) f->Get("hResdcaR_vs_phi_pT_Cside");
  hResdcaRTPCInner_vs_eta_pT_Aside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_eta_pT_Aside");    hResdcaRTPCInner_vs_eta_pT_Cside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_eta_pT_Cside");
  hResdcaRTPCInner_vs_phi_pT_Aside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_phi_pT_Aside");    hResdcaRTPCInner_vs_phi_pT_Cside   = (TH3D*) f->Get("hResdcaRTPCInner_vs_phi_pT_Cside");

  hphiPull_vs_phi_pT_Aside           = (TH3D*) f->Get("hphiPull_vs_phi_pT_Aside");            hphiPull_vs_phi_pT_Cside           = (TH3D*) f->Get("hphiPull_vs_phi_pT_Cside");
  hphiRes_vs_phi_pT_Aside            = (TH3D*) f->Get("hphiRes_vs_phi_pT_Aside");             hphiRes_vs_phi_pT_Cside            = (TH3D*) f->Get("hphiRes_vs_phi_pT_Cside");

  h1ptRes_vs_phi_pT_Aside            = (TH3D*) f->Get("h1ptRes_vs_phi_pT_Aside");             h1ptRes_vs_phi_pT_Cside            = (TH3D*) f->Get("h1ptRes_vs_phi_pT_Cside");
  h1ptRes_vs_mult_pT_Aside           = (TH3D*) f->Get("h1ptRes_vs_mult_pT_Aside");            h1ptRes_vs_mult_pT_Cside           = (TH3D*) f->Get("h1ptRes_vs_mult_pT_Cside");
  h1ptSigma_vs_phi_pT_Aside          = (TH3D*) f->Get("h1ptSigma_vs_phi_pT_Aside");           h1ptSigma_vs_phi_pT_Cside          = (TH3D*) f->Get("h1ptSigma_vs_phi_pT_Cside");
  h1ptSigma_vs_mult_pT_Aside         = (TH3D*) f->Get("h1ptSigma_vs_mult_pT_Aside");          h1ptSigma_vs_mult_pT_Cside         = (TH3D*) f->Get("h1ptSigma_vs_mult_pT_Cside");
  h1ptResTPCInnerC_vs_phi_pT_Aside   = (TH3D*) f->Get("h1ptResTPCInnerC_vs_phi_pT_Aside");    h1ptResTPCInnerC_vs_phi_pT_Cside   = (TH3D*) f->Get("h1ptResTPCInnerC_vs_phi_pT_Cside");
  h1ptResTPCInnerC_vs_mult_pT_Aside  = (TH3D*) f->Get("h1ptResTPCInnerC_vs_mult_pT_Aside");   h1ptResTPCInnerC_vs_mult_pT_Cside  = (TH3D*) f->Get("h1ptResTPCInnerC_vs_mult_pT_Cside");
  h1ptSigmaTPCInnerC_vs_phi_pT_Aside = (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_phi_pT_Aside");  h1ptSigmaTPCInnerC_vs_phi_pT_Cside = (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_phi_pT_Cside");
  h1ptSigmaTPCInnerC_vs_mult_pT_Aside= (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_mult_pT_Aside"); h1ptSigmaTPCInnerC_vs_mult_pT_Cside= (TH3D*) f->Get("h1ptSigmaTPCInnerC_vs_mult_pT_Cside");
  h1ptResTPCInner_vs_phi_pT_Aside    = (TH3D*) f->Get("h1ptResTPCInner_vs_phi_pT_Aside");     h1ptResTPCInner_vs_phi_pT_Cside    = (TH3D*) f->Get("h1ptResTPCInner_vs_phi_pT_Cside");
  h1ptResTPCInner_vs_mult_pT_Aside   = (TH3D*) f->Get("h1ptResTPCInner_vs_mult_pT_Aside");    h1ptResTPCInner_vs_mult_pT_Cside   = (TH3D*) f->Get("h1ptResTPCInner_vs_mult_pT_Cside");
  h1ptSigmaTPCInner_vs_phi_pT_Aside  = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_phi_pT_Aside");   h1ptSigmaTPCInner_vs_phi_pT_Cside  = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_phi_pT_Cside");
  h1ptSigmaTPCInner_vs_mult_pT_Aside = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_mult_pT_Aside");  h1ptSigmaTPCInner_vs_mult_pT_Cside = (TH3D*) f->Get("h1ptSigmaTPCInner_vs_mult_pT_Cside");

  hK0sPull_vs_alpha_1pT_pos = (TH3D*) f->Get("hK0sPull_vs_alpha_1pT_pos");
  hK0sRes_vs_alpha_1pT_pos  = (TH3D*) f->Get("hK0sRes_vs_alpha_1pT_pos");
  hK0sPull_vs_alpha_1pT_neg = (TH3D*) f->Get("hK0sPull_vs_alpha_1pT_neg");
  hK0sRes_vs_alpha_1pT_neg  = (TH3D*) f->Get("hK0sRes_vs_alpha_1pT_neg");

  hPulldcaRTPConly_vs_eta_1pT  = (TH3D*) f->Get("hPulldcaRTPConly_vs_eta_1pT");
  hPulldcaRcomb_vs_eta_1pT     = (TH3D*) f->Get("hPulldcaRcomb_vs_eta_1pT");
  hResdcaRTPConly_vs_eta_1pT   = (TH3D*) f->Get("hResdcaRTPConly_vs_eta_1pT");
  hResdcaRcomb_vs_eta_1pT      = (TH3D*) f->Get("hResdcaRcomb_vs_eta_1pT");

  hphiPull_vs_eta_1pT          = (TH3D*) f->Get("hphiPull_vs_eta_1pT");
  hphiRes_vs_eta_1pT           = (TH3D*) f->Get("hphiRes_vs_eta_1pT");

  if(hK0sPull_vs_alpha_1pT_pos && hK0sRes_vs_alpha_1pT_pos && hK0sPull_vs_alpha_1pT_neg && hK0sRes_vs_alpha_1pT_neg) fMakeV0s = kTRUE;

  return kTRUE;
}

void AliHighPtTreeAnalysis::SetMakePlots(Bool_t makeAllPlots){ fMakePlots = makeAllPlots; }

void AliHighPtTreeAnalysis::SetApplyCorrections( const char *correctionFile )
{
  //
  // 
  //
  fApplyCorrections = kTRUE;

  TFile *fCorr = TFile::Open( correctionFile );
  TTree *CorrectionTree;
  fCorr->GetObject("PeriodTree",CorrectionTree);


  fCorrectionAside          = new Double_t[18];
  fCorrectionCside          = new Double_t[18];
  fCorrectionAsideTPCInner  = new Double_t[18];
  fCorrectionCsideTPCInner  = new Double_t[18];
  fCorrectionAsideTPCInnerC = new Double_t[18];
  fCorrectionCsideTPCInnerC = new Double_t[18];

  TH1D *hCorrtmp = new TH1D("hname","htitle",100,-1.,1.);
  for(Int_t i = 0; i != 18; ++i){
    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAside.fElements[%d]",i));
    fCorrectionAside[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();

    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCside.fElements[%d]",i));
    fCorrectionCside[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();

    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAsideTPCInner.fElements[%d]",i));
    fCorrectionAsideTPCInner[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();

    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCsideTPCInner.fElements[%d]",i));
    fCorrectionCsideTPCInner[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();

    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCAsideTPCInnerC.fElements[%d]",i));
    fCorrectionAsideTPCInnerC[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();

    CorrectionTree->Project("hname",Form("vec1ptShiftFitSec_TPCCsideTPCInnerC.fElements[%d]",i));
    fCorrectionCsideTPCInnerC[i] = hCorrtmp->GetMean();
    hCorrtmp->Reset();     
  }
}
