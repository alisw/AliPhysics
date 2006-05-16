#include "AliESDtrackCuts.h"

#include <Riostream.h>

//____________________________________________________________________
ClassImp(AliESDtrackCuts)

//____________________________________________________________________
AliESDtrackCuts::AliESDtrackCuts() {

  //##############################################################################
  // setting default cuts

  SetMinNClustersTPC();
  SetMinNClustersITS();	    
  SetMaxChi2PerClusterTPC();
  SetMaxChi2PerClusterITS();  				    
  SetMaxCovDiagonalElements();  				    
  SetRequireTPCRefit();
  SetRequireITSRefit();
  SetAcceptKingDaughters();
  SetMinNsigmaToVertex();
  SetRequireSigmaToVertex();
  SetPRange();
  SetPtRange();
  SetPxRange();
  SetPyRange();
  SetPzRange();
  SetEtaRange();
  SetRapRange();

  SetHistogramsOn();

  // set the cut names
  fCutNames[0]  = "require TPC refit";
  fCutNames[1]  = "require ITS refit";  
  fCutNames[2]  = "n clusters TPC";
  fCutNames[3]  = "n clusters ITS";
  fCutNames[4]  = "#Chi^{2}/clusters TPC";
  fCutNames[5]  = "#Chi^{2}/clusters ITS";
  fCutNames[6]  = "cov 11";
  fCutNames[7]  = "cov 22";
  fCutNames[8]  = "cov 33";
  fCutNames[9]  = "cov 44";
  fCutNames[10] = "cov 55";
  fCutNames[11] = "trk-to-vtx";
  fCutNames[12] = "trk-to-vtx failed";
  fCutNames[13] = "kink daughters";

  fCutNames[14] = "p";
  fCutNames[15] = "p_{T}";
  fCutNames[16] = "p_{x}";
  fCutNames[17] = "p_{y}";
  fCutNames[18] = "p_{z}";
  fCutNames[19] = "y";
  fCutNames[20] = "eta";

}

//____________________________________________________________________
Bool_t 
AliESDtrackCuts::AcceptTrack(AliESDtrack* esdTrack) {
  // 
  // figure out if the tracks survives all the track cuts defined
  //

  UInt_t status = esdTrack->GetStatus();
  
  // getting quality parameters from the ESD track
  Int_t nClustersITS = esdTrack->GetITSclusters(fIdxInt);
  Int_t nClustersTPC = esdTrack->GetTPCclusters(fIdxInt);
  
  Float_t chi2PerClusterITS = -1;
  Float_t chi2PerClusterTPC = -1;
  if (nClustersITS!=0)
    chi2PerClusterITS = esdTrack->GetITSchi2()/Float_t(nClustersITS);
  if (nClustersTPC!=0)
    chi2PerClusterTPC = esdTrack->GetTPCchi2()/Float_t(nClustersTPC);  

  Double_t extCov[15];
  esdTrack->GetExternalCovariance(extCov);  

  // getting the track to vertex parameters
  Float_t b[2];
  Float_t bRes[2];
  Float_t bCov[3];
  esdTrack->GetImpactParameters(b,bCov);    
  if (bCov[0]<=0 || bCov[2]<=0) {
    AliDebug(1, "Estimated b resolution zero!");
    bCov[0]=0; bCov[1]=0;
  }
  bRes[0] = TMath::Sqrt(bCov[0]);
  bRes[1] = TMath::Sqrt(bCov[2]);

  // FIX !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  //
  // this is not correct - it will not give n sigma!!!
  // 
  Float_t nSigmaToVertex = -1;
  if (bRes[0]!=0 && bRes[1]!=0)
    nSigmaToVertex = TMath::Sqrt(TMath::Power(b[0]/bRes[0],2) + TMath::Power(b[1]/bRes[1],2));  

  // getting the kinematic variables of the track 
  // (assuming the mass is known)
  Double_t p[3];
  esdTrack->GetPxPyPz(p);
  Float_t momentum = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2) + TMath::Power(p[2],2));
  Float_t pt       = TMath::Sqrt(TMath::Power(p[0],2) + TMath::Power(p[1],2));
  Float_t energy   = TMath::Sqrt(TMath::Power(esdTrack->GetMass(),2) + TMath::Power(momentum,2));

  //y-eta related calculations
  Float_t eta = -100.;
  Float_t y   = -100.;
  if((momentum != TMath::Abs(p[2]))&&(momentum != 0))
    eta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));
  if((energy != TMath::Abs(p[2]))&&(momentum != 0))
    y = 0.5*TMath::Log((energy + p[2])/(energy - p[2]));

  
  //########################################################################
  // cut the track?
  
  Bool_t cuts[fNCuts];
  for (Int_t i=0; i<fNCuts; i++) cuts[i]=kFALSE;
  
  // track quality cuts
  if (fCut_RequireTPCRefit && (status&AliESDtrack::kTPCrefit)==0)
    cuts[0]=kTRUE;
  if (fCut_RequireITSRefit && (status&AliESDtrack::kITSrefit)==0)
    cuts[1]=kTRUE;
  if (nClustersTPC<fCut_MinNClusterTPC) 
    cuts[2]=kTRUE;
  if (nClustersITS<fCut_MinNClusterITS) 
    cuts[3]=kTRUE;
  if (chi2PerClusterTPC>fCut_MaxChi2PerClusterTPC) 
    cuts[4]=kTRUE; 
  if (chi2PerClusterITS>fCut_MaxChi2PerClusterITS) 
    cuts[5]=kTRUE;
  if (extCov[0]  > fCut_MaxC11) 
    cuts[6]=kTRUE;  
  if (extCov[2]  > fCut_MaxC22) 
    cuts[7]=kTRUE;  
  if (extCov[5]  > fCut_MaxC33) 
    cuts[8]=kTRUE;  
  if (extCov[9]  > fCut_MaxC44) 
    cuts[9]=kTRUE;  
  if (extCov[14]  > fCut_MaxC55) 
    cuts[10]=kTRUE;  
  if (nSigmaToVertex > fCut_NsigmaToVertex) 
    cuts[11] = kTRUE;
  // if n sigma could not be calculated
  if (nSigmaToVertex<0 && fCut_SigmaToVertexRequired)   
    cuts[12]=kTRUE;
  if (!fCut_AcceptKinkDaughters && esdTrack->GetKinkIndex(0)>0) 
    cuts[13]=kTRUE;
  // track kinematics cut
  if((momentum < fPMin) || (momentum > fPMax)) 
    cuts[14]=kTRUE;
  if((pt < fPtMin) || (pt > fPtMax)) 
    cuts[15] = kTRUE;
  if((p[0] < fPxMin) || (p[0] > fPxMax)) 
    cuts[16] = kTRUE;
  if((p[1] < fPyMin) || (p[1] > fPyMax)) 
    cuts[17] = kTRUE;
  if((p[2] < fPzMin) || (p[2] > fPzMax)) 
    cuts[18] = kTRUE;
  if((eta < fEtaMin) || (eta > fEtaMax)) 
    cuts[19] = kTRUE;
  if((y < fRapMin) || (y > fRapMax)) 
    cuts[20] = kTRUE;

  Bool_t cut=kFALSE;
  for (Int_t i=0; i<fNCuts; i++) 
    if (cuts[i]) cut = kTRUE;
  
  //########################################################################
  // filling histograms
  if (fHistogramsOn) {
    fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n tracks")));
    
    if (cut)
      fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin("n cut tracks")));
    
    for (Int_t i=0; i<fNCuts; i++) {
      if (cuts[i])
 	fhCutStatistics->Fill(fhCutStatistics->GetBinCenter(fhCutStatistics->GetXaxis()->FindBin(fCutNames[i])));
      
      for (Int_t j=i; j<fNCuts; j++) {
 	if (cuts[i] && cuts[j]) {
 	  Float_t x = fhCutCorrelation->GetXaxis()->GetBinCenter(fhCutCorrelation->GetXaxis()->FindBin(fCutNames[i]));
 	  Float_t y = fhCutCorrelation->GetYaxis()->GetBinCenter(fhCutCorrelation->GetYaxis()->FindBin(fCutNames[j]));
 	  fhCutCorrelation->Fill(x,y);
 	}
      }
    }
    

    fhNClustersITS[0]->Fill(nClustersITS);        
    fhNClustersTPC[0]->Fill(nClustersTPC);        
    fhChi2PerClusterITS[0]->Fill(chi2PerClusterITS);
    fhChi2PerClusterTPC[0]->Fill(chi2PerClusterTPC);   
    
    fhC11[0]->Fill(extCov[0]);                 
    fhC22[0]->Fill(extCov[2]);                 
    fhC33[0]->Fill(extCov[5]);                 
    fhC44[0]->Fill(extCov[9]);                                  
    fhC55[0]->Fill(extCov[14]);                                  
    
    fhDZ[0]->Fill(b[1]);     
    fhDXY[0]->Fill(b[0]);    
    fhDXYvsDZ[0]->Fill(b[1],b[0]);

    if (bRes[0]!=0 && bRes[1]!=0) {
      fhDZNormalized[0]->Fill(b[1]/bRes[1]);     
      fhDXYNormalized[0]->Fill(b[0]/bRes[0]);    
      fhDXYvsDZNormalized[0]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
    }
  }

  //########################################################################  
  // cut the track!
  if (cut) return kFALSE;

  //########################################################################  
  // filling histograms after cut
  if (fHistogramsOn) {
    fhNClustersITS[1]->Fill(nClustersITS);        
    fhNClustersTPC[1]->Fill(nClustersTPC);        
    fhChi2PerClusterITS[1]->Fill(chi2PerClusterITS);
    fhChi2PerClusterTPC[1]->Fill(chi2PerClusterTPC);   
    
    fhC11[1]->Fill(extCov[0]);                 
    fhC22[1]->Fill(extCov[2]);                 
    fhC33[1]->Fill(extCov[5]);                 
    fhC44[1]->Fill(extCov[9]);                                  
    fhC55[1]->Fill(extCov[14]);                                  
    
    fhDZ[1]->Fill(b[1]);     
    fhDXY[1]->Fill(b[0]);    
    fhDXYvsDZ[1]->Fill(b[1],b[0]);

    fhDZNormalized[1]->Fill(b[1]/bRes[1]);     
    fhDXYNormalized[1]->Fill(b[0]/bRes[0]);    
    fhDXYvsDZNormalized[1]->Fill(b[1]/bRes[1], b[0]/bRes[0]);
  }
  
  return kTRUE;
}

//____________________________________________________________________
TObjArray*
AliESDtrackCuts::GetAcceptedTracks(AliESD* esd) {
  
  // returns an array of all tracks that pass the cuts
  fAcceptedTracks->Clear();
  
  // loop over esd tracks
  for (Int_t iTrack = 0; iTrack < esd->GetNumberOfTracks(); iTrack++) {
    AliESDtrack* track = esd->GetTrack(iTrack);
    
    if(AcceptTrack(track)) fAcceptedTracks->Add(track);
  }

  return fAcceptedTracks;
}

//____________________________________________________________________
 void AliESDtrackCuts::DefineHistograms(Int_t color) {

   fHistogramsOn=kTRUE;

//   //###################################################################################
   // defining histograms

   fhCutStatistics = new TH1F("cut_statistics","cut statistics",fNCuts+4,-0.5,fNCuts+3.5);

   fhCutStatistics->GetXaxis()->SetBinLabel(1,"n tracks");
   fhCutStatistics->GetXaxis()->SetBinLabel(2,"n cut tracks");

   fhCutCorrelation = new TH2F("cut_correlation","cut correlation",fNCuts,-0.5,fNCuts-0.5,fNCuts,-0.5,fNCuts-0.5);;
  
   for (Int_t i=0; i<fNCuts; i++) {
     fhCutStatistics->GetXaxis()->SetBinLabel(i+4,fCutNames[i]);
     fhCutCorrelation->GetXaxis()->SetBinLabel(i+1,fCutNames[i]);
     fhCutCorrelation->GetYaxis()->SetBinLabel(i+1,fCutNames[i]);
   } 

  fhCutStatistics  ->SetLineColor(color);
  fhCutCorrelation ->SetLineColor(color);
  fhCutStatistics  ->SetLineWidth(2);
  fhCutCorrelation ->SetLineWidth(2);


  TH1F *fhNClustersITS  = new TH1F[2];
  TH1F *fhNClustersTPC  = new TH1F[2];
  TH1F *fhChi2PerClusterITS   = new TH1F[2];
  TH1F *fhChi2PerClusterTPC   = new TH1F[2];
  		       
  TH1F *fhC11                 = new TH1F[2];
  TH1F *fhC22                 = new TH1F[2];
  TH1F *fhC33                 = new TH1F[2];
  TH1F *fhC44                 = new TH1F[2];
  TH1F *fhC55                 = new TH1F[2];
  
  TH1F *fhDXY                 = new TH1F[2];
  TH1F *fhDZ		       = new TH1F[2];
  TH2F *fhDXYvsDZ	       = new TH2F[2];

  TH1F *fhDXYNormalized       = new TH1F[2];
  TH1F *fhDZNormalized        = new TH1F[2];
  TH2F *fhDXYvsDZNormalized   = new TH2F[2];


  Char_t str[256];
  for (Int_t i=0; i<2; i++) {
    if (i==0) sprintf(str," ");
    else sprintf(str,"_cut");

    fhNClustersITS[i]        = TH1F(Form("nClustersITS%s",str),"",8,-0.5,7.5);
    fhNClustersTPC[i]        = TH1F(Form("nClustersTPC%s",str),"",165,-0.5,164.5);
    fhChi2PerClusterITS[i]   = TH1F(Form("chi2PerClusterITS%s",str),"",500,0,10);
    fhChi2PerClusterTPC[i]   = TH1F(Form("chi2PerClusterTPC%s",str),"",500,0,10);
    
    fhC11[i]                 = TH1F(Form("covMatrixDiagonal11%s",str),"",1000,0,5);
    fhC22[i]                 =  TH1F(Form("covMatrixDiagonal22%s",str),"",1000,0,5);
    fhC33[i]                 =  TH1F(Form("covMatrixDiagonal33%s",str),"",1000,0,0.5);
    fhC44[i]                 =  TH1F(Form("covMatrixDiagonal44%s",str),"",1000,0,5);
    fhC55[i]                 =  TH1F(Form("covMatrixDiagonal55%s",str),"",1000,0,5);
    
    fhDXY[i]                 =  TH1F(Form("dXY%s",str),"",500,-10,10);
    fhDZ[i]                  =  TH1F(Form("dZ%s",str),"",500,-10,10);
    fhDXYvsDZ[i]             =  TH2F(Form("dXYvsDZ%s",str),"",200,-10,10,200,-10,10);

    fhDXYNormalized[i]       =  TH1F(Form("dXYNormalized%s",str),"",500,-10,10);
    fhDZNormalized[i]        =  TH1F(Form("dZNormalized%s",str),"",500,-10,10);
    fhDXYvsDZNormalized[i]   =  TH2F(Form("dXYvsDZNormalized%s",str),"",200,-10,10,200,-10,10);


    fhNClustersITS[i].SetXTitle("n ITS clusters");  
    fhNClustersTPC[i].SetXTitle("n TPC clusters"); 
    fhChi2PerClusterITS[i].SetXTitle("#Chi^{2} per ITS cluster"); 
    fhChi2PerClusterTPC[i].SetXTitle("#Chi^{2} per TPC cluster"); 
    
    fhC11[i].SetXTitle("cov 11 : #sigma_{y}^{2} [cm^{2}]"); 
    fhC22[i].SetXTitle("cov 22 : #sigma_{z}^{2} [cm^{2}]"); 
    fhC33[i].SetXTitle("cov 33 : #sigma_{sin(#phi)}^{2}"); 
    fhC44[i].SetXTitle("cov 44 : #sigma_{tan(#theta_{dip})}^{2}"); 
    fhC55[i].SetXTitle("cov 55 : #sigma_{1/p_{T}}^{2} [(c/GeV)^2]"); 
   
    fhDXY[i].SetXTitle("transverse impact parameter"); 
    fhDZ[i].SetXTitle("longitudinal impact parameter"); 
    fhDXYvsDZ[i].SetXTitle("longitudinal impact parameter"); 
    fhDXYvsDZ[i].SetYTitle("transverse impact parameter"); 

    fhDXYNormalized[i].SetXTitle("normalized trans impact par"); 
    fhDZNormalized[i].SetXTitle("normalized long impact par"); 
    fhDXYvsDZNormalized[i].SetXTitle("normalized long impact par"); 
    fhDXYvsDZNormalized[i].SetYTitle("normalized trans impact par"); 

    fhNClustersITS[i].SetLineColor(color);   fhNClustersITS[i].SetLineWidth(2);
    fhNClustersTPC[i].SetLineColor(color);   fhNClustersTPC[i].SetLineWidth(2);
    fhChi2PerClusterITS[i].SetLineColor(color);   fhChi2PerClusterITS[i].SetLineWidth(2);
    fhChi2PerClusterTPC[i].SetLineColor(color);   fhChi2PerClusterTPC[i].SetLineWidth(2);
    						     					      
    fhC11[i].SetLineColor(color);   fhC11[i].SetLineWidth(2);
    fhC22[i].SetLineColor(color);   fhC22[i].SetLineWidth(2);
    fhC33[i].SetLineColor(color);   fhC33[i].SetLineWidth(2);
    fhC44[i].SetLineColor(color);   fhC44[i].SetLineWidth(2);
    fhC55[i].SetLineColor(color);   fhC55[i].SetLineWidth(2);
   						     					      
    fhDXY[i].SetLineColor(color);   fhDXY[i].SetLineWidth(2);
    fhDZ[i].SetLineColor(color);   fhDZ[i].SetLineWidth(2);
						     
    fhDXYNormalized[i].SetLineColor(color);   fhDXYNormalized[i].SetLineWidth(2);
    fhDZNormalized[i].SetLineColor(color);   fhDZNormalized[i].SetLineWidth(2);

  }
}

//____________________________________________________________________
void 
AliESDtrackCuts::Print(const Option_t*) const {

  AliInfo("AliESDtrackCuts...");
}


//____________________________________________________________________
void AliESDtrackCuts::SaveHistograms(Char_t* dir) {
  
  if (!fHistogramsOn) {
    AliDebug(0, "Histograms not on - cannot save histograms!!!");
    return;
  }

  gDirectory->mkdir(dir);
  gDirectory->cd(dir);

  gDirectory->mkdir("before_cuts");
  gDirectory->mkdir("after_cuts");
 
  fhCutStatistics->Write();
  fhCutCorrelation->Write();

  for (Int_t i=0; i<2; i++) {
    if (i==0)
      gDirectory->cd("before_cuts");
    else
      gDirectory->cd("after_cuts");
    
    fhNClustersITS[i]        ->Write();
    fhNClustersTPC[i]        ->Write();
    fhChi2PerClusterITS[i]   ->Write();
    fhChi2PerClusterTPC[i]   ->Write();
    
    fhC11[i]                 ->Write();
    fhC22[i]                 ->Write();
    fhC33[i]                 ->Write();
    fhC44[i]                 ->Write();
    fhC55[i]                 ->Write();

    fhDXY[i]                 ->Write();
    fhDZ[i]                  ->Write();
    fhDXYvsDZ[i]             ->Write();
    
    fhDXYNormalized[i]       ->Write();
    fhDZNormalized[i]        ->Write();
    fhDXYvsDZNormalized[i]   ->Write();

    gDirectory->cd("../");
  }

  gDirectory->cd("../");
}



