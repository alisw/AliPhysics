// $Id$
/*
 * Modified version of $ALICE_ROOT/TPC/AliTPCComparison.C for evaluating the
 * performance of the CA tracker. 
 * 
 * Usage:
 * <pre>
 *   aliroot AliHLTTPCTrackerEvaluation.C
 * </pre>
 *
 * documentation to be filled in when the macro is finalized 
 *
 *
 * @ingroup alihlt_tpc
 */
#if defined(__CINT__) && !defined(__MAKECINT__)

//Int_t AliHLTTPCTrackerEvaluation(const Char_t *dir=".", const char* input="AliHLTESDs.root", Float_t ptLowerCut=0.1, Float_t ptUpperCut=10., Bool_t bDedxAndClus=kFALSE)
{
  gSystem->AddIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE_ROOT/TPC -I$ALICE/geant3/TGeant3");
  gROOT->LoadMacro("$ALICE_ROOT/HLT/TPCLib/macros/AliHLTTPCTrackerEvaluation.C++");
 
  AliHLTTPCTrackerEvaluation();
  // the attempt to pass the arguments with a self-called macro hasn't worked so far
  //AliHLTTPCTrackerEvaluation(dir, input, ptLowerCut, ptUpperCut, bDedxAndClus);
 
}
#else

#include <AliLog.h>
#include <TMath.h>
#include <TError.h>
#include <Riostream.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TTree.h>
#include <TParticle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TText.h>
#include <TBenchmark.h>
#include <TStyle.h>
#include <TFile.h>
#include <TROOT.h>

#include "AliStack.h"
#include "AliHeader.h"
#include "AliTrackReference.h"
#include "AliRunLoader.h"
#include "AliRun.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliTPCtrack.h"

#include "AliSimDigits.h"
#include "AliTPC.h"
#include "AliTPCParamSR.h"
#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"
#include "AliTPCcluster.h"
#include "AliTPCLoader.h"
#include "TParticlePDG.h"
#include "TDatabasePDG.h"

#include "AliCDBManager.h"
#include "AliGRPManager.h"
#include "AliTPCcalibDB.h"
#include "TProfile.h"

//_______________________________________________________________________________________


Int_t GoodTPCTracks(const Char_t *dir = ".");

extern AliRun     *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT      *gROOT;

static Int_t allgood     = 0;
static Int_t allselected = 0;
static Int_t allfound    = 0;



Int_t AliHLTTPCTrackerEvaluation(const Char_t *dir=".", const char* input="AliHLTESDs.root", Float_t ptLowerCut=0.1, Float_t ptUpperCut=10., Bool_t bDedxAndClus=kFALSE){
 
   if(!input){
       cerr << "please specify an input file" << endl;
       return 1;
   }

   gBenchmark->Start("AliHLTTPCTrackerEvaluation");

   AliCDBManager* man = AliCDBManager::Instance();
   man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
   man->SetRun(0);
   man->SetSpecificStorage("GRP/GRP/Data","local://$PWD");
   AliGRPManager grpman;
   if (!grpman.ReadGRPEntry() ||
       !grpman.SetMagField()) {
     cerr << "can not set magnetic field from GRP" << endl;
     return -1;
   }

   ::Info("AliHLTTPCTrackerEvaluation.C","Calculating reconstruction efficiency...");

   gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include -I$ALICE_ROOT -I$ALICE/geant3/TGeant3");

   gStyle->SetTitleOffset(1.2,"X");
   gStyle->SetTitleOffset(1.1,"Y");

   
   // --------------- Efficiency histograms ------------------------

   TH1F *hClus      = new TH1F("hClus","# clusters on reconstructed tracks",300,0,300);  
  
   TH1F *hGood      = new TH1F("hGood",     "P_{t} distribution of MC tracks (selected for efficiency calculation)", 30, 0, 6);
   TH1F *hFake      = new TH1F("hFake",     "P_{t} distribution of fake tracks",                                     30, 0, 6);
   TH1F *hFound     = new TH1F("hFound",    "P_{t} distribution of reconstructed tracks",                            30, 0, 6);
   TH1F *hFoundMult = new TH1F("hFoundMult","P_{t} distribution of reconstructed tracks with clones",                30, 0, 6);
   TH1F *hClone     = new TH1F("hClone",    "P_{t} distribution of clone tracks",                                    30, 0, 6);
 
   TH1F *hEff      = new TH1F("hEff",     "Reconstruction efficiency based on selected MC tracks", 30, 0, 6);
   TH1F *hFakeEff  = new TH1F("hFakeEff", "Efficiency for fake tracks",                            30, 0, 6);
   TH1F *hCloneEff = new TH1F("hCloneEff","Efficiency for clone tracks",                           30, 0, 6);
 
   hEff->SetLineWidth(2);
   hEff->SetMaximum(1.4);
   hEff->SetYTitle("efficiency"); 
   hEff->GetYaxis()->CenterTitle();
   hEff->SetXTitle("P_{t} (GeV/c)");
 
   hFakeEff->SetLineColor(kRed); 
   hFakeEff->SetLineWidth(2);
 
   hCloneEff->SetLineColor(kBlue); 
   hCloneEff->SetLineWidth(2);
 
   TH2F *hDedx = new TH2F("hDedx","dE/dx vs. momentum", 50,0.,2.,50,0.,400.);
   hDedx->SetMarkerStyle(8);
   hDedx->SetMarkerSize(0.4);
   hDedx->SetXTitle("P (Gev/c)"); 
   hDedx->SetYTitle("dE/dx (a.u.)");

   TH1F *hnhit_ref = new TH1F("hnhit_ref","# clusters on reco tracks, used in fit performance",300,0,300);
   TH1F *he = new TH1F("he","dE/dx for pions with 0.4<p<0.5 GeV/c",50,0.,100.);

   // --------------- Resolution histograms ---------------------------
   
   TH2F *hPt     = new TH2F("hPt",     "", 30, 0, 6, 50,  -3,  3); 
   TH2F *hPhi    = new TH2F("hPhi",    "", 30, 0, 6, 50, -16, 16); 
   TH2F *hLambda = new TH2F("hLambda", "", 30, 0, 6, 50, -10, 10);
   TH2F *hY	 = new TH2F("hY",      "", 30, 0, 6, 30,  -5,  5); 
   TH2F *hZ	 = new TH2F("hZ",      "", 30, 0, 250, 30,  -5,  5); // fix the x axis, if this is z dependent 
   
   TProfile *hResPt     = new TProfile("hResPt",     "P_{t} resolution (%)",  30, 0, 6);
   TProfile *hResPhi    = new TProfile("hResPhi",    "#phi resolution (mrad)",   30, 0, 6);    
   TProfile *hResLambda = new TProfile("hResLambda", "#lambda resolution (mrad)",30, 0, 6);
   TProfile *hResY      = new TProfile("hResY",      "Y resolution (mm)",      30, 0, 6); 
   TProfile *hResZ      = new TProfile("hResZ",      "Z resolution (mm)",      30, 0, 250); // select appropriate x axis, KKK

   hResPt->SetXTitle("P_{t} (GeV/c)");
   hResPt->SetYTitle("#sigma_{(P_{t}-P_{t_{MC}})/P_{t_{MC}}} (%)");
   hResPt->GetYaxis()->CenterTitle(true);
   hResPt->GetYaxis()->CenterTitle(true);
  
   hResPhi->SetXTitle("P_{t} (GeV/c)");
   hResPhi->SetYTitle("#sigma_{(#phi_{rec}-#phi_{MC})} (mrad)");
   hResPhi->GetYaxis()->CenterTitle(true);

   hResLambda->SetXTitle("P_{t} (GeV/c)");
   hResLambda->SetYTitle("#sigma_{(#lambda_{rec}-#lambda_{MC})} (mrad)"); // KKK perhaps we should add the definition on the histo pad
   hResLambda->GetYaxis()->CenterTitle(true);
  
   hResY->SetXTitle("P_{t} (GeV/c)");
   hResY->SetYTitle("#sigma_{(Y_{rec}-Y_{MC})} (mm)");
   hResY->GetYaxis()->CenterTitle(true);
   
   hResZ->SetXTitle("Z (mm)");
   hResZ->SetYTitle("#sigma_{(Z_{rec}-Z_{MC})} (mm)");
   hResZ->GetYaxis()->CenterTitle(true);

   // KKK I leave the Z coordinate up to you. You can make the hResZ Z dependent.
   
   /* 
     
   TH1D *hProjPt[15], *hProjPhi[15], *hProjLambda[15], *hProjY[15], *hProjZ[15]; // I am not sure we need these histos to be pointers
   Char_t name[15];  Char_t title[100];
   
   for(Int_t i=0; i<15; i++){
       sprintf(name,"hProjPt%i",i+1);
       sprintf(title,"(Pt_{MC}-Pt_{Rec})/Pt_{MC} @ Pt#in[%f, %f] GeV/c", i-0.5, i+0.5);
       hProjPt[i] = new TH1D(name, title, 50, -3, -3); 
       
       sprintf(name,"hProjPhi%i",i+1);
       sprintf(title,"(#phi_{MC}-#phi_{Rec}) @ #phi#in[%f, %f] GeV/c", i-0.5, i+0.5);
       hProjPhi[i] = new TH1D(name, title, 50, -16, -16);        
      
       sprintf(name,"hProjLambda%i",i+1);
       sprintf(title,"(#lambda_{MC}-#lambda_{Rec}) @ #lambda#in[%f, %f] GeV/c", i-0.5, i+0.5);
       hProjLambda[i] = new TH1D(name, title, 50, -16, -16);        
       
       sprintf(name,"hProjY%i",i+1);
       sprintf(title,"(Y_{MC}-Y_{Rec}) @ Y#in[%f, %f] GeV/c", i-0.5, i+0.5);
       hProjY[i] = new TH1D(name, title, 50, -16, -16);   
       
       // here you add the hProjZ[] histograms    KKK 
   }
   */ 
   
     
   // --------------- Pull variable histograms ------------------------

   TH1F *hpullPhi  = new TH1F("hpullPhi", "SinPhi pull", 30,-10.,10.); 
   TH1F *hpullY    = new TH1F("hpullY",   "Y pull",      30,-10.,10.); 
   TH1F *hpullZ    = new TH1F("hpullZ",   "Z pull",      30,-10.,10.); 
   TH1F *hpullDzds = new TH1F("hpullDzds","Dzds pull",   30,-10.,10.); 
   TH1F *hpullK    = new TH1F("hpullK",   "Kappa pull",  30,-10.,10.); 

//---------------------------------------------------------------------------------------

   Char_t fname[100];
   sprintf(fname,"%s/GoodTPCTracks.root",dir);

   TFile *refFile = TFile::Open(fname,"old");
   if(!refFile || !refFile->IsOpen()){
      ::Info("AliHLTTPCTrackerEvaluation.C","Marking good tracks (will take a while)...");
      if(GoodTPCTracks(dir)){
         ::Error("AliHLTTPCTrackerEvaluation.C","Cannot generate the reference file!");
         return 1;
     }
   }
  
   refFile = TFile::Open(fname,"old");
   if(!refFile || !refFile->IsOpen()){
      ::Error("AliHLTTPCTrackerEvaluation.C","Cannot open the reference file!");
      return 2;
   }   
   
   //------------ access the contents of the tree created by the function GoodTracks() -------------
  
   TTree *tpcTree = (TTree*)refFile->Get("tpcTree");
   if(!tpcTree){
      ::Error("AliHLTTPCTrackerEvaluation.C","Cannot get the reference tree!");
      return 3;
   }
   
   TBranch *branch = tpcTree->GetBranch("TPC");
   if(!branch){
      ::Error("AliHLTTPCTrackerEvaluation.C","Cannot get the TPC branch!");
      return 4;
   }
  
   TClonesArray dummy("AliTrackReference",1000), *refs = &dummy;
   branch->SetAddress(&refs);

   //------------ access the HLT output --------------------

   sprintf(fname,"%s/%s", dir, input);
   TFile *ef = TFile::Open(fname);
   
   if( !ef || !ef->IsOpen() ){
      ::Error("AliHLTTPCTrackerEvaluation.C","Cannot open [%s]!",fname);
      return 5;   
   }
   
   AliESDEvent *event = new AliESDEvent();
   TTree *esdTree = NULL;
   TString tsinput = input;
   
   if(tsinput.CompareTo("AliESDs.root")==0){
      esdTree = (TTree*)ef->Get("HLTesdTree");
      if(!esdTree){
         ::Error("AliHLTTPCTrackerEvaluation.C", "no HLTESD tree found");
         return 6;
      }
   event->ReadFromTree(esdTree);
   } else if(tsinput.CompareTo("AliHLTESDs.root")==0){
      esdTree = (TTree*)ef->Get("esdTree");
      if(!esdTree){
         ::Error("AliHLTTPCTrackerEvaluation.C", "no ESD tree found");
         return 7;
      }
   event->ReadFromTree(esdTree);
   } else return 8;
      

   //---------- Loop over events -------------------

   Int_t iEvent = 0;
   while(esdTree->GetEvent(iEvent)){
     
      Int_t nEntries = event->GetNumberOfTracks();      
      cout << "********* Processing event number: " << iEvent << ", " << nEntries << " reconstructed tracks *******\n" << endl;

      allfound += nEntries;

      if(tpcTree->GetEvent(iEvent++)==0){
	 cerr << "No reconstructable tracks !\n";
         continue;
      }

      Int_t ngood = refs->GetEntriesFast(); // access the track references
      cout << ngood << " good MC tracks" << endl;
      allgood += ngood;

      const Int_t MAX = 15000;
  
      Int_t track_notfound[MAX],                            itrack_notfound   = 0;
      Int_t track_fake[MAX],                                itrack_fake       = 0;
      Int_t track_multifound[MAX], track_multifound_n[MAX], itrack_multifound = 0;

      for(Int_t i=0; i<nEntries; i++){
	  
	  hClus->Fill(event->GetTrack(i)->GetTPCNcls()); // filled but not displayed
	  cout << "TPC label for track " << i << " = " << event->GetTrack(i)->GetTPCLabel() << ", # clusters " << event->GetTrack(i)->GetTPCNcls() << endl;
      }
      
      Int_t i;
      for(Int_t k=0; k<ngood; k++){ // read the MC information

  	  AliTrackReference *ref = (AliTrackReference*)refs->UncheckedAt(k); 
         
	  Int_t  label = ref->Label();
	  Int_t tlabel = -1;
          Float_t ptMC = TMath::Sqrt(ref->Px()*ref->Px() + ref->Py()*ref->Py());
	  Float_t pMC  = TMath::Sqrt(ref->Px()*ref->Px() + ref->Py()*ref->Py()+ref->Pz()*ref->Pz());
	
          if (ptMC<1e-33 || ptMC<ptLowerCut || ptMC>ptUpperCut) continue; // first condition is for tracks not crossing padrow 0

          allselected++;

          hGood->Fill(ptMC);
	  
	  for(i=0; i<nEntries; i++){		   
	      Int_t ttlabel = event->GetTrack(i)->GetTPCLabel();
	      if( label==TMath::Abs(ttlabel) && ttlabel<0 ){
	         track_fake[itrack_fake++] = label;
	         hFake->Fill(ptMC);
	      }	      
	  }
	
          AliESDtrack *esdtrack = 0;
          for(i=0; i<nEntries; i++){            
	    esdtrack = event->GetTrack(i);
	    tlabel   = esdtrack->GetTPCLabel();
	    if(label==tlabel) break;
          }
        
	  if(i==nEntries){
	    track_notfound[itrack_notfound++] = label;
	    cout << "MC track " << label << " not reconstructed" << endl;
	    continue;
          }
      
          Int_t micount = 0;
          Int_t mi;
          AliESDtrack *mitrack;

          for(mi=0; mi<nEntries; mi++){
	      mitrack = event->GetTrack(mi);	      
	      if(label==mitrack->GetTPCLabel()) micount++;
          }
	
          if(micount>1){
	     track_multifound[itrack_multifound]   = label;
	     track_multifound_n[itrack_multifound] = micount;
	     itrack_multifound++;
	     hClone->Fill(ptMC,micount-1); 	    
          }
        
	//if((mitrack->GetStatus()&AliESDtrack::kTPCrefit)==0) continue;
        //if((mitrack->GetStatus()&AliESDtrack::kTPCin)==0)    continue;
	  cout << "MC track " << label << " reconstructed "<<micount<<" times" << endl;
	  
        if(label==tlabel){
	   hFound->Fill(ptMC);
	   //if( micount<1 ) cout<<"ERROR!!!!"<<endl;
	   hFoundMult->Fill(ptMC,micount);
	}	
	
	AliExternalTrackParam tpctrack =*(esdtrack->GetInnerParam());

	if( TMath::Abs(tpctrack.GetSigned1Pt()) <1.e-8 ) continue;
	
	Double_t mcxyz[3]   = { ref->X(), ref->Y(), ref->Z() };
	Double_t mclocal[3] = { ref->X(), ref->Y(), ref->Z() };
	
	if(1){
	  Double_t c = TMath::Cos(tpctrack.GetAlpha());
	  Double_t s = TMath::Sin(tpctrack.GetAlpha());
	  Double_t x = mclocal[0];
	  mclocal[0] = x*c + mclocal[1]*s;
	  mclocal[1] =-x*s + mclocal[1]*c;
	}
	
	if(0){	  
	  Double_t local[3] = { tpctrack.GetX(),tpctrack.GetY(),tpctrack.GetZ() };
// 	  cout << "label: "      << label << endl;
// 	  cout << "Z diff = "    << local[2] - mclocal[2] << endl;
// 	  cout << "orig rec: "   << local[0]<<" "<<local[1]<<" "<<local[2] <<" "<<TMath::Sqrt(local[0]*local[0]+local[1]*local[1])<<endl;
// 	  cout << "rotated mc: " << mclocal[0]<<" "<<mclocal[1]<<" "<<mclocal[2] <<" "<<TMath::Sqrt(mclocal[0]*mclocal[0]+mclocal[1]*mclocal[1])<<endl;
// 	  cout << "orig mc: "    << mcxyz[0]<<" "<<mcxyz[1]<<" "<<mcxyz[2]<<endl;
// 	  cout << "alpha = "     << tpctrack.GetAlpha()/TMath::Pi()*180.<<" "<<ref->Alpha()/TMath::Pi()*180.<<endl;
	}

	//cout << "X = " << mclocal[0] << " / " << tpctrack.GetX() << endl;	
	tpctrack.AliExternalTrackParam::PropagateTo(mclocal[0],5.);
       
        Double_t xyz[3], pxpypz[3]; 

	tpctrack.GetXYZ(xyz);       // GetInnerXYZ
        tpctrack.GetPxPyPz(pxpypz); // GetInnerPxPyPz
	
	Double_t local[3] = { tpctrack.GetX(),tpctrack.GetY(),tpctrack.GetZ() };
        Float_t phiRec = TMath::ATan2(pxpypz[1],pxpypz[0]);
        
	if(phiRec <= TMath::Pi()) phiRec += 2*TMath::Pi();
        if(phiRec >= TMath::Pi()) phiRec -= 2*TMath::Pi();
        
	Double_t ptRec = TMath::Sqrt(pxpypz[0]*pxpypz[0]+pxpypz[1]*pxpypz[1]);
	Double_t pRec  = TMath::Sqrt(pxpypz[0]*pxpypz[0]+pxpypz[1]*pxpypz[1]+pxpypz[2]*pxpypz[2]);
	
        Float_t  lambdaRec = TMath::ATan2(pxpypz[2],ptRec); 
        Float_t pt_1 = 1./ptRec;
	//Float_t pts_1 = tpctrack->GetSigned1Pt();
	//Float_t ptMCs_1 = ref->/ptMC;

	

        Int_t pdg = (Int_t)ref->GetLength();  // particle PDG
	const Double_t kCLight = 0.000299792458;  
	//Double_t Bz = 5.;
	Double_t q  = 0.;
	
	if( TMath::Abs(pdg)<10000 ){
	   TParticlePDG *ppdg = TDatabasePDG::Instance()->GetParticle(pdg);
	   if( ppdg ) q = ppdg->Charge()/3.;
	}
	
	Double_t kappaMC = q/ptMC; // /Bz/kCLight;
	hnhit_ref->Fill(esdtrack->GetTPCNcls());    
 
	hPt->Fill(ptMC, (ptRec - ptMC)/(ptMC)*100.);	       
	
	Float_t phiMC = TMath::ATan2(ref->Py(),ref->Px());
        hPhi->Fill(ptMC, (phiRec - phiMC)/phiMC*1000.);

        Float_t lambdaMC = TMath::ATan2(ref->Pz(),ptMC);
        hLambda->Fill(ptMC, (lambdaRec - lambdaMC)/lambdaMC*1000.);    
	
	hY->Fill(ptMC, (local[1]-mclocal[1])*10.);
	hZ->Fill( fabs(mclocal[2]), (local[2]-mclocal[2] ) *10.); // Please check hY and hZ! KKK
	
	if( tpctrack.GetSigmaY2()>0   && finite(tpctrack.GetSigmaY2())   )  hpullY   ->Fill( (local[1]-mclocal[1])/TMath::Sqrt(TMath::Abs(tpctrack.GetSigmaY2()))  );
	if( tpctrack.GetSigmaZ2()>0   && finite(tpctrack.GetSigmaZ2())   )  hpullZ   ->Fill( (local[2]-mclocal[2])/TMath::Sqrt(TMath::Abs(tpctrack.GetSigmaZ2()))  );
	if( tpctrack.GetSigmaSnp2()>0 && finite(tpctrack.GetSigmaSnp2()) )  hpullPhi ->Fill( (phiRec-phiMC)/TMath::Sqrt(TMath::Abs(tpctrack.GetSigmaSnp2()))	   );
	if( tpctrack.GetSigmaTgl2()>0 && finite(tpctrack.GetSigmaTgl2()) )  hpullDzds->Fill( (lambdaRec-lambdaMC)/TMath::Sqrt(TMath::Abs(tpctrack.GetSigmaTgl2())) );
	
	if( tpctrack.GetSigma1Pt2()>0 && finite(tpctrack.GetSigma1Pt2()) )
	if( q!=0 ) hpullK->Fill((tpctrack.GetSigned1Pt()-kappaMC)/TMath::Sqrt(TMath::Abs(tpctrack.GetSigma1Pt2())));	
 	
        Float_t dedx = esdtrack->GetTPCsignal();	
        hDedx->Fill(pRec,dedx,1.);
        
	if(TMath::Abs(pdg)==211) // pions
	   if(pRec>0.4 && pRec<0.5){
              he->Fill(dedx,1.);
           }
      } // end of loop over "good" tracks
            
      /*
      cout<<"\nList of Not found tracks :\n";
      for ( i = 0; i< itrack_notfound; i++){
        cout<<track_notfound[i]<<"\t";
        if ((i%5)==4) cout<<"\n";
      }
      cout<<"\nList of fake  tracks :\n";
      for ( i = 0; i< itrack_fake; i++){
        cout<<track_fake[i]<<"\t";
        if ((i%5)==4) cout<<"\n";
      }
      cout<<"\nList of multiple found tracks :\n";
      for ( i=0; i<itrack_multifound; i++) {
          cout<<"id.   "<<track_multifound[i]
              <<"     found - "<<track_multifound_n[i]<<"times\n";
      }

      cout<<"Number of found tracks ="<<nEntries<<endl;
      cout<<"Number of \"good\" tracks ="<<ngood<<endl;
      */
      
  refs->Clear();
  } // end of the loop over events

   delete event;
   delete esdTree;
   ef->Close();

   delete tpcTree;
   refFile->Close();

   Stat_t nGoodMC  = hGood->GetEntries();
   Stat_t nFakes   = hFake->GetEntries();
   Stat_t nClones  = hClone->GetEntries();
   Stat_t nRec     = hFound->GetEntries();
   Stat_t nRecMult = hFoundMult->GetEntries();
   
   
    if(nGoodMC!=0)         ::Info("\n\nAliHLTTPCTrackerEvaluation","Integral efficiency in (%) is about : %f\n", nRec/nGoodMC*100.); 
    if(nRecMult+nFakes!=0) ::Info("AliHLTTPCTrackerEvaluation","Integral fake rate in (%) is about  : %f\n", nFakes/(nRecMult+nFakes)*100.);
    if(nRecMult+nFakes!=0) ::Info("AliHLTTPCTrackerEvaluation","Integral clone rate in (%) is about : %f\n", nClones/(nRecMult+nFakes)*100.);
    
    ::Info("AliHLTTPCTrackerEvaluation", "Total number of selected MC tracks                    : %d\n", allgood);
    ::Info("AliHLTTPCTrackerEvaluation", "Total number of selected MC tracks with momentum cuts : %d\n", allselected);
    ::Info("AliHLTTPCTrackerEvaluation", "Total number of reconstructed tracks                  : %d\n", allfound);


   // ---------- Plotting the histograms ----------------------------------------------------------------------------
   
   hFound    ->Sumw2(); 
   hFoundMult->Sumw2(); 
   hGood     ->Sumw2(); 
   hFake     ->Sumw2(); 
   hClone    ->Sumw2();
  
   hEff     ->Divide(hFound,hGood,     1,1.,"b");
   hFakeEff ->Divide(hFake, hGood,     1,1.,"b");
   hCloneEff->Divide(hClone,hFoundMult,1,1.,"b"); 
   // KKK are you sure this is the correct definition? When plotted, it looks like we have 70% clones in the lower momenta...
   // perhaps fakes and clones should be displayed in another histogram??

   gROOT->SetStyle("Plain");
   gStyle->SetOptStat(""); //gStyle->SetOptStat(111110);
   gStyle->SetOptFit(0);  //gStyle->SetOptFit(1);
   gStyle->SetPalette(1);
 

   //--------------- canvas 0 for EFFICIENCY -------------------
   
   TCanvas *c0 = new TCanvas("c0","reconstruction efficiency",500,450);
   
   c0->cd();
   hEff     ->Draw();
   hCloneEff->Draw("histsame");
   hFakeEff ->Draw("histsame");
   
   TLegend *leg = new TLegend(0.3,0.5,0.85,0.8);
   leg->SetFillColor(10);
   leg->SetLineColor(10);
   leg->AddEntry(hEff,      "reconstruction eff.", "l");
   leg->AddEntry(hCloneEff, "clone contribution", "l");
   leg->AddEntry(hFakeEff,  "fake contribution",  "l");
   leg->Draw("same");
      
   TLine *line1 = new TLine(0,1,6,1); 
   line1->SetLineStyle(4);
   line1->Draw("same");
   
   Float_t fakeData  = nFakes/(nRecMult+nFakes)*100.;
   Float_t cloneData = nClones/(nRecMult+nFakes)*100.;
   Char_t  fakeBuf[100];
   Char_t  cloneBuf[100];

   sprintf(fakeBuf,"%d", (Int_t)fakeData);
   TString fakeStr = fakeBuf;
   fakeStr.Append(" %");
   
   sprintf(cloneBuf,"%d", (Int_t)cloneData);
   TString cloneStr = cloneBuf;
   cloneStr.Append(" %");
     
   TText *text = new TText(1.4,0.1,fakeStr.Data());
   text->SetTextSize(0.05);
   text->SetTextColor(kRed);
   text->Draw();
   text = new TText(0.6,0.3,cloneStr.Data());
   text->SetTextSize(0.05);
   text->SetTextColor(kBlue);
   text->Draw();
   
   
   //--------------- canvas 1 for RESOLUTIONS ----------------

   
   gStyle->SetOptStat(0);
   TCanvas *c1 = new TCanvas("c1","resolutions",1100,900); // please add the Y and Z projections

   TF1 *fGauss = new TF1("gauss","gaus");   
   for(Int_t i=0; i<15; i++){              
     Float_t pMin = i*6/15.;// -.5;
     Float_t pMax = (i+1)*6/15.;// +.5;
     Int_t binx1 = hPt->GetXaxis()->FindBin(pMin);
     Int_t binx2 = hPt->GetXaxis()->FindBin(pMax);
     Float_t zMin = i*250./15.;
     Float_t zMax = (i+1)*250./15.;
     Int_t binz1 = hZ->GetXaxis()->FindBin(zMin);
     Int_t binz2 = hZ->GetXaxis()->FindBin(zMax);
     //if(i>0) binx1++;
            
     TH1D *p;
     p = (hPt    ->ProjectionY("hPtProj", binx1, binx2));  
     //cout<<i<<" "<<pMin<<" "<<pMax<<" "<<binx1<<" "<<binx2<<": "<<p->GetEntries()<<endl;
     if(p->GetEntries()>50){
       fGauss->SetParameter(2,0);
       p->Fit("gauss","Q"); 
       hResPt->Fill(pMin, fGauss->GetParameter(2));fGauss->GetParError(2);
     }

     //KKK I am resetting only the sigma and its error, perhaps we should do this with all the fit parameters ???
     // SG I don't know;

     p    = (hPhi   ->ProjectionY("hPtProj2", binx1, binx2));  
     if(p->GetEntries()>50){
       fGauss->SetParameter(2,0);
       p->Fit("gauss","Q"); 
       hResPhi->Fill(pMin, fGauss->GetParameter(2)); fGauss->GetParError(2);
     }
     
     p = (hLambda->ProjectionY("hLambdaProj", binx1, binx2));  
     if(p->GetEntries()>50){
       fGauss->SetParameter(2,0);
       p->Fit("gauss","Q"); 
       hResLambda->Fill(pMin, fGauss->GetParameter(2)); fGauss->GetParError(2);
     }
     
     p     = (hY     ->ProjectionY("hYProj", binx1, binx2));  
     if(p->GetEntries()>50){
       fGauss->SetParameter(2,0);
       p->Fit("gauss","Q"); 
       hResY->Fill(pMin, fGauss->GetParameter(2)); fGauss->GetParError(2);
     }
       
     p      = (hZ     ->ProjectionY("hZProj", binz1, binz2));  
     //cout<<i<<" "<<zMin<<" "<<zMax<<" "<<binz1<<" "<<binz2<<endl;
     //cout<<p->GetEntries()<<endl;
     if(p->GetEntries()>50){
       fGauss->SetParameter(2,0);
       p->Fit("gauss","Q"); 
       hResZ->Fill(zMin, fGauss->GetParameter(2)); fGauss->GetParError(2);
     }
   }   
   
   c1->Clear(); c1->Divide(3,2);  
  
   TVirtualPad *p = c1->cd(1);
   p->SetGridy();
   hResPt->SetMarkerStyle(20);
   hResPt->Draw("P"); // KKK I haven't filled the errors for the sigma, perhaps a TGraph would make them easier to plot
  
   p = c1->cd(2);
   p->SetGridy();
   hResPhi->SetMarkerStyle(20);
   hResPhi->Draw("P");
   
   p=c1->cd(3);
   p->SetGridy();
   hResLambda->SetMarkerStyle(20);
   hResLambda->Draw("P");
   
   p=c1->cd(4);
   p->SetGridy();
   hResY->SetMarkerStyle(20);
   hResY->Draw("P");
   
   p=c1->cd(5);
   p->SetGridy();
   hResZ->SetMarkerStyle(20);
   hResZ->Draw("P");
   c1->Update();
// pad 6 stays empty
  

   //----------------- optional canvases 2 and 3 for dE/dx and clusters ----------------

   TCanvas *c2 = NULL;
   TCanvas *c3 = NULL;
   if(bDedxAndClus){     
      c2 = new TCanvas("c2","dE/dx",500,450);
      c2->cd();
      hDedx->Draw();  
        
      c3 = new TCanvas("c","clusters",500,450);
      c3->cd();   
      hClus->Draw();
      hnhit_ref->SetLineColor(kRed);
      hnhit_ref->Draw("same");
   }
   
   
   //----------------- canvas 3 for PULL VARIABLES --------------------

   TCanvas *c4 = new TCanvas("c4","pull variables",1100,900);
 
   gStyle->SetOptFit(2);
   gStyle->SetOptStat("e");

   c4->Divide(3,2);
   
   c4->cd(1);
   hpullY->Draw();
   hpullY->Fit("gaus","Q");
   
   c4->cd(2);
   hpullZ->Draw();
   hpullZ->Fit("gaus","Q");
   
   c4->cd(3);
   hpullPhi->Draw();
   hpullPhi->Fit("gaus","Q");

   c4->cd(4);
   hpullDzds->Draw();
   hpullDzds->Fit("gaus","Q");
   
   c4->cd(5);
   hpullK->Draw();
   hpullK->Fit("gaus","Q");



//-------------------- OUTPUT FILE ---------------------------------------- 
   
   TFile fc("CAtrackerEvaluation.root","RECREATE");
   c0->Write();
   c1->Write();
   if (c2) c2->Write();
   if (c3) c3->Write();
   c4->Write();
   fc.Close();
   
   gBenchmark->Stop("AliHLTTPCTrackerEvaluation");
   gBenchmark->Show("AliHLTTPCTrackerEvaluation");

   return 0;
}


//============================================================================

Int_t GoodTPCTracks(const Char_t *dir){

  Char_t fname[100];
  sprintf(fname,"%s/galice.root",dir);

  AliRunLoader *runLoader = AliRunLoader::Open(fname,"COMPARISON");
  if (!runLoader) {
     ::Error("GoodTPCTracks","Cannot start session!");
     return 1;
  }
  
  // load MC information
  
  runLoader->LoadgAlice();
  runLoader->LoadHeader();
  runLoader->LoadKinematics();
  runLoader->LoadTrackRefs();

  AliTPCLoader *tpcLoader = (AliTPCLoader*)runLoader->GetLoader("TPCLoader");
  if (tpcLoader == 0x0) {
     ::Error("GoodTPCTracks","Cannot find TPCLoader!");
     delete runLoader;
     return 2;
  }
     
  AliTPC *TPC = (AliTPC*)runLoader->GetAliRun()->GetDetector("TPC");
  Int_t tpcVersion = TPC->IsVersion(); 
  
  //cout << "TPC version "<< tpcVersion << " has been found!\n";
   
  if      (tpcVersion==1) tpcLoader->LoadRecPoints();
  else if (tpcVersion==2) tpcLoader->LoadDigits();
  else {
     ::Error("GoodTPCTracks","Wrong TPC version!");
     delete runLoader;
     return 3;
  } 

  AliCDBManager *manager = AliCDBManager::Instance();
  manager->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  manager->SetRun(0);
  
  // loading the TPC parameters
  AliTPCParamSR *digp = (AliTPCParamSR*)(AliTPCcalibDB::Instance()->GetParameters());
  if (!digp) { 
    ::Error("AliHLTTPCTrackerEvaluation.C","TPC parameters have not been found!");
    delete runLoader;
    return 4; 
  }

  Int_t nrow_up     = digp->GetNRowUp();                           //get the number of pad rows in up sector
  Int_t nrows       = digp->GetNRowLow()+nrow_up;                  //get the number of pad rows in low sector and add 
  Int_t zeroSup     = digp->GetZeroSup();
  Int_t gap         = Int_t(0.125*nrows), shift = Int_t(0.5*gap);
  Int_t good_number = Int_t(0.4*nrows);                            // will be used for selecting tracks with hits in 40% of the rows

  Int_t nEvents = runLoader->GetNumberOfEvents();
  ::Info("GoodTPCTracks","Number of events : %d\n", nEvents);  

  sprintf(fname,"%s/GoodTPCTracks.root",dir);
  TFile *file = TFile::Open(fname,"recreate");

  TClonesArray dummy("AliTrackReference",1000), *refs = &dummy;
  TTree tpcTree("tpcTree","Tree with information about the reconstructable TPC tracks");
  tpcTree.Branch("TPC",&refs);

  
  
  //-------------- Loop over generated events ------------------------------
  
  for(Int_t iEvent=0; iEvent<nEvents; iEvent++){

      Int_t nt = 0;
      refs->Clear();

      Int_t i;

      runLoader->GetEvent(iEvent);  
      file->cd();

      Int_t nParticles = runLoader->GetHeader()->GetNtrack();
      //cout << "Event " << iEvent << ", Number of particles: " << nParticles << endl;

      // ------- Set the selection criteria for the MC sample that will be used for the efficiency -------
      
      Int_t *good = new Int_t[nParticles]; for(i=0; i<nParticles; i++) good[i] = 0;
      
      switch (tpcVersion){
      case 1: // not used any more
         {
           AliTPCClustersArray *pca = new AliTPCClustersArray, &ca=*pca;
           ca.Setup(digp);
           ca.SetClusterType("AliTPCcluster");
           ca.ConnectTree(tpcLoader->TreeR());
           Int_t nrows=Int_t(ca.GetTree()->GetEntries());
           for (Int_t n=0; n<nrows; n++) {
    	     AliSegmentID *s=ca.LoadEntry(n);
    	     Int_t sec,row;
    	     digp->AdjustSectorRow(s->GetID(),sec,row);
    	     AliTPCClustersRow &clrow = *ca.GetRow(sec,row);
    	     Int_t ncl=clrow.GetArray()->GetEntriesFast();
    	     while (ncl--) {
    	   	 AliTPCcluster *c=(AliTPCcluster*)clrow[ncl];
    	   	 Int_t label=c->GetLabel(0);
    	   	 if (label<0) continue; //noise cluster
    	   	 label=TMath::Abs(label);

    	   	 if (sec>=digp->GetNInnerSector())
    	   	    if (row==nrow_up-1) good[label]|=0x4000;
    	   	 if (sec>=digp->GetNInnerSector())
    	   	    if (row==nrow_up-1-gap) good[label]|=0x1000;

    	   	 if (sec>=digp->GetNInnerSector())
    	   	    if (row==nrow_up-1-shift) good[label]|=0x2000;
    	   	 if (sec>=digp->GetNInnerSector())
    	   	    if (row==nrow_up-1-gap-shift) good[label]|=0x800;

    	   	 good[label]++;
    	     }
    	     ca.ClearRow(sec,row);
           }
           delete pca;
           }
         break;
     
      case 2:
         {
           TTree *TD = tpcLoader->TreeD(); // get the digits tree
       
           AliSimDigits da, *digits = &da;
           TD->GetBranch("Segment")->SetAddress(&digits);

           Int_t *count = new Int_t[nParticles];
           Int_t i;
           for (i=0; i<nParticles; i++) count[i] = 0;

           Int_t sectors_by_rows = (Int_t)TD->GetEntries();
           
	   for(i=0; i<sectors_by_rows; i++){

    	     if (!TD->GetEvent(i)) continue;
    	     Int_t sec,row;
    	     
	     digp->AdjustSectorRow(digits->GetID(),sec,row);
             if(digits->First()){
               do { 
    	   	 Int_t it    = digits->CurrentRow();
	   	 Int_t ip    = digits->CurrentColumn(); 	 
    	   	 Short_t dig = digits->GetDigit(it,ip);
    	   	 Int_t idx0  = digits->GetTrackID(it,ip,0); 
    	   	 Int_t idx1  = digits->GetTrackID(it,ip,1);
    	   	 Int_t idx2  = digits->GetTrackID(it,ip,2);
    	   	 if(idx0>=0 && dig>=zeroSup && idx0<nParticles) count[idx0]+=1;
    	   	 if(idx1>=0 && dig>=zeroSup && idx1<nParticles) count[idx1]+=1;
    	   	 if(idx2>=0 && dig>=zeroSup && idx2<nParticles) count[idx2]+=1;
               } while (digits->Next());
             }
            
	     for(Int_t j=0; j<nParticles; j++){
    	   	
		 if(count[j]>1){
    	   	    if(sec>=digp->GetNInnerSector())
          	    if(row==nrow_up-1)  	     good[j]|=0x4000;
    	   	    if(sec>=digp->GetNInnerSector())
          	    if(row==nrow_up-1-gap)	     good[j]|=0x1000;

    	   	    if(sec>=digp->GetNInnerSector())
          	    if(row==nrow_up-1-shift)         good[j]|=0x2000;
    	   	    if(sec>=digp->GetNInnerSector())
          	    if(row==nrow_up-1-gap-shift)     good[j]|=0x800;

		    good[j]++;
    	   	 }
    	   	 count[j] = 0;
    	     } // end of loop over particles
           } // end of loop over sectors_by_rows
           delete [] count;
         } // end of case 2
         break;
      } // end of switch




    // ---------------- Select sample of MC particles that will be used for forming the efficiency ------------------
   
    AliStack *stack = runLoader->Stack();

    for(i=0; i<nParticles; i++){ // loop over stack

        if ((good[i]&0x5000) != 0x5000)//SG!!!
        if ((good[i]&0x2800) != 0x2800)     continue;
        if ((good[i]&0x7FF ) < good_number) continue;
	// N TPC rows with digits >= 64 => common checks of cluster finder & tracker
        // There are digits in rows (159&&139) || (149&&129)


        TParticle *part = (TParticle*)stack->Particle(i);
        if(part == 0x0){
	   cerr << "Cannot get particle "<< i << endl;
      	   continue;
        }
       
        if(part->Pt()<0.100) continue;                        // skip particles with pt below 0.1 GeV/c
        if(TMath::Abs(part->Pz()/part->Pt())>0.999) continue; // reject tracks outside the TPC acceptance, below 45 degrees 

        Double_t vx = part->Vx();
	Double_t vy = part->Vy();
	Double_t vz = part->Vz();
        //if (TMath::Sqrt(vx*vx+vy*vy)>3.5) continue; // vertex cuts
        //if (TMath::Abs(vz) > 50.) continue;

        AliTrackReference *ref = new((*refs)[nt])AliTrackReference();

        ref->SetLabel(i);
        Int_t pdg = part->GetPdgCode();
        ref->SetLength(pdg);  //This will the particle's PDG !
        ref->SetMomentum(0.,0.,0.);  
        ref->SetPosition(0.,0.,0.);

        nt++;
    } // end of loop over stack

    //----------- Check whether there is also information at the entrance of the TPC
    
    TTree   *TR     = runLoader->TreeTR();
    TBranch *branch = TR->GetBranch("TrackReferences");
    if(branch==0){
       ::Error("GoodTPCTracks","No track references!");
       delete runLoader;
       return 5;
    }
    
    TClonesArray tpcdummy("AliTrackReference",1000), *tpcRefs = &tpcdummy;
    branch->SetAddress(&tpcRefs);

    for(Int_t r=0; r<(Int_t)TR->GetEntries(); r++){ 

  	//cerr<<r<<' '<<(Int_t)TR->GetEntries()<<'\r';
  	TR->GetEvent(r);

	if(!tpcRefs->GetEntriesFast()) continue;
  	
	AliTrackReference *tpcRef = 0x0; 	        
	for(Int_t iref=0; iref<tpcRefs->GetEntriesFast(); ++iref){

            tpcRef = (AliTrackReference*)tpcRefs->UncheckedAt(iref);
            if(tpcRef->DetectorId() == AliTrackReference::kTPC) break;
            tpcRef = 0x0;
        }

        if(!tpcRef) continue;

  	Int_t j;
        AliTrackReference *ref = 0;
  	
	for(j=0; j<nt; j++){

          ref = (AliTrackReference*)refs->UncheckedAt(j);
  	  if(ref->Label()==tpcRef->Label()) break;
  	  ref = 0;
        }  
  	if(ref==0) continue;
        
  	ref->SetMomentum(tpcRef->Px(),tpcRef->Py(),tpcRef->Pz());
  	ref->SetPosition(tpcRef->X(), tpcRef->Y(), tpcRef->Z());

        tpcRefs->Clear();
    }

    tpcTree.Fill();

    delete [] good;
  } //---------- end of the loop over generated events
  
  tpcTree.Write();
  file->Close();

  delete runLoader;
  return 0;
}
#endif
