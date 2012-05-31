#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TRD/AliTRDcluster.h>
#include <AliRunLoader.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#endif

//#define EVENT_BY_EVENT      // Show amplitude event-by-event
//#define EVENT_BY_EVENT2     // Show tracks event-by-event
//#define PRINT_OTHER         // Display various additional histograms

// Number of noisy detectors
const Int_t numNoisyDets = 10;

// Holds the detector number of detectors that produce too much noise and
// are therefore excluded from the calculation
const Int_t exDets[numNoisyDets] = {25, 511, 515, 538, 523, 249, 294, 27, 24, 247};   

Bool_t isBadCluster(Int_t det, Int_t Q, Int_t PadCol, Int_t PadRow)
{  
  if (Q < 20) return 1;

  // Exclude bad cols
  if (PadCol == 142 || PadCol == 136 || PadCol == 66 || PadCol == 65 || PadCol <= 2) return 1;
  //if (PadCol == 142 || PadCol == 66 || PadCol == 65 || PadCol <= 2) return 1;
  //if (PadCol == 142 || PadCol <= 2) return 1;
  
  // Exclude noisy detectors
  for (Int_t indDet = 0; indDet < numNoisyDets; indDet++)
  {
    if (det == exDets[indDet])
    {
      // Exclude several regions of bad detectors
      if ((PadCol >= 134 && PadRow <= 5) || (PadCol == 141 && PadRow >= 10)) return 1;
      //if (PadRow == 15 || (PadCol >= 134 && PadRow <= 11) || (PadCol >= 119 && PadCol <= 123 && PadRow == 11)) return 1; 

      break;
    }
  }
  
  // Cluster seems to be clean
  return 0;
}


//=====================================

void AliTRDanalyseRecPoints (const char *filename="galice.root") 
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(1);

  // open run loader and load gAlice, kinematics and header
  AliRunLoader* runLoader = AliRunLoader::Open(filename);
  if (!runLoader) 
  {
    printf("Error: readKine:\nGetting run loader from file \"%s\" failed", filename);
    return;
  }

  runLoader->LoadHeader();
  runLoader->LoadRecPoints("TRD");
  TObjArray *module = new TObjArray();

  runLoader->CdGAFile();
  

  // Flag indicating whether the signal comes from a noisy detector
  Bool_t fromNoisyDet = 0;

  // Counts the number of clusters per event
  Int_t clusters = 0;

  // Threshold for the amplitude
  const Int_t kQThreshold = 30;

  // Selection criterions for clusters
  const Int_t kThresholdLower = 30;
  const Int_t kThresholdUpper = 500;

  // Number of events passing the criterion
  Int_t numPassed = 0;

  // hist
#ifdef PRINT_OTHER
  TH1D *mNCls = new TH1D("NCls", ";number of clusters", 400, 0, 400);
  TH1D *mResY = new TH1D("resY", ";sigma Y (cm)", 100, 0, 1);
  TH1D *mResZ = new TH1D("resZ", ";sigma Z (cm)", 100, 1, 3);

  TH1D *mNPads = new TH1D("NPads", ";number of pads", 10, -0.5, 9.5);
  TH1D *mCenter = new TH1D("center", ";center", 100, -1, 1);
  TH1D *mTime = new TH1D("timeBin", ";time bin", 25, -0.5, 24.5);
  TH1D *mPlane = new TH1D("plane", ";plane", 6, -0.5, 5.5);

  TH1D *mX = new TH1D("x", "; X (cm)", 100, 0, 3.5);
  TH1D *mY = new TH1D("y", "; y / chamber width (%)", 100, -1, 1);
  TH1D *mZ = new TH1D("z", "; Z (cm)", 200, -400, 400);

  TH2D *mPlaneY = new TH2D("planeY", ";plane;Y (cm);#Clusters", 6, -0.5, 5.5, 100, -60, 60);
  TH2D *mYX = new TH2D("yx", ";Y (cm);X (cm);#Clusters", 100, -60, 60, 120, 280, 400);
#endif

  TH1D *mChrg = new TH1D("chrg", ";Charge(no pad-filter)", 200, 0, 200);
  TH1D *mChrg2 = new TH1D("chrg2", ";Charge(pads<9)", 200, 0, 200);
  TH1D *mChrg3 = new TH1D("chrg3", ";Charge(2<pads<9)", 200, 0, 200);
  
  TH2D *mQvsPads = new TH2D("QvsPads", ";Number of pads;Q;#Clusters", 15, -0.5, 14.5, 100, 0, 1000);

  // There are 540 detectors
  TH2D *mQvsDetector = new TH2D("QvsDet", ";Number of detector;Q;#Clusters", 540, -0.5, 539.5, 100, 0, 1000);

  // Supermoduleid=int(detectorNum / 30)
  TH2D *mQvsSModId = new TH2D("QvsSModId", ";Supermodule id;Q;#Clusters", 18, -0.5, 17.5, 100, 0, 1000);
  // Histograms for the already installed supermodules
  TH1D *mSMod0 = new TH1D("SMod0", ";Charge", 200, 0, 200);
  TH1D *mSMod8 = new TH1D("SMod8", ";Charge", 200, 0, 200);
  TH1D *mSMod9 = new TH1D("SMod9", ";Charge", 200, 0, 200);
  TH1D *mSMod17 = new TH1D("SMod17", ";Charge", 200, 0, 200);

  // Index + 2 corresponds to the number of activated pads
  TH1D *mPads[14];
  char name[8];
  for (Int_t ind = 0; ind < 14; ind++)
  {
    sprintf(name,"%dPads",ind + 2);
    mPads[ind] = new TH1D(name, ";Charge", 200, 0, 200);
  }

  // Amplitude distribution will be displayed event by event for the first 100 events
  TH1D *mChrgEBE = new TH1D("chrgEBE", ";Charge", 200, 0, 200);

  // Keep track of the number of clusters per event
  TH1D *mClusters = new TH1D("Clusters", ";Number of clusters per event", 300, -0.5, 299.5);

  // Keep track of the number of clusters per chamber with Q < kQThreshold
  char clusTitle[55];
  sprintf(clusTitle,";Number of clusters with Q < %d vs detector", kQThreshold);
  TH1D *mClusCham = new TH1D("ClusCham", clusTitle, 540, -0.5, 539.5);
  // Same but this time noisy detectors excluded
  sprintf(clusTitle,";Number of clusters with Q < %d vs detector", kQThreshold);
  TH1D *mClusChamEx = new TH1D("ClusChamEx", clusTitle, 540, -0.5, 539.5);

  TH2D *mColRow = new TH2D("ColRow", ";Column;Row;#Clusters", 144, -0.5, 143.5, 18, -0.5, 17.5);
  TH2D *mColRowFiltered = new TH2D("ColRowFiltered", ";Column;Row;#Clusters", 144, -0.5, 143.5, 18, -0.5, 17.5);

  TH2D *mColTime = new TH2D("ColTime", ";TimeBin;Column", 18, -0.5, 17.5, 144, -0.5, 143.5);
  TH2D *mColTimeFiltered = new TH2D("ColTimeFiltered", ";TimeBin;Column;#Clusters", 18, -0.5, 17.5, 144, -0.5, 143.5);

  // YX histograms for the already installed supermodules
  TH2D *mTrSMod0 = new TH2D("TrSMod0", ";Col;time+layer*20;#Clusters", 144, -0.5, 143.5, 120, -0.5, 119.5);
  TH2D *mTrSMod8 = new TH2D("TrSMod8", ";Col;time+layer*20;#Clusters", 144, -0.5, 143.5, 120, -0.5, 119.5);
  TH2D *mTrSMod9 = new TH2D("TrSMod9", ";Col;time+layer*20;#Clusters", 144, -0.5, 143.5, 120, -0.5, 119.5);
  TH2D *mTrSMod17 = new TH2D("TrSMod17", ";Col;time+layer*20;#Clusters", 144, -0.5, 143.5, 120, -0.5, 119.5);


#ifdef EVENT_BY_EVENT  
  TCanvas* eventc = new TCanvas("eventc","Amplitude distribution for event 0",700,500);
  char title[40];
#endif

#ifdef EVENT_BY_EVENT2
  TCanvas* sTrc = new TCanvas("sTrc","Track -  Supermodules",800,800);
  sTrc->SetLogz();
  sTrc->Divide(2,2);
  char title2[40];
#endif

 
  int nEvents = runLoader->GetNumberOfEvents();
  
  AliTRDcluster *cls = 0;

  for(Int_t ev = 0; ev < nEvents - 1; ev++) 
  {
    clusters = 0;    

    TTree *tree = runLoader->GetTreeR("TRD", 0);
    tree->SetBranchAddress("TRDcluster", &module);

    int N = tree->GetEntries();

    // Check number of clusters
    for(Int_t ind = 0; ind < N; ind++)
    {
      tree->GetEntry(ind);
      Int_t m = module->GetEntries();
  
      for (Int_t j = 0; j < m; j++)
      {
        if (cls != 0) delete cls;
        cls = (AliTRDcluster*)module->At(j);
        if (!isBadCluster(cls->GetDetector(), (Int_t)cls->GetQ(), cls->GetPadCol(), cls->GetPadRow())) clusters++;
        else
	      {
          // Uncomment the following lines AND comment the next if-clause to display the histograms for the 
          // EXCLUDED clusters

	        //mQvsPads->Fill(cls->GetNPads(),cls->GetQ());
	        //mQvsDetector->Fill(cls->GetDetector(),cls->GetQ());
	        //mQvsSModId->Fill((Int_t)(cls->GetDetector() / 30.),cls->GetQ());	
	        //mColTime->Fill(cls->GetLocalTimeBin(), cls->GetPadCol());
          //mColRow->Fill(cls->GetPadCol(), cls->GetPadRow());    
	        //mClusCham->Fill(cls->GetDetector());
          //mYX->Fill(cls->GetY(), cls->GetX());
        }
      }  
    }

    // Process event, if number of clusters is in the correct range
    if (clusters >= kThresholdLower && clusters <= kThresholdUpper)
    {
      clusters = 0;
      numPassed++;

      for (Int_t i = 0; i < N; i++) 
      {
	      tree->GetEntry(i);
	      int m = module->GetEntries();

#ifdef PRINT_OTHER
	      mNCls->Fill(m);
#endif

	      for(Int_t j = 0; j < m; j++) 
	      {
          fromNoisyDet = 0;
	        clusters++;

          if (cls != 0) delete cls;
	        cls = (AliTRDcluster*)module->At(j);
	      
	        // Filter the clusters of the events, too
          if (cls->GetQ() < 20) continue;
           
	        mColTime->Fill(cls->GetLocalTimeBin(), cls->GetPadCol());

          if (cls->GetQ() <= kQThreshold) 
          {
            mClusCham->Fill(cls->GetDetector());
	        }
          
          // Exclude noisy detectors
          for (Int_t indDet = 0; indDet < numNoisyDets; indDet++)
	        {
	          if (cls->GetDetector() == exDets[indDet])
	          {
          	  fromNoisyDet = 1;
	            break;
	          }
	        }

          // Exclude I
          if (cls->GetPadCol() == 142 || cls->GetPadCol() == 136 || cls->GetPadCol() == 66 || cls->GetPadCol() == 65 || cls->GetPadCol() <= 2) continue;
          //if (cls->GetPadCol() == 142 || cls->GetPadCol() == 66 || cls->GetPadCol() == 65 || cls->GetPadCol() <= 2) continue;
          //if (cls->GetPadCol() == 142 || cls->GetPadCol() <= 2) continue;

	        if (fromNoisyDet) 
	        {
	          mColRow->Fill(cls->GetPadCol(), cls->GetPadRow());    
 
	          // Exclude II
	          //if (cls->GetPadRow() == 15 || (cls->GetPadCol() >= 134 && cls->GetPadRow() <= 11) || (cls->GetPadCol() >= 119 && cls->GetPadCol() <= 123 && cls->GetPadRow() == 11)) continue; 
            
            if ((cls->GetPadCol() >= 134 && cls->GetPadRow() <= 5) || (cls->GetPadCol() == 141 && cls->GetPadRow() >= 10)) continue; 
            
            mColRowFiltered->Fill(cls->GetPadCol(), cls->GetPadRow());
	        }

          mColTimeFiltered->Fill(cls->GetLocalTimeBin(), cls->GetPadCol());

          if (cls->GetQ() <= kQThreshold) 
          {
            mClusChamEx->Fill(cls->GetDetector());
	        }
	        
#ifdef PRINT_OTHER
	        mX->Fill(cls->GetX());
	        mZ->Fill(cls->GetZ());
          mYX->Fill(cls->GetY(), cls->GetX());
	        mResY->Fill(TMath::Sqrt(cls->GetSigmaY2()));
	        mResZ->Fill(TMath::Sqrt(cls->GetSigmaZ2()));
          mNPads->Fill(cls->GetNPads());
          mCenter->Fill(cls->GetCenter());
	        mTime->Fill(cls->GetLocalTimeBin());
#endif

	        if (cls->GetNPads() > 2 && cls->GetNPads() < 9) mChrg3->Fill(cls->GetQ());
          if (cls->GetNPads() < 9) mChrg2->Fill(cls->GetQ());
          mChrg->Fill(cls->GetQ());	         
	        
	        mQvsPads->Fill(cls->GetNPads(),cls->GetQ());
	        mQvsDetector->Fill(cls->GetDetector(),cls->GetQ());
	        mQvsSModId->Fill((Int_t)(cls->GetDetector() / 30.),cls->GetQ());	

          Int_t superModule = (Int_t)(cls->GetDetector() / 30.);

          // Position in supermodule = cls->GetDetector() - superModule * 30
          // There are 6 layers per supermodule:
          Int_t layer = (cls->GetDetector() - superModule * 30) % 6;
	        switch (superModule)
	        {
	          case 0:
	            mSMod0->Fill(cls->GetQ());
                    mTrSMod0->Fill(cls->GetPadCol(), cls->GetLocalTimeBin() + layer * 20);
	            break;
	          case 8:
	            mSMod8->Fill(cls->GetQ());
                    mTrSMod8->Fill(cls->GetPadCol(), cls->GetLocalTimeBin() + layer * 20);
	            break;
	          case 9:
	            mSMod9->Fill(cls->GetQ());
                    mTrSMod9->Fill(cls->GetPadCol(), cls->GetLocalTimeBin() + layer * 20);
	            break;
	          case 17:
	            mSMod17->Fill(cls->GetQ());
                    mTrSMod17->Fill(cls->GetPadCol(), cls->GetLocalTimeBin() + layer * 20);
	            break;
	          default:
	            break;
	        }

	        // Only consider events with maximal 15 pads here
	        if (cls->GetNPads() - 2 < 14 && cls->GetNPads() - 2 >= 0)
	        {
	          mPads[cls->GetNPads() - 2]->Fill(cls->GetQ());
	        }

	        mChrgEBE->Fill(cls->GetQ());
	      }
      }
      

#ifdef EVENT_BY_EVENT  
      // Event by event display for the first 100 events
      if (numPassed < 100)
      {
        sprintf(title,"Charge distribution for event %d", ev);
        eventc->SetTitle(title);
        mChrgEBE->Draw();
        eventc->Update();
        gSystem->Sleep(400);
        mChrgEBE->Reset();
      }
#endif

#ifdef EVENT_BY_EVENT2
      // Event by event display 
      sprintf(title2,"Supermodule tracks for event %d", ev);
      sTrc->SetTitle(title2);
      sTrc->cd(1);
      mTrSMod0->Draw("colz");
      sTrc->cd(2);
      mTrSMod8->Draw("colz");
      sTrc->cd(3);
      mTrSMod9->Draw("colz");
      sTrc->cd(4);
      mTrSMod17->Draw("colz");
      sTrc->Update();
      gSystem->Sleep(400);
      mTrSMod0->Reset();
      mTrSMod8->Reset();
      mTrSMod9->Reset();
      mTrSMod17->Reset();
#endif


      mClusters->Fill(clusters);    
    }
     
    runLoader->GetNextEvent();
  }

#ifdef PRINT_OTHER
  new TCanvas();
  gPad->SetLogy();
  mNCls->Draw();

  new TCanvas();
  gPad->SetLogy();
  mResY->Draw();

  new TCanvas();
  gPad->SetLogy();
  mResZ->Draw();
#endif

  new TCanvas("chrgc", "Charge(no pad-filter)", 700, 500);
  mChrg->Draw();
  // Try landau-fit
  TF1* lFit = new TF1("lFit", "landau", 20, 200);
  cout << "\n\nFit (no pad-filter):\n";
  mChrg->Fit(lFit, "", "", 20, 200);

  new TCanvas("chrg2c", "Charge(pads<9)", 700, 500);
  mChrg2->Draw();
  // Try landau-fit
  TF1* lFit2 = new TF1("lFit2", "landau", 20, 200);
  cout << "\n\nFit (pads<9):\n";
  mChrg2->Fit(lFit2, "", "", 20, 200);

  new TCanvas("chrg3c", "Charge(2<pads<9)", 700, 500);
  mChrg3->Draw();
  // Try landau-fit
  TF1* lFit3 = new TF1("lFit3", "landau", 20, 200);
  cout << "\n\nFit (2<pads<9):\n";
  mChrg3->Fit(lFit3, "", "", 20, 200);

#ifdef PRINT_OTHER
  new TCanvas();
  mNPads->Draw();

  new TCanvas();
  mTime->Draw();

  new TCanvas();
  mX->Draw();

  new TCanvas();
  mY->Draw();

  new TCanvas();
  mZ->Draw();

  new TCanvas();
  mCenter->Draw();

  new TCanvas();
  mYX->Draw("colz");

  new TCanvas();
  mPlane->Draw();

  new TCanvas();
  mPlaneY->Draw("colz");
#endif

  new TCanvas();
  gPad->SetLogz();
  mQvsPads->Draw("colz");  

  new TCanvas();
  gPad->SetLogz();
  mQvsDetector->Draw("colz");
  TLine* tl = 0;
  for (Int_t num = 30; num < 540; num += 30)
  {
    // Draw vertical lines every 30 bins to sketch the layer modules
    tl = new TLine(num - 0.5, 0, num - 0.5, 1000);
    tl->Draw();
  }

  new TCanvas();
  gPad->SetLogz();
  mQvsSModId->Draw("colz");

  TCanvas* sc = new TCanvas("sc","Supermodules", 800, 800);
  sc->SetLogz();
  sc->Divide(2,2);
  sc->cd(1);
  mSMod0->Draw();
  sc->cd(2);
  mSMod8->Draw();
  sc->cd(3);
  mSMod9->Draw();
  sc->cd(4);
  mSMod17->Draw();
 

  TCanvas* padc = new TCanvas("padc","Pads", 1000, 1000);
  padc->Divide(4,4);
  for (Int_t ind = 0; ind < 14; ind++)
  {
    padc->cd(ind + 1);
    mPads[ind]->Draw();
  }

  new TCanvas();
  mClusters->Draw();

  TCanvas* clusChamc = new TCanvas("clusChamc", "Clusters against detector with threshold", 800, 800);
  clusChamc->Divide(1,2);
  clusChamc->cd(1);
  mClusCham->Draw();
  clusChamc->cd(2);
  mClusChamEx->Draw();

  TCanvas* noisyc = new TCanvas("noisyc", "Col/Row of noisy detectors");
  noisyc->Divide(1,2);
  noisyc->cd(1);
  gPad->SetLogz();
  mColRow->Draw("colz");
  noisyc->cd(2);
  gPad->SetLogz();
  mColRowFiltered->Draw("colz");

  TCanvas* ColTimec = new TCanvas("ColTimec", "Col/Time");
  ColTimec->Divide(1,2);
  ColTimec->cd(1);
  gPad->SetLogz();
  mColTime->Draw("colz");
  ColTimec->cd(2);
  gPad->SetLogz();
  mColTimeFiltered->Draw("colz");

  cout << numPassed << " out of " << nEvents << " events passed the criterion (number of 'good' clusters >= ";
  cout << kThresholdLower << " && clusters <= " << kThresholdUpper << ").\n";
}
