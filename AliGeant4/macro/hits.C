// $Id$

#include <iostream.h>

void hits()
{
  hits("FMD","1");
  hits("ITS","5");
  hits("MUON","0");
  hits("PMD","0");
  hits("START","0");
  hits("TOF","0");
  hits("TPC","1");
  hits("PHOS","1");
}  

void hits(const TString& detName, const TString& detVersion)
{
  // labels
  TString g3 = "G3: ";
  TString g4 = "G4: ";

  // construct file names
  TString top = getenv("ALICE_ROOT");
  TString f1NameEnd = "test10.root";
  TString f2NameEnd = "test20.root";
  TString g3File1Name = top + "/test/" + f1NameEnd;
  TString g3File2Name = top + "/test/" + f2NameEnd;
  TString g4File1Name = top + "/AliGeant4/test/" + f1NameEnd;
  TString g4File2Name = top + "/AliGeant4/test/" + f2NameEnd;
  
  cout << "Files: " << endl;
  cout << g3File1Name << endl;
  cout << g3File2Name << endl;
  cout << g4File1Name << endl;
  cout << g4File2Name << endl;

  // link shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
  cout << "AliRoot libraries were loaded." << endl;
      
  // set histogram ranges
  Int_t nbin;
  Int_t g3xmin; Int_t g4xmin; Int_t g3xmax; Int_t g4xmax;
  Int_t g3ymin; Int_t g4ymin; Int_t g3ymax; Int_t g4ymax;
  Int_t g3zmin; Int_t g4zmin; Int_t g3zmax; Int_t g4zmax;
  SetHistogramRanges(detName, nbin, 
                     g3xmin, g3xmax, g4xmin, g4xmax,
                     g3ymin, g3ymax, g4ymin, g4ymax,
                     g3zmin, g3zmax, g4zmin, g4zmax);


  // create histograms
  TH1F* g3x  = new TH1F("g3x",  g3 + detName + " hits per x", nbin, g3xmin, g3xmax);
  TH1F* g3xh = new TH1F("g3xh", g3 + detName + " hits per x", nbin, g3xmin, g3xmax);
  TH1F* g4x  = new TH1F("g4x",  g4 + detName + " hits per x", nbin, g4xmin, g4xmax);
  TH1F* g4xh = new TH1F("g4xh", g4 + detName + " hits per x", nbin, g4xmin, g4xmax);

  TH1F* g3z  = new TH1F("g3z",  g3 + detName + " hits per z", nbin, g3zmin, g3zmax);
  TH1F* g3zh = new TH1F("g3zh", g3 + detName + " hits per z", nbin, g3zmin, g3zmax);
  TH1F* g4z  = new TH1F("g4z",  g4 + detName + " hits per z", nbin, g4zmin, g4zmax);
  TH1F* g4zh = new TH1F("g4zh", g4 + detName + " hits per z", nbin, g4zmin, g4zmax);

  TH1F* g3y  = 0;
  TH1F* g3yh = 0;
  TH1F* g4y  = 0;
  TH1F* g4yh = 0;
  if (detName == "PHOS") {
    TH1F* g3y  = new TH1F("g3y",  g3 + detName + " hits per y", nbin, g3ymin, g3ymax);
    TH1F* g3yh = new TH1F("g3yh", g3 + detName + " hits per y", nbin, g3ymin, g3ymax);
    TH1F* g4y  = new TH1F("g4y",  g4 + detName + " hits per y", nbin, g4ymin, g4ymax);
    TH1F* g4yh = new TH1F("g4yh", g4 + detName + " hits per y", nbin, g4ymin, g4ymax);
  }

  cout << "Histograms were created." << endl;

  // fill histograms
  AliDetector* detector;

  //detector = LoadDetector(g3File1Name, detName, g3);
  //FillHistogram(detector, g3, g3x, g3y, g3z); 

  //detector = LoadDetector(g3File2Name, detName,  g3);
  FillHistogram(detector, g3, g3xh, g3yh, g3zh); 

  detector = LoadDetector(g4File1Name, detName, g4);
  FillHistogram(detector, g4, g4x, g4y, g4z); 

  //detector = LoadDetector(g4File2Name, detName, g4);
  FillHistogram(detector, g4, g4xh, g4yh, g4zh); 

  // compose picture name
  TString gifNameBase =  "hits" + detName + "v" + detVersion;
  TString title = detName + " hits";

  // draw histohrams
  DrawHistograms(title, gifNameBase + "_x.gif", g3x, g4x, g3xh, g4xh);
  DrawHistograms(title, gifNameBase + "_y.gif", g3y, g4y, g3yh, g4yh);
  DrawHistograms(title, gifNameBase + "_z.gif", g3z, g4z, g3zh, g4zh);
}

void SetHistogramRanges(TString& detName, Int_t& nbin, 
                        Int_t& g3xmin, Int_t& g3xmax,
                        Int_t& g4xmin, Int_t& g4xmax,
                        Int_t& g3ymin, Int_t& g3ymax, 
                        Int_t& g4ymin, Int_t& g4ymax, 
                        Int_t& g3zmin, Int_t& g3zmax, 
                        Int_t& g4zmin, Int_t& g4zmax)
{			 
  nbin = 200;

  g3xmin = 0; g4xmin = g3xmin;
  g3xmax = 1; g4xmax = g3xmax;
  g3ymin = 0; g4ymin = g3xmin;
  g3ymax = 1; g4ymax = g3xmax;
  g3zmin = 0; g4zmin = g3zmin;
  g3zmax = 1; g4zmax = g3zmax;

  if (detName == "FMD") {
    g3xmin = -70; g4xmin = g3xmin;
    g3xmax =  70; g4xmax = g3xmax;
    g3zmin = -300; g4zmin = g3zmin;
    g3zmax =  300; g4zmax = g3zmax;
  }
  else if (detName == "ITS") {
    g3xmin = -60; g4xmin = g3xmin;
    g3xmax =  60; g4xmax = g3xmax;
    g3zmin = -60; g4zmin = g3zmin;
    g3zmax =  60; g4zmax = g3zmax;
  }  
  else if (detName == "MUON") {
    g3xmin = -450; g4xmin = g3xmin;
    g3xmax =  450; g4xmax = g3xmax;
    g3zmin =  400; g4zmin = g3zmin;
    g3zmax = 1800; g4zmax = g3zmax;
  }
  else if (detName == "PHOS") {
    g3xmin = -400; g4xmin = g3xmin;
    g3xmax =  400; g4xmax = g3xmax;
    g3ymin = -500; g4ymin = g3ymin;
    g3ymax = -250; g4ymax = g3ymax;
    g3zmin =  -80; g4zmin = g3zmin;
    g3zmax =   80; g4zmax = g3zmax;
  }
  else if (detName == "PMD") {
    g3xmin = -200; g4xmin = g3xmin*10;
    g3xmax =  200; g4xmax = g3xmax*10;
    g3zmin = -582; g4zmin = g3zmin*10;
    g3zmax = -578; g4zmax = g3zmax*10;
  }
  else if (detName == "START") {
    g3xmin = -10; g4xmin = g3xmin;
    g3xmax =  10; g4xmax = g3xmax;
    g3zmin = -80; g4zmin = g3zmin;
    g3zmax =  80; g4zmax = g3zmax;
  }
  else if (detName == "TOF") {
    g3xmin = -400; g4xmin = g3xmin;
    g3xmax =  400; g4xmax = g3xmax;
    g3zmin = -400; g4zmin = g3zmin;
    g3zmax =  400; g4zmax = g3zmax;
  }
  else if (detName == "TPC") {
    g3xmin = -300; g4xmin = g3xmin;
    g3xmax =  300; g4xmax = g3xmax;
    g3zmin = -300; g4zmin = g3zmin;
    g3zmax =  300; g4zmax = g3zmax;
  }
}  

AliDetector* LoadDetector(TString& fileName, TString& detName, TString& label)
{ 
  // connect the Root files
  TFile* file = new TFile(fileName);
  if (!file) {
    cerr << "Root file was not found." << endl;
    return 0;
  }

  // get AliRun object from file 
  if (gAlice) delete gAlice;
  gAlice = 0;
  gAlice = (AliRun*)file->Get("gAlice");
  if (!gAlice) {
    cerr << label << "AliRun object not found in a file." << endl;
    return 0;
  }
      
  // import the hits trees 
  Int_t nparticles = gAlice->GetEvent(0);
  if (nparticles <= 0) {
    cerr << label << "No particles in a file." << endl;
    return 0;
  }  
  cout << label << "got nparticles = " << nparticles << endl;

  // get detector
  AliDetector* detector = gAlice->GetDetector(detName);
  if (!detector) {
    cerr << label << "No detector " << detName << " in a file." << endl;
    return 0;
  }  
  
  return detector;
}  

void FillHistogram(AliDetector* detector, TString& label, 
                   TH1F* hx, TH1F* hy, TH1F* hz) 
{ 
  Int_t nofHits = 0;
  // get number of primary tracks
  Int_t ntracks = gAlice->TreeH()->GetEntries();
  cout << label << "got ntracks = " << ntracks << endl;

  // loop on tracks in the hits container
  for (Int_t i=0; i<ntracks; i++)
    // loop on hits  
    for(AliHit* hit= detector->FirstHit(i); hit; hit=detector->NextHit()) {
      Float_t x = hit->X();
      Float_t y = hit->Y();
      Float_t z = hit->Z();
      if (hx) hx->Fill(x);
      if (hy) hy->Fill(y);
      if (hz) hz->Fill(z);
      nofHits++;
    }

  cout << label << "filled " << nofHits << " hits" << endl;
}  

void DrawHistograms(TString& title, TString& gifName,
                    TH1F* h1, TH1F* h2, TH1F* h3, TH1F* h4)
{
// Create canvas, set the view range, show histograms.
// ---

  if (h1 && h2 && h3 && h4) {
    // create canvas
    //TCanvas* canvas = new TCanvas("c1",detName + " hits", 400, 10, 800, 600);
    TCanvas* canvas = new TCanvas("c1", title, 400, 10, 800, 600);
    canvas->Divide(2,2);
    
    // draw histograms
    canvas->cd(1); h1->Draw();
    canvas->cd(2); h2->Draw(); 
    canvas->cd(3); h3->Draw(); 
    canvas->cd(4); h4->Draw();
    
    // save gif
    //canvas->SaveAs(gifNameBase + "_x.gif"); 
    canvas->SaveAs(gifName); 
  }  
}
