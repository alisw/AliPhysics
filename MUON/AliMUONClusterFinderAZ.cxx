#include "AliMUONClusterFinderAZ.h"

#include <stdlib.h>
#include <fcntl.h>
#include <Riostream.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TTree.h>
#include <TH2.h>
#include <TView.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TMatrixD.h>

#include "AliHeader.h"
#include "AliRun.h"
#include "AliMUON.h"
#include "AliMUONChamber.h"
#include "AliMUONDigit.h"
#include "AliMUONHit.h"
#include "AliMUONChamber.h"
#include "AliMUONRawCluster.h"
#include "AliMUONClusterInput.h"
#include "AliMUONPixel.h"
#include "AliMC.h"

// Clusterizer class developped by Zitchenko (Dubna)
//
//
//


ClassImp(AliMUONClusterFinderAZ)
 
 const Double_t AliMUONClusterFinderAZ::fgkCouplMin = 1.e-3; // threshold on coupling 
 AliMUONClusterFinderAZ* AliMUONClusterFinderAZ::fgClusterFinder = 0x0;
 TMinuit* AliMUONClusterFinderAZ::fgMinuit = 0x0;


//_____________________________________________________________________________
AliMUONClusterFinderAZ::AliMUONClusterFinderAZ(Bool_t draw=0, Int_t iReco=0)
  : AliMUONClusterFinderVS()
{
// Constructor
  for (Int_t i=0; i<4; i++) {fHist[i] = 0;}
  fMuonDigits = 0;
  fSegmentation[1] = fSegmentation[0] = 0; 
  fgClusterFinder = 0x0;
  fgMinuit = 0x0; 
  if (!fgClusterFinder) fgClusterFinder = this;
  if (!fgMinuit) fgMinuit = new TMinuit(8);
  fDraw = draw;
  fReco = iReco;
  fPixArray = new TObjArray(20); 
  /*
  fPoints = 0;
  fPhits = 0;
  fRpoints = 0;
  fCanvas = 0;
  fNextCathode = kFALSE; 
  fColPad = 0;
  */
}

//_____________________________________________________________________________
AliMUONClusterFinderAZ::AliMUONClusterFinderAZ(const AliMUONClusterFinderAZ& rhs)
  : AliMUONClusterFinderVS(rhs)
{
// Protected copy constructor

  Fatal("AliMUONClusterFinderAZModule", "Not implemented.");
}

//_____________________________________________________________________________
AliMUONClusterFinderAZ::~AliMUONClusterFinderAZ()
{
  // Destructor
  delete fgMinuit; fgMinuit = 0; delete fPixArray; fPixArray = 0;
  /*
  // Delete space point structure
  if (fPoints) fPoints->Delete();
  delete fPoints;
  fPoints     = 0;
  //
  if (fPhits) fPhits->Delete();
  delete fPhits;
  fPhits     = 0;
  //
  if (fRpoints) fRpoints->Delete();
  delete fRpoints;
  fRpoints     = 0;
  */
}

//_____________________________________________________________________________
AliMUONClusterFinderAZ&  
AliMUONClusterFinderAZ::operator=(const AliMUONClusterFinderAZ& rhs)
{
// Protected assignement operator

  if (this == &rhs) return *this;

  Fatal("operator=", "Not implemented.");
    
  return *this;  
}    
          
//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindRawClusters()
{
// To provide the same interface as in AliMUONClusterFinderVS

  EventLoop (gAlice->GetHeader()->GetEvent(), AliMUONClusterInput::Instance()->Chamber());
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::EventLoop(Int_t nev=0, Int_t ch=0)
{
// Loop over events
  
  FILE *lun = 0;
  TCanvas *c1 = 0;
  TView *view = 0;
  TH2F *hist = 0;
  Double_t p1[3]={0}, p2[3];
  TTree *treeR = 0;
  if (fDraw) {
    // File
    lun = fopen("pool.dat","w");
    c1 = new TCanvas("c1","Clusters",0,0,600,700);
    c1->Divide(1,2);
    new TCanvas("c2","Mlem",700,0,600,350);
  }

newev:
  Int_t nparticles = 0, nent;

  //Loaders
  AliRunLoader * rl = AliRunLoader::GetRunLoader();
  AliLoader * gime  = rl->GetLoader("MUONLoader");

  if (!fReco) nparticles = rl->GetEvent(nev);
  else nparticles = gAlice->GetMCApp()->GetNtrack();
  cout << "nev         " << nev <<endl;
  cout << "nparticles  " << nparticles <<endl;
  if (nparticles <= 0) return;
  
  TTree *treeH = gime->TreeH();
  Int_t ntracks = (Int_t) treeH->GetEntries();
  cout<<"ntracks "<<ntracks<<endl;
    
  // Get pointers to Alice detectors and Digits containers
  AliMUON *muon  = (AliMUON*) gAlice->GetModule("MUON");
  if (!muon) return;
  //       TClonesArray *Particles = gAlice->Particles();
  if (!fReco) {
    treeR = gime->TreeR();
    if (treeR) {
      muon->ResetRawClusters();
      nent = (Int_t) treeR->GetEntries();
      if (nent != 1) {
	cout << "Error in MUONdrawClust" << endl;
	cout << " nent = " <<  nent << " not equal to 1" << endl;
	//exit(0);
      }
    } // if (treeR)
  } // if (!fReco)

  TTree *treeD = gime->TreeD();
  //muon->ResetDigits();

  TClonesArray *listMUONrawclust ;
  AliMUONChamber*       iChamber = 0;
  
  // As default draw the first cluster of the chamber #0
      
newchamber:
  if (ch > 9) {if (fReco) return; nev++; ch = 0; goto newev;}
  //gAlice->ResetDigits();
  fMuonDigits  = muon->GetMUONData()->Digits(ch);
  if (fMuonDigits == 0) return;
  iChamber = &(muon->Chamber(ch));
  fSegmentation[0] = iChamber->SegmentationModel(1);
  fSegmentation[1] = iChamber->SegmentationModel(2);
  fResponse = iChamber->ResponseModel();
    
  nent = 0;
 
  if (treeD) {
    nent = (Int_t) treeD->GetEntries();
    //printf(" entries %d \n", nent);
  }

  Int_t ndigits[2]={9,9}, nShown[2]={0};
  for (Int_t i=0; i<2; i++) {
    for (Int_t j=0; j<fgkDim; j++) {fUsed[i][j]=kFALSE;}
  }

next:
  if (ndigits[0] == nShown[0] && ndigits[1] == nShown[1]) {
    // No more clusters
    if (fReco) return;
    ch++;
    goto newchamber; // next chamber
  }
  Float_t xpad, ypad, zpad, zpad0;
  TLine *line[99]={0};
  Int_t nLine = 0;
  Bool_t first = kTRUE;
  cout << " *** Event # " << nev << " chamber: " << ch << endl;
  fnPads[0] = fnPads[1] = 0;
  for (Int_t i=0; i<fgkDim; i++) {fPadIJ[1][i] = 0;}
  //for (Int_t iii = 0; iii<999; iii++) { 
  for (Int_t iii = 0; iii<2; iii++) { 
    Int_t cath = TMath::Odd(iii);
    gAlice->ResetDigits();
    treeD->GetEvent(cath);
    fMuonDigits  = muon->GetMUONData()->Digits(ch);

    ndigits[cath] = fMuonDigits->GetEntriesFast();
    if (!ndigits[0] && !ndigits[1]) {if (fReco) return; ch++; goto newchamber;}
    if (ndigits[cath] == 0) continue;
    cout << " ndigits: " << ndigits[cath] << " " << cath << endl;

    AliMUONDigit  *mdig;
    Int_t digit;

    Bool_t eEOC = kTRUE; // end-of-cluster
    for (digit = 0; digit < ndigits[cath]; digit++) {
      mdig    = (AliMUONDigit*)fMuonDigits->UncheckedAt(digit);
      if (mdig->Cathode() != cath) continue;
      if (first) {
	// Find first unused pad
	if (fUsed[cath][digit]) continue;
	fSegmentation[cath]->GetPadC(mdig->PadX(),mdig->PadY(),xpad,ypad,zpad0);
      } else {
	if (fUsed[cath][digit]) continue;
	fSegmentation[cath]->GetPadC(mdig->PadX(),mdig->PadY(),xpad,ypad,zpad);
	if (TMath::Abs(zpad-zpad0)>0.1) continue; // different slats
	// Find a pad overlapping with the cluster
	if (!Overlap(cath,mdig)) continue;
      }
      // Add pad - recursive call
      AddPad(cath,digit);
      eEOC = kFALSE;
      if (digit >= 0) break;
    }
    if (first && eEOC) {
      // No more unused pads 
      if (cath == 0) continue; // on cathode #0 - check #1
      else {
	// No more clusters
	if (fReco) return;
	ch++;
	goto newchamber; // next chamber
      }
    }
    if (eEOC) break; // cluster found
    first = kFALSE;
    cout << " nPads: " << fnPads[cath] << " " << nShown[cath]+fnPads[cath] << " " << cath << endl;
  } // for (Int_t iii = 0;

  
  if (fReco) goto skip;
  char hName[4];
  for (Int_t cath = 0; cath<2; cath++) {
    // Build histograms
    if (fHist[cath*2]) {fHist[cath*2]->Delete(); fHist[cath*2] = 0;}
    if (fHist[cath*2+1]) {fHist[cath*2+1]->Delete(); fHist[cath*2+1] = 0;}
    if (fnPads[cath] == 0) continue; // cluster on one cathode only
    Float_t wxMin=999, wxMax=0, wyMin=999, wyMax=0; 
    Int_t minDx=0, maxDx=0, minDy=0, maxDy=0;
    for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
      if (fPadIJ[0][i] != cath) continue;
      if (fXyq[3][i] < wxMin) {wxMin = fXyq[3][i]; minDx = i;}
      if (fXyq[3][i] > wxMax) {wxMax = fXyq[3][i]; maxDx = i;}
      if (fXyq[4][i] < wyMin) {wyMin = fXyq[4][i]; minDy = i;}
      if (fXyq[4][i] > wyMax) {wyMax = fXyq[4][i]; maxDy = i;}
    }
    cout << minDx << maxDx << minDy << maxDy << endl;
    Int_t nx, ny, padSize;
    Float_t xmin=9999, xmax=-9999, ymin=9999, ymax=-9999;
    if (TMath::Nint(fXyq[3][minDx]*1000) == TMath::Nint(fXyq[3][maxDx]*1000) &&
	TMath::Nint(fXyq[4][minDy]*1000) == TMath::Nint(fXyq[4][maxDy]*1000)) {
      // the same segmentation
      cout << " Same" << endl;
      cout << fXyq[3][minDx] << " " << fXyq[3][maxDx] << " " << fXyq[4][minDy] << " " << fXyq[4][maxDy] << endl;
      for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
	if (fPadIJ[0][i] != cath) continue;
	if (fXyq[0][i] < xmin) xmin = fXyq[0][i];
	if (fXyq[0][i] > xmax) xmax = fXyq[0][i];
	if (fXyq[1][i] < ymin) ymin = fXyq[1][i];
	if (fXyq[1][i] > ymax) ymax = fXyq[1][i];
      }
      xmin -= fXyq[3][minDx]; xmax += fXyq[3][minDx];
      ymin -= fXyq[4][minDy]; ymax += fXyq[4][minDy];
      nx = TMath::Nint ((xmax-xmin)/wxMin/2);
      ny = TMath::Nint ((ymax-ymin)/wyMin/2);
      sprintf(hName,"h%d",cath*2);
      fHist[cath*2] = new TH2F(hName,"cluster",nx,xmin,xmax,ny,ymin,ymax);
      cout << fHist[cath*2] << " " << fnPads[cath] << endl;
      for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
	if (fPadIJ[0][i] != cath) continue;
	fHist[cath*2]->Fill(fXyq[0][i],fXyq[1][i],fXyq[2][i]);
	//cout << fXyq[0][i] << fXyq[1][i] << fXyq[2][i] << endl;
      }
    } else {
      // different segmentation in the cluster
      cout << " Different" << endl;
      cout << fXyq[3][minDx] << " " << fXyq[3][maxDx] << " " << fXyq[4][minDy] << " " << fXyq[4][maxDy] << endl;
      Int_t nOK = 0;
      Int_t indx, locMin, locMax;
      if (TMath::Nint(fXyq[3][minDx]*1000) != TMath::Nint(fXyq[3][maxDx]*1000)) {
	// different segmentation along x
	indx = 0;
	locMin = minDx;
	locMax = maxDx;
      } else {
	// different segmentation along y
	indx = 1;
	locMin = minDy;
	locMax = maxDy;
      }
      Int_t loc = locMin;
      for (Int_t i=0; i<2; i++) {
	// loop over different pad sizes
	if (i>0) loc = locMax;
	padSize = TMath::Nint(fXyq[indx+3][loc]*1000);
	xmin = 9999; xmax = -9999; ymin = 9999; ymax = -9999;
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	  if (fPadIJ[0][j] != cath) continue;
	  if (TMath::Nint(fXyq[indx+3][j]*1000) != padSize) continue;
	  nOK++;
	  xmin = TMath::Min (xmin,fXyq[0][j]);
	  xmax = TMath::Max (xmax,fXyq[0][j]);
	  ymin = TMath::Min (ymin,fXyq[1][j]);
	  ymax = TMath::Max (ymax,fXyq[1][j]);
	}
	xmin -= fXyq[3][loc]; xmax += fXyq[3][loc];
	ymin -= fXyq[4][loc]; ymax += fXyq[4][loc];
	nx = TMath::Nint ((xmax-xmin)/fXyq[3][loc]/2);
	ny = TMath::Nint ((ymax-ymin)/fXyq[4][loc]/2);
	sprintf(hName,"h%d",cath*2+i);
	fHist[cath*2+i] = new TH2F(hName,"cluster",nx,xmin,xmax,ny,ymin,ymax);
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	  if (fPadIJ[0][j] != cath) continue;
	  if (TMath::Nint(fXyq[indx+3][j]*1000) != padSize) continue;
	  fHist[cath*2+i]->Fill(fXyq[0][j],fXyq[1][j],fXyq[2][j]);
	}
      } // for (Int_t i=0;
      if (nOK != fnPads[cath]) cout << " *** Too many segmentations: nPads, nOK " << fnPads[cath] << " " << nOK << endl;
    } // if (TMath::Nint(fXyq[3][minDx]*1000)
  } // for (Int_t cath = 0;
	
  // Draw histograms and coordinates
  for (Int_t cath=0; cath<2; cath++) {
    if (cath == 0) ModifyHistos();
    if (fnPads[cath] == 0) continue; // cluster on one cathode only
    if (fDraw) {
      c1->cd(cath+1);
      gPad->SetTheta(55);
      gPad->SetPhi(30);
      Double_t x, y, x0, y0, r1=999, r2=0;
      if (fHist[cath*2+1]) {
	// 
	x0 = fHist[cath*2]->GetXaxis()->GetXmin() - 1000*TMath::Cos(30*TMath::Pi()/180);
	y0 = fHist[cath*2]->GetYaxis()->GetXmin() - 1000*TMath::Sin(30*TMath::Pi()/180);
	r1 = 0;
	Int_t ihist=cath*2;
	for (Int_t iy=1; iy<=fHist[ihist]->GetNbinsY(); iy++) {
	  y = fHist[ihist]->GetYaxis()->GetBinCenter(iy) 
	    + fHist[ihist]->GetYaxis()->GetBinWidth(iy);
	  for (Int_t ix=1; ix<=fHist[ihist]->GetNbinsX(); ix++) {
	    if (fHist[ihist]->GetCellContent(ix,iy) > 0.1) {
	      x = fHist[ihist]->GetXaxis()->GetBinCenter(ix)
		+ fHist[ihist]->GetXaxis()->GetBinWidth(ix);
	      r1 = TMath::Max (r1,TMath::Sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)));
	    }
	  }
	}
	ihist = cath*2 + 1 ;
	for (Int_t iy=1; iy<=fHist[ihist]->GetNbinsY(); iy++) {
	  y = fHist[ihist]->GetYaxis()->GetBinCenter(iy)
	    + fHist[ihist]->GetYaxis()->GetBinWidth(iy);
	  for (Int_t ix=1; ix<=fHist[ihist]->GetNbinsX(); ix++) {
	    if (fHist[ihist]->GetCellContent(ix,iy) > 0.1) {
	      x = fHist[ihist]->GetXaxis()->GetBinCenter(ix)
		+ fHist[ihist]->GetXaxis()->GetBinWidth(ix);
	      r2 = TMath::Max (r2,TMath::Sqrt((x-x0)*(x-x0)+(y-y0)*(y-y0)));
	    }
	  }
	}
	cout << r1 << " " << r2 << endl;
      } // if (fHist[cath*2+1])
      if (r1 > r2) {
	//fHist[cath*2]->Draw("lego1");
	fHist[cath*2]->Draw("lego1Fb");
	//if (fHist[cath*2+1]) fHist[cath*2+1]->Draw("lego1SameAxisBb");
	if (fHist[cath*2+1]) fHist[cath*2+1]->Draw("lego1SameAxisBbFb");
      } else {
	//fHist[cath*2+1]->Draw("lego1");
	fHist[cath*2+1]->Draw("lego1Fb");
	//fHist[cath*2]->Draw("lego1SameAxisBb");
	fHist[cath*2]->Draw("lego1SameAxisFbBb");
      }
      c1->Update();
    } // if (fDraw)
  } // for (Int_t cath = 0;

  // Draw generated hits
  Double_t xNDC[6];
  hist = fHist[0] ? fHist[0] : fHist[2];
  p2[2] = hist->GetMaximum();
  view = 0;
  if (c1) view = c1->Pad()->GetView();
  cout << " *** GEANT hits *** " << endl;
  fnMu = 0;
  Int_t ix, iy, iok;
  for (Int_t i=0; i<ntracks; i++) {
    treeH->GetEvent(i);
    for (AliMUONHit* mHit=(AliMUONHit*)muon->FirstHit(-1); 
	 mHit;
	 mHit=(AliMUONHit*)muon->NextHit()) {
      if (mHit->Chamber() != ch+1) continue;  // chamber number
      if (TMath::Abs(mHit->Z()-zpad0) > 1) continue; // different slat
      p2[0] = p1[0] = mHit->X();        // x-pos of hit
      p2[1] = p1[1] = mHit->Y();        // y-pos
      if (p1[0] < hist->GetXaxis()->GetXmin() || 
	  p1[0] > hist->GetXaxis()->GetXmax()) continue;
      if (p1[1] < hist->GetYaxis()->GetXmin() || 
	  p1[1] > hist->GetYaxis()->GetXmax()) continue;
      // Check if track comes thru pads with signal
      iok = 0;
      for (Int_t ihist=0; ihist<4; ihist++) {
	if (!fHist[ihist]) continue;
	ix = fHist[ihist]->GetXaxis()->FindBin(p1[0]);
	iy = fHist[ihist]->GetYaxis()->FindBin(p1[1]);
	if (fHist[ihist]->GetCellContent(ix,iy) > 0.5) {iok = 1; break;}
      }
      if (!iok) continue;
      gStyle->SetLineColor(1);
      if (TMath::Abs((Int_t)mHit->Particle()) == 13) {
	gStyle->SetLineColor(4);
	fnMu++;
	if (fnMu <= 2) {
	  fxyMu[fnMu-1][0] = p1[0];
	  fxyMu[fnMu-1][1] = p1[1];
	}
      }	    
      printf(" X=%10.4f, Y=%10.4f, Z=%10.4f\n",p1[0],p1[1],mHit->Z());
      if (view) {
	view->WCtoNDC(p1, &xNDC[0]);
	view->WCtoNDC(p2, &xNDC[3]);
	for (Int_t ipad=1; ipad<3; ipad++) {
	  c1->cd(ipad);
	  //c1->DrawLine(xpad[0],xpad[1],xpad[3],xpad[4]);
	  line[nLine] = new TLine(xNDC[0],xNDC[1],xNDC[3],xNDC[4]);
	  line[nLine++]->Draw();
	}
      }
    } // for (AliMUONHit* mHit=
  } // for (Int_t i=0; i<ntracks;

  // Draw reconstructed coordinates
  listMUONrawclust  = muon->GetMUONData()->RawClusters(ch);
  treeR->GetEvent(ch);
  //cout << listMUONrawclust  << " " << listMUONrawclust ->GetEntries() << endl;
  AliMUONRawCluster *mRaw;
  gStyle->SetLineColor(3);
  cout << " *** Reconstructed hits *** " << endl;
  for (Int_t i=0; i<listMUONrawclust ->GetEntries(); i++) {
    mRaw = (AliMUONRawCluster*)listMUONrawclust ->UncheckedAt(i);
    if (TMath::Abs(mRaw->GetZ(0)-zpad0) > 1) continue; // different slat
    p2[0] = p1[0] = mRaw->GetX(0);        // x-pos of hit
    p2[1] = p1[1] = mRaw->GetY(0);        // y-pos
    if (p1[0] < hist->GetXaxis()->GetXmin() || 
	p1[0] > hist->GetXaxis()->GetXmax()) continue;
    if (p1[1] < hist->GetYaxis()->GetXmin() || 
	p1[1] > hist->GetYaxis()->GetXmax()) continue;
    /*
      treeD->GetEvent(cath);
      cout << mRaw->fMultiplicity[0] << mRaw->fMultiplicity[1] << endl;
      for (Int_t j=0; j<mRaw->fMultiplicity[cath]; j++) {
      Int_t digit = mRaw->fIndexMap[j][cath];
      cout << ((AliMUONDigit*)fMuonDigits->UncheckedAt(digit))->Signal() << endl;
      }
    */
    // Check if track comes thru pads with signal
    iok = 0;
    for (Int_t ihist=0; ihist<4; ihist++) {
      if (!fHist[ihist]) continue;
      ix = fHist[ihist]->GetXaxis()->FindBin(p1[0]);
      iy = fHist[ihist]->GetYaxis()->FindBin(p1[1]);
      if (fHist[ihist]->GetCellContent(ix,iy) > 0.5) {iok = 1; break;}
    }
    if (!iok) continue;
    printf(" X=%10.4f, Y=%10.4f, Z=%10.4f\n",p1[0],p1[1],mRaw->GetZ(0));
    if (view) {
      view->WCtoNDC(p1, &xNDC[0]);
      view->WCtoNDC(p2, &xNDC[3]);
      for (Int_t ipad=1; ipad<3; ipad++) {
	c1->cd(ipad);
	line[nLine] = new TLine(xNDC[0],xNDC[1],xNDC[3],xNDC[4]);
	line[nLine++]->Draw();
      }
    }
  } // for (Int_t i=0; i<listMUONrawclust ->GetEntries();
  if (fDraw) c1->Update();

skip:
  // Use MLEM for cluster finder
  fZpad = zpad0;
  Int_t nMax = 1, localMax[100], maxPos[100];
  Double_t maxVal[100];
  
  if (CheckPrecluster(nShown)) {
    BuildPixArray();
    if (fnPads[0]+fnPads[1] > 50) nMax = FindLocalMaxima(localMax, maxVal);
    if (nMax > 1) TMath::Sort(nMax, maxVal, maxPos, kTRUE); // in decreasing order
    for (Int_t i=0; i<nMax; i++) {
      if (nMax > 1) FindCluster(localMax, maxPos[i]);
      if (!MainLoop()) cout << " MainLoop failed " << endl;
      if (i < nMax-1) {
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	  if (fPadIJ[1][j] == 0) continue; // pad charge was not modified
	  fPadIJ[1][j] = 0;
	  fXyq[2][j] = fXyq[5][j]; // use backup charge value
	}
      }
    }
  }
  if (fReco) goto next;

  for (Int_t i=0; i<fnMu; i++) {
    // Check again if muon come thru the used pads (due to extra splitting)
    for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
      if (TMath::Abs(fxyMu[i][0]-fXyq[0][j])<fXyq[3][j] && 
	  TMath::Abs(fxyMu[i][1]-fXyq[1][j])<fXyq[4][j]) {
	printf("%12.3e %12.3e %12.3e %12.3e\n",fxyMu[i][2],fxyMu[i][3],fxyMu[i][4],fxyMu[i][5]);
	if (lun) fprintf(lun,"%4d %2d %12.3e %12.3e %12.3e %12.3e\n",nev,ch,fxyMu[i][2],fxyMu[i][3],fxyMu[i][4],fxyMu[i][5]);
	break;
      }
    }
  } // for (Int_t i=0; i<fnMu;

  // What's next?
  char command[8];
  cout << " What is next? " << endl;
  command[0] = ' '; 
  if (fDraw) gets(command);
  if (command[0] == 'n' || command[0] == 'N') {nev++; goto newev;} // next event 
  else if (command[0] == 'q' || command[0] == 'Q') {fclose(lun); return;} // exit display 
  //else if (command[0] == 'r' || command[0] == 'R') goto redraw; // redraw points
  else if (command[0] == 'c' || command[0] == 'C') {
    // new chamber
    sscanf(command+1,"%d",&ch);
    goto newchamber;
  } 
  else if (command[0] == 'e' || command[0] == 'E') {
    // new event
    sscanf(command+1,"%d",&nev);
    goto newev;
  } 
  else goto next; // Next cluster
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::ModifyHistos(void)
{
  // Modify histograms to bring them to the same size
  Int_t nhist = 0;
  Float_t hlim[4][4], hbin[4][4]; // first index - xmin, xmax, ymin, ymax
  Float_t binMin[4] = {999,999,999,999};

  for (Int_t i=0; i<4; i++) {
    if (!fHist[i]) continue;
    hlim[0][nhist] = fHist[i]->GetXaxis()->GetXmin(); // xmin
    hlim[1][nhist] = fHist[i]->GetXaxis()->GetXmax(); // xmax
    hlim[2][nhist] = fHist[i]->GetYaxis()->GetXmin(); // ymin
    hlim[3][nhist] = fHist[i]->GetYaxis()->GetXmax(); // ymax
    hbin[0][nhist] = hbin[1][nhist] = fHist[i]->GetXaxis()->GetBinWidth(1);
    hbin[2][nhist] = hbin[3][nhist] = fHist[i]->GetYaxis()->GetBinWidth(1);
    binMin[0] = TMath::Min(binMin[0],hbin[0][nhist]);
    binMin[2] = TMath::Min(binMin[2],hbin[2][nhist]);
    nhist++;
  }
  binMin[1] = binMin[0];
  binMin[3] = binMin[2];
  cout << " Nhist: " << nhist << endl;

  Int_t imin, imax;
  for (Int_t lim=0; lim<4; lim++) {
    while (1) {
      imin = TMath::LocMin(nhist,hlim[lim]);
      imax = TMath::LocMax(nhist,hlim[lim]);
      if (TMath::Abs(hlim[lim][imin]-hlim[lim][imax])<0.01*binMin[lim]) break;
      if (lim == 0 || lim == 2) {
	// find lower limit
	hlim[lim][imax] -= hbin[lim][imax];
      } else {
	// find upper limit
	hlim[lim][imin] += hbin[lim][imin];
      }
    } // while (1)
  }
    
  // Rebuild histograms 
  nhist = 0;
  TH2F *hist = 0;
  Int_t nx, ny;
  Double_t x, y, cont, cmax=0;
  char hName[4];
  for (Int_t ihist=0; ihist<4; ihist++) {
    if (!fHist[ihist]) continue;
    nx = TMath::Nint((hlim[1][nhist]-hlim[0][nhist])/hbin[0][nhist]);
    ny = TMath::Nint((hlim[3][nhist]-hlim[2][nhist])/hbin[2][nhist]);
    //hist =  new TH2F("h","hist",nx,hlim[0][nhist],hlim[1][nhist],ny,hlim[2][nhist],hlim[3][nhist]);
    sprintf(hName,"hh%d",ihist);
    hist =  new TH2F(hName,"hist",nx,hlim[0][nhist],hlim[1][nhist],ny,hlim[2][nhist],hlim[3][nhist]);
    for (Int_t i=1; i<=fHist[ihist]->GetNbinsX(); i++) {
      x = fHist[ihist]->GetXaxis()->GetBinCenter(i);
      for (Int_t j=1; j<=fHist[ihist]->GetNbinsY(); j++) {
	y = fHist[ihist]->GetYaxis()->GetBinCenter(j);
	cont = fHist[ihist]->GetCellContent(i,j);
	hist->Fill(x,y,cont);
      }
    }
    cmax = TMath::Max (cmax,hist->GetMaximum());
    fHist[ihist]->Delete();
    fHist[ihist] = new TH2F(*hist);
    hist->Delete(); 
    nhist++;
  }
  printf("%f \n",cmax);

  for (Int_t ihist=0; ihist<4; ihist++) {
    if (!fHist[ihist]) continue;
    fHist[ihist]->SetMaximum(cmax);
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddPad(Int_t cath, Int_t digit)
{
  // Add pad to the cluster
  AliMUONDigit *mdig = (AliMUONDigit*)fMuonDigits->UncheckedAt(digit);

  Int_t charge = mdig->Signal();
  // get the center of the pad
  Float_t xpad, ypad, zpad;
  fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad);
  
  Int_t   isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());
  Int_t nPads = fnPads[0] + fnPads[1];
  fXyq[0][nPads] = xpad;
  fXyq[1][nPads] = ypad;
  fXyq[2][nPads] = charge;
  fXyq[3][nPads] = fSegmentation[cath]->Dpx(isec)/2;
  fXyq[4][nPads] = fSegmentation[cath]->Dpy(isec)/2;
  fXyq[5][nPads] = digit;
  fPadIJ[0][nPads] = cath;
  fPadIJ[1][nPads] = 0;
  fUsed[cath][digit] = kTRUE;
  //cout << " bbb " << fXyq[cath][2][nPads] << " " << fXyq[cath][0][nPads] << " " << fXyq[cath][1][nPads] << " " << fXyq[cath][3][nPads] << " " << fXyq[cath][4][nPads] << " " << zpad << " " << nPads << endl;
  fnPads[cath]++;

  // Check neighbours
  Int_t nn, ix, iy, xList[10], yList[10];
  AliMUONDigit  *mdig1;

  Int_t ndigits = fMuonDigits->GetEntriesFast();
  fSegmentation[cath]->Neighbours(mdig->PadX(),mdig->PadY(),&nn,xList,yList); 
  for (Int_t in=0; in<nn; in++) {
    ix=xList[in];
    iy=yList[in];
    for (Int_t digit1 = 0; digit1 < ndigits; digit1++) {
      if (digit1 == digit) continue;
      mdig1 = (AliMUONDigit*)fMuonDigits->UncheckedAt(digit1);
      if (mdig1->Cathode() != cath) continue;
      if (!fUsed[cath][digit1] && mdig1->PadX() == ix && mdig1->PadY() == iy) {
	fUsed[cath][digit1] = kTRUE;
	// Add pad - recursive call
	AddPad(cath,digit1);
      }
    } //for (Int_t digit1 = 0;
  } // for (Int_t in=0;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::Overlap(Int_t cath, TObject *dig)
{
  // Check if the pad from one cathode overlaps with a pad 
  // in the precluster on the other cathode

  AliMUONDigit *mdig = (AliMUONDigit*) dig;

  Float_t xpad, ypad, zpad;
  fSegmentation[cath]->GetPadC(mdig->PadX(), mdig->PadY(), xpad, ypad, zpad);
  Int_t   isec = fSegmentation[cath]->Sector(mdig->PadX(), mdig->PadY());

  Float_t xy1[4], xy12[4];
  xy1[0] = xpad - fSegmentation[cath]->Dpx(isec)/2;
  xy1[1] = xy1[0] + fSegmentation[cath]->Dpx(isec);
  xy1[2] = ypad - fSegmentation[cath]->Dpy(isec)/2;
  xy1[3] = xy1[2] + fSegmentation[cath]->Dpy(isec);
  //cout << " ok " << fnPads[0]+fnPads[1] << xy1[0] << xy1[1] << xy1[2] << xy1[3] << endl;

  Int_t cath1 = TMath::Even(cath);
  for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {
    if (fPadIJ[0][i] != cath1) continue;
    if (Overlap(xy1, i, xy12, 0)) return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::Overlap(Float_t *xy1, Int_t iPad, Float_t *xy12, Int_t iSkip)
{
  // Check if the pads xy1 and iPad overlap and return overlap area

  Float_t xy2[4];
  xy2[0] = fXyq[0][iPad] - fXyq[3][iPad];
  xy2[1] = fXyq[0][iPad] + fXyq[3][iPad];
  if (xy1[0] > xy2[1]-1.e-4 || xy1[1] < xy2[0]+1.e-4) return kFALSE;
  xy2[2] = fXyq[1][iPad] - fXyq[4][iPad];
  xy2[3] = fXyq[1][iPad] + fXyq[4][iPad];
  if (xy1[2] > xy2[3]-1.e-4 || xy1[3] < xy2[2]+1.e-4) return kFALSE;
  if (!iSkip) return kTRUE; // just check overlap (w/out computing the area)
  xy12[0] = TMath::Max (xy1[0],xy2[0]);
  xy12[1] = TMath::Min (xy1[1],xy2[1]);
  xy12[2] = TMath::Max (xy1[2],xy2[2]);
  xy12[3] = TMath::Min (xy1[3],xy2[3]);
  return kTRUE;
}

//_____________________________________________________________________________
/*
Bool_t AliMUONClusterFinderAZ::Overlap(Int_t i, Int_t j, Float_t *xy12, Int_t iSkip)
{
  // Check if the pads i and j overlap and return overlap area

  Float_t xy1[4], xy2[4];
  return Overlap(xy1, xy2, xy12, iSkip);
}
*/
//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::CheckPrecluster(Int_t *nShown)
{
  // Check precluster in order to attempt to simplify it (mostly for
  // two-cathode preclusters)

  Int_t i1, i2;
  Float_t xy1[4], xy12[4];
  
  Int_t npad = fnPads[0] + fnPads[1];

  // If pads have the same size take average of pads on both cathodes 
  Int_t sameSize = (fnPads[0] && fnPads[1]) ? 1 : 0;
  if (sameSize) {
    Double_t xSize = -1, ySize = 0;
    for (Int_t i=0; i<npad; i++) {
      if (fXyq[2][i] < 0) continue;
      if (xSize < 0) { xSize = fXyq[3][i]; ySize = fXyq[4][i]; }
      if (TMath::Abs(xSize-fXyq[3][i]) > 1.e-4 ||  TMath::Abs(ySize-fXyq[4][i]) > 1.e-4) { sameSize = 0; break; }
    }
  } // if (sameSize)
  if (sameSize && (fnPads[0] > 2 || fnPads[1] > 2)) {
    nShown[0] += fnPads[0];
    nShown[1] += fnPads[1];
    fnPads[0] = fnPads[1] = 0;
    Int_t div;
    for (Int_t i=0; i<npad; i++) {
      if (fXyq[2][i] < 0) continue; // used pad
      fXyq[2][fnPads[0]] = fXyq[2][i];
      div = 1;
      for (Int_t j=i+1; j<npad; j++) {
	if (fPadIJ[0][j] == fPadIJ[0][i]) continue; // same cathode
	if (TMath::Abs(fXyq[0][j]-fXyq[0][i]) > 1.e-4) continue;
	if (TMath::Abs(fXyq[1][j]-fXyq[1][i]) > 1.e-4) continue;
	fXyq[2][fnPads[0]] += fXyq[2][j];
	div = 2;
	fXyq[2][j] = -2;
	break;
      }
      fXyq[2][fnPads[0]] /= div;
      fXyq[0][fnPads[0]] = fXyq[0][i];
      fXyq[1][fnPads[0]] = fXyq[1][i];
      fPadIJ[0][fnPads[0]++] = 0;
    }
  } // if (sameSize)

  // Check if one-cathode precluster
  i1 = fnPads[0]!=0 ? 0 : 1;
  i2 = fnPads[1]!=0 ? 1 : 0;

  if (i1 != i2) { // two-cathode 

    Int_t *flags = new Int_t[npad];
    for (Int_t i=0; i<npad; i++) { flags[i] = 0; }

    // Check pad overlaps
    for (Int_t i=0; i<npad; i++) {
      if (fPadIJ[0][i] != i1) continue;
      xy1[0] = fXyq[0][i] - fXyq[3][i];
      xy1[1] = fXyq[0][i] + fXyq[3][i];
      xy1[2] = fXyq[1][i] - fXyq[4][i];
      xy1[3] = fXyq[1][i] + fXyq[4][i];
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[0][j] != i2) continue;
	if (!Overlap(xy1, j, xy12, 0)) continue;
	flags[i] = flags[j] = 1; // mark overlapped pads
      } // for (Int_t j=0;
    } // for (Int_t i=0;

    // Check if all pads overlap
    Int_t digit=0, cath, nFlags=0;
    for (Int_t i=0; i<npad; i++) {nFlags += !flags[i];}
    if (nFlags) cout << " nFlags = " << nFlags << endl;
    //if (nFlags > 2 || (Float_t)nFlags / npad > 0.2) { // why 2 ??? - empirical choice
    if (nFlags > 0) {
      for (Int_t i=0; i<npad; i++) {
	if (flags[i]) continue;
	digit = TMath::Nint (fXyq[5][i]);
	cath = fPadIJ[0][i];
	fUsed[cath][digit] = kFALSE; // release pad
	fXyq[2][i] = -2;
	fnPads[cath]--;
      }
    } // if (nFlags > 2)

    // Check correlations of cathode charges
    if (fnPads[0] && fnPads[1]) { // two-cathode
      Double_t sum[2]={0};
      Int_t over[2] = {1, 1};
      for (Int_t i=0; i<npad; i++) {
	cath = fPadIJ[0][i];
	if (fXyq[2][i] > 0) sum[cath] += fXyq[2][i];
	if (fXyq[2][i] > fResponse->MaxAdc()-1) over[cath] = 0;
      }
      cout << " Total charge: " << sum[0] << " " << sum[1] << endl;
      if ((over[0] || over[1]) && TMath::Abs(sum[0]-sum[1])/(sum[0]+sum[1])*2 > 1) { // 3 times difference
	cout << " Release " << endl;
	// Big difference
	cath = sum[0]>sum[1] ? 0 : 1;
	Int_t imax = 0;
	Double_t cmax=-1;
	Double_t *dist = new Double_t[npad];
	for (Int_t i=0; i<npad; i++) {
	  if (fPadIJ[0][i] != cath) continue;
	  if (fXyq[2][i] < cmax) continue;
	  cmax = fXyq[2][i];
	  imax = i;
	}
	// Arrange pads according to their distance to the max, 
	// normalized to the pad size
	for (Int_t i=0; i<npad; i++) {
	  dist[i] = 0;
	  if (fPadIJ[0][i] != cath) continue;
	  if (i == imax) continue; 
	  if (fXyq[2][i] < 0) continue;
	  dist[i] = (fXyq[0][i]-fXyq[0][imax])*(fXyq[0][i]-fXyq[0][imax])/
	             fXyq[3][imax]/fXyq[3][imax]/4;
	  dist[i] += (fXyq[1][i]-fXyq[1][imax])*(fXyq[1][i]-fXyq[1][imax])/
                      fXyq[4][imax]/fXyq[4][imax]/4;
	  dist[i] = TMath::Sqrt (dist[i]);
	}
	TMath::Sort(npad, dist, flags, kFALSE); // in increasing order
	Int_t indx;
	Double_t xmax = -1;
	for (Int_t i=0; i<npad; i++) {
	  indx = flags[i];
	  if (fPadIJ[0][indx] != cath) continue;
	  if (fXyq[2][indx] < 0) continue;
	  if (fXyq[2][indx] <= cmax || TMath::Abs(dist[indx]-xmax)<1.e-3) {
	    // Release pads
	    if (TMath::Abs(dist[indx]-xmax)<1.e-3) 
                cmax = TMath::Max((Double_t)(fXyq[2][indx]),cmax);
	    else cmax = fXyq[2][indx];
	    xmax = dist[indx];
	    digit = TMath::Nint (fXyq[5][indx]);
	    fUsed[cath][digit] = kFALSE; 
	    fXyq[2][indx] = -2;
	    fnPads[cath]--;
	    // xmax = dist[i]; // Bug?
	  }
	  else break;
	} 
	delete [] dist; dist = 0;
      } // TMath::Abs(sum[0]-sum[1])...
    } // if (fnPads[0] && fnPads[1])
    delete [] flags; flags = 0;
  } // if (i1 != i2) 

  if (!sameSize) { nShown[0] += fnPads[0]; nShown[1] += fnPads[1]; }

  // Move released pads to the right
  Int_t beg = 0, end = npad-1, padij;
  Double_t xyq;
  while (beg < end) {
    if (fXyq[2][beg] > 0) { beg++; continue; }
    for (Int_t j=end; j>beg; j--) {
      if (fXyq[2][j] < 0) continue;
      end = j - 1;
      for (Int_t j1=0; j1<2; j1++) {
	padij = fPadIJ[j1][beg]; 
	fPadIJ[j1][beg] = fPadIJ[j1][j];
	fPadIJ[j1][j] = padij;
      }
      for (Int_t j1=0; j1<6; j1++) {
	xyq = fXyq[j1][beg]; 
	fXyq[j1][beg] = fXyq[j1][j];
	fXyq[j1][j] = xyq;
      }
      break;
    } // for (Int_t j=end;
    beg++;
  } // while
  npad = fnPads[0] + fnPads[1];
  if (npad > 500) { cout << " ***** Too large cluster. Give up. " << npad << endl; return kFALSE; }
  // Back up charge value
  for (Int_t j=0; j<npad; j++) fXyq[5][j] = fXyq[2][j];

  return kTRUE;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::BuildPixArray()
{
  // Build pixel array for MLEM method
  
  Int_t nPix=0, i1, i2;
  Float_t xy1[4], xy12[4];
  AliMUONPixel *pixPtr=0;

  Int_t npad = fnPads[0] + fnPads[1];

  // One cathode is empty
  i1 = fnPads[0]!=0 ? 0 : 1;
  i2 = fnPads[1]!=0 ? 1 : 0;

  // Build array of pixels on anode plane
  if (i1 == i2) { // one-cathode precluster
    for (Int_t j=0; j<npad; j++) {
      pixPtr = new AliMUONPixel();
      for (Int_t i=0; i<2; i++) {
	pixPtr->SetCoord(i, fXyq[i][j]); // pixel coordinates
	pixPtr->SetSize(i, fXyq[i+3][j]); // pixel size
      }
      pixPtr->SetCharge(fXyq[2][j]); // charge
      fPixArray->Add((TObject*)pixPtr);
      nPix++;
    }
  } else { // two-cathode precluster    
    for (Int_t i=0; i<npad; i++) {
      if (fPadIJ[0][i] != i1) continue;
      xy1[0] = fXyq[0][i] - fXyq[3][i];
      xy1[1] = fXyq[0][i] + fXyq[3][i];
      xy1[2] = fXyq[1][i] - fXyq[4][i];
      xy1[3] = fXyq[1][i] + fXyq[4][i];
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[0][j] != i2) continue;
	if (!Overlap(xy1, j, xy12, 1)) continue;
	pixPtr = new AliMUONPixel();
	for (Int_t k=0; k<2; k++) {
	  pixPtr->SetCoord(k, (xy12[2*k]+xy12[2*k+1])/2); // pixel coordinates
	  pixPtr->SetSize(k, xy12[2*k+1]-pixPtr->Coord(k)); // size
	}
	pixPtr->SetCharge(TMath::Min (fXyq[2][i],fXyq[2][j])); //charge
	fPixArray->Add((TObject*)pixPtr);
	nPix++;
      } // for (Int_t j=0;
    } // for (Int_t i=0;
  } // else

  Float_t wxmin=999, wymin=999;
  for (Int_t i=0; i<npad; i++) {
    if (fPadIJ[0][i] == i1) wymin = TMath::Min (wymin,fXyq[4][i]);
    if (fPadIJ[0][i] == i2) wxmin = TMath::Min (wxmin,fXyq[3][i]);
  }
  cout << wxmin << " " << wymin << endl;

  // Check if small pixel X-size
  AjustPixel(wxmin, 0);
  // Check if small pixel Y-size
  AjustPixel(wymin, 1);
  // Check if large pixel size
  AjustPixel(wxmin, wymin);

  // Remove discarded pixels
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    //pixPtr->Print();
    if (pixPtr->Charge() < 1) { fPixArray->RemoveAt(i); delete pixPtr; }// discarded pixel
  }
  fPixArray->Compress();
  nPix = fPixArray->GetEntriesFast();

  if (nPix > npad) {
    cout << nPix << endl;
    // Too many pixels - sort and remove pixels with the lowest signal
    fPixArray->Sort();
    for (Int_t i=npad; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      //pixPtr->Print();
      fPixArray->RemoveAt(i);
      delete pixPtr;
    }
    nPix = npad;
  } // if (nPix > npad)

  // Set pixel charges to the same value (for MLEM)
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    //pixPtr->SetCharge(10);
    cout << i+1 << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << " " << pixPtr->Size(0) << " " << pixPtr->Size(1) << endl;
  }
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AjustPixel(Float_t width, Int_t ixy)
{
  // Check if some pixels have small size (ajust if necessary)

  AliMUONPixel *pixPtr, *pixPtr1 = 0;
  Int_t ixy1 = TMath::Even(ixy);
  Int_t nPix = fPixArray->GetEntriesFast();

  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(ixy)-width < -1.e-4) {
      // try to merge 
      cout << " Small X or Y: " << ixy << " " << pixPtr->Size(ixy) << " " << width << " " << pixPtr->Coord(0) << " " << pixPtr->Coord(1) << endl;
      for (Int_t j=i+1; j<nPix; j++) {
	pixPtr1 = (AliMUONPixel*) fPixArray->UncheckedAt(j);
	if (pixPtr1->Charge() < 1) continue; // discarded pixel
	if (TMath::Abs(pixPtr1->Size(ixy)-width) < 1.e-4) continue; // right size 
	if (TMath::Abs(pixPtr1->Coord(ixy1)-pixPtr->Coord(ixy1)) > 1.e-4) continue; // different rows/columns
	if (TMath::Abs(pixPtr1->Coord(ixy)-pixPtr->Coord(ixy)) < 2*width) {
	  // merge
	  pixPtr->SetSize(ixy, width);
	  pixPtr->SetCoord(ixy, (pixPtr->Coord(ixy)+pixPtr1->Coord(ixy))/2);
	  pixPtr->SetCharge(TMath::Min (pixPtr->Charge(),pixPtr1->Charge()));
	  pixPtr1->SetCharge(0);
	  pixPtr1 = 0;
	  break;
	}
      } // for (Int_t j=i+1;
      //if (!pixPtr1) { cout << " I am here!" << endl; pixPtr->SetSize(ixy, width); } // ???
      //else if (pixPtr1->Charge() > 0.5 || i == nPix-1) {
      if (pixPtr1 || i == nPix-1) {
	// edge pixel - just increase its size
	cout << " Edge ..." << endl; 
	for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
	  // ???if (fPadIJ[0][j] != i1) continue;
	  if (TMath::Abs(pixPtr->Coord(ixy1)-fXyq[ixy1][j]) > 1.e-4) continue;
	  if (pixPtr->Coord(ixy) < fXyq[ixy][j]) 
	    pixPtr->Shift(ixy, -pixPtr->Size(ixy));
	  else pixPtr->Shift(ixy, pixPtr->Size(ixy));
  	  pixPtr->SetSize(ixy, width);
	  break;
	}
      }
    } // if (pixPtr->Size(ixy)-width < -1.e-4)
  } // for (Int_t i=0; i<nPix;
  return;
}
  
//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AjustPixel(Float_t wxmin, Float_t wymin)
{
  // Check if some pixels have large size (ajust if necessary)

  Int_t nx, ny;
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr, *pixPtr1, pix;

  // Check if large pixel size
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 1) continue; // discarded pixel
    if (pixPtr->Size(0)-wxmin > 1.e-4 || pixPtr->Size(1)-wymin > 1.e-4) {
      cout << " Different " << pixPtr->Size(0) << " " << wxmin << " " << pixPtr->Size(1) << " " << wymin << endl;
      pix = *pixPtr;
      nx = TMath::Nint (pix.Size(0)/wxmin);
      ny = TMath::Nint (pix.Size(1)/wymin);
      pix.Shift(0, -pix.Size(0)-wxmin);
      pix.Shift(1, -pix.Size(1)-wymin);
      pix.SetSize(0, wxmin);
      pix.SetSize(1, wymin);
      for (Int_t ii=0; ii<nx; ii++) {
	pix.Shift(0, wxmin*2);
	for (Int_t jj=0; jj<ny; jj++) {
	  pix.Shift(1, wymin*2);
	  pixPtr1 = new AliMUONPixel(pix);
	  fPixArray->Add((TObject*)pixPtr1);
	}
      }
      pixPtr->SetCharge(0);
    }
  } // for (Int_t i=0; i<nPix;
  return;
}

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::MainLoop()
{
  // Repeat MLEM algorithm until pixel size becomes sufficiently small
  
  TH2D *mlem;

  Int_t ix, iy;
  //Int_t nn, xList[10], yList[10];
  Int_t nPix = fPixArray->GetEntriesFast();
  Int_t npadTot = fnPads[0] + fnPads[1], npadOK = 0;
  AliMUONPixel *pixPtr = 0;
  Double_t *coef = 0, *probi = 0; 
  for (Int_t i=0; i<npadTot; i++) if (fPadIJ[1][i] == 0) npadOK++;

  while (1) {

    mlem = (TH2D*) gROOT->FindObject("mlem");
    if (mlem) mlem->Delete();
    // Calculate coefficients
    cout << " nPix, npadTot, npadOK " << nPix << " " << npadTot << " " << npadOK << endl;

    // Calculate coefficients and pixel visibilities
    coef = new Double_t [npadTot*nPix];
    probi = new Double_t [nPix];
    Int_t indx = 0, cath;
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      probi[ipix] = 0;
      for (Int_t j=0; j<npadTot; j++) {
	if (fPadIJ[1][j] < 0) { coef[j*nPix+ipix] = 0; continue; }
	cath = fPadIJ[0][j];
	fSegmentation[cath]->GetPadI(fXyq[0][j],fXyq[1][j],fZpad,ix,iy);
	fSegmentation[cath]->SetPad(ix,iy);
	/*
	  fSegmentation[cath]->Neighbours(ix,iy,&nn,xList,yList); 
	  if (nn != 4) {
	  cout << nn << ": ";
	  for (Int_t i=0; i<nn; i++) {cout << xList[i] << " " << yList[i] << ", ";}
	  cout << endl;
	  }
	*/
	Double_t sum = 0;
	fSegmentation[cath]->SetHit(pixPtr->Coord(0),pixPtr->Coord(1),fZpad);
	sum += fResponse->IntXY(fSegmentation[cath]);
	indx = j*nPix + ipix;
	coef[indx] = sum; 
	probi[ipix] += coef[indx];
	//cout << j << " " << ipix << " " << coef[indx] << endl;
      } // for (Int_t j=0;
      //cout << " prob: " << probi[ipix] << endl;
      if (probi[ipix] < 0.01) pixPtr->SetCharge(0); // "invisible" pixel
    } // for (Int_t ipix=0;

    // MLEM algorithm
    Mlem(coef, probi);

    Double_t xylim[4] = {999, 999, 999, 999};
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      for (Int_t i=0; i<4; i++) 
	xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
      //cout << ipix+1; pixPtr->Print();
    }
    for (Int_t i=0; i<4; i++) {
      xylim[i] -= pixPtr->Size(i/2); cout << (i%2 ? -1 : 1)*xylim[i] << " "; }
    cout << endl;

    // Ajust histogram to approximately the same limits as for the pads
    // (for good presentation)
    //*
    Float_t xypads[4];
    if (fHist[0]) {
      xypads[0] = fHist[0]->GetXaxis()->GetXmin();
      xypads[1] = -fHist[0]->GetXaxis()->GetXmax();
      xypads[2] = fHist[0]->GetYaxis()->GetXmin();
      xypads[3] = -fHist[0]->GetYaxis()->GetXmax();
      for (Int_t i=0; i<4; i++) {
	while(1) {
	  if (xylim[i] < xypads[i]) break;
	  xylim[i] -= 2*pixPtr->Size(i/2);
	}
      }
    } // if (fHist[0])
    //*/

    Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
    Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);
    mlem = new TH2D("mlem","mlem",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      mlem->Fill(pixPtr->Coord(0),pixPtr->Coord(1),pixPtr->Charge());
    }
    //gPad->GetCanvas()->cd(3);
    if (fDraw) {
      ((TCanvas*)gROOT->FindObject("c2"))->cd();
      gPad->SetTheta(55);
      gPad->SetPhi(30);
      mlem->Draw("lego1Fb");
      gPad->Update();
      gets((char*)&ix);
    }

    // Check if the total charge of pixels is too low
    Double_t qTot = 0;
    for (Int_t i=0; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      qTot += pixPtr->Charge();
    }
    if (qTot < 1.e-4 || npadOK < 3 && qTot < 50) {
      delete [] coef; delete [] probi; coef = 0; probi = 0;
      fPixArray->Delete(); 
      return kFALSE; 
    }

    // Plot data - expectation
    /*
    Double_t x, y, cont;
    for (Int_t j=0; j<npadTot; j++) {
      Double_t sum1 = 0;
      for (Int_t i=0; i<nPix; i++) {
        // Caculate expectation
	pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
	sum1 += pixPtr->Charge()*coef[j*nPix+i];
      }
      sum1 = TMath::Min (sum1,(Double_t)fResponse->MaxAdc());
      x = fXyq[0][j];
      y = fXyq[1][j];
      cath = fPadIJ[0][j];
      Int_t ihist = cath*2;
      ix = fHist[ihist]->GetXaxis()->FindBin(x);
      iy = fHist[ihist]->GetYaxis()->FindBin(y);
      cont = fHist[ihist]->GetCellContent(ix,iy);
      if (cont == 0 && fHist[ihist+1]) {
        ihist += 1;
	ix = fHist[ihist]->GetXaxis()->FindBin(x);
	iy = fHist[ihist]->GetYaxis()->FindBin(y);
      }
      fHist[ihist]->SetBinContent(ix,iy,fXyq[2][j]-sum1);
    }
    ((TCanvas*)gROOT->FindObject("c1"))->cd(1);
    //gPad->SetTheta(55);
    //gPad->SetPhi(30);
    //mlem->Draw("lego1");
    gPad->Modified();
    ((TCanvas*)gROOT->FindObject("c1"))->cd(2);
    gPad->Modified();
    */

    // Calculate position of the center-of-gravity around the maximum pixel
    Double_t xyCOG[2];
    FindCOG(mlem, xyCOG);

    if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) < 0.07 && pixPtr->Size(0) > pixPtr->Size(1)) break;
    //if (TMath::Min(pixPtr->Size(0),pixPtr->Size(1)) >= 0.07 || pixPtr->Size(0) < pixPtr->Size(1)) {
    // Sort pixels according to the charge
    fPixArray->Sort();
    /*
    for (Int_t i=0; i<nPix; i++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
      cout << i+1; pixPtr->Print();
    }
    */
    Double_t pixMin = 0.01*((AliMUONPixel*)fPixArray->UncheckedAt(0))->Charge();
    pixMin = TMath::Min (pixMin,50.);

    // Decrease pixel size and shift pixels to make them centered at 
    // the maximum one
    indx = (pixPtr->Size(0)>pixPtr->Size(1)) ? 0 : 1;
    Double_t width = 0, shift[2]={0};
    ix = 1;
    for (Int_t i=0; i<4; i++) xylim[i] = 999;
    Int_t nPix1 = nPix; nPix = 0;
    for (Int_t ipix=0; ipix<nPix1; ipix++) {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      if (nPix >= npadOK) { // too many pixels already
	fPixArray->RemoveAt(ipix); 
	delete pixPtr; 
	continue;
      }
      if (pixPtr->Charge() < pixMin) { // low charge
	fPixArray->RemoveAt(ipix); 
	delete pixPtr; 
	continue;
      }
      for (Int_t i=0; i<2; i++) {
	if (!i) {
	  pixPtr->SetCharge(10);
	  pixPtr->SetSize(indx, pixPtr->Size(indx)/2);
	  width = -pixPtr->Size(indx);
	  pixPtr->Shift(indx, width);
	  // Shift pixel position
	  if (ix) {
	    ix = 0;
	    for (Int_t j=0; j<2; j++) {
	      shift[j] = pixPtr->Coord(j) - xyCOG[j];
	      shift[j] -= ((Int_t)(shift[j]/pixPtr->Size(j)/2))*pixPtr->Size(j)*2;
	    }
	    //cout << ipix << " " << i << " " << shift[0] << " " << shift[1] << endl;
	  } // if (ix)
	  pixPtr->Shift(0, -shift[0]);
	  pixPtr->Shift(1, -shift[1]);
	} else {
	  pixPtr = new AliMUONPixel(*pixPtr);
	  pixPtr->Shift(indx, -2*width);
	  fPixArray->Add((TObject*)pixPtr);
	} // else
	//pixPtr->Print();
	for (Int_t i=0; i<4; i++) 
	  xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
      } // for (Int_t i=0; i<2;
      nPix += 2;
    } // for (Int_t ipix=0;

    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();

    // Remove excessive pixels
    if (nPix > npadOK) {
      for (Int_t ipix=npadOK; ipix<nPix; ipix++) { 
	pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
	fPixArray->RemoveAt(ipix); 
	delete pixPtr;
      }
    } else {
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(0);
      // add pixels if the maximum is at the limit of pixel area
      // start from Y-direction
      Int_t j = 0;
      for (Int_t i=3; i>-1; i--) {
	if (nPix < npadOK && 
	    TMath::Abs((i%2 ? -1 : 1)*xylim[i]-xyCOG[i/2]) < pixPtr->Size(i/2)) {
	  pixPtr = new AliMUONPixel(*pixPtr);
	  pixPtr->SetCoord(i/2, xyCOG[i/2]+(i%2 ? 2:-2)*pixPtr->Size(i/2));
	  j = TMath::Even (i/2);
	  pixPtr->SetCoord(j, xyCOG[j]);
	  fPixArray->Add((TObject*)pixPtr);
	  nPix++;
	}
      }
    } // else    

    fPixArray->Compress();
    nPix = fPixArray->GetEntriesFast();
    delete [] coef; delete [] probi; coef = 0; probi = 0;
  } // while (1)

  // remove pixels with low signal or low visibility
  // Cuts are empirical !!!
  Double_t thresh = TMath::Max (mlem->GetMaximum()/100.,1.);
  thresh = TMath::Min (thresh,50.);
  Double_t cmax = -1, charge = 0;
  for (Int_t i=0; i<nPix; i++) cmax = TMath::Max (cmax,probi[i]); 
  // Mark pixels which should be removed
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    charge = pixPtr->Charge();
    if (charge < thresh) pixPtr->SetCharge(-charge);
    else if (cmax > 1.91) {
      if (probi[i] < 1.9) pixPtr->SetCharge(-charge);
    }
    else if (probi[i] < cmax*0.9) pixPtr->SetCharge(-charge);
  }
  // Move charge of removed pixels to their nearest neighbour (to keep total charge the same)
  Int_t near = 0;
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    charge = pixPtr->Charge();
    if (charge > 0) continue;
    near = FindNearest(pixPtr);
    pixPtr->SetCharge(0);
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(near);
    pixPtr->SetCharge(pixPtr->Charge() - charge);
  }
  // Update histogram
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    ix = mlem->GetXaxis()->FindBin(pixPtr->Coord(0));
    iy = mlem->GetYaxis()->FindBin(pixPtr->Coord(1));
    mlem->SetBinContent(ix, iy, pixPtr->Charge());
  }
  if (fDraw) {
    ((TCanvas*)gROOT->FindObject("c2"))->cd();
    gPad->SetTheta(55);
    gPad->SetPhi(30);
    mlem->Draw("lego1Fb");
    gPad->Update();
  }

  fxyMu[0][6] = fxyMu[1][6] = 9999;
  // Try to split into clusters
  Bool_t ok = kTRUE;
  if (mlem->GetSum() < 1) ok = kFALSE;
  else Split(mlem, coef);
  delete [] coef; delete [] probi; coef = 0; probi = 0;
  fPixArray->Delete(); 
  return ok;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Mlem(Double_t *coef, Double_t *probi)
{
  // Use MLEM to find pixel charges
  
  Int_t nPix = fPixArray->GetEntriesFast();
  Int_t npad = fnPads[0] + fnPads[1];
  Double_t *probi1 = new Double_t [nPix];
  Int_t indx, indx1;
  AliMUONPixel *pixPtr;

  for (Int_t iter=0; iter<15; iter++) {
    // Do iterations
    for (Int_t ipix=0; ipix<nPix; ipix++) {
      // Correct each pixel
      if (probi[ipix] < 0.01) continue; // skip "invisible" pixel
      Double_t sum = 0;
      probi1[ipix] = probi[ipix];
      for (Int_t j=0; j<npad; j++) {
	if (fPadIJ[1][j] < 0) continue; 
	Double_t sum1 = 0;
	indx1 = j*nPix;
	indx = indx1 + ipix;
	for (Int_t i=0; i<nPix; i++) {
	  // Caculate expectation
	  pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
	  sum1 += pixPtr->Charge()*coef[indx1+i];
	} // for (Int_t i=0;
	if (fXyq[2][j] > fResponse->MaxAdc()-1 && sum1 > fResponse->MaxAdc()) { probi1[ipix] -= coef[indx]; continue; } // correct for pad charge overflows
	//cout << sum1 << " " << fXyq[2][j] << " " << coef[j*nPix+ipix] << endl;
	if (coef[indx] > 1.e-6) sum += fXyq[2][j]*coef[indx]/sum1;
      } // for (Int_t j=0;
      pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
      if (probi1[ipix] > 1.e-6) pixPtr->SetCharge(pixPtr->Charge()*sum/probi1[ipix]);
    } // for (Int_t ipix=0;
  } // for (Int_t iter=0;
  delete [] probi1;
  return;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindCOG(TH2D *mlem, Double_t *xyc)
{
  // Calculate position of the center-of-gravity around the maximum pixel

  Int_t ixmax, iymax, ix, nsumx=0, nsumy=0, nsum=0;
  Int_t i1 = -9, j1 = -9;
  mlem->GetMaximumBin(ixmax,iymax,ix);
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t thresh = mlem->GetMaximum()/10;
  Double_t x, y, cont, xq=0, yq=0, qq=0;

  for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
    y = mlem->GetYaxis()->GetBinCenter(i);
    for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
      cont = mlem->GetCellContent(j,i);
      if (cont < thresh) continue;
      if (i != i1) {i1 = i; nsumy++;}
      if (j != j1) {j1 = j; nsumx++;}
      x = mlem->GetXaxis()->GetBinCenter(j);
      xq += x*cont;
      yq += y*cont;
      qq += cont;
      nsum++;
    }
  }

  Double_t cmax = 0;
  Int_t i2 = 0, j2 = 0;
  x = y = 0;
  if (nsumy == 1) {
    // one bin in Y - add one more (with the largest signal)
    for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
      if (i == iymax) continue;
      for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
	cont = mlem->GetCellContent(j,i);
	if (cont > cmax) {
	  cmax = cont;
	  x = mlem->GetXaxis()->GetBinCenter(j);
	  y = mlem->GetYaxis()->GetBinCenter(i);
	  i2 = i;
	  j2 = j;
	}
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumy == 1)

  if (nsumx == 1) {
    // one bin in X - add one more (with the largest signal)
    cmax = x = y = 0;
    for (Int_t j=TMath::Max(1,ixmax-1); j<=TMath::Min(nx,ixmax+1); j++) {
      if (j == ixmax) continue;
      for (Int_t i=TMath::Max(1,iymax-1); i<=TMath::Min(ny,iymax+1); i++) {
	cont = mlem->GetCellContent(j,i);
	if (cont > cmax) {
	  cmax = cont;
	  x = mlem->GetXaxis()->GetBinCenter(j);
	  y = mlem->GetYaxis()->GetBinCenter(i);
	  i2 = i;
	  j2 = j;
	}
      }
    }
    xq += x*cmax;
    yq += y*cmax;
    qq += cmax;
    if (i2 != i1) nsumy++;
    if (j2 != j1) nsumx++;
    nsum++;
  } // if (nsumx == 1)

  xyc[0] = xq/qq; xyc[1] = yq/qq;
  cout << xyc[0] << " " << xyc[1] << " " << qq << " " << nsum << " " << nsumx << " " << nsumy << endl;
  return;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::FindNearest(AliMUONPixel *pixPtr0)
{
  // Find the pixel nearest to the given one
  // (algorithm may be not very efficient)

  Int_t nPix = fPixArray->GetEntriesFast(), imin = 0;
  Double_t rmin = 99999, dx = 0, dy = 0, r = 0;
  Double_t xc = pixPtr0->Coord(0), yc = pixPtr0->Coord(1);
  AliMUONPixel *pixPtr;

  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 0.5) continue;
    dx = (xc - pixPtr->Coord(0)) / pixPtr->Size(0);
    dy = (yc - pixPtr->Coord(1)) / pixPtr->Size(1);
    r = dx *dx + dy * dy;
    if (r < rmin) { rmin = r; imin = i; }
  }
  return imin;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Split(TH2D *mlem, Double_t *coef)
{
  // The main steering function to work with clusters of pixels in anode
  // plane (find clusters, decouple them from each other, merge them (if
  // necessary), pick up coupled pads, call the fitting function)
  
  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Int_t nPix = fPixArray->GetEntriesFast();

  Bool_t *used = new Bool_t[ny*nx];
  Double_t cont;
  Int_t nclust = 0, indx, indx1;

  for (Int_t i=0; i<ny*nx; i++) used[i] = kFALSE; 

  TObjArray *clusters[200]={0};
  TObjArray *pix;

  // Find clusters of histogram bins (easier to work in 2-D space)
  for (Int_t i=1; i<=ny; i++) {
    for (Int_t j=1; j<=nx; j++) {
      indx = (i-1)*nx + j - 1;
      if (used[indx]) continue;
      cont = mlem->GetCellContent(j,i);
      if (cont < 0.5) continue;
      pix = new TObjArray(20);
      used[indx] = 1;
      pix->Add(BinToPix(mlem,j,i));
      AddBin(mlem, i, j, 0, used, pix); // recursive call
      clusters[nclust++] = pix;
      if (nclust > 200) { cout << " Too many clusters " << endl; ::exit(0); }
    } // for (Int_t j=1; j<=nx; j++) {
  } // for (Int_t i=1; i<=ny;
  cout << nclust << endl;
  delete [] used; used = 0;
  
  // Compute couplings between clusters and clusters to pads
  Int_t npad = fnPads[0] + fnPads[1];

  // Exclude pads with overflows
  for (Int_t j=0; j<npad; j++) {
    if (fXyq[2][j] > fResponse->MaxAdc()-1) fPadIJ[1][j] = -9;
    else fPadIJ[1][j] = 0;
  }

  // Compute couplings of clusters to pads
  TMatrixD *aijclupad = new TMatrixD(nclust,npad);
  *aijclupad = 0;
  Int_t npxclu;
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    pix = clusters[iclust];
    npxclu = pix->GetEntriesFast();
    for (Int_t i=0; i<npxclu; i++) {
      indx = fPixArray->IndexOf(pix->UncheckedAt(i));
      for (Int_t j=0; j<npad; j++) {
	// Exclude overflows
	if (fPadIJ[1][j] < 0) continue;
	if (coef[j*nPix+indx] < fgkCouplMin) continue;
	(*aijclupad)(iclust,j) += coef[j*nPix+indx];
      }
    }
  }
  // Compute couplings between clusters
  TMatrixD *aijcluclu = new TMatrixD(nclust,nclust);
  *aijcluclu = 0;
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    for (Int_t j=0; j<npad; j++) {
      // Exclude overflows
      if (fPadIJ[1][j] < 0) continue;
      if ((*aijclupad)(iclust,j) < fgkCouplMin) continue;
      for (Int_t iclust1=iclust+1; iclust1<nclust; iclust1++) {
	if ((*aijclupad)(iclust1,j) < fgkCouplMin) continue;
	(*aijcluclu)(iclust,iclust1) += 
	  TMath::Sqrt ((*aijclupad)(iclust,j)*(*aijclupad)(iclust1,j));
      }
    }
  }
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    for (Int_t iclust1=iclust+1; iclust1<nclust; iclust1++) {
      (*aijcluclu)(iclust1,iclust) = (*aijcluclu)(iclust,iclust1);
    }
  }

  if (nclust > 1) aijcluclu->Print();

  // Find groups of coupled clusters
  used = new Bool_t[nclust];
  for (Int_t i=0; i<nclust; i++) used[i] = kFALSE;
  Int_t *clustNumb = new Int_t[nclust];
  Int_t nCoupled, nForFit, minGroup[3], clustFit[3], nfit = 0;
  Double_t parOk[8];

  for (Int_t igroup=0; igroup<nclust; igroup++) {
    if (used[igroup]) continue;
    used[igroup] = kTRUE;
    clustNumb[0] = igroup;
    nCoupled = 1;
    // Find group of coupled clusters
    AddCluster(igroup, nclust, aijcluclu, used, clustNumb, nCoupled); // recursive
    cout << " nCoupled: " << nCoupled << endl;
    for (Int_t i=0; i<nCoupled; i++) cout << clustNumb[i] << " "; cout << endl; 

    while (nCoupled > 0) {

      if (nCoupled < 4) {
	nForFit = nCoupled;
	for (Int_t i=0; i<nCoupled; i++) clustFit[i] = clustNumb[i];
      } else {
	// Too many coupled clusters to fit - try to decouple them
	// Find the lowest coupling of 1, 2, min(3,nLinks/2) pixels with 
	// all the others in the group 
	for (Int_t j=0; j<3; j++) minGroup[j] = -1;
	Double_t coupl = MinGroupCoupl(nCoupled, clustNumb, aijcluclu, minGroup);

	// Flag clusters for fit
	nForFit = 0;
	while (minGroup[nForFit] >= 0 && nForFit < 3) {
	  cout << clustNumb[minGroup[nForFit]] << " ";
	  clustFit[nForFit] = clustNumb[minGroup[nForFit]];
	  clustNumb[minGroup[nForFit]] -= 999;
	  nForFit++;
	}
	cout << nForFit << " " << coupl << endl;
      } // else

      // Select pads for fit. 
      if (SelectPad(nCoupled, nForFit, clustNumb, clustFit, aijclupad) < 3 && nCoupled > 1) {
	// Deselect pads
	for (Int_t j=0; j<npad; j++) if (TMath::Abs(fPadIJ[1][j]) == 1) fPadIJ[1][j] = 0;
	// Merge the failed cluster candidates (with too few pads to fit) with 
	// the one with the strongest coupling
	Merge(nForFit, nCoupled, clustNumb, clustFit, clusters, aijcluclu, aijclupad);
      } else {
	// Do the fit
	nfit = Fit(nForFit, clustFit, clusters, parOk);
      }

      // Subtract the fitted charges from pads with strong coupling and/or
      // return pads for further use
      UpdatePads(nfit, parOk);

      // Mark used pads
      for (Int_t j=0; j<npad; j++) {if (fPadIJ[1][j] == 1) fPadIJ[1][j] = -1;}

      // Sort the clusters (move to the right the used ones)
      Int_t beg = 0, end = nCoupled - 1;
      while (beg < end) {
        if (clustNumb[beg] >= 0) { beg++; continue; }
        for (Int_t j=end; j>beg; j--) {
          if (clustNumb[j] < 0) continue;
          end = j - 1;
          indx = clustNumb[beg];
          clustNumb[beg] = clustNumb[j];
          clustNumb[j] = indx;
          break;
        }
	beg++;
      }

      nCoupled -= nForFit;
      if (nCoupled > 3) {
	// Remove couplings of used clusters
	for (Int_t iclust=nCoupled; iclust<nCoupled+nForFit; iclust++) {
	  indx = clustNumb[iclust] + 999;
	  for (Int_t iclust1=0; iclust1<nCoupled; iclust1++) {
	    indx1 = clustNumb[iclust1];
	    (*aijcluclu)(indx,indx1) = (*aijcluclu)(indx1,indx) = 0;
	  }
	}

	// Update the remaining clusters couplings (exclude couplings from 
	// the used pads)
	for (Int_t j=0; j<npad; j++) {
	  if (fPadIJ[1][j] != -1) continue;
	  for (Int_t iclust=0; iclust<nCoupled; iclust++) {
	    indx = clustNumb[iclust];
	    if ((*aijclupad)(indx,j) < fgkCouplMin) continue;
	    for (Int_t iclust1=iclust+1; iclust1<nCoupled; iclust1++) {
	      indx1 = clustNumb[iclust1];
	      if ((*aijclupad)(indx1,j) < fgkCouplMin) continue;
	      // Check this
	      (*aijcluclu)(indx,indx1) -= 
		TMath::Sqrt ((*aijclupad)(indx,j)*(*aijclupad)(indx1,j));
	      (*aijcluclu)(indx1,indx) = (*aijcluclu)(indx,indx1);
	    }
	  }
	  fPadIJ[1][j] = -9;
	} // for (Int_t j=0; j<npad;
      } // if (nCoupled > 3)
    } // while (nCoupled > 0)
  } // for (Int_t igroup=0; igroup<nclust;

  //delete aij_clu; aij_clu = 0; delete aijclupad; aijclupad = 0;
  aijcluclu->Delete(); aijclupad->Delete();
  for (Int_t iclust=0; iclust<nclust; iclust++) {
    pix = clusters[iclust]; 
    pix->Clear();
    delete pix; pix = 0;
  }
  delete [] clustNumb; clustNumb = 0; delete [] used; used = 0;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddBin(TH2D *mlem, Int_t ic, Int_t jc, Int_t mode, Bool_t *used, TObjArray *pix)
{
  // Add a bin to the cluster

  Int_t nx = mlem->GetNbinsX();
  Int_t ny = mlem->GetNbinsY();
  Double_t cont1, cont = mlem->GetCellContent(jc,ic);
  AliMUONPixel *pixPtr = 0;

  for (Int_t i=TMath::Max(ic-1,1); i<=TMath::Min(ic+1,ny); i++) {
    for (Int_t j=TMath::Max(jc-1,1); j<=TMath::Min(jc+1,nx); j++) {
      if (i != ic && j != jc) continue;
      if (used[(i-1)*nx+j-1]) continue;
      cont1 = mlem->GetCellContent(j,i);
      if (mode && cont1 > cont) continue;
      used[(i-1)*nx+j-1] = kTRUE;
      if (cont1 < 0.5) continue;
      if (pix) pix->Add(BinToPix(mlem,j,i)); 
      else {
	pixPtr = new AliMUONPixel (mlem->GetXaxis()->GetBinCenter(j), 
				   mlem->GetYaxis()->GetBinCenter(i), 0, 0, cont1);
	fPixArray->Add((TObject*)pixPtr);
      }
      AddBin(mlem, i, j, mode, used, pix); // recursive call
    }
  }
}

//_____________________________________________________________________________
TObject* AliMUONClusterFinderAZ::BinToPix(TH2D *mlem, Int_t jc, Int_t ic)
{
  // Translate histogram bin to pixel 
  
  Double_t yc = mlem->GetYaxis()->GetBinCenter(ic);
  Double_t xc = mlem->GetXaxis()->GetBinCenter(jc);
  
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr;

  // Compare pixel and bin positions
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    if (pixPtr->Charge() < 0.5) continue;
    if (TMath::Abs(pixPtr->Coord(0)-xc)<1.e-4 && TMath::Abs(pixPtr->Coord(1)-yc)<1.e-4) return (TObject*) pixPtr;
  }
  cout << " Something wrong ??? " << endl;
  return NULL;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddCluster(Int_t ic, Int_t nclust, TMatrixD *aijcluclu, Bool_t *used, Int_t *clustNumb, Int_t &nCoupled)
{
  // Add a cluster to the group of coupled clusters

  for (Int_t i=0; i<nclust; i++) {
    if (used[i]) continue;
    if ((*aijcluclu)(i,ic) < fgkCouplMin) continue;
    used[i] = kTRUE;
    clustNumb[nCoupled++] = i;
    AddCluster(i, nclust, aijcluclu, used, clustNumb, nCoupled);
  }
}

//_____________________________________________________________________________
Double_t AliMUONClusterFinderAZ::MinGroupCoupl(Int_t nCoupled, Int_t *clustNumb, TMatrixD *aijcluclu, Int_t *minGroup)
{
  // Find group of clusters with minimum coupling to all the others

  Int_t i123max = TMath::Min(3,nCoupled/2); 
  Int_t indx, indx1, indx2, indx3, nTot = 0;
  Double_t *coupl1 = 0, *coupl2 = 0, *coupl3 = 0;

  for (Int_t i123=1; i123<=i123max; i123++) {

    if (i123 == 1) {
      coupl1 = new Double_t [nCoupled];
      for (Int_t i=0; i<nCoupled; i++) coupl1[i] = 0;
    }
    else if (i123 == 2) {
      nTot = nCoupled*nCoupled;
      coupl2 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl2[i] = 9999;
    } else {
      nTot = nTot*nCoupled;
      coupl3 = new Double_t [nTot];
      for (Int_t i=0; i<nTot; i++) coupl3[i] = 9999;
    } // else

    for (Int_t i=0; i<nCoupled; i++) {
      indx1 = clustNumb[i];
      for (Int_t j=i+1; j<nCoupled; j++) {
	indx2 = clustNumb[j];
	if (i123 == 1) {
	  coupl1[i] += (*aijcluclu)(indx1,indx2);
	  coupl1[j] += (*aijcluclu)(indx1,indx2);
	} 
	else if (i123 == 2) {
	  indx = i*nCoupled + j;
	  coupl2[indx] = coupl1[i] + coupl1[j];
	  coupl2[indx] -= 2 * ((*aijcluclu)(indx1,indx2));
	} else {
	  for (Int_t k=j+1; k<nCoupled; k++) {
	    indx3 = clustNumb[k];
	    indx = i*nCoupled*nCoupled + j*nCoupled + k;
	    coupl3[indx] = coupl2[i*nCoupled+j] + coupl1[k];
	    coupl3[indx] -= 2 * ((*aijcluclu)(indx1,indx3)+(*aijcluclu)(indx2,indx3));
	  }
	} // else
      } // for (Int_t j=i+1;
    } // for (Int_t i=0;
  } // for (Int_t i123=1;

  // Find minimum coupling
  Double_t couplMin = 9999;
  Int_t locMin = 0;

  for (Int_t i123=1; i123<=i123max; i123++) {
    if (i123 == 1) {
      locMin = TMath::LocMin(nCoupled, coupl1);
      couplMin = coupl1[locMin];
      minGroup[0] = locMin;
      delete [] coupl1; coupl1 = 0;
    } 
    else if (i123 == 2) {
      locMin = TMath::LocMin(nCoupled*nCoupled, coupl2);
      if (coupl2[locMin] < couplMin) {
	couplMin = coupl2[locMin];
	minGroup[0] = locMin/nCoupled;
	minGroup[1] = locMin%nCoupled;
      }
      delete [] coupl2; coupl2 = 0;
    } else {
      locMin = TMath::LocMin(nTot, coupl3);
      if (coupl3[locMin] < couplMin) {
	couplMin = coupl3[locMin];
	minGroup[0] = locMin/nCoupled/nCoupled;
	minGroup[1] = locMin%(nCoupled*nCoupled)/nCoupled;
	minGroup[2] = locMin%nCoupled;
      }
      delete [] coupl3; coupl3 = 0;
    } // else
  } // for (Int_t i123=1;
  return couplMin;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::SelectPad(Int_t nCoupled, Int_t nForFit, Int_t *clustNumb, Int_t *clustFit, TMatrixD *aijclupad)
{
  // Select pads for fit. If too many coupled clusters, find pads giving 
  // the strongest coupling with the rest of clusters and exclude them from the fit.

  Int_t npad = fnPads[0] + fnPads[1];
  Double_t *padpix = 0;

  if (nCoupled > 3) {
    padpix = new Double_t[npad];
    for (Int_t i=0; i<npad; i++) padpix[i] = 0; 
  }

  Int_t nOK = 0, indx, indx1;
  for (Int_t iclust=0; iclust<nForFit; iclust++) {
    indx = clustFit[iclust];
    for (Int_t j=0; j<npad; j++) {
      if (fPadIJ[1][j] < 0) continue; // exclude overflows and used pads
      if ((*aijclupad)(indx,j) < fgkCouplMin) continue;
      fPadIJ[1][j] = 1; // pad to be used in fit
      nOK++;
      if (nCoupled > 3) {
	// Check other clusters
	for (Int_t iclust1=0; iclust1<nCoupled; iclust1++) {
	  indx1 = clustNumb[iclust1];
	  if (indx1 < 0) continue;
	  if ((*aijclupad)(indx1,j) < fgkCouplMin) continue;
	  padpix[j] += (*aijclupad)(indx1,j);
	}
      } // if (nCoupled > 3)
    } // for (Int_t j=0; j<npad;
  } // for (Int_t iclust=0; iclust<nForFit
  if (nCoupled < 4) return nOK;

  Double_t aaa = 0;
  for (Int_t j=0; j<npad; j++) {
    if (padpix[j] < fgkCouplMin) continue;
    cout << j << " " << padpix[j] << " "; 
    cout << fXyq[0][j] << " " << fXyq[1][j] << endl;
    aaa += padpix[j];
    fPadIJ[1][j] = -1; // exclude pads with strong coupling to the other clusters
    nOK--;
  }
  delete [] padpix; padpix = 0;
  return nOK;
}
  
//_____________________________________________________________________________
void AliMUONClusterFinderAZ::Merge(Int_t nForFit, Int_t nCoupled, Int_t *clustNumb, Int_t *clustFit, TObjArray **clusters, TMatrixD *aijcluclu, TMatrixD *aijclupad)
{
  // Merge the group of clusters with the one having the strongest coupling with them

  Int_t indx, indx1, npxclu, npxclu1, imax=0;
  TObjArray *pix, *pix1;
  Double_t couplMax;

  for (Int_t icl=0; icl<nForFit; icl++) {
    indx = clustFit[icl];
    pix = clusters[indx];
    npxclu = pix->GetEntriesFast();
    couplMax = -1;
    for (Int_t icl1=0; icl1<nCoupled; icl1++) {
      indx1 = clustNumb[icl1];
      if (indx1 < 0) continue;
      if ((*aijcluclu)(indx,indx1) > couplMax) {
	couplMax = (*aijcluclu)(indx,indx1);
	imax = indx1;
      }
    } // for (Int_t icl1=0;
    /*if (couplMax < fgkCouplMin) {
      cout << " Oops " << couplMax << endl;
      aijcluclu->Print();
      cout << icl << " " << indx << " " << npxclu << " " << nLinks << endl;
      ::exit(0);
      }*/
    // Add to it
    pix1 = clusters[imax];
    npxclu1 = pix1->GetEntriesFast();
    // Add pixels 
    for (Int_t i=0; i<npxclu; i++) { pix1->Add(pix->UncheckedAt(i)); pix->RemoveAt(i); }
    cout << " New number of pixels: " << npxclu1 << " " << pix1->GetEntriesFast() << endl;
    //Add cluster-to-cluster couplings
    //aijcluclu->Print();
    for (Int_t icl1=0; icl1<nCoupled; icl1++) {
      indx1 = clustNumb[icl1];
      if (indx1 < 0 || indx1 == imax) continue;
      (*aijcluclu)(indx1,imax) += (*aijcluclu)(indx,indx1);
      (*aijcluclu)(imax,indx1) = (*aijcluclu)(indx1,imax);
    }
    (*aijcluclu)(indx,imax) = (*aijcluclu)(imax,indx) = 0;
    //aijcluclu->Print();
    //Add cluster-to-pad couplings
    for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
      if (fPadIJ[1][j] < 0) continue; // exclude overflows and used pads
      (*aijclupad)(imax,j) += (*aijclupad)(indx,j);
      (*aijclupad)(indx,j) = 0;
    }
  } // for (Int_t icl=0; icl<nForFit;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::Fit(Int_t nfit, Int_t *clustFit, TObjArray **clusters, Double_t *parOk)
{
  // Find selected clusters to selected pad charges
  
  TH2D *mlem = (TH2D*) gROOT->FindObject("mlem");
  //Int_t nx = mlem->GetNbinsX();
  //Int_t ny = mlem->GetNbinsY();
  Double_t xmin = mlem->GetXaxis()->GetXmin() - mlem->GetXaxis()->GetBinWidth(1);
  Double_t xmax = mlem->GetXaxis()->GetXmax() + mlem->GetXaxis()->GetBinWidth(1);
  Double_t ymin = mlem->GetYaxis()->GetXmin() - mlem->GetYaxis()->GetBinWidth(1);
  Double_t ymax = mlem->GetYaxis()->GetXmax() + mlem->GetYaxis()->GetBinWidth(1);
  //Double_t qmin = 0, qmax = 1;
  Double_t step[3]={0.01,0.002,0.02};

  Double_t cont, cmax = 0, xseed = 0, yseed = 0, errOk[8];
  TObjArray *pix;
  Int_t npxclu;

  // Number of pads to use
  Int_t npads = 0;
  for (Int_t i=0; i<fnPads[0]+fnPads[1]; i++) {if (fPadIJ[1][i] == 1) npads++;}
  for (Int_t i=0; i<nfit; i++) {cout << i+1 << " " << clustFit[i] << " ";}
  cout << nfit << endl;
  cout << " Number of pads to fit: " << npads << endl;
  fNpar = 0;
  fQtot = 0;
  if (npads < 2) return 0; 

  // Take cluster maxima as fitting seeds
  AliMUONPixel *pixPtr;
  Double_t xyseed[3][2], qseed[3];
  for (Int_t ifit=1; ifit<=nfit; ifit++) {
    cmax = 0;
    pix = clusters[clustFit[ifit-1]];
    npxclu = pix->GetEntriesFast();
    for (Int_t clu=0; clu<npxclu; clu++) {
      pixPtr = (AliMUONPixel*) pix->UncheckedAt(clu);
      cont = pixPtr->Charge();
      fQtot += cont;
      if (cont > cmax) { 
	cmax = cont; 
	xseed = pixPtr->Coord(0);
	yseed = pixPtr->Coord(1);
      }
    }
    xyseed[ifit-1][0] = xseed;
    xyseed[ifit-1][1] = yseed;
    qseed[ifit-1] = cmax;
  } // for (Int_t ifit=1;

  Int_t nDof, maxSeed[3];
  Double_t fmin, chi2o = 9999, chi2n;

  // Try to fit with one-track hypothesis, then 2-track. If chi2/dof is 
  // lower, try 3-track (if number of pads is sufficient).
  
  TMath::Sort(nfit, qseed, maxSeed, kTRUE); // in decreasing order
  nfit = TMath::Min (nfit, (npads + 1) / 3);

  Double_t *gin = 0, func0, func1, param[8], param0[2][8], deriv[2][8], step0[8];
  Double_t shift[8], stepMax, derMax, parmin[8], parmax[8], func2[2], shift0;
  Double_t delta[8], scMax, dder[8], estim, shiftSave = 0;
  Int_t min, max, nCall = 0, memory[8] = {0}, nLoop, idMax = 0, iestMax = 0, nFail;

  for (Int_t iseed=0; iseed<nfit; iseed++) {

    for (Int_t j=0; j<3; j++) step0[fNpar+j] = shift[fNpar+j] = step[j];
    param[fNpar] = xyseed[maxSeed[iseed]][0];
    parmin[fNpar] = xmin; 
    parmax[fNpar++] = xmax; 
    param[fNpar] = xyseed[maxSeed[iseed]][1];
    parmin[fNpar] = ymin; 
    parmax[fNpar++] = ymax; 
    if (fNpar > 2) {
      param[fNpar] = fNpar == 4 ? 0.5 : 0.3;
      parmin[fNpar] = 0; 
      parmax[fNpar++] = 1; 
    }

    // Try new algorithm
    min = nLoop = 1; stepMax = func2[1] = derMax = 999999; nFail = 0;

    while (1) {
      max = !min;
      fcn1(fNpar, gin, func0, param, 1); nCall++;
      //cout << " Func: " << func0 << endl;

      func2[max] = func0;
      for (Int_t j=0; j<fNpar; j++) {
	param0[max][j] = param[j];
	delta[j] = step0[j];
	param[j] += delta[j] / 10;
	if (j > 0) param[j-1] -= delta[j-1] / 10;
	fcn1(fNpar, gin, func1, param, 1); nCall++;
	deriv[max][j] = (func1 - func0) / delta[j] * 10; // first derivative
	//cout << j << " " << deriv[max][j] << endl;
	dder[j] = param0[0][j] != param0[1][j] ? (deriv[0][j] - deriv[1][j]) / 
	                                         (param0[0][j] - param0[1][j]) : 0; // second derivative
      }
      param[fNpar-1] -= delta[fNpar-1] / 10;
      if (nCall > 2000) ::exit(0);

      min = func2[0] < func2[1] ? 0 : 1;
      nFail = min == max ? 0 : nFail + 1;

      stepMax = derMax = estim = 0;
      for (Int_t j=0; j<fNpar; j++) { 
	// Estimated distance to minimum
	shift0 = shift[j];
	if (nLoop == 1) shift[j] = TMath::Sign (step0[j], -deriv[max][j]); // first step
	else if (TMath::Abs(deriv[0][j]) < 1.e-3 && TMath::Abs(deriv[1][j]) < 1.e-3) shift[j] = 0;
	else if (deriv[min][j]*deriv[!min][j] > 0 && TMath::Abs(deriv[min][j]) > TMath::Abs(deriv[!min][j])
	      || TMath::Abs(deriv[0][j]-deriv[1][j]) < 1.e-3) {
	  shift[j] = -TMath::Sign (shift[j], (func2[0]-func2[1]) * (param0[0][j]-param0[1][j]));
	  if (min == max) { 
	    if (memory[j] > 1) { shift[j] *= 2; } //cout << " Memory " << memory[j] << " " << shift[j] << endl; }
	    memory[j]++;
	  }
	} else {
	  shift[j] = -deriv[min][j] / dder[j];
	  memory[j] = 0;
	}
	if (TMath::Abs(shift[j])/step0[j] > estim) { 
	  estim = TMath::Abs(shift[j])/step0[j];
	  iestMax = j;
	}

	// Too big step
	if (TMath::Abs(shift[j])/step0[j] > 10) shift[j] = TMath::Sign(10.,shift[j]) * step0[j]; // 

	// Failed to improve minimum
	if (min != max) {
	  memory[j] = 0;
	  param[j] = param0[min][j];
	  if (TMath::Abs(shift[j]+shift0) > 0.1*step0[j]) shift[j] = (shift[j] + shift0) / 2;
	  else shift[j] /= -2;
	} 

	// Too big step
	if (TMath::Abs(shift[j]*deriv[min][j]) > func2[min]) 
	  shift[j] = TMath::Sign (func2[min]/deriv[min][j], shift[j]);

	// Introduce step relaxation factor
	if (memory[j] < 3) {
	  scMax = 1 + 4 / TMath::Max(nLoop/2.,1.);
	  if (TMath::Abs(shift0) > 0 && TMath::Abs(shift[j]/shift0) > scMax) 
	    shift[j] = TMath::Sign (shift0*scMax, shift[j]);
	}
	param[j] += shift[j]; 
	  
	//cout << " xxx " << j << " " << shift[j] << " " << param[j] << endl;
	stepMax = TMath::Max (stepMax, TMath::Abs(shift[j]/step0[j]));
	if (TMath::Abs(deriv[min][j]) > derMax) {
	  idMax = j;
	  derMax = TMath::Abs (deriv[min][j]);
	}
      } // for (Int_t j=0; j<fNpar;
      //cout << max << " " << func2[min] << " " << derMax << " " << stepMax << " " << estim << " " << iestMax << " " << nCall << endl;
      if (estim < 1 && derMax < 2 || nLoop > 100) break; // minimum was found

      nLoop++;
      // Check for small step
      if (shift[idMax] == 0) { shift[idMax] = step0[idMax]/10; param[idMax] += shift[idMax]; continue; }
      if (!memory[idMax] && derMax > 0.5 && nLoop > 10) {
	//cout << " ok " << deriv[min][idMax] << " " << deriv[!min][idMax] << " " << dder[idMax]*shift[idMax] << " " << shift[idMax] << endl;
	if (dder[idMax] != 0 && TMath::Abs(deriv[min][idMax]/dder[idMax]/shift[idMax]) > 10) {
	  if (min == max) dder[idMax] = -dder[idMax];
	  shift[idMax] = -deriv[min][idMax] / dder[idMax] / 10; 
	  param[idMax] += shift[idMax];
	  stepMax = TMath::Max (stepMax, TMath::Abs(shift[idMax])/step0[idMax]);
	  //cout << shift[idMax] << " " << param[idMax] << endl;
	  if (min == max) shiftSave = shift[idMax];
	}
	if (nFail > 10) {
	  param[idMax] -= shift[idMax];
	  shift[idMax] = 4 * shiftSave * (gRandom->Rndm(0) - 0.5);
	  param[idMax] += shift[idMax];
	  //cout << shift[idMax] << endl;
	}
      }      
    } // while (1)
    fmin = func2[min];

    nDof = npads - fNpar;
    chi2n = nDof ? fmin/nDof : 0;

    if (chi2n*1.2+1.e-6 > chi2o ) { fNpar -= 3; break; }
    // Save parameters and errors
    for (Int_t i=0; i<fNpar; i++) {
      parOk[i] = param0[min][i];
      errOk[i] = fmin;
    }

    cout << chi2o << " " << chi2n << endl;
    chi2o = chi2n;
    if (fmin < 0.1) break; // !!!???
  } // for (Int_t iseed=0; 

  for (Int_t i=0; i<fNpar; i++) {
    if (i == 4 || i == 7) continue;
    cout << parOk[i] << " " << errOk[i] << endl;
  }
  nfit = (fNpar + 1) / 3;
  Double_t rad;
  Int_t indx, imax;
  if (fReco) {
    for (Int_t j=0; j<nfit; j++) {
      indx = j<2 ? j*2 : j*2+1;  
      AddRawCluster (parOk[indx], parOk[indx+1], errOk[indx]);
    }
    return nfit;
  } 
  for (Int_t i=0; i<fnMu; i++) {
    cmax = fxyMu[i][6];
    for (Int_t j=0; j<nfit; j++) {
      indx = j<2 ? j*2 : j*2+1;  
      rad = (fxyMu[i][0]-parOk[indx])*(fxyMu[i][0]-parOk[indx]) +
            (fxyMu[i][1]-parOk[indx+1])*(fxyMu[i][1]-parOk[indx+1]);
      if (rad < cmax) {
	cmax = rad; 
	imax = indx;
	fxyMu[i][6] = cmax;
	fxyMu[i][2] = parOk[imax] - fxyMu[i][0];
	fxyMu[i][4] = parOk[imax+1] - fxyMu[i][1];
	fxyMu[i][3] = errOk[imax];
	fxyMu[i][5] = errOk[imax+1];
      }
    }      
  }
  return nfit;
}  

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::fcn1(Int_t & /*npar*/, Double_t * /*gin*/, Double_t &f, Double_t *par, Int_t /*iflag*/)
{
  // Fit for one track
  AliMUONClusterFinderAZ& c = *(AliMUONClusterFinderAZ::fgClusterFinder);    
  
  Int_t cath, ix, iy, indx, npads=0;
  Double_t charge, delta, coef=0, chi2=0;
  for (Int_t j=0; j<c.fnPads[0]+c.fnPads[1]; j++) {
    if (c.fPadIJ[1][j] != 1) continue;
    cath = c.fPadIJ[0][j];
    npads++;
    c.fSegmentation[cath]->GetPadI(c.fXyq[0][j],c.fXyq[1][j],c.fZpad,ix,iy);
    c.fSegmentation[cath]->SetPad(ix,iy);
    charge = 0;
    for (Int_t i=c.fNpar/3; i>=0; i--) { // sum over tracks
      indx = i<2 ? 2*i : 2*i+1;
      c.fSegmentation[cath]->SetHit(par[indx],par[indx+1],c.fZpad);
      //charge += c.fResponse->IntXY(c.fSegmentation[cath])*par[icl*3+2];
      if (c.fNpar == 2) coef = 1;
      else coef = i==c.fNpar/3 ? par[indx+2] : 1-coef;
      //coef = TMath::Max (coef, 0.);
      if (c.fNpar == 8 && i < 2) coef = i==1 ? coef*par[indx+2] : coef - par[7];
      //coef = TMath::Max (coef, 0.);
      charge += c.fResponse->IntXY(c.fSegmentation[cath])*coef;
    }
    charge *= c.fQtot;
    //if (c.fXyq[2][j] > c.fResponse->MaxAdc()-1 && charge > 
    //	c.fResponse->MaxAdc()) charge = c.fResponse->MaxAdc(); 
    delta = charge - c.fXyq[2][j];
    delta /= TMath::Sqrt ((Double_t)c.fXyq[2][j]);
    //chi2 += TMath::Abs(delta);
    chi2 += delta*delta;
  } // for (Int_t j=0;
  f = chi2; 
  Double_t qAver = c.fQtot/npads; //(c.fnPads[0]+c.fnPads[1]);
  f = chi2/qAver;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::UpdatePads(Int_t /*nfit*/, Double_t *par)
{
  // Subtract the fitted charges from pads with strong coupling

  Int_t cath, ix, iy, indx;
  Double_t charge, coef=0;
  for (Int_t j=0; j<fnPads[0]+fnPads[1]; j++) {
    if (fPadIJ[1][j] != -1) continue;
    if (fNpar != 0) {
      cath = fPadIJ[0][j];
      fSegmentation[cath]->GetPadI(fXyq[0][j],fXyq[1][j],fZpad,ix,iy);
      fSegmentation[cath]->SetPad(ix,iy);
      charge = 0;
      for (Int_t i=fNpar/3; i>=0; i--) { // sum over tracks
	indx = i<2 ? 2*i : 2*i+1;
	fSegmentation[cath]->SetHit(par[indx],par[indx+1],fZpad);
	if (fNpar == 2) coef = 1;
	else coef = i==fNpar/3 ? par[indx+2] : 1-coef;
	if (fNpar == 8 && i < 2) coef = i==1 ? coef*par[indx+2] : coef - par[7];
	charge += fResponse->IntXY(fSegmentation[cath])*coef;
      }
      charge *= fQtot;
      fXyq[2][j] -= charge;
    } // if (fNpar != 0)
    if (fXyq[2][j] > fResponse->ZeroSuppression()) fPadIJ[1][j] = 0; // return pad for further using
  } // for (Int_t j=0;
}  

//_____________________________________________________________________________
Bool_t AliMUONClusterFinderAZ::TestTrack(Int_t /*t*/) {
// Test if track was user selected
  return kTRUE;
  /*
    if (fTrack[0]==-1 || fTrack[1]==-1) {
	return kTRUE;
    } else if (t==fTrack[0] || t==fTrack[1]) {
	return kTRUE;
    } else {
	return kFALSE;
    }
  */
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::AddRawCluster(Double_t x, Double_t y, Double_t fmin)
{
  //
  // Add a raw cluster copy to the list
  //
  AliMUONRawCluster cnew;
  AliMUON *pMUON=(AliMUON*)gAlice->GetModule("MUON");
  //pMUON->AddRawCluster(fInput->Chamber(),c); 

  Int_t cath;    
  for (cath=0; cath<2; cath++) {
    cnew.SetX(cath, x);
    cnew.SetY(cath, y);
    cnew.SetZ(cath, fZpad);
    cnew.SetCharge(cath, 100);
    cnew.SetPeakSignal(cath,20);
    cnew.SetMultiplicity(cath, 5);
    cnew.SetNcluster(cath, 1);
    cnew.SetChi2(cath, fmin); //0.1;
    /*
    cnew.fMultiplicity[cath]=c->fMultiplicity[cath];
    for (i=0; i<fMul[cath]; i++) {
      cnew.fIndexMap[i][cath]=c->fIndexMap[i][cath];
      fSeg[cath]->SetPad(fIx[i][cath], fIy[i][cath]);
    }
    fprintf(stderr,"\nRawCluster %d cath %d\n",ico,cath);
    fprintf(stderr,"mult_av %d\n",c->fMultiplicity[cath]);
    FillCluster(&cnew,cath);
    */
  } 
  //cnew.fClusterType=cnew.PhysicsContribution();
  pMUON->GetMUONData()->AddRawCluster(AliMUONClusterInput::Instance()->Chamber(),cnew); 
  //fNPeaks++;
}

//_____________________________________________________________________________
Int_t AliMUONClusterFinderAZ::FindLocalMaxima(Int_t *localMax, Double_t *maxVal)
{
  // Find local maxima in pixel space for large preclusters in order to
  // try to split them into smaller pieces (to speed up the MLEM procedure)

  TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  if (hist) hist->Delete();

  Double_t xylim[4] = {999, 999, 999, 999};
  Int_t nPix = fPixArray->GetEntriesFast();
  AliMUONPixel *pixPtr = 0;
  for (Int_t ipix=0; ipix<nPix; ipix++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
    for (Int_t i=0; i<4; i++) 
         xylim[i] = TMath::Min (xylim[i], (i%2 ? -1 : 1)*pixPtr->Coord(i/2));
  }
  for (Int_t i=0; i<4; i++) xylim[i] -= pixPtr->Size(i/2); 

  Int_t nx = TMath::Nint ((-xylim[1]-xylim[0])/pixPtr->Size(0)/2);
  Int_t ny = TMath::Nint ((-xylim[3]-xylim[2])/pixPtr->Size(1)/2);
  hist = new TH2D("anode","anode",nx,xylim[0],-xylim[1],ny,xylim[2],-xylim[3]);
  for (Int_t ipix=0; ipix<nPix; ipix++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(ipix);
    hist->Fill(pixPtr->Coord(0), pixPtr->Coord(1), pixPtr->Charge());
  }
  if (fDraw) {
    ((TCanvas*)gROOT->FindObject("c2"))->cd();
    gPad->SetTheta(55);
    gPad->SetPhi(30);
    hist->Draw("lego1Fb");
    gPad->Update();
    int ia;
    cin >> ia;
  }

  Int_t nMax = 0, indx;
  Int_t *isLocalMax = new Int_t[ny*nx];
  for (Int_t i=0; i<ny*nx; i++) isLocalMax[i] = 0;

  for (Int_t i=1; i<=ny; i++) {
    indx = (i-1) * nx;
    for (Int_t j=1; j<=nx; j++) {
      if (hist->GetCellContent(j,i) < 0.5) continue;
      //if (isLocalMax[indx+j-1] < 0) continue;
      if (isLocalMax[indx+j-1] != 0) continue;
      FlagLocalMax(hist, i, j, isLocalMax);
    }
  }

  for (Int_t i=1; i<=ny; i++) {
    indx = (i-1) * nx;
    for (Int_t j=1; j<=nx; j++) {
      if (isLocalMax[indx+j-1] > 0) { 
	localMax[nMax] = indx + j - 1; 
	maxVal[nMax++] = hist->GetCellContent(j,i);
      }
      if (nMax > 99) { cout << " Too many local maxima !!!" << endl; ::exit(0); }
    }
  }
  cout << " Local max: " << nMax << endl;
  delete [] isLocalMax; isLocalMax = 0;
  return nMax;
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FlagLocalMax(TH2D *hist, Int_t i, Int_t j, Int_t *isLocalMax)
{
  // Flag pixels (whether or not local maxima)

  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t cont = TMath::Nint (hist->GetCellContent(j,i));
  Int_t cont1 = 0;

  for (Int_t i1=i-1; i1<i+2; i1++) {
    if (i1 < 1 || i1 > ny) continue;
    for (Int_t j1=j-1; j1<j+2; j1++) {
      if (j1 < 1 || j1 > nx) continue;
      if (i == i1 && j == j1) continue;
      cont1 = TMath::Nint (hist->GetCellContent(j1,i1));
      if (cont < cont1) { isLocalMax[(i-1)*nx+j-1] = -1; return; }
      else if (cont > cont1) isLocalMax[(i1-1)*nx+j1-1] = -1;
      else { // the same charge
	isLocalMax[(i-1)*nx+j-1] = 1; 
	if (isLocalMax[(i1-1)*nx+j1-1] == 0) {
	  FlagLocalMax(hist, i1, j1, isLocalMax);
	  if (isLocalMax[(i1-1)*nx+j1-1] < 0) { isLocalMax[(i-1)*nx+j-1] = -1; return; }
	  else isLocalMax[(i1-1)*nx+j1-1] = -1;
	}
      } 
    }
  }
  isLocalMax[(i-1)*nx+j-1] = 1; // local maximum
}

//_____________________________________________________________________________
void AliMUONClusterFinderAZ::FindCluster(Int_t *localMax, Int_t iMax)
{
  // Find pixel cluster around local maximum #iMax and pick up pads
  // overlapping with it

  TH2D *hist = (TH2D*) gROOT->FindObject("anode");
  Int_t nx = hist->GetNbinsX();
  Int_t ny = hist->GetNbinsY();
  Int_t ic = localMax[iMax] / nx + 1;
  Int_t jc = localMax[iMax] % nx + 1;
  Bool_t *used = new Bool_t[ny*nx];
  for (Int_t i=0; i<ny*nx; i++) used[i] = kFALSE;

  // Drop all pixels from the array - pick up only the ones from the cluster
  fPixArray->Delete();

  Double_t wx = hist->GetXaxis()->GetBinWidth(1)/2; 
  Double_t wy = hist->GetYaxis()->GetBinWidth(1)/2;  
  Double_t yc = hist->GetYaxis()->GetBinCenter(ic);
  Double_t xc = hist->GetXaxis()->GetBinCenter(jc);
  Double_t cont = hist->GetCellContent(jc,ic);
  AliMUONPixel *pixPtr = new AliMUONPixel (xc, yc, wx, wy, cont);
  fPixArray->Add((TObject*)pixPtr);
  used[(ic-1)*nx+jc-1] = kTRUE;
  AddBin(hist, ic, jc, 1, used, (TObjArray*)0); // recursive call

  Int_t nPix = fPixArray->GetEntriesFast(), npad = fnPads[0] + fnPads[1];
  for (Int_t i=0; i<nPix; i++) {
    ((AliMUONPixel*)fPixArray->UncheckedAt(i))->SetSize(0,wx); 
    ((AliMUONPixel*)fPixArray->UncheckedAt(i))->SetSize(1,wy); 
  }
  cout << iMax << " " << nPix << endl;

  Float_t xy[4], xy12[4];
  // Pick up pads which overlap with found pixels
  for (Int_t i=0; i<npad; i++) fPadIJ[1][i] = -1;
  for (Int_t i=0; i<nPix; i++) {
    pixPtr = (AliMUONPixel*) fPixArray->UncheckedAt(i);
    for (Int_t j=0; j<4; j++) 
      xy[j] = pixPtr->Coord(j/2) + (j%2 ? 1 : -1)*pixPtr->Size(j/2);
    for (Int_t j=0; j<npad; j++) 
      if (Overlap(xy, j, xy12, 0)) fPadIJ[1][j] = 0; // flag for use
  }

  delete [] used; used = 0;
}
