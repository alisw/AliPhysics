#include "AliDisplacedVertexSelection.h"
// #include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMultiplicity.h"
// #include "AliMCEventHandler.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include <iostream>
#include <TROOT.h>
#include <TH1D.h>
#include <TH2D.h>
#include "AliESDEvent.h"
#include "AliESDZDC.h"
ClassImp(AliDisplacedVertexSelection)
#if 0
; // This is for Emacs
#endif

//____________________________________________________________________
AliDisplacedVertexSelection::AliDisplacedVertexSelection()
  : TObject(),
    fDeltaTdc(0),
    fSumTdc(0),
    fZdcEnergy(0),
    fZemEnergy(0),
    fCorrelationZemZdc(0),
    fCorrelationSumDelta(0), 
    fVertexZ(kInvalidVtxZ), 
    fCent(100), 
    fHVertexZ(0),
    fHCent(0),
    fMC(kFALSE)
{
}
//____________________________________________________________________
AliDisplacedVertexSelection::AliDisplacedVertexSelection(const AliDisplacedVertexSelection& o)
  : TObject(o),
    fDeltaTdc(o.fDeltaTdc),
    fSumTdc(o.fSumTdc),
    fZdcEnergy(o.fZdcEnergy),
    fZemEnergy(o.fZemEnergy),
    fCorrelationZemZdc(o.fCorrelationZemZdc),
    fCorrelationSumDelta(o.fCorrelationSumDelta), 
    fVertexZ(kInvalidVtxZ), 
    fCent(100), 
    fHVertexZ(0),
    fHCent(0),
    fMC(kFALSE)
{
}
//____________________________________________________________________
AliDisplacedVertexSelection&
AliDisplacedVertexSelection::operator=(const AliDisplacedVertexSelection& o)
{
  if (&o == this) return *this;
  
  fDeltaTdc  = o.fDeltaTdc;
  fSumTdc    = o.fSumTdc;
  fZdcEnergy = o.fZdcEnergy;
  fZemEnergy = o.fZemEnergy;
  fCorrelationZemZdc = o.fCorrelationZemZdc;
  fCorrelationSumDelta = o.fCorrelationSumDelta;
  fMC = o.fMC;
  
  return *this;
}

//____________________________________________________________________
void
AliDisplacedVertexSelection::SetupForData(TList* l, 
					  const char* /* name*/,
					  Bool_t mc)
{
  fMC = mc;
  
  TList* out = new TList;
  out->SetName("displacedVertex");
  out->SetOwner();
  l->Add(out);

  Double_t dVz   = 37.5;
  Double_t vzMin = (-kMaxK-.5) * dVz;
  Double_t vzMax = (+kMaxK+.5) * dVz;

  fHVertexZ = new TH1D("vertexZ", "Interaction point Z", 
		       2*kMaxK+1, vzMin, vzMax);
  fHVertexZ->SetXTitle("IP_{z} [cm]");
  fHVertexZ->SetYTitle("events");
  fHVertexZ->SetDirectory(0);
  fHVertexZ->SetFillColor(kRed+1);
  fHVertexZ->SetFillStyle(3001);
  out->Add(fHVertexZ);

  Int_t    nCent   = 6;
  Double_t bCent[] = { 0, 5, 10, 20, 30, 40, 100 };
  fHCent = new TH1D("cent", "Centrality", nCent, bCent);
  fHCent->SetXTitle("Centrality [%]");
  fHCent->SetYTitle("events");
  fHCent->SetDirectory(0);
  fHCent->SetFillColor(kBlue+1);
  fHCent->SetFillStyle(3001);
  out->Add(fHCent);

  Int_t    nDeltaT   = 2000;
  Double_t minDeltaT = -250;
  Double_t maxDeltaT = +250;
  fDeltaTdc          = new TH1D("DeltaTdc","#DeltaTDC",
				nDeltaT,minDeltaT,maxDeltaT);
  fDeltaTdc->SetXTitle("#DeltaTDC");
  fDeltaTdc->SetDirectory(0);
  out->Add(fDeltaTdc);

  Int_t    nSumT   = nDeltaT;
  Double_t minSumT = minDeltaT;
  Double_t maxSumT = maxDeltaT;
  fSumTdc          = new TH1D("SumTdc","#sumTDC",nSumT,minSumT,maxSumT);
  fSumTdc->SetXTitle("#sumTDC");
  fSumTdc->SetDirectory(0);
  out->Add(fSumTdc);
  
  Int_t    nZdcE   = 10000;
  Double_t minZdcE = 0;
  Double_t maxZdcE = 200000;
  fZdcEnergy           = new TH1D("ZDCEnergy","E_{ZDC}",nZdcE,minZdcE,maxZdcE);
  fZdcEnergy->SetXTitle("E_{ZDC}");
  fZdcEnergy->SetDirectory(0);
  out->Add(fZdcEnergy);
  
  Int_t    nZemE   = 1000;
  Double_t minZemE = -100;
  Double_t maxZemE = 2900;
  fZemEnergy           = new TH1D("ZEMEnergy","E_{ZEM}",nZemE,minZemE,maxZemE);
  fZemEnergy->SetXTitle("E_{ZEM}");
  fZemEnergy->SetDirectory(0);
  out->Add(fZemEnergy);

  fCorrelationZemZdc   = new TH2D("corrZEMZDC","E_{ZEM} vs E_{ZDC}",
				  nZemE, minZemE, maxZemE, 
				  nZdcE, minZdcE, maxZdcE);
  fCorrelationZemZdc->SetXTitle("E_{ZEM}");
  fCorrelationZemZdc->SetYTitle("E_{ZDC}");
  fCorrelationZemZdc->SetDirectory(0);
  out->Add(fCorrelationZemZdc);
  
  fCorrelationSumDelta =  new TH2D("corrSumDelta","#sumTDC vs #DeltaTDC",
				   nSumT, minSumT, maxSumT,
				   nDeltaT, minDeltaT, maxDeltaT);
  fCorrelationSumDelta->SetXTitle("#sum TDC");
  fCorrelationSumDelta->SetYTitle("#DeltaTDC");
  fCorrelationSumDelta->SetDirectory(0);  
  out->Add(fCorrelationSumDelta);

}

  
//____________________________________________________________________
void
AliDisplacedVertexSelection::Print(Option_t*) const
{
#if 0
  char ind[gROOT->GetDirLevel()+1];
  for (Int_t i = 0; i < gROOT->GetDirLevel(); i++) ind[i] = ' ';
  ind[gROOT->GetDirLevel()] = '\0';
  std::cout << std::boolalpha 
	    << std::noboolalpha << std::endl;
#endif
}

//____________________________________________________________________
Float_t
AliDisplacedVertexSelection::GetZemCorr(Int_t k, Bool_t minusminus) const
{
  if (-kMaxK > k || k > kMaxK) return 0;

  // Corrections for magnetic fields
  const Float_t kPlusPlus[21]   = { 0.6840,0.7879,0.8722,
				    0.9370,0.9837,1.0137,
				    1.0292,1.0327,1.0271,
				    1.0152,1.0000,0.9844,
				    0.9714,0.9634,0.9626,
				    0.9708,0.9891,1.0181,
				    1.0574,1.1060,1.1617};
  const Float_t kMoinsMoins[21] = { 0.7082,0.8012,0.8809,
				    0.9447,0.9916,1.0220,
				    1.0372,1.0395,1.0318,
				    1.0174,1.0000,0.9830,
				    0.9699,0.9635,0.9662,
				    0.9794,1.0034,1.0371,
				    1.0782,1.1224,1.1634};
  
  // MC specific code by Hans, is it used - why are ++ and -- the
  // same?
  const Float_t kPlusPlusMC[21]   = { 0.68400, 0.78790, 0.87220,
				      0.93700, 0.98370, 1.08982,
				      1.03518, 1.01534, 1.01840,
				      1.01493, 1.00000, 0.970083,
				      0.979705,0.960576,0.925394,
				      0.92167, 0.971241,1.08338,
				      1.23517, 1.39308, 1.53943 };
  const Float_t kMoinsMoinsMC[21] = { 0.68400, 0.78790, 0.87220,
				      0.9370,  0.98370, 1.08982,
				      1.03518, 1.01534, 1.0184,
				      1.01493, 1.00000, 0.970083,
				      0.979705,0.960576,0.925394,
				      0.92167, 0.971241,1.08338,
				      1.23517, 1.39308, 1.53943 };
  Int_t i = k+kMaxK;
  if (fMC) 
    return minusminus ? kMoinsMoinsMC[i] : kPlusPlusMC[i];
  return minusminus ? kMoinsMoins[i] : kPlusPlus[i];
}


//____________________________________________________________________
Bool_t
AliDisplacedVertexSelection::CheckOutlier(Int_t ivtx, 
					  const AliESDEvent* esd) const
{
  if (fMC) return false;
  if (ivtx <= 4) return false; // Not outliers

  // --- Find outliers -----------------------------------------------
  AliESDVZERO* esdV0 = esd->GetVZEROData(); 
  Double_t nClusters = esd->GetMultiplicity()->GetNumberOfITSClusters(1);

  // parameter from the fit of VZERO multiplicity vs N SPD cluster
  // distribution
  const Double_t slp[8][21] = {{0.000,0.000,0.000,0.000,0.000,0.000,0.000,//0
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.000,0.000,0.000,0.000,0.000,0.000,0.000},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//1
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.000,0.000,0.000,0.000,0.000,0.000,0.000},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//2
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.000,0.000,0.000,0.000,0.000,0.000,0.000},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//3
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.000,0.000,0.000,0.000,0.000,0.000,0.000},
			       {0.000,0.000,0.000,0.000,0.000,0.155,0.173,//4 
				0.238,0.348,0.384,0.209,0.362,0.433,0.435,
				0.421,0.400,0.403,0.398,0.407,0.394,0.369},
			       {0.000,0.000,0.000,0.000,0.000,0.356,0.402,//5
				0.504,0.650,0.643,0.331,0.543,0.640,0.642,
				0.614,0.591,0.587,0.564,0.572,0.524,0.576},
			       {0.000,0.000,0.000,0.000,0.000,0.386,0.500,//6
				0.692,0.856,0.810,0.390,0.619,0.713,0.708,
				0.670,0.633,0.622,0.593,0.587,0.538,0.587},
			       {0.000,0.000,0.000,0.000,0.000,0.353,0.454,//7
				0.697,0.944,0.941,0.444,0.670,0.746,0.736,
				0.683,0.642,0.619,0.583,0.576,0.535,0.581}};
  const Double_t cst[8][21] = {{0.000,0.000,0.000,0.000,0.000,0.000,0.000,//0
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.0000,0.000,0.000,0.000,0.000,0.000,0.00},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//1
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.0000,0.000,0.000,0.000,0.000,0.000,0.00},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//2
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.0000,0.000,0.000,0.000,0.000,0.000,0.00},
			       {0.000,0.000,0.000,0.000,0.000,0.000,0.000,//3
				0.000,0.000,0.000,0.000,0.000,0.000,0.000,
				0.0000,0.000,0.000,0.000,0.000,0.000,0.00},
			       {0.000,0.000,0.000,0.000,0.000,51.67,54.79,//4
				53.51,23.49,24.60,28.47,27.35,24.47,23.93,
				13.162,12.91,5.305,4.007,-14.34,-9.16,-9.43},
			       {0.000,0.000,0.000,0.000,0.000,136.34,88.84,
				84.85,24.31,27.92,35.64,33.15,32.25,23.07,
				13.746,-6.05,-11.95,1.995,-24.91,-22.55,-37.86},
			       {0.000,0.000,0.000,0.000,0.000,380.42,204.22,
				123.89,34.60,29.26,32.49,30.94,19.67,7.760,
				1.1010,-6.01,-15.66,-3.16,-18.34,-17.33,-22.17},
			       {0.000,0.000,0.000,0.000,0.000,128.00,105.04,
				114.05,15.96,21.67,30.11,30.69,27.66,0.363,-0.511,
				-17.67,-22.40,-11.84,-25.65,-24.29,-24.46}};

  Bool_t outlier = kFALSE;
  // loop over VZERO rings
  for (int iring = 4; iring < 8; iring++) { 
    Double_t multRing = 0;
    for (int iCh=0; iCh<8; iCh++) {
      Int_t idx = iCh+iring*8;
      multRing += esdV0->GetMultiplicity(idx);
    }
    
    // Tigher cut for z=0.  Looser cut for suffician statistics and
    // see saturation effects if present.
    Double_t upFac    = (ivtx == 10 ? .35 : .42);
    Double_t downFac  = (ivtx == 10 ? .20 : ivtx < 8 ? .42 : .26);
    Double_t slpUP    = slp[iring][ivtx]*(1.+upFac); //upper cut
    Double_t slpDOWN  = slp[iring][ivtx]*(1.-downFac); //lower cut
    Double_t constant = cst[iring][ivtx];
    //Cut is applied for rings and vertex of interest
    if ((slpDOWN * nClusters + constant) > multRing || 
	(slpUP   * nClusters + constant) < multRing ) { 
	outlier = kTRUE;
    }
  }

  return outlier;
  //--- End of outlier code ----------------------------------------- 
}
//____________________________________________________________________
Bool_t
AliDisplacedVertexSelection::ProcessMC(const AliMCEvent* mcevent)
{
  if (!fMC) return true;
  if (!mcevent) return false;

  AliMCEvent*        event     = const_cast<AliMCEvent*>(mcevent);
  AliHeader*         header    = event->Header();
  AliGenEventHeader* genHeader = header->GenEventHeader();
  TArrayF            vertex;
  genHeader->PrimaryVertex(vertex);
  Double_t zvtx           = vertex.At(2);
  Double_t intTime        = genHeader->InteractionTime();
  const Double_t kVtx[16] = {-187.5, -150., -112.5, -75., -37.5, 0.,
			     37.5,   75.,  112.5, 150., 187.5, 225.,
			     262.5,  300.,  337.5, 375.};
  Int_t vtxClass          = -1;
  Float_t SigmaVtx        = 5.3;
  Float_t SigmaTime       = 0.18;
  
  for(Int_t iVtx=0; iVtx<16; ++iVtx) {
    // Do not do something for nominal
    if (iVtx == 5) continue;
    if ((zvtx >= kVtx[iVtx] - 3.54 * SigmaVtx && 
	 zvtx <= kVtx[iVtx] + 3.54 * SigmaVtx) && 
	(intTime >= (iVtx * 1.25e-9 - 6.25e-9 - 4 * SigmaTime) && 
	 intTime <= (iVtx * 1.25e-9 - 6.25e-9 + 4 * SigmaTime)))  {
      vtxClass = iVtx;
      // break here to not do more 
      break;
    }
  }
  if(vtxClass > -1) {
    fVertexZ = kVtx[vtxClass];
    return true;
  }
  return false;
}
//____________________________________________________________________
Bool_t
AliDisplacedVertexSelection::Process(const AliESDEvent* esd)
{
  fVertexZ = kInvalidVtxZ; // Default vertex value 
  fCent    = 100;  // Default centrality value 

  // Some constants 
  const Float_t kZDCrefSum        = -66.5;
  const Float_t kZDCrefDelta      =  -2.10;
  const Float_t kZDCsigmaSum      =   3.25;
  const Float_t kZDCsigmaDelta    =   2.25;
  const Float_t kZDCsigmaSumSat   =   2.50;
  const Float_t kZDCsigmaDeltaSat =   1.35;  
  
  // --- Get the ESD object ------------------------------------------
  AliESDZDC* esdZDC = esd->GetESDZDC();
  if (!esdZDC) { 
    AliWarning("No ZDC ESD object!");
    return false; 
  }
  Double_t currentL3   = esd->GetCurrentL3();
  Double_t currentDipo = esd->GetCurrentDip();

  // --- ZDC and ZEM energy signal -----------------------------------
  Double_t zdcEn      = (esdZDC->GetZDCN1Energy()+
			 esdZDC->GetZDCP1Energy()+
			 esdZDC->GetZDCN2Energy()+
			 esdZDC->GetZDCP2Energy());
  Double_t zemEn      = (esdZDC->GetZDCEMEnergy(0)+
			 esdZDC->GetZDCEMEnergy(1))/8.;
 
  // HHD/Maxime inclusion - MC check!
  if (fMC) {
    zdcEn      = (2.9 * esdZDC->GetZDCN1Energy() + 
		  7.2 * esdZDC->GetZDCP1Energy() + 
		  3   * esdZDC->GetZDCN2Energy() +
		  8.7 * esdZDC->GetZDCP2Energy());
    zemEn      = 0.57 * (esdZDC->GetZDCEMEnergy(0) + 
			 esdZDC->GetZDCEMEnergy(1));
  }

  // --- Calculate DeltaT and sumT -----------------------------------
  Double_t deltaTdc = 0;
  Double_t sumTdc   = 0;
   
  for(Int_t i = 0; i < 4; ++i) {
    if(esdZDC->GetZDCTDCData(10,i) != 0) {
      Double_t  tdcCnoCorr = 0.025*(esdZDC->GetZDCTDCData(10,i)-
				    esdZDC->GetZDCTDCData(14,i)); 
      Double_t  tdcC       = esdZDC->GetZDCTDCCorrected(10,i); 
      for(Int_t j = 0; j < 4; ++j) {
	if(esdZDC->GetZDCTDCData(12,j) != 0) {
	  Double_t   tdcAnoCorr = 0.025*(esdZDC->GetZDCTDCData(12,j)-
					 esdZDC->GetZDCTDCData(14,j));
	  Double_t   tdcA       = esdZDC->GetZDCTDCCorrected(12,j);
	  if(esdZDC->TestBit(AliESDZDC::kCorrectedTDCFilled)) {	    
	    deltaTdc = tdcC-tdcA;
	    sumTdc   = tdcC+tdcA;
	  }
	  else {
	    deltaTdc = tdcCnoCorr-tdcAnoCorr;
	    sumTdc   = tdcCnoCorr+tdcAnoCorr;
	  }
	}
      }
    }
  }
  fDeltaTdc->Fill(deltaTdc);
  fSumTdc->Fill(sumTdc);
  fCorrelationSumDelta->Fill(sumTdc, deltaTdc);

  // --- Find the vertex ---------------------------------------------
  Int_t ivtx = -1;
  if(deltaTdc!=0. || sumTdc!=0.) {
    Double_t fillVz = kInvalidVtxZ;
    for (Int_t k = -kMaxK; k <= kMaxK; ++k) {
      Float_t zsat  = 2.55F * k;
      Float_t delta = (k == 0 ? kZDCsigmaDelta : kZDCsigmaDeltaSat);
      Float_t sum   = (k == 0 ? kZDCsigmaSum   : kZDCsigmaSumSat);
      Float_t dT    = deltaTdc - kZDCrefDelta - zsat;
      Float_t sT    = sumTdc  - kZDCrefSum - zsat;
      Float_t check = dT * dT / delta / delta + sT * sT / sum  / sum;
      if (check > 1.0) continue;
      if (k == 0) { 
	fillVz = 0;
	continue;
      }
      // Set the vertex 
      fVertexZ = 37.5 * k;
      fillVz   = fVertexZ;    
      ivtx     = k + 10; // Used for outlier calculations
      // Correct zem energy 
      if(currentDipo>0 && currentL3>0) zemEn /= GetZemCorr(k,false);
      if(currentDipo<0 && currentL3<0) zemEn /= GetZemCorr(k,true);
      // We need to break here, or we could possibly update zemEn again 
      break;
    }
    if (fillVz != kInvalidVtxZ) fHVertexZ->Fill(fillVz);
  }
  fZemEnergy->Fill(zemEn);
  fZdcEnergy->Fill(zdcEn);
  fCorrelationZemZdc->Fill(zemEn, zdcEn);

  // --- Check if this is an outlier event ---------------------------
  if (CheckOutlier(ivtx, esd)) {
    fVertexZ = kInvalidVtxZ;
    return false;
  }

  // --- Calculate the centrality ------------------------------------
  Float_t c1, c2, c3, c4, c5;
  Int_t runnumber = esd->GetRunNumber();
#if 1 // New centrality determination by HHD
  if (runnumber < 137165 ) {
    // ref 137161
    c1 = 16992;
    c2 = 345;
    c3 = 2.23902e+02;
    c4 = 1.56667;
    c5 = 9.49434e-05;
  }
  else if (runnumber <= 137848 && runnumber >= 137230) {
    // ref 137722
    c1  = 15000;
    c2  = 295;
    c3  = 2.23660e+02;
    c4  = 1.56664;
    c5  = 8.99571e-05;
  }
  else if (runnumber <= 138275 && runnumber >= 138125) {
    // ref 137161, yes that's the correct run!
    c1 = 16992;
    c2 = 345;
    c3 = 2.23902e+02;
    c4 = 1.56667;
    c5 = 9.49434e-05;
  }
  else { // if (runnumber <= 139517 && runnumber >= 138358) {
    // Ref run 139172
    c1 = 16992;
    c2 = 345;
    c3 = 2.04789e+02;
    c4 = 1.56629;
    c5 = 1.02768e-04;
  }
#else // Old centrality determination
  if (runnumber < 137722 && runnumber >= 137161 ) {
    c1 = 16992;
    c2 = 345;
    c3 = 2.23902e+02;
    c4 = 1.56667;
    c5 = 9.49434e-05;
  }
  else if (runnumber < 139172 && runnumber >= 137722) {
    c1  = 15000;
    c2  = 295;
    c3  = 2.23660e+02;
    c4 = 1.56664;
    c5 = 8.99571e-05;
  }
  else if (runnumber >= 139172) {
    c1 = 16992;
    c2 = 345;
    c3 = 2.04789e+02;
    c4 = 1.56629;
    c5 = 1.02768e-04;
  }
  else {
    c1 = 16992;
    c2 = 345;
    c3 = 2.04789e+02;
    c4 = 1.56629;
    c5 = 1.02768e-04;
  }
#endif
  if (zemEn > c2) { 
    Float_t slope   = (zdcEn + c1) / (zemEn - c2) + c3;
    Float_t zdcCent = (TMath::ATan(slope) - c4) / c5;
    if (zdcCent >= 0) fCent = zdcCent;
  }
  
  fHCent->Fill(fCent);

  return true;
}
//____________________________________________________________________
//
// EOF
//
