#include "AliAnaChargedJetResponseMaker.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "Riostream.h"
#include "TH1.h"
#include "TRandom.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "THnSparse.h"
#include "TArrayD.h"

#define round(x) ((x)>=0?(long)((x)+0.5):(long)((x)-0.5))

using std::cout;
using std::endl;

ClassImp(AliAnaChargedJetResponseMaker)

AliAnaChargedJetResponseMaker::AliAnaChargedJetResponseMaker(): 
  fDebug(kFALSE),
  fResolutionType(kParam),
  fDeltaPt(0x0),
  fhDeltaPt(0x0),
  fh1MeasuredTruncated(0x0),
  fh2DetectorResponse(0x0),
  fh2DetectorResponseRebin(0x0),
  fhEfficiencyDet(0x0),
  fh2ResponseMatrixCombinedFineFull(0x0),
  fh2ResponseMatrixCombinedFull(0x0),
  fh2ResponseMatrixCombined(0x0),
  fhEfficiencyCombined(0x0),
  fDimensions(1),
  fDimRec(0),
  fDimGen(1),
  fPtMin(-999),
  fPtMax(-999),
  fNbins(0),
  fBinArrayPtRec(0x0),
  fPtMeasured(0x0),
  fEffFlat(0),
  fEfficiency(0x0),
  fEfficiencyFine(0x0),
  fResponseMatrix(0x0),
  fResponseMatrixFine(0x0),
  fPtMinUnfolded(0.),
  fPtMaxUnfolded(0.),
  fPtMaxUnfoldedHigh(-1.),
  fBinWidthFactorUnfolded(2),
  fSkipBinsUnfolded(0),
  fExtraBinsUnfolded(5),
  fbVariableBinning(kFALSE),
  fPtMaxUnfVarBinning(0),
  f1MergeFunction(0x0),
  fFineFrac(10),
  fbCalcErrors(kFALSE)
{;}


//--------------------------------------------------------------------------------------------------------------------------------------------------
AliAnaChargedJetResponseMaker::AliAnaChargedJetResponseMaker(const AliAnaChargedJetResponseMaker& obj):
  fDebug(obj.fDebug),
  fResolutionType(obj.fResolutionType),
  fDeltaPt(obj.fDeltaPt),
  fhDeltaPt(obj.fhDeltaPt),
  fh1MeasuredTruncated(obj.fh1MeasuredTruncated),
  fh2DetectorResponse(obj.fh2DetectorResponse),
  fh2DetectorResponseRebin(obj.fh2DetectorResponseRebin),
  fhEfficiencyDet(obj.fhEfficiencyDet),
  fh2ResponseMatrixCombinedFineFull(obj.fh2ResponseMatrixCombinedFineFull),
  fh2ResponseMatrixCombinedFull(obj.fh2ResponseMatrixCombinedFull),
  fh2ResponseMatrixCombined(obj.fh2ResponseMatrixCombined),
  fhEfficiencyCombined(obj.fhEfficiencyCombined),
  fDimensions(obj.fDimensions),
  fDimRec(obj.fDimRec),
  fDimGen(obj.fDimGen),
  fPtMin(obj.fPtMin),
  fPtMax(obj.fPtMax),
  fNbins(obj.fNbins),
  fBinArrayPtRec(obj.fBinArrayPtRec),
  fPtMeasured(obj.fPtMeasured),
  fEffFlat(obj.fEffFlat),
  fEfficiency(obj.fEfficiency),
  fEfficiencyFine(obj.fEfficiencyFine),
  fResponseMatrix(obj.fResponseMatrix),
  fResponseMatrixFine(obj.fResponseMatrixFine),
  fPtMinUnfolded(obj.fPtMinUnfolded),
  fPtMaxUnfolded(obj.fPtMaxUnfolded),
  fPtMaxUnfoldedHigh(obj.fPtMaxUnfoldedHigh),
  fBinWidthFactorUnfolded(obj.fBinWidthFactorUnfolded),
  fSkipBinsUnfolded(obj.fSkipBinsUnfolded),
  fExtraBinsUnfolded(obj.fExtraBinsUnfolded),
  fbVariableBinning(obj.fbVariableBinning),
  fPtMaxUnfVarBinning(obj.fPtMaxUnfVarBinning),
  f1MergeFunction(obj.f1MergeFunction),
  fFineFrac(obj.fFineFrac),
  fbCalcErrors(obj.fbCalcErrors)
{;}

//--------------------------------------------------------------------------------------------------------------------------------------------------
AliAnaChargedJetResponseMaker& AliAnaChargedJetResponseMaker::operator=(const AliAnaChargedJetResponseMaker& other)
{
// Assignment
  if(&other == this) return *this;
  AliAnaChargedJetResponseMaker::operator=(other);
  fDebug                    = other.fDebug;
  fResolutionType           = other.fResolutionType;
  fDeltaPt                  = other.fDeltaPt;
  fhDeltaPt                 = other.fhDeltaPt;
  fh1MeasuredTruncated      = other.fh1MeasuredTruncated;
  fh2DetectorResponse       = other.fh2DetectorResponse;
  fh2DetectorResponseRebin  = other.fh2DetectorResponseRebin;
  fhEfficiencyDet           = other.fhEfficiencyDet;
  fh2ResponseMatrixCombinedFineFull = other.fh2ResponseMatrixCombinedFineFull;
  fh2ResponseMatrixCombinedFull = other.fh2ResponseMatrixCombinedFull;
  fh2ResponseMatrixCombined = other.fh2ResponseMatrixCombined;
  fhEfficiencyCombined      = other.fhEfficiencyCombined;
  fDimensions               = other.fDimensions;
  fDimRec                   = other.fDimRec;
  fDimGen                   = other.fDimGen;
  fPtMin                    = other.fPtMin;
  fPtMax                    = other.fPtMax;
  fNbins                    = other.fNbins;
  fBinArrayPtRec            = other.fBinArrayPtRec;
  fPtMeasured               = other.fPtMeasured;
  fEffFlat                  = other.fEffFlat;
  fEfficiency               = other.fEfficiency;
  fEfficiencyFine           = other.fEfficiencyFine;
  fResponseMatrix           = other.fResponseMatrix;
  fResponseMatrixFine       = other.fResponseMatrixFine;
  fPtMinUnfolded            = other.fPtMinUnfolded;
  fPtMaxUnfolded            = other.fPtMaxUnfolded;
  fPtMaxUnfoldedHigh        = other.fPtMaxUnfoldedHigh;
  fBinWidthFactorUnfolded   = other.fBinWidthFactorUnfolded;
  fSkipBinsUnfolded         = other.fSkipBinsUnfolded;
  fExtraBinsUnfolded        = other.fExtraBinsUnfolded;
  fbVariableBinning         = other.fbVariableBinning;
  fPtMaxUnfVarBinning       = other.fPtMaxUnfVarBinning;
  f1MergeFunction           = other.f1MergeFunction;
  fFineFrac                 = other.fFineFrac;
  fbCalcErrors              = other.fbCalcErrors;

  return *this;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::SetMeasuredSpectrum(TH1D *hPtMeasured)
{
  //
  // Set measured spectrum in THnSparse format
  // This defines the minimum and maximum pT on the reconstructed axis of the response matrix
  //
  if(fDebug) printf(">>AliAnaChargedJetResponseMaker::SetMeasuredSpectrum \n");

  fNbins = hPtMeasured->GetXaxis()->GetNbins();
  fPtMin = hPtMeasured->GetXaxis()->GetXmin();
  fPtMax = hPtMeasured->GetXaxis()->GetXmax();

  if(fDebug) printf("fNbins: %d  fPtMin: %f  fPtMax: %f \n",fNbins,fPtMin,fPtMax);
  
  if(fBinArrayPtRec) delete fBinArrayPtRec;
  fBinArrayPtRec = new Double_t[fNbins+1];
  for(int j = 0; j<fNbins; j++) {
    fBinArrayPtRec[j] = hPtMeasured->GetXaxis()->GetBinLowEdge(j+1);
  }
  fBinArrayPtRec[fNbins] = hPtMeasured->GetXaxis()->GetBinUpEdge(fNbins);
  

  Int_t nbins[fDimensions];
  Double_t xmin[fDimensions]; 
  Double_t xmax[fDimensions];
  for(int dim = 0; dim<fDimensions; dim++) {
    nbins[dim] = fNbins;
    xmin[dim]  = fPtMin;
    xmax[dim]  = fPtMax;
  }

  if(fPtMeasured) delete fPtMeasured;
  fPtMeasured = new THnSparseD("fPtMeasured","Measured pT spectrum; p_{T}^{rec}",fDimensions,nbins,xmin,xmax);
  fPtMeasured->Sumw2();
  fPtMeasured->GetAxis(0)->SetTitle("p_{T}^{rec}");
  fPtMeasured->SetBinEdges(0,fBinArrayPtRec);

  //Fill
  if(fDebug) printf("fill measured THnSparse\n");
  if(fNbins!=hPtMeasured->GetNbinsX()) 
    printf("WARNING: nbins not correct \t %d vs %d !!!\n",fNbins,hPtMeasured->GetNbinsX());

  int bin[1] = {0};
  bin[0] = 0;
  for(int i = hPtMeasured->FindBin(fPtMin); i<hPtMeasured->FindBin(fPtMax); i++) {
    fPtMeasured->SetBinContent(bin,(Double_t)hPtMeasured->GetBinContent(i));
    fPtMeasured->SetBinError(bin,(Double_t)hPtMeasured->GetBinError(i));
    bin[0]++;
  }
  
  if(fDebug) printf("fPtMeasured->GetNbins(): %d \n",(int)fPtMeasured->GetNbins());

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::SetFlatEfficiency(Double_t eff) {

  //
  // Set flat efficiency to value of eff
  //

  fEffFlat = eff;
  return;

  /*
  Int_t nbins[fDimensions];
  Double_t xmin[fDimensions]; 
  Double_t xmax[fDimensions];
  for(int dim = 0; dim<fDimensions; dim++) {
    nbins[dim] = fNbins;
    xmin[dim] = fPtMin;
    xmax[dim] = fPtMax;
  }

  if(fEfficiency) delete fEfficiency;
  fEfficiency = new THnSparseD("fEfficiency","Efficiency - no resolution effects; p_{T}^{gen}",fDimensions,nbins,xmin,xmax);
  fEfficiency->Sumw2();
  fEfficiency->GetAxis(0)->SetTitle("p_{T}^{gen}");
  fEfficiency->SetBinEdges(0,fBinArrayPtRec);

  for(int i=0; i<fNbins; i++) {
    int bin[1];
    bin[0] = i;
    fEfficiency->SetBinContent(bin,fEffFlat);
    fEfficiency->SetBinError(bin,0);
  }
  */
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::SetEfficiency(TGraphErrors *grEff)
{
  //
  // Fill fEfficiency (THnSparse) with values from grEff
  //

  Int_t nbins[fDimensions];
  Double_t xmin[fDimensions]; 
  Double_t xmax[fDimensions];
  for(int dim = 0; dim<fDimensions; dim++) {
    nbins[dim] = fNbins;
    xmin[dim] = fPtMin;
    xmax[dim] = fPtMax;
  }

  if(fEfficiency) delete fEfficiency;
  fEfficiency = new THnSparseD("fEfficiency","Efficiency - no resolution effects; p_{T}^{gen}",fDimensions,nbins,xmin,xmax);
  fEfficiency->Sumw2();
  fEfficiency->GetAxis(0)->SetTitle("p_{T}^{gen}");
  fEfficiency->SetBinEdges(0,fBinArrayPtRec);

  double pT[1]; 
  double yield = 0.;
  double error = 0.;
  double dummy = 0.;
  int nbinsgr = grEff->GetN();
  
  for(int i=0; i<nbinsgr; i++) {
    grEff->GetPoint(i,dummy,yield);
    pT[0] = dummy;
    error = grEff->GetErrorY(i);
    
    fEfficiency->Fill(pT,yield);
    fEfficiency->SetBinError(i,error);

  }
  
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::MakeResponseMatrixJetsFineMerged(Int_t skipBins, Int_t binWidthFactor, Int_t extraBins, Bool_t bVariableBinning, Double_t ptmax)
{
  //
  // Make jet response matrix
  //

  if(fDebug) printf(">>AliAnaChargedJetResponseMaker::MakeResponseMatrixJetsFineMerged\n");

  if(!fPtMeasured) {
    if(fDebug) printf("fPtMeasured does not exist. Aborting!!\n");
    return;
  }
  if(fResponseMatrix)     { delete fResponseMatrix; }
  if(fResponseMatrixFine) { delete fResponseMatrixFine; }

  SetSkipBinsUnfolded(skipBins);
  SetBinWidthFactorUnfolded(binWidthFactor);
  SetExtraBinsUnfolded(extraBins);
  SetVariableBinning(bVariableBinning,ptmax);

  InitializeResponseMatrix();
  InitializeEfficiency();

  InitializeResponseMatrixFine();
  InitializeEfficiencyFine();

  FillResponseMatrixFineAndMerge();

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::InitializeResponseMatrix() {
  //
  //Define bin edges of RM to be used for unfolding
  //

  Int_t nbins[fDimensions*2];
  nbins[fDimRec] = fNbins;
  nbins[fDimGen] = fNbins;

  double binWidthMeas = (double)((fPtMax-fPtMin)/fNbins);
  double binWidthUnf  = (double)fBinWidthFactorUnfolded*binWidthMeas;
  double binWidthUnfLowPt = -1.;
  if(fbVariableBinning) 
    binWidthUnfLowPt = binWidthUnf*0.5;

  fPtMaxUnfolded = fPtMax+(double)(fExtraBinsUnfolded)*binWidthUnf;
  nbins[fDimGen]+=fExtraBinsUnfolded;


  printf("fPtMinMeas: %f  fPtMaxMeas: %f\n",fPtMin,fPtMax);
  printf("binWidthMeas: %f  binWidthUnf: %f   fBinWidthFactorUnfolded: %d\n",binWidthMeas,binWidthUnf,fBinWidthFactorUnfolded);
  printf("binWidthUnfLowPt: %f\n",binWidthUnfLowPt);

  printf("fPtMinUnfolded: %f  fPtMaxUnfolded: %f\n",fPtMinUnfolded,fPtMaxUnfolded);


  if(fbVariableBinning) {
    // cout << "fPtMaxUnfVarBinning: " << fPtMaxUnfVarBinning << " \tfPtMinUnfolded: " << fPtMinUnfolded << "  binWidthUnfLowPt: " << binWidthUnfLowPt << endl;
    Int_t tmp = (int)((fPtMaxUnfVarBinning-fPtMinUnfolded)/binWidthUnfLowPt);
    tmp = tmp - fSkipBinsUnfolded;
    fPtMinUnfolded = fPtMaxUnfVarBinning-(double)(tmp)*binWidthUnfLowPt;  
    //cout << "tmp = " << tmp << "  fSkipBinsUnfolded = " << fSkipBinsUnfolded << " fPtMinUnfolded = " << fPtMinUnfolded << endl;
    //Redefine also number of bins on generated axis in case of variable binning
    nbins[fDimGen] = (int)((fPtMaxUnfVarBinning-fPtMinUnfolded)/binWidthUnfLowPt)+(int)((fPtMaxUnfolded-fPtMaxUnfVarBinning)/binWidthUnf);
  }
  else {
    int tmp = round((fPtMaxUnfolded-fPtMinUnfolded)/binWidthUnf); //nbins which fit between 0 and fPtMaxUnfolded
    tmp = tmp - fSkipBinsUnfolded;
    fPtMinUnfolded = fPtMaxUnfolded-(double)(tmp)*binWidthUnf; //set ptmin unfolded
    fPtMaxUnfolded = fPtMinUnfolded+(double)(tmp)*binWidthUnf; //set ptmax unfolded
    nbins[fDimGen] = (int)((fPtMaxUnfolded-fPtMinUnfolded)/binWidthUnf); //adjust nbins to bin width
  }

  printf("   nbins[fDimGen] = %d   nbins[fDimRec] = %d\n",nbins[fDimGen],nbins[fDimRec]);

  Double_t binWidth[2];
  binWidth[fDimRec] = (double)((fPtMax-fPtMin)/nbins[fDimRec]);
  binWidth[fDimGen] = (double)((fPtMaxUnfolded-fPtMinUnfolded)/nbins[fDimGen]);

  Double_t xmin[fDimensions*2]; 
  Double_t xmax[fDimensions*2];
  xmin[fDimRec] = fPtMin;
  xmax[fDimRec] = fPtMax;
  xmin[fDimGen] = fPtMinUnfolded;
  xmax[fDimGen] = fPtMaxUnfolded;

  printf("xmin[fDimRec]: %f  xmin[fDimGen]: %f \n",xmin[fDimRec],xmin[fDimGen]);
  printf("xmax[fDimRec]: %f  xmax[fDimGen]: %f \n",xmax[fDimRec],xmax[fDimGen]);
  printf("nbins[fDimRec]: %d  nbins[fDimGen]: %d \n",nbins[fDimRec],nbins[fDimGen]);
  if(!fbVariableBinning) printf("binWidth[fDimRec]: %f  binWidth[fDimGen]: %f \n",binWidth[fDimRec],binWidth[fDimGen]);

  Double_t binArrayPt0[nbins[fDimRec]+1];
  Double_t binArrayPt1[nbins[fDimGen]+1];
  Double_t xmaxGen = TMath::Max(xmax[fDimGen],fPtMaxUnfoldedHigh);

  //Define bin limits reconstructed/measured axis
  for(int i=0; i<nbins[fDimRec]; i++) {
    binArrayPt0[i] = xmin[fDimRec]+(double)i*binWidth[fDimRec];
  }
  binArrayPt0[nbins[fDimRec]]= xmax[fDimRec];

  //Define bin limits generated/unfolded axis
  binArrayPt1[0] = fPtMinUnfolded;
  for(int i=1; i<nbins[fDimGen]; i++) {
    if(fbVariableBinning) {
      double test = xmin[fDimGen]+(double)i*binWidthUnfLowPt;
      if(test<=fPtMaxUnfVarBinning) binArrayPt1[i] = test;
      else binArrayPt1[i] = binArrayPt1[i-1]+binWidthUnf;
    }
    else
      binArrayPt1[i] = xmin[fDimGen]+(double)i*binWidth[fDimGen];
    //printf("RM. i = %d \t binArrayPt[i] = %f \n",i,binArrayPt1[i]);
  }
  binArrayPt1[nbins[fDimGen]]= xmaxGen;


  // Response matrix : dimensions must be 2N = 2 x (number of variables)
  // dimensions 0 ->  N-1 must be filled with reconstructed values
  // dimensions N -> 2N-1 must be filled with generated values
  if(fResponseMatrix) delete fResponseMatrix;
  fResponseMatrix = new THnSparseD("fResponseMatrix","Response Matrix pTMC vs pTrec",fDimensions*2,nbins,xmin,xmax);
  fResponseMatrix->Sumw2();
  fResponseMatrix->GetAxis(fDimRec)->SetTitle("p_{T}^{rec}");
  fResponseMatrix->GetAxis(fDimGen)->SetTitle("p_{T}^{gen}");

  fResponseMatrix->SetBinEdges(fDimRec,binArrayPt0);
  fResponseMatrix->SetBinEdges(fDimGen,binArrayPt1);

  Int_t bin[2] = {1,1};
  for(int i=1; i<fResponseMatrix->GetAxis(0)->GetNbins(); i++) {
    bin[0]=i;
    for(int j=1; j<fResponseMatrix->GetAxis(1)->GetNbins(); j++) {
    bin[1]=j;
      fResponseMatrix->SetBinContent(bin,0.);
    }
  }



}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::InitializeEfficiency() {
  //
  // Define dimensions of efficiency THnSparse
  //

  if(!fResponseMatrix) {
    printf("AliAnaChargedJetResponseMaker::InitializeEfficiency()\n");
    printf("Can not define dimensions efficiency without fResponseMatrix dimensions. Aborting \n");
    return;
  }

  TAxis *genAxis = fResponseMatrix->GetAxis(fDimGen);

  Int_t nbinsEff[fDimensions];
  Double_t xminEff[fDimensions]; 
  Double_t xmaxEff[fDimensions];

  for(int dim = 0; dim<fDimensions; dim++) {
    nbinsEff[dim] = genAxis->GetNbins();
    xminEff[dim]  = genAxis->GetXmin();
    xmaxEff[dim]  = genAxis->GetXmax();
  }

  if(fEfficiency) delete fEfficiency;
  fEfficiency = new THnSparseD("fEfficiency","Efficiency - no resolution effects; p_{T}^{gen}",fDimensions,nbinsEff,xminEff,xmaxEff);
  fEfficiency->Sumw2();
  fEfficiency->GetAxis(0)->SetTitle("p_{T}^{gen}");

  const Double_t *binArrayPt = genAxis->GetXbins()->GetArray();
  fEfficiency->SetBinEdges(0,binArrayPt);

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::InitializeResponseMatrixFine() {
  //
  //Initialize fine response matrix
  //

  Int_t nbinsFine[fDimensions*2];
  Double_t xminFine[fDimensions*2];
  Double_t xmaxFine[fDimensions*2];
  // Double_t pTarrayFine[fDimensions*2];

  nbinsFine[fDimRec] = fResponseMatrix->GetAxis(fDimRec)->GetNbins()*fFineFrac; 
  nbinsFine[fDimGen] = fResponseMatrix->GetAxis(fDimRec)->GetNbins()*fFineFrac; 
  xminFine[fDimRec] = TMath::Min(fPtMin,0.);
  xminFine[fDimGen] = TMath::Min(fPtMin,0.);
  xmaxFine[fDimRec] = fResponseMatrix->GetAxis(fDimGen)->GetXmax();//+40.;
  xmaxFine[fDimGen] = fResponseMatrix->GetAxis(fDimGen)->GetXmax();//+40.;
  // pTarrayFine[fDimRec] = 0.;
  // pTarrayFine[fDimGen] = 0.;

  Double_t binWidth[2];
  binWidth[fDimRec] = fResponseMatrix->GetAxis(fDimRec)->GetBinWidth(1);

  Double_t binWidthFine[2];
  binWidthFine[fDimRec] = binWidth[fDimRec]/((double)fFineFrac);
  binWidthFine[fDimGen] = binWidthFine[fDimRec]*(double)fBinWidthFactorUnfolded;
  //  binWidthFine[fDimRec] = binWidthFine[fDimGen];

  nbinsFine[fDimRec] = (int)((xmaxFine[fDimRec]-xminFine[fDimRec])/binWidthFine[fDimRec]); //adjust nbins to bin width
  nbinsFine[fDimGen] = (int)((xmaxFine[fDimGen]-xminFine[fDimGen])/binWidthFine[fDimGen]); //adjust nbins to bin width

  printf("xminFine[fDimRec]: %f  xminFine[fDimGen]: %f \n",xminFine[fDimRec],xminFine[fDimGen]);
  printf("xmaxFine[fDimRec]: %f  xmaxFine[fDimGen]: %f \n",xmaxFine[fDimRec],xmaxFine[fDimGen]);
  printf("nbinsFine[fDimRec]: %d  nbinsFine[fDimGen]: %d \n",nbinsFine[fDimRec],nbinsFine[fDimGen]);
  printf("binWidthFine[fDimRec]: %f  binWidthFine[fDimGen]: %f \n",binWidthFine[fDimRec],binWidthFine[fDimGen]);


  Double_t binArrayPt0Fine[nbinsFine[fDimRec]+1];
  Double_t binArrayPt1Fine[nbinsFine[fDimGen]+1];

  for(int i=0; i<nbinsFine[fDimRec]; i++) {
    binArrayPt0Fine[i] = xminFine[fDimRec]+(double)i*binWidthFine[fDimRec];
    //    printf("RM. i = %d \t binArrayPtFine[i] = %f \n",i,binArrayPt0Fine[i]);
  }
  binArrayPt0Fine[nbinsFine[fDimRec]]= xmaxFine[fDimRec];

  for(int i=0; i<nbinsFine[fDimGen]; i++) {
    binArrayPt1Fine[i] = xminFine[fDimGen]+(double)i*binWidthFine[fDimGen];
    //    printf("RM. i = %d \t binArrayPtFine[i] = %f \n",i,binArrayPt1Fine[i]);
  }
  binArrayPt1Fine[nbinsFine[fDimGen]]= xmaxFine[fDimGen];

  // Response matrix : dimensions must be 2N = 2 x (number of variables)
  // dimensions 0 ->  N-1 must be filled with reconstructed values
  // dimensions N -> 2N-1 must be filled with generated values
  if(fResponseMatrixFine) delete fResponseMatrixFine;
  fResponseMatrixFine = new THnSparseD("fResponseMatrixFine","Response Matrix pTMC vs pTrec",fDimensions*2,nbinsFine,xminFine,xmaxFine);
  fResponseMatrixFine->Sumw2();
  fResponseMatrixFine->GetAxis(fDimRec)->SetTitle("p_{T}^{rec}");
  fResponseMatrixFine->GetAxis(fDimGen)->SetTitle("p_{T}^{gen}");

  fResponseMatrixFine->SetBinEdges(fDimRec,binArrayPt0Fine);
  fResponseMatrixFine->SetBinEdges(fDimGen,binArrayPt1Fine);

  Int_t bin[2] = {1,1};
  for(int i=1; i<fResponseMatrixFine->GetAxis(0)->GetNbins(); i++) {
    bin[0]=i;
    for(int j=1; j<fResponseMatrixFine->GetAxis(1)->GetNbins(); j++) {
    bin[1]=j;
      fResponseMatrixFine->SetBinContent(bin,0.);
    }
  }

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::InitializeEfficiencyFine() {
  //
  // Define dimensions of efficiency THnSparse
  //

  if(!fResponseMatrixFine) {
    printf("AliAnaChargedJetResponseMaker::InitializeEfficiencyFine()\n");
    printf("Can not define dimensions efficiency without fResponseMatrixFine dimensions. Aborting \n");
    return;
  }

  TAxis *genAxis = fResponseMatrixFine->GetAxis(fDimGen);

  Int_t nbinsEff[fDimensions];
  Double_t xminEff[fDimensions]; 
  Double_t xmaxEff[fDimensions];

  for(int dim = 0; dim<fDimensions; dim++) {
    nbinsEff[dim] = genAxis->GetNbins();
    xminEff[dim]  = genAxis->GetXmin();
    xmaxEff[dim]  = genAxis->GetXmax();
  }

  if(fEfficiencyFine) delete fEfficiencyFine;
  fEfficiencyFine = new THnSparseD("fEfficiencyFine","EfficiencyFine - no resolution effects; p_{T}^{gen}",fDimensions,nbinsEff,xminEff,xmaxEff);
  fEfficiencyFine->Sumw2();
  fEfficiencyFine->GetAxis(0)->SetTitle("p_{T}^{gen}");

  const Double_t *binArrayPt = genAxis->GetXbins()->GetArray();
  fEfficiencyFine->SetBinEdges(0,binArrayPt);

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::FillResponseMatrixFineAndMerge() {
  //
  // Fill fine response matrix
  //

  if(!fResponseMatrix) {
    printf("Dimensions of fResponseMatrix have to be defined first. Aborting!");
    return;
  }

  if(!fResponseMatrixFine) {
    printf("Dimensions of fResponseMatrixFine have to be defined first. Aborting!");
    return;
  }

  TAxis *genAxis = fResponseMatrixFine->GetAxis(fDimGen);
  TAxis *recAxis = fResponseMatrixFine->GetAxis(fDimRec);

  Int_t nbinsFine[fDimensions*2];
  //Double_t xminFine[fDimensions*2]; 
  //Double_t xmaxFine[fDimensions*2];
  Double_t pTarrayFine[fDimensions*2];

  nbinsFine[fDimGen] = genAxis->GetNbins();
  nbinsFine[fDimRec] = recAxis->GetNbins();
  //  xminFine[fDimGen]  = genAxis->GetXmin();
  //  xminFine[fDimRec]  = recAxis->GetXmin();
  //  xmaxFine[fDimGen]  = genAxis->GetXmax();
  //  xmaxFine[fDimRec]  = recAxis->GetXmax();
  pTarrayFine[fDimGen] = 0.;
  pTarrayFine[fDimRec] = 0.;

  double sumyield = 0.;
  double sumyielderror2 = 0.;

  int ipt[2]    = {0,0};
  int iptMerged[2]    = {0,0};
  int ierror[2] = {0,0};

  Int_t check = 0;
  Double_t pTgen, pTrec;
  Double_t yield = 0.;
  Double_t error = 0.;

  const int nng = fResponseMatrix->GetAxis(fDimGen)->GetNbins();
  const int nnr = fResponseMatrix->GetAxis(fDimRec)->GetNbins();
  Double_t errorArray[nng][nnr];
  for(int iig =0; iig<nng; iig++) {
    for(int iir =0; iir<nnr; iir++) {
      errorArray[iig][iir] = 0.;
    }
  }

  for(int iy=1; iy<=nbinsFine[fDimGen]; iy++) { //gen
    pTgen = fResponseMatrixFine->GetAxis(fDimGen)->GetBinCenter(iy);
    pTarrayFine[fDimGen] = pTgen;
    ierror[fDimGen]=iy;
    sumyield = 0.;
    check = 0;

    for(int ix=1; ix<=nbinsFine[fDimRec]; ix++) { //rec
      pTrec = fResponseMatrixFine->GetAxis(fDimRec)->GetBinCenter(ix);
      Double_t width = fResponseMatrixFine->GetAxis(fDimRec)->GetBinWidth(ix);
      if(fResolutionType==kParam) {
	yield = fDeltaPt->Eval(pTrec-pTgen);
	error = 0.;
      }
      else if(fResolutionType==kResiduals) {
	yield = fhDeltaPt->Interpolate(pTrec-pTgen);
	error = 0.;
      }
      else if(fResolutionType==kResidualsErr) {
	Double_t deltaPt = pTrec-pTgen;
	int bin = fhDeltaPt->FindBin(deltaPt);
	yield = fhDeltaPt->GetBinContent(bin);
	error = fhDeltaPt->GetBinError(bin);
      }
      yield=yield*width;
      error=error*width;
      //avoid trouble with empty bins in the high pT tail of deltaPt distribution
      if(check==0 && yield>0. && pTrec>pTgen) check=1;
      if(check==1 && yield==0.) ix=nbinsFine[fDimRec];

      sumyield+=yield;
      sumyielderror2 += error*error;

      pTarrayFine[fDimRec] = pTrec;
      ierror[fDimRec]  = ix;
      fResponseMatrixFine->Fill(pTarrayFine,yield);
      if(fbCalcErrors) fResponseMatrixFine->SetBinError(ierror,error);

    }// ix (dimRec) loop

    //Normalize to total number of counts =1
    if(sumyield>1) {
      ipt[fDimGen]=iy;
      for(int ix=1; ix<=nbinsFine[fDimRec]; ix++) {
	ipt[fDimRec]=ix;
	fResponseMatrixFine->SetBinContent(ipt,fResponseMatrixFine->GetBinContent(ipt)/sumyield);
	if(fResolutionType==kResidualsErr && fbCalcErrors) {
	  Double_t A = 1./sumyield*fResponseMatrix->GetBinError(ipt);
	  Double_t B = -1.*fResponseMatrix->GetBinContent(ipt)/sumyield/sumyield*TMath::Sqrt(sumyielderror2);
	  Double_t tmp2 = A*A + B*B;
	  fResponseMatrix->SetBinError(ipt,TMath::Sqrt(tmp2));
	}

      }
    }

    int bin[1];
    bin[0] = iy;
    fEfficiencyFine->SetBinContent(bin,sumyield);
    if(fbCalcErrors) fEfficiencyFine->SetBinError(bin,TMath::Sqrt(sumyielderror2));

    //fill merged response matrix

    //find bin in fine RM corresponding to mimimum/maximum bin of merged RM on rec axis
    int ixMin = fResponseMatrixFine->GetAxis(fDimRec)->FindBin(fResponseMatrix->GetAxis(fDimRec)->GetXmin()); 
    int ixMax = fResponseMatrixFine->GetAxis(fDimRec)->FindBin(fResponseMatrix->GetAxis(fDimRec)->GetXmax());

    if(fResponseMatrixFine->GetAxis(fDimGen)->GetBinLowEdge(iy) >= fResponseMatrix->GetAxis(fDimGen)->GetXmin()) { 
      ipt[fDimGen]=iy;
      iptMerged[fDimGen]=fResponseMatrix->GetAxis(fDimGen)->FindBin(pTgen);

      Double_t weight = 1.;
      if(f1MergeFunction) {
	Double_t loEdge = fResponseMatrix->GetAxis(fDimGen)->GetBinLowEdge(iptMerged[fDimGen]);
	Double_t upEdge = fResponseMatrix->GetAxis(fDimGen)->GetBinUpEdge(iptMerged[fDimGen]);
	Float_t powInteg = f1MergeFunction->Integral(loEdge,upEdge);
	//printf("loEdge = %f  upEdge = %f  powInteg = %f\n",loEdge,upEdge,powInteg);
	if(powInteg>0.)
	  weight = f1MergeFunction->Integral(fResponseMatrixFine->GetAxis(fDimGen)->GetBinLowEdge(iy),fResponseMatrixFine->GetAxis(fDimGen)->GetBinUpEdge(iy))/powInteg;
	//	printf("weight: %f \n", weight );
      } else {
	weight = 1./((double)fFineFrac);
	if(fbVariableBinning && pTgen<fPtMaxUnfVarBinning) weight=1./(0.5*(double)fFineFrac);
      }

      for(int ix=ixMin; ix<=ixMax; ix++) {                    //loop reconstructed axis
	pTrec = fResponseMatrixFine->GetAxis(fDimRec)->GetBinCenter(ix);
	ipt[fDimRec]=ix;
	iptMerged[fDimRec]=fResponseMatrix->GetAxis(fDimRec)->FindBin(pTrec);

	fResponseMatrix->AddBinContent(iptMerged,fResponseMatrixFine->GetBinContent(ipt)*weight);
	if(fbCalcErrors) errorArray[iptMerged[fDimGen]-1][iptMerged[fDimRec]-1] += fResponseMatrixFine->GetBinError(ipt)*fResponseMatrixFine->GetBinError(ipt)*weight*weight;
      }

   }//loop gen axis fine response matrix

  } // iy (dimGen) loop

  //Fill Efficiency corresponding to merged response matrix
  for(int iy=1; iy<=fResponseMatrix->GetAxis(fDimGen)->GetNbins(); iy++) { //gen
    sumyield = 0.;
    ipt[fDimGen] = iy;

    for(int ix=1; ix<=fResponseMatrix->GetAxis(fDimRec)->GetNbins(); ix++) { //rec
      ipt[fDimRec] = ix;
      sumyield += fResponseMatrix->GetBinContent(ipt);
      
      if(fbCalcErrors) fResponseMatrix->SetBinError(ipt,TMath::Sqrt(errorArray[iy][ix]));
    }
    printf("RM: pTgen: %f \t sumyield: %f \n",fResponseMatrix->GetAxis(fDimGen)->GetBinCenter(iy),sumyield);
    int bin[1];
    bin[0] = iy;
    fEfficiency->SetBinContent(bin,sumyield);
    if(fbCalcErrors) fEfficiency->SetBinError(bin,0);
  }
  
  if(fDebug) printf("fResponseMatrixFine->GetNbins(): %d \n",(int)fResponseMatrixFine->GetNbins());
  if(fDebug) printf("fResponseMatrix->GetNbins(): %d \n",(int)fResponseMatrix->GetNbins());

  if(fDebug) printf("Done constructing response matrix\n");

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::MakeResponseMatrixRebin(TH2 *hRMFine, TH2 *hRM, Bool_t useFunctionWeight) {

  //
  // Rebin matrix hRMFine to dimensions of hRM
  // function returns matrix in TH2D format (hRM2) with dimensions from hRM
  //

  TH2 *hRM2 = (TH2*)hRM->Clone("hRM2");
  hRM2->Reset();

  if(useFunctionWeight)  cout << "Use function to do weighting" << endl;

  //First normalize columns of input
  const Int_t nbinsNorm = hRM2->GetNbinsX();
  if(fDebug) cout << "nbinsNorm: " << nbinsNorm << endl;

  TArrayF *normVector = new TArrayF(nbinsNorm);

  for(int ix=1; ix<=hRM2->GetNbinsX(); ix++) {
    Double_t xLow = hRM2->GetXaxis()->GetBinLowEdge(ix);
    Double_t xUp = hRM2->GetXaxis()->GetBinUpEdge(ix);
    //cout << "xLow: " << xLow << " xUp: " << xUp << "\t center: " << hRM2->GetXaxis()->GetBinCenter(ix) << endl;
    Int_t binxLowFine = hRMFine->GetXaxis()->FindBin(xLow);
    Int_t binxUpFine = hRMFine->GetXaxis()->FindBin(xUp)-1;
    //cout << "xLowFine: " << hRMFine->GetXaxis()->GetBinLowEdge(binxLowFine) << "\txUpFine: " << hRMFine->GetXaxis()->GetBinUpEdge(binxUpFine) << endl;
    if(useFunctionWeight)
      normVector->SetAt(f1MergeFunction->Integral(xLow,xUp),ix-1);
    else
      normVector->SetAt(hRMFine->Integral(binxLowFine,binxUpFine,1,hRMFine->GetYaxis()->GetNbins()),ix-1);
    //if(fDebug) cout << "ix norm: " << normVector->At(ix-1) << endl;
  }

  Double_t content, oldcontent = 0.;
  Int_t ixNew = 0;
  Int_t iyNew = 0;
  Double_t xvalLo, xvalUp, yvalLo, yvalUp;
  Double_t xmin = hRM2->GetXaxis()->GetXmin();
  Double_t ymin = hRM2->GetYaxis()->GetXmin();
  Double_t xmax = hRM2->GetXaxis()->GetXmax();
  Double_t ymax = hRM2->GetYaxis()->GetXmax();
  for(int ix=1; ix<=hRMFine->GetXaxis()->GetNbins(); ix++) {
    xvalLo = hRMFine->GetXaxis()->GetBinLowEdge(ix);
    xvalUp = hRMFine->GetXaxis()->GetBinUpEdge(ix);
    if(xvalLo<xmin || xvalUp>xmax) continue;
    ixNew = hRM2->GetXaxis()->FindBin(hRMFine->GetXaxis()->GetBinCenter(ix));
    for(int iy=1; iy<=hRMFine->GetYaxis()->GetNbins(); iy++) {
      yvalLo = hRMFine->GetYaxis()->GetBinLowEdge(iy);
      yvalUp = hRMFine->GetYaxis()->GetBinUpEdge(iy);
      if(yvalLo<ymin || yvalUp>ymax) continue;
      content = hRMFine->GetBinContent(ix,iy);
      iyNew = hRM2->GetYaxis()->FindBin(hRMFine->GetYaxis()->GetBinCenter(iy));
      oldcontent = hRM2->GetBinContent(ixNew,iyNew);

      //if(fDebug) cout << "ixNew: " << ixNew << " " << xvalLo << " iyNew: " << iyNew << " " << yvalLo << " content: " << content << " oldContent: " << oldcontent << " newContent: " << oldcontent+content << endl;
      Double_t weight = 1.;
      if(normVector->At(ixNew-1)>0.) {
	if(useFunctionWeight)
	  weight = f1MergeFunction->Integral(xvalLo,xvalUp)/normVector->At(ixNew-1);
	else
	  weight = 1./normVector->At(ixNew-1);
      }
      hRM2->SetBinContent(ixNew,iyNew,oldcontent+content*weight);
    }
  }

  if(normVector) delete normVector;
  
  return hRM2;

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::CreateTruncated2DHisto(TH2 *h2, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {

  //
  // Limit axis range of 2D histogram
  //

  Int_t binMinXh2 = h2->GetXaxis()->FindBin(xmin);
  if(h2->GetXaxis()->GetBinLowEdge(binMinXh2) < xmin ) binMinXh2++;
  if(h2->GetXaxis()->GetBinLowEdge(binMinXh2) > xmin ) binMinXh2--;

  Int_t binMinYh2 = h2->GetYaxis()->FindBin(ymin);
  if(h2->GetYaxis()->GetBinLowEdge(binMinYh2) < ymin ) binMinYh2++;
  if(h2->GetYaxis()->GetBinLowEdge(binMinYh2) > ymin ) binMinYh2--;

  Int_t binMaxXh2 = h2->GetXaxis()->FindBin(xmax);
  if(h2->GetXaxis()->GetBinUpEdge(binMaxXh2) < xmax ) binMaxXh2++;
  if(h2->GetXaxis()->GetBinUpEdge(binMaxXh2) > xmax ) binMaxXh2--;

  Int_t binMaxYh2 = h2->GetYaxis()->FindBin(ymax);
  if(h2->GetYaxis()->GetBinUpEdge(binMaxYh2) < ymax ) binMaxYh2++;
  if(h2->GetYaxis()->GetBinUpEdge(binMaxYh2) > ymax ) binMaxYh2--;

  Int_t nbinsX = binMaxXh2-binMinXh2+1;
  Int_t nbinsY = binMaxYh2-binMinYh2+1;

  Double_t *binsX = new Double_t[nbinsX+1];
  Double_t *binsY = new Double_t[nbinsY+1];

  for(int ix=1; ix<=nbinsX; ix++)
    binsX[ix-1] = h2->GetXaxis()->GetBinLowEdge(binMinXh2+ix-1);
  binsX[nbinsX] = h2->GetXaxis()->GetBinUpEdge(binMaxXh2);

  for(int iy=1; iy<=nbinsY; iy++)
    binsY[iy-1] = h2->GetYaxis()->GetBinLowEdge(binMinYh2+iy-1);
  binsY[nbinsY] = h2->GetYaxis()->GetBinUpEdge(binMaxYh2);

  TH2 *h2Lim = new TH2D("h2Lim","h2Lim",nbinsX,binsX,nbinsY,binsY);

  for(int ix=1; ix<=nbinsX; ix++) {
    //    cout << "ix: " << ix << "  " << binsX[ix] << endl;
    for(int iy=1; iy<=nbinsY; iy++) {
      // cout << "ix: " << ix << "  " << binsX[ix] << "\tiy: " << iy << "  " << binsY[iy] << endl;

      double content = h2->GetBinContent(binMinXh2+ix-1,binMinYh2+iy-1);
      double error = h2->GetBinContent(binMinXh2+ix-1,binMinYh2+iy-1);
      h2Lim->SetBinContent(ix,iy,content);
      if(fbCalcErrors) h2Lim->SetBinError(ix,iy,error);

    }
  }



  return h2Lim;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::TruncateAxisRangeResponseMatrix(TH2 *hRMOrig, Double_t xmin, Double_t xmax, Double_t ymin, Double_t ymax) {

  //
  // Limit axis range of response matrix without changing normalization
  //

  //TH2 *hRMLimit
  //TH2 *hRMLimit2 = (TH2*)hRMLimit->Clone("hRMLimit2");

  TH2 *hRMLimit2 = CreateTruncated2DHisto(hRMOrig, xmin, xmax, ymin, ymax);
  hRMLimit2->Reset();

  double binCent[2] = {0.,0.}; 
  double content = 0.;
  double error = 0.;
  Int_t binOrig[2] = {0};
  for(int ix=1; ix<=hRMLimit2->GetXaxis()->GetNbins(); ix++) {
    binCent[0] = hRMLimit2->GetXaxis()->GetBinCenter(ix);
    binOrig[0] = hRMOrig->GetXaxis()->FindBin(binCent[0]);
    for(int iy=1; iy<=hRMLimit2->GetYaxis()->GetNbins(); iy++) {
      binCent[1] = hRMLimit2->GetYaxis()->GetBinCenter(iy);
      binOrig[1] = hRMOrig->GetYaxis()->FindBin(binCent[1]);

      content = hRMOrig->GetBinContent(binOrig[0],binOrig[1]);
      error = hRMOrig->GetBinError(binOrig[0],binOrig[1]);

      hRMLimit2->SetBinContent(ix,iy,content);
      if(fbCalcErrors) hRMLimit2->SetBinError(ix,iy,error);

    }
  }
  

  return hRMLimit2;

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
Bool_t AliAnaChargedJetResponseMaker::CheckInputForCombinedResponse() {

  Bool_t check = kTRUE;

  if(!fhDeltaPt)             check = kFALSE;
  if(!fh2DetectorResponse)   check = kFALSE;
  if(!fhEfficiencyDet)       check = kFALSE;
  if(!fh1MeasuredTruncated)  check = kFALSE;
  if(fPtMin<0. || fPtMax<0.) check = kFALSE;

  return check;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnaChargedJetResponseMaker::MakeResponseMatrixCombined(Int_t skipBins, Int_t binWidthFactor, Int_t extraBins, Bool_t bVariableBinning, Double_t ptmin) {

  //Check if all input histos are there otherwise break
  if(!CheckInputForCombinedResponse()) {
    printf("Not all input available..... \n");
    return;
  }
 
  // Make response bkg fluctuations
  MakeResponseMatrixJetsFineMerged(skipBins,binWidthFactor,extraBins,bVariableBinning,ptmin);

  //Get TH2D with dimensions we want final RM
  TH2D *h2ResponseMatrix = fResponseMatrix->Projection(0,1,"E");
  h2ResponseMatrix->Reset();

  //Get fine binned response matrix from bkg fluctuations
  TH2D *h2ResponseMatrixFine = fResponseMatrixFine->Projection(0,1,"E");

  // Rebin response detector effects
  const TArrayD *arrayX = h2ResponseMatrixFine->GetXaxis()->GetXbins();
  const Double_t *xbinsArray = arrayX->GetArray();

  TH2D *h2ResponseDetTmp = new TH2D("h2ResponseDetTmp","h2ResponseDetTmp",h2ResponseMatrixFine->GetNbinsX(),xbinsArray,h2ResponseMatrixFine->GetNbinsX(),xbinsArray);

  fh2DetectorResponseRebin = (TH2D*)MakeResponseMatrixRebin(fh2DetectorResponse,h2ResponseDetTmp,kFALSE);
  fh2DetectorResponseRebin->SetName("fh2DetectorResponseRebin");
  fh2DetectorResponseRebin->SetTitle("fh2DetectorResponseRebin");
  fh2DetectorResponseRebin->SetXTitle("p_{T}^{jet} gen");
  fh2DetectorResponseRebin->SetYTitle("p_{T}^{jet} rec");

  // Make full combined response matrix fine binning (background flucutuations + detector effects)
  fh2ResponseMatrixCombinedFineFull = (TH2D*)MultiplityResponseMatrices(h2ResponseMatrixFine,fh2DetectorResponseRebin);
  fh2ResponseMatrixCombinedFineFull->SetName("fh2ResponseMatrixCombinedFineFull");
  fh2ResponseMatrixCombinedFineFull->SetTitle("fh2ResponseMatrixCombinedFineFull");

  // Rebin full combined response matrix with weighting
  fh2ResponseMatrixCombinedFull = (TH2D*)MakeResponseMatrixRebin(fh2ResponseMatrixCombinedFineFull,h2ResponseMatrix,kTRUE); 
  fh2ResponseMatrixCombinedFull->SetName("fh2ResponseMatrixCombinedFull");
  fh2ResponseMatrixCombinedFull->SetTitle("fh2ResponseMatrixCombinedFull");
  fh2ResponseMatrixCombinedFull->SetXTitle("#it{p}_{T,gen}^{ch jet} (GeV/#it{c})");
  fh2ResponseMatrixCombinedFull->SetYTitle("#it{p}_{T,rec}^{ch jet} (GeV/#it{c})");

  fh2ResponseMatrixCombined = (TH2D*)TruncateAxisRangeResponseMatrix(fh2ResponseMatrixCombinedFull, h2ResponseMatrix->GetXaxis()->GetXmin(),h2ResponseMatrix->GetXaxis()->GetXmax(),fh1MeasuredTruncated->GetXaxis()->GetXmin(),fh1MeasuredTruncated->GetXaxis()->GetXmax());
  fh2ResponseMatrixCombined->SetTitle("fh2ResponseMatrixCombined");
  fh2ResponseMatrixCombined->SetName("fh2ResponseMatrixCombined");

  fhEfficiencyCombined = (TH1D*)fh2ResponseMatrixCombined->ProjectionX("fhEfficiencyCombined");
  fhEfficiencyCombined->Reset();

  for(int i=1; i<=fh2ResponseMatrixCombined->GetNbinsX(); i++) {
    Double_t sumyield = 0.;
    for(int j=1; j<=fh2ResponseMatrixCombined->GetYaxis()->GetNbins(); j++) {
      sumyield+=fh2ResponseMatrixCombined->GetBinContent(i,j);
    }
    Double_t kCounter = 0.;
    Double_t kPtLowEdge = fh2ResponseMatrixCombined->GetXaxis()->GetBinLowEdge(i);
    Double_t kPtUpEdge  = fh2ResponseMatrixCombined->GetXaxis()->GetBinUpEdge(i);
    Double_t kJetEffDet = 0.;

    for(int k=1; k<=fhEfficiencyDet->GetNbinsX(); k++) {
      if(fhEfficiencyDet->GetXaxis()->GetBinLowEdge(k)>=kPtLowEdge && fhEfficiencyDet->GetXaxis()->GetBinUpEdge(k)<=kPtUpEdge) {
	kJetEffDet+=fhEfficiencyDet->GetBinContent(k);
	kCounter+=1.;
      }
    }
    kJetEffDet = kJetEffDet/kCounter;

    if(kJetEffDet==0.) kJetEffDet=1.;

    fhEfficiencyCombined->SetBinContent(i,sumyield*kJetEffDet);

  }

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::MultiplityResponseMatrices(TH2 *h2RMDeltaPt, TH2 *h2RMDetector) {

  //
  // Make combined response matrix (background flucutuations + detector effects)
  //
  // hRMDeltaPt is the response matrix for background fluctuations from \delta(p_t) measurement
  // hRMDetector is the response matrix for detector effects: needs to be a squared matrix with on each axis the same bins as on the generated axis of the bkg fluctuations response matrix
  //
  // Function assumes that generated/unfolded axis is x-axis and reconstructed is on y-axis on both respone matrices


  TH2D *h2ResponseMatrixCombined = (TH2D*)h2RMDeltaPt->Clone("h2ResponseMatrixCombined"); //h2RMDeltaPt is the bkg fluctuations RM which has the dimensions we want for the combined response matrix
  h2ResponseMatrixCombined->SetTitle("h2ResponseMatrixCombined");
  h2ResponseMatrixCombined->SetName("h2ResponseMatrixCombined");

  // M = RM_deltaPt * RM_DetEffects * T   (M=measured T=truth)
  // Dimensions:
  // mx1 = mxd * dxt * tx1
  // m = measured bins
  // t = truth bins
  // d = rec for RM_DetEffects and gen for RM_deltaPt
  // RM_comb = RM_deltaPt * RM_DetEffects (dimensions mxt)
  // RM_comb(m,t) = Sum_d RM_deltaPt(m,d)*RM_DetEffects(d,t)

  if(fDebug) {
    printf("Nt=%d\n",h2ResponseMatrixCombined->GetNbinsX());
    printf("Nm=%d\n",h2ResponseMatrixCombined->GetNbinsY());
    printf("Nd=%d\n",h2RMDeltaPt->GetNbinsX());
  }


  for(Int_t t=1; t<=h2ResponseMatrixCombined->GetNbinsX();t++) { 
    for(Int_t m=1; m<=h2ResponseMatrixCombined->GetNbinsY();m++) { 
      Double_t valueSum = 0.;    
      for(Int_t d=1; d<=h2RMDeltaPt->GetNbinsX();d++) { 
	valueSum += h2RMDeltaPt->GetBinContent(d,m) * h2RMDetector->GetBinContent(t,d);
	//	if(t==10 && m==10) cout << "sum m,d=" << m << "," << d << endl;
      }//d-loop
      //  if(t==10) cout << "t,m = " << t << "," << m << "\tvalueSum: " << valueSum << endl; 
      h2ResponseMatrixCombined->SetBinContent(t,m,valueSum);
    } //m-loop
  }//t-loop

  return h2ResponseMatrixCombined;

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::GetTransposeResponsMatrix(TH2 *h2RM) {
  //
  // Transpose response matrix
  //

  Printf("AliAnaChargedJetResponseMaker::GetTransposeResponsMatrix");

  //Initialize transposed response matrix h2RMTranspose
  // TArrayD *arrayX = (TArrayD*)h2RM->GetXaxis()->GetXbins();
  // for(int i=0; i<arrayX->GetSize(); i++) Printf("i: %d  %f", i,arrayX->At(i));
  // Double_t *xbinsArray = arrayX->GetArray();

  // TArrayD *arrayY = (TArrayD*)h2RM->GetYaxis()->GetXbins();
  // for(int i=0; i<arrayY->GetSize(); i++) Printf("i: %d  %f", i,arrayY->At(i));
  // Double_t *ybinsArray = arrayY->GetArray();


  Double_t *ybinsArray = new Double_t[h2RM->GetNbinsY()+1];
  for(int i=1; i<=h2RM->GetNbinsY(); i++) {
    ybinsArray[i-1] = h2RM->GetYaxis()->GetBinLowEdge(i);
    Printf("i=%d  %f   %f",i,ybinsArray[i-1],h2RM->GetYaxis()->GetBinLowEdge(i));
  }
  ybinsArray[h2RM->GetNbinsY()] = h2RM->GetYaxis()->GetBinUpEdge(h2RM->GetNbinsY());

  Double_t *xbinsArray = new Double_t[h2RM->GetNbinsX()+1];
  for(int i=1; i<=h2RM->GetNbinsX(); i++) 
    xbinsArray[i-1] = h2RM->GetXaxis()->GetBinLowEdge(i);
  xbinsArray[h2RM->GetNbinsX()] = h2RM->GetXaxis()->GetBinUpEdge(h2RM->GetNbinsX());


  TH2D *h2RMTranspose = new TH2D("h2RMTranspose","h2RMTranspose",h2RM->GetNbinsY(),ybinsArray,h2RM->GetNbinsX(),xbinsArray);

  //Fill tranposed response matrix
  for (Int_t ibin = 1; ibin <= h2RMTranspose->GetNbinsX(); ibin++) {
    for (Int_t jbin = 1; jbin <= h2RMTranspose->GetNbinsY(); jbin++) {
      h2RMTranspose->SetBinContent(ibin,jbin, h2RM->GetBinContent(jbin,ibin));
    }
  }


  return h2RMTranspose;

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH2* AliAnaChargedJetResponseMaker::NormalizeResponsMatrixYaxisWithPrior(TH2 *h2RM, TH1 *hPrior) {
  //
  // Normalize such that the Y projection is the prior
  //

  //  TH1D *hProjRespY = (TH1D*)h2RM->ProjectionY("hProjRespY");
  double intPrior = hPrior->Integral();//"width");
  for (Int_t jbin = 1; jbin <= h2RM->GetNbinsY(); jbin++) {
    //    double corr = hPrior->GetBinContent(jbin)/hProjRespY->GetBinContent(jbin);
      for (Int_t ibin = 1; ibin <= h2RM->GetNbinsX(); ibin++) {
	double content = h2RM->GetBinContent(ibin,jbin);
	//	h2RM->SetBinContent(ibin,jbin,content*corr);
	h2RM->SetBinContent(ibin,jbin,hPrior->GetBinContent(jbin)/hPrior->GetBinWidth(jbin)/intPrior*content);
    }
  }

  return h2RM;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* AliAnaChargedJetResponseMaker::MultiplyResponseGenerated(TH1 *hGen, TH2 *hResponse,TH1 *hEfficiency,Bool_t bDrawSlices) {

  //
  // Multiply hGen with response matrix to obtain refolded spectrum
  // Efficiency must be given. In case efficiency is 1 use SetFlatEfficiency(1.) before calling function
  //

  if(!hEfficiency) {
    //    printf("Efficiency must be given. In case efficiency is 1 use SetFlatEfficiency(1.) before calling function. Aborting!");
    //    return 0;
    printf("Setting efficiency to 1 \n");
    hEfficiency = (TH1D*)hGen->Clone("hEfficiency");
    hEfficiency->Reset();
    for(int i=1; i<=hEfficiency->GetNbinsX(); i++) hEfficiency->SetBinContent(i,1.);
  }

  //For response
  //x-axis: generated
  //y-axis: reconstructed
  if(fDebug>0) cout << "nbins hGen: " << hGen->GetNbinsX() << "\t nbins hResponseGen: " << hResponse->GetXaxis()->GetNbins() << "\t nbins hResponseRec: " << hResponse->GetYaxis()->GetNbins()  << endl;

  TH1D *hRec = hResponse->ProjectionY("hRec");
  hRec->Sumw2();
  hRec->Reset();
  hRec->SetTitle("hRec");
  hRec->SetName("hRec");

  for(int irec=1; irec<=hRec->GetNbinsX(); irec++)
    hRec->SetBinContent(irec,0);

  TH1D *hRecGenBin = 0x0;
  TCanvas *cSlices = 0x0;
  if(bDrawSlices) {
    cSlices = new TCanvas("cSlices","cSlices:Slices",400,400);
    gPad->SetLogy();
  }

  Double_t yieldMC = 0.;
  Double_t yieldMCerror = 0.;
  Double_t sumYield = 0.;
  const Int_t nbinsRec = hRec->GetNbinsX()+1;
  Double_t sumError2[999] = {0.};
  Double_t eff = 0.;

  for(int igen=1; igen<=hGen->GetNbinsX(); igen++) {
    //get pTMC
    sumYield = 0.;
    if(fEffFlat>0.)
      eff = fEffFlat;
    else
      eff = hEfficiency->GetBinContent(igen);
    yieldMC = hGen->GetBinContent(igen)*eff;
    yieldMCerror = hGen->GetBinError(igen)*eff;

    if(bDrawSlices) {
      hRecGenBin = hResponse->ProjectionY(Form("hRecGenBin%d",igen));
      hRecGenBin->Sumw2();
      hRecGenBin->Reset();
      hRecGenBin->SetTitle(Form("hRecGenBin%d",igen));
      hRecGenBin->SetName(Form("hRecGenBin%d",igen));
    }

    for(int irec=1; irec<=hRec->GetNbinsX(); irec++) {
      hRec->AddBinContent(irec,yieldMC*hResponse->GetBinContent(igen,irec));
      sumYield+=hResponse->GetBinContent(igen,irec);
      //      Double_t A = yieldMC*hResponse->GetBinError(igen,irec);
      Double_t B = hResponse->GetBinContent(igen,irec)*yieldMCerror;
      //      Double_t tmp2 = A*A + B*B;
      //sumError2[irec-1] += tmp2 ;
      sumError2[irec-1] += B*B;

      if(bDrawSlices)
	hRecGenBin->SetBinContent(irec,yieldMC*hResponse->GetBinContent(igen,irec));

    }
    if(bDrawSlices) {
      cSlices->cd();
      hRecGenBin->SetLineColor(igen);
      if(igen==1) hRecGenBin->DrawCopy();      
      else hRecGenBin->DrawCopy("same");
    }

    if(hRecGenBin) delete hRecGenBin;
    
    cout << "igen: " << igen << "\tpTMC: " << hGen->GetXaxis()->GetBinCenter(igen) << "\teff:" << eff << "\tsumYield: " << sumYield << endl;
  }
  
  for(int i=0; i<nbinsRec; i++) {
    if(sumError2[i]>0.)
      hRec->SetBinError(i+1,TMath::Sqrt(sumError2[i]));
  }


  return hRec;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
TH1D* AliAnaChargedJetResponseMaker::MultiplyResponseGenerated(TF1 *fGen, TH2 *hResponse,TH1 *hEfficiency) {

  //
  // Multiply fGen function with response matrix to obtain (re)folded spectrum
  // Efficiency must be given. In case efficiency is 1 use SetFlatEfficiency(1.) before calling function
  //

  //For response
  //x-axis: generated
  //y-axis: reconstructed

  if(fDebug>0) printf(">>AliAnaChargedJetResponseMaker::MultiplyResponseGenerated(TF1 *fGen, TH2 *hResponse,TH1 *hEfficiency)\n");

  TH1D *hRec = hResponse->ProjectionY("hRec");
  hRec->Sumw2();
  hRec->Reset();
  hRec->SetTitle("hRec");
  hRec->SetName("hRec");

  //  TH1D *hRec = new TH1D("hRec","hRec",hResponse->GetNbinsY(),hResponse->GetYaxis()->GetXmin(),hResponse->GetYaxis()->GetXmax());
  
  for(int irec=1; irec<=hRec->GetNbinsX(); irec++)
    hRec->SetBinContent(irec,0);
  
  Double_t yieldMC = 0.;
  Double_t sumYield = 0.;
  Double_t eff = 0.;
  for(int igen=1; igen<=hResponse->GetNbinsX(); igen++) {
    //get pTMC
    sumYield = 0.;
    double pTMC = hResponse->GetXaxis()->GetBinCenter(igen);
    if(hEfficiency) { 
      int binEff = hEfficiency->FindBin(pTMC);
      eff = hEfficiency->GetBinContent(binEff);
    }
    else eff = 1.;
    // yieldMC = fGen->Eval(pTMC)*eff;
    yieldMC = fGen->Integral(hResponse->GetXaxis()->GetBinLowEdge(igen),hResponse->GetXaxis()->GetBinUpEdge(igen))*eff;
    for(int irec=1; irec<=hResponse->GetNbinsY(); irec++) {
      hRec->AddBinContent(irec,yieldMC*hResponse->GetBinContent(igen,irec));
      sumYield+=hResponse->GetBinContent(igen,irec);
    }
    //    cout << "igen: " << igen << "\tpTMC: " << pTMC << "\tsumYield: " << sumYield << endl;
    //    cout << "yieldMC: " << yieldMC << endl;
    
  }

  return hRec;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------
/* static */ Double_t AliAnaChargedJetResponseMaker::GetBetaPerDOFValue(Int_t betaColl, Int_t betaOpt) {

  Float_t betaPerDOFoptions[6] = {0.};
  //define beta per degree of freedom (DOF=nbins in measured)
  if(betaColl==0) {
    betaPerDOFoptions[0] = 9.1e-9;
    betaPerDOFoptions[1] = 2.25e-8;
    betaPerDOFoptions[2] = 4.55e-8;
    betaPerDOFoptions[3] = 9.1e-8;
    betaPerDOFoptions[4] = 2.25e-7;
    betaPerDOFoptions[5] = 4.55e-7;
  }
  if(betaColl==1) {
    betaPerDOFoptions[0] = 9.1e-7;
    betaPerDOFoptions[1] = 2.25e-6;
    betaPerDOFoptions[2] = 4.55e-6;
    betaPerDOFoptions[3] = 9.1e-6;
    betaPerDOFoptions[4] = 2.25e-5;
    betaPerDOFoptions[5] = 4.55e-5;
  }
  if(betaColl==2) {
    betaPerDOFoptions[0] = 9.1e-5;
    betaPerDOFoptions[1] = 2.25e-4;
    betaPerDOFoptions[2] = 4.55e-4;
    betaPerDOFoptions[3] = 9.1e-4;
    betaPerDOFoptions[4] = 2.25e-3;
    betaPerDOFoptions[5] = 4.55e-3;
  }
  if(betaColl==3) {
    betaPerDOFoptions[0] = 9.1e-3;
    betaPerDOFoptions[1] = 2.25e-2;
    betaPerDOFoptions[2] = 4.55e-2;
    betaPerDOFoptions[3] = 9.1e-2;
    betaPerDOFoptions[4] = 2.25e-1;
    betaPerDOFoptions[5] = 4.55e-1;
  }
  if(betaColl==4) {
    betaPerDOFoptions[0] = 9.1e-1;
    betaPerDOFoptions[1] = 2.25;
    betaPerDOFoptions[2] = 4.55;
    betaPerDOFoptions[3] = 9.1;
    betaPerDOFoptions[4] = 22.5;
    betaPerDOFoptions[5] = 45.5;
  }

  return betaPerDOFoptions[betaOpt];

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnaChargedJetResponseMaker::InterpolateFast(TGraph *gr, Double_t x) {

  Double_t ipmax = gr->GetN()-1.;
  Double_t x2,y2,xold,yold;

  Double_t xmin,ymin,xmax, ymax;
  gr->GetPoint(0,xmin,ymin);
  gr->GetPoint(gr->GetN()-1,xmax,ymax);
  if(x<xmin || x>xmax) return 0;

  Double_t ip = ipmax/2.;
  Double_t ipold = 0.;
  gr->GetPoint((int)(ip),x2,y2);

  Int_t i = 0;

  if(x2>x) {
    while(x2>x) { 
      if(i==0) ipold=0.;
      ip -= (ip)/2.;
      gr->GetPoint((int)(ip),x2,y2);
      if(x2>x){}
      else ipold = ip;
      i++;
      //      cout << "ipold: " << ipold << "\tip: " << ip << "\tx2: " << x2 << "\ty2: " << y2 << endl;
    }
  }
  else if(x2<x) {
    while(x2<x) {
      if(i==0) ipold=ipmax;
      ip += (ipold-ip)/2.;
      gr->GetPoint((int)(ip),x2,y2);
      if(x2>x) ipold = ip;
      else {}
      i++;
      //      cout << "ipold: " << ipold << "\tip: " << ip << "\tx2: " << x2 << "\ty2: " << y2 << endl;
    }
  }
  
  int ip2 = 0;
  if(x2>x) {
    ip2 = (int)(ip)-1;
    gr->GetPoint(ip2,x2,y2);
    while(x2>x) {
      ip2--;
      gr->GetPoint(ip2,x2,y2);
    }
    gr->GetPoint(ip2+1,xold,yold);
  }
  else {
    ip2 = (int)(ip)+1;
    gr->GetPoint(ip2,x2,y2);
    while(x2<x) {
      ip2++;
      gr->GetPoint(ip2,x2,y2);
    }
    gr->GetPoint(ip2-1,xold,yold);

  }
  //  cout << "For x=" << x << " interpolate between: " << xold << " and " << x2 << endl;
  return ((x-xold)*y2 + (x2-x)*yold) / (x2-xold);

}

//--------------------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnaChargedJetResponseMaker::InterpolateFast(TH1 *h, Double_t x) {

  // Double_t ipmax = gr->GetN()-1.;
  Double_t ipmax = (double)h->GetNbinsX();
  Double_t x2,y2,xold,yold;

  Double_t xmin = h->GetXaxis()->GetXmin();
  Double_t xmax = h->GetXaxis()->GetXmax();
  if(x<xmin || x>xmax) return 0;

  Double_t ip = ipmax/2.;
  Double_t ipold = 0.;

  x2 = h->GetBinCenter((int)ip);
  y2 = h->GetBinContent((int)ip);

  Int_t i = 0;

  if(x2>x) {
    while(x2>x) { 
      if(i==0) ipold=0.;
      ip -= (ip)/2.;
      x2 = h->GetBinCenter((int)ip);
      y2 = h->GetBinContent((int)ip);
      if(x2>x) {}
      else ipold = ip;
      i++;
      //      cout << "ipold: " << ipold << "\tip: " << ip << "\tx2: " << x2 << "\ty2: " << y2 << endl;
    }
  }
  else if(x2<x) {
    while(x2<x) {
      if(i==0) ipold=ipmax;
      ip += (ipold-ip)/2.;
      x2 = h->GetBinCenter((int)ip);
      y2 = h->GetBinContent((int)ip);
      if(x2>x) ipold = ip;
      else {}
      i++;
      //      cout << "ipold: " << ipold << "\tip: " << ip << "\tx2: " << x2 << "\ty2: " << y2 << endl;
    }
  }
  
  int ip2 = 0;
  if(x2>x) {
    ip2 = (int)(ip)-1;
    x2 = h->GetBinCenter(ip2);
    y2 = h->GetBinContent(ip2);
    while(x2>x) {
      ip2--;
      x2 = h->GetBinCenter(ip2);
      y2 = h->GetBinContent(ip2);
    }
    xold = h->GetBinCenter(ip2+1);
    yold = h->GetBinContent(ip2+1);
  }
  else {
    ip2 = (int)(ip)+1;
    x2 = h->GetBinCenter(ip2);
    y2 = h->GetBinContent(ip2);
    while(x2<x) {
      ip2++;
      x2 = h->GetBinCenter(ip2);
      y2 = h->GetBinContent(ip2);
    }
    xold = h->GetBinCenter(ip2-1);
    yold = h->GetBinContent(ip2-1);

  }
  //  cout << "For x=" << x << " interpolate between: " << xold << " and " << x2 << endl;
  return ((x-xold)*y2 + (x2-x)*yold) / (x2-xold);

}


