//------------------------------------------------------------------------------
// Implementation of AliTPCPerformanceSummary class. 
// It has only static member functions to extract some TPC Performance
// parameters and produce trend graphs.
// The function MakeReport is to be called for every run. It reads AliPerformanceTPC
// and AliPerformanceDEdx objects from file and produces a
// rootfile with the results stored in a TTree.
// The function MakeReport needs a list of these rootfiles as input
// and writes the output (tree and histograms) to another rootfile.
//
// by M.Knichel 15/10/2010
//
// Modifications:
//    by Marian Ivanov, m.ivanov@cern.ch;
//------------------------------------------------------------------------------

#include <fstream>

#include "TSystem.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TList.h"
#include "TFile.h"
#include "TGrid.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "THnSparse.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TError.h"
#include <iostream>
#include "AliGRPObject.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "TTreeStream.h"
#include "AliPerformanceTPC.h"
#include "AliPerformanceDEdx.h"
#include "AliPerformanceDCA.h"
#include "AliPerformanceMatch.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
//
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPManager.h"
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCdataQA.h"
#include "AliTPCCalROC.h"
#include "TGeoGlobalMagField.h"
#include "AliMathBase.h"
#include "AliTPCPerformanceSummary.h"
#include "AliDAQ.h"

using std::ifstream;

ClassImp(AliTPCPerformanceSummary)

Bool_t AliTPCPerformanceSummary::fgForceTHnSparse = kFALSE;

Bool_t AliTPCPerformanceSummary::GetStatInfo(TH1 * histo, TVectorF &statInfo, Int_t axis){
  //
  // fill basic statistical information
  // 
  if (!histo){
    return kFALSE;
  }
  statInfo[0]=histo->GetEntries();
  statInfo[1]=histo->GetMean(axis);
  statInfo[2]=histo->GetMeanError(axis);
  statInfo[3]=histo->GetRMS(axis);
  statInfo[4]=histo->GetRMSError(axis);
  return kTRUE;
}

Bool_t AliTPCPerformanceSummary::GetFitInfo(TF1 * fitFunction, TVectorF &fitInfo){
  //
  //  fill basic statistical information
  //  parameters, covariance, chi2
  if (!fitFunction) return kFALSE;
  Int_t npar=fitFunction->GetNpar();
  if (2*npar+1>fitInfo.GetNrows()) return kFALSE;

  for (Int_t ipar=0; ipar<npar; ipar++){
    fitInfo[ipar]=fitFunction->GetParameter(ipar);
    fitInfo[npar+ipar]=fitFunction->GetParError(ipar);
  }
  if (fitFunction->GetNDF()>0.) fitInfo[2*npar]=fitFunction->GetChisquare()/fitFunction->GetNDF();
  return kTRUE;
}

//_____________________________________________________________________________
void AliTPCPerformanceSummary::WriteToTTreeSRedirector(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const AliPerformanceMatch* pTPCMatch,const AliPerformanceMatch* pTPCPull, const AliPerformanceMatch* pConstrain, TTreeSRedirector* const pcstream, Int_t run)
{
  // 
  // Extracts performance parameters from pTPC and pTPCgain.
  // Output is written to pcstream.
  // The run number must be provided since it is not stored in 
  // AliPerformanceTPC or AliPerformanceDEdx.
  // Here we assume that 
  //
  if (run <= 0 ) {
    if (pTPCMatch) {run = pTPCMatch->GetRunNumber(); }
    if (pTPCgain) {run = pTPCgain->GetRunNumber(); }
    if (pTPC) { run = pTPC->GetRunNumber(); }
  }

  // check the presence of the detectors
  AliCDBEntry *entry=0x0;
  try {
    entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  } catch(...) {
    Info("AliTPCPerformanceSummary::WriteToTTreeSRedirector","No GRP entry found");
    entry = 0x0;
  }

  if (!entry) return;

  AliGRPObject* grpData = dynamic_cast<AliGRPObject*>(entry->GetObject());
  if (!grpData) {
    Info("AliTPCPerformanceSummary::WriteToTTreeSRedirector","Failed to get GRP data for run %d\n",run);
    return;
  }
  Int_t activeDetectors = grpData->GetDetectorMask();
  TString detStr = AliDAQ::ListOfTriggeredDetectors(activeDetectors);
  //printf("Detectors in the data:\n%s\n",detStr.Data());
  if ( detStr.Contains("TPC")==0){
    Info("AliTPCPerformanceSummary::WriteToTTreeSRedirector","TPC not present in run %d", run);
    return;
  }

  TObjString runType;
  
  Int_t startTimeGRP=0;
  Int_t stopTimeGRP=0;   
  Int_t time=0;
  Int_t duration=0;
  Float_t currentL3 =0;
  Int_t polarityL3 = 0;
  Float_t bz = 0;
  if (AliCDBManager::Instance()->GetRun()==run){
    AliTPCcalibDB     *calibDB=0;
    calibDB = AliTPCcalibDB::Instance();         
    if (calibDB->GetGRP(run)){
      startTimeGRP = AliTPCcalibDB::GetGRP(run)->GetTimeStart();
      stopTimeGRP  = AliTPCcalibDB::GetGRP(run)->GetTimeEnd();
      currentL3 = AliTPCcalibDB::GetL3Current(run);
      polarityL3 = AliTPCcalibDB::GetL3Polarity(run);
      bz = AliTPCcalibDB::GetBz(run);
      if (polarityL3>0) bz*=-1; 
      runType = AliTPCcalibDB::GetRunType(run).Data();  
    }    
  }
  
  time = startTimeGRP;
  duration = (stopTimeGRP-startTimeGRP);
  TObjString period(gSystem->Getenv("eperiod"));
  TObjString pass(gSystem->Getenv("epass"));
  TObjString dataType(gSystem->Getenv("edataType"));
  ::Info("AliTPCPerformanceSummary::WriteToTTreeSRedirector","%s/%s/%s",dataType.GetName(), period.GetName(), pass.GetName());
  Int_t year=0;
  if (gSystem->Getenv("eyear")) year=atoi(gSystem->Getenv("eyear"));
    if (!pcstream) return;
    (*pcstream)<<"tpcQA"<<      
      "run="<<run<<
      "time="<<time<<
      "year="<<year<<
      "period.="<<&period<<
      "pass.="<<&pass<<
      "dataType.="<<&dataType<<
      "startTimeGRP="<<startTimeGRP<<
      "stopTimeGRP="<<stopTimeGRP<<
      "duration="<<duration<<
      "bz="<<bz<<
      "runType.="<<&runType;
    if (pTPC) {
        if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(9)->SetRangeUser(0.5,1.5);
        if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0.25,10);
        if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-1,1);
        AnalyzeNCL(pTPC, pcstream);
        MakeRawOCDBQAPlot(pcstream);
        AnalyzeDrift(pTPC, pcstream);
        AnalyzeDriftPos(pTPC, pcstream);
        AnalyzeDriftNeg(pTPC, pcstream);
        AnalyzeDCARPhi(pTPC, pcstream);
        AnalyzeDCARPhiPos(pTPC, pcstream);
        AnalyzeDCARPhiNeg(pTPC, pcstream);
        AnalyzeEvent(pTPC, pcstream);
        AnalyzePt(pTPC,pcstream);
        AnalyzeChargeOverPt(pTPC,pcstream);
        AnalyzeQAPosNegDpT(pTPC,pcstream);
        AnalyzeQADCAFitParameter(pTPC,pcstream);
        AnalyzeOcc(pTPC,pcstream);
    }
    AnalyzeGain(pTPCgain, pcstream);
    AnalyzeMatch(pTPCMatch, pcstream);
    AnalyzePull(pTPCPull, pcstream);
    AnalyzeConstrain(pConstrain, pcstream);

    (*pcstream)<<"tpcQA"<<"\n";
    TTree * tree = ((*pcstream)<<"tpcQA").GetTree();
    tree->SetAlias("nEvents","entriesMult");
    
}

//_____________________________________________________________________________
void AliTPCPerformanceSummary::WriteToFile(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const AliPerformanceMatch* pMatch,  const AliPerformanceMatch* pPull, const AliPerformanceMatch* pConstrain, const Char_t* outfile, Int_t run)
{
    //
    // Extracts performance parameters from pTPC and pTPCgain.
    // Output is written to a TTree saved in outfile.
    // The run number must be provided since it is not stored in 
    // AliPerformanceTPC or AliPerformanceDEdx.
    //
    // The function creates a TTreeSRedirector and calls the 
    // function WriteToTTreeSRedirector.
    //
    
    if (!outfile) return;
    TTreeSRedirector* pcstream = 0;
    pcstream = new TTreeSRedirector(outfile,"recreate");
    if (!pcstream) return;
    WriteToTTreeSRedirector(pTPC, pTPCgain, pMatch, pPull, pConstrain, pcstream, run);
    if (pcstream) { delete pcstream; pcstream = 0; }    
    
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::MakeReport(const Char_t* infile, const Char_t* outfile, Int_t run)
{
    //
    // Reads QA information (AliPerformanceTPC and AliPerformanceDEdx) from
    // infile (this must be a rootfile) and writes the output to a TTree
    // stored in outfile.
    // The run number must be provided since it is not stored in 
    // AliPerformanceTPC or AliPerformanceDEdx.
    // 
    // The input objects must be named "AliPerformanceTPC" and 
    // "AliPerformanceDEdxTPCInner" and stored in a TList which name must
    // be one of the following: "TPC", "TPCQA", "TPC_PerformanceQA"
    // or "TPC_PerformanceQA/TPC" (with directory)
    //
    
    if (!infile) return -1;
    if (!outfile) return -1;
    TFile *f =0;
    f=TFile::Open(infile,"read");
    if (!f) {
        printf("File %s not available\n", infile);
        return -1;
    } 
    TList* list = 0;
    list = dynamic_cast<TList*>(f->Get("TPC")); 
    if (!list) { list = dynamic_cast<TList*>(f->Get("TPCQA")); }
    if (!list) { list = dynamic_cast<TList*>(f->Get("TPC_PerformanceQA/TPCQA")); }
    if (!list) { list = dynamic_cast<TList*>(f->Get("TPC_PerformanceQA")); }
    if (!list) { list = dynamic_cast<TList*>(f->Get("ITSTPCMatch")); }
    if (!list) {
            printf("QA %s not available\n", infile);
            return -1;
    } 
    AliPerformanceTPC* pTPC = 0;
    AliPerformanceDEdx* pTPCgain = 0; 
    AliPerformanceMatch* pTPCmatch = 0; 
    AliPerformanceMatch* pTPCPull = 0; 
    AliPerformanceMatch* pConstrain = 0;
    
    if (list) {  pTPC = dynamic_cast<AliPerformanceTPC*>(list->FindObject("AliPerformanceTPC")); }
    if (list) {  pTPCgain = dynamic_cast<AliPerformanceDEdx*>(list->FindObject("AliPerformanceDEdxTPCInner")); }
    if (list) {  pTPCmatch = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchTPCITS")); }
    if (list) {  pTPCPull = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchITSTPC")); }
    if (list) {  pConstrain = dynamic_cast<AliPerformanceMatch*>(list->FindObject("AliPerformanceMatchTPCConstrain")); }
    
    Int_t returncode = 0;
    WriteToFile(pTPC, pTPCgain, pTPCmatch , pTPCPull, pConstrain, outfile, run);
    if (f) { delete f; f=0; }
    return returncode;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::ProduceTrends(const Char_t* infilelist, const Char_t* outfile)
{
    //
    // Produces trend graphs.
    //
    // Input: infilelist is a textfile with one rootfile per line.
    // There should be one rootfile for each run, the rootfile must
    // contain the output of the MakeReport function
    // Output: the information for all runs is merged into a TTree
    // that is saved in outfile along with the trend graphs.
    // Trend graphs are stored as TCanvas objects to include axis labels etc.
    //
    
    if (!infilelist) return -1;
    if (!outfile) return -1;
     
    TChain* chain = new TChain("tpcQA");
    if(!chain) return -1;

    ifstream in;
    in.open(infilelist);

    TString currentFile;    
    while(in.good()) {
        in >> currentFile;

        if (!currentFile.Contains("root")) continue; // protection            
        chain->Add(currentFile.Data());
    }
    in.close();
    //TTree *tree = chain;
    TTree *tree = chain->CopyTree("1");
    if(!tree) return -1;
    if (chain) { delete chain; chain=0; }
    //TGraph* graph = dynamic_cast<TGraph*>(tree->DrawClone("run:run"));
    //TGraph *graph = (TGraph*)gPad->GetPrimitive("Graph");
    
    TFile* out = new TFile(outfile,"RECREATE");
    if(!out) return -1;

    out->cd();
    const Char_t* condition = "meanTPCncl>0";
    SaveGraph(tree,"meanTPCnclF","run",condition);
    SaveGraph(tree,"rmsTPCnclF","run",condition);
    SaveGraph(tree,"meanTPCChi2","run",condition);
    SaveGraph(tree,"rmsTPCChi2","run",condition);
    SaveGraph(tree,"slopeATPCnclF","run",condition);
    SaveGraph(tree,"slopeCTPCnclF","run",condition);
    SaveGraph(tree,"slopeATPCnclFErr","run",condition);
    SaveGraph(tree,"slopeCTPCnclFErr","run",condition);
    SaveGraph(tree,"meanTPCncl","run",condition);
    SaveGraph(tree,"rmsTPCncl","run",condition);
    SaveGraph(tree,"slopeATPCncl","run",condition);
    SaveGraph(tree,"slopeCTPCncl","run",condition);
    SaveGraph(tree,"slopeATPCnclErr","run",condition);
    SaveGraph(tree,"slopeCTPCnclErr","run",condition);
    
    SaveGraph(tree,"offsetdRA","run",condition);
    SaveGraph(tree,"slopedRA","run",condition);
    SaveGraph(tree,"offsetdRC","run",condition);
    SaveGraph(tree,"slopedRC","run",condition);
    SaveGraph(tree,"offsetdRAErr","run",condition);
    SaveGraph(tree,"slopedRAErr","run",condition);    
    SaveGraph(tree,"offsetdRCErr","run",condition);
    SaveGraph(tree,"slopedRCErr","run",condition);
    SaveGraph(tree,"offsetdRAchi2","run",condition);
    SaveGraph(tree,"slopedRAchi2","run",condition);
    SaveGraph(tree,"offsetdRCchi2","run",condition);    
    SaveGraph(tree,"slopedRCchi2","run",condition);    
    
    SaveGraph(tree,"offsetdRAPos","run",condition);
      
    SaveGraph(tree,"slopedRAPos","run",condition);
    SaveGraph(tree,"offsetdRCPos","run",condition);
    SaveGraph(tree,"slopedRCPos","run",condition);
    SaveGraph(tree,"offsetdRAErrPos","run",condition);
    SaveGraph(tree,"slopedRAErrPos","run",condition);
    SaveGraph(tree,"offsetdRCErrPos","run",condition); 
    SaveGraph(tree,"slopedRCErrPos","run",condition);
    SaveGraph(tree,"offsetdRAchi2Pos","run",condition);
    SaveGraph(tree,"slopedRAchi2Pos","run",condition);
    SaveGraph(tree,"offsetdRCchi2Pos","run",condition);
    SaveGraph(tree,"slopedRCchi2Pos","run",condition);
        
    SaveGraph(tree,"offsetdRANeg","run",condition);
    SaveGraph(tree,"slopedRANeg","run",condition);
    SaveGraph(tree,"offsetdRCNeg","run",condition);
    SaveGraph(tree,"slopedRCNeg","run",condition);
    SaveGraph(tree,"offsetdRAErrNeg","run",condition);
    SaveGraph(tree,"slopedRAErrNeg","run",condition);
    SaveGraph(tree,"offsetdRCErrNeg","run",condition);
    SaveGraph(tree,"slopedRCErrNeg","run",condition);
    SaveGraph(tree,"offsetdRAchi2Neg","run",condition);
    SaveGraph(tree,"slopedRAchi2Neg","run",condition);
    SaveGraph(tree,"offsetdRCchi2Neg","run",condition);
    SaveGraph(tree,"slopedRCchi2Neg","run",condition);
        
    SaveGraph(tree,"offsetdZAPos","run",condition);
    SaveGraph(tree,"slopedZAPos","run",condition);
    SaveGraph(tree,"offsetdZCPos","run",condition);
    SaveGraph(tree,"slopedZCPos","run",condition);
    SaveGraph(tree,"offsetdZAErrPos","run",condition);
    SaveGraph(tree,"slopedZAErrPos","run",condition);
    SaveGraph(tree,"offsetdZCErrPos","run",condition);
    SaveGraph(tree,"slopedZCErrPos","run",condition);
    SaveGraph(tree,"offsetdZAchi2Pos","run",condition);
    SaveGraph(tree,"slopedZAchi2Pos","run",condition);
    SaveGraph(tree,"offsetdZCchi2Pos","run",condition);
    SaveGraph(tree,"slopedZCchi2Pos","run",condition);
    
    SaveGraph(tree,"offsetdZANeg","run",condition);
    SaveGraph(tree,"slopedZANeg","run",condition);
    SaveGraph(tree,"offsetdZCNeg","run",condition);
    SaveGraph(tree,"slopedZCNeg","run",condition);
    SaveGraph(tree,"offsetdZAErrNeg","run",condition);
    SaveGraph(tree,"slopedZAErrNeg","run",condition);
    SaveGraph(tree,"offsetdZCErrNeg","run",condition);
    SaveGraph(tree,"slopedZCErrNeg","run",condition);
    SaveGraph(tree,"offsetdZAchi2Neg","run",condition);
    SaveGraph(tree,"slopedZAchi2Neg","run",condition);
    SaveGraph(tree,"offsetdZCchi2Neg","run",condition);
    SaveGraph(tree,"slopedZCchi2Neg","run",condition);    
    
    SaveGraph(tree,"offsetdZA","run",condition);
    SaveGraph(tree,"slopedZA","run",condition);
    SaveGraph(tree,"offsetdZC","run",condition);
    SaveGraph(tree,"slopedZC","run",condition);
    SaveGraph(tree,"offsetdZAErr","run",condition);
    SaveGraph(tree,"slopedZAErr","run",condition);
    SaveGraph(tree,"offsetdZCErr","run",condition);
    SaveGraph(tree,"slopedZCErr","run",condition);
    SaveGraph(tree,"offsetdZAchi2","run",condition);
    SaveGraph(tree,"slopedZAchi2","run",condition);
    SaveGraph(tree,"offsetdZCchi2","run",condition);
    SaveGraph(tree,"slopedZCchi2","run",condition);        

    SaveGraph(tree,"meanVertX","run",condition);
    SaveGraph(tree,"rmsVertX","run",condition);
    SaveGraph(tree,"meanVertY","run",condition);
    SaveGraph(tree,"rmsVertY","run",condition);
    SaveGraph(tree,"meanVertZ","run",condition);
    SaveGraph(tree,"rmsVertZ","run",condition);
    SaveGraph(tree,"vertStatus","run",condition);
    SaveGraph(tree,"meanMult","run",condition);
    SaveGraph(tree,"rmsMult","run",condition);
    SaveGraph(tree,"meanMultPos","run",condition);
    SaveGraph(tree,"rmsMultPos","run",condition);
    SaveGraph(tree,"errorMultPos","run",condition);
    SaveGraph(tree,"meanMultNeg","run",condition);
    SaveGraph(tree,"rmsMultNeg","run",condition);
    SaveGraph(tree,"errorMultNeg","run",condition);
    SaveGraph(tree,"vertAll","run",condition);
    SaveGraph(tree,"vertOK","run",condition);


    SaveGraph(tree,"meanPtAPos","run",condition);
    SaveGraph(tree,"mediumPtAPos","run",condition);
    SaveGraph(tree,"highPtAPos","run",condition);
    SaveGraph(tree,"meanPtCPos","run",condition);
    SaveGraph(tree,"mediumPtCPos","run",condition);
    SaveGraph(tree,"highPtCPos","run",condition);
    SaveGraph(tree,"meanPtANeg","run",condition);
    SaveGraph(tree,"mediumPtANeg","run",condition);
    SaveGraph(tree,"highPtANeg","run",condition);
    SaveGraph(tree,"meanPtCNeg","run",condition);
    SaveGraph(tree,"mediumPtCNeg","run",condition);
    SaveGraph(tree,"highPtCNeg","run",condition);
 
    SaveGraph(tree,"qOverPt","run",condition);
    SaveGraph(tree,"qOverPtA","run",condition);
    SaveGraph(tree,"qOverPtC","run",condition);

    SaveGraph(tree,"dcarAP0","run",condition);
    SaveGraph(tree,"dcarAP1","run",condition);
    SaveGraph(tree,"dcarCP0","run",condition);
    SaveGraph(tree,"dcarCP1","run",condition);

    condition = "";
    SaveGraph(tree,"tpcItsMatchA","run",condition);
    SaveGraph(tree,"tpcItsMatchHighPtA","run",condition);
    SaveGraph(tree,"tpcItsMatchC","run",condition);
    SaveGraph(tree,"tpcItsMatchHighPtC","run",condition);
    
    SaveGraph(tree,"phiPull","run",condition);
    SaveGraph(tree,"phiPullHighPt","run",condition);
    SaveGraph(tree,"ptPull","run",condition);
    SaveGraph(tree,"ptPullHighPt","run",condition);
    SaveGraph(tree,"yPull","run",condition);
    SaveGraph(tree,"yPullHighPt","run",condition);
    SaveGraph(tree,"zPull","run",condition);
    SaveGraph(tree,"zPullHighPt","run",condition);
    SaveGraph(tree,"lambdaPull","run",condition);
    SaveGraph(tree,"lambdaPullHighPt","run",condition);
    
    SaveGraph(tree,"tpcConstrainPhiA","run",condition);
    SaveGraph(tree,"tpcConstrainPhiC","run",condition);
     
    SaveGraph(tree,"deltaPt","run",condition);
    SaveGraph(tree,"deltaPtchi2","run",condition);
    SaveGraph(tree,"deltaPtA","run",condition);
    SaveGraph(tree,"deltaPtchi2A","run",condition);
    SaveGraph(tree,"deltaPtC","run",condition);
    SaveGraph(tree,"deltaPtchi2C","run",condition);
    SaveGraph(tree,"deltaPtA_Err","run",condition);
    SaveGraph(tree,"deltaPtA_Err","run",condition);
    SaveGraph(tree,"deltaPtC_Err","run",condition);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
 //save dca fit parameters
    SaveGraph(tree,"dcar_posA_0","run",condition);
    SaveGraph(tree,"dcar_posA_1","run",condition);
    SaveGraph(tree,"dcar_posA_2","run",condition);
    SaveGraph(tree,"dcar_posA_chi2","run",condition);
    SaveGraph(tree,"dcar_posA_0_Err","run",condition);
    SaveGraph(tree,"dcar_posA_1_Err","run",condition);
    SaveGraph(tree,"dcar_posA_2_Err","run",condition);

    SaveGraph(tree,"dcaz_posA_0","run",condition);
    SaveGraph(tree,"dcaz_posA_1","run",condition);
    SaveGraph(tree,"dcaz_posA_2","run",condition);
    SaveGraph(tree,"dcaz_posA_chi2","run",condition);
    SaveGraph(tree,"dcaz_posA_0_Err","run",condition);
    SaveGraph(tree,"dcaz_posA_1_Err","run",condition);
    SaveGraph(tree,"dcaz_posA_2_Err","run",condition);

    SaveGraph(tree,"dcaz_posC_0","run",condition);
    SaveGraph(tree,"dcaz_posC_1","run",condition);
    SaveGraph(tree,"dcaz_posC_2","run",condition);
    SaveGraph(tree,"dcaz_posC_chi2","run",condition);
    SaveGraph(tree,"dcaz_posC_0_Err","run",condition);
    SaveGraph(tree,"dcaz_posC_1_Err","run",condition);
    SaveGraph(tree,"dcaz_posC_2_Err","run",condition);

    SaveGraph(tree,"dcar_posC_0","run",condition);
    SaveGraph(tree,"dcar_posC_1","run",condition);
    SaveGraph(tree,"dcar_posC_2","run",condition);
    SaveGraph(tree,"dcar_posC_chi2","run",condition);
    SaveGraph(tree,"dcar_posC_0_Err","run",condition);
    SaveGraph(tree,"dcar_posC_1_Err","run",condition);
    SaveGraph(tree,"dcar_posC_2_Err","run",condition);

    SaveGraph(tree,"dcar_negA_0","run",condition);
    SaveGraph(tree,"dcar_negA_1","run",condition);
    SaveGraph(tree,"dcar_negA_2","run",condition);
    SaveGraph(tree,"dcar_negA_chi2","run",condition);
    SaveGraph(tree,"dcar_negA_0_Err","run",condition);
    SaveGraph(tree,"dcar_negA_1_Err","run",condition);
    SaveGraph(tree,"dcar_negA_2_Err","run",condition);

    SaveGraph(tree,"dcaz_negA_0","run",condition);
    SaveGraph(tree,"dcaz_negA_1","run",condition);
    SaveGraph(tree,"dcaz_negA_2","run",condition);
    SaveGraph(tree,"dcaz_negA_chi2","run",condition);
    SaveGraph(tree,"dcaz_negA_0_Err","run",condition);
    SaveGraph(tree,"dcaz_negA_1_Err","run",condition);
    SaveGraph(tree,"dcaz_negA_2_Err","run",condition);
    
    SaveGraph(tree,"dcaz_negC_0","run",condition);
    SaveGraph(tree,"dcaz_negC_1","run",condition);
    SaveGraph(tree,"dcaz_negC_2","run",condition);
    SaveGraph(tree,"dcaz_negC_chi2","run",condition);
    SaveGraph(tree,"dcaz_negC_0_Err","run",condition);
    SaveGraph(tree,"dcaz_negC_1_Err","run",condition);
    SaveGraph(tree,"dcaz_negC_2_Err","run",condition);
    
    SaveGraph(tree,"dcar_negC_0","run",condition);
    SaveGraph(tree,"dcar_negC_1","run",condition);
    SaveGraph(tree,"dcar_negC_2","run",condition);
    SaveGraph(tree,"dcar_negC_chi2","run",condition);
    SaveGraph(tree,"dcar_negC_0_Err","run",condition);
    SaveGraph(tree,"dcar_negC_1_Err","run",condition);
    SaveGraph(tree,"dcar_negC_2_Err","run",condition);

    SaveGraph(tree,"iroc_A_side","run",condition);    
    SaveGraph(tree,"iroc_C_side","run",condition);    
    SaveGraph(tree,"oroc_A_side","run",condition);    
    SaveGraph(tree,"oroc_C_side","run",condition);    

    //A/C side IROC                                                                                                                                                           
    SaveGraph(tree,"TPC_Occ_IROC.","run",condition);
    SaveGraph(tree,"TPC_Occ_OROC.","run",condition);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    tree->Write();
    out->Close();   
    if (tree) { delete tree; tree=0; }
    if (out) { delete out; out=0; }
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::SaveGraph(TTree* tree, const Char_t* y, const Char_t* x, const Char_t* condition)
{    
    //
    // Creates a Graph and writes the canvas to the current directory
    // called by ProduceTrends function.
    //
    
    TString s(y);
    s += ':';
    s += x;
    tree->Draw(s.Data(),condition,"goff");
    TCanvas* c = new TCanvas(s.Data(),s.Data());
    c->cd();
    TPad* p = new TPad("pad0","pad0",0,0,1,1);
    p->Draw();
    p->cd();
    if (tree->GetSelectedRows() > 0) {
      TGraph* graph = new TGraph(tree->GetSelectedRows(), tree->GetV2(), tree->GetV1());
      graph->Draw("AP");
      graph->GetXaxis()->SetTitle(x);
      graph->GetYaxis()->SetTitle(y);
      c->Write(s.Data());
      delete graph;
    }
    //graph->Write(s.Data());
    delete c;
    
return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDCARPhi(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA R imperfections
    //
    
    if (!pcstream) return 8;
    if (!pTPC) return 8;
        
    // variables:
    static Double_t offsetdRA=0;
    static Double_t slopedRA=0;
    static Double_t offsetdRC=0;
    static Double_t slopedRC=0;
    static Double_t offsetdRAErr=0;
    static Double_t slopedRAErr=0;
    static Double_t offsetdRCErr=0;
    static Double_t slopedRCErr=0;
    static Double_t offsetdRAchi2=0;
    static Double_t slopedRAchi2=0;
    static Double_t offsetdRCchi2=0;
    static Double_t slopedRCchi2=0;
    static Double_t dcarAP0 = 0;
    static Double_t dcarAP1 = 0;
    static Double_t dcarCP0 = 0;
    static Double_t dcarCP1 = 0;

    //AliPerformanceTPC* pTPC =  dynamic_cast<AliPerformanceTPC*>(pTPCObject);    
    
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_3_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_3_5_7"));
        if(!his3D) return 8;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
    }
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    if (his3D && !fgForceTHnSparse) { 
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        his2D = pTPC->GetTPCTrackHisto()->Projection(3,5);
    }            
    if(!his2D) return 8;

    
    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdRC=fpol1->GetParameter(0);
    slopedRC=fpol1->GetParameter(1);
    offsetdRCchi2=fpol1->GetChisquare();
    slopedRCchi2=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdRA=fpol1->GetParameter(0);
    slopedRA=fpol1->GetParameter(1);
    offsetdRAErr=fpol1->GetParError(0);
    slopedRAErr=fpol1->GetParError(1);
    offsetdRAchi2=fpol1->GetChisquare();
    slopedRAchi2=fpol1->GetChisquare();
    //
    printf("DCA R QA report\n");
    printf("offsetdRA\t%f\n",offsetdRA);
    printf("slopedRA\t%f\n",slopedRA);
    printf("offsetdRC\t%f\n",offsetdRC);
    printf("slopedRC\t%f\n",slopedRC);

    //
    //extraction of DCAr versus pt
    //
    TLinearFitter linearFit;
    linearFit.SetFormula("pol1");
    TObjArray arrayWidth;  
    TH1 *width;
    Int_t nXbins;
    Double_t x,y;
    Int_t pn = 1;

    if(!his3D)
      return 8;
    his3D->GetYaxis()->SetRangeUser(-1,1);

    //get his2D in A Side
    his3D->GetYaxis()->SetRangeUser(0,1);
    his3D->GetZaxis()->SetRangeUser(0.35,8);
    his2D  = dynamic_cast<TH2*>(his3D->Project3D("xz"));
    if (his2D)
      his2D->FitSlicesY(0,0,-1,0,"QNR",&arrayWidth);

    width =  dynamic_cast<TH1*>(arrayWidth.At(2));
    if (width) 
    {
      nXbins = width->GetNbinsX();
    for(Int_t i=2; i<nXbins; i++){
      x = width->GetBinCenter(i);
      if(x!=0)
	x = 1.0/(x*x);
      y = width->GetBinContent(i);
      y = y*y;
      linearFit.AddPoint(&x,y,1);
    }
    if(!linearFit.Eval()){
      
      dcarAP0 = linearFit.GetParameter(0);
      if(dcarAP0!=0)
	pn = Int_t(TMath::Abs(dcarAP0)/dcarAP0);
      dcarAP0 = pn*TMath::Sqrt(TMath::Abs(dcarAP0));

      dcarAP1 = linearFit.GetParameter(1);
      if(dcarAP1!=0)
	pn = Int_t(TMath::Abs(dcarAP1)/dcarAP1);
      dcarAP1 = pn*TMath::Sqrt(TMath::Abs(dcarAP1));
    }
    }

    linearFit.ClearPoints();
    
    //get his2D in C Side
    his3D->GetYaxis()->SetRangeUser(-1,-0.001);
    his2D  = dynamic_cast<TH2*>(his3D->Project3D("xz"));
    if (his2D)
      his2D->FitSlicesY(0,0,-1,0,"QNR",&arrayWidth);
      width =  dynamic_cast<TH1*>(arrayWidth.At(2));
    if (width) 
    {
      nXbins = width->GetNbinsX();
    for(Int_t i=2; i<nXbins; i++){
      x = width->GetBinCenter(i);
      if(x!=0)
	x = 1.0/(x*x);
      y = width->GetBinContent(i);
      y = y*y;
      linearFit.AddPoint(&x,y);
    }
    if(!linearFit.Eval()){
      dcarCP0 = linearFit.GetParameter(0);
      if(dcarCP0!=0)
	pn = Int_t(TMath::Abs(dcarCP0)/dcarCP0);
      dcarCP0 = pn*TMath::Sqrt(TMath::Abs(dcarCP0));

      dcarCP1 = linearFit.GetParameter(1);
      if(dcarCP1!=0)
	pn = Int_t(TMath::Abs(dcarCP1)/dcarCP1);
      dcarCP1 = pn*TMath::Sqrt(TMath::Abs(dcarCP1));
    }
    }
    his3D->GetYaxis()->SetRangeUser(-1,1);
    his3D->GetZaxis()->SetRangeUser(0,20);

    //
    // dump values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdRA="<< offsetdRA<<
        "slopedRA="<< slopedRA<<
        "offsetdRC="<< offsetdRC<<
        "slopedRC="<<slopedRC<<
        //
        "offsetdRAErr="<< offsetdRAErr<<
        "slopedRAErr="<< slopedRAErr<<
        "offsetdRCErr="<< offsetdRCErr<<
        "slopedRCErr="<<slopedRCErr<<
        //
        "offsetdRAchi2="<< offsetdRAchi2<<
        "slopedRAchi2="<< slopedRAchi2<<
        "offsetdRCchi2="<< offsetdRCchi2<<
        "slopedRCchi2="<<slopedRCchi2<<
        //
        "dcarAP0="<<dcarAP0<<
        "dcarAP1="<<dcarAP1<<
        "dcarCP0="<<dcarCP0<<
        "dcarCP1="<<dcarCP1;
        
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDCARPhiPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA R imperfections for positive particles
    //
    
    if (!pcstream) return 16;
    if (!pTPC) return 16;

    // variables:
    static Double_t offsetdRAPos=0;
    static Double_t slopedRAPos=0;
    static Double_t offsetdRCPos=0;
    static Double_t slopedRCPos=0;
    static Double_t offsetdRAErrPos=0;
    static Double_t slopedRAErrPos=0;
    static Double_t offsetdRCErrPos=0;
    static Double_t slopedRCErrPos=0;
    static Double_t offsetdRAchi2Pos=0;
    static Double_t slopedRAchi2Pos=0;
    static Double_t offsetdRCchi2Pos=0;
    static Double_t slopedRCchi2Pos=0;

    //AliPerformanceTPC* pTPC =  dynamic_cast<AliPerformanceTPC*>(pTPCObject);    
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_7"));
        if(!his3D) return 16;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
    }
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    if (his3D && !fgForceTHnSparse) { 
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        pTPC->GetTPCTrackHisto()->GetAxis(8)->SetRangeUser(0,1.5);        
        his2D = pTPC->GetTPCTrackHisto()->Projection(3,5);
        pTPC->GetTPCTrackHisto()->GetAxis(8)->SetRangeUser(-1.5,1.5);
    }            
    if(!his2D) return 16;
    
    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdRCPos=fpol1->GetParameter(0);
    slopedRCPos=fpol1->GetParameter(1);
    offsetdRCchi2Pos=fpol1->GetChisquare();
    slopedRCchi2Pos=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdRAPos=fpol1->GetParameter(0);
    slopedRAPos=fpol1->GetParameter(1);
    offsetdRAErrPos=fpol1->GetParError(0);
    slopedRAErrPos=fpol1->GetParError(1);
    offsetdRAchi2Pos=fpol1->GetChisquare();
    slopedRAchi2Pos=fpol1->GetChisquare();
    //
    printf("DCA R QA Pos report\n");
    printf("offsetdRAPos\t%f\n",offsetdRAPos);
    printf("slopedRAPos\t%f\n",slopedRAPos);
    printf("offsetdRCPos\t%f\n",offsetdRCPos);
    printf("slopedRCPos\t%f\n",slopedRCPos);
    //
    // dump values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdRAPos="<< offsetdRAPos<<
        "slopedRAPos="<< slopedRAPos<<
        "offsetdRCPos="<< offsetdRCPos<<
        "slopedRCPos="<<slopedRCPos<<
        //
        "offsetdRAErrPos="<< offsetdRAErrPos<<
        "slopedRAErrPos="<< slopedRAErrPos<<
        "offsetdRCErrPos="<< offsetdRCErrPos<<
        "slopedRCErrPos="<<slopedRCErrPos<<
        //
        "offsetdRAchi2Pos="<< offsetdRAchi2Pos<<
        "slopedRAchi2Pos="<< slopedRAchi2Pos<<
        "offsetdRCchi2Pos="<< offsetdRCchi2Pos<<
        "slopedRCchi2Pos="<<slopedRCchi2Pos;
        
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDCARPhiNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA R imperfections for negative particles
    //
    if (!pcstream) return 32;
    if (!pTPC) return 32;

    // variables:
    static Double_t offsetdRANeg=0;
    static Double_t slopedRANeg=0;
    static Double_t offsetdRCNeg=0;
    static Double_t slopedRCNeg=0;
    static Double_t offsetdRAErrNeg=0;
    static Double_t slopedRAErrNeg=0;
    static Double_t offsetdRCErrNeg=0;
    static Double_t slopedRCErrNeg=0;
    static Double_t offsetdRAchi2Neg=0;
    static Double_t slopedRAchi2Neg=0;
    static Double_t offsetdRCchi2Neg=0;
    static Double_t slopedRCchi2Neg=0;

    //AliPerformanceTPC* pTPC =  dynamic_cast<AliPerformanceTPC*>(pTPCObject);    
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_7"));
	if(!his3D) return 32;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
    }
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    if (his3D && !fgForceTHnSparse) {
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        pTPC->GetTPCTrackHisto()->GetAxis(8)->SetRangeUser(-1.5,0);        
        his2D = pTPC->GetTPCTrackHisto()->Projection(3,5);
        pTPC->GetTPCTrackHisto()->GetAxis(8)->SetRangeUser(-1.5,1.5);
    }            
    if(!his2D) return 32;

    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdRCNeg=fpol1->GetParameter(0);
    slopedRCNeg=fpol1->GetParameter(1);
    offsetdRCchi2Neg=fpol1->GetChisquare();
    slopedRCchi2Neg=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdRANeg=fpol1->GetParameter(0);
    slopedRANeg=fpol1->GetParameter(1);
    offsetdRAErrNeg=fpol1->GetParError(0);
    slopedRAErrNeg=fpol1->GetParError(1);
    offsetdRAchi2Neg=fpol1->GetChisquare();
    slopedRAchi2Neg=fpol1->GetChisquare();
    //
    printf("DCA R QA Neg report\n");
    printf("offsetdRANeg\t%f\n",offsetdRANeg);
    printf("slopedRANeg\t%f\n",slopedRANeg);
    printf("offsetdRCNeg\t%f\n",offsetdRCNeg);
    printf("slopedRCNeg\t%f\n",slopedRCNeg);
    //
    // dump drift QA values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdRANeg="<< offsetdRANeg<<
        "slopedRANeg="<< slopedRANeg<<
        "offsetdRCNeg="<< offsetdRCNeg<<
        "slopedRCNeg="<<slopedRCNeg<<
        //
        "offsetdRAErrNeg="<< offsetdRAErrNeg<<
        "slopedRAErrNeg="<< slopedRAErrNeg<<
        "offsetdRCErrNeg="<< offsetdRCErrNeg<<
        "slopedRCErrNeg="<<slopedRCErrNeg<<
        //
        "offsetdRAchi2Neg="<< offsetdRAchi2Neg<<
        "slopedRAchi2Neg="<< slopedRAchi2Neg<<
        "offsetdRCchi2Neg="<< offsetdRCchi2Neg<<
        "slopedRCchi2Neg="<<slopedRCchi2Neg;
        
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeNCL(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
  //
  // Analyse number of TPC clusters 
  // Standard NCl Ncl/Findable histograms 
  //   
  // Functionality for missing chambers detection.
  //    <NCL>med per phi|sector|TPC
  //    <Ntr>med per phi|sector|TPC

    if (!pcstream) return 1;
    if (!pTPC) return 1;
 
    // variables:
    static Double_t meanTPCnclF=0;
    static Double_t rmsTPCnclF=0;
    static Double_t meanTPCChi2=0;
    static Double_t rmsTPCChi2=0;  
    static Double_t slopeATPCnclF=0;
    static Double_t slopeCTPCnclF=0;
    static Double_t slopeATPCnclFErr=0;
    static Double_t slopeCTPCnclFErr=0;
    static Double_t meanTPCncl=0;
    static Double_t rmsTPCncl=0;
    static Double_t slopeATPCncl=0;
    static Double_t slopeCTPCncl=0;
    static Double_t slopeATPCnclErr=0;
    static Double_t slopeCTPCnclErr=0;  
    //
    TH1* his1D=0;
    TH3* his3D0=0;
    TH3* his3D1=0;
    TH3* his3D2=0;
    TProfile* hprof=0;
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    //
    // all clusters
    // only events with rec. vertex
    // eta cut - +-1
    // pt cut  - 0.250 GeV
    if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-1.,1.);
    if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0.25,10);
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_0_5_7")) {    
        his3D0 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_0_5_7"));
        if(!his3D0) return 1;
        his3D0->GetYaxis()->SetRangeUser(-1,1);
        his3D0->GetZaxis()->SetRangeUser(0.25,10);
    }
    if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_1_5_7")) {    
        his3D1 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_1_5_7"));
        if(!his3D1) return 1;
        his3D1->GetYaxis()->SetRangeUser(-1,1);
        his3D1->GetZaxis()->SetRangeUser(0.25,10);
    }
    if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_2_5_7")) {    
        his3D2 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_2_5_7"));
        if(!his3D2) return 1;
        his3D2->GetYaxis()->SetRangeUser(-1,1);
        his3D2->GetZaxis()->SetRangeUser(0.25,10);
        his3D2->GetXaxis()->SetRangeUser(0.4,1.1);        
    }    
    

    if (his3D0 && !fgForceTHnSparse) { 
         his1D = his3D0->Project3D("x"); 
    } else {
         his1D = pTPC->GetTPCTrackHisto()->Projection(0);
    }
 
    meanTPCncl= his1D->GetMean();
    rmsTPCncl= his1D->GetRMS();
    delete his1D;
    
    if (his3D1 && !fgForceTHnSparse) {
         his1D = his3D1->Project3D("x"); 
    } else {
         his1D = pTPC->GetTPCTrackHisto()->Projection(1);
    }
          
    meanTPCChi2= his1D->GetMean();
    rmsTPCChi2= his1D->GetRMS();
    delete his1D;  
    
   if (his3D0 && !fgForceTHnSparse) {
        hprof = (dynamic_cast<TH2*>(his3D0->Project3D("xy")))->ProfileX(); 
    } else {
        hprof = pTPC->GetTPCTrackHisto()->Projection(0,5)->ProfileX();
    }
    if(!hprof) return 1;
    
    hprof->Fit(fpol1,"QNR","QNR",0.1,0.8);
    slopeATPCncl= fpol1->GetParameter(1);
    slopeATPCnclErr= fpol1->GetParError(1);
    hprof->Fit(fpol1,"QNR","QNR",-0.8,-0.1);
    slopeCTPCncl= fpol1->GetParameter(1);
    slopeCTPCnclErr= fpol1->GetParameter(1);
    delete hprof;
    
    //
    // findable clusters
    //
    
   if (his3D2 && !fgForceTHnSparse) {
        his1D = his3D2->Project3D("x"); 
    } else {    
        pTPC->GetTPCTrackHisto()->GetAxis(2)->SetRangeUser(0.4,1.1);
        his1D = pTPC->GetTPCTrackHisto()->Projection(2);
    }    
        
    meanTPCnclF= his1D->GetMean();
    rmsTPCnclF= his1D->GetRMS();
    delete his1D;
    
   if (his3D2 && !fgForceTHnSparse) { 
         his1D = (dynamic_cast<TH2*>(his3D2->Project3D("xy")))->ProfileX(); 
    } else {    
        pTPC->GetTPCTrackHisto()->GetAxis(2)->SetRangeUser(0.4,1.1);
        his1D = pTPC->GetTPCTrackHisto()->Projection(2,5)->ProfileX();
    }      
    if(!his1D) return 1;
    
    his1D->Fit(fpol1,"QNR","QNR",0.1,0.8);
    slopeATPCnclF= fpol1->GetParameter(1);
    slopeATPCnclFErr= fpol1->GetParError(1);
    his1D->Fit(fpol1,"QNR","QNR",-0.8,-0.1);
    slopeCTPCnclF= fpol1->GetParameter(1);
    slopeCTPCnclFErr= fpol1->GetParameter(1);
    delete his1D;
        
    if(fgForceTHnSparse) pTPC->GetTPCTrackHisto()->GetAxis(2)->SetRangeUser(0,10);
    
    printf("Cluster QA report\n");
    printf("meanTPCnclF=\t%f\n",meanTPCnclF);
    printf("rmsTPCnclF=\t%f\n",rmsTPCnclF);
    printf("slopeATPCnclF=\t%f\n",slopeATPCnclF);
    printf("slopeCTPCnclF=\t%f\n",slopeCTPCnclF);
    printf("meanTPCncl=\t%f\n",meanTPCncl);
    printf("rmsTPCncl=\t%f\n",rmsTPCncl);
    printf("slopeATPCncl=\t%f\n",slopeATPCncl);
    printf("slopeCTPCncl=\t%f\n",slopeCTPCncl);
    printf("meanTPCChi2=\t%f\n",meanTPCChi2);
    printf("rmsTPCChi2=\t%f\n",rmsTPCChi2);
    //
    // dump results to the tree
    //
    (*pcstream)<<"tpcQA"<<
      "meanTPCnclF="<<meanTPCnclF <<   
      "rmsTPCnclF="<<rmsTPCnclF <<
      "meanTPCChi2="<<meanTPCChi2 <<
      "rmsTPCChi2="<<rmsTPCChi2 <<
      "slopeATPCnclF="<< slopeATPCnclF<<
      "slopeCTPCnclF="<< slopeCTPCnclF<<
      "slopeATPCnclFErr="<< slopeATPCnclFErr<<
      "slopeCTPCnclFErr="<< slopeCTPCnclFErr<<
      "meanTPCncl="<<meanTPCncl <<
      "rmsTPCncl="<< rmsTPCncl<<
      "slopeATPCncl="<< slopeATPCncl<<
      "slopeCTPCncl="<< slopeCTPCncl<<
      "slopeATPCnclErr="<< slopeATPCnclErr<<
      "slopeCTPCnclErr="<< slopeCTPCnclErr;
    //
    // Get ncl:sector map estimator without having ncl map
    Double_t bz = 0;
    if (TGeoGlobalMagField::Instance()->GetField()) {
      Int_t polarityL3 = AliTPCcalibDB::GetL3Polarity(pTPC->GetRunNumber());
      bz = AliTPCcalibDB::GetBz(pTPC->GetRunNumber());
      if (polarityL3>0) bz*=-1; 
    }
    
    TH3D * hisNclpos = dynamic_cast<TH3D*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_2_5_6"));
    if(!hisNclpos) return 1;
    TH3D * hisNclneg = dynamic_cast<TH3D*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_2_5_6"));
    if(!hisNclneg) return 1;
    Int_t nbins= hisNclpos->GetZaxis()->GetNbins();
    TAxis* paxisPhi=hisNclpos->GetZaxis();
    TH3D * hisNcl= 0;
    //
    // 0.) As an estimator most probable value of the ncl/nclfindable used
    //
    static TGraphErrors * graphNclMostProbPhi[8]={0};
    static TGraphErrors * graphNclMostProbPhiSector[8]={0};
    TVectorD vecMP(nbins);
    TVectorD vecPhi(nbins);
    TVectorD vecMostProb(nbins);
    TVectorD vecSector(nbins);
    TVectorD vecNclSector(nbins);
    TVectorD vecNcl(nbins);
    //
    //
    for (Int_t jgr=0; jgr<8; jgr++){
      Int_t igr=jgr%4;
      hisNcl= ((igr%2==0)) ?  new TH3D(*hisNclpos) : new TH3D(*hisNclneg);   // track sign
      if (igr<2) {
	hisNcl->GetYaxis()->SetRangeUser(0.3,0.9);   // A side
      }else{
	hisNcl->GetYaxis()->SetRangeUser(-0.9,-0.3); // C side
      }
      Int_t sign=((igr%2==0)) ? 1:-1;
      if (bz>0) sign*=-1;  // sign of the Bz
      TH2* his2=  (TH2*)hisNcl->Project3D("xz");
      for (Int_t ibin=1; ibin<=nbins; ibin++){
	TH1 * his1D = (TH1*)his2->ProjectionY("his1D",ibin,ibin);
	if (jgr<4) vecMP[ibin-1]=his1D->GetBinCenter(his1D->GetMaximumBin());
	if (jgr>=4) vecMP[ibin-1]=his1D->GetEntries();
	vecPhi[ibin-1]=9.*paxisPhi->GetBinCenter((ibin+sign*4+nbins)%nbins)/TMath::Pi();
	delete his1D;
      }
      graphNclMostProbPhi[jgr]=new TGraphErrors(nbins,vecPhi.GetMatrixArray(), vecMP.GetMatrixArray());          
      graphNclMostProbPhi[jgr]->SetMarkerStyle(21+igr);
      graphNclMostProbPhi[jgr]->SetMarkerColor(1+igr);
      graphNclMostProbPhi[jgr]->SetLineColor(1+igr);
      delete his2;
      delete hisNcl;
      hisNcl=0;
    }    
    //
    //
    for (Int_t igr=0;igr<8; igr++){
      for (Int_t isec=0; isec<18; isec++){
	Int_t bins=0;
	for (Int_t ibin=0;ibin<nbins; ibin++){    
	  if (TMath::Abs(graphNclMostProbPhi[igr]->GetX()[ibin]-isec-0.5)>0.4) continue;
	  vecNcl[bins]=graphNclMostProbPhi[igr]->GetY()[ibin];
	  bins++;
	}
	vecSector[isec]=isec;
	vecNclSector[isec]=TMath::Median(bins, vecNcl.GetMatrixArray());
      }
      graphNclMostProbPhiSector[igr]=new TGraphErrors(18, vecSector.GetMatrixArray(), vecNclSector.GetMatrixArray());      
    }
    static TVectorD normMedian(8);
    for (Int_t jgr=0;jgr<8; jgr++){
      Int_t igr=jgr%4;
      normMedian[igr] = TMath::Median(nbins, graphNclMostProbPhi[igr]->GetY());
      graphNclMostProbPhi[igr]->SetMarkerStyle(21+igr);
      graphNclMostProbPhi[igr]->SetMarkerColor(1+igr);
      graphNclMostProbPhiSector[igr]->SetMarkerStyle(21+igr);
      graphNclMostProbPhiSector[igr]->SetMarkerColor(1+igr);
    } 

    (*pcstream)<<"tpcQA"<<
      "grNclPhiMedian.="<<&normMedian<<            //  median value (144 phi bins)  of the number of clusters 
      "grNclPhiPosA.="<< graphNclMostProbPhi[0]<<  //  phi NCL/findable profile per phi bin - positive tracks A side
      "grNclPhiNegA.="<< graphNclMostProbPhi[1]<<  //  phi NCL/findable profile per phi bin - negative tracks A side
      "grNclPhiPosC.="<< graphNclMostProbPhi[2]<<  //  phi NCL/findable profile per phi bin - positive tracks C side
      "grNclPhiNegC.="<< graphNclMostProbPhi[3]<<  //  phi NCL/findable profile per phi bin - negative tracks C side
      "grNtrPhiPosA.="<< graphNclMostProbPhi[4]<<  //  phi entries per phi bin - positive tracks A side
      "grNtrPhiNegA.="<< graphNclMostProbPhi[5]<<  //  phi entries per phi bin - negative tracks A side
      "grNtrPhiPosC.="<< graphNclMostProbPhi[6]<<  //  phi entries per phi bin - positive tracks C side
      "grNtrPhiNegC.="<< graphNclMostProbPhi[7];   //  phi entries per phi bin - negative tracks C side

    (*pcstream)<<"tpcQA"<<
      "grNclSectorPosA.="<< graphNclMostProbPhiSector[0]<<  //  sector NCL/findable profile
      "grNclSectorNegA.="<< graphNclMostProbPhiSector[1]<<  // 
      "grNclSectorPosC.="<< graphNclMostProbPhiSector[2]<<  // 
      "grNclSectorNegC.="<< graphNclMostProbPhiSector[3]<<  // 
      "grNtrSectorPosA.="<< graphNclMostProbPhiSector[4]<<  //  sector entries per sector
      "grNtrSectorNegA.="<< graphNclMostProbPhiSector[5]<<  // 
      "grNtrSectorPosC.="<< graphNclMostProbPhiSector[6]<<  // 
      "grNtrSectorNegC.="<< graphNclMostProbPhiSector[7];   // 
 

    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDrift(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA Z imperferctions (drift velocity)
    //
    
    if (!pcstream) return 2;
    if (!pTPC) return 2;

    // variables:
    static Double_t offsetdZA=0;
    static Double_t slopedZA=0;
    static Double_t offsetdZC=0;
    static Double_t slopedZC=0;
    static Double_t offsetdZAErr=0;
    static Double_t slopedZAErr=0;
    static Double_t offsetdZCErr=0;
    static Double_t slopedZCErr=0;
    static Double_t offsetdZAchi2=0;
    static Double_t slopedZAchi2=0;
    static Double_t offsetdZCchi2=0;
    static Double_t slopedZCchi2=0;
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
   if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_4_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_4_5_7"));
        if(!his3D) return 2;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
   }

   if (his3D && !fgForceTHnSparse) { 
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        his2D = pTPC->GetTPCTrackHisto()->Projection(4,5);
    }        
    if(!his2D) return 2;
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdZC=fpol1->GetParameter(0);
    slopedZC=fpol1->GetParameter(1);
    offsetdZCErr=fpol1->GetParError(0);
    slopedZCErr=fpol1->GetParError(1);        
    offsetdZCchi2=fpol1->GetChisquare();
    slopedZCchi2=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdZA=fpol1->GetParameter(0);
    slopedZA=fpol1->GetParameter(1);
    offsetdZAErr=fpol1->GetParError(0);
    slopedZAErr=fpol1->GetParError(1);
    offsetdZAchi2=fpol1->GetChisquare();
    slopedZAchi2=fpol1->GetChisquare();
    //
    printf("Drift velocity QA report\n");
    printf("offsetdZA\t%f\n",offsetdZA);
    printf("slopedZA\t%f\n",slopedZA);
    printf("offsetdZC\t%f\n",offsetdZC);
    printf("slopedZC\t%f\n",slopedZC);
    //
    // dump drift QA values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdZA="<< offsetdZA<<
        "slopedZA="<< slopedZA<<
        "offsetdZC="<< offsetdZC<<
        "slopedZC="<<slopedZC<<
        //
        "offsetdZAErr="<< offsetdZAErr<<
        "slopedZAErr="<< slopedZAErr<<
        "offsetdZCErr="<< offsetdZCErr<<
        "slopedZCErr="<<slopedZCErr<<
        //
        "offsetdZAchi2="<< offsetdZAchi2<<
        "slopedZAchi2="<< slopedZAchi2<<
        "offsetdZCchi2="<< offsetdZCchi2<<
        "slopedZCchi2="<<slopedZCchi2;
    
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDriftPos(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA Z imperferctions (drift velocity)
    // for positive particles
    //
    if (!pcstream) return 64;
    if (!pTPC) return 64;
        
    // variables:
    static Double_t offsetdZAPos=0;
    static Double_t slopedZAPos=0;
    static Double_t offsetdZCPos=0;
    static Double_t slopedZCPos=0;
    static Double_t offsetdZAErrPos=0;
    static Double_t slopedZAErrPos=0;
    static Double_t offsetdZCErrPos=0;
    static Double_t slopedZCErrPos=0;
    static Double_t offsetdZAchi2Pos=0;
    static Double_t slopedZAchi2Pos=0;
    static Double_t offsetdZCchi2Pos=0;
    static Double_t slopedZCchi2Pos=0;
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
   if (pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_4_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_4_5_7"));
        if(!his3D) return 64;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
    } 

    if (his3D && !fgForceTHnSparse) { 
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        his2D = pTPC->GetTPCTrackHisto()->Projection(4,5);
    }            
    if(!his2D) return 64;
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;
    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdZCPos=fpol1->GetParameter(0);
    slopedZCPos=fpol1->GetParameter(1);
    offsetdZCErrPos=fpol1->GetParError(0);
    slopedZCErrPos=fpol1->GetParError(1);        
    offsetdZCchi2Pos=fpol1->GetChisquare();
    slopedZCchi2Pos=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdZAPos=fpol1->GetParameter(0);
    slopedZAPos=fpol1->GetParameter(1);
    offsetdZAErrPos=fpol1->GetParError(0);
    slopedZAErrPos=fpol1->GetParError(1);
    offsetdZAchi2Pos=fpol1->GetChisquare();
    slopedZAchi2Pos=fpol1->GetChisquare();
    //
    printf("Drift velocity QA report\n");
    printf("offsetdZAPos\t%f\n",offsetdZAPos);
    printf("slopedZAPos\t%f\n",slopedZAPos);
    printf("offsetdZCPos\t%f\n",offsetdZCPos);
    printf("slopedZCPos\t%f\n",slopedZCPos);
    //
    // dump drift QA values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdZAPos="<< offsetdZAPos<<
        "slopedZAPos="<< slopedZAPos<<
        "offsetdZCPos="<< offsetdZCPos<<
        "slopedZCPos="<<slopedZCPos<<
        //
        "offsetdZAErrPos="<< offsetdZAErrPos<<
        "slopedZAErrPos="<< slopedZAErrPos<<
        "offsetdZCErrPos="<< offsetdZCErrPos<<
        "slopedZCErrPos="<<slopedZCErrPos<<
        //
        "offsetdZAchi2Pos="<< offsetdZAchi2Pos<<
        "slopedZAchi2Pos="<< slopedZAchi2Pos<<
        "offsetdZCchi2Pos="<< offsetdZCchi2Pos<<
        "slopedZCchi2Pos="<<slopedZCchi2Pos;
        
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeDriftNeg(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA Z imperferctions (drift velocity)
    // for negative particles
    //
    if (!pcstream) return 128;
    if (!pTPC) return 128;
            
    // variables:
    static Double_t offsetdZANeg=0;
    static Double_t slopedZANeg=0;
    static Double_t offsetdZCNeg=0;
    static Double_t slopedZCNeg=0;
    static Double_t offsetdZAErrNeg=0;
    static Double_t slopedZAErrNeg=0;
    static Double_t offsetdZCErrNeg=0;
    static Double_t slopedZCErrNeg=0;
    static Double_t offsetdZAchi2Neg=0;
    static Double_t slopedZAchi2Neg=0;
    static Double_t offsetdZCchi2Neg=0;
    static Double_t slopedZCchi2Neg=0;
    TH1* his1D=0;
    TH2* his2D=0;
    TH3* his3D=0;
    
    
   if (pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_4_5_7")) {    
        his3D = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_4_5_7"));
        if(!his3D) return 128;
        his3D->GetYaxis()->SetRangeUser(-1,1);
        his3D->GetZaxis()->SetRangeUser(0.25,10);
    }
    if (his3D && !fgForceTHnSparse) { 
        his2D = dynamic_cast<TH2*>(his3D->Project3D("xy")); 
    } else {    
        his2D = pTPC->GetTPCTrackHisto()->Projection(4,5);
    }                
    if(!his2D) return 128;
    
    static TF1 *fpol1 = new TF1("fpol1","pol1");
    TObjArray arrayFit;
    his2D->FitSlicesY(0,0,-1,10,"QNR",&arrayFit);
    delete his2D;
    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",-0.8,-0.1);
    offsetdZCNeg=fpol1->GetParameter(0);
    slopedZCNeg=fpol1->GetParameter(1);
    offsetdZCErrNeg=fpol1->GetParError(0);
    slopedZCErrNeg=fpol1->GetParError(1);        
    offsetdZCchi2Neg=fpol1->GetChisquare();
    slopedZCchi2Neg=fpol1->GetChisquare();
    //
    his1D->Fit(fpol1,"QNRROB=0.8","QNR",0.1,0.8);
    offsetdZANeg=fpol1->GetParameter(0);
    slopedZANeg=fpol1->GetParameter(1);
    offsetdZAErrNeg=fpol1->GetParError(0);
    slopedZAErrNeg=fpol1->GetParError(1);
    offsetdZAchi2Neg=fpol1->GetChisquare();
    slopedZAchi2Neg=fpol1->GetChisquare();
    //
    printf("Drift velocity QA report\n");
    printf("offsetdZANeg\t%f\n",offsetdZANeg);
    printf("slopedZANeg\t%f\n",slopedZANeg);
    printf("offsetdZCNeg\t%f\n",offsetdZCNeg);
    printf("slopedZCNeg\t%f\n",slopedZCNeg);
    //
    // dump drift QA values
    //
    (*pcstream)<<"tpcQA"<<
        "offsetdZANeg="<< offsetdZANeg<<
        "slopedZANeg="<< slopedZANeg<<
        "offsetdZCNeg="<< offsetdZCNeg<<
        "slopedZCNeg="<<slopedZCNeg<<
        //
        "offsetdZAErrNeg="<< offsetdZAErrNeg<<
        "slopedZAErrNeg="<< slopedZAErrNeg<<
        "offsetdZCErrNeg="<< offsetdZCErrNeg<<
        "slopedZCErrNeg="<<slopedZCErrNeg<<
        //
        "offsetdZAchi2Neg="<< offsetdZAchi2Neg<<
        "slopedZAchi2Neg="<< slopedZAchi2Neg<<
        "offsetdZCchi2Neg="<< offsetdZCchi2Neg<<
        "slopedZCchi2Neg="<<slopedZCchi2Neg;
    
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeGain(const AliPerformanceDEdx* pTPCgain, TTreeSRedirector* const pcstream)
{
    //
    // Analyse Gain
    //
    
    if (!pcstream) return 4;
    if (!pTPCgain) return 4;

    static TVectorD meanMIPvsSector(36);
    static TVectorD sigmaRelMIPvsSector(36);
    static TVectorD sector(36);
    static Float_t meanMIP = 0;
    static Float_t resolutionMIP = 0;
    static Float_t attachSlopeC = 0;
    static Float_t attachSlopeA = 0;
    static Float_t meanMIPele = 0;
    static Float_t resolutionMIPele = 0;
    static Float_t electroMIPSeparation = 0;
    
    TH1 * his1D = 0;
    //TH1 * hisProj1D=0;
    TH2* his2D=0;
     
    meanMIPvsSector.Zero();
    sigmaRelMIPvsSector.Zero();
    //
    // select MIP particles
    //
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(7)->SetRangeUser(0.4,0.55);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(0)->SetRangeUser(35,60);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(6)->SetRangeUser(80,160);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-1,1);
    //
    // MIP position and resolution
    //    
    TF1 gausFit("gausFit","gaus");
   
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_0") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_0")->Clone());
    } else {
       his1D =  pTPCgain->GetDeDxHisto()->Projection(0);
    }
    if(!his1D) return 4;
    his1D->Fit(&gausFit,"QN","QN");

    meanMIP = gausFit.GetParameter(1);
    resolutionMIP = 0;
    if (meanMIP!=0) resolutionMIP = gausFit.GetParameter(2)/meanMIP;
    //removedtotest// delete his1D;
    
    //
    // MIP position vs. dip angle (attachment)
    //    
      if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-3,0); // C side
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_c_0_5") && !fgForceTHnSparse) {    
        his2D = dynamic_cast<TH2*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_c_0_5")->Clone());
    } else {
        his2D =  pTPCgain->GetDeDxHisto()->Projection(0,5);
    }        
    if(!his2D) return 4;

    TF1 * fpol = new TF1("fpol","pol1");
    TObjArray arrayFit;
    his2D->FitSlicesY(0,0,-1,10,"QN",&arrayFit);    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpol,"QNROB=0.8","QNR",-1,0);
    attachSlopeC = fpol->GetParameter(1);
    //removedtotest// delete his2D;
    //removedtotest// delete his1D;
    //
      if(fgForceTHnSparse)  pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(0,3); // A side
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_5") && !fgForceTHnSparse) {    
        his2D = dynamic_cast<TH2*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_5")->Clone());
    } else {
        his2D =  pTPCgain->GetDeDxHisto()->Projection(0,5);
    }         
    if(!his2D) return 4;

    TF1 * fpolA = new TF1("fpolA","pol1");
    TObjArray arrayFitA;
    //FitSlicesY(TF1* f1 = 0, Int_t firstxbin = 0, Int_t lastxbin = -1, Int_t cut = 0, Option_t* option = "QNR", TObjArray* arr = 0)   
    his2D->FitSlicesY(0,0,-1,10,"QN",&arrayFit); 
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpolA,"QNROB=0.8","QN",0,1);
    attachSlopeA = fpolA->GetParameter(1);
     //removedtotest// delete his2D;
     //removedtotest// delete his1D;
    //
    // MIP position vs. sector
    //
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-3,0); // C side
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_c_0_1") && !fgForceTHnSparse) {    
        his2D = dynamic_cast<TH2*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_c_0_1")->Clone());
    } else {
        his2D =  pTPCgain->GetDeDxHisto()->Projection(0,1);
    }
    if(!his2D) return 4;

    for(Int_t i = 0; i < 18; i++) { // loop over sectors; correct mapping to be checked!
        //TH1* his1D=0;
        Float_t phiLow = -TMath::Pi() + i*(20./360.)*(2*TMath::Pi());
        Float_t phiUp    = -TMath::Pi() + (i+1)*(20./360.)*(2*TMath::Pi());
        //pTPCgain->GetDeDxHisto()->GetAxis(1)->SetRangeUser(phiLow,phiUp);
        his2D->GetXaxis()->SetRangeUser(phiLow,phiUp);
        //his1D = pTPCgain->GetDeDxHisto()->Projection(0); 
        his1D = his2D->ProjectionY(); 
        TF1 gausFunc("gausFunc","gaus");
        his1D->Fit(&gausFunc, "QN");
        const Double_t mean=gausFunc.GetParameter(1);
        const Double_t res =gausFunc.GetParameter(2);
        meanMIPvsSector(i) = mean;
        sigmaRelMIPvsSector(i) = mean>0?res/mean:0.;
        sector(i)=i;
        //removedtotest// delete his1D;
    }
     //removedtotest// delete his2D;
    //
    if(fgForceTHnSparse)  pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(0,3); // A side
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_1") && !fgForceTHnSparse) {    
        his2D = dynamic_cast<TH2*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_1")->Clone());
    } else {
        his2D =  pTPCgain->GetDeDxHisto()->Projection(0,1);
    }    
    if(!his2D) return 4; 

    for(Int_t i = 0; i < 18; i++) { // loop over sectors; correct mapping to be checked!
        //TH1* his1D=0;
        Float_t phiLow = -TMath::Pi() + i*(20./360.)*(2*TMath::Pi());
        Float_t phiUp    = -TMath::Pi() + (i+1)*(20./360.)*(2*TMath::Pi());
        //pTPCgain->GetDeDxHisto()->GetAxis(1)->SetRangeUser(phiLow,phiUp);
        his2D->GetXaxis()->SetRangeUser(phiLow,phiUp);
        //his1D = pTPCgain->GetDeDxHisto()->Projection(0);
        his1D = his2D->ProjectionY();
        TF1 gausFunc("gausFunc","gaus");
        his1D->Fit(&gausFunc, "QN");
        const Double_t mean=gausFunc.GetParameter(1);
        const Double_t res =gausFunc.GetParameter(2);
        meanMIPvsSector(i+18) = mean;
        sigmaRelMIPvsSector(i+18) = mean>0?res/mean:0.;
        sector(i+18)=i+18;
        //removedtotest// delete his1D;
    }
     //removedtotest// delete his2D;

     //                                                         
    //  
    // select electrons                                                                                                               
    //                                                                                                                                                                           
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(7)->SetRangeUser(0.32,0.38); // momenta
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(0)->SetRangeUser(70,100); // dedx
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(6)->SetRangeUser(80,160); // nr clusters
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-1,1); // eta

    TF1 gausFitEle("gausFitEle","gaus");

    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mipsele_0") && !fgForceTHnSparse) {
      his1D = dynamic_cast<TH1*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mipsele_0")->Clone());
    } else {
      his1D =  pTPCgain->GetDeDxHisto()->Projection(0);
    }
    if(!his1D) return 4;
    his1D->Fit(&gausFitEle,"QN","QN");

    meanMIPele = gausFitEle.GetParameter(1);
    resolutionMIPele = 0;
    if (meanMIPele!=0) resolutionMIPele = gausFitEle.GetParameter(2)/meanMIPele;
    
    //restore cuts as before
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(7)->SetRangeUser(0.4,0.55);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(0)->SetRangeUser(35,60);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(6)->SetRangeUser(80,160);
    if(fgForceTHnSparse) pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-1,1);

    //                                                                                                                                                                        
    // separation between electrons and MIPs                                                                                                                                  
    // 
    electroMIPSeparation = TMath::Abs((meanMIP-meanMIPele));

    printf("Gain QA report\n");
    printf("MIP mean\t%f\n",meanMIP);
    printf("MIP resolution\t%f\n",resolutionMIP);
    printf("MIPslopeA\t%f\n",attachSlopeA);
    printf("MIPslopeC\t%f\n",attachSlopeC);
    printf("Electons energy loss MIP mean\t%f\n",meanMIPele);
    printf("Electons MIP resolution\t%f\n",resolutionMIPele);
    // 
    
    (*pcstream)<<"tpcQA"<<
      "MIPattachSlopeC="<<attachSlopeC<<
      "MIPattachSlopeA="<<attachSlopeA<<
      "resolutionMIP="<<resolutionMIP<<
      "meanMIPvsSector.="<<&meanMIPvsSector<<
      "sigmaRelMIPvsSector.="<<&sigmaRelMIPvsSector<<
      "sector.="<<&sector<<
      "meanMIP="<<meanMIP<<
      "meanMIPele="<<meanMIPele<<
      "resolutionMIPele="<<resolutionMIPele<<
      "electroMIPSeparation="<<electroMIPSeparation;
    
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzeEvent(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse Primary Vertex Distribution and Multiplicities
    //
  if (!pcstream) return 1;
    if (!pTPC) return 1;
    //
    // 
    //
    static TVectorF infoVertexX(5);
    static TVectorF infoVertexY(5);
    static TVectorF infoVertexZ(5);
    static TVectorF infoMult(5);
    static TVectorF infoMultPos(5);
    static TVectorF infoMultNeg(5);
    
    static Double_t entriesVertX=0;
    static Double_t meanVertX=0;
    static Double_t rmsVertX=0;
    static Double_t entriesVertY=0;
    static Double_t meanVertY=0;
    static Double_t rmsVertY=0;
    static Double_t entriesVertZ=0;
    static Double_t meanVertZ=0;
    static Double_t rmsVertZ=0;
    static Double_t vertStatus=0;
    static Double_t entriesMult=0;
    static Double_t meanMult=0;
    static Double_t rmsMult=0;
    static Double_t meanMultPos=0;
    static Double_t rmsMultPos=0;
    static Double_t errorMultPos=0;
    static Double_t meanMultNeg=0;
    static Double_t rmsMultNeg=0;
    static Double_t errorMultNeg=0;
    static Double_t vertAll = 0;
    static Double_t vertOK = 0;
    
    TH1* his1D=0;
    TH1* hc=0;
    if (pTPC->GetHistos()->FindObject("h_tpc_event_6") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_6")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(6);
    }
    if(!his1D) return 1;

    vertAll = his1D->GetEntries();
    vertOK  = his1D->GetBinContent(2);
    if (vertAll>=1) {
            vertStatus = vertOK / vertAll;
    }
    delete his1D;
    
    if(fgForceTHnSparse) pTPC->GetTPCEventHisto()->GetAxis(6)->SetRange(2,2);
   
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_0") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_0")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(0);
    }
    if(!his1D) return 1;
    entriesVertX = his1D->GetEntries(); 
    meanVertX = his1D->GetMean();    
    rmsVertX    = his1D->GetRMS();
    delete his1D;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_1") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_1")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(1);
    }
    if(!his1D) return 1;

    entriesVertY = his1D->GetEntries();
    meanVertY = his1D->GetMean();
    rmsVertY    = his1D->GetRMS();
    delete his1D;
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_2") && !fgForceTHnSparse) {    
        hc = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_2"));
        if(!hc) return 1;
        his1D = (TH1*)hc->Clone();
    }
    else {
       his1D = pTPC->GetTPCEventHisto()->Projection(2);
    }    
    if(!his1D) return 1;

    entriesVertZ = his1D->GetEntries();
    meanVertZ = his1D->GetMean();
    rmsVertZ    = his1D->GetRMS();
    delete his1D;
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_3") && !fgForceTHnSparse) {    
        hc = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_3"));
	if(!hc) return 1;
        his1D = (TH1*)hc->Clone();
    }
    else {
       his1D = pTPC->GetTPCEventHisto()->Projection(3);
    }
    if(!his1D) return 1;

    entriesMult = his1D->GetEntries();
    meanMult    = his1D->GetMean();
    rmsMult     = his1D->GetRMS();
    delete his1D;
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_4") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_4")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(4);
    }
    if(!his1D) return 1;

    meanMultPos    = his1D->GetMean();
    rmsMultPos     = his1D->GetRMS();
    errorMultPos   = his1D->GetRMS() / TMath::Sqrt(his1D->GetEntries());
    GetStatInfo(his1D,infoMultPos,0);

    delete his1D;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_5") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_5")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(5);
    }
    if(!his1D) return 1;

    meanMultNeg    = his1D->GetMean();
    rmsMultNeg     = his1D->GetRMS();
    errorMultNeg   = his1D->GetRMS() / TMath::Sqrt(his1D->GetEntries());
    GetStatInfo(his1D,infoMultNeg,0);

    delete his1D;
    //
    (*pcstream)<<"trending"<<
      "infoVertX.="<<&infoVertexX <<       // vertex X stat information
      "infoVertY.="<<&infoVertexY <<       // vertex Y stat information
      "infoVertZ.="<<&infoVertexZ <<       // vertex Z stat information
      "infoMult.="<<&infoMult <<         // multipicity stat information
      "infoMultPos.="<<&infoMultPos <<   // multiplicity stat information
      "infoMultNeg.="<<&infoMultPos <<   // multiplicity stat information
      
      "entriesVertX="<<entriesVertX<<
      "meanVertX="<<meanVertX<<
      "rmsVertX="<<rmsVertX<<
      "entriesVertY="<<entriesVertY<<
      "meanVertY="<<meanVertY<<
      "rmsVertY="<<rmsVertY<<
      "entriesVertZ="<<entriesVertZ<<
      "meanVertZ="<<meanVertZ<<
      "rmsVertZ="<<rmsVertZ<<
      "vertStatus="<<vertStatus<<
      "vertAll="<<vertAll<<
      "vertOK="<<vertOK<<
      "entriesMult="<<entriesMult<<
      "meanMult="<<meanMult<<
      "rmsMult="<<rmsMult<<
      "meanMultPos="<<meanMultPos<<
      "rmsMultPos="<<rmsMultPos<<
      "errorMultPos="<<errorMultPos<<
      "meanMultNeg="<<meanMultNeg<<
      "rmsMultNeg="<<rmsMultNeg<<          
      "errorMultNeg="<<errorMultNeg;
    
    if(fgForceTHnSparse) pTPC->GetTPCEventHisto()->GetAxis(6)->SetRange(1,2);
    //
    (*pcstream)<<"tpcQA"<<
        "entriesVertX="<<entriesVertX<<
        "meanVertX="<<meanVertX<<
        "rmsVertX="<<rmsVertX<<
        "entriesVertY="<<entriesVertY<<
        "meanVertY="<<meanVertY<<
        "rmsVertY="<<rmsVertY<<
        "entriesVertZ="<<entriesVertZ<<
        "meanVertZ="<<meanVertZ<<
        "rmsVertZ="<<rmsVertZ<<
        "vertStatus="<<vertStatus<<
        "vertAll="<<vertAll<<
        "vertOK="<<vertOK<<
        "entriesMult="<<entriesMult<<
        "meanMult="<<meanMult<<
        "rmsMult="<<rmsMult<<
        "meanMultPos="<<meanMultPos<<
        "rmsMultPos="<<rmsMultPos<<
        "meanMultNeg="<<meanMultNeg<<
        "rmsMultNeg="<<rmsMultNeg;     
     
    return 0;
}

//_____________________________________________________________________________
Int_t AliTPCPerformanceSummary::AnalyzePt(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
    //
    // Analyse DCA R imperfections for positive particles
    //
    
    if (!pcstream) return 256;
    if (!pTPC) return 256;

    // variables:
    static Double_t meanPtAPos = 0;
    static Double_t mediumPtAPos = 0;
    static Double_t highPtAPos = 0;
    static Double_t meanPtCPos = 0;
    static Double_t mediumPtCPos = 0;
    static Double_t highPtCPos = 0;

    static Double_t meanPtANeg = 0;
    static Double_t mediumPtANeg = 0;
    static Double_t highPtANeg = 0;
    static Double_t meanPtCNeg = 0;
    static Double_t mediumPtCNeg = 0;
    static Double_t highPtCNeg = 0;

    TH3* his3D1=0;
    TH3* his3D2=0;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_7")) {    

      his3D1 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_7"));
      if(!his3D1) return 256;
	
      his3D1->GetYaxis()->SetRangeUser(0.1,0.8);
      
      his3D1->GetZaxis()->SetRangeUser(0.25,10);
      meanPtAPos = his3D1->GetMean(3);
      his3D1->GetZaxis()->SetRangeUser(2,5);
      mediumPtAPos = his3D1->GetMean(3);
      his3D1->GetZaxis()->SetRangeUser(5,10);
      highPtAPos = his3D1->GetMean(3);
      
      his3D1->GetYaxis()->SetRangeUser(-0.8,-0.1);

      his3D1->GetZaxis()->SetRangeUser(0.25,10);
      meanPtCPos = his3D1->GetMean(3);
      his3D1->GetZaxis()->SetRangeUser(2,5);
      mediumPtCPos = his3D1->GetMean(3);
      his3D1->GetZaxis()->SetRangeUser(5,10);
      highPtCPos = his3D1->GetMean(3);

      his3D1->GetYaxis()->SetRangeUser(-1,1);
      his3D1->GetZaxis()->SetRangeUser(0.25,10);
    }


    if (pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_7")) {    

      his3D2 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_7"));
      if(!his3D2) return 256;
	
      his3D2->GetYaxis()->SetRangeUser(0.1,0.8);

      his3D2->GetZaxis()->SetRangeUser(0.25,10);
      meanPtANeg = his3D2->GetMean(3);
      his3D2->GetZaxis()->SetRangeUser(2,5);
      mediumPtANeg = his3D2->GetMean(3);
      his3D2->GetZaxis()->SetRangeUser(5,10);
      highPtANeg = his3D2->GetMean(3);
      
      his3D2->GetYaxis()->SetRangeUser(-0.8,-0.1);

      his3D2->GetZaxis()->SetRangeUser(0.25,10);
      meanPtCNeg = his3D2->GetMean(3);
      his3D2->GetZaxis()->SetRangeUser(2,5);
      mediumPtCNeg = his3D2->GetMean(3);
      his3D2->GetZaxis()->SetRangeUser(5,10);
      highPtCNeg = his3D2->GetMean(3);
      
      his3D2->GetYaxis()->SetRangeUser(-1,1);
      his3D2->GetZaxis()->SetRangeUser(0.25,10);
    }



    // dump values
    //
    (*pcstream)<<"tpcQA"<<
      "meanPtAPos="<< meanPtAPos<<
      "mediumPtAPos="<< mediumPtAPos<<
      "highPtAPos="<< highPtAPos<<
      //
      "meanPtCPos="<< meanPtCPos<<
      "mediumPtCPos="<< mediumPtCPos<<
      "highPtCPos="<< highPtCPos<<
      //
      "meanPtANeg="<< meanPtANeg<<
      "mediumPtANeg="<< mediumPtANeg<<
      "highPtANeg="<< highPtANeg<<
        //
      "meanPtCNeg="<< meanPtCNeg<<
      "mediumPtCNeg="<< mediumPtCNeg<<
      "highPtCNeg="<< highPtCNeg;

        
    return 0;
}

//_____________________________________________________________________________

Int_t AliTPCPerformanceSummary::AnalyzeChargeOverPt(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream){
    //
    // Analyse DCA R imperfections for positive particles
    //
    
    if (!pcstream) return 512;
    if (!pTPC) return 512;

    // variables:
    static Double_t qOverPt = 0;
    static Double_t qOverPtA = 0;
    static Double_t qOverPtC = 0;

    TH2* his2D=0;
    TH1* his1D1=0;
    TH1* his1D2=0;
    TH1* his1D3=0;
    TF1 *fp1 = new TF1("fp1","pol2",-1.0,1.0);
    TF1 *fp2 = new TF1("fp2","pol2",-1.0,1.0);
    TF1 *fp3 = new TF1("fp3","pol2",-1.0,1.0);
    
    if (pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_5_8")) {

      his2D = dynamic_cast<TH2*>(pTPC->GetHistos()->FindObject("h_tpc_track_all_recvertex_5_8"));
      if(!his2D) return 512;

      his1D1 = his2D->ProjectionX();
      his1D1->Fit(fp1,"R");
      if(fp1->GetParameter(2)!=0){
	qOverPt = (-1.0)*(fp1->GetParameter(1)/(2.0*fp1->GetParameter(2)));
       }
      delete fp1;
      delete his1D1;

       his2D->GetYaxis()->SetRangeUser(0.1,0.8);
       his1D2 = his2D->ProjectionX();
       his1D2->Fit(fp2,"R");
       if(fp2->GetParameter(2)!=0)
	 qOverPtA = (-1.0)*(fp2->GetParameter(1)/(2.0*fp2->GetParameter(2)));
       delete fp2;
       delete his1D2;
     
       his2D->GetYaxis()->SetRangeUser(-0.8,-0.1);       
       his1D3 = his2D->ProjectionX();
       his1D3->Fit(fp3,"R");
       if(fp3->GetParameter(2)!=0)
	 qOverPtC = (-1.0)*(fp3->GetParameter(1)/(2.0*fp3->GetParameter(2)));
       delete fp3;
       delete his1D3;
       
      his2D->GetYaxis()->SetRangeUser(-1.0,1.0);
    }
    
    
    (*pcstream)<<"tpcQA"<<
      "qOverPt="<< qOverPt<<
      "qOverPtA="<< qOverPtA<<
      "qOverPtC="<< qOverPtC;
    
    return 0;
}

Int_t AliTPCPerformanceSummary::AnalyzeMatch(const AliPerformanceMatch* pMatch, TTreeSRedirector* const pcstream)
{
  /* if ((pMatch == 0) or (0 == pcstream)) { printf("this will not work anyway..."); }
     printf("funtion not implemented");*/

  if (!pcstream) return 1024;
  if (!pMatch) return 1024;
  static Double_t tpcItsMatchA = 0;
  static Double_t tpcItsMatchHighPtA = 0; 
  static Double_t tpcItsMatchC = 0;
  static Double_t tpcItsMatchHighPtC = 0; 

  TH2 *h2D = 0;
  TH2 *h2D1 = 0;
  if(pMatch->GetHistos()->FindObject("h_tpc_match_trackingeff_all_2_3") &&
     pMatch->GetHistos()->FindObject("h_tpc_match_trackingeff_tpc_2_3")){
    h2D = dynamic_cast<TH2*>(pMatch->GetHistos()->FindObject("h_tpc_match_trackingeff_all_2_3"));
    h2D1 = dynamic_cast<TH2*>(pMatch->GetHistos()->FindObject("h_tpc_match_trackingeff_tpc_2_3"));
   
    if(!h2D) return 4;
    if(!h2D1) return 4;

    h2D->GetXaxis()->SetRangeUser(0,1.5);
    h2D1->GetXaxis()->SetRangeUser(0,1.5);

    Double_t entries,entries1;
    entries = h2D->GetEffectiveEntries();
    entries1 = h2D1->GetEffectiveEntries();
    if(entries > 0)
      tpcItsMatchA = entries1/entries;

    h2D->GetYaxis()->SetRangeUser(4.01,20.);
    h2D1->GetYaxis()->SetRangeUser(4.01,20.);
    entries = h2D->GetEffectiveEntries();
    entries1 = h2D1->GetEffectiveEntries();
    if(entries > 0)
    tpcItsMatchHighPtA = entries1/entries;


    h2D->GetXaxis()->SetRangeUser(-1.5,-0.01);
    h2D1->GetXaxis()->SetRangeUser(-1.5,-0.01);
    h2D->GetYaxis()->SetRangeUser(0.0,20.);
    h2D1->GetYaxis()->SetRangeUser(0.0,20.);

    entries = h2D->GetEffectiveEntries();
    entries1 = h2D1->GetEffectiveEntries();
    if(entries > 0)
      tpcItsMatchC = entries1/entries;

    h2D->GetYaxis()->SetRangeUser(4.01,20.);
    h2D1->GetYaxis()->SetRangeUser(4.01,20.);
    entries = h2D->GetEffectiveEntries();
    entries1 = h2D1->GetEffectiveEntries();
    if(entries > 0)
      tpcItsMatchHighPtC = entries1/entries;

    h2D->GetXaxis()->SetRangeUser(-1.5,1.5);
    h2D1->GetXaxis()->SetRangeUser(-1.5,1.5);
    h2D->GetYaxis()->SetRangeUser(0.0,20.);
    h2D1->GetYaxis()->SetRangeUser(0.0,20.);
    //    delete h2D;
    //    delete h2D1;
  }

  (*pcstream)<<"tpcQA"<<
    "tpcItsMatchA="<< tpcItsMatchA<<
    "tpcItsMatchHighPtA="<< tpcItsMatchHighPtA<<
    "tpcItsMatchC="<< tpcItsMatchC<<
    "tpcItsMatchHighPtC="<< tpcItsMatchHighPtC;

  return 0;
}

Int_t AliTPCPerformanceSummary::AnalyzePull(const AliPerformanceMatch* pPull, TTreeSRedirector* const pcstream)
{
  /* if ((pPull == 0) or (0 == pcstream)) { printf("this will not work anyway..."); }
     printf("funtion not implemented");*/

  if (!pcstream) return 2048;
  if (!pPull) return 2048;
  static Double_t phiPull = 0;
  static Double_t phiPullHighPt = 0; 
  static Double_t ptPull = 0;
  static Double_t ptPullHighPt = 0; 
  static Double_t yPull = 0;
  static Double_t yPullHighPt = 0; 
  static Double_t zPull = 0;
  static Double_t zPullHighPt = 0; 
  static Double_t lambdaPull = 0;
  static Double_t lambdaPullHighPt = 0; 

  TH2 *h2D1 = 0;
  if(pPull->GetHistos()->FindObject("h_tpc_match_pull_2_7")){
    h2D1 = dynamic_cast<TH2*>(pPull->GetHistos()->FindObject("h_tpc_match_pull_2_7"));
    if(!h2D1) return 4;
    phiPull = h2D1->GetMean(2);
    h2D1->SetAxisRange(0.0,1.0/5.0,"X");
    phiPullHighPt = h2D1->GetMean(2);
    h2D1->SetAxisRange(0.0,10.0,"X");
    //    delete h2D1;
  }

  TH2 *h2D2 = 0;
  if(pPull->GetHistos()->FindObject("h_tpc_match_pull_4_7")){
    h2D2 = dynamic_cast<TH2*>(pPull->GetHistos()->FindObject("h_tpc_match_pull_4_7"));
    if(!h2D2) return 4;
    ptPull = h2D2->GetMean(2);

    h2D2->SetAxisRange(0.0,1.0/5.0,"X");
    ptPullHighPt = h2D2->GetMean(2);
    h2D2->SetAxisRange(0.0,10.0,"X");
    //    delete h2D2;
  }

  TH2 *h2D3 = 0;
  if(pPull->GetHistos()->FindObject("h_tpc_match_pull_0_7")){
    h2D3 = dynamic_cast<TH2*>(pPull->GetHistos()->FindObject("h_tpc_match_pull_0_7"));
    if(!h2D3) return 4;
    yPull = h2D3->GetMean(2);

    h2D3->SetAxisRange(0.0,1.0/5.0,"X");
    yPullHighPt = h2D3->GetMean(2);
    h2D3->SetAxisRange(0.0,10.0,"X");
    //    delete h2D3;
  }

  TH2 *h2D4 = 0;
  if(pPull->GetHistos()->FindObject("h_tpc_match_pull_1_7")){
    h2D4 = dynamic_cast<TH2*>(pPull->GetHistos()->FindObject("h_tpc_match_pull_1_7"));
    if(!h2D4) return 4;
    zPull = h2D4->GetMean(2);

    h2D4->SetAxisRange(0.0,1.0/5.0,"X");
    zPullHighPt = h2D4->GetMean(2);
    h2D4->SetAxisRange(0.0,10.0,"X");
    //    delete h2D4;
 }

  TH2 *h2D5 = 0;
  if(pPull->GetHistos()->FindObject("h_tpc_match_pull_3_7")){
    h2D5 = dynamic_cast<TH2*>(pPull->GetHistos()->FindObject("h_tpc_match_pull_3_7"));
    if(!h2D5) return 4;
    lambdaPull = h2D5->GetMean(2);

    h2D5->SetAxisRange(0.0,1.0/5.0,"X");
    lambdaPullHighPt = h2D5->GetMean(2);
    h2D5->SetAxisRange(0.0,10.0,"X");
    //    delete h2D5;
}

  (*pcstream)<<"tpcQA"<<
    "phiPull="<< phiPull<<
    "phiPullHighPt="<< phiPullHighPt<<
    "ptPull="<< ptPull<<
    "ptPullHighPt="<< ptPullHighPt<<
    "yPull="<< yPull<<
    "yPullHighPt="<< yPullHighPt<<
    "zPull="<< zPull<<
    "zPullHighPt="<< zPullHighPt<<
    "lambdaPull="<< lambdaPull<<
    "lambdaPullHighPt="<< lambdaPullHighPt;
    
  return 0;
}
Int_t AliTPCPerformanceSummary::AnalyzeConstrain(const AliPerformanceMatch* pConstrain, TTreeSRedirector* pcstream)
{
  if (!pcstream) return 5126;
  if (!pConstrain) return 5126;

    TH3* his3D=0;
    static Double_t tpcConstrainPhiA = 0;
    static Double_t tpcConstrainPhiC = 0;
    
    if (pConstrain->GetHistos()->FindObject("h_tpc_constrain_tpc_0_2_3")) {    
      
      his3D = dynamic_cast<TH3*>(pConstrain->GetHistos()->FindObject("h_tpc_constrain_tpc_0_2_3"));//phi pull:pt:eta
      if(!his3D) return 5126;
      
      his3D->GetZaxis()->SetRangeUser(0.0,1.0);
      tpcConstrainPhiA = his3D->GetMean(1);
      his3D->GetZaxis()->SetRangeUser(-1.0,-0.001);
      tpcConstrainPhiC = his3D->GetMean(1);
    }

  (*pcstream)<<"tpcQA"<<
    "tpcConstrainPhiA="<<tpcConstrainPhiA <<
    "tpcConstrainPhiC="<< tpcConstrainPhiC;
  
  return 0;
}

//_____________________________________________________________________________
 Int_t AliTPCPerformanceSummary::AnalyzeQAPosNegDpT(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
  //function which plot 1/Pt for negative and 
  //positive particles
  
  if (!pcstream) return 512;
  if (!pTPC) return 512;
  
  TH3D* pos3=0;
  TH3D* neg3=0;
  TH1D* pos=0;
  TH1D* neg=0;
  TH1D* posC=0;
  TH1D* negC=0;
  TH1D* posA=0;
  TH1D* negA=0;
  static Double_t deltaPtC = 0;
  static Double_t deltaPtchi2C = 0;
  static Double_t slopeC = 0;
  static Double_t deltaPtA = 0;
  static Double_t deltaPtchi2A = 0;
  static Double_t slopeA = 0;
  static Double_t deltaPt = 0;
  static Double_t deltaPtchi2 = 0;
  static Double_t slope = 0;
  static Double_t deltaPt_Err = 0; 
  static Double_t deltaPtA_Err = 0;
  static Double_t deltaPtC_Err = 0;


//C side

  if(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_0_5_7"))
    {
    pos3 = dynamic_cast<TH3D*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_0_5_7")); 
    if(!pos3) return 512;
  
    pos = pos3->ProjectionZ("pos",71,-1,6,25);
    posC = pos3->ProjectionZ("posC",71,-1,6,15);
    posA = pos3->ProjectionZ("posA",71,-1,16,25);
    }

    if(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_0_5_7")){
    neg3 = dynamic_cast<TH3D*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_0_5_7")); 
    if(!neg3) return 512;
    
    neg = neg3->ProjectionZ("neg",71,-1,6,25);
    negC = neg3->ProjectionZ("negC",71,-1,6,15);
    negA = neg3->ProjectionZ("negA",71,-1,16,25);
}

if(!pos) return 512; 
if(!neg) return 512; 
if(!posA) return 512; 
if(!negA) return 512; 
if(!posC) return 512; 
if(!negC) return 512; 

pos->Sumw2();
neg->Sumw2();
posA->Sumw2();
negA->Sumw2();
posC->Sumw2();
negC->Sumw2();

pos->Scale(1.,"width");
neg->Scale(1.,"width");
posA->Scale(1.,"width");
negA->Scale(1.,"width");
posC->Scale(1.,"width");
negC->Scale(1.,"width");

//both sides

TF1 fpt("fpt","[1]*exp(-1/((1/x))*[0])",0.1,10);
TF1 fpt2("fpt2","[1]*exp(-1/((1/x))*[0])",0.1,10);
fpt.SetParameters(1,0.5);
fpt2.SetParameters(1,0.5);
pos->Fit(&fpt,"","",1,4); pos->Fit(&fpt,"","",1,4); pos->Fit(&fpt,"","",1,4);
neg->Fit(&fpt2,"","",1,4); neg->Fit(&fpt2,"","",1,4); neg->Fit(&fpt2,"","",1,4);

slope = (fpt.GetParameter(0)+fpt2.GetParameter(0))/2.;

TH1D* ratio = new TH1D(*pos); 
ratio->Divide(neg);

ratio->Draw();
TF1 fptRatio("fptratio","[2]*exp(-1/((1/x)+[1])*[0])/exp(-1/((1/x)-[1])*[0])",0.1,10);
fptRatio.SetParameters(0.5,0.006,1);
fptRatio.FixParameter(0,slope);
fptRatio.Draw();
ratio->Fit(&fptRatio,"","",1,4); ratio->Fit(&fptRatio,"","",1,4); 
ratio->Fit(&fptRatio,"","",1,4);

deltaPt = fptRatio.GetParameter(1);
deltaPtchi2 = fptRatio.GetChisquare();

//get the errors
deltaPt_Err = fptRatio.GetParError(1);


//A side

TF1 fptA("fptA","[1]*exp(-1/((1/x))*[0])",0.1,10);
TF1 fpt2A("fpt2A","[1]*exp(-1/((1/x))*[0])",0.1,10);
fptA.SetParameters(1,0.5);
fpt2A.SetParameters(1,0.5);
posA->Fit(&fptA,"","",1,4); posA->Fit(&fptA,"","",1,4); posA->Fit(&fptA,"","",1,4);
negA->Fit(&fpt2A,"","",1,4); negA->Fit(&fpt2A,"","",1,4); negA->Fit(&fpt2A,"","",1,4);

slopeA = (fptA.GetParameter(0)+fpt2A.GetParameter(0))/2.;

TH1D* ratioA = new TH1D(*posA); 
ratioA->Divide(negA);

ratioA->Draw();
TF1 fptRatioA("fptratioA","[2]*exp(-1/((1/x)+[1])*[0])/exp(-1/((1/x)-[1])*[0])",0.1,10);
fptRatioA.SetParameters(0.5,0.006,1);
fptRatioA.FixParameter(0,slopeA);
fptRatioA.Draw();
ratioA->Fit(&fptRatioA,"","",1,4); ratio->Fit(&fptRatioA,"","",1,4); 
ratioA->Fit(&fptRatioA,"","",1,4);

deltaPtA = fptRatioA.GetParameter(1);
deltaPtchi2A = fptRatioA.GetChisquare();

//get the errors
deltaPtA_Err = fptRatioA.GetParError(1);

 delete ratioA;
 delete pos;
 delete neg;


//C side
TF1 fptC("fptC","[1]*exp(-1/((1/x))*[0])",0.1,10);
TF1 fpt2C("fpt2C","[1]*exp(-1/((1/x))*[0])",0.1,10);
fptC.SetParameters(1,0.5);
fpt2C.SetParameters(1,0.5);
posC->Fit(&fptC,"","",1,4); posC->Fit(&fptC,"","",1,4); posC->Fit(&fptC,"","",1,4);
negC->Fit(&fpt2C,"","",1,4); negC->Fit(&fpt2C,"","",1,4); negC->Fit(&fpt2C,"","",1,4);

slopeC = (fptC.GetParameter(0)+fpt2C.GetParameter(0))/2.;

TH1D* ratioC = new TH1D(*posC); 
ratioC->Divide(negC);

ratioC->Draw();
TF1 fptRatioC("fptratioC","[2]*exp(-1/((1/x)+[1])*[0])/exp(-1/((1/x)-[1])*[0])",0.1,10);
fptRatioC.SetParameters(0.5,0.006,1);
fptRatioC.FixParameter(0,slopeC);
fptRatioC.Draw();
ratioC->Fit(&fptRatioC,"","",1,4); ratio->Fit(&fptRatioC,"","",1,4); 
ratioC->Fit(&fptRatioC,"","",1,4);

deltaPtC = fptRatioC.GetParameter(1);
deltaPtchi2C = fptRatioC.GetChisquare();

//get the errors
deltaPtC_Err = fptRatioC.GetParError(1);


 delete posC;
 delete negC;
 delete ratioC;
    
    (*pcstream)<<"tpcQA"<<      
      "deltaPt="<< deltaPt<<
      "deltaPtchi2="<< deltaPtchi2<<
      "deltaPtA="<< deltaPtA<<
      "deltaPtchi2A="<< deltaPtchi2A<<
      "deltaPtC="<< deltaPtC<<
      "deltaPtchi2C="<< deltaPtchi2C<<
      "deltaPt_Err="<< deltaPt_Err<<
      "deltaPtA_Err="<< deltaPtA_Err<<
      "deltaPtC_Err="<< deltaPtC_Err;    
      
    return 0;
}

//_____________________________________________________________________________
 Int_t AliTPCPerformanceSummary::AnalyzeQADCAFitParameter(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
  
  //
  //function which retrieve DCA fit parameters
  //
  
  if (!pcstream) return 16;
  if (!pTPC) return 16;
  
  TH3* dcar_pos3=0;
  TH3* dcaz_pos3=0;
  TH3* dcar_neg3=0;
  TH3* dcaz_neg3=0;
  static TGraphErrors * graphEtaProfile[4]={0};
  static TGraphErrors * graphPhiProfile[8]={0};
  static Double_t dcar_posA_0=0; 
  static Double_t dcar_posA_1=0; 
  static Double_t dcar_posA_2=0; 
  static Double_t dcar_posA_chi2=0; 
  static Double_t dcar_posA_0_Err=0; 
  static Double_t dcar_posA_1_Err=0; 
  static Double_t dcar_posA_2_Err=0; 
  
  static Double_t dcar_posC_0=0; 
  static Double_t dcar_posC_1=0; 
  static Double_t dcar_posC_2=0; 
  static Double_t dcar_posC_chi2=0; 
  static Double_t dcar_posC_0_Err=0; 
  static Double_t dcar_posC_1_Err=0; 
  static Double_t dcar_posC_2_Err=0; 

  static Double_t dcaz_posA_0=0; 
  static Double_t dcaz_posA_1=0; 
  static Double_t dcaz_posA_2=0; 
  static Double_t dcaz_posA_chi2=0; 
  static Double_t dcaz_posA_0_Err=0; 
  static Double_t dcaz_posA_1_Err=0; 
  static Double_t dcaz_posA_2_Err=0; 
  
  static Double_t dcaz_posC_0=0; 
  static Double_t dcaz_posC_1=0; 
  static Double_t dcaz_posC_2=0; 
  static Double_t dcaz_posC_chi2=0; 
  static Double_t dcaz_posC_0_Err=0; 
  static Double_t dcaz_posC_1_Err=0; 
  static Double_t dcaz_posC_2_Err=0; 
  
  static Double_t dcar_negA_0=0; 
  static Double_t dcar_negA_1=0; 
  static Double_t dcar_negA_2=0; 
  static Double_t dcar_negA_chi2=0; 
  static Double_t dcar_negA_0_Err=0; 
  static Double_t dcar_negA_1_Err=0; 
  static Double_t dcar_negA_2_Err=0; 
  
  static Double_t dcar_negC_0=0; 
  static Double_t dcar_negC_1=0; 
  static Double_t dcar_negC_2=0; 
  static Double_t dcar_negC_chi2=0; 
  static Double_t dcar_negC_0_Err=0; 
  static Double_t dcar_negC_1_Err=0; 
  static Double_t dcar_negC_2_Err=0; 

  static Double_t dcaz_negA_0=0; 
  static Double_t dcaz_negA_1=0; 
  static Double_t dcaz_negA_2=0; 
  static Double_t dcaz_negA_chi2=0; 
  static Double_t dcaz_negA_0_Err=0; 
  static Double_t dcaz_negA_1_Err=0; 
  static Double_t dcaz_negA_2_Err=0; 
 
  static Double_t dcaz_negC_0=0; 
  static Double_t dcaz_negC_1=0; 
  static Double_t dcaz_negC_2=0; 
  static Double_t dcaz_negC_chi2=0; 
  static Double_t dcaz_negC_0_Err=0; 
  static Double_t dcaz_negC_1_Err=0; 
  static Double_t dcaz_negC_2_Err=0;
 
  if (pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_6")) {  
    dcar_pos3 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_3_5_6"));
  }
  
  if (pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_4_5_6")) {
    dcaz_pos3 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_pos_recvertex_4_5_6"));
  }
  
  if (pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_6")) {
    dcar_neg3 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_3_5_6")); 
  }
  
  if (pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_4_5_6")) {
    dcaz_neg3 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_track_neg_recvertex_4_5_6"));
  }

  TF1 fit("fit","[0]+[1]*cos(x)+[2]*sin(x)",0,7); 

  dcar_pos3->GetYaxis()->SetRangeUser(0,0.99);    
  TH1* dcar_posA = (dynamic_cast<TH2*>(dcar_pos3->Project3D("xz_1")))->ProfileX(); 
  dcar_posA->Fit(&fit,"NQ");
  dcar_posA_0 = fit.GetParameter(0);
  dcar_posA_1 = fit.GetParameter(1);
  dcar_posA_2 = fit.GetParameter(2);
  dcar_posA_chi2 = fit.GetChisquare();  
  dcar_posA_0_Err = fit.GetParError(0);
  dcar_posA_1_Err = fit.GetParError(1);
  dcar_posA_2_Err = fit.GetParError(2);
  
  dcar_pos3->GetYaxis()->SetRangeUser(-1.0,-0.01);    
  TH1* dcar_posC = (dynamic_cast<TH2*>(dcar_pos3->Project3D("xz_2")))->ProfileX(); 
  dcar_posC->Fit(&fit,"NQ");
  dcar_posC_0 = fit.GetParameter(0);
  dcar_posC_1 = fit.GetParameter(1);
  dcar_posC_2 = fit.GetParameter(2);
  dcar_posC_chi2 = fit.GetChisquare();    
  dcar_posC_0_Err = fit.GetParError(0);
  dcar_posC_1_Err = fit.GetParError(1);
  dcar_posC_2_Err = fit.GetParError(2);

  dcaz_pos3->GetYaxis()->SetRangeUser(0,0.99);    
  TH1* dcaz_posA = (dynamic_cast<TH2*>(dcaz_pos3->Project3D("xz_3")))->ProfileX(); 
  dcaz_posA->Fit(&fit,"NQ");
  dcaz_posA_0 = fit.GetParameter(0);
  dcaz_posA_1 = fit.GetParameter(1);
  dcaz_posA_2 = fit.GetParameter(2);
  dcaz_posA_chi2 = fit.GetChisquare();      
  dcaz_posA_0_Err = fit.GetParError(0);
  dcaz_posA_1_Err = fit.GetParError(1);
  dcaz_posA_2_Err = fit.GetParError(2);  
  
  dcaz_pos3->GetYaxis()->SetRangeUser(-1.0,-0.01);    
  TH1* dcaz_posC = (dynamic_cast<TH2*>(dcaz_pos3->Project3D("xz_4")))->ProfileX(); 
  dcaz_posC->Fit(&fit,"NQ");
  dcaz_posC_0 = fit.GetParameter(0);
  dcaz_posC_1 = fit.GetParameter(1);
  dcaz_posC_2 = fit.GetParameter(2);
  dcaz_posC_chi2 = fit.GetChisquare();    
  dcaz_posC_0_Err = fit.GetParError(0);
  dcaz_posC_1_Err = fit.GetParError(1);
  dcaz_posC_2_Err = fit.GetParError(2);  
    

  
   dcar_neg3->GetYaxis()->SetRangeUser(0,0.99);    
  TH1* dcar_negA = (dynamic_cast<TH2*>(dcar_neg3->Project3D("xz_1")))->ProfileX(); 
  dcar_negA->Fit(&fit,"NQ");
  dcar_negA_0 = fit.GetParameter(0);
  dcar_negA_1 = fit.GetParameter(1);
  dcar_negA_2 = fit.GetParameter(2);
  dcar_negA_chi2 = fit.GetChisquare();  
  dcar_negA_0_Err = fit.GetParError(0);
  dcar_negA_1_Err = fit.GetParError(1);
  dcar_negA_2_Err = fit.GetParError(2);
  
  dcar_neg3->GetYaxis()->SetRangeUser(-1.0,-0.01);    
  TH1* dcar_negC = (dynamic_cast<TH2*>(dcar_neg3->Project3D("xz_2")))->ProfileX(); 
  dcar_negC->Fit(&fit,"NQ");
  dcar_negC_0 = fit.GetParameter(0);
  dcar_negC_1 = fit.GetParameter(1);
  dcar_negC_2 = fit.GetParameter(2);
  dcar_negC_chi2 = fit.GetChisquare();    
  dcar_negC_0_Err = fit.GetParError(0);
  dcar_negC_1_Err = fit.GetParError(1);
  dcar_negC_2_Err = fit.GetParError(2);

  dcaz_neg3->GetYaxis()->SetRangeUser(0,0.99);    
  TH1* dcaz_negA = (dynamic_cast<TH2*>(dcaz_neg3->Project3D("xz_3")))->ProfileX(); 
  dcaz_negA->Fit(&fit,"NQ");
  dcaz_negA_0 = fit.GetParameter(0);
  dcaz_negA_1 = fit.GetParameter(1);
  dcaz_negA_2 = fit.GetParameter(2);
  dcaz_negA_chi2 = fit.GetChisquare();      
  dcaz_negA_0_Err = fit.GetParError(0);
  dcaz_negA_1_Err = fit.GetParError(1);
  dcaz_negA_2_Err = fit.GetParError(2);  
  
  dcaz_neg3->GetYaxis()->SetRangeUser(-1.0,-0.01);    
  TH1* dcaz_negC = (dynamic_cast<TH2*>(dcaz_neg3->Project3D("xz_4")))->ProfileX(); 
  dcaz_negC->Fit(&fit,"NQ");
  dcaz_negC_0 = fit.GetParameter(0);
  dcaz_negC_1 = fit.GetParameter(1);
  dcaz_negC_2 = fit.GetParameter(2);
  dcaz_negC_chi2 = fit.GetChisquare();    
  dcaz_negC_0_Err = fit.GetParError(0);
  dcaz_negC_1_Err = fit.GetParError(1);
  dcaz_negC_2_Err = fit.GetParError(2);  
  //
  //
  //
  TH2 *hisTemp2D=0;
  TH1 *hisTemp1D=0;
  TH3 * hisEtaInput[4]={dcar_pos3, dcar_neg3, dcaz_pos3, dcaz_neg3};
  
  for (Int_t igr=0;igr<4; igr++){
    hisTemp2D = (dynamic_cast<TH2*>(hisEtaInput[igr]->Project3D(TString::Format("xy_%d",igr).Data())));
    hisTemp2D->GetXaxis()->SetRangeUser(-1,1);
    if (graphEtaProfile[igr]) delete graphEtaProfile[igr];
    graphEtaProfile[igr]=new TGraphErrors(hisTemp2D->ProfileX());
    delete hisTemp2D;
  }

  for (Int_t igr=0;igr<4; igr++){  
    // A side
    hisEtaInput[igr]->GetYaxis()->SetRangeUser(0,0.99);
    hisTemp2D = (dynamic_cast<TH2*>(hisEtaInput[igr]->Project3D(TString::Format("xz_%dAside",igr).Data())));
    if (graphPhiProfile[igr]) delete graphPhiProfile[igr];
    graphPhiProfile[igr]=new TGraphErrors(hisTemp2D->ProfileX());
    delete hisTemp2D;
    //
    // C side
    hisEtaInput[igr]->GetYaxis()->SetRangeUser(-1.00,-0.01);
    hisTemp2D = (dynamic_cast<TH2*>(hisEtaInput[igr]->Project3D(TString::Format("xz_%dCside",igr).Data())));
    if (graphPhiProfile[4+igr]) delete graphPhiProfile[4+igr];
    graphPhiProfile[4+igr]=new TGraphErrors(hisTemp2D->ProfileX());
    delete hisTemp2D;
    //
  }



// store results (shift in dca) in ttree
  (*pcstream)<<"tpcQA"<<      
    "grdcar_pos_Eta.="<<graphEtaProfile[0]<<
    "grdcar_neg_Eta.="<<graphEtaProfile[1]<<
    "grdcaz_pos_Eta.="<<graphEtaProfile[2]<<
    "grdcaz_neg_Eta.="<<graphEtaProfile[3]<<
    //
    "grdcar_pos_ASidePhi.="<<graphPhiProfile[0]<<
    "grdcar_neg_ASidePhi.="<<graphPhiProfile[1]<<
    "grdcaz_pos_ASidePhi.="<<graphPhiProfile[2]<<
    "grdcaz_neg_ASidePhi.="<<graphPhiProfile[3]<<
    "grdcar_pos_CSidePhi.="<<graphPhiProfile[4]<<
    "grdcar_neg_CSidePhi.="<<graphPhiProfile[5]<<
    "grdcaz_pos_CSidePhi.="<<graphPhiProfile[6]<<
    "grdcaz_neg_CSidePhi.="<<graphPhiProfile[7];
    
    
    (*pcstream)<<"tpcQA"<<      
      "dcar_posA_0="<< dcar_posA_0<<
      "dcar_posA_1="<< dcar_posA_1<<
      "dcar_posA_2="<< dcar_posA_2<<
      "dcar_posA_chi2="<< dcar_posA_chi2<<
      "dcar_posA_0_Err="<< dcar_posA_0_Err<<
      "dcar_posA_1_Err="<< dcar_posA_1_Err<<
      "dcar_posA_2_Err="<< dcar_posA_2_Err;    
      
      (*pcstream)<<"tpcQA"<<            
      "dcaz_posA_0="<< dcaz_posA_0<<
      "dcaz_posA_1="<< dcaz_posA_1<<
      "dcaz_posA_2="<< dcaz_posA_2<<
      "dcaz_posA_chi2="<< dcaz_posA_chi2<<
      "dcaz_posA_0_Err="<< dcaz_posA_0_Err<<
      "dcaz_posA_1_Err="<< dcaz_posA_1_Err<<
      "dcaz_posA_2_Err="<< dcaz_posA_2_Err;          
      
      (*pcstream)<<"tpcQA"<<            
      "dcaz_posC_0="<< dcaz_posC_0<<
      "dcaz_posC_1="<< dcaz_posC_1<<
      "dcaz_posC_2="<< dcaz_posC_2<<
      "dcaz_posC_chi2="<< dcaz_posC_chi2<<
      "dcaz_posC_0_Err="<< dcaz_posC_0_Err<<
      "dcaz_posC_1_Err="<< dcaz_posC_1_Err<<
      "dcaz_posC_2_Err="<< dcaz_posC_2_Err;           

      (*pcstream)<<"tpcQA"<<            
      "dcar_posC_0="<< dcar_posC_0<<
      "dcar_posC_1="<< dcar_posC_1<<
      "dcar_posC_2="<< dcar_posC_2<<
      "dcar_posC_chi2="<< dcar_posC_chi2<<
      "dcar_posC_0_Err="<< dcar_posC_0_Err<<
      "dcar_posC_1_Err="<< dcar_posC_1_Err<<
      "dcar_posC_2_Err="<< dcar_posC_2_Err;           
            
      
     (*pcstream)<<"tpcQA"<<      
      "dcar_negA_0="<< dcar_negA_0<<
      "dcar_negA_1="<< dcar_negA_1<<
      "dcar_negA_2="<< dcar_negA_2<<
      "dcar_negA_chi2="<< dcar_negA_chi2<<
      "dcar_negA_0_Err="<< dcar_negA_0_Err<<
      "dcar_negA_1_Err="<< dcar_negA_1_Err<<
      "dcar_negA_2_Err="<< dcar_negA_2_Err;    
      
      (*pcstream)<<"tpcQA"<<            
      "dcaz_negA_0="<< dcaz_negA_0<<
      "dcaz_negA_1="<< dcaz_negA_1<<
      "dcaz_negA_2="<< dcaz_negA_2<<
      "dcaz_negA_chi2="<< dcaz_negA_chi2<<
      "dcaz_negA_0_Err="<< dcaz_negA_0_Err<<
      "dcaz_negA_1_Err="<< dcaz_negA_1_Err<<
      "dcaz_negA_2_Err="<< dcaz_negA_2_Err;          
      
      (*pcstream)<<"tpcQA"<<            
      "dcaz_negC_0="<< dcaz_negC_0<<
      "dcaz_negC_1="<< dcaz_negC_1<<
      "dcaz_negC_2="<< dcaz_negC_2<<
      "dcaz_negC_chi2="<< dcaz_negC_chi2<<
      "dcaz_negC_0_Err="<< dcaz_negC_0_Err<<
      "dcaz_negC_1_Err="<< dcaz_negC_1_Err<<
      "dcaz_negC_2_Err="<< dcaz_negC_2_Err;           

      (*pcstream)<<"tpcQA"<<            
      "dcar_negC_0="<< dcar_negC_0<<
      "dcar_negC_1="<< dcar_negC_1<<
      "dcar_negC_2="<< dcar_negC_2<<
      "dcar_negC_chi2="<< dcar_negC_chi2<<
      "dcar_negC_0_Err="<< dcar_negC_0_Err<<
      "dcar_negC_1_Err="<< dcar_negC_1_Err<<
      "dcar_negC_2_Err="<< dcar_negC_2_Err;                 

      
      return 0;

}

//_____________________________________________________________________________                                                                                                  
Int_t AliTPCPerformanceSummary::AnalyzeOcc(const AliPerformanceTPC* pTPC, TTreeSRedirector* const pcstream)
{
  //                                                                                                                                                                      
  //function which make trending of occupany per side and IROC-OROC                                     
  //                                                                                                                                                                            

  if (!pcstream) return 16;
  if (!pTPC) return 16;

  TH3* h3D_1=0;
  TH2* his2D=0;
  TH1* his1D=0;

  static Double_t norm=0; 
  static Double_t mean_occ_chamber=0;                                                                                                                                         
  static Double_t rms_mean_occ_chamber=0;   
  static Float_t occ_chamber=0;
  static Double_t rmsNr   = 3.0;
  static Int_t n_chamber_lowOcc = 0;  
  static Double_t minOcc= 0;  
  
  //nr of chamber within the thresholds
  static Int_t iroc_A_side =0;
  static Int_t oroc_A_side=0;
  static Int_t iroc_C_side =0;
  static Int_t oroc_C_side =0;
  
  //occupancy for each chamber, normalized to the total occupancy  
  static TVectorF meanOccArray_iroc(36);
  static TVectorF meanOccArray_oroc(36);

  if (pTPC->GetHistos()->FindObject("h_tpc_clust_0_1_2")) {
      h3D_1 = dynamic_cast<TH3*>(pTPC->GetHistos()->FindObject("h_tpc_clust_0_1_2"));
  }
  else h3D_1 = pTPC->GetTPCClustHisto()->Projection(0,1,2);
  if(!h3D_1) {
    printf("E-AliTPCPerformanceSummary::AnalyzeOcc: h_tpc_clust_0_1_2 not found");
    return 4;
  }

  //////////////////////////////////////////
  // normalization
  h3D_1->GetZaxis()->SetRangeUser(0.2,0.99); //A side
  h3D_1->GetXaxis()->SetRangeUser(0,160); //IROC + OROC
  his2D  = dynamic_cast<TH2*>(h3D_1->Project3D("xy_A_norm"));
  if(!his2D) return 4;
  his1D = his2D->ProjectionX();
  norm = his1D->Integral();
  printf("normalization:  \t%f\n",norm);
  if (norm < 0.001) norm=0.00001;
  delete his2D;
  
  //////////////////////////////////////////
  // A_side IROC
  h3D_1->GetZaxis()->SetRangeUser(0.2,0.99); //A_side
  h3D_1->GetXaxis()->SetRangeUser(0,63); //IROC    

  his2D = dynamic_cast<TH2*>(h3D_1->Project3D("xy_A_side_IROC"));
  if(!his2D) return 4;

  printf("-------------- A_IROC occupancy \t\n");

  for(Int_t i = 0; i < 18; i++) { //fill IROC A_side         
    Float_t phiLow = i*(20./360.)*(2*TMath::Pi());
    Float_t phiUp  = (i+1)*(20./360.)*(2*TMath::Pi());
    his2D->GetXaxis()->SetRangeUser(phiLow,phiUp); 
    his1D = his2D->ProjectionX();
    occ_chamber = (his1D->Integral()) / norm;    
    printf("%d occ_chamber \t%f \t phiLow phiUp  \t%f  %f \n",i, occ_chamber, phiLow, phiUp);
    meanOccArray_iroc[i]= occ_chamber;//fill array with occupancy for each chamber
    mean_occ_chamber += occ_chamber;//compute the average occupancy        
    rms_mean_occ_chamber  += occ_chamber*occ_chamber;
    delete his1D;
  }
  delete his2D;

  mean_occ_chamber /= 18; //nr of sectors                                                                                                                              
  rms_mean_occ_chamber  /= 18; //nr of sectors                                            
  
  rms_mean_occ_chamber   =  TMath::Sqrt( TMath::Abs(rms_mean_occ_chamber - (mean_occ_chamber*mean_occ_chamber)) );                                         
  minOcc    = mean_occ_chamber - rmsNr*rms_mean_occ_chamber;  

  printf("mean_occ_chamber +- rms_mean_occ_chamber \t%f\t%f \n", mean_occ_chamber, rms_mean_occ_chamber);
  printf("min Occ allowed (threshold) \t%f \n", minOcc);

  for (Int_t i = 0; i<18; i++) {
    if (meanOccArray_iroc[i] < minOcc) {n_chamber_lowOcc++;}
  }
  iroc_A_side = (18 - n_chamber_lowOcc); //nr of iroc above the thr
  printf("Nr of iroc_A_side \t%i \n \n ",iroc_A_side);

  mean_occ_chamber=0;
  rms_mean_occ_chamber=0;
  occ_chamber=0.;
  n_chamber_lowOcc=0;
  minOcc=0.;
  ////////////////////////////////////////////
  // A_side OROC
  h3D_1->GetZaxis()->SetRangeUser(0.2,0.99); //A_side
  h3D_1->GetXaxis()->SetRangeUser(64,160); //OROC    

  his2D = dynamic_cast<TH2*>(h3D_1->Project3D("xy_A_side_OROC"));
  if(!his2D) return 4;

  printf("-------------- A_OROC occupancy \t\n");

  for(Int_t i = 0; i < 18; i++) {          
    Float_t phiLow = i*(20./360.)*(2*TMath::Pi());
    Float_t phiUp  = (i+1)*(20./360.)*(2*TMath::Pi());
    his2D->GetXaxis()->SetRangeUser(phiLow,phiUp); 
    his1D = his2D->ProjectionX();
    occ_chamber = (his1D->Integral()) / norm;    
    printf("%d occ_chamber \t%f \t phiLow phiUp  \t%f  %f \n",i, occ_chamber, phiLow, phiUp);
    meanOccArray_oroc[i]= occ_chamber;//fill array with occupancy for each chamber
    mean_occ_chamber += occ_chamber;//compute the average occupancy        
    rms_mean_occ_chamber  += occ_chamber*occ_chamber;
    delete his1D;
  }
  delete his2D;

  mean_occ_chamber /= 18; //nr of sectors                                                                                                                              
  rms_mean_occ_chamber  /= 18; //nr of sectors                                            
  
  rms_mean_occ_chamber   =  TMath::Sqrt( TMath::Abs(rms_mean_occ_chamber - (mean_occ_chamber*mean_occ_chamber)) );                                         
  minOcc    = mean_occ_chamber - rmsNr*rms_mean_occ_chamber;  

  printf("mean_occ_chamber +- rms_mean_occ_chamber \t%f\t%f \n", mean_occ_chamber, rms_mean_occ_chamber);
  printf("min Occ allowed (threshold) \t%f \n", minOcc);

  for (Int_t i = 0; i<18; i++) {
    if (meanOccArray_oroc[i] < minOcc) {n_chamber_lowOcc++;}
  }
  oroc_A_side = (18 - n_chamber_lowOcc); //variable stored in the trending
  printf("Nr of oroc_A_side \t%i \n \n ",oroc_A_side);

  mean_occ_chamber=0;
  rms_mean_occ_chamber=0;
  occ_chamber=0.;
  n_chamber_lowOcc=0;
  minOcc=0.;

  ////////////////////////////////////////////////////////////////////////////////
  // C side
  //////////////////////////////////////////
  
  // normalization
  h3D_1->GetZaxis()->SetRangeUser(1.1,2.1); //C side
  h3D_1->GetXaxis()->SetRangeUser(0,160); //IROC + OROC
  his2D  = dynamic_cast<TH2*>(h3D_1->Project3D("xy_C_norm"));
  if(!his2D) return 4;
  his1D = his2D->ProjectionX();
  norm = his1D->Integral();
  printf("normalization:  \t%f\n",norm);
  if (norm < 0.001) norm=0.00001;
  delete his2D;
  
  //////////////////////////////////////////
  // C_side IROC
  h3D_1->GetZaxis()->SetRangeUser(1.1,2.1); //C_side
  h3D_1->GetXaxis()->SetRangeUser(0,63); //IROC    

  his2D = dynamic_cast<TH2*>(h3D_1->Project3D("xy_C_side_IROC"));
  if(!his2D) return 4;

  printf("-------------- C_IROC occupancy \t\n");

  for(Int_t i = 0; i < 18; i++) {          
    Float_t phiLow = i*(20./360.)*(2*TMath::Pi());
    Float_t phiUp  = (i+1)*(20./360.)*(2*TMath::Pi());
    his2D->GetXaxis()->SetRangeUser(phiLow,phiUp); 
    his1D = his2D->ProjectionX();
    occ_chamber = (his1D->Integral()) / norm;    
    printf("%d occ_chamber \t%f \t phiLow phiUp  \t%f  %f \n",i, occ_chamber, phiLow, phiUp);
    meanOccArray_iroc[18+i]= occ_chamber;//fill array with occupancy for each chamber, C side, therefore index 18+i
    mean_occ_chamber += occ_chamber;//compute the average occupancy        
    rms_mean_occ_chamber  += occ_chamber*occ_chamber;
    delete his1D;
  }
  delete his2D;

  mean_occ_chamber /= 18; //nr of sectors                                                                                                                              
  rms_mean_occ_chamber  /= 18; //nr of sectors                                            
  
  rms_mean_occ_chamber   =  TMath::Sqrt( TMath::Abs(rms_mean_occ_chamber - (mean_occ_chamber*mean_occ_chamber)) );                                         
  minOcc    = mean_occ_chamber - rmsNr*rms_mean_occ_chamber;  

  printf("mean_occ_chamber +- rms_mean_occ_chamber \t%f\t%f \n", mean_occ_chamber, rms_mean_occ_chamber);
  printf("min Occ allowed (threshold) \t%f \n", minOcc);

  for (Int_t i = 18; i<36; i++) {
    if (meanOccArray_iroc[i] < minOcc) {n_chamber_lowOcc++;}
  }
  iroc_C_side = (18 - n_chamber_lowOcc); //variable stored in the trending
  printf("Nr of iroc_C_side \t%i \n \n ",iroc_C_side);

  mean_occ_chamber=0;
  rms_mean_occ_chamber=0;
  occ_chamber=0.;
  n_chamber_lowOcc=0;
  minOcc=0.;

  ////////////////////////////////////////////
  // C_side OROC
  h3D_1->GetZaxis()->SetRangeUser(1.1,2.1); //C_side
  h3D_1->GetXaxis()->SetRangeUser(64,160); //OROC    

  his2D = dynamic_cast<TH2*>(h3D_1->Project3D("xy_C_side_OROC"));
  if(!his2D) return 4;

  printf("-------------- C_OROC occupancy \t\n");

  for(Int_t i = 0; i < 18; i++) {          
    Float_t phiLow = i*(20./360.)*(2*TMath::Pi());
    Float_t phiUp  = (i+1)*(20./360.)*(2*TMath::Pi());
    his2D->GetXaxis()->SetRangeUser(phiLow,phiUp); 
    his1D = his2D->ProjectionX();
    occ_chamber = (his1D->Integral()) / norm;    
    printf("%d occ_chamber \t%f \t phiLow phiUp  \t%f  %f \n",i, occ_chamber, phiLow, phiUp);
    meanOccArray_oroc[18+i]= occ_chamber;//fill array with occupancy for each chamber
    mean_occ_chamber += occ_chamber;//compute the average occupancy        
    rms_mean_occ_chamber  += occ_chamber*occ_chamber;
    delete his1D;
  }
  delete his2D;

  mean_occ_chamber /= 18; //nr of sectors                                                                                                                              
  rms_mean_occ_chamber  /= 18; //nr of sectors                                            
  
  rms_mean_occ_chamber   =  TMath::Sqrt( TMath::Abs(rms_mean_occ_chamber - (mean_occ_chamber*mean_occ_chamber)) );                                         
  minOcc    = mean_occ_chamber - rmsNr*rms_mean_occ_chamber;  

  printf("mean_occ_chamber +- rms_mean_occ_chamber \t%f\t%f \n", mean_occ_chamber, rms_mean_occ_chamber);
  printf("min Occ allowed (threshold) \t%f \n", minOcc);

  for (Int_t i = 18; i<36; i++) {
    if (meanOccArray_oroc[i] < minOcc) {n_chamber_lowOcc++;}
  }
  oroc_C_side = (18 - n_chamber_lowOcc); //variable stored in the trending
  printf("Nr of oroc_C_side \t%i \n \n ",oroc_C_side);

  mean_occ_chamber=0;
  rms_mean_occ_chamber=0;
  occ_chamber=0.;
  n_chamber_lowOcc=0;
  minOcc=0.;

  (*pcstream)<<"tpcQA"<<      
   "iroc_A_side="<< iroc_A_side<<
   "oroc_A_side="<< oroc_A_side<<
   "iroc_C_side="<< iroc_C_side<<
   "oroc_C_side="<< oroc_C_side<<
   //A/C side IROC 
   "TPC_Occ_IROC.="<< &meanOccArray_iroc<< 
   //A/C side OROC
   "TPC_Occ_OROC.="<< &meanOccArray_oroc;   

  return 0;
}







void   AliTPCPerformanceSummary::MakeRawOCDBQAPlot(TTreeSRedirector *pcstream){
  //
  //  Draw RAW OCDB QA plot and store the summary graphs in the trending tree
  //    1.) Occupancy and cluster charge maps
  //    2.) HV information
  //
  gStyle->SetLabelSize(0.08,"XYZ");
  gStyle->SetTitleSize(0.08,"XYZ");
  gStyle->SetTitleOffset(0.5,"XYZ");
  AliTPCcalibDB *calibDB = AliTPCcalibDB::Instance();
  TVectorD   value(72), valueNorm(72),valueRMS(72), roc(72);  
  AliTPCdataQA*  dataQA= calibDB->GetDataQA();
  for (Int_t isec=0; isec<72; isec++) roc[isec]=isec;
  //
  //
  AliTPCCalPad * padActive=calibDB->GetPadGainFactor();
  AliTPCCalPad * padLocalMax=dataQA->GetNLocalMaxima();
  AliTPCCalPad * padNoThreshold=dataQA->GetNoThreshold();
  AliTPCCalPad * padMaxCharge=dataQA->GetMaxCharge();
  AliTPCCalPad * padMeanCharge=dataQA->GetMeanCharge();
  AliCDBManager *man = AliCDBManager::Instance();
  AliCDBEntry * entry=  man->Get("TPC/Calib/QA");
  static Bool_t hasRawQA = entry->GetId().GetFirstRun()== man->GetRun();


  const char   *typeNames[5]={"ActiveChannelMap", "LocaleMaxima", "NoThreshold",   "MaxCharge"  };
  AliTPCCalPad *padInput [5]={ padActive,          padLocalMax,    padNoThreshold,  padMaxCharge};
  TGraphErrors *grRaw[8]={0};
  TGraphErrors *grStatus[8]={0};
  const char * side[2]={"A","C"};
  Int_t kcolors[5]={1,2,4,3,6};
  Int_t kmarkers[5]={21,25,20,24,26};
   
  for (Int_t itype=0; itype<4; itype++){
    if (!padInput[itype]) {
      ::Error("AliTPCPerformanceSummary::MakeRawOCDBQAPlot","Could not get input for type %d: %s", itype, typeNames[itype]);
    }

    for (Int_t isec=0; isec<72; isec++) {
      value[isec]              =padInput[itype]&&padInput[itype]->GetCalROC(isec)?padInput[itype]->GetCalROC(isec)->GetMedian():-1;
      if (itype==0) value[isec]=padInput[itype]&&padInput[itype]->GetCalROC(isec)?padInput[itype]->GetCalROC(isec)->GetMean()  :-1;
      if (itype==2) value[isec]=padInput[itype]&&padInput[itype]->GetCalROC(isec)?padInput[itype]->GetCalROC(isec)->GetMedian():-1;
      valueRMS[isec]=padInput[itype]&&padInput[itype]->GetCalROC(isec)?padInput[itype]->GetCalROC(isec)->GetRMS():-1;
      if (value[isec]>0)valueRMS[isec]/=value[isec];
      if (padInput[itype]&&padInput[itype]->GetCalROC(isec)) valueRMS[isec]/=TMath::Sqrt( 16*padInput[itype]->GetCalROC(isec)->GetNchannels()/padInput[itype]->GetCalROC(isec)->GetNrows());
    }
    for (Int_t isec=0; isec<72; isec++) {
      Int_t offset=36*(isec/36);
      Double_t norm,rms;
      AliMathBase::EvaluateUni(36, &(value.GetMatrixArray()[offset]),norm,rms,33);
      valueNorm[isec]=1;
      if (norm>0) valueNorm[isec]=value[isec]/norm;  // if norm==0 means values not defined - we defineed normalized values in that case==1
    }
    grRaw[itype]= new TGraphErrors(72,roc.GetMatrixArray(),valueNorm.GetMatrixArray(),0,valueRMS.GetMatrixArray());
    grRaw[itype]->SetMarkerColor(kcolors[itype]);
    grRaw[itype]->SetMarkerStyle(kmarkers[itype]);
    grRaw[itype]->GetXaxis()->SetTitle("sector");
    grRaw[itype]->GetYaxis()->SetTitle("Norm. value");
    grRaw[itype]->SetMinimum(0); 
    grRaw[itype]->SetMaximum(2); 
  }


  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TCanvas *canvasRawQA0 = new TCanvas("canvasRawQA0","canvasRawQA0",1400,700);
  for (Int_t itype=0; itype<4; itype++){
    for (Int_t iplot=0; iplot<2; iplot++){
      canvasRawQA0->cd()->SetLogz();
      TPad *pad = new TPad("xxx","xxx",(itype+0)*0.25+0.005,(iplot+1)*0.33+0.005, (itype+1)*0.25-0.005,(iplot+2)*0.33-0.005);
      pad->Draw();
      pad->cd()->SetLogz();
      TH1* histo = 0;
      if (padInput[itype] && padInput[itype]->GetMedian()>0){
        padInput[itype]->Multiply(1./padInput[itype]->GetMedian());
	histo=padInput[itype]->MakeHisto2D((iplot+1)%2);
	if (histo->GetMaximum()>2) histo->SetMaximum(2.); 
	histo->SetName(TString::Format("%s  %s Side",padInput[itype]->GetTitle(),side[(iplot+1)%2]).Data());
	histo->SetTitle(TString::Format("%s %s Side",padInput[itype]->GetTitle(),side[(iplot+1)%2]).Data());
	histo->Draw("colz");
      } else {
        TLatex l;
        l.DrawLatexNDC(.2,.5,TString::Format("%s not available", typeNames[itype]));
      }
    }
  }
  canvasRawQA0->cd(); 
  gStyle->SetOptTitle(1);
  TPad *pad = new TPad("graph","graph",0,0,1,0.33);
  pad->SetRightMargin(0.01);
  pad->Draw();
  pad->cd();
  TLegend * legend = new TLegend(0.2,0.70,0.5,0.89,"Raw data QA (normalized to median))");
  legend->SetBorderSize(0);
  legend->SetNColumns(2);
  for (Int_t itype=0; itype<4; itype++){
    if (itype==0) grRaw[itype]->Draw("ap");
    grRaw[itype]->Draw("p");
    legend->AddEntry(grRaw[itype], typeNames[itype],"p");
  }
  legend->Draw();
  TLine *line = new TLine(0,0.70,72,0.70);
  line->SetLineStyle(2);
  line->SetLineColor(2);
  line->SetLineWidth(3);
  line->Draw();
  TLatex latex;
  latex.DrawLatex(73,0.70,"Threshold");
  canvasRawQA0->SaveAs("rawQAInformation.png");
  static Int_t clusterCounter= dataQA->GetClusterCounter();
  static Int_t signalCounter= dataQA->GetSignalCounter();
  //
  // add aso HV status
  //
  for (Int_t itype=0; itype<5; itype++){
    for (Int_t isec=0; isec<72; isec++){
      if (itype==0) valueNorm[isec]=calibDB->GetChamberHVStatus(isec);
      if (itype==1) valueNorm[isec]=calibDB->GetChamberGoodHighVoltageFraction(isec);
      if (itype==2) valueNorm[isec]=calibDB->GetChamberHighVoltageMedian(isec);
      if (itype==3) valueNorm[isec]=calibDB->GetChamberCurrentNominalHighVoltage(isec);  
      if (itype==4) valueNorm[isec]=calibDB->GetChamberCurrentNominalHighVoltage(isec)-50.;  
    }    
    grStatus[itype]= new TGraphErrors(72,roc.GetMatrixArray(),valueNorm.GetMatrixArray(),0,0);
    grStatus[itype]->SetMarkerColor(kcolors[itype]);
    grStatus[itype]->SetMarkerStyle(kmarkers[itype]);
    grStatus[itype]->GetXaxis()->SetTitle("sector");
    grStatus[itype]->GetYaxis()->SetTitle("Norm. value");
    grStatus[itype]->SetMinimum(0); 
  }
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(1);
  TCanvas *canvasROCStatusOCDB = new TCanvas("canvasROCStatusOCDB","canvasROCStatusOCDB",1000,500);
  canvasROCStatusOCDB->Divide(1,2);
  canvasROCStatusOCDB->SetRightMargin(0.02);
  TLegend * legendStatus = new TLegend(0.2,0.20,0.5,0.4,"Chamber status");
  TLegend * legendVoltage = new TLegend(0.2,0.20,0.5,0.4,"Chamber voltage");
  legendStatus->SetBorderSize(0);
  legendVoltage->SetBorderSize(0);
  {
    canvasROCStatusOCDB->cd(1)->SetRightMargin(0.02); 
    grStatus[0]->SetMinimum(0);
    grStatus[0]->Draw("ap");
    grStatus[1]->Draw("p");
    legendStatus->AddEntry(grStatus[0],"Chamber status","p");
    legendStatus->AddEntry(grStatus[1],"Time fraction in nominal condition","p");
    legendStatus->Draw();
    //
    canvasROCStatusOCDB->cd(2)->SetRightMargin(0.02);
    grStatus[3]->SetMinimum(0);
    grStatus[3]->Draw("ap");
    grStatus[2]->Draw("p");
    grStatus[4]->SetMarkerColor(2);
    grStatus[4]->SetMarkerStyle(kmarkers[2]);
    grStatus[4]->SetMarkerSize(0.75);
    grStatus[4]->SetLineColor(2);
    grStatus[4]->SetLineWidth(3);
    grStatus[4]->Draw("lp");
    legendVoltage->AddEntry(grStatus[2],"Chamber voltage","p");
    legendVoltage->AddEntry(grStatus[3],"Chamber nominal voltage (-TPC median)","p");
    legendVoltage->AddEntry(grStatus[4],"Chamber voltage Threshold","p");
    legendVoltage->Draw();
  }
  canvasROCStatusOCDB->SaveAs("canvasROCStatusOCDB.png");

  if (pcstream){
    (*pcstream)<<"tpcQA"<<
      "hasRawQA="<<hasRawQA<<                   // flag - Raw QA present
      "rawClusterCounter="<<clusterCounter<<    // absolute number of cluster  in Raw QA          -  calibDB->GetDataQA()->GetClusterCounter();
      "rawSignalCounter="<<signalCounter<<      // absolute number of signal above Thr  in Raw QA -  calibDB->GetDataQA()->GetSignalCounter()
      "grOCDBStatus.="<<grRaw[0]<<              // OCDB status as used in Reco         - LTM:calibDB->GetPadGainFactor();
      //
      "grRawLocalMax.="<<grRaw[1]<<             // RAW QA OCDB local cluster counter   - LTM:calibDB->GetDataQA()->GetNLocalMaxima();
      "grRawAboveThr.="<<grRaw[2]<<             // RAW QA OCDB above threshold counter - LTM:calibDB->GetDataQA->GetNoThreshold(); 
      "grRawQMax.="<<grRaw[3]<<                 // RAQ QA OCDB max charge              - LTM:calibDB->GetDataQA()->GetMaxCharge();
      //
      "grROCHVStatus.="<<grStatus[0]<<          // ROC HV status - enable/disable chambers beacus of LOW HV
      "grROCHVTimeFraction.="<<grStatus[1]<<    // ROC HV fraction of time disabled
      "grROCHVMedian.="<<grStatus[2]<<          // ROC median voltage
      "grROCHVNominal.="<<grStatus[3];          // ROC Nominal voltage corrected for common median shift
  }
}



void  AliTPCPerformanceSummary::MakeMissingChambersAliases(TTree * tree){
  //
  // Make default aliases for the "missing chamber selection"
  //
  tree->SetAlias("sectorNclMissing70","Sum$((max(grNclSectorNegA.fY,grNclSectorPosA.fY)/grNclPhiMedian.fElements[0])<0.7)+Sum$((max(grNclSectorNegC.fY,grNclSectorPosC.fY)/grNclPhiMedian.fElements[0])<0.7)");  
  // Missing chambers according StandardQA - using Ncl:phi  histogram.  
  //  Counter:  sum( Ncl_sector/<Ncl>_median < 70% ) 
  tree->SetAlias("sectorNtrMissing70","Sum$((18*grNtrSectorPosC.fY/Sum$(grNtrSectorPosC.fY))<0.7)+Sum$((18*grNtrSectorPosA.fY/Sum$(grNtrSectorPosA.fY))<0.7)");
  // Missing chambers according StandardQA - using Ncl:phi  track counter
  //  Counter sum( Ntr_sector/<Ntr>_median < 70% ) 
  tree->SetAlias("ocdbStatusCounter","hasRawQA*Sum$(grOCDBStatus.fY<0.75)");               // counter status for the OCDB
  tree->SetAlias("rawLowQMaxCounter75","Sum$((grRawQMax.fY)<0.75)");                       // counter small gain based on RAW QA  Qmax 
  tree->SetAlias("rawLowOccupancyCounter75","Sum$((grRawLocalMax.fY)<0.75)");              // counter small occupancy (either gain or fraction of time)
  tree->SetAlias("rawLowCounter75","Sum$(((grRawQMax.fY)<0.75||(grRawLocalMax.fY)<0.75))*(-1+2*(grRawQMax.fEY<0.03))");
  // counter outliers 75 (Qmax or occupancy based)  
  // for laser data Q is not reliable - disable decission for data with big RMS 
  tree->SetAlias("ocdbHVStatusCounter","Sum$(grROCHVStatus.fY<0.5)");                      // chambers not active  - according  HV decission
  tree->SetAlias("ocdbLowerHVCounter50","Sum$(grROCHVMedian.fY<(grROCHVNominal.fY-50))");    // chambers not active  - according  maximal diff of HV
  tree->SetAlias("qaClOccupancyCounter60","(Sum$((36*TPC_Occ_IROC.fElements)<0.60)+Sum$((36*TPC_Occ_OROC.fElements)<0.60))");   // missing chambers according cluster occupancy

  //
  // dEdx outlier alaiases
  //
  tree->SetAlias("normdEdxSector","(36*meanMIPvsSector.fElements/Sum$(meanMIPvsSector.fElements))");  
  // dEdx_{sector}/<dEdx>
  tree->SetAlias("wrongdEdxSectorCounter5","Sum$(abs(normdEdxSector-1)>0.05)&&rawLowCounter75==0&&sectorNtrMissing70<36");  
  tree->SetAlias("wrongdEdxSectorCounter2","Sum$(abs(normdEdxSector-1)>0.02)&&rawLowCounter75==0&&sectorNtrMissing70<36");  
  // dEdx_{sector}/<dEdx>
  //
  //
  tree->SetAlias("disabledGoodChambers","(sectorNclMissing70>0||sectorNtrMissing70>0)&&rawLowCounter75==0&&sectorNtrMissing70<36");
  //
  // disabled good chambers counter
  // 
  tree->SetAlias("emptyQA","(sectorNtrMissing70==36)");
  // empty QA counter

}
