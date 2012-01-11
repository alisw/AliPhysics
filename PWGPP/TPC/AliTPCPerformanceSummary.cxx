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

#include "AliGRPObject.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibDButil.h"
#include "TTreeStream.h"
#include "AliPerformanceTPC.h"
#include "AliPerformanceDEdx.h"
#include "AliPerformanceDCA.h"
#include "AliPerformanceMatch.h"

#include "AliTPCPerformanceSummary.h"

ClassImp(AliTPCPerformanceSummary)

Bool_t AliTPCPerformanceSummary::fgForceTHnSparse = kFALSE;


//_____________________________________________________________________________
void AliTPCPerformanceSummary::WriteToTTreeSRedirector(const AliPerformanceTPC* pTPC, const AliPerformanceDEdx* pTPCgain, const AliPerformanceMatch* pTPCMatch,const AliPerformanceMatch* pTPCPull, const AliPerformanceMatch* pConstrain, TTreeSRedirector* const pcstream, Int_t run)
{
   // 
    // Extracts performance parameters from pTPC and pTPCgain.
    // Output is written to pcstream.
    // The run number must be provided since it is not stored in 
    // AliPerformanceTPC or AliPerformanceDEdx.
    //
    if (run <= 0 ) {
        if (pTPCMatch) {run = pTPCMatch->GetRunNumber(); }
        if (pTPCgain) {run = pTPCgain->GetRunNumber(); }
        if (pTPC) { run = pTPC->GetRunNumber(); }
    }
    TObjString runType;

    //AliTPCcalibDB     *calibDB=0;

//     AliTPCcalibDButil *dbutil =0;
    Int_t startTimeGRP=0;
    Int_t stopTimeGRP=0;   
    Int_t time=0;
    Int_t duration=0;

    //Float_t currentL3 =0;
    //Int_t polarityL3 = 0;
    //Float_t bz = 0;

    //calibDB = AliTPCcalibDB::Instance();

//     dbutil= new AliTPCcalibDButil;   
        
    //printf("Processing run %d ...\n",run);
    //if (calibDB) { 
    //AliTPCcalibDB::Instance()->SetRun(run); 

//     dbutil->UpdateFromCalibDB();
//     dbutil->SetReferenceRun(run);
//     dbutil->UpdateRefDataFromOCDB();     
     
    //if (calibDB->GetGRP(run)){
    //startTimeGRP = AliTPCcalibDB::GetGRP(run)->GetTimeStart();
    //stopTimeGRP  = AliTPCcalibDB::GetGRP(run)->GetTimeEnd();
    //currentL3 = AliTPCcalibDB::GetL3Current(run);
    //polarityL3 = AliTPCcalibDB::GetL3Polarity(run);
    //bz = AliTPCcalibDB::GetBz(run);
    
    //}    
    //runType = AliTPCcalibDB::GetRunType(run).Data();  
    //}  
  time = (startTimeGRP+stopTimeGRP)/2;
  duration = (stopTimeGRP-startTimeGRP);
    
    if (!pcstream) return;
    (*pcstream)<<"tpcQA"<<      
      "run="<<run<<
      "time="<<time<<
      "startTimeGRP="<<startTimeGRP<<
      "stopTimeGRP="<<stopTimeGRP<<
      "duration="<<
      "runType.="<<&runType;
    if (pTPC) {
        pTPC->GetTPCTrackHisto()->GetAxis(9)->SetRangeUser(0.5,1.5);
        pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0.25,10);
        pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-1,1);    
        AnalyzeNCL(pTPC, pcstream);    
        AnalyzeDrift(pTPC, pcstream);
        AnalyzeDriftPos(pTPC, pcstream);
        AnalyzeDriftNeg(pTPC, pcstream);    
        AnalyzeDCARPhi(pTPC, pcstream);
        AnalyzeDCARPhiPos(pTPC, pcstream);
        AnalyzeDCARPhiNeg(pTPC, pcstream);
        AnalyzeEvent(pTPC, pcstream);         

	AnalyzePt(pTPC,pcstream);
	AnalyzeChargeOverPt(pTPC,pcstream); 
	

        pTPC->GetTPCTrackHisto()->GetAxis(9)->SetRangeUser(-10,10);
        pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0,100);
        pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-10,10); 
    }
    AnalyzeGain(pTPCgain, pcstream);
    AnalyzeMatch(pTPCMatch, pcstream);
    AnalyzePull(pTPCPull, pcstream);
    AnalyzeConstrain(pConstrain, pcstream);
   
    (*pcstream)<<"tpcQA"<<"\n";
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
    pcstream = new TTreeSRedirector(outfile);
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
    SaveGraph(tree,"meanMultNeg","run",condition);
    SaveGraph(tree,"rmsMultNeg","run",condition);
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
    his2D->FitSlicesY(0,0,-1,0,"QNR",&arrayWidth);
    width =  dynamic_cast<TH1*>(arrayWidth.At(2));
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

    linearFit.ClearPoints();
    
    //get his2D in C Side
    his3D->GetYaxis()->SetRangeUser(-1,-0.001);
    his2D  = dynamic_cast<TH2*>(his3D->Project3D("xz"));
    his2D->FitSlicesY(0,0,-1,0,"QNR",&arrayWidth);
    width =  dynamic_cast<TH1*>(arrayWidth.At(2));
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
    //
    
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
    pTPC->GetTPCTrackHisto()->GetAxis(5)->SetRangeUser(-1.,1.);
    pTPC->GetTPCTrackHisto()->GetAxis(7)->SetRangeUser(0.25,10);
    
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
        
    pTPC->GetTPCTrackHisto()->GetAxis(2)->SetRangeUser(0,10);
    
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
    static TVectorD sector(36);
    static Float_t meanMIP = 0;
    static Float_t resolutionMIP = 0;
    static Float_t attachSlopeC = 0;
    static Float_t attachSlopeA = 0;

    TH1 * his1D = 0;
    //TH1 * hisProj1D=0;
    TH2* his2D=0;
     

    meanMIPvsSector.Zero();
    //
    // select MIP particles
    //
    pTPCgain->GetDeDxHisto()->GetAxis(7)->SetRangeUser(0.4,0.55);
    pTPCgain->GetDeDxHisto()->GetAxis(0)->SetRangeUser(35,60);
    pTPCgain->GetDeDxHisto()->GetAxis(6)->SetRangeUser(80,160);
    pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-1,1);
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
    pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-3,0); // C side
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
    pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(0,3); // A side
    if (pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_5") && !fgForceTHnSparse) {    
        his2D = dynamic_cast<TH2*>(pTPCgain->GetHistos()->FindObject("h_tpc_dedx_mips_a_0_5")->Clone());
    } else {
        his2D =  pTPCgain->GetDeDxHisto()->Projection(0,5);
    }         
    if(!his2D) return 4;

    TF1 * fpolA = new TF1("fpolA","pol1");
    TObjArray arrayFitA;
    his2D->FitSlicesY(0,0,-1,10,"QN",&arrayFit);    
    his1D = (TH1*) arrayFit.At(1);
    his1D->Fit(fpolA,"QNROB=0.8","QN",0,1);
    attachSlopeA = fpolA->GetParameter(1);
     //removedtotest// delete his2D;
     //removedtotest// delete his1D;
    //
    // MIP position vs. sector
    //
    pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(-3,0); // C side
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
        meanMIPvsSector(i) = gausFunc.GetParameter(1);
        sector(i)=i;
        //removedtotest// delete his1D;
    }
     //removedtotest// delete his2D;
    //
    pTPCgain->GetDeDxHisto()->GetAxis(5)->SetRangeUser(0,3); // A side
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
        meanMIPvsSector(i+18) = gausFunc.GetParameter(1);
        sector(i+18)=i+18;
        //removedtotest// delete his1D;
    }
     //removedtotest// delete his2D;
    //
    printf("Gain QA report\n");
    printf("MIP mean\t%f\n",meanMIP);
    printf("MIP resolution\t%f\n",resolutionMIP);
    printf("MIPslopeA\t%f\n",attachSlopeA);
    printf("MIPslopeC\t%f\n",attachSlopeC);
    // 
    
    (*pcstream)<<"tpcQA"<<
        "MIPattachSlopeC="<<attachSlopeC<<
        "MIPattachSlopeA="<<attachSlopeA<<
        "resolutionMIP="<<resolutionMIP<<
        "meanMIPvsSector.="<<&meanMIPvsSector<<
        "sector.="<<&sector<<
        "meanMIP="<<meanMIP;

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
    static Double_t meanVertX=0;
    static Double_t rmsVertX=0;
    static Double_t meanVertY=0;
    static Double_t rmsVertY=0;
    static Double_t meanVertZ=0;
    static Double_t rmsVertZ=0;
    static Double_t vertStatus=0;
    static Double_t meanMult=0;
    static Double_t rmsMult=0;
    static Double_t meanMultPos=0;
    static Double_t rmsMultPos=0;
    static Double_t meanMultNeg=0;
    static Double_t rmsMultNeg=0;
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
    
    pTPC->GetTPCEventHisto()->GetAxis(6)->SetRange(2,2);
   
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_0") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_0")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(0);
    }
    if(!his1D) return 1;

    meanVertX = his1D->GetMean();    
    rmsVertX    = his1D->GetRMS();
    delete his1D;
    
    
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_1") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_1")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(1);
    }
    if(!his1D) return 1;

    meanVertY = his1D->GetMean();
    rmsVertY    = his1D->GetRMS();
    delete his1D;
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_2") && !fgForceTHnSparse) {    
        hc = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_2"));
	if(!hc) return 1;
        //his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_2")->Clone());
        his1D = (TH1*)hc->Clone();
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(2);
    }    
    if(!his1D) return 1;

    meanVertZ = his1D->GetMean();
    rmsVertZ    = his1D->GetRMS();
    delete his1D;
    
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_3") && !fgForceTHnSparse) {    
        hc = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_3"));
	if(!hc) return 1;
        //his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_3")->Clone());
        his1D = (TH1*)hc->Clone();
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(3);
    }
    if(!his1D) return 1;

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
    delete his1D;
    
    if (pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_5") && !fgForceTHnSparse) {    
        his1D = dynamic_cast<TH1*>(pTPC->GetHistos()->FindObject("h_tpc_event_recvertex_5")->Clone());
    } else {
       his1D = pTPC->GetTPCEventHisto()->Projection(5);
    }
    if(!his1D) return 1;

    meanMultNeg    = his1D->GetMean();
    rmsMultNeg     = his1D->GetRMS();
    delete his1D;
    
    pTPC->GetTPCEventHisto()->GetAxis(6)->SetRange(1,2);
    //
    (*pcstream)<<"tpcQA"<<
        "meanVertX="<<meanVertX<<
        "rmsVertX="<<rmsVertX<<
        "meanVertY="<<meanVertY<<
        "rmsVertY="<<rmsVertY<<
        "meanVertZ="<<meanVertZ<<
        "rmsVertZ="<<rmsVertZ<<
        "vertStatus="<<vertStatus<<
        "vertAll="<<vertAll<<
        "vertOK="<<vertOK<<
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
