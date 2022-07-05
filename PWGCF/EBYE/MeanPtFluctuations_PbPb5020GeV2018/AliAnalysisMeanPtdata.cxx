
/**************************************************************************
 Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *Class designed for Mean pt fluctuation Analysis in 5 TeV.
  Author: (Tulika Tripathy, Sadhana Dash),   IIT Bombay                                                                      *                                    *
 * Contributors are mentioned in the code where appropriate.              *
 **************************************************************************/



#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TList.h"
#include "TPDGCode.h"
#include "THnSparse.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliMultSelection.h"
#include "AliVTrack.h"

#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliAODPid.h"


#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"

#include "AliInputEventHandler.h"
#include "AliAODpidUtil.h"
#include "AliAnalysisMeanPtdata.h"
#include "TProfile.h"
#include "AliEventCuts.h"
#include "TProfile2D.h"

ClassImp(AliAnalysisMeanPtdata)

//________________________________________________________________________
AliAnalysisMeanPtdata::AliAnalysisMeanPtdata() // All data members should be initialised here
:AliAnalysisTaskSE(),
  fOutput(0),
  fPIDResponse(0),
  fAOD(0), 
  fPrimaryVtx(0),
  fEvents(0),
  fVtxZ(0),
  fVtxZDiff1(0),
  fVtxZDiff2(0),
  fVtxZDiff3(0),
  fVtxZCut(0),
  fVtxZCont(0),
  fVtxZCutDiff(0),
  fHistPt(0),
  fHistEta(0),
  fHistCentralityMultSelection(0),
//  _usePileupCut_PbPb5TeV(0),
  _hProfPileupCut(0),
 _profV0MvsTPCout(0),
  _fb              (0),
  _nClusterMin          ( 80),
  _dcaZMin              ( -3.2),
  _dcaZMax              (  3.2),
 _dcaXYMin             ( -2.4),
 _dcaXYMax             (  2.4),
 _chi2perTPC             (36),
 _chi2perITS             (36),
 _nTPCCrossRows (120),
  _pileUpEvent(0),
  _vZMin(0),
  _vZMax(10),
  _fhV0MvsTracksTPCout_before(0),
  _fhV0MvsTracksTPCout_after(0),
 _fV0MmultVsSpdTracklet_before(0),
 _fV0MmultVsSpdTracklet_after(0),
  _singlesOnly   ( 0),
//  fCorrDet_beforeCut( 0 ),
//  fCorrDet_afterCut( 0 ),
  fCount(0),
  fMultMeanQ(0),
  fMultTwopart(0),
  fMultThreepart(0),
  fMultFourpart(0),
  fMultQ1(0),
  fMultTwopartA(0),
  fMultThreepartA(0),
  fMultFourpartA(0),
  fMultA(0),
  fMultPairsA(0),
  fMultTripletsA(0),
  fMultQuadsA(0),
  fcentEvents(0),
 fcentEvents5(0),
 fMeanQcent0p5(0),
 fTwopartcent0p5(0),
 fMeanQcent(0),
  fTwopartcent(0),
  fMeanQcent5(0),
  fTwopartcent5(0),
  fMeanQ(0),
  fTwopart(0),
  fMeanQ10(0),
  fTwopart10(0),
  fMeanQ20(0),
  fTwopart20(0),
  fMeanQ50(0),
  fTwopart50(0),
  fMeanQ100(0),
  fTwopart100(0),
  fThreepart(0),
  fFourpart(0),
  fMeanQ_ss(0),
  fTwopart_ss(0),
  fThreepart_ss(0),
  fFourpart_ss(0),
  fMeanQcent_ss0p5(0),
  fTwopartcent_ss0p5(0),				     
  fMeanQcent_ss(0),                                                
  fTwopartcent_ss(0),                                              
  fMeanQcent_ss5(0),                                         
  fTwopartcent_ss5(0),
  fQ1A(0),
  fTwopartA(0),
  fThreepartA(0),
  fFourpartA(0),
  fA(0),
  fPairsA(0),
  fTripletsA(0),
  fQuadsA(0),
  fAliEventCuts(0x0),
  fcentnpart(0),
  h2d_tpc_ITSl1_before(0),
  h2d_tpc_ITSl2_before(0),
  h2d_ITSl1_ITSl2_before(0),
  h2d_tpc_ITSl1_after(0),
  h2d_tpc_ITSl2_after(0),
  h2d_ITSl1_ITSl2_after(0),
  eventcount(0),
  fEventSee(0),
  fNContributors(0),
  fMaxVertexZDiff1(0),
  fcentnpart2d(0),
  fV0MtoTrkTPCout_lower_PbPb5TeV(NULL),
  fV0MtoTrkTPCout_upper_PbPb5TeV(NULL),
  fNoOfTPCoutTracks(0),
  fV0Multiplicity(0),
  fV0AMultiplicity(0),
  fV0CMultiplicity(0),
  ftrack(0),
  fHdcaz(0),
  fHdcaxy(0),
  fTpcNCrossedRows(0),
  fTpcNCluster(0),
  fcentnpart2d_1(0)

  //  fEventCuts(0) 

{
    // Dummy constructor ALW

}

//________________________________________________________________________
AliAnalysisMeanPtdata::AliAnalysisMeanPtdata(const char *name) // All data members should be initialised here
:AliAnalysisTaskSE(name), 
fOutput(0),
fPIDResponse(0),
fAOD(0), 
fPrimaryVtx(0),
fEvents(0),
fVtxZ(0),
 fVtxZDiff1(0),
 fVtxZDiff2(0),
 fVtxZDiff3(0),
 fVtxZCut(0),
 fVtxZCont(0),
 fVtxZCutDiff(0),
 // _usePileupCut_PbPb5TeV(0),
 _hProfPileupCut(0),
_profV0MvsTPCout(0),
 _fb              (0),
 _nClusterMin          ( 80),
 _dcaZMin              ( -3.2),
 _dcaZMax              (  3.2),
 _dcaXYMin             ( -2.4),
 _dcaXYMax             (  2.4),
 _chi2perTPC             (36),
 _chi2perITS             (36),
 _nTPCCrossRows (120),
  _pileUpEvent(0),
 _vZMin(0),
 _vZMax(10),
 _fhV0MvsTracksTPCout_before(0),
 _fhV0MvsTracksTPCout_after(0),
 _fV0MmultVsSpdTracklet_before(0),
 _fV0MmultVsSpdTracklet_after(0),
 _singlesOnly   ( 0),
// fCorrDet_beforeCut( 0 ),
//  fCorrDet_afterCut( 0 ),
fHistPt(0),
 fHistEta(0),
  fHistCentralityMultSelection(0),
  fCount(0),
  fMultMeanQ(0),
  fMultTwopart(0),
  fMultThreepart(0),
  fMultFourpart(0),
  fMultQ1(0),
  fMultTwopartA(0),
  fMultThreepartA(0),
  fMultFourpartA(0),
  fMultA(0),
  fMultPairsA(0),
  fMultTripletsA(0),
  fMultQuadsA(0),
 fcentEvents(0),
 fcentEvents5(0),
 fMeanQcent0p5(0),
 fTwopartcent0p5(0),
 fMeanQcent(0),
  fTwopartcent(0),
  fMeanQcent5(0),
  fTwopartcent5(0),
  fMeanQ(0),
  fTwopart(0),
 fMeanQ10(0),
 fTwopart10(0),
 fMeanQ20(0),
 fTwopart20(0),
  fMeanQ50(0),
  fTwopart50(0),
  fMeanQ100(0),
  fTwopart100(0),
  fThreepart(0),
  fFourpart(0),
  fMeanQ_ss(0),
  fTwopart_ss(0),
  fThreepart_ss(0),
  fFourpart_ss(0),
  fMeanQcent_ss0p5(0),
  fTwopartcent_ss0p5(0),
  fMeanQcent_ss(0),
  fTwopartcent_ss(0),
  fMeanQcent_ss5(0),
  fTwopartcent_ss5(0),
  fQ1A(0),
  fTwopartA(0),
  fThreepartA(0),
  fFourpartA(0),
  fA(0),
  fPairsA(0),
  fTripletsA(0),
  fQuadsA(0),
  fAliEventCuts(0x0),
  fcentnpart(0),
  h2d_tpc_ITSl1_before(0),
  h2d_tpc_ITSl2_before(0),
  h2d_ITSl1_ITSl2_before(0),
  h2d_tpc_ITSl1_after(0),
  h2d_tpc_ITSl2_after(0),
  h2d_ITSl1_ITSl2_after(0),
  eventcount(0),
  fEventSee(0),
  fMaxVertexZDiff1(0),
  fNContributors(0),
  fcentnpart2d(0),
  fV0MtoTrkTPCout_lower_PbPb5TeV(NULL),
  fV0MtoTrkTPCout_upper_PbPb5TeV(NULL),
  fNoOfTPCoutTracks(0),
  fV0Multiplicity(0),
  fV0AMultiplicity(0),
  fV0CMultiplicity(0),
  ftrack(0),
  fHdcaz(0),
  fHdcaxy(0),
  fTpcNCrossedRows(0),
  fTpcNCluster(0),
  fcentnpart2d_1(0)

  //  fEventCuts(0)
  
{
  fAliEventCuts = new AliEventCuts();
DefineOutput(1, TList::Class());                            // for output list
}

//________________________________________________________________________
AliAnalysisMeanPtdata::~AliAnalysisMeanPtdata()
{

  if (fOutput) delete fOutput;    
    
}

//________________________________________________________________________
void AliAnalysisMeanPtdata::UserCreateOutputObjects()
{
fOutput = new TList();
fOutput->SetOwner(true);
fAliEventCuts->AddQAplotsToList(fOutput);   
 

 fEvents = new TH1D("fEvents","Check cuts",15,0,15);


  fVtxZ = new TH1D("fVtxZ","Vertex Z distribution before cuts",100,-20,20);
  fVtxZCut = new TH1F("fVtxZCut","Vertex Z distribution after vtxZ cut",110,-11,11);
  fVtxZCont = new TH1F("fVtxZCont","Vertex Z distribution after nCont cut",110,-11,11);
  fVtxZCutDiff = new TH1F("fVtxZCutDiff","Vertex Z distribution after cut on vtx Z Diff",110,-11,11);

  fVtxZDiff1 = new TH1F("fVtxZDiff1","Difference 1 between vertex Z distributions",100,-5,5);
  fVtxZDiff2 = new TH1F("fVtxZDiff2","Difference 2 between vertex Z distributions",100,-5,5);
  fVtxZDiff3 = new TH1F("fVtxZDiff3","Difference 3 between vertex Z distributions",100,-5,5);
  
        
  Int_t ptbins = 500;
  Float_t ptlow = 0.1, ptup = 5.0;
  fHistPt = new TH1D("fHistPt", "P_{T} distribution ", ptbins, ptlow, ptup);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  Int_t etabins = 40;
  Float_t etalow = -1.0, etaup = 1.0;
  fHistEta = new TH1D("fHistEta","#eta distribution",etabins, etalow, etaup);
  fHistEta->GetXaxis()->SetTitle("#eta");
  fHistEta->GetYaxis()->SetTitle("counts");


  fHistCentralityMultSelection = new TH1D("fHistCentralityMultSelection","Centrality Percentile ;Centrality;Entries",100,0.,100.);  //new sk
  fHistCentralityMultSelection->GetXaxis()->SetTitle("Centrality (%)");
  fHistCentralityMultSelection->GetYaxis()->SetTitle("counts");
  fHistCentralityMultSelection->SetMarkerStyle(kFullCircle);  


  fCount = new TH1D("fCount"," number of events",50,0.0,5000.0);

  fMultMeanQ = new TH1D("fMultMeanQ"," meanpt",50,0.0,5000.0);
  fMultTwopart = new TH1D("fMultTwopart"," ",50,0.0,5000.0);
  fMultThreepart = new TH1D("fMultThreepart"," ",50,0.0,5000.0);
  fMultFourpart = new TH1D("fMultFourpart"," ",50,0.0,5000.0);


  fMultQ1 = new TH1D("fMultQ1"," meanpt",50,0.0,5000.0);
  fMultTwopartA = new TH1D("fMultTwopartA"," ",50,0.0,5000.0);
  fMultThreepartA = new TH1D("fMultThreepartA"," ",50,0.0,5000.0);
  fMultFourpartA = new TH1D("fMultFourpartA"," ",50,0.0,5000.0);


  fMultA = new TH1D("fMultA"," meanpt",50,0.0,5000.0);
  fMultPairsA = new TH1D("fMultPairsA"," ",50,0.0,5000.0);
  fMultTripletsA = new TH1D("fMultTripletsA"," ",50,0.0,5000.0);
  fMultQuadsA = new TH1D("fMultQuadsA"," ",50,0.0,5000.0);


  fMeanQ100 = new TProfile("fMeanQ100"," meanpt",50,0.0,5000.0);
  fTwopart100 = new TProfile("fTwopart100"," ",50,0.0,5000.0);

 fMeanQ20 = new TProfile("fMeanQ20"," meanpt",250,0.0,5000.0);
  fTwopart20 = new TProfile("fTwopart20"," ",250,0.0,5000.0);

 fMeanQ10 = new TProfile("fMeanQ10"," meanpt",500,0.0,5000.0);
 fTwopart10 = new TProfile("fTwopart10"," ",500,0.0,5000.0);

  fMeanQ50 = new TProfile("fMeanQ50"," meanpt",100,0.0,5000.0);
  fTwopart50 = new TProfile("fTwopart50"," ",100,0.0,5000.0);


  fcentEvents = new TH1D("fcentEvents"," meanpt",100,0.0,100.0);
  fcentEvents5 = new TH1D("fcentEvents5"," meanpt",20,0.0,100.0);

  //        const Int_t NBINS = 35;
  //	Double_t edgestwo[NBINS + 1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0,70.0,90.0,120.0,140.0,160.0,180.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1200.0,1400.0,1600.0,1800.0,2000.0,2200.0,2400.0,2600.0,2800.0,3000.0,3400.0,3800.0,4200.0,4600.0,5000.0};

  //26+2+4+8+10+5
  //        const Int_t NBINS = 58;
  // Double_t edgestwo[NBINS + 1] = {0.0, 2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0, 20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,45.0,50.0, 55.0,60.0,65.0,70.0,75.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1200.0,1400.0,1600.0,1800.0,2000.0,2200.0,2400.0,2600.0,2800.0,3000.0,3400.0,3800.0,4200.0,4600.0,5000.0};

  const Int_t NBINS =52;
  Double_t edgestwo[NBINS + 1] = {0.0, 2.0,4.0,6.0,8.0,10.0,12.0,14.0,16.0,18.0, 20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,38.0,40.0,45.0,50.0, 55.0,60.0,65.0,70.0,75.0,80.0,90.0,100.0,120.0,140.0,160.0,180.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0,900.0,1000.0,1200.0,1400.0,1600.0,1800,2400.0,3000.0,3600.0,4600.0,5000.0};
     
   fMeanQ = new TProfile("fMeanQ","Two-particle correlator",/* number of bins */ NBINS,/* edge array */edgestwo );
   fTwopart = new TProfile("fTwopart","Two-particle correlator",/* number of bins */ NBINS,/* edge array */edgestwo );

   fMeanQcent0p5=new TProfile("fMeanQcent0p5"," meanpt",200,0.0,100.0);
   fTwopartcent0p5=new TProfile("fTwopartcent0p5"," meanpt",200,0.0,100.0);
   fMeanQcent=new TProfile("fMeanQcent"," meanpt",100,0.0,100.0);
   fTwopartcent=new TProfile("fTwopartcent"," meanpt",100,0.0,100.0);
   fMeanQcent5=new TProfile("fMeanQcent5"," meanpt",20,0.0,100.0);
   fTwopartcent5=new TProfile("fTwopartcent5"," meanpt",20,0.0,100.0);
   

   fThreepart = new TProfile("fThreepart"," ", NBINS,/* edge array */edgestwo );
   fFourpart = new TProfile("fFourpart"," ", NBINS,/* edge array */edgestwo );

   fMeanQ_ss = new TProfile2D("fMeanQ_ss","Check cuts",50,0,50,NBINS,/* edge array */edgestwo );
   fTwopart_ss = new TProfile2D("fTwopart_ss","Check cuts",50,0,50,NBINS,/* edge array */edgestwo );
   fThreepart_ss = new TProfile2D("fThreepart_ss","Check cuts",50,0,50,NBINS,/* edge array */edgestwo );
   fFourpart_ss = new TProfile2D("fFourpart_ss","Check cuts",50,0,50,NBINS,/* edge array */edgestwo );

   fMeanQcent_ss0p5=new TProfile2D("fMeanQcent_ss0p5","Check cuts",50,0,50,200,0.0,100.0);
   fTwopartcent_ss0p5=new TProfile2D("fTwopartcent_ss0p5","Check cuts",50,0,50,200,0.0,100.0);
   fMeanQcent_ss=new TProfile2D("fMeanQcent_ss","Check cuts",50,0,50,100,0.0,100.0);
   fTwopartcent_ss=new TProfile2D("fTwopartcent_ss","Check cuts",50,0,50,100,0.0,100.0);
   fMeanQcent_ss5=new TProfile2D("fMeanQcent_ss5","Check cuts",50,0,50,20,0.0,100.0);
   fTwopartcent_ss5=new TProfile2D("fTwopartcent_ss5","Check cuts",50,0,50,20,0.0,100.0);
   
   fEventSee = new TH1D("fEventSee","Check cuts",20,0,20);
   
   fQ1A = new TProfile("fQ1A"," meanpt",50,0.0,5000.0);
   fTwopartA = new TProfile("fTwopartA"," ",50,0.0,5000.0);
   fThreepartA = new TProfile("fThreepartA"," ",50,0.0,5000.0);
   fFourpartA = new TProfile("fFourpartA"," ",50,0.0,5000.0);
   
   fA = new TProfile("fA"," meanpt",50,0.0,5000.0);
   fPairsA = new TProfile("fPairsA"," ",50,0.0,5000.0);
   fTripletsA = new TProfile("fTripletsA"," ",50,0.0,5000.0);
   fQuadsA = new TProfile("fQuadsA"," ",50,0.0,5000.0);
   
   
   
   const Int_t NBINStwoSq = 13;
   Double_t edgestwoSq[NBINStwoSq + 1] = {0.0, 2.5, 5.0, 7.5, 10.0, 20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100.0};
     
   fcentnpart = new TH1D("fcentnpart","Two-particle correlator",/* number of bins */ NBINStwoSq,/* edge array */ edgestwoSq );
   fcentnpart2d = new TH2D("fcentnpart2d","Two-particle correlator",/* number of bins */ NBINStwoSq,/* edge array */ edgestwoSq,5000,0.,5000. );
   fcentnpart2d_1 = new TH2D("fcentnpart2d_1","Two-particle correlator",100,0.0,100.0,5000,0,5000);

   h2d_tpc_ITSl1_before=new TH2D("h2d_tpc_ITSl1_before","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);
   h2d_tpc_ITSl2_before=new TH2D("h2d_tpc_ITSl2_before","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);
   h2d_ITSl1_ITSl2_before=new TH2D("h2d_ITSl1_ITSl2_before","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);
   
   h2d_tpc_ITSl1_after=new TH2D("h2d_tpc_ITSl1_after","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);
   h2d_tpc_ITSl2_after=new TH2D("h2d_tpc_ITSl2_after","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);
   h2d_ITSl1_ITSl2_after=new TH2D("h2d_ITSl1_ITSl2_after","  ",5000,0.0,5000000.0,5000,0.0,5000000.0);

   _fV0MmultVsSpdTracklet_before =new TH2D("_fV0MmultVsSpdTracklet_before","  ",5000,0.0,50000.0,5000,0.0,50000.0);
   _fV0MmultVsSpdTracklet_before->GetXaxis()->SetTitle("SPD tracklet");
   _fV0MmultVsSpdTracklet_before->GetYaxis()->SetTitle("V0");

   _fV0MmultVsSpdTracklet_after =new TH2D("_fV0MmultVsSpdTracklet_after","  ",5000,0.0,50000.0,5000,0.0,50000.0);
   _fV0MmultVsSpdTracklet_after->GetXaxis()->SetTitle("SPD tracklet");
   _fV0MmultVsSpdTracklet_after->GetYaxis()->SetTitle("V0");

   _fhV0MvsTracksTPCout_before =new TH2D("_fhV0MvsTracksTPCout_before","  ",5000,0.0,50000.0,5000,0.0,50000.0);
   _fhV0MvsTracksTPCout_before->GetXaxis()->SetTitle("V0");
   _fhV0MvsTracksTPCout_before->GetYaxis()->SetTitle("TPC tracks");

   _fhV0MvsTracksTPCout_after =new TH2D("_fhV0MvsTracksTPCout_after","  ",5000,0.0,50000.0,5000,0.0,50000.0);
   _fhV0MvsTracksTPCout_after->GetXaxis()->SetTitle("V0");
   _fhV0MvsTracksTPCout_after->GetYaxis()->SetTitle("TPC tracks");


   _profV0MvsTPCout = new TProfile("_profV0MvsTPCout", "V0M vs TPCout profile",5000,0.0,50000.0,"s");

   ftrack= new TH1D("ftrack","Multiplicity histogram",1000,0.0,5000.0);
  
   // const Int_t ndims = 3; //number of dimensions 
     
 /*   TString sparseTitle_beforeCut[ndims] = {"TPCout_beforeCut",
					   "V0Mult_beforeCut",
					   "spdTracklet_beforeCut"};
     
     
   TString sparseTitle_afterCut[ndims] = {"TPCout_afterCut",
					  "V0Mult_afterCut",
					  "spdTracklet_afterCut"};
     
   Int_t    binVal = 2000;
   Double_t minVal = 0.0;
   Double_t maxVal = 2000.0;
     
   Int_t bins[ndims]; //number of bins in every dimension
     
     
   Double_t xmin[ndims];//minimum range
     
   Double_t xmax[ndims];//maximum range 
     
   TString addChar[ndims];
   for(Int_t iDim = 0; iDim < ndims; iDim++){
     bins[iDim] = binVal;
     xmin[iDim] = minVal;
     xmax[iDim] = maxVal;
     addChar[iDim] = "_";
   }
     
   fCorrDet_beforeCut = new THnSparseD("fCorrDet_beforeCut", "",
				       ndims, bins, xmin, xmax);
   fCorrDet_afterCut = new THnSparseD("fCorrDet_afterCut", "",
				      ndims, bins, xmin, xmax);
     
   for(Int_t iDim = 0; iDim < ndims; iDim++){
     fCorrDet_beforeCut->GetAxis(iDim)->SetTitle(sparseTitle_beforeCut[iDim]);
     fCorrDet_beforeCut->GetAxis(iDim)->SetName( sparseTitle_beforeCut[iDim]+ addChar[iDim] + std::to_string(iDim) );
         
     fCorrDet_afterCut->GetAxis(iDim)->SetTitle(sparseTitle_afterCut[iDim]);
     fCorrDet_afterCut->GetAxis(iDim)->SetName( sparseTitle_afterCut[iDim] + addChar[iDim] + std::to_string(iDim));
   }

 */

   fHdcaxy = new TH1F("fHdcaxy","fHdcaxy",1000,-3,3);

   fHdcaz = new TH1F("fHdcaz","fHdcaz",1000,-3,3);

   fTpcNCrossedRows = new TH1D("fTpcNCrossedRows","fTpcNCrossedRows",200,0,200);
   fTpcNCluster = new TH1D("fTpcNCluster","fTpcNCluster",200,0,200);




    //Fill histos here

   fOutput->Add(fVtxZ);
   fOutput->Add(fVtxZCut);
   fOutput->Add(fVtxZCont);
   fOutput->Add(fVtxZCutDiff);
   fOutput->Add(fVtxZDiff1);
   fOutput->Add(fVtxZDiff2);
   fOutput->Add(fVtxZDiff3);
   fOutput->Add(fEvents);
   fOutput->Add(fcentnpart);
   fOutput->Add(fcentnpart2d);
   fOutput->Add(fcentnpart2d_1);
   fOutput->Add(fHistPt);
   fOutput->Add(fHistEta);
   fOutput->Add(fHistCentralityMultSelection);
   fOutput->Add(fCount);
   fOutput->Add(fMultMeanQ);
   fOutput->Add(fMultTwopart);
   fOutput->Add(fMultThreepart);
   fOutput->Add(fMultFourpart);
   fOutput->Add(fMultQ1);
   fOutput->Add(fMultTwopartA);
   fOutput->Add(fMultThreepartA);
   fOutput->Add(fMultFourpartA);

    fOutput->Add(fMultA);
   fOutput->Add(fMultPairsA);
   fOutput->Add(fMultTripletsA);
   fOutput->Add(fMultQuadsA);

   fOutput->Add(fMeanQcent0p5);
   fOutput->Add(fTwopartcent0p5);
   fOutput->Add(fMeanQcent);
   fOutput->Add(fTwopartcent);
   fOutput->Add(fMeanQcent5);
   fOutput->Add(fTwopartcent5);
   
   fOutput->Add(fMeanQ);
   fOutput->Add(fTwopart);

   fOutput->Add(fMeanQ10);
   fOutput->Add(fTwopart10);

   fOutput->Add(fMeanQ20);
   fOutput->Add(fTwopart20);
   
   fOutput->Add(fMeanQ50);
   fOutput->Add(fTwopart50);
   
   fOutput->Add(fMeanQ100);
   fOutput->Add(fTwopart100);

   fOutput->Add(fThreepart);
   fOutput->Add(fFourpart);

   fOutput->Add(fMeanQ_ss);
   fOutput->Add(fTwopart_ss);
   fOutput->Add(fThreepart_ss);
   fOutput->Add(fFourpart_ss);

   fOutput->Add(fMeanQcent_ss0p5);
   fOutput->Add(fTwopartcent_ss0p5);
   fOutput->Add(fMeanQcent_ss);
   fOutput->Add(fTwopartcent_ss);
   fOutput->Add(fMeanQcent_ss5);
   fOutput->Add(fTwopartcent_ss5);

   fOutput->Add(fQ1A);
   fOutput->Add(fTwopartA);
   fOutput->Add(fThreepartA);
   fOutput->Add(fFourpartA);

   fOutput->Add(fA);
   fOutput->Add(fPairsA);
   fOutput->Add(fTripletsA);
   fOutput->Add(fQuadsA);
 
   /*   fOutput->Add(h2d_tpc_ITSl1_before);
   fOutput->Add(h2d_tpc_ITSl2_before);
   fOutput->Add(h2d_ITSl1_ITSl2_before);

   fOutput->Add(h2d_tpc_ITSl1_after);
   fOutput->Add(h2d_tpc_ITSl2_after);
   fOutput->Add(h2d_ITSl1_ITSl2_after);*/

   fOutput->Add(fcentEvents);
   fOutput->Add(fcentEvents5);

   fOutput->Add(fEventSee);
   fOutput->Add(_fV0MmultVsSpdTracklet_before);
   fOutput->Add(_fV0MmultVsSpdTracklet_after);
   fOutput->Add(_fhV0MvsTracksTPCout_before);
   fOutput->Add(_fhV0MvsTracksTPCout_after);
   fOutput->Add(_profV0MvsTPCout);
   fOutput-> Add(ftrack);
   //   fOutput->Add(fCorrDet_beforeCut);
   //   fOutput->Add(fCorrDet_afterCut);
   fOutput->Add(fHdcaxy);
   fOutput->Add(fHdcaz);
   fOutput->Add(fTpcNCrossedRows);
   fOutput->Add(fTpcNCluster);

  PostData(1, fOutput); 
}

//Track cuts
//_________________________________________________________________
static Bool_t AcceptTrack(const AliAODTrack *trk){



	if (!trk) return kFALSE;

	//	if(!trk->TestFilterBit(768)) return kFALSE;//Hybrid track cuts: filterbit (768) (ITSRefit required + max shared fraction of TPC clusters < 0.4)
//	if(!trk->TestFilterBit(8)) return kFALSE;	//	if(!trk->TestFilterBit(8)) return kFALSE;
	if (trk->Charge() == 0) return kFALSE;
	if (trk->Eta() < -0.8 || trk->Eta() > 0.8) return kFALSE;
	if (trk->Pt() < 0.15 || trk->Pt() > 2.0)  return kFALSE;
	//if (trk->GetTPCNcls() <=80)return kFALSE;
	//if (!((trk->GetStatus() & AliVTrack::kTPCout) && trk->GetID() > 0 ) )return kFALSE;

	return kTRUE;
	  }
//_________________________________________________________________

void AliAnalysisMeanPtdata::UserExec(Option_t *) 
{


  //  cout<< "I am here 1 "<<endl;
   fEvents->Fill(0);  //before all cuts
   
   fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
   if (!fAOD) {
     printf("ERROR: fAOD not available\n");
     return;
  }
   //    cout<< "I am here 2 "<<endl;
     fEvents->Fill(1);  //after AOD and before pile up cuts

      AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
   AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
   fPIDResponse = inputHandler->GetPIDResponse();


   if (!(((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & AliVEvent::kINT7)) return;    
   
   //   cout<< "I am here 3"<<endl;    

   fEvents->Fill(2);  //after pile up cuts


   // --- Vertex cuts ---
   fPrimaryVtx = fAOD->GetPrimaryVertex();
   if (!fPrimaryVtx){
     printf ("ERROR: no primary vertex\n");
     return;
   }
   fEvents->Fill(3);  //primery vertex                              

   Double_t vtxZGlobal=0., vtxZSPD=0., vtxZTPC=0.;
   Double_t vtxZ=0., vtxNCont=0.;
   Double_t vtxZdiff1=0., vtxZdiff2=0., vtxZdiff3=0.;
   Double_t vtxNContGlobal=0.;
   Double_t vtxZdiff=0.;
   Double_t xv=fPrimaryVtx->GetX();
   Double_t yv=fPrimaryVtx->GetY();
   Double_t zv=fPrimaryVtx->GetZ(); //vtxZGlobal -- in the alice code

   //   const AliESDVertex* vtxESD = fAOD->GetPrimaryVertexTracks();
   //   vtxZGlobal = vtxESD->GetZ();
   //   vtxNContGlobal = vtxESD->GetNContributors();
   // SPD vertex

   const AliAODVertex* vtxAODSPD = fAOD->GetPrimaryVertexSPD();
   vtxZSPD = vtxAODSPD->GetZ();
   //vtxNContSPD = vtxESDSPD->GetNContributors();

   // TPC vertex
   const AliAODVertex* vtxAODTPC = fAOD->GetPrimaryVertexTPC();
   vtxZTPC = vtxAODTPC->GetZ();


   vtxZ = zv; // vertex global
   vtxNContGlobal = fPrimaryVtx->GetNContributors();
   vtxNCont = vtxNContGlobal;


   fVtxZ->Fill(zv);    


   vtxZdiff1 = vtxZTPC - vtxZ;
   vtxZdiff2 = vtxZTPC - vtxZSPD;
   vtxZdiff3 = vtxZ - vtxZSPD;
   fVtxZDiff1->Fill(vtxZdiff1);
   fVtxZDiff2->Fill(vtxZdiff2);
   fVtxZDiff3->Fill(vtxZdiff3);

   /*    if (TMath::Abs(zv) > 10.0)
      {
	return;
	}*/

   //   cout<<" vz b4  "<<zv<<endl;

   if (!(TMath::Abs(zv) < _vZMax))
     {return;}

   //   cout<<" vz after  "<<zv<<endl;


   fEvents->Fill(4);  //primery vertex -10 to 10                              
   fVtxZCut->Fill(vtxZ); // VtxZ after cut on vtxZ

    if(vtxNCont < fNContributors) {
      //     Printf("Vertex has no contributors");
      return;
    }

   fEvents->Fill(5);  //vertex contributor
                            
    fVtxZCont->Fill(vtxZ); // VtxZ after cut on nContributors

    fEvents->Fill(6);  //vertex difference global and TPC
    fVtxZCutDiff->Fill(vtxZ); // VtxZ after cut on vtxZDiff

   
    Int_t event=0;
    event++;

    Double_t fcentrality=0;
    Int_t fCurrentEventCentralityBin ;
    Double_t fITSlayer1=0;
    Double_t fITSlayer2=0;
    Double_t fTPCCluster=0;

 
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");

    if(!MultSelection) {
      //            cout << "AliMultSelection object not found!" << endl;
      return;

    }
	

    else fcentrality = MultSelection->GetMultiplicityPercentile("V0M",kTRUE);
    //    if(fcentrality < 0 || fcentrality >= 100) return;	
    if(fcentrality < 0 || fcentrality >= 90) return;	

    fEvents->Fill(7);  //after pile up cuts

    if(fcentrality > 0.0 && fcentrality < 5.0) fCurrentEventCentralityBin = 0;
    else if(fcentrality >= 5.0 && fcentrality < 10.0 ) fCurrentEventCentralityBin = 1;
    else if(fcentrality >= 10.0 && fcentrality < 15.0 ) fCurrentEventCentralityBin =2;
    else if(fcentrality >= 15.0 && fcentrality < 20.0 ) fCurrentEventCentralityBin =3;
    else if(fcentrality >= 20.0 && fcentrality < 25.0 ) fCurrentEventCentralityBin = 4;
    else if(fcentrality >= 25.0 && fcentrality < 30.0 ) fCurrentEventCentralityBin = 5;
    else if(fcentrality >= 30.0 && fcentrality < 35.0 ) fCurrentEventCentralityBin = 6;
    else if(fcentrality >= 35.0 && fcentrality < 40.0 ) fCurrentEventCentralityBin = 7;
    else if(fcentrality >= 40.0 && fcentrality < 45.0 ) fCurrentEventCentralityBin = 8;
    else if(fcentrality >= 45.0 && fcentrality < 50.0 ) fCurrentEventCentralityBin = 9;
    else if(fcentrality >= 50 && fcentrality <55.0 ) fCurrentEventCentralityBin = 10;
    else if(fcentrality >= 55.0 && fcentrality < 60.0 ) fCurrentEventCentralityBin = 11;
    else if(fcentrality >= 60.0 && fcentrality <65.0 ) fCurrentEventCentralityBin = 12;
    else if(fcentrality >=65.0 && fcentrality < 70.0 ) fCurrentEventCentralityBin =13;
    else if(fcentrality >= 70.0 && fcentrality < 75.0 ) fCurrentEventCentralityBin = 14;
    else if(fcentrality >= 75.0 && fcentrality < 80.0 ) fCurrentEventCentralityBin = 15;
    else if(fcentrality >= 80.0 && fcentrality < 85.0 ) fCurrentEventCentralityBin = 16;
    else if(fcentrality >= 85.0 && fcentrality < 90.0 ) fCurrentEventCentralityBin = 17;
    else if(fcentrality >= 90.0 && fcentrality < 95.0 ) fCurrentEventCentralityBin = 18;
    else if(fcentrality >= 95.0 && fcentrality < 100.0 ) fCurrentEventCentralityBin = 19;
    else  fCurrentEventCentralityBin = -1;
    Int_t centbin =  fCurrentEventCentralityBin;

    //    vsparse_beforeCut[7] = multSelection->GetMultiplicityPercentile("SPDClusters" , kFALSE);
    //  vsparse_beforeCut[8] = multSelection->GetMultiplicityPercentile("SPDTracklets", kFALSE);
    fITSlayer1 = fAOD->GetMultiplicity()->GetNumberOfITSClusters(0);//SPD 1st layer
    fITSlayer2 = fAOD->GetMultiplicity()->GetNumberOfITSClusters(1);//SPD 2nd layer

    fTPCCluster = fAOD->GetNumberOfTPCClusters();//SPD 1st layer
   
    //    h2d_tpc_ITSl1_before->Fill(fITSlayer1,fTPCCluster);
    // h2d_tpc_ITSl2_before->Fill(fITSlayer2,fTPCCluster);
    // h2d_ITSl1_ITSl2_before->Fill(fITSlayer1,fITSlayer2);

    //    cout<<"########fcentrality#######   before "<<fcentrality<<"  "<<fITSlayer1<<"  "<<fITSlayer2<<"   "<<fTPCCluster<<endl;
    fHistCentralityMultSelection->Fill(fcentrality);

    //    cout<< "I am here 5"<<endl;    


    /*        AliVEvent *ev = InputEvent();
     //cout<<"pass cuts: ----------------- "<<fAliEventCuts->AcceptEvent(ev)<<endl;    
     if(!fAliEventCuts->AcceptEvent(ev)) 
       {PostData(1,fOutput);
       return;}
 
       fAliEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);  */
     fEvents->Fill(8);  //after pile up cuts
     //       cout<< "I am here 6"<<endl;    

    //    h2d_tpc_ITSl1_after->Fill(fITSlayer1,fTPCCluster);
    //    h2d_tpc_ITSl2_after->Fill(fITSlayer2,fTPCCluster);
    //    h2d_ITSl1_ITSl2_after->Fill(fITSlayer1,fITSlayer2);

    //     fAliEventCuts->fUseITSTPCCluCorrelationCut = true;
     fEvents->Fill(9);  //after pile up cuts

     //         cout<<"########fcentrality#######   after  "<<fcentrality<<"  "<<fITSlayer1<<"  "<<fITSlayer2<<"   "<<fTPCCluster<<endl;
    Int_t totTrack=0;
    totTrack=fAOD->GetNumberOfTracks();
    //    cout<<"I m here"<<totTrack<<endl;
  
    Int_t nParts = 0 ;
    Double_t Q1=0.0, Q2 = 0.0, Q3 = 0.0, Q4=0.0;
    Double_t Pt=0.;
    Double_t Ptp=0.;
    Double_t Ptn=0.;


    
    Double_t MeanQ1=0., twopart=0., threepart=0., fourpart=0.;
    Double_t  twopart1=0., threepart1=0., fourpart1=0.;
    float spdTracklet = -999.;
    Double_t tracks[100000]= {0.};	// array of track Pts, needed for the two-particle correlator

	    if( centbin == 0){fEventSee->Fill(1);}
	    if( centbin == 1){fEventSee->Fill(2);}
	    if( centbin == 2){fEventSee->Fill(3);}
	    if( centbin == 3){fEventSee->Fill(4);}
	    if( centbin == 4){fEventSee->Fill(5);}
	    if( centbin == 5){fEventSee->Fill(6);}
	    if( centbin == 6){fEventSee->Fill(7);}
	    if( centbin == 7){fEventSee->Fill(8);}
	    if( centbin == 8){fEventSee->Fill(9);}
	    if( centbin == 9){fEventSee->Fill(10);}
	    if( centbin == 10){fEventSee->Fill(11);}
	    if( centbin == 11){fEventSee->Fill(12);}
	    if( centbin == 12){fEventSee->Fill(13);}
	    if( centbin == 13){fEventSee->Fill(14);}
	    if( centbin == 14){fEventSee->Fill(15);}
	    if( centbin == 15){fEventSee->Fill(16);}
	    if( centbin == 16){fEventSee->Fill(17);}
	    if( centbin == 17){fEventSee->Fill(18);}
	    if( centbin == 18){fEventSee->Fill(19);} 
	    if( centbin == 19){fEventSee->Fill(20);} 

    if( fAOD )
      {
	if(_pileUpEvent) // 0=goes out, 1=goes in-pile up cut applied
	  {
	    if (!fAliEventCuts->AcceptEvent(fAOD))//fAliEventCuts //pile up cut applied
	      {PostData(1,fOutput);
		return;}
	    fAliEventCuts->SetRejectTPCPileupWithITSTPCnCluCorr(kTRUE,1);
	    fAliEventCuts->fUseITSTPCCluCorrelationCut = true;
	  }
	    StoreEventMultiplicities( fAOD );

	    spdTracklet = MultSelection->GetMultiplicityPercentile("SPDTracklets");
	    _fhV0MvsTracksTPCout_before->Fill(fV0Multiplicity, fNoOfTPCoutTracks);
	    _fV0MmultVsSpdTracklet_before->Fill(spdTracklet, fV0Multiplicity);

	    if ( _singlesOnly )
	      {
		Double_t vsparse_beforeCut[3];
		vsparse_beforeCut[0] = fNoOfTPCoutTracks;
		vsparse_beforeCut[1] = fV0Multiplicity;
		vsparse_beforeCut[2] = spdTracklet;
		//		fCorrDet_beforeCut->Fill(vsparse_beforeCut);
	      }



    //Tulika -- this is something I will get after extraction of parameters 
	    /*    if(!_usePileupCut_PbPb5TeV){
      fV0MtoTrkTPCout_lower_PbPb5TeV = new TFormula(Form("fV0MtoTrkTPCout_lower_%s", ""), "-42.71378+10.677351*x");//1sigma cut

      fV0MtoTrkTPCout_upper_PbPb5TeV
	= new TFormula(Form("fV0MtoTrkTPCout_upper_%s", ""), "85.77193 +20.053612*x");//2sigma cut
      
      
      if( (fV0Multiplicity < fV0MtoTrkTPCout_lower_PbPb5TeV->Eval(fNoOfTPCoutTracks)) || (fV0Multiplicity > fV0MtoTrkTPCout_upper_PbPb5TeV->Eval(fNoOfTPCoutTracks))) return;
    }


    cout<<fV0Multiplicity<<"   fV0Multiplicity      "<<endl;

       if(_usePileupCut_PbPb5TeV){
     
         if(fNoOfTPCoutTracks >(_hProfPileupCut->GetBinContent(_hProfPileupCut->FindBin(fV0Multiplicity))+.5*_hProfPileupCut->GetBinError(_hProfPileupCut->FindBin(fV0Multiplicity))) || fNoOfTPCoutTracks<(_hProfPileupCut->GetBinContent(_hProfPileupCut->FindBin(fV0Multiplicity))-1.5*_hProfPileupCut->GetBinError(_hProfPileupCut->FindBin(fV0Multiplicity)))) return;
    
	 }*/

       _fhV0MvsTracksTPCout_after->Fill(fV0Multiplicity, fNoOfTPCoutTracks);
       _fV0MmultVsSpdTracklet_after->Fill(spdTracklet, fV0Multiplicity);
       _profV0MvsTPCout->Fill(fV0Multiplicity, fNoOfTPCoutTracks);    //This is the making of the profile

       if ( _singlesOnly )
	 {
	   Double_t vsparse_afterCut[3];
	   vsparse_afterCut[0] = fNoOfTPCoutTracks;
	   vsparse_afterCut[1] = fV0Multiplicity;
	   vsparse_afterCut[2] = spdTracklet;
	   //	   fCorrDet_afterCut->Fill(vsparse_afterCut);
	 }

       /*       AliVEvent *ev = InputEvent();
       if(!fAliEventCuts->AcceptEvent(ev))
	 {PostData(1,fOutput);
	 return;}*/


       Float_t fMultGlobal  = 0;  // global track multiplicity
       Float_t trkDCAz=0.0,trkDCAxy=0.0,dcaxys=0.0,trkDCAxyA=0.0, trkDCAzA=0.0;
       Float_t nCrossedRowsTPC = -1;
       Float_t nTPCNcls = -1;


 for (Int_t iTracks = 0; iTracks <totTrack; iTracks++) {
 

   AliAODTrack* track =(AliAODTrack*)fAOD->GetTrack(iTracks);
  
   if (!track) {
     printf("ERROR: Could not receive track %d\n", iTracks);
     continue;
   }
 


   if (!AcceptTrack(track)) continue;
   if (!(track->TestFilterBit(_fb)))continue;
   fHistEta->Fill(track->Eta());
   Pt =track->Pt();
   fHistPt->Fill(Pt);

   trkDCAxy=  track->DCA();
   trkDCAz=   track->ZAtDCA();

  
     if ((TMath::Abs(trkDCAxy)==999)||(TMath::Abs(trkDCAz)==999))
       {
	 Double_t    bval[2] = {-99., -99.};
	 Double_t    bCov[3] = {-99., -99., -99.};
	 AliAODTrack copy(*track);
	 if(copy.PropagateToDCA(fAOD->GetPrimaryVertex(), fAOD->GetMagneticField(), 100., bval, bCov) && TMath::Abs(bval[0]) < 0.3 &&  TMath::Abs(bval[1]) < 0.3){
	   fMultGlobal++;
	 }

	 trkDCAzA = bval[1];
	 trkDCAxyA = bval[0];
       }
     else
       {
	 trkDCAzA=trkDCAz;
	 trkDCAxyA=trkDCAxy;
       }

     dcaxys=7*(0.0026+(0.005/(TMath::Power(Pt,1.01))));

     //     cout<<"B4 DCA Z    "<< trkDCAzA<<"DCA XY  " <<trkDCAxyA<<endl;
   //   if(!(track->GetTPCchi2perCluster()<_chi2perTPC))continue;

   //  cout<<"Before Dca z  "<<TMath::Abs(trkDCAzA)<<" Dca z max  "<<_dcaZMax<<" dca xy "<<TMath::Abs(trkDCAxyA)<<" dca xy max  "<<_dcaXYMax<<" Tpc crossrows "<<nCrossedRowsTPC<<" Tpc crossrows max  "<<_nTPCCrossRows<<endl;
   

     if (!(TMath::Abs(trkDCAzA)<_dcaZMax && TMath::Abs(trkDCAxyA)<_dcaXYMax))continue; 
     nCrossedRowsTPC = track->GetTPCCrossedRows();
     nTPCNcls=track->GetTPCNcls();

       if(nCrossedRowsTPC < _nTPCCrossRows)continue;
       if(nTPCNcls < _nClusterMin)continue;

     

     //cout<<"After Dca z  "<<TMath::Abs(trkDCAzA)<<" Dca z max  "<<_dcaZMax<<" dca xy "<<TMath::Abs(trkDCAxyA)<<" dca xy max  "<<_dcaXYMax<<" Tpc crossrows "<<nCrossedRowsTPC<<" Tpc crossrows max  "<<_nTPCCrossRows<<endl;

     fHdcaz->Fill(trkDCAzA);
     fHdcaxy->Fill(trkDCAxyA);
     fTpcNCrossedRows->Fill(nCrossedRowsTPC);
     fTpcNCluster->Fill(nTPCNcls);

   
   Int_t nClustersITS = 0;
   nClustersITS =track->GetITSNcls();
   //   if(nClustersITS==0)continue;
   //   Float_t fMinNClustersITS = -1;   
   //if (nClustersITS<fMinNClustersITS) return kFALSE; //cut on minimum number of ITS clusters

   Float_t chi2PerClusterITS = -1;
   chi2PerClusterITS =track->GetITSchi2()/Float_t(nClustersITS);
   //   cout<<" chisq per ITS before  "<<chi2PerClusterITS<<endl;  
   // if(!(chi2PerClusterITS<_chi2perITS))continue;
   //   cout<<" chisq per ITS after  "<<chi2PerClusterITS<<endl;  



   tracks[nParts] = Pt;
   Q1 += Pt;
   Q2 +=  Pt*Pt;
   Q3 += Pt*Pt*Pt;
   Q4 += Pt*Pt*Pt*Pt;
   nParts++;
  
 }
 
 
 ftrack->Fill(nParts);
   Double_t nPairs=0., nTriplets=0., nQuads=0., meanQ1 = 0.0;
   
   
      if(nParts>1){

	meanQ1 = Q1/nParts;
	nPairs = nParts*(nParts-1);//It takes the mean of a particlular
	twopart = ( (Q1*Q1) - Q2);	
	twopart1 = ( (Q1*Q1) - Q2)/(nPairs);
      }

   
       if(nParts > 2) { 
	 nTriplets = nParts*(nParts-1)*(nParts-2);
	 threepart = ((Q1*Q1*Q1) - (3.0*Q2*Q1) + 2*Q3);
	 threepart1 = ((Q1*Q1*Q1) - (3.0*Q2*Q1) + 2*Q3)/(nTriplets);
	 
       }

       
        if(nParts > 3) { 

	  nQuads = nParts*(nParts-1)*(nParts-2)*(nParts - 3);
	  
	  fourpart =  ((Q1*Q1*Q1*Q1) - (6*Q2*Q1*Q1) + (3*Q2*Q2) + (8*Q3*Q1) - 6*Q4) ;
	  
	  fourpart1 =  ((Q1*Q1*Q1*Q1) - (6*Q2*Q1*Q1) + (3*Q2*Q2) + (8*Q3*Q1) - 6*Q4) / (nQuads) ;
	}
	

	eventcount++;

     fCount->Fill(nParts,1);
     fMultMeanQ->Fill(nParts,meanQ1);
     fMultTwopart->Fill(nParts,twopart1);
     fMultThreepart->Fill(nParts,threepart1);
     fMultFourpart->Fill(nParts,fourpart1);

     fMeanQ->Fill(nParts,meanQ1);
     fTwopart->Fill(nParts,twopart1);

     fMeanQ10->Fill(nParts,meanQ1);
     fTwopart10->Fill(nParts,twopart1);

     fMeanQ20->Fill(nParts,meanQ1);
     fTwopart20->Fill(nParts,twopart1);

     fMeanQ50->Fill(nParts,meanQ1);
     fTwopart50->Fill(nParts,twopart1);

     fMeanQ100->Fill(nParts,meanQ1);
     fTwopart100->Fill(nParts,twopart1);

     fThreepart->Fill(nParts,threepart1);
     fFourpart->Fill(nParts,fourpart1);

     fMeanQ_ss->Fill(eventcount%30,nParts,meanQ1);
     fTwopart_ss->Fill(eventcount%30,nParts,twopart1);
     fThreepart_ss->Fill(eventcount%30,nParts,threepart1);
     fFourpart_ss->Fill(eventcount%30,nParts,fourpart1);


     fMeanQcent_ss0p5->Fill(eventcount%30,fcentrality,meanQ1);
     fTwopartcent_ss0p5->Fill(eventcount%30,fcentrality,twopart1);
     fMeanQcent_ss->Fill(eventcount%30,fcentrality,meanQ1);
     fTwopartcent_ss->Fill(eventcount%30,fcentrality,twopart1);
     fMeanQcent_ss5->Fill(eventcount%30,fcentrality,meanQ1);
     fTwopartcent_ss5->Fill(eventcount%30,fcentrality,twopart1);


     fMultQ1->Fill(nParts,Q1);
     fMultTwopartA->Fill(nParts,twopart);
     fMultThreepartA->Fill(nParts,threepart);
     fMultFourpartA->Fill(nParts,fourpart);

      fMultA->Fill(nParts,nParts);
      fMultPairsA->Fill(nParts,nPairs);
      fMultTripletsA->Fill(nParts,nTriplets);
      fMultQuadsA->Fill(nParts,nQuads);

      fQ1A->Fill(nParts,Q1);
      fTwopartA->Fill(nParts,twopart);
      fThreepartA->Fill(nParts,threepart);
      fFourpartA->Fill(nParts,fourpart);

      fA->Fill(nParts,nParts);
      fPairsA->Fill(nParts,nPairs);
      fTripletsA->Fill(nParts,nTriplets);
      fQuadsA->Fill(nParts,nQuads);
     
     fcentnpart->Fill(fcentrality,nParts);
     fcentnpart2d->Fill(fcentrality,nParts);
     fcentnpart2d_1->Fill(fcentrality,nParts);


     fMeanQcent0p5->Fill(fcentrality,meanQ1);
     fTwopartcent0p5->Fill(fcentrality,twopart1);
     fMeanQcent->Fill(fcentrality,meanQ1);
     fTwopartcent->Fill(fcentrality,twopart1);
     fMeanQcent5->Fill(fcentrality,meanQ1);
     fTwopartcent5->Fill(fcentrality,twopart1);
     
     
     fcentEvents->Fill(fcentrality,event);
     fcentEvents5->Fill(fcentrality,event);
      }
    PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisMeanPtdata::Terminate(Option_t *) 
{
	// Draw result to screen, or perform fitting, normalizations
	// Called once at the end of the query
	
	fOutput = dynamic_cast<TList*> (GetOutputData(1));
	if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
	
	fHistPt = dynamic_cast<TH1D*> (fOutput->FindObject("fHistPt"));
	if (!fHistPt) { Printf("ERROR: could not retrieve fHistPt"); return;}
	


}


//________________________________________________________________________
Bool_t AliAnalysisMeanPtdata::StoreEventMultiplicities(AliVEvent *ev)
{  
  AliAODEvent *aodEvent=dynamic_cast<AliAODEvent*>(ev);

  fNoOfTPCoutTracks = 0;
  //fV0Multiplicity_Victor = aodEvent->GetVZEROData()->GetMTotV0A()+event->GetVZEROData()->GetMTotV0C(); //Victor's way

  AliVVZERO *vzero = (AliVVZERO*)ev->GetVZEROData();
  if(vzero)
    {
      fV0Multiplicity = 0;
      fV0AMultiplicity = 0;
      fV0CMultiplicity = 0;

      for(int ich=0; ich < 64; ich++){
        fV0Multiplicity += vzero->GetMultiplicity(ich);
	//        fV0AMultiplicity += vzero->GetMultiplicityV0A(ich);
        //fV0CMultiplicity += vzero->GetMultiplicityV0C(ich);
      } //ich loop   

    } // AliEventCuts
  
  Int_t nTracks = 0;
  nTracks = aodEvent->GetNumberOfTracks();
  
  for (Int_t itrk = 0; itrk < nTracks; itrk++)
    {
      AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(aodEvent->GetTrack(itrk));

      //      cout << "_max_eta_1: "<<_max_eta_1<<"\t" << "_nClusterMin: "<<_nClusterMin<<endl;
      //      if (/*(TMath::Abs(aodt->Eta()) < 0.8) &&*/ (aodt->GetTPCNcls() >= _nClusterMin)/* && (aodt->Pt() >= 0.15) && (aodt->Pt() < 2.)*/)
	      if(AcceptTrack(aodt))
		if (aodt->TestFilterBit(_fb))
		

	{
    	  /*  if ((aodt->GetStatus() & AliVTrack::kTPCout) && aodt->GetID() > 0 )  */ fNoOfTPCoutTracks++;
	}
    }
  return kTRUE;
}
//______________________________________________________________________
//_______________________________________________________________________
