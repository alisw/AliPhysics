/**************************************************************************
 * Authors: Andrey Ivanov, Igor Altsybeev.                                                 *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// This class creates TH2D histogramms for Nch - Nch , Nch - Pt , Pt - Pt 
// dirtributions for given ETA windows and some supplementary data for Long Range Correlation (LRC) analysis .  
// Class is designid to work with AliAnalysisTaskLRC

// Authors : Andrey Ivanov, Igor Altsybeev, St.Peterburg State University
// Email: Igor.Altsybeev@cern.ch

#include "Riostream.h"
#include "AliLRCProcess.h"
//#include "THnSparse.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"
#include "TString.h"
#include "TMath.h"
//#include "TClonesArray.h"
using std::cout;
using std::endl;
using std::cerr;

//#include <AliPID.h> //for particle mass only
ClassImp(AliLRCProcess)

//const bool useSparse = false;//false;
//const bool useAccumulatingHist = true;//false;

AliLRCProcess::AliLRCProcess():fIsEventOpend(kFALSE), fIsOnline(kFALSE), fDisplayInitOnDemandWarning(kTRUE)
  ,fUseSparse(false)
  ,fUseAccumulatingHist(true)
  ,fEventCount(0),fStartForwardETA(0), fEndForwardETA(0)
  ,fStartForwardPhi(0)
  ,fEndForwardPhi(0)
  ,fStartBackwardETA(0)
  ,fEndBackwardETA(0)
  ,fStartBackwardPhi(0)
  ,fEndBackwardPhi(0)
  ,fDoubleSidedBackwardPhiWindow(kFALSE)
  ,fHiPt(0)
  ,fLoPt(0)
  ,fHiMultHor(0)
  ,fLowMultHor(0)
  ,fHiMultVert(0)
  ,fLowMultVert(0)
  ,fMultBinsHor(0)
  ,fMultBinsVert(0)
  ,fPtBins(0)
  ,fPtHistXaxisRebinFactor(1)
  ,fSumPtFw(0)
  ,fSumPtBw(0)
  ,fSumPtBw2(0)
  ,fNchFw(0)
  ,fNchBw(0)
  ,fNchFwPlus(0)
  ,fNchBwPlus(0)
  ,fNchFwMinus(0)
  ,fNchBwMinus(0)
  ,fOutList(0), fShortDef(0),fOutputSlot(0), fHistPt(0),fHistEta(0),fHistClouds(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fProfNberr(0),fProfNberrPtPt(0),fProfdPtB(0),fProfTestLRC(0),fHistSparseDimensionLabeling(0),fHistSparsePIDblocksLabeling(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistNchForwardPtPt(0)
  ,fHistPhiForward(0)
  ,fHistTracksChargeForward(0)
  ,fHistPtBackward(0)
  ,fHistEtaBackward(0)
  ,fHistNchBackward(0)
  ,fHistPhiBackward(0)
  ,fHistTracksChargeBackward(0)
  ,fArrAccumulatedValues(0)
  //  ,fHistArrayLabeling(0)
  ,fEventCentrality(0)
  ,fHistNfCentrality(0)
  ,fHistTestPIDForward(0)
  ,fHistTestPIDBackward(0)
  ,fHistDifferenceNf(0)
  ,fHistDifferenceNb(0)

  ,fPidForward(-1)//kLRCany)
  ,fPidBackward(-1)//kLRCany)
  //,fWhichParticleToProcess(kLRCany)//kLRCany)
  //,fPidFillCondition(kLRCpidIgnored)
  //,fNumberOfSectors(1)
  //,fNeedToRotateSector(kFALSE)
{
    //fWhichParticleToProcess = kLRCany;  //default - all particle types
    //fPidFillCondition = kLRCpidIgnored;
}

AliLRCProcess::AliLRCProcess(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA ): fIsEventOpend(kFALSE), fIsOnline(kFALSE), fDisplayInitOnDemandWarning(kTRUE)
  ,fUseSparse(false)
  ,fUseAccumulatingHist(true)
  ,fEventCount(0),fStartForwardETA(0), fEndForwardETA(0), fStartForwardPhi(0),fEndForwardPhi(0),fStartBackwardETA(0), fEndBackwardETA(0),fStartBackwardPhi(0)
  ,fEndBackwardPhi(0)
  ,fDoubleSidedBackwardPhiWindow(kFALSE)
  ,fHiPt(0)
  ,fLoPt(0)
  ,fHiMultHor(0)
  ,fLowMultHor(0)
  ,fHiMultVert(0)
  ,fLowMultVert(0)
  ,fMultBinsHor(0)
  ,fMultBinsVert(0)
  ,fPtBins(0)
  ,fPtHistXaxisRebinFactor(1)
  ,fSumPtFw(0),  fSumPtBw(0), fSumPtBw2(0),fNchFw(0)
  ,/*fNchFwPtPt(0),*/ fNchBw(0)
  ,fNchFwPlus(0)
  ,fNchBwPlus(0)
  ,fNchFwMinus(0)
  ,fNchBwMinus(0)
  ,fOutList(0), fShortDef(0),fOutputSlot(0), fHistPt(0),fHistEta(0),fHistClouds(0),fHistNN(0),fHistPtN(0),fHistPtPt(0),fProfNberr(0),fProfNberrPtPt(0),fProfdPtB(0),fProfTestLRC(0),fHistSparseDimensionLabeling(0),fHistSparsePIDblocksLabeling(0),fHistPtForward(0),fHistEtaForward(0),fHistNchForward(0),fHistNchForwardPtPt(0)
  ,fHistPhiForward(0)
  ,fHistTracksChargeForward(0)
  ,fHistPtBackward(0)
  ,fHistEtaBackward(0)
  ,fHistNchBackward(0)
  ,fHistPhiBackward(0)
  ,fHistTracksChargeBackward(0)
  ,fArrAccumulatedValues(0)
  //  ,fHistArrayLabeling(0)
  ,fEventCentrality(0)
  ,fHistNfCentrality(0)
  ,fHistTestPIDForward(0)
  ,fHistTestPIDBackward(0)
  ,fHistDifferenceNf(0)
  ,fHistDifferenceNb(0)
  ,fPidForward(-1)//kLRCany)
  ,fPidBackward(-1)//kLRCany)
  //,fWhichParticleToProcess(kLRCany)//kLRCany)
  //,fPidFillCondition(kLRCpidIgnored)
  //,fNumberOfSectors(1)
  //,fNeedToRotateSector(kFALSE)
{
    // Constructor with window setup makes ready-to-run processor
    fEventCount=0;

    //fWhichParticleToProcess = kLRCany;  //default - all particle types
    //fPidFillCondition = kLRCpidIgnored;
    
    //cout << "TEST" << endl;
    SetETAWindows( _StartForwardETA, _EndForwardETA,_StartBakwardETA,_EndBakwardETA);
    SetHistPtRange( 0.15, 1.5, 270 );
    SetHistMultRange( 0, 0, 100 );
    SetForwardWindowPhi( 0, 2*TMath::Pi() );
    SetBackwardWindowPhi( 0, 2*TMath::Pi() );
    

}

Bool_t AliLRCProcess::InitDataMembers()
{
    //Printf("INITDATAMEMBERS");
    // This method is actualy creating output histogramms
    // Thist method  is to be called in CreateOutputObjects method of AliAnalysisTask
    //cout<<" # Init for "<<fShortDef<<" this="<<this<<"\n";
    if( fIsOnline )
    {
        Printf("Can't init data members more then one time! \n");
        return kFALSE;
    }
    fEventCount=0;
    fOutList = new TList();
    fOutList->SetOwner();  // IMPORTANT!

    fOutList->SetName(fShortDef);

    Double_t lowMultHor, hiMultHor;
    Double_t lowMultVert, hiMultVert;

    lowMultHor = fLowMultHor - 0.5;
    hiMultHor = fHiMultHor + 0.5;

    lowMultVert = fLowMultVert - 0.5;
    hiMultVert = fHiMultVert + 0.5;


    //TArray to accumulate data, with names hist
    //26.01.2013: array with accumulated values
    //fArrAccumulatedValues = new TClonesArray("Float_t", en_arr_labels_total );//TArrayF(en_arr_labels_total);
    if ( fUseAccumulatingHist )
    {
        fArrAccumulatedValues = new TH1D( "fArrAccumulatedValues", "Accumulating hist with labeling", en_arr_labels_total,-0.5,en_arr_labels_total-0.5);
        TString gArrayMemberNames[en_arr_labels_total];
        gArrayMemberNames[ en_arr_labels_NN_Nevents  ] = "NN_Nevents"       ;
        gArrayMemberNames[ en_arr_labels_NN_Nf       ] = "NN_Nf"            ;
        gArrayMemberNames[ en_arr_labels_NN_Nb       ] = "NN_Nb"            ;
        gArrayMemberNames[ en_arr_labels_NN_N2_f     ] = "NN_N2_f"          ;
        gArrayMemberNames[ en_arr_labels_NN_Nf_Nb    ] = "NN_Nf_Nb"         ;

        gArrayMemberNames[ en_arr_labels_PtN_Nevents ] = "PtN_Nevents"      ;
        gArrayMemberNames[ en_arr_labels_PtN_Nf      ] = "PtN_Nf"           ;
        gArrayMemberNames[ en_arr_labels_PtN_PtB     ] = "PtN_PtB"          ;
        gArrayMemberNames[ en_arr_labels_PtN_N2_f    ] = "PtN_N2_f"         ;
        gArrayMemberNames[ en_arr_labels_PtN_Ptb_Nf  ] = "PtN_Ptb_Nf"       ;

        gArrayMemberNames[ en_arr_labels_PtPt_Nevents] = "PtPt_Nevents"     ;
        gArrayMemberNames[ en_arr_labels_PtPt_PtF    ] = "PtPt_PtF"         ;
        gArrayMemberNames[ en_arr_labels_PtPt_PtB    ] = "PtPt_PtB"         ;
        gArrayMemberNames[ en_arr_labels_PtPt_Pt2_f  ] = "PtPt_Pt2_f"       ;
        gArrayMemberNames[ en_arr_labels_PtPt_Ptf_Ptb] = "PtPt_Ptf_Ptb"     ;

        for( Int_t i = 1; i <= en_arr_labels_total; i++ )
            fArrAccumulatedValues->GetXaxis()->SetBinLabel(i,gArrayMemberNames[i-1].Data());
        //fOutList->Add(fArrAccumulatedValues);
        fOutList->Add(fArrAccumulatedValues);
    }



    
    // Window statistics histograms

    // Forward

    fHistPtForward = new TH1D("fHistPtForward", "P_{T} distribution in Forward window", 100, 0.0, 5);
    fHistPtForward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtForward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtForward->SetMarkerStyle(kFullCircle);


    fHistEtaForward = new TH1D("fEtaForward", "#eta distribution in Forward window", 200, -2, 2);
    fHistEtaForward->GetXaxis()->SetTitle("ETA");
    fHistEtaForward->GetYaxis()->SetTitle("dN/ETA");
    fHistEtaForward->SetMarkerStyle(kFullCircle);

    
    fHistNchForward = new TH1D("fHistNchForward", "N_{ch} distribution in Forward window", fMultBinsHor, lowMultHor, hiMultHor);
    fHistNchForward->GetXaxis()->SetTitle("N_{ch}");
    fHistNchForward->GetYaxis()->SetTitle("dN/dN_{ch}");
    fHistNchForward->SetMarkerStyle(kFullCircle);

    fHistNchForwardPtPt = new TH1D("fHistNchForwardPtPt", "N_{ch} distribution in Forward window for PtPt accept conditions", fMultBinsHor, lowMultHor, hiMultHor);
    fHistNchForwardPtPt->GetXaxis()->SetTitle("N_{ch}");
    fHistNchForwardPtPt->GetYaxis()->SetTitle("dN/dN_{ch}");
    fHistNchForwardPtPt->SetMarkerStyle(kFullCircle);

    fHistPhiForward = new TH1D("fPhiForward", "#phi distribution in Forward window", 144, 0, 2*TMath::Pi());
    fHistPhiForward->GetXaxis()->SetTitle("Phi");
    fHistPhiForward->GetYaxis()->SetTitle("dN/Phi");
    fHistPhiForward->SetMarkerStyle(kFullCircle);

    fHistTestPIDForward = new TH1D("fHistTestPIDForward","PID distribution in Forward window;PID;N",5,-0.5,4.5);
    TString gBinParticleNames[5] = {/*"Other",*/"Electron","Muon","Pion","Kaon", "Proton"};
    for(Int_t i = 1; i <= 5; i++)
        fHistTestPIDForward->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());

    //15/12/2012: charge hist
    fHistTracksChargeForward = new TH1D("fHistTracksChargeForward","Accepted tracks charge;charge;Entries",3,-1.5,1.5);

    // Backward

    fHistPtBackward = new TH1D("fHistPtBakward", "P_{T} distribution in Backward window", 100, 0.0, 5);
    fHistPtBackward->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPtBackward->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPtBackward->SetMarkerStyle(kFullCircle);


    fHistEtaBackward = new TH1D("fEtaBakward", "#eta distribution in Backward window", 200, -2, 2);
    fHistEtaBackward->GetXaxis()->SetTitle("ETA");
    fHistEtaBackward->GetYaxis()->SetTitle("dN/ETA");
    fHistEtaBackward->SetMarkerStyle(kFullCircle);


    fHistNchBackward = new TH1D("fHistNchBakward", "N_{ch} distribution in Backward window", fMultBinsVert, lowMultVert, hiMultVert);
    fHistNchBackward->GetXaxis()->SetTitle("N_{ch}");
    fHistNchBackward->GetYaxis()->SetTitle("dN/dN_{ch}");
    fHistNchBackward->SetMarkerStyle(kFullCircle);

    fHistPhiBackward = new TH1D("fPhiBakward", "#phi distribution in Backward window", 144, 0, 2*TMath::Pi());
    fHistPhiBackward->GetXaxis()->SetTitle("Phi");
    fHistPhiBackward->GetYaxis()->SetTitle("dN/Phi");
    fHistPhiBackward->SetMarkerStyle(kFullCircle);

    fHistTestPIDBackward = new TH1D("fHistTestPIDBackward","PID distribution in Backward window;PID;N",5,-0.5,4.5);
    for(Int_t i = 1; i <= 5; i++)
        fHistTestPIDBackward->GetXaxis()->SetBinLabel(i,gBinParticleNames[i-1].Data());


    //15/12/2012: charge hist
    fHistTracksChargeBackward = new TH1D("fHistTracksChargeBackward","Accepted tracks charge;charge;Entries",3,-1.5,1.5);






    //Overal statistics histograms

    fHistPt = new TH1F("fHistPt", "P_{T} distribution", 100, 0.0, 5.0);
    fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
    fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
    fHistPt->SetMarkerStyle(kFullCircle);


    fHistEta = new TH1F("fHistEta", "#eta distribution", 200, -2, 2);
    fHistEta->GetXaxis()->SetTitle("ETA");
    fHistEta->GetYaxis()->SetTitle("dN/ETA");
    fHistEta->SetMarkerStyle(kFullCircle);
    
    
    
    // -------- LRC histograms
    //new cloud implementation
    //const int lSparseDim = en_sparse_total;
    //const int nSparseBins = 1000;
    fCorrespondanceWithAliROOTpid[kSparsePIDany] = kSparsePIDany;
    fCorrespondanceWithAliROOTpid[kSparsePIDdefined] = -1000;
    fCorrespondanceWithAliROOTpid[kSparsePIDpion] = 2;
    fCorrespondanceWithAliROOTpid[kSparsePIDkaon] = 3;
    fCorrespondanceWithAliROOTpid[kSparsePIDproton] = 4;

    /* from AliROOT //deprecated!
        enum EParticleType {
            kElectron = 0,
            kMuon = 1,
            kPion = 2,
            kKaon = 3,
            kProton = 4,
            kPhoton = 5,
            kPi0 = 6,
            kNeutron = 7,
            kKaon0 = 8,
            kEleCon = 9,
            kDeuteron = 10,
            kTriton = 11,
            kHe3 = 12,
            kAlpha = 13,
            kUnknown = 14
        };*/


    if ( fUseSparse )
    {

        //correspondance of PID blocks
        //it's a way to unlink THnSparse data dimenstion from enum
        fHistSparsePIDblocksLabeling = new TH1D("fHistSparsePIDblocksLabeling","THnSparse PID blocks labeling", kSparsePIDtotal,-0.5,kSparsePIDtotal-0.5);
        TString gEventCutBinPIDblocksNames[kSparsePIDtotal];    // = {"Total","No trigger","Centrality","No vertex","Bad vertex position","HighMult cut","LowMult cut","Analyzed"};
        gEventCutBinPIDblocksNames[kSparsePIDany]               = "any";
        gEventCutBinPIDblocksNames[kSparsePIDdefined]       = "defined";
        gEventCutBinPIDblocksNames[kSparsePIDpion]      = "pion";
        gEventCutBinPIDblocksNames[kSparsePIDkaon]      = "kaon";
        gEventCutBinPIDblocksNames[kSparsePIDproton]    = "proton";


        for(Int_t i = 1; i <= kSparsePIDtotal; i++)fHistSparsePIDblocksLabeling->GetXaxis()->SetBinLabel(i,gEventCutBinPIDblocksNames[i-1].Data());
        //for(Int_t i = 0; i < nEnumBins; i++)fHistSparseDimensionLabeling->Fill( i );

        //dimensions labelling

        fHistSparseDimensionLabeling = new TH1D("fHistSparseDimensionLabeling","THnSparse labeling", kSparseTotal,-0.5,kSparseTotal-0.5);
        TString gSparseDimensionsNames[kSparseTotal];    // = {"Total","No trigger","Centrality","No vertex","Bad vertex position","HighMult cut","LowMult cut","Analyzed"};
        gSparseDimensionsNames[kSparseNf] = "N_f";
        gSparseDimensionsNames[kSparseNb] = "N_b";
        gSparseDimensionsNames[kSparsePtF] = "Pt_f";
        gSparseDimensionsNames[kSparsePtB] = "Pt_b";
        gSparseDimensionsNames[en_sparse_N2_f] = "N2_f";
        gSparseDimensionsNames[en_sparse_Nf_Nb] = "Nf_Nb";
        gSparseDimensionsNames[en_sparse_Ptb_Nf] = "Ptb_Nf";
        gSparseDimensionsNames[en_sparse_Pt2_f] = "Pt2_f";
        gSparseDimensionsNames[en_sparse_Ptf_Ptb] = "Ptf_Ptb";

        for( Int_t i = 1; i <= kSparseTotal; i++ )
            fHistSparseDimensionLabeling->GetXaxis()->SetBinLabel(i,gSparseDimensionsNames[i-1].Data());



        Int_t* lSparseBins 	= new Int_t[kSparseTotal*kSparsePIDtotal];
        Double_t *lSparseXmin 	= new Double_t[kSparseTotal*kSparsePIDtotal];
        Double_t *lSparseXmax 	= new Double_t[kSparseTotal*kSparsePIDtotal];
        TString *lSparseAxisNames 	= new TString[kSparseTotal*kSparsePIDtotal];

        TString *lPIDNames 	= new TString[kSparsePIDtotal];
        lPIDNames[ kSparsePIDany      ] = Form( "any" );
        lPIDNames[ kSparsePIDdefined  ] = Form( "defined" );
        lPIDNames[ kSparsePIDpion     ] = Form( "pion" );
        lPIDNames[ kSparsePIDkaon     ] = Form( "kaon" );
        lPIDNames[ kSparsePIDproton   ] = Form( "proton" );



        for ( Int_t d = 0; d < kSparsePIDtotal; ++d )
        {
            Int_t binShift = kSparseTotal*d;

            lSparseAxisNames[kSparseNf + binShift] = Form( "axisNf_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[kSparseNb + binShift] = Form( "axisNb_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[kSparsePtF + binShift] = Form( "axisPtf_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[kSparsePtB + binShift] = Form( "axisPtb_%s", lPIDNames[ d ].Data() );

            lSparseAxisNames[en_sparse_N2_f + binShift]     = Form( "axisN2_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[en_sparse_Nf_Nb + binShift]    = Form( "axisNf_Nb_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[en_sparse_Ptb_Nf + binShift]   = Form( "axisPtb_Nf_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[en_sparse_Pt2_f + binShift]    = Form( "axisPt2_f_%s", lPIDNames[ d ].Data() );
            lSparseAxisNames[en_sparse_Ptf_Ptb + binShift]  = Form( "axisPtf_Ptb_%s", lPIDNames[ d ].Data() );


            lSparseBins[kSparseNf + binShift ] = fMultBinsHor;
            lSparseXmin[kSparseNf + binShift ] = lowMultHor;
            lSparseXmax[kSparseNf + binShift ] = hiMultHor;
            lSparseBins[kSparseNb + binShift ] = fMultBinsVert;
            lSparseXmin[kSparseNb + binShift ] = lowMultVert;
            lSparseXmax[kSparseNb + binShift ] = hiMultVert;
            //}
            //for (Int_t d = 2; d < lSparseDim; ++d) {
            lSparseBins[kSparsePtF + binShift ] = fPtBins;
            lSparseXmin[kSparsePtF + binShift ] = fLoPt;
            lSparseXmax[kSparsePtF + binShift ] = fHiPt;
            lSparseBins[kSparsePtB + binShift ] = fPtBins;
            lSparseXmin[kSparsePtB + binShift ] = fLoPt;
            lSparseXmax[kSparsePtB + binShift ] = fHiPt;





            //}
            /*
                lSparseBins[en_sparse_Et_f + binShift ] = 500;
            lSparseXmin[en_sparse_Et_f + binShift ] = 0.2;
            lSparseXmax[en_sparse_Et_f + binShift ] = 2.7;
            lSparseBins[en_sparse_Et_b + binShift ] = 500;
            lSparseXmin[en_sparse_Et_b + binShift ] = 0.2;
            lSparseXmax[en_sparse_Et_b + binShift ] = 2.7;
*/

            lSparseBins[en_sparse_N2_f + binShift ] = (fMultBinsHor-1)*(fMultBinsHor-1)+1;
            lSparseXmin[en_sparse_N2_f + binShift ] = fLowMultHor*fLowMultHor-0.5; // ! use global mult without shift
            lSparseXmax[en_sparse_N2_f + binShift ] = fHiMultHor*fHiMultHor+0.5;
            //lSparseBins[en_sparse_N_b] = fMultBinsVert;
            //lSparseXmin[en_sparse_N_b] = lowMultVert;
            //lSparseXmax[en_sparse_N_b] = hiMultVert;
            lSparseBins[en_sparse_Nf_Nb + binShift ] = (fMultBinsHor-1)*(fMultBinsVert-1)+1;
            lSparseXmin[en_sparse_Nf_Nb + binShift ] = fLowMultHor*fLowMultVert-0.5;
            lSparseXmax[en_sparse_Nf_Nb + binShift ] = fHiMultHor*fHiMultVert+0.5;



            lSparseBins[en_sparse_Ptb_Nf + binShift ] = /*(fMultBinsHor-1)**/(10*fPtBins);//+1;
            lSparseXmin[en_sparse_Ptb_Nf + binShift ] = fLowMultHor*fLoPt-0.5;
            lSparseXmax[en_sparse_Ptb_Nf + binShift ] = fHiMultHor*fHiPt+0.5;

            lSparseBins[en_sparse_Pt2_f + binShift ] = (10*fPtBins);//*(fPtBins);//+1;
            lSparseXmin[en_sparse_Pt2_f + binShift ] = fLoPt*fLoPt;
            lSparseXmax[en_sparse_Pt2_f + binShift ] = fHiPt*fHiPt;

            lSparseBins[en_sparse_Ptf_Ptb + binShift ] = (10*fPtBins);//*(fPtBins)+1;
            lSparseXmin[en_sparse_Ptf_Ptb + binShift ] = fLoPt*fLoPt;
            lSparseXmax[en_sparse_Ptf_Ptb + binShift ] = fHiPt*fHiPt;

            /*
                lSparseBins[en_sparse_Nf_plus + binShift ] = fMultBinsHor;
            lSparseXmin[en_sparse_Nf_plus + binShift ] = lowMultHor;
            lSparseXmax[en_sparse_Nf_plus + binShift ] = hiMultHor;
                lSparseBins[en_sparse_Nf_minus + binShift ] = fMultBinsHor;
            lSparseXmin[en_sparse_Nf_minus + binShift ] = lowMultHor;
            lSparseXmax[en_sparse_Nf_minus + binShift ] = hiMultHor;

                lSparseBins[en_sparse_Nb_plus + binShift ] = fMultBinsVert;
            lSparseXmin[en_sparse_Nb_plus + binShift ] = lowMultVert;
            lSparseXmax[en_sparse_Nb_plus + binShift ] = hiMultVert;
                lSparseBins[en_sparse_Nb_minus + binShift ] = fMultBinsVert;
            lSparseXmin[en_sparse_Nb_minus + binShift ] = lowMultVert;
            lSparseXmax[en_sparse_Nb_minus + binShift ] = hiMultVert;

                lSparseBins[en_sparse_Nf_plus_Nb_minus + binShift ] = (fMultBinsHor-1)*(fMultBinsVert-1)+1;
            lSparseXmin[en_sparse_Nf_plus_Nb_minus + binShift ] = fLowMultHor*fLowMultVert-0.5;
            lSparseXmax[en_sparse_Nf_plus_Nb_minus + binShift ] = fHiMultHor*fHiMultVert+0.5;

                lSparseBins[en_sparse_Nb_plus_Nf_minus + binShift ] = (fMultBinsHor-1)*(fMultBinsVert-1);
            lSparseXmin[en_sparse_Nb_plus_Nf_minus + binShift ] = fLowMultHor*fLowMultVert-0.5;
            lSparseXmax[en_sparse_Nb_plus_Nf_minus + binShift ] = fHiMultHor*fHiMultVert+0.5;
        */
        }

        fHistClouds = new THnSparseD("cloudLRC", "cloudLRC", kSparseTotal*kSparsePIDtotal, lSparseBins, lSparseXmin, lSparseXmax);
        //end of cloud implementation

        //set axis names
        TAxis *lSparseAxis = 0x0;
        for ( Int_t d = 0; d < kSparseTotal*kSparsePIDtotal; ++d )
        {
            lSparseAxis = fHistClouds->GetAxis( d );
            lSparseAxis->SetNameTitle( lSparseAxisNames[d], lSparseAxisNames[d] );
        }


        delete [] lSparseBins;
        delete [] lSparseXmin;
        delete [] lSparseXmax;
        delete [] lSparseAxisNames;


        fOutList->Add(fHistSparseDimensionLabeling);
        fOutList->Add(fHistSparsePIDblocksLabeling);

        //    bool useSparse = false;

        // !!!!!! temp comment!
        fOutList->Add(fHistClouds);
    }


    fHistNN = new TH2D("NN","NN",fMultBinsHor, lowMultHor, hiMultHor,fMultBinsVert, lowMultVert, hiMultVert );
    fHistPtN = new TH2D("PtN","PtN",fMultBinsHor, lowMultHor, hiMultHor,fPtBins/* /fPtHistXaxisRebinFactor*/,fLoPt,fHiPt);
    fHistPtPt = new TH2D("PtPt","PtPt",fPtBins/fPtHistXaxisRebinFactor,fLoPt,fHiPt,fPtBins,fLoPt,fHiPt);
    fProfNberr = new TProfile("nber","nber",fMultBinsHor, lowMultHor, hiMultHor);
    fProfNberrPtPt = new TProfile("nberPtPt","nberPtPt",fPtBins/fPtHistXaxisRebinFactor,fLoPt,fHiPt);
    fProfdPtB = new TProfile("dPtB","Overal multievent Pt_Backward (first bin) Pt_Backward^2 (sec. bin) ",16,0.5,16.5);
    fProfTestLRC = new TProfile("TestLRC","Test LRC calculaion via TProfile",fMultBinsHor, lowMultHor, hiMultHor);

    fHistNfCentrality = new TH2D("NfCentrality","NfCentrality",fMultBinsHor, lowMultHor, hiMultHor,101,-1.01,100.01);
    fHistDifferenceNf = new TH2D("fHistDifferenceNf","Hist nF-nB;nF;nF-nB",fMultBinsHor, lowMultHor, hiMultHor,fMultBinsHor,-hiMultHor,hiMultHor);
    fHistDifferenceNb = new TH2D("fHistDifferenceNb","Hist nB-nF;nB;nB-nF",fMultBinsVert, lowMultVert, hiMultVert,fMultBinsVert,-hiMultVert,hiMultVert);

    // ---------- Adding data members to output list

    // Adding overal statistics

    //commented: to save memory
    //fOutList->Add(fHistPt);
    //fOutList->Add(fHistEta);

    //Adding LRC hists

    
    if (1)
    {
        fOutList->Add(fHistNN);
        fOutList->Add(fHistPtN);
        fOutList->Add(fHistPtPt);
    }
    fOutList->Add(fProfNberr);
    fOutList->Add(fProfNberrPtPt);
    fOutList->Add(fProfdPtB);
    fOutList->Add(fProfTestLRC);


    //Adding window statistics
    


    fOutList->Add(fHistNchForward);
    fOutList->Add(fHistNchBackward);
    fOutList->Add(fHistNchForwardPtPt);

    fOutList->Add(fHistPtForward);
    fOutList->Add(fHistPtBackward);

    fOutList->Add(fHistEtaForward);
    fOutList->Add(fHistEtaBackward);

    fOutList->Add(fHistPhiForward);
    fOutList->Add(fHistPhiBackward);

    fOutList->Add(fHistTracksChargeForward);
    fOutList->Add(fHistTracksChargeBackward);

    fOutList->Add(fHistTestPIDForward);
    fOutList->Add(fHistTestPIDBackward);

    //    fOutList->Add(fHistNfCentrality);
    
    
    //fOutList->Add(fHistDifferenceNf);
    //fOutList->Add(fHistDifferenceNb);

    // Adding status to dPtB

    fProfdPtB->Fill(3 , fStartForwardETA);
    fProfdPtB->Fill(4 , fEndForwardETA);
    fProfdPtB->Fill(5 , fStartBackwardETA);
    fProfdPtB->Fill(6 , fEndBackwardETA);
    fProfdPtB->Fill(7 , fStartForwardPhi);
    fProfdPtB->Fill(8 , fEndForwardPhi);
    fProfdPtB->Fill(9 , fStartBackwardPhi);
    fProfdPtB->Fill(10 , fEndBackwardPhi);




    fIsOnline = kTRUE;
    return kTRUE;
}
AliLRCProcess::~AliLRCProcess()
{
    //Destructor

}

// ---------------------------------------  Setters ------------------
void AliLRCProcess::SetShortDef()
{
    // Creating task and output container name
    char str[200];
    snprintf(str,200, "TaskLRCw%3.1fto%3.1fvs%3.1fto%3.1f",fStartForwardETA,fEndForwardETA,fStartBackwardETA,fEndBackwardETA);
    /*if ( fWhichParticleToProcess != kLRCany
        && (int)fWhichParticleToProcess > 0 && (int)fWhichParticleToProcess <= 6 ) //to avoid program falling
    {
        char str2[80];
        TString gBinParticleNames[6] = {"Other","Electron","Muon","Pion","Kaon", "Proton"};
        snprintf(str2,80, "%s_%s",str,gBinParticleNames[(int)fWhichParticleToProcess].Data());
        fShortDef= str2;
    }
    else*/
    fShortDef= str;

}

void AliLRCProcess::SetForwardWindow(Double_t StartETA,Double_t EndETA)
{
    //setter for the forward eta window
    fStartForwardETA=StartETA;
    fEndForwardETA=EndETA;
    SetShortDef();
}
void AliLRCProcess::SetBackwardWindow(Double_t StartETA,Double_t EndETA)
{
    //setter for the backward eta window
    fStartBackwardETA=StartETA;
    fEndBackwardETA=EndETA;
    SetShortDef();
}
void AliLRCProcess::SetETAWindows(Double_t _StartForwardETA,Double_t _EndForwardETA,Double_t _StartBakwardETA,Double_t _EndBakwardETA)
{
    //setter for the eta windows
    fStartForwardETA=_StartForwardETA;
    fEndForwardETA=_EndForwardETA;
    fStartBackwardETA=_StartBakwardETA;
    fEndBackwardETA=_EndBakwardETA;
    SetShortDef();
}
void AliLRCProcess::GetETAWindows(Double_t &_StartForwardETA,Double_t &_EndForwardETA,Double_t &_StartBakwardETA,Double_t &_EndBakwardETA)
{
    //getter for the eta windows
    _StartForwardETA    = fStartForwardETA;
    _EndForwardETA      = fEndForwardETA;
    _StartBakwardETA    = fStartBackwardETA;
    _EndBakwardETA      = fEndBackwardETA;
}

void AliLRCProcess::GetPhiWindows(Double_t &_StartForwardPhi,Double_t &_EndForwardPhi,Double_t &_StartBakwardPhi,Double_t &_EndBakwardPhi)
{
    //getter for the eta windows
    _StartForwardPhi    = fStartForwardPhi;
    _EndForwardPhi      = fEndForwardPhi;
    _StartBakwardPhi    = fStartBackwardPhi;
    _EndBakwardPhi      = fEndBackwardPhi;
}

void AliLRCProcess::SetParticleType( char* strForwardOrBackward, char* strPid )
{
    //cout << "hm! strForwardOrBackward = " << strForwardOrBackward
    //	<< ", strPid = " << strPid << endl;
    //cout << "before ae! fPidForward = " << fPidForward << ", fPidBackward = " << fPidBackward << endl;

    int lPid = -1;//kLRCany;
    if ( !strcmp( strPid, "pion") )
        lPid = 2;//kLRCpion;
    else if ( !strcmp( strPid, "kaon") )
        lPid = 3;//kLRCkaon;
    else if ( !strcmp( strPid, "proton") )
        lPid = 4;//kLRCproton;
    else if ( !strcmp( strPid, "knownpid") )
        lPid = 100;//will will histos if we KNOW PID! (not important which)

    //set pid for window
    if ( !strcmp( strForwardOrBackward, "fwd") )
        fPidForward = lPid;
    else if ( !strcmp( strForwardOrBackward, "bkwd") )
        fPidBackward = lPid;
    //cout << "ae! lPid = " << lPid << ", fPidForward = " << fPidForward << ", fPidBackward = " << fPidBackward << endl;
    //int aaaa;
    //cin>> aaaa;
}


void AliLRCProcess::SetHistPtRangeForwardWindowRebinFactor( Int_t ptHistXaxisRebinFactor )
{
    // Rebining for Pt histograms X-axis
    if(fIsOnline)
    {
        Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
        return ;
    }
    fPtHistXaxisRebinFactor = ptHistXaxisRebinFactor;
}

void AliLRCProcess::SetHistPtRange(Double_t LoPt,Double_t HiPt,Int_t PtBins)
{
    // Sets Pt range and number of bins for Pt axis of histos
    if(fIsOnline)
    {
        Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
        return ;
    }
    // Setter for Pt range and N bins in histos
    fLoPt=LoPt;
    fHiPt=HiPt;
    fPtBins=PtBins;
}
void AliLRCProcess::SetHistMultRange( Int_t whichWindow, Int_t LoMult,Int_t HiMult,Int_t MultBins)
{
    // Setter for multiplicity range and N bins in histos
    if ( whichWindow == 0 ) //set range for both windows
    {
        SetHistMultRangeHor( LoMult, HiMult, MultBins) ;
        SetHistMultRangeVert( LoMult, HiMult, MultBins) ;
    }
    else if ( whichWindow == 1 ) //for fwd
        SetHistMultRangeHor( LoMult, HiMult, MultBins) ;
    else if ( whichWindow == 2 ) //for bwd
        SetHistMultRangeVert( LoMult, HiMult, MultBins) ;
    /*

    if(fIsOnline)
    {
        Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
        return ;
    }
    fLoMult=LoMult;
    fHiMult=HiMult;
    if(!MultBins)
    {
    fMultBins=fHiMult-fLoMult+1;
    }else
    {
    fMultBins=MultBins;
    }*/
}

void AliLRCProcess::SetHistMultRangeHor(Int_t LoMult,Int_t HiMult,Int_t MultBins)
{
    // Setter for multiplicity range and N bins in histos
    if(fIsOnline)
    {
        Printf("Can't change histos paramiters after InitDataMembers() was called! \n");
        return ;
    }
    fLowMultHor	= LoMult;
    fHiMultHor 	= HiMult;
    if(!MultBins)
    {
        fMultBinsHor = fHiMultHor-fLowMultHor+1;
    }else
    {
        fMultBinsHor = MultBins;
    }
}

void AliLRCProcess::SetHistMultRangeVert(Int_t LoMult,Int_t HiMult,Int_t MultBins)
{
    // Setter for multiplicity range and N bins in histos
    if(fIsOnline)
    {
        Printf("Can't change histos parameters after InitDataMembers() was called! \n");
        return ;
    }
    fLowMultVert	= LoMult;
    fHiMultVert 	= HiMult;
    if(!MultBins)
    {
        fMultBinsVert = fHiMultVert-fLowMultVert+1;
    }else
    {
        fMultBinsVert = MultBins;
    }
}

void AliLRCProcess::SetOutputSlotNumber(Int_t SlotNumber)
{
    //Sets number of output slot for LRCProcessor
    fOutputSlot=SlotNumber;
}

//________________________________________________________________________



TList*   AliLRCProcess::CreateOutput() const
{
    // Creates a link to output data TList
    return fOutList;
}

TString  AliLRCProcess::GetShortDef() const
{
    return fShortDef;
}

Int_t AliLRCProcess::GetOutputSlotNumber() const
{
    // Returns number of output slot for LRCProcessor
    return fOutputSlot;
}

void AliLRCProcess::StartEvent()
{
    // Open new Event for track by track event import
    if(fIsEventOpend)                     // Check if trying to open event more than once !
    {
        Printf("Event is already opened! Auto finishing ! \n");
        cout<<fShortDef<<": event count = "<<fEventCount<<" ";
        Printf("NchF = %i,NchB = %i \n",fNchFw,fNchBw);

        FinishEvent();
    }
    if(!fIsOnline)                        // Autocreating histos if InitDataMembers was not called by hand
    {
        Printf("InitDataMembers was not called by hand ! Autocreating histos...\n");
        InitDataMembers();
    }

    fNchFw=0;
    fSumPtFw=0;
    fNchBw=0;
    fSumPtBw=0;
    fSumPtBw2=0;
    
    fNchFwPlus = 0;
    fNchBwPlus = 0;
    fNchFwMinus = 0;
    fNchBwMinus = 0;
    
    //added 23.03
    for ( int pid = 0; pid < kSparsePIDtotal; pid++ )
    {
        fNchFwPID[pid] = 0;
        fNchFwPlusPID[pid] = 0;
        fNchFwMinusPID[pid] = 0;
        fSumPtFwPID[pid] = 0;
        fSumEtFwPID[pid] = 0;
        
        fNchBwPID[pid] = 0;
        fNchBwPlusPID[pid] = 0;
        fNchBwMinusPID[pid] = 0;
        fSumPtBwPID[pid] = 0;
        fSumEtBwPID[pid] = 0;
    }
    
    //fNchFwPtPt = 0;

    fIsEventOpend=kTRUE;
}
void AliLRCProcess::AddTrackForward(Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType  )
{
    // Imports track to the event directly to Forward window
    if(!fIsEventOpend)
    {Printf("Event is not opened!\n");
        return;}

    //Bool_t lAddDecision = kFALSE;

    //if ( fPidForward == 100 ) //add particle if we know pid (i.e. particleType != -1)
    //{
    //if ( particleType != -1 )
    //{
    //	lAddDecision = kTRUE;
    //cout << "fill fwd with pid " << particleType << endl;
    //}
    //}
    //else if ( fPidForward != -1 ) //if we specify pid for this window - just check it
    //{
    //if ( particleType == fPidForward )
    //{
    //lAddDecision = kTRUE;
    //cout << "fill fwd with pid " << particleType << endl;
    //}
    //}
    //else
    //lAddDecision = kTRUE;

    //if ( !lAddDecision )
    //return;

    fHistTestPIDForward->Fill( particleType );

    fNchFw++;
    Charge > 0 ? fNchFwPlus++ : fNchFwMinus++;
    fSumPtFw+=Pt;
    fHistPtForward->Fill(Pt);
    fHistEtaForward->Fill(Eta);
    fHistPhiForward->Fill(Phi);

    //added 15.12.12
    fHistTracksChargeForward->Fill(Charge);
    
    //added 23.03
    for ( int pid = 0; pid < kSparsePIDtotal; pid++ )
    {
        if (   pid == kSparsePIDany //write ALL pid types
               ||
               ( pid == kSparsePIDdefined && particleType != -1 ) //write defined particles
               ||
               ( fCorrespondanceWithAliROOTpid[pid] == particleType ) //write not defined particles
               )
        {
            fNchFwPID[pid]++;
            Charge > 0 ? fNchFwPlusPID[pid]++ : fNchFwMinusPID[pid]++;
            fSumPtFwPID[pid] += Pt;
            if ( pid != kSparsePIDany )
            {
                Double_t lMass = 0;//AliPID::ParticleMass( particleType );
                fSumEtFwPID[pid] += sqrt( Pt*Pt + lMass*lMass ) ;
            }
        }
    }

}
void AliLRCProcess::AddTrackBackward(Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType  )
{
    // Imports track to the event directly to Backward window
    if(!fIsEventOpend)
    {Printf("Event is not opened!\n");
        return;}

    /*Bool_t lAddDecision = kFALSE;

    if ( fPidBackward == 100 ) //add particle if we know pid (i.e. particleType != -1)
    {
        if ( particleType != -1 )
        {
            lAddDecision = kTRUE;
            //cout << "fill fwd with pid " << particleType << endl;
        }
    }
    else if ( fPidBackward != -1 ) //if we specify pid for this window - just check it
    {
        if ( particleType == fPidBackward )
        {
            lAddDecision = kTRUE;
            //cout << "fill fwd with pid " << particleType << endl;
        }
    }
    else
        lAddDecision = kTRUE;

    if ( !lAddDecision )
        return;
        */
    fHistTestPIDBackward->Fill( particleType );

    fNchBw++;
    Charge > 0 ? fNchBwPlus++ : fNchBwMinus++;
    fSumPtBw += Pt;
    fSumPtBw2 += Pt*Pt;
    fProfdPtB->Fill( 1, Pt );
    fProfdPtB->Fill( 2, Pt*Pt );
    fHistPtBackward->Fill( Pt );
    fHistEtaBackward->Fill( Eta );
    fHistPhiBackward->Fill( Phi );

    //added 15.12.12
    fHistTracksChargeBackward->Fill(Charge);

    //added 23.03
    for ( int pid = 0; pid < kSparsePIDtotal; pid++ )
    {
        if (   pid == kSparsePIDany //write ALL pid types
               ||
               ( pid == kSparsePIDdefined && particleType != -1 ) //write defined particles
               ||
               ( fCorrespondanceWithAliROOTpid[pid] == particleType )
               )
        {
            fNchBwPID[pid]++;
            Charge > 0 ? fNchBwPlusPID[pid]++ : fNchBwMinusPID[pid]++;
            fSumPtBwPID[pid] += Pt;
            if ( pid != kSparsePIDany )
            {
                Double_t lMass = 0;//AliPID::ParticleMass( particleType );
                fSumEtBwPID[pid] += sqrt( Pt*Pt + lMass*lMass ) ;
            }
            
        }
    }

}



void AliLRCProcess::AddTrackPtEta(Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType  )
{
    //cout << Pt << endl;


    //if particle type is different - ignore the track
    //if ( fWhichParticleToFill != -1 && ParticleType != fWhichParticleToFill )
    //	return;
    //Track by track event import :  Imports track to the event

    //cout << "fWhichParticleToProcess = " << fWhichParticleToProcess
    //	<< ", particleType = " << particleType << endl;
    if(!fIsEventOpend)
    {Printf("Event is not opened!\n");
        return;}

    //  Global track data
    fHistPt->Fill(Pt);
    fHistEta->Fill(Eta);

    //Bool_t lAddForwardDecision = kFALSE;
    //Bool_t lAddBackwardDecision = kFALSE;

    //Forward window
    if( ( fStartForwardETA < Eta ) && ( Eta < fEndForwardETA ) )
        if( IsPhiInRange( Phi, fStartForwardPhi, fEndForwardPhi) )//( fStartForwardPhi < Phi ) && ( Phi < fEndForwardPhi ) )
        {
            AddTrackForward( Pt, Eta, Phi, Charge, particleType );
            //if ( fPidFillCondition == kLRCpidForForwardOnly
            //	|| fPidFillCondition == kLRCpidForBoth )
            //if ( fPidForward != -1 )//kLRCany )
            //{
            //	if ( particleType == fPidForward )//kLRCpion )//particleType )
            //	{
            //		lAddForwardDecision = kTRUE;//AddTrackForward( Pt, Eta, Phi );
            //cout << "fill fwd with pid " << particleType << endl;
            //	}
            //}
            //else
            //	lAddForwardDecision = kTRUE;//AddTrackForward( Pt, Eta, Phi );
        }
    //if ( lAddForwardDecision )
    //{
    //	AddTrackForward( Pt, Eta, Phi, particleType );
    //	fHistTestPIDForward->Fill( particleType );
    //}

    //Backward window
    if( ( fStartBackwardETA < Eta ) && ( Eta < fEndBackwardETA ) )
        if (
                (
                    IsPhiInRange( Phi, fStartBackwardPhi, fEndBackwardPhi)  //( fStartBackwardPhi < Phi ) && ( Phi < fEndBackwardPhi )
                    )
                ||
                (
                    fDoubleSidedBackwardPhiWindow //if this option is true
                    && IsPhiInRange( Phi, fStartBackwardPhi + TMath::Pi(), fEndBackwardPhi + TMath::Pi() )  //
                    //&& ( fStartBackwardPhi + TMath::Pi() < Phi )
                    //&& ( Phi < fEndBackwardPhi + TMath::Pi() )
                    )
                )
        {
            AddTrackBackward( Pt, Eta, Phi, Charge, particleType );
            //if ( fPidFillCondition == kLRCpidForBackwardOnly
            //	|| fPidFillCondition == kLRCpidForBoth )
            //if ( fPidBackward != -1 )//kLRCany )
            //{
            //	if ( particleType == fPidBackward )//kLRCpion )//particleType )
            //	{
            //		lAddBackwardDecision = kTRUE;//AddTrackBackward( Pt, Eta, Phi );
            //cout << "fill bkwd with pid " << particleType << endl;
            //	}

            //}
            //else
            //	lAddBackwardDecision = kTRUE;//AddTrackBackward( Pt, Eta, Phi );
        }
    //if ( lAddBackwardDecision )
    //{
    //	AddTrackBackward( Pt, Eta, Phi, particleType );
    //	fHistTestPIDBackward->Fill( particleType );
    //}

}



void AliLRCProcess::AddTrackPtEtaMixing( Int_t winFB, Double_t Pt, Double_t Eta ,Double_t Phi, Short_t Charge, Int_t particleType  )
{
    // put track in F or B window using varible winFB
    if(!fIsEventOpend)
    {
        Printf("Event is not opened!\n");
        return;
    }

    //  Global track data
    fHistPt->Fill(Pt);
    fHistEta->Fill(Eta);


    //Forward window
    if( winFB == 0
            && ( fStartForwardETA < Eta ) && ( Eta < fEndForwardETA ) )
        if( IsPhiInRange( Phi, fStartForwardPhi, fEndForwardPhi) ) // (fStartForwardPhi < Phi ) && ( Phi < fEndForwardPhi ) )
        {
            AddTrackForward( Pt, Eta, Phi, Charge, particleType );
        }

    //Backward window
    if( winFB == 1
            && ( fStartBackwardETA < Eta ) && ( Eta < fEndBackwardETA ) )
        if (
                (
                    IsPhiInRange( Phi, fStartBackwardPhi, fEndBackwardPhi)  //( fStartBackwardPhi < Phi ) && ( Phi < fEndBackwardPhi )
                    )
                ||
                (
                    fDoubleSidedBackwardPhiWindow //if this option is true
                    && IsPhiInRange( Phi, fStartBackwardPhi + TMath::Pi(), fEndBackwardPhi + TMath::Pi() )
                    //&& ( fStartBackwardPhi + TMath::Pi() < Phi )
                    //&& ( Phi < fEndBackwardPhi + TMath::Pi() )
                    )
                )
        {
            AddTrackBackward( Pt, Eta, Phi, Charge, particleType );
        }


}

void AliLRCProcess::FinishEvent(Bool_t kDontCount)
{
    // Track by track event import : Close opened event and fill event summary histos

    if(!fIsEventOpend)
    {
        Printf("Event is not opened!\n");
        return;
    }
    if ( kDontCount ) //don't count this event! just ignore it
    {
        fIsEventOpend = kFALSE;
        return;
    }
    //fHistSparseDimensionLabeling->Fill(1);
    //Filling even-total data
    //cout << "filling" << endl;
    /*Double_t lCloudData[en_sparse_total*en_sparse_PID_total];
    lCloudData[en_sparse_N_f] = fNchFw;     //write Nf
    lCloudData[en_sparse_N_b] = fNchBw;     //write Nb
    lCloudData[en_sparse_N2_f] = fNchFw*fNchFw;     //write Nf^2
    lCloudData[en_sparse_Nf_Nb] = fNchFw*fNchBw;     //write Nb
    
    lCloudData[en_sparse_Pt_f] = 0; //fill bin 0, if don't have appropriate PtSum
    lCloudData[en_sparse_Pt_b] = 0; //fill bin 0, if don't have appropriate PtSum

    lCloudData[en_sparse_Nf_plus] = fNchFwPlus;
    lCloudData[en_sparse_Nf_minus] = fNchFwMinus;
    lCloudData[en_sparse_Nb_plus] = fNchBwPlus;
    lCloudData[en_sparse_Nb_minus] = fNchBwMinus;
    lCloudData[en_sparse_Nf_plus_Nb_minus] = fNchFwPlus * fNchBwMinus;
    lCloudData[en_sparse_Nb_plus_Nf_minus] = fNchBwPlus * fNchFwMinus; */
    
    fHistNN->Fill(fNchFw,fNchBw);

    if ( fUseAccumulatingHist )
    {

        fArrAccumulatedValues->Fill( en_arr_labels_NN_Nevents,   1                  );
        fArrAccumulatedValues->Fill( en_arr_labels_NN_Nf     ,   fNchFw             );
        fArrAccumulatedValues->Fill( en_arr_labels_NN_Nb     ,   fNchBw             );
        fArrAccumulatedValues->Fill( en_arr_labels_NN_N2_f   ,   fNchFw*fNchFw      );
        fArrAccumulatedValues->Fill( en_arr_labels_NN_Nf_Nb  ,   fNchFw*fNchBw      );
    }

    if( fNchBw != 0 )
    {
        fSumPtBw = fSumPtBw / fNchBw;
        //lCloudData[en_sparse_Pt_b] = fSumPtBw; //write <PtB>
        fProfTestLRC->Fill( fNchFw, fSumPtBw );
        fHistPtN->Fill( fNchFw, fSumPtBw );
        //cout << "fill PtN: fNchFw = " << fNchFw << ", fSumPtBw=" << fSumPtBw <<  endl;
        fProfNberr->Fill(fNchFw, 1.0 / fNchBw);

        if ( fUseAccumulatingHist )
        {
            fArrAccumulatedValues->Fill( en_arr_labels_PtN_Nevents,  1                  );
            fArrAccumulatedValues->Fill( en_arr_labels_PtN_Nf     ,  fNchFw             );
            fArrAccumulatedValues->Fill( en_arr_labels_PtN_PtB    ,  fSumPtBw           );
            fArrAccumulatedValues->Fill( en_arr_labels_PtN_N2_f   ,  fNchFw*fNchFw      );
            fArrAccumulatedValues->Fill( en_arr_labels_PtN_Ptb_Nf ,  fSumPtBw*fNchFw    );
        }

        if( fNchFw != 0 )
        {
            fSumPtFw = fSumPtFw / fNchFw;
            //lCloudData[en_sparse_Pt_f] = fSumPtFw; //write <PtF>
            fHistPtPt->Fill( fSumPtFw, fSumPtBw );
            fProfNberrPtPt->Fill( fSumPtFw, 1.0 / fNchBw );
            // dPtB for PtPt
            fProfdPtB->Fill( 15, fSumPtBw, fNchBw );
            fProfdPtB->Fill( 16, fSumPtBw2 / fNchBw, fNchBw );
            fHistNchForwardPtPt->Fill(fNchFw);

            if ( fUseAccumulatingHist )
            {
                fArrAccumulatedValues->Fill( en_arr_labels_PtPt_Nevents,   1                   );
                fArrAccumulatedValues->Fill( en_arr_labels_PtPt_PtF    ,   fSumPtFw            );
                fArrAccumulatedValues->Fill( en_arr_labels_PtPt_PtB    ,   fSumPtBw            );
                fArrAccumulatedValues->Fill( en_arr_labels_PtPt_Pt2_f  ,   fSumPtFw*fSumPtFw   );
                fArrAccumulatedValues->Fill( en_arr_labels_PtPt_Ptf_Ptb,   fSumPtBw*fSumPtFw   );
            }

        }
    }





    if ( fUseSparse )
    {
        Double_t lCloudData[kSparseTotal*kSparsePIDtotal];

        for (Int_t d = 0; d < kSparsePIDtotal; ++d)
        {
            Int_t binShift = kSparseTotal*d; //step over dimension set

            lCloudData[kSparseNf + binShift ] = fNchFwPID[d];     //write Nf
            lCloudData[kSparseNb + binShift ] = fNchBwPID[d];     //write Nb
            lCloudData[en_sparse_N2_f + binShift ] = fNchFwPID[d]*fNchFwPID[d];     //write Nf^2
            lCloudData[en_sparse_Nf_Nb + binShift ] = fNchFwPID[d]*fNchBwPID[d];     //write Nb

            lCloudData[kSparsePtF + binShift ] = 0; //fill bin 0, if don't have appropriate PtSum
            lCloudData[kSparsePtB + binShift ] = 0; //fill bin 0, if don't have appropriate PtSum

            lCloudData[en_sparse_Ptb_Nf + binShift ] = 0;
            lCloudData[en_sparse_Pt2_f + binShift ] = 0;
            lCloudData[en_sparse_Ptf_Ptb + binShift ] = 0;


            if( fNchBwPID[d] != 0 )
            {
                double lSumPtBwPID = fSumPtBwPID[d] / fNchBwPID[d];
                lCloudData[kSparsePtB + binShift ] = lSumPtBwPID; //write <PtB>

                fSumEtBwPID[d] = fSumEtBwPID[d] / fNchBwPID[d];
                //lCloudData[en_sparse_Et_b + binShift ] = fSumEtBwPID[d]; //write <PtB>
                lCloudData[en_sparse_Ptb_Nf + binShift ] = lSumPtBwPID * fNchFwPID[d];

                //fProfTestLRC->Fill( fNchFw, fSumPtBw );
                //fHistPtN->Fill( fNchFw, fSumPtBw );
                //cout << "fill PtN: fNchFw = " << fNchFw << ", fSumPtBw=" << fSumPtBw <<  endl;

                //fProfNberr->Fill(fNchFw, 1.0 / fNchBw);


                if( fNchFwPID[d] != 0 )
                {
                    double lSumPtFwPID = fSumPtFwPID[d] / fNchFwPID[d];
                    lCloudData[kSparsePtF + binShift ] = lSumPtFwPID; //write <PtF>

                    fSumEtFwPID[d] = fSumEtFwPID[d] / fNchFwPID[d];

                    lCloudData[en_sparse_Pt2_f + binShift ] = lSumPtFwPID * lSumPtFwPID;
                    lCloudData[en_sparse_Ptf_Ptb + binShift ] = lSumPtFwPID*lSumPtBwPID;

                    //lCloudData[en_sparse_Et_f + binShift ] = fSumEtFwPID[d]; //write <PtF>
                    //fHistPtPt->Fill( fSumPtFw, fSumPtBw );
                    //fProfNberrPtPt->Fill( fSumPtFw, 1.0 / fNchBw );
                    // dPtB for PtPt
                    //fProfdPtB->Fill( 15, fSumPtBw, fNchBw );
                    //fProfdPtB->Fill( 16, fSumPtBw2 / fNchBw, fNchBw );
                    //fHistNchForwardPtPt->Fill(fNchFw);

                }
            }
            /*
            lCloudData[en_sparse_Nf_plus + binShift ] = fNchFwPlusPID[d];
            lCloudData[en_sparse_Nf_minus + binShift ] = fNchFwMinusPID[d];
            lCloudData[en_sparse_Nb_plus + binShift ] = fNchBwPlusPID[d];
            lCloudData[en_sparse_Nb_minus + binShift ] = fNchBwMinusPID[d];
            lCloudData[en_sparse_Nf_plus_Nb_minus + binShift ] = fNchFwPlusPID[d] * fNchBwMinusPID[d];
            lCloudData[en_sparse_Nb_plus_Nf_minus + binShift ] = fNchBwPlusPID[d] * fNchFwMinusPID[d];
             */
        }

        //tmp (22.03): fill pid with fignya data
        /* for (Int_t d = 1; d < en_sparse_PID_total; ++d)
        {
            Int_t binShift = en_sparse_total*d;
            lCloudData[en_sparse_N_f + binShift ] = 1+d;     //write Nf
            lCloudData[en_sparse_N_b + binShift ] = 2+d;     //write Nb
            lCloudData[en_sparse_N2_f + binShift ] = 3+d;     //write Nf^2
            lCloudData[en_sparse_Nf_Nb + binShift ] = 4+d;     //write Nb

            lCloudData[en_sparse_Pt_f + binShift ] = 5+d; //fill bin 0, if don't have appropriate PtSum
            lCloudData[en_sparse_Pt_b + binShift ] = 6+d; //fill bin 0, if don't have appropriate PtSum

            lCloudData[en_sparse_Nf_plus + binShift ] = 7+d;
            lCloudData[en_sparse_Nf_minus + binShift ] = 8+d;
            lCloudData[en_sparse_Nb_plus + binShift ] = 9+d;
            lCloudData[en_sparse_Nb_minus + binShift ] = 10+d;
            lCloudData[en_sparse_Nf_plus_Nb_minus + binShift ] = 11+d;
            lCloudData[en_sparse_Nb_plus_Nf_minus + binShift ] = 12+d;
        }
    */

        //    cout << "before filling" << endl;
        fHistClouds->Fill( lCloudData ); //fill sparse hist with Nf, Nb, <PtF>, <PtB>
        //        cout << "filled." << endl;
    }




    //additional info-hist
    if ( fNchFw > 0 || fNchBw > 0 )
    {
        Double_t lAvMult = ( fNchFw + fNchBw ) / 2.;
        fHistDifferenceNf->Fill( fNchFw, ( fNchFw-fNchBw ) / lAvMult );
        fHistDifferenceNb->Fill( fNchBw, ( fNchBw-fNchFw ) / lAvMult);
    }

    //cout << "n particles: " << fNchFw << " , Back = " << fNchBw << endl;
    //cout << "fHistNN: " << fHistNN->GetEntries() <<  endl;
    //cout << "mean= " << fHistNN->GetMean() <<  endl;


    fHistNchForward->Fill(fNchFw);
    fHistNchBackward->Fill(fNchBw);

    fEventCount++;
    fIsEventOpend = kFALSE;
    
    //fill nf-centr plot
    fHistNfCentrality->Fill( fNchFw, fEventCentrality );
    
    //cout<<fShortDef<<": event count = "<<fEventCount<<" ";
    //	 Printf("NchF = %i,NchB = %i",fNchFw,fNchBw);
}

Bool_t AliLRCProcess::IsPhiInRange( Double_t phi, Double_t phiBoundMin, Double_t phiBoundMax )
{
    if ( ( phiBoundMin < phi ) && ( phi < phiBoundMax ) )
        return kTRUE;

    //when bound is more than 2pi - check phi+2pi!
    phi += 2 * TMath::Pi();
    if ( ( phiBoundMin < phi ) && ( phi < phiBoundMax ) )
        return kTRUE;

    return kFALSE; //phi not in range
}
