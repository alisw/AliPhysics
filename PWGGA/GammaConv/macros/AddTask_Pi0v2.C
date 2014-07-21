#include <exception>                                             

class AliAnalysisTaskPi0v2;
class AliV0ReaderV1;

////////////////////CURRENTLY NOT WORKING ///////////////////////////////


// Settings
Int_t nBinsPhi=6;
Int_t epselectionmask[4]={1,1,1,1};// TPC,TPCEtaGap,V0A,V0C
const Int_t nCentralityBins=5;
Double_t fCentralityBins[nCentralityBins+1]={0,5,10,20,30,40};
Double_t fInvMassRange[2]={0.0,0.3};
//const Int_t fNRadialBins=9;
//Float_t fRadialBins[fNRadialBins+1]={0,13,20,26,35,40,55,70,90,200};

Bool_t fFillQA=kTRUE;

//Bool_t fWeightMult=kFALSE;  // cut number for mult =9

AliV0ReaderV1 *fV0Reader=NULL;
AliAnalysisManager *mgr=NULL;

const Int_t numberOfCuts=24;
TString cutarray[numberOfCuts];
TString mesoncutarray[numberOfCuts];


// Standard Cuts
                                                                  01525065009000
cutarray[0] = "1080000042092970023220000000"; mesoncutarray[0] = "01522045000000";  //standard cut Pi0 PbPb 00-100

// TPC PID
cutarray[1] = "1080001042093970023220000000"; mesoncutarray[1] = "01522045000000";
cutarray[2] = "1080001042096970023220000000"; mesoncutarray[2] = "01522045000000";
cutarray[3] = "1080001042092470023220000000"; mesoncutarray[3] = "01522045000000";
cutarray[4] = "1080001042092770023220000000"; mesoncutarray[4] = "01522045000000";
cutarray[5] = "1080001042092950023220000000"; mesoncutarray[5] = "01522045000000";

// TOF PID
cutarray[6] = "1080001042092970033220000000"; mesoncutarray[6] = "01522045000000";
cutarray[7] = "1080001042092970043220000000"; mesoncutarray[7] = "01522045000000";

// Qt max
cutarray[8] = "1080001042092970024220000000"; mesoncutarray[8] = "01522045000000";
cutarray[9] = "1080001042092970022220000000"; mesoncutarray[9] = "01522045000000";

//  Chi2 Gamma
cutarray[10] = "1080001042092970023120000000"; mesoncutarray[10] = "01522045000000";
cutarray[11] = "1080001042092970023820000000"; mesoncutarray[11] = "01522045000000";
                                        // Psi Pair
cutarray[12] = "1080001042092970023210000000"; mesoncutarray[12] = "01522045000000";
cutarray[13] = "1080001042092970023230000000"; mesoncutarray[13] = "01522045000000";

//  R Cut
cutarray[14] = "1080001044092970023220000000"; mesoncutarray[14] = "01522045000000";   //5-70
cutarray[15] = "1080001045092970023220000000"; mesoncutarray[15] = "01522045000000";   //10-180
cutarray[16] = "1080001046092970023220000000"; mesoncutarray[16] = "01522045000000";   //20
cutarray[17] = "1080001047092970023220000000"; mesoncutarray[17] = "01522045000000";   //26
cutarray[18] = "1080001048092970023220000000"; mesoncutarray[18] = "01522045000000";   //35
cutarray[19] = "1080001045092970023220000000"; mesoncutarray[19] = "01522045000000";  //60

// Single Pt
cutarray[20] = "1080001042492970023220000000"; mesoncutarray[20] = "01522045000000";
cutarray[21] = "1080001042192970023220000000"; mesoncutarray[21] = "01522045000000";

// Alpha
cutarray[22] = "1080001042092970023220000000"; mesoncutarray[22] = "01022085000000";
cutarray[23] = "1080001042092970023220000000"; mesoncutarray[23] = "01022005000000";

AliAnalysisTask *AddTask_Pi0v2(Int_t harmonic=2,Bool_t IsHeavyIon=kTRUE,Bool_t doSys=kTRUE){

    // standard with task
    printf("========================================================================================\n");
    printf("Pi0v2Analysis: Initialising AliAnalysisTaskPi0v2\n");
    printf("========================================================================================\n");

    //get the current analysis manager

    mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
	Error("AddTask_dlohner_Pi0v2", "No analysis manager found.");
	return 0;
    }

    Bool_t isMC=kFALSE;
    AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>(mgr->GetMCtruthEventHandler());
    if(mcH){
	isMC=kTRUE;
    }

    // For 2011 data
    /*
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskVZEROEPSelection.C");
    AddTaskVZEROEPSelection();
    gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskEventplane.C");
    AddTaskEventplane();
    */

    TString fV0ReaderCut="";

    if(IsHeavyIon){
	fV0ReaderCut = "1080000002084001001500000000";
    }
    else{
	fV0ReaderCut = "0000000002084001001500000000";
    }
    fV0Reader=new AliV0ReaderV1(Form("PhotonPi0v%d",harmonic));
    mgr->AddTask(fV0Reader);
    ConfigV0Reader(fV0Reader,fV0ReaderCut.Data(),IsHeavyIon);
    mgr->ConnectInput(fV0Reader,  0, mgr->GetCommonInputContainer());
    fV0Reader->GetConversionCuts()->SetFillCutHistograms("V0ReaderCuts");
   
    // Setup Task

    //========= Add task to the ANALYSIS manager =====

    AliAnalysisTaskPi0v2 *task = new AliAnalysisTaskPi0v2(Form("dlohnerTask_Pi0v%d",harmonic),harmonic);
    if(doSys)SetupPi0v2(task,IsHeavyIon,isMC,Form("dlohner_Pi0v%d",harmonic),numberOfCuts);
    else SetupPi0v2(task,IsHeavyIon,isMC,Form("dlohner_Pi0v%d",harmonic),1);

    return task;
}

void ConfigV0Reader(AliV0ReaderV1 *fV0Reader,TString analysiscut="",Bool_t IsHeavyIon=kTRUE){

    fV0Reader->SetUseOwnXYZCalculation(kTRUE);

    // Set AnalysisCut Number
    AliConversionCuts *fCuts=NULL;
    if(analysiscut!=""){
	fCuts= new AliConversionCuts(analysiscut.Data(),analysiscut.Data());
	if(fCuts->InitializeCutsFromCutString(analysiscut.Data())){
	    fV0Reader->SetConversionCuts(fCuts);
	}
    }
    else{
	// Init standard cuts
	if(IsHeavyIon){fCuts=AliConversionCuts::GetStandardCuts2010PbPb();}
	else{fCuts=AliConversionCuts::GetStandardCuts2010pp();}
	fV0Reader->SetConversionCuts(fCuts);
    }
    // Initialize
    fV0Reader->Init();
}


void SetupPi0v2(AliAnalysisTaskPi0v2 *task,Bool_t IsHeavyIon,Bool_t IsMC=kFALSE,TString outputname,Int_t ncuts=0){
 
    cout<<"Settings for Task : "<<outputname.Data()<<endl;

    task->SetV0Reader(fV0Reader);
    task->SetInvMassRange(fInvMassRange);
    task->SetNBinsPhi(nBinsPhi);
    task->SetFillQA(fFillQA);
    //task->SetEPSelectionMask(epselectionmask);

    if(IsHeavyIon){
	task->SetCentralityBins(fCentralityBins,nCentralityBins);
       // task->SetWeightMultiplicity(fWeightMult);

	// Set Cuts
	if(ncuts==0)ncuts=numberOfCuts;
        if(IsMC)ncuts=1;

	AliConversionSelection **selection=new AliConversionSelection*[ncuts];
	for(Int_t ii=0;ii<ncuts;ii++){
	    cout<<"AddingCut: "<<cutarray[ii]<<" "<<mesoncutarray[ii]<<endl;
	    selection[ii]=new AliConversionSelection(cutarray[ii],mesoncutarray[ii]);
	    selection[ii]->SetInvMassRange(fInvMassRange);
	}
	task->SetCuts(selection,ncuts);
	task->SetEtaGap(1);
    }
    else{
	// no cuts defined
    }

    //SetFlattening(task,"LHC10h");

    mgr->AddTask(task);

    AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();

    AliAnalysisDataContainer *coutput0 =
	mgr->CreateContainer("dlohner_tree",TTree::Class(),AliAnalysisManager::kExchangeContainer,"dlohner_default");

    AliAnalysisDataContainer *coutput1 =
	mgr->CreateContainer(outputname.Data(), TList::Class(),
			     AliAnalysisManager::kOutputContainer,Form("%s.root",outputname.Data()));
    //connect containers
    mgr->ConnectInput  (task,  0, cinput );
    //  mgr->ConnectOutput (task,  0, coutput0);
    mgr->ConnectOutput (task,  1, coutput1);
}

void SetFlattening(AliAnalysisTaskPi0v2 *task,TString period){

    const Int_t nCent=5;

    Int_t harmonic=task->GetHarmonic();

    Int_t periodindex=task->GetPeriodIndex(period);

    if(periodindex==0){
	// TPC EP
	if(harmonic==2){
	    Double_t cc2[nCent]={0.00904396,0.00472483,0.00306154,0.00218462,0.00167447};
	    Double_t cs2[nCent]={0.00885519,0.00516223,0.00411065,0.00380145,0.00324424};
	    Double_t cc4[nCent]={-0.00110933,-0.00110521,-0.00124342,0.00104131,0.000651779};
	    Double_t cs4[nCent]={0.00163869,-0.00053565,0.000878745,-0.000563657,-0.000604021};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPC,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}

	if(harmonic==3){
	    Double_t cc2[5]={0.0116542,0.0103631,0.00897965,0.00707409,0.00605151};
	    Double_t cs2[5]={-0.0171191,-0.013024,-0.0114752,-0.0086613,-0.00706863};
	    Double_t cc4[5]={-0.000602948,0.00144836,-0.000193641,0.000108773,-0.000518333};
	    Double_t cs4[5]={-0.00164769,0.00134327,-0.00106369,7.96546e-06,-0.000261517};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPC,0,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}

	//TPC ETA A
	if(harmonic==2){
	    Double_t cc2[5]={0.00529447,0.00278029,0.00315325,0.00173634,0.000763168};
	    Double_t cs2[5]={0.00314285,0.00170173,0.00263333,0.0018509,0.00223784};
	    Double_t cc4[5]={-0.000737254,-0.00037845,-0.000492715,0.000775897,0.000768656};
	    Double_t cs4[5]={0.000347583,3.79872e-05,0.000387037,-0.000186129,0.000432698};

	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPCEtaA,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}
	if(harmonic==3){
	    Double_t cc2[5]={0.000386277,0.000119225,0.00111969,0.000534801,0.000642703};
	    Double_t cs2[5]={-0.00581604,-0.00607255,-0.00443819,-0.00268834,-0.00299961};
	    Double_t cc4[5]={0.00051635,0.00036326,-0.000221272,4.66775e-05,-3.05784e-06};
	    Double_t cs4[5]={1.43285e-05,0.000514099,0.000619339,0.00106466,0.000344196};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPCEtaA,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}

	//TPC ETA C

	if(harmonic==2){
	    Double_t cc2[5]={-0.00562282,-0.00456735,-0.00306068,-0.0027173,-0.00172432};
	    Double_t cs2[5]={0.0101804,0.00430782,0.00394715,0.00350156,0.00302749};
	    Double_t cc4[5]={0.00150831,-0.00159271,-0.000964157,0.000525894,9.93172e-05};
	    Double_t cs4[5]={0.00119279,-4.74629e-05,0.000118845,0.000278554,3.20868e-05};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPCEtaC,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}
	if(harmonic==3){
	    Double_t cc2[5]={0.0116475,0.0102385,0.00801121,0.00552336,0.00423273};
	    Double_t cs2[5]={-0.0112722,-0.00796059,-0.00683678,-0.00531097,-0.00430716};
	    Double_t cc4[5]={-0.000609051,1.36573e-08,-0.000464961,-0.000387943,-2.28363e-05};
	    Double_t cs4[5]={0.00125449,0.00168484,-0.000390491,-0.000219447,8.11997e-07};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPTPCEtaC,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}

	//V0A

	if(harmonic==2){
	    Double_t cc2[5]={0.046427,0.0105401,-0.000152992,-0.00578274,-0.0108038};
	    Double_t cs2[5]={0.00551503,0.0158159,0.00965148,0.00135414,-0.00548846};
	    Double_t cc4[5]={0.00362833,0.00170777,0.000152998,0.00223823,0.00215164};
	    Double_t cs4[5]={0.00349056,0.00142802,0.00123298,0.00207995,0.00145625};

	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPV0A,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}
	if(harmonic==3){
	    Double_t cc2[5]={-0.0057427,-0.00482728,-0.00565919,-0.000717094,-0.00933233};
	    Double_t cs2[5]={0.0306554,-0.0144675,-0.0159243,-0.0120465,-0.00814124};
	    Double_t cc4[5]={-0.002868,0.00159533,0.00754171,0.00683898,0.00689441};
	    Double_t cs4[5]={0.00083196,0.00198133,4.68307e-05,-0.00018187,-0.0014258};

	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPV0A,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}

	// V0 C
	if(harmonic==2){
	    Double_t cc2[5]={-0.00473277,-0.000371313,0.000857122,-1.54263e-05,-0.000686139};
	    Double_t cs2[5]={0.00408304,-0.00208615,-0.00149018,-0.000853616,-2.78855e-05};
	    Double_t cc4[5]={-0.00451741,-0.00399036,-0.00318784,-0.00186472,-0.00106299};
	    Double_t cs4[5]={0.00188045,-0.00713956,-0.00484254,-0.00448149,-0.00482164};
	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPV0C,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}
	if(harmonic==3){
	    Double_t cc2[5]={-0.00259141,-0.00115826,-0.000738658,-4.96667e-05,-0.000346694};
	    Double_t cs2[5]={-0.0111001,0.00258109,0.00110959,-0.000147296,-0.000199817};
	    Double_t cc4[5]={0.000968742,0.00157903,0.000206157,0.000444206,-0.00046573};
	    Double_t cs4[5]={-0.00307319,-0.0047952,-0.00412117,-0.00320344,-0.00386629};

	    task->SetFlatteningCoeff(AliAnalysisTaskPi0v2::kEPV0C,periodindex,nCent,&cc2[0],&cs2[0],&cc4[0],&cs4[0]);
	}
    }

}


