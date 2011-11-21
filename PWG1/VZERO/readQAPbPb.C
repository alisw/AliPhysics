void readQAPbPb(const char * period ="LHC11h", char * files = "list.txt"){
	gStyle->SetPalette(1);
	 TGrid::Connect("alien://");
	// gSystem->Exec("alien_find /alice/data/2011/LHC11a/* ESDs/pass1/QA*/QAresults.root > list.txt");
	// gSystem->Exec("sed '$d' < list.txt > tmplist.txt ; mv tmplist.txt list.txt");

   // TGridResult *res = gGrid->Query(Form("/alice/data/2011/%s/*",period),"ESDs/pass1/QAresults.root");
   // const Int_t nFiles = res->GetEntries();
   // if (nFiles ==0) {
   //   Error("QA","No QA files found");
   //   delete res;
   //   return;
   // }

	TArrayI valid(1);
	TArrayI runs(1);
	AliCDBManager *man = AliCDBManager::Instance();

  man->SetDefaultStorage("raw://");

Int_t minRun = 1000000;
Int_t maxRun = 0;
Int_t nValid =0;

FILE * fin = fopen(files,"r");
Int_t runNumber, nfiles=0;
while(EOF!=fscanf(fin,"%d, ",&runNumber)){
	valid.Set(nfiles+1);
	runs.Set(nfiles+1);
	valid.SetAt(kTRUE,nfiles);
	runs.SetAt(runNumber,nfiles);


 gEnv->SetValue("XNet.ConnectTimeout",10);
 gEnv->SetValue("XNet.RequestTimeout",10);
 gEnv->SetValue("XNet.MaxRedirectCount",2);
 gEnv->SetValue("XNet.ReconnectTimeout",10);
 gEnv->SetValue("XNet.FirstConnectMaxCnt",1);	

  	man->SetRun(runNumber);
	AliCDBEntry *entry2=0;
	entry2 = man->Get("GRP/GRP/Data");
	AliGRPObject* fGRPData=0;
  	if (entry2) {
      	printf("Found an AliGRPObject in GRP/GRP/Data, reading it\n");
      	fGRPData = dynamic_cast<AliGRPObject*>(entry2->GetObject());  // new GRP entry
      	entry2->SetOwner(0);
	}
	TString activeDetList(AliDAQ::ListOfTriggeredDetectors(fGRPData->GetDetectorMask()));
	TString runType(fGRPData->GetRunType());
	TString beamType(fGRPData->GetBeamType());
	TString machineMode(fGRPData->GetMachineMode());
	TString lhcState(fGRPData->GetLHCState());
	cout<<"beamType "<<beamType<<endl;
	cout<<"machineMode "<<machineMode<<endl;
	cout<<"lhcState "<<lhcState<<endl;
	
	time_t duration = fGRPData->GetTimeEnd() - fGRPData->GetTimeStart();

	if(!lhcState.Contains("STABLE BEAMS")){ // Remove no BEAM runs
		valid.SetAt(kFALSE,nfiles);
		continue;
	}
	if(!activeDetList.Contains("VZERO")){ // Remove Runs where VZERO is not active
		valid.SetAt(kFALSE,nfiles);
		continue;
	}
	if(!runType.Contains("PHYSICS")){ // Remove no PHYSICS runs
		valid.SetAt(kFALSE,nfiles);
		continue;
	}
	if(duration<600){ // Remove Runs shorter than 10 min
		valid.SetAt(kFALSE,nfiles);
		continue;
	}
	nfiles++;
	nValid++;
	
	if(runNumber>maxRun) maxRun = runNumber;
	if(runNumber<minRun) minRun = runNumber;
}

TH1F * hTimeA = new TH1F("hTimeA","BB Leading time;;Time (ns)",nValid,-0.5,nValid-0.5);
TH1F * hTimeC = new TH1F("hTimeC","BB Leading time;;Time (ns)",nValid,-0.5,nValid-0.5);
TH1F * hBB_BG = new TH1F("hBB_BG","Trigger ratio",nValid,-0.5,nValid-0.5);
TH1F * hBB_EE = new TH1F("hBB_EE","Trigger ratio",nValid,-0.5,nValid-0.5);
TH1F * hAdcA = new TH1F("hAdcA","Average Charge",nValid,-0.5,nValid-0.5);
TH1F * hAdcC = new TH1F("hAdcC","Average Charge",nValid,-0.5,nValid-0.5);
TH1F * hMultA = new TH1F("hMultA","Average number of Fired cell",nValid,-0.5,nValid-0.5);
TH1F * hMultC = new TH1F("hMultC","Average number of Fired cell",nValid,-0.5,nValid-0.5);
TH1F * hTriggerEff_CVLN = new TH1F("hTriggerEff_CVLN","CVLN / CVBN",nValid,-0.5,nValid-0.5);
TH1F * hTriggerEff_CVHN = new TH1F("hTriggerEff_CVHN","CVHN / CVBN",nValid,-0.5,nValid-0.5);
TH1F * hTriggerEff_CVHN2 = new TH1F("hTriggerEff_CVHN2","CVHN / CVLN",nValid,-0.5,nValid-0.5);
TH1F * hPMTEdges[64];
for(int i = 0; i < 64; ++i){
	hPMTEdges[i] = new  TH1F(Form("hPMTEdges%d",i),Form("Multiplicity edge Cell %d",i),nValid,-0.5,nValid-0.5);
}

int nEntries=0;
TString trigMB, trigCVLN, trigCVHN;
for(int ifile = nfiles-1; ifile >= 0; ifile--){

	if(!valid.At(ifile)) continue;
	runNumber = runs.At(ifile);
	
   	if(runNumber>168171){
		trigMB   = "CPBI2_B1";
		trigCVLN = "CVLN_R1";
		trigCVHN = "CVHN_R2";
	} else if(runNumber>167693){
		trigMB   = "CPBI2_B1";
		trigCVLN = "CVLN_B2";
		trigCVHN = "CVHN_R2";
	} else if(runNumber>166532)	{
		trigMB   = "CPBI1";
		trigCVLN = "CVLN";
		trigCVHN = "CVHN";
	} else {
		trigMB   = "CPBI1";
		trigCVLN = "CVLN";
		trigCVHN = "CVHN";
	}
		
	TGridResult *res;
	if(strcmp(period,"LHC11h")==0){
		res = gGrid->Query(Form("/alice/data/2011/%s/000%d/*",period,runNumber),"ESDs/pass1_HLT/QAresults.root");
	}


if (res->GetEntries() ==0) {
     Error("QA",Form("No QA files found for run %d\n",runNumber));
     delete res;
	 nEntries++;
     continue;
   }

	   man->SetRun(runNumber);
 	
	AliCDBEntry *entryCTP = man->Get("GRP/CTP/Config");
   AliTriggerConfiguration *configCTP = (AliTriggerConfiguration*)entryCTP->GetObject();
	TObjArray  inputsArray = configCTP->GetInputs();
	
	Double_t rnd1=1., rnd2=1., bc1=1., bc2=1.;
	AliTriggerInput * input;
	input = (AliTriggerInput*)(inputsArray.FindObject("RND1"));
	if(input) rnd1 =  (input->GetSignature())/(double)(0x7fffffff );
	if(input)cout<<Form("RND1 = %d",input->GetSignature())<<endl;

	input = (AliTriggerInput*)(inputsArray.FindObject("RND2"));
	if(input) rnd2 =  (input->GetSignature())/(double)(0x7fffffff );
	if(input)cout<<Form("RND2 = %d",input->GetSignature())<<endl;
	
	input = (AliTriggerInput*)(inputsArray.FindObject("BC1"));
	if(input) bc1 =  1./(input->GetSignature()+1.);
	if(input)cout<<Form("BC1 = %d",input->GetSignature())<<endl;
	
	input = (AliTriggerInput*)(inputsArray.FindObject("BC2"));
	if(input) bc2 =  1./(input->GetSignature()+1.);
	if(input)cout<<Form("BC2 = %d",input->GetSignature())<<endl;


	
   	TString filename = res->GetKey(0, "turl");
   	TObjArray* tmp = filename.Tokenize("/");
	if(filename == "") continue;
   	TFile *fQA = TFile::Open(filename.Data());
   	if (!fQA) {
     	Error("QA",Form("Can not open QA file found for run %d\n",runNumber));
		nEntries++;
     	continue;
   	}
	TList *list = (TList*)fQA->Get("VZERO_PbPb_Performance/PbPbVZEROHists");
	TH2F *hTriggerDecision = (TH2F*)list->FindObject("hTriggerDecision");
	TH1F *hAdcNoTimeA = (TH1F*)list->FindObject("hAdcNoTimeV0A");
	TH1F *hAdcWithTimeA = (TH1F*)list->FindObject("hAdcWithTimeV0A");
	TH1F *hAdcNoTimeC = (TH1F*)list->FindObject("hAdcNoTimeV0C");
	TH1F *hAdcWithTimeC = (TH1F*)list->FindObject("hAdcWithTimeV0C");
	TH2F *hadcpmtwithtime = (TH2F*)list->FindObject("hadcpmtwithtime");	
	TH1F *htimepmtA = (TH1F*)list->FindObject("htimepmtV0A");
	TH1F *htimepmtC = (TH1F*)list->FindObject("htimepmtV0C");
	TH1F *hwidthA = (TH1F*)list->FindObject("hwidthV0A");
	TH1F *hwidthC = (TH1F*)list->FindObject("hwidthV0C");
	TH1F *hV0ampl = (TH1F*)list->FindObject("hV0ampl");
	TH2F *htimepmt = (TH2F*)list->FindObject("htimepmt");	
	TH2F *hwidthpmt = (TH2F*)list->FindObject("hwidthpmt");	
	TH2F *hadcwidthA = (TH2F*)list->FindObject("hadcwidthV0A");	
	TH2F *hadcwidthC = (TH2F*)list->FindObject("hadcwidthV0C");	
	TH2F *hAdcTimeA = (TH2F*)list->FindObject("hAdcTimeV0A");	
	TH2F *hAdcTimeC = (TH2F*)list->FindObject("hAdcTimeV0C");	
	TH2F *htimecorr = (TH2F*)list->FindObject("htimecorr");	
	TH2F *hNFlags   = (TH2F*)list->FindObject("hNFlags");
	TH1F *hV0A = (TH1F*)hNFlags->ProjectionX("hV0A",1,hNFlags->GetNbinsY());
	TH1F *hV0C = (TH1F*)hNFlags->ProjectionY("hV0C",1,hNFlags->GetNbinsX());
	TH2F* hVtxXYBB  =(TH2F*) list->FindObject("fhVtxXYBB");
	TH1F* hVtxZBB   =(TH1F*) list->FindObject("fhVtxZBB");
	TH2F* hVtxXYBGA =(TH2F*) list->FindObject("fhVtxXYBGA");
	TH1F* hVtxZBGA  =(TH1F*) list->FindObject("fhVtxZBGA");
	TH2F* hVtxXYBGC =(TH2F*) list->FindObject("fhVtxXYBGC");
	TH1F* hVtxZBGC  =(TH1F*) list->FindObject("fhVtxZBGC");
	
	TH2F* hRecoMult = (TH2F*) list->FindObject(Form("hRecoMult_%s-",trigMB.Data()));
	TH2F* hRecoMultPMT = (TH2F*) list->FindObject(Form("hRecoMultPMT_%s-",trigMB.Data()));
	TH1F* hTotRecoMult = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigMB.Data()));
	TH1F* hTotRecoMult_CVLN = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigCVLN.Data()));
	TH1F* hTotRecoMult_CVHN = (TH1F*) list->FindObject(Form("hTotRecoMult_%s-",trigCVHN.Data()));
	TH2F* hEqualizedMult = (TH2F*) list->FindObject(Form("hEqualizedMult_%s-",trigMB.Data()));
	
	Double_t BB  = hTriggerDecision->GetBinContent(2,2);
	Double_t EE  = hTriggerDecision->GetBinContent(1,1);
	Double_t BGA = hTriggerDecision->GetBinContent(3,2);
	Double_t BGC = hTriggerDecision->GetBinContent(2,3);
		
	hTimeA->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hTimeC->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hBB_BG->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hBB_EE->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hAdcA->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hAdcC->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hMultA->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hMultC->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hTriggerEff_CVLN->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hTriggerEff_CVHN->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
	hTriggerEff_CVHN2->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
    for(int i = 0; i < 64; ++i) {
        hPMTEdges[i]->GetXaxis()->SetBinLabel(nEntries+1,Form("%d",runNumber));
    }
	
	if(hAdcWithTimeA->GetEntries()==0) {
		delete fQA;	
		nEntries++;		
		continue;
	}

	Double_t cVLN = hTotRecoMult_CVLN->GetEntries();
	Double_t cVHN = hTotRecoMult_CVHN->GetEntries();
	Double_t cVBN = hTotRecoMult->GetEntries();
	
	// Correct for downscaling factor
	Double_t scaleVBN = 1., scaleVLN = 1., scaleVHN = 1.;
	if(trigMB.Contains("B1")) {
		if(bc1>0.) scaleVBN=1./bc1;
	} else 	if(trigMB.Contains("B2")) {
		if(bc2>0.) scaleVBN=1./bc2;
	} else if(trigMB.Contains("R1")) {
		if(rnd1>0.) scaleVBN=1./rnd1;
	} else if(trigMB.Contains("R2")) {
		if(rnd2>0.) scaleVBN=1./rnd2;
	}
	if(trigCVLN.Contains("B1")) {
		if(bc1>0.) scaleVLN=1./bc1;
	} else 	if(trigCVLN.Contains("B2")) {
		if(bc2>0.) scaleVLN=1./bc2;
	} else if(trigCVLN.Contains("R1")) {
		if(rnd1>0.) scaleVLN=1./rnd1;
	} else if(trigCVLN.Contains("R2")) {
		if(rnd2>0.) scaleVLN=1./rnd2;
	}
	
	if(trigCVHN.Contains("B1")) {
		if(bc1>0.) scaleVHN=1./bc1;
	} else 	if(trigCVHN.Contains("B2")) {
		if(bc2>0.) scaleVHN=1./bc2;
	} else if(trigCVHN.Contains("R1")) {
		if(rnd1>0.) scaleVHN=1./rnd1;
	} else if(trigCVHN.Contains("R2")) {
		if(rnd2>0.) scaleVHN=1./rnd2;
	}

	cout<<Form("CVBN = %lf \nCVLN = %lf \nCVHN = %lf\n",cVBN,cVLN,cVHN);
	hTotRecoMult->Scale(scaleVBN);
	hTotRecoMult_CVLN->Scale(scaleVLN);
	hTotRecoMult_CVHN->Scale(scaleVHN);

	cVBN *= scaleVBN/100.;
	cVLN *= scaleVLN;
	cVHN *= scaleVHN;
	
	cout<<Form("CVBN = %lf \nCVLN = %lf \nCVHN = %lf\n",cVBN,cVLN,cVHN);

	if(cVBN >0.){
		hTriggerEff_CVLN->SetBinContent(nEntries+1,cVBN/cVLN);
		if(cVLN>0.) hTriggerEff_CVLN->SetBinError(nEntries+1,cVBN/cVLN*(TMath::Sqrt(1./cVLN+1./cVBN)));
		hTriggerEff_CVHN->SetBinContent(nEntries+1,cVBN/cVHN);
		if(cVHN>0.) hTriggerEff_CVHN->SetBinError(nEntries+1,cVBN/cVHN*(TMath::Sqrt(1./(cVHN/scaleVHN)+1./(cVBN/scaleVBN*100.))));
		if(cVLN>0.) hTriggerEff_CVHN2->SetBinContent(nEntries+1,cVLN/cVHN);
		if(cVHN>0. && cVLN>0.) hTriggerEff_CVHN2->SetBinError(nEntries+1,cVLN/cVHN*(TMath::Sqrt(1./(cVHN/scaleVHN)+1./(cVLN/scaleVLN))));
	}

    Double_t beta1[64], beta2[64];
    Double_t q = 1. - 1.e-4;
    Double_t q2 = 1. - 2.e-4;
    for(int i = 0; i < 64; ++i) {
    	((TH1D*)hRecoMultPMT->ProjectionY(Form("hRecoMultPMT%d",i),i+1,i+1))->GetQuantiles(1,&beta1[i],&q);
    	((TH1D*)hRecoMultPMT->ProjectionY(Form("hRecoMultPMT%d",i),i+1,i+1))->GetQuantiles(1,&beta2[i],&q2);
        hPMTEdges[i]->SetBinContent(nEntries+1,(beta1[i]+beta2[i])/2.);
		hPMTEdges[i]->SetBinError(nEntries+1,(beta1[i] - beta2[i]));
    }

    Double_t betaSide1[2];
    Double_t betaSide2[2];
    for(int i = 0; i < 2; ++i) {
		if(i==0){
			((TH1D*)hRecoMult->ProjectionY(Form("hRecoMult1%d",i)))->GetQuantiles(1,&betaSide1[i],&q);
			((TH1D*)hRecoMult->ProjectionY(Form("hRecoMult2%d",i)))->GetQuantiles(1,&betaSide2[i],&q2);
	    }else{
			((TH1D*)hRecoMult->ProjectionX(Form("hRecoMult1%d",i)))->GetQuantiles(1,&betaSide1[i],&q);
			((TH1D*)hRecoMult->ProjectionX(Form("hRecoMult2%d",i)))->GetQuantiles(1,&betaSide2[i],&q2);
		}
    }
	hAdcA->SetBinContent(nEntries+1,(betaSide1[1] + betaSide2[1])/2.);
	hAdcC->SetBinContent(nEntries+1,(betaSide1[0] + betaSide2[0])/2.);
	hAdcA->SetBinError(nEntries+1,(betaSide1[1] - betaSide2[1]));
	hAdcC->SetBinError(nEntries+1,(betaSide1[0] - betaSide2[0]));
	
	

	
		TSpectrum s;
		int nPeaksFound;
		float * peaks;
		float max;
		float shiftA = 8.;

		nPeaksFound = s.Search(htimepmtA); peaks = s.GetPositionX(); max = -25.;
		for(int i=0;i<nPeaksFound;i++) if(peaks[i]>max) max = peaks[i];	
		htimepmtA->Fit("gaus","","",max-1.,max+1.);
		hTimeA->Fill(nEntries,gaus->GetParameter(1)-shiftA);

		nPeaksFound = s.Search(htimepmtC); peaks = s.GetPositionX(); max = -25.;
		for(int i=0;i<nPeaksFound;i++) if(peaks[i]>max) max = peaks[i];	
		htimepmtC->Fit("gaus","","",max-1.,max+1.);
		hTimeC->Fill(nEntries,gaus->GetParameter(1));

		if(BB) {
			hBB_BG->SetBinContent(nEntries,(BGA+BGC)/BB);
			hBB_EE->SetBinContent(nEntries,EE/BB);
		} else {
			hBB_BG->SetBinContent(nEntries,0);
			hBB_EE->SetBinContent(nEntries,0);
		}

		hMultA->SetBinContent(nEntries,hV0A->GetMean());
		hMultC->SetBinContent(nEntries,hV0C->GetMean());

		hAdcA->SetBinContent(nEntries,hAdcWithTimeA->GetMean());
		hAdcC->SetBinContent(nEntries,hAdcWithTimeC->GetMean());

//-------------
	TCanvas * cOut = new TCanvas("cOut",Form("Run %d",runNumber));
	cOut->Divide(2,2);
	cOut->cd(1); cOut->GetPad(1)->SetLogy();
	hAdcNoTimeA->Draw("l");
	hAdcWithTimeA->Draw("same"); hAdcWithTimeA->SetLineColor(2);
	
	cOut->cd(2); cOut->GetPad(2)->SetLogy();
	hAdcNoTimeC->Draw("l");
	hAdcWithTimeC->Draw("same"); hAdcWithTimeC->SetLineColor(2);
	
	cOut->cd(3); cOut->GetPad(3)->SetLogz();
	hadcpmtwithtime->Draw("colz");


	cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz(1);
	hEqualizedMult->Draw("colz");  
	
	cOut->Print(Form("QA_Run_%d.pdf(",runNumber));
//-------------

	cOut->Clear();
	cOut->cd(1); cOut->GetPad(1)->SetLogy();
	htimepmtA->GetXaxis()->SetRangeUser(-25.,25.); htimepmtA->Draw();

	cOut->cd(2); cOut->GetPad(2)->SetLogy();
	htimepmtC->GetXaxis()->SetRangeUser(-25.,25.); htimepmtC->Draw();
	
	cOut->cd(3); cOut->GetPad(3)->SetLogy();cOut->GetPad(3)->SetLogz(0);
	//hwidthA->GetXaxis()->SetRangeUser(0.,50.); 
	hwidthA->Draw();
	
	cOut->cd(4); cOut->GetPad(4)->SetLogy();cOut->GetPad(4)->SetLogz(0);
	//hwidthC->GetXaxis()->SetRangeUser(0.,50.); 
	hwidthC->Draw();

	cOut->Print(Form("QA_Run_%d.pdf",runNumber));
//-------------

	cOut->Clear();
	cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
	htimepmt->Draw("colz");

	cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
	//hwidthpmt->GetYaxis()->SetRangeUser(0.,50.); 
	hwidthpmt->Draw("colz");
	
	cOut->cd(3); cOut->GetPad(3)->SetLogy(0);cOut->GetPad(3)->SetLogz();
	//hadcwidthA->GetYaxis()->SetRangeUser(0.,50.); 
	hadcwidthA->Draw("colz");
	
	cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz();
	//hadcwidthC->GetYaxis()->SetRangeUser(0.,50.); 
	hadcwidthC->Draw("colz");

	cOut->Print(Form("QA_Run_%d.pdf",runNumber));
//-------------

	cOut->Clear();
	cOut->cd(1); cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz();
	hAdcTimeA->Draw("colz");

	cOut->cd(2); cOut->GetPad(2)->SetLogy(0);cOut->GetPad(2)->SetLogz();
	hAdcTimeC->Draw("colz");
	
	cOut->cd(3); cOut->GetPad(3)->SetLogz();
	hTriggerDecision->Draw("text");
	
	cOut->cd(4); cOut->GetPad(4)->SetLogz(); cOut->GetPad(4)->SetLogz(0);
	htimecorr->Draw("colz");

	cOut->Print(Form("QA_Run_%d.pdf",runNumber));
//-------------

	cOut->Clear();
	cOut->cd(1);  cOut->GetPad(1)->SetLogy(1);cOut->GetPad(1)->SetLogz(0);
	hV0A->GetXaxis()->SetRangeUser(0.,33.);hV0A->Draw();

	cOut->cd(2); cOut->GetPad(2)->SetLogy(1);cOut->GetPad(2)->SetLogz(0);
	hV0C->GetXaxis()->SetRangeUser(0.,33.);hV0C->Draw();

	cOut->cd(3); cOut->GetPad(3)->SetLogy(0);cOut->GetPad(3)->SetLogz(0); cOut->GetPad(3)->SetGridy(1);
	hTotRecoMult->Sumw2();
	hTotRecoMult_CVLN->Sumw2();
	hTotRecoMult_CVHN->Sumw2();
	hTotRecoMult_CVLN->Divide(hTotRecoMult);
	hTotRecoMult_CVHN->Divide(hTotRecoMult);
	hTotRecoMult_CVLN->Draw("e"); hTotRecoMult_CVLN->SetLineColor(4);hTotRecoMult_CVLN->SetMarkerStyle(0);hTotRecoMult_CVLN->SetMaximum(1.2);
	hTotRecoMult_CVLN->SetTitle("Multiplicity efficiency CVLN (blue) and CHVN (red)");
	hTotRecoMult_CVHN->Draw("esame"); hTotRecoMult_CVHN->SetLineColor(2);hTotRecoMult_CVHN->SetMarkerStyle(0);
	
	cOut->cd(4); cOut->GetPad(4)->SetLogy(0);cOut->GetPad(4)->SetLogz(1); cOut->GetPad(4)->SetGridx(1); cOut->GetPad(4)->SetGridy(1);
	hRecoMult->Draw("colz");
	
	
	cOut->Print(Form("QA_Run_%d.pdf",runNumber));
//-------------

	cOut->Clear();
	cOut->Divide(2,3);

	cOut->cd(1);  cOut->GetPad(1)->SetLogy(0);cOut->GetPad(1)->SetLogz(1);
	hVtxXYBB->Draw("colz");

	cOut->cd(2); cOut->GetPad(2)->SetLogy(1);cOut->GetPad(2)->SetLogz(0);
	hVtxZBB->Draw();
	
	cOut->cd(3); cOut->GetPad(3)->SetLogy(0); cOut->GetPad(3)->SetLogz(1);
	hVtxXYBGA->Draw("colz");
	
	cOut->cd(4); cOut->GetPad(4)->SetLogy(); cOut->GetPad(3)->SetLogz(0);
	hVtxZBGA->Draw();
	
	cOut->cd(5); cOut->GetPad(5)->SetLogy(0); cOut->GetPad(5)->SetLogz(1);
	hVtxXYBGC->Draw("colz");
	
	cOut->cd(6); cOut->GetPad(6)->SetLogy(); cOut->GetPad(6)->SetLogz(0);
	hVtxZBGC->Draw();

	cOut->Print(Form("QA_Run_%d.pdf)",runNumber));
//-------------

			delete cOut;
			delete fQA;			

			nEntries++;

		}
		delete res;

gStyle->SetOptStat(0);
hTimeA->SetMarkerStyle(20);
hTimeA->SetMarkerColor(2);

hTimeC->SetMarkerStyle(20);
hTimeC->SetMarkerColor(4);

TCanvas * c = new TCanvas("c","Leading time versus run number");
c->SetGridy();
hTimeA->Draw("P");
hTimeA->SetMinimum(TMath::Min(hTimeA->GetMinimum(),hTimeC->GetMinimum())-1.);
hTimeA->SetMaximum(TMath::Max(hTimeA->GetMaximum(),hTimeC->GetMaximum())+1.);
// hTimeA->GetXaxis()->SetLabelOptions("v");
// hTimeC->GetXaxis()->SetLabelOptions("v");

hTimeC->Draw("Psame");
TLegend * lg = new TLegend(0.8,0.9,1,1);
lg->AddEntry(hTimeA,"V0A - 8 ns","p");
lg->AddEntry(hTimeC,"V0C","p");
lg->Draw("same");
// TPave * pavA = new TPave(-0.5,TMath::Max(hTimeA->GetMinimum(),1.5-shiftA),nValid-0.5,TMath::Min(hTimeA->GetMaximum(),33.5-shiftA),0);
// pavA->SetFillStyle(3004);
// pavA->SetFillColor(2);
// TPave * pavC = new TPave(-0.5,TMath::Max(hTimeA->GetMinimum(),0.5),nValid-0.5,TMath::Min(hTimeA->GetMaximum(),25.5),0);
// pavC->SetFillStyle(3005);
// pavC->SetFillColor(4);
// 
// pavA->Draw("same");
// pavC->Draw("same");
TFile * fout = TFile::Open(Form("QA_Resume_%d_%d.root",minRun,maxRun),"RECREATE");

c->Print(Form("QA_Resume_%d_%d.pdf(",minRun,maxRun));
c->Write();

TCanvas * c2 = new TCanvas("c2","Trigger ratios");
c2->SetGridy();

hBB_BG->SetMarkerStyle(20);
hBB_BG->SetMarkerColor(2);
hBB_EE->SetMarkerStyle(20);
hBB_EE->SetMarkerColor(4);

hBB_BG->Draw("P");hBB_BG->SetTitle("Beam-Gas / Beam-Beam");
// hBB_EE->Draw("Psame");
// TLegend * lg2 = new TLegend(0.8,0.9,1,1);
// lg2->AddEntry(hBB_BG,"BG / BB","p");
// lg2->AddEntry(hBB_EE,"EE / BB","p");
// lg2->Draw("same");

c2->Print(Form("QA_Resume_%d_%d.pdf",minRun,maxRun));
c2->Write();


TCanvas * c3 = new TCanvas("c3","Multiplicity Eddes V0A/V0C");
c3->SetGridy();

hAdcA->SetMarkerStyle(20);
hAdcA->SetMarkerColor(2);
hAdcC->SetMarkerStyle(20);
hAdcC->SetMarkerColor(4);
hAdcA->SetMinimum(0);
hAdcA->SetMaximum(100);

hAdcA->Draw("P");
hAdcC->Draw("Psame");
TLegend * lg3 = new TLegend(0.8,0.9,1,1);
lg3->AddEntry(hAdcA,"V0A","p");
lg3->AddEntry(hAdcC,"V0C","p");
lg3->Draw("same");

c3->Print(Form("QA_Resume_%d_%d.pdf",minRun,maxRun));
c3->Write();


TCanvas * c4 = new TCanvas("c4","Average number of cell");
c4->SetGridy();

hMultA->SetMarkerStyle(20);
hMultA->SetMarkerColor(2);
hMultC->SetMarkerStyle(20);
hMultC->SetMarkerColor(4);
hMultA->SetMinimum(0);
hMultA->SetMaximum(32);

hMultA->Draw("P");
hMultC->Draw("Psame");
TLegend * lg4 = new TLegend(0.8,0.9,1,1);
lg4->AddEntry(hMultA,"V0A","p");
lg4->AddEntry(hMultC,"V0C","p");
lg4->Draw("same");

c4->Print(Form("QA_Resume_%d_%d.pdf",minRun,maxRun));
c4->Write();

TCanvas * c5 = new TCanvas("c5","Trigger Efficiency");
c5->SetGridy();

hTriggerEff_CVLN->SetMarkerStyle(20);
hTriggerEff_CVLN->SetMarkerColor(2);
hTriggerEff_CVHN->SetMarkerStyle(20);
hTriggerEff_CVHN->SetMarkerColor(4);
hTriggerEff_CVHN2->SetMarkerStyle(4);
hTriggerEff_CVHN2->SetMarkerColor(1);
hTriggerEff_CVLN->SetMinimum(0);
hTriggerEff_CVLN->SetMaximum(15.);
hTriggerEff_CVLN->SetTitle("Centrality triggers fraction");

hTriggerEff_CVLN->Draw("P");
hTriggerEff_CVHN->Draw("Psame");
hTriggerEff_CVHN2->Draw("Psame");
TLegend * lg5 = new TLegend(0.7,0.8,1,1);
lg5->AddEntry(hTriggerEff_CVLN,"(CPBI2/100) / CVLN","p");
lg5->AddEntry(hTriggerEff_CVHN,"(CPBI2/100) / CVHN","p");
lg5->AddEntry(hTriggerEff_CVHN2,"CVLN / CVHN","p");
lg5->Draw("same");

c5->Print(Form("QA_Resume_%d_%d.pdf",minRun,maxRun));
c5->Write();


TCanvas * cedge[8];
for(int i = 0; i < 8; ++i){
	cedge[i] = new TCanvas(Form("cedge%d",i),Form("Edge Ring %d",i));
	cedge[i]->SetGridy();
	cedge[i]->Divide(3,3);
	for(int iCh = 0; iCh < 8; ++iCh){
		cedge[i]->cd(iCh+1);
		hPMTEdges[iCh+i*8]->SetMarkerStyle(20);
		hPMTEdges[iCh+i*8]->Draw("P");
	}
	if(i==7) cedge[i]->Print(Form("QA_Resume_%d_%d.pdf)",minRun,maxRun));
	else cedge[i]->Print(Form("QA_Resume_%d_%d.pdf",minRun,maxRun));
	cedge[i]->Write();
}





fout->Close();

gSystem->Exec(Form("tar cvf QA_Runs_%d_%d.tar QA_Run*.pdf QA_Resume_%d_%d.pdf QA_Resume_%d_%d.root",minRun,maxRun,minRun,maxRun,minRun,maxRun));
gSystem->Exec(Form("rm -f QA_Run*.pdf QA_Resume_%d_%d.pdf  QA_Resume_%d_%d.root",minRun,maxRun,minRun,maxRun));


}

