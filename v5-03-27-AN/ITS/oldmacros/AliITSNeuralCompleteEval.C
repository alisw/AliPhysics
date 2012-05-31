// Neural performance evaluation
//
// before using this macro, you should have make run the macro
// 'ITSstoreFindableTracks.C' in order to have a file named following
// the rule to append '_fnd' to the name 
// (e.g.   galice.root --> galice_fnd.root)
// and have in it a tree with some tracks statistics useful for evaluation
// of performance with the same event without loosing too much time...

void AliITSNeuralCompleteEval
(Int_t nsecs, const char *filename, Bool_t low = kFALSE, 
 Bool_t draw = kTRUE, const char *save = 0)
{
	Int_t N, M;
	Float_t *edge = 0;
	if (!low) {
		N = 7;
		edge = new Float_t[8];
		edge[0] = 0.0;
		edge[1] = 0.5;
		edge[2] = 1.0;
		edge[3] = 1.5;
		edge[4] = 2.0;
		edge[5] = 3.0;
		edge[6] = 4.0;
		edge[7] = 5.0;
	}
	else {
		N = 5;
		edge = new Float_t[6];
		edge[0] = 0.0;
		edge[1] = 0.2;
		edge[2] = 0.4;
		edge[3] = 0.6;
		edge[4] = 0.8;
		edge[5] = 1.0;
	}
	
	Float_t find[2] = {0.0, 0.0}, find1[2] = {0.0, 0.0};
	Float_t good[2] = {0.0, 0.0}, good1[2] = {0.0, 0.0};
	Float_t fake[2] = {0.0, 0.0}, fake1[2] = {0.0, 0.0};
	
	// histos filled with neural results
	TH1D *hgood[2], *hfake[2], *hfind[2];
	TH1D *hgood1[2], *hfake1[2], *hfind1[2];
	TH1D *hg[2], *hf[2];
	if (draw) {
		hgood[0] = new TH1D("hgood0", "Good found tracks", N, edge);
		hfake[0] = new TH1D("hfake0", "Fake found tracks", N, edge);
		hfind[0] = new TH1D("hfound0", "Findable tracks", N, edge);
		hgood[1] = new TH1D("hgood1", "Good found tracks", N, edge);
		hfake[1] = new TH1D("hfake1", "Fake found tracks", N, edge);
		hfind[1] = new TH1D("hfound1", "Findable tracks", N, edge);
	
		hgood[0]->Sumw2(); 
		hfake[0]->Sumw2(); 
		hfind[0]->Sumw2(); 
		hgood[1]->Sumw2(); 
		hfake[1]->Sumw2(); 
		hfind[1]->Sumw2(); 
		
		// histos for evaluating percentual efficiency
		hg[0] = new TH1D("hg0", "Efficiency (%) for #geq 5 right pts.", N, edge);
		hf[0] = new TH1D("hf0", "Fake probability (%) for #geq 5 right pts.", N, edge);
		hg[1] = new TH1D("hg1", "Efficiency (%) for 6 right pts.", N, edge);
		hf[1] = new TH1D("hf1", "Fake probability (%) for 6 right pts.", N, edge);
	}
	
	TFile *ffind;
	TTree *tree;
	
	if (low) {
		ffind = new TFile(Form("%s.root", filename), "READ");
		tree  = (TTree*)ffind->Get("TreeF");
	}
	else {
		ffind = new TFile(Form("%s_fnd.root", filename), "READ");
		tree  = (TTree*)ffind->Get("tree");
	}
	
	TFile *ftracks = new TFile(Form("%s_%d.root", filename, nsecs), "READ");
	
	
	Double_t pt;
	Int_t i, j, count, prim, *none = 0, div;
	Int_t entries = tree->GetEntries(), label, min[] = {5,6};
	tree->SetBranchAddress("pt", &pt);
	tree->SetBranchAddress("count", &count);
	tree->SetBranchAddress("prim", &prim);
	
	AliITSneuralTrack *trk = 0;
	div = low ? 50 : 500;
	for(i = 1;;i++) {
		trk = (AliITSneuralTrack*)ftracks->Get(Form("AliITSneuralTrack;%d", i));
		if (i%div == 0) cout << "\rEvaluating found track " << i << flush;
		if (!trk) break;
		for (j = 0; j < 2; j++) {
			label = trk->EvaluateTrack(0, min[j], none);
			tree->GetEntry(abs(label));
			if (count >= min[j]) {
				if (label > 0) {
					if (draw) hgood[j]->Fill(pt);
					good[j]++;
					if (pt >= 1.0) good1[j]++;
				}
				else {
					if (draw) hfake[j]->Fill(pt);
					fake[j]++;
					if (pt >= 1.0) fake1[j]++;
				}
			}
		}
	}
	cout << endl;
		
	div = low ? 200 : 20000;
	
	for (i = 0; i < entries; i++) {
		if (i%div == 0) cout << "\rEvaluating findable track no. " << i << flush;
		tree->GetEntry(i);
		for (j = 0; j < 2; j++) {
			if (count >= min[j]) {
				find[j]++;
				if (draw) hfind[j]->Fill(pt);
				if (pt >= 1.0) find1[j]++;
			}
		}
	}
	cout << endl;
	cout << hgood[0]->GetEntries() << " " << hgood[1]->GetEntries() << endl;
	cout << hfake[0]->GetEntries() << " " << hfake[1]->GetEntries() << endl;
	cout << hfind[0]->GetEntries() << " " << hfind[1]->GetEntries() << endl << endl;
	
	if (draw) {
		TCanvas *canv[2];
		canv[0] = new TCanvas("c_0", "Tracking efficiency (soft)", 0, 0, 600, 500);
		canv[1] = new TCanvas("c_1", "Tracking efficiency (hard)", 630, 0, 600, 500);
		
		TLine *line1 = new TLine(1,100.0,edge[N],100.0); line1->SetLineStyle(4);
		TLine *line2 = new TLine(1,90,edge[N],90); line2->SetLineStyle(4);
			
		Bool_t good_drawn;
		for (i = 0; i < 2; i++) {
			canv[i]->cd();
			good_drawn = kFALSE;
			if (hgood[i]->GetEntries() > 0.0) {
				good_drawn = kTRUE;
				hg[i]->Divide(hgood[i], hfind[i], 100.0, 1.0);
				hg[i]->SetMaximum(120);
				hg[i]->SetMinimum(0);
				hg[i]->SetMarkerStyle(21);
				hg[i]->SetMarkerSize(1);
				hg[i]->SetStats(kFALSE);
				hg[i]->GetXaxis()->SetTitle("pt (GeV/c)");
				hg[i]->Draw("PE1");
			}
			if (hfake[i]->GetEntries() > 0.0) {
				hf[i]->Divide(hfake[i], hfind[i], 100.0, 1.0);
				hf[i]->SetMaximum(120);
				hf[i]->SetMinimum(0);
				hf[i]->SetMarkerStyle(25);
				hf[i]->SetMarkerSize(1);
				hf[i]->SetStats(kFALSE);
				if (good_drawn)
					hf[i]->Draw("PE1SAME");
				else
					hf[i]->Draw("PE1");
			}
			line1->Draw("histosame");
			line2->Draw("histosame");
			canv[i]->Update();
		}
		canv[0]->SaveAs(Form("%s_soft.eps", filename));
		canv[1]->SaveAs(Form("%s_hard.eps", filename));
		cout << endl;
	}
	
	Float_t sgood[2] = {0.0, 0.0}, sgood1[2] = {0.0, 0.0};
	Float_t sfake[2] = {0.0, 0.0}, sfake1[2] = {0.0, 0.0};
	for (i = 0; i < 2; i++) {
		sgood[i] = error(good[i], find[i]);
		sgood1[i] = error(good1[i], find1[i]);
		sfake[i] = error(fake[i], find[i]);
		sfake1[i] = error(fake1[i], find1[i]);
		
		good[i] = good[i] * 100.0 / find[i];
		fake[i] = fake[i] * 100.0 / find[i];
		good1[i] = good1[i] * 100.0 / find1[i];
		fake1[i] = fake1[i] * 100.0 / find1[i];
	}
	
	if (save) {
		fstream data(save, ios::app);
		data.setf(ios::fixed);
		data.precision(1);
		data << good1[0] << " " << fake1[0] << " ";
		data << good1[1] << " " << fake1[1] << endl;
		data.close();
	}
	else {
		cout.setf(ios::fixed);
		cout.precision(1);
		cout << "*****************************************" << endl;
		cout << "* Tracks with at least 5 correct points *" << endl;
		cout << "*****************************************" << endl;
		cout << "(all particles)" << endl;
		cout << "Efficiency: " << good[0] << " +/- " << sgood[0] << "%" << endl;
		cout << "Fake prob.: " << fake[0] << " +/- " << sfake[0] << "%" << endl;
		if (!low) {
			cout << "(pt >= 1 GeV/c)" << endl;
			cout << "Efficiency: " << good1[0] << " +/- " << sgood1[0] << "%" << endl;
			cout << "Fake prob.: " << fake1[0] << " +/- " << sfake1[0] << "%" << endl;
		}
		cout << endl;
		cout << "************************************" << endl;
		cout << "* Tracks with all 6 correct points *" << endl;
		cout << "************************************" << endl;
		cout << "(all particles)" << endl;
		cout << "Efficiency: " << good[1] << " +/- " << sgood[1] << "%" << endl;
		cout << "Fake prob.: " << fake[1] << " +/- " << sfake[1] << "%" << endl;
		if (!low) {
			cout << "(pt >= 1 GeV/c)" << endl;
			cout << "Efficiency: " << good1[1] << " +/- " << sgood1[1] << "%" << endl;
			cout << "Fake prob.: " << fake1[1] << " +/- " << sfake1[1] << "%" << endl;
		}
		cout << endl;
	}
}

Double_t error(Double_t x, Double_t y) {
	Double_t radq = x + x * x / y;
	return 100.0 * sqrt(radq) / y;
}
