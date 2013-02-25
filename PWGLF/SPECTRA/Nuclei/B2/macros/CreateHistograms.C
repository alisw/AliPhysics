/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
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

// macro for creating histograms
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

AliLnHistoMap* CreateHistograms(const TString& species, Bool_t simulation, Double_t maxDCAxy, Double_t maxEta, Double_t maxY, Bool_t heavyIons)
{
//
// Define an create the histograms for the analysis
//
	AliLnHistoMap* hMap = new AliLnHistoMap();
	
	Int_t A = 1;
	if(species == "Deuteron")    A = 2;
	else if(species == "Triton") A = 3;
	else if(species == "He3")    A = 3;
	else if(species == "Alpha")  A = 4;
	
	// pt
	Int_t ptBins = 100;
	Double_t ptMin = 0.*A;
	Double_t ptMax = 10.*A;
	
	// eta and rapidity
	Int_t etaBins = 300;
	Double_t etaMin = -1.5;
	Double_t etaMax = 1.5;
	
	// phi
	Int_t phiBins = 180;
	Double_t phiMin = 0.;
	Double_t phiMax = 2.*TMath::Pi();
	
	// DCA distribution bin width = 0.001 cm
	Int_t dcaxyBins = 2.*maxDCAxy/0.001;
	Double_t dcaxyMin = -maxDCAxy;
	Double_t dcaxyMax = maxDCAxy;
	
	Int_t dcazBins = 1300;
	Double_t dcazMin = -3.25;
	Double_t dcazMax = 3.25;
	
	// m2 bins
	Int_t m2Bins = 400;
	Double_t m2Min = 0;
	Double_t m2Max = 20;
	
	// track multiplicity
	Int_t ntrkBins = 200;
	Double_t ntrkMin = 0;
	Double_t ntrkMax = 200;
	
	if(heavyIons)
	{
		ntrkBins = 2000;
		ntrkMin  = 0;
		ntrkMax  = 2000;
	}
	
	// stats
	
	TH1D* hStats = new TH1D(species + "_Stats", "Stats", 8, 0, 8);
	
	hStats->GetXaxis()->SetBinLabel(1,"Events");
	hStats->GetXaxis()->SetBinLabel(2,"Trig");
	hStats->GetXaxis()->SetBinLabel(3,"Ana");
	hStats->GetXaxis()->SetBinLabel(4,"Vtx");
	hStats->GetXaxis()->SetBinLabel(5,"Vz");
	hStats->GetXaxis()->SetBinLabel(6,"Vxy");
	hStats->GetXaxis()->SetBinLabel(7,"VertexerZ");
	hStats->GetXaxis()->SetBinLabel(8,"NoPileUp");
	
	hStats->SetStats(0);
	hMap->Add( species + "_Stats", (TObject*)hStats);
	
	if(simulation)
	{
		// PID table
		
		TH2D* hPidTable = new TH2D(species + "_Stats_PID_Table", "Particle table (Row: PID, Col: MC)", 10,0,10,10,0,10);
		
		hPidTable->GetXaxis()->SetBinLabel(1,"e");
		hPidTable->GetXaxis()->SetBinLabel(2,"#mu");
		hPidTable->GetXaxis()->SetBinLabel(3,"#pi");
		hPidTable->GetXaxis()->SetBinLabel(4,"K");
		hPidTable->GetXaxis()->SetBinLabel(5,"p");
		hPidTable->GetXaxis()->SetBinLabel(6,"d");
		hPidTable->GetXaxis()->SetBinLabel(7,"t");
		hPidTable->GetXaxis()->SetBinLabel(8,"h");
		hPidTable->GetXaxis()->SetBinLabel(9,"#alpha");
		hPidTable->GetXaxis()->SetBinLabel(10,"sum");
		
		hPidTable->GetYaxis()->SetBinLabel(1,"e");
		hPidTable->GetYaxis()->SetBinLabel(2,"#mu");
		hPidTable->GetYaxis()->SetBinLabel(3,"#pi");
		hPidTable->GetYaxis()->SetBinLabel(4,"K");
		hPidTable->GetYaxis()->SetBinLabel(5,"p");
		hPidTable->GetYaxis()->SetBinLabel(6,"d");
		hPidTable->GetYaxis()->SetBinLabel(7,"t");
		hPidTable->GetYaxis()->SetBinLabel(8,"h");
		hPidTable->GetYaxis()->SetBinLabel(9,"#alpha");
		hPidTable->GetYaxis()->SetBinLabel(10,"sum");
		
		hPidTable->SetStats(0);
		hMap->Add( species + "_Stats_PID_Table", (TObject*)hPidTable);
	}
	
	TH1D* hStatsPid = new TH1D(species + "_Stats_PID", "Stats PID", 9, 0, 9);
	hStatsPid->GetXaxis()->SetBinLabel(1,"e");
	hStatsPid->GetXaxis()->SetBinLabel(2,"\\mu");
	hStatsPid->GetXaxis()->SetBinLabel(3,"\\pi");
	hStatsPid->GetXaxis()->SetBinLabel(4,"K");
	hStatsPid->GetXaxis()->SetBinLabel(5,"p");
	hStatsPid->GetXaxis()->SetBinLabel(6,"d");
	hStatsPid->GetXaxis()->SetBinLabel(7,"t");
	hStatsPid->GetXaxis()->SetBinLabel(8,"h");
	hStatsPid->GetXaxis()->SetBinLabel(9,"\\alpha");
	hStatsPid->SetStats(0);
	hMap->Add( species + "_Stats_PID", (TObject*)hStatsPid);
	
	// detector signals
	
	hMap->Add( species + "_ITS_dEdx_P", 1000, 0.001, 10., 1500, 0.001, 1500., "", "p/Z (GeV/c)", "dE/dx");
	hMap->Add( species + "_TPC_dEdx_P", 1000, 0.001, 10., 1500, 0.001, 1500., "", "p/Z (GeV/c)", "dE/dx");
	hMap->Add( species + "_TOF_Beta_P", 1000, 0.001, 10., 1200, 0.001, 1.2, "", "p/Z (GeV/c)", "#beta");
	hMap->Add( species + "_TOF_Mass_P", 1000, 0.001, 10., 500, 0.01, 5., "", "p/Z (GeV/c)", "Mass (GeV/c^{2})");
	
	// TrackCuts
	
	hMap->Add( species + "_TrackCuts_DCAxy", 200, -4, 4, "Track cuts", "DCA_{xy} (cm)");
	hMap->Add( species + "_TrackCuts_DCAz", 200, -6, 6, "Track cuts", "DCA_{z} (cm)");
	hMap->Add( species + "_TrackCuts_NSigma", 200, 0, 10, "Track cuts", "N_{#sigma}");
	
	hMap->Add( species + "_TrackCuts_ITSchi2PerCls", 200, 0, 40, "Track cuts", "ITS #chi^{2} / Cluster");
	
	hMap->Add( species + "_TrackCuts_TPCncls", 200, 0, 200, "Track cuts", "TPC clusters");
	hMap->Add( species + "_TrackCuts_TPCclsOverF", 200, 0, 2, "Track cuts", "TPC clusters/Findable");
	hMap->Add( species + "_TrackCuts_TPCxRows", 200, 0, 200, "Track cuts", "TPC crossed rows");
	hMap->Add( species + "_TrackCuts_TPCchi2PerCls", 200, 0, 20, "Track cuts", "TPC #chi^{2} / Cluster");
	hMap->Add( species + "_TrackCuts_TPCchi2Global", 200, -2, 38, "Track cuts", "#chi^{2} of constrained TPC vs global track");
	
	hMap->Add( species + "_Before_Phi_Theta", 180, 0, TMath::Pi(), 360, 0, 2.*TMath::Pi(), "Global tracks (before track cuts)", "#theta (rad)", "#phi (rad)");
	
	hMap->Add( species + "_After_Phi_Theta", 180, 0, TMath::Pi(), 360, 0, 2.*TMath::Pi(), "Global tracks (after track cuts)", "#theta (rad)", "#phi (rad)");
	
	// multiplicity
	
	hMap->Add( species + "_Event_Ntrk", ntrkBins, ntrkMin, ntrkMax, "Track multiplicity (all events)", "N_{trk}", "Events");
	hMap->Add( species + "_Ana_Event_Ntrk", ntrkBins, ntrkMin, ntrkMax, "Track multiplicity", "N_{trk}", "Events");
	
	if(simulation)
	{
		hMap->Add( species + "_Ana_Event_Nch_Ntrk", 200, 0., 200., 200, 0., 200., "Track multiplicity", "N_{trk}", "N_{ch}");
	}
	
	hMap->Add( species + "_Vertex_XZ", 500, -50., 50., 800, -1., 1., "Primary Vertex", "Z (cm)", "X (cm)");
	hMap->Add( species + "_Vertex_YZ", 500, -50., 50., 800, -1., 1., "Primary Vertex", "Z (cm)", "Y (cm)");
	hMap->Add( species + "_Vertex_YX", 800, -1., 1., 800, -1., 1., "Primary Vertex", "X (cm)", "Y (cm)");
	
	hMap->Add( species + "_Ana_Vertex_XZ", 500, -50., 50., 800, -1., 1., "Primary Vertex", "Z (cm)", "X (cm)");
	hMap->Add( species + "_Ana_Vertex_YZ", 500, -50., 50., 800, -1., 1., "Primary Vertex", "Z (cm)", "Y (cm)");
	hMap->Add( species + "_Ana_Vertex_YX", 800, -1., 1., 800, -1., 1., "Primary Vertex", "X (cm)", "Y (cm)");
	
	if(heavyIons)
	{
		hMap->Add( species + "_V0_Mult", 1000, 0, 100000., "Corrected V0 Multiplicity", "V0M", "Events");
		hMap->Add( species + "_V0_Scaled_Mult", 1000, 0, 10000., "Corrected Scaled V0 Multiplicity", "V0SM", "Events");
	
		hMap->Add( species + "_Ana_V0_Mult", 1000, 0, 100000., "Corrected V0 Multiplicity (after centrality sel.)", "V0M", "Events");
		hMap->Add( species + "_Ana_V0_Scaled_Mult", 1000, 0, 10000., "Corrected Scaled V0 Multiplicity (after centrality sel.)", "V0SM", "Events");
	}
	
	// positive and negative
	
	TString particle[] = { Form("Anti%s",species.Data()), species };
	
	for(Int_t i=0; i<2; ++i)
	{
		// pid with TOF
		
		hMap->Add( particle[i] + "_PID_M2_Pt", ptBins, ptMin, ptMax, m2Bins, m2Min, m2Max, Form("%s candidate",particle[i].Data()), "p_{T} (GeV/c)", "m^{2} (GeV^{2}/c^{4})");
		
		// as a function of momentum
		hMap->Add( particle[i] + "_PID_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("%s candidate",particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		
		hMap->Add( particle[i] + "_PID_DCAz_Pt", ptBins, ptMin, ptMax, dcazBins, dcazMin, dcazMax, Form("%s candidate",particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{z} (cm)");
		
		hMap->Add( particle[i] + "_PID_NSigma_Pt", ptBins, ptMin, ptMax, 200, 0., 10., Form("%s candidate",particle[i].Data()), "p_{T} (GeV/c)", "N\\sigma");
		
		if(simulation)
		{
			hMap->Add( particle[i] + "_Sim_PID_M2_Pt", ptBins, ptMin, ptMax, m2Bins, m2Min, m2Max, Form("%s after PID", particle[i].Data()), "p_{T} (GeV/c)", "m^{2} (GeV^{2}/c^{4})");
			
			hMap->Add( particle[i] + "_Sim_Prim_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Primary %s", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		
			hMap->Add( particle[i] + "_Sim_Fdwn_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Feed-down %s", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
			
			hMap->Add( particle[i] + "_Sim_Mat_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("%s from materials", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
			
			// primaries
			hMap->Add( particle[i] + "_Sim_PID_Prim_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Primary %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Prim_DCAz_Pt", ptBins, ptMin, ptMax, dcazBins, dcazMin, dcazMax, Form("Primary %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{z} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Prim_NSigma_Pt", ptBins, ptMin, ptMax, 200, 0., 10., Form("Primary %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "N\\sigma");
			
			// feed-down
			hMap->Add( particle[i] + "_Sim_PID_Fdwn_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Feed-down %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Fdwn_DCAz_Pt", ptBins, ptMin, ptMax, dcazBins, dcazMin, dcazMax, Form("Feed-down %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{z} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Fdwn_NSigma_Pt", ptBins, ptMin, ptMax, 200, 0., 10., Form("Feed-down %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "N\\sigma");
			
			// secondaries from materials
			hMap->Add( particle[i] + "_Sim_PID_Mat_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("%s from materials after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Mat_DCAz_Pt", ptBins, ptMin, ptMax, dcazBins, dcazMin, dcazMax, Form("%s from materials after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{z} (cm)");
		
			hMap->Add( particle[i] + "_Sim_PID_Mat_NSigma_Pt", ptBins, ptMin, ptMax, 200, 0., 10., Form("%s from materials after PID", particle[i].Data()), "p_{T} (GeV/c)", "N\\sigma");
			
			// fake tracks
			hMap->Add( particle[i] + "_Sim_PID_Fake_Prim_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Fake Primary %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
			
			hMap->Add( particle[i] + "_Sim_PID_Fake_Fdwn_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Fake feed-down %s after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
			
			hMap->Add( particle[i] + "_Sim_PID_Fake_Mat_DCAxy_Pt", ptBins, ptMin, ptMax, dcaxyBins, dcaxyMin, dcaxyMax, Form("Fake %s from materials after PID", particle[i].Data()), "p_{T} (GeV/c)", "sign #times DCA_{xy} (cm)");
		}
		
		// identification
		
		hMap->Add( particle[i] + "_PID_ITSdEdx_P", 1000, 1.e-3, 10., 1500, 1.e-3, 1500., Form("%s candidate",particle[i].Data()), "p (GeV/c)", "dE/dx");
		
		hMap->Add( particle[i] + "_PID_TPCdEdx_P", 1000, 1.e-3, 10., 1500, 1.e-3, 1500., Form("%s candidate",particle[i].Data()), "p_{TPC} (GeV/c)", "dE/dx");
		
		hMap->Add( particle[i] + "_PID_Beta_P", 1000, 1.e-3, 10., 1200, 1.e-3, 1.2, Form("%s candidate",particle[i].Data()), "p_{TOF} (GeV/c)", "\\beta");
		
		hMap->Add( particle[i] + "_PID_Mass_P", 1000, 1.e-3, 10., 500, 0.01, 5., Form("%s candidate",particle[i].Data()), "p_{TOF} (GeV/c)", "m (GeV/c^{2})");
		
		if(simulation)
		{
			hMap->Add( particle[i] + "_Sim_PID_Mass", 500, 0.01, 5., Form("%s after PID", particle[i].Data()), "m (GeV/c^{2})");
		}
		
		// results
		
		hMap->Add( particle[i] + "_PID_Pt_Y", etaBins, etaMin, etaMax, ptBins, ptMin, ptMax, Form("%s candidate",particle[i].Data()), "y", "p_{T} (GeV/c)");
		
		hMap->Add( particle[i] + "_PID_Pt", ptBins, ptMin, ptMax, Form("%s candidate (|y| < %0.1f, 0 < \\phi < 2\\pi)", particle[i].Data(),maxY), "p_{T} (GeV/c)");
		
		hMap->Add( Form("%s_PID_Y",particle[i].Data()), etaBins, etaMin, etaMax, Form("%s candidate (p_{T} > 0, 0 < \\phi < 2\\pi)",particle[i].Data()), "y");
		
		hMap->Add( particle[i] + "_PID_Phi", phiBins, phiMin, phiMax, Form("%s candidate (p_{T} > 0, |y| < %0.1f)",particle[i].Data(),maxY), "\\phi (rad)");
		
		if(simulation)
		{
			hMap->Add( particle[i] + "_Sim_Pt", ptBins, ptMin, ptMax, Form("%s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Prim_Phi", phiBins, phiMin, phiMax, Form("Primary %s (p_{T} > 0, |y| < %0.1f)",particle[i].Data(),maxY), "\\phi (rad)");
			
			hMap->Add( particle[i] + "_Sim_Prim_Y", etaBins, etaMin, etaMax, Form("Primary %s (p_{T} > 0, 0 < \\phi < 2\\pi)",particle[i].Data()), "y");
			
			hMap->Add( particle[i] + "_Sim_Fdwn_Pt", ptBins, ptMin, ptMax, Form("Feed-down %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Mat_Pt", ptBins, ptMin, ptMax, Form("%s from materials (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_PID_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s after PID (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
		
			hMap->Add( particle[i] + "_Sim_PID_Pt", ptBins, ptMin, ptMax, Form("%s after PID (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
		
			hMap->Add( particle[i] + "_Sim_PID_Prim_R", 200, 0., 200., Form("%s after PID (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "R_{xyz} (cm)");
			
			hMap->Add( particle[i] + "_Sim_PtY", etaBins, etaMin, etaMax, ptBins, ptMin, ptMax, particle[i].Data(), "y", "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Prim_M2_P", ptBins, ptMin, ptMax, m2Bins, m2Min, m2Max, Form("Primary %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p (GeV/c)", "m^{2} (GeV^{2}/c^{4})");
			
			hMap->Add( particle[i] + "_Sim_Prim_M2_Pt", ptBins, ptMin, ptMax, m2Bins, m2Min, m2Max, Form("Primary %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)", "m^{2} (GeV^{2}/c^{4})");
			
			// unfolding
			
			hMap->Add( particle[i] + "_Response_Matrix",  ptBins, ptMin, ptMax,  ptBins, ptMin, ptMax, particle[i].Data(), "Measured p_{T} (GeV/c)", "True p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Prim_Response_Matrix",  ptBins, ptMin, ptMax,  ptBins, ptMin, ptMax, Form("Primary %s",particle[i].Data()), "Measured p_{T} (GeV/c)", "True p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Prim_DiffPt_RecPt",  ptBins, ptMin, ptMax,  1000, -0.5, 0.5, Form("Primary %s",particle[i].Data()), "p_{T}^{rec} (GeV/c)", "p_{T}^{gen}-p_{T}^{rec} (GeV/c)");
			
			// fake tracks
			hMap->Add( particle[i] + "_Sim_Fake_Pt", ptBins, ptMin, ptMax, Form("Fake %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Fake_Prim_Pt", ptBins, ptMin, ptMax, Form("Fake Primary %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Fake_Fdwn_Pt", ptBins, ptMin, ptMax, Form("Fake Feed-down %s (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_Fake_Mat_Pt", ptBins, ptMin, ptMax, Form("Fake %s from materials (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Sim_PID_Fake_Pt", ptBins, ptMin, ptMax, Form("Fake %s after PID (|y| < %0.1f, 0 < \\phi < 2\\pi)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
		}
	}
	
	// generation
	
	if(simulation)
	{
		for(Int_t i=0; i<2; ++i)
		{
			hMap->Add( particle[i] + "_Gen_Prim_P", ptBins, ptMin, ptMax, Form("Primary %s", particle[i].Data()), "p (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s", particle[i].Data()), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_Prim_Y", 600, -12, 12, Form("Primary %s", particle[i].Data()), "y");
			
			hMap->Add( particle[i] + "_Gen_Prim_Eta", 600, -12, 12, Form("Primary %s", particle[i].Data()), "#eta");
			
			hMap->Add( particle[i] + "_Gen_Prim_Phi", phiBins, phiMin, phiMax, Form("Primary %s", particle[i].Data()), "\\phi (rad)");
			
			hMap->Add( particle[i] + "_Gen_Prim_PtY", etaBins, etaMin, etaMax, ptBins, ptMin, ptMax, Form("Primary %s", particle[i].Data()), "y", "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_Prim_EtaY", etaBins, etaMin, etaMax, etaBins, etaMin, etaMax, Form("Primary %s", particle[i].Data()), "y", "#eta");
			
			// after multiplicity cut
			
			hMap->Add( particle[i] + "_Gen_Mult_Prim_PtY", etaBins, etaMin, etaMax, ptBins, ptMin, ptMax, Form("Primary %s", particle[i].Data()), "y", "p_{T} (GeV/c)");
			
			// in the phase space within the acceptance
			
			hMap->Add( particle[i] + "_Gen_PhS_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (all events, |y| < %0.1f)",particle[i].Data(),maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_Trig_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (triggering events, |y| < %0.1f)", particle[i].Data(), maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_Vtx_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (good vertex, |y| < %0.1f)", particle[i].Data(), maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_NoTrig_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (non-triggering events, |y| < %0.1f)", particle[i].Data(), maxY), "p_{T} (GeV/c)");
			
			hMap->Add( particle[i] + "_Gen_NoVtx_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (no reconstructed vertex, |y| < %0.1f)", particle[i].Data(), maxY), "p_{T} (GeV/c)");
			
			// within the acceptance
			
			hMap->Add( particle[i] + "_Gen_Acc_Prim_P", ptBins, ptMin, ptMax, Form("Primary %s (|y| < %0.1f and |\\eta| < %0.1f)",particle[i].Data(),maxY,maxEta), "p (GeV/c)");
		
			hMap->Add( particle[i] + "_Gen_Acc_Prim_Pt", ptBins, ptMin, ptMax, Form("Primary %s (|y| < %0.1f and |\\eta| < %0.1f)",particle[i].Data(),maxY,maxEta), "p_{T} (GeV/c)");
		
			hMap->Add( particle[i] + "_Gen_Acc_Prim_Y", etaBins, etaMin, etaMax, Form("Primary %s (|\\eta| < %0.1f)",particle[i].Data(),maxEta), "y");
		
			hMap->Add( particle[i] + "_Gen_Acc_Prim_Phi", phiBins, phiMin, phiMax, Form("Primary %s (|y|< %0.1f and |\\eta| < %0.1f)",particle[i].Data(),maxY,maxEta), "\\phi (rad)");
			
			hMap->Add( particle[i] + "_Gen_Acc_Prim_PtY", etaBins, etaMin, etaMax, ptBins, ptMin, ptMax, Form("Primary %s",particle[i].Data()), "y", "p_{T} (GeV/c)");
		}
	}
	
	return hMap;
}
