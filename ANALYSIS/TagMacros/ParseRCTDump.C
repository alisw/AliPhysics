// Macro to parse the text dump of the RCT and convert it into
// a RunTag object. This object needs to be then merged with the
// proper RunTag from the Run.
//
// Parameters:
// Name of the text file with the dump from MonaLisa RCT page
//
// Author: Adam.Kisiel@cern.ch

int getposition(const char *token, int poscount)
{
  int position = -1;

  if (strstr(token, "#raw_run")) {
    position = 0;
    cout << "Run number at " << poscount << endl;
  }
  else if (strstr(token, "partition")) {
    position = 1;
    cout << "LHC Period at " << poscount << endl;
  }
  else if (strstr(token, "fillno")) {
    position = 4;
    cout << "Fill number at " << poscount << endl;
  }
  else if (strstr(token, "energy")) {
    position = 5;
    cout << "Beam Energy at " << poscount << endl;
  }
  else if (strstr(token, "intensity_per_bunch")) {
    position = 6;
    cout << "Bunch Intensity at " << poscount << endl;
  }
  else if (strstr(token, "cbeamb_abce_nopf_all")) {
    position = 7;
    cout << "Beam Triggers at " << poscount << endl;
  }
  else if (strstr(token, "cint1b_abce_nopf_all")) {
    position = 8;
    cout << "Collision Triggers at " << poscount << endl;
  }
  else if (strstr(token, "rate")) {
    position = 9;
    cout << "Collision rate at " << poscount << endl;
  }
  else if (strstr(token, "cint1_e_nopf_all")) {
    position = 10;
    cout << "Empty Triggers at " << poscount << endl;
  }
  else if (strstr(token, "cint1a_abce_nopf_all")) {
    position = 11;
    cout << "ASide Triggers at " << poscount << endl;
  }
  else if (strstr(token, "cint1c_abce_nopf_all")) {
    position = 12;
    cout << "CSide Triggers at " << poscount << endl;
  }
  else if (strstr(token, "mean_vertex_xyz")) {
    position = 13;
    cout << "Mean Vertex at " << poscount << endl;
  }
  else if (strstr(token, "vertex_quality")) {
    position = 14;
    cout << "Vertex Quality at " << poscount << endl;
  }
  else if (strstr(token, "field")) {
    position = 15;
    cout << "Magnetic Field at " << poscount << endl;
  }
  else if (strstr(token, "det_aco")) {
    position = 16;
    cout << "ACORDE at " << poscount << endl;
  }
  else if (strstr(token, "det_emc")) {
    position = 17;
    cout << "EMC at " << poscount << endl;
  }
  else if (strstr(token, "det_fmd")) {
    position = 18;
    cout << "FMD at " << poscount << endl;
  }
  else if (strstr(token, "det_hlt")) {
    position = 19;
    cout << "HLT at " << poscount << endl;
  }
  else if (strstr(token, "det_hmp")) {
    position = 20;
    cout << "HMPID at " << poscount << endl;
  }
  else if (strstr(token, "det_mch")) {
    position = 21;
    cout << "MCH at " << poscount << endl;
  }
  else if (strstr(token, "det_mtr")) {
    position = 22;
    cout << "MTR at " << poscount << endl;
  }
  else if (strstr(token, "det_phs")) {
    position = 23;
    cout << "PHOS at " << poscount << endl;
  }
  else if (strstr(token, "det_pmd")) {
    position = 24;
    cout << "pmd at " << poscount << endl;
  }
  else if (strstr(token, "det_spd")) {
    position = 25;
    cout << "SPD at " << poscount << endl;
  }
  else if (strstr(token, "det_sdd")) {
    position = 26;
    cout << "SDD at " << poscount << endl;
  }
  else if (strstr(token, "det_ssd")) {
    position = 27;
    cout << "SSD at " << poscount << endl;
  }
  else if (strstr(token, "det_tof")) {
    position = 28;
    cout << "TOF at " << poscount << endl;
  }
  else if (strstr(token, "det_tpc")) {
    position = 29;
    cout << "TPC at " << poscount << endl;
  }
  else if (strstr(token, "det_trd")) {
    position = 30;
    cout << "TRD at " << poscount << endl;
  }
  else if (strstr(token, "det_t00")) {
    position = 31;
    cout << "TZERO at " << poscount << endl;
  }
  else if (strstr(token, "det_v00")) {
    position = 32;
    cout << "VZERO at " << poscount << endl;
  }
  else if (strstr(token, "det_zdc")) {
    position = 33;
    cout << "ZDC at " << poscount << endl;
  }
  else if (strstr(token, "filling_config")) {
    position = 39;
    cout << "Filling Scheme at " << poscount << endl;
  }
  else if (strstr(token, "quality")) {
    position = 40;
    cout << "Quality at " << poscount << endl;
  }
  else if (strstr(token, "muon_trigger")) {
    position = 41;
    cout << "Muon Trigger at " << poscount << endl;
  }
  else if (strstr(token, "high_multiplicity_trigger")) {
    position = 42;
    cout << "High Multiplicity trigger at " << poscount << endl;
  }

  return position;
}

void ParseRCTDump(const char *insum = "rctdump.txt")
{
  ifstream *istr = new ifstream(insum);

  char buf[10000];
  istr->getline(buf, 10000);
  
  if (strstr(buf, "#raw_run"))
    cout << "Found raw run " << endl;

  char *nval;
  char *nlin = strstr(buf, "|");
  char *nbeg = buf-1;
  char ntok[1000];
  char nstr[1000];

  // Parse field name list

  int fieldpos[200];
  int poscount = 0;
  
  do {
    if (nlin == NULL) {
      strcpy(ntok, nbeg+1);
      //      ntok[nlin-nbeg-1] = '\0';
    }
    else {
      strncpy(ntok, nbeg+1, nlin-nbeg-1);
      ntok[nlin-nbeg-1] = '\0';
    }

    cout << "Token is " << ntok << endl;

    fieldpos[poscount] = getposition(ntok, poscount);
    
    nbeg = nlin;
    nlin = strstr(nlin+1, "|");
    poscount++;
  }
  while (nlin != NULL);

  strcpy(ntok, nbeg+1);

  cout << "Token is " << ntok << endl;
  
  fieldpos[poscount] = getposition(ntok, poscount);
  
  TFile* ftag = TFile::Open("RunTableTagSummary.root", "recreate");

  AliRunTag *tag = new AliRunTag();
  TTree * ttag = new TTree("T","A Tree with event tags");
  TBranch * btag = ttag->Branch("AliTAG", &tag);
  btag->SetCompressionLevel(9);

  int nfield = 0;


  while (!istr->eof()) {
    istr->getline(buf, 10000);
    nlin = strstr(buf, "|");
    nbeg = buf-1;

    nfield = 0;

    cout << "Line " << nlin << endl;
    
    int endline = 0;

    while (!endline) {
      
      if (nlin == NULL) {
	strcpy(ntok, nbeg+1);
	endline = 1;
      }
      else {
	strncpy(ntok, nbeg+1, nlin-nbeg-1);
	ntok[nlin-nbeg-1] = '\0';
      }
      if (ntok[0] == '"') {
	strncpy(nstr, ntok+1, strlen(ntok)-2); 
	nstr[strlen(ntok)-2] = '\0';
      }
      else if (ntok[0] == '(') {
	char *com = strstr(ntok, ";");
	cout << ntok << " "<< com << endl;
	strncpy(nstr, ntok+1, com-ntok-1);
	nstr[com-ntok-1] = '\0';
      }
      else {
	nstr[0] = '\0';
	nstr[1] = '\0';
	nstr[2] = '\0';
      }

      if (nlin != NULL) {
	nbeg = nlin;
	nlin = strstr(nlin+1, "|");
      }

      switch (fieldpos[nfield]) {
      case 0:
	tag->SetRunId(atoi(ntok));
	break;
      case 1:
	tag->SetLHCPeriod(nstr);
	break;
      case 2:
	break;
      case 3:
	break;
      case 4:
	tag->GetLHCTag()->SetFillNo(atoi(ntok));
	break;
      case 5:
	tag->GetLHCTag()->SetBeamEnergy(atof(ntok));
	break;
      case 6:
	tag->GetLHCTag()->SetBunchIntensity(atof(ntok));
	break;
      case 7:
	tag->SetBeamTriggers(atol(ntok));
	break;
      case 8:
	tag->SetCollisionTriggers(atol(ntok));
	break;
      case 9:
	tag->SetCollisionRate(atof(ntok));
	break;
      case 10:
	tag->SetEmptyTriggers(atol(ntok));
	break;
      case 11:
	tag->SetASideTriggers(atol(ntok));
	break;
      case 12:
	tag->SetCSideTriggers(atol(ntok));
	break;
      case 13:
	tag->SetMeanVertex(atof(ntok));
	break;
      case 14:
	tag->SetVertexQuality(atof(ntok));
	break;
      case 15:
	tag->SetMagneticField(atoi(ntok+1));
	break;
      case 16:
	tag->GetDetectorTags()->SetDetectorStatus(17, nstr);
	break;
      case 17:
	tag->GetDetectorTags()->SetDetectorStatus(18, nstr);
	break;
      case 18:
	tag->GetDetectorTags()->SetDetectorStatus(12, nstr);
	break;
      case 19:
	tag->GetDetectorTags()->SetDetectorStatus(21, nstr);
	break;
      case 20:
	tag->GetDetectorTags()->SetDetectorStatus(6, nstr);
	break;
      case 21:
	tag->GetDetectorTags()->SetDetectorStatus(10, nstr);
	break;
      case 22:
	tag->GetDetectorTags()->SetDetectorStatus(11, nstr);
	break;
      case 23:
	tag->GetDetectorTags()->SetDetectorStatus(7, nstr);
	break;
      case 24:
	tag->GetDetectorTags()->SetDetectorStatus(9, nstr);
	break;
      case 25:
	tag->GetDetectorTags()->SetDetectorStatus(0, nstr);
	break;
      case 26:
	tag->GetDetectorTags()->SetDetectorStatus(1, nstr);
	break;
      case 27:
	tag->GetDetectorTags()->SetDetectorStatus(2, nstr);
	break;
      case 28:
	tag->GetDetectorTags()->SetDetectorStatus(5, nstr);
	break;
      case 29:
	tag->GetDetectorTags()->SetDetectorStatus(3, nstr);
	break;
      case 30:
	tag->GetDetectorTags()->SetDetectorStatus(4, nstr);
	break;
      case 31:
	tag->GetDetectorTags()->SetDetectorStatus(13, nstr);
	break;
      case 32:
	tag->GetDetectorTags()->SetDetectorStatus(14, nstr);
	break;
      case 33:
	tag->GetDetectorTags()->SetDetectorStatus(15, nstr);
	break;
      case 34:
	// comment
	break;
      case 35:
	// detectors
	break;
      case 36:
	// changedby
	break;
      case 37:
	// changedon
	break;
      case 38:
	// run duration
	break;
      case 39:
	tag->GetLHCTag()->SetFillingScheme(nstr);
	break;
      case 40:
	cout << "Setting Quality ]" << nstr << "[" << atoi(nstr) << endl;
	if (strlen(nstr) == 0)
	  tag->SetRunQuality(0);
	else
	  tag->SetRunQuality(atoi(nstr));
	break;
      case 41:
	tag->SetMuonTriggers(atol(ntok));
	break;
      case 42:
	tag->SetHMTriggers(atol(ntok));
	break;
      }

      nfield++;
    }


    if (nfield > 1) {
      ttag->Fill();
      tag->Clear("");
    }
  }

  ftag->cd();
  tag->Clear();
  ttag->Write();
  ftag->Close();

}
