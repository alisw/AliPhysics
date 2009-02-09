  void UpdateCDBGRPEntryMC() {
	// produce the GRP entry in CDB 
	// reading MC parameter from Config.C
	
	AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(0);
	AliCDBId id(AliQA::GetGRPPath(), 0, AliCDBRunRange::Infinity());
  AliCDBMetaData *metadata= new AliCDBMetaData();

  // Get root version
  const char* rootv = gROOT->GetVersion();
	metadata->SetResponsible("Yves.Schutz@cern.ch");
	metadata->SetComment("GRP parameters for Monte Carlo");
 
		
	TMap *mappp = GetGRPList();
		
  //  Printf("Storing in CDB the default values for the GRP %d parameters produced with root %s and AliRoot version %s",list->GetEntries(),rootv,alirootv);
	
  man->Put(mappp,id,metadata);
}

TString ParseConfig(char * option)
{
	// Parse Config file to retrieve usefull information for the GRP
	
	ifstream in("Config.C", ios::in) ; 
  if (!in)
		AliError("Config.C file not found in current path") ; 	

	TString rv("") ; 
	
	TString soption(option) ;
	
	string line ; 
	TString sline ; 

	// search for the system p-p or Pb-Pb or ..... 
	if (soption.Contains("system") || soption.Contains("fAliceBeamType")){
		while (getline(in, line)) { 
			sline = line ;		
			if (sline.Contains("gAlice->SetTriggerDescriptor")) 
				break ;
		}
		TString sarg1 = sline(sline.Index("[")+1, sline.Index("]") - sline.Index("[")-1) ; 
		sarg1 = sarg1.Strip(TString::kLeading) ; 
		sarg1 = sarg1.Strip(TString::kTrailing) ; 
		TString sarg2 = sline(sline.Index("(")+1, sline.Index("[") - sline.Index("(")-1) ; 
		sarg2 = sarg2.Strip(TString::kLeading) ; 
		sarg2 = sarg2.Strip(TString::kTrailing) ; 
		in.seekg(0) ;
		while (getline(in, line)) {
			sline = line ;
			if (sline.Contains(sarg1)) {
				sarg1 = sline(sline.Index("=")+1, sline.Index(";")-sline.Index("=")-1) ; 
				sarg1 = sarg1.Strip(TString::kLeading) ; 
				sarg1 = sarg1.Strip(TString::kTrailing) ; 
				break ; 
			}
		}
		in.seekg(0) ;
		Int_t index = -1 ; 
		while (getline(in, line)) {
			sline = line ;
			if (sline.Contains(sarg1)) {
				sline = sline.Strip(TString::kLeading) ; 				
				Ssiz_t deb = sline.Index(sarg1) ;
				Ssiz_t kom = sline.Index(",") ; 
				if ( deb == 0 || deb < kom)
					index = 1 ; 
				else 
					index = 2 ; 
				continue ; 
			}
			if (sline.Contains(sarg2)) {
				Ssiz_t deb = sline.First("\"") ; 
				if ( deb == -1 ) {
					getline(in, line) ; 
					sline = line ; 
				}
				sline.ReplaceAll("\"", "") ;
				sline = sline.Strip(TString::kLeading) ; 
				if (index == 1 ) 
					rv = sline(0, sline.Index(",")) ;
				else if (index == 2) 
					rv = sline(sline.Index(",")+1, sline.Length()-sline.Index(",")) ; 
				break ;
			}
		}
	} 
	// get the beam energy
	else if (soption.Contains("fAliceBeamEnergy")) {
		while (getline(in, line)) { 
			sline = line ;		
			if ( sline.Contains("->SetEnergyCMS(") ) {
				rv = sline(sline.Index("(")+1, sline.Index(")")-sline.Index("(")-1) ; 
				break ;		
			}
		}
	}
	// get the number of active detectors
	else if (soption.Contains("fNumberOfDetectors")) {
		UShort_t det = 0 ; 
		while (getline(in, line)) {
			sline = line ; 
			sline.ReplaceAll(" ", "") ; 
			if (sline.Contains("iITS=1") ||
					sline.Contains("iTPC=1") ||
					sline.Contains("iTRD=1") ||
					sline.Contains("iTOF=1") ||			
					sline.Contains("iPHOS=1") ||
					sline.Contains("iHMPID=1") ||
					sline.Contains("iEMCAL=1") ||
					sline.Contains("iMUON=1") ||
					sline.Contains("iFMD=1") ||
					sline.Contains("iPMD=1") ||
					sline.Contains("iT0=1") ||
					sline.Contains("iVZERO=1") ||
					sline.Contains("iZDC=1") ||
					sline.Contains("iACORDE=1"))
				det++ ; 
		}
		rv = Form("%d", det) ; 
	}
	else if (soption.Contains("fL3Current")) {
		TString sarg ;
		while (getline(in, line)) {
			sline = line ; 
			sline.ReplaceAll(" ", "") ; 
			if (sline.Contains("AliMagF*field=newAliMagF(")) {
				sarg = sline(sline.Last(',')+1, sline.Last(')')-sline.Last(',')-1) ; 
				break ; 
			}
		}
		in.seekg(0) ;
		while (getline(in, line)) {
			sline = line ; 
			sline.ReplaceAll(" ", "") ; 
			if (sline.Contains(Form("%s=",sarg.Data()))) {
				sarg = sline(sline.Index("=")+1, sline.Index(";")-sline.Index("=")-1) ;
				if (sarg == "k5kG") 
					rv = "30000" ;
				else if (sarg == "k2kG")
					rv = "12000" ;
				break ; 
			}
		}		
	}
		in.close() ;
		rv.Strip(TString::kLeading) ; 
		rv.Strip() ; 
	return rv ; 
}
	
//_______________________________________//
TMap *GetGRPList() {
 
  TString system = ParseConfig("fAliceBeamType") ;  
	TString fSystem = system;
  TMap *map = new TMap();
  map->SetName("MONTECARLO");

  //DAQ
	map->Add(new TObjString("fRunType"),new TObjString(AliQA::GetRunTypeName(AliQA::kPHYSICS)));
	map->Add(new TObjString("fAliceStartTime"),new TObjString("0"));
  map->Add(new TObjString("fAliceStopTime"),new TObjString("9999"));
	map->Add(new TObjString("fAliceBeamEnergy"),new TObjString(ParseConfig("fAliceBeamEnergy")));
  map->Add(new TObjString("fAliceBeamType"),new TObjString(system));
  map->Add(new TObjString("fNumberOfDetectors"),new TObjString(ParseConfig("fNumberOfDetectors")));
  map->Add(new TObjString("fDetectorMask"),new TObjString("1074790399"));
  map->Add(new TObjString("fLHCPeriod"),new TObjString("LHC08c"));

  //DCS
  map->Add(new TObjString("fLHCState"),new TObjString("0"));
  map->Add(new TObjString("fLHCCondition"),new TObjString("0"));
  map->Add(new TObjString("fLHCLuminosity"),new TObjString("0"));
  map->Add(new TObjString("fBeamIntensity"),new TObjString("0"));
  map->Add(new TObjString("fL3Current"),new TObjString(ParseConfig("fL3Current")));
  map->Add(new TObjString("fL3Polarity"),new TObjString("0"));
  map->Add(new TObjString("fDipoleCurrent"),new TObjString("6000"));
  map->Add(new TObjString("fDipolePolarity"),new TObjString("0"));
  map->Add(new TObjString("fCavernTemperature"),new TObjString("0"));
  map->Add(new TObjString("fCavernPressure"),new TObjString("0"));

  return map;
}
