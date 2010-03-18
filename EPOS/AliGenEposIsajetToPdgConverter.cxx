//
// AliGenEposIsajetToPdgConverter.cxx
//
//  Helper class used by TEpos to insert EPOS internal objects and higher
//  resonances to PDG database.
//  Quark clusters has has unaltered code, unknown objects have code of
//  form 6xxxxxxxx, and higher resonances have codes from Particle Physics
//  Review '96
//
//  Created on: Aug 03, 2009
//      Author: Piotr Ostrowski, postrow@if.pw.edu.pl
//


#include "TDatabasePDG.h"
#include "AliGenEposIsajetToPdgConverter.h"

ClassImp(AliGenEposIsajetToPdgConverter)

AliGenEposIsajetToPdgConverter::AliGenEposIsajetToPdgConverter() {

}

AliGenEposIsajetToPdgConverter::~AliGenEposIsajetToPdgConverter() { }

Int_t AliGenEposIsajetToPdgConverter::ConvertIsajetToPdg(Int_t isajet) const{
  // Perform code conversion, init PDG if first time called
	if (!fgParticlesAdded) {
		AddHigherResonances();
	}	
	TDatabasePDG *pdgDb = TDatabasePDG::Instance();
	Int_t pdg = pdgDb->ConvertIsajetToPdg(isajet);
	if(pdg == 0) {
		pdg = ExtendedMapping(isajet);
	}
	if(pdg != 0){
		return pdg;
	}
	if(isajet / 100000000 == 8 || isajet / 100000000 == 7) {
		AddQuarkCluster(isajet);
		return isajet;
	}
	if(isajet % 100 == 99) {
		pdg = AddUnknownObject(isajet);
	}
	
	if (pdg == 0) {
		Printf("TEpos: Warning, unknown particle, ISAJET: %d\n",isajet);
	}
	return pdg;
}

void AliGenEposIsajetToPdgConverter::AddQuarkCluster(Int_t clusterCode) const {
  // Add EPOS internal object to the PDG db
	TDatabasePDG *pdgDb = TDatabasePDG::Instance();
	if(!pdgDb->GetParticle(clusterCode)) {
		pdgDb->AddParticle("Quark cluster", "Quark cluster", 0.0, kFALSE, 0, 0, "EPOS Quark cluster", clusterCode);
	}
}

Int_t AliGenEposIsajetToPdgConverter::AddUnknownObject(Int_t code) const {
  // Add EPOS internal object to the PDG db
	Int_t newCode = 600000000 + code;
	TDatabasePDG *pdgDb = TDatabasePDG::Instance();
	if(!pdgDb->GetParticle(newCode)) {
		pdgDb->AddParticle("Unknown EPOS particle", "Unknown EPOS particle", 0.0, kFALSE, 0, 0, "Unknown", newCode);
	}
	return newCode;
}

Int_t AliGenEposIsajetToPdgConverter::ExtendedMapping(Int_t isajet) const {
  // Mappings for some internal EPOS objects and resonaces
	Int_t sign = isajet < 0 ? -1 : 1;
	isajet*=sign;
	Int_t retVal = 0;
	switch(isajet) {
		//diquarks
		case 1100: retVal=2203; break;
		case 1200: retVal=2101; break;
		case 1300: retVal=3201; break;
		case 2200: retVal=1103; break;
		case 2300: retVal=3101; break;
		case 3300: retVal=3303; break;

		case 1112: retVal=2222; break;
		case 1113: retVal=12224; break;
		case 1114: retVal=12226; break; //delta++
		case 2222: retVal=1112; break;
		case 2223: retVal=11114; break;
		case 2224: retVal=11116; break; //delta-
		case 1122: retVal=12212; break;
		case 1123: retVal=22212; break;
		case 1125: retVal=32212; break;
		case 1127: retVal=42212; break; //p*
		case 1222: retVal=12112; break;
		case 1223: retVal=22112; break;
		case 1225: retVal=32112; break;
		case 1227: retVal=42112; break; //n*
		case 1124: retVal=2122; break;
		case 1126: retVal=12214; break;
		case 1128: retVal=12126; break; //delta+
		case 1224: retVal=1212; break;
		case 1226: retVal=12114; break;
		case 1228: retVal=11216; break; //delta0
		case 1132: retVal=13222; break;
		case 1133: retVal=3226; break;
		case 1134: retVal=23224; break; //sigma+
		case 1236: retVal=13212; break;
		case 1237: retVal=3216; break;
		case 1239: retVal=23214; break; //sigma0
		case 2232: retVal=13112; break;
		case 2233: retVal=3116; break;
		case 2234: retVal=23114; break; //sigma-
		case 1233: retVal=13122; break;
		case 1234: retVal=3124; break;
		case 1235: retVal=33122; break;
		case 1238: retVal=13126; break; //lambda0
	}
	return sign*retVal;
}

Bool_t AliGenEposIsajetToPdgConverter::fgParticlesAdded = kFALSE;

void AliGenEposIsajetToPdgConverter::AddHigherResonances() {
  // Adds higher resonances known by EPOS but not to the PDG
	if(fgParticlesAdded)
		return;
	fgParticlesAdded=kTRUE;
	TDatabasePDG *fPdgDb = TDatabasePDG::Instance();
	//| 332 | F0(975)     |
	fPdgDb->AddAntiParticle("F0_bar", -10221);
	//| 112 | A0(980)     |
	//| 122 | A+(980)     |
	//| 1112 | DL++(1620) |
	fPdgDb->AddParticle("delta++", "delta++", 1.630, kFALSE, 0.145, 6, "Baryon", 2222);
	fPdgDb->AddAntiParticle("delta++_bar", -2222);
	//| 1113 | DL++(1700) |
	fPdgDb->AddParticle("delta++", "delta++", 1.700, kFALSE, 0.300, 6, "Baryon", 12224);
	fPdgDb->AddAntiParticle("delta++_bar", -12224);
	//| 1114 | DL++(1925) |
	fPdgDb->AddParticle("delta++", "delta++", 1.960, kFALSE, 0.360, 6, "Baryon", 12226);
	fPdgDb->AddAntiParticle("delta++_bar", -12226);
	//| 2222 | DL-(1620) |
	fPdgDb->AddParticle("delta-", "delta-", 1.630, kFALSE, 0.145, -3, "Baryon", 1112);
	fPdgDb->AddAntiParticle("delta-_bar", -1112);
	//| 2223 | DL-(1700) |
	fPdgDb->AddParticle("delta-", "delta-", 1.700, kFALSE, 0.300, -3, "Baryon", 11114);
	fPdgDb->AddAntiParticle("delta-_bar", -11114);
	//| 2224 | DL-(1925) |
	fPdgDb->AddParticle("delta-", "delta-", 1.960, kFALSE, 0.360, -3, "Baryon", 11116);
	fPdgDb->AddAntiParticle("delta-_bar", -11116);
	//| 1122 | N*+(1440) |
	fPdgDb->AddParticle("p*", "proton (excited)", 1.440, kFALSE, 0.3, 3, "Baryon", 12212);
	fPdgDb->AddAntiParticle("p*_bar", -12212);
	//| 1123 | N*+(1530) |
	fPdgDb->AddParticle("p*", "proton (excited)", 1.535, kFALSE, 0.150, 3, "Baryon", 22212);
	fPdgDb->AddAntiParticle("p*_bar", -22212);
	//| 1125 | N*+(1665) |
	fPdgDb->AddParticle("p*", "proton (excited)", 1.655, kFALSE, 0.165, 3, "Baryon", 32212);
	fPdgDb->AddAntiParticle("p*_bar", -32212);
	//| 1127 | N*+(1710) |
	fPdgDb->AddParticle("p*", "proton (excited)", 1.710, kFALSE, 0.1, 3, "Baryon", 42212);
	fPdgDb->AddAntiParticle("p*_bar", -42212);
	//| 1222 | N*0(1440) |
	fPdgDb->AddParticle("n*", "neutron (excited)", 1.440, kFALSE, 0.3, 3, "Baryon", 12112);
	fPdgDb->AddAntiParticle("n*_bar", -12112);
	//| 1223 | N*0(1530) |
	fPdgDb->AddParticle("n*", "neutron (excited)", 1.535, kFALSE, 0.150, 3, "Baryon", 22112);
	fPdgDb->AddAntiParticle("n*_bar", -22112);
	//| 1225 | N*0(1665) |
	fPdgDb->AddParticle("n*", "neutron (excited)", 1.655, kFALSE, 0.165, 3, "Baryon", 32112);
	fPdgDb->AddAntiParticle("n*_bar", -32112);
	//| 1227 | N*0(1710) |
	fPdgDb->AddParticle("n*", "neutron (excited)", 1.710, kFALSE, 0.1, 3, "Baryon", 42112);
	fPdgDb->AddAntiParticle("n*_bar", -42112);
	//| 1124 | DL+(1620) |
	fPdgDb->AddParticle("delta+", "delta+", 1.630, kFALSE, 0.145, 3, "Baryon", 2122);
	fPdgDb->AddAntiParticle("delta+_bar", -2122);
	//| 1126 | DL+(1700) |
	fPdgDb->AddParticle("delta+", "delta+", 1.700, kFALSE, 0.300, 3, "Baryon", 12214);
	fPdgDb->AddAntiParticle("delta+_bar", -12214);
	//| 1128 | DL+(1925) |
	fPdgDb->AddParticle("delta+", "delta+", 1.960, kFALSE, 0.360, 3, "Baryon", 12126);
	fPdgDb->AddAntiParticle("delta+_bar", -12126);
	//| 1224 | DL0(1620) |
	fPdgDb->AddParticle("delta0", "delta0", 1.630, kFALSE, 0.145, 0, "Baryon", 1212);
	fPdgDb->AddAntiParticle("delta0_bar", -1212);
	//| 1226 | DL0(1700) |
	fPdgDb->AddParticle("delta0", "delta0", 1.700, kFALSE, 0.300, 0, "Baryon", 12114);
	fPdgDb->AddAntiParticle("delta0_bar", -12114);
	//| 1228 | DL0(1925) |
	fPdgDb->AddParticle("delta0", "delta0", 1.960, kFALSE, 0.360, 0, "Baryon", 11216);
	fPdgDb->AddAntiParticle("delta0_bar", -11216);
	//| 1233 | L(1405)    |
	fPdgDb->AddParticle("lambda0", "lambda0", 1.406, kFALSE, 0.050, 0, "Baryon", 13122);
	fPdgDb->AddAntiParticle("lambda0_bar", -13122);
	//| 1234 | L(1520)    |
	//| 1235 | L(1645)    |
	fPdgDb->AddParticle("lambda0", "lambda0", 1.670, kFALSE, 0.035, 0, "Baryon", 33122);
	fPdgDb->AddAntiParticle("lambda0_bar", -33122);
	//| 1238 | L(1845)    |
	fPdgDb->AddParticle("lambda0", "lambda0", 1.830, kFALSE, 0.095, 0, "Baryon", 13126);
	fPdgDb->AddAntiParticle("lambda0_bar", -13126);
	//| 1236 | S0(1665)   |
	fPdgDb->AddParticle("sigma0", "sigma0", 1.660, kFALSE, 0.100, 0, "Baryon", 13212);
	fPdgDb->AddAntiParticle("sigma0_bar", -13212);
	//| 1237 | S0(1776)   |
	fPdgDb->AddParticle("sigma0", "sigma0", 1.775, kFALSE, 0.120, 0, "Baryon", 3216);
	fPdgDb->AddAntiParticle("sigma0_bar", -3216);
	//| 1239 | S0(1930)   |
	fPdgDb->AddParticle("sigma0", "sigma0", 1.940, kFALSE, 0.220, 0, "Baryon", 23214);
	fPdgDb->AddAntiParticle("sigma0_bar", -23214);
	//| 1132 | S+(1665)   |
	fPdgDb->AddParticle("sigma+", "sigma+", 1.660, kFALSE, 0.100, 0, "Baryon", 13222);
	fPdgDb->AddAntiParticle("sigma+_bar", -13222);
	//| 1133 | S+(1776)   |
	fPdgDb->AddParticle("sigma+", "sigma+", 1.775, kFALSE, 0.120, 0, "Baryon", 3226);
	fPdgDb->AddAntiParticle("sigma+_bar", -3226);
	//| 1134 | S+(1930)   |
	fPdgDb->AddParticle("sigma+", "sigma+", 1.940, kFALSE, 0.220, 0, "Baryon", 23224);
	fPdgDb->AddAntiParticle("sigma+_bar", -23224);
	//| 2232 | S-(1665)   |
	fPdgDb->AddParticle("sigma-", "sigma-", 1.660, kFALSE, 0.100, 0, "Baryon", 13112);
	fPdgDb->AddAntiParticle("sigma-_bar", -13112);
	//| 2233 | S-(1776)   |
	fPdgDb->AddParticle("sigma-", "sigma-", 1.775, kFALSE, 0.120, 0, "Baryon", 3116);
	fPdgDb->AddAntiParticle("sigma-_bar", -3116);
	//| 2234 | S-(1930)   |
	fPdgDb->AddParticle("sigma-", "sigma-", 1.940, kFALSE, 0.220, 0, "Baryon", 23114);
	fPdgDb->AddAntiParticle("sigma-_bar", -23114);
	//+------+------------+

}
