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

// average track multiplicity
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TString.h>
#endif

Double_t GetMeanNtrk08(const TString& period);
Double_t GetMeanNtrk05(const TString& period);
Double_t GetV0ANDMeanNtrk08(const TString& period);
Double_t GetV0ANDMeanNtrk05(const TString& period);

Double_t MeanNtrk(const TString& period, Double_t eta, Bool_t v0and)
{
//
// average track multiplicity for the given period and speudorapidity
//
	Double_t mntrk = 1;
	
	TString periodl = period;
	periodl.ToLower();
	
	if(TMath::Abs(eta) > 0.51)
	{
		mntrk = (v0and) ? GetV0ANDMeanNtrk08(periodl) : GetMeanNtrk08(periodl);
	}
	else
	{
		mntrk = (v0and) ? GetV0ANDMeanNtrk05(periodl) : GetMeanNtrk05(periodl);
	}
	
	if(mntrk == 1)
	{
		std::cerr << "Warning in MeanNtrk: no <Ntrk> for period " << period << " and |eta| < " << eta << std::endl;
	}
	
	return mntrk;
}

Double_t GetMeanNtrk08(const TString& period)
{
//
// average track multiplicity for |eta| < 0.8
//
	if(period == "lhc10b")   return 9.68887; // pass3
	if(period == "lhc10c")   return 9.66970; // pass3
	if(period == "lhc10d")   return 9.47466; // pass2
	if(period == "lhc10e")   return 9.55678; // pass2
	
	// MC
	if(period == "lhc10f6a") return 7.15259;
	if(period == "lhc10e21") return 7.69483;
	
	return 1;
}

Double_t GetMeanNtrk05(const TString& period)
{
//
// average track multiplicity for |eta| < 0.5
//
	if(period == "lhc10c900")    return 3.70158; // pass3
	if(period == "lhc10b_pass2") return 9.86258;
	if(period == "lhc10c_pass2") return 9.61402;
	if(period == "lhc10b")       return 5.96104; // pass3
	if(period == "lhc10c")       return 5.94719; // pass3
	if(period == "lhc10d")       return 5.82333; // pass2
	if(period == "lhc10e")       return 5.89367; // pass2
	if(period == "lhc11a_wosdd") return 4.28597; // pass3
	if(period == "lhc11a_wsdd")  return 4.69927; // pass4
	
	// MC
	if(period == "lhc10e13")            return 3.13712;
	if(period == "lhc10f6a")            return 4.41362;
	if(period == "lhc10e21")            return 4.74991;
	if(period == "lhc11e3a_plus_wosdd") return 3.37669;
	if(period == "lhc11e3a_plus_wsdd")  return 3.47885;
	
	if(period == "lhc12a5a")            return 29.264;
	if(period == "lhc12a5bb")           return 31.0288;
	if(period == "lhc12a5bc")           return 30.6888;
	if(period == "lhc12a5bd")           return 30.3528;
	if(period == "lhc12a5be")           return 29.9859;
	if(period == "lhc12a5c_wsdd")       return 27.5981;
	
	return 1;
}

Double_t GetV0ANDMeanNtrk08(const TString& period)
{
//
// average track multiplicity for |eta| < 0.8 (V0AND)
//
	if(period == "lhc10b")   return 10.0630; // pass3
	if(period == "lhc10c")   return 10.0292; // pass3
	if(period == "lhc10d")   return 9.77129; // pass2
	if(period == "lhc10e")   return 9.90511; // pass2
	
	// MC
	if(period == "lhc10f6a") return 7.5087;
	if(period == "lhc10e21") return 7.91423;
	
	return 1;
}

Double_t GetV0ANDMeanNtrk05(const TString& period)
{
//
// average track multiplicity for |eta| < 0.5 (V0AND)
//
	if(period =="lhc10c900")    return 3.84362; // pass3
	if(period =="lhc10b")       return 6.18470; // pass3
	if(period =="lhc10c")       return 6.16175; // pass3
	if(period =="lhc10d")       return 6.03108; // pass2
	if(period =="lhc10e")       return 6.10384; // pass2
	if(period =="lhc11a_wosdd") return 4.40312; // pass3
	if(period =="lhc11a_wsdd")  return 4.87609; // pass4
	
	// MC
	if(period =="lhc10e13")            return 3.33273;
	if(period =="lhc10f6a")            return 4.80771;
	if(period =="lhc10e21")            return 4.91967;
	if(period =="lhc11e3a_plus_wosdd") return 3.4774;
	if(period =="lhc11e3a_plus_wsdd")  return 3.6467;
	
	if(period =="lhc12a5a")            return 29.3414;
	if(period =="lhc12a5bb")           return 31.1514;
	if(period =="lhc12a5bc")           return 30.7877;
	if(period =="lhc12a5bd")           return 30.4706;
	if(period =="lhc12a5be")           return 30.1013;
	if(period =="lhc12a5c_wsdd")       return 27.6921;
	
	return 1;
}
