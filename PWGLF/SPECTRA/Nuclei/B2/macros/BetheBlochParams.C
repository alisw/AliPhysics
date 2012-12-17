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

// macro for Bethe-Bloch parameters
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

void SetParameters(Double_t* param, Double_t c0, Double_t c1, Double_t c2, Double_t c3, Double_t c4)
{
//
// assign values to an array
//
	param[0] = c0;
	param[1] = c1;
	param[2] = c2;
	param[3] = c3;
	param[4] = c4;
}

void BetheBlochParams(Double_t* param, const TString& periodname)
{
//
// TPC Bethe Bloch ALEPH paramaters for different beam periods
//
	TString period = periodname;
	period.ToLower();
	
	if(period == "lhc10b_pass2")
	{
		SetParameters(param, 3.22422, 10.9345, 1.26309e-05, 2.26343, 2.43587);
	}
	else if(period == "lhc10b" || period == "lhc10b_pass3")
	{
		SetParameters(param, 5.24531, 5.82813, 0.000431522, 2.47198, 1.38539);
		
	}
	else if(period == "lhc10c900") // pass3
	{
		SetParameters(param, 1.4857, 22.9345, 1.77678e-11, 2.26456, 4.44306);
	}
	else if(period == "lhc10c_pass2")
	{
		SetParameters(param, 1.49726, 24.5879, 2.76442e-11, 2.15661, 4.91248);
	}
	else if(period == "lhc10c" || period == "lhc10c_pass3")
	{
		SetParameters(param, 7.41249, 5.13753, 0.000738319, 2.55495, 1.33519);
	}
	else if(period == "lhc10d" || period == "lhc10e" || period == "lhc10d_pass2" || period == "lhc10e_pass2")
	{
		SetParameters(param, 1.59526, 24.6438, 3.5082e-11, 2.18753, 3.7487);
	}
	else if(period == "lhc10h") // heavy ions
	{
		SetParameters(param, 2.77047, 14.6777, 5.62959e-08, 2.30422, 2.35623);
	}
	else if(period == "lhc11a" || period == "lhc11a_wsdd" || period == "lhc11a_wosdd") // pp 2.76 TeV
	{
		SetParameters(param, 4.94865, 8.29784, 9.95186e-05, 2.22417, 1.51139);
	}
	
	// simulation
	
	else if(period == "lhc10e12" || period == "lhc10e13" || period == "lhc10f6a" || period == "lhc10f6" || period == "lhc10e21")
	{
		SetParameters(param, 1.98509, 16.9132, 2.27954e-10, 2.1544, 3.94486);
	}
	else
	{
		SetParameters(param, 4.4194, 7.50931, 1.34e-05, 2.22085, 1.80461);
	}
}
