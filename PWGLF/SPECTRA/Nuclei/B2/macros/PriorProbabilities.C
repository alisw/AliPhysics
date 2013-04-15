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

// macro for prior probabilities
// author: Eulogio Serradilla <eulogio.serradilla@cern.ch>

void SetParameters(Double_t* prob, Double_t c0, Double_t c1, Double_t c2, Double_t c3, Double_t c4, Double_t c5, Double_t c6, Double_t c7, Double_t c8)
{
//
// assign values to an array
//
	prob[0] = c0;
	prob[1] = c1;
	prob[2] = c2;
	prob[3] = c3;
	prob[4] = c4;
	prob[5] = c5;
	prob[6] = c6;
	prob[7] = c7;
	prob[8] = c8;
}

Bool_t PriorProbabilities(Double_t* prob, const TString& periodname, const TString& tracksel, const TString& ztag="")
{
//
// Prior probabilities for the different samples and track selection criteria.
// Obtained from several iterations of the analysis program
//
	TString period = periodname;
	period.ToLower();
	
	TString trksel = tracksel;
	trksel.ToLower();
	
	if(period == "lhc10c900") // pp 900 GeV
	{
		if (trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0229385, 2.13852e-06, 0.868889, 0.0596281, 0.0484369, 0.00010232, 2.54977e-06, 0, 0);
		}
		else if (trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00670787, 1.03214e-06, 0.887492, 0.0627365, 0.0429736, 8.72157e-05, 2.06428e-06, 0, 0);
		}
		else if(trksel=="its_tpc_nsigma_spd")
		{
			SetParameters(prob, 0.00598811, 0.000283397, 0.889894, 0.0655886, 0.0382077, 3.76161e-05, 8.51043e-08, 0, 0);
		}
		else if(trksel == "its_tpc_tof_nsigma")
		{
			SetParameters(prob, 0.00968656, 0.00230198, 0.844375, 0.084224, 0.059334, 7.8127e-05, 7.89161e-07, 0, 0);
		}
		else if(trksel == "its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0123784, 0.00368129, 0.830011, 0.0795686, 0.0742296, 0.000127461, 2.99204e-06, 1.19682e-07, 0);
		}
		else // its_tpc_nsigma
		{
			SetParameters(prob, 0.0181447, 1.09966e-06, 0.880063, 0.0635145, 0.0382359, 4.06874e-05, 1.83277e-07, 0, 0);
		}
	}
	
	else if(   period.Contains("lhc10b")
	        || period.Contains("lhc10c")
	        || period.Contains("lhc10d")
	        || period.Contains("lhc10e") ) // pp 7 TeV
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0250199, 0.00330477, 0.856905, 0.0621491, 0.0524821, 0.000136071, 3.44981e-06, 6.86137e-08, 0);
		}
		else if(trksel=="its_tpc_dca_spd")
		{
			if(ztag=="ntrk0002")
			{
				SetParameters(prob, 0.00572332, 0.00235974, 0.888793, 0.060643, 0.0424185, 6.09002e-05, 1.29575e-06, 0, 0);
			}
			else if(ztag=="ntrk0204")
			{
				SetParameters(prob, 0.00605528, 0.00176134, 0.888974, 0.0614645, 0.0416846, 5.91396e-05, 1.27997e-06, 0, 0);
			}
			else if(ztag=="ntrk0408")
			{
				SetParameters(prob, 0.00679582, 0.00129216, 0.883082, 0.0641175, 0.044637, 7.44268e-05, 9.96841e-07, 4.80405e-08, 0);
			}
			else if(ztag=="ntrk0811")
			{
				SetParameters(prob, 0.00741765, 0.00105688, 0.87915, 0.0661869, 0.046099, 8.77512e-05, 1.61561e-06, 9.07646e-08, 0);
			}
			else if(ztag=="ntrk1120")
			{
				SetParameters(prob, 0.00806164, 0.000955738, 0.877048, 0.0674902, 0.0463403, 0.00010188, 1.96668e-06, 7.13305e-08, 0);
			}
			else if(ztag=="ntrk20xx")
			{
				SetParameters(prob, 0.00899995, 0.000887649, 0.876102, 0.0681763, 0.0457112, 0.00011968, 2.87651e-06, 1.16222e-07, 0);
			}
			else // MB
			{
				SetParameters(prob, 0.00883871, 0.00306025, 0.877528, 0.0646921, 0.0457766, 0.000102032, 2.23315e-06, 3.71013e-08, 0);
			}
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			if(ztag=="ntrk0002")
			{
				SetParameters(prob, 0.00354158, 0.0163535, 0.853117, 0.0703677, 0.0565492, 6.95016e-05, 1.68052e-06, 0, 0);
			}
			else if(ztag=="ntrk0204")
			{
				SetParameters(prob, 0.00383418, 0.0106071, 0.853307, 0.0739765, 0.0582004, 7.27388e-05, 1.58847e-06, 0, 0);
			}
			else if(ztag=="ntrk0408")
			{
				SetParameters(prob, 0.00457824, 0.00613158, 0.846351, 0.0795435, 0.0632952, 9.89517e-05, 1.19219e-06, 8.72334e-08, 0);
			}
			else if(ztag=="ntrk0811")
			{
				SetParameters(prob, 0.00517035, 0.00408579, 0.841406, 0.0833815, 0.0658357, 0.000118306, 2.28488e-06, 2.11563e-07, 0);
			}
			else if(ztag=="ntrk1120")
			{
				SetParameters(prob, 0.00577945, 0.00321126, 0.838286, 0.0861044, 0.0664752, 0.000140696, 3.25918e-06, 1.37712e-07, 0);
			}
			else if(ztag=="ntrk20xx")
			{
				SetParameters(prob, 0.00666829, 0.00269162, 0.836384, 0.0883033, 0.0657773, 0.000170996, 4.77082e-06, 1.88322e-07, 0);
			}
			else // MB
			{
				SetParameters(prob, 0.0155752, 0.0130319, 0.813626, 0.0808387, 0.0767456, 0.000177436, 4.32771e-06, 3.93428e-07, 0);
			}
		}
		else // its_tpc_nsigma
		{
			SetParameters(prob, 0.0207349, 0.000587565, 0.871052, 0.0673557, 0.0402021, 6.74649e-05, 0, 0, 0);
		}
	}
	
	else if(period == "lhc10h") // heavy ions 2.76 TeV, 20% centrality
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.032994, 0.012492, 0.824044, 0.0809969, 0.0492257, 0.000242368, 5.01681e-06, 3.34454e-07, 0);
		}
		else if(trksel=="its_tpc_dca_spd")
		{
			SetParameters(prob, 0.0262206, 0.00398172, 0.845586, 0.081878, 0.0421142, 0.000215343, 3.92153e-06, 0, 0);
		}
		else if(trksel=="its_tpc_nsigma")
		{
			SetParameters(prob, 0.0279673, 0.00329778, 0.845802, 0.085995, 0.0367557, 0.000179264, 2.62155e-06, 2.49671e-07, 1.24836e-07);
		}
		else if(trksel=="its_tpc_tof_dca")
		{
			SetParameters(prob, 0.017615, 0.00752739, 0.784779, 0.109672, 0.080031, 0.000370991, 3.95272e-06, 8.47012e-07, 0);
		}
		else if(trksel=="its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0130896, 0.0079917, 0.797236, 0.112833, 0.0685078, 0.000338538, 3.07762e-06, 0, 0);
		}
		else if(trksel=="its_tpc_tof_nsigma")
		{
			SetParameters(prob, 0.0157654, 0.00721323, 0.79879, 0.117158, 0.0607624, 0.000308557, 1.26199e-06, 6.30996e-07, 0);
		}
		else
		{
			SetParameters(prob, 0.032994, 0.012492, 0.824044, 0.0809969, 0.0492257, 0.000242368, 5.01681e-06, 3.34454e-07, 3.34454e-07);
		}
	}
	
	else if(period == "lhc11a_wosdd") // pp 2.76 TeV without SDD pass3
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0172658, 0.0179012, 0.857385, 0.0590938, 0.0482499, 0.000104021, 8.05072e-07, 0, 0);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00808143, 0.016243, 0.873576, 0.0593675, 0.0426588, 7.27021e-05, 6.45822e-07, 0, 0);
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0113244, 0.0275349, 0.790324, 0.0899881, 0.0806648, 0.000161212, 2.45189e-06, 0, 0);
		}
		else
		{
			SetParameters(prob, 0.0132164, 0.0209664, 0.864317, 0.0624972, 0.0389603, 4.30391e-05, 0, 0, 0);
		}
	}
	
	else if(period == "lhc11a_wsdd") // pp 2.76 TeV without SDD pass3
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0198007, 0.0130987, 0.857079, 0.0604184, 0.0494965, 0.000104874, 1.5572e-06, 0, 0);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00808143, 0.016243, 0.873576, 0.0593675, 0.0426588, 7.27021e-05, 6.45822e-07, 0, 0);
		}
		else if(trksel == "its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0206518, 0.00922775, 0.859646, 0.0595307, 0.0508038, 0.000137122, 3.02079e-06, 0, 0);
		}
		else
		{
			SetParameters(prob, 0.0150916, 0.0122592, 0.868485, 0.0649622, 0.0391601, 4.18627e-05, 0, 0, 0);
		}
	}
	
	// --------- simulation ----------
	
	else if(period == "lhc12a5a")
	{
		if(trksel == "its_tpc_dca")
		{
			SetParameters(prob, 0.00269435, 0.00142106, 0.0831012, 0.00522487, 0.371055, 0.332289, 0.0602906, 0.0830945, 0.0608297);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.000903954, 0.000740461, 0.0787924, 0.0052352, 0.377698, 0.331812, 0.0598412, 0.0834821, 0.0614951);
		}
		else if(trksel == "its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.000699209, 0.000628928, 0.0488789, 0.00274458, 0.378142, 0.35414, 0.0643921, 0.0830996, 0.0672741);
		}
		else
		{
			SetParameters(prob, 0.00177397, 0.000490646, 0.0745835, 0.00552979, 0.383725, 0.331756, 0.0592882, 0.0827403, 0.0601128);
		}
	}
	
	else if(period == "lhc12a5bb")
	{
		if(trksel == "its_tpc_dca")
		{
			SetParameters(prob, 0.00370442, 0.00158498, 0.108206, 0.00792265, 0.375296, 0.33667, 0.0610851, 0.0715099, 0.034021);
		}
		else if(trksel == "its_tpc_tof_dca")
		{
			SetParameters(prob, 0.00112885, 0.000832592, 0.0672591, 0.00469105, 0.385956, 0.361319, 0.0655949, 0.0767322, 0.0364865);
		}
	}
	
	else if(period == "lhc12a5bc")
	{
		if(trksel == "its_tpc_dca")
		{
			SetParameters(prob, 0.00370442, 0.00158498, 0.108206, 0.00792265, 0.375296, 0.33667, 0.0610851, 0.0715099, 0.034021);
		}
		else if(trksel == "its_tpc_tof_dca")
		{
			SetParameters(prob, 0.00112885, 0.000832592, 0.0672591, 0.00469105, 0.385956, 0.361319, 0.0655949, 0.0767322, 0.0364865);
		}
	}
	
	else if(period == "lhc12a5bd")
	{
		if(trksel == "its_tpc_dca")
		{
			SetParameters(prob, 0.00389124, 0.00181114, 0.10725, 0.00753775, 0.361361, 0.322078, 0.0584554, 0.0785471, 0.0590672);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00129393, 0.000963743, 0.102737, 0.00760901, 0.36841, 0.322256, 0.0581189, 0.0791033, 0.0595082);
		}
		else if(trksel == "its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.00122986, 0.000908924, 0.0696766, 0.00462655, 0.385381, 0.360127, 0.0653381, 0.0762508, 0.0364611);
		}
		else
		{
			SetParameters(prob, 0.00265834, 0.000648637, 0.0985374, 0.00801093, 0.373786, 0.322076, 0.0576089, 0.0782477, 0.0584268);
		}
	}
	
	else if(period == "lhc12a5be")
	{
		if(trksel == "its_tpc_tof_dca")
		{
			SetParameters(prob, 0.00114523, 0.000795738, 0.0666132, 0.00465889, 0.386032, 0.36171, 0.0656093, 0.0767938, 0.0366419);
		}
		else
		{
			SetParameters(prob, 0.00114523, 0.000795738, 0.0666132, 0.00465889, 0.386032, 0.36171, 0.0656093, 0.0767938, 0.0366419);
		}
	}
	
	else if(period == "lhc12a5c_wosdd" || period == "lhc12a5c_wsdd")
	{
		if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00140969, 0.00100911, 0.0901532, 0.00627527, 0.375878, 0.333126, 0.0601881, 0.0752193, 0.0567413);
		}
		else if(trksel == "its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.00254534, 0.00138188, 0.0920481, 0.00642957, 0.381436, 0.344603, 0.0627622, 0.0736729, 0.0351211);
		}
		else
		{
			SetParameters(prob, 0.00254534, 0.00138188, 0.0920481, 0.00642957, 0.381436, 0.344603, 0.0627622, 0.0736729, 0.0351211);
		}
	}
	
	else if(period == "lhc10e12")
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0274572, 0.0132276, 0.868083, 0.0579849, 0.0332055, 2.94442e-05, 1.1225e-05, 0, 8.63465e-07);
		}
		else if(trksel=="its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00902306, 0.00879859, 0.891372, 0.0600273, 0.0307466, 2.43343e-05, 7.44928e-06, 0, 4.96619e-07);
		}
		else if(trksel=="its_tpc_tof_dca")
		{
			SetParameters(prob, 0.0113533, 0.0104639, 0.877447, 0.0577014, 0.0429803, 3.77285e-05, 1.40853e-05, 0, 2.51523e-06);
		}
		else
		{
			SetParameters(prob, 0.0228473, 0.00644153, 0.881383, 0.061429, 0.0278978, 1.21828e-06, 2.81142e-07, 0, 0);
		}
	}
	
	else if(period == "lhc10e13")
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0275037, 0.0129471, 0.864088, 0.0566767, 0.03875, 2.46576e-05, 9.52441e-06, 0, 5.29134e-07);
		}
		else if(trksel=="its_tpc_dca_spd")
		{
			SetParameters(prob, 0.00881856, 0.00853512, 0.887761, 0.0587557, 0.03611, 1.46288e-05, 5.01096e-06, 0, 4.0411e-07);
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0118565, 0.0102582, 0.869917, 0.0579405, 0.0499923, 2.71412e-05, 8.44393e-06, 0, 0);
		}
		else
		{
			SetParameters(prob, 0.0227488, 0.00644463, 0.877399, 0.0598542, 0.033552, 1.08815e-06, 3.43625e-07, 0, 0);
		}
	}
	
	else if(period == "lhc10f6a")
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0301784, 0.0130575, 0.849464, 0.0631854, 0.0440839, 1.65739e-05, 8.28697e-06, 0, 5.52465e-06);
		}
		else if(trksel=="its_tpc_dca_spd")
		{
			SetParameters(prob, 0.0101971, 0.00853634, 0.875408, 0.0657018, 0.0401396, 1.18943e-05, 4.23147e-06, 0, 8.53405e-07);
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0149055, 0.0105359, 0.852053, 0.0663013, 0.0561727, 2.17209e-05, 9.18959e-06, 0, 8.35418e-07);
		}
		else
		{
			SetParameters(prob, 0.0248513, 0.00652446, 0.863878, 0.0670812, 0.0376654, 3.27516e-05, 2.45637e-05, 0, 5.4586e-06);
		}
	}
	
	else if(period == "lhc11e3a_wsdd") // from the grid without SDD
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0288269, 0.0143436, 0.85568, 0.0600814, 0.0410309, 2.04805e-05, 1.56616e-05, 0, 1.20474e-06);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.0104567, 0.00998477, 0.878832, 0.0622996, 0.0384064, 1.43839e-05, 5.60782e-06, 0, 7.60383e-07);
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0129161, 0.0117642, 0.852733, 0.0671011, 0.0554566, 2.17717e-05, 5.44293e-06, 0, 1.36073e-06);
		}
		else
		{
			SetParameters(prob, 0.0249858, 0.00773001, 0.87308, 0.0601418, 0.0340623, 0, 0, 0, 0);
		}
	}
	
	else if(period == "lhc11e3a_wosdd") // from the grid without SDD
	{
		if(trksel=="its_tpc_dca")
		{
			SetParameters(prob, 0.0222021, 0.0131403, 0.862641, 0.0613467, 0.0406294, 2.3342e-05, 1.55613e-05, 0, 1.29678e-06);
		}
		else if(trksel == "its_tpc_dca_spd")
		{
			SetParameters(prob, 0.0104567, 0.00998477, 0.878832, 0.0622996, 0.0384064, 1.43839e-05, 5.60782e-06, 0, 7.60383e-07);
		}
		else if(trksel=="its_tpc_tof_dca" || trksel == "its_tpc_tof_dca_spd")
		{
			SetParameters(prob, 0.0118237, 0.0110965, 0.855035, 0.0672175, 0.0547956, 2.49547e-05, 6.23867e-06, 0, 6.93185e-07);
		}
		else
		{
			SetParameters(prob, 0.0186686, 0.00736997, 0.875823, 0.062905, 0.0352281, 4.98982e-06, 0, 0, 0);
		}
	}
	
	else
	{
		std::cout << "Setting equal probabilities to all particle species: " << periodname << ", " << trksel << std::endl;
		
		SetParameters(prob, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11);
		
		return kFALSE;
	}
	
	return kTRUE;
}
