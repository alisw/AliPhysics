#ifndef __CINT__
# include "AliFMDParameters.h"
# include "AliFMDCorrNoiseGain.h"
# include "AliForwardCorrectionManager.h"
# include "AliForwardUtil.h"
# include "AliCDBManager.h"
# include "AliCDBEntry.h"
# include "AliCDBStorage.h"
# include "AliGRPObject.h"
# include "AliMagF.h"
# include "AliLog.h"
# include <TSystem.h>
# include <TROOT.h>
# include <TError.h>
# include <TString.h>
#endif

/** 
 * Do one one 
 * 
 * @param runNo Run number 
 */
void 
ExtractForRun(Int_t runNo) 
{ 
  // --- Figure out the year --------------------------------------
  UShort_t year = 0;
  if      (runNo <= 99999)  year = 2009;
  else if (runNo <= 139667) year = 2010;
  else if (runNo <= 170718) year = 2011;
  else if (runNo <= 194306) year = 2012;
  else if (runNo <= 197709) year = 2013;
  if (year <= 0) { 
    Error("", "Couldn't deduce the year from the run number");
    return;
  }

  // --- Connect to OCDB ---------------------------------------------
  AliCDBManager* cdb = AliCDBManager::Instance();
  cdb->SetRun(runNo);
  cdb->SetDefaultStorageFromRun(runNo);

#if 0
  // --- Get the general run parameters ------------------------------
  // AliLog::SetModuleDebugLevel("STEER", 3);
  AliCDBId grpId("GRP/GRP/Data", runNo + 100, runNo - 100);
  AliCDBEntry* grpE = cdb->GetDefaultStorage()->GetEntry(grpId);
  if (!grpE) { 
    Warning("ExtractForRun", "No GRP entry found for run %d", runNo);
    return;
  }
  AliGRPObject* grp = static_cast<AliGRPObject*>(grpE->GetObject());
  if (!grp) { 
    Warning("ExtractForRun", "No GRP object found for run %d", runNo);
    return;
  }
  Float_t  beamE = grp->GetBeamEnergy();
  TString  beamT = grp->GetBeamType();
#if 0
  // This isn't really needed as the acceptance map is indifferent to
  // the field settings.
  Float_t  l3cur = grp->GetL3Current(AliGRPObject::kMean);
  Char_t   l3pol = grp->GetL3Polarity();
  Bool_t   l3lhc = grp->IsPolarityConventionLHC();
  Bool_t   l3uni = grp->IsUniformBMap();
  AliMagF* fldM  = 
    AliMagF::CreateFieldMap(TMath::Abs(l3cur) * (l3pol ? -1:1), 0, 
			    (l3lhc ? 0 : 1), l3uni, beamE, beamT.Data());
  Float_t  l3fld = fldM->SolenoidField();
#endif
  Printf("=== From GRP: Beam: E=%f T=%s", beamE, beamT.Data());
  if (beamE > 14000) beamE = 450;
  if (beamT.IsNull()) beamT = "pp";

  UShort_t sys = AliForwardUtil::ParseCollisionSystem(beamT);
  UShort_t sNN = AliForwardUtil::ParseCenterOfMassEnergy(sys, 2 * beamE);
  Short_t  fld = +999; // AliForwardUtil::ParseMagneticField(l3fld);
  Printf("=== Run=%d, year=%d, sys=%d, sNN=%d, fld=%d", 
	 runNo, year, sys, sNN, fld);
#endif
  UShort_t sys = +999;
  UShort_t sNN = +999;
  Short_t  fld = +999;

  // --- Get our parameters ------------------------------------------
  AliFMDParameters* param = AliFMDParameters::Instance();
  param->Init(true, AliFMDParameters::kPulseGain|AliFMDParameters::kPedestal);

  // --- Get the object to store -------------------------------------
  AliFMDCorrNoiseGain* ret   = new AliFMDCorrNoiseGain();
  Float_t              konst = param->GetDACPerMIP();

  // --- Loop over all strips ----------------------------------------
  for (UShort_t d = 1; d <= 3; d++) { 
    UShort_t nQ = (d == 1 ? 1 : 2);
    for (UShort_t q = 0; q < nQ; q++) { 
      Char_t   r  = (q == 0 ? 'I' : 'O');
      UShort_t nS = (q == 0 ?  20 :  40);
      UShort_t nT = (q == 0 ? 512 : 256);
      for (UShort_t s = 0; s < nS; s++) { 
	for (UShort_t t = 0; t < nT; t++) { 
	  Float_t noise = param->GetPedestalWidth(d,r,s,t);
	  Float_t gain  = param->GetPulseGain(d,r,s,t);
	  Float_t corr  = 0;
	  if (noise > .5 && gain > .5) corr = noise / (gain * konst);
	  if (corr > 1 || corr < 0) { 
	    Warning("", "FMD%d%c[%2d,%3d] corr= %f (=%f/(%f*%f))",
		    d, r, s, t, corr, noise, gain, konst);
	    corr = 0;
	  }
	  ret->Set(d,r,s,t,corr);
	}
      }
    }
  }

  // --- Write to a file ---------------------------------------------
  Printf("=== Writing to disk");
  AliForwardCorrectionManager& cm = AliForwardCorrectionManager::Instance();
  if (!cm.Store(ret, runNo, sys, sNN, fld, false, false, 
		"fmd_corrections.root", "OLDER")) { 
    Error("", "Failed to store acceptance correction in local file");
    return;
  }
  
}

void ExtractAll() {
  // We need to get a list of runs.  We should make an entry for every
  // pedestal and gain run.
  // 
  //  for y in `seq 2009 2013` ; do \
  //    for c in PulseGain Pedestal ; do \
  //      alien_ls /alice/data/${y}/OCDB/FMD/Calib/${c}/ | sed -e 's/Run//' -e 's/_.*//' ; \
  //    done ; \
  //  done | sort -u -n | grep -v ^0  | sed 's/\([0-9]*\)/    \1,/'
  // 
  Int_t runs[] = {
    // 58360,
    // 61383,
    // <-- Start of 2009
    75201,  75238,  75311,  75330,  75383,  75384,  75630,  75631,
    77901,  80649,  80650,  80738,  82055,  82301,  85947,  85948,
    87560,  87561,  91794,  91795,  92436,  92441,  93271,  93273,
    93587,  93588,  93595,  94793,  96747,  96945,  96947,  96949,
    96956,  96962,  96963,  97226,  97228,  97593,  97996,  98780,
    98782,  98974,  98980,  99033,  99084,  99085,  99414,  99726,
    100075, 100273, 100594, 100595, 100868, 100967, 101808, 101809, 
    102034, 102036, 102043, 102238, 102240, 103639, 103641, 103984, 
    103989, 104108, 104109, 104167, 104169, 104526, 104529, 104901, 
    104904, 104914, 105111, 105112,
    // <-- Start of 2010
    105827, 105834, 105963, 107718, 110372, 110373, 113268, 113331,
    113650, 113651, 114594, 114596, 115126, 115127, 115136, 115160,
    115171, 115244, 115334, 115355, 115360, 115477, 115480, 115528,
    115530, 115538, 115545, 115651, 115660, 116294, 116335, 116438,
    116439, 116440, 116613, 116614, 116655, 117255, 117257, 117258,
    117395, 117396, 117397, 117398, 117399, 117400, 117401, 117402,
    117403, 117404, 117405, 117406, 117407, 117408, 117409, 117411,
    117413, 117415, 117416, 117417, 117418, 117420, 117421, 117422,
    117423, 117424, 117425, 117426, 117427, 117428, 117429, 117431,
    117432, 117433, 117434, 117435, 117437, 117438, 117439, 117441,
    117442, 117443, 117444, 117445, 117446, 117447, 117448, 117450,
    117452, 117453, 117458, 117459, 117461, 117464, 117471, 117476,
    117478, 117483, 117492, 117495, 117496, 117498, 117512, 117524,
    117561, 117571, 117588, 117606, 117612, 117616, 117624, 117630,
    117648, 117650, 117664, 117677, 117678, 117680, 117681, 117682,
    117684, 117687, 117688, 117689, 117690, 117693, 117695, 117698,
    117699, 117700, 117701, 117702, 117703, 117704, 117706, 117707,
    117708, 117709, 117712, 117714, 117717, 117719, 117720, 117722,
    117723, 117724, 117725, 117726, 117728, 117730, 117731, 117733,
    117735, 117736, 117737, 117738, 117739, 117741, 117742, 117744,
    117745, 117746, 117747, 117763, 117772, 117776, 117779, 117784,
    119187, 119200, 120542, 120543, 121526, 121527, 121528, 121532,
    121550, 121554, 121557, 121621, 121638, 121651, 121656, 121662,
    121667, 121677, 121686, 121691, 121697, 121698, 121700, 121702,
    121704, 121706, 121708, 121710, 121712, 121713, 121714, 121715,
    121717, 121719, 121722, 121724, 121726, 121728, 121730, 121732,
    121734, 121736, 121738, 121740, 121742, 121743, 121963, 121965,
    121968, 121970, 121980, 121982, 121983, 121985, 121986, 121987,
    121988, 121989, 121990, 121991, 121992, 121993, 121994, 121995,
    121996, 121997, 121998, 121999, 122000, 122001, 122002, 122003,
    122004, 122005, 122006, 122007, 122008, 122009, 122010, 122011,
    122012, 122013, 124437, 124438, 125892, 125894, 125896, 125897,
    126043, 126046, 126939, 126940, 128138, 128140, 128511, 128512,
    129793, 129796, 129804, 130290, 130291, 131336, 131337, 131383,
    131386, 131387, 131389, 131392, 131396, 131397, 131398, 131399,
    131400, 132637, 132638, 132789, 132799, 135478, 135479, 136733,
    136738, 136740, 138243, 138244, 138935,
    //  <-- Start of 2011
    144429, 144751, 144837, 144856, 145167, 145169, 145170, 146647,
    146650, 147281, 147285, 147289, 147299, 147301, 147307, 147318,
    149252, 149256, 154053, 154054, 154057, 155422, 155424, 155430,
    155933, 155934, 156030, 156032, 156037, 156038, 156201, 156202,
    156204, 156207, 156217, 156222, 156778, 156780, 157800, 157801,
    157802, 157804, 157805, 157808, 158993, 158994, 165623, 165636,
    165637, 166817, 166819, 167219, 167221, 169438, 169443, 169446,
    169448, 169449, 169450, 169451, 169484, 169486,
    // <-- Start of 2012
    172968, 172970, 175768, 175769, 176615, 176616, 176617, 177999,
    178000, 178562, 178571, 179962, 179963, 182380, 182381, 183249,
    183253, 183594, 183596, 184904, 184911, 184912, 185266, 185267,
    185269, 185270, 185873, 186408, 186471, 187021, 187022, 187023,
    187259, 187263, 187803, 187804, 188300, 188301, 188851, 188856,
    189256, 189257, 191769, 191770,
    // <-- Start of 2013
    194507, 194522, 194523, 194526, 194589, 194590, 195027, 195030,
    // <-- End marker 
    -1 };
  Int_t* pRun = runs;
  Int_t  skipped = 0;
  Int_t  total   = 0;
  Int_t  last    = 0;
  while (*pRun > 0) { 
    Int_t next = *(pRun+1);
    Int_t dist = next - *pRun;
    total++;
    if (next > 0 && dist <= 20) {
      skipped++;
      pRun++;
      continue;
    }
#if 0
    if (last > 0) {
      dist = *pRun - last; 
      Printf("%-6d %s%d,%d,",dist, url,last,*pRun);
    }
#endif
    last = *pRun;
    ExtractForRun(last);
    pRun++;
  }
  Info("", "Skipped %d of %d", skipped, total);
}
// 
// EOF
// 
