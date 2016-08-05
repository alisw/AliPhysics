#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliInputEventHandler.h"
#include "AliESDpid.h"
#include "AliVEvent.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODHeader.h"
#include "AliAODpidUtil.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliCentrality.h"

#include "AliAODForwardMult.h"
#include "AliForwardUtil.h"

#include "AliAnaCoulWave.h"
#include "AliAnalysisTaskPsi3HBTsystFitrange.h"

static const float PI = acos(-1.0);
static const float mpi = 0.139;
static const double M2Pion = 0.01947983612222;
static const float mK = 0.494;
static const float mp = 0.938272;


ClassImp(AliAnalysisTaskPsi3HBTsystFitrange)

	AliAnalysisTaskPsi3HBTsystFitrange::AliAnalysisTaskPsi3HBTsystFitrange()
  : AliAnalysisTaskSE()

	,fAOD						(0x0)
	,fOutputList		(0x0)

	,fHListEPCalib	(0x0)
	,fFileEPCalib		(0x0)

	,fCentralAna		(kFALSE)
	,fTrigSel				(-1)

	,fRunNumber							(0)
	,fRunNumberForPrevEvent	(-1)

	,fNEntries			(0x0)
	,fNumberE				(0x0)

	,cent_V0M				(0x0)
	,cent_Calib			(0x0)

	,fV0flag				(0x0)
	,fTPCflag				(0x0)
	,fFMDflag				(0x0)
	,fMagsign				(0x0)

	,fH1Event				(0x0)
	,fH1Nch					(0x0)

	,fH1Phi					(0x0)
	,fH1Eta					(0x0)
	,fH1Pt					(0x0)

	,fPrfV0Sig			(0x0)
	,fH2V0Sig				(0x0)

	,fFMD						(0x0)
	,fH2FMD2D				(0x0)

	,fH1kt					(0x0)
	,fH1kt_all			(0x0)

	,fH2DCA					(0x0)
	,fH2DCAxyPt			(0x0)
	,fH2DCAzPt 			(0x0)
	,fH2PhiEta			(0x0)

	,fH2TOF					(0x0)

	,fCoulomb				(0x0)

{
  // Constructor

	for (int ipart = 0; ipart < 3; ipart++) {
		fH2NsigTPC[ipart] = (0x0)	;
		fH2NsigTOF[ipart] = (0x0)	;
	}
	for (int itrig = 0; itrig < 4; itrig++) {
			fH1Cent[itrig]	= (0x0)	;
			fH1Zvtx[itrig]	= (0x0)	;
	}
  
	for (int iside = 0; iside < fTPCs; iside++) {
			fH1EtaD[iside]	= (0x0)	;
	}

	for (int ikt = 0; ikt < 4; ikt++) {
			fH1kt_div[ikt]	= (0x0)	;
	}
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfV0Qvcos		[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfV0Qvsin		[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfTPCQvcos	[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfTPCQvsin	[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfFMDQvcos	[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfFMDQvsin	[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop

	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1V0Qvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1V0Qvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1TPCQvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1TPCQvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1FMDQvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1FMDQvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPV0		[icor][itrig][iside][iharm][icent] = (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPTPC	[icor][itrig][iside][iharm][icent] = (0x0);
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPFMD	[icor][itrig][iside][iharm][icent] = (0x0);
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibV0Qvx		[itrig][iside][iharm] = (0x0)	;
				fEPCalibV0Qvy		[itrig][iside][iharm] = (0x0)	;
			}
		}
	}
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fTPCs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibTPCQvx	[itrig][iside][iharm] = (0x0)	;
				fEPCalibTPCQvy	[itrig][iside][iharm] = (0x0)	;
			}
		}
	}
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fFMDs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibFMDQvx	[itrig][iside][iharm] = (0x0)	;
				fEPCalibFMDQvy	[itrig][iside][iharm] = (0x0)	;
			}
		}
	}

	for (int idet = 0; idet < fDet; idet++) {
		for (int iharm = 0; iharm < fHarm; iharm++) {
			fPrfCos	[idet][iharm] = (0x0)	;
			fPrfSin	[idet][iharm] = (0x0)	;
		}
	}

	for (int idet = 0; idet < fVnDet; idet++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {

				fPrfvncos_cent	[idet][iside][iharm] = (0x0);
				fPrfvnsin_cent	[idet][iside][iharm] = (0x0);
				for (int icent = 0; icent < fVnCent; icent++) {
					fPrfvncos_pt	[idet][iside][iharm][icent] = (0x0);
					fPrfvnsin_pt	[idet][iside][iharm][icent] = (0x0);
				}
			}
		}
	}

	for (int icharge = 0; icharge < fChargeH; icharge++) {
		for (int icent = 0; icent < fCentH; icent++) {
			for (int iep = 0; iep < fPsi3H; iep++) {
				for (int iqbin = 0; iqbin < fQbin; iqbin++) {
					fH1Qinv			[icharge][icent][iep][iqbin]	= (0x0);
					fH1CQinv		[icharge][icent][iep][iqbin]	= (0x0);
					fH1Qinv_mix	[icharge][icent][iep][iqbin]	= (0x0);
					fH3Q				[icharge][icent][iep][iqbin]	= (0x0);
					fH3Q_mix		[icharge][icent][iep][iqbin]	= (0x0);
					fPrfConv		[icharge][icent][iep][iqbin]	= (0x0);
				}	//iqbin
			}	//iep
		}	//icent
	}	//icharge
}



//________________________________________________________________________
AliAnalysisTaskPsi3HBTsystFitrange::AliAnalysisTaskPsi3HBTsystFitrange(const char *name)
  : AliAnalysisTaskSE(name)

	,fAOD						(0x0)
	,fOutputList		(0x0)

	,fHListEPCalib	(0x0)
	,fFileEPCalib		(0x0)

	,fCentralAna		(kFALSE)
	,fTrigSel				(-1)

	,fRunNumber							(0)
	,fRunNumberForPrevEvent	(-1)

	,fNEntries			(0x0)
	,fNumberE				(0x0)

	,cent_V0M				(0x0)
	,cent_Calib			(0x0)

	,fV0flag				(0x0)
	,fTPCflag				(0x0)
	,fFMDflag				(0x0)
	,fMagsign				(0x0)

	,fH1Event				(0x0)
	,fH1Nch					(0x0)

	,fH1Phi					(0x0)
	,fH1Eta					(0x0)
	,fH1Pt					(0x0)

	,fPrfV0Sig			(0x0)
	,fH2V0Sig				(0x0)

	,fFMD						(0x0)
	,fH2FMD2D				(0x0)

	,fH1kt					(0x0)
	,fH1kt_all			(0x0)

	,fH2DCA					(0x0)
	,fH2DCAxyPt			(0x0)
	,fH2DCAzPt 			(0x0)
	,fH2PhiEta			(0x0)

	,fH2TOF					(0x0)

	,fCoulomb				(0x0)

{
  // Constructor

	for (int ipart = 0; ipart < 3; ipart++) {
		fH2NsigTPC[ipart] = (0x0)	;
		fH2NsigTOF[ipart] = (0x0)	;
	}

	for (int itrig = 0; itrig < 4; itrig++) {
			fH1Cent[itrig]	= (0x0)	;
			fH1Zvtx[itrig]	= (0x0)	;
	}
  
	for (int iside = 0; iside < fTPCs; iside++) {
			fH1EtaD[iside]	= (0x0)	;
	}

	for (int ikt = 0; ikt < 4; ikt++) {
			fH1kt_div[ikt]	= (0x0)	;
	}
  
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfV0Qvcos		[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfV0Qvsin		[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfTPCQvcos	[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfTPCQvsin	[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfFMDQvcos	[icor][itrig][iside][iharm]	= (0x0)	;
					fPrfFMDQvsin	[icor][itrig][iside][iharm]	= (0x0)	;
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop

	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1V0Qvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1V0Qvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1TPCQvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1TPCQvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1FMDQvcos	[icor][itrig][iside][iharm][icent]	= (0x0)	;
						fH1FMDQvsin	[icor][itrig][iside][iharm][icent]	= (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPV0		[icor][itrig][iside][iharm][icent] = (0x0)	;
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPTPC	[icor][itrig][iside][iharm][icent] = (0x0);
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						fH1EPFMD	[icor][itrig][iside][iharm][icent] = (0x0);
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibV0Qvx		[itrig][iside][iharm] = (0x0)	;
				fEPCalibV0Qvy		[itrig][iside][iharm] = (0x0)	;
			}
		}
	}
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fTPCs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibTPCQvx	[itrig][iside][iharm] = (0x0)	;
				fEPCalibTPCQvy	[itrig][iside][iharm] = (0x0)	;
			}
		}
	}
  
	for (int itrig = 0; itrig < fTrigH; itrig++) {
		for (int iside = 0; iside < fFMDs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fEPCalibFMDQvx	[itrig][iside][iharm] = (0x0)	;
				fEPCalibFMDQvy	[itrig][iside][iharm] = (0x0)	;
			}
		}
	}

	for (int idet = 0; idet < fDet; idet++) {
		for (int iharm = 0; iharm < fHarm; iharm++) {
			fPrfCos	[idet][iharm] = (0x0);
			fPrfSin	[idet][iharm] = (0x0);
		}
	}

	for (int idet = 0; idet < fVnDet; idet++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {

				fPrfvncos_cent	[idet][iside][iharm] = (0x0);
				fPrfvnsin_cent	[idet][iside][iharm] = (0x0);
				for (int icent = 0; icent < fVnCent; icent++) {
					fPrfvncos_pt	[idet][iside][iharm][icent] = (0x0);
					fPrfvnsin_pt	[idet][iside][iharm][icent] = (0x0);
				}
			}
		}
	}
	for (int icharge = 0; icharge < fChargeH; icharge++) {
		for (int icent = 0; icent < fCentH; icent++) {
			for (int iep = 0; iep < fPsi3H; iep++) {
				for (int iqbin = 0; iqbin < fQbin; iqbin++) {
					fH1Qinv			[icharge][icent][iep][iqbin]	= (0x0);
					fH1CQinv		[icharge][icent][iep][iqbin]	= (0x0);
					fH1Qinv_mix	[icharge][icent][iep][iqbin]	= (0x0);
					fH3Q				[icharge][icent][iep][iqbin]	= (0x0);
					fH3Q_mix		[icharge][icent][iep][iqbin]	= (0x0);
					fPrfConv		[icharge][icent][iep][iqbin]	= (0x0);
				}	//iqbin
			}	//iep
		}	//icent
	}	//icharge

  // Define input and output slots here
  // Input slot #0 works with a TChain
  //DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());

}


//________________________________________________________________________
void AliAnalysisTaskPsi3HBTsystFitrange::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

	srand((unsigned) time(NULL));
	cout << "start" << endl;

	for (int icharge = 0; icharge < fChargeH; icharge++) {
		for (int icent = 0; icent < fCentM; icent++ ) {
			for (int izvtx = 0; izvtx < fZvtxM; izvtx++ ) {
				for (int iep3 = 0; iep3 < fPsi3M; iep3++) {
					mix_flag[icharge][icent][izvtx][iep3] = 0;
				}
			}
		}
	}
	cout << "Mixing flag" << endl;

	const int fReserve[12] = {1000, 1000, 521, 435, 367, 308, 256, 211, 172, 138, 111, 89};

	for (int icent = 0; icent < fCentM; icent++ ) {
		for (int izvtx = 0; izvtx < fZvtxM; izvtx++ ) {
			for (int iep3 = 0; iep3 < fPsi3M; iep3++) {
				for (int imix = 0; imix < fMix; imix++) {
					pion_ppx	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mpx	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_ppy	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mpy	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_ppz	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mpz	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_ppt	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mpt	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_pe		[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_me		[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_pphi	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mphi	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_peta	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_meta	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_pclu	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_mclu	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_psha	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_msha	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
					pion_pax	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;   pion_max	[icent][izvtx][iep3][imix].reserve(fReserve[icent])		;
				}
			}
		}
	}

	fOutputList	= new TList();

	fH1Event		= new TH1F("fH1Eve", "TotalNevent", 5, 0., 5.);
	fH1Nch			= new TH1F("fH1Nch", "Number of charged particle", 200, 0.0, 5000.0);

	fH1Phi			= new TH1F("fH1Phi", "Phi dist", 100, 0.0, 10.0);
	fH1Eta			= new TH1F("fH1Eta", "Eta dist", 100, -2.0, 2.0);
	fH1Pt				= new TH1F("fH1Pt", "Pt_hist", 500, 0., 5.);

	fPrfV0Sig		= new TProfile("fPrfV0Sig_st2", "V0 Signal", 64, 0, 64, 0., 1000.);
	fH2V0Sig		= new TH2F("fH2V0Sig_st2", "V0 Signal", 64, 0, 64, 200, 0., 1000.);

	fH2FMD2D		= new TH2F("fH2FMD2D", "FMD 2D hist", 200, -4., 6., 20, 0., TMath::TwoPi());

	fH1kt				= new TH1F("fH1kt", "kt_hist", 300, 0., 3.)															;
	fH1kt_all		= new TH1F("fH1kt_all", "All kt_hist", 300, 0., 3.)											;

	fH2DCA			= new TH2F("fH2DCA", "DCA dist", 500, -5.0, 5.0, 500, -5.0, 5.0)				;
	fH2DCAxyPt	= new TH2F("fH2DCAxyPt",	"DCA vs Pt dist", 200, 0.0, 5.0, 500, -5.0, 5.0)	;
	fH2DCAzPt 	= new TH2F("fH2DCAzPt",		"DCA vs Pt dist", 200, 0.0, 5.0, 500, -5.0, 5.0)	;
	fH2PhiEta		= new TH2F("fH2PhiEta", "Phi v.s. Eta", 100, 0.0, 10.0, 100, -2.0, 2.0)	;

	fH2TOF			= new TH2F("fH2TOF", "Pt v.s. TOFbeta", 500, 0.0, 5.0, 1000, -0.4, 1.1)	;

	fCoulomb		= new AliAnaCoulWave();
	fCoulomb		-> Set_Qtype(0);

	Char_t det_name		[2]	[20]	= {"Before", "After"};
	Char_t side_name	[3]	[20]	= {"C", "A", "A+C"};
	Char_t side_name2	[5]	[30]	= {"C side(-1, -0.5)", "C side(-0.5, 0.0)", "A side(0.0, 0.5)", "A side(0.5, 1.0)"};
	Char_t side_name3	[13][30]	= {"C1", "C2", "C3", "C", "A1", "A2", "A3", "A4", "A5", "A", "A123", "A45", "AC"};
	Char_t trig_name	[4]	[20]	= {"MB", "Semicent_trig", "Central_trig", "ALL"};
	Char_t cent_name	[20][20]	= {"0-5%", "5-10%", "10-15%", "15-20%", "20-25%", "25-30%", "30-35%", "35-40%", "40-45%", "40-45%", "45-50%", "50-55%", "55-60%", "60-65%", "65-70%", "70-75%", "80-85%", "85-90%", "90-95%", "95-100%" };
	Char_t cent_name2	[9]	[20]	= {"0-10%", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "40-95%"};
	Char_t zvt_name		[8]	[20]	= {"(-8,-6)", "(-6,-4)", "(-4,-2)", "(-2, 0)", "(0, 2)", "(2, 4)", "(4, 6)", "(6, 8)"};
	Char_t part_name	[3]	[20]	= {"pion nsigma", "kaon nsigma", "proton nsigma"};
	Char_t h_name			[20][100]	;

	for (int ipart = 0; ipart < 3; ipart++) {
		sprintf(h_name[0], "fH2NsigTPC%d", ipart+1);
		sprintf(h_name[1], "fH2NsigTOF%d", ipart+1);

		fH2NsigTPC[ipart] = new TH2F(h_name[0], part_name[ipart], 500, 0.0, 5.0, 500, -10., 10.);
		fH2NsigTOF[ipart] = new TH2F(h_name[1], part_name[ipart], 500, 0.0, 5.0, 500, -10., 10.);
	}

	for (int itrig = 0; itrig < 4; itrig++) {
			sprintf(h_name[0], "fH1Cent%d", itrig+1);
			sprintf(h_name[1], "fH1ZVtx%d", itrig+1);
			fH1Cent[itrig]	= new TH1F(h_name[0], trig_name[itrig], 100, 0., 100.);
			fH1Zvtx[itrig]	= new TH1F(h_name[1], trig_name[itrig], 100, -10.0, 10.0)		;
	}
  
	for (int iside = 0; iside < fTPCs; iside++) {
			sprintf(h_name[0], "fH1EtaD%d", iside+1);
			fH1EtaD[iside]	= new TH1F(h_name[0], "Eta dist", 100, -2.0, 2.0);
	}

	for (int ikt = 0; ikt < 4; ikt++) {
		sprintf(h_name[0], "fH1kt_div%d",		ikt+1);
		fH1kt_div			[ikt]		= new TH1F(h_name[0], "kt_hist", 300, 0., 3.);
	}
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
	
					sprintf(h_name[2],	"fPrfV0Qvcos%d%d%d%d",	icor+1, itrig+1, iside+1, iharm+1);
					sprintf(h_name[3],	"fPrfV0Qvsin%d%d%d%d",	icor+1, itrig+1, iside+1, iharm+1);
	
					sprintf(h_name[8],	"%s %s QvxV0 %s side Psi%d",	det_name[icor], trig_name[itrig], side_name[iside], iharm+1);
					sprintf(h_name[9],	"%s %s QvyV0 %s side Psi%d",	det_name[icor], trig_name[itrig], side_name[iside], iharm+1);
	
					fPrfV0Qvcos	[icor][itrig][iside][iharm]	= new TProfile(h_name[2], h_name[8],	20, 0, 100, "S");
					fPrfV0Qvsin	[icor][itrig][iside][iharm]	= new TProfile(h_name[3], h_name[9],	20, 0, 100, "S");
  
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
		
					sprintf(h_name[6], "fPrfTPCQvcos%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1);
					sprintf(h_name[7], "fPrfTPCQvsin%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1);
	
					sprintf(h_name[8], "%s %s PrfQvcosTPC%s Psi%d", det_name[icor], trig_name[itrig], side_name2[iside], iharm+1);
					sprintf(h_name[9], "%s %s PrfQvsinTPC%s Psi%d", det_name[icor], trig_name[itrig], side_name2[iside], iharm+1);
	
					fPrfTPCQvcos[icor][itrig][iside][iharm]	= new TProfile(h_name[6], h_name[8], 20, 0, 100, "S");
					fPrfTPCQvsin[icor][itrig][iside][iharm]	= new TProfile(h_name[7], h_name[9], 20, 0, 100, "S");
  
				}	//iside loop
			}	//iharm loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
	
					sprintf(h_name[6], "fPrfFMDQvcos%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1);
					sprintf(h_name[7], "fPrfFMDQvsin%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1);
	
					sprintf(h_name[8], "%s %s PrfQvcosFMD%s Psi%d", det_name[icor], trig_name[itrig], side_name3[iside], iharm+1);
					sprintf(h_name[9], "%s %s PrfQvsinFMD%s Psi%d", det_name[icor], trig_name[itrig], side_name3[iside], iharm+1);
	
					fPrfFMDQvcos[icor][itrig][iside][iharm]	= new TProfile(h_name[6], h_name[8], 20, 0, 100, "S");
					fPrfFMDQvsin[icor][itrig][iside][iharm]	= new TProfile(h_name[7], h_name[9], 20, 0, 100, "S");

				}	//iside loop
			}	//iharm loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
  
						sprintf(h_name[2], "fH1V0Qvcos%d%d%d%d%d",	icor+1, itrig+1, iside+1, iharm+1, icent+1);
						sprintf(h_name[3], "fH1V0Qvsin%d%d%d%d%d",	icor+1, itrig+1, iside+1, iharm+1, icent+1);
  
						sprintf(h_name[8],	"%s QvcosV0%s Psi%d %s %s",		trig_name[itrig], side_name[iside], iharm+1, det_name[icor], cent_name2[icent]);
						sprintf(h_name[9],	"%s QvsinV0%s Psi%d %s %s",		trig_name[itrig], side_name[iside], iharm+1, det_name[icor], cent_name2[icent]);
  
						fH1V0Qvcos	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[2], h_name[8],	100, -10., 10.);
						fH1V0Qvsin	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[3], h_name[9],	100, -10., 10.);

					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
  
						sprintf(h_name[0], "fH1TPCQvcos%d%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1, icent+1);
						sprintf(h_name[1], "fH1TPCQvsin%d%d%d%d%d", icor+1, itrig+1, iside+1, iharm+1, icent+1);
  
						sprintf(h_name[2], "%s QvcosTPC%s Psi%d %s %s", trig_name[itrig], side_name2[iside], iharm+1, det_name[icor], cent_name2[icent]);
						sprintf(h_name[3], "%s QvsinTPC%s Psi%d %s %s", trig_name[itrig], side_name2[iside], iharm+1, det_name[icor], cent_name2[icent]);
  
						fH1TPCQvcos	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[0], h_name[2], 100, -10., 10.);
						fH1TPCQvsin	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[1], h_name[3], 100, -10., 10.);

					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
  
						sprintf(h_name[0], "fH1FMDQvcos%d_%d_%d_%d_%d", icor+1, itrig+1, iside+1, iharm+1, icent+1);
						sprintf(h_name[1], "fH1FMDQvsin%d_%d_%d_%d_%d", icor+1, itrig+1, iside+1, iharm+1, icent+1);
  
						sprintf(h_name[2], "%s QvcosFMD%s Psi%d %s %s", trig_name[itrig], side_name3[iside], iharm+1, det_name[icor], cent_name2[icent]);
						sprintf(h_name[3], "%s QvsinFMD%s Psi%d %s %s", trig_name[itrig], side_name3[iside], iharm+1, det_name[icor], cent_name2[icent]);
  
						fH1FMDQvcos	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[0], h_name[2], 100, -10., 10.);
						fH1FMDQvsin	[icor][itrig][iside][iharm][icent]	= new TH1F(h_name[1], h_name[3], 100, -10., 10.);
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	}	//icor loop
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						sprintf(h_name[1], "fH1EPV0%d%d%d%d%d",		icor+1, itrig+1, iside+1, iharm+1, icent+1);
  
						fH1EPV0		[icor][itrig][iside][iharm][icent] = new TH1F(h_name[1], h_name[1], 50, -TMath::Pi(), TMath::Pi());
					}
				}
			}
		}
	}
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						sprintf(h_name[0], "fH1EPTPC%d%d%d%d%d",	icor+1, itrig+1, iside+1, iharm+1, icent+1);
						fH1EPTPC	[icor][itrig][iside][iharm][icent] = new TH1F(h_name[0], h_name[0], 50, -TMath::Pi(), TMath::Pi());
					}
				}
			}
		}
	}
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fQvCent; icent++) {
						sprintf(h_name[0], "fH1EPFMD%d_%d_%d_%d_%d",	icor+1, itrig+1, iside+1, iharm+1, icent+1);
						fH1EPFMD	[icor][itrig][iside][iharm][icent] = new TH1F(h_name[0], h_name[0], 50, -TMath::Pi(), TMath::Pi());
					}
				}
			}
		}
	}

	for (int idet = 0; idet < fDet; idet++) {
		for (int iharm = 0; iharm < fHarm; iharm++) {
			sprintf(h_name[1], "fPrfCos%d_%d",		idet+1, iharm+1);
			sprintf(h_name[2], "fPrfSin%d_%d",		idet+1, iharm+1);

			fPrfCos		[idet][iharm] = new TProfile(h_name[1], h_name[1], 18, 0.0, 90.0);
			fPrfSin		[idet][iharm] = new TProfile(h_name[2], h_name[2], 18, 0.0, 90.0);
		}	//iharm loop
	}	//idet loop

	Char_t det_name2[3][10] = {"TPC", "V0", "FMD"};
	int nxbin = 31;
	Double_t xbins[32] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.4, 2.8, 3.2, 3.6, 4.0, 4.4, 4.8, 5.2, 5.6, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 14.0}; 

	for (int idet = 0; idet < fVnDet; idet++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				sprintf(h_name[0], "fPrfvncos_cent%d%d%d", idet+1, iside+1, iharm+1);
				sprintf(h_name[1], "fPrfvnsin_cent%d%d%d", idet+1, iside+1, iharm+1);

				sprintf(h_name[2], "v%d cos %s %sside ",	iharm+1, det_name2[idet], side_name[iside]);
				sprintf(h_name[3], "v%d sin %s %sside ",	iharm+1, det_name2[idet], side_name[iside]);

				fPrfvncos_cent	[idet][iside][iharm] = new TProfile(h_name[0], h_name[2], 18, 0.0, 90.0);
				fPrfvnsin_cent	[idet][iside][iharm] = new TProfile(h_name[1], h_name[3], 18, 0.0, 90.0);
				for (int icent = 0; icent < fVnCent; icent++) {
					sprintf(h_name[4], "fPrfvncos_pt%d%d%d%d", idet+1, iside+1, iharm+1, icent+1);
					sprintf(h_name[5], "fPrfvnsin_pt%d%d%d%d", idet+1, iside+1, iharm+1, icent+1);

					sprintf(h_name[6], "v%d cos %s %sside %s",	iharm+1, det_name2[idet], side_name[iside], cent_name[icent]);
					sprintf(h_name[7], "v%d sin %s %sside %s",	iharm+1, det_name2[idet], side_name[iside], cent_name[icent]);
					fPrfvncos_pt	[idet][iside][iharm][icent] = new TProfile(h_name[4], h_name[6], nxbin, xbins, "");
					fPrfvnsin_pt	[idet][iside][iharm][icent] = new TProfile(h_name[5], h_name[7], nxbin, xbins, "");
				}
			}
		}
	}

	for (int icharge = 0; icharge < fChargeH; icharge++) {
		for (int icent = 0; icent < fCentH; icent++) {
			for (int iep = 0; iep < fPsi3H; iep++) {
				for (int iqbin = 0; iqbin < fQbin; iqbin++) {

					sprintf(h_name[0], "fH1Qinv%d%d%d%d"			,icharge+1 ,icent+1, iep+1, iqbin+1);
					sprintf(h_name[1], "fH1CQinv%d%d%d%d"			,icharge+1 ,icent+1, iep+1, iqbin+1);
					sprintf(h_name[2], "fH1Qinv_mix%d%d%d%d"	,icharge+1 ,icent+1, iep+1, iqbin+1);
					sprintf(h_name[3], "fH3Q%d%d%d%d"					,icharge+1 ,icent+1, iep+1, iqbin+1);
					sprintf(h_name[4], "fH3Q_mix%d%d%d%d"			,icharge+1 ,icent+1, iep+1, iqbin+1);
					sprintf(h_name[5], "fPrfConv%d%d%d%d"			,icharge+1 ,icent+1, iep+1, iqbin+1);

					if (!iqbin)	{
						fH1Qinv			[icharge][icent][iep][iqbin]	= new TH1D(h_name[0], h_name[0], 34, 0.0, 0.34);
						fH1CQinv		[icharge][icent][iep][iqbin]	= new TH1D(h_name[1], h_name[1], 34, 0.0, 0.34);
						fH1Qinv_mix	[icharge][icent][iep][iqbin]	= new TH1D(h_name[2], h_name[2], 34, 0.0, 0.34);
						fH3Q				[icharge][icent][iep][iqbin]	= new TH3D(h_name[3], h_name[3], 60, -0.3, 0.3, 60, -0.3, 0.3, 60, -0.3, 0.3);
						fH3Q_mix		[icharge][icent][iep][iqbin]	= new TH3D(h_name[4], h_name[4], 60, -0.3, 0.3, 60, -0.3, 0.3, 60, -0.3, 0.3);
						fPrfConv		[icharge][icent][iep][iqbin]	= new TProfile(h_name[5],h_name[5],64000, 0, 64000, 0., 1.0);
					}
					else if (iqbin==2) {
						fH1Qinv			[icharge][icent][iep][iqbin]	= new TH1D(h_name[0], h_name[0], 17, 0.0, 0.34);
						fH1CQinv		[icharge][icent][iep][iqbin]	= new TH1D(h_name[1], h_name[1], 17, 0.0, 0.34);
						fH1Qinv_mix	[icharge][icent][iep][iqbin]	= new TH1D(h_name[2], h_name[2], 17, 0.0, 0.34);
						fH3Q				[icharge][icent][iep][iqbin]	= new TH3D(h_name[3], h_name[3], 30, -0.3, 0.3, 30, -0.3, 0.3, 30, -0.3, 0.3);
						fH3Q_mix		[icharge][icent][iep][iqbin]	= new TH3D(h_name[4], h_name[4], 30, -0.3, 0.3, 30, -0.3, 0.3, 30, -0.3, 0.3);
						fPrfConv		[icharge][icent][iep][iqbin]	= new TProfile(h_name[5],h_name[5],8000, 0, 8000, 0., 1.0);
					}
					else {
						fH1Qinv			[icharge][icent][iep][iqbin]	= new TH1D(h_name[0], h_name[0], 68, 0.0, 0.34);
						fH1CQinv		[icharge][icent][iep][iqbin]	= new TH1D(h_name[1], h_name[1], 68, 0.0, 0.34);
						fH1Qinv_mix	[icharge][icent][iep][iqbin]	= new TH1D(h_name[2], h_name[2], 68, 0.0, 0.34);
						//fH3Q				[icharge][icent][iep][iqbin]	= new TH3D(h_name[3], h_name[3], 80, -0.20, 0.20, 80, -0.20, 0.20, 80, -0.20, 0.20);
						//fH3Q_mix		[icharge][icent][iep][iqbin]	= new TH3D(h_name[4], h_name[4], 80, -0.20, 0.20, 80, -0.20, 0.20, 80, -0.20, 0.20);
						//fPrfConv		[icharge][icent][iep][iqbin]	= new TProfile(h_name[5],h_name[5],512000, 0, 512000, 0., 1.0);
						fH3Q				[icharge][icent][iep][iqbin]	= new TH3D(h_name[3], h_name[3], 60, -0.15, 0.15, 60, -0.15, 0.15, 60, -0.15, 0.15);
						fH3Q_mix		[icharge][icent][iep][iqbin]	= new TH3D(h_name[4], h_name[4], 60, -0.15, 0.15, 60, -0.15, 0.15, 60, -0.15, 0.15);
						fPrfConv		[icharge][icent][iep][iqbin]	= new TProfile(h_name[5],h_name[5],216000, 0, 216000, 0., 1.0);
					}
				}	//iqbin
			}	//iep
		}	//icent
	}	//icharge

	fOutputList->Add(	fH1Event	)					;
	fOutputList->Add(	fH1Nch		)					;
  
	fOutputList->Add(	fH1Phi		)					;
	fOutputList->Add(	fH1Eta		)					;
	fOutputList->Add(	fH1Pt			)					;
  
	fOutputList->Add(	fPrfV0Sig	)					;
	fOutputList->Add(	fH2V0Sig	)					;
  
	fOutputList->Add(	fH2FMD2D	)					;

	fOutputList->Add(	fH1kt				)				;
	fOutputList->Add(	fH1kt_all		)				;
  
	fOutputList->Add(	fH2DCA			)				;
	fOutputList->Add(	fH2DCAxyPt	)				;
	fOutputList->Add(	fH2DCAzPt		)				;
	fOutputList->Add(	fH2TOF			)				;
	fOutputList->Add(	fH2PhiEta		)				;
  
  
	for (int itrig = 0; itrig < 4; itrig++) {
		fOutputList->Add(	fH1Cent[itrig]	);
		fOutputList->Add(	fH1Zvtx[itrig]	);
	}
  //
	//for (int iside = 0; iside < fTPCs; iside++) {
	//	fOutputList->Add(	fH1EtaD[iside]	)	;
	//}
  //
	//for (int ikt = 0; ikt < 4; ikt++) {
	//	fOutputList->Add(	fH1kt_div	[ikt]	)	;
	//}
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fV0s; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				fOutputList -> Add(fPrfV0Qvcos	[icor][itrig][iside][iharm]);
	//				fOutputList -> Add(fPrfV0Qvsin	[icor][itrig][iside][iharm]);
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fTPCs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				fOutputList -> Add(fPrfTPCQvcos[icor][itrig][iside][iharm]);
	//				fOutputList -> Add(fPrfTPCQvsin[icor][itrig][iside][iharm]);
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fFMDs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				fOutputList -> Add(fPrfFMDQvcos[icor][itrig][iside][iharm]);
	//				fOutputList -> Add(fPrfFMDQvsin[icor][itrig][iside][iharm]);
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fV0s; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList -> Add(fH1V0Qvcos	[icor][itrig][iside][iharm][icent]);
	//					fOutputList -> Add(fH1V0Qvsin	[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fTPCs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList -> Add(fH1TPCQvcos[icor][itrig][iside][iharm][icent]);
	//					fOutputList -> Add(fH1TPCQvsin[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fFMDs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList -> Add(fH1FMDQvcos[icor][itrig][iside][iharm][icent]);
	//					fOutputList -> Add(fH1FMDQvsin[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iharm loop
	//		}	//iside loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fV0s; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList->Add(fH1EPV0	[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iside loop
	//		}	//iharm loop
	//	}	//itrig loop
	//}	//icor loop
  //
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fTPCs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList->Add(fH1EPTPC	[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iside loop
	//		}	//iharm loop
	//	}	//itrig loop
	//}	//icor loop
  
	//for (int icor = 0; icor < fCorr; icor++) {
	//	for (int itrig = 0; itrig < fTrigH; itrig++) {
	//		for (int iside = 0; iside < fFMDs; iside++) {
	//			for (int iharm = 0; iharm < fHarm; iharm++) {
	//				for (int icent = 0; icent < fQvCent; icent++) {
	//					fOutputList->Add(fH1EPFMD	[icor][itrig][iside][iharm][icent]);
	//				}	//icent loop
	//			}	//iside loop
	//		}	//iharm loop
	//	}	//itrig loop
	//}	//icor loop

	for (int ipart = 0; ipart < 3; ipart++) {
		fOutputList	->	Add(fH2NsigTPC[ipart]);
		fOutputList	->	Add(fH2NsigTOF[ipart]);
	}

	for (int idet = 0; idet < fDet; idet++) {
		for (int iharm = 0; iharm < fHarm; iharm++) {

			fOutputList	->	Add(fPrfCos		[idet][iharm] );
			fOutputList	->	Add(fPrfSin		[idet][iharm] );
		}	//iharm loop
	}	//idet loop

	for (int idet = 0; idet < fVnDet; idet++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				fOutputList	->	Add(fPrfvncos_cent	[idet][iside][iharm]);
				fOutputList	->	Add(fPrfvnsin_cent	[idet][iside][iharm]);
				for (int icent = 0; icent < fVnCent; icent++) {
					fOutputList	->	Add(fPrfvncos_pt	[idet][iside][iharm][icent]);
					fOutputList	->	Add(fPrfvnsin_pt	[idet][iside][iharm][icent]);
				}	//icent loop
			}	//iharm loop
		}	//iside loop
	}	//idet loop

	for (int icharge = 0; icharge < fChargeH; icharge++) {
		for (int icent = 0; icent < fCentH; icent++) {
			for (int iep = 0; iep < fPsi3H; iep++) {
				for (int iqbin = 0; iqbin < fQbin; iqbin++) {

					fOutputList -> Add(fH1Qinv		[icharge][icent][iep][iqbin]);
					fOutputList -> Add(fH1CQinv		[icharge][icent][iep][iqbin]);
					fOutputList -> Add(fH1Qinv_mix[icharge][icent][iep][iqbin]);
					fOutputList -> Add(fH3Q				[icharge][icent][iep][iqbin]);
					fOutputList -> Add(fH3Q_mix		[icharge][icent][iep][iqbin]);
					fOutputList -> Add(fPrfConv		[icharge][icent][iep][iqbin]);
				}	//iqbin
			}	//iep
		}	//icent
	}	//icharge
}

//________________________________________________________________________
void AliAnalysisTaskPsi3HBTsystFitrange::UserExec(Option_t *)
{

	//Coulom source size Input
	Float_t r_inv		[6]	= {11, 10, 9, 8, 7, 6};
	Float_t r_side	[6]	= {11, 10, 9, 8, 7, 6};
	Float_t r_out		[6]	= {11, 10, 9, 8, 7, 6};
	Float_t r_long	[6]	= {11, 10, 9, 8, 7, 6};
	//Float_t r_inv		[5]	= {9, 8, 7, 6, 5};
	//Float_t r_side	[5]	= {9, 8, 7, 6, 5};
	//Float_t r_out		[5]	= {9, 8, 7, 6, 5};
	//Float_t r_long	[5]	= {9, 8, 7, 6, 5};

	// Main loop
	// Called for each event
  
	// Post output data.
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  
	if (!fAOD) {
					Printf("ERROR: fAOD not available");
					PairCalc_Final();
					return;
	}
	fH1Event -> Fill(0);
  
	fRunNumber = fAOD->GetRunNumber();
  
	TTree *fTree_temp = AliAnalysisManager::GetAnalysisManager()->GetTree();                                                                                                                                                                                                      
	fNEntries = fTree_temp->GetEntries();
	//delete fTree_temp;
  
	fNumberE = AliAnalysisManager::GetAnalysisManager()->GetNcalls();

	cout << "-------------------------------------" << endl;
	cout << "Event : " << fNumberE << " / " << fNEntries << endl;
	cout << "-------------------------------------" << endl;



	Char_t h_name			[20][100]	;
  
	//Load Event Plane Calibration Parameter
	//if Run Number is not same for previous LOOP
	if( fRunNumberForPrevEvent!=fRunNumber) {
  
		fRunNumberForPrevEvent=fRunNumber;
  
		TString path = Form("alien:///alice/cern.ch/user/n/ntanaka/run11/EPCalib_AOD145/Recentering/Param/000%d/EPCalib.AOD.root",fRunNumber);
  
		if(fFileEPCalib){
			fFileEPCalib->Close()	;
			delete fFileEPCalib		;
		}
  
		fFileEPCalib=0x0;
		fFileEPCalib = TFile::Open(path.Data(),"READ");
		//fFileEPCalib = TFile::Open("EPCalib.AOD.167915.root","READ");
  
		if ( !fFileEPCalib ) {
			printf("Error : %s is NOT Exist!!! \n",path.Data());
			cout << endl << endl << endl;
			PairCalc_Final();
			return;
		}
		cout << "Calibration file read !" << endl;
  
		fHListEPCalib=0x0;
		fHListEPCalib = (TList*)fFileEPCalib->Get("EPcaliblist_st1");
		if ( !fHListEPCalib ) {
			printf("Error : TList EPcaliblist_CalibStage2 is NOT Exist in %s!!! \n", path.Data());
			cout << endl << endl << endl;
			PairCalc_Final();
			return;
		}
  
		//fHListEPCalib->Print();
  
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s ; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
	
					fEPCalibV0Qvx	[itrig][iside][iharm]	= 0x0;	fEPCalibV0Qvy	[itrig][iside][iharm]	= 0x0;
	
					sprintf(h_name[2], "fPrfV0Qvcos%d%d%d",		itrig+1,	iside+1, iharm+1);
					sprintf(h_name[3], "fPrfV0Qvsin%d%d%d",		itrig+1,	iside+1, iharm+1);
	
					fEPCalibV0Qvx	[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[2]);
					fEPCalibV0Qvy	[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[3]);
	
					if (!fEPCalibV0Qvx	[itrig][iside][iharm])		{printf("Error : Couldn't Get Recentering Parameters(x) for harm=%d VZERO\n",	iharm+1);}
					if (!fEPCalibV0Qvy	[itrig][iside][iharm])		{printf("Error : Couldn't Get Recentering Parameters(y) for harm=%d VZERO\n",	iharm+1);}

				}	//iharm loop
			}	//iside loop
		}	//itrig loop
	
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs ; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
	
					fEPCalibTPCQvx[itrig][iside][iharm]	= 0x0;	fEPCalibTPCQvy[itrig][iside][iharm]	= 0x0;
	
					sprintf(h_name[0], "fPrfTPCQvcos%d%d%d", itrig+1, iside+1, iharm+1);
					sprintf(h_name[1], "fPrfTPCQvsin%d%d%d", itrig+1, iside+1, iharm+1);
	
					fEPCalibTPCQvx[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[0]);
					fEPCalibTPCQvy[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[1]);
	
					if (!fEPCalibTPCQvx[itrig][iside][iharm])	{printf("Error : Couldn't Get Recentering Parameters(x) for harm=%d TPC\n",iharm+1);}
					if (!fEPCalibTPCQvy[itrig][iside][iharm])	{printf("Error : Couldn't Get Recentering Parameters(y) for harm=%d TPC\n",iharm+1);}
	
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
  
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs ; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
		
					fEPCalibFMDQvx[itrig][iside][iharm]	= 0x0;	fEPCalibFMDQvy[itrig][iside][iharm]	= 0x0;
	
					sprintf(h_name[0], "fPrfFMDQvcos%d%d%d", itrig+1, iside+1, iharm+1);
					sprintf(h_name[1], "fPrfFMDQvsin%d%d%d", itrig+1, iside+1, iharm+1);
	
					fEPCalibFMDQvx[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[0]);
					fEPCalibFMDQvy[itrig][iside][iharm]	= (TProfile*)fHListEPCalib->FindObject(h_name[1]);
	
					if (!fEPCalibFMDQvx[itrig][iside][iharm])	{printf("Error : Couldn't Get Recentering Parameters(x) for harm=%d FMD\n",iharm+1);}
					if (!fEPCalibFMDQvy[itrig][iside][iharm])	{printf("Error : Couldn't Get Recentering Parameters(y) for harm=%d FMD\n",iharm+1);}
	
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
  
		cout << "Calib Param Loading Finish !" << endl;
  
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fVnCent; icent++) {
	
						//Initialize recentering param
						MeanQxV0	[itrig][iside][iharm][icent] = 0.0;		MeanQyV0	[itrig][iside][iharm][icent] = 0.0;
						RMSQxV0		[itrig][iside][iharm][icent] = 0.0;		RMSQyV0		[itrig][iside][iharm][icent] = 0.0;
						
						// Load Recentering param
						MeanQxV0	[itrig][iside][iharm][icent] = fEPCalibV0Qvx	[itrig][iside][iharm]->GetBinContent(icent+1)	;
						MeanQyV0	[itrig][iside][iharm][icent] = fEPCalibV0Qvy	[itrig][iside][iharm]->GetBinContent(icent+1)	;
						RMSQxV0		[itrig][iside][iharm][icent] = fEPCalibV0Qvx	[itrig][iside][iharm]->GetBinError(icent+1)		;
						RMSQyV0		[itrig][iside][iharm][icent] = fEPCalibV0Qvy	[itrig][iside][iharm]->GetBinError(icent+1)		;

					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
  
		cout << "Load V0Recentering Param ! " << endl;
  
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fVnCent; icent++) {
	
						//Initialize recentering param
						MeanQxTPC	[itrig][iside][iharm][icent] = 0.0;	MeanQyTPC	[itrig][iside][iharm][icent] = 0.0;
						RMSQxTPC	[itrig][iside][iharm][icent] = 0.0;	RMSQyTPC	[itrig][iside][iharm][icent] = 0.0;
						
						// Load Recentering param
						MeanQxTPC	[itrig][iside][iharm][icent] = fEPCalibTPCQvx	[itrig][iside][iharm]->GetBinContent(icent+1);
						MeanQyTPC	[itrig][iside][iharm][icent] = fEPCalibTPCQvy	[itrig][iside][iharm]->GetBinContent(icent+1);
						RMSQxTPC	[itrig][iside][iharm][icent] = fEPCalibTPCQvx	[itrig][iside][iharm]->GetBinError(icent+1)	;
						RMSQyTPC	[itrig][iside][iharm][icent] = fEPCalibTPCQvy	[itrig][iside][iharm]->GetBinError(icent+1)	;
	
					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
  
		for (int itrig = 0; itrig < fTrigH; itrig++) {
			for (int iside = 0; iside < fFMDs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					for (int icent = 0; icent < fVnCent; icent++) {
		
						//Initialize recentering param
						MeanQxFMD	[itrig][iside][iharm][icent] = 0.0;	MeanQyFMD	[itrig][iside][iharm][icent] = 0.0;
						RMSQxFMD	[itrig][iside][iharm][icent] = 0.0;	RMSQyFMD	[itrig][iside][iharm][icent] = 0.0;
																									
						// Load Recentering param             
						MeanQxFMD	[itrig][iside][iharm][icent] = fEPCalibFMDQvx	[itrig][iside][iharm]->GetBinContent(icent+1);
						MeanQyFMD	[itrig][iside][iharm][icent] = fEPCalibFMDQvy	[itrig][iside][iharm]->GetBinContent(icent+1);
						RMSQxFMD	[itrig][iside][iharm][icent] = fEPCalibFMDQvx	[itrig][iside][iharm]->GetBinError(icent+1)	;
						RMSQyFMD	[itrig][iside][iharm][icent] = fEPCalibFMDQvy	[itrig][iside][iharm]->GetBinError(icent+1)	;

					}	//icent loop
				}	//iharm loop
			}	//iside loop
		}	//itrig loop
		cout << "Load TPC and FMD Recentering Param ! " << endl;
		cout << "Calib Param Loading Finish !" << endl;
  
	}	//EPCalib Load
  
	//Minimum Bias Trigger
	Bool_t itrig_kMB = (((AliInputEventHandler*)
		(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
		& AliVEvent::kMB);
  
	Bool_t itrig_kCentral = (((AliInputEventHandler*)
		(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
		& AliVEvent::kCentral);
  
	Bool_t itrig_kSemiCentral = (((AliInputEventHandler*)
		(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()
		& AliVEvent::kSemiCentral);
  
	Bool_t itrig_Accept = itrig_kMB || itrig_kCentral || itrig_kSemiCentral;
  
	if(!(itrig_Accept)){
		//if (fDebug > 1 ) Printf(" Trigger Selection: event REJECTED ... ");
		cout<<" Trigger Selection: event REJECTED ... "<<endl;
		PairCalc_Final();
		return;
	}
  
	Int_t fTrig = trig_dec(itrig_kMB, itrig_kSemiCentral, itrig_kCentral);
	if (fTrig < 0) {
		PairCalc_Final();
		return;
	}

	if (!fTrigSel) {
		if (fTrig>0) {
			PairCalc_Final();
			return;
		}
	}
	else if (fTrigSel==1) {
		if (fTrig!=1) {
			PairCalc_Final();
			return;
		}
	}
	else if (fTrigSel==2) {
		if (fTrig!=2) {
			PairCalc_Final();
			return;
		}
	}

	fH1Event -> Fill(1);

	AliAODInputHandler *aodH = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
	if(!aodH) {
		PairCalc_Final();
	}

	AliAODpidUtil * pidres = aodH->GetAODpidUtil();
	pidres -> SetUseTPCMultiplicityCorrection(kTRUE);
	pidres -> SetOADBPath("$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponse_special.root");
	pidres -> SetCustomTPCetaMaps("$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_special.root");
	pidres -> SetUseTPCEtaCorrection(kTRUE);

	AliTOFPIDResponse *fTOFResponse = 0x0;
	AliTPCPIDResponse *fTPCResponse = 0x0;
  
	if (pidres) {
		fTOFResponse = &(pidres->GetTOFResponse());
		fTPCResponse = &(pidres->GetTPCResponse());
	}
	fH1Event -> Fill(2);
  
  
	//Vertexing
	AliAODVertex *PrimaryVertex_AOD = fAOD->GetPrimaryVertex();
	Double_t Vertex[3];
	Vertex[0]=PrimaryVertex_AOD->GetX();
	Vertex[1]=PrimaryVertex_AOD->GetY();
	Vertex[2]=PrimaryVertex_AOD->GetZ();
	if(Vertex[0]<10e-5 && Vertex[1]<10e-5 &&  Vertex[2]<10e-5) return;
	if(fabs(Vertex[2]) > 8.0) {
		cout << "Vertex Cut---- " << endl;
		PairCalc_Final();
		return; // Z-Vertex Cut
	}
	fH1Zvtx[fTrig]	-> Fill(Vertex[2]);
	fH1Zvtx[3]			-> Fill(Vertex[2]);
	fH1Event -> Fill(3);

  Int_t fZvtx = (Int_t) ( ( Vertex[2] + 8.0 ) / 2.0 );
	if (fZvtx < 0 || fZvtx > 7) {
		cout << "Z vertex class has unexpected value +++++++ " << endl;
		PairCalc_Final();
		return;
	}
	fH1Event -> Fill(4);

	Double_t fMagField = fAOD->GetMagneticField();
	fMagsign = MagDec(fMagField);
  
	//Centrality AOD
	AliAODHeader *aodheader = (AliAODHeader*) fAOD->GetHeader();
	cent_V0M		= aodheader->GetCentrality();
	cent_Calib	= aodheader->GetCentrality();
	Int_t fCent				= (Int_t)(cent_V0M / 5.0)	;
	Int_t fCentCalib	= (Int_t)(cent_V0M / 5.0)	;
	Int_t fCent10			= (Int_t)(cent_V0M / 10.0);
	Int_t fCent_H			= HistCent(cent_V0M)			;

	if (cent_V0M < 0.0 || cent_V0M >100.0) {
		PairCalc_Final();
		return;
	}
	if (fCent < 0 || fCent > 18) {
		PairCalc_Final();
		return;
	}
	if (fCentCalib < 0 || fCentCalib > 18) {
		PairCalc_Final();
		return;
	}
  
	if (fTrig==2 && fCentCalib == 2)	{
		cent_Calib = 7.5	;
		fCentCalib = 1		;
		fCent10 = 0				;
	}
	if (fTrig==1 && fCent < 2)	{
		cent_Calib	= 12.5	;
		fCentCalib = 2			;
		fCent10 = 1					;
	}
	if (fTrig==1 && fCent > 9){
		cent_Calib = 45.5	;
		fCentCalib = 9		;
		fCent10	= 4				;
	}
	if (fTrig==0 && fCent == 18) {
		cent_V0M = 87.5;
		fCent		= 17;
		fCent10 = 8	;
	}
	fH1Event -> Fill(5);
  
	////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////
	for (int icor = 0; icor < fCorr; icor++) {
		for (int iside = 0; iside < fV0s; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				QxV0	[icor][iside][iharm] = 0.0;	QyV0	[icor][iside][iharm] = 0.0;
			}	//iharm loop
		}	//iside loop
	}	//icor loop
  
	Double_t	V0MC, V0MA, V0MCA;
	//---------------------------------------	
	V0MC	= 0.0;	V0MA	= 0.0;	V0MCA = 0.0;
  
	//******************************************************************************************************
	//VZERO EP Calibration
	//******************************************************************************************************
  
	AliAODVZERO *aodV0 = fAOD -> GetVZEROData();
  
	if (!aodV0) {
		cout << "ERROR: aodV0 not available" << endl;
		PairCalc_Final();
		return;
	}
  
	for (int iseg = 0; iseg < 64; iseg++) {
  
		Float_t v0_amp	= aodV0 -> GetMultiplicity(iseg)	;
  
		fPrfV0Sig		-> Fill(iseg, v0_amp)		;
		fH2V0Sig		-> Fill(iseg, v0_amp)		;
  
		Float_t v0phi = 22.5 + 45.0 * (iseg%8);
		v0phi = v0phi * PI / 180.0;
  
		if (iseg<32) V0MC += v0_amp;
		else				 V0MA += v0_amp;
		V0MCA += v0_amp;
  
		//[correction][side][harmonics]
		for (int iharm = 0; iharm < fHarm; iharm++) {
			if (iseg<32) {	//C side
				QxV0[0][0][iharm] += (v0_amp)*TMath::Cos((iharm+1)*v0phi)	;
				QyV0[0][0][iharm] += (v0_amp)*TMath::Sin((iharm+1)*v0phi)	;
			}
			else { //A side
				QxV0[0][1][iharm] += (v0_amp)*TMath::Cos((iharm+1)*v0phi)	;
				QyV0[0][1][iharm] += (v0_amp)*TMath::Sin((iharm+1)*v0phi)	;
			}
			//A + C side
			QxV0[0][2][iharm] += (v0_amp)*TMath::Cos((iharm+1)*v0phi)		;
			QyV0[0][2][iharm] += (v0_amp)*TMath::Sin((iharm+1)*v0phi)		;
		}
	}
  
	fV0flag	= 0;
	//V0Multiplicity Normalization
	for (int iharm = 0; iharm < fHarm; iharm++) {
		if ( V0MC > 0.0 && V0MA > 0.0 ) {
			QxV0[0][0][iharm] = QxV0[0][0][iharm] / sqrt( V0MC	);		QyV0[0][0][iharm] = QyV0[0][0][iharm] / sqrt( V0MC	);
			QxV0[0][1][iharm] = QxV0[0][1][iharm] / sqrt( V0MA	);		QyV0[0][1][iharm] = QyV0[0][1][iharm] / sqrt( V0MA	);
			QxV0[0][2][iharm] = QxV0[0][2][iharm] / sqrt( V0MCA	);		QyV0[0][2][iharm] = QyV0[0][2][iharm] / sqrt( V0MCA );
		}
		else {
			fV0flag++;
			cout << "**********************************"	<< endl;
			cout << " ERROR: VZERO multiplicity is 0 "		<< endl;
			cout << "**********************************"	<< endl;
		}
	}
  
	for (int iside = 0; iside < fV0s; iside++) {
		for (int iharm = 0; iharm < fHarm; iharm++) {
			if (RMSQxV0 [fTrig][iside][iharm][fCentCalib]*RMSQyV0 [fTrig][iside][iharm][fCentCalib] == 0.0) {
				fV0flag++;
				continue;
			}
			// Recentering
			QxV0[1][iside][iharm] = ( QxV0[0][iside][iharm] - MeanQxV0[fTrig][iside][iharm][fCentCalib] ) / RMSQxV0 [fTrig][iside][iharm][fCentCalib];
			QyV0[1][iside][iharm] = ( QyV0[0][iside][iharm] - MeanQyV0[fTrig][iside][iharm][fCentCalib] ) / RMSQyV0 [fTrig][iside][iharm][fCentCalib];
		}	//iharm loop
	}	//iside loop
	cout << "V0 recentering finish !" << endl;
  
	if (fV0flag == 0) {
		for (int icor = 0; icor < fCorr; icor++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfV0Qvcos	[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QxV0[icor][iside][iharm]);
					fPrfV0Qvsin	[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QyV0[icor][iside][iharm]);
  
					fH1V0Qvcos	[icor][fTrig][iside][iharm][fCent10] -> Fill(QxV0[icor][iside][iharm]);
					fH1V0Qvsin	[icor][fTrig][iside][iharm][fCent10] -> Fill(QyV0[icor][iside][iharm]);
				}	//iharm loop
			}	//iside loop
		}	//icor loop
		cout << "V0 Qvector Fill !" << endl;
  
		for (int icor = 0; icor < fCorr; icor++) {
			for (int iside = 0; iside < fV0s; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fV0psi	[icor][iside][iharm]= (1.0/(iharm+1.0))*TMath::ATan2(QyV0	[icor][iside][iharm], QxV0[icor][iside][iharm]);
					fH1EPV0	[icor][fTrig][iside][iharm][fCent10] -> Fill((iharm+1.0)*fV0psi[icor][iside][iharm]);
				}	//iharm	 loop
			}	//iside loop
		}	//icor loop
	}	//V0 recentering param request
  
  
	//======================================================================================================
	//TPC EP Calibration
	//======================================================================================================
	//Initialization of TPC Qvector  & Ncharge counter
	Int_t			fNch	[7] = {}		;
	Double_t	fQwgt	[7] = {0.0}	;
	Double_t	fWgt			= 0.0		;
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int iside = 0; iside < fTPCs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				QxTPC[icor][iside][iharm] = 0.0;	QyTPC[icor][iside][iharm] = 0.0;
			}
		}
	}
  
	cout << "TPC Calibration" << endl;
	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
		AliAODTrack *ttrack = (AliAODTrack*) (fAOD->GetTracks())->At(iTracks);
		if (!ttrack->TestFilterBit(1<<(7))) continue;
		if (!ttrack) {
			printf("ERROR: Could not receive ttrack %d\n", iTracks);
			continue;
		}
		if (!ttrack->IsPrimaryCandidate()) continue;										//Kink rejection
		if (!(ttrack->GetType()== AliAODTrack::kPrimary)) continue;
		if (ttrack->GetTPCNcls() < 70) continue;												//# of TPC cluster
		if (TMath::Abs(ttrack->Eta()) > 1.0) continue;
		if (ttrack->Chi2perNDF() > 4.0) continue;
		if (ttrack->GetLabel() < 0) continue;
  
		AliAODTrack* trk_clone = (AliAODTrack*)ttrack->Clone("trk_clone");
		Double_t  impact		[2];
		Double_t  covimpact	[3];
		Double_t	fImpactXY	= -999.;
		Double_t	fImpactZ	= -999.;
		if (trk_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(), 20.,impact,covimpact)) {
			fImpactXY = impact[0];
			fImpactZ  = fabs(impact[1]) - fabs(Vertex[2]);
		}
		else {
			fImpactXY = -1000.0;
			fImpactZ = -1000.0;
		}
		delete trk_clone;

		//fImpactXY	= ttrack->DCA();
		//fImpactZ	= ttrack->ZAtDCA();
  
		if (TMath::Abs(fImpactXY)> 2.4) continue;
		if (TMath::Abs(fImpactZ) > 3.0) continue;
  
		Double_t	fPhi		= ttrack->Phi();
		Double_t	fEta		= ttrack->Eta();
		Double_t	fPt			= ttrack->Pt();
		Int_t		fnEta2	= TPCeta2(fEta);
		Int_t		fnEta 	= TPCeta(fEta);

		if (fPt < 2.0)	fWgt	= fPt	;
		else						fWgt	= 2.0	;
  
		if (fnEta<0)	continue;
		if (fnEta2<0)	continue;
  
		//fH1Phi	->	Fill(fPhi);
		//fH1Eta	->	Fill(fEta);
		//fH1Pt		->	Fill(fPt);
  
		fH1EtaD[fnEta]	->	Fill(fEta);
		fH1EtaD[fnEta2]	->	Fill(fEta);
		fH1EtaD[6]			->	Fill(fEta);
  
		fNch[fnEta]++					;
		fNch[fnEta2]++				;
		fNch[6]++							;
  
		fQwgt[fnEta]	+= fWgt	;
		fQwgt[fnEta2]	+= fWgt	;
		fQwgt[6]			+= fWgt	;
  
		for (int iharm = 0; iharm < fHarm; iharm++) {
			QxTPC[0][fnEta][iharm]	+= fWgt * TMath::Cos((iharm+1)*fPhi);
			QyTPC[0][fnEta][iharm]	+= fWgt * TMath::Sin((iharm+1)*fPhi);
  
			QxTPC[0][fnEta2][iharm]	+= fWgt * TMath::Cos((iharm+1)*fPhi);
			QyTPC[0][fnEta2][iharm]	+= fWgt * TMath::Sin((iharm+1)*fPhi);
  
			QxTPC[0][6][iharm]			+= fWgt * TMath::Cos((iharm+1)*fPhi);
			QyTPC[0][6][iharm]			+= fWgt * TMath::Sin((iharm+1)*fPhi);
		}	//iharm loop
	}	//iTracks loop
  
	fH1Nch -> Fill(fNch[6]);
  
	fTPCflag	= 0;
	//if (fQwgt[0] > 0.0 && fQwgt[1] > 0.0 && fQwgt[2] > 0.0 && fQwgt[3] > 0.0) {
	if (fNch[0] > 0 && fNch[1] > 0 && fNch[2] > 0 && fNch[3] > 0) {
		//TPC Pt weight Normalization
		for (int iside = 0; iside < fTPCs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				QxTPC[0][iside][iharm] = QxTPC[0][iside][iharm] / TMath::Sqrt(	fQwgt[iside]	);
				QyTPC[0][iside][iharm] = QyTPC[0][iside][iharm] / TMath::Sqrt(	fQwgt[iside]	);
			}
		}
		for (int iside = 0; iside < fTPCs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				if (RMSQxTPC[fTrig][iside][iharm][fCentCalib] * RMSQyTPC[fTrig][iside][iharm][fCentCalib] == 0.0) {
					fTPCflag++	;
					continue		;
				}
				//Recentering
				QxTPC[1][iside][iharm] = ( QxTPC[0][iside][iharm] - MeanQxTPC[fTrig][iside][iharm][fCentCalib] ) / RMSQxTPC[fTrig][iside][iharm][fCentCalib]	;
				QyTPC[1][iside][iharm] = ( QyTPC[0][iside][iharm] - MeanQyTPC[fTrig][iside][iharm][fCentCalib] ) / RMSQyTPC[fTrig][iside][iharm][fCentCalib]	;
			}	//iside loop
		}	//iharm loop
	}	//fTPC flag request
	else {
		fTPCflag++;
	}
	cout << "TPC Recentering Finish ! " << endl;
  
	if (fTPCflag == 0) {
		for (int icor = 0; icor < fCorr; icor++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fPrfTPCQvcos[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QxTPC[icor][iside][iharm]);
					fPrfTPCQvsin[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QyTPC[icor][iside][iharm]);
  
					fH1TPCQvcos	[icor][fTrig][iside][iharm][fCent10] -> Fill(QxTPC[icor][iside][iharm]);
					fH1TPCQvsin	[icor][fTrig][iside][iharm][fCent10] -> Fill(QyTPC[icor][iside][iharm]);
				}	//iharm loop
			}	//iside loop
		}	//icor loop
  
		for (int icor = 0; icor < fCorr; icor++) {
			for (int iside = 0; iside < fTPCs; iside++) {
				for (int iharm = 0; iharm < fHarm; iharm++) {
					fTPCpsi	[icor][iside][iharm]  = (1.0/(iharm+1))*TMath::ATan2(QyTPC[icor][iside][iharm], QxTPC[icor][iside][iharm]);
					fH1EPTPC[icor][fTrig][iside][iharm][fCent10] -> Fill((iharm+1)*fTPCpsi[icor][iside][iharm]);
				}	//iharm	 loop
			}	//iside loop
		}	//icor loop
	}	//TPC flag
  
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
	//FMD EP Calibration
	//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  
	AliAODEvent *fAOD_FMD = AliForwardUtil::GetAODEvent(this);
  
	AliAODForwardMult* aodfmult =
		static_cast<AliAODForwardMult*>(fAOD_FMD->FindListObject("Forward"));
  
	if (!aodfmult) {
		cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"	<< endl;
		cout << "ForwardMult is not available !!!!!!!"	<< endl;
		cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"	<< endl;
		PairCalc_Final();
		return;
	}
  
	if (!aodfmult->IsTriggerBits(AliAODForwardMult::kOffline)) {
		cout << "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
		PairCalc_Final();
		return;
	}
  
	for (int icor = 0; icor < fCorr; icor++) {
		for (int iside = 0; iside < fFMDs; iside++) {
			for (int iharm = 0; iharm < fHarm; iharm++) {
				QxFMD[icor][iside][iharm] = 0.0;	QyFMD[icor][iside][iharm] = 0.0;
			}
		}
	}
  
	fFMD = &aodfmult->GetHistogram();
	fH2FMD2D -> Add(fFMD);
  
	Double_t Sval = 0.0;
	TH1F	*fFMDphi;
  
  fFMDflag = 0;
	for (int iharm = 0; iharm < fHarm; iharm++) {
		for (int iside = 0; iside < fFMDs; iside++) {
			fFMDphi = ProjectPhiFMD(fFMD, EtaLimitFMDAC1(iside),EtaLimitFMDAC2(iside));
			Sval = 0.;
			for(int i=1;i<=fFMDphi->GetNbinsX();i++){
				Double_t val=fFMDphi->GetBinContent(i)							;
				Double_t phi=fFMDphi->GetBinCenter (i)							;
				QxFMD[0][iside][iharm] += val*cos((iharm+1.0)*phi)	;
				QyFMD[0][iside][iharm] += val*sin((iharm+1.0)*phi)	;
				Sval += val																					;
			}
			if (Sval) {
				QxFMD[0][iside][iharm] /= TMath::Sqrt(Sval);
				QyFMD[0][iside][iharm] /= TMath::Sqrt(Sval);
  
				if (RMSQxFMD[fTrig][iside][iharm][fCentCalib] * RMSQyFMD[fTrig][iside][iharm][fCentCalib] == 0.0) {
					fFMDflag++;
					cout << "FMD Calibration parameter is ++++ ZERO +++" << endl;
					continue	;
				}
				//Recentering
				QxFMD[1][iside][iharm] = ( QxFMD[0][iside][iharm] - MeanQxFMD[fTrig][iside][iharm][fCentCalib] ) / RMSQxFMD[fTrig][iside][iharm][fCentCalib]	;
				QyFMD[1][iside][iharm] = ( QyFMD[0][iside][iharm] - MeanQyFMD[fTrig][iside][iharm][fCentCalib] ) / RMSQyFMD[fTrig][iside][iharm][fCentCalib]	;
  
				for (int icor = 0; icor < fCorr; icor++) {
					fFMDpsi	[icor][iside][iharm]  = (1.0/(iharm+1))*TMath::ATan2(QyFMD[icor][iside][iharm], QxFMD[icor][iside][iharm]);
					fH1EPFMD[icor][fTrig][iside][iharm][fCent10] -> Fill((iharm+1)*fFMDpsi[icor][iside][iharm]);
  
					fPrfFMDQvcos[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QxFMD[icor][iside][iharm]);
					fPrfFMDQvsin[icor][fTrig][iside][iharm] -> Fill(cent_Calib, QyFMD[icor][iside][iharm]);
  
					fH1FMDQvcos	[icor][fTrig][iside][iharm][fCent10] -> Fill(QxFMD[icor][iside][iharm]);
					fH1FMDQvsin	[icor][fTrig][iside][iharm][fCent10] -> Fill(QyFMD[icor][iside][iharm]);
				}	//icor loop
			}	//Sval
			else {
				cout << "FMD Multiplicity is ++++ ZERO +++" << endl;
				fFMDflag++;
			}	//Sval = 0
		}	//iside loop
	}	//iharm loop
  
	cout << "FMD Flattening Param Calc ! " << endl;
  
	//##############################################################
	//##############################################################
  
	Float_t fEP3FMD	= fFMDpsi[1][12][2];
	Int_t		fnEP3FMD	= (Int_t)( ( fFMDpsi[1][12][2] + ( TMath::Pi() ) / 3.0 ) / ( TMath::Pi() / 30.0) );
	//fEP3V0	= fV0psi[2][2][2];
	//fnEP3V0	= (Int_t)( ( fV0psi[2][2][2] + ( TMath::Pi() ) / 3.0 ) / ( TMath::Pi() / 24.0) );
	//if (fnEP3V0 < 0) {
	//	PairCalc_Final();
	//	return;
	//}
	//##############################################################
	//##############################################################
  
	FillEPCorr(fV0flag, fTPCflag, fFMDflag);

	if (fnEP3FMD < 0) {
		PairCalc_Final();
		return;
	}

	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
		AliAODTrack *fvntrack = (AliAODTrack*) (fAOD->GetTracks())->At(iTracks);
		if (!fvntrack->TestFilterBit(1<<(7))) continue;
		if (!fvntrack) {
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
		if (!fvntrack->IsPrimaryCandidate()) continue;										//Kink rejection
		if (!(fvntrack->GetType()== AliAODTrack::kPrimary)) continue;
		if (fvntrack->GetTPCNcls() < 70) continue;												//# of TPC cluster
		if (TMath::Abs(fvntrack->Eta()) > 0.8) continue;
		if (fvntrack->Chi2perNDF() > 4.0) continue;
		if (fvntrack->GetLabel() < 0) continue;


		AliAODTrack* trk_clone = (AliAODTrack*)fvntrack->Clone("trk_clone");
		Double_t  impact		[2];
		Double_t  covimpact	[3];
		Double_t	fImpactXY	= -999.;
		Double_t	fImpactZ	= -999.;
		if (trk_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),20.,impact,covimpact)) {
			fImpactXY = impact[0];
			fImpactZ  = fabs(impact[1]) - fabs(Vertex[2]);
		}

		if (TMath::Abs(fImpactXY)> 2.4) continue;
		if (TMath::Abs(fImpactZ) > 3.0) continue;

		delete trk_clone;

		Double_t	fPhi	= fvntrack->Phi()	;
		Double_t	fPt		=	fvntrack->Pt()	;
		Double_t	fEta	= fvntrack->Eta()	;
		Int_t			fnEta = TPCeta(fEta)		;
  
		Double_t fvncos	= 0.0;
		Double_t fvnsin	= 0.0;

		for (int iharm = 0; iharm < fHarm; iharm++) {
			for (int iside = 0; iside < fV0s; iside++) {

				if (fV0flag == 0) {
					fvncos	= TMath::Cos( (iharm+1.0) * (fPhi - fV0psi	[1][iside][iharm] ) );
					fvnsin	= TMath::Sin( (iharm+1.0) * (fPhi - fV0psi	[1][iside][iharm] ) );

					fPrfvncos_pt		[1][iside][iharm][fCent]	-> Fill(fPt,			fvncos);
					fPrfvnsin_pt		[1][iside][iharm][fCent]	-> Fill(fPt,			fvnsin);
					if ( ( fPt > 0.2 ) && ( fPt < 5.0 ) ) {
						fPrfvncos_cent	[1][iside][iharm]					-> Fill(cent_V0M,	fvncos);
						fPrfvnsin_cent	[1][iside][iharm]					-> Fill(cent_V0M,	fvnsin);
					}
				}
				if (fFMDflag == 0) {

					Int_t fConv_FMD[3] = {3, 9, 12};
					fvncos	= TMath::Cos( (iharm+1.0) * (fPhi - fFMDpsi	[1][fConv_FMD[iside]][iharm] ) );
					fvnsin	= TMath::Sin( (iharm+1.0) * (fPhi - fFMDpsi	[1][fConv_FMD[iside]][iharm] ) );

					fPrfvncos_pt		[2][iside][iharm][fCent]	-> Fill(fPt,			fvncos);
					fPrfvnsin_pt		[2][iside][iharm][fCent]	-> Fill(fPt,			fvnsin);
					if ( ( fPt > 0.2 ) && ( fPt < 5.0 ) ) {
						fPrfvncos_cent	[2][iside][iharm]					-> Fill(cent_V0M,	fvncos);
						fPrfvnsin_cent	[2][iside][iharm]					-> Fill(cent_V0M,	fvnsin);
					}
				}	
			}	//iside loop

			if (fTPCflag == 0) {
				if ((fnEta == 0) || (fnEta == 3)) {

					Int_t fConv_TPCE[4] = {3, -1, -1, 0};
					Int_t fConv_TPCV[4] = {0, -1, -1, 1};

					if (fConv_TPCE[fnEta] < 0) continue;
					if (fConv_TPCV[fnEta] < 0) continue;

					fvncos	= TMath::Cos( (iharm+1.0) * (fPhi - fTPCpsi	[1][fConv_TPCE[fnEta]][iharm] ) );
					fvnsin	= TMath::Sin( (iharm+1.0) * (fPhi - fTPCpsi	[1][fConv_TPCE[fnEta]][iharm] ) );

					fPrfvncos_pt		[0][fConv_TPCV[fnEta]][iharm][fCent]	-> Fill(fPt,			fvncos);
					fPrfvnsin_pt		[0][fConv_TPCV[fnEta]][iharm][fCent]	-> Fill(fPt,			fvnsin);
					if ( ( fPt > 0.2 ) && ( fPt < 5.0 ) ) {
						fPrfvncos_cent	[0][fConv_TPCV[fnEta]][iharm]				-> Fill(cent_V0M,	fvncos);
						fPrfvnsin_cent	[0][fConv_TPCV[fnEta]][iharm]				-> Fill(cent_V0M,	fvnsin);
					}
				}
			}
		}	//iharm loop
	}	//Track loop

	if (fCent > 9) {
		PairCalc_Final();
		return;
	}

	if (fCent_H < 0 || fCent_H > 5) {
		PairCalc_Final();
		return;
	}

	if (fFMDflag > 0) {
		PairCalc_Final();
		return;
	}

	if (!fCentralAna) {
		if (fCent < 2) {
			PairCalc_Final();
			return;
		}
	}
	else {
		if (fCent > 1) {
			PairCalc_Final();
			return;
		}
	}

	if (!fCentralAna) fCent_H = fCent_H-2	;

	fH1Cent[fTrig]	-> Fill(cent_V0M);
	fH1Cent[3]			-> Fill(cent_V0M);
  
	Int_t Nch_n1 = 0;	Int_t Nch_n2 = 0;
  
	Int_t labels[25000];
	for (Int_t il=0; il<25000; il++) labels[il] = -1;
  
	Int_t Nch_label;
  
	Nch_label = 0;
  
	cout << "CentMix : " << fCent << ", ZVtx : " << fZvtx << ", Psi3 : " << fnEP3FMD << ", Mixflag : " << mix_flag[0][fCent][fZvtx][fnEP3FMD] << endl;
  
	// looking for global tracks and saving their numbers to copy from them PID information to TPC-only tracks in the main loop over tracks
	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
		const AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(iTracks));
		if (!aodtrack->TestFilterBit(1<<(7))) {
			if(aodtrack->GetID() < 0) {
				continue;
			}
			else {
				labels[aodtrack->GetID()] = iTracks;
				Nch_label++;
			}
		}
	}
  
	for (Int_t iTracks = 0; iTracks < fAOD->GetNumberOfTracks(); iTracks++) {
		Int_t fPID = 0;
		AliAODTrack *track = (AliAODTrack*) (fAOD->GetTracks())->At(iTracks);
  
		if (!track->TestFilterBit(1<<(7))) continue;
		if (!track) {
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
  
		AliAODTrack* trk_clone = (AliAODTrack*)track->Clone("trk_clone");
		Double_t  impact		[2];
		Double_t  covimpact	[3];
		Double_t	fImpactXY	= -999.;
		Double_t	fImpactZ	= -999.;

		if (trk_clone->PropagateToDCA(fAOD->GetPrimaryVertex(),fAOD->GetMagneticField(),20.,impact,covimpact)) {
			fImpactXY = impact[0];
			fImpactZ  = fabs(impact[1]) - fabs(Vertex[2]);
		}
  
		if (!(track->GetType()== AliAODTrack::kPrimary)) continue	;
		if (track->GetTPCNcls() < 80) continue										;
		if (TMath::Abs(track->Eta()) > 0.8) continue							;
		if (track->Chi2perNDF() > 4.0) continue										;
		if (track->GetLabel() < 0) continue												;

		fH2DCA		-> Fill(fImpactXY, fImpactZ);
		fH2DCAxyPt-> Fill(track->Pt(), fImpactXY);
		fH2DCAzPt -> Fill(track->Pt(), fImpactZ);

		if (TMath::Abs(fImpactXY) > 2.4) continue									;
		if (TMath::Abs(fImpactZ) > 3.0) continue									;

		delete trk_clone;
  
		AliAODTrack *aodtrack = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(labels[-1-fAOD->GetTrack(iTracks)->GetID()]));
  
		if (!aodtrack) {
			printf("ERROR: Could not receive track %d\n", iTracks);
			continue;
		}
  
		Int_t fCharge = track->Charge();
		Double_t	fP	= track->P()	;
		Double_t	fPx	= track->Px()	;
		Double_t	fPy	= track->Py()	;
		Double_t	fPz	= track->Pz()	;
		Double_t	fPt = track->Pt()	;
  
		Double_t	fE_Pi	= TMath::Sqrt(fP*fP+M2Pion);
		Double_t	fPhi  = track->Phi();
		Double_t	fEta	= track->Eta();

		fH1Phi	->	Fill(fPhi);
		fH1Eta	->	Fill(fEta);
		fH1Pt		->	Fill(fPt);
  
  
		fH2PhiEta -> Fill(fPhi, fEta);
  
		/************************************************/
		//PID parameters
		/************************************************/
		Double_t	TOFtof		= aodtrack->GetTOFsignal();
		Double_t	TPCeloss	= aodtrack->GetTPCsignal();
		Double_t	fTPCP			= aodtrack->GetTPCmomentum();

		Float_t		TOFexp[3], TOFnsig[3];
		Float_t		TPCexp[3], TPCnsig[3];
		Double_t	TOFsig[3];
		Double_t	TPCsig[3];

		Double_t times[AliPID::kSPECIES];
		/************************************************/
    
		if (aodtrack->GetStatus() & AliESDtrack::kTOFpid) {
			Float_t TimeZeroTOF = fTOFResponse->GetStartTime(fP);
			aodtrack -> GetIntegratedTimes(times);
		}
    
		float vp = -1000.;
		for (int ipart = 0; ipart < 3; ipart++) {
			//PID parameters
			TOFexp[ipart] = -100000.;
			TOFsig[ipart] = -1000.;
			TOFnsig[ipart] = -1000.;
    
			if ((aodtrack->GetStatus() & AliESDtrack::kTOFpid) &&
					(aodtrack->GetStatus() & AliESDtrack::kTOFout) &&
					(aodtrack->GetStatus() & AliESDtrack::kTIME) ) {
				if (aodtrack->IsOn(AliESDtrack::kTOFpid)) {

					//Double_t len = 200; // esdtrack->GetIntegratedLength(); !!!!!
					Double_t len = aodtrack->GetIntegratedLength(); // esdtrack->GetIntegratedLength(); !!!!!
					Double_t tof = aodtrack->GetTOFsignal();
					if (tof > 0.) vp = len / tof / 0.03;
					fH2TOF -> Fill(fP, vp);

					TOFexp[ipart] = fTOFResponse->GetExpectedSignal(aodtrack, AliPID::EParticleType(ipart+2));
					TOFsig[ipart] = fTOFResponse->GetExpectedSigma(fTPCP,times[ipart+2],AliPID::ParticleMassZ(ipart+2));
					TOFnsig[ipart] = pidres->NumberOfSigmasTOF(aodtrack, AliPID::EParticleType(ipart+2));

					fH2NsigTOF[ipart] -> Fill(fTPCP, TOFnsig[ipart]);

				}
			}
		}
    
		for (int ipart = 0; ipart < 3; ipart++) {
			//PID parameters
			TPCexp	[ipart] = -100000.;
			TPCsig	[ipart] = -1000.;
			TPCnsig	[ipart] = -1000.;
    
			if (aodtrack->IsOn(AliESDtrack::kTPCpid)) {
				TPCexp		[ipart] = fTPCResponse->GetExpectedSignal(fTPCP, AliPID::EParticleType(ipart+2));
				TPCsig		[ipart] = fTPCResponse->GetExpectedSigma(fTPCP, aodtrack->GetTPCNcls(), AliPID::EParticleType(ipart+2));
				TPCnsig		[ipart] = fTPCResponse->GetNumberOfSigmas(fTPCP, TPCeloss, aodtrack->GetTPCNcls(), AliPID::EParticleType(ipart+2));
				fH2NsigTPC[ipart] -> Fill(fTPCP, TPCnsig[ipart]);
			}
		}
  
		if			(IsPionNSigma(fP, TPCnsig[0], TOFnsig[0]) == true)		fPID = fCharge*1;
		else if (IsKaonNSigma(fP, TPCnsig[1], TOFnsig[1]) == true)		fPID = fCharge*2;
		else if (IsProtonNSigma(fP, TPCnsig[2], TOFnsig[2]) == true)	fPID = fCharge*3;

		//pT cut for HBT
		if (track->Pt() > 1.5)	continue;
		if (track->Pt() < 0.15) continue;
  
  
		if (fPID == 1) {
			pion_ppx		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fPx)		;
			pion_ppy		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fPy)		;
			pion_ppz		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fPz)		;
			pion_ppt		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fPt)		;
			pion_pe			[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fE_Pi)	;
			pion_pphi		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fPhi)		;
			pion_peta		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fEta)		;
			pion_pclu		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(track -> GetTPCClusterMap())	;
			pion_psha		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(track -> GetTPCSharedMap())		;
			pion_pax		[fCent][fZvtx][fnEP3FMD][mix_flag[0][fCent][fZvtx][fnEP3FMD]].push_back(fEP3FMD)	;
			Nch_n1++;
		}
		if (fPID == -1) {
			pion_mpx		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fPx)		;
			pion_mpy		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fPy)		;
			pion_mpz		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fPz)		;
			pion_mpt		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fPt)		;
			pion_me			[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fE_Pi)	;
			pion_mphi		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fPhi)		;
			pion_meta		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fEta)		;
			pion_mclu		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(track -> GetTPCClusterMap())	;
			pion_msha		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(track -> GetTPCSharedMap())		;
			pion_max		[fCent][fZvtx][fnEP3FMD][mix_flag[1][fCent][fZvtx][fnEP3FMD]].push_back(fEP3FMD)	;
			Nch_n2++;
		}
	}
  
	if (Nch_n1!=0) {
		mix_flag[0][fCent][fZvtx][fnEP3FMD]++;
	}
	if (Nch_n2!=0) {
		mix_flag[1][fCent][fZvtx][fnEP3FMD]++;
	}
  
	if (mix_flag[0][fCent][fZvtx][fnEP3FMD] == fMix) {
		cout << "*****************************" << endl;
		cout << "* ++ Real Pair Calculation  *" << endl;
		cout << "*****************************" << endl;
		for (int ieve = 0; ieve < fMix; ieve++) {
			for (ULong_t i_loop = 0; i_loop < pion_ppx [fCent][fZvtx][fnEP3FMD][ieve].size()-1; i_loop++) {
				for (ULong_t j_loop = i_loop+1; j_loop < pion_ppx [fCent][fZvtx][fnEP3FMD][ieve].size(); j_loop++) {
  
					Double_t pxi		= pion_ppx	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pxj		= pion_ppx	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pyi		= pion_ppy	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pyj		= pion_ppy	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pzi		= pion_ppz	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pzj		= pion_ppz	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pei		= pion_pe		[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pej		= pion_pe		[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pti		= pion_ppt	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t ptj		= pion_ppt	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t phii		= pion_pphi	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t phij		= pion_pphi	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t etai		= pion_peta	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t etaj		= pion_peta	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					TBits		 clu1		= pion_pclu	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		 clu2		= pion_pclu	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					TBits		 sha1		= pion_psha	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		 sha2		= pion_psha	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t EPi		= pion_pax	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);
  
					Double_t dphi	= phii - phij;
					Double_t deta	= etai - etaj;
  
					Double_t kx = 0.5*(pxi + pxj);
					Double_t ky = 0.5*(pyi + pyj);
					Double_t kt = sqrt(kx*kx + ky*ky);
					fH1kt_all -> Fill(kt);
  
					Int_t	fkT = (int)(kt*10.0) - 2;
					if (fkT < 0 || fkT > 17) continue;

					fH1kt	-> Fill(kt);
  
					Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
					Double_t dphi_E	= dphi_p - EPi;
  
					dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

					Int_t	fDphi_E	= PiEPAngle(dphi_E);
					if (fDphi_E < 0) continue;

					Int_t	fOrder;
					if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;
  
					//Calculate Phi* ( distance = r[m] )
					Double_t dphi_s[2];
					for (int idist = 0; idist < 2; idist++) {
						Double_t dist = dec_dist(idist);
						if (dist < 0) continue;
						dphi_s[idist]	= dphi + TMath::ASin(-0.075*fMagsign*dist/pti) - TMath::ASin(-0.075*fMagsign*dist/ptj);
						dphi_s[idist] = fOrder*dphi_s[idist];
						dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
					}
					dphi = fOrder*dphi;
  
					//3 sigma cut
					if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;
  
					Float_t	sharity;
					if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
					else						sharity = Sharity(clu2, clu1, sha2, sha1);
  
					if (sharity >= 0.05) continue;
  
					Double_t qinv, c_weight;
					Double_t qside, qout, qlong;
					if (fOrder > 0) {
						qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
						qout 		= qout_cms(pxi, pxj, pyi, pyj);
						qside 	= qside_cms(pxi, pxj, pyi, pyj);
						qlong 	= qlong_cms(pzi, pzj, pei, pej);
					}
					else {
						qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
						qout 		= qout_cms(pxj, pxi, pyj, pyi);
						qside 	= qside_cms(pxj, pxi, pyj, pyi);
						qlong 	= qlong_cms(pzj, pzi, pej, pei);
					}
  
					if (fabs(qinv) <= 0.34) {
						if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
						else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
						for (int iqbin = 0; iqbin < fQbin; iqbin++) {
							fH1Qinv	[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
							fH1CQinv[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
						}
					}
  
					if (conv_over(qout, qside, qlong) == 1) continue;
  
					for (int iqbin = 0; iqbin < fQbin; iqbin++) {
						fH3Q[0][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
					}
  
					Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
					Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
					Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
					if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
						fPrfConv	[0][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql, qinv);
						//fPrfConv	[0][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2 + 400*conv_ql2, qinv);
					}
					if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
						fPrfConv	[0][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3 + 3600*conv_ql3, qinv);
					}
				}	//j_loop
			}	//i_loop
		}	//ieve
  
		cout << "***************************" << endl;
		cout << "* ++ Mix Pair Calculation *" << endl;
		cout << "***************************" << endl;
  
		for (int ieve = 0; ieve < fMix-1; ieve++) {
			for (int jeve = ieve+1; jeve < fMix; jeve++) {
				for (ULong_t i_loop = 0; i_loop < pion_ppx [fCent][fZvtx][fnEP3FMD][ieve].size(); i_loop++) {
					for (ULong_t j_loop = 0; j_loop < pion_ppx [fCent][fZvtx][fnEP3FMD][jeve].size(); j_loop++) {
  
						Double_t pxi		= pion_ppx	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pxj		= pion_ppx	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pyi		= pion_ppy	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pyj		= pion_ppy	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pzi		= pion_ppz	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pzj		= pion_ppz	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pei		= pion_pe		[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pej		= pion_pe		[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pti		= pion_ppt	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t ptj		= pion_ppt	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t phii		= pion_pphi	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t phij		= pion_pphi	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t etai		= pion_peta	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t etaj		= pion_peta	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						TBits		clu1		= pion_pclu	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		clu2		= pion_pclu	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						TBits		sha1		= pion_psha	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		sha2		= pion_psha	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t EPi		= pion_pax	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);
  
						Double_t dphi	= phii - phij;
						Double_t deta	= etai - etaj;
  
						Double_t kx = 0.5*(pxi + pxj);
						Double_t ky = 0.5*(pyi + pyj);
						Double_t kt = sqrt(kx*kx + ky*ky);
		
						Int_t	fkT = (int)(kt*10.0) - 2;
						if (fkT < 0 || fkT > 17) continue;
  
						Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
						Double_t dphi_E	= dphi_p - EPi;
		
						dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

						Int_t	fDphi_E	= PiEPAngle(dphi_E);
						if (fDphi_E < 0) continue;
  
						Int_t	fOrder;
						if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;
  
						//Calculate Phi* ( distance = r[m] )
						Double_t dphi_s[2];
						for (int idist = 0; idist < 2; idist++) {
							Double_t dist = dec_dist(idist);
							if (dist < 0) continue;
							dphi_s[idist]	= dphi + TMath::ASin(-0.075*fMagsign*dist/pti) - TMath::ASin(-0.075*fMagsign*dist/ptj);
							dphi_s[idist] = fOrder*dphi_s[idist];
							dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
						}
						dphi = fOrder*dphi;
  
						//3 sigma cut
						if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;
  
						Float_t	sharity;
						if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
						else						sharity = Sharity(clu2, clu1, sha2, sha1);
  
						if (sharity >= 0.05) continue;
  
						Double_t qinv, c_weight;
						Double_t qside, qout, qlong;
						if (fOrder > 0) {
							qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
							qout 		= qout_cms(pxi, pxj, pyi, pyj);
							qside 	= qside_cms(pxi, pxj, pyi, pyj);
							qlong 	= qlong_cms(pzi, pzj, pei, pej);
						}
						else {
							qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
							qout 		= qout_cms(pxj, pxi, pyj, pyi);
							qside 	= qside_cms(pxj, pxi, pyj, pyi);
							qlong 	= qlong_cms(pzj, pzi, pej, pei);
						}

						if (fabs(qinv) <= 0.34) {
							for (int iqbin = 0; iqbin < fQbin; iqbin++) {
								fH1Qinv_mix	[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
							}
						}
						if (conv_over(qout, qside, qlong) == 1) continue;
  
						for (int iqbin = 0; iqbin < fQbin; iqbin++) {
							fH3Q_mix		[0][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
						}
  
					} //j_loop
				}	//i_loop
			}	//jeve
		}	//ieve
  
		for (int ieve = 0; ieve < fMix; ieve++) {
			pion_ppx	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_ppy	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_ppz	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_ppt	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_pe		[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_pphi	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_peta	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_pclu	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_psha	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_pax	[fCent][fZvtx][fnEP3FMD][ieve].clear();
		}
		mix_flag[0][fCent][fZvtx][fnEP3FMD] = 0;
	} //Pos charge Mixing
  
	if (mix_flag[1][fCent][fZvtx][fnEP3FMD] == fMix) {
		cout << "*****************************" << endl;
		cout << "* -- Real Pair Calculation  *" << endl;
		cout << "*****************************" << endl;
		for (int ieve = 0; ieve < fMix; ieve++) {
			for (ULong_t i_loop = 0; i_loop < pion_mpx [fCent][fZvtx][fnEP3FMD][ieve].size()-1; i_loop++) {
				for (ULong_t j_loop = i_loop+1; j_loop < pion_mpx [fCent][fZvtx][fnEP3FMD][ieve].size(); j_loop++) {
  
					Double_t pxi		= pion_mpx	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pxj		= pion_mpx	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pyi		= pion_mpy	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pyj		= pion_mpy	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pzi		= pion_mpz	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pzj		= pion_mpz	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pei		= pion_me		[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pej		= pion_me		[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t pti		= pion_mpt	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t ptj		= pion_mpt	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t phii		= pion_mphi	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t phij		= pion_mphi	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t etai		= pion_meta	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t etaj		= pion_meta	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					TBits		clu1		= pion_mclu	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		clu2		= pion_mclu	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					TBits		sha1		= pion_msha	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		sha2		= pion_msha	[fCent][fZvtx][fnEP3FMD][ieve].at(j_loop);
					Double_t EPi		= pion_max	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);
  
					Double_t dphi	= phii - phij;
					Double_t deta	= etai - etaj;
  
					Double_t kx = 0.5*(pxi + pxj);
					Double_t ky = 0.5*(pyi + pyj);
					Double_t kt = sqrt(kx*kx + ky*ky);
					fH1kt_all -> Fill(kt);
	
					Int_t	fkT = (int)(kt*10.0) - 2;
					if (fkT < 0 || fkT > 17) continue;
					fH1kt -> Fill(kt);
  
					Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
					Double_t dphi_E	= dphi_p - EPi;

					dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

					Int_t	fDphi_E	= PiEPAngle(dphi_E);
					if (fDphi_E < 0) continue;
  
					Int_t	fOrder;
					if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;
  
					//Calculate Phi* ( distance = r[m] )
					Double_t dphi_s[2];
					for (int idist = 0; idist < 2; idist++) {
						Double_t dist = dec_dist(idist);
						if (dist < 0) continue;
						dphi_s[idist]	= dphi + TMath::ASin(0.075*fMagsign*dist/pti) - TMath::ASin(0.075*fMagsign*dist/ptj);
						dphi_s[idist] = fOrder*dphi_s[idist];
						dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
					}
					dphi = fOrder*dphi;
  
					//3 sigma cut
					if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;
  
					Float_t	sharity;
					if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
					else						sharity = Sharity(clu2, clu1, sha2, sha1);
  
					if (sharity >= 0.05) continue;
  
					Double_t qinv, c_weight;
					Double_t qside, qout, qlong;
					if (fOrder > 0) {
						qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
						qout 		= qout_cms(pxi, pxj, pyi, pyj);
						qside 	= qside_cms(pxi, pxj, pyi, pyj);
						qlong 	= qlong_cms(pzi, pzj, pei, pej);
					}
					else {
						qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
						qout 		= qout_cms(pxj, pxi, pyj, pyi);
						qside 	= qside_cms(pxj, pxi, pyj, pyi);
						qlong 	= qlong_cms(pzj, pzi, pej, pei);
					}
  
					if (fabs(qinv) <= 0.34) {
						if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
						else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
						for (int iqbin = 0; iqbin < fQbin; iqbin++) {
							fH1Qinv	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
							fH1CQinv[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
						}
					}
  
					if (conv_over(qout, qside, qlong) == 1) continue;
  
					for (int iqbin = 0; iqbin < fQbin; iqbin++) {
						fH3Q[1][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
					}
  
					Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
					Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
					Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
					if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
						fPrfConv	[1][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql,	qinv);
						//fPrfConv	[1][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2	+ 400*conv_ql2,	qinv);
					}
					if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
						fPrfConv	[1][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3	+ 3600*conv_ql3,	qinv);
					}
				}	//j_loop
			}	//i_loop
		}	//ieve
  
		cout << "***************************" << endl;
		cout << "* -- Mix Pair Calculation *" << endl;
		cout << "***************************" << endl;
  
		for (int ieve = 0; ieve < fMix; ieve++) {
			for (int jeve = ieve+1; jeve < fMix; jeve++) {
				for (ULong_t i_loop = 0; i_loop < pion_mpx [fCent][fZvtx][fnEP3FMD][ieve].size(); i_loop++) {
					for (ULong_t j_loop = 0; j_loop < pion_mpx [fCent][fZvtx][fnEP3FMD][jeve].size(); j_loop++) {
  
						Double_t pxi		= pion_mpx	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pxj		= pion_mpx	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pyi		= pion_mpy	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pyj		= pion_mpy	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pzi		= pion_mpz	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pzj		= pion_mpz	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pei		= pion_me		[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t pej		= pion_me		[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t pti		= pion_mpt	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t ptj		= pion_mpt	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t phii		= pion_mphi	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t phij		= pion_mphi	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t etai		= pion_meta	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	Double_t etaj		= pion_meta	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						TBits		clu1		= pion_mclu	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		clu2		= pion_mclu	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						TBits		sha1		= pion_msha	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);	TBits		sha2		= pion_msha	[fCent][fZvtx][fnEP3FMD][jeve].at(j_loop);
						Double_t EPi		= pion_max	[fCent][fZvtx][fnEP3FMD][ieve].at(i_loop);
  
						Double_t dphi	= phii - phij;
						Double_t deta	= etai - etaj;
  
						Double_t kx = 0.5*(pxi + pxj);
						Double_t ky = 0.5*(pyi + pyj);
						Double_t kt = sqrt(kx*kx + ky*ky);
						//fH1kt_all -> Fill(kt);
		
						Int_t	fkT = (int)(kt*10.0) - 2;
						if (fkT < 0 || fkT > 17) continue;
  
						Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
						Double_t dphi_E	= dphi_p - EPi;

						dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

						Int_t	fDphi_E	= PiEPAngle(dphi_E);
						if (fDphi_E < 0) continue;
  
						Int_t	fOrder;
						if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;
  
						//Calculate Phi* ( distance = r[m] )
						Double_t dphi_s[2];
						for (int idist = 0; idist < 2; idist++) {
							Double_t dist = dec_dist(idist);
							if (dist < 0) continue;
							dphi_s[idist]	= dphi + TMath::ASin(0.075*fMagsign*dist/pti) - TMath::ASin(0.075*fMagsign*dist/ptj);
							dphi_s[idist] = fOrder*dphi_s[idist];
							dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
						}
						dphi = fOrder*dphi;
  
						//3 sigma cut
						if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;
  
						Float_t	sharity;
						if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
						else						sharity = Sharity(clu2, clu1, sha2, sha1);
  
						if (sharity >= 0.05) continue;
  
						Double_t qinv, c_weight;
						Double_t qside, qout, qlong;
						if (fOrder > 0) {
							qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
							qout 		= qout_cms(pxi, pxj, pyi, pyj);
							qside 	= qside_cms(pxi, pxj, pyi, pyj);
							qlong 	= qlong_cms(pzi, pzj, pei, pej);
						}
						else {
							qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
							qout 		= qout_cms(pxj, pxi, pyj, pyi);
							qside 	= qside_cms(pxj, pxi, pyj, pyi);
							qlong 	= qlong_cms(pzj, pzi, pej, pei);
						}

						if (fabs(qinv) <= 0.34) {
							for (int iqbin = 0; iqbin < fQbin; iqbin++) {
								fH1Qinv_mix	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
							}
						}
  
						if (conv_over(qout, qside, qlong) == 1) continue;
  
						for (int iqbin = 0; iqbin < fQbin; iqbin++) {
							fH3Q_mix	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
						}
					} //j_loop
				}	//i_loop
			}	//jeve
		}	//ieve
  
		for (int ieve = 0; ieve < fMix; ieve++) {
			pion_mpx	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_mpy	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_mpz	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_mpt	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_me		[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_mphi	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_meta	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_mclu	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_msha	[fCent][fZvtx][fnEP3FMD][ieve].clear();
			pion_max	[fCent][fZvtx][fnEP3FMD][ieve].clear();
		}
		mix_flag[1][fCent][fZvtx][fnEP3FMD] = 0;
	} //neg charge Mixing

	if (fNumberE == fNEntries)	{
		PairCalc_Final()					;
		return;
	}
	else {
		PostData(1, fOutputList)	;
	}
}

//________________________________________________________________________
void AliAnalysisTaskPsi3HBTsystFitrange::Terminate(Option_t *)
{

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}


Int_t AliAnalysisTaskPsi3HBTsystFitrange::TPCeta(Float_t ieta){
	//Eta	for TPC
	if				(ieta >= -1.0	&& ieta < -0.5) return	0;
	else	if	(ieta >= -0.5 && ieta < 0.0)	return	1;
	else	if	(ieta >= 0.0	&& ieta < 0.5)	return	2;
	else	if	(ieta >= 0.5	&& ieta <= 1.0)	return	3;
	else {
		printf("Eta is Out of Range\n");
		return -999;
	}
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::TPCeta2(Float_t ieta){
	//Eta	for TPC
	if				(ieta >= -1.0	&& ieta <	0.0) 	return	4;
	else	if	(ieta >=  0.0 && ieta <= 1.0)	return	5;
	else {
		printf("Eta is Out of Range\n");
		return -999;
	}
}
Float_t AliAnalysisTaskPsi3HBTsystFitrange::PhiZDC(Int_t itow){
	//Phi for ZDC
	Int_t sec[8]={2,3,1,0,3,2,0,1};
	if(itow>=0 && itow<8) return TMath::Pi()*(45.+ 90*sec[itow])/180;
	else printf("Error in PhiZDC() === Out of Range\n");
	return -999;
}
Int_t	AliAnalysisTaskPsi3HBTsystFitrange::trig_dec(Bool_t ikMB, Bool_t ikSemi, Bool_t ikCent) {
	if      (ikMB == true)    return 0;
	else if (ikSemi == true)  return 1;
	else if (ikCent == true)  return 2;
	else                      return -999;
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::HistCent(Double_t fV0Cent) {
	if (fV0Cent < 50.0) {
		if			(fV0Cent <=	10.0)	return	(Int_t)(cent_V0M / 5.0	)		;
		else if (fV0Cent > 10.0)	return	(Int_t)(cent_V0M / 10.0	)+1	;
		else											return -1;
	}
	else if		(fV0Cent < 80.0)	return	6;
	else												return -1;
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::CoulCent(Int_t ftV0Cent) {
	if			(ftV0Cent < 2)		return ftV0Cent;
	else if (ftV0Cent < 10)		return ( (Int_t) (ftV0Cent/2.0) + 1);
	else											return 6;
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::PiEPAngle(Double_t fPlane) {
	if      ( (fPlane > -PI/24.0		)		&& (fPlane <= PI/24.0			) )																														return 0;
	else if ( (fPlane > PI/24.0			)		&& (fPlane <= 3.0*PI/24.0	) )																														return 1;
	else if ( (fPlane > 3.0*PI/24.0	)		&& (fPlane <= 5.0*PI/24.0	) )																														return 2;
	else if ( (fPlane > 5.0*PI/24.0	)		&& (fPlane <= 7.0*PI/24.0	) )																														return 3;
	else if ( ( (fPlane > 7.0*PI/24.0	) && (fPlane <= PI/3.0	) )	|| ( (fPlane > -PI/3.0	)   && (fPlane <= -7.0*PI/24.0	) ) )	return 4;
	else if ( (fPlane > -7.0*PI/24.0)		&& (fPlane <= -5.0*PI/24.0) )																														return 5;
	else if ( (fPlane > -5.0*PI/24.0)		&& (fPlane <= -3.0*PI/24.0) )																														return 6;
	else if ( (fPlane > -3.0*PI/24.0)		&& (fPlane <= -PI/24.0		) )																														return 7;
	else																																																												return -1;
}
Double_t AliAnalysisTaskPsi3HBTsystFitrange::qinv_cms(double& px1, double& px2, double& py1, double& py2, double& pz1, double& pz2, double& pe1, double& pe2)
{
	double qx = px1 - px2;
	double qy = py1 - py2;
	double qz = pz1 - pz2;
	double qe = pe1 - pe2;
	//float beta = (pzi + pzj)/(pei + pej);
	//float gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));
	//float ql = gamma*(qz - beta*qe);
	//float qe_b = gamma*(qe - beta*qz);

	double tmp_qinv = sqrt(fabs(qx*qx + qy*qy + qz*qz - qe*qe));

	return tmp_qinv;
}


Double_t AliAnalysisTaskPsi3HBTsystFitrange::qout_cms(double& x1, double& x2, double& y1, double& y2)
{
	double dx = x1 - x2;		double xt = x1 + x2;
	double dy = y1 - y2;		double yt = y1 + y2;

	double k1 = (TMath::Sqrt(xt*xt + yt*yt));
	double k2 = (dx*xt + dy*yt);

	double tmp_qout = k2/k1;

	return tmp_qout;
}


Double_t AliAnalysisTaskPsi3HBTsystFitrange::qside_cms(double& x1, double& x2, double& y1, double& y2)
{
	double xt = x1 + x2;		double yt = y1 + y2;

	double k1 = TMath::Sqrt(xt*xt + yt*yt);

	double tmp_qside = 2.0*(x2*y1-x1*y2)/k1;

	return tmp_qside;
}

Double_t AliAnalysisTaskPsi3HBTsystFitrange::qlong_cms(double& pz1, double& pz2, double& pe1, double& pe2)
{
	double beta = (pz1 + pz2)/(pe1 + pe2);
	double gamma = 1.0/TMath::Sqrt((1.-beta)*(1.+beta));
	double qz = pz1 - pz2;
	double qe = pe1 - pe2;

	double tmp_qlong = gamma*(qz - beta*qe);

	return (tmp_qlong);
}

int AliAnalysisTaskPsi3HBTsystFitrange::MagDec(Float_t fMag_t)
{
	if			(fMag_t > 1.0)	return	1;
	else if	(fMag_t < 1.0)	return	-1;
	else										return	0;
}

float AliAnalysisTaskPsi3HBTsystFitrange::dec_dist(float fDist) {
	if			( fDist == 0 )	return 1.1		;
	else if ( fDist == 1 )	return 2.45		;
	else										return -999.0	;
}

int AliAnalysisTaskPsi3HBTsystFitrange::conv_if(int qoIn, int qsIn, int qlIn) {
	//if (qoIn >= 0 && qoIn < 30 && qsIn >= 0 && qsIn < 30 && qlIn >= 0 && qlIn < 30)	return 0		;
	if (qoIn >= 0 && qoIn < 40 && qsIn >= 0 && qsIn < 40 && qlIn >= 0 && qlIn < 40)	return 0		;
	else																																						return 1		;
}

int AliAnalysisTaskPsi3HBTsystFitrange::conv_if2(int qoIn, int qsIn, int qlIn) {
	if (qoIn >= 0 && qoIn < 60 && qsIn >= 0 && qsIn < 60 && qlIn >= 0 && qlIn < 60)	return 0		;
	else																																						return 1		;
}

int AliAnalysisTaskPsi3HBTsystFitrange::conv_over(double qoIn, double qsIn, double qlIn) {
	//if (fabs(qinvIn) > 1.0 && fabs(qoIn) > 0.3 && fabs(qsIn) > 0.3 && fabs(qlIn) > 0.3) return 1;
	//if (fabs(qoIn) > 0.3 || fabs(qsIn) > 0.3 || fabs(qlIn) > 0.3) return 1;
	if (fabs(qoIn) > 0.3 || fabs(qsIn) > 0.3 || fabs(qlIn) > 0.3) return 1;
	else return 0;
}

Int_t AliAnalysisTaskPsi3HBTsystFitrange::conv_qosl(double qIn) {
	//return (Int_t) (qIn*100.0 + 15.0);
	return (Int_t) (qIn*100.0 + 20.0);
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::conv_qosl2(double qIn) {
	//return (Int_t) (qIn*100.0 + 15.0);
	return (Int_t) (qIn*50.0 + 10.0);
}
Int_t AliAnalysisTaskPsi3HBTsystFitrange::conv_qosl3(double qIn) {
	//return (Int_t) (qIn*100.0 + 15.0);
	return (Int_t) (qIn*200.0 + 30.0);
}

float AliAnalysisTaskPsi3HBTsystFitrange::Sharity(TBits& cl1,  TBits& cl2, TBits& sh1, TBits& sh2) {

	int ncl1 = cl1.GetNbits();
	int ncl2 = cl2.GetNbits();
	int sumCls = 0; int sumSha = 0; int sumQ = 0;
	double shfrac = 0;
	for(int ich = 0; ich < ncl1 && ich < ncl2; ich++) {
		if (cl1.TestBitNumber(ich) && cl2.TestBitNumber(ich)) { // Both clusters
			if (sh1.TestBitNumber(ich) && sh2.TestBitNumber(ich)) { // Shared
				sumQ++;
				sumCls+=2;
				sumSha+=2;
			}
			else {sumQ--; sumCls+=2;}
		}
		else if (cl1.TestBitNumber(ich) || cl2.TestBitNumber(ich)) { // Non shared
			sumQ++;
			sumCls++;
		}
	}
	if (sumCls>0) {
		shfrac = sumSha*1.0/sumCls;
	}

	return shfrac;
}

float AliAnalysisTaskPsi3HBTsystFitrange::Qfac(TBits& cl1,  TBits& cl2, TBits& sh1, TBits& sh2) {

	int ncl1 = cl1.GetNbits();
	int ncl2 = cl2.GetNbits();
	int sumCls = 0; int sumSha = 0; int sumQ = 0;
	double qfactor = 0;
	for(int ipad = 0; ipad < ncl1 && ipad < ncl2; ipad++) {
		if (cl1.TestBitNumber(ipad) && cl2.TestBitNumber(ipad)) { // Both clusters
			if (sh1.TestBitNumber(ipad) && sh2.TestBitNumber(ipad)) { // Shared
				sumQ++;
				sumCls+=2;
				sumSha+=2;
			}
			else {sumQ--; sumCls+=2;}
		}
		else if (cl1.TestBitNumber(ipad) || cl2.TestBitNumber(ipad)) { // Non shared
			sumQ++;
			sumCls++;
		}
	}
	if (sumCls>0) {
		qfactor = sumQ*1.0/sumCls;
	}

	return qfactor;
}


bool AliAnalysisTaskPsi3HBTsystFitrange::IsPionNSigma(Double_t& mom, Float_t& nsigmaTPCPi, Float_t& nsigmaTOFPi)
{
	if(mom<0.65)
	{
		if(nsigmaTOFPi<-999.)
		{
			if(mom<0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
			else if(mom<0.5 && mom>=0.35 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
			else if(mom>=0.5 && TMath::Abs(nsigmaTPCPi)<2.0) return true;
			else return false;
		}
		else if(TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<3.0) return true;
	}
	else if(nsigmaTOFPi<-999.)
	{
		return false;
	}
	else if(mom<1.5 && TMath::Abs(nsigmaTOFPi)<3.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;
	else if(mom>=1.5 && TMath::Abs(nsigmaTOFPi)<2.0 && TMath::Abs(nsigmaTPCPi)<5.0) return true;

	return false;
}

bool AliAnalysisTaskPsi3HBTsystFitrange::IsKaonNSigma(Double_t& mom, Float_t& nsigmaTPCK, Float_t& nsigmaTOFK)
{

	if(mom<0.4)
	{
		if(nsigmaTOFK<-999.)
		{
			if(TMath::Abs(nsigmaTPCK)<2.0) return true;
		}
		else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
	}
	else if(mom>=0.4 && mom<=0.6)
	{
		if(nsigmaTOFK<-999.)
		{
			if(TMath::Abs(nsigmaTPCK)<2.0) return true;
		}
		else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;
	}
	else if(nsigmaTOFK<-999.)
	{
		return false;
	}
	else if(TMath::Abs(nsigmaTOFK)<3.0 && TMath::Abs(nsigmaTPCK)<3.0) return true;

	return false;
}

bool AliAnalysisTaskPsi3HBTsystFitrange::IsProtonNSigma(Double_t& mom, Float_t& nsigmaTPCP, Float_t& nsigmaTOFP)
{
	if (mom > 0.8) {
		if (TMath::Hypot( nsigmaTOFP, nsigmaTPCP )/TMath::Sqrt(2) < 3.0)
			return true;
	}
	else {
		if (TMath::Abs(nsigmaTPCP) < 3.0)
			return true;
	}

	return false;
}

Float_t  AliAnalysisTaskPsi3HBTsystFitrange::EtaLimitFMDC(Int_t i){Double_t eta[]={-3.4, -2.8, -2.2, -1.7};             return eta[i];}
Float_t  AliAnalysisTaskPsi3HBTsystFitrange::EtaLimitFMDA(Int_t i){Double_t eta[]={ 1.7,  2.2,  2.8,  3.4,  4.0,  5.0}; return eta[i];}

//"C1", "C2", "C3", "C", "A1", "A2", "A3", "A4", "A5", "A", "A123", "A45", "AC"
Float_t  AliAnalysisTaskPsi3HBTsystFitrange::EtaLimitFMDAC1(Int_t i){Double_t eta[]={ -3.4,  -2.8,  -2.2,  -3.4,  1.7,  2.2,  2.8,	3.4,	4.0,	1.7,	1.7,	3.4,	-10.0}; return eta[i];}
Float_t  AliAnalysisTaskPsi3HBTsystFitrange::EtaLimitFMDAC2(Int_t i){Double_t eta[]={ -2.8,  -2.2,  -1.7,  -1.7,  2.2,  2.8,  3.4,	4.0,	5.0,	5.0,	3.4,	5.0,	10.0};	return eta[i];}

TH1F* AliAnalysisTaskPsi3HBTsystFitrange::ProjectPhiFMD(TH2D* fmd, Float_t eta1, Float_t eta2){   //don't forget delete projected hist
	Int_t bin1,bin2;
	bin1 = fmd->GetXaxis()->FindBin(eta1);
	bin2 = fmd->GetXaxis()->FindBin(eta2)+1;
	return (TH1F*)fmd->ProjectionY(Form("fFMDphi_bin%d_bin%d",bin1,bin2),bin1,bin2);
}
void AliAnalysisTaskPsi3HBTsystFitrange::FillEPCorr(Int_t fV0flag, Int_t fTPCflag, Int_t fFMDflag){   //don't forget delete projected hist
	for (int iharm = 0; iharm < fHarm; iharm++) {
		//-------------------------------------------------------------------------------------------------------------
		//[correct][side][harmonics]	
		//-------------------------------------------------------------------------------------------------------------
		//	===VZERO================
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[0]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0A		- TPC          
		if (fV0flag	== 0)							 		 		fPrfCos	[1]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fV0psi	[1][0][iharm]		)	)	);      //			+ V0A		- V0C          
		if (fV0flag	== 0 && fFMDflag	== 0)		fPrfCos	[2]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fFMDpsi	[1][3][iharm]		)	)	);      //			+ V0A		- FMDC         
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[3]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0C		- TPC          
    if (fV0flag	== 0 && fFMDflag	== 0)		fPrfCos	[4]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fFMDpsi	[1][9][iharm]		)	)	);      //			+ V0C		- FMDA         
    if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[5]	[iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0AC	- TPC          
		//-------------------------------------------------------------------------------------------------------------
		//	===FMD==================                                                              
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[6][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDA	- TPC      
		if (fFMDflag	== 0)							 			fPrfCos	[7][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fFMDpsi	[1][3][iharm]		)	)	);      //			+ FMDA	- FMDC     
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[8][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDC	- TPC      
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[9][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][12]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDAC	- TPC      
		//	===TPC Divided=========
		//	=== with V0A  =========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[10][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0A		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[11][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ V0A		- TPCC1(-0.5, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[12][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ V0A		- TPCA1( 0.0, 0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[13][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0A		- TPCA2( 0.5, 1.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[14][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ V0A		- TPCC (-1.0, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[15][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ V0A		- TPCA ( 0.0, 1.0)          
		//	=== with V0C  =========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[16][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0C		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[17][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ V0C		- TPCC1(-0.5, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[18][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ V0C		- TPCA1( 0.0, 0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[19][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0C		- TPCA2( 0.5, 1.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[20][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ V0C		- TPCC (-1.0, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[21][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ V0C		- TPCA ( 0.0, 1.0)          
		//	=== with V0A+C=========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[22][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0C		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfCos	[23][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0C		- TPCA2( 0.5, 1.0)          
		//	=== with FMDA  =========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[24][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDA	- TPCC2(-1.0,-0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[25][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ FMDA	- TPCC1(-0.5, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[26][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ FMDA	- TPCA1( 0.0, 0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[27][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDA	- TPCA2( 0.5, 1.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[28][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ FMDA	- TPCC (-1.0, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[29][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ FMDA	- TPCA ( 0.0, 1.0)          
		//	=== with FMDC  =========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[30][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDC	- TPCC2(-1.0,-0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[31][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ FMDC	- TPCC1(-0.5, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[32][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ FMDC	- TPCA1( 0.0, 0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[33][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDC	- TPCA2( 0.5, 1.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[34][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ FMDC	- TPCC (-1.0, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[35][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ FMDC	- TPCA ( 0.0, 1.0)          
		//	=== with FMDA+C=========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[36][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][12][iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDAC	- TPCC2(-1.0,-0.5)     
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfCos	[37][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fFMDpsi	[1][12][iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDAC	- TPCA2( 0.5, 1.0)          
		//	=== inside TPC =========
		if (fTPCflag	== 0)							 		 	fPrfCos	[38][iharm] -> Fill(	cent_V0M,	TMath::Cos(	(iharm+1.0)*	(	fTPCpsi	[1][0]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ TPCC2	- TPCA2          
	}

	for (int iharm = 0; iharm < fHarm; iharm++) {
		//-------------------------------------------------------------------------------------------------------------
		//[correct][side][harmonics]	
		//-------------------------------------------------------------------------------------------------------------
		//	===VZERO================
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[0]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0A		- TPC          
		if (fV0flag	== 0)							 		 		fPrfSin	[1]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fV0psi	[1][0][iharm]		)	)	);      //			+ V0A		- V0C          
		if (fV0flag	== 0 && fFMDflag	== 0)		fPrfSin	[2]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fFMDpsi	[1][3][iharm]		)	)	);      //			+ V0A		- FMDC         
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[3]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0C		- TPC          
    if (fV0flag	== 0 && fFMDflag	== 0)		fPrfSin	[4]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fFMDpsi	[1][9][iharm]		)	)	);      //			+ V0C		- FMDA         
    if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[5]	[iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ V0AC	- TPC          
		//-------------------------------------------------------------------------------------------------------------
		//	===FMD==================                                                               	  
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[6][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDA	- TPC      
		if (fFMDflag	== 0)							 			fPrfSin	[7][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fFMDpsi	[1][3][iharm]		)	)	);      //			+ FMDA	- FMDC     
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[8][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDC	- TPC      
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[9][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][12]	[iharm]	-	fTPCpsi	[1][6][iharm]		)	)	);      //			+ FMDAC	- TPC      
		//-------------------------------------------------------------------------------------------------------------
		//	===TPC Divided=========
		//	=== with V0A  =========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[10][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0A		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[11][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ V0A		- TPCC1(-0.5, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[12][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ V0A		- TPCA1( 0.0, 0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[13][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0A		- TPCA2( 0.5, 1.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[14][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ V0A		- TPCC (-1.0, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[15][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][1]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ V0A		- TPCA ( 0.0, 1.0)          
		//	=== with V0C  =========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[16][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0C		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[17][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ V0C		- TPCC1(-0.5, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[18][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ V0C		- TPCA1( 0.0, 0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[19][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0C		- TPCA2( 0.5, 1.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[20][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ V0C		- TPCC (-1.0, 0.0)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[21][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][0]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ V0C		- TPCA ( 0.0, 1.0)          
		//	=== with V0A+C=========
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[22][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ V0C		- TPCC2(-1.0,-0.5)          
		if (fV0flag	== 0 && fTPCflag	== 0)		fPrfSin	[23][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fV0psi	[1][2]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ V0C		- TPCA2( 0.5, 1.0)          
		//	=== with FMDA  =========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[24][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDA	- TPCC2(-1.0,-0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[25][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ FMDA	- TPCC1(-0.5, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[26][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ FMDA	- TPCA1( 0.0, 0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[27][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDA	- TPCA2( 0.5, 1.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[28][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ FMDA	- TPCC (-1.0, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[29][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][9]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ FMDA	- TPCA ( 0.0, 1.0)          
		//	=== with FMDC  =========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[30][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDC	- TPCC2(-1.0,-0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[31][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][1][iharm]		)	)	);      //			+ FMDC	- TPCC1(-0.5, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[32][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][2][iharm]		)	)	);      //			+ FMDC	- TPCA1( 0.0, 0.5)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[33][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDC	- TPCA2( 0.5, 1.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[34][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][4][iharm]		)	)	);      //			+ FMDC	- TPCC (-1.0, 0.0)          
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[35][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][3]	[iharm]	-	fTPCpsi	[1][5][iharm]		)	)	);      //			+ FMDC	- TPCA ( 0.0, 1.0)          
		//	=== with FMDA+C=========
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[36][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][12][iharm]	-	fTPCpsi	[1][0][iharm]		)	)	);      //			+ FMDAC	- TPCC2(-1.0,-0.5)     
		if (fFMDflag	== 0 && fTPCflag	== 0)	fPrfSin	[37][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fFMDpsi	[1][12][iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ FMDAC	- TPCA2( 0.5, 1.0)          
		//	=== inside TPC =========
		if (fTPCflag	== 0)							 		 	fPrfSin	[38][iharm] -> Fill(	cent_V0M,	TMath::Sin(	(iharm+1.0)*	(	fTPCpsi	[1][0]	[iharm]	-	fTPCpsi	[1][3][iharm]		)	)	);      //			+ TPCC2	- TPCA2          
	}
}

void AliAnalysisTaskPsi3HBTsystFitrange::PairCalc_Final(){
	//cout << "$$$$$$$$$$$$$$$$$$$$$$" << endl;
	//cout << "$ Pair Calc Final !! $" << endl;
	//cout << "$$$$$$$$$$$$$$$$$$$$$$" << endl;
	//cout << endl << endl << endl << endl << endl;
	//Coulomb source size Input
	Float_t r_inv		[7]	= {11, 10, 9, 8, 7, 6, 5};
	Float_t r_side	[7]	= {11, 10, 9, 8, 7, 6, 5};
	Float_t r_out		[7]	= {11, 10, 9, 8, 7, 6, 5};
	Float_t r_long	[7]	= {11, 10, 9, 8, 7, 6, 5};

	if (fNumberE == fNEntries) {
		cout << "$$$$$$$$$$$$$$$$$$$$$$" << endl;
		for (int tcharge = 0; tcharge < fChargeH; tcharge++) {
			for (int ftempCent = 0; ftempCent < fCentM; ftempCent++) {
				for (int ftZvtx = 0; ftZvtx < fZvtxM; ftZvtx++) {
					for (int ftEP3V0 = 0; ftEP3V0 < fPsi3M; ftEP3V0++) {
						if (mix_flag[tcharge][ftempCent][ftZvtx][ftEP3V0]) {
							cout << "Mixing Flag[" << tcharge << "][" << ftempCent << "][" << ftZvtx << "][" << ftEP3V0 << "] : " << mix_flag[tcharge][ftempCent][ftZvtx][ftEP3V0] << endl;
						}	//flag != 0
					}	//Psi3 loop
				}	//Zvtx loop
			}	//Centrality loop
		}	//charge loop
		cout << "$$$$$$$$$$$$$$$$$$$$$$" << endl;
		cout << endl << endl << endl << endl << endl;

		for (int ftempCent = 0; ftempCent < fCentM; ftempCent++) {
			Int_t fCent_H = CoulCent(ftempCent);
			//if (!fCentralAna)	fCent_H = fCent_H-2;
			if (!fCentralAna) {
				fCent_H = fCent_H-2;
				cout << "Centrality Bin : " << fCent_H << endl;
				cout << "Coulomb Correction : " << r_inv[fCent_H+2] << endl;
			}

			for (int ftZvtx = 0; ftZvtx < fZvtxM; ftZvtx++) {
				for (int ftEP3V0 = 0; ftEP3V0 < fPsi3M; ftEP3V0++) {
					if (mix_flag[0][ftempCent][ftZvtx][ftEP3V0]>0 && mix_flag[1][ftempCent][ftZvtx][ftEP3V0]>0) {
						cout << "GGGGGGGGGGGGGGGGGGGGGG" << endl;
						cout << "Mixing Flag[0][" << ftempCent << "][" << ftZvtx << "][" << ftEP3V0 << "] : " <<mix_flag[0][ftempCent][ftZvtx][ftEP3V0] << endl;
						cout << "Mixing Flag[1][" << ftempCent << "][" << ftZvtx << "][" << ftEP3V0 << "] : " <<mix_flag[1][ftempCent][ftZvtx][ftEP3V0] << endl;
						cout << endl << endl;
					}
					if ( ( mix_flag[0][ftempCent][ftZvtx][ftEP3V0] < 10 ) && ( mix_flag[0][ftempCent][ftZvtx][ftEP3V0] > 1 )) {
						cout << "*****************************" << endl;
						cout << "* ++ Real Pair Calculation  *" << endl;
						cout << "*****************************" << endl;
						for (int ieve = 0; ieve < mix_flag[0][ftempCent][ftZvtx][ftEP3V0]; ieve++) {
							for (ULong_t i_loop = 0; i_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][ieve].size()-1; i_loop++) {
								for (ULong_t j_loop = i_loop+1; j_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); j_loop++) {

									Double_t pxi		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pyi		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pzi		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pei		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pti		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t phii		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t etai		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		 clu1		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		 clu2		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		 sha1		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		 sha2		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t EPi		= pion_pax	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

									Double_t dphi	= phii - phij;
									Double_t deta	= etai - etaj;

									Double_t kx = 0.5*(pxi + pxj);
									Double_t ky = 0.5*(pyi + pyj);
									Double_t kt = sqrt(kx*kx + ky*ky);
									fH1kt_all -> Fill(kt);
				
									Int_t	fkT = (int)(kt*10.0) - 2;
									if (fkT < 0 || fkT > 17) continue;

									fH1kt -> Fill(kt);

									Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
									Double_t dphi_E	= dphi_p - EPi;

									dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

									Int_t	fDphi_E	= PiEPAngle(dphi_E);
									if (fDphi_E < 0) continue;

									Int_t	fOrder;
									if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

									//Calculate Phi* ( distance = r[m] )
									Double_t dphi_s[2];
									for (int idist = 0; idist < 2; idist++) {
										Double_t dist = dec_dist(idist);
										if (dist < 0) continue;
										dphi_s[idist]	= dphi + TMath::ASin(-0.075*fMagsign*dist/pti) - TMath::ASin(-0.075*fMagsign*dist/ptj);
										dphi_s[idist] = fOrder*dphi_s[idist];
										dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
									}
									dphi = fOrder*dphi;

									//3 sigma cut
									if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

									Float_t	sharity;
									if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
									else						sharity = Sharity(clu2, clu1, sha2, sha1);

									if (sharity >= 0.05) continue;

									Double_t qinv, c_weight;
									Double_t qside, qout, qlong;
									if (fOrder > 0) {
										qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
										qout 		= qout_cms(pxi, pxj, pyi, pyj);
										qside 	= qside_cms(pxi, pxj, pyi, pyj);
										qlong 	= qlong_cms(pzi, pzj, pei, pej);
									}
									else {
										qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
										qout 		= qout_cms(pxj, pxi, pyj, pyi);
										qside 	= qside_cms(pxj, pxi, pyj, pyi);
										qlong 	= qlong_cms(pzj, pzi, pej, pei);
									}

									if (fabs(qinv) <= 0.34) {
										if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
                    else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH1Qinv	[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											fH1CQinv[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
										}
									}

									if (conv_over(qout, qside, qlong) == 1) continue;

									for (int iqbin = 0; iqbin < fQbin; iqbin++) {
										fH3Q[0][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
									}

									Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
									Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
									Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
									if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
										fPrfConv	[0][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql, qinv);
										//fPrfConv	[0][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2 + 400*conv_ql2, qinv);
									}
									if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
										fPrfConv	[0][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3 + 3600*conv_ql3, qinv);
									}
								}	//j_loop
							}	//i_loop
						}	//ieve

						cout << "***************************" << endl;
						cout << "* ++ Mix Pair Calculation *" << endl;
						cout << "***************************" << endl;

						for (int ieve = 0; ieve < mix_flag[0][ftempCent][ftZvtx][ftEP3V0]-1; ieve++) {
							for (int jeve = ieve+1; jeve < mix_flag[0][ftempCent][ftZvtx][ftEP3V0]; jeve++) {
								for (ULong_t i_loop = 0; i_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); i_loop++) {
									for (ULong_t j_loop = 0; j_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][jeve].size(); j_loop++) {

										Double_t pxi		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pyi		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pzi		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pei		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pti		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t phii		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t etai		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										TBits		clu1		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		clu2		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										TBits		sha1		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		sha2		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t EPi		= pion_pax	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

										Double_t dphi	= phii - phij;
										Double_t deta	= etai - etaj;

										Double_t kx = 0.5*(pxi + pxj);
										Double_t ky = 0.5*(pyi + pyj);
										Double_t kt = sqrt(kx*kx + ky*ky);

										Int_t	fkT = (int)(kt*10.0) - 2;
										if (fkT < 0 || fkT > 17) continue;

										Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
										Double_t dphi_E	= dphi_p - EPi;

										dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

										Int_t	fOrder;
										if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

										Int_t	fDphi_E	= PiEPAngle(dphi_E);
										if (fDphi_E < 0) continue;

										//Calculate Phi* ( distance = r[m] )
										Double_t dphi_s[2];
										for (int idist = 0; idist < 2; idist++) {
											Double_t dist = dec_dist(idist);
											if (dist < 0) continue;
											dphi_s[idist]	= dphi + TMath::ASin(-0.075*fMagsign*dist/pti) - TMath::ASin(-0.075*fMagsign*dist/ptj);
											dphi_s[idist] = fOrder*dphi_s[idist];
											dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
										}
										dphi = fOrder*dphi;

										//3 sigma cut
										if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

										Float_t	sharity;
										if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
										else						sharity = Sharity(clu2, clu1, sha2, sha1);

										if (sharity >= 0.05) continue;

										Double_t qinv, c_weight;
										Double_t qside, qout, qlong;
										if (fOrder > 0) {
											qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
											qout 		= qout_cms(pxi, pxj, pyi, pyj);
											qside 	= qside_cms(pxi, pxj, pyi, pyj);
											qlong 	= qlong_cms(pzi, pzj, pei, pej);
										}
										else {
											qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
											qout 		= qout_cms(pxj, pxi, pyj, pyi);
											qside 	= qside_cms(pxj, pxi, pyj, pyi);
											qlong 	= qlong_cms(pzj, pzi, pej, pei);
										}
										if (fabs(qinv) <= 0.34) {
											for (int iqbin = 0; iqbin < fQbin; iqbin++) {
												fH1Qinv_mix	[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											}
										}

										if (conv_over(qout, qside, qlong) == 1) continue;

										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH3Q_mix		[0][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
										}
									} //j_loop
								}	//i_loop
							}	//jeve
						}	//ieve

						for (int ieve = 0; ieve < mix_flag[0][ftempCent][ftZvtx][ftEP3V0]; ieve++) {
							pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pax	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
						}
						mix_flag[0][ftempCent][ftZvtx][ftEP3V0]	= 0;
					} //Pos charge Mixing
					else if (mix_flag[0][ftempCent][ftZvtx][ftEP3V0] == 1) {
						cout << "*****************************" << endl;
						cout << "* ++ Real Pair Calculation  *" << endl;
						cout << "*****************************" << endl;
						for (int ieve = 0; ieve < 1; ieve++) {
							for (ULong_t i_loop = 0; i_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][ieve].size()-1; i_loop++) {
								for (ULong_t j_loop = i_loop+1; j_loop < pion_ppx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); j_loop++) {

									Double_t pxi		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pyi		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pzi		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pei		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pti		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t phii		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t etai		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		 clu1		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		 clu2		= pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		 sha1		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		 sha2		= pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t EPi		= pion_pax	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

									Double_t dphi	= phii - phij;
									Double_t deta	= etai - etaj;

									Double_t kx = 0.5*(pxi + pxj);
									Double_t ky = 0.5*(pyi + pyj);
									Double_t kt = sqrt(kx*kx + ky*ky);
									fH1kt_all -> Fill(kt);

									Int_t	fkT = (int)(kt*10.0) - 2;
									if (fkT < 0 || fkT > 17) continue;

									fH1kt						-> Fill(kt);

									Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
									Double_t dphi_E	= dphi_p - EPi;

									dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

									Int_t	fDphi_E	= PiEPAngle(dphi_E);
									if (fDphi_E < 0) continue;

									Int_t	fOrder;
									if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

									//Calculate Phi* ( distance = r[m] )
									Double_t dphi_s[2];
									for (int idist = 0; idist < 2; idist++) {
										Double_t dist = dec_dist(idist);
										if (dist < 0) continue;
										dphi_s[idist]	= dphi + TMath::ASin(-0.075*fMagsign*dist/pti) - TMath::ASin(-0.075*fMagsign*dist/ptj);
										dphi_s[idist] = fOrder*dphi_s[idist];
										dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
									}
									dphi = fOrder*dphi;

									//3 sigma cut
									if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

									Float_t	sharity;
									if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
									else						sharity = Sharity(clu2, clu1, sha2, sha1);

									if (sharity >= 0.05) continue;

									Double_t qinv, c_weight;
									Double_t qside, qout, qlong;
									if (fOrder > 0) {
										qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
										qout 		= qout_cms(pxi, pxj, pyi, pyj);
										qside 	= qside_cms(pxi, pxj, pyi, pyj);
										qlong 	= qlong_cms(pzi, pzj, pei, pej);
									}
									else {
										qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
										qout 		= qout_cms(pxj, pxi, pyj, pyi);
										qside 	= qside_cms(pxj, pxi, pyj, pyi);
										qlong 	= qlong_cms(pzj, pzi, pej, pei);
									}

									if (fabs(qinv) <= 0.34) {
										if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
                    else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH1Qinv	[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											fH1CQinv[0][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
										}
									}

									if (conv_over(qout, qside, qlong) == 1) continue;

									for (int iqbin = 0; iqbin < fQbin; iqbin++) {
										fH3Q[0][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
									}

									Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
									Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
									Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
									if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
										fPrfConv	[0][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql, qinv);
										//fPrfConv	[0][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2 + 400*conv_ql2, qinv);
									}
									if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
										fPrfConv	[0][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3	+ 3600*conv_ql3, qinv);
									}
								}	//j_loop
							}	//i_loop
						}	//ieve
						for (int ieve = 0; ieve < 1; ieve++) {
							pion_ppx	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppy	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppz	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_ppt	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pe		[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pphi	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_peta	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pclu	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_psha	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_pax	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
						}
						mix_flag[0][ftempCent][ftZvtx][ftEP3V0] = 0;
					}
					if (( mix_flag[1][ftempCent][ftZvtx][ftEP3V0] < 10 ) && ( mix_flag[1][ftempCent][ftZvtx][ftEP3V0] > 1 )) {
						cout << "*****************************" << endl;
						cout << "* -- Real Pair Calculation  *" << endl;
						cout << "*****************************" << endl;
						for (int ieve = 0; ieve < mix_flag[1][ftempCent][ftZvtx][ftEP3V0]; ieve++) {
							for (ULong_t i_loop = 0; i_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][ieve].size()-1; i_loop++) {
								for (ULong_t j_loop = i_loop+1; j_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); j_loop++) {

									Double_t pxi		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pyi		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pzi		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pei		= pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pti		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t phii		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t etai		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		clu1		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		clu2		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		sha1		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		sha2		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t EPi		= pion_max	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

									Double_t dphi	= phii - phij;
									Double_t deta	= etai - etaj;

									Double_t kx = 0.5*(pxi + pxj);
									Double_t ky = 0.5*(pyi + pyj);
									Double_t kt = sqrt(kx*kx + ky*ky);
									fH1kt_all -> Fill(kt);

									Int_t	fkT = (int)(kt*10.0) - 2;
									if (fkT < 0 || fkT > 17) continue;

									fH1kt						-> Fill(kt);

									Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
									Double_t dphi_E	= dphi_p - EPi;

									dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

									Int_t	fDphi_E	= PiEPAngle(dphi_E);
									if (fDphi_E < 0) continue;

									Int_t	fOrder;
									if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

									//Calculate Phi* ( distance = r[m] )
									Double_t dphi_s[2];
									for (int idist = 0; idist < 2; idist++) {
										Double_t dist = dec_dist(idist);
										if (dist < 0) continue;
										dphi_s[idist]	= dphi + TMath::ASin(0.075*fMagsign*dist/pti) - TMath::ASin(0.075*fMagsign*dist/ptj);
										dphi_s[idist] = fOrder*dphi_s[idist];
										dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
									}
									dphi = fOrder*dphi;

									//3 sigma cut
									if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

									Float_t	sharity;
									if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
									else						sharity = Sharity(clu2, clu1, sha2, sha1);

									if (sharity >= 0.05) continue;

									Double_t qinv, c_weight;
									Double_t qside, qout, qlong;
									if (fOrder > 0) {
										qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
										qout 		= qout_cms(pxi, pxj, pyi, pyj);
										qside 	= qside_cms(pxi, pxj, pyi, pyj);
										qlong 	= qlong_cms(pzi, pzj, pei, pej);
									}
									else {
										qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
										qout 		= qout_cms(pxj, pxi, pyj, pyi);
										qside 	= qside_cms(pxj, pxi, pyj, pyi);
										qlong 	= qlong_cms(pzj, pzi, pej, pei);
									}

									if (fabs(qinv) <= 0.34) {
										if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
                    else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH1Qinv	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											fH1CQinv[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
										}
									}

									if (conv_over(qout, qside, qlong) == 1) continue;

									for (int iqbin = 0; iqbin < fQbin; iqbin++) {
										fH3Q[1][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
									}

									Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
									Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
									Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
									if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
										fPrfConv	[1][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql, qinv);
										//fPrfConv	[1][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2 + 400*conv_ql2, qinv);
									}
									if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
										fPrfConv	[1][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3	+ 3600*conv_ql3, qinv);
									}
								}	//j_loop
							}	//i_loop
						}	//ieve

						cout << "***************************" << endl;
						cout << "* -- Mix Pair Calculation *" << endl;
						cout << "***************************" << endl;

						for (int ieve = 0; ieve < mix_flag[1][ftempCent][ftZvtx][ftEP3V0]-1; ieve++) {
							for (int jeve = ieve+1; jeve < mix_flag[1][ftempCent][ftZvtx][ftEP3V0]; jeve++) {
								for (ULong_t i_loop = 0; i_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); i_loop++) {
									for (ULong_t j_loop = 0; j_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][jeve].size(); j_loop++) {

										Double_t pxi		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pyi		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pzi		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pei		= pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_me		[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t pti		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t phii		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t etai		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										TBits		clu1		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		clu2		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										TBits		sha1		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		sha2		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][jeve].at(j_loop);
										Double_t EPi		= pion_max	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

										Double_t dphi	= phii - phij;
										Double_t deta	= etai - etaj;

										Double_t kx = 0.5*(pxi + pxj)			;
										Double_t ky = 0.5*(pyi + pyj)			;
										Double_t kt = sqrt(kx*kx + ky*ky)	;

										Int_t	fkT = (int)(kt*10.0) - 2;
										if (fkT < 0 || fkT > 17) continue;

										Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
										Double_t dphi_E	= dphi_p - EPi;

										dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

										Int_t	fDphi_E	= PiEPAngle(dphi_E);
										if (fDphi_E < 0) continue;

										Int_t	fOrder;
										if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

										//Calculate Phi* ( distance = r[m] )
										Double_t dphi_s[2];
										for (int idist = 0; idist < 2; idist++) {
											Double_t dist = dec_dist(idist);
											if (dist < 0) continue;
											dphi_s[idist]	= dphi + TMath::ASin(0.075*fMagsign*dist/pti) - TMath::ASin(0.075*fMagsign*dist/ptj);
											dphi_s[idist] = fOrder*dphi_s[idist];
											dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
										}
										dphi = fOrder*dphi;

										//3 sigma cut
										if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

										Float_t	sharity;
										if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
										else						sharity = Sharity(clu2, clu1, sha2, sha1);

										if (sharity >= 0.05) continue;

										Double_t qinv, c_weight;
										Double_t qside, qout, qlong;
										if (fOrder > 0) {
											qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
											qout 		= qout_cms(pxi, pxj, pyi, pyj);
											qside 	= qside_cms(pxi, pxj, pyi, pyj);
											qlong 	= qlong_cms(pzi, pzj, pei, pej);
										}
										else {
											qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
											qout 		= qout_cms(pxj, pxi, pyj, pyi);
											qside 	= qside_cms(pxj, pxi, pyj, pyi);
											qlong 	= qlong_cms(pzj, pzi, pej, pei);
										}

										if (fabs(qinv) <= 0.34) {
											for (int iqbin = 0; iqbin < fQbin; iqbin++) {
												fH1Qinv_mix	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											}
										}

										if (conv_over(qout, qside, qlong) == 1) continue;

										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH3Q_mix	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
										}
									} //j_loop
								}	//i_loop
							}	//jeve
						}	//ieve

						for (int ieve = 0; ieve < mix_flag[1][ftempCent][ftZvtx][ftEP3V0]; ieve++) {
							pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_max	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
						}
						mix_flag[1][ftempCent][ftZvtx][ftEP3V0] = 0;
					} //neg charge Mixing
					else if (mix_flag[1][ftempCent][ftZvtx][ftEP3V0] == 1) {
						cout << "*****************************" << endl;
						cout << "* -- Real Pair Calculation  *" << endl;
						cout << "*****************************" << endl;
						for (int ieve = 0; ieve < 1; ieve++) {
							for (ULong_t i_loop = 0; i_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][ieve].size()-1; i_loop++) {
								for (ULong_t j_loop = i_loop+1; j_loop < pion_mpx [ftempCent][ftZvtx][ftEP3V0][ieve].size(); j_loop++) {

									Double_t pxi		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pxj		= pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pyi		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pyj		= pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pzi		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pzj		= pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pei		= pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t pej		= pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t pti		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t ptj		= pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t phii		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t phij		= pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t etai		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	Double_t etaj		= pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		clu1		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		clu2		= pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									TBits		sha1		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);	TBits		sha2		= pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].at(j_loop);
									Double_t EPi		= pion_max	[ftempCent][ftZvtx][ftEP3V0][ieve].at(i_loop);

									Double_t dphi	= phii - phij;
									Double_t deta	= etai - etaj;

									Double_t kx = 0.5*(pxi + pxj);
									Double_t ky = 0.5*(pyi + pyj);
									Double_t kt = sqrt(kx*kx + ky*ky);
									fH1kt_all -> Fill(kt);

									Int_t	fkT = (int)(kt*10.0) - 2;
									if (fkT < 0 || fkT > 17) continue;

									fH1kt						-> Fill(kt);

									Double_t dphi_p	= TMath::ATan2(pyi+pyj, pxi+pxj);
									Double_t dphi_E	= dphi_p - EPi;
									dphi_E	= ( 1.0 / 3.0 ) * TMath::ATan2( sin(3.0 * dphi_E), cos( 3.0 * dphi_E ) );

									Int_t	fDphi_E	= PiEPAngle(dphi_E);
									if (fDphi_E < 0) continue;

									Int_t	fOrder;
									if (rand()/(double)RAND_MAX < 0.50) fOrder = -1;	else	fOrder = 1;

									//Calculate Phi* ( distance = r[m] )
									Double_t dphi_s[2];
									for (int idist = 0; idist < 2; idist++) {
										Double_t dist = dec_dist(idist);
										if (dist < 0) continue;
										dphi_s[idist]	= dphi + TMath::ASin(0.075*fMagsign*dist/pti) - TMath::ASin(0.075*fMagsign*dist/ptj);
										dphi_s[idist] = fOrder*dphi_s[idist];
										dphi_s[idist]	= TMath::ATan2(sin(dphi_s[idist]), cos(dphi_s[idist]));
									}
									dphi = fOrder*dphi;

									//3 sigma cut
									if (fabs(deta) < 0.018 && fabs(dphi_s[0]) < 0.066) continue;

									Float_t	sharity;
									if (fOrder > 0)	sharity = Sharity(clu1, clu2, sha1, sha2);
									else						sharity = Sharity(clu2, clu1, sha2, sha1);

									if (sharity >= 0.05) continue;

									Double_t qinv, c_weight;
									Double_t qside, qout, qlong;
									if (fOrder > 0) {
										qinv  	= qinv_cms(pxi, pxj, pyi, pyj, pzi, pzj, pei, pej);
										qout 		= qout_cms(pxi, pxj, pyi, pyj);
										qside 	= qside_cms(pxi, pxj, pyi, pyj);
										qlong 	= qlong_cms(pzi, pzj, pei, pej);
									}
									else {
										qinv  	= qinv_cms(pxj, pxi, pyj, pyi, pzj, pzi, pej, pei);
										qout 		= qout_cms(pxj, pxi, pyj, pyi);
										qside 	= qside_cms(pxj, pxi, pyj, pyi);
										qlong 	= qlong_cms(pzj, pzi, pej, pei);
									}

									if (fabs(qinv) <= 0.34) {
										if (!fCentralAna)	c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H+2], r_side[fCent_H+2],	r_out[fCent_H+2], r_long[fCent_H+2],	qinv, qside, qout, qlong);
                    else							c_weight = fCoulomb->Coul_Wave(mpi, r_inv[fCent_H],		r_side[fCent_H],		r_out[fCent_H],		r_long[fCent_H],		qinv, qside, qout, qlong);
										for (int iqbin = 0; iqbin < fQbin; iqbin++) {
											fH1Qinv	[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv);
											fH1CQinv[1][fCent_H][fDphi_E][iqbin]	-> Fill(qinv, c_weight);
										}
									}

									if (conv_over(qout, qside, qlong) == 1) continue;

									for (int iqbin = 0; iqbin < fQbin; iqbin++) {
										fH3Q[1][fCent_H][fDphi_E][iqbin]	-> Fill(qout, qside, qlong);
									}

									Int_t	conv_qo		= conv_qosl(qout);	Int_t	conv_qs		= conv_qosl(qside);		Int_t	conv_ql		= conv_qosl(qlong);
									Int_t	conv_qo2	= conv_qosl2(qout);	Int_t	conv_qs2	= conv_qosl2(qside);	Int_t	conv_ql2	= conv_qosl2(qlong);
									Int_t	conv_qo3	= conv_qosl3(qout);	Int_t	conv_qs3	= conv_qosl3(qside);	Int_t	conv_ql3	= conv_qosl3(qlong);
									if (conv_if(conv_qo, conv_qs, conv_ql) == 0) {
										//fPrfConv	[1][fCent_H][fDphi_E]	-> Fill(conv_qo + 30*conv_qs + 900*conv_ql, qinv);
										fPrfConv	[1][fCent_H][fDphi_E][0]	-> Fill(conv_qo		+ 40*conv_qs	+ 1600*conv_ql, qinv);
										//fPrfConv	[1][fCent_H][fDphi_E][2]	-> Fill(conv_qo2	+ 20*conv_qs2	+ 400*conv_ql2, qinv);
									}
									if (conv_if2(conv_qo3, conv_qs3, conv_ql3) == 0) {
										fPrfConv	[1][fCent_H][fDphi_E][1]	-> Fill(conv_qo3	+ 60*conv_qs3	+ 3600*conv_ql3, qinv);
									}
								}	//j_loop
							}	//i_loop
						}	//ieve
						for (int ieve = 0; ieve < 1; ieve++) {
							pion_mpx	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpy	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpz	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mpt	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_me		[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mphi	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_meta	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_mclu	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_msha	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
							pion_max	[ftempCent][ftZvtx][ftEP3V0][ieve].clear();
						}	//ieve loop
						mix_flag[1][ftempCent][ftZvtx][ftEP3V0] = 0;
					}	//event pool==1 loop
				}	//Psi3 loop
			}	//Zvertex loop
		}	//Centrality loop
		PostData(1, fOutputList);
	}	//Last event of this file
}
