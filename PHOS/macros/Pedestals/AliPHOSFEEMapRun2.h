#ifndef AliPHOSFEEMapRun2_h
#define AliPHOSFEEMapRun2_h

//This class is developed for calibration of PHOS in Run-2 (2015-2017).
//Author : Daiki Sekihata (Hiroshima University), 12.November.2015
//daiki.sekihata@cern.ch

#include <TObject.h>

class AliPHOSFEEMapRun2 : public TObject {
	public:
		AliPHOSFEEMapRun2();//default constructor
		AliPHOSFEEMapRun2(Int_t module,Int_t cellx,Int_t cellz);

		~AliPHOSFEEMapRun2();//destructor

	public:
		Int_t CellToSRUID(Int_t cellx);//this returns SRU ID [0-3].
		Int_t CellToFEEID(Int_t cellz);//this returns FEE ID [1-14] on branch 0, [21-34] on branch 1
		Int_t CellToALTRO(Int_t cellx,Int_t cellz);//this returns ALTRO ID [2,3,0,4] on 1 FEE card.
		Int_t CellToCSPID(Int_t module,Int_t cellx,Int_t cellz);//this returns CSP ID [0-15] on 1 FEE card.
		Int_t CSPToHVID(Int_t csp);//this returns HVBIAS register [0x60-0x7f] in decimal.
		Int_t CSPToALTROChannel(Int_t csp,TString gain);//this returns ALTRO Hi/Lo channel [0-15] for pedestal calculation.

    Int_t CellToTRUID(Int_t module,Int_t cellx,Int_t cellz);//this returns TRU ID [0-27];
    Int_t CellToTRUChannel(Int_t cellx,Int_t cellz);//this returns TRU channel [0-111];
    void TRUHWToCellID(Int_t ddl, Int_t hwaddress, Int_t &cellx, Int_t &cellz);// this returns minimum cell ID in 1 fastOR = 2x2 crystals.

		void  GetElectronicsMap(Int_t &sru, Int_t &fee, Int_t &altro, Int_t &csp, Int_t &hvch, Int_t &altlg, Int_t &althg){
			sru   = fSRUID;
			fee   = fFEEID;
			altro = fALTRO;
			csp   = fCSPID;
			hvch  = fHVID;
      altlg = fALTCH_LG;
      althg = fALTCH_HG;
		}

    void Print(Option_t *option="") const;

	private:
		Int_t fSRUID;
		Int_t fFEEID;
		Int_t fALTRO;
		Int_t fCSPID;
		Int_t fHVID;
		Int_t fALTCH_LG;
		Int_t fALTCH_HG;
		Int_t fTRUID;
		Int_t fTRUCH;

	ClassDef(AliPHOSFEEMapRun2,5);
};
#endif
