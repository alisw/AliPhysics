#include <Riostream.h>
#include "TVirtualMCApplication.h"
#include "TVirtualMCStack.h"
#include "TFluka.h"
#include "TFlukaCodes.h"

// Fluka includes
#include "Fdimpar.h"  //(DIMPAR) 
#include "Fdblprc.h"  //(DBLPRC) 
#include "Ftrackr.h"  //(TRACKR) 
#include "Fopphst.h"  //(OPPHST) 
#include "Fflkstk.h"  //(FLKSTK) 
#include "Fltclcm.h"  //(LTCLCM) 
#include "Fpaprop.h"  //(PAPROP) 
#include "Falldlt.h"  //(ALLDLT) 
#include "Fdpdxcm.h"  //(DPDXCM) 
#include "Fflkmat.h"  //(FLKMAT) 

#include "TGeoManager.h"

#ifndef WIN32
# define mgdraw mgdraw_
#else
# define mgdraw MGDRAW
#endif

extern "C" {
void mgdraw(Int_t& icode, Int_t& mreg)
{
    TFluka* fluka =  (TFluka*) gMC;
    if (mreg == fluka->GetDummyRegion()) return;
//
//  Make sure that stack has currrent track Id
//
    Int_t trackId = -1;
    TVirtualMCStack* cppstack = fluka->GetStack();
    
    if (TRACKR.jtrack == -1) {
        trackId = OPPHST.louopp[OPPHST.lstopp];
        if (trackId == 0) {
            trackId = FLKSTK.ispark[FLKSTK.npflka][mkbmx2-1];
        }
    } else {
        trackId = TRACKR.ispusr[mkbmx2-1];
    }
    
    Int_t verbosityLevel = fluka->GetVerbosityLevel();

    if (TRACKR.jtrack < -6) {
       // from -7 to -12 = "heavy" fragment
       // assing parent id
       // id < -6 was skipped in stuprf =>   if (kpart < -6) return;
       if (verbosityLevel >= 3) {
          cout << "mgdraw: (heavy fragment) jtrack < -6 =" << TRACKR.jtrack
               << " assign parent pdg=" << fluka->PDGFromId(TRACKR.ispusr[mkbmx2 - 3]) << endl;
       }
       TRACKR.jtrack = TRACKR.ispusr[mkbmx2 - 3];
    }
    
    cppstack->SetCurrentTrack(trackId);
//
//    
    Int_t mlttc = TRACKR.lt1trk; // LTCLCM.mlatm1;
    fluka->SetMreg(mreg, mlttc);
    fluka->SetIcode((FlukaProcessCode_t) icode);
    fluka->SetCaller(kMGDRAW);

    Int_t nodeId;
    Int_t volId   = fluka->CurrentVolID(nodeId);
    Int_t crtlttc = gGeoManager->GetCurrentNodeId()+1;

    // check region lattice consistency (debug Ernesto)
    // *****************************************************
    if( mreg != volId  && !gGeoManager->IsOutside() ) {
       cout << "  mgdraw:   track=" << trackId << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
            << " icode=" << icode << " gNstep=" << fluka->GetNstep() << endl
            << "               fluka   mreg=" << mreg << " mlttc=" << mlttc << endl
            << "               TGeo   volId=" << volId << " crtlttc=" << crtlttc << endl
            << "     common TRACKR   lt1trk=" << TRACKR.lt1trk << " lt2trk=" << TRACKR.lt2trk << endl
            << "     common LTCLCM   newlat=" << LTCLCM.newlat << " mlatld=" <<  LTCLCM.mlatld << endl
            << "                     mlatm1=" << LTCLCM.mlatm1 << " mltsen=" <<  LTCLCM.mltsen << endl
            << "                     mltsm1=" << LTCLCM.mltsm1 << " mlattc=" << LTCLCM.mlattc << endl;
        if( mlttc == crtlttc ) cout << "   *************************************************************" << endl;
    }
    // *****************************************************

    
    Int_t med   = FLKMAT.medium[mreg - 1];  // Medium
    Int_t msd   = DPDXCM.msdpdx[med  - 1];  // Iionisation model

    if (!TRACKR.ispusr[mkbmx2 - 2]) {
	if (verbosityLevel >= 3) {
	    cout << endl << "mgdraw: energy deposition for:" << trackId
		 << " icode=" << icode
		 << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
		 << " flukaid="<< TRACKR.jtrack
		 << " mreg=" << mreg
		 << " np  =" << ALLDLT.nalldl
		 << endl;
	         
	}

	if (msd > 0) {
//
// Primary ionsiations
	    Int_t nprim = ALLDLT.nalldl;
// Protection against nprim > mxalld
	    if (nprim >= mxalld) {
		nprim = mxalld;
		Warning("mgdraw", "nprim > mxalld, nprim: %6d  pdg: %6d mreg %6d p %13.3f step %13.3f\n", 
		       ALLDLT.nalldl, 
		       fluka->PDGFromId(TRACKR.jtrack), 
		       mreg, 
		       TRACKR.ptrack, 
		       TRACKR.ctrack);
		
	    }
// Multiple steps for nprim > 0
	    if (nprim > 0) {
		for (Int_t i = 0; i < nprim; i++) {
		    fluka->SetCurrentPrimaryElectronIndex(i);
		    (TVirtualMCApplication::Instance())->Stepping();
		    if (i == 0) fluka->SetTrackIsNew(kFALSE);
		}
	    } else {
		// No primary electron ionisation
		// Call Stepping anyway but flag nprim = 0 as index = -2
		fluka->SetCurrentPrimaryElectronIndex(-2);
		(TVirtualMCApplication::Instance())->Stepping();
	    }
	    // Reset the index
	    fluka->SetCurrentPrimaryElectronIndex(-1);
	} else {
	    // Single step
	    (TVirtualMCApplication::Instance())->Stepping();
	    fluka->SetTrackIsNew(kFALSE);
	}
	
    } else {
        //
        // Tracking is being resumed after secondary tracking
        //
        if (verbosityLevel >= 3) {
            cout << endl << "mgdraw: resuming Stepping(): " << trackId << endl;
        }

        fluka->SetTrackIsNew(kTRUE);
        fluka->SetCaller(kMGResumedTrack);
        (TVirtualMCApplication::Instance())->Stepping();

        // Reset flag and stored values
        TRACKR.ispusr[mkbmx2 - 2] = 0;
        for (Int_t i = 0; i < 9; i++) TRACKR.spausr[i] = -1.;


        if (verbosityLevel >= 3) {
            cout << endl << " !!! I am in mgdraw - first Stepping() after resume: " << icode << endl;
            cout << " Track= " << trackId << " region = " << mreg << endl;
        }

        fluka->SetTrackIsNew(kFALSE);
        fluka->SetCaller(kMGDRAW);
	if (msd > 0) fluka->SetCurrentPrimaryElectronIndex(-2);
        (TVirtualMCApplication::Instance())->Stepping();
	if (msd > 0) fluka->SetCurrentPrimaryElectronIndex(-1);
    }
} // end of mgdraw
} // end of extern "C"

