/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Macro to visualise rootified raw-data from AD.
//

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TStyle.h>
#include <TEveManager.h>

#include <AliESDAD.h>
#include <AliEveEventManager.h>
#include <AliEveADModule.h>
#else
class AliESDAD;
class AliEveADModule;
#endif

void ad_esd(Int_t maxCharge = 300, Bool_t showLegend = kFALSE)
{
    printf("\n*** ESD AD ***\n");
    
    gStyle->SetPalette(1, 0);
    
    AliESDAD *adESD = AliEveEventManager::Instance()->AssertESD()->GetADData();
    
    gEve->DisableRedraw();
    
    AliEveADModule* esdA = new AliEveADModule("AD_ESD_A", kTRUE, maxCharge, showLegend);
    esdA->LoadEsd(adESD);
    
    
    AliEveADModule* esdC = new AliEveADModule("AD_ESD_C", kFALSE, maxCharge, showLegend);
    esdC->LoadEsd(adESD);
    
    
    gEve->EnableRedraw();
}
