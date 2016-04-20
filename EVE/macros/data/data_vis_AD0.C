/// \file data_visAD0.C
/// \brief Data visualisation macro for AD
///
/// This macro is a starting point for visualisation of AD data. Available data types: Raw, Hits, ESD.
///
/// Parameters which can be modified in AliEVE config file:
/// AD.maxCharge.Raw - maximum expected charge in Raw data, used for colouring
/// AD.maxCharge.ESD - maximum expected charge in ESD, used for colouring
/// AD.showLegend - switch on/off legend's drawing
///
///
/// \author Matevz Tadel <Matevz.Tadel@cern.ch>, Univ. of California San Diego
/// \author Alja Mrak-Tadel <Alja.Mrak.Tadel@cern.ch>, Univ. of California San Diego
/// \author Michal Broz <Michal.Broz@cern.ch>, Czech Technical University
/// \author Jeremi Niedziela <jeremi.niedziela@cern.ch>, Warsaw University of Technology
/// \date Feb 23, 2016


/// Main function of the macro. Call data-specific function according to requested data type.
void data_vis_AD0(AliEveEventManager::EDataType type)
{

    switch (type)
    {
        case AliEveEventManager::kRaw:
            data_vis_AD0_raw();
            break;
        case AliEveEventManager::kHits:
            data_vis_AD0_hits();
            break;
        case AliEveEventManager::kESD:
            data_vis_AD0_esd();
            break;
        default:break;
    }
}

/// AD raw data visualisation
void data_vis_AD0_raw()
{
    AliEveEventManager *manager = AliEveEventManager::Instance();
    
    // get Raw Reader from the Event Manager
    AliRawReader *reader = manager->GetRawReader();
    // always check if an object exists!
    if(!reader)
    {
        Error("data_vis_AD0_raw", "No Raw Reader");
        return;
    }
    reader->Reset();
    
    // read settings from AliEVE's config file
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    Int_t maxCharge = settings.GetValue("AD.maxCharge.Raw",1023);
    Bool_t showLegend = settings.GetValue("AD.showLegend",false);
    
    gStyle->SetPalette(1, 0);
    
    // Ask Event Manager to block redrawing
    manager->DisableRedraw();
    
    // create detector-specific visualisation object
    AliEveADModule* rawA = new AliEveADModule("AD_RAW_A", kTRUE, maxCharge, showLegend);
    AliEveADModule* rawC = new AliEveADModule("AD_RAW_C", kFALSE, maxCharge, showLegend);
    
    // load data from the Raw Reader
    rawA->LoadRaw(reader);
    rawC->LoadRaw(reader);
    
    // Tell Event Manager to draw modules
    manager->AddElement(rawA);
    manager->AddElement(rawC);
    manager->Redraw3D();
    manager->EnableRedraw();
}

/// AD hits visualisation
void data_vis_AD0_hits()
{
    AliEveEventManager *maanger = AliEveEventManager::Instance();
    
    // ask Event Manager for Run Loader
    AliRunLoader* runLoader =  manager->GetRunLoader();

    // check if Run Loader exists
    if(!runLoader)
    {
        Error("data_vis_AD0_hits", "No Run Loader");
        return;
    }
    // load hits from the Run Loader
    runLoader->LoadHits("AD");
    
    // Get hits tree from Run Loader
    TTree* hitsTree = runLoader->GetTreeH("AD", false);

    // check if tree exists
    if(!hitsTree)
    {
        Error("data_vis_AD0_hits", "Hits tree not found");
        return;
    }
    
    // create visualisation object (set of points in this case) and fill it with data from a tree
    // in general it would be better to put this part in a dedicated class!
    const char  *selection = "";
    TEvePointSet* points = new TEvePointSet(Form("AD Hits '%s'", selection));
    TEvePointSelector pointSelector(hitsTree, points, "AD.fX:AD.fY:AD.fZ", selection);
    pointSelector.Select();
    
    if(points->Size() == 0)
    {
        Warning("data_vis_AD0_hits", "No hits match '%s'", selection);
        delete points;
        return;
    }
    
    points->SetName(Form("AD Hits"));
    const TString viz_tag("SIM Hits AD");
    points->ApplyVizTag(viz_tag, "Hits");
    
    points->SetTitle(Form("N=%d", points->Size()));
    points->SetMarkerSize(.5);
    points->SetMarkerColor(2);
    
    // register visualisation object in Event Manager and ask to redraw
    manager->AddElement(points, 0);
    manager->Redraw3D();
}

/// AD data visualisation from ESD
void data_vis_AD0_esd()
{
    AliEveEventManager *manager = AliEveEventManager::Instance();
    
    // read settings from AliEVE's config file
    TEnv settings;
    AliEveInit::GetConfig(&settings);
    Int_t maxCharge = settings.GetValue("AD.maxCharge.Raw",300);
    Bool_t showLegend = settings.GetValue("AD.showLegend",false);
    
    gStyle->SetPalette(1, 0);
    
    // get ESD from the Event Manager
    AliESDEvent *esdEvent = manager->GetESD();
    if(!esdEvent)
    {
        cout<<"data_vis_AD0_esd -- No ESD event"<<endl;
        return;
    }
    AliESDAD *adESD = esdEvent->GetADData();
    
    // check if ESD exists
    if(!adESD)
    {
        Error("data_vis_AD0_esd", "No AD ESD");
        return;
    }
    
    // Ask Event Manager to block redrawing
    manager->DisableRedraw();
    
    // create detector-specific visualisation object
    AliEveADModule* esdA = new AliEveADModule("AD_ESD_A", kTRUE, maxCharge, showLegend);
    AliEveADModule* esdC = new AliEveADModule("AD_ESD_C", kFALSE, maxCharge, showLegend);

    // load data from ESD
    esdA->LoadEsd(adESD);
    esdC->LoadEsd(adESD);
    
    // Tell Event Manager to draw modules
    manager->AddElement(esdA);
    manager->AddElement(esdC);
    manager->Redraw3D();
    manager->EnableRedraw();
}