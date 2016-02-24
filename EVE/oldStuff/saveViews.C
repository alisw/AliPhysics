/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
// Author: Mihai Niculescu 2013


/*
 * This script creates a collage containing all OpenGL views from a running AliEve
 *
 * Given Collage size (width, height), the size for all OpenGL
 * views are computed using the same aspect ratio (width/height) as the main 3D View
 */

#include <TASImage.h>
#include <TGLViewer.h>
#include <TEveViewer.h>
#include <TEveManager.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTimeStamp.h>
#include <TEveElement.h>
#include <TIterator.h>
#include <TList.h>
#include <TROOT.h>

#include <STEER/ESD/AliESDEvent.h>
#include <STEER/ESD/AliESDRun.h>
#include <STEER/STEER/AliGRPObject.h>
#include <STEER/CDB/AliCDBEntry.h>
#include <STEER/CDB/AliCDBManager.h>
#include <RAW/RAWDatarec/AliRawReader.h>
#include <RAW/RAWDatabase/AliRawEventHeaderBase.h>
#include <EVE/EveBase/AliEveEventManager.h>
#include <STEER/STEERBase/AliDAQ.h>
#include <MONITOR/alionlinereco/AliOnlineReconstructionUtil.h>


#include <TEnv.h>
#include <TString.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TSystem.h>

#include <iostream>

TString getEventInfo();

/***********Save all OpenGL views into one picture
	compositeImgFileName - save final image to this file
 showLiveBar - whether to show the LIVE bar, useful when not online (using offline)
	width - of the collage image
	height -of the collage image
 */
void saveViews(const char* compositeImgFileName="views.png", Bool_t showLiveBar=kTRUE, Int_t width = 1440, Int_t height= 900)
{
    Info("saveViews.C", "saving views to [%s]", compositeImgFileName);
    
    Int_t heightInfoBar = 65; // hold height of the Information bar
    
    TASImage* compositeImg=0; // this holds the final image
    TASImage* tempImg=0; // temporary used for loading images
    
    TEveViewerList* viewers = gEve->GetViewers();
    Int_t Nviewers = viewers->NumChildren()-2; // remark: 3D view is counted twice
    
    compositeImg = new TASImage(width, height);
    
    // 3D View size
    Int_t width3DView = TMath::FloorNint((float)Nviewers*width/(float)(Nviewers+1)); // the width of the 3D view
    Int_t height3DView= height-heightInfoBar; // the height of the 3D view
    Float_t aspectRatio = (float)width3DView/(float)height3DView; // 3D View aspect ratio
    
    // Children View Size
    Int_t heightChildView = TMath::FloorNint((float)height3DView/(float)Nviewers);
    Int_t widthChildView  = TMath::FloorNint(aspectRatio*heightChildView); // has the same aspect ratio as the 3D view
    
    int index=0; // iteration counter
    int x = width3DView; // x position of the child view
    int y = 0;// y position of the child view
    TString viewFilename; // save view to this file
    
    for(TEveElement::List_i i = (++viewers->BeginChildren()); i != viewers->EndChildren(); i++)
    { // NB: this skips the first children (first 3D View)
        TEveViewer* view = ((TEveViewer*)*i);
        viewFilename = Form("view-%d.png", index);
        
        // Save OpenGL views in files
        if(index==0){
            view->GetGLViewer()->SavePictureUsingBB(viewFilename);//, width3DView, height3DView);
        }
        else {
            view->GetGLViewer()->SavePictureUsingBB(viewFilename);//, widthChildView, heightChildView);
        }
        
        tempImg = new TASImage(viewFilename);
        
        // copy view image in the composite image
        if(index==0){
            if(tempImg->GetWidth()>width3DView)
            {
                tempImg->Crop((tempImg->GetWidth()-width3DView)/2,
                               0,width3DView,height3DView);
            }
            if(tempImg->GetHeight()>height3DView)
            {
                tempImg->Crop(0,(tempImg->GetHeight()-height3DView)/2,
                              width3DView,height3DView);
            }
            
            tempImg->CopyArea(compositeImg, 0,0, width3DView, height3DView);
        }
        else {
            tempImg->Crop((tempImg->GetWidth()-widthChildView)/2,
                          (tempImg->GetHeight()-heightChildView)/2,
                          widthChildView,heightChildView);
            tempImg->CopyArea(compositeImg, 0,0, widthChildView, heightChildView, x,y);
            
            // draw a border around child views
            compositeImg->DrawRectangle(x,y, widthChildView, heightChildView, "#C0C0C0");
        }
        
        /*
         final touches inside loop
         */
        delete tempImg;
        if(index>0) // skip 3D View
            y+=heightChildView;
        
        index++;
    }
    
    // Create a glow (bloom) effect
    tempImg = (TASImage*)compositeImg->Clone("tempImg");
    tempImg->Blur(10.0,10.0);
    compositeImg->Merge(tempImg, "lighten");
    delete tempImg; tempImg = 0;
    
    
    // show LIVE bar
    if(showLiveBar){
        TTimeStamp ts;
        TString tNow = ts.AsString("s"); // display date & time
        
        compositeImg->Gradient( 90, "#EAEAEA #D2D2D2 #FFFFFF", 0, 75, 0, 239, 95);
        compositeImg->Gradient( 90, "#D6D6D6 #242424 #000000", 0, 155, 60, 152, 26);
        compositeImg->BeginPaint();
        compositeImg->DrawRectangle(50,0, 264, 94);
        compositeImg->DrawText(162, 6, "LIVE", 70, "#FF2D00", "FreeSansBold.otf");
        compositeImg->DrawText(162, 65, tNow, 16, "#FFFFFF", "arial.ttf");
        compositeImg->EndPaint();
        //include ALICE Logo
//        tempImg = new TASImage(Form("%s/../src/picts/alice_logo_2009_blue_80x80.png", gSystem->Getenv("ALICE_ROOT")) );
        
        
        tempImg = new TASImage("alice_logo.png");
        
        tempImg->Scale(64,87);
        //tempImg->CopyArea(compositeImg, 0,0, 236, 319, 59, 4);
        compositeImg->Merge(tempImg, "alphablend", 82, 4);
        delete tempImg; tempImg = 0;
    }
    
    //---------------------------------------------------
    //              show Information bar
    //---------------------------------------------------
    
    // read crededentials to logbook from config file:
    TEnv settings;
    settings.ReadFile(AliOnlineReconstructionUtil::GetPathToServerConf(), kEnvUser);
    const char *dbHost = settings.GetValue("logbook.host", "");
    Int_t   dbPort =  settings.GetValue("logbook.port", 0);
    const char *dbName =  settings.GetValue("logbook.db", "");
    const char *user =  settings.GetValue("logbook.user", "");
    const char *password = settings.GetValue("logbook.pass", "");
    
//    cout<<dbHost<<"\t"<<dbPort<<"\t"<<dbName<<"\t"<<user<<"\t"<<password<<endl;
    
    // get entry from logbook:
    cout<<"creating sql server...";
    TSQLServer* server = TSQLServer::Connect(Form("mysql://%s:%d/%s", dbHost, dbPort, dbName), user, password);
    cout<<"created"<<endl;
    
    AliESDEvent* esd =  AliEveEventManager::Instance()->AssertESD();
    int run = esd->GetRunNumber();
    TString sqlQuery;
    sqlQuery.Form("SELECT * FROM logbook_trigger_clusters WHERE run = %d",run);
    
    TSQLResult* result = server->Query(sqlQuery);
    if (!result){
        Printf("ERROR: Can't execute query <%s>!", sqlQuery.Data());
        return;
    }
    
    if (result->GetRowCount() == 0){
        Printf("ERROR: Run %d not found", run);
        delete result;
        return;
    }

    // read values from logbook entry:
    vector<int> cluster;
    vector<ULong64_t> inputDetectorMask;
    vector<ULong64_t> triggerClassMask0;
    vector<ULong64_t> triggerClassMask1;
    
    TSQLRow* row;
    int i=0;
    while(row=result->Next())
    {
        cluster.push_back(atoi(row->GetField(2)));
        inputDetectorMask.push_back(atol(row->GetField(4)));
        triggerClassMask0.push_back(atol(row->GetField(5)));
        triggerClassMask1.push_back(atol(row->GetField(6)));
        i++;
    }

    ULong64_t mask = 1;
    
    // 3 lines of trigger classes
    TString triggerClasses1 = "";
    TString triggerClasses2 = "";
    TString triggerClasses3 = "";
    int sw=0;
    
    vector<TString> clustersDescription;
    TString detectorsInCluster;

    ULong64_t triggerMask = esd->GetTriggerMask();
    ULong64_t triggerMaskNext50 = esd->GetTriggerMaskNext50();
    
    for(int k=0;k<i;k++)//loop over all clusters in run
    {
        detectorsInCluster="";
        
        // get trigger classes for given cluster
        mask=1;
        for(int i=0;i<50;i++)
        {
            if(mask & triggerMask)
            {
                switch (sw) {
                    case 0:
                        triggerClasses1+=esd->GetESDRun()->GetTriggerClass(i);
                        triggerClasses1+=Form("(%d)",cluster[k]);
                        triggerClasses1+="   ";
                        sw=1;
                        break;
                    case 1:
                        triggerClasses2+=esd->GetESDRun()->GetTriggerClass(i);
                        triggerClasses2+=Form("(%d)",cluster[k]);
                        triggerClasses2+="   ";
                        sw=2;
                        break;
                    case 2:
                        triggerClasses3+=esd->GetESDRun()->GetTriggerClass(i);
                        triggerClasses3+=Form("(%d)",cluster[k]);
                        triggerClasses3+="   ";
                        sw=0;
                        break;
                    default:break;
                }
            }
            if(mask & triggerMaskNext50)
            {
                switch (sw) {
                    case 0:
                        triggerClasses1+=esd->GetESDRun()->GetTriggerClass(i+50);
                        triggerClasses1+=Form("(%d)",cluster[k]);
                        triggerClasses1+="   ";
                        sw=1;
                        break;
                    case 1:
                        triggerClasses2+=esd->GetESDRun()->GetTriggerClass(i+50);
                        triggerClasses2+=Form("(%d)",cluster[k]);
                        triggerClasses2+="   ";
                        sw=2;
                        break;
                    case 2:
                        triggerClasses3+=esd->GetESDRun()->GetTriggerClass(i+50);
                        triggerClasses3+=Form("(%d)",cluster[k]);
                        triggerClasses3+="   ";
                        sw=0;
                        break;
                    default:break;
                }
            }
            mask = mask<<1;
        }
        
        // get trigger detectors in given cluster
        mask=1;
        for(int ii=0;ii<20;ii++)
        {
            if(inputDetectorMask[k] & mask)
            {
                detectorsInCluster+=AliDAQ::fgkDetectorName[ii];
                detectorsInCluster+=", ";
            }
            
            mask=mask<<1;
        }
        clustersDescription.push_back(detectorsInCluster);
    }

    // put trigger classes in blue bar on the bottom:
    TString stringInfo;
    stringInfo = getEventInfo();
    compositeImg->Gradient( 90, "#1B58BF #1D5CDF #0194FF", 0, 0, height-1.33*heightInfoBar, width, heightInfoBar);
    compositeImg->BeginPaint();
    compositeImg->DrawText(10, height-1.33*heightInfoBar+15, stringInfo, 28, "#FFFFFF", "FreeSansBold.otf");
    compositeImg->DrawText(750, height-1.33*heightInfoBar+4,triggerClasses1, 16, "#FFFFFF", "FreeSansBold.otf");
    compositeImg->DrawText(750, height-1.33*heightInfoBar+24,triggerClasses2, 16, "#FFFFFF", "FreeSansBold.otf");
    compositeImg->DrawText(750, height-1.33*heightInfoBar+44,triggerClasses3, 16, "#FFFFFF", "FreeSansBold.otf");
    compositeImg->EndPaint();

    
    // put clusters description in green bar on the bottom:
    compositeImg->Gradient( 90, "#1BDD1B #1DDD1D #01DD01", 0, 0, height-0.33*heightInfoBar, width, 0.33*heightInfoBar);
    compositeImg->BeginPaint();
    
    TString cd="";
    i=0;
    for (TString desc : clustersDescription) {
        cd+="Cluster ";
        cd+=cluster[i];
        cd+=":";
        cd+=desc;
        cd+="\t";
//        compositeImg->DrawText(10, height-heightInfoBar+4+20*i,cd, 16, "#FFFFFF", "FreeSansBold.otf");
        i++;
    }
    compositeImg->DrawText(10, height-0.33*heightInfoBar+2,cd, 16, "#000000", "FreeSansBold.otf");
    compositeImg->EndPaint();
    
    // write composite image to disk
    compositeImg->CopyArea(compositeImg, 0,0, width, height);
    compositeImg->WriteImage(compositeImgFileName);
    
    delete compositeImg;
    
    return;
}

// This function retrieves a string containing some information regarding the current event
TString getEventInfo()
{
    // For general public please show as less or technical information as possible
    
    TString rawInfo, esdInfo;
    
    if (!AliEveEventManager::HasRawReader())
    {
        rawInfo = "";
    }
    else
    {
        AliRawReader* rawReader = AliEveEventManager::AssertRawReader();
        if(!rawReader) return "";
        rawInfo.Form("Run: %d  Event#: %d (%s)",
                     rawReader->GetRunNumber(),
                     AliEveEventManager::CurrentEventId(),
                     AliRawEventHeaderBase::GetTypeName(rawReader->GetType())
                     );
        
        return rawInfo;
    }
    
    if (!AliEveEventManager::HasESD())
    {
        esdInfo = "";
    }
    else
    {
        AliESDEvent* esd =  AliEveEventManager::Instance()->AssertESD();
       
        
        
        esdInfo.Form("Colliding: %s Run: %d  Event: %d (%s)",
                     esd->GetESDRun()->GetBeamType(),
                     esd->GetRunNumber(),
                     esd->GetEventNumberInFile(),
                     AliRawEventHeaderBase::GetTypeName(esd->GetEventType())
                     );
    }
    
    return esdInfo;
}
