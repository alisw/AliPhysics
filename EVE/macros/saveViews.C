/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
 // Author: Mihai Niculescu 2013
 
 
 /*
 * This script creates a single picture containing all OpenGL views from a running AliEve
 *
 * Considering a given size (width, height) for the final picture, the size for all OpenGL 
 * views are computed using the same aspect ratio (width/height) as the main 3D View
 */

#include <TASImage.h>
#include <TGLViewer.h>
#include <TEveViewer.h>
#include <TEveManager.h>
#include <TMath.h>
#include <TSystem.h>
#include <TTimeStamp.h>

#include <STEER/ESD/AliESDEvent.h>
#include <STEER/ESD/AliESDRun.h>
#include <STEER/STEER/AliGRPObject.h>
#include <STEER/CDB/AliCDBEntry.h>
#include <STEER/CDB/AliCDBManager.h>
#include <RAW/AliRawReader.h>
#include <RAW/AliRawEventHeaderBase.h>
#include <EVE/EveBase/AliEveEventManager.h>

TString getEventInfo();

/* composite picture final size
	compositeImgFileName - save final image to this file	
	width - of the final image
	height -of the final image
*/
void saveViews(const char* compositeImgFileName="views.png", Int_t width = 1440, Int_t height= 900)
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
	for(TEveElement::List_i i = (++viewers->BeginChildren()); i != viewers->EndChildren(); i++){ // NB: this skips the first children (first 3D View)
		TEveViewer* view = ((TEveViewer*)*i);
		viewFilename = Form("view-%d.png", index);

    // Save OpenGL views in files
    if(index==0){
			view->GetGLViewer()->SavePictureUsingFBO(viewFilename, width3DView, height3DView);
		}
		else {
			view->GetGLViewer()->SavePictureUsingFBO(viewFilename, widthChildView, heightChildView);
		}
		
		tempImg = new TASImage(viewFilename);
		
		// copy view image in the composite image
		if(index==0){
			tempImg->CopyArea(compositeImg, 0,0, width3DView, height3DView);
		}
		else {
			tempImg->CopyArea(compositeImg, 0,0, widthChildView, heightChildView, x,y);
		    
    // draw a border around child views
    compositeImg->DrawRectangle(x,y, widthChildView, heightChildView);
    //compositeImg->FillRectangle("#CFEDF0", x+1,y+1, widthChildView-1, 40);
    //compositeImg->DrawText(x+5,y+10, view->GetName(), 28, "#000000", "verdana.ttf", TImage::kOutlineBelow);
    
    }
    
    /*
     final touches inside loop
    */
    delete tempImg;
    if(index>0) // skip 3D View
    	y+=heightChildView;
    	
    index++;
   }

 
 // show LIVE bar
	TTimeStamp ts;
	TString tNow = ts.AsString("s"); // display date & time
 
	compositeImg->Gradient( 90, "#EAEAEA #D2D2D2 #FFFFFF", 0, 30, 0, 234, 75);
	compositeImg->Gradient( 90, "#D6D6D6 #242424 #000000", 0, 35, 40, 224, 30);
	compositeImg->BeginPaint();
	compositeImg->DrawRectangle(30,0, 234, 75);
	compositeImg->DrawText(35, 4, "ALICE", 40, "#000000", "FreeSans.otf");
	compositeImg->DrawText(162, 6, "LIVE", 40, "#FF2D00", "FreeSansBold.otf");
	compositeImg->DrawText(59, 48, tNow, 20, "#FFFFFF", "arial.ttf");
	compositeImg->EndPaint();
	
	// show Information bar
	TString stringInfo;
	stringInfo = getEventInfo();
	compositeImg->Gradient( 90, "#1B58BF #1D5CDF #0194FF", 0, 0, height-heightInfoBar, width, heightInfoBar);
	compositeImg->BeginPaint();
	compositeImg->DrawText(10, height-heightInfoBar+15, stringInfo, 28, "#FFFFFF", "FreeSansBold.otf");
	compositeImg->EndPaint();
		
	
	// write composite image to disk
	compositeImg->CopyArea(compositeImg, 0,0, width, height);
	compositeImg->WriteImage(compositeImgFileName);
	
	delete compositeImg;
	
	return;
}

TString getEventInfo()
{
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
		AliESDEvent* esd =  AliEveEventManager::AssertESD();

		esdInfo.Form("Colliding: %s Run: %d  Event: %d (%s)",
		esd->GetESDRun()->GetBeamType(),
		esd->GetRunNumber(),
		AliEveEventManager::CurrentEventId(),
		AliRawEventHeaderBase::GetTypeName(esd->GetEventType())
		);
  }

  return esdInfo;
}
