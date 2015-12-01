/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

// Author: Pawel Debski 2010

#include <AliEveMomentumVectors.h>

#include <TPolyMarker3D.h>
#include <TString.h>
#include <TEveLine.h>
#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveElement.h>
#include <TEveUtil.h>

#include <AliESDtrack.h>
#include <AliEveTrack.h>
#include <AliEveMultiView.h>

#include <iostream>

using namespace std;

AliEveMomentumVectors::AliEveMomentumVectors()
{
    
    TString str1;
    TString str2;
    
    Bool_t draw = kFALSE;
    
    TEveElement::List_i i = gEve->GetEventScene()->FirstChild()->BeginChildren();
    TEveElement::List_i j = gEve->GetEventScene()->FirstChild()->EndChildren();
    TEveElement::List_i k;
    
    TEveElementList* momentumVectorList = new TEveElementList("Momentum Vectors");
    
    Double_t maxMomentum = 0;
    Double_t vectorLength = 600.0;
    
    Double_t x1 = 0;
    Double_t y1 = 0;
    Double_t z1 = 0;
    
    Double_t x2 = 0;
    Double_t y2 = 0;
    Double_t z2 = 0;
    
    //==============================================
    // find highest momentum (to normalize)
    //==============================================
    
    for(k = i; k != j; k++)
    {
        TEveElement* element = (TEveElement*) *k;
        
        str1 = element->GetElementName();
        
        if(str1.Contains("Tracks") || str1.Contains("tracks"))
        {
            
            TEveElement::List_i m = element->BeginChildren();
            TEveElement::List_i n = element->EndChildren();
            TEveElement::List_i l;
            
            for(l = m; l != n; l++)
            {
                TEveElement* trackType = (TEveElement*) *l;
                str2 = trackType->GetElementName();
                
                if(str2.Contains("Sigma < 3"))
                {
                    if(trackType->HasChildren())
                    {
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                                maxMomentum = trackSingle1->GetESDTrack()->P();
                            
                        }
                    }
                }
                
                
                if(str2.Contains("3 < Sigma < 5"))
                {
                    
                    if(trackType->HasChildren())
                    {
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                                maxMomentum = trackSingle1->GetESDTrack()->P();
                            
                        }
                    }
                }
                
                if(str2.Contains("5 < Sigma"))
                {
                    
                    if(trackType->HasChildren())
                    {
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            if(trackSingle1->GetESDTrack()->P() > maxMomentum)
                                maxMomentum = trackSingle1->GetESDTrack()->P();
                            
                            
                        }
                    }
                }
            }
        }
    }
    
    //==============================================
    // clean the display
    //==============================================
    /*
     if(!drawWithTracks)
     
     for(k = i; k != j; k++)
     {
     TEveElement* element = (TEveElement*) *k;
     
     str1 = element->GetElementName();
     
     element->SetRnrSelf(kFALSE);
     
     if(element->HasChildren())
     element->SetRnrChildren(kFALSE);
     
     }
     
     }
     */
    //==============================================
    // draw momentum vectors
    //==============================================
    
    if(maxMomentum)
        vectorLength = vectorLength/maxMomentum;
    //      vectorLength = vectorLength/TMath::Log(maxMomentum);
    
    for(k = i; k != j; k++)
    {
        TEveElement* element = (TEveElement*) *k;
        
        str1 = element->GetElementName();
        
        if(str1.Contains("Tracks") || str1.Contains("tracks"))
        {
            
            TEveElement::List_i m = element->BeginChildren();
            TEveElement::List_i n = element->EndChildren();
            TEveElement::List_i l;
            
            for(l = m; l != n; l++)
            {
                TEveElement* trackType = (TEveElement*) *l;
                str2 = trackType->GetElementName();
                
                trackType->SetRnrSelf(kFALSE);
                
                if(trackType->HasChildren())
                    trackType->SetRnrChildren(kFALSE);
                
                if(str2.Contains("Sigma < 3"))
                {
                    
                    if(trackType->HasChildren())
                    {
                        
                        TEveElementList* momentumVectorList1 = new TEveElementList("sigma < 3");
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));
                            
                            x1 = trackSingle1->GetESDTrack()->Xv();
                            y1 = trackSingle1->GetESDTrack()->Yv();
                            z1 = trackSingle1->GetESDTrack()->Zv();
                            
                            momentumVector->SetPoint(0, x1, y1, z1);
                            
                            x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                            y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                            z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                            
                            momentumVector->SetPoint(1, x2, y2, z2);
                            
                            /*
                             if(trackSingle1->GetESDTrack()->Charge() == -1)
                             momentumVector->SetLineColor(kGreen);
                             else
                             momentumVector->SetLineColor(kRed);
                             */
                            momentumVector->SetLineColor(kRed);
                            
                            momentumVector->SetLineWidth(1);
                            momentumVector->SetLineStyle(0);
                            momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));
                            
                            momentumVectorList1->AddElement(momentumVector);
                            
                        }
                        
                        //                  gEve->AddElement(momentumVectorList1);
                        momentumVectorList->AddElement(momentumVectorList1);
                        
                        draw = kTRUE;
                        
                    }
                }
                
                
                if(str2.Contains("3 < Sigma < 5"))
                {
                    
                    if(trackType->HasChildren())
                    {
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        TEveElementList* momentumVectorList2 = new TEveElementList("3 < sigma < 5");
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));
                            
                            x1 = trackSingle1->GetESDTrack()->Xv();
                            y1 = trackSingle1->GetESDTrack()->Yv();
                            z1 = trackSingle1->GetESDTrack()->Zv();
                            
                            momentumVector->SetPoint(0, x1, y1, z1);
                            
                            x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                            y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                            z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                            
                            momentumVector->SetPoint(1, x2, y2, z2);
                            /*
                             if(trackSingle1->GetESDTrack()->Charge() == -1)
                             momentumVector->SetLineColor(kGreen+2);
                             else
                             momentumVector->SetLineColor(kRed+2);
                             */
                            momentumVector->SetLineColor(kRed+2);
                            
                            momentumVector->SetLineWidth(1);
                            momentumVector->SetLineStyle(0);
                            momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));
                            
                            momentumVectorList2->AddElement(momentumVector);
                            
                        }
                        
                        //                  gEve->AddElement(momentumVectorList2);
                        momentumVectorList->AddElement(momentumVectorList2);
                        
                        draw = kTRUE;
                        
                    }
                }
                
                if(str2.Contains("5 < Sigma"))
                {
                    
                    if(trackType->HasChildren())
                    {
                        
                        TEveElementList* momentumVectorList3 = new TEveElementList("5 < sigma");
                        
                        TEveElement::List_i x = trackType->BeginChildren();
                        TEveElement::List_i y = trackType->EndChildren();
                        TEveElement::List_i z;
                        
                        for(z = x; z != y; z++)
                        {
                            
                            AliEveTrack* trackSingle1 = dynamic_cast<AliEveTrack*>((TEveElement*) *z);
                            
                            TEveLine* momentumVector = new TEveLine(TString::Format("Momentum Vector"));
                            
                            x1 = trackSingle1->GetESDTrack()->Xv();
                            y1 = trackSingle1->GetESDTrack()->Yv();
                            z1 = trackSingle1->GetESDTrack()->Zv();
                            
                            momentumVector->SetPoint(0, x1, y1, z1);
                            
                            x2 = x1+vectorLength*trackSingle1->GetESDTrack()->Px();
                            y2 = y1+vectorLength*trackSingle1->GetESDTrack()->Py();
                            z2 = z1+vectorLength*trackSingle1->GetESDTrack()->Pz();
                            
                            momentumVector->SetPoint(1, x2, y2, z2);
                            /*
                             if(trackSingle1->GetESDTrack()->Charge() == -1)
                             momentumVector->SetLineColor(kGreen+3);
                             else
                             momentumVector->SetLineColor(kRed+3);
                             */
                            momentumVector->SetLineColor(kRed+3);
                            
                            momentumVector->SetLineWidth(1);
                            momentumVector->SetLineStyle(0);
                            momentumVector->SetTitle(Form("%f GeV/c", trackSingle1->GetESDTrack()->P()));
                            
                            momentumVectorList3->AddElement(momentumVector);
                            
                        }
                        
                        //gEve->AddElement(momentumVectorList3);
                        momentumVectorList->AddElement(momentumVectorList3);
                        
                        draw = kTRUE;
                        
                    }
                }
            }
        }
    }
    
    gEve->AddElement(momentumVectorList);
    
    TEveElement* top = gEve->GetCurrentEvent();
    
    AliEveMultiView *mv = AliEveMultiView::Instance();
    
    mv->DestroyEventRPhi();
    mv->DestroyEventRhoZ();
    
    mv->ImportEventRPhi(top);
    mv->ImportEventRhoZ(top);
    
    gEve->FullRedraw3D(kFALSE, kTRUE);
    
}
