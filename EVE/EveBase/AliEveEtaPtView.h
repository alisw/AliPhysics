// $Id$
// Author: Jeremi Niedziela 2015

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveEtaPtView_H
#define AliEveEtaPtView_H

#include <TEveManager.h>

#include <TEveViewer.h>
#include <TGLViewer.h>

#include <TEveScene.h>
#include <TEveGeoShape.h>

#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>

#include <TEveBrowser.h>
#include <TEveWindow.h>

class AliEveEtaPtView
{
public:
    AliEveEtaPtView();
    ~AliEveEtaPtView();
    
    static AliEveEtaPtView* Instance();
    
    void InitGeom();
    
    TEveViewer* Get3DView()   { return f3DView; }
private:
    TEveViewer *f3DView;
    TEveWindowPack *fPack;
    TEveGeoShape *fGeom;
    
    static AliEveEtaPtView* fInstance;
    
    AliEveEtaPtView(const AliEveEtaPtView&);
    AliEveEtaPtView& operator=(const AliEveEtaPtView&);
    
    ClassDef(AliEveEtaPtView, 0);
};

#endif
