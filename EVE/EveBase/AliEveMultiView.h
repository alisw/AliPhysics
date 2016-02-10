// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMultiView_H
#define AliEveMultiView_H

#include <TEveManager.h>

#include <TEveViewer.h>
#include <TGLViewer.h>

#include <TEveScene.h>
#include <TEveGeoShape.h>

#include <TEveProjectionManager.h>
#include <TEveProjectionAxes.h>

#include <TEveBrowser.h>
#include <TEveWindow.h>

#include <vector>

//______________________________________________________________________________
// Short description of AliEveMultiView
//

class AliEveMultiView
{
public:
    AliEveMultiView();
    virtual ~AliEveMultiView();
    
    static AliEveMultiView* Instance();
    
    void InitGeomGentle(TEveGeoShape* g3d, TEveGeoShape* grphi, TEveGeoShape* grhoz);
    void InitGeomGentleTrd(TEveGeoShape* gtrd);
    void InitSimpleGeom(TEveGeoShape* geom, bool rPhi=true, bool rhoZ=true);
    
    //-------------------------------------------------------------------------
    
    void SetDepth(Float_t d);
    
    //-------------------------------------------------------------------------
    
    void ImportGeomRPhi(TEveElement* el);
    void ImportGeomRhoZ(TEveElement* el);
    void ImportEventRPhi(TEveElement* el);
    void ImportEventRhoZ(TEveElement* el);
    
    void DestroyEventRPhi();
    void DestroyEventRhoZ();
    
    void SetCenterRPhi(Double_t x, Double_t y, Double_t z);
    void SetCenterRhoZ(Double_t x, Double_t y, Double_t z);
    
    void DestroyAllGeometries();
    
    //-------------------------------------------------------------------------
    
    TEveViewer* Get3DView()   { return f3DView; }
    TEveViewer* GetRPhiView() { return fRPhiView; }
    TEveViewer* GetRhoZView() { return fRhoZView; }
    
    TEveScene* GetRPhiScene() {return fRPhiGeomScene;}
    TEveScene* GetRhoZScene() {return fRhoZGeomScene;}
    
    TEveProjectionManager* GetRhoZMgr() { return fRhoZMgr;}
    
    
protected:
    TEveProjectionManager *fRPhiMgr; // Obvious meaning.
    TEveProjectionManager *fRhoZMgr; // Obvious meaning.
    
    TEveViewer            *f3DView;   // Obvious meaning.
    TEveViewer            *fRPhiView; // Obvious meaning.
    TEveViewer            *fRhoZView; // Obvious meaning.
    
    TEveWindowPack *fPack;
    
    TEveScene             *fRPhiGeomScene;  // Obvious meaning.
    TEveScene             *fRhoZGeomScene;  // Obvious meaning.
    TEveScene             *fRPhiEventScene; // Obvious meaning.
    TEveScene             *fRhoZEventScene; // Obvious meaning.
    
    std::vector<TEveGeoShape*> fGeomVector;
    
    TEveGeoShape          *fGeomGentle;     // Obvious meaning.
    TEveGeoShape          *fGeomGentleRPhi; // Obvious meaning.
    TEveGeoShape          *fGeomGentleRhoZ; // Obvious meaning.
    TEveGeoShape          *fGeomGentleTrd;  // Obvious meaning.
        
    static AliEveMultiView* fgInstance;     // Obvious meaning.
    
private:
    AliEveMultiView(const AliEveMultiView&);            // Not implemented
    AliEveMultiView& operator=(const AliEveMultiView&); // Not implemented
    
    ClassDef(AliEveMultiView, 0); // Multiple-views.
};

#endif
