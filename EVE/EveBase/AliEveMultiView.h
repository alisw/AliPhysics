// $Id$
// Author: Matevz Tadel 2009

/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

#ifndef AliEveMultiView_H
#define AliEveMultiView_H

#include <TEveGeoShape.h>
#include <TEveScene.h>
#include <TEveViewer.h>
#include <TGLViewer.h>

#include <vector>

class AliEveMultiView
{
public:
    static AliEveMultiView* Instance(){return fgInstance;}
    
    AliEveMultiView();
    ~AliEveMultiView();
    
    void InitSimpleGeom(TEveGeoShape* geom,bool threeD=true, bool rPhi=true, bool rhoZ=true);
    
    void ImportEvent(TEveElement *el);
    void DestroyAllEvents();
    void DestroyAllGeometries();
    
    void ImportEvent3D(TEveElement* el);
    void ImportEventRPhi(TEveElement* el);
    void ImportEventRhoZ(TEveElement* el);
    
    void DestroyEvent3D();
    void DestroyEventRPhi();
    void DestroyEventRhoZ();
    
    TEveViewer* Get3DView()   { return f3DView; }
    TEveViewer* GetRPhiView() { return fRPhiView; }
    TEveViewer* GetRhoZView() { return fRhoZView; }
    
    TEveScene* GetRPhiScene() {return fRPhiGeomScene;}
    TEveScene* GetRhoZScene() {return fRhoZGeomScene;}
    
private:
    static AliEveMultiView* fgInstance;     // Single instance of the multiview
    
    TEveProjectionManager *fRPhiMgr; // Manager of R-Phi projection
    TEveProjectionManager *fRhoZMgr; // Manager of Rho-Z projection
    
    TEveViewer            *f3DView;   // 3D view
    TEveViewer            *fRPhiView; // R-Phi view
    TEveViewer            *fRhoZView; // Pho-Z view
    
    TEveScene             *f3DGeomScene;    // 3D geometry scene
    TEveScene             *fRPhiGeomScene;  // R-Phi geometry scene
    TEveScene             *fRhoZGeomScene;  // Rho-Z geometry scene
    TEveScene             *f3DEventScene;   // 3D event scene
    TEveScene             *fRPhiEventScene; // R-Phi event scene
    TEveScene             *fRhoZEventScene; // Rho-Z event scene
    
    std::vector<TEveGeoShape*> fGeomVector;
    
    
    AliEveMultiView(const AliEveMultiView&);            // Not implemented
    AliEveMultiView& operator=(const AliEveMultiView&); // Not implemented
    
    ClassDef(AliEveMultiView, 0);
};

#endif
