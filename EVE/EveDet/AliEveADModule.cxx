/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The drawing module for the AD detector                            //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "AliEveADModule.h"
#include "AliADConst.h"

#include <TEveManager.h>
#include <TEveTrans.h>
#include <TEveRGBAPaletteOverlay.h>
#include <TGLViewer.h>
#include <TGLAnnotation.h>

static const Float_t xSize = 18.1;
static const Float_t ySize = 21.6;
static const Float_t zSize = 2.54;
static const Float_t xGap = 0.27;
static const Float_t yGap = 0.05;
static const Float_t zGap = 0.23+2.54;
static const Int_t xModule[4] = {1,1,-1,-1};
static const Int_t yModule[4] = {1,-1,-1,1};
static const Int_t zLayer[2] = {-1,1};
static const Float_t ACHole[2] = {6.2,3.7};

ClassImp(AliEveADModule)

/******************************************************************************/
AliEveADModule::AliEveADModule(const Text_t* name, Bool_t side, Int_t maxCharge, Bool_t showLegend)
: TEveQuadSet(name),
fStream(NULL),
fIsASide(side),
fShowLegend(showLegend)
{
    TEveRGBAPalette* rawPalette  = new TEveRGBAPalette(0, maxCharge);
    rawPalette->SetLimits(0, maxCharge);
    SetPalette(rawPalette);
    Reset(TEveQuadSet::kQT_FreeQuad, kFALSE, 16);
}

/******************************************************************************/
AliEveADModule::~AliEveADModule()
{
    delete fStream;
}
/******************************************************************************/
void AliEveADModule::LoadEsd(AliESDAD *adESD)
{
    //
    // Load AD ESD
    //
    if (!adESD) return;
    
    for (Int_t iChannel=0; iChannel < 16; ++iChannel)
    {
        Int_t iLayer; Int_t iModule; Float_t rHole;
        if (fIsASide) {
            if (iChannel < 8) continue;
            iLayer = (iChannel-8)/4;
            iModule = (iChannel-8)%4;
            rHole = ACHole[0];
            
        }
        else {
            if (iChannel >= 8) continue;
            iLayer = iChannel/4;
            iModule = iChannel%4;
            rHole = ACHole[1];
            
        }
        Float_t v[12];
        
        v[ 0] = (rHole + xGap)*xModule[iModule]; v[ 1] = yGap*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = (rHole + xGap)*xModule[iModule]; v[ 4] = (yGap + ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (xGap + xSize)*xModule[iModule]; v[ 7] = (yGap + ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (xGap + xSize)*xModule[iModule]; v[10] = yGap*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(adESD->GetAdc(iChannel));
        
        v[ 0] = xGap*xModule[iModule]; v[ 1] = (rHole + yGap)*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = xGap*xModule[iModule]; v[ 4] = (yGap + ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (rHole/2 + xGap)*xModule[iModule]; v[ 7] = (yGap+ ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (rHole/2 + xGap)*xModule[iModule]; v[10] = (rHole + yGap)*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(adESD->GetAdc(iChannel));
        
        v[ 0] = (rHole/2 + xGap)*xModule[iModule]; v[ 1] = (rHole + yGap)*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = (rHole/2 + xGap)*xModule[iModule]; v[ 4] = (yGap+ ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (rHole + xGap)*xModule[iModule]; v[ 7] = (yGap + ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (rHole + xGap)*xModule[iModule]; v[10] = (rHole/2 + yGap)*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(adESD->GetAdc(iChannel));
    }
    
    if (fIsASide)
        RefMainTrans().SetPos(0, 0, 1699.7);
    else
        RefMainTrans().SetPos(0, 0, -1954.4);
    
    if(fShowLegend){
        TEveRGBAPalette* pal = this->GetPalette();
        TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(pal, 0.69, 0.1, 0.3, 0.05);
        TGLViewer* v = gEve->GetDefaultGLViewer();
        v->AddOverlayElement(po);
        TGLAnnotation *ann = new TGLAnnotation(v,"Integrated charge [ADC counts]",0.69, 0.2);
        ann->SetTextSize(0.04);
    }
}

/******************************************************************************/
void AliEveADModule::LoadRaw(AliRawReader *rawReader)
{
    //
    // Load AD raw-data
    //
    if (fStream) delete fStream;
    fStream = new AliADRawStream(rawReader);
    if (!fStream->Next()) {
        delete fStream;
        fStream = NULL;
        return;
    }
    
    for (Int_t iChannel=0; iChannel < AliADRawStream::kNChannels; ++iChannel) {
        Int_t offChannel = kOfflineChannel[iChannel];
        Int_t iLayer; Int_t iModule; Float_t rHole;
        if (fIsASide) {
            if (offChannel < 8) continue;
            iLayer = (offChannel-8)/4;
            iModule = (offChannel-8)%4;
            rHole = ACHole[0];
            
        }
        else {
            if (offChannel >= 8) continue;
            iLayer = offChannel/4;
            iModule = offChannel%4;
            rHole = ACHole[1];
            
        }
        Float_t v[12];
        
        v[ 0] = (rHole + xGap)*xModule[iModule]; v[ 1] = yGap*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = (rHole + xGap)*xModule[iModule]; v[ 4] = (yGap + ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (xGap + xSize)*xModule[iModule]; v[ 7] = (yGap + ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (xGap + xSize)*xModule[iModule]; v[10] = yGap*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(fStream->GetADC(iChannel));
        
        v[ 0] = xGap*xModule[iModule]; v[ 1] = (rHole + yGap)*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = xGap*xModule[iModule]; v[ 4] = (yGap + ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (rHole/2 + xGap)*xModule[iModule]; v[ 7] = (yGap+ ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (rHole/2 + xGap)*xModule[iModule]; v[10] = (rHole + yGap)*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(fStream->GetADC(iChannel));
        
        v[ 0] = (rHole/2 + xGap)*xModule[iModule]; v[ 1] = (rHole + yGap)*yModule[iModule]; v[ 2] = zGap*zLayer[iLayer];
        v[ 3] = (rHole/2 + xGap)*xModule[iModule]; v[ 4] = (yGap+ ySize)*yModule[iModule]; v[ 5] = zGap*zLayer[iLayer];
        v[ 6] = (rHole + xGap)*xModule[iModule]; v[ 7] = (yGap + ySize)*yModule[iModule]; v[ 8] = zGap*zLayer[iLayer];
        v[ 9] = (rHole + xGap)*xModule[iModule]; v[10] = (rHole/2 + yGap)*yModule[iModule]; v[11] = zGap*zLayer[iLayer];
        
        AddQuad(v);
        QuadValue(fStream->GetADC(iChannel));
    }
    
    if (fIsASide) 
        RefMainTrans().SetPos(0, 0, 1699.7);
    else
        RefMainTrans().SetPos(0, 0, -1954.4);
    
    if(fShowLegend){
        TEveRGBAPalette* pal = this->GetPalette();
        TEveRGBAPaletteOverlay *po = new TEveRGBAPaletteOverlay(pal, 0.69, 0.1, 0.3, 0.05);
        TGLViewer* v = gEve->GetDefaultGLViewer();
        v->AddOverlayElement(po);
        TGLAnnotation *ann = new TGLAnnotation(v,"Amplitude of measured charge [ADC counts]",0.69, 0.2);
        ann->SetTextSize(0.04);
    }
}


