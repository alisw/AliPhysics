/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, all rights reserved. *
 * See http://aliceinfo.cern.ch/Offline/AliRoot/License.html for          *
 * full copyright notice.                                                 *
 **************************************************************************/
#ifndef AliEveADModule_H
#define AliEveADModule_H

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// The drawing module for the AD detector                               //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <TEveQuadSet.h>

#include <AliRawReader.h>
#include <AliESDAD.h>
#include <AliADRawStream.h>

class AliEveADModule : public TEveQuadSet
{
public:
    AliEveADModule(const Text_t* name="AliEveADModule",
                   Bool_t side = kTRUE,
                   Int_t maxCharge = 1023,
                   Bool_t showLegend = kFALSE);
    ~AliEveADModule();
    
    void LoadRaw(AliRawReader *rawReader);
    void LoadEsd(AliESDAD *adESD);
    
private:
    
    AliADRawStream *fStream;      // Raw-stream
    Bool_t          fIsASide;     // A or C side module
    Bool_t	        fShowLegend;  // Legend
    
    AliEveADModule(const AliEveADModule&);
    AliEveADModule& operator=(const AliEveADModule&);
    
    ClassDef(AliEveADModule,1) // Representation of a AD module
};

#endif
