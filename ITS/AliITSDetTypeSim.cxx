/***************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
 $Id$
*/

/*********************************************************************
 * This class defines the "Standard" reconstruction for the ITS 
 * detector.
 **********************************************************************/
#include "AliITSDetTypeSim.h"

ClassImp(AliITSDetTypeSim)

//----------------------------------------------------------------------
AliITSDetTypeSim::AliITSDetTypeSim():
TObject(),
fGeom(),         //
fSimulation(),   // [NDet]
fSegmentation(), // [NDet]
fResponse(),     // [NMod]
fPreProcess(),   // [] e.g. Fill fHitModule with hits
fPostProcess(),  // [] e.g. Wright Raw data
fHitModule(),    //! [NMod][Nhits]
fNhits(0),       //! number of hits
fHits(),         //! local pointer
fNSDigits(0),    //! number of SDigits
fSDigits(),      //! [NMod][NSDigits]
fNDigits(0),     //! number of Digits
fDigits(),       //! [NMod][NDigits]
fHitClassName(), // String with Hit class name.
fSDigClassName(),// String with SDigit class name.
fDigClassName(){ // String with digit class name.
    // Default Constructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    A properly zero-ed AliITSDetTypeSim class.
}
//----------------------------------------------------------------------
AliITSDetTypeSim::~AliITSDetTypeSim(){
    // Destructor
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    Nothing.

    delete fGeom;
    delete fSimulation;
    delete fSegmentation;
    delete fResponse;
    delete fPreProcess;
    delete fPostProcess;
    delete fHitModule;
    // Do not delete fSDigits Not owned by this class see AliITS
    // Do not delete fDigits Not owned by this class see AliITS
}
//----------------------------------------------------------------------
AliITSDetTypeSim::AliITSDetTypeSim(const AliITSDetTypeSim &s) : TObject(s){
    // Copy Constructor for object AliITSDetTypeSim
    // Inputs:
    //   AliITSDetTypeSim &s  The copy Sourse
    // Outputs:
    //   none.
    // Return:
    //   A new AliITSDetTypeSim Object with the same data as that in the
    //   source s.

    if(&s==this) return;
    *this = source;
    return;
}
//----------------------------------------------------------------------
AliITSDetTypeSim::operator=(const AliITSDetTypeSim &s) : TObject(s){
    // The = operator for object AliITSDetTypeSim
    // Inputs:
    //   AliITSDetTypeSim &s  The copy Sourse
    // Outputs:
    //   none.
    // Return:
    //   A new AliITSDetTypeSim Object with the same data as that in the
    //   source s.

    if(&s==this) return;
    // Make copies of the objects in the arrays as well
    this->fGeom            = new AliITSgeom(s.fGeom);// Create a new instance
    if(this->fSimulation!=0) delete this->fSimulation;
    this->fSimulation      = new TObjArray(s.fSimulation->GetSize(),
                                           s.fSimulation->LowerBound());
    for(i=0;i<s.fSimulation->GetSize();i++) if(s.fSimulation[i]!=0)
        this->fSimulation[i] = new AliITSsimulation(*(s.fSimulation[i]));
        else this->fSimulation[i] = 0;
    if(this->fSegmentation!=0) delete this->fSegentation;
    this->fSegmentation = new TObjArray(s.fSegmentation->GetSize(),
                                        s.fSegmentation->GetLowerBound());
    for(i=0;i<s.fSegmentation->GetSize();i++) if(s.fSegmentation[i]!=0)
        this->fSegmentation[i] = new AliITSsegmentation(*(s.fSegmentation[i]));
        else this->fSegmentation[i] = 0;
    if(this->fResponse!=0) delete fResponse;
    this->fResponse = new TObjArray(s.fResponse->GetSize(),
                                    s.fResponse->GetLowerBound());
    for(i=0;i<s.fResponse->GetSize();i++) if(s.Response[i]!=0)
        this->fResponse[i] = new AliITSresponse(*(s.fResponse[i]));
        else this->fResponse[i] = 0;
    this->fPreProcess      = s.fPreProcess;  // Improper copy
    this->fPostProcess     = s.fPostProcess; // Improper copy
    this->fHitModule       = s.fHitModule;   // Improper copy
    this->fNhits           = s.fNhits;       //
    this->fHits            = s.fHits;        // copy pointer address only
    this->fNSDigits        = s.fNSDigits;    //
    this->fSDigits         = s.fSDigits;     // copy pointer address only
    this->fNDigits         = s.fNDigits;     //
    this->fDigits          = s.fDigits;      // copy pointer address only
    this->fHitClassName    = s.fHitClassName;
    this->fSDigitClassName = s.fSDigitClassName;
    this->fDigitClassName  = s.FDigitClassName;
    return *this;
}
//______________________________________________________________________
void AliITSDetTypeSim::InitModules(Int_t size,Int_t &nmodules){
    // Initialize the modules array.
    // Inputs:
    //      Int_t size  Size of array of the number of modules to be
    //                  created. If size <=0 then the number of modules
    //                  is gotten from AliITSgeom class kept in fGeom.
    // Outputs:
    //      Int_t &nmodules The number of modules existing.
    // Return:
    //      none.

    if(fHitModule){ 
        fHitModule->Delete();
        delete fHitModule;
    } // end fir fITSmoudles

    Int_t nl,indexMAX,index;

    if(size<=0){ // default to using data stored in AliITSgeom
        if(fGeom==0) {
            Error("InitModules","fGeom not defined");
            return;
        } // end if fGeom==0
        nl = fGeom->GetNlayers();
        indexMAX = fGeom->GetIndexMax();
        nmodules = indexMAX;
        fHitModule = new TObjArray(indexMAX);
        for(index=0;index<indexMAX;index++){
            fHitModule->AddAt( new AliITSmodule(index),index);
        } // end for index
    }else{
        fHitModule = new TObjArray(size);
        for(index=0;index<size;index++) {
            fHitModule->AddAt( new AliITSmodule(index),index);
        } // end for index

        nmodules = size;
    } // end i size<=0
}
//______________________________________________________________________
void AliITSDetTypeSim::FillModules(TTree *treeH, Int_t mask) {
    // fill the modules with the sorted by module hits; 
    // can be called many times to do a merging
    // Inputs:
    //      TTree *treeH  The tree containing the hits to be copied into
    //                    the modules.
    //      Int_t mask    The track number mask to indecate which file
    //                    this hits came from.
    // Outputs:
    //      none.
    // Return:
    //      none.

    if (treeH == 0x0){Error("FillModules","Tree is NULL");return;}

    Int_t lay,lad,det,index;
    AliITShit *itsHit=0;
    AliITSmodule *mod=0;
    TBranch *branch = treeH->GetBranch(fHitClassName.Data());
    if (!branch){Error("FillModules","%s branch in TreeH not found",
                       fHitClassName.Data());return;} // end if !branch
    branch->SetAddress(&fHits);
    Int_t nTracks =(Int_t) treeH->GetEntries();
    Int_t iPrimTrack,h;
    for(iPrimTrack=0; iPrimTrack<nTracks; iPrimTrack++){
        ResetHits();
        Int_t nBytes = treeH->GetEvent(iPrimTrack);
        if (nBytes <= 0) continue;
        Int_t nHits = fHits->GetEntriesFast();
        for(h=0; h<nHits; h++){
            itsHit = (AliITShit *)fHits->UncheckedAt(h);
            itsHit->GetDetectorID(lay,lad,det);
            if (fGeom) index = fGeom->GetModuleIndex(lay,lad,det);
            else index=det-1; // This should not be used.
            mod = GetModule(index);
            itsHit->SetTrack(itsHit->GetTrack()+mask); // Set track mask.
            mod->AddHit(itsHit,iPrimTrack,h);
        } // end loop over hits 
    } // end loop over tracks
}
