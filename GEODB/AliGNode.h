#ifndef ALIGNODE_H
#define ALIGNODE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TObjArray.h>
#include <TNamed.h>
#include <TArrayI.h>
#include <TFile.h>
#include "AliGTransform.h"
#include "AliGShape.h"
#include "AliGBox.h"
#include "AliGSphere.h"
#include "AliGTube.h"
#include "AliGCone.h"
#include "AliGPCone.h"
#include "AliGTRD1.h"
#include "AliGMaterial.h"
#include "AliGConfig.h"

class AliGNode : public TNamed {

protected:
    AliGConfig*   fConfig;
    Int_t         fId;        /* node id number                                    */
    AliGMaterial* fMaterial;  /* material                                          */
    TArrayI*      fNcopy;     /* copy numbers list (associated with the sons list) */
    TString       fNname;     /* node 'short' name                                 */
    TArrayI*      fNnode;     /* Sons list                                         */
    TObjArray*    fNode;      /* list of all nodes below                           */
    Int_t         fNsons;     /* Number of elements in fNnode                      */
    TArrayI*      fNtrans;    /* transf list (associated with the sons list)       */
    AliGNode*     fParent;    /* Pointer to the parent node                        */
    AliGShape*    fShape;     /* shape                                             */
    TObjArray*    fTrans;     /* list of all the transformations involved          */
    Int_t         fVis;	      /* Visibility	                                   */
    

public:
    AliGNode( AliGNode* node=NULL ); /* Copy or Default Constructor */

    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGBox* box, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for AliGBox shape */

    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGSphere* sphere, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for AliGSphere shape */

    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGTube* tube, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for AliGTube shape */

    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGCone* cone, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for AliGCone shape */
    
    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGPCone* pcone,AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for AliGPCone shape */
    
    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGShape* shape, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL ); /* AliGNode Constructor for an AliGShape */

    AliGNode( Text_t* name, Int_t id, Text_t* title, AliGTRD1* trd1, AliGMaterial* material=NULL, const Text_t* matrixname="", AliGConfig* config=NULL );

    virtual ~AliGNode(); /* Destructor */
    AliGNode* operator=( const AliGNode* node ); /* Operator = */
    virtual AliGMaterial *Material() const { return fMaterial;}
 
            void           Add(AliGNode *son,AliGTransform *tran);
            void           AddConf(AliGConfig* configname) { fConfig = configname;}
            void           AddConfig(Text_t* name, Text_t* title, Text_t* detail="",Int_t beg=0, Int_t end=0);
            void           AddMaterial(AliGMaterial* materialname) { fMaterial = materialname;}
            void           AddShape( AliGShape* shapename) { fShape = shapename;}
            Int_t          DistancetoPrimitive( Int_t px, Int_t py);
            void           Draw( Option_t* option);
            void           DrawShape( Option_t* option); /* MENU   */
            void           EndConstructor( Text_t* name, Int_t id, Text_t* title, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config );
            AliGConfig*    GetConfig() {return fConfig;}
	    
	    TArrayI*       GetfNnode() {return fNnode;}
	    TObjArray*     GetfNode()  {return fNode;}
            Int_t          GetfNsons() {return fNsons;}
	    TArrayI*       GetfNtrans(){return fNtrans;}
	    TArrayI*       GetfNcopy(){return fNcopy;}
	    TObjArray*     GetfTrans() {return fTrans;}
	    AliGMaterial*  GetMaterial() {return fMaterial;}
            AliGNode*      GetNodeFromfNnode(int position) {return (AliGNode*)
	    this->fNode->At(fNnode->At(position));}
	    const Text_t*        GetFamName() {return fNname;}
	    AliGNode*      GetNodeFromfNode(int position) {return (AliGNode*)
	    this->fNode->At(position);}
            AliGNode*      GetParent() {return fParent;}
            Text_t*        GetPath();
            AliGShape*     GetShape()  {return fShape;}
            AliGTransform* GetTransFromfTrans(int position) {return (AliGTransform*) this->fTrans->At(position);}
            //void           IncreaseSize(const Int_t total, TArrayI** array);
            void	   IncreaseSize(Int_t total,TArrayI &array);
	    void           Paint( Option_t *option);
            void           Save(TFile* file);
            void           SaveAll(TFile* file);
            Int_t          SizefNode() {return this->fNode->GetSize();}
            //void           Streamer(TBuffer &b);
	    void           SetVis(Int_t val) {fVis = val;}

    ClassDef(AliGNode,1) //Node class
};

//R__EXTERN TArrayF* gMatrix;

#endif
