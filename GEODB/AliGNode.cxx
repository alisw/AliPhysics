/**************************************************************************
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
$Log$
*/

/**********************************************/
/*                                            */
/* FILE: AliGNode.cxx                         */
/* PURPOSE: Tree elemnt of the geometry       */
/* LANGUAGE: C++                              */
/* COMPILER: CC for HP-UX 9.x and 10.         */
/* AUTHOR: Joana && David                     */
/* DATE: May 28, 1999                         */
/* ADDRESS: jesanto@cern.ch, dcollado@cern.ch */
/*                                            */
/**********************************************/

#include <TCanvas.h>
#include <TView.h>
#include <TPadView3D.h>
#include <iostream.h>
#include <TFile.h>
#include <TROOT.h>
#include "AliGeometry.h"
#include "AliGBox.h"
#include "AliGCone.h"
#include "AliGPCone.h"
#include "AliGSphere.h"
#include "AliGNode.h"

const Int_t kSonsInvisible = BIT(17);
Text_t* ROOT_FILE = "AliGeoDB.root";

ClassImp(AliGNode)

//----------------------------------------------------------------------

AliGNode::AliGNode( AliGNode* node ) 
{

    if( node ) {
        /* Copy Constructor */
        fConfig   = new AliGConfig( node->fConfig );
        fId       = node->fId;
        fMaterial = new AliGMaterial( node->fMaterial );
        fName     = node->fName;
        fNcopy    = new TArrayI( *(node->fNcopy) );
        fNname    = node->fNname;
        fNnode    = new TArrayI( *(node->fNnode) );
        fNode     = new TObjArray( *(node->fNode) );
        fNsons    = node->fNsons;
        fNtrans   = new TArrayI( *(node->fNtrans) );
	fParent   = node->fParent;
	fVis      = node->fVis;
        fTitle    = node->fTitle;
        fTrans    = new TObjArray(*(node->fTrans));
	//        fParent   = new AliGNode(node->fParent);

        if (node->fShape) {
            if( !strcmp(node->fShape->ClassName(), "AliGBox") ) {
                fShape = new AliGBox( (AliGBox*) node->fShape );
                //fShape = new AliGBox(*(AliGBox*)node->fShape);
            }
            else {
                if( !strcmp(node->fShape->ClassName(), "AliGSphere") )
                    fShape = new AliGSphere( (AliGSphere*) node->fShape );
                else {
                    if( !strcmp(node->fShape->ClassName(), "AliGTube") )
                        fShape = new AliGTube( (AliGTube*) node->fShape );
                    else {
                        if( !strcmp(node->fShape->ClassName(), "AliGCone") )
                            fShape = new AliGCone( (AliGCone*) node->fShape );
			else {
                            if( !strcmp(node->fShape->ClassName(), "AliGPCone") )
                                fShape = new AliGPCone( (AliGPCone*) node->fShape );
			    else {
                                if( !strcmp(node->fShape->ClassName(), "AliGTRD1") ) 
                                    fShape = new AliGTRD1( (AliGTRD1*) node->fShape );
			    }
			}
                    }
                }
            }
	}
        else
            fShape = NULL;

        /*cout << endl << " Comparacion." << endl;
        cout << " El que paso por parametro:" << endl;
        node->Dump();
        cout << " El resultado de la copia:" << endl;
        this->Dump();*/

    }
    else {
        /* Default Constructor */
        fConfig   = NULL;
        fId       = 0;
        fMaterial = NULL;
        fName     = "";
        fNcopy    = NULL;
        fNname    = "";
        fNnode    = NULL;
        fNode     = new TObjArray();
        fNsons    = 0;
        fNtrans   = NULL;
        fParent   = NULL;
        fShape    = NULL;
        fTitle    = "";
        fTrans    = new TObjArray();
	fVis      = 0;
    }
}

//----------------------------------------------------------------------

AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGBox* box, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGBox shape */

    if( box )
        fShape = new AliGBox(box); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------

AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGSphere* sphere, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGSphere shape */

    if( sphere )
        fShape = new AliGSphere(sphere); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------

AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGTube* tube, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGTube shape */

    if( tube )
        fShape = new AliGTube(tube); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------

AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGCone* cone, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGCone shape */

    if( cone )
        fShape = new AliGCone(cone); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------


AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGPCone* pcone, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGPCone shape */

    if( pcone )
        fShape = new AliGPCone(pcone); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------


AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGTRD1* trd1, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* AliGNode Constructor for AliGTRD1 shape */

    if( trd1 )
	fShape = new AliGTRD1(trd1); // Shape inside node
    else
        fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}

//----------------------------------------------------------------------

void AliGNode::EndConstructor( Text_t* name, Int_t id, Text_t* title, AliGMaterial *material, const Text_t *matrixname, AliGConfig* config )
{
    /* Finish the construction of AliGNodes */
    
    fConfig   = new AliGConfig(config);
    fId       = id;                         // Node id number
    fMaterial = new AliGMaterial(material); // Material inside shape

    char*  TmpName  = new char[strlen(name)+int(TMath::Log10((double) id)+3)];
    sprintf(TmpName,"%s_%d",name,id);
    fName     = TmpName;             // Complete node name
    delete [] TmpName;

    fNcopy    = new TArrayI(50);            // Array of integers containing the copy number (by family name)
    fNname    = name;                       // Node family name
    fNnode    = new TArrayI(50);            // Array of integers containing every node son index
    fNode     = new TObjArray(100,0);       // Array containing every (different)node son   
    fNsons    = 0;
    fNtrans   = new TArrayI(50);            // Array of integers containing every transformation index
    fParent   = NULL;
    fTitle    = title;
    fTrans    = new TObjArray(100,0);       // Array containing every (different)transformation used
    fVis      = 1;

}

//-----------------------------------------------------------------------

AliGNode::AliGNode( Text_t* name, Int_t id, Text_t* title, AliGShape* shape, AliGMaterial* material, const Text_t* matrixname, AliGConfig* config ) 
{
    /* AliGNode Constructor for an AliGShape */

    fShape = NULL;

    EndConstructor( name, id, title, material, matrixname, config );
}


//----------------------------------------------------------------------

AliGNode::~AliGNode()
{
    /* Destructor */

    if(fConfig)   delete fConfig;
    if(fMaterial) delete fMaterial;
    if(fNcopy)    delete fNcopy;
    if(fNnode)    delete fNnode;
    if(fNode)     delete fNode;
    if(fNtrans)   delete fNtrans;
    if(fShape)    delete fShape;
    if(fTrans)    delete fTrans;
}

//----------------------------------------------------------------------

AliGNode* AliGNode::operator=( const AliGNode* node )
{
    if( this == node ) return this; // special case.

    fConfig   = node->fConfig;
    fId       = node->fId;
    fMaterial = node->fMaterial;
    fName     = node->fName;
    fNcopy    = node->fNcopy;
    fNname    = node->fNname;
    fNnode    = node->fNnode;

    if( fNode ) fNode->Delete();

    if( node->fNode ) {
      //for( int i=0; i<node->fNode->GetSize(); i++ )
        for( int i=0; i<node->fNode->GetSize(); i++ )
            fNode->AddLast(node->fNode->At(i));
    }

    fNsons    = node->fNsons;
    fNtrans   = node->fNtrans;
    fParent   = node->fParent;
    fShape    = node->fShape;
    fTitle    = node->fTitle;

    if( fTrans ) fTrans->Delete();
    if( node->fTrans ) {
        for( int i=0; i<node->fTrans->GetSize(); i++ )
            fTrans->AddLast(node->fTrans->At(i));
    }
    //fTrans    = node->fTrans;

    return this;
}

//----------------------------------------------------------------------

void AliGNode::Add(AliGNode* son, AliGTransform* tran)
{
    
    Int_t index = fNode->IndexOf((TObject*)son);
    
    if( index == -1 ) {
        fNode->Add(son);  // fNode is TObjArray
        index = fNode->IndexOf((TObject*)son);
        son->fParent = this;
    }

    Int_t tranumber = fTrans->IndexOf((TObject*)tran); 
    if( tranumber == -1 ) {
        fTrans->Add( tran );
        tranumber = fTrans->IndexOf((TObject*)tran);
    }

    Int_t Nsons = this->GetfNsons();
    Int_t copy = 1;
   
    for( Int_t i = 0; i<Nsons; i++ ) {
        if( index == fNnode->At(i) ) {
            if( fNtrans->At(i) == tranumber ) {
	      cout << "Error: Node already exists!";
	        return;
            }
        }
	
	
        if( ((AliGNode*)fNode->At(fNnode->At(i)))->fNname == son->fNname ) 
	    copy++;
	    
    }

    if( fNnode->GetSize() == Nsons ) {
    
        IncreaseSize( Nsons+100, *fNnode  );
        IncreaseSize( Nsons+100, *fNtrans );
        IncreaseSize( Nsons+100, *fNcopy  );
	
    }
    fNnode->AddAt( index, Nsons ); // Add value of index at position Nsons
    
    fNtrans->AddAt( tranumber, Nsons ); 
    fNcopy->AddAt( copy, Nsons );
    fNsons++; // Because we have added a new node
   
}

//----------------------------------------------------------------------   

void AliGNode::IncreaseSize(Int_t total,TArrayI &array)
{
    TArrayI *newarray = new TArrayI(array);
    array.Set(total);    

    for( int i=0; i<newarray->GetSize(); i++ )
        array[i] = (*newarray)[i];
        
    delete newarray;
}

//----------------------------------------------------------------------

Text_t* AliGNode::GetPath()
{
    const int kMAXDEPTH = 128;
    const AliGNode* d[kMAXDEPTH];
    AliGNode* node = this;
    int depth = 0, len = 0;

    d[depth++] = node;
    len = strlen(node->GetName()) + 1;

    while (node->fParent && depth < kMAXDEPTH) {
        node = node->fParent;     
        d[depth++] = node;
        len += strlen (node->GetName()) + 1;
    }

    char* path = new char[len + 2];

    for (int i = depth-1; i>=0; i--) {
        if (i == depth-1) {   // file or TROOT name 
            strcpy(path, d[i]->GetName());
            if (i == 0) strcat(path, "/");
        } 
        else {
            strcat(path, "/");
            strcat(path, d[i]->GetName());
        }
    }
 
    return path;  
} 


//----------------------------------------------------------------------

void AliGNode::Save( TFile* file )
{
    if( fParent ) { // it's not the root node
        // So I must check if the father node is in disk

        if( file->cd(fParent->GetPath()) ) {

            if( !gDirectory->Get(GetName()) )
                file->mkdir( GetName() );

            file->cd( GetPath() );

            if( !gDirectory->Get( fShape->GetName() ) )
                if ( fShape ) {
/*                    const Text_t* shapetype = fShape->ClassName();
                    if( !strcmp(shapetype,"AliGBox") ) {
                        fShape->Write();
                    }
                    else {
                        if (!strcmp(shapetype,"AliGSphere")) {
                            fShape->Write();
                        }
                        else {
	                    if (!strcmp(shapetype,"AliGTube")) {
                                fShape->Write();
	                    }
                            else {
	                        if (!strcmp(shapetype,"AliGCone")) {
                                    fShape->Write();
                                }
                                else {
                                    if (!strcmp(shapetype,"AliGPCone")) {
                                        cout << " Saving an AliGPCone:" << endl;
                                        fShape->Dump();
                                        cout << endl;
                                        fShape->Write();
                                    }
                                    else {
                                        if (!strcmp(shapetype,"AliGTRD1")) {
                                            fShape->Write();
                                        }
                                    }
                                }
                            }
                        }
                    }*/
                    fShape->Write();
                }


            if( !gDirectory->Get( fMaterial->GetName() ) )
                if ( fMaterial )
                    fMaterial->Write();

            if ( fConfig ) {
                /*cout << " Con capacidad: " << fConfig->GetFormula().Capacity() << " voy a grabar esta configuracion: " << fConfig->GetFormula().Data() << endl;*/
                fConfig->Write();
            }

            if ( fTrans )
                fTrans->Write();
        }
        else
	    cout << " ERROR : Father exists but I can't find it. See AliGNode::Save, please..." << endl;
    }
    else { // if it's the root node
        file->mkdir( GetName() );
        file->cd( GetPath() );

        if ( fShape )
            fShape->Write();

        if ( fMaterial )
            fMaterial->Write();

        if ( fConfig ) {
            /*cout << " Con capacidad: " << fConfig->GetFormula().Capacity() << " voy a grabar esta configuracion: " << fConfig->GetFormula().Data() << endl;*/
            fConfig->Write();
        }

        if ( fTrans )
            fTrans->Write();
    }

    fflush(stdout);
    file->Write();
}

//----------------------------------------------------------------------

void AliGNode::SaveAll( TFile* file )
{
/*    if( fParent )
        file->cd( fParent->GetPath() );
    else
        file->cd();

    Save( file );

    for( int i=0; i<fNsons; i++ ) {
        cout << " Saving node: " << GetPath() << " with copy number " << << endl;
        GetNodeFromfNnode(i)->SaveAll( file );
    }
*/

    if( !gDirectory->Get(GetName()) ) {
        Save( file );

        for( int i=0; i<fNsons; i++ )
	    GetNodeFromfNnode(i)->SaveAll( file );
    }

/*    if (!file->cd(GetPath()))  {
        Save( file );

        for( int i=0; i<fNsons; i++ )
	    GetNodeFromfNnode(i)->SaveAll( file );
    }
*/
}

//----------------------------------------------------------------------

void AliGNode::AddConfig( Text_t* name, Text_t* title, Text_t* detail, Int_t beg, Int_t end )
{
    Int_t Nsons = GetfNsons();

    TStringLong formula = "";

    for( int i=0; i<Nsons; i++ ) {
        Int_t index = fNnode->At(i);  // Gets every node son's index
        const Text_t* NodeName = ((AliGNode*)fNode->At(index))->GetName();
        Int_t tranumber = fNtrans->At(i); // And the correspondent transform
        const Text_t* TransfName = ((AliGTransform*)fTrans->At(tranumber))->GetName();

        formula.Append(NodeName);
        formula.Append(":");
        formula.Append(TransfName);

        if( i != fNsons-1 )
            formula.Append("+");
    }

    //cout << " Acabo de crear la formula: " << formula << endl << endl;

    const Text_t* shapetype    = fShape->ClassName();
    const Text_t* shapename    = fShape->GetName();
    const Text_t* materialname = fMaterial->GetName();

    AliGConfig* tmp = new AliGConfig( name, title, formula, detail, shapetype, shapename, materialname, beg, end );

    //if(fConfig) delete fConfig;

    fConfig = tmp;
}

//----------------------------------------------------------------------

/*void AliGNode::AddConfig( Text_t* name, Text_t* title, Text_t* detail, Int_t beg, Int_t end )
{
    // Creates the configuration of a node using the information passed as arguments
    // and building the formula of the level below by building the list of its sons
    // recursively.

    Int_t Nsons = GetfNsons();

    char* formula = NULL;

    cout << " Nsons = " << Nsons << endl;

    for( int i=0; i<Nsons; i++ ) {
        Int_t index = fNnode->At(i);  // Gets every node son's index
        const Text_t* NodeName = ((AliGNode*)fNode->At(index))->GetName();
        Int_t tranumber = fNtrans->At(i); // And the correspondent transform
        const Text_t* TransfName = ((AliGTransform*)fTrans->At(tranumber))->GetName();

        int lenF, lenN, lenT;
        lenF = lenN = lenT = 0;

        cout << " i = " << i << " and fNsons = " << fNsons << endl;

        if ( formula ) {
            lenF = strlen(formula);
            cout << " formula = #" << formula << "#" << endl;
        }
        
        if ( NodeName ) {
            lenN = strlen(NodeName);
            cout << " NodeName = #" << NodeName << "#" << endl;
        }
        
        if ( TransfName ) {
            lenT = strlen(TransfName);
            cout << " TransfName = #" << TransfName << "#" << endl;
        }
        
        cout << " lenF: " << lenF << " - lenN: " << lenN << " - lenT: " << lenT << endl;

        char* formula2;

        if( i != fNsons-1 ) {
            cout << " Arriba" << endl;
            formula2 = new char[lenF+lenN+lenT+2];
        }
        else {
            cout << " Embaixo" << endl;
            formula2 = new char[lenF+lenN+lenT+1];
        }
        
        formula2[0] = '\x0';

        cout << "@@@ :) @@@" << endl;

        if ( formula ) {
            strcpy( formula2, formula );
            cout << "1 formula2 = #" << formula2 << "#" << endl;
        }

        if( NodeName ) {
            strcat( formula2, NodeName );
            cout << "2 formula2 = #" << formula2 << "#" << endl;
        }
      
        strcat( formula2, ":" );

        if( TransfName ) {
            strcat( formula2, TransfName );
            cout << "3 formula2 = #" << formula2 << "#" << endl;
        }

        if( i != fNsons-1 ) {
            strcat( formula2, "+" );
            cout << "4 formula2 = #" << formula2 << "#" << endl;
        }

        //strcat( formula2, "\x0" );
        //cout << "5 formula2 = #" << formula2 << "#" << endl;        

        cout << " Y formula2 tiene un tamanyo de: " << strlen(formula2) << " caracteres." << endl;

        //if( formula )  delete [] formula;
        formula = new char[strlen(formula2)];
        strcpy( formula, formula2 );
        if( formula2 ) delete [] formula2;
    }

    const Text_t* shapetype    = fShape->ClassName();
    const Text_t* shapename    = fShape->GetName();
    const Text_t* materialname = fMaterial->GetName();

    AliGConfig* tmp = new AliGConfig( name, title, formula, detail, shapetype, shapename, materialname, beg, end );

    if(fConfig) delete fConfig;

    fConfig = tmp;
}
*/

//-------------------------------------------------------------------------

void AliGNode::DrawShape(Option_t *option)
{
    Draw(option);
    gPad->Update();
}

//----------------------------------------------------------------------

void AliGNode::Draw(Option_t* option)
{

   TString opt = option;
   opt.ToLower();
//*-*- Clear pad if option "same" not given
   if (!gPad) {
      if (!gROOT->GetMakeDefCanvas()) return;
      (gROOT->GetMakeDefCanvas())();
   }
   if (!opt.Contains("same")) gPad->Clear();

//*-*- Draw Referenced node
   //gGeometry->SetGeomLevel();
   //gGeometry->UpdateTempMatrix();

   AppendPad(option);

//*-*- Create a 3-D View
   TView *view = gPad->GetView();
   if (!view) {
      view = new TView(1);
      view->SetAutoRange(kTRUE);
      Paint(option);
      view->SetAutoRange(kFALSE);
   }
}

//------------------------------------------------------------------

void AliGNode::Paint(Option_t *option)
{ 
  //     Print();
     TPadView3D *view3D=gPad->GetView3D();

   //Int_t level = gGeometry->GeomLevel();
// restrict the levels for "range" option
   //if (level > 3 && strcmp(option,"range") == 0) return;
//*-*- Update translation vector and rotation matrix for new level
   //if (level) {
   //   gGeometry->UpdateTempMatrix(fX,fY,fZ,fMatrix->GetMatrix(),fMatrix->IsReflection());
   //   if (view3D)
   //       view3D->UpdateNodeMatrix(this,option);
  // }
  
//*-*- Paint Referenced shape
//   Int_t vis = fShape->GetVisibility();
//   if ( vis == -1) return;
//   if (vis == -3) {
//     if (nsons == 0) vis = 1;
//     else            vis = 0;
//   }

  // TAttLine::Modify();
  // TAttFill::Modify();
   //if (fVisibility && fShape->GetVisibility()) {
      gNode = this;
      //fShape->SetLineColor(GetLineColor());
      //fShape->SetLineStyle(GetLineStyle());
      //fShape->SetLineWidth(GetLineWidth());
      //fShape->SetFillColor(GetFillColor());
      //fShape->SetFillStyle(GetFillStyle());
      
   Int_t nsons = 0;
   if (fNnode) nsons = GetfNsons();
      if (view3D)
         view3D->SetAttNode((TNode*)this,option);
      
     
    if(fVis) {
        if( fShape )
            fShape->Paint(option);
        else
            cout << " Intente dibujar donde no habia shape!!! " << endl;
   }
   //}
   if ( TestBit(kSonsInvisible) ) return;

//*-*- Paint all sons
   if(!nsons) return;
    
   //gGeometry->PushLevel();
   for( int i=0; i<nsons; i++ ) {
      
	gAliGeometry->PushMatrix(gMatrix);

        AliGNode* node = (AliGNode*) fNode->At(fNnode->At(i));
	
	AliGTransform* trans = (AliGTransform*) fTrans->At(fNtrans->At(i));
        gAliGeometry->UpdateMatrix(trans);
        node->Paint(option);
	gAliGeometry->PopMatrix();
	
    }
   //gGeometry->PopLevel();
   
}
   
   
//--------------------------------------------------------------------

Int_t AliGNode::DistancetoPrimitive(Int_t px, Int_t py)
{
    TView *view = gPad->GetView();
    gPad->SetSelected(view);
    return 0;
}

//--------------------------------------------------------------------

