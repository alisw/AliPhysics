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
/* FILE: AliGeometry.cxx                      */
/* PURPOSE: To store and retrieve a complete  */
/*          geometry in/from a root file.     */
/* LANGUAGE: C++                              */
/* COMPILER: CC for HP-UX 9.x and 10.         */
/* AUTHOR: Joana && David                     */
/* DATE: May 28, 1999                         */
/* ADDRESS: jesanto@cern.ch, dcollado@cern.ch */
/*                                            */
/**********************************************/

#include <iostream.h>
#include <TFile.h>
#include <TROOT.h>
#include "AliGConfig.h"
#include "AliGeometry.h"
#include "AliGBox.h"
#include "AliGSphere.h"
#include "AliGCone.h"
#include "AliGTube.h"
#include "AliGTransform.h"

#define FormLeng 80

//float matrix[16] = {1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.};
//TArrayF* gMatrix = new TArrayF( 16, matrix );
TVector *gMatrix = new TVector(0,15,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1., "END");
AliGeometry *gAliGeometry = 0;

ClassImp(AliGeometry)


// ---------------------------------------------------------------------------

AliGeometry::AliGeometry( AliGeometry* Geom )
{
    if( Geom ) {
        /* Copy Constructor */
        fBomb           = Geom->fBomb; //Bomb factor for exploded geometry
        fGeomLevel      = Geom->fGeomLevel;
        fMaxDepth       = Geom->fMaxDepth;
        fName           = Geom->fName;

        fRules          = new TList();
        if( Geom->fRules ) {
            for( int i=0; i<Geom->fRules->GetSize(); i++ )
                fRules->AddLast(Geom->fRules->At(i));
        }

        fTitle          = Geom->fTitle;
        fTopNode        = Geom->fTopNode;

        gAliGeometry    = this;
        fTransformation = new TList();
        if( Geom->fTransformation ) {
            for( int i=0; i<Geom->fTransformation->GetSize(); i++ )
                fTransformation->AddLast(Geom->fTransformation->At(i));
        }
    }
    else {
        /* Default constructor */
        fBomb           = 1.; //Bomb factor for exploded geometry
        fGeomLevel      = 0;
        fMaxDepth       = 0; // 0 = No Depth limitation
        fName           = "";
        fRules          = new TList(); // Antes estaba a NULL
        fTitle          = "";
        fTopNode        = "";
        fTransformation = new TList(); // Antes estaba a NULL
	gAliGeometry = this;
    }
}

// ---------------------------------------------------------------------------

AliGeometry::~AliGeometry()
{
    /* Destructor */
    
    if( fRules ) {
        fRules->Delete();
        fRules = NULL;
    }

    if( fTransformation ) {
        fTransformation->Delete();
        fTransformation = NULL;
    }
}

// ---------------------------------------------------------------------------

AliGeometry* AliGeometry::operator=( const AliGeometry* Geom )
{
    /* Operator = */
    if( this == Geom ) return this; // special case.

    fBomb           = Geom->fBomb; //Bomb factor for exploded geometry
    fGeomLevel      = Geom->fGeomLevel;
    fMaxDepth       = Geom->fMaxDepth;
    fName           = Geom->fName;
    gAliGeometry    = this;

    if( fRules )
        fRules->Delete();

    if( Geom->fRules ) {
        for( int i=0; i<Geom->fRules->GetSize(); i++ )
            fRules->AddLast(Geom->fRules->At(i));
    }

    fTitle   = Geom->fTitle;
    fTopNode = Geom->fTopNode;

    if( fTransformation )
        fTransformation->Delete();

    fTransformation = new TList();

    if( Geom->fTransformation ) {
        for( int i=0; i<Geom->fTransformation->GetSize(); i++ )
            fTransformation->AddLast(Geom->fTransformation->At(i));
    }

    return this;
}

// ---------------------------------------------------------------------------

AliGNode* AliGeometry::FileToMemTree( TFile* file )
{
    // Geometry constructor that retrieves a whole geometry from the root file
    // and stores it in memory as a tree structure

    // Read the name of my geometry's first node
    file->cd();
    const char* node_name = (const char*) GetTop().Data();

    // Call recursion
    AliGNode* new_node = new AliGNode( RecurToMem( node_name, file, NULL, NULL ) );

    return new_node;
}

// ---------------------------------------------------------------------------

AliGNode* AliGeometry::RecurToMem( const char* node_name, TFile* file, AliGNode* father, char *subtransf ) 
{
    // Recursive function used to retrieve the root_file structure into memory

    const int LenName = strlen(node_name);

    int i;
    for( i=0; i<LenName; i++ )
        if( node_name[i] == '_' )
            break;
    
    char* name = new char[i+1];
    strncpy( name, node_name, i );
    name[i] = '\x0';

    char* ids = new char[LenName-i];
    for( int g1=0, g2=i+1; g2<LenName; g1++, g2++ )
        ids[g1] = node_name[g2];

    ids[LenName-i-1] = '\x0';

    const int id = atoi( ids );
    if( ids )
        delete [] ids;

    AliGNode* node = NULL;

    int new_node = 1;

    if (father) {
        AliGTransform* trans = NULL;

        int l;
        for( l=0; father->GetNodeFromfNode(l); l++ ) {
            // check if the node has been already added
            if ( !strcmp( node_name, (father->GetNodeFromfNode(l))->GetName()) ) {
                new_node = 0;
                break;
            }
            else
                node = new AliGNode( (Text_t*)name, id, (Text_t*)name, (AliGBox*) NULL );
        }

        int new_trans = 1;
        int m;
        for( m=0; father->GetTransFromfTrans(m); m++ ) {
            // check if this tranformations has already been used with this father
            if (!strcmp(subtransf,(father->GetTransFromfTrans(m))->GetName()))  {
                new_trans=0;
                break;
            }
            else {
                trans = new AliGTransform();
                file->cd(father->GetPath());
                trans->Read(subtransf);
            }
        }

        if ( l == 0 ) node = new AliGNode( (Text_t*)name, id, (Text_t*)name, (AliGBox*) NULL );
        if ( m == 0 ) {
            trans = new AliGTransform();
            trans->Read(subtransf);
        }

        if (new_node) {
            if(!new_trans) {
                father->Add(node,father->GetTransFromfTrans(m));
                file->cd(node->GetPath());
            }
            else {
                father->Add(node,trans);
                file->cd(node->GetPath());
             }
        }
        else {
            if(!new_trans) {
                father->Add(father->GetNodeFromfNode(l),father->GetTransFromfTrans(m));
                file->cd(father->GetPath());
            }
            else {
                father->Add(father->GetNodeFromfNode(l),trans);
                file->cd(father->GetPath());
            }
        }

        //if( trans ) delete trans;
    }
    else { // There is no father
        node = new AliGNode( (Text_t*)name, id, (Text_t*)name, (AliGBox*) NULL );
        file->cd(node->GetPath());
    }

    if( name )
        delete [] name;

    if( new_node ) {
        // Read the configuration name
        char* configs = NULL;

        if( fRules ) {
            for( int k=0; k<fRules->GetSize(); k++ ) {
                char* rule = new char[strlen((char*)fRules->At(k)->GetName())];
                strcpy( rule, (char*)fRules->At(k)->GetName() );

                int w;
                for( w=0; w<strlen(rule); w++ )
                    if( rule[w] == '/' )
                        break;
    
                configs = new char[w+1];
                strncpy( configs, rule, w );
                configs[w] = '\x0';

                if( !strcmp(node->GetName(), configs) ) {
                    delete [] rule;
                    break;
                }
                else {
                    delete [] configs;
                    configs = NULL;
                }
                
                delete [] rule;
            }
        }
        else {
            printf( " ERROR: I couldn't find the fRule inside RecurToMem.\n" );
            return NULL;
        }

        // Build Configuration
        AliGConfig* Config = new AliGConfig();
        
        Config->Read(configs);
        node->AddConf(Config);

        // Read Shape, material and formula from configuration
        const Text_t* shapetype = Config->GetShapeType().Data();
        const Text_t* shapename = Config->GetShapeName().Data();

        // Build Shape

        if( !strcmp(shapetype,"AliGBox") ) {
            AliGBox* shape = new AliGBox();
            shape->Read(shapename);
            node->AddShape(shape);
        }
        else {
            if (!strcmp(shapetype,"AliGSphere")) {
                AliGSphere* shape = new AliGSphere();
                shape->Read(shapename);
                node->AddShape(shape);
            }
            else {
	        if (!strcmp(shapetype,"AliGTube")) {
                    AliGTube* shape = new AliGTube();
                    shape->Read(shapename);
                    node->AddShape(shape);
	        }
                else {
	            if (!strcmp(shapetype,"AliGCone")) {
                        AliGCone* shape = new AliGCone();
                        shape->Read(shapename);
                        node->AddShape(shape);
                    }
                    else {
                        if (!strcmp(shapetype,"AliGPCone")) {
                            AliGPCone* shape = new AliGPCone();
                            shape->Read(shapename);
                            node->AddShape(shape);
                        }    
                        else {
                            if (!strcmp(shapetype,"AliGTRD1")) {
                                AliGTRD1* shape = new AliGTRD1();
                                shape->Read(shapename);
                                node->AddShape(shape);
                            }
                        }
                    }
                }
            }
        }

        // Build Material
        const Text_t* materialname = Config->GetMaterialName().Data();
        AliGMaterial *material = new AliGMaterial();
        material->Read(materialname);
        node->AddMaterial(material);

        // Read formula
        const char* formula = Config->GetFormula();

        int Len = strlen(formula);
        int j = 0;

        while( formula[j] ) {
            char* subnode   = new char[Len];
            char* subtransf = new char[Len];
            int k;

            for( k=0; formula[j] != ':'; k++, j++ )
                subnode[k] = formula[j];

            subnode[k] = '\x0';
            if( formula[j] ) j++;

            for( k=0; formula[j] && formula[j] != '+'; k++, j++ )
                subtransf[k] = formula[j];

            subtransf[k] = '\x0';
            if( formula[j] ) j++;

            AliGNode* node1 = this->RecurToMem( subnode, file, node, subtransf );
       }
    }

    return node;
}

// ---------------------------------------------------------------------------

AliGeometry::AliGeometry( Text_t* name, Text_t* title, AliGNode* topNode, int maxDepth ): TNamed(name, title)
{
    /* Geometry constructor that stores in the root file the geometry file that corresponds to one geometry */

    gAliGeometry    = this;
    fBomb           = 1.; //Bomb factor for exploded geometry
    fGeomLevel      = 0;
    fMaxDepth       = maxDepth;
    fRules          = new TList();
    fTopNode        = topNode->GetName();
    fTransformation = new TList();

    RulesList( fRules, topNode, 0 );

}

// ---------------------------------------------------------------------------

void AliGeometry::RulesList( TList *fRules, AliGNode* node, int position ) 
{
    /* Recursive function used to create the list of rules for the Geometry file */

    int len = strlen(node->GetName()) + strlen(node->GetConfig()->GetName());

    char* path = new char[len+2];
    strcpy( path, node->GetName() );
    strcat( path, "/" );
    strcat( path, node->GetConfig()->GetName());
    TNamed* rule = new TNamed( path, "" );
    delete [] path;

    fRules->AddLast( rule );
    //cout << " Rule: " << rule->GetName() << " added." << endl;
    //if( rule ) delete rule;

    while( (position < node->SizefNode()) && (node->GetNodeFromfNode(position) != NULL) ) {
        this->RulesList( fRules, node->GetNodeFromfNode(position), 0 );
        position++;
    }

}

// ---------------------------------------------------------------------------

void AliGeometry::UpdateMatrix(AliGTransform* trans)
{
//    Update global rotation matrix/translation vector for this node
//   this function must be called before invoking Local2Master

        TVector* matrix = trans->GetMatrix(); // Get father's transf. matrix
	TVector* temp = new TVector( 0,15,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1.,0.,0.,0.,0.,1., "END");
      

	//	printf("Inside updatematrix, this is the son's matrix\n");
	//	matrix->Print();
	//	printf("Inside updatematrix, this is the global matrix before\n");
	//	gMatrix->Print();
        /* Multiplication of matrices*/
        for( int i=0; i<4; i++ )
            for( int j=0; j<4; j++ ) {
	       double tmp = 0;
	        for( int k=0; k<4; k++ ) 
		  //        tmp += (*matrix)(i*4+k) * (*gMatrix)(j + k*4);
	            tmp += (*gMatrix)(i*4+k) * (*matrix)(j + k*4);
		(*temp)(i*4+j)= tmp;
            }
	//if(gMatrix) delete gMatrix;
	    gMatrix = temp;
	
	    //	printf("Inside updatematrix, this is the global matrix after\n");
	    //	gMatrix->Print();
	
	    //	for (int k = 0; k < 16; k++) {
	    //cout << "gMatrix[" << k << "]" << (*gMatrix)(k) << endl;
	    //cout << "matrix[" << k << "]" << (*matrix)(k) << endl;
	    //}

	    //if (temp)   delete [] temp;
        //if (matrix) delete [] matrix;
 
}

// ---------------------------------------------------------------------------

void AliGeometry::Local2Master(Float_t *local, Float_t *master)
{
// *-*-*-*-*Convert one point from local system to master reference system*-*-*
// *-*      ==============================================================
//  Note that before invoking this function, the  matrix of the transformation
//  for this node must have been computed.
//  AliGeometry::UpdateMatrix should be called before.


  //    gMatrix->Print();
    if( gNode->GetParent() && gMatrix) {
        
        Double_t tmp;
	Float_t loc[4];

	for(int i=0;i<3;i++) loc[i]=local[i];
	loc[3]=1;
	
	TVector &mat = *gMatrix;

        for (int i=0; i<3; i++) {
	    tmp=0;
            for (int j=0; j<4; j++) tmp+= mat(i*4+j) * loc[j];
	    master[i] = tmp;
        }
    }
    else
        memcpy( master, local, sizeof(Double_t)* kVectorSize );	 
}

//----------------------------------------------------------------------------

void AliGeometry::PushMatrix(TVector* matrix)
{
    //if( fTransformation == NULL )
        //fTransformation = new TList();

    //cout << " Va a dar error " << endl;
    if (matrix) fTransformation->AddLast((TObject*)matrix);
    //cout << " Ves como aqui no llega? :-) " << endl;
}

//----------------------------------------------------------------------------

TVector* AliGeometry::PopMatrix()
{
  
  gMatrix = (TVector*)fTransformation->Last();
  fTransformation->Remove(fTransformation->Last());
  return gMatrix;
}

//----------------------------------------------------------------------------



