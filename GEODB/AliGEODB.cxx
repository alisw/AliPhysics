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
Revision 1.3  1999/11/03 13:17:05  fca
Have ProdProcess return const char*

Revision 1.2  1999/09/29 09:24:19  fca
Introduction of the Copyright and cvs Log

*/

//////////////////////////////////////////////////////
//  C++ dummy interface to Geant3 basic routines    //
//////////////////////////////////////////////////////

#include <iostream.h>
#include <TBrowser.h>
#include <TCanvas.h>
#include "AliGEODB.h"
#include "AliGNode.h"
#include "AliGShape.h"
#include "AliGMaterial.h"
#include "AliGTransform.h"
#include "AliGTube.h"
#include "AliGBox.h"
#include "AliGTRD1.h"
#include "AliGPCone.h"
#include "AliGeometry.h"

//---------------------------------------------------------

#define MAXMAT 1000

struct Mater{
    char*    fName;
    Float_t  fA;
    Float_t  fZ;
    Float_t  fDens;
    Float_t  fRadl;
    Float_t  fAbsl;
    Float_t* fBuf;
    Int_t    fNwbuf;
};

Mater* AliMaterials[MAXMAT]; // Array of Mater structs containing the AliMaterial parameters.

struct Med{
    char*    fName;
    Int_t    fNmat;
    Int_t    fIsvol;
    Int_t    fIfield;
    Float_t  fFieldm;
    Float_t  fTmaxfd;
    Float_t  fStemax;
    Float_t  fDeemax;
    Float_t  fEpsil;
    Float_t  fStmin;
    Float_t* fUbuf;
    Int_t    fNbuf;
};

Med* AliMediums[MAXMAT]; // Array of Med structs containing the AliMedium parameters.

struct NodePosp{
    char* name;
    char* shape;
    Int_t nmed;
};

NodePosp* Posp[MAXMAT];

AliGMaterial* Mat[MAXMAT]; // List of Pointers to my AliGMaterials created

Int_t  MateCount   = 0;
Int_t  MedCount    = 0;
Int_t  MatrCount   = 0;
Int_t  PospCount   = 0;
Int_t  Color       = 2;

TList* listNodes     = new TList(); // List of Nodes
TList* listTransf    = new TList();
TList* listShapes    = new TList();
TList* listMaterials = new TList();

AliGNode* TopNode    = NULL;

/***********************************************************************/

ClassImp(AliGEODB)

/***********************************************************************/

AliGEODB::AliGEODB(const char *title, Int_t) : AliMC("AliGEODB",title)
{
    //cout << " AliGEODB" << endl;
}

/***********************************************************************/

void AliGEODB::DefaultRange() 
{
    cout << " DefaultRange" << endl;
}
 
//=======================functions from GBASE

/***********************************************************************/

void    AliGEODB::Gfile(const char*, const char*) 
{
    cout << " Gfile" << endl;
}

/***********************************************************************/

void    AliGEODB::GeomIter() 
{
    cout << " GeomIter" << endl;
}

/***********************************************************************/

Int_t   AliGEODB::CurrentMaterial(Float_t &, Float_t &, Float_t &, Float_t &, Float_t &) const 
{
    cout << " CurrentMaterial" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::NextVolUp(Text_t*, Int_t&) 
{
    cout << " NextVolUp" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::CurrentVol(Text_t*, Int_t&) const 
{
    cout << " CurrentVol" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::NofVolumes() const 
{
    cout << " NofVolumes" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::CurrentVolOff(Int_t, Text_t*, Int_t&) const 
{
    cout << " CurrentVolOff" << endl;
    return 0;
}

/***********************************************************************/

void    AliGEODB::TrackPosition(Float_t*) const 
{
    cout << " TrackPosition" << endl;
}

/***********************************************************************/

void    AliGEODB::TrackMomentum(Float_t*) const 
{
    cout << " TrackMomentum" << endl;
}

/***********************************************************************/

Int_t   AliGEODB::VolId( Text_t* name ) const 
{
    cout << " VolId" << endl;
    return 0;
}

/***********************************************************************/

const char* AliGEODB::VolName(Int_t ) const 
{
    cout << " VolName" << endl;
    return 0;
}

/***********************************************************************/
    
Float_t AliGEODB::TrackCharge() const 
{
    cout << " TrackCharge" << endl;
    return 0;
}

/***********************************************************************/

Float_t AliGEODB::TrackMass() const 
{
    cout << " TrackMass" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackInside() const 
{
    cout << " TrackInside" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackEntering() const 
{
    cout << " TrackEntering" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackExiting() const 
{
    cout << " TrackExiting" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackOut() const 
{
    cout << " TrackOut" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackDisappear() const 
{
    cout << " TrackDisappear" << endl;
    return 0;
}

/***********************************************************************/

Bool_t  AliGEODB::TrackStop() const 
{
    cout << " TrackStop" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::NSecondaries() const 
{
    cout << " NSecondaries" << endl;
    return 0;
}

/***********************************************************************/

AliMCProcess AliGEODB::ProdProcess() const 
{
    cout << " ProdProcess" << endl;
    return 0;
}

/***********************************************************************/

void    AliGEODB::GetSecondary(Int_t, Int_t&, Float_t*, Float_t*)
{
    cout << " GetSecondary" << endl;
}

/***********************************************************************/

Float_t AliGEODB::MaxStep() const 
{
    cout << " MaxStep" << endl;
    return 0;
}

/***********************************************************************/

void    AliGEODB::SetMaxStep(Float_t ) 
{
    cout << " SetMaxStep" << endl;
}

/***********************************************************************/

void    AliGEODB::GetParticle(const Int_t, char*, Float_t&) const 
{
    cout << " GetParticle" << endl;
}

/***********************************************************************/

Int_t   AliGEODB::CurrentEvent() const 
{
    cout << " CurrentEvent" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::GetMedium() const 
{
    cout << " GetMedium" << endl;
    return 0;
}

/***********************************************************************/

Float_t AliGEODB::Edep() const 
{
    cout << " Edep" << endl;
    return 0;
}

/***********************************************************************/
    
Float_t AliGEODB::Etot() const 
{
    cout << " Etot" << endl;
    return 0;
}

/***********************************************************************/

void    AliGEODB::Rndm(Float_t*, const Int_t) const 
{
    cout << " Rndm" << endl;
}

/***********************************************************************/

Float_t AliGEODB::TrackStep() const 
{
    cout << " TrackStep" << endl;
    return 0;
}

/***********************************************************************/
    
Float_t AliGEODB::TrackLength() const 
{
    cout << " TrackLength" << endl;
    return 0;
}

/***********************************************************************/

Float_t AliGEODB::TrackTime() const 
{
    cout << " TrackTime" << endl;
    return 0;
}

/***********************************************************************/

Int_t   AliGEODB::TrackPid() const 
{
    cout << " TrackPid" << endl;
    return 0;
}

/***********************************************************************/
    
Bool_t  AliGEODB::TrackAlive() const 
{
    cout << " TrackAlive" << endl;
    return 0;
}

/***********************************************************************/

void    AliGEODB::StopTrack() 
{
    cout << " StopTrack" << endl;
}

/***********************************************************************/

void    AliGEODB::StopEvent() 
{
    cout << " StopEvent" << endl;
}

/***********************************************************************/

void    AliGEODB::SetMaxNStep(Int_t) 
{
    cout << " SetMaxNStep" << endl;
}

/***********************************************************************/

void    AliGEODB::SetColors() 
{
    cout << " SetColors" << endl;
}

/***********************************************************************/

Int_t   AliGEODB::GetMaxNStep() const 
{
    cout << " GetMaxNStep" << endl;
    return 0;
}

/***********************************************************************/

void AliGEODB::Material( Int_t& kmat, const char* name, Float_t a, Float_t z, Float_t dens, Float_t radl, Float_t absl,Float_t* buf, Int_t nwbuf )
{
    if( MateCount >= MAXMAT ) {
        printf( " ERROR: Material number: %d bigger that MAXMAT: %d => Out of array\n", MateCount, MAXMAT );
        exit(0);
    }

    Mater*   AliMat = new Mater;
    char*    fname  = new char[strlen(name)+1];
    Float_t* fbuf   = new Float_t[nwbuf];

    strcpy( fname, name );
    fname[strlen(name)] = '\x0';

    AliMat->fName  = fname;
    AliMat->fA     = a;
    AliMat->fZ     = z;
    AliMat->fDens  = dens;
    AliMat->fRadl  = radl;
    AliMat->fAbsl  = absl;
    AliMat->fNwbuf = nwbuf;

    for( int i=0; i<nwbuf; i++ )
        fbuf[i] = buf[i];

    AliMat->fBuf = fbuf;
    
    AliMaterials[MateCount] = AliMat;

    kmat = MateCount++;
}

/***********************************************************************/

void  AliGEODB::Mixture( Int_t& imat, const char* name, Float_t* a, Float_t* z, Float_t dens, Int_t nlmat, Float_t* wmat )
{
    if( MateCount >= MAXMAT ) {
        printf( " ERROR: Mixture number: %d bigger that MAXMAT: %d => Out of array\n", MateCount, MAXMAT );
        exit(0);
    }

    Mater* AliMat   = new Mater;
    char*    fname  = new char[strlen(name)+1];

    strcpy( fname, name );
    fname[strlen(name)] = '\x0';

    AliMat->fName   = fname;
    Float_t MidA    = 0.;
    Float_t MidZ    = 0.;
    Float_t SUMwmat = 0.;

    for( int i=0; i<nlmat; i++ ) {
        MidA    += wmat[i]*a[i];
        MidZ    += wmat[i]*z[i];
        SUMwmat += wmat[i];
    }

    MidA            = MidA/SUMwmat;
    MidZ            = MidZ/SUMwmat;

    AliMat->fA      = MidA;
    AliMat->fZ      = MidZ;
    AliMat->fDens   = dens;
    AliMat->fRadl   = 0.;
    AliMat->fAbsl   = 0.;
    AliMat->fNwbuf  = 0;

    AliMat->fBuf    = NULL;

    AliMaterials[MateCount] = AliMat;

    imat = MateCount++;
}

/***********************************************************************/

void  AliGEODB::Medium( Int_t& numed, const char* name, Int_t nmat, Int_t isvol, Int_t ifield, Float_t fieldm, Float_t tmaxfd, Float_t stemax, Float_t deemax, Float_t epsil, Float_t stmin, Float_t* ubuf, Int_t nbuf )
{
    if( MedCount >= MAXMAT ) {
        printf( " ERROR: Medium number: %d bigger that MAXMAT: %d => Out of array\n", MedCount, MAXMAT );
        exit(0);
    }

    Med*     AliMed = new Med;
    char*    fname  = new char[strlen(name)+1];
    Float_t* fUbuf  = new Float_t[nbuf];

    strcpy( fname, name );
    fname[strlen(name)] = '\x0';
    
    AliMed->fName   = fname;
    AliMed->fNmat   = nmat;
    AliMed->fIsvol  = isvol;
    AliMed->fIfield = ifield;
    AliMed->fFieldm = fieldm;
    AliMed->fTmaxfd = tmaxfd;
    AliMed->fStemax = stemax;
    AliMed->fDeemax = deemax;
    AliMed->fEpsil  = epsil;
    AliMed->fStmin  = stmin;
    AliMed->fNbuf   = nbuf;
    
    for( int i=0; i<nbuf; i++ )
        fUbuf[i] = ubuf[i];

    AliMed->fUbuf = fUbuf;

    AliMediums[MedCount] = AliMed;

    numed = MedCount++;
}

/***********************************************************************/

void  AliGEODB::Matrix( Int_t &nmate , Float_t theta1, Float_t phi1, Float_t theta2, Float_t phi2, Float_t theta3, Float_t phi3 ) 
{
    char* transname = new char[10];
    sprintf(transname, "%s%d", "tra", MatrCount);

    AliGTransform* tra = new AliGTransform( transname, transname, theta1, phi1, theta2, phi2, theta3, phi3 );

    //printf( " Created new AliGTransform( %s, %s, %f, %f, %f, %f, %f, %f );\n", transname, transname, theta1, phi1, theta2, phi2, theta3, phi3 );

    listTransf->AddAt( tra, MatrCount );

    nmate = ++MatrCount;

}

/***********************************************************************/

void  AliGEODB::WriteEuclid(const char*, const char*, Int_t, Int_t) 
{
    cout << " WriteEuclid" << endl;
}

/***********************************************************************/

void  AliGEODB::Gpcxyz() 
{
    cout << " Gpcxyz " << endl;
}

/***********************************************************************/

void  AliGEODB::Ggclos() 
{

    for( int i=0; i<listNodes->GetSize(); i++ ) {
        AliGNode* node = (AliGNode*) listNodes->At(i);
        char* name = new char[strlen(node->GetName())+1];
        strcpy(name, node->GetName());
        name[strlen(node->GetName())] = '\x0';
        node->AddConfig(name, name,"detail",981201,991201);
    }

    printf( " Total number of nodes = %d\n", listNodes->GetSize() );

    SetVisibility( "TPC_1",   1 );
    SetVisibility( "TGAS_1",  1 );
    SetVisibility( "TPSG_1",  1 );
    SetVisibility( "TPHV_1",  1 );
    SetVisibility( "TRCS_1",  1 );
    SetVisibility( "TSGA_1",  1 );
    SetVisibility( "TRCL_1",  1 );
    SetVisibility( "TLGA_1",  1 );
    SetVisibility( "TSST_1",  1 );
    SetVisibility( "TSST_2",  1 );
    SetVisibility( "TSST_3",  1 );
    SetVisibility( "TSST_4",  1 );
    SetVisibility( "TSST_5",  1 );
    SetVisibility( "TSST_6",  1 );
    SetVisibility( "TSST_7",  1 );
    SetVisibility( "TSST_8",  1 );
    SetVisibility( "TSST_9",  1 );
    SetVisibility( "TSST_10", 1 );
    SetVisibility( "TSST_11", 1 );
    SetVisibility( "TSST_12", 1 );
    SetVisibility( "TSST_13", 1 );
    SetVisibility( "TSST_14", 1 );
    SetVisibility( "TSST_15", 1 );
    SetVisibility( "TSST_16", 1 );
    SetVisibility( "TSST_17", 1 );
    SetVisibility( "TSST_18", 1 );
    SetVisibility( "TSST_19", 1 );
    SetVisibility( "TSST_20", 1 );
    SetVisibility( "TSST_21", 1 );
    SetVisibility( "TSST_22", 1 );
    SetVisibility( "TSST_23", 1 );
    SetVisibility( "TLST_1",  1 );
    SetVisibility( "TLST_2",  1 );
    SetVisibility( "TLST_3",  1 );
    SetVisibility( "TLST_4",  1 );
    SetVisibility( "TLST_5",  1 );
    SetVisibility( "TLST_6",  1 );
    SetVisibility( "TLST_7",  1 );
    SetVisibility( "TLST_8",  1 );
    SetVisibility( "TLST_9",  1 );
    SetVisibility( "TLST_10", 1 );
    SetVisibility( "TLST_11", 1 );
    SetVisibility( "TLST_12", 1 );
    SetVisibility( "TLST_13", 1 );
    SetVisibility( "TLST_14", 1 );
    SetVisibility( "TLST_15", 1 );
    SetVisibility( "TLST_16", 1 );
    SetVisibility( "TLST_17", 1 );
    SetVisibility( "TLST_18", 1 );
    SetVisibility( "TLST_19", 1 );
    SetVisibility( "TLST_20", 1 );
    SetVisibility( "TLST_21", 1 );
    SetVisibility( "TLST_22", 1 );
    SetVisibility( "TLST_23", 1 );
    SetVisibility( "TLST_24", 1 );
    SetVisibility( "TLST_25", 1 );
    SetVisibility( "TLST_26", 1 );
    SetVisibility( "TLST_27", 1 );
    SetVisibility( "TLST_28", 1 );
    SetVisibility( "TLST_29", 1 );
    SetVisibility( "TLST_30", 1 );
    SetVisibility( "TLST_31", 1 );
    SetVisibility( "TLST_32", 1 );
    SetVisibility( "TLST_33", 1 );
    SetVisibility( "TLST_34", 1 );
    SetVisibility( "TLST_35", 1 );
    SetVisibility( "TLST_36", 1 );
    SetVisibility( "TLST_37", 1 );
    SetVisibility( "TLST_38", 1 );
    SetVisibility( "TLST_39", 1 );
    SetVisibility( "TLST_40", 1 );
    SetVisibility( "TLST_41", 1 );
    SetVisibility( "TLST_42", 1 );
    SetVisibility( "TLST_43", 1 );
    SetVisibility( "TLST_44", 1 );
    SetVisibility( "TLST_45", 1 );
    SetVisibility( "TLST_46", 1 );
    SetVisibility( "TLST_47", 1 );
    SetVisibility( "TLST_48", 1 );
    SetVisibility( "TLST_49", 1 );
    SetVisibility( "TLST_50", 1 );
    SetVisibility( "TLST_51", 1 );
    SetVisibility( "TLST_52", 1 );
    SetVisibility( "TSWS_1",  1 );
    SetVisibility( "TSWS_2",  1 );
    SetVisibility( "TSWS_3",  1 );
    SetVisibility( "TPW1_1",  1 );
    SetVisibility( "TPW2_1",  1 );
    SetVisibility( "TPW3_1",  1 );
    SetVisibility( "TPW4_1",  1 );
    SetVisibility( "TSPI_1",  1 );
    SetVisibility( "TSP1_1",  1 );
    SetVisibility( "TSPO_1",  1 );
    SetVisibility( "TSP2_1",  1 );
    SetVisibility( "TSWH_1",  1 );
    SetVisibility( "TSW1_1",  1 );
    SetVisibility( "TCOV_1",  1 );
    SetVisibility( "TPOI_1",  1 );
    SetVisibility( "TSWS_4",  1 );
    SetVisibility( "TSWS_5",  1 );
    SetVisibility( "TSWS_6",  1 );
    SetVisibility( "TSWS_7",  1 );
    SetVisibility( "TSWS_8",  1 );
    SetVisibility( "TSWS_9",  1 );
    SetVisibility( "TSWS_10", 1 );
    SetVisibility( "TSWS_11", 1 );
    SetVisibility( "TSWS_12", 1 );
    SetVisibility( "TSWS_13", 1 );
    SetVisibility( "TPIV_1",  1 );
    SetVisibility( "TPVD_1",  1 );

    printf( " Visibilities activated \n" );

    AliGeometry* Geom1 = new AliGeometry( "Geom1", "Geom1", TopNode, 0 );

    TCanvas* c1 = new TCanvas( "TPC", "Geometry Shapes", 100, 100, 750, 750 );
    c1->Range(0,0,1,1);
    c1->SetFillColor(32); // Light Green
    c1->SetBorderSize(3);
    c1->SetBorderMode(0); // -1 (down) 0 (no) 1 (up)

    TopNode->Draw("same");

    TFile* file1 = new TFile( "GeoDB.root", "RECREATE" );
    cout << endl << " Storing GeoDB.root file.  Please, be patient...    :-)" << endl;
    TopNode->SaveAll(file1);

    file1->cd();

    if( Geom1 ) {
        Geom1->Write();
        cout << " Geometry saved in disk " << endl;
    }

    file1->Write();
    file1->Close();

    TFile* file2 = new TFile( "GeoDB.root", "READ" );
    file2->cd();
    AliGeometry* Geom2  = new AliGeometry( (AliGeometry*) file2->Get("Geom1") );

    // Retrieve the memory tree structure from AliGeoDB.root and stores it in memory, being tree_root the top node
    AliGNode* tree_root = new AliGNode( Geom2->FileToMemTree( file2 ) );
    //AliGNode* tree_root = new AliGNode( Geom2->FileToMemTree( "GeoDB.root", "Geom1" ) );

    TFile* file3 = new TFile( "GeoDB2.root", "RECREATE" );
    cout << " Storing GeoDB2.root file. Please, be patient again...   :-)" << endl;
    tree_root->SaveAll( file3 ); // Memory tree structure saved in disk (second time)

    file3->cd();

    if( Geom2 ) {
        Geom2->Write();
        cout << " Geometry saved in disk " << endl;
    }

    file3->Write();
    file3->Close();

    file2->Close();

    cout << " The two databases are stored in disk." << endl;

    TFile* file4 = new TFile( "GeoDB.root" , "READ" );
    TFile* file5 = new TFile( "GeoDB2.root", "READ" );

    TBrowser* a = new TBrowser();

    //DrawTree( TopNode, 0 );
}

/***********************************************************************/

void AliGEODB::SetVisibility(Text_t* name, Int_t val)
{
   AliGNode* node = (AliGNode*) listNodes->FindObject( name );

   if( node )
       node->SetVis(val);
   else
       printf( " ERROR!!! I couldn't find the node %s in the listNodes.\n", name );
}

/***********************************************************************/

void AliGEODB::DrawTree( AliGNode* topnode, Int_t tabs )
{

    for( int i=0; topnode->GetNodeFromfNode(i); i++ ) {
        for( int j=0; j<tabs; j++ )
            cout << " ";
        cout << topnode->GetNodeFromfNode(i)->GetName() << endl;
        tabs+=5;
        DrawTree( topnode->GetNodeFromfNode(i), tabs );
        tabs-=5;
    }

}

/***********************************************************************/

void  AliGEODB::Glast() 
{
    cout << " Glast " << endl;
}

/***********************************************************************/

void  AliGEODB::Gprint(const char*) 
{
    cout << " Gprint " << endl;
}

/***********************************************************************/

void  AliGEODB::Grun() 
{
    cout << " Grun " << endl;
}

/***********************************************************************/

void  AliGEODB::Gtrig() 
{
    cout << " Gtrig " << endl;
}

/***********************************************************************/

void  AliGEODB::Gtrigc() 
{
    cout << " Gtrigc " << endl;
}

/***********************************************************************/

void  AliGEODB::Gtrigi() 
{
    cout << " Gtrigi " << endl;
}

/***********************************************************************/

void  AliGEODB::Gwork(Int_t) 
{
    cout << " Gwork " << endl;
}

/***********************************************************************/

void  AliGEODB::Gzinit() 
{
    cout << " Gzinit " << endl;
}

//=======================functions from GCONS

/***********************************************************************/

void  AliGEODB::Gfmate(Int_t, char*, Float_t&, Float_t&, Float_t&, Float_t&, Float_t&, Float_t*, Int_t&)
{
    cout << " Gfmate " << endl;
}

/***********************************************************************/

void  AliGEODB::Gfpart(Int_t, char*, Int_t&, Float_t&, Float_t&, Float_t&)
{
    cout << " Gfpart " << endl;
}

/***********************************************************************/

void  AliGEODB::Gftmed(Int_t, char*, Int_t&, Int_t&, Int_t&, Float_t&, Float_t&,Float_t&, Float_t&, Float_t&, Float_t&, Float_t*, Int_t*)
{
    cout << " Gftmed " << endl;
}

/***********************************************************************/

void  AliGEODB::Gmate() {
    cout << " Gmate " << endl;
}

/***********************************************************************/

void  AliGEODB::Gpart() 
{
    AliGTransform* Identity = new AliGTransform("Identity","Identity","rot 0. 0. 0.");
    //printf( " Created new AliGTransform( Identity, Identity, rot 0. 0. 0. );\n" );

    listTransf->AddAt( Identity, MatrCount );
    MatrCount++;

    AliGBox*  box  = new AliGBox("box","box", 400,400,400);
    AliGNode* ALIC = new AliGNode( "ALIC", 1,"ALIC", box);
    TopNode = ALIC;

    listNodes->Add(ALIC);
}

/***********************************************************************/

void  AliGEODB::Gsdk(Int_t, Float_t*, Int_t*) {
    cout << " Gsdk " << endl;
}

/***********************************************************************/

void  AliGEODB::Gsmate(Int_t, const char*, Float_t, Float_t, Float_t, Float_t, Float_t) {
    cout << " Gsmate " << endl;
}

/***********************************************************************/

void  AliGEODB::Gsmixt(Int_t, const char*, Float_t*, Float_t*, Float_t, Int_t, Float_t*) {
    cout << " Gsmixt " << endl;
}

/***********************************************************************/

void  AliGEODB::Gspart(Int_t, const char*, Int_t,   Float_t, Float_t, Float_t) {
    cout << " Gspart " << endl;
}

/***********************************************************************/

void  AliGEODB::Gstmed(Int_t, const char*, Int_t, Int_t, Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) {
    cout << " Gstmed " << endl;
}

/***********************************************************************/

void  AliGEODB::Gstpar(Int_t, const char*, Float_t) 
{
    //cout << " Gstpar " << endl;
}

/***********************************************************************/

void  AliGEODB::Gsckov(Int_t, Int_t, Float_t *, Float_t *, Float_t *, Float_t *) 
{
    cout << " Gsclov " << endl;
}

/***********************************************************************/

//=======================functions from GKINE

void  AliGEODB::Gfkine(Int_t, Float_t*, Float_t*, Int_t&, Int_t&) 
{
    cout << " Gfkine" << endl;
}

/***********************************************************************/

void  AliGEODB::Gfvert(Int_t, Float_t*, Int_t&, Int_t&, Float_t&) 
{
    cout << " Gfvert" << endl;
}

/***********************************************************************/

Int_t AliGEODB::Gskine(Float_t*, Int_t, Int_t, Float_t*, Int_t)
{
    cout << " Gskine" << endl;
    return 0;
}

/***********************************************************************/

Int_t AliGEODB::Gsvert(Float_t*, Int_t, Int_t, Float_t*, Int_t)
{
    cout << " Gsvert" << endl;
    return 0;
}

/***********************************************************************/

//=======================functions from GPHYS

void  AliGEODB::Gphysi() 
{
    cout << " Gphysi" << endl;
}

/***********************************************************************/

//=======================functions from GTRAK

void  AliGEODB::Gdebug() 
{
    cout << " Gdebug" << endl;
}

/***********************************************************************/

void  AliGEODB::Gekbin() 
{
    cout << " Gekbin" << endl;
}

/***********************************************************************/

void  AliGEODB::Gfinds() 
{
    cout << " Gfinds" << endl;
}

/***********************************************************************/

void  AliGEODB::Gsking(Int_t) 
{
    cout << " Gsking" << endl;
}

/***********************************************************************/

void  AliGEODB::Gskpho(Int_t) 
{
    cout << " Gskpho" << endl;
}

/***********************************************************************/

void  AliGEODB::Gsstak(Int_t) 
{
    cout << " Gsstak" << endl;
}

/***********************************************************************/

void  AliGEODB::Gsxyz() 
{
    cout << " Gsxyz" << endl;
}

/***********************************************************************/

void  AliGEODB::Gtrack() 
{
    cout << " Gtrack" << endl;
}

/***********************************************************************/

void  AliGEODB::Gtreve() 
{
    cout << " Gtreve" << endl;
}

/***********************************************************************/

void  AliGEODB::Grndm(Float_t*, const Int_t) const 
{
    cout << " Grndm" << endl;
}

/***********************************************************************/

void  AliGEODB::Grndmq(Int_t&, Int_t&, const Int_t, const Text_t*) 
{
    cout << " Grndmq" << endl;
}

/***********************************************************************/

//=======================functions from GDRAW

void  AliGEODB::Gdxyz(Int_t ) 
{
    cout << " Gdxyz" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdcxyz() 
{
    cout << " Gdcxyz" << endl;
}

/***********************************************************************/

//=======================functions from GGEOM

void  AliGEODB::Gdtom(Float_t*, Float_t*, Int_t) 
{
  printf( " Gdtom.\n" );
}

/***********************************************************************/

void  AliGEODB::Glmoth(const char*, Int_t, Int_t&, Int_t*, Int_t*) 
{
  printf( " Glmoth.\n" );

}

/***********************************************************************/

void  AliGEODB::Gmedia(Float_t*, Int_t&) 
{
  printf( " Gmedia.\n" );
}

/***********************************************************************/

void  AliGEODB::Gmtod(Float_t*, Float_t*, Int_t) 
{
  printf( " Gmtod.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsdvn(const char* name, const char* parentname, Int_t ndiv, Int_t iaxis) 
{
    cout << " Inside Gsdvn " << endl;
    /* Divides one node in Ndiv nodes in the iaxis direction*/

    char* transfname = new char[10];
    char* nodename   = new char[10];

    AliGNode* parent = (AliGNode*) listNodes->FindObject(parentname);

    /* Find the shape that will result from the division*/

    /* BOX */
    if( parent->GetShape()->ClassName() == "AliGBox" ) {
        AliGMaterial* material = new AliGMaterial(parent->GetMaterial());

        Float_t Dx, Dy, Dz;

        switch (iaxis) {
            case 1 : Dx = ((AliGBox*)(parent->GetShape()))->GetX() / ndiv;
	             Dy = ((AliGBox*)(parent->GetShape()))->GetY();
	             Dz = ((AliGBox*)(parent->GetShape()))->GetZ();
                     break;
                    
            case 2 : Dx = ((AliGBox*)(parent->GetShape()))->GetX() ;
                     Dy = ((AliGBox*)(parent->GetShape()))->GetY()/ ndiv;
                     Dz = ((AliGBox*)(parent->GetShape()))->GetZ();
                     break;

            case 3 : Dx = ((AliGBox*)(parent->GetShape()))->GetX() ;
                     Dy = ((AliGBox*)(parent->GetShape()))->GetY();
                     Dz = ((AliGBox*)(parent->GetShape()))->GetZ()/ ndiv;
                     break;
                     
            default: Dx = Dy = Dz = 0.;
        };

        AliGBox* box = new AliGBox("box","box", Dx, Dy, Dz);

        /* Create the nodes son*/

        Float_t x, y, z;
        Text_t* expression;
        
        for( int i=0; i<ndiv; i++ ) {
            //nodename = "Node";
            nodename[4] = '\x0';
            sprintf( nodename, "%s%d", nodename, i );
            transfname = "trans";
            transfname[5] = '\x0';
            sprintf( transfname, "%s%d", transfname, i );

            switch (iaxis) {
                case 'X' :
                case 'x' : x = ndiv *i;
	                   y = 0;
                           z = 0;
	                   break;
                case 'Y' :
                case 'y' : x = 0;
	                   y = ndiv *i;
	                   z = 0;
	                   break;

                case 'Z' :
                case 'z' : x = 0;
	                   y = 0;
	                   z = ndiv *i;
	                   break;
            };

            /* Positioning the sons (only translations) and adding them to the parent*/
    
            sprintf( expression, "%s %f %f %f", "tra", x, y, z );
            AliGTransform* tran = new AliGTransform( transfname, transfname, expression );

            char* nodname = new char[strlen(name)];
            strcpy( nodname, name );
            AliGNode* son = new AliGNode( nodname, i, nodname, box, material );
            if( nodname ) delete [] nodname;
            
            parent->Add( son, tran );
            parent->AddConfig("","");
            listNodes->AddLast(son);
        }
    }
}

/***********************************************************************/

void  AliGEODB::Gsdvn2(const char*, const char*, Int_t, Int_t, Float_t, Int_t) 
{
  printf( " Gsdvn2.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsdvs(const char*, const char*, Float_t, Int_t, Int_t) 
{
  printf( " Gsdvs.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsdvs2(const char*, const char*, Float_t, Int_t, Float_t, Int_t) 
{
  printf( " Gsdvs2.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsdvt(const char*, const char*, Float_t, Int_t, Int_t, Int_t) 
{
  printf( " Gsdvt.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsdvt2(const char *, const char *, Float_t, Int_t, Float_t, Int_t, Int_t) 
{
  printf( " Gsdvt2.\n" );
}

/***********************************************************************/

void  AliGEODB::Gsord(const char*, Int_t) 
{
  // Aqui entra pero dice Federico que Dummy.
}

/***********************************************************************/

void  AliGEODB::Gspos( const char *name, Int_t nr, const char *mother, Float_t x, Float_t y, Float_t z, Int_t irot, const char* konly ) 
{

    if( !strcmp(mother, "ALIC") )
        mother = "ALIC_1";

    AliGNode* son    = NULL;
    AliGNode* parent = NULL;

    char* nodename   = new char[strlen(name)+1];
    strcpy( nodename, name );
    nodename[strlen(name)] = '\x0';

    TObjLink *lnk = listNodes->FirstLink();
    int quit=0;

    while( (lnk) && (!quit) ) {
        TObject *obj = lnk->GetObject();

        if( obj->GetName() ) {
            char* NodeName = new char[strlen(obj->GetName())+1];
            strcpy( NodeName, obj->GetName() );
            NodeName[strlen(obj->GetName())] = '\x0';

            int i;
            for( i=0; NodeName[i]!='_'; i++ );

            if( !strncmp(nodename, NodeName, i) )
                son = (AliGNode*) obj;
            else
                if( !strncmp(mother, NodeName, i) )
                    parent = (AliGNode*) obj;

            if( NodeName ) delete [] NodeName;

            if( (son != NULL) && (parent != NULL) )
                quit = 1;
        }

        lnk = lnk->Next();
    }

    if( son == NULL ) {
        cout << " ERROR in Gspos. Couldn't find node " << nodename << " inside listNodes." << endl;
        exit(1);
    }
    else
        if( parent == NULL) {
            cout << " ERROR in Gspos. Couldn't find node " << mother << " inside listNodes." << endl;
            exit(1);
        }

    AliGTransform* tra = (AliGTransform*)listTransf->At(MatrCount - 1);

    TVector* matrix = tra->GetMatrix();
    char* transname = new char[10];
    sprintf(transname, "%s%d", "tra", MatrCount );

    TVector &mat = *matrix;
    AliGTransform* newtra = new AliGTransform(transname, transname, mat(0), mat(1),  mat(2), mat(4), mat(5), mat(6),mat(8), mat(9), mat(10), x, y, z);

    //printf( "\n Created new AliGTransform( %s, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f );\n", transname, transname, mat(0), mat(1), mat(2), mat(4), mat(5), mat(6), mat(8), mat(9), mat(10), x, y, z );

    listTransf->AddAt(newtra, MatrCount);
    MatrCount++;

    parent->Add(son, newtra);
    //printf( " %s->Add( %s, %s );\n", parent->GetName(), son->GetName(), newtra->GetName() );
}

/***********************************************************************/

void  AliGEODB::Gsposp( const char* nodename, Int_t nr, const char* mother, Float_t x, Float_t y, Float_t z, Int_t irot, const char* konly, Float_t* upar, Int_t np ) 
{
    //printf( "\n" );

    Int_t i;
    
    for( i=0; i<PospCount; i++ )
        if( !strcmp(Posp[i]->name, nodename) )
	    break;

    char* shapename  = Posp[i]->shape;
    Int_t nmed       = Posp[i]->nmed;

    char* fnodename  = new char[strlen(nodename)+1];
    strcpy(fnodename, nodename);
    fnodename[strlen(nodename)] = '\x0';

    char* fmother    = new char[strlen(mother)+1];
    strcpy(fmother, mother);
    fmother[strlen(mother)] = '\x0';
    
    AliGNode* parent = NULL;

    TObjLink *lnk = listNodes->FirstLink();
    int quit=0;

    while( (lnk) && (!quit) ) {
        TObject *obj = lnk->GetObject();

        if( obj->GetName() ) {
            char* NodeName = new char[strlen(obj->GetName())];
            strcpy( NodeName, obj->GetName() );

            int i;
            for( i=0; NodeName[i]!='_'; i++ );

            if( !strncmp(fmother, NodeName, i) )
                    parent = (AliGNode*) obj ;
 
            if( NodeName ) delete [] NodeName;

            if( parent != NULL )
                quit = 1;
        }

        lnk = lnk->Next();
    }

    if( parent == NULL) {
        cout << " ERROR in Gsposp. Couldn't find node " << mother << " inside listNodes." << endl;
        exit(1);
    }

    Int_t nmat = AliMediums[nmed]->fNmat;

    AliGMaterial* MyMaterial = NULL;
    
    for( int j=0; j<listMaterials->GetSize(); j++ )
        if( ((AliGMaterial*) listMaterials->At(j))->GetfImat() == nmat ) {
            MyMaterial = (AliGMaterial*) listMaterials->At(j);
            /*printf( " AliGMaterial( %d, %s, %s, %d, %d, %f, %f, %f, %f, %f, %f, *fUbuf, %d, %f, %f, %f, %f, %f, *fBuf, %d ); already existed.\n", nmat,AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fNbuf, AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]->fNwbuf );*/
            break;
        }
    
    if( MyMaterial == NULL ) {
        MyMaterial = new AliGMaterial( nmat, AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fUbuf, AliMediums[nmed]->fNbuf,  AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]
->fBuf, AliMaterials[nmat]->fNwbuf );

        /*printf( " Created new AliGMaterial( %d, %s, %s, %d, %d, %f, %f, %f, %f, %f, %f, *fUbuf, %d, %f, %f, %f, %f, %f, *fBuf, %d );\n", nmat,AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fNbuf, AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]->fNwbuf );*/

        listMaterials->Add( MyMaterial );
    }

    AliGNode* son = NULL;

    // Different constructor depending on the shape
    if( !strncmp(shapename, "TUBE", 4) ) {

        AliGTube* Shape = NULL;
        
        for( int j=0; j<listShapes->GetSize(); j++ )
            if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                Shape = (AliGTube*) listShapes->At(j);
                //printf( " AliGTube( %s, %s, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2] );
                break;
            }

        if( Shape == NULL ) {
            Shape = new AliGTube( shapename, shapename, upar[0], upar[1], upar[2] );
            //printf( " Created new AliGTube( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
            listShapes->Add(Shape);
        }

        Shape->SetCol(Color);
        son = new AliGNode( fnodename, nr, fnodename, Shape, MyMaterial );
    }
    else
        if( !strncmp(shapename, "BOX", 3) ) {

            AliGBox* Shape = NULL;
        
            for( int j=0; j<listShapes->GetSize(); j++ )
                if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                    Shape = (AliGBox*) listShapes->At(j);
                    //printf( " AliGBox( %s, %s, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2] );
                    break;
                }

            if( Shape == NULL ) {
                Shape = new AliGBox( shapename, shapename, upar[0], upar[1], upar[2] );
                //printf( " Created new AliGBox( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
                listShapes->Add(Shape);
            }
            
            Shape->SetCol(Color);
            son = new AliGNode( fnodename, nr, fnodename, Shape, MyMaterial );
        }
        else
            if( !strncmp(shapename, "TRD1", 4) ) {

                AliGTRD1* Shape = NULL;
        
                for( int j=0; j<listShapes->GetSize(); j++ )
                    if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                        Shape = (AliGTRD1*) listShapes->At(j);
                        //printf( " AliGTRD1( %s, %s, %f, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2], upar[3] );
                        break;
                    }

                if( Shape == NULL ) {
                    Shape = new AliGTRD1( shapename, shapename, upar[0], upar[1], upar[2], upar[3] );
                    //printf( " Created new AliGTRD1( %s, %s, %f, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2], upar[3] );
                    listShapes->Add(Shape);
                }

                Shape->SetCol(Color);
                son = new AliGNode( fnodename, nr, fnodename, Shape, MyMaterial );
            }
            else
                if( !strncmp(shapename, "PCON", 4) ) {

                    AliGPCone* Shape = NULL;
        
                    for( int j=0; j<listShapes->GetSize(); j++ )
                        if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                            Shape = (AliGPCone*) listShapes->At(j);
                            //printf( " AliGPCone( %s, %s, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2] );
                            break;
                        }

                    if( Shape == NULL ) {
                        Shape = new AliGPCone( shapename, shapename, upar, np );
                        //printf( " Created new AliGPCone( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
                        listShapes->Add(Shape);
                    }

                    Shape->SetCol(Color);
                    son = new AliGNode( fnodename, nr, fnodename, Shape, MyMaterial );
                }

    Color++;
    if( Color == 50 )
        Color = 0;

    AliGTransform* tra = (AliGTransform*) listTransf->At(0);

    TVector* matrix = tra->GetMatrix();

    char* transname = new char[10];
    sprintf( transname, "%s%d", "tra", MatrCount );
    
    TVector &mat = *matrix;
    AliGTransform* newtra = new AliGTransform(transname, transname, mat(0), mat(1),  mat(2), mat(4), mat(5), mat(6), mat(8), mat(9), mat(10), x, y, z);

    //printf( " Created new AliGTransform( %s, %s, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f );\n", transname, transname, mat(0), mat(1), mat(2), mat(4), mat(5), mat(6), mat(8), mat(9), mat(10), x, y, z );

    listTransf->AddAt(newtra, MatrCount);
    MatrCount++;

    parent->Add(son, newtra);

    //printf( " %s->Add( %s, %s );\n", parent->GetName(), son->GetName(), newtra->GetName() );

    listNodes->Add(son);
}

/***********************************************************************/

void  AliGEODB::Gsrotm(Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t) 
{
  printf( " Gsrotm.\n" );
}

/***********************************************************************/

void  AliGEODB::Gprotm(Int_t)
{
  printf( " Gprotm.\n" );
}

/***********************************************************************/

Int_t AliGEODB::Gsvolu(const char* vname, const char* shape, Int_t nmed, Float_t* upar, Int_t np)
{
    Int_t lenvname = strlen(vname);
    
    char* nodename = new char[lenvname+1];
    strcpy( nodename, vname );
    nodename[lenvname] = '\x0';

    int pos;
    for( pos=0; pos<lenvname; pos++ )
        if( (nodename[pos] == '\x0') || (nodename[pos] == ' ') )
            break;

    if( pos < lenvname ) {
        if( nodename ) delete [] nodename;
        nodename = new char[pos+1];
        strncpy( nodename, vname, pos );
        nodename[pos] = '\x0';
    }

    char* shapename = new char[strlen(shape)+1];
    strcpy( shapename, shape );
    shapename[strlen(shape)] = '\x0';

    //printf( "\n Trying to build node: %s with shape: %s, nmed: %d, upar[0]: %f and np: %d\n", nodename, shapename, nmed, upar[0], np );

    if( !np ) {
        /* if number of parameters is 0, stores the nodename, shapename and nmed to be used in GSPOSP */
        NodePosp* AliPosp = new NodePosp;
	AliPosp->name     = nodename;
	AliPosp->shape    = shapename;
	AliPosp->nmed     = nmed;
        Posp[PospCount]   = AliPosp;
        PospCount++;
        return 0;
    }

    if( (listNodes) && (listNodes->FindObject(nodename)) ) { // Trying to find the node inside listNodes
        // The node existed in listNodes.
        printf( "\n ERROR!!! Node %s already existed in listNodes\n", nodename );
        exit(1);
    }
    else {
        // It's a new node.
        Int_t nmat = AliMediums[nmed]->fNmat;

        AliGMaterial* MyMaterial = NULL;
    
        for( int j=0; j<listMaterials->GetSize(); j++ )
            if( ((AliGMaterial*) listMaterials->At(j))->GetfImat() == nmat ) {
                MyMaterial = (AliGMaterial*) listMaterials->At(j);
                /*printf( " AliGMaterial( %d, %s, %s, %d, %d, %f, %f, %f, %f, %f, %f, *fUbuf, %d, %f, %f, %f, %f, %f, *fBuf, %d ); already existed.\n", nmat,AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fNbuf, AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]->fNwbuf );*/
                break;
            }
    
        if( MyMaterial == NULL ) {
            MyMaterial = new AliGMaterial( nmat, AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fUbuf, AliMediums[nmed]->fNbuf,  AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]->fBuf, AliMaterials[nmat]->fNwbuf );

            listMaterials->Add( MyMaterial );

            /*printf( " Created new AliGMaterial( %d, %s, %s, %d, %d, %f, %f, %f, %f, %f, %f, *fUbuf, %d, %f, %f, %f, %f, %f, *fBuf, %d );\n", nmat,AliMaterials[nmat]->fName, AliMaterials[nmat]->fName, AliMediums[nmed]->fIsvol, AliMediums[nmed]->fIfield, AliMediums[nmed]->fFieldm, AliMediums[nmed]->fTmaxfd, AliMediums[nmed]->fStemax, AliMediums[nmed]->fDeemax, AliMediums[nmed]->fEpsil, AliMediums[nmed]->fStmin, AliMediums[nmed]->fNbuf, AliMaterials[nmat]->fA, AliMaterials[nmat]->fZ, AliMaterials[nmat]->fDens, AliMaterials[nmat]->fRadl, AliMaterials[nmat]->fAbsl, AliMaterials[nmat]->fNwbuf );*/
        }

        // Different constructor depending on the shape
        AliGNode* node = NULL;

        if( !strncmp(shapename, "TUBE", 4) ) {

            AliGTube* Shape = NULL;
        
            for( int j=0; j<listShapes->GetSize(); j++ )
                if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                    Shape = (AliGTube*) listShapes->At(j);
                    printf( " AliGTube( %s, %s, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2] );
                    break;
                }

            if( Shape == NULL ) {
                Shape = new AliGTube( shapename, shapename, upar[0], upar[1], upar[2] );
                //printf( " Created new AliGTube( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
                listShapes->Add(Shape);
            }

            Shape->SetCol(Color);
            node = new AliGNode( nodename, 1, nodename, Shape, MyMaterial );
        }
        else
            if( !strncmp(shapename, "BOX", 3) ) {

                AliGBox* Shape = NULL;
        
                for( int j=0; j<listShapes->GetSize(); j++ )
                    if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                        Shape = (AliGBox*) listShapes->At(j);
                        printf( " AliGBox( %s, %s, %f, %f, %f ); already existed.\n", shapename, shapename, upar[0], upar[1], upar[2] );
                        break;
                    }

                if( Shape == NULL ) {
                    Shape = new AliGBox( shapename, shapename, upar[0], upar[1], upar[2] );
                    //printf( " Created new AliGBox( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
                    listShapes->Add(Shape);
                }

                Shape->SetCol(Color);
                node = new AliGNode( nodename, 1, nodename, Shape, MyMaterial );
            }
            else
                if( !strncmp(shapename, "TRD1", 4) ) {

                    AliGTRD1* Shape = NULL;

                    for( int j=0; j<listShapes->GetSize(); j++ )
                        if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                            Shape = (AliGTRD1*) listShapes->At(j);

                            break;
                        }

                    if( Shape == NULL ) {
                        Shape = new AliGTRD1( shapename, shapename, upar[0], upar[1], upar[2], upar[3] );
                        //printf( " Created new AliGTRD1( %s, %s, %f, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2], upar[3] );
                        listShapes->Add(Shape);
                    }

                    Shape->SetCol(Color);
                    node = new AliGNode( nodename, 1, nodename, Shape, MyMaterial );
                }
                else
                    if( !strncmp(shapename, "PCON", 4) ) {

                        AliGPCone* Shape = NULL;

                        for( int j=0; j<listShapes->GetSize(); j++ )
                            if( ((AliGShape*) listShapes->At(j))->GetName() == shapename ) {
                                Shape = (AliGPCone*) listShapes->At(j);
                                break;
                            }

                        if( Shape == NULL ) {
                            Shape = new AliGPCone( shapename, shapename, upar, np );
                            //printf( " Created new AliGPCone( %s, %s, %f, %f, %f );\n", shapename, shapename, upar[0], upar[1], upar[2] );
                            listShapes->Add(Shape);
                        }

                        Shape->SetCol(Color);
                        node = new AliGNode( nodename, 1, nodename, Shape, MyMaterial );
                    }

        //printf( " Created new AliGNode( %s, %d, %s, %s, %s );\n", nodename, 1+listNodes->GetSize(), nodename, shapename, MyMaterial->GetName() );

        listNodes->Add(node);

        //printf( " node->fShape->Color = %d\n", node->GetShape()->GetCol() );

        Color++;
        if( Color == 50 )
            Color = 1;
    }

    return 0;
}

/***********************************************************************/

void  AliGEODB::Gsatt(const char*, const char*, Int_t) 
{
  printf( " Gsatt.\n" );
}

/***********************************************************************/

void  AliGEODB::Gfpara(const char*, Int_t, Int_t, Int_t&, Int_t&, Float_t*, Float_t*) 
{
  printf( " Gfpara.\n" );
}

/***********************************************************************/

void  AliGEODB::Gckpar(Int_t, Int_t, Float_t*) 
{
  printf( " Gckpar.\n" );
}

/***********************************************************************/

void  AliGEODB::Gckmat(Int_t, char*) 
{
  printf( " Gckmat.\n" );
}

/***********************************************************************/

//=======================DRAW functions

void  AliGEODB::InitHIGZ() 
{
    printf( "InitHIGZ.\n" );
}

/***********************************************************************/

void  AliGEODB::Gdopen(Int_t) 
{
    printf( "Gdopen.\n" );
}

/***********************************************************************/

void  AliGEODB::Gdclose() 
{
    printf( "Gdclose.\n" );
}

/***********************************************************************/

void  AliGEODB::Gdelete(Int_t) 
{
    cout << " Gdelete" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdshow(Int_t) 
{
    cout << " Gdshow" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdopt(const char *,const char *) 
{
    cout << " Gdopt" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdraw(const char *,Float_t, Float_t, Float_t,Float_t,Float_t,Float_t,Float_t) 
{
    cout << " Gdraw" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdrawc(const char *,Int_t, Float_t,Float_t,Float_t,Float_t,Float_t) 
{
    cout << " Gdrawc" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdrawx(const char *,Float_t, Float_t, Float_t, Float_t, Float_t,Float_t,Float_t,Float_t,Float_t) 
{
    cout << " Gdrawx" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdhead(Int_t, const char *, Float_t) 
{
    cout << " Gdhead" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdman(Float_t, Float_t, const char *) 
{
    cout << " Gdman" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdspec(const char *) 
{
    cout << " Gdspec" << endl;
}

/***********************************************************************/

void  AliGEODB::DrawOneSpec(const char *) 
{
    cout << " DrawOneSpec" << endl;
}

/***********************************************************************/

void  AliGEODB::Gdtree(const char *,Int_t,Int_t) 
{
    cout << " Gdtree" << endl;
}

/***********************************************************************/

void  AliGEODB::GdtreeParent(const char *,Int_t,Int_t) 
{
    cout << " GdtreeParent" << endl;
}

/***********************************************************************/

//=======================Set functions

void AliGEODB::SetABAN(Int_t)
{
    cout << " SetABAN" << endl;
}

/***********************************************************************/

void AliGEODB::SetANNI(Int_t)
{
    cout << " SetANNI" << endl;
}

/***********************************************************************/

void AliGEODB::SetAUTO(Int_t)
{
    cout << " SetAUTO" << endl;
}

/***********************************************************************/

void AliGEODB::SetBOMB(Float_t)
{
    cout << " SetBOMB" << endl;
}

/***********************************************************************/

void AliGEODB::SetBREM(Int_t)
{
    cout << " SetBREM" << endl;
}

/***********************************************************************/

void AliGEODB::SetCKOV(Int_t)
{
    cout << " SetCKOV" << endl;
}

/***********************************************************************/

void AliGEODB::SetClipBox(const char *,Float_t,Float_t, Float_t,Float_t,Float_t,Float_t) 
{
    cout << " SetClipBox" << endl;
}

/***********************************************************************/

void AliGEODB::SetCOMP(Int_t)
{
    cout << " SetCOMP" << endl;
}

/***********************************************************************/

void AliGEODB::SetCUTS( Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t )
{
    cout << " SetCUTS" << endl;
}

/***********************************************************************/

void AliGEODB::SetDCAY(Int_t)
{
    cout << " SetDCAY" << endl;
}

/***********************************************************************/

void AliGEODB::SetDEBU(Int_t, Int_t, Int_t)
{
    cout << " SetDEBU" << endl;
}

/***********************************************************************/

void AliGEODB::SetDRAY(Int_t)
{
    cout << " SetDRAY" << endl;
}

/***********************************************************************/

void AliGEODB::SetHADR(Int_t)
{
    cout << " SetHADR" << endl;
}

/***********************************************************************/

void AliGEODB::SetKINE(Int_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t, Float_t,Float_t ) 
{
    cout << " SetKINE" << endl;
}

/***********************************************************************/

void AliGEODB::SetLOSS(Int_t)
{
    cout << " SetLOSS" << endl;
}

/***********************************************************************/

void AliGEODB::SetMULS(Int_t)
{
    cout << " SetMULS" << endl;
}

/***********************************************************************/

void AliGEODB::SetMUNU(Int_t)
{
    cout << " SetMUNU" << endl;
}

/***********************************************************************/

void AliGEODB::SetOPTI(Int_t)
{
    cout << " SetOPTI" << endl;
}

/***********************************************************************/

void AliGEODB::SetPAIR(Int_t)
{
    cout << " SetPAIR" << endl;
}

/***********************************************************************/

void AliGEODB::SetPFIS(Int_t)
{
    cout << " SetPFIS" << endl;
}

/***********************************************************************/

void AliGEODB::SetPHOT(Int_t)
{
    cout << " SetPHOT" << endl;
}

/***********************************************************************/

void AliGEODB::SetRAYL(Int_t)
{
    cout << " SetRAYL" << endl;
}

/***********************************************************************/

void AliGEODB::SetSWIT(Int_t , Int_t)
{
    cout << " SetSWIT" << endl;
}

/***********************************************************************/

void AliGEODB::SetTRIG(Int_t)
{
    cout << " SetTRIG" << endl;
}

/***********************************************************************/

void AliGEODB::Vname(const char *, char *) 
{
    cout << " Vname" << endl;
}

/***********************************************************************/

void AliGEODB::InitLego() 
{
    cout << " InitLego" << endl;
}

/***********************************************************************/

extern "C" void sxpart_(){}
