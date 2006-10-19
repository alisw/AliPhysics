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
 
/* $Id$ */
///////////////////////////////////////////////////////////////////
//Class for reconstruction of clusters V2                        //
//                                                               //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TParticle.h>

#include "AliRun.h"
#include "AliHeader.h"

#include "AliITSRecPoint.h"
#include "AliITSFindClustersV2.h"
#include "AliITSclusterV2.h"
#include "AliITSgeom.h"
#include "AliMC.h"

ClassImp(AliITSFindClustersV2)

//______________________________________________________________________
AliITSFindClustersV2::AliITSFindClustersV2():
fAr(0),
fDeletfAr(kFALSE),
fGeom(0),
fInFileName(0),
fOutFileName(0),
fIn(0),
fOut(0),
fInit(kFALSE),
fSlowFast(kFALSE){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A zero-ed constructed AliITSFindClustersV2 class.

}
//______________________________________________________________________
AliITSFindClustersV2::AliITSFindClustersV2(const TString infile,
                                           const TString outfile):
fAr(0),
fDeletfAr(kFALSE),
fGeom(0),
fInFileName(0),
fOutFileName(0),
fIn(0),
fOut(0),
fInit(kFALSE),
fSlowFast(kFALSE){
    // Standard constructor.
    // Inputs:
    //  const TString infile   Input file name where the RecPoints are
    //                         to be read from.
    //  const TString outfile  Output file where V2 tracking clulsters
    //                         are to be written. if =="" writen to the
    //                         same file as input.
    // Outputs:
    //   none.
    // Return:
    //    A standardly constructed AliITSFindClustersV2 class.

    fInFileName = new TString(infile);
    if(outfile.CompareTo("")!=0){
	fOutFileName = new TString(outfile);
    } // end if outfile.CompareeTo("")!=0
    
    if(fOutFileName!=0){
	fIn = new TFile(fInFileName->Data(),"READ");
	fOut = new TFile(fOutFileName->Data(),"UPDATE");
    }else{ // open fIn file for update
	fIn = new TFile(fInFileName->Data(),"UPDATE");
    } // end if

    fAr  = (AliRun*) fIn->Get("gAlice");
    if(!fAr){
	Warning("AliITSFindClusterV2",
		"Can't fine gAlice in file %s. Aborting.",fIn->GetName());
	return;
    } // end if !fAr
    fDeletfAr = kTRUE; // Since gAlice was read in, delete it.
    /*
    AliITS *its = (AliITS*) fAr->GetModule("ITS");
    if(!its){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS in gAlice. Aborting.");
	return;
    } // end if !its
    fGeom = its->GetITSgeom();
    if(!fGeom){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS geometry in gAlice. Aborting.");
	return;
    } // end if !fGeom
    */
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    fGeom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
    if(!fGeom){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS geometry in gAlice. Aborting.");
	return;
    } // end if !fGeom
    delete rl;

    if(fOut) fOut->cd();
    fInit = kTRUE;
}
//______________________________________________________________________
AliITSFindClustersV2::AliITSFindClustersV2(TFile *in,
					   TFile *out):
fAr(0),
fDeletfAr(kFALSE),
fGeom(0),
fInFileName(0),
fOutFileName(0),
fIn(0),
fOut(0),
fInit(kFALSE),
fSlowFast(kFALSE){
    // Standard constructor.
    // Inputs:
    //  const TFile *in   Input file where the RecPoints are
    //                         to be read from.
    //  const Tfile *out  Output file where V2 tracking clulsters
    //                         are to be written. if ==0 writen to the
    //                         same file as input.
    // Outputs:
    //   none.
    // Return:
    //    A standardly constructed AliITSFindClustersV2 class.

    fIn  = in;
    fOut = out;
    fAr  = (AliRun*) fIn->Get("gAlice");
    if(!fAr){
	Warning("AliITSFindClusterV2",
		"Can't fine gAlice in file %s. Aborting.",fIn->GetName());
	return;
    } // end if !fAr
    fDeletfAr = kTRUE; // Since gAlice was read in, delete it.
    /*
    AliITS *its = (AliITS*) fAr->GetModule("ITS");
    if(!its){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS in gAlice. Aborting.");
	return;
    } // end if !its
    fGeom = its->GetITSgeom();
    */
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    fGeom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
    if(!fGeom){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS geometry in gAlice. Aborting.");
	return;
    } // end if !fGeom
    delete rl;
 
    if(fOut) fOut->cd();
    fInit = kTRUE;
}
//______________________________________________________________________
AliITSFindClustersV2::AliITSFindClustersV2(AliRun *ar,
					   const TString outfile):
fAr(0),
fDeletfAr(kFALSE),
fGeom(0),
fInFileName(0),
fOutFileName(0),
fIn(0),
fOut(0),
fInit(kFALSE),
fSlowFast(kFALSE){
    // Standard constructor.
    // Inputs:
    //  const TString infile   Input file name where the RecPoints are
    //                         to be read from.
    //  const TString outfile  Output file where V2 tracking clulsters
    //                         are to be written. if =="" writen to the
    //                         same file as input.
    // Outputs:
    //   none.
    // Return:
    //    A standardly constructed AliITSFindClustersV2 class.

    if(outfile.CompareTo("")!=0){
	fOutFileName = new TString(outfile);
    } // end if outfile.CompareeTo("")!=0
    
    if(fOutFileName!=0){
	fOut = new TFile(fOutFileName->Data(),"UPDATE");
    } // end if

    fAr  = ar;
    if(!fAr){
	Warning("AliITSFindClusterV2",
		"ar==0. Aborting.");
	return;
    } // end if !fAr

    /*   AliITS *its = (AliITS*) fAr->GetModule("ITS");
    if(!its){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS in gAlice. Aborting.");
	return;
    } // end if !its
    fGeom = its->GetITSgeom();
    */
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    rl->CdGAFile();
    fGeom = (AliITSgeom*)gDirectory->Get("AliITSgeom");
    if(!fGeom){
	Warning("AliITSFindClusterV2",
		"Can't fine the ITS geometry in gAlice. Aborting.");
	return;
    } // end if !fGeom
    delete rl;

    if(fOut) fOut->cd();
    fInit = kTRUE;
}
//______________________________________________________________________
AliITSFindClustersV2::AliITSFindClustersV2(const AliITSFindClustersV2 &rec):TTask(rec),
fAr(rec.fAr),
fDeletfAr(rec.fDeletfAr),
fGeom(rec.fGeom),
fInFileName(rec.fInFileName),
fOutFileName(rec.fOutFileName),
fIn(rec.fIn),
fOut(rec.fOut),
fInit(rec.fInit),
fSlowFast(rec.fSlowFast){
    // Copy constructor. 

  
}
//______________________________________________________________________
AliITSFindClustersV2& AliITSFindClustersV2::operator=(const AliITSFindClustersV2& source){
    // Assignment operator. 
  this->~AliITSFindClustersV2();
  new(this) AliITSFindClustersV2(source);
  return *this;
}

//______________________________________________________________________
AliITSFindClustersV2::~AliITSFindClustersV2(){
    // Default constructor.
    // Inputs:
    //  none.
    // Outputs:
    //   none.
    // Return:
    //    A destroyed AliITSFindclustersV2 class.

    fGeom = 0; // Deleted when AliRun/ITS is deleted.
    if(fInFileName!=0){ // input file opened by AliITSFindClustersV2
	if(fIn!=0){
	    if(fIn->IsOpen()) fIn->Close();
	    delete fIn;
	    fIn = 0;
	} // end if
	delete fInFileName;
	fInFileName = 0;
    } // end if

    if(fOutFileName!=0){ // input file opened by AliITSFindClustersV2
	if(fOut!=0){
	    if(fOut->IsOpen()) fOut->Close();
	    delete fOut;
	    fOut = 0;
	} // end if
	delete fOutFileName;
	fOutFileName = 0;
    } // end if
    if(fDeletfAr && !fAr){
	cout << "deleting fAr."<< endl;
	delete fAr;
	fAr = 0;
	cout << "fAr deleted OK."<< endl;
    } // end if fDeletfAr
}
//______________________________________________________________________ 
void AliITSFindClustersV2::Exec(const Option_t *opt){
    // Main FindclustersV2 function.
    // Inputs:
    //      Option_t * opt   list of subdetector to digitize. =0 all.
    // Outputs:
    //      none.
    // Return:
    //      none.
    Char_t name[50];

    if(!fInit){
	Warning("Exec","(opt=%s) Initilization not succesfull. Aborting.\n",opt);
	return;
    } // end if !fInit

    fGeom->Write();

    fAr->GetEvent(0);
    TTree *pTree = fAr->TreeR();
    if(!pTree){
	Warning("Exec","Error getting TreeR. TreeR=%p",pTree);
	return;
    } // end if !pTree
    TBranch *branch = 0;
    if(fSlowFast) sprintf(name,"ITSRecPointsF");
    else sprintf(name,"ITSRecPoints");
    branch = pTree->GetBranch(name);
    if(!branch){
	Warning("Exec","Can't find branch %s in TreeR fSlowFast=%d",
		name,fSlowFast);
	return;
    } // end if !branch
    TClonesArray *points = new TClonesArray("AliITSRecPoint",10000);
    branch->SetAddress(&points);
    Int_t nentr = (Int_t) branch->GetEntries();

    if(fOut) fOut->cd();
    TClonesArray *cluster = new TClonesArray("AliITSclusterV2",10000);
    sprintf(name,"TreeC_ITS_%d",fAr->GetHeader()->GetEvent());
    TTree *cTree = new TTree(name,"ITS clusters");
    cTree->Branch("Clusters",&cluster);
    TClonesArray &cl = *cluster;

    Float_t lp[5];
    Int_t lab[6];
    Int_t i,j,lay,lad,det,nclusters=0,ncl;
    Float_t kmip,x,y,zshift,yshift;
    Double_t rot[9];
    AliITSRecPoint *p;
    TParticle *part;

    for(i=0;i<nentr;i++){
	points->Clear();
	branch->GetEvent(i);
	ncl = points->GetEntriesFast();
	if(ncl<=0) {
	    cTree->Fill();
	    continue;
	} // end if ncl<=0
	nclusters += ncl;
	fGeom->GetModuleId(i,lay,lad,det);
	fGeom->GetTrans(i,x,y,zshift);
	fGeom->GetRotMatrix(i,rot);
	yshift = x*rot[0] + y*rot[1];
	kmip = 1.0; // fGeom->GetModuletype(i)==kSPD
	if(fGeom->GetModuleType(i)==kSDD) kmip = 280.0;
	if(fGeom->GetModuleType(i)==kSSD) kmip = 38.0;
	for(j=0;j<ncl;j++){
	    p = (AliITSRecPoint*)(points->UncheckedAt(j));
	    lp[0] = - (p->GetDetLocalX() + yshift);
	    if(lay==1) lp[0] = -lp[0];
	    lp[1] = p->GetZ() + zshift;
	    lp[2] = p->GetSigmaDetLocX2();
	    lp[3] = p->GetSigmaZ2();
	    lp[4] = p->GetQ()/kmip;
	    lab[0] = p->GetLabel(0);
	    lab[1] = p->GetLabel(1);
	    lab[2] = p->GetLabel(2);
	    lab[3] = i;
	    lad = lab[0];
	    if(lad>=0){
		part = (TParticle*) fAr->GetMCApp()->Particle(lad);
		lad = -3;
		while(part->P() < 0.005){
		    if(part->GetFirstMother()<0){
			Warning("Exec","Primary momentum: %e",part->P());
			break;
		    } // end if part->GetFirstMother()<0
		    lad = part->GetFirstMother();
		    part = (TParticle*) fAr->GetMCApp()->Particle(lad);
		} // end while part->P() < 0.005
		if(lab[1]<0) lab[1] = lad;
		else if(lab[2]<0) lab[2] = lad;
		else Warning("Exec","No empty lables!");
	    } // end if lab>=0
	    Int_t info[3]={0,0,0};
	    new(cl[j]) AliITSclusterV2(lab,lp,info);
	} // end for j
	cTree->Fill();
	cluster->Delete();
	points->Delete();
    } // end for i
    cTree->Write();
    
    delete cTree;
    delete cluster;
    delete points;
}
