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


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  This class is a container class for the laser beam positions of          //
//  the TPC laser system (ALICE-INT-2002-022  v.1 ).                         //
//  The system consits of 6 Laser Rods, each contains four mirror bundels    //
//  on each side of the TPC. Each bundle has seven mirrors.                  //
//                                                                           //
//  The TPC side 0 corresponds to the "Shaft Side (A)", side 1 corresponds   //
//  to the "Muon Side (C)".
//  The laser rods are counted counter clockwise starting with the one       //
//  with the smalles phi angle.                                              //
//  The mirror bundles in one rod for each side are counted starting from    //
//  form the greatest z value to the smalles.                                //
//  The beams in each bundel are counted counter clockwise. Having all beams //
//  pointing downward, the first beam will be the leftmost.                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*
 Create File with design positions

 .L AliTPCLaserTrack.cxx+
 AliTPCLaserTracks ltracks
 ltracks.WriteTreeDesignData()
 TFile f("LaserTracksDesing.root");
 TTree * tree =(TTree*)f.Get("LaserTracks");
 
*/

#include <iostream.h>
#include <TString.h>
#include <TPolyLine3D.h>
#include <TMath.h>
#include <TTree.h>
#include <TFile.h>
#include <TEventList.h>
#include <TVector2.h>
#include <TVector3.h>
#include <TROOT.h>

#include <TGraph.h>


#include "AliTPCLaserTracks.h"

////////////////////////////////////////////////////////////////////////
//              Class AliTPCLaserTracks
////////////////////////////////////////////////////////////////////////

ClassImp(AliTPCLaserTracks)

AliTPCLaserTracks::AliTPCLaserTracks():
    fId(-1),
    fSide(-1),
    fRod(-1),
    fBundle(-1),
    fBeam(-1),
    fX(0),
    fY(0),
    fZ(0),
    fPhi(0),
    fTheta(0),
    fMaxSize(0),
    fNpoints(0),
    fXarr(0x0),
    fYarr(0x0),
    fZarr(0x0)

{
    //
    //  AliTPCLaserTracks default constructor
    //
}

//_______________________________________________________________________
AliTPCLaserTracks::AliTPCLaserTracks(Int_t npoints):
    fId(-1),
    fSide(-1),
    fRod(-1),
    fBundle(-1),
    fBeam(-1),
    fX(0),
    fY(0),
    fZ(0),
    fPhi(0),
    fTheta(0),
    fMaxSize(0),
    fNpoints(0),
    fXarr(0x0),
    fYarr(0x0),
    fZarr(0x0)

{
    //
    //  AliTPCLaserTracks constructor to initialise array of points
    //
    fNpoints = npoints;
    InitPoints();
}

//_______________________________________________________________________
Double_t AliTPCLaserTracks::FindBeamLength(TVector3 vF, TVector3 vP)
{
    //
    //  Find distance between mirror and the intersection between
    //  inner and outer field cage respectively
    //
    //  use the pq formula to find the intersection point between
    //  the beam and the inner and outer fieldcage cylinders
    //  The formulae solved are
    //
    //                        line = zylinder
    //  xPoint = vF.X() + n*vP.X() = r*cos(phi)
    //  yPoint = vF.Y() + n*vP.Y() = r*sin(phi)
    //  zPoint = vF.Y() + n*vP.Y()
    //
    //  get n from the first two equations, calculate x,y,zPoint
    //
    //  vF is the mirror position and vP the beam direction
    Double_t r[2];  //smalles r[0] and largest r[1] pad radius
    r[0] = 80.0;  //smalles pad radius
    r[1] = 260.0;   //largest pad radius

    Double_t fxpxfypy = vF.X()*vP.X()+vF.Y()*vP.Y();
    Double_t px2py2   = vP.X()*vP.X()+vP.Y()*vP.Y();

    for (Int_t i=0; i<2; i++){
	Double_t rad = (fxpxfypy/px2py2)*(fxpxfypy/px2py2)-
	              (vF.X()*vF.X()+vF.Y()*vF.Y()-r[i]*r[i])/px2py2;

	if ( rad >= 0 ){
	    Double_t n1   = -(fxpxfypy)/px2py2+TMath::Sqrt(rad);
	    TVector3 vI1(vF.X()+n1*vP.X(),
			 vF.Y()+n1*vP.Y(),
			 vF.Z()+n1*vP.Z());

	    Double_t n2   = -(fxpxfypy)/px2py2-TMath::Sqrt(rad);
	    TVector3 vI2(vF.X()+n2*vP.X(),
			 vF.Y()+n2*vP.Y(),
			 vF.Z()+n2*vP.Z());


	    //if we cross two boarders on our way return the closer boarder
            if ( n1>0 && n2>0 )
		if ( (vF-vI1).Mag() <= (vF-vI2).Mag() )
		    return (vF-vI1).Mag();
		else
		    return (vF-vI2).Mag();

	    if ( n1>0 )
		return (vF-vI1).Mag();

	    if ( n2>0 )
		return (vF-vI2).Mag();

	}
    }

    return 3.4e38;  //not found!!
}

//_______________________________________________________________________
void AliTPCLaserTracks::WriteTreeDesignData()
{
    //
    //   Write a tree with the design data of the laser track positions
    //

    TFile *f = TFile::Open("LaserTracksDesing.root","recreate");

    AliTPCLaserTracks *ltp = new AliTPCLaserTracks(2);

    TTree *t = new TTree("LaserTracks","LaserTracks");
    t->Branch("LaserTracks","AliTPCLaserTracks",&ltp);
    t->AutoSave();

    Double_t phiBeam[7]         = {-31.8,-16.,-9.2,2.5,9.2,16.,31.8};
    Double_t phiRod[6]          = {40.,100.,160.,220.,280.,340.};     //-20. for the muon side
//    Double_t zBundleCorse[2][4] = {{10,79,163,241},{13,85,169,247}};
    Double_t zBundleCorse[2][4] = {{241.,163.,79.,10.},{247.,169.,85.,13.}};
    Double_t zBeamFine[7]       = {1.45,1.3,1.15,1.3,1.15,1.3,1.45};

    Int_t id=0;
    //loop over TPC side -- 0:A (Shaft) side -- 1:C (Muon) side
    for (Int_t side=0; side<2; side++){
	//loop over laser rods -- counterclockwise
	for (Int_t rod=0; rod<6; rod++){
            //center of the rod
	    TVector2 vRod;

            Double_t phiRod2 = phiRod[rod]-side*20;
	    vRod.SetMagPhi(254.25,                   //center of the rod at 254.25cm
			   TMath::DegToRad()*
			   (phiRod2)  //-20 deg on C-Side
			  );

	    //loop over bundle -- counted from large z to small z
	    for (Int_t bundle=0; bundle<4; bundle++){
		//center of the bundle; not yet rotated
                Int_t bundleS = bundle;
		//		if ( side == 1 ) bundleS = 4-bundle;
		if ( side == 1 ) bundleS = 3-bundle;

		TVector2 vBundle;
		vBundle.SetMagPhi(.65,
				  TMath::DegToRad()*
				//?  TMath::Power(-1,side)*
				  (phiRod2+180.-90.*bundleS)
				 );

		//loop over beam
                Int_t i=0;
		for (Int_t beam=0; beam<7; beam++){
		    TVector2 vBeam;
		    if ( beam == 3 )
			vBeam.Set(0.,0.);
		    else{
			vBeam.SetMagPhi(1.,TMath::DegToRad()*(phiRod2+i*60));
                        i++;
		    }

		    TVector2 xyMirror = vRod+vBundle+vBeam;
		    Double_t zMirror   = (zBundleCorse[rod%2][bundleS]+zBeamFine[beam])*TMath::Power(-1,side);
                    TVector3 v3Mirror(xyMirror.X(),xyMirror.Y(),zMirror);

		    Double_t phi      = TMath::DegToRad()*(phiRod2+180+phiBeam[beam]*TMath::Power(-1,side));
                    Double_t theta    = TMath::DegToRad()*90;

		    TVector3 v3Beam;
		    v3Beam.SetMagThetaPhi(1,theta,phi);  //Direction Vector
                    v3Beam.SetMagThetaPhi(FindBeamLength(v3Mirror,v3Beam), theta, phi);



		    ltp->SetId(id);
		    ltp->SetSide(side);
		    ltp->SetRod(rod);
		    ltp->SetBundle(bundleS);
		    ltp->SetBeam(beam);

		    ltp->SetX(xyMirror.X());
		    ltp->SetY(xyMirror.Y());
		    ltp->SetZ(zMirror);

		    ltp->SetPhi(phi);
		    ltp->SetTheta(theta);

//		    ltp->InitPoints(2);
		    ltp->SetPoint(0,xyMirror.X(),xyMirror.Y(),zMirror);
		    ltp->SetPoint(1,xyMirror.X()+v3Beam.X(),xyMirror.Y()+v3Beam.Y(),zMirror+v3Beam.Z());

		    cout<< "Id: " << id
			<< "  Side: " << side
			<< "  Rod: " << rod
			<< "  Bundle: " << bundleS
			<< "  Beam: " << beam
			<< endl;


		    t->Fill();
                    //delete line;
                    id++;
//		    cout << "Filled!!!" << endl;
		}
	    }
	}
    }

    t->Write();
    delete f;
//    delete t;
//    delete ltp;
//    delete line;
}

//_______________________________________________________________________
TPolyLine3D* AliTPCLaserTracks::GetLine()
{
   if ( fNpoints ) return new TPolyLine3D(fNpoints,fXarr,fYarr,fZarr);

    return 0x0;
}

//_______________________________________________________________________
TObjArray* AliTPCLaserTracks::GetLines(Char_t* file, Char_t *cuts)
{

    TObjArray *array = new TObjArray;

    TFile *f = TFile::Open(file,"read");
    TTree *tree = (TTree*)f->Get("LaserTracks");

    AliTPCLaserTracks *ltp = new AliTPCLaserTracks();
//    TEventList *evList = new TEventList("evList");
    TEventList evList("evList");

    tree->SetBranchAddress("LaserTracks",&ltp);

    tree->Draw(">>evList",cuts,"goff");

//    cout << "N: " << evList.GetN() << endl;
    gROOT->cd();
    for (Int_t ev = 0; ev < evList.GetN(); ev++){
//        cout << ev << endl;
	tree->GetEntry( evList.GetEntry(ev) );
        array->Add(ltp->GetLine());
    }
//    cout << "delte f"<< endl;
//    delete f;

//    cout << "evlist" << endl;

    return array;
}

//_______________________________________________________________________
void AliTPCLaserTracks::InitPoints()
{
    //
    //  Init point arrays
    //

    fMaxSize = fNpoints;
    fXarr = new Double_t[fMaxSize];
    fYarr = new Double_t[fMaxSize];
    fZarr = new Double_t[fMaxSize];
}

//_______________________________________________________________________
Int_t AliTPCLaserTracks::SetPoint(Int_t point, Double_t x, Double_t y, Double_t z)
{
    //
    // Set point to point arrays
    //

    if ( !fXarr || !fYarr || !fZarr ) return 1;
    if ( point > fNpoints-1 ) return 1;

    fXarr[point] = x;
    fYarr[point] = y;
    fZarr[point] = z;
    return 0;
}

//_______________________________________________________________________
Int_t AliTPCLaserTracks::FindMirror(Char_t *file, Double_t x, Double_t y, Double_t z, Double_t phi)
{
    //
    //  return the id of the mirror in 'file' that maches best the given parameters
    //

    TFile *f = TFile::Open(file,"read");
    TTree *tree = (TTree*)f->Get("LaserTracks");

    AliTPCLaserTracks *ltp = new AliTPCLaserTracks();
    TEventList evList("evList");

    tree->SetBranchAddress("LaserTracks",&ltp);

    TString s;
    s = "abs(fX-"; s+= x;
    s+= ")<3&&abs(fY-"; s+=y;
    s+= ")<3&&abs(fZ-"; s+= z;
    s+= ")<1.5";
//    s+= "&&((abs(fPhi-";s+= phi;
//    s+= ")<.06)||(abs(fPhi-";s+= (phi-TMath::DegToRad()*180);
//    s+= ")<.06))";
    cout << s.Data() << endl;

    tree->Draw(">>evList",s.Data());

    Double_t dphiMin=TMath::Pi();
    Int_t id=-1;

    cout << "nev: " << evList.GetN() << endl;
    for (Int_t i=0; i<evList.GetN(); i++){
	tree->GetEntry( evList.GetEntry(i) );
	Double_t dphi = TMath::Abs(ltp->GetPhi() - phi);
	if ( dphi > TMath::DegToRad()*180 ) dphi-=TMath::DegToRad()*180;
	if ( dphi<dphiMin ) {
            dphiMin = dphi;
	    id = ltp->GetId();
	}
    }

    return id;
    delete f;
}

//_______________________________________________________________________
AliTPCLaserTracks::~AliTPCLaserTracks()
{
    //
    //  standard destructor
    //

//    if ( fP ) delete fP;
    if ( fXarr ) delete fXarr;
    if ( fYarr ) delete fYarr;
    if ( fZarr ) delete fZarr;
}
