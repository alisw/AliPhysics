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
  $Id:
 */

#include <TObject.h>
#include <TObjArray.h>

#include <TGeoMaterial.h>
#include <TGeoMatrix.h>
#include <TGeoMedium.h>
#include <TGeoPgon.h>
#include "AliITSBaseGeometry.h"


ClassImp(AliITSMixture)

AliITSMixture::AliITSMixture(const char *name,Int_t N,Double_t *w,TObjArray *m,
			     Double_t rho,Double_t radlen,Double_t intleng)
    :TGeoMixture(name,1,rho){
    // Defines a new mixture from a number of Mixtures, and put the
    // resulting mixture into this object. This will compute avarage
    // isotopic value between different elements.
    // Inputs:
    //   Int_t    N   The number of mixtures in te TObjArray
    //   Double_t *w  The array of weights of each mixture
    //   TObjArray *m The array of AliITSMixture (TGeoMixture)s 
    //                to be mixed.
    // Output:
    //   none.
    // Return:
    //   none.
    Int_t i,z=0,j,Nel;
    Double_t tw,*nw,wel[110],Ael[110],el[110];
    TGeoMixture *mix;

    if(N>m->GetEntries()){ // Error not enough mixtures defined
	Error("Mixing","There are more weight defined than mixtures");
	return;
    } // end if
    // First normilize the weights just in case.
    tw = 0.0;
    for(i=0;i<N;i++) if(w[i]>0.0) tw += w[i];
    nw = new Double_t[N];
    for(i=0;i<N;i++) {if(w[i]>0.0) nw[i] = w[i]/tw;else nw[i] = 0.0;}
    //
    Nel=0;
    for(i=0;i<110;i++) {el[i] = wel[i] = Ael[i] = 0.0;}
    for(i=0;i<N;i++)if(w[i]>0.0) {
	mix = (TGeoMixture*) (m->At(i));
	for(j=0;j<mix->GetNelements();j++) {
	    z = (Int_t) ((mix->GetZmixt())[j]);
	    wel[z] += nw[i]*((mix->GetWmixt())[j]);
	    el[z]  += wel[z]*((mix->GetZmixt())[j]);
	    Ael[z] += wel[z]*((mix->GetAmixt())[j]);
	} // end for j
    } // end for i
    tw = 0.0;
    for(i=1;i<110;i++) if(wel[i]>0.0){
	Nel++;
	tw += wel[i];
    } // end for
    if(tw<=0.0) {  // Error no elements defined.
	Error("Mixing","Total weight of this mixture is zero");
	delete[] nw;
	return;
    } // end if
    // setup TGeoMixture data members.
    fNelements = Nel;
    if(fZmixture!=0) delete[] fZmixture;
    if(fAmixture!=0) delete[] fAmixture;
    if(fWeights!=0)  delete[] fWeights;
    fZmixture = new Double_t[Nel];
    fAmixture = new Double_t[Nel];
    fWeights  = new Double_t[Nel];
    if(rho>0.) fDensity = rho;
    else { // try to compute density form mixture.
	rho = 0.0;
	for(i=0;i<N;i++) if(nw[i]>0.0) {
	    mix = (TGeoMixture*) (m->At(i));
	    rho += nw[i]*(mix->GetDensity());
	} // end for i
	fDensity = rho;
    } // end if
    if(radlen>0.) fRadLen = radlen;
    else { // try to compute radiation form mixture.
	// From "Review of Particle Physics" Particle Data Group Section
	// 26.4.1 equation 26.21 (2002).
	radlen = 0.0;
	for(i=0;i<N;i++) if(nw[i]>0.0) {
	    mix = (TGeoMixture*) (m->At(i));
	    if(mix->GetRadLen()>0.0) rho += 1.0/(nw[i]*(mix->GetRadLen()));
	} // end for i
	fRadLen = 1.0/radlen;
    } // end if
    if(intleng>0.) fIntLen = intleng;
    else { // try to compute interaction form mixture.
	intleng = 0.0;
	for(i=0;i<N;i++) if(nw[i]>0.0) {
	    mix = (TGeoMixture*) (m->At(i));
	    if(mix->GetIntLen()>0.0) intleng += 1.0/(nw[i]*(mix->GetIntLen()));
	} // end for i
	fIntLen = 1.0/intleng;
    } // end if
    j = 0;
    for(z=1;z<110;z++){
	wel[z] /= tw;
	el[z]  /= tw;
	Ael[z] /= tw;
	if(wel[z]>0.0) this->DefineElement(j++,Ael[z],el[z],wel[z]);
    } // end for i
    delete[] nw;
}
//======================================================================
ClassImp(AliITSGeoCable)
;
AliITSGeoCable::AliITSGeoCable(){
    //

    fRmin = fRmax = 0.0;
    fNs.SetXYZ(0.0,0.0,0.0);
    fNe.SetXYZ(0.0,0.0,0.0);
    fTubes    = 0;
    fTranRot  = 0;
}
//----------------------------------------------------------------------
AliITSGeoCable::AliITSGeoCable(const char *name,const TObjArray *vect,
                               const Double_t Rmin,const Double_t Rmax,
                               const TVector3 ns,const TVector3 ne){
    //
    // Inputs:
    //    char      *name   Name of this compound object
    //    TObjArray *vect   Array of TVector3's of points representing the
    //                      path of the cable
    //    TVector3  ns=0    Normal vector representing the angle of the 
    //                      starting surface, default perpendicular
    //    TVector3  ne=0    Normal vector representing the angle of the 
    //                      ending surface, default perpendicular
    // Outputs:
    //    none.
    // Return:
    //    A fully initilized and created AliITSGeoCable class.
    Char_t nam[500];
    Int_t i,n;
    Double_t s,th,ph;
    TVector3 x0,x1,x2,d,t,n0,n1;

    fRmin = Rmin;
    fRmax = Rmax;
    fRmin = fRmax = 0.0;
    fNs.SetXYZ(0.0,0.0,-1.0);
    fNe.SetXYZ(0.0,0.0,1.0);
    n = vect->GetEntries();
    fTubes = new TObjArray(n-1);
    fTranRot  = new TObjArray(n-1);
    fNs = ns*(1./ns.Mag());
    fNe = ne*(1./ne.Mag());
    //
    x0 = *((TVector3 *)(vect->At(0)));
    n0 = ns;
    for(i=1;i<n;i++){
        x1 = *((TVector3 *)(vect->At(i)));
        d  = x1 - x0;
        if(i<n-1) {
            x2 = *((TVector3 *)(vect->At(i+1)));
            n1 = d + (x2-x1);
            n1 *= 1./n1.Mag();
        }else{
            n1 = fNe;
        } // end if
        t  = 0.5*(x1 + x0);
        th =  TMath::ATan2(TMath::Sqrt(2.*d.Mag2()-d.z()-2.*d.Mag()*d.y()-
                                        2.*d.Mag()*d.x()),
                           TMath::Sqrt(d.z()*d.z()-2.*d.Mag()*d.z()+d.Mag2()));
        th *= TMath::RadToDeg();
        ph =  TMath::ATan2(d.y()-d.Mag(),d.x()-d.Mag());
        ph *= TMath::RadToDeg();
        sprintf(nam,"%sCombiTrans%dCable",name,i-1);
        fTranRot->AddAt(new TGeoCombiTrans(nam,t.x(),t.y(),t.z(),
					 new TGeoRotation("",ph,th,0.0)),i-1);
        s  = d.Mag();
        sprintf(nam,"%sPart%dCable",name,i-1);
        fTubes->AddAt( new TGeoCtub(nam,fRmin,fRmax,0.5*s,0.0,360.0,n0.x(),
				    n0.y(),n0.z(),n1.x(),n1.y(),n1.z()),i-1);
        n0 = -n1;
        x0 = x1;
    } // end for i
}
//----------------------------------------------------------------------
AliITSGeoCable::~AliITSGeoCable(){
    //
    // Inputs:
    //    none.
    // Outputs:
    //    none.
    // Return:
    //    none.
    Int_t i;

    if(fTubes){
        for(i=0;i<fTubes->GetEntries();i++) 
            delete (TGeoCtub*)(fTubes->At(i));
        delete fTubes;
    } // end if
    fTubes = 0;
    if(fTranRot){
        for(i=0;i<fTranRot->GetEntries();i++) 
            delete (TGeoCombiTrans*)(fTranRot->At(i));
        delete fTranRot;
    } // end if
    fTranRot  = 0;
}
//----------------------------------------------------------------------
