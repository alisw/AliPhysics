// $Id: AliPhJBaseTrack.h,v 1.5 2008/05/08 15:19:52 djkim Exp $

////////////////////////////////////////////////////
/*!
  \file AliPhJBaseTrack.h
  \brief
  \author J. Rak, D.J.Kim, R.Diaz (University of Jyvaskyla)
  \email: djkim@jyu.fi
  \version $Revision: 1.5 $
  \date $Date: 2008/05/08 15:19:52 $
 */
////////////////////////////////////////////////////

#ifndef ALIPHJBASETRACK_H
#define ALIPHJBASETRACK_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include  "JConst.h"


class AliPhJBaseTrack : public TObject {

    public:
        AliPhJBaseTrack();
        AliPhJBaseTrack(float pt,float theta, float eta, float phi, float e, short charge, float tof,float ptot); // constructor
        AliPhJBaseTrack(const AliPhJBaseTrack& a);

        virtual ~AliPhJBaseTrack(){;}		//destructor

        float  GetPt()       const {return fBasePt;}
        float  GetTheta()    const {return fBaseTheta;}
        float  GetEta()      const {return fBaseEta;}
        float  GetPhi()      const {return fBasePhi;}
        float  GetTwoPiPhi() const {return fBasePhi>-kJPi/3 ? fBasePhi : kJTwoPi+fBasePhi;}
        float  GetE()        const {return fBaseE;}
	short  GetCharge()   const {return fBaseCharge;}
        int    GetTrackID()  const {return fBaseID;}
        particleType  GetFlavor() const  {return fBaseFlavor;}
        float  GetTof()     const  {return fBaseTof;}
        float  GetPtot()    const  {return fBasePtot;}
        float  GetPx()       const {return fBasePt*cos(fBasePhi);}
        float  GetPy()       const {return fBasePt*sin(fBasePhi);}
        float  GetPz()       const {return fBasePt*sinh(fBaseEta);}

        void SetPt(float inpt) {fBasePt=inpt;}
        void SetTheta(float inth) {fBaseTheta=inth;}
        void SetEta(float inEta) {fBaseEta=inEta;}
        void SetPhi(float inphi) {fBasePhi=inphi;}
        void SetE(float inE) {fBaseE=inE;}
	void SetCharge(short q) {fBaseCharge=q;}
        void SetFlavor(particleType ptye) {fBaseFlavor=ptye;}
        void SetTrackID(int inID){fBaseID=inID;}
        void SetTof(float tof) {fBaseTof=tof;}
        void SetPtot(float ptot){ fBasePtot = ptot;}

        void PrintOut(const char *message) const;
        void PrintJetInput(char *message) const;

       AliPhJBaseTrack& operator=(const AliPhJBaseTrack& trk);

    protected:

        float fBasePt, fBaseTheta, fBaseEta, fBasePhi, fBaseE;// track pT, theta, eta, phi, E
        int fBaseID;// unique track id
        short fBaseCharge;//track charge
        particleType fBaseFlavor;//track flavour
        float  fBaseTof, fBasePtot; //track tof and totoal momentum

        ClassDef(AliPhJBaseTrack,1)
};

#endif

