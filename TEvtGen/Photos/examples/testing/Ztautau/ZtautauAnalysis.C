//#include "UserTreeAnalysis.H"   // remove if copied to user working directory
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "MC4Vector.H"
#include "HEPParticle.H"
#include "TH1.h"
#include "Setup.H"
#include "TObjArray.h"
#include "TMath.h"
#define PI 3.141592653

inline bool ifSofty(int Id,int nparams, double *params){
  if (nparams<5 && Id==22) return true;   // to remove photons only
  for (int i=nparams-1; i>3; i--)
    if (Id==params[i]) return true;       // to remove all what is in params from nparams down to 4
  return false;
}

// very similar to  MC_FillUserHistogram from Generate.cxx
inline void fillUserHisto(char *name,double val, double weight=1.0, 
                          double min=Setup::bin_min[0][0], 
                          double max=Setup::bin_max[0][0]){

    TH1D *h=(TH1D*)(Setup::user_histograms->FindObject(name));
    if(!h){
      h=new TH1D(name,name,Setup::nbins[0][0],min,max);
      if(!h) return;
      Setup::user_histograms->Add(h);
      //      printf("user histogram created %s\n", name);
}
    h->Fill(val,weight);

}
double angle(double X, double Y){

 
    double  an=0.0; 
    double R=sqrt(X*X+Y*Y); 
    //    if(R<pow(10,-20)) printf(" angle this time X %f\n", 10000*X);
    if(R<pow(10,-20)) return an;

    if(TMath::Abs(X)/R<0.8) 
    { 
        an=acos(X/R);
	if(Y<0 && an>0) an=-an;
        if(Y>0 && an<0) an=-an; 
    }
    else 
    { 
      an=asin(Y/R);
      if(X<0 && an>=0.)  an=PI-an; 
      else if(X<0.)      an=-PI-an; 
       
    } 
  return an;
}

int counter=0;

int ZtautauAnalysis(HEPParticle *mother,HEPParticleList *stableDaughters, int nparams, double *params)
{
  // THIS IS EXAMPLE of userTreeAnalysis. It acts on list of particle X (under our MC-test) final
  // daughters. Of course in real life many options may need to be introduced by the user.
    // PARAMETERS:
    // params[0] threshold on Energy (or p_T) of particle expressed as a fraction of mother's 
    //		 Energy (or p_T) in mother's frame. If not specified - default is 0.05 
    //           Note that setting it to zero is possible.
    // params[1] maximum number of left soft-suspected particles. 
    //           0: (default) all listed particles are removed, even if hard
    //           1: 2:  one two etc  removable particles are kept (at most)
    //          -1: this option is off, all hard particles stay.
    // params[2] control tag on discriminating particles, 
    //           0: (default) Energy in decaying particle rest frame
    //           1: Energy in lab frame
    //           2: p_T with respect to beam direction in lab frame.
    // params[3] type of action to be applied on removed daughters 
    //           0: (default) nothing is done, they are lost
    //           1:  algorithm as in studies on PHOTOS is used, see papers P. Golonka and Z. Was.
    // params[4] from this paramter on PDG id-s of particles to be removed can be added. 
    //           Default is that only photons are removed.
    // To get to origin of our particle (mother) one need to go after UserEventAnalysis,
    // instructive example will be given later. 
    assert(mother!=0);
    assert(stableDaughters!=0);


    double threshold=0.05, leftmax=0.0, selector=0.0, actLost=0.0; 
    if (nparams>0 && params==0)  return -1;
    if (nparams>0) threshold=params[0]; 
    if (nparams>1) leftmax  =params[1]; 
    if (nparams>2) selector =params[2]; 
    if (nparams>3) actLost  =params[3]; 

    
    HEPParticleList *lostDaughters=new HEPParticleList;    
    double pt_limit=threshold*mother->GetM();

    HEPParticleListIterator daughters (*stableDaughters);
    int nphot=0;
    double ephot=pow(10,22);
    bool redo=false;
    for (HEPParticle *part=daughters.first(); part!=0; (redo)? part=part:part=daughters.next() ) {
      if(redo) redo=false;
	MC4Vector d4(part->GetE(),part->GetPx(),part->GetPy(),part->GetPz(),part->GetM());
	// boost to mother's frame:
	d4.Boost(mother->GetPx(),mother->GetPy(),mother->GetPz(),mother->GetE(),mother->GetM());



	double p_t=d4.GetX0();  // default
	switch ((int) selector)
        {
	case 1: p_t=part->GetE(); break;
	case 2: p_t=d4.Xt(); break;
	}


        if(ifSofty(part->GetPDGId(),nparams,params)) nphot++;
	if ( ifSofty(part->GetPDGId(),nparams,params) && p_t < pt_limit) {
	  //	  printf("Odrzucamy! %i\n",count++);
          nphot=0;
            lostDaughters->push_back(part);
	    stableDaughters->remove(part);
	    part=daughters.first();
            redo=true;
            continue;
	}
        if( ifSofty(part->GetPDGId(),nparams,params) && ephot>p_t) ephot=p_t;
    }
    while(nphot>(int) leftmax)
    {
      double ephot1=pow(10,22);
      redo=false;
      for (HEPParticle *part=daughters.first(); part!=0; (redo)? part=part:part=daughters.next() ) {
        if(redo) redo=false;
	MC4Vector d4(part->GetE(),part->GetPx(),part->GetPy(),part->GetPz(),part->GetM());
	// boost to mother's frame:
	d4.Boost(mother->GetPx(),mother->GetPy(),mother->GetPz(),mother->GetE(),mother->GetM());


	double p_t=d4.GetX0();  // default
	switch ((int) selector)
        {
	case 1: p_t=part->GetE(); break;
	case 2: p_t=d4.Xt(); break;
	}


	if ( ifSofty(part->GetPDGId(),nparams,params) && p_t == ephot) {

	  nphot--;
            lostDaughters->push_back(part);
	    stableDaughters->remove(part);
            ephot1=pow(10,22);
	     part=daughters.first();
             redo=true;
             continue;
	}
        if( ifSofty(part->GetPDGId(),nparams,params) && ephot1>p_t) ephot1=p_t;
      }
      ephot=ephot1;

    }

    // ##############################################################
    // here functionality of removig some daughters is finished
    // now we reconsider what to do with them
    // ##############################################################
    // delete lostDaughters;
    // return 1; 


    // ##############################################################
    // Now: What to do with lost daughters?
    // ##############################################################
    //    lostDaughters->ls();
    int version=(int) actLost;
        
    switch(version)
      {
      case 1:  // add lost to charged
	{
          HEPParticleListIterator lost (*lostDaughters);
	  for (HEPParticle *lostpart=lost.first(); lostpart!=0; lostpart=lost.next() ) {
            HEPParticle* Best=0;
	    double searchvirt=pow(10.0,22);
            MC4Vector VV;
	    for( HEPParticle *part=daughters.first(); part!=0; part=daughters.next() ){
	      if(part->GetCharge()==0) continue;
              VV=lostpart->GetP4()+part->GetP4();
	      VV.AdjustM();
	      if (VV.GetM()<searchvirt) {searchvirt=VV.GetM(); Best=part;}
	    }
	    if(Best) {
              Best->SetPx(Best->GetPx()+lostpart->GetPx());
              Best->SetPy(Best->GetPy()+lostpart->GetPy());
              Best->SetPz(Best->GetPz()+lostpart->GetPz());
              Best->SetE (Best->GetE ()+lostpart->GetE ());
	    }
	  }
	break;
	}
      default: break; // do nothing
      }

    delete lostDaughters;


	/*
		Create histogram
	*/
	double px=0,py=0,pz=0,e=0;
	double px1=0,py1=0,pz1=0,e1=0;

	//loop over all daughters
	for(HEPParticle *part=daughters.first();part!=0;part=daughters.next())
	{
		if(abs(part->GetPDGId())==211)
		{
			px1+=part->GetPx();
			py1+=part->GetPy();
			pz1+=part->GetPz();
			e1+=part->GetE();
		}
		if(abs(part->GetPDGId())==211) continue;
		//Sum e+ e- 4-vectors
		px+=part->GetPx();
		py+=part->GetPy();
		pz+=part->GetPz();
		e+=part->GetE();
	}
	double mass=e*e-px*px-py*py-pz*pz;
	double mass1=(e+e1)*(e+e1)-(px+px1)*(px+px1)-(py+py1)*(py+py1)-(pz+pz1)*(pz+pz1);
	mass=1-mass/mass1;
	char mm[] = "pi-energy";
	fillUserHisto(mm,mass,1.0,0.0,1.1);

    //For ZtautauAnalysis - other histograms are not needed
    bool activateUserHist=false;
    if(!activateUserHist) return 1;

    // segmet of user defined histograms

    double X=mother->GetPx(), Y=mother->GetPy(), Z=mother->GetPz(); 
    // double E=mother->GetE(),  MM=mother->GetM();


    double pt=sqrt(X*X+Y*Y);
    double eta=log((sqrt(pt*pt+Z*Z)+TMath::Abs(Z))/pt);
    if(Z<0 && eta>0) eta=-eta;
    if(Z>0 && eta<0) eta=-eta;
    double phi=angle(X,Y);
    char hist1[] = "mother-PT";
    char hist2[] = "mother-eta";
    char hist3[] = "mother-phi";
    fillUserHisto(hist1,pt,1.0,0.0,100.0);
    fillUserHisto(hist2,eta,1.0,-8.0,8.0);
    fillUserHisto(hist3,phi,1.0,-PI,PI);
    return 1;
}
