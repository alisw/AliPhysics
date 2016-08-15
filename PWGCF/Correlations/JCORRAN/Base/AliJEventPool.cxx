#include <TRandom.h>
#include  "AliJEventPool.h"

#include <TH1D.h>

#include  "AliJTrack.h"
#include  "AliJPhoton.h"
#include <TClonesArray.h>

#include  "AliJBaseTrack.h"
#include  "AliJPhoton.h"


#include  "AliJPiZero.h"

#include  "AliJCard.h"
#include  "AliJCorrelationInterface.h"
#include  "AliJHistogramInterface.h"


// event mixing
// blah
// blah
// blah
// blah

AliJEventPool::AliJEventPool(AliJCard *cardin, AliJHistogramInterface *histosin, AliJCorrelationInterface *coin, particleType particle ) :
  fcard(cardin),
  fcorrelations(coin),
  fhistos(histosin),
  //ftk(NULL),
  //ftk1(NULL),
  //ftk2(NULL),
  fthisPoolType(particle),
  fpoolList(NULL)
{       
  // constructor
  
  if(fcard->GetNoOfBins(kCentrType) > kMaxNoCentrBin ){
    cout<<"ERROR: No of Centrality bins exceed max dim in AliJEventPool.h"<<endl;
    exit(0);
  }
  for(int ic=0;ic<fcard->GetNoOfBins(kCentrType);ic++){
    cout<<"Mixing pool depth for icbin = "<< ic << " is " << fcard->GetEventPoolDepth(ic) <<" for <"<<kParticleTypeStrName[particle] <<"> and prototype "<<kParticleProtoType[particle]<<endl;
    if(fcard->GetEventPoolDepth(ic) > MAXNOEVENT ){
      cout<<"ERROR: Event pool depth exeed max="<<MAXNOEVENT<<" in AliJEventPool.h"<<endl;
      exit(0);
    }
  } cout <<endl; 

  for(int ic=0;ic<fcard->GetNoOfBins(kCentrType);ic++){
    for(int ie=0;ie<fcard->GetEventPoolDepth(ic); ie++){ 
      fLists[ic][ie]  = new TClonesArray(kParticleProtoType[particle],1500);
    }
    flastAccepted[ic] = -1; //to start from 0
    fwhereToStore[ic] = -1; //to start from 0
    fnoMix[ic]    = 0;
    fnoMixCut[ic] = 0;
  }

  //ftk  = new AliJBaseTrack;
  //ftk1 = new AliJBaseTrack;
  //ftk2 = new AliJBaseTrack;

}

AliJEventPool::~AliJEventPool( ){
  // destructor
  //delete ftk;
  //delete ftk1;
  //delete ftk2;
}  

AliJEventPool::AliJEventPool(const AliJEventPool& obj) :
  fcard(obj.fcard),
  fcorrelations(obj.fcorrelations),
  fhistos(obj.fhistos),
  //ftk(obj.ftk),
  //ftk1(obj.ftk1),
  //ftk2(obj.ftk2),
  fthisPoolType(obj.fthisPoolType),
  fpoolList(obj.fpoolList)
{
  // copy constructor
  JUNUSED(obj);
}

AliJEventPool& AliJEventPool::operator=(const AliJEventPool& obj){
  // equal sign operator
  JUNUSED(obj);
  return *this;
}
  

//______________________________________________________________________________
void AliJEventPool::Mix( TClonesArray *triggList, 
        corrFillType cFTyp, 
        float cent, float Z, float thisMult, int iev, bool leadingParticle){
  // mixer
    int cBin = fcard->GetBin(kCentrType, cent);
    int zBin = fcard->GetBin(kZVertType, Z);
    int noTrigg=triggList->GetEntriesFast();
    int noAssoc=0;

    if ( cBin< 0 ) return;
   
//     cout << "c: " << cBin << endl;


    for(int backCounter=0; backCounter <= flastAccepted[cBin]; backCounter++){
        fpoolList = fLists [cBin] [backCounter];
        noAssoc = fpoolList->GetEntries();

        if(noAssoc<=0) continue;

        //mixit=======
        fnoMix[cBin]++;
        
        if( 
                fcard->SimilarCentrality(fcentrality[cBin][backCounter], cent, cBin) &&
                fcard->SimilarMultiplicity(fmult[cBin][backCounter], thisMult) &&
                fcard->GetBin(kZVertType, fZVertex[cBin][backCounter])==zBin      &&
                fevent[cBin][backCounter] != iev )
        {
            fnoMixCut[cBin]++;
            //=================================================
            // try to use only one track from each fevent
            //=================================================
            for(int ii=0;ii<noTrigg;ii++){
                AliJBaseTrack *ftk1 = (AliJBaseTrack*)triggList->At(ii);        
                //fhistos->fhTriggPtBin[kMixed][cBin][iptt]->Fill(ptt); //who needs that?
                for(int jj=0;jj<noAssoc ;jj++){
                    AliJBaseTrack *ftk2 = (AliJBaseTrack*)fpoolList->At(jj);
                    if(leadingParticle && ftk1->Pt() < ftk2->Pt()) continue; // In leading particle correlations, accept only those associated particles whose pT is lower than that of the trigger
                    fcorrelations->FillHisto(cFTyp,kMixed, cBin, zBin, ftk1, ftk2);
                } //inner loop mixing
            }//outer loop mixing
        }//if good for mix
    }//mixed fevent loop
}

//______________________________________________________________________________
void AliJEventPool::AcceptList(TClonesArray *inList, float cent, float Z, float inMult, int iev){
    //////////////////////////////////////////////////////////////
    // circular buffer. New fevent added after the previous one, 
    // mixing goes backwards 
    //////////////////////////////////////////////////////////////
    int cBin = fcard->GetBin(kCentrType, cent);
    if (cBin <0 ) return;
    flastAccepted[cBin]++;
    fwhereToStore[cBin]++;
    if( flastAccepted[cBin] >= fcard->GetEventPoolDepth(cBin) ) flastAccepted[cBin] = fcard->GetEventPoolDepth(cBin)-1;
    if( fwhereToStore[cBin] >= fcard->GetEventPoolDepth(cBin) ) fwhereToStore[cBin] = 0;
    fevent     [cBin][fwhereToStore[cBin]] = iev;
    fZVertex   [cBin][fwhereToStore[cBin]] = Z;
    fcentrality[cBin][fwhereToStore[cBin]] = cent;
    fmult      [cBin][fwhereToStore[cBin]] = inMult;

    fLists[cBin][fwhereToStore[cBin]]->Clear();
    for(int i=0;i<inList->GetEntriesFast();i++){
				if( fthisPoolType == kJPhoton || fthisPoolType == kJDecayphoton ){
					AliJPhoton *tkp = (AliJPhoton*)inList->At(i);
					new ((*fLists[cBin][fwhereToStore[cBin]])[i]) AliJPhoton(*tkp);
				}
				else if( fthisPoolType == kJPizero || fthisPoolType == kJEta ){
					AliJPiZero *tkpz = (AliJPiZero*)inList->At(i);
					new ((*fLists[cBin][fwhereToStore[cBin]])[i]) AliJPiZero(*tkpz);
				}
				else{
					AliJTrack *tk3 = (AliJTrack*)inList->At(i);
					new ((*fLists[cBin][fwhereToStore[cBin]])[i]) AliJTrack(*tk3);
				}
    }

}





//==================== Sampling ===========================
void AliJEventPool::Mysample(TH1D *fromh, TH1D *toh )
{
  // sampler
    if(fabs(toh->GetBinWidth(1)-fromh->GetBinWidth(1))>1e-5){
        cout<<" Attempt to integrate histograms with non-compatible binning"<<endl;
        printf("%16.14f %16.14f %16.14f \n", fromh->GetBinWidth(1), toh->GetBinWidth(1), fabs(toh->GetBinWidth(1)-fromh->GetBinWidth(1)) );
        return;
    }

    for(int idphi=1; idphi<=toh->GetNbinsX(); idphi++){
        double sum=0,err=0, a,b,e;
        for (int i = 1;i <= fromh->GetNbinsX();i++)
        {
            double phi1 = fromh->GetBinCenter(i);
            double fdphi = toh->GetBinCenter(idphi);
            double phi2 = atan2(sin(phi1+fdphi),cos(phi1+fdphi));
            a=fromh->GetBinContent(i);
            b=fromh->GetBinContent(fromh->FindBin(phi2));
            sum += a*b;
            err += a*b*(a+b);
        }
        if( idphi==1 || idphi== toh->GetNbinsX() ) sum /= 2.0;
        toh->SetBinContent(idphi,toh->GetBinContent(idphi)+sum);
        e=toh->GetBinError(idphi);
        toh->SetBinError(idphi,sqrt(err+e*e));
    }
}

//______________________________________________________________________________





















