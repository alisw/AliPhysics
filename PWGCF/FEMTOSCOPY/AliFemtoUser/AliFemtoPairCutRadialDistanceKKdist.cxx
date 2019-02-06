/////////////////////////////////////////////////////////////////////////////////
//                                                                             //
// AliFemtoPairCutRadialDistanceKKdist - a pair cut which checks                     //
// for some pair qualities that attempt to identify slit/doubly                //
// reconstructed tracks and also selects pairs based on their separation       //
// at the entrance to the TPC                                                  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////
/********************************************************************************
 *
 * Author: Johanna Gramling, University of Heidelberg, jgramlin@cern.ch
 *         Malgorzata Janik, Warsaw University of Technology, majanik@cern.ch
 *         Lukasz Graczykowski, Warsaw University of Technology, lgraczyk@cern.ch
 *
 ********************************************************************************/

#include "AliFemtoPairCutRadialDistanceKKdist.h"
#include <string>
#include <cstdio>

#ifdef __ROOT__
ClassImp(AliFemtoPairCutRadialDistanceKKdist)
#endif

//__________________
AliFemtoPairCutRadialDistanceKKdist::AliFemtoPairCutRadialDistanceKKdist():
AliFemtoPairCutAntiGamma(),
  //fDPhiStarMin(0),
  //fEtaMin(0),
  fMinRad(0.0)
 // fMaxRad(2.5),
 // fMagSign(1),
  //fPhistarmin(kTRUE)
{
}
//__________________
AliFemtoPairCutRadialDistanceKKdist::AliFemtoPairCutRadialDistanceKKdist(const AliFemtoPairCutRadialDistanceKKdist& c) :
  AliFemtoPairCutAntiGamma(c),
 // fDPhiStarMin(0),
//  fEtaMin(0),
  fMinRad(0.0)
  //fMaxRad(2.5),
  //fMagSign(1),
  //fPhistarmin(kTRUE)
{
 // fDPhiStarMin = c.fDPhiStarMin;
 // fEtaMin = c.fEtaMin;
  fMinRad = c.fMinRad;
  //fMaxRad = c.fMaxRad;
 // fMagSign = c.fMagSign;
 // fPhistarmin = c.fPhistarmin;
}

//__________________
AliFemtoPairCutRadialDistanceKKdist::~AliFemtoPairCutRadialDistanceKKdist(){
  /* no-op */
}
AliFemtoPairCutRadialDistanceKKdist& AliFemtoPairCutRadialDistanceKKdist::operator=(const AliFemtoPairCutRadialDistanceKKdist& c)
{
  if (this != &c) {
//    fDPhiStarMin = c.fDPhiStarMin;
//    fEtaMin = c.fEtaMin;
    fMinRad = c.fMinRad;
//    fMaxRad = c.fMaxRad;
//    fMagSign = c.fMagSign;
//    fPhistarmin = c.fPhistarmin;

  }

  return *this;
}
//__________________
bool AliFemtoPairCutRadialDistanceKKdist::Pass(const AliFemtoPair* pair){
  // Accept pairs based on their TPC AVERAGE separation and
  // quality and sharity
   // The radii at which we get the global positions
      // Float_t Rwanted[9] = {85., 105., 125., 145., 165., 185., 205., 225., 245.};

   double avgSep = 0;
   int count = 0;

// loop through the 8 points of the 'NominalTpcPoint' methods
  for (int i = 0; i < 8; i++) {
     const AliFemtoThreeVector &point_1 = pair->Track1()->Track()->NominalTpcPoint(i),
                               &point_2 = pair->Track2()->Track()->NominalTpcPoint(i);

            if (point_1.x() < -9000 || point_1.y() < -9000 || point_1.z() < -9000 ||
                point_2.x() < -9000 || point_2.y() < -9000 || point_2.z() < -9000  )
                 {
                   break;
                  }

        avgSep += (point_1 - point_2).Mag();
                  count++;
                          }
        if (count != 0) {
                        //    output->Fill(avgSep / count);
        avgSep = avgSep/count;
                         }

        Bool_t pass5 = kTRUE;

         if(avgSep < fMinRad)  pass5 = kFALSE;

         //cout<<" avgSep "<<avgSep<<" fMinRad " <<fMinRad<<" pass5 "<<pass5<<endl;

 
  if (pass5) {
    pass5 = AliFemtoPairCutAntiGamma::Pass(pair);
  }
  else {
    fNPairsFailed++;
  }


  return pass5;
}
//__________________
AliFemtoString AliFemtoPairCutRadialDistanceKKdist::Report(){
  // Prepare a report from the execution
  string stemp = "AliFemtoRadialDistance Pair Cut - remove shared and split pairs and pairs with small average separation at the specified radius\n";  char ctemp[100];
  snprintf(ctemp , 100, "Accept pair with average separation more that %f",fMinRad);
  stemp += ctemp;
  snprintf(ctemp , 100, "Number of pairs which passed:\t%ld  Number which failed:\t%ld\n",fNPairsPassed,fNPairsFailed);
  stemp += ctemp;
  AliFemtoString returnThis = stemp;
  return returnThis;}
//__________________

TList *AliFemtoPairCutRadialDistanceKKdist::ListSettings()
{
  // return a list of settings in a writable form
  TList *tListSetttings =  AliFemtoPairCutAntiGamma::ListSettings();
  char buf[200];
  snprintf(buf, 200, "AliFemtoPairCutRadialDistanceKKdist.average.separation=%f", fMinRad);
  tListSetttings->AddLast(new TObjString(buf));

  return tListSetttings;
}


void AliFemtoPairCutRadialDistanceKKdist::SetAverageSeparation(double minrad)
{
  fMinRad = minrad;
}

