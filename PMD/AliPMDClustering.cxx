/***************************************************************************
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

//-----------------------------------------------------//
//                                                     //
//  Source File : PMDClustering.cxx, Version 00        //
//                                                     //
//  Date   : September 26 2002                         //
//                                                     //
//  clustering code for alice pmd                      //
//                                                     //
//-----------------------------------------------------//

/* 
   --------------------------------------------------------------------
   Code developed by S. C. Phatak, Institute of Physics, 
   Bhubaneswar 751 005 ( phatak@iopb.res.in ) Given the energy deposited
   ( or ADC value ) in each cell of supermodule ( pmd or cpv ), the code
   builds up superclusters and breaks them into clusters. The input is 
   in array d[ndimx][ndimy] and cluster information is in array
   clusters[5][5000]. integer clno gives total number of clusters in the 
   supermodule.

   d, clno  and clusters are the only global ( public ) variables. Others 
   are local ( private ) to the code. 

   At the moment, the data is read for whole detector ( all supermodules
   and pmd as well as cpv. This will have to be modify later )

   LAST UPDATE  :  October 23, 2002
-----------------------------------------------------------------------
*/



#include <TNtuple.h>
#include <TObjArray.h>
#include "AliPMDcluster.h"
#include "AliPMDClustering.h"
#include <stdio.h>

ClassImp(AliPMDClustering)

const double AliPMDClustering::pi=3.141593;
const double AliPMDClustering::sqrth=0.8660254;  // sqrth = sqrt(3.)/2.


AliPMDClustering::AliPMDClustering()
{
  fDebug  = 0;
  fCutoff = 0;
  for(int i = 0; i < ndimx; i++)
    {
      for(int j = 0; j < ndimy; j++)
	{
	  coord[0][i][j] = i+j/2.;
	  coord[1][i][j] = sqrth*j;
	}
    }
}
AliPMDClustering::~AliPMDClustering()
{

}

void AliPMDClustering::DoClust(double celladc[48][96], TObjArray *pmdcont)
{

  AliPMDcluster *pmdcl = 0;

  int i, i1, i2, j, nmx1, incr;
  double  cutoff, ave;
  Float_t clusdata[5];

  const float twobysqrt3 = 1.1547; // 2./sqrt(3.)


  for (i = 0; i < ndimx; i++)
    {
      for (j = 0; j < ndimy; j++)
	{
	  d[i][j] = celladc[i][j];
	}
    }
  order(); // order the data
  //  cutoff=400.; // cutoff used to discard cells having ener. dep. 
  cutoff = fCutoff; // cutoff used to discard cells having ener. dep. 
  ave=0.; 
  nmx1=-1;
  
  for(j=0;j<nmx; j++)
    {
      i1 = iord[0][j];
      i2 = iord[1][j];
      if (d[i1][i2] > 0.) {ave=ave+d[i1][i2];}
      if (d[i1][i2] >= cutoff ) nmx1 = nmx1 + 1;
    }
  // nmx1 --- number of cells having ener dep >= cutoff
  if (fDebug == 1)
    {
      cout << " nmx1 " << nmx1 << endl;
    }
  ave=ave/nmx1;
  if (fDebug == 1)
    {
      cout <<"nmx " << nmx << " nmx1 " << nmx1<< " ave "<<ave<<
	" cutoff " << cutoff << endl;
    }
  
  incr = crclust(ave, cutoff, nmx1);
  
  refclust(incr);
  
  if (fDebug == 1)
    {
      cout << "clno " << clno << endl;
    }
  
  for(i1=0; i1<clno; i1++)
    {
      float clu_xc    = (float) clusters[0][i1];
      float clu_yc    = (float) clusters[1][i1];
      float clu_adc   = (float) clusters[2][i1];
      float clu_cells = (float) clusters[3][i1];
      float clu_rad   = (float) clusters[4][i1];
      
      float clu_y0 = twobysqrt3*clu_yc;
      float clu_x0 = clu_xc - clu_y0/2.;

      clusdata[0] = clu_cells;
      clusdata[1] = clu_x0;
      clusdata[2] = clu_y0;
      clusdata[3] = clu_adc;
      clusdata[4] = clu_rad;
      
      pmdcl = new AliPMDcluster(clusdata);
      pmdcont->Add(pmdcl);
    }
  delete pmdcl;
}

void AliPMDClustering::order()
{
  // using simple sort
  double dd[nmx], adum;// matrix d converted into 
  // one dimensional array dd. adum a place holder for double
  int i, j, i1, i2, iord1[nmx], itst, idum; // information of 
  // ordering is stored in iord1, original array not ordered
  //
  // define arrays dd and iord1
  for(i1=0; i1 < ndimx; i1++){
    for(i2=0; i2 < ndimy; i2++){
      i=i1+i2*ndimx;
      iord1[i]=i; dd[i]=d[i1][i2];
    }
  }
  // sort and store sorting information in iord1 
  for(j=1; j < nmx; j++){
    itst=0; adum=dd[j]; idum=iord1[j];
    for(i1=0; i1 < j ; i1++){
      if(adum > dd[i1] && itst == 0){
	itst=1;
	for(i2=j-1; i2 >= i1 ; i2=i2--){
	  dd[i2+1]=dd[i2]; 
	  iord1[i2+1]=iord1[i2];
	}
	dd[i1]=adum; iord1[i1]=idum;
      }
    }
  }
  // store the sorted information in iord for later use
  for(i=0; i<nmx; i++){
    j  = iord1[i];
    i2 = j/ndimx; 
    i1 = j-i2*ndimx; 
    iord[0][i]=i1; 
    iord[1][i]=i2;
  }
}


  
int AliPMDClustering::crclust(double ave, double cutoff, int nmx1)
{
  int i,j,k,id1,id2,icl, numcell, clust[2][5000];
  int jd1,jd2, icell, cellcount;
  static int neibx[6]={1,0,-1,-1,0,1}, neiby[6]={0,1,1,0,-1,-1};
  // neibx and neiby define ( incremental ) (i,j) for the neighbours of a 
  // cell. There are six neighbours.
  // cellcount --- total number of cells having nonzero ener dep
  // numcell --- number of cells in a given supercluster
  //ofstream ofl0("cells_loc",ios::out);
  // initialize infocl[2][ndimx][ndimy] 

  if (fDebug == 1)
    {
      printf(" *** Inside crclust **  nmx = %d nmx1 = %d ndimx = %d ndimy = %d ave = %f cutoff = %f\n",
	     nmx,nmx1,ndimx,ndimy,ave,cutoff);
    }
  for (j=0; j < ndimx; j++){
    for(k=0; k < ndimy; k++){
      infocl[0][j][k] = 0; 
      infocl[1][j][k] = 0;
    }
  }
  for(i=0; i < nmx; i++){
    infcl[0][i] = -1;
    id1=iord[0][i]; 
    id2=iord[1][i];
    if(d[id1][id2] <= cutoff){infocl[0][id1][id2]=-1;}
  }
  // ---------------------------------------------------------------
  // crude clustering begins. Start with cell having largest adc 
  // count and loop over the cells in descending order of adc count
  // ---------------------------------------------------------------
  icl=-1;
  cellcount=-1;
  for(icell=0; icell <= nmx1; icell++){
    id1=iord[0][icell]; 
    id2=iord[1][icell]; 
    if(infocl[0][id1][id2] == 0 ){
      // ---------------------------------------------------------------
      // icl -- cluster #, numcell -- # of cells in it, clust -- stores 
      // coordinates of the cells in a cluster, infocl[0][i1][i2] is 1 for 
      // primary and 2 for secondary cells, 
      // infocl[1][i1][i2] stores cluster #
      // ---------------------------------------------------------------
      icl=icl+1; 
      numcell=0; 
      cellcount=cellcount+1;
      infocl[0][id1][id2]=1; 
      infocl[1][id1][id2]=icl;
      infcl[0][cellcount]=icl; 
      infcl[1][cellcount]=id1; 
      infcl[2][cellcount]=id2;


      clust[0][numcell]=id1;
      clust[1][numcell]=id2;
      for(i=1; i<5000; i++)clust[0][i]=0;
      // ---------------------------------------------------------------
      // check for adc count in neib. cells. If ne 0 put it in this clust
      // ---------------------------------------------------------------
      for(i=0; i<6; i++){
	jd1=id1+neibx[i]; 
	jd2=id2+neiby[i];
	//if( (jd1 >= 0 && jd1 < 72) && (jd2 >= 0 && jd2 < 72) && 
	if( (jd1 >= 0 && jd1 < ndimx) && (jd2 >= 0 && jd2 < ndimy) && 
	    infocl[0][jd1][jd2] == 0){
	  numcell=numcell+1;
	  infocl[0][jd1][jd2]=2; 
	  infocl[1][jd1][jd2]=icl;
	  clust[0][numcell]=jd1;
	  clust[1][numcell]=jd2;
	  cellcount=cellcount+1;
	  infcl[0][cellcount]=icl; 
	  infcl[1][cellcount]=jd1; 
	  infcl[2][cellcount]=jd2;
	}
      }
      // ---------------------------------------------------------------
      // check adc count for neighbour's neighbours recursively and 
      // if nonzero, add these to the cluster.
      // ---------------------------------------------------------------
      for(i=1;i < 5000;i++){
	if(clust[0][i] != 0){
	  id1=clust[0][i]; 
	  id2=clust[1][i];
	  for(j=0; j<6 ; j++){
	    jd1=id1+neibx[j]; 
	    jd2=id2+neiby[j];
	    //if( (jd1 >= 0 && jd1 < 72) && (jd2 >= 0 && jd2 < 72) && 
	    if( (jd1 >= 0 && jd1 < ndimx) && (jd2 >= 0 && jd2 < ndimy) && 
		infocl[0][jd1][jd2] == 0 ){
	      infocl[0][jd1][jd2]=2; 
	      infocl[1][jd1][jd2]=icl;
	      numcell=numcell+1; 
	      clust[0][numcell]=jd1;
	      clust[1][numcell]=jd2;
	      cellcount=cellcount+1;
	      infcl[0][cellcount]=icl; 
	      infcl[1][cellcount]=jd1; 
	      infcl[2][cellcount]=jd2;
	    }
	  }
	}
      }
    }
  }
  //  for(icell=0; icell<=cellcount; icell++){
  //    ofl0 << infcl[0][icell] << " " << infcl[1][icell] << " " << 
  //      infcl[2][icell] << endl;
  //  }
  return cellcount;
}

void AliPMDClustering::refclust(int incr)
{
  int i, j, k, i1, i2, id, icl, ncl[4500], iord[4500], itest; 
  int ihld;
  int ig, nsupcl, lev1[20], lev2[20];
  double x[4500], y[4500], z[4500], x1, y1, z1, x2, y2, z2, dist;
  double xc[4500], yc[4500], zc[4500], cells[4500], sum, rc[4500], rr;
  // clno counts the final clusters
  // nsupcl =  # of superclusters; ncl[i]= # of cells in supercluster i
  // x, y and z store (x,y) coordinates of and energy deposited in a cell
  // xc, yc store (x,y) coordinates of the cluster center
  // zc stores the energy deposited in a cluster
  // rc is cluster radius
  // finally the cluster information is put in 2-dimensional array clusters
  //  ofstream ofl1("checking.5",ios::app);
  clno=-1;
  nsupcl=-1;
  for(i=0; i<4500; i++){ncl[i]=-1;}
  for(i=0; i<incr; i++){
    if(infcl[0][i] != nsupcl){ nsupcl=nsupcl+1; }
    ncl[nsupcl]=ncl[nsupcl]+1;
  }
  if (fDebug == 1)
    {
      cout << " # of cells " <<incr+1 << " # of superclusters " << nsupcl+1
	   << endl;
    }
  id=-1;
  icl=-1;
  for(i=0; i<nsupcl; i++){
    if(ncl[i] == 0){ 
      id=id+1; 
      icl=icl+1;
      // one  cell super-clusters --> single cluster
      // cluster center at the centyer of the cell
      // cluster radius = half cell dimension
      clno=clno+1; 
      i1=infcl[1][id]; 
      i2=infcl[2][id];
      clusters[0][clno]=coord[0][i1][i2]; 
      clusters[1][clno]=coord[1][i1][i2];
      clusters[2][clno]=d[i1][i2]; 
      clusters[3][clno]=1.; 
      clusters[4][clno]=0.5;
      //ofl1 << icl << " " << coord[0][i1][i2] << " " << coord[1][i1][i2] << 
      //" " << d[i1][i2] << " " << clusters[3][clno] <<endl; 
    }else if(ncl[i] == 1){
      // two cell super-cluster --> single cluster
      // cluster center is at ener. dep.-weighted mean of two cells
      // cluster radius == half cell dimension
      id=id+1; 
      icl=icl+1;
      clno=clno+1; 
      i1=infcl[1][id]; 
      i2=infcl[2][id]; 
      x1=coord[0][i1][i2];
      y1=coord[1][i1][i2]; 
      z1=d[i1][i2];
      id=id+1; 
      i1=infcl[1][id]; 
      i2=infcl[2][id];
      x2=coord[0][i1][i2]; 
      y2=coord[1][i1][i2]; 
      z2=d[i1][i2];
      clusters[0][clno]=(x1*z1+x2*z2)/(z1+z2); 
      clusters[1][clno]=(y1*z1+y2*z2)/(z1+z2);
      clusters[2][clno]=z1+z2; 
      clusters[3][clno]=2.; 
      clusters[4][clno]=0.5;
      //ofl1 << icl << " " << clusters[0][clno] << " " << clusters[1][clno]
      //   << " " << clusters[2][clno] << " " <<clusters[3][clno] <<endl; 
    }else{

      id=id+1; 
      iord[0]=0;
      // super-cluster of more than two cells - broken up into smaller 
      // clusters gaussian centers computed. (peaks separated by > 1 cell) 
      // Begin from cell having largest energy deposited This is first
      // cluster center
      i1=infcl[1][id]; 
      i2=infcl[2][id];
      x[0]=coord[0][i1][i2]; 
      y[0]=coord[1][i1][i2]; 
      z[0]=d[i1][i2];
      iord[0]=0;
      for(j=1;j<=ncl[i];j++){

	id=id+1;
	i1=infcl[1][id]; 
	i2=infcl[2][id];
	iord[j]=j;
	x[j]=coord[0][i1][i2]; 
	y[j]=coord[1][i1][i2]; 
	z[j]=d[i1][i2];
      }
      // arranging cells within supercluster in decreasing order 
      for(j=1;j<=ncl[i];j++){
	itest=0; ihld=iord[j];
	for(i1=0;i1<j;i1++){
	  if(itest == 0 && z[iord[i1]] < z[ihld]){
	    itest=1;
	    for(i2=j-1;i2>=i1;i2--){
	      iord[i2+1]=iord[i2];
	    }
	    iord[i1]=ihld;
	  }
	}
      }


      // compute the number of Gaussians and their centers ( first 
      // guess ) 
      // centers must be separated by cells having smaller ener. dep.
      // neighbouring centers should be either strong or well-separated
      ig=0;
      xc[ig]=x[iord[0]]; 
      yc[ig]=y[iord[0]]; 
      zc[ig]=z[iord[0]];
      for(j=1;j<=ncl[i];j++){
	itest=-1; 
	x1=x[iord[j]]; 
	y1=y[iord[j]];
	for(k=0;k<=ig;k++){
	  x2=xc[k]; y2=yc[k]; 
	  rr=Dist(x1,y1,x2,y2);
	  if( rr >= 1.1 && rr < 1.8 && z[iord[j]] > zc[k]/4.)
	    itest=itest+1;
	  if( rr >= 1.8 && rr < 2.1 && z[iord[j]] > zc[k]/10.)
	    itest=itest+1;
	  if( rr >= 2.1)itest=itest+1;
	} 
	if(itest == ig){
	  ig=ig+1; 
	  xc[ig]=x1; 
	  yc[ig]=y1; 
	  zc[ig]=z[iord[j]];
	}
      }
      // for(j=0; j<=ig; j++){
      //ofl1 << icl+j+1 << " " << xc[j] << " " <<yc[j] <<" "<<zc[j]<<endl;
      //}
      // gaussfit to adjust cluster parameters to minimize
      gaussfit(ncl[i], ig, x[0], y[0] ,z[0], xc[0], yc[0], zc[0], rc[0]);
      icl=icl+ig+1;
      // compute the number of cells belonging to each cluster.
      // cell is shared between several clusters ( if they are equidistant 
      // from it ) in the ratio of cluster energy deposition
      for(j=0; j<=ig; j++){
	cells[j]=0.;
      }
      if(ig > 0){
	for(j=0; j<=ncl[i]; j++){
	  lev1[0]=0; 
	  lev2[0]=0;
	  for(k=0; k<=ig; k++){
	    dist=Dist(x[j], y[j], xc[k], yc[k]);
	    if(dist < sqrt(3.) ){
	      lev1[0]++; 
	      i1=lev1[0]; 
	      lev1[i1]=k;
	    }else{
	      if(dist < 2.1){
		lev2[0]++; 
		i1=lev2[0]; 
		lev2[i1]=k;
	      }
	    }
	  }
	  if(lev1[0] != 0){
	    if(lev1[0] == 1){cells[lev1[1]]=cells[lev1[1]]+1.;}
	    else{
	      sum=0.;
	      for(k=1; k<=lev1[0]; k++){
		sum=sum+zc[lev1[k]];
	      }
	      for(k=1; k<=lev1[0]; k++){
		cells[lev1[k]]=cells[lev1[k]]+zc[lev1[k]]/sum;
	      }
	    }
	  }else{
	    if(lev2[0] == 0){cells[lev2[1]]=cells[lev2[1]]+1.;}
	    else{
	      sum=0.;
	      for(k=1; k<=lev2[0]; k++){
		sum=sum+zc[lev2[k]];
	      }
	      for(k=1; k<=lev2[0]; k++){
		cells[lev2[k]]=cells[lev2[k]]+zc[lev2[k]]/sum;
	      }
	    }
	  }
	}
      }
      for(j=0; j<=ig; j++){
	clno=clno+1; 
	clusters[0][clno]=xc[j]; 
	clusters[1][clno]=yc[j]; 
	clusters[2][clno]=zc[j];
	clusters[4][clno]=rc[j];
	if(ig == 0){
	  clusters[3][clno]=ncl[i];
	}else{
	  clusters[3][clno]=cells[j];
	}
      }
    }
  }

  cout << " COMING OUT of refclust" << endl;

}

void AliPMDClustering::gaussfit(int ncell, int nclust, double &x, double &y ,double &z, double &xc, double &yc, double &zc, double &rc)
{
  int i, j, i1, i2, jmax, novar, idd, jj;
  double xx[4500], yy[4500], zz[4500], xxc[4500], yyc[4500]; 
  double a[4500], b[4500], c[4500], d[4500], ha[4500], hb[4500];
  double hc[4500], hd[4500], zzc[4500], rrc[4500];
  int neib[4500][50];
  double sum, dx, dy, str, str1, aint, sum1, rr, dum;
  double x1, x2, y1, y2;
  str=0.; 
  str1=0.; 
  rr=0.3; 
  novar=0;

  j = 0;  // Just put not to see the compiler warning, BKN


  for(i=0; i<=ncell; i++){
    xx[i]=*(&x+i); 
    yy[i]=*(&y+i); 
    zz[i]=*(&z+i);
    str=str+zz[i];
  }
  for(i=0; i<=nclust; i++){
    xxc[i]=*(&xc+i); 
    yyc[i]=*(&yc+i); 
    zzc[i]=*(&zc+i); 
    str1=str1+zzc[i]; 
    rrc[i]=0.5;

  }
  for(i=0; i<=nclust; i++){
    zzc[i]=str/str1*zzc[i];
    ha[i]=xxc[i]; 
    hb[i]=yyc[i]; 
    hc[i]=zzc[i]; 
    hd[i]=rrc[i];
    x1=xxc[i]; 
    y1=yyc[i];
  }
  for(i=0; i<=ncell; i++){
    idd=0; 
    x1=xx[i]; 
    y1=yy[i];
    for(j=0; j<=nclust; j++){
      x2=xxc[j]; 
      y2=yyc[j];
      if(Dist(x1,y1,x2,y2) <= 3.){ idd=idd+1; neib[i][idd]=j; }
    }

    neib[i][0]=idd;
  }
  sum=0.;
  for(i1=0; i1<=ncell; i1++){
    aint=0.; 
    idd=neib[i1][0];
    for(i2=1; i2<=idd; i2++){
      jj=neib[i1][i2];
      dx=xx[i1]-xxc[jj]; 
      dy=yy[i1]-yyc[jj]; 
      dum=rrc[j]*rrc[jj]+rr*rr;
      aint=aint+exp(-(dx*dx+dy*dy)/dum)*zzc[idd]*rr*rr/dum;
    }
    sum=sum+(aint-zz[i1])*(aint-zz[i1])/str;
  }
  jmax=nclust*1000; 
  if(nclust > 20)jmax=20000;
  for(j=0; j<jmax; j++){
    str1=0.;
    for(i=0; i<=nclust; i++){
      a[i]=xxc[i]+0.6*(ranmar()-0.5); 
      b[i]=yyc[i]+0.6*(ranmar()-0.5);
      c[i]=zzc[i]*(1.+(ranmar()-0.5)*0.2); 
      str1=str1+zzc[i];
      d[i]=rrc[i]*(1.+(ranmar()-0.5)*0.1);
      if(d[i] < 0.25)d[i]=0.25;
    }
    for(i=0; i<=nclust; i++){ c[i]=c[i]*str/str1; }
    sum1=0.;
    for(i1=0; i1<=ncell; i1++){
      aint=0.; 
      idd=neib[i1][0];
      for(i2=1; i2<=idd; i2++){
	jj=neib[i1][i2];
	dx=xx[i1]-a[jj]; 
	dy=yy[i1]-b[jj]; 
	dum=d[jj]*d[jj]+rr*rr;
	aint=aint+exp(-(dx*dx+dy*dy)/dum)*c[i2]*rr*rr/dum;
      }
      sum1=sum1+(aint-zz[i1])*(aint-zz[i1])/str;
    }

    if(sum1 < sum){
      for(i2=0; i2<=nclust; i2++){
	xxc[i2]=a[i2]; 
	yyc[i2]=b[i2]; 
	zzc[i2]=c[i2]; 
	rrc[i2]=d[i2]; 
	sum=sum1;

      }
    }
  }
  for(j=0; j<=nclust; j++){
    *(&xc+j)=xxc[j]; 
    *(&yc+j)=yyc[j]; 
    *(&zc+j)=zzc[j]; 
    *(&rc+j)=rrc[j];
  }
}


double AliPMDClustering::Dist(double x1, double y1, double x2, double y2)
{
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double AliPMDClustering::ranmar()
{
  /*                                   C==========================C*/
  /*===================================C==========================*/
  /*  Universal random number generator proposed by Marsaglia and Zaman
      in report FSU-SCRI-87-50 */

  //  clock_t start;
  int ii, jj;
  static int i=96, j=32, itest=0, i1, i2, i3, i4, i5;
  static double u[97], c, cd, cm, s, t;
  static double uni;
  int count1,count2,idum;
  /*    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$  */
  if (itest == 0) {
    //*******************************************************
    // following three lines if the seed to be provided by computer
    // start = time(NULL);
    // ii=start;
    // jj=start;
    //*******************************************************
    //following two lines for fixed seed ( during testing only. Else
    //use preceeing three lines
    ii=8263;
    jj=5726;
    if(ii > 31328 ) ii = ii - ( ii / 31328 ) * 31328;
    if(jj > 30081 ) jj = jj - ( jj / 30081 ) * 30081;
    itest=itest+1;
    if((( ii > 0 ) &&  ( ii <= 31328 )) && (( jj > 0 ) && 
					    ( jj <= 30081 ))){
      i1=ii/177+2; i2=ii-(i1-2)*177+2; i3=jj/169+1; i4=jj-(i3-1)*169; 
      i4 = jj - (i3-1)*169;
      count1=0;
      while ( count1 < 97 ){
	s=0.;
	t=0.5;
	count2=0;
	while( count2 < 24 ){
	  idum=i1*i2/179;
	  idum=( i1*i2 - (i1*i2/179)*179 ) * i3;
	  i5=idum-(idum/179)*179;
	  i1=i2; i2=i3; i3=i5; idum=53*i4+1; i4=idum-(idum/169)*169;
	  if( i4*i5-((i4*i5)/64)*64 >= 32 ) s=s+t;
	  t=0.5*t;
	  count2=count2+1;
	}
	u[count1] = s;
	count1 = count1 +1;
      }
      c = 362436./16777216.;  cd = 7654321./16777216.; 
      cm = 16777213./16777216.;
    }
    else{
      cout << " wrong initialization " << endl;
    }
  }
  else{
    uni = u[i] - u[j]; if( uni < 0.) uni = uni + 1; u[i] = uni; 
    i = i -1;
    if( i < 0 ) i = 96; j = j - 1; if ( j < 0 ) j = 96; c = c - cd;
    if( c < 0. ) c = c+cm; uni = uni-c ; if( uni < 0. )uni = uni+1.;
    //    return uni;
  }
  return uni;

}   

void AliPMDClustering::SetEdepCut(Float_t decut)
{
  fCutoff = decut;
}
void AliPMDClustering::SetDebug(Int_t idebug)
{
  fDebug = idebug;
}
