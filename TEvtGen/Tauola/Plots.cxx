#include "Plots.h"
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

namespace Tauolapp
{

Plots::Plots():
    m_incoming_pdg_id(1),
    m_cosTheta       (-0.2),
    m_n_plot_points  (1000)
{
}

void Plots::SANCtest1(){

  cout<<"SANC plot 1 (short)..."<<endl;

  double smin = log(6.*6.)+0.0001;
  double smax = log(17000.*17000.);
  double step = (smax-smin)/(m_n_plot_points-1);

  ofstream f1,f2,f3;
  f1.open("f-sanc.txt");
  f2.open("f-born.txt");
  f3.open("f-plzap0.txt");
  for(int i=0; i<m_n_plot_points; i++)
  {
    double s = exp(smin+i*step);

    // Write SANC value
    t_pair.recalculateRij(m_incoming_pdg_id,15,s,m_cosTheta);
    f1<<sqrt(s)<<" "<<t_pair.m_R[0][3]<<endl;

    // Write Born-level value
    // (assuming table 11 is filled with born-level data)
    t_pair.recalculateRij(11,15,s,m_cosTheta);
    f2<<sqrt(s)<<" "<<t_pair.m_R[0][3]<<endl;

    int outgoing_pdg_id = 15;

    // Write PLZAP0 value
    double pz = 1-plzap0_(&m_incoming_pdg_id,&outgoing_pdg_id,&s, &m_cosTheta);
    t_pair.m_R[0][3]=2*pz-1;
    f3<<sqrt(s)<<" "<<t_pair.m_R[0][3]<<endl;
  }
  f1.close();
  f2.close();
  f3.close();
}

void Plots::SANCtest2(){

  cout<<"SANC plot 2 (short)..."<<endl;

  double smin = log(6.*6.)+0.0001;
  double smax = log(17000.*17000.);
  double step = (smax-smin)/(m_n_plot_points-1);

  ofstream f1,f2,f3;
  f1.open("f-w-single-point.txt");
  f2.open("f-w0-single-point.txt");
  f3.open("f-ww0-single-point.txt");

  for(int i=0; i<m_n_plot_points; i++){

    double s=exp(smin+i*step);
    t_pair.recalculateRij(m_incoming_pdg_id,15,s,m_cosTheta);

    // Write w, w0 and w/w0
    f1<<sqrt(s)<<" "<<Tauola::getEWwt()<<endl;
    f2<<sqrt(s)<<" "<<Tauola::getEWwt0()<<endl;
    f3<<sqrt(s)<<" "<<Tauola::getEWwt()/Tauola::getEWwt0()<<endl;
  }
  f1.close();
  f2.close();
  f3.close();
}

void Plots::SANCtest3(){

  cout<<"SANC plot 3 (long)..."<<endl;

  double smin = log(6.*6.)+0.0001;
  double smax = log(17000.*17000.);
  double step = (smax-smin)/(m_n_plot_points-1);

  ofstream f1;
  f1.open("f-err.txt");
  double costheta=-1.;

  for(int i=0; i<m_n_plot_points; i++){

    double buf=0.,err=0.;

    for(int j=0; j<m_n_plot_points; j++){

      double s = exp(smin+j*step);

      // Get value from SANC table
      t_pair.recalculateRij(m_incoming_pdg_id,15,s,costheta);
      buf = t_pair.m_R[0][3];
      t_pair.recalculateRij(11,15,s,costheta);

      // Calculate error between SANC and Born-level
      err += (buf-t_pair.m_R[0][3])*(buf-t_pair.m_R[0][3]);
    }

    f1<<costheta<<" "<<err/m_n_plot_points<<endl;
    err=0.;
    costheta+=2./m_n_plot_points;
  }

  f1.close();
}

void Plots::SANCtest4(){

  cout<<"SANC plot 4 (medium)..."<<endl;

  double smin = log(6.*6.);
  double smax = log(17000.*17000.);
  double step = (smax-smin)/(m_n_plot_points-1);

  ofstream f1,f2,f3;
  f1.open("f-cross.txt");
  f2.open("f-w.txt");
  f3.open("f-w0.txt");

  for(int i=0; i<m_n_plot_points; i++){

    double s        =  exp(smin+i*step);
    double sumEWwt  =  0.;
    double sumEWwt0 =  0.;
    double costheta = -1.;
    int    NCOS     =  21;

    // Calculate sum of w/w0 for all cosTheta
    for(int j=0; j<NCOS; j++){

      costheta = -1. + 1.0/NCOS + j*2./NCOS;
      
      t_pair.recalculateRij(m_incoming_pdg_id,15,s,costheta);

      sumEWwt +=Tauola::getEWwt();
      sumEWwt0+=Tauola::getEWwt0();
    }

    f1<<sqrt(s)<<" "<<sumEWwt/sumEWwt0/m_n_plot_points<<endl;
    f2<<sqrt(s)<<" "<< 2./NCOS * sumEWwt  <<endl;
    f3<<sqrt(s)<<" "<< 2./NCOS * sumEWwt0 <<endl;
  }

  f1.close();
  f2.close();
  f3.close();
}

void Plots::setSancVariables(int incoming, double cosTheta) {
  m_incoming_pdg_id = incoming;
  m_cosTheta        = cosTheta;
}

} // namespace Tauolapp
