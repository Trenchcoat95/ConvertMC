#pragma once
#ifndef HELIXOLD_H_
#define HELIXOLD_H_
#include <iostream>
#include "TVectorD.h"
#include "TMatrix.h"
#include "TMath.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TF2.h"
#include "Math/Vector3D.h"
#include <algorithm>
#include <vector>
#include <ctime>
using namespace ROOT::Math;

void sort_TPCClusters_along_track(const std::vector<XYZVector>  &TPCClusters,
                                            std::vector<int> &hlf,
                                            std::vector<int> &hlb,
                                            int printlevel,
                                            float &lengthforwards,
                                            float &lengthbackwards,
                                            float /* sorttransweight */,
                                            float /* sortdistback */)
{

  // this sorting code appears in the patrec stage too, if only to make tracks that can be drawn
  // on the event display

  // find candidate endpoints and sort TPCClusters

  float cmin[3];  // min x, y, and z coordinates over all TPCClusters
  float cmax[3];  // max x, y, and z coordinates over all TPCClusters
  size_t ihex[6];  // index of TPCCluster which gave the min or max ("extreme") 0-2: (min xyz)  3-5 (max xyz)


  for (size_t iTPCCluster=0; iTPCCluster < TPCClusters.size(); ++iTPCCluster)
    {
      for (int i=0; i<3; ++i)
        {
          float c ;
          
          if(i==0) c = TPCClusters[iTPCCluster].X();
          else if (i==1) c = TPCClusters[iTPCCluster].Y();
          else if (i==2) c = TPCClusters[iTPCCluster].Z();
          if(printlevel>3) std::cout<<"i: "<<i<<" c: "<<c<<std::endl;
  
          
          
          if (iTPCCluster==0)
            {
              cmin[i] = c;
              cmax[i] = c;
              ihex[i] = 0;
              ihex[i+3] = 0;
            }
          else
            {
              if (c<cmin[i])
                {
                  cmin[i] = c;
                  ihex[i] = iTPCCluster;
                }
              if (c>cmax[i])
                {
                  cmax[i] = c;
                  ihex[i+3] = iTPCCluster;
                }
            }
        }
    }

  if(printlevel>3)
  {
    std::cout<<"ihex: ";
    for(int ki=0;ki<6;ki++) std::cout<<ihex[ki]<<" ";
    std::cout<<std::endl;
  }


  // now we have six TPCClusters that have the min and max x, y, and z values.  Find out which of these six
  // TPCClusters has the biggest sum of distances to all the other TPCClusters (the most extreme)
  float sumdmax = 0;
  size_t imax = 0;
  for (size_t i=0; i<6; ++i)
    {
      float sumd = 0;
      XYZVector poshc = TPCClusters[ihex[i]];
      for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
          XYZVector hp = TPCClusters[iTPCCluster]; 
          sumd += sqrt((poshc - hp).Mag2());
        }
      if (sumd > sumdmax)
        {
          sumdmax = sumd;
          imax = i;
        }
    }
  
  if(printlevel>3)
  {
    std::cout<<"imax: "<<imax<<std::endl;
  }


  //  Use this TPCCluster, indexed by ihex[imax] as a starting point
  //  For cases with less than dmindir length, sort TPCClusters in order of how
  //  far they are from the first TPCCluster.  If longer, calculate a curdir vector
  //  and add hits that are a minimum step along curdir

  hlf.clear();
  lengthforwards = 0;
  hlf.push_back(ihex[imax]);
  XYZVector lpos = TPCClusters[hlf[0]];
  XYZVector lastpos=lpos;
  XYZVector cdpos=lpos;  // position of the beginning of the direction vector
  size_t cdposindex=0;  // and its index in the hlf vector
  XYZVector hpadd;
  XYZVector curdir(0,0,0);
  bool havecurdir = false;
  float dmindir=10;       // promote to a fcl parameter someday.  Distance over which to compute curdir

  for (size_t inh=1;inh<TPCClusters.size();++inh)
    {
      float dmin=0;
      int jmin=-1;
      for (size_t jh=0;jh<TPCClusters.size();++jh)
        {
          bool found = false;
          for (size_t kh=0;kh<hlf.size();++kh)
            {
              if (hlf[kh] == (int) jh)
                {
                  found = true;
                  break;
                }
            }
          if (found) continue;   // skip if we've already assigned this TPCCluster on this track
          XYZVector hpos = TPCClusters[jh];

          float d=0;
          if (havecurdir)
            {
              float dtmp = (hpos-lastpos).Dot(curdir);
              if (dtmp<-2)
                {
                  dtmp = -4 - dtmp;
                }
              d= dtmp + 0.1 * sqrt(((hpos-lastpos).Cross(curdir)).Mag2());
            }
          else
            {
              d=sqrt((hpos-lpos).Mag2());
            }
          if (jmin == -1)
            {
              jmin = jh;
              dmin = d;
              hpadd = hpos;
            }
          else
            {
              if (d<dmin)
                {
                  jmin = jh;
                  dmin = d;
                  hpadd = hpos;
                }
            }
        }
      //  std::cout << "dmin: " << dmin << std::endl;
      hlf.push_back(jmin);
      lengthforwards += sqrt((hpadd-lastpos).Mag2());   // add up track length
      lastpos = hpadd;
      if ( sqrt((hpadd-cdpos).Mag2()) > dmindir)
        {
          XYZVector cdcand = TPCClusters[hlf[cdposindex]];
          size_t itctmp = cdposindex;
          for (size_t itc=cdposindex+1; itc<hlf.size(); ++itc)
            {
              XYZVector cdcandt=TPCClusters[hlf[itc]];
              if (  sqrt((hpadd-cdcandt).Mag2())  < dmindir  )
                {
                  break;
                }
              else
                {
                  itctmp = itc;
                  cdcand = cdcandt;
                }
            }
          //std::cout << "Updated curdir: " << cdposindex << " " << itctmp << " " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << " ";
          cdposindex = itctmp;
          cdpos = cdcand;
          curdir = hpadd - cdcand;
          curdir *= (1.0/sqrt(curdir.Mag2()));
          //std::cout << "To: " << curdir.X() << " " << curdir.Y() << " " << curdir.Z() << std::endl;
          havecurdir = true;
        }
    }

  if (printlevel>2)
    {
      for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
          printf("Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
                 (int) iTPCCluster,
                 TPCClusters[iTPCCluster].X(),
                 TPCClusters[iTPCCluster].Y(),
                 TPCClusters[iTPCCluster].Z(),
                 hlf[iTPCCluster],
                 TPCClusters[hlf[iTPCCluster]].X(),
                 TPCClusters[hlf[iTPCCluster]].Y(),
                 TPCClusters[hlf[iTPCCluster]].Z());
        }
    }

  // now go backwards -- just invert the sort order

  hlb.clear();
  for (size_t i=0; i< hlf.size(); ++i)
    {
      hlb.push_back(hlf[hlf.size()-1-i]);  // just invert the order for the backward sort
    }
  lengthbackwards = lengthforwards;

  if (printlevel>2)
    {
      for (size_t iTPCCluster=0; iTPCCluster<TPCClusters.size(); ++iTPCCluster)
        {
          printf("Backward Sort compare: %5d %10.3f %10.3f %10.3f  %5d %10.3f %10.3f %10.3f\n",
                 (int) iTPCCluster,
                 TPCClusters[iTPCCluster].X(),
                 TPCClusters[iTPCCluster].Y(),
                 TPCClusters[iTPCCluster].Z(),
                 hlb[iTPCCluster],
                 TPCClusters[hlb[iTPCCluster]].X(),
                 TPCClusters[hlb[iTPCCluster]].Y(),
                 TPCClusters[hlb[iTPCCluster]].Z());
        }
    }
}


double capprox(double x1,double y1,
                        double x2,double y2,
                        double x3,double y3,
                        double &xc, double &yc)
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature -- copied from ALICE
  // here x is y and y is z for us
  //-----------------------------------------------------------------
  x3 -=x1;
  x2 -=x1;
  y3 -=y1;
  y2 -=y1;
  //
  double det = x3*y2-x2*y3;
  if (TMath::Abs(det)<1e-10){
    return 100;
  }
  //
  double u = 0.5* (x2*(x2-x3)+y2*(y2-y3))/det;
  double x0 = x3*0.5-y3*u;
  double y0 = y3*0.5+x3*u;
  double c2 = 1/TMath::Sqrt(x0*x0+y0*y0);
  xc = x0 + x1;
  yc = y0 + y1;
  if (det<0) c2*=-1;
  return c2;
}

void ouchef(double *x, double *y, int np, double &xcirccent, double &ycirccent,
                      double &rcirc, double &chisq, int &ifail)
{

  // converted from the original Fortran by T. Junk and de-OPALified -- apologies for the gotos

  //      SUBROUTINE OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *
  // *...OUCHEF  circle fit by linear regression method
  // *.
  // *.            OUCHEF(X,Y,NP,XCIRCCENT,YCIRCCENT,RCIRC,CHISQ,IFAIL)
  // *.
  // *.  INPUT
  // *.       X(.)      x coordinates of points
  // *.       Y(.)      y coordinates of points
  // *.       NP        number of points
  // *.
  // *.  OUTPUT
  // *.       XMED, YMED   center of circle
  // *.       R            radius of circle
  // *.       CHISQ     residual sum of squares : to be a chi**2 ,it should
  // *.                 be divided by (NP-3)*(mean residual)**2
  // *.
  // *.       IFAIL     error flag : 0   OK
  // *.                              1   less than 3 points
  // *.                              2   no convergence
  // *.                              5   number of points exceeding
  // *.                                  arrays dimensions -> fit performed
  // *.                                  with allowed maximum (NPTMAX)
  // *.
  // *. AUTHOR :  N.Chernov, G.Ososkov
  // *.
  // *.     written by N. Chernov and G. Ososkov, Dubna,
  // *.     reference:  computer physics communications vol 33, p329
  // *.     adapted slightly by F. James, March, 1986
  // *.     further adapted slightly by M.Hansroul May,1989
  // *.
  // *. CALLED   : any
  // *.
  // *. CREATED  : 31 May 1989
  // *. LAST MOD : 11 January 1996
  // *.
  // *. Modification Log.:
  // *. 11-jan-96  M.Hansroul  preset chisq to 9999.
  // *. 12-oct-95  M.Hansroul  clean up
  // *. 12-jul-93  M.Hansroul  protection against huge d0
  // *.                        put OUATG2 in line
  // *. 31-May-89  M.Hansroul  use OPAL output parameters
  // *.*********************************************************************

  const double eps = 1.0e-10;
  const int itmax=20;
  const int nptmax=1000;
  //const double twopi=8.0*atan(1.0);
  //const double piby2=2.0*atan(1.0);

  int nptlim=0;
  int npf=0;
  int n=0;
  int ihalf=0;
  int i=0;
  int iter=0;
  int iswap=0;
  //double sfi0=0;
  //double cfi0=0;
  double xmid=0;
  double ymid=0;
  double *xf;
  double *yf;
  double vlf[2];
  double vmf[2];
  //double vlm[2];

  double alf=0;
  double alm=0;
  double a0=0;
  double a1=0;
  double a2=0;
  double a22=0;
  double bem=0;
  double bet=0;
  double cur=0;
  double cu2=0;
  double dd=0;
  double den=0;
  double det=0;
  double dy=0;
  double d2=0;
  double f=0;
  double fact=0;
  double fg=0;
  double f1=0;
  double g=0;
  double gam=0;
  double gam0=0;
  double gmm=0;
  double g1=0;
  double h=0;
  double h2=0;
  double p2=0;
  double q2=0;
  double rm=0;
  double rn=0;
  double xa=0;
  double xb=0;
  double xd=0;
  double xi=0;
  double xm=0;
  double xx=0;
  double xy=0;
  double x1=0;
  double x2=0;
  double ya=0;
  double yb=0;
  double yd=0;
  double yi=0;
  double ym=0;
  double yy=0;
  double y1=0;
  double y2=0;
  double xcd=0;
  double ycd=0;
  //double xcycsq=0;
  //double rcd;

  //*     -----------------------------------------------------------------
  //*     xmid,ymid       = centre of fitted circle of radius radfit
  //*     chisq           = residual sum of squares per point
  //*     alf,bet,gam,cur = internal parameters of the circle

  ifail = 0;
  chisq = 9999.;

  nptlim = np;
  if (nptlim>nptmax)
    {
      nptlim = nptmax;
      ifail  = 5;
      return;
    }

  npf = nptlim + 1;
 label_12:
  npf = npf - 1;
  if (npf < 2)
    {
      ifail = 1;
      return;
    }
  vlf[0] = x[npf-1] - x[0];
  vlf[1] = y[npf-1] - y[0];

  //*                get a point near the middle

  ihalf  = (npf+1)/2;

  vmf[0] = x[ihalf-1] - x[0];
  vmf[1] = y[ihalf-1] - y[0];
  //vlm[0] = x[npf-1] - x[ihalf-1];
  //vlm[1] = y[npf-1] - y[ihalf-1];

  //*            check for more than half a cicle of points

  if (vlf[0]*vmf[0] + vlf[1]*vmf[1] < 0.) goto label_12;

  //*            if the track is nearer to the y axis
  //*            interchange x and y

  if (fabs(y[npf-1] - y[0]) > fabs(x[npf-1]-x[0]))
    {
      yf = x;
      xf = y;
      iswap = 1;
    }
  else
    {
      xf = x;
      yf = y;
      iswap = 0;
    }

  n  = npf;
  rn = 1./((double) n);
  xm = 0.;
  ym = 0.;
  for (i=0;i<n;++i)
    {
      xm = xm + xf[i];
      ym = ym + yf[i];
    }

  xm = xm * rn;
  ym = ym * rn;
  x2 = 0.;
  y2 = 0.;
  xy = 0.;
  xd = 0.;
  yd = 0.;
  d2 = 0.;
  for (i=0; i<n; ++i)
    {
      xi = xf[i] - xm;
      yi = yf[i] - ym;
      xx = xi*xi;
      yy = yi*yi;
      x2 = x2 + xx;
      y2 = y2 + yy;
      xy = xy + xi*yi;
      dd = xx + yy;
      xd = xd + xi*dd;
      yd = yd + yi*dd;
      d2 = d2 + dd*dd;
    }

  x2   = x2*rn;
  y2   = y2*rn;
  xy   = xy*rn;
  d2   = d2*rn;
  xd   = xd*rn;
  yd   = yd*rn;
  f    = 3.*x2 + y2;
  g    = 3.*y2 + x2;
  fg   = f*g;
  h    = xy + xy;
  h2   = h*h;
  p2   = xd*xd;
  q2   = yd*yd;
  gam0 = x2 + y2;
  fact = gam0*gam0;
  a2   = (fg-h2-d2)/fact;
  fact = fact*gam0;
  a1   = (d2*(f+g) - 2.*(p2+q2))/fact;
  fact = fact*gam0;
  a0   = (d2*(h2-fg) + 2.*(p2*g + q2*f) - 4.*xd*yd*h)/fact;
  a22  = a2 + a2;
  yb   = 1.0e30;
  iter = 0;
  xa   = 1.;
  //*                main iteration
 label_3:
  ya = a0 + xa*(a1 + xa*(a2 + xa*(xa-4.0)));
  if (fabs(ya) > fabs(yb))  goto label_4;
  if (iter >= itmax)  goto label_4;
  dy = a1 + xa*(a22 + xa*(4.*xa - 12.));
  xb = xa - ya/dy;
  if (fabs(xa-xb) <= eps)  goto label_5;
  xa = xb;
  yb = ya;
  iter = iter + 1;
  goto label_3;

 label_4:
  ifail = 2;
  //*      print 99, iter
  //*   99 format ('  circle - no convergence after',i4,' iterations.')
  xb = 1.;

 label_5:
  gam   = gam0*xb;
  f1    = f - gam;
  g1    = g - gam;
  x1    = xd*g1 - yd*h;
  y1    = yd*f1 - xd*h;
  det   = f1*g1 - h2;
  den   = 1./sqrt(x1*x1 + y1*y1 + gam*det*det);
  cur   = det*den;
  alf   = -(xm*det + x1)*den;
  bet   = -(ym*det + y1)*den;
  rm    = xm*xm + ym*ym;
  gam   = ((rm-gam)*det + 2.*(xm*x1 + ym*y1))*den*0.5;
  alm   = alf + cur*xm;
  bem   = bet + cur*ym;

  gmm   = gam + alf*xm + bet*ym + 0.5*cur*rm;
  chisq = ((0.5*cur)*(0.5*cur)*d2 + cur*(alm*xd + bem*yd + gmm*gam0)
           + alm*alm*x2 + bem*bem*y2 + 2.*alm*bem*xy + gmm*gmm) /rn;
  cu2   = 0.;
  if (cur != 0.)
    {
      cu2    =  1./cur;
      xcd    = -alf*cu2;
      ycd    = -bet*cu2;
      //rcd    =  fabs(cu2);
      xmid   =  xcd;
      ymid   =  ycd;
      xcirccent = xmid;
      ycirccent = ymid;
      if (iswap == 1)
        {
          xcirccent = ymid;
          ycirccent = xmid;
        }
      rcirc = cu2;
    }
  else
    {
      rcirc = 1e6;
    }

  return;
}


int Helix_Fit(const std::vector<XYZVector>  TPCClusters,
              XYZVector  &TPCClustersSeed,
              double &curvature_init,
              double &lambda_init,
              double &phi_init,
              double dir,
              int printlevel)
{

  size_t InitialTPNTPCClusters = 100;
  size_t nTPCClusters = TPCClusters.size();
  size_t firstTPCCluster = 0;
  size_t farTPCCluster = TMath::Min(nTPCClusters-1, InitialTPNTPCClusters);
  size_t intTPCCluster = farTPCCluster/2;
  size_t lastTPCCluster = nTPCClusters-1;
  

  double trackbeg[3] = {TPCClusters.at(firstTPCCluster).X(),
                       TPCClusters.at(firstTPCCluster).Y(),
                       TPCClusters.at(firstTPCCluster).Z()};

  double tp1[3] = {TPCClusters.at(intTPCCluster).X(),
                  TPCClusters.at(intTPCCluster).Y(),
                  TPCClusters.at(intTPCCluster).Z()};

  double tp2[3] = {TPCClusters.at(farTPCCluster).X(),
                  TPCClusters.at(farTPCCluster).Y(),
                  TPCClusters.at(farTPCCluster).Z()};
    
  if (printlevel>1)
    {
      std::cout << "TPCCluster Dump in initial_trackpar_estimate: " << std::endl;
      for (size_t i=0;i<nTPCClusters;++i)
        {
          size_t ihf = i;
          std::cout << i << " : " <<
            TPCClusters.at(ihf).X() << " " <<
            TPCClusters.at(ihf).Y() << " " <<
            TPCClusters.at(ihf).Z() << std::endl;
        }
    }
  if (printlevel>0)
    {
      std::cout << "First TPCCluster x, y, z: " << trackbeg[0] << " " << trackbeg[1] << " " << trackbeg[2] << std::endl;
      std::cout << "Inter TPCCluster x, y, z: " << tp1[0] << " " << tp1[1] << " " << tp1[2] << std::endl;
      std::cout << "Far   TPCCluster x, y, z: " << tp2[0] << " " << tp2[1] << " " << tp2[2] << std::endl;
    }

  TPCClustersSeed.SetXYZ(trackbeg[0],trackbeg[1],trackbeg[2]);
  double x_other_end = TPCClusters.at(lastTPCCluster).X();


  double ycc=0;
  double zcc=0;
  curvature_init = capprox(trackbeg[1],trackbeg[2],tp1[1],tp1[2],tp2[1],tp2[2],ycc,zcc);
  //std::cout << " inputs to trackpar circ fit (y,z): " << trackbeg[1] << " " << trackbeg[2] << " : "
  //        << tp1[1] << " " << tp1[2] << " : " << tp2[1] << " " << tp2[2] << std::endl;
  //std::cout << "curvature output: " << curvature_init << std::endl;

  // redo calc with a circle fit to all points, not just three

  double dycc=0;
  double dzcc=0;
  double drcc=0;
  double dchisq=0;
  int ifail=0;
  std::vector<double> dtpcclusy;
  std::vector<double> dtpcclusz;
  for (size_t i=0;i<nTPCClusters; ++i)
    {
      dtpcclusy.push_back(TPCClusters.at(i).Y());
      dtpcclusz.push_back(TPCClusters.at(i).Z());
    }
  int ict = nTPCClusters;
  ouchef(dtpcclusy.data(),dtpcclusz.data(),ict,dycc,dzcc,drcc,dchisq,ifail);
  if (ifail == 0 && drcc != 0)
    {
      ycc = dycc;
      zcc = dzcc;
      if (curvature_init < 0) drcc = -std::abs(drcc);
      //if (ict > 3)
      //  {
      //    std::cout << "updating curvature " << ict << " " << curvature_init << " with " << 1.0/drcc << std::endl;
      //  }
      curvature_init = 1.0/drcc;
    }


  phi_init = TMath::ATan2( trackbeg[2] - zcc, ycc - trackbeg[1] );
  double phi2 = phi_init;
  if (curvature_init<0) phi_init += TMath::Pi();
  double radius_init = 10000;
  if (curvature_init != 0) radius_init = 1.0/curvature_init;

  double dx1 = tp2[0] - TPCClustersSeed.X();
  if (dx1 != 0)
    {
      double dphi2 = TMath::ATan2(tp2[2]-zcc,ycc-tp2[1])-phi2;
      if (dphi2 > TMath::Pi()) dphi2 -= 2.0*TMath::Pi();
      if (dphi2 < -TMath::Pi()) dphi2 += 2.0*TMath::Pi();
      lambda_init = TMath::ATan(1.0/((radius_init/dx1)*dphi2));
    }
  else
    {
      //std::cout << "initial track par estimate failure" << std::endl;
      lambda_init = 0;
      return 0;
    } // got fMinNumTPCClusters all at exactly the same value of x (they were sorted).  Reject track.

  if (printlevel>0)
    {
      std::cout << "phi calc: dz, dy " << tp2[2]-trackbeg[2] << " " <<  tp2[1]-trackbeg[1] << std::endl;
      std::cout << "initial curvature, phi, lambda: " << curvature_init << " " << phi_init << " " << lambda_init << std::endl;
    }
  return 1;

}



void CalculateCovariance(const std::vector<XYZVector>  TPCClusters,
              TMatrixD &P,
              double curvature_init,
              double lambda_init,
              double phi_init,
              double dir,
              int printlevel,
              double sxy
              )
{
  size_t InitialTPNTPCClusters = 100;
  size_t nTPCClusters = TPCClusters.size();
  size_t firstTPCCluster = 0;
  size_t farTPCCluster = TMath::Min(nTPCClusters-1, InitialTPNTPCClusters);
  size_t intTPCCluster = farTPCCluster/2;
  size_t lastTPCCluster = nTPCClusters-1;
  

  std::vector<XYZVector> xyz3points;
  XYZVector xyzpoint;

  

  xyzpoint.SetXYZ(TPCClusters.at(firstTPCCluster).X(),
                       TPCClusters.at(firstTPCCluster).Y(),
                       TPCClusters.at(firstTPCCluster).Z());

  xyz3points.push_back(xyzpoint);                    

  xyzpoint.SetXYZ(TPCClusters.at(intTPCCluster).X(),
                  TPCClusters.at(intTPCCluster).Y(),
                  TPCClusters.at(intTPCCluster).Z());

  xyz3points.push_back(xyzpoint);

  xyzpoint.SetXYZ(TPCClusters.at(farTPCCluster).X(),
                  TPCClusters.at(farTPCCluster).Y(),
                  TPCClusters.at(farTPCCluster).Z());

  xyz3points.push_back(xyzpoint);

  
  double curvature_smear, lambda_smear, phi_smear;

  /////
  xyz3points[2].SetY(xyz3points[2].Y()+sxy);

  Helix_Fit(xyz3points,xyzpoint,curvature_smear,lambda_smear,phi_smear,dir,printlevel);

  Double_t f40=(curvature_smear-curvature_init)/sxy;
  Double_t f23=(phi_smear-phi_init)/sxy;
  Double_t f30=(lambda_smear-lambda_init)/sxy;

  xyz3points[2].SetY(xyz3points[2].Y()-sxy);
  //////



  /////
  xyz3points[1].SetY(xyz3points[1].Y()+sxy);

  Helix_Fit(xyz3points,xyzpoint,curvature_smear,lambda_smear,phi_smear,dir,printlevel);

  Double_t f42=(curvature_smear-curvature_init)/sxy;
  Double_t f22=(phi_smear-phi_init)/sxy;
  Double_t f32=(lambda_smear-lambda_init)/sxy;

  xyz3points[1].SetY(xyz3points[1].Y()-sxy);
  //////


  /////
  xyz3points[0].SetY(xyz3points[0].Y()+sxy);

  Helix_Fit(xyz3points,xyzpoint,curvature_smear,lambda_smear,phi_smear,dir,printlevel);

  Double_t f43=(curvature_smear-curvature_init)/sxy;
  Double_t f20=(phi_smear-phi_init)/sxy;

  xyz3points[0].SetY(xyz3points[0].Y()-sxy);
  //////

  /////
  xyz3points[2].SetX(xyz3points[2].X()+sxy);

  Helix_Fit(xyz3points,xyzpoint,curvature_smear,lambda_smear,phi_smear,dir,printlevel);

  Double_t f31=(lambda_smear-lambda_init)/sxy;

  xyz3points[2].SetY(xyz3points[2].X()-sxy);
  //////

  /////
  xyz3points[1].SetX(xyz3points[1].X()+sxy);

  Helix_Fit(xyz3points,xyzpoint,curvature_smear,lambda_smear,phi_smear,dir,printlevel);

  Double_t f34=(lambda_smear-lambda_init)/sxy;

  xyz3points[1].SetY(xyz3points[1].X()-sxy);
  //////


  double sxy2 = sxy*sxy;
  P[0][0]=sxy2;
  P[1][1]=sxy2;
  P[2][2]=f20*sxy2*f20+f22*sxy2*f22+f23*sxy2*f23;
  P[3][3]=f30*sxy2*f30+f31*sxy2*f31+f32*sxy2*f32+f34*sxy2*f34; 
  P[4][4]=f40*sxy2*f40+f42*sxy2*f42+f43*sxy2*f43;

  P[2][2]*=cos(phi_init)*cos(phi_init);//=sin(sqrt(P[2][2]))*sin(sqrt(P[2][2]));
  P[3][3]/=cos(lambda_init)*cos(lambda_init);//=tan(sqrt(P[3][3]))*tan(sqrt(P[3][3]));
  P[4][4]/=(0.5*0.3e-2)*(0.5*0.3e-2); // transform to 1/pt
  //P[4][4]*=0.2*0.2;
  

}


#endif