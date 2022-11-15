#ifndef garutils_H
#define garutils_H

#include "TTreeStream.h"
#include "TStopwatch.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TPad.h"
#include "TCanvas.h"
#include "AliPID.h"
#include <iostream>
#include <string>
#include <fstream>


void idx2ijsort2(size_t n, size_t k, size_t &i, size_t &j)
{
  if (n<2)
    {
      i=0;
      j=0;
      return;
    }
  if (k > n*(n-1)/2)
    {
      throw std::runtime_error(Form("idx2ijsort2: k too big: %lu %lu",k,n)) ;
    }
  i = n - 2 - TMath::Floor(TMath::Sqrt( (double) (-8*k + 4*n*(n-1)-7) )/2.0 - 0.5);
  j = k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2;
}

size_t ij2idxsort2(size_t n, size_t i_in, size_t j_in)
{
  if (i_in >= n)
    {
      throw std::runtime_error("ij2idxsort2 i_in >= n");
    }
  if (j_in >= n)
    {
      throw std::runtime_error("ij2idxsort2 j_in >= n");
    }
  size_t i = TMath::Min(i_in,j_in);
  size_t j = TMath::Max(i_in,j_in);
  size_t k = (n*(n-1)/2) - (n-i)*((n-i)-1)/2 + j - i - 1;
  return k;
}


bool sort2_check_cyclic(std::vector<int> &link1, std::vector<int> &link2, int i, int j)
{
  if (i == j) return true;

  int prev = i;
  int next = (link1.at(i) == -1) ? link2.at(i) : link1.at(i);
  while (next != -1)
    {
      if (next == j) return true;
      int ptmp = next;
      next = (link1.at(next) == prev) ? link2.at(next) : link1.at(next);
      prev = ptmp;
    }
  return false;
}


void sort_TPCClusters_along_track2(const std::vector<TVector3>  &TPCClusters,
                                        std::vector<int> &hlf,
                                        std::vector<int> &hlb,
                                        int /* printlevel */,
                                        float &lengthforwards,
                                        float &lengthbackwards,
					     float dcut)
{
  size_t nclus = TPCClusters.size();
  hlf.clear();
  hlb.clear();
  lengthforwards = 0;
  lengthbackwards = 0;
  if (nclus == 0) return;
  if (nclus == 1)
    {
      hlf.push_back(0);
      hlb.push_back(0);
      return;
    }

  size_t ndists = nclus*(nclus-1)/2;
  std::vector<float> dists(ndists);
  std::vector<int> dsi(ndists);
  std::vector<int> link1(nclus,-1);
  std::vector<int> link2(nclus,-1);
  for (size_t i=0; i<nclus-1; ++i)
    {
      TVector3 iclus(TPCClusters.at(i));
      for (size_t j=i+1; j<nclus; ++j)
        {
          TVector3 jclus(TPCClusters.at(j));
          float d = (iclus-jclus).Mag();
          dists.at(ij2idxsort2(nclus,i,j)) = d;
        }
    }
  TMath::Sort( (int) ndists, dists.data(), dsi.data(), false );
  //std::cout << "Dists sorting: " << std::endl;
  //for (size_t i=0; i<ndists; ++i)
  //  {
  //    std::cout << i << " " << dists[i] << " " << dists[dsi[i]] << std::endl;
  //  }

  // greedy pairwise association -- may end up leaving some hits stranded

  for (size_t k=0; k<ndists; ++k)
    {
      size_t i = 0;
      size_t j = 0;
      size_t k2 = dsi.at(k);
      if (dists.at(k2) > dcut) break;
      idx2ijsort2(nclus,k2,i,j);
      TVector3 iclus(TPCClusters.at(i));
      TVector3 jclus(TPCClusters.at(j));

      if (link1.at(i) == -1)
        {
          if (link1.at(j) == -1)
            {
              link1.at(i) = j;
              link1.at(j) = i;
            }
          else  // already have a link on j'th cluster, add a second one only if it's on the other side
            {
              if (link2.at(j) == -1)
                {
                  TVector3 j2clus(TPCClusters.at(link1.at(j)));

		  // old test on angle instead of cyclicness
                  //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::Pi()/2)

		  if (!sort2_check_cyclic(link1,link2,j,i))
                    {
                      link1.at(i) = j;
                      link2.at(j) = i;
                    }
                }
            }
        }
      else
        {
          if (link2.at(i) == -1)
            {
              if (link1.at(j) == -1)
                {
                  link2.at(i) = j;
                  link1.at(j) = i;
                }
              else  // already have a link on j'th cluster, add a second one only if it's on the other side
                {
                  if (link2.at(j) == -1)
                    {
                      TVector3 j2clus(TPCClusters.at(link1.at(j)));

		      // old test on angle instead of cyclicness
                      //if ( (j2clus-jclus).Angle(iclus-jclus) > TMath::Pi()/2)

		      if (!sort2_check_cyclic(link1,link2,j,i))
                        {
                          link2.at(i) = j;
                          link2.at(j) = i;
                        }
                    }
                }
            }
        }
    }

  //std::cout << "sorting clusters and links" << std::endl;
  //for (size_t i=0; i<nclus; ++i)
  //  {
  //    std::cout << i << " " << link1.at(i) << " " << link2.at(i) << " " << 
  //	TPCClusters.at(i).Position()[0] << " " <<
  // 	TPCClusters.at(i).Position()[1] << " " <<
  //	TPCClusters.at(i).Position()[2] << std::endl;
  //}

  std::vector<int> used(nclus,-1);
  std::vector<std::vector<int> > clistv;

  // find the first tpc cluster with just one link

  int ifirst = -1;
  for (size_t i=0; i<nclus; ++i)
    {
      int il1 = link1.at(i);
      int il2 = link2.at(i);

      if ( (il1 == -1 && il2 != -1) || (il1 != -1 && il2 == -1) ) 
	{
	  ifirst = i;
	  break;
	}
    }
  // this can happen if all TPC clusters are singletons -- no links anywhere.  Happens if
  // there are no TPC clusters at all or just one.
  if (ifirst == -1)
    {
      return;
    }

  // follow the chains

  for (size_t i=ifirst; i<nclus; ++i)
    {
      if (used.at(i) == 1) continue;
      std::vector<int> clist;
      clist.push_back(i);
      used.at(i) = 1;

      int idx = i;
      int il1 = link1.at(idx);
      int il2 = link2.at(idx);
      
      int u1 = -2;
      if (il1 != -1) u1 = used.at(il1);
      int u2 = -2;
      if (il2 != -1) u2 = used.at(il2);

      while (u1 == -1 || u2 == -1)
	{
	  idx =  (u1 == -1) ?  il1 : il2;
	  clist.push_back(idx);
	  used.at(idx) = 1;
          il1 = link1.at(idx);
          il2 = link2.at(idx);
          u1 = -2;
          if (il1 != -1) u1 = used.at(il1);
          u2 = -2;
          if (il2 != -1) u2 = used.at(il2);
	}
      clistv.push_back(clist);
    } 

  // temporary solution -- think about what to do.  Report the longest chain

  size_t msize=0;
  int ibest=-1;
  for (size_t ichain=0; ichain<clistv.size(); ++ichain)
    {
      size_t cursize = clistv.at(ichain).size();
      if (cursize > msize)
	{
	  ibest = ichain;
	  msize = cursize;  
	}
    }
  if (ibest == -1) return;

  for (size_t i=0; i<clistv.at(ibest).size(); ++i)
    {
      hlf.push_back(clistv[ibest].at(i));
      if (i>0)
	{
	  TVector3 lastpoint(TPCClusters.at(hlf.at(i-1)));
	  TVector3 curpoint(TPCClusters.at(hlf.at(i)));
	  lengthforwards += (curpoint-lastpoint).Mag();
	}
    } 
  for (size_t i=0; i< hlf.size(); ++i)
    {
      hlb.push_back(hlf[hlf.size()-1-i]);  // just invert the order for the backward sort
    }
  lengthbackwards = lengthforwards;
 
}


#endif