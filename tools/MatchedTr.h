#ifndef _MATCHEDTR_
#define _MATCHEDTR_

#include "includeROOT.h"

using namespace std;

class MatchedTr {
    int matchid;
    bool MatchedSt;
    float Length;
    float Ionization;
    TVector3 xyz;
    TVector3 pxyz;
    int PDG=0;
    int Q=0;
  public:
    MatchedTr (int,TVector3,TVector3,bool);
    MatchedTr (TVector3,TVector3,float,float);

    int GetID(void) {return matchid;};
    TVector3 GetXYZ(void) {return xyz;};
    TVector3 GetP(void) {return pxyz;};
    float GetLength(void) {return Length;};
    float GetIon(void) {return Ionization;};
    bool GetMatchedSt(void) {return MatchedSt;};
    float GetPDG(void) {return PDG;};
    float GetQ(void) {return Q;};
    void SetPDG(int);
    void SetID(int);
    void SetQ(int);
};

MatchedTr::MatchedTr (int id,TVector3 xyzt,TVector3 pxyzt, bool MatchedStt) {
  xyz = xyzt;
  pxyz = pxyzt;
  matchid = id;
  MatchedSt = MatchedStt;
}

MatchedTr::MatchedTr (TVector3 xyzt,TVector3 pxyzt, float lengtht, float iont) {
  xyz = xyzt;
  pxyz = pxyzt;
  Length=lengtht;
  Ionization=iont;
}

void MatchedTr::SetPDG(int PDGt){PDG=PDGt;};
void MatchedTr::SetID(int IDt){matchid=IDt;};
void MatchedTr::SetQ(int Qt){Q=Qt;};

#endif