// A selection of includes & inlines for analysis



#ifndef includeROOT_h
#define includeROOT_h



#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <sys/stat.h>
#include <ctime>
#include <string>
using std::cout;		using std::endl;
using std::string;		using std::vector;
#define  COUT   std::cout<<
#define  TABL   <<"\t"<<
#define  ENDL   <<std::endl



// Includes for common ROOT things
#include "Rtypes.h"         // Look in $ROOT_DIR/source/root-<version>/core/base/inc
#include "TMath.h"          // Lacks useful chisquared_cdf_c function
#include "Math/ProbFuncMathCore.h"
#include "TVector3.h"
#include "TCut.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TMinuit.h"
#include "TLorentzVector.h"
#include "TRandom3.h"



// Math made E-Z.  Too E-Z
#define BIG       1.79769313486231570814527423731704e+308
#define BIGINT    2147483647
#define BIGSIZE   18446744073709551615
#define BIGFLOAT  3.40282346638528859811704183484517e+38
#define TINY      2.22044604925031308084726333618164e-16
#define TINYFLOAT 1.1920928955078125e-07

#define MIN(X,Y) ( (X)<(Y) ? (X) : (Y) )
#define MAX(X,Y) ( (X)>(Y) ? (X) : (Y) )

// Why isn't this in cmath?
#define M_PI 3.14159265358979323846
// std::min and std::max are in <algorithm>
#define NINT(a)	 ( int(floor((a)+0.5)) )
#define SGN(a)   ( (a)<0 ? -1 : ( (a)>0 ? +1 : 0 ) )

inline int			SQR(int b)			{return b*b;}
inline float		SQR(float b)		{return b*b;}
inline double		SQR(double b)		{return b*b;}
inline long double	SQR(long double b)	{return b*b;}

inline int			CUBE(int b)			{return b*b*b;}
inline float		CUBE(float b)		{return b*b*b;}
inline double		CUBE(double b)		{return b*b*b;}
inline long double	CUBE(long double b)	{return b*b*b;}

inline int			QUAD(int b)			{int		 p = b*b; return p*p;}
inline float		QUAD(float b)		{float		 p = b*b; return p*p;}
inline double		QUAD(double b)		{double		 p = b*b; return p*p;}
inline long double	QUAD(long double b)	{long double p = b*b; return p*p;}



// std::hypot appeared in C++17.
inline double Qadd(double const x, double const y){
	return std::hypot(x,y);
}
inline double Qadd(double const x, double const y, double const z){
    double absX,absY,absZ;
    absX = fabs(x);     absY = fabs(y);     absZ = fabs(z);
    if ( absX>absY && absX>absZ ) return std::hypot(x,std::hypot(y,z));
    if ( absY>absX && absY>absZ ) return std::hypot(y,std::hypot(x,z));
    return Qadd(z,Qadd(x,y));
}
inline double Qsub(double const x, double const y) {
	if (std::fabs(x) < std::fabs(y)) {
		std::cout << "Fatality in Qsub(double,double): Invalid arguments";
		std::cout << std::endl;
		std::exit(1);
	}
	long double x2 = x*x;		long double y2 = y*y;
	return double(std::sqrt( x2 - y2 )); 
}



// Primitive earth algorithm for uncertainty in binomial-distributed variables
inline double binerr(double const passed, double const tested) {
    // Garbage input
    if ((passed<0)   || (tested<=0) || (passed>tested)) {
		std::cout << "Fatality in binerr(double,double): Invalid arguments ";
		std::cout << passed << "/" << tested << std::endl;
		std::exit(1);
	}
	double err;
	// Fairly ad-hoc statistical uncertainties for low stats
	if ( tested==1.0 ) {
		err = 0.393;
	} else if ( passed==0.0 || passed==tested ) {
		err = 0.421 / pow(tested,0.94);
	} else {
		err = std::sqrt( passed*(tested-passed)/tested );
	}
	err /= tested;
    return err;
}



// A string for std::cout that is a ratio and its binomial uncertainty
char* stringEffErr(int num, int den) {
	// Yea, its a memory leak.  Fight me.
	char* temp = new char[32];
	double rat = 1000.0*double(num)/double(den);
	double err = 1000.0*binerr(double(num),double(den));
// Assuming UTF-8 encoding on your terminal or whatever
	sprintf(temp,"%5.2f %c%c %5.2f m%c%c", rat, 0xC2,0xB1, err, 0xC2,0xBA);
	return temp;
}
char* stringEffErrMicro(int num, int den) {
	// Yea, its a memory leak.  Fight me.
	char* temp = new char[32];
	double rat = 1000000.0*double(num)/double(den);
	double err = 1000000.0*binerr(double(num),double(den));
// Assuming UTF-8 encoding on your terminal or whatever
	sprintf(temp,"%5.2f %c%c %5.2f %c%c%c%c", rat, 0xC2,0xB1, err, 0xCE,0xBC, 0xC2,0xBA);
	return temp;
}



// Yea does that file even exist?  Probably doesn't work for xroot access.
inline bool filehamna(const std::string& filename) {
    struct stat buf;
    int retval = stat(filename.c_str(), &buf);
    if (retval==0 && !S_ISREG(buf.st_mode)) {
        cout << "filehamna(string) : Not a regular file" << endl;
    }
    // -1 is the error condition, meaning that the file isn't there
    return (retval == -1);
}



#endif
