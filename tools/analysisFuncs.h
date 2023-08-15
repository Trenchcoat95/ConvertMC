// A selection of routines to compute useful analysis quantities in ../CoherentAna/????.cxx etc.



#ifndef analysisFuncs_h
#define analysisFuncs_h



#include "PDGutils.h"



// In a fiducial corresponding to 1 metric ton of Ar nuclei?  (Corresponds to
// about 35cm from edge of the TPC).
bool inFiducial(double x, double y, double z) {
	double const Rcut = 445.0 /2.0;		double const Zcut = 430.0 /2.0;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};



// Is this particle likely to be seen in the TPC as from vertex?
// Means of course only stable charged particles
bool visibleTPC(int PDGcode, double E, double distFromVert) {

	double EthreshChgd = 0.003;		// In GeV
	double visibleDcut = 2.0;		// in cm

	if ( chargedPDG(PDGcode) && E<EthreshChgd ) return false;
	if ( distFromVert > visibleDcut ) return false;

	switch (PDGcode){
		case +electronPDG:	case -electronPDG:
		case +muonPDG:		case -muonPDG:
		case +pichPDG:		case -pichPDG:
		case +KchPDG:		case -KchPDG:
		case +protonPDG:	case -protonPDG:
			return true;

		case  KlongPDG:		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
		case +neutronPDG:	case -neutronPDG:
			return false;

		default:
			return false;			// Nuclei, esp Ar40, not considered visible
	}
};



//  Are we inside a (12-sided) figure?
bool inDodecagon(double X, double Y, double Apothem) {
	X = fabs(X);	Y = fabs(Y);
	double theta = atan(Y/X);
	if (theta<    M_PI/12.0 && X>Apothem) return false;
	if (theta>5.0*M_PI/12.0 && Y>Apothem) return false;
	double dot;
	if (theta>    M_PI/ 4.0) {
		dot = (X + sqrt(3)*Y) / 2.0;
		if (dot>Apothem) return false;
	} else {
		dot = (sqrt(3)*X + Y) / 2.0;
		if (dot>Apothem) return false;		
	}
	return true;
};
bool insideECAL(double X, double Y, double Z) {
	// Inferred from gdml files for 42l w/ SPY_v3.  Could stand a better check.
	double const insideApothem     = 278.01;
	double const outsideApothem    = 308.21;
	double const insideHalfLength  = 329.25;
	double const outsideHalfLength = 359.25;
	
	X = fabs(X);	Y = fabs(Y);	Z = fabs(Z);
	
	if ( outsideHalfLength<=X) return false;
	if ( insideHalfLength <=X && X<outsideHalfLength ) {
		// Could be inside an endcap
		return inDodecagon(Z,Y,outsideApothem);
	}
	return (inDodecagon(Z,Y,outsideApothem) && !inDodecagon(Z,Y,insideApothem));
};
bool beyondECAL(double X, double Y, double Z) {
	// Inferred from gdml files for 42l w/ SPY_v3.  Could stand a better check.
	double const outsideApothem    = 308.21;
	double const outsideHalfLength = 359.25;
	
	X = fabs(X);	Y = fabs(Y);	Z = fabs(Z);

	if ( outsideHalfLength <=X ) return true;
	return (!inDodecagon(Z,Y,outsideApothem));
};



//  Are we inside the (cylindrical) TPC volume?
bool insideTPC(double x, double y, double z) {
	double const Rcut = 246.6;		double const Zcut = 249.7;
	double r = sqrt( z*z + y*y );
	if ( r > Rcut ) return false;
	if (fabs(x) > Zcut ) return false;
	return true;
};



// Is this particle likely to be seen in the MPD as from vertex?
// Means of course only stable particles - e.g. for a pi0 will be false, 
// but will be true for the two photons.
bool visibleMPD(int PDGcode, double E, double distFromVert) {

	double EthreshChgd = 0.003;		// In GeV
	double visibleDcut = 2.0;		// in cm

	if ( chargedPDG(PDGcode) && E<EthreshChgd ) return false;
	if ( distFromVert > visibleDcut ) return false;

	switch (PDGcode){	
		case +electronPDG:	case -electronPDG:
		case +muonPDG:		case -muonPDG:
		case  gammaPDG:
		case +pichPDG:		case -pichPDG:
		case +KchPDG:		case -KchPDG:
		case +protonPDG:	case -protonPDG:
			return true;

		case  KlongPDG:		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
		case +neutronPDG:	case -neutronPDG:
			return false;

		default:
			return false;
	}
};



// Is this particle likely to be seen in the ND LAr as from vertex?
bool visibleLAr(int PDGcode, double E, double distFromVert) {

	double EthreshChgd = 0.030;		// In GeV
	double visibleDcut = 2.0;		// in cm

	if ( chargedPDG(PDGcode) && E<EthreshChgd ) return false;
	if ( distFromVert > visibleDcut ) return false;

	switch (PDGcode){
		case +electronPDG:	case -electronPDG:
		case +muonPDG:		case -muonPDG:
		case  gammaPDG:
		case +pichPDG:		case -pichPDG:
		case +KchPDG:		case -KchPDG:
		case +protonPDG:	case -protonPDG:
			return true;

		case  KlongPDG:		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
		case +neutronPDG:	case -neutronPDG:
			return false;

		default:
			return false;
	}
};



// To kill protons using dE/dx
bool NoTaProton(double I,double p) {
	// I is average ionization
	// Units are those of TPCCluster.Signal().
	// p is momentum in GeV
	if ( p > 1.400 ) return true;
	if ( p > 1.000 && p < 1.400 && I >  1.75e3 ) return false;
	if ( p > 0.500 && p < 1.000 && I >  3.25e3 - 1.5e3*p ) return false;
	if ( p > 0.000 && p < 0.500 && I >  8.00e3 -11.0e3*p ) return false;
	return true;
};



// String for interaction code.  Must be fixed if MCNeutrino.h in
// nusimdata/SimulationBase ever changes.  Which is unlikely.  Is
// somewhat redundant to gar::sim::TruthInteractionTypeName, but has
// shorter output and handles 1000 correctly, at least as of GENIE 2.12.10
std::string interactionName(int interactionType) {
	switch (interactionType){
	case   -1: return "KnownUnknownInteraction";

	// These are Mode rather than the InteractionType; unusual values omitted
	case    0: return "QE";
	case    1: return "RES";
	case    2: return "DIS";
	case    3: return "COH";
	case    5: return "ENU";

	case 1000: return "Mystery Resonance";
	case 1001: return "CCQE";
	case 1002: return "NCQE";
	case 1003: return "ResCCNuProtonPi+";
	case 1004: return "ResCCNuNeutronPi0";
	case 1005: return "ResCCNuNeutronPi+";
	case 1006: return "ResNCNuProtonPi0";
	case 1007: return "ResNCNuProtonPi+";
	case 1008: return "ResNCNuNeutronPi0";
	case 1009: return "ResNCNuNeutronPi-";
	case 1010: return "ResCCNuBarNeutronPi-";
	case 1011: return "ResCCNuBarProtonPi0";
	case 1012: return "ResCCNuBarProtonPi-";
	case 1013: return "ResNCNuBarProtonPi0";
	case 1014: return "ResNCNuBarProtonPi+";
	case 1015: return "ResNCNuBarNeutronPi0";
	case 1016: return "ResNCNuBarNeutronPi-";
	case 1017: return "ResCCNuDelta+Pi+";
	case 1021: return "ResCCNuDelta2+Pi-";
	case 1028: return "ResCCNuBarDelta0Pi-";
	case 1032: return "ResCCNuBarDelta-Pi+";
	case 1039: return "ResCCNuProtonRho+";
	case 1041: return "ResCCNuNeutronRho+";
	case 1046: return "ResCCNuBarNeutronRho-";
	case 1048: return "ResCCNuBarNeutronRho0";
	case 1053: return "ResCCNuSigma+K+";
	case 1055: return "ResCCNuSigma+K0";
	case 1060: return "ResCCNuBarSigma-K0";
	case 1062: return "ResCCNuBarSigma0K0";
	case 1067: return "ResCCNuProtonEta";
	case 1070: return "ResCCNuBarNeutronEta";
	case 1073: return "ResCCNuK+Lambda0";
	case 1076: return "ResCCNuBarK0Lambda0";
	case 1079: return "ResCCNuProtonPi+Pi-";
	case 1080: return "ResCCNuProtonPi0Pi0";
	case 1085: return "ResCCNuBarNeutronPi+Pi-";
	case 1086: return "ResCCNuBarNeutronPi0Pi0";
	case 1090: return "ResCCNuBarProtonPi0Pi0";
	case 1091: return "CCDIS";
	case 1092: return "NCDIS";
	case 1093: return "UnUsed1";
	case 1094: return "UnUsed2";
	case 1095: return "CCQEHyperon";
	case 1096: return "NCCOH";
	case 1097: return "CCCOH";
	case 1098: return "NuElectronElastic";
	case 1099: return "InverseMuDecay";
	}          return "Idunno";
};



// A simplified PDG code for plotting
int simpleCode(int PDGcode) {
	int thisPDG;
	switch (PDGcode) {
	case (+electronPDG):
		thisPDG = -1;	break;
	case (-electronPDG):
		thisPDG = +1;	break;
	case (-muonPDG)	:
		thisPDG = -2;	break;
	case (+muonPDG)	:
		thisPDG = +2;	break;
	case (+pichPDG)	:
		thisPDG = +3;	break;
	case (-pichPDG)	:
		thisPDG = -3;	break;
	case (+protonPDG):
		thisPDG = +4;	break;
	case (-protonPDG):
		thisPDG = -4;	break;
	case (gammaPDG):
		thisPDG = -5;	break;
	default:
		thisPDG = +5;
	}
	return thisPDG;
};



#endif
