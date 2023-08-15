// A selection of codebits for analysis related to particles & their PID PDG codes
// Eventually, you should switch over to ROOT's tool for this:
//
//			TDatabasePDG* pdgInstance = TDatabasePDG::Instance();
//			TParticlePDG* particle = pdgInstance->GetParticle(2212);
//			particle->Mass();		// returns 0.938272
//			particle->Charge();		// returns 3.000000 ('cuz quarks)
//			particle->GetTitle();	// returns "proton"
//



#ifndef PDGutils_h
#define PDGutils_h

using namespace std;


int const nothingPDG   =    0;
int const protonPDG    = 2212;
int const neutronPDG   = 2112;
int const electronPDG  =   11;	// Means an electron not a positron
int const muonPDG      =   13;
int const tauPDG       =   15;
int const nuePDG       =   12;
int const numuPDG      =   14;
int const nutauPDG     =   16;
int const pichPDG      =  211;	// Means a pi+; pi- is -211
int const pi0PDG       =  111;
int const etaPDG       =  221;
int const omegaPDG     =  223;
int const etaPRPDG     =  331;
int const rho0PDG      =  113;
int const rhoCHPDG     =  213;
int const KchPDG       =  321;	// Means a K+; K- is -321
int const KlongPDG     =  130;
int const KshortPDG    =  310;
int const KzeroPDG     =  311;
int const gammaPDG     =   22;

int const ArgonPDG     = 1000180400;
int const CarbonPDG    = 1000060120;
int const DeuteronPDG  = 1000010020;
int const HeliumPDG    = 1000020040;
int const HydrogenPDG  = 1000010010;
int const OxygenPDG    = 1000080160;

int const lambdaPDG    = 3122;
int const sigmaM_PDG   = 3112;
int const sigma0_PDG   = 3212;
int const sigmaP_PDG   = 3222;
int const cascadeM_PDG = 3312;
int const cascade0_PDG = 3322;
int const OmegaM_PDG   = 3334;


double const Melectron = 0.00051100;	// in GeV
double const Mmuon     = 0.10565837;
double const Mtau      = 1.77686;
double const Mpi0      = 0.13497770;
double const Meta      = 0.547862;
double const Momega    = 0.78265;
double const MetaPR    = 0.95778;
double const Mrho      = 0.77526;
double const MpiCH     = 0.13957061;
double const MchK      = 0.493677;
double const MneuK     = 0.497611;
double const Mprot     = 0.93827208;
double const Mneut     = 0.93956541;
double const Mlambda   = 1.115683;
double const MsigmaNeg = 1.197449;
double const Msigma0   = 1.192642;
double const MsigmaPos = 1.18937;
double const Mcasc0    = 1.31486;
double const McascNeg  = 1.32171;
double const MOmegaNeg = 1.67245;


vector<vector<int>> generateCombinations(int n) {
    vector<vector<int>> result;
    vector<int> combination(n, 0); // Initialize combination with all 0s

    // Helper function to recursively generate combinations
    std::function<void(int)> generate = [&](int index) {
        // Base case: If we have reached the end of the combination vector
        if (index == n) {
            // Add the generated combination to the result vector
            result.push_back(combination);
            return;
        }

        // Generate combinations for current position with 0
        combination[index] = 0;
        generate(index + 1);

        // Generate combinations for current position with 1
        combination[index] = 1;
        generate(index + 1);
    };

    // Start generating combinations from index 0
    generate(0);

    return result;
}

double calculateMean(const std::vector<TVector3>& vectors) {
    TVector3 sum(0.0, 0.0, 0.0);

    for (const TVector3& vector : vectors) {
        sum += vector;
    }

    return vectors.empty() ? 0.0 : sum.Mag() / vectors.size();
}

double calculateStandardDeviation(const std::vector<TVector3>& vectors, double mean) {
    double sumSquaredDifference = 0.0;

    for (const TVector3& vector : vectors) {
        double diff = (vector.Mag() - mean);
        sumSquaredDifference += diff * diff;
    }

    return vectors.empty() ? 0.0 : std::sqrt(sumSquaredDifference / vectors.size());
}
// Mass of this particle, in GeV, from PDG code.  Negative value is err
double massPDG(int PDGcode) {

	if (PDGcode > 1000000000) {
		int A = (PDGcode -1000000000)/10;
		A = A%1000;
		return A * (39.948 / 40.0);
	}

	switch (PDGcode){
		case +electronPDG:	case -electronPDG:
			return Melectron;

		case +muonPDG:		case -muonPDG:
			return Mmuon;

		case +tauPDG:		case -tauPDG:
			return Mtau;

		case  gammaPDG:
		case +nuePDG:		case -nuePDG:
		case +numuPDG:		case -numuPDG:
		case +nutauPDG:		case -nutauPDG:
			return 0.000000;

		case +pichPDG:		case -pichPDG:
			return Mpi0;

		case  rho0PDG:
		case +rhoCHPDG:		case -rhoCHPDG:
			return Mrho;

		case +KchPDG:		case -KchPDG:
			return MpiCH;

		case +protonPDG:	case -protonPDG:
			return Mprot;

		case  pi0PDG:
			return Mpi0;

		case  etaPDG:
			return Meta;

		case  omegaPDG:
			return Momega;

		case  etaPRPDG:
			return MetaPR;

		case  KlongPDG:
			return MneuK;

		case +neutronPDG:	case -neutronPDG:
			return Mneut;

		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
			return MneuK;

		case +lambdaPDG:	case -lambdaPDG:
			return Mlambda;

		case +sigmaM_PDG:	case -sigmaM_PDG:
			return MsigmaNeg;

		case +sigma0_PDG:	case -sigma0_PDG:
			return  Msigma0;

		case +sigmaP_PDG:	case -sigmaP_PDG:
			return MsigmaPos;

		case +cascade0_PDG:	case -cascade0_PDG:
			return Mcasc0;

		case +cascadeM_PDG:	case -cascadeM_PDG:
			return McascNeg;

		case +OmegaM_PDG:	case -OmegaM_PDG:
			return MOmegaNeg;

		default:
			cout << "PDGutils balk: " << PDGcode << endl;
			return -999;
	}
}



// Charge of this particle, from the PDG code
int chargePDG(int PDGcode) {

	// Atoms
	if (PDGcode > 1000000000) {
		int retval = PDGcode / 10000;
		retval = retval % 1000;
		return retval;
	}

	switch (PDGcode){
		case -electronPDG:
		case -muonPDG:
		case -tauPDG:
		case +pichPDG:
		case +rhoCHPDG:
		case +KchPDG:
		case +protonPDG:
		case -sigmaM_PDG:
		case +sigmaP_PDG:
		case -cascadeM_PDG:
		case +OmegaM_PDG:
			return +1;

		case +electronPDG:
		case +muonPDG:
		case +tauPDG:
		case -pichPDG:
		case -rhoCHPDG:
		case -KchPDG:
		case -protonPDG:
		case +sigmaM_PDG:
		case -sigmaP_PDG:
		case +cascadeM_PDG:
		case -OmegaM_PDG:
			return -1;

		case  gammaPDG:
		case +nuePDG:		case -nuePDG:
		case +numuPDG:		case -numuPDG:
		case +nutauPDG:		case -nutauPDG:
		case  pi0PDG:		case  rho0PDG:
		case  etaPDG:		case etaPRPDG:
		case  omegaPDG:
		case  KlongPDG:		case  KshortPDG:
		case +KzeroPDG:		case -KzeroPDG:
		case +neutronPDG:	case -neutronPDG:
		case +lambdaPDG:	case -lambdaPDG:
		case +sigma0_PDG:	case -sigma0_PDG:
		case +cascade0_PDG:	case -cascade0_PDG:
			return false;

		default:
			cout << "PDGutils balk: " << PDGcode << endl;
			return -999;
	}
};
bool chargedPDG(int PDGcode) {
	int charge = chargePDG(PDGcode);
	if (charge==0 || charge==-999) return false;
	return true;
}

Int_t GetZ_PDG(int PDG){

  Int_t count_digits = int(log10(PDG) + 1);
  if(count_digits!=10) {
	cout<< "Not an atomic PDGcode: "<<PDG<<endl;
	return -999;
  }
  
  Int_t Z = (PDG % 10000000) /10000;
  return Z;

}

Int_t GetA_PDG(int PDG){

  Int_t count_digits = int(log10(PDG) + 1);
  if(count_digits!=10) {
	cout<< "Not an atomic PDGcode: "<<PDG<<endl;
	return -999;
  }
  
  Int_t A = (PDG % 10000) /10;
  return A;

}


#endif
