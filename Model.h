#ifndef MODEL_H_
#define MODEL_H_

#include "Coordinate.h"
#include "LatticeSite.h"

//#include "CellPopulation.h" //Commented out by Marjan 26/4/12

#include <vector>
#include <fstream> 

using namespace std;

#define HEIGHT 121
#define WIDTH 121

enum LatticeType {T, T1};

class Model
{
public:
	Model(const double glucoseBoundaryConcentration, const double oxygenBoundaryConcentration,
		const ChemicalConsumptionRates &metabolicThresholds, const ChemicalConsumptionRates &scaledConsumptionRates,
		const int mitoticCycleTime, 
		const PorosityType porosityEqu, const int modelPrecisionDP, const double latticeSiteLength,
		const double D0g, const double D0o, const double ifp, const int initialTumourRadius, 
		const double hoursPerCellUpdate, const double maxSimDays, ofstream &logFile);

	virtual ~Model();
	
	// This routine can be substituted for others
	static double ScaledChemicalConsumptionRate(const double chemicalConsumptionRate,
		const double chemicalBoundaryConcentration) {
			return chemicalConsumptionRate/chemicalBoundaryConcentration; }
	
	// Rounds a value to a certain precions for comparing doubles		
	double RoundToModelPrecision(const double value) const {
		double val1 = value * pow(10.0, m_modelPrecisionDP);
		val1 += 0.5;
		int val2 = (int)val1;
		return ((double)val2)/pow(10.0, m_modelPrecisionDP);
	}
	
	void LogModelParameters();
	
	// The execution file handle for logging model updates
	ofstream *m_logFile;
	
	// The lattice sites at the edge of the tumour growth region
	vector<Coordinate> *m_tumourGrowthEdgeSites;
	// The lattice sites making up the tumour growth region
	vector<Coordinate> *m_tumourGrowthSites;
	// Added by Marjan
	vector<Coordinate> *m_nonTumourGrowthRegion;
	
	// The lattice
	LatticeSite *m_latticeT[HEIGHT][WIDTH];
	LatticeSite *m_latticeT1[HEIGHT][WIDTH];
	double m_doubleLattice[2*HEIGHT][2*WIDTH];
	
	// Model precision for comparing doubles
	const int m_modelPrecisionDP;
	
	// A lattice site's length [100 um]
	const double m_latticeSiteLength;
	
	// The chemical boundary conditions
	const double m_glucoseBoundaryConcentration;
	const double m_oxygenBoundaryConcentration;
	// Cell chemical consumption rates, which are universal to all tumour sites
	const ChemicalConsumptionRates m_metabolicThresholds;
	// Scaled chemical consumption rates
	const ChemicalConsumptionRates m_scaledConsumptionRates;
	
	// Cell doubling time, which are universal to all tumour cells
	const int m_mitoticCycleTime;
	
	// ECM properties
	// Glucose diffusion coefficient [(100 um)^2/s]
	const double m_D0g;
	// Oxygen diffusion coefficient [(100 um)^2/s]  
	const double m_D0o;
	// Interstitium Fluid Pressure [kg/(100 um)(s)^2]
	const double m_ifp;
	
	// The equation for calculating a site's porosity
	const PorosityType m_porosityEqu;
	
	// The initial tumour radius [100 um]
	const int m_initialTumourRadius;
	
	// The cellular update interval [h]
	const double m_hoursPerCellUpdate;
	// The maximum number of simulation days [d]
	const double m_maxSimDays;
};

#endif /*MODEL_H_*/
