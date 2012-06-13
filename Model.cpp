#include "Model.h"


Model::Model(const double glucoseBoundaryConcentration, const double oxygenBoundaryConcentration,
	const ChemicalConsumptionRates &metabolicThresholds, const ChemicalConsumptionRates &scaledConsumptionRates, 
	const int mitoticCycleTime, const PorosityType porosityEqu,
	const int modelPrecisionDP, const double latticeSiteLength, 
	const double D0g, const double D0o, const double ifp, const int initialTumourRadius, 
	const double hoursPerCellUpdate, const double maxSimDays, ofstream &logFile): 
	m_tumourGrowthSites(NULL), m_tumourGrowthEdgeSites(NULL),
	m_glucoseBoundaryConcentration(glucoseBoundaryConcentration), 
	m_oxygenBoundaryConcentration(oxygenBoundaryConcentration), 
	m_metabolicThresholds(metabolicThresholds), m_scaledConsumptionRates(scaledConsumptionRates),
	m_mitoticCycleTime(mitoticCycleTime), m_porosityEqu(porosityEqu), 
	m_modelPrecisionDP(modelPrecisionDP), m_latticeSiteLength(latticeSiteLength),
	m_D0g(D0g), m_D0o(D0o), m_ifp(ifp), m_initialTumourRadius(initialTumourRadius), 
	m_hoursPerCellUpdate(hoursPerCellUpdate), m_maxSimDays(maxSimDays) {
	
	m_logFile = &logFile;
		
	for (int y = 0; y < HEIGHT; ++y) {
		for (int x = 0; x < WIDTH; ++x) {
			m_latticeT[y][x] = new LatticeSite(m_porosityEqu, m_D0g, m_D0o, m_ifp, 
				m_glucoseBoundaryConcentration, m_oxygenBoundaryConcentration, 
				m_metabolicThresholds, m_scaledConsumptionRates, m_mitoticCycleTime);
			m_latticeT1[y][x] = new LatticeSite(m_porosityEqu, m_D0g, m_D0o, m_ifp,
				m_glucoseBoundaryConcentration, m_oxygenBoundaryConcentration,
				m_metabolicThresholds, m_scaledConsumptionRates, m_mitoticCycleTime);
		}
	}

	// Added by Marjan 23/3/12
	for (int y = 0; y < HEIGHT * 2; ++y) {
		for (int x = 0; x < WIDTH * 2; ++x) {
			m_doubleLattice[y][x] =  m_glucoseBoundaryConcentration;
			m_doubleLattice[y][x] =  m_oxygenBoundaryConcentration;
		}
	}
	// Added by Marjan 23/3/12

	m_tumourGrowthEdgeSites = new vector<Coordinate>();
	m_tumourGrowthSites = new vector<Coordinate>();
	//Added by Marjan 8/3/12
	m_nonTumourGrowthRegion = new vector<Coordinate>();
	// Log Model parameters
	LogModelParameters();
}

Model::~Model()
{
	for (int y = 0; y < HEIGHT; ++y) {
		for (int x = 0; x < WIDTH; ++x) {
			delete m_latticeT[y][x];
			delete m_latticeT1[y][x];
		}
	}
	delete m_tumourGrowthEdgeSites;	
	delete m_tumourGrowthSites;
}

/*
 * Log model parameters
 */
void Model::LogModelParameters() {
	*m_logFile << "Model parameter: " << endl;
	*m_logFile << "   Maximum simulation days: " << m_maxSimDays << " [d]" << endl;
	*m_logFile << "   Cellular update (time) interval: " << m_hoursPerCellUpdate << " [h]" << endl; 
	*m_logFile << "Lattice dimensions: " << endl;
	*m_logFile << "   Width: " << WIDTH << " [100 um]" << endl; 
	*m_logFile << "   Height: " << HEIGHT << " [100 um]" << endl;
	*m_logFile << "   Lattice site length: " << m_latticeSiteLength << " [100 um]" << endl;
	*m_logFile << "Chemical parameters: " << endl;
	*m_logFile << "   Glucose boundary concentration: " << m_glucoseBoundaryConcentration << " [(mol)/(100 um)^3]" << endl;
	*m_logFile << "   Oxygen boundary concentration: " << m_oxygenBoundaryConcentration << " [(mol)/(100 um)^3]" << endl; 
	*m_logFile << "   Chemical consumption rates: " << endl;
	*m_logFile << "      Glucose metabolic thresholds: " << endl;
	*m_logFile << "         pat: " << m_metabolicThresholds.patGlucose << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         pant: " << m_metabolicThresholds.pantGlucose << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         qat: " << m_metabolicThresholds.qatGlucose << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         qant: " << m_metabolicThresholds.qantGlucose << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         dead: " << m_metabolicThresholds.dGlucose << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "      Oxygen metabolic thresholds: " << endl;
	*m_logFile << "         pat: " << m_metabolicThresholds.patOxygen << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         pant: " << m_metabolicThresholds.pantOxygen << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         qat: " << m_metabolicThresholds.qatOxygen << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         qant: " << m_metabolicThresholds.qantOxygen << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "         dead: " << m_metabolicThresholds.dOxygen << " [(mol)/(cell)(s)]" << endl;
	*m_logFile << "      Scaled glucose: " << endl;
	*m_logFile << "         pat: " << m_scaledConsumptionRates.patGlucose << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         pant: " << m_scaledConsumptionRates.pantGlucose << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         qat: " << m_scaledConsumptionRates.qatGlucose << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         qant: " << m_scaledConsumptionRates.qantGlucose << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         dead: " << m_scaledConsumptionRates.dGlucose << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "      Scaled oxygen: " << endl;
	*m_logFile << "         pat: " << m_scaledConsumptionRates.patOxygen << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         pant: " << m_scaledConsumptionRates.pantOxygen << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         qat: " << m_scaledConsumptionRates.qatOxygen << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         qant: " << m_scaledConsumptionRates.qantOxygen << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "         dead: " << m_scaledConsumptionRates.dOxygen << " [(100 um)^3/(cell)(s)]" << endl;
	*m_logFile << "Tumour parameters: " << endl;
	*m_logFile << "   Initial tumour radius: " << m_initialTumourRadius << " [100 um]" << endl;
	*m_logFile << "   Mitotic cycle time: " << m_mitoticCycleTime << " [h]" << endl;
	*m_logFile << "   Porosity equation: " << m_porosityEqu << endl;
	*m_logFile << "ECM parameters: " << endl;
	*m_logFile << "   IFP: " << m_ifp << " [kg/(100 um)(s)^2]" << endl;
	*m_logFile << "   Glucose diffusion coefficient: " << m_D0g << " [(100 um)^2/s]" << endl;
	*m_logFile << "   Oxygen diffusion coefficient: " << m_D0o << " [(100 um)^2/s]" << endl;
	*m_logFile << endl;
	
}
