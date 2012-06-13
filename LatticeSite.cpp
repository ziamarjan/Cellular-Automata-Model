#include "LatticeSite.h"
#include <cstdlib>
#include <sstream> 

/*
 * Default constructor
 */
 LatticeSite::LatticeSite(const PorosityType porosityEqu, const double D0g, const double D0o, const double ifp,
 	const double ECMGlucoseConcentration, const double ECMOxygenConcentration,
 	const ChemicalConsumptionRates &metabolicThresholds, 
 	const ChemicalConsumptionRates &scaledConsumptionRates,
	const int mitoticCycleTime):
	// Initialise the chemical concentrations 
	m_glucoseConcentration(ECMGlucoseConcentration), 
	m_oxygenConcentration(ECMOxygenConcentration),
	// Initialise the chemical diffusion coefficients 
	m_D0g(D0g), m_D0o(D0o), m_Dg(D0g), m_Do(D0o),
	// Initialise the metabolic thresholds
	m_metabolicThresholds(metabolicThresholds),
	// Initialise the chemical consumption rates
	m_scaledConsumptionRates(scaledConsumptionRates),	
	// Initialise the cell doubling time
	m_mitoticCycleTime(mitoticCycleTime),
	// Initialise the physical properties
	m_ifp(ifp), m_porosity(1.0), m_alpha(IntrinsicCellVolume() * ifp), m_pressure(ifp),
	m_porosityEqu(porosityEqu),
	// Initialise the cell collections
	m_patCells(0), m_pantCells(0), m_qatCells(0), m_qantCells(0), m_dCells(0), m_boneCells(0) {
}
 
/*
 * Copy constructor
 */
LatticeSite::LatticeSite(const LatticeSite &latticeSite): 
	// Initialise the chemical concentrations 
	m_glucoseConcentration(latticeSite.m_glucoseConcentration), 
	m_oxygenConcentration(latticeSite.m_oxygenConcentration),
	// Initialise the chemical diffusion coefficients 
	m_D0g(latticeSite.m_D0g), m_D0o(latticeSite.m_D0o), 
	m_Dg(latticeSite.m_Dg), m_Do(latticeSite.m_Do), 
	// Initialise the metabolic thresholds
	m_metabolicThresholds(latticeSite.m_metabolicThresholds),
	// Initialise the chemical consumption rates
	m_scaledConsumptionRates(latticeSite.m_scaledConsumptionRates),	
	// Initialise the physical properties
	m_ifp(latticeSite.m_ifp), m_porosity(latticeSite.m_porosity), 
	m_alpha(IntrinsicCellVolume() * latticeSite.m_ifp), 
	m_pressure(latticeSite.m_pressure), m_porosityEqu(latticeSite.m_porosityEqu),
	// Initialise the cell densities
	m_patCells(latticeSite.m_patCells), m_pantCells(latticeSite.m_pantCells),
	m_qatCells(latticeSite.m_qatCells), m_qantCells(latticeSite.m_qantCells),
	m_dCells(latticeSite.m_dCells), m_boneCells(latticeSite.m_boneCells),
	// Initialise the mitotic cycle time
	m_mitoticCycleTime(latticeSite.m_mitoticCycleTime) {
}

/*
 * Destructor
 */
LatticeSite::~LatticeSite() {
}

/*
 * Update the cell's metabolic state and proliferate new cells
 */
void LatticeSite::CellularUpdate(const double hoursPerCellUpdate) {
	
	// Find the total number of alive cells in this site
	double aliveCells = m_patCells + m_pantCells + m_qatCells + m_qantCells;
	// Re-set the cells in preparation of a state transition
	m_patCells = m_pantCells = m_qatCells = m_qantCells = 0.0;
	
	if (aliveCells > 0) {
	
		// Is there sufficient oxygen for this cell to transition to a PAT state
		if ( m_oxygenConcentration >= m_metabolicThresholds.patOxygen) {
			// Is there sufficient glucose for this cell to transition to a PAT state
			if (m_glucoseConcentration >= m_metabolicThresholds.patGlucose) {
				// PAT state
				m_patCells = aliveCells;
			}
			else {
				// QAT state
				m_qatCells = aliveCells;
			}
		}
		// Is there sufficient glucose for this cell to transition to a PANT state
		else if (m_glucoseConcentration >= m_metabolicThresholds.pantGlucose) {
			// PANT state
			m_pantCells = aliveCells;
		}
		// Is there sufficient glucose for this cell to transition to a QANT state
		else if (m_glucoseConcentration >= m_metabolicThresholds.qantGlucose) {
			// QANT state
			m_qantCells = aliveCells;
		}
		else {
			// DEAD state
			m_dCells += aliveCells;
		}
	}

	// Find this site's maximum cellular density
	int maxDensity = 0;
	switch (m_porosityEqu) {
		case HYP_POROSITY_106701:
			maxDensity = 106701;
			break;
		case HYP_POROSITY_1224:
			maxDensity = 1224;
			break;
		case HYP_POROSITY_691:
			maxDensity = 691;
			break;
		case HYP_POROSITY_300:
			maxDensity = 300;
			break;
		default:
			throw string("Error: porosity equation not recognised");
	}
	
	//if (m_porosity > MIN_POROSITY) {
	if (GetCellularDensity() < maxDensity) {
		if (m_mitoticCycleTime <= 0) {
			stringstream ss;
			ss << "Error: the mitotic cycle time is invalid." << endl;
			ss << "   Mitotic cycle time: " << m_mitoticCycleTime << endl; 
			throw ss.str(); 
		} 
		if (hoursPerCellUpdate <= 0.0) {
			stringstream ss;
			ss << "Error: the hours per cell update is invalid." << endl;
			ss << "   Hours per cell update: " << hoursPerCellUpdate << endl;
		}		
		
		// Proliferate new cells
		if (m_patCells > 0) {
			m_patCells += hoursPerCellUpdate * (log(2.0)/m_mitoticCycleTime) * m_patCells;
		}
		if (m_pantCells > 0) {
			m_pantCells += hoursPerCellUpdate * (log(2.0)/m_mitoticCycleTime) * m_pantCells;
		}
	// The cellular density has changed - update the porosity
	UpdatePorosity();	
	}
}

/*
 * The cells at this site consume glucose
 */
void LatticeSite::ConsumeGlucose() {
	double gConsumption = GetTotalGlucoseConsumption();
	if (m_glucoseConcentration - gConsumption  >= 0.0) 
		m_glucoseConcentration -= gConsumption;
	else 
		m_glucoseConcentration = 0.0;
}

/*
 * The cells at this site consume oxygen
 */
void LatticeSite::ConsumeOxygen() {
	double oConsumption = GetTotalOxygenConsumption();
	if (m_oxygenConcentration - oConsumption >= 0.0)
		m_oxygenConcentration -= oConsumption;
	else 
		m_oxygenConcentration = 0.0;
}

/*
 * Calculate the cell volume from the simple cell volume pressure relationship 
 */
double LatticeSite::SimpleCellVolume(double cellPressure) {
	return (m_alpha / cellPressure);
}

/*
 * Updates this site's porosity using a specified equation, called after this site's cellular composition
 * has changed.
 */
void LatticeSite::UpdatePorosity() {
		
	if (GetCellularDensity() < CriticalCellularDensity()) {
		m_porosity = 1.0 - (IntrinsicCellVolume() * GetCellularDensity());
	}
	// Calculate the porosity using the hypothetical porosity equation
	else if (m_porosityEqu == HYP_POROSITY_106701) {
		m_porosity = pow(1.001, (1 - (IntrinsicCellVolume() * GetCellularDensity()))) 
			- (pow(1.001, (1 - (IntrinsicCellVolume() * CriticalCellularDensity()))) - (1.0 - MAX_PACKING_DENSITY));   
	}
	else if (m_porosityEqu == HYP_POROSITY_691) {
		m_porosity = pow(1.2, (1 - (IntrinsicCellVolume() * GetCellularDensity())))
			- (pow(1.2, (1 - (IntrinsicCellVolume() * CriticalCellularDensity()))) - (1.0 - MAX_PACKING_DENSITY));
	}
	else if (m_porosityEqu == HYP_POROSITY_1224) {
		m_porosity = pow(1.1, (1 - (IntrinsicCellVolume() * GetCellularDensity())))
			- (pow(1.1, (1 - (IntrinsicCellVolume() * CriticalCellularDensity()))) - (1.0 - MAX_PACKING_DENSITY));
	}
	else if (m_porosityEqu == HYP_POROSITY_300) {
		m_porosity = pow(1.75951, (1 - (IntrinsicCellVolume() * GetCellularDensity())))
			- (pow(1.75951, (1 - (IntrinsicCellVolume() * CriticalCellularDensity()))) - (1.0 - MAX_PACKING_DENSITY));		
	}
	else {
		stringstream ss;
		ss << "Error: unrecognised porosity equation." << endl;
		ss << "   Porosity equation: " << m_porosityEqu << endl;
		throw ss.str();
	}
	
	//if (m_porosity < MIN_POROSITY)
	//	m_porosity = MIN_POROSITY;
	if (m_porosity < 0)
		m_porosity = 0;
	// Update this site's diffusion coefficients
	UpdateDiffusionCoefficients();
	// Update this site's pressure
	UpdatePressure();
}

/*
 * Updates this site's diffusion coefficients, called after the porosity has been updated
 */
void LatticeSite::UpdateDiffusionCoefficients() {
	double tortuosity = sqrt(m_porosity);
	m_Dg = tortuosity * m_D0g;
	m_Do = tortuosity * m_D0o;
}

/*
 * Update this site's pressure, called after the porosity has been updated
 */
void LatticeSite::UpdatePressure() {
	if (m_alpha <= 0.0) {
		stringstream ss;
		ss << "Error: alpha has not been initialised." << endl;
		ss << "   Alpha: " << m_alpha << endl;
		throw ss.str(); 
	}
	if (m_ifp <= 0.0) {
		stringstream ss;
		ss << "Error: the IFP has not been initialised." << endl;
		ss << "   IFP: " << m_ifp << endl;
		throw ss.str();
	}	 
	
	double density = GetCellularDensity();
	if (density < CriticalCellularDensity()) 
		m_pressure = m_ifp;
	else if (density >= CriticalCellularDensity()) 
		m_pressure = (m_alpha * density)/(1 - m_porosity);
	else {
		stringstream ss;
		ss << "Error: negative cellular density." << endl;
		ss << "   Cellular density: " << density << endl;
		throw ss.str();
	}	
}

/*******************************************************************************************
 * Mutators
 * *****************************************************************************************/

/*
 * Set this site's glucose concentration
 */
void LatticeSite::SetGlucoseConcentration(const double glucoseConcentration) {
	m_glucoseConcentration = glucoseConcentration;
}

/*
 * Set this site's oxygen concentration
 */
void LatticeSite::SetOxygenConcentration(const double oxygenConcentration) {
	m_oxygenConcentration = oxygenConcentration;
}

/*
 * Set this site's initial cellular composition
 */
void LatticeSite::SetCells(const int patCellPopulation, const int pantCellPopulation,
	const int qatCellPopulation, const int qantCellPopulation, const int dCellPopulation,
	const int boneCellPopulation) {		
	
	// Create this site's PAT population
	m_patCells = patCellPopulation;
	
	// Create this site's PANT population
	m_pantCells = pantCellPopulation;
	
	// Create this site's QAT population
	m_qatCells = qatCellPopulation;
	
	// Create this site's QANT population
	m_qantCells = qantCellPopulation;
	
	// Create this site's D population
	m_dCells = dCellPopulation;

	// Create this site's bone population
	m_boneCells = boneCellPopulation;
	
	// Update this site's porosity and diffusion ceofficients
	UpdatePorosity();	
}

/*
 * Set the number of PAT cells after migration
 */
void LatticeSite::SetPATCellPopulation(const double patCells) {
	m_patCells = patCells;
	// Update this site's physical properties
	UpdatePorosity();
}

/*
 * Set the number of PANT cells after migration
 */ 
void LatticeSite::SetPANTCellPopulation(const double pantCells) {
	m_pantCells = pantCells;
	UpdatePorosity();
}

/*
 * Set the number of QAT cells after migration
 */
void LatticeSite::SetQATCellPopulation(const double qatCells) {
	m_qatCells = qatCells;
	UpdatePorosity();
}

/*
 * Set the number of QANT cells after migration
 */
void LatticeSite::SetQANTCellPopulation(const double qantCells) {
	m_qantCells = qantCells;
	UpdatePorosity();
}

/*
 * Set the number of dead cells after migration
 */
void LatticeSite::SetDCellPopulation(const double dCells) {
	m_dCells = dCells;
	// Update this site's physical properties
	UpdatePorosity();
}

/*
 * Set the number of bone cells after migration
 */
void LatticeSite::SetBoneCellPopulation(const double boneCells) {
	m_boneCells = boneCells;
	UpdatePorosity();
}

/*******************************************************************************************
 * Accessors
 * *****************************************************************************************/
 
/*
 * Return the total cellular density
 */
double LatticeSite::GetCellularDensity() const {
	return GetPATCellPopulation() + GetPANTCellPopulation() + GetQATCellPopulation() 
		+ GetQANTCellPopulation() + GetDCellPopulation() + GetBoneCellPopulation();
}

/*
 * Return the glucose consumption at this site
 */
double LatticeSite::GetTotalGlucoseConsumption() const {
	double glucoseConsumption = 
		(GetPATCellPopulation() * m_scaledConsumptionRates.patGlucose * m_glucoseConcentration) +
		(GetPANTCellPopulation() * m_scaledConsumptionRates.pantGlucose * m_glucoseConcentration) +
		(GetQATCellPopulation() * m_scaledConsumptionRates.qatGlucose * m_glucoseConcentration) +
		(GetQANTCellPopulation() * m_scaledConsumptionRates.qantGlucose * m_glucoseConcentration) +
		(GetDCellPopulation() * m_scaledConsumptionRates.dGlucose * m_glucoseConcentration);
	return glucoseConsumption;
}

/*
 * Return the oxygen consumption at this site
 */ 
double LatticeSite::GetTotalOxygenConsumption() const {
	double oxygenConsumption = 
		(GetPATCellPopulation() * m_scaledConsumptionRates.patOxygen * m_oxygenConcentration) +
		(GetPANTCellPopulation() * m_scaledConsumptionRates.pantOxygen * m_oxygenConcentration) +
		(GetQATCellPopulation() * m_scaledConsumptionRates.qatOxygen * m_oxygenConcentration) +
		(GetQANTCellPopulation() * m_scaledConsumptionRates.qantOxygen * m_oxygenConcentration) +
		(GetDCellPopulation() * m_scaledConsumptionRates.dOxygen * m_oxygenConcentration);
	return oxygenConsumption;
}
