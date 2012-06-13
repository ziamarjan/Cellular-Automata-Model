#ifndef LATTICESITE_H_
#define LATTICESITE_H_

#include <string>
#include <cmath>
#include <list> 
#include "CellPopulation.h"


using namespace std;


#define PI 3.14159

// The maximum irregular packing density of spheres
#define MAX_PACKING_DENSITY 0.64

// The minimum porosity that corresponds to a porosity of zero
#define MIN_POROSITY 1e-4

// The radius of a cell under normal conditions (IFP) [100 um]
const double intrinsicCellRadius = 0.1;

// Calculates the volume of a cell under normal conditions (IFP) [(100 um)^3]
inline double IntrinsicCellVolume() { return ((4 * PI * pow(intrinsicCellRadius, 3)) / 3); }

// Calculates the density at which all the cells come into contact [cell/(100 um)^3]
inline int CriticalCellularDensity() { return ceil(MAX_PACKING_DENSITY / IntrinsicCellVolume()); }

// A type for the different methods of calculating a site's porosity. 
// The trailing number indicates the cellular density at which the porosity becomes zero.
enum PorosityType { HYP_POROSITY_106701, HYP_POROSITY_691, HYP_POROSITY_1224, HYP_POROSITY_300 };

// Cell chemical consumption rates
struct ChemicalConsumptionRates {
	double patGlucose;
	double pantGlucose;
	double qatGlucose;
	double qantGlucose;
	double dGlucose;
	double patOxygen;
	double pantOxygen;
	double qatOxygen;
	double qantOxygen;
	double dOxygen;
};

/*
 * A structure for storing lattice site attributes
 */
class LatticeSite
{
public:
	LatticeSite(const PorosityType porosityEqu, const double D0g, const double D0o, const double ifp,
 	const double ECMGlucoseConcentration, const double ECMOxygenConcentration,
 	const ChemicalConsumptionRates &metabolicThresholds, 
 	const ChemicalConsumptionRates &scaledConsumptionRates,
	const int mitoticCycleTime);
	LatticeSite(const LatticeSite &latticeSite);
	virtual ~LatticeSite();
	
	// Perform a cellular update
	void CellularUpdate(const double hoursPerCellUpdate);
	
	// Consume chemicals
	void ConsumeGlucose();
	void ConsumeOxygen();
	
	// Chemical concentrations
	void SetGlucoseConcentration(const double glucoseConcentration); 	
	void SetOxygenConcentration(const double oxygenConcentration);

	// Set chemical diffusion coefficients
	void SetApparentGlucoseDiffusionCoefficient(const double Dg) { m_Dg = Dg; }
	void SetApparentOxygenDiffusionCoefficient(const double Do) { m_Do = Do; }
	
	void SetCells(const int patCellPopulation, const int pantCellPopulation,
		const int qatCellPopulation, const int qantCellPopulation, const int dCellPopulation,
		const int boneCellPopulation);

		
	void SetPATCellPopulation(const double patCells);
	void SetPANTCellPopulation(const double pantCells);
	void SetQATCellPopulation(const double qatCells);
	void SetQANTCellPopulation(const double qantCells);
	void SetDCellPopulation(const double dCells);
	void SetBoneCellPopulation(const double boneCells);
	
	// Chemical diffusion coefficients
	double GetApparentGlucoseDiffusionCoefficient() const { return m_Dg; }
	double GetApparentOxygenDiffusionCoefficient() const { return m_Do; }
	// Chemical concentrations
	double GetGlucoseConcentration() const { return m_glucoseConcentration; }
	double GetOxygenConcentration() const { return m_oxygenConcentration; }
	
	// Chemical consumption 
	double GetTotalGlucoseConsumption() const;
	double GetTotalOxygenConsumption() const;
	
	// Scaled chemical consumption rates
	double GetScaledPATGlucoseConsumptionRate() const { return m_scaledConsumptionRates.patGlucose; }
	double GetScaledPANTGlucoseConsumptionRate() const { return m_scaledConsumptionRates.pantGlucose; }
	double GetScaledQATGlucoseConsumptionRate() const { return m_scaledConsumptionRates.qatGlucose; }
	double GetScaledQANTGlucoseConsumptionRate() const { return m_scaledConsumptionRates.qantGlucose; }
	double GetScaledDGlucoseConsumptionRate() const { return m_scaledConsumptionRates.dGlucose; }
	
	double GetScaledPATOxygenConsumptionRate() const { return m_scaledConsumptionRates.patOxygen; }
	double GetScaledPANTOxygenConsumptionRate() const { return m_scaledConsumptionRates.pantOxygen; }
	double GetScaledQATOxygenConsumptionRate() const { return m_scaledConsumptionRates.qatOxygen; }
	double GetScaledQANTOxygenConsumptionRate() const { return m_scaledConsumptionRates.qantOxygen; }
	double GetScaledDOxygenConsumptionRate() const { return m_scaledConsumptionRates.dOxygen; }
	
	// Cell numbers
	double GetPATCellPopulation() const { return m_patCells; }
	double GetPANTCellPopulation() const { return m_pantCells; }
	double GetQATCellPopulation() const { return m_qatCells; }
	double GetQANTCellPopulation() const { return m_qantCells; }
	double GetDCellPopulation() const { return m_dCells; }
	double GetBoneCellPopulation() const { return m_boneCells; }
	double GetCellularDensity() const;	 
	
	double GetPressure() const { return m_pressure; }
	

private:
	// Calculate the cell volume from the simple cell volume pressure relationship
	double SimpleCellVolume(double cellPressure);
	
	void UpdatePorosity();
	void UpdateDiffusionCoefficients();
	void UpdatePressure();	
	
	// Cell collections
	double m_patCells;
	double m_pantCells;
	double m_qatCells;
	double m_qantCells;
	double m_dCells;
	double m_boneCells;
	
	// The cell doubling time
	const int m_mitoticCycleTime;
	
	// Chemical attributes
	// Default: 2.8 * 10 ^ -12 [mol / (100 um)^3]
	double m_glucoseConcentration;
	// Default: 3.5 * 10 ^ -14 [mol / (100 um)^3]
	double m_oxygenConcentration;
	// The glucose diffusion coefficient of the extra-cellular matrix, default: 0.91 [(100 um)^2 / s]
	const double m_D0g;
	// The oxygen diffusion coefficient of the extra-cellular matrix, default: 0.182 [(100 um)^2 / s]
	const double m_D0o;
	// The apparent glucose diffusion ceofficient
	double m_Dg;
	// The apparent oxygen diffusion ceofficient
	double m_Do;
	
	// Average glucose consumption rate [(mol)/(cell)(s)]
	const ChemicalConsumptionRates m_metabolicThresholds;
	// Scaled glucose consumption rate [(100 um)^3/(cell)(s)]
	const ChemicalConsumptionRates m_scaledConsumptionRates;
		
	// Extra-cellular matrix pressure
	// Default: Interstitium Fluid Pressure (IFP) 0.0933257 [kg/(100 um)(s)^2]
	const double m_ifp;
	// The equation for calculating porosity
	const PorosityType m_porosityEqu;
	// Cellular Porosity 
	double m_porosity;
	// Simple cell volume pressure relationship constant
	const double m_alpha;
	// Site pressure [kg/(100 um)(s)^2]
	double m_pressure;
};

#endif /*LATTICESITE_H_*/
