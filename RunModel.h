#ifndef RUNMODEL_H_
#define RUNMODEL_H_

#include "Model.h"
#include <string>
#include <vector>

using namespace std;

enum DisplayType {GLUCOSE_CONCENTRATION, OXYGEN_CONCENTRATION, TGFBETA_CONCENTRATION,GLUCOSE_DIFFCOEF, OXYGEN_DIFFCOEF,
	GLUCOSE_CONSUMPTION, OXYGEN_CONSUMPTION,
	// Cell population by type 
	PAT_POPULATION, PANT_POPULATION, QAT_POPULATION, QANT_POPULATION, D_POPULATION, CELL_POPULATION, BONE_POPULATION};
enum ChemicalType {GLUCOSE, OXYGEN, TGFBETA};

class RunModel
{
public:
	RunModel();
	virtual ~RunModel();
	
	// Find the time dependent steady states
	int FindTimeDependentHomogeneousGlucoseSteadyState(const double percentageChange, Model *model);	
	int FindTimeDependentHomogeneousOxygenSteadyState(const double percentageChange, Model *model);
	int FindTimeDependentInhomogeneousGlucoseSteadyState(const double percentageChange, Model *model);
	int FindTimeDependentInhomogeneousOxygenSteadyState(const double percentageChange, Model *model);
	
	// Find the time independent chemical steady states
	int FindTimeIndependentHomogeneousGlucoseSteadyState(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model);
	int FindTimeIndependentHomogeneousOxygenSteadyState(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model);
	int FindTimeIndependentInhomogeneousGlucoseSteadyState(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model);
	int FindTimeIndependentInhomogeneousOxygenSteadyState(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model, CellPopulation *celPop,
		Coordinate source, double sourceConcentration);
	int FindTimeIndependentInhomogeneousTGFBetaSteadyState(const int maxIts,
		const double terminationThreshold, const LatticeType lattice, Model *model, CellPopulation *celPop,
		Coordinate source, double sourceConcentration);
	int FindTimeIndependentInhomogeneousGlucoseSteadyStateWithAvgConsumption(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model);
	int FindTimeIndependentInhomogeneousOxygenSteadyStateWithAvgConsumption(const int maxIts, 
		const double terminationThreshold, const LatticeType lattice, Model *model);
	
	// Lagacy methods, from implement the SOR algorithm
	int FindChemicalsSteadyState(const int maxIts, const double terminationThreshold, Model *model);
	// Added by Marjan 11/3/12
	int FindChemicalsSteadyStateIgnoreConsumption(const int maxIts, const double terminationThreshold,
			Model *model, Coordinate glucoseSource, double sourceConcentration,
			const LatticeType lattice);
	// Added by Marjan 11/3/12

	// Added by Marjan 11/4/12
	int FindChemicalsSteadyStateIgnoreConsumption(const int maxIts, const double terminationThreshold,
			Model *model, Coordinate glucoseSource, Coordinate secondGlucoseSource, double firstSourceConcentration, double secondSourceConcentration,
			const LatticeType lattice);

	// Added by Marjan 11/4/12

	// Added by Marjan 17/3/12
	double chebechev(double previousOmega, int iteration);
	// Added by Marjan 17/3/12

	//Added by Marjan 24/3/12
	void fillInDoubleLatticeForMigration(Model * model);
	//Added by Marjan 24/3/12

	// Added by Marjan 28/3/12
	void migrateTumourCellsUsingDoubleLattice(CellPopulation *cellPop, Model *model);
	// Added by Marjan 28/3/12

	// Added by Marjan 6/4/12
	void migrateCellCalculatingVelocityOfCells(CellPopulation *cellPop, Model *model, double deltaT, int timeStep);
	// Added by Marjan 6/4/12

	// Added by Marjan 25/4/12
	void cellEvolution(CellPopulation *cellPop, Model *model, double chemicalSource, double cellEvolutionParameter);

	// Added 8/5/12
	void pressureInduced(Model *model, CellPopulation *cellPop, Coordinate chemical, double prositiy, int timeStep);
	// Added 8/5/12

	int TestSORAlgorithm(const int maxIts, const double terminationThreshold, Model *model);
	
	void CellularUpdate(const double hoursPerCellUpdate, const LatticeType lattice, Model *model);
	
	// Cell migreation methods
	double AdaptivelyMigrateAllCells(const double maxCellMigrationDistance, const double lambda, Model *model);
	double MigrateAllCells(const double lambda, const double deltaT, Model *model);
	double MigrateCells(const int y, const int x, const double lambda, const double deltaT,
	const double maxCellMigrationDistance, double &cellsMigrated, Model *model);
	
	// Old homogeneous diffusion method
	void HomogeneousDiffuseChemicals(Model *model);
	// Old chemical consumption method
	void ConsumeChemicals(Model *model);
	
	void DisplayModel(const DisplayType display, const LatticeType lattice, Model *model);
	
	void SaveModel(const char *postfix, const int t, const DisplayType display, const LatticeType lattice,
		Model *model);
	void SaveModel(const int t, const DisplayType display, const LatticeType lattice, Model *model);
	void SaveResults(const LatticeType lattice, ofstream &radiusFile, ofstream &maxDensityFile, 
		ofstream &totalPopulationFile, ofstream &patPopulationFile, ofstream &pantPopulationFile,
		ofstream &qatPopulationFile, ofstream &qantPopulationFile, ofstream &dPopulationFile, Model *model);
	void SaveResults(const LatticeType lattice, ofstream &widthFile, ofstream &heightFile,
		ofstream &widthHeightRatioFile, Model *model);
		
	// Has the tumour reached the edge of the tumour growth region? 
	bool TumourAtGrowthRegionEdge(const LatticeType lattice, Model *model);		
	
	//Added by Marjan 21/3/12
	vector<long int>findCellKeysWithinEachLatticeSite(Model *model, CellPopulation *cellPop, Coordinate coord, int deltaX, int deltaY);
	//Added by Marjan 21/3/12

private:

	// Model to file methods
	string GenerateFileNamePrefix(const int t, const DisplayType display, const LatticeType lattice);
	void Save(const string fileName, const DisplayType display, const LatticeType lattice, Model *model);
	void Save(ofstream &fileHandle, const DisplayType display, const LatticeType lattice, Model *model);

	// Lagacy method, from the implementation of the SOR algorithm
	double FindResidWithBC(int i, int j, ChemicalType chemical, Model *model);

	void EmptyLattice(const LatticeType lattice, Model *model);
	void SwapLattices(Model *model);
	
	void TimeStep(Model *model);
	bool TimeStep(const double percentageChange, ChemicalType chemical, Model *model);
	
	// Old homogeneous diffusion methods
	double OrthogonalSum(int y, int x, ChemicalType chemical, Model *model);
	double OrthogonalSumTopLeft(int y, int x, ChemicalType chemical, Model *model);
	double OrthogonalSumTopRight(int y, int x, ChemicalType chemical, Model *model);
	double OrthogonalSumBottomLeft(int y, int x, ChemicalType chemical, Model *model);
	double OrthogonalSumBottomRight(int y, int x, ChemicalType chemical, Model *model);	
	
	// The time interval over which migration takes place, initially 30 [s]
	double m_migrationInterval;		
};

#endif /*RUNMODEL_H_*/
