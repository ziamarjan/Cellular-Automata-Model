#include "Model.h"
#include "InitialiseModel.h"
#include "RunModel.h"
#include "CellPopulation.h"

#include <cmath>
#include <iostream>
#include <vector>
#include <ctime>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cstring>

#define HOURS_PER_DAY 24
#define MODEL_PRECISION_DP 5

using namespace std;

// Chemical constants
// The initial chemical concentrations [mol/(100 um)^3]

// Changed by Marjan 2/3/12
const double glucoseBoundaryConcentration = /*1.01 * pow(10.0, -12);*/ 0.0;
//Changed by Marjan 2/3/12

//const double glucoseBoundaryConcentration = 2.8 * pow(10, -12)/pow(2.0,7.0);
const double oxygenBoundaryConcentration = /*3.5 * pow(10.0,-14)*/0.0;
//const double oxygenBoundaryConcentration = 3.5 * pow(10,-14)/pow(2.0,2.0);
// The boundary chemical concentration at which a fingering morphology is observed [mol/(100 um)^3]
const double fingeringGlucoseBoundaryConcentration = 2.8 * pow(10.0, -12)/pow(2.0,7.0);
const double fingeringOxygenBoundaryConcentration = 3.5 * pow(10.0,-14)/pow(2.0,2.0);

// The glucose consumption rate [mol/(cell)(s)]
const double glucoseConsumptionRate = 52e-17;
// The oxygen consumption rate [mol/(cell)(s)]
const double oxygenConsumptionRate = 8.3e-17;

// ECM Diffusion coefficients [(100 um)^2 / s]
const double D0g = 0.91;
const double D0o = 0.182;

// Interstitium Fluid Pressure [kg/(100 um)(s)^2]
const double ifp = 0.0933257;

// Cell migration parameters
// Viscour drag coefficient [kg/s]
const double lambda = 10.0;
// A lattice site's length [100 um] 
const double latticeSiteLength = 1.0;
// The maximum distance cells can migrate [100 um]
const double maxCellMigrationDistance = 0.5 * latticeSiteLength;

// Tumour parameters
// Initial tumour radius [100 um]
const int initialTumourRadius = 0;
// The cell doubling time [h]
const int mitoticCycleTime = 22;

// Model parameter 
// Cellular update time interval [h]
const double hoursPerCellUpdate = 0.1875;
// The maximum number of simulation days [d]
const double maxSimDays = 30.0;

//Added by Marjan 23/4/12
//Number of days for the cells to evolve
#define EVOLUTIONDAYS 8

//const double maxSimDays = 120.0;

// A method for replacing a character within a string
string &replaceAll(string &context, const string &from, const string &to) {
	size_t lookHere = 0;
	size_t foundHere;
	while((foundHere = context.find(from, lookHere)) != string::npos) {
		context.replace(foundHere, from.size(), to);
		lookHere = foundHere + to.size();
	}
	return context;
}

int main(int argc, char *argv[]) {

	//Added by Marjan 1/3/12
	CellPopulation *celPop = new CellPopulation(2, 4, 10);
	//cout<< "Helllo"<< endl;
	//Added by Marjan 1/3/12

	Model *model = NULL;
	InitialiseModel *initialiseModel = NULL;
	RunModel *runModel = NULL;

	// Generate the log file time
	time_t startTime = time(NULL);
	tm *t = localtime(&startTime);	
	//char *timeStr = asctime(t);
	char timeStr[40];
	strftime(timeStr, 40, "Time_%H_%M_%S_Date_%d_%m", t);
	string runFilePostfix(timeStr);
	//logFilePostfix = replaceAll(logFilePostfix, string(" "), string("_"));
	stringstream logFileName;
	logFileName << "log_" << runFilePostfix << ".txt";
	// Open the log file  
	ofstream logFile(logFileName.str().c_str(), ios::out);
	
	// Generate the tumour composition file
	stringstream tumourCompositionFileName;
	tumourCompositionFileName << "tumourComposition_" << runFilePostfix;
	
	try {
		// Check the log file was opened correctly
		if (!logFile)
			throw string("Error: failed to open the log file.");
		// Error checking */
		if (HEIGHT != WIDTH)
			throw string("Error: The lattice is not square in shape.");
			
		// Set the metabolic thresholds [mol/(cell)(s)]				
		ChemicalConsumptionRates metabolicThresholds;
		metabolicThresholds.patGlucose = 1.8e-16;
		metabolicThresholds.pantGlucose = 5.2e-16;
		metabolicThresholds.qatGlucose = 0;
		//metabolicThresholds.qantGlucose = 4.3e-16;
		metabolicThresholds.qantGlucose = 5.2e-17;
		metabolicThresholds.dGlucose = 0;
		metabolicThresholds.patOxygen = 8.3e-17;
		metabolicThresholds.pantOxygen = 0; 
		metabolicThresholds.qatOxygen = 0;
		metabolicThresholds.qantOxygen = 0;		
		metabolicThresholds.dOxygen = 0;
		
		// Calculate the scaled consumption rate [(100 um)^3/(cell)(s)]
		ChemicalConsumptionRates scaledConsumptionRates;
		scaledConsumptionRates.patGlucose = 
			Model::ScaledChemicalConsumptionRate(glucoseConsumptionRate, fingeringGlucoseBoundaryConcentration);
		scaledConsumptionRates.pantGlucose = 
			Model::ScaledChemicalConsumptionRate(glucoseConsumptionRate, fingeringGlucoseBoundaryConcentration);
		scaledConsumptionRates.qatGlucose = 
			Model::ScaledChemicalConsumptionRate(glucoseConsumptionRate, fingeringGlucoseBoundaryConcentration);
		scaledConsumptionRates.qantGlucose = 
			Model::ScaledChemicalConsumptionRate(glucoseConsumptionRate, fingeringGlucoseBoundaryConcentration);
		scaledConsumptionRates.dGlucose = 
			Model::ScaledChemicalConsumptionRate(0, fingeringGlucoseBoundaryConcentration);
		scaledConsumptionRates.patOxygen = 
			Model::ScaledChemicalConsumptionRate(oxygenConsumptionRate, fingeringOxygenBoundaryConcentration);	
		scaledConsumptionRates.pantOxygen = 
			Model::ScaledChemicalConsumptionRate(oxygenConsumptionRate, fingeringOxygenBoundaryConcentration);
		scaledConsumptionRates.qatOxygen = 
			Model::ScaledChemicalConsumptionRate(oxygenConsumptionRate, fingeringOxygenBoundaryConcentration);
		scaledConsumptionRates.qantOxygen = 
			Model::ScaledChemicalConsumptionRate(oxygenConsumptionRate, fingeringOxygenBoundaryConcentration);
		scaledConsumptionRates.dOxygen = 
			Model::ScaledChemicalConsumptionRate(0, fingeringOxygenBoundaryConcentration);		

		// Initialise the model chemical consumption rates
		model = new Model(glucoseBoundaryConcentration, oxygenBoundaryConcentration,
			metabolicThresholds, scaledConsumptionRates, mitoticCycleTime, HYP_POROSITY_300, MODEL_PRECISION_DP, 
			latticeSiteLength, D0g, D0o, ifp, initialTumourRadius, hoursPerCellUpdate, 
			maxSimDays, logFile);
		
		initialiseModel = new InitialiseModel();
			
		// Find a circular tumour growth region
		int tumourGrowthRadius = (WIDTH / 2) - 1;
		Coordinate tumourGrowthRegionCentre;
		tumourGrowthRegionCentre.x = WIDTH / 2;
		tumourGrowthRegionCentre.y = HEIGHT / 2;
		initialiseModel->FindTumourGrowthRegion(tumourGrowthRegionCentre, tumourGrowthRadius, model);

		//Added by Marjan 8/2/12
		//Find noncircular region in which the tumour deos not grow within
		initialiseModel->FindTumourNonGrowthRegion(model);
		//Added by Marjan 8/2/12

		// Initialse a bone within the growth region 
		Coordinate topLeftCorner, bottomRightCorner;
		//topLeftCorner.y = 51; 
		topLeftCorner.y = 45;
		topLeftCorner.x = 45;
		//bottomRightCorner.y = 53; 
		bottomRightCorner.y = 63;
		bottomRightCorner.x = 55;
		initialiseModel->CreateBone(topLeftCorner, bottomRightCorner, T, model); 


		// Added by Marjan 6/4/12
		/*
		 * Initialize the cells with
		 * random positions within the
		 * growth region
		 */
		celPop->cells[0].cellPosition[0] = 82.3; celPop->cells[0].cellPosition[1] = 80.6;
		celPop->cells[1].cellPosition[0] = 50.1; celPop->cells[1].cellPosition[1] = 50.3;
		/*celPop->cells[2].cellPosition[0] = 55.2; celPop->cells[1].cellPosition[1] = 55.2;
		celPop->cells[3].cellPosition[0] = 55.2; celPop->cells[1].cellPosition[1] = 55.2;
		celPop->cells[4].cellPosition[0] = 55.2; celPop->cells[1].cellPosition[1] = 55.2;*/
		//initialiseModel->SetTTheumourComposition(celPop);
		//Added by Marjan 1/3/12
		
		runModel = new RunModel();

		//Added by Marjan
		/*
		 * Perform the box counting
		 */
		/*vector<long int> cellKeys;
		Coordinate cellsCoord;
		cellsCoord.x = 57;
		cellsCoord.y = 57;
		cellKeys = runModel->findCellKeysWithinEachLatticeSite(model, celPop, cellsCoord, 1, 1);
		vector<long int>::iterator hh;
		for(hh = cellKeys.begin(); hh != cellKeys.end(); hh++)
			cout<< *hh<< endl;
		 */
		//Added by Marjan 7/3/12
		Coordinate tgfBetaSource;
		tgfBetaSource.y = HEIGHT/2;
		tgfBetaSource.x = WIDTH/2;
		double tgfBetaSourceConcentration = 1.3*pow(10.0, -12);
		//Added by Marjan 7/3/12

		// Added by Marjan 11/4/12
		Coordinate secondGlucoseSource;
		secondGlucoseSource.y = 90;
		secondGlucoseSource.x = 90;
		double secondGlucoseSourceConcentration = 1.3*pow(10.0, -12);

		const int maxIts = 1000;
		const double terminationThreshold = 1.0e-20;

		// Simulate tumour growth while the number of simulated days is less than the maximum number of simulated days
		// and the tumour has not reached the edge of the growth region

		bool tumourAtGrowthRegionEdge = false;
		const int updatesPerDay = HOURS_PER_DAY/hoursPerCellUpdate;
		for (int t = 1; t <=/*= (maxSimDays*updatesPerDay)*/ 36 && !tumourAtGrowthRegionEdge; ++t) {
			
			// Diffuse chemical
			int numIts = 0;	
			/*if ((numIts = runModel->FindChemicalsSteadyStateIgnoreConsumption(maxIts, terminationThreshold, model,
					                                                          tgfBetaSource, tgfBetaSourceConcentration, T)) < maxIts) {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, GLUCOSE_CONCENTRATION, T, model);
			}
			else {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, GLUCOSE_CONCENTRATION, T, model);
				stringstream ss;
				ss << "Error: Maximum iteration exceeded finding glucose steady state." << endl;
				ss << "   Number of iterations: " << numIts << endl;			
				throw ss.str();
			}*/

			/*if ((numIts = runModel->FindTimeIndependentInhomogeneousTGFBetaSteadyState(maxIts, terminationThreshold, T,
					model, celPop, tgfBetaSource, tgfBetaSourceConcentration)) < maxIts) {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, TGFBETA_CONCENTRATION, T, model);
			}
			else {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, TGFBETA_CONCENTRATION, T, model);
				stringstream ss;
				ss << "Error: Maximum iteration exceeded finding tgf beta steady state." << endl;
				ss << "   Number of iterations: " << numIts << endl;
				throw ss.str();
			}*/
			if ((numIts = runModel->FindChemicalsSteadyStateIgnoreConsumption(maxIts, terminationThreshold, model, tgfBetaSource, secondGlucoseSource,
					tgfBetaSourceConcentration, secondGlucoseSourceConcentration, T)) < maxIts) {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, TGFBETA_CONCENTRATION, T, model);
			}
			else {
				runModel->SaveModel(tumourCompositionFileName.str().c_str(), t, TGFBETA_CONCENTRATION, T, model);
				stringstream ss;
				ss << "Error: Maximum iteration exceeded finding tgf beta steady state." << endl;
				ss << "   Number of iterations: " << numIts << endl;			
				throw ss.str();
			}


			//Evolve the system
			runModel->cellEvolution(celPop, model, tgfBetaSourceConcentration, 0.34);
			
			// Migrate the cells
			runModel->migrateCellCalculatingVelocityOfCells(celPop, model, 2000000000000, t);

			// Distribute the cells
			//runModel->pressureInduced(model, celPop, tgfBetaSource, 20, t);
			// Has the tumour reached the edge of the growth region (ECM)?
			//tumourAtGrowthRegionEdge = runModel->TumourAtGrowthRegionEdge(T, model);
		} 
			
		// Clean up
		delete model;
		delete initialiseModel;
		delete runModel;				
	}
	catch (string &errStr) {
		cerr << errStr << endl;
		logFile << errStr << endl;
		// Clean up
		delete model;
		delete initialiseModel;
		delete runModel;		
	}
	catch (bad_alloc &e) {
		cerr << "Error allocating memory: " << e.what() << endl;
		logFile << "Error allocating memory: " << e.what() << endl;
		// Clean up
		delete model;
		delete initialiseModel;
		delete runModel;
	}
	// Catch all errors
	catch (...) {
		cerr << "Error: unkown." << endl;
		logFile << "Error: unkown." << endl;
		// Clean up
		delete model;
		delete initialiseModel;
		delete runModel;
	}
	
	// Save the execution time and close the log file 
	if (logFile) { 
		time_t endTime = time(NULL);
		double diff = difftime(endTime, startTime);
		logFile << "Execution time: " << diff << " [s]" << endl;
		logFile.close();
	}
}
