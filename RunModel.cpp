#include "RunModel.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <algorithm>

/*
 * Default constructor 
 */
RunModel::RunModel(): m_migrationInterval(30.0) {
}

/*
 * Destructor
 */
RunModel::~RunModel() {
}

/*
 * Calculate the contribution from the orthogonal sites on the lattice corners
 */
double RunModel::OrthogonalSumTopLeft(int y, int x, ChemicalType chemical, Model *model) {
	double sum = 0.0;
	if (chemical == GLUCOSE)
		sum = model->m_latticeT[y][x+1]->GetGlucoseConcentration() + model->m_latticeT[y+1][x]->GetGlucoseConcentration();
	else if (chemical == OXYGEN)
		sum = model->m_latticeT[y][x+1]->GetOxygenConcentration() + model->m_latticeT[y+1][x]->GetOxygenConcentration();	
	return sum;
}

double RunModel::OrthogonalSumTopRight(int y, int x, ChemicalType chemical, Model *model) {
	double sum = 0.0;
	if (chemical == GLUCOSE)
		sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y+1][x]->GetGlucoseConcentration();
	else if (chemical == OXYGEN)
		sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y+1][x]->GetOxygenConcentration();	
	return sum;
}

double RunModel::OrthogonalSumBottomLeft(int y, int x, ChemicalType chemical, Model *model) {
	double sum = 0.0;
	if (chemical == GLUCOSE) 
		sum = model->m_latticeT[y-1][x]->GetGlucoseConcentration() + model->m_latticeT[y][x+1]->GetGlucoseConcentration();
	else if (chemical == OXYGEN)
		sum = model->m_latticeT[y-1][x]->GetOxygenConcentration() + model->m_latticeT[y][x+1]->GetOxygenConcentration();
	return sum;
}

double RunModel::OrthogonalSumBottomRight(int y, int x, ChemicalType chemical, Model *model) {
	double sum = 0;
	if (chemical == GLUCOSE)
		sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y-1][x]->GetGlucoseConcentration();
	else if (chemical == OXYGEN)
		 sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y-1][x]->GetOxygenConcentration();
	return sum;
}


/*
 * Calculate the sum of a site's orthogonal neighbours
 */
double RunModel::OrthogonalSum(int y, int x, ChemicalType chemical, Model *model) {
	
	double sum = 0.0;
	// Is the site on the lattice edge?
	if (x > 0 && y > 0 && x < WIDTH-1 && y < HEIGHT-1) {
		// No
		if (chemical == GLUCOSE)
			sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y][x+1]->GetGlucoseConcentration() + 
				model->m_latticeT[y-1][x]->GetGlucoseConcentration() + model->m_latticeT[y+1][x]->GetGlucoseConcentration();
		else if (chemical == OXYGEN)
			sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y][x+1]->GetOxygenConcentration() + 
				model->m_latticeT[y-1][x]->GetOxygenConcentration() + model->m_latticeT[y+1][x]->GetOxygenConcentration();			
	} 
	else {
		// Is the site on the left edge?
		if (x == 0) {
			// Is the site in the corner?
			if (y == 0 || y == HEIGHT-1) {
				if (y == 0) {
					sum = OrthogonalSumTopLeft(y,x, chemical, model);						
				}
				else if (y == HEIGHT-1) {
					sum = OrthogonalSumBottomLeft(y,x, chemical, model);
				}
			}
			else {
				if (chemical == GLUCOSE)
					sum = model->m_latticeT[y][x+1]->GetGlucoseConcentration() + model->m_latticeT[y-1][x]->GetGlucoseConcentration() + 
						model->m_latticeT[y+1][x]->GetGlucoseConcentration();
				else if (chemical == OXYGEN)
					sum = model->m_latticeT[y][x+1]->GetOxygenConcentration() + model->m_latticeT[y-1][x]->GetOxygenConcentration() + 
						model->m_latticeT[y+1][x]->GetOxygenConcentration();				
			}
		}
		// Is the site on the right edge?
		else if (x == WIDTH-1) {
			// Is the site in the corner?
			if (y == 0 || y == HEIGHT-1) {
				if (y == 0) {
					sum = OrthogonalSumTopRight(y, x, chemical, model);
				}
				else if (y == HEIGHT-1) {
					sum = OrthogonalSumBottomRight(y, x, chemical, model);
				}
			}
			else {
				if (chemical == GLUCOSE)
					sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y-1][x]->GetGlucoseConcentration() + 
						model->m_latticeT[y+1][x]->GetGlucoseConcentration();
				else if (chemical == OXYGEN) 
					sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y-1][x]->GetOxygenConcentration() + 
						model->m_latticeT[y+1][x]->GetOxygenConcentration();
			}
		}	
		// Is the site on the top edge?
		else if (y == 0) {
			if (x == 0 || x == WIDTH-1) {
				if (x == 0) {
					sum = OrthogonalSumTopLeft(y,x, chemical, model);
				}
				else if (x == WIDTH-1) {
					sum = OrthogonalSumTopRight(y,x, chemical, model);
				}
			}
			else {
				if (chemical == GLUCOSE)
					sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y][x+1]->GetGlucoseConcentration() + 
						model->m_latticeT[y+1][x]->GetGlucoseConcentration();
				else if (chemical == OXYGEN)
					sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y][x+1]->GetOxygenConcentration() + 
						model->m_latticeT[y+1][x]->GetOxygenConcentration();				
			}
		} 	
		// Is the site on the bottom edge?
		else if (y == HEIGHT-1) {
			if (x == 0 || x == WIDTH-1) {
				if (x == 0) {
					sum = OrthogonalSumBottomLeft(y,x, chemical, model);
				}
				else if (x == WIDTH-1) {
					sum = OrthogonalSumBottomRight(y,x, chemical, model);
				}
			}
			else {
				if (chemical == GLUCOSE)
					sum = model->m_latticeT[y][x-1]->GetGlucoseConcentration() + model->m_latticeT[y][x+1]->GetGlucoseConcentration() + 
						model->m_latticeT[y-1][x]->GetGlucoseConcentration();
				else if (chemical == OXYGEN)
					sum = model->m_latticeT[y][x-1]->GetOxygenConcentration() + model->m_latticeT[y][x+1]->GetOxygenConcentration() + 
						model->m_latticeT[y-1][x]->GetOxygenConcentration();					
			}
		}
	}
	return sum;
}

/*
 * Copy the lattice at t+1 to t, effectively incrementing the time step 
 */
void RunModel::TimeStep(Model *model) {	
	for (int y = 0; y < HEIGHT; ++y) {
		for (int x = 0; x < WIDTH; ++x) {
			double g = model->m_latticeT1[y][x]->GetGlucoseConcentration();
			model->m_latticeT[y][x]->SetGlucoseConcentration(g);
			double o = model->m_latticeT1[y][x]->GetOxygenConcentration();
			model->m_latticeT[y][x]->SetOxygenConcentration(o);
		}
	}
}

/*
 * Copy the lattice at t+1 to t, effectively incrementing the time step 
 * 
 * Return: true is the concentration field is in a steady state
 */
bool RunModel::TimeStep(const double percentageChange, ChemicalType chemical, Model *model) {
	
	bool inSteadyState = true;
	for (int y = 0; y < HEIGHT; ++y) {
		for (int x = 0; x < WIDTH; ++x) {
			if (chemical == GLUCOSE) {
				// Update the glucose concentration
				double g1 = model->m_latticeT1[y][x]->GetGlucoseConcentration();				
				if (inSteadyState) {
					// Find the change in glucose concentration
					double g = model->m_latticeT[y][x]->GetGlucoseConcentration();
					double onePercent = g / 100.0;
					double threshold = onePercent * percentageChange;
					double change = abs(g - g1);
					if (change > threshold)
						inSteadyState = false;
				}
				model->m_latticeT[y][x]->SetGlucoseConcentration(g1);
			}
			else if (chemical == OXYGEN) {
				// Update the oxygen concentration
				double o1 = model->m_latticeT1[y][x]->GetOxygenConcentration();
				if (inSteadyState) {
					// Find the change in oxygen concentration
					double o = model->m_latticeT[y][x]->GetOxygenConcentration();
					double onePercent = o / 100.0;
					double threshold = onePercent * percentageChange;
					double change = abs(o - o1);
					if (change > threshold)
						inSteadyState = false;
				}
				model->m_latticeT[y][x]->SetOxygenConcentration(o1);
			}	
		}
	} 
	return inSteadyState;
}

/*
 * Implements the Successive Over-relaxation (SOR) algorithm to find the glucose field's steady state for the 
 * homogeneous diffusion of glucose.
 * 
 * maxIts: the maximum number of iteration to find the steady state in.
 * terminationThreshold: the stopping criteria for the chemical concentration
 * lattice: the lattice (T or T1) to calculate the steady state on.
 * model: the model object.
 * Return: the number of iterations.
 */
 int RunModel::FindTimeIndependentHomogeneousGlucoseSteadyState(const int maxIts, const double terminationThreshold,
 	const LatticeType lattice, Model *model) {

	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1); 		
 		
 	double gResidVector, normGlucoseResidVector, omega = 1.9;
 	
 	for (int n = 0; n < maxIts; ++n) {
 		normGlucoseResidVector = 0.0;
 		
 		vector<Coordinate>::const_iterator it;
 		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 			int x = (*it).x;
 			int y = (*it).y;
 			
			// Solve for the glucose concentration
			double g = (*l)[y][x]->GetGlucoseConcentration();
			double Dg = (*l)[y][x]->GetApparentGlucoseDiffusionCoefficient();
			double gE = (-4 
				// PAT consumption
				- (((*l)[y][x]->GetScaledPATGlucoseConsumptionRate() * (*l)[y][x]->GetPATCellPopulation()) / Dg)
				// PANT consumption
				- (((*l)[y][x]->GetScaledPANTGlucoseConsumptionRate() * (*l)[y][x]->GetPANTCellPopulation()) / Dg)
				// QAT consumption
				- (((*l)[y][x]->GetScaledQATGlucoseConsumptionRate() * (*l)[y][x]->GetQATCellPopulation()) / Dg)
				// QANT consumption
				- (((*l)[y][x]->GetScaledQANTGlucoseConsumptionRate() * (*l)[y][x]->GetQANTCellPopulation()) / Dg)
				// Dead cell consumption
				- (((*l)[y][x]->GetScaledDGlucoseConsumptionRate() * (*l)[y][x]->GetDCellPopulation()) / Dg));
			
			gResidVector = (*l)[y+1][x]->GetGlucoseConcentration() 
				+ (*l)[y-1][x]->GetGlucoseConcentration() 
				+ (*l)[y][x+1]->GetGlucoseConcentration()
				+ (*l)[y][x-1]->GetGlucoseConcentration()
				+ (gE * g); 

			normGlucoseResidVector += fabs(gResidVector);
			double gIts = gResidVector / gE;
			g -= omega * gIts;
			(*l)[y][x]->SetGlucoseConcentration(g); 			
 		}
 		// Termination check
 		if (normGlucoseResidVector < terminationThreshold) {
 			return n;
 		}
 	}
 	return maxIts;
 }
 
/*
 * Implements the SOR algorithm to model the non-linear diffusion of glucose and find the field's steady state
 * 
 * maxIts: the maximum number of iteration to find the steady state in.
 * terminationThreshold: the stopping criteria for the chemical concentration
 * lattice: the lattice (T or T1) to calculate the steady state on.
 * model: the model object.
 * Return: the number of iterations.
 */
int RunModel::FindTimeIndependentInhomogeneousGlucoseSteadyState(const int maxIts, 
	const double terminationThreshold, const LatticeType lattice,  Model *model) {
		
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);

	// Find the max site density for checking when the site has a diffusion coefficient of zero
	int maxDensity = 0;
	switch (model->m_porosityEqu) {
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
			throw string("Error: porosity equation not recognised.");				
	} 
		
 	double gResidVector, normGlucoseResidVector, omega = 1.9;
 	
 	for (int n = 0; n < maxIts; ++n) {
 		normGlucoseResidVector = 0.0;
 		
 		vector<Coordinate>::const_iterator it;
 		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 			int x = (*it).x;
 			int y = (*it).y;

			// Solve for the glucose concentration
			double g = (*l)[y][x]->GetGlucoseConcentration();
			double Dg = (*l)[y][x]->GetApparentGlucoseDiffusionCoefficient();
			double a = (0.5 * ((*l)[y][x+1]->GetApparentGlucoseDiffusionCoefficient() + Dg));
			double b = (0.5 * (Dg + (*l)[y][x-1]->GetApparentGlucoseDiffusionCoefficient()));
			double c = (0.5 * ((*l)[y+1][x]->GetApparentGlucoseDiffusionCoefficient() + Dg));
			double d = (0.5 * (Dg + (*l)[y-1][x]->GetApparentGlucoseDiffusionCoefficient())); 


			// Does an adjacent lattice site have the max density
			// Is this site on the left border of a wall (the site to the right is a wall)
			if ((*l)[y][x+1]->GetCellularDensity() >= maxDensity) {
				a = 0;
			}
			// Is this site on the right border of a wall (the site to the left is a wall)
			if ((*l)[y][x-1]->GetCellularDensity() >= maxDensity) {
				b = 0;
			}
			// Is this site on the bottom border of a wall (the site below a wall)
			if ((*l)[y+1][x]->GetCellularDensity() >= maxDensity) {
				c = 0;
			}
			// Is this site on the top border of a wall (the site above a wall)
			if ((*l)[y-1][x]->GetCellularDensity() >= maxDensity) {
				d = 0;
			}

			double consumption = 
				// PAT consumption
				((*l)[y][x]->GetScaledPATGlucoseConsumptionRate() * (*l)[y][x]->GetPATCellPopulation())
				// PANT consumption
				+ ((*l)[y][x]->GetScaledPANTGlucoseConsumptionRate() * (*l)[y][x]->GetPANTCellPopulation())
				// QAT consumption
				+ ((*l)[y][x]->GetScaledQATGlucoseConsumptionRate() * (*l)[y][x]->GetQATCellPopulation())
				// QANT consumption
				+ ((*l)[y][x]->GetScaledQANTGlucoseConsumptionRate() * (*l)[y][x]->GetQANTCellPopulation())
				// Dead cell consumption
				+ ((*l)[y][x]->GetScaledDGlucoseConsumptionRate() * (*l)[y][x]->GetDCellPopulation()); 
			double e = - a - b - c - d - consumption; 

			// Is this site inside a wall?
			if ((*l)[y][x]->GetCellularDensity() >= maxDensity) {
				a = b = c = d = 0;
				e = 1;
				// Yes - the concentration does not change
				//(*l)[y][x]->SetGlucoseConcentration(g);
			}
			
			gResidVector = (a * (*l)[y][x+1]->GetGlucoseConcentration()) 
				+ (b * (*l)[y][x-1]->GetGlucoseConcentration()) 
				+ (c * (*l)[y+1][x]->GetGlucoseConcentration())
				+ (d * (*l)[y-1][x]->GetGlucoseConcentration())
				+ (e * g); 

			normGlucoseResidVector += fabs(gResidVector);
			double gIts = gResidVector / e;
			g -= omega * gIts;
			(*l)[y][x]->SetGlucoseConcentration(g);
 		}
 		// Termination check
 		if (normGlucoseResidVector < terminationThreshold) {
 			return n;
 		}
 	}
 	return maxIts;	
}
	 
/*
 * Implements the Successive Over-relaxation (SOR) algorithm to find the oxygen field's steady state
 * 
 * maxIts: the maximum number of iteration to find the steady state in.
 * terminationThreshold: the stopping criteria for the chemical concentration
 * lattice: the lattice (T or T1) to calculate the steady state on.
 * model: the model object.
 * Return: the number of iterations.
 */
int RunModel::FindTimeIndependentHomogeneousOxygenSteadyState(const int maxIts, 
	const double terminationThreshold, const LatticeType lattice, Model *model) {

	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1); 		
	
	double oResidVector, normOxygenResidVector, omega = 1.9;
	
	for (int n = 0; n < maxIts; ++n) {
		normOxygenResidVector = 0.0;
		
		vector<Coordinate>::const_iterator it;
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;
			
			// Solve for the oxygen concentration
			double o = (*l)[y][x]->GetOxygenConcentration();
			double Do = (*l)[y][x]->GetApparentOxygenDiffusionCoefficient();
			double oE = (-4
				// PAT consumption 
				- (((*l)[y][x]->GetScaledPATOxygenConsumptionRate() * (*l)[y][x]->GetPATCellPopulation()) / Do)
				// PANT consumption
				- (((*l)[y][x]->GetScaledPANTOxygenConsumptionRate() * (*l)[y][x]->GetPANTCellPopulation()) / Do)
				// QAT consumption
				- (((*l)[y][x]->GetScaledQATOxygenConsumptionRate() * (*l)[y][x]->GetQATCellPopulation()) / Do)
				// QANT consumption
				- (((*l)[y][x]->GetScaledQANTOxygenConsumptionRate() * (*l)[y][x]->GetQANTCellPopulation()) / Do)
				// Dead cell consumption
				- (((*l)[y][x]->GetScaledDOxygenConsumptionRate() * (*l)[y][x]->GetDCellPopulation()) / Do));
						
			oResidVector = (*l)[y+1][x]->GetOxygenConcentration() 
				+ (*l)[y-1][x]->GetOxygenConcentration() 
				+ (*l)[y][x+1]->GetOxygenConcentration()
				+ (*l)[y][x-1]->GetOxygenConcentration()
				+ (oE * o);
				
			normOxygenResidVector += fabs(oResidVector);
			double oIts = oResidVector / oE;
			o -= omega * oIts;
			(*l)[y][x]->SetOxygenConcentration(o);
		}
		// Termination check
		if (normOxygenResidVector < terminationThreshold) {
			return n;
		}
	}
	return maxIts;
} 

/*
 * Implements the SOR algorithm to model the non-linear diffusion of glucose and find the field's steady state
 * 
 * maxIts: the maximum number of iteration to find the steady state in.
 * terminationThreshold: the stopping criteria for the chemical concentration
 * lattice: the lattice (T or T1) to calculate the steady state on.
 * model: the model object.
 * Return: the number of iterations.
 */
int RunModel::FindTimeIndependentInhomogeneousOxygenSteadyState(const int maxIts, 
	const double terminationThreshold, const LatticeType lattice, Model *model, CellPopulation *celPop,
	Coordinate source, double sourceConcentration) {
		
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1); 	

	// Find the max site density for checking when the site has a diffusion coefficient of zero
	int maxDensity = 0;
	switch (model->m_porosityEqu) {
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
			throw string("Error: porosity equation not recognised.");				
	} 	
		
 	double oResidVector, normOxygenResidVector, omega = 0.0;

	map<long int, Cell>::iterator kk;
	long int numberOfCell = 0;

	(*l)[source.y][source.x]->SetGlucoseConcentration(sourceConcentration);
	for (int n = 0; n < maxIts; ++n) {
 		normOxygenResidVector = 0.0;
 		omega = chebechev(omega, n);
 		
 		vector<Coordinate>::const_iterator it;
 		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 			int x = (*it).x;
 			int y = (*it).y;
 			numberOfCell = 0;
 			
 			for(kk = (celPop->cells).begin(); kk != (celPop->cells).end(); kk++)
 			{
 				int thisY = (kk->second).cellPosition[0];
 				int thisX = (kk->second).cellPosition[1];
 				if((thisY == y) && (thisX == x))
 				{
 					numberOfCell++;
 				}
 			}

			// Solve for the glucose concentration
			double o = (*l)[y][x]->GetOxygenConcentration();
			double Do = (*l)[y][x]->GetApparentOxygenDiffusionCoefficient();
			double a = (0.5 * ((*l)[y][x+1]->GetApparentOxygenDiffusionCoefficient() + Do));
			double b = (0.5 * (Do + (*l)[y][x-1]->GetApparentOxygenDiffusionCoefficient()));
			double c = (0.5 * ((*l)[y+1][x]->GetApparentOxygenDiffusionCoefficient() + Do));
			double d = (0.5 * (Do + (*l)[y-1][x]->GetApparentOxygenDiffusionCoefficient())); 

			// Is this site on the left border of a wall (the site to the right of a wall)
			/*if ((*l)[y][x+1]->GetCellularDensity() >= maxDensity) {
				a = 0;
			}
			// Is this site on the right border of a wall (the site to the left of a wall)
			if ((*l)[y][x-1]->GetCellularDensity() >= maxDensity) {
				b = 0;
			} 
			// Is this site on the bottom border of a wall (the site below a wall)
			if ((*l)[y+1][x]->GetCellularDensity() >= maxDensity) {
				c = 0;
			}
			// Is this site on the top border of a wall (the site above a wall)
			if ((*l)[y-1][x]->GetCellularDensity() >= maxDensity) {
				d = 0;
			}*/

			// Hard coded now
			double consumption = (52 * pow(10, -17.0)) * numberOfCell;
			double e = - a - b - c - d - consumption; 

			// Is this site within a wall?
			/*if ((*l)[y][x]->GetCellularDensity() >= maxDensity) {
				a = b = c = d = 0;
				e = 1;
				// Yes - the concentration does not change
				//(*l)[y][x]->SetOxygenConcentration(o);
			}*/
	
			if((x == source.x) && (y == source.y)){
			oResidVector = (a * (*l)[y][x+1]->GetOxygenConcentration()) 
				+ (b * (*l)[y][x-1]->GetOxygenConcentration()) 
				+ (c * (*l)[y+1][x]->GetOxygenConcentration())
				+ (d * (*l)[y-1][x]->GetOxygenConcentration())
				+ (e * o) + sourceConcentration;
			}
			else
			{
				oResidVector = (a * (*l)[y][x+1]->GetOxygenConcentration())
				+ (b * (*l)[y][x-1]->GetOxygenConcentration())
				+ (c * (*l)[y+1][x]->GetOxygenConcentration())
				+ (d * (*l)[y-1][x]->GetOxygenConcentration())
				+ (e * o);
			}
			normOxygenResidVector += fabs(oResidVector);
			double oIts = oResidVector / e;
			o -= omega * oIts;
			(*l)[y][x]->SetOxygenConcentration(o); 			
 		}
 		// Termination check
 		if (normOxygenResidVector < terminationThreshold) {
 			return n;
 		}
 	}
 	return maxIts;	
}
// Added by Marjan for TGFBeta consumption
int RunModel::FindTimeIndependentInhomogeneousTGFBetaSteadyState(const int maxIts,
		const double terminationThreshold, const LatticeType lattice, Model *model, CellPopulation *celPop,
		Coordinate source, double sourceConcentration){

	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);

 	double oResidVector, normOxygenResidVector, omega = 0.0;

	map<long int, Cell>::iterator kk;
	long int numberOfCell = 0;

	(*l)[source.y][source.x]->SetOxygenConcentration(sourceConcentration);
	for (int n = 0; n < maxIts; ++n) {
 		normOxygenResidVector = 0.0;
 		omega = chebechev(omega, n);

 		vector<Coordinate>::const_iterator it;
 		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 			int x = (*it).x;
 			int y = (*it).y;
 			numberOfCell = 0;

 			for(kk = (celPop->cells).begin(); kk != (celPop->cells).end(); kk++)
 			{
 				int thisY = (kk->second).cellPosition[0];
 				int thisX = (kk->second).cellPosition[1];
 				if((thisY == y) && (thisX == x))
 				{
 					numberOfCell++;
 				}
 			}

			// Solve for the glucose concentration
			double o = (*l)[y][x]->GetOxygenConcentration();
			double Do = (*l)[y][x]->GetApparentOxygenDiffusionCoefficient();
			double a = (0.5 * ((*l)[y][x+1]->GetApparentOxygenDiffusionCoefficient() + Do));
			double b = (0.5 * (Do + (*l)[y][x-1]->GetApparentOxygenDiffusionCoefficient()));
			double c = (0.5 * ((*l)[y+1][x]->GetApparentOxygenDiffusionCoefficient() + Do));
			double d = (0.5 * (Do + (*l)[y-1][x]->GetApparentOxygenDiffusionCoefficient()));

			// Hard coded now
			double consumption = (52 * pow(10, -14.0)) * numberOfCell;
			double e = - a - b - c - d - consumption;

			if((x == source.x) && (y == source.y)){
			oResidVector = (a * (*l)[y][x+1]->GetOxygenConcentration())
				+ (b * (*l)[y][x-1]->GetOxygenConcentration())
				+ (c * (*l)[y+1][x]->GetOxygenConcentration())
				+ (d * (*l)[y-1][x]->GetOxygenConcentration())
				+ (e * o) + sourceConcentration;
			}
			else
			{
				oResidVector = (a * (*l)[y][x+1]->GetOxygenConcentration())
				+ (b * (*l)[y][x-1]->GetOxygenConcentration())
				+ (c * (*l)[y+1][x]->GetOxygenConcentration())
				+ (d * (*l)[y-1][x]->GetOxygenConcentration())
				+ (e * o);
			}
			normOxygenResidVector += fabs(oResidVector);
			double oIts = oResidVector / e;
			o -= omega * oIts;
			(*l)[y][x]->SetOxygenConcentration(o);
 		}
 		// Termination check
 		if (normOxygenResidVector < terminationThreshold) {
 			return n;
 		}
 	}
 	return maxIts;
}

/*
 * Implements the Successive Over-relaxation (SOR) algorithm to find the chemical
 * fields' steady states.
 * 
 * Note: This is legacy code used in the implementation of the SOR algorithm, please ignore.
 */
int RunModel::FindChemicalsSteadyState(const int maxIts, const double terminationThreshold, Model *model) {
	
	double gResidVector, normGlucoseResidVector, oResidVector, normOxygenResidVector, omega = 1.9;
	
	for (int n = 0; n < maxIts; ++n) {
		normGlucoseResidVector = 0.0;
		normOxygenResidVector = 0.0;
		
		vector<Coordinate>::const_iterator it;
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {		
			int x = (*it).x;
			int y = (*it).y; 
			
			// Solve for the glucose concentration
			double g = model->m_latticeT1[y][x]->GetGlucoseConcentration();
			double Dg = model->m_latticeT1[y][x]->GetApparentGlucoseDiffusionCoefficient();
			double gE = (-4 
				// PAT consumption
				- ((model->m_latticeT1[y][x]->GetScaledPATGlucoseConsumptionRate() * model->m_latticeT1[y][x]->GetPATCellPopulation()) / Dg)
				// PANT consumption
				- ((model->m_latticeT1[y][x]->GetScaledPANTGlucoseConsumptionRate() * model->m_latticeT1[y][x]->GetPANTCellPopulation()) / Dg)
				// QAT consumption
				- ((model->m_latticeT1[y][x]->GetScaledQATGlucoseConsumptionRate() * model->m_latticeT1[y][x]->GetQATCellPopulation()) / Dg)
				// QANT consumption
				- ((model->m_latticeT1[y][x]->GetScaledQANTGlucoseConsumptionRate() * model->m_latticeT1[y][x]->GetQANTCellPopulation()) / Dg)
				// Dead cell consumption
				- ((model->m_latticeT1[y][x]->GetScaledDGlucoseConsumptionRate() * model->m_latticeT1[y][x]->GetDCellPopulation()) / Dg));
			
			gResidVector = model->m_latticeT1[y+1][x]->GetGlucoseConcentration() 
				+ model->m_latticeT1[y-1][x]->GetGlucoseConcentration() 
				+ model->m_latticeT1[y][x+1]->GetGlucoseConcentration()
				+ model->m_latticeT1[y][x-1]->GetGlucoseConcentration()
				+ (gE * g); 

			normGlucoseResidVector += fabs(gResidVector);
			double gIts = gResidVector / gE;
			//u[i][j] -= omega * resid / e[i][j];
			g -= omega * gIts;
			model->m_latticeT1[y][x]->SetGlucoseConcentration(g);
				
			// Solve for the oxygen concentration
			double o = model->m_latticeT1[y][x]->GetOxygenConcentration();
			double Do = model->m_latticeT1[y][x]->GetApparentOxygenDiffusionCoefficient();
			double oE = (-4
				// PAT consumption 
				- ((model->m_latticeT1[y][x]->GetScaledPATOxygenConsumptionRate() * model->m_latticeT1[y][x]->GetPATCellPopulation()) / Do)
				// PANT consumption
				- ((model->m_latticeT1[y][x]->GetScaledPANTOxygenConsumptionRate() * model->m_latticeT1[y][x]->GetPANTCellPopulation()) / Do)
				// QAT consumption
				- ((model->m_latticeT1[y][x]->GetScaledQATOxygenConsumptionRate() * model->m_latticeT1[y][x]->GetQATCellPopulation()) / Do)
				// QANT consumption
				- ((model->m_latticeT1[y][x]->GetScaledQANTOxygenConsumptionRate() * model->m_latticeT1[y][x]->GetQANTCellPopulation()) / Do)
				// Dead cell consumption
				- ((model->m_latticeT1[y][x]->GetScaledDOxygenConsumptionRate() * model->m_latticeT1[y][x]->GetDCellPopulation()) / Do));
					
			oResidVector = model->m_latticeT1[y+1][x]->GetOxygenConcentration() 
				+ model->m_latticeT1[y-1][x]->GetOxygenConcentration() 
				+ model->m_latticeT1[y][x+1]->GetOxygenConcentration()
				+ model->m_latticeT1[y][x-1]->GetOxygenConcentration()
				+ (oE * o);
				
			normOxygenResidVector += fabs(oResidVector);
			double oIts = oResidVector / oE;
			o -= omega * oIts;
			model->m_latticeT1[y][x]->SetOxygenConcentration(o);
						
		}						
		if (normGlucoseResidVector < terminationThreshold && 
			normOxygenResidVector < terminationThreshold) {
			return n;
		}
	}
	return maxIts;		
}
// Added by Marjan 11/4/12
/*
 * In the below function the bone is considered for diffustion
 * for implementing SOR
 */
int RunModel::FindChemicalsSteadyStateIgnoreConsumption(const int maxIts, const double terminationThreshold,
		Model *model, Coordinate glucoseSource, Coordinate secondGlucoseSource, double firstSourceConcentration, double secondSourceConcentration,
		const LatticeType lattice)
{
	LatticeSite *(*l)[HEIGHT][WIDTH];
	LatticeSite *(*l1)[HEIGHT][WIDTH];
	l = &(model->m_latticeT);
	l1 = &(model->m_latticeT1);

	// Find the nonGrowth region and and the boundry to it
	vector<Coordinate> *nonGrowthAreaPlusBoundry = new vector<Coordinate>();
	vector<Coordinate>::const_iterator ii, kk;
	Coordinate coord;

	// Set the glucose concentration 0 everywhere for the new iteration
	for(ii = model->m_tumourGrowthEdgeSites->begin(); ii != model->m_tumourGrowthEdgeSites->end(); ii++){
		int x = (*ii).x;
		int y = (*ii).y;
		(*l)[y][x]->SetGlucoseConcentration(0);
	}
	for(ii = model->m_tumourGrowthEdgeSites->begin(); ii != model->m_tumourGrowthEdgeSites->end(); ii++)
	{
		coord.x = ii->x;
		coord.y = ii->y;
		nonGrowthAreaPlusBoundry->push_back(coord);
	}

	for(kk = model->m_nonTumourGrowthRegion->begin(); kk != model->m_nonTumourGrowthRegion->end(); kk++)
	{
		coord.x = kk->x;
		coord.y = kk->y;
		nonGrowthAreaPlusBoundry->push_back(coord);
	}

	// Make sure that boundary and outside the boundary have been initialised with 0 for glucose
	vector<Coordinate>::const_iterator il;
	for(il = nonGrowthAreaPlusBoundry->begin(); il != nonGrowthAreaPlusBoundry->end(); il++){
		(*l)[il->y][il->x]->SetGlucoseConcentration(0);
	}

	// Set the glucose source
	(*l)[glucoseSource.y][glucoseSource.x]->SetOxygenConcentration(firstSourceConcentration);
	(*l)[secondGlucoseSource.y][secondGlucoseSource.x]->SetOxygenConcentration(secondSourceConcentration);

	double gResidVector, normGlucoseResidVector, omega = 0.0; bool flag = false;
	for (int n = 0; n < maxIts; ++n) {

		normGlucoseResidVector = 0.0;
		omega = chebechev(omega, n);

		vector<Coordinate>::const_iterator it;
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;

			// Solve for the glucose concentration
			double g = (*l)[y][x]->GetOxygenConcentration();
			//cout<< x<< endl<< y<< endl;
			double Dg = (*l)[y][x]->GetApparentGlucoseDiffusionCoefficient();
			double gE = (-4 * Dg);
			double a = Dg * (*l)[y+1][x]->GetOxygenConcentration();
			double b = Dg * (*l)[y-1][x]->GetOxygenConcentration();
			double c = Dg * (*l)[y][x+1]->GetOxygenConcentration();
			double d = Dg * (*l)[y][x-1]->GetOxygenConcentration();
			double e = gE * g;

			// The below greyed out peice of code consideres one source of chemical concentration
			if((x == glucoseSource.x) && (y == glucoseSource.y))
			{
				gResidVector = a + b + c + d + e + firstSourceConcentration; //note that sourceConcentration is the constant variable (f in numerical recepeis)
				//cout<< gResidVector<< endl;
			}
			else if((x == secondGlucoseSource.x) && (y == secondGlucoseSource.y))
			{
				gResidVector = a + b + c + d + e + secondSourceConcentration;
			}
			else{
				gResidVector = a + b + c + d + e;
			}
			normGlucoseResidVector += fabs(gResidVector);
			double gIts = gResidVector / gE;
			g -= omega * gIts;

			(*l)[y][x]->SetOxygenConcentration(g);
			//(*l1)[y][x]->SetGlucoseConcentration(g);
			//model->m_latticeT[y][x]->SetGlucoseConcentration(g);
		}
		if (normGlucoseResidVector < terminationThreshold) {
			//cout<< "normGlucoseResidVector is: "<< normGlucoseResidVector<< '\t'<< gResidVector<< endl;
			return n;
		}
	}

	//cout<< "normGlucoseResidVector is: "<< normGlucoseResidVector<< '\t'<< gResidVector<< endl;
	return maxIts;
}

// Added by Marjan 17/3/12

// Added by Marjan 11/3/12
/*
 * In the below function the bone is considered for diffustion
 * for implementing SOR
 */
int RunModel::FindChemicalsSteadyStateIgnoreConsumption(const int maxIts, const double terminationThreshold,
		Model *model, Coordinate glucoseSource, double sourceConcentration,
		const LatticeType lattice)
{
	LatticeSite *(*l)[HEIGHT][WIDTH];
	LatticeSite *(*l1)[HEIGHT][WIDTH];
	l = &(model->m_latticeT);
	l1 = &(model->m_latticeT1);

	// Find the nonGrowth region and and the boundry to it
	vector<Coordinate> *nonGrowthAreaPlusBoundry = new vector<Coordinate>();
	vector<Coordinate>::const_iterator ii, kk;
	Coordinate coord;

	// Set the glucose concentration 0 everywhere for the new iteration
	for(ii = model->m_tumourGrowthEdgeSites->begin(); ii != model->m_tumourGrowthEdgeSites->end(); ii++){
		int x = (*ii).x;
		int y = (*ii).y;
		(*l)[y][x]->SetGlucoseConcentration(0);
	}
	for(ii = model->m_tumourGrowthEdgeSites->begin(); ii != model->m_tumourGrowthEdgeSites->end(); ii++)
	{
		coord.x = ii->x;
		coord.y = ii->y;
		nonGrowthAreaPlusBoundry->push_back(coord);
	}

	for(kk = model->m_nonTumourGrowthRegion->begin(); kk != model->m_nonTumourGrowthRegion->end(); kk++)
	{
		coord.x = kk->x;
		coord.y = kk->y;
		nonGrowthAreaPlusBoundry->push_back(coord);
	}

	// Make sure that boundary and outside the boundary have been initialised with 0 for glucose
	vector<Coordinate>::const_iterator il;
	for(il = nonGrowthAreaPlusBoundry->begin(); il != nonGrowthAreaPlusBoundry->end(); il++){
		(*l)[il->y][il->x]->SetGlucoseConcentration(0);
	}

	// Set the glucose source
	(*l)[glucoseSource.y][glucoseSource.x]->SetGlucoseConcentration(sourceConcentration);

	double gResidVector, normGlucoseResidVector, omega = 0.0; bool flag = false;
	for (int n = 0; n < maxIts; ++n) {

		normGlucoseResidVector = 0.0;
		omega = chebechev(omega, n);

		vector<Coordinate>::const_iterator it;
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;

			// Solve for the glucose concentration
			double g = (*l)[y][x]->GetGlucoseConcentration();
			//cout<< x<< endl<< y<< endl;
			double Dg = (*l)[y][x]->GetApparentGlucoseDiffusionCoefficient();
			double gE = (-4 * Dg);
			double a = Dg * (*l)[y+1][x]->GetGlucoseConcentration();
			double b = Dg * (*l)[y-1][x]->GetGlucoseConcentration();
			double c = Dg * (*l)[y][x+1]->GetGlucoseConcentration();
			double d = Dg * (*l)[y][x-1]->GetGlucoseConcentration();
			double e = gE * g;

			// The below greyed out peice of code consideres one source of chemical concentration
			if((x == glucoseSource.x) && (y == glucoseSource.y))
			{
				gResidVector = a + b + c + d + e + sourceConcentration; //note that sourceConcentration is the constant variable (f in numerical recepeis)
				//cout<< gResidVector<< endl;
			}
			else
			{
				gResidVector = a + b + c + d + e;
			}
			normGlucoseResidVector += fabs(gResidVector);
			double gIts = gResidVector / gE;
			g -= omega * gIts;

			//This part changed by Marjan
			/*
			 * Keep the updated in another lattice (T1)
			 */
			(*l)[y][x]->SetGlucoseConcentration(g);
			//(*l1)[y][x]->SetGlucoseConcentration(g);
			//model->m_latticeT[y][x]->SetGlucoseConcentration(g);
		}
/*		for(int thisy = 0; thisy < HEIGHT; thisy++){
			for(int thisx = 0; thisx < WIDTH; thisx++){
				double glu = (*l1)[thisy][thisx]->GetGlucoseConcentration();
					(*l)[thisy][thisx]->SetGlucoseConcentration(glu);
				}
		}*/
		if (normGlucoseResidVector < terminationThreshold) {
			//cout<< "normGlucoseResidVector is: "<< normGlucoseResidVector<< '\t'<< gResidVector<< endl;
			return n;
		}
	}

	//cout<< "normGlucoseResidVector is: "<< normGlucoseResidVector<< '\t'<< gResidVector<< endl;
	return maxIts;
}

// Added by Marjan 17/3/12
double RunModel::chebechev(double previousOmega, int iteration)
{
	double newOmega;
	const double jCobi = (cos(PI/HEIGHT) + cos(PI/WIDTH))/2;
	if( iteration == 0)
	{
		newOmega = 1;
	}
	else
	{
		newOmega = 1/(1-((jCobi * jCobi * previousOmega)/4));
	}
	return newOmega;
}
// Added by Marjan 17/3/12

// Added by Marjan 24/3/12
void RunModel::fillInDoubleLatticeForMigration(Model *model)
{
	double temp[HEIGHT][WIDTH];
	double tempDoubleWIDTH[HEIGHT][WIDTH * 2];
	double tempDoubleHEIGHT[HEIGHT * 2][WIDTH];
	//double tempDoubleWIDTHHEIGHT[HEIGHT * 2][WIDTH * 2];

	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			double g = model->m_latticeT[y][x]->GetGlucoseConcentration();
			temp[y][x] = g;
		}
	}

	/*Repeat each site of the lattice once in the
	 * horizontal direction(Make the HEIGHT doubled)
	 */
	int repeatX = 0;
	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH; x++){
			tempDoubleWIDTH[y][repeatX] = temp[y][x];
			repeatX++;
			tempDoubleWIDTH[y][repeatX] = temp[y][x];
			repeatX++;
		}
		repeatX = 0;
	}

	//Swap the column with the row
	for(int y = 0; y < HEIGHT; y++){
		for(int x = 0; x < WIDTH * 2; x++){
			tempDoubleHEIGHT[x][y] = tempDoubleWIDTH[y][x];
		}
	}

	/*Now make a lattice with the new WIDTH of WIDTH*2
	 * and the new HEIGHT of HEIGHT *2
	 */
	repeatX = 0;
	for(int y = 0; y < HEIGHT * 2; y++){
		for(int x = 0; x < WIDTH; x++){
			model->m_doubleLattice[y][repeatX] = tempDoubleHEIGHT[y][x];
			repeatX++;
			model->m_doubleLattice[y][repeatX] = tempDoubleHEIGHT[y][x];
			repeatX++;
		}
		repeatX = 0;
	}
}
// Added by Marjan 24/3/12

// Added by Marjan 28/3/12
/*
void RunModel::migrateTumourCellsUsingDoubleLattice(CellPopulation *cells, Model *model)
{
	vector<Cell>::iterator ii, jj, kk;
	vector<Cell> newCellContainer;
	vector<Cell>::iterator ss;
	vector<long int>::iterator con;
	vector<double> temp;

	for(kk = cells->cellPop->begin(); kk != cells->cellPop->end(); kk++)
	{
		kk->migrationFlag = true;
	}

	for(ii = cells->cellPop->begin(); ii != cells->cellPop->end(); ii++)
	{
		double yPos = ii->cellPosition[0]; double xPos = ii->cellPosition[1];

		//External movement
		for(jj = ii; jj != cells->cellPop->end(); jj++)
		{
			//cout<< jj->cellPosition[0]<< " "<< jj->cellPosition[1]<< " "<< jj->migrationFlag<< endl;
			if((jj->cellPosition[0] == yPos) && (jj->cellPosition[1] == xPos) && ((jj->migrationFlag) == true))
			{
				// Just add the cell keys
				newCellContainer.push_back(*jj);
			}
		}
		// In case the cell container is not empty do as follow
		if(newCellContainer.size() != 0)
		{
			int ySite = yPos; int xSite = xPos;
			map<string, double> nonRep;
			map<string, double>::iterator it;

			double currentChem = model->m_doubleLattice[ySite][xSite];
			double chemUp = model->m_doubleLattice[ySite][xSite - 1];
			double chemDown = model->m_doubleLattice[ySite][xSite + 1];
			double chemRight = model->m_doubleLattice[ySite + 1][xSite];
			double chemLeft = model->m_doubleLattice[ySite - 1][xSite];

			// Add the non repeatetive chems to nonRep
			if((currentChem != chemUp) && ((currentChem < chemUp)))
			{
				nonRep.insert(pair<string, double>("ChemUp", chemUp));
			}
			if((currentChem != chemDown) && ((currentChem < chemDown)))
			{
				nonRep.insert(pair<string, double>("ChemDown", chemDown));
			}
			if((currentChem != chemRight) && ((currentChem < chemRight)))
			{
				nonRep.insert(pair<string, double>("ChemRight", chemRight));
			}
			if((currentChem != chemLeft) && ((currentChem < chemLeft)))
			{
				nonRep.insert(pair<string, double>("ChemLeft", chemLeft));
			}

			if(nonRep.size() != 0)
			{
				double content = 0;
				double up, down, right, left;
				bool flagUp = false; bool flagDown = false; bool flagRight = false; bool flagLeft = false;
				bool flagUpDown = false; bool flagUpRight = false; bool flagUpLeft = false;
				bool flagDownRight = false; bool flagDownLeft = false; bool flagRightLeft = false;
				bool flagJustUp = true; bool flagJustDown = true; bool flagJustRight = true; bool flagJustLeft = true;
				int numberOfCellsMovingUp, numberOfCellsMovingDown, numberOfCellsMovingRight, numberOfCellsMovingLeft;

				//for all the chems in the table
				for(it = nonRep.begin(); it != nonRep.end(); it++)
				{
					if((it->first) == "ChemUp")
					{
						content += it->second;
						flagUp = true;
						up = it->second;
					}
					else if((it->first) == "ChemDown")
					{
						content += it->second;
						flagDown = true;
						down = it->second;
					}
					else if((it->first) == "ChemRight")
					{
						content += it->second;
						flagRight = true;
						right = it->second;
					}
					else if((it->first) == "ChemLeft")
					{
						content += it->second;
						flagLeft = true;
						left = it->second;
					}
				}// end of for
				if(flagUp){
					if(flagUp && flagDown){
						flagUpDown = true;
						numberOfCellsMovingUp = (up/content) * newCellContainer.size();
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();

						int countUp = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countUp < numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
							countUp++;
							if((countUp > numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagUp && flagRight){
						flagUpRight = true;
						numberOfCellsMovingUp = (up/content)* newCellContainer.size();
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();

						int countUp = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countUp < numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
							countUp++;
							if((countUp > numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagUp && flagLeft){
						flagUpLeft = true;
						numberOfCellsMovingUp = (up/content)* newCellContainer.size();
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();

						int countUp = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countUp < numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
							countUp++;
							if((countUp > numberOfCellsMovingUp) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagUpDown || flagUpRight || flagUpLeft){
						flagJustUp = false;
						flagJustDown = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
					if(flagJustUp){
						flagJustDown = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
				}
				if(flagDown){
					if(flagDown && flagUp){
						flagUpDown = true;
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();
						numberOfCellsMovingUp = (up/content)* newCellContainer.size();

						int countDown = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countDown < numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
							countDown++;
							if((countDown > numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagDown && flagRight){
						flagDownRight = true;
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();

						int countDown = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countDown < numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
							countDown++;
							if((countDown > numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagDown && flagLeft){
						flagDownLeft = true;
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();

						int countDown = 0;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countDown < numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
							countDown++;
							if((countDown > numberOfCellsMovingDown) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagDownRight || flagDownLeft || flagUpDown){
						flagJustDown = false;
						flagJustUp = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
					if(flagJustDown){
						flagJustUp = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
				}
				if(flagRight){
					if(flagRight && flagUp){
						flagUpRight = true;
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();
						numberOfCellsMovingUp = (up/content)* newCellContainer.size();

						int countRight;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countRight < numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countRight++;
							if((countRight > numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagRight && flagDown){
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();
						flagDownRight = true;

						int countRight;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countRight < numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countRight++;
							if((countRight > numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagRight && flagLeft){
						flagRightLeft = true;
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();

						int countRight;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countRight < numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countRight++;
							if((countRight > numberOfCellsMovingRight) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagUpRight || flagDownRight || flagDownLeft){
						flagJustUp = false;
						flagJustDown = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
					if(flagJustRight){
						flagJustLeft = false;
						flagJustUp = false;
						flagJustDown = false;
					}
				}
				if(flagLeft){
					if(flagLeft && flagRight){
						flagRightLeft = true;
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();
						numberOfCellsMovingRight = (right/content)* newCellContainer.size();

						int countLeft;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countLeft < numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countLeft++;
							if((countLeft > numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
						}
					}
					if(flagLeft && flagDown){
						flagDownLeft = true;
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();
						numberOfCellsMovingDown = (down/content)* newCellContainer.size();

						int countLeft;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countLeft < numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countLeft++;
							if((countLeft > numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagLeft && flagUp){
						flagUpLeft = true;
						numberOfCellsMovingLeft = (left/content)* newCellContainer.size();
						numberOfCellsMovingUp = (up/content)* newCellContainer.size();

						int countLeft;
						for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
							if((countLeft < numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
								ss->migrationFlag = false;
							}
							countLeft++;
							if((countLeft > numberOfCellsMovingLeft) && (ss->migrationFlag == true)){
								ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
								ss->migrationFlag = false;
							}
						}
					}
					if(flagRightLeft || flagDownLeft || flagUpLeft){
						flagJustUp = false;
						flagJustDown = false;
						flagJustRight = false;
						flagJustLeft = false;
					}
					if(flagJustLeft){
						flagJustUp = false;
						flagJustDown = false;
						flagJustRight = false;
					}
				}
				// Note: Should move proportionally based on the source concentration
				// Just move up
				if(flagJustUp){
					numberOfCellsMovingUp = newCellContainer.size();
					for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
						ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] - 1;
						ss->migrationFlag = false;
					}
				}
				// Just move down
				if(flagJustDown){
					numberOfCellsMovingDown = newCellContainer.size();
					for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
						ss->cellPosition[0] = ss->cellPosition[0]; ss->cellPosition[1] = ss->cellPosition[1] + 1;
						ss->migrationFlag = false;
					}
				}
				// Just move right
				if(flagJustRight){
					numberOfCellsMovingRight = newCellContainer.size();
					for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
						ss->cellPosition[0] = ss->cellPosition[0] + 1; ss->cellPosition[1] = ss->cellPosition[1];
						ss->migrationFlag = false;
					}
				}
				// Just move left
				if(flagJustLeft){
					numberOfCellsMovingLeft = newCellContainer.size();
					for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
						ss->cellPosition[0] = ss->cellPosition[0] - 1; ss->cellPosition[1] = ss->cellPosition[1];
						ss->migrationFlag = false;
					}
				}
			}// end of if
		}
		// Make the migration for the already migrated cells false
		for(kk = cells->cellPop->begin(); kk != cells->cellPop->end(); kk++)
		{
			if((kk->cellPosition[0] == yPos) && (kk->cellPosition[1] == xPos) && ((kk->migrationFlag) == true))
			{
				kk->migrationFlag = false;
			}
		}
		for(ss = newCellContainer.begin(); ss != newCellContainer.end(); ss++){
			cout<< "key: "<< ss->m_cellKey<< " y: "<< ss->cellPosition[0]<< " x: "<< ss->cellPosition[1]<< endl;
		}

		newCellContainer.clear();
	}
}
*/
//Added by Marjan 28/3/12

// Added by Marjan 6/4/12

void RunModel::migrateCellCalculatingVelocityOfCells(CellPopulation *cellPop, Model *model, double deltaT, int timeStep)
{
	ofstream migration("/Users/marjan/Desktop/migration20.txt", ios::app);

	map<long int, Cell>::iterator kk,ll;
	for(kk = cellPop->cells.begin(); kk != cellPop->cells.end(); kk++)
	{
		int y = (kk->second).cellPosition[0]; int x = (kk->second).cellPosition[1];
		//double current = model->m_latticeT[y][x]->GetGlucoseConcentration();
		//v0 = -vx
		double upX = -1 * (model->m_latticeT[y][x-1]->GetOxygenConcentration());
		//v1x = -v1sin45, v1y = v1cos45
		double upRightX = -1 * (model->m_latticeT[y+1][x-1]->GetOxygenConcentration()) * sin(PI/4);
		double upRightY = (model->m_latticeT[y+1][x-1]->GetOxygenConcentration()) * cos(PI/4);
		//v2 = vy
		double rightY = model->m_latticeT[y+1][x]->GetOxygenConcentration();
		//v3x = v3sin45, v3y = v3cos45
		double downRightX = (model->m_latticeT[y+1][x+1]->GetOxygenConcentration()) * sin(PI/4);
		double downRightY = (model->m_latticeT[y+1][x+1]->GetOxygenConcentration()) * cos(PI/4);
		//v4 = vx
		double downX = model->m_latticeT[y][x+1]->GetOxygenConcentration();
		//v5x = v5sin45, v5y = -v5cos45
		double downLeftX = (model->m_latticeT[y-1][x+1]->GetOxygenConcentration()) * sin(PI/4);
		double downLeftY = -1 * (model->m_latticeT[y-1][x+1]->GetOxygenConcentration()) * cos(PI/4);
		//v6 = -vy
		double leftY = -1* (model->m_latticeT[y-1][x]->GetOxygenConcentration());
		// v7x = -v7sin45, v7y = -v7cos45
		double upLeftX = -1 * (model->m_latticeT[y-1][x-1]->GetOxygenConcentration()) * sin(PI/4);
		double upLeftY = -1 * (model->m_latticeT[y-1][x-1]->GetOxygenConcentration()) * cos(PI/4);

		// Calculate the value of the resultant vector
		double zigmaX = upX + upRightX + downRightX + downX + downLeftX + upLeftX;
		double zigmaY = upRightY + rightY + downRightY + downLeftY + leftY + upLeftY;
		//double resultantVector = sqrt(pow(zigmaX, 2.0) + pow(zigmaY, 2.0));

		// Calculate the angle
/*		double tanAngel = tan(zigmaX/zigmaY);
		double angle = atan(tanAngel);
		//cout << "angel is: "<< angle<< endl;

		// Calculate the distance the cell should move
		double l = resultantVector * deltaT;

		// Cell's new position is
		double yDistanceToMove = l * cos(angle);
		double xDistanceToMove = l * sin(angle);*/

		//Another x distances to move x = v(x)t + x0 and y = v(y)t + y0

		double newY =  zigmaY * deltaT + (kk->second).cellPosition[0];
		double newX =  zigmaX * deltaT + (kk->second).cellPosition[1];


		/*double newY = kk->cellPosition[0] + yDistanceToMove;
		double newX = kk->cellPosition[1] + xDistanceToMove;*/

		(kk->second).cellPosition[0] = newY;
		(kk->second).cellPosition[1] = newX;

	}
	if((timeStep % 5) == 0)
	{
		//cout<< timeStep<< endl;
		migration<< "t="<< timeStep<< endl;
		for(ll = cellPop->cells.begin(); ll != cellPop->cells.end(); ll++)
			{
				//migration<<cellPop->cells.size()<< endl;
				//migration<< (ll->second).m_cellKey<< '\t'<< (ll->second).cellPosition[0]<< '\t'<< (ll->second).cellPosition[1]<< endl;
			migration<< (ll->second).cellPosition[0]<< "\t"<< (ll->second).cellPosition[1]<< endl;
			}
		migration<<"----"<< endl;
	}

}

/*
 * Pressure Induced cell movement
 * Distributes cells within lattice
 * sites
 */
void RunModel::pressureInduced(Model *model, CellPopulation *cellPop, Coordinate chemicalCoordinate, double prosity, int timeStep)
{
	ofstream migration("/Users/marjan/Desktop/migration_11.txt", ios::app);
	map<long int, Cell>::iterator kk, mm, ll, ss;

	int sourceLocationY = chemicalCoordinate.y;
	int sourceLocationX = chemicalCoordinate.x;



	for(mm = (cellPop->cells).begin(); mm != (cellPop->cells).end(); mm++)
	{
		(mm->second).pressureInduced = true;
		(mm->second).chosenForPressure = false;
	}

	vector<Coordinate>::const_iterator it;
 	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it)
 	{
 		int x = (*it).x;
 		int y = (*it).y;
 		long int cellCounter = 0;

 		for(kk = (cellPop->cells).begin(); kk != (cellPop->cells).end(); kk++)
 		{
 			int cellY = (kk->second).cellPosition[0];
 			int cellX = (kk->second).cellPosition[1];

 			if(((y == cellY) && (x == cellX)) && (((kk->second).pressureInduced) == true))
 			{
 				cellCounter++;
 				(kk->second).chosenForPressure = true;
 			}
 		}// end of for
 		// Take the correct action
 		if(cellCounter > prosity)
 		{
 			// first determine the distance between the and the cell site
 			vector<int> xDistancesFromTheSource;
 			vector<int> yDistancesFromTheSource;
 			vector<int>::iterator ii, jj;

 			int xDisctaneFromNorthEast = x - 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromNorthEast);
 			int xDisctaneFromEast = x;
 			xDistancesFromTheSource.push_back(xDisctaneFromEast);
 			int xDisctaneFromSouthEast = x + 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromSouthEast);
 			int xDisctaneFromSouth = x + 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromSouth);
 			int xDisctaneFromSouthWest = x + 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromSouthWest);
 			int xDisctaneFromWest = x;
 			xDistancesFromTheSource.push_back(xDisctaneFromWest);
 			int xDisctaneFromWestNorth = x - 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromWestNorth);
 			int xDisctaneFromNorth = x - 1;
 			xDistancesFromTheSource.push_back(xDisctaneFromNorth);

 			int yDisctaneFromNorthEast = y + 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromNorthEast);
 			int yDisctaneFromEast = y + 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromEast);
 			int yDisctaneFromSouthEast = y + 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromSouthEast);
 			int yDisctaneFromSouth = y;
 			yDistancesFromTheSource.push_back(yDisctaneFromSouth);
 			int yDisctaneFromSouthWest = y - 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromSouthWest);
 			int yDisctaneFromWest = y - 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromWest);
 			int yDisctaneFromWestNorth = y - 1;
 			yDistancesFromTheSource.push_back(yDisctaneFromWestNorth);
 			int yDisctaneFromNorth = y;
 			yDistancesFromTheSource.push_back(yDisctaneFromNorth);

 			sort(xDistancesFromTheSource.begin(), xDistancesFromTheSource.end());
 			sort(yDistancesFromTheSource.begin(), yDistancesFromTheSource.end());

 			int largestXDistance = 0;
 			for(ii = xDistancesFromTheSource.begin(); ii != xDistancesFromTheSource.end(); ii++)
 			{
 				if(abs(((*ii) - sourceLocationX)) > largestXDistance)
 				{
 					largestXDistance = *ii;
 				}
 			}// end of for

 			int largestYDistance = 0;
 			for(jj = yDistancesFromTheSource.begin(); jj != yDistancesFromTheSource.end(); jj++)
 			{
 				if(abs(((*jj) - sourceLocationY)) > largestYDistance)
 				{
 					largestYDistance = *jj;
 				}
 			}// end of for
 			//largestYDistance = y + 1;
 			//largestXDistance = x + 1;
 			long int howManyCellWouldMove = cellCounter - prosity;
 			for(ll = (cellPop->cells).begin(); ll != (cellPop->cells).end(); ll++)
 			{
 				if(((ll->second).chosenForPressure == true))
 				{
 					double currentY = (ll->second).cellPosition[0];
 					int takeYInteger = (ll->second).cellPosition[0];
 					double remainY = currentY - takeYInteger;
 					(ll->second).cellPosition[0] = largestYDistance + remainY;

 					double currentX = (ll->second).cellPosition[1];
 					int takeXInteger = (ll->second).cellPosition[1];
 					double remainX = currentX - takeXInteger;
 					(ll->second).cellPosition[1] = largestXDistance + remainX;

 					howManyCellWouldMove--;
 				}
				if(howManyCellWouldMove < 0)
					break;
 			}// end of for
 		}// end of big if

 		/*for(ll = (cellPop->cells).begin(); ll != (cellPop->cells).end(); ll++)
 		{
 			if((ll->second).chosenForPressure == true)
 				(ll->second).pressureInduced = false;
 		}*/

 		// Make the whole cells have the pressureforce of false
 		for(ll = (cellPop->cells).begin(); ll != (cellPop->cells).end(); ll++)
 		{
 			(ll->second).chosenForPressure = false;
 		}
 	}
	if((timeStep % 2) == 0)
	{
		//cout<< timeStep<< endl;
		migration<< "t="<<timeStep<< endl;
		for(ll = cellPop->cells.begin(); ll != cellPop->cells.end(); ll++)
			{
				//migration<<cellPop->cells.size()<< endl;
				//migration<< (ll->second).m_cellKey<< '\t'<< (ll->second).cellPosition[0]<< '\t'<< (ll->second).cellPosition[1]<< endl;
			migration<< (ll->second).cellPosition[0]<< "\t"<< (ll->second).cellPosition[1]<< endl;
			}
		migration<<"--"<< endl;
	}
}

/* 
 * Note: This is legacy code used in the implementation of the SOR algorithm, please ignore.
 */
double RunModel::FindResidWithBC(int i, int j, ChemicalType chemical, Model *model) {
	double resid = 0.0;
	if (chemical == GLUCOSE) {
		double Dg = model->m_latticeT1[i][j]->GetApparentGlucoseDiffusionCoefficient();
		if (i == 0) {
			if (j == 0)
				resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())	
					+ (-2 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration()); 
			else if (j == WIDTH-1)
				resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
					+ (-2 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
			else 
				resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
					+ (-3 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
		}
		else if (i == HEIGHT-1) {
			if (j == 0)
				resid = (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())
					+ (-2 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
			else if (j == WIDTH-1)
				resid = (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
					+ (-2 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
			else
				resid = (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration()) 
					+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())
					+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
					+ (-3 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
		}
		else if (j == 0) {
			resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration()) 
				+ (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration())
				+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())
				+ (-3 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());			
		}
		else if (j == WIDTH-1) {
			resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration()) 
				+ (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration())
				+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
				+ (-3 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
		}
		else
			resid = (Dg * model->m_latticeT1[i+1][j]->GetGlucoseConcentration()) 
				+ (Dg * model->m_latticeT1[i-1][j]->GetGlucoseConcentration())
				+ (Dg * model->m_latticeT1[i][j+1]->GetGlucoseConcentration())
				+ (Dg * model->m_latticeT1[i][j-1]->GetGlucoseConcentration())
				+ (-4 * Dg * model->m_latticeT1[i][j]->GetGlucoseConcentration());
	}
	else if (chemical == OXYGEN) {
		double Do = model->m_latticeT1[i][j]->GetApparentOxygenDiffusionCoefficient();
		if (i == 0) {
			if (j == 0)
				resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())	
					+ (-2 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration()); 
			else if (j == WIDTH-1)
				resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
					+ (-2 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
			else 
				resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
					+ (-3 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
		}
		else if (i == HEIGHT-1) {
			if (j == 0)
				resid = (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())
					+ (-2 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
			else if (j == WIDTH-1)
				resid = (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
					+ (-2 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
			else
				resid = (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration()) 
					+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())
					+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
					+ (-3 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
		}
		else if (j == 0) {
			resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration()) 
				+ (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration())
				+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())
				+ (-3 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());			
		}
		else if (j == WIDTH-1) {
			resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration()) 
				+ (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration())
				+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
				+ (-3 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());
		}
		else
			resid = (Do * model->m_latticeT1[i+1][j]->GetOxygenConcentration()) 
				+ (Do * model->m_latticeT1[i-1][j]->GetOxygenConcentration())
				+ (Do * model->m_latticeT1[i][j+1]->GetOxygenConcentration())
				+ (Do * model->m_latticeT1[i][j-1]->GetOxygenConcentration())
				+ (-4 * Do * model->m_latticeT1[i][j]->GetOxygenConcentration());		
	}
	return resid;	
}

/*
 * Note: This is legacy code used in the implementation of the SOR algorithm, please ignore.
 */ 
int RunModel::TestSORAlgorithm(const int maxIts, const double terminationThreshold, Model *model) {
	
	double gResidVector, normGlucoseResidVector, oResidVector, normOxygenResidVector, omega = 1.5;
	
	for (int n = 0; n < maxIts; ++n) {
		normGlucoseResidVector = 0.0;
		normOxygenResidVector = 0.0;
		
		for (int i = 0; i < HEIGHT; ++i) {
			for (int j = 0; j < WIDTH; ++j) {
				
				double gE = 0.0;
				double oE = 0.0;
				double Dg = model->m_latticeT1[i][j]->GetApparentGlucoseDiffusionCoefficient();
				double Do = model->m_latticeT1[i][j]->GetApparentOxygenDiffusionCoefficient();
				if (i == 0) {
					if (j == 0 || j == WIDTH-1) {
						gE = (-2 * Dg);
						oE = (-2 * Do);
					}
					else {
						gE = (-3 * Dg);
						oE = (-3 * Do);
					} 
				}
				else if (i == HEIGHT-1) {
					if (j == 0 || j == WIDTH-1) {
						gE = (-2 * Dg);
						oE = (-2 * Do);
					}
					else {
						gE = (-3 * Dg);
						oE = (-3 * Do);
					}
				}
				else if (j == 0) {
					gE = (-3 * Dg);
					oE = (-3 * Do);
				}
				else if (j == WIDTH-1) {
					gE = (-3 * Dg);
					oE = (-3 * Do);
				}
				else {
					gE = (-4 * Dg);
					oE = (-4 * Do);
				}
							
				// Solve for the glucose concentration
				double gT = model->m_latticeT1[i][j]->GetGlucoseConcentration();
							
				gResidVector = FindResidWithBC(i , j, GLUCOSE, model);
				normGlucoseResidVector += fabs(gResidVector);
				
				double gIts = gResidVector / gE;
				gT -= omega * gIts;
				model->m_latticeT1[i][j]->SetGlucoseConcentration(gT);
				
				// Solve for the oxygen concentration
				double oT = model->m_latticeT1[i][j]->GetOxygenConcentration();

				oResidVector = FindResidWithBC(i, j, OXYGEN, model); 
				normOxygenResidVector += fabs(oResidVector);
				
				double oIts = oResidVector / oE;
				oT -= omega * oIts;
				model->m_latticeT1[i][j]->SetOxygenConcentration(oT);
			}		
		}						
		if (normGlucoseResidVector < terminationThreshold && 
			normOxygenResidVector < terminationThreshold) {
			return n;
		}
	}
	return maxIts;		
}

/*
 * Diffuse glucose across the lattice until the glucose is in a steady state.
 * 
 * percentageChange: the percentage by which the glucose concentration at a site has to change
 * for it not to be in a steady state.
 * Returns: the number of iterations for the glucose concentration field to reach a steady state. 
 */
int RunModel::FindTimeDependentHomogeneousGlucoseSteadyState(const double percentageChange, Model *model) {
	int n = 0;
	vector<Coordinate>::const_iterator it;

	while (!TimeStep(percentageChange, GLUCOSE, model)) {		
		
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;
			
			// Diffuse glucose
			double gT = model->m_latticeT[y][x]->GetGlucoseConcentration();			
			double gT1 =  gT + (model->m_latticeT[y][x]->GetApparentGlucoseDiffusionCoefficient()
				* (OrthogonalSum(y, x, GLUCOSE, model) - (4.0 * gT)));
			model->m_latticeT1[y][x]->SetGlucoseConcentration(gT1);
			
			// Consume glucose
			model->m_latticeT1[y][x]->ConsumeGlucose();			
		}
		// Increment time
		++n;
	}
	return n;
}

/*
 * Model the inhomogeneous diffusion of glucose across the lattice until the concentration field is in a steady state.
 * 
 * percentageChange: the percentage by which the oxygen concentration at a site has to change
 * for it not to be in a steady state.
 * model: the model object.
 * Return: the number of iterations for the glucose concentration field to reach a steady state.
 */
int RunModel::FindTimeDependentInhomogeneousGlucoseSteadyState(const double percentageChange, Model *model) {
	int n = 0;
	vector<Coordinate>::const_iterator it;
	
	while (!TimeStep(percentageChange, GLUCOSE, model)) {	
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;
			// Diffuse glucose
			double gT = model->m_latticeT[y][x]->GetGlucoseConcentration();
			double Dg = model->m_latticeT[y][x]->GetApparentGlucoseDiffusionCoefficient();	
			double gT1 = gT 
				+ (((0.5 *(model->m_latticeT[y][x+1]->GetApparentGlucoseDiffusionCoefficient() + Dg))
					*(model->m_latticeT[y][x+1]->GetGlucoseConcentration() - gT)) 
				 - ((0.5 *(Dg + model->m_latticeT[y][x-1]->GetApparentGlucoseDiffusionCoefficient()))
				 	*(gT - model->m_latticeT[y][x-1]->GetGlucoseConcentration())))
				+ (((0.5 *(model->m_latticeT[y+1][x]->GetApparentGlucoseDiffusionCoefficient() + Dg))
					*(model->m_latticeT[y+1][x]->GetGlucoseConcentration() - gT))
				 - ((0.5 *(Dg + model->m_latticeT[y-1][x]->GetApparentGlucoseDiffusionCoefficient()))
				 	*(gT - model->m_latticeT[y-1][x]->GetGlucoseConcentration())));
			model->m_latticeT1[y][x]->SetGlucoseConcentration(gT1);
			// Consume glucose
			model->m_latticeT1[y][x]->ConsumeGlucose();
		}
		// Increment time
		++n;
	}
	return n;			
}

/*
 * Diffuse oxygen across the lattice until the oxygen is in a steady state.
 * 
 * percentageChange: the percentage by which the oxygen concentration at a site has to change
 * for it not to be in a steady state.
 * model: the model object.
 * Returns: the number of iterations for the oxygen concentration field to reach a steady state.  
 */
int RunModel::FindTimeDependentHomogeneousOxygenSteadyState(const double percentageChange, Model *model) {
	int n = 0;
	vector<Coordinate>::const_iterator it; 
	while (!TimeStep(percentageChange, OXYGEN, model)) {
		
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;
			
			// Diffuse oxygen
			double oT = model->m_latticeT[y][x]->GetOxygenConcentration();			
			double oT1 = oT + (model->m_latticeT[y][x]->GetApparentOxygenDiffusionCoefficient()
				* (OrthogonalSum(y, x, OXYGEN, model) - (4.0 * oT)));
			model->m_latticeT1[y][x]->SetOxygenConcentration(oT1);
			
			// Consume oxygen
			model->m_latticeT1[y][x]->ConsumeOxygen();	
		}
		// Increment time
		++n;
	}
	return n;
}

/*
 * Model the inhomogeneous diffusion of oxygen across the lattice until the concentration field is in a 
 * steady state
 * 
 * percentageChange: the percentage by which the oxygen concentration at a site has to change
 * for it not to be in a steady state.
 * model: the model object.
 * Return: the number of iterations for the oxygen concentration field to reach a steady state.
 */
int RunModel::FindTimeDependentInhomogeneousOxygenSteadyState(const double percentageChange, Model *model) {
	int n = 0;
	vector<Coordinate>::const_iterator it;
	
	while (!TimeStep(percentageChange, OXYGEN, model)) {	
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
			int x = (*it).x;
			int y = (*it).y;
			// Diffuse oxygen
			double oT = model->m_latticeT[y][x]->GetOxygenConcentration();
			double Do = model->m_latticeT[y][x]->GetApparentOxygenDiffusionCoefficient();	
			double oT1 = oT 
				+ (((0.5 *(model->m_latticeT[y][x+1]->GetApparentOxygenDiffusionCoefficient() + Do))
					*(model->m_latticeT[y][x+1]->GetOxygenConcentration() - oT)) 
				 - ((0.5 *(Do + model->m_latticeT[y][x-1]->GetApparentOxygenDiffusionCoefficient()))
				 	*(oT - model->m_latticeT[y][x-1]->GetOxygenConcentration())))
				+ (((0.5 *(model->m_latticeT[y+1][x]->GetApparentOxygenDiffusionCoefficient() + Do))
					*(model->m_latticeT[y+1][x]->GetOxygenConcentration() - oT))
				 - ((0.5 *(Do + model->m_latticeT[y-1][x]->GetApparentOxygenDiffusionCoefficient()))
				 	*(oT - model->m_latticeT[y-1][x]->GetOxygenConcentration())));
			model->m_latticeT1[y][x]->SetOxygenConcentration(oT1);
			// Consume oxygen
			model->m_latticeT1[y][x]->ConsumeOxygen();
		}
		// Increment time
		++n;
	}
	return n;		
}

/*
 * Diffuse chemical across the lattice using a time dependent implementation.
 */
void RunModel::HomogeneousDiffuseChemicals(Model *model) {
	TimeStep(model);
	vector<Coordinate>::const_iterator it;

	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
		int x = (*it).x;
		int y = (*it).y;
		
		double gT = model->m_latticeT[y][x]->GetGlucoseConcentration();
		double oT = model->m_latticeT[y][x]->GetOxygenConcentration();
		
		double gT1 =  gT + (model->m_latticeT[y][x]->GetApparentGlucoseDiffusionCoefficient()
				* (OrthogonalSum(y, x, GLUCOSE, model) - (4 * gT)));
		model->m_latticeT1[y][x]->SetGlucoseConcentration(gT1);
		
		double oT1 = oT + (model->m_latticeT[y][x]->GetApparentOxygenDiffusionCoefficient()
				* (OrthogonalSum(y, x, OXYGEN, model) - (4 * oT)));
		model->m_latticeT1[y][x]->SetOxygenConcentration(oT1);			
	}
}

/*
 * The consumption of chemical by cells within the tumour growth region using a time dependent implementation.
 */
void RunModel::ConsumeChemicals(Model *model) {
	vector<Coordinate>::const_iterator it;
	
	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
		int x = (*it).x;
		int y = (*it).y;
		
		model->m_latticeT1[y][x]->ConsumeGlucose();
		model->m_latticeT1[y][x]->ConsumeOxygen();		
	}
}

/*
 * Update the state of cells in the tumour growth region.
 * 
 * hoursPerCellUpdate: the length of a cellular update interval in hours.
 * lattice: the lattice (T or T1) to update.
 * model: the model object.
 */
void RunModel::CellularUpdate(const double hoursPerCellUpdate, const LatticeType lattice, Model *model) {
	
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);	
	
	vector<Coordinate>::const_iterator it;
	
	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
		int x = (*it).x;
		int y = (*it).y;
			
		(*l)[y][x]->CellularUpdate(hoursPerCellUpdate);
	}
}

/*
 * Set to zero the cellular population of the tumour growth region.
 * 
 * lattice: the lattice (T or T1) to set to zero.
 * model: the model object.
 */
void RunModel::EmptyLattice(const LatticeType lattice, Model *model) {
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);
		
	vector<Coordinate>::const_iterator it;
	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
		int x = (*it).x;
		int y = (*it).y;
		
		(*l)[y][x]->SetPATCellPopulation(0.0);
		(*l)[y][x]->SetPANTCellPopulation(0.0);
		(*l)[y][x]->SetQATCellPopulation(0.0);
		(*l)[y][x]->SetQANTCellPopulation(0.0);
		(*l)[y][x]->SetDCellPopulation(0.0);		
	}	
}

/*
 * Set the cellular population of T to T1 and set to zero the cellular population of T1
 * 
 * model: the model object.
 */
void RunModel::SwapLattices(Model *model) {
	
	vector<Coordinate>::const_iterator it;
	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
		int x = (*it).x;
		int y = (*it).y;
			
		model->m_latticeT[y][x]->SetPATCellPopulation(model->m_latticeT1[y][x]->GetPATCellPopulation());		
		model->m_latticeT[y][x]->SetPANTCellPopulation(model->m_latticeT1[y][x]->GetPANTCellPopulation());		
		model->m_latticeT[y][x]->SetQATCellPopulation(model->m_latticeT1[y][x]->GetQATCellPopulation());
		model->m_latticeT[y][x]->SetQANTCellPopulation(model->m_latticeT1[y][x]->GetQANTCellPopulation());
		model->m_latticeT[y][x]->SetDCellPopulation(model->m_latticeT1[y][x]->GetDCellPopulation());
		
		model->m_latticeT1[y][x]->SetPATCellPopulation(0.0);
		model->m_latticeT1[y][x]->SetPANTCellPopulation(0.0);
		model->m_latticeT1[y][x]->SetQATCellPopulation(0.0);
		model->m_latticeT1[y][x]->SetQANTCellPopulation(0.0);
		model->m_latticeT1[y][x]->SetDCellPopulation(0.0);
	}
}

/*
 * Adapt the migration interval to dynamically reduce the distance cells migrate to below the
 * maximum cell migration threshold.
 * 
 * maxCellMigrationDistance: the maximum distance a cell can migrate in um.
 * lambda: the viscous drag coefficient.
 * model: the model object.
 * Return: the maximum number of cells migrated.
 */
double RunModel::AdaptivelyMigrateAllCells(const double maxCellMigrationDistance, const double lambda, Model *model) {
	
	double maxCellsMigrated = 0.0;
	double cellMigrationDistance = 0.0;
	// Loop while the cell migration distance is greater than a specified maximum
	do {
		maxCellsMigrated = 0.0;
		cellMigrationDistance = 0.0;
		EmptyLattice(T1, model);
		// For each lattice site in the tumour growth region
		int x, y;
		vector<Coordinate>::const_iterator it;	
		for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end() && 
			cellMigrationDistance < maxCellMigrationDistance; ++it) {			
			x = (*it).x;
			y = (*it).y;
			
			if (x == 47 && y == 51) {
				bool breakPoint = true;
			}
			
			cellMigrationDistance = 0.0;
			double cellsMigrated = 0.0;
			cellMigrationDistance = MigrateCells( y, x, lambda, m_migrationInterval, maxCellMigrationDistance, cellsMigrated, model);
			// Find the maximum number of cells migrated if the migration distance is not too large
			if (cellsMigrated > maxCellsMigrated && cellMigrationDistance < maxCellMigrationDistance)
				maxCellsMigrated = cellsMigrated; 
		}
		if (cellMigrationDistance >= maxCellMigrationDistance) {
			*model->m_logFile << "Warning: a cell's migration distance exceeded the maximum." << endl;						 
			*model->m_logFile << "   Cell migration distance: " << cellMigrationDistance << " [100 um]" << endl;
			*model->m_logFile << "   Maximum migration distance: " << maxCellMigrationDistance << " [100 um]" << endl;
			*model->m_logFile << "   Migration interval: " << m_migrationInterval << " [s]" << endl;
			*model->m_logFile << "   Lattice coordinate: (" << x << ", " << y << ")" << endl;
			*model->m_logFile << "   Lambda: " << lambda << " [kg/s]" << endl;
			
			// Reduce the migration interval
			m_migrationInterval /= 2.0; 
			
			*model->m_logFile << "Correction: the migration interval will be reduced and migration re-calculated." << endl;
			*model->m_logFile << "   New migration interval: " <<  m_migrationInterval << " [s]" << endl;
		} 
	} while (cellMigrationDistance >= maxCellMigrationDistance);
	
	SwapLattices(model);
	
	return maxCellsMigrated;
}

/*
 * Migrate all the cells in the tumour growth region.
 * 
 * lambda: the viscous drag coefficient.
 * deltaT: the cell migration time interval.
 * model: the model object.
 * Return: the maximum number of cells migrated.
 */
double RunModel::MigrateAllCells(const double lambda, const double deltaT, Model *model) {
 	
 	double maxCellsMigrated = 0.0;
 	vector<Coordinate>::const_iterator it;
 	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 		int x = (*it).x;
 		int y = (*it).y;
		
		double cellMigrationDistance = 0.0;
 		double cellsMigrated = 0.0;
 		cellMigrationDistance = MigrateCells( y, x, lambda, deltaT, model->m_latticeSiteLength/2, cellsMigrated, model);
 		
 		if (cellMigrationDistance >= model->m_latticeSiteLength/2) {
 			stringstream ss;
 			ss << "Error: cells have migrated by more than half the length of a lattice site, ";
 			ss << "try the adaptive migreation method.";
 			ss << "   Lattice site: (" << x << ", " << y << ")" << endl;
 			ss << "   Lambda: " << lambda << " [kg/s]" << endl;
 			ss << "   Migration interval: " << deltaT << " [s]" << endl;
 			ss << "   Maximum migration distance: " << model->m_latticeSiteLength/2 << " [100 um]" << endl;
 			throw ss.str();
 		}
 		
 		if (cellsMigrated > maxCellsMigrated)
 			maxCellsMigrated = cellsMigrated; 
 	}
 	
 	SwapLattices(model);
 	
 	return maxCellsMigrated;
 }

/*
 * Migrate the cells in a square of sites.
 * 
 * y: the y coordinate of the bottom left-hand corner of the square.
 * x: the x coordinate of the bottom right-hand corner of the square.
 * lambda: the viscous drag coefficient.
 * deltaT: the cell migration time interval.
 * maxCellMigrationDistance: the maximum distance a cell can migrate.
 * cellMigrated: the number of cells that moved.
 * model: the model object.
 * Return: the migration distance.
 */
double RunModel::MigrateCells(const int y, const int x, const double lambda, const double deltaT, 
	const double maxCellMigrationDistance, double &cellsMigrated,  Model *model) {
	// The lattice coordinates
	int X = 0;
	int Y = 0;
	int X1 = 0;
	int Y1 = 0;
	// Are the migrating sites within the lattice
	if ((x < WIDTH - 1) && (y > 0)) {		
		// The migrating sites are within the lattice
		X = x;
		Y = y;
		X1 = x+1;
		Y1 = y-1;
	}
	// The lattice sites wrap around the lattice
	else if (x == WIDTH - 1 && y == 0) {
		// The lattice site in the top RHS corner
		X = x;
		Y = y;
		X1 = 0;
		Y1 = HEIGHT-1;
	}
	else if (x < WIDTH - 1 && y == 0) {
		// The lattice sites in the first row
		X = x;
		Y = y;
		X1 = x+1;
		Y1 = HEIGHT - 1;
	} 
	else if (x == WIDTH - 1 && y > 0) {
		// The lattice sites in the last column
		X = x;
		Y = y;
		X1 = 0;
		Y1 = y-1; 
	} 
	
	if (X == 0 && Y == 0 && X1 == 0 && Y1 == 0) {
		stringstream ss;
		ss << "Error: failed to initialise the lattice coordinates." << endl;
		throw ss.str();
	}
	
	// Find the pressures
	double pXY = model->m_latticeT[Y][X]->GetPressure();
	double pX1Y = model->m_latticeT[Y][X1]->GetPressure();
	double pXY1 = model->m_latticeT[Y1][X]->GetPressure();
	double pX1Y1 = model->m_latticeT[Y1][X1]->GetPressure();	
	
	// Lattice site length [100 um]
	const double height = 1.0;
	// The distance by which our square of cells moves [100 um]
	double deltaX = 0.0;
	double deltaY = 0.0;
	
	// Calculate the amount by which the cells moved in the x direction
	deltaX = -(height/(2*lambda))*((pX1Y - pXY)+(pX1Y1 - pXY1))*deltaT;
	// Calculate the amount by which the cells moved in the y direction
	deltaY = -(height/(2*lambda))*((pXY1 - pXY)+(pX1Y1 - pX1Y))*deltaT;
	
	// Check the time step is not too great
	if (abs(deltaX) >= maxCellMigrationDistance || abs(deltaY) >= maxCellMigrationDistance) {
		double migrationDistance = 0.0;
		if (abs(deltaX) > abs(deltaY))
			migrationDistance = abs(deltaX);
		else 
			migrationDistance = abs(deltaY);
		return migrationDistance;
	}
	
	// The change in the cellular density at each of the four sites
	double deltaRhoXYPat = 0.0;
	double deltaRhoXYPant = 0.0;
	double deltaRhoXYQat = 0.0;
	double deltaRhoXYQant = 0.0;
	double deltaRhoXYD = 0.0;
	
	double deltaRhoX1YPat = 0.0;
	double deltaRhoX1YPant = 0.0;
	double deltaRhoX1YQat = 0.0;
	double deltaRhoX1YQant = 0.0;
	double deltaRhoX1YD = 0.0;
	
	double deltaRhoXY1Pat = 0.0;
	double deltaRhoXY1Pant = 0.0;
	double deltaRhoXY1Qat = 0.0;
	double deltaRhoXY1Qant = 0.0;
	double deltaRhoXY1D = 0.0;
	
	double deltaRhoX1Y1Pat = 0.0;
	double deltaRhoX1Y1Pant = 0.0;
	double deltaRhoX1Y1Qat = 0.0;
	double deltaRhoX1Y1Qant = 0.0;
	double deltaRhoX1Y1D = 0.0;
	// The box cell densities
	double aPat = 0.0;
	double aPant = 0.0;
	double aQat = 0.0;
	double aQant = 0.0;
	double aD = 0.0;
	
	double bPat = 0.0;
	double bPant = 0.0;
	double bQat = 0.0;
	double bQant = 0.0;
	double bD = 0.0;
	
	double cPat = 0.0;
	double cPant = 0.0;
	double cQat = 0.0;
	double cQant = 0.0;
	double cD = 0.0;
	
	double dPat = 0.0;
	double dPant = 0.0;
	double dQat = 0.0;
	double dQant = 0.0;
	double dD = 0.0;
	
	double ePat = 0.0;
	double ePant = 0.0;
	double eQat = 0.0;
	double eQant = 0.0;
	double eD = 0.0;
	
	double fPat = 0.0;
	double fPant = 0.0;
	double fQat = 0.0;
	double fQant = 0.0;
	double fD = 0.0;
	
	double gPat = 0.0;
	double gPant = 0.0;
	double gQat = 0.0;
	double gQant = 0.0;
	double gD = 0.0;
	
	double hPat = 0.0;
	double hPant = 0.0;
	double hQat = 0.0;
	double hQant = 0.0;
	double hD = 0.0;
	
	double iPat = 0.0;
	double iPant = 0.0;
	double iQat = 0.0;
	double iQant = 0.0;
	double iD = 0.0;
	
	// Calculate the magnitude of the x and y vector
	double xVectorMagnitude = abs(deltaX);
	double yVectorMagnitude = abs(deltaY);
	// Calculate the area of the boxes of cells moved
	const double aCoef = xVectorMagnitude * yVectorMagnitude;
	const double bCoef = (height/2 - xVectorMagnitude) * yVectorMagnitude;
	const double cCoef = height/2 * yVectorMagnitude;
	const double dCoef = xVectorMagnitude * (height/2 - yVectorMagnitude);
	const double eCoef = xVectorMagnitude * height/2;
	const double fCoef = height/2 * yVectorMagnitude;
	const double gCoef = height/2 * yVectorMagnitude;
	const double hCoef = xVectorMagnitude * height/2;
	const double iCoef = xVectorMagnitude * height/2;
	// In which direction does the square of cells move
	if (deltaX > 0.0 && deltaY > 0.0) {
		// The square of cells moves to the top right
		aPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * aCoef;
		aPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * aCoef;
		aQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * aCoef;
		aQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * aCoef;
		aD = model->m_latticeT[Y][X]->GetDCellPopulation() * aCoef;
		
		bPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * bCoef;
		bPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * bCoef;
		bQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * bCoef;
		bQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * bCoef;
		bD = model->m_latticeT[Y][X]->GetDCellPopulation() * bCoef;
		
		cPat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * cCoef;
		cPant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * cCoef;
		cQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * cCoef;
		cQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * cCoef;
		cD = model->m_latticeT[Y][X1]->GetDCellPopulation() * cCoef;
		
		dPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * dCoef;
		dPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * dCoef;
		dQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * dCoef;
		dQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * dCoef;
		dD = model->m_latticeT[Y][X]->GetDCellPopulation() * dCoef;
		
		ePat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * eCoef;
		ePant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * eCoef;
		eQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * eCoef;
		eQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * eCoef;
		eD = model->m_latticeT[Y1][X]->GetDCellPopulation() * eCoef;
		
		deltaRhoXYPat = - aPat - bPat - dPat;
		deltaRhoXYPant = - aPant - bPant - dPant;
		deltaRhoXYQat = - aQat - bQat - dQat;
		deltaRhoXYQant = - aQant - bQant - dQant;
		deltaRhoXYD = - aD - bD - dD;
		
		deltaRhoX1YPat = - cPat + dPat;
		deltaRhoX1YPant = - cPant + dPant;
		deltaRhoX1YQat = - cQat + dQat;
		deltaRhoX1YQant = - cQant + dQant;
		deltaRhoX1YD = - cD + dD;
		
		deltaRhoXY1Pat = - ePat + bPat;
		deltaRhoXY1Pant = - ePant + bPant;
		deltaRhoXY1Qat = - eQat + bQat;
		deltaRhoXY1Qant = - eQant + bQant;
		deltaRhoXY1D = - eD + bD;
		
		deltaRhoX1Y1Pat = aPat + cPat + ePat;
		deltaRhoX1Y1Pant = aPant + cPant + ePant;
		deltaRhoX1Y1Qat = aQat + cQat + eQat;
		deltaRhoX1Y1Qant = aQant + cQant + eQant;
		deltaRhoX1Y1D = aD + cD + eD;
	}
	else if (deltaX < 0.0 && deltaY > 0.0) {
		// The square of cells moves to the top left
		aPat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * aCoef;
		aPant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * aCoef;
		aQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * aCoef;
		aQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * aCoef;
		aD = model->m_latticeT[Y][X1]->GetDCellPopulation() * aCoef;
		
		bPat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * bCoef;
		bPant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * bCoef;
		bQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * bCoef;
		bQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * bCoef;
		bD = model->m_latticeT[Y][X1]->GetDCellPopulation() * bCoef;
		
		cPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * cCoef;
		cPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * cCoef;
		cQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * cCoef;
		cQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * cCoef;
		cD = model->m_latticeT[Y][X]->GetDCellPopulation() * cCoef;
		
		dPat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * dCoef;
		dPant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * dCoef;
		dQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * dCoef;
		dQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * dCoef;
		dD = model->m_latticeT[Y][X1]->GetDCellPopulation() * dCoef;
		
		ePat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * eCoef;
		ePant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * eCoef;
		eQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * eCoef;
		eQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * eCoef;
		eD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * eCoef;
		
		deltaRhoXYPat = - cPat + dPat;
		deltaRhoXYPant = - cPant + dPant;
		deltaRhoXYQat = - cQat + dQat;
		deltaRhoXYQant = - cQant + dQant;
		deltaRhoXYD = - cD + dD;
		
		deltaRhoX1YPat = - aPat - bPat - dPat;
		deltaRhoX1YPant = - aPant - bPant - dPant;
		deltaRhoX1YQat = - aQat - bQat - dQat;
		deltaRhoX1YQant = - aQant - bQant - dQant;
		deltaRhoX1YD = - aD - bD - dD;
		
		deltaRhoXY1Pat = aPat + ePat + cPat;
		deltaRhoXY1Pant = aPant + ePant + cPant;
		deltaRhoXY1Qat = aQat + eQat + cQat;
		deltaRhoXY1Qant = aQant + eQant + cQant;
		deltaRhoXY1D = aD + eD + cD;
		
		deltaRhoX1Y1Pat = - ePat + bPat;
		deltaRhoX1Y1Pant = - ePant + bPant;
		deltaRhoX1Y1Qat = - eQat + bQat;
		deltaRhoX1Y1Qant = - eQant + bQant;
		deltaRhoX1Y1D = - eD + bD;
	}
	else if (deltaX > 0.0 && deltaY < 0.0) {
		// The square of cells moves to bottom right
		aPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * aCoef;
		aPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * aCoef;
		aQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * aCoef;
		aQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * aCoef;
		aD = model->m_latticeT[Y1][X]->GetDCellPopulation() * aCoef;
		
		bPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * bCoef;
		bPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * bCoef;
		bQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * bCoef;
		bQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * bCoef;
		bD = model->m_latticeT[Y1][X]->GetDCellPopulation() * bCoef;
		
		cPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * cCoef;
		cPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * cCoef;
		cQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * cCoef;
		cQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * cCoef;
		cD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * cCoef;
		
		dPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * dCoef;
		dPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * dCoef;
		dQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * dCoef;
		dQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * dCoef;
		dD = model->m_latticeT[Y1][X]->GetDCellPopulation() * dCoef;
		
		ePat = model->m_latticeT[Y][X]->GetPATCellPopulation() * eCoef;
		ePant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * eCoef;
		eQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * eCoef;
		eQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * eCoef;
		eD = model->m_latticeT[Y][X]->GetDCellPopulation() * eCoef;
		
		deltaRhoXYPat = - ePat + bPat;
		deltaRhoXYPant = - ePant + bPant;
		deltaRhoXYQat = - eQat + bQat;
		deltaRhoXYQant = - eQant + bQant;
		deltaRhoXYD = - eD + bD;
		
		deltaRhoX1YPat = aPat + cPat + ePat;
		deltaRhoX1YPant = aPant + cPant + ePant;
		deltaRhoX1YQat = aQat + cQat + eQat;
		deltaRhoX1YQant = aQant + cQant + eQant;
		deltaRhoX1YD = aD + cD + eD;
		
		deltaRhoXY1Pat = - aPat - bPat - dPat;
		deltaRhoXY1Pant = - aPant - bPant - dPant;
		deltaRhoXY1Qat = - aQat - bQat - dQat;
		deltaRhoXY1Qant = - aQant - bQant - dQant;
		deltaRhoXY1D = - aD - bD - dD;
		
		deltaRhoX1Y1Pat = - cPat + dPat;
		deltaRhoX1Y1Pant = - cPant + dPant;
		deltaRhoX1Y1Qat = - cQat + dQat;
		deltaRhoX1Y1Qant = - cQant + dQant;
		deltaRhoX1Y1D = - cD + dD;
	}
	else if (deltaX < 0.0 && deltaY < 0.0) {
		// The square of cells moves to the bottom left
		aPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * aCoef;
		aPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * aCoef;
		aQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * aCoef;
		aQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * aCoef;
		aD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * aCoef;
		
		bPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * bCoef;
		bPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * bCoef;
		bQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * bCoef;
		bQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * bCoef;
		bD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * bCoef;
		
		cPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * cCoef;
		cPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * cCoef;
		cQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * cCoef;
		cQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * cCoef;
		cD = model->m_latticeT[Y1][X]->GetDCellPopulation() * cCoef;
		
		dPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * dCoef;
		dPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * dCoef;
		dQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * dCoef;
		dQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * dCoef;
		dD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * dCoef;
		
		ePat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * eCoef;
		ePant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * eCoef;
		eQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * eCoef;
		eQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * eCoef;
		eD = model->m_latticeT[Y][X1]->GetDCellPopulation() * eCoef;
		
		deltaRhoXYPat = aPat + cPat + ePat;
		deltaRhoXYPant = aPant + cPant + ePant;
		deltaRhoXYQat = aQat + cQat + eQat;
		deltaRhoXYQant = aQant + cQant + eQant;
		deltaRhoXYD = aD + cD + eD;
		
		deltaRhoX1YPat = - ePat + bPat;
		deltaRhoX1YPant = - ePant + bPant;
		deltaRhoX1YQat = - eQat + bQat;
		deltaRhoX1YQant = - eQant + bQant;
		deltaRhoX1YD = - eD + bD;
		
		deltaRhoXY1Pat = - cPat + dPat;
		deltaRhoXY1Pant = - cPant + dPant;
		deltaRhoXY1Qat = - cQat + dQat;
		deltaRhoXY1Qant = - cQant + dQant;
		deltaRhoXY1D = - cD + dD;
		
		deltaRhoX1Y1Pat = - aPat - bPat - dPat;
		deltaRhoX1Y1Pant = - aPant - bPant - dPant;
		deltaRhoX1Y1Qat = - aQat - bQat - dQat;
		deltaRhoX1Y1Qant = - aQant - bQant - dQant;
		deltaRhoX1Y1D = - aD - bD - dD;
	}
	else if (deltaX == 0.0 && deltaY > 0.0) {
		// The square of cells moves up
		fPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * fCoef;
		fPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * fCoef;
		fQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * fCoef;
		fQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * fCoef;
		fD = model->m_latticeT[Y][X]->GetDCellPopulation() * fCoef;
		
		gPat = model->m_latticeT[y][X1]->GetPATCellPopulation() * fCoef;
		gPant = model->m_latticeT[y][X1]->GetPANTCellPopulation() * fCoef;
		gQat = model->m_latticeT[y][X1]->GetQATCellPopulation() * fCoef;
		gQant = model->m_latticeT[y][X1]->GetQANTCellPopulation() * fCoef;
		gD = model->m_latticeT[y][X1]->GetDCellPopulation() * fCoef;
		
		deltaRhoXYPat = -fPat;
		deltaRhoXYPant = -fPant;
		deltaRhoXYQat = -fQat;
		deltaRhoXYQant = -fQant;
		deltaRhoXYD = -fD;
		
		deltaRhoX1YPat = -gPat;
		deltaRhoX1YPant = -gPant;
		deltaRhoX1YQat = -gQat;
		deltaRhoX1YQant = -gQant;
		deltaRhoX1YD = -gD;
		
		deltaRhoXY1Pat = fPat;
		deltaRhoXY1Pant = fPant;
		deltaRhoXY1Qat = fQat;
		deltaRhoXY1Qant = fQant;
		deltaRhoXY1D = fD;
		
		deltaRhoX1Y1Pat = gPat;
		deltaRhoX1Y1Pant = gPant;
		deltaRhoX1Y1Qat = gQat;
		deltaRhoX1Y1Qant = gQant;
		deltaRhoX1Y1D = gD;
	}
	else if (deltaX == 0.0 && deltaY < 0.0) {
		// The square of cells moves down
		fPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * fCoef;
		fPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * fCoef;
		fQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * fCoef;
		fQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * fCoef;
		fD = model->m_latticeT[Y1][X]->GetDCellPopulation() * fCoef;
		
		gPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * gCoef;
		gPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * gCoef;
		gQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * gCoef;
		gQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * gCoef;
		gD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * gCoef;
		
		deltaRhoXYPat = fPat;
		deltaRhoXYPant = fPant;
		deltaRhoXYQat = fQat;
		deltaRhoXYQant = fQant;
		deltaRhoXYD = fD;
		
		deltaRhoX1YPat = gPat;
		deltaRhoX1YPant = gPant;
		deltaRhoX1YQat = gQat;
		deltaRhoX1YQant = gQant;
		deltaRhoX1YD = gD;
		
		deltaRhoXY1Pat = -fPat;
		deltaRhoXY1Pant = -fPant;
		deltaRhoXY1Qat = -fQat;
		deltaRhoXY1Qant = -fQant;
		deltaRhoXY1D = -fD;
		
		deltaRhoX1Y1Pat = -gPat;
		deltaRhoX1Y1Pant = -gPant;
		deltaRhoX1Y1Qat = -gQat;
		deltaRhoX1Y1Qant = -gQant;
		deltaRhoX1Y1D = -gD;
	}
	else if (deltaX > 0.0 && deltaY == 0.0) {
		// The squares of cells moves to the right
		hPat = model->m_latticeT[Y][X]->GetPATCellPopulation() * hCoef;
		hPant = model->m_latticeT[Y][X]->GetPANTCellPopulation() * hCoef;
		hQat = model->m_latticeT[Y][X]->GetQATCellPopulation() * hCoef;
		hQant = model->m_latticeT[Y][X]->GetQANTCellPopulation() * hCoef;
		hD = model->m_latticeT[Y][X]->GetDCellPopulation() * hCoef;
		
		iPat = model->m_latticeT[Y1][X]->GetPATCellPopulation() * iCoef;
		iPant = model->m_latticeT[Y1][X]->GetPANTCellPopulation() * iCoef;
		iQat = model->m_latticeT[Y1][X]->GetQATCellPopulation() * iCoef;
		iQant = model->m_latticeT[Y1][X]->GetQANTCellPopulation() * iCoef;
		iD = model->m_latticeT[Y1][X]->GetDCellPopulation() * iCoef;
		
		deltaRhoXYPat = -hPat;
		deltaRhoXYPant = -hPant;
		deltaRhoXYQat = -hQat;
		deltaRhoXYQant = -hQant;
		deltaRhoXYD = -hD;
		
		deltaRhoX1YPat = hPat;
		deltaRhoX1YPant = hPant;
		deltaRhoX1YQat = hQat;
		deltaRhoX1YQant = hQant;
		deltaRhoX1YD = hD;
		
		deltaRhoXY1Pat = -iPat;
		deltaRhoXY1Pant = -iPant;
		deltaRhoXY1Qat = -iQat;
		deltaRhoXY1Qant = -iQant;
		deltaRhoXY1D = -iD;
		
		deltaRhoX1Y1Pat = iPat;
		deltaRhoX1Y1Pant = iPant;
		deltaRhoX1Y1Qat = iQat;
		deltaRhoX1Y1Qant = iQant;
		deltaRhoX1Y1D = iD;
	}
	else if (deltaX < 0.0 && deltaY == 0.0) {
		// The square of cells moves to the left
		hPat = model->m_latticeT[Y][X1]->GetPATCellPopulation() * hCoef;
		hPant = model->m_latticeT[Y][X1]->GetPANTCellPopulation() * hCoef;
		hQat = model->m_latticeT[Y][X1]->GetQATCellPopulation() * hCoef;
		hQant = model->m_latticeT[Y][X1]->GetQANTCellPopulation() * hCoef;
		hD = model->m_latticeT[Y][X1]->GetDCellPopulation() * hCoef;
		
		iPat = model->m_latticeT[Y1][X1]->GetPATCellPopulation() * iCoef;
		iPant = model->m_latticeT[Y1][X1]->GetPANTCellPopulation() * iCoef;
		iQat = model->m_latticeT[Y1][X1]->GetQATCellPopulation() * iCoef;
		iQant = model->m_latticeT[Y1][X1]->GetQANTCellPopulation() * iCoef;
		iD = model->m_latticeT[Y1][X1]->GetDCellPopulation() * iCoef;
		
		deltaRhoXYPat = hPat;
		deltaRhoXYPant = hPant;
		deltaRhoXYQat = hQat;
		deltaRhoXYQant = hQant;
		deltaRhoXYD = hD;
		
		deltaRhoX1YPat = -hPat;
		deltaRhoX1YPant = -hPant;
		deltaRhoX1YQat = -hQat;
		deltaRhoX1YQant = -hQant;
		deltaRhoX1YD = -hD;
		
		deltaRhoXY1Pat = iPat;
		deltaRhoXY1Pant = iPant;
		deltaRhoXY1Qat = iQat;
		deltaRhoXY1Qant = iQant;
		deltaRhoXY1D = iD;
		
		deltaRhoX1Y1Pat = -iPat;
		deltaRhoX1Y1Pant = -iPant;
		deltaRhoX1Y1Qat = -iQat;
		deltaRhoX1Y1Qant = -iQant;
		deltaRhoX1Y1D = -iD;
	}
	
	// Error check
	if (model->RoundToModelPrecision(deltaRhoXYPat + deltaRhoX1YPat + deltaRhoXY1Pat + deltaRhoX1Y1Pat
		+ deltaRhoXYPant + deltaRhoX1YPant + deltaRhoXY1Pant + deltaRhoX1Y1Pant
		+ deltaRhoXYQat + deltaRhoX1YQat + deltaRhoXY1Qat + deltaRhoX1Y1Qat
		+ deltaRhoXYQant + deltaRhoX1YQant + deltaRhoXY1Qant + deltaRhoX1Y1Qant
		+ deltaRhoXYD + deltaRhoX1YD + deltaRhoXY1D + deltaRhoX1Y1D) != 0.0) {
		stringstream ss;
		ss << "Error: The sum of the change in cellular density is non-zero: " << endl;
		ss << "deltaRhoXYPat: " << deltaRhoXYPat << endl;
		ss << "deltaRhoXYPant: " << deltaRhoXYPant << endl;
		ss << "deltaRhoXYQat: " << deltaRhoXYQat << endl;
		ss << "deltaRhoXYQant: " << deltaRhoXYQant << endl;
		ss << "deltaRhoXYD: " << deltaRhoXYD << endl;
		ss << "deltaRhoX1YPat: " << deltaRhoX1YPat << endl;
		ss << "deltaRhoX1YPant: " << deltaRhoX1YPant << endl;
		ss << "deltaRhoX1YQat: " << deltaRhoX1YQat << endl;
		ss << "deltaRhoX1YQant: " << deltaRhoX1YQant << endl;
		ss << "deltaRhoX1YD: " << deltaRhoX1YD << endl;
		ss << "deltaRhoXY1Pat: " << deltaRhoXY1Pat << endl;
		ss << "deltaRhoXY1Pant: " << deltaRhoXY1Pant << endl;
		ss << "deltaRhoXY1Qat: " << deltaRhoXY1Qat << endl;
		ss << "deltaRhoXY1Qant: " << deltaRhoXY1Qant << endl;
		ss << "deltaRhoXY1D: " << deltaRhoXY1D << endl;
		ss << "deltaRhoX1Y1Pat: " << deltaRhoX1Y1Pat << endl;
		ss << "deltaRhoX1Y1Pant: " << deltaRhoX1Y1Pant << endl;
		ss << "deltaRhoX1Y1Qat: " << deltaRhoX1Y1Qat << endl;
		ss << "deltaRhoX1Y1Qant: " << deltaRhoX1Y1Qant << endl;
		ss << "deltaRhoX1Y1D: " << deltaRhoX1Y1D << endl;
		ss << "Sum: " << deltaRhoXYPat + deltaRhoX1YPat + deltaRhoXY1Pat + deltaRhoX1Y1Pat
		+ deltaRhoXYPant + deltaRhoX1YPant + deltaRhoXY1Pant + deltaRhoX1Y1Pant
		+ deltaRhoXYQat + deltaRhoX1YQat + deltaRhoXY1Qat + deltaRhoX1Y1Qat
		+ deltaRhoXYQant + deltaRhoX1YQant + deltaRhoXY1Qant + deltaRhoX1Y1Qant
		+ deltaRhoXYD + deltaRhoX1YD + deltaRhoXY1D + deltaRhoX1Y1D << endl;
		throw ss.str();
	}		
	
	// Update the number of cells in the next time step
	// Has the site already been partially updated
	if (model->m_latticeT1[Y][X]->GetPATCellPopulation() == 0.0)
		model->m_latticeT1[Y][X]->SetPATCellPopulation(model->m_latticeT[Y][X]->GetPATCellPopulation() + deltaRhoXYPat);
	else 
		model->m_latticeT1[Y][X]->SetPATCellPopulation(model->m_latticeT1[Y][X]->GetPATCellPopulation() + deltaRhoXYPat);
	if (model->m_latticeT1[Y][X]->GetPANTCellPopulation() == 0.0) 
		model->m_latticeT1[Y][X]->SetPANTCellPopulation(model->m_latticeT[Y][X]->GetPANTCellPopulation() + deltaRhoXYPant);
	else 
		model->m_latticeT1[Y][X]->SetPANTCellPopulation(model->m_latticeT1[Y][X]->GetPANTCellPopulation() + deltaRhoXYPant);
	if (model->m_latticeT1[Y][X]->GetQATCellPopulation() == 0.0)
		model->m_latticeT1[Y][X]->SetQATCellPopulation(model->m_latticeT[Y][X]->GetQATCellPopulation() + deltaRhoXYQat);
	else
		model->m_latticeT1[Y][X]->SetQATCellPopulation(model->m_latticeT1[Y][X]->GetQATCellPopulation() + deltaRhoXYQat); 
	if (model->m_latticeT1[Y][X]->GetQANTCellPopulation() == 0.0)
		model->m_latticeT1[Y][X]->SetQANTCellPopulation(model->m_latticeT[Y][X]->GetQANTCellPopulation() + deltaRhoXYQant);
	else
		model->m_latticeT1[Y][X]->SetQANTCellPopulation(model->m_latticeT1[Y][X]->GetQANTCellPopulation() + deltaRhoXYQant);
	if (model->m_latticeT1[Y][X]->GetDCellPopulation() == 0.0)
		model->m_latticeT1[Y][X]->SetDCellPopulation(model->m_latticeT[Y][X]->GetDCellPopulation() + deltaRhoXYD);
	else
		model->m_latticeT1[Y][X]->SetDCellPopulation(model->m_latticeT1[Y][X]->GetDCellPopulation() + deltaRhoXYD);				
		
	// Update the bottom right lattice site
	if (model->m_latticeT1[Y][X1]->GetPATCellPopulation() == 0.0)
		model->m_latticeT1[Y][X1]->SetPATCellPopulation(model->m_latticeT[Y][X1]->GetPATCellPopulation() + deltaRhoX1YPat);
	else 
		model->m_latticeT1[Y][X1]->SetPATCellPopulation(model->m_latticeT1[Y][X1]->GetPATCellPopulation() + deltaRhoX1YPat);
	if (model->m_latticeT1[Y][X1]->GetPANTCellPopulation() == 0.0)
		model->m_latticeT1[Y][X1]->SetPANTCellPopulation(model->m_latticeT[Y][X1]->GetPANTCellPopulation() + deltaRhoX1YPant);
	else 
		model->m_latticeT1[Y][X1]->SetPANTCellPopulation(model->m_latticeT1[Y][X1]->GetPANTCellPopulation() + deltaRhoX1YPant);
	if (model->m_latticeT1[Y][X1]->GetQATCellPopulation() == 0.0)
		model->m_latticeT1[Y][X1]->SetQATCellPopulation(model->m_latticeT[Y][X1]->GetQATCellPopulation() + deltaRhoX1YQat);
	else 
		model->m_latticeT1[Y][X1]->SetQATCellPopulation(model->m_latticeT1[Y][X1]->GetQATCellPopulation() + deltaRhoX1YQat);
	if (model->m_latticeT1[Y][X1]->GetQANTCellPopulation() == 0.0)
		model->m_latticeT1[Y][X1]->SetQANTCellPopulation(model->m_latticeT[Y][X1]->GetQANTCellPopulation() + deltaRhoX1YQant);
	else 
		model->m_latticeT1[Y][X1]->SetQANTCellPopulation(model->m_latticeT1[Y][X1]->GetQANTCellPopulation() + deltaRhoX1YQant);
	if (model->m_latticeT1[Y][X1]->GetDCellPopulation() == 0.0)
		model->m_latticeT1[Y][X1]->SetDCellPopulation(model->m_latticeT[Y][X1]->GetDCellPopulation() + deltaRhoX1YD);
	else 
		model->m_latticeT1[Y][X1]->SetDCellPopulation(model->m_latticeT1[Y][X1]->GetDCellPopulation() + deltaRhoX1YD);							
		
	// Update the top left lattice site
	if (model->m_latticeT1[Y1][X]->GetPATCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X]->SetPATCellPopulation(model->m_latticeT[Y1][X]->GetPATCellPopulation() + deltaRhoXY1Pat);
	else 
		model->m_latticeT1[Y1][X]->SetPATCellPopulation(model->m_latticeT1[Y1][X]->GetPATCellPopulation() + deltaRhoXY1Pat);
	if (model->m_latticeT1[Y1][X]->GetPANTCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X]->SetPANTCellPopulation(model->m_latticeT[Y1][X]->GetPANTCellPopulation() + deltaRhoXY1Pant);
	else 
		model->m_latticeT1[Y1][X]->SetPANTCellPopulation(model->m_latticeT1[Y1][X]->GetPANTCellPopulation() + deltaRhoXY1Pant);
	if (model->m_latticeT1[Y1][X]->GetQATCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X]->SetQATCellPopulation(model->m_latticeT[Y1][X]->GetQATCellPopulation() + deltaRhoXY1Qat);
	else 
		model->m_latticeT1[Y1][X]->SetQATCellPopulation(model->m_latticeT1[Y1][X]->GetQATCellPopulation() + deltaRhoXY1Qat);
	if (model->m_latticeT1[Y1][X]->GetQANTCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X]->SetQANTCellPopulation(model->m_latticeT[Y1][X]->GetQANTCellPopulation() + deltaRhoXY1Qant);
	else 
		model->m_latticeT1[Y1][X]->SetQANTCellPopulation(model->m_latticeT1[Y1][X]->GetQANTCellPopulation() + deltaRhoXY1Qant);
	if (model->m_latticeT1[Y1][X]->GetDCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X]->SetDCellPopulation(model->m_latticeT[Y1][X]->GetDCellPopulation() + deltaRhoXY1D);
	else 
		model->m_latticeT1[Y1][X]->SetDCellPopulation(model->m_latticeT1[Y1][X]->GetDCellPopulation() + deltaRhoXY1D);				
		
	// Update the top right lattice site
	if (model->m_latticeT1[Y1][X1]->GetPATCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X1]->SetPATCellPopulation(model->m_latticeT[Y1][X1]->GetPATCellPopulation() + deltaRhoX1Y1Pat);
	else 
		model->m_latticeT1[Y1][X1]->SetPATCellPopulation(model->m_latticeT1[Y1][X1]->GetPATCellPopulation() + deltaRhoX1Y1Pat);
	if (model->m_latticeT1[Y1][X1]->GetPANTCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X1]->SetPANTCellPopulation(model->m_latticeT[Y1][X1]->GetPANTCellPopulation() + deltaRhoX1Y1Pant);
	else 
		model->m_latticeT1[Y1][X1]->SetPANTCellPopulation(model->m_latticeT1[Y1][X1]->GetPANTCellPopulation() + deltaRhoX1Y1Pant);
	if (model->m_latticeT1[Y1][X1]->GetQATCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X1]->SetQATCellPopulation(model->m_latticeT[Y1][X1]->GetQATCellPopulation() + deltaRhoX1Y1Qat);
	else 
		model->m_latticeT1[Y1][X1]->SetQATCellPopulation(model->m_latticeT1[Y1][X1]->GetQATCellPopulation() + deltaRhoX1Y1Qat);
	if (model->m_latticeT1[Y1][X1]->GetQANTCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X1]->SetQANTCellPopulation(model->m_latticeT[Y1][X1]->GetQANTCellPopulation() + deltaRhoX1Y1Qant);
	else 
		model->m_latticeT1[Y1][X1]->SetQANTCellPopulation(model->m_latticeT1[Y1][X1]->GetQANTCellPopulation() + deltaRhoX1Y1Qant);
	if (model->m_latticeT1[Y1][X1]->GetDCellPopulation() == 0.0) 
		model->m_latticeT1[Y1][X1]->SetDCellPopulation(model->m_latticeT[Y1][X1]->GetDCellPopulation() + deltaRhoX1Y1D);
	else 
		model->m_latticeT1[Y1][X1]->SetDCellPopulation(model->m_latticeT1[Y1][X1]->GetDCellPopulation() + deltaRhoX1Y1D);

	cellsMigrated = aPat + aPant + aQat + aQant + aD + bPat + bPant + bQat + bQant + bD
		+ cPat + cPant + cQat + cQant + cD + dPat + dPant + dQat + dQant + dD
		+ ePat + ePant + eQat + eQant + eD + fPat + fPant + fQat + fQant + fD 
		+ gPat + gPant + gQat + gQant + gD + hPat + hPant + hQat + hQant + hD
		+ iPat + iPant + iQat + iQant + iD;

	double migrationDistance = 0.0;
	if (abs(deltaX) > abs(deltaY))
		migrationDistance = abs(deltaX);
	else 
		migrationDistance = abs(deltaY);
	return migrationDistance;	
}

/*
 * Has the tumour grown to the edge of the tumour growth region?
 * 
 * Return: true - the tumour has grown to the edge of the tumour growth region
 *  false - otherwise
 */
bool RunModel::TumourAtGrowthRegionEdge(const LatticeType lattice, Model *model) {
	bool tumourAtEdge = false;
	
	// From which lattice will the results be saved
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);	
		
	// Foreach tumour growth site	
 	vector<Coordinate>::const_iterator it;
 	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end() && !tumourAtEdge; ++it) {
 		int x = (*it).x;
 		int y = (*it).y; 		
 		if ((*l)[y][x]->GetCellularDensity() > 1.0) {
 			if (find(model->m_tumourGrowthEdgeSites->begin(), model->m_tumourGrowthEdgeSites->end(), *it) 
 				!= model->m_tumourGrowthEdgeSites->end()) {
 				tumourAtEdge = true;
 			}
 		}	  
 	}
 	return tumourAtEdge;
}

/*
 * Save results regarding the tumour's radius, maximum density and metabolic populations to file.
 * 
 * lattice: the lattice to inspect, either T or T1. 
 * ...File: the file handle to save the respective results to. 
 * model: the model object.
 */
void RunModel::SaveResults(const LatticeType lattice, ofstream &radiusFile, ofstream &maxDensityFile, 
		ofstream &totalPopulationFile, ofstream &patPopulationFile, ofstream &pantPopulationFile,
		ofstream &qatPopulationFile, ofstream &qantPopulationFile, ofstream &dPopulationFile, Model *model) {
	
	double maxRadius, maxDensity, totalPopulation, patPop, pantPop, qatPop, qantPop, dPop;
	maxRadius = maxDensity = totalPopulation = patPop = pantPop = qatPop = qantPop = dPop = 0.0;	
		
	// The circle is assumed to be centred on the lattice centre
	int latticeCentreX = WIDTH / 2;
	int latticeCentreY = HEIGHT / 2;		
	
	// From which lattice will the results be saved
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);	
	
	// Foreach tumour growth site	
 	vector<Coordinate>::const_iterator it;
 	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 		int x = (*it).x;
 		int y = (*it).y;
 		
 		// Update the total population
 		totalPopulation += (*l)[y][x]->GetCellularDensity();
 		// Update the sub-populations
 		patPop += (*l)[y][x]->GetPATCellPopulation();
 		pantPop += (*l)[y][x]->GetPANTCellPopulation();
 		qatPop += (*l)[y][x]->GetQATCellPopulation();
 		qantPop += (*l)[y][x]->GetQANTCellPopulation();
 		dPop += (*l)[y][x]->GetDCellPopulation();
 		
 		// Update the maximum density
 		if ((*l)[y][x]->GetCellularDensity() > maxDensity)
 			maxDensity = (*l)[y][x]->GetCellularDensity();
 		
 		// Only consider sites with a cellular density greater than one as having a tumour cell population
 		// for the purposes of calculating maximum radius
 		if ((*l)[y][x]->GetCellularDensity() > 1.0) {
 			double radius = sqrt(pow(x-latticeCentreX, 2.0) + pow(y-latticeCentreY, 2.0));
 			if (radius > maxRadius)
 				maxRadius = radius;
 		}
 	}
 	// Save the results to file
 	radiusFile << maxRadius << "," << endl;
 	maxDensityFile << maxDensity << "," << endl;
 	totalPopulationFile << totalPopulation << "," << endl;
 	patPopulationFile << patPop << "," << endl;
 	pantPopulationFile << pantPop << "," << endl;
 	qatPopulationFile << qatPop << "," << endl;
 	qantPopulationFile << qantPop << "," << endl;
 	dPopulationFile << dPop << "," << endl; 	
}

/*
 * Save results regarding the tumour's width and height to file.
 * 
 * lattice: the lattice to inspect, either T or T1. 
 * ...File: the file handle to save the respective results to. 
 * model: the model object.
 */
void RunModel::SaveResults(const LatticeType lattice, ofstream &widthFile, ofstream &heightFile,
	ofstream &widthHeightRatioFile, Model *model) {
	
	// The lattice site containing cells furthest to the left of the lattice
	double lowestWidth;
	lowestWidth = WIDTH-1;
	// The lattice site containing cells furthest to the right of the lattice
	double greatestWidth;
	greatestWidth = 0;
	// The lattice site containing cells furthest to the top of the lattice
	double lowestHeight;
	lowestHeight = HEIGHT-1;
	// The lattice site containing cells furthest to the bottom of the lattice
	double greatestHeight;
	greatestHeight = 0;
	
	// From which lattice will the results be saved
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);	
	
	// Foreach tumour growth site	
 	vector<Coordinate>::const_iterator it;
 	for (it = model->m_tumourGrowthSites->begin(); it != model->m_tumourGrowthSites->end(); ++it) {
 		int x = (*it).x;
 		int y = (*it).y; 		
 		
 		// Only consider sites with a cellular density greater than one as having a tumour cell population
 		// for the purposes of calculating the tumour's width and height
 		if ((*l)[y][x]->GetCellularDensity() > 1.0) {

			if (x < lowestWidth)
				lowestWidth = x;
			if (x > greatestWidth)
				greatestWidth = x;
			if (y < lowestHeight)
				lowestHeight = y;				
			if (y > greatestHeight)
				greatestHeight = y;				
 		}
 	}
 	// Calculate the tumour's width and height
 	int tumourWidth = (greatestWidth - lowestWidth) + 1;
 	int tumourHeight = (greatestHeight - lowestHeight) + 1;
 	   
 	// Save the results to file
 	widthFile << tumourWidth << "," << endl;
 	heightFile << tumourHeight << "," << endl;
 	widthHeightRatioFile << (double)tumourHeight/(double)tumourWidth << "," << endl; 	
}

/*
 * Print the model to the stdout.
 * 
 * display: the data to display.
 * lattice: the lattice to inspect, either T or T1. 
 * model: the model object.
 */ 
void RunModel::DisplayModel(const DisplayType display, const LatticeType lattice, Model *model) {
	
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);
	
	for (int i = 0; i < HEIGHT; ++i) {
		for (int j = 0; j < WIDTH; ++j) {
			switch (display) {
				case GLUCOSE_CONCENTRATION:
					cout << setw(10) << (*l)[i][j]->GetGlucoseConcentration();
					break;
				case OXYGEN_CONCENTRATION:
					cout << setw(10) << (*l)[i][j]->GetOxygenConcentration();
					break;
				case GLUCOSE_DIFFCOEF:
					cout << setw(10) << (*l)[i][j]->GetApparentGlucoseDiffusionCoefficient();
					break;
				case OXYGEN_DIFFCOEF:
					cout << setw(10) << (*l)[i][j]->GetApparentOxygenDiffusionCoefficient();
					break;
				case GLUCOSE_CONSUMPTION:
					cout << setw(10) << (*l)[i][j]->GetTotalGlucoseConsumption();
					break;
				case OXYGEN_CONSUMPTION:
					cout << setw(10) << (*l)[i][j]->GetTotalOxygenConsumption();
					break;
				case PAT_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetPATCellPopulation();
					break;
				case PANT_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetPANTCellPopulation();
					break;
				case QAT_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetQATCellPopulation();
					break;
				case QANT_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetQANTCellPopulation();
					break;
				case D_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetDCellPopulation();
					break;
				case CELL_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetCellularDensity();
					break;
				case BONE_POPULATION:
					cout << setw(10) << (*l)[i][j]->GetBoneCellPopulation();
					break;
				default:					
					throw string("Error: display type not recognised");				
			}
			// Print the appropriate EOL character			
			if (j < WIDTH - 1)
				cout << ", ";
			else 
				cout << endl;
		}
	}	
}

/*
 * Save the model to file and close the file.
 * 
 * fileName: the file name.
 * display: the data to display.
 * lattice: the lattice to inspect, either T or T1.
 * model: the model object.
 */ 
void RunModel::Save(const string fileName, const DisplayType display, const LatticeType lattice, Model *model) {
	// Open the model file with the supplied file name
	ofstream modelFile(fileName.c_str(), ios::out);
	Save(modelFile, display, lattice, model);
	modelFile.close();
}

/*
 * Save the model to file.
 * 
 * fileHandle: the file to save the model to.
 * display: the data to display.
 * lattice: the lattice to inspect, either T or T1.
 * model: the model object.
 */	
void RunModel::Save(ofstream &fileHandle, const DisplayType display, const LatticeType lattice, Model *model) {
	
	LatticeSite *(*l)[HEIGHT][WIDTH];
	if (lattice == T)
		l = &(model->m_latticeT);
	else if (lattice == T1)
		l = &(model->m_latticeT1);
	
	if (!fileHandle) {
		throw string("Error: model file could not be opened");
	}
	else {
		// Save the model to file
		for (int y = 0; y < HEIGHT; ++y) {
			for (int x = 0; x < WIDTH; ++x) {
				switch (display) { 
					case GLUCOSE_CONCENTRATION:
						fileHandle << (*l)[y][x]->GetGlucoseConcentration();
						//fileHandle<< model->m_doubleLattice[y][x];
						break;
					case OXYGEN_CONCENTRATION:
						fileHandle << (*l)[y][x]->GetOxygenConcentration();
						break;
					case TGFBETA:
						fileHandle << (*l)[y][x]->GetOxygenConcentration();
						//fileHandle<< model->m_doubleLattice[y][x];
						break;
					case GLUCOSE_DIFFCOEF:
						fileHandle << (*l)[y][x]->GetApparentGlucoseDiffusionCoefficient();
						break;
					case OXYGEN_DIFFCOEF:
						fileHandle << (*l)[y][x]->GetApparentOxygenDiffusionCoefficient();
						break;
					case GLUCOSE_CONSUMPTION:
						fileHandle << (*l)[y][x]->GetTotalGlucoseConsumption();
						break;
					case OXYGEN_CONSUMPTION:
						fileHandle << (*l)[y][x]->GetTotalOxygenConsumption();
						break;
					case PAT_POPULATION:
						fileHandle << (*l)[y][x]->GetPATCellPopulation();
						break;
					case PANT_POPULATION:
						fileHandle << (*l)[y][x]->GetPANTCellPopulation();
						break;
					case QAT_POPULATION:
						fileHandle << (*l)[y][x]->GetQATCellPopulation();
						break;
					case QANT_POPULATION:
						fileHandle << (*l)[y][x]->GetQANTCellPopulation();
						break;
					case D_POPULATION:
						fileHandle << (*l)[y][x]->GetDCellPopulation();
						break;
					case CELL_POPULATION:
						fileHandle << (*l)[y][x]->GetCellularDensity();						
						break;
					case BONE_POPULATION:
						fileHandle << (*l)[y][x]->GetBoneCellPopulation();
						break;
					default:
						throw string("Error: display type not recognised");
				}				
				// Print the appropriate EOL character
				if (x != WIDTH - 1)
					fileHandle << ", ";
				else
					fileHandle << endl;
			}
		}
	}	
}

/*
 * Generate a model file name in the standard format
 * 
 * t: the iteration.
 * display: the data type to display.
 * lattice: the lattice to inspect, either T or T1.
 * Return: the file name.
 */ 
string RunModel::GenerateFileNamePrefix(const int t, const DisplayType display, const LatticeType lattice) {
	// Create the file name from the display type
	stringstream fileName;

	switch (display) {
		case GLUCOSE_CONCENTRATION:
			fileName << "glucoseConcen_";
			break;
		case TGFBETA_CONCENTRATION:
			fileName << "tgfBetaConcen_";
			break;
		case OXYGEN_CONCENTRATION:
			fileName << "oxygenConcen_";
			break;
		case GLUCOSE_DIFFCOEF:
			fileName << "glucoseDiffCoef_";
			break;
		case OXYGEN_DIFFCOEF:
			fileName << "oxygenDiffCoef_";
			break;
		case GLUCOSE_CONSUMPTION:
			fileName << "glucoseConsump_";
			break;
		case OXYGEN_CONSUMPTION:
			fileName << "oxygenConsump_";
			break;
		case PAT_POPULATION:
			fileName << "patPop_";
			break;
		case PANT_POPULATION:
			fileName << "pantPop_";
			break;
		case QAT_POPULATION:
			fileName << "qatPop_";
			break;
		case QANT_POPULATION:
			fileName << "qantPop_";
			break;
		case D_POPULATION:
			fileName << "dPop_";
			break;
		case CELL_POPULATION:
			fileName << "cellPop_";
			break;
		case BONE_POPULATION:
			fileName << "bonePop_";
			break;
		default:
			throw string("Error: display type not recognised");				
	}
	
	if (lattice == T)
		fileName << "T_";
	else if (lattice == T1)
		fileName << "T1_"; 
	
	fileName << t << "_";
	return fileName.str();	
}

/*
 * Wrappers for saving the model to file
 */
void RunModel::SaveModel(const char *postfix, const int t, const DisplayType display, const LatticeType lattice, Model *model) {
	stringstream fileName;
	fileName << GenerateFileNamePrefix(t, display, lattice) << postfix << ".csv";
	Save(fileName.str(), display, lattice, model);
}

void RunModel::SaveModel(const int t, const DisplayType display, const LatticeType lattice, Model *model) {
	stringstream fileName;
	fileName << GenerateFileNamePrefix(t, display, lattice) << ".csv";
	Save(fileName.str(), display, lattice, model);
}

// Modified by Marjan 25/4/12
vector<long int>RunModel::findCellKeysWithinEachLatticeSite(Model *model, CellPopulation *cellPop, Coordinate coord, int deltaX, int deltaY)
{
	vector<long int> cellKeys;
	//vector<Cell>::iterator ii;
	map<long int, Cell>::iterator ii;

	// Determine the whole coordinates
	int x0 = coord.x, y0 = coord.y;
	int x1 = x0 + 0, y1 = y0 + 1;
	int x2 = x0 + 1, y2 = y0 + 1;

	for(ii = (cellPop->cells).begin(); ii != (cellPop->cells).end(); ii++)
	{
		double xCell = (ii->second).cellPosition[0];
		double yCell = (ii->second).cellPosition[1];
		if( (xCell > x0) && (xCell < x2) && (yCell > y0) && yCell < y2){
			cellKeys.push_back((ii->second).m_cellKey);
		}
	}

	return cellKeys;
}

// Added by Marjan 25/4/12
void RunModel::cellEvolution(CellPopulation * cellPop, Model *model, double chemicalSource, double cellEvolutionParameter)
{
	LatticeSite *(*l)[HEIGHT][WIDTH];
	l = &(model->m_latticeT);

	double copyLattice[HEIGHT][WIDTH];
	for(int j = 0; j < HEIGHT; j++){
		for(int i = 0; i < WIDTH; i++){
			double chemical = (*l)[j][i]->GetOxygenConcentration();
			copyLattice[j][i] = chemical;
		}
	}
	cellPop->evolve(cellPop, copyLattice, chemicalSource, cellEvolutionParameter);
}

