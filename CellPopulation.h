/*
 * CellPopulation.h
 *
 *  Created on: Oct 4, 2011
 *      Author: marjan
 */

#ifndef CELLPOPULATION_H_
#define CELLPOPULATION_H_

#include "Cell.h"

#include <map>
#include <iostream>
#include <vector>

using namespace std;

class CellPopulation {
public:
	CellPopulation(int initialNumberOfCells, int timeBetweenDivision, int lifeSpan);
	virtual ~CellPopulation();

	//Start evolving the cells
	void evolve(CellPopulation *celPop, double copyLattice[121][121], double chemicalSource, double cellEvolutionParameter);
	void step(CellPopulation *celPop, double lattice[121][121], double chemicalSource, double cellEvolutionParameter);

	int getInitialNumberOfCells(void);
	int getTimeBetweenDivision(void);
	int getLifeSpan(void);
	long int generateNewCellKey(void);
	void generateRandomAge(map <long int, Cell> *cells);
	vector<double> returnCellKeys(void);


	void logFile(int simulationTime, long int numberOfAliveCells, long int numberOfDeadCells);
	void saveInFile(int simulationTime, long int totalNumberOfCells);

	map<long int, Cell> cells;
	//vector<Cell> *cellPop;// added feb 15th
	double m_position[2];
	int m_initialNumberOfCells;
	int m_timeBetweenDivision;
	int m_lifeSpan;
};

#endif /* CELLPOPULATION_H_ */
