/*
 * Cell.h
 *
 *  Created on: Oct 4, 2011
 *      Author: marjan
 */

#ifndef CELL_H_
#define CELL_H_

#include <vector>

using namespace std;

class Cell {
public:
	Cell();

	//Cell divide(long int cellKey); //Gives birth to a daughter cell
	Cell divide(long int cellKey, double y, double x); // Gives birth to a daughter cell

	void updateRealAge(void);
	void setRealAge(int);
	int getRealAge(void);

	void updatePregnancyAge(void);
	int getPregnancyAge(void);

	int getLifeSpan(void);
	void setLifeSpan(int);

	void resetPregnancyAgeAfterGivingBirth(void);

	bool getLifeStatusOfCell(void);
	void setLifeStatusOfCell(bool);

	bool getJustBorn(void);
	void setJustBorn(bool);

	void setCellKey(long int);
	long int getCellKey(void);
	long int m_cellKey;
	double cellPosition[2];
	bool divideImmedialtely;
	bool pressureInduced;
	bool chosenForPressure;
	//bool migrationFlag;

	virtual ~Cell();

private:
	short int timeBetweenDivision;
	short int lifeSpan;
	int pregnancyAge;
	int realAge;
	bool alive;
	bool justBorn;
};

#endif /* CELL_H_ */
