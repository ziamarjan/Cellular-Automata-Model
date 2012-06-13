#ifndef COORDINATE_H_
#define COORDINATE_H_

class Coordinate 
{
public:
	Coordinate();
	virtual ~Coordinate();
	
	bool operator==(const Coordinate &rhs) const {
		return (this->x == rhs.x) && (this->y == rhs.y);
	} 
	
	int x;
	int y;
};

#endif /*COORDINATE_H_*/
