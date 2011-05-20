class Point
{
public:
	Point(){x=0.;y=0.;z=0.;};
	// copy operator
	Point operator=(const Point pt) ;
	Point operator+(const Point pt) const;
	//Point operator-(const Point pt) const;
	Point operator*(double m) const;
	Point operator/(double m) const;
	double x,y,z;
};

Point Point::operator=(const Point pt)
{
	x = pt.x;
	y = pt.y;
	z = pt.z;
	return *this;
}
Point Point::operator+(const Point pt) const
{
	Point temp;
	temp.x = x + pt.x;
	temp.y = y + pt.y;
	temp.z = z + pt.z;
	return temp;
}
Point Point::operator*(double m) const
{
	Point temp;
	temp.x = x*m;
	temp.y = y*m;
	temp.z = z*m;
	return temp;
}
Point Point::operator/(double m) const
{
	Point temp;
	temp.x = x/m;
	temp.y = y/m;
	temp.z = z/m;
	return temp;
}
