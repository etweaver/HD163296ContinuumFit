/*
 *  geometry.h
 *  Created by Erik on 11/9/14.
 *  
 *  This header defines the geometry involved in the integral calculation,
 *  including vectors, planes, lines, cones, and spheres, and their intersections.
 *  Path length through a cell is included here, possibly temporarily.
 *  
 *  Note that for the spherical coordinates, r is radial component,
 *  theta is the polar angle, and phi is the azimuthal angle.
 *  
 */

#include <cmath>
#include <cstddef>
#include <iostream>
#include <vector>
#include <stdexcept>
//#include <type_traits>

class vect{
public:
	double x;
	double y;
	double z;
	vect() : x(0.0), y(0.0), z(0.0) {	}
	vect(const vect& v) : x(v.x), y(v.y), z(v.z) {	}
	vect(double x1, double y1, double z1) : x(x1), y(y1), z(z1) {	}
	
	friend std::ostream& operator << (std::ostream& s, const vect& v);
	
	double mag() const{
		return sqrt(x*x + y*y + z*z);
	}
	vect cross(const vect& other) const{
		return vect(y*other.z-z*other.y, -(x*other.z-z*other.x), x*other.y-y*other.x);
	}
	
	void rotateX(double angle){
		double newY=y*cos(angle)-z*sin(angle);
		double newZ=y*sin(angle)+z*cos(angle);
		y=newY;
		z=newZ;
	}
	
	void rotateY(double angle){
		double newX=x*cos(angle)+z*sin(angle);
		double newZ=-x*sin(angle)+z*cos(angle);
		x=newX;
		z=newZ;
	}
	
	void rotateZ(double angle){
		double newX=x*cos(angle)-y*sin(angle);
		double newY=x*sin(angle)+y*cos(angle);
		x=newX;
		y=newY;
	}
	
	double r() const{
		return sqrt(x*x + y*y + z*z); //this is the same as mag, but I'm leaving them separate
	}								  //so I don't confuse myself later
	double theta() const{
		if(this->r()==0)
			return 0;
		return acos(z/(this->r()));
	}
	double phi() const{
		double pi=4*atan(1);
		double p=atan2(y,x);
		if(p<0) p+=2*pi;
		if(p>2*pi) p-=2*pi;
		return p;
	}
	
	vect& operator += (const vect& other){
		x += other.x;
		y += other.y;
		z += other.z;
		return(*this);
	}
	
	vect operator + (const vect& other) const{
		vect v1 = *this;
		v1 += other;
		return v1;
	}
	
	vect& operator -= (const vect& other){
		x -= other.x;
		y -= other.y;
		z -= other.z;
		return(*this);
	}
	
	vect operator - (const vect& other) const{
		vect v1 = *this;
		v1 -= other;
		return v1;
	}
	
	double operator * (const vect& other) const{	//dot product
		return (x*other.x + y*other.y + z*other.z);
	}
	
	vect& operator *= (const double d){
		x*=d;
		y*=d;
		z*=d;
		return *this;
	}
	
	vect operator * (const double d) const{
		return vect(*this)*=d;
	}
	
	
	vect& operator /= (const double d){
		x/=d;
		y/=d;
		z/=d;
		return *this;
	}
	
	
	vect operator / (const double d) const{
		return vect(*this)/=d;
	}
	
	bool operator == (const vect& other){
		return (x==other.x && y==other.y && z==other.z);
	}
	
	bool operator != (const vect& other){
		return (!(*this==other));
	}
	
	vect operator - (){
		vect v(-x,-y,-z);
		return v;	
	}
};

vect operator * (double d, const vect& v);
std::ostream& operator << (std::ostream& s, const vect& v);

//a dyad is the outer product of two vectors
//need these for line/cone intersection
class dyadic{
public:
	double contents[3][3];
	dyadic(){
		for(int i=0;i<3;i++){
			for(int j=0;j<3;j++){
				contents[i][j]=0;
			}
		}
	}
	
	//matrix right multiplied by vector
	vect operator * (const vect v){
		vect v2;
		v2.x = v.x*contents[0][0]+v.y*contents[0][1]+v.z*contents[0][2];
		v2.y = v.x*contents[1][0]+v.y*contents[1][1]+v.z*contents[1][2];
		v2.z = v.x*contents[2][0]+v.y*contents[2][1]+v.z*contents[2][2];
		return v2;
	}
};

class line{
public:
	vect position;
	vect direction;
	line() : position(), direction() {	}
	line(const line& l) : position(l.position), direction(l.direction) {	}
	line(vect vpos, vect vdir) : position(vpos), direction(vdir) {	}
	bool operator == (const line& other){
		return (position==other.position && direction==other.direction);
	}
};

class plane{
public:
	vect normal;
	vect position;
	plane() : normal(), position() {	}
	plane(const plane& p) : normal(p.normal), position(p.position) {	}
	plane(vect vpos, vect vnorm) : normal(vnorm), position(vpos) {	}
	bool operator == (const plane& other){
		return (position==other.position && normal==other.normal);
	}
};

class cone{
public:
	vect vertex; //vertex position
	vect axis; //axis direction
	double angle; //opening angle
	cone() : vertex(), axis(), angle(0) {	}
	cone(const cone& c) : vertex(c.vertex), axis(c.axis), angle(c.angle) {	}
	cone(vect vvertex, vect vaxis, double angle1) : vertex(vvertex), axis(vaxis), angle(angle1) {	}
	bool operator == (const cone& other){
		return (vertex==other.vertex && axis==other.axis && angle==other.angle);
	}
};

class sphere{
public:
	double r;
	vect position;
	sphere() : r(0.0), position() {	}
	sphere(const sphere& s) : r(s.r), position(s.position)	{	}
	sphere(vect p1, double r1) : r(r1), position(p1)	{	}
	bool operator == (const sphere& other){
		return (position==other.position && r==other.r);
	}
};

enum objectType {
	PLANE, SPHERE, CONE, ERROR
};

class intersection{
public:
	vect location;
	objectType type;
	uint8_t side;
	intersection() : location(), type(ERROR)	{	}
	intersection(vect v, objectType t) : location(v), type(t)	{	}
	intersection(vect v, objectType t, uint8_t s) : location(v), type(t), side(s)	{	}
};

//So I used to do all the intersection calculations with std::vectors, but it turns out
//that allocating the memory was wasting a lot of time. So instead, this is basically a
//small vector of predefined capacity. It's always allocated on the stack.
//Adding this actually makes the main code 50% faster
template<uint8_t Capacity>
struct shortIntersectionList{
private:
	//static_assert(std::is_default_constructible<intersection>::value,"Intersections are assumed to be default constructible");
	//static_assert(std::is_trivially_destructible<intersection>::value,"Intersections are assumed to be trivially destructible");
	intersection isects[Capacity];
	uint8_t s; //current size
public:
	//Here I duplicate the standard container interface so we can use this like a vector
	using value_type=intersection;	//basically a typedef
	using reference=intersection&;
	using const_reference=const intersection&;
	using pointer = intersection*;
	using const_pointer = const intersection*;
	using iterator = intersection*;
	using const_iterator = const intersection*;
	using difference_type = ptrdiff_t;
	using size_type = uint8_t;
	
	shortIntersectionList() : s(0) {	}
	template<uint8_t otherCapacity>shortIntersectionList(const shortIntersectionList<otherCapacity>& other): s(0) {
		static_assert(Capacity>=otherCapacity, "Cannot construct an intersection list from a potentially longer one.");
		for(const intersection& isect: other)
			isects[s++]=isect;
	}
	template<uint8_t otherCapacity> shortIntersectionList<Capacity>& operator = (const shortIntersectionList<otherCapacity>& other){
		static_assert(Capacity>=otherCapacity, "Cannot assign a potentially longer intersect list to a shorter one.");
		s=0;
		for(const intersection& isect: other)
			isects[s++]=isect;
		return *this;
	}
	bool empty() const{
		return (s==0);
	}
	size_type size() const{
		return s;
	}
	size_type capacity() const{
		return Capacity;
	}
	reference operator [] (size_t i){
		return isects[i];
	}
	const_reference operator [] (size_t i) const {
		return isects[i];
	}
	void push_back(const intersection& i){
		if(s==Capacity)
			throw std::logic_error("push_back called on full shortIntersectionList");
		isects[s++] = i;
	}
	void insert(iterator pos, const intersection& i){
		if(s==Capacity)
			throw std::logic_error("insert called on full shortIntersectionList");
		for(iterator it=end();it!=pos;it--)
			*it=*(it-1);
		*pos=i;
		s++;
	}
	iterator begin(){
		return isects;
	}
	const_iterator begin()const{
		return isects;
	}
	iterator end(){
		return isects+s;
	}
	const_iterator end()const{
		return isects+s;
	}
	reference front(){
		return isects[0];
	}
	const_reference front()const{
		return isects[0];
	}
	reference back(){
		return isects[s-1];
	}
	const_reference back()const{
		return isects[s-1];
	}
};

shortIntersectionList<2> linePlaneIntersect(line l, plane p);

shortIntersectionList<2> lineSphereIntersect(line l, sphere s);

//take outer product of two vectors
dyadic outerProduct(vect v1, vect v2);


//This process is explained well here:
//http://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
shortIntersectionList<2> lineConeIntersect(line l, cone c);

//calculate the path length through a given cell
//a cell is specified by 2 spheres, 2 cones, and 2 planes
//Typically, we expect two intersections, but complications happen
//when a line enters through a plane, then hits a sphere, and exits through the opposite plane
//or exits through the sphere, re-enters through the sphere, and then exits through the
//oppostie plane. These will give three or four intersections, respectively
double lengthThroughCell(line l, plane p1, plane p2, sphere s1, sphere s2, cone c1, cone c2);

const vect origin(0,0,0);