/*
 *  geometry.cpp
 *  Created by Erik on 11/9/14.
 *  
 *  This file defines the geometry involved in the integral calculation,
 *  including vectors, planes, lines, cones, and spheres, and their intersections.
 *  Path length through a cell is included here, possibly temporarily.
 *  
 *  Note that for the spherical coordinates, r is radial component,
 *  theta is the polar angle, and phi is the azimuthal angle.
 *  
 */

#include "geometry.h"
//#include "diskPhysics.h"
const double pi=4*atan(1);

vect operator * (double d, const vect& v){
	return v*d;
}

std::ostream& operator << (std::ostream& s, const vect& v){
	s << v.x << ", " << v.y << ", " << v.z;
	return s;
}

shortIntersectionList<2> linePlaneIntersect(line l, plane p){
	shortIntersectionList<2> results;
	if(l.direction*p.normal == 0){ //line parallel to plane
		return results;
	}
	double d = ((p.position-l.position)*p.normal)/(l.direction*p.normal);
	if (d==0) { //line is inside plane
		return results;
	}
	vect intersect = l.direction*d + l.position;
	intersection i(intersect, PLANE);
	results.push_back(i);
	return results;
}

shortIntersectionList<2> lineSphereIntersect(line l, sphere s){
	double d1, d2;
	vect p1,p2;
	shortIntersectionList<2> results;
	vect temp= l.position-s.position;
	vect ldir= l.direction;
	double root = (ldir*temp)*(ldir*temp)- (ldir*ldir)*(temp*temp - s.r* s.r);
	if(root < 0){ //no intersection
		return results;
	}
	d1=(-(ldir*temp)+sqrt(root))/(ldir*ldir);
	d2=(-(ldir*temp)-sqrt(root))/(ldir*ldir);
	p1 = l.direction*d1 + l.position;
	//std::cout << p1 << std::endl;
	intersection i1(p1,SPHERE);
	//std::cout << p1 << std::endl;
	results.push_back(i1);
	if(d1 != d2){
		p2 = l.direction*d2 + l.position;
		intersection i2(p2,SPHERE);
		results.push_back(i2);
	}
	//std::cout << p2 << std::endl;
	return results;
}

//take outer product of two vectors
dyadic outerProduct(vect v1, vect v2){
	dyadic d;
	d.contents[0][0]=v2.x*v1.x; d.contents[1][0]=v2.x*v1.y; d.contents[2][0]=v2.x*v1.z;
	d.contents[0][1]=v2.y*v1.x; d.contents[1][1]=v2.y*v1.y; d.contents[2][1]=v2.y*v1.z;
	d.contents[0][2]=v2.z*v1.x; d.contents[1][2]=v2.z*v1.y; d.contents[2][2]=v2.z*v1.z;
	return d;
}

//This process is explained well here:
//http://www.geometrictools.com/Documentation/IntersectionLineCone.pdf
shortIntersectionList<2> lineConeIntersect(line l, cone c){
	shortIntersectionList<2> results;
	double c0,c1,c2,d1,d2;
	vect temp;
	vect delta=l.position-c.vertex;
	dyadic M=outerProduct(c.axis,c.axis);
	M.contents[0][0]-=(cos(c.angle)*cos(c.angle));
	M.contents[1][1]-=(cos(c.angle)*cos(c.angle));
	M.contents[2][2]-=(cos(c.angle)*cos(c.angle));
	temp = M*delta;
	c0=delta*temp;
	c1=l.direction*temp;
	temp = M*l.direction;
	c2= l.direction*temp;
	//std::cout<<"c0: "<<c0<<", c1: "<<c1<<", c2: "<<c2<<std::endl;
	
	if(c2 !=0 && (c1*c1-c0*c2) <0) //no intersection
		return results;
	if(c1*c1-c0*c2 == 0){//repeated roots (one intersection)
		//std::cout << "  one intersection" << std::endl;
		d1=-c1/c2;
		//std::cout << d1 << std::endl;
		vect v=l.direction*d1+l.position;
		temp=v-c.vertex;
		if(c.axis*temp>=0){	//only include if it intersects the right half of the cone
			intersection i(v,CONE);
			results.push_back(i);
		}
		return results;
	}
	//two intersections
	//std::cout << "  two intersections" << std::endl;
	d1= (-c1-sqrt(c1*c1-c0*c2))/c2;
	d2= (-c1+sqrt(c1*c1-c0*c2))/c2;
	//std::cout << d1 << "\t" << d2 << std::endl;
	vect v1=l.direction*d1+l.position;
	vect v2=l.direction*d2+l.position;
	temp=v1-c.vertex;
	if(c.axis*temp>=0){	//only include if it intersects the right half of the cone
		intersection i1(v1,CONE);
		results.push_back(i1);
	}
	temp=v2-c.vertex;
	if(c.axis*temp>=0){
		intersection i2(v2,CONE);
		results.push_back(i2);
	}
	return results;
}

//calculate the path length through a given cell
//a cell is specified by 2 spheres, 2 cones, and 2 planes
//Typically, we expect two intersections, but complications happen
//when a line enters through a plane, then hits a sphere, and exits through the opposite plane
//or exits through the sphere, re-enters through the sphere, and then exits through the
//oppostie plane. These will give three or four intersections, respectively
double lengthThroughCell(line l, plane p1, plane p2, sphere s1, sphere s2, cone c1, cone c2){
	double r1,r2,theta1,theta2,phi1,phi2;
	double length;
	vect origin(0,0,0);
	r1=s1.r;
	r2=s2.r;
	theta1=c1.angle;
	theta2=c2.angle;
	phi1=p1.normal.phi();
	phi2=p2.normal.phi();
	
	//all cones are assumed to be convex. That means that the polar angle theta will never be
	//more than pi/2. To find the absolute angle, check if the cone is pointing down, and if it is
	//add pi/2 to the angle. We can do this because the cone must be pointing either straight up or
	//down in order to be part of the coordinate system
	if (c1.axis.z < 0)
		theta1+=(pi/2);
	if (c2.axis.z < 0)
		theta2+=(pi/2);
	
	//check to make sure the cell is well specified
	//All boundaries of the cell must be centered at the origin, for one thing
	//all cone axes must point either straight up or straight down as well
	//each pair must also be distinct
	if(p1.position!=origin || p2.position!=origin || s1.position!=origin || s2.position!=origin || c1.vertex!=origin || c2.vertex!=origin){
		std::cout << "ERROR: cell not well specified. All elements must be centered at the origin" << std::endl;
		exit(1);
	}
	if(r1==r2 || theta1==theta2 || phi1==phi2){
		std::cout << "ERROR: cell not well specified. Cell is too thin." << std::endl;
		exit(1);
	}
	if(!(c1.axis.z==1 || c1.axis.z==-1) && !(c2.axis.z==1 || c2.axis.z==-1)){
		std::cout << "ERROR: cell not well specified. Cones not aligned with Z axis." << std::endl;
		exit(1);
	}
	
	std::vector<intersection> intersects; //store all intersections here.
	
	//first get the shpere intersection points, and check to make
	//sure that they are within the correct theta and phi boundaries
	shortIntersectionList<2> temp; //use this for each stage of intersection checking
	temp=lineSphereIntersect(l, s1);
	if(!temp.empty()){
		//to be in the allowed region, the point must not be greater than both phi values,
		//less than both phi values, greater than both theta values, or less than both
		//theta values.
		for(const auto &isect: temp){
			if(!((isect.location.phi()<phi1&&isect.location.phi()<phi2) || (isect.location.phi()>phi1&&isect.location.phi()>phi2)||
				 (isect.location.theta()<theta1&&isect.location.theta()<theta2) || (isect.location.theta()<theta1&&isect.location.theta()<theta2))){
				vect v(isect.location-l.position);
				if(v*l.direction >=0){//make sure the point is past the start of the line
					intersects.push_back(isect);
					std::cout << "line/sphere1" << std::endl;
				}
			}
		}
	}
	temp=lineSphereIntersect(l, s2);
	if(!temp.empty()){
		for(int i=0;i<temp.size();i++){
			if(!((temp[i].location.phi()<phi1&&temp[i].location.phi()<phi2) || (temp[i].location.phi()>phi1&&temp[i].location.phi()>phi2)||
				 (temp[i].location.theta()<theta1&&temp[i].location.theta()<theta2) || (temp[i].location.theta()<theta1&&temp[i].location.theta()<theta2))){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					std::cout << "line/sphere2" << std::endl;
				}
			}
		}
	}
	
	//next find the intersections between the line and the cones,
	//and check that they're within the r/phi boundaries
	temp=lineConeIntersect(l, c1);
	if(!temp.empty()){
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<r1 && temp[i].location.r()<r2) || (temp[i].location.r()>r1 && temp[i].location.r()>r2) || 
				 (temp[i].location.phi()<phi1 && temp[i].location.phi()<phi2)||(temp[i].location.phi()>phi1&&temp[i].location.phi()>phi2) )){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					std::cout << "line/cone1" << std::endl;
				}
			}
		}
	}
	temp=lineConeIntersect(l, c2);
	
	if(!temp.empty()){
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<r1 && temp[i].location.r()<r2) || (temp[i].location.r()>r1 && temp[i].location.r()>r2) || 
				 (temp[i].location.phi()<phi1 && temp[i].location.phi()<phi2)||(temp[i].location.phi()>phi1&&temp[i].location.phi()>phi2) )){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					std::cout << "line/cone2" << std::endl;
				}
			}
		}
	}
	
	//then find intersections between the line and the planes,
	//and check that they're within r/theta boundaries
	temp=linePlaneIntersect(l, p1);
	if(!temp.empty()){
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<r1 && temp[i].location.r()<r2) || (temp[i].location.r()>r1 && temp[i].location.r()>r2) || 
				 (temp[i].location.theta()<theta1&&temp[i].location.theta()<theta2) || (temp[i].location.theta()>theta1&&temp[i].location.theta()>theta2) )){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					std::cout << "line/plane1" << std::endl;
				}
			}
			
		}
	}
	temp=linePlaneIntersect(l, p2);
	if(!temp.empty()){
		for(int i=0;i<temp.size();i++){
			if(!( (temp[i].location.r()<r1 && temp[i].location.r()<r2) || (temp[i].location.r()>r1 && temp[i].location.r()>r2) || 
				 (temp[i].location.theta()<theta1&&temp[i].location.theta()<theta2) || (temp[i].location.theta()>theta1&&temp[i].location.theta()>theta2) )){
				vect v(temp[i].location-l.position);
				if(v*l.direction >=0){
					intersects.push_back(temp[i]);
					std::cout << "line/plane2" << std::endl;
				}
			}
		}
	}
	
	//now for the special cases:
	length=0; //default case
	//if, by some small chance, there is only one intersection, we just leave it at zero
	//if there are two intersections, we compute the distance between them
	if(intersects.size()==2){
		vect dist(intersects[1].location-intersects[0].location);
		length=dist.mag();
	}
	//if there are three intersections, it means that the line just grazes the inner sphere
	//or one of the cones, and we can ignore the middle intersection
	if(intersects.size()==3){
		vect dist(intersects[2].location-intersects[0].location);
		length=dist.mag();
	}
	//if there are four intersections, the line enters the cell, leaves again, reenters, and
	//then leaves. We only care about the path between the first and second, and third and 
	//fourth intersections
	if(intersects.size()==4){
		vect dist1(intersects[1].location-intersects[0].location);
		vect dist2(intersects[3].location-intersects[2].location);
		length=dist1.mag()+dist2.mag();
	}
	
	return length;
}
