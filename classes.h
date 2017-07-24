#ifndef classes_h
#define classes_h

#include<vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include<string>
#include "macros.h"

using namespace std;

void tic();
void toc();
#include "globalvar.h"
#define PI 3.14159265
using namespace std;

//class Gears
//{
//	Gears(double wd, int sections)
//	{
//		depth = wd;
//		m = sections;
//		delX = depth / m;
//	}
//	Gears() {}
//public:
//	double depth; //gear whole depth
//	int m; //# of sections of TSV
//	double delX; //depth of each TSV slice
//};


class Mtable
{
public:
	Mtable(const string& TXTfile1, const string& TXTfile2);
	Mtable(){}
	
	
	vector<vector<double> > matTable1, matTable2; // for negative, positive gauge pressure, repsectively
	double pstep1, pstep2;
	double RefRho; // [kg/m**3] reference density
	
	double BatP(double p);
	double RhoatP(double p);
	double NuatP(double p);

private:
};

class Gtable
{
public:
	Gtable(){}
	vector<double> GeomAngle;
	vector<vector<double> > GeomTable;

	int length;
	int NumberOfFeatures;
	double angleStep()
	{
		return 360.00 / (length - 1);
	}

	vector<vector<double> > Vderiv;

private:
};


class single_Chamber
{
	// Private
protected:
	int neq_;
	int gindex_;

public:
	typeName(single_Chamber);
	single_Chamber();
	virtual ~single_Chamber();
	// Construct from Index and Dead Volume
	single_Chamber(int,double,Gtable*,Mtable*);

	virtual void setTime(double t);
	virtual void initialize(double* y,int& i);
	virtual void update(double* ,double);
	virtual void getydot(double* );
	
	virtual inline double BatP() const
	{
		return materialTable->BatP(p);
	}
	inline double RhoatP() const
	{
		return materialTable->RhoatP(p);
	}
	inline double RhoCalculated() const
	{
		return mass/volume;
	}
	inline double NuatP() const
	{
		return materialTable->NuatP(p);
	}
	
	double RefDensity() const
	{
		return materialTable->RefRho;
	}

	inline int neq() const
	{
		return neq_;
	}
	
	virtual double getF(); // evaluete RHS of pressure build-up
	//inline double getP() const
	//{
	//	return last_p;
	//}// get the pressure of last time step
	//inline double getdPdt() const
	//{
	//	return last_dpdt;
	//}// get the diffierential pressure of last time step

	double getMomentumfluxforconnection()
	{
		double velocity = (vleft + vright) / 2;
		double vgrad = (vright - vleft) / delX;
		double mu = RhoatP()*NuatP();
//		cout << endl << p << "  " << velocity*velocity*RhoatP() << "  " << 4.0 / 3.0*mu*vgrad << endl;
		return (p + velocity*velocity*RhoatP() - 4.0 / 3.0*mu*vgrad)*Area;
		//return (p + velocity*velocity*RhoatP() )*Area;
	}
	double time;
	double initialPosition; // [deg]
	double angularV;
	int index; // volume drive or slave, 1 for drive, 2 for slave; also the volume index in geometric table.
	double p;
	//double dp;
	double angularPosition;
	//vector<double> output_angularPosition;
	/*double last_time;
	double last_p;
	double last_dpdt;*/
	
	double mass;
	double mflux;

	//vector<double> recordTime;
	//vector<double> recordP;
	//vector<double> recordDpdt;

	double volume;
	double volumeDeriv;
	double bulkModulus;
	double deadVolume; // non-zero for inlet/outlet volume.
	vector<double> dm; // all the mass flux at specific time

	Gtable *Geometry;
	Mtable *materialTable;
	//double pstep; // p step of material look-up table

	double Area; //axial area of TSV
	double vleft;
	double vright;//axial velocity in TSV
	double last_p,last_time;
	
private:

};


class single_linkage
{
	protected:
		int inlet_index;
		int outlet_index;
public:

	typeName(single_linkage);

	single_linkage(single_Chamber& inletC, single_Chamber& outletC)
		:inlet(&inletC),
		outlet(&outletC),
		inlet_index(0),
		outlet_index(0)
	{
		inlet->dm.push_back(0.00);
		outlet->dm.push_back(0.00);
		inlet_index = inlet->dm.size()-1;
		outlet_index = outlet->dm.size()-1;
//		inlet = &inletC;
//		outlet = &outletC;
	}
	// Definition: inlet is the chamber from which linkage reads angular position from:

	single_linkage(){}

	//virtual void setTime(double t){ std::cout << "This is wrong..." << std::endl; cin.get();}
	virtual void setTime(double) = 0;
	virtual void getFlux()=0;
	//virtual void getFlux(){ std::cout << "This is wrong..." << std::endl; cin.get();}
	virtual void update(double);
	double BatP();
	double RhoatP();
	double NuatP();
	inline double RefDensity() const
	{
		return materialTable->RefRho;
	}

	single_Chamber* inlet;
	single_Chamber* outlet;
	double flux;
	double time;
	

	double p; // the average pressure btw inlet and outlet
	//vector<double> output_mflux;
	double last_time;
	vector<double> recordTime;
	//vector<double> output_orificeA;
	//vector<double> output_orificeHD;

	Gtable *Geometry;
	Mtable *materialTable;
	//double pstep; // p step of material look-up table
	int indexA,indexHD;

private:
	//single_Chamber Null;
};

class single_orifice : public single_linkage
{
public:
	typeName(single_orifice);
	single_orifice():single_linkage(){}
	single_orifice(single_Chamber& inletC, single_Chamber& outletC):single_linkage(inletC,outletC)
	{
		inlet = &inletC;
		outlet = &outletC;
	}

	void getFlux() override;
	
	
	void setTime(double t);
	int indexA; // index in geometric table
	int indexHD;
	double orificeA;
	double orificeHD;
	int side; //-1 for left, 1 for right applicable only to HG/LG
	
private:
	double lambda_crit = 1000.0;  // critical flow number (larminar -> turbulence)
	double cqmax = 0.7; // maximum discharge coefficient
	double cq; // discharge coefficient
	double lambda; // flow number
};

class single_port : public single_linkage
{
public:
	typeName(single_port);
	single_port(){}
	single_port(double pressure, single_Chamber& outletC, double d, double cq_max);
	void getFlux();
	void setPortA(double);
	void setTime(double t)
	{
		time = t;
	}

	double portA; // m^2
	double portHD; // m
	double portP;
	
	double lambda_crit = 1000.0;  // critical flow number (larminar -> turbulence)
	double cqmax; // maximum discharge coefficient
	double cq; // discharge coefficient
	double lambda; // flow number
private:
};




class single_drainLeakage : public single_linkage
{
public:
	typeName(single_drainLeakage);
	single_drainLeakage(){}
	single_drainLeakage(double pressure, single_Chamber& outletC, double, double, double, int);
	void getFlux();
	inline void setTime(double t)
	{
		time = t;
	}
	double portP;
	double L;
	double h;
	double b;
	

private:
};


class single_LateralLeakage : public single_linkage
{
public:
	single_LateralLeakage() {}
	single_LateralLeakage(single_Chamber& inletC, single_Chamber& outletC, double, double, double, double, int);
	void getFlux();
	inline void setTime(double t)
	{
		time = t;
	}
	double u;
	double L;
	double h;
	double b;


private:
};

class single_RadialLeakage : public single_linkage
{
public:
	single_RadialLeakage() {}
	single_RadialLeakage(single_Chamber& inletC, single_Chamber& outletC, double, double, double, double, double, double, double, int);
	void getFlux();
	inline void setTime(double t)
	{
		time = t;
	}
	double u;
	double L;
	double h;
	double b;
	double casing_start;
	double casing_end;


private:
};






class Chambers
{
public:
	typeName(Chambers);
	Chambers(int n, int m=1);
	Chambers()
	{}

	int number_of_chambers;
	int number_of_cells;
	vector<vector<single_Chamber> > single_cell; 
	inline int size() const
	{
		return number_of_chambers;
	}

private:
	//double angle;
};


class Orifices
{
public:
	typeName(Orifices);
	//Orifices(int n)
	//{
	//	number_of_orifices = n;
	//	single_array.resize(n);
	//}
	//Orifices(int n);
	Orifices(int n, int m=1);

	Orifices()
	{}

	//int number_of_orifices;
	int number_of_cells;
	//vector<single_orifice> single_array;
	vector<vector<single_orifice> >single_cell;

private:
};

class Ports
{
public:
//	typeName(Ports);
//	Ports(int n)
//	{
//		single_array.resize(n);
//	}
	Ports(int n, int m);
	Ports()
	{}

	//vector<single_port> single_array;
	vector<vector<single_port> >single_cell;

private:
};

class DrainLeakage
{
public:
	DrainLeakage(int n, int m=1);
	DrainLeakage()
	{}
	
	vector<vector<single_drainLeakage> > single_cell;
	
private:
};

class LateralLeakage
{
public:
	LateralLeakage(int n, int m=1);

	LateralLeakage()
	{}

	vector<vector<single_LateralLeakage> > single_cell;

private:
};

class RadialLeakage
{
public:
	RadialLeakage(int n, int m);

	RadialLeakage()
	{}

	vector<vector<single_RadialLeakage> > single_cell;

private:
};





class single_cell_connection: public single_linkage
{
public:
	single_cell_connection() {}
	single_cell_connection(single_Chamber& inletC, single_Chamber& outletC) :single_linkage(inletC, outletC)
	{
		//velocity = 0;
		momentum = 0;
		gindex =-1;

	}
	void getFlux() override  //gets flux from momentum
	{
		flux = momentum / delX;  // mdot = mv/delX   ... need to define delX properly.
		(*inlet).vright = flux / ((*inlet).RhoatP()*(*inlet).Area);
		(*outlet).vleft = flux / ((*outlet).RhoatP()*(*outlet).Area);  
			// update the dm term at inlet/outlet chamber
		inlet->dm[inlet_index] = -flux;
		outlet->dm[outlet_index] = flux;
		return;
	}

//	double getMomentumFlux() 
//	{
//		return (*inlet).getMomentumfluxforconnection() - (*outlet).getMomentumfluxforconnection();
//	}

	void initialize(double* y, int& i)
	{
		gindex=i;
		y[i+1] = momentum;
		i += neq_;
	}
	void update(double* y,double t)
	{
		if (gindex==-1)
		{
			cout << "single_cell_connection not initialized by multistep" << endl;
			exit(0);
		}
		momentum = y[gindex];
		this->getFlux();
	
	}
	void update(double t)
	{
		return;
	}

	void setTime(double t)
	{
		return;
	}
	double getwallshearforce()
	{
		if (momentum == 0)
			return 0;
		p = 0.5 * ((*inlet).p + (*outlet).p);
		double rho = RhoatP();
		double nu = NuatP();
		double Area = 0.5 * ((*inlet).Area + (*outlet).Area);
		double vel = momentum / (rho*Area*delX);
		double diameter = sqrt(4.0 * Area / PI);
		double Re = vel*diameter / nu;
		double f;
		double e = 5e-6; //roughness = 5 microns
		if (Re < 2000)
			f = 64.0 / Re;
		else
			f = pow(1.0 / (-1.8*log(pow(e / diameter / 3.7, 1.11) + 6.9 / Re)), 2);
		double l = 4 * Area / diameter;//wettedperimeter
		double shearforce = f*l*delX * 1 / 8.0*rho*vel*vel;
		if (vel < 0)
			shearforce = -shearforce;
		return shearforce;
	}
	void getydot(double* ydot)
	{
		ydot[gindex] = (*inlet).getMomentumfluxforconnection() - (*outlet).getMomentumfluxforconnection() -getwallshearforce();
	}

	double neq()
	{
		return neq_;
	}
	

	int gindex,neq_=1;
//	double velocity;
	double momentum;  //direction is from inlet identified chamber to outlet identified chamber
	
private:
	
};

class cell_connections
{
public:
	cell_connections(int n, int m)
	{
		
		single_cell.resize(n);
		if(m>1)
			for (int i = 0; i < n; i++)
				single_cell[i].resize(m-1);
	}

	cell_connections()
	{}

	
	vector<vector<single_cell_connection> > single_cell;

private:
};


#endif
