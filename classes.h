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
	virtual void outerUpdate();
	
	virtual inline double BatP() const
	{
		return materialTable->BatP(p);
	}
	double RhoatP() const
	{
		return materialTable->RhoatP(p);
	}
	inline double RhoCalculated() const
	{
		return mass/volume;
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
	inline double getP() const
	{
		return last_p;
	}// get the pressure of last time step
	inline double getdPdt() const
	{
		return last_dpdt;
	}// get the diffierential pressure of last time step

	double time;
	double initialPosition; // [deg]
	double angularV;
	int index; // volume drive or slave, 1 for drive, 2 for slave; also the volume index in geometric table.
	double p;
	double dp;
	double angularPosition;
	//vector<double> output_angularPosition;
	double last_time;
	double last_p;
	double last_dpdt;
	
	double mass;
	double mflux;

	vector<double> recordTime;
	vector<double> recordP;
	vector<double> recordDpdt;

	double volume;
	double volumeDeriv;
	double bulkModulus;
	double deadVolume; // non-zero for inlet/outlet volume.

	vector<double> dm; // all the mass flux at specific time

	Gtable *Geometry;
	Mtable *materialTable;
	//double pstep; // p step of material look-up table

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
	int indexA; // index in geometric table
	int indexHD;
	double orificeA;
	double orificeHD;

	double p; // the average pressure btw inlet and outlet
	//vector<double> output_mflux;
	double last_time;
	vector<double> recordTime;
	//vector<double> output_orificeA;
	//vector<double> output_orificeHD;

	Gtable *Geometry;
	Mtable *materialTable;
	//double pstep; // p step of material look-up table

private:
	//single_Chamber Null;
};

class single_hyd_joint: public single_linkage
{
public:
	typeName(single_hyd_joint);

	single_hyd_joint(){}
	
	single_hyd_joint(single_Chamber& inletC, single_Chamber& outletC, double area, double hydD, double length)
	:single_linkage(inletC,outletC)
	,A(area)
	,HD(hydD)
	,L(length)
	{
	}

	void getFlux()
	{
		p = 0.5 * ((*inlet).p + (*outlet).p);
		
		int sign;
		if ((*inlet).p - (*outlet).p  > 0)
			sign = 1;
		else
			sign = -1;
		
		
		flux = sign * A * fabs((*inlet).p - (*outlet).p) * HD * HD/32.0/L/NuatP();
		
		// update the dm term at outlet chamber only... (inlet is pressure boundary condition)
//		(*inlet).dm.push_back(-flux);
//		(*outlet).dm.push_back(flux);
		inlet->dm[inlet_index] = -flux;
		outlet->dm[outlet_index] = flux;
		return;
	}

	void setTime(double t)
	{
		time = t;
	}

	void update(double t)
	{
		this->setTime(t);
		this->getFlux();
	}

private:
	double A;
	double HD;
	double L;
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
		
//		inlet->dm.push_back(0.00);
//		outlet->dm.push_back(0.00);
//		inlet_index=inlet->dm.size()-1;
//		outlet_index=outlet->dm.size()-1;
	}

	void getFlux() override;
	
	
	void setTime(double t);
	
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
	single_drainLeakage(double pressure, single_Chamber& outletC, double La, double Lb, double Lc);
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


class single_Leakage : public single_linkage
{
public:
	typeName(single_Leakage);
	single_Leakage(){}
	single_Leakage(single_Chamber& inletC, single_Chamber& outletC, double La, double Lb, double Lc, double inputU);
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





class Chambers
{
public:
	typeName(Chambers);
	Chambers(int n);
	//Rituraj
	Chambers(int n, int m);
	Chambers()
	{}

	int number_of_chambers;
	int number_of_cells;
	vector<single_Chamber> single_array;
	vector<vector<single_Chamber> > single_cell;  //Rituraj
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
	Orifices(int n)
	{
		number_of_orifices = n;
		single_array.resize(n);
	}

	Orifices()
	{}

	int number_of_orifices;
	vector<single_orifice> single_array;

private:
};

class Ports
{
public:
	typeName(Ports);
	Ports(int n)
	{
		single_array.resize(n);
	}

	Ports()
	{}

	vector<single_port> single_array;

private:
};

class DrainLeakage
{
public:
	typeName(DrainLeakage);
	DrainLeakage(int n)
	{
		single_array.resize(n);
	}
	
	DrainLeakage()
	{}
	
	vector<single_drainLeakage> single_array;
	
private:
};

class Leakage
{
public:
	typeName(Leakage);
	Leakage(int n)
	{
		single_array.resize(n);
	}
	
	Leakage()
	{}
	
	vector<single_Leakage> single_array;
	
private:
};


class RHS
{
public:
	RHS(){}
	void operator () (double t, vector<double> &tempP, vector<double> &k, vector<single_Chamber*> allChambers, vector<single_linkage*> allLinkages);
	void operator () (double t, double* P, double* Dpdt, vector<single_Chamber*> allChambers, vector<single_linkage*> allLinkages);
	
private:
};









#endif
