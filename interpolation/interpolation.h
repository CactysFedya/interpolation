#include <iostream>
#include <vector>
#include <iomanip>

#ifndef INTERPOLATION_H    
#define INTERPOLATION_H

class Interpolation {
private:
	std::vector <double> m_Ox;
	std::vector <double> m_Oy;
	double** m_matrix;
	double* m_k;
	double* m_b;
	double* m_C;
	double* m_q;
	double* m_L;
	double* m_y;

public:
	Interpolation(std::vector <double>, std::vector <double>);

	void linear();

	void canonical();

	void formulaLagrange();

	void formulaNewton();

	~Interpolation();

private:
	void outputPoints();

	void outputCanonical();

	void outputFormulaLagrange();

	void outputFormulaNewton();

	void gauss();

	void outputLinear();

	void creatingMatrix();
};

#endif INTERPOLATION_H