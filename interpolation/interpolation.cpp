#include "Interpolation.h"

Interpolation::Interpolation(std::vector <double> Ox, std::vector <double> Oy) {
	m_Ox = Ox;
	m_Oy = Oy;
	m_k = new double[m_Ox.size()];
	m_b = new double[m_Ox.size()];
	m_C = new double[m_Ox.size()];
	m_q = new double[m_Ox.size()];
	m_L = new double[m_Ox.size()];
	m_matrix = new double* [m_Ox.size()];
	for (int i = 0; i < m_Ox.size(); i++) {
		m_matrix[i] = new double[m_Ox.size() + 1];
	}
	m_y = new double[m_Ox.size()];
	outputPoints();
}

void Interpolation::canonical() {
	creatingMatrix();
	gauss();
	outputCanonical();
}

void Interpolation::outputPoints() {
	std::cout << "Координаты по оси Ox: ";
	for (int i = 0; i < m_Ox.size(); i++) {
		std::cout << std::setw(8) << m_Ox[i];
	}
	std::cout << std::endl;

	std::cout << "Координаты по оси Oy: ";
	for (int i = 0; i < m_Oy.size(); i++) {
		std::cout << std::setw(8) << m_Oy[i];
	}

	std::cout << std::endl << std::endl;
}

void Interpolation::formulaLagrange() {
	for (int i = 0; i < m_Ox.size() / 2 + 1; i++) {
		m_q[i] = 1;
		for (int j = 0; j < m_Ox.size() / 2 + 1; j++) {
			if (i != j) {
				m_q[i] *= (m_Ox[i * 2] - m_Ox[j * 2]);
			}
		}
		m_q[i] = m_Oy[i * 2] / m_q[i];
	}
	
	for (int i = 0; i < m_Ox.size() / 2 + 1; i++) {
		m_q[i] = round(m_q[i] * 1000) / 1000;
	}
	
	for (int h = 0; h < m_Ox.size() / 2; h++) {
		m_Oy[h * 2 + 1] = 0;
		for (int i = 0; i < m_Ox.size() / 2 + 1; i++) {
			m_L[i] = 1;
			for (int j = 0; j < m_Ox.size() / 2 + 1; j++) {
				if (i != j) {
					m_L[i] *= (m_Ox[h * 2 + 1] - m_Ox[j * 2]);
				}
			}
			m_L[i] *= m_q[i];
			m_Oy[h * 2 + 1] += m_L[i];
		}
	}
	
	for (int i = 0; i < m_Ox.size() / 2; i++) {
		m_Oy[i * 2 + 1] = round(m_Oy[i * 2 + 1] * 1000) / 1000;
	}

	outputFormulaLagrange();
}

void Interpolation::linear() {
	for (int i = 0; i < (m_Ox.size() + 1) / 2 - 1; i++) {
		m_k[i] = (m_Oy[2 * i + 2] - m_Oy[2 * i]) / (m_Ox[2 * i + 2] - m_Ox[2 * i]);
		m_b[i] = m_Oy[2 * i] - m_Ox[2 * i] * m_k[i];
	}
	outputLinear();
}

void Interpolation::outputLinear() {
	std::cout << "Линейная интерполяция. Задается кусочно." << "\n";
	for (int i = 0; i < (m_Ox.size() + 1) / 2 - 1; i++) {
		std::cout << "y = " << std::setw(3) << m_k[i] << " * x ";
		if (m_b < 0) {
			std::cout << "- ";
		}
		else {
			std::cout << "+ ";
		}
		std::cout << std::setw(3) << abs(m_b[i]);
		std::cout << " , при x [" << std::setw(4) << m_Ox[2 * i] << " ; " <<
			std::setw(4) << m_Ox[2 * i + 2] << "] ";
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Interpolation::outputFormulaLagrange() {
	std::cout << "формула Лагранжа. Координаты точек." << "\n";
	outputPoints();
}

void Interpolation::creatingMatrix() {
	for (int i = 0; i < m_Ox.size(); i++) {
		for (int j = 0; j < m_Ox.size(); j++) {
			m_matrix[i][j] = pow(m_Ox[i], j);
		}
	}

	for (int i = 0; i < m_Ox.size(); i++) {
		m_matrix[i][m_Ox.size()] = m_Oy[i];
	}
}

void Interpolation::gauss() {
	double k;
	// Приведение матрицы к трапециевидному виду
	for (int h = 0; h < m_Ox.size() - 1; h++) {
		for (int i = 1 + h; i < m_Ox.size(); i++) {
			k = m_matrix[h][h] / m_matrix[i][h];
			for (int j = h; j < m_Ox.size() + 1; j++) {
				m_matrix[i][j] = m_matrix[i][j] * k - m_matrix[h][j];
			}
		}
	}
	// Решение системы получившихся уравнений
	for (int i = m_Ox.size() - 1; i >= 0; i--) {
		m_C[i] = 0;
		for (int j = m_Ox.size() - 1; j > i; j--)
			m_C[i] += m_matrix[i][j] * m_C[j];
		m_C[i] = (m_matrix[i][m_Ox.size()] - m_C[i]) / m_matrix[i][i];
	}
	// Округление результатов
	for (int i = 0; i < m_Ox.size(); i++) {
		m_C[i] = round(m_C[i] * 1000) / 1000;
	}
}

void Interpolation::outputCanonical() {
	std::cout << "Интерполяция каноническим полиномом." << "\n";
	std::cout << "ф(x) = " << m_C[0];
	for (int i = 1; i < m_Ox.size(); i++) {
		if (m_C[i] < 0)
			std::cout << " - ";
		else
			std::cout << " + ";
		std::cout << abs(m_C[i]) << " * x^" << i;
	}
	std::cout << "\n" << "\n";
}

void Interpolation::formulaNewton() {
	double** matrix = new double* [m_Ox.size()];
	for (int i = 0; i < m_Ox.size(); i++) {
		matrix[i] = new double[m_Ox.size()];
	}
	
	for (int i = 0; i < m_Ox.size(); i++) {
		matrix[i][0] = m_Oy[i];
	}
	
	for (int j = 1; j < m_Ox.size(); j++) {
		for (int i = 0; i < m_Ox.size(); i++) {
			if (m_Ox.size() - j > i) {
				matrix[i][j] = (matrix[i + 1][j - 1] - matrix[i][j - 1]) / (m_Ox[j + i] - m_Ox[i]);
			}
		}
	}

	for (int i = 0; i < m_Ox.size(); i++) {
		m_y[i] = round(matrix[0][i] * 1000) / 1000;
	}

	for (int i = 0; i < m_Ox.size(); i++) {
		delete[] matrix[i];
	}
	delete[] matrix;

	outputFormulaNewton();
}

void Interpolation::outputFormulaNewton() {
	std::cout << "Интерполяционные полином Ньютона." << "\n";
	std::cout << "P(x) = " << m_y[0];
	for (int i = 1; i < m_Ox.size(); i++) {
		if (m_y[i] < 0)
			std::cout << " - ";
		else
			std::cout << " + ";
		std::cout << abs(m_y[i]);
		for (int j = 0; j < i; j++) {
			std::cout << " * (x";
			if (m_Ox[j] < 0)
				std::cout << " + ";
			else
				std::cout << " - ";
			std::cout << abs(m_Ox[j]) << ")";
		}
	}
	std::cout << "\n";
}

Interpolation::~Interpolation() {
	for (int i = 0; i < m_Ox.size(); i++) {
		delete[] m_matrix[i];
	}
	delete[] m_matrix;
	delete[] m_b;
	delete[] m_C;
	delete[] m_k;
	delete[] m_q;
	delete[] m_y;
	delete[] m_L;
}