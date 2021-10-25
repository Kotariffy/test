#include <iostream>
#include <cmath>
#include "math.h"
using namespace std;


int main()
{
	setlocale(LC_ALL, "rus");

	int const SIZE = 160;
	double const PI = 3.14159265;
	double
		koeff_d, F, Rv, Rn,
		b, E, sigma, n, alfa_vi, alfa_zapr, h,
		m_rad, J, delta_x, delta_y,
		a1, b1, max_M_krut, min_M_krut,
		N, Q, S_pop,
		eps_max_n, eps_min_v,
		alfa_v, alfa_n, alfa_ni,
		Lvo, Lno, divo, dino,
		s1n, s2n, s3n,
		s1v, s2v, s3v;

	double
		alfa[SIZE],
		gamma[SIZE],
		x[SIZE],
		y[SIZE],
		z[SIZE],
		u[SIZE],
		ro[SIZE],
		Rv_ellipse[SIZE],
		Rn_ellipse[SIZE],
		rad[SIZE],
		M_krut[SIZE],
		Nm[SIZE],
		Qm[SIZE],
		sigma_v[SIZE],
		sigma_n[SIZE],
		sigma_norm[SIZE],
		sigma_sum_v[SIZE],
		sigma_sum_n[SIZE],
		eps_v[SIZE], 
		eps_n[SIZE];

	// ---------------------------------- //
	// ������� ������� ������� ���������� //
	// ---------------------------------- //
	F = 1000;	// ���� F 1000 [��?]
	// �������������� ���������:
	E = 205.94 * 1000000000;	// ������ ��������� ��������� E
	sigma = 1.079 * 1000000000; // ������ ��������� sigma
	n = 4;						// ����������� ������ ��������� n
	// �������������� �������:
	alfa_vi = 60;			    // ��������� ���� �����. ����������� [��.]
	alfa_zapr = 30;			    // ����������� ���� [��.]
	Rn = 0.03;					// ��������� ������ [�]
	Rv = 0.025;					// ���������� ������ [�]
	b = 0.01;					// ������� ������ [�]

	cout << "----------------------------------" << endl;
	cout << "������� ������� ������� ����������" << endl;
	cout << "----------------------------------" << endl;
	cout << "���� F " << F << " [��] " << endl;
	cout << endl << "�������������� ���������: " << endl;
	cout << "������ ��������� ��������� E " << E << " [��] " << endl;
	cout << "������ ��������� sigma " << sigma << " [��] " << endl;
	cout << "�����. ������ ��������� n " << n << endl;
	cout << endl << "�������������� �������: " << endl;
	cout << "���. ���� �����. ����������� alfa_v " << alfa_vi << " [��.] " << endl;
	cout << "����������� ���� alfa_zapr " << alfa_zapr << " [��.] " << endl;
	cout << "���������� ������ Rv " << Rv << " [�] " << endl;
	cout << "��������� ������ Rn " << Rn << " [�] " << endl;
	cout << "������� ������ b " << b << " [�] " << endl;
	cout << "----------------------------------" << endl;
	
	F /= 10;				    //���� � [�]
	
	h = Rn - Rv;			    // ������� ������ h
					
	m_rad = (Rn + Rv) / 2;		//������� ������ ������ m_rad	
	
	J = b * pow(h, 3) / 12;		// ������ ������ ������� ����������� ������� J
	
	delta_x = 0.137 / 2 * (F * pow(m_rad, 3)) / (E * J); // ��������� ������� ������ � ����������� x
	
	delta_y = -0.149 / 2 * (F * pow(m_rad, 3)) / (E * J);

	// ������� ������� �������
	a1 = (m_rad + delta_x);
	b1 = (m_rad + delta_y);

	double temp_i = 0;
	int count = 0;

	for (int i = 0; i < SIZE; i++)
	{

		alfa[i] = temp_i;
		gamma[i] = (PI / 2) - alfa[i];
		temp_i += 0.01;
		count++;
		//cout << "count = " << count << endl;
		//cout << "Size =  "<< i << endl;
		if (temp_i > (PI / 2))
		{
		//	cout << "break! " << endl;
			break;
		}
		
		//cout << i << ": alfa = " << alfa[i] << endl;
		//cout << i << ": gamma = " << gamma[i] << endl;
	}

	// ������� ����������� ������� ������
	S_pop = b * h;

	for (int i = 0; i < count; i++)
	{
		x[i] = a1 * cos(alfa[i]); // ���������� �� X ������� ����� �������
		y[i] = b1 * sin(alfa[i]); // ���������� �� Y ������� ����� �������

		z[i] = sin(gamma[i]);
		u[i] = cos(gamma[i]);

		// ������� ������ �������� ������� 
		ro[i] = pow(a1, 2) * pow(b1, 2) * pow((pow(x[i], 2) / pow(a1, 4)) + (pow(y[i], 2) / pow(b1, 4)), 1.5);

		
		

		// ������� ����� �������� ����������� �������
		Rn_ellipse[i] = ro[i] + h / 2;
	    // ������� ����� ���������� ����������� �������
		Rv_ellipse[i] = ro[i] - h / 2;
		// ������ �������� ������������ ����
		rad[i] = h / log(Rn_ellipse[i] / Rv_ellipse[i]);

		

		// �������� ������ ���. ������� ����� �� ��. ����� ������
		M_krut[i] = F * m_rad * (0.3183 - z[i] / 2);
		// ������������� ���� ���. ������� ����� �� ��. ����� ������
		Nm[i] = -1 / 2 * F * z[i];
		// ���� ������ � ������� �����
		Qm[i] = -1 / 2 * F * u[i];

		// ����������, ����������� �������� M � ������ ����������� ������� ������
		// �� �����. �����������:
		sigma_v[i] = (M_krut[i] * (Rv_ellipse[i] - rad[i])) / (S_pop * Rv_ellipse[i] * (ro[i] - rad[i]));
		// �� �������� �����������:
		sigma_n[i] = (M_krut[i] * (Rn_ellipse[i] - rad[i])) / (S_pop * Rn_ellipse[i] * (ro[i] - rad[i]));

		// ������ ����������� ������� sigma[i] < sigma, ����� break

		// ����������, ����������� ����� ���������� �� ���� ������ ����������� ������� ������
		sigma_norm[i] = -(1/2) * (F * z[i] / S_pop);

		// ��������� ���������� �� ���������� ������������ ����������� ������� ������
		sigma_sum_v[i] = sigma_v[i] + sigma_norm[i];
		// ��������� ���������� �� �������� ������������ ����������� ������� ������
		sigma_sum_n[i] = sigma_n[i] + sigma_norm[i];

		if ((sigma_sum_v[i] >= n * sigma) || (sigma_sum_n[i] >= n * sigma))
		{
			cout << "�������� ������ ����������. ������ ������������ ���. ������! ";
			continue;
		}

		//���������� ���������� ����������� ������
		eps_v[i] = sigma_sum_v [i] / E;
		//���������� �������� ����������� ������
		eps_n[i] = sigma_sum_n[i] / E;

	//  cout << "(" << alfa[i] * (180 / PI) <<";"<< eps_n[i] <<") " << endl;
	//	cout << "(" << alfa[i] * (180 / PI) << ";" << eps_v[i] << ") " << endl;

		eps_max_n = eps_n[0];
		eps_min_v = eps_v[count - 1];

		max_M_krut = M_krut[0];
		min_M_krut = M_krut[count - 1];
		//!!!!!!! ����� ����� ��������� ������� eps_n(alfa[i]) � eps_v(alfa[i]) !!!!!
		
		
	}

	alfa_vi = alfa_vi * (PI / 180); // ���������� ���� [���.]
	
	alfa_v = (PI / 2) - alfa_vi;				// ���� ������� ���������� ����������� ��� ������������� [���.]
	alfa_n = (alfa_v * Rv) / Rn;				// ���� ������� ������� ����������� ��� ������������� [���.]
	alfa_ni = (PI / 2 - alfa_n) * (180 / PI);	// �������� ����, ��� ������� Ln = Lv [��.]
	Lvo = Rv * alfa_v;							// ����� �� ���������� ����������� [�]
	Lno = Rn * alfa_n;							// ����� �� �������� ����������� [�]

	/* !!!!!! �������� � �������� eps_n(alfa[i]) � eps_v(alfa[i]) ������������ ����� (�����., ����., � ����. ����) 
	for (int i = -round(F); i < round(eps_n[count-1]*1000000); i++)
	{
		alfa_vi * (180 / PI) �� i
		alfa_ni �� i
		90 - alfa_zapr �� i
	}

	*/
	//--------
	// ������ ��������� �����. � �����. ������������ ������ ���������� �� ������ ��������.
	//--------
	/*
	s1n = eps_n[0] + eps_n[count - 1];
	s2n = 0; 
	s3n = 0;

	int k = 0;
	while (k <= count - 3)
	{
		k += 2;
		s2n += eps_n[k];
	}
	
	
	k = 1;
	while (k <= count - 2)
	{
		s3n += eps_n[k];
		k += 2;
	}

	// ��������� ����� �������� ����������� ������. 
	double Integral_n = (alfa[count - 1] - alfa[1]) * (s1n + 2 + s2n + 4 * s3n) / (6 * 79);
	*/
	s1n = eps_n[0] + eps_n[count - 1];
	s2n = 0;
	s3n = 0;

	int k = 0;
	while (k <= count - 3)
	{
		k += 2;
		s2n += eps_n[k];
	}


	k = 1;
	while (k <= count - 2)
	{
		s3n += eps_n[k];
		k += 2;
	}

	// ��������� ����� �������� ����������� ������. 
	double Integral_n = 2 * (alfa[count - 1] - alfa[1]) * (s1n + 2 * s2n + 4 * s3n) / (6 * 79);

	s1v = eps_v[0] + eps_v[count - 1];
	s2v = 0; 
	s3v = 0;
	
	k = 0;
	while (k <= count - 3)
	{
		k += 2;
		s2v += eps_v[k];
	}

	k = 1;
	while (k <= count - 2)
	{
		s3v += eps_v[k];
		k += 2;
	}

	// ��������� ����� ���������� ����������� ������. �������� �� ������ ��������. 

	double Integral_v =  2 * (alfa[count - 1] - alfa[1]) * (s1v + 2 * s2v + 4 * s3v) / (6 * 79); 
	
	double LvF = 2 * (Lvo + Integral_v);		// ����� ���������� ����������� ������ ��� �������� ���� 
	double LnF = 2 * (Lno + Integral_n);		// ����� �������� ����������� ������ ��� �������� ���� 

	alfa_vi *= 180 / PI;
	

	cout << endl;
	cout << "----------------------------------" << endl;
	cout << "���������� ������������� ������" << endl;
	cout << "----------------------------------" << endl;
	cout << " ������ ������ " << h << " [�] " << endl;
	cout << " ������� ������ ������ " << m_rad << " [�] " << endl;
	cout << " ������ ������ ������� " << J << " [�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ��������� �������� ������ � ����������� x " << delta_x << " [�] " << endl;
	cout << " ��������� ������� ������ � ����������� y " << delta_y << " [�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ������� ������� ������� a1 " << a1 << " [�] " << endl;
	cout << " ������� ������� ������� b1 " << b1 << " [�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ������� ����������� ������� " << S_pop << " [�^2] " << endl;
	cout << " ������������ �������� ������ " << max_M_krut << " [�*�] " << endl;
	cout << " ����������� �������� ������ " << min_M_krut <<" [�*�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ���� ��� ������������� �� �����. ����������� " << alfa_n  <<" [���] " << endl;
	cout << " ���� ��� ������������� �� �����. ����������� " << alfa_v  <<" [���] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ����� ������� �� �����. ����������� " << Lno << " [�] " << endl;
	cout << " ����� ������� �� �����. ����������� " << Lvo << " [�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ��������� ����� �� �����. ����������� " << Integral_n << " [�] " << endl;
	cout << " ��������� ����� �� �����. ����������� " << Integral_v << " [�] " << endl;
	cout << "----------------------------------" << endl;
	cout << " ����� �������� ����������� ������ ��� �������� ���� " << LnF << " [�] " << endl;
	cout << " ����� ���������� ����������� ������ ��� �������� ���� " <<LvF << " [�] " << endl;
	cout << "----------------------------------" << endl;

	cout << "----------------------------------" << endl;
	cout << " ���������� ����. ��������������" << endl;
	cout << "----------------------------------" << endl;
	cout << "     ���������� ����������� " << endl;
	cout << "----------------------------------" << endl;
	cout << "  ���� [��.]\t\t���������� 10^5 " << endl;

	

	for (int i = 0; i < count; i+=10)
	{
		cout << fixed;
		cout.precision(4);
		cout << "  " << alfa[i] * 180 / PI << "\t\t";
		cout << fixed;
		cout.precision(4);
		cout << eps_v[i] * 100000 << " " << endl;
	}

	
	cout << "----------------------------------" << endl;
	cout << "       �������� ����������� " << endl;
	cout << "----------------------------------" << endl;
	cout << "  ���� [��.]\t\t���������� 10^5 " << endl;

	for (int i = 0; i < count; i+=10)
	{
		cout << fixed;
		cout.precision(4);
		cout << "  " << alfa[i] * 180 / PI << "\t\t";
		cout << fixed;
		cout.precision(4); 
		cout<< eps_n[i] * 100000 << " " << endl;
	}
	
	return 0;
}
	/*
	int i; // �������

	double Integral;
	double a = 0.0, b = alfa_n; // ����� ������� �������������� = 0,436775
	double h = 0.01;// ����� ��� ��������������
	double n; // ����� ����� ��������� n

	n = (b - a) / h;

	// ��������� �������� �� ������� ��������
	
	Integral = h * (f(a) + f(b)) / 6.0;

	for (i = 1; i <= n; i++)

		Integral = Integral + 4.0 / 6.0 * h * f(a + h * (i - 0.5));

	for (i = 1; i <= n - 1; i++)
		Integral = Integral + 2.0 / 6.0 * h * f(a + h * i);
	*/

