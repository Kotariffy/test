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
	// Типовой вариант входных параметров //
	// ---------------------------------- //
	F = 1000;	// Сила F 1000 [кг?]
	// Характеристики материала:
	E = 205.94 * 1000000000;	// Модуль упругости материала E
	sigma = 1.079 * 1000000000; // Предел прочности sigma
	n = 4;						// Коэффициент запаса прочности n
	// Геометрические размеры:
	alfa_vi = 60;			    // Начальный угол внутр. поверхности [гр.]
	alfa_zapr = 30;			    // Запрещенный угол [гр.]
	Rn = 0.03;					// Наруженый радиус [м]
	Rv = 0.025;					// Внутренний радиус [м]
	b = 0.01;					// Толщина кольца [м]

	cout << "----------------------------------" << endl;
	cout << "ТИПОВОЙ ВАРИАНТ ВХОДНЫХ ПАРАМЕТРОВ" << endl;
	cout << "----------------------------------" << endl;
	cout << "Сила F " << F << " [кГ] " << endl;
	cout << endl << "Характеристики материала: " << endl;
	cout << "Модуль упругости материала E " << E << " [Па] " << endl;
	cout << "Предел прочности sigma " << sigma << " [Па] " << endl;
	cout << "Коэфф. запаса просночти n " << n << endl;
	cout << endl << "Геометрические размеры: " << endl;
	cout << "Нач. угол внутр. поверхности alfa_v " << alfa_vi << " [Гр.] " << endl;
	cout << "Запрещенный угол alfa_zapr " << alfa_zapr << " [Гр.] " << endl;
	cout << "Внутренний радиус Rv " << Rv << " [м] " << endl;
	cout << "Наруженый радиус Rn " << Rn << " [м] " << endl;
	cout << "Толщина кольца b " << b << " [м] " << endl;
	cout << "----------------------------------" << endl;
	
	F /= 10;				    //Сила в [Н]
	
	h = Rn - Rv;			    // Толщина кольца h
					
	m_rad = (Rn + Rv) / 2;		//Средний радиус кольца m_rad	
	
	J = b * pow(h, 3) / 12;		// Осевой момент инерции поперечного сечения J
	
	delta_x = 0.137 / 2 * (F * pow(m_rad, 3)) / (E * J); // Изменение радиуса кольца в направлении x
	
	delta_y = -0.149 / 2 * (F * pow(m_rad, 3)) / (E * J);

	// Средние полуоси эллипса
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

	// Площадь поперечного сечения кольца
	S_pop = b * h;

	for (int i = 0; i < count; i++)
	{
		x[i] = a1 * cos(alfa[i]); // Координаты по X текущих точек эллипса
		y[i] = b1 * sin(alfa[i]); // Координаты по Y текущих точек эллипса

		z[i] = sin(gamma[i]);
		u[i] = cos(gamma[i]);

		// Средний радиус кривизны эллипса 
		ro[i] = pow(a1, 2) * pow(b1, 2) * pow((pow(x[i], 2) / pow(a1, 4)) + (pow(y[i], 2) / pow(b1, 4)), 1.5);

		
		

		// Текущие точки наружной поверхности эллипса
		Rn_ellipse[i] = ro[i] + h / 2;
	    // Текущие точки внутренней поверхности эллипса
		Rv_ellipse[i] = ro[i] - h / 2;
		// Радиус кривизны нейтрального слоя
		rad[i] = h / log(Rn_ellipse[i] / Rv_ellipse[i]);

		

		// Крутящий момент отн. текущей точки на ср. линии кольца
		M_krut[i] = F * m_rad * (0.3183 - z[i] / 2);
		// Растягивающая сила отн. текущей точки на ср. линии кольца
		Nm[i] = -1 / 2 * F * z[i];
		// Сила сдвига в текущей точке
		Qm[i] = -1 / 2 * F * u[i];

		// Напряжение, создаваемое моментом M в точках поперечного сечения кольца
		// на внутр. поверхности:
		sigma_v[i] = (M_krut[i] * (Rv_ellipse[i] - rad[i])) / (S_pop * Rv_ellipse[i] * (ro[i] - rad[i]));
		// на наружной поверхности:
		sigma_n[i] = (M_krut[i] * (Rn_ellipse[i] - rad[i])) / (S_pop * Rn_ellipse[i] * (ro[i] - rad[i]));

		// Должно выполняться условие sigma[i] < sigma, иначе break

		// Напряжение, создаваемое силой растяжения во всех точках поперечного сечения кольца
		sigma_norm[i] = -(1/2) * (F * z[i] / S_pop);

		// Суммарное напряжение на внутренних поверхностях поперечного сечения кольца
		sigma_sum_v[i] = sigma_v[i] + sigma_norm[i];
		// Суммарное напряжение на наружных поверхностях поперечного сечения кольца
		sigma_sum_n[i] = sigma_n[i] + sigma_norm[i];

		if ((sigma_sum_v[i] >= n * sigma) || (sigma_sum_n[i] >= n * sigma))
		{
			cout << "Превышен предел деформации. Ошибка корректности исх. данных! ";
			continue;
		}

		//Деформация внутренней поверхности кольца
		eps_v[i] = sigma_sum_v [i] / E;
		//Деформация наружней поверхности кольца
		eps_n[i] = sigma_sum_n[i] / E;

	//  cout << "(" << alfa[i] * (180 / PI) <<";"<< eps_n[i] <<") " << endl;
	//	cout << "(" << alfa[i] * (180 / PI) << ";" << eps_v[i] << ") " << endl;

		eps_max_n = eps_n[0];
		eps_min_v = eps_v[count - 1];

		max_M_krut = M_krut[0];
		min_M_krut = M_krut[count - 1];
		//!!!!!!! Здесь нужно начертить графики eps_n(alfa[i]) и eps_v(alfa[i]) !!!!!
		
		
	}

	alfa_vi = alfa_vi * (PI / 180); // Внутренний угол [рад.]
	
	alfa_v = (PI / 2) - alfa_vi;				// Угол участка внутренней поверхности для тензоэлемента [рад.]
	alfa_n = (alfa_v * Rv) / Rn;				// Угол участка внешней поверхности для тензоэлемента [рад.]
	alfa_ni = (PI / 2 - alfa_n) * (180 / PI);	// Наружный угол, при котором Ln = Lv [гр.]
	Lvo = Rv * alfa_v;							// Длина на внутренней поверхности [м]
	Lno = Rn * alfa_n;							// Длина на наружной поверхности [м]

	/* !!!!!! Добавить к графикам eps_n(alfa[i]) и eps_v(alfa[i]) вертикальные линии (внутр., внеш., и запр. углы) 
	for (int i = -round(F); i < round(eps_n[count-1]*1000000); i++)
	{
		alfa_vi * (180 / PI) от i
		alfa_ni от i
		90 - alfa_zapr от i
	}

	*/
	//--------
	// Найдем изменение внутр. и наруж. поверхностей кольца интегралом по методу Симпсона.
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

	// Изменение длины наружней поверхности кольца. 
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

	// Изменение длины наружней поверхности кольца. 
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

	// Изменение длины внутренней поверхности кольца. Интеграл по методу Симпсона. 

	double Integral_v =  2 * (alfa[count - 1] - alfa[1]) * (s1v + 2 * s2v + 4 * s3v) / (6 * 79); 
	
	double LvF = 2 * (Lvo + Integral_v);		// Длина внутренней повехрности кольца при действии силы 
	double LnF = 2 * (Lno + Integral_n);		// Длина наружней повехрности кольца при действии силы 

	alfa_vi *= 180 / PI;
	

	cout << endl;
	cout << "----------------------------------" << endl;
	cout << "РАСПЕЧАТКА ПРОМЕЖУТОЧНЫХ ДАННЫХ" << endl;
	cout << "----------------------------------" << endl;
	cout << " Ширина кольца " << h << " [м] " << endl;
	cout << " Средний радиус кольца " << m_rad << " [м] " << endl;
	cout << " Осевой момент инерции " << J << " [м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Изменение раадиуса кольца в направлении x " << delta_x << " [м] " << endl;
	cout << " Изменение радиуса кольца в направлении y " << delta_y << " [м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Средняя полуось эллипса a1 " << a1 << " [м] " << endl;
	cout << " Средняя полуось эллипса b1 " << b1 << " [м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Площадь поперечного сечения " << S_pop << " [м^2] " << endl;
	cout << " Максимальный крутящий момент " << max_M_krut << " [Н*м] " << endl;
	cout << " Минимальный крутящий момент " << min_M_krut <<" [Н*м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Угол для тензоэлемента на наруж. поверхности " << alfa_n  <<" [рад] " << endl;
	cout << " Угол для тензоэлемента на внутр. поверхности " << alfa_v  <<" [рад] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Длина участка на наруж. поверхности " << Lno << " [м] " << endl;
	cout << " Длина участка на внутр. поверхности " << Lvo << " [м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Изменение длины на наруж. поверхности " << Integral_n << " [м] " << endl;
	cout << " Изменение длины на внутр. поверхности " << Integral_v << " [м] " << endl;
	cout << "----------------------------------" << endl;
	cout << " Длина наружней поверхности кольца при действии силы " << LnF << " [м] " << endl;
	cout << " Длина внутренней поверхности кольца при действии силы " <<LvF << " [м] " << endl;
	cout << "----------------------------------" << endl;

	cout << "----------------------------------" << endl;
	cout << " РАСПЕЧАТКА СТАТ. ХАРАКТЕРИСТИКИ" << endl;
	cout << "----------------------------------" << endl;
	cout << "     Внутренняя поверхность " << endl;
	cout << "----------------------------------" << endl;
	cout << "  Угол [гр.]\t\tДеформация 10^5 " << endl;

	

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
	cout << "       Наружняя поверхность " << endl;
	cout << "----------------------------------" << endl;
	cout << "  Угол [гр.]\t\tДеформация 10^5 " << endl;

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
	int i; // счётчик

	double Integral;
	double a = 0.0, b = alfa_n; // задаём отрезок интегрирования = 0,436775
	double h = 0.01;// задаём шаг интегрирования
	double n; // задаём число разбиений n

	n = (b - a) / h;

	// вычисляем интеграл по формуле Симпсона
	
	Integral = h * (f(a) + f(b)) / 6.0;

	for (i = 1; i <= n; i++)

		Integral = Integral + 4.0 / 6.0 * h * f(a + h * (i - 0.5));

	for (i = 1; i <= n - 1; i++)
		Integral = Integral + 2.0 / 6.0 * h * f(a + h * i);
	*/

