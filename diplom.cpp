#include "stdlib.h"
#include <stdio.h>
#include <time.h> 
#include <chrono>
#include <iostream>
#include <fstream>
#include <string.h>
#include "math.h"
#include <boost/algorithm/string/replace.hpp>
//автор Бабкин Михаил 2020г
using namespace std;

#include "consts.impl"
#include "otherFuncs.impl"

/**-------------------------------------------------------------------=========================================------
//строка для запуска программы через консоль linux: DIPLOM
//g++ -std=c++14 -fopenmp diplom.cpp && ./a.out
--------------------------------------------------------------------=========================================-----**/
double array_A[Nr] = {0};

void processing(int NOW, int MAX){
	if (NOW <= 0)
		cout << "00%";

	int proc = 0;
	proc = 100 * NOW / MAX;
	char b = 8;
	cout << b << b << b;
	if (proc < 10)
		cout << "0" << proc << "%";
	else
		cout << proc << "%";
	if (proc == 100)
		cout << " Complite " << endl;
	fflush(stdout);
}

double decubeD(int t){
	return  pow((-0.0002*t*t + 0.0592*t + 3.8005),1.0/3.0);	
}

bool file_input(double **f, int param, string file_name){	
	ofstream file_stream(file_name);	
	if(file_stream.is_open()){
		file_stream << "r\t";
		for(int t=0; t < Nt; t++){ //формируем шапку
			file_stream << "time" << t << "\t";
		}
		file_stream << endl;		
		for(int r = 0; r< Nr; r++){
			for(int t=-1; t < Nt; t++){	
				if(t<0)	{ file_stream << 2*(delta_r*r+r_min)*pow(10,6)*r0 << "\t";}
				else 	{ file_stream << f[t][r] << "\t";}
			}
			file_stream << endl;
		}
		file_stream.close();
		return true;
	}
	file_stream.close();
	return false;
}

double FF (int i){
	return delta_r*i+r_min;
}

int analog_f(double X){
	int start_point = 0;
	
	if(X> FF(Nr/2)){ start_point = (X <(FF(Nr/2+Nr/4)))?(-1+Nr/2):(Nr/2+Nr/4 -1);
	}else{ start_point = (X <(FF(Nr/4)))?(0):(Nr/4 -1);}
	
	for(int i = start_point;i<Nr;i++)
		if(X<=FF(i))return i;
	
	return Nr;

}

long double V(int r) {
	return (4.0 / 3.0)*3.141596*pow(2*(delta_r*r+r_min), 3.0); //по радиусу безразмерный
}

double f_start(double x){
	
	x = (delta_r*x+r_min)*2.0*r0*pow(10.0, 6.0);
	
	double mass_corrector = 0.2; //снизить для снижения экстремума по ординате, не влияет на средний размер при введении своих значений стоит отрегулировать им массу, так что бы она не превышала 1
	double sigma = 10;
	double math_oj=averStartSize;
	return mass_corrector*exp((-pow(x-math_oj,2))/(2.0*sigma*sigma))/(sigma*pow(2*3.141598, 0.5));
	
	/*  //старое распределение в виде горба, крайне не достоверное((
		x = (delta_r*x+r_min+0.01);
		const int N_p = 4;
		double p[N_p] = { 2.72737678945911e+22,	-4.53079771200336e+18,	146886712364303,	-36751521.1442473};
		double f = 0;
		for(int i = 0; i < N_p;i++)
			f += pow(x*r0,i) * p[N_p - i - 1];  /*в микронах для радиуса 
		
		f /= f0; //обезразмериваем
		
		if (f > f_max) f_max = f;
		if (f < 0 or x > 50) return 0;
		return f;
	*/
	
}

double get_Gamma(){
	//type Al2O3drinding -1 , Dry grinding - 0, SAS - 1, Isopropyl - 2 , Ethanol - 3, Al2O3 - 4
	double Z1 = gammaBabk;
	switch (typePAV) {
		case -1:
			Z1 = gammaBabk; //model by Babkin
			break;
		case 0: 			//models by Ivannikov
			Z1 = surfaceEnergy;
			break;
		case 1:
			//c /= 0.001;
			Z1 = surfaceEnergy + p3*c+p33*c*c/1000;
			break;
		case 2:
			Z1 = surfaceEnergy + p1;
			break;
		case 3:
			Z1 = surfaceEnergy + p2;
			break;
		case 4:
			Z1 = gammaBabk;
			break;
		default:
			break;
	}
	return Z1;
}

double getEpsIv(){
    return b0 + b1*x1 + b11*x1*x1 + b2*x2 + b22*x2*x2;
}

double getEpsBa(){
    return a0 + a1 * z1 + a12 * z1 * z2 + a2 * z2;
}

void init_A (){
	
	double eps = isBabkin?getEpsBa():getEpsIv();
	
	double L0 = 1 /t0; 	
	for(int r =0;r<Nr;r++){
		double We = (ps) * pow(2*(delta_r*r+r_min)*r0, 5.0 / 3.0) * pow(eps, 2.0 / 5.0) / get_Gamma(); //размерный [1]
		array_A[r] = (r < Barrier_A)?(0):(We * L/L0); //если обнулять дробление не нужно для данного размера частиц обезразмериваем А
	}
}

std::string printRadius(double **f, int t){
	stringstream strStream;
	
	for(int r =0; r < Nr; r++){
		strStream << r0*(r*delta_r + r_min) << "," ;
	}
	strStream << "end";
	string out(strStream.str());
	boost::replace_all(out, ",end", "");
	return out;
}

std::string printSizes(double **f, int t){
	stringstream strStream;
	
	for(int r =0; r < Nr; r++){
		strStream << f[t][r] << "," ;
	}
	strStream << "end";
	string out(strStream.str());
	boost::replace_all(out, ",end", "");
	return out;
}

double A(int r) {
	return  array_A[r];
}

double Betta(double x, double y){
	
	double summ = 0.0;
	double h = 0.000001;
	double t_start = 0.00001;
	for(int t = 0; t < int(1.0/h); t++){		
		summ += (pow(t_start+double(t)*h,x-1.0) * pow(1.0-t_start+double(t)*h,y-1.0)) + (pow(double(t+1.0)*h,x-1.0) * pow(1.0-double(t+1.0)*h,y-1.0));
	}
	cout << summ * h/2.0 << endl;
	return summ * h/2.0;
	
}

double B(double i, double k) {
	long double Vk = V(k);		//родительская часица
	long double Vi = V(i);		//дочерняя частица
	double _B;
	double C = 1.2;
	//_B = (Vk>Vi)?(0.5*( C/(delta_r*k+r_min) + ((1.0-C/2.0)/(delta_r*k+r_min))*(24.0*(pow(((delta_r*i+r_min)/(delta_r*k+r_min)),2)) - 24.0*(delta_r*i+r_min)/(delta_r*k+r_min) + 6.0))):(0);return P*_B; 		//параболическое ядро
	//_B = (Vk>Vi)?(0.5*( C/(Vk) + ((1.0-C/2.0)/(Vk))*(24.0*(pow(((Vi)/(Vk)),2)) - 24.0*(Vi)/(Vk) + 6.0))):(0); return P*_B;		//параболическое ядро
	
	
	_B = (Vk>Vi)?((30.0 / (Vk)) * (Vi / Vk) * (Vi / Vk) * (1.0 - Vi / Vk)* (1.0 - Vi / Vk)):(0); return P*_B;
	/*_B = (Vk>Vi)?((30.0 / (delta_r*k+r_min)) * ((delta_r*i+r_min) / (delta_r*k+r_min)) * 
		((delta_r*i+r_min)  / (delta_r*k+r_min)) * 
		(1.0 - (delta_r*i+r_min)  / (delta_r*k+r_min))* (1.0 - (delta_r*i+r_min)  / (delta_r*k+r_min))):(0);
		return P*_B;*/
	
	
	//уравнения Хилла-Нгу
	double z = Vi/Vk;
	double 	w1 = 1.0, //режим "истирание"
			w2 = 0.5,
			p1 = 2.0,
			p2 = 2.0,
			q1 = 0.0001,
			q2 = 1,
			k1 = 1,
			k2 = 0.0001;	
	double Q =  (Vk>Vi)?(w1*p1* (pow(z,q1-1.0)*pow(1-z, k1-1.0))/Betta(q1, k1)  + 
				w2*p2* (pow(z,q2-1.0)*pow(1-z, k2-1.0))/Betta(q2, k2)):(0); return  Q/Vk;
	
	

}

double SumInegral(double **f, int t, int r) {
	double SumI = 0.0;

	#pragma omp parallel for reduction(+:SumI) num_threads(Num_T)
	for (int i = r; i < Nr - 1; i++){
		SumI += f[t][i] * A(i) * B(r, i) + f[t][i + 1] * A(i + 1) * B(r, i + 1);
	}

	SumI *= delta_r / 2.0;
	return  SumI;

}

double summ_f(double **f, int t){
	double result_sum =0;
	#pragma omp parallel for reduction(+:result_sum) num_threads(Num_T)
	for(int i = 0; i < Nr; i++)
		result_sum += f[t][i];
	return result_sum;
}

double RadiusMean(double **f, int t) {
	
	double a=0;
	double b=0;
	
	for (int i = 0; i < Nr - 1; i++){
		a += f[t][i] * (delta_r*double(i) + r_min)+ f[t][i+1] * (delta_r*double(i+1) + r_min);
		b += f[t][i] + f[t][i + 1];
	}
	a *= delta_r / 2.0;
	b *= delta_r / 2.0;
	return a / (b);
}

double Deviance(double x, double **f, int t){
	double sm = 0;
	double sf = 0;
	for(int i = 0; i < Nr-1; i++){
		sm += f[t][i]*f0 * (pow((delta_r*double(i)*r0 + r_min*r0 - x*r0),2.0));
		sf = sf + f[t][i]*f0;
	}
	double _Dev = sm/sf;
	return pow(_Dev, 0.5);

}

double Mass_Of_All_Particle(double **f, int t) {
	double Volume = 0;
	for (int r = 0; r < Nr-1; r++) {
		Volume += (V0*V(r)*f[t][r]*f0 + V0*V(r+1)*f[t][r+1]*f0) * delta_r/2.0;
	}
	return ps*Volume;
}

void magic_func(double ** f, int t,  double target_mass){
	double bad_mass = Mass_Of_All_Particle(f, t);
	
	for(int i = 0; i< Nr; i++)
		f[t][i] *= (target_mass/bad_mass);
	
}

double SumInegral_with_V(double **f, int t, int r) {
	double SumI = 0.0;

	#pragma omp parallel for reduction(+:SumI) num_threads(Num_T)
	
	for (int i = r; i < Nr - 1; i++){
		SumI += V(i)*f[t][i] * A(i) * B(r, i) + V(i+1)*f[t][i + 1] * A(i + 1) * B(r, i + 1);
	}

	SumI *= delta_r / 2.0;
	return  SumI;

}

double verification_balance(double **f, int t){
	double left_summ = 0;
	double right_summ = 0;
	#pragma omp parallel for reduction(+:left_summ,right_summ) shared(f,t) num_threads(Num_T)
	for(int i = 0; i<Nr-1; i++){
		left_summ += (V(i)*f[t][i]*A(i) + V(i+1)*f[t][i+1]*A(i+1));
		right_summ += (SumInegral_with_V(f,t, i) + SumInegral_with_V(f,t,i+1));
	}
	double balance = (right_summ * delta_r) / 2.0  - (left_summ * delta_r) / 2.0;	
	
	return balance;
}

void log_out(double dM, double r_0, double r_f, double time_proc, double Deviance){
	ofstream file0("log.txt", ios_base::app);
	time_t t;
    time(&t);
    tm* local = localtime(&t);
    char Hours[10] = {0};
    char Minutes[10] = {0};
    char Date [20] = {0};
    strftime(Hours, 	sizeof(Hours)/sizeof(Hours[0]), 	"%H", local);
    strftime(Minutes, 	sizeof(Minutes)/sizeof(Minutes[0]), "%M", local);
    strftime(Date, 		sizeof(Date)/sizeof(Date[0]), 	"%d %B %Y", local);	
	if (file0.is_open()) {
		file0<< " \nВывод_параметров_во_время_вычисления_на " << Hours << ":" << Minutes << " " << Date << endl;
		file0<< "f0 " 		<< f0 	<< endl ;
		file0<< "r0 " 		<< r0 	<< endl;
		file0<< "V0 " 		<< V0 	<< endl;
		file0<< "A0 " 		<< A0 	<< endl ;
		file0<< "B0" 		<< B0 	<< endl ;
		file0<< "m0 " 		<< m0 	<< endl ;
		file0<< "Nr "		<< Nr	<< endl ;
		file0<< "Nt " 		<< Nt 	<< endl ;
		file0<< "T_max " 	<< T_max*t0/60.0<< endl ;
		file0<< "ps " 		<< ps 	<< endl ;	
		file0<< "m_sphere " << m_sphere << endl ;
		file0<< "m_poroh "  << m_poroh  << endl ;
		file0<< "d_sphere " << d_sphere << " mm" << endl ;	
		file0<< "delta_r "  << delta_r 	<< endl ; // шаг по радиусу обезразмеренный
		file0<< "delta_t " 	<< delta_t	<< endl;  //обезрамеренная
		file0<< "L " 		<< L 	<< endl ;	
		file0<< "P " 		<< P	 << endl ;	
		file0<< "Barrier A (micro) " << 2*pow(10,6)*(Barrier_A*delta_r+r_min)*r0	 << endl ;	
		file0<< "z1 " 		<< (m_sphere / m_poroh) << endl;
		file0<< "z2 " 		<< 1.0 / d_sphere << endl;
		file0<< "f_max " 	<< f_max << endl;
		file0<< "T_max(мин) " 	<< T_max*t0<< endl ;
		file0 << "M0/M_fin = " 	<< dM << endl;
		file0 << "d_исходная micro = " << r_0*2*pow(10, 6) << endl;
		file0 << "d_конечная micro = " << r_f*2*pow(10, 6) << endl;
		file0 << "SquareMeanDeviance = " << Deviance*2 << endl;
		file0 << "время_расчёта_в_минутах = " << time_proc << endl;
	}else{
		cout<< "Ошибка вывода файла логов"<<endl;
	}
	file0.close();
}

int main(int argc, char* argv[]) {
	setlocale(LC_ALL, "rus"); // корректное отображение Кириллицы
	
	if(argc > 1){
		for(int i=0;i < argc; i++ ){
			//cout << argv[i] << endl;
			if(std::string(argv[i]).compare("-h") == 0 || std::string(argv[i]).compare("--help") == 0){
				
				cout  	<< "-L [dbl] - феноменологический коэффициент" << endStr
						<< "-typePAV [int] - тип пов. активного вещ-ва  Al2O3drinding: -1 , Dry grinding: 0, SAS: 1, Isopropyl: 2 , Ethanol: 3 " << endStr
						<< "-Con [dbl] масса SiC г." << endStr
						<< "-oborot [dbl] 250 или 400 оборотов в минуту" << endStr
						<< "-densBalls [dbl] плотность шаров кг/м^3" << endStr
						<< "-densParticle [dbl] плотность частиц кг/м^3" << endStr
						<< "-typeMill [int] тип мельницы, 1 если модель Иванникова, 0 - модель Бабкина" << endStr
						<< "-P [dbl] прединтегральный коэффициент" << endStr
						<< "-massBallsAndParticles [dbl] [dbl] масса шаров и масса частиц" << endStr
						<< "-massRatio [dbl] отношение массы шаров и массы частиц" << endStr
						<< "-sizeBall [dbl] размер мелющих шаров, мм" << endStr
						<< "-searchSize [dbl] размер который необходимо получить (в микронах)" << endStr
						<< "-avStart [dbl] average particle size (23,7)micron " << endStr
						<< "-Tmax [dbl] время дробления в секундах" << endl;
				return 0;				 
			}
			if(std::string(argv[i]).compare("-typeMill") == 0){
				typeMill = std::atoi(argv[i+1]); //cout << "Type of mill init (" << ((typeMill==0)?("PlanetaryMill"):("CentralForceMill")) << ")" << endl;
				isBabkin = (typeMill==0);
				/*if(isBabkin){
					L = 7.72;
					P = 0.17;
				}else
				{ 
					L=19.6971;
					P=10.1343;
				}*/			 
			}
			if(std::string(argv[i]).compare("-searchSize") == 0){
				size_search = std::atof(argv[i+1])/std::atof(argv[i+2]); //cout << "Search size init ( " << size_search << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-massBallsAndParticles") == 0){
				z1 = std::atof(argv[i+1])/std::atof(argv[i+2]);// cout << "mass ratio init ( " << z1 << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-massRatio") == 0 || std::string(argv[i]).compare("-z1") == 0 ){
				z1 = std::atof(argv[i+1]); //cout << "mass ratio init ( " << z1 << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-z2") == 0 ){
				z2 = std::atof(argv[i+1]); //cout << "z2 init ( " << z2 << ")"<<endl;
				d_sphere = pow( std::atof(argv[i+1]),-1.0);
			}
			if(std::string(argv[i]).compare("-sizeBall") == 0 ){
				z2 = 1.0 / std::atof(argv[i+1]); //cout << "size balls init ( " << z2*pow(10.0, -3.0) << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-P") == 0){
				P = std::atof(argv[i+1]); //cout << "P init ( " << P << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-L") == 0){
				L = std::atof(argv[i+1]); //cout << "L init ( " << L << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-avStart") == 0){
				averStartSize = std::atof(argv[i+1]); 
			}
			if(std::string(argv[i]).compare("-typePAV") == 0){
				typePAV = std::atoi(argv[i+1]); //cout << "typePAV init (" << typePAV << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-con") == 0){
				c = std::atof(argv[i+1]); //cout << "Concentration (weight in terms of dry matter) init (" << c << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-oborot") == 0){
				x1 = std::atof(argv[i+1])/x10; //cout << "oborotov in sec init (" << x1*x10 << ")"<<endl;
			}
			if(std::string(argv[i]).compare("-densParticle") == 0){
				ps = std::atof(argv[i+1]); //cout << "density of particle init (" << ps << ")"<<endl;				
			}
			if(std::string(argv[i]).compare("-densBalls") == 0){
				x2 = std::atof(argv[i+1])/x20; //cout << "density of milling balls init (" << x2*x20 << ")"<<endl; //params.x2 = Double.parseDouble(args[i+1])/params.x20;
			}
			
			if(std::string(argv[i]).compare("-tMax") == 0){
				T_max = std::atof(argv[i+1])/t0; //cout << "Tmax init (" << T_max << ")"<<endl;
				delta_t = T_max / (double(Nt)); 
			}
		}
	}
	
	time_t start, end;
	r_max /= r0;
	r_min /= r0;
	
	init_A();
	if(d_sphere >= 5){
		P = 0.16;
		K_vol = pow(3.0, 1.0/3.0);
	}
	if(d_sphere <= 2){
		P = 0.18;
		K_vol = pow(4.0, 1.0/3.0);
	}

	/*
		cout << "Вывод параметров во время вычисления: " << endl ;	
		cout<< "r_max " << r_max << endl ;
		cout<< "r_min " << r_min << endl ;
		cout<< "t0 " 	<< t0 << endl ;
		cout<< "f0 " 	<< f0 << endl ;
		cout<< "g0 " 	<< g0 << endl ;
		cout<< "r0 " 	<< r0 << endl;
		cout<< "V0 " 	<< V0 << endl;
		cout<< "A0 " 	<< A0 << endl ;
		cout<< "B0" 	<< B0 << endl ;
		cout<< "m0 " 	<< m0 << endl ;
		cout<< "Nr "	<< Nr<< endl ;
		cout<< "Nt " 	<< Nt << endl ;
		cout<< "N_particle " 	<< N_particle << endl ;
		cout<< "T_max (мин)" << T_max*t0/60.0<< endl ;
		cout<< "L " 	<< L 	<< endl ;	
		cout<< "ps " 	<< ps 	<< endl ;	
		cout<< "m_sphere " 	<< m_sphere << endl ;
		cout<< "m_poroh "  	<< m_poroh  << endl ;
		cout<< "d_sphere " << d_sphere << " mm" << endl ;	
		cout<< "delta_r "  	<< delta_r << endl ; 	// шаг по радиусу обезразмеренный
		cout<< "delta_t " 	<< delta_t<< endl;  	//обезрамеренная
		cout<< "z1 " << z1 << endl;
		cout<< "z2 " << z2 << endl;
		cout<< "stable size " << 2*pow(10,6)*(Barrier_A*delta_r+r_min)*r0	 << endl ;	
		
	*/
	double **f = new double*[Nt];
	//cout << "delta_r " << delta_r << endl;
	
	for (int i = 0; i < Nt; i++) {//выделяем память под массив данных
		f[i] = new double[Nr];
	}
	
	double *f_half = new double[Nr];
	//cout << "Data array is reserved in memory" << endl;
	for (int t = 0; t < Nt; t++)
		for (int r = 0; r < Nr; r++) {
			f[t][r] = 0;
		}
		
	
	double u = 0, fff = 0;
	for (int i = 0; i < Nr; i++) {
		f_half[i] = 0;
		f[0][i] =  f_start(i);
	}
	//cout << "f_max " << f_max << endl;
	//cout << "Initialization complate!" << endl;
	//cout << "--koeff P-----  " << P << "   --------" << endl;
	double diametr_mean_0 = r0*RadiusMean(f, 0)*2;
	//cout << "средни ДИАМЕТР в микронах (t=0) =  \x1b[5;32m" << diametr_mean_0*pow(10, 6) << " \x1b[0m" << endl;
	//cout << "deviace (t=0)  = " << Deviance( RadiusMean(f, 0) ,f, 0) << endl;
	long double dt_integral = 0; 	//буфер для расчёта итеграла умноженного на дельта_тэ
	long double one_plus_dt_A = 0; 	// буфер для расчёт знаменателя разносностной схемы на шаге n+1
	//char * file_name = new char[100];
	//sprintf(file_name, "statestic for %.1f %.2f %.3f %.3f %.3f.txt", m_sphere/m_poroh, d_sphere, L, P, 2*pow(10,6)*(Barrier_A*delta_r+r_min)*r0);
	//cout << file_name << endl;
	//ofstream file_f_mean(file_name);
	//file_f_mean << "time\tdiametr\tmass\tparticle\tSquareMeanDeviance\n";
	
	double time_of_search_size = 0;
	double eps_r = 0.1*pow(10, -6)/r0;
	double diam_temp =0;
	double min_mass = 100;
	double mass_temp = 0;
	
	
	time(&start); //засекаем время расчёта
	//cout << "Прогресс расчёта:    "<<endl ;
	
	
	/*
	//------------------------------------------------------------------
	char * file_name_A = new char[255];
	sprintf(file_name_A, "A_distridution");
	ofstream file_A(file_name_A);
	file_A << "r\tA";
	for(int r = 0; r<Nr;r++){
		
		file_A << delta_r*r+r_min << "\t" << A(r) <<"\t" << B(r, 1500) << endl;
		
	}
	file_A.close();
	//------------------------------------------------------------------------
	return 0;
	*/
	double startSize;

	for (int t = 0; t < Nt - 1; t++) {
		if(t == 2) startSize = RadiusMean(f,t); //save mean first particle's raduis 
		
		if(d_sphere >= 5){
			//P = 0.087;
			if(t == 0) K_vol = pow(3.0, 1.0/3.0);
			if(t == 40) K_vol = pow(4.0, 1.0/3.0);
		}
		if(d_sphere <= 2){
			//P = 0.1;
			if(t == 0) K_vol = pow(4.0, 1.0/3.0);
			if(t == 27) K_vol = pow(5.0, 1.0/3.0);
			if(t == 45) K_vol = pow(6.0, 1.0/3.0);
		}
		
		for (int r = 0; r < Nr; r++) {
			dt_integral = 4*3.1415*(delta_r*r+r_min)*delta_t * SumInegral(f, t, analog_f((delta_r*r+r_min)*K_vol/*decubeD(t)*/));
			f_half[r] = (f[t][r]) + P*dt_integral; //считаем на половинном шаге
			one_plus_dt_A = 1.0 + delta_t*A(r)*A0;
			f[t + 1][r] = (f_half[r]) / (one_plus_dt_A); //считаем на шаге "n + 1"
		}
		//magic_func(f, t+1, Mass_Of_All_Particle(f, t));
		if((RadiusMean(f, t) < (size_search+eps_r)) && (RadiusMean(f, t) > (size_search-eps_r)))
			time_of_search_size = delta_t*t*t0/60;
		//processing(t+1, (Nt - 1));
		diam_temp = r0*RadiusMean(f, t)*2/pow(10,-6);
		mass_temp = Mass_Of_All_Particle(f, t);
		cout <<  "\n";
		//if(t > 15 && t%5 ==0){
			
		cout << "{" <<
				"\"L,P\":" <<"["<<L<<","<<P<<"]"<< "," <<
				"\"typeMill\":" << ((isBabkin)?("\"AluminiunOxide\""):("\"SiliconCarbide\"")) << "," <<
				"\"roShar\":" << ((isBabkin)?(5.68):(x2)) << "," <<
				"\"velocity\":" << ((isBabkin)?(2.5):(x1)) << ","<<
				"\"density\":" << ps << ","<<
				"\"avStSize\":"<< averStartSize<<","<<
				"\"SurfPowerShar\":" << ((isBabkin)?(7.480):((x2>6)?(15):(7.48))) << ","<<
				"\"SurfPowerMater\":" << ((isBabkin)?(gammaBabk):(get_Gamma())) << ","<<
				"\"massRatio\":" << ((isBabkin)?(z1):(4.0)) << ","<<
				"\"sizeShar\":" << ((isBabkin)?((1.0/z2)):(10.0)) << ","<<
				"\"Sizes\": [[" << (printRadius(f,t)) << "],[" << (printSizes(f,t)) << "]]," <<
				"\"t\":" << ((t+1)*delta_t *t0) << ","<<
				"\"averageSize\":" << diam_temp  << 
				"}\n";	
		//}

		if (mass_temp < min_mass) min_mass = mass_temp;
	}
	//cout << "Size_fin " << diam_temp << "\r\n";
	//file_f_mean.close();
	time(&end); //останавливаем счётчик времени
	/*double	mass_delta = Mass_Of_All_Particle(f, 0) / Mass_Of_All_Particle(f, Nt - 1);
	cout << "Расчёт окончен и занял времени (мин): " << difftime(end, start)/60<< endl;
	cout << "Deviance fin = " << Deviance(RadiusMean(f, Nt - 1),f, Nt - 1)/pow(10, -6) << endl;
	cout << "Size fin = " << "\x1b[32m"<<r0*RadiusMean(f, Nt - 1)*2 / pow(10, -6)<< " \x1b[0m" << endl;
	cout << "размер равный \x1b[33m"<< size_search*r0/pow(10, -6) << " \x1b[0m достигается в момент времени \x1b[5;32m"<< time_of_search_size << " \x1b[0m" << endl;
	cout << "M0/M_fin = " <<mass_delta<< endl;
	cout << "M0 - M_fin = " <<  Mass_Of_All_Particle(f, 0) - Mass_Of_All_Particle(f, Nt - 1)  << endl;
	cout << "M0/M_MIN = " << Mass_Of_All_Particle(f, 0) / min_mass << endl;
	
	//вывод в файл 
	//cout << "Начинаем вывод в файл Кол-во частиц в единице объема -> ";
	FILE *file_f_integr;
	//cout << "Начинаем вывод в файл интегрального распределения -> ";
	file_name = new char[100];
	
	sprintf(file_name, "results %.1f %.2f %f %f %f.txt", m_sphere/m_poroh, d_sphere, L, P, 2*pow(10,6)*(Barrier_A*delta_r+r_min)*r0);
	
	file_input(f, 1, file_name);
	
	log_out(mass_delta, diametr_mean_0, r0*RadiusMean(f, Nt - 1), difftime(end, start)/60, Deviance(RadiusMean(f, Nt - 1),f, Nt - 1));
	
	
	cout << "Закончили вывод в файлы" << endl;*/
	return 0;
}
