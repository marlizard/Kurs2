#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

typedef double real;

//Номер используемого теста
#define function 2

//Используется ли генератор сетки (0 - не используется, 1 - используется)
#define generator 1

class Method {
public:
	//Количество конечных элементов в сетке области
	int n;

	//Количество узлов в сетке по времени
	int time;

	//Тип используемых базисных функций
	int type;

	//Глобальная матрица
	vector<real> di;
	vector<real> au;
	vector<int> ia;

	//Матрица жесткости
	vector<real> gdi;
	vector<real> gau;

	//Матрица массы
	vector<real> mdi;
	vector<real> mau;

	//Массив координат узлов области
	vector<real> cord;

	//Массив координат узлов по времени
	vector<real> time_cord;

	//Глобальный вектор правой части
	vector<real> b;

	//Новый вектор весов базисных функций
	vector<real> new_q;

	//Старый вектор весов базисных функций
	vector<real> old_q;

	//Массив значений лямбда
	vector<real> lambda;

	//Массив значений гамма
	vector<real> gamma;

	//Массив значений сигма
	vector<real> sigma;

	//Флаг ошибки
	int flag;

	//Номера краевых условий
	int left, right;

	//Коэффициент для третьих краевых
	real beta;

	//Искомая функция
	real fun_u(real r, real t);

	//Производная искомой функции
	real derfun_u(real r, real t);

	//Функция вектора правой части
	real func(real r, int el, real t);

	//Производная функции вектора правой части
	real der_func(real r, int el, real t);

	//Расчет относительной погрешности
	real relative_error(real t);

	//Нахождение произведения матрицы массы на старый вектор весов
	vector<real> multply_additive_component();

	//Ввод данных
	void input();

	//Генератор сетки по области
	void area_grid_generator();

	//Генератор сетки по времени
	void time_grid_generator();

	//Выявление типа базисных функций
	void basis_type();

	//Создание портрета матрицы
	void matrix_portrait();

	//Инициализация пространства глобальной матрицы
	void matrix_clean();

	//Посчет локальной матрицы жесткости
	void local_G(int i);
	//Сборка матрицы жесткости
	void stiffness_matrix();

	//Подсчет локальной матрицы масс
	void local_M(int i);
	//Сборка матрицы массы
	void mass_matrix();

	//Подсчет локального вектора правой части
	void local_b(int i, real t);
	//Сборка вектора правой части
	void right_part_vector(int time_index);

	//Подсчет глобальной матрицы
	void global_matrix(int time_index);

	//Сохранение нового вектора весов как старого
	void save_as_previous();

	//Учет первого краевого условия
	void dirichlet_condition_left(real t);
	void dirichlet_condition_right(real t);

	//Учет второго краевого условия
	void neyman_condition_left(real t);
	void neyman_condition_right(real t);

	//Учет третьего краевого условия
	void newton_condition_left(real t);
	void newton_condition_right(real t);

	//Учет краевых условий
	void boundary_conditions(real t);

	//Разложение глобальной матрицы в LLt-формат
	void LLt();

	//Прямой и обратный ход для вычисления q
	void direct_reverse();

	//Основной алгоритм
	void iteration_process();

	//Вывод
	void output(real t, ofstream &out);
};


real Method::fun_u(real r, real t) {
	switch (function) {
	case 0: {
		return 5;
		break;
	}
	case 1: {
		return 5 * t;
		break;
	}
	case 2: {
		return 5 * t * t;
		break;
	}
	case 3: {
		return r * r;
		break;
	}
	case 4: {
		return r * r * t;
		break;
	}
	case 5: {
		return r * r * t * t;
		break;
	}
	case 6: {
		return t * t * cos(r);
		break;
	}
	}
}


real Method::derfun_u(real r, real t) {
	switch (function) {
	case 0: {
		return 0;
		break;
	}
	case 1: {
		return 0;
		break;
	}
	case 2: {
		return 0;
		break;
	}
	case 3: {
		return 2 * r;
		break;
	}
	case 4: {
		return 2 * r * t;
		break;
	}
	case 5: {
		return 2 * r * t * t;
		break;
	}
	case 6: {
		return - t * t * sin(r);
		break;
	}
	}
}


real Method::func(real r, int el, real t) {
	switch (function) {
	case 0: {
		return 0;
		break;
	}
	case 1: {
		return 5 * sigma[el];
		break;
	}
	case 2: {
		return 10 * t * sigma[el];
		break;
	}
	case 3: {
		return -6 * lambda[el];
		break;
	}
	case 4: {
		return -6 * lambda[el] * t + sigma[el] * r * r;
		break;
	}
	case 5: {
		return -6 * lambda[el] * t * t + 2 * sigma[el] * r * r * t;
		break;
	}
	case 6: {
		return 2 * lambda[el] * t * t * sin(r) / r + lambda[el] * t * t * cos(r) + 2 * t * sigma[el] * cos(r);
		break;
	}
	}
}


real Method::der_func(real r, int el, real t) {
	switch (function) {
	case 0: {
		return 0;
		break;
	}
	case 1: {
		return 0;
		break;
	}
	case 2: {
		return 0;
		break;
	}
	case 3: {
		return 0;
		break;
	}
	case 4: {
		return 2 * sigma[el] * r;
		break;
	}
	case 5: {
		return 4 * sigma[el] * r * t;
		break;
	}
	case 6: {
		return - 2 * t * sigma[el] * sin(r) - lambda[el] * t * t * sin(r) + 2 * lambda[el] * t * t * r * cos(r) / r / r - 2 * lambda[el] * t * t * sin(r) / r / r;
		break;
	}
	}
}


real Method::relative_error(real t) {
	real error_sum = 0, real_sum = 0;

	int b = 0;

	for (int i = 0; i <= n; i++) {
		error_sum += (new_q[b] - fun_u(cord[i], t)) * (new_q[b] - fun_u(cord[i], t));
		real_sum += fun_u(cord[i], t) * fun_u(cord[i], t);
		if (type == 1) b++;
		else b += 2;
	}
	return sqrt(error_sum / real_sum);
}


vector<real> Method::multply_additive_component() {
	vector<real> m(old_q.size(), 0);
	if (type == 1) {
		m[0] = mdi[0] * old_q[0] + mau[0] * old_q[1];

		for (int i = 1; i < n; i++)
			m[i] = mau[i - 1] * old_q[i - 1] + mdi[i] * old_q[i] + mau[i] * old_q[i + 1];

		m[n] = mau[n - 1] * old_q[n - 1] + mdi[n] * old_q[n];
		//для линейного базиса
	}
	else {
		int block = b.size() / 2 - 1;

		m[0] = mdi[0] * old_q[0] +
			   mau[0] * old_q[1] +
			   mau[1] * old_q[2] + 
			   mau[3] * old_q[3];

		m[1] = mau[0] * old_q[0] +
			   mdi[1] * old_q[1] +
			   mau[2] * old_q[2] +
			   mau[4] * old_q[3];

		for (int i = 1; i < block; i++){
			m[2 * i] = mau[5 * (i - 1) + 1] * old_q[2 * i - 2] +
					   mau[5 * (i - 1) + 2] * old_q[2 * i - 1] +
					   mdi[2 * i] * old_q[2 * i] +
					   mau[5 * (i - 1) + 5] * old_q[2 * i + 1] + 
					   mau[5 * (i - 1) + 6] * old_q[2 * i + 2] + 
					   mau[5 * (i - 1) + 8] * old_q[2 * i + 3];

			m[2 * i + 1] = mau[5 * (i - 1) + 3] * old_q[2 * i - 2] +
						   mau[5 * (i - 1) + 4] * old_q[2 * i - 1] +
						   mau[5 * (i - 1) + 5] * old_q[2 * i] +
						   mdi[2 * i + 1] * old_q[2 * i + 1] +
						   mau[5 * (i - 1) + 7] * old_q[2 * i + 2] +
						   mau[5 * (i - 1) + 9] * old_q[2 * i + 3];
		}

		m[2 * block] = mau[5 * (block - 1) + 1] * old_q[2 * block - 2] +
					   mau[5 * (block - 1) + 2] * old_q[2 * block - 1] +
					   mdi[2 * block] * old_q[2 * block] +
					   mau[5 * (block - 1) + 5] * old_q[2 * block + 1];

		m[2 * block + 1] = mau[5 * (block - 1) + 3] * old_q[2 * block - 2] +
						   mau[5 * (block - 1) + 4] * old_q[2 * block - 1] +
						   mau[5 * (block - 1) + 5] * old_q[2 * block] +
						   mdi[2 * block + 1] * old_q[2 * block + 1];
	}
	return m;
}


void Method::input() {
	beta = 1;
	flag = 0;

	if (generator == 0) {
		ifstream elem("elem.txt");
		elem >> n;

		lambda.resize(n, 0);
		gamma.resize(n, 0);
		sigma.resize(n, 0);

		int start, finish;

		for (int i = 0; i < n; i++) {
			elem >> start;
			elem >> finish;

			if (i != start || i != finish - 1) {
				cout << "Lattice creation error";
				flag = 1;
				return;
			}

			elem >> lambda[i];
			elem >> gamma[i];
			elem >> sigma[i];
		}
		elem.close();

		int s;
		ifstream nodes("nodes.txt");
		nodes >> s;
		if (s != n + 1) {
			cout << "Not correct number of nodes";
			flag = 1;
			return;
		}

		cord.resize(n + 1, 0);
		for (int i = 0; i <= n; i++) nodes >> cord[i];
		nodes.close();

		ifstream time_file("time_cord.txt");
		time_file >> time;
		time_cord.resize(time, 0);

		for (int i = 0; i < time; i++)
			time_file >> time_cord[i];

		time_file.close();
	}
	else {
		area_grid_generator();
		time_grid_generator();
	}

	ifstream bound_cond("bound_cond.txt");
	bound_cond >> left;
	bound_cond >> right;
	bound_cond.close();
	if (left < 1 || left > 3 || right < 1 || right > 3) {
		cout << "Wrong type of boundary conditions";
		flag = 1;
		return;
	}
}


void Method::area_grid_generator() {
	real start, finish, k, h, r;
	ifstream area("area_grid_gen.txt");
	area >> n >> start >> finish >> k;
	area.close();
	lambda.resize(n, 0);
	gamma.resize(n, 0);
	sigma.resize(n, 0);
	cord.resize(n + 1, 0);
	for (int i = 0; i < lambda.size(); i++) {
		lambda[i] = 1;
		gamma[i] = 1;
		sigma[i] = 1;
	}
	if (k == 1)
		h = (finish - start) / n;
	else
		h = (finish - start) * (k - 1) / (pow(k, n) - 1);
	r = start;
	for (int i = 0; i < cord.size(); i++) {
		cord[i] = r;
		r += pow(k, i) * h;
	}
}


void Method::time_grid_generator() {
	real start, finish, k, h, t;
	ifstream time_file("time_grid_gen.txt");
	time_file >> time >> start >> finish >> k;
	time_file.close();
	time_cord.resize(time, 0);
	if (k == 1)
		h = (finish - start) / (time - 1);
	else
		h = (finish - start) * (k - 1) / (pow(k, (time - 1)) - 1);
	t = start;
	for (int i = 0; i < time_cord.size(); i++) {
		time_cord[i] = t;
		t += pow(k, i) * h;
	}
}


void Method::basis_type() {
	if (flag) return;
	cout << "Choose the type of basis functions:" << endl;
	cout << "1. Linear functions" << endl;
	cout << "2. Cubic hermitian functions" << endl;
	cin >> type;
	if (type != 1 && type != 2) {
		cout << "Wrong type of basis functions";
		flag = 1;
	}
}


void Method::matrix_portrait() {
	basis_type();

	if (flag) return;

	if (type == 1) {
		ia.resize(n + 2, 0);
		ia[0] = 0;
		for (int i = 0; i < n; i++) ia[i + 1] = i;
		ia[n + 1] = n;

		di.resize(n + 1, 0);
		au.resize(n, 0);

		gdi.resize(n + 1, 0);
		gau.resize(n, 0);

		mdi.resize(n + 1, 0);
		mau.resize(n, 0);

		b.resize(n + 1, 0);
		new_q.resize(n + 1, 0);
		old_q.resize(n + 1, 0);
	}
	else {
		ia.resize(2 * n + 3, 0);
		ia[0] = 0;
		ia[1] = 0;
		int index = 1;
		for (int i = 0; i < n; i++) {
			ia[2 * i + 2] = index;
			index += 2;
			ia[2 * i + 3] = index;
			index += 3;
		}
		ia[2 * n + 2] = index;

		di.resize(2 * n + 2, 0);
		au.resize(index, 0);

		gdi.resize(2 * n + 2, 0);
		gau.resize(index, 0);

		mdi.resize(2 * n + 2, 0);
		mau.resize(index, 0);

		b.resize(2 * n + 2, 0);
		new_q.resize(2 * n + 2, 0);
		old_q.resize(2 * n + 2, 0); 
	}
}


void Method::matrix_clean() {
	for (int i = 0; i < di.size(); i++) di[i] = 0;
	for (int i = 0; i < au.size(); i++) au[i] = 0;
	for (int i = 0; i < b.size(); i++) b[i] = 0;
	for (int i = 0; i < new_q.size(); i++) new_q[i] = 0;
}


void Method::local_G(int i) {
	real h = cord[i + 1] - cord[i];
	if (type == 1) {
		real c = lambda[i] * (cord[i] * cord[i] + cord[i] * h + h * h / 3) / h;
		gdi[i] += c;
		gdi[i + 1] += c;
		gau[i] -= c;
	}
	else {
		//Диагональные элементы
		gdi[2 * i] += lambda[i] * (6 * cord[i] * cord[i] / 5 + 6 * h * cord[i] / 5 + 12 * h * h / 35) / h;
		gdi[2 * i + 1] += lambda[i] * h * (2 * cord[i] * cord[i] / 15 + h * cord[i] / 15 + 2 * h * h / 105);
		gdi[2 * i + 2] += lambda[i] * (6 * cord[i] * cord[i] / 5 + 6 * h * cord[i] / 5 + 12 * h * h / 35) / h;
		gdi[2 * i + 3] += lambda[i] * h * (2 * cord[i] * cord[i] / 15 + h * cord[i] / 5 + 3 * h * h / 35);

		//Элементы нижнего треугольника
		gau[5 * i] += lambda[i] * (cord[i] * cord[i] / 10 + h * cord[i] / 5 + h * h / 14);
		gau[5 * i + 1] -= lambda[i] * (6 * cord[i] * cord[i] / 5 + 6 * h * cord[i] / 5 + 12 * h * h / 35) / h;
		gau[5 * i + 2] -= lambda[i] * (cord[i] * cord[i] / 10 + h * cord[i] / 5 + h * h / 14);
		gau[5 * i + 3] += lambda[i] * (cord[i] * cord[i] / 10 - h * h / 35);
		gau[5 * i + 4] -= lambda[i] * h * (cord[i] * cord[i] / 30 + h * cord[i] / 30 + h * h / 70);
		gau[5 * i + 5] -= lambda[i] * (cord[i] * cord[i] / 10 - h * h / 35);
	}
}


void Method::stiffness_matrix() {
	for (int i = 0; i < n; i++) {
		local_G(i);
	}
}


void Method::local_M(int i) {
	real h = cord[i + 1] - cord[i];
	real c = gamma[i] * h;
	if (type == 1) {
		mdi[i] += c * (cord[i] * cord[i] / 3 + h * cord[i] / 6 + h * h / 30);
		mdi[i + 1] += c * (cord[i] * cord[i] / 3 + h * cord[i] / 2 + h * h / 5);
		mau[i] += c * (cord[i] * cord[i] / 6 + h * cord[i] / 6 + h * h / 20);
	}
	else {
		//Диагональные элементы
		mdi[2 * i] += c * (13 * cord[i] * cord[i] / 35 + 6 * h * cord[i] / 35 + 19 * h * h / 630);
		mdi[2 * i + 1] += c * h * h * (cord[i] * cord[i] / 105 + h * cord[i] / 140 + h * h / 630);
		mdi[2 * i + 2] += c * (13 * cord[i] * cord[i] / 35 + 4 * h * cord[i] / 7 + 29 * h * h / 126);
		mdi[2 * i + 3] += c * h * h * (cord[i] * cord[i] / 105 + h * cord[i] / 84 + h * h / 252);

		//Элементы нижнего треугольника
		mau[5 * i] += c * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 30 + 17 * h * h / 2520);
		mau[5 * i + 1] += c * (9 * cord[i] * cord[i] / 70 + 9 * h * cord[i] / 70 + 23 * h * h / 630);
		mau[5 * i + 2] += c * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 30 + 5 * h * h / 504);
		mau[5 * i + 3] -= c * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 35 + 19 * h * h / 2520);
		mau[5 * i + 4] -= c * h * h * (cord[i] * cord[i] / 140 + h * cord[i] / 140 + h * h / 504);
		mau[5 * i + 5] -= c * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 14 + 13 * h * h / 504);
	}
}


void Method::mass_matrix() {
	for (int i = 0; i < n; i++) {
		local_M(i);
	}
}


void Method::local_b(int i, real t) {
	real h = cord[i + 1] - cord[i];
	if (type == 1) {
		b[i] += func(cord[i], i, t) * h * (cord[i] * cord[i] / 3 + h * cord[i] / 6 + h * h / 30) +
			func(cord[i + 1], i, t) * h * (cord[i] * cord[i] / 6 + h * cord[i] / 6 + h * h / 20);
		b[i + 1] += func(cord[i], i, t) * h * (cord[i] * cord[i] / 6 + h * cord[i] / 6 + h * h / 20) +
			func(cord[i + 1], i, t) * h * (cord[i] * cord[i] / 3 + h * cord[i] / 2 + h * h / 5);
	}
	else {
		real f1 = func(cord[i], i, t);
		real f2 = der_func(cord[i], i, t);
		real f3 = func(cord[i + 1], i, t);
		real f4 = der_func(cord[i + 1], i, t);

		b[2 * i] += f1 * h * (13 * cord[i] * cord[i] / 35 + 6 * h * cord[i] / 35 + 19 * h * h / 630) +
			f2 * h * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 30 + 17 * h * h / 2520) +
			f3 * h * (9 * cord[i] * cord[i] / 70 + 9 * h * cord[i] / 70 + 23 * h * h / 630) -
			f4 * h * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 35 + 19 * h * h / 2520);

		b[2 * i + 1] += f1 * h * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 30 + 17 * h * h / 2520) +
			f2 * h * h * h * (cord[i] * cord[i] / 105 + h * cord[i] / 140 + h * h / 630) +
			f3 * h * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 30 + 5 * h * h / 504) -
			f4 * h * h * h * (cord[i] * cord[i] / 140 + h * cord[i] / 140 + h * h / 504);

		b[2 * i + 2] += f1 * h * (9 * cord[i] * cord[i] / 70 + 9 * h * cord[i] / 70 + 23 * h * h / 630) +
			f2 * h * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 30 + 5 * h * h / 504) +
			f3 * h * (13 * cord[i] * cord[i] / 35 + 4 * h * cord[i] / 7 + 29 * h * h / 126) -
			f4 * h * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 14 + 13 * h * h / 504);

		b[2 * i + 3] -= f1 * h * h * (13 * cord[i] * cord[i] / 420 + h * cord[i] / 35 + 19 * h * h / 2520) +
			f2 * h * h * h * (cord[i] * cord[i] / 140 + h * cord[i] / 140 + h * h / 504) +
			f3 * h * h * (11 * cord[i] * cord[i] / 210 + h * cord[i] / 14 + 13 * h * h / 504) -
			f4 * h * h * h * (cord[i] * cord[i] / 105 + h * cord[i] / 84 + h * h / 252);
	}
}


void Method::right_part_vector(int time_index) {
	real dt = time_cord[time_index] - time_cord[time_index - 1];
	vector<real> m = multply_additive_component();

	for (int i = 0; i < n; i++)
		local_b(i, time_cord[time_index]);
	for(int i = 0; i < b.size(); i++)
		b[i] += m[i] / dt;
}


void Method::global_matrix(int time_index) {
	if (flag) return;
	real dt = time_cord[time_index] - time_cord[time_index - 1];

	for (int i = 0; i < di.size(); i++)
		di[i] = gdi[i] + mdi[i] / dt;

	for(int i = 0; i < au.size(); i++)
		au[i] = gau[i] + mau[i] / dt;
}


void Method::save_as_previous() {
	for (int i = 0; i < old_q.size(); i++)
		old_q[i] = new_q[i];
}


void Method::dirichlet_condition_left(real t) {
	di[0] = 1e+70;
	b[0] = 1e+70 * fun_u(cord[0], t);

	if (type == 2) {
		di[1] = 1e+70;
		b[1] = 1e+70 * derfun_u(cord[0], t);
	}
}


void Method::dirichlet_condition_right(real t) {
	int m;
	if (type == 1) m = n;
	else m = 2 * n;

	di[m] = 1e+70;
	b[m] = 1e+70 * fun_u(cord[n], t);

	if (type == 2) {
		di[m + 1] = 1e+70;
		b[m + 1] = 1e+70 * derfun_u(cord[n], t);
	}
}


void Method::neyman_condition_left(real t) {
	b[0] -= lambda[0] * derfun_u(cord[0], t) * cord[0] * cord[0];
}


void Method::neyman_condition_right(real t) {
	b[n] += lambda[n - 1] * derfun_u(cord[n], t) * cord[n] * cord[n];
}


void Method::boundary_conditions(real t) {
	if (flag) return;

	if (left == 2) neyman_condition_left(t);
	if (right == 2) neyman_condition_right(t);

	if (left == 1) dirichlet_condition_left(t);
	if (right == 1) dirichlet_condition_right(t);
}


void Method::LLt() {
	if (flag) return;

	int m;
	if (type == 1) m = n;
	else m = 2 * n + 1;

	for (int i = 0; i <= m; i++) {
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);

		real d = 0;
		for (int ij = i0; ij < i1; ij++, j++) {
			int j0 = ia[j];
			int j1 = ia[j + 1];
			int k, ik, jk;
			int jz = j - (j1 - j0);
			int iz = i - (i1 - i0);
			real s = 0;

			if (jz > iz) {
				jk = j0;
				ik = i0 + (jz - iz);
				k = jz;
			}
			else {
				ik = i0;
				jk = j0 + (iz - jz);
				k = iz;
			}

			for (; ik < ij; ik++, jk++) {
				s += au[ik] * au[jk];
			}
			au[ij] = (au[ij] - s) / di[j];
			d += au[ij] * au[ij];
		}
		di[i] -= d;
		di[i] = sqrt(di[i]);
	}
}


void Method::direct_reverse() {
	if (flag) return;

	int m;
	if (type == 1) m = n;
	else m = 2 * n + 1;

	real* z;
	z = new real[m + 1];
	for (int i = 0; i <= m; i++) z[i] = 0;

	for (int i = 0; i <= m; i++) {
		real s = 0;
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int m = i0; m < i1; j++, m++)
			s += z[j] * au[m];
		z[i] = b[i] - s;
		z[i] /= di[i];
	}

	for (int i = m; i >= 0; i--) {
		real s = z[i] / di[i];
		int i0 = ia[i];
		int i1 = ia[i + 1];
		int j = i - (i1 - i0);
		for (int k = i0; k < i1; k++, j++)
			z[j] -= au[k] * s;
		new_q[i] = s;
	}
}


void Method::output(real t, ofstream &out) {
	if (flag) return;

	cout << "Time = " << t << endl;
	out << "Time = " << t << endl;
	cout << "Relative error = " << setprecision(14) << relative_error(t) << endl;

	if (type == 1) {

		for (int i = 0; i <= n; i++) {
			cout << fixed;
			cout << setprecision(14) << "q[" << i + 1 << "] = " << new_q[i];
			cout << "    u_r[" << i + 1 << "] = " << fun_u(cord[i], t);
			cout << "    |q[" << i + 1 << "] - u_r[" << i + 1 << "]| = " << abs(new_q[i] - fun_u(cord[i], t)) << endl;
			out << fixed;
			out << setprecision(14) << new_q[i] << "     |     ";
			out << fun_u(cord[i], t) << "     |     ";
			out << abs(new_q[i] - fun_u(cord[i], t)) << endl;
		};
	}
	else {

		for (int i = 0, j = 0; i <= n; i++, j += 2) {
			cout << fixed;
			cout << setprecision(14) << "q[" << i + 1 << "] = " << new_q[j];
			cout << "    u_r[" << i + 1 << "] = " << fun_u(cord[i], t);
			cout << "    |q[" << i + 1 << "] - u_r[" << i + 1 << "]| = " << abs(new_q[j] - fun_u(cord[i], t)) << endl;
			out << fixed;
			out << setprecision(14) << new_q[j] << "     |     ";
			out << fun_u(cord[i], t) << "     |     ";
			out << abs(new_q[j] - fun_u(cord[i], t)) << endl;
		}
	}
	cout << endl;
	out << endl;
}


void Method::iteration_process() {
	input();
	matrix_portrait();
	stiffness_matrix();
	mass_matrix();

	if (type == 1) {
		for (int i = 0; i < old_q.size(); i++)
			old_q[i] = fun_u(cord[i], time_cord[0]);
	}
	else {
		for (int i = 0; i < old_q.size() / 2; i++) {
			old_q[2 * i] = fun_u(cord[i], time_cord[0]);
			old_q[2 * i + 1] = derfun_u(cord[i], time_cord[0]);
		}
	}

	ofstream out("output.txt");
	for (int i = 1; i < time; i++) {
		//сборка глобальной матрицы
		global_matrix(i);
		//сборка вектора правой части
		right_part_vector(i);
		//учет краевых
		boundary_conditions(time_cord[i]);
		//решение СЛАУ
		LLt();
		direct_reverse();
		//сохранение вектора весов
		save_as_previous();
		//Вывод результата
		output(time_cord[i], out);
		//чистка глобальной матрицы и вектора правой части
		matrix_clean();
	}
	out.close();
}


int main() {
	Method m;
	m.iteration_process();
}