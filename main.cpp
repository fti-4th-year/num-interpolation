#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <functional>
#include <vector>

struct vec2 {
	double x, y;
	vec2() = default;
	vec2(double x, double y) {
		this->x = x;
		this->y = y;
	}
};

class Polynomial {
protected:
	int size;
	std::vector<vec2> points;
	
	void setPoints(const vec2 *pts, int s) {
		points.resize(s);
		for(int i = 0; i < s; ++i) {
			points[i] = pts[i];
		}
		size = s;
	}
	
public:
	Polynomial(const vec2 *points, int size) {
		setPoints(points, size);
	}
	virtual ~Polynomial() = default;
	
	virtual double sample(double x) = 0;
};

class NewtonPolynomial : public Polynomial {
private:
	std::vector<double> factors;
	
	double getFactor(int first, int last) {
		if(first == last)
			return points[first].y;
		return (getFactor(first + 1, last) - getFactor(first, last - 1))/(points[last].x - points[first].x);
		
	}
	
	void findFactors() {
		factors.resize(size);
		for(int i = 0; i < size; ++i) {
			factors[i] = getFactor(0, i);
		}
	}
	
public:
	NewtonPolynomial(const vec2 *points, int size) 
	    : Polynomial(points, size) 
	{
		findFactors();
	}
	
	virtual double sample(double x) override {
		double sum = 0.0;
		double base = 1.0;
		for(int i = 0; i < size; ++i) {
			sum += factors[i]*base;
			base *= (x - points[i].x);
		}
		return sum;
	}
};

class LagrangePolynomial : public Polynomial {
private:
	std::vector<double> factors;
	
	void findFactors() {
		factors.resize(size);
		for(int i = 0; i < size; ++i) {
			double mul = points[i].y;
			for(int j = 0; j < size; ++j) {
				mul /= (j != i) ? (points[i].x - points[j].x) : 1.0;
			}
			factors[i] = mul;
		}
	}
	
public:
	LagrangePolynomial(const vec2 *points, int size) 
	    : Polynomial(points, size) 
	{
		findFactors();
	}
	
	virtual double sample(double x) override {
		double sum = 0.0;
		for(int i = 0; i < size; ++i) {
			double base = 1.0;
			for(int j = 0; j < size; ++j) {
				base *= (j != i) ? (x - points[j].x) : 1.0;
			}
			sum += factors[i]*base;
		}
		return sum;
	}
};

double majorant(const std::vector<vec2> &v, double x) {
	double fact = 1.0;
	for(int i = 2; i <= int(v.size()) + 1; ++i) {
		fact *= i;
	}
	double omega = 1.0;
	for(int i = 0; i < int(v.size()); ++i) {
		omega *= (x - v[i].x);
	}
	return omega/fact;
}

//#define DEP

int main(int argc, char *argv[]) {
	
	std::function<double(double)> func = [](double x) -> double {return sin(2.0*M_PI*x);};
	
#ifndef DEP
	
	FILE *points, *line, *error;
	points = fopen("points.txt", "w");
	line = fopen("line.txt", "w");
	error = fopen("error.txt", "w");
	
	std::vector<vec2> v;
	const int n = 6;
#ifdef RAND_PTS
	srand(0xABCDEF);
#endif
	for(int i = 0; i < n; ++i) {
		double x = (double) i/(n - 1);
#ifdef RAND_PTS
		double y = 2.0*((double) random())/RAND_MAX - 1.0;
#else
		double y = func(x);
#endif
		fprintf(points, "%.32lf %.32lf\n", x, y);
		v.push_back(vec2(x, y));
	}
	
	// NewtonPolynomial p(v.data(), v.size());
	LagrangePolynomial p(v.data(), v.size());
	
	//fprintf(line, "%lf %lf\n", v[0].x, v[0].y);
	const double extra = 0.2;
	const double step = 1e-3;
	for(double x = step - extra; x < 1.0 + extra; x += step) {
		double value = p.sample(x);
		fprintf(line, "%.32lf %.32lf\n", x, value);
		fprintf(error, "%.32lf %.32lf %lf\n", x, fabs(value - func(x)), fabs(pow(2*M_PI, n + 1)*majorant(v, x)));
	}
	//fprintf(line, "%lf %lf\n", v[v.size() - 1].x, v[v.size() - 1].y);
	
	fclose(error);
	fclose(line);
	fclose(points);
	
#else // DEP
	
	const int n = 48;
	const double step = 1e-3;
	
	FILE *errdep = fopen("errdep.txt", "w");
	for(int i = 2; i <= n; ++i) {
		std::vector<vec2> lv;
		for(int j = 0; j < i; ++j) {
			double x = (double) j/(i - 1);
			double y = func(x);
			lv.push_back(vec2(x, y));
		}
		
		//NewtonPolynomial lp(lv.data(), lv.size());
		LagrangePolynomial lp(lv.data(), lv.size());
		
		double max_err = 0.0;
		for(double x = step; x < 1.0; x += step) {
			double err = fabs(lp.sample(x) - func(x));
			if(err > max_err)
				max_err = err;
		}
		fprintf(errdep, "%d %.32lf\n", i, max_err);
	}
	fclose(errdep);
	
#endif // DEP
}
