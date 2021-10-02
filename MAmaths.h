#pragma once

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <array>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <assert.h>

#define PI      3.14159265359
#define TWO_PI  6.28318530718
#define HALF_PI 1.57079632679

namespace MAmaths {

#pragma region __HELPER
	static std::unordered_map<double, std::pair<double, double>> WrongDegrees
	{
		std::pair<double, std::pair<double, double>>(90.0,  std::pair<double, double>(1.0, 0.0)),
		std::pair<double, std::pair<double, double>>(180.0,  std::pair<double, double>(0.0, -1.0)),
		std::pair<double, std::pair<double, double>>(270.0,  std::pair<double, double>(-1.0, 0.0))
	};

	//Helper function declarations
	class Complex;
	std::string RemoveWhitespaces(const std::string& str);
	size_t GetLastOf(const std::string& str, const std::string& chars);
	void RightTrim(std::string& number, char c);
	size_t GetDecimalCount(double num);
	double GetGreatestCommonDivisor(double a, double b);
	double GetLeastCommonMultiple(double a, double b);
#pragma endregion

	struct Fraction
	{
	public:
		Fraction(double a = 0, double b = 1)
			: A(a / GetGreatestCommonDivisor(a, b)), B(b / GetGreatestCommonDivisor(a, b)) { }

		Fraction(const std::string& number)
		{
			if (number.find('/') != std::string::npos)
			{
				int pos = number.find('/');
				std::string a = number.substr(0, pos);
				std::string b = number.substr(pos + 1, number.size() - pos - 1);
				A = stod(a); B = stod(b);
			}
			else *this = ToFraction(stod(number));
		}

		double GetValue() const { return A / B; }
		double ToDecimal() const { return A / B; }

	public:
		Fraction operator+(const Fraction& other) const
		{
			double denom = GetLeastCommonMultiple(B, other.B);
			double numer = (A * denom / B) + (other.A * denom / other.B);
			double gcd = GetGreatestCommonDivisor(numer, denom);
			return Fraction(numer / gcd, denom / gcd);
		}
		Fraction operator+() const { return Fraction(A, B); }
		Fraction& operator+=(const Fraction& other)
		{
			*this = *this + other;
			return *this;
		}
		Fraction operator-(const Fraction& other) const
		{
			return *this + (-other);
		}
		Fraction operator-() const { return Fraction(-A, B); }
		Fraction& operator-=(const Fraction& other)
		{
			*this = *this - other;
			return *this;
		}
		Fraction operator*(const Fraction& other) const
		{
			double rA = A * other.A;
			double rB = B * other.B;
			double gcd = GetGreatestCommonDivisor(rA, rB);
			return Fraction(rA / gcd, rB / gcd);
		}
		Fraction& operator*=(const Fraction& other)
		{
			*this = *this * other;
			return *this;
		}
		Fraction& NotSimplifyingMultiplication(double value)
		{
			A *= value;
			B *= value;
			return *this;
		}
		Fraction operator/(const Fraction& other) const { return *this * other.GetReciprocal(); }
		Fraction& operator/=(const Fraction& other) { *this = *this / other; return *this; }
		bool operator==(const Fraction& other) const { return A == other.A && B == other.B; }
		bool operator==(double value) const { return *this == ToFraction(value); }
		bool operator==(int value) const { return *this == Fraction(value, 1); }
		bool operator!=(const Fraction& other) const { return A != other.A || B != other.B; }
		bool operator!=(double value) const { return *this != ToFraction(value); }
		bool operator!=(int value) const { return *this != Fraction(value, 1); }
		bool operator<(const Fraction& other) const { return ToDecimal() < other.ToDecimal(); }
		bool operator<(double value) const { return ToDecimal() < value; }
		bool operator<(int value) const { return ToDecimal() < value; }
		bool operator<=(const Fraction& other) const { return ToDecimal() <= other.ToDecimal(); }
		bool operator<=(double value) const { return ToDecimal() <= value; }
		bool operator<=(int value) const { return ToDecimal() <= value; }
		bool operator>(const Fraction& other) const { return ToDecimal() > other.ToDecimal(); }
		bool operator>(double value) const { return ToDecimal() > value; }
		bool operator>(int value) const { return ToDecimal() > value; }
		bool operator>=(const Fraction& other) const { return ToDecimal() >= other.ToDecimal(); }
		bool operator>=(double value) const { return ToDecimal() >= value; }
		bool operator>=(int value) const { return ToDecimal() >= value; }
		operator double() const { return this->ToDecimal(); }
		friend std::ostream& operator<<(std::ostream& stream, const Fraction& fraction)
		{
			return stream << fraction.A << '/' << fraction.B;
		}
		Fraction GetReciprocal() const
		{
			if (A == 0) return *this;
			else if (A < 0) return Fraction(-B, -A);
			else return Fraction(B, A);
		}
		std::string	ToString() const
		{
			std::stringstream ss;
			ss << std::setprecision(10) << std::noshowpoint;
			if (A == 0)
				ss << 0;
			else
			{
				if (B == 1) ss << A;
				else if (B == -1) ss << -A;
				else
				{
					if (B < 0)
						ss << -A << "/" << -B;
					else
						ss << A << "/" << B;
				}
			}
			return ss.str();
		}

	public:
		double A;
		double B;

	public:
		static Fraction ToFraction(double a)
		{
			size_t decimalCount = GetDecimalCount(a);

			double denom = pow(10, decimalCount);
			double numer = round(a * denom);
			double gcd = GetGreatestCommonDivisor(numer, denom);
			if (gcd == 0.0) gcd = 1.0;

			return Fraction(numer / gcd, denom / gcd);
		}
	};

	class Complex
	{
	public:
		Complex()
			: A(0), B(0) { }
		Complex(const Fraction& a, const Fraction& b = { 0, 1 })
			: A(a), B(b) { }
		Complex(const std::string& complex)
		{
			std::string number = RemoveWhitespaces(complex);
			size_t pos = GetLastOf(number, "+-");
			std::replace(number.begin(), number.end(), ',', '.');
			//Only real or imaginary part
			if (pos == 0 || pos == -1)
			{
				size_t iPos = GetLastOf(number, "i");
				//Only real part
				if (iPos == -1)
				{
					A = Fraction(number);
					B = { 0, 1 };
				}
				//Only imaginary part
				else
				{
					A = { 0, 1 };
					if (iPos == 0 || number[iPos - 1] == '+') B = { 1, 1 };
					else if (number[iPos - 1] == '-') B = { -1, 1 };
					else
					{
						if (number[0] == '+') B = Fraction(number.substr(1, iPos));
						if (number[0] == '-') B = -Fraction(number.substr(1, iPos));
						else B = Fraction(number.substr(0, iPos));
					}
				}
			}
			//Both real and imaginary parts
			else
			{
				std::string aStr = number.substr(0, pos);
				A = Fraction(aStr);

				size_t iPos = GetLastOf(number, "i");
				if (number[iPos - 1] == '+') B = { 1, 1 };
				else if (number[iPos - 1] == '-') B = { -1, 1 };
				else
				{
					std::string bStr = number.substr(pos, number.length() - pos - 1);
					B = Fraction(bStr);
				}
			}
		}

	public:
		std::vector<Complex> GetPower(Fraction& power) const
		{
			//Whole power
			if (power.B == 1)
			{
				std::vector<Complex> roots;
				if (power >= 0)
				{
					double r = sqrt(A * A + B * B);
					double alpha;
					if (A == 0)
					{
						if (B > 0) alpha = HALF_PI;
						else if (B == 0) alpha = 0.0;
						else alpha = -HALF_PI;
					}
					else
					{
						alpha = atan(B / A);
						if (A < 0) alpha += PI;
					}
					double rRoot = pow(r, power);
					double alphaRoot = alpha * power;
					double aRoot, bRoot;

					if (WrongDegrees.find(round(alphaRoot * 180.0 / PI)) != WrongDegrees.end())
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * WrongDegrees[round(alphaRoot * 180.0 / PI)].second),
							Fraction::ToFraction(rRoot * WrongDegrees[round(alphaRoot * 180.0 / PI)].first));
						return roots;
					}
					else
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * cos(alphaRoot)), Fraction::ToFraction(rRoot * sin(alphaRoot)));
						return roots;
					}
				}
				else
				{
					double r = sqrt(A * A + B * B);
					double alpha;
					power *= Fraction(-1, 1);
					if (A == 0)
					{
						if (B > 0) alpha = HALF_PI;
						else if (B == 0) alpha = 0.0;
						else alpha = -HALF_PI;
					}
					else
					{
						alpha = atan(B / A);
						if (A < 0) alpha += PI;
					}
					double rRoot = pow(r, power);
					double alphaRoot = alpha * power;
					double aRoot, bRoot;

					if (WrongDegrees.find(round(alphaRoot * 180.0 / PI)) != WrongDegrees.end())
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * WrongDegrees[round(alphaRoot * 180.0 / PI)].second).GetReciprocal(),
							Fraction::ToFraction(rRoot * WrongDegrees[round(alphaRoot * 180.0 / PI)].first).GetReciprocal());
						return roots;
					}
					else
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * cos(alphaRoot)).GetReciprocal(), Fraction::ToFraction(rRoot * sin(alphaRoot)).GetReciprocal());
						return roots;
					}
				}
			}
			//Fraction power
			else
			{
				std::vector<MAmaths::Complex> roots = GetRoots(power.B);
				for (uint32_t i = 0; i < roots.size(); i++)
					roots[i] = roots[i].GetPower(MAmaths::Fraction(power.A))[0];
				return roots;
			}
		}
		Complex GetComplexPower(const Complex& power) const
		{
			double r = sqrt(A * A + B * B);
			double alpha;
			if (A == 0)
			{
				if (B > 0) alpha = HALF_PI;
				else if (B == 0) alpha = 0.0;
				else alpha = -HALF_PI;
			}
			else
			{
				alpha = atan(B / A);
				if (A < 0) alpha += PI;
			}

			double a = log(r);
			double alphaResult = a * power.B + alpha * power.A;
			double rResult = exp(a * power.A - alpha * power.B);
			double aResult = rResult * cos(alphaResult);
			double bResult = rResult * sin(alphaResult);

			return Complex(Fraction::ToFraction(aResult), Fraction::ToFraction(bResult));
		}
		std::vector<Complex> GetRoots(const Fraction& n) const
		{
			std::vector<Complex> roots;
			if (n >= 0)
			{
				double power = 1.0 / (double)n;
				double r = sqrt(A * A + B * B);
				double alpha;
				if (A == 0)
				{
					if (B > 0) alpha = HALF_PI;
					else if (B == 0) alpha = 0.0;
					else alpha = -HALF_PI;
				}
				else
				{
					alpha = atan(B / A);
					if (A < 0) alpha += PI;
				}
				double rRoot = pow(r, power);
				double alphaRoot = alpha * power;
				double alphaI, aRoot, bRoot;

				for (uint32_t k = 0; k < n; k++)
				{
					alphaI = alphaRoot + TWO_PI * k * power;
					if (WrongDegrees.find(round(alphaI * 180.0 / PI)) != WrongDegrees.end())
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * WrongDegrees[round(alphaI * 180.0 / PI)].second),
							Fraction::ToFraction(rRoot * WrongDegrees[round(alphaI * 180.0 / PI)].first));
					}
					else
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * cos(alphaI)),
							Fraction::ToFraction(rRoot * sin(alphaI)));
					}
				}
			}
			else
			{
				double power = -1.0 / (double)n;
				double r = sqrt(A * A + B * B);
				double alpha;
				if (A == 0)
				{
					if (B > 0) alpha = HALF_PI;
					else if (B == 0) alpha = 0.0;
					else alpha = -HALF_PI;
				}
				else
				{
					alpha = atan(B / A);
					if (A < 0) alpha += PI;
				}
				double rRoot = pow(r, power);
				double alphaRoot = alpha * power;
				double alphaI, aRoot, bRoot;

				for (uint32_t k = 0; k < n; k++)
				{
					alphaI = alphaRoot + TWO_PI * k * power;
					if (WrongDegrees.find(round(alphaI * 180.0 / PI)) != WrongDegrees.end())
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * WrongDegrees[round(alphaI * 180.0 / PI)].second).GetReciprocal(),
							Fraction::ToFraction(rRoot * WrongDegrees[round(alphaI * 180.0 / PI)].first).GetReciprocal());
					}
					else
					{
						roots.emplace_back(Fraction::ToFraction(rRoot * cos(alphaI)).GetReciprocal(),
							Fraction::ToFraction(rRoot * sin(alphaI)).GetReciprocal());
					}
				}
			}
			return roots;
		}
		std::string ToString() const
		{
			if (B < 0)
			{
				std::string a = std::to_string(A.ToDecimal()); RightTrim(a, '0');
				std::string b = std::to_string(-B.ToDecimal()); RightTrim(b, '0');
				if (a[a.length() - 1] == '.') a.erase(a.length() - 1, 1);
				if (b[b.length() - 1] == '.') b.erase(b.length() - 1, 1);
				return a + " - " + b + "i";
			}
			else
			{
				std::string a = std::to_string(A.ToDecimal()); RightTrim(a, '0');
				std::string b = std::to_string(abs(B.ToDecimal())); RightTrim(b, '0');
				if (a[a.length() - 1] == '.') a.erase(a.length() - 1, 1);
				if (b[b.length() - 1] == '.') b.erase(b.length() - 1, 1);
				return a + " + " + b + "i";
			}
		}

	public:
		static std::array<Complex, 2> SolveQuadratic(const Complex& a, const Complex& b, const Complex& c)
		{
			std::vector<Complex> discriminantRoots = ((b * b) - (4 * a * c)).GetRoots(Fraction(2));

			MAmaths::Complex r1 = (-b + discriminantRoots[0]) / (2 * a);
			MAmaths::Complex r2 = (-b + discriminantRoots[1]) / (2 * a);
			return { r1, r2 };
		}
		static std::array<Complex, 3> SolveCubic(const Complex& a, const Complex& b, const Complex& c, const Complex& d)
		{
			std::array<Complex, 3> roots;

			//First root
			MAmaths::Complex q = (d / a) + (2 * b.GetPower(Fraction(3))[0]) / (27 * a.GetPower(Fraction(3))[0]) - (b * c) / (3 * a.GetPower(Fraction(2))[0]);
			MAmaths::Complex p = (c / a) - (b.GetPower(Fraction(2))[0]) / (3.0 * a.GetPower(Fraction(2))[0]);
			std::vector<Complex> z = (q.GetPower(Fraction(2))[0] + (4.0 * p.GetPower(Fraction(3))[0] / 27.0)).GetRoots(Fraction(2)); z[0] /= 2.0; z[1] /= 2.0;
			Complex A = (-q / 2 + z[0]).GetRoots(3)[0]; //Choose the first cube root
			std::vector<Complex> BRoots = (-q / 2 + z[1]).GetRoots(3); //Get all cube roots of B
			//Check all cube roots of B whether it works with A (one of the will always work)
			for (uint32_t i = 0; i < 3; i++)
			{
				Complex temp = 3 * A * BRoots[i];
				if (Complex(-round(p.A * Fraction(1000)), -round(p.B * Fraction(1000))) == Complex(round(temp.A * Fraction(1000)), round(temp.B * Fraction(1000))))
				{
					roots[0] = (A + BRoots[i]) - (b / (3 * a)); //y - b/3a
					break;
				}
			}

			//Two roots from quadratic
			Complex a1 = a;
			Complex b1 = b + roots[0] * a1;
			Complex c1 = c + roots[0] * b1;
			std::array<MAmaths::Complex, 2> quadraticRoots = MAmaths::Complex::SolveQuadratic(a1, b1, c1);
			roots[1] = quadraticRoots[0];
			roots[2] = quadraticRoots[1];

			return roots;
		}
		static std::array<Complex, 4> SolveQuartic(Complex& a, Complex& b, Complex& c, Complex& d, Complex& e)
		{
			std::array<Complex, 4> roots;

			//Make the main coefficient 1
			b /= a;
			c /= a;
			d /= a;
			e /= a;
			a.A = 1; a.B = 0;

			//Calculate helper variables
			MAmaths::Complex p = (-3.0 * b.GetPower(Fraction(2))[0] / 8.0) + (c);
			MAmaths::Complex q = (b.GetPower(Fraction(3))[0] / 8.0) - (b * c / 2.0) + d;
			MAmaths::Complex r = (-3.0 * b.GetPower(Fraction(4))[0] / 256.0) + (b.GetPower(Fraction(2))[0] * c / 16.0) - (d * b / 4.0) + (e);

			//lSquareCubic contains the three possible values for l^2 (turns out that all give the same roots in the end)
			std::array<MAmaths::Complex, 3> lSquareCubic = MAmaths::Complex::SolveCubic({ 1 }, (2 * p), (p.GetPower(Fraction(2))[0] - 4 * r), (-1 * q.GetPower(Fraction(2))[0]));
			MAmaths::Complex l = lSquareCubic[0].GetRoots(Fraction(2))[0]; //The value of l

			//Calculate n and m using p, q and l
			MAmaths::Complex n = (p + l.GetPower(Fraction(2))[0] + (q / l)) / 2.0;
			MAmaths::Complex m = (p + l.GetPower(Fraction(2))[0] - (q / l)) / 2.0;

			//Get the roots from both quadratic equations
			std::array<MAmaths::Complex, 2> roots1 = MAmaths::Complex::SolveQuadratic({ 1 }, l, m); for (uint32_t i = 0; i < roots1.size(); i++) roots1[i] -= (b / 4.0);
			std::array<MAmaths::Complex, 2> roots2 = MAmaths::Complex::SolveQuadratic({ 1 }, -l, n); for (uint32_t i = 0; i < roots2.size(); i++) roots2[i] -= (b / 4.0);

			//Fill the roots array
			roots[0] = roots1[0];
			roots[1] = roots1[1];
			roots[2] = roots2[0];
			roots[3] = roots2[1];

			return roots;
		}

	public:
		Complex operator+(const Complex& other) const { return Complex(A + other.A, B + other.B); }
		Complex operator+(const Fraction& value) const { return Complex(A + value, B); }
		Complex operator+() const { return Complex(A, B); }
		Complex& operator+=(const Complex& other) { A += other.A; B += other.B; return *this; }
		Complex& operator+=(const Fraction& value) { A += value; B += value; return *this; }
		friend Complex operator+(const Fraction& value, const Complex& complex) { return complex + value; }
		Complex operator-(const Complex& other) const { return Complex(A - other.A, B - other.B); }
		Complex operator-(const Fraction& value) const { return Complex(A - value, B); }
		Complex operator-() const { return Complex(-A, -B); }
		Complex& operator-=(const Complex& other) { A -= other.A; B -= other.B; return *this; }
		Complex& operator-=(const Fraction& value) { A -= value; B -= value; return *this; }
		friend Complex operator-(const Fraction& value, const Complex& complex) { return complex - value; }
		Complex operator*(const Complex& other) const { return Complex(A * other.A - B * other.B, A * other.B + B * other.A); }
		Complex operator*(const Fraction& value) const { return Complex(A * value, B * value); }
		Complex& operator*=(const Complex& other) { *this = *this * other; return *this; }
		Complex& operator*=(const Fraction& value) { A* value; B* value; return *this; }
		friend Complex operator*(const Fraction& value, const Complex& complex) { return complex * value; }
		Complex operator/(const Complex& other) const
		{
			Complex otherConjugate(other.A, -other.B);
			Complex numerator = (*this) * otherConjugate;
			Complex denominator = other * otherConjugate;
			return Complex(numerator.A / denominator.A, numerator.B / denominator.A);
		}
		Complex operator/(const Fraction& value) const { return Complex(A / value, B / value); }
		Complex& operator/=(const Complex& other) { *this = *this / other; return *this; }
		Complex& operator/=(const Fraction& value) { A /= value; B /= value; return *this; }
		friend Complex operator/(const Fraction& value, const Complex& complex)
		{
			Complex complexConjugate(complex.A, -complex.B);
			Complex numerator = value * complexConjugate;
			Complex denominator = complex * complexConjugate;
			return Complex(numerator.A / denominator.A, numerator.B / denominator.A);
		}
		bool operator==(const Complex& other) const { return A == other.A && B == other.B; }
		bool operator==(const Fraction& value) const { return A == value && B == 0; }
		friend bool operator==(const Fraction& value, Complex complex) { return value == complex.A && complex.B == 0; }
		friend std::ostream& operator<<(std::ostream& stream, const Complex& complex)
		{
			stream << complex.ToString();
			return stream;
		}

	public:
		Fraction A;
		Fraction B;
	};

	class Matrix
	{
	public:
		class Row
		{
		public:
			Row(uint32_t length = 0)
				: m_Length(length)
			{
				m_Data.resize(length);
			}
			Row(const std::vector<Fraction>& data)
				: m_Length(data.size()), m_Data(data) { }
			Row(const Row& other)
				: m_Length(other.m_Length), m_Data(other.m_Data) { }

		public:
			uint32_t GetSize() const { return m_Length; }
			std::vector<Fraction>& GetData() { return m_Data; }
			const std::vector<Fraction>& GetData() const { return m_Data; }

		public:
			operator std::vector<Fraction>& () { return m_Data; }
			operator const std::vector<Fraction>& () const { return m_Data; }

			Row operator+(const Row& other) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] + other[i];;
				return result;
			}
			Row operator+(Fraction scalar) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] + scalar;
				return result;
			}
			Row& operator+=(const Row& other)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] += other[i];
				return *this;
			}
			Row& operator+=(Fraction scalar)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] += scalar;
				return *this;
			}
			Row operator-(const Row& other) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] - other[i];
				return result;
			}
			Row operator-(Fraction scalar) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] - scalar;
				return result;
			}
			Row& operator-=(const Row& other)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] -= other[i];
				return *this;
			}
			Row& operator-=(Fraction scalar)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] -= scalar;
				return *this;
			}
			Row operator*(const Row& other) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] * other[i];
				return result;
			}
			Row operator*(Fraction scalar) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] * scalar;
				return result;
			}
			Row& operator*=(const Row& other)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] *= other[i];
				return *this;
			}
			Row& operator*=(Fraction scalar)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] *= scalar;
				return *this;
			}
			Row operator/(const Row& other) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] / other[i];
				return result;
			}
			Row operator/(Fraction scalar) const
			{
				Row result(m_Length);
				for (int i = 0; i < m_Length; i++)
					result[i] = m_Data[i] / scalar;
				return result;
			}
			Row& operator/=(const Row& other)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] /= other[i];
				return *this;
			}
			Row& operator/=(Fraction scalar)
			{
				for (int i = 0; i < m_Length; i++)
					m_Data[i] /= scalar;
				return *this;
			}
			Fraction& operator[](uint32_t idx)
			{
				assert(idx < m_Length, "Index out of range on row");
				return m_Data[idx];
			}
			const Fraction& operator[](uint32_t idx) const
			{
				assert(idx < m_Length, "Index out of range on row");
				return m_Data[idx];
			}

		private:
			uint32_t m_Length = -1;
			std::vector<Fraction> m_Data;
		};

	public:
		Matrix(uint32_t size)
			: M(size), N(size)
		{
			m_Data.reserve(M);
			for (uint32_t i = 0; i < M; i++)
				m_Data.emplace_back(N);
		}
		Matrix(uint32_t rows, uint32_t columns)
			: M(rows), N(columns)
		{
			m_Data.reserve(M);
			for (uint32_t i = 0; i < M; i++)
				m_Data.emplace_back(N);
		}
		Matrix(const std::vector<Row>& data)
			: M(data.size()), N(data[0].GetSize()), m_Data(data) { }
		Matrix(const Matrix& other)
			: M(other.M), N(other.N)
		{
			m_Data = other.m_Data;
		}

	public:
		uint32_t GetRowCount() const { return M; }
		uint32_t GetColumnCount() const { return N; }
		uint32_t GetM() const { return M; }
		uint32_t GetN() const { return N; }

	public:
		Matrix operator+(const Matrix& other)
		{
			Matrix result(M, N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] + other[m][n];
				}
			}
			return result;
		}
		Matrix& operator+=(const Matrix& other)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] += other[m][n];
				}
			}
			return *this;
		}
		Matrix operator+(const Fraction& scaler)
		{
			Matrix result(M, N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] + scaler;
				}
			}
			return result;
		}
		Matrix& operator+=(const Fraction& scaler)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] += scaler;
				}
			}
			return *this;
		}
		Matrix operator-(const Matrix& other)
		{
			Matrix result(M, N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] - other[m][n];
				}
			}
			return result;
		}
		Matrix& operator-=(const Matrix& other)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] -= other[m][n];
				}
			}
			return *this;
		}
		Matrix operator-(const Fraction& scaler)
		{
			Matrix result(M, N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] - scaler;
				}
			}
			return result;
		}
		Matrix& operator-=(const Fraction& scaler)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] -= scaler;
				}
			}
			return *this;
		}
		Matrix operator*(const Matrix& other)
		{
			assert(N == other.M, "Cannot multiply the matrices");
			Matrix result(M, other.N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < other.N; n++)
				{
					for (uint32_t p = 0; p < N; p++)
					{
						result[m][n] += m_Data[m][p] * other.m_Data[p][n];
					}
				}
			}
			return result;
		}
		Matrix operator*(const Fraction& scaler)
		{
			Matrix result(M * N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] * scaler;
				}
			}
			return result;
		}
		Matrix& operator*=(const Fraction& scaler)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] *= scaler;
				}
			}
			return *this;
		}
		Matrix operator/(const Fraction& scaler)
		{
			Matrix result(M * N);
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					result[m][n] = m_Data[m][n] / scaler;
				}
			}
			return result;
		}
		Matrix& operator/=(const Fraction& scaler)
		{
			for (uint32_t m = 0; m < M; m++)
			{
				for (uint32_t n = 0; n < N; n++)
				{
					m_Data[m][n] /= scaler;
				}
			}
			return *this;
		}
		Row& operator[](uint32_t row)
		{
			assert(row < M, "Index out of range on matrix");
			return m_Data[row];
		}
		const Row& operator[](uint32_t row) const
		{
			assert(row < M, "Index out of range on matrix");
			return m_Data[row];
		}
		void SetData(const std::vector<Row>& data) { m_Data = data; }
		std::string GetString()
		{
			std::string str = "{";
			for (int m = 0; m < M; m++)
			{
				str += "{";
				for (int n = 0; n < N; n++)
				{
					str += m_Data[m][n].ToString();
					if (n < N - 1) str += ",";
				}
				str += "}";
				if (m < M - 1) str += ",";
			}
			str += "}";
			return str;
		}
		static Matrix GetFromString(const std::string& str)
		{
			std::string temp = str;
			RemoveWhitespaces(temp);
			temp.erase(0, 1);
			temp.erase(temp.size() - 1, 1);

			std::vector<Row> rows;
			int m = -1;
			int n = 0;
			int colCount = 0;
			for (int i = temp.find('{'); i < temp.size(); i++)
			{
				if (temp[i] == ',') colCount++;
				else if (temp[i] == '}')
				{
					colCount++;
					break;
				}
			}

			int i = 0;
			std::string num = "";
			while (++i < temp.size())
			{
				rows.emplace_back(colCount);
				m++;
				n = 0;
				num = "";
				while (temp[i] != '}')
				{
					if (temp[i] == ',')
					{
						if (num.find('/') != std::string::npos)
						{
							int pos = num.find('/');
							std::string a = num.substr(0, pos);
							std::string b = num.substr(pos + 1, num.size() - pos - 1);
							rows[m][n] = Fraction(stod(a), stod(b));
						}
						else rows[m][n] = Fraction::ToFraction(stod(num));
						num = ""; n++;
					}
					else num += temp[i];
					i++;
				}
				if (num.find('/') != std::string::npos)
				{
					int pos = num.find('/');
					std::string a = num.substr(0, pos);
					std::string b = num.substr(pos + 1, num.size() - pos - 1);
					rows[m][n] = Fraction(stod(a), stod(b));
				}
				else rows[m][n] = Fraction::ToFraction(stod(num));
				i += 2;
			}
			if (num.find('/') != std::string::npos)
			{
				int pos = num.find('/');
				std::string a = num.substr(0, pos);
				std::string b = num.substr(pos + 1, num.size() - pos - 1);
				rows[m][n] = Fraction(stod(a), stod(b));
			}
			else rows[m][n] = Fraction::ToFraction(stod(num));

			return { rows };
		}
		Matrix GetTransposed()
		{
			Matrix result(N, M);

			for (int m = 0; m < N; m++)
				for (int n = 0; n < M; n++)
					result[m][n] = m_Data[n][m];

			return result;
		}
		Fraction GetDeterminant()
		{
			assert(m_Data.size() == m_Data[0].GetSize(), "Determinant does not exist, because this is not a square matrix");

			if (m_Data.size() == 1) return m_Data[0][0];
			else if (m_Data.size() == 2) return m_Data[0][0] * m_Data[1][1] - m_Data[0][1] * m_Data[1][0];

			Fraction determinant(0, 1);
			for (int i = 0; i < m_Data.size(); i++)
			{
				Matrix subMatrix = __GetSubMatrix(0, i);
				determinant += (i & 1 ? -1 : 1) * m_Data[0][i] * subMatrix.GetDeterminant();
			}
			return determinant;
		}
		Fraction GetDeterminant2()
		{
			Matrix matrix(m_Data);

			//Make everything 0 below the main diagonal
			for (int n = 0; n < N - 1; n++)
			{
				for (int m = n + 1; m < M; m++)
				{
					if (matrix[m][n] == 0) continue;
					Fraction mul = matrix[n][n] / matrix[m][n];
					matrix[m] -= matrix[n] / mul;
				}
			}
			return __GetMainDiagonalMultiple(matrix);
		}
		Matrix GetInverse()
		{
			Fraction determinant = GetDeterminant();
			assert(determinant.A != 0, "Inverse does not exist, because the determinant is 0");

			//Setup the matrices
			Matrix matrix(m_Data);
			Matrix inverse(M, N);
			for (int m = 0; m < M; m++)
			{
				for (int n = 0; n < N; n++)
				{
					if (m == n) inverse[m][n] = Fraction(1, 1);
					else inverse[m][n] = Fraction(0, 1);
				}
			}

			//Make everything 0 below the main diagonal
			for (int n = 0; n < matrix.GetN() - 1; n++)
			{
				for (int m = n + 1; m < matrix.GetM(); m++)
				{
					if (matrix[m][n] == 0) continue;
					Fraction mul = matrix[n][n] / matrix[m][n];
					matrix[m] *= mul;
					inverse[m] *= mul;
					matrix[m] -= matrix[n];
					inverse[m] -= inverse[n];
				}
				inverse[n] /= matrix[n][n];
				matrix[n] /= matrix[n][n];
			}
			inverse[matrix.GetM() - 1] /= matrix[matrix.GetM() - 1][matrix.GetN() - 1];
			matrix[matrix.GetM() - 1] /= matrix[matrix.GetM() - 1][matrix.GetN() - 1];

			//Make everything 0 above the main diagonal
			for (int n = matrix.GetN() - 1; n >= 0; n--)
			{
				inverse[n] /= matrix[n][n];
				matrix[n] /= matrix[n][n];

				for (int m = n - 1; m >= 0; m--)
				{
					if (matrix[m][n] == 0) continue;
					Fraction mul = matrix[n][n] / matrix[m][n];
					matrix[m] *= mul;
					inverse[m] *= mul;
					matrix[m] -= matrix[n];
					inverse[m] -= inverse[n];
				}
			}

			if (inverse.GetM() == 1) return inverse;

			double lcm = determinant.ToDecimal();

			for (int m = 0; m < M; m++)
			{
				for (int n = 0; n < N; n++)
				{
					double mul = lcm / inverse[m][n].B;
					inverse[m][n].NotSimplifyingMultiplication(mul);
				}
			}

			return inverse;
		}
		Matrix GetInverse2()
		{
			Fraction determinant = GetDeterminant();
			assert(determinant.A != 0, "Inverse does not exist, because the determinant is 0");

			Fraction determinantReciprocal(1, GetDeterminant());
			Matrix inverse(N, M);

			for (int m = 0; m < M; m++)
			{
				for (int n = 0; n < N; n++)
				{
					Fraction r = __GetSubMatrix(n, m).GetDeterminant() * determinantReciprocal * ((m + n) & 1 ? MAmaths::Fraction(-1, 1) : MAmaths::Fraction(1, 1));
					Fraction mul = determinant.A / r.B;
					inverse[m][n] = r.NotSimplifyingMultiplication(mul);
				}
			}

			return inverse;
		}

	private:
		Matrix __GetSubMatrix(uint32_t row, uint32_t column)
		{
			Matrix subMatrix(m_Data.size() - 1);
			int rowOffs = 0, colOffs = 0;
			for (int m = 0; m < m_Data.size(); m++)
			{
				if (m == row)
				{
					rowOffs = -1;
					continue;
				}
				colOffs = 0;
				for (int n = 0; n < m_Data.size(); n++)
				{
					if (n == column)
					{
						colOffs = -1;
						continue;
					}
					subMatrix[m + rowOffs][n + colOffs] = m_Data[m][n];
				}
			}
			return subMatrix;
		}
		void __SwapLines(Matrix& matrix, int m1, int m2)
		{
			Row temp = matrix.m_Data[m2];
			matrix.m_Data[m2] = matrix.m_Data[m1];
			matrix.m_Data[m1] = temp;
		}
		Fraction __GetCommonMultiple(std::vector<Fraction> numbers)
		{
			Fraction max = numbers[0];
			for (int i = 1; i < numbers.size(); i++)
			{
				if (numbers[i] > max)
					max = numbers[i];
			}
			max.B /= abs(max.B);
			return max;
		}
		Matrix __MatrixMix(const Matrix& m1, const Matrix& m2)
		{
			Matrix result(m1.GetM(), m1.GetN() + m2.GetN());
			for (int m = 0; m < result.GetM(); m++)
			{
				for (int n = 0; n < result.GetN(); n++)
				{
					if (n < m1.GetN()) result[m][n] = m1[m][n];
					else result[m][n] = m2[m][n % m1.GetN()];
				}
			}
			return result;
		}
		Fraction __GetMainDiagonalMultiple(const Matrix& matrix)
		{
			Fraction multiple(1, 1);
			for (int m = 0; m < matrix.GetM(); m++)
				multiple *= matrix[m][m];
			return multiple;
		}

	private:
		uint32_t M;
		uint32_t N;

	private:
		std::vector<Row> m_Data;
	};

	//Function definitions
	std::string RemoveWhitespaces(const std::string& str)
	{
		std::string num = str;
		std::string::iterator end_pos = std::remove(num.begin(), num.end(), ' ');
		num.erase(end_pos, num.end());
		return num;
	}
	size_t GetLastOf(const std::string& str, const std::string& chars)
	{
		for (int i = str.length() - 1; i >= 0; i--)
		{
			for (int j = 0; j < chars.length(); j++)
				if (str[i] == chars[j]) return i;
		}
		return -1;
	}
	void RightTrim(std::string& number, char c)
	{
		int index;
		for (index = number.length() - 1; index >= 0; index--)
			if (number[index] != c) break;
		number.erase(index + 1, number.length() - index - 1);
	}
	size_t GetDecimalCount(double num)
	{
		std::string number = std::to_string(num);
		size_t count = 0;
		bool was = false;
		for (uint32_t i = number.length() - 1; i >= 0; i--)
		{
			if (number[i] == '.') break;
			else if (number[i] != '0' || was)
			{
				count++;
				was = true;
			}
		}
		return count;
	}
	double GetGreatestCommonDivisor(double a, double b)
	{
		if (b == 0)
			return a;
		return GetGreatestCommonDivisor(b, std::fmod(a, b));
	}
	double GetLeastCommonMultiple(double a, double b)
	{
		return a / GetGreatestCommonDivisor(a, b) * b;
	}
}