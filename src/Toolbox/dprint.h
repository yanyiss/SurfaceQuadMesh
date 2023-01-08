#pragma once
//一系列输出信息函数的重载和定义
#include <iostream>
#include <qstring.h>
#include <vector>
#include <map>
#include <Eigen/Dense>
std::ostream& operator<<(std::ostream& os, const QString &qs);
std::ostream& operator << (std::ostream& os, const Eigen::MatrixXd &m);

template <class T>
std::ostream& operator << (std::ostream& os, const std::vector<T> &v) {
	for (auto itr : v) os << itr << " ";
	return os;
}

template <class T>
inline void dprinterror(const T &error) { std::cerr << error << std::endl; }

inline void dprint() { std::cout << "\n"; }

template <class T, class... Args>
static void dprint(const T info, const Args... rest)
{
	std::cout << info << " ";
	dprint(rest...);
}

template <class T, class... Args>
static void dprintwithprecision(const int pr, const T info, const Args... rest)
{
	std::cout << std::setprecision(pr);
	dprint(info, rest...);
}

#include <ctime>
class timeRecorder
{
public:
	timeRecorder() { timeKnots.push_back(clock()); }
	timeRecorder(const timeRecorder& tr) = delete;
	~timeRecorder() {}
private:
	std::vector<clock_t> timeKnots;
	std::map<std::string, clock_t> timeMark;
public:
	inline void tog(bool if_discard_before = false)
	{
		if (if_discard_before)
			timeKnots.back() = clock();
		else
			timeKnots.push_back(clock());
	}
	void out(const std::string &info = "time:")
	{
		clock_t presentTime = clock();
		dprint(info, presentTime - timeKnots.back(), "ms");
		timeKnots.push_back(presentTime);
	}
	void sum(const std::string &info = "sum time:")
	{
		clock_t presentTime = clock();
		dprint(info, presentTime - timeKnots.front(), "ms");
		timeKnots.push_back(presentTime);
	}
	void begin(const std::string &info)
	{
		if (timeMark.find(info) == timeMark.end())
			timeMark.insert(std::make_pair(info, clock_t(0)));
		tog();
	}
	void end(const std::string &info)
	{
		clock_t presentTime = clock();
		timeMark[info] += presentTime - timeKnots.back();
		timeKnots.push_back(presentTime);
	}
	void mark(const std::string &info)
	{
		dprint(info, timeMark[info], "ms");
	}
};
