//#pragma once
////一系列输出信息函数的重载和定义
//#include <iostream>
//#include <qstring.h>
//#include <vector>
//#include <Eigen/Dense>
//std::ostream& operator<<(std::ostream& os, const QString &qs);
//std::ostream& operator << (std::ostream& os, const Eigen::MatrixXd &m);
//
//template <class T>
//std::ostream& operator << (std::ostream& os, const std::vector<T> &v) {
//	for (auto itr : v) os << itr << " ";
//	return os;
//}
//
//template <class T>
//inline void dprinterror(const T &error) { std::cerr << error << std::endl; }
//
//inline void dprint() { std::cout << "\n"; }
//
//template <class T, class... Args>
//static void dprint(const T info, const Args... rest)
//{
//	std::cout << info << " ";
//	dprint(rest...);
//}
//
//template <class T, class... Args>
//static void dprintwithprecision(const int pr, const T info, const Args... rest)
//{
//	std::cout << std::setprecision(pr);
//	dprint(info, rest...);
//}
//
//#include <ctime>
//class timeRecorder
//{
//public:
//	timeRecorder() { timeKnots.push_back(clock()); }
//	timeRecorder(const timeRecorder& tr) = delete;
//	~timeRecorder() {}
//private:
//	std::vector<clock_t> timeKnots;
//public:
//	inline void tog(bool if_discard_before = false)
//	{
//		if (if_discard_before)
//			timeKnots.back() = clock();
//		else
//			timeKnots.push_back(clock());
//	}
//	inline void out(const std::string &info = "time:")
//	{
//		clock_t presentTime = clock();
//		dprint(info, presentTime - timeKnots.back(), "ms");
//		timeKnots.push_back(presentTime);
//	}
//	inline void sum(const std::string &info = "sum time:")
//	{
//		clock_t presentTime = clock();
//		dprint(info, presentTime - timeKnots.front(), "ms");
//		timeKnots.push_back(presentTime);
//	}
//};
