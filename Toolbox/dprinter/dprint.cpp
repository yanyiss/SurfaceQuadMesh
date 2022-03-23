#include "dprint.h"
std::ostream& operator<<(std::ostream& os, const QString &qs) {
	os << qs.toLatin1().data();
	return os;
}

std::ostream& operator << (std::ostream& os, const Eigen::MatrixXd &m) {
	int r = m.rows();
	int c = m.cols();
	for (int i = 0; i < r; i++) {
		os << i;
		for (int j = 0; j < c; j++) {
			os << " " << m(i, j);
		}
		os << "\n";
	}
	return os;
}


