/**	\file pdbutils.hpp
 *
 *	\brief This file desperately needs comments.
 *
 */

#include "core.hpp"
#include "exceptions.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <set>

namespace tpp {

class BondMatrix {
private:
	std::set<int> curset;
	boost::numeric::ublas::matrix<bool> mtx;
	void rec_check(int i);
	bool check();
public:
	BondMatrix(const OpenBabel::OBMol&);
	~BondMatrix() {
	}
	const boost::numeric::ublas::matrix<bool> &operator *() {
		return mtx;
	}
	void print(std::ostream &os);
};

extern std::pair<boost::numeric::ublas::matrix<int>,
		boost::numeric::ublas::matrix<int> > top_neig(OpenBabel::OBMol &,
		boost::numeric::ublas::matrix<int>, int, int, int, int);
extern boost::numeric::ublas::vector<int> generate_long_tail(
		OpenBabel::OBMol &);
extern boost::numeric::ublas::vector<int> generate_long_tail1(
		OpenBabel::OBMol &);
extern AtomArray mol_renum(OpenBabel::OBMol &, AtomArray &,
		boost::numeric::ublas::vector<int>);
extern AtomArray mol_renum1(OpenBabel::OBMol &, AtomArray &,
		boost::numeric::ublas::vector<int>);
}

