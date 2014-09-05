#include "global.hpp"
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <set>


namespace tpp {


  class BondMatrix {
    private:
      std::set<int> curset;
      boost::numeric::ublas::matrix<bool> mtx;
      void rec_check(int i);
      bool check() throw (t_exception);
    public:
      BondMatrix(const OBMol&) throw (t_exception);
      ~BondMatrix() {;}
      const boost::numeric::ublas::matrix<bool> &operator *() {
        return mtx;
      }
      void print (std::ostream &os);
  };

  extern std::pair<boost::numeric::ublas::matrix<int>, 
         boost::numeric::ublas::matrix<int> > top_neig
   (OBMol &, boost::numeric::ublas::matrix<int>, int, int, int, int) throw (t_exception);
  extern boost::numeric::ublas::vector<int> generate_long_tail(OBMol &) throw (t_exception);
  extern boost::numeric::ublas::vector<int> generate_long_tail1(OBMol &) throw (t_exception);
  extern t_atom_array mol_renum(OBMol &, t_atom_array &, boost::numeric::ublas::vector<int>) throw (t_exception);
  extern t_atom_array mol_renum1(OBMol &, t_atom_array &, boost::numeric::ublas::vector<int>) throw (t_exception);
}

