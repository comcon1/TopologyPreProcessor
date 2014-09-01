#include "global.hpp"
#include "calc.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "ublasadd.hpp"
#include "boost/numeric/bindings/lapack/gesv.hpp"
#include "boost/numeric/bindings/lapack/gesvd.hpp"
#include "hessian_accept.hpp"
#include "set"

namespace tpp {

  using ublas::slice;
  using ublas::matrix;
  using ublas::row;
  namespace lapack = boost::numeric::bindings::lapack;

static const double kbond = 973890.53;
static const double kang  = 2633.40;
static const double BOHR  = 0.52917706;

bool accept_hessian(t_topology &tp, const matrix<double> &hessian) throw (t_exception) {
  // simple correctness testing
  BOOST_REQUIRE( hessian.size1() == tp.atoms.size()*3 );
  BOOST_REQUIRE( hessian.size1() == hessian.size2() );

  // generate Jacobi matrix
  matrix<double> jacobi(hessian.size1(), hessian.size2());
  jacobi *= 0.0;
  int jacrow = 0;
  // map for exclude from jacobi
  map<string, vector<int> > coord_map;
  map<string, vector<double> > coord_values;
  // applying all bonds into new basis
  int _Counter = 0;
  for (t_top_map::nth_index<2>::type::iterator typit = tp.parameters.get<2>().begin();
       typit != tp.parameters.get<2>().end(); ++typit) 
  if (typit->type == TPP_TTYPE_BON) {
    if (typit->f == -1) { // bond is lack
    cout << "Processing " << typit->defname << endl;    
      coord_map.insert(std::pair<string, vector<int> >(typit->defname, vector<int>() ));
      coord_values.insert(std::pair<string, vector<double> >(typit->defname, vector<double>() ) );
    }
    for (t_top_array::nth_index<1>::type::iterator elit = tp.elements.get<1>().lower_bound(typit->defname);
        elit != tp.elements.get<1>().upper_bound(typit->defname); ++elit) {
      int i = elit->i, j = elit->j;
      t_point v1 = tp.atoms.find(i)->coord / BOHR, 
              v2 = tp.atoms.find(j)->coord / BOHR;
      int rn = svd_rank(jacobi);
      runtime.log_write("Adding into jacobi matrix bond (" 
          + lexical_cast<string>(i) + ", "
          + lexical_cast<string>(j) + ")\n");
      jacobi(jacrow, i*3-3) = __CALC_BOND(1, v1, v2);
      jacobi(jacrow, i*3-2) = __CALC_BOND(2, v1, v2);
      jacobi(jacrow, i*3-1) = __CALC_BOND(3, v1, v2);
      jacobi(jacrow, j*3-3) = __CALC_BOND(4, v1, v2);
      jacobi(jacrow, j*3-2) = __CALC_BOND(5, v1, v2);
      jacobi(jacrow, j*3-1) = __CALC_BOND(6, v1, v2);
      BOOST_CHECK(svd_rank(jacobi) == rn+1);
      if (coord_map.find(typit->defname) != coord_map.end())
        coord_map.find(typit->defname)->second.push_back(jacrow);
      if (coord_values.find(typit->defname) != coord_values.end())
        coord_values.find(typit->defname)->second.push_back(__CALC_BOND(0, v1, v2) );      
      jacrow++;
    } // end-for (tp.elements)
    _Counter++;
  } // if-BOND
  
  cout << jacrow << " bonds added to jacobi matrix." << endl;
  {
  // proceed single value decomposition:
  ublas::matrix<double, boost::numeric::ublas::column_major > _m(jacobi); 
  ublas::vector<double> S(_m.size1() < _m.size2() ? _m.size1() : _m.size2()); 
  lapack::gesvd(_m, S); 
  cout << S << endl;
  }

  
  // applaying bad angles into new basis
  for (t_top_map::nth_index<2>::type::iterator typit = tp.parameters.get<2>().begin();
       typit != tp.parameters.get<2>().end(); ++typit) 
  if (typit->type == TPP_TTYPE_ANG) {
    if (typit->f == -1) {
      cout << "Processing " << typit->defname;
      coord_map.insert(std::pair<string, vector<int> >(typit->defname, vector<int>() ));
      coord_values.insert(std::pair<string, vector<double> >(typit->defname, vector<double>() ) );
    } else continue;
    for (t_top_array::nth_index<1>::type::iterator elit = tp.elements.get<1>().lower_bound(typit->defname);
        elit != tp.elements.get<1>().upper_bound(typit->defname); ++elit) {
      int i = elit->i, j = elit->j, k = elit->k;
      t_point v1 = tp.atoms.find(i)->coord / BOHR, 
              v2 = tp.atoms.find(j)->coord / BOHR,
              v3 = tp.atoms.find(k)->coord / BOHR;
      int rn = svd_rank(jacobi);
      runtime.log_write("Adding into jacobi matrix angle (" 
          + lexical_cast<string>(i) + ", "
          + lexical_cast<string>(j) + ", "
          + lexical_cast<string>(k) + ")\n");
      jacobi(jacrow, i*3-3) = __CALC_ANG(1, v1, v2, v3);
      jacobi(jacrow, i*3-2) = __CALC_ANG(2, v1, v2, v3);
      jacobi(jacrow, i*3-1) = __CALC_ANG(3, v1, v2, v3);
      jacobi(jacrow, j*3-3) = __CALC_ANG(4, v1, v2, v3);
      jacobi(jacrow, j*3-2) = __CALC_ANG(5, v1, v2, v3);
      jacobi(jacrow, j*3-1) = __CALC_ANG(6, v1, v2, v3);
      jacobi(jacrow, k*3-3) = __CALC_ANG(7, v1, v2, v3);
      jacobi(jacrow, k*3-2) = __CALC_ANG(8, v1, v2, v3);
      jacobi(jacrow, k*3-1) = __CALC_ANG(9, v1, v2, v3);
      if (svd_rank(jacobi) != rn+1) {
        cout << 'x';
        row(jacobi, jacrow) *= 0.0;
      } else {
        cout << '.' << flush;

        if (coord_map.find(typit->defname) != coord_map.end())
          coord_map.find(typit->defname)->second.push_back(jacrow);
        if (coord_values.find(typit->defname) != coord_values.end())
          coord_values.find(typit->defname)->second.push_back(__CALC_ANG(0, v1, v2, v3) );      

        jacrow++;
      }
    } // end-for
    cout << endl;
  } // end-ANG
  // temporary write jacobi
  {
  ostringstream os;
  os << format("Current basis set is %1$3d, need to be updated to %2$3d.\n") % jacrow % hessian.size1();
  os << "Jacobi matrix after bond/angle completion:" << endl;
  for(int i=0; i<jacobi.size1(); ++i)
   os << row(jacobi,i) << endl;
  runtime.log_write(os.str());
  // proceed single value decomposition:
  ublas::matrix<double, boost::numeric::ublas::column_major > _m(jacobi); 
  ublas::vector<double> S(_m.size1() < _m.size2() ? _m.size1() : _m.size2()); 
  lapack::gesvd(_m, S); 
  cout << S << endl;
  }
  // *********************************************
  // Starting construction of orthogonal expansion
  //      according to comcon1 algorythm
  // *********************************************
  { // memory protection block
    do { // jacrow cycle --
    using std::set;
    set<int> goodcols;
    set<int> badcols;
    int crnk = 0;
    ublas::matrix<double> mt(0,0);
    for (int i=0; i<jacobi.size1(); ++i) {
      if (badcols.count(i) || goodcols.count(i)) continue;
      mt.resize(jacrow,mt.size2()+1);
      subrange(mt, 0, jacrow, mt.size2()-1, mt.size2()) = subrange(jacobi, 0, jacrow, i, i+1);

      if (svd_rank(mt) == crnk + 1) {
        goodcols.insert(i);
        crnk++;
        continue;
      } else {
        badcols.insert(i);
        mt.resize(jacrow,mt.size2()-1);
      }
    }
/*      // temporary output mt
      for (int xx=0; xx<mt.size1(); ++xx) {
       for (int yy=0; yy<mt.size2(); ++yy)
        cout << format("%1$6.3f ") % mt(xx,yy);
       cout << endl; 
      }
      cout << "-------------------------------------" << endl;
      // temporary output mt */
    BOOST_CHECK(mt.size1() == jacrow && mt.size2() == jacrow); // check for corectness
    matrix<double,ublas::column_major> col_(jacrow,1);
    col_ = ublas::zero_matrix<double>(jacrow,1);
    matrix<double,ublas::column_major> lins(mt);
    for (set<int>::iterator it = badcols.begin(); it != badcols.end(); ++it)
      col_ -= 1.0*subrange(jacobi,0,jacrow,*it,*it+1);
//    cout << "TEST COL 0: " << endl << col_ << endl;
    lapack::gesv(lins, col_); // solution now in col_
//    cout << "TEST COL 1: " << endl << col_ << endl;
    int ii = 0;
    for (int xx=0; xx<jacobi.size1(); xx++)
      jacobi(jacrow, xx) = goodcols.count(xx) ? col_(ii++,0) : 1.0;
// TEST FOR ORTHOGONALITY
    for (ii = 0; ii < jacrow; ++ii)
      BOOST_CHECK(fabs(ublas::inner_prod(ublas::row(jacobi,jacrow),ublas::row(jacobi,ii))) < 1e-10 );
    } while(++jacrow < jacobi.size1());
    // jacrow cycle --
  }

  {
  // proceed single value decomposition:
  ublas::matrix<double, boost::numeric::ublas::column_major > _m(jacobi); 
  ublas::vector<double> S(_m.size1() < _m.size2() ? _m.size1() : _m.size2()); 
  lapack::gesvd(_m, S); 
  cout << S << endl;
  }

// temporary write
  {
  ostringstream os;
  os << "Jacobi matrix after bond/angle completion:" << endl;
  for(int i=0; i<jacobi.size1(); ++i) {
   for (int j=0; j<jacobi.size2(); ++j)
    os << format("%1$6.3f ") % jacobi(i,j);
   os << endl;
  }
  os << format("DETERMINANT: %1$8.3e") % lu_det(jacobi);
  os << endl;
  runtime.log_write(os.str());
  }
  // JACOBIAN FINALLY CONVERTED
  // STARTING CONVERTING HESSIAN
  ublas::matrix<double> invjac(jacobi.size1(), jacobi.size2());
  svd_inv(jacobi, invjac); // invert jacobi matrix

  {
  ostringstream os;
  os << "Invert Jacobi matrix:" << endl;
  for(int i=0; i<jacobi.size1(); ++i) {
   for (int j=0; j<jacobi.size2(); ++j)
    os << format("%1$6.3f ") % invjac(i,j);
   os << endl;
  }
  os << format("DETERMINANT: %1$8.3e") % lu_det(invjac);
  os << endl;
  runtime.log_write(os.str());
  }  
  
  
  
  // cycle all lack parameters
  for (t_top_map::nth_index<2>::type::iterator typit = tp.parameters.get<2>().lower_bound(-1);
       typit != tp.parameters.get<2>().upper_bound(-1); ++typit) {
    runtime.log_write("Applying hessian for " + typit->defname + "\n");
    double sumh = 0.00, sumv = 0.00;
    BOOST_CHECK( coord_map.find(typit->defname)->second.size() ==
                 coord_values.find(typit->defname)->second.size() 
                );
    runtime.log_write("HESSE elements: ");
    for (vector<int>::iterator elit = coord_map.find(typit->defname)->second.begin();
        elit != coord_map.find(typit->defname)->second.end(); ++elit) {
      for (int i=0; i<invjac.size1(); ++i)
        for (int j=0; j<invjac.size2(); ++j)
           sumh += hessian(i,j)*invjac(i,*elit)*invjac(j,*elit);
      runtime.log_write(lexical_cast<string>(sumh)+" ");
    } // end for coord_map
    runtime.log_write("\n");    
    runtime.log_write("Values: ");
    for (vector<double>::iterator elit = coord_values.find(typit->defname)->second.begin();
      elit != coord_values.find(typit->defname)->second.end(); ++elit) {
      runtime.log_write(boost::lexical_cast<string>(*elit)+ " ");
      sumv += *elit;
    }
    runtime.log_write("\n");
    sumh /= coord_map.find(typit->defname)->second.size();
    sumv /= coord_values.find(typit->defname)->second.size();
    t_top_coord ttc;
    ttc.defname = typit->defname;
    ttc.type = typit->type;  
    ttc.f = -1;  
    switch (typit->type) {
        case TPP_TTYPE_BON:
          ttc.c0 = sumv * 0.052;
          ttc.c1 = sumh*kbond; 
          break;
        case TPP_TTYPE_ANG:
          ttc.c0 = sumv / M_PI * 180;
          ttc.c1 = sumh*kang;
          break;
        case TPP_TTYPE_IMPDIH:
          break;
    };  // end switch
    // correctly updating
    tp.parameters.get<2>().replace(typit, ttc);
  } // end.for (t_top_map)*/
  return 0;
} // end function hessian_accept

} //end.namespace

