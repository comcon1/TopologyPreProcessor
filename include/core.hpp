/** \file core.hpp
 *
 *  \brief This file contains core definitions for TPP project.
 *
 *  Those include atoms, topologies and their respective containers.
 *
 */

#ifndef TPP_CORE_HEADER
#define TPP_CORE_HEADER

#include <string>

#include <boost/array.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/sequenced_index.hpp>

#include <boost/numeric/ublas/vector.hpp>

#include <openbabel/mol.h>

// boost serialization features
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

///  configurational definitions
#define TPP_GMX_EXTERNAL
#define TPP_MAX_ATOM_NUM 200
#define TPP_MAX_FREQ_NUM 600
#define TPP_MAX_BONDS 4
#define TPP_MAX_LONGTAIL_DEEP 10

/// definitions of atom definer
#define TPP_SMART_COEFF 100
#define TPP_ZNUC_COEFF  100000
#define TPP_BOND_COEFF  1
#define TPP_ANGLE_COEFF 2
#define TPP_DIHED_COEFF 3

/** \namespace tpp is used for scoping TPP specific variables.
 *
 */
namespace tpp {

  // namespace aliases

  namespace bub = boost::numeric::ublas;
  namespace bmi = boost::multi_index;

  // typedef aliases

  typedef bub::bounded_vector<double, 3> Point; //!< type of space coordinate
        typedef unsigned TppIndex; //!< alias for type of atom indexing in molecule array

  /**
   *  \brief This class describes atom properties in the topology or
   *  molecular structure.
   *
   */
  struct Atom {
    TppIndex index; //!< place in vector array
    TppIndex oldindex; //!< index in its initial array (make sense when renumbering)
    std::string atom_type; //!< non-bonded atom type (e.g. opls_136)
    std::string atom_type2; //!< bonded atom type (e.g. CT)
    std::string atom_name; //!< atom name unique for the residue
    std::string old_aname; //!< atom name in the original molecule (make sense when renumbering)
    unsigned char ncharge; //!< nuclear charge
    std::string res_name; //!< residue name
    std::string comment; //!< free string
    unsigned char mol_id; //!< TODO: i do not remember if we really need this field ..
    Point coord; //!< atom coordinates
    boost::array<TppIndex, TPP_MAX_BONDS> connects; //!< array of atom connectivity
    unsigned char num_connects; //!< atom valency
    double charge; //!< atom partial charge
    double mass; //!< atom mass (in a.u.)
                TppIndex c_gnr; //!< charge group!
    std::string qmname; //!< This surely means something!
  };

  /** \brief Container of atoms indexed by index, name or charge group.
    *
    *   It is important to locate atom by name and by index. Locating by
    *   charge group if questionable (may be should be removed).
    */
  typedef bmi::multi_index_container<Atom,
      bmi::indexed_by<bmi::ordered_unique<bmi::member<Atom, TppIndex, &Atom::index> >,
                        // key index (like array)
      bmi::ordered_non_unique<bmi::member<Atom, std::string, &Atom::atom_name> >,
                        // key by name
      bmi::ordered_non_unique<bmi::member<Atom, TppIndex, &Atom::c_gnr> >
                        // key by charge group
      > > AtomArray;

  /**
   * \brief This enumeration specifies format of the input data.
   *
   * Input data can be in coordinate format (like PDB) or can be read
   * from output of QM programs (GAMOPT).
   */
  typedef enum {
    TPP_IF_PDB = 0,       //!< Input in Protein Data Bank (PDB) format
    TPP_IF_GRO = 1,       //!< Input in GROMOS short format (GRO)
    TPP_IF_G96 = 2,       //!< Input in GROMOS "double precision" format (G96)
    TPP_IF_GAMOPT = 3,    //!< GAMESS optimization output as an input
    TPP_IF_GAMHESS = 4,   //!< GAMESS hessian output as an input
    TPP_IF_GAMSP = 5      //!< GAMESS singlepoint output as an input
  } InputFormat;

  /// Returns human-radable description of the input format.
  std::string in_fmt_descr( InputFormat fmt);

  /**
   * \brief This enumeration the output format
   *
   *  Note that OutputFormat is used for structure only.
   *  Topology output is not set by this enum.
   */
  typedef enum {
    TPP_OF_PDB = 0,       //!< Output in Protein Data Bank (PDB) format
    TPP_OF_GRO = 2,       //!< Output in GROMOS short format (GRO)
    TPP_OF_G96 = 3,       //!< Output in GROMOS "double precision" format (G96)
    TPP_OF_GAMIN = 4      //!< Output in GAMESS input format for further QM
  } OutputFormat;

  /// Returns human-radable description of the output format.
  std::string out_fmt_descr(OutputFormat fmt);

  /// Type of internal coordinate
  typedef enum {
    TPP_IC_BON  = 1,      //!< Valence bond
                TPP_IC_ANG  = 2,      //!< Valence angle
                TPP_IC_DIH  = 3,      //!< Dihedral angle
                TPP_IC_CHIR = 4       //!< TODO: what it is? Chirality definition?
  } IntCoordType;

  /** \brief Separate internal coordinate includes
   *
   * Includes information about atoms forming coordinate, type of the
   * coordinate and reference to corresponding tpp::TopCoord element
   * associated with this coordinate.
   */
  typedef struct {
    IntCoordType type;    //!< type of the coordinate
    TppIndex i,  //!< atom index No.1 (for e.g. bond type only i and j make sense)
                         j,  //!< atom index No.2
                         k,  //!< atom index No.3
                         l;  //!< atom index No.4
    std::string defname;  //!< reference to top::TopCoord element
  } IntCoord;

        /// Array of internal coordinates
  typedef std::vector<IntCoord> InternalsArray;

  /** \brief Parameters of potential function associated with internal coordinate
   *
   * Here we do not use all possible potentials but only that which we
   * operate.
   * TODO: rename to TopCoordType
   */
  typedef enum {
    TPP_TTYPE_BON = 1,     //!< bond with square potential
    TPP_TTYPE_ANG = 2,     //!< angle with square potential
    TPP_TTYPE_RBDIH = 3,   //!< dihedral with Ryckaert-Belleman potential
    TPP_TTYPE_IMPDIH = 4,  //!< dihedral with square potential
    TPP_TTYPE_SYMDIH = 5,  //!< dihedral with proper potential (gmx proper)
    TPP_TTYPE_PAIR = 6,    //!< turn pair on
    TPP_TTYPE_EXCL = 7,    //!< turn exclusion on (not used yet)
    TPP_TTYPE_SPECIMP = 8  //!< special improper type (TODO: do we need it?)
  } TopologyType;

  /** \brief Potential parameter set for similar internal coordinates.
    *
    *  This structure defines parameter set for some interaction types. It
    *  corresponds to a separate definition in .prm file (CHARMM) or in ffbonded.itp file (GMX).
    *
    * TODO: rename in some way
    */
  struct TopCoord {
    long int dbid;         //!< ID in corresponding SQL table
    std::string defname;   //!< define that will be reference topology in future .itp
    TopologyType type;     //!< type of the parameter set
    short f; //!< f = -1 (means parameter lack). Other correspond to ftype in GMX notation

    /** \brief c0..c5 values parametrize the potential interaction
     * function.
     *
     * c0-5 values have different meaning for different TopCoord::type values:
     *
     * TPP_TTYPE_BON: equilibrium bond length (c0), bond spring constant (c11)
     * TPP_TTYPE_ANG: equilibrium angle (c0), angle spring constant (c1)
     * TPP_TTYPE_RBDIH: c0..c5 means coefficients in RB
     * TPP_TTYPE_IMPDIH: equilibrium angle (c0), angle spring constant (c1)
     * TPP_TTYPE_SYMDIH: equilibrium angle (c0), potential
     *          repeating coeff. (c1), potential amplitude (c2)
     * other: c0..c5 means nothing
     */
    double c0, c1, c2, c3, c4, c5;

                /*
     t_top_coord(): defname(""), type(1),f(-1),c0(0),c1(0),c2(0),c3(0),c4(0),c5(0) {;}
     t_top_coord(const t_top_coord &_t): defname(_t.defname), type(_t.type),f(_t.f),c0(_t.c0),c1(_t.c1),c2(_t.c2),c3(_t.c3),c4(_t.c4),c5(_t.c5) {;}
                 */
  };

  /** \brief Single record combining interaction function for exact atoms
    * of the molecule.
    *
    * TopElement::defname defines interaction described in TopCoord
    * i,j,k,l defines atom types. For pair interaction only i and j make
    * sense. For tri-centered interactions only i,j and k make sense.
    *
    */
  struct TopElement {
    std::string defname;
    TppIndex i, j, k, l;
  };

  /** \brief Container of interaction potential definitions
    *
    * Array can be ordered by interaction type (alphabetically), type of
    * interaction (TopCoord::type) and by GMX interaction type f.
    */
  typedef boost::multi_index::multi_index_container<
      TopCoord,
      boost::multi_index::indexed_by<
       boost::multi_index::ordered_unique<
         boost::multi_index::member<TopCoord, std::string, &TopCoord::defname>
                         >, // key by defname
       boost::multi_index::ordered_non_unique<
         boost::multi_index::member<TopCoord, TopologyType, &TopCoord::type>
                         >, // key by directive
       boost::multi_index::ordered_non_unique<
                           boost::multi_index::member<TopCoord, short, &TopCoord::f>
                         > // key by function
      >
                > TopMap;

  /** \brief Container of exact interactions
   *
   * Array can easily be ordered by interaction types (defname)
   *
   */
  typedef boost::multi_index::multi_index_container<
      TopElement,
      boost::multi_index::indexed_by<
        boost::multi_index::sequenced<>, // array key
        boost::multi_index::ordered_non_unique<
          boost::multi_index::member<TopElement, std::string, &TopElement::defname>
                          >
                        > // key by defname
           > TopArray;

  /** \brief Internal molecule topology definition.
   *
   * This class describes both molecular structure and mechanical
   * structure (topology). Now supports single residue.
   *
   */
  class Topology {
  public:
    TopMap parameters;     //!< interaction potentials definitions
    TopArray elements;     //!< map internal coordinates to definitions
    AtomArray atoms;       //!< atoms (nb params & structure)
    OpenBabel::OBMol mol;  //!< corresponding OBMol object

    unsigned short nrexcl; //!< number of exclusions (see gmx manual)
    std::string ffinclude; //!< string for mother parameters include for TOP-file
    std::string ffdefaults;//!< string for force field defaults (required for expanded mode)
    std::string ffinfo;    //!< information about the force field
    std::string ffcheck;   //!< define that should be defined in force field ITP
    std::string res_name;  //!< residue name (short & caps)
    std::string name;      //!< molecule name (long)

    /* serialization support */

    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & BOOST_SERIALIZATION_NVP(parameters);
      ar & BOOST_SERIALIZATION_NVP(elements);
      ar & BOOST_SERIALIZATION_NVP(atoms);
      ar & BOOST_SERIALIZATION_NVP(nrexcl);
      ar & BOOST_SERIALIZATION_NVP(ffinclude);
      ar & BOOST_SERIALIZATION_NVP(res_name);
    }
  };

  //TODO: DO WE NEED IT?
  typedef std::vector<double*> t_optimize_list; /// optimize constant list

}

#endif
