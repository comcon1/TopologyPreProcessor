/** \file core.hpp
 *
 *	\brief This file contains core definitions for TPP project.
 *
 *	Those include atoms, topologies and their respective containers.
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

//
//	Those are quite mysterious
//
#define TPP_GMX_EXTERNAL
#define TPP_MAX_ATOM_NUM 200
#define TPP_INDEX unsigned
#define TPP_MAX_FREQ_NUM 600
#define TPP_MAX_BONDS 4

namespace tpp {
	// atom and atom_array definition
	typedef boost::numeric::ublas::bounded_vector<double, 3> Point;

	/**
	 *	\brief This class describes atoms!
	 *
	 */
	struct Atom {
		TPP_INDEX index; //!< place in vector array
		TPP_INDEX oldindex; //!< This surely means something!
		std::string atom_type; //!< This surely means something!
		std::string atom_type2; //!< This surely means something!
		std::string atom_name; //!< This surely means something!
		std::string old_aname; //!< This surely means something!
		unsigned char ncharge; //!< nuclear charge
		std::string res_name; //!< This surely means something!
		std::string comment; //!< This surely means something!
		unsigned char mol_id; //!< This surely means something!
		Point coord; //!< This surely means something!
		boost::array<TPP_INDEX, TPP_MAX_BONDS> connects; //!< This surely means something!
		unsigned char num_connects; //!< This surely means something!
		double charge; //!< This surely means something!
		double mass;TPP_INDEX c_gnr; //!< This surely means something!
		std::string qmname; //!< This surely means something!
	};

	typedef boost::multi_index::multi_index_container<Atom,
			boost::multi_index::indexed_by<boost::multi_index::ordered_unique<boost::multi_index::member<Atom, TPP_INDEX, &Atom::index> >, // key index (like array)
			boost::multi_index::ordered_non_unique<boost::multi_index::member<Atom, std::string, &Atom::atom_name> >, // key by name
			boost::multi_index::ordered_non_unique<boost::multi_index::member<Atom, TPP_INDEX, &Atom::c_gnr> > // key by charge group
			> > AtomArray;

	/// This enumeration surely relates to input somehow
	typedef enum {
		TPP_IF_PDB = 0,
		TPP_IF_GRO = 1,
		TPP_IF_G96 = 2,
		TPP_IF_GAMOPT = 3,
		TPP_IF_GAMHESS = 4,
		TPP_IF_GAMSP = 5
	} InputFormat;

	/// This enumeration surely relates to output somehow
	typedef enum {
		TPP_OF_PDB = 0, TPP_OF_GRO = 2, TPP_OF_G96 = 3, TPP_OF_GAMIN = 4
	} OutputFormat;

	/// internal coordinates
	typedef enum {
		TPP_IC_BON = 1, TPP_IC_ANG = 2, TPP_IC_DIH = 3, TPP_IC_CHIR = 4
	} IntCoordType;

	/// M... some other coordinates?..
	typedef struct {
		IntCoordType type;
		TPP_INDEX i, j, k, l;
		std::string defname;
	} IntCoord;

	typedef std::vector<IntCoord> InternalsArray;

	/// topology parameters definition
	typedef enum {
		TPP_TTYPE_BON = 1,
		TPP_TTYPE_ANG = 2,
		TPP_TTYPE_RBDIH = 3,
		TPP_TTYPE_IMPDIH = 4,
		TPP_TTYPE_SYMDIH = 5,
		TPP_TTYPE_PAIR = 6,
		TPP_TTYPE_EXCL = 7,
		TPP_TTYPE_SPECIMP = 8
	} TopologyType;

	/// Some coordinates, it seems like
	struct TopCoord {
		long int dbid;
		std::string defname;
		TopologyType type;
		short f; //!< f = -1 (means parameter lack)
		double c0, c1, c2, c3, c4, c5;/*
		 t_top_coord(): defname(""), type(1),f(-1),c0(0),c1(0),c2(0),c3(0),c4(0),c5(0) {;}
		 t_top_coord(const t_top_coord &_t): defname(_t.defname), type(_t.type),f(_t.f),c0(_t.c0),c1(_t.c1),c2(_t.c2),c3(_t.c3),c4(_t.c4),c5(_t.c5) {;}  */
	};

	struct TopElement {
		std::string defname;TPP_INDEX i, j, k, l;
	};

	/// DEFINITIONS CONTAINER
	typedef boost::multi_index::multi_index_container<
			TopCoord,
			boost::multi_index::indexed_by<
			boost::multi_index::ordered_unique<
			boost::multi_index::member<TopCoord, std::string, &TopCoord::defname> >, // key by defname
			boost::multi_index::ordered_non_unique<
			boost::multi_index::member<TopCoord, TopologyType, &TopCoord::type> >, // key by directive
			boost::multi_index::ordered_non_unique<boost::multi_index::member<TopCoord, short, &TopCoord::f> > // key by function
			> > TopMap;

	typedef boost::multi_index::multi_index_container< /* TOPOLOGY CONTAINER */
			TopElement,
			boost::multi_index::indexed_by<
			boost::multi_index::sequenced<>, // array key
			boost::multi_index::ordered_non_unique<
			boost::multi_index::member<TopElement, std::string, &TopElement::defname> > > // key by defname
	> TopArray;

	// internal topology definition
	class Topology {
	public:
		TopMap parameters;
		TopArray elements;
		AtomArray atoms;
		OpenBabel::OBMol mol;

		unsigned short nrexcl;
		std::string ffinclude;
		std::string ffinfo;
		std::string res_name;
		std::string name;

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


	typedef std::vector<double*> t_optimize_list; /// optimize constant list

}

#endif
