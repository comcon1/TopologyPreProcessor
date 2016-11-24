#ifndef TPP_BONDDEFINER_H
#define PP_BONDDEFINER_H

#include "db_base.hpp"

namespace tpp
{
/**
 * \brief A class that defines bonds between atoms, I guess?
 *
 */
class BondDefiner: public DbBase {
public:

	BondDefiner(Parameters, Topology &);
	virtual ~BondDefiner();
	void bond_align();
	void log_needed_bonds();

private:
	InternalsArray bonds;
	Topology &tp;
	std::map<std::string, std::string> namemap; //! map of uname -> name
	short ffid;
	bool genpairs;
	std::ofstream qalcfile;

	// methods
	bool connect_db() override;
	void fill_bonds();
	void fill_angles();
	void fill_dihedrals();
	void fill_special();
	void fill_impropers();
	void fill_pairs();

};
} // of namespace tpp

#endif
