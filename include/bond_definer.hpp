#ifndef TPP_BONDDEFINER_H
#define PP_BONDDEFINER_H

#include "db_base.hpp"

namespace tpp
{
/**
 * \brief A class that defines bonds between atoms, I guess?
 *
 */
class bond_definer: public db_base {
  private:
    InternalsArray bonds;
    Topology &tp;
    std::map<std::string, std::string> namemap; // map of uname -> name ))
    short ffid;
    bool  genpairs;
    std::ofstream qalcfile;

    // methods
    virtual bool connect_db();
    void fill_bonds();
    void fill_angles();
    void fill_dihedrals();
    void fill_special();
    void fill_impropers();
    void fill_pairs();
  public:

    bond_definer(Parameters, Topology &);
    virtual ~bond_definer();
    void bond_align();
    void log_needed_bonds();
};
}

#endif
