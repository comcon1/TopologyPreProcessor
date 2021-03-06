#ifndef TPP_GAMESSCALL_H
#define TPP_GAMESSCALL_H

#include "global.hpp"


// estimation default control
//       - GAMESS -
#define TPP_DEFAULT_NICE 19
#define TPP_GAMESS_EXEC "pcgamess"
#define TPP_GAMESS_WORKDIR "GMS_ESTIM"
#define TPP_GAMESS_PUNCH "PUNCH"
// input control
#define TPP_GAMESS_DIRSCF true

// --
typedef ublas::bounded_vector<double, TPP_MAX_FREQ_NUM> tpp_frequencies;
typedef ublas::bounded_matrix<double, TPP_MAX_FREQ_NUM, TPP_MAX_FREQ_NUM>  tpp_hessian;


template<typename Type>
class gms_abstract {
  protected:
    // estimation restores from cash
    bool cash_only; // if no data in cash exception throwes
    // execution parameters
    string     gmsexec;
    string       punch;
    string     workdir;
    char          nice;
    // input paramters
    map<string,string> params;

    void    execute() throw (tpp_exception);
    void buildinput() throw (tpp_exception);
    virtual Type    read_value() throw (tpp_exception) = 0;
    virtual init(vector<string> params);
  public:
    gms_abstract(string a, string b, string c, char n, map<string,string> &p):
      gmsexec(a), punch(b), workdir(c), nice(n), params(p) {;}
    gms_abstract(map<string,string> &p):
      gmsexec(TPP_GAMESS_EXEC), punch(TPP_GAMESS_PUNCH), 
      workdir(TPP_GAMESS_WORKDIR), nice(TPP_DEFAULT_NICE), params(p) {;}
    Type operator(tpp_atom_array &) throw (tpp_exception) = 0;
    virtual ~gms_abstract();
};

class gms_energy: public gms_abstract<double> {
  private:
    virtual void buildinput() throw (tpp_exception);
    virtual double read_value() throw (tpp_exception);
  public:
    gms_energy(): cash_only(true) {;};
    gms_energy(tpp_input_params &p):
      internals(i) {
      gms_abstract::gms_abstract(p);
      p.insert(tpp_input_param("runtyp","energy"));
    };
    virtual ~gms_energy();
};

class gms_optimize: public gms_abstract<tpp_atom_array> {
  private:
    tpp_internals_array internals;
    virtual tpp_atom_array read_value() throw (tpp_exception);
  public:
    gms_optimize(tpp_input_params &p, tpp_internals_array &i):
      internals(i) {
      gms_abstract::gms_abstract(p);
      p.insert(tpp_input_param("runtyp","optimize"));
      p.insert(tpp_input_param("nzvar",lexical_cast<string>(i.size()) ));
    };
    virtual ~gms_optimise();
};

class gms_hessian: public gms_abstract<tpp_hessian> {
  private:
    virtual tpp_hessian read_value() throw (tpp_exception);
  public:
    gms_hessian(): cash_only(true) {;};
    gms_hessian(tpp_input_params &p, tpp_internals_array &t) {
      gms_abstract::gms_abstract(p);
      p.insert(tpp_input_param("runtyp","hessian"));
    };
    virtual ~gms_hessian();
};

class gms_fconstants: public gms_abstract<tpp_top_array> {
  private:
    tpp_internals_array internals;
    virtual tpp_top_array read_value() throw (tpp_exception);
  public:
    gms_hessian(tpp_input_params &p, tpp_internals_array &i):
      internals(i) {
      gms_abstract::gms_abstract(p);
      p.insert(tpp_input_param("runtyp","hessian"));
      p.insert(tpp_input_param("nzvar",lexical_cast<string>(i.size()) ));
      p.insert(tpp_input_param("prtifc",".t."));
      p.insert(tpp_input_param("decomp",".t."));
    };
    virtual ~gms_hessian();
};

class gms_frequencies: public gms_abstract<tpp_frequencies> {
  private:
    virtual tpp_frequencies read_value() throw (tpp_exception);
  public:
    gms_frequencies(): cash_only(true) {;};
    gms_frequencies(tpp_input_params &p) {
      gms_abstract::gms_abstract(p);
      p.insert(tpp_input_param("runtyp","hessian"));
      p.insert(tpp_input_param("vibanl",".t."));
    }
    virtual ~gms_hessian();
};

#endif

   