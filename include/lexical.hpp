/** \file lexical.hpp
 *
 * \brief This file desperately needs clarification.
 *
 * Also, this file was originally in cp1251! MS Windows user detected!
 *
 * LEX namespace includes cast and parse classes for parametric files read.
 * Also it has some useful overrides of classes of boost::spirit::actor.
 */
#ifndef LEXICAL_SPIRIT_H
#define LEXICAL_SPIRIT_H

// for some spirit features??
#ifdef WIN32
#undef WIN32
#endif

#include "global.hpp"
#include <cassert>
#include "logger.hpp"

#include "core.hpp"

#include <boost/type_traits.hpp>


#if HAVE_BOOST_SPIRIT_CORE_HPP

	#include <boost/spirit/core.hpp>
	#include <boost/spirit/iterator/file_iterator.hpp>
	#include <boost/spirit/dynamic/stored_rule.hpp>
	#include <boost/spirit/actor/ref_value_actor.hpp>
	#include <boost/spirit/actor/assign_actor.hpp>
	#include <boost/spirit/actor/increment_actor.hpp>
	#include <boost/spirit/actor/push_back_actor.hpp>
	#include <boost/spirit/actor/insert_key_actor.hpp>
	#include <boost/spirit/actor/clear_actor.hpp>
	#include <boost/spirit/error_handling/exceptions.hpp>
	#include <boost/spirit/utility.hpp>
	#include <boost/spirit/utility/chset.hpp>
	#include <boost/spirit/utility/chset_operators.hpp>

#elif HAVE_BOOST_SPIRIT_HOME_CLASSIC_CORE_HPP

	#include <boost/spirit/home/classic/core.hpp>
	#include <boost/spirit/home/classic/iterator/file_iterator.hpp>
	#include <boost/spirit/home/classic/dynamic/stored_rule.hpp>
	#include <boost/spirit/home/classic/actor/ref_value_actor.hpp>
	#include <boost/spirit/home/classic/actor/assign_actor.hpp>
	#include <boost/spirit/home/classic/actor/increment_actor.hpp>
	#include <boost/spirit/home/classic/actor/push_back_actor.hpp>
	#include <boost/spirit/home/classic/actor/insert_key_actor.hpp>
	#include <boost/spirit/home/classic/actor/clear_actor.hpp>
	#include <boost/spirit/home/classic/error_handling/exceptions.hpp>
	#include <boost/spirit/home/classic/utility.hpp>
	#include <boost/spirit/home/classic/utility/chset.hpp>
	#include <boost/spirit/home/classic/utility/chset_operators.hpp>
#endif

#include <set>
#include <boost/lexical_cast.hpp>

using boost::lexical_cast;
using boost::bad_lexical_cast;
using namespace boost::spirit::classic;

namespace lex {
    using tpp::AtomArray;
    using tpp::Atom;
    struct Tcoord {
        double x, y, z;
    };

    typedef int flag;
    typedef file_iterator<char> iterator;
    typedef scanner<iterator> scanner;
    typedef stored_rule<scanner> rule_t;
    typedef rule<scanner> rule_r;
    typedef enum {
        endblock = 0, endfile = 1
    } manipulator;

    static std::ostringstream outs;

    const int eolbefore_flag = 1; // 000001
    const int multistring_flag = 2; // 000010

    static int offset;
    static unsigned incrementor;
    static unsigned nullVal = 0;
    static bool scan_str_flag;

    std::set<void *> BLOCKED_VARIABLES;

    // namespace for debug output
    namespace echo {
        static void val_start(iterator f, iterator l) {
            std::string v;
            std::copy(f, l, std::inserter(v, v.end()));
            for (int i = 0; i < offset; ++i)
                outs.put(' ');
            outs << std::setw(20) << std::left << std::string("[ ") + v
                    << (scan_str_flag ? "" : "= ");
        }

        static void block_start(iterator f, iterator l) {
            std::string v;
            std::copy(f + 1, l - 1, std::inserter(v, v.end()));
            for (int i = 0; i < offset; ++i)
                outs.put('-');
            outs << "|" << std::setw(20) << std::left << v << "\n";
            offset += 4;
        }

        static void block_stop(iterator f, iterator l) {
            offset -= 4;
            for (int i = 0; i < offset; ++i)
                outs.put('-');
            outs << "block end." << std::endl;
        }

        static void val_print(iterator f, iterator l) {
            std::string v;
            std::copy(f, l, std::inserter(v, v.end()));
            outs << ((incrementor) ? ", " : "") << v << std::endl;
        }

        static void print(iterator f, iterator l) {
            std::string v;
            std::copy(f, l, std::inserter(v, v.end()));
            std::cout << v << std::endl;
        }

        static void val_stop(iterator, iterator) {
            outs << " ]" << std::endl;
        }
    } // end echo namespace

    // Класс для актора boost::spirit
    struct assign_sp {
    public:
        template<typename T, typename ValueT>
        void act(T& ref, ValueT const& value) const {
            ref = value;
        }

        template<typename T, typename IteratorT>
        void act(std::set<T>& ref, IteratorT const& first,
                IteratorT const& last) const {
            std::string s;
            std::copy(first, last, std::inserter(s, s.end()));
            try {
                ref.insert(lexical_cast<T>(s));
            } catch (boost::bad_lexical_cast &e) {
                outs << "*\n";
                throw_(first, 1);
            }
        }

        template<typename T, typename IteratorT>
        void act(std::vector<T>& ref, IteratorT const& first,
                IteratorT const& last) const {
            std::string s;
            std::copy(first, last, std::inserter(s, s.end()));
            try {
                ref.push_back(lexical_cast<T>(s));
            } catch (boost::bad_lexical_cast &e) {
                outs << "*\n";
                throw_(first, 1);
            }
        }

        template<typename T, typename IteratorT>
        void act(T& ref, IteratorT const& first, IteratorT const& last) const {
            void *ptr = &ref;
            std::string s;
            std::copy(first, last, std::inserter(s, s.end()));
            try {
                ref = lexical_cast<T>(s);
            } catch (boost::bad_lexical_cast &e) {
                outs << "*\n";
                throw_(first, 1);
            }
            if (!(BLOCKED_VARIABLES.insert(ptr)).second && !incrementor) {
                outs << "*\n";
                throw_(first, 2); // повторное включение | неверный парс
            }
        }
    };

    // *********************************************************************
    // **
    // **********************************************************************
    struct insert_mi_action {
        template<typename T, typename ValueT, typename ReferentT>
        void act(T& ref_, ValueT const& value_, ReferentT const& key_) const {
            ref_.insert(value_);
        }

        template<typename T, typename ValueT, typename IteratorT>
        void act(T& ref_, ValueT const& value_, IteratorT const& first_,
                IteratorT const& last_) const {
            ref_.insert(value_);
        }
    };

    template<typename T, typename ValueT>
    inline ref_const_ref_value_actor<T, ValueT, insert_mi_action> insert_mi_a(
            T& ref_, ValueT const& value_) {
        return ref_const_ref_value_actor<T, ValueT, insert_mi_action>(ref_, value_);
    }
    // ***********************************************************************
    // **
    //  T - AtomArray
    //  ValueT - Atom
    struct update_mi_action {
        template<typename ReferentT>
        void act(AtomArray& ref_, Atom const& value_, ReferentT const& key_) const {
            AtomArray::iterator it = ref_.find(value_.index);
            assert(it != ref_.end());
            Atom atom0 = *it;
            atom0.atom_name = value_.atom_name;
            atom0.charge = value_.charge;
            atom0.atom_type = value_.atom_type;
            atom0.comment = value_.comment;
            atom0.c_gnr = value_.c_gnr;
            ref_.replace(it, atom0);
        }

        template<typename IteratorT>
        void act(AtomArray& ref_, Atom const& value_, IteratorT const& first_,
                IteratorT const& last_) const {
            AtomArray::iterator it = ref_.find(value_.index);
            assert(it != ref_.end());
            Atom atom0 = *it;
            atom0.atom_name = value_.atom_name;
            atom0.charge = value_.charge;
            atom0.atom_type = value_.atom_type;
            atom0.comment = value_.comment;
            atom0.c_gnr = value_.c_gnr;
            ref_.replace(it, atom0);
        }
    };

    inline ref_const_ref_value_actor<tpp::AtomArray, Atom, update_mi_action> update_mi_a(
            AtomArray& ref_, Atom const& value_) {
        return ref_const_ref_value_actor<tpp::AtomArray, Atom, update_mi_action>(
                ref_, value_);
    }

    // ***********************************************************************
    // ***
    struct assign_with_cast_action {
        template<typename T, typename ValueT>
        void act(T& ref_, ValueT const& value_) const {
            ref_ = static_cast<T>(value_);
        }
        template<typename T, typename IteratorT>
        void act(T& ref_, IteratorT const& first_, IteratorT const& last_) const {
            typedef T value_type;
    #ifndef BOOST_NO_TEMPLATED_ITERATOR_CONSTRUCTORS
            value_type value(first_, last_);
    #else
            value_type value;
            std::copy(first_, last_, std::inserter(value, value.end()));
    #endif
            ref_ = static_cast<T>(value);
        }
    };

    template<typename T>
    inline ref_value_actor<T, assign_with_cast_action> assign_cast_a(T& ref_) {
        return ref_value_actor<T, assign_with_cast_action>(ref_);
    }

    template<typename T, typename ValueT>
    inline ref_const_ref_actor<T, ValueT, assign_with_cast_action> assign_cast_a(
            T& ref_, ValueT const& value_) {
        return ref_const_ref_actor<T, ValueT, assign_with_cast_action>(ref_, value_);
    }

    // Структура для определения блока параметров (или блоков :)
    struct single_block {
        rule_t *rules;
        unsigned maxnumrules;
        unsigned numrules;
    public:
        rule_t value;
        std::string header, footer;
        single_block(const char *s, int num) :
                maxnumrules(num), numrules(0) {
            header = std::string("<") + std::string(s) + std::string(">");
            footer = std::string("<") + std::string(s) + ".end>";
            rules = new rule_t[num];
        }
        // добавление нового правила к цепочке правил
        friend single_block& operator <<(single_block& c, rule_t nxt_rule) {
            c.rules[c.numrules] = nxt_rule;
            c.numrules++;
            assert(c.numrules <= c.maxnumrules);
            return c;
        }
        // добавление концевого манипулятора
        friend single_block& operator <<(single_block& c, manipulator man) {
            if (man == endblock) {
                c.value = *eol_p >> *blank_p
                        >> (str_p(c.header.c_str()))[&echo::block_start] >> *blank_p
                        >> eol_p;
                rule_t sin = c.rules[0];
                for (int i = 1; i < c.numrules; ++i)
                    sin = sin.copy() | c.rules[i];
                c.value = c.value.copy() >> repeat_p(1, c.maxnumrules)[sin.copy()]
                        >> *blank_p >> (str_p(c.footer.c_str()))[&echo::block_stop]
                        >> *blank_p >> *eol_p;
            }
            return c;
        }
        // удаление цепочки правил
        ~single_block() {
            delete[] rules;
        }
    };

    // Переделка стандартного актора boost::spirit::assign_a
    template<typename T>
    inline ref_value_actor<T, assign_sp> assign_at(T &ref) {
        return ref_value_actor<T, assign_sp>(ref);
    }

    // добавления правил для параметра типа vector
    template<typename VectorT>
    static rule_t custom_param(const char* nam, std::vector<VectorT> &vec, flag fl,
            unsigned num) {
        rule_t result;
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result =
                result.copy() >> ch_p('[')[assign_a(incrementor, nullVal)]
                        >> *blank_p
                        >> (((*(graph_p - ',' - ']'))[assign_at(vec)]))[&echo::val_print]
                        >> *blank_p
                        >> !repeat_p(1, num)[(ch_p(','))[increment_a(incrementor)]
                                >> *space_p
                                >> ((*(graph_p - ',' - ']'))[assign_at(vec)])[&echo::val_print]
                                >> *blank_p] >> ch_p(']') >> *blank_p
                        >> ((eol_p)[&echo::val_stop])[assign_a(incrementor, nullVal)];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // добавления правил для параметра типа pair
    template<typename PairT1, typename PairT2>
    static rule_t custom_param(const char* nam, std::pair<PairT1, PairT2> &vec,
            flag fl, unsigned num) {
        rule_t result;
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result =
                result.copy() >> ch_p('[')[assign_a(incrementor, nullVal)]
                        >> *blank_p
                        >> (((*(graph_p - ',' - ']'))[assign_at(vec.first)]))[&echo::val_print]
                        >> *blank_p
                        >> ((ch_p(','))[increment_a(incrementor)] >> *blank_p
                                >> ((*(graph_p - ',' - ']'))[assign_at(vec.second)])[&echo::val_print]
                                >> *blank_p) >> ch_p(']') >> *blank_p
                        >> ((eol_p)[&echo::val_stop])[assign_a(incrementor, nullVal)];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // for Tcoord type special
    static rule_t custom_param(const char* nam, Tcoord &vec, flag fl,
            unsigned num) {
        rule_t result;
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result = result.copy() >> ch_p('[')[assign_a(incrementor, nullVal)]
                >> *blank_p
                >> (((*(graph_p - ',' - ']'))[assign_at(vec.x)]))[&echo::val_print]
                >> *blank_p >> (ch_p(','))[increment_a(incrementor)] >> *blank_p
                >> ((*(graph_p - ',' - ']'))[assign_at(vec.y)])[&echo::val_print]
                >> *blank_p >> (ch_p(','))[increment_a(incrementor)] >> *blank_p
                >> ((*(graph_p - ',' - ']'))[assign_at(vec.z)])[&echo::val_print]
                >> *blank_p >> ch_p(']') >> *blank_p
                >> ((eol_p)[&echo::val_stop])[assign_a(incrementor, nullVal)];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // правила для коротких массивов типа (Type *)
    template<typename PointerT>
    static rule_t custom_param(const char* nam, PointerT *vec, flag fl,
            unsigned num) {
        rule_t result;
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result =
                result.copy() >> ch_p('[')[assign_a(incrementor, nullVal)]
                        >> *blank_p
                        >> (((*(graph_p - ',' - ']'))[assign_at(*vec)]))[&echo::val_print]
                        >> *blank_p
                        >> !repeat_p(num - 1)[(ch_p(','))[increment_a(incrementor)]
                                >> *blank_p
                                >> ((*(graph_p - ',' - ']'))[assign_at(
                                        vec[incrementor])])[&echo::val_print]
                                >> *blank_p] >> ch_p(']') >> *blank_p
                        >> ((eol_p)[&echo::val_stop])[assign_a(incrementor, nullVal)];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // Добавление цепочек правил для параметров типа std::set
    template<typename SetT>
    static rule_t custom_param(const char* nam, std::set<SetT> &vec, flag fl,
            unsigned num) {
        rule_t result;
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result =
                result.copy() >> ch_p('[')[assign_a(incrementor, nullVal)]
                        >> *blank_p
                        >> (((*(graph_p - ',' - ']'))[assign_at(vec)]))[&echo::val_print]
                        >> *blank_p
                        >> !repeat_p(1, num)[(ch_p(','))[increment_a(incrementor)]
                                >> *space_p
                                >> ((*(graph_p - ',' - ']'))[assign_at(vec)])[&echo::val_print]
                                >> *blank_p]
                        >> (ch_p(']'))[assign_a(incrementor, nullVal)]
                        >> (*blank_p >> eol_p)[&echo::val_stop];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // добавление правил для остальных типов (enum, bool, char, int, double, ...)
    // (... функция частично специализирована выше ...)
    template<typename T>
    static rule_t custom_param(const char *nam, T &val, flag fl, unsigned num) {
        rule_t result, empty;
        assert(!::boost::is_pointer<T>::value);
        result = *blank_p >> (str_p(nam))[&echo::val_start] >> *blank_p >> ch_p('=')
                >> *blank_p;
        result = result.copy() >> ((*graph_p)[assign_at(val)])[&echo::val_print];
        result = result.copy() >> *blank_p >> (eol_p)[&echo::val_stop];
        if (fl & eolbefore_flag)
            result = *eol_p >> result.copy();
        return result;
    }

    // добавление правила для блока без параметров.
    static rule_t scan_strings(std::vector<std::string> &vec) {
        rule_t result = (*(anychar_p - eol_p - ch_p('<')))[assign_a(scan_str_flag,
                1)][push_back_a(vec)][&echo::val_start][&echo::val_stop][assign_a(
                scan_str_flag, 0)] >> (eol_p);
        return result;
    }

}
#endif

