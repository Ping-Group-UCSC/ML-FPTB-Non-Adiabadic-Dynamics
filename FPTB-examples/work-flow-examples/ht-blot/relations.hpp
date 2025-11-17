#ifndef RELATIONS
#define RELATIONS


#include <boost/hana.hpp>
#include <variant>
#include <set>
#include <map>
#include <iostream>
#include <typeinfo>

using namespace std;

// solve_formula( r( a, _1 , b) && s(c, _2))
// should return a list of (_1,_2) that satisfy the formula. 
// should be lazy evaluated.
// should be templated


// each relation should implement a get function
// for instance r.get(a,_1,b) should return a list of _1 that satisfy this relation.
// The implementation depends on the relation, could use a hash table or a function. 
// How to deal with infinite relations?
//
// Some relations can have a r.set(..) function which manually extends the relation
//
// Should run at run-time since the relations may only be built at run-time, although there can be some compile-time operations to make coding easier.
// indexing of hana tuples must be done at compile time.

//implement binary symmetric relations first
//
//is there some kind of multi hash table for the implementation of m-ary relations?
//one option is to use sqlite
//
//auto return type of virtual member functions not allowed in c++
//templated virtual member functions not allowed in c++ (but virtual members allowed in class templates)

class Variable {};

//vector<hana::tuple<Xs...>> get (Ys...) // where Xs, Ys are template parameter packs for the search and constraint types.
//but how to treat the case where there are no Xs ?

template<class T1, class T2> 
class binary_relation{
    //set<T1> universe1; 
    //set<T2> universe2; 
    map<T1,set<T2>> hash1;
    map<T2,set<T1>> hash2;
    
  public:
    void setr(T1 x1, T2 x2){ 
      hash1[x1].insert(x2);
      hash2[x2].insert(x1);
    };

    template<class TT1, class TT2>
    auto operator() (TT1 x1, TT2 x2);
}

//we want:
//vector<T1> r( _ , x2)
//vector<T2> r( x1 , _)
//vector< pair<T1,T2>> r( _ , _)
//bool r(x1, x2)
//
//This can be done with overloading for simple cases, but maybe we should use hana for more complex ones

template<class T1, class T2> 
template<class TT1, class TT2>
auto binary_relation::operator() (TT1

#endif
