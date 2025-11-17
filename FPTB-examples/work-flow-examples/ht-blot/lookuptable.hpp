#ifndef LOOKUPTABLE
#define LOOKUPTABLE


#include <set>
#include <tuple>
#include <map>
#include <iostream>

//
// Implementation of lookuptables
// J = type of index column
// T... = type of other columns, each must be unique
// Unique column label types can be defined by luttype<A> 
// Each row must also have a unique value of the index
//
// Basic usage:
//
// class G_lut{};
// typedef luttype<G_lut,int> G;
// //same for H1, H2, H3 ...
// 
// lut.insert<G>(J(j), G(g)); //inserts g into row j, column G
// lut.lookup_idx<H1,H2,H3>(H1(h1),H2(h2),H3(h3)) //looks up the index j that appear in the same row as h1, h2, h3
// J const j = *lut.lookup_idx<H1,H2,H3>(H1(h1),H2(h2),H3(h3)).begin() //same as above, but when only one result is expected
// lut.lookup_col<G>(J(j)); //looks  up column G from index j
// lut.lookup_col<G,H1,H2,H3>(H1(h1), H2(h2), H3(h3)) //looks up the G that appear in the same row as h1, h2, h3
//
// Some template arguments can also be omitted:
// lut.insert(J(j), G(g)); //inserts g into row j, column G
// lut.lookup_idx(H1(h1),H2(h2),H3(h3)) //looks up the index j that appear in the same row as h1, h2, h3
// J const j = *lut.lookup_idx(H1(h1),H2(h2),H3(h3)).begin() //same as above, but when only one result is expected
// lut.lookup_col<G>(J(j)); //looks  up column G from index j
// lut.lookup_col<G>(H1(h1), H2(h2), H3(h3)) //looks up the G that appear in the same row as h1, h2, h3
//
// lookuptable is zero-initialized

template<class J, class ...T>
struct lookuptable{
  std::map<J, std::tuple<T...>> row_lut;
  std::tuple< std::multimap<T,J>... > idx_lut;

  lookuptable();

  template<class G>
  void insert(J j, G g);

  template<class G, class ...S>
  void insert(J j, G g, S... s);

  template<class H>
  std::set<J> lookup_idx (H h) const;

  template<class H, class ...S>
  std::set<J> lookup_idx (H h, S... s) const;

  template<class G>
  G lookup_col (J j) const;

  template<class G, class H>
  std::set<G> lookup_col (H h) const;

  template<class G, class H, class ...S>
  std::set<G> lookup_col (H h, S... s) const;

};


template<class J, class ...T>
lookuptable<J,T...>::lookuptable(){

  std::map<J, std::tuple<T...>> row_lut;
  std::tuple< std::multimap<T,J>... > idx_lut;
}

template<class J, class ...T>
template<class G>
void lookuptable<J,T...>::insert(J j, G g){
  if (row_lut.count(j)==0){
    std::tuple<T...> tup;
    std::get<G>(tup) = g;
    row_lut.insert({j,tup});
  }
  else{
    std::get<G>(row_lut.at(j)) = g;
  }

  auto &mm = std::get<std::multimap<G,J>>(idx_lut);
  mm.insert(std::pair<G,J>(g,j));

}

template<class J, class ...T>
template<class G, class ...S>
void lookuptable<J,T...>::insert(J j, G g, S... s){
  lookuptable<J,T...>::insert(j, g);
  lookuptable<J,T...>::insert(j, s...);
}

template<class J, class ...T>
template<class H>
std::set<J> lookuptable<J,T...>::lookup_idx(H h) const {

  std::set<J> ans;
  auto j_iters = std::get<std::multimap<H,J>>(idx_lut).equal_range(h);
  for (auto j_iter=j_iters.first; j_iter!=j_iters.second; ++j_iter){
    auto j = j_iter->second;
    auto tup = row_lut.at(j);
    ans.insert(j);
  }
  return ans;
}



template<class J, class ...T>
template<class H, class ...S>
std::set<J> lookuptable<J,T...>::lookup_idx(H h, S... s) const {

  std::set<J> j1 = lookuptable<J,T...>::lookup_idx<H>(h);
  std::set<J> j2 = lookuptable<J,T...>::lookup_idx<S...>(s...);

  std::set<J> j3;
  std::set_intersection(j1.begin(), j1.end(), j2.begin(), j2.end(), std::inserter(j3, j3.begin()));

  return j3;
}

template<class J, class ...T>
template<class G>
G lookuptable<J,T...>::lookup_col(J j) const {

  auto tup = row_lut.at(j);
  return std::get<G>(tup);
}
  

template<class J, class ...T>
template<class G, class H>
std::set<G> lookuptable<J,T...>::lookup_col(H h) const {

  auto js = lookuptable<J,T...>::lookup_idx<H>(h);
  std::set<G> ans;
  for (auto j : js){
    auto g = lookuptable<J,T...>::lookup_col<G>(j);
    ans.insert(g);
  }

  return ans;
}

template<class J, class ...T>
template<class G, class H, class ...S>
std::set<G> lookuptable<J,T...>::lookup_col(H h, S... s) const {

  auto js = lookuptable<J,T...>::lookup_idx<H,S...>(h, s...);
  std::set<G> ans;
  for (auto j : js){
    auto g = lookuptable<J,T...>::lookup_col<G>(j);
    ans.insert(g);
  }

  return ans;
}

//
// Label templates
//

template<class T, class A>
struct luttype{
  T value;
  luttype(T v) : value(v) {};
  luttype() {};
  //luttype( luttype<T,A> & x){ value = x.value;};
  luttype(const luttype<T,A> & x){ value = x.value;};
};

template<class T, class A>
bool operator< (luttype<T,A> a, luttype<T,A> b){
  return a.value < b.value;
}
// ==
template<class T, class A, class B>
bool operator== (luttype<T,A> a, luttype<T,B> b){
  return a.value == b.value;
}

template<class T, class A>
bool operator== (luttype<T,A> a, T b){
  return a.value == b;
}

template<class T, class A>
bool operator== (T a, luttype<T,A> b){
  return a == b.value;
}

template<class T1, class T2, class A>
bool operator== (luttype<T1,A> a, T2 b){
  return a.value == T1(b);
}

template<class T1, class T2, class A>
bool operator== (T1 a, luttype<T2,A> b){
  return T2(a) == b.value;
}
// !=
template<class T, class A, class B>
bool operator!= (luttype<T,A> a, luttype<T,B> b){
  return a.value != b.value;
}

template<class T, class A>
bool operator!= (luttype<T,A> a, T b){
  return a.value != b;
}

template<class T, class A>
bool operator!= (T a, luttype<T,A> b){
  return a != b.value;
}

template<class T1, class T2, class A>
bool operator!= (luttype<T1,A> a, T2 b){
  return a.value != T1(b);
}

template<class T1, class T2, class A>
bool operator!= (T1 a, luttype<T2,A> b){
  return T2(a) != b.value;
}


//template<class A>
//struct lutstring{
//  std::string value;
//  lutstring(std::string v) : value(v) {};
//  lutstring() {};
//};

//template<class A>
//bool operator< (lutstring<A> a, lutstring<A> b){
//  return a.value < b.value;
//}



#endif
