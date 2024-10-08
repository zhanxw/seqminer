#ifndef _ORDEREDMAP_H_
#define _ORDEREDMAP_H_

#include <cstdio>
#include <cassert>
#include <vector>
#include <map>
#include "Exception.h"
#include "R.h"

/**
 * a simple OrderedMap class
 * use operator[] to insert elements
 * use at() to access elements
 * NOTE: KEY will be stored twice.
 */
template <class KEY, class TYPE>
class OrderedMap{
 public:
  bool find(const KEY& key) const {
    if (this->keyTypeMap.find(key) == this->keyTypeMap.end()){
      return false;
    }
    return true;
  }
  TYPE& operator[] (const KEY& key) {
    if (!this->find(key)){
      this->keyVec.push_back(key);
    }
    return this->keyTypeMap[key];
  }
  const TYPE& operator[] (const KEY& key) const{
    if (!this->find(key)){
      throw "key not found in OrderedMap";
    }
    return this->keyTypeMap.find(key)->second;
  }
  void front(KEY* k, TYPE* v) {
    *k = this->keyVec.front();
    *v = this->keyTypeMap[(*k)];
  }
  bool at(unsigned int idx, KEY* k, TYPE* v) {
    if (idx >= this->size()) return false;
    *k = this->keyVec[idx];
    *v = this->keyTypeMap[(*k)];
    return true;
  }
  bool at(unsigned int idx, KEY* k, TYPE* v) const {
    if (idx >= this->size()) return false;
    *k = this->keyVec[idx];
    if (this->keyTypeMap.find(*k) == this->keyTypeMap.end()){
      v =NULL;
    } else {
      *v = this->keyTypeMap.find(*k)->second;
    }
    return true;
  }
  const KEY& keyAt(unsigned int idx) const {
    if (idx >= this->size()) {
      log_error("Index out of bound, now quitting...");
    }
    return this->keyVec[idx];
  }
  const TYPE& valueAt(unsigned int idx) const {
    if (idx >= this->size()) {
      log_error("Index out of bound, now quitting...");
    }
    const KEY& k = this->keyVec[idx];
    if (this->keyTypeMap.find(k) == this->keyTypeMap.end()){
      REprintf("Cannot find KEY in valueAt()\n");
      REprintf("Critical error happening!\n");
      //abort();
      // just return the first
      return this->keyTypeMap.begin()->second;
    } else {
      return this->keyTypeMap.find(k)->second;
    }
  }
  bool value(const KEY& key, TYPE* val) const{
    typename std::map < KEY, TYPE >::const_iterator iter;
    iter = this->keyTypeMap.find(key);
    if (iter != this->keyTypeMap.end()){
      *val = iter->second;
      return true;
    } else {
      return false;
    }
  }    /**
        * compare
        */
  void compareKey(const OrderedMap<KEY, TYPE>& other, int* overlap, int* thisUniqueKeys, int* otherUniqueKeys) const{
    (*overlap) = (*thisUniqueKeys) = 0;
    //assert(overlap && thisUniqueKeys && otherUniqueKeys);
    for (unsigned int i = 0; i != this->size(); i ++ ){
      KEY& k = this->keyAt(i);
      if (other.find(k))
        (*overlap)++;
      else
        (*thisUniqueKeys)++;
    }
    *otherUniqueKeys = other.size() - *overlap;
  }
  unsigned int size() const { return this->keyVec.size();} ;
  void clear() {
    this->keyVec.clear();
    this->keyTypeMap.clear();
  };
  const std::vector<KEY>& getKey() const{
    return this->keyVec;
  };
 private:
  std::vector < KEY > keyVec;
  std::map < KEY, TYPE > keyTypeMap;
};

/**
 * For each element in @param input, give its order
 * e.g. input = {a, b, c, d, b}
 *      output = {a:0, b:1, c:2, d:3}
 */
template<class KEY>
inline int numberVectorAsMap(const std::vector<KEY>& input,
                             OrderedMap<KEY, int>* output) {
  output->clear();
  int n = input.size();
  for(int i = 0; i < n; ++i){
    if (output->find(input[i])){
      continue;
    }
    const int s = output->size();
    output[input[i]] = s;
  }
  return 0;
}

#endif /* _ORDEREDMAP_H_ */
