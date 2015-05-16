#ifndef _MAXIMUM_H_
#define _MAXIMUM_H_

#include <limits.h>

class Maximum{
public:
Maximum(): maxValue(INT_MIN) {}
  Maximum& add(int i){
    // maxValue  =  std::max(maxValue, i);
    // the above won't work for SolarisStudio 12.3
    maxValue = maxValue > i ? maxValue : i;
    return *this;
  }
  int max() {
    return this->maxValue;
  }
  void reset() {
    this->maxValue = INT_MIN;
  }
private:
  // int counter; // how many time counted
  int maxValue; // maximum value so far
};

#endif /* _MAXIMUM_H_ */
