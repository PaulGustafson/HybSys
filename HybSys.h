/*
 * Copyright (C) 2020 Paul Gustafson 
 * License: MIT  (see the file LICENSE)
 */

#ifndef HybSys_h
#define HybSys_h
#include <stdint.h>
#include <vector>
#include <functional>

using namespace std;

class Mode {    
 public:
  int dim;
  function<vector<double>(vector<double>)> vectorField;
  
  Mode(int dim, function<vector<double>(vector<double>)> vectorField);
  Mode();
  
  static Mode parallel(Mode M, Mode N);  
};


class Reset {
 public:
  function<bool(vector<double>)> guard;
  function<vector<double>(vector<double>)> reset;
  
  Reset(function<bool(vector<double>)> guard,
	function<vector<double>(vector<double>)> reset);

  static Reset either(Reset r1, Reset r2);
  
  Reset();
};

class HybSys {
 public:
  vector<Mode> modes;
  vector<vector<Reset>> resets;
  
  HybSys(vector<Mode> modes, vector<vector<Reset>> resets);
  HybSys();

  static HybSys parallel(HybSys H, HybSys K);

  static HybSys sequential(HybSys H, HybSys K);

  static HybSys loop(HybSys H);

};

#endif


