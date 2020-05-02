/*
 * Copyright (C) 2020 Paul Gustafson - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 */

#ifndef HybSys_h
#define HybSys_h

#include <stdint.h>
#include <vector>

using namespace std;

class Mode {    
 public:
  int dim;
  function<vector<double>(vector<double>)> vectorField;
  
  Mode(int dim, function<vector<double>(vector<double>)> vectorField);
    
  static Mode parallel(Mode M, Mode N);
};


class Reset {
 public:
  function<bool(vector<double>)> guard;
  function<vector<double>(vector<double>)> reset;
  
  Reset(function<bool(vector<double>)> guard,
	function<vector<double>(vector<double>)> reset);
};

class HybSys {
 public:
  vector<Mode> modes;
  vector<vector<Reset>> resets;
  
  HybSys(vector<Mode> modes, vector<vector<Reset>> resets);

  static HybSys parallel(HybSys H, HybSys K);
};



#endif
