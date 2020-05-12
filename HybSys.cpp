/*
 * Copyright (C) 2020 Paul Gustafson - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 */


#include "HybSys.h"

using namespace std;


template <typename T> vector<T> concat(vector<T> a, vector<T> b) {
    vector<T> ret = vector<T>();
    copy(a.begin(), a.end(), back_inserter(ret));
    copy(b.begin(), b.end(), back_inserter(ret));
    return ret;
}

Mode::Mode(int dim, function<vector<double>(vector<double>)> vectorField) {
  this->dim = dim;
  this->vectorField = vectorField;  
}


Mode Mode::parallel(Mode M, Mode N) {
  function<vector<double>(vector<double>)> vectorField = [&M, &N](vector<double> coords) {
      vector<double> mCoords = vector<double>(coords.begin(), coords.begin() + M.dim);
      vector<double> nCoords = vector<double>(coords.begin() + M.dim,
					      coords.begin() + M.dim + N.dim);
    
      vector<double> mAns = M.vectorField(mCoords);
      vector<double> nAns = N.vectorField(nCoords);
    
      return concat(mAns, nAns);
    };

  return Mode(M.dim + N.dim, vectorField);
}

Reset::Reset(function<bool(vector<double>)> guard,
	     function<vector<double>(vector<double>)> reset) {
  this->guard = guard;
  this->reset = reset;
}

HybSys::HybSys(vector<Mode> modes, vector<vector<Reset>> resets) {
  this->modes = modes;
  this->resets = resets;
}


HybSys HybSys::parallel(HybSys H, HybSys K) {
  vector<Mode> modes (H.modes.size() * K.modes.size());
  vector<vector<Reset>> resets (H.modes.size() * K.modes.size());
  int n = K.modes.size();

  for (int i = 0; i < H.modes.size(); i++) {
    for (int j = 0; j < K.modes.size(); j++) {
      modes[i*n + j] = Mode::parallel(H.modes[i], K.modes[j]);
    }
  }

  
  for (int i = 0; i < H.modes.size(); i++) {
    for (int j = 0; j < K.modes.size(); j++) {
      resets[i*n + j]  = vector<Reset> (H.modes.size() * K.modes.size());
      for (int k = 0; k < H.modes.size(); k++) {
	for (int l = 0; l < K.modes.size(); l++) {
	  int offset = H.modes[i].dim;
	  Reset rh = H.resets[i][k];
	  Reset rk = K.resets[j][l];

	  function<bool(vector<double>)> guard1, guard2, guard3;
	  function<vector<double>(vector<double>)> reset1, reset2, reset3;
	  
	  guard1 = [&offset, &rh, &rk](vector<double> x) {
	     vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	     vector<double> xk = vector<double>(x.begin() + offset, x.end());
	     return rh.guard(xh) && !rk.guard(xk); 
	   };
	  reset1 = [&offset, &rh, &rk](vector<double> x) {
	    vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	    vector<double> xk = vector<double>(x.begin() + offset, x.end());
	    return concat(rh.reset(xh), xk);
	  };
	  resets[i*n + j][k*n + j] = Reset(guard1, reset1);
	   
	   
	  guard2 = [&offset, &rh, &rk](vector<double> x) {
	    vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	    vector<double> xk = vector<double>(x.begin() + offset, x.end());
	    return !rh.guard(xh) && rk.guard(xk);
	  };
	  reset2 = [&offset, &rh, &rk](vector<double> x){
	     vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	     vector<double> xk = vector<double>(x.begin() + offset, x.end());
	     return concat(xh, rk.reset(xk));
	  };
	  resets[i*n + j][i*n +l] = Reset(guard2, reset2);
	   
	  guard3 = [&offset, &rh, &rk](vector<double> x){
	     vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	     vector<double> xk = vector<double>(x.begin() + offset, x.end());
	     return rh.guard(xh) && rk.guard(xk);
	  };
	  reset3 = [&offset, &rh, &rk](vector<double> x){
	    vector<double> xh = vector<double>(x.begin(), x.begin() + offset);
	    vector<double> xk = vector<double>(x.begin() + offset, x.end());
	    return concat(rh.reset(xh), rk.reset(xk));
	  };
	  resets[i*n + j][k*n + l] = Reset(guard3, reset3);
	}
      }
    }
  }
  return HybSys(modes, resets);
}

Mode::Mode() {
}

Reset::Reset() {
}

HybSys::HybSys() {
}
