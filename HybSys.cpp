/*
 * Copyright (C) 2020 Paul Gustafson - All Rights Reserved
 * Unauthorized copying of this file, via any medium is strictly prohibited
 * Proprietary and confidential
 */


#include "HybSys.h"

using namespace std;


template <typename T> vector<T> concat(vector<T> &a, vector<T> &b) {
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
    function<vector<double>(vector<double>)> vectorField = [](vector<double> coords) {
      vector<double> mCoords = vector<double>(coords.begin(), coords.begin() + M.dim);
      vector<double> nCoords = vector<double>(coords.begin() + M.dim,
					    coords.begin() + M.dim + N.dim)
    
      mAns = M.vectorField(mCoords);
      nAns = N.vectorField(nCoords);
    
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
  n = K.modes.size();

  for (int i = 0; i < H.modes.size(); i++) {
    for (int j = 0; j < K.modes.size(); j++) {
      modes[i*n + j] = parallel(H.modes[i], K.modes[j]);
    }
  }

  
  for (int i = 0; i < H.modes.size(); i++) {
    for (int j = 0; j < K.modes.size(); j++) {
      resets[i*n + j]  = vector<Reset> (H.modes.size() * K.modes.size());
      for (int k = 0; k < H.modes.size(); k++) {
	for (int l = 0; l < K.modes.size(); l++) {
	  Reset rh = H.resets[i][k];
	  Mode domRh = H.modes[i];
	  Guard gh = rh.guard;
	  Resets rk = K.resets[j][l];
	  Guard gk = rk.guard;

	  function<bool>(vector<double>) guard1, guard2, guard3;
	  function<vector<double>(vector<double>) reset1, reset2, reset3;
	  
	  guard1 = [](vector<double> x) {
	     xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	     xk = vector<double>(x.begin() + domRh.dim, x.end());
	     return gh(xh) && !gk(xk); 
	   };
	  reset1 = [](vector<double> x) {
	    xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	    xk = vector<double>(x.begin() + domRh.dim, x.end());
	    return concat(rh(xh), xk);
	  };
	  resets[i*n + j][k*n + j] = Reset(guard1, reset1);
	   
	   
	  guard2 = [](vector<double> x) {
	    xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	    xk = vector<double>(x.begin() + domRh.dim, x.end());
	    return !gh(xh) && gk(xk);
	  };
	  reset2 = [](vector<double> x){
	     xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	     xk = vector<double>(x.begin() + domRh.dim, x.end());
	     return concat(xh, rk(xk));
	  };
	  resets[i*n + j][i*n +l] = Reset(guard2, reset2);
	   
	  guard3 = [](vector<double> x){
	     xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	     xk = vector<double>(x.begin() + domRh.dim, x.end());
	     return gh(xh) && gk(xk);
	  };
	  reset3 = [](vector<double> x){
	    xh = vector<double>(x.begin(), x.begin() + domRh.dim);
	    xk = vector<double>(x.begin() + domRh.dim, x.end());
	    return concat(rh(xh), rk(xk));
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
