/*
 * Copyright (C) 2020 Paul Gustafson
 * License: MIT  (see the file LICENSE)
 */


#include "HybSys.h"
#include <Eigen/Dense>

using namespace std;


Eigen::VectorXd concatE(Eigen::VectorXd a, Eigen::VectorXd b) {
    Eigen::VectorXd vec_joined(a.size() + b.size());
    vec_joined << a, b;
    return vec_joined;
}

template <typename T> vector<T> concat(vector<T> a, vector<T> b) {
    vector<T> ret = vector<T>();
    copy(a.begin(), a.end(), back_inserter(ret));
    copy(b.begin(), b.end(), back_inserter(ret));
    return ret;
}


Mode::Mode(int dim, function<Eigen::VectorXd(Eigen::VectorXd)> vectorField) {
  this->dim = dim;
  this->vectorField = vectorField;  
}


Mode Mode::parallel(Mode M, Mode N) {
  function<Eigen::VectorXd(Eigen::VectorXd)> vectorField = [&M, &N](Eigen::VectorXd coords) {
      Eigen::VectorXd mCoords = Eigen::VectorXd(0, M.dim);
      Eigen::VectorXd nCoords = Eigen::VectorXd(M.dim, M.dim + N.dim);
    
      Eigen::VectorXd mAns = M.vectorField(mCoords);
      Eigen::VectorXd nAns = N.vectorField(nCoords);
    
      return concatE(mAns, nAns);
    };

  return Mode(M.dim + N.dim, vectorField);
}

Reset::Reset(function<bool(Eigen::VectorXd)> guard,
	     function<Eigen::VectorXd(Eigen::VectorXd)> reset) {
  this->guard = guard;
  this->reset = reset;
}

HybSys::HybSys(vector<Mode> modes, vector<vector<Reset>> resets) {
  this->modes = modes;
  this->resets = resets;
}

Semiconjugacy::Semiconjugacy(HybSys dom, HybSys cod, vector<int> nodeMap, vector<function<Eigen::VectorXd(Eigen::VectorXd)>> manifoldMap) {
  this->dom = dom;
  this->cod = cod;
  this->modeMap = modeMap;
  this->manifoldMap = manifoldMap;
}

TAPair::TAPair(HybSys anchor, HybSys temp, Semiconjugacy asymptoticPhase) { 
  this->anchor = anchor;
  this->temp = temp;
  this->asymptoticPhase = asymptoticPhase;
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

	  function<bool(Eigen::VectorXd)> guard1, guard2, guard3;
	  function<Eigen::VectorXd(Eigen::VectorXd)> reset1, reset2, reset3;
	  
	  guard1 = [&offset, &rh, &rk](Eigen::VectorXd x) {
	     Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	     Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	     return rh.guard(xh) && !rk.guard(xk); 
	   };
	  reset1 = [&offset, &rh, &rk](Eigen::VectorXd x) {
	    Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	    Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	    return concatE(rh.reset(xh), xk);
	  };
	  resets[i*n + j][k*n + j] = Reset(guard1, reset1);
	   
	   
	  guard2 = [&offset, &rh, &rk](Eigen::VectorXd x) {
	    Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	    Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	    return !rh.guard(xh) && rk.guard(xk);
	  };
	  reset2 = [&offset, &rh, &rk](Eigen::VectorXd x){
	     Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	     Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	     return concatE(xh, rk.reset(xk));
	  };
	  resets[i*n + j][i*n +l] = Reset(guard2, reset2);
	   
	  guard3 = [&offset, &rh, &rk](Eigen::VectorXd x){
	     Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	     Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	     return rh.guard(xh) && rk.guard(xk);
	  };
	  reset3 = [&offset, &rh, &rk](Eigen::VectorXd x){
	    Eigen::VectorXd xh = Eigen::VectorXd(0, offset);
	    Eigen::VectorXd xk = Eigen::VectorXd(offset, x.size());
	    return concatE(rh.reset(xh), rk.reset(xk));
	  };
	  resets[i*n + j][k*n + l] = Reset(guard3, reset3);
	}
      }
    }
  }
  return HybSys(modes, resets);
}

/*
 * Replace the final mode of H with the first mode of K
 */
HybSys HybSys::sequential(HybSys H, HybSys K) {
  H.modes.pop_back();
  vector<Mode> modes = concat(H.modes, K.modes);
  vector<vector<Reset>> resets (H.modes.size() + K.modes.size() - 1);

  for (int i = 0; i < resets.size(); i++) {
    resets[i] = vector<Reset>  (H.modes.size() + K.modes.size() - 1);
    if (i < H.modes.size() - 1) {
      resets[i] = H.resets[i];
    }
    else { 
      resets[i + H.modes.size() - 1] = K.resets[i - H.modes.size() + 1];
    }
  }
  
  return HybSys(modes, resets);  
}

Reset Reset::either(Reset r1, Reset r2) {
  function<bool(Eigen::VectorXd)> guard = [r1, r2](Eigen::VectorXd x) {
    return r1.guard(x) || r2.guard(x);
  };
  function<Eigen::VectorXd(Eigen::VectorXd)> reset = [r1, r2](Eigen::VectorXd x) {
    if (r1.guard(x)) {
      return r1.reset(x);
    }
    return r2.reset(x);
  };
  return Reset(guard, reset);
}

/*
 * Identify the final mode and the first mode.  Use the first mode's guards/resets
 */
HybSys HybSys::loop(HybSys H) {
  H.modes.pop_back();
  H.resets.pop_back();
  int last = H.modes.size() - 1;

  // If you would reset into the final mode, reset into the first mode instead
  for (int i = 0; i < H.modes.size(); i++) {
    H.resets[i][0] = Reset::either(H.resets[i][0], H.resets[i][last]);
    H.resets[i].pop_back();
  }
  return H;
}

Mode::Mode() {
}

Reset::Reset() {
}

HybSys::HybSys() {
}

Semiconjugacy::Semiconjugacy() {
}

TAPair::TAPair() {
}

