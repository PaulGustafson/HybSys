/*
 * Copyright (C) 2020 Paul Gustafson 
 * License: MIT  (see the file LICENSE)
 */

#define EIGEN_NO_DEBUG
#include <Eigen/Dense>
#include <Eigen/Core>

#ifndef HybSys_h
#define HybSys_h
#include <stdint.h>
#include <vector>
#include <functional>

using namespace std;


class Mode {    
 public:
  int dim;
  function<Eigen::VectorXd(Eigen::VectorXd)> vectorField;
  
  Mode(int dim, function<Eigen::VectorXd(Eigen::VectorXd)> vectorField);
  Mode();
  
  static Mode parallel(Mode M, Mode N);  
};


class Reset {
 public:
  function<bool(Eigen::VectorXd)> guard;
  function<Eigen::VectorXd(Eigen::VectorXd)> reset;
  
  Reset(function<bool(Eigen::VectorXd)> guard,
	function<Eigen::VectorXd(Eigen::VectorXd)> reset);

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

class Semiconjugacy {
  public:
    HybSys dom, cod;
    vector<int> modeMap;
    vector<function<Eigen::VectorXd(Eigen::VectorXd)>> manifoldMap;

    Semiconjugacy(HybSys dom, HybSys cod, vector<int> nodeMap, vector<function<Eigen::VectorXd(Eigen::VectorXd)>> manifoldMap);
    Semiconjugacy();
};


class TAPair {
 public:
  HybSys anchor;
  HybSys temp;

  Semiconjugacy asymptoticPhase;

  TAPair(HybSys anchor, HybSys temp, Semiconjugacy asymptoticPhase);
  TAPair();
};

#endif


