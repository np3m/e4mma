/***********************************************************************
 * Copyright (c) 2016 Luke F. Roberts.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 ***********************************************************************/
/// \file Tensor.hpp
/// \author lroberts
/// \since Apr 02, 2016
///
/// \brief
///
///

#ifndef TENSOR_HPP_
#define TENSOR_HPP_

#include <iostream>
#include <vector> 

namespace nuopac { 

  template<class T> 
  struct Tensor; 

  template<class T> 
  T FullContract(Tensor<T> a, Tensor<T> b);

  template<class T> 
  Tensor<T> Add(Tensor<T> a, Tensor<T> b);

  template<class T> 
  Tensor<T> PartialContract(Tensor<T> a, Tensor<T> b, int idxA=1,
                            int idxB=0);

  /// Class for storing tensors decomposed relative to the energy momentum 
  /// transfer vector
  template<class T> 
  struct Tensor {
    T L, Q, Tp, Tm, Mp, Mm;
  
    Tensor() 
      : L(0.0), Q(0.0), Tp(0.0), Tm(0.0), Mp(0.0), Mm(0.0) {}
    Tensor(T l, T q, T tp, T tm, T mp, T mm) 
      : L(l), Q(q), Tp(tp), Tm(tm), Mp(mp), Mm(mm) {}
  
    void Multiply(T f) {L*=f; Q*=f; Tp*=f; Tm*=f; Mp*=f; Mm*=f;} 
    std::vector<T> GetVector() {return {L, Q, Tp, Tm, Mp, Mm};}
    std::vector<T> GetLVector() {return {L, Q, Mp, Mm};}
    std::vector<T> GetTVector() {return {Tp, Tm};}
  
    friend T FullContract<T>(Tensor a, Tensor b); 
    friend Tensor<T> Add<T>(Tensor a, Tensor b); 
    friend Tensor<T> PartialContract<T>(Tensor<T> a, Tensor<T> b, 
                                        int idxA, int idxB); 

  };

  // Contract all indices of two tensors
  template<class T> 
  T FullContract(Tensor<T> a, Tensor<T> b) {
    return a.L*b.L + a.Q*b.Q + 2.0*a.Tp*b.Tp + 2.0*a.Tm*b.Tm 
      -2.0*a.Mp*b.Mp - 2.0*a.Mm*b.Mm;
  }

  template<class T> 
  Tensor<T> Add(Tensor<T> a, Tensor<T> b) {
    return Tensor<T>(a.L + b.L, a.Q + b.Q, a.Tp + b.Tp, a.Tm + b.Tm, 
                     a.Mp + b.Mp, a.Mm + b.Mm);
  }

  // Perform the contraction C^{\alpha \beta} = A^{\alpha \mu} B_{\mu}^{\beta}
  // idxA and idxB choose which indices of A and B to contract, by default 
  // 1 and 0, respectively.  
  template<class T> 
  Tensor<T> PartialContract(Tensor<T> a, Tensor<T> b, int idxA, int idxB) {
    
    // This accounts for changing sign on assymetric parts when indices are 
    // switched 
    double sgna = 1.0;
    double sgnb = 1.0;
    
    if (idxA == 0) sgna = -1.0;
    if (idxB == 1) sgnb = -1.0;
    T cQ = a.Q*b.Q - a.Mp*b.Mp + sgnb*a.Mp*b.Mm 
      - sgna*a.Mm*b.Mp + sgna*sgnb*a.Mm*b.Mm; 
    T cL = a.L*b.L - a.Mp*b.Mp - sgnb*a.Mp*b.Mm 
      + sgna*a.Mm*b.Mp + sgna*sgnb*a.Mm*b.Mm; 
    T cTp = a.Tp*b.Tp - sgna*sgnb*a.Tm*b.Tm; 
    T cTm = sgnb*a.Tp*b.Tm - sgna*a.Tm*b.Tp; 
    T cMp = 0.5*(a.Q*(b.Mp + sgnb*b.Mm) + a.L*(b.Mp - sgna*b.Mm)
                 + b.Q*(a.Mp - a.Mm) + b.L*(a.Mp + a.Mm));
    T cMm = 0.5*(a.Q*(b.Mp + sgnb*b.Mm) - a.L*(b.Mp - sgnb*b.Mm)
                 + b.Q*(a.Mp - sgna*a.Mm) - b.L*(a.Mp + sgna*a.Mm));
    return Tensor<T>(cL, cQ, cTp, cTm, cMp, cMm);
  }

}

#endif // TENSOR_HPP_
