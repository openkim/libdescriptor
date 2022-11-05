Brief Introduction to Enzyme AD
================================

## Fundamentals of Automatic Differentiation (AD)
For a \f$ \mathbb{R}_n \rightarrow \mathbb{R}_m \f$ vector function 
$$
F_{\mathbb{R}_n \rightarrow \mathbb{R}_m}\mathbf{x} = \mathbf{y}
$$
the derivative at \f$ \mathbf{x} \f$ of can be written in terms of Jacobian of the function and Jacobian-vector product (or vector-Jacobian product), i.e.
$$
\frac{ d \mathbf{y}}{d \mathbf{x}}|_\mathbf{x} =  \mathbf{y}^T\mathbf{J}^{\mathbf{F}^{-1}}_{m \times n}
$$
This approach can be translated more efficiently using directional derivatives in direction of \f$ \mathbf{y} \f$. 
Such an approach sidesteps the complexity of computing the full Jacobian matrix. 

## Automatic differentiation 
In Automatic differentiation one can compute the gradients using the Jacobian-vector product of forward function, or using vector-Jacobian product of the inverse function. 
The AD framework such as Enzyme provides means of generating functions which can take input as the resultant vector \f$\mathbf{y}\f$ and return the vector-Jacobian product to yield gradients w.r.t \f$\mathbf{x}\f$. 

The major advantage of AD is the flexibility to take derivatives against any parameters, including hyperparameters, which sometimes are not trivial to analytically estimate. 
This provides means to rapidly iterate and optimize descriptors for more accurate results.
The major disadvantage of AD frameworks is the lost performance due to "trackers" that are used to invert the function.
Enzyme generates highly optimized gradient functions, that offer same order of performance of analytical functions. 
See benchmarks below for evaluation of different AD frameworks, calculating a 64 atom Si configuration using Stillinger-Weber potential for energy, and AD for forces. 
One reason why TorchScript come of somewhat worse is that to avoid any cashing, these benchmarks contain average of 10 one-shot computations. 
Once several warm up runs are allowed, the TorchScript time does come down, but only by a factor of 2 in most optimized cases. So still it is nowhere viable option for large scale computations.

| High lvl frameworks | C++ frameworks |
|---------------------|----------------|
| <img src="./E_F.png" width=400> | <img src="./E_F_C++Updated.png" width=400>|

It is clear that as opposed to analytical function (implemented in OpenKIM, called through ASE) popular frameworks like Pytorch suffers from 4 orders of magnitude extra computation time. While Enzyme is only about 4-6 times slower then pure analytical function. When accounted for bookkeeping in OpenKIM, it is a competitive alternative to analytical functions.

## Enzyme simple examples
Although enzyme offers forward differentiation methods as well, we will focus on reverse diff, as it is currently implemented in `libdescriptor`.
### Enzyme Components
Enzyme works at compile time by generating the gradient functions. 
For identifying which functions to differentiate, Enzyme searches for certain strings. Below is a bare minimal valid Enzyme code
```Cpp
#include <iostream>

// Enzyme arg kinds
int enzyme_dup, enzyme_const, enzyme_out;

// function to diff
double pow(int x, double y){
    double z = 1.0;
    for (int i = 0; i < x; i++){
        z *= y;
    }
    return z;
}

// declaration of diff
double __enzyme_autodiff_d_pow(double (*)(int , double) /* pointer to function to diff */ , 
                                int /* kind of arg */, int /* x */, 
                                int /* kind of arg */, double /* y */);

int main(){
    int x = 3; double y = 4.0;
    // call to gradient
    double d_pow_y =  __enzyme_autodiff_d_pow(pow, enzyme_const, x, enzyme_out, y);
    std::cout << "Function: " << pow(x, y) << "\n";
    std::cout << "Derivative: " << d_pow_y << "\n";
    return 0;
}
```

Enzyme recognizes any function that starts with the identifier `__enzyme_autodiff`.
The arguments of the function includes 

1. Pointer to the function to be differentiated
2. Integer indicating what king of arguments will come next
3. Arguments to be passed on to function, and to be used for diff

Most used integer arguments are of three kind, namely `enzyme_const`, `enzyme_dup`, `enzyme_out`. 
You just need to declare them in your file, and Enzyme will fill out required values.
Their meaning are,

1. `enzyme_const`: The argument that follows is to be considered as a constant, and to be differentiated against. Usually it is used for integers.
2. `enzyme_out`: The argument that follows next is to differentiated against, and results will be returned as a function value.
3. `enzyme_dup`: The following argument will accompany a second "shadow argument", which either will save the output of the derivatives, or contain the vector for vector-jacobian product.

Now to compile the program you need to load the compiler library as
```shell
$ clang++ pow.cpp -Xclang -load -Xclang /usr/local/ClangEnzyme-12.so -O2 -o pow.x
$ ./pow.x
$ Function: 64
Derivative: 48
```

## More Examples
Here is another example of vector functions to demonstrate vector-jacobian products,

```Cpp

#include<vector>
#include <cmath>

int enzyme_dup, enzyme_out, enzyme_const;

void vxp(std::vector<double>& vector_in, double p, int n, std::vector<double>& vector_out){
    for (int i = 0; i < vector_in.size(); i++){
            vector_out[i] = std::pow(vector_in[i], n) * p;
        }
}

void __enzyme_autodiff(
        void (*) (std::vector<double>&, double, int, std::vector<double>& ), //*tr to vec func 
        int, std::vector<double>&, std::vector<double>&, // input vector, grad_vector
        int, double, int, int, // const arguments
        int, std::vector<double>&, std::vector<double>&); //resultant vector for vjp

int main(){
    // Initialize values
    double p = 2.0; int n = 3;
    std::vector<double> vector_in(5), d_vector_in(5), vector_out(5), d_vector_out(5);
    for (int i = 0; i < 5 ; i++){
        vector_in[i] = i; // input vector
        d_vector_in[i] = 0.0; // vector to save gradient w.r.t input
        vector_out[i] = 0.0; // output vector
        d_vector_out[i] = 1.0; // vector to left multiply to jacobian
    }
    __enzyme_autodiff(vxp, 
            enzyme_dup, vector_in, d_vector_in, 
            enzyme_const, p, enzyme_const, n, 
            enzyme_dup, vector_out, d_vector_out);

    std::cout << vector_in[0] << " " << vector_in[4] << "\n";
    std::cout << vector_out[0] << " " << vector_out[4] << "\n";
    std::cout << d_vector_in[0] << " " << d_vector_in[4] << "\n";
    std::cout << d_vector_out[0] << " " << d_vector_out[4] << "\n";
    return 0;
}
```

The main idea to notice in above example is the initialization of `d_vector_out` for unity and passing it along for differentiation. It the the vector that will be left multiplied to the jacobian, or along which direction the directional derivative will be computed. 
In practical usage, it will be the derivative of energy with respect to the descriptors, which will form input to the descriptor gradient function.
Both expected gradient and supplied vector are indicated by shadow argument indicator `enzyme_dup`.