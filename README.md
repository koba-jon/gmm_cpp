# GMM C++ Samples
These are sample programs of Gaussian Mixture Model from scratch in C++.

## 1. Requirement

### Boost

This is used for command line arguments, etc. <br>
~~~
$ sudo apt install libboost-dev libboost-all-dev
~~~

## 2. Preparation

### Git Clone
~~~
$ git clone https://github.com/koba-jon/gmm_cpp.git
$ cd svm_cpp
~~~

## 3. Usage

### 3.1. Build
Please build the source file according to the procedure.
~~~
$ mkdir build
$ cd build
$ cmake ..
$ make
$ cd ..
~~~

### 3.2. Dataset Setting

The following hierarchical relationships are recommended.


### 3.3. Execution

The following is an example for Toy Dataset.

#### Setting
Please set the shell for executable file.
~~~
$ vi scripts/toy.sh
~~~
If you want to view specific examples of command line arguments, please view "src/main.cpp" or add "--help" to the argument.
~~~
#!/bin/bash

DATA='toy'

./GMM \
    --dataset ${DATA} \
    --D 3 \
    --K 8
~~~

#### Run
Please execute the following to start the program.
~~~
$ sh scripts/toy.sh
~~~


## 4. License

This repository: [MIT License](LICENSE)

### 3rd-Party Libraries
- Boost <br>
Official : https://www.boost.org/ <br>
License : https://www.boost.org/users/license.html <br>

## References
- https://qiita.com/kenmatsu4/items/59ea3e5dfa3d4c161efb (Japanese)
- https://qiita.com/mnanri/items/f0e9b20395545dd674c9 (Japanese)
