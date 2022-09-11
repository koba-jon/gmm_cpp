#ifndef PARAMETER_HPP
#define PARAMETER_HPP

#include <vector>
#include <cstring>


// -------------------
// struct{Parameter}
// -------------------
struct Parameter{

private:

    // Member variable
    bool created = false;
    long long *step = nullptr;

    // Function (read only)
    void alloc_memory(Parameter &param){
        this->free_memory(param);
        param.created = true;
        param.data = new double[param.numel];
        param.step = new long long[param.numst];
        return;
    }
    /********************************************************************************/
    void free_memory(Parameter &param){
        if (param.created){
            param.created = false;
            delete[] param.data;
            delete[] param.step;
        }
        return;
    }
    /********************************************************************************/
    void deep_copy(const Parameter &paramI, Parameter &paramO){

        // Free memory
        this->free_memory(paramO);

        // Copy for variable
        paramO.numel = paramI.numel;
        paramO.numst = paramI.numst;
        paramO.shape = paramI.shape;

        // Copy for array
        if (paramI.created){
            this->alloc_memory(paramO);
            std::memcpy(paramO.data, paramI.data, sizeof(double) * paramO.numel);
            std::memcpy(paramO.step, paramI.step, sizeof(long long) * paramO.numst);
        }

        return;

    }

public:

    double *data = nullptr;
    size_t numel = 0;
    size_t numst = 0;
    std::vector<long long> shape;

    // Constructor
    Parameter() = default;  // default constructor
    Parameter(const Parameter &param){  // copy constructor
        this->deep_copy(param, *this);
    }

    // Operator
    template <class... Args>
    inline Parameter operator()(Args... args_){

        std::vector<long long> args{args_...};
        /*********************************/
        long long idx = 0;
        for (size_t i = 0; i < args.size(); i++){
            idx += args[i] * this->step[i];
        }

        Parameter param;
        std::vector<long long> shape_(this->shape.size() - args.size());
        /*********************************/
        for (size_t i = 0; i < shape_.size(); i++){
            shape_[i] = this->shape[args.size() + i];
        }
        /*********************************/
        param.create({shape_});
        /*********************************/
        for (size_t i = 0; i < param.numel; i++){
            param.data[i] = this->data[idx + i];
        }

        return param;

    }
    /********************************************************************************/
    Parameter operator-(const Parameter &param){  // subtraction operator
        Parameter paramO;
        paramO.create(this->shape);
        for (size_t i = 0; i < this->numel; i++){
            paramO.data[i] = this->data[i] - param.data[i];
        }
        return paramO;
    }
    /********************************************************************************/
    Parameter operator*(const double value){  // multiplication operator
        Parameter paramO;
        paramO.create(this->shape);
        for (size_t i = 0; i < this->numel; i++){
            paramO.data[i] = this->data[i] * value;
        }
        return paramO;
    }
    /********************************************************************************/
    Parameter operator*(const Parameter &param){  // multiplication operator
        Parameter paramO;
        paramO.create(this->shape);
        for (size_t i = 0; i < this->numel; i++){
            paramO.data[i] = this->data[i] * param.data[i];
        }
        return paramO;
    }
    /********************************************************************************/
    Parameter operator/(const double divisor){  // division operator
        Parameter paramO;
        paramO.create(this->shape);
        for (size_t i = 0; i < this->numel; i++){
            paramO.data[i] = this->data[i] / divisor;
        }
        return paramO;
    }
    /********************************************************************************/
    Parameter operator=(const double value){  // assignment operator
        for (size_t i = 0; i < this->numel; i++){
            this->data[i] = value;
        }
        return *this;
    }
    /********************************************************************************/
    Parameter operator=(const Parameter &param){  // assignment operator
        this->deep_copy(param, *this);
        return *this;
    }
    /********************************************************************************/
    Parameter &operator+=(const Parameter &param){  // addition assignment operator
        for (size_t i = 0; i < this->numel; i++){
            this->data[i] += param.data[i];
        }
        return *this;
    }
    /********************************************************************************/
    Parameter &operator/=(const double divisor){  // division assignment operator
        for (size_t i = 0; i < this->numel; i++){
            this->data[i] /= divisor;
        }
        return *this;
    }

    // Function (rewrite)
    void create(const std::vector<long long> &shape_){

        this->numel = 1;
        this->shape = shape_;
        for (long long size : this->shape){
            this->numel *= size;
        }
        this->numst = this->shape.size() - 1;
        this->alloc_memory(*this);

        long long s = 1;
        for (size_t i = this->numst; i > 0; i--){
            s *= this->shape[i];
            this->step[i - 1] = s;
        }

        return;

    }
    /********************************************************************************/
    void create(const std::vector<long long> &shape_, const double value){

        this->numel = 1;
        this->shape = shape_;
        for (long long size : this->shape){
            this->numel *= size;
        }
        this->numst = this->shape.size() - 1;
        this->alloc_memory(*this);

        for (size_t i = 0; i < this->numel; i++){
            this->data[i] = value;
        }

        long long s = 1;
        for (size_t i = this->numst; i > 0; i--){
            s *= this->shape[i];
            this->step[i - 1] = s;
        }

        return;

    }
    /********************************************************************************/
    template <class... Args>
    inline double &at(Args... args_){

        std::vector<long long> args{args_...};
        /*********************************/
        long long idx = 0;
        for (size_t i = 0; i < this->numst; i++){
            idx += args[i] * this->step[i];
        }
        idx += args[this->numst];

        return this->data[idx];

    }
    /********************************************************************************/
    void inplace(const std::vector<long long> &args, const Parameter &param){

        long long idx = 0;
        for (size_t i = 0; i < args.size(); i++){
            idx += args[i] * this->step[i];
        }

        for (size_t i = 0; i < param.numel; i++){
            this->data[idx + i] = param.data[i];
        }

        return;

    }


    // Destructor
    ~Parameter(){
        this->free_memory(*this);
    }

};


#endif