#include "linalg.h"
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <bitset>
#include <cmath>

namespace nm{

    bool isZero(double entry){
        return entry == (double)0;
    }

    bool checkIfValid(double entry){

        //we use 0x7FF for our mask (equal to 0b11111111111, which represents NaN and +/- infinity in exponent)
        //and use fixed-size uint64_t for 64-bit double-precision floating point 
        //(as per IEEE 754 Double-Precision standard)
        const uint64_t exponentMask{0x7FF};

        //technically undefined behaviour, but we leverage the same size
        //to convert into an integer type, where we can use bitwise operations
        uint64_t* uint64Representation {reinterpret_cast<uint64_t*>(&entry)}; 

        //no masking necessary, just read the sign bit
        std::bitset<1> sign         (( *uint64Representation) >> 63); 

        //mask with our exponent mask from above
        std::bitset<11> exponent    (((*uint64Representation) >> 52) & exponentMask);

        //mask with 52 bits of 1s 
        std::bitset<52> mantissa    (((*uint64Representation) >> 0)  & 0xfffffffffffff);

        //if we have nan or inf, return false
        //otherwise true
        return exponent != std::bitset<11>(exponentMask);
    }

    namespace linalg{

        template <std::size_t N> vector<N>::vector(){
            //size is known
            //default constructor
            if(N == 0){
                throw std::invalid_argument("The dimension cannot be 0.");
            }
            if(checkIfValid() == false){
                throw std::invalid_argument("One or more entries contains an infinity or NaN.");
            }
        }

        template <std::size_t N> vector<N>::vector(double* const entries, std::size_t capacity){
            if(N == 0){
                throw std::invalid_argument("The dimension cannot be 0.");
            }
            if(N != capacity){
                throw std::invalid_argument("Capacity of the provided array does not match the vector's size. Provided capacity: " + std::to_string(capacity) + ", vector size: " + std::to_string(N));
            }
            if(entries == nullptr){
                throw std::invalid_argument("The pointer passed may not be nullptr.");
            }
            for(std::size_t i{0}; i < N; ++i){
                if(nm::checkIfValid(entries[i]) == false){
                    throw std::invalid_argument("One or more entries contains an infinity or NaN.");
                }
            }
            for(std::size_t i{0}; i < N; ++i){
                this->_entries[i] = entries[i];
                if(this->_entries[i] != 0){
                    _density++;
                }
            }
        }

        template <std::size_t N> bool vector<N>::checkIfValid() const{
            for(std::size_t i {0}; i < N; ++i){
                double tempDouble {this->_entries[i]};

                //if any value is invalid, the whole vector is invalid
                if(nm::checkIfValid(tempDouble) == false){return false;}
            }
            return true;
        }

        template <std::size_t N> bool vector<N>::isZeroVector() const{
            return _density == N;
            /*for(std::size_t i {0}; i < N; ++i){
                double tempDouble {this->_entries[i]};

                //if any value is non-empty, the whole vector is non-empty
                if(nm::isZero(tempDouble) == false){return false;}
            }
            return true;*/
        }

        template <std::size_t N> void vector<N>::setValue(double entry, std::size_t index){
            if(index >= N){
                throw std::invalid_argument("Index was outside the bounds of the vector capacity. Provided index: " + std::to_string(index) + ", maximum vector index: " + std::to_string(N-1));
            }
            if(nm::checkIfValid(entry) == false){
                throw std::invalid_argument("Entry is NaN or infinity.");
            }
            if(isZero(this->_entries[index]) && isZero(entry) == false){
                //existing zero replaced with non-zero -> density increases
                _density++;
            }else if(isZero(this->_entries[index]) == false && isZero(entry)){
                //existing non-zero replaced with zero -> density decreases
                _density--;
            }
            this->_entries[index] = entry;
        }

        template <std::size_t N> void vector<N>::setValues(double* const entries, std::size_t offset, std::size_t capacity){

            if(entries == nullptr){
                throw std::invalid_argument("The pointer passed may not be nullptr.");
            }

            std::size_t density{0};
            for(std::size_t i{0}; i < capacity; ++i){

                if(nm::checkIfValid(entries[i]) == false){
                    throw std::invalid_argument("Entry at index " + std::to_string(i) + " (of the supplied array) is NaN or infinity.");
                }
                if(nm::isZero(entries[i]) == false){
                    density++;
                }
            }
            
            std::size_t subarrayDensity{0};
            for(std::size_t i{0}; i < capacity; ++i){
                std::size_t offsettedIndex { i + offset - 1};
                if(nm::isZero(this->_entries[offsettedIndex]) == false){
                    subarrayDensity++;
                }
                this->_entries[offsettedIndex] = entries[i];
            }

            //calculate difference (should be int to avoid underflow)
            int delta = subarrayDensity - density;
            _density += delta;

            /*
            hypothetical:
            entries: 1, 5, 3, 4, 5

            entries is of size 5, this->_entries is of size 8, offset is 3
            i = 0 -> i + offset - 1 = 2

             1   3   2   1   0   0   7   3
            [0] [1] [2] [3] [4] [5] [6] [7]

             1   3   1   5   3   4   5   3

            i = 4 -> i + offset - 1 = 4 + 2 = 6

            */
        }

        template <std::size_t N> double vector<N>::getValue(std::size_t index) const{
            if(index >= N){
                throw std::invalid_argument("Index was outside the bounds of the vector capacity. Provided index: " + std::to_string(index) + ", maximum vector index: " + std::to_string(N-1));
            }
            return this->_entries[index];
        }

        template <std::size_t N> const double* const vector<N>::getValues() const{
            return this->_entries;
        }

        template <std::size_t N> std::size_t vector<N>::getDensity() const{
            return _density;
        }

        template <std::size_t N> double vector<N>::getDensityNormalized() const{
            return (double)_density/(double)N;
        }

        template <std::size_t N> vector<N> vector<N>::operator*(double scalar) const{
            double entries[N]{0};
            if(nm::checkIfValid(scalar) == false){
                throw std::invalid_argument("Scalar provided was NaN or infinity.");
            }
            for(std::size_t i{0}; i < N; ++i){
                entries[i] = this->_entries[i] * scalar;
            }
            return vector<N>{entries, N};
        }

        template <std::size_t N> double vector<N>::operator*(vector<N>& other) const{
            double dotProduct{0};
            for(std::size_t i{0}; i < N; ++i){
                dotProduct += (this->_entries[i] * other._entries[i]);
            }
            return dotProduct;
        }

        template <std::size_t N> vector<N> vector<N>::operator+(vector<N>& other) const{
            double entries[N]{0};
            for(std::size_t i{0}; i < N; ++i){
                entries[i] = this->_entries[i] + other._entries[i];
            }
            return vector<N>{entries, N};
        }

        template <std::size_t N> std::size_t vector<N>::dimension() const{
            return N;
        }

        template <std::size_t N> double vector<N>::norm2() const{

            switch(N){
                case 1:
                    return this->_entries[0];
                case 2:
                    return std::hypot(this->_entries[0], this->_entries[1]);
                case 3:
                    return std::hypot(this->_entries[0], this->_entries[1], this->_entries[2]);
                default:double squareSums{0};
                    for(std::size_t i{0}; i < N; ++i){
                        squareSums += (this->_entries[i])*(this->_entries[i]);
                    }
                    return std::sqrt(squareSums);
            }
        }

        template <std::size_t N> double vector<N>::leadingEntry() const{
            for(std::size_t i{0}; i < N; ++i){
                if(isZero(this->_entries[i]) == false){
                    return this->_entries[i];
                }
            }
            return 0;
        }

        template <std::size_t N> std::size_t vector<N>::leadingEntryIndex() const{
            for(std::size_t i{0}; i < N; ++i){
                if(isZero(this->_entries[i]) == false){
                    return i;
                }
            }
            return 0;
        }


        /** MATRIX CLASS **/
        
        template <std::size_t M, std::size_t N> double matrix<M,N>::_square2Determinant(std::size_t row1, std::size_t row2, std::size_t column1, std::size_t column2) const{
            return search(row1,column1) * search(row2,column2) - search(row1,column2) * search(row2,column1);
        }

        template <std::size_t M, std::size_t N> matrix<M,N>::matrix(double* entries, std::size_t capacity, bool columnMajorOrder){
            
            if(entries == nullptr){
                throw std::invalid_argument("The pointer passed may not be nullptr.");
            }
            if(capacity == 0 || M == 0 || N == 0){
                throw std::invalid_argument("Matrix cannot have M = 0, N = 0, or capacity of 0.");
            }


            if(capacity != N*M){
                throw std::invalid_argument("Capacity of the provided array does not match the matrix dimensions. Given capacity: " + std::to_string(capacity) + ", dimensions (N*M): " + std::to_string(N*M));
            }

            this->_entries = entries;
            this->_columnMajorOrder = columnMajorOrder;
            
            //columnMajorOrder means that the entries in the matrix are defined by going down the columns
            //then moving to the right.
            if(columnMajorOrder){

                //we put the entries into column vectors, and put those column vectors in a row
                //we know the amount of column vectors we will need (N)
                //we also know that we will need the column vectors to be of size M
                //we need i to be less then N, and never equal, otherwise we overflow
                for(std::size_t i {0}; i < N; ++i){
                    vector<M> tempColumnVector {};

                    //remember that we rely on truncation to the vector size to avoid overflow
                    //pointer magic to offset
                    //e.g. if M = 3, then we move the array pointer by 3 each time
                    tempColumnVector.setValues(&this->_entries[i*M]); 

                    this->_variables[i] = tempColumnVector;
                }

            }else{

                //now we put the entries into row vectors, and put those row vectors in a column
                //we know the amount of row vectors we will need (M)
                //we know the row vectors will be of size N
                //we need i to be less than M, and never equal, lest we overflow
                for(std::size_t i {0}; i < M; ++i){
                    vector<N> tempRowVector {};

                    //remember that we rely on truncation to the vector size to avoid overflow
                    //pointer magic to offset
                    //e.g. if N = 3, then we move the array pointer by 3 each time
                    tempRowVector.setValues(&this->_entries[i*N]);

                    this->_systems[i] = tempRowVector;

                }
            }
            convert();
        }

        template <std::size_t M, std::size_t N> double matrix<M,N>::search(std::size_t i, std::size_t j) const{
            if(this->_columnMajorOrder){
                return this->_variables[j].getValue(i);
            }else{
                return this->_systems[i].getValue(j);
            }
        }

        template <std::size_t M, std::size_t N> void matrix<M,N>::print() const{
            if(this->_columnMajorOrder){
                //vector<M> _variables[N] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have column vectors, we hold an i constant and iterate j
                //once we finished j-iteration, iterate i and in turn iterate j again
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        std::cout << this->_variables[j].getValue(i) << "\t";
                    }
                    std::cout << "\n";
                }

            }else{
                //vector<N> _systems[M] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have row vectors, we hold a j constant and iterate i
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        std::cout << this->_systems[i].getValue(j) << "\t";
                    }
                    std::cout << "\n";
                }
            }
        }

        template <std::size_t M, std::size_t N> void matrix<M,N>::convert(){
            for(std::size_t i {0}; i < M; ++i){
                for(std::size_t j {0}; j < N; ++j){
                    if(this->_columnMajorOrder){
                        //vector<M> _variables[N] is populated
                        //_variables -> _systems
                        this->_systems[i].setValue(search(i,j),j);
                    }else{
                        //vector<N> _systems[M] is populated
                        //_systems -> _variables
                        this->_variables[j].setValue(search(i,j),i);
                    }
                }
            }
            /*
            std::cout << "original:\n";
            print();
            std::cout << "converted:\n";
            if(!this->_columnMajorOrder){
                //vector<M> _variables[N] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have column vectors, we hold an i constant and iterate j
                //once we finished j-iteration, iterate i and in turn iterate j again
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        std::cout << this->_variables[j].getValue(i) << "\t";
                    }
                    std::cout << "\n";
                }

            }else{
                //vector<N> _systems[M] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have row vectors, we hold a j constant and iterate i
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        std::cout << this->_systems[i].getValue(j) << "\t";
                    }
                    std::cout << "\n";
                }
            }
            */
        }

        template <std::size_t M, std::size_t N> bool matrix<M,N>::checkIfValidMatrix() const{
            if(this->_columnMajorOrder){
                for(std::size_t i{0}; i < N; ++i){
                    if(this->_variables[i].checkIfValid() == false){
                        return false;
                    }
                }
            }else{
                for(std::size_t i{0}; i < M; ++i){
                    if(this->_systems[i].checkIfValid() == false){
                        return false;
                    }
                }
            }
            return true;
        }

        template <std::size_t M, std::size_t N> bool matrix<M,N>::isRowEmpty(std::size_t row) const{
            if(row >= M){
                throw std::invalid_argument("Row specified is outside the matrix's bounds.");
            }
            return this->_systems[row].isZeroVector();
        }

        template <std::size_t M, std::size_t N> bool matrix<M,N>::isColumnEmpty(std::size_t column) const{
            if(column >= N){
                throw std::invalid_argument("Column specified is outside the matrix's bounds.");
            }
            return this->_variables[column].isZeroVector();
        }

        template <std::size_t M, std::size_t N> double matrix<M,N>::leadingEntry(std::size_t row) const{
            if(row >= M){
                throw std::invalid_argument("Row specified is outside the matrix's bounds.");
            }
            return _systems[row].leadingEntry();
        }

        /*template <std::size_t M, std::size_t N> double matrix<M,N>::determinant() const{

            if(N != M){
                throw std::invalid_argument("The matrix must be square (i.e. # of rows (M) is equal to # of columns (N)). Given matrix of size: " + std::to_string(M) + " (M), " + std::to_string(N) + " (N)");
            }
            
            //simplest case: 1x1 matrix
            if(N == 1){
                return this->_entries[0];
            }else if(N == 2){
                return _square2Determinant(0,1,0,1);
            }

            return 0;
        }*/

        template <std::size_t M, std::size_t N> matrix<M,N> matrix<M,N>::operator+(matrix<M,N>& other) const{
            double zeros[M*N]{0};
            if(other._columnMajorOrder && this->_columnMajorOrder){
                for(std::size_t i{0}; i < M*N; ++i){
                    zeros[i] = other._entries[i] + this->_entries[i];
                }
            }
            matrix<M,N> output {zeros, M*N, true};
            return output;
        }
        
        template <std::size_t M, std::size_t N> template <std::size_t K> matrix<M,K> matrix<M,N>::operator*(matrix<N,K>& other) const{
		//algorithm: (naive)
		//this matrix (A) should be accessed in row-major order
		//other matrix (B) should be accessed in column-major order
		//transpose B -> B^T (now KxN matrix)
		matrix<K,N> B_transpose = ~other;
		//now we can access B in row-major order
		
		
	}

        template <std::size_t M, std::size_t N> matrix<N,M> matrix<M,N>::operator~() const{
            matrix<N,M> transpose {this->_entries, N*M, !this->_columnMajorOrder};
            return transpose;
        }
        

    }
}

bool vectorTests(){
    std::size_t score{0};

    double array3[3] = {1,2.5,5};
    
    std::cout << "Testing correct array size\n";
    try{
        nm::linalg::vector<3> test{array3, 3};
        std::cout << "Success!\n";
        score++;
    }catch(std::invalid_argument message){
        std::cout << message.what();
    }
    std::cout << "\n";

    std::cout << "Testing incorrect declared array size\n";
    try{
        nm::linalg::vector<3> test{array3, 4};
    }catch(std::invalid_argument message){
        std::cout << "Success!\n";
        std::cout << message.what();
        score++;
    }
    std::cout << "\n";


    return score == 2;
}




int main(){

    //vectorTests();
    /*
    int blah{0};
    //std::cin >> blah;
    double entries[12] {1,2,3,4,5,6,7,8,9,10,11,12};
    nm::linalg::matrix<3,4> ma {entries, 12, (blah == 0 ? true : false)};
    ma.print();
    std::cout << "\n";
    nm::linalg::matrix<4,3> tr {~ma};
    tr.print();
    std::size_t i{0}; std::size_t j{0};

    ma.convert();
*/
    double entries1 [12] {1,2,3,4,5,6,7,8,9,10,11,12};
    double entries2 [12] {12,11,10,9,8,4,6,5,4,3,2,1};
    double zeros [12] {0};
    /*
    nm::linalg::matrix<3,4> ma1 {entries1, 12, true};
    nm::linalg::matrix<3,4> ma2 {entries2, 12, true};
    nm::linalg::matrix<3,4> sum {zeros, 12, true};
    sum = ma1 + ma2;
    sum.print();*/

    double entries3 [3] {1,2,3};
    nm::linalg::vector<3> initial {entries3, 3};
    nm::linalg::vector<3> final {initial*3};

    for(std::size_t i{0}; i < 3; ++i){
        std::cout << "init dim " << i << ": " << initial.getValue(i) << ", final dim " << i << ": " << final.getValue(i) << "\n"; 
    }
    //initial*final;
    std::cout << "dot product: " << initial*final << ", vector sum: [ " << (initial+final).getValue(0) << " " << (initial+final).getValue(1) << " " << (initial+final).getValue(2) << " ]\n";
    std::cout << "dimension of init: " << initial.dimension() << ", 2-norm of initial: " << initial.norm2() << ", 2-norm of zero vector " << (nm::linalg::vector<12>{zeros, 12}).norm2() << ", norm of [1,2,3,4,5,6,7,8,9,10,11,12] " << (nm::linalg::vector<12>{entries1, 12}).norm2() << "\n";
    //std::cerr << "pingas\n";

    nm::linalg::matrix<3,4> ma3 {entries1, 12, true};

    std::cout << sizeof(ma3) << "\n";
}