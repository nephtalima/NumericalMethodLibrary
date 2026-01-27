#include "linalg.h"
#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <cstdint>
#include <bitset>
#include <cmath>

namespace nm{



    

    bool isZero(double entry, double tolerance){
        return isEqual(entry, 0, tolerance);
    }

    bool isEqual(double a, double b, double tolerance){
        if(checkIfValid(a) == false){
            throw std::invalid_argument("The provided argument a is invalid.");
        }
        if(checkIfValid(b) == false){
            throw std::invalid_argument("The provided argument b is invalid.");
        }
        if(checkIfValid(tolerance) == false){
            throw std::invalid_argument("The provided argument tolerance is invalid.");
        }
        //std::cout << "a: " << a << ", b: " << b << "\n" << "a + tolerance: " << a + tolerance << ", b - tolerance: " << b - tolerance << "\n";
        if(a + tolerance > b && a - tolerance < b){
            return true;
        }
        return false;
    }

    bool checkIfValid(double entry){

        //Inspired by code from Douglas Wilhelm Harder (cool prof)
        //https://ece.uwaterloo.ca/~dwharder/nm/Algorithms/Formated_printing/src/float_rep.tpp

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

    void printBitsD(double entry){
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

        std::cout << "Sign: " << sign.to_string() << '\n' << "Exponent: " << exponent.to_string() << '\n' << "Mantissa: " << mantissa.to_string() << '\n';
        return;

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

        template <std::size_t N> vector<N>::vector(double* const entries, std::size_t capacity, double tolerance){
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
                if(nm::isZero(this->_entries[i], tolerance) == false){
                    _density++;
                    //std::cout << "density: " << _density << "\n";
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

        template <std::size_t N> void vector<N>::print() const{
            std::cout << " [";
            for(std::size_t i {0}; i < N; ++i){
                std::cout << " " << this->_entries[i];
            }
            std::cout << " ]\n";
        }

        template <std::size_t N> bool vector<N>::isZeroVector() const{
            //std::cout << "density: " << _density << "\n";
            return _density == 0;
            /*for(std::size_t i {0}; i < N; ++i){
                double tempDouble {this->_entries[i]};

                //if any value is non-empty, the whole vector is non-empty
                if(nm::isZero(tempDouble) == false){return false;}
            }
            return true;*/
        }

        template <std::size_t N> void vector<N>::setValue(double entry, std::size_t index, double tolerance){
            if(index >= N){
                throw std::invalid_argument("Index was outside the bounds of the vector capacity. Provided index: " + std::to_string(index) + ", maximum vector index: " + std::to_string(N-1));
            }
            if(nm::checkIfValid(entry) == false){
                throw std::invalid_argument("Entry is NaN or infinity.");
            }
            //std::cout << "is existing zero? " << isZero(this->_entries[index], tolerance) << ", is new entry zero? " << isZero(entry, tolerance) << "\n";
            //std::cout << "old density: " << _density << "\n";
            if(isZero(this->_entries[index], tolerance) && isZero(entry, tolerance) == false){
                //existing zero replaced with non-zero -> density increases
                _density++;
            }else if(isZero(this->_entries[index], tolerance) == false && isZero(entry, tolerance)){
                //existing non-zero replaced with zero -> density decreases
                _density--;
            }
            //std::cout << "new density: " << _density << "\n";
            this->_entries[index] = entry;
        }

        template <std::size_t N> void vector<N>::setValues(double* const entries, std::size_t offset, std::size_t capacity, double tolerance){

            if(entries == nullptr){
                throw std::invalid_argument("The pointer passed may not be nullptr.");
            }

            std::size_t overlap = std::min(N - offset, capacity);
            //std::cout << "overlap: " << overlap << "\n";
            std::size_t subarrayDensity {0};
            for(std::size_t i{offset}; i < overlap + offset; ++i){
                std::size_t offsetIndex {i-offset};
                if(isZero(entries[offsetIndex]) == false) subarrayDensity++;
                this->_entries[i] = entries[offsetIndex];
            }

            //std::cout << "subarray density: " << subarrayDensity << ", _density: " << _density << "\n";

            _density = N - overlap + subarrayDensity;
            //std::cout << "subarray density: " << subarrayDensity << ", _density: " << _density << "\n";

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

        template <std::size_t N> vector<N> vector<N>::operator-(vector<N>& other) const{
            double entries[N]{0};
            //std::cout << "other vector: \n"; 
            //other.print();
            for(std::size_t i{0}; i < N; ++i){
                entries[i] = this->_entries[i] - other._entries[i];
            }
            return vector<N>{entries, N};
        }

        template <std::size_t N> bool vector<N>::operator==(vector<N>& other) const{
            return this->isEqual(other);
        }

        template <std::size_t N> bool vector<N>::operator!=(vector<N>& other) const{
            return this->isEqual(other) == false;
        }

        template <std::size_t N> bool vector<N>::isEqual(vector<N>& other, double tolerance) const {
            for(std::size_t i {0}; i < N; ++i){
                if(nm::isEqual(this->_entries[i],other._entries[i]) == false)return false;
            }
            return true;
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

        template <std::size_t M, std::size_t N> std::size_t matrix<M,N>::calculateDensity() const{
            /*std::size_t density{0};
            
            // shamelessly stolen from print()
            if(this->_columnMajorOrder){
                //vector<M> _variables[N] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have column vectors, we hold an i constant and iterate j
                //once we finished j-iteration, iterate i and in turn iterate j again
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        if(nm::isZero(this->_variables[j].getValue(i)) == false) density++;
                    }
                    //std::cout << "\n";
                }

            }else{
                //vector<N> _systems[M] is populated
                //i,j for M,N => i iterates M, j iterates N
                //since we have row vectors, we hold a j constant and iterate i
                for(std::size_t i{0}; i < M; ++i){
                    for(std::size_t j{0}; j < N; ++j){
                        if(nm::isZero(this->_systems[i].getValue(j)) == false) density++;
                    }
                    //std::cout << "\n";
                }
            }

            std::cout << "testing if vector density can be used instead of the above code:\n";
            std::size_t densityA{0};
            std::size_t densityB{0};
            for(std::size_t i {0}; i < M; ++i){
                densityA += this->_systems[i].getDensity();
            } 
            for(std::size_t j{0}; j < N; ++j){
                densityB += this->_variables[j].getDensity();
            }
            //std::cout << "original: " << density << '\n';
            std::cout << "alt calc (rows): " << densityA << "\n";
            std::cout << "alt calc (columns): " << densityB << "\n";
            */
            
            if(this->_columnMajorOrder){
                //use _variables
                std::size_t density{0};
                for(std::size_t i {0}; i < N; ++i){
                    density += this->_variables[N].getDensity();
                }
                return density;
            }else{
                //use _systems
                std::size_t density{0};
                for(std::size_t i{0}; i < M; ++i){
                    density += this->_systems[M].getDensity();
                }
                return density;
            }


            return 0;
        }

        //what if supplied array is changed after initialization?
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
                    //remember that we rely on truncation to the vector size to avoid overflow
                    //pointer magic to offset
                    //e.g. if M = 3, then we move the array pointer by 3 each time
                    this->_variables[i].setValues(&this->_entries[i*M]);
                }

            }else{

                //now we put the entries into row vectors, and put those row vectors in a column
                //we know the amount of row vectors we will need (M)
                //we know the row vectors will be of size N
                //we need i to be less than M, and never equal, lest we overflow
                for(std::size_t i {0}; i < M; ++i){
                    //remember that we rely on truncation to the vector size to avoid overflow
                    //pointer magic to offset
                    //e.g. if N = 3, then we move the array pointer by 3 each time
                    this->_systems[i].setValues(&this->_entries[i*N]);

                }
            }

            convert();
        }

        template <std::size_t M, std::size_t N> std::size_t matrix<M,N>::getCapacity() const{
            return M*N;
        }

        template <std::size_t M, std::size_t N> double matrix<M,N>::search(std::size_t rowIndex, std::size_t columnIndex) const{

            
            if(rowIndex >= M){
                throw std::invalid_argument("Index rowIndex must be strictly less than " + std::to_string(M));
            }
            if(columnIndex >= N){
                throw std::invalid_argument("Index columnIndex must be strictly less than " + std::to_string(N));
            }

            if(this->_columnMajorOrder){
                return this->_variables[columnIndex].getValue(rowIndex);
            }else{
                return this->_systems[rowIndex].getValue(columnIndex);
            }
        }

        template <std::size_t M, std::size_t N> void matrix<M,N>::setEntry(std::size_t rowIndex, std::size_t columnIndex, double value){
            if(rowIndex >= M){
                throw std::invalid_argument("Index rowIndex must be strictly less than " + std::to_string(M));
            }
            if(columnIndex >= N){
                throw std::invalid_argument("Index columnIndex must be strictly less than " + std::to_string(N));
            }

            //do not replace .setValue with raw array accesses; density calculations
            if(this->_columnMajorOrder){
                this->_variables[columnIndex].setValue(value, rowIndex);
            }else{
                this->_systems[rowIndex].setValue(value,columnIndex);
            }
            this->convert();
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

        template <std::size_t M, std::size_t N> vector<N> matrix<M, N>::getSystemByIndex(std::size_t index) const{
            if(index >= M){
                throw std::invalid_argument("Index must be strictly less than " + std::to_string(M));
            }
            return _systems[index];
        }

        template <std::size_t M, std::size_t N> vector<N> matrix<M, N>::getRowByIndex(std::size_t index) const{
            return this->getSystemByIndex(index);
        }

        template <std::size_t M, std::size_t N> vector<M> matrix<M, N>::getVariableByIndex(std::size_t index) const{
            if(index >= N){
                throw std::invalid_argument("Index must be strictly less than " + std::to_string(N));
            }
            return _variables[index];
        }

        template <std::size_t M, std::size_t N> vector<M> matrix<M, N>::getColumnByIndex(std::size_t index) const{
            return this->getVariableByIndex(index);
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
        }

        template <std::size_t M, std::size_t N> void matrix<M, N>::convert(bool columnMajorOrderTarget){
            
            
            
            std::cout << "Not implemented.\n";
            return;
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

        template <std::size_t M, std::size_t N> bool matrix<M, N>::isZeroMatrix() const{
            if(this->_columnMajorOrder){
                for(std::size_t i{0}; i < N; ++i){
                    if(this->_variables[i].isZeroVector() == false)return false;
                }
            }else{
                for(std::size_t i{0}; i < M; ++i){
                    if(this->_systems[i].isZeroVector() == false) return false;
                }
            }
            return true;
        }

        template <std::size_t M, std::size_t N> double matrix<M,N>::leadingEntry(std::size_t row) const{
            if(row >= M){
                throw std::invalid_argument("Row specified is outside the matrix's bounds.");
            }
            return _systems[row].leadingEntry();
        }

        template <std::size_t M, std::size_t N> bool matrix<M,N>::isIdentityMatrix() const{
            /*  Algorithm:
             *  First check if square matrix. If no, return false. (Square check)
             *  Check if the k,k index of the matrix is 1. If no, return false. (Diagonality check)
             *  Check if density of each row/column = 1. (Off-diagonal check)
             * 
             *  IGNORE THIS
             *  Add the rows going down (if row-major), or add the columsn going right (if column-major)
             *  If the magnitude of the summed vector is not equal to the dimension, we have non-zero entries in the off-diagonal region.
             */
            
            if(M != N)return false;
            
            for(std::size_t i {0}; i < M; ++i){
                if(_columnMajorOrder){
                    //_variables
                    if(_variables[i].getValue(i) != 1)return false;
                }else{
                    //_systems
                    if(_systems[i].getValue(i) != 1)return false;
                }
            }

            for(std::size_t i {0}; i < M; ++i){
                if(_columnMajorOrder){
                    if(_variables[i].getDensity() != 1) return false;
                }else{
                    if(_systems[i].getDensity() != 1) return false;
                }
            }

            return true;



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
            double outputArray[M*N]{0};
            matrix<M,N> output {outputArray, M*N, false};
            for(std::size_t i{0}; i < M; ++i){
                for(std::size_t j{0}; j < N; ++j){
                    output.setEntry(i, j, 
                        this->search(i, j) + other.search(i, j)
                    );
                }
            }
            return output;
        }

        template <std::size_t M, std::size_t N> matrix<M,N> matrix<M, N>::operator-(matrix<M,N>& other) const{
            matrix<M,N> output = other*(-1);
            return this + output;
        }

        template <std::size_t M, std::size_t N> matrix<M, N> matrix<M, N>::operator*(double scalar) const{
            double outputArray[M*N]{0};
            matrix<M,N> output{outputArray, M*N, false};
            for(std::size_t i{0}; i < M; ++i){
                for(std::size_t j{0}; j < N; ++j){
                    output.setEntry(i, j, 
                        this->search(i, j) * scalar
                    );
                }
            }
            return output;
        }
        
        template <std::size_t M, std::size_t N> template <std::size_t K> matrix<M,K> matrix<M,N>::operator*(matrix<N,K>& other) const{
            //algorithm: (naive)
            //this matrix (A) should be accessed in row-major order
            //other matrix (B) should be accessed in column-major order
            //transpose B -> B^T (now KxN matrix)
            std::size_t newCapacity{M*K};
            double newMatrixArray[newCapacity];
            matrix<K,N> B_transpose { ~other };

            //now we can access B in row-major order
            for(std::size_t i{0}; i < M; ++i){
                for(std::size_t j{0}; j < K; ++j){
                    vector<N> system {B_transpose.getSystemByIndex(j)};
                    newMatrixArray[i * K + j]= this->_systems[i] * system;
                }
            }

            return matrix<M,K>{newMatrixArray, newCapacity, false};
            
	    }

        template <std::size_t M, std::size_t N> matrix<N,M> matrix<M,N>::operator~() const{
            matrix<N,M> transpose {this->_entries, N*M, !this->_columnMajorOrder};
            return transpose;
        }

        template <std::size_t M, std::size_t N> bool matrix<M, N>::operator==(matrix<M, N>& other) const{
            return this->isEqual(other);
        }

        template <std::size_t M, std::size_t N> bool matrix<M, N>::operator!=(matrix<M, N>& other) const{
            return this->isEqual(other) == false;
        }

        template <std::size_t M, std::size_t N> bool matrix<M, N>::isEqual(matrix<M,N>& other, double tolerance) const{
            for(std::size_t i{0}; i < M; ++i){
                for(std::size_t j{0}; j < N; ++j){
                    if(nm::isEqual(this->_systems[i].getValue(j), other.getSystemByIndex(i).getValue(j), tolerance) == false)return false;
                }
            }
            return true;
        }

        template <std::size_t M, std::size_t N> void matrix<M, N>::swapRows(std::size_t row1Index, std::size_t row2Index){

            if(row1Index >= M){
                throw std::invalid_argument("Index row1Index must be strictly less than " + std::to_string(M));
            }
            if(row2Index >= M){
                throw std::invalid_argument("Index row2Index must be strictly less than " + std::to_string(M));
            }
            vector<N> system {this->_systems[row1Index]}; //row1 is saved
            this->_systems[row1Index] = this->_systems[row2Index];
            this->_systems[row2Index] = system;
        }

        template <std::size_t M,std::size_t N> void matrix<M, N>::multiplyRowByScalar(std::size_t rowIndex, double scalar){
            if(rowIndex >= M){
                throw std::invalid_argument("Index rowIndex must be strictly less than " + std::to_string(M));
            }

            this->_systems[rowIndex] = this->_systems[rowIndex] * scalar;
        }

        template <std::size_t M,std::size_t N> void matrix<M, N>::addRow(std::size_t targetRowIndex, std::size_t sourceRowIndex, double scalar){
            if(targetRowIndex >= M){
                throw std::invalid_argument("Index targetRowIndex must be strictly less than " + std::to_string(M));
            }
            if(sourceRowIndex >= M){
                throw std::invalid_argument("Index sourceRowIndex must be strictly less than " + std::to_string(M));
            }

            this->_systems[targetRowIndex] = this->_systems[sourceRowIndex] * scalar + this->_systems[targetRowIndex];
        }

        template <std::size_t M, std::size_t N> void matrix<M,N>::gaussJordanElimination(vector<M> b) const {

            const std::size_t capacity { M*(N+1)};
            double augmentedSystemArray [capacity] = {0};

            //for(std::size_t i{0}; i < M; ++i)std::cout << "index " << i << " density: " << _systems[i].getDensity() << "\n";
            //for(std::size_t j{0}; j < N; ++j)std::cout << "index " << j << " density: " << _variables[j].getDensity() << "\n";



            //populate augmented system (populate b vector first)
            for(std::size_t a{1}; a <= M; ++a){
                augmentedSystemArray[a*(N+1)-1] = b.getValue(a-1);
            }


            for(std::size_t a{0}; a < capacity; ++a){

                if((a+1) % (N+1) == 0)continue; //skip the entries with the b vector entries
                std::size_t row {(a+1) / (N+1)};
                augmentedSystemArray[a] = this->_systems[row].getValue((a) % (N+1));
                std::cout << "index " << a << ": " << augmentedSystemArray[a] << ", row: " << row << "\n" ;
            }

            //don't remove the for loop, the c++ compiler is smart enough to optimize it out :)
            for(std::size_t a{0}; a<capacity; ++a){
                std::cout << "index " << a << ": " << augmentedSystemArray[a] << "\n";
            }
            matrix<M, N+1> augmentedSystem {augmentedSystemArray, capacity, false }; //row-major by default

            for(std::size_t i{0}; i < N+1; ++i){
                std::cout << "variable " << i << " density: " << _variables[i].getDensity() << "\n";
            }
            for(std::size_t j{0}; j < M; ++j){
                std::cout << "system " << j << " density: " << _systems[j].getDensity() << "\n";
            }


            augmentedSystem.print();

            for(std::size_t j {0}; j < N+1; ++j){
                //if this column vector has 0 leading entries, skip to the next one
                std::cout << "j: " << j << "\n";
                augmentedSystem.getVariableByIndex(j).print();
                if(augmentedSystem.getVariableByIndex(j).isZeroVector())continue;

                //get the max leading entry
                std::size_t maxIndex{0};
                double max{augmentedSystem.getSystemByIndex(0).getValue(j)};
                for(std::size_t l {0}; l < M; ++l){
                    double candidate{augmentedSystem.getSystemByIndex(l).getValue(j)};
                    if(candidate > max) {
                        max = candidate;
                        maxIndex = l;
                    }
                }
                std::cout << "max index: " << maxIndex << ", max: " << max << "\n";

                for(std::size_t i {0}; i < M; ++i){
                    std::cout << "max index: " << maxIndex << ", i: " << i << "\n"; 
                    //swap the rows
                    augmentedSystem.swapRows(maxIndex,i);
                    
                    for(std::size_t k{0}; k < M; ++k){
                        if(k == i)continue;
                        
                        double a_kj = augmentedSystem.getSystemByIndex(k).getValue(j);
                        double a_ij = augmentedSystem.getSystemByIndex(i).getValue(j);
                        std::cout << "a_kj: " << a_kj << ", a_ij: " << a_ij << "\n";

                        augmentedSystem.print();
                        augmentedSystem.addRow(k, i, (a_kj/a_ij)*(-1));
                        std::cout << "j: " << j << ", i: " << i << ", k: " << k << "\n"; 
                        augmentedSystem.print();
                    }
                    
                }

            }

            //matrix<M, N+1> augmentedSystem {};

            augmentedSystem.print();

        }



    }
}


bool basicTests(){
    std::size_t score {0};
    std::size_t numberOfTests {0};

    std::cout << "\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nSTART BASIC TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";
    
    std::cout << "Testing NaN detection against C++ implementation: ";
    numberOfTests++;
    double nan_number = 0.0/0.0;
    //std::cout << "C++ Built-in: " << std::isnan(num) << ", NM: " << !nm::checkIfValid(num) << "\n";
    if(std::isnan(nan_number) == !nm::checkIfValid(nan_number)){
        //std::cout << "Both functions identified the number as not-a-number\n";
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    
    std::cout << "Testing +infinity detection: ";
    numberOfTests++;
    double positive_infinity {1.0/0.0};
    if(!nm::checkIfValid(positive_infinity)){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "Testing -infinity detection: ";
    numberOfTests++;
    double negative_infinity {-1.0/0.0};
    if(!nm::checkIfValid(negative_infinity)){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing isEqual() (for equal numbers) with global tolerance " << nm::globalTolerance << ": ";
    numberOfTests++;
    double isEqualTestNumber1 {0.0002};
    double isEqualTestNumber2 {0.0005};
    std::cout << "(" << isEqualTestNumber1 << "," << isEqualTestNumber2 << ") ";
    if(nm::isEqual(isEqualTestNumber1, isEqualTestNumber2, nm::globalTolerance)){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing isEqual() (for non-equal numbers) with global tolerance " << nm::globalTolerance << ": ";
    numberOfTests++;
    double isEqualTestNumber3 {0.0306};
    double isEqualTestNumber4 {0.0005};
    std::cout << "(" << isEqualTestNumber3 << "," << isEqualTestNumber4 << ") ";
    if(nm::isEqual(isEqualTestNumber3, isEqualTestNumber4, nm::globalTolerance) == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing isZero() (for equal numbers) with global tolerance " << nm::globalTolerance << ": ";
    numberOfTests++;
    double isZeroTestNumber1{0.00004};
    std::cout << "(" << isZeroTestNumber1 << ") ";
    if(nm::isZero(isZeroTestNumber1, nm::globalTolerance)){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing isZero() (for non-equal numbers) with global tolerance " << nm::globalTolerance << ": ";
    numberOfTests++;
    double isZeroTestNumber2{0.04045};
    std::cout << "(" << isZeroTestNumber2 << ") ";
    if(nm::isZero(isZeroTestNumber2, nm::globalTolerance) == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nEND BASIC TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";


    return numberOfTests == score;

}

bool vectorTests(){
    using namespace nm::linalg;
    std::size_t score{0};
    std::size_t numberOfTests = 0;

    std::cout << "\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nSTART VECTOR TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";

    //Usual vector array
    double threeDimensionalVectorArray1 [3] = {1,2.5,6};

    //Vector array with density of 2, with zero in the 0th position
    double threeDimensionalVectorArray2 [3] = {0 , 4, 3.5};

    double negativeThreeDimensionalVectorArray1 [3] = {0, 0, 0};
    for(std::size_t i {0}; i < 3; ++i){
        negativeThreeDimensionalVectorArray1[i] = (-1)*threeDimensionalVectorArray1[i];
    }

    double threeDimensionalZeroVectorArray [3] = {0,0,0};


    std::cout << "Testing correct array size: ";
    numberOfTests++;
    try{
        vector<3> threeDimensionalVector{threeDimensionalVectorArray1, 3};
        std::cout << "Success!";
        score++;
    }catch(std::invalid_argument message){
        std::cout << "Failure\n";
        std::cout << message.what();
    }
    std::cout << "\n";



    std::cout << "Testing incorrect declared array size: ";
    numberOfTests++;
    try{
        vector<3> threeDimensionalVector{threeDimensionalVectorArray1, 4};
        std::cout << "Failure\n";
    }catch(std::invalid_argument message){
        std::cout << "Success!\n";
        std::cout << message.what();
        score++;
    }
    std::cout << "\n";



    std::cout << "Testing zero vector detection: ";
    numberOfTests++;
    vector<3> zeroVector {threeDimensionalZeroVectorArray, 3};
    if(zeroVector.isZeroVector()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "Testing vector equality (for equal vectors): ";
    numberOfTests++;
    double equalTestVectorArray1[3] = {2,3,4};
    vector<3> equalTestVector1 {equalTestVectorArray1, 3};
    if(equalTestVector1 == equalTestVector1){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing vector equality (for non-equal vectors): ";
    numberOfTests++;
    double equalTestVectorArray2[3] = {4,3,2};
    vector<3> equalTestVector2 {equalTestVectorArray2, 3};
    if(equalTestVector2 == equalTestVector1){
        std::cout << "Failure\n";
    }else{
        std::cout << "Success!\n";
        score++;
    }
    
    
    std::cout << "Testing vector inequality (for non-equal vectors): ";
    numberOfTests++;
    if(equalTestVector2 != equalTestVector1){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing vector inequality (for equal vectors): ";
    numberOfTests++;
    if(equalTestVector1 != equalTestVector1){
        std::cout << "Failure\n";
    }else{
        std::cout << "Success!\n";
       score++;
    }
    

    std::cout << "Testing scalar multiplication by 0 (uses zero vector detection): ";
    numberOfTests++;
    vector<3> scalarMultiplicationTestVector {threeDimensionalVectorArray1, 3};
    if((scalarMultiplicationTestVector*0).isZeroVector()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing vector subtraction 1 (additive inverse property, uses zero vector detection): ";
    numberOfTests++;
    vector<3> additiveInverseTestVector1 {threeDimensionalVectorArray1, 3};
    if((additiveInverseTestVector1 - additiveInverseTestVector1).isZeroVector()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing vector subtraction 2: ( 0 - a = -a )";
    numberOfTests++;
    vector<3> additiveInverseTestVector2 {negativeThreeDimensionalVectorArray1, 3};
    if((zeroVector - additiveInverseTestVector1) == additiveInverseTestVector2){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing vector addition. (Depends on scalar multiplication by 2): ";
    numberOfTests++;
    vector<3> additiveTestVector {threeDimensionalVectorArray1, 3};
    vector<3> additiveTestVectorScaledBy2 = additiveTestVector*2;
    if((additiveTestVector + additiveTestVector) == additiveTestVectorScaledBy2){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing dot product: ";
    numberOfTests++;
    double dotProductTestArray1[3] = {0,0,1};
    double dotProductTestArray2[3] = {1,0,0};
    vector<3> dotProductTestVector1 {dotProductTestArray1,3};
    vector<3> dotProductTestVector2 {dotProductTestArray2, 3};
    if((dotProductTestVector1 * dotProductTestVector2) == 0){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "Testing norm (depends on dot product): ";
    numberOfTests++;
    vector<3> normTestVector {threeDimensionalVectorArray1, 3};
    if(normTestVector.norm2()*normTestVector.norm2() == normTestVector*normTestVector){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "Leading entry test: ";
    numberOfTests++;
    vector<3> leadingEntryValueTestVector { threeDimensionalVectorArray2, 3};
    if(leadingEntryValueTestVector.leadingEntry() == 4){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Leading entry index test: ";
    numberOfTests++;
    vector<3> leadingEntryValueIndexVector { threeDimensionalVectorArray2, 3};
    if(leadingEntryValueIndexVector.leadingEntryIndex() == 1){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing set value: ";
    numberOfTests++;
    double setValueVectorArray1[3] = {0,2,0};
    double setValueVectorArray2[3] = {0,0,0};
    vector<3> setValueVector1 {setValueVectorArray1, 3};
    vector<3> setValueVector2 {setValueVectorArray2, 3};
    setValueVector2.setValue(2, 1);
    if ((setValueVector1 - setValueVector2).isZeroVector()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }
    std::cout << "Testing set value density: ";
    numberOfTests++;
    if(setValueVector1.getDensity() == setValueVector2.getDensity() && setValueVector1.getDensity() == 1){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing set values: ";
    numberOfTests++;

    double getsetValuesVectorArray1[3] = {4,0,1};
    double getsetValuesVectorArray2[3] = {1,2,3};
    double getsetValuesVectorArray3[3] = {4,1,2};
    //this vector gets changed
    vector<3> getsetValuesVector1 {getsetValuesVectorArray1, 3};
    vector<3> getsetValuesVector2 {getsetValuesVectorArray2, 3};
    vector<3> getsetValuesVector3 {getsetValuesVectorArray3, 3};
    std::size_t sv_offset{1};
    getsetValuesVector1.setValues(getsetValuesVectorArray2, sv_offset, 3);

    if(getsetValuesVector1.getValue(sv_offset) == 1 && getsetValuesVector1 == getsetValuesVector3){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }
    //setValuesVector1.print();
    //setValuesVector2.print();
    //setValuesVector3.print();
    //std::cout << setValuesVector1.getDensity() << ", " << setValuesVector2.getDensity() << ", " << setValuesVector3.getDensity() << "\n"; 


    std::cout << "Testing set values density: ";
    numberOfTests++;
    if(getsetValuesVector1.getDensity() == getsetValuesVector3.getDensity()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing get value: ";
    numberOfTests++;
    std::size_t gv_index {0};
    //std::cout << setValuesVector1.getValue(gv_index) << " " << setValuesVectorArray1[gv_index] << "\n";
    if(getsetValuesVector1.getValue(gv_index) == getsetValuesVectorArray1[gv_index]){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing get values: ";
    numberOfTests++;
    std::size_t getValueCounter {0};
    for(std::size_t i {0}; i < 3; ++i){
        if(getsetValuesVector2.getValues()[i] == getsetValuesVectorArray2[i]) getValueCounter++;
    }
    if(getValueCounter == 3){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing get density: ";
    numberOfTests++;
    double densityArray [3] = {4,0,2};
    vector<3> getDensityTest {densityArray, 3};
    if(getDensityTest.getDensity() == 2){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing get density normalized: ";
    numberOfTests++;
    vector<3> getDensityNormalizedTest {densityArray, 3};
    if(nm::isEqual(getDensityNormalizedTest.getDensityNormalized(),(double)2/(double)3, nm::globalTolerance)){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing get dimension: ";
    numberOfTests++;
    double dimensionTestVectorArray[6] = {0,0,0,0,0,4};
    vector<6> dimensionTestVector{dimensionTestVectorArray, 6};
    if(dimensionTestVector.dimension() == 6){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nEND VECTOR TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";

    return score == numberOfTests;

}

bool matrixTests(){

    using namespace nm::linalg;
    std::size_t score {0};
    std::size_t numberOfTests {0};

    std::cout << "\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nSTART MATRIX TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";
    
    // 1 x N matrix test
    // M x 1 matrix test
    // 2 x 2 matrix, CM vs RM test
    // 3 x 3 matrix, CM vs RM test
    // 1 x N matrix density test
    // M x 1 matrix density test
    // 2 x 2 matrix density test, CM vs RM
    // 3 x 3 matrix density test, CM vs RM
    // 1 x 1 search test
    // 2 x 2 search test
    // 2 x 2 setEntry test
    // 3 x 3 setEntry test
    // 1 x 1 convert test
    // 2 x 2 convert test
    // 3 x 3 convert test
    // 3 x 2 convert test
    // 2 x 3 convert test
    // 4 x 2 convert test
    // checkIfValidMatrix test
    // 3 x 3 isRowEmpty
    // 3 x 3 isColumnEmpty
    // 3 x 3 leadingEntryTest (use REF)
    // 2 x 2 isIdentityMatrix test
    // 3 x 3 isIdentityMatrix test
    // 3 x 3 scalar multiplication test
    // 3 x 3 matrix addition test 1 (additive inverse)
    // 3 x 3 matrix addition test 2 ()
    // 3 x 4 matrix multiplied by 4 x 2 matrix test
    // 3 x 3 matrix transpose test
    // 4 x 3 -> 3 x 4 matrix transpose test
    // 3 x 3 matrix equality test
    // 3 x 3 swapRows test
    // 3 x 3 multiplyRowByScalar test
    // 3 x 3 addRow test

    std::cout << "Testing 1 X N matrix size: ";
    numberOfTests++;
    double m1x5Array [5] = {1,2,3,4,5};
    matrix<1, 5> m1x5 {m1x5Array, 5, false};
    if(m1x5.getSystemByIndex(0).dimension() == 5 && m1x5.getVariableByIndex(0).dimension() == 1){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 1 x N matrix equality vs vector<N>: ";
    numberOfTests++;
    vector<5> v5 {m1x5Array, 5};
    if(m1x5.getSystemByIndex(0) == v5){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing M x 1 matrix size: ";
    numberOfTests++;
    double m6x1Array [6] = {3,5,4,6,5,7};
    matrix<6, 1> m6x1 {m6x1Array, 6, true};
    if(m6x1.getSystemByIndex(0).dimension() == 1 && m6x1.getVariableByIndex(0).dimension() == 6){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing M x 1 matrix equality vs vector<M>: ";
    numberOfTests++;
    vector<6> v6 {m6x1Array, 6};
    if(m6x1.getVariableByIndex(0) == v6){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 2 x 2 matrix equality between column-major and row-major forms: ";
    numberOfTests++;
    double m2x2equalityTestArray [4] = {3,6,4,8};
    matrix<2, 2> m2x2equalityTestMatrixCM {m2x2equalityTestArray, 4, true};
    matrix<2, 2> m2x2equalityTestMatrixRM {m2x2equalityTestArray, 4, false};
    //std::cout << "\n";
    //m2x2equalityTestMatrixCM.print();
    //m2x2equalityTestMatrixRM.print();

    vector<2> m2x2equalityTestMatrixCMColumn0 = m2x2equalityTestMatrixCM.getVariableByIndex(0);
    vector<2> m2x2equalityTestMatrixCMColumn1 = m2x2equalityTestMatrixCM.getVariableByIndex(1);
    vector<2> m2x2equalityTestMatrixCMRow0 = m2x2equalityTestMatrixCM.getSystemByIndex(0);
    vector<2> m2x2equalityTestMatrixCMRow1 = m2x2equalityTestMatrixCM.getSystemByIndex(1);

    vector<2> m2x2equalityTestMatrixRMColumn0 = m2x2equalityTestMatrixRM.getVariableByIndex(0);
    vector<2> m2x2equalityTestMatrixRMColumn1 = m2x2equalityTestMatrixRM.getVariableByIndex(1);
    vector<2> m2x2equalityTestMatrixRMRow0 = m2x2equalityTestMatrixRM.getSystemByIndex(0);
    vector<2> m2x2equalityTestMatrixRMRow1 = m2x2equalityTestMatrixRM.getSystemByIndex(1);

    if(
        m2x2equalityTestMatrixCMColumn0 == m2x2equalityTestMatrixRMRow0 &&
        m2x2equalityTestMatrixCMColumn1 == m2x2equalityTestMatrixRMRow1 /*&&
        m2x2equalityTestMatrixRMColumn0 == m2x2equalityTestMatrixCMRow0 &&
        m2x2equalityTestMatrixRMColumn1 == m2x2equalityTestMatrixCMRow1*/
        // we DO NOT want RM columns or CM rows, those are computed by convert(); separate test
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix equality between column-major and row-major forms: ";
    numberOfTests++;
    double m3x3equalityTestArray [9] = {2,4,3,5,4,6,5,7,6};
    matrix<3, 3> m3x3equalityTestMatrixCM {m3x3equalityTestArray, 9, true};
    matrix<3, 3> m3x3equalityTestMatrixRM {m3x3equalityTestArray, 9, false};

    //std::cout << "\n";
    //m3x3equalityTestMatrixCM.print();
    //m3x3equalityTestMatrixRM.print();

    vector<3> m3x3equalityTestMatrixCMColumn0 = m3x3equalityTestMatrixCM.getVariableByIndex(0);
    vector<3> m3x3equalityTestMatrixCMColumn1 = m3x3equalityTestMatrixCM.getVariableByIndex(1);
    vector<3> m3x3equalityTestMatrixCMColumn2 = m3x3equalityTestMatrixCM.getVariableByIndex(2);
    vector<3> m3x3equalityTestMatrixCMRow0 = m3x3equalityTestMatrixCM.getSystemByIndex(0);
    vector<3> m3x3equalityTestMatrixCMRow1 = m3x3equalityTestMatrixCM.getSystemByIndex(1);
    vector<3> m3x3equalityTestMatrixCMRow2 = m3x3equalityTestMatrixCM.getSystemByIndex(2);

    vector<3> m3x3equalityTestMatrixRMColumn0 = m3x3equalityTestMatrixRM.getVariableByIndex(0);
    vector<3> m3x3equalityTestMatrixRMColumn1 = m3x3equalityTestMatrixRM.getVariableByIndex(1);
    vector<3> m3x3equalityTestMatrixRMColumn2 = m3x3equalityTestMatrixRM.getVariableByIndex(2);
    vector<3> m3x3equalityTestMatrixRMRow0 = m3x3equalityTestMatrixRM.getSystemByIndex(0);
    vector<3> m3x3equalityTestMatrixRMRow1 = m3x3equalityTestMatrixRM.getSystemByIndex(1);
    vector<3> m3x3equalityTestMatrixRMRow2 = m3x3equalityTestMatrixRM.getSystemByIndex(2);

    if(
        m3x3equalityTestMatrixCMColumn0 == m3x3equalityTestMatrixRMRow0 &&
        m3x3equalityTestMatrixCMColumn1 == m3x3equalityTestMatrixRMRow1 &&
        m3x3equalityTestMatrixCMColumn2 == m3x3equalityTestMatrixRMRow2 /*&&
        m3x3equalityTestMatrixCMRow0 == m3x3equalityTestMatrixRMColumn0 &&
        m3x3equalityTestMatrixCMRow1 == m3x3equalityTestMatrixRMColumn1 &&
        m3x3equalityTestMatrixCMRow2 == m3x3equalityTestMatrixRMColumn2*/
        // we DO NOT want RM columns or CM rows, those are computed by convert(); separate test
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 1 x N matrix density vs vector<N>: ";
    numberOfTests++;
    double m1x5densityTestArray [5] = {0,3,2,0,1};
    vector<5> m1x5densityTestVector {m1x5densityTestArray, 5};
    matrix<1, 5> m1x5densityTestMatrix {m1x5densityTestArray, 5, false};
    if(m1x5densityTestMatrix.getSystemByIndex(0).getDensity() == m1x5densityTestVector.getDensity() && m1x5densityTestMatrix.getSystemByIndex(0).getDensity() == 3){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing M x 1 matrix density vs vector<M>: ";
    numberOfTests++;
    double m6x1densityTestArray [6] = {0,0,1,0,3.5,8};
    vector<6> m6x1densityTestVector {m6x1densityTestArray, 6};
    matrix<6, 1> m6x1densityTestMatrix {m6x1densityTestArray, 6, true};
    if(m6x1densityTestMatrix.getVariableByIndex(0).getDensity() == m6x1densityTestVector.getDensity() && m6x1densityTestMatrix.getVariableByIndex(0).getDensity() == 3){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 2 x 2 matrix density equality between row-major and column-major forms: ";
    numberOfTests++;
    double m2x2densityTestArray [4] = {0,2,3,5};
    matrix<2, 2> m2x2densityTestMatrixCM = {m2x2densityTestArray, 4, true};
    matrix<2, 2> m2x2densityTestMatrixRM = {m2x2densityTestArray, 4, false};
    if(
        m2x2densityTestMatrixCM.getSystemByIndex(0).getDensity() == m2x2densityTestMatrixRM.getVariableByIndex(0).getDensity() &&
        m2x2densityTestMatrixCM.getSystemByIndex(0).getDensity() == 1 &&
        m2x2densityTestMatrixCM.getSystemByIndex(1).getDensity() == m2x2densityTestMatrixRM.getVariableByIndex(1).getDensity() &&
        m2x2densityTestMatrixCM.getSystemByIndex(1).getDensity() == 2
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix density equality between row-major and column-major forms: ";
    numberOfTests++;
    double m3x3densityTestArray [9] = {1,2,3,0,4,5,0,0,6};
    matrix<3, 3> m3x3densityTestMatrixCM = {m3x3densityTestArray, 9, true};
    matrix<3, 3> m3x3densityTestMatrixRM = {m3x3densityTestArray, 9, false};
    if(
        m3x3densityTestMatrixCM.getSystemByIndex(0).getDensity() == m3x3densityTestMatrixRM.getVariableByIndex(0).getDensity() &&
        m3x3densityTestMatrixCM.getSystemByIndex(0).getDensity() == 1 &&
        m3x3densityTestMatrixCM.getSystemByIndex(1).getDensity() == m3x3densityTestMatrixRM.getVariableByIndex(1).getDensity() &&
        m3x3densityTestMatrixCM.getSystemByIndex(1).getDensity() == 2 &&
        m3x3densityTestMatrixCM.getSystemByIndex(2).getDensity() == m3x3densityTestMatrixRM.getVariableByIndex(2).getDensity() &&
        m3x3densityTestMatrixCM.getSystemByIndex(2).getDensity() == 3
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 1 x 1 matrix search(): ";
    numberOfTests++;
    double m1x1searchTestArray [1] = {4};
    matrix<1, 1> m1x1searchTestMatrix {m1x1searchTestArray, 1, false};
    if(m1x1searchTestMatrix.search(0, 0) == m1x1searchTestArray[0]){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 2 x 2 matrix search() (input array vs both row-major and column-major forms): ";
    numberOfTests++;
    double m2x2searchTestArray [4] = {0,2,1,0};
    matrix<2, 2> m2x2searchTestMatrixCM {m2x2searchTestArray, 4, true};
    matrix<2, 2> m2x2searchTestMatrixRM {m2x2searchTestArray, 4, false};
    if(
        m2x2searchTestMatrixRM.search(0, 0) == m2x2searchTestArray[0] &&
        m2x2searchTestMatrixRM.search(0, 1) == m2x2searchTestArray[1] &&
        m2x2searchTestMatrixRM.search(1, 0) == m2x2searchTestArray[2] &&
        m2x2searchTestMatrixRM.search(1, 1) == m2x2searchTestArray[3] &&
        m2x2searchTestMatrixCM.search(0, 0) == m2x2searchTestArray[0] &&
        m2x2searchTestMatrixCM.search(0, 1) == m2x2searchTestArray[2] &&
        m2x2searchTestMatrixCM.search(1, 0) == m2x2searchTestArray[1]
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 2 x 2 matrix setEntry(): ";
    numberOfTests++;
    double m2x2setEntryTestArray [4] = {0,0,2,1};
    matrix<2, 2> m2x2setEntryTestMatrix {m2x2setEntryTestArray, 4, false}; //RM vs CM is arbitrary
    m2x2setEntryTestMatrix.setEntry(0, 1, 5.5);
    if(m2x2setEntryTestMatrix.search(0, 1) == 5.5){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix setEntry(): ";
    numberOfTests++;
    double m3x3setEntryTestArray [9] = {0,0,2,1,0,0,0,0,0};
    matrix<3, 3> m3x3setEntryTestMatrix {m3x3setEntryTestArray, 9, false}; //RM vs CM is arbitrary
    m2x2setEntryTestMatrix.setEntry(0, 1, 5.5);
    if(m2x2setEntryTestMatrix.search(0, 1) == 5.5){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    //RM <-> systems
    //CM <-> variables

    std::cout << "Testing 1 x 1 convert() equality between row-major and column-major forms: ";
    numberOfTests++;
    double m1x1convertTestArray [1] = {3};
    matrix<1, 1> m1x1convertTestMatrixRM {m1x1convertTestArray, 1, false};
    matrix<1, 1> m1x1convertTestMatrixCM {m1x1convertTestArray, 1, true};
    vector<1> m1x1convertTestMatrixCMRow0 {m1x1convertTestMatrixCM.getSystemByIndex(0)};
    vector<1> m1x1convertTestMatrixRMColumn0 {m1x1convertTestMatrixRM.getVariableByIndex(0)};
    if(m1x1convertTestMatrixRMColumn0 == m1x1convertTestMatrixCMRow0){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 2 x 2 convert() equality between row-major and column-major forms: ";
    numberOfTests++;
    double m2x2convertTestArray [4] = {3,4,2,7};
    matrix<2, 2> m2x2convertTestMatrixRM {m2x2convertTestArray, 4, false};
    matrix<2, 2> m2x2convertTestMatrixCM {m2x2convertTestArray, 4, true};
    //std::cout << '\n';
    //m2x2convertTestMatrixRM.print();m2x2convertTestMatrixCM.print();

    vector<2> m2x2convertTestMatrixCMRow0 {m2x2convertTestMatrixCM.getSystemByIndex(0)};
    vector<2> m2x2convertTestMatrixCMRow1 {m2x2convertTestMatrixCM.getSystemByIndex(1)};
    vector<2> m2x2convertTestMatrixRMColumn0 {m2x2convertTestMatrixRM.getVariableByIndex(0)};
    vector<2> m2x2convertTestMatrixRMColumn1 {m2x2convertTestMatrixRM.getVariableByIndex(1)};

    //m2x2convertTestMatrixCMRow0.print(); m2x2convertTestMatrixCMRow1.print();
    //m2x2convertTestMatrixRMColumn0.print(); m2x2convertTestMatrixRMColumn1.print();

    if(
        m2x2convertTestMatrixCMRow0 == m2x2convertTestMatrixRMColumn0 &&
        m2x2convertTestMatrixCMRow1 == m2x2convertTestMatrixRMColumn1
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 convert() equality between row-major and column-major forms: ";
    numberOfTests++;
    double m3x3convertTestArray [9] = {5.5,3.5,7.8,9.2,4.2,5,1,-5,-7.3};
    matrix<3, 3> m3x3convertTestMatrixRM {m3x3convertTestArray, 9, false};
    matrix<3, 3> m3x3convertTestMatrixCM {m3x3convertTestArray, 9, true};
    //std::cout << '\n';
    //m3x3convertTestMatrixRM.print();m3x3convertTestMatrixCM.print();
    vector<3> m3x3convertTestMatrixCMRow0 {m3x3convertTestMatrixCM.getSystemByIndex(0)};
    vector<3> m3x3convertTestMatrixCMRow1 {m3x3convertTestMatrixCM.getSystemByIndex(1)};
    vector<3> m3x3convertTestMatrixCMRow2 {m3x3convertTestMatrixCM.getSystemByIndex(2)};
    vector<3> m3x3convertTestMatrixRMColumn0 {m3x3convertTestMatrixRM.getVariableByIndex(0)};
    vector<3> m3x3convertTestMatrixRMColumn1 {m3x3convertTestMatrixRM.getVariableByIndex(1)};
    vector<3> m3x3convertTestMatrixRMColumn2 {m3x3convertTestMatrixRM.getVariableByIndex(2)};

    //m3x3convertTestMatrixCMRow0.print(); m3x3convertTestMatrixCMRow1.print(); m3x3convertTestMatrixCMRow0.print();
    //m3x3convertTestMatrixRMColumn0.print(); m3x3convertTestMatrixRMColumn1.print(); m3x3convertTestMatrixRMColumn2.print();

    if(
        m3x3convertTestMatrixCMRow0 == m3x3convertTestMatrixRMColumn0 &&
        m3x3convertTestMatrixCMRow1 == m3x3convertTestMatrixRMColumn1 &&
        m3x3convertTestMatrixCMRow2 == m3x3convertTestMatrixRMColumn2
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    //for 4x2:
    //add vector padding and trimming


    double m3x3isRowColumnEmptyTestArray [9] = {0,2.2,3,0,0,4.7,0,0,0};
    matrix<3, 3> m3x3isRowColumnEmptyTestMatrix {m3x3isRowColumnEmptyTestArray, 9, false};
    

    std::cout << "Testing 3 x 3 isRowEmpty(): ";
    numberOfTests++;
    if(
        m3x3isRowColumnEmptyTestMatrix.isRowEmpty(0) == false &&
        m3x3isRowColumnEmptyTestMatrix.isRowEmpty(1) == false &&
        m3x3isRowColumnEmptyTestMatrix.isRowEmpty(2)
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 isColumnEmpty(): ";
    numberOfTests++;
    if(
        m3x3isRowColumnEmptyTestMatrix.isColumnEmpty(0) &&
        m3x3isRowColumnEmptyTestMatrix.isColumnEmpty(1) == false &&
        m3x3isRowColumnEmptyTestMatrix.isColumnEmpty(2) == false
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 leadingEntry() (REF): ";
    numberOfTests++;
    double m3x3leadingEntryREFTestArray [9] = {1,2,3,0,4,5,0,0,6};
    matrix<3, 3> m3x3leadingEntryREFTestMatrix {m3x3leadingEntryREFTestArray, 9, false};
    if(
        m3x3leadingEntryREFTestMatrix.leadingEntry(0) == 1 &&
        m3x3leadingEntryREFTestMatrix.leadingEntry(1) == 4 &&
        m3x3leadingEntryREFTestMatrix.leadingEntry(2) == 6
    ){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 2 x 2 isIdentityMatrix() (given non-identity matrix): ";
    numberOfTests++;
    double m2x2notIdentityTestArray [4] = {3,4.5,-9,1};
    matrix<2, 2> m2x2notIdentityTestMatrix {m2x2notIdentityTestArray, 4, false};
    if(m2x2notIdentityTestMatrix.isIdentityMatrix() == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 2 x 2 isIdentityMatrix() (given identity matrix): ";
    numberOfTests++;
    double m2x2identityTestArray [4] = {1,0,0,1};
    matrix<2, 2> m2x2identityTestMatrix {m2x2identityTestArray, 4, false};
    if(m2x2identityTestMatrix.isIdentityMatrix()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }



    std::cout << "Testing 3 x 3 isIdentityMatrix() (given non-identity matrix): ";
    numberOfTests++;
    double m3x3notIdentityTestArray [9] = {3,4.5,-9,1,7.6, -10.2, 3.2, 1.3,0.4};
    matrix<3, 3> m3x3notIdentityTestMatrix {m3x3notIdentityTestArray, 9, false};
    if(m3x3notIdentityTestMatrix.isIdentityMatrix() == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 3 x 3 isIdentityMatrix() (given identity matrix): ";
    numberOfTests++;
    double m3x3identityTestArray [9] = {1,0,0,0,1,0,0,0,1};
    matrix<3, 3> m3x3identityTestMatrix {m3x3identityTestArray, 9, false};
    if(m3x3identityTestMatrix.isIdentityMatrix()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 isZeroVector() (given non-zero matrix): ";
    numberOfTests++;
    double m3x3zeroVectorTestArray1 [9] = {0,0,0,0,0,0,0,0,1};
    matrix<3, 3> m3x3zeroVectorTestMatrix1 {m3x3zeroVectorTestArray1, 9, false};
    if(m3x3zeroVectorTestMatrix1.isZeroMatrix() == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 isZeroVector() (given zero matrix): ";
    numberOfTests++;
    double m3x3zeroVectorTestArray2 [9] = {0,0,0,0,0,0,0,0,0};
    matrix<3, 3> m3x3zeroVectorTestMatrix2 {m3x3zeroVectorTestArray2, 9, false};
    if(m3x3zeroVectorTestMatrix2.isZeroMatrix()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix scalar multiplication: ";
    numberOfTests++;
    double m3x3scalarMultiplicationTestScalar {3};
    double m3x3scalarMultiplicationTestArray [9] = {0,-3.5,2,9,5.4,-0.2,3,-2.79,0};
    double m3x3scalarMultiplicationTestArrayScaled [9] = {0,0,0,0,0,0,0,0,0};
    for(std::size_t i {0}; i < 9; ++i){
        m3x3scalarMultiplicationTestArrayScaled[i] = m3x3scalarMultiplicationTestScalar * m3x3scalarMultiplicationTestArray[i];
    }
    matrix<3, 3> m3x3scalarMultiplicationTestMatrix {m3x3scalarMultiplicationTestArray, 9, false};
    matrix<3, 3> m3x3scalarMultiplicationTestMatrixScaled {m3x3scalarMultiplicationTestArrayScaled, 9, false};
    if(m3x3scalarMultiplicationTestMatrix * m3x3scalarMultiplicationTestScalar == m3x3scalarMultiplicationTestMatrixScaled){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix addition 1 (additive inverse): ";
    numberOfTests++;
    double m3x3additionTest1Array [9] = {0,-9.7,0.746,4.32,3.5,-5.3,0,-1.2,0};
    double m3x3additionTest1InverseArray [9] = {0,0,0,0,0,0,0,0,0};
    for(std::size_t i{0}; i < 9; ++i){
        m3x3additionTest1InverseArray[i] = (-1)*m3x3additionTest1Array[i];
    }
    matrix<3, 3> m3x3additionTest1Matrix{m3x3additionTest1Array, 9, false};
    matrix<3, 3> m3x3additionTest1InverseMatrix{m3x3additionTest1InverseArray, 9, false};
    if((m3x3additionTest1Matrix + m3x3additionTest1InverseMatrix).isZeroMatrix()){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix addition 2 (depends on scalar multiplication): ";
    numberOfTests++;
    double m3x3additionTest2Array [9] = {1,0.79,-4.367,0,6,-4.3,1.43,2.3,-7};
    matrix<3, 3> m3x3additionTest2Matrix {m3x3additionTest2Array, 9, false};
    matrix<3, 3> m3x3additionTest2MatrixScaled = m3x3additionTest2Matrix * 2;
    if(m3x3additionTest2Matrix + m3x3additionTest2Matrix == m3x3additionTest2MatrixScaled){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 matrix transpose (using convert() and matrix equality operator): ";
    numberOfTests++;
    double m3x3transposeTestArray [9] = {3.4,-9,0,0.793,-5.32,8,1,7.5,0};
    matrix<3, 3> m3x3transposeTestMatrixRM {m3x3transposeTestArray, 9, false};
    matrix<3, 3> m3x3transposeTestMatrixCM {m3x3transposeTestArray, 9, true};
    if(~m3x3transposeTestMatrixRM == m3x3transposeTestMatrixCM){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 full matrix equality (with equal matrices): ";
    numberOfTests++;
    double m3x3fullEqualityTest1Array [9] = {3,-0.79,-7,5.326,3.69,0,-1,1,4.5555};
    matrix<3, 3> m3x3fullEqualityTest1Matrix{m3x3fullEqualityTest1Array, 9, false};
    if(m3x3fullEqualityTest1Matrix == m3x3fullEqualityTest1Matrix){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 full matrix equality (with unequal matrices): ";
    numberOfTests++;
    double m3x3fullEqualityTest2Array [9] = {4,-5,-2,7.953,0,2.445,-8.3364,1,0};
    matrix<3, 3> m3x3fullEqualityTest2Matrix1 {m3x3fullEqualityTest2Array, 9, false};
    matrix<3, 3> m3x3fullEqualityTest2Matrix2 = m3x3fullEqualityTest2Matrix1*(2.76);
    if((m3x3fullEqualityTest2Matrix1 == m3x3fullEqualityTest2Matrix2) == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 full matrix inequality (with equal matrices): ";
    numberOfTests++;
    double m3x3fullInequalityTest1Array [9] = {4,-3,-2.798,0,6,-7.8,4.420,0,-1};
    matrix<3, 3> m3x3fullInequalityTest1Matrix {m3x3fullInequalityTest1Array, 9, false};
    if((m3x3fullInequalityTest1Matrix != m3x3fullInequalityTest1Matrix) == false){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 full matrix inequality (with unequal matrices): ";
    numberOfTests++;
    double m3x3fullInequalityTest2Array [9] = {3,-2,7.4346,6.3,-5.4,0,-3.24,4.66,-1};
    matrix<3, 3> m3x3fullInequalityTest2Matrix1 { m3x3fullInequalityTest2Array, 9, false};
    matrix<3, 3> m3x3fullInequalityTest2Matrix2 = m3x3fullInequalityTest2Matrix1 * (-9.4);
    if(m3x3fullInequalityTest2Matrix1 != m3x3fullInequalityTest2Matrix2){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing matrix multiplication (3 x 4 * 4 x 2) (depends on matrix equality): ";
    numberOfTests++;
    double m3x4matrixMultiplicationTestArray [12] {1,2,3,4,5,6,7,8,9,10,11,12};
    double m4x2matrixMultiplicationTestArray [8] {1,2,3,4,5,6,7,8};
    double m3x2matrixMultiplicationTestArray [6] = {50,60,114,140,178,220}; //verified
    nm::linalg::matrix<3, 4> m3x4matrixMultiplicationTestMatrix {m3x4matrixMultiplicationTestArray, 12, false};
    nm::linalg::matrix<4, 2> m4x2matrixMultiplicationTestMatrix {m4x2matrixMultiplicationTestArray, 8, false};
    nm::linalg::matrix<3, 2> m3x2matrixMultiplicationTestMatrix {m3x2matrixMultiplicationTestArray, 6, false};
    if(m3x4matrixMultiplicationTestMatrix * m4x2matrixMultiplicationTestMatrix == m3x2matrixMultiplicationTestMatrix){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }

    std::cout << "Testing 3 x 3 swapRows(): ";
    numberOfTests++;
    double m3x3swapRowsTestArray[9] = {0,3,2,4.5,-6.2,1,0,-4,-3};
    std::size_t rowSwapIndex1 {2};
    std::size_t rowSwapIndex2 {0};
    double m3x3swapRowsTestArraySwapped[9] = {0,3,2,4.5,-6.2,1,0,-4,-3};
    std::size_t rowSwap1Position {rowSwapIndex1 * 3};
    std::size_t rowSwap2Position {rowSwapIndex2 * 3};
    for(std::size_t i {0}; i < 3; ++i){
        m3x3swapRowsTestArraySwapped[i + rowSwap2Position] = m3x3swapRowsTestArray[i + rowSwap1Position];
        m3x3swapRowsTestArraySwapped[i + rowSwap1Position] = m3x3swapRowsTestArray[i + rowSwap2Position];
    }
    /*
    for(std::size_t j {0}; j < 9; ++j)std::cout << m3x3swapRowsTestArray[j] << " ";
    std::cout << "\n";
    for(std::size_t j {0}; j < 9; ++j)std::cout << m3x3swapRowsTestArraySwapped[j] << " ";
    std::cout << "\n";
    */
    matrix<3, 3> m3x3swapRowsTestMatrix {m3x3swapRowsTestArray, 9, false};
    matrix<3, 3> m3x3swapRowsTestMatrixSwapped {m3x3swapRowsTestArraySwapped, 9, false};
    m3x3swapRowsTestMatrix.swapRows(0, 2);
    if(m3x3swapRowsTestMatrix == m3x3swapRowsTestMatrixSwapped){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 3 x 3 multiplyRowByScalar(): ";
    numberOfTests++;
    double m3x3multiplyRowByScalarTestArray [9] = {5,-3.54,0.3,9,2,0,-1,-4.675,2.3};
    double m3x3multiplyRowByScalarTestArrayMultiplied [9] = {5,-3.54,0.3,9,2,0,-1,-4.675,2.3};
    std::size_t rowScalarMultipliedIndex {1};
    double testScalarMultiple {-4.53};
    for(std::size_t i{0}; i < 3; ++i){
        std::size_t offset{3 * rowScalarMultipliedIndex};
        m3x3multiplyRowByScalarTestArrayMultiplied[i + offset] = m3x3multiplyRowByScalarTestArray[i + offset] * testScalarMultiple;
    }
    matrix<3, 3> m3x3multiplyRowByScalarTestMatrix {m3x3multiplyRowByScalarTestArray, 9, false};
    matrix<3, 3> m3x3multiplyRowByScalarTestMatrixMultiplied {m3x3multiplyRowByScalarTestArrayMultiplied, 9, false};
    m3x3multiplyRowByScalarTestMatrix.multiplyRowByScalar(rowScalarMultipliedIndex, testScalarMultiple);
    if(m3x3multiplyRowByScalarTestMatrix == m3x3multiplyRowByScalarTestMatrixMultiplied){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }


    std::cout << "Testing 3 x 3 addRow(): ";
    numberOfTests++;
    double m3x3addRowTestArray [9] = {0,4,-3.4,5.56,-9.542,2,-1,0,8.5};
    double m3x3addRowTestArrayAdded [9] = {0,4,-3.4,5.56,-9.542,2,-1,0,8.5};
    std::size_t rowAddTargetIndex {0};
    std::size_t rowAddSourceIndex {2};
    for(std::size_t i {0}; i < 3; ++i){
        std::size_t offsetT = rowAddTargetIndex * 3;
        std::size_t offsetS = rowAddSourceIndex * 3;
        m3x3addRowTestArrayAdded[i + offsetT] = m3x3addRowTestArray[i + offsetS] + m3x3addRowTestArrayAdded[i + offsetT];
    }
    matrix<3, 3> m3x3addRowTestMatrix {m3x3addRowTestArray, 9, false};
    matrix<3, 3> m3x3addRowTestMatrixAdded {m3x3addRowTestArrayAdded, 9, false};
    m3x3addRowTestMatrix.addRow(rowAddTargetIndex, rowAddSourceIndex); //scalar is 1, no need for complication as scalar multiplication covers that case
    if(m3x3addRowTestMatrix == m3x3addRowTestMatrixAdded){
        std::cout << "Success!\n";
        score++;
    }else{
        std::cout << "Failure\n";
    }




    





    std::cout << "\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\nEND MATRIX TESTS.\n――――――――――――――――――――――――――――――――――――――――――――――――――――――――\n";


    return score == numberOfTests;
}

void ece206mobius1helper(){
    using namespace nm::linalg;

    std::cout << "first 3d vector (A):\n";
    double a[3] = {0,0,0};
    std::cin >> a[0];
    std::cin >> a[1];
    std::cin >> a[2];
    vector<3> v_a {a, 3};
    v_a.print();

    std::cout << "second 3d vector (B):\n";
    double b[3] = {0,0,0};
    std::cin >> b[0];
    std::cin >> b[1];
    std::cin >> b[2];
    vector<3> v_b {b, 3};
    v_b.print();


    vector<3> difference = v_b - v_a;
    std::cout << "(b - a):\n";
    difference.print();

    double v = a[2];
    double u = a[0]/std::cos(v);


    double n[3] = {std::sin(v),
                        -1*(std::cos(v)),
                        u};

    vector<3> normal {n, 3};


    std::cout << "_________________________________\n" << "dot product: " << difference*normal << "\n";
}




int main(){

    

    std::cout << basicTests() << "\n";
    std::cout << vectorTests() << "\n";
    std::cout << matrixTests() << "\n";

    //std::size_t l{static_cast<std::size_t>((1-2))};
    //std::cout << l + 2 << "\n";

    //nm::printBitsD(1000.000000045);
    //nm::printBitsD(1.0/0.0);
    //nm::printBitsD(1.0/(-0.0));

    

    



	std::cout << "sizeof of std::size_t: " << sizeof(std::size_t) << "\n";
	std::cout << "sizeof of char: " << sizeof(char) << ", sizeof of float: " << sizeof(float) << ", sizeof of double: " << sizeof(double) << "\n";
	std::cout << "sizeof of uin8_t: " << sizeof(uint8_t) << ", sizeof of uint16_t: " << sizeof(uint16_t) << ", sizeof of uint32_t: " << sizeof(uint32_t) << "\n";   
    std::cout << (char)(((uint8_t)(0-1)) ^ 0x80) << "\n";
	
}
