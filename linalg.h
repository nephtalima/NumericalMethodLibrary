#include <cstddef>

//209 893 2936 (some spam number)

namespace nm{


    constexpr double globalTolerance {0.001};

    /*
     *  Checks if the double is equal to 0. 
     *  Both -0 and +0 evaluate to true.
     *  Returns false if the entry is non-zero, true otherwise.
     */
    //static bool isZero(double entry);
    static bool isZero(double entry, double tolerance = globalTolerance);

    /*
     *  Checks if two doubles are equal to each other within a certain tolerance.
     *  Returns true if they are equal, false otherwise.
     */
    static bool isEqual(double a, double b, double tolerance = globalTolerance);


    /*
     *  Checks the given value to check if it is NaN or infinity as defined
     *  by the IEEE 754 Double-Precision Floating-Point standard.
     *  Returns false if the value is invalid, and true otherwise.
     */
    static bool checkIfValid(double entry);

    /*
     *  Prints the sign, exponent, and matissa of a given IEEE 754 Double-Precision Floating-Point number
     */
    static void printBitsD(double entry);

    namespace linalg{

        template <std::size_t N> 
        class vector{
            private:
                /*
                *  The underlying entries of the vector.
                *  This may be changed to type std::array or std::vector at a later date.
                */
                double _entries[N]{0};

                std::size_t _density{0};

            public:
                /*
                *  Default constructor. Generates the class with all entries being 0.
                */
                vector();

                /*
                *  Constructor with pre-defined values.
                *  This function will throw an error if capacity is not equal to N,
                *  but will truncate the data down to the vector's size should you 
                *  choose to supply a bigger array.
                *  If entries is nullptr, the vector will be a zero-vector (i.e. all
                *  all entries are 0).
                */
                vector(double* const entries, std::size_t capacity);

                /*
                *  Checks the entire vector to check if any values are NaN or infinity
                *  as defined by the IEEE 754 Double-Precision Floating-Point standard.
                *  Returns false if any values are invalid, and true otherwise.
                */
                bool checkIfValid() const;

                /*
                *  Prints the vector to the console on a single line.
                */
                void print() const;

                
                /*
                *  Checks if the vector is the N-dimensional zero vector.
                *  Both -0 and +0 are counted as 0.
                *  Returns false if any entry is non-zero, true otherwise.
                */
                bool isZeroVector() const;
                

                /*
                *  Sets an entry to the provided value at the specified index.
                *  If the entry is invalid, no data is written and this function
                *  throws std::invalid_exception.
                */
                void setValue(double entry, std::size_t index);

                /*
                *  Sets entries with an array full of values, where offset
                *  specifies the amount of entries to offset from 0 (default is 0)
                *  e.g.: offset of 4 means that the array will write to position 3
                *  (0-based index) onwards.
                *  The provided array must be at least size capacity = N - offset,
                *  otherwise this function will read into forbidden memory.
                *  Moreover, if the pointer passed is nullptr or any of the values
                *  are NaN or infinity, this function will throw std::invalid_exception.
                */
                void setValues(double* const entries, std::size_t offset = 0, std::size_t capacity = N); //untested

                /*
                *  Gets a value from a specified index.
                *  The index must be from 0 to N-1, otherwise this function will
                *  throw std::invalid_argument.
                *  Returns a double from the vector at the index specified.
                */
                double getValue(std::size_t index) const;

                /*
                *  Gets all values from the vector.
                *  Returns a pointer pointing to the array of size N-1 (0-based index).
                */
                const double* const getValues() const;

                /*
                *  Returns the amount of non-zero entries in the vector.
                */
                std::size_t getDensity() const;

                /*
                *  Gets the normalized density of the vector.
                *  If M is the number of non-zero entries in the vector, this function
                *  returns M/N (and the value is between 0 and 1 inclusive).
                */
                double getDensityNormalized() const;



                /** Vector Operations **/

                /*
                *  Scalar multiplcation by scalar.
                *  Returns an N-dimensional vector scaled by the scalar.
                */
                vector<N> operator*(double scalar) const;

                /*
                *  Dot product by another N-dimensional vector.
                *  Returns the dot product of the two vectors.
                */
                double operator*(vector<N>& other) const;

                /*
                *  Vector addition by another N-dimensional vector.
                *  Returns the sum of the two vectors.
                */
                vector<N> operator+(vector<N>& other) const;

                /*
                *  Vector subtraction by another N-dimensional vector.
                *  Returns the sum of the two vectors.
                */
                vector<N> operator-(vector<N>& other) const;

                /*
                *   Tests if two vectors are equal to each other.
                */
                bool operator==(vector<N>& other) const;

                bool isEqual(vector<N>& other, double tolerance = nm::globalTolerance) const;

                /*
                *  Returns the dimension of the vector (always N).
                */
                std::size_t dimension() const;

                /*
                 *  Returns the 2-norm of the vector.
                 */
                double norm2() const;

                /*
                 *  Returns the leading entry of the vector.
                 */
                double leadingEntry() const;

                /*
                 *  Returns the index of the leading entry of the vector.
                 */
                std::size_t leadingEntryIndex() const;


        };


        template <std::size_t M, std::size_t N>
        class matrix {
        private:
            /*
             * M rows of N-dimensional row vectors
             * a.k.a. M systems of equations with N variables
             * Default when _columnMajorOrder == false
             */
            vector<N> _systems[M]{};

            /*
             *  N columns of M-dimensional column vectors
             *  a.k.a. N variables with M systems of equations
             *  Default when _columnMajorOrder == true
             */
            vector<M> _variables[N]{};

            /*
             *  List of the entries as provided in the constructor.
             */
            double* _entries;

            /*
             *  Determines how the matrix stores the entries from
             *  _entries into either _systems or _variables.
             *  If true, by default populates _variables, and
             *  if false, by default populates _systems.
             *  Note that the constructor automatically populates
             *  both _variables and _systems, but is designed for
             *  the above statement. 
             * 
             *  Example: [1,2,3,4] in a 2x2 matrix will be stored as:
             * 
             *  │ 1 3 │                │ 1 2 │
             *  │ 2 4 │                │ 3 4 │
             *  
             *  (when true)            (when false)
             *  (column-major)         (row-major)
             */
            bool _columnMajorOrder;

            /*
             *  Gets the determinant from a 2x2 matrix or 2x2 submatrix.
             *  Returns canonical ad - bc for a 2x2 matrix as follows:
             *  │ a b │ 
             *  │ c d │
             *  Note that row1, row2, column1, column2 are 0-based indices.
             */
            double _square2Determinant(std::size_t row1, std::size_t row2, std::size_t column1, std::size_t column2) const;  //untested

        public:

            /*
             *  Constructor.
             *  Populates the matrix with the array of doubles supplied.
             *  The entries will be stored in column-major order if
             *  columnMajorOrder is true, and row-major order otherwise.
             *  Example: [1,2,3,4] in a 2x2 matrix will be stored as:
             * 
             *  │ 1 3 │                │ 1 2 │
             *  │ 2 4 │                │ 3 4 │
             *  
             *  (when true)            (when false)
             *  (column-major)         (row-major)
             * 
             *  Throws std::invalid_argument if
             *  - entries is nullptr
             *  - capacity is 0
             *  - either M or N are 0
             *  - capacity is not equal to M*N
             *  
             */
            matrix(double* entries, std::size_t capacity, bool columnMajorOrder); 

            /** Low-level Operations **/

            /*
             *  Returns the value located at the i-th row and j-th column.
             * (i and j are 0-based indices)
             */
            double search(std::size_t rowIndex, std::size_t columnIndex) const;  //untested

            void setEntry(std::size_t rowIndex, std::size_t columnIndex, double value);   //untested

            /*
             *  Prints the matrix to the standard output stream.
             */
            void print() const;

            /*
             *  Returns the system (or row) specified by the zero-based index.
             */
            vector<N> getSystemByIndex(std::size_t index) const;

            /*
             *  Returns the system (or row) specified by the zero-based index. Alias of matrix<M,N>::getSystemByIndex
             */
            vector<N> getRowByIndex(std::size_t index) const;

            /*
             *  Returns the variable (or column) specified by the zero-based index.
             */
            vector<M> getVariableByIndex(std::size_t index) const;
            
            /*
             *  Returns the system (or row) specified by the zero-based index. Alias of matrix<M,N>::getSystemByIndex
             */
            vector<M> getColumnByIndex(std::size_t index) const;

            /*
             *  Converts between row-major and column-major forms.
             *  Automatically called in the constructor and thus
             *  both row-major and column-major arrays are populated.
             */
            void convert();  //untested

            /*
             *  Checks if all entries are valid.
             *  Returns false if any entries are invalid, true otherwise.
             */
            bool checkIfValidMatrix() const;  //untested
            
            /*
             *  Checks if the i-th (0-based) row is non-zero.
             *  Returns false if any entries are non-zero, true otherwise.
             */
            bool isRowEmpty(std::size_t row) const;  //untested

            /*
             *  Checks if the j-th (0-based) column is non-zero.
             *  Returns false if any entries are non-zero, true otherwise.
             */
            bool isColumnEmpty(std::size_t column) const;  //untested
            double leadingEntry(std::size_t row) const; //untested
            bool isIdentityMatrix() const; //untested
            //bool isEmptyMatrix() const; //untested
            //bool isElementaryMatrix() const; //untested
            
            

            /** Basic Matrix Operations **/

            //by default makes row-major matrix
            matrix<M,N> operator+(matrix<M,N>& other) const; //untested


            /*matrix multiplication*/
            template <std::size_t K> matrix<M,K> operator*(matrix<N,K>& other) const; //untested

            /*
             *  Transposes the matrix into an NxM matrix.
             *  Example:
             *  (Original)
             *  │ 1 3 │    
             *  │ 2 4 │ 
             *  
             *  (Transposed)
             *  │ 1 2 │
             *  │ 3 4 │
             */
            matrix<N,M> operator~() const; //untested


            bool operator==(matrix<M,N>& other) const; //untested

            bool isEqual(matrix<M,N>& other, double tolerance = nm::globalTolerance) const;

            /** Elementary Row Operations **/

            void swapRows(std::size_t row1Index, std::size_t row2Index); //untested
            void multiplyRowByScalar(std::size_t rowIndex, double scalar); //untested
            void addRow(std::size_t targetRowIndex, std::size_t sourceRowIndex, double scalar = 1); //untested

            


            void gaussJordanElimination(vector<M> b) const; //untested
            //double determinant() const; //untested

            

        };


    }
}