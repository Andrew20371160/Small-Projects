# 0) Quick start

## Include
```cpp
#include "matrix_algebra/matrix_algebra.h"
```

## Main types
- `complex` (custom complex number)
- `template<typename DataType> class matrix`

## Key notes / conventions
- Many error paths print to `cout` and return **“garbage default”** values like:
  - `matrix<DataType>(1,1,-1)` for matrix-returning functions
  - `-1` for scalar-returning functions
- Numeric thresholds:
  - `tolerance = 1e-32` (for computations: gauss/inverse/gram-schmidt)
  - `check_tolerance = 1e-6` (for `is_*` checks)

---

# 1) Constants, enums, and global helpers

## 1.1 Sentinel
- `#define no_init INT_MIN`  
  Used in `matrix(rows, cols, initialization_value, ...)`: if the initialization value equals `no_init`, the constructor allocates but **does not fill**.

## 1.2 Global error strings
- `shape_error`
- `square_error`
- `uninit_error`

Printed via `cout` when operations fail.

## 1.3 Enums (public API)

### Placement / splitting
- `upper_left, lower_left, lower_right, upper_right`  
  Used by `matrix::at_quarter(...)`.

- `lower_half, upper_half`  
  Used by `matrix::split(...)` and FFT helper logic.

### Pivot mapping mode
- `old_locations, new_locations`  
  Used by `matrix::gauss_down(...)` to decide whether pivot indices are returned in:
  - *post-swap* row order (`new_locations`)
  - or tracking *original row locations* via row swaps (`old_locations`) (used by `det()` sign logic)

### Matrix feature types
`general, utri, ltri, diagonal, symmetric, anti_symmetric, constant, iden, orthonormal`

### Compression types
`general_compress, utri_compress, ltri_compress, diagonal_compress, symmetric_compress, anti_symmetric_compress, constant_compress, iden_compress, orthonormal_compress`

Important: internally, the library compares `compression_type` and the requested `matrix_type`; if they don’t match it may **decompress and fall back to general access**.

---

# 2) Complex number module (`complex.h/.cpp`)

## 2.1 Class `complex`
Represents a complex number `(re + im*j)` with:
- `long double re, im`

### Constructors / destructor
- `complex(long double _re=0, long double _im=0)`
- `complex(int)` (real = int, imag = 0)
- copy constructor `complex(const complex&)`
- destructor `~complex()`

### Accessors
- `long double get_re() const`
- `long double get_im() const`

### Arithmetic
- `operator+`, `operator-`, `operator*`, `operator/` (complex-complex)
- `operator*(long double)`, `operator/(long double)`
- compound: `+=, -=, *=, /=`
- power:
  - `complex operator^(const long double&) const`
  - `void operator^=(const long double&)`

### Assignment
- `void operator=(const complex&)`
- `void operator=(const long double&)` (sets imag to 0)

### Comparisons
- `>`, `<` compare by magnitude (`abs`)
- `==` compares componentwise within `tolerance`
- `!=` is negation of `==`

### Angle
- `long double theta() const` uses `atan2(im, re)`

### Output
- `friend ostream& operator<<(ostream&, const complex&)` prints `re` and `±imj` if `|im| > tolerance`.

## 2.2 Free functions (helpers)
- `long double abs(const complex&)`
- `string to_string(complex&)`  
  Used by `matrix::mat_to_string()` for `DataType=complex`.

### Conjugation overload family (very important)
- `complex conjugate(const complex&)` → `(re, -im)`
- For numeric primitive types (`int`, `double`, etc.) there are overloads returning the same value.

This design makes:
- `matrix::transpose()` behave as **conjugate-transpose** for complex matrices
- `matrix::dot()` accumulate `mat(i,j) * conjugate(this(i,j))` (a Hermitian-style dot for complex)

### `pow` wrapper
- `complex pow(const complex& c, int power)` calls `c ^ power`

---

# 3) Core matrix type (`matrix_algebra.h/.cpp`)

## 3.1 Class template `matrix<DataType>`

### Storage model
- General matrices store all elements row-major in `vec`.
- Compressed matrices store only required elements (triangular / packed symmetric / etc.) and use:
  - `row_start[]`, `pindex[]` (for `utri`/`ltri` mapping)
  - packed indexing for symmetric / anti-symmetric
  - singleton `vec[0]` for constant and identity
- Element access is dispatched through member-function pointers:
  - `at_ptr` and `at_ptr_c` for mutable/const access
  - `start_ptr` and `end_ptr` for “nonzero range” in each row, used to accelerate loops in multiplication and some vector ops.

### “Feature” vs “Compression”
- `matrix_type` = what the matrix *represents* (general, diagonal, utri, …)
- `compression_type` = how it is *currently stored* (general_compress, diagonal_compress, …)
- `is_compressed` = storage state

The library sometimes sets `matrix_type` even when not compressed; you can then call `compress()` to pack storage.

---

# 4) Construction, destruction, basic access (Foundation section)

## 4.1 Constructors / destructor
- `matrix()`  
  Creates an empty matrix with `rows=0, cols=0`, no storage.

- `matrix(int rows, int cols, DataType initialization_value = no_init, bool compressed=false)`  
  Allocates a matrix. Behavior:
  - If `compressed==true`: allocates `1` element and marks matrix as `constant` storage (stores `initialization_value` if not `no_init`)
  - Else: allocates full `rows*cols`, fills only if init value is not `no_init`

- `matrix(int rows, int cols, DataType* initialization_array, int size_of_array)`  
  Fills in row-major from the provided array up to `size_of_array`, pads the rest with zero.

- copy constructor `matrix(const matrix&)`
- destructor `~matrix()` frees `vec`, `row_start`, `pindex`

## 4.2 Basic getters
- `int get_rows() const`
- `int get_cols() const`
- `int get_size() const`  
  Returns `actual_size` (stored elements, not always `rows*cols`).
- `int get_type() const`  
  Returns `matrix_type` enum.

## 4.3 Indexing and validation
- `bool is_valid_index(int row, int col) const`  
  Checks bounds.

- `DataType& at(int r, int c)`
- `const DataType& at(int r, int c) const`  
  Dispatches to the correct accessor for the current matrix feature/compression.

### Important caveat (design)
For compressed triangular/diagonal/identity/constant matrices, out-of-structure elements return a reference to `empty_vec[0]` (often set to 0). Writing through that reference does **not** store into the matrix (it only overwrites the dummy).

## 4.4 Display / conversion
- `void show() const` prints rows like `[a , b , c]`
- `string mat_to_string() const`  
  Uses `to_string(at(i,j))`. For custom `DataType`, you must provide `to_string(DataType&)` overload.

## 4.5 Basic mutators
- `void fill(DataType value)`  
  Fills every entry with `value` and sets `matrix_type = constant` (note: not automatically compressed).
- `void set_identity()`  
  Fills identity (square only) and sets `matrix_type = iden` (not automatically compressed).

---

# 5) Algebraic operations (Arithmetic section)

## 5.1 Shape checks
- `bool same_shape(const matrix&) const`

## 5.2 Operators
- `matrix operator+(const matrix&) const`
- `matrix operator-(const matrix&) const`
- `matrix operator*(const matrix&) const`  
  Matrix multiplication. Requires `this->cols == mat.rows`.  
  Uses `start(row)`/`end(row)` for the left matrix to reduce loops on structured matrices.

- `matrix operator/(const matrix&) const`  
  Defined as `A / B = A * inverse(B)`.

- `matrix operator*(DataType scalar) const`  
  Scalar multiply.

- `matrix operator^(const DataType& power) const`  
  Attempts diagonalization-like power:
  1. compute eigenvalues via QR iteration (`eigen_values`)
  2. compute eigenvectors (`eigen_vectors`)
  3. return `S * (Λ^power) * S^-1` (implemented as scaling columns by diagonal entries)

  Warning: depends on eigen routines that may diverge or be inaccurate.

## 5.3 Vector-style helpers
- `DataType dot(const matrix&) const`  
  Requires same shape. Accumulates:
  `sum += mat(i,j) * conjugate(this(i,j))`.

- `DataType axpy(DataType alpha, const matrix& y) const`  
  Computes a scalar accumulator of `alpha*this + y` (note: **returns a scalar**, not a matrix).  
  Iterates from `max(start(i), y.start(i))` to `min(end(i), y.end(i))`.

- `DataType norm2() const` / `DataType length() const`  
  `sqrt(abs(dot(*this)))`.

- `DataType theta(matrix& other) const`  
  Uses cosine formula with `acos(abs(adotb/(|a||b|)))` and returns degrees.

---

# 6) Structural predicates and special matrix properties (Matrix classification section)

All use `check_tolerance` for numeric decisions unless otherwise indicated.

## 6.1 Shape-based
- `bool is_square() const`

## 6.2 Equality / zero / identity
- `bool operator==(const matrix&) const`
- `bool is_zero() const`  
  If true, sets `matrix_type = constant`.

- `bool is_identity() const`  
  If true, sets feature to `iden`.

## 6.3 Triangular / diagonal / scalar / symmetry
- `bool is_upper_tri() const`  
  If true, sets feature to `utri`.

- `bool is_lower_tri() const`  
  If true, sets feature to `ltri`.

- `bool is_diagonal() const`  
  Implemented as `is_upper_tri() && is_lower_tri()`. If true sets `diagonal`.

- `bool is_scalar() const`  
  Checks diagonal all equal and off-diagonal all zero.

- `bool is_symmetric() const`  
  Checks `a(i,j) == conjugate(a(j,i))`. If true sets `symmetric`.

- `bool is_anti_symmetric() const`  
  Checks `a(i,j) == -conjugate(a(j,i))`. If true sets `anti_symmetric`.

## 6.4 Orthogonality / parallelism
- `bool is_orthogonal() const`  
  Checks `A * A^T` is identity, sets `orthonormal` on success.

- `bool is_orthogonal(const matrix&) const`  
  Checks `transpose() * mat` is zero (requires same number of rows).

- `bool is_parallel(const matrix&) const`  
  True if `theta≈0` or `theta≈180`.

## 6.5 Polynomial-like properties
- `bool is_idempotent() const` checks `A*A == A`
- `bool is_involutory() const` checks `A*A` is identity
- `bool is_nilpotent() const`  
  Repeatedly squares `mat_pow = mat_pow * mat_pow` up to `rows` times and checks if it becomes zero. (This is **not** the usual A^k iteration; it squares each time: A, A^2, A^4, …)

---

# 7) Gaussian elimination, pivots, rref, rank (Elimination & solving section)

## 7.1 Helper row operation
- `void row_axpy(DataType alpha, int x_row, int y_row)`  
  Performs `row_y += alpha * row_x` (over columns from `start(y_row)` to `end(y_row)`).

## 7.2 Pivot detection helpers
- `int is_pivot(int row_index, int col_index)`  
  Searches **downward in the same column** for a nonzero entry.
  - If found at a lower row, it **switches rows** and returns that pivot row index.
  - If none found, returns `-1`.

- `int is_pivot_up(int row_index, int col_index)`  
  Searches **upward** in the same column for a nonzero entry.
  - Does **not** switch rows; returns the row index where pivot was found or `-1`.

## 7.3 Row swap
- `bool switch_rows(int row1, int row2)`  
  Swaps full rows in-place.

## 7.4 `gauss_down` (forward elimination)
- `matrix gauss_down(matrix<int>* pivots_indices = NULL, int pivots_locations = new_locations) const`

Produces an upper-triangular-like form:
- If called on a matrix already marked `utri`, it returns itself and optionally copies `pindex`.
- Otherwise:
  - performs forward elimination
  - returns a matrix whose feature is set to `utri`
  - optionally returns pivot indices:
    - `new_locations`: `pivots_indices(i,0)` = pivot column index for row i, or -1
    - `old_locations`: `pivots_indices` is initialized as `[0,1,2,...]` and row swaps are recorded by swapping entries in that vector

## 7.5 `gauss_up` (backward elimination)
- `matrix gauss_up(matrix<int>* pivots_indices = NULL) const`

Produces a lower-triangular-like form:
- If matrix already `ltri`, returns itself (optionally copies `pindex`)
- Else performs backward elimination (from bottom) and returns feature `ltri`

## 7.6 Back/forward substitution helpers
These assume you are solving an appended system-like matrix (augmented with RHS in last column).

- `DataType back_sub(int selected_row, const matrix& solution_matrix) const`  
  Uses last column as RHS and divides by diagonal element.

- `DataType fwd_sub(int selected_row, const matrix& solution_matrix) const`  
  Forward-sub variant.

## 7.7 `solve()` (augmented system solver)
- `matrix solve() const`

Assumes caller is an augmented matrix `[A | b]` with `rows` equations and `cols = rows + 1` (typical usage).
- Runs `gauss_down` and checks every row has a pivot.
- Performs back substitution to return `rows x 1` solution.

If pivots insufficient, returns `matrix(1,1,-1)`.

## 7.8 `rref()` (reduced row echelon form)
- `matrix rref(matrix<int>* pivots_indices = NULL) const`

Implementation:
1. `gauss_down` to get pivot locations by row
2. eliminates above pivots
3. normalizes pivot rows to pivot = 1

Returns feature set to `utri` (even though it’s rref-like).

## 7.9 `rank()`
- `int rank() const`
Uses `gauss_down` to count pivots until a `-1` pivot index occurs.

---

# 8) Four Fundamental Subspaces

This section groups:
- column space (range)
- row space
- null space (kernel)
- left null space

And documents all helper functions used for these.

## 8.1 Dimension helpers
- `int dim() const`  
  Returns `rank()` (dimension of column space).

- `int dim_null_cols() const`  
  Returns `cols - rank()` (nullity).

- `int dim_null_rows() const`  
  Returns `rows - rank()` (left null space dimension).

## 8.2 Independence / basis checks
- `bool is_independent() const`  
  True if `rank() == cols`.

- `bool is_basis(int dimension) const`  
  True if `cols == dimension && is_independent()`.

## 8.3 Basis for column space
- `matrix basis_cols() const`

Algorithm:
1. compute pivot column indices from `gauss_down`
2. take the corresponding pivot columns from the **original matrix**
3. return a matrix whose columns form a basis of `Col(A)`

Return shape: `rows x rank(A)`.

## 8.4 Basis for row space
- `matrix basis_rows() const`

Algorithm:
1. compute `rref()` and count pivot rows
2. return the first `rank(A)` rows of rref (zeroed where very small)

Return shape: `rank(A) x cols`.

## 8.5 Null space of columns: `Null(A)`
- `matrix null_cols() const`

Uses `rref` and then constructs special-solution vectors by choosing free variables.
Return shape: `cols x nullity`.

- If full rank (no free vars), returns `cols x 1` zero vector.

### Helper used here
- Internally relies on pivot bookkeeping from `rref(...)`.

## 8.6 Left null space: `Null(A^T)` (null space of rows)
- `matrix null_rows(matrix* elementary = NULL) const`

Despite the comment “null space of row space”, the returned vectors correspond to the left-null-space basis (vectors `y` such that `y^T A = 0`), built from elimination operations tracked on an elementary matrix.

High-level algorithm (as coded):
1. Copy matrix and call `fix_pivots()` to stabilize pivot row ordering
2. Perform elimination while simultaneously updating an `elementary` matrix (starts as identity)
3. After reaching rref-like form, rows of `elementary` corresponding to zero rows in the transformed `A` form left-null-space vectors

Return shape: `(rows - rank(A)) x rows` (basis row vectors).
- If no left null space dimension, returns `matrix(1, rows, 0)`.

### Optional output helper
- If `elementary != NULL`, writes the computed elementary matrix into `*elementary`.

## 8.7 Pivot-row ordering helper
- `void fix_pivots()`
Uses `gauss_down(... old_locations)` and then performs row switches to move pivot rows upward in a consistent way, to avoid issues with elementary-matrix tracking.

---

# 9) Determinant, inverse, trace, cofactors (Square-matrix analytics section)

## 9.1 `trace()`
- `DataType trace() const`
Sums diagonal (implementation uses an `int res`, which is a type bug for non-int DataType).

## 9.2 `det()`
- `DataType det() const`
Uses `gauss_down(old_locations)`:
- multiplies diagonal entries of the eliminated matrix
- attempts to correct sign based on how many rows moved (counts `original_pivots_indices(i) != i`)

## 9.3 `inverse()`
- `matrix inverse() const`

- If matrix is marked `orthonormal`, it returns `transpose()` directly (fast path).
- Else performs Gauss-Jordan-like elimination while applying the same row operations to an identity matrix.

Returns `matrix(1,1,-1)` on failure.

## 9.4 Cofactors
- `DataType cofactor(int row_i, int col_i) const`
Builds minor matrix (size `n-1`) and returns `det(minor) * (-1)^(row_i+col_i)`.

- `matrix cofactors() const`
Returns cofactor matrix (same shape).

---

# 10) Factorizations and decompositions (LU / QR / SVD section)

## 10.1 LU factorization with permutation
- `void lu_fact(matrix& lower_fact, matrix& permutation, matrix& upper_fact) const`

Behavior:
- If caller is already `utri`: returns `L=I`, `P=I`, `U=this`.
- Else (square only):
  - performs elimination on `upper_fact`
  - stores multipliers in `lower_fact`
  - records row swaps into `permutation` (by switching identity’s rows)
  - sets features: `lower_fact` as `ltri`, `upper_fact` as `utri`

## 10.2 Gram–Schmidt (orthonormalization)
- `matrix gram_shmidt() const`

Treats columns as vectors. Produces an orthonormal set (as implemented).
Returns feature `orthonormal`.

Notes from code comments:
- recommended to use `double` and very low tolerance; results can be unstable.

## 10.3 QR factorization
- `void qr_fact(matrix& q, matrix& r) const`

Implementation:
- `q = gram_shmidt()` (feature set orthonormal)
- `r = q.transpose() * (*this)` (feature set `utri`)

## 10.4 Eigenvalues via QR iteration
- `matrix eigen_values(const int& max_iteration, const double& min_diff = check_tolerance) const`

Iterates:
1. save current diagonal
2. `qr_fact(q,r)`
3. `A = r*q`
4. after halfway iterations, checks diagonal stability vs previous diagonal (min_diff)

Returns a diagonal matrix `rows x rows` with eigenvalues on diagonal. Feature set to `diagonal`.
Code comment warns it “doesn’t generate correct results all the time due to divergence issues”.

## 10.5 Eigenvectors from eigenvalues
- `matrix eigen_vectors(const matrix& eigen_values) const`

For each eigenvalue λ:
1. compute `A - λI`
2. compute `null_cols()` of that matrix
3. take the first null-space vector as the eigenvector column

Return shape: `rows x rows`.

## 10.6 SVD (Singular Value Decomposition)
- `void svd(matrix& u, matrix& s, matrix& vt) const`

Implements a rough approach:
- Compute `AAT = A * A^T`
- `s = eigen_values(AAT)` (treated as Σ² initially)
- `u = eigen_vectors(s)`
- normalize columns of `u` by `col_length`
- `vt = eigen_vectors( (A^T * A) )`-like path (code uses `(transpose()*(*this)).eigen_vectors(s)` then transposes)
- set `s(i,i) = sqrt(s(i,i))`
- features set: `u` orthonormal, `s` diagonal, `vt` orthonormal

---

# 11) Least squares, projection, and linear systems (Applied linear algebra section)

## 11.1 Projection matrix
- `matrix projection() const`

Returns:
\[
P = A (A^T A)^{-1} A^T
\]
Sets feature to `symmetric`.

## 11.2 Least-squares fit
- `matrix fit_least_squares(const matrix& output) const`

Given `A` (caller) and `b` (`output`), returns:
\[
x^* = (A^T A)^{-1} A^T b
\]
Implementation builds `[A^T A | A^T b]` and calls `solve()`.

## 11.3 Transformer (map input vectors to output vectors)
- `matrix transformer(const matrix& output) const`

Goal: find matrix `T` such that:
\[
T \cdot \text{input} = \text{output}
\]
where “input” are columns of the caller, “output” are columns of the parameter.

Implementation:
- builds an augmented system from `transpose().append_cols(output.transpose())`
- `gauss_down` once
- then back-substitutes for each output column to get columns of the transformer
- returns `solution.transpose()`

Requires that pivot indices match identity positions (expects independent columns); otherwise returns error matrix.

---

# 12) FFT / Fourier utilities (Signal processing section)

## 12.1 Fourier matrix
- `matrix<complex> fourier_mat(int dimension)`
Builds a symmetric-filled Fourier matrix using:
- `w = cos(2π/N) + j sin(2π/N)`
- `ret(i,j) = w^(i*j)` filled for j>=i and mirrored

## 12.2 Fourier diagonal vector
- `matrix<complex> fourier_diagonal(int dimension, int n)`
Returns column vector `[w^0, w^1, ..., w^(n-1)]^T`.

## 12.3 Rearrangement helper
- `matrix arrange(const matrix<int>& seq) const`
Permutes rows according to `seq(i,0)`.

Used by FFT to reorder even/odd indices.

## 12.4 Quarter placement helper
- `void at_quarter(int quarter, const matrix& src)`
Copies `src` into one quarter of caller (upper-left/lower-left/upper-right/lower-right), if it fits within half dimensions.

## 12.5 Split
- `matrix split(int half) const`
Returns half the rows (top or bottom half, per code’s `half` logic).

## 12.6 FFT of a column
- `matrix fft_col(int dimension) const`
Recursive Cooley–Tukey style (expects length power-of-2; caller `fft()` pads columns).

## 12.7 FFT of a matrix (column-wise)
- `matrix fft() const`
For each column:
- pads to next power of 2 (if needed)
- runs `fft_col`
- writes back first `rows` values

---

# 13) Reshaping, extraction, filtering (Utilities section)

## 13.1 Append utilities
- `matrix append_cols(const matrix& src) const`  → `[A | B]`
- `matrix append_rows(const matrix& src) const`  → stacked rows

## 13.2 Extract column
- `matrix extract_col(int index) const`

## 13.3 Resize with padding
- `matrix resize(int wanted_rows=get_rows(), int wanted_cols=get_cols(), DataType padding_value=0) const`

Copies overlap region and fills remaining entries with padding value.

## 13.4 Filter small values
- `void filter(double filter_tolerance = check_tolerance)`
Sets entries with `abs(value) <= filter_tolerance` to 0.

## 13.5 Column length helper
- `DataType col_length(const int& col_i) const`
Computes Euclidean length of a column (note: uses `at(i,col_i)*at(i,col_i)` without conjugate; works best for real types).

---

# 14) Compression & sparse-structure acceleration (Storage/Compression section)

This library supports compressing certain structured matrices to reduce memory and speed up loops via `start(row)` / `end(row)`.

## 14.1 Public API
- `void fill_features(double check_tol=1e-6)`  
  Detects features (utri/ltri/symmetric/anti_symmetric/constant/iden/diagonal) efficiently.

- `void compress()`  
  Calls `fill_features()` then compresses accordingly:
  - symmetric / anti-symmetric: packed upper triangle
  - utri / ltri: packed by row with pivot indices
  - diagonal: main diagonal only
  - identity / constant: one element

- `void decompress()`
Re-expands into full `rows*cols` row-major storage.

- `int start(int row_i) const`, `int end(int row_i) const`
Return computation bounds per row depending on the structure:
- general/constant/symmetric: `[0, cols)`
- diagonal/identity: `[row_i, row_i+1)`
- utri: `[pindex[row_i], cols)` or empty if no pivot
- ltri: `[0, pindex[row_i]+1)` or empty if no pivot

## 14.2 Compression helpers (private but behavior matters)
- `compress_utri(const matrix<int>& pivot_index_by_row)`
- `compress_ltri(const matrix<int>& pivot_index_by_row)`
- `compress_symmetric()`
- `compress_anti_symmetric()`
- `compress_diagonal()`
- `compress_identity()`
- `compress_const()`

## 14.3 Pivot index extractor
- `matrix<int> get_pindex()`
Returns pivot positions per row used for compressing triangular matrices.

---

# 15) Free functions outside the class (creation helpers)

## 15.1 Identity matrix
- `template<typename DataType> matrix<DataType> identity(int dimension)`

## 15.2 Random matrix
- `template<typename DataType> matrix<DataType> rand(int rows, int cols, int max_val=INT_MAX)`
Fills with `rand() % max_val` (uses `srand(time(0))`).

## 15.3 Fourier helpers
- `matrix<complex> fourier_diagonal(int dimension, int n)`
- `matrix<complex> fourier_mat(int dimension)`

---

# 16) Demo reference
- There is a demo on youtube the link is in DEMO.md

---
