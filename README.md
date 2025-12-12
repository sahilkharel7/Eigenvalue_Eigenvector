project:
  name: Eigenvalues and Eigenvectors Calculator
  platform: TI-84 Plus CE Python
  purpose: >
    Compute integer eigenvalues and eigenvectors of square matrices
    using exact arithmetic on the TI-84 Plus CE Python calculator.

overview:
  description: >
    This program calculates integer eigenvalues and their corresponding
    eigenvectors for 3x3 and 4x4 matrices. It avoids floating-point
    arithmetic to ensure exact results that match textbook solutions.
  supported_matrix_sizes:
    - 3x3
    - 4x4
  input_type: Integer matrices only
  output_type: Integer eigenvalues and integer eigenvector bases

features:
  - Exact integer arithmetic
  - No floating-point rounding errors
  - No complex-number artifacts
  - TI-84 Python compatible
  - Prints eigenvector bases
  - Handles repeated eigenvalues

mathematical_background:
  eigenvalue_definition: >
    A scalar lambda is an eigenvalue of matrix A if det(A - lambda I) = 0.
  approach:
    - Test integer lambda values in a user-defined range
    - Compute determinant using the Bareiss algorithm
    - Solve (A - lambda I)x = 0 using exact rational row reduction
    - Output a basis for the eigenvector space
  algorithms_used:
    determinant: Bareiss algorithm (fraction-free)
    linear_system: Rational row-reduced echelon form

usage:
  steps:
    - Open the Python app on the TI-84 Plus CE
    - Create a new script
    - Paste the program code
    - Run the script
  inputs:
    matrix_dimension:
      description: Enter the size of the matrix
      allowed_values: [3, 4]
    matrix_entries:
      description: >
        Enter each matrix entry when prompted.
        All entries must be integers.
    lambda_search_range:
      description: >
        Range of integer values to test as possible eigenvalues.
      recommended_default:
        from: -10
        to: 10

example:
  matrix:
    - [3, -2, 4]
    - [-2, 6, 2]
    - [4, 2, 3]
  lambda_range:
    from: -10
    to: 10
  output:
    - eigenvalue: -2
      eigenvectors:
        - [1, 1, -2]
    - eigenvalue: 7
      eigenvectors:
        - [1, -1, 0]
        - [2, 0, -1]

output_format:
  description: >
    For each eigenvalue found, the program prints the eigenvalue
    followed by one or more eigenvectors forming a basis.
  notes:
    - Multiple eigenvectors indicate higher geometric multiplicity
    - Vectors are scaled to integers and reduced

limitations:
  - Only integer eigenvalues are detected
  - Eigenvalues must lie within the chosen search range
  - Non-integer and complex eigenvalues are not shown
  - Designed for educational use, not symbolic computation

recommended_use_cases:
  - Linear Algebra coursework
  - Exam preparation
  - Homework verification
  - TI-84 Python environments without CAS

author_notes:
  design_goals:
    - Exact results
    - Calculator compatibility
    - Conceptual clarity
  rationale: >
    Floating-point methods introduce rounding and complex-number noise.
    This implementation ensures outputs align with classroom expectations.

license:
  type: Educational
  permissions:
    - Free to use
    - Free to modify
    - Free to extend
