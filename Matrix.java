/**
 * @author BATOSALEM, Angelika && CORTEZ, Louise && RIVERA, Sophia
 * @section S17
 */

import java.util.ArrayList;
import java.util.List;

public class Matrix {
    // The usage of an Array/List-like structure to store Matrix data 
    // as a list of Vectors. You may also store the Matrix as a 2d array.
    private List<Vector> matrix;

    // The usage of immutable Integer variables to hold values for 
    // the number of rows/columns.
    private int rows;
    private int columns;
    
    private int dimension; // ?? not sure

    // A proper implementation of a default constructor that initializes 
    // the matrix as an identity matrix of a given dimension.
    public Matrix(int dimension) {
        rows = dimension;
        this.dimension = dimension;
        double[] row = new double[dimension];
        matrix = new ArrayList<Vector>();

        for(int i = 0; i < dimension; i++) {
            for(int j = 0; j < dimension; j++) {
                if(i == j) {
                    row[j] = 1;
                } else {
                    row[j] = 0;
                }
            }

            Vector v = new Vector(row, dimension);
            matrix.add(v);
        }
    }

    // A proper implementation of a constructor, converting an 
    // already-existing array/list of data from a rudimentary data 
    // structure into the vector class.
    public Matrix (List<Vector> list, int dimension) {
        columns = list.size();
        rows = dimension;
        this.dimension = dimension;
    }

    // An implementation of function for matrix multiplication.
    public double[][] times (Matrix other) {
        // Errors for size mismatches when multiplying matrices must also be handled.
        if(matrix.getNumCols() != other.getNumRows()) {
            return null;
        }
    }

    // An implementation of a function that performs Gauss-Jordan 
    // Elimination to find the determinant of the matrix.
    public double det () {
        return 0.0;
    }

    // An implementation of a function that finds the inverse of the matrix.
    public Matrix inverse () {
        // The function must return a null value if the matrix has no inverse.
        return null;
    }

    public int getNumRows() {
        return rows;
    }

    public int getNumCols() {
        return columns;
    }

    // For testing
    public static void main(String[] args) {
        Matrix m = new Matrix(3);
    }
}