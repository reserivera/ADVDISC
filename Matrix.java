/**
 * @author BATOSALEM, Angelika && CORTEZ, Louise && RIVERA, Sophia
 * @section S17
 */

import java.util.ArrayList;

public class Matrix {
    // The usage of an Array/List-like structure to store Matrix data 
    // as a list of Vectors. You may also store the Matrix as a 2d array.
    private ArrayList<Vector> vectors;

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
    public static Matrix times (Matrix other) {
        // Errors for size mismatches when multiplying matrices must also be handled.
        return null;
    }

    // An implementation of a function that performs Gauss-Jordan 
    // Elimination to find the determinant of the matrix.
    public static double det () {
        return 0.0;
    }

    // An implementation of a function that finds the inverse of the matrix.
    public static Matrix inverse () {
        // The function must return a null value if the matrix has no inverse.
        return null;
    }
}